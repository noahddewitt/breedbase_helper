#Dependencies
#library(BrAPI)
#library(openxlsx)
#library(breedbase)
#library(sommer)
#library(emmeans)

#Need to break up into reasonable chunks

make_lsuDB_connection <- function() {
  lsuDB <- createBrAPIConnection(host = "sgbreedbase.agcenter.lsu.edu", protocol = "http")  
}


convert_trait_if_numeric <- function(phenoVec) {
  #It's okay to have NAs. But do we generate any *new* NAs?
  if (any(!is.na(suppressWarnings(as.numeric(phenoVec[which(!is.na(phenoVec))]))))) {
    phenoVec <- as.numeric(phenoVec)
  } 
  return(phenoVec)
}

get_all_trials <- function(brapiDB) {
  allTrials <- brapiDB$get("/studies", page = "all", query = list(Location = "Baton Rouge"), pageSize = 10)
  
  allTrialIds <- c()
  for (page in allTrials$content) { 
    for (trial in page$result$data) {
      allTrialIds <- c(allTrialIds, trial$studyDbId)
      names(allTrialIds)[length(allTrialIds)] <- trial$studyName
    }
  }
  
  return(allTrialIds)
}

parse_BrAPI_JSON <- function(brapiResponse, selectedCols, selectedPlotMetadata = NULL) {
  #JSON parsed as list of lists, one per row, but with variables in random order so have to sort
  brapiResponse_orderedSubset <- lapply(brapiResponse$combined_data, 
                                        function(obsList) {unlist(obsList)[selectedCols]})
  trialPhenotypes_DF <- data.frame(do.call(rbind, brapiResponse_orderedSubset))
  
  if (!is.null(selectedPlotMetadata)) {
    for (dataLevel in selectedPlotMetadata) {
      #Order of the levelCode/Order/Name randomized within plot, but level type order consistent
      obsLevelOrder <- brapiResponse$combined_data[[1]]$observationUnitPosition$observationLevelRelationships
      obsLevelOrder_DF <- do.call(rbind, lapply(obsLevelOrder, data.frame))
      obsLevelIndex <- which(obsLevelOrder_DF$levelName == dataLevel)
      
      if (length(obsLevelIndex) == 0) {
        print("Error: selected plot metadata not found")
        break
      } else {
        brapiResponse_metaData <- lapply(brapiResponse$combined_data, 
                                         function(obsList) {obsList$observationUnitPosition$observationLevelRelationships[[obsLevelIndex]]$levelCode})
        trialPhenotypes_DF[, dataLevel] <- unlist(brapiResponse_metaData)
      }
    }
  }
  return(trialPhenotypes_DF)
}

#Different than function for trial design in the label maker file
get_trial_design <- function(trialId, brapiDB, pageSize = 40) {
  trialPlots_JSON <- brapiDB$get("/observationunits", query = list(studyDbId = trialId), page = "all", pageSize = pageSize) #, verbose = FALSE)
  
  trialPlots_DF <- parse_BrAPI_JSON(trialPlots_JSON, selectedCols =  c("observationUnitName", 
                                                                       "observationUnitDbId", 
                                                                       "trialName", 
                                                                       "observationUnitPosition.positionCoordinateX",
                                                                       "observationUnitPosition.positionCoordinateY",
                                                                       "germplasmName", 
                                                                       "observationUnitPosition.observationLevel.levelCode"),
                                    selectedPlotMetadata = c("block")
  )
  
  trialPlots_DF <- unique.data.frame(trialPlots_DF)
  
  colnames(trialPlots_DF) <- c("PLOTNAME", "PLOTID", "TRIAL", "RANGE", "ROW", "DESIG", "PLOTNUM", "REP")
  return(trialPlots_DF)
}

get_trial_phenotypes <- function(trialId, brapiDB, pageSize = 20) {
  trialPhenotypes_JSON <- brapiDB$get("/observations", query = list(studyDbId = trialId), page = "all", pageSize = pageSize)
  
  
  trialPhenotypes_DF <- parse_BrAPI_JSON(trialPhenotypes_JSON, selectedCols = c("observationUnitName",
                                                                                "studyDbId",
                                                                                "germplasmName",
                                                                                "observationVariableName", 
                                                                                "value"))
  if (length(trialPhenotypes_DF) == 0) {
    return(NULL)
  } else {
    trialPhenotypes_wideDF <- reshape(trialPhenotypes_DF, idvar = c("observationUnitName", "studyDbId", "germplasmName"), 
                                      timevar = "observationVariableName", 
                                      direction = "wide")
    
    for (colIndex in which(grepl("^value\\.", colnames(trialPhenotypes_wideDF)))) {
      trialPhenotypes_wideDF[, colIndex] <- convert_trait_if_numeric(trialPhenotypes_wideDF[,colIndex])
    }
    
    colnames(trialPhenotypes_wideDF) <- gsub("value\\.", "", colnames(trialPhenotypes_wideDF))
    
    #temporary -- fix once breedbase update fixes perl NA/0 bug
    trialPhenotypes_wideDF[is.na(trialPhenotypes_wideDF)] <- 0
    
    #Should this live here?
    coKey <- c("STAND09" = "Plant stand - 0-9 density scale", 
               "PHE09" = "Agronomic Merit - 1-9 Agronomic Merit Rating", 
               "STRIPE09" = "Stripe rust plant response - 0-9 Mc Neal scale",
               "LFRUST09" = "Leaf rust plant response - 0-9 Mc Neal scale|CO_321:0501123",
               "HD_JUL" = "Heading time - day",
               "HT_IN" = "HT_IN",
               "BYDV09" = "BYDV09",
               "BYDV09" = "BYDV09")
    
    colSubIndices <- which(grepl(paste(coKey, collapse = "|"), colnames(trialPhenotypes_wideDF)))
    
    for (i in colSubIndices) { 
      traitName <- names(coKey)[which(grepl(colnames(trialPhenotypes_wideDF)[i], coKey))]
      
      #Composed trait
      if (grepl("\\|day", colnames(trialPhenotypes_wideDF)[i])) {
        dayInfo <- regmatches(colnames(trialPhenotypes_wideDF)[i], regexpr("\\|day \\d*", colnames(trialPhenotypes_wideDF)[i])) 
        traitName <- paste0(traitName, dayInfo)
      }
      
      colnames(trialPhenotypes_wideDF)[i] <- traitName
    }
    
    
    colnames(trialPhenotypes_wideDF)[1:3] <- c("PLOTNAME", "STUDYID", "DESIG")
    
    return(trialPhenotypes_wideDF)
  }
}

get_trial_genotype_info <- function(trialId, brapiDB, pageSize = 20) {
  #Honestly surprised this works, was going to try a more complicated solution
  trialGenotypes_JSON <- brapiDB$get("/germplasm", query = list(studyDbId = trialId), page = "all", pageSize = pageSize)
  
  trialGenotypes_DF <- parse_BrAPI_JSON(trialGenotypes_JSON, selectedCols = c("germplasmName",
                                                                              "pedigree",
                                                                              "additionalInfo.additionalProps.purdy pedigree",
                                                                              "synonyms.synonym",
                                                                              "seedSource"))
  
  colnames(trialGenotypes_DF) <- c("DESIG", "PARS", "PED", "SYNS", "SRC")
  
  trialGenotypes_DF <- unique.data.frame(trialGenotypes_DF)
  trialGenotypes_DF <- trialGenotypes_DF[order(trialGenotypes_DF$DESIG), ]
  
  return(trialGenotypes_DF)
}

get_fieldbook_DF <- function(trialId, brapiDB, pageSize = 20) {
  phenoDF <- get_trial_phenotypes(trialId, brapiDB, pageSize)
  genoDF <- get_trial_genotype_info(trialId, brapiDB, pageSize)
  trialDF <- get_trial_design(trialId, brapiDB, pageSize)
  
  phenoNames <- colnames(phenoDF)[-c(1:3)]
  
  combDF <- merge(trialDF, genoDF, by.x = "DESIG")
  combDF$PLOTNUM <-  as.numeric(gsub("^.*_", "", combDF$PLOTNAME))
  combDF$LOC <-  gsub("_.*$", "", gsub("^[^_]*\\d{2}", "", combDF$PLOTNAME))
  combDF$EXPT <-  gsub("\\d{2}.*$", "", combDF$PLOTNAME)
  
  if (!is.null(phenoDF)) { 
    combDF <- merge(combDF, phenoDF[, -which(colnames(phenoDF) %in% c("DESIG", "STUDYID"))], by.x = "PLOTNAME", by.y = "PLOTNAME", all.x = T)
  }
  
  combDF <- combDF[, c("PLOTNAME", "TRIAL", "LOC", "EXPT", "RANGE", "ROW", "REP", "PLOTNUM", "DESIG", "PED", "SYNS", phenoNames)]
  combDF <- combDF[order(combDF$PLOTNAME), ]
  return(combDF)
}

format_fieldbook_DF <- function(fieldbook_DF, extraCols = 0, lessInfo = TRUE) {
  trialName <- gsub("_.*$", "", fieldbook_DF$PLOTNAME[1])
  
  if(lessInfo) {
    fieldbook_DF <- fieldbook_DF[, -which(colnames(fieldbook_DF) %in% c("PLOTNAME", "LOC", "EXPT", "SYNS"))]
  }
  ######Code from Jeanette
  workbook_name <- paste0(trialName, "_pheno_book", format(Sys.Date(),  "%b%d%y"), ".xlsx") #generate a workbook name or add a name here
  
  wb <- createWorkbook() #create a new workbook
  
  leftAlignedColNames <- c("PLOTNAME", "TRIAL", "LOC", "EXPT", "DESIG", "PED", "SYNS")
  leftAlignedColInds <- which(colnames(fieldbook_DF) %in% leftAlignedColNames)
  centAlignedColInds <- which(!colnames(fieldbook_DF) %in% leftAlignedColNames)
  
  # Styles for the data table 
  #title_style <- createStyle(textDecoration = "bold")
  colnames_style <- createStyle(wrapText = TRUE, textDecoration = "bold", fgFill = "#d9ead3",textRotation = 45,
                                border = "TopBottomLeftRight", borderStyle = "dotted", borderColour = "black", fontSize = 10, halign = "center") 
  body_style <- createStyle(halign = "left", border = "TopBottomLeftRight", borderStyle = "dotted", borderColour = "grey", fontSize = 9)
  
  addWorksheet(wb=wb, sheetName = "RawData")
  writeData(wb=wb, headerStyle = colnames_style , sheet = "RawData", fieldbook_DF, keepNA = T, na.string = "NA")
  
  
  addStyle(wb, "RawData", style = body_style, rows = 2:(nrow(fieldbook_DF)+1) , cols = 1:(ncol(fieldbook_DF) + extraCols), gridExpand = TRUE)
  addStyle(wb, "RawData", style = colnames_style, rows = 1 , cols = (ncol(fieldbook_DF)):(ncol(fieldbook_DF) + extraCols), gridExpand = TRUE)
  
  #Fix to left or right aligned
  addStyle(wb, "RawData", style = createStyle(halign = "left"), 
           rows = 2:(nrow(fieldbook_DF)+1) , cols = leftAlignedColInds, gridExpand = TRUE, stack = TRUE)
  addStyle(wb, "RawData", style = createStyle(halign = "center"), 
           rows = 2:(nrow(fieldbook_DF)+1), cols = centAlignedColInds, gridExpand = TRUE, stack = TRUE)
  
  #Adjust PED and DESIG columns
  addStyle(wb, "RawData", style = createStyle(fontSize = 7), rows = 2:(nrow(fieldbook_DF)+1), 
           cols = which(colnames(fieldbook_DF) %in% c("PED", "SYNS")), gridExpand = TRUE, stack = TRUE)
  addStyle(wb, "RawData", style = createStyle(textDecoration = "bold", fontColour = "#073763"), rows = 2:(nrow(fieldbook_DF)+1) ,
           cols = which(colnames(fieldbook_DF) == "DESIG"), gridExpand = TRUE, stack = TRUE)
  
  setRowHeights(wb, "RawData", rows = 1, heights = 55)
  
  for (colInd in c(1:ncol(fieldbook_DF))) {
    print(colInd)
    charWidths <- nchar(fieldbook_DF[, colInd])
    charWidths <- charWidths[!is.na(charWidths)]
    charWidths <- charWidths[order(charWidths)]
    colWidth <- charWidths[ceiling(length(charWidths) * 0.9)] #Don't want the longest. Like there can be some really long ones.
    if (colWidth < 4.5) {
      colWidth <- 4.5
    } else if (colnames(fieldbook_DF)[colInd] == "PED") {
      colWidth <- 15
    } else if (colnames(fieldbook_DF)[colInd] == "SYNS") {
      colWidth <- 5
    }
    setColWidths(wb, "RawData", cols = colInd, widths = colWidth + 1) #A bit of spacing.
  }
  
  
  #save the workbook
  saveWorkbook(wb, workbook_name, overwrite = T) 
}

#This relies on your phenotype column names following the SunGrains standards
jam_fieldbook_into_template <- function(fieldbook_DF, phenoSummary_DF, templatePath) {
  filledTemplateName <- paste0(gsub("\\.xl.*$", "", templatePath), "_pheno_book", format(Sys.Date(),  "%b%d%y"), "_", fieldbook_DF$LOC[1], ".xlsx")
  
  wb <- loadWorkbook(templatePath)
  
  ###Fill raw data first
  locData <- read.xlsx(wb, sheet = "Location Data Form")
  
  #Need LOC RANGE	ROW	REP	PLOT 
  #Plotname vs plot
  colnames(fieldbook_DF)[which(colnames(fieldbook_DF) == "PLOTNAME")] <- "PLOT " #Ask jeanette to fix
  
  #Fragile part here -- maybe look for phenotypes in reference list to figure cols?
  phenoNameVec <- colnames(fieldbook_DF)[c((which(colnames(fieldbook_DF) == "PED") + 1):ncol(fieldbook_DF))]
  phenoNameVec <- phenoNameVec[phenoNameVec != "SYNS"]
  
  newDataIndex <- which(locData == "NOTES", arr.ind = T)
  for (colName in c("LOC", "RANGE", "ROW", "PLOT ", "DESIG", "REP", "DESIG", "PED", phenoNameVec)) {
    phenoVec <- fieldbook_DF[, colName]
    headerIndex <- which(locData == colName, arr.ind = T)
    
    #Will probably have to extend this for converting julian DOY
    #If anything that's not NA already becomes NA after as.numeric, probly not numeric data
    #UPDATE - should be obsoleted by conversion in get_phenotype
    #phenoVec <- convert_trait_if_numeric(phenoVec)
    
    #Begins writing one row before specified
    if (length(headerIndex) == 0) {
      #Append to end. Include column name so don't add 1 to row posn.
      newDataIndex[2] <- newDataIndex[2] + 1
      writeData(wb, sheet = "Location Data Form", x = c(colName, phenoVec),
                startRow = newDataIndex[1] + 1, startCol = newDataIndex[2],
                colNames = T)
    } else {
      writeData(wb, sheet = "Location Data Form", x = phenoVec, 
                startRow = headerIndex[1] + 2, startCol = headerIndex[2],
                colNames = F)
    }
  }
  
  if (!is.null(phenoSummary_DF)) {
    ### Fill phenotype summary data
    sumData <- read.xlsx(wb, sheet = "Summary Data Form")
    
    newDataIndex_Sum <- which(sumData == "NOTES", arr.ind = T)
    
    ##Read in genotype order so we can order our genotypes correctly.
    germplasmIndex <- which(sumData == "ID", arr.ind = T)
    
    #Add 1 to avoid reading header
    germOrder <- sumData[(1 + germplasmIndex[1]):(nrow(sumData)-4), germplasmIndex[2]]
    
    ####resort fieldbook_DF here
    config <- list(database_terms = list(name = TRUE, synonyms = TRUE, accession_numbers = TRUE),
                   search_routines = list(name = TRUE, punctuation = TRUE, substring = TRUE, edit_distance = FALSE, max_edit_distance = 2),
                   return_records = FALSE)
    accDB <- createAccessionSearchDB("LSU DB", "http://sgbreedbase.agcenter.lsu.edu/brapi/v1", "v1.3")
    updateAccessionSearchCache(accDB)
    
    germMatched <- matchAccessions(accDB, germOrder, config)
    if (any(is.na(germMatched))) {
      print("Error: Not all germplasm matched")
      break
    } else if (length(germMatched[!germMatched %in% phenoSummary_DF$DESIG]) > 0) {
      #This will also trigger in cases where line is missing in phenotype data.
      print("Check duplicate accessions for:")
      print(germMatched[!germMatched %in% phenoSummary_DF$DESIG])
    }
    
    #Make sure NA rows are PRESERVED
    phenoSummaryEntries_DF <- phenoSummary_DF[match(c(germMatched, "Mu", "CV", "LSD", "H"), phenoSummary_DF$DESIG), ]
    
    #REFACTOR THIS WITH ABOVE FOR LOOP FOR LOC DATA
    for (colName in c("DESIG", phenoNameVec)) {
      phenoVec <- phenoSummaryEntries_DF[, colName]
      headerIndex <- which(sumData == colName, arr.ind = T)
      
      #Will probably have to extend this for converting julian DOY
      #If anything that's not NA already becomes NA after as.numeric, probly not numeric data
      #Convert to character first as may be factor...
      if (any(!is.na(suppressWarnings(as.numeric(as.character(phenoVec[which(!is.na(phenoVec))])))))) {
        phenoVec <- as.numeric(phenoVec)
        print(colName)
      } 
      
      #Begins writing one row before specified
      if (length(headerIndex) == 0) {
        #Append to end. Include column name so don't add 1 to row posn.
        newDataIndex_Sum[2] <- newDataIndex_Sum[2] + 1
        writeData(wb, sheet = "Summary Data Form", x = c(colName, phenoVec),
                  startRow = newDataIndex_Sum[1] + 1, startCol = newDataIndex_Sum[2],
                  colNames = T)
      } else {
        writeData(wb, sheet = "Summary Data Form", x = phenoVec, 
                  startRow = headerIndex[1] + 2, startCol = headerIndex[2],
                  colNames = F)
      }
    }
  }
  
  saveWorkbook(wb, filledTemplateName, overwrite = T)
  
}


calculate_stats <- function(traitDF, phenoName, testBlock = TRUE) {

  #NEED TO CLEAN THIS UP -- DONT NEED TO BE RUNNING ANOVA HERE
  traitDF <- traitDF[, c("DESIG", "REP", phenoName)]
  traitDF$REP <- as.factor(traitDF$REP)
  traitDF$DESIG <- as.factor(traitDF$DESIG)
  colnames(traitDF)[3] <- "trait"
  traitDF$trait <- as.numeric(traitDF$trait)
  traitMod <- lm(data = traitDF, trait ~ DESIG + REP)
  traitAnova <- anova(traitMod)
  
  traitRan <- mmer(data = data.frame(traitDF), fixed = trait ~ REP,
                   random = ~ DESIG)
  
  #Wait, have to run this in reverse... 
  if (testBlock) {
    if ((traitAnova$`Pr(>F)`[2] > .05)) { #Does block have *any* effect? If not, don't include it.
      traitMod <- lm(data = traitDF, trait ~ DESIG) 
      traitAnova <- anova(traitMod)
      traitRan <- mmer(data = traitDF, fixed = trait ~ 1,
                       random = ~ DESIG)
    }
  }
  
  traitMMeans <- as.data.frame(emmeans(traitMod, specs = "DESIG"))
  traitMMeans <- traitMMeans[,c(1:2)]
  traitMMeans[,2] <- round(traitMMeans[,2], 2)
  colnames(traitMMeans)[2] <- phenoName
  
  ranSum <- summary(traitRan)
  modH <- ranSum$varcomp["DESIG.trait-trait", "VarComp"] / ( ranSum$varcomp["DESIG.trait-trait", "VarComp"] +  ranSum$varcomp["units.trait-trait", "VarComp"])
  #modR2 <- meanReliabilityBLUPs(traitRan, "Name")
  
  modDf <- traitAnova$Df[length(traitAnova$Df)]

  geno_emmeans <- emmeans(m2, specs = "DESIG")
  avg_se_genoMeans <- mean(data.frame(geno_emmeans)$SE)
  
  #Alpha /2 for two-sided test
  modT <- qt(.05/2, df = modDf, lower.tail = F)
  modLSD <- modT * sqrt(2) * avg_se_genoMeans
  
  modMu <- mean(traitDF$trait, na.rm = T)
  modCV <- 100 * (sd(traitDF$trait, na.rm = T) / mean(traitDF$trait, na.rm = T))
  
  statsVec <- c(round(c(modMu, modCV, modLSD, modH), 2))
  names(statsVec) <- c("Mu", "CV", "LSD", "H")
  
  statsDF <- data.frame(DESIG = names(statsVec), statsVec)
  colnames(statsDF)[2] <- phenoName
  
  return(rbind(traitMMeans, statsDF))
  
}

get_lsmeans <- function(traitDF, testBlock = TRUE) {
  #Fragile piece here
  phenoNames <- colnames(traitDF)[(which(colnames(traitDF) == "SYNS")+1):ncol(traitDF)]
  
  statFrame <- NULL
  for (phenoName in phenoNames) {
    newStatFrame <- calculate_stats(traitDF, phenoName, testBlock)
    if (is.null(statFrame)) {
      statFrame <- newStatFrame
    } else {
      statFrame <- merge(statFrame, newStatFrame, by = "DESIG", sort = F)
    }
  }
  
  return(statFrame)
}
