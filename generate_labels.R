#In this case, we assume we walk along "X"
segmentPlotStakes <- function(trialDesign_DF, colMax = 25) {
  
  #If we're too long, want multiple stakes per X
  subStart_index <- 1
  stakeIndices <- c()
  stakeIndices_sub <- c(subStart_index)
  
  #Keep track of when you increase and decrease for next part...
  #Iterate through plots, adding stake at change in y. At end of each subtrial, add in the halfway stakes if greater than colmax.
  for (plot_index in c(2:nrow(trialDesign_DF))) {
    yChange <- as.numeric(trialDesign_DF$yCoord[plot_index]) - as.numeric(trialDesign_DF$yCoord[plot_index - 1])
    
    if (yChange == 1) {
      stakeIndices_sub <- c(stakeIndices_sub, plot_index)
      
      #This will happen when the trial is composed of two sub-trials
    } else if (yChange < 0 | plot_index == nrow(trialDesign_DF)) {
      #Have to subtract by one, otherwise will get first plot of next trial.
      if (plot_index == nrow(trialDesign_DF)) {
        subEnd_index <- plot_index
      } else {
        subEnd_index <- plot_index - 1
      }
      trialDesign_DF_sub <- trialDesign_DF[subStart_index:subEnd_index, ]
      
      #1 make it inclusive
      trialWidth <- 1 + max(as.numeric(trialDesign_DF_sub$xCoord)) - min(as.numeric(trialDesign_DF_sub$xCoord))
      
      if (trialWidth > colMax) {
        halfPoint <-  floor(trialWidth / 2)
        
        #Each stakeIndex gets a new halfway stake index buddy
        for (plot_index_sub in stakeIndices_sub) {
          plot_x <- as.numeric(trialDesign_DF_sub$xCoord[plot_index_sub - subStart_index + 1])
          #Adding 1 going forward and subtracking backwards gets us to our intuition of "halfway"
          if (plot_x == max(as.numeric(trialDesign_DF_sub$xCoord))) {
            new_x <- (plot_x + 1) - halfPoint
          } else if (plot_x == min(as.numeric(trialDesign_DF_sub$xCoord))) {
            new_x <- (plot_x - 1) + halfPoint
          } else {
            print("Error: Check that trial edges are regular")
            # break
          }
          
          newStakeIndex <- which(as.numeric(trialDesign_DF_sub$xCoord == new_x) & 
                                   as.numeric(trialDesign_DF_sub$yCoord) == as.numeric(trialDesign_DF_sub$yCoord[plot_index_sub - subStart_index + 1]))
          
          stakeIndices_sub <- c(stakeIndices_sub, newStakeIndex - 1 + subStart_index)
        }
        stakeIndices_sub <- stakeIndices_sub[order(stakeIndices_sub)]
      }
      #Go from sub-df index to overall index
      stakeIndices <- c(stakeIndices, stakeIndices_sub)
      #Restart everything at new starting plot
      subStart_index <- plot_index
      stakeIndices_sub <- c(subStart_index)
    }
  }
  
  return(trialDesign_DF[stakeIndices, ])
}

#Makes stake labels with Fieldbook-readable barcodes
makeStakeLabels <- function(trialStakes_DF, stakeOrder = "columnWise", includeBarcode = TRUE, labelSize = "1x3") {
  stake_labels <- ""
  #if stakeOrder == plotWise, do nothing. We sorted previously.
  if (stakeOrder == "columnWise") {
    
  } else if (stakeOrder %in% c("plotWise", "columnWise")) {
    return("Error: stake order option not recognized. Choose from 'plotWise', 'columnWise'")
  }
  if (labelSize == "1x3") {
    for (plot_index in c(1:nrow(trialStakes_DF))) {
      trialString <- gsub("\\w{3}$", "", trialStakes_DF$studyName[plot_index])
      locString <- gsub("^.*\\d{2}", "", trialStakes_DF$studyName[plot_index])
      
      #Technically I should be able to pull this from breedbase
      locNames <- c(LAB = "Baton Rouge", LAW = "Winnsboro", LAA = "Alexandria", MSS = "Stoneville", ARM = "Marianna")
      if (locString %in% names(locNames)) {
        locString <- locNames[locString]
      }
      
      stake_label_string <- createStakeLabelString(plotIdString = trialStakes_DF$observationUnitDbId[plot_index], 
                                                   trialString = trialString,
                                                   plotNumString = trialStakes_DF$plotNum[plot_index],
                                                   accessionString = trialStakes_DF$germplasmName[plot_index],
                                                   xString = trialStakes_DF$xCoord[plot_index],
                                                   yString = trialStakes_DF$yCoord[plot_index],
                                                   locString = locString,
                                                   includeBarcode = TRUE, labelSize = "1x3")
      
      stake_labels <- paste0(stake_labels, "\n", stake_label_string)
      
    }
  }
  return(stake_labels)
  
}

createStakeLabelString <- function(plotIdString, trialString, plotNumString, accessionString,
                                   xString, yString, locString, includeBarcode = TRUE, labelSize = "1x3") {
  paste0(
    "^XA\n",
    "^FWR,0\n",
    "^CF0,45\n",
    "^FO540,20^FD", trialString, "^FS\n",
    "^FO490,20^FD", plotNumString, "^FS\n",
    "^FO485,20^GB2,160,2^FS\n",
    "^CF0,30\n",
    "^TBR,195, 100^FO370,20^FD", accessionString, "^FS\n",
    "^FO350,20^GB2,160,2^FS\n",
    "^FO300,20^FDX: ", xString, "^FS\n",
    "^FO260,20^FDY: ", yString, "^FS\n",
    "^FO210,20^FD", locString, "^FS\n",
    "^FO30,10^BQN,2,7,Q,7^FDQA,", plotIdString,"^FS\n",
    "^XZ") #R write function includes terminal newline
}


createStockLabelString <- function(stockString, accessionString, dateString, defSize = "20") {
  paste0(
    "^XA\n",
    "^FO15,10^BQN,2,7,Q,7^FDQA,", stockString, "^FS\n",
    "^CF0,", defSize, "\n",
    "^TBN,200,80^FO200,20^FD", accessionString, "^FS\n",
    "^FO200,95^GB200,2,2^FS\n",
    "^TBN,200,80^FO200,105^FD", stockString, "^FS\n",
    "^FO200,163^GB200,2,2^FS\n",
    "^CF0,16\n",
    "^FO200,173^FD", dateString, "^FS\n",
    "^XZ") #R write function includes terminal newline
}