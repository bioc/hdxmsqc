#' Quality Control table function. Generate a table that collates quality
#' control metrics
#' @param object An object of class Qfeatures, with the data used for the analysis
#' @param massError The output of the `computeMassError` function
#' @param intensityOutlier The output of the `intensityOutliers` function
#' @param retentionOutlier The output of the `rTimeOutliers` function
#' @param monotonicityStat The output of the `computeMonotoneStats` function
#' @param mobilityOutlier The output of the `imTimeOutliers` function
#' @param chargeCorrelation The output of the `chargeCorrelationsHdx` function
#' @param replicateCorrelation The output of the `replicateCorrelation` function
#' @param replicateOutlier The output of the `replicateOutlier` function
#' @param sequenceCheck The output of the `compatibleUptake` function
#' @param spectraCheck The output of the `spectraSimiarity` function
#' @param experiment The experimental conditions.
#' @param timepoints The timepoints used in the analysis, must include repeat
#' for replicates
#' @param undeuterated A logical indicating whether only the undeuterated data 
#' should be exported 
#' @return An object of class `DataFrame` containing a summary of the quality
#' control results.
#' @md
#' @author Oliver Crook
#' @export 
qualityControl <- function(object, 
                           massError = NULL,
                           intensityOutlier = NULL,
                           retentionOutlier = NULL,
                           monotonicityStat = NULL,
                           mobilityOutlier = NULL,
                           chargeCorrelation = NULL,
                           replicateCorrelation = NULL,
                           replicateOutlier = NULL,
                           sequenceCheck = NULL,
                           spectraCheck = NULL,
                           experiment = NULL,
                           timepoints = NULL,
                           undeuterated = FALSE){
    
    
    QCtable <- DataFrame(sequence = massError$sequence)
    QCtable$timepoints <- rep(rep(timepoints,
                                  times = length(experiment)),
                                  times = length(unique(massError$sequence)))
    QCtable$experiment <- rep(rep(experiment,
                                  each = length(timepoints)),
                                  times = length(unique(massError$sequence)))
    QCtable$replicate <- as_tibble(QCtable) |>  group_by(sequence, timepoints, experiment) |>
        mutate(replicate = row_number()) |> pull(replicate)
    
    ## Get start and end
    i <- grep(pattern = "Start", x = rowDataNames(object)[[1]])
    j <- grep(pattern = "End", x = rowDataNames(object)[[1]])
    startend <- rowData(object[[1]])[,c(i[1], j[1])]
    QCtable$Start <- 0
    QCtable$End <-0
    QCtable[, c("Start", "End")] <- startend[QCtable$sequence,]
    
    ## Get score from software (e.g. HDExaminer)
    i <- grep(pattern = "Score", x = rowDataNames(object)[[1]])
    HDscore <- data.frame(rowData(object[[1]])[, i])
    QCtable$HDScore <- 0 
    QCtable$HDScore <- c(t(HDscore))
    
    
    
    
    if (is.null(massError)){
        QCtable$massCheck <- 0
    } else{
        QCtable$massCheck <- 1 * (abs(massError$y) > 300)
    }    

    if (is.null(intensityOutlier)){
        QCtable$intensityCheck <- 0 
    } else{
        iOut <- which(QCtable$sequence %in%
                      intensityOutlier$x[intensityOutlier$outlier == 1])
    
        QCtable$intensityCheck <- 0
        QCtable$intensityCheck[iOut] <- 1
    }

    if (is.null(retentionOutlier)){
        QCtable$retentionOutlier <- 0
    } else{
        QCtable$retentionOutlier <- apply(rbind(retentionOutlier$leftRT$outlier,
                                    retentionOutlier$rightRT$outlier), 2, max)
    }
    

    if (is.null(monotonicityStat)){
        QCtable$monotonicityStat <- 0
    } else{
            QCtable$monotonicityStat <- 0
    for (i in seq_along(experiment)){
          monoOut <- which(QCtable$sequence %in% 
                         monotonicityStat[[i]]$x[monotonicityStat[[i]]$outlier == 1])
          wh <- which(QCtable$experiment %in% experiment[[i]])
          
    n <- length(QCtable$monotonicityStat[intersect(wh, monoOut)])      
    QCtable$monotonicityStat[intersect(wh, monoOut)] <- rep(1, n)
    }
    }
    
    if (is.null(mobilityOutlier)){
        QCtable$mobilityOutlier <- 0
    } else {
            
        QCtable$mobilityOutlier <- apply(rbind(mobilityOutlier$leftIMS$outlier,
                                    mobilityOutlier$rightIMS$outlier), 2, 
                                   max)
    }

    if (is.null(chargeCorrelation)){
        QCtable$chargecorrelation <- 0
    } else{
    QCtable$chargecorrelation <- 0
    for (i in seq_along(experiment)){
        cRes <- colSums(chargeCorrelation[[i]] < 0.9, na.rm = TRUE) >= 1
        wC <- cRes[cRes]
        seq <- sapply(strsplit(QCtable$sequence, split = ""),
                      function(x) paste(x[seq.int(length(x)) - 1],
                                        sep = "", collapse = ""))
        sub1 <- QCtable$experiment %in% experiment[[i]]
        sub2 <- seq %in% names(wC)
        QCtable$chargecorrelation[sub1 & sub2] <- 1
    }
    }
    
    # Replicate-based check 1
    
    if(is.null(replicateCorrelation)){
        QCtable$replicatecorrelation <- 0
    } else{
        QCtable$replicatecorrelation <- 0
        
        for (i in seq_along(experiment)){
            
            rcRes <- as.numeric(replicateCorrelation[[i]]$outlier) == 1
            rcC <- replicateCorrelation[[i]]$x[rcRes]
            sub1 <- QCtable$experiment %in% experiment[[i]]
            sub2 <- QCtable$sequence %in% rcC
            
            QCtable$replicatecorrelation[sub1 & sub2] <- 1 
            
        }
    }
    
    # Replicated-based check 2
    if(is.null(replicateOutlier)){
        QCtable$replicateoutlier <- 0
    } else{
        QCtable$replicateoutlier <- 0
        
        for (i in seq_along(experiment)){
            
            rcRes <- as.numeric(replicateOutlier[[i]]$outlier) == 1
            rcC <- replicateOutlier[[i]]$x[rcRes]
            sub1 <- QCtable$experiment %in% experiment[[i]]
            sub2 <- QCtable$sequence %in% rcC
            
            QCtable$replicateOutlier[sub1 & sub2] <- 1 
        }
    }
    
    
    
    
    if (is.null(sequenceCheck)){
        QCtable$sequenceCheck <- 0
    } else{
            QCtable$sequenceCheck <- 0
        for (j in seq_along(sequenceCheck)){
            ex <- gsub("([0-9]+).*$","", sequenceCheck[[j]])
            time <- as.numeric(gsub(".*?([0-9]+).*", "\\1", sequenceCheck[[j]]))
        
            whseqCheck <- which((QCtable$sequence %in% names(sequenceCheck[[j]])) &
                                (QCtable$timepoints %in% time) &
                                (QCtable$experiment %in% ex))
        
            QCtable$sequenceCheck[whseqCheck] <- 1 
        }
    }



    # spit out fully deuterated samples
    spectraCheck$observedSpectra <- spectraCheck$observedSpectra[spectraCheck$observedSpectra$DeutTime != "FD",]

    # convert times seconds, currently as character:
    spectraCheck$observedSpectra$DeutTime <- vapply(strsplit(spectraCheck$observedSpectra$DeutTime, "s"),
                                       function(x) as.numeric(x),
                                       FUN.VALUE = numeric(1))
    
    if (is.null(spectraCheck)){
        QCtable$score <- 0
        
    } else {
        for (j in seq_along(experiment)){
            i <- grep(experiment[j], spectraCheck$observedSpectra$experiment)
            spectraCheck$observedSpectra$experiment[i] <- experiment[j]
        }
    
    charge <- sapply(strsplit(QCtable$sequence, split = ""),
                                    function(x) paste(x[length(x)],
                                                      sep = "", collapse = ""))
    QCtable$score <- NA
    
        for (j in seq.int(nrow(QCtable))){
            wh <- spectraCheck$observedSpectra$Sequence == seq[j] &
                spectraCheck$observedSpectra$Charge == charge[j]
            wh2 <- spectraCheck$observedSpectra$DeutTime == QCtable$timepoints[j] &
                spectraCheck$observedSpectra$experiment == QCtable$experiment[j] &
                spectraCheck$observedSpectra$replicate == QCtable$replicate[j]
            QCtable$score[j] <- max(spectraCheck$observedSpectra$score[which(wh & wh2)],
                                0) # give 0 to missing spectra
        }
    }

    
    QCtable$flagged <- rowSums(as.matrix(QCtable[, -(c(seq.int(6), which(colnames(QCtable) == "score")))]))
    
    if(isTRUE(undeuterated)){
        QCtable <- QCtable[QCtable$timepoints == 0,]
    }
    
    return(QCtable)
}