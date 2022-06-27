#' Function to curate and HDExaminer file so that in contains all the information
#' in a sensible format.
#' 
#' 
#' @param HDExaminerFile an object of class data.frame containing an HDExaminer
#' data
#' @param proteinStates a character vector indicating the protein states

processHDE <- function(HDExaminerFile, proteinStates = NULL){
    
    stopifnot("Not a data.frame"=is(HDExaminerFile, "data.frame"))
    
    numSeq <- length(unique(HDExaminerFile$Sequence))
    numTimepoints <- length(unique(HDExaminerFile$Deut.Time))
    numStates <- length(unique(HDExaminerFile$Protein.State))
    
    message("Number of peptide sequence: ", numSeq, 
            "\nNumber of timepoints: ", numTimepoints,
            "\nNumber of Protein States: ", numStates)
    
    # processing steps
    # Convert n/a s to NA
    HDExaminerFile[HDExaminerFile == "n/a"] <- NA
    
    # Making RT and IMS into reasonable left and right windows
    HDExaminerFile$leftRT <- 0
    HDExaminerFile$rightRT <- 0
    HDExaminerFile$leftIMS <- 0
    HDExaminerFile$rightIMS <- 0

    # assumes data are stored as "x1-x2"
    left_right_Rt <- t(vapply(strsplit(HDExaminerFile$Actual.RT,
                                       fixed = TRUE, split = "-"),
                              function(x) as.numeric(x[1:2]),
                              FUN.VALUE = numeric(2))) 
    left_right_ims <-  t(vapply(strsplit(HDExaminerFile$IMS.Range,
                                         fixed = TRUE, split = "-"),
                                function(x) as.numeric(x[1:2]),
                                FUN.VALUE = numeric(2))) 

    # Put in dataframe
    HDExaminerFile[, c("leftRT", "rightRT")] <- left_right_Rt
    HDExaminerFile[, c("leftIMS", "rightIMS")] <- left_right_ims
    
    # spit out fully deuterated samples
    HDExaminerFile_fd <- HDExaminerFile[HDExaminerFile$Deut.Time == "FD",]
    HDExaminerFile <- HDExaminerFile[HDExaminerFile$Deut.Time != "FD",]
    
    # convert times seconds, currently as character:
    HDExaminerFile$Deut.Time <- vapply(strsplit(HDExaminerFile$Deut.Time, "s"),
                                      function(x) as.numeric(x),
                                      FUN.VALUE = numeric(1))

    # add in repliate numbers
    HDExaminerFile <- HDExaminerFile |> 
        group_by(Deut.Time, Sequence, Protein.State, Charge) |>
        mutate(replicate = row_number())

    # remove annoying spaces in files and replace with 0
    # also convert to numeric
    HDExaminerFile$X..Deut[HDExaminerFile$Deut.. == ""] <- 0
    HDExaminerFile$Deut..[HDExaminerFile$Deut.. == ""] <- 0
    HDExaminerFile$Deut.. <- as.numeric(HDExaminerFile$Deut..)
    HDExaminerFile$X..Deut <- as.numeric(HDExaminerFile$X..Deut)

    if (is.null(proteinStates)){
        proteinStates <- paste0(rep("Condition ", numStates), seq.int(numStates))
        proteinStatesCurrent <- unique(HDExaminerFile$Protein.State)
        for (j in seq_along(proteinStatesCurrent)){
            HDExaminerFile$Protein.State[HDExaminerFile$Protein.State 
                                         == proteinStatesCurrent[j]] <- proteinStates[j]

        }
    } else{
        stopifnot("proteinStates does not match
                  number of states"=numStates==length(proteinStates))
        proteinStatesCurrent <- unique(HDExaminerFile$Protein.State)
        for (j in seq_along(proteinStatesCurrent)){
            HDExaminerFile$Protein.State[HDExaminerFile$Protein.State 
                                         == proteinStatesCurrent[j]] <- proteinStates[j]
            
        }
    }    

    # convert to wide format
    HDExaminerFile_wide <- pivot_wider(data = HDExaminerFile,
                                       id_cols = c("Sequence",
                                                   "Charge"),
                                       names_from = c("Protein.State",
                                                      "Deut.Time",
                                                      "replicate"),
                                       values_from = c("X..Deut",
                                                       "Search.RT",
                                                       "Actual.RT",
                                                       "X..Spectra",
                                                       "Search.IMS",
                                                       "IMS.Range",
                                                       "Max.Inty",
                                                       "Exp.Cent",
                                                       "Theor.Cent",
                                                       "Score",
                                                       "Confidence",
                                                       "leftRT",
                                                       "rightRT",
                                                       "leftIMS",
                                                       "rightIMS",
                                                       "X..Spectra",
                                                       "Start",
                                                       "End"))
    # make feature names
    HDExaminerFile_wide$fnames <- paste0(HDExaminerFile_wide$Sequence,
                                         HDExaminerFile_wide$Charge)
    
    return(HDExaminerFile_wide)

}

#' missing value plot
#' 
#' @param object An object of class `QFeatures`
#' 
#' 
#' 
#' 
#' 
#' 
#' 
plotMissing <- function(object, ...){
    
    stopifnot("Object is not a QFeatures object"=is(object, "QFeatures"))
    
    na_mat <- 1*is.na(assay(object))
    pheatmap(na_mat, 
             cluster_rows = FALSE,
             cluster_cols = FALSE,
             color = brewer.pal(n = 3, name = "Greys")[c(1, 3)],
             legend_breaks = c(0,1),
             legend_labels = c("Not Missing", "Missing"),
             main = "Missing value plot",
             fontsize = 12,...)
    
}
#'Missing at random versus missing not at random
#'
#'@param object An object of class `QFeatures`
#'@param threshold A threshold indicated how many missing values indicate
#'whether missingness is not at random. Default is NULL, which means leads to a
#'threshold which is half the number of columns.
#'@param filter A logial indicating whether to filter out data that is deemed
#' missing not at random
#' 
#' 
isMissingAtRandom <- function(object, threshold = NULL, filter = TRUE){
    
    stopifnot("Object is not a QFeatures object"=is(object, "QFeatures"))
    
    na_mat <- 1*is.na(assay(object))
    if (is.null(threshold)){
        threshold <- ncol(na_mat)/2
    }
    
    to_filter_missing <- 1*(rowSums(na_mat) > threshold)
    rowData(object)[[1]]$mnar <- to_filter_missing
    
    if (isTRUE(filter)){
        object <- filterFeatures(object, ~ to_filter_missing != 1) 
        message("Number of peptides filtered:", sum(to_filter_missing))   
    }

    
    return(object)
    
}
#'Empirical versus theoretical mass errors
#'
#'@param object An object of class `QFeatures`
#'@param eCentroid character string indicating column identifier for 
#'experimental centroid
#'@param tCentroid character string indicating column identifier for 
#'theoretical centroid
#'
computeMassError <- function(object,
                             eCentroid = "Exp.Cent",
                             tCentroid = "Theor.Cent"){
    
    j <- grep(pattern = eCentroid, x = rowDataNames(object)[[1]])
    k <- grep(pattern = tCentroid, x = rowDataNames(object)[[1]])
    
    deltaPPM <- ((as.matrix(rowData(object)[[1]][, j]) - 
                      as.matrix(rowData(object)[[1]][,k]))/as.matrix(rowData(object)[[1]][,k])) * 10^6

    ppmerror <- data.frame(x = c(t(as.matrix(rowData(object)[[1]][,k]))),
                           y = c(t(deltaPPM)), 
                           sequence = rep(rownames(object)[[1]],
                                          each = ncol(assay(object))))
    return(ppmerror)
}

#' Mass error plot
#' 
#'@param object An object of class `QFeatures`
#'@param eCentroid character string indicating column identifier for 
#'experimental centroid
#'@param tCentroid character string indicating column identifier for 
#'theoretical centroid
#'
plotMassError <- function(object,
                          eCentroid = "Exp.Cent",
                          tCentroid = "Theor.Cent"){
    
    ppmerror <- computeMassError(object = object,
                                 eCentroid = eCentroid,
                                 tCentroid = tCentroid)
    gg <- ppmerror |> ggplot(aes(x = x, y = y, col = sequence)) +
        geom_point(size = 2, alpha = 0.8) + 
        scale_color_manual(values = 
                               colorRampPalette(brewer.pal(n  = 11, 
                                                           name = "Set3"))(169)) +
        theme_classic() + 
        theme(legend.position = "none") + 
        xlab("Theoretical Centroid") + 
        ylab("Empirical Error")
    
    return(gg)
}    
    