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