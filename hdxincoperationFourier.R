## fourier transform approach to computing isotopic distribution




fourierIsotope <- function(elements,
                           incorp = 0,
                           num_exch_sites = 0,
                           charge = 1,
                           isotopes = NULL){
    
    ## add in charge
    el <- c(elements, charge)
    el["H"] <- el["H"] - num_exch_sites
    el["D"] <- num_exch_sites
    el <- el[names(el) != "P"]
    
    # infilling with zeros to make matricies match
    
    ## carbon details
    cmass <- c(12.0000000, 13.0033548378, 0, 0)
    cprob <- c(0.9893, 0.0107, 0, 0)


    ## expected mass of unexchanged hydrogens
    hmass <- c(1.0078250321, 2.0141017780, 0, 0)
    hprob <- c(0.999885, 0.000115, 0, 0)

    ## expected mass of exchanged hydrogens
    dmass <- c(1.0078250321, 2.0141017780, 0, 0)
    dprob <- c(1 - incorp, incorp, 0, 0)

    ## expected mass of nitrogens
    nmass <- c(14.0030740052, 15.0001088984, 0, 0)
    nprob <- c(0.99632, 0.00368, 0, 0)


    ## expected mass of oxygens 
    omass <- c(15.9949146221, 16.99913150, 17.9991604, 0)
    oprob <- c(0.99757, 0.00038, 0.00205, 0)

    ## expected mass of sulphers
    smass <- c(31.97207069, 32.97145850, 33.96786683, 35.96708088)
    sprob <- c(0.9493, 0.0076, 0.0429, 0.0002)

    ## expected mass of charge
    chmass <- c(1.0072764522, 2.0135531981, 0, 0)
    chprob <- c(0.999885, 0.000115, 0 , 0)

    # Must be in the same order as elments
    prob <- rbind(cprob, hprob, nprob, oprob, sprob, chprob, dprob)
    mass <- rbind(cmass, hmass, nmass, omass, smass, chmass, dmass)
    

    # guess sensible isotopes to compute
    if(is.null(isotopes)){
        isotopes <- max(4,
                        0.5 * sqrt(1 + sum(el * apply(mass,1, function(x) var(x[x!=0])))))
        isotopes <- floor(isotopes)
    }
    # isotopic expansion
    isotopes <- isotopes + num_exch_sites
    
    # make to size of isotopes required
    finalProb <- matrix(0, nrow = nrow(prob), ncol = isotopes)
    finalMass <- matrix(0, nrow = nrow(prob), ncol = isotopes)
    finalProb[,1:4] <- prob
    finalMass[,1:4] <- mass
    
    # isotope distribution of elements
    # fourier transform of probability
    fa <- t(apply(finalProb, 1, fft))
    # expsumlog trick for stability
    a <- exp(colSums(el * log(fa)))
    # Compute probability distribution of isotopes, correct normalising const 
    pout <- fft(a, inverse = TRUE)/length(a) 
    # Compute weighted mass matrix
    b <- finalProb * finalMass
    # fourier transform
    fb <- t(apply(b, 1, fft))
    #Compute ratio
    fr <- fb/fa
    # mass distribution of el, correct normalising const
    bout <- fft(colSums(el * fr) * a, inverse = TRUE)/length(pout)
    m <- bout/pout 
    
    if(charge != 0 ) mz <- m/charge
    else mz <- m
    
    return(list(intensity = Re(pout), mz = Re(mz)))
}    

isotopicDistributionHDXfourier <- function(sequence,
                                           incorp = 0,
                                           charge = 1, 
                                           custom = list(code = NULL, elements = NULL)) {
    
    if(length(custom$elements != 0)) {
        custom_elements <- c(C = 0, H = 0, N = 0, O = 0, S = 0, P = 0)
        custom_elements[names(custom$elements)] <- custom$elements
    }
    
    if(charge < 0 | charge > 8) stop("charge must be between 1 and 8")
    
    # create vector for sequences and holder for element annotations
    seq_vector <- strsplit(sequence, split = "")[[1]]
    x <- c(C = 0, H = 0, N = 0, O = 0, S = 0, P = 0)
    
    # Update x with amino acid information
    for(i in seq.int(length(seq_vector))) {
        if(seq_vector[i] == "A") x <- x + c(C = 3, H = 5, N = 1, O = 1, S = 0, P = 0)
        if(seq_vector[i] == "R") x <- x + c(C = 6, H =12, N = 4, O = 1, S = 0, P = 0)
        if(seq_vector[i] == "N") x <- x + c(C = 4, H = 6, N = 2, O = 2, S = 0, P = 0)
        if(seq_vector[i] == "D") x <- x + c(C = 4, H = 5, N = 1, O = 3, S = 0, P = 0)
        if(seq_vector[i] == "C") x <- x + c(C = 3, H = 5, N = 1, O = 1, S = 1, P = 0) 
        if(seq_vector[i] == "E") x <- x + c(C = 5, H = 7, N = 1, O = 3, S = 0, P = 0)
        if(seq_vector[i] == "Q") x <- x + c(C = 5, H = 8, N = 2, O = 2, S = 0, P = 0)
        if(seq_vector[i] == "G") x <- x + c(C = 2, H = 3, N = 1, O = 1, S = 0, P = 0)
        if(seq_vector[i] == "H") x <- x + c(C = 6, H = 7, N = 3, O = 1, S = 0, P = 0)
        if(seq_vector[i] == "I") x <- x + c(C = 6, H =11, N = 1, O = 1, S = 0, P = 0)
        if(seq_vector[i] == "L") x <- x + c(C = 6, H =11, N = 1, O = 1, S = 0, P = 0)
        if(seq_vector[i] == "K") x <- x + c(C = 6, H =12, N = 2, O = 1, S = 0, P = 0)
        if(seq_vector[i] == "M") x <- x + c(C = 5, H = 9, N = 1, O = 1, S = 1, P = 0)
        if(seq_vector[i] == "F") x <- x + c(C = 9, H = 9, N = 1, O = 1, S = 0, P = 0)
        if(seq_vector[i] == "P") x <- x + c(C = 5, H = 7, N = 1, O = 1, S = 0, P = 0)
        if(seq_vector[i] == "S") x <- x + c(C = 3, H = 5, N = 1, O = 2, S = 0, P = 0)
        if(seq_vector[i] == "T") x <- x + c(C = 4, H = 7, N = 1, O = 2, S = 0, P = 0)
        if(seq_vector[i] == "W") x <- x + c(C =11, H =10, N = 2, O = 1, S = 0, P = 0)
        if(seq_vector[i] == "Y") x <- x + c(C = 9, H = 9, N = 1, O = 2, S = 0, P = 0)
        if(seq_vector[i] == "V") x <- x + c(C = 5, H = 9, N = 1, O = 1, S = 0, P = 0)
        
        if(length(custom$elements != 0))
            if(seq_vector[i] == custom$code) x <- x + custom_elements    
    }
    
    ## add N-terminal H and C-terminal OH
    elements <- x + c(C = 0, H = 2, N = 0, O = 1, S = 0, P = 0) 
    
    # get the number of exchangeable amides
    num_exch_sites <- exchangeableAmides(sequence)
    
    ## compute distribution using fourier method 
    res  <- fourierIsotope(elements = elements,
                           incorp = incorp,
                           num_exch_sites = num_exch_sites,
                           charge = charge,
                           isotopes = NULL)
    
    intensity <- as.numeric(res$intensity)
    mz <- as.numeric(res$mz)
    mz <- mz[!(intensity < 10^{-8})]
    intensity <- intensity[!(intensity < 10^{-8})]
    
    # generate spectra object
    spec <- DataFrame(
        msLevel = c(1L),
        charge = charge,
        sequence = sequence)
    
    spec$mz <- list(mz)
    spec$intensity <-  list(intensity)
    
    # construct spectra object
    sps <- Spectra(spec)
    
    return(sps)
    
}

##' Computes the number of exchangeable amides based on the sequnece
##' @title Compute exchangeable amides.
##' @param sequence The sequence of the peptide
##' @return Returns a numeric indicating the number of exchangeable amides
##' @md
##' 
##' @rdname hdx-distributions
exchangeableAmides <- function(sequence) {
    
    n <- length(sequence)
    x <- vector(mode = "numeric", length = n)
    
    for(i in 1:n) {
        seq_vector <- strsplit(as.character(sequence[i]), split = "")[[1]]
        x[i] <- length(na.omit(sub("P", NA, seq_vector))) - 2
    }
    
    return(x)	
}
