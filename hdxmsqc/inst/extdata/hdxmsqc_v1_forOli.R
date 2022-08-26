#hdxmsqc v1
#Nathan Gittens, 20th July 2022
#hdxmsqc v0.99.0 - Author: Oliver Crook
#######################################

####Raw Data Files for PROTEIN####
hdexaminer_file <- "Olly_multiplechargestates_v2.csv"
protein_states <- c("apo_PROTEIN","mAb_1","mAb_2","mAb_3","mAb_4","mAb_5","mAb_6")
hdx_timepoints <- c(0, 15, 150, 1500)

####################################################################


#Overwrite default options based on values found in "hdxstats_parameters.txt". 
#if(file.exists("hdxstats_parameters_Bcl.txt")){
#  source("hdxstats_parameters_Bcl.txt")
#} else {
#  stop("ERROR! hdxstats_parameters.txt not found")
#}

# packages

#The packages required are the following.
require(hdxmsqc)
require(S4Vectors)
require(dplyr)
require(tidyr)
require(QFeatures)
require(RColorBrewer)
require(ggplot2)
require(MASS)
require(pheatmap)
require(Spectra)
require(patchwork)
# Data

#We first load the data, as exported from HDExaminer. 
HDXuncurated <- data.frame(read.csv(hdexaminer_file))

#The following code chunk tidies dataset, which improves the formatting and converts
#to wide format. It will also note the number of states, timepoints and peptides.

HDXuncurated_wide <- processHDE(HDExaminerFile = HDXuncurated,
                                 proteinStates = protein_states)

#The next code chunk extracts the columns with the quantitative data. 

i <- grep(pattern = "X..Deut",
          x = names(HDXuncurated_wide))

#We now parse the object into an object of class `Qfeatures`. This standardises
#the formatting of the data.

HDXdf <- readQFeatures(table = HDXuncurated_wide,
                        ecol = i,
                        names = "Deuteration",
                        fnames = "fnames")

# Visualisation

#A simple heatmap of our data can give us a sense of it.
pheatmap(assay(HDXdf), cluster_rows = FALSE, scale = "row")

# Examining missing values
#Here, we can plot where the missing values are:
plotMissing(object = HDXdf)

#Here, we can filter data that is not missing at random:
HDXdf_filtered <- isMissingAtRandom(object = HDXdf)

#We can then replot missing-ness:
#plotMissing(object = HDXdf_filtered)

#The values that are missing are all at the zero time-points where deuterium
#uptake should be 0, we can simply impute these values.
HDXdf_filtered_imputed <- impute(HDXdf_filtered, method = "zero")

# Empirical vs Theoretical errors
massError <- computeMassError(object = HDXdf_filtered_imputed, )
plotMassError(object = HDXdf_filtered_imputed)

# Intensity based outlier detection

#Using linear-model based outlier detection we see whether there
#are Spectra that have variable intensity based on their mean intensity. A linear
#model is fitted to the log-mean and log-variance of the intensities. These
#should follow a linear trend. Cook distance is used to determine outliers are
#consider if their distance is greater than 2/$\sqrt(n)$, where $n$ is the 
#number of peptides.
intensityOutlier <- intensityOutliers(object = HDXdf_filtered_imputed)
plotIntensityOutliers(object = HDXdf_filtered_imputed)


# Retention time analysis

#Retention time outlier detection looks at the
#usual variability of retention time search window and the
#actual left/right windows of the retention time. Outliers are flagged
#if their retention time falls outside 1.5 * interquartile range.
dfrt <- rTimeOutliers(object = HDXdf_filtered_imputed)
plotrTimeOutliers(object = HDXdf_filtered_imputed)

# Monotonicity statistics 

#This uses a statistic to detect differences from monotonic behavior. First,
#we need to specify the experimental design and the timepoints used. 
experiment <- protein_states
timepoints <- hdx_timepoints

#The monotonicity statistic measure the deviation from monotoncity. Whilst
#some deviation is expected from random fluctuations, it is worth double
#checking those that are strong deviates compare to the rest of the data.
monoStat <- computeMonotoneStats(object = HDXdf_filtered_imputed,
                                 experiment = experiment, 
                                 timepoints = timepoints)
out <- plotMonotoneStat(object = HDXdf_filtered_imputed,
                                 experiment = experiment, 
                                 timepoints = timepoints)
out

# Ion Mobility Time analysis

#In a similar analysis to the retention time analysis, for ion mobility time
#we can also see whether there are random deviation in the ion mobility windows.
#Again, we define outliers that deviate outside the typical 1.5 * IQR.

imTimeOut <- imTimeOutlier(object = HDXdf_filtered_imputed)
plotImTimeOutlier(object = HDXdf_filtered_imputed)

# Charge state correlation

#We check that charge states are correlated. Whilst we do not expect exactly
#the same before - low correlation maybe concerning.

csCor <- chargeCorrelationHdx(object = HDXdf_filtered_imputed,
                              experiment = experiment,
                              timepoints = timepoints)
csCor

# Using sequence overlap information are uptake values compatible

#We can also check whether uptakes are compatible with overlapping peptides.
#The difference in uptake cannot be more different than the difference
#in the number of exchangeable amides. The default methodology only checks
#whether sequence with up-to 5 different exchangeable amides are compatible
#to keep run-times lower. Larger difference may indicate different 
#chemical changes or back-exchange properties. 

tocheck <- compatibleUptake(object = HDXdf_filtered_imputed,
                            experiment = experiment,
                            timepoints = timepoints)

# Comparison of Spectra

#In this section, we can directly examine the differences between the 
#theoretical spectra one would expect from the computed deuterium uptake and
#the actual observed spectra. Deviations observed in the spectra could 
#suggest contamination, false identifications or poor quality spectra.
#A score is generated using the cosine similarity between the spectra - which
#is equivalent to the normalized dot product. The spectra pairs can be 
#also be visualized. 


#Load in some Spectra from HDsite which should match those of HDExaminer

hdxsite <- data.frame(read.csv(system.file("extdata", "HDX_RowChecked_20220628_HDsite.csv",
                                           package = "hdxmsqc", mustWork = TRUE),
                               header = TRUE, fileEncoding = 'UTF-8-BOM'))
HDXmatched <- read.csv(system.file("extdata", "HDX_RowChecked_20220628_HDE.csv",
                                    package = "hdxmsqc", mustWork = TRUE),
                        header = TRUE, fileEncoding = 'UTF-8-BOM')

hdxsite <- data.frame(read.csv("BRD4_RowChecked_20220628_HDsite.csv",
                               header = TRUE, fileEncoding = 'UTF-8-BOM'))
HDXmatched <- read.csv("BRD4_RowChecked_20220628_HDE.csv",
                       header = TRUE, fileEncoding = 'UTF-8-BOM')

spectraCompare <- spectraSimilarity(peaks = hdxsite,
                                    object = HDXmatched, 
                                    experiment = experiment,
                                    numSpectra = NULL)

#The scores can be accesses as follows:

head(spectraCompare$observedSpectra$score)

#To visualise these spectra we can use the following function

plotSpectraMirror(spectraCompare$observedSpectra[1,], spectraCompare$matchedSpectra[1,], ppm = 300)

qctable <- qualityControl(object = HDXdf_filtered_imputed, 
                          massError = massError,
                          intensityOutlier = intensityOutlier,
                          retentionOutlier = dfrt,
                          monotonicityStat = monoStat,
                          mobilityOutlier = imTimeOut,
                          chargeCorrelation = csCor,
                          sequenceCheck = tocheck,
                          spectraCheck = spectraCompare,
                          experiment = experiment,
                          timepoints = timepoints )

