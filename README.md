[![Project Status: Active - The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)

# The `hdxmsqc` package

The [**hdxmsqc**] (Hydrogen-Deuterium Exchange Mass spectrometry quality control)
package is an R/Bioconductor package for HDX-MS quality control. The package
allows the enconding of experimental design, dedicated storage, quality
control visualisation and metrics. Whilst many tools metrics are difficult
to understand, unintpretable or rely on using a particular preprocessing software.
`hdxmsqc` can examine any data that can be formatted correctly into a `QFeatures`
object.

# Installation requirements

Users will require a working version of R, currently at least version >4.
It is recommended to use [RStudio](https://www.rstudio.com). The package can 
then be installed using the `devtools` package. The package should take a few 
minutes to install on a regular desktop or laptop. The package will need to be
loaded using `library(hdxmsqc)`

# Basic ideas and concepts

The [`hdxmsqc`] package implements quality control metrics for hdx-ms data.
Data from form such experiments most commonly yield a matrix of measurements 
where we have peptides along the rows, and samples (timepoints/conditions)
along the columns. Assessing the quality of `hdxmsqc` data can be challenging
and laborious. This package seeks to streamline this analysis and provide
a uniform way to present quality control of data.


To use `hdxmsqc` the data must be stored as a `QFeatures`/`Spectra`, as
implemented in the Bioconductor 
[`RForMassSpectrometry`](https://www.rformassspectrometry.org/) packages. There
are functions within the package to coherse the data into the correct format.

# Vignettes

There is currently one vignette in the package which runs through
the analysis. Analysis is designed to be quick. Though assessing the quality
of Spectra can take longer.

# Documentation

Functions are documented and be founded in the `Help`

# Contribution

Contributions are welcome, please open an issue so ew can discuss any contribution
in advance

# Feature requests

This package is actively being developed and maintained, please open Github
issue if you would like to request or discuss a particular feature.

# Licensing

The package is free to use for academic and commerical users at no cost.
Commerical users who need help installing on internal systems or more elaborate
requests designed for their infrastructure may contact us.



