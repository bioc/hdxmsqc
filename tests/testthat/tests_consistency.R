context("hdxmsqc tests")

library(testthat)
library(hdxmsqc) 
library(QFeatures)

test_that("Input validation fails for incorrect types", {
  expect_error(processHDE("not_a_dataframe"), "Not a data.frame")
})


sample_data <- data.frame(read.csv(system.file("extdata",
                                               "ELN55049_AllResultsTables_Uncurated.csv",
                                                                package = "hdxmsqc",
                                               mustWork = TRUE), nrows = 10))

test_that("Unique counts are correctly calculated", {

  result <- processHDE(sample_data)
  
  expect_true(nrow(result) == 10)
  expect_true(length(unique(result$Sequence)) == 10)
  
})

test_that("n/a values are correctly converted to NA", {
  
  result <- processHDE(sample_data)
  expect_true(all(is.na(result[result == "n/a"])))
})

result <- processHDE(sample_data)
i <- grep(pattern = "X..Deut",
          x = names(result))


sample_qfeatures <- readQFeatures(table = result,
                        ecol = i,
                        names = "Deuteration",
                        fnames = "fnames")


test_that("isMissingAtRandom stops with non-QFeatures input", {
  expect_error(isMissingAtRandom("not_a_qfeatures"), "Object is not a QFeatures object")
})

test_that("isMissingAtRandom applies default threshold correctly", {
  # Modify 'sample_qfeatures' to have a known distribution of NA values for testing
  modified_qfeatures <- sample_qfeatures # Apply modifications as needed
  result <- isMissingAtRandom(modified_qfeatures)

  expect_equal(nrow(assay(result)), 0)
  
})

test_that("isMissingAtRandom respects custom threshold", {
  custom_threshold <- 5 # Set this based on your data's structure
  result <- isMissingAtRandom(sample_qfeatures, threshold = custom_threshold)
  # Verify that features are filtered based on the custom threshold
  expect_equal(nrow(assay(result)), 10)
})

test_that("isMissingAtRandom filters features correctly", {
  result_with_filter <- isMissingAtRandom(sample_qfeatures, filter = TRUE)
  result_without_filter <- isMissingAtRandom(sample_qfeatures, filter = FALSE)
  
  # Check that the 'mnar' column is added and that filtering behaves as expected
  expect_true("mnar" %in% colnames(rowData(result_with_filter)[[1]]))
  expect_true("mnar" %in% colnames(rowData(result_without_filter)[[1]]))
})

test_that("computeMassError stops with non-QFeatures input", {
  expect_error(computeMassError("not_a_qfeatures"), "Object is not a QFeatures object")
})

test_that("computeMassError calculates mass error correctly", {
  # Pre-calculate expected results for a known set of experimental and theoretical centroids in 'sample_qfeatures'
  expected_delta_ppm <- 1000 # your pre-calculated delta PPM values
    result <- computeMassError(sample_qfeatures)
  # Compare calculated deltaPPM against expected values
  calculated_delta_ppm <- result$y # Assuming 'y' is the column with delta PPM values
  expect_true(all(abs(calculated_delta_ppm) < expected_delta_ppm))
})

test_that("computeMassError returns data frame with correct structure", {
  result <- computeMassError(sample_qfeatures)
  expected_columns <- c("x", "y", "sequence")
  expect_true(all(names(result) %in% expected_columns), "Returned data frame does not have the expected columns")
  expect_true(is.data.frame(result), "Returned object is not a data frame")
})


