# Testing data files
metannotate_data_filepath <- "../test-data/inputs/rpoB_annotations.tsv"
hmm_info_filepath <- "../test-data/inputs/hmm_info.tsv"
dataset_info_filepath <- "../test-data/inputs/dataset_info.tsv"

test_read_1_md5 <- read_metannotate_data(metannotate_data_filepath) %>%
  digest::digest("md5")

test_map_1_md5 <- read_metannotate_data(metannotate_data_filepath) %>%
  map_naming_information(hmm_info_filepath, dataset_info_filepath) %>%
  digest::digest("md5")

test_that("load functions work end to-end", {
  expect_identical(test_read_1_md5, "7ece44a61dfa500f015c6f3788154c24")
  expect_identical(test_map_1_md5, "237d22808bd5b7612a53a5b62c8a2f8f")
})
