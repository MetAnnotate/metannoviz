# Testing data files
# TODO - is there a better place to store these?
metannotate_data_filepath <- "../test-data/inputs/rpoB_0_MetagenomeTest_0_annotations_5z4KAl762541689.tsv"
hmm_info_filepath <- "../test-data/inputs/hmm_info_template_FILLED.tsv"
dataset_info_filepath <- "../test-data/inputs/dataset_info_template_FILLED.tsv"

# TODO - these functions rely on the load functions working. Longer-term, find a way to separate
# e.g., by having the mapped table pre-loaded somehow
# TODO - maybe can eliminate this one longer term
test_norm_1_md5 <- read_metannotate_data(metannotate_data_filepath) %>%
  map_naming_information(hmm_info_filepath, dataset_info_filepath) %>%
  explore(dump_raw_data = TRUE, quietly = TRUE) %>%
  digest::digest("md5")

test_plot_1 <- read_metannotate_data(metannotate_data_filepath) %>%
  map_naming_information(hmm_info_filepath, dataset_info_filepath) %>%
  explore(plot_type = "bar", quietly = TRUE)
test_plot_1_md5 <- test_plot_1$data %>%
  digest::digest("md5")

test_that("end-to-end test works", {
  expect_identical(test_norm_1_md5, "e6fe48672b27efef5489d13c781cb6ae")
  expect_identical(test_plot_1_md5, "e6fe48672b27efef5489d13c781cb6ae")
})
