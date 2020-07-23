# Load the mapped metannotate data generated from test_map_1_md5 (before the md5)
metannotate_data_mapped <- readRDS("../test-data/intermediates/metannotate_data_mapped.rds")

# TODO - maybe can eliminate this one longer term
test_norm_1_md5 <- explore(metannotate_data_mapped,
                           dump_raw_data = TRUE, quietly = TRUE) %>%
  digest::digest("md5")

test_plot_1 <- explore(metannotate_data_mapped,
                       plot_type = "bar", quietly = TRUE)
test_plot_1_md5 <- test_plot_1$data %>%
  digest::digest("md5")

test_that("explore works end-to-end", {
  expect_identical(test_norm_1_md5, "ce37837e31bbad812e80aef31ad02eca")
  expect_identical(test_plot_1_md5, "392eec495df3f0d484ba0922091745cd")
})
