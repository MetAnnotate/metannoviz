metannotate_data_filepath <- "../inputs/rpoB_0_MetagenomeTest_0_annotations_5z4KAl762541689.tsv"
hmm_info_filepath <- "../inputs/hmm_info_template_FILLED.tsv"
dataset_info_filepath <- "../inputs/dataset_info_template_FILLED.tsv"

metannotate_data_mapped <- read_metannotate_data(metannotate_data_filepath) %>%
  map_naming_information(hmm_info_filepath, dataset_info_filepath)

saveRDS(metannotate_data_mapped, file = "metannotate_data_mapped.rds")