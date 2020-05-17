# Code to generate metannotate_data_mapped.rds
metannotate_data_filepath <- "../inputs/rpoB_annotations.tsv"
hmm_info_filepath <- "../inputs/hmm_info.tsv"
dataset_info_filepath <- "../inputs/dataset_info.tsv"

metannotate_data_mapped <- read_metannotate_data(metannotate_data_filepath) %>%
  map_naming_information(hmm_info_filepath, dataset_info_filepath)

saveRDS(metannotate_data_mapped, file = "metannotate_data_mapped.rds")