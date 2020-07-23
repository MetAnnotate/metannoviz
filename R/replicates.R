# replicates.R
# Combine replicates
# Copyright Jackson M. Tsuji, 2020 (Neufeld Lab)

#' Combine replicates
#'
#' @aliases dereplicate
#' @description: Combines replicates in a normalized metannotate table into means with standard deviations
#' @param metannotate_data_normalized_list List output of \code{\link{normalize}}
#' @return List of two:
#' - Tibble of metannotate data, with 'percent_abundance' and 'percent_abundance_sd'
#' - Tibble summarizing total normalized counts for genes
#' @export
combine_replicates <- function(metannotate_data_normalized_list) {
  # # Example column names of the plotting table, if collapsed to family
  # [1] "Dataset"       "replicate"         "HMM.Family"                   "Closest.Homolog.Superkingdom"
  # [4] "Closest.Homolog.Phylum"       "Closest.Homolog.Class"        "Closest.Homolog.Order"
  # [7] "Closest.Homolog.Family"       "percent_abundance"

  # Check that the first input is actually a list
  if (class(metannotate_data_normalized_list)[1] != "list") {
    stop("Input object is not a list as expected. Are you sure that you have already run normalize()?")
  }

  # Extract list components
  metannotate_data <- metannotate_data_normalized_list$metannotate_data
  hit_totals <- tidyr::pivot_longer(metannotate_data_normalized_list$total_normalized_hits, -c('Dataset', 'replicate'),
                                    names_to = "HMM.Family", values_to = "percent_abundance")
  hit_totals$HMM.Family <- factor(hit_totals$HMM.Family, levels = unique(hit_totals$HMM.Family), ordered = TRUE)

  # Check metannotate data has been normalized
  # TODO - consider making a more exhaustive data check
  if ("percent_abundance" %in% colnames(metannotate_data) == FALSE) {
    stop(paste0("Provided metannotate_data table does not contain the expected 'percent_abundance' column. "),
         "Are you sure that you have already run normalize()?")
  }

  # Detect the taxonomy that the data has been collapsed to
  plotting_taxon_colname <- TAXONOMY_NAMING$metannotate_colnames[
    TAXONOMY_NAMING$metannotate_colnames %in% colnames(metannotate_data)] %>%
    tail(n = 1)
  plotting_taxon <- TAXONOMY_NAMING$taxonomy[match(plotting_taxon_colname,
                                                   TAXONOMY_NAMING$metannotate_colnames)]
  futile.logger::flog.debug(paste0("Plotting input dataframe has been collapsed to the '", plotting_taxon, "' level."))

  # Get the mean and standard deviation of replicates
  grouping_cols <- c("Dataset", "HMM.Family", TAXONOMY_NAMING$metannotate_colnames[
    TAXONOMY_NAMING$metannotate_colnames %in% colnames(metannotate_data)])
  metannotate_data <- dplyr::group_by_at(metannotate_data, grouping_cols) %>%
    dplyr::summarise(percent_abundance_mean = mean(percent_abundance), percent_abundance_sd = sd(percent_abundance)) %>%
    dplyr::rename(percent_abundance = percent_abundance_mean)

  # Also get means of the total hits
  hit_totals <- dplyr::group_by(hit_totals, Dataset, HMM.Family) %>%
    dplyr::summarise(percent_abundance_mean = mean(percent_abundance), percent_abundance_sd = sd(percent_abundance)) %>%
    dplyr::rename(percent_abundance = percent_abundance_mean) %>%
    dplyr::select(-percent_abundance_sd) %>%
    dplyr::ungroup() %>%
    tidyr::pivot_wider(names_from = HMM.Family, values_from = percent_abundance)

  # Output the list in the same format as for normalize()
  output_list <- list(metannotate_data, hit_totals)
  names(output_list) <- c("metannotate_data", "total_normalized_hits")

  return(output_list)
}
