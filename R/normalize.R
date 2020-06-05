# normalize.R
# Normalize metannotate data
# Copyright Jackson M. Tsuji, 2020 (Neufeld Lab)

#' Normalize metannotate data
#'
#' @aliases normalize_collapsed_metannotate_data
#' @description Double-normlizes metannotate data by HMM length and a single-copy gene marker. The input data has to
#' be pre-mapped and collapsed by taxonomy
#' @param metannotate_data Tibble of taxonomy-collapsed and mapped metannotate data -
#' see \code{\link{collapse_by_taxon}}
#' @param normalizing_HMM Character vector (length 1) giving the plotting name of an HMM.Family in your table
#' that you want to normalize all genes within each sample to. This should be a reliable single-copy taxonomic marker gene, 
#' e.g., rpoB, dnaK, and so on.
#' @return List of two:
#' - Tibble of metannotate data; 'hits' has now been changed to 'percent_abundance' and is double-normalized
#' - Tibble summarizing total normalized counts for genes
#' @export
normalize <- function(metannotate_data, normalizing_HMM = "rpoB") {

  # Check metannotate data has been mapped
  if ("HMM_length" %in% colnames(metannotate_data) == FALSE) {
    stop(paste0("Provided metannotate_data table does not contain the expected 'HMM_length' column. "),
         "Are you sure that you have already run map_naming_information()?")
  }

  # Check metannotate data has been collapsed by taxonomy
  # TODO - consider making a more exhaustive data check
  if ("hits" %in% colnames(metannotate_data) == FALSE) {
    stop(paste0("Provided metannotate_data table does not contain the expected 'hits' column. "),
         "Are you sure that you have already run collapse_by_taxonomy()?")
  }

  ### First normalize by HMM length (for between-HMM comparison)
  metannotate_data$hits <- metannotate_data$hits /
    metannotate_data$HMM_length
  
  ### Then normalize relative to the normalizing HMM (allows for comparison between samples)
  # Start by getting the total HMM-length normalized counts for the normalizing HMM
  total_hits <- summarize_total_reads_all_genes(metannotate_data, format = "wide", collapsed = TRUE)
  
  # Join the normalizing HMM column onto the main table
  total_hits_join <- dplyr::select(total_hits, Dataset, all_of(normalizing_HMM)) %>%
    dplyr::rename(normalizing_HMM_total_hits = normalizing_HMM)
  metannotate_data <- dplyr::left_join(metannotate_data, total_hits_join, by = "Dataset")
  
  # Normalize metannotate table
  metannotate_data$hits <- metannotate_data$hits /
    metannotate_data$normalizing_HMM_total_hits * 100
  
  ### Clean up the output
  metannotate_data <- dplyr::ungroup(metannotate_data) %>%
    dplyr::select(-HMM_length, -normalizing_HMM_total_hits) %>%
    dplyr::rename(percent_abundance = hits)
  
  ### Collect overall summary stats for the normalization
  total_hits_normalized <- (dplyr::select(total_hits, -Dataset) /
    dplyr::pull(total_hits, normalizing_HMM) * 100) %>%
    tibble::add_column(Dataset = dplyr::pull(total_hits, Dataset), .before = 1)
  
  output_list <- list(metannotate_data, total_hits_normalized)
  names(output_list) <- c("metannotate_data", "total_normalized_hits")

  return(output_list)
}

#' Filter normalized metannotate data by relative abundance
#'
#' @aliases subset subset_normalized_metannotate_data filter
#' @description Filters out low relative abundance taxa from normalized metannotate data based on a desired threshold
#' @param metannotate_data Tibble of normalized metannotate data - see \code{\link{normalize}}
#' @param top_x Numeric vector (length 1) giving the subsetting amount you desire.
#' If top_x >=1, the script will return the "top_x most abundant taxa" for each Dataset/HMM.Family
#' If top_x <1, the script will return "all taxa of (top_x * 100%) abundance or greater for each Dataset/HMM.Family -
#' but see below.
#' @param percent_mode If top_x <1, there are two different methods for keeping the most abundant organisms:
#' - "within_sample" -- the normalized % abundance relative to rpoB is used
#' - "within_HMM" -- the percent abundance of that taxon within the specific HMM gene hits is used.
#' You won't notice much of a different between these modes unless one of your HMMs has very few hits and you want to
#' show some of the taxa that were hit. This would be a good time to use 'within_HMM'.
#' @return List of two:
#' - Tibble of metannotate data; 'hits' has now been changed to 'percent_abundance' and is double-normalized
#' - Tibble summarizing total normalized counts for genes
#' @export
filter_by_abundance <- function(metannotate_data, top_x, percent_mode = "within_sample") {
  # Parse the top_x value to determine what the user wants
  if (is.na(top_x) == TRUE) {
    futile.logger::flog.info("Not filtering by abundance because 'top_x' was set as NA")

  } else if (top_x >= 1) {
    futile.logger::flog.info(paste0("Subsetting to the top ", top_x, " most abundant taxa per sample"))
    
    metannotate_data <- dplyr::group_by(metannotate_data, Dataset, HMM.Family) %>%
      dplyr::top_n(n = top_x, wt = percent_abundance)
    
  } else if (top_x < 1 && top_x > 0) {
    if (percent_mode == "within_sample") {
      futile.logger::flog.info(paste0("Subsetting all taxa of at least ", top_x * 100, "% normalized relative abundance"))
      
      metannotate_data <- dplyr::filter(metannotate_data, percent_abundance >= (top_x * 100))
      
    } else if (percent_mode == "within_HMM") {
      futile.logger::flog.info(paste0("Subsetting all taxa of at least ", top_x * 100, "% relative abundance per HMM"))
      
      # Re-calculate the total hits for each HMM
      # TODO - generalize the summary function so I don't have to rename so many columns here
      total_hits <- dplyr::rename(metannotate_data, hits = percent_abundance) %>%
        summarize_total_reads_all_genes(format = "long", collapsed = TRUE) %>%
        dplyr::rename(HMM.Family = gene, total_percent_abundance_within_HMM = hits)
      
      # Calculate the relative % contribution that the HMM makes to its own total counts
      metannotate_data <- dplyr::left_join(metannotate_data, total_hits,
                                                      by = c("Dataset", "HMM.Family"))
      metannotate_data$percent_abudance_within_HMM <- metannotate_data$percent_abundance /
        metannotate_data$total_percent_abundance_within_HMM * 100
      
      # Filter
      metannotate_data <- dplyr::filter(metannotate_data, percent_abudance_within_HMM >= (top_x * 100))
      
      # Clean up tibble
      metannotate_data <- dplyr::select(metannotate_data, -total_percent_abundance_within_HMM,
                                             -percent_abudance_within_HMM)
    }
  }
  
  return(metannotate_data)
}
