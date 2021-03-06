# explore.R
# High-level data exploration of metannotate data
# Copyright Jackson M. Tsuji, 2020 (Neufeld Lab)

#' Explore MetAnnotate data
#'
#' @aliases explore_metannotate_data
#' @description high-level exploration function for examining MetAnnotate data
#' @param metannotate_data The mapped metannotate tibble output by \code{\link{map_naming_information}}
#' @param evalue E-value cutoff for HMM hits
#' @param taxon Character vector (length 1) giving the taxon name to collapse to
#' Can be: domain, phylum, class, order, family, genus, species (case insensitive)
#' @param normalizing_HMM Name of the normalizing HMM (e.g., 'rpoB')]; specify 'auto' to attempt auto-detection
#' @param top_x Numeric vector (length 1) giving the subsetting amount you desire.
#' If top_x >=1, the script will return the "top_x most abundant taxa" for each Dataset/HMM.Family
#' If top_x <1, the script will return "all taxa of (top_x * 100%) abundance or greater for each Dataset/HMM.Family -
#' but see below.
#' @param percent_mode If top_x <1, there are two different methods for keeping the most abundant organisms:
#' - "within_sample" -- the normalized % abundance relative to rpoB is used
#' - "within_HMM" -- the percent abundance of that taxon within the specific HMM gene hits is used.
#' You won't notice much of a different between these modes unless one of your HMMs has very few hits and you want to
#' show some of the taxa that were hit. This would be a good time to use 'within_HMM'.
#' @param colouring_template_filename Filename of the colouring template you want to load
#' If the file does not exist, then this function will write a template to that file
#' If 'NA' is entered, then the function will auto-generate colours and continue on
#' @param quietly logical (TRUE/FALSE); if TRUE, only reports warnings and errors
#' @param ... Other fine-tuned plotting options controlled by \code{\link{visualize}} and the underlying
#' \code{\link{generate_ggplot}}. Highlights include plot_type, which can be "bar" or "bubble"
#' @return A ggplot of MetAnnotate data
#' @export
explore <- function(metannotate_data, evalue = 1e-10, taxon = "Family",
                    normalizing_HMM = "rpoB", top_x = 0.02, percent_mode = "within_sample",
                    colouring_template_filename = NA, quietly = FALSE, ...) {

  if (quietly == TRUE) {
    futile.logger::flog.threshold(WARN)
  }

  # Check metannotate data has been mapped
  # TODO - consider making a more exhaustive data check
  if ("HMM_length" %in% colnames(metannotate_data) == FALSE) {
    stop(paste0("Provided metannotate_data table does not contain the expected 'HMM_length' column. "),
         "Are you sure that you have already run map_naming_information()?")
  }

  # Filter by e-value cutoff and report stats to user
  futile.logger::flog.info(paste0("Filtering by e-value cutoff of ", evalue))
  metannotate_data_filtered <- filter_by_evalue(metannotate_data, evalue = evalue)
  metannotate_data <- metannotate_data_filtered$metannotate_data
  futile.logger::flog.info("Percent change from e-value filtration:")
  if (quietly == FALSE) {
    print(metannotate_data_filtered$read_counts$percent_change)
  }
  # TODO - optionally output the info to the user
  
  # Collapse the table to the desired taxonomic rank
  futile.logger::flog.info(paste0("Collapsing table to taxonomic rank '", taxon, "'"))
  metannotate_data_collapsed <- collapse_by_taxon(metannotate_data, taxon = taxon)
  
  # Normalize the data by HMM length
  futile.logger::flog.info("Normalizing data")
  metannotate_data_normalized_list <- normalize(metannotate_data_collapsed, normalizing_HMM = normalizing_HMM)

  futile.logger::flog.info("Combining any replicates")
  metannotate_data_normalized_list <- combine_replicates(metannotate_data_normalized_list)
  if (quietly == FALSE) {
    futile.logger::flog.info("Total normalized % abundance of analyzed genes compared to the marker gene:")
    print(metannotate_data_normalized_list$total_normalized_hits)
  }

  # Make plots
  futile.logger::flog.info("Plotting data")
  metannotate_plot <- visualize(metannotate_data_normalized_list = metannotate_data_normalized_list,
                                colouring_template_filename = colouring_template_filename,
                                top_x = top_x,
                                percent_mode = percent_mode,
                                normalizing_HMM = normalizing_HMM,
                                ...)

  return(metannotate_plot)
}
