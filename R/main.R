#' Encode a five nucleotide sequence
#'
#' Split each sequence into categorical nucleotide positions. This internal
#' helper uses the factor levels required by the supplied random forest.
#'
#' @param dna_strings Character vector of equal-length DNA sequences.
#' @return A data frame with one factor column per nucleotide position.
dna_encoding <- function(dna_strings) {
  nn <- nchar(dna_strings[1])
  seq_m <- matrix(unlist(strsplit(dna_strings, "")), ncol = nn, byrow = TRUE)
  colnames(seq_m) <- paste0("nt_pos", seq_len(nn))
  seq_df <- as.data.frame(seq_m)
  seq_df[] <- lapply(seq_df, factor, levels = c("A", "T", "C", "G"))
  seq_df
}

#' Predict m6A status for several candidate sites
#'
#' Encodes each five-base context, obtains the model's Positive class
#' probability, and labels scores strictly greater than the threshold Positive.
#' The course's exon length and junction distance inputs are already log2
#' transformed. The bundled example data demonstrate the expected format.
#'
#' @param ml_fit Fitted randomForest classification model, such as the bundled
#'   `rf_fit.rds` resource.
#' @param feature_df Data frame with columns `gc_content`, `RNA_type`,
#'   `RNA_region`, `exon_length`, `distance_to_junction`,
#'   `evolutionary_conservation`, and `DNA_5mer`. RNA types are mRNA, lincRNA,
#'   lncRNA, or pseudogene; regions are CDS, intron, 3'UTR, or 5'UTR.
#' @param positive_threshold Numeric cutoff; only scores greater than it are
#'   classified Positive. Default is 0.5.
#' @return The input data frame (with RNA columns converted to factors) plus
#'   numeric `predicted_m6A_prob` and character `predicted_m6A_status` columns.
#' @import randomForest
#' @importFrom stats predict
#' @export
#' @examples
#' ml_fit <- readRDS(system.file("extdata", "rf_fit.rds",
#'                              package = "m6APrediction", mustWork = TRUE))
#' feature_df <- read.csv(system.file("extdata", "m6A_input_example.csv",
#'                                    package = "m6APrediction", mustWork = TRUE))
#' prediction_multiple(ml_fit, head(feature_df, 3), positive_threshold = 0.6)
prediction_multiple <- function(ml_fit, feature_df, positive_threshold = 0.5) {
  stopifnot(all(c("gc_content", "RNA_type", "RNA_region", "exon_length",
                  "distance_to_junction", "evolutionary_conservation",
                  "DNA_5mer") %in% colnames(feature_df)))
  feature_df$RNA_type <- factor(feature_df$RNA_type,
                                levels = c("mRNA", "lincRNA", "lncRNA", "pseudogene"))
  feature_df$RNA_region <- factor(feature_df$RNA_region,
                                  levels = c("CDS", "intron", "3'UTR", "5'UTR"))
  inferred_prob <- predict(ml_fit,
                           newdata = cbind(feature_df, dna_encoding(feature_df$DNA_5mer)),
                           type = "prob")
  feature_df$predicted_m6A_prob <- inferred_prob[, "Positive"]
  feature_df$predicted_m6A_status <- "Negative"
  feature_df$predicted_m6A_status[
    feature_df$predicted_m6A_prob > positive_threshold
  ] <- "Positive"
  feature_df
}

#' Predict m6A status for one candidate site
#'
#' Creates a one-row feature table and calls [prediction_multiple()].
#' Scores equal to the threshold receive Negative status.
#'
#' @param ml_fit Fitted randomForest classification model.
#' @param gc_content Numeric GC content feature.
#' @param RNA_type One of mRNA, lincRNA, lncRNA, or pseudogene.
#' @param RNA_region One of CDS, intron, 3'UTR, or 5'UTR.
#' @param exon_length Numeric, already log2 transformed as in the training data.
#' @param distance_to_junction Numeric, already log2 transformed.
#' @param evolutionary_conservation Numeric conservation feature.
#' @param DNA_5mer Character string of five DNA bases with A at the center.
#' @param positive_threshold Numeric Positive class cutoff (default 0.5).
#' @return Named character vector containing `predicted_m6A_prob` and
#'   `predicted_m6A_status`. The probability is coerced to text by `c()`.
#' @export
#' @examples
#' ml_fit <- readRDS(system.file("extdata", "rf_fit.rds",
#'                              package = "m6APrediction", mustWork = TRUE))
#' prediction_single(ml_fit, gc_content = 0.6, RNA_type = "mRNA",
#'                   RNA_region = "CDS", exon_length = 12,
#'                   distance_to_junction = 5, evolutionary_conservation = 0.8,
#'                   DNA_5mer = "GGACA", positive_threshold = 0.5)
prediction_single <- function(ml_fit, gc_content, RNA_type, RNA_region,
                              exon_length, distance_to_junction,
                              evolutionary_conservation, DNA_5mer,
                              positive_threshold = 0.5) {
  pred_outcome <- prediction_multiple(
    ml_fit = ml_fit,
    feature_df = data.frame(
      gc_content = gc_content,
      RNA_type = factor(RNA_type,
                        levels = c("mRNA", "lincRNA", "lncRNA", "pseudogene")),
      RNA_region = factor(RNA_region,
                          levels = c("CDS", "intron", "3'UTR", "5'UTR")),
      exon_length = exon_length,
      distance_to_junction = distance_to_junction,
      evolutionary_conservation = evolutionary_conservation,
      DNA_5mer = DNA_5mer
    ),
    positive_threshold = positive_threshold
  )
  returned_vector <- c(pred_outcome$predicted_m6A_prob,
                       pred_outcome$predicted_m6A_status)
  names(returned_vector) <- c("predicted_m6A_prob", "predicted_m6A_status")
  returned_vector
}
