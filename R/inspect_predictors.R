#' Diagnostic plots of predictor distributions
#'
#' For each predictor, draws a distribution panel intended to help the user
#' visually pick filtering boundaries before calling [fit_beinf()]:
#' \itemize{
#'   \item Predictors named in `categorical` (default: `"n_isoforms"`) are
#'         shown as bar plots, with each bar annotated by the percent of
#'         reads at that level.
#'   \item All other predictors are shown as densities with dashed lines
#'         marking `quantile_probs` (default 5\% and 95\%), labelled with
#'         the cutpoint values.
#' }
#'
#' Quantile values and per-level proportions are returned invisibly so the
#' user can use them when picking `max_n_isoforms` or
#' `pct_covered_winsor_bounds` for [clean_features()].
#'
#' Typical workflow:
#' \preformatted{
#'   features <- extract_features(bam, gtf, fasta)
#'   info     <- inspect_predictors(features)   # inspect raw distributions
#'   clean    <- clean_features(features, max_n_isoforms = 7)
#'   fit      <- fit_beinf(clean)
#' }
#'
#' @param features Data frame from [extract_features()] or
#'   [clean_features()]. Safe to pass the raw output of
#'   [extract_features()] before any cleaning.
#' @param predictors Character vector of columns to plot. Default is the
#'   four canonical predictors.
#' @param categorical Character vector of predictors to draw as bar plots
#'   (one bar per unique value, annotated with percent of reads). Default
#'   `"n_isoforms"`.
#' @param primary_only Logical; subset to primary alignments
#'   (`aln_type == "P"`) before plotting. Defaults to `TRUE` so the
#'   panels match the data [fit_beinf()] actually trains on.
#' @param quantile_probs Length-2 numeric of quantile probabilities to
#'   mark on the numeric panels. Default `c(0.05, 0.95)`.
#'
#' @return A `patchwork` object (printed as a side effect). Invisibly
#'   returns a named list with one element per predictor: for categorical
#'   predictors, a data frame of `level` / `prop`; for numeric
#'   predictors, the two-element quantile vector.
#'
#' @examples
#' ex       <- calmer_example_data()
#' features <- extract_features(ex$bam, ex$gtf, ex$fasta)
#' info     <- inspect_predictors(features)   # inspect raw primaries
#' info$tx_length    # 5th and 95th percentile cutpoints
#'
#' @importFrom ggplot2 ggplot aes geom_bar geom_text geom_density
#'   geom_vline annotate labs theme_minimal scale_y_continuous expansion
#' @importFrom patchwork wrap_plots
#' @importFrom rlang .data
#' @importFrom stats quantile
#' @export
inspect_predictors <- function(features,
                               predictors     = c("n_isoforms", "tx_length",
                                                  "gc_content", "pct_covered"),
                               categorical    = "n_isoforms",
                               primary_only   = TRUE,
                               quantile_probs = c(0.05, 0.95)) {
  stopifnot(
    is.data.frame(features),
    all(predictors %in% names(features)),
    is.logical(primary_only), length(primary_only) == 1L,
    is.numeric(quantile_probs),
    length(quantile_probs) == 2L,
    quantile_probs[1] < quantile_probs[2],
    quantile_probs[1] >= 0,
    quantile_probs[2] <= 1
  )

  if (primary_only) {
    if (!"aln_type" %in% names(features)) {
      stop("primary_only = TRUE but `aln_type` column is missing.",
           call. = FALSE)
    }
    features <- features[features$aln_type == "P", , drop = FALSE]
  }

  plots        <- list()
  summary_info <- list()

  for (p in predictors) {
    x <- features[[p]]
    x <- x[!is.na(x)]

    if (p %in% categorical) {
      tab <- prop.table(table(x))
      df  <- data.frame(level = names(tab),
                        prop  = as.numeric(tab),
                        stringsAsFactors = FALSE)
      num_order <- suppressWarnings(as.numeric(df$level))
      if (!any(is.na(num_order))) {
        df$level <- factor(df$level, levels = df$level[order(num_order)])
      } else {
        df$level <- factor(df$level, levels = df$level)
      }
      summary_info[[p]] <- df

      plots[[p]] <- ggplot2::ggplot(df, ggplot2::aes(x = .data$level,
                                                     y = .data$prop)) +
        ggplot2::geom_bar(stat = "identity", fill = "steelblue") +
        ggplot2::geom_text(
          ggplot2::aes(label = sprintf("%.1f%%", 100 * .data$prop)),
          vjust = -0.3, size = 2.2
        ) +
        ggplot2::scale_y_continuous(
          labels = function(v) sprintf("%.0f%%", 100 * v),
          expand = ggplot2::expansion(mult = c(0, 0.12))
        ) +
        ggplot2::labs(title = p, x = p, y = "% of reads") +
        ggplot2::theme_minimal()

    } else {
      qs <- stats::quantile(x, probs = quantile_probs, na.rm = TRUE)
      summary_info[[p]] <- qs

      plots[[p]] <- ggplot2::ggplot(data.frame(x = x),
                                    ggplot2::aes(x = .data$x)) +
        ggplot2::geom_density(fill = "steelblue", alpha = 0.3) +
        ggplot2::geom_vline(xintercept = qs, linetype = "dashed",
                            colour = "firebrick") +
        ggplot2::annotate(
          "text", x = qs, y = Inf,
          label = sprintf("%s: %.2f", names(qs), qs),
          vjust = 1.5, hjust = -0.05, size = 2.2, colour = "firebrick"
        ) +
        ggplot2::labs(title = p, x = p, y = "density") +
        ggplot2::theme_minimal()
    }
  }

  print(patchwork::wrap_plots(plots, ncol = 2))
  invisible(summary_info)
}
