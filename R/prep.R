#' Prepare responses and response times for model fitting
#'
#' Assembles item responses and item-level response times into the format
#' expected by [grmrt_fit()]. The function reorders items into their order of
#' presentation, recodes reverse-keyed items, converts times to log seconds,
#' and sets implausible timings to missing. With \code{drop_first = TRUE}, the
#' time of the first presented item is set to missing because, in web
#' questionnaires, it typically includes reading the instructions and page
#' orientation. Respondents are never removed; timings outside the plausible
#' window are treated as missing at the person-item level.
#'
#' @param data A data frame.
#' @param items Names of the response columns.
#' @param times Names of the time columns, in the same order as \code{items}.
#' @param order Optional character vector giving the presentation order of
#'   \code{items}. Defaults to the order of \code{items}.
#' @param time_unit Unit of the time columns: \code{"ms"} or \code{"s"}.
#' @param reverse Optional names of reverse-keyed items to recode.
#' @param scale_range Minimum and maximum of the response scale, used for
#'   reverse coding. Defaults to the observed range of each item.
#' @param missing_codes Response codes to treat as missing (for example, 0).
#' @param min_time,max_time Plausible window for item times, in seconds.
#' @param drop_first Logical; set the first presented item's time to missing.
#' @return An object of class \code{grmrt_data} with elements
#'   \code{responses}, \code{logrt}, and \code{screening} (a per-item summary).
#' @examples
#' sim <- simulate_grmrt(n_persons = 100, n_items = 5)
#' df <- data.frame(sim$responses, exp(sim$logrt) * 1000)
#' names(df) <- c(paste0("i", 1:5), paste0("i", 1:5, "_ms"))
#' d <- prep_rt(df, items = paste0("i", 1:5), times = paste0("i", 1:5, "_ms"))
#' d
#' @export
prep_rt <- function(data, items, times, order = NULL, time_unit = c("ms", "s"),
                    reverse = NULL, scale_range = NULL, missing_codes = NULL,
                    min_time = 0.5, max_time = 60, drop_first = TRUE) {
  time_unit <- match.arg(time_unit)
  if (length(items) != length(times)) stop("'items' and 'times' must have equal length.")
  miss <- setdiff(c(items, times), names(data))
  if (length(miss)) stop("Columns not found: ", paste(miss, collapse = ", "))
  if (!is.null(order)) {
    if (!setequal(order, items)) stop("'order' must contain exactly the names in 'items'.")
    idx <- match(order, items); items <- items[idx]; times <- times[idx]
  }
  X <- as.matrix(data[, items, drop = FALSE])
  if (!is.null(missing_codes)) X[X %in% missing_codes] <- NA
  if (!is.null(reverse)) {
    bad <- setdiff(reverse, items)
    if (length(bad)) stop("Unknown item(s) in 'reverse': ", paste(bad, collapse = ", "))
    for (r in reverse) {
      rg <- if (is.null(scale_range)) range(X[, r], na.rm = TRUE) else scale_range
      X[, r] <- rg[1] + rg[2] - X[, r]
    }
  }
  Tm <- as.matrix(data[, times, drop = FALSE])
  storage.mode(Tm) <- "double"
  if (time_unit == "ms") Tm <- Tm / 1000
  fast <- !is.na(Tm) & Tm < min_time
  slow <- !is.na(Tm) & Tm > max_time
  Tm[fast | slow] <- NA
  if (drop_first) Tm[, 1L] <- NA
  Y <- log(Tm)
  colnames(X) <- colnames(Y) <- items
  scr <- data.frame(item = items, position = seq_along(items),
                    missing_response = colSums(is.na(X)),
                    too_fast = colSums(fast), too_slow = colSums(slow),
                    missing_time = colSums(is.na(Y)), row.names = NULL)
  structure(list(responses = X, logrt = Y, screening = scr,
                 window = c(min_time, max_time), drop_first = drop_first),
            class = "grmrt_data")
}

#' @export
print.grmrt_data <- function(x, ...) {
  cat("grmrt data:", nrow(x$responses), "respondents,", ncol(x$responses), "items\n")
  cat(sprintf("Plausible time window: %.2f to %.0f s; first item time dropped: %s\n",
              x$window[1], x$window[2], x$drop_first))
  cat(sprintf("Timings set to missing: %d too fast, %d too slow (%.2f%% of all)\n",
              sum(x$screening$too_fast), sum(x$screening$too_slow),
              100 * mean(is.na(x$logrt[, -1L]))))
  invisible(x)
}
