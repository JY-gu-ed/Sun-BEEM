registry_end_date <- function(registry) {
  setting <- switch(registry, hospital = "SUN_BEEM_HOSPITAL_END_DATE",
                    cancer = "SUN_BEEM_CANCER_END_DATE", stop("Unknown registry."))
  value <- Sys.getenv(setting)
  if (!nzchar(value)) stop("Set ", setting, " to the documented registry follow-up end date (YYYY-MM-DD).")
  date <- as.Date(value, format = "%Y-%m-%d")
  if (is.na(date)) stop("Invalid registry end date: ", setting)
  date
}

apply_incidence_followup <- function(data, status_col, event_date_col, cutoff,
                                     baseline_col = "date_assess") {
  required <- c(status_col, event_date_col, baseline_col, "death_date", "loss_to_followup_date")
  if (!all(required %in% names(data))) {
    stop("Rebuild the step4 incidence cohort with follow-up dates. Missing: ",
         paste(setdiff(required, names(data)), collapse = ", "))
  }
  baseline <- as.Date(data[[baseline_col]])
  event_date <- as.Date(data[[event_date_col]])
  death_date <- as.Date(data$death_date)
  loss_date <- as.Date(data$loss_to_followup_date)
  end <- pmin(rep(as.Date(cutoff), nrow(data)), death_date, loss_date, na.rm = TRUE)
  reported <- !is.na(data[[status_col]]) & as.character(data[[status_col]]) %in% c("yes", "1")
  complete <- !is.na(baseline) & !(reported & is.na(event_date))
  event <- reported & !is.na(event_date) & event_date > baseline & event_date <= end
  event[is.na(event)] <- FALSE
  observed_end <- end
  observed_end[event] <- event_date[event]
  time <- as.numeric(observed_end - baseline)
  data[[status_col]] <- ifelse(event, "yes", "no")
  data$diag_time <- time
  data[complete & is.finite(time) & time > 0, , drop = FALSE]
}
