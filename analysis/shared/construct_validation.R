fit_construct_comparison <- function(data, response) {
  variables <- c(response, "BEEM", "age", "sex", "bmi", "season", "assessment_centre")
  if (!all(variables %in% names(data))) stop("Missing construct-validation variables.")
  data <- data[complete.cases(data[, variables]), variables, drop = FALSE]
  data$BEEM <- factor(data$BEEM, levels = c("low", "medium", "high"))
  for (v in c("age", "sex", "bmi", "season", "assessment_centre")) data[[v]] <- factor(data[[v]])
  lm(reformulate(variables[-1], response = response), data = data)
}

baseline_season <- function(date) {
  month <- as.integer(format(as.Date(date), "%m"))
  labels <- c("winter", "winter", "spring", "spring", "spring", "summer",
              "summer", "summer", "autumn", "autumn", "autumn", "winter")
  factor(labels[month], levels = c("winter", "spring", "summer", "autumn"))
}
