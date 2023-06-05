unidim <- function(start_con = 1, end_con = 120, start_set = 1, end_set = 10) {
  library(tidyverse)
  library(semTools)
  library(reliacoef)
  library(psych)
  library(Lambda4)
  library(tictoc)
  library(misty)
  # specify simulation conditions
  n <- c(50, 100, 250, 500, 1000)
  k <- c(3, 5, 7, 9)
  fskew = c(0, 2.5)
  fkurt = c(-1.2, 0, 14)
  conditions <- tidyr::crossing(n, k, fskew, fkurt) 
  colnames(conditions) <- c("n", "k", "fskew", "fkurt")
  condition_numbers <- start_con:end_con
  rep_sets <- start_set:end_set
  rep_per_set <- 1:1000
  #========================================================================
  # Loop
  #========================================================================
  for (condition_number in condition_numbers) {
    condition <- conditions[condition_number, ]
    print(condition)
    for (rep_set in rep_sets) {
      tictoc::tic()
      print(paste("Starting condition number", condition_number, "rep", rep_set))
      print(condition)
      filename <- paste0("unidim", condition_number, "-", rep_set, ".csv")
      if (!file.exists(filename)) {
        for (rep in rep_per_set) {
          cat("unidim: ", condition_number, "rep set: ", rep_set, "rep: ", rep)
          data <- generate(conditions, condition_number, rep_set, rep)
          if (rep == 1) {
            temp <- analyze(conditions, condition_number, rep_set, rep, data)
          } else {
            temp <- rbind(temp,
                         analyze(conditions, condition_number, rep_set, rep, data))
          }
        } # end of for (rep in rep_per_set)
        out <- temp
        readr::write_csv(data.frame(out), file = filename)
        print(out)
        tictoc::toc()
      } # end of if (!file.exists(filename))
    } # end of  for (rep_set in rep_sets)
  } # end of for (condition_number in condition_numbers)
} # end of function