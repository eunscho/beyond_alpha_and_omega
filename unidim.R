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
  # generate
  #========================================================================
  generate <- function(conditions, condition_number, rep_set, rep) {
    #========================================================================
    # seed and variable names
    #========================================================================
    set.seed(10000 * condition_number + 1000 * rep_set + rep)
    n <- as.integer(conditions[condition_number, 1])
    k <- as.integer(conditions[condition_number, 2])
    fskew <- as.double(conditions[condition_number, 3])
    fkurt <- as.double(conditions[condition_number, 4])
    # ========================================================================
    # true reliability and omega h
    #=========================================================================
    loadings <- .3 + .6 * runif(k) # range between .3 and .9
    avg <- mean(loadings)
    true_score <- loadings%*%t(loadings)
    error_score <- 1 - diag(true_score)
    rel <- sum(true_score) / (sum(true_score) + sum(error_score))
    dev <- sum(diag(true_score)) - sum(true_score) / k
    rel_dev <- (k / (k - 1)) * (dev / (sum(true_score) + sum(error_score)))
    sqrt_error_cov <- matrix(rep(0, k^2), nrow = k)
    diag(sqrt_error_cov) <- sqrt(error_score)
    ##############################################################
    # Random data
    #############################################################
    Skewness <- c(fskew, rep(0, k))
    Kurtosis <- c(fkurt, rep(0, k))
    Sigma <- matrix(rep(0, (k + 1)^2), nrow = k + 1)
    diag(Sigma) <- 1
    for (i in 1:n) {
      randoms <- mvrnonnorm(n = n, mu = rep(0, k + 1), Sigma = Sigma,
                            skewness = Skewness, kurtosis = Kurtosis)
      factor_score <- randoms[1]
      error_score <- randoms[2:(k + 1)]
      if (i == 1) {
        sample_obs <- factor_score %*% loadings + error_score %*% sqrt_error_cov
      } else {
        sample_obs <- rbind(sample_obs,
                            factor_score %*% loadings + error_score %*% sqrt_error_cov)
      }
    }
    
    out <- list(sample_obs = sample_obs,
                rel = rel,
                avg = avg,
                rel_dev = rel_dev)
    
    return(out)
  }
  #========================================================================
  # analyze
  #========================================================================
  analyze <- function(conditions, condition_number, rep_set, rep, data) {
    n <- as.integer(conditions[condition_number, 1])
    k <- as.integer(conditions[condition_number, 2])
    fskew <- as.double(conditions[condition_number, 3])
    fkurt <- as.double(conditions[condition_number, 4])
    m <- var(data$sample_obs)
    r <- cov2cor(m)
    rel <- data$rel
    avg <- data$avg
    rel_dev <- data$rel_dev
    #===========================================================================
    # reliability
    #===========================================================================
    alpha <- reliacoef::alpha(m, print = F)
    lambda2 <- reliacoef::mu1(m, print = F)
    mu2 <- reliacoef::mu2(m, print = F)
    std_alpha <- reliacoef::alpha(r, print = F)
    std_lambda2 <- reliacoef::mu1(r, print = F)
    std_mu2 <- reliacoef::mu2(r, print = F)
    glb <- psych::glb.fa(m)$glb
    joreskog <- reliacoef::joreskog(m, print = F )
    posr <- uni_cfa(r)
    joreskog_r <- sum(posr$lambda)^2/(sum(posr$lambda)^2 + sum(posr$theta))
    safeomega <- purrr::safely(psych::omega)
    minres <- safeomega(r, nfactors = 1, fm = "minres")$result$omega.tot
    kaiser <- reliacoef::kaisercaffrey(r, print = F)
    gilmer <- reliacoef::gilmer(m)
    out <- tibble::tibble(condition_number, n, k, fskew, fkurt, rep_set, rep, rel,
                          avg, rel_dev, alpha, std_alpha, lambda2, std_lambda2,
                          mu2, std_mu2, glb, joreskog, joreskog_r, minres,
                          gilmer, kaiser)
    return(out)
  }
  #========================================================================
  # Main Loop
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
