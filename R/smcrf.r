#' Approximate Bayesian Computation sequential Monte Carlo via random forests
#'
#' `smcrf()` uses random forests to find the posterior distribution(s) for one or more parameters in a model.
#' It implements the sequential Monte Carlo framework, where each iteration
#' uses either ABC-RF (functions `regAbcrf` and `predict` in R package `abcrf`)
#' or ABC-DRF (functions `drf` and `predict` in R package `drf`) to update the posterior distribution(s).
#'
#' @param method Random forest method to implement in each iteration ("smcrf-single-param" by default).
#' method = "smcrf-single-param": implements ABC-RF for each parameter and results in their marginal posterior distributions.
#' method = "smcrf-multi-param": implements ABC-DRF for all parameters and results in the joint posterior distribution.
#' @param statistics_target A dataframe containing statistics from data.
#' Column names are the statistics IDs.
#' `smcrf()` only supports one row of statistics.
#' If there are multiple observations, we recommend applying `smcrf()` to each row individually.
#' @param smcrf_results An existing ABC-SMC-RF result.
#' If provided, smcrf will continue ABC-SMC-RF from the last iteration of the previous run.
#' @param model Model for the statistics.
#' The function must take two inputs: a data frame parameters and logic variable parallel.
#' @param perturb blah blah blah...
#' @param bounds ...
#' @param parameters_initial ...
#' @param nParticles ...
#' @param parallel ...
#' @param n_cores ...
#' @param ... Additional arguments to be passed to `abcrf` or `drf`.
#' @return An object `smcrf_results` containing the results of the inference.
#' If the posterior distributions have not converged to a satisfactory level,
#' the user may continue with `smcrf(smcrf_results = smcrf_results, ...)`,
#' in which case ABC-SMC-(D)RF will continue iterating from the last run in `smcrf_results`.
#' @export
#' @examples
#' library(abcsmcrf)
#' #---------------------------ABC-SMC-RF for a model with one parameter
#' #    Data to be fitted, consisting of two statistics s1 and s2:
#' statistics_target <- data.frame(s1 = 0, s2 = 2)
#' #    Parametrized model for the statistics:
#' model <- function(parameters) {
#'     statistics <- data.frame(
#'         s1 = parameters$theta - 1 + runif(nrow(parameters)),
#'         s2 = parameters$theta + 1 + runif(nrow(parameters))
#'     )
#'     cbind(parameters, statistics)
#' }
#' #    Function to perturb the parameters with random noise:
#' perturb <- function(parameters) {
#'     parameters$theta <- parameters$theta + runif(nrow(parameters), min = -0.1, max = 0.1)
#'     return(parameters)
#' }
#' #    Initial guesses for the parameter theta:
#' parameters_initial <- data.frame(theta = runif(100000, -10, 10))
#' #    Make sure that theta stays within bounds of prior distribution:
#' bounds <- data.frame(
#'     parameter = c("theta"),
#'     min = c(-10),
#'     max = c(10)
#' )
#' #    Run ABC-SMC-RF:
#' smcrf_results <- smcrf(
#'     method = "smcrf-single-param",
#'     statistics_target = statistics_target,
#'     model = model,
#'     perturb = perturb,
#'     bounds = bounds,
#'     parameters_initial = parameters_initial,
#'     nParticles = c(100, 100, 100),
#' )
smcrf <- function(method = "smcrf-single-param",
                  statistics_target = NULL,
                  smcrf_results = NULL,
                  model,
                  perturb,
                  bounds = NULL,
                  parameters_initial = NULL,
                  nParticles,
                  parallel = FALSE,
                  n_cores = NULL,
                  ...) {
    if (method == "smcrf-single-param") {
        return(smcrf_single_param(
            statistics_target = statistics_target,
            model = model,
            perturb = perturb,
            bounds = bounds,
            parameters_initial = parameters_initial,
            smcrf_single_param_results = smcrf_results,
            nParticles = nParticles,
            parallel = parallel,
            n_cores = n_cores,
            ...
        ))
    } else if (method == "smcrf-multi-param") {
        return(smcrf_multi_param(
            statistics_target = statistics_target,
            model = model,
            perturb = perturb,
            bounds = bounds,
            parameters_initial = parameters_initial,
            smcrf_multi_param_results = smcrf_results,
            nParticles = nParticles,
            parallel = parallel,
            n_cores = n_cores,
            ...
        ))
    } else {
        stop("Invalid method for smcrf")
    }
}

#' .......
#' @export
smcrf_single_param <- function(statistics_target = NULL,
                               model,
                               perturb,
                               bounds = NULL,
                               parameters_initial = NULL,
                               nParticles,
                               parallel,
                               n_cores = NULL,
                               save_model = TRUE,
                               save_rda = FALSE,
                               filename_rda = "ABCSMCRF.rda",
                               smcrf_single_param_results = NULL,
                               ...) {
    library(abcrf)
    library(Hmisc)
    nSimulations <<- 0
    cat("\n\n==============================================================================================================================================\n")
    #---Obtain information from previous SMC-RF
    if (!is.null(smcrf_single_param_results)) {
        SMCRF <- smcrf_single_param_results
        old_nIterations <- SMCRF[["nIterations"]]
        begin_iteration <- old_nIterations + 1
        nIterations <- length(nParticles) + old_nIterations
        nParticles <- c(SMCRF[["nParticles"]], nParticles)
        statistics_target <- SMCRF[["statistics_target"]]
        parameters_ids <- colnames(SMCRF[["Iteration_1"]]$parameters)
    } else {
        #---Initialize information for SMC-RF
        SMCRF <- list()
        parameters_ids <- colnames(parameters_initial)
        nIterations <- length(nParticles)
        begin_iteration <- 1
    }
    SMCRF[["method"]] <- "smcrf-single-param"
    SMCRF[["nIterations"]] <- nIterations
    SMCRF[["nParticles"]] <- nParticles
    SMCRF[["statistics_target"]] <- statistics_target
    SMCRF[["parameters_labels"]] <- data.frame(parameter = parameters_ids)
    SMCRF[["statistics_labels"]] <- data.frame(ID = colnames(statistics_target))
    for (iteration in begin_iteration:(nIterations + 1)) {
        if (iteration == (nIterations + 1)) {
            cat(paste0("\n\n++++++++++++++++++++++++\nSimulation count: ", nSimulations, "\n++++++++++++++++++++++++\n\n\n"))
            cat("\n\nSIMULATING STATISTICS FROM FINAL POSTERIOR DISTRIBUTION...\n")
        } else {
            cat(paste0("SMC-RF FOR SINGLE PARAMETERS: iteration ", iteration, "...\n"))
        }
        #---Sample prior parameters for this round of iteration...
        if (!is.null(smcrf_single_param_results)) {
            if (iteration == (old_nIterations + 1)) {
                SMCRF[[paste0("Iteration_", (old_nIterations + 1))]] <- c()
                ABCRF_weights <- SMCRF[[paste0("Iteration_", old_nIterations)]]$weights
                parameters <- SMCRF[[paste0("Iteration_", old_nIterations)]]$parameters
            }
        }
        if (iteration == 1) {
            #   ... For iteration 1: sample from initial parameters
            parameters <- data.frame(parameters_initial[1:nParticles[iteration], ])
            colnames(parameters) <- parameters_ids
            parameters_unperturbed <- parameters
        } else {
            #   ... For later iterations:
            ifelse(iteration == (nIterations + 1), nrow <- nParticles[nIterations], nrow <- nParticles[iteration])
            parameters_unperturbed <- data.frame(matrix(NA, nrow = nrow, ncol = length(parameters_ids)))
            colnames(parameters_unperturbed) <- parameters_ids
            parameters_next <- data.frame(matrix(NA, nrow = nrow, ncol = length(parameters_ids)))
            colnames(parameters_next) <- parameters_ids
            for (parameter_id in parameters_ids) {
                invalid_indices <- 1:nrow
                while (length(invalid_indices) > 0) {
                    #   Sample parameters from previous posterior distribution
                    parameter_replace <- data.frame(sample(parameters[, parameter_id], size = length(invalid_indices), prob = ABCRF_weights[, parameter_id], replace = TRUE))
                    colnames(parameter_replace) <- parameter_id
                    if (length(invalid_indices) == nrow) parameters_unperturbed[[parameter_id]][invalid_indices] <- parameter_replace[[parameter_id]]
                    #   Perturb parameters
                    if (iteration < (nIterations + 1)) {
                        if (is.function(perturb)) {
                            parameter_replace <- perturb(parameters = parameter_replace)
                        } else if (perturb == "Beaumont") {
                            wtd.var <- var(data.frame(sample(parameters[, parameter_id], size = 10000, prob = ABCRF_weights[, parameter_id], replace = TRUE)))
                            parameter_replace[[parameter_id]] <- rnorm(
                                n = nrow(parameter_replace),
                                mean = parameter_replace[[parameter_id]],
                                sd = sqrt(2 * wtd.var)
                            )
                        }
                    }
                    #   Check if parameters are within bounds, otherwise redo the failed parameters
                    parameters_next[[parameter_id]][invalid_indices] <- parameter_replace[[parameter_id]]
                    if (is.null(bounds)) {
                        invalid_indices <- which(is.na(parameters_next[, parameter_id]))
                    } else {
                        invalid_indices <- which(
                            is.na(parameters_next[, parameter_id]) |
                                parameters_next[, parameter_id] < bounds$min[which(bounds$parameter == parameter_id)] |
                                parameters_next[, parameter_id] > bounds$max[which(bounds$parameter == parameter_id)]
                        )
                    }
                }
            }
            parameters <- parameters_next
        }
        #---Simulate statistics
        cat("Making model simulations...\n")
        reference <- model(parameters = parameters)
        #---Re-sample and perturb the parameters if there's NA in the reference table
        invalid_rows <- rowSums(is.na(reference)) == ncol(reference)
        while (any(invalid_rows)) {
            if (iteration == 1) {
                parameters <- data.frame(parameters_initial[nrow(parameters) + 1:(nrow(parameters) + sum(invalid_rows)), ])
                colnames(parameters) <- parameters_ids
                reference[invalid_rows, ] <- model(parameters = parameters)
                invalid_rows <- rowSums(is.na(reference)) == ncol(reference)
            } else {
                #   ... For later iterations:
                parameters_next <- data.frame(matrix(NA, nrow = sum(invalid_rows), ncol = length(parameters_ids)))
                colnames(parameters_next) <- parameters_ids
                parameters_tmp <- parameters
                for (parameter_id in parameters_ids) {
                    invalid_indices <- 1:sum(invalid_rows)
                    while (length(invalid_indices) > 0) {
                        #   Sample parameters from previous posterior distribution
                        parameter_replace <- data.frame(parameters_tmp[sample(nrow(parameters_tmp), size = length(invalid_indices), prob = ABCRF_weights[, parameter_id], replace = T), parameter_id])
                        colnames(parameter_replace) <- parameter_id
                        #   Perturb parameters
                        if (iteration < (nIterations + 1)) {
                            if (is.function(perturb)) {
                                parameter_replace <- perturb(parameters = parameter_replace)
                            } else if (perturb == "Beaumont") {
                                wtd.var <- var(data.frame(parameters_tmp[sample(nrow(parameters_tmp), size = 10000, prob = ABCRF_weights[, parameter_id], replace = T), parameter_id]))
                                parameter_replace[[parameter_id]] <- rnorm(
                                    n = nrow(parameter_replace),
                                    mean = parameter_replace[[parameter_id]],
                                    sd = sqrt(2 * wtd.var)
                                )
                            }
                        }
                        #   Check if parameters are within bounds, otherwise redo the failed parameters
                        parameters_next[[parameter_id]][invalid_indices] <- parameter_replace[[parameter_id]]
                        if (is.null(bounds)) {
                            invalid_indices <- which(is.na(parameters_next[, parameter_id]))
                        } else {
                            invalid_indices <- which(
                                is.na(parameters_next[, parameter_id]) |
                                    parameters_next[, parameter_id] < bounds$min[which(bounds$parameter == parameter_id)] |
                                    parameters_next[, parameter_id] > bounds$max[which(bounds$parameter == parameter_id)]
                            )
                        }
                    }
                }
                parameters_tmp <- parameters_next
                colnames(parameters_tmp) <- parameters_ids
                reference[invalid_rows, ] <- model(parameters = parameters_tmp)
                invalid_rows <- rowSums(is.na(reference)) == ncol(reference)
            }
        }
        statistics <- data.frame(reference[, colnames(reference)[!colnames(reference) %in% parameters_ids]])
        parameters <- data.frame(reference[, colnames(reference)[colnames(reference) %in% parameters_ids]])
        colnames(parameters) <- parameters_ids
        colnames(statistics) <- colnames(reference)[!colnames(reference) %in% parameters_ids]
        #---Save SMC-RF results from the last iteration
        if (iteration == (nIterations + 1)) {
            SMCRF_iteration <- list()
            SMCRF_iteration$reference <- reference
            SMCRF_iteration$parameters <- parameters
            SMCRF_iteration$parameters_unperturbed <- parameters_unperturbed
            SMCRF_iteration$statistics <- statistics
            SMCRF[[paste0("Iteration_", iteration)]] <- SMCRF_iteration
            if (save_rda == TRUE) {
                save(SMCRF, file = filename_rda)
            }
            break
        }
        #---Run ABCRF for each parameter
        cat("Performing Random Forest prediction...\n")
        ABCRF_weights <- data.frame(matrix(NA, nrow = nParticles[iteration], ncol = 0))
        RFmodels <- list()
        posterior_gamma_RFs <- list()
        for (parameter_id in parameters_ids) {
            mini_reference <- reference[, c(parameter_id, colnames(reference)[!colnames(reference) %in% parameters_ids])]
            colnames(mini_reference)[1] <- "para"
            f <- as.formula("para ~.")
            RFmodel <- regAbcrf(
                formula = f,
                data = mini_reference,
                paral = parallel,
                ...
            )
            posterior_gamma_RF <- predict(
                object = RFmodel,
                obs = statistics_target,
                training = mini_reference,
                paral = parallel,
                rf.weights = T
            )
            ABCRF_weights[, parameter_id] <- posterior_gamma_RF$weights
            if (save_model == TRUE) {
                RFmodels[[parameter_id]] <- RFmodel
                posterior_gamma_RFs[[parameter_id]] <- posterior_gamma_RF
            }
        }
        cat("\n\n")
        #---Save SMC-RF results from this iteration
        SMCRF_iteration <- list()
        SMCRF_iteration$reference <- reference
        SMCRF_iteration$parameters <- parameters
        SMCRF_iteration$parameters_unperturbed <- parameters_unperturbed
        SMCRF_iteration$statistics <- statistics
        SMCRF_iteration$weights <- ABCRF_weights
        if (save_model == TRUE) {
            SMCRF_iteration$rf_model <- RFmodels
            SMCRF_iteration$rf_predict <- posterior_gamma_RFs
        }
        SMCRF[[paste0("Iteration_", iteration)]] <- SMCRF_iteration
        if (save_rda == TRUE) {
            save(SMCRF, file = filename_rda)
        }
    }
    return(SMCRF)
}

#' .......
#' @export
smcrf_multi_param <- function(statistics_target = NULL,
                              model,
                              perturb,
                              bounds = NULL,
                              parameters_initial = NULL,
                              nParticles,
                              parallel,
                              n_cores = NULL,
                              save_model = TRUE,
                              save_rda = FALSE,
                              filename_rda = "ABCSMCDRF.rda",
                              splitting.rule = "CART",
                              smcrf_multi_param_results = NULL,
                              ...) {
    library(drf)
    library(matrixStats)
    library(Hmisc)
    nSimulations <<- 0
    #---Obtain information from previous SMC-DRF
    if (!is.null(smcrf_multi_param_results)) {
        SMCDRF <- smcrf_multi_param_results
        old_nIterations <- SMCDRF[["nIterations"]]
        begin_iteration <- old_nIterations + 1
        nIterations <- length(nParticles) + old_nIterations
        nParticles <- c(SMCDRF[["nParticles"]], nParticles)
        statistics_target <- SMCDRF[["statistics_target"]]
        parameters_ids <- colnames(SMCDRF[["Iteration_1"]]$parameters)
    } else {
        #---Initialize information for SMC-DRF
        SMCDRF <- list()
        parameters_ids <- colnames(parameters_initial)
        nIterations <- length(nParticles)
        begin_iteration <- 1
    }
    SMCDRF[["method"]] <- "smcrf-multi-param"
    SMCDRF[["nIterations"]] <- nIterations
    SMCDRF[["nParticles"]] <- nParticles
    SMCDRF[["statistics_target"]] <- statistics_target
    SMCDRF[["parameters_labels"]] <- data.frame(parameter = parameters_ids)
    SMCDRF[["statistics_labels"]] <- data.frame(ID = colnames(statistics_target))
    for (iteration in begin_iteration:(nIterations + 1)) {
        if (iteration == (nIterations + 1)) {
            cat(paste0("\n\n++++++++++++++++++++++++\nSimulation count: ", nSimulations, "\n++++++++++++++++++++++++\n\n\n"))
            cat("\n\nSIMULATING STATISTICS FROM FINAL POSTERIOR DISTRIBUTION...\n")
        } else {
            cat(paste0("SMC-RF FOR MULTIPLE PARAMETERS: iteration ", iteration, "...\n"))
        }
        #---Sample prior parameters for this round of iteration...
        #---Sample prior parameters for this round of iteration...
        if (!is.null(smcrf_multi_param_results)) {
            if (iteration == (old_nIterations + 1)) {
                SMCDRF[[paste0("Iteration_", (old_nIterations + 1))]] <- c()
                DRF_weights <- SMCDRF[[paste0("Iteration_", old_nIterations)]]$weights
                parameters <- SMCDRF[[paste0("Iteration_", old_nIterations)]]$parameters
            }
        }
        if (iteration == 1) {
            #   ... For iteration 1: sample from initial parameters
            parameters <- data.frame(parameters_initial[1:nParticles[iteration], ])
            colnames(parameters) <- parameters_ids
            parameters_unperturbed <- parameters
        } else {
            #   ... For later iterations:
            ifelse(iteration == (nIterations + 1), nrow <- nParticles[nIterations], nrow <- nParticles[iteration])
            parameters_unperturbed <- data.frame(matrix(NA, nrow = nrow, ncol = length(parameters_ids)))
            colnames(parameters_unperturbed) <- parameters_ids
            parameters_next <- data.frame(matrix(NA, nrow = nrow, ncol = length(parameters_ids)))
            colnames(parameters_next) <- parameters_ids
            invalid_indices <- 1:nrow
            while (length(invalid_indices) > 0) {
                #   Sample parameters from previous posterior distribution
                parameter_replace <- data.frame(parameters[sample(nrow(parameters), size = length(invalid_indices), prob = DRF_weights[, 1], replace = T), ])
                colnames(parameter_replace) <- parameters_ids
                if (length(invalid_indices) == nrow) parameters_unperturbed <- parameter_replace
                #   Perturb parameters
                if (iteration < (nIterations + 1)) {
                    if (is.function(perturb)) {
                        parameter_replace <- perturb(parameters = parameter_replace)
                    } else if (perturb == "Beaumont") {
                        for (parameter_id in parameters_ids) {
                            if (length(invalid_indices) == nrow) wtd.var <- var(data.frame(parameters[sample(nrow(parameters), size = 10000, prob = DRF_weights[, 1], replace = T), parameter_id]))
                            parameter_replace[[parameter_id]] <- rnorm(
                                n = nrow(parameter_replace),
                                mean = parameter_replace[[parameter_id]],
                                sd = sqrt(2 * wtd.var)
                            )
                        }
                    }
                }
                #   Check if parameters are within bounds, otherwise redo the failed parameters
                parameters_next[invalid_indices, ] <- parameter_replace
                if (is.null(bounds)) {
                    invalid_indices <- which(apply(parameters_next, 1, function(x) any(is.na(x))))
                } else {
                    invalid_indices <- c()
                    for (parameter_id in parameters_ids) {
                        invalid_indices <- union(invalid_indices, which(
                            is.na(parameters_next[, parameter_id]) |
                                parameters_next[, parameter_id] < bounds$min[which(bounds$parameter == parameter_id)] |
                                parameters_next[, parameter_id] > bounds$max[which(bounds$parameter == parameter_id)]
                        ))
                    }
                }
            }
            parameters <- parameters_next
        }
        #---Simulate statistics
        cat("Making model simulations...\n")
        reference <- model(parameters = parameters)
        #---Re-sample and perturb the parameters if there's NA in the reference table
        invalid_rows <- rowSums(is.na(reference)) == ncol(reference)
        while (any(invalid_rows)) {
            if (iteration == 1) {
                parameters <- data.frame(parameters_initial[(nrow(parameters) + 1):(nrow(parameters) + sum(invalid_rows)), ])
                colnames(parameters) <- parameters_ids
                reference[invalid_rows, ] <- model(parameters = parameters)
                invalid_rows <- rowSums(is.na(reference)) == ncol(reference)
            } else {
                parameters_tmp <- parameters
                #   ... For later iterations:
                parameters_next <- data.frame(matrix(NA, nrow = sum(invalid_rows), ncol = length(parameters_ids)))
                colnames(parameters_next) <- parameters_ids
                invalid_indices <- 1:sum(invalid_rows)
                while (length(invalid_indices) > 0) {
                    #   Sample parameters from previous posterior distribution
                    parameter_replace <- data.frame(parameters_tmp[sample(nrow(parameters_tmp), size = length(invalid_indices), prob = DRF_weights[, 1], replace = T), ])
                    colnames(parameter_replace) <- parameters_ids
                    #   Perturb parameters
                    if (iteration < (nIterations + 1)) {
                        if (is.function(perturb)) {
                            parameter_replace <- perturb(parameters = parameter_replace)
                        } else if (perturb == "Beaumont") {
                            for (parameter_id in parameters_ids) {
                                wtd.var <- var(data.frame(parameters_tmp[sample(nrow(parameters_tmp), size = 10000, prob = DRF_weights[, 1], replace = T), parameter_id]))
                                parameter_replace[[parameter_id]] <- rnorm(
                                    n = nrow(parameter_replace),
                                    mean = parameter_replace[[parameter_id]],
                                    sd = sqrt(2 * wtd.var)
                                )
                            }
                        }
                    }
                    #   Check if parameters are within bounds, otherwise redo the failed parameters
                    parameters_next[invalid_indices, ] <- parameter_replace
                    if (is.null(bounds)) {
                        invalid_indices <- which(apply(parameters_next, 1, function(x) any(is.na(x))))
                    } else {
                        invalid_indices <- c()
                        for (parameter_id in parameters_ids) {
                            invalid_indices <- union(invalid_indices, which(
                                is.na(parameters_next[, parameter_id]) |
                                    parameters_next[, parameter_id] < bounds$min[which(bounds$parameter == parameter_id)] |
                                    parameters_next[, parameter_id] > bounds$max[which(bounds$parameter == parameter_id)]
                            ))
                        }
                    }
                }
                parameters_tmp <- parameters_next
                colnames(parameters_tmp) <- parameters_ids
                reference[invalid_rows, ] <- model(parameters = parameters_tmp)
                invalid_rows <- rowSums(is.na(reference)) == ncol(reference)
            }
        }
        statistics <- data.frame(reference[, colnames(reference)[!colnames(reference) %in% parameters_ids]])
        parameters <- data.frame(reference[, colnames(reference)[colnames(reference) %in% parameters_ids]])
        colnames(parameters) <- parameters_ids
        colnames(statistics) <- colnames(reference)[!colnames(reference) %in% parameters_ids]
        #---Save SMC-DRF results for the last iteration
        if (iteration == (nIterations + 1)) {
            SMCDRF_iteration <- list()
            SMCDRF_iteration$reference <- reference
            SMCDRF_iteration$parameters <- parameters
            SMCDRF_iteration$parameters_unperturbed <- parameters_unperturbed
            SMCDRF_iteration$statistics <- statistics
            SMCDRF[[paste0("Iteration_", iteration)]] <- SMCDRF_iteration
            if (save_rda == TRUE) {
                save(SMCDRF, file = filename_rda)
            }
            break
        }
        #---Run DRF for all parameter
        cat("Performing Random Forest prediction...\n")
        Xdrf <- statistics
        Ydrf <- reference[, parameters_ids]
        drfmodel <- drf(Xdrf, Ydrf, splitting.rule = splitting.rule, ...)
        def_pred <- predict(drfmodel, statistics_target)
        DRF_weights <- as.vector(get_sample_weights(drfmodel, statistics_target))
        DRF_weights <- data.frame(matrix(rep(DRF_weights, length(parameters_ids)), ncol = length(parameters_ids)))
        colnames(DRF_weights) <- parameters_ids
        # ---Save SMC-DRF results from this iteration
        SMCDRF_iteration <- list()
        SMCDRF_iteration$reference <- reference
        SMCDRF_iteration$parameters <- parameters
        SMCDRF_iteration$parameters_unperturbed <- parameters_unperturbed
        SMCDRF_iteration$statistics <- statistics
        SMCDRF_iteration$weights <- DRF_weights
        if (save_model == TRUE) {
            SMCDRF_iteration$rf_model <- drfmodel
            SMCDRF_iteration$rf_predict <- def_pred
        }
        SMCDRF[[paste0("Iteration_", iteration)]] <- SMCDRF_iteration
        if (save_rda == TRUE) {
            save(SMCDRF, file = filename_rda)
        }
    }
    return(SMCDRF)
}
