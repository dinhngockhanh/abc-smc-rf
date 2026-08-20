#' Approximate Bayesian Computation sequential Monte Carlo via random forests
#'
#' \code{\link{smcrf}} uses random forests to find the posterior distribution(s) for one or more parameters in a model.
#' It implements the sequential Monte Carlo framework, where each iteration
#' uses either ABC-RF (functions \code{regAbcrf} and \code{predict} in R package \code{abcrf})
#' or ABC-DRF (functions \code{drf} and \code{predict} in R package \code{drf}) to update the posterior distribution(s).
#'
#' @param method Random forest method to implement in each iteration (\code{"smcrf-single-param"} by default).
#' method = \code{"smcrf-single-param"}: implements ABC-RF for each parameter and results in their marginal posterior distributions.
#' method = \code{"smcrf-multi-param"}: implements ABC-DRF for all parameters and results in the joint posterior distribution.
#' @param statistics_target A dataframe containing statistics from data.
#' Column names are the statistics IDs.
#' \code{\link{smcrf}} only supports one row of statistics.
#' If there are multiple observations, we recommend applying \code{\link{smcrf}} to each row individually.
#' @param statistics_selection A dataframe indicating selection of statistics for fitting individual parameters (only works for method \code{"smcrf-single-param"}; \code{NULL} by default).
#' Each column's name matches one statistic ID, and each row's name matches one parameter ID.
#' The value is 1 if the statistic is used for the parameter, 0 otherwise.
#' @param smcrf_results An existing ABC-SMC-RF result.
#' If provided, smcrf will continue ABC-SMC-RF from the last iteration of the previous run.
#' @param model Model for the statistics.
#' The function must take two inputs: a dataframe \code{parameters} and logic variable \code{parallel}.
#' The model must output a reference table, where each row contains parameters for each simulation and corresponding statistics.
#' The column names of the reference table must match the parameter and statistics IDs.
#' @param rprior Function to generate particles from the prior distribution.
#' The function must take one input: \code{Nparameters}, the number of particles to generate.
#' The output is a dataframe where column names match parameter IDs,
#' and each row contains one parameter set.
#' @param dprior Function to compute the prior density.
#' The function must take two inputs: \code{parameters} and \code{parameter_id}.
#' The dataframe \code{parameters} contains parameter sets in each row, with column names as parameter IDs.
#' The \code{parameter_id} is either \code{"all"} or one of the parameter IDs.
#' The output is a vector of prior probabilities corresponding to rows in \code{parameters},
#' either for the parameter indicated by \code{parameter_id} or jointly for all parameters (if \code{parameter_id} = \code{"all"}).
#' @param rperturb Function to generate perturbed particles (\code{NULL} by default).
#' If \code{rperturb} and \code{dperturb} are both \code{NULL}, \code{\link{smcrf}} uses the
#' independent Beaumont kernel from \code{\link{make_beaumont_kernel}}.
#' If supplied, both \code{rperturb} and \code{dperturb} must be provided.
#' The function must take three inputs: \code{parameters_unperturbed}, \code{parameters_previous_sampled}, and \code{iteration}.
#' The dataframe \code{parameters_unperturbed} contains unperturbed parameter sets in each row.
#' The dataframe \code{parameters_previous_sampled} contains parameter sets sampled from the previous iteration in each row.
#' In both dataframes, the column names match parameter IDs.
#' The integer \code{iteration} indicates the index for current iteration.
#' The output is a dataframe where column names match parameter IDs,
#' and each row contains one perturbed parameter set corresponding to each row in \code{parameters_unperturbed}.
#' @param dperturb Function to compute the perturbation density (\code{NULL} by default).
#' If \code{rperturb} and \code{dperturb} are both \code{NULL}, \code{\link{smcrf}} uses the
#' independent Beaumont kernel from \code{\link{make_beaumont_kernel}}.
#' If supplied, both \code{rperturb} and \code{dperturb} must be provided.
#' The function must take five inputs: \code{parameters}, \code{parameters_previous}, \code{parameters_previous_sampled}, \code{iteration}, and \code{parameter_id}.
#' The dataframe \code{parameters} contains parameter sets in each row.
#' The dataframe \code{parameters_previous} contains one parameter set from the previous iteration.
#' The dataframe \code{parameters_previous_sampled} contains parameter sets sampled from the previous iteration in each row.
#' In all three dataframes, the column names match parameter IDs.
#' The integer \code{iteration} indicates the index for current iteration.
#' The \code{parameter_id} is either \code{"all"} or one of the parameter IDs.
#' The output is a vector of perturbation probabilities from \code{parameters_previous} to rows in \code{parameters},
#' either for the parameter indicated by \code{parameter_id} or jointly for all parameters (if \code{parameter_id} = \code{"all"}).
#' @param nParticles A vector of particle counts.
#' Each entry indicates the number of simulations (e.g. particles) in the corresponding iteration.
#' @param final_sample A logic variable (\code{TRUE} by default).
#' If \code{final_sample} = \code{TRUE}, parameters will be sampled from the final iteration and saved in iteration \code{length(nParticles)+1}.
#' @param model_redo_if_NA A logic variable (\code{FALSE} by default).
#' If \code{model_redo_if_NA} = \code{TRUE}, the particles where \code{model} returns \code{NA} will be simulated again.
#' @param verbose A logic variable (\code{TRUE} by default).
#' If \code{verbose} = \code{TRUE}, progress messages will be printed during the fitting process.
#' @param parallel A logic variable (\code{FALSE} by default).
#' If \code{parallel} = \code{TRUE}, the ABC-RF functions will be computed in parallel.
#' @param save_model A logic variable (\code{FALSE} by default).
#' If \code{save_model} = \code{TRUE}, the random forest will be saved in \code{\link{smcrf}}'s output.
#' @param save_rds A logic variable (\code{FALSE} by default).
#' If \code{save_rds} = \code{TRUE}, the ABC-SMC-RF results will be saved in an rds file.
#' @param filename_rds A string (\code{"ABCSMCDRF.rds"} by default).
#' If \code{save_rds} = \code{TRUE}, the output from ABC-SMC-(D)RF will be saved in a file with this name.
#' @param parameter_bounds Optional dataframe with columns \code{parameter}, \code{min}, and \code{max}
#' giving truncation intervals for the default Beaumont kernel.
#' Ignored when custom \code{rperturb} and \code{dperturb} are supplied.
#' Parameters not listed are treated as unbounded.
#' See \code{\link{make_beaumont_kernel}}.
#' @param ... Additional arguments to be passed to \code{abcrf} or \code{drf}.
#' @seealso \code{\link{make_beaumont_kernel}}
#' @return An object \code{smcrf_results} containing the results of the inference.
#' If the posterior distributions have not converged to a satisfactory level,
#' the user may continue with \code{smcrf(smcrf_results = smcrf_results, ...)},
#' in which case ABC-SMC-(D)RF will continue iterating from the last run in \code{smcrf_results}.
#' @export
smcrf <- function(method = "smcrf-single-param",
                  statistics_target = NULL,
                  statistics_selection = NULL,
                  smcrf_results = NULL,
                  model,
                  rprior,
                  dprior,
                  rperturb = NULL,
                  dperturb = NULL,
                  nParticles,
                  final_sample = TRUE,
                  model_redo_if_NA = FALSE,
                  verbose = TRUE,
                  parallel = FALSE,
                  save_model = TRUE,
                  save_rds = FALSE,
                  filename_rds = "ABCSMCDRF.rds",
                  parameter_bounds = NULL,
                  ...) {
    suppressPackageStartupMessages(library(matrixStats))
    suppressPackageStartupMessages(library(Hmisc))
    suppressPackageStartupMessages(library(crayon))
    kernel <- resolve_perturbation(
        rperturb = rperturb,
        dperturb = dperturb,
        parameter_bounds = parameter_bounds
    )
    rperturb <- kernel$rperturb
    dperturb <- kernel$dperturb
    if (method == "smcrf-multi-param" & !is.null(statistics_selection)) stop("statistics_selection is only available for method 'smcrf-single-param'")
    if (method == "smcrf-single-param") {
        return(smcrf_single_param(
            statistics_target = statistics_target,
            statistics_selection = statistics_selection,
            model = model,
            rprior = rprior,
            dprior = dprior,
            rperturb = rperturb,
            dperturb = dperturb,
            nParticles = nParticles,
            final_sample = final_sample,
            model_redo_if_NA = model_redo_if_NA,
            verbose = verbose,
            parallel = parallel,
            save_model = save_model,
            save_rds = save_rds,
            filename_rds = filename_rds,
            smcrf_single_param_results = smcrf_results,
            ...
        ))
    } else if (method == "smcrf-multi-param") {
        if (!is.null(statistics_selection)) stop("statistics_selection is only available for method 'smcrf-single-param'")
        return(smcrf_multi_param(
            statistics_target = statistics_target,
            model = model,
            rprior = rprior,
            dprior = dprior,
            rperturb = rperturb,
            dperturb = dperturb,
            nParticles = nParticles,
            final_sample = final_sample,
            model_redo_if_NA = model_redo_if_NA,
            verbose = verbose,
            parallel = parallel,
            save_model = save_model,
            save_rds = save_rds,
            filename_rds = filename_rds,
            smcrf_multi_param_results = smcrf_results,
            ...
        ))
    } else {
        stop("Invalid method for smcrf")
    }
}

smcrf_single_param <- function(statistics_target,
                               statistics_selection = NULL,
                               model,
                               rprior,
                               dprior,
                               rperturb,
                               dperturb,
                               nParticles,
                               final_sample,
                               model_redo_if_NA,
                               verbose,
                               parallel,
                               save_model = TRUE,
                               save_rds = FALSE,
                               filename_rds,
                               smcrf_single_param_results,
                               ...) {
    suppressPackageStartupMessages(library(abcrf))
    nSimulations <<- 0
    #---Obtain information from previous SMC-RF
    if (!is.null(smcrf_single_param_results)) {
        #   ... If continuing from a previous ABC-SMC-DRF chain:
        SMCRF <- smcrf_single_param_results
        parameters_ids <- colnames(SMCRF[["Iteration_1"]]$parameters)
        statistics_ids <- colnames(SMCRF[["Iteration_1"]]$statistics)
        nIterations <- length(nParticles) + SMCRF[["nIterations"]]
        iteration_start <- 1 + SMCRF[["nIterations"]]
        if (is.null(statistics_target)) statistics_target <- SMCRF[["statistics_target"]]
        nParticles <- c(SMCRF[["nParticles"]], nParticles)
        SMCRF[[paste0("Iteration_", iteration_start)]] <- c()
        ABCRF_weights <- SMCRF[[paste0("Iteration_", iteration_start - 1)]]$weights
        parameters <- SMCRF[[paste0("Iteration_", iteration_start - 1)]]$parameters
    } else {
        #   ... If starting a new ABC-SMC-DRF chain:
        SMCRF <- list()
        parameters_ids <- colnames(rprior(Nparameters = 1))
        statistics_ids <- colnames(statistics_target)
        nIterations <- length(nParticles)
        iteration_start <- 1
    }
    reference_ids <- colnames(model(parameters = rprior(Nparameters = 1)))
    iteration_end <- ifelse(final_sample, nIterations + 1, nIterations)
    SMCRF[["method"]] <- "smcrf-single-param"
    SMCRF[["nIterations"]] <- nIterations
    SMCRF[["nParticles"]] <- nParticles
    SMCRF[["statistics_target"]] <- statistics_target
    SMCRF[["parameters_labels"]] <- data.frame(parameter = parameters_ids)
    SMCRF[["statistics_labels"]] <- data.frame(ID = colnames(statistics_target))
    if (!verbose) {
        pb <- txtProgressBar(
            min = 0,
            max = iteration_end - iteration_start + 1,
            style = 3,
            width = 50,
            char = "+"
        )
    }
    for (iteration in iteration_start:iteration_end) {
        if (verbose) {
            if (iteration == (nIterations + 1)) {
                cat(bold(red("ABC-SMC-DRF FOR SINGLE PARAMETERS:")), paste0(bold(yellow("final posterior distribution", "\n"))))
            } else {
                cat(bold(red("ABC-SMC-DRF FOR SINGLE PARAMETERS:")), paste0(bold(yellow("iteration", iteration, "\n"))))
            }
        }
        #---Prepare sampled particles from the previous iteration
        if (iteration > 1) {
            parameters_previous_sampled <- data.frame(matrix(NA, nrow = 10000, ncol = length(parameters_ids)))
            colnames(parameters_previous_sampled) <- parameters_ids
            for (parameter_id in parameters_ids) {
                parameters_previous_sampled[, parameter_id] <- parameters[sample(nrow(parameters), size = 10000, prob = ABCRF_weights[, parameter_id], replace = TRUE), parameter_id]
            }
        }
        #---Create training set
        if (verbose) cat(blue("Sampling parameters and computing model simulations...\n"))
        nrow <- ifelse(iteration == (nIterations + 1), nParticles[nIterations], nParticles[iteration])
        invalid_indices <- 1:nrow
        parameters_unperturbed <- data.frame(matrix(NA, nrow = nrow, ncol = length(parameters_ids)))
        colnames(parameters_unperturbed) <- parameters_ids
        parameters_next <- data.frame(matrix(NA, nrow = nrow, ncol = length(parameters_ids)))
        colnames(parameters_next) <- parameters_ids
        reference_next <- data.frame(matrix(NA, nrow = nrow, ncol = length(parameters_ids) + ncol(statistics_target)))
        colnames(reference_next) <- reference_ids
        while (length(invalid_indices) > 0) {
            #   Generate new particles...
            if (iteration == 1) {
                #   ... for iteration 1:
                #   Sample from the prior distribution
                parameters_unperturbed[invalid_indices, ] <- rprior(Nparameters = length(invalid_indices))
                parameters_next[invalid_indices, ] <- parameters_unperturbed[invalid_indices, ]
            } else {
                #   ... for later iterations:
                #   Sample from the previous posterior distribution
                parameters_previous <- parameters
                weights_previous <- ABCRF_weights
                parameter_replace <- data.frame(matrix(NA, nrow = length(invalid_indices), ncol = length(parameters_ids)))
                colnames(parameter_replace) <- parameters_ids
                for (parameter_id in parameters_ids) {
                    parameter_replace[, parameter_id] <- parameters_previous[sample(nrow(parameters_previous), size = length(invalid_indices), prob = weights_previous[, parameter_id], replace = T), parameter_id]
                }
                parameters_unperturbed[invalid_indices, ] <- parameter_replace
                #   Perturb parameters
                if (iteration < (nIterations + 1)) {
                    parameter_replace <- rperturb(
                        parameters_unperturbed = parameter_replace,
                        parameters_previous_sampled = parameters_previous_sampled,
                        iteration = iteration
                    )
                }
                parameters_next[invalid_indices, ] <- parameter_replace
            }
            invalid_indices_next <- which(dprior(parameters_next, parameter_id = "all") <= 0)
            #   Generate statistics for the new particles
            ids <- setdiff(invalid_indices, invalid_indices_next)
            tmp <- parameters_next[ids, ]
            if (!is.data.frame(tmp)) {
                tmp <- as.data.frame(tmp, ncol = length(parameters_ids))
                colnames(tmp) <- parameters_ids
            }
            if (length(ids) > 0) reference_next[ids, ] <- model(parameters = tmp)
            #   Find particles that need to be regenerated
            invalid_indices <- invalid_indices_next
            if (model_redo_if_NA) invalid_indices <- unique(c(invalid_indices, which(apply(reference_next, 1, function(x) any(is.na(x))))))
        }
        reference <- reference_next
        parameters <- reference[, parameters_ids, drop = FALSE]
        statistics <- reference[, statistics_ids, drop = FALSE]
        #   Finish the last iteration
        if (iteration == (nIterations + 1)) {
            SMCRF_iteration <- list()
            SMCRF_iteration$reference <- reference
            SMCRF_iteration$parameters <- parameters
            SMCRF_iteration$parameters_unperturbed <- parameters_unperturbed
            SMCRF_iteration$statistics <- statistics
            SMCRF[[paste0("Iteration_", iteration)]] <- SMCRF_iteration
            if (save_rds == TRUE) {
                saveRDS(SMCRF, file = filename_rds)
            }
            if (!verbose) setTxtProgressBar(pb, iteration - iteration_start + 1)
            if (!verbose) cat("\n")
            return(SMCRF)
        }
        #---Run ABCRF for all parameters
        if (verbose) cat(blue("Performing Random Forest prediction...\n"))
        ABCRF_weights <- data.frame(matrix(NA, nrow = nParticles[iteration], ncol = 0))
        if (save_model == TRUE) {
            RFmodels <- list()
            posterior_gamma_RFs <- list()
        }
        for (parameter_id in parameters_ids) {
            #   Select statistics for the parameter
            if (is.null(statistics_selection)) {
                mini_reference <- reference[, c(parameter_id, colnames(reference)[!colnames(reference) %in% parameters_ids])]
            } else {
                mini_reference <- reference[, c(parameter_id, colnames(statistics_selection)[which(statistics_selection[rownames(statistics_selection) == parameter_id, ] == 1)])]
            }
            colnames(mini_reference)[1] <- "para"
            RFmodel <- regAbcrf(
                formula = as.formula("para ~."),
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
        #---Modify ABCRF weights
        if (verbose) cat(blue("Recalibrating Random Forest weights...\n"))
        if (iteration > 1) {
            for (parameter_id in parameters_ids) {
                #   Compute numerators for weight recalibration
                weight_modifiers_numerator <- dprior(parameters, parameter_id = parameter_id)
                #   Compute denominators for weight recalibration
                weight_modifiers_denominator <- rep(0, nrow(parameters))
                for (i in 1:nrow(parameters_previous)) {
                    weight_modifiers_denominator_i <- rep(weights_previous[i, parameter_id], nrow(parameters)) *
                        dperturb(
                            parameters = parameters,
                            parameters_previous = parameters_previous[i, , drop = FALSE],
                            parameters_previous_sampled = parameters_previous_sampled,
                            iteration = iteration,
                            parameter_id = parameter_id
                        )
                    weight_modifiers_denominator <- weight_modifiers_denominator + weight_modifiers_denominator_i
                }
                #   Modify weights for new particles
                ABCRF_weights[, parameter_id] <- as.vector(ABCRF_weights[, parameter_id]) * weight_modifiers_numerator / weight_modifiers_denominator
                ABCRF_weights[, parameter_id] <- ABCRF_weights[, parameter_id] / sum(ABCRF_weights[, parameter_id])
            }
        }
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
        if (save_rds == TRUE) {
            saveRDS(SMCRF, file = filename_rds)
        }
        if (!verbose) setTxtProgressBar(pb, iteration - iteration_start + 1)
    }
    if (save_rds == TRUE) {
        saveRDS(SMCRF, file = filename_rds)
    }
    if (!verbose) cat("\n")
    return(SMCRF)
}

smcrf_multi_param <- function(statistics_target,
                              model,
                              rprior,
                              dprior,
                              rperturb,
                              dperturb,
                              nParticles,
                              final_sample,
                              model_redo_if_NA,
                              verbose,
                              parallel,
                              save_model,
                              save_rds,
                              filename_rds,
                              splitting.rule = "CART",
                              smcrf_multi_param_results,
                              ...) {
    suppressPackageStartupMessages(library(drf))
    nSimulations <<- 0
    #---Obtain information from previous SMC-DRF
    if (!is.null(smcrf_multi_param_results)) {
        #   ... If continuing from a previous ABC-SMC-DRF chain:
        SMCDRF <- smcrf_multi_param_results
        parameters_ids <- colnames(SMCDRF[["Iteration_1"]]$parameters)
        statistics_ids <- colnames(SMCDRF[["Iteration_1"]]$statistics)
        nIterations <- length(nParticles) + SMCDRF[["nIterations"]]
        iteration_start <- 1 + SMCDRF[["nIterations"]]
        if (is.null(statistics_target)) statistics_target <- SMCDRF[["statistics_target"]]
        nParticles <- c(SMCDRF[["nParticles"]], nParticles)
        SMCDRF[[paste0("Iteration_", iteration_start)]] <- c()
        DRF_weights <- SMCDRF[[paste0("Iteration_", iteration_start - 1)]]$weights
        parameters <- SMCDRF[[paste0("Iteration_", iteration_start - 1)]]$parameters
    } else {
        #   ... If starting a new ABC-SMC-DRF chain:
        SMCDRF <- list()
        parameters_ids <- colnames(rprior(Nparameters = 1))
        statistics_ids <- colnames(statistics_target)
        nIterations <- length(nParticles)
        iteration_start <- 1
    }
    reference_ids <- colnames(model(parameters = rprior(Nparameters = 1)))
    iteration_end <- ifelse(final_sample, nIterations + 1, nIterations)
    SMCDRF[["method"]] <- "smcrf-multi-param"
    SMCDRF[["nIterations"]] <- nIterations
    SMCDRF[["nParticles"]] <- nParticles
    SMCDRF[["statistics_target"]] <- statistics_target
    SMCDRF[["parameters_labels"]] <- data.frame(parameter = parameters_ids)
    SMCDRF[["statistics_labels"]] <- data.frame(ID = statistics_ids)
    if (!verbose) {
        pb <- txtProgressBar(
            min = 0,
            max = iteration_end - iteration_start + 1,
            style = 3,
            width = 50,
            char = "+"
        )
    }
    for (iteration in iteration_start:iteration_end) {
        if (verbose) {
            if (iteration == (nIterations + 1)) {
                cat(bold(red("ABC-SMC-DRF FOR MULTIPLE PARAMETERS:")), paste0(bold(yellow("final posterior distribution", "\n"))))
            } else {
                cat(bold(red("ABC-SMC-DRF FOR MULTIPLE PARAMETERS:")), paste0(bold(yellow("iteration", iteration, "\n"))))
            }
        }
        #---Prepare sampled particles from the previous iteration
        if (iteration > 1) {
            parameters_previous_sampled <- as.data.frame(parameters[sample(nrow(parameters), size = 10000, prob = DRF_weights[, 1], replace = TRUE), , drop = FALSE])
        }
        #---Create training set
        if (verbose) cat(blue("Sampling parameters and computing model simulations...\n"))
        nrow <- ifelse(iteration == (nIterations + 1), nParticles[nIterations], nParticles[iteration])
        invalid_indices <- 1:nrow
        parameters_unperturbed <- data.frame(matrix(NA, nrow = nrow, ncol = length(parameters_ids)))
        colnames(parameters_unperturbed) <- parameters_ids
        parameters_next <- data.frame(matrix(NA, nrow = nrow, ncol = length(parameters_ids)))
        colnames(parameters_next) <- parameters_ids
        reference_next <- data.frame(matrix(NA, nrow = nrow, ncol = length(parameters_ids) + ncol(statistics_target)))
        colnames(reference_next) <- reference_ids
        while (length(invalid_indices) > 0) {
            #   Generate new particles...
            if (iteration == 1) {
                #   ... for iteration 1:
                #   Sample from the prior distribution
                parameters_unperturbed[invalid_indices, ] <- rprior(Nparameters = length(invalid_indices))
                parameters_next[invalid_indices, ] <- parameters_unperturbed[invalid_indices, ]
            } else {
                #   ... for later iterations:
                #   Sample from the previous posterior distribution
                parameters_previous <- parameters
                weights_previous <- DRF_weights
                parameter_replace <- parameters[sample(nrow(parameters), size = length(invalid_indices), prob = DRF_weights[, 1], replace = T), ]
                parameters_unperturbed[invalid_indices, ] <- parameter_replace
                #   Perturb parameters
                if (iteration < (nIterations + 1)) {
                    parameter_replace <- rperturb(
                        parameters_unperturbed = parameter_replace,
                        parameters_previous_sampled = parameters_previous_sampled,
                        iteration = iteration
                    )
                }
                parameters_next[invalid_indices, ] <- parameter_replace
            }
            invalid_indices_next <- which(dprior(parameters_next, parameter_id = "all") <= 0)
            #   Generate statistics for the new particles
            ids <- setdiff(invalid_indices, invalid_indices_next)
            tmp <- parameters_next[ids, ]
            if (!is.data.frame(tmp)) {
                tmp <- as.data.frame(tmp, ncol = length(parameters_ids))
                colnames(tmp) <- parameters_ids
            }
            if (length(ids) > 0) reference_next[ids, ] <- model(parameters = tmp)
            #   Find particles that need to be regenerated
            invalid_indices <- invalid_indices_next
            if (model_redo_if_NA) invalid_indices <- unique(c(invalid_indices, which(apply(reference_next, 1, function(x) any(is.na(x))))))
        }
        reference <- reference_next
        parameters <- reference[, parameters_ids, drop = FALSE]
        statistics <- reference[, statistics_ids, drop = FALSE]
        #   Finish the last iteration
        if (iteration == (nIterations + 1)) {
            SMCDRF_iteration <- list()
            SMCDRF_iteration$reference <- reference
            SMCDRF_iteration$parameters <- parameters
            SMCDRF_iteration$parameters_unperturbed <- parameters_unperturbed
            SMCDRF_iteration$statistics <- statistics
            SMCDRF[[paste0("Iteration_", iteration)]] <- SMCDRF_iteration
            if (save_rds == TRUE) {
                saveRDS(SMCDRF, file = filename_rds)
            }
            if (!verbose) setTxtProgressBar(pb, iteration - iteration_start + 1)
            if (!verbose) cat("\n")
            return(SMCDRF)
        }
        #---Run DRF for all parameters
        if (verbose) cat(blue("Performing Random Forest prediction...\n"))
        drfmodel <- drf(
            X = statistics,
            Y = parameters,
            splitting.rule = splitting.rule,
            ...
        )
        def_pred <- predict(
            object = drfmodel,
            newdata = statistics_target
        )
        DRF_weights <- as.vector(get_sample_weights(drfmodel, statistics_target))
        #---Modify DRF weights
        if (verbose) cat(blue("Recalibrating Random Forest weights...\n"))
        if (iteration > 1) {
            #   Compute numerators for weight recalibration
            weight_modifiers_numerator <- dprior(parameters, parameter_id = "all")
            #   Compute denominators for weight recalibration
            weight_modifiers_denominator <- rep(0, nrow(parameters))
            for (i in 1:nrow(parameters_previous)) {
                weight_modifiers_denominator_i <- rep(weights_previous[i, 1], nrow(parameters)) *
                    dperturb(
                        parameters = parameters,
                        parameters_previous = parameters_previous[i, , drop = FALSE],
                        parameters_previous_sampled = parameters_previous_sampled,
                        iteration = iteration,
                        parameter_id = "all"
                    )
                weight_modifiers_denominator <- weight_modifiers_denominator + weight_modifiers_denominator_i
            }
            #   Modify weights for new particles
            DRF_weights <- DRF_weights * weight_modifiers_numerator / weight_modifiers_denominator
            DRF_weights <- DRF_weights / sum(DRF_weights)
        }
        DRF_weights <- data.frame(matrix(rep(DRF_weights, length(parameters_ids)), ncol = length(parameters_ids)))
        colnames(DRF_weights) <- parameters_ids
        #---Save SMC-DRF results from this iteration
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
        if (save_rds == TRUE) {
            saveRDS(SMCDRF, file = filename_rds)
        }
        if (!verbose) setTxtProgressBar(pb, iteration - iteration_start + 1)
    }
    if (save_rds == TRUE) {
        saveRDS(SMCDRF, file = filename_rds)
    }
    if (!verbose) cat("\n")
    return(SMCDRF)
}
