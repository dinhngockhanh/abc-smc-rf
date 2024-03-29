# # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Zijin - Macbook
# R_workplace <- "/Users/xiangzijin/Documents/ABC_SMCRF/0327_sfs"
# R_libPaths <- ""
# R_libPaths_extra <- "/Users/xiangzijin/SMC-RF/R"
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Zijin - Macmini
R_workplace <- "/Users/khanhngocdinh/Documents/Zijin/0328_sfs_coala"
R_libPaths <- ""
R_libPaths_extra <- "/Users/khanhngocdinh/Documents/Zijin/SMC-RF/R"
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Zhihan - Macbook
# R_workplace <- "/Users/lexie/Documents/DNA/SMC-RF/vignettes"
# R_libPaths <- ""
# R_libPaths_extra <- "/Users/lexie/Documents/DNA/SMC-RF/R"
# =======================================SET UP FOLDER PATHS & LIBRARIES
.libPaths(R_libPaths)
library(ggplot2)
library(gridExtra)
library(grid)
library(invgamma)
setwd(R_libPaths_extra)
files_sources <- list.files(pattern = "\\.[rR]$")
sapply(files_sources, source)
setwd(R_workplace)
# =========================Model for the Site Frequency Spectrum (SFS)
#   Input:  data frame of parameters, each row is one set of parameters
#   Output: data frame of parameters & statistics, each row contains statistics for one set of parameters:
#           first columns = input parameters
#           next columns = summary statistics
model <- function(parameters, parallel = TRUE) {
    nSamples <- 1000
    nNoise <- 0
    if (exists("nSimulations")) nSimulations <<- nSimulations + nrow(parameters)
    #   Make simulations & compute summary statistics (allele count)
    if (parallel) {
        library(parallel)
        library(pbapply)
        library(data.table)
        cl <- makePSOCKcluster(detectCores() - 1)
        clusterExport(cl, varlist = c("SFS_model"))
        stats <- pblapply(
            cl = cl, X = 1:nrow(parameters),
            FUN = function(i) {
                SFS_model(theta = parameters$theta[i], n = nSamples)
            }
        )
        stopCluster(cl)
        stats <- rbindlist(stats)
        class(stats) <- "data.frame"
    } else {
        stats <- c()
        for (i in 1:nrow(parameters)) {
            stats <- rbind(stats, SFS_model(theta = parameters$theta[i], n = nSamples))
        }
    }
    #   Add noise statistics
    noise <- matrix(runif(nrow(parameters) * nNoise), nrow(parameters), nNoise)
    #   Add column names
    data <- data.frame(cbind(stats, noise))
    if (nNoise > 0) {
        colnames(data) <- c(
            colnames(stats),
            paste0("noise_", c(1:nNoise))
        )
    }
    return(data)
}

SFS_model <- function(theta, n) {
    library(coala)
    model <- coal_model(n,20) +
        feat_mutation(par_const(theta))+
        sumstat_sfs()
    sim_data <- simulate(model, nsim = 1)
    sfs <- create_abc_sumstat(sim_data, model)
    sval <- sum(sfs)
    lvec <- floor(sqrt(n))
    # stats <- data.frame(matrix(c(theta, sval), nrow = 1))
    # colnames(stats) <- c("theta", "Mutation_count_S")
    stats <- data.frame(matrix(c(theta, sval, sfs[1:lvec]), nrow = 1))
    colnames(stats) <- c("theta", "Mutation_count_S", paste0("SFS_", 1:lvec))
    return(stats)
}
# =====================================================Target statistics
set.seed(1)
theta <- runif(1, 1, 20)
parameters_ground_truth <- data.frame(
    theta = theta
)
statistics_target <- model(parameters = parameters_ground_truth, parallel = FALSE)[-c(1:ncol(parameters_ground_truth))]
# ======================================Model for parameter perturbation
#   Input:  data frame of parameters, each row is one set of parameters
#   Output: data frame of parameters, after perturbation
perturb <- function(parameters) {
    for (i in 1:ncol(parameters)) parameters[[i]] <- parameters[[i]] + runif(nrow(parameters), min = -1, max = 1)
    return(parameters)
}
# ======================================Define ranges for the parameters
range <- data.frame(
    parameter = c("theta"),
    min = c(1),
    max = c(20)
)
# ========================================Initial guesses for parameters
# ====================================(sampled from prior distributions)
theta <- runif(10000, 1, 20)
parameters_initial <- data.frame(
    theta = theta
)
# ====================================Labels for parameters in the plots
parameters_labels <- data.frame(
    parameter = c("theta"),
    label = c(deparse(expression(theta)))
)
# ================================================================ABC-RF
#---Run ABC-RF
abcrf_results <- smcrf(
    method = "smcrf-single-param",
    statistics_target = statistics_target,
    parameters_initial = parameters_initial,
    model = model,
    perturb = perturb,
    range = range,
    nParticles = rep(8000, 1),
    parallel = TRUE
)
#---Plot posterior marginal distributions against other methods
plots <- plot_compare_marginal(
    # plots = plots,
    parameters_truth = parameters_ground_truth,
    abc_results = abcrf_results,
    parameters_labels = parameters_labels,
    plot_statistics = TRUE
)
# ========================================
#   Plot the out-of-bag estimates (equivalent to cross-validation)
abcrf_results[["Iteration_1"]][["rf_model"]][["model.rf"]]$predictions
abcrf_results$Iteration_1$rf_model$model.rf$predictions
abcrf_results$Iteration_1$parameters$theta
png(paste0("NEUTRAL_abcrf_theta_out_of_bag.png"))
plot(abcrf_results$Iteration_1$parameters$theta,
    abcrf_results$Iteration_1$rf_model$model.rf$predictions,
    xlab = "True value",
    ylab = "Out-of-bag estimate"
) + abline(a = 0, b = 1, col = "red")
dev.off()
#   Can the error be lowered by increasing the number of trees?
library(abcrf)
oob_error <- err.regAbcrf(abcrf_results$Iteration_1$rf_model, training = abcrf_results$Iteration_1$reference, paral = T)
png(paste0("NEUTRAL_abcrf_theta_error_by_ntree.png"))
plot(oob_error[, "ntree"], oob_error[, "oob_mse"], type = "l", xlab = "Number of trees", ylab = "Out-of-bag MSE")
dev.off()
#   Variance Importance of each statistic in inferring gamma
png(paste0("NEUTRAL_abcrf_theta_variance_importance.png"), width = 1500, height = 800, res = 150)
n.var <- min(30, length(abcrf_results$Iteration_1$rf_model$model.rf$variable.importance))
imp <- abcrf_results$Iteration_1$rf_model$model.rf$variable.importance
names(imp) <- colnames(statistics_target)
ord <- rev(order(imp, decreasing = TRUE)[1:n.var])
xmin <- 0
xlim <- c(xmin, max(imp) + 1)
dotchart(imp[ord], pch = 19, xlab = "Variable Importance", ylab = "", xlim = xlim, main = NULL, bg = "white", cex = 0.7)
dev.off()






# library(coala)
# model <- coal_model(10, 50) +
#     feat_mutation(par_const(4)) +
#     sumstat_sfs()

# sim_data <- simulate(model, nsim = 1)
# sfs <- create_abc_sumstat(sim_data, model)
# sim_data

# library(coala)
# model <- coal_model(10, 50) +
#     feat_mutation(par_named("theta")) +
#     sumstat_sfs()

# sim_data <- simulate(model, nsim = 1, pars = c(theta = 5))

# sim_data


# sim_param <- create_abc_param(sim_data, model)
# sim_param

# length(sim_data)

# sim_sumstat <- create_abc_sumstat(sim_data, model)

# nrow(sim_sumstat)
# row_sums <- apply(sim_sumstat, 1, sum)






















# row_sums
# activate_ms(priority = 500)
# model <- coal_model(20,500) +
#     feat_mutation(2) +
#     sumstat_sfs()

# stats <- simulate(model)
# barplot(stats$sfs)


# model <- coal_model(20, 500) + feat_mutation(2) + sumstat_sfs() stats <- simulate(model)
# barplot(stats$sfs)
# model <- coal_model(5, 1, 10) +
#     feat_mutation(5, model = "GTR", gtr_rates = rep(1, 6)) +
#     sumstat_dna()

# simulate(model)$dna


# model <- coal_model(c(10, 1), 1, 25) + feat_mutation(7.5, model = "GTR", gtr_rates = c(1, 1, 1, 1, 1, 1) / 6) +
#     feat_outgroup(2) + feat_pop_merge(1.0, 2, 1) + sumstat_dna()
# simulate(model)$dna