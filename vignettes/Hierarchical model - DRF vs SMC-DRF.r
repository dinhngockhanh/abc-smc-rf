library(abcsmcrf)
library(invgamma)
set.seed(1)
# ===================================Function for plotting extrame cases
plot_hierarchical_extreme <- function(drf_results,
                                      smcdrf_results,
                                      parameters_truth,
                                      limits) {
    library(ggplot2)
    library(patchwork)
    color_drf <- "#FFD320"
    color_smcdrf <- "#C03728"
    color_truth <- "#2E2A2B"
    #---Get parameters from each method
    Y_drf <- drf_results[[paste0("Iteration_", drf_results$nIterations + 1)]]$parameters_unperturbed
    Y_drf <- Y_drf[sample(1:nrow(Y_drf), size = 400, replace = T), ]
    Y_smcdrf <- smcdrf_results[[paste0("Iteration_", smcdrf_results$nIterations + 1)]]$parameters_unperturbed
    Y_smcdrf <- Y_smcdrf[sample(1:nrow(Y_smcdrf), size = 400, replace = T), ]
    Y_truth <- parameters_truth
    #---TOP plot = marginal distribution for theta_1
    drf_density <- density(Y_drf$theta1)
    drf_df <- data.frame(x = drf_density$x, y = drf_density$y)
    p_top <- ggplot() +
        geom_histogram(
            data = Y_truth, aes(x = theta1, y = ..density..),
            fill = color_truth, color = color_truth, alpha = 1
        ) +
        geom_line(
            data = drf_df,
            aes(x = x, y = y),
            color = color_drf, size = 3
        ) +
        theme_void() +
        xlim(c(limits$min[which(limits$parameter == "theta1")], limits$max[which(limits$parameter == "theta1")]))
    smcdrf_density <- density(Y_smcdrf$theta1)
    smcdrf_df <- data.frame(x = smcdrf_density$x, y = smcdrf_density$y)
    p_top <- p_top +
        geom_line(
            data = smcdrf_df,
            aes(x = x, y = y),
            color = color_smcdrf, size = 3
        )
    #---RIGHT plot = marginal distribution for theta_2
    drf_density <- density(Y_drf$theta2)
    drf_df <- data.frame(x = drf_density$x, y = drf_density$y)
    p_right <- ggplot() +
        geom_histogram(
            data = Y_truth, aes(x = theta2, y = ..density..),
            fill = color_truth, color = color_truth, alpha = 1
        ) +
        geom_line(
            data = drf_df,
            aes(x = x, y = y),
            color = color_drf, size = 3
        ) +
        theme_void() +
        coord_flip() +
        xlim(c(limits$min[which(limits$parameter == "theta2")], limits$max[which(limits$parameter == "theta2")]))
    smcdrf_density <- density(Y_smcdrf$theta2)
    smcdrf_df <- data.frame(x = smcdrf_density$x, y = smcdrf_density$y)
    p_right <- p_right +
        geom_line(
            data = smcdrf_df,
            aes(x = x, y = y),
            color = color_smcdrf, size = 3
        )
    #---MAIN plot = joint distribution from DRF, TRUE and SMCRF
    p_main <- ggplot() +
        geom_density_2d_filled(data = Y_truth, aes(x = theta1, y = theta2)) +
        geom_density_2d(data = Y_drf, aes(x = theta1, y = theta2), colour = color_drf, linewidth = 2, bins = 6) +
        geom_density_2d(data = Y_smcdrf, aes(x = theta1, y = theta2), colour = color_smcdrf, linewidth = 2, bins = 6) +
        theme_minimal() +
        scale_fill_grey() +
        xlim(c(limits$min[which(limits$parameter == "theta1")], limits$max[which(limits$parameter == "theta1")])) +
        ylim(c(limits$min[which(limits$parameter == "theta2")], limits$max[which(limits$parameter == "theta2")])) +
        labs(x = expression(theta[1]), y = expression(theta[2])) +
        theme(
            text = element_text(size = 30),
            panel.background = element_rect(fill = "white", colour = "white"),
            plot.margin = unit(c(0, 0, 0, 0), "cm"),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            legend.position = "none"
        )
    # Combine the plots using patchwork
    combined_plot <- p_top + plot_spacer() + p_main + p_right +
        plot_layout(ncol = 2, nrow = 2, widths = c(4, 1), heights = c(1, 4))
    # Save the plot to a file
    file_name <- paste0("NORMAL_drf_contour_with_densities_theta1_theta2.png")
    ggsave(file_name, plot = combined_plot, width = 17, height = 17, units = "in", dpi = 150)
}
# ==============================Model for the hierarchical model example
#   Input:  data frame of parameters, each row is one set of parameters
#   Output: data frame of parameters & statistics, each row contains statistics for one set of parameters:
#           first columns = input parameters
#           next columns = summary statistics
model <- function(parameters) {
    theta1 <- parameters$theta1
    theta2 <- parameters$theta2
    nSamples <- 10
    nNoise <- 50
    #   Make simulations
    y <- matrix(NA, length(theta1), nSamples)
    for (i in 1:length(theta1)) {
        y[i, ] <- rnorm(nSamples, theta1[i], sqrt(theta2[i]))
    }
    #   Compute some summary statistics (mean, variance, MAD)
    summary <- matrix(NA, length(
        theta1
    ), 3)
    for (i in 1:length(theta1)) {
        summary[i, ] <- c(mean(y[i, ]), var(y[i, ]), mad(y[i, ]))
    }
    data <- cbind(parameters, summary)
    #   Add some others summary statistics
    x <- data[, -c(1:2)]
    x <- cbind(
        x[, 1] + x[, 2], x[, 1] + x[, 3], x[, 2] + x[, 3], x[, 1] + x[, 2] + x[, 3],
        x[, 1] * x[, 2], x[, 1] * x[, 3], x[, 2] * x[, 3], x[, 1] * x[, 2] * x[, 3]
    )
    data <- cbind(data, x)
    #   Add noise statistics
    noise <- matrix(runif(nrow(parameters) * nNoise), nrow(parameters), nNoise)
    data <- cbind(data, noise)
    #   Add column names
    data <- data.frame(data)
    colnames(data) <- c(
        "theta1", "theta2",
        "mean", "variance", "mad",
        "sum_mean_var",
        "sum_mean_mad",
        "sum_var_mad",
        "sum_mean_var_mad",
        "prod_mean_var",
        "prod_mean_mad",
        "prod_var_mad",
        "prod_mean_var_mad",
        paste0("noise_", c(1:nNoise))
    )
    return(data)
}
# =====================================================Target statistics
alpha <- 4
beta <- 5
theta2 <- 1 / rgamma(1, shape = alpha, rate = beta)
theta1 <- rnorm(1, 0, sqrt(theta2))
parameters_target <- data.frame(
    theta1 = theta1,
    theta2 = theta2
)
statistics_target <- model(parameters = parameters_target)[-c(1:ncol(parameters_target))]
# ====================================================Prior distribution
dprior <- function(parameters, parameter_id = "all") {
    vec_NAN <- which(parameters$theta2 <= 0)
    vec_notNAN <- which(parameters$theta2 > 0)
    probs <- rep(1, nrow(parameters))
    probs[vec_NAN] <- 0
    if (parameter_id %in% c("all", "theta2")) {
        probs[vec_notNAN] <- probs[vec_notNAN] * dinvgamma(parameters$theta2[vec_notNAN], shape = alpha, rate = beta)
    }
    if (parameter_id %in% c("all", "theta1")) {
        probs[vec_notNAN] <- probs[vec_notNAN] * dnorm(parameters$theta1[vec_notNAN], 0, sqrt(parameters$theta2[vec_notNAN]))
    }
    return(probs)
}
rprior <- function(Nparameters) {
    theta2 <- rinvgamma(Nparameters, shape = alpha, rate = beta)
    theta1 <- rnorm(Nparameters, 0, sqrt(theta2))
    return(data.frame(theta1 = theta1, theta2 = theta2))
}
# ======================Perturbation bounds (default Beaumont kernel)
parameter_bounds <- data.frame(
    parameter = c("theta1", "theta2"),
    min = c(-Inf, 0),
    max = c(Inf, Inf)
)
# ====================================Labels for parameters in the plots
parameters_labels <- data.frame(
    parameter = c("theta1", "theta2"),
    label = c(deparse(expression(theta[1])), deparse(expression(theta[2])))
)
# ==============================================True marginal posteriors
nSamples <- 10
s_2 <- statistics_target[, "variance"] * (nSamples - 1)
ybar <- statistics_target[, "mean"]
shape_theta2 <- nSamples / 2 + alpha
scale_theta2 <- 0.5 * (s_2 + 2 * beta + nSamples * ybar^2 / (nSamples + 1))
theta2_true <- rinvgamma(10000, shape = shape_theta2, scale = 1 / scale_theta2)
theta1_true <- rnorm(10000, mean = nSamples * ybar / (nSamples + 1), sd = sqrt(2 * theta2_true / (nSamples + 1)))
parameters_truth <- data.frame(
    theta1 = theta1_true,
    theta2 = theta2_true
)
# ===================================================================DRF
drf_results <- smcrf(
    method = "smcrf-multi-param",
    statistics_target = statistics_target,
    model = model,
    rprior = rprior,
    dprior = dprior,
    nParticles = rep(20000, 1),
    parallel = TRUE
)
# ===========================================================ABC-SMC-DRF
smcdrf_results <- smcrf(
    method = "smcrf-multi-param",
    statistics_target = statistics_target,
    model = model,
    rprior = rprior,
    dprior = dprior,
    parameter_bounds = parameter_bounds,
    nParticles = rep(20000, 20),
    parallel = TRUE
)
# =============================================Plot marginal joint plots
#---Plot marginal joint plots
plot_hierarchical_extreme(
    drf_results = drf_results,
    smcdrf_results = smcdrf_results,
    limits = data.frame(parameter = c("theta1", "theta2"), min = c(0, 0), max = c(5, 8)),
    parameters_truth = parameters_truth
)
