#' Plot distribution(s) in each iteration from ABC-SMC-(D)RF result
#'
#' \code{\link{plot_smcrf_marginal}} plot the marginal distribution(s) for each iteration from an
#' Approximate Bayesian Computation sequential Monte Carlo via random forest result.
#'
#' @param smcrf_results An ABC-SMC-(D)RF result containing the inference distributions of parameters from each iteration.
#' @param parameters_truth A dataframe containing true values of parameters from the ground-truth distributions.
#' If provided, the function will plot the true values or distributions of parameters.
#' @param parameters_labels A dataframe containing labels in the plots for corresponding parameters.
#' If provided, parameter labels will be exhibited on the plots' axes.
#' @param statistics_labels A dataframe containing labels in the plots for corresponding statistics.
#' If provided, statistics labels will exhibit on the plots' axes.
#' @param plot_statistics A logic variable (plot_statistics = FALSE by default).
#' If plot_statistics = TRUE, the marginal distributions in each iteration for corresponding statistics will also be output.
#' @param xlimit A dataframe containing the maximum and minimum bounds for parameters.
#' If provided, the x-axis will be scaled by them.
#' @param alpha The numeric number to modify transparency. Default is 0.3.
#' @param plot_hist A logic variable (plot_hist = FALSE by default).
#' If plot_hist = TRUE, marginal distributions will be plotted in histograms.
#' @seealso
#' \code{\link{smcrf}}
#'
#' @export
#' @examples
#' #    Dataframe containing the true parameters
#' parameters_truth <- data.frame(
#'     theta = 2 # this will plot a vertical line
#' )
#' # if you prefer to have a histogram of true posterior
#' parameters_truth <- data.frame(
#'     theta = rnorm(10000, mean = 2, 1)
#' )
#' #    Dataframe containing the parameter labels
#' parameters_labels <- data.frame(
#'     parameter = c("theta"),
#'     label = c(deparse(expression(theta)))
#' )
#' #    Dataframe containing the x-axis bounds
#' xlimit <- data.frame(
#'     parameter = c("theta"),
#'     min = c(1),
#'     max = c(20)
#' )
plot_smcrf_marginal <- function(smcrf_results,
                                parameters_truth = NULL,
                                parameters_labels = NULL,
                                statistics_labels = NULL,
                                plot_statistics = FALSE,
                                xlimit = NULL,
                                alpha = 0.3,
                                plot_hist = FALSE) {
    library(ggplot2)
    nIterations <- smcrf_results[["nIterations"]]
    if (is.null(parameters_labels)) parameters_labels <- smcrf_results[["parameters_labels"]]
    if (is.null(statistics_labels)) statistics_labels <- smcrf_results[["statistics_labels"]]
    #---Set up color scheme for plotting
    color_scheme <- c(
        "True Posterior Distribution" = "#666666",
        "Data statistic" = "#666666",
        "Prior Distribution" = "#E5E5E5",
        setNames(hcl.colors(max(nIterations, 2), palette = "Earth"), paste0("Iter. ", 1:nIterations))
    )
    #---Set up legend order for plotting
    legend_order <- c("True Posterior Distribution", "Data statistic", "Prior Distribution", paste0("Iter. ", 1:nIterations))
    #---Plot marginal distributions for each parameter
    for (i in 1:length(parameters_labels$parameter)) {
        parameter_id <- parameters_labels$parameter[i]
        #   Begin plot
        p <- ggplot()
        #   Plot True Posterior distribution (if provided)
        if (!is.null(parameters_truth)) {
            true_posterior_df <- data.frame(value = parameters_truth[[parameter_id]], legend = "True Posterior Distribution")
            if (nrow(parameters_truth) == 1) {
                p <- p + geom_vline(data = true_posterior_df, aes(xintercept = value), linetype = "solid", linewidth = 5)
            } else {
                if (plot_hist) {
                    p <- p + geom_histogram(data = true_posterior_df, aes(x = value, y = ..density.., fill = legend, color = legend), alpha = 1)
                } else {
                    p <- p + geom_density(data = true_posterior_df, aes(x = value, fill = legend, color = legend), alpha = 1, linewidth = 2)
                }
            }
        }
        #   Plot prior distribution
        prior_df <- data.frame(value = smcrf_results[["Iteration_1"]]$parameters_unperturbed[[parameter_id]], legend = "Prior Distribution")
        if (plot_hist) {
            p <- p + geom_histogram(data = prior_df, aes(x = value, y = ..density.., fill = legend, color = legend), alpha = alpha)
        } else {
            p <- p + geom_density(data = prior_df, aes(x = value, fill = legend, color = legend), alpha = alpha, linewidth = 2)
        }
        #   Plot posterior distribution for each iteration
        for (iteration in 1:nIterations) {
            posterior_df <- data.frame(
                value = smcrf_results[[paste0("Iteration_", iteration + 1)]]$parameters_unperturbed[[parameter_id]],
                legend = paste0("Iter. ", iteration)
            )
            if (plot_hist) {
                p <- p + geom_histogram(data = posterior_df, aes(x = value, y = ..density.., fill = legend, color = legend), alpha = alpha)
            } else {
                p <- p + geom_density(data = posterior_df, aes(x = value, fill = legend, color = legend), alpha = alpha, linewidth = 2)
            }
            if (!is.null(xlimit)) {
                p <- p +
                    xlim(xlimit$min[which(xlimit$parameter == parameter_id)], xlimit$max[which(xlimit$parameter == parameter_id)])
            }
        }
        #   Add label for parameter
        if ("label" %in% colnames(parameters_labels)) {
            p <- p + labs(x = eval(parse(text = parameters_labels$label[i])))
        } else {
            p <- p + labs(x = parameter_id)
        }
        #   Beautify plot
        p <- p +
            scale_fill_manual(values = color_scheme, name = "", breaks = legend_order) +
            scale_color_manual(values = color_scheme, name = "", breaks = legend_order) +
            guides(fill = guide_legend(nrow = 1, keywidth = 2.5, keyheight = 1)) +
            theme(
                text = element_text(size = 50),
                panel.background = element_rect(fill = "white", colour = "white"),
                panel.grid.major = element_line(colour = "white"),
                panel.grid.minor = element_line(colour = "white"),
                legend.position = "top",
                legend.justification = c(0, 0.5)
            )
        #   Print marginal distribution plot
        file_name <- paste0(smcrf_results[["method"]], "-marginal-parameter=", parameter_id, ".png")
        png(file_name, res = 150, width = 30, height = 15, units = "in", pointsize = 12)
        print(p)
        dev.off()
    }
    #---Plot marginal distributions for each statistic
    if (plot_statistics) {
        for (i in 1:length(statistics_labels$ID)) {
            statistic_id <- statistics_labels$ID[i]
            p <- ggplot()
            #   Plot statistic value from target data
            data_df <- data.frame(value = smcrf_results$statistics_target[[statistic_id]], legend = "Data statistic")
            p <- p + geom_vline(data = data_df, aes(xintercept = value), linetype = "solid", linewidth = 5)
            #   Plot prior distribution
            prior_df <- data.frame(value = smcrf_results[["Iteration_1"]]$statistics[[statistic_id]], legend = "Prior Distribution")
            if (plot_hist) {
                p <- p + geom_histogram(data = prior_df, aes(x = value, y = ..density.., fill = legend, color = legend), alpha = alpha)
            } else {
                p <- p + geom_density(data = prior_df, aes(x = value, fill = legend, color = legend), alpha = alpha, linewidth = 2)
            }

            #   Plot posterior distribution for each iteration
            for (iteration in 1:nIterations) {
                posterior_df <- data.frame(
                    value = smcrf_results[[paste0("Iteration_", iteration + 1)]]$statistics[[statistic_id]],
                    legend = paste0("Iter. ", iteration)
                )
                if (plot_hist) {
                    p <- p + geom_histogram(data = posterior_df, aes(x = value, y = ..density.., fill = legend, color = legend), alpha = alpha)
                } else {
                    p <- p + geom_density(data = posterior_df, aes(x = value, fill = legend, color = legend), alpha = alpha, linewidth = 2)
                }
            }
            #   Add label for parameter
            if ("label" %in% colnames(statistics_labels)) {
                p <- p + labs(x = eval(parse(text = statistics_labels$label[i])))
            } else {
                p <- p + labs(x = statistic_id)
            }
            #   Beautify plot
            p <- p +
                scale_fill_manual(values = color_scheme, name = "", breaks = legend_order) +
                scale_color_manual(values = color_scheme, name = "", breaks = legend_order) +
                guides(fill = guide_legend(nrow = 1, keywidth = 2.5, keyheight = 1)) +
                theme(
                    text = element_text(size = 50),
                    panel.background = element_rect(fill = "white", colour = "white"),
                    panel.grid.major = element_line(colour = "white"),
                    panel.grid.minor = element_line(colour = "white"),
                    legend.position = "top",
                    legend.justification = c(0, 0.5)
                )
            #   Print marginal distribution plot
            file_name <- paste0(smcrf_results[["method"]], "-marginal-statistic=", statistic_id, ".png")
            png(file_name, res = 150, width = 30, height = 15, units = "in", pointsize = 12)
            print(p)
            dev.off()
        }
    }
}
#' Plot and compare marginal posterior distribution(s) from ABC-SMC-(D)RF result
#'
#' \code{\link{plot_compare_marginal}} plots the marginal posterior distribution(s) for the provided ABC-SMC-(D)RF result.
#' It can also compare the marginal posterior distributions for the provided ABC-SMC-(D)RF result with several ABC-SMC-(D)RF plots results.
#'
#' @param plots An existed ABC-SMC-(D)RF marginal plots result.
#' If provided, \code{\link{plot_compare_marginal}} will plot the new ABC-SMC-(D)RF result and compare it with provided plots results.
#' If plots = NULL, \code{\link{plot_compare_marginal}} will make a new plot for the ABC-SMC-(D)RF result.
#' @param abc_results An ABC-SMC-(D)RF result.
#' Will be plotted out by the function.
#' @param parameters_truth A dataframe containing true values of parameters from the ground-truth distributions.
#' If provided, the function will plot the true values or distributions of parameters.
#' @param parameters_labels A dataframe containing labels in the plots for corresponding parameters.
#' If provided, parameter labels will be exhibited on the plots' axes.
#' @param statistics_labels A dataframe containing labels in the plots for corresponding statistics.
#' If provided, statistics labels will exhibit on the plots' axes.
#' @param xlimit A dataframe containing the maximum and minimum bounds for parameters.
#' If provided, the x-axis will be scaled by them.
#' @param plot_statistics A logic variable (plot_statistics = FALSE by default).
#' If plot_statistics = TRUE, the marginal distributions in each iteration for corresponding statistics will also be output.
#' @param plot_truth_hist A logic variable (TRUE by default).
#' If plot_truth_hist = TRUE and \code{parameters_truth} contains more than one row,
#' the true posterior is drawn as a histogram; otherwise it is drawn as a density or line.
#' @param plot_hist A logic variable (plot_hist = FALSE by default).
#' If plot_hist = TRUE, marginal distributions will be plotted in histograms.
#' @param plot_hist_point A logic variable (plot_hist_point = FALSE by default).
#' If plot_hist_point = TRUE, marginal distributions will be plotted in histograms with points in the middle.
#' @param breaks Optional numeric vector of histogram break points, passed to \code{geom_histogram} and \code{hist}.
#' If \code{NULL}, ggplot2 / \code{hist} defaults are used.
#' @param alpha_truth Transparency for the true posterior distribution. Default is 0.8.
#' @param alpha The numeric number to modify transparency. Default is 0.3.
#' @param plot_prior A logic variable (plot_prior = FALSE by default)
#' If plot_prior = TRUE, the prior distribution will be plotted out.
#' @return A list of \code{ggplot2} objects containing the marginal plots results of posterior distributions.
#' The user can use the function to compare ABC-SMC-(D)RF marginal plots with the marginal posterior distribution(s)
#' of other ABC-SMC-(D)RF result(s).
#' @seealso
#' \code{\link{smcrf}}
#' @export
plot_compare_marginal <- function(plots = NULL,
                                  abc_results,
                                  parameters_truth = NULL,
                                  parameters_labels = NULL,
                                  statistics_labels = NULL,
                                  xlimit = NULL,
                                  plot_statistics = FALSE,
                                  plot_truth_hist = TRUE,
                                  plot_hist = FALSE,
                                  plot_hist_point = FALSE,
                                  breaks = NULL,
                                  alpha_truth = 0.8,
                                  alpha = 0.3,
                                  plot_prior = FALSE) {
    library(ggplot2)
    if (is.null(parameters_labels)) parameters_labels <- abc_results[["parameters_labels"]]
    if (is.null(statistics_labels)) statistics_labels <- abc_results[["statistics_labels"]]
    method <- abc_results[["method"]]
    #---Set up color scheme for plotting
    color_scheme <- c(
        "Prior Distribution" = "#E5E5E5",
        "True Posterior Distribution" = "#666666",
        "ABC-REJ" = "#20854E",
        "ABC-RF" = "#7876B1",
        "ABC-DRF" = "#0072B5",
        "MCMC" = "#E18727",
        "ABC-MCMC" = "#E18727",
        "ABC-SMC" = "#E18727",
        "ABC-SMC-RF" = "#BC3C29",
        "ABC-SMC-DRF" = "#BC3C29"
    )
    #---Set up legend order for plotting
    legend_order <- c(
        "Prior Distribution",
        "True Posterior Distribution",
        "ABC-REJ",
        "ABC-MCMC",
        "ABC-SMC",
        "ABC-RF",
        "ABC-DRF",
        "MCMC",
        "ABC-SMC-RF",
        "ABC-SMC-DRF"
    )
    #---Begin plots
    if (is.null(plots)) {
        plots <- list()
        for (parameter_id in parameters_labels$parameter) {
            plots$parameters[[parameter_id]] <- ggplot()
        }
        if (plot_statistics) {
            for (statistic_id in statistics_labels$ID) {
                plots$statistics[[statistic_id]] <- ggplot()
            }
        }
        new_plot <- TRUE
    } else {
        new_plot <- FALSE
    }
    #---Plot prior distributions (if provided)
    if (plot_prior) {
        for (parameter_id in parameters_labels$parameter) {
            prior_df <- data.frame(value = abc_results[["Iteration_1"]]$parameters[[parameter_id]], legend = "Prior Distribution")
            if (plot_hist) {
                plots$parameters[[parameter_id]] <- plots$parameters[[parameter_id]] +
                    geom_histogram(
                        data = prior_df,
                        aes(x = value, y = ..density.., fill = legend, color = legend),
                        breaks = breaks, alpha = alpha
                    )
            } else {
                plots$parameters[[parameter_id]] <- plots$parameters[[parameter_id]] +
                    geom_density(
                        data = prior_df,
                        aes(x = value, fill = legend, color = legend),
                        alpha = alpha, linewidth = 2
                    )
            }
        }
    }
    #---Plot True Posterior Distribution parameter distributions (if provided)
    if (!is.null(parameters_truth)) {
        for (parameter_id in parameters_labels$parameter) {
            if (any(grepl("density", colnames(parameters_truth)))) {
                true_posterior_df <- data.frame(value = parameters_truth[[parameter_id]], density = parameters_truth[["density"]], legend = "True Posterior Distribution")
                if (plot_truth_hist) {
                    plots$parameters[[parameter_id]] <- plots$parameters[[parameter_id]] +
                        geom_histogram(
                            data = true_posterior_df,
                            aes(x = value, weight = density, y = ..density.., fill = legend, color = legend),
                            breaks = breaks, alpha = alpha_truth
                        )
                } else {
                    plots$parameters[[parameter_id]] <- plots$parameters[[parameter_id]] +
                        geom_line(
                            data = true_posterior_df,
                            aes(x = value, y = density, color = legend),
                            alpha = alpha_truth, linetype = "solid", linewidth = 5
                        ) +
                        geom_area(alpha = 0.5, fill = legend)
                }
            } else {
                true_posterior_df <- data.frame(value = parameters_truth[[parameter_id]], legend = "True Posterior Distribution")
                if (nrow(parameters_truth) == 1) {
                    plots$parameters[[parameter_id]] <- plots$parameters[[parameter_id]] +
                        geom_vline(data = true_posterior_df, aes(xintercept = value), linetype = "solid", linewidth = 5)
                } else {
                    if (plot_truth_hist) {
                        plots$parameters[[parameter_id]] <- plots$parameters[[parameter_id]] +
                            geom_histogram(
                                data = true_posterior_df,
                                aes(x = value, y = ..density.., fill = legend, color = legend),
                                breaks = breaks, alpha = alpha_truth
                            )
                    } else {
                        plots$parameters[[parameter_id]] <- plots$parameters[[parameter_id]] +
                            geom_density(
                                data = true_posterior_df,
                                aes(x = value, fill = legend, color = legend),
                                alpha = alpha_truth, linewidth = 2
                            )
                    }
                }
            }
        }
    }
    #---Plot statistic value from target data
    if (plot_statistics & new_plot) {
        for (statistic_id in statistics_labels$ID) {
            data_df <- data.frame(value = abc_results$statistics_target[[statistic_id]], legend = "Data statistic")
            plots$statistics[[statistic_id]] <- plots$statistics[[statistic_id]] +
                geom_vline(data = data_df, aes(xintercept = value), linetype = "solid", linewidth = 5)
        }
    }
    #---Extract final posterior distributions
    if (method == "smcrf-single-param") {
        nIterations <- abc_results[["nIterations"]]
        if (nIterations == 1) {
            legend_label <- "ABC-RF"
        } else {
            legend_label <- "ABC-SMC-RF"
        }
        parameters_values <- abc_results[[paste0("Iteration_", nIterations + 1)]]$parameters_unperturbed
        if (plot_statistics) statistics_values <- abc_results[[paste0("Iteration_", nIterations + 1)]]$statistics
    } else if (method == "smcrf-multi-param") {
        nIterations <- abc_results[["nIterations"]]
        if (nIterations == 1) {
            legend_label <- "ABC-DRF"
        } else {
            legend_label <- "ABC-SMC-DRF"
        }
        parameters_values <- abc_results[[paste0("Iteration_", nIterations + 1)]]$parameters_unperturbed
        if (plot_statistics) statistics_values <- abc_results[[paste0("Iteration_", nIterations + 1)]]$statistics
    } else if (method == "abc-rejection") {
        legend_label <- "ABC-REJ"
        parameters_values <- abc_results[["Iteration_2"]]$parameters
        if (plot_statistics) statistics_values <- abc_results[["Iteration_2"]]$statistics
    } else if (method == "abc-smc") {
        legend_label <- "ABC-SMC"
        parameters_values <- abc_results[["Iteration_2"]]$parameters
        if (plot_statistics) statistics_values <- abc_results[["Iteration_2"]]$statistics
    } else if (method == "abc-mcmc") {
        legend_label <- "ABC-MCMC"
        parameters_values <- abc_results[["Iteration_2"]]$parameters
        if (plot_statistics) statistics_values <- abc_results[["Iteration_2"]]$statistics
    } else if (method == "mcmc") {
        legend_label <- "MCMC"
        parameters_values <- abc_results$selected_params
        if (plot_statistics) statistics_values <- abc_results$statistics
    }
    #---Plot marginal distributions for each parameter
    for (parameter_id in parameters_labels$parameter) {
        posterior_df <- data.frame(
            value = parameters_values[[parameter_id]],
            legend = legend_label
        )
        if (plot_hist) {
            plots$parameters[[parameter_id]] <- plots$parameters[[parameter_id]] +
                geom_histogram(
                    data = posterior_df,
                    aes(x = value, y = ..density.., fill = legend, color = legend),
                    breaks = breaks, alpha = alpha
                )
        } else if (plot_hist_point) {
            plots$parameters[[parameter_id]] <- plots$parameters[[parameter_id]] +
                stat_bin(data = posterior_df, aes(x = value, y = ..density.., fill = legend, color = legend), geom = "point", bins = 30, position = "identity", size = 10) +
                stat_bin(data = posterior_df, aes(x = value, y = ..density.., fill = legend, color = legend), geom = "line", bins = 30, position = "identity", size = 3, show.legend = FALSE)
        } else {
            rows_remove <- which((posterior_df$value < min(breaks)) | (posterior_df$value > max(breaks)))
            if (length(rows_remove) > 0) posterior_df <- posterior_df[-rows_remove, ]
            tmp <- hist(posterior_df$value, breaks = breaks, plot = FALSE)
            tmp <- data.frame(mids = tmp$mids, density = tmp$density, legend_label = legend_label)
            plots$parameters[[parameter_id]] <- plots$parameters[[parameter_id]] +
                geom_line(
                    data = tmp,
                    aes(x = mids, y = density, color = legend_label, fill = legend_label),
                    size = 3
                )
            # geom_density(data = posterior_df[sample(1:nrow(posterior_df), size = 600, replace = T), ], aes(x = value, fill = legend, color = legend), alpha = 0, linewidth = 3)
        }
        if (!is.null(xlimit)) {
            plots$parameters[[parameter_id]] <- plots$parameters[[parameter_id]] +
                xlim(xlimit$min[which(xlimit$parameter == parameter_id)], xlimit$max[which(xlimit$parameter == parameter_id)])
        }
    }
    if (plot_statistics) {
        for (statistic_id in statistics_labels$ID) {
            posterior_df <- data.frame(
                value = statistics_values[[statistic_id]],
                legend = legend_label
            )
            if (plot_hist) {
                plots$statistics[[statistic_id]] <- plots$statistics[[statistic_id]] +
                    geom_histogram(data = posterior_df, aes(x = value, y = ..density.., fill = legend, color = legend), alpha = alpha)
            } else if (plot_hist_point) {
                plots$statistics[[statistic_id]] <- plots$statistics[[statistic_id]] +
                    stat_bin(data = posterior_df, aes(x = value, y = ..density.., color = legend), geom = "point", bins = 30, position = "identity", size = 10) +
                    stat_bin(data = posterior_df, aes(x = value, y = ..density.., color = legend), geom = "line", bins = 30, position = "identity", size = 3)
            } else {
                plots$statistics[[statistic_id]] <- plots$statistics[[statistic_id]] +
                    geom_density(data = posterior_df, aes(x = value, fill = legend, color = legend), alpha = alpha, linewidth = 2)
            }
        }
    }
    #---Add label for parameters and statistics
    if (new_plot == TRUE) {
        if ("label" %in% colnames(parameters_labels)) {
            for (parameter_id in parameters_labels$parameter) {
                plots$parameters[[parameter_id]] <- plots$parameters[[parameter_id]] +
                    labs(x = eval(parse(text = parameters_labels$label[which(parameters_labels$parameter == parameter_id)])))
            }
        } else {
            for (parameter_id in parameters_labels$parameter) {
                plots$parameters[[parameter_id]] <- plots$parameters[[parameter_id]] +
                    labs(x = parameter_id)
            }
        }
        if (plot_statistics) {
            if ("label" %in% colnames(statistics_labels)) {
                for (statistic_id in statistics_labels$ID) {
                    plots$statistics[[statistic_id]] <- plots$statistics[[statistic_id]] +
                        labs(x = eval(parse(text = statistics_labels$label[which(statistics_labels$ID == statistic_id)])))
                }
            } else {
                for (statistic_id in statistics_labels$ID) {
                    plots$statistics[[statistic_id]] <- plots$statistics[[statistic_id]] +
                        labs(x = statistic_id)
                }
            }
        }
    }
    #---Beautify plots
    if (new_plot) {
        for (parameter_id in parameters_labels$parameter) {
            plots$parameters[[parameter_id]] <- plots$parameters[[parameter_id]] +
                scale_fill_manual(values = color_scheme, name = "", breaks = legend_order) +
                scale_color_manual(values = color_scheme, name = "", breaks = legend_order) +
                guides(fill = guide_legend(nrow = 1, keywidth = 2.5, keyheight = 1)) +
                theme(
                    text = element_text(size = 50),
                    panel.background = element_rect(fill = "white", colour = "white"),
                    panel.grid.major = element_line(colour = "white"),
                    panel.grid.minor = element_line(colour = "white"),
                    legend.position = "top",
                    legend.justification = c(0, 0.5)
                )
        }
        if (plot_statistics) {
            for (statistic_id in statistics_labels$ID) {
                plots$statistics[[statistic_id]] <- plots$statistics[[statistic_id]] +
                    scale_fill_manual(values = color_scheme, name = "", breaks = legend_order) +
                    scale_color_manual(values = color_scheme, name = "", breaks = legend_order) +
                    guides(fill = guide_legend(nrow = 1, keywidth = 2.5, keyheight = 1)) +
                    theme(
                        text = element_text(size = 50),
                        panel.background = element_rect(fill = "white", colour = "white"),
                        panel.grid.major = element_line(colour = "white"),
                        panel.grid.minor = element_line(colour = "white"),
                        legend.position = "top",
                        legend.justification = c(0, 0.5)
                    )
            }
        }
    }
    #---Print marginal distribution plots
    for (parameter_id in parameters_labels$parameter) {
        file_name <- paste0("comparison-marginal-parameter=", parameter_id, ".png")
        png(file_name, res = 150, width = 30, height = 15, units = "in", pointsize = 12)
        print(plots$parameters[[parameter_id]])
        dev.off()
    }
    if (plot_statistics) {
        for (statistic_id in statistics_labels$ID) {
            file_name <- paste0("comparison-marginal-statistic=", statistic_id, ".png")
            png(file_name, res = 150, width = 30, height = 15, units = "in", pointsize = 12)
            print(plots$statistics[[statistic_id]])
            dev.off()
        }
    }
    return(plots)
}
#' Plot joint distribution(s) of each iteration from ABC-SMC-(D)RF result
#'
#' \code{\link{plot_smcrf_joint}} plots the joint distribution(s) of two parameters for each iteration from an
#' Approximate Bayesian Computation sequential Monte Carlo via random forest result.
#'
#' @param smcrf_results An ABC-SMC-(D)RF result containing the inference distributions of parameters from each iteration.
#' @param parameters_truth A dataframe containing true values of parameters from the ground-truth distributions.
#' If provided, the function will plot the true values or distributions of parameters.
#' @param parameters_labels A dataframe containing labels in the plots for corresponding parameters.
#' If provided, parameter labels will be exhibited on the plots' axes.
#' @param lims A dataframe containing the maximum and minimum bounds for parameters.
#' If provided, x-axis and y-axis will be scaled by them.
#' @param nBins Number of contour bins shown in the plot. Default is 5.
#' @seealso
#' \code{\link{smcrf}}
#' @export
#' @examples
#' lims <- data.frame(
#'     parameter = c("theta", "mu"),
#'     min = c(0, 0),
#'     max = c(10, 10)
#' )
#' parameters_truth <- data.frame(
#'          theta = rnorm(10000, mean = 2, sd = 1),
#'          mu = rnorm(10000, mean = 1, sd = 2)
#')
plot_smcrf_joint <- function(smcrf_results,
                             parameters_truth = NULL,
                             parameters_labels = NULL,
                             lims = NULL,
                             nBins = 5) {
    library(MASS)
    library(dplyr)
    library(ggplot2)
    nIterations <- smcrf_results[["nIterations"]]
    method <- smcrf_results[["method"]]
    if (is.null(parameters_labels)) parameters_labels <- smcrf_results[["parameters_labels"]]
    if (nrow(parameters_labels) != 2) stop("ERROR: plot_smcrf_joint only works for two parameters. Please check parameters_labels.")
    #---Set up color scheme for plotting
    color_scheme <- c(
        "Prior Distribution" = "#E5E5E5",
        setNames(rev(rainbow(nIterations)), paste0("Iter. ", 1:nIterations))
    )
    #---Set up legend order for plotting
    legend_order <- c("True Posterior Distribution", "Prior Distribution", paste0("Iter. ", 1:nIterations))
    #---Begin plot
    p <- ggplot()
    #---Plot True Posterior distribution (if provided)
    if (!is.null(parameters_truth)) {
        true_posterior_df <- data.frame(
            x = parameters_truth[[parameters_labels$parameter[1]]],
            y = parameters_truth[[parameters_labels$parameter[2]]]
        )
        # if (!is.null(lims)) true_posterior_df <- apply_lims(true_posterior_df)
        p <- p + geom_density_2d_filled(data = true_posterior_df, aes(x = x, y = y), show.legend = FALSE)
    }
    #---Plot prior distribution
    prior_df <- data.frame(
        x = smcrf_results[["Iteration_1"]]$parameters_unperturbed[[parameters_labels$parameter[1]]],
        y = smcrf_results[["Iteration_1"]]$parameters_unperturbed[[parameters_labels$parameter[2]]],
        legend = "Prior Distribution"
    )
    # if (!is.null(lims)) prior_df <- apply_lims(prior_df)
    p <- p + geom_density_2d(data = prior_df, aes(x = x, y = y, color = legend), linewidth = 3, bins = nBins)
    #---Plot posterior distribution for each iteration
    for (iteration in 1:nIterations) {
        posterior_df <- data.frame(
            x = smcrf_results[[paste0("Iteration_", iteration + 1)]]$parameters_unperturbed[[parameters_labels$parameter[1]]],
            y = smcrf_results[[paste0("Iteration_", iteration + 1)]]$parameters_unperturbed[[parameters_labels$parameter[2]]],
            legend = paste0("Iter. ", iteration)
        )
        p <- p + geom_density_2d(data = posterior_df, aes(x = x, y = y, color = legend), linewidth = 3, bins = nBins)
    }
    #---Add label for parameter
    if ("label" %in% colnames(parameters_labels)) {
        p <- p + labs(x = eval(parse(text = parameters_labels$label[1])), y = eval(parse(text = parameters_labels$label[2])))
    } else {
        p <- p + labs(x = parameter_id[1], y = parameter_id[2])
    }
    #---Beautify plot
    p <- p +
        scale_color_manual(values = color_scheme, name = "", breaks = legend_order) +
        guides(color = guide_legend(nrow = 1, keywidth = 5, override.aes = list(lwd = 10))) +
        theme(
            text = element_text(size = 50),
            panel.background = element_rect(fill = "white", colour = "white"),
            panel.grid.major = element_line(colour = "white"),
            panel.grid.minor = element_line(colour = "white"),
            legend.position = "top",
            legend.justification = c(0, 0.5)
        )
    #---Plot within the limits if provided
    if (!is.null(lims)) {
        p <- p +
            xlim(c(lims$min[which(lims$parameter == parameters_labels$parameter[1])], lims$max[which(lims$parameter == parameters_labels$parameter[1])])) +
            ylim(c(lims$min[which(lims$parameter == parameters_labels$parameter[2])], lims$max[which(lims$parameter == parameters_labels$parameter[2])]))
    }
    #---Print joint distribution plot
    file_name <- paste0(smcrf_results[["method"]], "-joint-parameters=", parameters_labels$parameter[1], "-vs-", parameters_labels$parameter[2], ".png")
    png(file_name, res = 150, width = 30, height = 31, units = "in", pointsize = 12)
    print(p)
    dev.off()
}
#' Plot and compare joint posterior distribution(s) from ABC-SMC-(D)RF result
#'
#' \code{\link{plot_compare_joint}} plots the joint posterior distribution(s) for the provided ABC-SMC-(D)RF result.
#' It can also compare the joint posterior distributions for the provided ABC-SMC-(D)RF result with several ABC-SMC-(D)RF plots results.
#'
#' @param plots An existed ABC-SMC-(D)RF joint plots result.
#' If provided, \code{\link{plot_compare_joint}} will plot the new ABC-SMC-(D)RF result and compare it with provided plots results.
#' If plots = NULL, \code{\link{plot_compare_joint}} will make a new plot for the ABC-SMC-(D)RF result.
#' @param abc_results An ABC-SMC-(D)RF result.
#' Will be plotted out by the function.
#' @param parameters_truth A dataframe containing true values of parameters from the ground-truth distributions.
#' If provided, the function will plot the true values or distributions of parameters.
#' @param parameters_labels A dataframe containing labels in the plots for corresponding parameters.
#' If provided, parameter labels will be exhibited on the plots' axes.
#' @param lims A dataframe containing the maximum and minimum bounds for parameters.
#' If provided, x-axis and y-axis will be scaled by them.
#' @param nBins Number of contour bins shown in the plot. Default is 5.
#' @return A list of \code{ggplot2} objects containing the joint plots results of posterior distributions.
#' The user can use the function to compare ABC-SMC-(D)RF joint plots with the joint posterior distribution(s)
#' of other ABC-SMC-(D)RF result(s).
#'
#' @export
plot_compare_joint <- function(plots = NULL,
                               abc_results,
                               parameters_truth = NULL,
                               parameters_labels = NULL,
                               lims = NULL,
                               nBins = 4) {
    library(MASS)
    library(dplyr)
    library(RColorBrewer)
    library(ggplot2)
    method <- abc_results[["method"]]

    if (is.null(parameters_labels)) parameters_labels <- abc_results[["parameters_labels"]]
    if (nrow(parameters_labels) != 2) stop("ERROR: plot_compare_joint only works for two parameters. Please check parameters_labels.")
    #---Set up color scheme for plotting
    color_scheme <- c(
        "ABC-REJ" = "#20854E",
        "ABC-RF" = "#7876B1",
        "ABC-DRF" = "#0072B5",
        "MCMC" = "#E18727",
        "ABC-MCMC" = "#E18727",
        "ABC-SMC" = "#E18727",
        "ABC-SMC-RF" = "#BC3C29",
        "ABC-SMC-DRF" = "#BC3C29"
    )
    #---Set up legend order for plotting
    legend_order <- c(
        "ABC-REJ",
        "MCMC",
        "ABC-MCMC",
        "ABC-SMC",
        "ABC-RF",
        "ABC-DRF",
        "ABC-SMC-RF",
        "ABC-SMC-DRF"
    )
    #---Begin plots
    if (is.null(plots)) {
        plots <- ggplot()
        new_plot <- TRUE
    } else {
        new_plot <- FALSE
    }
    #---Plot True Posterior parameter distributions (if provided)
    if (!is.null(parameters_truth)) {
        true_posterior_df <- data.frame(
            x = parameters_truth[[parameters_labels$parameter[1]]],
            y = parameters_truth[[parameters_labels$parameter[2]]]
        )
        plots <- plots + geom_density_2d_filled(data = true_posterior_df, aes(x = x, y = y), show.legend = FALSE)
    }
    #---Extract final posterior distributions
    if (method == "smcrf-single-param") {
        nIterations <- abc_results[["nIterations"]]
        if (nIterations == 1) {
            legend_label <- "ABC-RF"
        } else {
            legend_label <- "ABC-SMC-RF"
        }
        posterior_df <- data.frame(
            x = abc_results[[paste0("Iteration_", nIterations + 1)]]$parameters_unperturbed[[parameters_labels$parameter[1]]],
            y = abc_results[[paste0("Iteration_", nIterations + 1)]]$parameters_unperturbed[[parameters_labels$parameter[2]]],
            legend = legend_label
        )
    } else if (method == "smcrf-multi-param") {
        nIterations <- abc_results[["nIterations"]]
        if (nIterations == 1) {
            legend_label <- "ABC-DRF"
        } else {
            legend_label <- "ABC-SMC-DRF"
        }
        posterior_df <- data.frame(
            x = abc_results[[paste0("Iteration_", nIterations + 1)]]$parameters_unperturbed[[parameters_labels$parameter[1]]],
            y = abc_results[[paste0("Iteration_", nIterations + 1)]]$parameters_unperturbed[[parameters_labels$parameter[2]]],
            legend = legend_label
        )
    } else if (method == "abc-rejection") {
        legend_label <- "ABC-REJ"
        posterior_df <- data.frame(
            x = abc_results[["Iteration_2"]]$parameters_unperturbed[[parameters_labels$parameter[1]]],
            y = abc_results[["Iteration_2"]]$parameters_unperturbed[[parameters_labels$parameter[2]]],
            legend = legend_label
        )
    } else if (method == "abc-smc") {
        legend_label <- "ABC-SMC"
        posterior_df <- data.frame(
            x = abc_results[["Iteration_2"]]$parameters_unperturbed[[parameters_labels$parameter[1]]],
            y = abc_results[["Iteration_2"]]$parameters_unperturbed[[parameters_labels$parameter[2]]],
            legend = legend_label
        )
    } else if (method == "abc-mcmc") {
        legend_label <- "ABC-MCMC"
        posterior_df <- data.frame(
            x = abc_results[["Iteration_2"]]$parameters_unperturbed[[parameters_labels$parameter[1]]],
            y = abc_results[["Iteration_2"]]$parameters_unperturbed[[parameters_labels$parameter[2]]],
            legend = legend_label
        )
    } else if (method == "mcmc") {
        legend_label <- "MCMC"
        posterior_df <- data.frame(
            x = abc_results$selected_params[[parameters_labels$parameter[1]]],
            y = abc_results$selected_params[[parameters_labels$parameter[2]]],
            legend = legend_label
        )
    } else if (method == "true-joint") {
        legend_label <- "TRUE-JOINT"
        posterior_df <- data.frame(
            x = abc_results$parameters_truth[[parameters_labels$parameter[1]]],
            y = abc_results$parameters_truth[[parameters_labels$parameter[2]]],
            legend = legend_label
        )
    }
    #---Plot joint distribution
    if (method == "true-joint") {
        plots <- plots +
            geom_density_2d_filled(data = posterior_df, aes(x = x, y = y), show.legend = FALSE) +
            scale_fill_grey()
    } else {
        posterior_df <- posterior_df[sample(1:nrow(posterior_df), size = 400, replace = T), ]
        plots <- plots +
            geom_density_2d(data = posterior_df, aes(x = x, y = y, color = legend), linewidth = 3, bins = nBins)
    }
    #---Add label for parameter
    if (new_plot == TRUE) {
        if ("label" %in% colnames(parameters_labels)) {
            plots <- plots + labs(x = eval(parse(text = parameters_labels$label[1])), y = eval(parse(text = parameters_labels$label[2])))
        } else {
            plots <- plots + labs(x = parameters_labels$parameter[1], y = parameters_labels$parameter[2])
        }
    }
    #---Beautify plot
    plots <- plots +
        scale_color_manual(values = color_scheme, name = "", breaks = legend_order) +
        guides(color = guide_legend(nrow = 1, keywidth = 5, override.aes = list(lwd = 10))) +
        theme(
            text = element_text(size = 50),
            panel.background = element_rect(fill = "white", colour = "white"),
            panel.grid.major = element_line(colour = "white"),
            panel.grid.minor = element_line(colour = "white"),
            legend.position = "top",
            legend.justification = c(0, 0.5)
        )
    #---Plot within the limits if provided
    if (!is.null(lims)) {
        plots <- plots +
            xlim(c(lims$min[which(lims$parameter == parameters_labels$parameter[1])], lims$max[which(lims$parameter == parameters_labels$parameter[1])])) +
            ylim(c(lims$min[which(lims$parameter == parameters_labels$parameter[2])], lims$max[which(lims$parameter == parameters_labels$parameter[2])]))
    }
    #---Print joint distribution plot
    file_name <- paste0("comparison-joint-parameters=", parameters_labels$parameter[1], "-vs-", parameters_labels$parameter[2], ".png")
    png(file_name, res = 150, width = 30, height = 31, units = "in", pointsize = 12)
    print(plots)
    dev.off()
    return(plots)
}
#' Plot and compare marginal quantile-quantile plots from ABC-SMC-(D)RF result
#'
#' \code{\link{plot_compare_qqplot}} plots the marginal quantile-quantile plots for inferred parameters and parameters from ground-truth distributions.
#'
#' @param plots An existed ABC-SMC-(D)RF quantile-quantile plots result.
#' If provided, \code{\link{plot_compare_qqplot}} will plot the quantile-qantile plot for inferred parameters in the new ABC-SMC-(D)RF result and compare it with provided quantile-quantile plots result.
#' If plots = NULL, \code{\link{plot_compare_qqplot}} will make a new quantile-quantile plot for inferred parameters in the ABC-SMC-(D)RF result and true parameters.
#' @param abc_results An ABC-SMC-(D)RF result.
#' The function will plot the quantile-quantile plot between the inferred parameters from the ABC-SMC-(D)RF result and true parameters.
#' @param parameters_truth A dataframe containing true values of parameters from the ground-truth distributions.
#' @param parameters_labels A dataframe containing labels in the plots for corresponding parameters.
#' If provided, parameter labels will be exhibited on the plots' axes.
#' @param lims A dataframe containing the maximum and minimum bounds for parameters.
#' If provided, x-axis and y-axis will be scaled by them.
#' @return A list of \code{ggplot2} objects containing the quantile-quantile plots results.
#' The user can use the function to compare ABC-SMC-(D)RF quantile-quantile plots with the quantile-quantile plots
#' of other ABC-SMC-(D)RF result(s).
#' @export
#' @examples
#' #    parameters_truth is optional here!
#' parameters_truth <- data.frame(
#'     theta = rnorm(10000, mean = 2, 1)
#' )
#' #    Dataframe containing the parameter labels
#' parameters_labels <- data.frame(
#'     parameter = c("theta"),
#'     label = c(deparse(expression(theta)))
#' )
#' #    Dataframe containing the x-axis bounds
#' xlimit <- data.frame(
#'     parameter = c("theta"),
#'     min = c(1),
#'     max = c(20)
#' )
plot_compare_qqplot <- function(plots = NULL,
                                abc_results,
                                parameters_truth,
                                parameters_labels = NULL,
                                lims = NULL) {
    library(ggplot2)
    if (is.null(parameters_labels)) parameters_labels <- abc_results[["parameters_labels"]]
    method <- abc_results[["method"]]
    #---Set up color scheme for plotting
    color_scheme <- c(
        "ABC-REJ" = "#20854E",
        "ABC-RF" = "#7876B1",
        "ABC-DRF" = "#0072B5",
        "MCMC" = "#E18727",
        "ABC-MCMC" = "#E18727",
        "ABC-SMC" = "#E18727",
        "ABC-SMC-RF" = "#BC3C29",
        "ABC-SMC-DRF" = "#BC3C29"
    )
    #---Set up legend order for plotting
    legend_order <- c(
        "ABC-REJ",
        "ABC-MCMC",
        "ABC-SMC",
        "ABC-RF",
        "ABC-DRF",
        "MCMC",
        "ABC-SMC-RF",
        "ABC-SMC-DRF"
    )
    #---Begin plots
    if (is.null(plots)) {
        plots <- list()
        for (parameter_id in parameters_labels$parameter) {
            plots$parameters[[parameter_id]] <- ggplot()
        }
        new_plot <- TRUE
    } else {
        new_plot <- FALSE
    }
    #---Extract final posterior distributions
    if (method == "smcrf-single-param") {
        nIterations <- abc_results[["nIterations"]]
        if (nIterations == 1) {
            legend_label <- "ABC-RF"
        } else {
            legend_label <- "ABC-SMC-RF"
        }
        parameters_values <- abc_results[[paste0("Iteration_", nIterations + 1)]]$parameters_unperturbed
    } else if (method == "smcrf-multi-param") {
        nIterations <- abc_results[["nIterations"]]
        if (nIterations == 1) {
            legend_label <- "ABC-DRF"
        } else {
            legend_label <- "ABC-SMC-DRF"
        }
        parameters_values <- abc_results[[paste0("Iteration_", nIterations + 1)]]$parameters_unperturbed
    } else if (method == "abc-rejection") {
        legend_label <- "ABC-REJ"
        parameters_values <- abc_results[["Iteration_2"]]$parameters_unperturbed
    } else if (method == "abc-smc") {
        legend_label <- "ABC-SMC"
        parameters_values <- abc_results[["Iteration_2"]]$parameters_unperturbed
    } else if (method == "abc-mcmc") {
        legend_label <- "ABC-MCMC"
        parameters_values <- abc_results[["Iteration_2"]]$parameters_unperturbed
    } else if (method == "mcmc") {
        legend_label <- "MCMC"
        parameters_values <- abc_results$selected_params
    }

    if (new_plot) {
        for (parameter_id in parameters_labels$parameter) {
            plots$parameters[[parameter_id]] <- plots$parameters[[parameter_id]] +
                geom_abline(intercept = 0, slope = 1, color = "black", size = 1.5)
        }
    }

    for (parameter_id in parameters_labels$parameter) {
        true_para <- qqplot(x = parameters_truth[[parameter_id]], y = parameters_values[[parameter_id]], plot.it = FALSE)$x
        fitted_para <- qqplot(x = parameters_truth[[parameter_id]], y = parameters_values[[parameter_id]], plot.it = FALSE)$y
        qq_data <- data.frame(True_parameter = true_para, Fitted_parameter = fitted_para, legend = legend_label)
        plots$parameters[[parameter_id]] <- plots$parameters[[parameter_id]] +
            geom_point(data = qq_data, aes(x = True_parameter, y = Fitted_parameter, color = legend), size = 10)
    }
    if (new_plot == TRUE) {
        if ("label" %in% colnames(parameters_labels)) {
            for (parameter_id in parameters_labels$parameter) {
                label_for_id <- parameters_labels$label[parameters_labels$parameter == parameter_id]
                label_string <- gsub("expression\\((.*)\\)", "\\1", label_for_id)
                label_expression <- parse(text = label_string)[[1]]
                plots$parameters[[parameter_id]] <- plots$parameters[[parameter_id]] +
                    labs(
                        x = bquote("Likelihood-based " ~ .(label_expression)),
                        y = bquote("Inferred " ~ .(label_expression)),
                    )
            }
        } else {
            for (parameter_id in parameters_labels$parameter) {
                plots$parameters[[parameter_id]] <- plots$parameters[[parameter_id]] +
                    labs(
                        x = paste0(
                            "Likelihood ", parameter_id
                        ),
                        y = paste0(
                            "Inferred ", parameter_id
                        )
                    )
            }
        }
    }
    #---Beautify plots
    if (new_plot) {
        for (parameter_id in parameters_labels$parameter) {
            plots$parameters[[parameter_id]] <- plots$parameters[[parameter_id]] +
                scale_color_manual(values = color_scheme, name = "", breaks = legend_order) +
                guides(color = guide_legend(override.aes = list(size = 10))) +
                theme(
                    text = element_text(size = 50),
                    panel.background = element_rect(fill = "white", colour = "white"),
                    panel.grid.major = element_line(colour = "white"),
                    panel.grid.minor = element_line(colour = "white"),
                    legend.key.size = unit(10, "points"),
                    legend.position = "top",
                    legend.justification = c(0, 0.5)
                )
        }
    }
    if (!is.null(lims)) {
        for (parameter_id in parameters_labels$parameter) {
            plots$parameters[[parameter_id]] <- plots$parameters[[parameter_id]] +
                xlim(c(lims$min[which(lims$parameter == parameter_id)], lims$max[which(lims$parameter == parameter_id)])) +
                ylim(c(lims$min[which(lims$parameter == parameter_id)], lims$max[which(lims$parameter == parameter_id)]))
        }
    }
    #---Print marginal distribution plots
    for (parameter_id in parameters_labels$parameter) {
        file_name <- paste0("qq-plot-compare-parameter=", parameter_id, ".png")
        png(file_name, res = 150, width = 30, height = 15, units = "in", pointsize = 12)
        print(plots$parameters[[parameter_id]])
        dev.off()
    }
    return(plots)
}
