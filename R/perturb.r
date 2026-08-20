#' Independent Beaumont perturbation kernel
#'
#' Builds the default ABC-SMC perturbation kernel used by \code{\link{smcrf}}:
#' an independent truncated normal in each parameter, with variance equal to
#' twice the empirical variance of particles resampled from the previous
#' iteration (Beaumont et al. 2009). This is the same scheme used in DECODE.
#'
#' Parameters listed in \code{parameter_bounds} are sampled from a truncated
#' normal on \code{[min, max]}. Parameters not listed are treated as unbounded
#' and use an ordinary normal kernel. Half-bounded intervals such as
#' \code{[0, Inf)} are allowed.
#'
#' @param parameter_bounds Optional dataframe with columns \code{parameter},
#'   \code{min}, and \code{max}. Each row gives the truncation interval for one
#'   parameter. Missing parameters default to \code{(-Inf, Inf)}.
#' @return A list with functions \code{rperturb} and \code{dperturb} that match
#'   the signatures required by \code{\link{smcrf}}.
#' @seealso \code{\link{smcrf}}
#' @export
#' @examples
#' kernel <- make_beaumont_kernel(
#'     parameter_bounds = data.frame(
#'         parameter = c("theta1", "theta2"),
#'         min = c(-Inf, 0),
#'         max = c(Inf, Inf)
#'     )
#' )
make_beaumont_kernel <- function(parameter_bounds = NULL) {
    parameter_bounds <- normalize_parameter_bounds(parameter_bounds)
    list(
        rperturb = function(parameters_unperturbed, parameters_previous_sampled, iteration) {
            beaumont_rperturb(
                parameters_unperturbed = parameters_unperturbed,
                parameters_previous_sampled = parameters_previous_sampled,
                parameter_bounds = parameter_bounds
            )
        },
        dperturb = function(parameters, parameters_previous, parameters_previous_sampled, iteration, parameter_id = "all") {
            beaumont_dperturb(
                parameters = parameters,
                parameters_previous = parameters_previous,
                parameters_previous_sampled = parameters_previous_sampled,
                parameter_id = parameter_id,
                parameter_bounds = parameter_bounds
            )
        }
    )
}

normalize_parameter_bounds <- function(parameter_bounds) {
    if (is.null(parameter_bounds)) {
        return(NULL)
    }
    if (!is.data.frame(parameter_bounds) ||
        !all(c("parameter", "min", "max") %in% names(parameter_bounds))) {
        stop("parameter_bounds must be a data.frame with columns 'parameter', 'min', and 'max'.")
    }
    bounds <- data.frame(
        parameter = as.character(parameter_bounds$parameter),
        min = as.numeric(parameter_bounds$min),
        max = as.numeric(parameter_bounds$max),
        stringsAsFactors = FALSE
    )
    if (anyDuplicated(bounds$parameter)) {
        stop("parameter_bounds must contain unique parameter names.")
    }
    if (anyNA(bounds$parameter) || any(bounds$parameter == "")) {
        stop("parameter_bounds$parameter must contain non-empty parameter IDs.")
    }
    if (anyNA(bounds$min) || anyNA(bounds$max) || any(bounds$min >= bounds$max)) {
        stop("Each row of parameter_bounds must have min < max (finite or infinite).")
    }
    bounds
}

bounds_for_parameters <- function(parameter_ids, parameter_bounds) {
    mins <- stats::setNames(rep(-Inf, length(parameter_ids)), parameter_ids)
    maxs <- stats::setNames(rep(Inf, length(parameter_ids)), parameter_ids)
    if (is.null(parameter_bounds) || length(parameter_ids) == 0) {
        return(list(min = mins, max = maxs))
    }
    for (i in seq_len(nrow(parameter_bounds))) {
        id <- parameter_bounds$parameter[i]
        if (!id %in% parameter_ids) next
        mins[[id]] <- parameter_bounds$min[i]
        maxs[[id]] <- parameter_bounds$max[i]
    }
    list(min = mins, max = maxs)
}

beaumont_variances <- function(parameters_previous_sampled) {
    sampled <- as.data.frame(parameters_previous_sampled)
    if (ncol(sampled) == 0) {
        return(numeric(0))
    }
    vars <- vapply(sampled, function(x) stats::var(as.numeric(x), na.rm = TRUE), numeric(1))
    vars[is.na(vars)] <- 1e-10
    2 * pmax(vars, 1e-10)
}

is_unbounded_interval <- function(lower, upper) {
    identical(lower, -Inf) && identical(upper, Inf)
}

perturb_ids <- function(parameters, parameter_id) {
    ids <- if (identical(parameter_id, "all")) colnames(parameters) else parameter_id
    ids[ids %in% colnames(parameters)]
}

beaumont_rperturb <- function(parameters_unperturbed, parameters_previous_sampled, parameter_bounds) {
    parameters_perturbed <- as.data.frame(parameters_unperturbed)
    variances <- beaumont_variances(parameters_previous_sampled)
    bounds <- bounds_for_parameters(colnames(parameters_perturbed), parameter_bounds)
    n <- nrow(parameters_perturbed)
    for (id in colnames(parameters_perturbed)) {
        sd <- sqrt(if (is.null(variances[[id]]) || !is.finite(variances[[id]])) 1e-10 else variances[[id]])
        lower <- bounds$min[[id]]
        upper <- bounds$max[[id]]
        if (is_unbounded_interval(lower, upper)) {
            parameters_perturbed[[id]] <- stats::rnorm(
                n = n,
                mean = parameters_perturbed[[id]],
                sd = sd
            )
        } else {
            parameters_perturbed[[id]] <- truncnorm::rtruncnorm(
                n = n,
                a = lower,
                b = upper,
                mean = parameters_perturbed[[id]],
                sd = sd
            )
        }
    }
    parameters_perturbed
}

beaumont_dperturb <- function(parameters, parameters_previous, parameters_previous_sampled, parameter_id, parameter_bounds) {
    parameters <- as.data.frame(parameters)
    parameters_previous <- as.data.frame(parameters_previous)
    variances <- beaumont_variances(parameters_previous_sampled)
    ids <- perturb_ids(parameters, parameter_id)
    bounds <- bounds_for_parameters(ids, parameter_bounds)
    probs <- rep(1, nrow(parameters))
    for (id in ids) {
        sd <- sqrt(if (is.null(variances[[id]]) || !is.finite(variances[[id]])) 1e-10 else variances[[id]])
        lower <- bounds$min[[id]]
        upper <- bounds$max[[id]]
        mean <- parameters_previous[[id]]
        if (is_unbounded_interval(lower, upper)) {
            probs <- probs * stats::dnorm(
                parameters[[id]],
                mean = mean,
                sd = sd
            )
        } else {
            probs <- probs * truncnorm::dtruncnorm(
                parameters[[id]],
                a = lower,
                b = upper,
                mean = mean,
                sd = sd
            )
        }
    }
    probs
}

resolve_perturbation <- function(rperturb, dperturb, parameter_bounds) {
    supplied_r <- !is.null(rperturb)
    supplied_d <- !is.null(dperturb)
    if (xor(supplied_r, supplied_d)) {
        stop("rperturb and dperturb must both be NULL (to use the default Beaumont kernel) or both supplied.")
    }
    if (!supplied_r) {
        return(make_beaumont_kernel(parameter_bounds = parameter_bounds))
    }
    list(rperturb = rperturb, dperturb = dperturb)
}
