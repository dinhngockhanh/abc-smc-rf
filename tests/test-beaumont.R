# Comprehensive tests for the default Beaumont perturbation kernel.
# Run from the package root:
#   Rscript tests/test-beaumont.R

Sys.setenv(R_TESTS = "")
options(warn = 1)

if (!file.exists("R/perturb.r")) {
    args <- commandArgs(trailingOnly = FALSE)
    file_arg <- grep("^--file=", args, value = TRUE)
    if (length(file_arg)) {
        setwd(dirname(dirname(normalizePath(sub("^--file=", "", file_arg)))))
    }
}
if (!file.exists("R/perturb.r") && file.exists(file.path("..", "R", "perturb.r"))) {
    setwd("..")
}

stopifnot(file.exists("R/perturb.r"), file.exists("R/smcrf.r"))
sys.source("R/perturb.r", envir = environment())
sys.source("R/smcrf.r", envir = environment())

n_run <- 0L
n_pass <- 0L
failures <- character()

assert <- function(cond, msg) {
    n_run <<- n_run + 1L
    ok <- isTRUE(cond)
    if (length(cond) > 1) ok <- isTRUE(all(cond))
    if (ok) {
        n_pass <<- n_pass + 1L
        cat("  PASS ", msg, "\n", sep = "")
    } else {
        failures <<- c(failures, msg)
        cat("  FAIL ", msg, "\n", sep = "")
    }
    invisible(ok)
}

assert_equal <- function(a, b, msg, tolerance = sqrt(.Machine$double.eps)) {
    assert(isTRUE(all.equal(a, b, tolerance = tolerance, check.attributes = FALSE)), msg)
}

assert_error <- function(expr, msg, pattern = NULL) {
    err <- tryCatch(
        {
            force(expr)
            NULL
        },
        error = function(e) e
    )
    if (is.null(err)) {
        assert(FALSE, paste0(msg, " (no error raised)"))
        return(invisible(FALSE))
    }
    if (!is.null(pattern) && !grepl(pattern, conditionMessage(err))) {
        assert(FALSE, paste0(msg, " (message: ", conditionMessage(err), ")"))
        return(invisible(FALSE))
    }
    assert(TRUE, msg)
}

section <- function(title) cat("\n== ", title, " ==\n", sep = "")

decode_style_rperturb <- function(parameters_unperturbed, parameters_previous_sampled, bounds) {
    Beaumont_variances <- 2 * pmax(sapply(as.data.frame(parameters_previous_sampled), var), 1e-10)
    parameters_perturbed <- as.data.frame(parameters_unperturbed)
    for (id in names(parameters_perturbed)) {
        lower <- bounds$min[bounds$parameter == id]
        upper <- bounds$max[bounds$parameter == id]
        if (!length(lower)) {
            lower <- -Inf
            upper <- Inf
        }
        if (is.infinite(lower) && is.infinite(upper) && lower < 0) {
            parameters_perturbed[[id]] <- rnorm(
                n = nrow(parameters_perturbed),
                mean = parameters_perturbed[[id]],
                sd = sqrt(Beaumont_variances[[id]])
            )
        } else {
            parameters_perturbed[[id]] <- truncnorm::rtruncnorm(
                n = nrow(parameters_perturbed),
                a = lower,
                b = upper,
                mean = parameters_perturbed[[id]],
                sd = sqrt(Beaumont_variances[[id]])
            )
        }
    }
    parameters_perturbed
}

decode_style_dperturb <- function(parameters, parameters_previous, parameters_previous_sampled, parameter_id, bounds) {
    Beaumont_variances <- 2 * pmax(sapply(as.data.frame(parameters_previous_sampled), var), 1e-10)
    probs <- rep(1, nrow(parameters))
    ids <- if (identical(parameter_id, "all")) names(parameters) else parameter_id
    for (id in ids) {
        if (is.null(parameters[[id]])) next
        lower <- bounds$min[bounds$parameter == id]
        upper <- bounds$max[bounds$parameter == id]
        if (!length(lower)) {
            lower <- -Inf
            upper <- Inf
        }
        if (is.infinite(lower) && is.infinite(upper) && lower < 0) {
            probs <- probs * dnorm(
                parameters[[id]],
                mean = parameters_previous[[id]],
                sd = sqrt(Beaumont_variances[[id]])
            )
        } else {
            probs <- probs * truncnorm::dtruncnorm(
                parameters[[id]],
                a = lower,
                b = upper,
                mean = parameters_previous[[id]],
                sd = sqrt(Beaumont_variances[[id]])
            )
        }
    }
    probs
}

toy_model <- function(parameters, parallel = TRUE) {
    n <- nrow(parameters)
    out <- parameters
    out$s_a <- parameters$a + rnorm(n, 0, 0.05)
    if ("b" %in% names(parameters)) {
        out$s_b <- parameters$b + rnorm(n, 0, 0.05)
    }
    out
}

# -----------------------------------------------------------------------------
section("parameter_bounds validation")

assert(is.null(normalize_parameter_bounds(NULL)), "NULL bounds stay NULL")

ok_bounds <- data.frame(parameter = c("a", "b"), min = c(0, -Inf), max = c(1, Inf))
normed <- normalize_parameter_bounds(ok_bounds)
assert_equal(normed$parameter, c("a", "b"), "normalize keeps parameter names")
assert_equal(normed$min, c(0, -Inf), "normalize keeps min")
assert_equal(normed$max, c(1, Inf), "normalize keeps max")

assert_error(
    normalize_parameter_bounds(data.frame(min = 0, max = 1)),
    "rejects bounds missing parameter column",
    pattern = "parameter"
)
assert_error(
    normalize_parameter_bounds(data.frame(parameter = c("a", "a"), min = c(0, 0), max = c(1, 2))),
    "rejects duplicate parameter names",
    pattern = "unique"
)
assert_error(
    normalize_parameter_bounds(data.frame(parameter = "a", min = 1, max = 1)),
    "rejects min >= max",
    pattern = "min < max"
)
assert_error(
    normalize_parameter_bounds(data.frame(parameter = "a", min = NA_real_, max = 1)),
    "rejects NA bounds",
    pattern = "min < max"
)
assert_error(
    normalize_parameter_bounds(data.frame(parameter = "", min = 0, max = 1)),
    "rejects empty parameter ID",
    pattern = "non-empty"
)

# snapshot: later mutation of the input must not change the kernel
pb <- data.frame(parameter = "x", min = 0, max = 1)
kernel_snap <- make_beaumont_kernel(pb)
pb$max <- 100
set.seed(1)
sampled <- data.frame(x = runif(400, 0, 1))
set.seed(2)
out_snap <- kernel_snap$rperturb(data.frame(x = rep(0.5, 200)), sampled, 2)
assert(all(out_snap$x >= 0 & out_snap$x <= 1 + 1e-10), "kernel snapshots bounds (ignores later mutation)")

# -----------------------------------------------------------------------------
section("Beaumont variance")

sampled <- data.frame(a = c(1, 2, 3, 4, 5), b = c(0, 0, 0, 0, 0))
vars <- beaumont_variances(sampled)
assert_equal(unname(vars["a"]), 2 * var(sampled$a), "variance is twice the empirical variance")
assert_equal(unname(vars["b"]), 2e-10, "zero empirical variance is floored then doubled")
one_row <- data.frame(a = 1)
vars_one <- beaumont_variances(one_row)
assert(is.finite(vars_one[["a"]]) && vars_one[["a"]] >= 2e-10, "NA variance from n=1 is floored")

# -----------------------------------------------------------------------------
section("resolve_perturbation")

assert_error(
    resolve_perturbation(rperturb = identity, dperturb = NULL, parameter_bounds = NULL),
    "error if only rperturb is supplied",
    pattern = "both"
)
assert_error(
    resolve_perturbation(rperturb = NULL, dperturb = identity, parameter_bounds = NULL),
    "error if only dperturb is supplied",
    pattern = "both"
)
default_kernel <- resolve_perturbation(NULL, NULL, NULL)
assert(is.function(default_kernel$rperturb) && is.function(default_kernel$dperturb), "NULL/NULL yields default kernel")
custom <- resolve_perturbation(rperturb = identity, dperturb = abs, parameter_bounds = ok_bounds)
assert(identical(custom$rperturb, identity) && identical(custom$dperturb, abs), "custom pair is passed through")
assert(
    identical(custom$rperturb, identity),
    "custom kernel wins even if parameter_bounds is also supplied"
)

# -----------------------------------------------------------------------------
section("parity with DECODE / hierarchical handwritten kernel")

hier_bounds <- data.frame(
    parameter = c("theta1", "theta2"),
    min = c(-Inf, 0),
    max = c(Inf, Inf)
)
kernel_hier <- make_beaumont_kernel(hier_bounds)
set.seed(10)
sampled_h <- data.frame(theta1 = rnorm(300), theta2 = rexp(300))
unpert_h <- data.frame(theta1 = c(-0.4, 1.2, 0), theta2 = c(0.3, 2, 1.1))
prev_h <- data.frame(theta1 = 0.15, theta2 = 0.9)

set.seed(11)
got_r <- kernel_hier$rperturb(unpert_h, sampled_h, iteration = 2)
set.seed(11)
exp_r <- decode_style_rperturb(unpert_h, sampled_h, hier_bounds)
assert_equal(got_r, exp_r, "hierarchical rperturb matches handwritten Beaumont")
assert(all(got_r$theta2 > 0), "theta2 samples stay strictly positive")

got_d <- kernel_hier$dperturb(unpert_h, prev_h, sampled_h, 2, "all")
exp_d <- decode_style_dperturb(unpert_h, prev_h, sampled_h, "all", hier_bounds)
assert_equal(got_d, exp_d, "hierarchical dperturb matches handwritten Beaumont")
got_d1 <- kernel_hier$dperturb(unpert_h, prev_h, sampled_h, 2, "theta1")
exp_d1 <- decode_style_dperturb(unpert_h, prev_h, sampled_h, "theta1", hier_bounds)
assert_equal(got_d1, exp_d1, "single-parameter dperturb (theta1) matches reference")
got_d2 <- kernel_hier$dperturb(unpert_h, prev_h, sampled_h, 2, "theta2")
exp_d2 <- decode_style_dperturb(unpert_h, prev_h, sampled_h, "theta2", hier_bounds)
assert_equal(got_d2, exp_d2, "single-parameter dperturb (theta2) matches reference")
assert_equal(got_d, got_d1 * got_d2, "joint density is the product of independent margins")

outside <- data.frame(theta1 = 0, theta2 = -0.2)
assert_equal(
    kernel_hier$dperturb(outside, prev_h, sampled_h, 2, "theta2"),
    0,
    "truncated density is 0 outside [0, Inf)"
)

decode_bounds <- data.frame(
    parameter = c("alpha", "omega_inference_A_0", "p_1", "omega_inference_A_1"),
    min = c(0.9, 0, 0.01, 0),
    max = c(5, 200, 1, 200)
)
kernel_dec <- make_beaumont_kernel(decode_bounds)
set.seed(20)
sampled_d <- data.frame(
    alpha = runif(250, 0.9, 5),
    omega_inference_A_0 = runif(250, 0, 200),
    p_1 = runif(250, 0.01, 1),
    omega_inference_A_1 = runif(250, 0, 200)
)
unpert_d <- sampled_d[1:5, ]
prev_d <- sampled_d[6, , drop = FALSE]
set.seed(21)
got_dec <- kernel_dec$rperturb(unpert_d, sampled_d, 3)
set.seed(21)
exp_dec <- decode_style_rperturb(unpert_d, sampled_d, decode_bounds)
assert_equal(got_dec, exp_dec, "DECODE-style 4-parameter rperturb matches")
assert(
    all(got_dec$alpha >= 0.9 & got_dec$alpha <= 5) &&
        all(got_dec$p_1 >= 0.01 & got_dec$p_1 <= 1) &&
        all(got_dec$omega_inference_A_0 >= 0 & got_dec$omega_inference_A_0 <= 200),
    "DECODE-style samples remain inside box bounds"
)
assert_equal(
    kernel_dec$dperturb(unpert_d, prev_d, sampled_d, 3, "all"),
    decode_style_dperturb(unpert_d, prev_d, sampled_d, "all", decode_bounds),
    "DECODE-style dperturb matches"
)

# extra bounds for unused names are ignored; missing names stay unbounded
partial <- make_beaumont_kernel(data.frame(parameter = c("b", "unused"), min = c(0, 0), max = c(1, 1)))
set.seed(3)
sampled_ab <- data.frame(a = rnorm(100), b = runif(100))
set.seed(4)
out_ab <- partial$rperturb(data.frame(a = 0, b = 0.5), sampled_ab, 2)
assert(is.finite(out_ab$a), "unlisted parameter a is treated as unbounded")
assert(out_ab$b >= 0 && out_ab$b <= 1, "listed parameter b is truncated")

# iteration is unused (kernel is not annealed); same seed => same sample
set.seed(99)
r_i2 <- kernel_hier$rperturb(unpert_h, sampled_h, 2)
set.seed(99)
r_i9 <- kernel_hier$rperturb(unpert_h, sampled_h, 9)
assert_equal(r_i2, r_i9, "iteration argument does not change Beaumont samples")

# -----------------------------------------------------------------------------
section("smcrf argument checks")

rprior_ab <- function(Nparameters) {
    data.frame(a = runif(Nparameters, 0, 1), b = runif(Nparameters, 0, 1))
}
dprior_ab <- function(parameters, parameter_id = "all") {
    probs <- rep(1, nrow(parameters))
    if (parameter_id %in% c("all", "a")) probs <- probs * dunif(parameters$a, 0, 1)
    if (parameter_id %in% c("all", "b")) probs <- probs * dunif(parameters$b, 0, 1)
    probs
}
stats_ab <- data.frame(s_a = 0.4, s_b = 0.6)
bounds_ab <- data.frame(parameter = c("a", "b"), min = c(0, 0), max = c(1, 1))

assert_error(
    smcrf(
        method = "nope",
        statistics_target = stats_ab,
        model = toy_model,
        rprior = rprior_ab,
        dprior = dprior_ab,
        nParticles = 10
    ),
    "invalid method is rejected",
    pattern = "Invalid method"
)
assert_error(
    smcrf(
        method = "smcrf-multi-param",
        statistics_target = stats_ab,
        statistics_selection = data.frame(s_a = 1),
        model = toy_model,
        rprior = rprior_ab,
        dprior = dprior_ab,
        nParticles = 10
    ),
    "statistics_selection is rejected for multi-param",
    pattern = "statistics_selection"
)
assert_error(
    smcrf(
        statistics_target = stats_ab,
        model = toy_model,
        rprior = rprior_ab,
        dprior = dprior_ab,
        rperturb = function(...) NULL,
        nParticles = 10
    ),
    "smcrf errors when only rperturb is supplied",
    pattern = "both"
)

# -----------------------------------------------------------------------------
section("smcrf integration: default Beaumont vs custom identity kernel")

npart <- c(80, 80)

run_smcrf <- function(...) {
    args <- list(
        statistics_target = stats_ab,
        model = toy_model,
        rprior = rprior_ab,
        dprior = dprior_ab,
        nParticles = npart,
        final_sample = FALSE,
        verbose = FALSE,
        save_model = FALSE,
        parallel = FALSE
    )
    extra <- list(...)
    args[names(extra)] <- extra
    suppressWarnings(do.call(smcrf, args))
}

set.seed(101)
res_default_single <- run_smcrf(
    method = "smcrf-single-param",
    parameter_bounds = bounds_ab,
    ntree = 20
)
assert(identical(res_default_single$method, "smcrf-single-param"), "single-param run completes with default kernel")
assert_equal(res_default_single$nIterations, 2, "single-param records two iterations")
assert(
    isTRUE(all.equal(res_default_single$Iteration_1$parameters, res_default_single$Iteration_1$parameters_unperturbed)),
    "iteration 1 does not perturb (parameters == unperturbed)"
)
assert(
    !isTRUE(all.equal(res_default_single$Iteration_2$parameters, res_default_single$Iteration_2$parameters_unperturbed)),
    "iteration 2 default kernel actually perturbs particles"
)
assert(
    all(res_default_single$Iteration_2$parameters$a >= 0 & res_default_single$Iteration_2$parameters$a <= 1) &&
        all(res_default_single$Iteration_2$parameters$b >= 0 & res_default_single$Iteration_2$parameters$b <= 1),
    "default truncated Beaumont keeps particles inside parameter_bounds"
)
w2 <- res_default_single$Iteration_2$weights
assert(isTRUE(all.equal(colSums(w2), setNames(rep(1, ncol(w2)), names(w2)))), "recalibrated weights sum to 1")

r_calls <- 0L
d_calls <- 0L
identity_rperturb <- function(parameters_unperturbed, parameters_previous_sampled, iteration) {
    r_calls <<- r_calls + 1L
    as.data.frame(parameters_unperturbed)
}
identity_dperturb <- function(parameters, parameters_previous, parameters_previous_sampled, iteration, parameter_id = "all") {
    d_calls <<- d_calls + 1L
    rep(1, nrow(as.data.frame(parameters)))
}

set.seed(101)
res_identity <- run_smcrf(
    method = "smcrf-single-param",
    rperturb = identity_rperturb,
    dperturb = identity_dperturb,
    parameter_bounds = bounds_ab,
    ntree = 20
)
assert(r_calls >= 1L, "custom rperturb is called on later iterations")
assert(d_calls >= 1L, "custom dperturb is used for weight recalibration")
assert(
    isTRUE(all.equal(res_identity$Iteration_2$parameters, res_identity$Iteration_2$parameters_unperturbed)),
    "custom identity kernel is used (iteration 2 not Beaumont-perturbed)"
)

# 1-iteration ABC-RF (vignette style) does not need rperturb
set.seed(3)
res_one <- run_smcrf(
    method = "smcrf-single-param",
    nParticles = 80,
    final_sample = FALSE,
    ntree = 20
)
assert_equal(res_one$nIterations, 1, "1-iteration run works without supplying a kernel")

# final_sample should not perturb
r_iters <- integer(0)
count_r <- function(parameters_unperturbed, parameters_previous_sampled, iteration) {
    r_iters <<- c(r_iters, iteration)
    as.data.frame(parameters_unperturbed)
}
count_d <- function(parameters, parameters_previous, parameters_previous_sampled, iteration, parameter_id = "all") {
    rep(1, nrow(as.data.frame(parameters)))
}
set.seed(12)
res_final <- run_smcrf(
    method = "smcrf-single-param",
    nParticles = c(60, 60),
    final_sample = TRUE,
    rperturb = count_r,
    dperturb = count_d,
    ntree = 20
)
assert_equal(res_final$nIterations, 2, "final_sample stores nIterations as length(nParticles)")
assert(!is.null(res_final$Iteration_3$parameters), "final_sample writes Iteration_3")
assert(
    isTRUE(all.equal(res_final$Iteration_3$parameters, res_final$Iteration_3$parameters_unperturbed)),
    "final sampling iteration does not perturb"
)
assert(length(r_iters) >= 1L && all(r_iters == 2L), "rperturb is used only in iteration 2, not the final sample")

set.seed(202)
res_default_multi <- run_smcrf(
    method = "smcrf-multi-param",
    nParticles = c(200, 200),
    parameter_bounds = bounds_ab,
    num.trees = 50
)
assert(identical(res_default_multi$method, "smcrf-multi-param"), "multi-param run completes with default kernel")
assert(
    !isTRUE(all.equal(res_default_multi$Iteration_2$parameters, res_default_multi$Iteration_2$parameters_unperturbed)),
    "multi-param iteration 2 is perturbed by default Beaumont"
)
assert(
    all(res_default_multi$Iteration_2$parameters$a >= 0 & res_default_multi$Iteration_2$parameters$a <= 1),
    "multi-param default kernel respects bounds"
)
w_multi <- res_default_multi$Iteration_2$weights
assert(isTRUE(all.equal(sum(w_multi[, 1]), 1)), "multi-param weights sum to 1")

# continue from previous results with default kernel
set.seed(303)
res_cont <- run_smcrf(
    method = "smcrf-multi-param",
    smcrf_results = res_default_multi,
    nParticles = 200,
    parameter_bounds = bounds_ab,
    final_sample = FALSE,
    num.trees = 50
)
assert_equal(res_cont$nIterations, 3, "continuation adds iterations onto existing results")
assert(!is.null(res_cont$Iteration_3$weights), "continuation writes a new weighted iteration")

# -----------------------------------------------------------------------------
section("summary")
cat(sprintf("\n%d / %d checks passed\n", n_pass, n_run))
if (length(failures)) {
    cat("Failures:\n")
    cat(paste0("  - ", failures, collapse = "\n"), "\n")
    quit(status = 1)
}
cat("All checks passed.\n")
