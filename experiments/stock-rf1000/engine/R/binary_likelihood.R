# Unpenalized Bernoulli likelihood supremum with explicit facial reduction.
# Public numerical helper. Exact feasibility and high-precision objectives
# do not certify the separate nominal chi-squared reference.
# A has signed predictor rows: (2*y-1) * cbind(1,X).
.binary_softplus <- function(z) pmax(z, 0) + log1p(exp(-abs(z)))
.binary_entropy <- function(q) {
  a <- numeric(length(q))
  i <- q > 0 & q < 1
  a[i] <- -q[i] * log(q[i]) - (1 - q[i]) * log1p(-q[i])
  sum(a)
}
.binary_basis <- function(A) {
  decomposition <- qr(A, tol = 1e-10, LAPACK = FALSE)
  if (!decomposition$rank) return(matrix(0, nrow(A), 0L))
  qr.Q(decomposition)[, seq_len(decomposition$rank), drop = FALSE]
}
# gmp's native solver can reject a nonsingular matrix with a zero leading
# pivot. Preserve its fast path, but verify the exact residual and fall back to
# exact Gaussian elimination with row pivoting when necessary.
.binary_solve_exact <- function(M, rhs) {
  n <- nrow(M)
  good <- function(x) !inherits(x, 'error') && length(x) == n &&
    all(gmp::`%*%`(M, x) == rhs)
  answer <- tryCatch(gmp::solve.bigq(M, rhs), error = function(e) e)
  if (good(answer)) return(answer)
  original <- M
  original_rhs <- rhs
  for (k in seq_len(n)) {
    available <- which(M[k:n, k] != 0)
    if (!length(available)) stop('Exact pivoted system is singular')
    pivot <- k + available[1L] - 1L
    if (pivot != k) {
      temporary <- M[k, ]
      M[k, ] <- M[pivot, ]
      M[pivot, ] <- temporary
      temporary <- rhs[k]
      rhs[k] <- rhs[pivot]
      rhs[pivot] <- temporary
    }
    if (k < n) for (i in (k + 1L):n) if (M[i, k] != 0) {
      ratio <- M[i, k] / M[k, k]
      M[i, k:n] <- M[i, k:n] - ratio * M[k, k:n]
      rhs[i] <- rhs[i] - ratio * rhs[k]
    }
  }
  answer <- gmp::as.bigq(rep(0, n))
  for (k in n:1L) {
    rest <- if (k < n) sum(M[k, (k + 1L):n] * answer[(k + 1L):n]) else gmp::as.bigq(0)
    answer[k] <- (rhs[k] - rest) / M[k, k]
  }
  if (any(gmp::`%*%`(original, answer) != original_rhs)) stop('Exact pivoted solve residual failed')
  answer
}
# Repair selected coordinates in exact rational arithmetic. Every input double
# is treated as its exact binary rational value; all equations are rechecked.
.binary_exact_null <- function(C, v, interior_pivots = FALSE) {
  q <- gmp::as.bigq(v)
  if (!nrow(C)) return(q)
  r <- qr(C, tol = 1e-10)$rank
  if (!r) {
    if (any(C != 0)) stop('Unresolved exact rank')
    return(q)
  }
  rows <- qr(t(C), tol = 1e-10)$pivot[seq_len(r)]
  pivot_matrix <- C[rows, , drop = FALSE]
  if (interior_pivots) pivot_matrix <- sweep(pivot_matrix, 2, sqrt(pmax(v * (1 - v), 0)), '*')
  cols <- qr(pivot_matrix, LAPACK = TRUE)$pivot[seq_len(r)]
  free <- setdiff(seq_len(ncol(C)), cols)
  Q <- gmp::as.bigq(C)
  rhs <- if (length(free)) -gmp::`%*%`(Q[rows, free, drop = FALSE], q[free]) else gmp::as.bigq(rep(0, r))
  q[cols] <- .binary_solve_exact(Q[rows, cols, drop = FALSE], rhs)
  if (any(gmp::`%*%`(Q, q) != 0)) stop('Exact rational moment/rank verification failed')
  q
}
.binary_exact_objectives <- function(A, eta, q, w) {
  qexact <- .binary_exact_null(t(A), q, interior_pivots = TRUE)
  if (any(qexact <= 0 | qexact >= 1)) stop('Exact repaired dual is not interior')
  fit <- qr(A, tol = 1e-10)
  beta <- qr.coef(fit, eta)
  beta[is.na(beta)] <- 0
  etaexact <- gmp::`%*%`(gmp::as.bigq(A), gmp::as.bigq(beta))
  evaluate <- function(bits) {
    z <- Rmpfr::mpfr(etaexact, precBits = bits)
    qq <- Rmpfr::mpfr(qexact, precBits = bits)
    primal <- sum(log1p(exp(-abs(z))) + (abs(z) - z) / 2)
    dual <- sum(-qq * log(qq) - (1 - qq) * log1p(-qq))
    c(primal = as.numeric(primal), dual = as.numeric(dual))
  }
  low <- evaluate(128L)
  high <- evaluate(256L)
  if (max(abs(low - high)) > 1e-11) stop('High precision objective agreement failed')
  # Explicit numerical rounding envelope. Feasibility is exact; the entropy
  # calculation is high-precision numerical evidence, not interval arithmetic.
  envelope <- 1e-10
  list(upper = high['primal'] + envelope, lower = high['dual'] - envelope,
    primal = unname(high['primal']), dual_value = unname(high['dual']),
    qexact = qexact, beta = beta, envelope = envelope, precision_bits = 256L,
    precision_difference = max(abs(low - high)), exact_moments = TRUE)
}
binary_face <- function(A, feasibility_tol = 1e-7) {
  n <- nrow(A)
  d <- ncol(A)
  # Free v is vplus-vminus. The capped t objective exposes maximal support:
  # any rows positive under any recession direction can simultaneously reach 1.
  fit <- lpSolve::lp('max', c(rep(0, 2 * d), rep(1, n)),
    rbind(cbind(A, -A, -diag(n)), cbind(matrix(0, n, 2 * d), diag(n))),
    c(rep('>=', n), rep('<=', n)), c(rep(0, n), rep(1, n)))
  if (fit$status != 0L) stop('Facial LP failed: status ', fit$status)
  direction <- fit$solution[seq_len(d)] - fit$solution[d + seq_len(d)]
  margin <- drop(A %*% direction)
  t <- tail(fit$solution, n)
  removed <- which(t > .5)
  remaining <- setdiff(seq_len(n), removed)
  if (any(margin < -feasibility_tol) ||
    (length(removed) && min(margin[removed]) < 1 - feasibility_tol) ||
    (length(remaining) && max(abs(margin[remaining])) > feasibility_tol) ||
    abs(fit$objval - length(removed)) > feasibility_tol)
    stop('Facial LP certificate failed numerical checks')
  exact_direction <- .binary_exact_null(A[remaining, , drop = FALSE], direction)
  exact_margin <- gmp::`%*%`(gmp::as.bigq(A), exact_direction)
  if ((length(removed) && any(exact_margin[removed] <= 0)) ||
    (length(remaining) && any(exact_margin[remaining] != 0))) stop('Exact facial certificate failed')
  list(removed = removed, remaining = remaining, direction = direction, margin = margin,
    exact_direction = exact_direction, exact_margin = exact_margin,
    objective = fit$objval, lp_status = fit$status,
    state = if (!length(removed)) 'finite' else if (!length(remaining)) 'strict_separation' else 'quasi_separation')
}
.binary_interior_dual <- function(B, tolerance = 1e-10) {
  n <- nrow(B)
  r <- ncol(B)
  # Positive weights annihilating B prove no further nonconstant recession
  # direction remains. This independently checks maximality of the exposed face.
  fit <- lpSolve::lp('max', c(rep(0, n), 1),
    rbind(cbind(t(B), rep(0, r)), c(rep(1, n), 0), cbind(diag(n), rep(-1, n))),
    c(rep('=', r), '=', rep('>=', n)), c(rep(0, r), 1, rep(0, n)))
  if (fit$status != 0L) stop('Interior-dual LP failed')
  w <- fit$solution[seq_len(n)]
  if (r) w <- w - drop(B %*% crossprod(B, w))
  w <- w / sum(w)
  residual <- if (r) max(abs(crossprod(B, w))) else 0
  if (min(w) <= tolerance || residual > tolerance) stop('No verified positive interior dual vector')
  list(weights = w, residual = residual, minimum = min(w), lp_minimum = tail(fit$solution, 1))
}
.binary_dual_bound <- function(B, q, w) {
  r <- ncol(B)
  projected <- if (r) q - drop(B %*% crossprod(B, q)) else q
  # Mix two nullspace vectors just enough to stay strictly inside [0,1].
  epsilon <- 1e-14
  amount <- 0
  if (any(projected < epsilon)) amount <- max(amount,
    max((epsilon - projected[projected < epsilon]) / (w[projected < epsilon] - projected[projected < epsilon])))
  if (any(projected > 1 - epsilon)) amount <- max(amount,
    max((projected[projected > 1 - epsilon] - (1 - epsilon)) / (projected[projected > 1 - epsilon] - w[projected > 1 - epsilon])))
  amount <- min(1, max(0, amount) * (1 + 1e-8))
  feasible <- (1 - amount) * projected + amount * w
  if (any(feasible < 0 | feasible > 1)) stop('Dual probabilities outside [0,1]')
  list(q = feasible, lower = .binary_entropy(feasible), mixture = amount,
    residual = if (r) max(abs(crossprod(B, feasible))) else 0)
}
# design includes an intercept column; y uses numeric 0/1 coding.
# Return the fitted face, likelihood bounds, gap and feasibility certificates.
binary_loglik_supremum <- function(design, y, gap_tolerance = 1e-8, max_iterations = 150L) {
  if (!is.matrix(design) || !is.numeric(design) || any(!is.finite(design)) ||
    nrow(design) != length(y) || !all(y %in% c(0, 1))) stop('Invalid binary design/response')
  for (package in c('lpSolve',
    'gmp',
    'Rmpfr')) if (!requireNamespace(package,
    quietly = TRUE)) stop(package,
    ' is required')
  A <- design * (2 * y - 1)
  face <- binary_face(A)
  result <- list(face = face, design_rank = qr(design, tol = 1e-10)$rank,
    n = nrow(A), gap_tolerance = gap_tolerance, method = 'facial-unpenalized-logistic-v2')
  if (!length(face$remaining)) {
    result$nll_upper <- result$nll_lower <- result$nll_value <- result$gap <- 0
    result$dual <- rep(0, nrow(A))
    result$dual_residual <- 0
    result$iterations <- 0L
    result$linear_predictor_remaining <- numeric()
    result$optimality <- 'verified_numeric'
    return(result)
  }
  active <- A[face$remaining, , drop = FALSE]
  B <- .binary_basis(active)
  interior <- .binary_interior_dual(B)
  coefficients <- numeric(ncol(B))
  converged <- FALSE
  for (iteration in 0:max_iterations) {
    eta <- drop(B %*% coefficients)
    q <- stats::plogis(-eta)
    upper <- sum(.binary_softplus(-eta))
    dual <- .binary_dual_bound(B, q, interior$weights)
    gap <- upper - dual$lower
    if (gap < -1e-9) stop('Negative primal-dual gap beyond roundoff')
    # Approximate projected entropy is only a progress diagnostic. Acceptance
    # uses exact original-coordinate feasibility and the unchanged MPFR gap.
    if (max(abs(crossprod(B, q))) <= 1e-7) {
      candidate_certificate <- tryCatch(.binary_exact_objectives(active,
        eta,
        q,
        interior$weights),
      error = function(e) NULL)
      if (!is.null(candidate_certificate) && candidate_certificate$upper - candidate_certificate$lower >= 0 &&
        candidate_certificate$upper - candidate_certificate$lower <= gap_tolerance) {
        certificate <- candidate_certificate
        converged <- TRUE
        break
      }
    }
    if (iteration == max_iterations) break
    gradient <- -drop(crossprod(B, q))
    hessian <- crossprod(B, B * as.numeric(q * (1 - q)))
    eigenfit <- eigen(hessian, symmetric = TRUE)
    # Eigenvalue floor stabilizes a Newton direction; no penalty enters objective.
    direction <- -drop(eigenfit$vectors %*% (crossprod(eigenfit$vectors, gradient) / pmax(eigenfit$values, 1e-14)))
    slope <- sum(gradient * direction)
    if (!is.finite(slope) || slope >= 0) {
      direction <- -gradient
      slope <- -sum(gradient^2)
    }
    step <- 1
    accepted <- FALSE
    for (backtrack in 0:60) {
      candidate <- coefficients + step * direction
      delta <- drop(B %*% (step * direction))
      # Exact algebra for the likelihood difference avoids subtracting two
      # nearly equal sums near the optimum. Large steps use stable softplus.
      change <- if (max(abs(delta)) < .5) sum(log1p(q * expm1(-delta))) else
        sum(.binary_softplus(-(eta + delta)) - .binary_softplus(-eta))
      if (is.finite(change) && change <= 1e-4 * step * slope) {
        accepted <- TRUE
        break
      }
      step <- step / 2
    }
    if (!accepted) stop('Monotone line search failed')
    coefficients <- candidate
  }
  if (!converged) stop('Primal-dual tolerance not reached')
  certificate <- .binary_exact_objectives(active, eta, q, interior$weights)
  upper <- certificate$upper
  lower <- certificate$lower
  gap <- upper - lower
  if (gap < 0 || gap > gap_tolerance) stop('Exact-feasible high-precision objective gap failed')
  result$nll_upper <- unname(upper)
  result$nll_lower <- unname(lower)
  result$gap <- unname(gap)
  result$nll_value <- certificate$primal
  result$exact_certificate <- certificate
  result$dual <- numeric(nrow(A))
  result$dual[face$remaining] <- as.numeric(certificate$qexact)
  result$dual_residual <- max(abs(crossprod(A, result$dual)))
  if (result$dual_residual > 1e-8) stop('Original-coordinate dual feasibility failed')
  result$iterations <- iteration
  result$linear_predictor_remaining <- eta
  result$interior_dual <- interior
  result$optimality <- 'exact_feasibility_high_precision_objective'
  result
}
# X excludes the intercept; target contains one-based tested column indices.
# Return nested likelihood differences, nominal degrees of freedom and tail.
binary_lrt_supremum <- function(X, y, target, gap_tolerance = 1e-8) {
  if (!is.matrix(X) || !is.numeric(X) || any(!is.finite(X)) || nrow(X) != length(y) ||
    !length(target) || !is.numeric(target) || anyNA(target) || anyDuplicated(target) ||
    any(target != as.integer(target)) || any(target < 1 | target > ncol(X))) stop('Invalid X or target')
  full <- binary_loglik_supremum(cbind(1, X), y, gap_tolerance)
  reduced <- binary_loglik_supremum(cbind(1, X[, setdiff(seq_len(ncol(X)), target), drop = FALSE]), y, gap_tolerance)
  statistic <- 2 * (reduced$nll_value - full$nll_value)
  lower <- 2 * (reduced$nll_lower - full$nll_upper)
  upper <- 2 * (reduced$nll_upper - full$nll_lower)
  if (statistic < 0 || upper < 0) stop('Negative nested likelihood-supremum difference; do not truncate')
  df <- full$design_rank - reduced$design_rank
  if (df <= 0) stop('Nonpositive original nominal degrees of freedom')
  list(pvalue = stats::pchisq(statistic, df, lower.tail = FALSE), statistic = statistic,
    statistic_lower = lower, statistic_upper = upper, df = df, full = full, reduced = reduced,
    calibration = 'original nominal chi-square; no new boundary calibration')
}
