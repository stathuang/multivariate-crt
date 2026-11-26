# Unpenalized multinomial likelihood with exact feasibility checks.
# Numerical objective certification does not certify the chi-squared reference.
bcr_qsolve_slow <- function(A, b) {
  m <- nrow(A)
  for (j in seq_len(m)) {
    ii <- (j:m)[which(A[j:m, j] != 0)[1L]]
    if (is.na(ii)) stop("Exact rational system singular")
    if (ii != j) {
      tmp <- A[j, ]
      A[j, ] <- A[ii, ]
      A[ii, ] <- tmp
      tmp <- b[j]
      b[j] <- b[ii]
      b[ii] <- tmp
    }
    div <- A[j, j]
    A[j, ] <- A[j, ] / div
    b[j] <- b[j] / div
    for (i in setdiff(seq_len(m), j)) if (A[i, j] != 0) {
      f <- A[i, j]
      A[i, ] <- A[i, ] - f * A[j, ]
      b[i] <- b[i] - f * b[j]
    }
  }
  b
}
bcr_qsolve <- function(A, b) {
  answer <- tryCatch({
    N <- matrix(as.numeric(Rmpfr::mpfr(A, 128)), nrow(A), ncol(A))
    m <- nrow(N)
    perm <- seq_len(m)
    for (j in seq_len(m)) {
      ii <- (j:m)[which.max(abs(N[j:m, j]))]
      if (ii != j) {
        tmp <- N[j, ]
        N[j, ] <- N[ii, ]
        N[ii, ] <- tmp
        tmp <- perm[j]
        perm[j] <- perm[ii]
        perm[ii] <- tmp
      }
      if (!is.finite(N[j, j]) || N[j, j] == 0) stop("Numerical pivot proposal failed")
      if (j < m) for (i in (j + 1L):m) N[i, ] <- N[i, ] - N[i, j] / N[j, j] * N[j, ]
    }
    value <- gmp::solve.bigq(A[perm, , drop = FALSE], b[perm])
    if (!all(gmp::`%*%`(A, value) == b)) stop("Compiled exact residual check failed")
    attr(value, "solve_backend") <- "compiled_gmp_with_numeric_row_permutation"
    value
  }, error = function(e) e)
  if (inherits(answer, "error")) {
    reason <- conditionMessage(answer)
    answer <- bcr_qsolve_slow(A, b)
    if (!all(gmp::`%*%`(A, answer) == b)) stop("Fallback exact residual check failed")
    attr(answer, "solve_backend") <- "exact_R_GaussJordan_fallback"
    attr(answer, "fast_path_error") <- reason
  }
  answer
}
# Multinomial recession constraints use last class as reference. An exposed
# alternative is removed only from that observation's log-sum-exp support.
# X excludes the intercept and y supplies categorical responses.
# Return exposed support, objective gap and exact feasibility evidence.
bcr_multinomial_supremum <- function(X, y, gap_tolerance = 1e-5,
                                     feasibility_tolerance = 1e-7) {
  X <- as.matrix(X)
  if (!is.numeric(X) || any(!is.finite(X)) || length(y) != nrow(X) || anyNA(y)) stop("Invalid X/y")
  y <- droplevels(as.factor(y))
  D <- cbind(intercept = 1, X)
  y <- as.integer(y)
  n <- nrow(D)
  d <- ncol(D)
  K <- max(y)
  q <- d * (K - 1L)
  if (K < 2L) stop("At least two observed classes required")
  A <- matrix(0, n * (K - 1L), q)
  ij <- matrix(0L, nrow(A), 2L)
  r <- 0L
  for (i in seq_len(n)) for (k in setdiff(seq_len(K), y[i])) {
    r <- r + 1L
    ij[r, ] <- c(i, k)
    if (y[i] < K) A[r, (y[i] - 1L) * d + seq_len(d)] <- D[i, ]
    if (k < K) A[r, (k - 1L) * d + seq_len(d)] <- A[r, (k - 1L) * d + seq_len(d)] - D[i, ]
  }
  nr <- nrow(A)
  face <- lpSolve::lp("max", c(rep(0, 2 * q), rep(1, nr)),
    rbind(cbind(A, -A, -diag(nr)), cbind(matrix(0, nr, 2 * q), diag(nr))),
    c(rep(">=", nr), rep("<=", nr)), c(rep(0, nr), rep(1, nr)))
  face_solver <- "lpSolve"
  original_face_status <- face$status
  if (face$status != 0L) {
    V <- CVXR::Variable(q)
    T <- CVXR::Variable(nr)
    proposal <- CVXR::solve(CVXR::Problem(CVXR::Maximize(CVXR::sum_entries(T)),
      list(A %*% V >= T, T >= 0, T <= 1)), solver = "ECOS",
    abstol = 1e-9, reltol = 1e-9, feastol = 1e-9, max_iters = 1000L)
    if (!proposal$status %in% c("optimal", "optimal_inaccurate")) stop("Both face LP solvers failed")
    v <- as.numeric(proposal$getValue(V))
    t <- as.numeric(proposal$getValue(T))
    face <- list(status = 0L, solution = c(pmax(v, 0), pmax(-v, 0), t), objval = sum(t))
    face_solver <- paste0("ECOS_after_lpSolve_", original_face_status)
  }
  direction <- face$solution[seq_len(q)] - face$solution[q + seq_len(q)]
  margin <- as.numeric(A %*% direction)
  if (min(margin) < -1e-7) stop("Face recession certificate violates constraints")
  slack <- face$solution[2 * q + seq_len(nr)]
  removed <- which(slack > .5)
  # LP maximality is a proposal diagnostic. Acceptance below depends on exact
  # recession feasibility and the final primal/dual objective certificate.
  face_objective_discrepancy <- abs(face$objval - length(removed))
  retained_numeric_margin <- max(c(0, abs(margin[setdiff(seq_len(nr), removed)])))
  if (length(removed) && min(margin[removed]) < .99) stop("Face removal certificate failed")
  # Exact rational repair of recession direction on retained constraints.
  retained0 <- setdiff(seq_len(nr), removed)
  Aq <- gmp::as.bigq(A)
  vq <- gmp::as.bigq(direction)
  direction_solve_backend <- "not_required"
  if (length(removed) && length(retained0)) {
    C <- A[retained0, , drop = FALSE]
    rk <- qr(C, tol = 1e-10)$rank
    rr <- qr(t(C), tol = 1e-10)$pivot[seq_len(rk)]
    cc <- qr(C[rr, , drop = FALSE], tol = 1e-10)$pivot[seq_len(rk)]
    if (rk) {
      correction <- bcr_qsolve(gmp::as.bigq(C[rr, cc, drop = FALSE]),
        -gmp::`%*%`(gmp::as.bigq(C[rr, , drop = FALSE]), vq))
      direction_solve_backend <- attr(correction, "solve_backend")
      vq[cc] <- vq[cc] + correction
    }
  }
  exact_margin <- gmp::`%*%`(Aq, vq)
  if (!all(exact_margin >= 0) || (length(removed) && !all(exact_margin[removed] > 0)) ||
    (length(retained0) && !all(exact_margin[retained0] == 0))) stop("Exact face certificate failed")
  support <- matrix(TRUE, n, K)
  if (length(removed)) support[ij[removed, , drop = FALSE]] <- FALSE
  objective <- function(b, gradient = FALSE) {
    eta <- cbind(D %*% matrix(b, d, K - 1L), 0)
    eta[!support] <- -Inf
    m <- apply(eta, 1L, max)
    e <- exp(eta - m)
    p <- e / rowSums(e)
    if (gradient) {
      p[cbind(seq_len(n), y)] <- p[cbind(seq_len(n), y)] - 1
      return(as.numeric(crossprod(D, p[, seq_len(K - 1L), drop = FALSE])))
    }
    sum(m + log(rowSums(e)) - eta[cbind(seq_len(n), y)])
  }
  opt <- optim(rep(0, q), objective, function(b) objective(b, TRUE), method = "BFGS",
    control = list(maxit = 5000L, reltol = 1e-13))
  # Independent entropy dual on the same exposed face. The objective and
  # equality residual are retained; acceptance is numerical, at stated tolerances.
  P <- CVXR::Variable(n, K)
  Y <- matrix(0, n, K)
  Y[cbind(seq_len(n), y)] <- 1
  constraints <- list(P >= 0, P %*% rep(1, K) == 1,
    t(D) %*% (P - Y) == 0)
  for (k in seq_len(K)) if (any(!support[, k])) constraints <- c(constraints, list(P[which(!support[, k]), k] == 0))
  problem <- CVXR::Problem(CVXR::Maximize(CVXR::sum_entries(CVXR::entr(P))), constraints)
  dual <- CVXR::solve(problem, solver = "ECOS", abstol = 1e-9, reltol = 1e-9, feastol = 1e-9,
    max_iters = 1000L)
  pp <- dual$getValue(P)
  residual <- max(abs(rowSums(pp) - 1), abs(crossprod(D, pp - Y)), pmax(-pp, 0),
    if (any(!support)) abs(pp[!support]) else 0)
  approximate_entropy <- -sum(ifelse(pp > 0, pp * log(pmax(pp, .Machine$double.xmin)), 0))
  # Enforce moment equations over exact rationals after eliminating row sums.
  # z contains probabilities of retained nonobserved alternatives; observed
  # probability is 1 minus their row sum. Thus A' z = 0 is the dual constraint.
  retained <- setdiff(seq_len(nr), removed)
  project_dual <- function(weighted = FALSE) {
    B <- t(A[retained, , drop = FALSE])
    raw_z <- pp[ij[retained, , drop = FALSE]]
    z <- if (weighted) pmax(raw_z, 0) else raw_z
    rank <- qr(B, tol = 1e-10)$rank
    rows <- qr(t(B), tol = 1e-10)$pivot[seq_len(rank)]
    if (weighted) {
      weights <- pmax(pmin(z, pp[cbind(ij[retained, 1L], y[ij[retained, 1L]])]), 1e-14)
      pivots <- qr(sweep(B[rows, , drop = FALSE], 2L, weights, "*"), LAPACK = TRUE)$pivot[seq_len(rank)]
    } else pivots <- qr(B[rows, , drop = FALSE], tol = 1e-10)$pivot[seq_len(rank)]
    Bq <- gmp::as.bigq(B)
    zq <- gmp::as.bigq(z)
    projection_backend <- "not_required"
    if (rank) {
      rhs <- -gmp::`%*%`(Bq[rows, , drop = FALSE], zq)
      correction <- bcr_qsolve(gmp::as.bigq(B[rows, pivots, drop = FALSE]), rhs)
      projection_backend <- attr(correction, "solve_backend")
      zq[pivots] <- zq[pivots] + correction
    }
    exact_moments <- all(gmp::`%*%`(Bq, zq) == 0)
    obsq <- gmp::as.bigq(rep(1, n))
    for (j in seq_along(retained)) obsq[ij[retained[j], 1L]] <- obsq[ij[retained[j], 1L]] - zq[j]
    exact_feasible <- exact_moments && all(zq >= 0) && all(obsq >= 0)
    entropy <- NA_real_
    if (exact_feasible) {
      pq <- c(zq, obsq)
      pq <- pq[pq > 0]
      pm <- Rmpfr::mpfr(pq, precBits = 256)
      # MPFR rounding allowance is conservative for this fixed <=627-term sum;
      # this is a recorded numerical envelope, not an interval-arithmetic proof.
      entropy <- as.numeric(-sum(pm * log(pm))) - 1e-10
    }
    list(zq = zq, obsq = obsq, exact_feasible = exact_feasible, entropy = entropy, backend = projection_backend)
  }
  safe_project <- function(weighted) tryCatch(project_dual(weighted), error = function(e)
    list(zq = NULL,
      obsq = NULL,
      exact_feasible = FALSE,
      entropy = NA_real_,
      backend = "error",
      error = conditionMessage(e)))
  candidate <- safe_project(FALSE)
  zq <- candidate$zq
  obsq <- candidate$obsq
  exact_feasible <- candidate$exact_feasible
  entropy <- candidate$entropy
  # Independently evaluate the attained primal objective in MPFR, using the
  # same exact binary-double design/coefficient inputs as the rational checks.
  dm <- Rmpfr::mpfr(D, 256)
  bm <- Rmpfr::mpfr(matrix(opt$par, d, K - 1L), 256)
  primal_mp <- Rmpfr::mpfr(0, 256)
  for (i in seq_len(n)) {
    eta <- Rmpfr::mpfr(rep(0, K), 256)
    for (k in seq_len(K - 1L)) eta[k] <- sum(dm[i, ] * bm[, k])
    keep <- which(support[i, ])
    maximum <- max(eta[keep])
    primal_mp <- primal_mp + maximum + log(sum(exp(eta[keep] - maximum))) - eta[y[i]]
  }
  primal_upper <- as.numeric(primal_mp) + 1e-10
  gap <- primal_upper - entropy
  certified <- exact_feasible && is.finite(gap) && gap >= 0 && gap <= gap_tolerance
  certificate_route <- "original_exact_projection"
  dual_attempts <- list(original = list(exact_feasible = exact_feasible,
    gap = gap,
    backend = candidate$backend,
    error = candidate$error))
  if (!certified) {
    candidate <- safe_project(TRUE)
    candidate_gap <- primal_upper - candidate$entropy
    candidate_ok <- candidate$exact_feasible && is.finite(candidate_gap) && candidate_gap >= 0 && candidate_gap <= gap_tolerance
    dual_attempts$weighted <- list(exact_feasible = candidate$exact_feasible,
      gap = candidate_gap,
      backend = candidate$backend,
      error = candidate$error)
    if (candidate_ok) {
      zq <- candidate$zq
      obsq <- candidate$obsq
      exact_feasible <- candidate$exact_feasible
      entropy <- candidate$entropy
      gap <- candidate_gap
      certified <- TRUE
      certificate_route <- "interior_weighted_exact_projection"
    }
  }
  list(nll = primal_upper,
    double_nll = opt$value,
    certified = certified,
    status = if (certified) "numerically_certified" else "uncertified",
    direction_solve_backend = direction_solve_backend,
    face_objective_discrepancy = face_objective_discrepancy,
    retained_numeric_margin = retained_numeric_margin,
    face_solver = face_solver,
    certificate_route = certificate_route,
    dual_attempts = dual_attempts,
    exact_dual_alternatives = zq,
    exact_observed_probabilities = obsq,
    retained_indices = retained,
    alternative_map = ij,
    dual_entropy = entropy,
    approximate_dual_entropy = approximate_entropy,
    exact_dual_feasible = exact_feasible,
    gap = gap,
    dual_residual = residual,
    optimizer_code = opt$convergence,
    dual_status = dual$status, removed_alternatives = length(removed), support = support,
    direction = direction,
    exact_direction = vq,
    exact_face_certificate = TRUE,
    recession_margin = margin,
    coefficients = opt$par,
    gradient_max = max(abs(objective(opt$par, TRUE))), gap_tolerance = gap_tolerance,
    feasibility_tolerance = feasibility_tolerance)
}
# Compare full and reduced working likelihoods for specified target columns.
# The nominal tail is distinct from numerical objective certification.
bcr_supremum_lrt <- function(X, y, target, ...) {
  if (anyDuplicated(target) || any(target < 1L | target > ncol(X))) stop("Invalid target")
  full <- bcr_multinomial_supremum(X, y, ...)
  reduced <- bcr_multinomial_supremum(X[, setdiff(seq_len(ncol(X)), target), drop = FALSE], y, ...)
  statistic <- 2 * (reduced$nll - full$nll)
  rounding_bound <- 2 * (max(0, full$gap) + max(0, reduced$gap) + 2e-10)
  ok <- full$certified && reduced$certified && statistic >= -rounding_bound
  df <- length(target) * (nlevels(as.factor(y)) - 1L)
  list(pvalue = if (ok) pchisq(max(0, statistic), df, lower.tail = FALSE) else NA_real_,
    statistic = statistic, rounding_bound = rounding_bound, df = df, full = full, reduced = reduced,
    objective_status = if (ok) "numerically_certified" else "uncertified",
    inference_status = "working_chisquare_reference_calibration_not_proved")
}
