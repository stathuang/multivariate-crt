# Gaussian graphical co-sufficient hub-and-spoke sampling.
# All randomness belongs to the caller's RNG stream. No seed is set here.

.ms_matrix <- function(X, name = "X") {
  if (!is.matrix(X) || !is.numeric(X) || any(!is.finite(X)) ||
      nrow(X) < 1L || ncol(X) < 1L) stop(name, " must be a finite numeric matrix")
  invisible(X)
}

.ms_target <- function(target, p) {
  if (!is.numeric(target) || !length(target) || anyNA(target) ||
      any(target != as.integer(target)) || any(target < 1L | target > p) ||
      anyDuplicated(target)) stop("target must contain unique column indices")
  as.integer(target)
}

.ms_count <- function(x, name) {
  if (length(x) != 1L || !is.numeric(x) || !is.finite(x) ||
      x < 0 || x != floor(x) || x > .Machine$integer.max)
    stop(name, " must be a nonnegative integer")
  as.integer(x)
}

# Preserve the historical QR direction and random sign on ordinary full-rank
# inputs. Rank-deficient designs use the independent columns of the QR basis.
.ms_orthogonal <- function(A) {
  q <- qr(A)
  if (q$rank >= nrow(A)) return(NULL)
  if (q$rank < ncol(A)) A <- A[, q$pivot[seq_len(q$rank)], drop = FALSE]
  qr.Q(qr(cbind(A, stats::rnorm(nrow(A)))))[, ncol(A) + 1L] *
    sample(c(-1, 1), 1L)
}

.ms_rotate <- function(X, i, graph) {
  neighbors <- setdiff(which(graph[i, ] != 0), i)
  ones <- matrix(1, nrow(X), 1L)
  if (!length(neighbors)) {
    mu <- mean(X[, i])
    radius <- sqrt(sum((X[, i] - mu)^2))
    direction <- .ms_orthogonal(ones)
    if (is.null(direction)) return(X[, i])
    return(mu + radius * direction)
  }
  # lm is deliberate: its QR arithmetic matches the original implementation.
  fit <- stats::lm(X[, i] ~ X[, neighbors, drop = FALSE])
  direction <- .ms_orthogonal(cbind(ones, X[, neighbors, drop = FALSE]))
  if (is.null(direction)) return(X[, i])
  fit$fitted.values + sqrt(sum(fit$residuals^2)) * direction
}

mscrt_sample <- function(X, graph, target = seq_len(ncol(X)), M = 100L,
                         L = 3L, reduced = TRUE) {
  .ms_matrix(X)
  target <- .ms_target(target, ncol(X))
  M <- .ms_count(M, "M"); L <- .ms_count(L, "L")
  if (!is.matrix(graph) || !(is.numeric(graph) || is.logical(graph)) ||
      !identical(dim(graph), c(ncol(X), ncol(X))) || any(!is.finite(graph)) ||
      any((graph != 0) != t(graph != 0))) stop("graph must be a finite symmetric adjacency matrix")
  if (length(reduced) != 1L || !is.logical(reduced) || is.na(reduced))
    stop("reduced must be TRUE or FALSE")
  sweep <- function(Z, order) {
    for (i in order) Z[, i] <- .ms_rotate(Z, i, graph)
    Z
  }
  hub <- X
  for (l in seq_len(L)) hub <- sweep(hub, target)
  copies <- vector("list", M)
  for (m in seq_len(M)) {
    candidate <- hub
    for (l in seq_len(L)) candidate <- sweep(candidate, rev(target))
    copies[[m]] <- if (reduced) candidate[, target, drop = FALSE] else candidate
  }
  list(X_hub = hub, X_copies = copies, reduced = reduced, target = target)
}

# Upper-tail rank. Exact floating-point equality defines ties; no tolerance is
# applied because changing a tie group changes the randomization distribution.
rank_pvalue <- function(observed, randomized, rule = "discrete") {
  if (length(observed) != 1L || !is.numeric(observed) || !is.finite(observed) ||
      !is.numeric(randomized) || !is.null(dim(randomized)) ||
      any(!is.finite(randomized))) stop("statistics must be finite numeric scalars/vector")
  rule <- match.arg(rule, c("discrete", "continuous", "conservative"))
  greater <- sum(randomized > observed)
  equal <- 1L + sum(randomized == observed)
  offset <- switch(rule, discrete = sample.int(equal, 1L),
                   continuous = stats::runif(1L, 0, equal), conservative = equal)
  (greater + offset) / (length(randomized) + 1L)
}
