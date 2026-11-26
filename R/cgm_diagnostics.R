# Version-pinned comparator diagnostics; no package namespace mutation.
# It preserves SIHR's initial fit and random stream. Only the zero-retained
# boundary gains a typed error before an undefined projection is attempted.
cgm_lf_diagnostic <- function(...) {
  ns <- asNamespace('SIHR')
  expected <- c(LF = 'd92d868be0bb7bedfd6d6e0407119c0386bf8253d3af4d83df5a18fdc30c8550',
    relevant.funs = '68ec2047353184c560537c870f9a1e8648d8b5a5b5ad168ea828eae426fdc277')
  if (as.character(utils::packageVersion('SIHR')) != '2.0.1') stop('CGM diagnostic adapter requires SIHR 2.0.1')
  for (nm in names(expected)) {
    f <- get(nm, ns)
    if (!identical(digest::digest(list(formals(f), body(f)), algo = 'sha256'), unname(expected[nm])))
      stop('CGM diagnostic adapter runtime body mismatch: ', nm)

  }
  state <- new.env(parent = emptyenv())
  state$lambda <- NULL
  private <- new.env(parent = ns)
  original_relevant <- get('relevant.funs', ns)
  private$relevant.funs <- function(...) {
    answer <- original_relevant(...)
    # A child environment overrides only the CV function lookup. The original
    # training closure, arguments, selected coefficient rule and RNG stay intact.
    train_env <- new.env(parent = environment(answer$train.fun))
    original_cv <- get('cv.glmnet', ns)
    train_env$cv.glmnet <- function(...) {
      fit <- original_cv(...)
      state$lambda <- c(lambda.min = fit$lambda.min, lambda.1se = fit$lambda.1se)
      fit

    }
    environment(answer$train.fun) <- train_env
    answer

  }
  private$empty_filter_error <- function(n, n_retained, threshold, lambda, beta_supplied) {
    selected <- if (beta_supplied) NA_real_ else
    if (is.null(lambda) || identical(lambda, 'CV.min')) state$lambda['lambda.min'] else
    if (identical(lambda, 'CV')) state$lambda['lambda.1se'] else lambda
    structure(list(message = 'CGM probability filter retained zero observations; projection is undefined', call = NULL,
      n = as.integer(n), n_retained = as.integer(n_retained), prob_filter = threshold,
      selected_lambda = unname(selected), initial_beta_supplied = beta_supplied),
    class = c('CGM_NO_RETAINED_OBSERVATIONS', 'error', 'condition'))

  }
  fit <- get('LF', ns)
  statements <- as.list(body(fit))
  stopifnot(identical(statements[[33]], quote({
    n.filter = nrow(X.filter)

  })[[2L]]))
  # Capture supplied-beta status before LF replaces NULL with its initial fit.
  statements <- append(statements, list(quote(.beta_supplied <- !is.null(beta.init))), after = 1L)
  statements <- append(statements,
    list(quote(if (n.filter == 0L)
      stop(empty_filter_error(n, n.filter, prob.filter, lambda, .beta_supplied)))),
    after = 34L)
  body(fit) <- as.call(statements)
  environment(fit) <- private
  fit(...)

}
