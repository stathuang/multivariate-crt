# Reconstructed historical simulations. The caller owns the replicate seed.
simulation_configs <- function(path = NULL) {
  if (is.null(path)) {
    candidates <- file.path(c(getOption("mscrt.study_root", "."),
                             "experiments/reproducible-study"), "config/simulations.json")
    path <- candidates[file.exists(candidates)][1L]
  }
  if (is.na(path) || !file.exists(path)) stop("Cannot locate config/simulations.json")
  cfg <- jsonlite::fromJSON(path, simplifyVector = FALSE)$scenarios
  vectors <- c("target", "theta", "statistics", "baselines", "graph_bandwidths",
               "t3_lag_coefficients", "adjustments", "fixed_beta", "graph_ids", "methods")
  cfg <- lapply(cfg, function(x) {
    for (nm in intersect(vectors, names(x))) x[[nm]] <- unlist(x[[nm]], use.names = FALSE)
    x
  })
  names(cfg) <- vapply(cfg, `[[`, "", "id")
  cfg
}

simulation_band_graph <- function(p, bandwidth) {
  d <- abs(outer(seq_len(p), seq_len(p), `-`))
  1L * (d > 0 & d <= bandwidth)
}

simulation_t3 <- function(n, p, coefficients = c(.2,.2,-.2,-.2,.2,.2),
                          df = 3, noise_scale = .7) {
  X <- matrix(0, n, p)
  X[,1L] <- stats::rt(n, df)
  if (p > 1L) for (j in 2:p) {
    lag <- seq_len(min(j - 1L, length(coefficients)))
    # Lag 1 is the immediately preceding column, including j=2,...,6.
    X[,j] <- drop(X[, j-lag, drop=FALSE] %*% coefficients[lag]) +
      noise_scale * stats::rt(n, df)
  }
  X
}

simulation_design <- function(config) {
  p <- config$p
  graph <- simulation_band_graph(p, config$bandwidth)
  precision <- diag(p) + config$precision_offdiag * graph
  covariance <- solve(precision)
  covariance <- stats::cov2cor(covariance)
  if (config$family == "t3") {
    X <- simulation_t3(config$n, p, config$t3_lag_coefficients,
                       config$t3_df, config$t3_noise_scale)
    permutation <- sample.int(p)
    X <- X[,permutation,drop=FALSE]
  } else {
    permutation <- sample.int(p)
    covariance <- covariance[permutation,permutation,drop=FALSE]
    # MASS::mvrnorm matches the historical Gaussian generator.
    X <- MASS::mvrnorm(config$n, rep(0,p), covariance)
  }
  graphs <- lapply(config$graph_bandwidths, function(k)
    simulation_band_graph(p,k)[permutation,permutation,drop=FALSE])
  names(graphs) <- paste0("K=", config$graph_bandwidths)
  beta <- numeric(p)
  if (config$coefficient_policy == "historical_fixed_seed_2025") {
    beta[seq_along(config$fixed_beta)] <- config$fixed_beta
  } else {
    active <- seq_len(min(p,20L))
    beta[active] <- stats::runif(length(active),1,2)/sqrt(length(active))
  }
  list(X=X, covariance=covariance, graphs=graphs, beta=beta,
       noise_uniform=stats::runif(config$n), permutation=permutation)
}

simulation_response <- function(X, setting, target, theta, beta, noise_uniform) {
  model <- substring(setting,3L)
  if (model %in% c("l_w","gl_w")) {
    beta[target] <- beta[target] * theta
    predictor <- drop(X %*% beta)
  } else if (model == "l_m") {
    basis <- function(z) cbind(z[,1]^2/2, 1/(1+z[,2]^2), cos(pi*z[,3]),
      sin(pi*z[,4]), z[,5]*z[,6], sin(z[,7]*z[,8]), sin(pi*z[,9])/(4+z[,10]^2))
    predictor <- drop(cbind(basis(X[,1:10]), basis(X[,11:20])) %*% c(rep(theta,7),rep(1,7)))
  } else if (model == "gl_m") {
    z <- X[,1:20,drop=FALSE]
    b <- ifelse(z >= stats::qnorm(1/3) & z <= stats::qnorm(2/3), -.5, .25)
    predictor <- drop(b %*% c(rep(theta,10),rep(1,10)))
  } else stop("Unknown simulation setting: ", setting)
  if (startsWith(model,"gl")) as.numeric(noise_uniform < stats::plogis(predictor))
  else predictor + stats::qnorm(noise_uniform)
}

simulation_capture <- function(expr) {
  warning_text <- character()
  start <- proc.time()[["elapsed"]]
  error_text <- ""; error_condition <- NULL
  value <- tryCatch(withCallingHandlers(force(expr), warning=function(w) {
    warning_text <<- c(warning_text, conditionMessage(w)); invokeRestart("muffleWarning")
  }), error=function(e) { error_text <<- conditionMessage(e); error_condition <<- e; NULL })
  list(value=value, condition=error_condition, status=if (nzchar(error_text)) "failed" else "ok",
       error=error_text, warnings=paste(unique(warning_text),collapse=" | "),
       elapsed_seconds=unname(proc.time()[["elapsed"]]-start))
}

simulation_pvalue <- function(x) {
  if (length(x)!=1L || !is.numeric(x) || !is.finite(x) || x<0 || x>1)
    stop("Method returned an invalid p-value")
  as.numeric(x)
}

simulation_row <- function(config, replicate, theta, graph, method, result) {
  row <- data.frame(scenario=config$id, replicate=as.integer(replicate), theta=theta,
    graph=graph, method=method, pvalue=if(result$status=="ok") result$value else NA_real_,
    status=result$status, error=result$error, warnings=result$warnings,
    elapsed_seconds=result$elapsed_seconds, stringsAsFactors=FALSE)
  if (exists("add_result_diagnostics", mode="function")) row <- add_result_diagnostics(row, result$diagnostics, result$condition)
  row
}

simulation_replicate <- function(config, replicate) {
  if (length(replicate)!=1L || !is.finite(replicate) || replicate<1 || replicate!=as.integer(replicate))
    stop("replicate must be a positive integer")
  design_result <- simulation_capture(simulation_design(config))
  if (design_result$status=="failed") {
    cells <- expected_simulation_rows(config)
    return(do.call(rbind,lapply(seq_len(nrow(cells)),function(i)
      if (config$family=="group_selection")
        simulation_group_metric(config,replicate,cells$theta[i],cells$method[i],design_result,integer())
      else simulation_row(config,replicate,cells$theta[i],cells$graph[i],cells$method[i],design_result))))
  }
  design <- design_result$value
  if (config$family == "group_selection") {
    result <- simulation_group_replicate(config,replicate,design)
    result$warnings <- vapply(result$warnings,function(w) paste(c(design_result$warnings,w)[nzchar(c(design_result$warnings,w))],collapse=" | "),"")
    return(result)
  }
  X <- design$X
  cache <- simulation_prepare_de_lasso(config,replicate,X)
  rows <- list()
  for (graph_name in names(design$graphs)) {
    sampling <- simulation_capture(mscrt_sample(X,design$graphs[[graph_name]],config$target,
                                               M=config$M,L=config$L,reduced=TRUE))
    for (theta in config$theta) {
      response <- simulation_capture(simulation_response(X,config$setting,config$target,theta,design$beta,design$noise_uniform))
      y <- response$value
      for (statistic in config$statistics) {
        if (response$status=="failed") result <- response else if (sampling$status == "failed") result <- sampling else {
          result <- simulation_capture(simulation_pvalue(gcrt_test(X,y,config$target,
            sampling$value,statistic,config$response_family,config$p_rule)))
          result$warnings <- paste(c(sampling$warnings,result$warnings)[nzchar(c(sampling$warnings,result$warnings))],collapse=" | ")
          result$elapsed_seconds <- result$elapsed_seconds + sampling$elapsed_seconds
        }
        result$warnings <- paste(c(response$warnings,result$warnings)[nzchar(c(response$warnings,result$warnings))],collapse=" | ")
        rows[[length(rows)+1L]] <- simulation_row(config,replicate,theta,graph_name,statistic,result)
      }
      for (method in config$baselines) {
        if (response$status=="failed") result <- response else if(method=="De-Lasso" && cache$enabled) {
          result <- cache$result
          if(result$status=="ok") {
            response_seed <- simulation_de_lasso_seed(config,replicate,"response",match(theta,config$theta))
            result <- simulation_capture(simulation_with_seed(response_seed,
              simulation_pvalue(baseline_test_cached_candidate(X,y,config$target,method,
                config$response_family,cache$result$value,
                suppress_unused_group_testing=config$de_lasso_suppress_unused_group_testing))))
            result$warnings <- paste(c(cache$result$warnings,result$warnings)[nzchar(c(cache$result$warnings,result$warnings))],collapse=" | ")
            result$elapsed_seconds <- result$elapsed_seconds + cache$result$elapsed_seconds
          }
        } else {
          result <- simulation_capture(if (exists("baseline_test_detail", mode="function"))
            baseline_test_detail(X,y,config$target,method,config$response_family) else
            list(pvalue=baseline_test(X,y,config$target,method,config$response_family)))
          if (result$status == "ok") {
            result$diagnostics <- result$value$diagnostics
            result$value <- simulation_pvalue(result$value$pvalue)
          }
        }
        result$warnings <- paste(c(response$warnings,result$warnings)[nzchar(c(response$warnings,result$warnings))],collapse=" | ")
        rows[[length(rows)+1L]] <- simulation_row(config,replicate,theta,graph_name,method,result)
      }
    }
  }
  result <- if (exists("bind_diagnostic_rows", mode="function")) bind_diagnostic_rows(rows) else do.call(rbind,rows)
  if (cache$enabled) {
    is_de_lasso <- result$method=="De-Lasso"
    result$de_lasso_policy <- config$de_lasso_policy
    result$de_lasso_seed_policy <- config$de_lasso_seed_policy
    result$de_lasso_nodewise_seed <- ifelse(is_de_lasso,cache$nodewise_seed,NA_integer_)
    result$de_lasso_response_seed <- vapply(seq_len(nrow(result)),function(i)
      if(is_de_lasso[i]) simulation_de_lasso_seed(config,replicate,"response",match(result$theta[i],config$theta)) else NA_integer_,0L)
  }
  result$warnings <- vapply(result$warnings,function(w) paste(c(design_result$warnings,w)[nzchar(c(design_result$warnings,w))],collapse=" | "),"")
  result
}

simulation_group_select <- function(pvalues, method, q=.1) {
  if (!length(pvalues) || any(!is.finite(pvalues) | pvalues<0 | pvalues>1))
    stop("Group selection requires a valid p-value for every group")
  if (method=="e-BH") {
    g <- length(pvalues)
    cutoffs <- cumsum(q/g/(1+seq_len(g)))
    index <- rowSums(outer(pvalues,cutoffs,`>`))
    evalues <- g/q/(index+1)
    which(stats::p.adjust(1/evalues,"BH") <= q)
  } else {
    if (!method %in% c("BH","BY")) stop("Unknown adjustment")
    which(stats::p.adjust(pvalues,method) <= q)
  }
}

# Julia initialization never downloads a runtime or installs Julia packages.
simulation_knockoff_setup <- function() {
  if (!requireNamespace("knockoffsr",quietly=TRUE)) stop("Required package knockoffsr is unavailable")
  executable <- Sys.getenv("MSCRT_JULIA",unset=Sys.which("julia"))
  if (!nzchar(executable) || !file.exists(executable))
    stop("Group knockoff blocked: no Julia executable on PATH; set MSCRT_JULIA to an existing runtime (no automatic installation)")
  JuliaCall::julia_setup(JULIA_HOME=dirname(normalizePath(executable)), installJulia=FALSE, install=FALSE)
  knockoffsr::knockoff_setup(julia_dir=dirname(normalizePath(executable)),pkg_check=FALSE,
                           installJulia=FALSE,install=FALSE)
}

simulation_group_metric <- function(config,replicate,theta,method,result,true_groups) {
  selected <- result$value
  if (result$status=="ok" && (any(!is.finite(selected)) || any(selected!=as.integer(selected)) ||
      any(selected<1 | selected>config$p/config$group_size) || anyDuplicated(selected))) {
    result$status <- "failed"; result$error <- "Invalid selected group indices"; selected <- NULL
  }
  result$value <- NA_real_
  row <- simulation_row(config,replicate,theta,config$graph_ids[1L],method,result)
  row$record_type <- "group_metric"
  row$selected <- if (result$status=="ok") paste(selected,collapse=";") else NA_character_
  row$TP <- if(result$status=="ok") sum(selected %in% true_groups) else NA_integer_
  row$FP <- if(result$status=="ok") sum(!selected %in% true_groups) else NA_integer_
  row$FDP <- if(result$status=="ok") row$FP/max(1,length(selected)) else NA_real_
  row$power <- if(result$status=="ok") row$TP/length(true_groups) else NA_real_
  row
}

simulation_group_replicate <- function(config,replicate,design) {
  X <- design$X
  groups <- split(seq_len(config$p),ceiling(seq_len(config$p)/config$group_size))
  true_groups <- which(vapply(groups,function(g) any(design$beta[g]!=0),logical(1)))
  ko <- simulation_capture({
    # Draw before setup; isolate any cold-start R draws from the replicate stream.
    julia_seed <- sample.int(.Machine$integer.max,1L)
    api <- simulation_with_seed(julia_seed,simulation_knockoff_setup())
    JuliaCall::julia_command(sprintf("using Random; Random.seed!(%d)",julia_seed))
    object <- api$modelX_gaussian_group_knockoffs(X,"maxent",as.integer(ceiling(seq_len(config$p)/config$group_size)),
                                                 rep(0,config$p),design$covariance,m=1,verbose=FALSE)
    list(api=api,object=object)
  })
  copies <- lapply(groups,function(g) simulation_capture(mscrt_sample(X,design$graphs[[1]],g,
                                                  M=config$M,L=config$L,reduced=TRUE)))
  rows <- list()
  for(theta in config$theta) {
    response <- simulation_capture(simulation_response(X,config$setting,which(design$beta!=0),theta,design$beta,design$noise_uniform))
    y <- response$value
    for(statistic in config$statistics) {
      pv <- if(response$status=="failed") response else simulation_capture(vapply(seq_along(groups),function(j) {
        if(copies[[j]]$status=="failed") stop("Group ",j," sampler: ",copies[[j]]$error)
        simulation_pvalue(gcrt_test(X,y,groups[[j]],copies[[j]]$value,statistic,config$response_family,config$p_rule))
      },numeric(1)))
      pv$warnings <- paste(unique(c(response$warnings,pv$warnings,vapply(copies,`[[`,"","warnings"))),collapse=" | ")
      pv$elapsed_seconds <- pv$elapsed_seconds + sum(vapply(copies,`[[`,0,"elapsed_seconds"))
      for (adjustment in config$adjustments) {
        result <- pv
        if(result$status=="ok") result$value <- simulation_group_select(pv$value,adjustment,config$q_fdr)
        rows[[length(rows)+1L]] <- simulation_group_metric(config,replicate,theta,adjustment,result,true_groups)
      }
    }
    result <- if(response$status=="failed") response else ko
    if(response$status=="ok" && ko$status=="ok") {
      result <- simulation_capture({
        fit <- ko$value$api$fit_lasso(y,ko$value$object)
        # fit$fdr_target[j] is the target FDR for the group indices in fit$selected[[j]].
        # Match the configured level (0.1 here), not its position in the default grid.
        fdr <- as.numeric(fit$fdr_target)
        # Allow floating-point rounding; exactly one target must match.
        index <- which(abs(fdr-config$q_fdr)<1e-12)
        stopifnot(length(index)==1L)
        # Flatten the R/Julia result to integer group indices; empty means no selection.
        # The surrounding simulation_capture records any extraction error as a failure.
        as.integer(unlist(fit$selected[[index]],use.names=FALSE))
      })
      result$elapsed_seconds <- result$elapsed_seconds + ko$elapsed_seconds
    }
    result$warnings <- paste(c(response$warnings,result$warnings)[nzchar(c(response$warnings,result$warnings))],collapse=" | ")
    rows[[length(rows)+1L]] <- simulation_group_metric(config,replicate,theta,"Second-Order-Group",result,true_groups)
  }
  do.call(rbind,rows)
}

expected_simulation_rows <- function(config) {
  expand.grid(theta=config$theta, graph=config$graph_ids, method=config$methods,
              KEEP.OUT.ATTRS=FALSE, stringsAsFactors=FALSE)
}

# Explicit campaign-v2 scoped streams. These restore the caller's RNG state;
# they apply only to the opt-in Gaussian De-Lasso cache, never to generators.
simulation_with_seed <- function(seed, expr) {
  had_seed <- exists(".Random.seed",envir=.GlobalEnv,inherits=FALSE)
  if(had_seed) previous <- get(".Random.seed",envir=.GlobalEnv,inherits=FALSE)
  on.exit(if(had_seed) assign(".Random.seed",previous,envir=.GlobalEnv)
          else if(exists(".Random.seed",envir=.GlobalEnv,inherits=FALSE)) rm(".Random.seed",envir=.GlobalEnv),add=TRUE)
  set.seed(seed)
  force(expr)
}

simulation_de_lasso_seed <- function(config, replicate, stage, theta_index=0L) {
  offset <- switch(stage,nodewise=1357911,response=2468022+104729*theta_index,
                   stop("Unknown De-Lasso RNG stage"))
  as.integer((as.double(config$seed_base)+replicate+offset) %% 2147483646 + 1)
}

simulation_prepare_de_lasso <- function(config, replicate, X) {
  policy <- config$de_lasso_policy
  if(is.null(policy) || policy=="uncached-ordered-v1") return(list(enabled=FALSE))
  if(policy!="one-nodewise-fit-per-design-v1") stop("Unknown De-Lasso cache policy")
  if(config$response_family!="gaussian" || !"De-Lasso" %in% config$baselines)
    stop("De-Lasso cache policy requires Gaussian De-Lasso in the method list")
  if(!identical(config$de_lasso_seed_policy,"separate-nodewise-response-seeds-v1"))
    stop("The Gaussian cache requires its explicit separate-stream policy")
  if(!is.logical(config$de_lasso_suppress_unused_group_testing) ||
      length(config$de_lasso_suppress_unused_group_testing)!=1L || is.na(config$de_lasso_suppress_unused_group_testing))
    stop("Specify the unused-group-testing policy explicitly")
  seed <- simulation_de_lasso_seed(config,replicate,"nodewise")
  result <- simulation_capture(simulation_with_seed(seed,hdi_gaussian_cache(X)))
  list(enabled=TRUE,nodewise_seed=seed,result=result)
}
