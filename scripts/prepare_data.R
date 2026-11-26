# Explicit import or preparation; input directories are never written.
args <- commandArgs(TRUE)
if (!length(args) %in% c(2L, 3L) || !args[1L] %in% c('import', 'raw'))
  stop('Usage: Rscript scripts/prepare_data.R raw|import INPUT [all|stock|bcr]')
dataset <- if (length(args) == 3L) args[3L] else 'all'
if (!dataset %in% c('all', 'stock', 'bcr')) stop('Dataset must be all, stock, or bcr')
selected <- if (dataset == 'all') c('stock', 'bcr') else dataset
script <- sub('^--file=', '', grep('^--file=', commandArgs(), value = TRUE)[1L])
root <- normalizePath(file.path(dirname(script), '..'))
input <- normalizePath(args[2L], mustWork = TRUE)
local <- file.path(root, 'data/processed')
source(file.path(root, 'R/run_helpers.R'))
source(file.path(root, 'R/comparators.R'))
source(file.path(root, 'R/data.R'))
source(file.path(root, 'R/stock_data.R'))
source(file.path(root, 'R/bcr_data.R'))
RNGkind('Mersenne-Twister', 'Inversion', 'Rejection')
objects <- list()
if (args[1L] == 'raw') {
  expected <- jsonlite::fromJSON(file.path(root, 'data/required-inputs.json'))
  prefixes <- c(stock = 'Stock', bcr = 'Breast')[selected]
  expected <- expected[dirname(expected$file) %in% prefixes, , drop = FALSE]
  missing <- expected$file[!file.exists(file.path(input, expected$file))]
  if (length(missing)) stop('Missing input files: ', paste(missing, collapse = ', '))
  for (i in seq_len(nrow(expected))) {
    actual <- digest::digest(file = file.path(input, expected$file[i]), algo = 'sha256', serialize = FALSE)
    if (actual != expected$sha256[i]) stop('Input snapshot hash differs: ', expected$file[i])
  }
  if ('stock' %in% selected) {
    stock <- prepare_stock(file.path(input, 'Stock'))
    # Replace the incomplete historical final DIA window using the full 461-date snapshot.
    chart <- jsonlite::fromJSON(file.path(input, 'Stock/DIA.json'), simplifyVector = FALSE)$chart$result[[1L]]
    prices <- as.numeric(unlist(chart$indicators$adjclose[[1L]]$adjclose))
    dates <- as.Date(as.POSIXct(as.numeric(unlist(chart$timestamp)), origin = '1970-01-01', tz = 'UTC'))
    stopifnot(length(prices) == 461L, length(dates) == 461L, !anyDuplicated(dates), all(is.finite(prices)), all(prices > 0))
    returns <- exp(diff(log(prices))) - 1
    groups <- (seq_along(returns) - 1L) %/% 5L + 1L
    stopifnot(all(tabulate(groups) == 5L))
    stock$y <- as.numeric(tapply(returns, groups, mean)) * 100
    stock$nulls <- lapply(stock$targets, function(target) with_seed(2023L, {
      outside <- setdiff(seq_len(ncol(stock$X)), target)
      fit <- glmnet::cv.glmnet(stock$X[, outside, drop = FALSE], stock$y)
      predicted <- predict(fit, newx = stock$X[, outside, drop = FALSE], s = fit$lambda.min)
      list(coef = as.numeric(coef(fit, s = fit$lambda.min)), sd = sd(stock$y - predicted), lambda = fit$lambda.min, outside = outside, seed = 2023L)
    }))
    stock$starting_point <- 'Original 103-stock matrix; corrected 92 complete five-day DIA windows.'
    objects[['stock.rds']] <- stock
  }
  if ('bcr' %in% selected)
    objects[['bcr_erplus60.rds']] <- prepare_bcr(file.path(input, 'Breast'), selection_seed = 2023L, null_seed = 2023L)
} else {
  filenames <- c(stock = 'stock.rds', bcr = 'bcr_erplus60.rds')[selected]
  for (name in filenames) objects[[name]] <- readRDS(file.path(input, name))
}
if ('stock' %in% selected) {
  stock <- objects[['stock.rds']]
  stopifnot(identical(dim(stock$X), c(92L, 103L)), length(stock$y) == 92L)
}
if ('bcr' %in% selected) {
  bcr <- objects[['bcr_erplus60.rds']]
  stopifnot(
    identical(dim(bcr$X), c(209L, 21L)), length(bcr$target) == 18L, length(bcr$condition_probes) == 3L,
    identical(bcr$panel, 'wang_erplus60')
  )
}
destinations <- file.path(local, names(objects))
if (any(file.exists(destinations)))
  stop('Refusing to overwrite: ', paste(destinations[file.exists(destinations)], collapse = ', '))
dir.create(local, recursive = TRUE, showWarnings = FALSE)
for (name in names(objects)) {
  saveRDS(objects[[name]], file.path(local, name), version = 3L)
}
manifest <- data.frame(name = names(objects), input_mode = args[1L], sha256 = vapply(names(objects), function(n)
  digest::digest(file = file.path(local, n), algo = 'sha256', serialize = FALSE), ''))
manifest_path <- file.path(local, 'IMPORT.csv')
if (file.exists(manifest_path)) {
  previous <- read.csv(manifest_path, stringsAsFactors = FALSE)
  manifest <- rbind(previous[!previous$name %in% manifest$name, , drop = FALSE], manifest)
}
write.csv(manifest[order(manifest$name), ], manifest_path, row.names = FALSE)
cat('Prepared: ', paste(names(objects), collapse = ', '), '\n', sep = '')
