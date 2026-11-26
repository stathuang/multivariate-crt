# Export the two named functions without sourcing unrelated legacy code.
args <- commandArgs(TRUE)
stopifnot(length(args) == 2L, !file.exists(args[2L]))
definitions <- list()
for (expr in parse(args[1L])) {
  if (is.call(expr) && as.character(expr[[1L]]) %in% c('<-', '=') &&
      is.symbol(expr[[2L]]) && as.character(expr[[2L]]) %in%
        c('dCRT', 'd_CRT', 'Creat_condition_gaussian')) {
    name <- as.character(expr[[2L]])
    stopifnot(is.call(expr[[3L]]), identical(expr[[3L]][[1L]], as.name('function')))
    definitions[[name]] <- paste(deparse(expr, width.cutoff=500L), collapse='\n')
  }
}
stopifnot(length(definitions) == 2L, 'Creat_condition_gaussian' %in% names(definitions))
dir.create(dirname(args[2L]), recursive=TRUE, showWarnings=FALSE)
writeLines(unlist(definitions, use.names=FALSE), args[2L])
