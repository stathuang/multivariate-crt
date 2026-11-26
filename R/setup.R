# Entry scripts call this from the repository root.
options(stringsAsFactors = FALSE)
root <- normalizePath('.', mustWork = TRUE)
args <- commandArgs(TRUE)
quick <- '--quick' %in% args
argument <- function(flag, default = NULL) {
  i <- match(flag, args)
  if (is.na(i)) return(default)
  if (i == length(args)) stop('Missing value for ', flag)
  args[i + 1L]
}
load_modules <- function(names, envir = .GlobalEnv) {
  for (name in names) sys.source(file.path(root, 'R', paste0(name, '.R')), envir)
}
load_modules('run_helpers')
output_directory <- function(study, part) {
  path <- file.path(root, 'results', if (quick) 'quick' else 'generated', study, part)
  if (dir.exists(path) && length(list.files(path))) stop(
    'Output already exists: ', path,
    '. Move or remove this generated directory before rerunning.'
  )
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  writeLines(capture.output(sessionInfo()), file.path(path, 'sessionInfo.txt'))
  path
}
save_settings <- function(config, out) jsonlite::write_json(list(mode = if (quick) 'quick' else 'full', settings = config), file.path(out, 'settings.json'),
  pretty = TRUE, auto_unbox = TRUE, digits = 17
)
RNGkind('Mersenne-Twister', 'Inversion', 'Rejection')
