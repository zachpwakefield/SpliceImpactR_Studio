bootstrap_source <- tryCatch(sys.frame(1)$ofile, error = function(error) NULL)
app_dir <- if (!is.null(bootstrap_source) && nzchar(bootstrap_source)) {
  dirname(normalizePath(bootstrap_source, winslash = "/", mustWork = TRUE))
} else {
  normalizePath(getwd(), winslash = "/", mustWork = TRUE)
}

if (!isTRUE(getOption("spliceimpactr.v4.loaded"))) {
  options(
    spliceimpactr.app_dir = app_dir,
    spliceimpactr.v4.loaded = TRUE
  )
  source(file.path(app_dir, "src", "00_config.R"), local = FALSE)
  source(file.path(app_dir, "src", "05_resources.R"), local = FALSE)
  source(file.path(app_dir, "src", "10_utils.R"), local = FALSE)
  source(file.path(app_dir, "src", "20_io.R"), local = FALSE)
  source(file.path(app_dir, "src", "30_state.R"), local = FALSE)
  source(file.path(app_dir, "src", "40_pipeline.R"), local = FALSE)
  source(file.path(app_dir, "src", "50_ui.R"), local = FALSE)
  source(file.path(app_dir, "src", "60_server.R"), local = FALSE)
}
