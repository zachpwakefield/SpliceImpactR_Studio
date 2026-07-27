app_source <- tryCatch(sys.frame(1)$ofile, error = function(error) NULL)
app_dir <- if (!is.null(app_source) && nzchar(app_source)) {
  dirname(normalizePath(app_source, winslash = "/", mustWork = TRUE))
} else {
  normalizePath(getwd(), winslash = "/", mustWork = TRUE)
}
options(spliceimpactr.app_dir = app_dir)

source(file.path(app_dir, "src", "00_config.R"), local = FALSE)
source(file.path(app_dir, "src", "05_resources.R"), local = FALSE)
source(file.path(app_dir, "src", "10_utils.R"), local = FALSE)
source(file.path(app_dir, "src", "20_io.R"), local = FALSE)
source(file.path(app_dir, "src", "30_state.R"), local = FALSE)
source(file.path(app_dir, "src", "40_pipeline.R"), local = FALSE)
source(file.path(app_dir, "src", "50_ui.R"), local = FALSE)
source(file.path(app_dir, "src", "60_server.R"), local = FALSE)

shiny::shinyApp(
  ui = sir_ui(),
  server = sir_server
)
