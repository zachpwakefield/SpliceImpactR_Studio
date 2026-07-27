suppressPackageStartupMessages({
  library(shiny)
  library(bslib)
  library(DT)
  library(data.table)
  library(ggplot2)
  library(SpliceImpactR)
})

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0L) y else x
}

SIR_APP_VERSION <- "4.0.0"
SIR_PACKAGE_SHA <- "5436813f56606af40cc33e917f2b8fb75bb90c12"
SIR_SUPPORTED_EVENTS <- c("ALE", "AFE", "MXE", "SE", "A3SS", "A5SS", "RI")
SIR_EVENT_LABELS <- c(
  ALE = "ALE · Alternative last exon",
  AFE = "AFE · Alternative first exon",
  MXE = "MXE · Mutually exclusive exons",
  SE = "SE · Skipped exon",
  A3SS = "A3SS · Alternative 3′ splice site",
  A5SS = "A5SS · Alternative 5′ splice site",
  RI = "RI · Retained intron"
)
SIR_FEATURE_SOURCES <- c(
  "interpro", "pfam", "gene3d", "signalp", "tmhmm",
  "ncoils", "seg", "mobidblite", "elm"
)
SIR_ENRICHMENT_SOURCES <- c(
  "GO:BP",
  "GO:MF",
  "GO:CC",
  "MSigDB:H",
  "MSigDB:C2:CP:REACTOME",
  "MSigDB:C5:GO:BP"
)

sir_env_number <- function(name, default, min = 1, max = Inf) {
  raw <- trimws(Sys.getenv(name, unset = ""))
  value <- suppressWarnings(as.numeric(raw))
  if (!nzchar(raw) || !is.finite(value) || value < min || value > max) {
    return(default)
  }
  value
}

SIR_MAX_UPLOAD_MB <- sir_env_number(
  "SPLICEIMPACTR_MAX_UPLOAD_MB",
  default = 200,
  min = 10,
  max = 2048
)
SIR_MAX_ARCHIVE_EXPANDED_MB <- sir_env_number(
  "SPLICEIMPACTR_MAX_ARCHIVE_EXPANDED_MB",
  default = 1500,
  min = 50,
  max = 10240
)
SIR_MAX_ARCHIVE_FILES <- as.integer(sir_env_number(
  "SPLICEIMPACTR_MAX_ARCHIVE_FILES",
  default = 5000,
  min = 10,
  max = 100000
))

options(shiny.maxRequestSize = SIR_MAX_UPLOAD_MB * 1024^2)

SIR_APP_DIR <- getOption(
  "spliceimpactr.app_dir",
  normalizePath(".", winslash = "/", mustWork = TRUE)
)

.sir_process_cache <- new.env(parent = emptyenv())

sir_is_hosted_runtime <- function() {
  override <- tolower(trimws(Sys.getenv(
    "SPLICEIMPACTR_HOSTED",
    unset = ""
  )))
  if (override %in% c("1", "true", "yes", "on")) {
    return(TRUE)
  }
  if (override %in% c("0", "false", "no", "off")) {
    return(FALSE)
  }

  hosted_markers <- c(
    "SHINY_PORT",
    "RSCONNECT_USER",
    "CONNECT_SERVER",
    "SHINY_SERVER_VERSION",
    "SHINYPROXY_USERNAME",
    "DYNO",
    "K_SERVICE",
    "WEBSITE_INSTANCE_ID"
  )
  hosted_environment <- any(nzchar(Sys.getenv(
    hosted_markers,
    unset = ""
  )))
  hosted_path <- grepl(
    "^/(srv/(connect/apps|shiny-server)|opt/shiny-server)(/|$)",
    normalizePath(SIR_APP_DIR, winslash = "/", mustWork = FALSE)
  )

  hosted_environment || hosted_path
}

sir_load_demo_features <- function() {
  cache_key <- "demo_features"
  if (exists(cache_key, envir = .sir_process_cache, inherits = FALSE)) {
    return(get(cache_key, envir = .sir_process_cache, inherits = FALSE))
  }
  ex <- SpliceImpactR::load_example_data(c("protein_feature_total", "exon_features"))
  value <- list(
    protein_feature_total = data.table::as.data.table(ex$protein_feature_total),
    exon_features = data.table::as.data.table(ex$exon_features)
  )
  assign(cache_key, value, envir = .sir_process_cache)
  value
}

sir_load_default_ppi <- function() {
  cache_key <- "default_ppi"
  if (exists(cache_key, envir = .sir_process_cache, inherits = FALSE)) {
    return(get(cache_key, envir = .sir_process_cache, inherits = FALSE))
  }
  value <- data.table::as.data.table(SpliceImpactR::get_ppi_interactions())
  assign(cache_key, value, envir = .sir_process_cache)
  value
}

sir_reference_index <- function(reference_id, reference) {
  cache_key <- paste0("reference_index__", reference_id)
  if (exists(cache_key, envir = .sir_process_cache, inherits = FALSE)) {
    return(get(cache_key, envir = .sir_process_cache, inherits = FALSE))
  }
  annotations <- data.table::as.data.table(reference$annotations)
  columns <- intersect(
    c("gene_id", "gene_name", "transcript_id", "protein_id", "type", "chr", "strand"),
    names(annotations)
  )
  index <- unique(annotations[
    (is.na(type) | type == "transcript") &
      !is.na(transcript_id) &
      nzchar(transcript_id),
    ..columns
  ])
  assign(cache_key, index, envir = .sir_process_cache)
  index
}
