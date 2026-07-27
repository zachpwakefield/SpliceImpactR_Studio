SIR_STATE_FIELDS <- c(
  "reference", "reference_info", "sample_frame", "raw_events", "di",
  "matched", "hits_sequences", "pairs", "sequence_results",
  "protein_features", "exon_features", "ppi", "domain_results",
  "final_results", "background", "integrated", "enrichment",
  "probe", "browser", "hit_compare", "psi_overview"
)

SIR_INVALIDATION <- list(
  reference = c(
    "protein_features", "exon_features", "ppi", "matched", "hits_sequences",
    "pairs", "sequence_results", "domain_results", "final_results",
    "background", "integrated", "enrichment", "browser"
  ),
  events = c(
    "di", "matched", "hits_sequences", "pairs", "sequence_results",
    "domain_results", "final_results", "background", "integrated",
    "enrichment", "probe", "hit_compare", "psi_overview"
  ),
  di = c(
    "matched", "hits_sequences", "pairs", "sequence_results",
    "domain_results", "final_results", "integrated", "enrichment"
  ),
  thresholds = c(
    "matched", "hits_sequences", "pairs", "sequence_results",
    "domain_results", "final_results", "integrated", "enrichment"
  ),
  pairing = c(
    "hits_sequences", "pairs", "sequence_results", "domain_results",
    "final_results", "integrated", "enrichment"
  ),
  features = c(
    "exon_features", "domain_results", "final_results", "background",
    "integrated", "enrichment", "browser"
  ),
  ppi = c("final_results", "integrated", "enrichment"),
  background = c("enrichment")
)

sir_invalidation_targets <- function(cause) {
  targets <- SIR_INVALIDATION[[cause]]
  if (is.null(targets)) {
    stop("Unknown invalidation cause: ", cause, call. = FALSE)
  }
  targets
}

sir_invalidate_state <- function(rv, cause) {
  targets <- sir_invalidation_targets(cause)
  changed <- character()
  for (field in targets) {
    if (!is.null(rv[[field]])) {
      changed <- c(changed, field)
    }
    rv[[field]] <- NULL
  }
  provenance <- rv$provenance %||% list()
  provenance[intersect(names(provenance), targets)] <- NULL
  rv$provenance <- provenance
  invisible(changed)
}

sir_record_stage <- function(rv, stage, signature, details = list()) {
  provenance <- rv$provenance %||% list()
  provenance[[stage]] <- c(
    list(
      signature = signature,
      completed_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
    ),
    details
  )
  rv$provenance <- provenance
  invisible(provenance[[stage]])
}

sir_stage_signature <- function(rv, stage) {
  (rv$provenance %||% list())[[stage]]$signature %||% NULL
}

sir_workspace_manifest <- function(rv) {
  table_count <- function(x) {
    if (is.null(x) || !is.data.frame(x)) NULL else nrow(x)
  }
  unique_count <- function(x, column) {
    if (is.null(x) || !is.data.frame(x) || !column %in% names(x)) {
      return(NULL)
    }
    data.table::uniqueN(data.table::as.data.table(x)[[column]], na.rm = TRUE)
  }

  reference_info <- rv$reference_info %||% NULL
  if (is.list(reference_info)) {
    reference_info$paths <- NULL
  }

  list(
    app = list(
      name = "SpliceImpactR Studio",
      version = SIR_APP_VERSION,
      created_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
      spliceimpactr_version = as.character(utils::packageVersion("SpliceImpactR")),
      spliceimpactr_expected_sha = SIR_PACKAGE_SHA
    ),
    source = rv$source_meta %||% NULL,
    reference = reference_info,
    conditions = rv$condition_meta %||% NULL,
    counts = list(
      samples = table_count(rv$sample_frame),
      raw_rows = table_count(rv$raw_events),
      raw_events = unique_count(rv$raw_events, "event_id"),
      di_rows = table_count(rv$di),
      di_events = unique_count(rv$di, "event_id"),
      matched_rows = table_count(rv$matched),
      matched_events = unique_count(rv$matched, "event_id"),
      sequence_pairs = table_count(rv$sequence_results),
      domain_pairs = table_count(rv$domain_results),
      final_pairs = table_count(rv$final_results),
      feature_rows = table_count(rv$protein_features),
      ppi_rows = table_count(rv$ppi),
      background_genes = table_count(rv$background)
    ),
    provenance = rv$provenance %||% list()
  )
}
