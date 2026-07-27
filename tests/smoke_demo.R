source("src/00_config.R")
source("src/05_resources.R")
source("src/10_utils.R")
source("src/20_io.R")
source("src/30_state.R")
source("src/40_pipeline.R")
source("src/50_ui.R")
source("src/60_server.R")

shiny::testServer(sir_server, {
  session$setInputs(annotation_mode = "library", input_mode = "project")
  session$flushReact()
  library_controls <- paste(output$reference_source_controls, collapse = "")
  full_references <- setdiff(
    names(sir_reference_catalog()),
    "demo_human_v45"
  )
  stopifnot(
    if (length(full_references)) {
      grepl("reference_library_id", library_controls, fixed = TRUE)
    } else {
      grepl("No server-library reference found", library_controls, fixed = TRUE)
    },
    grepl(
      "project_archive",
      paste(output$input_mode_controls, collapse = ""),
      fixed = TRUE
    )
  )

  session$setInputs(annotation_mode = "cached")
  session$flushReact()
  stopifnot(grepl(
    "annotation_species",
    paste(output$reference_source_controls, collapse = ""),
    fixed = TRUE
  ))

  session$setInputs(annotation_mode = "upload")
  session$flushReact()
  upload_controls <- paste(output$reference_source_controls, collapse = "")
  stopifnot(
    grepl("annotation_gtf_upload", upload_controls, fixed = TRUE),
    grepl("annotation_transcripts_upload", upload_controls, fixed = TRUE),
    grepl("annotation_translations_upload", upload_controls, fixed = TRUE),
    grepl("gencode.v45.annotation.gtf.gz", upload_controls, fixed = TRUE)
  )

  hosted_reference_modes <- unname(sir_reference_mode_choices(hosted = TRUE))
  if ("bundled" %in% hosted_reference_modes) {
    session$setInputs(annotation_mode = "bundled", load_reference = 1)
    session$flushReact()
    stopifnot(
      identical(rv$reference_info$kind, "preprocessed"),
      identical(rv$reference_info$gencode_release, 45L),
      nrow(rv$reference$annotations) > 800000L,
      nrow(rv$reference$sequences) > 100000L
    )
  }

  session$setInputs(annotation_mode = "test", input_mode = "demo")
  session$flushReact()
  session$setInputs(load_demo_workspace = 1)
  session$flushReact()

  stopifnot(
    !is.null(rv$reference),
    nrow(rv$raw_events) == 5863L,
    nrow(rv$protein_features) == 1449L
  )

  stopifnot(
    identical(resolved_conditions()$reference, "control"),
    identical(resolved_conditions()$case, "case")
  )
  session$setInputs(
    min_reads = 10,
    min_prop = 0.5,
    cooks_cutoff = "Inf",
    padj_threshold = 0.05,
    dpsi_threshold = 0.1,
    pair_mode = "multi",
    show_protein_domains = FALSE,
    preview_rows = 25,
    run_all_pipeline = 1
  )
  session$flushReact()
  stopifnot(
    nrow(rv$di) == 533L,
    identical(rv$condition_meta$reference, "control"),
    identical(rv$condition_meta$case, "case"),
    data.table::uniqueN(rv$matched$event_id) == 40L,
    nrow(rv$sequence_results) == 36L,
    nrow(rv$domain_results) == 36L,
    sir_nonempty_table(rv$ppi),
    "n_ppi" %in% names(rv$final_results),
    !is.null(rv$integrated),
    identical(rv$final_kind, "sequence + domains + PPI")
  )

  stopifnot(identical(overview_event_types(), SIR_SUPPORTED_EVENTS))

  a3ss_gene <- rv$raw_events[event_type == "A3SS", gene_id][[1]]
  session$setInputs(
    probe_gene = a3ss_gene,
    probe_event_type = "A3SS"
  )
  session$flushReact()
  a3ss_event <- probe_event_ids()[[1]]
  session$setInputs(probe_event = a3ss_event)
  session$flushReact()
  session$setInputs(run_probe = 1)
  session$flushReact()
  stopifnot(!is.null(rv$probe), nrow(rv$probe$data) > 0L)

  browser_index <- annotation_index()
  browser_gene <- browser_index[
    ,
    .(transcripts = data.table::uniqueN(transcript_id)),
    by = gene_id
  ][transcripts >= 2L, gene_id][[1]]
  browser_transcripts <- sort(unique(
    browser_index[gene_id == browser_gene, transcript_id]
  ))
  session$setInputs(
    browser_gene = browser_gene,
    transcript_a = browser_transcripts[[1]],
    transcript_b = browser_transcripts[[2]],
    browser_view = "protein",
    browser_feature_databases = unique(rv$protein_features$database),
    browser_combine_domains = TRUE,
    run_browser = 1
  )
  session$flushReact()
  stopifnot(
    !is.null(rv$browser),
    nrow(rv$browser$sequence) > 0L,
    inherits(rv$browser$plot, "ggplot")
  )

  session$setInputs(
    reference_condition = "control",
    case_condition = "case",
    overview_event_type = "A3SS"
  )
  session$flushReact()
  session$setInputs(run_overview = 1)
  session$flushReact()
  stopifnot(!is.null(rv$hit_compare), !is.null(rv$psi_overview))

  session$setInputs(overview_event_type = "ALE")
  session$flushReact()
  session$setInputs(run_overview = 2)
  session$flushReact()
  stopifnot(!is.null(rv$hit_compare), !is.null(rv$psi_overview))

  first_mapping_signature <- sir_stage_signature(rv, "matched")
  session$setInputs(dpsi_threshold = 0.5)
  session$flushReact()

  stopifnot(
    is.null(current_matched()),
    is.null(current_final()),
    data.table::uniqueN(significant_di()$event_id) == 17L
  )

  session$setInputs(run_consequences = 1)
  session$flushReact()

  stopifnot(
    !is.null(current_matched()),
    !identical(first_mapping_signature, sir_stage_signature(rv, "matched")),
    data.table::uniqueN(rv$matched$event_id) == 16L,
    nrow(rv$final_results) == 23L
  )

  imported <- data.table::copy(rv$di)
  replace_di_evidence(imported, list(type = "smoke_import", label = "smoke test"))
  stopifnot(is.null(rv$raw_events), identical(rv$di_origin, "imported"))

  session$setInputs(run_all_pipeline = 2)
  session$flushReact()
  stopifnot(!is.null(current_matched()), !is.null(current_final()))
})

message("SpliceImpactR v4 demo smoke test passed.")
