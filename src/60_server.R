sir_server <- function(input, output, session) {
  session_root <- sir_ensure_directory(SIR_SESSION_ROOT, "Session workspace root")
  session_dir <- tempfile("spliceimpactr-session-", tmpdir = session_root)
  dir.create(session_dir, recursive = TRUE, showWarnings = FALSE)
  session$onSessionEnded(function() {
    unlink(session_dir, recursive = TRUE, force = TRUE)
  })

  rv <- reactiveValues(
    reference = NULL,
    reference_info = NULL,
    sample_frame = NULL,
    raw_events = NULL,
    di = NULL,
    matched = NULL,
    hits_sequences = NULL,
    pairs = NULL,
    sequence_results = NULL,
    protein_features = NULL,
    exon_features = NULL,
    ppi = NULL,
    domain_results = NULL,
    final_results = NULL,
    background = NULL,
    integrated = NULL,
    enrichment = NULL,
    probe = NULL,
    browser = NULL,
    hit_compare = NULL,
    psi_overview = NULL,
    source_meta = NULL,
    condition_meta = NULL,
    di_origin = NULL,
    reference_token = NULL,
    events_token = NULL,
    feature_token = NULL,
    ppi_token = NULL,
    background_token = NULL,
    final_kind = NULL,
    provenance = list(),
    logs = data.table::data.table(
      time = character(),
      level = character(),
      message = character()
    )
  )

  append_log <- function(level, message, notify = TRUE) {
    level <- toupper(level)
    message <- as.character(message)
    entry <- data.table::data.table(
      time = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
      level = level,
      message = message
    )
    rv$logs <- data.table::rbindlist(list(rv$logs, entry), use.names = TRUE, fill = TRUE)
    if (isTRUE(notify)) {
      type <- switch(
        level,
        ERROR = "error",
        WARNING = "warning",
        SUCCESS = "message",
        INFO = "default",
        "default"
      )
      showNotification(message, type = type, duration = if (level == "ERROR") 9 else 6)
    }
    invisible(entry)
  }

  run_safely <- function(label, expression) {
    tryCatch(
      force(expression),
      error = function(error) {
        append_log(
          "error",
          paste0(label, " failed: ", sir_error_text(error))
        )
        NULL
      }
    )
  }

  new_token <- function(kind, ...) {
    sir_signature(kind, format(Sys.time(), "%Y-%m-%dT%H:%M:%OS6%z"), runif(1), ...)
  }

  preview_rows <- reactive({
    value <- suppressWarnings(as.integer(input$preview_rows %||% 25L))
    max(1L, min(value, 250L))
  })

  set_button_state <- function(id, enabled) {
    session$sendCustomMessage(
      "sir-button-state",
      list(id = id, disabled = !isTRUE(enabled))
    )
  }

  reset_workspace <- function(clear_log = TRUE) {
    for (field in SIR_STATE_FIELDS) {
      rv[[field]] <- NULL
    }
    rv$source_meta <- NULL
    rv$condition_meta <- NULL
    rv$di_origin <- NULL
    rv$reference_token <- NULL
    rv$events_token <- NULL
    rv$feature_token <- NULL
    rv$ppi_token <- NULL
    rv$background_token <- NULL
    rv$final_kind <- NULL
    rv$provenance <- list()
    if (isTRUE(clear_log)) {
      rv$logs <- data.table::data.table(
        time = character(),
        level = character(),
        message = character()
      )
    }
    invisible(TRUE)
  }

  uploaded_reference_request <- function(prepare_files = TRUE) {
    filter_tsl <- input$annotation_tsl %||% c("1", "2", "3")
    request <- sir_reference_request(
      mode = "cached",
      species = "human",
      release = 45L,
      filter_tsl = filter_tsl
    )
    request$id <- "uploaded_human_v45"
    request$label <- "Uploaded human GENCODE v45 reference"
    request$short_label <- "Uploaded human GENCODE v45"
    request$kind <- "uploaded_raw"
    request$assembly <- "GRCh38.p14"
    request$note <- paste0(
      "Official GTF and FASTA assets uploaded through the browser. ",
      "Files and processed objects remain isolated to this session."
    )
    if (!isTRUE(prepare_files)) {
      return(request)
    }

    uploads <- list(
      gtf = input$annotation_gtf_upload,
      transcripts = input$annotation_transcripts_upload,
      translations = input$annotation_translations_upload
    )
    missing <- names(uploads)[vapply(
      uploads,
      function(fileinfo) is.null(fileinfo) || !nrow(fileinfo),
      logical(1)
    )]
    if (length(missing)) {
      stop(
        "Upload all three GENCODE files: annotation GTF, transcript FASTA, ",
        "and protein FASTA.",
        call. = FALSE
      )
    }
    total_size <- sum(vapply(
      uploads,
      function(fileinfo) as.numeric(fileinfo$size[[1]] %||% 0),
      numeric(1)
    ))
    if (
      is.finite(total_size) &&
      total_size > SIR_MAX_UPLOAD_MB * 1024^2
    ) {
      stop(
        "The three reference files total more than the configured ",
        sir_fmt_num(SIR_MAX_UPLOAD_MB, 0),
        " MB upload limit.",
        call. = FALSE
      )
    }

    upload_root <- tempfile("reference-upload-", tmpdir = session_dir)
    source_dir <- file.path(upload_root, "source")
    paths <- c(
      gtf = sir_copy_reference_upload(
        uploads$gtf,
        source_dir,
        "gtf"
      ),
      transcripts = sir_copy_reference_upload(
        uploads$transcripts,
        source_dir,
        "transcripts"
      ),
      translations = sir_copy_reference_upload(
        uploads$translations,
        source_dir,
        "translations"
      )
    )
    upload_signature <- sir_signature(lapply(
      uploads,
      function(fileinfo) {
        list(
          name = basename(fileinfo$name[[1]]),
          size = fileinfo$size[[1]],
          type = fileinfo$type[[1]],
          datapath = fileinfo$datapath[[1]]
        )
      }
    ))
    request$id <- paste0("uploaded_human_v45_", upload_signature)
    request$paths <- paths
    request$workspace <- file.path(upload_root, "processed")
    request
  }

  reference_request_from_input <- function(prepare_upload = FALSE) {
    mode <- input$annotation_mode %||% sir_default_reference_mode()
    if (identical(mode, "bundled")) {
      return(sir_bundled_human_reference())
    }
    if (identical(mode, "upload")) {
      return(uploaded_reference_request(prepare_files = prepare_upload))
    }
    sir_reference_request(
      mode = mode,
      library_id = input$reference_library_id %||% NULL,
      species = input$annotation_species %||% "human",
      release = input$annotation_release %||% 45L,
      filter_tsl = input$annotation_tsl %||% c("1", "2", "3")
    )
  }

  set_reference <- function(request) {
    if (is.character(request) && length(request) == 1L) {
      request <- sir_reference_info(request)
    }
    reference <- sir_load_reference(request)
    info <- request
    sir_invalidate_state(rv, "reference")
    rv$reference <- reference
    rv$reference_info <- sir_reference_public_info(info)
    rv$reference_token <- new_token(
      "reference",
      info$id,
      nrow(reference$annotations),
      nrow(reference$sequences)
    )
    sir_record_stage(
      rv,
      "reference",
      rv$reference_token,
      list(
        id = info$id,
        kind = info$kind,
        species = info$species,
        gencode_release = info$gencode_release,
        ensembl_release = info$ensembl_release,
        assembly = info$assembly,
        annotation_rows = nrow(reference$annotations),
        sequence_rows = nrow(reference$sequences)
      )
    )
    updateRadioButtons(
      session,
      "feature_mode",
      selected = if (identical(info$kind, "demo")) "demo" else "query"
    )
    if (length(info$ensembl_release) && is.finite(info$ensembl_release)) {
      updateNumericInput(
        session,
        "biomart_release",
        value = as.integer(info$ensembl_release)
      )
    }
    append_log("success", paste0(info$label, " loaded."))
    invisible(reference)
  }

  replace_raw_evidence <- function(raw_events, sample_frame = NULL, source_meta = list()) {
    raw <- data.table::as.data.table(raw_events)
    if (!nrow(raw)) {
      stop("No raw event rows were loaded.", call. = FALSE)
    }
    counts <- sir_validate_raw_conditions(raw)
    sir_invalidate_state(rv, "events")
    rv$sample_frame <- if (is.null(sample_frame)) NULL else as.data.frame(sample_frame)
    rv$raw_events <- raw
    rv$di <- NULL
    rv$di_origin <- NULL
    rv$source_meta <- source_meta
    rv$condition_meta <- NULL
    rv$events_token <- new_token(
      "events",
      source_meta,
      nrow(raw),
      data.table::uniqueN(raw$event_id),
      counts
    )
    sir_record_stage(
      rv,
      "raw_events",
      rv$events_token,
      list(
        rows = nrow(raw),
        events = data.table::uniqueN(raw$event_id),
        conditions = as.list(stats::setNames(counts$samples, counts$condition))
      )
    )
    append_log(
      "success",
      paste0(
        "Loaded ", sir_fmt_int(nrow(raw)), " raw rows across ",
        sir_fmt_int(data.table::uniqueN(raw$event_id)), " events and ",
        sir_fmt_int(nrow(counts)), " conditions."
      )
    )
    invisible(raw)
  }

  replace_di_evidence <- function(di, source_meta = list()) {
    result <- sir_normalize_di(di)
    sir_invalidate_state(rv, "events")
    rv$sample_frame <- NULL
    rv$raw_events <- NULL
    rv$di <- result
    rv$di_origin <- "imported"
    rv$source_meta <- source_meta
    rv$condition_meta <- NULL
    rv$events_token <- NULL
    signature <- new_token(
      "imported-di",
      source_meta,
      nrow(result),
      data.table::uniqueN(result$event_id)
    )
    sir_record_stage(
      rv,
      "di",
      signature,
      list(
        origin = "imported",
        rows = nrow(result),
        events = data.table::uniqueN(result$event_id)
      )
    )
    append_log(
      "success",
      paste0(
        "Imported ", sir_fmt_int(nrow(result)), " differential-inclusion rows across ",
        sir_fmt_int(data.table::uniqueN(result$event_id)), " events."
      )
    )
    invisible(result)
  }

  set_features <- function(features, exon_features, origin) {
    features <- data.table::as.data.table(features)
    exon_features <- data.table::as.data.table(exon_features)
    if (!nrow(features)) {
      stop("The feature source returned no rows.", call. = FALSE)
    }
    if (!nrow(exon_features)) {
      stop(
        "No protein features overlapped coding exons in the loaded reference.",
        call. = FALSE
      )
    }
    sir_invalidate_state(rv, "features")
    rv$protein_features <- features
    rv$exon_features <- exon_features
    rv$feature_token <- new_token("features", origin, nrow(features), nrow(exon_features))
    coverage <- sir_feature_coverage(features, rv$reference)
    sir_record_stage(
      rv,
      "protein_features",
      rv$feature_token,
      list(
        origin = origin,
        rows = nrow(features),
        exon_feature_rows = nrow(exon_features),
        covered_transcripts = coverage$feature_transcripts,
        reference_transcripts = coverage$reference_transcripts,
        coverage = coverage$proportion
      )
    )
    append_log(
      "success",
      paste0(
        "Loaded ", sir_fmt_int(nrow(features)), " feature rows covering ",
        sir_fmt_int(coverage$feature_transcripts), " reference transcripts."
      )
    )
    if (
      !identical(rv$reference_info$kind, "demo") &&
      is.finite(coverage$proportion) &&
      coverage$proportion < 0.2
    ) {
      append_log(
        "warning",
        paste0(
          "Feature coverage is only ",
          sir_fmt_num(100 * coverage$proportion, 1, "%"),
          " of reference transcripts; domain results may be sparse."
        )
      )
    }
    invisible(features)
  }

  set_ppi <- function(ppi, origin, source_signature) {
    ppi <- sir_validate_ppi_table(ppi)
    sir_invalidate_state(rv, "ppi")
    rv$ppi <- ppi
    rv$ppi_token <- new_token(
      "ppi",
      origin,
      source_signature,
      nrow(ppi)
    )
    sir_record_stage(
      rv,
      "ppi",
      rv$ppi_token,
      list(
        origin = origin,
        source_signature = source_signature,
        rows = nrow(ppi)
      )
    )
    append_log(
      "success",
      paste0(
        "Loaded ", sir_fmt_int(nrow(ppi)), " PPI rows from ", origin, "."
      )
    )
    invisible(ppi)
  }

  load_selected_ppi <- function(force = FALSE) {
    mode <- input$ppi_mode %||% "default"
    if (!mode %in% c("default", "upload")) {
      stop("Unknown PPI source.", call. = FALSE)
    }

    if (identical(mode, "upload")) {
      upload <- input$ppi_upload
      if (is.null(upload) || !nrow(upload)) {
        stop(
          "The uploaded PPI source is selected, but no PPI table was provided.",
          call. = FALSE
        )
      }
      origin <- paste0("uploaded network: ", basename(upload$name[[1]]))
      source_signature <- sir_signature(
        "ppi-upload-v1",
        basename(upload$name[[1]]),
        upload$size[[1]],
        upload$type[[1]],
        upload$datapath[[1]]
      )
    } else {
      origin <- "bundled SpliceImpactR network"
      source_signature <- sir_signature("ppi-default-v1")
    }

    current_signature <- (
      (rv$provenance %||% list())$ppi$source_signature %||% ""
    )
    if (
      !isTRUE(force) &&
      sir_nonempty_table(rv$ppi) &&
      identical(current_signature, source_signature)
    ) {
      return(rv$ppi)
    }

    ppi <- if (identical(mode, "upload")) {
      sir_read_uploaded_table(input$ppi_upload, session_dir)
    } else {
      sir_load_default_ppi()
    }
    set_ppi(ppi, origin, source_signature)
  }

  raw_condition_counts <- reactive({
    if (!sir_nonempty_table(rv$raw_events)) {
      return(data.table::data.table())
    }
    sir_validate_raw_conditions(rv$raw_events)
  })

  resolved_conditions <- reactive({
    if (!sir_nonempty_table(rv$raw_events)) {
      return(NULL)
    }
    sir_resolve_condition_selection(
      rv$raw_events,
      input$reference_condition,
      input$case_condition
    )
  })

  expected_di_signature <- reactive({
    if (identical(rv$di_origin, "imported") && sir_nonempty_table(rv$di)) {
      return(sir_stage_signature(rv, "di"))
    }
    if (!sir_nonempty_table(rv$raw_events) || is.null(rv$events_token)) {
      return(NULL)
    }
    conditions <- resolved_conditions()
    if (is.null(conditions)) {
      return(NULL)
    }
    sir_signature(
      "computed-di-v1",
      rv$events_token,
      conditions$reference,
      conditions$case,
      sir_di_parameters(input)
    )
  })

  current_di <- reactive({
    if (!sir_nonempty_table(rv$di)) {
      return(NULL)
    }
    expected <- expected_di_signature()
    actual <- sir_stage_signature(rv, "di")
    if (is.null(expected) || !identical(expected, actual)) {
      return(NULL)
    }
    rv$di
  })

  significant_di <- reactive({
    di <- current_di()
    if (is.null(di)) {
      return(data.table::data.table())
    }
    sir_significant_di(
      di,
      padj_threshold = input$padj_threshold,
      dpsi_threshold = input$dpsi_threshold
    )
  })

  expected_matched_signature <- reactive({
    if (is.null(rv$reference_token) || is.null(current_di())) {
      return(NULL)
    }
    sir_signature(
      "matched-v1",
      rv$reference_token,
      sir_stage_signature(rv, "di"),
      sir_threshold_parameters(input)
    )
  })

  current_matched <- reactive({
    expected <- expected_matched_signature()
    if (
      !sir_nonempty_table(rv$matched) ||
      is.null(expected) ||
      !identical(sir_stage_signature(rv, "matched"), expected)
    ) {
      return(NULL)
    }
    rv$matched
  })

  expected_sequence_signature <- reactive({
    if (is.null(current_matched())) {
      return(NULL)
    }
    sir_signature(
      "sequence-v1",
      sir_stage_signature(rv, "matched"),
      input$pair_mode %||% "multi"
    )
  })

  current_sequence <- reactive({
    expected <- expected_sequence_signature()
    if (
      !sir_nonempty_table(rv$sequence_results) ||
      is.null(expected) ||
      !identical(sir_stage_signature(rv, "sequence_results"), expected)
    ) {
      return(NULL)
    }
    rv$sequence_results
  })

  expected_domain_signature <- reactive({
    if (is.null(current_sequence()) || is.null(rv$feature_token)) {
      return(NULL)
    }
    sir_signature(
      "domains-v1",
      sir_stage_signature(rv, "sequence_results"),
      rv$feature_token,
      isTRUE(input$show_protein_domains)
    )
  })

  current_domains <- reactive({
    expected <- expected_domain_signature()
    if (
      !sir_nonempty_table(rv$domain_results) ||
      is.null(expected) ||
      !identical(sir_stage_signature(rv, "domain_results"), expected)
    ) {
      return(NULL)
    }
    rv$domain_results
  })

  expected_final_signature <- reactive({
    if (!is.null(current_domains())) {
      if (!is.null(rv$ppi_token) && sir_nonempty_table(rv$ppi)) {
        return(sir_signature(
          "final-v1",
          sir_stage_signature(rv, "domain_results"),
          rv$ppi_token
        ))
      }
      return(sir_signature(
        "final-v1",
        sir_stage_signature(rv, "domain_results"),
        "no-ppi"
      ))
    }
    if (!is.null(current_sequence())) {
      return(sir_signature(
        "final-v1",
        sir_stage_signature(rv, "sequence_results"),
        "sequence-only"
      ))
    }
    NULL
  })

  current_final <- reactive({
    expected <- expected_final_signature()
    if (
      !sir_nonempty_table(rv$final_results) ||
      is.null(expected) ||
      !identical(sir_stage_signature(rv, "final_results"), expected)
    ) {
      return(NULL)
    }
    rv$final_results
  })

  parse_background_transcripts <- reactive({
    text <- input$background_transcripts %||% ""
    values <- trimws(unlist(strsplit(text, "[,;\\n\\r\\t]+", perl = TRUE)))
    unique(values[nzchar(values)])
  })

  expected_background_signature <- reactive({
    if (is.null(rv$reference_token) || is.null(rv$feature_token)) {
      return(NULL)
    }
    source <- input$background_source %||% "annotated"
    source_token <- switch(
      source,
      annotated = "all-annotated",
      hit_index = {
        conditions <- resolved_conditions()
        if (is.null(conditions)) {
          NULL
        } else {
          list(
            events = rv$events_token,
            reference_condition = conditions$reference,
            case_condition = conditions$case
          )
        }
      },
      `user-given` = parse_background_transcripts(),
      NULL
    )
    if (is.null(source_token)) {
      return(NULL)
    }
    sir_signature(
      "background-v1",
      source,
      source_token,
      rv$reference_token,
      rv$feature_token
    )
  })

  current_background <- reactive({
    expected <- expected_background_signature()
    if (
      !sir_nonempty_table(rv$background) ||
      is.null(expected) ||
      !identical(sir_stage_signature(rv, "background"), expected)
    ) {
      return(NULL)
    }
    rv$background
  })

  current_integrated <- reactive({
    final <- current_final()
    if (is.null(final) || is.null(rv$integrated)) {
      return(NULL)
    }
    expected <- sir_signature(
      "integrated-v1",
      sir_stage_signature(rv, "final_results"),
      sir_stage_signature(rv, "di")
    )
    if (!identical(sir_stage_signature(rv, "integrated"), expected)) {
      return(NULL)
    }
    rv$integrated
  })

  compute_di <- function(force = FALSE) {
    if (identical(rv$di_origin, "imported") && sir_nonempty_table(rv$di)) {
      return(rv$di)
    }
    if (!sir_nonempty_table(rv$raw_events)) {
      stop(
        "Load raw events or import a differential-inclusion table first.",
        call. = FALSE
      )
    }
    if (!isTRUE(force) && !is.null(current_di())) {
      return(current_di())
    }

    conditions <- resolved_conditions()
    prepared <- sir_prepare_conditions(
      rv$raw_events,
      rv$sample_frame,
      conditions$reference,
      conditions$case
    )
    parameters <- sir_di_parameters(input)
    result <- SpliceImpactR::get_differential_inclusion(
      DT = data.table::copy(prepared$raw_events),
      min_total_reads = parameters$min_total_reads,
      minimum_proportion_containing_event = parameters$minimum_proportion,
      terminal_fill = parameters$terminal_fill,
      cooks_cutoff = parameters$cooks_cutoff,
      parallel_glm = FALSE,
      BPPARAM = BiocParallel::SerialParam(),
      verbose = FALSE
    )
    result <- sir_normalize_di(result)

    sir_invalidate_state(rv, "di")
    rv$di <- result
    rv$di_origin <- "computed"
    rv$condition_meta <- prepared$metadata
    signature <- expected_di_signature()
    sir_record_stage(
      rv,
      "di",
      signature,
      list(
        origin = "computed",
        parameters = parameters,
        conditions = prepared$metadata,
        rows = nrow(result),
        events = data.table::uniqueN(result$event_id)
      )
    )
    append_log(
      "success",
      paste0(
        "Differential inclusion completed for ",
        sir_fmt_int(nrow(result)), " event-form rows (",
        prepared$metadata$case, " versus ", prepared$metadata$reference, ")."
      )
    )
    result
  }

  compute_mapping <- function(force = FALSE) {
    if (is.null(rv$reference)) {
      stop("Load a reference before matching events.", call. = FALSE)
    }
    compute_di(force = FALSE)
    if (!isTRUE(force) && !is.null(current_matched())) {
      return(current_matched())
    }
    significant <- significant_di()
    if (!nrow(significant)) {
      stop(
        "No events pass the current adjusted p-value and ΔPSI cutoffs.",
        call. = FALSE
      )
    }
    signature <- expected_matched_signature()
    sir_invalidate_state(rv, "thresholds")
    result <- SpliceImpactR::get_matched_events_chunked(
      events = data.table::copy(significant),
      annotations = rv$reference$annotations,
      chunk_size = 2000L
    )
    result <- data.table::as.data.table(result)
    if (!nrow(result)) {
      stop(
        "Significant events did not map to the loaded reference. Check genome assembly and reference release.",
        call. = FALSE
      )
    }
    rv$matched <- result
    sir_record_stage(
      rv,
      "matched",
      signature,
      list(
        thresholds = sir_threshold_parameters(input),
        rows = nrow(result),
        events = data.table::uniqueN(result$event_id),
        transcripts = data.table::uniqueN(result$transcript_id)
      )
    )
    append_log(
      "success",
      paste0(
        "Mapped ", sir_fmt_int(data.table::uniqueN(result$event_id)),
        " significant events to ",
        sir_fmt_int(data.table::uniqueN(result$transcript_id)),
        " transcripts."
      )
    )
    result
  }

  compute_consequences <- function(force = FALSE) {
    matched <- compute_mapping(force = FALSE)
    sequence_signature <- expected_sequence_signature()
    sequence <- current_sequence()

    if (isTRUE(force) || is.null(sequence)) {
      sir_invalidate_state(rv, "pairing")
      hits_sequences <- SpliceImpactR::attach_sequences(
        data.table::copy(matched),
        rv$reference$sequences
      )
      pairs <- SpliceImpactR::get_pairs(
        data.table::copy(data.table::as.data.table(hits_sequences)),
        source = input$pair_mode
      )
      pairs <- data.table::as.data.table(pairs)
      if (!nrow(pairs)) {
        stop(
          "No compatible transcript pairs were produced. Try the multi pairing mode or inspect the mapped rows.",
          call. = FALSE
        )
      }
      sequence <- SpliceImpactR::compare_sequence_frame(
        data.table::copy(pairs),
        rv$reference$annotations
      )
      sequence <- data.table::as.data.table(sequence)
      if (!nrow(sequence)) {
        stop("Sequence comparison returned no transcript pairs.", call. = FALSE)
      }
      rv$hits_sequences <- data.table::as.data.table(hits_sequences)
      rv$pairs <- pairs
      rv$sequence_results <- sequence
      sir_record_stage(
        rv,
        "sequence_results",
        sequence_signature,
        list(
          pairing = input$pair_mode,
          rows = nrow(sequence),
          events = data.table::uniqueN(sequence$event_id)
        )
      )
      append_log(
        "success",
        paste0("Sequence comparison completed for ", sir_fmt_int(nrow(sequence)), " transcript pairs.")
      )
    }

    domains <- NULL
    if (sir_nonempty_table(rv$exon_features) && !is.null(rv$feature_token)) {
      domain_signature <- expected_domain_signature()
      domains <- current_domains()
      if (isTRUE(force) || is.null(domains)) {
        rv$domain_results <- NULL
        rv$final_results <- NULL
        rv$integrated <- NULL
        rv$enrichment <- NULL
        domains <- SpliceImpactR::get_domains(
          hits = data.table::copy(sequence),
          exon_features = rv$exon_features,
          show_protein_domains = isTRUE(input$show_protein_domains)
        )
        domains <- data.table::as.data.table(domains)
        rv$domain_results <- domains
        sir_record_stage(
          rv,
          "domain_results",
          domain_signature,
          list(
            rows = nrow(domains),
            rows_with_changes = if ("diff_n" %in% names(domains)) {
              sum(domains$diff_n > 0, na.rm = TRUE)
            } else {
              0L
            }
          )
        )
        append_log(
          "success",
          paste0(
            "Domain comparison completed; ",
            sir_fmt_int(sum(domains$diff_n > 0, na.rm = TRUE)),
            " pairs have a predicted domain change."
          )
        )
      }
    }

    final_signature <- expected_final_signature()
    final <- current_final()
    if (isTRUE(force) || is.null(final)) {
      if (!is.null(domains) && sir_nonempty_table(rv$ppi) && !is.null(rv$ppi_token)) {
        hit_genes <- unique(na.omit(as.character(domains$gene_id)))
        ppi_for_hits <- data.table::as.data.table(rv$ppi)[
          geneA %in% hit_genes | geneB %in% hit_genes
        ]
        final <- SpliceImpactR::get_ppi_switches(
          hits_domain = data.table::copy(domains),
          ppi = ppi_for_hits,
          protein_feature_total = rv$protein_features
        )
        final <- data.table::as.data.table(final)
        rv$final_kind <- "sequence + domains + PPI"
      } else if (!is.null(domains)) {
        final <- domains
        rv$final_kind <- "sequence + domains"
      } else {
        final <- sequence
        rv$final_kind <- "sequence only"
      }
      rv$final_results <- final
      rv$integrated <- NULL
      rv$enrichment <- NULL
      sir_record_stage(
        rv,
        "final_results",
        final_signature,
        list(kind = rv$final_kind, rows = nrow(final))
      )
    }

    required_integrated <- c(
      "summary_classification", "case_only_n", "control_only_n",
      "n_ppi", "n_case_ppi", "n_control_ppi"
    )
    if (all(required_integrated %in% names(final))) {
      integrated_signature <- sir_signature(
        "integrated-v1",
        sir_stage_signature(rv, "final_results"),
        sir_stage_signature(rv, "di")
      )
      if (
        isTRUE(force) ||
        is.null(rv$integrated) ||
        !identical(sir_stage_signature(rv, "integrated"), integrated_signature)
      ) {
        rv$integrated <- SpliceImpactR::integrated_event_summary(final, rv$di)
        sir_record_stage(
          rv,
          "integrated",
          integrated_signature,
          list(kind = rv$final_kind)
        )
      }
    }

    if (identical(rv$final_kind, "sequence only")) {
      append_log(
        "warning",
        "Sequence results are current. Load protein features to add domain consequences."
      )
    } else if (identical(rv$final_kind, "sequence + domains")) {
      append_log(
        "warning",
        "Sequence and domain results are current. Load a PPI network to add interaction changes."
      )
    } else {
      append_log(
        "success",
        paste0(
          "All consequence layers are current for ",
          sir_fmt_int(nrow(final)), " transcript pairs."
        )
      )
    }
    final
  }

  observe({
    raw_ready <- sir_nonempty_table(rv$raw_events)
    di_ready <- !is.null(current_di())
    evidence_ready <- raw_ready || sir_nonempty_table(rv$di)
    reference_ready <- !is.null(rv$reference)
    features_ready <- sir_nonempty_table(rv$protein_features)
    background_ready <- !is.null(current_background())
    overview_ready <- raw_ready && !is.null(rv$sample_frame)
    background_input_ready <- switch(
      input$background_source %||% "annotated",
      annotated = TRUE,
      hit_index = raw_ready && !is.null(rv$sample_frame),
      `user-given` = length(parse_background_transcripts()) > 0L,
      FALSE
    )

    set_button_state("run_di", raw_ready)
    set_button_state("run_analysis", reference_ready && evidence_ready)
    set_button_state(
      "run_all_pipeline",
      reference_ready && evidence_ready && features_ready
    )
    set_button_state("map_events", reference_ready && evidence_ready)
    set_button_state("run_consequences", reference_ready && evidence_ready)
    set_button_state("load_features", reference_ready)
    set_button_state(
      "build_background",
      reference_ready && features_ready && background_input_ready
    )
    set_button_state("run_probe", raw_ready)
    set_button_state("run_browser", reference_ready && features_ready)
    set_button_state("run_overview", overview_ready)
    set_button_state("run_enrichment", (di_ready || raw_ready) && background_ready)
    set_button_state("download_raw", raw_ready)
    set_button_state("download_di", di_ready)
    set_button_state("download_sig", di_ready)
    set_button_state("download_matched", !is.null(current_matched()))
    set_button_state("download_sequence", !is.null(current_sequence()))
    set_button_state("download_domains", !is.null(current_domains()))
    set_button_state(
      "download_ppi_results",
      !is.null(current_final()) && "n_ppi" %in% names(current_final())
    )
    set_button_state("download_final", !is.null(current_final()))
    set_button_state("download_enrichment", !is.null(current_enrichment()))
  })

  observeEvent(input$new_workspace, {
    reset_workspace(clear_log = TRUE)
    append_log("success", "Started a new empty workspace.")
    updateRadioButtons(
      session,
      "annotation_mode",
      selected = sir_default_reference_mode()
    )
    updateRadioButtons(session, "input_mode", selected = "demo")
  })

  observeEvent(input$clear_analysis, {
    imported_di <- identical(rv$di_origin, "imported")
    if (!imported_di) {
      rv$di <- NULL
      rv$di_origin <- NULL
      rv$condition_meta <- NULL
      provenance <- rv$provenance
      provenance$di <- NULL
      rv$provenance <- provenance
    }
    sir_invalidate_state(rv, "di")
    append_log(
      "success",
      if (imported_di) {
        "Cleared derived analysis while retaining the imported DI evidence."
      } else {
        "Cleared differential and downstream analysis while retaining loaded inputs."
      }
    )
  })

  observeEvent(input$load_reference, {
    withProgress(message = "Loading reference", value = 0, {
      run_safely("Reference load", {
        request <- reference_request_from_input(prepare_upload = TRUE)
        detail <- switch(
          request$kind,
          download = "Downloading into the managed cache",
          cached = "Reading the managed processed cache",
          raw = "Processing trusted local GTF and FASTA assets",
          uploaded_raw = "Processing the uploaded GTF and FASTA assets",
          "Reading or reusing the server reference"
        )
        incProgress(0.25, detail = detail)
        set_reference(request)
        incProgress(1, detail = "Reference ready")
      })
    })
  })

  observeEvent(input$load_demo_workspace, {
    withProgress(message = "Loading guided demo", value = 0, {
      run_safely("Demo workspace load", {
        reset_workspace(clear_log = TRUE)
        updateRadioButtons(session, "annotation_mode", selected = "test")
        updateRadioButtons(session, "input_mode", selected = "demo")
        incProgress(0.15, detail = "Reference")
        set_reference("demo_human_v45")
        incProgress(0.35, detail = "Sample-level splice events")
        sample_frame <- SpliceImpactR::load_example_data("sample_frame")$sample_frame
        raw <- SpliceImpactR::get_rmats_hit(
          sample_frame = sample_frame,
          event_types = SIR_SUPPORTED_EVENTS,
          use = "JCEC",
          keep_annotated_first_last = TRUE
        )
        replace_raw_evidence(
          raw,
          sample_frame,
          list(
            type = "demo",
            label = "Packaged SpliceImpactR demo samples",
            event_types = SIR_SUPPORTED_EVENTS,
            count_mode = "JCEC"
          )
        )
        incProgress(0.72, detail = "Protein features")
        features <- sir_load_demo_features()
        set_features(
          features$protein_feature_total,
          features$exon_features,
          origin = "packaged demo subset"
        )
        incProgress(1, detail = "Ready to analyze")
        append_log(
          "success",
          "Guided demo is ready. Review condition meaning, then run the current analysis."
        )
      })
    })
  })

  observeEvent(input$load_input, {
    withProgress(message = "Loading evidence", value = 0, {
      run_safely("Evidence load", {
        mode <- input$input_mode %||% "demo"
        event_types <- input$event_types %||% SIR_SUPPORTED_EVENTS
        if (!length(event_types) && mode %in% c("demo", "project")) {
          stop("Choose at least one event type.", call. = FALSE)
        }

        if (identical(mode, "demo")) {
          incProgress(0.2, detail = "Packaged sample manifest")
          sample_frame <- SpliceImpactR::load_example_data("sample_frame")$sample_frame
          raw <- SpliceImpactR::get_rmats_hit(
            sample_frame = sample_frame,
            event_types = event_types,
            use = "JCEC",
            keep_annotated_first_last = isTRUE(input$keep_first_last)
          )
          replace_raw_evidence(
            raw,
            sample_frame,
            list(
              type = "demo",
              label = "Packaged SpliceImpactR demo samples",
              event_types = event_types,
              count_mode = "JCEC"
            )
          )
        } else if (identical(mode, "project")) {
          incProgress(0.15, detail = "Validating project archive")
          project <- sir_extract_project_archive(
            input$project_archive,
            session_dir,
            event_types = event_types,
            count_mode = input$rmats_use
          )
          incProgress(0.45, detail = "Reading rMATS outputs")
          raw <- SpliceImpactR::get_rmats_hit(
            sample_frame = project$manifest,
            event_types = event_types,
            use = input$rmats_use,
            keep_annotated_first_last = isTRUE(input$keep_first_last)
          )
          replace_raw_evidence(
            raw,
            project$manifest,
            list(
              type = "project_zip",
              label = project$archive_name,
              event_types = event_types,
              count_mode = input$rmats_use
            )
          )
        } else if (identical(mode, "raw")) {
          incProgress(0.25, detail = "Reading normalized raw table")
          uploaded <- sir_read_uploaded_table(input$raw_upload, session_dir)
          raw <- SpliceImpactR::get_user_data(data.table::copy(uploaded))
          replace_raw_evidence(
            raw,
            sample_frame = NULL,
            source_meta = list(
              type = "raw_upload",
              label = basename(input$raw_upload$name[[1]])
            )
          )
        } else if (identical(mode, "di")) {
          incProgress(0.25, detail = "Reading differential-inclusion table")
          uploaded <- sir_read_uploaded_table(input$di_upload, session_dir)
          if (identical(input$di_import_mode, "native")) {
            parsed <- SpliceImpactR::get_user_data_post_di(data.table::copy(uploaded))
          } else {
            event_type_column <- trimws(input$generic_event_type_column %||% "")
            column_map <- list(
              gene_id = trimws(input$generic_gene_column),
              chr = trimws(input$generic_chr_column),
              strand = trimws(input$generic_strand_column),
              inc = trimws(input$generic_inc_column),
              exc = trimws(input$generic_exc_column),
              delta_psi = trimws(input$generic_dpsi_column),
              pvalue = trimws(input$generic_pvalue_column),
              event_type = if (nzchar(event_type_column)) event_type_column else NULL
            )
            parsed <- SpliceImpactR::import_di_table(
              df = data.table::copy(uploaded),
              colmap = column_map,
              default_event_type = input$generic_default_event_type,
              add_chr_prefix = isTRUE(input$generic_add_chr)
            )
          }
          replace_di_evidence(
            parsed,
            list(
              type = "di_upload",
              label = basename(input$di_upload$name[[1]]),
              schema = input$di_import_mode
            )
          )
        } else {
          stop("Unknown evidence input mode.", call. = FALSE)
        }
        incProgress(1, detail = "Evidence ready")
      })
    })
  })

  observeEvent(input$load_features, {
    withProgress(message = "Loading protein features", value = 0, {
      run_safely("Protein feature load", {
        if (is.null(rv$reference)) {
          stop("Load a reference first.", call. = FALSE)
        }
        mode <- input$feature_mode %||% "demo"
        if (identical(mode, "demo")) {
          if (!identical(rv$reference_info$kind, "demo")) {
            stop(
              "The packaged demo features cover only the demo reference. Query or upload features for the full reference.",
              call. = FALSE
            )
          }
          incProgress(0.45, detail = "Packaged feature subset")
          feature_bundle <- sir_load_demo_features()
          set_features(
            feature_bundle$protein_feature_total,
            feature_bundle$exon_features,
            origin = "packaged demo subset"
          )
        } else {
          if (identical(mode, "query")) {
            sources <- input$feature_sources %||% character()
            if (!length(sources)) {
              stop("Choose at least one feature database.", call. = FALSE)
            }
            ensembl_release <- suppressWarnings(as.integer(
              input$biomart_release %||% rv$reference_info$ensembl_release
            ))
            if (!length(ensembl_release) || !is.finite(ensembl_release) || ensembl_release < 1L) {
              stop(
                "Enter the Ensembl release corresponding to the loaded GENCODE reference.",
                call. = FALSE
              )
            }
            feature_cache <- sir_ensure_directory(
              SIR_FEATURE_CACHE_DIR,
              "Managed protein-feature cache"
            )
            incProgress(0.25, detail = "Querying reference-matched feature sources")
            features <- SpliceImpactR::get_protein_features(
              biomaRt_databases = sources,
              gtf_df = rv$reference$annotations,
              sequences = rv$reference$sequences,
              base_dir = feature_cache,
              use_cache = TRUE,
              species = rv$reference_info$species,
              release = ensembl_release,
              test = FALSE,
              combine_overlaps = isTRUE(input$combine_feature_overlaps)
            )
            origin <- paste0(
              "Ensembl ", ensembl_release, ": ",
              paste(sources, collapse = ", ")
            )
          } else if (identical(mode, "upload")) {
            incProgress(0.2, detail = "Reading feature table")
            uploaded <- sir_read_uploaded_table(input$feature_upload, session_dir)
            if (identical(input$feature_upload_mode, "manual")) {
              features <- SpliceImpactR::get_manual_features(
                data.table::copy(uploaded),
                gtf_df = rv$reference$annotations
              )
            } else {
              features <- sir_validate_normalized_features(uploaded)
            }
            origin <- paste0(
              "uploaded ", input$feature_upload_mode, " table: ",
              basename(input$feature_upload$name[[1]])
            )
          } else {
            stop("Unknown protein feature mode.", call. = FALSE)
          }
          incProgress(0.72, detail = "Mapping features to coding exons")
          features <- data.table::as.data.table(
            SpliceImpactR::get_comprehensive_annotations(list(features))
          )
          exon_features <- SpliceImpactR::get_exon_features(
            rv$reference$annotations,
            features
          )
          set_features(features, exon_features, origin = origin)
        }
        incProgress(1, detail = "Features ready")
      })
    })
  })

  observeEvent(input$load_ppi, {
    withProgress(message = "Loading interaction network", value = 0, {
      run_safely("PPI load", {
        incProgress(
          0.35,
          detail = if (identical(input$ppi_mode, "upload")) {
            "Reading uploaded network"
          } else {
            "Reading or reusing bundled network"
          }
        )
        load_selected_ppi(force = TRUE)
        incProgress(1, detail = "Network ready")
      })
    })
  })

  observeEvent(input$build_background, {
    withProgress(message = "Building enrichment universe", value = 0, {
      run_safely("Background build", {
        if (is.null(rv$reference) || !sir_nonempty_table(rv$protein_features)) {
          stop("Load a reference and protein features first.", call. = FALSE)
        }
        source <- input$background_source
        background_input <- switch(
          source,
          annotated = character(),
          hit_index = {
            if (is.null(rv$sample_frame) || !sir_nonempty_table(rv$raw_events)) {
              stop(
                "The hit-index background requires a loaded rMATS project or demo manifest.",
                call. = FALSE
              )
            }
            conditions <- resolved_conditions()
            prepared <- sir_prepare_conditions(
              rv$raw_events,
              rv$sample_frame,
              conditions$reference,
              conditions$case
            )
            prepared$sample_frame
          },
          `user-given` = {
            transcripts <- parse_background_transcripts()
            if (!length(transcripts)) {
              stop("Enter at least one transcript ID.", call. = FALSE)
            }
            transcripts
          },
          stop("Unknown background source.", call. = FALSE)
        )
        incProgress(0.3, detail = "Matching genes and feature-bearing transcripts")
        result <- SpliceImpactR::get_background(
          source = source,
          input = background_input,
          annotations = rv$reference$annotations,
          protein_features = rv$protein_features,
          keep_annotated_first_last = isTRUE(input$keep_first_last),
          BPPARAM = BiocParallel::SerialParam()
        )
        result <- data.table::as.data.table(result)
        if (!nrow(result)) {
          stop("The selected background produced no genes.", call. = FALSE)
        }
        sir_invalidate_state(rv, "background")
        rv$background <- result
        rv$background_token <- new_token("background", source, nrow(result))
        signature <- expected_background_signature()
        sir_record_stage(
          rv,
          "background",
          signature,
          list(source = source, rows = nrow(result))
        )
        incProgress(1, detail = "Universe ready")
        append_log(
          "success",
          paste0(
            "Built a ", source, " enrichment universe with ",
            sir_fmt_int(nrow(result)), " rows."
          )
        )
      })
    })
  })

  observeEvent(input$run_di, {
    withProgress(message = "Differential inclusion", value = 0, {
      run_safely("Differential inclusion", {
        incProgress(0.15, detail = "Validating conditions and coverage")
        compute_di(force = TRUE)
        incProgress(1, detail = "Differential inclusion ready")
      })
    })
  })

  observeEvent(input$map_events, {
    withProgress(message = "Matching significant events", value = 0, {
      run_safely("Annotation matching", {
        incProgress(0.2, detail = "Checking differential inclusion")
        compute_di(force = FALSE)
        incProgress(0.55, detail = "Matching events to transcripts")
        compute_mapping(force = TRUE)
        incProgress(1, detail = "Annotation matching ready")
      })
    })
  })

  observeEvent(input$run_all_pipeline, {
    withProgress(message = "Running full pipeline", value = 0, {
      run_safely("Full pipeline", {
        if (
          !sir_nonempty_table(rv$protein_features) ||
          !sir_nonempty_table(rv$exon_features)
        ) {
          stop(
            "Run Full Pipeline requires loaded protein features. ",
            "Load packaged, BioMart, or uploaded features first.",
            call. = FALSE
          )
        }
        incProgress(0.08, detail = "Loading the selected PPI network")
        load_selected_ppi(force = FALSE)
        incProgress(0.18, detail = "Checking differential inclusion")
        compute_di(force = FALSE)
        incProgress(0.38, detail = "Matching events to transcripts")
        compute_mapping(force = FALSE)
        incProgress(0.65, detail = "Comparing sequences, domains, and interactions")
        compute_consequences(force = FALSE)
        if (
          !identical(rv$final_kind, "sequence + domains + PPI") ||
          !"n_ppi" %in% names(rv$final_results)
        ) {
          stop(
            "The full pipeline did not produce PPI-switch results.",
            call. = FALSE
          )
        }
        incProgress(1, detail = "Sequence, domain, and PPI stages complete")
      })
    })
  })

  observeEvent(input$run_analysis, {
    withProgress(message = "Running current analysis", value = 0, {
      run_safely("Analysis", {
        incProgress(0.1, detail = "Checking differential inclusion")
        compute_di(force = FALSE)
        incProgress(0.35, detail = "Matching events to transcripts")
        compute_mapping(force = FALSE)
        incProgress(0.62, detail = "Comparing transcript sequences")
        compute_consequences(force = FALSE)
        incProgress(1, detail = "Current stages complete")
      })
    })
  })

  observeEvent(input$run_consequences, {
    withProgress(message = "Refreshing consequences", value = 0, {
      run_safely("Consequence refresh", {
        incProgress(0.2, detail = "Checking upstream provenance")
        compute_consequences(force = FALSE)
        incProgress(1, detail = "Consequences current")
      })
    })
  })

  expected_enrichment_signature <- reactive({
    background <- current_background()
    if (is.null(background)) {
      return(NULL)
    }
    mode <- input$enrichment_mode %||% "di"
    foreground_signature <- switch(
      mode,
      di = if (!is.null(current_di())) {
        sir_signature(
          sir_stage_signature(rv, "di"),
          sir_threshold_parameters(input)
        )
      } else {
        NULL
      },
      domain = if (!is.null(current_domains())) {
        sir_stage_signature(rv, "domain_results")
      } else {
        NULL
      },
      ppi = if (
        !is.null(current_final()) &&
        "n_ppi" %in% names(current_final())
      ) {
        sir_stage_signature(rv, "final_results")
      } else {
        NULL
      },
      NULL
    )
    if (is.null(foreground_signature)) {
      return(NULL)
    }
    sir_signature(
      "enrichment-v1",
      mode,
      foreground_signature,
      sir_stage_signature(rv, "background"),
      input$enrichment_gene_id_type,
      sort(input$enrichment_sources %||% character()),
      as.integer(input$enrichment_min_size),
      as.integer(input$enrichment_max_size),
      as.numeric(input$enrichment_fdr),
      input$enrichment_plot_type
    )
  })

  current_enrichment <- reactive({
    expected <- expected_enrichment_signature()
    if (
      is.null(rv$enrichment) ||
      is.null(expected) ||
      !identical(sir_stage_signature(rv, "enrichment"), expected)
    ) {
      return(NULL)
    }
    rv$enrichment
  })

  observeEvent(input$run_enrichment, {
    withProgress(message = "Running enrichment", value = 0, {
      run_safely("Enrichment", {
        background <- current_background()
        if (is.null(background)) {
          stop(
            "Build or refresh the enrichment background before running enrichment.",
            call. = FALSE
          )
        }
        if (!length(input$enrichment_sources)) {
          stop("Choose at least one enrichment collection.", call. = FALSE)
        }
        if (input$enrichment_min_size >= input$enrichment_max_size) {
          stop("Maximum set size must be larger than minimum set size.", call. = FALSE)
        }

        mode <- input$enrichment_mode
        incProgress(0.2, detail = "Preparing foreground")
        if (identical(mode, "di")) {
          compute_di(force = FALSE)
          foreground <- SpliceImpactR::get_gene_enrichment(
            mode = "di",
            res = current_di(),
            padj_threshold = input$padj_threshold,
            delta_psi_threshold = input$dpsi_threshold
          )
        } else {
          compute_consequences(force = FALSE)
          if (identical(mode, "domain")) {
            domains <- current_domains()
            if (is.null(domains)) {
              stop(
                "Domain enrichment requires loaded protein features and current domain results.",
                call. = FALSE
              )
            }
            foreground <- SpliceImpactR::get_gene_enrichment(
              mode = "domain",
              hits = domains
            )
          } else {
            final <- current_final()
            if (is.null(final) || !"n_ppi" %in% names(final)) {
              stop(
                "PPI enrichment requires a loaded interaction network and current PPI results.",
                call. = FALSE
              )
            }
            foreground <- SpliceImpactR::get_gene_enrichment(
              mode = "ppi",
              hits = final
            )
          }
        }
        foreground <- unique(na.omit(as.character(foreground)))
        if (!length(foreground)) {
          stop("The selected foreground contains no genes.", call. = FALSE)
        }
        if (!"gene_id" %in% names(background)) {
          stop("The current background does not contain a gene_id column.", call. = FALSE)
        }

        incProgress(0.58, detail = "Querying gene-set collections")
        result <- SpliceImpactR::get_enrichment(
          foreground = foreground,
          background = unique(na.omit(as.character(background$gene_id))),
          species = rv$reference_info$species,
          gene_id_type = input$enrichment_gene_id_type,
          sources = input$enrichment_sources,
          min_size = as.integer(input$enrichment_min_size),
          max_size = as.integer(input$enrichment_max_size),
          p_adjust_cutoff = as.numeric(input$enrichment_fdr),
          plot_type = input$enrichment_plot_type
        )
        rv$enrichment <- result
        signature <- expected_enrichment_signature()
        significant_rows <- tryCatch(
          nrow(data.table::as.data.table(result$results_signif)),
          error = function(error) 0L
        )
        sir_record_stage(
          rv,
          "enrichment",
          signature,
          list(
            mode = mode,
            foreground_genes = length(foreground),
            background_genes = data.table::uniqueN(background$gene_id),
            significant_sets = significant_rows,
            sources = input$enrichment_sources,
            adjusted_p_cutoff = input$enrichment_fdr
          )
        )
        incProgress(1, detail = "Enrichment ready")
        append_log(
          "success",
          paste0(
            "Enrichment completed for ", sir_fmt_int(length(foreground)),
            " foreground genes; ", sir_fmt_int(significant_rows),
            " gene sets pass the selected cutoff."
          )
        )
      })
    })
  })

  observeEvent(input$reference_condition, {
    counts <- raw_condition_counts()
    if (!nrow(counts)) {
      return()
    }
    choices <- setdiff(counts$condition, input$reference_condition %||% "")
    if (!length(choices)) {
      return()
    }
    selection <- sir_resolve_condition_selection(
      rv$raw_events,
      input$reference_condition,
      isolate(input$case_condition)
    )
    updateSelectInput(
      session,
      "case_condition",
      choices = choices,
      selected = selection$case
    )
  }, ignoreInit = TRUE)

  annotation_index <- reactive({
    if (is.null(rv$reference) || is.null(rv$reference_info$id)) {
      return(data.table::data.table())
    }
    sir_reference_index(rv$reference_info$id, rv$reference)
  })

  observeEvent(rv$raw_events, {
    genes <- if (sir_nonempty_table(rv$raw_events)) {
      sort(unique(na.omit(as.character(rv$raw_events$gene_id))))
    } else {
      character()
    }
    updateSelectizeInput(
      session,
      "probe_gene",
      choices = genes,
      selected = if (length(genes)) genes[[1]] else character(),
      server = TRUE
    )
  }, ignoreInit = FALSE)

  probe_event_types <- reactive({
    if (!sir_nonempty_table(rv$raw_events) || !nzchar(input$probe_gene %||% "")) {
      return(character())
    }
    sort(unique(as.character(
      data.table::as.data.table(rv$raw_events)[gene_id == input$probe_gene, event_type]
    )))
  })

  probe_event_ids <- reactive({
    if (!sir_nonempty_table(rv$raw_events) || !nzchar(input$probe_gene %||% "")) {
      return(character())
    }
    dt <- data.table::as.data.table(rv$raw_events)[gene_id == input$probe_gene]
    selected_event_type <- input$probe_event_type %||% ""
    if (nzchar(selected_event_type)) {
      dt <- sir_filter_event_type(dt, selected_event_type)
    }
    sort(unique(as.character(dt$event_id)))
  })

  overview_event_types <- reactive({
    sir_present_event_types(rv$raw_events)
  })

  observeEvent(overview_event_types(), {
    choices <- overview_event_types()
    current <- isolate(input$overview_event_type %||% "")
    selected <- if (length(current) && current %in% choices) {
      current
    } else if (length(choices)) {
      choices[[1]]
    } else {
      character()
    }
    updateSelectInput(
      session,
      "overview_event_type",
      choices = sir_event_type_choices(choices),
      selected = selected
    )
  }, ignoreInit = FALSE)

  observeEvent(
    list(input$probe_gene, input$probe_event_type, input$probe_event),
    {
      rv$probe <- NULL
    },
    ignoreInit = TRUE
  )

  observeEvent(
    list(
      input$overview_event_type,
      input$reference_condition,
      input$case_condition
    ),
    {
      rv$hit_compare <- NULL
      rv$psi_overview <- NULL
    },
    ignoreInit = TRUE
  )

  observeEvent(annotation_index(), {
    index <- annotation_index()
    genes <- if (nrow(index)) {
      sort(unique(na.omit(as.character(index$gene_id))))
    } else {
      character()
    }
    updateSelectizeInput(
      session,
      "browser_gene",
      choices = genes,
      selected = if (length(genes)) genes[[1]] else character(),
      server = TRUE
    )
  }, ignoreInit = FALSE)

  observeEvent(list(annotation_index(), input$browser_gene), {
    index <- annotation_index()
    gene <- input$browser_gene %||% ""
    transcripts <- if (nrow(index) && nzchar(gene)) {
      sort(unique(na.omit(as.character(index[gene_id == gene, transcript_id]))))
    } else {
      character()
    }
    selected_a <- if (length(transcripts)) transcripts[[1]] else character()
    selected_b <- if (length(transcripts) > 1L) transcripts[[2]] else selected_a
    updateSelectInput(session, "transcript_a", choices = transcripts, selected = selected_a)
    updateSelectInput(session, "transcript_b", choices = transcripts, selected = selected_b)
  }, ignoreInit = FALSE)

  observeEvent(rv$protein_features, {
    databases <- if (
      sir_nonempty_table(rv$protein_features) &&
      "database" %in% names(rv$protein_features)
    ) {
      sort(unique(na.omit(as.character(rv$protein_features$database))))
    } else {
      character()
    }
    updateSelectizeInput(
      session,
      "browser_feature_databases",
      choices = databases,
      selected = head(databases, 4L),
      server = TRUE
    )
  }, ignoreInit = FALSE)

  observeEvent(input$run_probe, {
    withProgress(message = "Probing event", value = 0, {
      run_safely("Event probe", {
        if (!sir_nonempty_table(rv$raw_events)) {
          stop("Load raw events first.", call. = FALSE)
        }
        event <- input$probe_event %||% ""
        if (!nzchar(event)) {
          stop("Choose an event to probe.", call. = FALSE)
        }
        rv$probe <- SpliceImpactR::probe_individual_event(
          rv$raw_events,
          event = event
        )
        incProgress(1)
        append_log("success", paste0("Event probe ready for ", event, "."))
      })
    })
  })

  observeEvent(input$run_browser, {
    withProgress(message = "Comparing transcripts", value = 0, {
      run_safely("Transcript comparison", {
        if (is.null(rv$reference) || !sir_nonempty_table(rv$protein_features)) {
          stop("Load a reference and protein features first.", call. = FALSE)
        }
        transcript_a <- input$transcript_a %||% ""
        transcript_b <- input$transcript_b %||% ""
        if (!nzchar(transcript_a) || !nzchar(transcript_b)) {
          stop("Choose two transcripts.", call. = FALSE)
        }
        if (identical(transcript_a, transcript_b)) {
          stop("Choose two different transcripts.", call. = FALSE)
        }
        request <- data.frame(
          transcript1 = transcript_a,
          transcript2 = transcript_b,
          stringsAsFactors = FALSE
        )
        matched <- SpliceImpactR::compare_transcript_pairs(
          request,
          rv$reference$annotations
        )
        if (!nrow(matched)) {
          stop("The selected transcript pair could not be built.", call. = FALSE)
        }
        sequences <- SpliceImpactR::attach_sequences(
          matched,
          rv$reference$sequences
        )
        pairs <- SpliceImpactR::get_pairs(sequences, source = "multi")
        sequence <- SpliceImpactR::compare_sequence_frame(
          pairs,
          rv$reference$annotations
        )
        domains <- SpliceImpactR::get_domains(
          sequence,
          rv$exon_features
        )
        plot <- SpliceImpactR::plot_two_transcripts_with_domains_unified(
          transcripts = c(transcript_a, transcript_b),
          gtf_df = rv$reference$annotations,
          protein_features = rv$protein_features,
          feature_db = if (length(input$browser_feature_databases)) {
            input$browser_feature_databases
          } else {
            NULL
          },
          combine_domains = isTRUE(input$browser_combine_domains),
          view = input$browser_view
        )
        rv$browser <- list(
          sequence = data.table::as.data.table(sequence),
          domains = data.table::as.data.table(domains),
          plot = plot,
          transcripts = c(transcript_a, transcript_b)
        )
        incProgress(1)
        append_log(
          "success",
          paste0("Compared ", transcript_a, " with ", transcript_b, ".")
        )
      })
    })
  })

  observeEvent(input$run_overview, {
    rv$hit_compare <- NULL
    rv$psi_overview <- NULL
    withProgress(message = "Building sample overview", value = 0, {
      run_safely("Sample overview", {
        if (is.null(rv$sample_frame) || !sir_nonempty_table(rv$raw_events)) {
          stop(
            "A sample-level overview requires a demo or project manifest.",
            call. = FALSE
          )
        }
        conditions <- resolved_conditions()
        prepared <- sir_prepare_conditions(
          rv$raw_events,
          rv$sample_frame,
          conditions$reference,
          conditions$case
        )
        hit_compare <- SpliceImpactR::compare_hit_index(
          prepared$sample_frame,
          condition_map = c(control = "control", test = "case")
        )
        selected_event_type <- input$overview_event_type %||% ""
        event_rows <- sir_filter_event_type(
          prepared$raw_events,
          selected_event_type
        )
        if (!nrow(event_rows)) {
          stop("No raw rows exist for the selected event type.", call. = FALSE)
        }
        psi_overview <- SpliceImpactR::overview_spicing_comparison(
          events = event_rows,
          sample_df = prepared$sample_frame,
          depth_norm = "exon_files",
          event_type = selected_event_type,
          conditions = c(control = "control", experimental = "case")
        )
        rv$hit_compare <- hit_compare
        rv$psi_overview <- psi_overview
        incProgress(1)
        append_log(
          "success",
          paste0(
            "Sample-level HIT and ", selected_event_type,
            " PSI overview is ready."
          )
        )
      })
    })
  })

  output$runtime_banner <- renderUI({
    description <- tryCatch(
      utils::packageDescription("SpliceImpactR"),
      error = function(error) NULL
    )
    installed_sha <- description$RemoteSha %||% ""
    if (nzchar(installed_sha) && !identical(installed_sha, SIR_PACKAGE_SHA)) {
      return(div(
        class = "sir-runtime-banner",
        "Deployment warning: the installed SpliceImpactR commit does not match the app's expected commit."
      ))
    }
    NULL
  })

  output$reference_source_controls <- renderUI({
    mode <- input$annotation_mode %||% sir_default_reference_mode()
    if (identical(mode, "upload")) {
      urls <- sir_gencode_urls("human", 45L)
      return(tagList(
        sir_panel_note(
          "Load from this computer",
          paste0(
            "Upload the three official human GENCODE v45 source files. ",
            "They are copied into an isolated session workspace and removed ",
            "when the session ends."
          ),
          "info"
        ),
        fileInput(
          "annotation_gtf_upload",
          "Annotation GTF · gencode.v45.annotation.gtf.gz",
          accept = c(".gz", "application/gzip", "application/x-gzip"),
          buttonLabel = "Choose GTF"
        ),
        fileInput(
          "annotation_transcripts_upload",
          "Transcript FASTA · gencode.v45.pc_transcripts.fa.gz",
          accept = c(".gz", "application/gzip", "application/x-gzip"),
          buttonLabel = "Choose transcript FASTA"
        ),
        fileInput(
          "annotation_translations_upload",
          "Protein FASTA · gencode.v45.pc_translations.fa.gz",
          accept = c(".gz", "application/gzip", "application/x-gzip"),
          buttonLabel = "Choose protein FASTA"
        ),
        checkboxGroupInput(
          "annotation_tsl",
          "Transcript support levels",
          choices = stats::setNames(
            as.character(1:5),
            paste("TSL", 1:5)
          ),
          selected = c("1", "2", "3"),
          inline = TRUE
        ),
        tags$div(
          class = "reference-download-links",
          tags$strong("Need the official files first?"),
          tags$p(
            class = "help-note",
            "Download them directly from GENCODE to this computer, then upload all three above."
          ),
          tags$a(
            class = "btn btn-default btn-sm",
            href = urls[["gtf"]],
            target = "_blank",
            rel = "noopener noreferrer",
            "Download annotation GTF"
          ),
          tags$a(
            class = "btn btn-default btn-sm",
            href = urls[["transcripts"]],
            target = "_blank",
            rel = "noopener noreferrer",
            "Download transcript FASTA"
          ),
          tags$a(
            class = "btn btn-default btn-sm",
            href = urls[["translations"]],
            target = "_blank",
            rel = "noopener noreferrer",
            "Download protein FASTA"
          )
        )
      ))
    }
    if (identical(mode, "library")) {
      catalog <- sir_reference_catalog()
      catalog <- catalog[setdiff(names(catalog), "demo_human_v45")]
      if (!length(catalog)) {
        return(sir_panel_note(
          "No server-library reference found",
          "Set SPLICEIMPACTR_REFERENCE_DIRS to one or more trusted directories containing a complete preprocessed bundle or matching GTF/FASTA triplet.",
          "warning"
        ))
      }
      choices <- stats::setNames(
        names(catalog),
        vapply(catalog, `[[`, character(1), "label")
      )
      return(selectInput(
        "reference_library_id",
        "Available server reference",
        choices = choices,
        selected = names(choices)[[1]]
      ))
    }
    if (mode %in% c("cached", "download")) {
      return(tagList(
        selectInput(
          "annotation_species",
          "Species",
          choices = c("Human" = "human", "Mouse" = "mouse"),
          selected = "human"
        ),
        numericInput(
          "annotation_release",
          "GENCODE release",
          value = 45,
          min = 1,
          step = 1
        ),
        checkboxGroupInput(
          "annotation_tsl",
          "Transcript support levels",
          choices = stats::setNames(as.character(1:5), paste("TSL", 1:5)),
          selected = c("1", "2", "3"),
          inline = TRUE
        )
      ))
    }
    NULL
  })

  output$reference_description <- renderUI({
    info <- tryCatch(
      reference_request_from_input(),
      error = function(error) NULL
    )
    if (is.null(info)) {
      return(NULL)
    }
    ensembl_text <- if (
      length(info$ensembl_release) &&
      is.finite(info$ensembl_release)
    ) {
      paste0(" maps to Ensembl ", info$ensembl_release)
    } else {
      "; enter the matching Ensembl release under Protein Features"
    }
    sir_panel_note(
      paste0(info$species_label, " · ", info$assembly),
      paste0(
        "GENCODE ", sir_release_tag(info$species, info$gencode_release),
        ensembl_text, ". ",
        info$note
      ),
      tone = if (info$kind %in% c("demo", "download")) "warning" else "info"
    )
  })

  output$reference_storage_status <- renderUI({
    storage <- sir_reference_storage_summary()
    if (isTRUE(storage$hosted)) {
      return(NULL)
    }
    display_path <- function(path) {
      if (isTRUE(storage$hosted)) basename(path) else path
    }
    roots <- vapply(storage$roots, display_path, character(1))
    sir_panel_note(
      "Declared storage",
      paste0(
        "Library roots: ",
        if (length(roots)) paste(roots, collapse = ", ") else "none",
        ". Managed cache: ", display_path(storage$cache),
        ". GENCODE downloads: ",
        if (storage$downloads_enabled) "enabled" else "disabled by server policy",
        "."
      ),
      "info"
    )
  })

  output$input_mode_controls <- renderUI({
    mode <- input$input_mode %||% "demo"
    if (mode %in% c("demo", "project")) {
      controls <- tagList(
        checkboxGroupInput(
          "event_types",
          "Event types",
          choices = SIR_SUPPORTED_EVENTS,
          selected = SIR_SUPPORTED_EVENTS,
          inline = TRUE
        ),
        checkboxInput(
          "keep_first_last",
          "Keep annotated first/last exon events",
          value = TRUE
        )
      )
      if (identical(mode, "demo")) {
        return(tagList(
          sir_panel_note(
            "Self-contained example",
            "Eight packaged samples (four control, four case) are loaded with JCEC counts.",
            "info"
          ),
          controls
        ))
      }
      return(tagList(
        fileInput(
          "project_archive",
          "Project ZIP",
          accept = ".zip",
          buttonLabel = "Choose ZIP"
        ),
        radioButtons(
          "rmats_use",
          "rMATS count table",
          choices = c("Junction + exon counts" = "JCEC", "Junction counts" = "JC"),
          selected = "JCEC",
          inline = TRUE
        ),
        controls,
        sir_help_text(
          "The archive must contain one manifest.csv/tsv/txt with relative path, sample_name, and condition columns, plus every sample directory."
        )
      ))
    }
    if (identical(mode, "raw")) {
      return(tagList(
        fileInput(
          "raw_upload",
          "Normalized raw-event table",
          accept = c(".csv", ".tsv", ".txt", ".csv.gz", ".tsv.gz", ".txt.gz")
        ),
        sir_help_text(
          "Required: event_id, form, gene_id, chr, strand, inc, exc, inclusion_reads, exclusion_reads, sample, and condition. Each selected condition needs at least two distinct samples."
        )
      ))
    }
    tagList(
      fileInput(
        "di_upload",
        "Differential-inclusion table",
        accept = c(".csv", ".tsv", ".txt", ".csv.gz", ".tsv.gz", ".txt.gz")
      ),
      radioButtons(
        "di_import_mode",
        "Schema",
        choices = c(
          "SpliceImpactR / native event form" = "native",
          "Generic site-level table" = "generic"
        ),
        selected = "native"
      ),
      conditionalPanel(
        "input.di_import_mode === 'generic'",
        div(
          class = "sir-generic-grid",
          textInput("generic_gene_column", "Gene column", "gene_id"),
          textInput("generic_chr_column", "Chromosome column", "chr"),
          textInput("generic_strand_column", "Strand column", "strand"),
          textInput("generic_inc_column", "Included coordinates", "inc"),
          textInput("generic_exc_column", "Excluded coordinates", "exc"),
          textInput("generic_dpsi_column", "ΔPSI column", "delta_psi"),
          textInput("generic_pvalue_column", "p-value column", "p.value"),
          textInput("generic_event_type_column", "Event-type column (optional)", ""),
          selectInput(
            "generic_default_event_type",
            "Default event type",
            choices = c(SIR_SUPPORTED_EVENTS, "SITE"),
            selected = "SITE"
          ),
          checkboxInput("generic_add_chr", "Add chr prefix", FALSE)
        )
      ),
      sir_help_text(
        "Imported DI is treated as evidence. Downstream and enrichment runs reuse it; they never attempt to rerun a missing raw-event stage."
      )
    )
  })

  output$condition_controls <- renderUI({
    counts <- raw_condition_counts()
    if (!nrow(counts)) {
      return(sir_panel_note(
        "Not required for imported DI",
        if (sir_nonempty_table(rv$di)) {
          "The uploaded DI table already encodes the comparison. Its ΔPSI direction is preserved."
        } else {
          "Load raw events to choose the reference and case levels."
        },
        "info"
      ))
    }
    conditions <- counts$condition
    selection <- sir_resolve_condition_selection(
      rv$raw_events,
      isolate(input$reference_condition),
      isolate(input$case_condition)
    )

    tagList(
      div(
        class = "sir-condition-grid",
        selectInput(
          "reference_condition",
          "Reference / baseline",
          choices = conditions,
          selected = selection$reference
        ),
        selectInput(
          "case_condition",
          "Case / comparison",
          choices = conditions,
          selected = selection$case
        )
      ),
      div(
        class = "sir-condition-counts",
        lapply(seq_len(nrow(counts)), function(index) {
          span(
            class = "sir-condition-chip",
            paste0(counts$condition[[index]], " · ", counts$samples[[index]], " samples")
          )
        })
      ),
      if (length(conditions) > 2L) {
        sir_help_text("Only the two selected levels enter the model; all other conditions are excluded explicitly.")
      }
    )
  })

  output$workspace_status <- renderUI({
    di_ready <- !is.null(current_di())
    matched_ready <- !is.null(current_matched())
    consequence_ready <- !is.null(current_final())
    tagList(
      sir_status_row(
        "Reference",
        !is.null(rv$reference),
        if (is.null(rv$reference)) "Choose and load a deployed bundle" else rv$reference_info$label
      ),
      sir_status_row(
        "Evidence",
        sir_nonempty_table(rv$raw_events) || sir_nonempty_table(rv$di),
        if (sir_nonempty_table(rv$raw_events)) {
          paste0(
            sir_fmt_int(data.table::uniqueN(rv$raw_events$event_id)),
            " raw events"
          )
        } else if (sir_nonempty_table(rv$di)) {
          paste0(sir_fmt_int(data.table::uniqueN(rv$di$event_id)), " imported DI events")
        } else {
          "Load demo, project, raw, or DI evidence"
        }
      ),
      sir_status_row(
        "Differential inclusion",
        di_ready,
        if (di_ready) {
          paste0(
            sir_fmt_int(nrow(current_di())), " rows · ",
            rv$di_origin
          )
        } else if (sir_nonempty_table(rv$di)) {
          "Parameters or condition meaning changed"
        } else {
          "Not run"
        },
        dirty = sir_nonempty_table(rv$di) && !di_ready
      ),
      sir_status_row(
        "Transcript matching",
        matched_ready,
        if (matched_ready) {
          paste0(
            sir_fmt_int(data.table::uniqueN(current_matched()$event_id)),
            " events mapped"
          )
        } else if (sir_nonempty_table(rv$matched)) {
          "Thresholds or upstream data changed"
        } else {
          "Not run"
        },
        dirty = sir_nonempty_table(rv$matched) && !matched_ready
      ),
      sir_status_row(
        "Consequences",
        consequence_ready,
        if (consequence_ready) rv$final_kind else if (sir_nonempty_table(rv$final_results)) {
          "Upstream resource or pairing changed"
        } else {
          "Not run"
        },
        dirty = sir_nonempty_table(rv$final_results) && !consequence_ready
      )
    )
  })

  output$analysis_freshness <- renderUI({
    if (sir_nonempty_table(rv$di) && is.null(current_di())) {
      return(sir_panel_note(
        "Differential inclusion is out of date",
        "Condition meaning or model settings changed. Run the current analysis to recompute it before viewing or exporting results.",
        "warning"
      ))
    }
    if (sir_nonempty_table(rv$matched) && is.null(current_matched())) {
      return(sir_panel_note(
        "Downstream stages are out of date",
        "The significance thresholds changed. The next analysis run will remap significant events before any consequence step.",
        "warning"
      ))
    }
    if (!is.null(current_final())) {
      return(sir_panel_note(
        "Current",
        paste0("All completed results match the active parameters (", rv$final_kind, ")."),
        "info"
      ))
    }
    NULL
  })

  output$analysis_message <- renderUI({
    if (is.null(rv$reference)) {
      return(sir_panel_note("Reference required", "Load a reference on the Set up page.", "warning"))
    }
    if (!sir_nonempty_table(rv$raw_events) && !sir_nonempty_table(rv$di)) {
      return(sir_panel_note("Evidence required", "Load raw events or a DI table on the Set up page.", "warning"))
    }
    if (!sir_nonempty_table(rv$protein_features)) {
      return(sir_panel_note(
        "Sequence analysis is available",
        "Protein features are not loaded, so the run will stop after sequence consequences. Add them on Functional data for domains.",
        "info"
      ))
    }
    if (!sir_nonempty_table(rv$ppi)) {
      return(sir_panel_note(
        "Domain analysis is available",
        "No PPI network is loaded. Interaction consequences will be skipped without silently loading a large network.",
        "info"
      ))
    }
    NULL
  })

  output$consequence_message <- renderUI({
    if (is.null(current_sequence())) {
      return(sir_panel_note(
        "Run the analysis first",
        "Current sequence-pair results are required before domain or PPI consequences can be shown.",
        "warning"
      ))
    }
    if (is.null(current_domains())) {
      return(sir_panel_note(
        "No current domain layer",
        "Load reference-matched protein features, then refresh consequences.",
        "warning"
      ))
    }
    if (is.null(rv$ppi)) {
      return(sir_panel_note(
        "PPI layer not loaded",
        "Domain results are current; load an interaction network to calculate PPI changes.",
        "info"
      ))
    }
    NULL
  })

  output$metric_cards <- renderUI({
    annotation_genes <- if (!is.null(rv$reference)) {
      data.table::uniqueN(
        data.table::as.data.table(rv$reference$annotations)$gene_id,
        na.rm = TRUE
      )
    } else {
      0L
    }
    raw_events <- if (sir_nonempty_table(rv$raw_events)) {
      data.table::uniqueN(rv$raw_events$event_id)
    } else {
      0L
    }
    significant_events <- if (nrow(significant_di())) {
      data.table::uniqueN(significant_di()$event_id)
    } else {
      0L
    }
    final_pairs <- if (!is.null(current_final())) nrow(current_final()) else 0L
    div(
      class = "sir-metric-grid",
      sir_metric_card(
        "Reference genes",
        sir_fmt_int(annotation_genes),
        if (is.null(rv$reference)) "Reference not loaded" else rv$reference_info$short_label,
        "teal"
      ),
      sir_metric_card(
        "Raw events",
        sir_fmt_int(raw_events),
        if (raw_events) "Unique event_id values" else "Raw events not loaded",
        "blue"
      ),
      sir_metric_card(
        "Significant events",
        sir_fmt_int(significant_events),
        paste0("padj ≤ ", input$padj_threshold, " · |ΔPSI| ≥ ", input$dpsi_threshold),
        "rust"
      ),
      sir_metric_card(
        "Current pairs",
        sir_fmt_int(final_pairs),
        if (!is.null(current_final())) {
          rv$final_kind
        } else if (sir_nonempty_table(rv$final_results)) {
          "Out of date · rerun analysis"
        } else {
          "Consequences not run"
        },
        "green"
      )
    )
  })

  input_preview_data <- reactive({
    query <- input$input_search %||% ""
    if (sir_nonempty_table(rv$raw_events)) {
      columns <- intersect(
        c(
          "event_id", "event_type", "form", "gene_id", "chr", "strand",
          "inc", "exc", "inclusion_reads", "exclusion_reads", "psi",
          "sample", "condition"
        ),
        names(rv$raw_events)
      )
      return(sir_preview_data(
        rv$raw_events,
        preview_rows(),
        query,
        columns = columns,
        search_columns = c("event_id", "event_type", "gene_id", "sample", "condition")
      ))
    }
    if (sir_nonempty_table(rv$di)) {
      return(sir_preview_data(
        rv$di,
        preview_rows(),
        query,
        search_columns = c("event_id", "event_type", "gene_id", "chr")
      ))
    }
    if (!is.null(rv$sample_frame)) {
      manifest <- data.table::as.data.table(rv$sample_frame)
      if ("path" %in% names(manifest)) {
        manifest <- data.table::copy(manifest)
        manifest[, path := basename(sub("/+$", "", path))]
      }
      return(sir_preview_data(manifest, preview_rows(), query))
    }
    sir_preview_data(NULL)
  })

  output$input_preview <- DT::renderDT({
    sir_dt_widget(input_preview_data(), page_length = preview_rows())
  })

  output$input_preview_caption <- renderUI({
    preview <- input_preview_data()
    if (!preview$total) "No input is loaded." else sir_preview_caption(preview)
  })

  di_preview_data <- reactive({
    sir_preview_data(
      significant_di(),
      preview_rows(),
      input$analysis_search %||% "",
      search_columns = c("event_id", "event_type", "gene_id", "chr")
    )
  })

  output$di_table <- DT::renderDT({
    sir_dt_widget(di_preview_data(), page_length = preview_rows())
  })

  output$di_table_caption <- renderUI({
    preview <- di_preview_data()
    if (!preview$total) "No significant, current DI rows." else sir_preview_caption(preview)
  })

  output$di_plot <- renderPlot({
    di <- current_di()
    if (is.null(di)) {
      return(NULL)
    }
    SpliceImpactR::plot_di_volcano_dt(
      di,
      padj_thr = input$padj_threshold,
      dpsi_thr = input$dpsi_threshold
    )
  })

  output$di_plot_note <- renderUI({
    di <- current_di()
    if (is.null(di)) {
      return(sir_panel_note("Waiting for current DI", "Run differential inclusion or import a DI table.", "info"))
    }
    significant <- significant_di()
    sir_panel_note(
      "How to read this",
      paste0(
        "Each point is one event form. ",
        sir_fmt_int(data.table::uniqueN(significant$event_id)),
        " unique events pass the active cutoffs out of ",
        sir_fmt_int(data.table::uniqueN(di$event_id)), "."
      ),
      "info"
    )
  })

  matched_preview_data <- reactive({
    sir_preview_data(
      current_matched(),
      preview_rows(),
      input$analysis_search %||% "",
      search_columns = c("event_id", "event_type", "gene_id", "transcript_id", "protein_id")
    )
  })

  output$matched_table <- DT::renderDT({
    sir_dt_widget(matched_preview_data(), page_length = preview_rows())
  })

  output$matched_table_caption <- renderUI({
    preview <- matched_preview_data()
    if (!preview$total) "No current event mapping." else sir_preview_caption(preview)
  })

  output$mapping_summary <- renderUI({
    matched <- current_matched()
    if (is.null(matched)) {
      return(sir_panel_note(
        "Not available",
        "Run the current analysis after loading a reference and evidence.",
        "info"
      ))
    }
    summary <- data.table::as.data.table(matched)[
      , .(
        events = data.table::uniqueN(event_id),
        transcripts = data.table::uniqueN(transcript_id),
        rows = .N
      ),
      by = event_type
    ][order(-events)]
    tags$table(
      class = "table table-sm sir-summary-table",
      tags$thead(tags$tr(lapply(names(summary), tags$th))),
      tags$tbody(lapply(seq_len(nrow(summary)), function(index) {
        tags$tr(lapply(summary[index], function(value) tags$td(as.character(value))))
      }))
    )
  })

  output$alignment_plot <- renderPlot({
    sequence <- current_sequence()
    if (is.null(sequence)) return(NULL)
    SpliceImpactR::plot_alignment_summary(sequence, mode = "protein")
  })

  output$length_plot <- renderPlot({
    sequence <- current_sequence()
    if (is.null(sequence)) return(NULL)
    SpliceImpactR::plot_length_comparison(sequence, mode = "protein")
  })

  output$sequence_note <- renderUI({
    sequence <- current_sequence()
    if (is.null(sequence)) {
      return(sir_panel_note("No current sequence results", "Run the current analysis.", "info"))
    }
    protein_identity <- if ("prot_pid" %in% names(sequence)) {
      suppressWarnings(stats::median(as.numeric(sequence$prot_pid), na.rm = TRUE))
    } else {
      NA_real_
    }
    sir_panel_note(
      "Current sequence layer",
      paste0(
        sir_fmt_int(nrow(sequence)), " transcript pairs across ",
        sir_fmt_int(data.table::uniqueN(sequence$event_id)), " events",
        if (is.finite(protein_identity)) {
          paste0("; median protein identity is ", sir_fmt_num(protein_identity, 1, "%"), ".")
        } else {
          "."
        }
      ),
      "info"
    )
  })

  sequence_preview_data <- reactive({
    sir_preview_data(
      current_sequence(),
      preview_rows(),
      input$analysis_search %||% "",
      search_columns = c(
        "event_id", "event_type", "gene_id", "transcript_id_control",
        "transcript_id_case", "protein_id_control", "protein_id_case"
      )
    )
  })

  output$sequence_table <- DT::renderDT({
    sir_dt_widget(sequence_preview_data(), page_length = preview_rows())
  })

  output$sequence_table_caption <- renderUI({
    preview <- sequence_preview_data()
    if (!preview$total) "No current sequence results." else sir_preview_caption(preview)
  })

  output$feature_coverage <- renderUI({
    if (!sir_nonempty_table(rv$protein_features)) {
      return(sir_panel_note("No feature set loaded", "Domain analysis is not yet available.", "info"))
    }
    coverage <- sir_feature_coverage(rv$protein_features, rv$reference)
    sir_panel_note(
      "Reference coverage",
      paste0(
        sir_fmt_int(coverage$feature_transcripts), " of ",
        sir_fmt_int(coverage$reference_transcripts), " transcripts (",
        sir_fmt_num(100 * coverage$proportion, 1, "%"), ")."
      ),
      if (coverage$proportion < 0.2 && !identical(rv$reference_info$kind, "demo")) "warning" else "info"
    )
  })

  output$ppi_status <- renderUI({
    if (!sir_nonempty_table(rv$ppi)) {
      return(sir_panel_note("No network loaded", "PPI consequences will be skipped.", "info"))
    }
    sir_panel_note(
      "Network ready",
      paste0(sir_fmt_int(nrow(rv$ppi)), " interaction rows are available to this session."),
      "info"
    )
  })

  output$background_status <- renderUI({
    background <- current_background()
    if (!is.null(background)) {
      return(sir_panel_note(
        "Universe current",
        paste0(sir_fmt_int(nrow(background)), " rows from the ", input$background_source, " source."),
        "info"
      ))
    }
    if (sir_nonempty_table(rv$background)) {
      return(sir_panel_note(
        "Universe out of date",
        "The source or upstream resource changed. Rebuild before enrichment.",
        "warning"
      ))
    }
    sir_panel_note("Not built", "Enrichment remains disabled until an explicit universe is ready.", "info")
  })

  feature_preview_data <- reactive({
    columns <- if (sir_nonempty_table(rv$protein_features)) {
      intersect(
        c(
          "database", "feature_id", "name", "alt_name",
          "ensembl_transcript_id", "ensembl_peptide_id",
          "chr", "strand", "start", "stop", "method"
        ),
        names(rv$protein_features)
      )
    } else {
      NULL
    }
    sir_preview_data(
      rv$protein_features,
      preview_rows(),
      input$feature_search %||% "",
      columns = columns,
      search_columns = c(
        "database", "feature_id", "name", "alt_name",
        "ensembl_transcript_id", "ensembl_peptide_id"
      )
    )
  })

  output$feature_table <- DT::renderDT({
    sir_dt_widget(feature_preview_data(), page_length = preview_rows())
  })

  output$feature_table_caption <- renderUI({
    preview <- feature_preview_data()
    if (!preview$total) "No feature set loaded." else sir_preview_caption(preview)
  })

  reference_preview_data <- reactive({
    annotations <- if (!is.null(rv$reference)) rv$reference$annotations else NULL
    columns <- if (!is.null(annotations)) {
      intersect(
        c(
          "type", "gene_id", "gene_name", "transcript_id", "protein_id",
          "chr", "start", "end", "strand", "exon_id", "exon_number"
        ),
        names(annotations)
      )
    } else {
      NULL
    }
    sir_preview_data(
      annotations,
      preview_rows(),
      input$reference_search %||% "",
      columns = columns,
      search_columns = c("gene_id", "gene_name", "transcript_id", "protein_id", "chr")
    )
  })

  output$reference_table <- DT::renderDT({
    sir_dt_widget(reference_preview_data(), page_length = preview_rows())
  })

  output$reference_table_caption <- renderUI({
    preview <- reference_preview_data()
    if (!preview$total) "No reference loaded." else sir_preview_caption(preview)
  })

  resource_preview_data <- reactive({
    object <- if (sir_nonempty_table(rv$background)) {
      rv$background
    } else if (sir_nonempty_table(rv$ppi)) {
      rv$ppi
    } else {
      NULL
    }
    sir_preview_data(
      object,
      preview_rows(),
      input$resource_search %||% ""
    )
  })

  output$resource_table <- DT::renderDT({
    sir_dt_widget(resource_preview_data(), page_length = preview_rows())
  })

  output$resource_table_caption <- renderUI({
    preview <- resource_preview_data()
    if (!preview$total) {
      return("Load a PPI network or build an enrichment background to preview it.")
    }
    source <- if (sir_nonempty_table(rv$background)) "Enrichment background" else "PPI network"
    paste0(source, " · ", sir_preview_caption(preview))
  })

  domain_long <- reactive({
    sir_build_domain_long(current_domains())
  })

  domain_preview_data <- reactive({
    sir_preview_data(
      domain_long(),
      preview_rows(),
      input$domain_search %||% "",
      search_columns = c(
        "event_id", "event_type", "gene_id", "transcript_id_control",
        "transcript_id_case", "database", "domain", "direction"
      )
    )
  })

  output$domain_plot <- renderPlot({
    long <- domain_long()
    if (!nrow(long)) return(NULL)
    aggregate <- long[, .N, by = .(event_type, direction, database)]
    ggplot2::ggplot(
      aggregate,
      ggplot2::aes(x = event_type, y = N, fill = database)
    ) +
      ggplot2::geom_col() +
      ggplot2::facet_wrap(~direction, ncol = 1, scales = "free_y") +
      ggplot2::labs(x = "Event type", y = "Domain calls", fill = "Database") +
      ggplot2::theme_minimal(base_size = 12) +
      ggplot2::theme(legend.position = "bottom")
  })

  output$domain_note <- renderUI({
    long <- domain_long()
    if (!nrow(long)) {
      return(sir_panel_note(
        "No domain changes to show",
        "Load protein features and refresh consequences. A valid run may still yield zero changes.",
        "info"
      ))
    }
    sir_panel_note(
      "Domain calls",
      paste0(
        sir_fmt_int(nrow(long)), " calls across ",
        sir_fmt_int(data.table::uniqueN(long$event_id)),
        " events. Shared calls are shown separately from directional gains or losses."
      ),
      "info"
    )
  })

  output$domain_table <- DT::renderDT({
    sir_dt_widget(domain_preview_data(), page_length = preview_rows())
  })

  output$domain_table_caption <- renderUI({
    preview <- domain_preview_data()
    if (!preview$total) "No current domain-change rows." else sir_preview_caption(preview)
  })

  ppi_results <- reactive({
    final <- current_final()
    if (is.null(final) || !"n_ppi" %in% names(final)) {
      return(data.table::data.table())
    }
    tryCatch(
      data.table::as.data.table(SpliceImpactR::get_hits_ppi(final)),
      error = function(error) data.table::data.table()
    )
  })

  ppi_preview_data <- reactive({
    sir_preview_data(
      ppi_results(),
      preview_rows(),
      input$ppi_search %||% "",
      search_columns = c(
        "event_id", "event_type", "gene_id", "transcript_id_control",
        "transcript_id_case", "case_ppi", "control_ppi"
      )
    )
  })

  output$ppi_plot <- renderPlot({
    final <- current_final()
    if (is.null(final) || !"n_ppi" %in% names(final)) return(NULL)
    SpliceImpactR::plot_ppi_summary(final)
  })

  output$ppi_note <- renderUI({
    final <- current_final()
    if (is.null(final) || !"n_ppi" %in% names(final)) {
      return(sir_panel_note(
        "No current PPI layer",
        "Load a network and refresh consequences.",
        "info"
      ))
    }
    changed <- sum(final$n_ppi > 0, na.rm = TRUE)
    sir_panel_note(
      "Predicted interaction changes",
      paste0(
        sir_fmt_int(changed), " of ", sir_fmt_int(nrow(final)),
        " transcript pairs have at least one changed PPI partner."
      ),
      "info"
    )
  })

  output$ppi_table <- DT::renderDT({
    sir_dt_widget(ppi_preview_data(), page_length = preview_rows())
  })

  output$ppi_table_caption <- renderUI({
    preview <- ppi_preview_data()
    if (!preview$total) "No current PPI-change rows." else sir_preview_caption(preview)
  })

  output$integrated_plot <- renderPlot({
    integrated <- current_integrated()
    if (is.null(integrated)) return(NULL)
    integrated$plot
  })

  output$integrated_summary <- renderUI({
    integrated <- current_integrated()
    if (is.null(integrated)) {
      return(sir_panel_note(
        "Complete functional layers required",
        "The integrated summary becomes available when sequence, domain, and PPI results are all current.",
        "info"
      ))
    }
    by_type <- data.table::as.data.table(integrated$summaries$by_type)
    tagList(
      sir_panel_note(
        "Integrated view current",
        paste0(
          "The summary combines ", sir_fmt_int(nrow(current_final())),
          " transcript pairs across ", sir_fmt_int(nrow(by_type)), " event classes."
        ),
        "info"
      ),
      tags$table(
        class = "table table-sm sir-summary-table",
        tags$thead(tags$tr(lapply(names(by_type), tags$th))),
        tags$tbody(lapply(seq_len(nrow(by_type)), function(index) {
          tags$tr(lapply(by_type[index], function(value) tags$td(as.character(value))))
        }))
      )
    )
  })

  output$probe_event_type_ui <- renderUI({
    choices <- probe_event_types()
    selectInput(
      "probe_event_type",
      "Event type",
      choices = choices,
      selected = if (length(choices)) choices[[1]] else character()
    )
  })

  output$probe_event_ui <- renderUI({
    choices <- probe_event_ids()
    selectizeInput(
      "probe_event",
      "Event",
      choices = choices,
      selected = if (length(choices)) choices[[1]] else character(),
      options = list(placeholder = "Choose an event")
    )
  })

  output$overview_event_types_note <- renderUI({
    choices <- overview_event_types()
    if (!length(choices)) {
      return(sir_help_text(
        "Load sample-level events to populate the event-class selector."
      ))
    }
    sir_help_text(
      paste0(
        "Available in the loaded data: ",
        paste(choices, collapse = ", "),
        ". The selector controls the PSI panel for any listed class; ",
        "the HIT panel remains the terminal-exon HIT-index comparison."
      )
    )
  })

  output$probe_plot <- renderPlot({
    if (is.null(rv$probe)) return(NULL)
    rv$probe$plot
  })

  output$probe_note <- renderUI({
    if (is.null(rv$probe)) return(NULL)
    data <- data.table::as.data.table(rv$probe$data)
    sir_panel_note(
      "Selected event",
      paste0(
        sir_fmt_int(nrow(data)), " sample/form rows are shown for ",
        input$probe_event, "."
      ),
      "info"
    )
  })

  output$browser_plot <- renderPlot({
    if (is.null(rv$browser)) return(NULL)
    rv$browser$plot
  })

  output$browser_note <- renderUI({
    if (is.null(rv$browser)) return(NULL)
    sequence <- rv$browser$sequence
    classification <- if ("summary_classification" %in% names(sequence)) {
      paste(unique(na.omit(sequence$summary_classification)), collapse = ", ")
    } else {
      "not available"
    }
    sir_panel_note(
      "Direct transcript comparison",
      paste0(
        paste(rv$browser$transcripts, collapse = " versus "),
        " · classification: ", classification, "."
      ),
      "info"
    )
  })

  output$hit_plot <- renderPlot({
    if (is.null(rv$hit_compare)) return(NULL)
    rv$hit_compare$plot
  })

  output$psi_plot <- renderPlot({
    if (is.null(rv$psi_overview)) return(NULL)
    rv$psi_overview
  })

  output$enrichment_context <- renderUI({
    background <- current_background()
    if (is.null(background)) {
      return(sir_panel_note(
        "Background required",
        "Build a current enrichment universe on Functional data.",
        "warning"
      ))
    }
    sir_panel_note(
      "Explicit universe",
      paste0(
        sir_fmt_int(data.table::uniqueN(background$gene_id)),
        " unique genes from the ", input$background_source,
        " background are current."
      ),
      "info"
    )
  })

  output$enrichment_plot <- renderPlot({
    enrichment <- current_enrichment()
    if (is.null(enrichment)) return(NULL)
    enrichment$plot
  })

  enrichment_results <- reactive({
    enrichment <- current_enrichment()
    if (is.null(enrichment)) {
      return(data.table::data.table())
    }
    data.table::as.data.table(enrichment$results_signif)
  })

  enrichment_preview_data <- reactive({
    sir_preview_data(
      enrichment_results(),
      preview_rows(),
      input$enrichment_search %||% "",
      search_columns = c("source", "description", "term_name", "name", "ID")
    )
  })

  output$enrichment_note <- renderUI({
    enrichment <- current_enrichment()
    if (is.null(enrichment)) {
      if (!is.null(rv$enrichment)) {
        return(sir_panel_note(
          "Enrichment out of date",
          "A foreground, background, collection, or cutoff changed. Run enrichment again.",
          "warning"
        ))
      }
      return(NULL)
    }
    rows <- enrichment_results()
    sir_panel_note(
      "Current enrichment",
      if (nrow(rows)) {
        paste0(sir_fmt_int(nrow(rows)), " gene sets pass the adjusted p-value cutoff.")
      } else {
        "The run completed, but no gene sets pass the adjusted p-value cutoff."
      },
      "info"
    )
  })

  output$enrichment_table <- DT::renderDT({
    sir_dt_widget(enrichment_preview_data(), page_length = preview_rows())
  })

  output$enrichment_table_caption <- renderUI({
    preview <- enrichment_preview_data()
    if (!preview$total) "No current significant gene sets." else sir_preview_caption(preview)
  })

  log_text <- reactive({
    if (!nrow(rv$logs)) {
      return("No actions have been recorded in this session.")
    }
    paste(
      sprintf("[%s] %-7s %s", rv$logs$time, rv$logs$level, rv$logs$message),
      collapse = "\n"
    )
  })

  output$activity_log <- renderText(log_text())

  provenance_json <- reactive({
    jsonlite::toJSON(
      sir_workspace_manifest(rv),
      auto_unbox = TRUE,
      pretty = TRUE,
      null = "null",
      na = "null",
      digits = NA
    )
  })

  output$provenance_preview <- renderText(provenance_json())

  safe_raw_export <- reactive({
    if (!sir_nonempty_table(rv$raw_events)) {
      return(NULL)
    }
    result <- data.table::copy(data.table::as.data.table(rv$raw_events))
    if ("source_file" %in% names(result)) {
      result[, source_file := basename(source_file)]
    }
    result
  })

  table_download <- function(getter, stem) {
    downloadHandler(
      filename = function() sir_download_name(stem, "csv"),
      content = function(file) {
        value <- getter()
        if (is.null(value) || !is.data.frame(value)) {
          stop("This result is not currently available.", call. = FALSE)
        }
        data.table::fwrite(data.table::as.data.table(value), file)
      }
    )
  }

  output$download_raw <- table_download(function() safe_raw_export(), "spliceimpactr_raw_events")
  output$download_di <- table_download(function() current_di(), "spliceimpactr_di")
  output$download_sig <- table_download(function() significant_di(), "spliceimpactr_significant_di")
  output$download_matched <- table_download(function() current_matched(), "spliceimpactr_matched_events")
  output$download_sequence <- table_download(function() current_sequence(), "spliceimpactr_sequence")
  output$download_domains <- table_download(function() domain_long(), "spliceimpactr_domain_changes")
  output$download_ppi_results <- table_download(function() ppi_results(), "spliceimpactr_ppi_changes")
  output$download_final <- table_download(function() current_final(), "spliceimpactr_final_hits")
  output$download_enrichment <- table_download(
    function() {
      enrichment <- current_enrichment()
      if (is.null(enrichment)) return(NULL)
      data.table::as.data.table(enrichment$results_combined)
    },
    "spliceimpactr_enrichment"
  )

  output$download_log <- downloadHandler(
    filename = function() sir_download_name("spliceimpactr_activity", "txt"),
    content = function(file) {
      writeLines(log_text(), con = file, useBytes = TRUE)
    }
  )

  output$download_provenance <- downloadHandler(
    filename = function() sir_download_name("spliceimpactr_provenance", "json"),
    content = function(file) {
      writeLines(provenance_json(), con = file, useBytes = TRUE)
    }
  )
}
