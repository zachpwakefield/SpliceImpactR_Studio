sir_normalize_di <- function(object) {
  if (is.null(object) || !is.data.frame(object)) {
    stop("Differential inclusion did not return a table.", call. = FALSE)
  }
  out <- data.table::copy(data.table::as.data.table(object))

  rename_first <- function(target, candidates) {
    if (target %in% names(out)) {
      return(invisible(NULL))
    }
    found <- intersect(candidates, names(out))
    if (length(found)) {
      data.table::setnames(out, found[[1]], target)
    }
    invisible(NULL)
  }
  rename_first("delta_psi", c("delta.psi", "delta_psi.x", "delta.psi.x"))
  rename_first("padj", c("padj.x", "FDR", "fdr"))
  rename_first("p.value", c("pvalue", "p_value"))

  required <- c("delta_psi", "padj", "gene_id", "chr", "strand", "inc", "exc")
  missing <- setdiff(required, names(out))
  if (length(missing)) {
    stop(
      "The differential inclusion table is missing required columns: ",
      paste(missing, collapse = ", "), ".",
      call. = FALSE
    )
  }
  if (!"event_type" %in% names(out)) {
    out[, event_type := "SITE"]
  }
  if (!"form" %in% names(out)) {
    out[, form := "SITE"]
  }
  if (!"p.value" %in% names(out)) {
    out[, p.value := padj]
  }
  if (!"event_id" %in% names(out)) {
    if ("site_id" %in% names(out)) {
      out[, event_id := as.character(site_id)]
    } else {
      out[, event_id := paste(event_type, gene_id, chr, inc, exc, sep = "|")]
    }
  }

  out[, `:=`(
    event_id = as.character(event_id),
    event_type = as.character(event_type),
    gene_id = as.character(gene_id),
    chr = as.character(chr),
    strand = as.character(strand),
    inc = as.character(inc),
    exc = as.character(exc),
    form = as.character(form),
    delta_psi = suppressWarnings(as.numeric(delta_psi)),
    padj = suppressWarnings(as.numeric(padj)),
    p.value = suppressWarnings(as.numeric(p.value))
  )]
  if (!nrow(out)) {
    stop("The differential inclusion table is empty.", call. = FALSE)
  }
  out
}

sir_validate_raw_conditions <- function(raw_events) {
  dt <- data.table::as.data.table(raw_events)
  required <- c("condition", "sample")
  missing <- setdiff(required, names(dt))
  if (length(missing)) {
    stop(
      "The raw event table is missing required columns: ",
      paste(missing, collapse = ", "), ".",
      call. = FALSE
    )
  }
  conditions <- sort(unique(trimws(as.character(dt$condition))))
  conditions <- conditions[nzchar(conditions)]
  if (length(conditions) < 2L) {
    stop("Raw-event analysis requires at least two conditions.", call. = FALSE)
  }
  counts <- unique(dt[, .(condition = as.character(condition), sample = as.character(sample))])[
    , .(samples = uniqueN(sample)), by = condition
  ][order(condition)]
  counts
}

sir_resolve_condition_selection <- function(
  raw_events,
  reference_condition = NULL,
  case_condition = NULL
) {
  counts <- sir_validate_raw_conditions(raw_events)
  conditions <- counts$condition

  scalar_condition <- function(value) {
    value <- as.character(value %||% "")
    if (!length(value) || is.na(value[[1]])) {
      return("")
    }
    trimws(value[[1]])
  }

  reference_condition <- scalar_condition(reference_condition)
  case_condition <- scalar_condition(case_condition)

  if (!reference_condition %in% conditions) {
    reference_condition <- if ("control" %in% conditions) {
      "control"
    } else {
      conditions[[1]]
    }
  }

  case_choices <- setdiff(conditions, reference_condition)
  if (!case_condition %in% case_choices) {
    case_condition <- if ("case" %in% case_choices) {
      "case"
    } else {
      case_choices[[1]]
    }
  }

  list(
    reference = reference_condition,
    case = case_condition,
    counts = counts
  )
}

sir_prepare_conditions <- function(raw_events, sample_frame, reference_condition, case_condition) {
  reference_condition <- trimws(reference_condition %||% "")
  case_condition <- trimws(case_condition %||% "")
  if (!nzchar(reference_condition) || !nzchar(case_condition)) {
    stop("Choose both a reference condition and a case condition.", call. = FALSE)
  }
  if (identical(reference_condition, case_condition)) {
    stop("Reference and case conditions must be different.", call. = FALSE)
  }

  counts <- sir_validate_raw_conditions(raw_events)
  selected_counts <- counts[condition %in% c(reference_condition, case_condition)]
  if (nrow(selected_counts) != 2L) {
    stop("One or both selected conditions are not present in the raw event table.", call. = FALSE)
  }
  if (any(selected_counts$samples < 2L)) {
    detail <- paste0(selected_counts$condition, "=", selected_counts$samples, collapse = ", ")
    stop(
      "Each selected condition needs at least two distinct samples; found ",
      detail, ".",
      call. = FALSE
    )
  }

  raw <- data.table::copy(data.table::as.data.table(raw_events))[
    condition %in% c(reference_condition, case_condition)
  ]
  raw[, condition := data.table::fcase(
    condition == reference_condition, "control",
    condition == case_condition, "case",
    default = NA_character_
  )]

  manifest <- NULL
  if (!is.null(sample_frame) && is.data.frame(sample_frame)) {
    manifest <- as.data.frame(sample_frame, stringsAsFactors = FALSE)
    manifest <- manifest[manifest$condition %in% c(reference_condition, case_condition), , drop = FALSE]
    manifest$condition <- ifelse(
      manifest$condition == reference_condition,
      "control",
      "case"
    )
  }

  list(
    raw_events = raw,
    sample_frame = manifest,
    metadata = list(
      reference = reference_condition,
      case = case_condition,
      excluded = setdiff(counts$condition, c(reference_condition, case_condition)),
      sample_counts = stats::setNames(selected_counts$samples, selected_counts$condition)
    )
  )
}

sir_significant_di <- function(di, padj_threshold, dpsi_threshold) {
  if (!sir_nonempty_table(di)) {
    return(data.table::data.table())
  }
  data.table::as.data.table(SpliceImpactR::keep_sig_pairs(
    data.table::copy(data.table::as.data.table(di)),
    padj_thr = padj_threshold,
    dpsi_thr = dpsi_threshold
  ))
}

sir_feature_coverage <- function(features, reference) {
  if (!sir_nonempty_table(features) || is.null(reference)) {
    return(list(feature_transcripts = 0L, reference_transcripts = 0L, proportion = 0))
  }
  feature_ids <- unique(na.omit(as.character(features$ensembl_transcript_id)))
  annotation <- data.table::as.data.table(reference$annotations)
  reference_ids <- unique(na.omit(as.character(annotation$transcript_id)))
  list(
    feature_transcripts = length(intersect(feature_ids, reference_ids)),
    reference_transcripts = length(reference_ids),
    proportion = if (length(reference_ids)) {
      length(intersect(feature_ids, reference_ids)) / length(reference_ids)
    } else {
      0
    }
  )
}

sir_build_domain_long <- function(hits) {
  if (!sir_nonempty_table(hits)) {
    return(data.table::data.table())
  }
  dt <- tryCatch(
    data.table::as.data.table(SpliceImpactR::get_hits_domain(hits)),
    error = function(error) data.table::data.table()
  )
  if (!nrow(dt)) {
    return(dt)
  }

  list_columns <- c(
    case_only = "case_only_domains_list",
    control_only = "control_only_domains_list",
    shared = "either_domains_list"
  )
  pieces <- lapply(names(list_columns), function(direction) {
    column <- list_columns[[direction]]
    if (!column %in% names(dt)) {
      return(data.table::data.table())
    }
    base_columns <- intersect(
      c(
        "event_id", "event_type", "gene_id",
        "transcript_id_control", "transcript_id_case"
      ),
      names(dt)
    )
    rows <- lapply(seq_len(nrow(dt)), function(index) {
      domains <- as.character(unlist(dt[[column]][[index]], use.names = FALSE))
      domains <- domains[!is.na(domains) & nzchar(domains)]
      if (!length(domains)) {
        return(NULL)
      }
      base <- dt[index, ..base_columns]
      base[rep(1L, length(domains))][
        , `:=`(domain_id = domains, direction = direction)
      ]
    })
    data.table::rbindlist(rows, use.names = TRUE, fill = TRUE)
  })
  long <- data.table::rbindlist(pieces, use.names = TRUE, fill = TRUE)
  if (!nrow(long)) {
    return(long)
  }
  parsed <- data.table::tstrsplit(long$domain_id, ";", fixed = TRUE)
  long[, `:=`(
    database = if (length(parsed)) parsed[[1]] else "unlabeled",
    domain = if (length(parsed) >= 2L) parsed[[2]] else domain_id
  )]
  long
}

sir_di_parameters <- function(input) {
  list(
    min_total_reads = as.integer(input$min_reads),
    minimum_proportion = as.numeric(input$min_prop),
    terminal_fill = "event_max",
    cooks_cutoff = trimws(input$cooks_cutoff %||% "Inf")
  )
}

sir_threshold_parameters <- function(input) {
  list(
    padj = as.numeric(input$padj_threshold),
    delta_psi = as.numeric(input$dpsi_threshold)
  )
}
