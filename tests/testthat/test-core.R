test_that("bounded previews never materialize more than requested", {
  input <- data.table::data.table(
    gene_id = paste0("ENSG", seq_len(1000)),
    label = paste0("value[", seq_len(1000))
  )

  preview <- sir_preview_data(
    input,
    n = 17,
    query = "[",
    search_columns = "label"
  )

  expect_equal(preview$total, 1000L)
  expect_equal(preview$matching, 1000L)
  expect_equal(nrow(preview$data), 17L)
})

test_that("event-type helpers expose and filter every loaded class safely", {
  events <- data.table::data.table(
    event_id = paste0("event_", seq_along(SIR_SUPPORTED_EVENTS)),
    event_type = SIR_SUPPORTED_EVENTS,
    payload = seq_along(SIR_SUPPORTED_EVENTS)
  )

  expect_equal(sir_present_event_types(events), SIR_SUPPORTED_EVENTS)
  expect_equal(unname(sir_event_type_choices(SIR_SUPPORTED_EVENTS)), SIR_SUPPORTED_EVENTS)

  filtered <- sir_filter_event_type(events, "A3SS")
  expect_equal(nrow(filtered), 1L)
  expect_equal(filtered$event_type, "A3SS")
  expect_equal(filtered$event_id, events[event_type == "A3SS", event_id])
})

test_that("terminal-exon filling is fixed to event maximum", {
  parameters <- sir_di_parameters(list(
    min_reads = 10,
    min_prop = 0.5,
    terminal_fill = "none",
    cooks_cutoff = "Inf"
  ))

  expect_identical(parameters$terminal_fill, "event_max")
})

test_that("archive path policy accepts relative paths and rejects traversal", {
  expect_true(sir_archive_paths_are_safe(c(
    "manifest.csv",
    "case/sample_1/SE.MATS.JCEC.txt",
    "control/sample_2/"
  )))
  expect_false(sir_archive_paths_are_safe("../outside.txt"))
  expect_false(sir_archive_paths_are_safe("sample/../../outside.txt"))
  expect_false(sir_archive_paths_are_safe("/absolute/path.txt"))
  expect_false(sir_archive_paths_are_safe("C:/absolute/path.txt"))
  expect_false(sir_archive_paths_are_safe("sample//file.txt"))
})

test_that("manifest paths stay inside an extracted project", {
  root <- tempfile("sir-project-")
  dir.create(file.path(root, "samples", "case_1"), recursive = TRUE)
  manifest_path <- file.path(root, "manifest.csv")
  manifest <- data.frame(
    path = "samples/case_1",
    sample_name = "case_1",
    condition = "case"
  )

  resolved <- sir_resolve_manifest_paths(manifest, manifest_path, root)

  expect_true(dir.exists(resolved$path[[1]]))
  expect_true(sir_path_within(resolved$path[[1]], root))
})

test_that("condition preparation is explicit and excludes unselected levels", {
  raw <- data.table::data.table(
    event_id = rep("event_1", 8),
    sample = rep(c("A1", "A2", "B1", "B2"), each = 2),
    condition = rep(c("baseline", "baseline", "treated", "treated"), each = 2),
    form = rep(c("INC", "EXC"), 4)
  )
  extra <- data.table::copy(raw)
  extra[, `:=`(
    sample = sub("^[AB]", "C", sample),
    condition = "other"
  )]
  raw <- data.table::rbindlist(list(raw, extra))

  prepared <- sir_prepare_conditions(
    raw,
    sample_frame = NULL,
    reference_condition = "baseline",
    case_condition = "treated"
  )

  expect_setequal(unique(prepared$raw_events$condition), c("control", "case"))
  expect_equal(prepared$metadata$reference, "baseline")
  expect_equal(prepared$metadata$case, "treated")
  expect_equal(prepared$metadata$excluded, "other")
})

test_that("condition defaults prefer control as reference and case as comparison", {
  raw <- data.table::data.table(
    event_id = rep("event_1", 8),
    sample = rep(c("case_1", "case_2", "control_1", "control_2"), each = 2),
    condition = rep(c("case", "case", "control", "control"), each = 2),
    form = rep(c("INC", "EXC"), 4)
  )

  selection <- sir_resolve_condition_selection(raw)

  expect_identical(selection$reference, "control")
  expect_identical(selection$case, "case")
})

test_that("condition preparation requires distinct replicated groups", {
  raw <- data.table::data.table(
    event_id = rep("event_1", 4),
    sample = rep(c("A1", "B1"), each = 2),
    condition = rep(c("baseline", "treated"), each = 2),
    form = rep(c("INC", "EXC"), 2)
  )

  expect_error(
    sir_prepare_conditions(raw, NULL, "baseline", "treated"),
    "at least two distinct samples"
  )
  expect_error(
    sir_prepare_conditions(raw, NULL, "baseline", "baseline"),
    "must be different"
  )
})

test_that("generic DI normalization creates stable mapping columns", {
  input <- data.table::data.table(
    gene_id = "ENSG1",
    chr = "chr1",
    strand = "+",
    inc = "10-20",
    exc = "10-15",
    delta.psi = 0.3,
    FDR = 0.01
  )

  normalized <- sir_normalize_di(input)

  expect_true(all(c("event_id", "event_type", "form", "delta_psi", "padj") %in% names(normalized)))
  expect_equal(normalized$event_type, "SITE")
  expect_equal(normalized$form, "SITE")
  expect_equal(normalized$p.value, normalized$padj)
})

test_that("state invalidation follows the declared dependency graph", {
  state <- new.env(parent = emptyenv())
  for (field in SIR_STATE_FIELDS) {
    state[[field]] <- data.frame(value = 1)
  }
  state$provenance <- stats::setNames(
    lapply(SIR_STATE_FIELDS, function(field) list(signature = field)),
    SIR_STATE_FIELDS
  )

  sir_invalidate_state(state, "thresholds")

  expect_null(state$matched)
  expect_null(state$sequence_results)
  expect_null(state$final_results)
  expect_null(state$enrichment)
  expect_true(is.data.frame(state$di))
  expect_true(is.data.frame(state$protein_features))
})

test_that("deployed reference catalog never advertises absent bundles", {
  catalog <- sir_reference_catalog()
  expect_true("demo_human_v45" %in% names(catalog))
  expect_true(all(vapply(catalog, `[[`, logical(1), "available")))
  expect_true(all(vapply(catalog, function(info) {
    is.null(info$paths) || all(file.exists(info$paths))
  }, logical(1))))
})

test_that("server reference discovery requires a complete trusted bundle", {
  root <- tempfile("sir-reference-root-")
  dir.create(root, recursive = TRUE)

  annotations <- data.table::data.table(
    type = "transcript",
    gene_id = "ENSG_TEST",
    transcript_id = "ENST_TEST"
  )
  sequences <- data.table::data.table(
    transcript_id = "ENST_TEST",
    transcript_seq = "ATG",
    protein_seq = "M"
  )
  hybrids <- list(
    first_hybrids = data.table::data.table(),
    last_hybrids = data.table::data.table()
  )

  saveRDS(annotations, file.path(root, "human_gencode_v99.gtf.rds"))
  saveRDS(sequences, file.path(root, "human_gencode_v99_sequences.rds"))
  expect_length(sir_discover_preprocessed_references(root), 0L)

  saveRDS(hybrids, file.path(root, "human_gencode_v99_hybrids.rds"))
  discovered <- sir_discover_preprocessed_references(root)
  expect_length(discovered, 1L)

  loaded <- sir_load_reference(discovered[[1]])
  expect_equal(loaded$annotations$transcript_id, "ENST_TEST")
  expect_equal(loaded$sequences$protein_seq, "M")
})

test_that("preprocessed references take precedence over duplicate raw assets", {
  root <- tempfile("sir-reference-priority-")
  dir.create(root, recursive = TRUE)

  saveRDS(
    data.table::data.table(
      type = "transcript",
      gene_id = "ENSG_TEST",
      transcript_id = "ENST_TEST"
    ),
    file.path(root, "human_gencode_v77.gtf.rds")
  )
  saveRDS(
    data.table::data.table(
      transcript_id = "ENST_TEST",
      transcript_seq = "ATG",
      protein_seq = "M"
    ),
    file.path(root, "human_gencode_v77_sequences.rds")
  )
  saveRDS(
    list(
      first_hybrids = data.table::data.table(),
      last_hybrids = data.table::data.table()
    ),
    file.path(root, "human_gencode_v77_hybrids.rds")
  )
  file.create(file.path(root, c(
    "gencode.v77.annotation.gtf.gz",
    "gencode.v77.pc_transcripts.fa.gz",
    "gencode.v77.pc_translations.fa.gz"
  )))

  catalog <- sir_reference_catalog(roots = root)
  full <- catalog[setdiff(names(catalog), "demo_human_v45")]

  expect_length(full, 1L)
  expect_identical(full[[1]]$kind, "preprocessed")
})

test_that("hosted reference choices expose bundled and browser-safe sources only", {
  catalog <- list(
    demo_human_v45 = list(
      id = "demo_human_v45",
      kind = "demo",
      species = "human",
      gencode_release = 45L
    ),
    bundled_human_v45 = list(
      id = "bundled_human_v45",
      label = "Internal full reference",
      short_label = "Human v45",
      kind = "preprocessed",
      species = "human",
      species_label = "Human",
      gencode_release = 45L,
      ensembl_release = 111L,
      assembly = "Server reference library",
      paths = c(
        annotations = "/srv/human_gencode_v45.gtf.rds",
        sequences = "/srv/human_gencode_v45_sequences.rds",
        hybrids = "/srv/human_gencode_v45_hybrids.rds"
      ),
      source_root = "/srv",
      note = "Internal storage detail"
    )
  )

  choices <- sir_reference_mode_choices(hosted = TRUE, catalog = catalog)
  expect_identical(unname(choices), c("test", "bundled", "upload"))
  expect_false(any(c("library", "cached", "download") %in% unname(choices)))
  expect_identical(
    sir_default_reference_mode(hosted = TRUE, catalog = catalog),
    "bundled"
  )

  bundled <- sir_bundled_human_reference(catalog)
  expect_identical(bundled$kind, "preprocessed")
  expect_match(bundled$label, "Bundled human GENCODE 45")
  expect_false(grepl("/srv", bundled$note, fixed = TRUE))
})

test_that("GENCODE browser links and reference uploads use raw safe formats", {
  urls <- sir_gencode_urls("human", 45)
  expect_match(urls[["gtf"]], "gencode.v45.annotation.gtf.gz", fixed = TRUE)
  expect_match(
    urls[["transcripts"]],
    "gencode.v45.pc_transcripts.fa.gz",
    fixed = TRUE
  )
  expect_match(
    urls[["translations"]],
    "gencode.v45.pc_translations.fa.gz",
    fixed = TRUE
  )

  source <- tempfile(fileext = ".gtf.gz")
  writeBin(charToRaw("reference"), source)
  destination <- tempfile("reference-upload-")
  fileinfo <- data.frame(
    name = "gencode.v45.annotation.gtf.gz",
    size = file.info(source)$size,
    type = "application/gzip",
    datapath = source,
    stringsAsFactors = FALSE
  )

  copied <- sir_copy_reference_upload(fileinfo, destination, "gtf")
  expect_true(file.exists(copied))
  expect_error(
    sir_copy_reference_upload(
      transform(fileinfo, name = "unsafe-reference.rds"),
      destination,
      "gtf"
    ),
    "\\.gtf\\.gz"
  )
})

test_that("managed-cache requests validate release, TSL, and download policy", {
  cached <- sir_reference_request(
    "cached",
    species = "mouse",
    release = 36,
    filter_tsl = c("1", "3")
  )

  expect_identical(cached$kind, "cached")
  expect_identical(cached$species, "mouse")
  expect_identical(cached$gencode_release, 36L)
  expect_identical(cached$filter_tsl, c("1", "3"))
  expect_error(
    sir_reference_request("cached", release = 0),
    "positive integer"
  )
  expect_error(
    sir_reference_request("cached", filter_tsl = character()),
    "at least one transcript-support level"
  )
  expect_error(
    sir_reference_request("download", allow_download = FALSE),
    "downloads are disabled"
  )
})

test_that("hosted-runtime policy supports an explicit deployment override", {
  previous <- Sys.getenv("SPLICEIMPACTR_HOSTED", unset = NA_character_)
  on.exit({
    if (is.na(previous)) {
      Sys.unsetenv("SPLICEIMPACTR_HOSTED")
    } else {
      Sys.setenv(SPLICEIMPACTR_HOSTED = previous)
    }
  }, add = TRUE)

  Sys.setenv(SPLICEIMPACTR_HOSTED = "true")
  expect_true(sir_is_hosted_runtime())
  Sys.setenv(SPLICEIMPACTR_HOSTED = "false")
  expect_false(sir_is_hosted_runtime())
})

test_that("PPI uploads require switch evidence rather than bare edges", {
  expect_error(
    sir_validate_ppi_table(data.table::data.table(geneA = "A", geneB = "B")),
    "DDI/DDI_A/DDI_B"
  )
  expect_silent(sir_validate_ppi_table(data.table::data.table(
    geneA = "A",
    geneB = "B",
    DDI = TRUE,
    DDI_A = "PF0001",
    DDI_B = "PF0002"
  )))
})
