sir_env_flag <- function(name, default = FALSE) {
  value <- tolower(trimws(Sys.getenv(name, unset = "")))
  if (!nzchar(value)) {
    return(isTRUE(default))
  }
  if (value %in% c("1", "true", "yes", "on")) {
    return(TRUE)
  }
  if (value %in% c("0", "false", "no", "off")) {
    return(FALSE)
  }
  warning(
    name, " must be true/false, yes/no, on/off, or 1/0; using the default.",
    call. = FALSE
  )
  isTRUE(default)
}

sir_config_path <- function(path, base_dir = SIR_APP_DIR) {
  path <- path.expand(trimws(path %||% ""))
  if (!nzchar(path)) {
    return(NULL)
  }
  if (!grepl("^(/|[A-Za-z]:[/\\\\])", path)) {
    path <- file.path(base_dir, path)
  }
  normalizePath(path, winslash = "/", mustWork = FALSE)
}

sir_env_paths <- function(name, base_dir = SIR_APP_DIR) {
  raw <- trimws(Sys.getenv(name, unset = ""))
  if (!nzchar(raw)) {
    return(character())
  }
  values <- trimws(strsplit(raw, .Platform$path.sep, fixed = TRUE)[[1]])
  values <- values[nzchar(values)]
  unique(vapply(values, sir_config_path, character(1), base_dir = base_dir))
}

sir_reference_roots <- function(app_dir = SIR_APP_DIR) {
  configured <- sir_env_paths("SPLICEIMPACTR_REFERENCE_DIRS", app_dir)
  defaults <- c(
    file.path(app_dir, "references"),
    file.path(app_dir, "annotation_cache"),
    app_dir
  )
  candidates <- unique(c(configured, defaults))
  candidates <- vapply(candidates, sir_config_path, character(1), base_dir = app_dir)
  candidates[dir.exists(candidates)]
}

sir_default_cache_dir <- function(app_dir = SIR_APP_DIR) {
  configured <- sir_config_path(
    Sys.getenv("SPLICEIMPACTR_CACHE_DIR", unset = ""),
    base_dir = app_dir
  )
  if (!is.null(configured)) {
    return(configured)
  }
  if (sir_is_hosted_runtime()) {
    return(normalizePath(
      tools::R_user_dir("SpliceImpactRStudio", which = "cache"),
      winslash = "/",
      mustWork = FALSE
    ))
  }
  normalizePath(
    file.path(app_dir, "runtime_cache"),
    winslash = "/",
    mustWork = FALSE
  )
}

sir_default_session_root <- function() {
  configured <- sir_config_path(
    Sys.getenv("SPLICEIMPACTR_SESSION_DIR", unset = ""),
    base_dir = SIR_APP_DIR
  )
  if (!is.null(configured)) {
    return(configured)
  }
  normalizePath(
    file.path(tempdir(), "spliceimpactr-studio-v4", "sessions"),
    winslash = "/",
    mustWork = FALSE
  )
}

SIR_CACHE_DIR <- sir_default_cache_dir()
SIR_ANNOTATION_CACHE_DIR <- file.path(SIR_CACHE_DIR, "annotations")
SIR_FEATURE_CACHE_DIR <- file.path(SIR_CACHE_DIR, "features")
SIR_SESSION_ROOT <- sir_default_session_root()
SIR_ALLOW_REFERENCE_DOWNLOADS <- sir_env_flag(
  "SPLICEIMPACTR_ALLOW_REFERENCE_DOWNLOADS",
  default = !sir_is_hosted_runtime()
)

sir_ensure_directory <- function(path, label = "Directory") {
  if (!dir.exists(path)) {
    created <- dir.create(path, recursive = TRUE, showWarnings = FALSE)
    if (!isTRUE(created) && !dir.exists(path)) {
      stop(label, " could not be created: ", path, call. = FALSE)
    }
  }
  probe <- tempfile(".spliceimpactr-write-test-", tmpdir = path)
  writable <- isTRUE(file.create(probe, showWarnings = FALSE))
  if (writable) {
    unlink(probe, force = TRUE)
  }
  if (!writable) {
    stop(
      label, " is not writable: ", path,
      ". Set SPLICEIMPACTR_CACHE_DIR to a writable persistent location.",
      call. = FALSE
    )
  }
  normalizePath(path, winslash = "/", mustWork = TRUE)
}

sir_species_label <- function(species) {
  switch(species, human = "Human", mouse = "Mouse", species)
}

sir_release_tag <- function(species, release) {
  if (identical(species, "mouse")) {
    paste0("M", as.integer(release))
  } else {
    as.character(as.integer(release))
  }
}

sir_gencode_ensembl_release <- function(species, release) {
  release <- as.integer(release)
  known <- if (identical(species, "human")) {
    c(`45` = 111L)
  } else {
    c(`36` = 111L)
  }
  value <- unname(known[as.character(release)])
  if (!length(value) || is.na(value)) NA_integer_ else as.integer(value)
}

sir_reference_id <- function(kind, paths, species, release) {
  identity <- paste(
    kind,
    species,
    as.integer(release),
    paste(normalizePath(paths, winslash = "/", mustWork = FALSE), collapse = "|"),
    sep = "::"
  )
  paste0(
    "library_",
    species,
    "_v",
    sir_release_tag(species, release),
    "_",
    substr(digest::digest(identity, algo = "xxhash64", serialize = FALSE), 1L, 10L)
  )
}

sir_discover_preprocessed_references <- function(roots = sir_reference_roots()) {
  found <- list()
  for (root in roots) {
    annotation_files <- sort(list.files(
      root,
      pattern = "^(human|mouse)_gencode_v(M?[0-9]+)\\.gtf\\.rds$",
      full.names = TRUE,
      ignore.case = TRUE
    ))
    for (annotation_path in annotation_files) {
      filename <- basename(annotation_path)
      match <- regexec(
        "^(human|mouse)_gencode_v(M?)([0-9]+)\\.gtf\\.rds$",
        filename,
        ignore.case = TRUE
      )
      parts <- regmatches(filename, match)[[1]]
      if (length(parts) != 4L) {
        next
      }
      species <- tolower(parts[[2]])
      release <- as.integer(parts[[4]])
      stem <- sub("\\.gtf\\.rds$", "", filename, ignore.case = TRUE)
      paths <- c(
        annotations = annotation_path,
        sequences = file.path(root, paste0(stem, "_sequences.rds")),
        hybrids = file.path(root, paste0(stem, "_hybrids.rds"))
      )
      if (!all(file.exists(paths))) {
        next
      }
      id <- sir_reference_id("preprocessed", paths, species, release)
      found[[id]] <- list(
        id = id,
        label = paste0(
          sir_species_label(species), " full reference · GENCODE ",
          sir_release_tag(species, release), " · preprocessed"
        ),
        short_label = paste0(
          sir_species_label(species), " GENCODE ",
          sir_release_tag(species, release)
        ),
        kind = "preprocessed",
        species = species,
        species_label = sir_species_label(species),
        gencode_release = release,
        ensembl_release = sir_gencode_ensembl_release(species, release),
        assembly = "Server reference library",
        available = TRUE,
        paths = paths,
        source_root = root,
        note = paste0(
          "Trusted preprocessed RDS bundle discovered in ",
          basename(root), ". It is read in place and never rewritten."
        )
      )
    }
  }
  found
}

sir_discover_raw_references <- function(roots = sir_reference_roots()) {
  found <- list()
  for (root in roots) {
    gtf_files <- sort(list.files(
      root,
      pattern = "^gencode\\.v(M?[0-9]+)\\.annotation\\.gtf(?:\\.gz)?$",
      full.names = TRUE,
      ignore.case = TRUE
    ))
    for (gtf_path in gtf_files) {
      filename <- basename(gtf_path)
      match <- regexec(
        "^gencode\\.v(M?)([0-9]+)\\.annotation\\.gtf(?:\\.gz)?$",
        filename,
        ignore.case = TRUE
      )
      parts <- regmatches(filename, match)[[1]]
      if (length(parts) != 3L) {
        next
      }
      species <- if (nzchar(parts[[2]])) "mouse" else "human"
      release <- as.integer(parts[[3]])
      tag <- paste0("v", parts[[2]], parts[[3]])
      paths <- c(
        gtf = gtf_path,
        transcripts = file.path(root, paste0("gencode.", tag, ".pc_transcripts.fa.gz")),
        translations = file.path(root, paste0("gencode.", tag, ".pc_translations.fa.gz"))
      )
      if (!all(file.exists(paths))) {
        next
      }
      id <- sir_reference_id("raw", paths, species, release)
      found[[id]] <- list(
        id = id,
        label = paste0(
          sir_species_label(species), " full reference · GENCODE ",
          sir_release_tag(species, release), " · local GTF/FASTA"
        ),
        short_label = paste0(
          sir_species_label(species), " GENCODE ",
          sir_release_tag(species, release)
        ),
        kind = "raw",
        species = species,
        species_label = sir_species_label(species),
        gencode_release = release,
        ensembl_release = sir_gencode_ensembl_release(species, release),
        assembly = "Server reference library",
        available = TRUE,
        paths = paths,
        source_root = root,
        note = paste0(
          "Trusted local GTF and FASTA files discovered in ", basename(root),
          ". Their processed objects are written only to the managed cache."
        )
      )
    }
  }
  found
}

sir_reference_catalog <- function(
  app_dir = SIR_APP_DIR,
  roots = sir_reference_roots(app_dir)
) {
  catalog <- list(
    demo_human_v45 = list(
      id = "demo_human_v45",
      label = "Human demo subset · GENCODE 45",
      short_label = "Demo subset",
      kind = "demo",
      species = "human",
      species_label = "Human",
      gencode_release = 45L,
      ensembl_release = 111L,
      assembly = "GRCh38.p14",
      available = TRUE,
      paths = NULL,
      source_root = NULL,
      note = "Small packaged subset for learning and validation; not a genome-wide reference."
    )
  )

  preprocessed <- sir_discover_preprocessed_references(roots)
  raw <- sir_discover_raw_references(roots)
  preprocessed_keys <- vapply(
    preprocessed,
    function(info) paste(info$species, info$gencode_release, sep = "::"),
    character(1)
  )
  if (length(raw) && length(preprocessed_keys)) {
    raw_keys <- vapply(
      raw,
      function(info) paste(info$species, info$gencode_release, sep = "::"),
      character(1)
    )
    raw <- raw[!raw_keys %in% preprocessed_keys]
  }

  entries <- c(preprocessed, raw)
  if (length(entries)) {
    entry_keys <- vapply(
      entries,
      function(info) paste(info$kind, info$species, info$gencode_release, sep = "::"),
      character(1)
    )
    entries <- entries[!duplicated(entry_keys)]
  }
  c(catalog, entries)
}

sir_reference_info <- function(reference_id, catalog = sir_reference_catalog()) {
  if (!length(reference_id) || !reference_id %in% names(catalog)) {
    stop("Choose an available server reference before continuing.", call. = FALSE)
  }
  catalog[[reference_id]]
}

sir_bundled_human_reference <- function(
  catalog = sir_reference_catalog(),
  release = 45L
) {
  matches <- Filter(
    function(info) {
      identical(info$kind, "preprocessed") &&
        identical(info$species, "human") &&
        identical(as.integer(info$gencode_release), as.integer(release))
    },
    catalog
  )
  if (!length(matches)) {
    stop(
      "The bundled human GENCODE ",
      sir_release_tag("human", release),
      " reference is not available in this deployment.",
      call. = FALSE
    )
  }
  info <- matches[[1]]
  info$label <- paste0(
    "Bundled human GENCODE ",
    sir_release_tag("human", release),
    " reference"
  )
  info$short_label <- paste0(
    "Human GENCODE ",
    sir_release_tag("human", release)
  )
  info$assembly <- "GRCh38.p14"
  info$note <- paste0(
    "Full preprocessed reference bundled with this application. ",
    "It is loaded read-only and reused by the current R worker."
  )
  info
}

sir_reference_mode_choices <- function(
  hosted = sir_is_hosted_runtime(),
  catalog = sir_reference_catalog()
) {
  if (isTRUE(hosted)) {
    choices <- c("Use bundled test set" = "test")
    bundled <- tryCatch(
      sir_bundled_human_reference(catalog),
      error = function(error) NULL
    )
    if (!is.null(bundled)) {
      choices <- c(
        choices,
        "Bundled human GENCODE 45 reference · recommended" = "bundled"
      )
    }
    return(c(
      choices,
      "Load GENCODE 45 files from this computer" = "upload"
    ))
  }

  choices <- c("Bundled test data" = "test")
  if (length(setdiff(names(catalog), "demo_human_v45"))) {
    choices <- c(choices, "Server reference library" = "library")
  }
  choices <- c(choices, "Managed processed cache" = "cached")
  if (isTRUE(SIR_ALLOW_REFERENCE_DOWNLOADS)) {
    choices <- c(choices, "Download from GENCODE to managed cache" = "download")
  }
  choices
}

sir_default_reference_mode <- function(
  hosted = sir_is_hosted_runtime(),
  catalog = sir_reference_catalog()
) {
  choices <- unname(sir_reference_mode_choices(hosted, catalog))
  if (isTRUE(hosted) && "bundled" %in% choices) "bundled" else "test"
}

sir_gencode_urls <- function(species = "human", release = 45L) {
  species <- match.arg(species, c("human", "mouse"))
  release <- suppressWarnings(as.integer(release))
  if (!is.finite(release) || release < 1L) {
    stop("GENCODE release must be a positive integer.", call. = FALSE)
  }
  if (identical(species, "human")) {
    tag <- paste0("v", release)
    base <- paste0(
      "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_",
      release
    )
  } else {
    tag <- paste0("vM", release)
    base <- paste0(
      "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M",
      release
    )
  }
  c(
    gtf = paste0(base, "/gencode.", tag, ".annotation.gtf.gz"),
    transcripts = paste0(
      base,
      "/gencode.", tag, ".pc_transcripts.fa.gz"
    ),
    translations = paste0(
      base,
      "/gencode.", tag, ".pc_translations.fa.gz"
    )
  )
}

sir_reference_request <- function(
  mode,
  library_id = NULL,
  species = "human",
  release = 45L,
  filter_tsl = c("1", "2", "3"),
  allow_download = SIR_ALLOW_REFERENCE_DOWNLOADS
) {
  mode <- match.arg(mode, c("test", "library", "cached", "download"))
  if (identical(mode, "test")) {
    return(sir_reference_info("demo_human_v45"))
  }
  if (identical(mode, "library")) {
    info <- sir_reference_info(library_id)
    if (identical(info$kind, "demo")) {
      stop("Choose a full server-library reference.", call. = FALSE)
    }
    return(info)
  }

  species <- match.arg(species, c("human", "mouse"))
  release <- suppressWarnings(as.integer(release))
  if (!is.finite(release) || release < 1L) {
    stop("GENCODE release must be a positive integer.", call. = FALSE)
  }
  filter_tsl <- sort(unique(as.character(filter_tsl)))
  if (!length(filter_tsl) || !all(filter_tsl %in% as.character(1:5))) {
    stop("Choose at least one transcript-support level from 1 through 5.", call. = FALSE)
  }
  if (identical(mode, "download") && !isTRUE(allow_download)) {
    stop(
      "Reference downloads are disabled on this server. Configure a local reference root ",
      "or set SPLICEIMPACTR_ALLOW_REFERENCE_DOWNLOADS=true.",
      call. = FALSE
    )
  }

  tag <- sir_release_tag(species, release)
  list(
    id = paste(mode, species, tag, paste(filter_tsl, collapse = "-"), sep = "__"),
    label = paste0(
      sir_species_label(species), " GENCODE ", tag,
      if (identical(mode, "cached")) " · managed cache" else " · GENCODE download"
    ),
    short_label = paste0(sir_species_label(species), " GENCODE ", tag),
    kind = mode,
    species = species,
    species_label = sir_species_label(species),
    gencode_release = release,
    ensembl_release = sir_gencode_ensembl_release(species, release),
    assembly = "GENCODE reference",
    available = TRUE,
    paths = NULL,
    source_root = NULL,
    filter_tsl = filter_tsl,
    note = if (identical(mode, "cached")) {
      paste0(
        "Loads processed objects from the managed annotation cache at ",
        SIR_ANNOTATION_CACHE_DIR, "."
      )
    } else {
      paste0(
        "Downloads the three official GENCODE assets into ",
        SIR_ANNOTATION_CACHE_DIR,
        " and stores processed objects there for reuse."
      )
    }
  )
}

sir_reference_fingerprint <- function(info) {
  if (is.null(info$paths)) {
    return(info$id)
  }
  metadata <- file.info(info$paths)
  digest::digest(
    list(
      id = info$id,
      paths = normalizePath(info$paths, winslash = "/", mustWork = TRUE),
      size = metadata$size,
      modified = as.numeric(metadata$mtime)
    ),
    algo = "xxhash64"
  )
}

sir_validate_reference_bundle <- function(value) {
  required <- c("annotations", "sequences", "hybrids")
  if (!is.list(value) || !all(required %in% names(value))) {
    stop(
      "The reference bundle must contain annotations, sequences, and hybrids.",
      call. = FALSE
    )
  }
  if (!is.data.frame(value$annotations) || !nrow(value$annotations)) {
    stop("The reference annotation table is empty or invalid.", call. = FALSE)
  }
  annotation_columns <- c("gene_id", "transcript_id", "type")
  missing_annotations <- setdiff(annotation_columns, names(value$annotations))
  if (length(missing_annotations)) {
    stop(
      "The reference annotation table is missing: ",
      paste(missing_annotations, collapse = ", "), ".",
      call. = FALSE
    )
  }
  if (!is.data.frame(value$sequences) || !nrow(value$sequences)) {
    stop("The reference sequence table is empty or invalid.", call. = FALSE)
  }
  missing_sequences <- setdiff(
    c("transcript_id", "transcript_seq", "protein_seq"),
    names(value$sequences)
  )
  if (length(missing_sequences)) {
    stop(
      "The reference sequence table is missing: ",
      paste(missing_sequences, collapse = ", "), ".",
      call. = FALSE
    )
  }
  value
}

sir_load_reference <- function(request) {
  if (is.character(request) && length(request) == 1L) {
    request <- sir_reference_info(request)
  }
  if (!is.list(request) || is.null(request$id) || is.null(request$kind)) {
    stop("The reference request is invalid.", call. = FALSE)
  }

  cacheable <- !identical(request$kind, "uploaded_raw")
  cache_key <- paste0("reference__", sir_reference_fingerprint(request))
  if (
    isTRUE(cacheable) &&
    exists(cache_key, envir = .sir_process_cache, inherits = FALSE)
  ) {
    return(get(cache_key, envir = .sir_process_cache, inherits = FALSE))
  }

  value <- switch(
    request$kind,
    demo = SpliceImpactR::load_example_data("annotation_df")$annotation_df,
    preprocessed = list(
      annotations = readRDS(request$paths[["annotations"]]),
      sequences = readRDS(request$paths[["sequences"]]),
      hybrids = readRDS(request$paths[["hybrids"]])
    ),
    raw = {
      cache_dir <- sir_ensure_directory(
        SIR_ANNOTATION_CACHE_DIR,
        "Managed annotation cache"
      )
      SpliceImpactR::get_annotation(
        load = "path",
        base_dir = cache_dir,
        species = request$species,
        release = request$gencode_release,
        gtf_path = request$paths[["gtf"]],
        transcript_path = request$paths[["transcripts"]],
        translation_path = request$paths[["translations"]],
        filter_tsl = request$filter_tsl %||% c("1", "2", "3")
      )
    },
    uploaded_raw = {
      cache_dir <- sir_ensure_directory(
        request$workspace,
        "Uploaded reference workspace"
      )
      SpliceImpactR::get_annotation(
        load = "path",
        base_dir = cache_dir,
        species = request$species,
        release = request$gencode_release,
        gtf_path = request$paths[["gtf"]],
        transcript_path = request$paths[["transcripts"]],
        translation_path = request$paths[["translations"]],
        filter_tsl = request$filter_tsl %||% c("1", "2", "3")
      )
    },
    cached = {
      cache_dir <- sir_ensure_directory(
        SIR_ANNOTATION_CACHE_DIR,
        "Managed annotation cache"
      )
      SpliceImpactR::get_annotation(
        load = "cached",
        base_dir = cache_dir,
        species = request$species,
        release = request$gencode_release,
        filter_tsl = request$filter_tsl
      )
    },
    download = {
      if (!isTRUE(SIR_ALLOW_REFERENCE_DOWNLOADS)) {
        stop("Reference downloads are disabled on this server.", call. = FALSE)
      }
      cache_dir <- sir_ensure_directory(
        SIR_ANNOTATION_CACHE_DIR,
        "Managed annotation cache"
      )
      SpliceImpactR::get_annotation(
        load = "link",
        base_dir = cache_dir,
        species = request$species,
        release = request$gencode_release,
        filter_tsl = request$filter_tsl
      )
    },
    stop("Unknown reference kind: ", request$kind, call. = FALSE)
  )

  value <- sir_validate_reference_bundle(value)
  if (isTRUE(cacheable)) {
    assign(cache_key, value, envir = .sir_process_cache)
  }
  value
}

sir_reference_public_info <- function(info) {
  info$paths <- NULL
  info$workspace <- NULL
  if (!is.null(info$source_root)) {
    info$source_root <- basename(info$source_root)
  }
  info
}

sir_reference_storage_summary <- function() {
  list(
    roots = sir_reference_roots(),
    cache = SIR_ANNOTATION_CACHE_DIR,
    feature_cache = SIR_FEATURE_CACHE_DIR,
    downloads_enabled = isTRUE(SIR_ALLOW_REFERENCE_DOWNLOADS),
    hosted = sir_is_hosted_runtime()
  )
}
