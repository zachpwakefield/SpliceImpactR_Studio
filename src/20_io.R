sir_upload_suffix <- function(name) {
  lower <- tolower(basename(name %||% ""))
  known <- c(".csv.gz", ".tsv.gz", ".txt.gz", ".csv", ".tsv", ".txt", ".zip")
  hit <- known[endsWith(lower, known)]
  if (!length(hit)) "" else hit[[1]]
}

sir_copy_upload <- function(fileinfo, destination_dir) {
  if (is.null(fileinfo) || !nrow(fileinfo)) {
    stop("Choose a file to upload.", call. = FALSE)
  }
  if (nrow(fileinfo) != 1L) {
    stop("Upload one file at a time.", call. = FALSE)
  }
  if (is.finite(fileinfo$size[[1]]) && fileinfo$size[[1]] > SIR_MAX_UPLOAD_MB * 1024^2) {
    stop(
      "The upload is larger than the configured ",
      sir_fmt_num(SIR_MAX_UPLOAD_MB, 0), " MB limit.",
      call. = FALSE
    )
  }

  suffix <- sir_upload_suffix(fileinfo$name[[1]])
  if (!nzchar(suffix)) {
    stop(
      "Unsupported file type. Use CSV, TSV, TXT, their gzip variants, or a project ZIP.",
      call. = FALSE
    )
  }

  dir.create(destination_dir, recursive = TRUE, showWarnings = FALSE)
  safe_stem <- gsub("[^A-Za-z0-9._-]", "_", basename(fileinfo$name[[1]]))
  target <- file.path(destination_dir, paste0(format(Sys.time(), "%H%M%OS6"), "_", safe_stem))
  if (!file.copy(fileinfo$datapath[[1]], target, overwrite = TRUE)) {
    stop("The uploaded file could not be copied into the session workspace.", call. = FALSE)
  }
  target
}

sir_copy_reference_upload <- function(fileinfo, destination_dir, asset) {
  asset <- match.arg(asset, c("gtf", "transcripts", "translations"))
  if (is.null(fileinfo) || !nrow(fileinfo)) {
    stop(
      "Choose the ",
      switch(
        asset,
        gtf = "annotation GTF",
        transcripts = "transcript FASTA",
        translations = "protein FASTA"
      ),
      " file.",
      call. = FALSE
    )
  }
  if (nrow(fileinfo) != 1L) {
    stop("Upload one reference file in each field.", call. = FALSE)
  }
  if (
    is.finite(fileinfo$size[[1]]) &&
    fileinfo$size[[1]] > SIR_MAX_UPLOAD_MB * 1024^2
  ) {
    stop(
      "The reference file is larger than the configured ",
      sir_fmt_num(SIR_MAX_UPLOAD_MB, 0),
      " MB upload limit.",
      call. = FALSE
    )
  }

  filename <- basename(fileinfo$name[[1]])
  valid_name <- switch(
    asset,
    gtf = grepl("\\.gtf\\.gz$", filename, ignore.case = TRUE),
    transcripts = grepl(
      "\\.(fa|fasta)\\.gz$",
      filename,
      ignore.case = TRUE
    ),
    translations = grepl(
      "\\.(fa|fasta)\\.gz$",
      filename,
      ignore.case = TRUE
    )
  )
  if (!isTRUE(valid_name)) {
    stop(
      switch(
        asset,
        gtf = "The annotation file must end in .gtf.gz.",
        transcripts = "The transcript FASTA must end in .fa.gz or .fasta.gz.",
        translations = "The protein FASTA must end in .fa.gz or .fasta.gz."
      ),
      call. = FALSE
    )
  }

  dir.create(destination_dir, recursive = TRUE, showWarnings = FALSE)
  safe_name <- gsub("[^A-Za-z0-9._-]", "_", filename)
  target <- file.path(destination_dir, safe_name)
  if (!file.copy(fileinfo$datapath[[1]], target, overwrite = TRUE)) {
    stop(
      "The uploaded reference file could not be copied into this session.",
      call. = FALSE
    )
  }
  target
}

sir_read_delimited_path <- function(path) {
  suffix <- sir_upload_suffix(path)
  if (!suffix %in% c(".csv", ".csv.gz", ".tsv", ".tsv.gz", ".txt", ".txt.gz")) {
    stop("Expected a CSV, TSV, TXT, or gzip-compressed delimited table.", call. = FALSE)
  }
  separator <- if (suffix %in% c(".csv", ".csv.gz")) "," else "\t"
  result <- data.table::fread(
    path,
    sep = separator,
    header = TRUE,
    data.table = TRUE,
    showProgress = FALSE,
    na.strings = c("", "NA", "NaN")
  )
  if (!ncol(result)) {
    stop("The uploaded table has no columns.", call. = FALSE)
  }
  result
}

sir_read_uploaded_table <- function(fileinfo, destination_dir) {
  path <- sir_copy_upload(fileinfo, destination_dir)
  sir_read_delimited_path(path)
}

sir_path_is_absolute <- function(path) {
  grepl("^(/|~|[A-Za-z]:[/\\\\])", path)
}

sir_path_within <- function(path, root) {
  normalized_path <- normalizePath(path, winslash = "/", mustWork = FALSE)
  normalized_root <- sub(
    "/+$",
    "",
    normalizePath(root, winslash = "/", mustWork = TRUE)
  )
  identical(normalized_path, normalized_root) ||
    startsWith(normalized_path, paste0(normalized_root, "/"))
}

sir_archive_paths_are_safe <- function(names) {
  names <- gsub("\\\\", "/", as.character(names))
  path_names <- sub("/+$", "", names)
  valid <- nzchar(path_names) &
    !startsWith(path_names, "/") &
    !grepl("^[A-Za-z]:/", path_names) &
    vapply(
      strsplit(path_names, "/", fixed = TRUE),
      function(parts) !any(parts %in% c("..", "")),
      logical(1)
    )
  all(valid)
}

sir_archive_entries <- function(zip_path) {
  listing <- utils::unzip(zip_path, list = TRUE)
  if (!nrow(listing)) {
    stop("The project ZIP is empty.", call. = FALSE)
  }
  if (nrow(listing) > SIR_MAX_ARCHIVE_FILES) {
    stop(
      "The project ZIP contains too many files (limit: ",
      sir_fmt_int(SIR_MAX_ARCHIVE_FILES), ").",
      call. = FALSE
    )
  }

  if (!sir_archive_paths_are_safe(listing$Name)) {
    stop("The project ZIP contains an unsafe or invalid path.", call. = FALSE)
  }

  expanded_bytes <- sum(as.numeric(listing$Length), na.rm = TRUE)
  if (expanded_bytes > SIR_MAX_ARCHIVE_EXPANDED_MB * 1024^2) {
    stop(
      "The expanded project ZIP exceeds the configured ",
      sir_fmt_num(SIR_MAX_ARCHIVE_EXPANDED_MB, 0), " MB limit.",
      call. = FALSE
    )
  }
  listing
}

sir_find_manifest <- function(root) {
  files <- list.files(root, recursive = TRUE, full.names = TRUE, all.files = FALSE)
  names <- tolower(basename(files))
  allowed <- c("manifest.csv", "manifest.tsv", "manifest.txt")
  manifests <- files[names %in% allowed]
  if (!length(manifests)) {
    stop("The project ZIP must contain manifest.csv, manifest.tsv, or manifest.txt.", call. = FALSE)
  }
  if (length(manifests) > 1L) {
    stop("The project ZIP contains more than one manifest file.", call. = FALSE)
  }
  manifests[[1]]
}

sir_resolve_manifest_paths <- function(manifest, manifest_path, archive_root) {
  required <- c("path", "sample_name", "condition")
  missing <- setdiff(required, names(manifest))
  if (length(missing)) {
    stop(
      "The manifest is missing required columns: ",
      paste(missing, collapse = ", "), ".",
      call. = FALSE
    )
  }

  out <- as.data.frame(manifest, stringsAsFactors = FALSE)
  out$path <- trimws(as.character(out$path))
  out$sample_name <- trimws(as.character(out$sample_name))
  out$condition <- trimws(as.character(out$condition))

  if (any(!nzchar(out$path)) || any(!nzchar(out$sample_name)) || any(!nzchar(out$condition))) {
    stop("Manifest path, sample_name, and condition values cannot be blank.", call. = FALSE)
  }
  if (anyDuplicated(out$sample_name)) {
    stop("Manifest sample_name values must be unique.", call. = FALSE)
  }
  if (any(vapply(out$path, sir_path_is_absolute, logical(1)))) {
    stop("Project ZIP manifest paths must be relative to the archive.", call. = FALSE)
  }

  manifest_dir <- dirname(manifest_path)
  resolved <- vapply(
    out$path,
    function(relative_path) {
      candidates <- c(
        file.path(manifest_dir, relative_path),
        file.path(archive_root, relative_path)
      )
      candidate <- candidates[dir.exists(candidates)][1]
      if (!length(candidate) || is.na(candidate)) {
        stop("Sample directory not found in project ZIP: ", relative_path, call. = FALSE)
      }
      candidate <- normalizePath(candidate, winslash = "/", mustWork = TRUE)
      if (!sir_path_within(candidate, archive_root)) {
        stop("A manifest path resolves outside the extracted project.", call. = FALSE)
      }
      candidate
    },
    character(1)
  )
  out$path <- resolved
  out
}

sir_expected_event_file <- function(sample_dir, event_type, count_mode) {
  files <- list.files(sample_dir, full.names = FALSE)
  if (event_type %in% c("AFE", "ALE")) {
    pattern <- paste0("\\.", event_type, "PSI(?:\\.gz)?$")
  } else {
    pattern <- paste0("^", event_type, "\\.MATS\\.", count_mode, "\\.txt(?:\\.gz)?$")
  }
  any(grepl(pattern, files, ignore.case = FALSE, perl = TRUE))
}

sir_validate_manifest_event_files <- function(manifest, event_types, count_mode) {
  if (!length(event_types)) {
    stop("Choose at least one event type.", call. = FALSE)
  }
  missing <- character()
  for (row in seq_len(nrow(manifest))) {
    for (event_type in event_types) {
      if (!sir_expected_event_file(manifest$path[[row]], event_type, count_mode)) {
        missing <- c(missing, paste0(manifest$sample_name[[row]], " / ", event_type))
      }
    }
  }
  if (length(missing)) {
    example <- paste(head(missing, 8L), collapse = ", ")
    more <- if (length(missing) > 8L) paste0(" and ", length(missing) - 8L, " more") else ""
    stop("Missing expected event files for: ", example, more, ".", call. = FALSE)
  }
  invisible(TRUE)
}

sir_extract_project_archive <- function(fileinfo, session_dir, event_types, count_mode) {
  zip_path <- sir_copy_upload(fileinfo, session_dir)
  if (!identical(sir_upload_suffix(zip_path), ".zip")) {
    stop("Project ingestion requires a ZIP archive.", call. = FALSE)
  }
  sir_archive_entries(zip_path)

  extract_root <- tempfile("project-", tmpdir = session_dir)
  dir.create(extract_root, recursive = TRUE, showWarnings = FALSE)
  utils::unzip(zip_path, exdir = extract_root)

  extracted <- list.files(
    extract_root,
    recursive = TRUE,
    full.names = TRUE,
    all.files = TRUE,
    include.dirs = TRUE,
    no.. = TRUE
  )
  if (length(extracted) && any(nzchar(Sys.readlink(extracted)))) {
    stop("Symbolic links are not allowed in project ZIP files.", call. = FALSE)
  }
  if (length(extracted) && any(!vapply(extracted, sir_path_within, logical(1), root = extract_root))) {
    stop("The extracted project contains a path outside its session workspace.", call. = FALSE)
  }

  manifest_path <- sir_find_manifest(extract_root)
  manifest <- sir_read_delimited_path(manifest_path)
  manifest <- sir_resolve_manifest_paths(manifest, manifest_path, extract_root)
  sir_validate_manifest_event_files(manifest, event_types, count_mode)

  list(
    manifest = manifest,
    archive_name = basename(fileinfo$name[[1]]),
    extract_root = extract_root
  )
}

sir_validate_ppi_table <- function(object) {
  dt <- data.table::as.data.table(object)
  missing <- setdiff(c("geneA", "geneB"), names(dt))
  if (length(missing)) {
    stop(
      "The PPI table is missing required columns: ",
      paste(missing, collapse = ", "), ".",
      call. = FALSE
    )
  }
  if (!nrow(dt)) {
    stop("The PPI table is empty.", call. = FALSE)
  }
  has_ddi <- all(c("DDI", "DDI_A", "DDI_B") %in% names(dt))
  has_dmi <- all(c("DMI", "DMI_A", "DMI_B") %in% names(dt))
  if (!has_ddi && !has_dmi) {
    stop(
      "A PPI switch table needs either DDI/DDI_A/DDI_B or DMI/DMI_A/DMI_B evidence columns.",
      call. = FALSE
    )
  }
  dt
}

sir_validate_normalized_features <- function(object) {
  dt <- data.table::as.data.table(object)
  required <- c(
    "ensembl_transcript_id", "start", "stop", "database",
    "feature_id", "name", "alt_name"
  )
  missing <- setdiff(required, names(dt))
  if (length(missing)) {
    stop(
      "The normalized feature table is missing required columns: ",
      paste(missing, collapse = ", "), ".",
      call. = FALSE
    )
  }
  dt
}
