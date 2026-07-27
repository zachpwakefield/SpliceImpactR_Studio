expected_sha <- "5436813f56606af40cc33e917f2b8fb75bb90c12"

deploy_source <- tryCatch(sys.frame(1)$ofile, error = function(error) NULL)
if (is.null(deploy_source) || !nzchar(deploy_source)) {
  file_argument <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  deploy_source <- if (length(file_argument)) {
    sub("^--file=", "", file_argument[[1]])
  } else {
    NULL
  }
}
deploy_dir <- if (!is.null(deploy_source) && nzchar(deploy_source)) {
  dirname(normalizePath(deploy_source, winslash = "/", mustWork = TRUE))
} else {
  normalizePath(".", winslash = "/", mustWork = TRUE)
}

app_files <- c(
  "app.R",
  "src",
  "www",
  "DESCRIPTION",
  "README.md",
  "renv.lock"
)

reference_files <- c(
  "human_gencode_v45.gtf.rds",
  "human_gencode_v45_sequences.rds",
  "human_gencode_v45_hybrids.rds"
)

missing <- app_files[!file.exists(file.path(deploy_dir, app_files))]
if (length(missing)) {
  stop("v4 deployment files are missing: ", paste(missing, collapse = ", "))
}

reference_present <- file.exists(file.path(deploy_dir, reference_files))
if (any(reference_present) && !all(reference_present)) {
  stop(
    "The optional GENCODE reference bundle is incomplete. Provide all three ",
    "preprocessed files or remove the partial bundle: ",
    paste(reference_files[!reference_present], collapse = ", "),
    "."
  )
}
if (all(reference_present)) {
  app_files <- c(app_files, reference_files)
  message("Including the local preprocessed human GENCODE 45 reference.")
} else {
  message(
    "No local preprocessed reference found; deploying the source-only ",
    "test and browser-upload workflows."
  )
}

lockfile <- jsonlite::read_json(
  file.path(deploy_dir, "renv.lock"),
  simplifyVector = FALSE
)
locked_sha <- lockfile$Packages$SpliceImpactR$RemoteSha
if (!identical(locked_sha, expected_sha)) {
  found_sha <- if (is.null(locked_sha) || !length(locked_sha)) "<none>" else locked_sha
  stop(
    "renv.lock does not pin the expected SpliceImpactR commit. Expected ",
    expected_sha, "; found ", found_sha, "."
  )
}

installed <- utils::packageDescription("SpliceImpactR")
if (!identical(installed$RemoteSha, expected_sha)) {
  stop(
    "The local deployment runtime is not using the pinned SpliceImpactR commit. ",
    "Install fiszbein-lab/SpliceImpactR@", expected_sha, " before deploying."
  )
}

dry_run <- tolower(trimws(Sys.getenv(
  "SPLICEIMPACTR_DEPLOY_DRY_RUN",
  unset = "false"
))) %in% c("1", "true", "yes", "on")
deploy_account <- trimws(Sys.getenv(
  "SPLICEIMPACTR_RSCONNECT_ACCOUNT",
  unset = "fiszbein-lab"
))
deploy_server <- trimws(Sys.getenv(
  "SPLICEIMPACTR_RSCONNECT_SERVER",
  unset = "shinyapps.io"
))
deploy_app_name <- trimws(Sys.getenv(
  "SPLICEIMPACTR_RSCONNECT_APP_NAME",
  unset = "SpliceImpactR_Studio"
))

message(
  "Deployment target: ",
  deploy_account,
  " / ",
  deploy_server,
  " / ",
  deploy_app_name,
  "."
)

if (isTRUE(dry_run)) {
  paths <- rsconnect::listDeploymentFiles(
    appDir = deploy_dir,
    appFiles = app_files
  )
  message(
    "Deployment validation passed for ",
    length(paths),
    " source-only or locally bundled files. No upload was performed."
  )
} else {
  rsconnect::deployApp(
    appDir = deploy_dir,
    appPrimaryDoc = "app.R",
    appName = deploy_app_name,
    appFiles = app_files,
    account = deploy_account,
    server = deploy_server,
    forceUpdate = TRUE,
    launch.browser = FALSE
  )
}
