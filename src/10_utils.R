sir_fmt_int <- function(x) {
  if (!length(x) || is.na(x[[1]])) {
    return("0")
  }
  format(as.integer(x[[1]]), big.mark = ",", scientific = FALSE, trim = TRUE)
}

sir_fmt_num <- function(x, digits = 2, suffix = "") {
  if (!length(x) || !is.finite(x[[1]])) {
    return("NA")
  }
  paste0(format(round(x[[1]], digits), nsmall = digits, trim = TRUE), suffix)
}

sir_nonempty_table <- function(x) {
  !is.null(x) && is.data.frame(x) && nrow(x) > 0L
}

sir_present_event_types <- function(object) {
  if (
    is.null(object) ||
    !is.data.frame(object) ||
    !"event_type" %in% names(object)
  ) {
    return(character())
  }

  present <- trimws(as.character(object[["event_type"]]))
  present <- unique(present[!is.na(present) & nzchar(present)])
  c(
    SIR_SUPPORTED_EVENTS[SIR_SUPPORTED_EVENTS %in% present],
    sort(setdiff(present, SIR_SUPPORTED_EVENTS))
  )
}

sir_event_type_choices <- function(event_types) {
  event_types <- unique(as.character(event_types))
  event_types <- event_types[!is.na(event_types) & nzchar(event_types)]
  labels <- unname(SIR_EVENT_LABELS[event_types])
  labels[is.na(labels) | !nzchar(labels)] <- event_types[is.na(labels) | !nzchar(labels)]
  stats::setNames(event_types, labels)
}

sir_filter_event_type <- function(object, selected_event_type) {
  if (is.null(object) || !is.data.frame(object)) {
    stop("Event data must be a table.", call. = FALSE)
  }
  if (!"event_type" %in% names(object)) {
    stop("Event data is missing the event_type column.", call. = FALSE)
  }

  selected_event_type <- trimws(as.character(selected_event_type %||% ""))
  if (length(selected_event_type) != 1L || !nzchar(selected_event_type)) {
    stop("Choose an event type.", call. = FALSE)
  }

  dt <- data.table::as.data.table(object)
  event_values <- trimws(as.character(dt[["event_type"]]))
  row_index <- which(!is.na(event_values) & event_values == selected_event_type)
  data.table::copy(dt[row_index])
}

sir_signature <- function(...) {
  digest::digest(list(...), algo = "xxhash64", serialize = TRUE)
}

sir_error_text <- function(error, fallback = "The operation could not be completed.") {
  message <- trimws(conditionMessage(error) %||% "")
  if (!nzchar(message)) fallback else message
}

sir_stringify_value <- function(value, max_chars = 180L) {
  if (is.null(value) || !length(value)) {
    return("")
  }
  text <- paste(head(as.character(unlist(value, use.names = FALSE)), 8L), collapse = "; ")
  if (nchar(text) > max_chars) {
    paste0(substr(text, 1L, max_chars - 3L), "...")
  } else {
    text
  }
}

sir_preview_data <- function(
  object,
  n = 25L,
  query = "",
  columns = NULL,
  search_columns = NULL
) {
  if (is.null(object) || !is.data.frame(object)) {
    return(list(data = data.table::data.table(), total = 0L, matching = 0L))
  }

  dt <- data.table::as.data.table(object)
  total <- nrow(dt)
  n <- max(1L, min(as.integer(n), 250L))
  query <- trimws(query %||% "")
  columns <- intersect(columns %||% names(dt), names(dt))
  if (!length(columns)) {
    columns <- names(dt)
  }

  if (nzchar(query)) {
    candidates <- intersect(search_columns %||% columns, names(dt))
    keep <- rep(FALSE, total)
    for (column in candidates) {
      values <- dt[[column]]
      if (is.list(values)) {
        values <- vapply(values, sir_stringify_value, character(1))
      } else {
        values <- as.character(values)
      }
      hits <- grepl(tolower(query), tolower(values), fixed = TRUE)
      hits[is.na(hits)] <- FALSE
      keep <- keep | hits
    }
    row_index <- which(keep)
  } else {
    row_index <- seq_len(total)
  }

  matching <- length(row_index)
  row_index <- head(row_index, n)
  out <- data.table::copy(dt[row_index, ..columns])
  for (column in names(out)) {
    if (is.list(out[[column]])) {
      out[, (column) := vapply(get(column), sir_stringify_value, character(1))]
    } else if (is.character(out[[column]])) {
      out[, (column) := vapply(get(column), sir_stringify_value, character(1))]
    }
  }

  list(data = out, total = total, matching = matching)
}

sir_preview_caption <- function(preview) {
  if (preview$total == preview$matching) {
    paste0(
      "Showing ", sir_fmt_int(nrow(preview$data)), " of ",
      sir_fmt_int(preview$total), " rows."
    )
  } else {
    paste0(
      "Showing ", sir_fmt_int(nrow(preview$data)), " of ",
      sir_fmt_int(preview$matching), " matching rows (",
      sir_fmt_int(preview$total), " total)."
    )
  }
}

sir_dt_widget <- function(preview, page_length = NULL) {
  if (!nrow(preview$data)) {
    return(NULL)
  }
  DT::datatable(
    preview$data,
    rownames = FALSE,
    filter = "none",
    escape = TRUE,
    options = list(
      dom = "tip",
      pageLength = page_length %||% min(25L, nrow(preview$data)),
      lengthChange = FALSE,
      scrollX = TRUE,
      autoWidth = TRUE
    )
  )
}

sir_metric_card <- function(label, value, detail, tone = "teal") {
  div(
    class = paste("sir-metric", paste0("sir-metric--", tone)),
    div(class = "sir-metric__label", label),
    div(class = "sir-metric__value", value),
    div(class = "sir-metric__detail", detail)
  )
}

sir_status_row <- function(label, ready, detail, dirty = FALSE) {
  state <- if (dirty) "outdated" else if (ready) "ready" else "waiting"
  div(
    class = "sir-status-row",
    div(
      class = paste("sir-status-dot", paste0("sir-status-dot--", state)),
      `aria-hidden` = "true"
    ),
    div(
      class = "sir-status-copy",
      div(class = "sir-status-label", label),
      div(class = "sir-status-detail", detail)
    ),
    span(class = paste("sir-status-badge", paste0("sir-status-badge--", state)), state)
  )
}

sir_panel_note <- function(title, text, tone = "info") {
  div(
    class = paste("sir-note", paste0("sir-note--", tone)),
    div(class = "sir-note__title", title),
    div(class = "sir-note__text", text)
  )
}

sir_download_name <- function(stem, extension) {
  paste0(stem, "_", format(Sys.Date(), "%Y%m%d"), ".", extension)
}
