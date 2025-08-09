# R/audit.R
ensure_audit_file <- function(path) {
  if (!file.exists(path)) {
    dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
    readr::write_csv(tibble::tibble(
      timestamp = character(), user = character(), role = character(),
      event = character(), details = character()
    ), path)
  }
}

log_event <- function(path, user_info, event, details = list()) {
  ensure_audit_file(path)
  u <- tryCatch(if (is.function(user_info$user)) user_info$user() else user_info$user, error = function(e) NULL)
  r <- tryCatch(if (is.function(user_info$role)) user_info$role() else user_info$role, error = function(e) NULL)
  df <- tibble::tibble(
    timestamp = as.character(Sys.time()),
    user = u %||% "guest",
    role = r %||% "guest",
    event = event,
    details = jsonlite::toJSON(details, auto_unbox = TRUE)
  )
  readr::write_csv(df, path, append = TRUE)
}


# --- Hash-chained Audit Log ---
audit_append_hashchain <- function(file = "log/audit.csv", actor, action, payload = list()) {
  if (!requireNamespace("digest", quietly = TRUE)) stop("Bitte Paket 'digest' installieren.")
  dir.create(dirname(file), showWarnings = FALSE, recursive = TRUE)
  prev_hash <- "GENESIS"
  if (file.exists(file)) {
    tb <- try(readr::read_csv(file, show_col_types = FALSE, progress = FALSE), silent = TRUE)
    if (!inherits(tb, "try-error") && nrow(tb) > 0) prev_hash <- tail(tb$hash, 1)
  }
  ts <- format(Sys.time(), tz = "UTC", usetz = TRUE)
  payload_json <- jsonlite::toJSON(payload, auto_unbox = TRUE, null = "null")
  chain_input <- paste(prev_hash, ts, actor, action, payload_json, sep = "|")
  h <- digest::digest(chain_input, algo = "sha256")
  df <- data.frame(ts = ts, actor = actor, action = action, payload = payload_json, prev_hash = prev_hash, hash = h)
  readr::write_csv(df, file, append = file.exists(file))
  invisible(h)
}

audit_verify_chain <- function(file = "log/audit.csv") {
  if (!requireNamespace("digest", quietly = TRUE)) stop("Bitte Paket 'digest' installieren.")
  if (!file.exists(file)) return(TRUE)
  tb <- readr::read_csv(file, show_col_types = FALSE, progress = FALSE)
  prev <- "GENESIS"
  for (i in seq_len(nrow(tb))) {
    chain_input <- paste(prev, tb$ts[i], tb$actor[i], tb$action[i], tb$payload[i], sep = "|")
    h <- digest::digest(chain_input, algo = "sha256")
    if (!identical(h, tb$hash[i])) return(FALSE)
    prev <- tb$hash[i]
  }
  TRUE
}


# ---- Central audit_event ----
audit_event <- function(action, payload = list(), session = shiny::getDefaultReactiveDomain(), require_reason = FALSE, require_esign = FALSE) {
  actor <- try(session$userData$user(), silent = TRUE); if (inherits(actor, "try-error") || is.null(actor)) actor <- "anonymous"
  if (require_reason) {
    # Show modal to capture reason; store in payload$reason via input$aud_reason
    shiny::showModal(shiny::modalDialog(
      title = paste("Begründung erforderlich:", action),
      shiny::textAreaInput("aud_reason", "Begründung", rows = 3),
      footer = shiny::tagList(
        shiny::modalButton("Abbrechen"),
        shiny::actionButton("aud_reason_ok", "OK", class = "btn-primary")
      ), easyClose = FALSE
    ))
    # Wait for aud_reason_ok via observeEvent in server; here we append asynchronously in app wiring.
  } else {
    .audit_do_append(actor, action, payload)
  }
}

.aud_db_write <- function(entry) {
  # Optional DB sink
  if (!requireNamespace("DBI", quietly = TRUE) || !requireNamespace("RPostgres", quietly = TRUE)) return(invisible(FALSE))
  dsn_ok <- nzchar(Sys.getenv("PGDATABASE",""))
  if (!dsn_ok) return(invisible(FALSE))
  con <- try(connect_pg(), silent = TRUE); if (inherits(con, "try-error")) return(invisible(FALSE))
  on.exit(try(DBI::dbDisconnect(con), silent = TRUE))
  prev <- try(DBI::dbGetQuery(con, "SELECT hash FROM audit_log ORDER BY id DESC LIMIT 1"), silent = TRUE)
  prev_hash <- if (!inherits(prev, "try-error") && nrow(prev) > 0) prev$hash[1] else "GENESIS"
  ts <- format(Sys.time(), tz = "UTC", usetz = TRUE)
  payload_json <- jsonlite::toJSON(entry$payload, auto_unbox = TRUE, null = "null")
  chain_input <- paste(prev_hash, ts, entry$actor, entry$action, payload_json, sep = "|")
  h <- digest::digest(chain_input, algo = "sha256")
  sql <- "INSERT INTO audit_log(ts, actor, action, payload, prev_hash, hash) VALUES ($1,$2,$3,$4,$5,$6)"
  DBI::dbExecute(con, sql, params = list(ts, entry$actor, entry$action, payload_json, prev_hash, h))
  invisible(TRUE)
}

.audit_do_append <- function(actor, action, payload) {
  entry <- list(actor = actor, action = action, payload = payload)
  try(audit_append_hashchain(actor = actor, action = action, payload = payload), silent = TRUE)
  try(.aud_db_write(entry), silent = TRUE)
  invisible(TRUE)
}
