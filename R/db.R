# R/db.R
# Lightweight DB helpers for Postgres; falls back to local CSV if no DSN provided
get_db_con <- function() {
  dsn <- Sys.getenv("PG_DSN", "")
  if (dsn == "") return(NULL)
  if (!requireNamespace("DBI", quietly = TRUE) || !requireNamespace("RPostgres", quietly = TRUE)) {
    warning("DB packages fehlen, kehre zu NULL (kein DB) zurück.")
    return(NULL)
  }
  con <- tryCatch(DBI::dbConnect(RPostgres::Postgres(), dsn = dsn), error = function(e) NULL)
  con
}

db_write_audit <- function(user, role, event, details = list()) {
  con <- get_db_con()
  if (is.null(con)) return(invisible(FALSE))
  try({
    DBI::dbExecute(con, "INSERT INTO audit_log(user_name, role, event, details) VALUES($1,$2,$3,$4::jsonb)",
                   params = list(user %||% "guest", role %||% "guest", event, jsonlite::toJSON(details, auto_unbox = TRUE)))
    DBI::dbDisconnect(con)
  }, silent = TRUE)
  invisible(TRUE)
}


# R/db.R - Postgres helpers for antibiogram & dataset_versions

connect_pg <- function() {
  if (!requireNamespace("DBI", quietly = TRUE) || !requireNamespace("RPostgres", quietly = TRUE)) {
    stop("Bitte Pakete 'DBI' und 'RPostgres' installieren.")
  }
  host <- Sys.getenv("PGHOST", "localhost")
  port <- as.integer(Sys.getenv("PGPORT", "5432"))
  db   <- Sys.getenv("PGDATABASE", "tdmx")
  user <- Sys.getenv("PGUSER", "tdmx")
  pass <- Sys.getenv("PGPASSWORD", "")
  DBI::dbConnect(RPostgres::Postgres(), host = host, port = port, dbname = db, user = user, password = pass)
}

db_write_dataset_version <- function(con, kind, version, checksum, meta = list()) {
  stopifnot(!is.null(kind), !is.null(version), !is.null(checksum))
  sql <- "INSERT INTO dataset_versions(kind, version, checksum, meta) VALUES ($1,$2,$3,$4)"
  DBI::dbExecute(con, sql, params = list(kind, version, checksum, jsonlite::toJSON(meta, auto_unbox = TRUE)))
}

db_import_antibiogram <- function(con, df, source = "upload", version = NULL) {
  stopifnot(all(c("drug","mic","prob") %in% colnames(df)))
  # Normalize per drug (ensure sum(prob)=1)
  df <- dplyr::group_by(df, drug) |> dplyr::mutate(prob = prob / sum(prob)) |> dplyr::ungroup()
  # Insert rows
  df$source <- source
  DBI::dbWriteTable(con, "antibiogram", df, append = TRUE, row.names = FALSE)
  # Version entry
  if (!requireNamespace("digest", quietly = TRUE)) stop("Bitte Paket 'digest' installieren.")
  checksum <- digest::digest(jsonlite::toJSON(df, digits = NA), algo = "sha256")
  version <- version %||% format(Sys.time(), "%Y%m%d%H%M%S")
  try(db_write_dataset_version(con, "antibiogram", version, checksum, meta = list(source = source)), silent = TRUE)
  invisible(list(version = version, checksum = checksum, n = nrow(df)))
}

db_list_antibiogram_drugs <- function(con) {
  tb <- DBI::dbGetQuery(con, "SELECT DISTINCT drug FROM antibiogram ORDER BY drug ASC")
  tb$drug
}

db_get_antibiogram <- function(con, drug) {
  stopifnot(nzchar(drug))
  DBI::dbGetQuery(con, "SELECT drug, mic, prob FROM antibiogram WHERE drug = $1 ORDER BY mic ASC", params = list(drug))
}
