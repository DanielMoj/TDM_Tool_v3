
# R/auth.R (hardened)
# Uses libsodium for password hashing (argon2) and verifies stored hashes.
# Falls back to plaintext only if 'password_hash' missing (discouraged).

credentials_load <- function(path = "config/users.yaml") {
  if (!file.exists(path)) stop("users.yaml fehlt")
  if (!requireNamespace("yaml", quietly = TRUE)) stop("Bitte Paket 'yaml' installieren.")
  yaml::read_yaml(path)
}

password_hash <- function(p) {
  if (!requireNamespace("sodium", quietly = TRUE)) stop("Bitte Paket 'sodium' installieren.")
  sodium::password_store(charToRaw(p))
}

password_verify <- function(p, hash) {
  if (!requireNamespace("sodium", quietly = TRUE)) stop("Bitte Paket 'sodium' installieren.")
  tryCatch(sodium::password_verify(hash, charToRaw(p)), error = function(e) FALSE)
}

auth_check <- function(username, password, path = "config/users.yaml") {
  users <- credentials_load(path)
  if (is.null(users$users)) stop("users.yaml: Abschnitt 'users' fehlt")
  u <- NULL
  for (usr in users$users) if (identical(usr$username, username)) { u <- usr; break }
  if (is.null(u)) return(FALSE)
  if (!is.null(u$password_hash)) return(isTRUE(password_verify(password, u$password_hash)))
  # Fallback to plaintext (discouraged)
  if (!is.null(u$password)) return(identical(as.character(u$password), as.character(password)))
  FALSE
}

# Utility: upgrade users.yaml by hashing plaintext passwords (one-time action)
auth_upgrade_hashes <- function(path = "config/users.yaml", backup = TRUE) {
  if (!requireNamespace("yaml", quietly = TRUE)) stop("Bitte Paket 'yaml' installieren.")
  users <- credentials_load(path)
  changed <- FALSE
  for (i in seq_along(users$users)) {
    u <- users$users[[i]]
    if (!is.null(u$password) && is.null(u$password_hash)) {
      u$password_hash <- as.character(password_hash(u$password))
      u$password <- NULL
      users$users[[i]] <- u; changed <- TRUE
    }
  }
  if (changed) {
    if (backup) file.copy(path, paste0(path, ".bak"), overwrite = TRUE)
    yaml::write_yaml(users, path)
  }
  invisible(changed)
}


# ---- Roles & Policies ----
policies_load <- function(path = "config/policies.yaml") {
  if (!requireNamespace("yaml", quietly = TRUE)) stop("Bitte Paket 'yaml' installieren.")
  if (!file.exists(path)) stop("policies.yaml fehlt")
  yaml::read_yaml(path)
}

policy_allow <- function(role, action, policies = NULL) {
  if (is.null(policies)) policies <- policies_load()
  if (is.null(role) || !nzchar(role)) return(FALSE)
  allow <- try(policies$roles[[role]]$allow, silent = TRUE)
  if (inherits(allow, "try-error") || is.null(allow)) return(FALSE)
  action %in% allow
}

# ---- Session / TTL ----
session_config <- function(path = "config/session.yml") {
  if (!requireNamespace("yaml", quietly = TRUE)) return(list(ttl_minutes = 30L))
  if (!file.exists(path)) return(list(ttl_minutes = 30L))
  yaml::read_yaml(path)
}

auth_is_authenticated <- function(session) {
  !is.null(session$userData$user()) && nzchar(session$userData$user())
}

auth_user_role <- function(session) {
  session$userData$role()
}

auth_set_user <- function(session, username, role) {
  session$userData$user <- reactiveVal(username)
  session$userData$role <- reactiveVal(role)
  session$userData$login_time <- reactiveVal(Sys.time())
  cfg <- session_config()
  ttl <- as.numeric(cfg$ttl_minutes %||% 30)
  session$userData$expire_at <- reactiveVal(Sys.time() + ttl*60)
}

auth_logout <- function(session, reason = "logout") {
  audit_event("logout", list(reason = reason), session = session)
  session$userData$user <- reactiveVal(NULL)
  session$userData$role <- reactiveVal(NULL)
  session$userData$login_time <- reactiveVal(NULL)
  session$userData$expire_at <- reactiveVal(NULL)
}

login_modal <- function(title = "Login") {
  shiny::modalDialog(
    title = title,
    shiny::textInput("auth_user", "Benutzername"),
    shiny::passwordInput("auth_pass", "Passwort"),
    footer = shiny::tagList(
      shiny::modalButton("Abbrechen"),
      shiny::actionButton("auth_do_login", "Anmelden", class = "btn-primary")
    ),
    easyClose = FALSE
  )
}

# ---- E-Sign ----
esign_modal <- function(action = "sign", msg = "Bitte bestätigen Sie mit Passwort") {
  shiny::modalDialog(
    title = paste("E-Signatur:", action),
    shiny::passwordInput("esign_pass", "Passwort"),
    footer = shiny::tagList(
      shiny::modalButton("Abbrechen"),
      shiny::actionButton("esign_confirm", "Signieren", class = "btn-primary")
    ),
    easyClose = FALSE
  )
}

auth_esign_verify <- function(session, password) {
  u <- session$userData$user(); if (is.null(u) || !nzchar(u)) return(FALSE)
  users <- credentials_load()
  for (usr in users$users) {
    if (identical(usr$username, u)) {
      if (!is.null(usr$password_hash)) return(password_verify(password, usr$password_hash))
      if (!is.null(usr$password)) return(identical(as.character(usr$password), as.character(password)))
    }
  }
  FALSE
}
