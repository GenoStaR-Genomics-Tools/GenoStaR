# environment to cache the snapshot in memory
.allele_cache <- new.env(parent = emptyenv())

# Load the snapshot into cache
load_bundled_snapshot <- function(version = NULL, pkg = "GenoStaR") {
  datadir <- system.file("data", package = pkg)
  files <- list.files(datadir, pattern = "^allele_definitions_snapshot_.*\\.rda$", full.names = TRUE)
  
  if (length(files) == 0) stop("No snapshot .rda files found in package data/")
  
  if (is.null(version)) {
    # pick the latest by date in filename
    version <- sub(".*_(\\d{4}-\\d{2}-\\d{2})\\.rda$", "\\1", basename(files))
    version <- sort(version, decreasing = TRUE)[1]
  }
  
  filepath <- files[grepl(version, files)]
  if (length(filepath) == 0) stop("Snapshot not found: ", version)
  
  e <- new.env(parent = emptyenv())
  load(filepath, envir = e)
  assign("allele_definitions_snapshot", get("allele_definitions_snapshot", envir = e), envir = .allele_cache)
  assign("snapshot_version", version, envir = .allele_cache)
  
  invisible(version)
}

# Accessor to get the snapshot version
current_snapshot_version <- function() {
  if (exists("snapshot_version", envir = .allele_cache)) {
    get("snapshot_version", envir = .allele_cache)
  } else {
    NA_character_
  }
}