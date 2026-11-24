# Download only published inputs and verify the exact snapshots used.
source("R/lib_inputs.R")
options(timeout = max(180, getOption("timeout")))
inputs <- input_manifest()
for (i in seq_len(nrow(inputs))) {
  x <- inputs[i, ]
  if (file.exists(x$path)) {
    verify_input_files(x)
    next
  }
  dir.create(dirname(x$path), recursive = TRUE, showWarnings = FALSE)
  temporary <- paste0(x$path, ".part")
  tryCatch({
    status <- download.file(x$source, temporary, mode = "wb", quiet = TRUE)
    if (status != 0L || !file.exists(temporary) ||
        file.size(temporary) != x$bytes || unname(tools::md5sum(temporary)) != x$md5)
      stop("Download differs from the pinned snapshot: ", x$path)
    if (!file.rename(temporary, x$path)) stop("Cannot install input: ", x$path)
  }, finally = if (file.exists(temporary)) unlink(temporary))
  cat("Verified", x$path, "\n")
}
verify_input_files(inputs)
