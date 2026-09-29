fixture_path <- function(...) {
  parts <- c("extdata", ...)
  installed <- do.call(system.file, c(as.list(parts), list(package = "immunread")))
  if (nzchar(installed)) return(installed)

  path <- do.call(testthat::test_path, as.list(c("..", "..", "inst", parts[-1L])))
  if (!file.exists(path)) stop("Fixture not found: ", path)
  path
}

read_fixture <- function(reader, fixture, schema = c("v_call", "junction_aa"), ...) {
  reader(do.call(fixture_path, as.list(fixture)), schema = schema, output_folder = tempfile(), ...)
}

options(immundata.verbose = FALSE)

annotations <- function(idata) dplyr::collect(idata)
