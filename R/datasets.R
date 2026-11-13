# Resolve against the active namespace, not the caller's system.file shim.
# pkgload records the source directory as its namespace path; installations
# record the installed package directory. The model is read only on access.
flexbayes_model_dataset = function(name) {
  root = getNamespaceInfo(asNamespace("flexbayes"), "path")
  filename = paste0(name, ".rda")
  directory = if (file.exists(file.path(root, "DESCRIPTION")) &&
                  !file.exists(file.path(root, "Meta", "package.rds"))) "inst/extdata" else "extdata"
  fixture = file.path(root, directory, filename)
  if (!file.exists(fixture)) stop("Fitted dataset file not found: ", name, call. = FALSE)
  environment = new.env(parent = baseenv())
  load(fixture, envir = environment)
  get(name, envir = environment, inherits = FALSE)
}
