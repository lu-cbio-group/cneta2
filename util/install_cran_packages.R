#!/usr/bin/env Rscript

# Install the R packages used by util/ that are on CRAN but not on conda,
# into the active environment (run after `conda activate cneta-util`):
#
#   Rscript util/install_cran_packages.R
#
#   Quartet  compare_trees.R
#   rwty     check_convergence.R
#
# Packages already installed are skipped, so it is safe to rerun.

pkgs = c("Quartet", "rwty")
repo = "https://cloud.r-project.org"

missing = pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing) > 0) {
  cat("Installing from CRAN:", paste(missing, collapse = ", "), "\n")
  install.packages(missing, repos = repo, lib = .libPaths()[1])
}

ok = vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)
for (p in pkgs) {
  if (ok[[p]]) {
    cat(sprintf("OK      %s %s\n", p, as.character(packageVersion(p))))
  } else {
    cat(sprintf("FAILED  %s (see the messages above)\n", p))
  }
}
if (!all(ok)) quit(status = 1)
