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
  # Inside conda, install only into the environment. A cneta-util created
  # before environment-util.yml set R_LIBS_USER still puts the personal R
  # folder first, and packages installed there may not load with this R.
  lib = .libPaths()[1]
  prefix = Sys.getenv("CONDA_PREFIX")
  if (nzchar(prefix) &&
      !startsWith(normalizePath(lib), paste0(normalizePath(prefix, mustWork = FALSE), "/"))) {
    cat("Stopped: R would install into ", lib, "\n", sep = "")
    cat("which is outside the active conda environment ", prefix, "\n", sep = "")
    cat("Packages there may not load with this environment's R.\n")
    cat("If cneta-util is active, update it so that R skips that folder,\n")
    cat("then run this script again:\n")
    cat("  conda env update -f envs/environment-util.yml\n")
    cat("  conda deactivate\n")
    cat("  conda activate cneta-util\n")
    quit(status = 1)
  }
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
