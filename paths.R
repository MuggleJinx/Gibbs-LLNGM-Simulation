# Path helpers for running scripts either from this directory as a standalone
# repository or from the parent manuscript repository.

simulation_dir <- function() {
  candidates <- c(".", "Simulation")
  for (candidate in candidates) {
    if (file.exists(file.path(candidate, "gibbs.R"))) {
      return(normalizePath(candidate, mustWork = TRUE))
    }
  }
  stop("Could not locate the simulation directory. Run scripts from this repository root or from its parent directory.")
}

sim_path <- function(...) {
  file.path(simulation_dir(), ...)
}
