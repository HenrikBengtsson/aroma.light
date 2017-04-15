library("devtools")

availableCores <- function() {
  getenv <- function(name) {
    as.integer(Sys.getenv(name, NA_character_))
  }
  getopt <- function(name) {
    as.integer(getOption(name, NA_integer_))
  }
  if (is.finite(n <- getopt("mc.cores") + 1L)) return(n)
  if (is.finite(n <- getopt("Ncpus") + 1L)) return(n)
  if (is.finite(n <- getenv("PBS_NUM_PPN"))) return(n)
  if (is.finite(n <- getenv("SLURM_CPUS_PER_TASK"))) return(n)
  if (is.finite(n <- getenv("NSLOTS"))) return(n)
  1L
}

## Some Bioconductor packages test with more than the two-parallel processes
## accepted by CRAN / --as-cran.
env_vars = c("_R_CHECK_LIMIT_CORES_" = "false")
revdep_check(bioconductor = TRUE, env_vars = env_vars, recursive = TRUE, threads = availableCores())
revdep_check_save_summary()
revdep_check_print_problems()
