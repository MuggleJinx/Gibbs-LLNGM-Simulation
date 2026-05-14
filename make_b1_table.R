# Generate LaTeX table for B1 summary from experiment_summary.csv
# Outputs: Table \ref{tab:sim-b1-summary}

source(if (file.exists("paths.R")) "paths.R" else file.path("Simulation", "paths.R"))

in_path <- sim_path("experiment_summary.csv")

if (!file.exists(in_path)) {
  stop("Missing input: ", in_path)
}

x <- read.csv(in_path, stringsAsFactors = FALSE)

need <- c(
  "regime", "point",
  "iact_S_plus", "iact_S_minus", "iact_S_log",
  "ess_sec_S_plus", "ess_sec_S_minus", "ess_sec_S_log"
)
miss <- setdiff(need, names(x))
if (length(miss) > 0) {
  stop("Missing columns in CSV: ", paste(miss, collapse = ", "))
}

fmt2 <- function(v) sprintf("%.2f", as.numeric(v))

rows <- data.frame(
  Regime = x$regime,
  Point = x$point,
  IACT_S_plus = fmt2(x$iact_S_plus),
  IACT_S_minus = fmt2(x$iact_S_minus),
  IACT_S_log = fmt2(x$iact_S_log),
  ESS_sec_S_plus = fmt2(x$ess_sec_S_plus),
  ESS_sec_S_minus = fmt2(x$ess_sec_S_minus),
  ESS_sec_S_log = fmt2(x$ess_sec_S_log),
  stringsAsFactors = FALSE
)

cat("\\begin{table}[t]\n")
cat("\\centering\n")
cat("\\caption{B1 simulation summary (IACT and ESS/sec).}\n")
cat("\\label{tab:sim-b1-summary}\n")
cat("\\begin{tabular}{llrrrrrr}\n")
cat("\\toprule\n")
cat("Regime & Point & IACT $S_+$ & IACT $S_-$ & IACT $S_{\\log}$ & ESS/sec $S_+$ & ESS/sec $S_-$ & ESS/sec $S_{\\log}$ \\\\\n")
cat("\\midrule\n")

for (i in seq_len(nrow(rows))) {
  r <- rows[i, ]
  cat(sprintf(
    "%s & %s & %s & %s & %s & %s & %s & %s \\\\\n",
    r$Regime, r$Point, r$IACT_S_plus, r$IACT_S_minus, r$IACT_S_log,
    r$ESS_sec_S_plus, r$ESS_sec_S_minus, r$ESS_sec_S_log
  ))
}

cat("\\bottomrule\n")
cat("\\end{tabular}\n")
cat("\\end{table}\n")
