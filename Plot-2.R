library(ggplot2)
library(tidyr)
library(dplyr)

source(if (file.exists("paths.R")) "paths.R" else file.path("Simulation", "paths.R"))

# ============================================================
# 1) Read data
# ============================================================
file_path <- sim_path("experiment2_summary_scanmu_4stats.csv")
stopifnot(file.exists(file_path))

df <- read.csv(file_path) %>%
  arrange(gamma_ns)

# ============================================================
# 2) Build long-format table for ESS/sec and IACT
# ============================================================
# Select the 4 statistics for ESS/sec and IACT
df_long <- df %>%
  select(
    gamma_ns,
    ess_sec_S_plus,  ess_sec_S_minus,  ess_sec_S_log,  ess_sec_T_null,
    iact_S_plus,     iact_S_minus,     iact_S_log,     iact_T_null
  ) %>%
  pivot_longer(
    cols = -gamma_ns,
    names_to = c("Metric", "Statistic"),
    names_pattern = "^(ess_sec|iact)_(S_plus|S_minus|S_log|T_null)$",
    values_to = "Value"
  ) %>%
  mutate(
    Metric = recode(Metric,
                    "ess_sec" = "ESS/sec",
                    "iact"    = "IACT")
  ) %>%
  filter(Metric == "IACT")

# Optional: order legend as you like
stat_labels <- c(
  "S_plus"  = expression(S["+"]),
  "S_minus" = expression(S["-"]),
  "S_log"   = expression(S[log]),
  "T_null"  = expression(T[null])
)

df_long$Statistic <- factor(
  df_long$Statistic,
  levels = c("S_plus", "S_minus", "S_log", "T_null")
)

# ============================================================
# 3) Plot: 2 panels (ESS/sec and IACT), each with 4 curves
# ============================================================
use_log_ess <- FALSE  # set FALSE if you don't want log-scale on ESS/sec

p <- ggplot(df_long, aes(x = gamma_ns, y = Value, color = Statistic, shape = Statistic)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2.8) +
  # facet_wrap(~ Metric, scales = "free_y", ncol = 2) +
  labs(
    x = expression(gamma[ns](mu)),
    y = "IACT",
    # title = "B2: Mixing efficiency (IACT)",
    # subtitle = "Geometry fixed; drift parameter mu scanned"
  ) +
  theme_bw() +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(size = 15),
    legend.key.size = unit(2, "lines"),
    strip.text = element_text(face = "bold")
  ) +
  scale_color_discrete(labels = stat_labels) +
  scale_shape_discrete(labels = stat_labels)

# Apply log-scale only to the ESS/sec panel (keep IACT linear)
# Easiest and cleanest way: build two plots and combine? We avoid extra packages.
# We'll do it by transforming ESS/sec values in the data when use_log_ess=TRUE.

if (use_log_ess) {
  df_long2 <- df_long %>%
    mutate(Value_plot = ifelse(Metric == "ESS/sec", log10(Value), Value),
           ylab = ifelse(Metric == "ESS/sec", "log10(ESS/sec)", "IACT"))

  p <- ggplot(df_long2, aes(x = gamma_ns, y = Value_plot, color = Statistic, shape = Statistic)) +
    geom_line(linewidth = 1) +
    geom_point(size = 2.8) +
    facet_wrap(~ Metric, scales = "free_y", ncol = 2) +
    labs(
      x = expression(gamma[ns](mu)),
      y = NULL,
      # title = "B2: Mixing efficiency along the null-smallness scan",
      # subtitle = "Left: log10(ESS/sec). Right: IACT. Geometry fixed; mu scanned."
    ) +
    theme_bw() +
    theme(
      legend.position = "bottom",
      legend.title = element_blank(),
      legend.text = element_text(size = 15),
      legend.key.size = unit(2, "lines"),
      strip.text = element_text(face = "bold")
    ) +
    scale_color_discrete(labels = stat_labels) +
    scale_shape_discrete(labels = stat_labels)
}

print(p)

# ============================================================
# 4) Save (single combined figure)
# ============================================================
ggsave(
  filename = sim_path("B2_IACT.png"),
  plot = p,
  width = 8,
  height = 5,
  dpi = 300
)

cat("Saved: ", sim_path("B2_IACT.png"), "\n", sep = "")
