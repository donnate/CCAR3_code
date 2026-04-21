repo_root <- "~/Documents/CCAR3_code"

write_plot_file <- function(plot_obj, relative_path, width = 11, height = 7) {
  out_path <- file.path(repo_root, relative_path)
  dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(out_path, plot = plot_obj, width = width, height = height, dpi = 300)
  out_path
}

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 1L) {
  stop("Usage: Rscript plot_rrr_effect_nnzeros_compare.R <results_csv>", call. = FALSE)
}

results_csv <- args[1]
name_stub <- if (length(args) >= 2L) args[2] else tools::file_path_sans_ext(basename(results_csv))

suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(dplyr))

results <- read.csv(results_csv, stringsAsFactors = FALSE)

summary_results <- results %>%
  dplyr::filter(status == "ok", method %in% c("ccar3", "FIT_SAR_CV")) %>%
  dplyr::group_by(method, p2, nnzeros) %>%
  dplyr::summarise(
    distance_tot_mean = mean(distance_tot, na.rm = TRUE),
    distance_tot_sd = stats::sd(distance_tot, na.rm = TRUE),
    distance_tot_se = distance_tot_sd / sqrt(dplyr::n()),
    distance_tot_ci_lo = distance_tot_mean - 1.96 * distance_tot_se,
    distance_tot_ci_hi = distance_tot_mean + 1.96 * distance_tot_se,
    counts = dplyr::n(),
    .groups = "drop"
  )

method_colors <- c("ccar3" = "#111111", "FIT_SAR_CV" = "#0B8F5A")
method_labels <- c("ccar3" = "ccar3", "FIT_SAR_CV" = "FIT_SAR_CV")

plot_obj <- ggplot(
  summary_results,
  aes(x = nnzeros, y = distance_tot_mean, color = method, group = method)
) +
  geom_errorbar(
    aes(ymin = distance_tot_ci_lo, ymax = distance_tot_ci_hi),
    width = 0.8,
    alpha = 0.8
  ) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2) +
  facet_wrap(
    ~p2,
    nrow = 2,
    scales = "free_y",
    labeller = as_labeller(c(`5` = "q = 5", `10` = "q = 10", `20` = "q = 20", `30` = "q = 30"))
  ) +
  scale_color_manual(values = method_colors, labels = method_labels) +
  labs(
    title = "Mean Subspace Distance vs Number of Nonzeros (error bars = 95% CI)",
    x = "Number of nonzero rows in U",
    y = "Mean subspace distance",
    color = "Method"
  ) +
  theme_bw(base_size = 13) +
  theme(
    legend.position = "top",
    panel.grid.minor = element_blank()
  )

png_path <- write_plot_file(
  plot_obj,
  file.path("experiments", "simulations", "results", paste0(name_stub, "_facetted.png"))
)
pdf_path <- write_plot_file(
  plot_obj,
  file.path("experiments", "simulations", "results", paste0(name_stub, "_facetted.pdf"))
)

cat("Wrote plot PNG to: ", png_path, "\n", sep = "")
cat("Wrote plot PDF to: ", pdf_path, "\n", sep = "")
