#!/usr/bin/env Rscript
# More views of chain length vs the IFNB response (chain length on the x axis throughout).
# Run 34_lipidomics_chain_length.R first - this reads its lipids_parsed.csv and applies the same filters
# (no standards, no oxidised species, no rows with missing values).
#
# Writes to <project>/lipidomics_chain_length/ (PNG + PDF):
#   A_core_trend_lines          one plot, a fitted line per core lysosomal class (compare slopes directly)
#   B_mean_logFC_by_chain_length  bars of mean log2FC at each chain length, per class (shows WHERE it changes)
#   C_heatmap_class_x_chain_length  every class x chain length, coloured by mean log2FC
#   D_lysosomal_vs_other_trends  pooled trends (carbons per chain) for core / general / other lipids
#   E_chain_length_x_unsaturation  chain length vs double bonds coloured by log2FC (is it length or saturation?)
# log2FC > 0 (red) = higher in IFNB; < 0 (blue) = lower.

suppressPackageStartupMessages({ library(readr); library(dplyr); library(tidyr); library(ggplot2) })

MIN_SPECIES <- 5
PROJ <- Sys.getenv("PROJ", unset = file.path(
  Sys.getenv("HOME"),
  "Documents/aviti_rnaseq/2025/aviti",
  "20251009-100725_Maury_25161_WivfFacICOTHFAUJwSPeHBEFfCfNtCRpUSgVGUkz",
  "results/Maury_25161"))
OUT <- file.path(PROJ, "lipidomics_chain_length")
pf <- file.path(OUT, "lipids_parsed.csv")
if (!file.exists(pf)) stop("Run 34_lipidomics_chain_length.R first (missing ", pf, ")")
x <- read_csv(pf, show_col_types = FALSE) |>
  filter(!standard, !oxidized, !miss, !is.na(total_C), !is.na(lfc), class != "Co")
message("Lipids used: ", nrow(x))

INK <- "#1F2933"; UP <- "#C53030"; DOWN <- "#2B6CB0"
theme_f <- function() theme_bw(base_size = 11) + theme(panel.grid.minor = element_blank(),
  panel.border = element_rect(colour = "grey70"), axis.text = element_text(colour = INK),
  strip.background = element_rect(fill = "grey95", colour = "grey70"), strip.text = element_text(colour = INK, face = "bold"),
  plot.title = element_text(colour = INK, face = "bold"), plot.subtitle = element_text(colour = "grey35", size = 9),
  legend.position = "bottom")
save2 <- function(p, name, w, h) {
  ggsave(file.path(OUT, paste0(name, ".png")), p, width = w, height = h, dpi = 300, limitsize = FALSE)
  ggsave(file.path(OUT, paste0(name, ".pdf")), p, width = w, height = h, limitsize = FALSE)
}
big <- x |> count(class) |> filter(n >= MIN_SPECIES) |> pull(class)
lyso_big <- x |> filter(tier != "Other", class %in% big)
core_big <- x |> filter(tier == "Core lysosomal", class %in% big)

# ---- A: fitted lines per core class ---------------------------------------------------------------------------
if (nrow(core_big)) {
  cols <- c("#D55E00", "#0072B2", "#009E73", "#CC79A7", "#E69F00", "#56B4E9", "#7A7A7A")
  cl <- sort(unique(core_big$class)); pal <- setNames(cols[seq_along(cl)], cl)
  pA <- ggplot(core_big, aes(total_C, lfc, colour = class)) +
    geom_hline(yintercept = 0, colour = "grey50", linewidth = 0.3) +
    geom_point(size = 1.4, alpha = 0.35) +
    geom_smooth(method = "lm", formula = y ~ x, se = FALSE, linewidth = 1.1) +
    scale_colour_manual(values = pal, name = NULL) +
    labs(title = "Core lysosomal-membrane lipids: log2FC vs chain length",
         subtitle = "one fitted line per class; a downward line = shorter chains relatively higher in IFNB",
         x = "total acyl carbons", y = "log2 fold-change (IFNB vs untreated)") + theme_f()
  save2(pA, "A_core_trend_lines", 8, 5.5)
}

# ---- B: mean log2FC at each chain length ------------------------------------------------------------------------
if (nrow(lyso_big)) {
  prof <- lyso_big |> group_by(class, tier, total_C) |>
    summarise(mean_lfc = mean(lfc), se = ifelse(n() > 1, sd(lfc) / sqrt(n()), NA_real_), n = n(), .groups = "drop") |>
    mutate(dir = ifelse(mean_lfc >= 0, "higher in IFNB", "lower in IFNB"))
  pB <- ggplot(prof, aes(total_C, mean_lfc)) +
    geom_hline(yintercept = 0, colour = "grey40", linewidth = 0.3) +
    geom_col(aes(fill = dir), width = 0.8) +
    geom_errorbar(aes(ymin = mean_lfc - se, ymax = mean_lfc + se), width = 0.3, colour = INK, linewidth = 0.3, na.rm = TRUE) +
    geom_smooth(data = lyso_big, aes(total_C, lfc), method = "lm", formula = y ~ x, se = FALSE,
                colour = INK, linetype = "dashed", linewidth = 0.6, inherit.aes = FALSE) +
    facet_wrap(~class, scales = "free_x") +
    scale_fill_manual(values = c(`higher in IFNB` = UP, `lower in IFNB` = DOWN), name = NULL) +
    labs(title = "Mean log2FC at each chain length", x = "total acyl carbons",
         subtitle = "bars = mean of the species with that chain length (+/- SE where more than one); dashed line = linear trend",
         y = "mean log2 fold-change (IFNB vs untreated)") + theme_f()
  nf <- n_distinct(prof$class)
  save2(pB, "B_mean_logFC_by_chain_length", 3.3 * min(nf, 4) + 1, 3 * ceiling(nf / 4) + 1.6)
}

# ---- C: heatmap, every class x chain length ---------------------------------------------------------------------------
tile <- x |> filter(class %in% big) |> group_by(class, tier, total_C) |>
  summarise(mean_lfc = mean(lfc), n = n(), .groups = "drop")
if (nrow(tile)) {
  ord <- tile |> group_by(tier, class) |> summarise(tot = sum(n), .groups = "drop") |> arrange(tier, desc(tot))
  tile$class <- factor(tile$class, levels = rev(ord$class))
  tile$tier <- factor(tile$tier, levels = c("Core lysosomal", "General membrane", "Other"))
  pC <- ggplot(tile, aes(total_C, class, fill = mean_lfc)) +
    geom_tile(colour = "white", linewidth = 0.4) +
    facet_grid(tier ~ ., scales = "free_y", space = "free_y") +
    scale_fill_gradient2(low = "#2C5F9E", mid = "#F7F7F7", high = "#B2382C", midpoint = 0,
                         limits = c(-1.5, 1.5), oob = scales::squish, name = "mean log2FC\n(IFNB vs untreated)") +
    labs(title = "Where along the chain-length range do lipids change?",
         subtitle = "each tile = mean log2FC of the species with that class and total carbons (blank = none measured; colour capped at +/-1.5)",
         x = "total acyl carbons", y = NULL) + theme_f() +
    theme(panel.background = element_rect(fill = "grey95"), panel.grid = element_blank())
  save2(pC, "C_heatmap_class_x_chain_length", 11, 0.32 * n_distinct(tile$class) + 3.5)
}

# ---- D: pooled trends, lysosomal vs other -----------------------------------------------------------------------------
pool <- x |> filter(class %in% big, !is.na(carbons_per_chain))
if (nrow(pool)) {
  pool$tier <- factor(pool$tier, levels = c("Core lysosomal", "General membrane", "Other"))
  pD <- ggplot(pool, aes(carbons_per_chain, lfc, colour = tier, fill = tier)) +
    geom_hline(yintercept = 0, colour = "grey50", linewidth = 0.3) +
    geom_point(size = 1, alpha = 0.18) +
    geom_smooth(method = "lm", formula = y ~ x, linewidth = 1, alpha = 0.15) +
    scale_colour_manual(values = c(`Core lysosomal` = UP, `General membrane` = "#6B7280", Other = "#2F855A"), name = NULL) +
    scale_fill_manual(values = c(`Core lysosomal` = UP, `General membrane` = "#6B7280", Other = "#2F855A"), name = NULL) +
    coord_cartesian(ylim = c(-2, 2)) +
    labs(title = "Chain-length trend: lysosomal-membrane lipids vs the rest",
         subtitle = "all species pooled by group (view clipped to +/-2 log2FC); x = carbons per chain (total carbons / number of chains)",
         x = "carbons per chain", y = "log2 fold-change (IFNB vs untreated)") + theme_f()
  save2(pD, "D_lysosomal_vs_other_trends", 8, 5.5)
}

# ---- E: chain length x unsaturation ------------------------------------------------------------------------------------
eset <- x |> filter(class %in% big, tier != "Other", !is.na(total_DB))
if (nrow(eset)) {
  pE <- ggplot(eset, aes(total_C, total_DB, fill = lfc)) +
    geom_point(shape = 21, size = 3, colour = "white", stroke = 0.4) +
    facet_wrap(~class, scales = "free") +
    scale_fill_gradient2(low = "#2C5F9E", mid = "#F7F7F7", high = "#B2382C", midpoint = 0,
                         limits = c(-1.5, 1.5), oob = scales::squish, name = "log2FC") +
    labs(title = "Chain length and unsaturation together",
         subtitle = "each dot is one species; red = higher, blue = lower in IFNB. Separates a chain-length effect from a double-bond effect",
         x = "total acyl carbons", y = "total double bonds") + theme_f() +
    theme(panel.background = element_rect(fill = "grey93"), panel.grid = element_line(colour = "white"))
  nf <- n_distinct(eset$class)
  save2(pE, "E_chain_length_x_unsaturation", 3.3 * min(nf, 4) + 1, 3 * ceiling(nf / 4) + 1.6)
}
message("\nFigures written to: ", OUT)
