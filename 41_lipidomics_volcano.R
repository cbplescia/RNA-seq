#!/usr/bin/env Rscript
# Differential-abundance plot of the pilot lipidomics, IFNB vs untreated, with the top 20 increased and top 20
# decreased lipids labelled.
#
# The file has no p-values or replicate-level data, so a classic volcano (-log10 p) is not possible. The y axis is the
# MEAN ABUNDANCE of each lipid instead (log10 of the mean of the untreated and IFNB normalised intensities), i.e. an
# MA-style plot: it shows whether a big change is in a well-measured lipid (high) or a barely detected one (low).
# Dashed lines mark |log2FC| = 1 as a guide only; they are not a significance cut-off.
#
# Filters (switch below):
#   EXCLUDE_STANDARDS  drop internal standards (names with +D7:(s), -D7 ...)                         default TRUE
#   EXCLUDE_MISSING    drop lipids with any missing value in untreated / IFNB (log2FC inflated by 0) default TRUE
#   EXCLUDE_OXIDIZED   drop oxidised species (COOH / OOH in the name); FALSE keeps them as triangles  default FALSE
#   MIN_MEAN_INTENSITY drop lipids whose mean intensity is below this (0 = keep all)                  default 0
# Top 20 up / down = largest positive / most negative log2FC among the kept lipids.
#
# Usage: Rscript 41_lipidomics_volcano.R [path/to/01_lipidomics_ranked_IFNB_vs_Untreated.csv]
# Writes to <project>/lipidomics_chain_length/: lipidomics_volcano_IFNB_vs_untreated.png/.pdf and
# lipidomics_volcano_labelled_lipids.csv

suppressPackageStartupMessages({ library(readr); library(dplyr); library(ggplot2) })
if (!requireNamespace("ggrepel", quietly = TRUE)) stop("Install ggrepel (conda install -c conda-forge r-ggrepel)")

EXCLUDE_STANDARDS <- TRUE; EXCLUDE_MISSING <- TRUE; EXCLUDE_OXIDIZED <- FALSE; MIN_MEAN_INTENSITY <- 0
N_LABEL <- 20; LFC_GUIDE <- 1
INK <- "#1F2933"; UP <- "#C53030"; DOWN <- "#2B6CB0"; NS <- "#B8B8B8"

PROJ <- Sys.getenv("PROJ", unset = file.path(
  Sys.getenv("HOME"),
  "Documents/aviti_rnaseq/2025/aviti",
  "20251009-100725_Maury_25161_WivfFacICOTHFAUJwSPeHBEFfCfNtCRpUSgVGUkz",
  "results/Maury_25161"))
OUT <- file.path(PROJ, "lipidomics_chain_length"); dir.create(OUT, showWarnings = FALSE)

a <- commandArgs(trailingOnly = TRUE)
fname <- if (length(a)) a[1] else "01_lipidomics_ranked_IFNB_vs_Untreated.csv"
cand <- c(fname, file.path(PROJ, fname), file.path(getwd(), fname)); f <- cand[file.exists(cand)][1]
if (is.na(f)) stop("Cannot find ", fname, " - give its full path as the first argument")
raw <- read_csv(f, show_col_types = FALSE)
need <- c("Lipid_ID", "Class", "Untreated_Norm", "B_Norm", "log2FC_IFNB_vs_Untreated")
if (!all(need %in% names(raw))) stop("Missing columns: ", paste(setdiff(need, names(raw)), collapse = ", "))

d <- raw |> transmute(lipid = Lipid_ID, class = Class, lfc = log2FC_IFNB_vs_Untreated,
                      mean_intensity = (Untreated_Norm + B_Norm) / 2,
                      missing = (if ("any_missing_Untreated" %in% names(raw)) any_missing_Untreated else FALSE) |
                                (if ("any_missing_IFNB" %in% names(raw)) any_missing_IFNB else FALSE))
d$standard <- grepl("\\+D[0-9]|-D[0-9]|\\(s\\)", d$lipid)
d$oxidized <- grepl("COOH|OOH", d$lipid)
message(sprintf("Lipids: %d  (standards %d, missing-value rows %d, oxidised %d)", nrow(d), sum(d$standard), sum(d$missing), sum(d$oxidized)))
keep <- !is.na(d$lfc) & is.finite(d$lfc) & d$mean_intensity > 0
if (EXCLUDE_STANDARDS) keep <- keep & !d$standard
if (EXCLUDE_MISSING) keep <- keep & !d$missing
if (EXCLUDE_OXIDIZED) keep <- keep & !d$oxidized
keep <- keep & d$mean_intensity >= MIN_MEAN_INTENSITY
d <- d[keep, ]
message("Plotted: ", nrow(d), " lipids")
d$y <- log10(d$mean_intensity)
d$status <- ifelse(d$lfc >= LFC_GUIDE, "increased (log2FC >= 1)", ifelse(d$lfc <= -LFC_GUIDE, "decreased (log2FC <= -1)", "within +/-1"))
lab <- bind_rows(d |> filter(lfc > 0) |> slice_max(lfc, n = N_LABEL, with_ties = FALSE),
                 d |> filter(lfc < 0) |> slice_min(lfc, n = N_LABEL, with_ties = FALSE))
write_csv(lab |> mutate(direction = ifelse(lfc > 0, "up", "down")) |> arrange(desc(lfc)) |>
            select(lipid, class, direction, lfc, mean_intensity, oxidized), file.path(OUT, "lipidomics_volcano_labelled_lipids.csv"))
message("Top ", N_LABEL, " increased: ", paste(lab$lipid[lab$lfc > 0], collapse = ", "))
message("Top ", N_LABEL, " decreased: ", paste(lab$lipid[lab$lfc < 0], collapse = ", "))

xr <- max(abs(d$lfc)) * 1.08
p <- ggplot(d, aes(lfc, y)) +
  geom_vline(xintercept = c(-LFC_GUIDE, LFC_GUIDE), linetype = "dashed", colour = "grey60", linewidth = 0.3) +
  geom_vline(xintercept = 0, colour = "grey75", linewidth = 0.3) +
  geom_point(aes(colour = status, shape = oxidized), size = 1.8, alpha = 0.7) +
  geom_point(data = lab, aes(colour = status, shape = oxidized), size = 2.4) +
  ggrepel::geom_text_repel(data = lab, aes(label = lipid), colour = INK, size = 2.7, max.overlaps = Inf,
                           min.segment.length = 0, segment.colour = "grey60", segment.size = 0.25,
                           box.padding = 0.35, seed = 1) +
  scale_colour_manual(values = c(`increased (log2FC >= 1)` = UP, `decreased (log2FC <= -1)` = DOWN, `within +/-1` = NS), name = NULL) +
  scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 17), labels = c(`FALSE` = "lipid", `TRUE` = "oxidised species (COOH/OOH)"), name = NULL) +
  scale_x_continuous(limits = c(-xr, xr)) +
  labs(title = "Pilot lipidomics: IFNB vs untreated",
       subtitle = sprintf("each point is one lipid; top %d increased and top %d decreased by log2FC labelled; no p-values in the data, so y = abundance",
                          N_LABEL, N_LABEL),
       x = "log2 fold-change (IFNB vs untreated)", y = "log10 mean normalised intensity") +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_line(colour = "grey94", linewidth = 0.3),
        panel.border = element_rect(colour = "grey70"), axis.text = element_text(colour = INK), axis.title = element_text(colour = INK),
        plot.title = element_text(colour = INK, face = "bold"), plot.subtitle = element_text(colour = "grey35", size = 9),
        legend.position = "bottom") + guides(colour = guide_legend(nrow = 1), shape = guide_legend(nrow = 1))
ggsave(file.path(OUT, "lipidomics_volcano_IFNB_vs_untreated.png"), p, width = 10, height = 8, dpi = 300)
ggsave(file.path(OUT, "lipidomics_volcano_IFNB_vs_untreated.pdf"), p, width = 10, height = 8)
message("\nSaved to: ", OUT)
