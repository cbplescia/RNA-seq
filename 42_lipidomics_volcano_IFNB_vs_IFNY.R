#!/usr/bin/env Rscript
# Differential-abundance plot of the pilot lipidomics, IFNB vs IFNy:
#   lipids MORE abundant in IFNB are dots on the RIGHT (log2FC > 0), MORE abundant in IFNy on the LEFT (< 0).
# Uses the file's own log2FC_IFNB_vs_IFNY column (positive = higher in IFNB; checked against B_Norm / Y_Norm).
# The top 20 on each side are labelled.
#
# The file has no p-values or replicates, so a classic volcano (-log10 p) is not possible: the y axis is the MEAN
# ABUNDANCE (log10 of the mean of the IFNB and IFNy normalised intensities), an MA-style plot. Dashed lines at
# log2FC = +/-1 are a guide, not a significance cut-off.
#
# Filters (switch below):
#   EXCLUDE_STANDARDS  internal standards (+D7:(s), -D7 ...)                                      default TRUE
#   EXCLUDE_MISSING    rows flagged with missing values in IFNB or IFNy                           default TRUE
#   MIN_GROUP_INTENSITY lipid must reach this normalised intensity in BOTH IFNB and IFNy; lipids that do not are
#                      "detected in one condition only" - their log2FC is set by a near-zero value (e.g. 0-26 vs
#                      millions gives +18) - and are written to lipidomics_IFNB_vs_IFNY_one_sided.csv instead   default 10000 (0 = off)
#   EXCLUDE_OXIDIZED   oxidised species (COOH / OOH); FALSE keeps them as triangles                  default FALSE
#
# Usage: Rscript 42_lipidomics_volcano_IFNB_vs_IFNY.R [path/to/01_lipidomics_ranked_IFNB_vs_Untreated.csv]
# Writes to <project>/lipidomics_chain_length/: lipidomics_volcano_IFNB_vs_IFNY.png/.pdf,
# lipidomics_volcano_IFNB_vs_IFNY_labelled.csv, lipidomics_IFNB_vs_IFNY_one_sided.csv

suppressPackageStartupMessages({ library(readr); library(dplyr); library(ggplot2) })
if (!requireNamespace("ggrepel", quietly = TRUE)) stop("Install ggrepel (conda install -c conda-forge r-ggrepel)")

EXCLUDE_STANDARDS <- TRUE; EXCLUDE_MISSING <- TRUE; EXCLUDE_OXIDIZED <- FALSE; MIN_GROUP_INTENSITY <- 10000
N_LABEL <- 20; LFC_GUIDE <- 1
INK <- "#1F2933"; B_COL <- "#2B6CB0"; Y_COL <- "#D69E2E"; NS <- "#B8B8B8"

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
need <- c("Lipid_ID", "Class", "B_Norm", "Y_Norm", "log2FC_IFNB_vs_IFNY")
if (!all(need %in% names(raw))) stop("Missing columns: ", paste(setdiff(need, names(raw)), collapse = ", "))

d <- raw |> transmute(lipid = Lipid_ID, class = Class, lfc = log2FC_IFNB_vs_IFNY, ifnb = B_Norm, ifny = Y_Norm,
                      mean_intensity = (B_Norm + Y_Norm) / 2,
                      missing = (if ("any_missing_IFNB" %in% names(raw)) any_missing_IFNB else FALSE) |
                                (if ("any_missing_IFNY" %in% names(raw)) any_missing_IFNY else FALSE))
d$standard <- grepl("\\+D[0-9]|-D[0-9]|\\(s\\)", d$lipid)
d$oxidized <- grepl("COOH|OOH", d$lipid)
message(sprintf("Lipids: %d  (standards %d, missing-value rows %d, oxidised %d)", nrow(d), sum(d$standard), sum(d$missing), sum(d$oxidized)))
ok <- !is.na(d$lfc) & is.finite(d$lfc)
if (EXCLUDE_STANDARDS) ok <- ok & !d$standard
if (EXCLUDE_MISSING) ok <- ok & !d$missing
if (EXCLUDE_OXIDIZED) ok <- ok & !d$oxidized
d <- d[ok, ]
both <- d$ifnb >= MIN_GROUP_INTENSITY & d$ifny >= MIN_GROUP_INTENSITY
one_sided <- d[!both, ] |> arrange(desc(lfc))
write_csv(one_sided |> select(lipid, class, lfc, ifnb, ifny), file.path(OUT, "lipidomics_IFNB_vs_IFNY_one_sided.csv"))
d <- d[both, ]
message(sprintf("Below the detection floor (%d) in IFNB or IFNy, saved separately: %d;  plotted: %d", MIN_GROUP_INTENSITY, nrow(one_sided), nrow(d)))

d$y <- log10(d$mean_intensity)
d$status <- ifelse(d$lfc >= LFC_GUIDE, "higher in IFNB", ifelse(d$lfc <= -LFC_GUIDE, "higher in IFNγ", "within +/-1"))
lab <- bind_rows(d |> filter(lfc > 0) |> slice_max(lfc, n = N_LABEL, with_ties = FALSE),
                 d |> filter(lfc < 0) |> slice_min(lfc, n = N_LABEL, with_ties = FALSE))
write_csv(lab |> mutate(higher_in = ifelse(lfc > 0, "IFNB", "IFNy")) |> arrange(desc(lfc)) |>
            select(lipid, class, higher_in, lfc, ifnb, ifny, oxidized), file.path(OUT, "lipidomics_volcano_IFNB_vs_IFNY_labelled.csv"))
message("Top ", N_LABEL, " higher in IFNB (right): ", paste(lab$lipid[lab$lfc > 0], collapse = ", "))
message("Top ", N_LABEL, " higher in IFNy (left): ", paste(lab$lipid[lab$lfc < 0], collapse = ", "))
message(sprintf("Beyond +/-1: %d higher in IFNB, %d higher in IFNy (of %d plotted)", sum(d$lfc >= 1), sum(d$lfc <= -1), nrow(d)))

xr <- max(abs(d$lfc)) * 1.08
p <- ggplot(d, aes(lfc, y)) +
  geom_vline(xintercept = c(-LFC_GUIDE, LFC_GUIDE), linetype = "dashed", colour = "grey60", linewidth = 0.3) +
  geom_vline(xintercept = 0, colour = "grey75", linewidth = 0.3) +
  geom_point(aes(colour = status, shape = oxidized), size = 1.8, alpha = 0.7) +
  geom_point(data = lab, aes(colour = status, shape = oxidized), size = 2.4) +
  ggrepel::geom_text_repel(data = lab, aes(label = lipid), colour = INK, size = 2.7, max.overlaps = Inf,
                           min.segment.length = 0, segment.colour = "grey60", segment.size = 0.25, box.padding = 0.35, seed = 1) +
  scale_colour_manual(values = c(`higher in IFNB` = B_COL, `higher in IFNγ` = Y_COL, `within +/-1` = NS), name = NULL,
                      breaks = c("higher in IFNγ", "within +/-1", "higher in IFNB")) +
  scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 17), labels = c(`FALSE` = "lipid", `TRUE` = "oxidised species (COOH/OOH)"), name = NULL) +
  scale_x_continuous(limits = c(-xr, xr)) +
  labs(title = "Pilot lipidomics: IFNB vs IFNγ",
       subtitle = sprintf("each point is one lipid; top %d on each side labelled; lipids below the detection floor (%d) in either condition left out; no p-values, so y = abundance",
                          N_LABEL, MIN_GROUP_INTENSITY),
       x = "log2 fold-change (IFNB / IFNγ)     ← more abundant in IFNγ  |  more abundant in IFNB →",
       y = "log10 mean normalised intensity") +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_line(colour = "grey94", linewidth = 0.3),
        panel.border = element_rect(colour = "grey70"), axis.text = element_text(colour = INK), axis.title = element_text(colour = INK),
        plot.title = element_text(colour = INK, face = "bold"), plot.subtitle = element_text(colour = "grey35", size = 8.5),
        legend.position = "bottom") + guides(colour = guide_legend(nrow = 1), shape = guide_legend(nrow = 1))
ggsave(file.path(OUT, "lipidomics_volcano_IFNB_vs_IFNY.png"), p, width = 10, height = 8, dpi = 300)
ggsave(file.path(OUT, "lipidomics_volcano_IFNB_vs_IFNY.pdf"), p, width = 10, height = 8)
message("\nSaved to: ", OUT)
