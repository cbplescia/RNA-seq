#!/usr/bin/env Rscript
# BMP (LBPA) species from the pilot lipidomics next to the transcripts of the enzymes that make,
# degrade or acyl-modify BMP, IFNB vs untreated/Mock.
#
# Lipid side  (01_lipidomics_ranked_IFNB_vs_Untreated.csv; LBPA / BMP class)
#   M1_BMP_relative_abundance      each species as % of the BMP signal, untreated vs IFNB
#   M2_BMP_log2FC_ranked           species ranked by log2FC, coloured by chain type (22:6, 18:1, other)
#   M3_BMP_composition_shift       share of the BMP pool carried by the main species, untreated vs IFNB
# RNA side  (vsd_v2.rds from 04 + DE_v2_IFNB_vs_Mock_protein_coding.csv)
#   R1_enzyme_transcript_abundance VST expression Mock vs IFNB, every sample a dot, one panel per enzyme
#   R2_enzyme_transcript_log2FC    log2FC (+/- lfcSE) grouped by function; coloured when padj < 0.05
#   N_summary                      lipid composition shift beside enzyme log2FC (needs the patchwork package)
# Relative abundance = a species' intensity / the sum over all BMP species in that condition (semi-quantitative:
# intensities of different species are not directly comparable). One lipid measurement per condition, so the
# lipid panels are descriptive. Enzyme transcript labels name the proposed role; edit PANEL below to change it.
# Usage: Rscript 38_bmp_lipids_vs_enzyme_transcripts.R [path/to/lipidomics.csv]

suppressPackageStartupMessages({ library(readr); library(dplyr); library(tidyr); library(ggplot2); library(DESeq2) })
select <- dplyr::select; first <- dplyr::first; rename <- dplyr::rename
slice <- dplyr::slice; filter <- dplyr::filter; desc <- dplyr::desc; count <- dplyr::count

PROJ <- Sys.getenv("PROJ", unset = file.path(
  Sys.getenv("HOME"),
  "Documents/aviti_rnaseq/2025/aviti",
  "20251009-100725_Maury_25161_WivfFacICOTHFAUJwSPeHBEFfCfNtCRpUSgVGUkz",
  "results/Maury_25161"))
RNA_OUT <- file.path(PROJ, "requant_v2_outputs")
OUT <- file.path(PROJ, "lipidomics_chain_length"); dir.create(OUT, showWarnings = FALSE)
INK <- "#1F2933"; UP <- "#C53030"; DOWN <- "#2B6CB0"; NS <- "#B8B8B8"
TOP_SPECIES <- 6

PANEL <- tibble::tribble(
  ~group,                          ~gene,     ~role,
  "BMP synthesis & turnover",      "Cln5",    "BMP synthase",
  "BMP synthesis & turnover",      "Pld3",    "lysosomal PLD",
  "BMP synthesis & turnover",      "Pld4",    "lysosomal PLD",
  "BMP synthesis & turnover",      "Pla2g15", "LPLA2 (BMP turnover)",
  "BMP synthesis & turnover",      "Pgs1",    "PG synthase (precursor)",
  "BMP synthesis & turnover",      "Cds1",    "CDP-DAG synthase",
  "BMP synthesis & turnover",      "Cds2",    "CDP-DAG synthase",
  "BMP synthesis & turnover",      "Ptpmt1",  "PGP phosphatase",
  "18:1 (monounsaturated) chains", "Scd1",    "stearoyl-CoA desaturase",
  "18:1 (monounsaturated) chains", "Scd2",    "stearoyl-CoA desaturase",
  "PUFA / DHA supply",             "Fads1",   "fatty acid desaturase",
  "PUFA / DHA supply",             "Fads2",   "fatty acid desaturase",
  "PUFA / DHA supply",             "Elovl5",  "PUFA elongase",
  "PUFA / DHA supply",             "Elovl6",  "C16 elongase",
  "PUFA / DHA supply",             "Acsl1",   "acyl-CoA synthetase",
  "PUFA / DHA supply",             "Acsl3",   "acyl-CoA synthetase",
  "PUFA / DHA supply",             "Acsl4",   "PUFA acyl-CoA synthetase",
  "Acyl remodeling",               "Lpcat3",  "PUFA into PC/PE",
  "Acyl remodeling",               "Agpat3",  "acyltransferase",
  "Acyl remodeling",               "Mboat1",  "acyltransferase",
  "Acyl remodeling",               "Mboat7",  "acyltransferase",
  "Acyl remodeling",               "Pla2g4a", "cPLA2 alpha",
  "Acyl remodeling",               "Pla2g6",  "iPLA2 beta")

theme_f <- function() theme_bw(base_size = 11) + theme(panel.grid.minor = element_blank(),
  panel.border = element_rect(colour = "grey70"), axis.text = element_text(colour = INK),
  strip.background = element_rect(fill = "grey95", colour = "grey70"), strip.text = element_text(colour = INK, face = "bold"),
  plot.title = element_text(colour = INK, face = "bold"), plot.subtitle = element_text(colour = "grey35", size = 9),
  legend.position = "bottom")
save2 <- function(p, name, w, h) {
  ggsave(file.path(OUT, paste0(name, ".png")), p, width = w, height = h, dpi = 300, limitsize = FALSE)
  ggsave(file.path(OUT, paste0(name, ".pdf")), p, width = w, height = h, limitsize = FALSE)
}

# =========================== lipids ======================================================================
a <- commandArgs(trailingOnly = TRUE)
fname <- if (length(a)) a[1] else "01_lipidomics_ranked_IFNB_vs_Untreated.csv"
cand <- c(fname, file.path(PROJ, fname), file.path(getwd(), fname)); rf <- cand[file.exists(cand)][1]
if (is.na(rf)) stop("Cannot find ", fname, " - give its full path as the first argument")
raw <- read_csv(rf, show_col_types = FALSE)
miss <- (if ("any_missing_Untreated" %in% names(raw)) raw$any_missing_Untreated else FALSE) |
        (if ("any_missing_IFNB" %in% names(raw)) raw$any_missing_IFNB else FALSE)
b <- raw[raw$Class %in% c("LBPA", "BMP") & !grepl("\\+D[0-9]|-D[0-9]|\\(s\\)", raw$Lipid_ID) & !miss, ] |>
  transmute(lipid = Lipid_ID, untreated = Untreated_Norm, ifnb = B_Norm, lfc = log2FC_IFNB_vs_Untreated)
if (!nrow(b)) stop("No LBPA/BMP species found in the lipidomics file")
message("BMP/LBPA species: ", nrow(b))

chains <- regmatches(b$lipid, gregexpr("(?<![0-9])[0-9]{1,2}:[0-9]{1,2}", b$lipid, perl = TRUE))
sorted <- lapply(chains, function(z) { p <- do.call(rbind, strsplit(z, ":", fixed = TRUE)); z[order(as.numeric(p[, 1]), as.numeric(p[, 2]))] })
b$label <- make.unique(vapply(sorted, paste, character(1), collapse = "/"))
b$dha <- vapply(chains, function(z) "22:6" %in% z, logical(1))
b$c181 <- vapply(chains, function(z) "18:1" %in% z, logical(1))
b$type <- case_when(b$dha & b$c181 ~ "18:1 + 22:6", b$dha ~ "contains 22:6 (DHA)",
                    b$c181 ~ "contains 18:1 (no 22:6)", TRUE ~ "other chains")
b$pct_untreated <- 100 * b$untreated / sum(b$untreated)
b$pct_ifnb <- 100 * b$ifnb / sum(b$ifnb)
write_csv(b, file.path(OUT, "bmp_species_relative_abundance.csv"))
type_cols <- c(`contains 22:6 (DHA)` = "#D55E00", `18:1 + 22:6` = "#CC79A7",
               `contains 18:1 (no 22:6)` = "#0072B2", `other chains` = "#9CA3AF")
tshare <- b |> group_by(type) |> summarise(untreated = sum(pct_untreated), ifnb = sum(pct_ifnb), .groups = "drop")
message("\nShare of the BMP signal by chain type (%):"); print(as.data.frame(mutate(tshare, across(where(is.numeric), ~ round(.x, 1)))), row.names = FALSE)

# M1: paired relative abundance
m1 <- b |> arrange(desc(pct_untreated + pct_ifnb)) |>
  pivot_longer(c(pct_untreated, pct_ifnb), names_to = "cond", values_to = "pct") |>
  mutate(cond = factor(ifelse(cond == "pct_untreated", "untreated", "IFNB"), levels = c("untreated", "IFNB")))
m1$label <- factor(m1$label, levels = rev(unique(m1$label)))
pM1 <- ggplot(m1, aes(pct, label, fill = cond)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.75) +
  scale_fill_manual(values = c(untreated = "#8A8A8A", IFNB = DOWN), name = NULL) +
  labs(title = "BMP (LBPA) species: relative abundance", x = "% of total BMP signal in that condition", y = NULL,
       subtitle = "label = the two acyl chains; ordered by combined abundance; intensities are semi-quantitative") + theme_f()
save2(pM1, "M1_BMP_relative_abundance", 7, 0.3 * nrow(b) + 2.2)

# M2: ranked log2FC
m2 <- b |> arrange(lfc); m2$label <- factor(m2$label, levels = m2$label)
pM2 <- ggplot(m2, aes(lfc, label, fill = type)) +
  geom_vline(xintercept = 0, colour = "grey40", linewidth = 0.3) + geom_col(width = 0.75) +
  scale_fill_manual(values = type_cols, name = NULL) +
  labs(title = "BMP (LBPA) species ranked by log2FC", x = "log2 fold-change (IFNB vs untreated)", y = NULL,
       subtitle = "colour = chain type; DHA = 22:6") + theme_f() + guides(fill = guide_legend(nrow = 2))
save2(pM2, "M2_BMP_log2FC_ranked", 7, 0.26 * nrow(b) + 2.4)

# M3: composition shift
top <- b |> arrange(desc(pct_untreated)) |> head(TOP_SPECIES) |> pull(label)
m3 <- b |> mutate(grp = ifelse(label %in% top, label, "all other species")) |>
  group_by(grp) |> summarise(untreated = sum(pct_untreated), IFNB = sum(pct_ifnb), .groups = "drop") |>
  pivot_longer(c(untreated, IFNB), names_to = "cond", values_to = "pct") |>
  mutate(cond = factor(cond, levels = c("untreated", "IFNB")),
         grp = factor(grp, levels = c("all other species", rev(top))))
pal <- c("#9CA3AF", "#56B4E9", "#E69F00", "#009E73", "#CC79A7", "#0072B2", "#D55E00")
names(pal) <- c("all other species", rev(top))[seq_along(pal)]
pM3 <- ggplot(m3, aes(cond, pct, fill = grp)) +
  geom_col(width = 0.7, colour = "white", linewidth = 0.4) +
  geom_text(aes(label = ifelse(pct >= 4, sprintf("%.0f%%", pct), "")), position = position_stack(vjust = 0.5), size = 3, colour = INK) +
  scale_fill_manual(values = pal, name = NULL) +
  labs(title = "How the BMP pool is composed", x = NULL, y = "% of total BMP signal",
       subtitle = sprintf("the %d most abundant species and the rest; labels shown for slices >= 4%%", TOP_SPECIES)) + theme_f()
save2(pM3, "M3_BMP_composition_shift", 6, 6)

# =========================== transcripts ==================================================================
de <- read_csv(file.path(RNA_OUT, "DE_v2_IFNB_vs_Mock_protein_coding.csv"), show_col_types = FALSE)
pn <- PANEL |> left_join(de |> select(gene = gene_name, gene_id, baseMean, log2FoldChange, lfcSE, padj), by = "gene")
if (any(is.na(pn$gene_id))) message("\nNot in the RNA-seq data (left out): ", paste(pn$gene[is.na(pn$gene_id)], collapse = ", "))
pn <- pn |> filter(!is.na(gene_id)) |>
  mutate(group = factor(group, levels = unique(PANEL$group)), gene = factor(gene, levels = unique(PANEL$gene)),
         direction = case_when(!is.na(padj) & padj < 0.05 & log2FoldChange > 0 ~ "Up", !is.na(padj) & padj < 0.05 & log2FoldChange < 0 ~ "Down",
                               TRUE ~ "Not significant"),
         stars = case_when(is.na(padj) ~ "", padj < 0.001 ~ "***", padj < 0.01 ~ "**", padj < 0.05 ~ "*", TRUE ~ "ns"))
write_csv(pn |> select(group, gene, role, baseMean, log2FoldChange, lfcSE, padj, stars), file.path(OUT, "bmp_enzyme_transcripts.csv"))

# R2: log2FC by group
pn$gene_f <- factor(pn$gene, levels = rev(levels(pn$gene)))
pR2 <- ggplot(pn, aes(log2FoldChange, gene_f, fill = direction)) +
  geom_vline(xintercept = 0, colour = "grey40", linewidth = 0.3) +
  geom_col(width = 0.7) +
  geom_errorbar(aes(xmin = log2FoldChange - lfcSE, xmax = log2FoldChange + lfcSE), width = 0.25, colour = INK, linewidth = 0.3, orientation = "y") +
  geom_text(aes(x = log2FoldChange + sign(log2FoldChange) * (lfcSE + 0.04), label = stars), size = 3, colour = INK,
            hjust = ifelse(pn$log2FoldChange >= 0, 0, 1)) +
  facet_grid(group ~ ., scales = "free_y", space = "free_y") +
  scale_fill_manual(values = c(Up = UP, Down = DOWN, `Not significant` = NS), name = NULL) +
  scale_x_continuous(expand = expansion(mult = 0.12)) +
  labs(title = "Transcripts of BMP-related enzymes: IFNB vs Mock", x = "log2 fold-change (shrunk, +/- lfcSE); * padj<0.05, ** <0.01, *** <0.001", y = NULL,
       subtitle = "mRNA level, not enzyme activity") + theme_f() + theme(axis.text.y = element_text(face = "italic"), strip.text.y = element_text(angle = 0, hjust = 0))
save2(pR2, "R2_enzyme_transcript_log2FC", 7.5, 0.26 * nrow(pn) + 3)

# R1: VST abundance
vf <- file.path(RNA_OUT, "vsd_v2.rds")
if (file.exists(vf)) {
  vsd <- readRDS(vf); vm <- assay(vsd)
  smp <- tibble(sample_id = colnames(vm), group = as.character(vsd$group)) |> filter(group %in% c("Mock_Mock", "IFNB_Mock")) |>
    mutate(treatment = factor(ifelse(group == "Mock_Mock", "Mock", "IFNB"), levels = c("Mock", "IFNB")))
  ids <- pn$gene_id[pn$gene_id %in% rownames(vm)]
  lg <- as.data.frame(t(vm[ids, smp$sample_id, drop = FALSE])) |> tibble::rownames_to_column("sample_id") |>
    pivot_longer(-sample_id, names_to = "gene_id", values_to = "vst") |>
    left_join(smp, by = "sample_id") |> left_join(pn |> select(gene_id, gene, role, stars), by = "gene_id") |>
    mutate(panel = factor(sprintf("%s  %s\n%s", gene, stars, role), levels = sprintf("%s  %s\n%s", pn$gene, pn$stars, pn$role)))
  write_csv(lg |> select(-panel), file.path(OUT, "bmp_enzyme_vst_values.csv"))
  pR1 <- ggplot(lg, aes(treatment, vst, fill = treatment, colour = treatment)) +
    geom_boxplot(alpha = 0.25, outlier.shape = NA, width = 0.55, linewidth = 0.5) +
    geom_jitter(width = 0.1, height = 0, size = 1.6, alpha = 0.9) +
    facet_wrap(~panel, ncol = 6, scales = "free_y") +
    scale_fill_manual(values = c(Mock = "#8A8A8A", IFNB = DOWN), guide = "none") +
    scale_colour_manual(values = c(Mock = "#8A8A8A", IFNB = DOWN), guide = "none") +
    labs(title = "Transcript abundance of BMP-related enzymes: Mock vs IFNB",
         subtitle = "VST expression (log2 scale), one dot per sample (n = 3); stars = padj of the IFNB vs Mock test", x = NULL, y = "VST expression") + theme_f() +
    theme(strip.text = element_text(size = 8))
  save2(pR1, "R1_enzyme_transcript_abundance", 13, 2.9 * ceiling(nlevels(lg$panel) / 6) + 1.2)
} else message("vsd_v2.rds not found - R1 skipped")

# combined summary
if (requireNamespace("patchwork", quietly = TRUE)) {
  save2(patchwork::wrap_plots(pM3 + theme(legend.position = "right"), pR2, widths = c(1, 1.4)) +
          patchwork::plot_annotation(title = "BMP pool composition and the transcripts of its enzymes"),
        "N_summary", 14, 7)
} else message("patchwork not installed - skipping the combined N_summary figure (conda install -c conda-forge r-patchwork)")
message("\nFigures in: ", OUT)
