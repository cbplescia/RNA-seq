#!/usr/bin/env Rscript
# PCA of Mock, IFNB, EBOV and IFNB+EBOV samples (Mock 1-3, IFNB 10-12, EBOV 4-6, IFNB+EBOV 16-18) on three gene sets:
#   1. ALL protein-coding genes
#   2. ALL IFNB ISGs   (IFNB vs Mock: padj < 0.05, log2FC > 1)
#   3. ALL IFNB IDGs   (IFNB vs Mock: padj < 0.05, log2FC < -1)
# No gene filtering beyond the set definition (every gene in the set with non-zero variance is used; NOT just the most variable).
# Input = VST expression. PCA centres each gene across the samples, so using "log2 change vs Mock" instead would give exactly the
# same PCA (subtracting each gene's Mock mean only shifts the gene by a constant). Set SCALE <- TRUE to give every gene equal weight.
# Solid circles coloured by group (Mock grey, IFNB blue, EBOV dark green, IFNB+EBOV light green), replicate number inside the dot.
# Writes to <project>/requant_v2_outputs/figures_mock_ifnb_ebov/pca_gene_sets/ (PNG + PDF + CSV):
#   PCA_<set>.png/.pdf, PCA_<set>_coordinates.csv, PCA_<set>_top_loadings.csv, PCA_gene_sets_combined.png (needs patchwork),
#   PCA_gene_sets_variance.csv (percent variance of PC1-PC4 for each set)
# Needs DE_v2_IFNB_vs_Mock_protein_coding.csv and vsd_v2.rds. Usage: Rscript 53_pca_gene_sets.R   (PROJ from 'source setup_session.sh')

suppressPackageStartupMessages({ library(DESeq2); library(readr); library(dplyr); library(tidyr); library(ggplot2) })
select <- dplyr::select; first <- dplyr::first; rename <- dplyr::rename
slice <- dplyr::slice; filter <- dplyr::filter; desc <- dplyr::desc; count <- dplyr::count

LFC_CUT <- 1; PADJ_CUT <- 0.05; SCALE <- FALSE
SETS <- list(Mock = c("1", "2", "3"), IFNB = c("10", "11", "12"), `IFNB+EBOV` = c("16", "17", "18"), EBOV = c("4", "5", "6"))
EXPECTED <- c(Mock = "Mock_Mock", IFNB = "IFNB_Mock", `IFNB+EBOV` = "IFNB_EBOV", EBOV = "Mock_EBOV")
GCOL <- c(Mock = "#7F7F7F", IFNB = "#2B6CB0", `IFNB+EBOV` = "#9BD08E", EBOV = "#1B7837")
INK <- "#1F2933"

PROJ <- Sys.getenv("PROJ", unset = file.path(
  Sys.getenv("HOME"),
  "Documents/aviti_rnaseq/2025/aviti",
  "20251009-100725_Maury_25161_WivfFacICOTHFAUJwSPeHBEFfCfNtCRpUSgVGUkz",
  "results/Maury_25161"))
RNA <- file.path(PROJ, "requant_v2_outputs")
OUT <- file.path(RNA, "figures_mock_ifnb_ebov", "pca_gene_sets"); dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

vsd <- readRDS(file.path(RNA, "vsd_v2.rds")); vm <- assay(vsd)
sel <- tibble(sample = unlist(SETS), group = rep(names(SETS), lengths(SETS)), replicate = unlist(lapply(SETS, seq_along)))
sel$data_group <- as.character(vsd$group)[match(sel$sample, colnames(vm))]
if (any(is.na(sel$data_group))) stop("Samples not in the VST object: ", paste(sel$sample[is.na(sel$data_group)], collapse = ", "))
message("Check - chosen group vs group recorded in the data:")
print(as.data.frame(sel |> transmute(sample, chosen_as = group, replicate, group_in_data = data_group,
                                     match = ifelse(data_group == EXPECTED[group], "ok", "DIFFERENT"))), row.names = FALSE)
if (any(sel$data_group != EXPECTED[sel$group])) message("WARNING: some samples are recorded under a different group. Using YOUR sample numbers.")

de <- read_csv(file.path(RNA, "DE_v2_IFNB_vs_Mock_protein_coding.csv"), show_col_types = FALSE) |>
  filter(gene_id %in% rownames(vm), !is.na(padj))
GENESETS <- list(
  all_protein_coding = list(title = "All protein-coding genes", ids = de$gene_id),
  ISGs = list(title = "IFNB ISGs (log2FC > 1)", ids = de$gene_id[de$padj < PADJ_CUT & de$log2FoldChange > LFC_CUT]),
  IDGs = list(title = "IFNB IDGs (log2FC < -1)", ids = de$gene_id[de$padj < PADJ_CUT & de$log2FoldChange < -LFC_CUT]))
sym <- setNames(de$gene_name, de$gene_id)
X <- vm[, sel$sample, drop = FALSE]

plots <- list(); var_rows <- list()
for (id in names(GENESETS)) {
  gs <- GENESETS[[id]]
  x <- X[intersect(gs$ids, rownames(X)), , drop = FALSE]
  x <- x[apply(x, 1, var) > 0, , drop = FALSE]
  p <- prcomp(t(x), center = TRUE, scale. = SCALE)
  pct <- 100 * p$sdev^2 / sum(p$sdev^2)
  message(sprintf("\n%s: %d genes, 12 samples; variance explained PC1 %.1f%%, PC2 %.1f%%, PC3 %.1f%%", gs$title, nrow(x), pct[1], pct[2], pct[3]))
  var_rows[[id]] <- tibble(gene_set = gs$title, genes = nrow(x), PC1 = pct[1], PC2 = pct[2], PC3 = pct[3], PC4 = pct[4])
  co <- as.data.frame(p$x[, 1:4]) |> tibble::rownames_to_column("sample") |> left_join(sel |> select(sample, group, replicate), by = "sample") |>
    mutate(group = factor(group, levels = names(SETS)), txt = ifelse(group == "IFNB+EBOV", INK, "white"))
  write_csv(co, file.path(OUT, paste0("PCA_", id, "_coordinates.csv")))
  ld <- as.data.frame(p$rotation[, 1:2]) |> tibble::rownames_to_column("gene_id") |> mutate(gene_name = unname(sym[gene_id]))
  write_csv(bind_rows(ld |> slice_max(abs(PC1), n = 20) |> mutate(driver_of = "PC1"), ld |> slice_max(abs(PC2), n = 20) |> mutate(driver_of = "PC2")),
            file.path(OUT, paste0("PCA_", id, "_top_loadings.csv")))
  hull <- co |> group_by(group) |> slice(chull(PC1, PC2)) |> ungroup()
  g <- ggplot(co, aes(PC1, PC2)) +
    geom_hline(yintercept = 0, colour = "grey85", linewidth = 0.3) + geom_vline(xintercept = 0, colour = "grey85", linewidth = 0.3) +
    geom_polygon(data = hull, aes(fill = group, group = group), alpha = 0.13, colour = NA) +
    geom_point(aes(fill = group), shape = 21, size = 8, colour = "white", stroke = 0.9) +
    geom_text(aes(label = replicate, colour = txt), size = 3.6, fontface = "bold") +
    scale_colour_identity() + scale_fill_manual(values = GCOL, name = "group") +
    guides(fill = guide_legend(override.aes = list(shape = 21, size = 5, alpha = 1))) +
    scale_x_continuous(expand = expansion(mult = 0.12)) + scale_y_continuous(expand = expansion(mult = 0.12)) +
    labs(title = paste0("PCA: ", gs$title), subtitle = sprintf("%d genes (VST), 12 samples; number inside each dot = replicate%s", nrow(x), if (SCALE) "; genes scaled" else ""),
         x = sprintf("PC1 (%.1f%% of variance)", pct[1]), y = sprintf("PC2 (%.1f%% of variance)", pct[2])) +
    theme_bw(base_size = 12) + theme(panel.grid.minor = element_blank(), panel.border = element_rect(colour = "grey70"),
          axis.text = element_text(colour = INK), plot.title = element_text(colour = INK, face = "bold"),
          plot.subtitle = element_text(colour = "grey35", size = 9), legend.position = "right")
  ggsave(file.path(OUT, paste0("PCA_", id, ".png")), g, width = 8, height = 6.5, dpi = 300)
  ggsave(file.path(OUT, paste0("PCA_", id, ".pdf")), g, width = 8, height = 6.5)
  plots[[id]] <- g
}
write_csv(bind_rows(var_rows), file.path(OUT, "PCA_gene_sets_variance.csv"))
if (requireNamespace("patchwork", quietly = TRUE)) {
  comb <- patchwork::wrap_plots(lapply(plots, function(p) p + theme(legend.position = "none", plot.subtitle = element_text(size = 8))), nrow = 1) +
    patchwork::plot_annotation(title = "PCA by gene set: Mock, IFNB, EBOV and IFNB+EBOV")
  ggsave(file.path(OUT, "PCA_gene_sets_combined.png"), comb, width = 18, height = 6.2, dpi = 300)
} else message("patchwork not installed - skipping the combined 3-panel figure (conda install -c conda-forge r-patchwork)")
message("\nFigures and tables in: ", OUT)
