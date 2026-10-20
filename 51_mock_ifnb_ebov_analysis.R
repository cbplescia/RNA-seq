#!/usr/bin/env Rscript
# Mock / IFNB / EBOV / IFNB+EBOV only.  Samples: Mock 1-3, IFNB 10-12, EBOV 4-6, IFNB+EBOV 16-18.
# Heatmap column order, left to right: Mock | IFNB | IFNB+EBOV | EBOV.
#
#   1. PCA_Mock_IFNB_EBOV_IFNBEBOV   PCA of the 12 samples (500 most variable protein-coding genes, VST); solid circles coloured by
#                                    group, replicate number (1-3) inside each dot
#   2. heatmap_all_protein_coding_organised   every protein-coding gene, rows organised into labelled blocks:
#        IFNB-UP genes      -> EBOV suppresses / partly suppresses / does not suppress the induction
#        IFNB-DOWN genes    -> EBOV reverses / partly reverses / does not reverse the repression
#        not IFNB-responsive -> altered by EBOV alone / unchanged
#      Two colour bars on the rows give the IFNB response and the EBOV effect; within a block genes are ordered by IFNB log2FC.
#   3. heatmap_top30_ISGs_4conditions   the 30 most upregulated IFNB-vs-Mock genes (padj < 0.05, log2FC > 1), ranked, gene names shown
#   4. heatmap_top50_ISGs_Mock_IFNB (+ a top-30 version)   the same ISGs, Mock and IFNB columns only
# Colour = log2 change of each sample's VST vs the mean of the Mock samples (blue lower, red higher), capped at the 98th percentile.
# EBOV effect on the IFNB response (for IFNB-responsive genes: |log2FC| > LFC_CUT, padj < PADJ_CUT):
#   retention = log2FC(IFNB+EBOV vs Mock) / log2FC(IFNB vs Mock)
#   suppressed  retention <= R_STRONG (0.25) and the IFNB+EBOV-vs-IFNB contrast moves the gene back toward Mock with padj < PADJ_CUT
#   not suppressed  retention >= R_FAIL (0.75)       partly suppressed  everything in between
#   "altered by EBOV alone": not IFNB-responsive (|log2FC_B| < LFC_NOT) but |log2FC(EBOV vs Mock)| > LFC_CUT with padj < PADJ_CUT
# Needs DESeq2 objects from 04 and pheatmap. Usage: Rscript 51_mock_ifnb_ebov_analysis.R   (PROJ from 'source setup_session.sh')
# Writes to <project>/requant_v2_outputs/figures_mock_ifnb_ebov/ (PNG + PDF + CSV).

suppressPackageStartupMessages({ library(DESeq2); library(readr); library(dplyr); library(tidyr); library(ggplot2) })
if (!requireNamespace("pheatmap", quietly = TRUE)) stop("Install pheatmap (conda install -c conda-forge r-pheatmap)")
select <- dplyr::select; first <- dplyr::first; rename <- dplyr::rename
slice <- dplyr::slice; filter <- dplyr::filter; desc <- dplyr::desc; count <- dplyr::count

LFC_CUT <- 1; PADJ_CUT <- 0.05; LFC_NOT <- 0.5; R_STRONG <- 0.25; R_FAIL <- 0.75
N_ISG_4 <- 30; N_ISG_2 <- 50; NTOP_PCA <- 500
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
OUT <- file.path(RNA, "figures_mock_ifnb_ebov"); dir.create(OUT, showWarnings = FALSE)

vsd <- readRDS(file.path(RNA, "vsd_v2.rds")); vm <- assay(vsd)
sel <- tibble(sample = unlist(SETS), group = rep(names(SETS), lengths(SETS)), replicate = unlist(lapply(SETS, seq_along)))
sel$data_group <- as.character(vsd$group)[match(sel$sample, colnames(vm))]
if (any(is.na(sel$data_group))) stop("Samples not in the VST object: ", paste(sel$sample[is.na(sel$data_group)], collapse = ", "))
message("Check - chosen group vs group recorded in the data:")
print(as.data.frame(sel |> transmute(sample, chosen_as = group, replicate, group_in_data = data_group,
                                     match = ifelse(data_group == EXPECTED[group], "ok", "DIFFERENT"))), row.names = FALSE)
if (any(sel$data_group != EXPECTED[sel$group])) message("WARNING: some samples are recorded under a different group. Using YOUR sample numbers.")

pc <- read_csv(file.path(RNA, "DE_v2_IFNB_vs_Mock_protein_coding.csv"), show_col_types = FALSE)
pc <- pc |> filter(gene_id %in% rownames(vm), !is.na(padj))
X <- vm[, sel$sample, drop = FALSE]
mock_cols <- sel$sample[sel$group == "Mock"]
save_gt <- function(h, name, w, ht) for (ext in c("png", "pdf")) {
  f <- file.path(OUT, paste0(name, ".", ext))
  if (ext == "png") png(f, width = w, height = ht, units = "in", res = 300) else pdf(f, width = w, height = ht)
  grid::grid.newpage(); grid::grid.draw(h$gtable); dev.off()
}
pal <- colorRampPalette(c("#2C5F9E", "#F7F7F7", "#B2382C"))(101)
cap_of <- function(m) min(8, max(2, ceiling(quantile(abs(m), 0.98, na.rm = TRUE))))
change_matrix <- function(ids, cols) { m <- X[ids, cols, drop = FALSE]; m - rowMeans(X[ids, mock_cols, drop = FALSE]) }
col_ann <- function(cols) {
  grp <- sel$group[match(cols, sel$sample)]
  list(df = data.frame(Group = factor(grp, levels = names(SETS)), row.names = cols),
       gaps = unname(cumsum(table(factor(grp, levels = unique(grp))))[-length(unique(grp))]))
}

# ======================= 1. PCA ==========================================================================================
x <- vm[rownames(vm) %in% pc$gene_id, sel$sample, drop = FALSE]
v <- apply(x, 1, var); x <- x[order(v, decreasing = TRUE)[seq_len(min(NTOP_PCA, sum(v > 0)))], , drop = FALSE]
p <- prcomp(t(x), center = TRUE, scale. = FALSE); pct <- 100 * p$sdev^2 / sum(p$sdev^2)
message(sprintf("\nPCA: 12 samples, %d genes; variance explained PC1 %.1f%%, PC2 %.1f%%, PC3 %.1f%%", nrow(x), pct[1], pct[2], pct[3]))
co <- as.data.frame(p$x[, 1:5]) |> tibble::rownames_to_column("sample") |> left_join(sel |> select(sample, group, replicate), by = "sample") |>
  mutate(group = factor(group, levels = names(SETS)), txt = ifelse(group == "IFNB+EBOV", INK, "white"))
write_csv(co, file.path(OUT, "PCA_coordinates.csv"))
write_csv(as.data.frame(p$rotation[, 1:3]) |> tibble::rownames_to_column("gene_id") |> left_join(pc |> select(gene_id, gene_name), by = "gene_id"),
          file.path(OUT, "PCA_loadings.csv"))
for (pp in list(c("PC1", "PC2"), c("PC1", "PC3"))) {
  hull <- co |> group_by(group) |> slice(chull(.data[[pp[1]]], .data[[pp[2]]])) |> ungroup()
  ix <- as.integer(sub("PC", "", pp[1])); iy <- as.integer(sub("PC", "", pp[2]))
  g <- ggplot(co, aes(.data[[pp[1]]], .data[[pp[2]]])) +
    geom_hline(yintercept = 0, colour = "grey85", linewidth = 0.3) + geom_vline(xintercept = 0, colour = "grey85", linewidth = 0.3) +
    geom_polygon(data = hull, aes(fill = group, group = group), alpha = 0.13, colour = NA) +
    geom_point(aes(fill = group), shape = 21, size = 8, colour = "white", stroke = 0.9) +
    geom_text(aes(label = replicate, colour = txt), size = 3.6, fontface = "bold") +
    scale_colour_identity() + scale_fill_manual(values = GCOL, name = "group") +
    guides(fill = guide_legend(override.aes = list(shape = 21, size = 5, alpha = 1))) +
    scale_x_continuous(expand = expansion(mult = 0.12)) + scale_y_continuous(expand = expansion(mult = 0.12)) +
    labs(title = "PCA: Mock, IFNB, EBOV and IFNB+EBOV",
         subtitle = sprintf("top %d most variable protein-coding genes (VST), 12 samples; number inside each dot = replicate", nrow(x)),
         x = sprintf("%s (%.1f%% of variance)", pp[1], pct[ix]), y = sprintf("%s (%.1f%% of variance)", pp[2], pct[iy])) +
    theme_bw(base_size = 12) + theme(panel.grid.minor = element_blank(), panel.border = element_rect(colour = "grey70"),
          axis.text = element_text(colour = INK), plot.title = element_text(colour = INK, face = "bold"),
          plot.subtitle = element_text(colour = "grey35", size = 9), legend.position = "right")
  stub <- paste0("PCA_Mock_IFNB_EBOV_IFNBEBOV_", pp[1], "_", pp[2])
  ggsave(file.path(OUT, paste0(stub, ".png")), g, width = 8.5, height = 6.5, dpi = 300)
  ggsave(file.path(OUT, paste0(stub, ".pdf")), g, width = 8.5, height = 6.5)
}

# ======================= contrasts for the organised transcriptome heatmap ================================================
dds <- readRDS(file.path(RNA, "dds_v2_decoyAware_hostNormalised.rds"))
dds <- dds[rownames(dds) %in% pc$gene_id, ]
use_ashr <- requireNamespace("ashr", quietly = TRUE)
contrast_tab <- function(nm, a, b) {
  ct <- c("group", a, b); res <- results(dds, contrast = ct)
  shr <- if (use_ashr) lfcShrink(dds, contrast = ct, res = res, type = "ashr") else lfcShrink(dds, contrast = ct, type = "normal")
  out <- tibble(gene_id = rownames(res), l = shr$log2FoldChange, p = res$padj); names(out)[2:3] <- paste0(c("lfc_", "padj_"), nm); out
}
message("\nComputing contrasts (IFNB, EBOV, IFNB+EBOV vs Mock; IFNB+EBOV vs IFNB)...")
GT <- Reduce(function(a, b) left_join(a, b, by = "gene_id"),
            list(contrast_tab("B", "IFNB_Mock", "Mock_Mock"), contrast_tab("E", "Mock_EBOV", "Mock_Mock"),
                 contrast_tab("BE", "IFNB_EBOV", "Mock_Mock"), contrast_tab("BE_B", "IFNB_EBOV", "IFNB_Mock"))) |>
  left_join(pc |> select(gene_id, gene_name), by = "gene_id") |> relocate(gene_name, .after = gene_id)
GT <- GT |> mutate(
  resp = case_when(!is.na(padj_B) & padj_B < PADJ_CUT & lfc_B > LFC_CUT ~ "Up", !is.na(padj_B) & padj_B < PADJ_CUT & lfc_B < -LFC_CUT ~ "Down", TRUE ~ "None"),
  retain = ifelse(abs(lfc_B) > 0, lfc_BE / lfc_B, NA_real_),
  suppressed = resp != "None" & !is.na(retain) & retain <= R_STRONG & sign(lfc_BE_B) == -sign(lfc_B) & !is.na(padj_BE_B) & padj_BE_B < PADJ_CUT,
  not_supp = resp != "None" & !is.na(retain) & retain >= R_FAIL,
  ebov_alone = resp == "None" & abs(lfc_B) < LFC_NOT & !is.na(padj_E) & padj_E < PADJ_CUT & abs(lfc_E) > LFC_CUT,
  ebov_effect = case_when(resp != "None" & suppressed ~ "Suppressed", resp != "None" & not_supp ~ "Not suppressed", resp != "None" ~ "Partly suppressed",
                          ebov_alone ~ "EBOV-altered (no IFNB response)", TRUE ~ "Unchanged"))
BLOCKS <- tibble::tribble(~resp, ~eff, ~title,
  "Up", "Suppressed", "IFNB-up, suppressed by EBOV", "Up", "Partly suppressed", "IFNB-up, partly suppressed", "Up", "Not suppressed", "IFNB-up, NOT suppressed",
  "Down", "Suppressed", "IFNB-down, reversed by EBOV", "Down", "Partly suppressed", "IFNB-down, partly reversed", "Down", "Not suppressed", "IFNB-down, NOT reversed",
  "None", "EBOV-altered (no IFNB response)", "EBOV-altered only", "None", "Unchanged", "Unchanged")
GT$block <- NA_character_
for (i in seq_len(nrow(BLOCKS))) GT$block[GT$resp == BLOCKS$resp[i] & GT$ebov_effect == BLOCKS$eff[i]] <- BLOCKS$title[i]
write_csv(GT, file.path(OUT, "transcriptome_blocks.csv"))
message("\nGenes per block:"); print(as.data.frame(GT |> count(block) |> arrange(match(block, BLOCKS$title))), row.names = FALSE)

# ======================= 2. organised whole-transcriptome heatmap =============================================================
ord_cols <- unlist(SETS[c("Mock", "IFNB", "IFNB+EBOV", "EBOV")])
ordered_genes <- bind_rows(lapply(seq_len(nrow(BLOCKS)), function(i) {
  b <- GT |> filter(block == BLOCKS$title[i])
  b[order(if (BLOCKS$resp[i] == "Down") b$lfc_B else if (BLOCKS$resp[i] == "None" && BLOCKS$eff[i] != "Unchanged") -b$lfc_E else -b$lfc_B), ]
})) |> filter(gene_id %in% rownames(X))
m <- change_matrix(ordered_genes$gene_id, ord_cols)
ca <- col_ann(ord_cols)
ann_row <- data.frame(IFNB_response = factor(ordered_genes$resp, levels = c("Up", "Down", "None")),
                      EBOV_effect = factor(ordered_genes$ebov_effect, levels = c("Suppressed", "Partly suppressed", "Not suppressed", "EBOV-altered (no IFNB response)", "Unchanged")),
                      row.names = ordered_genes$gene_id)
blk <- factor(ordered_genes$block, levels = BLOCKS$title)
ann_cols <- list(Group = GCOL[names(SETS)],
                 IFNB_response = c(Up = "#C53030", Down = "#2B6CB0", None = "#D9D9D9"),
                 EBOV_effect = c(Suppressed = "#762A83", `Partly suppressed` = "#C2A5CF", `Not suppressed` = "#F1A340",
                                 `EBOV-altered (no IFNB response)` = "#1B7837", Unchanged = "#EFEFEF"))
h2 <- pheatmap::pheatmap(m, color = pal, breaks = seq(-cap_of(m), cap_of(m), length.out = 102), cluster_rows = FALSE, cluster_cols = FALSE,
                         annotation_col = ca$df, annotation_row = ann_row, annotation_colors = ann_cols,
                         gaps_col = ca$gaps, gaps_row = unname(cumsum(table(blk)))[-nlevels(blk)][table(blk)[-nlevels(blk)] > 0],
                         show_rownames = FALSE, show_colnames = TRUE, border_color = NA, fontsize = 9, silent = TRUE,
                         main = sprintf("All protein-coding genes (%d), organised by IFNB response and EBOV effect\n(suppressed / reversed = EBOV pulls the IFNB change back toward Mock)", nrow(m)))
save_gt(h2, "heatmap_all_protein_coding_organised", 9, 12)
message("Block order top to bottom: ", paste(BLOCKS$title, collapse = " | "))

# ======================= 3 and 4. ISG heatmaps ===================================================================================
isg <- pc |> filter(padj < PADJ_CUT, log2FoldChange > LFC_CUT) |> arrange(desc(log2FoldChange))
message(sprintf("\nISGs (padj<%g, log2FC>%d): %d", PADJ_CUT, LFC_CUT, nrow(isg)))
isg_heatmap <- function(n, cond_names, name, title) {
  g <- head(isg, n); g <- g[g$gene_id %in% rownames(X), ]
  cols <- unlist(SETS[cond_names]); mm <- change_matrix(g$gene_id, cols); cc <- col_ann(cols)
  h <- pheatmap::pheatmap(mm, color = pal, breaks = seq(-cap_of(mm), cap_of(mm), length.out = 102), cluster_rows = FALSE, cluster_cols = FALSE,
                          annotation_col = cc$df, annotation_colors = list(Group = GCOL[cond_names]), gaps_col = cc$gaps,
                          labels_row = make.unique(sprintf("%s  (%+.1f)", g$gene_name, g$log2FoldChange)), fontsize_row = 8.5, fontsize = 9,
                          show_colnames = TRUE, border_color = NA, silent = TRUE, main = sprintf("%s (label = gene (IFNB log2FC))", title))
  save_gt(h, name, if (length(cond_names) > 2) 7 else 5.5, 0.21 * nrow(mm) + 3.2)
  write_csv(g |> transmute(rank = row_number(), gene_id, gene_name, log2FoldChange, padj), file.path(OUT, paste0(name, "_genes.csv")))
}
isg_heatmap(N_ISG_4, c("Mock", "IFNB", "IFNB+EBOV", "EBOV"), "heatmap_top30_ISGs_4conditions", sprintf("Top %d IFNB ISGs: Mock, IFNB, IFNB+EBOV, EBOV", N_ISG_4))
isg_heatmap(N_ISG_2, c("Mock", "IFNB"), "heatmap_top50_ISGs_Mock_IFNB", sprintf("Top %d IFNB ISGs: Mock and IFNB", N_ISG_2))
isg_heatmap(N_ISG_4, c("Mock", "IFNB"), "heatmap_top30_ISGs_Mock_IFNB", sprintf("Top %d IFNB ISGs: Mock and IFNB", N_ISG_4))
message("\nFigures and tables in: ", OUT)
