#!/usr/bin/env Rscript
# How differently do EBOV and MARV suppress the IFNB response, gene by gene?
# Needs 48_patterns_virus_vs_ifnb.R to have been run (reads figures_patterns/master_table.csv).
#
# For every IFNB-responsive protein-coding gene (|log2FC_B| > LFC_CUT, padj_B < PADJ_CUT):
#   suppression by virus V  S_V = -sign(log2FC_B) x log2FC(IFNB+V vs IFNB)      (log2 units; > 0 = the virus pushes the gene
#                                                                                back toward Mock, < 0 = it enhances the response)
#   difference               dS  = S_E - S_M   (> 0 EBOV suppresses more, < 0 MARV suppresses more)
# and a direct statistical test of IFNB+EBOV vs IFNB+MARV (DESeq2, shrunk log2FC, padj) tells whether the two viruses
# really differ. Categories (genes up and down with IFNB are kept apart):
#   EBOV suppresses more     padj < PADJ_CUT, |log2FC(IFNB+EBOV vs IFNB+MARV)| >= DIFF_LFC, dS > 0
#   MARV suppresses more     same, dS < 0
#   No difference - both suppress      not significant and |dS| < NO_DIFF, retention <= 0.5 for both viruses
#   No difference - neither suppresses  not significant and |dS| < NO_DIFF, retention >= 0.75 for both
#   No difference - intermediate        not significant and |dS| < NO_DIFF, anything else
#   Inconclusive                        everything else
# retention = log2FC(IFNB+virus vs Mock) / log2FC(IFNB vs Mock) (1 = full response kept, 0 = gone).
#
# Writes to <project>/requant_v2_outputs/figures_patterns/virus_comparison/ (PNG + PDF + CSV):
#   EV_A_suppression_scatter   S_MARV against S_EBOV, one dot per gene, diagonal = equal suppression
#   EV_B_difference_vs_mean    dS against mean suppression (is the difference large where suppression is large?)
#   EV_C_category_counts       how many genes fall in each category
#   EV_D_heatmap_IFNB_up / EV_D_heatmap_IFNB_down   representative genes: strongest EBOV-biased, strongest MARV-biased, and
#                              no-difference genes that both viruses suppress strongly (columns as in the earlier heatmaps)
#   genes_<category>.csv, virus_comparison_table.csv
# Usage: Rscript 50_ebov_vs_marv_suppression.R   (PROJ from 'source setup_session.sh')

suppressPackageStartupMessages({ library(DESeq2); library(readr); library(dplyr); library(tidyr); library(ggplot2) })
for (pkg in c("ggrepel", "pheatmap")) if (!requireNamespace(pkg, quietly = TRUE)) stop("Install ", pkg)
select <- dplyr::select; first <- dplyr::first; rename <- dplyr::rename
slice <- dplyr::slice; filter <- dplyr::filter; desc <- dplyr::desc; count <- dplyr::count

LFC_CUT <- 1; PADJ_CUT <- 0.05; DIFF_LFC <- 0.58; NO_DIFF <- 0.3; N_EACH <- 20; N_LABEL <- 8
INK <- "#1F2933"
SETS <- list(Mock = c("1", "2", "3"), IFNB = c("10", "11", "12"), `IFNB+EBOV` = c("16", "17", "18"),
             EBOV = c("4", "5", "6"), `IFNB+MARV` = c("22", "23", "24"), MARV = c("7", "8", "9"))
GCOL <- c(Mock = "#7F7F7F", IFNB = "#2B6CB0", `IFNB+EBOV` = "#9BD08E", EBOV = "#1B7837", `IFNB+MARV` = "#C9A6DB", MARV = "#762A83")
CATS <- c("EBOV suppresses more", "MARV suppresses more", "No difference: both suppress",
          "No difference: neither suppresses", "No difference: intermediate", "Inconclusive")
CCOL <- c(`EBOV suppresses more` = "#1B7837", `MARV suppresses more` = "#762A83", `No difference: both suppress` = "#2B6CB0",
          `No difference: neither suppresses` = "#9E9E9E", `No difference: intermediate` = "#C8C8C8", Inconclusive = "#E6E6E6")

PROJ <- Sys.getenv("PROJ", unset = file.path(
  Sys.getenv("HOME"),
  "Documents/aviti_rnaseq/2025/aviti",
  "20251009-100725_Maury_25161_WivfFacICOTHFAUJwSPeHBEFfCfNtCRpUSgVGUkz",
  "results/Maury_25161"))
RNA <- file.path(PROJ, "requant_v2_outputs"); PAT <- file.path(RNA, "figures_patterns")
OUT <- file.path(PAT, "virus_comparison"); dir.create(OUT, showWarnings = FALSE)
save2 <- function(p, name, w, h) {
  ggsave(file.path(OUT, paste0(name, ".png")), p, width = w, height = h, dpi = 300, limitsize = FALSE)
  ggsave(file.path(OUT, paste0(name, ".pdf")), p, width = w, height = h, limitsize = FALSE)
}
theme_f <- function() theme_bw(base_size = 11) + theme(panel.grid.minor = element_blank(), panel.border = element_rect(colour = "grey70"),
  axis.text = element_text(colour = INK), strip.background = element_rect(fill = "grey95", colour = "grey70"),
  strip.text = element_text(colour = INK, face = "bold"), plot.title = element_text(colour = INK, face = "bold"),
  plot.subtitle = element_text(colour = "grey35", size = 8.5), legend.position = "bottom")

mt <- file.path(PAT, "master_table.csv")
if (!file.exists(mt)) stop("Run 48_patterns_virus_vs_ifnb.R first (missing ", mt, ")")
M <- read_csv(mt, show_col_types = FALSE)

# ---- direct IFNB+EBOV vs IFNB+MARV test ----------------------------------------------------------------------------------------
dds <- readRDS(file.path(RNA, "dds_v2_decoyAware_hostNormalised.rds"))
dds <- dds[rownames(dds) %in% M$gene_id, ]
ct <- c("group", "IFNB_EBOV", "IFNB_MARV")
res <- results(dds, contrast = ct)
shr <- if (requireNamespace("ashr", quietly = TRUE)) lfcShrink(dds, contrast = ct, res = res, type = "ashr") else lfcShrink(dds, contrast = ct, type = "normal")
M <- M |> left_join(tibble(gene_id = rownames(res), lfc_BEvsBM = shr$log2FoldChange, padj_BEvsBM = res$padj), by = "gene_id")

R <- M |> filter(!is.na(padj_B), padj_B < PADJ_CUT, abs(lfc_B) > LFC_CUT, !is.na(lfc_BE_B), !is.na(lfc_BM_B)) |>
  mutate(direction = factor(ifelse(lfc_B > 0, "IFNB-up genes", "IFNB-down genes"), levels = c("IFNB-up genes", "IFNB-down genes")),
         S_E = -sign(lfc_B) * lfc_BE_B, S_M = -sign(lfc_B) * lfc_BM_B, dS = S_E - S_M, S_mean = (S_E + S_M) / 2,
         diff_sig = !is.na(padj_BEvsBM) & padj_BEvsBM < PADJ_CUT & abs(lfc_BEvsBM) >= DIFF_LFC,
         category = case_when(diff_sig & dS > 0 ~ CATS[1], diff_sig & dS < 0 ~ CATS[2],
                              !diff_sig & abs(dS) < NO_DIFF & retain_E <= 0.5 & retain_M <= 0.5 ~ CATS[3],
                              !diff_sig & abs(dS) < NO_DIFF & retain_E >= 0.75 & retain_M >= 0.75 ~ CATS[4],
                              !diff_sig & abs(dS) < NO_DIFF ~ CATS[5], TRUE ~ CATS[6]),
         category = factor(category, levels = CATS))
write_csv(R, file.path(OUT, "virus_comparison_table.csv"))
for (cn in CATS) for (dn in levels(R$direction)) {
  g <- R |> filter(category == cn, direction == dn); if (nrow(g))
    write_csv(g, file.path(OUT, paste0("genes_", gsub("[^A-Za-z]+", "_", cn), "__", sub("-", "_", sub(" ", "_", dn)), ".csv")))
}
message(sprintf("IFNB responders: %d (%d up, %d down)", nrow(R), sum(R$direction == "IFNB-up genes"), sum(R$direction == "IFNB-down genes")))
print(as.data.frame(R |> count(direction, category) |> pivot_wider(names_from = direction, values_from = n, values_fill = 0)), row.names = FALSE)
message(sprintf("Median suppression: EBOV %.2f, MARV %.2f log2 units; correlation of S_E and S_M = %.2f", median(R$S_E), median(R$S_M), cor(R$S_E, R$S_M)))

# ---- A: S_MARV against S_EBOV --------------------------------------------------------------------------------------------------
lim <- quantile(c(R$S_E, R$S_M), c(0.005, 0.995)); lo <- floor(lim[[1]]); hi <- ceiling(lim[[2]])
lab <- R |> filter(category %in% CATS[1:2]) |> group_by(direction, category) |> slice_max(abs(dS), n = N_LABEL, with_ties = FALSE) |> ungroup()
pA <- ggplot(R, aes(S_E, S_M)) +
  geom_hline(yintercept = 0, colour = "grey80", linewidth = 0.3) + geom_vline(xintercept = 0, colour = "grey80", linewidth = 0.3) +
  geom_abline(slope = 1, intercept = 0, colour = "grey40", linewidth = 0.4) +
  geom_point(data = R |> filter(category %in% CATS[3:6]), aes(colour = category), size = 1.2, alpha = 0.6) +
  geom_point(data = R |> filter(category %in% CATS[1:2]), aes(colour = category), size = 1.9, alpha = 0.9) +
  ggrepel::geom_text_repel(data = lab, aes(label = gene_name), size = 2.6, colour = INK, max.overlaps = Inf, min.segment.length = 0,
                           segment.colour = "grey60", segment.size = 0.25, box.padding = 0.3, seed = 1) +
  facet_wrap(~direction) + coord_fixed(xlim = c(lo, hi), ylim = c(lo, hi)) +
  scale_colour_manual(values = CCOL, name = NULL, drop = FALSE) +
  labs(title = "Suppression of the IFNB response: EBOV vs MARV, gene by gene",
       subtitle = "suppression = log2 units the virus pulls a gene back toward Mock (negative = enhances); diagonal = equal suppression; labelled genes differ most (significant)",
       x = "suppression by EBOV (IFNB+EBOV vs IFNB)", y = "suppression by MARV (IFNB+MARV vs IFNB)") + theme_f() + guides(colour = guide_legend(nrow = 2))
save2(pA, "EV_A_suppression_scatter", 12, 7.2)

# ---- B: difference against mean ---------------------------------------------------------------------------------------------------
pB <- ggplot(R, aes(S_mean, dS, colour = category)) +
  geom_hline(yintercept = 0, colour = "grey40", linewidth = 0.4) +
  geom_hline(yintercept = c(-DIFF_LFC, DIFF_LFC), linetype = "dashed", colour = "grey65", linewidth = 0.3) +
  geom_point(size = 1.4, alpha = 0.75) + facet_wrap(~direction) +
  scale_colour_manual(values = CCOL, name = NULL, drop = FALSE) +
  labs(title = "Where do EBOV and MARV differ?", subtitle = sprintf("y > 0: EBOV suppresses more, y < 0: MARV suppresses more; dashed = +/-%.2f log2; x = how strongly the two suppress on average", DIFF_LFC),
       x = "mean suppression by the two viruses (log2 units)", y = "EBOV minus MARV suppression (log2 units)") + theme_f() + guides(colour = guide_legend(nrow = 2))
save2(pB, "EV_B_difference_vs_mean", 12, 6)

# ---- C: counts ------------------------------------------------------------------------------------------------------------------------
cn <- R |> count(direction, category)
pC <- ggplot(cn, aes(direction, n, fill = category)) + geom_col(width = 0.65, colour = "white") +
  geom_text(aes(label = n), position = position_stack(vjust = 0.5), size = 3, colour = INK) +
  scale_fill_manual(values = CCOL, name = NULL, drop = FALSE) +
  labs(title = "How many IFNB-responsive genes differ between the two viruses?", x = NULL, y = "genes") + theme_f() + guides(fill = guide_legend(nrow = 3))
save2(pC, "EV_C_category_counts", 7, 6)

# ---- D: heatmaps of representative genes --------------------------------------------------------------------------------------
vsd <- readRDS(file.path(RNA, "vsd_v2.rds")); vm <- assay(vsd)
sel <- tibble(sample = unlist(SETS), group = rep(names(SETS), lengths(SETS)))
X <- vm[, sel$sample, drop = FALSE]; mock_cols <- sel$sample[sel$group == "Mock"]
ann_col <- data.frame(Group = factor(sel$group, levels = names(SETS)), row.names = sel$sample)
gaps_col <- unname(cumsum(lengths(SETS))[-length(SETS)])
pal <- colorRampPalette(c("#2C5F9E", "#F7F7F7", "#B2382C"))(101)
for (dn in levels(R$direction)) {
  d <- R |> filter(direction == dn)
  blocks <- list(`EBOV suppresses more` = d |> filter(category == CATS[1]) |> slice_max(dS, n = N_EACH, with_ties = FALSE),
                 `MARV suppresses more` = d |> filter(category == CATS[2]) |> slice_min(dS, n = N_EACH, with_ties = FALSE),
                 `No difference, both suppress` = d |> filter(category == CATS[3]) |> slice_max(S_mean, n = N_EACH, with_ties = FALSE))
  blocks <- blocks[vapply(blocks, nrow, integer(1)) > 0]
  if (!length(blocks)) next
  g <- bind_rows(lapply(names(blocks), function(b) blocks[[b]] |> mutate(block = b)))
  g <- g[g$gene_id %in% rownames(X), ]
  m <- X[g$gene_id, , drop = FALSE] - rowMeans(X[g$gene_id, mock_cols, drop = FALSE])
  cap <- min(8, max(2, ceiling(quantile(abs(m), 0.98, na.rm = TRUE))))
  rownames(m) <- make.unique(g$gene_id)
  ann_row <- data.frame(Comparison = factor(g$block, levels = names(blocks)), row.names = rownames(m))
  bc <- c(`EBOV suppresses more` = "#1B7837", `MARV suppresses more` = "#762A83", `No difference, both suppress` = "#2B6CB0")
  h <- pheatmap::pheatmap(m, color = pal, breaks = seq(-cap, cap, length.out = 102), cluster_rows = FALSE, cluster_cols = FALSE,
                          annotation_col = ann_col, annotation_row = ann_row, annotation_colors = list(Group = GCOL, Comparison = bc[names(blocks)]),
                          gaps_col = gaps_col,
                          gaps_row = if (length(blocks) > 1) unname(cumsum(table(factor(g$block, levels = names(blocks)))))[-length(blocks)] else NULL,
                          labels_row = make.unique(sprintf("%s  (%+.1f)", g$gene_name, g$lfc_B)), fontsize_row = 8, fontsize = 9,
                          border_color = NA, silent = TRUE,
                          main = sprintf("%s: EBOV- vs MARV-specific suppression and genes where they agree (label = gene (IFNB log2FC))", dn))
  stub <- paste0("EV_D_heatmap_", sub(" genes", "", sub("-", "_", dn)))
  ht <- 0.19 * nrow(m) + 3.5
  for (ext in c("png", "pdf")) {
    f <- file.path(OUT, paste0(stub, ".", ext))
    if (ext == "png") png(f, width = 9, height = ht, units = "in", res = 300) else pdf(f, width = 9, height = ht)
    grid::grid.newpage(); grid::grid.draw(h$gtable); dev.off()
  }
}
message("\nFigures and tables in: ", OUT)
