#!/usr/bin/env Rscript
# Genes that EBOV UPREGULATES (EBOV vs Mock) but IFNB DOWNREGULATES (IFNB vs Mock), by any fold-change.
#   EBOV up   = padj < PADJ_CUT and shrunk log2FC > LFC_E (default 0, i.e. any increase)
#   IFNB down = padj < PADJ_CUT and shrunk log2FC < -LFC_B (default 0, i.e. any decrease)         => class "both"
#   a second, looser class "IFNB trend" = EBOV up (significant) and IFNB log2FC < 0 but not significant, shown lighter so you can judge the threshold yourself
# Figures (PNG + PDF) in requant_v2_outputs/figures_ebov_up_ifnb_down/:
#   scatter_EBOV_vs_IFNB_log2FC      every gene: x = log2FC IFNB vs Mock, y = log2FC EBOV vs Mock; the upper-left quadrant is the question; hits coloured and the top ones labelled
#   venn_EBOV_up_IFNB_down           how many genes are EBOV-up, IFNB-down, and in both (hypergeometric p of the overlap is in the subtitle)
#   bars_EBOV_up_IFNB_down           top N_BAR hits: log2FC in IFNB vs Mock (blue) and EBOV vs Mock (green), with padj stars
#   heatmap_EBOV_up_IFNB_down_*      top N_HEAT hits across Untreated | IFNB | EBOV | IFNB+EBOV replicates: _mockgrey (each replicate vs the Untreated mean, Mock = grey)
#                                    and _rowZ (row z-score); rows ordered by hit score, not clustered
#   genes_EBOV_up_IFNB_down.csv      every hit with log2FC / padj for IFNB, EBOV and IFNB+EBOV vs Mock (so you can see whether EBOV also wins in the co-treatment)
# hit score = log2FC(EBOV) - log2FC(IFNB): large when EBOV pushes the gene up and IFNB pushes it down. Uses all genes with an Ensembl ID; set PROTEIN_CODING=1 for protein-coding only.
# Usage: Rscript 79_ebov_up_ifnb_down.R   (PROJ from 'source setup_session.sh')

suppressPackageStartupMessages({ library(readr); library(dplyr); library(tidyr); library(ggplot2); library(ggrepel); library(pheatmap); library(grid) })
select <- dplyr::select; filter <- dplyr::filter

PADJ_CUT <- 0.05; LFC_E <- 0; LFC_B <- 0; N_LAB <- 20; N_BAR <- 30; N_HEAT <- 50; CAP_Q <- 0.98
PROTEIN_CODING <- identical(Sys.getenv("PROTEIN_CODING"), "1")
INK <- "#1F2933"; C_B <- "#2B6CB0"; C_E <- "#1B7837"; C_HIT <- "#7B3294"; C_TREND <- "#C2A5CF"
G <- c(Untreated = "Mock_Mock", IFNB = "IFNB_Mock", EBOV = "Mock_EBOV", `IFNB+EBOV` = "IFNB_EBOV")
PROJ <- Sys.getenv("PROJ", unset = file.path(
  Sys.getenv("HOME"),
  "Documents/aviti_rnaseq/2025/aviti",
  "20251009-100725_Maury_25161_WivfFacICOTHFAUJwSPeHBEFfCfNtCRpUSgVGUkz",
  "results/Maury_25161"))
RNA <- file.path(PROJ, "requant_v2_outputs")
OUT <- file.path(RNA, "figures_ebov_up_ifnb_down"); dir.create(OUT, showWarnings = FALSE)
save2 <- function(p, name, w, h) { ggsave(file.path(OUT, paste0(name, ".png")), p, width = w, height = h, dpi = 300); ggsave(file.path(OUT, paste0(name, ".pdf")), p, width = w, height = h) }
stars <- function(p) ifelse(is.na(p), "", ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", "ns"))))
rd <- function(f) read_csv(file.path(RNA, f), show_col_types = FALSE)

E <- rd("DE_v2_EBOV_vs_Mock_annotated.csv") |> filter(grepl("^ENSMUSG", gene_id)) |> transmute(gene_id, gene_name, gene_biotype, baseMean, lfc_E = log2FoldChange, padj_E = padj)
B <- rd("DE_v2_IFNB_vs_Mock_annotated.csv") |> transmute(gene_id, lfc_B = log2FoldChange, padj_B = padj)
BE <- rd("DE_v2_IFNB_EBOV_vs_Mock_annotated.csv") |> transmute(gene_id, lfc_BE = log2FoldChange, padj_BE = padj)
A <- E |> inner_join(B, by = "gene_id") |> left_join(BE, by = "gene_id") |> filter(!is.na(padj_E), !is.na(padj_B))
if (PROTEIN_CODING) A <- A |> filter(gene_biotype == "protein_coding")
A <- A |> mutate(label = ifelse(is.na(gene_name) | gene_name == "", gene_id, gene_name), score = lfc_E - lfc_B,
                 e_up = padj_E < PADJ_CUT & lfc_E > LFC_E, b_down = padj_B < PADJ_CUT & lfc_B < -LFC_B,
                 class = case_when(e_up & b_down ~ "both", e_up & lfc_B < 0 ~ "IFNB trend", TRUE ~ "other"))
nE <- sum(A$e_up); nB <- sum(A$b_down); nO <- sum(A$e_up & A$b_down)
ph <- phyper(nO - 1, nB, nrow(A) - nB, nE, lower.tail = FALSE)
message(sprintf("%d genes tested in both comparisons | EBOV up: %d | IFNB down: %d | both: %d (expected by chance %.1f; hypergeometric p = %.2g) | EBOV up with a non-significant IFNB decrease: %d",
                nrow(A), nE, nB, nO, nE * nB / nrow(A), ph, sum(A$class == "IFNB trend")))
hits <- A |> filter(class == "both") |> arrange(desc(score))
write_csv(A |> filter(class != "other") |> arrange(class, desc(score)) |> select(gene_id, gene_name, gene_biotype, baseMean, class, score, lfc_E, padj_E, lfc_B, padj_B, lfc_BE, padj_BE), file.path(OUT, "genes_EBOV_up_IFNB_down.csv"))
if (!nrow(hits)) stop("No gene is significantly up in EBOV and down in IFNB at these cut-offs - see genes_EBOV_up_IFNB_down.csv for the trend class")
print(as.data.frame(hits |> slice_head(n = 15) |> transmute(gene = label, IFNB = round(lfc_B, 2), EBOV = round(lfc_E, 2), IFNB_EBOV = round(lfc_BE, 2), padj_IFNB = signif(padj_B, 2), padj_EBOV = signif(padj_E, 2))), row.names = FALSE)

base_theme <- theme_bw(base_size = 12) + theme(panel.grid.minor = element_blank(), panel.grid.major = element_line(colour = "grey94", linewidth = 0.3), panel.border = element_rect(colour = "grey70"),
  axis.text = element_text(colour = INK), axis.title = element_text(colour = INK), plot.title = element_text(colour = INK, face = "bold"), plot.subtitle = element_text(colour = "grey35", size = 9))

# ---- 1. scatter ----
lab <- hits |> slice_head(n = N_LAB)
lim <- quantile(c(abs(A$lfc_B), abs(A$lfc_E)), 0.999)
ps <- ggplot(A, aes(lfc_B, lfc_E)) +
  geom_point(data = A |> filter(class == "other"), colour = "grey70", size = 0.5, alpha = 0.4) +
  geom_point(data = A |> filter(class == "IFNB trend"), colour = C_TREND, size = 1, alpha = 0.8) +
  geom_point(data = hits, colour = C_HIT, size = 1.6, alpha = 0.95) +
  geom_text_repel(data = lab, aes(label = label), colour = C_HIT, size = 3, fontface = "italic", max.overlaps = Inf, min.segment.length = 0, segment.size = 0.25, segment.colour = "grey60", box.padding = 0.4, seed = 1) +
  geom_hline(yintercept = 0, colour = INK, linewidth = 0.35) + geom_vline(xintercept = 0, colour = INK, linewidth = 0.35) + geom_abline(slope = -1, intercept = 0, linetype = "dotted", colour = "grey55") +
  coord_cartesian(xlim = c(-lim, lim), ylim = c(-lim, lim)) +
  annotate("text", x = -lim * 0.98, y = lim * 0.98, label = sprintf("EBOV up, IFNB down\n%d genes (purple)", nrow(hits)), colour = C_HIT, hjust = 0, vjust = 1, fontface = "bold", size = 4) +
  labs(title = "Genes EBOV raises that IFNB lowers", subtitle = sprintf("purple = EBOV up and IFNB down (both padj < %g); light purple = EBOV up, IFNB log2FC < 0 but not significant; grey = other genes; dotted line = equal and opposite", PADJ_CUT),
       x = "log2FC IFNB vs Mock", y = "log2FC EBOV vs Mock") + base_theme
save2(ps, "scatter_EBOV_vs_IFNB_log2FC", 7.5, 7)

# ---- 2. Venn ----
circ <- function(cx, r = 1, n = 200) tibble(x = cx + r * cos(seq(0, 2 * pi, length.out = n)), y = r * sin(seq(0, 2 * pi, length.out = n)))
pv <- ggplot() + geom_polygon(data = circ(-0.6), aes(x, y), fill = C_E, alpha = 0.35, colour = C_E) + geom_polygon(data = circ(0.6), aes(x, y), fill = C_B, alpha = 0.35, colour = C_B) +
  annotate("text", x = -1.25, y = 0, label = nE - nO, size = 7, fontface = "bold", colour = INK) + annotate("text", x = 0, y = 0, label = nO, size = 8, fontface = "bold", colour = C_HIT) +
  annotate("text", x = 1.25, y = 0, label = nB - nO, size = 7, fontface = "bold", colour = INK) +
  annotate("text", x = -0.9, y = 1.25, label = sprintf("EBOV up (%d)", nE), colour = C_E, fontface = "bold", size = 4.2) + annotate("text", x = 0.9, y = 1.25, label = sprintf("IFNB down (%d)", nB), colour = C_B, fontface = "bold", size = 4.2) +
  coord_equal(xlim = c(-2, 2), ylim = c(-1.2, 1.5)) + theme_void(base_size = 12) +
  labs(title = "EBOV-up genes that are IFNB-down", subtitle = sprintf("padj < %g, any fold-change; overlap expected by chance %.1f, hypergeometric p = %.2g; %d genes tested in both", PADJ_CUT, nE * nB / nrow(A), ph, nrow(A))) +
  theme(plot.title = element_text(colour = INK, face = "bold", hjust = 0.5), plot.subtitle = element_text(colour = "grey35", size = 9, hjust = 0.5), plot.margin = margin(10, 10, 10, 10))
save2(pv, "venn_EBOV_up_IFNB_down", 6, 4.2)

# ---- 3. paired bars ----
tb <- hits |> slice_head(n = N_BAR)
L <- bind_rows(tb |> transmute(label, cond = "IFNB vs Mock", lfc = lfc_B, padj = padj_B), tb |> transmute(label, cond = "EBOV vs Mock", lfc = lfc_E, padj = padj_E)) |>
  mutate(cond = factor(cond, levels = c("IFNB vs Mock", "EBOV vs Mock")), label = factor(label, levels = rev(tb$label)), hj = ifelse(lfc >= 0, 0, 1), xs = lfc + ifelse(lfc >= 0, 0.04, -0.04) * diff(range(c(tb$lfc_B, tb$lfc_E))))
pb <- ggplot(L, aes(lfc, label, fill = cond)) + geom_col(position = position_dodge(width = 0.8), width = 0.75) +
  geom_text(aes(x = xs, label = stars(padj), hjust = hj), position = position_dodge(width = 0.8), size = 2.6, vjust = 0.45, colour = INK) +
  geom_vline(xintercept = 0, colour = INK, linewidth = 0.4) + scale_fill_manual(values = c(`IFNB vs Mock` = C_B, `EBOV vs Mock` = C_E), name = NULL) + scale_x_continuous(expand = expansion(mult = 0.1)) +
  labs(title = sprintf("Top %d genes EBOV raises and IFNB lowers", nrow(tb)), subtitle = "shrunk log2FC vs Mock, ranked by log2FC(EBOV) - log2FC(IFNB); stars = padj", x = "log2 fold-change vs Mock", y = NULL) + base_theme +
  theme(axis.text.y = element_text(face = "italic"), legend.position = "top", panel.grid.major.y = element_blank())
save2(pb, "bars_EBOV_up_IFNB_down", 7.5, 0.28 * nrow(tb) + 2)

# ---- 4. heatmaps ----
vsd <- readRDS(file.path(RNA, "vsd_v2.rds")); vm <- SummarizedExperiment::assay(vsd); grp <- as.character(vsd$group)
cols <- lapply(G, function(g) colnames(vm)[grp == g]); cols <- cols[lengths(cols) > 0]
th <- hits |> slice_head(n = N_HEAT) |> mutate(label = make.unique(label)); th <- th[th$gene_id %in% rownames(vm), ]
nrep <- lengths(cols)
X <- vm[th$gene_id, unlist(cols), drop = FALSE]; rownames(X) <- th$label
colnames(X) <- paste0(rep(names(cols), nrep), " ", unlist(lapply(nrep, seq_len)))
W <- 0.5 * ncol(X) + 2.5; H <- 0.17 * nrow(X) + 2.2
save_ph <- function(mat, name, pal, brk, at, text, title) {
  p <- pheatmap(mat, color = pal, breaks = brk, cluster_rows = FALSE, cluster_cols = FALSE, border_color = NA, fontsize_row = 8, fontsize_col = 9, fontsize = 9, angle_col = 45, cellwidth = 22, cellheight = 10,
                gaps_col = cumsum(nrep)[-length(nrep)], silent = TRUE, legend_breaks = at, legend_labels = text, main = title)
  png(file.path(OUT, paste0(name, ".png")), width = W, height = H, units = "in", res = 300); grid.newpage(); grid.draw(p$gtable); invisible(dev.off())
  pdf(file.path(OUT, paste0(name, ".pdf")), width = W, height = H); grid.newpage(); grid.draw(p$gtable); invisible(dev.off())
}
Xm <- X - rowMeans(X[, seq_len(nrep[[1]]), drop = FALSE]); lm <- as.numeric(quantile(abs(Xm), CAP_Q))
save_ph(Xm, "heatmap_EBOV_up_IFNB_down_mockgrey", colorRampPalette(c(C_B, "#F2F2F2", "#C53030"))(100), seq(-lm, lm, length.out = 101), c(-lm, 0, lm), c(sprintf("%.1f", -lm), "0 (Mock mean)", sprintf("%.1f", lm)),
        sprintf("Top %d EBOV-up / IFNB-down genes: each replicate vs Untreated mean", nrow(Xm)))
Z <- t(scale(t(X))); lz <- as.numeric(quantile(abs(Z), CAP_Q))
save_ph(Z, "heatmap_EBOV_up_IFNB_down_rowZ", colorRampPalette(c(C_B, "white", "#C53030"))(100), seq(-lz, lz, length.out = 101), c(-lz, 0, lz), c(sprintf("%.1f", -lz), "0", sprintf("%.1f", lz)),
        sprintf("Top %d EBOV-up / IFNB-down genes: expression, row z-score", nrow(Z)))
message("Saved to: ", OUT)
