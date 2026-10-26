#!/usr/bin/env Rscript
# EBOV vs Mock (untreated): a volcano plot of every gene tested and a heatmap of the top TOP_N differentially expressed genes.
#   differentially expressed (DE) = padj < PADJ_CUT and |shrunk log2FC| > LFC_CUT (red = up, blue = down, black = not DE); the volcano shows ALL tested genes
#   (set PROTEIN_CODING=1 to restrict to protein-coding genes); the TOP_N_LAB most changed up and down genes are labelled.
#   heatmap: the TOP_N DE genes with the smallest padj (set TOP_BY=lfc to take the largest |log2FC| instead), rows in order of log2FC (up genes first), not clustered,
#   columns Untreated 1-3 | EBOV 1-3 (every replicate). Two versions:
#     heatmap_top50_EBOV_vs_Mock_rowZ        row z-score of VST expression (blue below the gene's mean across the 6 samples, red above)
#     heatmap_top50_EBOV_vs_Mock_mockgrey    each replicate's VST minus the gene's Untreated mean: Mock level = light grey / white, higher = red, lower = blue
# Reads DE_v2_EBOV_vs_Mock_annotated.csv and vsd_v2.rds. Writes to requant_v2_outputs/figures_ebov_vs_mock/ (PNG + PDF) plus EBOV_vs_Mock_DE_genes.csv and heatmap_top50_genes.csv
# Usage: Rscript 76_ebov_vs_mock_volcano_heatmap.R   (PROJ from 'source setup_session.sh';  PROTEIN_CODING=1 TOP_BY=lfc optional)

suppressPackageStartupMessages({ library(readr); library(dplyr); library(ggplot2); library(ggrepel); library(pheatmap); library(DESeq2); library(grid) })
select <- dplyr::select; filter <- dplyr::filter

TOP_N <- 50; TOP_N_LAB <- 10; PADJ_CUT <- 0.05; LFC_CUT <- 1; CAP_Q <- 0.98
PROTEIN_CODING <- identical(Sys.getenv("PROTEIN_CODING"), "1"); TOP_BY <- Sys.getenv("TOP_BY", unset = "padj")
INK <- "#1F2933"; RED <- "#D7191C"; BLUE <- "#2B6CB0"
G <- c(Untreated = "Mock_Mock", EBOV = "Mock_EBOV")
PROJ <- Sys.getenv("PROJ", unset = file.path(
  Sys.getenv("HOME"),
  "Documents/aviti_rnaseq/2025/aviti",
  "20251009-100725_Maury_25161_WivfFacICOTHFAUJwSPeHBEFfCfNtCRpUSgVGUkz",
  "results/Maury_25161"))
RNA <- file.path(PROJ, "requant_v2_outputs")
OUT <- file.path(RNA, "figures_ebov_vs_mock"); dir.create(OUT, showWarnings = FALSE)
save2 <- function(p, name, w, h) { ggsave(file.path(OUT, paste0(name, ".png")), p, width = w, height = h, dpi = 300); ggsave(file.path(OUT, paste0(name, ".pdf")), p, width = w, height = h) }

d <- read_csv(file.path(RNA, "DE_v2_EBOV_vs_Mock_annotated.csv"), show_col_types = FALSE) |> filter(!is.na(padj), grepl("^ENSMUSG", gene_id))
if (PROTEIN_CODING) d <- d |> filter(gene_biotype == "protein_coding")
d <- d |> mutate(label = ifelse(is.na(gene_name) | gene_name == "", gene_id, gene_name), y = -log10(pmax(padj, 1e-300)),
                 class = case_when(padj < PADJ_CUT & log2FoldChange >  LFC_CUT ~ "up", padj < PADJ_CUT & log2FoldChange < -LFC_CUT ~ "down", TRUE ~ "other"))
nU <- sum(d$class == "up"); nD <- sum(d$class == "down")
message(sprintf("EBOV vs Mock: %d genes tested%s | DE (padj < %g, |log2FC| > %g): %d up, %d down | padj < %g at any fold-change: %d", nrow(d), if (PROTEIN_CODING) " (protein-coding)" else "",
                PADJ_CUT, LFC_CUT, nU, nD, PADJ_CUT, sum(d$padj < PADJ_CUT)))
write_csv(d |> filter(class != "other") |> arrange(class, padj) |> select(gene_id, gene_name, gene_biotype, baseMean, log2FoldChange, lfcSE, padj, class), file.path(OUT, "EBOV_vs_Mock_DE_genes.csv"))

# ---------------- volcano ----------------
lab <- bind_rows(d |> filter(class == "up") |> slice_max(log2FoldChange, n = TOP_N_LAB, with_ties = FALSE), d |> filter(class == "down") |> slice_min(log2FoldChange, n = TOP_N_LAB, with_ties = FALSE))
xr <- max(abs(d$log2FoldChange)) * 1.05
pv <- ggplot(d, aes(log2FoldChange, y)) +
  geom_point(data = d |> filter(class == "other"), colour = "black", size = 0.6, alpha = 0.4) +
  geom_point(data = d |> filter(class != "other"), aes(colour = class), size = 1, alpha = 0.85) +
  geom_text_repel(data = lab, aes(label = label, colour = class), size = 2.8, fontface = "italic", max.overlaps = Inf, min.segment.length = 0, segment.size = 0.25, segment.colour = "grey60", box.padding = 0.35, seed = 1) +
  geom_hline(yintercept = -log10(PADJ_CUT), linetype = "dashed", colour = "grey50", linewidth = 0.35) +
  geom_vline(xintercept = c(-LFC_CUT, LFC_CUT), linetype = "dashed", colour = "grey50", linewidth = 0.35) +
  scale_colour_manual(values = c(up = RED, down = BLUE), guide = "none") + scale_x_continuous(limits = c(-xr, xr)) +
  annotate("text", x = xr * 0.97, y = Inf, label = sprintf("up: %d", nU), colour = RED, hjust = 1, vjust = 1.6, fontface = "bold", size = 4.2) +
  annotate("text", x = -xr * 0.97, y = Inf, label = sprintf("down: %d", nD), colour = BLUE, hjust = 0, vjust = 1.6, fontface = "bold", size = 4.2) +
  labs(title = "EBOV vs Mock: differential expression", subtitle = sprintf("red = up, blue = down (padj < %g and |log2FC| > %g); black = not differentially expressed; %d genes%s",
                                                                           PADJ_CUT, LFC_CUT, nrow(d), if (PROTEIN_CODING) ", protein-coding" else ""),
       x = "log2 fold-change (EBOV vs Mock)", y = expression(-log[10](padj))) +
  theme_bw(base_size = 12) + theme(panel.grid.minor = element_blank(), panel.grid.major = element_line(colour = "grey94", linewidth = 0.3), panel.border = element_rect(colour = "grey70"),
        axis.text = element_text(colour = INK), axis.title = element_text(colour = INK), plot.title = element_text(colour = INK, face = "bold"), plot.subtitle = element_text(colour = "grey35", size = 9))
save2(pv, "volcano_EBOV_vs_Mock", 8, 5.5)

# ---------------- heatmap of the top genes ----------------
de <- d |> filter(class != "other")
if (!nrow(de)) stop("No differentially expressed genes at these cut-offs")
top <- (if (TOP_BY == "lfc") de |> arrange(desc(abs(log2FoldChange))) else de |> arrange(padj, desc(abs(log2FoldChange)))) |> slice_head(n = TOP_N) |>
  arrange(desc(log2FoldChange)) |> mutate(label = make.unique(label))
write_csv(top |> select(gene_id, gene_name, baseMean, log2FoldChange, padj, class), file.path(OUT, "heatmap_top50_genes.csv"))
vsd <- readRDS(file.path(RNA, "vsd_v2.rds")); vm <- SummarizedExperiment::assay(vsd); grp <- as.character(vsd$group)
cols <- lapply(G, function(g) colnames(vm)[grp == g])
if (any(lengths(cols) == 0)) stop("Missing groups in the VST object: ", paste(names(G)[lengths(cols) == 0], collapse = ", "))
message("Samples used: ", paste(sprintf("%s = %s", names(cols), vapply(cols, paste, "", collapse = ",")), collapse = " | "))
nrep <- lengths(cols)
X <- vm[top$gene_id, unlist(cols), drop = FALSE]; rownames(X) <- top$label
colnames(X) <- paste0(rep(names(cols), nrep), " ", unlist(lapply(nrep, seq_len)))
W <- 0.5 * ncol(X) + 2.5; H <- 0.17 * nrow(X) + 2.2
save_ph <- function(mat, name, pal, brk, at, text, title) {
  ph <- pheatmap(mat, color = pal, breaks = brk, cluster_rows = FALSE, cluster_cols = FALSE, border_color = NA, fontsize_row = 8, fontsize_col = 9, fontsize = 9, angle_col = 45,
                 cellwidth = 26, cellheight = 10, gaps_col = cumsum(nrep)[-length(nrep)], silent = TRUE, legend_breaks = at, legend_labels = text, main = title)
  png(file.path(OUT, paste0(name, ".png")), width = W, height = H, units = "in", res = 300); grid.newpage(); grid.draw(ph$gtable); invisible(dev.off())
  pdf(file.path(OUT, paste0(name, ".pdf")), width = W, height = H); grid.newpage(); grid.draw(ph$gtable); invisible(dev.off())
}
Z <- t(scale(t(X))); lz <- as.numeric(quantile(abs(Z), CAP_Q))
save_ph(Z, "heatmap_top50_EBOV_vs_Mock_rowZ", colorRampPalette(c(BLUE, "white", "#C53030"))(100), seq(-lz, lz, length.out = 101), c(-lz, 0, lz), c(sprintf("%.1f", -lz), "0", sprintf("%.1f", lz)),
        sprintf("Top %d DE genes, EBOV vs Mock: expression, row z-score", nrow(Z)))
Xm <- X - rowMeans(X[, seq_len(nrep[[1]]), drop = FALSE]); lm <- as.numeric(quantile(abs(Xm), CAP_Q))
save_ph(Xm, "heatmap_top50_EBOV_vs_Mock_mockgrey", colorRampPalette(c(BLUE, "#F2F2F2", "#C53030"))(100), seq(-lm, lm, length.out = 101), c(-lm, 0, lm), c(sprintf("%.1f", -lm), "0 (Mock mean)", sprintf("%.1f", lm)),
        sprintf("Top %d DE genes, EBOV vs Mock: each replicate vs Untreated mean (Mock = grey)", nrow(Xm)))
message("Saved to: ", OUT)
