#!/usr/bin/env Rscript
# Volcano plots and heatmaps (protein-coding genes only) for
#   IFNB vs Mock, IFNY vs Mock, and IFNB vs IFNY.
# Needs 07, 08 and 09 to have been run (they create the *_protein_coding.csv
# tables) and the saved vsd_v2.rds from 04.
#
# Per contrast, in <requant_v2_outputs>/figures_protein_coding/:
#   *_volcano_up     genes with log2FC >= 0; red = log2FC > 1 & padj < 0.05;
#                    top 20 by highest log2FC labelled
#   *_volcano_down   genes with log2FC <= 0; red = log2FC < -1 & padj < 0.05;
#                    top 20 by lowest log2FC labelled
#   *_heatmap_all    every protein-coding gene, rows ordered by log2FC
#   *_heatmap_top50_up / _top50_down   top 50 significant up / down genes
# Heatmaps show VST expression z-scored per gene across the samples shown
# (only the two groups being compared). Each figure is written as PNG and PDF.

suppressPackageStartupMessages({
  library(DESeq2); library(readr); library(dplyr); library(ggplot2)
})
for (pkg in c("ggrepel", "pheatmap"))
  if (!requireNamespace(pkg, quietly = TRUE))
    stop("Package '", pkg, "' is needed. Install with: conda install -c conda-forge r-",
         tolower(pkg), "   (or install.packages('", pkg, "'))")

PROJ <- Sys.getenv("PROJ", unset = file.path(
  Sys.getenv("HOME"),
  "Documents/aviti_rnaseq/2025/aviti",
  "20251009-100725_Maury_25161_WivfFacICOTHFAUJwSPeHBEFfCfNtCRpUSgVGUkz",
  "results/Maury_25161"))
OUT  <- file.path(PROJ, "requant_v2_outputs")
FIG  <- file.path(OUT, "figures_protein_coding")
dir.create(FIG, showWarnings = FALSE)

LFC_CUT <- 1; PADJ_CUT <- 0.05; N_LABEL <- 20; N_HEAT <- 50
RED <- "#C53030"; INK <- "#1F2933"

vsd <- readRDS(file.path(OUT, "vsd_v2.rds"))
vmat <- assay(vsd)

# ---- what to plot ------------------------------------------------------------
# groups = (reference, tested); heatmap columns are shown in this order.
specs <- list(
  list(nm = "IFNB_vs_Mock", file = "DE_v2_IFNB_vs_Mock_protein_coding.csv",
       lfc = "log2FoldChange", padj = "padj",
       groups = c("Mock_Mock", "IFNB_Mock"),
       up = "Up in IFNB (ISGs)", down = "Down in IFNB (IDGs)",
       title = "IFNB vs Mock"),
  list(nm = "IFNY_vs_Mock", file = "DE_v2_IFNY_vs_Mock_protein_coding.csv",
       lfc = "log2FoldChange", padj = "padj",
       groups = c("Mock_Mock", "IFNY_Mock"),
       up = "Up in IFNγ (ISGs)", down = "Down in IFNγ (IDGs)",
       title = "IFNγ vs Mock"),
  list(nm = "IFNB_vs_IFNY", file = "DE_v2_IFNB_vs_IFNY_protein_coding.csv",
       lfc = "log2FC_IFNB_vs_IFNY", padj = "padj_diff",
       groups = c("IFNY_Mock", "IFNB_Mock"),
       up = "Higher in IFNB", down = "Higher in IFNγ",
       title = "IFNB vs IFNγ")
)

# ---- helpers -----------------------------------------------------------------
theme_fig <- function() {
  theme_bw(base_size = 12) +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_line(colour = "grey92", linewidth = 0.3),
          panel.border = element_rect(colour = "grey70"),
          axis.text = element_text(colour = INK), axis.title = element_text(colour = INK),
          plot.title = element_text(colour = INK, face = "bold"),
          plot.subtitle = element_text(colour = "grey35"),
          legend.position = "bottom")
}

save_gg <- function(p, base, w, h) {
  ggsave(paste0(base, ".png"), p, width = w, height = h, dpi = 300)
  ggsave(paste0(base, ".pdf"), p, width = w, height = h)
}

volcano <- function(d, direction, s) {
  up <- direction == "up"
  d <- d |> filter(!is.na(padj), if (up) lfc >= 0 else lfc <= 0) |>
    mutate(y = -log10(pmax(padj, 1e-300)),
           hit = padj < PADJ_CUT & (if (up) lfc > LFC_CUT else lfc < -LFC_CUT))
  lab <- d |> filter(hit) |> arrange(if (up) desc(lfc) else lfc) |> head(N_LABEL)
  cls <- if (up) s$up else s$down
  d$class <- factor(ifelse(d$hit, cls, "Other"), levels = c("Other", cls))
  xlab <- paste0("log2 fold-change, ", s$title)
  p <- ggplot(d, aes(lfc, y)) +
    geom_hline(yintercept = -log10(PADJ_CUT), linetype = "dashed", colour = "grey60", linewidth = 0.3) +
    geom_vline(xintercept = if (up) LFC_CUT else -LFC_CUT, linetype = "dashed",
               colour = "grey60", linewidth = 0.3) +
    geom_point(aes(colour = class), size = 0.9, alpha = 0.6) +
    geom_point(data = lab, aes(colour = class), size = 1.4) +
    ggrepel::geom_text_repel(data = lab, aes(label = gene_name), colour = INK, size = 3,
                             max.overlaps = Inf, min.segment.length = 0,
                             segment.colour = "grey60", segment.size = 0.25,
                             box.padding = 0.35, seed = 1) +
    scale_colour_manual(values = setNames(c("black", RED), c("Other", cls)), name = NULL) +
    scale_x_continuous(expand = expansion(mult = if (up) c(0.01, 0.06) else c(0.06, 0.01))) +
    labs(title = paste0(s$title, ": ", if (up) "upregulated" else "downregulated", " protein-coding genes"),
         subtitle = sprintf("%d genes with log2FC %s 0 shown; %d meet log2FC %s %d and padj < %g; top %d labelled",
                            nrow(d), if (up) ">=" else "<=", sum(d$hit),
                            if (up) ">" else "<", if (up) LFC_CUT else -LFC_CUT, PADJ_CUT, min(N_LABEL, nrow(lab))),
         x = xlab, y = expression(-log[10](padj))) +
    theme_fig()
  save_gg(p, file.path(FIG, paste0(s$nm, "_volcano_", direction)), 7, 6.5)
  invisible(nrow(d))
}

group_cols <- c(Mock_Mock = "#8A8A8A", IFNB_Mock = "#2B6CB0", IFNY_Mock = "#D69E2E")
heat_cols  <- colorRampPalette(c("#2C5F9E", "#F7F7F7", "#B2382C"))(101)
heat_breaks <- seq(-2.5, 2.5, length.out = 102)

heatmap_fig <- function(ids, labels, s, tag, w, h, show_names) {
  keep_s <- as.character(vsd$group) %in% s$groups
  cols <- colnames(vmat)[keep_s]
  grp  <- as.character(vsd$group)[keep_s]
  o <- order(match(grp, s$groups)); cols <- cols[o]; grp <- grp[o]
  ids2 <- ids[ids %in% rownames(vmat)]
  labels <- labels[match(ids2, ids)]
  m <- vmat[ids2, cols, drop = FALSE]
  z <- t(scale(t(m)))
  ok <- stats::complete.cases(z)               # drop zero-variance genes
  z <- z[ok, , drop = FALSE]; labels <- labels[ok]
  rownames(z) <- make.unique(ifelse(is.na(labels), rownames(z), labels))
  colnames(z) <- cols
  ann <- data.frame(Group = factor(grp, levels = s$groups), row.names = cols)
  p <- pheatmap::pheatmap(
    z, color = heat_cols, breaks = heat_breaks,
    cluster_rows = FALSE, cluster_cols = FALSE,
    annotation_col = ann, annotation_colors = list(Group = group_cols[s$groups]),
    gaps_col = cumsum(table(factor(grp, levels = s$groups)))[1],
    show_rownames = show_names, show_colnames = TRUE, fontsize = 9, fontsize_row = 7,
    border_color = NA, main = paste0(s$title, " - ", tag, " (", nrow(z), " genes, z-score of VST)"),
    silent = TRUE)
  for (ext in c("png", "pdf")) {
    f <- file.path(FIG, paste0(s$nm, "_heatmap_", gsub("[^A-Za-z0-9]+", "_", tag), ".", ext))
    if (ext == "png") png(f, width = w, height = h, units = "in", res = 300) else pdf(f, width = w, height = h)
    grid::grid.newpage(); grid::grid.draw(p$gtable); dev.off()
  }
  nrow(z)
}

# ---- run ----------------------------------------------------------------------
for (s in specs) {
  raw <- read_csv(file.path(OUT, s$file), show_col_types = FALSE)
  d <- tibble(gene_id = raw$gene_id, gene_name = raw$gene_name,
              lfc = raw[[s$lfc]], padj = raw[[s$padj]]) |> filter(!is.na(lfc))
  message("\n=== ", s$title, " : ", nrow(d), " protein-coding genes ===")

  volcano(d, "up", s); volcano(d, "down", s)

  by_lfc <- d |> arrange(desc(lfc))
  heatmap_fig(by_lfc$gene_id, by_lfc$gene_name, s, "all protein-coding genes", 5.5, 11, FALSE)

  sig <- d |> filter(!is.na(padj), padj < PADJ_CUT)
  top_up   <- sig |> filter(lfc >  LFC_CUT) |> arrange(desc(lfc)) |> head(N_HEAT)
  top_down <- sig |> filter(lfc < -LFC_CUT) |> arrange(lfc)       |> head(N_HEAT)
  message("significant up: ", sum(sig$lfc > LFC_CUT), "   significant down: ", sum(sig$lfc < -LFC_CUT))
  if (nrow(top_up))   heatmap_fig(top_up$gene_id,   top_up$gene_name,   s, "top 50 up",   5.5, 9, TRUE)
  if (nrow(top_down)) heatmap_fig(top_down$gene_id, top_down$gene_name, s, "top 50 down", 5.5, 9, TRUE)
}
message("\nFigures written to: ", FIG)
