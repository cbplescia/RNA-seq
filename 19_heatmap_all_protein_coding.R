#!/usr/bin/env Rscript
# One heatmap of all protein-coding genes: IFNB samples on the left, Mock in the
# centre, IFNY on the right.
#
# Values are VST expression z-scored per gene across the 9 samples shown. Rows are
# grouped into K k-means clusters; clusters run top to bottom from "highest in IFN
# vs Mock" to "lowest", and genes are hierarchically ordered inside each cluster.
# (Clustering all ~12,000 genes hierarchically in one go would be very slow, and a
# single log2FC sort would favour one interferon over the other.)
# Needs 07_add_gene_symbols.R and 08_split_protein_coding.R to have been run, and
# the pheatmap package.

suppressPackageStartupMessages({ library(DESeq2); library(readr); library(dplyr) })
if (!requireNamespace("pheatmap", quietly = TRUE)) stop("Install pheatmap (conda install -c conda-forge r-pheatmap)")

K <- 10
PROJ <- Sys.getenv("PROJ", unset = file.path(
  Sys.getenv("HOME"),
  "Documents/aviti_rnaseq/2025/aviti",
  "20251009-100725_Maury_25161_WivfFacICOTHFAUJwSPeHBEFfCfNtCRpUSgVGUkz",
  "results/Maury_25161"))
OUT <- file.path(PROJ, "requant_v2_outputs")
FIG <- file.path(OUT, "figures_protein_coding")
dir.create(FIG, showWarnings = FALSE)

group_order <- c("IFNB_Mock", "Mock_Mock", "IFNY_Mock")      # left, centre, right
group_cols  <- c(IFNB_Mock = "#2B6CB0", Mock_Mock = "#8A8A8A", IFNY_Mock = "#D69E2E")

vsd  <- readRDS(file.path(OUT, "vsd_v2.rds"))
vmat <- assay(vsd)
pc   <- read_csv(file.path(OUT, "DE_v2_IFNB_vs_Mock_protein_coding.csv"), show_col_types = FALSE)
ids  <- intersect(pc$gene_id, rownames(vmat))

grp  <- as.character(vsd$group); names(grp) <- colnames(vmat)
cols <- unlist(lapply(group_order, function(g) names(grp)[grp == g]))
gcol <- grp[cols]
message("Samples (left to right): ", paste(paste0(cols, " [", gcol, "]"), collapse = ", "))

z <- t(scale(t(vmat[ids, cols])))
z <- z[stats::complete.cases(z), , drop = FALSE]          # drops zero-variance genes
message("Protein-coding genes in heatmap: ", nrow(z))

# ---- cluster and order rows ---------------------------------------------------------
set.seed(1)
km <- kmeans(z, centers = K, nstart = 10, iter.max = 100)
gm <- sapply(group_order, function(g) rowMeans(z[, gcol == g, drop = FALSE]))   # genes x 3
score <- sapply(seq_len(K), function(k) {
  r <- km$cluster == k
  mean((gm[r, "IFNB_Mock"] + gm[r, "IFNY_Mock"]) / 2 - gm[r, "Mock_Mock"])
})
new_id <- match(seq_len(K), order(score, decreasing = TRUE))     # 1 = most induced by IFN
cl <- new_id[km$cluster]; names(cl) <- rownames(z)

ord <- unlist(lapply(seq_len(K), function(k) {
  g <- names(cl)[cl == k]
  if (length(g) > 2) g[hclust(dist(z[g, , drop = FALSE]), method = "ward.D2")$order] else g
}))
z <- z[ord, , drop = FALSE]
cl_f <- factor(paste0("C", cl[ord]), levels = paste0("C", seq_len(K)))

# cluster summary + gene table
summ <- data.frame(cluster = levels(cl_f), n = as.integer(table(cl_f)),
                   IFNB = tapply(gm[ord, "IFNB_Mock"], cl_f, mean),
                   Mock = tapply(gm[ord, "Mock_Mock"], cl_f, mean),
                   IFNY = tapply(gm[ord, "IFNY_Mock"], cl_f, mean))
message("\nCluster summary (mean z-score per group):")
print(format(summ, digits = 2), row.names = FALSE)
write_csv(tibble(gene_id = rownames(z), cluster = as.character(cl_f),
                 mean_z_IFNB = gm[rownames(z), "IFNB_Mock"], mean_z_Mock = gm[rownames(z), "Mock_Mock"],
                 mean_z_IFNY = gm[rownames(z), "IFNY_Mock"]) |>
            left_join(pc |> select(gene_id, gene_name), by = "gene_id") |>
            relocate(gene_name, .after = gene_id),
          file.path(FIG, "heatmap_all_protein_coding_clusters.csv"))

# ---- draw -----------------------------------------------------------------------
heat_cols   <- colorRampPalette(c("#2C5F9E", "#F7F7F7", "#B2382C"))(101)
heat_breaks <- seq(-2.5, 2.5, length.out = 102)
ann_col <- data.frame(Group = factor(gcol, levels = group_order), row.names = cols)
ann_row <- data.frame(Cluster = cl_f, row.names = rownames(z))
cl_cols <- setNames(grDevices::hcl.colors(K, "Dark 3"), levels(cl_f))

p <- pheatmap::pheatmap(
  z, color = heat_cols, breaks = heat_breaks,
  cluster_rows = FALSE, cluster_cols = FALSE,
  annotation_col = ann_col, annotation_row = ann_row,
  annotation_colors = list(Group = group_cols, Cluster = cl_cols),
  gaps_col = unname(cumsum(table(factor(gcol, levels = group_order)))[1:2]),
  gaps_row = unname(cumsum(table(cl_f))[-K]),
  show_rownames = FALSE, show_colnames = TRUE, fontsize = 9, border_color = NA,
  main = sprintf("Protein-coding genes (%d): IFNB | Mock | IFNy  (z-score of VST)", nrow(z)),
  silent = TRUE)
for (ext in c("png", "pdf")) {
  f <- file.path(FIG, paste0("heatmap_all_protein_coding_IFNB_Mock_IFNY.", ext))
  if (ext == "png") png(f, width = 7.5, height = 11, units = "in", res = 300) else pdf(f, width = 7.5, height = 11)
  grid::grid.newpage(); grid::grid.draw(p$gtable); dev.off()
}
message("\nSaved to: ", FIG)
