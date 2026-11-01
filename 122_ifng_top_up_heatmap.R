#!/usr/bin/env Rscript
# Heatmap of the top TOP_N (default 50) UPREGULATED protein-coding genes in IFNG (IFNY) vs Mock, rows ordered by log2FC (highest first), columns Untreated | IFNG | EBOV | IFNG+EBOV (every replicate; EBOV groups via EXTRA), not clustered. Two versions:
#   heatmap_top50_IFNG_vs_Mock_rowZ        row z-score of VST expression (blue below the gene's mean across the 6 samples, red above)
#   heatmap_top50_IFNG_vs_Mock_mockgrey    each replicate's VST minus the gene's Untreated mean: Mock level = light grey, higher = red, lower = blue
# Candidates: protein-coding, padj < PADJ_CUT (0.05) and shrunk log2FC > LFC_CUT (1). Ranked by the largest shrunk log2FC (default TOP_BY=lfc) or by smallest padj (TOP_BY=padj). Rows are genes (gene-level counts, not individual transcript isoforms).
# Env: EXTRA (default Mock_EBOV,IFNY_EBOV = EBOV and IFNG+EBOV columns added after Untreated and IFNG; EXTRA=none for the two-group version; the genes are always the top IFNG-vs-Mock ones), TOP_N (50), TOP_BY (lfc|padj), LFC_CUT, PADJ_CUT, CONTRAST (default IFNY_vs_Mock) and GROUP (default IFNY_Mock; the treated group name in the VST object), e.g. CONTRAST=IFNB_vs_Mock GROUP=IFNB_Mock for IFNB.
# Reads DE_v2_<CONTRAST>_annotated.csv and vsd_v2.rds. Writes to requant_v2_outputs/figures_top_up_<name>/ (PNG + PDF) plus heatmap_top50_genes.csv.   Usage: Rscript 122_ifng_top_up_heatmap.R   (PROJ from 'source setup_session.sh')

suppressPackageStartupMessages({ library(readr); library(dplyr); library(pheatmap); library(DESeq2); library(grid) })
select <- dplyr::select; filter <- dplyr::filter

TOP_N <- as.integer(Sys.getenv("TOP_N", "50")); TOP_BY <- Sys.getenv("TOP_BY", "lfc"); PADJ_CUT <- as.numeric(Sys.getenv("PADJ_CUT", "0.05")); LFC_CUT <- as.numeric(Sys.getenv("LFC_CUT", "1")); CAP_Q <- 0.98
CONTRAST <- Sys.getenv("CONTRAST", "IFNY_vs_Mock"); GROUP <- Sys.getenv("GROUP", "IFNY_Mock"); NAME <- sub("_vs_Mock", "", CONTRAST); NAME <- ifelse(NAME == "IFNY", "IFNG", NAME); BLUE <- "#2B6CB0"
LAB <- c(Mock_EBOV = "EBOV", IFNY_EBOV = "IFNG+EBOV", IFNB_EBOV = "IFNB+EBOV", Mock_MARV = "MARV", IFNY_MARV = "IFNG+MARV", IFNB_MARV = "IFNB+MARV", IFNB_Mock = "IFNB", IFNY_Mock = "IFNG")
EXTRA <- Sys.getenv("EXTRA", "Mock_EBOV,IFNY_EBOV"); EXTRA <- if (EXTRA %in% c("", "none")) character(0) else trimws(strsplit(EXTRA, ",")[[1]])   # extra groups drawn after Untreated and the treated group; EXTRA=none for the two-group version
G <- c(Untreated = "Mock_Mock", setNames(GROUP, NAME), setNames(EXTRA, ifelse(EXTRA %in% names(LAB), LAB[EXTRA], EXTRA))); SUF <- if (length(EXTRA)) paste0("_plus_", paste(names(G)[-(1:2)], collapse = "_")) else ""
SUF <- gsub("\\+", "", SUF)
PROJ <- Sys.getenv("PROJ", unset = file.path(
  Sys.getenv("HOME"),
  "Documents/aviti_rnaseq/2025/aviti",
  "20251009-100725_Maury_25161_WivfFacICOTHFAUJwSPeHBEFfCfNtCRpUSgVGUkz",
  "results/Maury_25161"))
RNA <- file.path(PROJ, "requant_v2_outputs"); OUT <- file.path(RNA, paste0("figures_top_up_", NAME)); dir.create(OUT, showWarnings = FALSE)

d <- read_csv(file.path(RNA, paste0("DE_v2_", CONTRAST, "_annotated.csv")), show_col_types = FALSE) |> filter(!is.na(padj), grepl("^ENSMUSG", gene_id), gene_biotype == "protein_coding")
d <- d |> mutate(label = ifelse(is.na(gene_name) | gene_name == "", gene_id, gene_name))
up <- d |> filter(padj < PADJ_CUT, log2FoldChange > LFC_CUT); message(nrow(up), " protein-coding genes are up (padj < ", PADJ_CUT, ", log2FC > ", LFC_CUT, ") in ", CONTRAST)
if (!nrow(up)) stop("No upregulated protein-coding genes at these cut-offs")
top <- (if (TOP_BY == "padj") up |> arrange(padj, desc(log2FoldChange)) else up |> arrange(desc(log2FoldChange))) |> slice_head(n = TOP_N) |> arrange(desc(log2FoldChange)) |> mutate(label = make.unique(label))
write_csv(top |> select(gene_id, gene_name, baseMean, log2FoldChange, lfcSE, padj), file.path(OUT, "heatmap_top50_genes.csv"))
vsd <- readRDS(file.path(RNA, "vsd_v2.rds")); vm <- SummarizedExperiment::assay(vsd); grp <- as.character(vsd$group)
cols <- lapply(G, function(g) colnames(vm)[grp == g])
if (any(lengths(cols) == 0)) stop("Missing groups in the VST object: ", paste(names(G)[lengths(cols) == 0], collapse = ", "), " (groups present: ", paste(unique(grp), collapse = ", "), ")")
message("Samples used: ", paste(sprintf("%s = %s", names(cols), vapply(cols, paste, "", collapse = ",")), collapse = " | "))
nrep <- lengths(cols); X <- vm[top$gene_id, unlist(cols), drop = FALSE]; rownames(X) <- top$label
colnames(X) <- paste0(rep(names(cols), nrep), " ", unlist(lapply(nrep, seq_len)))
W <- 0.5 * ncol(X) + 2.5; H <- 0.17 * nrow(X) + 2.2
save_ph <- function(mat, name, pal, brk, at, text, title) {
  ph <- pheatmap(mat, color = pal, breaks = brk, cluster_rows = FALSE, cluster_cols = FALSE, border_color = NA, fontsize_row = 8, fontsize_col = 9, fontsize = 9, angle_col = 45,
                 cellwidth = 26, cellheight = 10, gaps_col = cumsum(nrep)[-length(nrep)], silent = TRUE, legend_breaks = at, legend_labels = text, main = title)
  png(file.path(OUT, paste0(name, ".png")), width = W, height = H, units = "in", res = 300); grid.newpage(); grid.draw(ph$gtable); invisible(dev.off())
  pdf(file.path(OUT, paste0(name, ".pdf")), width = W, height = H); grid.newpage(); grid.draw(ph$gtable); invisible(dev.off())
}
Z <- t(scale(t(X))); lz <- as.numeric(quantile(abs(Z), CAP_Q))
save_ph(Z, paste0("heatmap_top", nrow(Z), "_", NAME, "_vs_Mock", SUF, "_rowZ"), colorRampPalette(c(BLUE, "white", "#C53030"))(100), seq(-lz, lz, length.out = 101), c(-lz, 0, lz), c(sprintf("%.1f", -lz), "0", sprintf("%.1f", lz)),
        sprintf("Top %d upregulated protein-coding genes (%s vs Mock), all conditions: row z-score", nrow(Z), NAME))
Xm <- X - rowMeans(X[, seq_len(nrep[[1]]), drop = FALSE]); lm <- as.numeric(quantile(abs(Xm), CAP_Q))
save_ph(Xm, paste0("heatmap_top", nrow(Xm), "_", NAME, "_vs_Mock", SUF, "_mockgrey"), colorRampPalette(c(BLUE, "#F2F2F2", "#C53030"))(100), seq(-lm, lm, length.out = 101), c(-lm, 0, lm), c(sprintf("%.1f", -lm), "0 (Mock mean)", sprintf("%.1f", lm)),
        sprintf("Top %d upregulated genes (%s vs Mock), all conditions: each replicate vs Untreated mean (Mock = grey)", nrow(Xm), NAME))
message("Saved to: ", OUT)
