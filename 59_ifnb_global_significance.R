#!/usr/bin/env Rscript
# Is the IFNB-induced change to the transcriptome (IFNB vs Mock) statistically significant, in the context of the whole transcriptome?
# Reports, for protein-coding genes (Mock_Mock vs IFNB_Mock samples):
#   1. Per-gene Wald test (DESeq2) with Benjamini-Hochberg padj across all genes tested: number significant / up / down, plus the
#      same with a proper test of |log2FC| > LFC_CUT (results(lfcThreshold, altHypothesis = "greaterAbs")) instead of a post-hoc cutoff
#   2. Storey pi1 = estimated fraction of genes truly changed (1 - pi0, from the p-value distribution); p-value histogram
#   3. Global separation of the samples: exact-permutation PERMANOVA on VST values (pseudo-F, R2) and a t-test on PC1
#   4. Permutation null for the number of DE genes: DESeq2 re-run on every distinct re-labelling of the 6 samples (RUN_PERM)
#   5. Hallmark gene-set context (msigdbr): Fisher over-representation of Hallmark sets among ISGs and IDGs (background = all tested genes),
#      rank-sum test and GSEA on the whole ranked transcriptome (Wald statistic) for Interferon Alpha / Gamma Response
# With 3 vs 3 samples there are only 20 label splits (10 distinct after swapping the group names), so the smallest possible permutation
# p-value is 0.1 - treat 3 and 4 as supporting evidence; the gene-level padj (1) and the gene-set tests (5) carry the statistics.
# Reads requant_v2_outputs/dds_v2_decoyAware_hostNormalised.rds, vsd_v2.rds, DE_v2_IFNB_vs_Mock_protein_coding.csv
# Writes to requant_v2_outputs/figures_significance/: significance_summary.csv (every number reported), per-gene results, hallmark tables,
#   pvalue_histogram, permutation_null, gsea_hallmark_IFN (PNG + PDF)
# Optional packages: msigdbr + clusterProfiler/enrichplot (section 5), qvalue (pi0; a simple estimate is used otherwise)
# Usage: Rscript 59_ifnb_global_significance.R   (PROJ from 'source setup_session.sh')

suppressPackageStartupMessages({ library(readr); library(dplyr); library(tidyr); library(ggplot2); library(DESeq2) })
select <- dplyr::select; filter <- dplyr::filter; rename <- dplyr::rename; slice <- dplyr::slice

PADJ_CUT <- 0.05; LFC_CUT <- 1; RUN_PERM <- TRUE
INK <- "#1F2933"; BLUE <- "#2B6CB0"
PROJ <- Sys.getenv("PROJ", unset = file.path(
  Sys.getenv("HOME"),
  "Documents/aviti_rnaseq/2025/aviti",
  "20251009-100725_Maury_25161_WivfFacICOTHFAUJwSPeHBEFfCfNtCRpUSgVGUkz",
  "results/Maury_25161"))
RNA <- file.path(PROJ, "requant_v2_outputs")
OUT <- file.path(RNA, "figures_significance"); dir.create(OUT, showWarnings = FALSE)
save2 <- function(p, name, w, h) {
  ggsave(file.path(OUT, paste0(name, ".png")), p, width = w, height = h, dpi = 300)
  ggsave(file.path(OUT, paste0(name, ".pdf")), p, width = w, height = h)
}
base_theme <- theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_line(colour = "grey94", linewidth = 0.3), panel.border = element_rect(colour = "grey70"),
        axis.text = element_text(colour = INK), axis.title = element_text(colour = INK), plot.title = element_text(colour = INK, face = "bold"),
        plot.subtitle = element_text(colour = "grey35", size = 9))
S <- list(); add <- function(test, metric, value, note = "") S[[length(S) + 1]] <<- tibble(test = test, metric = metric, value = as.character(value), note = note)

# ---------------- data ----------------
pcg <- read_csv(file.path(RNA, "DE_v2_IFNB_vs_Mock_protein_coding.csv"), show_col_types = FALSE)    # shrunk log2FC, as used in the figures
dds <- readRDS(file.path(RNA, "dds_v2_decoyAware_hostNormalised.rds"))
vsd <- readRDS(file.path(RNA, "vsd_v2.rds"))
dds <- dds[rownames(dds) %in% pcg$gene_id, ]
g <- as.character(dds$group)
idx <- which(g %in% c("Mock_Mock", "IFNB_Mock"))
message(sprintf("Protein-coding genes: %d | samples: Mock %d, IFNB %d", nrow(dds), sum(g == "Mock_Mock"), sum(g == "IFNB_Mock")))
CT <- c("group", "IFNB_Mock", "Mock_Mock")

# ---------------- 1. per-gene Wald test, BH across the transcriptome ----------------
res  <- results(dds, contrast = CT)
resT <- results(dds, contrast = CT, lfcThreshold = LFC_CUT, altHypothesis = "greaterAbs")
R <- as.data.frame(res) |> tibble::rownames_to_column("gene_id") |>
  left_join(pcg |> select(gene_id, gene_name, lfc_shrunk = log2FoldChange, padj_shrunk_file = padj), by = "gene_id") |>
  left_join(as.data.frame(resT) |> tibble::rownames_to_column("gene_id") |> select(gene_id, pvalue_vs_threshold = pvalue, padj_vs_threshold = padj), by = "gene_id")
write_csv(R |> arrange(padj), file.path(OUT, "IFNB_vs_Mock_gene_tests.csv"))
nt <- sum(!is.na(R$padj))
add("1 Wald test", "protein-coding genes tested (non-NA padj)", nt, "BH correction is applied across these genes")
add("1 Wald test", sprintf("genes with padj < %g", PADJ_CUT), sum(R$padj < PADJ_CUT, na.rm = TRUE), sprintf("%.1f%% of tested genes", 100 * mean(R$padj < PADJ_CUT, na.rm = TRUE)))
add("1 Wald test", sprintf("  up (padj < %g, log2FC > 0)", PADJ_CUT), sum(R$padj < PADJ_CUT & R$log2FoldChange > 0, na.rm = TRUE))
add("1 Wald test", sprintf("  down (padj < %g, log2FC < 0)", PADJ_CUT), sum(R$padj < PADJ_CUT & R$log2FoldChange < 0, na.rm = TRUE))
nISG <- sum(R$padj < PADJ_CUT & R$lfc_shrunk > LFC_CUT, na.rm = TRUE); nIDG <- sum(R$padj < PADJ_CUT & R$lfc_shrunk < -LFC_CUT, na.rm = TRUE)
add("1 Wald test", sprintf("ISGs (padj < %g and shrunk log2FC > %g)", PADJ_CUT, LFC_CUT), nISG, "definition used in the figures")
add("1 Wald test", sprintf("IDGs (padj < %g and shrunk log2FC < -%g)", PADJ_CUT, LFC_CUT), nIDG)
add("1 Wald test", sprintf("genes with |log2FC| > %g by test (lfcThreshold, padj < %g)", LFC_CUT, PADJ_CUT), sum(R$padj_vs_threshold < PADJ_CUT, na.rm = TRUE),
    "tests H0: |true log2FC| <= threshold, so the fold-change cutoff itself is statistically supported")
add("1 Wald test", "  of these, up", sum(R$padj_vs_threshold < PADJ_CUT & R$log2FoldChange > 0, na.rm = TRUE))
add("1 Wald test", "  of these, down", sum(R$padj_vs_threshold < PADJ_CUT & R$log2FoldChange < 0, na.rm = TRUE))

# ---------------- 2. Storey pi0 / pi1 and p-value histogram ----------------
pv <- R$pvalue[!is.na(R$padj)]
pi0 <- if (requireNamespace("qvalue", quietly = TRUE)) qvalue::pi0est(pv)$pi0 else min(1, mean(pv > 0.5) / 0.5)
add("2 pi1", "pi0 (fraction of genes not changed)", sprintf("%.3f", pi0), if (requireNamespace("qvalue", quietly = TRUE)) "qvalue::pi0est" else "simple estimate at lambda = 0.5 (install r-qvalue for the smoothed one)")
add("2 pi1", "pi1 = 1 - pi0 (estimated fraction of genes truly changed)", sprintf("%.3f", 1 - pi0), sprintf("about %d of %d genes", round((1 - pi0) * length(pv)), length(pv)))
p_hist <- ggplot(tibble(p = pv), aes(p)) + geom_histogram(breaks = seq(0, 1, 0.025), fill = BLUE, colour = "white", linewidth = 0.2) +
  geom_hline(yintercept = pi0 * length(pv) * 0.025, linetype = "dashed", colour = "grey35") +
  labs(title = "IFNB vs Mock: p-value distribution (protein-coding genes)",
       subtitle = sprintf("spike near 0 = real changes; dashed = flat level expected from unchanged genes (pi0 = %.2f, pi1 = %.2f)", pi0, 1 - pi0),
       x = "Wald p-value", y = "genes") + base_theme
save2(p_hist, "pvalue_histogram", 7, 4.5)

# ---------------- 3. global separation of the samples: exact-permutation PERMANOVA + PC1 ----------------
cn <- colnames(dds)[idx]
X <- t(SummarizedExperiment::assay(vsd)[rownames(dds), cn]); n <- nrow(X); lab0 <- factor(g[idx])
sumsq <- function(M) sum(sweep(M, 2, colMeans(M))^2)
pseudoF <- function(lab) {
  sst <- sumsq(X); ssw <- sum(vapply(split(seq_len(n), lab), function(i) sumsq(X[i, , drop = FALSE]), numeric(1)))
  k <- nlevels(lab); c(F = ((sst - ssw) / (k - 1)) / (ssw / (n - k)), R2 = (sst - ssw) / sst)
}
obs <- pseudoF(lab0)
splits <- combn(n, sum(lab0 == "IFNB_Mock"))
Fall <- apply(splits, 2, function(s) { l <- factor(rep("A", n), levels = c("A", "B")); l[s] <- "B"; pseudoF(l)[["F"]] })
p_perm <- mean(Fall >= obs[["F"]] - 1e-9)
add("3 global separation", "PERMANOVA pseudo-F (Euclidean on VST, all genes)", sprintf("%.2f", obs[["F"]]))
add("3 global separation", "PERMANOVA R2 (fraction of variance explained by IFNB)", sprintf("%.3f", obs[["R2"]]))
add("3 global separation", "PERMANOVA exact permutation p", sprintf("%.3f", p_perm), sprintf("%d label splits; the observed split and its mirror image tie, so the minimum possible p is %.2f", ncol(splits), 2 / ncol(splits)))
pca <- prcomp(X); pc1 <- pca$x[, 1]
tt <- t.test(pc1[lab0 == "IFNB_Mock"], pc1[lab0 == "Mock_Mock"])
add("3 global separation", "PC1 variance explained (Mock vs IFNB samples only)", sprintf("%.1f%%", 100 * pca$sdev[1]^2 / sum(pca$sdev^2)))
add("3 global separation", "PC1 Mock vs IFNB, Welch t-test p", signif(tt$p.value, 3))

# ---------------- 4. permutation null for the number of DE genes ----------------
n_obs <- sum(R$padj < PADJ_CUT & abs(R$log2FoldChange) > LFC_CUT, na.rm = TRUE)
add("4 permutation null", sprintf("observed DE genes (padj < %g, |log2FC| > %g, unshrunk)", PADJ_CUT, LFC_CUT), n_obs)
if (RUN_PERM) {
  message("Permutation null: re-running DESeq2 on re-labelled samples (about a minute per split)...")
  cnt <- round(counts(dds)[, cn]); nf <- normalizationFactors(dds)[, cn]
  key <- function(s) paste(sort(s), collapse = "-")
  true_s <- which(lab0 == "IFNB_Mock")
  canon <- function(s) min(key(s), key(setdiff(seq_len(n), s)))
  ks <- apply(splits, 2, canon); keep <- !duplicated(ks) & ks != canon(true_s)
  nulls <- vapply(which(keep), function(j) {
    l <- factor(rep("A", n), levels = c("A", "B")); l[splits[, j]] <- "B"
    dd <- DESeqDataSetFromMatrix(cnt, data.frame(lab = l, row.names = cn), design = ~ lab)
    normalizationFactors(dd) <- nf
    r <- results(DESeq(dd, quiet = TRUE), contrast = c("lab", "B", "A"))
    sum(r$padj < PADJ_CUT & abs(r$log2FoldChange) > LFC_CUT, na.rm = TRUE)
  }, numeric(1))
  p_null <- (1 + sum(nulls >= n_obs)) / (1 + length(nulls))
  add("4 permutation null", "DE genes under re-labelled samples: median [range]", sprintf("%g [%g - %g]", median(nulls), min(nulls), max(nulls)), sprintf("%d distinct re-labellings", length(nulls)))
  add("4 permutation null", "permutation p (observed count >= null)", sprintf("%.3f", p_null), sprintf("minimum possible with %d permutations is %.2f", length(nulls), 1 / (1 + length(nulls))))
  p_nl <- ggplot(tibble(n = nulls), aes(n)) + geom_dotplot(binwidth = max(1, diff(range(c(nulls, n_obs))) / 60), fill = "grey60", colour = "white") +
    geom_vline(xintercept = n_obs, colour = BLUE, linewidth = 1) +
    annotate("text", x = n_obs, y = 0.95, label = sprintf("observed: %d", n_obs), colour = BLUE, hjust = 1.05, fontface = "bold") +
    scale_y_continuous(NULL, breaks = NULL) +
    labs(title = "Number of DE genes: observed vs re-labelled samples", subtitle = sprintf("grey = %d re-labellings of the Mock / IFNB samples; blue = real labels", length(nulls)),
         x = sprintf("genes with padj < %g and |log2FC| > %g", PADJ_CUT, LFC_CUT)) + base_theme
  save2(p_nl, "permutation_null", 7, 4)
}

# ---------------- 5. Hallmark gene-set context ----------------
have <- requireNamespace("msigdbr", quietly = TRUE)
if (!have) {
  message("msigdbr not installed - skipping section 5 (conda install -c conda-forge -c bioconda r-msigdbr)")
} else {
  msig <- tryCatch(msigdbr::msigdbr(species = "Mus musculus", collection = "H"), error = function(e) msigdbr::msigdbr(species = "Mus musculus", category = "H"))
  sets <- split(msig$gene_symbol, msig$gs_name)
  RS <- R |> filter(!is.na(padj), !is.na(gene_name), gene_name != "") |> distinct(gene_name, .keep_all = TRUE)
  uni <- RS$gene_name
  fisher_tab <- function(genes, label) bind_rows(lapply(names(sets), function(s) {
    st <- intersect(sets[[s]], uni); a <- length(intersect(genes, st))
    m <- matrix(c(a, length(genes) - a, length(st) - a, length(uni) - length(genes) - length(st) + a), 2)
    tibble(gene_set = s, list = label, set_size = length(st), overlap = a, expected = length(genes) * length(st) / length(uni),
           odds_ratio = unname(fisher.test(m, alternative = "greater")$estimate), p = fisher.test(m, alternative = "greater")$p.value)
  })) |> mutate(padj = p.adjust(p, "BH")) |> arrange(p)
  isg <- RS$gene_name[RS$padj < PADJ_CUT & RS$lfc_shrunk >  LFC_CUT]
  idg <- RS$gene_name[RS$padj < PADJ_CUT & RS$lfc_shrunk < -LFC_CUT]
  FT <- bind_rows(fisher_tab(isg, "ISGs"), fisher_tab(idg, "IDGs"))
  write_csv(FT, file.path(OUT, "hallmark_overrepresentation.csv"))
  for (s in c("HALLMARK_INTERFERON_ALPHA_RESPONSE", "HALLMARK_INTERFERON_GAMMA_RESPONSE")) {
    r <- FT |> filter(list == "ISGs", gene_set == s)
    add("5 Hallmark", sprintf("Fisher: %s among ISGs", s), sprintf("%d / %d set genes are ISGs (%.1f expected); OR %.1f; p %.2g; BH padj %.2g", r$overlap, r$set_size, r$expected, r$odds_ratio, r$p, r$padj),
        sprintf("background = %d tested protein-coding genes with a symbol", length(uni)))
  }
  top_idg <- FT |> filter(list == "IDGs") |> slice(1)
  add("5 Hallmark", "most over-represented Hallmark set among IDGs", sprintf("%s (overlap %d, BH padj %.2g)", top_idg$gene_set, top_idg$overlap, top_idg$padj))
  rk <- RS |> filter(!is.na(stat)) |> arrange(desc(stat)); rkv <- setNames(rk$stat, rk$gene_name)
  for (s in c("HALLMARK_INTERFERON_ALPHA_RESPONSE", "HALLMARK_INTERFERON_GAMMA_RESPONSE")) {
    inn <- rk$gene_name %in% sets[[s]]
    w <- wilcox.test(rk$stat[inn], rk$stat[!inn])
    add("5 Hallmark", sprintf("rank-sum: Wald statistic of %s genes vs all other genes", s), sprintf("median %.1f vs %.1f; p %.2g", median(rk$stat[inn]), median(rk$stat[!inn]), w$p.value),
        "whole-transcriptome shift of the set, no cutoff")
  }
  if (requireNamespace("clusterProfiler", quietly = TRUE)) {
    gs <- tryCatch(clusterProfiler::GSEA(rkv, TERM2GENE = msig[, c("gs_name", "gene_symbol")], pvalueCutoff = 1, minGSSize = 10, maxGSSize = 500, eps = 0, seed = TRUE, verbose = FALSE), error = function(e) { message("GSEA failed: ", conditionMessage(e)); NULL })
    if (!is.null(gs)) {
      write_csv(as.data.frame(gs@result) |> select(-core_enrichment), file.path(OUT, "hallmark_gsea.csv"))
      write_csv(as.data.frame(gs@result) |> select(ID, setSize, enrichmentScore, NES, pvalue, p.adjust, core_enrichment), file.path(OUT, "hallmark_gsea_with_leading_edge.csv"))
      for (s in c("HALLMARK_INTERFERON_ALPHA_RESPONSE", "HALLMARK_INTERFERON_GAMMA_RESPONSE")) {
        r <- gs@result[gs@result$ID == s, ]
        if (nrow(r)) add("5 Hallmark", sprintf("GSEA: %s", s), sprintf("NES %.2f; p %.2g; BH padj %.2g; %d genes", r$NES, r$pvalue, r$p.adjust, r$setSize), "ranked by Wald statistic, whole transcriptome")
      }
      sig_n <- sum(gs@result$p.adjust < PADJ_CUT); add("5 Hallmark", sprintf("Hallmark sets with GSEA padj < %g", PADJ_CUT), sig_n, sprintf("of %d tested", nrow(gs@result)))
      ids <- intersect(c("HALLMARK_INTERFERON_ALPHA_RESPONSE", "HALLMARK_INTERFERON_GAMMA_RESPONSE"), gs@result$ID)
      if (length(ids) && requireNamespace("enrichplot", quietly = TRUE)) tryCatch(save2(enrichplot::gseaplot2(gs, geneSetID = ids, pvalue_table = TRUE, title = "Hallmark interferon responses, IFNB vs Mock (whole transcriptome ranked by Wald statistic)"), "gsea_hallmark_IFN", 8, 6),
                                                                                  error = function(e) message("GSEA plot failed: ", conditionMessage(e)))
    }
  } else message("clusterProfiler not installed - GSEA skipped")
}

# ---------------- summary ----------------
SUM <- bind_rows(S)
write_csv(SUM, file.path(OUT, "significance_summary.csv"))
message("\n================ IFNB vs Mock: significance summary ================")
for (i in seq_len(nrow(SUM))) message(sprintf("[%s] %s: %s%s", SUM$test[i], SUM$metric[i], SUM$value[i], if (nzchar(SUM$note[i])) paste0("   (", SUM$note[i], ")") else ""))
message("\nSaved to: ", OUT)
