#!/usr/bin/env Rscript
# Deep dive into IFNB-response patterns and how EBOV / MARV change them (protein-coding genes).
# Rule-based, so the groups do not depend on k-means cluster numbers.
#
# Contrasts (DESeq2, shrunk log2FC + padj, from dds_v2_decoyAware_hostNormalised.rds):
#   B  = IFNB vs Mock        E = EBOV vs Mock        M = MARV vs Mock
#   BE = IFNB+EBOV vs Mock   BM = IFNB+MARV vs Mock
#   BE_B = IFNB+EBOV vs IFNB   BM_B = IFNB+MARV vs IFNB   (what the virus does ON TOP of IFNB)
# retention  = log2FC(IFNB+virus vs Mock) / log2FC(IFNB vs Mock): 1 = full IFNB response kept, 0 = response gone,
#              negative = reversed.
# IFNB responder: |log2FC_B| > LFC_CUT and padj_B < PADJ_CUT.
#   "suppressed"      retention <= R_STRONG (default 0.25, i.e. >= 75% of the response lost) AND the virus-on-IFNB
#                     contrast moves the gene back toward Mock with padj < PADJ_CUT
#   "fails to suppress" retention >= R_FAIL (default 0.75, i.e. >= 75% of the response kept)
#   "IFNB does not alter" |log2FC_B| < LFC_NOT (default 0.5)
# Patterns written (CSV + labelled heatmap + enrichment):
#   P1  IFNB-UP genes suppressed by BOTH EBOV and MARV
#   P2  IFNB-DOWN genes whose repression is relieved by BOTH viruses
#   P3u/P3d  IFNB-up / IFNB-down genes EBOV FAILS to suppress        P4u/P4d  same for MARV
#   P5  genes EBOV alters (|log2FC_E| > LFC_CUT, padj < PADJ_CUT) that IFNB does not
#   P6  same for MARV         P7  genes BOTH viruses alter in the same direction that IFNB does not
# plus master_table.csv with every gene and every statistic, so you can try other thresholds.
# Heatmap columns: Mock | IFNB | IFNB+EBOV | EBOV | IFNB+MARV | MARV; colour = log2 change of VST vs Mock mean.
# Enrichment (GO BP, KEGG) uses clusterProfiler if installed; background = all protein-coding genes tested.
# Usage: Rscript 48_patterns_virus_vs_ifnb.R   (PROJ from 'source setup_session.sh')
# Writes to <project>/requant_v2_outputs/figures_patterns/

suppressPackageStartupMessages({ library(DESeq2); library(readr); library(dplyr); library(tidyr) })
if (!requireNamespace("pheatmap", quietly = TRUE)) stop("Install pheatmap")
select <- dplyr::select; first <- dplyr::first; rename <- dplyr::rename
slice <- dplyr::slice; filter <- dplyr::filter; desc <- dplyr::desc; count <- dplyr::count

LFC_CUT <- 1; PADJ_CUT <- 0.05; LFC_NOT <- 0.5; R_STRONG <- 0.25; R_FAIL <- 0.75
MAX_LABEL <- 60       # label genes on a heatmap only when it has at most this many rows
MIN_ENRICH <- 10      # run enrichment only for patterns with at least this many genes

SETS <- list(Mock = c("1", "2", "3"), IFNB = c("10", "11", "12"), `IFNB+EBOV` = c("16", "17", "18"),
             EBOV = c("4", "5", "6"), `IFNB+MARV` = c("22", "23", "24"), MARV = c("7", "8", "9"))
EXPECTED <- c(Mock = "Mock_Mock", IFNB = "IFNB_Mock", `IFNB+EBOV` = "IFNB_EBOV", EBOV = "Mock_EBOV",
              `IFNB+MARV` = "IFNB_MARV", MARV = "Mock_MARV")
GCOL <- c(Mock = "#7F7F7F", IFNB = "#2B6CB0", `IFNB+EBOV` = "#9BD08E", EBOV = "#1B7837", `IFNB+MARV` = "#C9A6DB", MARV = "#762A83")

PROJ <- Sys.getenv("PROJ", unset = file.path(
  Sys.getenv("HOME"),
  "Documents/aviti_rnaseq/2025/aviti",
  "20251009-100725_Maury_25161_WivfFacICOTHFAUJwSPeHBEFfCfNtCRpUSgVGUkz",
  "results/Maury_25161"))
RNA <- file.path(PROJ, "requant_v2_outputs")
OUT <- file.path(RNA, "figures_patterns"); dir.create(OUT, showWarnings = FALSE)

# ---- data -------------------------------------------------------------------------------------------------------------
dds <- readRDS(file.path(RNA, "dds_v2_decoyAware_hostNormalised.rds"))
vsd <- readRDS(file.path(RNA, "vsd_v2.rds")); vm <- assay(vsd)
pc <- read_csv(file.path(RNA, "DE_v2_IFNB_vs_Mock_protein_coding.csv"), show_col_types = FALSE)
dds <- dds[rownames(dds) %in% pc$gene_id, ]
message("Protein-coding genes in the model: ", nrow(dds))
need <- c("Mock_Mock", "IFNB_Mock", "Mock_EBOV", "Mock_MARV", "IFNB_EBOV", "IFNB_MARV")
if (!all(need %in% levels(dds$group))) stop("Missing groups in dds: ", paste(setdiff(need, levels(dds$group)), collapse = ", "))

CONTR <- list(B = c("IFNB_Mock", "Mock_Mock"), E = c("Mock_EBOV", "Mock_Mock"), M = c("Mock_MARV", "Mock_Mock"),
              BE = c("IFNB_EBOV", "Mock_Mock"), BM = c("IFNB_MARV", "Mock_Mock"),
              BE_B = c("IFNB_EBOV", "IFNB_Mock"), BM_B = c("IFNB_MARV", "IFNB_Mock"))
use_ashr <- requireNamespace("ashr", quietly = TRUE)
get_contrast <- function(nm) {
  ct <- c("group", CONTR[[nm]])
  res <- results(dds, contrast = ct)
  shr <- if (use_ashr) lfcShrink(dds, contrast = ct, res = res, type = "ashr") else lfcShrink(dds, contrast = ct, type = "normal")
  out <- tibble(gene_id = rownames(res), lfc = shr$log2FoldChange, padj = res$padj)
  names(out)[2:3] <- paste0(c("lfc_", "padj_"), nm); out
}
message("Computing ", length(CONTR), " contrasts (", if (use_ashr) "ashr" else "normal", " shrinkage)...")
M <- Reduce(function(a, b) left_join(a, b, by = "gene_id"), lapply(names(CONTR), get_contrast))
M <- M |> left_join(pc |> select(gene_id, gene_name), by = "gene_id") |> relocate(gene_name, .after = gene_id) |>
  mutate(retain_E = ifelse(abs(lfc_B) > 0, lfc_BE / lfc_B, NA_real_), retain_M = ifelse(abs(lfc_B) > 0, lfc_BM / lfc_B, NA_real_))
write_csv(M, file.path(OUT, "master_table.csv"))

# ---- pattern definitions ----------------------------------------------------------------------------------------------
resp_up <- !is.na(M$padj_B) & M$padj_B < PADJ_CUT & M$lfc_B > LFC_CUT
resp_dn <- !is.na(M$padj_B) & M$padj_B < PADJ_CUT & M$lfc_B < -LFC_CUT
supp <- function(retain, lfc_on, padj_on) !is.na(retain) & retain <= R_STRONG & sign(lfc_on) == -sign(M$lfc_B) & !is.na(padj_on) & padj_on < PADJ_CUT
sup_E <- supp(M$retain_E, M$lfc_BE_B, M$padj_BE_B); sup_M <- supp(M$retain_M, M$lfc_BM_B, M$padj_BM_B)
fail_E <- !is.na(M$retain_E) & M$retain_E >= R_FAIL; fail_M <- !is.na(M$retain_M) & M$retain_M >= R_FAIL
not_B <- abs(M$lfc_B) < LFC_NOT
alt_E <- !is.na(M$padj_E) & M$padj_E < PADJ_CUT & abs(M$lfc_E) > LFC_CUT
alt_M <- !is.na(M$padj_M) & M$padj_M < PADJ_CUT & abs(M$lfc_M) > LFC_CUT

PATTERNS <- list(
  list(id = "P1_IFNB_up_suppressed_by_EBOV_and_MARV", title = "IFNB-induced genes suppressed by both EBOV and MARV",
       sel = resp_up & sup_E & sup_M, order = quote(-lfc_B), lab = "lfc_B", labname = "IFNB log2FC"),
  list(id = "P2_IFNB_down_relieved_by_EBOV_and_MARV", title = "IFNB-repressed genes relieved by both EBOV and MARV",
       sel = resp_dn & sup_E & sup_M, order = quote(lfc_B), lab = "lfc_B", labname = "IFNB log2FC"),
  list(id = "P3u_IFNB_up_EBOV_fails_to_suppress", title = "IFNB-induced genes EBOV fails to suppress",
       sel = resp_up & fail_E, order = quote(-lfc_B), lab = "lfc_B", labname = "IFNB log2FC"),
  list(id = "P3d_IFNB_down_EBOV_fails_to_relieve", title = "IFNB-repressed genes EBOV fails to relieve",
       sel = resp_dn & fail_E, order = quote(lfc_B), lab = "lfc_B", labname = "IFNB log2FC"),
  list(id = "P4u_IFNB_up_MARV_fails_to_suppress", title = "IFNB-induced genes MARV fails to suppress",
       sel = resp_up & fail_M, order = quote(-lfc_B), lab = "lfc_B", labname = "IFNB log2FC"),
  list(id = "P4d_IFNB_down_MARV_fails_to_relieve", title = "IFNB-repressed genes MARV fails to relieve",
       sel = resp_dn & fail_M, order = quote(lfc_B), lab = "lfc_B", labname = "IFNB log2FC"),
  list(id = "P5_altered_by_EBOV_not_IFNB", title = "Genes EBOV alters that IFNB does not",
       sel = alt_E & not_B, order = quote(-lfc_E), lab = "lfc_E", labname = "EBOV log2FC"),
  list(id = "P6_altered_by_MARV_not_IFNB", title = "Genes MARV alters that IFNB does not",
       sel = alt_M & not_B, order = quote(-lfc_M), lab = "lfc_M", labname = "MARV log2FC"),
  list(id = "P7_altered_by_both_viruses_not_IFNB", title = "Genes both EBOV and MARV alter (same direction) that IFNB does not",
       sel = alt_E & alt_M & not_B & sign(M$lfc_E) == sign(M$lfc_M), order = quote(-(lfc_E + lfc_M) / 2), lab = "lfc_E", labname = "EBOV log2FC"))

message("\nPattern sizes (thresholds: |log2FC| > ", LFC_CUT, ", padj < ", PADJ_CUT, ", strong suppression retention <= ", R_STRONG,
        ", failure retention >= ", R_FAIL, "):")
message(sprintf("  IFNB responders: %d up, %d down", sum(resp_up), sum(resp_dn)))
message(sprintf("  of the %d up: suppressed by EBOV %d, by MARV %d, by both %d | fail to suppress: EBOV %d, MARV %d, both %d",
                sum(resp_up), sum(resp_up & sup_E), sum(resp_up & sup_M), sum(resp_up & sup_E & sup_M),
                sum(resp_up & fail_E), sum(resp_up & fail_M), sum(resp_up & fail_E & fail_M)))
message(sprintf("  of the %d down: relieved by EBOV %d, by MARV %d, by both %d | fail to relieve: EBOV %d, MARV %d, both %d",
                sum(resp_dn), sum(resp_dn & sup_E), sum(resp_dn & sup_M), sum(resp_dn & sup_E & sup_M),
                sum(resp_dn & fail_E), sum(resp_dn & fail_M), sum(resp_dn & fail_E & fail_M)))

# ---- heatmap + enrichment helpers -----------------------------------------------------------------------------------------
sel_s <- tibble(sample = unlist(SETS), group = rep(names(SETS), lengths(SETS)))
sel_s$data_group <- as.character(vsd$group)[match(sel_s$sample, colnames(vm))]
if (any(is.na(sel_s$data_group))) stop("Samples not in the VST object: ", paste(sel_s$sample[is.na(sel_s$data_group)], collapse = ", "))
if (any(sel_s$data_group != EXPECTED[sel_s$group])) message("\nWARNING: some chosen samples are recorded under a different group than assumed (check the group table in the data). Using your sample numbers.")
X <- vm[, sel_s$sample, drop = FALSE]
mock_cols <- sel_s$sample[sel_s$group == "Mock"]
ann_col <- data.frame(Group = factor(sel_s$group, levels = names(SETS)), row.names = sel_s$sample)
gaps <- unname(cumsum(lengths(SETS))[-length(SETS)])
pal <- colorRampPalette(c("#2C5F9E", "#F7F7F7", "#B2382C"))(101)

draw_heatmap <- function(g, p) {
  g <- g[g$gene_id %in% rownames(X), ]
  m <- X[g$gene_id, , drop = FALSE] - rowMeans(X[g$gene_id, mock_cols, drop = FALSE])
  cap <- min(8, max(2, ceiling(quantile(abs(m), 0.98, na.rm = TRUE))))
  lab <- nrow(m) <= MAX_LABEL
  args <- list(m, color = pal, breaks = seq(-cap, cap, length.out = 102), cluster_rows = FALSE, cluster_cols = FALSE,
               annotation_col = ann_col, annotation_colors = list(Group = GCOL), gaps_col = gaps, show_colnames = TRUE,
               border_color = NA, silent = TRUE, fontsize = 9, show_rownames = lab,
               main = sprintf("%s (%d genes; colour = log2 change vs Mock mean)", p$title, nrow(m)))
  if (lab) { args$labels_row <- make.unique(sprintf("%s  (%+.1f)", g$gene_name, g[[p$lab]])); args$fontsize_row <- 8 }
  h <- do.call(pheatmap::pheatmap, args)
  ht <- if (lab) 0.2 * nrow(m) + 3.2 else 8
  for (ext in c("png", "pdf")) {
    f <- file.path(OUT, paste0("heatmap_", p$id, ".", ext))
    if (ext == "png") png(f, width = 8, height = ht, units = "in", res = 300) else pdf(f, width = 8, height = ht)
    grid::grid.newpage(); grid::grid.draw(h$gtable); dev.off()
  }
}

have_cp <- requireNamespace("clusterProfiler", quietly = TRUE) && requireNamespace("org.Mm.eg.db", quietly = TRUE) &&
  requireNamespace("enrichplot", quietly = TRUE)
if (!have_cp) message("\nclusterProfiler / enrichplot / org.Mm.eg.db not all installed - enrichment skipped")
uni_ens <- M$gene_id
enrich <- function(ids, id) {
  if (!have_cp || length(ids) < MIN_ENRICH) return(invisible())
  map <- suppressMessages(suppressWarnings(clusterProfiler::bitr(unique(c(ids, uni_ens)), "ENSEMBL", "ENTREZID", org.Mm.eg.db::org.Mm.eg.db))) |>
    distinct(ENSEMBL, .keep_all = TRUE)
  g <- unique(na.omit(map$ENTREZID[match(ids, map$ENSEMBL)])); uni <- unique(na.omit(map$ENTREZID[match(uni_ens, map$ENSEMBL)]))
  run <- function(label, fn) {
    r <- tryCatch(fn(), error = function(e) { message("  ", id, " ", label, ": not run (", conditionMessage(e), ")"); NULL })
    if (is.null(r) || !nrow(r@result)) return(invisible())
    r <- tryCatch(clusterProfiler::setReadable(r, org.Mm.eg.db::org.Mm.eg.db, keyType = "ENTREZID"), error = function(e) r)
    write_csv(r@result, file.path(OUT, paste0("enrichment_", id, "_", label, ".csv")))
    sig <- sum(r@result$p.adjust < 0.05)
    message(sprintf("  %s %s: %d terms with p.adjust < 0.05%s", id, label, sig,
                    if (sig) paste0("  (", paste(head(r@result$Description[r@result$p.adjust < 0.05], 3), collapse = " | "), ")") else ""))
    top <- r@result[order(r@result$pvalue), ][seq_len(min(10, nrow(r@result))), ]
    r@result <- top
    pl <- tryCatch(enrichplot::dotplot(r, showCategory = 10, label_format = 45) +
                     ggplot2::labs(title = paste(id, "-", label),
                                   subtitle = if (sig) "terms with the lowest p-values (some p.adjust < 0.05)" else "NO term reached p.adjust < 0.05 - lowest raw p-values shown (exploratory)"),
                   error = function(e) NULL)
    if (!is.null(pl)) ggplot2::ggsave(file.path(OUT, paste0("enrichment_", id, "_", label, ".png")), pl, width = 8.5, height = 6, dpi = 300)
  }
  run("GO_BP", function() clusterProfiler::enrichGO(g, OrgDb = org.Mm.eg.db::org.Mm.eg.db, keyType = "ENTREZID", ont = "BP", universe = uni,
                                                    pvalueCutoff = 1, qvalueCutoff = 1, minGSSize = 5, maxGSSize = 500))
  run("KEGG", function() clusterProfiler::enrichKEGG(g, organism = "mmu", universe = uni, pvalueCutoff = 1, qvalueCutoff = 1, minGSSize = 5, maxGSSize = 500))
}

# ---- run every pattern ------------------------------------------------------------------------------------------------
summ <- list()
for (p in PATTERNS) {
  g <- M[p$sel & !is.na(p$sel), ]
  g <- g[order(eval(p$order, g)), ]
  message(sprintf("\n%s: %d genes", p$id, nrow(g)))
  summ[[p$id]] <- tibble(pattern = p$id, description = p$title, genes = nrow(g))
  if (!nrow(g)) next
  write_csv(g, file.path(OUT, paste0("genes_", p$id, ".csv")))
  message("  top: ", paste(head(g$gene_name, 12), collapse = ", "))
  if (nrow(g) >= 2) draw_heatmap(g, p)
  enrich(g$gene_id, p$id)
}
write_csv(bind_rows(summ), file.path(OUT, "pattern_summary.csv"))
message("\nAll files in: ", OUT)
