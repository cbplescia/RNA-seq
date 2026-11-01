#!/usr/bin/env Rscript
# Step 3 of the intron / exon check: does the change in EBOV vs Mock come with more NASCENT (intronic) signal, or only more mature (exonic) signal?
#   * transcription up      -> intronic reads (pre-mRNA) rise together with exonic reads  (points on the diagonal)
#   * mRNA stabilised / better processed (post-transcriptional) -> exonic reads rise, intronic reads do not (points above the diagonal)
# Counts come from salmon_quant_introns (130_intron_salmon.sh): exons = all transcripts of the gene, introns = the intron features of 129 (your gene + N control genes). Both are tested with DESeq2 (Mock vs EBOV), with the SAME size factors (taken from the exon counts).
# The control genes show what 'exon and intron move together' looks like in your data; the query gene (default Top1) is judged against that cloud. delta = log2FC(exon) - log2FC(intron); z = delta / sqrt(se_exon^2 + se_intron^2) is an approximate screening statistic, not an exact test.
# Outputs (requant_v2_outputs/figures_<gene>_intron_exon/): exon_vs_intron_scatter, query_gene_exon_intron (log2FC with SE + per-replicate counts), delta_distribution (.png/.pdf), exon_intron_all_genes.tsv.   Env: GENE (Top1), MIN_EXON (20: min total exon reads across samples for a control gene), MIN_INTRON (10), PROJ.
# CAUTION: 3 vs 3; intronic reads also come from unspliced transcripts and retained-intron isoforms, and a gene with few intronic reads has a noisy intron log2FC (look at its SE). Treat the result as evidence for or against transcription-driven change, not as proof.
# Usage: Rscript 131_intron_exon_analysis.R   (PROJ from 'source setup_session.sh')

suppressPackageStartupMessages({ library(readr); library(dplyr); library(tidyr); library(ggplot2); library(DESeq2) })
select <- dplyr::select; filter <- dplyr::filter
GENE <- Sys.getenv("GENE", "Top1"); MIN_EXON <- as.numeric(Sys.getenv("MIN_EXON", "20")); MIN_INTRON <- as.numeric(Sys.getenv("MIN_INTRON", "10")); INK <- "#1F2933"
PROJ <- Sys.getenv("PROJ", unset = file.path(
  Sys.getenv("HOME"),
  "Documents/aviti_rnaseq/2025/aviti",
  "20251009-100725_Maury_25161_WivfFacICOTHFAUJwSPeHBEFfCfNtCRpUSgVGUkz",
  "results/Maury_25161"))
RNA <- file.path(PROJ, "requant_v2_outputs"); OUT <- file.path(RNA, paste0("figures_", tolower(GENE), "_intron_exon")); dir.create(OUT, showWarnings = FALSE); QD <- file.path(PROJ, "salmon_quant_introns")
save2 <- function(p, name, w, h) { ggsave(file.path(OUT, paste0(name, ".png")), p, width = w, height = h, dpi = 300); ggsave(file.path(OUT, paste0(name, ".pdf")), p, width = w, height = h) }
theme_t <- theme_bw(base_size = 11) + theme(panel.grid.minor = element_blank(), panel.border = element_rect(colour = "grey70"), axis.text = element_text(colour = INK), plot.title = element_text(colour = INK, face = "bold", size = 11), plot.subtitle = element_text(colour = "grey35", size = 8.5), strip.background = element_rect(fill = "grey95", colour = "grey70"))
nov <- function(x) sub("\\..*$", "", x)

smp <- read_tsv(file.path(PROJ, "intron_samples.txt"), col_names = c("sample", "group"), col_types = "cc"); smp <- smp |> filter(file.exists(file.path(QD, sample, "quant.sf")))
if (n_distinct(smp$group) < 2 || min(table(smp$group)) < 2) stop("need >= 2 quantified samples in each of the two groups (found: ", paste(smp$sample, "(", smp$group, ")", collapse = " ", sep = ""), ") - run 130_intron_salmon.sh")
tx2g <- read_tsv(file.path(PROJ, "tx2gene_GRCm39.109_plus_model.tsv"), show_col_types = FALSE, progress = FALSE) |> transmute(tx = nov(TXNAME), gene = nov(GENEID))
imap <- read_tsv(file.path(PROJ, "ref", "introns_tx2gene.tsv"), show_col_types = FALSE) |> transmute(tx = TXNAME, gene = nov(GENEID))
gl <- read_tsv(file.path(PROJ, "ref", "intron_genes.tsv"), show_col_types = FALSE) |> mutate(gene_id = nov(gene_id)); qid <- gl$gene_id[gl$role == "query"][1]; message("query gene: ", GENE, " = ", qid, "; ", nrow(gl) - 1, " control genes")
count_by <- function(q, map, genes) { m <- map |> filter(gene %in% genes); d <- q |> mutate(tx = if (identical(map, imap)) Name else nov(Name)) |> inner_join(m, by = "tx"); tapply(d$NumReads, d$gene, sum) }
EX <- IN <- NULL
for (i in seq_len(nrow(smp))) { q <- read_tsv(file.path(QD, smp$sample[i], "quant.sf"), col_types = "cdddd", progress = FALSE)
  e <- count_by(q, tx2g, gl$gene_id); n <- count_by(q, imap, gl$gene_id); ev <- setNames(rep(0, nrow(gl)), gl$gene_id); iv <- ev; ev[names(e)] <- e; iv[names(n)] <- n
  EX <- cbind(EX, ev); IN <- cbind(IN, iv) }
colnames(EX) <- colnames(IN) <- smp$sample
cdat <- data.frame(row.names = smp$sample, condition = factor(ifelse(smp$group == "Mock_EBOV", "EBOV", "Mock"), levels = c("Mock", "EBOV")))
keep <- (rowSums(EX) >= MIN_EXON & rowSums(IN) >= MIN_INTRON) | rownames(EX) == qid; message(sum(keep), " genes with enough exonic and intronic reads (", sum(!keep), " dropped)")
if (!qid %in% rownames(EX)[keep]) stop("query gene missing")
if (rowSums(IN)[qid] < MIN_INTRON) message("WARNING: ", GENE, " has only ", rowSums(IN)[qid], " intronic reads in total - its intron log2FC will be very uncertain")
EXk <- round(EX[keep, , drop = FALSE]); INk <- round(IN[keep, , drop = FALSE])
dE <- DESeqDataSetFromMatrix(EXk, cdat, ~ condition); dE <- estimateSizeFactors(dE); sfE <- sizeFactors(dE); dE <- DESeq(dE, quiet = TRUE)
dI <- DESeqDataSetFromMatrix(INk, cdat, ~ condition); sizeFactors(dI) <- sfE; dI <- estimateDispersions(dI, quiet = TRUE); dI <- nbinomWaldTest(dI, quiet = TRUE)
rE <- as.data.frame(results(dE, contrast = c("condition", "EBOV", "Mock"))); rI <- as.data.frame(results(dI, contrast = c("condition", "EBOV", "Mock")))
T <- tibble(gene_id = rownames(rE), lfc_exon = rE$log2FoldChange, se_exon = rE$lfcSE, padj_exon = rE$padj, lfc_intron = rI$log2FoldChange, se_intron = rI$lfcSE, padj_intron = rI$padj, exon_reads = rowSums(EXk), intron_reads = rowSums(INk)) |>
  left_join(gl |> select(gene_id, gene_name, role), by = "gene_id") |> mutate(delta = lfc_exon - lfc_intron, z = delta / sqrt(se_exon^2 + se_intron^2), p_delta = 2 * pnorm(-abs(z)), is_q = gene_id == qid)
# GLOBAL shift: if most genes lose intronic signal relative to exonic signal (less nascent RNA in the infected samples, or a sample-prep difference), every gene's delta is shifted; judge the query gene against the controls' median, not against 0
med_delta <- median(T$delta[!T$is_q & is.finite(T$delta)]); T <- T |> mutate(delta_adj = delta - med_delta, z_adj = delta_adj / sqrt(se_exon^2 + se_intron^2), p_delta_adj = 2 * pnorm(-abs(z_adj)))
write_tsv(T |> arrange(desc(abs(z_adj))), file.path(OUT, "exon_intron_all_genes.tsv"))
ctl <- T |> filter(!is_q, is.finite(lfc_exon), is.finite(lfc_intron)); q <- T |> filter(is_q)
fit <- lm(lfc_exon ~ lfc_intron, data = ctl); r_c <- suppressWarnings(cor(ctl$lfc_exon, ctl$lfc_intron, method = "spearman")); pct <- 100 * mean(ctl$delta <= q$delta, na.rm = TRUE)
message(sprintf("\n%s: exon log2FC %+.2f (SE %.2f), intron log2FC %+.2f (SE %.2f); delta (exon - intron) %+.2f, z %.2f; delta is higher than %.0f%% of the %d control genes", GENE, q$lfc_exon, q$se_exon, q$lfc_intron, q$se_intron, q$delta, q$z, pct, nrow(ctl)))
message(sprintf("controls: median exon log2FC %+.2f, median intron log2FC %+.2f, so the GLOBAL exon-minus-intron shift is %+.2f; Spearman correlation of exon vs intron log2FC = %.2f", median(ctl$lfc_exon, na.rm = TRUE), median(ctl$lfc_intron, na.rm = TRUE), med_delta, r_c))
if (abs(med_delta) > 0.3) message("NOTE: the intron/exon balance shifts across ALL genes (", sprintf("%+.2f", med_delta), " log2). That is a global feature of these samples (less nascent signal relative to mature RNA in one group, a real drop in host transcription, or a library-prep difference) and not specific to ", GENE, "; the gene-level reading below is therefore made relative to the control median.")
message(sprintf("%s relative to the control median: delta %+.2f, z %.2f", GENE, q$delta_adj, q$z_adj))
verdict <- if (!is.finite(q$z_adj)) "no verdict (missing values)" else if (q$z_adj > 2) "exon rises MORE than intron -> points to a post-transcriptional contribution (mRNA stability / processing), not just more transcription" else if (q$z_adj < -2) "intron rises MORE than exon -> nascent signal up without matching mature RNA (slower processing / export, or intron retention)" else "after allowing for the global shift, exon and intron move like a typical control gene -> no gene-specific evidence of extra transcription or of a separate post-transcriptional effect (an intron signal that did not rise argues against a transcriptional increase)"
message("reading: ", verdict)
lim <- c(-2, 2); sc <- ggplot(ctl, aes(lfc_intron, lfc_exon)) + geom_hline(yintercept = 0, colour = "grey80") + geom_vline(xintercept = 0, colour = "grey80") + geom_abline(slope = 1, intercept = 0, colour = "#9AA5B1", linetype = "dashed") + geom_abline(slope = 1, intercept = med_delta, colour = "#2B6CB0", linetype = "dotted") + geom_point(colour = "#C9CDD3", size = 0.9, alpha = 0.7) +
  geom_point(data = ctl |> filter(is.finite(padj_exon), padj_exon < 0.05), colour = "#1B7837", size = 1.1) + geom_errorbar(data = q, aes(ymin = lfc_exon - se_exon, ymax = lfc_exon + se_exon), width = 0, colour = "#C53030") + geom_errorbarh(data = q, aes(xmin = lfc_intron - se_intron, xmax = lfc_intron + se_intron), height = 0, colour = "#C53030") +
  geom_point(data = q, colour = "#C53030", size = 3) + annotate("text", x = q$lfc_intron, y = q$lfc_exon, label = paste0("  ", GENE), colour = "#C53030", hjust = 0, vjust = -0.8, fontface = "italic") + coord_cartesian(xlim = lim, ylim = lim) +
  labs(title = paste0("Exonic vs intronic change, EBOV vs Mock: ", GENE, " among ", nrow(ctl), " control genes"), subtitle = paste0("dashed = exon and intron change equally (transcription-like); dotted blue = the typical control gene (global shift); above = exon rises more (post-transcriptional); green = control genes significant in exons (padj < 0.05); red bars = +/- 1 SE; control Spearman rho = ", sprintf("%.2f", r_c)), x = "intron log2FC (EBOV / Mock)", y = "exon log2FC (EBOV / Mock)") + theme_t
save2(sc, "exon_vs_intron_scatter", 7.2, 6.6)
nb <- tibble(type = factor(c("exon (mature)", "intron (nascent)"), levels = c("exon (mature)", "intron (nascent)")), lfc = c(q$lfc_exon, q$lfc_intron), se = c(q$se_exon, q$se_intron))
pa <- ggplot(nb, aes(type, lfc, fill = type)) + geom_col(width = 0.55) + geom_errorbar(aes(ymin = lfc - se, ymax = lfc + se), width = 0.2, colour = INK) + geom_hline(yintercept = 0, colour = INK, linewidth = 0.4) + scale_fill_manual(values = c("#1B7837", "#9BD08E"), guide = "none") + labs(title = paste0(GENE, ": exon vs intron log2FC"), x = NULL, y = "log2FC EBOV / Mock (+/- SE)") + theme_t
cn <- bind_rows(tibble(sample = colnames(EX), type = "exon (mature)", n = as.numeric(EX[qid, ]) / sfE[colnames(EX)]), tibble(sample = colnames(IN), type = "intron (nascent)", n = as.numeric(IN[qid, ]) / sfE[colnames(IN)])) |> left_join(smp, by = "sample") |> mutate(type = factor(type, levels = levels(nb$type)), group = factor(ifelse(group == "Mock_EBOV", "EBOV", "Mock"), levels = c("Mock", "EBOV")))
pb <- ggplot(cn, aes(group, n, fill = group)) + geom_col(data = cn |> group_by(type, group) |> summarise(n = mean(n), .groups = "drop"), alpha = 0.5, width = 0.6) + geom_point(position = position_jitter(width = 0.08, height = 0, seed = 1), shape = 21, size = 2.4, colour = INK, fill = "white") + facet_wrap(~type, scales = "free_y") + scale_fill_manual(values = c(Mock = "#7F7F7F", EBOV = "#1B7837"), guide = "none") +
  labs(title = paste0(GENE, ": normalised reads per replicate"), subtitle = "same size factors for exons and introns", x = NULL, y = "normalised reads") + theme_t
save2(pa, "query_gene_exon_intron_lfc", 4.2, 4.2); save2(pb, "query_gene_exon_intron_counts", 6.2, 4.2)
pd <- ggplot(ctl, aes(delta)) + geom_histogram(bins = 60, fill = "#C9CDD3", colour = "white") + geom_vline(xintercept = q$delta, colour = "#C53030", linewidth = 0.8) + annotate("text", x = q$delta, y = Inf, label = paste0(" ", GENE, "  (", sprintf("%.0f", pct), "th percentile)"), colour = "#C53030", hjust = 0, vjust = 1.6, fontface = "italic") +
  labs(title = "exon log2FC - intron log2FC across the control genes", subtitle = "0 = exon and intron change equally; positive = exon rises more", x = "delta (exon - intron log2FC)", y = "genes") + theme_t
save2(pd, "delta_distribution", 7, 4)
message("\nSaved to: ", OUT)
