#!/usr/bin/env Rscript
# First-pass checks on HOW a gene's transcripts go up (default Top1), using data you already have:
#   1  isoform / transcript level   which Salmon transcripts carry the change; does isoform usage (share of the gene's TPM) shift?      -> tx_composition_TPM, tx_isoform_fraction, tx_normalised_counts (.png/.pdf) + top1_transcripts.tsv
#   2  what does the gene track?    Top1 VST vs EBOV read fraction, an interferon-stimulated-gene score, a cell-cycle score and Top2a across all samples (Spearman)   -> top1_correlations (.png/.pdf, .tsv)
#   3  TF ACTIVITY (not motifs)     decoupleR (ULM) with the CollecTRI mouse regulons on each DE contrast: which TFs look active, and what the curated regulators of the gene do     -> tf_activity_* (.png/.pdf, .tsv)
#   4  artefact checks              (a) share of the gene's reads that are multi-mapping (Salmon unique vs ambiguous counts) and where that sits among all transcripts; (b) log2FC under different normalisations      -> artefact_* (.png/.pdf, .tsv)
# Paralog / retrocopy check by k-mer sharing: python 127_kmer_paralog_check.py   (companion script).
# Reads salmon_quant_v2/<sample>/quant.sf + aux_info/ambig_info.tsv, dds_v2_decoyAware_hostNormalised.rds, vsd_v2.rds, DE_v2_*_annotated.csv, construct_features_per_sample.csv, the GTF.
# Section 3 needs the Bioconductor package decoupleR (and OmnipathR + internet for the CollecTRI download):  conda install -c conda-forge -c bioconda bioconductor-decoupler bioconductor-omnipathr   (or BiocManager::install(c("decoupleR","OmnipathR"))). Without it that section is skipped.
# Env: GENE (default Top1), GROUPS (groups shown in the transcript plots; default Mock_Mock,Mock_EBOV,IFNB_Mock,IFNB_EBOV,IFNY_Mock,IFNY_EBOV), SKIP_TF=1 to skip section 3.   Writes to requant_v2_outputs/figures_<gene>_mechanism/.
# CAUTION: 3 replicates per group; the tests here rank and describe, they do not prove a mechanism. Transcript-level quantification of isoforms that share most of their sequence is the least certain part (see the ambiguity panel).
# Usage: Rscript 126_top1_mechanism_checks.R   (PROJ from 'source setup_session.sh')

suppressPackageStartupMessages({ library(readr); library(dplyr); library(tidyr); library(ggplot2); library(DESeq2) })
select <- dplyr::select; filter <- dplyr::filter

GENE <- Sys.getenv("GENE", "Top1"); SHOW <- trimws(strsplit(Sys.getenv("GROUPS", "Mock_Mock,Mock_EBOV,IFNB_Mock,IFNB_EBOV,IFNY_Mock,IFNY_EBOV"), ",")[[1]]); INK <- "#1F2933"
ORD <- c("Mock_Mock", "IFNB_Mock", "IFNY_Mock", "Mock_EBOV", "IFNB_EBOV", "IFNY_EBOV", "Mock_MARV", "IFNB_MARV", "IFNY_MARV", "Mock_EVLP", "Mock_MVLP")
LAB <- c(Mock_Mock = "Mock", IFNB_Mock = "IFNB", IFNY_Mock = "IFNG", Mock_EBOV = "EBOV", IFNB_EBOV = "IFNB+EBOV", IFNY_EBOV = "IFNG+EBOV", Mock_MARV = "MARV", IFNB_MARV = "IFNB+MARV", IFNY_MARV = "IFNG+MARV", Mock_EVLP = "EVLP", Mock_MVLP = "MVLP")
COL <- c(Mock = "#7F7F7F", EBOV = "#1B7837", IFNB = "#2B6CB0", `IFNB+EBOV` = "#9BD08E", IFNG = "#D69E2E", `IFNG+EBOV` = "#E8C77A", MARV = "#8C510A", `IFNB+MARV` = "#C7A27C", `IFNG+MARV` = "#DCC9A8", EVLP = "#6A51A3", MVLP = "#CBA4D8")
PROJ <- Sys.getenv("PROJ", unset = file.path(
  Sys.getenv("HOME"),
  "Documents/aviti_rnaseq/2025/aviti",
  "20251009-100725_Maury_25161_WivfFacICOTHFAUJwSPeHBEFfCfNtCRpUSgVGUkz",
  "results/Maury_25161"))
RNA <- file.path(PROJ, "requant_v2_outputs"); QUANT <- file.path(PROJ, "salmon_quant_v2"); OUT <- file.path(RNA, paste0("figures_", tolower(GENE), "_mechanism")); dir.create(OUT, showWarnings = FALSE)
save2 <- function(p, name, w, h) { ggsave(file.path(OUT, paste0(name, ".png")), p, width = w, height = h, dpi = 300); ggsave(file.path(OUT, paste0(name, ".pdf")), p, width = w, height = h) }
theme_t <- theme_bw(base_size = 11) + theme(panel.grid.minor = element_blank(), panel.border = element_rect(colour = "grey70"), axis.text = element_text(colour = INK), plot.title = element_text(colour = INK, face = "bold", size = 11), plot.subtitle = element_text(colour = "grey35", size = 8.5), strip.background = element_rect(fill = "grey95", colour = "grey70"))
glab <- function(g) ifelse(g %in% names(LAB), LAB[g], g)
stars <- function(p) ifelse(is.na(p), "", ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", "ns"))))

dds <- readRDS(file.path(RNA, "dds_v2_decoyAware_hostNormalised.rds")); vsd <- readRDS(file.path(RNA, "vsd_v2.rds"))
cd <- as.data.frame(colData(dds)); cd$sample_id <- rownames(cd); cd$group <- as.character(cd$group); smp <- cd$sample_id
de0 <- read_csv(file.path(RNA, "DE_v2_EBOV_vs_Mock_annotated.csv"), show_col_types = FALSE)
gid <- de0$gene_id[tolower(de0$gene_name) == tolower(GENE)][1]; if (is.na(gid)) stop(GENE, " is not in the DE table (removed by the pre-filter in 04?)")
message(GENE, " = ", gid)
cf <- file.path(RNA, "construct_features_per_sample.csv"); viral_ids <- if (file.exists(cf)) unique(read_csv(cf, show_col_types = FALSE)$feature) else character(0)
ncts <- counts(dds); vir <- intersect(viral_ids, rownames(ncts)); viral_reads <- if (length(vir)) colSums(ncts[vir, , drop = FALSE]) else rep(0, ncol(ncts)); host_reads <- colSums(ncts[setdiff(rownames(ncts), vir), , drop = FALSE])
viral_frac <- viral_reads / (viral_reads + host_reads)
# per-sample size factor. A tximport-based DESeq2 object stores normalizationFactors (a gene x sample matrix incl. transcript-length offsets) and has NO sizeFactors, so derive one per sample from it
sf <- sizeFactors(dds)
if (is.null(sf)) { nf <- normalizationFactors(dds); if (is.null(nf)) stop("the DESeq2 object has neither sizeFactors nor normalizationFactors"); sf <- exp(colMeans(log(nf[setdiff(rownames(nf), vir), , drop = FALSE]))); sf <- sf / exp(mean(log(sf))); message("size factors taken from the geometric mean of the normalizationFactors (tximport object)") }
names(sf) <- colnames(dds)

# ================= 1  transcript level =================
tx2gene <- read_tsv(file.path(PROJ, "tx2gene_GRCm39.109_plus_model.tsv"), show_col_types = FALSE, progress = FALSE); nov <- function(x) sub("\\..*$", "", x)
tx_ids <- tx2gene$TXNAME[nov(tx2gene$GENEID) == nov(gid)]; if (!length(tx_ids)) stop("no transcripts for ", GENE, " in the tx2gene file"); message(length(tx_ids), " transcripts in the index")
gtf <- file.path(PROJ, "ref", "Mus_musculus.GRCm39.109.gtf")
ann <- tryCatch({
  ln <- readLines(pipe(sprintf("grep -F 'gene_name \"%s\"' '%s'", GENE, gtf))); f <- strsplit(ln, "\t"); f <- f[vapply(f, length, 1L) >= 9]; f <- f[vapply(f, function(x) x[3] == "transcript", TRUE)]
  gv <- function(x, k) if (grepl(paste0(k, ' "'), x)) sub(paste0('.*', k, ' "([^"]+)".*'), "\\1", x) else NA_character_
  tibble(tid = vapply(f, function(x) gv(x[9], "transcript_id"), ""), tname = vapply(f, function(x) gv(x[9], "transcript_name"), ""), biotype = vapply(f, function(x) gv(x[9], "transcript_biotype"), ""), canonical = vapply(f, function(x) grepl("Ensembl_canonical", x[9]), TRUE)) }, error = function(e) NULL)
rd <- function(s) {
  q <- read_tsv(file.path(QUANT, s, "quant.sf"), col_types = "cdddd", progress = FALSE); ai <- file.path(QUANT, s, "aux_info", "ambig_info.tsv")
  a <- if (file.exists(ai)) read_tsv(ai, col_types = "dd", progress = FALSE) else tibble(UniqueCount = NA_real_, AmbigCount = NA_real_)[rep(1, nrow(q)), ]
  if (nrow(a) != nrow(q)) a <- tibble(UniqueCount = rep(NA_real_, nrow(q)), AmbigCount = NA_real_); q$unique <- a$UniqueCount; q$ambig <- a$AmbigCount; q }
message("reading ", length(smp), " quant.sf files ...")
Q <- lapply(smp, rd); names(Q) <- smp
TX <- bind_rows(lapply(smp, function(s) { q <- Q[[s]]; q <- q[nov(q$Name) %in% nov(tx_ids), ]; tibble(sample = s, tid = q$Name, TPM = q$TPM, reads = q$NumReads, unique = q$unique, ambig = q$ambig) }))
TX <- TX |> left_join(cd |> select(sample = sample_id, group), by = "sample") |> mutate(norm = reads / sf[sample], tidn = nov(tid))
lab_tx <- function(t) { a <- if (!is.null(ann)) ann[match(nov(t), nov(ann$tid)), ] else NULL
  if (is.null(a)) return(t); paste0(ifelse(is.na(a$tname) | a$tname == "", t, a$tname), " (", ifelse(is.na(a$biotype), "", gsub("_", " ", a$biotype)), ifelse(!is.na(a$canonical) & a$canonical, ", canonical", ""), ")") }
TX$label <- lab_tx(TX$tid)
tot <- TX |> group_by(sample) |> summarise(gene_TPM = sum(TPM), .groups = "drop"); TX <- TX |> left_join(tot, by = "sample") |> mutate(frac = ifelse(gene_TPM > 0, TPM / gene_TPM, NA_real_))
keep_tx <- TX |> group_by(label) |> summarise(m = max(tapply(TPM, group, mean)), .groups = "drop") |> filter(m >= 0.2) |> pull(label); TXs <- TX |> filter(group %in% SHOW) |> mutate(label = ifelse(label %in% keep_tx, label, "other (low)"), group = factor(glab(group), levels = glab(intersect(ORD, SHOW))))
gm <- TXs |> group_by(group, label) |> summarise(TPM = mean(TPM), frac = mean(frac, na.rm = TRUE), .groups = "drop")
p1 <- ggplot(gm, aes(group, TPM, fill = label)) + geom_col(width = 0.7, colour = "white", linewidth = 0.2) + labs(title = paste0(GENE, ": transcript composition (mean TPM per group)"), subtitle = "stacked Salmon transcripts of the gene; each colour is one transcript", x = NULL, y = paste0(GENE, " TPM"), fill = NULL) + theme_t + theme(axis.text.x = element_text(angle = 40, hjust = 1), legend.position = "right")
save2(p1, "tx_composition_TPM", 8, 4.8)
p2 <- ggplot(gm, aes(group, frac, fill = label)) + geom_col(width = 0.7, colour = "white", linewidth = 0.2) + labs(title = paste0(GENE, ": isoform usage (share of the gene's TPM)"), subtitle = "a shift in colours between groups = isoform switching (promoter choice, splicing or 3' end processing)", x = NULL, y = "share of gene TPM", fill = NULL) + theme_t + theme(axis.text.x = element_text(angle = 40, hjust = 1), legend.position = "right")
save2(p2, "tx_isoform_fraction", 8, 4.8)
pn <- TXs |> filter(label != "other (low)") |> mutate(l2 = log2(norm + 1))
if (nrow(pn)) { p3 <- ggplot(pn, aes(group, l2, fill = group)) + geom_boxplot(outlier.shape = NA, alpha = 0.75, width = 0.6, colour = INK, linewidth = 0.3) + geom_point(position = position_jitter(width = 0.08, height = 0, seed = 1), size = 1.6, shape = 21, fill = "white", colour = INK, stroke = 0.4) +
    facet_wrap(~label, scales = "free_y") + scale_fill_manual(values = COL, guide = "none") + labs(title = paste0(GENE, ": each transcript, host-normalised counts"), subtitle = "log2(normalised reads + 1); every replicate shown", x = NULL, y = "log2(normalised counts + 1)") + theme_t + theme(axis.text.x = element_text(angle = 40, hjust = 1))
  save2(p3, "tx_normalised_counts", 3 + 3 * min(3, n_distinct(pn$label)), 2.5 + 2.6 * ceiling(n_distinct(pn$label) / 3)) }
tt <- TX |> filter(group %in% c("Mock_Mock", "Mock_EBOV")) |> group_by(tid, label) |> group_modify(~ { a <- .x$norm[.x$group == "Mock_EBOV"]; b <- .x$norm[.x$group == "Mock_Mock"]; fa <- .x$frac[.x$group == "Mock_EBOV"]; fb <- .x$frac[.x$group == "Mock_Mock"]
  tibble(mean_mock = mean(b), mean_ebov = mean(a), log2FC_EBOV_vs_Mock = log2((mean(a) + 1) / (mean(b) + 1)), p_counts = tryCatch(t.test(log2(a + 1), log2(b + 1))$p.value, error = function(e) NA_real_),
         frac_mock = mean(fb, na.rm = TRUE), frac_ebov = mean(fa, na.rm = TRUE), p_isoform_share = tryCatch(t.test(fa, fb)$p.value, error = function(e) NA_real_)) }) |> ungroup() |> mutate(padj_counts = p.adjust(p_counts, "BH"), padj_isoform_share = p.adjust(p_isoform_share, "BH"))
write_tsv(tt |> arrange(desc(mean_ebov)), file.path(OUT, paste0(tolower(GENE), "_transcripts.tsv")))
message("\n--- ", GENE, " transcripts, EBOV vs Mock (host-normalised reads; Welch t on log2(x+1), 3 vs 3) ---"); print(as.data.frame(tt |> mutate(across(where(is.numeric), ~ signif(.x, 3))) |> arrange(desc(mean_ebov))), row.names = FALSE)

# ================= 2  what does the gene track across samples =================
vm <- SummarizedExperiment::assay(vsd); nm2id <- setNames(de0$gene_id, tolower(de0$gene_name)); nm2id <- nm2id[!duplicated(names(nm2id))]
score <- function(genes) { ids <- nm2id[tolower(genes)]; ids <- ids[!is.na(ids) & ids %in% rownames(vm)]; if (length(ids) < 3) return(rep(NA_real_, ncol(vm))); z <- t(scale(t(vm[ids, , drop = FALSE]))); colMeans(z) }
ISG <- c("Isg15", "Ifit1", "Ifit2", "Ifit3", "Mx1", "Mx2", "Oas1a", "Oasl2", "Irf7", "Usp18", "Rsad2", "Stat1", "Stat2", "Irf9", "Bst2", "Ifitm3", "Cxcl10", "Ddx58", "Ifih1", "Gbp2")
CC <- c("Mki67", "Pcna", "Mcm2", "Mcm3", "Mcm4", "Mcm5", "Mcm6", "Mcm7", "Ccnb1", "Ccna2", "Cdk1", "Cdc20", "Aurkb", "Plk1", "Ccne1", "Rrm2", "Tyms")
P <- tibble(sample = colnames(vm), group = cd$group[match(colnames(vm), cd$sample_id)], gene = as.numeric(vm[gid, ]), viral = log10(viral_frac[colnames(vm)] + 1e-6), ISG_score = score(ISG), cellcycle_score = score(CC),
            Top2a = if ("top2a" %in% names(nm2id) && nm2id["top2a"] %in% rownames(vm)) as.numeric(vm[nm2id["top2a"], ]) else NA_real_)
Pl <- P |> pivot_longer(c(viral, ISG_score, cellcycle_score, Top2a), names_to = "predictor", values_to = "x") |> filter(!is.na(x)) |>
  mutate(predictor = factor(predictor, levels = c("viral", "ISG_score", "cellcycle_score", "Top2a"), labels = c("log10 construct (EBOV) read fraction", "ISG score (mean z of 20 ISGs)", "cell-cycle score (mean z of 17 genes)", "Top2a VST")), glabel = factor(glab(group), levels = glab(ORD)))
isE <- P$group %in% c("Mock_EBOV", "IFNB_EBOV", "IFNY_EBOV"); rc <- function(x, y) { ok <- is.finite(x) & is.finite(y); if (sum(ok) < 4) return(c(NA, NA)); ct <- suppressWarnings(cor.test(x[ok], y[ok], method = "spearman", exact = FALSE)); c(unname(ct$estimate), ct$p.value) }
cors <- bind_rows(lapply(c("viral", "ISG_score", "cellcycle_score", "Top2a"), function(v) { a <- rc(P$gene, P[[v]]); b <- rc(P$gene[isE], P[[v]][isE]); tibble(predictor = v, rho_all = a[1], p_all = a[2], n_all = sum(is.finite(P[[v]])), rho_EBOV_samples = b[1], p_EBOV_samples = b[2], n_EBOV = sum(is.finite(P[[v]][isE]))) }))
write_tsv(cors, file.path(OUT, paste0(tolower(GENE), "_correlations.tsv"))); message("\n--- Spearman correlation of ", GENE, " VST with ---"); print(as.data.frame(cors |> mutate(across(where(is.numeric), ~ signif(.x, 3)))), row.names = FALSE)
lab_df <- cors |> mutate(predictor = factor(predictor, levels = c("viral", "ISG_score", "cellcycle_score", "Top2a"), labels = levels(Pl$predictor)), txt = sprintf("all: rho %.2f (p %.2g)\nEBOV samples: rho %.2f (p %.2g)", rho_all, p_all, rho_EBOV_samples, p_EBOV_samples))
p4 <- ggplot(Pl, aes(x, gene)) + geom_smooth(aes(group = 1), method = "lm", se = FALSE, colour = "grey50", linewidth = 0.4, formula = y ~ x) + geom_point(aes(fill = glabel), shape = 21, size = 2.4, colour = INK, stroke = 0.3) +
  geom_text(data = lab_df, aes(x = -Inf, y = Inf, label = txt), inherit.aes = FALSE, hjust = -0.05, vjust = 1.15, size = 2.6, colour = INK) + facet_wrap(~predictor, scales = "free_x", nrow = 1) + scale_fill_manual(values = COL, name = NULL) +
  labs(title = paste0(GENE, " expression across samples vs what might drive it"), subtitle = "each point = one sample; Spearman correlation; EBOV read fraction is only informative among the infected samples", x = NULL, y = paste0(GENE, " VST")) + theme_t + theme(legend.position = "bottom")
save2(p4, paste0(tolower(GENE), "_correlations"), 13, 4.6)

# ================= 3  TF activity (decoupleR + CollecTRI) =================
if (Sys.getenv("SKIP_TF") != "1") tryCatch({
  if (!requireNamespace("decoupleR", quietly = TRUE)) stop("decoupleR not installed")
  CON <- c(IFNB_vs_Mock = "IFNB", IFNY_vs_Mock = "IFNG", EBOV_vs_Mock = "EBOV", IFNB_EBOV_vs_Mock = "IFNB+EBOV", IFNY_EBOV_vs_Mock = "IFNG+EBOV")
  mats <- lapply(names(CON), function(k) { f <- file.path(RNA, paste0("DE_v2_", k, "_annotated.csv")); if (!file.exists(f)) return(NULL); d <- read_csv(f, show_col_types = FALSE) |> filter(!is.na(gene_name), gene_name != "")
    st <- if ("stat" %in% names(d) && any(is.finite(d$stat))) d$stat else sign(d$log2FoldChange) * -log10(pmax(if ("pvalue" %in% names(d)) d$pvalue else d$padj, 1e-300)); tibble(g = d$gene_name, v = st) |> filter(is.finite(v)) |> group_by(g) |> slice_max(abs(v), n = 1, with_ties = FALSE) |> ungroup() |> mutate(con = CON[[k]]) })
  mm <- bind_rows(mats); mat <- mm |> pivot_wider(names_from = con, values_from = v) |> as.data.frame(); rownames(mat) <- mat$g; mat$g <- NULL; mat <- as.matrix(mat); mat[is.na(mat)] <- 0
  message("\ngetting the CollecTRI mouse network (needs internet the first time) ..."); net <- decoupleR::get_collectri(organism = "mouse", split_complexes = FALSE)
  res <- decoupleR::run_ulm(mat, net, .source = "source", .target = "target", .mor = "mor", minsize = 5); res <- as_tibble(res); write_tsv(res, file.path(OUT, "tf_activity_all.tsv"))
  reg <- net |> filter(tolower(target) == tolower(GENE)) |> distinct(source, mor); message("CollecTRI regulators of ", GENE, ": ", if (nrow(reg)) paste(reg$source, collapse = ", ") else "none listed")
  if (nrow(reg)) { rr <- res |> filter(source %in% reg$source) |> left_join(reg, by = "source") |> mutate(condition = factor(condition, levels = unname(CON)), sig = stars(p_value))
    write_tsv(rr, file.path(OUT, paste0("tf_activity_regulators_of_", tolower(GENE), ".tsv")))
    p5 <- ggplot(rr, aes(condition, source, fill = score)) + geom_tile(colour = "white") + geom_text(aes(label = sig), size = 3, colour = INK) + scale_fill_gradient2(low = "#2B6CB0", mid = "#F2F2F2", high = "#C53030", midpoint = 0, name = "TF activity\n(ULM score)") +
      labs(title = paste0("Activity of the curated regulators of ", GENE, " (CollecTRI), inferred from each contrast"), subtitle = "red = target genes of the TF go up, blue = down; stars = ULM p; mode of regulation (activator/repressor) is in the table (mor)", x = NULL, y = NULL) + theme_t
    save2(p5, paste0("tf_activity_regulators_of_", tolower(GENE)), 7.5, 1.8 + 0.28 * n_distinct(rr$source)) }
  top <- res |> group_by(condition) |> slice_max(abs(score), n = 12, with_ties = FALSE) |> ungroup() |> mutate(condition = factor(condition, levels = unname(CON)))
  tfs <- unique(top$source); hm <- res |> filter(source %in% tfs) |> mutate(condition = factor(condition, levels = unname(CON)), sig = stars(p_value)); ord <- hm |> filter(condition == "EBOV") |> arrange(score) |> pull(source); hm$source <- factor(hm$source, levels = unique(c(ord, tfs)))
  p6 <- ggplot(hm, aes(condition, source, fill = score)) + geom_tile(colour = "white") + geom_text(aes(label = sig), size = 2.6, colour = INK) + scale_fill_gradient2(low = "#2B6CB0", mid = "#F2F2F2", high = "#C53030", midpoint = 0, name = "TF activity\n(ULM score)") +
    labs(title = "Inferred TF activity: the 12 strongest TFs of each contrast vs Mock", subtitle = "CollecTRI mouse regulons, univariate linear model on the DESeq2 statistic; stars = ULM p", x = NULL, y = NULL) + theme_t
  save2(p6, "tf_activity_top", 7.5, 2 + 0.2 * n_distinct(hm$source))
}, error = function(e) message("\nSection 3 (TF activity) skipped: ", conditionMessage(e), "\n  install:  conda install -c conda-forge -c bioconda bioconductor-decoupler bioconductor-omnipathr   (and check the internet connection)"))

# ================= 4  artefact checks =================
amb <- bind_rows(lapply(smp, function(s) { q <- Q[[s]]; tot <- q$unique + q$ambig; ok <- is.finite(tot) & tot >= 50; fr <- ifelse(ok, q$ambig / tot, NA_real_); isg <- nov(q$Name) %in% nov(tx_ids)
  tu <- sum(q$unique[isg], na.rm = TRUE); ta <- sum(q$ambig[isg], na.rm = TRUE); f <- if (tu + ta > 0) ta / (tu + ta) else NA_real_
  tibble(sample = s, gene_ambig_fraction = f, gene_reads_unique = tu, gene_reads_ambiguous = ta, percentile_among_transcripts = if (is.finite(f) && any(ok)) 100 * mean(fr[ok] <= f) else NA_real_) })) |> left_join(cd |> select(sample = sample_id, group), by = "sample")
write_tsv(amb, file.path(OUT, paste0(tolower(GENE), "_ambiguity.tsv")))
if (any(is.finite(amb$gene_ambig_fraction))) { pa <- amb |> filter(is.finite(gene_ambig_fraction)) |> mutate(g = factor(glab(group), levels = glab(ORD)))
  p7 <- ggplot(pa, aes(g, 100 * gene_ambig_fraction, fill = g)) + geom_boxplot(outlier.shape = NA, alpha = 0.75, width = 0.6, colour = INK, linewidth = 0.3) + geom_point(position = position_jitter(width = 0.08, height = 0, seed = 1), size = 1.6, shape = 21, fill = "white", colour = INK, stroke = 0.4) + scale_fill_manual(values = COL, guide = "none") +
    labs(title = paste0(GENE, ": share of its reads that are multi-mapping (Salmon ambiguous / (unique + ambiguous))"), subtitle = "includes reads shared between isoforms of the same gene (harmless) as well as with other genes (paralogs / retrocopies); a high value at gene level that differs between groups would be a warning; see 127_kmer_paralog_check.py", x = NULL, y = "% ambiguous reads") + theme_t + theme(axis.text.x = element_text(angle = 40, hjust = 1))
  save2(p7, paste0("artefact_", tolower(GENE), "_ambiguity"), 8, 4.2); message("\n", GENE, " ambiguous-read share, mean per group (%): "); print(as.data.frame(amb |> group_by(group) |> summarise(pct_ambiguous = round(100 * mean(gene_ambig_fraction, na.rm = TRUE), 1), .groups = "drop")), row.names = FALSE) }
raw <- ncts[gid, ]; hostN <- sf; libN <- host_reads / mean(host_reads); hk <- intersect(nm2id[tolower(c("Gapdh", "Actb", "B2m", "Rplp0", "Rpl13a", "Ppia", "Hprt", "Tbp"))], rownames(ncts))
hkN <- if (length(hk) >= 3) { gm_ <- exp(colMeans(log(ncts[hk, , drop = FALSE] + 1))); gm_ / mean(gm_) } else rep(NA_real_, ncol(ncts))
meth <- list(`DESeq2 size factors (host genes; the analysis default)` = hostN, `library size (host reads)` = libN, `housekeeping genes (Gapdh, Actb, B2m, ... may themselves change)` = hkN)
fc <- bind_rows(lapply(names(meth), function(m) { nrm <- raw / meth[[m]]; bind_rows(lapply(setdiff(unique(cd$group), "Mock_Mock"), function(g) tibble(method = m, group = g, log2FC = log2((mean(nrm[cd$group == g]) + 0.5) / (mean(nrm[cd$group == "Mock_Mock"]) + 0.5))))) })) |> filter(group %in% c("Mock_EBOV", "IFNB_Mock", "IFNB_EBOV", "IFNY_Mock", "IFNY_EBOV", "Mock_MARV"))
write_tsv(fc, file.path(OUT, paste0(tolower(GENE), "_log2FC_by_normalisation.tsv")))
p8 <- ggplot(fc |> mutate(g = factor(glab(group), levels = glab(ORD))), aes(g, log2FC, fill = method)) + geom_col(position = position_dodge(width = 0.8), width = 0.75) + geom_hline(yintercept = 0, colour = INK, linewidth = 0.4) + scale_fill_manual(values = c("#1F2933", "#7F7F7F", "#C9CDD3"), name = NULL) +
  labs(title = paste0(GENE, ": is the change robust to how the samples are normalised?"), subtitle = "log2 of group mean / Mock mean of normalised reads (+0.5); similar bars across methods = not a normalisation artefact", x = NULL, y = "log2 fold-change vs Mock") + theme_t + theme(legend.position = "bottom", axis.text.x = element_text(angle = 35, hjust = 1)) + guides(fill = guide_legend(ncol = 1))
save2(p8, paste0("artefact_", tolower(GENE), "_normalisation"), 8.5, 5)
message("\nSaved to: ", OUT)
