#!/usr/bin/env Rscript
# Is Top1's small rise in EBOV part of an NRF2 (oxidative-stress / detox) response?   EBOV vs Mock; NRF2 = Nfe2l2.
#   A  bar graph (EBOV vs Mock shrunk log2FC +/- 1 SE, padj stars): Top1 | canonical NRF2 targets | NRF2 pathway regulators (Nfe2l2, Keap1, Bach1, Maff)                                           -> bars_nrf2_top1_EBOV
#   B  is the NRF2 program up as a whole?  Wald statistics of the canonical NRF2 targets vs all other genes (Wilcoxon), shown on the ranked gene list ("barcode") with Top1 marked          -> nrf2_barcode_EBOV
#   C  NRF2 score per sample (mean z-score of the target genes' VST) in every group, and Top1 VST against that score across samples (Spearman)                                              -> nrf2_score_by_group, top1_vs_nrf2_score
#   D  TF activity with decoupleR (ULM) and the CollecTRI mouse regulons on the EBOV vs Mock Wald statistic: the strongest TFs, NRF2's activity, whether Top1 is a curated NRF2 target, and the curated regulators of Top1 with their activity   -> tf_activity_EBOV_top, tf_activity_regulators_of_top1, nrf2_regulon_targets, + .tsv
# The NRF2 target list below is curated from the literature by hand (edit NRF2_TARGETS); the Wald statistic is the DESeq2 statistic from the stored object (unshrunk).  D needs decoupleR + OmnipathR + internet:  conda install -c conda-forge -c bioconda bioconductor-decoupler bioconductor-omnipathr
# Without them D is skipped and A-C still run.   Env: GENE (Top1), TREAT (Mock_EBOV), REF (Mock_Mock), SKIP_TF=1, PROJ.   Outputs: requant_v2_outputs/figures_nrf2_<gene>_EBOV/.   Usage: Rscript 132_nrf2_top1_checks.R
# CAUTION: 3 vs 3 replicates; NRF2 activity is controlled mainly at the protein level (KEAP1 degradation, nuclear NRF2), so mRNA targets are a proxy; an enriched target set shows a program, not that NRF2 causes the Top1 change.

suppressPackageStartupMessages({ library(readr); library(dplyr); library(tidyr); library(ggplot2); library(DESeq2) })
select <- dplyr::select; filter <- dplyr::filter; slice <- dplyr::slice
GENE <- Sys.getenv("GENE", "Top1"); TREAT <- Sys.getenv("TREAT", "Mock_EBOV"); REF <- Sys.getenv("REF", "Mock_Mock"); INK <- "#1F2933"; LOW_EXPR <- 50
NRF2_TARGETS <- c("Hmox1", "Nqo1", "Gclc", "Gclm", "Srxn1", "Txnrd1", "Slc7a11", "Gsta3", "Gstm1", "Gstm4", "Gstp1", "Fth1", "Ftl1", "Prdx1", "Txn1", "Sqstm1", "Abcc1", "Me1", "G6pdx", "Pgd", "Cbr1", "Aldh1a1")
CLASSES <- list(Query = GENE, `Canonical NRF2 targets` = c("Hmox1", "Nqo1", "Gclc", "Gclm", "Srxn1", "Txnrd1", "Slc7a11"), `NRF2 pathway` = c("Nfe2l2", "Keap1", "Bach1", "Maff"))
ORD <- c("Mock_Mock", "IFNB_Mock", "IFNY_Mock", "Mock_EBOV", "IFNB_EBOV", "IFNY_EBOV", "Mock_MARV", "IFNB_MARV", "IFNY_MARV", "Mock_EVLP", "Mock_MVLP")
LAB <- c(Mock_Mock = "Mock", IFNB_Mock = "IFNB", IFNY_Mock = "IFNG", Mock_EBOV = "EBOV", IFNB_EBOV = "IFNB+EBOV", IFNY_EBOV = "IFNG+EBOV", Mock_MARV = "MARV", IFNB_MARV = "IFNB+MARV", IFNY_MARV = "IFNG+MARV", Mock_EVLP = "EVLP", Mock_MVLP = "MVLP")
COL <- c(Mock = "#7F7F7F", EBOV = "#1B7837", IFNB = "#2B6CB0", `IFNB+EBOV` = "#9BD08E", IFNG = "#D69E2E", `IFNG+EBOV` = "#E8C77A", MARV = "#8C510A", `IFNB+MARV` = "#C7A27C", `IFNG+MARV` = "#DCC9A8", EVLP = "#6A51A3", MVLP = "#CBA4D8")
PROJ <- Sys.getenv("PROJ", unset = file.path(
  Sys.getenv("HOME"),
  "Documents/aviti_rnaseq/2025/aviti",
  "20251009-100725_Maury_25161_WivfFacICOTHFAUJwSPeHBEFfCfNtCRpUSgVGUkz",
  "results/Maury_25161"))
RNA <- file.path(PROJ, "requant_v2_outputs"); OUT <- file.path(RNA, paste0("figures_nrf2_", tolower(GENE), "_", sub("^Mock_", "", TREAT))); dir.create(OUT, showWarnings = FALSE)
save2 <- function(p, name, w, h) { ggsave(file.path(OUT, paste0(name, ".png")), p, width = w, height = h, dpi = 300); ggsave(file.path(OUT, paste0(name, ".pdf")), p, width = w, height = h) }
theme_t <- theme_bw(base_size = 11) + theme(panel.grid.minor = element_blank(), panel.border = element_rect(colour = "grey70"), axis.text = element_text(colour = INK), plot.title = element_text(colour = INK, face = "bold", size = 11), plot.subtitle = element_text(colour = "grey35", size = 8.5), strip.background = element_rect(fill = "grey95", colour = "grey70"))
stars <- function(p) ifelse(is.na(p), "", ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", "ns"))))
glab <- function(g) ifelse(g %in% names(LAB), LAB[g], g)

dds <- readRDS(file.path(RNA, "dds_v2_decoyAware_hostNormalised.rds")); vsd <- readRDS(file.path(RNA, "vsd_v2.rds")); vm <- SummarizedExperiment::assay(vsd); grp <- as.character(vsd$group)
ann <- read_csv(file.path(RNA, "DE_v2_EBOV_vs_Mock_annotated.csv"), show_col_types = FALSE)
res <- as.data.frame(results(dds, contrast = c("group", TREAT, REF))); res$gene_id <- rownames(res)
rk <- res |> inner_join(ann |> select(gene_id, gene_name) |> filter(!is.na(gene_name), gene_name != ""), by = "gene_id") |> filter(is.finite(stat)) |> group_by(gene_name) |> slice_max(abs(stat), n = 1, with_ties = FALSE) |> ungroup()
ranks <- sort(setNames(rk$stat, rk$gene_name), decreasing = TRUE); lc <- tolower(names(ranks)); message(length(ranks), " genes ranked")
nm2id <- setNames(ann$gene_id, tolower(ann$gene_name)); nm2id <- nm2id[!duplicated(names(nm2id))]

# ---- A: bar graph ----
want <- bind_rows(lapply(names(CLASSES), function(k) tibble(class = k, query = CLASSES[[k]]))) |> mutate(pos = row_number(), key = tolower(query))
D <- want |> left_join(ann |> mutate(key = tolower(gene_name)) |> distinct(key, .keep_all = TRUE) |> select(key, gene = gene_name, baseMean, lfc = log2FoldChange, se = lfcSE, padj), by = "key")
miss <- D$query[is.na(D$gene)]; if (length(miss)) message("Not in the DE table (removed by the pre-filter, or a different symbol): ", paste(miss, collapse = ", "))
D <- D |> filter(!is.na(gene)) |> mutate(class = factor(class, levels = names(CLASSES)), lab = stars(padj), vj = ifelse(lfc >= 0, 0, 1), ys = lfc + ifelse(lfc >= 0, 1, -1) * (se + 0.03), xlab = paste0(gene, ifelse(baseMean < LOW_EXPR, "†", "")))
D$xlab <- factor(D$xlab, levels = unique(D$xlab[order(D$pos)])); write_csv(D |> select(class, gene, baseMean, log2FC = lfc, lfcSE = se, padj), file.path(OUT, "nrf2_targets_log2FC.tsv"))
pA <- ggplot(D, aes(xlab, lfc, fill = class)) + geom_col(width = 0.7) + geom_errorbar(aes(ymin = lfc - se, ymax = lfc + se), width = 0.25, colour = INK, linewidth = 0.35) + geom_text(aes(y = ys, label = lab, vjust = vj), colour = INK, size = 3) + geom_hline(yintercept = 0, colour = INK, linewidth = 0.4) +
  facet_grid(. ~ class, scales = "free_x", space = "free_x") + scale_fill_manual(values = c("#1B7837", "#9BD08E", "#C9CDD3"), guide = "none") + scale_y_continuous(expand = expansion(mult = 0.12)) +
  labs(title = paste0(GENE, " next to canonical NRF2 targets: ", glab(TREAT), " vs Mock"), subtitle = "shrunk log2FC +/- 1 SE; stars = padj (ns = not significant); dagger = baseMean < 50", x = NULL, y = "log2 fold-change vs Mock") + theme_t + theme(axis.text.x = element_text(face = "italic", angle = 45, hjust = 1), panel.grid.major.x = element_blank())
save2(pA, paste0("bars_nrf2_", tolower(GENE), "_", sub("^Mock_", "", TREAT)), max(6, 0.8 * nrow(D) + 2.5), 4.8)

# ---- B: is the NRF2 program up as a whole ----
inT <- lc %in% tolower(NRF2_TARGETS); message(sum(inT), " of ", length(NRF2_TARGETS), " NRF2 target genes are in the ranked list: ", paste(names(ranks)[inT], collapse = ", "))
wt <- tryCatch(wilcox.test(ranks[inT], ranks[!inT])$p.value, error = function(e) NA_real_)
message(sprintf("NRF2 targets: median Wald statistic %.2f vs %.2f for all other genes (Wilcoxon p = %.3g)", median(ranks[inT]), median(ranks[!inT]), wt))
bc <- tibble(rank = seq_along(ranks), gene = names(ranks), stat = as.numeric(ranks), nrf2 = inT, q = lc == tolower(GENE))
pB <- ggplot(bc, aes(rank, stat)) + geom_area(fill = "#E5E7EB") + geom_hline(yintercept = 0, colour = INK, linewidth = 0.3) + geom_segment(data = bc |> filter(nrf2), aes(xend = rank, y = min(bc$stat) * 0.98, yend = min(bc$stat) * 0.8), colour = "#C53030", linewidth = 0.5) +
  geom_vline(data = bc |> filter(q), aes(xintercept = rank), colour = "#1B7837", linewidth = 0.7) + geom_text(data = bc |> filter(nrf2), aes(label = gene, y = stat), size = 2.4, angle = 90, hjust = -0.1, colour = "#C53030") + geom_text(data = bc |> filter(q), aes(label = paste0(gene, " (stat ", sprintf("%.2f", stat), ")"), y = max(bc$stat) * 0.85), colour = "#1B7837", hjust = -0.05, fontface = "italic", size = 3.4) +
  labs(title = paste0("Canonical NRF2 targets on the ranked gene list, ", glab(TREAT), " vs Mock"), subtitle = sprintf("red ticks = NRF2 targets (median stat %.2f vs %.2f for all other genes; Wilcoxon p = %.2g); green line = %s", median(ranks[inT]), median(ranks[!inT]), wt, GENE), x = "gene rank (by Wald statistic; high = up in the infected samples)", y = "Wald statistic") + theme_t
save2(pB, paste0("nrf2_barcode_", sub("^Mock_", "", TREAT)), 10, 4.2)

# ---- C: NRF2 score per sample ----
ids <- nm2id[tolower(NRF2_TARGETS)]; ids <- ids[!is.na(ids) & ids %in% rownames(vm)]
if (length(ids) >= 3) {
  z <- t(scale(t(vm[ids, , drop = FALSE]))); sc <- colMeans(z); gid <- nm2id[tolower(GENE)]
  P <- tibble(sample = colnames(vm), group = grp, nrf2 = sc, gene = as.numeric(vm[gid, ])) |> mutate(g = factor(glab(group), levels = glab(ORD)))
  pC <- ggplot(P, aes(g, nrf2, fill = g)) + geom_boxplot(outlier.shape = NA, alpha = 0.75, width = 0.6, colour = INK, linewidth = 0.3) + geom_point(position = position_jitter(width = 0.08, height = 0, seed = 1), size = 1.8, shape = 21, fill = "white", colour = INK, stroke = 0.4) + scale_fill_manual(values = COL, guide = "none") +
    labs(title = "NRF2 target score per sample", subtitle = paste0("mean z-score of ", length(ids), " canonical NRF2 target genes (VST); every sample shown"), x = NULL, y = "NRF2 target score") + theme_t + theme(axis.text.x = element_text(angle = 40, hjust = 1))
  save2(pC, "nrf2_score_by_group", 8, 4.2)
  rc <- function(a, b) { ok <- is.finite(a) & is.finite(b); if (sum(ok) < 4) return(c(NA, NA)); t <- suppressWarnings(cor.test(a[ok], b[ok], method = "spearman", exact = FALSE)); c(unname(t$estimate), t$p.value) }
  six <- P$group %in% c(TREAT, REF); a <- rc(P$gene, P$nrf2); b <- rc(P$gene[six], P$nrf2[six])
  pD <- ggplot(P, aes(nrf2, gene)) + geom_smooth(aes(group = 1), method = "lm", se = FALSE, colour = "grey50", linewidth = 0.4, formula = y ~ x) + geom_point(aes(fill = g), shape = 21, size = 2.6, colour = INK, stroke = 0.3) + scale_fill_manual(values = COL, name = NULL) +
    annotate("text", x = -Inf, y = Inf, hjust = -0.05, vjust = 1.3, size = 3, colour = INK, label = sprintf("all samples: rho %.2f (p %.2g)\nMock + %s only (n = %d): rho %.2f (p %.2g)", a[1], a[2], glab(TREAT), sum(six), b[1], b[2])) +
    labs(title = paste0(GENE, " vs the NRF2 target score across samples"), subtitle = "each point = one sample; Spearman correlation", x = "NRF2 target score", y = paste0(GENE, " VST")) + theme_t
  save2(pD, paste0(tolower(GENE), "_vs_nrf2_score"), 7.5, 5)
  message(sprintf("%s vs NRF2 score: all samples rho %.2f (p %.2g); Mock + %s only rho %.2f (p %.2g)", GENE, a[1], a[2], glab(TREAT), b[1], b[2]))
} else message("fewer than 3 NRF2 target genes in the VST object - sections C skipped")

# ---- D: decoupleR ----
if (Sys.getenv("SKIP_TF") != "1") tryCatch({
  if (!requireNamespace("decoupleR", quietly = TRUE)) stop("decoupleR is not installed")
  mat <- matrix(as.numeric(ranks), ncol = 1, dimnames = list(names(ranks), glab(TREAT)))
  message("\ngetting the CollecTRI mouse network (needs internet the first time) ..."); net <- decoupleR::get_collectri(organism = "mouse", split_complexes = FALSE)
  r <- as_tibble(as.data.frame(decoupleR::run_ulm(mat, net, .source = "source", .target = "target", .mor = "mor", minsize = 5))); write_tsv(r, file.path(OUT, "tf_activity_EBOV_all.tsv"))
  top <- r |> slice_max(abs(score), n = 25, with_ties = FALSE) |> mutate(src = factor(source, levels = source[order(score)]), sig = stars(p_value)); nr <- r |> filter(tolower(source) == "nfe2l2")
  if (nrow(nr)) message(sprintf("NRF2 (Nfe2l2) activity in %s vs Mock: ULM score %.2f, p = %.3g (rank %d of %d TFs by |score|)", glab(TREAT), nr$score, nr$p_value, match("Nfe2l2", r$source[order(-abs(r$score))]), nrow(r)))
  pt <- ggplot(top, aes(score, src, fill = score > 0)) + geom_col(width = 0.7) + geom_text(aes(label = sig, x = score + sign(score) * 0.15), size = 3, colour = INK) + scale_fill_manual(values = c(`TRUE` = "#C53030", `FALSE` = "#2B6CB0"), guide = "none") +
    labs(title = paste0("Inferred TF activity, ", glab(TREAT), " vs Mock: the 25 strongest"), subtitle = "CollecTRI mouse regulons, univariate linear model on the Wald statistic; red = target genes go up; stars = ULM p", x = "ULM score", y = NULL) + theme_t
  save2(pt, "tf_activity_EBOV_top", 7, 6.2)
  reg <- net |> filter(tolower(target) == tolower(GENE)) |> distinct(source, mor); message("CollecTRI regulators of ", GENE, ": ", if (nrow(reg)) paste(reg$source, collapse = ", ") else "none listed")
  isn <- any(tolower(net$source) == "nfe2l2" & tolower(net$target) == tolower(GENE)); message(GENE, if (isn) " IS" else " is NOT", " listed as an NRF2 (Nfe2l2) target in CollecTRI")
  if (nrow(reg)) { rr <- r |> filter(source %in% reg$source) |> left_join(reg, by = "source") |> arrange(desc(abs(score))); write_tsv(rr, file.path(OUT, paste0("tf_activity_regulators_of_", tolower(GENE), ".tsv"))); message("\nregulators of ", GENE, " and their inferred activity:"); print(as.data.frame(rr |> transmute(TF = source, mor, score = round(score, 2), p = signif(p_value, 2))), row.names = FALSE)
    pr <- ggplot(rr |> mutate(src = factor(source, levels = source[order(score)]), sig = stars(p_value), mode = ifelse(mor > 0, "activator", "repressor")), aes(score, src, fill = score > 0)) + geom_col(width = 0.7) + geom_text(aes(label = paste0(sig, "  (", mode, ")"), x = score + sign(score) * 0.1, hjust = ifelse(score > 0, 0, 1)), size = 2.8, colour = INK) +
      scale_fill_manual(values = c(`TRUE` = "#C53030", `FALSE` = "#2B6CB0"), guide = "none") + labs(title = paste0("Curated regulators of ", GENE, " (CollecTRI): inferred activity, ", glab(TREAT), " vs Mock"), subtitle = "ULM score on the Wald statistic; red = their target genes go up; (activator/repressor) = how CollecTRI lists them acting on the gene", x = "ULM score", y = NULL) + theme_t
    save2(pr, paste0("tf_activity_regulators_of_", tolower(GENE)), 7.5, 1.8 + 0.3 * nrow(rr)) }
  nt <- net |> filter(tolower(source) == "nfe2l2") |> inner_join(tibble(target = names(ranks), stat = as.numeric(ranks)), by = "target")
  if (nrow(nt)) { write_tsv(nt, file.path(OUT, "nrf2_regulon_targets_stat.tsv")); lab <- nt |> slice_max(abs(stat), n = 12, with_ties = FALSE) |> bind_rows(nt |> filter(tolower(target) == tolower(GENE))) |> distinct(target, .keep_all = TRUE)
    pn <- ggplot(nt, aes(stat, fill = stat > 0)) + geom_histogram(bins = 30, colour = "white") + geom_vline(xintercept = 0, colour = INK, linewidth = 0.3) + scale_fill_manual(values = c(`TRUE` = "#C53030", `FALSE` = "#2B6CB0"), guide = "none") +
      labs(title = paste0("Targets of NRF2 (CollecTRI, ", nrow(nt), " genes): Wald statistic, ", glab(TREAT), " vs Mock"), subtitle = paste0("strongest: ", paste(lab$target, collapse = ", ")), x = "Wald statistic (positive = up in the infected samples)", y = "NRF2 target genes") + theme_t
    save2(pn, "nrf2_regulon_targets", 7.5, 4) }
}, error = function(e) message("\nSection D (decoupleR) skipped: ", conditionMessage(e), "\n  install:  conda install -c conda-forge -c bioconda bioconductor-decoupler bioconductor-omnipathr   (and check the internet connection)"))
message("\nSaved to: ", OUT)
