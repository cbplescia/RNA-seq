#!/usr/bin/env Rscript
# Is Top1 part of a program that EBOV turns up (or down)?  Gene-set enrichment analysis (GSEA) of EBOV vs Mock over ALL genes ranked by the DESeq2 Wald statistic (not a cut-off list, so a modest change like Top1's is not lost),
# then every enriched set that CONTAINS Top1 is listed, with whether Top1 sits in the set's leading edge (the genes that drive the enrichment).
# Collections (msigdbr, mouse): Hallmark, GO biological process, GO cellular component (nucleolus, nuclear speckles ...), Reactome.  Method: clusterProfiler::GSEA (permutation / fgsea based), sets of 10-500 genes, BH-adjusted.
# Outputs (requant_v2_outputs/figures_gsea_<gene>_EBOV/): gsea_<collection>.tsv (all sets), gsea_sets_containing_<gene>.tsv, bars_sets_containing_<gene>.png/.pdf (NES of the enriched sets that contain the gene, * = gene is in the leading edge),
#   top_sets_<collection>.png/.pdf (strongest up / down sets per collection, for context), enrichment_<gene>_set<N>.png/.pdf (running-score plots of the top sets containing the gene), rank_position_<gene>.png/.pdf (where the gene sits in the ranked list).
# Needs: clusterProfiler, msigdbr (and enrichplot for the running-score plots):  conda install -c conda-forge -c bioconda bioconductor-clusterprofiler bioconductor-enrichplot r-msigdbr   (msigdbr needs a one-time data download).
# Env: GENE (Top1), TREAT (Mock_EBOV), REF (Mock_Mock), PADJ (0.25: the GSEA convention for 'worth a look'), PROJ.   Usage: Rscript 128_ebov_gsea_top1.R
# CAUTION: with 3 vs 3 replicates and a ~1.2-fold change in Top1, enrichment shows co-movement, not cause. Sets are annotations, and a set that contains Top1 and is enriched is a lead to follow up (is the program transcriptional? which TF?), not a mechanism.

suppressPackageStartupMessages({ library(readr); library(dplyr); library(tidyr); library(ggplot2); library(DESeq2) })
for (pk in c("clusterProfiler", "msigdbr")) if (!requireNamespace(pk, quietly = TRUE)) stop("Install ", pk, ":  conda install -c conda-forge -c bioconda bioconductor-clusterprofiler bioconductor-enrichplot r-msigdbr")
select <- dplyr::select; filter <- dplyr::filter; slice <- dplyr::slice
GENE <- Sys.getenv("GENE", "Top1"); TREAT <- Sys.getenv("TREAT", "Mock_EBOV"); REF <- Sys.getenv("REF", "Mock_Mock"); PADJ <- as.numeric(Sys.getenv("PADJ", "0.25")); INK <- "#1F2933"
PROJ <- Sys.getenv("PROJ", unset = file.path(
  Sys.getenv("HOME"),
  "Documents/aviti_rnaseq/2025/aviti",
  "20251009-100725_Maury_25161_WivfFacICOTHFAUJwSPeHBEFfCfNtCRpUSgVGUkz",
  "results/Maury_25161"))
RNA <- file.path(PROJ, "requant_v2_outputs"); OUT <- file.path(RNA, paste0("figures_gsea_", tolower(GENE), "_", sub("^Mock_", "", TREAT))); dir.create(OUT, showWarnings = FALSE)
save2 <- function(p, name, w, h) { ggsave(file.path(OUT, paste0(name, ".png")), p, width = w, height = h, dpi = 300); ggsave(file.path(OUT, paste0(name, ".pdf")), p, width = w, height = h) }
theme_t <- theme_bw(base_size = 10) + theme(panel.grid.minor = element_blank(), panel.border = element_rect(colour = "grey70"), axis.text = element_text(colour = INK), plot.title = element_text(colour = INK, face = "bold", size = 11), plot.subtitle = element_text(colour = "grey35", size = 8.5), strip.background = element_rect(fill = "grey95", colour = "grey70"))
wrap <- function(x, w = 55) vapply(strwrap(gsub("_", " ", x), w, simplify = FALSE), paste, "", collapse = "\n")

# ---- ranked list: Wald statistic of TREAT vs REF over every tested gene, one value per gene symbol ----
dds <- readRDS(file.path(RNA, "dds_v2_decoyAware_hostNormalised.rds")); res <- as.data.frame(results(dds, contrast = c("group", TREAT, REF))); res$gene_id <- rownames(res)
ann <- read_csv(file.path(RNA, "DE_v2_EBOV_vs_Mock_annotated.csv"), show_col_types = FALSE) |> select(gene_id, gene_name) |> filter(!is.na(gene_name), gene_name != "")
rk <- res |> inner_join(ann, by = "gene_id") |> filter(is.finite(stat)) |> group_by(gene_name) |> slice_max(abs(stat), n = 1, with_ties = FALSE) |> ungroup()
ranks <- sort(setNames(rk$stat, rk$gene_name), decreasing = TRUE); message(length(ranks), " genes ranked (", TREAT, " vs ", REF, ")")
me <- rk |> filter(tolower(gene_name) == tolower(GENE)); if (!nrow(me)) stop(GENE, " is not in the ranked list")
rpos <- match(me$gene_name, names(ranks)); message(sprintf("%s: stat %.2f, log2FC (unshrunk) %.2f, padj %.3g, rank %d of %d (%.1f percentile from the top)", me$gene_name, me$stat, me$log2FoldChange, me$padj, rpos, length(ranks), 100 * rpos / length(ranks)))
gsym <- me$gene_name
pr <- ggplot(tibble(rank = seq_along(ranks), stat = as.numeric(ranks)), aes(rank, stat)) + geom_area(fill = "#C9CDD3") + geom_hline(yintercept = 0, colour = INK, linewidth = 0.3) + geom_vline(xintercept = rpos, colour = "#C53030", linewidth = 0.5) +
  annotate("text", x = rpos, y = max(ranks) * 0.9, label = sprintf(" %s\n rank %d / %d\n stat %.2f", gsym, rpos, length(ranks), me$stat), hjust = -0.05, colour = "#C53030", size = 3.2) +
  labs(title = paste0("Where ", gsym, " sits among all genes: ", sub("^Mock_", "", TREAT), " vs Mock"), subtitle = "genes ranked by the DESeq2 Wald statistic (positive = higher in the infected samples)", x = "gene rank", y = "Wald statistic") + theme_t
save2(pr, paste0("rank_position_", tolower(GENE)), 7.5, 3.6)

# ---- gene sets ----
get_sets <- function(coll, sub = NULL) {
  r <- tryCatch(suppressWarnings(if (is.null(sub)) msigdbr::msigdbr(species = "Mus musculus", collection = coll) else msigdbr::msigdbr(species = "Mus musculus", collection = coll, subcollection = sub)), error = function(e) NULL)
  if (is.null(r) || !nrow(r)) r <- tryCatch(suppressWarnings(if (is.null(sub)) msigdbr::msigdbr(species = "Mus musculus", category = coll) else msigdbr::msigdbr(species = "Mus musculus", category = coll, subcategory = sub)), error = function(e) NULL)
  if (is.null(r) || !nrow(r)) return(NULL); r |> distinct(gs_name, gene_symbol) |> transmute(term = gs_name, gene = gene_symbol) }
COLS <- list(Hallmark = list("H", NULL), GO_BP = list("C5", "GO:BP"), GO_CC = list("C5", "GO:CC"), Reactome = list("C2", "CP:REACTOME"))
allres <- list(); objs <- list()
for (nm in names(COLS)) tryCatch({
  t2g <- get_sets(COLS[[nm]][[1]], COLS[[nm]][[2]]); if (is.null(t2g)) stop("could not load the msigdbr sets"); message("\n", nm, ": ", n_distinct(t2g$term), " sets")
  set.seed(1); g <- clusterProfiler::GSEA(geneList = ranks, TERM2GENE = t2g, minGSSize = 10, maxGSSize = 500, pvalueCutoff = 1, eps = 0, seed = TRUE, verbose = FALSE)
  d <- as_tibble(as.data.frame(g)); if (!nrow(d)) stop("no sets tested")
  inset <- t2g |> filter(tolower(gene) == tolower(gsym)) |> pull(term)
  d <- d |> mutate(collection = nm, contains_gene = ID %in% inset, in_leading_edge = grepl(paste0("(^|/)", gsym, "(/|$)"), core_enrichment, ignore.case = TRUE)) |> arrange(p.adjust)
  allres[[nm]] <- d; objs[[nm]] <- g; write_tsv(d |> select(-any_of("Description")), file.path(OUT, paste0("gsea_", nm, ".tsv")))
  message(sprintf("  %d sets with padj < %.2f (%d up, %d down); %d sets contain %s, of which %d enriched", sum(d$p.adjust < PADJ), PADJ, sum(d$p.adjust < PADJ & d$NES > 0), sum(d$p.adjust < PADJ & d$NES < 0), sum(d$contains_gene), gsym, sum(d$contains_gene & d$p.adjust < PADJ)))
  tp <- bind_rows(d |> filter(NES > 0) |> slice_head(n = 10), d |> filter(NES < 0) |> slice_head(n = 10)) |> mutate(lab = wrap(ID), dir = ifelse(NES > 0, "up in infected", "down in infected"), lab = factor(lab, levels = unique(lab[order(NES)])))
  pt <- ggplot(tp, aes(NES, lab, fill = dir)) + geom_col(width = 0.7) + geom_vline(xintercept = 0, colour = INK, linewidth = 0.3) + geom_point(data = tp |> filter(contains_gene), mapping = aes(x = 0, y = lab), shape = 8, colour = "#C53030", size = 2, inherit.aes = FALSE) +
    scale_fill_manual(values = c(`up in infected` = "#C53030", `down in infected` = "#2B6CB0"), name = NULL) + labs(title = paste0(nm, ": strongest gene sets, ", sub("^Mock_", "", TREAT), " vs Mock"), subtitle = paste0("normalised enrichment score; red star = the set contains ", gsym, "; the 10 strongest in each direction"), x = "NES", y = NULL) + theme_t + theme(axis.text.y = element_text(size = 6.5), legend.position = "bottom")
  save2(pt, paste0("top_sets_", nm), 9, 2.5 + 0.2 * nrow(tp))
}, error = function(e) message("\n", nm, " skipped: ", conditionMessage(e)))
if (!length(allres)) stop("no collection could be analysed")
A <- bind_rows(allres); S <- A |> filter(contains_gene) |> arrange(p.adjust); write_tsv(S |> select(collection, ID, NES, pvalue, p.adjust, setSize, in_leading_edge, core_enrichment), file.path(OUT, paste0("gsea_sets_containing_", tolower(GENE), ".tsv")))
message("\n--- enriched sets that contain ", gsym, " (padj < ", PADJ, "), strongest first ---"); sh <- S |> filter(p.adjust < PADJ); if (nrow(sh)) print(as.data.frame(sh |> transmute(collection, set = substr(ID, 1, 60), NES = round(NES, 2), padj = signif(p.adjust, 2), leading_edge = in_leading_edge) |> slice_head(n = 25)), row.names = FALSE) else message("none - ", gsym, " is not in any enriched set at this threshold (see the .tsv for the nearest)")
pl <- (if (nrow(sh)) sh else S) |> slice_head(n = 20) |> mutate(lab = paste0(wrap(ID, 60), ifelse(in_leading_edge, "  *", "")), lab = factor(lab, levels = rev(unique(lab))), sig = p.adjust < PADJ)
if (nrow(pl)) { pb <- ggplot(pl, aes(NES, lab, fill = collection, alpha = sig)) + geom_col(width = 0.7) + geom_vline(xintercept = 0, colour = INK, linewidth = 0.3) + scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = 0.35), guide = "none") + scale_fill_brewer(palette = "Set2", name = NULL) +
    labs(title = paste0("Gene sets that contain ", gsym, ": ", sub("^Mock_", "", TREAT), " vs Mock"), subtitle = paste0("NES of the 20 best sets containing the gene (faded = padj >= ", PADJ, "); * = ", gsym, " is in the leading edge (a gene that drives the enrichment)"), x = "normalised enrichment score", y = NULL) + theme_t + theme(axis.text.y = element_text(size = 7), legend.position = "bottom")
  save2(pb, paste0("bars_sets_containing_", tolower(GENE)), 9.5, 2.6 + 0.26 * nrow(pl)) }
if (requireNamespace("enrichplot", quietly = TRUE) && nrow(sh)) for (k in seq_len(min(3, nrow(sh)))) tryCatch({ r <- sh[k, ]; p <- enrichplot::gseaplot2(objs[[r$collection]], geneSetID = r$ID, title = paste0(r$collection, ": ", substr(r$ID, 1, 70), "  (NES ", round(r$NES, 2), ", padj ", signif(r$p.adjust, 2), ")"), pvalue_table = FALSE)
  save2(p, paste0("enrichment_", tolower(GENE), "_set", k), 7.5, 5) }, error = function(e) message("running-score plot ", k, " skipped: ", conditionMessage(e)))
message("\nSaved to: ", OUT)
