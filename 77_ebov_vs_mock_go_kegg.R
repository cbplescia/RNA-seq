#!/usr/bin/env Rscript
# GO (BP, MF, CC) and KEGG enrichment of the genes UP and DOWN in EBOV vs Mock (untreated).
#   up   = padj < PADJ_CUT and shrunk log2FC >  LFC_CUT        down = padj < PADJ_CUT and shrunk log2FC < -LFC_CUT   (protein-coding genes)
#   background (universe) = every protein-coding gene tested in the EBOV vs Mock comparison. If a direction has fewer than MIN_GENES genes at the log2FC cut-off, the cut-off is
#   dropped for that direction (padj only) and the message says so; set LFC_CUT=0 to use padj only for both. KEGG needs internet; GO BP is also reduced with simplify(0.7) when possible.
# Writes to requant_v2_outputs/figures_ebov_vs_mock/enrichment/ (PNG + PDF):
#   enrichment_EBOV_vs_Mock_all_terms.csv      every term tested (set, ontology, ID, Description, GeneRatio, BgRatio, Count, pvalue, p.adjust, geneID)
#   enrichment_EBOV_vs_Mock_summary.csv        terms tested and terms at p.adjust < 0.05 per set and ontology
#   genes_EBOV_vs_Mock_up.csv / _down.csv      the gene lists used
#   enrichment_EBOV_vs_Mock_up.png / _down.png top TOP_PER_ONT terms per ontology (bar = -log10 p, shade = significance, number = genes in the term)
#   enrichment_EBOV_vs_Mock_compare_<ontology>.png   the top terms of up and down side by side (tile = -log10 p.adjust, * = p.adjust < 0.05)
# Always read the shading: a set with no term at p.adjust < 0.05 is drawn with its lowest raw p-values and marked exploratory.
# Usage: Rscript 77_ebov_vs_mock_go_kegg.R   (PROJ from 'source setup_session.sh';  LFC_CUT=0 optional)

suppressPackageStartupMessages({ library(readr); library(dplyr); library(tidyr); library(ggplot2) })
for (pkg in c("clusterProfiler", "org.Mm.eg.db"))
  if (!requireNamespace(pkg, quietly = TRUE)) stop("Install ", pkg, " (conda install -c conda-forge -c bioconda bioconductor-", tolower(pkg), ")")
select <- dplyr::select; filter <- dplyr::filter; slice <- dplyr::slice; rename <- dplyr::rename; first <- dplyr::first
OrgDb <- org.Mm.eg.db::org.Mm.eg.db

PADJ_CUT <- 0.05; LFC_CUT <- as.numeric(Sys.getenv("LFC_CUT", unset = "1")); MIN_GENES <- 10; PADJ_SIG <- 0.05; TOP_PER_ONT <- 12; TOP_COMPARE <- 6; TILE_CAP <- 6
INK <- "#1F2933"
PROJ <- Sys.getenv("PROJ", unset = file.path(
  Sys.getenv("HOME"),
  "Documents/aviti_rnaseq/2025/aviti",
  "20251009-100725_Maury_25161_WivfFacICOTHFAUJwSPeHBEFfCfNtCRpUSgVGUkz",
  "results/Maury_25161"))
RNA <- file.path(PROJ, "requant_v2_outputs")
OUT <- file.path(RNA, "figures_ebov_vs_mock", "enrichment"); dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
save2 <- function(p, name, w, h) { ggsave(file.path(OUT, paste0(name, ".png")), p, width = w, height = h, dpi = 300, limitsize = FALSE); ggsave(file.path(OUT, paste0(name, ".pdf")), p, width = w, height = h, limitsize = FALSE) }
wrap <- function(x, w = 52) vapply(strwrap(x, w, simplify = FALSE), paste, "", collapse = "\n")

d <- read_csv(file.path(RNA, "DE_v2_EBOV_vs_Mock_annotated.csv"), show_col_types = FALSE) |> filter(!is.na(padj), grepl("^ENSMUSG", gene_id), gene_biotype == "protein_coding")
d$ens <- sub("\\..*$", "", d$gene_id)
map <- suppressMessages(suppressWarnings(clusterProfiler::bitr(unique(d$ens), "ENSEMBL", "ENTREZID", OrgDb))) |> distinct(ENSEMBL, .keep_all = TRUE)
ent <- function(ens) unique(na.omit(map$ENTREZID[match(ens, map$ENSEMBL)]))
uni <- ent(d$ens)
message(sprintf("EBOV vs Mock: %d protein-coding genes tested (%d with Entrez IDs, used as the background)", nrow(d), length(uni)))

pick <- function(dir) {
  f <- function(cut) if (dir == "up") d |> filter(padj < PADJ_CUT, log2FoldChange > cut) else d |> filter(padj < PADJ_CUT, log2FoldChange < -cut)
  g <- f(LFC_CUT)
  if (nrow(g) < MIN_GENES && LFC_CUT > 0) { message(sprintf("   %s: only %d genes with |log2FC| > %g - using padj < %g alone", dir, nrow(g), LFC_CUT, PADJ_CUT)); g <- f(0) }
  g
}
sets <- list(up = pick("up"), down = pick("down"))
for (s in names(sets)) { message(sprintf("%s: %d genes", toupper(s), nrow(sets[[s]]))); write_csv(sets[[s]] |> arrange(padj) |> select(gene_id, gene_name, baseMean, log2FoldChange, padj), file.path(OUT, paste0("genes_EBOV_vs_Mock_", s, ".csv"))) }

run_set <- function(g) {
  out <- list()
  add <- function(label, fn) {
    r <- tryCatch(fn(), error = function(e) { message("   ", label, ": not run (", conditionMessage(e), ")"); NULL })
    if (!is.null(r) && nrow(r@result)) {
      if (label == "GO_BP") r <- tryCatch(clusterProfiler::simplify(r, cutoff = 0.7, by = "p.adjust", select_fun = min), error = function(e) r)
      r <- tryCatch(clusterProfiler::setReadable(r, OrgDb, keyType = "ENTREZID"), error = function(e) r)
      out[[label]] <<- r
    }
  }
  for (ont in c("BP", "MF", "CC"))
    add(paste0("GO_", ont), function() clusterProfiler::enrichGO(g, OrgDb = OrgDb, keyType = "ENTREZID", ont = ont, universe = uni, pvalueCutoff = 1, qvalueCutoff = 1, minGSSize = 5, maxGSSize = 500))
  add("KEGG", function() clusterProfiler::enrichKEGG(g, organism = "mmu", universe = uni, pvalueCutoff = 1, qvalueCutoff = 1, minGSSize = 5, maxGSSize = 500))
  out
}
res <- list()
for (s in names(sets)) {
  message("\n", toupper(s))
  o <- run_set(ent(sets[[s]]$ens))
  for (ont in names(o)) {
    r <- o[[ont]]@result |> mutate(set = s, ontology = ont); res[[paste(s, ont)]] <- r
    message(sprintf("   %-7s %4d terms tested, %d with p.adjust < %g%s", ont, nrow(r), sum(r$p.adjust < PADJ_SIG), PADJ_SIG,
                    if (any(r$p.adjust < PADJ_SIG)) paste0("  (", paste(head(r$Description[order(r$pvalue)][r$p.adjust[order(r$pvalue)] < PADJ_SIG], 2), collapse = " | "), ")") else ""))
  }
}
if (!length(res)) stop("No enrichment results (check the gene lists and, for KEGG, the internet connection)")
R <- bind_rows(res)
write_csv(R |> select(set, ontology, ID, Description, GeneRatio, BgRatio, Count, pvalue, p.adjust, geneID), file.path(OUT, "enrichment_EBOV_vs_Mock_all_terms.csv"))
write_csv(R |> group_by(set, ontology) |> summarise(terms_tested = n(), terms_padj_lt_0.05 = sum(p.adjust < PADJ_SIG), best_p = min(pvalue), .groups = "drop"), file.path(OUT, "enrichment_EBOV_vs_Mock_summary.csv"))

theme_f <- function() theme_bw(base_size = 10) + theme(panel.grid.minor = element_blank(), panel.border = element_rect(colour = "grey70"), axis.text = element_text(colour = INK),
  strip.background = element_rect(fill = "grey95", colour = "grey70"), strip.text = element_text(colour = INK, face = "bold"), plot.title = element_text(colour = INK, face = "bold"),
  plot.subtitle = element_text(colour = "grey35", size = 8.5), legend.position = "bottom")
SIGCOL <- c(`p.adjust < 0.05` = "#08519C", `p.adjust < 0.25` = "#6BAED6", `nominal p only` = "#C6C6C6")
for (s in unique(R$set)) {
  r <- R |> filter(set == s)
  top <- r |> group_by(ontology) |> slice_min(pvalue, n = TOP_PER_ONT, with_ties = FALSE) |> ungroup() |>
    mutate(sig = factor(case_when(p.adjust < PADJ_SIG ~ "p.adjust < 0.05", p.adjust < 0.25 ~ "p.adjust < 0.25", TRUE ~ "nominal p only"), levels = names(SIGCOL)),
           y = factor(paste(ontology, ID), levels = rev(paste(ontology, ID))), ontology = factor(ontology, levels = unique(ontology)))
  nm <- setNames(top$Description, as.character(top$y)); n_sig <- sum(r$p.adjust < PADJ_SIG)
  pl <- ggplot(top, aes(-log10(pvalue), y, fill = sig)) + geom_col(width = 0.75) + geom_text(aes(label = Count), hjust = -0.3, size = 2.8, colour = INK) +
    facet_grid(ontology ~ ., scales = "free_y", space = "free_y") + scale_y_discrete(labels = function(x) wrap(nm[x])) +
    scale_fill_manual(values = SIGCOL, name = NULL, drop = FALSE) + scale_x_continuous(expand = expansion(mult = c(0, 0.12))) +
    labs(title = sprintf("EBOV vs Mock: genes %s-regulated", s), subtitle = sprintf("%d genes; top %d terms per ontology by p-value; number = genes in the term; %s", nrow(sets[[s]]), TOP_PER_ONT,
                                                                                 if (n_sig) sprintf("%d terms reach p.adjust < 0.05", n_sig) else "NO term reaches p.adjust < 0.05 (exploratory)"), x = expression(-log[10](p)), y = NULL) + theme_f()
  save2(pl, paste0("enrichment_EBOV_vs_Mock_", s), 10, 0.2 * nrow(top) + 2.8)
}
for (ont in unique(R$ontology)) {
  r <- R |> filter(ontology == ont)
  if (n_distinct(r$set) < 2) next
  keep <- r |> group_by(set) |> slice_min(pvalue, n = TOP_COMPARE, with_ties = FALSE) |> ungroup() |> distinct(ID)
  tl <- r |> filter(ID %in% keep$ID) |> mutate(score = pmin(-log10(p.adjust), TILE_CAP), star = ifelse(p.adjust < PADJ_SIG, "*", ""), set = factor(set, levels = c("up", "down")))
  ord <- tl |> group_by(ID, Description) |> summarise(best = min(pvalue), .groups = "drop") |> arrange(desc(best))
  tl <- tl |> mutate(term = factor(Description, levels = unique(ord$Description)))
  pc <- ggplot(tl, aes(set, term, fill = score)) + geom_tile(colour = "white") + geom_text(aes(label = star), size = 5, colour = "white") +
    scale_fill_gradient(low = "#EAF1F8", high = "#08306B", limits = c(0, TILE_CAP), name = expression(-log[10](p.adjust))) + scale_y_discrete(labels = function(x) wrap(x, 55)) +
    labs(title = sprintf("EBOV vs Mock: %s, up vs down", ont), subtitle = "tile = -log10 p.adjust (capped), * = p.adjust < 0.05; top terms of each direction", x = NULL, y = NULL) + theme_f()
  save2(pc, paste0("enrichment_EBOV_vs_Mock_compare_", ont), 8, 0.28 * n_distinct(tl$term) + 2.2)
}
message("\nSaved to: ", OUT)
