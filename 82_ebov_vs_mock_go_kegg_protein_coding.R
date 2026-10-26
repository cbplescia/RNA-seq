#!/usr/bin/env Rscript
# GO (BP, MF, CC) and KEGG enrichment of the PROTEIN-CODING genes (1) up and (2) down in EBOV vs Mock - a rewrite of 77 that cannot silently hang.
# What changed from 77: (a) no simplify() on the full GO result (with p-value cut-off 1 that compares thousands of terms and can run for a very long time) - it is applied only
# to terms with p.adjust < 0.05, and only if SIMPLIFY=1; (b) every step prints a timed progress line; (c) results are written after each set, so a failure later does not lose earlier output;
# (d) KEGG gets a download timeout and is skipped (GO continues) if it cannot connect - KEGG needs internet; (e) the input is the protein-coding table DE_v2_EBOV_vs_Mock_protein_coding.csv.
#   up   = padj < PADJ_CUT and shrunk log2FC >  LFC_CUT        down = padj < PADJ_CUT and shrunk log2FC < -LFC_CUT   (default LFC_CUT = 1; LFC_CUT=0 for padj alone)
#   a direction with fewer than MIN_GENES genes at the log2FC cut-off falls back to padj alone (the message says so)
#   background (universe) = every protein-coding gene tested in EBOV vs Mock (with an Entrez ID)
# Writes to requant_v2_outputs/figures_ebov_vs_mock/enrichment_protein_coding/ (PNG + PDF):
#   enrichment_all_terms.csv, enrichment_summary.csv, genes_up.csv, genes_down.csv
#   enrichment_up.png / enrichment_down.png      top terms per ontology (bar = -log10 p, shade = significance, number = genes in the term)
#   compare_<ontology>.png                       up and down side by side (tile = -log10 p.adjust, * = p.adjust < 0.05)
# Usage: Rscript 82_ebov_vs_mock_go_kegg_protein_coding.R   (PROJ from 'source setup_session.sh';  LFC_CUT=0  SIMPLIFY=1  SKIP_KEGG=1 optional)

t0 <- Sys.time(); say <- function(...) { message(sprintf("[%5.1f min] ", as.numeric(difftime(Sys.time(), t0, units = "mins"))), ...); flush.console() }
suppressPackageStartupMessages({ library(readr); library(dplyr); library(tidyr); library(ggplot2) })
for (pkg in c("clusterProfiler", "org.Mm.eg.db"))
  if (!requireNamespace(pkg, quietly = TRUE)) stop("Install ", pkg, " (conda install -c conda-forge -c bioconda bioconductor-", tolower(pkg), ")")
select <- dplyr::select; filter <- dplyr::filter; slice <- dplyr::slice; rename <- dplyr::rename; first <- dplyr::first
OrgDb <- org.Mm.eg.db::org.Mm.eg.db
options(timeout = 90)                                                      # seconds, for the KEGG downloads

PADJ_CUT <- 0.05; LFC_CUT <- as.numeric(Sys.getenv("LFC_CUT", unset = "1")); MIN_GENES <- 10; PADJ_SIG <- 0.05; TOP_PER_ONT <- 12; TOP_COMPARE <- 6; TILE_CAP <- 6
SIMPLIFY <- identical(Sys.getenv("SIMPLIFY"), "1"); SKIP_KEGG <- identical(Sys.getenv("SKIP_KEGG"), "1")
INK <- "#1F2933"
PROJ <- Sys.getenv("PROJ", unset = file.path(
  Sys.getenv("HOME"),
  "Documents/aviti_rnaseq/2025/aviti",
  "20251009-100725_Maury_25161_WivfFacICOTHFAUJwSPeHBEFfCfNtCRpUSgVGUkz",
  "results/Maury_25161"))
RNA <- file.path(PROJ, "requant_v2_outputs")
OUT <- file.path(RNA, "figures_ebov_vs_mock", "enrichment_protein_coding"); dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
save2 <- function(p, name, w, h) { ggsave(file.path(OUT, paste0(name, ".png")), p, width = w, height = h, dpi = 300, limitsize = FALSE); ggsave(file.path(OUT, paste0(name, ".pdf")), p, width = w, height = h, limitsize = FALSE) }
wrap <- function(x, w = 52) vapply(strwrap(x, w, simplify = FALSE), paste, "", collapse = "\n")

# ---- input: protein-coding genes only ----
pcf <- file.path(RNA, "DE_v2_EBOV_vs_Mock_protein_coding.csv")
d <- if (file.exists(pcf)) read_csv(pcf, show_col_types = FALSE) else
  read_csv(file.path(RNA, "DE_v2_EBOV_vs_Mock_annotated.csv"), show_col_types = FALSE) |> filter(gene_biotype == "protein_coding")
d <- d |> filter(!is.na(padj), grepl("^ENSMUSG", gene_id)); d$ens <- sub("\\..*$", "", d$gene_id)
say(sprintf("input: %d protein-coding genes tested in EBOV vs Mock (%s)", nrow(d), if (file.exists(pcf)) basename(pcf) else "annotated table filtered to protein_coding"))
map <- suppressMessages(suppressWarnings(clusterProfiler::bitr(unique(d$ens), "ENSEMBL", "ENTREZID", OrgDb))) |> distinct(ENSEMBL, .keep_all = TRUE)
ent <- function(ens) unique(na.omit(map$ENTREZID[match(ens, map$ENSEMBL)]))
uni <- ent(d$ens); say(sprintf("%d genes have Entrez IDs and form the background", length(uni)))

pick <- function(dir) {
  f <- function(cut) if (dir == "up") d |> filter(padj < PADJ_CUT, log2FoldChange > cut) else d |> filter(padj < PADJ_CUT, log2FoldChange < -cut)
  g <- f(LFC_CUT)
  if (nrow(g) < MIN_GENES && LFC_CUT > 0) { say(sprintf("%s: only %d genes with |log2FC| > %g - using padj < %g alone", dir, nrow(g), LFC_CUT, PADJ_CUT)); g <- f(0) }
  g
}
sets <- list(up = pick("up"), down = pick("down"))
for (s in names(sets)) { say(sprintf("%s: %d protein-coding genes", toupper(s), nrow(sets[[s]]))); write_csv(sets[[s]] |> arrange(padj) |> select(gene_id, gene_name, baseMean, log2FoldChange, padj), file.path(OUT, paste0("genes_", s, ".csv"))) }
if (min(vapply(sets, nrow, 1L)) < 3) stop("A direction has fewer than 3 genes - nothing to enrich (try LFC_CUT=0)")

one <- function(label, fn) {
  say("   ", label, " ...")
  r <- tryCatch(fn(), error = function(e) { say("   ", label, " FAILED: ", conditionMessage(e)); NULL })
  if (is.null(r)) return(NULL)
  n <- tryCatch(nrow(r@result), error = function(e) 0L)
  if (!n) { say("   ", label, ": no terms returned"); return(NULL) }
  r
}
run_set <- function(g) {
  out <- list()
  for (ont in c("BP", "MF", "CC")) {
    r <- one(paste0("GO_", ont), function() clusterProfiler::enrichGO(g, OrgDb = OrgDb, keyType = "ENTREZID", ont = ont, universe = uni, pvalueCutoff = 1, qvalueCutoff = 1, minGSSize = 5, maxGSSize = 500, readable = TRUE))
    if (!is.null(r)) {
      if (ont == "BP" && SIMPLIFY) { r2 <- r; r2@result <- r@result[r@result$p.adjust < PADJ_SIG, ]
        if (nrow(r2@result) > 1) { s <- tryCatch(clusterProfiler::simplify(r2, cutoff = 0.7, by = "p.adjust", select_fun = min), error = function(e) NULL); if (!is.null(s)) r@result <- r@result[!r@result$ID %in% setdiff(r2@result$ID, s@result$ID), ] } }
      out[[paste0("GO_", ont)]] <- r
    }
  }
  if (!SKIP_KEGG) {
    r <- one("KEGG", function() clusterProfiler::enrichKEGG(g, organism = "mmu", universe = uni, pvalueCutoff = 1, qvalueCutoff = 1, minGSSize = 5, maxGSSize = 500))
    if (!is.null(r)) { r <- tryCatch(clusterProfiler::setReadable(r, OrgDb, keyType = "ENTREZID"), error = function(e) r); out[["KEGG"]] <- r }
  }
  out
}
res <- list()
for (s in names(sets)) {
  say(toupper(s), " enrichment")
  o <- run_set(ent(sets[[s]]$ens))
  for (ont in names(o)) {
    r <- o[[ont]]@result |> mutate(set = s, ontology = ont); res[[paste(s, ont)]] <- r
    say(sprintf("   %-7s %4d terms tested, %d with p.adjust < %g%s", ont, nrow(r), sum(r$p.adjust < PADJ_SIG), PADJ_SIG,
                if (any(r$p.adjust < PADJ_SIG)) paste0("  (", paste(head(r$Description[order(r$pvalue)][r$p.adjust[order(r$pvalue)] < PADJ_SIG], 2), collapse = " | "), ")") else ""))
  }
  if (length(res)) write_csv(bind_rows(res) |> select(set, ontology, ID, Description, GeneRatio, BgRatio, Count, pvalue, p.adjust, geneID), file.path(OUT, "enrichment_all_terms.csv"))   # saved after each set
}
if (!length(res)) stop("No enrichment results were produced - see the FAILED lines above (for KEGG, check the internet connection; GO needs only the local org.Mm.eg.db)")
R <- bind_rows(res)
write_csv(R |> group_by(set, ontology) |> summarise(terms_tested = n(), terms_padj_lt_0.05 = sum(p.adjust < PADJ_SIG), best_p = min(pvalue), .groups = "drop"), file.path(OUT, "enrichment_summary.csv"))

say("drawing figures")
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
    labs(title = sprintf("EBOV vs Mock: protein-coding genes %s-regulated", s), subtitle = sprintf("%d genes; top %d terms per ontology by p-value; number = genes in the term; %s", nrow(sets[[s]]), TOP_PER_ONT,
         if (n_sig) sprintf("%d terms reach p.adjust < 0.05", n_sig) else "NO term reaches p.adjust < 0.05 (exploratory)"), x = expression(-log[10](p)), y = NULL) + theme_f()
  save2(pl, paste0("enrichment_", s), 10, 0.2 * nrow(top) + 2.8)
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
  save2(pc, paste0("compare_", ont), 8, 0.28 * n_distinct(tl$term) + 2.2)
}
say("done. Saved to: ", OUT)
