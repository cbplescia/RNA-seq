#!/usr/bin/env Rscript
# GO (BP, MF, CC) and KEGG enrichment of every pattern from 48_patterns_virus_vs_ifnb.R, with visualisations.
# Reads figures_patterns/genes_P*.csv and master_table.csv (run 48 first). Background (universe) = all protein-coding genes
# tested. Patterns need at least MIN_GENES genes. KEGG (and Reactome, if ReactomePA is installed) need internet.
#
# Writes to <project>/requant_v2_outputs/figures_patterns/enrichment/ (PNG + PDF):
#   enrichment_all_patterns.csv                every term tested in every pattern (pvalue, p.adjust, genes, ontology)
#   overview_<pattern>.png                     top 12 terms per ontology: bar = -log10 p, shade = significance, number = genes
#   compare_GO_BP.png / compare_KEGG.png (...) the top terms of every pattern side by side (tile = -log10 p.adjust, * = p.adjust < 0.05)
#   network_<pattern>.png                      gene-concept network for patterns with significant GO BP terms
#   enrichment_summary.csv                     how many terms reach p.adjust < 0.05 per pattern and ontology
# Always read the shading: a pattern with no term at p.adjust < 0.05 gets its lowest raw p-values drawn, marked exploratory.
# Usage: Rscript 49_pattern_enrichment.R   (PROJ from 'source setup_session.sh')

suppressPackageStartupMessages({ library(readr); library(dplyr); library(tidyr); library(ggplot2) })
for (pkg in c("clusterProfiler", "enrichplot", "org.Mm.eg.db"))
  if (!requireNamespace(pkg, quietly = TRUE)) stop("Install ", pkg, " (conda install -c conda-forge -c bioconda bioconductor-", tolower(pkg), ")")
select <- dplyr::select; first <- dplyr::first; rename <- dplyr::rename
slice <- dplyr::slice; filter <- dplyr::filter; desc <- dplyr::desc; count <- dplyr::count
OrgDb <- org.Mm.eg.db::org.Mm.eg.db

MIN_GENES <- 10; PADJ_SIG <- 0.05; TOP_PER_PATTERN <- 5; TOP_OVERVIEW <- 12; TILE_CAP <- 6
INK <- "#1F2933"
LABELS <- c(P1 = "IFNB-up, suppressed by EBOV and MARV", P2 = "IFNB-down, relieved by EBOV and MARV",
            P3u = "IFNB-up, EBOV fails to suppress", P3d = "IFNB-down, EBOV fails to relieve",
            P4u = "IFNB-up, MARV fails to suppress", P4d = "IFNB-down, MARV fails to relieve",
            P5 = "altered by EBOV, not IFNB", P6 = "altered by MARV, not IFNB", P7 = "altered by both viruses, not IFNB")

PROJ <- Sys.getenv("PROJ", unset = file.path(
  Sys.getenv("HOME"),
  "Documents/aviti_rnaseq/2025/aviti",
  "20251009-100725_Maury_25161_WivfFacICOTHFAUJwSPeHBEFfCfNtCRpUSgVGUkz",
  "results/Maury_25161"))
PAT <- file.path(PROJ, "requant_v2_outputs", "figures_patterns")
OUT <- file.path(PAT, "enrichment"); dir.create(OUT, showWarnings = FALSE)
stringr_wrap <- function(x, width) vapply(x, function(s) paste(strwrap(s, width), collapse = "\n"), character(1), USE.NAMES = FALSE)
save2 <- function(p, name, w, h) {
  ggsave(file.path(OUT, paste0(name, ".png")), p, width = w, height = h, dpi = 300, limitsize = FALSE)
  ggsave(file.path(OUT, paste0(name, ".pdf")), p, width = w, height = h, limitsize = FALSE)
}

mt <- file.path(PAT, "master_table.csv")
if (!file.exists(mt)) stop("Run 48_patterns_virus_vs_ifnb.R first (missing ", mt, ")")
M <- read_csv(mt, show_col_types = FALSE)
files <- list.files(PAT, pattern = "^genes_P.*\\.csv$", full.names = TRUE)
if (!length(files)) stop("No genes_P*.csv files found in ", PAT)
pats <- lapply(files, function(f) { d <- read_csv(f, show_col_types = FALSE)
  id <- sub("^genes_|\\.csv$", "", basename(f)); key <- sub("_.*$", "", id)
  list(id = id, key = key, label = paste0(key, ": ", ifelse(key %in% names(LABELS), LABELS[key], id)), genes = d) })
pats <- Filter(function(p) nrow(p$genes) >= MIN_GENES, pats)
if (!length(pats)) stop("No pattern has at least ", MIN_GENES, " genes")
message("Patterns analysed: ", paste(vapply(pats, function(p) sprintf("%s (%d)", p$key, nrow(p$genes)), character(1)), collapse = ", "))

map <- suppressMessages(suppressWarnings(clusterProfiler::bitr(unique(M$gene_id), "ENSEMBL", "ENTREZID", OrgDb))) |> distinct(ENSEMBL, .keep_all = TRUE)
ent <- function(ens) unique(na.omit(map$ENTREZID[match(ens, map$ENSEMBL)]))
uni <- ent(M$gene_id)
message("Background: ", length(uni), " protein-coding genes with Entrez IDs")

# ---- enrichment ---------------------------------------------------------------------------------------------------------
run_set <- function(g) {
  out <- list()
  add <- function(label, fn) {
    r <- tryCatch(fn(), error = function(e) { message("   ", label, ": not run (", conditionMessage(e), ")"); NULL })
    if (!is.null(r) && nrow(r@result)) {
      r <- tryCatch(clusterProfiler::setReadable(r, OrgDb, keyType = "ENTREZID"), error = function(e) r)
      out[[label]] <<- r
    }
  }
  for (ont in c("BP", "MF", "CC"))
    add(paste0("GO_", ont), function() clusterProfiler::enrichGO(g, OrgDb = OrgDb, keyType = "ENTREZID", ont = ont, universe = uni,
                                                                  pvalueCutoff = 1, qvalueCutoff = 1, minGSSize = 5, maxGSSize = 500))
  add("KEGG", function() clusterProfiler::enrichKEGG(g, organism = "mmu", universe = uni, pvalueCutoff = 1, qvalueCutoff = 1, minGSSize = 5, maxGSSize = 500))
  if (requireNamespace("ReactomePA", quietly = TRUE))
    add("Reactome", function() ReactomePA::enrichPathway(g, organism = "mouse", universe = uni, pvalueCutoff = 1, qvalueCutoff = 1, minGSSize = 5, maxGSSize = 500))
  out
}
res <- list(); objs <- list()
for (p in pats) {
  message("\n", p$label, "  (", nrow(p$genes), " genes)")
  objs[[p$id]] <- run_set(ent(p$genes$gene_id))
  for (ont in names(objs[[p$id]])) {
    r <- objs[[p$id]][[ont]]@result |> mutate(pattern = p$id, key = p$key, pattern_label = p$label, ontology = ont)
    res[[paste(p$id, ont)]] <- r
    message(sprintf("   %-9s %4d terms tested, %d with p.adjust < %g%s", ont, nrow(r), sum(r$p.adjust < PADJ_SIG), PADJ_SIG,
                    if (any(r$p.adjust < PADJ_SIG)) paste0("  (", paste(head(r$Description[order(r$pvalue)][r$p.adjust[order(r$pvalue)] < PADJ_SIG], 2), collapse = " | "), ")") else ""))
  }
}
R <- bind_rows(res)
write_csv(R |> select(pattern, ontology, ID, Description, GeneRatio, BgRatio, Count, pvalue, p.adjust, geneID), file.path(OUT, "enrichment_all_patterns.csv"))
write_csv(R |> group_by(pattern, ontology) |> summarise(terms_tested = n(), terms_padj_lt_0.05 = sum(p.adjust < PADJ_SIG), best_p = min(pvalue), .groups = "drop"),
          file.path(OUT, "enrichment_summary.csv"))

theme_f <- function() theme_bw(base_size = 10) + theme(panel.grid.minor = element_blank(), panel.border = element_rect(colour = "grey70"),
  axis.text = element_text(colour = INK), strip.background = element_rect(fill = "grey95", colour = "grey70"),
  strip.text = element_text(colour = INK, face = "bold"), plot.title = element_text(colour = INK, face = "bold"),
  plot.subtitle = element_text(colour = "grey35", size = 8.5), legend.position = "bottom")

# ---- per-pattern overview -------------------------------------------------------------------------------------------------
SIGCOL <- c(`p.adjust < 0.05` = "#08519C", `p.adjust < 0.25` = "#6BAED6", `nominal p only` = "#C6C6C6")
for (p in pats) {
  r <- R |> filter(pattern == p$id)
  if (!nrow(r)) next
  top <- r |> group_by(ontology) |> slice_min(pvalue, n = TOP_OVERVIEW, with_ties = FALSE) |> ungroup() |>
    mutate(sig = factor(case_when(p.adjust < PADJ_SIG ~ "p.adjust < 0.05", p.adjust < 0.25 ~ "p.adjust < 0.25", TRUE ~ "nominal p only"), levels = names(SIGCOL)),
           y = factor(paste(ontology, ID), levels = rev(paste(ontology, ID))), ontology = factor(ontology, levels = unique(ontology)))
  names_map <- setNames(top$Description, as.character(top$y))
  n_sig <- sum(r$p.adjust < PADJ_SIG)
  pl <- ggplot(top, aes(-log10(pvalue), y, fill = sig)) +
    geom_col(width = 0.75) + geom_text(aes(label = Count), hjust = -0.3, size = 2.8, colour = INK) +
    facet_grid(ontology ~ ., scales = "free_y", space = "free_y") +
    scale_y_discrete(labels = function(x) stringr_wrap(names_map[x], 52)) +
    scale_fill_manual(values = SIGCOL, name = NULL, drop = FALSE) + scale_x_continuous(expand = expansion(mult = c(0, 0.12))) +
    labs(title = p$label, subtitle = sprintf("%d genes; top %d terms per ontology by p-value; number = genes in the term; %s", nrow(p$genes), TOP_OVERVIEW,
                                             if (n_sig) sprintf("%d terms reach p.adjust < 0.05", n_sig) else "NO term reaches p.adjust < 0.05 (exploratory)"),
         x = expression(-log[10](p)), y = NULL) + theme_f()
  save2(pl, paste0("overview_", p$id), 10, 0.2 * nrow(top) + 2.8)
}

# ---- comparison across patterns ---------------------------------------------------------------------------------------------
for (ont in unique(R$ontology)) {
  r <- R |> filter(ontology == ont)
  if (n_distinct(r$pattern) < 2) next
  keep <- r |> group_by(pattern) |> slice_min(pvalue, n = TOP_PER_PATTERN, with_ties = FALSE) |> ungroup() |> distinct(ID)
  tl <- r |> filter(ID %in% keep$ID) |> mutate(score = pmin(-log10(p.adjust), TILE_CAP), star = ifelse(p.adjust < PADJ_SIG, "*", ""))
  ord <- tl |> group_by(ID, Description) |> summarise(best = min(pvalue), .groups = "drop") |> arrange(desc(best))
  tl$term <- factor(tl$Description, levels = unique(ord$Description))
  tl$key <- factor(tl$key, levels = unique(vapply(pats, function(p) p$key, character(1))))
  pc <- ggplot(tl, aes(key, term, fill = score)) +
    geom_tile(colour = "white", linewidth = 0.4) + geom_text(aes(label = star), size = 4, colour = "white", vjust = 0.75) +
    scale_fill_gradient(low = "#EEF3FA", high = "#08306B", limits = c(0, TILE_CAP), name = expression(-log[10](p.adjust))) +
    scale_y_discrete(labels = function(x) stringr_wrap(x, 55)) +
    labs(title = paste0("Enrichment compared across patterns: ", ont),
         subtitle = sprintf("the %d lowest-p terms of each pattern; tile = -log10 p.adjust (capped at %d); * = p.adjust < %g; grey = term not tested for that pattern",
                            TOP_PER_PATTERN, TILE_CAP, PADJ_SIG), x = NULL, y = NULL) + theme_f() +
    theme(panel.background = element_rect(fill = "grey90"), panel.grid = element_blank())
  save2(pc, paste0("compare_", ont), 9, 0.22 * nlevels(tl$term) + 3)
}

# ---- gene-concept networks ------------------------------------------------------------------------------------------------
fc_col <- c(P1 = "lfc_B", P2 = "lfc_B", P3u = "lfc_B", P3d = "lfc_B", P4u = "lfc_B", P4d = "lfc_B", P5 = "lfc_E", P6 = "lfc_M", P7 = "lfc_E")
for (p in pats) {
  ego <- objs[[p$id]][["GO_BP"]]
  if (is.null(ego) || sum(ego@result$p.adjust < PADJ_SIG) < 1) next
  ego@result <- ego@result[ego@result$p.adjust < PADJ_SIG, ]
  col <- fc_col[[p$key]]; fc <- setNames(p$genes[[col]], p$genes$gene_name); fc <- fc[!is.na(names(fc)) & !duplicated(names(fc))]
  pn <- tryCatch(enrichplot::cnetplot(ego, showCategory = 5, foldChange = fc, colorEdge = TRUE) + ggtitle(paste0(p$label, " - GO BP (colour = ", sub("lfc_", "", col), " log2FC)")),
                 error = function(e) { message("  network for ", p$key, " not drawn (", conditionMessage(e), ")"); NULL })
  if (!is.null(pn)) save2(pn, paste0("network_", p$id), 10, 8)
}
message("\nEnrichment tables and figures in: ", OUT)
