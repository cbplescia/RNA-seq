#!/usr/bin/env Rscript
# Visuals of the signalling / molecular functions that are enriched among the genes UP in EBOV vs Mock. Reads the output of 82_ebov_vs_mock_go_kegg_protein_coding.R
# (enrichment_all_terms.csv - run 82 first) and the EBOV vs Mock DE table, and makes:
#   1. GO_MF_up_dotplot             the enriched GO molecular-function terms (x = gene ratio, size = genes, colour = p.adjust)
#   2. GO_MF_up_term_gene_map       terms (columns) x genes (rows): a tile = the gene is in the term, coloured by the gene's EBOV log2FC; shows which genes drive which function and which are shared
#   3. GO_MF_up_term_network        terms as nodes (size = genes), edges = shared genes (Jaccard >= JACCARD) - needs the igraph package; skipped if it is missing
#   4. KEGG pathway maps            the KEGG signalling pathways enriched among the up genes drawn as the real KEGG diagrams with every gene coloured by EBOV log2FC (pathview;
#                                   needs the pathview package and internet) - PNG per pathway in pathview/
#   GO_MF_up_terms_genes.csv        term -> genes with log2FC / padj, for the figures above
# GO_MF terms are molecular FUNCTIONS (e.g. kinase activity, cytokine receptor binding), not pathways: 1-3 show which functions / genes they are, 4 places the genes on the actual signalling pathways.
# Terms with p.adjust < PADJ_SIG are used; if none reach it, the TOP_FALLBACK lowest-p terms are shown and labelled exploratory.
# Writes to requant_v2_outputs/figures_ebov_vs_mock/signalling_visuals/ (PNG + PDF).   Usage: Rscript 84_ebov_up_signalling_visuals.R   (PROJ from 'source setup_session.sh')

suppressPackageStartupMessages({ library(readr); library(dplyr); library(tidyr); library(ggplot2) })
select <- dplyr::select; filter <- dplyr::filter; slice <- dplyr::slice

PADJ_SIG <- 0.05; TOP_FALLBACK <- 15; MAX_TERMS <- 25; MAX_GENES <- 70; JACCARD <- 0.25; N_PATH <- 8
INK <- "#1F2933"
PROJ <- Sys.getenv("PROJ", unset = file.path(
  Sys.getenv("HOME"),
  "Documents/aviti_rnaseq/2025/aviti",
  "20251009-100725_Maury_25161_WivfFacICOTHFAUJwSPeHBEFfCfNtCRpUSgVGUkz",
  "results/Maury_25161"))
RNA <- file.path(PROJ, "requant_v2_outputs")
IN <- file.path(RNA, "figures_ebov_vs_mock", "enrichment_protein_coding", "enrichment_all_terms.csv")
OUT <- file.path(RNA, "figures_ebov_vs_mock", "signalling_visuals"); dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
if (!file.exists(IN)) stop("Run 82_ebov_vs_mock_go_kegg_protein_coding.R first (missing ", IN, ")")
save2 <- function(p, name, w, h) { ggsave(file.path(OUT, paste0(name, ".png")), p, width = w, height = h, dpi = 300, limitsize = FALSE); ggsave(file.path(OUT, paste0(name, ".pdf")), p, width = w, height = h, limitsize = FALSE) }
wrap <- function(x, w = 45) vapply(strwrap(x, w, simplify = FALSE), paste, "", collapse = "\n")
theme_f <- function() theme_bw(base_size = 10) + theme(panel.grid.minor = element_blank(), panel.border = element_rect(colour = "grey70"), axis.text = element_text(colour = INK),
  plot.title = element_text(colour = INK, face = "bold"), plot.subtitle = element_text(colour = "grey35", size = 8.5))

E <- read_csv(file.path(RNA, "DE_v2_EBOV_vs_Mock_annotated.csv"), show_col_types = FALSE) |> transmute(gene_name, gene_id, lfc = log2FoldChange, padj)
R <- read_csv(IN, show_col_types = FALSE) |> filter(set == "up")
mf <- R |> filter(ontology == "GO_MF")
if (!nrow(mf)) message("No GO_MF terms for the up set in ", IN, " (GO_MF failed or nothing was tested) - skipping figures 1-3")
exploratory <- FALSE
if (nrow(mf)) {
  sig <- mf |> filter(p.adjust < PADJ_SIG)
  if (!nrow(sig)) { sig <- mf |> arrange(pvalue) |> slice_head(n = TOP_FALLBACK); exploratory <- TRUE; message("No GO_MF term reaches p.adjust < ", PADJ_SIG, " - showing the ", nrow(sig), " lowest p-values (exploratory)") }
  sig <- sig |> arrange(pvalue) |> slice_head(n = MAX_TERMS) |> mutate(ratio = vapply(strsplit(GeneRatio, "/"), function(x) as.numeric(x[1]) / as.numeric(x[2]), numeric(1)))
  message(sprintf("GO_MF terms used: %d%s", nrow(sig), if (exploratory) " (exploratory)" else ""))
  expl <- if (exploratory) "exploratory: no term reaches p.adjust < 0.05" else sprintf("p.adjust < %g", PADJ_SIG)
  # ---- 1. dotplot ----
  sig_dot <- sig |> mutate(Description = factor(Description, levels = rev(Description)))
  p1 <- ggplot(sig_dot, aes(ratio, Description, size = Count, colour = p.adjust)) + geom_point() + scale_colour_gradient(low = "#C53030", high = "#9ECAE1", name = "p.adjust") +
    scale_size_continuous(range = c(2, 7), name = "genes") + scale_y_discrete(labels = function(x) wrap(x, 50)) +
    labs(title = "EBOV vs Mock, up-regulated genes: enriched GO molecular functions", subtitle = expl, x = "gene ratio", y = NULL) + theme_f()
  save2(p1, "GO_MF_up_dotplot", 8.5, 0.32 * nrow(sig) + 2)
  # ---- term -> gene table ----
  TG <- sig |> select(ID, Description, pvalue, p.adjust, geneID) |> separate_rows(geneID, sep = "/") |> rename(gene_name = geneID) |> left_join(E |> distinct(gene_name, .keep_all = TRUE), by = "gene_name")
  write_csv(TG |> arrange(Description, desc(lfc)) |> select(ID, Description, p.adjust, gene_name, gene_id, log2FC = lfc, padj), file.path(OUT, "GO_MF_up_terms_genes.csv"))
  # ---- 2. term x gene map ----
  gn <- TG |> count(gene_name, name = "n_terms") |> left_join(E |> distinct(gene_name, .keep_all = TRUE), by = "gene_name") |> arrange(desc(n_terms), desc(lfc)) |> slice_head(n = MAX_GENES)
  M <- TG |> filter(gene_name %in% gn$gene_name) |> mutate(gene_name = factor(gene_name, levels = rev(gn$gene_name)), Description = factor(Description, levels = sig$Description[order(sig$pvalue)]))
  p2 <- ggplot(M, aes(Description, gene_name, fill = lfc)) + geom_tile(colour = "white", linewidth = 0.4) +
    scale_fill_gradient2(low = "#2B6CB0", mid = "white", high = "#C53030", midpoint = 0, name = "log2FC\nEBOV vs Mock") + scale_x_discrete(labels = function(x) wrap(x, 28)) +
    labs(title = "Which up-regulated genes carry each molecular function", subtitle = sprintf("%s; a tile = the gene is in the term, colour = its EBOV log2FC; genes in most terms first (top %d)", expl, nrow(gn)), x = NULL, y = NULL) +
    theme_f() + theme(axis.text.x = element_text(angle = 50, hjust = 1, size = 8), axis.text.y = element_text(face = "italic", size = 7))
  save2(p2, "GO_MF_up_term_gene_map", 0.35 * n_distinct(M$Description) + 4.5, 0.16 * nrow(gn) + 4)
  # ---- 3. term network ----
  if (requireNamespace("igraph", quietly = TRUE) && nrow(sig) > 2) {
    gs <- split(TG$gene_name, TG$Description); nm <- names(gs)
    J <- outer(seq_along(gs), seq_along(gs), Vectorize(function(i, j) length(intersect(gs[[i]], gs[[j]])) / length(union(gs[[i]], gs[[j]]))))
    ed <- which(upper.tri(J) & J >= JACCARD, arr.ind = TRUE)
    g <- igraph::make_empty_graph(n = length(nm), directed = FALSE) |> igraph::set_vertex_attr("name", value = nm)
    if (nrow(ed)) g <- igraph::add_edges(g, as.vector(t(ed)), weight = J[ed])
    set.seed(1); lay <- igraph::layout_with_fr(g)
    V <- tibble(term = nm, x = lay[, 1], y = lay[, 2], n = lengths(gs)) |> left_join(sig |> select(Description, p.adjust), by = c("term" = "Description"))
    ES <- if (nrow(ed)) tibble(x = lay[ed[, 1], 1], y = lay[ed[, 1], 2], xend = lay[ed[, 2], 1], yend = lay[ed[, 2], 2], w = J[ed]) else tibble(x = numeric(), y = numeric(), xend = numeric(), yend = numeric(), w = numeric())
    p3 <- ggplot() + geom_segment(data = ES, aes(x = x, y = y, xend = xend, yend = yend, linewidth = w), colour = "grey70", alpha = 0.8) + scale_linewidth(range = c(0.3, 2.5), name = "shared genes\n(Jaccard)") +
      geom_point(data = V, aes(x, y, size = n, colour = p.adjust)) + scale_colour_gradient(low = "#C53030", high = "#9ECAE1", name = "p.adjust") + scale_size_continuous(range = c(3, 11), name = "genes") +
      ggrepel::geom_text_repel(data = V, aes(x, y, label = wrap(term, 26)), size = 2.8, colour = INK, lineheight = 0.9, max.overlaps = Inf, box.padding = 0.5, seed = 1) +
      labs(title = "Network of enriched molecular functions (EBOV up genes)", subtitle = sprintf("%s; terms linked when they share >= %.0f%% of their genes", expl, 100 * JACCARD)) + theme_void() +
      theme(plot.title = element_text(colour = INK, face = "bold"), plot.subtitle = element_text(colour = "grey35", size = 8.5), plot.margin = margin(10, 10, 10, 10))
    save2(p3, "GO_MF_up_term_network", 9, 7.5)
  } else message("igraph not installed (or < 3 terms) - network skipped (conda install -c conda-forge r-igraph)")
}

# ---- 4. KEGG signalling pathway maps ----
kg <- R |> filter(ontology == "KEGG") |> arrange(pvalue)
if (!nrow(kg)) message("No KEGG terms for the up set in the enrichment table - pathway maps skipped (KEGG needs internet in script 82)") else {
  kp <- kg |> filter(p.adjust < PADJ_SIG); lab <- "p.adjust < 0.05"
  if (!nrow(kp)) { kp <- kg |> slice_head(n = min(5, nrow(kg))); lab <- "exploratory (no pathway reaches p.adjust < 0.05)" }
  kp <- kp |> slice_head(n = N_PATH)
  message("KEGG pathways to draw (", lab, "):\n  ", paste(sprintf("%s %s", kp$ID, kp$Description), collapse = "\n  "))
  write_csv(kp |> select(ID, Description, GeneRatio, Count, pvalue, p.adjust, geneID), file.path(OUT, "KEGG_up_pathways_drawn.csv"))
  if (!requireNamespace("pathview", quietly = TRUE) || !requireNamespace("org.Mm.eg.db", quietly = TRUE)) {
    message("pathview not installed - KEGG maps skipped (conda install -c conda-forge -c bioconda bioconductor-pathview). The pathway list above is in KEGG_up_pathways_drawn.csv")
  } else {
    pcf <- file.path(RNA, "DE_v2_EBOV_vs_Mock_protein_coding.csv")
    d <- (if (file.exists(pcf)) read_csv(pcf, show_col_types = FALSE) else read_csv(file.path(RNA, "DE_v2_EBOV_vs_Mock_annotated.csv"), show_col_types = FALSE)) |> filter(!is.na(log2FoldChange), grepl("^ENSMUSG", gene_id)) |> mutate(ens = sub("\\..*$", "", gene_id))
    map <- suppressMessages(suppressWarnings(clusterProfiler::bitr(unique(d$ens), "ENSEMBL", "ENTREZID", org.Mm.eg.db::org.Mm.eg.db))) |> distinct(ENSEMBL, .keep_all = TRUE)
    gd <- d |> inner_join(map, by = c("ens" = "ENSEMBL")) |> group_by(ENTREZID) |> slice_max(abs(log2FoldChange), n = 1, with_ties = FALSE) |> ungroup()
    gene.data <- setNames(gd$log2FoldChange, gd$ENTREZID)
    pv <- file.path(OUT, "pathview"); dir.create(pv, showWarnings = FALSE); old <- setwd(pv); on.exit(setwd(old), add = TRUE)
    for (i in seq_len(nrow(kp))) {
      pid <- sub("^mmu", "", kp$ID[i]); message("   pathview ", kp$ID[i], " ", kp$Description[i])
      tryCatch(suppressMessages(pathview::pathview(gene.data = gene.data, pathway.id = pid, species = "mmu", kegg.dir = pv, limit = list(gene = 2, cpd = 1), bins = list(gene = 20, cpd = 10),
                                                   low = list(gene = "#2B6CB0", cpd = "blue"), mid = list(gene = "#F2F2F2", cpd = "gray"), high = list(gene = "#C53030", cpd = "red"), out.suffix = "EBOV_vs_Mock")),
               error = function(e) message("      failed: ", conditionMessage(e)))
    }
    setwd(old); message("KEGG maps written to ", pv, " (files mmuXXXXX.EBOV_vs_Mock.png; red = up, blue = down, colour range +/- 2 log2FC)")
  }
}
message("\nSaved to: ", OUT)
