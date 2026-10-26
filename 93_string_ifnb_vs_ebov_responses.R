#!/usr/bin/env Rscript
# STRING analysis of how the IFNB response and the EBOV response differ, gene class by gene class (protein-coding genes; IFNB vs Mock and EBOV vs Mock from the DE_v2_* tables):
#   Note: log2FC(EBOV vs Mock) - log2FC(IFNB vs Mock) is exactly the EBOV-vs-IFNB contrast, so the genes whose two responses DIFFER significantly are the sets of 91_string_ebov_vs_ifnb.R. This script instead
#   splits genes by how each response behaves, so shared and opposite behaviour can be read directly. Significant = padj < PADJ_CUT and |shrunk log2FC| > LFC_CUT (default 0.58 = 1.5-fold; "not significant"
#   in the other comparison = padj >= PADJ_CUT). Classes:
#     shared_up            up in IFNB and up in EBOV                      shared_down        down in both
#     IFNB_up_EBOV_down    up in IFNB, down in EBOV (opposite)            EBOV_up_IFNB_down  up in EBOV, down in IFNB (opposite)
#     IFNB_only_up         up in IFNB, not significant in EBOV            EBOV_only_up       up in EBOV, not significant in IFNB
#     IFNB_only_down       down in IFNB, not significant in EBOV          EBOV_only_down     down in EBOV, not significant in IFNB
#   Classes with fewer than MIN_SET genes (default 15) are listed but not analysed. Background = every protein-coding gene tested in both comparisons. Needs STRINGdb and internet (mouse data cached in <project>/stringdb_cache).
# For each analysed class it writes to requant_v2_outputs/figures_string_IFNB_vs_EBOV_responses/<class>/ (PNG + PDF): network_clusters (clusters named after their most specific enriched STRING term), network_log2FC,
# enrichment_top_terms, and cluster_names.csv, interactions.csv, hubs.csv, clusters.csv, enrichment_all.csv, ppi_stats.csv, genes_in_set.csv; plus class_sizes.csv and gene_classes.csv in the parent folder.
# Usage: Rscript 93_string_ifnb_vs_ebov_responses.R   (PROJ from 'source setup_session.sh';  LFC_CUT=0  MIN_SET=10  STRING_SCORE=400  STRING_VERSION=11.5 optional)

t0 <- Sys.time(); say <- function(...) { message(sprintf("[%5.1f min] ", as.numeric(difftime(Sys.time(), t0, units = "mins"))), ...); flush.console() }
suppressPackageStartupMessages({ library(readr); library(dplyr); library(tidyr); library(ggplot2); library(ggrepel) })
for (pkg in c("STRINGdb", "igraph")) if (!requireNamespace(pkg, quietly = TRUE)) stop("Install ", pkg, " (conda install -c conda-forge -c bioconda ", if (pkg == "igraph") "r-igraph" else "bioconductor-stringdb", ")")
select <- dplyr::select; filter <- dplyr::filter; slice <- dplyr::slice; rename <- dplyr::rename
count <- dplyr::count; mutate <- dplyr::mutate; summarise <- dplyr::summarise; arrange <- dplyr::arrange; desc <- dplyr::desc   # STRINGdb loads plyr, whose versions of these verbs would otherwise take over

PADJ_CUT <- 0.05; LFC_CUT <- as.numeric(Sys.getenv("LFC_CUT", unset = "0.58")); MIN_SET <- as.integer(Sys.getenv("MIN_SET", unset = "15")); STRING_SCORE <- as.integer(Sys.getenv("STRING_SCORE", unset = "700")); STRING_VERSION <- Sys.getenv("STRING_VERSION", unset = "12.0")
N_LAB <- 40; N_CLUST_COL <- 8; TOP_TERMS <- 12; MIN_CLUSTER <- 3
INK <- "#1F2933"
PROJ <- Sys.getenv("PROJ", unset = file.path(
  Sys.getenv("HOME"),
  "Documents/aviti_rnaseq/2025/aviti",
  "20251009-100725_Maury_25161_WivfFacICOTHFAUJwSPeHBEFfCfNtCRpUSgVGUkz",
  "results/Maury_25161"))
RNA <- file.path(PROJ, "requant_v2_outputs")
OUTR <- file.path(RNA, "figures_string_IFNB_vs_EBOV_responses"); dir.create(OUTR, recursive = TRUE, showWarnings = FALSE)
cache <- file.path(PROJ, "stringdb_cache"); dir.create(cache, showWarnings = FALSE)
wrap <- function(x, w = 48) vapply(strwrap(x, w, simplify = FALSE), paste, "", collapse = "\n")
options(timeout = 600)

say("connecting to STRING v", STRING_VERSION, " (mouse); the first run downloads the species data into ", cache)
sdb <- STRINGdb::STRINGdb$new(version = STRING_VERSION, species = 10090, score_threshold = 0L, network_type = "full", input_directory = cache)

rd <- function(f) read_csv(file.path(RNA, f), show_col_types = FALSE)
B <- rd("DE_v2_IFNB_vs_Mock_annotated.csv") |> filter(grepl("^ENSMUSG", gene_id), gene_biotype == "protein_coding", !is.na(gene_name), gene_name != "") |> transmute(gene_id, gene_name, baseMean, lfc_B = log2FoldChange, padj_B = padj)
E <- rd("DE_v2_EBOV_vs_Mock_annotated.csv") |> transmute(gene_id, lfc_E = log2FoldChange, padj_E = padj)
d <- B |> inner_join(E, by = "gene_id") |> filter(!is.na(padj_B), !is.na(padj_E)) |> distinct(gene_name, .keep_all = TRUE)
sgn <- function(l, p, s) !is.na(p) & p < PADJ_CUT & (if (s > 0) l > LFC_CUT else l < -LFC_CUT)
d <- d |> mutate(B_up = sgn(lfc_B, padj_B, 1), B_dn = sgn(lfc_B, padj_B, -1), E_up = sgn(lfc_E, padj_E, 1), E_dn = sgn(lfc_E, padj_E, -1), B_ns = padj_B >= PADJ_CUT, E_ns = padj_E >= PADJ_CUT,
                 class = case_when(B_up & E_up ~ "shared_up", B_dn & E_dn ~ "shared_down", B_up & E_dn ~ "IFNB_up_EBOV_down", E_up & B_dn ~ "EBOV_up_IFNB_down",
                                   B_up & E_ns ~ "IFNB_only_up", E_up & B_ns ~ "EBOV_only_up", B_dn & E_ns ~ "IFNB_only_down", E_dn & B_ns ~ "EBOV_only_down", TRUE ~ "other"))
CLS <- c("shared_up", "shared_down", "IFNB_up_EBOV_down", "EBOV_up_IFNB_down", "IFNB_only_up", "EBOV_only_up", "IFNB_only_down", "EBOV_only_down")
DESC <- c(shared_up = "up in both IFNB and EBOV", shared_down = "down in both IFNB and EBOV", IFNB_up_EBOV_down = "up in IFNB, down in EBOV", EBOV_up_IFNB_down = "up in EBOV, down in IFNB",
          IFNB_only_up = "up in IFNB only (not significant in EBOV)", EBOV_only_up = "up in EBOV only (not significant in IFNB)", IFNB_only_down = "down in IFNB only (not significant in EBOV)", EBOV_only_down = "down in EBOV only (not significant in IFNB)")
LEG <- c(shared_up = "mean log2FC", shared_down = "mean |log2FC|", IFNB_up_EBOV_down = "log2FC in IFNB", EBOV_up_IFNB_down = "log2FC in EBOV", IFNB_only_up = "log2FC in IFNB", EBOV_only_up = "log2FC in EBOV", IFNB_only_down = "|log2FC| in IFNB", EBOV_only_down = "|log2FC| in EBOV")
cs <- tibble(class = CLS, description = DESC[CLS], n_genes = vapply(CLS, function(k) sum(d$class == k), 1L), analysed = n_genes >= MIN_SET)
write_csv(cs, file.path(OUTR, "class_sizes.csv")); write_csv(d |> filter(class != "other") |> arrange(class) |> select(class, gene_id, gene_name, baseMean, lfc_B, padj_B, lfc_E, padj_E), file.path(OUTR, "gene_classes.csv"))
say(sprintf("%d protein-coding genes tested in both comparisons (|log2FC| > %g, padj < %g)", nrow(d), LFC_CUT, PADJ_CUT)); print(as.data.frame(cs), row.names = FALSE)
use <- cs$class[cs$analysed]; if (!length(use)) stop("No class has at least ", MIN_SET, " genes - try LFC_CUT=0 or MIN_SET=8")
ups <- setNames(lapply(use, function(k) { g <- d |> filter(class == k); lf <- switch(k, shared_up = (g$lfc_B + g$lfc_E) / 2, shared_down = -(g$lfc_B + g$lfc_E) / 2, IFNB_up_EBOV_down = g$lfc_B, EBOV_up_IFNB_down = g$lfc_E,
                                                  IFNB_only_up = g$lfc_B, EBOV_only_up = g$lfc_E, IFNB_only_down = -g$lfc_B, EBOV_only_down = -g$lfc_E)
                                  list(all = d, up = g |> mutate(log2FoldChange = lf, padj = pmin(padj_B, padj_E)) |> select(gene_id, gene_name, baseMean, log2FoldChange, padj)) }), use)
for (k in names(ups)) { dir.create(file.path(OUTR, k), showWarnings = FALSE); write_csv(ups[[k]]$up |> arrange(padj), file.path(OUTR, k, "genes_in_set.csv")) }
if (any(vapply(ups, function(u) nrow(u$up), 1L) < 3)) stop("A set has fewer than 3 genes - nothing to build a network from (try LFC_CUT=0)")

for (nm in names(ups)) {
  OUT <- file.path(OUTR, nm); dir.create(OUT, showWarnings = FALSE)
  if (identical(Sys.getenv("SKIP_DONE"), "1") && file.exists(file.path(OUT, "enrichment_all.csv")) && file.exists(file.path(OUT, "network_clusters.png"))) { say("skipping ", nm, " (already done; SKIP_DONE=1)"); next }
  save2 <- function(p, name, w, h) { ggsave(file.path(OUT, paste0(name, ".png")), p, width = w, height = h, dpi = 300, limitsize = FALSE); ggsave(file.path(OUT, paste0(name, ".pdf")), p, width = w, height = h, limitsize = FALSE) }
  say("==== genes ", DESC[[nm]])
  hits <- sdb$map(as.data.frame(ups[[nm]]$up |> select(gene_name, log2FoldChange)), "gene_name", removeUnmappedRows = TRUE) |> as_tibble() |> distinct(STRING_id, .keep_all = TRUE)
  say(sprintf("%d of %d genes mapped to STRING", nrow(hits), nrow(ups[[nm]]$up)))
  bg <- tryCatch(sdb$map(as.data.frame(ups[[nm]]$all |> select(gene_name)), "gene_name", removeUnmappedRows = TRUE) |> as_tibble() |> distinct(STRING_id), error = function(e) NULL)
  if (!is.null(bg)) { sdb$set_background(bg$STRING_id); say(sprintf("background set: %d proteins", nrow(bg))) }

  # --- PPI statistics ---
  ppi <- tryCatch(sdb$get_ppi_enrichment(hits$STRING_id), error = function(e) { say("PPI enrichment failed: ", conditionMessage(e)); NULL })
  if (!is.null(ppi)) {
    ppi_df <- tryCatch(as.data.frame(lapply(as.list(ppi), function(x) if (length(x)) as.character(x[[1]]) else NA_character_), stringsAsFactors = FALSE), error = function(e) NULL)   # elements can be empty
    if (!is.null(ppi_df) && nrow(ppi_df)) { write_csv(ppi_df, file.path(OUT, "ppi_stats.csv")); say("PPI: ", paste(names(ppi_df), unlist(ppi_df[1, ]), sep = " = ", collapse = "; ")) } else say("PPI enrichment returned no usable values - skipped")
  }

  # --- network ---
  ed <- sdb$get_interactions(hits$STRING_id) |> as_tibble() |> filter(combined_score >= STRING_SCORE, from != to) |> mutate(a = pmin(from, to), b = pmax(from, to)) |> distinct(a, b, .keep_all = TRUE) |> select(from = a, to = b, combined_score)
  nmap <- hits |> select(STRING_id, gene_name, lfc = log2FoldChange)
  write_csv(ed |> left_join(nmap |> rename(from = STRING_id, gene_from = gene_name) |> select(from, gene_from), by = "from") |> left_join(nmap |> rename(to = STRING_id, gene_to = gene_name) |> select(to, gene_to), by = "to") |> arrange(desc(combined_score)),
            file.path(OUT, "interactions.csv"))
  say(sprintf("%d interactions with combined score >= %d", nrow(ed), STRING_SCORE))
  if (nrow(ed) >= 3) {
    g <- igraph::graph_from_data_frame(ed, directed = FALSE, vertices = nmap |> filter(STRING_id %in% c(ed$from, ed$to)))
    igraph::E(g)$weight <- ed$combined_score / 1000
    cl <- igraph::cluster_louvain(g); mem <- igraph::membership(cl)
    sz <- sort(table(mem), decreasing = TRUE); keep <- as.integer(names(sz)[sz >= MIN_CLUSTER]); top <- head(keep, N_CLUST_COL)
    clab <- setNames(seq_along(top), top)
    set.seed(1); lay <- igraph::layout_with_fr(g, weights = igraph::E(g)$weight)
    V <- tibble(STRING_id = igraph::V(g)$name, gene_name = igraph::V(g)$gene_name, lfc = igraph::V(g)$lfc, x = lay[, 1], y = lay[, 2], deg = igraph::degree(g), cl_raw = as.integer(mem[igraph::V(g)$name])) |>
      mutate(cluster = ifelse(cl_raw %in% top, paste0("cluster ", clab[as.character(cl_raw)]), "other"))
    # --- name every cluster after its most specific enriched term (STRING enrichment of that cluster's genes against the background); fall back to its hub genes ---
    name_cluster <- function(ids) {
      en <- tryCatch(as_tibble(sdb$get_enrichment(ids)), error = function(e) NULL)
      if (is.null(en) || !nrow(en)) return(tibble(name = NA_character_, top = NA_character_))
      en <- en |> mutate(fdr = as.numeric(fdr), nbg = as.numeric(number_of_genes_in_background), ng = as.numeric(number_of_genes)) |>
        filter(category %in% c("Process", "KEGG", "RCTM", "Function", "WikiPathways"), fdr < 0.05) |> arrange(fdr, desc(ng))
      sp <- en |> filter(nbg <= 600); if (!nrow(sp)) sp <- en
      if (!nrow(sp)) return(tibble(name = NA_character_, top = NA_character_))
      tibble(name = sp$description[1], top = paste(head(unique(sp$description), 3), collapse = " | "))
    }
    cinfo <- bind_rows(lapply(sort(setdiff(unique(V$cluster), "other")), function(cl) {
      ids <- V$STRING_id[V$cluster == cl]; r <- name_cluster(ids)
      hubs3 <- paste(head(V$gene_name[V$cluster == cl][order(-V$deg[V$cluster == cl])], 3), collapse = ", ")
      tibble(cluster = cl, n_genes = length(ids), cluster_name = if (is.na(r$name)) paste0("no enriched function (hubs: ", hubs3, ")") else r$name, top_terms = r$top, hubs = hubs3)
    })) |> mutate(short = ifelse(nchar(cluster_name) > 52, paste0(substr(cluster_name, 1, 49), "..."), cluster_name))
    for (i in seq_len(nrow(cinfo))) say(sprintf("   %s: %s  (%d genes; hubs %s)", cinfo$cluster[i], cinfo$cluster_name[i], cinfo$n_genes[i], cinfo$hubs[i]))
    write_csv(cinfo |> select(cluster, cluster_name, n_genes, top_terms, hubs), file.path(OUT, "cluster_names.csv"))
    cn <- V |> count(cluster) |> left_join(cinfo |> select(cluster, short), by = "cluster") |> mutate(lbl = ifelse(is.na(short), sprintf("%s (%d genes)", cluster, n), sprintf("%s. %s (%d genes)", sub("cluster ", "", cluster), short, n)))
    V <- V |> left_join(cinfo |> select(cluster, cluster_name), by = "cluster")
    V$cluster <- factor(V$cluster, levels = c(sort(setdiff(unique(V$cluster), "other")), "other"))
    write_csv(V |> select(gene_name, STRING_id, degree = deg, cluster, cluster_name, log2FC = lfc) |> arrange(cluster, desc(degree)), file.path(OUT, "hubs.csv"))
    write_csv(V |> filter(cluster != "other") |> group_by(cluster, cluster_name) |> summarise(n_genes = n(), genes = paste(gene_name[order(-deg)], collapse = ", "), .groups = "drop"), file.path(OUT, "clusters.csv"))
    EDp <- ed |> left_join(V |> select(STRING_id, x, y), by = c("from" = "STRING_id")) |> rename(x1 = x, y1 = y) |> left_join(V |> select(STRING_id, x, y), by = c("to" = "STRING_id")) |> rename(x2 = x, y2 = y)
    ccol <- setNames(c(scales::hue_pal(l = 55)(length(levels(V$cluster)) - 1), "grey70")[seq_along(levels(V$cluster))], levels(V$cluster))
    lab <- V |> slice_max(deg, n = N_LAB, with_ties = FALSE)
    ttl <- sprintf("STRING network of genes %s", DESC[[nm]]); sub <- sprintf("%d genes, %d interactions (combined score >= %d); node size = number of partners; %d genes have no interaction at this score", nrow(V), nrow(ed), STRING_SCORE, nrow(hits) - nrow(V))
    p1 <- ggplot() + geom_segment(data = EDp, aes(x = x1, y = y1, xend = x2, yend = y2, alpha = combined_score), colour = "grey55", linewidth = 0.3) + scale_alpha_continuous(range = c(0.2, 0.7), guide = "none") +
      geom_point(data = V, aes(x, y, size = deg, colour = cluster)) + scale_colour_manual(values = ccol, name = NULL, labels = function(x) vapply(x, function(z) if (z %in% cn$cluster) cn$lbl[cn$cluster == z] else z, "")) + scale_size_continuous(range = c(1, 7), guide = "none") +
      geom_text_repel(data = lab, aes(x, y, label = gene_name), size = 2.6, fontface = "italic", colour = INK, max.overlaps = Inf, min.segment.length = 0, segment.size = 0.2, segment.colour = "grey60", box.padding = 0.25, seed = 1) +
      geom_label_repel(data = V |> filter(cluster != "other") |> group_by(cluster) |> summarise(x = median(x), y = median(y), .groups = "drop") |> left_join(cinfo |> select(cluster, short), by = "cluster") |> mutate(txt = wrap(short, 24)),
                       aes(x, y, label = txt, colour = cluster), size = 3, fontface = "bold", fill = alpha("white", 0.85), label.size = 0.25, lineheight = 0.9, max.overlaps = Inf, box.padding = 0.6, show.legend = FALSE, seed = 2) +
      labs(title = ttl, subtitle = sub) + theme_void(base_size = 11) + theme(plot.title = element_text(colour = INK, face = "bold"), plot.subtitle = element_text(colour = "grey35", size = 8.5), legend.position = "right", plot.margin = margin(10, 10, 10, 10))
    save2(p1, "network_clusters", 11, 8.5)
    lim <- max(abs(V$lfc), na.rm = TRUE)
    p2 <- ggplot() + geom_segment(data = EDp, aes(x = x1, y = y1, xend = x2, yend = y2, alpha = combined_score), colour = "grey55", linewidth = 0.3) + scale_alpha_continuous(range = c(0.2, 0.7), guide = "none") +
      geom_point(data = V, aes(x, y, size = deg, fill = lfc), shape = 21, colour = "grey30", stroke = 0.3) + scale_fill_gradient(low = "#FEE0D2", high = "#67000D", limits = c(min(V$lfc, na.rm = TRUE), lim), name = gsub(" ", "\n", LEG[[nm]], fixed = TRUE)) + scale_size_continuous(range = c(1, 7), guide = "none") +
      geom_text_repel(data = lab, aes(x, y, label = gene_name), size = 2.6, fontface = "italic", colour = INK, max.overlaps = Inf, min.segment.length = 0, segment.size = 0.2, segment.colour = "grey60", box.padding = 0.25, seed = 1) +
      labs(title = ttl, subtitle = paste0(sub, "; colour = log2FC")) + theme_void(base_size = 11) + theme(plot.title = element_text(colour = INK, face = "bold"), plot.subtitle = element_text(colour = "grey35", size = 8.5), legend.position = "right", plot.margin = margin(10, 10, 10, 10))
    save2(p2, "network_log2FC", 11, 8.5)
  } else say("fewer than 3 interactions at this score - network figures skipped (try STRING_SCORE=400)")

  # --- functional enrichment ---
  en <- tryCatch(sdb$get_enrichment(hits$STRING_id) |> as_tibble(), error = function(e) { say("enrichment failed: ", conditionMessage(e)); NULL })
  if (!is.null(en) && nrow(en)) {
    en <- en |> mutate(fdr = as.numeric(fdr), p_value = as.numeric(p_value)) |> arrange(fdr)
    write_csv(en |> select(any_of(c("category", "term", "description", "number_of_genes", "number_of_genes_in_background", "p_value", "fdr", "preferredNames"))), file.path(OUT, "enrichment_all.csv"))
    say("enrichment categories returned: ", paste(sprintf("%s (%d)", names(table(en$category)), as.integer(table(en$category))), collapse = ", "))
    cats <- intersect(c("Process", "Function", "Component", "KEGG", "RCTM", "WikiPathways"), unique(en$category))
    if (!length(cats)) cats <- head(names(sort(table(en$category), decreasing = TRUE)), 5)      # fall back to whatever categories STRING returned
    say("enrichment: ", paste(sprintf("%s %d (FDR < 0.05: %d)", cats, vapply(cats, function(k) sum(en$category == k), 1L), vapply(cats, function(k) sum(en$category == k & en$fdr < 0.05), 1L)), collapse = " | "))
    tp <- en |> filter(category %in% cats) |> group_by(category) |> slice_min(fdr, n = TOP_TERMS, with_ties = FALSE) |> ungroup() |>
      mutate(category = factor(category, levels = cats), y = factor(paste(category, term), levels = rev(paste(category, term))), sigf = ifelse(fdr < 0.05, "FDR < 0.05", "FDR >= 0.05"))
    if (nrow(tp)) {
      nmv <- setNames(tp$description, as.character(tp$y))
      p3 <- ggplot(tp, aes(-log10(fdr), y, fill = sigf)) + geom_col(width = 0.75) + geom_text(aes(label = number_of_genes), hjust = -0.3, size = 2.8, colour = INK) +
        facet_grid(category ~ ., scales = "free_y", space = "free_y") + scale_y_discrete(labels = function(x) wrap(nmv[x])) + scale_fill_manual(values = c(`FDR < 0.05` = "#08519C", `FDR >= 0.05` = "#C6C6C6"), name = NULL) +
        scale_x_continuous(expand = expansion(mult = c(0, 0.12))) + labs(title = sprintf("STRING functional enrichment: genes %s", DESC[[nm]]), subtitle = "top terms per category by FDR; number = genes in the term; background = all protein-coding genes tested", x = expression(-log[10](FDR)), y = NULL) +
        theme_bw(base_size = 10) + theme(panel.grid.minor = element_blank(), panel.border = element_rect(colour = "grey70"), axis.text = element_text(colour = INK), strip.background = element_rect(fill = "grey95", colour = "grey70"),
                                         strip.text = element_text(colour = INK, face = "bold"), plot.title = element_text(colour = INK, face = "bold"), plot.subtitle = element_text(colour = "grey35", size = 8.5), legend.position = "bottom")
      save2(p3, "enrichment_top_terms", 10, 0.2 * nrow(tp) + 2.8)
    } else say("no enrichment terms to draw for this set (see enrichment_all.csv)")
  }
}
say("done. Saved to: ", OUTR)
