#!/usr/bin/env Rscript
# STRING protein-protein interaction analysis of the protein-coding genes UP in (1) IFNB vs Mock and, separately, (2) EBOV vs Mock.
#   up = padj < PADJ_CUT and shrunk log2FC > LFC_CUT (default 1);  background = every protein-coding gene tested in that comparison (set in STRING for the enrichment).
#   Uses the STRINGdb package (mouse, taxon 10090); it downloads the STRING data for mouse the first time (a few hundred MB, cached in <project>/stringdb_cache) and needs internet.
# For each set it writes to requant_v2_outputs/figures_string/<IFNB|EBOV>/ (PNG + PDF):
#   cluster_names.csv          each cluster named after its most specific enriched STRING term (GO Process / Function, KEGG, Reactome, WikiPathways; FDR < 0.05, small background term size),
#                              with its top 3 terms and hub genes; clusters without a significant term are named after their hubs. The names are on the network and in the legend
#   network_clusters           interactions with combined score >= STRING_SCORE (default 700 = high confidence), isolated genes dropped; node colour = Louvain cluster, size = number of
#                              partners, the top N_LAB hubs labelled; clusters are numbered and listed in the legend with their size
#   network_log2FC             the same network coloured by the gene's log2FC
#   enrichment_top_terms       STRING functional enrichment of the set against the background (GO Process, Function, Component, KEGG, Reactome): top terms per category by FDR
#   interactions.csv, hubs.csv (degree, cluster, log2FC), clusters.csv (genes per cluster), enrichment_all.csv, ppi_stats.csv (nodes, edges, expected edges, PPI enrichment p-value)
# and shared_up_genes_IFNB_EBOV.csv (genes up in both) in the parent folder.
# Install if needed:  conda install -c conda-forge -c bioconda bioconductor-stringdb r-igraph
# Usage: Rscript 88_string_ifnb_ebov_up.R   (PROJ from 'source setup_session.sh';  STRING_SCORE=400  LFC_CUT=0  STRING_VERSION=11.5 optional)

t0 <- Sys.time(); say <- function(...) { message(sprintf("[%5.1f min] ", as.numeric(difftime(Sys.time(), t0, units = "mins"))), ...); flush.console() }
suppressPackageStartupMessages({ library(readr); library(dplyr); library(tidyr); library(ggplot2); library(ggrepel) })
for (pkg in c("STRINGdb", "igraph")) if (!requireNamespace(pkg, quietly = TRUE)) stop("Install ", pkg, " (conda install -c conda-forge -c bioconda ", if (pkg == "igraph") "r-igraph" else "bioconductor-stringdb", ")")
select <- dplyr::select; filter <- dplyr::filter; slice <- dplyr::slice; rename <- dplyr::rename
count <- dplyr::count; mutate <- dplyr::mutate; summarise <- dplyr::summarise; arrange <- dplyr::arrange; desc <- dplyr::desc   # STRINGdb loads plyr, whose versions of these verbs would otherwise take over

PADJ_CUT <- 0.05; LFC_CUT <- as.numeric(Sys.getenv("LFC_CUT", unset = "1")); STRING_SCORE <- as.integer(Sys.getenv("STRING_SCORE", unset = "700")); STRING_VERSION <- Sys.getenv("STRING_VERSION", unset = "12.0")
N_LAB <- 40; N_CLUST_COL <- 8; TOP_TERMS <- 12; MIN_CLUSTER <- 3
INK <- "#1F2933"
PROJ <- Sys.getenv("PROJ", unset = file.path(
  Sys.getenv("HOME"),
  "Documents/aviti_rnaseq/2025/aviti",
  "20251009-100725_Maury_25161_WivfFacICOTHFAUJwSPeHBEFfCfNtCRpUSgVGUkz",
  "results/Maury_25161"))
RNA <- file.path(PROJ, "requant_v2_outputs")
OUTR <- file.path(RNA, "figures_string"); dir.create(OUTR, recursive = TRUE, showWarnings = FALSE)
cache <- file.path(PROJ, "stringdb_cache"); dir.create(cache, showWarnings = FALSE)
wrap <- function(x, w = 48) vapply(strwrap(x, w, simplify = FALSE), paste, "", collapse = "\n")
options(timeout = 600)

say("connecting to STRING v", STRING_VERSION, " (mouse); the first run downloads the species data into ", cache)
sdb <- STRINGdb::STRINGdb$new(version = STRING_VERSION, species = 10090, score_threshold = 0L, network_type = "full", input_directory = cache)

CMP <- c(IFNB = "DE_v2_IFNB_vs_Mock", EBOV = "DE_v2_EBOV_vs_Mock")
load_de <- function(base) {
  f <- file.path(RNA, paste0(base, "_protein_coding.csv"))
  d <- if (file.exists(f)) read_csv(f, show_col_types = FALSE) else read_csv(file.path(RNA, paste0(base, "_annotated.csv")), show_col_types = FALSE) |> filter(gene_biotype == "protein_coding")
  d |> filter(!is.na(padj), grepl("^ENSMUSG", gene_id), !is.na(gene_name), gene_name != "") |> distinct(gene_name, .keep_all = TRUE)
}
ups <- list()
for (nm in names(CMP)) {
  d <- load_de(CMP[[nm]])
  g <- d |> filter(padj < PADJ_CUT, log2FoldChange > LFC_CUT)
  if (nrow(g) < 10 && LFC_CUT > 0) { say(nm, ": only ", nrow(g), " genes with log2FC > ", LFC_CUT, " - using padj alone"); g <- d |> filter(padj < PADJ_CUT, log2FoldChange > 0) }
  ups[[nm]] <- list(all = d, up = g); say(sprintf("%s: %d protein-coding genes tested, %d up", nm, nrow(d), nrow(g)))
}
sh <- inner_join(ups$IFNB$up |> select(gene_id, gene_name, lfc_IFNB = log2FoldChange, padj_IFNB = padj), ups$EBOV$up |> select(gene_id, lfc_EBOV = log2FoldChange, padj_EBOV = padj), by = "gene_id")
write_csv(sh, file.path(OUTR, "shared_up_genes_IFNB_EBOV.csv"))
say(sprintf("genes up in both IFNB and EBOV: %d (of %d IFNB-up and %d EBOV-up)", nrow(sh), nrow(ups$IFNB$up), nrow(ups$EBOV$up)))

for (nm in names(ups)) {
  OUT <- file.path(OUTR, nm); dir.create(OUT, showWarnings = FALSE)
  if (identical(Sys.getenv("SKIP_DONE"), "1") && file.exists(file.path(OUT, "enrichment_all.csv")) && file.exists(file.path(OUT, "network_clusters.png"))) { say("skipping ", nm, " (already done; SKIP_DONE=1)"); next }
  save2 <- function(p, name, w, h) { ggsave(file.path(OUT, paste0(name, ".png")), p, width = w, height = h, dpi = 300, limitsize = FALSE); ggsave(file.path(OUT, paste0(name, ".pdf")), p, width = w, height = h, limitsize = FALSE) }
  say("==== ", nm, " vs Mock: up genes")
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
    ttl <- sprintf("STRING network of genes up in %s vs Mock", nm); sub <- sprintf("%d genes, %d interactions (combined score >= %d); node size = number of partners; %d genes have no interaction at this score", nrow(V), nrow(ed), STRING_SCORE, nrow(hits) - nrow(V))
    p1 <- ggplot() + geom_segment(data = EDp, aes(x = x1, y = y1, xend = x2, yend = y2, alpha = combined_score), colour = "grey55", linewidth = 0.3) + scale_alpha_continuous(range = c(0.2, 0.7), guide = "none") +
      geom_point(data = V, aes(x, y, size = deg, colour = cluster)) + scale_colour_manual(values = ccol, name = NULL, labels = function(x) vapply(x, function(z) if (z %in% cn$cluster) cn$lbl[cn$cluster == z] else z, "")) + scale_size_continuous(range = c(1, 7), guide = "none") +
      geom_text_repel(data = lab, aes(x, y, label = gene_name), size = 2.6, fontface = "italic", colour = INK, max.overlaps = Inf, min.segment.length = 0, segment.size = 0.2, segment.colour = "grey60", box.padding = 0.25, seed = 1) +
      geom_label_repel(data = V |> filter(cluster != "other") |> group_by(cluster) |> summarise(x = median(x), y = median(y), .groups = "drop") |> left_join(cinfo |> select(cluster, short), by = "cluster") |> mutate(txt = wrap(short, 24)),
                       aes(x, y, label = txt, colour = cluster), size = 3, fontface = "bold", fill = alpha("white", 0.85), label.size = 0.25, lineheight = 0.9, max.overlaps = Inf, box.padding = 0.6, show.legend = FALSE, seed = 2) +
      labs(title = ttl, subtitle = sub) + theme_void(base_size = 11) + theme(plot.title = element_text(colour = INK, face = "bold"), plot.subtitle = element_text(colour = "grey35", size = 8.5), legend.position = "right", plot.margin = margin(10, 10, 10, 10))
    save2(p1, "network_clusters", 11, 8.5)
    lim <- max(abs(V$lfc), na.rm = TRUE)
    p2 <- ggplot() + geom_segment(data = EDp, aes(x = x1, y = y1, xend = x2, yend = y2, alpha = combined_score), colour = "grey55", linewidth = 0.3) + scale_alpha_continuous(range = c(0.2, 0.7), guide = "none") +
      geom_point(data = V, aes(x, y, size = deg, fill = lfc), shape = 21, colour = "grey30", stroke = 0.3) + scale_fill_gradient(low = "#FEE0D2", high = "#67000D", limits = c(min(V$lfc, na.rm = TRUE), lim), name = "log2FC") + scale_size_continuous(range = c(1, 7), guide = "none") +
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
        scale_x_continuous(expand = expansion(mult = c(0, 0.12))) + labs(title = sprintf("STRING functional enrichment: genes up in %s vs Mock", nm), subtitle = "top terms per category by FDR; number = genes in the term; background = all protein-coding genes tested", x = expression(-log[10](FDR)), y = NULL) +
        theme_bw(base_size = 10) + theme(panel.grid.minor = element_blank(), panel.border = element_rect(colour = "grey70"), axis.text = element_text(colour = INK), strip.background = element_rect(fill = "grey95", colour = "grey70"),
                                         strip.text = element_text(colour = INK, face = "bold"), plot.title = element_text(colour = INK, face = "bold"), plot.subtitle = element_text(colour = "grey35", size = 8.5), legend.position = "bottom")
      save2(p3, "enrichment_top_terms", 10, 0.2 * nrow(tp) + 2.8)
    } else say("no enrichment terms to draw for this set (see enrichment_all.csv)")
  }
}
say("done. Saved to: ", OUTR)
