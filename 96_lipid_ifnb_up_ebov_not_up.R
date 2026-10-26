#!/usr/bin/env Rscript
# Lipid genes that IFNB induces but EBOV does not: of all host genes for lipid synthesis, metabolism, modification and transport, which are (1) UP in IFNB vs Mock and (2) DOWN or UNCHANGED (not
# significantly up) in EBOV vs Mock. "Lipid gene" = annotated in GO (org.Mm.eg.db, children included) to lipid metabolic process (GO:0006629), lipid biosynthetic process (GO:0008610), lipid catabolic
# process (GO:0016042), lipid transport (GO:0006869), lipid localization (GO:0010876), lipid modification (GO:0030258), lipid translocation (GO:0034204) or flippase / floppase / scramblase / intramembrane
# lipid transporter activity (GO:0140327, GO:0140328, GO:0017128, GO:0140326). Protein-coding genes only. IFNB-up = padj < PADJ_CUT and log2FC > LFC_UP (default 0.58 = 1.5-fold; LFC_UP=0 for any increase);
# EBOV down = padj < PADJ_CUT and log2FC < 0; EBOV unchanged = padj >= PADJ_CUT.
# Writes to requant_v2_outputs/figures_lipid_ifnb_up_ebov_not_up/ (PNG + PDF unless noted):
#   scatter_lipid_IFNB_vs_EBOV_log2FC        every lipid gene, x = log2FC IFNB, y = log2FC EBOV; hits coloured by EBOV status (down / unchanged), IFNB-induced genes that EBOV also raises in grey-red
#   bars_hits_p<n>                           vertical bars for every hit, left to right IFNB vs Mock, EBOV vs Mock, IFNB+EBOV vs Mock (shrunk log2FC +/- 1 SE, stars = padj), grouped by lipid function, paged
#   lipid_IFNB_up_EBOV_not_up.csv            the hit list with GO function groups, log2FC / padj in the three comparisons, and the EBOV status;  lipid_genes_all.csv = every lipid gene tested
#   IFNB_up_EBOV_not_up/                     STRING network (clusters named after their most specific enriched term), log2FC network, functional enrichment, tables (needs STRINGdb + internet; the figures
#                                            above are written first, so they exist even if STRING fails)
# Usage: Rscript 96_lipid_ifnb_up_ebov_not_up.R   (PROJ from 'source setup_session.sh';  LFC_UP=0  STRING_SCORE=400  STRING_VERSION=11.5  SKIP_STRING=1 optional)

t0 <- Sys.time(); say <- function(...) { message(sprintf("[%5.1f min] ", as.numeric(difftime(Sys.time(), t0, units = "mins"))), ...); flush.console() }
suppressPackageStartupMessages({ library(readr); library(dplyr); library(tidyr); library(ggplot2); library(ggrepel) })
for (pkg in c("clusterProfiler", "org.Mm.eg.db")) if (!requireNamespace(pkg, quietly = TRUE)) stop("Install ", pkg, " (conda install -c conda-forge -c bioconda bioconductor-", tolower(pkg), ")")
HAVE_STRING <- requireNamespace("STRINGdb", quietly = TRUE) && requireNamespace("igraph", quietly = TRUE) && !identical(Sys.getenv("SKIP_STRING"), "1")
if (!HAVE_STRING) message("STRING part skipped (STRINGdb / igraph not installed, or SKIP_STRING=1) - the lipid figures are still made")
select <- dplyr::select; filter <- dplyr::filter; slice <- dplyr::slice; rename <- dplyr::rename
count <- dplyr::count; mutate <- dplyr::mutate; summarise <- dplyr::summarise; arrange <- dplyr::arrange; desc <- dplyr::desc   # STRINGdb loads plyr, whose versions of these verbs would otherwise take over

PADJ_CUT <- 0.05; LFC_UP <- as.numeric(Sys.getenv("LFC_UP", unset = "0.58")); PER_PAGE <- 40; LOW_EXPR <- 50; STRING_SCORE <- as.integer(Sys.getenv("STRING_SCORE", unset = "700")); STRING_VERSION <- Sys.getenv("STRING_VERSION", unset = "12.0")
N_LAB <- 40; N_CLUST_COL <- 8; TOP_TERMS <- 12; MIN_CLUSTER <- 3
INK <- "#1F2933"
PROJ <- Sys.getenv("PROJ", unset = file.path(
  Sys.getenv("HOME"),
  "Documents/aviti_rnaseq/2025/aviti",
  "20251009-100725_Maury_25161_WivfFacICOTHFAUJwSPeHBEFfCfNtCRpUSgVGUkz",
  "results/Maury_25161"))
RNA <- file.path(PROJ, "requant_v2_outputs")
OUTR <- file.path(RNA, "figures_lipid_ifnb_up_ebov_not_up"); dir.create(OUTR, recursive = TRUE, showWarnings = FALSE)
cache <- file.path(PROJ, "stringdb_cache"); dir.create(cache, showWarnings = FALSE)
wrap <- function(x, w = 48) vapply(strwrap(x, w, simplify = FALSE), paste, "", collapse = "\n")
options(timeout = 600)

save2 <- function(p, name, w, h) { ggsave(file.path(OUTR, paste0(name, ".png")), p, width = w, height = h, dpi = 300, limitsize = FALSE); ggsave(file.path(OUTR, paste0(name, ".pdf")), p, width = w, height = h, limitsize = FALSE) }
stars <- function(p) ifelse(is.na(p), "", ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", "ns"))))
rd <- function(f) read_csv(file.path(RNA, f), show_col_types = FALSE)

# ---- the lipid gene universe from GO ----
GO <- list(synthesis = "GO:0008610", catabolism = "GO:0016042", modification = "GO:0030258",
           transport = c("GO:0006869", "GO:0010876", "GO:0034204", "GO:0140327", "GO:0140328", "GO:0017128", "GO:0140326"), metabolism = "GO:0006629")
g2a <- as.list(org.Mm.eg.db::org.Mm.egGO2ALLEGS)
ens_of <- function(go) { e <- unique(as.character(unlist(g2a[go]))); e <- e[!is.na(e)]; if (!length(e)) return(character()); unique(suppressMessages(suppressWarnings(clusterProfiler::bitr(e, "ENTREZID", "ENSEMBL", org.Mm.eg.db::org.Mm.eg.db)))$ENSEMBL) }
mem <- lapply(GO, ens_of); say(sprintf("GO lipid sets (genes): %s", paste(sprintf("%s %d", names(mem), lengths(mem)), collapse = " | ")))
lipid_ens <- unique(unlist(mem))
FUN_ORDER <- c("transport / translocation", "modification", "synthesis", "catabolism", "other lipid metabolism")
fn_of <- function(e) ifelse(e %in% mem$transport, FUN_ORDER[1], ifelse(e %in% mem$modification, FUN_ORDER[2], ifelse(e %in% mem$synthesis, FUN_ORDER[3], ifelse(e %in% mem$catabolism, FUN_ORDER[4], FUN_ORDER[5]))))
fn_all <- function(e) paste(c(if (e %in% mem$transport) "transport", if (e %in% mem$modification) "modification", if (e %in% mem$synthesis) "synthesis", if (e %in% mem$catabolism) "catabolism", if (e %in% mem$metabolism) "metabolism"), collapse = "; ")

# ---- the three comparisons ----
B <- rd("DE_v2_IFNB_vs_Mock_annotated.csv") |> filter(grepl("^ENSMUSG", gene_id), gene_biotype == "protein_coding", !is.na(gene_name), gene_name != "") |> transmute(gene_id, ens = sub("\\..*$", "", gene_id), gene_name, baseMean, lfc_B = log2FoldChange, se_B = lfcSE, padj_B = padj)
E <- rd("DE_v2_EBOV_vs_Mock_annotated.csv") |> transmute(gene_id, lfc_E = log2FoldChange, se_E = lfcSE, padj_E = padj)
BE <- rd("DE_v2_IFNB_EBOV_vs_Mock_annotated.csv") |> transmute(gene_id, lfc_BE = log2FoldChange, se_BE = lfcSE, padj_BE = padj)
d <- B |> inner_join(E, by = "gene_id") |> left_join(BE, by = "gene_id") |> filter(!is.na(padj_B), !is.na(padj_E)) |> distinct(gene_name, .keep_all = TRUE)
d <- d |> mutate(lipid = ens %in% lipid_ens, ifnb_up = padj_B < PADJ_CUT & lfc_B > LFC_UP,
                 ebov = case_when(padj_E < PADJ_CUT & lfc_E > 0 ~ "up", padj_E < PADJ_CUT & lfc_E < 0 ~ "down", TRUE ~ "unchanged"), hit = lipid & ifnb_up & ebov != "up")
L <- d |> filter(lipid); H <- d |> filter(hit)
say(sprintf("%d protein-coding genes tested; %d lipid genes; %d of them up in IFNB (log2FC > %g); %d also up in EBOV; %d HITS = IFNB-up and EBOV not up (%d EBOV down, %d EBOV unchanged)",
            nrow(d), nrow(L), sum(L$ifnb_up), sum(L$ifnb_up & L$ebov == "up"), nrow(H), sum(H$ebov == "down"), sum(H$ebov == "unchanged")))
write_csv(L |> mutate(go_functions = vapply(ens, fn_all, ""), group = fn_of(ens)) |> arrange(desc(lfc_B)) |> select(gene_name, gene_id, baseMean, group, go_functions, lfc_IFNB = lfc_B, padj_IFNB = padj_B, lfc_EBOV = lfc_E, padj_EBOV = padj_E, lfc_IFNB_EBOV = lfc_BE, padj_IFNB_EBOV = padj_BE), file.path(OUTR, "lipid_genes_all.csv"))
if (!nrow(H)) stop("No lipid gene is up in IFNB and not up in EBOV at these cut-offs (try LFC_UP=0)")
write_csv(H |> mutate(go_functions = vapply(ens, fn_all, ""), group = fn_of(ens), ebov_status = ebov) |> arrange(group, desc(lfc_B)) |> select(gene_name, gene_id, baseMean, group, go_functions, ebov_status, lfc_IFNB = lfc_B, padj_IFNB = padj_B, lfc_EBOV = lfc_E, padj_EBOV = padj_E, lfc_IFNB_EBOV = lfc_BE, padj_IFNB_EBOV = padj_BE),
          file.path(OUTR, "lipid_IFNB_up_EBOV_not_up.csv"))

# ---- scatter ----
CS <- c(`EBOV down` = "#2B6CB0", `EBOV unchanged` = "#E08214", `also up in EBOV` = "#C5A0A0", `other lipid genes` = "grey75")
L <- L |> mutate(cls = factor(case_when(hit & ebov == "down" ~ "EBOV down", hit ~ "EBOV unchanged", ifnb_up & ebov == "up" ~ "also up in EBOV", TRUE ~ "other lipid genes"), levels = names(CS)))
lim <- as.numeric(quantile(c(abs(L$lfc_B), abs(L$lfc_E)), 0.995)); lab <- H |> arrange(desc(lfc_B)) |> slice_head(n = 30)
ps <- ggplot(L, aes(lfc_B, lfc_E)) + geom_point(data = L |> filter(cls == "other lipid genes"), colour = CS[["other lipid genes"]], size = 0.9, alpha = 0.6) +
  geom_point(data = L |> filter(cls == "also up in EBOV"), colour = CS[["also up in EBOV"]], size = 1.5) + geom_point(data = L |> filter(hit), aes(colour = cls), size = 2) +
  geom_text_repel(data = lab, aes(label = gene_name, colour = ifelse(ebov == "down", "EBOV down", "EBOV unchanged")), size = 3, fontface = "italic", max.overlaps = Inf, min.segment.length = 0, segment.size = 0.2, segment.colour = "grey60", box.padding = 0.3, seed = 1, show.legend = FALSE) +
  geom_hline(yintercept = 0, colour = INK, linewidth = 0.35) + geom_vline(xintercept = c(0, LFC_UP), colour = c(INK, "grey55"), linetype = c("solid", "dashed"), linewidth = 0.35) +
  scale_colour_manual(values = CS, breaks = names(CS), name = NULL) + coord_cartesian(xlim = c(-lim, lim), ylim = c(-lim, lim)) +
  labs(title = "Lipid genes: induced by IFNB, but not by EBOV", subtitle = sprintf("%d lipid genes tested; coloured = up in IFNB (padj < %g, log2FC > %g) and down (blue) or not significantly changed (orange) in EBOV; pink = IFNB-up and also EBOV-up", nrow(L), PADJ_CUT, LFC_UP),
       x = "log2FC IFNB vs Mock", y = "log2FC EBOV vs Mock") + theme_bw(base_size = 12) + theme(panel.grid.minor = element_blank(), panel.grid.major = element_line(colour = "grey94", linewidth = 0.3), panel.border = element_rect(colour = "grey70"),
       axis.text = element_text(colour = INK), plot.title = element_text(colour = INK, face = "bold"), plot.subtitle = element_text(colour = "grey35", size = 8.5), legend.position = "bottom")
save2(ps, "scatter_lipid_IFNB_vs_EBOV_log2FC", 8, 7.5)

# ---- vertical bars, three comparisons, grouped by lipid function ----
CMPN <- c(`IFNB vs Mock` = "#2B6CB0", `EBOV vs Mock` = "#1B7837", `IFNB+EBOV vs Mock` = "#9BD08E")
HL <- H |> mutate(group = factor(fn_of(ens), levels = FUN_ORDER)) |> arrange(group, desc(lfc_B)) |> mutate(idx = row_number(), page = ceiling(idx / PER_PAGE), xl = paste0(gene_name, ifelse(baseMean < LOW_EXPR, "†", "")))
long <- bind_rows(HL |> transmute(gene_name, xl, group, page, idx, comparison = "IFNB vs Mock", lfc = lfc_B, se = se_B, padj = padj_B), HL |> transmute(gene_name, xl, group, page, idx, comparison = "EBOV vs Mock", lfc = lfc_E, se = se_E, padj = padj_E),
                  HL |> transmute(gene_name, xl, group, page, idx, comparison = "IFNB+EBOV vs Mock", lfc = lfc_BE, se = se_BE, padj = padj_BE)) |>
  mutate(comparison = factor(comparison, levels = names(CMPN)), lab = stars(padj), vj = ifelse(lfc >= 0, 0, 1), ys = lfc + ifelse(lfc >= 0, 1, -1) * (se + 0.03), key = paste(group, xl, sep = "||"))
dg <- position_dodge(width = 0.85); np <- max(HL$page)
for (pg in seq_len(np)) {
  dd <- long |> filter(page == pg) |> mutate(group = droplevels(group)); dd$key <- factor(dd$key, levels = unique(dd$key[order(dd$idx)]))
  p <- ggplot(dd, aes(key, lfc, fill = comparison)) + geom_col(position = dg, width = 0.8) + geom_errorbar(aes(ymin = lfc - se, ymax = lfc + se), position = dg, width = 0.25, colour = INK, linewidth = 0.3) +
    geom_text(aes(y = ys, label = lab, vjust = vj), position = dg, colour = INK, size = 2.1) + geom_hline(yintercept = 0, colour = INK, linewidth = 0.4) +
    facet_grid(. ~ group, scales = "free_x", space = "free_x", labeller = labeller(group = label_wrap_gen(16))) + scale_fill_manual(values = CMPN, name = NULL) +
    scale_x_discrete(labels = function(x) sub("^.*\\|\\|", "", x)) + scale_y_continuous(expand = expansion(mult = 0.12)) +
    labs(title = "Lipid genes up in IFNB and not up in EBOV", subtitle = sprintf("%d genes%s; left to right: IFNB, EBOV, IFNB+EBOV vs Mock; shrunk log2FC +/- 1 SE; stars = padj of each comparison; dagger = baseMean < %d", nrow(HL), if (np > 1) sprintf(", page %d of %d", pg, np) else "", LOW_EXPR),
         x = NULL, y = "log2 fold-change vs Mock") + theme_bw(base_size = 10) + theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(), panel.grid.major.y = element_line(colour = "grey92", linewidth = 0.3), panel.border = element_rect(colour = "grey70"),
          axis.text = element_text(colour = INK), axis.text.x = element_text(face = "italic", angle = 50, hjust = 1, size = 8), plot.title = element_text(colour = INK, face = "bold", size = 11), plot.subtitle = element_text(colour = "grey35", size = 8),
          strip.background = element_rect(fill = "grey95", colour = "grey70"), strip.text = element_text(colour = INK, face = "bold", size = 8), legend.position = "top", panel.spacing.x = unit(0.15, "lines"))
  save2(p, if (np > 1) sprintf("bars_hits_p%d", pg) else "bars_hits", max(6, 0.5 * n_distinct(dd$key) + 1.4 * n_distinct(dd$group) + 2), 5.8)
}
say(sprintf("figures written: scatter + %d bar page(s)", np))

# ---- STRING on the hit list ----
DESC <- c(IFNB_up_EBOV_not_up = "lipid genes up in IFNB and not up in EBOV"); LEG <- c(IFNB_up_EBOV_not_up = "log2FC in IFNB")
ups <- list(IFNB_up_EBOV_not_up = list(all = d, up = H |> transmute(gene_id, gene_name, baseMean, log2FoldChange = lfc_B, padj = padj_B)))
if (HAVE_STRING) {
say("connecting to STRING v", STRING_VERSION, " (mouse); the first run downloads the species data into ", cache)
sdb <- STRINGdb::STRINGdb$new(version = STRING_VERSION, species = 10090, score_threshold = 0L, network_type = "full", input_directory = cache)
for (k in names(ups)) { dir.create(file.path(OUTR, k), showWarnings = FALSE); write_csv(ups[[k]]$up |> arrange(padj), file.path(OUTR, k, "genes_in_set.csv")) }
if (any(vapply(ups, function(u) nrow(u$up), 1L) < 3)) stop("A set has fewer than 3 genes - nothing to build a network from (try LFC_UP=0)")

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
  ed_all <- sdb$get_interactions(hits$STRING_id) |> as_tibble() |> filter(from != to) |> mutate(a = pmin(from, to), b = pmax(from, to)) |> distinct(a, b, .keep_all = TRUE) |> select(from = a, to = b, combined_score)
  SC <- STRING_SCORE; ed <- ed_all |> filter(combined_score >= SC)
  for (lower in c(400L, 150L)) if (nrow(ed) < 3 && SC > lower) { say(sprintf("only %d interactions at score >= %d - lowering the threshold to %d for this set", nrow(ed), SC, lower)); SC <- lower; ed <- ed_all |> filter(combined_score >= SC) }
  nmap <- hits |> select(STRING_id, gene_name, lfc = log2FoldChange)
  write_csv(ed |> left_join(nmap |> rename(from = STRING_id, gene_from = gene_name) |> select(from, gene_from), by = "from") |> left_join(nmap |> rename(to = STRING_id, gene_to = gene_name) |> select(to, gene_to), by = "to") |> arrange(desc(combined_score)),
            file.path(OUT, "interactions.csv"))
  say(sprintf("%d interactions with combined score >= %d", nrow(ed), SC))
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
    ttl <- sprintf("STRING network of genes %s", DESC[[nm]]); sub <- sprintf("%d genes, %d interactions (combined score >= %d); node size = number of partners; %d genes have no interaction at this score", nrow(V), nrow(ed), SC, nrow(hits) - nrow(V))
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
}
