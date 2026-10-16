#!/usr/bin/env Rscript
# What do the up- and downregulated lipid genes (script 27) do?
#
# For each input table (all lipid genes; endocytic-vesicle lipid enzymes) this writes,
# into requant_v2_outputs/lipid_function_analysis/:
#   <tag>_gene_annotations.csv   one row per gene: full name, EC numbers (reactions),
#                                KEGG pathways, GO biological process / molecular
#                                function / cellular component terms, and TRUE/FALSE
#                                flags for major compartments
#   <tag>_enrichment_<GO_BP|GO_MF|GO_CC|KEGG|Reactome>.csv   Up vs Down over-
#                                representation results (ALL terms; see p.adjust)
#   <tag>_dotplot_<...>.png      dotplot of terms with p.adjust < ENRICH_PADJ
#   <tag>_STRING_edges.csv       protein-protein interactions among the genes (STRING)
#   <tag>_genes_up.txt / _down.txt   symbols, to paste into string-db.org or KEGG Mapper
#
# BACKGROUND: enrichment uses the parent list (all lipid-metabolism genes, or all
# endocytic lipid enzymes) as the universe, NOT the whole genome. All input genes are
# already lipid genes, so a genome-wide background would simply return "lipid
# metabolism". This way a hit means "over-represented among the responding genes
# compared with other lipid genes". Gene sets this small have little statistical power,
# so treat enrichment as hypothesis-generating and lean on the per-gene table.
#
# Needs (Bioconductor): clusterProfiler, enrichplot, org.Mm.eg.db, GO.db
# Optional: ReactomePA (Reactome), STRINGdb (interactions)
# Internet is needed for KEGG and STRING; the GO parts work offline.

suppressPackageStartupMessages({ library(readr); library(dplyr); library(tidyr); library(ggplot2) })
for (pkg in c("clusterProfiler", "enrichplot", "org.Mm.eg.db", "GO.db", "AnnotationDbi"))
  if (!requireNamespace(pkg, quietly = TRUE))
    stop("Package '", pkg, "' is needed:  BiocManager::install('", pkg, "')")
OrgDb <- org.Mm.eg.db::org.Mm.eg.db

# Bioconductor packages mask several dplyr verbs (S4Vectors::first, AnnotationDbi::select,
# ...). Pin the dplyr versions so the script behaves the same however packages load.
select <- dplyr::select; first <- dplyr::first; rename <- dplyr::rename
slice <- dplyr::slice; filter <- dplyr::filter; desc <- dplyr::desc; count <- dplyr::count

ENRICH_PADJ <- 0.05     # terms shown in the dotplots (CSVs keep everything)
EXCLUDE_IEA <- FALSE    # TRUE: GO terms with experimental/curated evidence only
STRING_SCORE <- 400     # STRING confidence cutoff (400 medium, 700 high)

PROJ <- Sys.getenv("PROJ", unset = file.path(
  Sys.getenv("HOME"),
  "Documents/aviti_rnaseq/2025/aviti",
  "20251009-100725_Maury_25161_WivfFacICOTHFAUJwSPeHBEFfCfNtCRpUSgVGUkz",
  "results/Maury_25161"))
OUT <- file.path(PROJ, "requant_v2_outputs")
RES <- file.path(OUT, "lipid_function_analysis")
dir.create(RES, showWarnings = FALSE)

inputs <- list(
  list(tag = "all_lipid",
       file = "DE_v2_IFNB_vs_Mock_lipid_genes_log2FC_gt0.7_padj.csv",
       universe = "DE_v2_IFNB_vs_Mock_protein_coding_lipid_metabolism.csv"),
  list(tag = "endocytic_enzymes",
       file = "DE_v2_IFNB_vs_Mock_lipid_endocytic_enzymes_log2FC_gt0.7_padj.csv",
       universe = "DE_v2_IFNB_vs_Mock_lipid_enzymes_endocytic_vesicle.csv"))

compartments <- c(ER = "GO:0005783", Golgi = "GO:0005794", mitochondrion = "GO:0005739",
                  peroxisome = "GO:0005777", lysosome = "GO:0005764", endosome = "GO:0005768",
                  plasma_membrane = "GO:0005886", lipid_droplet = "GO:0005811",
                  cytosol = "GO:0005829", nucleus = "GO:0005634", extracellular = "GO:0005576")

strip <- function(x) sub("\\.[0-9]+$", "", x)
collapse <- function(x) paste(sort(unique(na.omit(x))), collapse = "; ")

go_names <- function(ids) suppressMessages(AnnotationDbi::select(GO.db::GO.db, keys = unique(ids),
                                                                 keytype = "GOID", columns = "TERM"))

for (job in inputs) {
  f <- file.path(OUT, job$file); fu <- file.path(OUT, job$universe)
  if (!file.exists(f) || !file.exists(fu)) { message("Skipping ", job$tag, ": input not found (run 13/16/27 first)"); next }
  d <- read_csv(f, show_col_types = FALSE)
  uni_tbl <- read_csv(fu, show_col_types = FALSE)
  if (!nrow(d)) { message("Skipping ", job$tag, ": no genes passed the filter"); next }
  d$.key <- strip(d$gene_id); uni_tbl$.key <- strip(uni_tbl$gene_id)
  if (!"direction" %in% names(d)) d$direction <- ifelse(d$log2FoldChange > 0, "up", "down")
  tag <- job$tag
  message("\n=========== ", tag, " : ", nrow(d), " genes (", sum(d$direction == "up"), " up, ",
          sum(d$direction == "down"), " down); background ", nrow(uni_tbl), " genes ===========")

  # ---- Ensembl -> Entrez ------------------------------------------------------------
  map <- suppressMessages(suppressWarnings(
    clusterProfiler::bitr(unique(c(d$.key, uni_tbl$.key)), fromType = "ENSEMBL",
                          toType = "ENTREZID", OrgDb = OrgDb))) |> distinct(ENSEMBL, .keep_all = TRUE)
  ent <- function(keys) unique(map$ENTREZID[match(keys, map$ENSEMBL)] |> na.omit())
  unmapped <- setdiff(d$.key, map$ENSEMBL)
  if (length(unmapped)) message(length(unmapped), " genes have no Entrez ID and are left out of enrichment: ",
                                paste(d$gene_name[d$.key %in% unmapped], collapse = ", "))

  # ---- per-gene annotation ------------------------------------------------------------
  basic <- suppressMessages(AnnotationDbi::select(OrgDb, keys = d$.key, keytype = "ENSEMBL",
                                                  columns = c("GENENAME", "ENZYME"))) |>
    group_by(ENSEMBL) |> summarise(gene_full_name = first(na.omit(GENENAME)), EC_numbers = collapse(ENZYME), .groups = "drop")

  go <- suppressMessages(AnnotationDbi::select(OrgDb, keys = d$.key, keytype = "ENSEMBL",
                                               columns = c("GO", "ONTOLOGY", "EVIDENCE"))) |> filter(!is.na(GO))
  if (EXCLUDE_IEA) go <- go |> filter(EVIDENCE != "IEA")
  go <- go |> left_join(go_names(go$GO), by = c("GO" = "GOID")) |>
    group_by(ENSEMBL, ONTOLOGY) |> summarise(terms = collapse(TERM), .groups = "drop") |>
    mutate(ONTOLOGY = paste0("GO_", ONTOLOGY)) |> pivot_wider(names_from = ONTOLOGY, values_from = terms)
  for (cn in c("GO_BP", "GO_MF", "GO_CC")) if (!cn %in% names(go)) go[[cn]] <- NA_character_

  goall <- suppressMessages(AnnotationDbi::select(OrgDb, keys = d$.key, keytype = "ENSEMBL",
                                                  columns = c("GOALL", "ONTOLOGYALL"))) |>
    filter(ONTOLOGYALL == "CC")
  loc <- tibble(ENSEMBL = d$.key)
  for (nm in names(compartments)) loc[[nm]] <- loc$ENSEMBL %in% goall$ENSEMBL[goall$GOALL == compartments[[nm]]]

  kegg_map <- tibble(ENTREZID = character(), KEGG_pathways = character())
  ek <- tryCatch(clusterProfiler::enrichKEGG(gene = ent(d$.key), organism = "mmu", pvalueCutoff = 1,
                                             qvalueCutoff = 1, minGSSize = 1, maxGSSize = 5000),
                 error = function(e) { message("KEGG not reachable (", conditionMessage(e), ") - KEGG columns left empty"); NULL })
  if (!is.null(ek) && nrow(ek@result)) {
    kegg_map <- ek@result |> select(Description, geneID) |>
      tidyr::separate_rows(geneID, sep = "/") |>
      mutate(Description = sub(" - Mus musculus.*$", "", Description)) |>
      group_by(ENTREZID = geneID) |> summarise(KEGG_pathways = collapse(Description), .groups = "drop")
  }

  ann <- d |> select(gene_id, gene_name, direction, log2FoldChange, padj, .key) |>
    left_join(basic, by = c(".key" = "ENSEMBL")) |>
    left_join(map |> select(ENSEMBL, ENTREZID), by = c(".key" = "ENSEMBL")) |>
    left_join(kegg_map, by = "ENTREZID") |>
    left_join(go, by = c(".key" = "ENSEMBL")) |>
    left_join(loc, by = c(".key" = "ENSEMBL")) |>
    select(-.key) |> relocate(gene_full_name, .after = gene_name) |> arrange(direction, desc(abs(log2FoldChange)))
  write_csv(ann, file.path(RES, paste0(tag, "_gene_annotations.csv")))

  message("Compartment counts (genes annotated to each):")
  print(as.data.frame(ann |> group_by(direction) |> summarise(n = n(), across(all_of(names(compartments)), sum), .groups = "drop")),
        row.names = FALSE)

  # symbol lists for web tools
  for (dir in c("up", "down"))
    writeLines(ann$gene_name[ann$direction == dir], file.path(RES, paste0(tag, "_genes_", dir, ".txt")))

  # ---- enrichment: Up vs Down ----------------------------------------------------------
  gl <- list(Up = ent(d$.key[d$direction == "up"]), Down = ent(d$.key[d$direction == "down"]))
  gl <- gl[lengths(gl) > 0]
  uni <- ent(uni_tbl$.key)
  if (length(gl)) {
    enrich <- function(label, fun, ...) {
      cc <- tryCatch(clusterProfiler::compareCluster(geneClusters = gl, fun = fun, universe = uni,
                                                     pvalueCutoff = 1, qvalueCutoff = 1, ...),
                     error = function(e) { message(label, ": not run (", conditionMessage(e), ")"); NULL })
      if (is.null(cc) || !nrow(as.data.frame(cc))) { message(label, ": no terms returned"); return(invisible()) }
      cc <- tryCatch(clusterProfiler::setReadable(cc, OrgDb, keyType = "ENTREZID"), error = function(e) cc)
      res <- as.data.frame(cc)
      write_csv(res, file.path(RES, paste0(tag, "_enrichment_", label, ".csv")))
      sig <- res |> filter(p.adjust < ENRICH_PADJ)
      message(sprintf("%-9s %d terms tested, %d with p.adjust < %g", label, nrow(res), nrow(sig), ENRICH_PADJ))
      if (nrow(sig)) {
        for (g in unique(sig$Cluster)) message("   ", g, ": ", paste(head(sig$Description[sig$Cluster == g], 5), collapse = " | "))
        cc@compareClusterResult <- sig
        p <- tryCatch(enrichplot::dotplot(cc, showCategory = 10, label_format = 45) +
                        labs(title = paste(tag, "-", label)), error = function(e) NULL)
        if (!is.null(p)) ggsave(file.path(RES, paste0(tag, "_dotplot_", label, ".png")), p, width = 8, height = 6, dpi = 300)
      }
    }
    for (ont in c("BP", "MF", "CC"))
      enrich(paste0("GO_", ont), "enrichGO", OrgDb = OrgDb, keyType = "ENTREZID", ont = ont, minGSSize = 3, maxGSSize = 500)
    enrich("KEGG", "enrichKEGG", organism = "mmu", minGSSize = 3, maxGSSize = 500)
    if (requireNamespace("ReactomePA", quietly = TRUE))
      enrich("Reactome", "enrichPathway", organism = "mouse", minGSSize = 3, maxGSSize = 500)
    else message("Reactome skipped (install ReactomePA to add it)")
  }

  # ---- interactions (STRING) ---------------------------------------------------------------
  if (requireNamespace("STRINGdb", quietly = TRUE)) {
    tryCatch({
      sdb <- STRINGdb::STRINGdb$new(species = 10090, score_threshold = STRING_SCORE,
                                    network_type = "full", input_directory = tempdir())
      mp <- sdb$map(as.data.frame(ann |> select(gene_name, direction)), "gene_name", removeUnmappedRows = TRUE)
      inter <- sdb$get_interactions(unique(mp$STRING_id))
      sym <- setNames(mp$gene_name, mp$STRING_id); dirn <- setNames(mp$direction, mp$STRING_id)
      edges <- inter |> distinct(from, to, .keep_all = TRUE) |>
        transmute(gene_a = sym[from], gene_b = sym[to], dir_a = dirn[from], dir_b = dirn[to], combined_score)
      write_csv(edges, file.path(RES, paste0(tag, "_STRING_edges.csv")))
      deg <- sort(table(c(edges$gene_a, edges$gene_b)), decreasing = TRUE)
      message("STRING: ", nrow(edges), " interactions (score >= ", STRING_SCORE, "); most connected: ",
              paste(sprintf("%s (%d)", names(head(deg, 8)), head(deg, 8)), collapse = ", "))
    }, error = function(e) message("STRING not run (", conditionMessage(e), ")"))
  } else message("STRING skipped (install STRINGdb, or paste the _genes_up/down.txt lists into string-db.org, organism Mus musculus)")
}
message("\nResults in: ", RES)
