#!/usr/bin/env Rscript
# From the lipid-metabolism genes (script 13), pick the ENZYMES that are known or
# predicted to sit in the membrane or lumen of endocytic vesicles / endosomes,
# write a CSV, and draw a log2FC bar graph (IFNB vs Mock).
#
# How the list is compiled (all from Gene Ontology via org.Mm.eg.db, which carries
# every gene's parent terms, so a gene annotated to e.g. "early endosome membrane"
# is also counted under "endosome"):
#   enzyme       : annotated to Molecular Function GO:0003824 "catalytic activity"
#                  (or any child of it)
#   vesicle      : Cellular Component annotation to endocytic vesicle (GO:0030139),
#                  endosome (GO:0005768; includes early, late and recycling endosome),
#                  phagocytic vesicle (GO:0045335) or multivesicular body (GO:0005771)
#   lysosome     : GO:0005764 - the terminal compartment of the endocytic pathway.
#                  Included by default (INCLUDE_LYSOSOME <- TRUE); the lysosome_only
#                  column marks genes with no earlier endosome/vesicle annotation, so
#                  you can drop them again, or set INCLUDE_LYSOSOME <- FALSE.
#   EXTRA_GENES  : genes to consider even if GO does not call them lipid-metabolism
#                  genes (e.g. NAGLU, a lysosomal heparan-sulfate enzyme). They still
#                  have to pass the enzyme + endocytic/lysosome checks; the
#                  in_lipid_go_list column says which genes came in this way.
#   location     : "membrane" / "lumen" from specific GO terms (listed below);
#                  "membrane and lumen" if both; "not specified" if the gene is
#                  annotated to the compartment but not to a membrane/lumen term
# Annotations are a mix of experimental and computationally inferred (IEA); set
# EXCLUDE_IEA <- TRUE for curated/experimental evidence only. GO describes where a
# protein has been reported or predicted, not where it is exclusively - many lipid
# enzymes also work in the ER or Golgi.
#
# Input : DE_v2_IFNB_vs_Mock_protein_coding_lipid_metabolism.csv  (from 13)
# Output: DE_v2_IFNB_vs_Mock_lipid_enzymes_endocytic_vesicle.csv
#         DE_v2_IFNB_vs_Mock_lipid_endocytic_vesicle_non_enzymes.csv (for reference)
#         figures_lipid_pathways/lipid_enzymes_endocytic_vesicle_log2FC.png

suppressPackageStartupMessages({ library(readr); library(dplyr); library(ggplot2) })
for (pkg in c("org.Mm.eg.db", "AnnotationDbi"))
  if (!requireNamespace(pkg, quietly = TRUE))
    stop("Package '", pkg, "' is needed:  BiocManager::install('org.Mm.eg.db')")

INCLUDE_LYSOSOME <- TRUE
EXTRA_GENES      <- c("Naglu")   # symbols; add others here
EXCLUDE_IEA      <- FALSE
LFC_CUT <- 1; PADJ_CUT <- 0.05
INK <- "#1F2933"
COL <- c(Up = "#C53030", Down = "#2B6CB0", `Not significant` = "#B8B8B8")

PROJ <- Sys.getenv("PROJ", unset = file.path(
  Sys.getenv("HOME"),
  "Documents/aviti_rnaseq/2025/aviti",
  "20251009-100725_Maury_25161_WivfFacICOTHFAUJwSPeHBEFfCfNtCRpUSgVGUkz",
  "results/Maury_25161"))
OUT <- file.path(PROJ, "requant_v2_outputs")
FIG <- file.path(OUT, "figures_lipid_pathways")
dir.create(FIG, showWarnings = FALSE)

lip <- read_csv(file.path(OUT, "DE_v2_IFNB_vs_Mock_protein_coding_lipid_metabolism.csv"), show_col_types = FALSE)
lip$in_lipid_go_list <- TRUE
if (length(EXTRA_GENES)) {
  pc <- read_csv(file.path(OUT, "DE_v2_IFNB_vs_Mock_protein_coding.csv"), show_col_types = FALSE)
  add <- pc |> filter(toupper(gene_name) %in% toupper(EXTRA_GENES), !gene_id %in% lip$gene_id)
  nf <- setdiff(toupper(EXTRA_GENES), toupper(c(lip$gene_name, add$gene_name)))
  if (length(nf)) message("EXTRA_GENES not found among protein-coding genes in the data: ", paste(nf, collapse = ", "))
  if (nrow(add)) { add$in_lipid_go_list <- FALSE; lip <- bind_rows(lip, add)
                   message("Added from EXTRA_GENES: ", paste(add$gene_name, collapse = ", ")) }
}
lip$.key <- sub("\\.[0-9]+$", "", lip$gene_id)

go <- suppressMessages(AnnotationDbi::select(
  org.Mm.eg.db::org.Mm.eg.db, keys = unique(lip$.key), keytype = "ENSEMBL",
  columns = c("GOALL", "ONTOLOGYALL", "EVIDENCEALL"))) |> filter(!is.na(GOALL))
if (EXCLUDE_IEA) go <- go |> filter(EVIDENCEALL != "IEA")
has <- function(ont, ids) unique(go$ENSEMBL[go$ONTOLOGYALL == ont & go$GOALL %in% ids])

compartments <- c(
  endocytic_vesicle    = "GO:0030139",
  endosome             = "GO:0005768",
  early_endosome       = "GO:0005769",
  late_endosome        = "GO:0005770",
  recycling_endosome   = "GO:0055037",
  multivesicular_body  = "GO:0005771",
  phagocytic_vesicle   = "GO:0045335",
  lysosome             = "GO:0005764")
MEMBRANE <- c("GO:0010008", "GO:0030666", "GO:0031901", "GO:0031902", "GO:0055038",
              "GO:0032585", "GO:0030670", "GO:0005765")   # endosome/endocytic/early/late/recycling/MVB/phagocytic vesicle/lysosomal membrane
LUMEN    <- c("GO:0031904", "GO:0071682", "GO:0097013", "GO:0043202")  # endosome/endocytic vesicle/phagocytic vesicle/lysosomal lumen

for (nm in names(compartments)) lip[[nm]] <- lip$.key %in% has("CC", compartments[[nm]])
lip$is_enzyme  <- lip$.key %in% has("MF", "GO:0003824")
lip$in_membrane <- lip$.key %in% has("CC", MEMBRANE)
lip$in_lumen    <- lip$.key %in% has("CC", LUMEN)

vesicle_cols <- c("endocytic_vesicle", "endosome", "multivesicular_body", "phagocytic_vesicle")
lip$vesicle <- Reduce(`|`, lapply(vesicle_cols, function(n) lip[[n]]))
lip$lysosome_only <- lip$lysosome & !lip$vesicle
lip$in_scope <- lip$vesicle | (INCLUDE_LYSOSOME & lip$lysosome)

lip$location <- case_when(lip$in_membrane & lip$in_lumen ~ "membrane and lumen",
                          lip$in_membrane ~ "membrane", lip$in_lumen ~ "lumen",
                          TRUE ~ "not specified")
comp_cols <- names(compartments)
lip$compartments <- apply(lip[, comp_cols], 1, function(r) paste(comp_cols[as.logical(r)], collapse = ", "))
lip$direction <- case_when(!is.na(lip$padj) & lip$padj < PADJ_CUT & lip$log2FoldChange >  LFC_CUT ~ "Up",
                           !is.na(lip$padj) & lip$padj < PADJ_CUT & lip$log2FoldChange < -LFC_CUT ~ "Down",
                           TRUE ~ "Not significant")

keep_cols <- function(d) d |> select(-.key) |> arrange(padj, desc(abs(log2FoldChange)))
enz <- lip |> filter(in_scope, is_enzyme)
non <- lip |> filter(in_scope, !is_enzyme)
write_csv(keep_cols(enz), file.path(OUT, "DE_v2_IFNB_vs_Mock_lipid_enzymes_endocytic_vesicle.csv"))
write_csv(keep_cols(non), file.path(OUT, "DE_v2_IFNB_vs_Mock_lipid_endocytic_vesicle_non_enzymes.csv"))

message("Lipid-metabolism genes: ", nrow(lip))
message("  in endocytic vesicle / endosome / phagosome / MVB", if (INCLUDE_LYSOSOME) " / lysosome" else "", ": ", sum(lip$in_scope))
message("  of these, enzymes: ", nrow(enz), "   (non-enzymes, written separately: ", nrow(non), ")")
print(as.data.frame(enz |> count(location)), row.names = FALSE)
message("  significant: up ", sum(enz$direction == "Up"), ", down ", sum(enz$direction == "Down"))
if (!nrow(enz)) stop("No genes passed the filters.")

# ---- bar graph, split by membrane / lumen ---------------------------------------
d <- enz |> arrange(log2FoldChange) |>
  mutate(label = make.unique(ifelse(is.na(gene_name), gene_id, gene_name)),
         location = factor(location, levels = c("membrane", "lumen", "membrane and lumen", "not specified")),
         direction = factor(direction, levels = names(COL)))
d$label <- factor(d$label, levels = d$label)
has_se <- "lfcSE" %in% names(d)
p <- ggplot(d, aes(log2FoldChange, label, fill = direction)) +
  geom_vline(xintercept = 0, colour = "grey40", linewidth = 0.4) +
  geom_col(width = 0.7) +
  { if (has_se) geom_errorbarh(aes(xmin = log2FoldChange - lfcSE, xmax = log2FoldChange + lfcSE),
                               height = 0.3, colour = INK, linewidth = 0.3) } +
  facet_grid(location ~ ., scales = "free_y", space = "free_y", drop = TRUE) +
  scale_fill_manual(values = COL, drop = FALSE, name = NULL) +
  scale_x_continuous(expand = expansion(mult = 0.08)) +
  labs(title = "Lipid-metabolism enzymes at endocytic vesicles / endosomes",
       subtitle = sprintf("IFNB vs Mock, shrunk log2FC%s; %d enzymes, %d significant (padj<%g & |log2FC|>%g); panels = GO membrane / lumen annotation",
                          if (has_se) " ± lfcSE" else "", nrow(d), sum(d$direction != "Not significant"), PADJ_CUT, LFC_CUT),
       x = "log2 fold-change (IFNB vs Mock)", y = NULL) +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_line(colour = "grey92", linewidth = 0.3),
        panel.border = element_rect(colour = "grey70"),
        strip.background = element_rect(fill = "grey95", colour = "grey70"),
        strip.text = element_text(colour = INK, face = "bold"),
        axis.text = element_text(colour = INK), axis.text.y = element_text(face = "italic"),
        axis.title = element_text(colour = INK),
        plot.title = element_text(colour = INK, face = "bold"),
        plot.subtitle = element_text(colour = "grey35", size = 8.5),
        legend.position = "bottom")
ggsave(file.path(FIG, "lipid_enzymes_endocytic_vesicle_log2FC.png"), p,
       width = 7.5, height = max(4, 0.22 * nrow(d) + 3), dpi = 300, limitsize = FALSE)
message("\nCSV and PNG written to: ", OUT, "  and  ", FIG)
