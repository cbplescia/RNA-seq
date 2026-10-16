#!/usr/bin/env Rscript
# Pull the enzyme-encoding genes out of the IFNB vs Mock results (protein-coding
# genes) and split the responders into up / down, ready for 28_lipid_function_analysis.R.
#
# What counts as an enzyme (ENZYME_DEF below):
#   "EC"     the gene has an Enzyme Commission number in org.Mm.eg.db (a specific reaction)
#   "GO"     annotated (directly or via a child term) to GO:0003824 "catalytic activity"
#   "either" (default) EC OR GO. GO catalytic activity is broader: it also picks up
#            kinases, phosphatases, ubiquitin ligases, ATPases and similar.
# Responders: padj < PADJ_CUT and |log2FC| > LFC_CUT (shrunk IFNB vs Mock).
#
# Writes to requant_v2_outputs/:
#   DE_v2_IFNB_vs_Mock_protein_coding_enzymes.csv   all enzyme genes tested (the background)
#   DE_v2_IFNB_vs_Mock_enzymes_up_down.csv          the responding enzymes, with direction
# and a figure of responders by EC class in lipid_function_analysis/.
# Then run:  Rscript 28_lipid_function_analysis.R all_enzymes

suppressPackageStartupMessages({ library(readr); library(dplyr); library(tidyr); library(ggplot2) })
for (pkg in c("org.Mm.eg.db", "AnnotationDbi"))
  if (!requireNamespace(pkg, quietly = TRUE)) stop("Install ", pkg, " (BiocManager::install)")
OrgDb <- org.Mm.eg.db::org.Mm.eg.db
select <- dplyr::select; first <- dplyr::first; rename <- dplyr::rename
slice <- dplyr::slice; filter <- dplyr::filter; desc <- dplyr::desc; count <- dplyr::count

ENZYME_DEF <- "either"          # "EC", "GO" or "either"
LFC_CUT <- 1; PADJ_CUT <- 0.05  # use 0.7 to match the lipid analysis

PROJ <- Sys.getenv("PROJ", unset = file.path(
  Sys.getenv("HOME"),
  "Documents/aviti_rnaseq/2025/aviti",
  "20251009-100725_Maury_25161_WivfFacICOTHFAUJwSPeHBEFfCfNtCRpUSgVGUkz",
  "results/Maury_25161"))
OUT <- file.path(PROJ, "requant_v2_outputs")
RES <- file.path(OUT, "lipid_function_analysis")
dir.create(RES, showWarnings = FALSE)

pc <- read_csv(file.path(OUT, "DE_v2_IFNB_vs_Mock_protein_coding.csv"), show_col_types = FALSE)
pc$.key <- sub("\\.[0-9]+$", "", pc$gene_id)
message("Protein-coding genes tested: ", nrow(pc))

ec <- suppressMessages(AnnotationDbi::select(OrgDb, keys = unique(pc$.key), keytype = "ENSEMBL", columns = "ENZYME")) |>
  filter(!is.na(ENZYME)) |> distinct(ENSEMBL, ENZYME)
ec_tab <- ec |> group_by(ENSEMBL) |>
  summarise(EC_numbers = paste(sort(ENZYME), collapse = "; "),
            ec_class_n = list(sort(unique(as.integer(sub("\\..*$", "", ENZYME))))), .groups = "drop")
ec_names <- c("1" = "Oxidoreductase", "2" = "Transferase", "3" = "Hydrolase", "4" = "Lyase",
              "5" = "Isomerase", "6" = "Ligase", "7" = "Translocase")
ec_tab$EC_class <- vapply(ec_tab$ec_class_n, function(v) paste(ec_names[as.character(v)], collapse = " / "), character(1))
ec_tab$ec_class_n <- NULL

go_cat <- suppressMessages(AnnotationDbi::select(OrgDb, keys = unique(pc$.key), keytype = "ENSEMBL",
                                                 columns = c("GOALL", "ONTOLOGYALL"))) |>
  filter(ONTOLOGYALL == "MF", GOALL == "GO:0003824") |> distinct(ENSEMBL)

pc <- pc |> left_join(ec_tab, by = c(".key" = "ENSEMBL")) |>
  mutate(has_EC = !is.na(EC_numbers), has_GO_catalytic = .key %in% go_cat$ENSEMBL,
         is_enzyme = switch(ENZYME_DEF, EC = has_EC, GO = has_GO_catalytic, either = has_EC | has_GO_catalytic))
message(sprintf("Enzyme genes (%s): %d   [EC number: %d, GO catalytic activity: %d, both: %d]",
                ENZYME_DEF, sum(pc$is_enzyme), sum(pc$has_EC), sum(pc$has_GO_catalytic),
                sum(pc$has_EC & pc$has_GO_catalytic)))

enz <- pc |> filter(is_enzyme) |> select(-.key) |> arrange(padj)
write_csv(enz, file.path(OUT, "DE_v2_IFNB_vs_Mock_protein_coding_enzymes.csv"))

resp <- enz |> filter(!is.na(padj), padj < PADJ_CUT, abs(log2FoldChange) > LFC_CUT) |>
  mutate(direction = ifelse(log2FoldChange > 0, "up", "down")) |>
  relocate(direction, .after = padj) |> arrange(desc(log2FoldChange))
write_csv(resp, file.path(OUT, "DE_v2_IFNB_vs_Mock_enzymes_up_down.csv"))
message(sprintf("Responding enzymes (padj<%g, |log2FC|>%g): %d  (up %d, down %d)", PADJ_CUT, LFC_CUT,
                nrow(resp), sum(resp$direction == "up"), sum(resp$direction == "down")))

# ---- responders by EC class ---------------------------------------------------------------
cls <- resp |> mutate(EC_class = ifelse(is.na(EC_class), "No EC number (GO only)", EC_class)) |>
  count(direction, EC_class)
if (nrow(cls)) {
  print(as.data.frame(cls |> tidyr::pivot_wider(names_from = direction, values_from = n, values_fill = 0)), row.names = FALSE)
  p <- ggplot(cls, aes(reorder(EC_class, n, FUN = sum), n, fill = direction)) +
    geom_col(position = position_dodge2(preserve = "single"), width = 0.7) + coord_flip() +
    scale_fill_manual(values = c(up = "#C53030", down = "#2B6CB0"), name = NULL) +
    labs(title = "Responding enzymes by reaction class (EC)",
         subtitle = sprintf("IFNB vs Mock, padj<%g & |log2FC|>%g; genes with several EC classes are labelled with all", PADJ_CUT, LFC_CUT),
         x = NULL, y = "genes") +
    theme_bw(base_size = 11) + theme(panel.grid.minor = element_blank(), panel.border = element_rect(colour = "grey70"),
                                     legend.position = "bottom")
  ggsave(file.path(RES, "all_enzymes_EC_class_counts.png"), p, width = 8, height = 5, dpi = 300)
}
message("\nNext:  Rscript 28_lipid_function_analysis.R all_enzymes")
