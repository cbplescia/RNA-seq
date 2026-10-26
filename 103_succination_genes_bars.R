#!/usr/bin/env Rscript
# SUCCINATION (not succinylation): the covalent modification of protein CYSTEINE residues by fumarate, giving S-(2-succinyl)cysteine. It is a non-enzymatic reaction that rises when fumarate accumulates, so there are
# no "succinating enzymes" - the relevant genes are those that make or remove fumarate, plus proteins classically found succinated. Panels (my curated list; edit SUCC):
#   fumarate sources   Asl (urea cycle), Adsl (purine nucleotide cycle), Fah (tyrosine catabolism), the succinate dehydrogenase complex Sdha / Sdhb / Sdhc / Sdhd (TCA)
#   fumarate sink      Fh1 (fumarate hydratase; loss of it is the classic cause of succination); Mdh2 / Got2 as the neighbouring TCA steps
#   classic succinated proteins   Keap1 (succination activates NRF2), Gapdh, Aco2, Pgk1, Gsr (their expression only - it does not measure the modification)
#   NRF2 response      Nfe2l2, Hmox1, Nqo1 (what Keap1 succination switches on)
# Each gene is also tested against GO for an ENDOLYSOSOMAL link (endosome / lysosome compartment or endocytic / endolysosomal process, as in 102_endolysosomal_succ_ub_genes_bars.R); the endolysosomal figure only
# contains those that qualify. These genes are mostly mitochondrial or cytosolic, so that figure may be empty - the script says so. The evidence column of the csv shows which rule each gene met.
# Writes bars_succination_genes.png/.pdf (all genes found), bars_succination_genes_endolysosomal.png/.pdf (if any gene qualifies) and succination_genes_log2fc.csv to requant_v2_outputs/figures_gene_bars/.
# Usage: Rscript 103_succination_genes_bars.R   (PROJ from 'source setup_session.sh')

suppressPackageStartupMessages({ library(readr); library(dplyr); library(tidyr); library(ggplot2) })
select <- dplyr::select; filter <- dplyr::filter

SUCC <- list(`fumarate sources` = c("Asl", "Adsl", "Fah", "Sdha", "Sdhb", "Sdhc", "Sdhd"), `fumarate sink` = c("Fh1", "Mdh2", "Got2"),
             `classic succinated proteins` = c("Keap1", "Gapdh", "Aco2", "Pgk1", "Gsr"), `NRF2 response` = c("Nfe2l2", "Hmox1", "Nqo1"))
ISG <- c("Uba7", "Ube2l6", "Herc6", "Trim25", "Arih1")
PADJ_CUT <- 0.05; LFC_MIN <- as.numeric(Sys.getenv("LFC_MIN", unset = "0.58")); PER_PAGE <- as.integer(Sys.getenv("PER_PAGE", unset = "36")); UB_ALL <- !identical(Sys.getenv("SIG_ONLY"), "1"); LOW_EXPR <- 50
COMP <- c(IFNB_vs_Mock = "IFNB vs Mock", EBOV_vs_Mock = "EBOV vs Mock"); CCOL <- c(`IFNB vs Mock` = "#2B6CB0", `EBOV vs Mock` = "#1B7837"); INK <- "#1F2933"
PROJ <- Sys.getenv("PROJ", unset = file.path(
  Sys.getenv("HOME"),
  "Documents/aviti_rnaseq/2025/aviti",
  "20251009-100725_Maury_25161_WivfFacICOTHFAUJwSPeHBEFfCfNtCRpUSgVGUkz",
  "results/Maury_25161"))
RNA <- file.path(PROJ, "requant_v2_outputs")
OUT <- file.path(RNA, "figures_gene_bars"); dir.create(OUT, showWarnings = FALSE)
stars <- function(p) ifelse(is.na(p), "", ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", "ns"))))
tabs <- lapply(names(COMP), function(nm) read_csv(file.path(RNA, paste0("DE_v2_", nm, "_annotated.csv")), show_col_types = FALSE) |> mutate(key = tolower(gene_name), gid = sub("\\..*$", "", gene_id))); names(tabs) <- names(COMP)

if (!requireNamespace("clusterProfiler", quietly = TRUE) || !requireNamespace("org.Mm.eg.db", quietly = TRUE)) stop("Install clusterProfiler and org.Mm.eg.db")
g2a <- as.list(org.Mm.eg.db::org.Mm.egGO2ALLEGS)
ens_of <- function(go) { e <- unique(as.character(unlist(g2a[go]))); e <- e[!is.na(e)]; if (!length(e)) return(character()); unique(suppressMessages(suppressWarnings(clusterProfiler::bitr(e, "ENTREZID", "ENSEMBL", org.Mm.eg.db::org.Mm.eg.db)))$ENSEMBL) }
loc_ens <- ens_of(c("GO:0005768", "GO:0005769", "GO:0005770", "GO:0005771", "GO:0010008", "GO:0036019", "GO:0005764", "GO:0005765", "GO:0043202", "GO:0030139", "GO:0045334", "GO:0031901", "GO:0031902"))
proc_ens <- ens_of(c("GO:0006897", "GO:0016197", "GO:0007041", "GO:0007040", "GO:0006622", "GO:0030100"))
CUR <- c("Nedd4", "Nedd4l", "Itch", "Smurf1", "Smurf2", "Wwp1", "Wwp2", "Cbl", "Cblb", "Mib1", "Rnf13", "Rnf26", "Rnf152", "Rnf167", "Mgrn1", "Znrf1", "Znrf2", "Znrf3", "Rnf43", "Stub1", "Marchf1", "Marchf2", "Marchf3", "Marchf8", "Marchf9", "Marchf11", "Hecw1", "Hecw2")
message(sprintf("endolysosomal GO sets: compartment %d genes, process %d genes", length(loc_ens), length(proc_ens)))
evidence <- function(ens, name) { v <- c(if (ens %in% loc_ens) "GO compartment", if (ens %in% proc_ens) "GO process", if (tolower(name) %in% tolower(CUR)) "curated"); paste(v, collapse = " + ") }
# data frame of the three-column results for a gene table GT (columns panel, gene_id, pos)
collect <- function(GT) bind_rows(lapply(names(COMP), function(nm) tabs[[nm]] |> inner_join(GT, by = "gene_id") |>
  transmute(panel, gene = gene_name, gene_id, pos, baseMean, comparison = COMP[[nm]], lfc = log2FoldChange, se = lfcSE, padj))) |>
  mutate(comparison = factor(comparison, levels = COMP), lab = stars(padj), vj = ifelse(lfc >= 0, 0, 1), ys = lfc + ifelse(lfc >= 0, 1, -1) * (se + 0.03), xl = paste0(gene, ifelse(baseMean < LOW_EXPR, "†", "")))
draw <- function(D, title, fname, extra = "") {
  D$panel <- factor(D$panel, levels = unique(D$panel[order(D$pos)])); D$key <- paste(D$panel, D$xl, sep = "||"); D$key <- factor(D$key, levels = unique(D$key[order(D$pos)]))
  dg <- position_dodge(width = 0.85)
  p <- ggplot(D, aes(key, lfc, fill = comparison)) + geom_col(position = dg, width = 0.8) + geom_errorbar(aes(ymin = lfc - se, ymax = lfc + se), position = dg, width = 0.25, colour = INK, linewidth = 0.3) +
    geom_text(aes(y = ys, label = lab, vjust = vj), position = dg, colour = INK, size = 2.2) + geom_hline(yintercept = 0, colour = INK, linewidth = 0.4) +
    facet_grid(. ~ panel, scales = "free_x", space = "free_x", labeller = labeller(panel = label_wrap_gen(15))) + scale_fill_manual(values = CCOL, name = NULL) +
    scale_x_discrete(labels = function(x) sub("^.*\\|\\|", "", x)) + scale_y_continuous(expand = expansion(mult = 0.12)) +
    labs(title = title, subtitle = sprintf("IFNB vs Mock and EBOV vs Mock side by side; shrunk log2FC +/- 1 SE; stars = padj of each comparison (ns = not significant); dagger = baseMean < %d%s", LOW_EXPR, extra), x = NULL, y = "log2 fold-change vs Mock") +
    theme_bw(base_size = 10) + theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(), panel.grid.major.y = element_line(colour = "grey92", linewidth = 0.3), panel.border = element_rect(colour = "grey70"),
          axis.text = element_text(colour = INK), axis.text.x = element_text(face = "italic", angle = 50, hjust = 1, size = 7.5), plot.title = element_text(colour = INK, face = "bold", size = 11), plot.subtitle = element_text(colour = "grey35", size = 7.5),
          strip.background = element_rect(fill = "grey95", colour = "grey70"), strip.text = element_text(colour = INK, face = "bold", size = 7), legend.position = "top", panel.spacing.x = unit(0.15, "lines"))
  w <- max(6.5, 0.48 * nlevels(D$key) + 1.2 * nlevels(D$panel) + 2)
  ggsave(file.path(OUT, paste0(fname, ".png")), p, width = w, height = 5.5, dpi = 300, limitsize = FALSE); ggsave(file.path(OUT, paste0(fname, ".pdf")), p, width = w, height = 5.5, limitsize = FALSE)
}
find_gene <- function(q) { h <- tabs[[1]] |> filter(key == tolower(q)) |> distinct(gene_id, .keep_all = TRUE); if (nrow(h)) h$gene_id[1] else NA_character_ }

# ================= succination genes =================
GT <- bind_rows(lapply(names(SUCC), function(p) tibble(panel = p, query = SUCC[[p]]))) |> distinct(query, .keep_all = TRUE); GT$gene_id <- vapply(GT$query, find_gene, "")
miss <- GT$query[is.na(GT$gene_id)]; if (length(miss)) message("Not found in the DE tables (removed by the pre-filter, or a different symbol): ", paste(miss, collapse = ", "))
GT <- GT |> filter(!is.na(gene_id)) |> mutate(ev = vapply(seq_len(n()), function(i) evidence(sub("\\..*$", "", gene_id[i]), query[i]), ""), pos = row_number())
message("succination genes detected: ", nrow(GT), "; with an endolysosomal link: ", sum(nzchar(GT$ev)), if (any(nzchar(GT$ev))) paste0(" (", paste(GT$query[nzchar(GT$ev)], collapse = ", "), ")") else "")
if (!nrow(GT)) stop("No succination gene found in the DE tables")
D1 <- collect(GT)
write_csv(D1 |> arrange(pos, comparison) |> left_join(GT |> select(gene_id, endolysosomal_evidence = ev), by = "gene_id") |> select(panel, gene, gene_id, endolysosomal_evidence, baseMean, comparison, log2FC = lfc, lfcSE = se, padj), file.path(OUT, "succination_genes_log2fc.csv"))
draw(D1, "Succination-related genes (fumarate metabolism): IFNB vs Mock and EBOV vs Mock", "bars_succination_genes", sprintf("; %d genes found, %d not in the data", nrow(GT), length(miss)))
GE <- GT |> filter(nzchar(ev)) |> mutate(pos = row_number())
if (nrow(GE)) {
  draw(collect(GE), "Succination-related genes with an endolysosomal link: IFNB vs Mock and EBOV vs Mock", "bars_succination_genes_endolysosomal", sprintf("; %d genes", nrow(GE)))
  message(sprintf("endolysosomal succination figure: %d genes", nrow(GE)))
} else message("None of the succination-related genes has an endolysosomal annotation in GO - no endolysosomal figure (the full set is in bars_succination_genes)")
message("Saved to: ", OUT)
