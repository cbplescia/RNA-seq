#!/usr/bin/env Rscript
# The genes of script 101 (lysine succinylation genes, and the GO-defined ubiquitin conjugation / ligase genes) filtered to those KNOWN to be localized to, or to act on, the endocytic / endolysosomal system, then
# the same vertical bar graphs (IFNB vs Mock next to EBOV vs Mock). A gene is kept if it is
#   - annotated (GO, org.Mm.eg.db, children included) to an endosomal / lysosomal / endocytic-vesicle compartment: endosome GO:0005768, early / late endosome, multivesicular body, endosome membrane, endolysosome,
#     lysosome GO:0005764, lysosomal membrane / lumen, endocytic and clathrin-coated endocytic vesicle; OR
#   - annotated to an endocytic / endolysosomal PROCESS: endocytosis GO:0006897, endosomal transport GO:0016197, lysosomal transport GO:0007041, lysosome organization GO:0007040, protein targeting to lysosome
#     GO:0006622, regulation of endocytosis GO:0030100; OR
#   - on my short curated list of ubiquitin ligases with established endocytic / lysosomal roles that GO can miss (Nedd4, Nedd4l, Itch, Smurf1/2, Wwp1/2, Cbl, Cblb, Mib1, Rnf13, Rnf26, Rnf152, Rnf167, Mgrn1,
#     Znrf1/2/3, Rnf43, Stub1, Marchf1/2/3/8/9/11, Hecw1/2).
# The evidence for each gene (which of the three) is in the csv files. GO annotation is broad ("affects" includes regulation), so read the list as "has a reported endolysosomal link", not as a verdict.
# Writes bars_succinylation_genes_endolysosomal.png/.pdf, bars_ubiquitination_genes_endolysosomal(_p<n>).png/.pdf and the csv files to requant_v2_outputs/figures_gene_bars/.
# All kept ubiquitin genes are drawn by default (SIG_ONLY=1 keeps only those with padj < PADJ_CUT and |log2FC| > LFC_MIN in IFNB and / or EBOV); PER_PAGE genes per figure.
# Usage: Rscript 102_endolysosomal_succ_ub_genes_bars.R   (PROJ from 'source setup_session.sh';  SIG_ONLY=1  LFC_MIN=0.58  PER_PAGE=36 optional)

suppressPackageStartupMessages({ library(readr); library(dplyr); library(tidyr); library(ggplot2) })
select <- dplyr::select; filter <- dplyr::filter

SUCC <- list(`succinyltransferases (reported)` = c("Kat2a", "Hat1", "Cpt1a", "Ep300"), `alpha-KGDH complex (succinyl-CoA)` = c("Ogdh", "Ogdhl", "Dlst", "Dld"),
             `succinyl-CoA synthetase` = c("Sucla2", "Suclg1", "Suclg2"), `desuccinylases (sirtuins)` = c("Sirt5", "Sirt7", "Sirt3"))
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

# ================= (1) succinylation =================
GT <- bind_rows(lapply(names(SUCC), function(p) tibble(panel = p, query = SUCC[[p]]))) |> distinct(query, .keep_all = TRUE); GT$gene_id <- vapply(GT$query, find_gene, "")
miss <- GT$query[is.na(GT$gene_id)]; if (length(miss)) message("Succinylation genes not found in the DE tables (removed by the pre-filter, or a different symbol): ", paste(miss, collapse = ", "))
GT <- GT |> filter(!is.na(gene_id)) |> mutate(ev = vapply(seq_len(n()), function(i) evidence(sub("\\..*$", "", gene_id[i]), query[i]), ""))
message("succinylation genes detected: ", nrow(GT), "; with an endolysosomal link: ", sum(nzchar(GT$ev)), if (any(nzchar(GT$ev))) paste0(" (", paste(GT$query[nzchar(GT$ev)], collapse = ", "), ")") else "")
GT <- GT |> filter(nzchar(ev)) |> mutate(pos = row_number())
if (nrow(GT)) {
  D1 <- collect(GT); write_csv(D1 |> arrange(pos, comparison) |> left_join(GT |> select(gene_id, evidence = ev), by = "gene_id") |> select(panel, gene, gene_id, evidence, baseMean, comparison, log2FC = lfc, lfcSE = se, padj), file.path(OUT, "succinylation_genes_endolysosomal_log2fc.csv"))
  draw(D1, "Lysine succinylation genes with an endolysosomal link: IFNB vs Mock and EBOV vs Mock", "bars_succinylation_genes_endolysosomal", sprintf("; %d genes found, %d not in the data", nrow(GT), length(miss)))
  message(sprintf("succinylation: %d genes drawn", nrow(GT)))
} else message("No succinylation gene has an endolysosomal annotation - no succinylation figure")

# ================= (2) ubiquitination =================
e1 <- ens_of("GO:0004839"); e2 <- ens_of("GO:0061631"); e3 <- ens_of(c("GO:0061630", "GO:0004842")); all_ub <- unique(c(e1, e2, e3))
message(sprintf("GO ubiquitin transferase genes: E1 %d, E2 %d, E3 / other transferase %d (union %d)", length(e1), length(e2), length(e3), length(all_ub)))
U <- tabs[[1]] |> filter(gid %in% all_ub) |> distinct(gene_id, .keep_all = TRUE) |> transmute(gene_id, gene = gene_name, gid)
cur_ub <- tabs[[1]] |> filter(key %in% tolower(CUR)) |> distinct(gene_id, .keep_all = TRUE) |> transmute(gene_id, gene = gene_name, gid)
U <- bind_rows(U, cur_ub |> filter(!gene_id %in% U$gene_id))
U$ev <- vapply(seq_len(nrow(U)), function(i) evidence(U$gid[i], U$gene[i]), "")
message(sprintf("ubiquitin genes detected (GO + curated): %d; with an endolysosomal link: %d", nrow(U), sum(nzchar(U$ev))))
U <- U |> filter(nzchar(ev))
fam <- function(g, e) {
  if (g %in% ISG) return("ISG15 conjugation (ISGylation)")
  if (e %in% e1) return("E1 (activating)")
  if (grepl("^(Arih|Rnf31$|Rnf144|Rnf19|Prkn$|Rnf216$|Ankib)", g)) return("E3: RBR")
  if (grepl("^Trim[0-9]", g)) return("E3: TRIM")
  if (grepl("^(Herc|Hect|Huwe|Hace|Smurf|Wwp|Nedd4|Itch$|Ube3|Hecw|Trip12$)", g)) return("E3: HECT")
  if (grepl("^(Fbx|Cul[0-9]|Rbx1$|Ddb1$|Klhl|Kctd|Socs|Dcaf|Spop|Btrc$|Skp|Vhl$)", g)) return("E3: cullin-RING / adaptors")
  if (grepl("^Rnf[0-9]", g)) return("E3: RNF / other RING")
  if (e %in% e2 && grepl("^Ube2|^Ube[0-9]", g)) return("E2 (conjugating)")
  if (e %in% e2 && !(e %in% e3)) return("E2 (conjugating)")
  "E3: other"
}
ORDER <- c("ISG15 conjugation (ISGylation)", "E1 (activating)", "E2 (conjugating)", "E3: TRIM", "E3: RNF / other RING", "E3: RBR", "E3: HECT", "E3: cullin-RING / adaptors", "E3: other")
U$panel <- vapply(seq_len(nrow(U)), function(i) fam(U$gene[i], U$gid[i]), ""); U$panel <- factor(U$panel, levels = ORDER)
U0 <- collect(U |> arrange(panel, gene) |> mutate(pos = row_number()) |> mutate(panel = as.character(panel)))
write_csv(U0 |> arrange(pos, comparison) |> left_join(U |> select(gene_id, evidence = ev), by = "gene_id") |> select(panel, gene, gene_id, evidence, baseMean, comparison, log2FC = lfc, lfcSE = se, padj), file.path(OUT, "ubiquitination_genes_endolysosomal_log2fc.csv"))
sig <- U0 |> group_by(gene_id) |> summarise(keep = any(!is.na(padj) & padj < PADJ_CUT & abs(lfc) > LFC_MIN), .groups = "drop")
nsig <- sum(sig$keep); message(sprintf("ubiquitination: %d GO genes detected in your data; %d significantly changed (padj < %g and |log2FC| > %g) in IFNB and / or EBOV", nrow(sig), nsig, PADJ_CUT, LFC_MIN))
U1 <- if (UB_ALL) U0 else U0 |> filter(gene_id %in% sig$gene_id[sig$keep])
if (!nrow(U1)) stop("No endolysosomal ubiquitination gene passes the filter - use LFC_MIN=0 or leave SIG_ONLY unset")
ug <- U1 |> distinct(gene_id, panel, pos) |> arrange(pos) |> mutate(page = ceiling(row_number() / PER_PAGE)); np <- max(ug$page)
for (pg in seq_len(np)) {
  d <- U1 |> inner_join(ug |> select(gene_id, page), by = "gene_id") |> filter(page == pg) |> mutate(pos = match(gene_id, unique(gene_id[order(pos)])))
  draw(d, "Ubiquitin ligase / conjugation genes with an endolysosomal link: IFNB vs Mock and EBOV vs Mock", if (np > 1) sprintf("bars_ubiquitination_genes_endolysosomal_p%d", pg) else "bars_ubiquitination_genes_endolysosomal",
       sprintf("; %s%s", if (UB_ALL) sprintf("all %d detected genes", nrow(sig)) else sprintf("%d of %d detected genes (padj < %g and |log2FC| > %g in IFNB and / or EBOV)", nsig, nrow(sig), PADJ_CUT, LFC_MIN), if (np > 1) sprintf("; page %d of %d", pg, np) else ""))
}
message("Saved to: ", OUT)
