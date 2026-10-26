#!/usr/bin/env Rscript
# Everything STRING (mouse, taxon 10090) links to Gba1 (direct partners, any evidence channel: experiments, databases, text mining, co-expression, neighbourhood, fusion, co-occurrence), then the three-bar graph
# (left to right IFNB, EBOV, IFNB+EBOV vs Mock; shrunk log2FC +/- 1 lfcSE; stars = padj vs Mock; dagger = baseMean < 50) for every partner found in the DE tables, ordered by STRING combined score (strongest first), paged.
# Env: STRING_SCORE (default 400 = medium confidence; 150 = low, more partners; 700 = high), QUERY (default "Gba1", falls back to "Gba"), PER_PAGE (default 20), ONLY_SIG=1 (keep partners with padj < 0.05 in at least one comparison), STRING_VERSION (default 12.0).
# Needs STRINGdb + internet the first time (mouse data cached in <project>/stringdb_cache). A STRING link means "functional association" (physical or not), not necessarily direct binding.
# Writes bars_gba1_partners_pageN.png/.pdf, gba1_string_partners.csv (partner, score, whether found in the DE tables) and gba1_partners_three_comparisons.csv to requant_v2_outputs/figures_gba1_string/.   Usage: Rscript 111_gba1_string_partners_bars.R

suppressPackageStartupMessages({ library(readr); library(dplyr); library(tidyr); library(ggplot2) })
if (!requireNamespace("STRINGdb", quietly = TRUE)) stop("Install STRINGdb: BiocManager::install('STRINGdb')")
select <- dplyr::select; filter <- dplyr::filter; count <- dplyr::count; mutate <- dplyr::mutate; summarise <- dplyr::summarise; arrange <- dplyr::arrange; desc <- dplyr::desc   # STRINGdb loads plyr

QUERY <- Sys.getenv("QUERY", "Gba1"); STRING_SCORE <- as.integer(Sys.getenv("STRING_SCORE", "400")); STRING_VERSION <- Sys.getenv("STRING_VERSION", "12.0")
PER_PAGE <- as.integer(Sys.getenv("PER_PAGE", "20")); ONLY_SIG <- Sys.getenv("ONLY_SIG", "0") == "1"; LOW_EXPR <- 50
COMP <- c(IFNB_vs_Mock = "IFNB vs Mock", EBOV_vs_Mock = "EBOV vs Mock", IFNB_EBOV_vs_Mock = "IFNB+EBOV vs Mock"); CCOL <- c(`IFNB vs Mock` = "#2B6CB0", `EBOV vs Mock` = "#1B7837", `IFNB+EBOV vs Mock` = "#9BD08E")
INK <- "#1F2933"
PROJ <- Sys.getenv("PROJ", unset = file.path(
  Sys.getenv("HOME"),
  "Documents/aviti_rnaseq/2025/aviti",
  "20251009-100725_Maury_25161_WivfFacICOTHFAUJwSPeHBEFfCfNtCRpUSgVGUkz",
  "results/Maury_25161"))
RNA <- file.path(PROJ, "requant_v2_outputs"); OUT <- file.path(RNA, "figures_gba1_string"); dir.create(OUT, showWarnings = FALSE)
cache <- file.path(PROJ, "stringdb_cache"); dir.create(cache, showWarnings = FALSE); options(timeout = 600)
stars <- function(p) ifelse(is.na(p), "", ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", "ns"))))

message("connecting to STRING v", STRING_VERSION, " (mouse)")
sdb <- STRINGdb::STRINGdb$new(version = STRING_VERSION, species = 10090, score_threshold = 0L, network_type = "full", input_directory = cache)
prot <- sdb$get_proteins() |> as_tibble() |> select(STRING_id = protein_external_id, name = preferred_name)
qid <- NA_character_
for (q in unique(c(QUERY, "Gba1", "Gba"))) { h <- prot |> filter(tolower(name) == tolower(q)); if (nrow(h)) { qid <- h$STRING_id[1]; message("query '", q, "' -> ", qid); break } }
if (is.na(qid)) stop("Could not find ", QUERY, " in STRING mouse")
nb <- sdb$get_neighbors(qid); message(length(nb), " STRING neighbours (any score)")
ed <- sdb$get_interactions(c(qid, nb)) |> as_tibble() |> filter(from == qid | to == qid) |> mutate(STRING_id = ifelse(from == qid, to, from)) |> filter(STRING_id != qid) |>
  group_by(STRING_id) |> summarise(score = max(combined_score), .groups = "drop") |> left_join(prot, by = "STRING_id") |> filter(score >= STRING_SCORE) |> arrange(desc(score))
message(nrow(ed), " partners with combined score >= ", STRING_SCORE, " (", sum(ed$score >= 700), " high-confidence, >= 700)")

tabs <- lapply(names(COMP), function(nm) read_csv(file.path(RNA, paste0("DE_v2_", nm, "_annotated.csv")), show_col_types = FALSE) |> mutate(key = tolower(gene_name))); names(tabs) <- names(COMP)
ed <- ed |> mutate(key = tolower(name)) |> left_join(tabs[[1]] |> filter(!is.na(gene_name)) |> distinct(key, .keep_all = TRUE) |> select(key, gene_id), by = "key") |> mutate(in_DE = !is.na(gene_id))
write_csv(ed |> select(partner = name, STRING_id, combined_score = score, in_DE_tables = in_DE, gene_id), file.path(OUT, "gba1_string_partners.csv"))
message(sum(ed$in_DE), " of them are in the DE tables; not shown (removed by the pre-filter, or no matching symbol): ", paste(ed$name[!ed$in_DE], collapse = ", "))
ed <- ed |> filter(in_DE) |> mutate(pos = row_number())
if (ONLY_SIG) { sg <- unique(unlist(lapply(tabs, function(t) t$gene_id[!is.na(t$padj) & t$padj < 0.05]))); ed <- ed |> filter(gene_id %in% sg) |> mutate(pos = row_number()) }
if (!nrow(ed)) stop("No partners to plot")

D <- bind_rows(lapply(names(COMP), function(nm) tabs[[nm]] |> inner_join(ed |> select(gene_id, pos, score), by = "gene_id") |> transmute(gene = gene_name, gene_id, pos, score, baseMean, comparison = COMP[[nm]], lfc = log2FoldChange, se = lfcSE, padj)))
D <- D |> mutate(comparison = factor(comparison, levels = COMP), lab = stars(padj), vj = ifelse(lfc >= 0, 0, 1), ys = lfc + ifelse(lfc >= 0, 1, -1) * (se + 0.03), xlab = paste0(gene, ifelse(baseMean < LOW_EXPR, "†", "")), page = (pos - 1) %/% PER_PAGE + 1)
write_csv(D |> select(gene, gene_id, STRING_score = score, baseMean, comparison, log2FC = lfc, lfcSE = se, padj, page), file.path(OUT, "gba1_partners_three_comparisons.csv"))
print(as.data.frame(D |> transmute(gene, STRING = score, baseMean = round(baseMean, 1), comparison, log2FC = round(lfc, 2), padj = signif(padj, 2)) |> tidyr::pivot_wider(names_from = comparison, values_from = c(log2FC, padj))), row.names = FALSE)
npg <- max(D$page); dg <- position_dodge(width = 0.85)
for (pg in seq_len(npg)) {
  Dp <- D |> filter(page == pg); Dp$xlab <- factor(Dp$xlab, levels = unique(Dp$xlab[order(Dp$pos)]))
  p <- ggplot(Dp, aes(xlab, lfc, fill = comparison)) + geom_col(position = dg, width = 0.8) +
    geom_errorbar(aes(ymin = lfc - se, ymax = lfc + se), position = dg, width = 0.25, colour = INK, linewidth = 0.35) +
    geom_text(aes(y = ys, label = lab, vjust = vj), position = dg, colour = INK, size = 2.8) + geom_hline(yintercept = 0, colour = INK, linewidth = 0.4) +
    scale_fill_manual(values = CCOL, name = NULL) + scale_y_continuous(expand = expansion(mult = 0.12)) +
    labs(title = paste0("STRING partners of ", QUERY, " (mouse, score >= ", STRING_SCORE, "): IFNB, EBOV and IFNB+EBOV vs Mock", if (npg > 1) paste0("  [page ", pg, " of ", npg, "]") else ""),
         subtitle = "genes ordered by STRING combined score with Gba1 (strongest first); shrunk log2FC +/- 1 SE; stars = padj of each comparison (ns = not significant); dagger = baseMean < 50", x = NULL, y = "log2 fold-change vs Mock") +
    theme_bw(base_size = 11) + theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(), panel.grid.major.y = element_line(colour = "grey92", linewidth = 0.3), panel.border = element_rect(colour = "grey70"),
          axis.text = element_text(colour = INK), axis.text.x = element_text(face = "italic", angle = 45, hjust = 1), plot.title = element_text(colour = INK, face = "bold", size = 11), plot.subtitle = element_text(colour = "grey35", size = 8.5), legend.position = "top")
  w <- max(6, 0.9 * nlevels(Dp$xlab) + 2)
  ggsave(file.path(OUT, paste0("bars_gba1_partners_page", pg, ".png")), p, width = w, height = 5, dpi = 300); ggsave(file.path(OUT, paste0("bars_gba1_partners_page", pg, ".pdf")), p, width = w, height = 5)
}
message("\nSaved ", npg, " page(s) to: ", OUT)
