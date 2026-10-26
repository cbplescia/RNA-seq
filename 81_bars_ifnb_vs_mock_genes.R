#!/usr/bin/env Rscript
# Bar graph of the IFNB vs Mock differential expression of selected genes: shrunk log2FC from DE_v2_IFNB_vs_Mock_annotated.csv, error bars +/- 1 lfcSE,
# stars = padj (* < 0.05, ** < 0.01, *** < 0.001, ns), red = up, blue = down. Default genes (protein names -> mouse genes):
#   Gba1 = glucocerebrosidase (GCase) | Lrrk2 | Trip12 | c-Abl = Abl1 | LIMP-2 = Scarb2 | Saposin C = Psap (prosaposin gene; saposin C is one of its cleavage products,
#   so the RNA-seq value is for the whole Psap transcript)
# Other genes without editing the file:  GENES=Axl,Mertk OUT_NAME=my_genes Rscript 81_bars_ifnb_vs_mock_genes.R   (symbols or the aliases below; matched ignoring case)
# Writes bars_IFNB_vs_Mock_<OUT_NAME>.png / .pdf and a .csv to requant_v2_outputs/figures_gene_bars/
# Usage: Rscript 81_bars_ifnb_vs_mock_genes.R   (PROJ from 'source setup_session.sh')

suppressPackageStartupMessages({ library(readr); library(dplyr); library(ggplot2) })
select <- dplyr::select; filter <- dplyr::filter

QUERY <- c("Gba1", "Lrrk2", "Trip12", "c-Abl", "LIMP-2", "Saposin C")
if (nzchar(Sys.getenv("GENES"))) QUERY <- trimws(strsplit(Sys.getenv("GENES"), ",")[[1]])
OUT_NAME <- Sys.getenv("OUT_NAME", unset = "Gba1_Lrrk2_Trip12_Abl1_Scarb2_Psap")
ALIAS <- c(`gba1` = "Gba1|Gba", `c-abl` = "Abl1", `cabl` = "Abl1", `abl` = "Abl1", `limp-2` = "Scarb2", `limp2` = "Scarb2", `saposin c` = "Psap", `saposin` = "Psap",
           `prosaposin` = "Psap", `gcase` = "Gba1|Gba", `cert` = "Cert1|Col4a3bp")         # query (lower case) -> mouse symbol(s), "|" separates alternatives
INK <- "#1F2933"; RED <- "#C53030"; BLUE <- "#2B6CB0"
PROJ <- Sys.getenv("PROJ", unset = file.path(
  Sys.getenv("HOME"),
  "Documents/aviti_rnaseq/2025/aviti",
  "20251009-100725_Maury_25161_WivfFacICOTHFAUJwSPeHBEFfCfNtCRpUSgVGUkz",
  "results/Maury_25161"))
RNA <- file.path(PROJ, "requant_v2_outputs")
OUT <- file.path(RNA, "figures_gene_bars"); dir.create(OUT, showWarnings = FALSE)
stars <- function(p) ifelse(is.na(p), "", ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", "ns"))))

d <- read_csv(file.path(RNA, "DE_v2_IFNB_vs_Mock_annotated.csv"), show_col_types = FALSE) |> mutate(key = tolower(gene_name))
rows <- list()
for (q in QUERY) {
  syms <- tolower(strsplit(if (!is.na(ALIAS[tolower(q)])) ALIAS[[tolower(q)]] else q, "\\|")[[1]])
  h <- d |> filter(key %in% syms | tolower(sub("\\..*$", "", gene_id)) %in% syms) |> distinct(gene_id, .keep_all = TRUE)
  if (!nrow(h)) { message("Not found in the DE table (removed by the pre-filter, or a different symbol): ", q); next }
  rows[[length(rows) + 1]] <- h |> transmute(query = q, gene = gene_name, gene_id, baseMean, lfc = log2FoldChange, se = lfcSE, padj)
}
if (!length(rows)) stop("None of the genes were found")
D <- bind_rows(rows) |> mutate(lab = stars(padj), dir = ifelse(lfc >= 0, "up", "down"), ypos = ifelse(lfc >= 0, lfc + se, lfc - se),
                               name = ifelse(tolower(query) == tolower(gene), gene, sprintf("%s\n(%s)", gene, query)))
D$name <- factor(D$name, levels = unique(D$name))
write_csv(D |> select(query, gene, gene_id, baseMean, log2FC = lfc, lfcSE = se, padj), file.path(OUT, paste0("bars_IFNB_vs_Mock_", OUT_NAME, ".csv")))
print(as.data.frame(D |> transmute(query, gene, baseMean = round(baseMean, 1), log2FC = round(lfc, 2), padj = signif(padj, 3), stars = lab)), row.names = FALSE)
pad <- 0.05 * diff(range(c(D$lfc - D$se, D$lfc + D$se, 0)))
p <- ggplot(D, aes(name, lfc, fill = dir)) + geom_col(width = 0.65) +
  geom_errorbar(aes(ymin = lfc - se, ymax = lfc + se), width = 0.18, colour = INK, linewidth = 0.4) +
  geom_text(aes(y = ypos + ifelse(lfc >= 0, pad, -pad), label = lab), colour = INK, size = 4, vjust = 0.5) +
  geom_hline(yintercept = 0, colour = INK, linewidth = 0.4) + scale_fill_manual(values = c(up = RED, down = BLUE), guide = "none") + scale_y_continuous(expand = expansion(mult = 0.12)) +
  labs(title = "IFNB vs Mock: log2 fold-change of selected genes", subtitle = "shrunk log2FC +/- 1 SE; stars = padj (ns = not significant)", x = NULL, y = "log2 fold-change (IFNB vs Mock)") +
  theme_bw(base_size = 12) + theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(), panel.grid.major.y = element_line(colour = "grey92", linewidth = 0.3), panel.border = element_rect(colour = "grey70"),
        axis.text = element_text(colour = INK), axis.title = element_text(colour = INK), plot.title = element_text(colour = INK, face = "bold"), plot.subtitle = element_text(colour = "grey35", size = 9))
w <- max(4.5, 0.9 * nrow(D) + 2)
ggsave(file.path(OUT, paste0("bars_IFNB_vs_Mock_", OUT_NAME, ".png")), p, width = w, height = 4.8, dpi = 300)
ggsave(file.path(OUT, paste0("bars_IFNB_vs_Mock_", OUT_NAME, ".pdf")), p, width = w, height = 4.8)
message("Saved to: ", OUT)
