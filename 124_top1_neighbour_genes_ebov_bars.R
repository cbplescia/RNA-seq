#!/usr/bin/env Rscript
# Genes that sit next to Top1 in the mouse genome and the genes that share its promoter, with a bar graph of their EBOV vs Mock differential expression.
#   shared promoter = a gene on the OPPOSITE strand whose transcription start site lies within PROM_KB (default 2 kb) of the Top1 TSS (a head-to-head / bidirectional promoter), plus any gene whose TSS is on the same side within that distance
#   neighbours      = every other gene whose body lies within WINDOW_KB (default 250 kb) of the Top1 gene; split into <= NEAR_KB (default 50 kb) and farther. Distances are gaps between gene bodies in kb (0 = overlapping).
# Gene positions come from ref/Mus_musculus.GRCm39.109.gtf; expression from DE_v2_EBOV_vs_Mock_annotated.csv (shrunk log2FC +/- 1 SE, stars = padj, dagger = baseMean < 50). Protein-coding genes only unless ALL=1. Genes not in the DE table (removed by the pre-filter in 04) are listed as "Not in DE table".
# Genes in each panel are ordered by distance from Top1 (closest first); Top1 itself is the first panel. Writes bars_top1_neighbours_EBOV_vs_Mock.png / .pdf and top1_neighbours_EBOV_vs_Mock.csv to requant_v2_outputs/figures_gene_bars/.
# Env: GENE (default Top1; any gene symbol works, e.g. GENE=Top2a), WINDOW_KB, NEAR_KB, PROM_KB, ALL=1, COMP (default EBOV_vs_Mock; any DE_v2_<COMP>_annotated.csv).   Usage: Rscript 124_top1_neighbour_genes_ebov_bars.R   (PROJ from 'source setup_session.sh')
# This finds genomic neighbours and shared promoters. Genes that merely share transcription-factor motifs are a different question (needs the motif scan, 106 / 108).

suppressPackageStartupMessages({ library(readr); library(dplyr); library(ggplot2) })
select <- dplyr::select; filter <- dplyr::filter

GENE <- Sys.getenv("GENE", "Top1"); WINDOW_KB <- as.numeric(Sys.getenv("WINDOW_KB", "250")); NEAR_KB <- as.numeric(Sys.getenv("NEAR_KB", "50")); PROM_KB <- as.numeric(Sys.getenv("PROM_KB", "2")); ALL <- Sys.getenv("ALL", "0") == "1"
COMP <- Sys.getenv("COMP", "EBOV_vs_Mock"); COMP_LAB <- gsub("_", " ", gsub("_vs_", " vs ", COMP)); LOW_EXPR <- 50; INK <- "#1F2933"
FILL <- c(`EBOV_vs_Mock` = "#1B7837", `IFNB_vs_Mock` = "#2B6CB0", `IFNB_EBOV_vs_Mock` = "#9BD08E")[COMP]; if (is.na(FILL)) FILL <- "#1B7837"
PROJ <- Sys.getenv("PROJ", unset = file.path(
  Sys.getenv("HOME"),
  "Documents/aviti_rnaseq/2025/aviti",
  "20251009-100725_Maury_25161_WivfFacICOTHFAUJwSPeHBEFfCfNtCRpUSgVGUkz",
  "results/Maury_25161"))
RNA <- file.path(PROJ, "requant_v2_outputs"); OUT <- file.path(RNA, "figures_gene_bars"); dir.create(OUT, showWarnings = FALSE)
GTF <- file.path(PROJ, "ref", "Mus_musculus.GRCm39.109.gtf")
stars <- function(p) ifelse(is.na(p), "", ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", "ns"))))

message("reading gene records from ", GTF)
g <- read_tsv(GTF, comment = "#", col_names = c("chr", "src", "feature", "start", "end", "score", "strand", "frame", "attr"), col_types = "ccciiccci", progress = FALSE) |> filter(feature == "gene")
attr <- function(x, k) sub(paste0('.*', k, ' "([^"]*)".*'), "\\1", ifelse(grepl(paste0(k, ' "'), x), x, NA_character_))
g <- g |> mutate(gene_id = attr(attr, "gene_id"), gene_name = attr(attr, "gene_name"), biotype = attr(attr, "gene_biotype"), tss = ifelse(strand == "+", start, end)) |> select(-attr, -score, -frame, -src, -feature)
g$gene_name <- ifelse(is.na(g$gene_name), g$gene_id, g$gene_name)
tg <- g |> filter(tolower(gene_name) == tolower(GENE)); if (!nrow(tg)) stop("'", GENE, "' not found in the GTF")
tg <- tg[1, ]; message(GENE, ": chr", tg$chr, ":", tg$start, "-", tg$end, " (", tg$strand, " strand), TSS ", tg$tss)
nb <- g |> filter(chr == tg$chr, gene_id != tg$gene_id, end >= tg$start - WINDOW_KB * 1000, start <= tg$end + WINDOW_KB * 1000)
if (!ALL) nb <- nb |> filter(biotype == "protein_coding")
nb <- nb |> mutate(gap_kb = pmax(0, pmax(start - tg$end, tg$start - end)) / 1000, tss_dist_kb = abs(tss - tg$tss) / 1000, side = ifelse(end < tg$start, "upstream (lower coordinate)", ifelse(start > tg$end, "downstream (higher coordinate)", "overlaps Top1"))) |>
  mutate(class = ifelse(tss_dist_kb <= PROM_KB, paste0("Shared promoter (TSS within ", PROM_KB, " kb)"), ifelse(gap_kb <= NEAR_KB, paste0("Neighbours within ", NEAR_KB, " kb"), paste0("Neighbours ", NEAR_KB, "-", WINDOW_KB, " kb")))) |>
  mutate(head_to_head = strand != tg$strand & tss_dist_kb <= PROM_KB) |> arrange(gap_kb, tss_dist_kb)
message(nrow(nb), " genes in the window (", sum(nb$class == paste0("Shared promoter (TSS within ", PROM_KB, " kb)")), " share the promoter)")
print(as.data.frame(nb |> transmute(gene = gene_name, biotype, strand, side, gap_kb = round(gap_kb, 1), tss_dist_kb = round(tss_dist_kb, 1), class)), row.names = FALSE)

de <- read_csv(file.path(RNA, paste0("DE_v2_", COMP, "_annotated.csv")), show_col_types = FALSE) |> mutate(gid = sub("\\..*$", "", gene_id))
want <- bind_rows(tg |> mutate(gap_kb = 0, tss_dist_kb = 0, class = paste0(GENE, " (query)"), head_to_head = FALSE, side = "query"), nb) |> mutate(order = row_number())
D <- want |> left_join(de |> select(gid, baseMean, lfc = log2FoldChange, se = lfcSE, padj), by = c("gene_id" = "gid"))
nd <- D |> filter(is.na(lfc)); if (nrow(nd)) message("Not in DE table (removed by the pre-filter, or not expressed): ", paste(nd$gene_name, collapse = ", "))
D <- D |> filter(!is.na(lfc)) |> mutate(class = factor(class, levels = unique(want$class[order(want$class != paste0(GENE, " (query)"), !grepl("^Shared", want$class), grepl("-", want$class))])),
                  lab = stars(padj), vj = ifelse(lfc >= 0, 0, 1), ys = lfc + ifelse(lfc >= 0, 1, -1) * (se + 0.03), xlab = paste0(gene_name, ifelse(strand == "+", " (+)", " (-)"), ifelse(baseMean < LOW_EXPR, "†", "")))
D <- D |> arrange(class, order); D$xlab <- factor(D$xlab, levels = unique(D$xlab))
write_csv(D |> select(class, gene = gene_name, gene_id, chr, start, end, strand, biotype, gap_kb, tss_dist_kb, head_to_head, baseMean, log2FC = lfc, lfcSE = se, padj), file.path(OUT, paste0("top1_neighbours_", COMP, ".csv")))
print(as.data.frame(D |> transmute(class, gene = gene_name, baseMean = round(baseMean, 1), log2FC = round(lfc, 2), padj = signif(padj, 2))), row.names = FALSE)

p <- ggplot(D, aes(xlab, lfc)) + geom_col(fill = FILL, width = 0.7) + geom_errorbar(aes(ymin = lfc - se, ymax = lfc + se), width = 0.25, colour = INK, linewidth = 0.35) +
  geom_text(aes(y = ys, label = lab, vjust = vj), colour = INK, size = 3) + geom_hline(yintercept = 0, colour = INK, linewidth = 0.4) +
  facet_grid(. ~ class, scales = "free_x", space = "free_x", labeller = labeller(class = label_wrap_gen(22))) + scale_y_continuous(expand = expansion(mult = 0.12)) +
  labs(title = sprintf("Genes near / sharing the promoter of %s (chr%s): %s", GENE, tg$chr, COMP_LAB), subtitle = paste0("within ", WINDOW_KB, " kb of the gene; ordered by distance from ", GENE, " (closest first); (+)/(-) = strand; shrunk log2FC +/- 1 SE; stars = padj (ns = not significant); dagger = baseMean < 50"),
       x = NULL, y = paste0("log2 fold-change (", COMP_LAB, ")")) +
  theme_bw(base_size = 11) + theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(), panel.grid.major.y = element_line(colour = "grey92", linewidth = 0.3), panel.border = element_rect(colour = "grey70"),
        axis.text = element_text(colour = INK), axis.text.x = element_text(face = "italic", angle = 55, hjust = 1), plot.title = element_text(colour = INK, face = "bold", size = 11), plot.subtitle = element_text(colour = "grey35", size = 8),
        strip.background = element_rect(fill = "grey93", colour = "grey70"), strip.text = element_text(colour = INK, face = "bold", size = 8))
w <- max(7, 0.55 * nrow(D) + 3)
ggsave(file.path(OUT, paste0("bars_top1_neighbours_", COMP, ".png")), p, width = w, height = 5.6, dpi = 300); ggsave(file.path(OUT, paste0("bars_top1_neighbours_", COMP, ".pdf")), p, width = w, height = 5.6)
message("\nSaved to: ", OUT)
