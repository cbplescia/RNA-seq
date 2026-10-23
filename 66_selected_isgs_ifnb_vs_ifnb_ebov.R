#!/usr/bin/env Rscript
# Selected genes (edit GENES): for each gene two bars: log2FC of IFNB vs Mock (blue) and IFNB+EBOV vs Mock (light green), shrunk, error bars +/- 1 lfcSE;
# genes are ordered by their IFNB log2FC. A bracket joins the two bars and gives the stat of IFNB+EBOV vs IFNB: DESeq2 Wald test on the difference between the
# two log2FCs (contrast IFNB_EBOV vs IFNB_Mock, BH padj across all genes), i.e. whether EBOV significantly changes the IFNB induction of that gene
# (* < 0.05, ** < 0.01, *** < 0.001, ns). Small grey stars above each bar = padj vs Mock.
# Writes selected_ISGs_IFNB_vs_IFNB_EBOV.png / .pdf and .csv to requant_v2_outputs/figures_gene_bars/
# Usage: Rscript 66_selected_isgs_ifnb_vs_ifnb_ebov.R   (PROJ from 'source setup_session.sh')

suppressPackageStartupMessages({ library(readr); library(dplyr); library(ggplot2); library(DESeq2) })
select <- dplyr::select; filter <- dplyr::filter

GENES <- c("Zbp1", "Isg15", "Xaf1", "Irf7", "Mx2", "Srms", "Ifi208", "Mx1", "Usp18", "Gm43302", "Cfb")     # matched ignoring case
# a trailing * matches a whole family:  GENES='Gbp*' ... (all genes whose symbol starts with Gbp)
# any other gene(s) without editing the file:  GENES=Axl OUT_NAME=Axl_IFNB_vs_IFNB_EBOV Rscript 66_selected_isgs_ifnb_vs_ifnb_ebov.R   (comma-separate several genes)
if (nzchar(Sys.getenv("GENES"))) GENES <- trimws(strsplit(Sys.getenv("GENES"), ",")[[1]])
OUT_NAME <- Sys.getenv("OUT_NAME", unset = "selected_ISGs_IFNB_vs_IFNB_EBOV")
INK <- "#1F2933"; C_B <- "#2B6CB0"; C_BE <- "#9BD08E"
PROJ <- Sys.getenv("PROJ", unset = file.path(
  Sys.getenv("HOME"),
  "Documents/aviti_rnaseq/2025/aviti",
  "20251009-100725_Maury_25161_WivfFacICOTHFAUJwSPeHBEFfCfNtCRpUSgVGUkz",
  "results/Maury_25161"))
RNA <- file.path(PROJ, "requant_v2_outputs")
OUT <- file.path(RNA, "figures_gene_bars"); dir.create(OUT, showWarnings = FALSE)
stars <- function(p) ifelse(is.na(p), "", ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", "ns"))))

allg <- read_csv(file.path(RNA, "DE_v2_IFNB_vs_Mock_annotated.csv"), show_col_types = FALSE) |> mutate(key = tolower(gene_name))
wild <- tolower(sub("\\*$", "", grep("\\*$", GENES, value = TRUE)))                   # "Gbp*" = every gene whose symbol starts with Gbp
GENES <- GENES[!grepl("\\*$", GENES)]
want <- tolower(GENES)                                    # gene symbols or Ensembl gene IDs (ENSMUSG...) both work
is_wild <- if (length(wild)) Reduce(`|`, lapply(wild, function(p) startsWith(allg$key, p))) else rep(FALSE, nrow(allg))
ann <- allg |> filter(key %in% want | tolower(gene_id) %in% want | tolower(sub("\\..*$", "", gene_id)) %in% want | is_wild) |> distinct(gene_id, .keep_all = TRUE)
if (length(wild)) message("Wildcard matches: ", paste(ann$gene_name[startsWith(tolower(ann$gene_name), wild[1])], collapse = ", "))
miss <- GENES[!vapply(want, function(w) any(allg$key == w | tolower(allg$gene_id) == w | tolower(sub("\\..*$", "", allg$gene_id)) == w), logical(1))]
if (length(miss)) {
  message("Not found in the DE table: ", paste(miss, collapse = ", "))
  for (m in miss) {                                       # similar symbols, so a renamed gene can be spotted
    cand <- allg |> filter(startsWith(key, substr(tolower(m), 1, min(4, nchar(m))))) |> arrange(desc(baseMean)) |> slice_head(n = 12)
    if (nrow(cand)) { message("  symbols starting like '", m, "':"); print(as.data.frame(cand |> select(gene_name, gene_id, baseMean)), row.names = FALSE) }
  }
}
if (!nrow(ann)) stop("None of the genes were found")
top <- ann |> arrange(desc(log2FoldChange)) |> transmute(gene_id, gene_name, lfc_B = log2FoldChange, se_B = lfcSE, padj_B = padj)
be <- read_csv(file.path(RNA, "DE_v2_IFNB_EBOV_vs_Mock_annotated.csv"), show_col_types = FALSE) |> transmute(gene_id, lfc_BE = log2FoldChange, se_BE = lfcSE, padj_BE = padj)
dds <- readRDS(file.path(RNA, "dds_v2_decoyAware_hostNormalised.rds"))
df <- as.data.frame(results(dds, contrast = c("group", "IFNB_EBOV", "IFNB_Mock"))) |> tibble::rownames_to_column("gene_id") |> transmute(gene_id, diff_lfc = log2FoldChange, diff_padj = padj)
T <- top |> left_join(be, by = "gene_id") |> left_join(df, by = "gene_id") |> mutate(gene_name = factor(gene_name, levels = gene_name))
write_csv(T |> transmute(gene_name, gene_id, log2FC_IFNB_vs_Mock = lfc_B, padj_IFNB_vs_Mock = padj_B, log2FC_IFNB_EBOV_vs_Mock = lfc_BE, padj_IFNB_EBOV_vs_Mock = padj_BE,
                         log2FC_diff_IFNB_EBOV_vs_IFNB_unshrunk = diff_lfc, padj_IFNB_EBOV_vs_IFNB = diff_padj), file.path(OUT, paste0(OUT_NAME, ".csv")))
print(as.data.frame(T |> transmute(gene = gene_name, IFNB = round(lfc_B, 2), IFNB_EBOV = round(lfc_BE, 2), padj_diff = signif(diff_padj, 3), stars = stars(diff_padj))), row.names = FALSE)

L <- bind_rows(T |> transmute(gene_name, cond = "IFNB vs Mock", lfc = lfc_B, se = se_B, padj = padj_B),
               T |> transmute(gene_name, cond = "IFNB+EBOV vs Mock", lfc = lfc_BE, se = se_BE, padj = padj_BE)) |>
  mutate(cond = factor(cond, levels = c("IFNB vs Mock", "IFNB+EBOV vs Mock")), top = lfc + se)
pad <- 0.03 * max(L$top, na.rm = TRUE)
bk <- L |> group_by(gene_name) |> summarise(y = max(top, na.rm = TRUE) + 2.2 * pad, .groups = "drop") |> left_join(T |> select(gene_name, diff_padj), by = "gene_name") |>
  mutate(x = as.numeric(gene_name), lab = stars(diff_padj))
dg <- position_dodge(width = 0.8)
p <- ggplot(L, aes(gene_name, lfc, fill = cond)) + geom_col(position = dg, width = 0.8) +
  geom_errorbar(aes(ymin = lfc - se, ymax = lfc + se), position = dg, width = 0.22, colour = INK, linewidth = 0.35) +
  geom_text(aes(y = top + 0.4 * pad, label = stars(padj)), position = dg, colour = "grey40", size = 2.4, vjust = 0) +
  geom_segment(data = bk, aes(x = x - 0.2, xend = x + 0.2, y = y, yend = y), inherit.aes = FALSE, colour = INK, linewidth = 0.4) +
  geom_segment(data = bk, aes(x = x - 0.2, xend = x - 0.2, y = y, yend = y - pad), inherit.aes = FALSE, colour = INK, linewidth = 0.4) +
  geom_segment(data = bk, aes(x = x + 0.2, xend = x + 0.2, y = y, yend = y - pad), inherit.aes = FALSE, colour = INK, linewidth = 0.4) +
  geom_text(data = bk, aes(x = x, y = y + 0.4 * pad, label = lab), inherit.aes = FALSE, colour = INK, size = 4, vjust = 0, fontface = "bold") +
  geom_hline(yintercept = 0, colour = INK, linewidth = 0.4) + scale_fill_manual(values = c(`IFNB vs Mock` = C_B, `IFNB+EBOV vs Mock` = C_BE), name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.12))) +
  labs(title = if (nrow(T) == 1) sprintf("%s: IFNB vs IFNB+EBOV", T$gene_name[1]) else "Selected IFNB-induced genes: IFNB vs IFNB+EBOV",
       subtitle = "bars = shrunk log2FC vs Mock +/- 1 SE (small grey stars = padj vs Mock); bracket = IFNB+EBOV vs IFNB (does EBOV change the induction?)", x = NULL, y = "log2 fold-change vs Mock") +
  theme_bw(base_size = 12) + theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(), panel.grid.major.y = element_line(colour = "grey92", linewidth = 0.3),
        panel.border = element_rect(colour = "grey70"), axis.text = element_text(colour = INK), axis.text.x = element_text(angle = 45, hjust = 1, face = "italic"),
        axis.title = element_text(colour = INK), plot.title = element_text(colour = INK, face = "bold"), plot.subtitle = element_text(colour = "grey35", size = 8.5), legend.position = "top")
W <- max(4.5, 0.8 * nrow(T) + 2)
ggsave(file.path(OUT, paste0(OUT_NAME, ".png")), p, width = W, height = 5.5, dpi = 300)
ggsave(file.path(OUT, paste0(OUT_NAME, ".pdf")), p, width = W, height = 5.5)
message("Saved to: ", OUT)
