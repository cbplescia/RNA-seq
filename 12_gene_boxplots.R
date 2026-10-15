#!/usr/bin/env Rscript
# VST expression boxplots (Mock vs IFNB, every sample as a point) for chosen genes.
# Edit GENES below to plot others (case-insensitive mouse symbols).
# Each gene gets its own panel with its own y-axis, because VST levels differ a lot
# between genes; the panel title also gives the IFNB vs Mock log2FC and padj.
# Needs 07_add_gene_symbols.R to have been run.

suppressPackageStartupMessages({
  library(DESeq2); library(readr); library(dplyr); library(ggplot2)
})

# Mouse symbols. Human-style names with no mouse gene of that name (BTN3A3, IFI6, MICB
# are the likely ones) are reported as not found, with similarly named genes listed.
GENES <- c("Akt3", "Odc1", "Flt1", "Pfkfb3", "Map3k5", "Gbp2", "Dcp1a", "Ccdc92", "Ifi6", "Btn3a3")
NCOL <- 5      # panels are drawn left to right, row by row, in the order of GENES above

PROJ <- Sys.getenv("PROJ", unset = file.path(
  Sys.getenv("HOME"),
  "Documents/aviti_rnaseq/2025/aviti",
  "20251009-100725_Maury_25161_WivfFacICOTHFAUJwSPeHBEFfCfNtCRpUSgVGUkz",
  "results/Maury_25161"))
OUT <- file.path(PROJ, "requant_v2_outputs")
FIG <- file.path(OUT, "figures_genes")
dir.create(FIG, showWarnings = FALSE)

INK <- "#1F2933"
COL <- c(Mock = "#8A8A8A", IFNB = "#2B6CB0")

de <- read_csv(file.path(OUT, "DE_v2_IFNB_vs_Mock_annotated.csv"), show_col_types = FALSE)
hit <- de |> filter(toupper(gene_name) %in% toupper(GENES))
miss <- setdiff(toupper(GENES), toupper(hit$gene_name))
if (length(miss)) {
  message("\nNot found in the DE table: ", paste(miss, collapse = ", "),
          "\n  (either no mouse gene has that symbol, or it was too lowly expressed and removed by the filter in 04)")
  for (m in miss) {
    near <- sort(unique(de$gene_name[grepl(paste0("^", substr(m, 1, 4)), toupper(de$gene_name))]))
    if (length(near)) message("  genes starting '", substr(m, 1, 4), "' that ARE present: ",
                              paste(head(near, 15), collapse = ", "))
  }
  message("")
}
if (!nrow(hit)) stop("None of the requested genes were found.")

vsd  <- readRDS(file.path(OUT, "vsd_v2.rds"))
vmat <- assay(vsd)
samp <- tibble(sample_id = colnames(vmat), group = as.character(vsd$group)) |>
  filter(group %in% c("Mock_Mock", "IFNB_Mock")) |>
  mutate(treatment = factor(ifelse(group == "Mock_Mock", "Mock", "IFNB"), levels = c("Mock", "IFNB")))

fmt_p <- function(p) ifelse(is.na(p), "padj n/a", ifelse(p < 0.001, "padj < 0.001", sprintf("padj = %.3g", p)))
hit <- hit |>
  mutate(panel = sprintf("%s\nlog2FC %+.2f, %s", gene_name, log2FoldChange, fmt_p(padj)))
hit$panel <- factor(hit$panel, levels = hit$panel[order(match(toupper(hit$gene_name), toupper(GENES)))])

ids <- hit$gene_id[hit$gene_id %in% rownames(vmat)]
long <- as.data.frame(t(vmat[ids, samp$sample_id, drop = FALSE])) |>
  tibble::rownames_to_column("sample_id") |>
  tidyr::pivot_longer(-sample_id, names_to = "gene_id", values_to = "vst") |>
  left_join(samp, by = "sample_id") |>
  left_join(hit |> select(gene_id, gene_name, panel), by = "gene_id")
write_csv(long |> select(-panel), file.path(FIG, "selected_genes_vst_values.csv"))

p <- ggplot(long, aes(treatment, vst, fill = treatment, colour = treatment)) +
  geom_boxplot(alpha = 0.25, outlier.shape = NA, width = 0.55, linewidth = 0.5) +
  geom_jitter(width = 0.1, height = 0, size = 2, alpha = 0.9) +
  facet_wrap(~panel, ncol = NCOL, scales = "free_y") +
  scale_fill_manual(values = COL, guide = "none") +
  scale_colour_manual(values = COL, guide = "none") +
  labs(title = "Mock vs IFNB", subtitle = "VST-normalised expression; each point is one sample (n = 3 per group)",
       x = NULL, y = "VST expression (log2 scale)") +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(colour = "grey92", linewidth = 0.3),
        panel.border = element_rect(colour = "grey70"),
        strip.background = element_rect(fill = "grey95", colour = "grey70"),
        strip.text = element_text(colour = INK, face = "bold"),
        axis.text = element_text(colour = INK), axis.title = element_text(colour = INK),
        plot.title = element_text(colour = INK, face = "bold"),
        plot.subtitle = element_text(colour = "grey35"))
n <- length(unique(long$gene_id))
w <- 2.5 * min(n, NCOL) + 0.8
h <- 3.3 * ceiling(n / NCOL) + 1.2
ggsave(file.path(FIG, "selected_genes_vst_boxplot_Mock_vs_IFNB_ordered.png"), p, width = w, height = h, dpi = 300)
ggsave(file.path(FIG, "selected_genes_vst_boxplot_Mock_vs_IFNB_ordered.pdf"), p, width = w, height = h)
message("Plotted: ", paste(unique(long$gene_name), collapse = ", "), "\nSaved in: ", FIG)
