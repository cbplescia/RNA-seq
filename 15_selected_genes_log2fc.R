#!/usr/bin/env Rscript
# One PNG: IFNB vs Mock log2 fold-change bar plot for a chosen gene list, ordered
# in the order given in REQ, left to right. Shrunk log2FC +/- lfcSE; stars = padj.
# Genes not present in the data (no mouse gene with that symbol, or removed by the
# expression filter in 04) keep their place on the axis, marked "not in data", so
# the requested list is visibly complete.
# Needs 07_add_gene_symbols.R to have been run.

suppressPackageStartupMessages({ library(readr); library(dplyr); library(ggplot2) })

# lookup = mouse symbol searched for; label = how it is shown on the axis
REQ <- tibble::tribble(
  ~lookup,   ~label,
  "Akt3",    "Akt3",
  "Odc1",    "Odc1",
  "Flt1",    "Flt1",
  "Pfkfb3",  "Pfkfb3",
  "Map3k5",  "Map3k5",
  "Gbp2",    "Gbp2",
  "Dcp1a",   "Dcp1a",
  "Ccdc92",  "Ccdc92",
  "Ifi6",    "Ifi6",
  "Btn3a3",  "Btn3a3")
# Order on the x-axis is the order above (left to right). Edit this list to change it.

LFC_CUT <- 1; PADJ_CUT <- 0.05
INK <- "#1F2933"
COL <- c(Up = "#C53030", Down = "#2B6CB0", `Not significant` = "#B8B8B8")

PROJ <- Sys.getenv("PROJ", unset = file.path(
  Sys.getenv("HOME"),
  "Documents/aviti_rnaseq/2025/aviti",
  "20251009-100725_Maury_25161_WivfFacICOTHFAUJwSPeHBEFfCfNtCRpUSgVGUkz",
  "results/Maury_25161"))
OUT <- file.path(PROJ, "requant_v2_outputs")
FIG <- file.path(OUT, "figures_genes")
dir.create(FIG, showWarnings = FALSE)

de <- read_csv(file.path(OUT, "DE_v2_IFNB_vs_Mock_annotated.csv"), show_col_types = FALSE)
if (!"lfcSE" %in% names(de)) de$lfcSE <- NA_real_
de <- de |> mutate(key = toupper(gene_name)) |>
  distinct(key, .keep_all = TRUE) |>
  select(key, gene_id, log2FoldChange, lfcSE, padj)

d <- REQ |> mutate(key = toupper(lookup)) |> left_join(de, by = "key") |>
  mutate(found = !is.na(log2FoldChange),
         direction = case_when(!found ~ NA_character_,
                               !is.na(padj) & padj < PADJ_CUT & log2FoldChange >  LFC_CUT ~ "Up",
                               !is.na(padj) & padj < PADJ_CUT & log2FoldChange < -LFC_CUT ~ "Down",
                               TRUE ~ "Not significant"),
         stars = case_when(!found ~ "", is.na(padj) ~ "", padj < 0.001 ~ "***",
                           padj < 0.01 ~ "**", padj < 0.05 ~ "*", TRUE ~ "ns"))
d$label <- factor(d$label, levels = REQ$label)           # keep the order given in REQ
d$direction <- factor(d$direction, levels = names(COL))

miss <- d$label[!d$found]
if (length(miss)) message("Not in the data (shown as 'not in data'): ", paste(miss, collapse = ", "))

df <- d |> filter(found) |> mutate(
  top = ifelse(log2FoldChange >= 0, log2FoldChange + coalesce(lfcSE, 0), log2FoldChange - coalesce(lfcSE, 0)))
ymax <- max(c(df$log2FoldChange, 0)); ymin <- min(c(df$log2FoldChange, 0))
pad <- 0.06 * (ymax - ymin + 1)

p <- ggplot(d, aes(label, log2FoldChange)) +
  geom_hline(yintercept = 0, colour = "grey40", linewidth = 0.4) +
  geom_col(aes(fill = direction), width = 0.65, na.rm = TRUE) +
  geom_errorbar(data = df, aes(ymin = log2FoldChange - lfcSE, ymax = log2FoldChange + lfcSE),
                width = 0.2, colour = INK, linewidth = 0.4, na.rm = TRUE) +
  geom_text(data = df, aes(y = top, label = stars, vjust = ifelse(log2FoldChange >= 0, -0.4, 1.3)),
            colour = INK, size = 3.4) +
  geom_text(data = d |> filter(!found), aes(y = 0, label = "not in data"),
            angle = 90, hjust = -0.1, colour = "grey50", size = 3, fontface = "italic") +
  scale_x_discrete(drop = FALSE) +
  scale_fill_manual(values = COL, drop = FALSE, name = NULL, na.translate = FALSE) +
  scale_y_continuous(expand = expansion(add = c(pad, pad))) +
  labs(title = "IFNB vs Mock: selected genes",
       subtitle = sprintf("shrunk log2 fold-change ± lfcSE; coloured when padj<%g & |log2FC|>%g; * padj<0.05, ** <0.01, *** <0.001",
                          PADJ_CUT, LFC_CUT),
       x = NULL, y = "log2 fold-change (IFNB vs Mock)") +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(colour = "grey92", linewidth = 0.3),
        panel.border = element_rect(colour = "grey70"),
        axis.text = element_text(colour = INK), axis.text.x = element_text(angle = 45, hjust = 1, face = "italic"),
        axis.title = element_text(colour = INK),
        plot.title = element_text(colour = INK, face = "bold"),
        plot.subtitle = element_text(colour = "grey35", size = 9),
        legend.position = "bottom")
ggsave(file.path(FIG, "selected_genes_log2FC_IFNB_vs_Mock_ordered.png"), p, width = 9, height = 5.5, dpi = 300)
write_csv(d |> select(label, gene_id, log2FoldChange, lfcSE, padj, direction, found),
          file.path(FIG, "selected_genes_log2FC_IFNB_vs_Mock_ordered.csv"))
message("Saved: ", file.path(FIG, "selected_genes_log2FC_IFNB_vs_Mock_ordered.png"))
