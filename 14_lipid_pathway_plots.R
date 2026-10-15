#!/usr/bin/env Rscript
# One bar-graph PNG per lipid pathway / substrate group, showing the IFNB vs Mock
# log2 fold-change of every lipid-metabolism gene in that group.
#
# Grouping: each gene is put in the FIRST group below whose GO terms it is annotated
# to (using lipid_metabolism_gene_GO_terms.csv from 13). Order runs from specific
# to general, so e.g. a cholesterol-synthesis gene lands in "Cholesterol & sterol
# synthesis" and not in the broader steroid group. Groups with fewer than MIN_GENES
# genes are folded into "Other lipid metabolism", and groups with more than
# MAX_PER_PLOT genes are split over several PNGs.
#
# To move a gene yourself, make a file lipid_group_overrides.csv in the project
# folder with columns gene_name,group (group = one of the names in GROUPS below).
# Writes lipid_gene_groups.csv so you can review every assignment.

suppressPackageStartupMessages({ library(readr); library(dplyr); library(ggplot2) })

MIN_GENES <- 4; MAX_PER_PLOT <- 40
LFC_CUT <- 1; PADJ_CUT <- 0.05
INK <- "#1F2933"
COL <- c(Up = "#C53030", Down = "#2B6CB0", `Not significant` = "#B8B8B8")

GROUPS <- list(
  "Cholesterol & sterol synthesis"            = c("GO:0006695", "GO:0016126"),
  "Bile acid metabolism"                      = "GO:0008206",
  "Other cholesterol, sterol & steroid metabolism" = "GO:0008202",
  "Fatty acid oxidation"                      = c("GO:0019395", "GO:0006635"),
  "Fatty acid synthesis, elongation & desaturation" = c("GO:0006633", "GO:0030497", "GO:0006636"),
  "Eicosanoid & prostaglandin metabolism"     = c("GO:0006690", "GO:0019369"),
  "Sphingolipid, ceramide & glycolipid metabolism" = c("GO:0006665", "GO:0006664"),
  "Phosphoinositide metabolism"               = "GO:0046488",
  "Phospholipid breakdown (phospholipases)"   = "GO:0009395",
  "Phospholipid & glycerophospholipid metabolism" = c("GO:0006644", "GO:0006650"),
  "Triglyceride, glycerolipid & neutral lipid metabolism" = c("GO:0006641", "GO:0006638", "GO:0046486"),
  "Other fatty acid metabolism"               = "GO:0006631",
  "Other lipid catabolism"                    = "GO:0016042"
)
OTHER <- "Other lipid metabolism"

PROJ <- Sys.getenv("PROJ", unset = file.path(
  Sys.getenv("HOME"),
  "Documents/aviti_rnaseq/2025/aviti",
  "20251009-100725_Maury_25161_WivfFacICOTHFAUJwSPeHBEFfCfNtCRpUSgVGUkz",
  "results/Maury_25161"))
OUT <- file.path(PROJ, "requant_v2_outputs")
FIG <- file.path(OUT, "figures_lipid_pathways")
dir.create(FIG, showWarnings = FALSE)

lip <- read_csv(file.path(OUT, "DE_v2_IFNB_vs_Mock_protein_coding_lipid_metabolism.csv"), show_col_types = FALSE)
trm <- read_csv(file.path(OUT, "lipid_metabolism_gene_GO_terms.csv"), show_col_types = FALSE)
terms_by_gene <- split(trm$go_id, trm$gene_id)

assign_group <- function(gid) {
  t <- terms_by_gene[[gid]]
  for (g in names(GROUPS)) if (any(GROUPS[[g]] %in% t)) return(g)
  OTHER
}
lip$group <- vapply(lip$gene_id, assign_group, character(1))

ov <- file.path(PROJ, "lipid_group_overrides.csv")
if (file.exists(ov)) {
  o <- read_csv(ov, show_col_types = FALSE)
  lip$group[match(toupper(o$gene_name), toupper(lip$gene_name)) |> na.omit()] <-
    o$group[!is.na(match(toupper(o$gene_name), toupper(lip$gene_name)))]
  message("Applied ", nrow(o), " overrides from ", ov)
}

sizes <- table(lip$group)
small <- setdiff(names(sizes)[sizes < MIN_GENES], OTHER)
lip$group[lip$group %in% small] <- OTHER
if (length(small)) message("Folded into '", OTHER, "' (fewer than ", MIN_GENES, " genes): ", paste(small, collapse = "; "))

lip <- lip |> mutate(
  direction = case_when(!is.na(padj) & padj < PADJ_CUT & log2FoldChange >  LFC_CUT ~ "Up",
                        !is.na(padj) & padj < PADJ_CUT & log2FoldChange < -LFC_CUT ~ "Down",
                        TRUE ~ "Not significant"))
write_csv(lip |> select(gene_id, gene_name, group, log2FoldChange, padj, direction) |>
            arrange(group, desc(log2FoldChange)), file.path(OUT, "lipid_gene_groups.csv"))

has_se <- "lfcSE" %in% names(lip)
order_groups <- c(names(GROUPS), OTHER)
order_groups <- order_groups[order_groups %in% lip$group]

plot_group <- function(d, title, file_base) {
  d <- d |> arrange(log2FoldChange) |> mutate(
    label = make.unique(ifelse(is.na(gene_name), gene_id, gene_name)))
  d$label <- factor(d$label, levels = d$label)
  d$direction <- factor(d$direction, levels = names(COL))
  nsig <- sum(d$direction != "Not significant")
  p <- ggplot(d, aes(log2FoldChange, label, fill = direction)) +
    geom_vline(xintercept = 0, colour = "grey40", linewidth = 0.4) +
    geom_col(width = 0.7) +
    { if (has_se) geom_errorbarh(aes(xmin = log2FoldChange - lfcSE, xmax = log2FoldChange + lfcSE),
                                 height = 0.3, colour = INK, linewidth = 0.3) } +
    scale_fill_manual(values = COL, drop = FALSE, name = NULL) +
    scale_x_continuous(expand = expansion(mult = 0.08)) +
    labs(title = title,
         subtitle = sprintf("IFNB vs Mock, shrunk log2FC%s; %d genes, %d significant (padj<%g & |log2FC|>%g)",
                            if (has_se) " ± lfcSE" else "", nrow(d), nsig, PADJ_CUT, LFC_CUT),
         x = "log2 fold-change (IFNB vs Mock)", y = NULL) +
    theme_bw(base_size = 11) +
    theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(),
          panel.grid.major.x = element_line(colour = "grey92", linewidth = 0.3),
          panel.border = element_rect(colour = "grey70"),
          axis.text = element_text(colour = INK), axis.text.y = element_text(face = "italic"),
          axis.title = element_text(colour = INK),
          plot.title = element_text(colour = INK, face = "bold"),
          plot.subtitle = element_text(colour = "grey35", size = 9),
          legend.position = "bottom")
  h <- max(3.2, 0.24 * nrow(d) + 2)
  ggsave(paste0(file_base, ".png"), p, width = 7.5, height = h, dpi = 300)
}

slug <- function(x) gsub("^_|_$", "", gsub("[^a-z0-9]+", "_", tolower(x)))
msg <- list()
for (i in seq_along(order_groups)) {
  g <- order_groups[i]
  d <- lip |> filter(group == g)
  pages <- split(seq_len(nrow(d)), ceiling(rank(-d$log2FoldChange, ties.method = "first") / MAX_PER_PLOT))
  for (pg in seq_along(pages)) {
    dd <- d[pages[[pg]], ]
    suffix <- if (length(pages) > 1) sprintf("_part%d", pg) else ""
    ttl <- if (length(pages) > 1) sprintf("%s (part %d of %d)", g, pg, length(pages)) else g
    plot_group(dd, ttl, file.path(FIG, sprintf("%02d_%s%s", i, slug(g), suffix)))
  }
  msg[[g]] <- sprintf("%-55s %3d genes  (%d plot%s)", g, nrow(d), length(pages), if (length(pages) > 1) "s" else "")
}
message("\n", paste(unlist(msg), collapse = "\n"))
message("\nPNGs in: ", FIG, "\nGene-to-group table: ", file.path(OUT, "lipid_gene_groups.csv"))
