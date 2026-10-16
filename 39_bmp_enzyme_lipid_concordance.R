#!/usr/bin/env Rscript
# BMP/LBPA-related enzyme transcripts (RNA-seq, IFNB vs Mock) set against the matching BMP lipid readouts
# (pilot lipidomics, IFNB vs untreated), module by module, on one log2FC axis.
#
# Each module = enzymes + a testable prediction about lipid chain composition:
#   A  18:1 chains       Scd1, Scd2 (make 18:1)          enzyme DOWN  -> species with 18:1 should go DOWN
#   B  PUFA / DHA        Fads1, Fads2, Elovl5, Acsl4,    enzymes DOWN -> PUFA / DHA species should go DOWN
#                        Lpcat3, Agpat3
#   C  BMP synthesis &   Cln5, Pld3, Pld4, Pla2g15,      mixed roles  -> not scored; shown against the total
#      turnover          Pgs1, Cds1, Cds2, Ptpmt1        BMP signal
# Lipid readouts compare the mean log2FC of species that HAVE a chain feature with those that do not (Wilcoxon p),
# for BMP and, as a check on specificity, for the other diacyl phospholipids (PC/PE/PS/PG). Only species whose name
# lists every chain are used; standards, oxidised species and rows with missing values are excluded.
# "Concordant" = the sign of the mean enzyme log2FC matches the sign of the lipid difference (with minus without).
# This is a consistency check on direction, not proof of mechanism; mRNA is not enzyme activity.
#
# Writes to <project>/lipidomics_chain_length/:
#   P_enzyme_transcripts_vs_BMP_lipids.png/.pdf   the combined bar graph
#   bmp_enzyme_lipid_concordance.csv              one row per module x lipid readout
# Usage: Rscript 39_bmp_enzyme_lipid_concordance.R [path/to/lipidomics.csv]

suppressPackageStartupMessages({ library(readr); library(dplyr); library(tidyr); library(ggplot2) })
select <- dplyr::select; first <- dplyr::first; rename <- dplyr::rename
slice <- dplyr::slice; filter <- dplyr::filter; desc <- dplyr::desc; count <- dplyr::count

PROJ <- Sys.getenv("PROJ", unset = file.path(
  Sys.getenv("HOME"),
  "Documents/aviti_rnaseq/2025/aviti",
  "20251009-100725_Maury_25161_WivfFacICOTHFAUJwSPeHBEFfCfNtCRpUSgVGUkz",
  "results/Maury_25161"))
RNA_OUT <- file.path(PROJ, "requant_v2_outputs")
OUT <- file.path(PROJ, "lipidomics_chain_length"); dir.create(OUT, showWarnings = FALSE)
INK <- "#1F2933"

MODULES <- list(
  A = list(name = "A. 18:1 (monounsaturated) chains", genes = c("Scd1", "Scd2"), scored = TRUE),
  B = list(name = "B. PUFA / DHA supply and incorporation", genes = c("Fads1", "Fads2", "Elovl5", "Acsl4", "Lpcat3", "Agpat3"), scored = TRUE),
  C = list(name = "C. BMP synthesis and turnover", genes = c("Cln5", "Pld3", "Pld4", "Pla2g15", "Pgs1", "Cds1", "Cds2", "Ptpmt1"), scored = FALSE))

# ---- RNA ---------------------------------------------------------------------------------------------------
de <- read_csv(file.path(RNA_OUT, "DE_v2_IFNB_vs_Mock_protein_coding.csv"), show_col_types = FALSE)
enz <- bind_rows(lapply(names(MODULES), function(m) tibble(module = m, gene = MODULES[[m]]$genes))) |>
  left_join(de |> select(gene = gene_name, log2FoldChange, lfcSE, padj), by = "gene")
if (any(is.na(enz$log2FoldChange))) message("Not in the RNA-seq data (left out): ", paste(enz$gene[is.na(enz$log2FoldChange)], collapse = ", "))
enz <- enz |> filter(!is.na(log2FoldChange))

# ---- lipids ------------------------------------------------------------------------------------------------
a <- commandArgs(trailingOnly = TRUE)
fname <- if (length(a)) a[1] else "01_lipidomics_ranked_IFNB_vs_Untreated.csv"
cand <- c(fname, file.path(PROJ, fname), file.path(getwd(), fname)); rf <- cand[file.exists(cand)][1]
if (is.na(rf)) stop("Cannot find ", fname, " - give its full path as the first argument")
raw <- read_csv(rf, show_col_types = FALSE)
miss <- (if ("any_missing_Untreated" %in% names(raw)) raw$any_missing_Untreated else FALSE) |
        (if ("any_missing_IFNB" %in% names(raw)) raw$any_missing_IFNB else FALSE)
L <- raw[!grepl("\\+D[0-9]|-D[0-9]|\\(s\\)|COOH|OOH", raw$Lipid_ID) & !miss, ] |>
  transmute(lipid = Lipid_ID, class = Class, lfc = log2FC_IFNB_vs_Untreated, untreated = Untreated_Norm, ifnb = B_Norm)
ch <- regmatches(L$lipid, gregexpr("(?<![0-9])[0-9]{1,2}:[0-9]{1,2}", L$lipid, perl = TRUE))
L$chains <- ch
Lx <- L[lengths(ch) >= 2, ]                       # species whose name lists every chain
bmp <- Lx[Lx$class %in% c("LBPA", "BMP"), ]
dia <- Lx[Lx$class %in% c("PC", "PE", "PS", "PG"), ]
message("BMP species with explicit chains: ", nrow(bmp), "; PC/PE/PS/PG species: ", nrow(dia))

has_181 <- function(z) "18:1" %in% z
has_pufa <- function(z) any(as.numeric(sub("^.*:", "", z)) >= 4)
has_dha <- function(z) "22:6" %in% z
readout <- function(d, fn, module, label) {
  m <- vapply(d$chains, fn, logical(1))
  w <- d$lfc[m]; wo <- d$lfc[!m]
  if (length(w) < 3 || length(wo) < 3) return(NULL)
  tibble(module = module, readout = label, n_with = length(w), mean_with = mean(w), se_with = sd(w) / sqrt(length(w)),
         n_without = length(wo), mean_without = mean(wo), se_without = sd(wo) / sqrt(length(wo)),
         difference = mean(w) - mean(wo), p_wilcoxon = suppressWarnings(wilcox.test(w, wo)$p.value))
}
ro <- bind_rows(
  readout(bmp, has_181, "A", "BMP: species with 18:1"),
  readout(dia, has_181, "A", "PC/PE/PS/PG: species with 18:1"),
  readout(bmp, has_pufa, "B", "BMP: species with a PUFA chain (>=4 DB)"),
  readout(bmp, has_dha, "B", "BMP: species with 22:6 (DHA)"),
  readout(dia, has_pufa, "B", "PC/PE/PS/PG: species with a PUFA chain"))
bmp_all <- L[L$class %in% c("LBPA", "BMP"), ]
bmp_total <- tibble(module = "C", readout = "BMP: total signal (all species)",
                    n_with = nrow(bmp_all), mean_with = log2(sum(bmp_all$ifnb) / sum(bmp_all$untreated)), se_with = NA_real_,
                    n_without = NA_integer_, mean_without = NA_real_, se_without = NA_real_, difference = NA_real_, p_wilcoxon = NA_real_)
ro <- bind_rows(ro, bmp_total)

# ---- concordance table ----------------------------------------------------------------------------------------
em <- enz |> group_by(module) |> summarise(enzyme_mean_log2FC = mean(log2FoldChange), n_enzymes = n(),
                                          n_sig_down = sum(padj < 0.05 & log2FoldChange < 0), n_sig_up = sum(padj < 0.05 & log2FoldChange > 0), .groups = "drop")
tab <- ro |> left_join(em, by = "module") |>
  mutate(scored = MODULES[module] |> vapply(function(x) x$scored, logical(1)) |> unname(),
         prediction = ifelse(scored, ifelse(enzyme_mean_log2FC < 0, "lipids with this feature DOWN", "lipids with this feature UP"), "not scored"),
         observed = ifelse(is.na(difference), NA_character_, ifelse(difference < 0, "DOWN", "UP")),
         verdict = case_when(!scored ~ "not scored", sign(enzyme_mean_log2FC) == sign(difference) ~ "concordant", TRUE ~ "DISCORDANT"))
write_csv(tab, file.path(OUT, "bmp_enzyme_lipid_concordance.csv"))
message("\nEnzyme transcripts vs lipid readouts:")
print(as.data.frame(tab |> transmute(module, readout, enzyme_mean = round(enzyme_mean_log2FC, 2) , with = round(mean_with, 2),
                                     without = round(mean_without, 2), diff = round(difference, 2), p = signif(p_wilcoxon, 2), verdict)), row.names = FALSE)

# ---- combined bar graph -----------------------------------------------------------------------------------------
KIND <- c(`Enzyme transcript: up` = "#C53030", `Enzyme transcript: down` = "#2B6CB0", `Enzyme transcript: n.s.` = "#B8B8B8",
          `Lipid species WITH feature` = "#6A3D9A", `Lipid species WITHOUT feature` = "#CAB2D6", `Lipid total` = "#6A3D9A")
e <- enz |> transmute(module, label = gene, value = log2FoldChange, err = lfcSE,
        kind = case_when(padj < 0.05 & log2FoldChange > 0 ~ "Enzyme transcript: up", padj < 0.05 & log2FoldChange < 0 ~ "Enzyme transcript: down",
                         TRUE ~ "Enzyme transcript: n.s."),
        note = case_when(padj < 0.001 ~ "***", padj < 0.01 ~ "**", padj < 0.05 ~ "*", TRUE ~ "ns"), order = 1)
l <- bind_rows(
  tab |> filter(!is.na(mean_without)) |> transmute(module, label = sprintf("[lipid] %s (n=%d)", readout, n_with), value = mean_with, err = se_with,
        kind = "Lipid species WITH feature", note = sprintf("diff %+.2f, p=%.2g", difference, p_wilcoxon), order = 2),
  tab |> filter(!is.na(mean_without)) |> transmute(module, label = sprintf("[lipid] %s - without (n=%d)", sub(": species with", ": species lacking", readout), n_without), value = mean_without, err = se_without,
        kind = "Lipid species WITHOUT feature", note = "", order = 2.5),
  tab |> filter(is.na(mean_without)) |> transmute(module, label = sprintf("[lipid] %s", readout), value = mean_with, err = NA_real_,
        kind = "Lipid total", note = "", order = 2))
pd <- bind_rows(e, l) |> mutate(mod_name = vapply(module, function(m) MODULES[[m]]$name, character(1)))
pd <- pd |> arrange(module, order) |> mutate(row = row_number())
pd$label_f <- factor(pd$label, levels = rev(unique(pd$label)))
pd$kind <- factor(pd$kind, levels = names(KIND))
pd$mod_name <- factor(pd$mod_name, levels = vapply(MODULES, function(x) x$name, character(1)))
p <- ggplot(pd, aes(value, label_f, fill = kind)) +
  geom_vline(xintercept = 0, colour = "grey40", linewidth = 0.3) +
  geom_col(width = 0.7) +
  geom_errorbar(aes(xmin = value - err, xmax = value + err), width = 0.25, colour = INK, linewidth = 0.3, orientation = "y", na.rm = TRUE) +
  geom_text(aes(x = ifelse(value >= 0, value + ifelse(is.na(err), 0, err), value - ifelse(is.na(err), 0, err)), label = note),
            hjust = ifelse(pd$value >= 0, -0.1, 1.1), size = 2.7, colour = INK) +
  facet_grid(mod_name ~ ., scales = "free_y", space = "free_y") +
  scale_fill_manual(values = KIND, name = NULL, breaks = names(KIND)[1:5]) +
  scale_x_continuous(expand = expansion(mult = 0.25)) +
  labs(title = "BMP-related enzyme transcripts and the matching BMP lipid changes (IFNB)",
       subtitle = "blue/red = enzyme transcript log2FC vs Mock (stars = padj); purple = mean lipid log2FC vs untreated (+/- SE) of species with / without the chain feature",
       x = "log2 fold-change", y = NULL) +
  theme_bw(base_size = 10) + theme(panel.grid.minor = element_blank(), panel.border = element_rect(colour = "grey70"),
        axis.text = element_text(colour = INK), strip.background = element_rect(fill = "grey95", colour = "grey70"),
        strip.text.y = element_text(angle = 0, hjust = 0, colour = INK, face = "bold"),
        plot.title = element_text(colour = INK, face = "bold"), plot.subtitle = element_text(colour = "grey35", size = 8),
        legend.position = "bottom") + guides(fill = guide_legend(nrow = 2))
ggsave(file.path(OUT, "P_enzyme_transcripts_vs_BMP_lipids.png"), p, width = 11, height = 0.27 * nrow(pd) + 3.5, dpi = 300)
ggsave(file.path(OUT, "P_enzyme_transcripts_vs_BMP_lipids.pdf"), p, width = 11, height = 0.27 * nrow(pd) + 3.5)
message("\nFigure and table in: ", OUT)
