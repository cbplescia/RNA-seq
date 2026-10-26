#!/usr/bin/env Rscript
# For each target gene, find the genes that ACTIVATE or INHIBIT it according to STRING's curated "actions" data (pathway / literature-derived directed interactions), then draw vertical bar graphs of
# those regulators: three bars per regulator, left to right IFNB vs Mock, EBOV vs Mock, IFNB+EBOV vs Mock (shrunk log2FC +/- 1 lfcSE, stars = padj against Mock), grouped in panels by the gene they regulate.
# Targets (edit TARGETS): Suox, Sptlc1, Tst, Smpd2, Aox1, Aldh2, Cbs, Asah2, Aldh1b1, Cerk, Akr1b7, Sqor, Gm2a, Hexa, Pdxp, Glb1, Mocs3, Limk1.
# How regulators are found: STRING (mouse, taxon 10090) protein.actions file - rows with mode = activation / inhibition (add mode = expression with INCLUDE_EXPRESSION=1, which means one gene regulates the
# expression of the other without a stated sign) where the regulator acts on the target (directional, a_is_acting). Kept if the STRING score >= SCORE_MIN (default 400); at most N_REG regulators per target
# (highest score first) are drawn, but all are saved. Label "(+)" = activates the target, "(-)" = inhibits it. This is STRING's curated knowledge, which is sparse for some targets (several of yours are
# poorly annotated): a target with no directed regulator at this score is listed as having none - lower SCORE_MIN or tell me which regulators you know and I will add them by hand.
# It needs internet the first time (downloads two STRING files, ~100 MB, into <project>/stringdb_cache; the info file is already there if you ran script 88).
# Writes to requant_v2_outputs/figures_regulators/: regulators_of_targets.csv (every regulator found), regulator_log2fc.csv, bars_regulators_p<n>.png / .pdf
# Usage: Rscript 95_regulators_of_genes_bars.R   (PROJ from 'source setup_session.sh';  SCORE_MIN=700 N_REG=8 INCLUDE_EXPRESSION=1 STRING_VERSION=11.5 optional)

suppressPackageStartupMessages({ library(readr); library(dplyr); library(tidyr); library(ggplot2) })
select <- dplyr::select; filter <- dplyr::filter; count <- dplyr::count; mutate <- dplyr::mutate; summarise <- dplyr::summarise; arrange <- dplyr::arrange; desc <- dplyr::desc

TARGETS <- c("Suox", "Sptlc1", "Tst", "Smpd2", "Aox1", "Aldh2", "Cbs", "Asah2", "Aldh1b1", "Cerk", "Akr1b7", "Sqor", "Gm2a", "Hexa", "Pdxp", "Glb1", "Mocs3", "Limk1")
SCORE_MIN <- as.integer(Sys.getenv("SCORE_MIN", unset = "400")); N_REG <- as.integer(Sys.getenv("N_REG", unset = "12")); TARGETS_PER_PAGE <- 5; LOW_EXPR <- 50
INCL_EXPR <- identical(Sys.getenv("INCLUDE_EXPRESSION"), "1"); VER <- Sys.getenv("STRING_VERSION", unset = "12.0")
INK <- "#1F2933"
COMP <- c(IFNB_vs_Mock = "IFNB vs Mock", EBOV_vs_Mock = "EBOV vs Mock", IFNB_EBOV_vs_Mock = "IFNB+EBOV vs Mock"); CCOL <- c(`IFNB vs Mock` = "#2B6CB0", `EBOV vs Mock` = "#1B7837", `IFNB+EBOV vs Mock` = "#9BD08E")
PROJ <- Sys.getenv("PROJ", unset = file.path(
  Sys.getenv("HOME"),
  "Documents/aviti_rnaseq/2025/aviti",
  "20251009-100725_Maury_25161_WivfFacICOTHFAUJwSPeHBEFfCfNtCRpUSgVGUkz",
  "results/Maury_25161"))
RNA <- file.path(PROJ, "requant_v2_outputs")
OUT <- file.path(RNA, "figures_regulators"); dir.create(OUT, showWarnings = FALSE)
cache <- file.path(PROJ, "stringdb_cache"); dir.create(cache, showWarnings = FALSE)
options(timeout = 900)
stars <- function(p) ifelse(is.na(p), "", ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", "ns"))))
getf <- function(sub, name) { f <- file.path(cache, name); if (!file.exists(f)) { message("downloading ", name); download.file(sprintf("https://stringdb-downloads.org/download/%s.v%s/%s", sub, VER, name), f, mode = "wb", quiet = TRUE) }; f }

# ---- STRING files ----
inf <- read_tsv(getf("protein.info", sprintf("10090.protein.info.v%s.txt.gz", VER)), show_col_types = FALSE, progress = FALSE); names(inf)[1:2] <- c("id", "name")
act <- read_tsv(getf("protein.actions", sprintf("10090.protein.actions.v%s.txt.gz", VER)), show_col_types = FALSE, progress = FALSE)
names(act) <- sub("^#", "", names(act)); act <- act |> mutate(across(c(is_directional, a_is_acting), ~ as.character(.x) %in% c("t", "TRUE", "true", "1")))
message(sprintf("STRING actions: %d rows; modes: %s", nrow(act), paste(names(table(act$mode)), as.integer(table(act$mode)), sep = "=", collapse = ", ")))
modes <- c("activation", "inhibition", if (INCL_EXPR) "expression")

# ---- regulators of each target ----
tid <- inf |> mutate(key = tolower(name)); tt <- tibble(target = TARGETS, key = tolower(TARGETS)) |> left_join(tid |> select(key, target_id = id), by = "key")
if (any(is.na(tt$target_id))) message("Targets not found in STRING: ", paste(tt$target[is.na(tt$target_id)], collapse = ", "))
R <- bind_rows(lapply(seq_len(nrow(tt)), function(i) {
  if (is.na(tt$target_id[i])) return(NULL); t <- tt$target_id[i]
  a <- act |> filter(mode %in% modes, score >= SCORE_MIN, is_directional)
  r1 <- a |> filter(item_id_b == t, a_is_acting) |> transmute(regulator_id = item_id_a, mode, action, score)
  r2 <- a |> filter(item_id_a == t, !a_is_acting) |> transmute(regulator_id = item_id_b, mode, action, score)
  bind_rows(r1, r2) |> mutate(target = tt$target[i])
})) |> filter(regulator_id != "") |> left_join(inf |> select(regulator_id = id, regulator = name), by = "regulator_id")
R <- R |> mutate(role = case_when(mode == "activation" ~ "activates", mode == "inhibition" ~ "inhibits", TRUE ~ "regulates expression"), sign = case_when(role == "activates" ~ "(+)", role == "inhibits" ~ "(-)", TRUE ~ "(~)")) |>
  group_by(target, regulator, role, sign) |> summarise(score = max(score), .groups = "drop") |> arrange(target, desc(score))
per <- R |> count(target) ; message("Directed regulators per target (score >= ", SCORE_MIN, "):"); print(as.data.frame(tt |> select(target) |> left_join(per, by = "target") |> mutate(n = ifelse(is.na(n), 0L, n))), row.names = FALSE)

# ---- expression of the regulators ----
tabs <- lapply(names(COMP), function(nm) read_csv(file.path(RNA, paste0("DE_v2_", nm, "_annotated.csv")), show_col_types = FALSE) |> mutate(key = tolower(gene_name))); names(tabs) <- names(COMP)
D <- bind_rows(lapply(names(COMP), function(nm) tabs[[nm]] |> filter(key %in% tolower(R$regulator)) |> distinct(key, .keep_all = TRUE) |> transmute(key, gene_id, baseMean, comparison = COMP[[nm]], lfc = log2FoldChange, se = lfcSE, padj)))
R <- R |> mutate(key = tolower(regulator)); R$detected <- R$key %in% D$key
write_csv(R |> select(target, regulator, role, string_score = score, detected_in_data = detected), file.path(OUT, "regulators_of_targets.csv"))
message(sprintf("%d target-regulator pairs; %d with the regulator detected in your data", nrow(R), sum(R$detected)))
if (!any(R$detected)) stop("None of the regulators is in the DE tables - lower SCORE_MIN or add regulators by hand")
P <- R |> filter(detected) |> group_by(target) |> slice_head(n = N_REG) |> ungroup() |> inner_join(D, by = "key", relationship = "many-to-many") |>
  mutate(comparison = factor(comparison, levels = COMP), lab = stars(padj), vj = ifelse(lfc >= 0, 0, 1), ys = lfc + ifelse(lfc >= 0, 1, -1) * (se + 0.03),
         xl = paste0(regulator, ifelse(baseMean < LOW_EXPR, "†", ""), " ", sign), target = factor(target, levels = TARGETS))
write_csv(P |> arrange(target, desc(score), comparison) |> select(target, regulator, role, string_score = score, baseMean, comparison, log2FC = lfc, lfcSE = se, padj), file.path(OUT, "regulator_log2fc.csv"))

# ---- bar graphs: panels = target, pages of TARGETS_PER_PAGE targets ----
tg <- levels(P$target)[levels(P$target) %in% P$target]; pages <- split(tg, ceiling(seq_along(tg) / TARGETS_PER_PAGE)); dg <- position_dodge(width = 0.85)
for (pg in seq_along(pages)) {
  d <- P |> filter(target %in% pages[[pg]]) |> mutate(target = droplevels(target), panel = paste0("regulators of ", target))
  d$panel <- factor(d$panel, levels = paste0("regulators of ", levels(d$target))); d$key <- paste(d$target, d$xl, sep = "||")
  ord <- d |> distinct(target, key, score) |> arrange(target, desc(score)); d$key <- factor(d$key, levels = ord$key)
  p <- ggplot(d, aes(key, lfc, fill = comparison)) + geom_col(position = dg, width = 0.8) + geom_errorbar(aes(ymin = lfc - se, ymax = lfc + se), position = dg, width = 0.25, colour = INK, linewidth = 0.3) +
    geom_text(aes(y = ys, label = lab, vjust = vj), position = dg, colour = INK, size = 2.2) + geom_hline(yintercept = 0, colour = INK, linewidth = 0.4) +
    facet_grid(. ~ panel, scales = "free_x", space = "free_x", labeller = labeller(panel = label_wrap_gen(16))) + scale_fill_manual(values = CCOL, name = NULL) +
    scale_x_discrete(labels = function(x) sub("^.*\\|\\|", "", x)) + scale_y_continuous(expand = expansion(mult = 0.12)) +
    labs(title = "Regulators of the target genes: expression in IFNB, EBOV and IFNB+EBOV vs Mock", subtitle = sprintf("regulators from STRING actions (score >= %d, top %d per target); (+) activates, (-) inhibits%s the target; left to right: IFNB, EBOV, IFNB+EBOV; shrunk log2FC +/- 1 SE; stars = padj; dagger = baseMean < %d%s",
                                                                                                                      SCORE_MIN, N_REG, if (INCL_EXPR) ", (~) regulates expression of" else "", LOW_EXPR, if (length(pages) > 1) sprintf("; page %d of %d", pg, length(pages)) else ""),
         x = NULL, y = "log2 fold-change vs Mock") +
    theme_bw(base_size = 10) + theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(), panel.grid.major.y = element_line(colour = "grey92", linewidth = 0.3), panel.border = element_rect(colour = "grey70"),
          axis.text = element_text(colour = INK), axis.text.x = element_text(face = "italic", angle = 50, hjust = 1, size = 7.5), plot.title = element_text(colour = INK, face = "bold", size = 11),
          plot.subtitle = element_text(colour = "grey35", size = 7.5), strip.background = element_rect(fill = "grey95", colour = "grey70"), strip.text = element_text(colour = INK, face = "bold", size = 8), legend.position = "top", panel.spacing.x = unit(0.2, "lines"))
  w <- max(7, 0.55 * n_distinct(d$key) + 1.5 * n_distinct(d$target) + 2)
  ggsave(file.path(OUT, sprintf("bars_regulators_p%d.png", pg)), p, width = w, height = 5.6, dpi = 300, limitsize = FALSE); ggsave(file.path(OUT, sprintf("bars_regulators_p%d.pdf", pg)), p, width = w, height = 5.6, limitsize = FALSE)
}
none <- setdiff(TARGETS, unique(R$target[R$detected])); if (length(none)) message("Targets with no detected directed regulator at this score (nothing to plot): ", paste(none, collapse = ", "))
message("Saved to: ", OUT)
