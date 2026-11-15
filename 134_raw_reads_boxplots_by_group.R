#!/usr/bin/env Rscript
# Box plots of RAW READ COUNTS (Salmon NumReads summed over each gene's transcripts, from the UNFILTERED quantification, so genes removed by the pre-filter are included) in the untreated (Mock) samples,
# ALL groups in ONE figure (one panel per group; panel width follows the number of genes). Each gene's box = its replicates; the mean number of reads is printed above the box; the x-axis shows the mouse gene symbol and, underneath it, the name you gave when it differs.
# Groups and name mapping are defined in SETS below (mouse symbol = name as given):  edit freely.   Assumptions: B-tubulin = Tubb5; clec2 = Clec1b (CLEC-2); F4/80 = Adgre1; CD64 = Fcgr1; IBA1 = Aif1; Tim4 = Timd4; CRIg = Vsig4; CD11c = Itgax; Gpnmb1 = Gpnmb; CXCR1 = Cxcr1. CX3CR1 is in BOTH the capsule and the central vein group (a gene may appear in more than one group).
# Outputs (requant_v2_outputs/figures_raw_reads/): rawreads_all_groups.png / .pdf (everything in one figure; rawreads_<group>.png per group only with INDIVIDUAL=1), raw_reads_per_replicate.tsv (every gene x sample), raw_reads_summary.tsv (mean / min / max reads, TPM, kept by the pre-filter?).
# Env: ROWS (2 = two rows of panels, needs the R package patchwork: conda install -c conda-forge r-patchwork; 1 = one long row; default 2 if patchwork is installed), ROW1 (groups on the first row), INDIVIDUAL=1, GROUPS (default Mock_Mock = untreated; e.g. GROUPS=Mock_Mock,IFNB_Mock,IFNY_Mock to dodge several groups side by side), PROJ.   Usage: Rscript 134_raw_reads_boxplots_by_group.R   (PROJ from 'source setup_session.sh')

suppressPackageStartupMessages({ library(readr); library(dplyr); library(tidyr); library(ggplot2); library(DESeq2) })
select <- dplyr::select; filter <- dplyr::filter
SETS <- list(
  `Controls`     = c(Tubb5 = "B-tubulin", Alb = "Alb", Ins2 = "Ins2", Krt14 = "Krt14", Snap25 = "Snap25"),
  `General`      = c(Adgre1 = "F4/80", Fcgr1 = "CD64", Mertk = "MerTK", Aif1 = "IBA1"),
  `Kupffer`      = c(Timd4 = "Tim4", Clec4f = "Clec4f", Clec1b = "Clec2", Vsig4 = "CRIg/VSIG4", Folr2 = "FOLR2", Cd163 = "CD163", Cd5l = "CD5L", Marco = "MARCO"),
  `Capsule`      = c(Cx3cr1 = "CX3CR1", Itgax = "CD11c"),
  `Central vein` = c(Cx3cr1 = "CX3CR1", Cxcr1 = "CXCR1"),
  `Bile duct`    = c(Trem2 = "Trem2", Spp1 = "Spp1", Gpnmb = "Gpnmb1", Cd9 = "CD9"))
ORD <- c("Mock_Mock", "IFNB_Mock", "IFNY_Mock", "Mock_EBOV", "IFNB_EBOV", "IFNY_EBOV", "Mock_MARV", "IFNB_MARV", "IFNY_MARV", "Mock_EVLP", "Mock_MVLP")
LAB <- c(Mock_Mock = "Mock (untreated)", IFNB_Mock = "IFNB", IFNY_Mock = "IFNG", Mock_EBOV = "EBOV", IFNB_EBOV = "IFNB+EBOV", IFNY_EBOV = "IFNG+EBOV", Mock_MARV = "MARV", IFNB_MARV = "IFNB+MARV", IFNY_MARV = "IFNG+MARV", Mock_EVLP = "EVLP", Mock_MVLP = "MVLP")
COL <- c(`Mock (untreated)` = "#7F7F7F", IFNB = "#2B6CB0", IFNG = "#D69E2E", EBOV = "#1B7837", `IFNB+EBOV` = "#9BD08E", `IFNG+EBOV` = "#E8C77A", MARV = "#8C510A", `IFNB+MARV` = "#C7A27C", `IFNG+MARV` = "#DCC9A8", EVLP = "#6A51A3", MVLP = "#CBA4D8")
GRP <- trimws(strsplit(Sys.getenv("GROUPS", "Mock_Mock"), ",")[[1]]); INK <- "#1F2933"
PROJ <- Sys.getenv("PROJ", unset = file.path(
  Sys.getenv("HOME"),
  "Documents/aviti_rnaseq/2025/aviti",
  "20251009-100725_Maury_25161_WivfFacICOTHFAUJwSPeHBEFfCfNtCRpUSgVGUkz",
  "results/Maury_25161"))
RNA <- file.path(PROJ, "requant_v2_outputs"); QUANT <- file.path(PROJ, "salmon_quant_v2"); GTF <- file.path(PROJ, "ref", "Mus_musculus.GRCm39.109.gtf"); OUT <- file.path(RNA, "figures_raw_reads"); dir.create(OUT, showWarnings = FALSE)
nov <- function(x) sub("\\..*$", "", x)
dds <- readRDS(file.path(RNA, "dds_v2_decoyAware_hostNormalised.rds")); cd <- as.data.frame(colData(dds)); cd$sample <- rownames(cd); cd$group <- as.character(cd$group); kept <- nov(rownames(dds))
if (!all(GRP %in% cd$group)) stop("unknown group(s): ", paste(setdiff(GRP, cd$group), collapse = ", ")); cd <- cd |> filter(group %in% GRP)
allg <- unlist(lapply(SETS, names), use.names = FALSE)
pat <- tempfile(); writeLines(sprintf('gene_name "%s"', allg), pat); message("looking the genes up in the GTF (one pass) ...")
ln <- tryCatch(readLines(pipe(sprintf("grep -F -f '%s' '%s'", pat, GTF))), error = function(e) character(0), warning = function(w) character(0)); f <- strsplit(ln, "\t"); f <- f[vapply(f, length, 1L) >= 9]; f <- f[vapply(f, function(x) x[3] == "transcript", TRUE)]
gv <- function(x, k) sub(paste0('.*', k, ' "([^"]+)".*'), "\\1", x)
tab <- tibble(gene = vapply(f, function(x) gv(x[9], "gene_name"), ""), gene_id = nov(vapply(f, function(x) gv(x[9], "gene_id"), "")), tx = nov(vapply(f, function(x) gv(x[9], "transcript_id"), ""))) |> filter(tolower(gene) %in% tolower(allg))
miss <- allg[!tolower(allg) %in% tolower(tab$gene)]; if (length(miss)) message("NOT FOUND in the GTF under that symbol: ", paste(miss, collapse = ", "))
message("reading ", nrow(cd), " quant.sf files ..."); rows <- list()
for (s in cd$sample) { q <- read_tsv(file.path(QUANT, s, "quant.sf"), col_types = "cdddd", progress = FALSE); q$tx <- nov(q$Name)
  m <- q |> inner_join(tab, by = "tx") |> group_by(gene, gene_id) |> summarise(reads = sum(NumReads), TPM = sum(TPM), .groups = "drop"); m$sample <- s; rows[[s]] <- m }
genes_found <- tab |> distinct(gene, gene_id)
L <- bind_rows(rows) |> right_join(expand_grid(genes_found, sample = cd$sample), by = c("gene", "gene_id", "sample")) |> mutate(reads = coalesce(reads, 0), TPM = coalesce(TPM, 0)) |> left_join(cd |> select(sample, group), by = "sample") |>
  mutate(kept = gene_id %in% kept, group = factor(ifelse(group %in% names(LAB), LAB[group], group), levels = unname(LAB[GRP])), key = tolower(gene))
write_tsv(L |> select(gene, gene_id, sample, group, reads, TPM, passes_prefilter = kept) |> arrange(gene, group, sample), file.path(OUT, "raw_reads_per_replicate.tsv"))
S <- L |> group_by(gene, gene_id, kept, group) |> summarise(n = n(), mean_reads = mean(reads), min_reads = min(reads), max_reads = max(reads), mean_TPM = mean(TPM), .groups = "drop"); write_tsv(S, file.path(OUT, "raw_reads_summary.tsv"))
brk <- c(0, 1, 3, 10, 30, 100, 300, 1000, 3000, 10000, 30000, 100000, 300000)
# one data frame for all sets (a gene in two sets appears in both panels)
mk <- function(nm) {
  sy <- names(SETS[[nm]]); sy <- sy[tolower(sy) %in% L$key]; if (!length(sy)) { message("group '", nm, "': none of its genes were found"); return(NULL) }
  lab <- vapply(sy, function(g) { gg <- L$gene[match(tolower(g), L$key)]; gv_ <- SETS[[nm]][[g]]; if (tolower(gv_) == tolower(gg)) gg else paste0(gg, "\n(", gv_, ")") }, "")
  L |> filter(key %in% tolower(sy)) |> mutate(set = nm, xl = lab[match(key, tolower(sy))], xo = match(key, tolower(sy))) }
D <- bind_rows(lapply(names(SETS), mk)); if (!nrow(D)) stop("no gene was found")
D$set <- factor(D$set, levels = names(SETS)); nset <- D |> distinct(set, key) |> count(set); D$setlab <- factor(paste0(D$set, "  (", nset$n[match(D$set, nset$set)], " gene", ifelse(nset$n[match(D$set, nset$set)] > 1, "s", ""), ")"), levels = paste0(nset$set, "  (", nset$n, " gene", ifelse(nset$n > 1, "s", ""), ")"))
D$xl <- factor(D$xl, levels = unique(D$xl[order(D$set, D$xo)]))
sm <- D |> group_by(setlab, xl, group) |> summarise(mean_reads = mean(reads), top = max(reads), .groups = "drop") |> mutate(txt = format(round(mean_reads), big.mark = ","))
panel <- function(dd, ss, show_legend) ggplot(dd, aes(xl, reads, fill = group)) + geom_boxplot(position = position_dodge(width = 0.8), width = 0.6, outlier.shape = NA, alpha = 0.8, colour = INK, linewidth = 0.4) +
  geom_point(aes(group = group), position = position_jitterdodge(jitter.width = 0.08, dodge.width = 0.8, seed = 1), size = 1.9, shape = 21, colour = INK, fill = "white", stroke = 0.5) +
  geom_text(data = ss, mapping = aes(x = xl, y = top, label = txt, group = group), position = position_dodge(width = 0.8), vjust = -0.9, size = 2.7, colour = INK, inherit.aes = FALSE) +
  facet_grid(. ~ setlab, scales = "free_x", space = "free_x") + scale_fill_manual(values = COL, name = NULL) + scale_y_continuous(trans = scales::pseudo_log_trans(base = 10), breaks = brk, labels = scales::label_comma(), expand = expansion(mult = c(0.05, 0.18))) +
  labs(x = NULL, y = "raw reads") + theme_bw(base_size = 10.5) + theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(), panel.border = element_rect(colour = "grey70"), axis.text = element_text(colour = INK), axis.text.x = element_text(size = 8),
        strip.background = element_rect(fill = "grey93", colour = "grey70"), strip.text = element_text(face = "bold", colour = INK), legend.position = if (show_legend && length(GRP) > 1) "top" else "none")
ttl <- paste0("Raw read counts (Salmon, unfiltered), ", paste(unname(LAB[GRP]), collapse = " / "), ": all gene groups"); subt <- "every dot = one sample; number above each box = mean reads; log scale; mouse symbol with the name you gave in brackets where it differs"
ROWS <- Sys.getenv("ROWS", if (requireNamespace("patchwork", quietly = TRUE)) "2" else "1"); r1 <- trimws(strsplit(Sys.getenv("ROW1", "Controls,General,Capsule,Central vein"), ",")[[1]])
if (ROWS == "2" && requireNamespace("patchwork", quietly = TRUE)) {
  s1 <- levels(D$setlab)[sub("  \\(.*$", "", levels(D$setlab)) %in% r1]; s2 <- setdiff(levels(D$setlab), s1)
  if (length(s1) && length(s2)) { d1 <- D |> filter(setlab %in% s1); d2 <- D |> filter(setlab %in% s2)
    p <- patchwork::wrap_plots(panel(d1, sm |> filter(setlab %in% s1), TRUE), panel(d2, sm |> filter(setlab %in% s2), FALSE), ncol = 1) + patchwork::plot_annotation(title = ttl, subtitle = subt, theme = theme(plot.title = element_text(colour = INK, face = "bold", size = 13), plot.subtitle = element_text(colour = "grey35", size = 9)))
    wmax <- max(nrow(distinct(d1, xl, set)), nrow(distinct(d2, xl, set))); w <- max(8, 2.2 + 0.95 * wmax * max(1, length(GRP) * 0.6)); h <- 9.2 }
  else ROWS <- "1" }
if (ROWS != "2" || !exists("p")) { if (ROWS == "2") message("patchwork is not installed - drawing everything in one row (install it with: conda install -c conda-forge r-patchwork)"); nx <- nrow(distinct(D, xl, set)) + 0; p <- panel(D, sm, TRUE) + labs(title = ttl, subtitle = subt) + theme(plot.title = element_text(colour = INK, face = "bold", size = 13), plot.subtitle = element_text(colour = "grey35", size = 9)); w <- max(8, 2 + 0.95 * nx * max(1, length(GRP) * 0.6)); h <- 5.4 }
ggsave(file.path(OUT, "rawreads_all_groups.png"), p, width = w, height = h, dpi = 300, limitsize = FALSE); ggsave(file.path(OUT, "rawreads_all_groups.pdf"), p, width = w, height = h, limitsize = FALSE)
if (Sys.getenv("INDIVIDUAL") == "1") for (nm in names(SETS)) { dd <- D |> filter(set == nm); if (!nrow(dd)) next; ss <- sm |> filter(setlab %in% unique(dd$setlab)); pp <- panel(dd, ss, TRUE) + labs(title = nm); fn <- paste0("rawreads_", gsub("[^A-Za-z0-9]+", "_", tolower(nm))); ggsave(file.path(OUT, paste0(fn, ".png")), pp, width = max(4.5, 1.4 + 0.95 * n_distinct(dd$xl) * max(1, length(GRP) * 0.6)), height = 4.8, dpi = 300) }
for (nm in names(SETS)) { dd <- sm |> filter(grepl(paste0("^", nm, "  \\("), setlab)); if (nrow(dd)) message(sprintf("%-12s %s", nm, paste(sprintf("%s %s", gsub("\n.*", "", dd$xl), dd$txt), collapse = " | "))) }
message("\nSaved to: ", OUT)
