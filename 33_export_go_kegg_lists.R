#!/usr/bin/env Rscript
# Export full reference lists to CSV (written to requant_v2_outputs/reference_lists/):
#   go_terms_all.csv          every GO term: id, ontology (BP / MF / CC), name, definition   (offline, GO.db)
#   go_term_parents.csv       parent -> child links, so you can see the hierarchy               (offline, GO.db)
#   kegg_mmu_pathways.csv     every mouse KEGG pathway with its KEGG BRITE classes (A / B)     (needs internet)
# KEGG's REST API is free for academic use; see kegg.jp for licensing of bulk data.

suppressPackageStartupMessages({ library(readr); library(dplyr) })
select <- dplyr::select; filter <- dplyr::filter

PROJ <- Sys.getenv("PROJ", unset = file.path(
  Sys.getenv("HOME"),
  "Documents/aviti_rnaseq/2025/aviti",
  "20251009-100725_Maury_25161_WivfFacICOTHFAUJwSPeHBEFfCfNtCRpUSgVGUkz",
  "results/Maury_25161"))
OUT <- file.path(PROJ, "requant_v2_outputs", "reference_lists")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

# ---- GO (offline) ------------------------------------------------------------------------------
if (requireNamespace("GO.db", quietly = TRUE)) {
  ids <- AnnotationDbi::keys(GO.db::GO.db, keytype = "GOID")
  go <- suppressMessages(AnnotationDbi::select(GO.db::GO.db, keys = ids, keytype = "GOID",
                                               columns = c("ONTOLOGY", "TERM", "DEFINITION")))
  write_csv(go, file.path(OUT, "go_terms_all.csv"))
  message("GO terms: ", nrow(go), "  (", paste(names(table(go$ONTOLOGY)), table(go$ONTOLOGY), collapse = ", "), ")")
  par <- do.call(rbind, lapply(list(BP = GO.db::GOBPPARENTS, MF = GO.db::GOMFPARENTS, CC = GO.db::GOCCPARENTS), function(m) {
    l <- AnnotationDbi::as.list(m)
    data.frame(child = rep(names(l), lengths(l)), parent = unlist(l, use.names = FALSE),
               relation = unlist(lapply(l, names), use.names = FALSE), stringsAsFactors = FALSE)
  }))
  par <- par[par$parent != "all", ]
  write_csv(par, file.path(OUT, "go_term_parents.csv"))
  message("GO parent-child links: ", nrow(par))
} else message("GO.db not installed - skipping GO (BiocManager::install('GO.db'))")

# ---- KEGG (needs internet) -------------------------------------------------------------------------
tryCatch({
  pw <- readLines("https://rest.kegg.jp/list/pathway/mmu")
  pw <- do.call(rbind, strsplit(pw, "\t", fixed = TRUE))
  pw <- tibble(pathway_id = sub("^path:", "", pw[, 1]),
               pathway_name = sub(" - Mus musculus.*$", "", pw[, 2]))

  cls <- tibble(pathway_id = character(), class_A = character(), class_B = character())
  if (requireNamespace("jsonlite", quietly = TRUE)) {
    js <- jsonlite::fromJSON(paste(readLines("https://rest.kegg.jp/get/br:br08901/json", warn = FALSE), collapse = ""),
                             simplifyVector = FALSE)
    rows <- list()
    walk <- function(node, a = NA, b = NA, depth = 0) {
      kids <- node$children
      if (is.null(kids)) {   # leaf: "00010 Glycolysis / Gluconeogenesis"
        num <- sub("^([0-9]{5}).*$", "\\1", node$name)
        rows[[length(rows) + 1]] <<- data.frame(pathway_id = paste0("mmu", num), class_A = a, class_B = b)
      } else for (k in kids) walk(k, if (depth == 0) sub("^[0-9]+ ", "", k$name) else a,
                                  if (depth == 1) sub("^[0-9]+ ", "", k$name) else b, depth + 1)
    }
    walk(js)
    cls <- unique(do.call(rbind, rows))
  } else message("jsonlite not installed - KEGG classes (BRITE) skipped; install.packages('jsonlite')")

  out <- left_join(pw, cls, by = "pathway_id")
  write_csv(out, file.path(OUT, "kegg_mmu_pathways.csv"))
  message("Mouse KEGG pathways: ", nrow(out), "  (with a class: ", sum(!is.na(out$class_A)), ")")
}, error = function(e) message("KEGG not reachable (", conditionMessage(e), ") - use the links in the chat instead"))
message("\nWritten to: ", OUT)
