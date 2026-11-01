#!/usr/bin/env python3
"""Genes whose promoters share transcription-factor motifs with Top1 (not just genes that sit next to it), and their EBOV vs Mock differential expression.
How it works (needs 108_promoter_tf_groups.py in the same folder - its promoter / motif code is reused):
  1. the Top1 promoter (canonical transcript, UP bp upstream .. DOWN bp downstream of the TSS; default 1500 / 500) is scanned with the JASPAR vertebrate matrices (relative score >= REL, default 0.88); a motif is kept if the TF is
     expressed in your cells (baseMean >= MIN_TF_EXPR, default 10, matched on the TF name) - capped at M_MAX (default 120) motifs, best Top1 score first (set TFS=Sp1,Nfya,Usf1,... to use only those TFs)
  2. the promoters of every expressed protein-coding gene (baseMean >= MIN_BASEMEAN, default 10) are scanned with those motifs
  3. each gene gets a SHARED-MOTIF SCORE = sum over the Top1 motifs present in its promoter of weight = -log2(fraction of all promoters carrying that motif), so a motif found in half of all promoters (GC-rich / SP-like) counts little
     and a rare motif counts a lot; n_shared motifs and the TF names are reported; the top TOP_N (default 30) genes are plotted
Outputs (requant_v2_outputs/top1_shared_tf/): top1_shared_tf_genes.tsv (every gene, score, n shared, percentile, shared TFs), top1_motifs_used.tsv, bars_top1_shared_tf_<COMP>.png/.pdf (score panel on top, log2FC bars with padj stars below, genes in score order),
  top1_shared_motif_map.png/.pdf (genes x shared motifs, tile = best relative score).   Bars: shrunk log2FC +/- 1 SE; stars = padj (ns = not significant); dagger = baseMean < 50.
CAUTION: sharing motifs means a gene COULD be regulated by the same TFs, not that it is; promoters rich in GC share many motifs by chance (the weighting reduces, not removes, this); enhancers are not scanned; mouse macrophage ChIP-seq / ATAC would be the real test.
Env: GENE=Top1  COMP=EBOV_vs_Mock  TOP_N=30  REL=0.88  UP=1500  DOWN=500  M_MAX=120  MIN_BASEMEAN=10  MIN_TF_EXPR=10  TFS=comma list  JASPAR_FILE=path  PROJ=path       Usage: python 125_top1_shared_tf_genes.py       (python 3, numpy, scipy, matplotlib; ~10-20 min)
"""
import os, sys, re, csv, importlib.util
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

here = os.path.dirname(os.path.abspath(__file__)); p108 = os.path.join(here, "108_promoter_tf_groups.py")
if not os.path.exists(p108): sys.exit("Put 108_promoter_tf_groups.py in the same folder as this script (its promoter / motif code is reused)")
spec = importlib.util.spec_from_file_location("p108", p108); M = importlib.util.module_from_spec(spec); spec.loader.exec_module(M)
RNA = M.RNA; GENE = os.environ.get("GENE", "Top1"); COMP = os.environ.get("COMP", "EBOV_vs_Mock"); TOP_N = int(os.environ.get("TOP_N", "30")); REL = M.REL; M_MAX = int(os.environ.get("M_MAX", "120"))
MIN_BASEMEAN = float(os.environ.get("MIN_BASEMEAN", "10")); MIN_TF_EXPR = float(os.environ.get("MIN_TF_EXPR", "10")); TFS = [t.strip().lower() for t in os.environ.get("TFS", "").split(",") if t.strip()]
OUT = os.path.join(RNA, "top1_shared_tf" if GENE == "Top1" else "%s_shared_tf" % GENE.lower()); os.makedirs(OUT, exist_ok=True); INK = "#1F2933"; say = M.say
FILL = {"EBOV_vs_Mock": "#1B7837", "IFNB_vs_Mock": "#2B6CB0", "IFNB_EBOV_vs_Mock": "#9BD08E"}.get(COMP, "#1B7837")

de = list(csv.DictReader(open(os.path.join(RNA, "DE_v2_%s_annotated.csv" % COMP))))
def f(x):
    try: return float(x)
    except Exception: return np.nan
D = {}
for r in de:
    nm = (r.get("gene_name") or "").strip()
    if nm and nm.lower() not in D: D[nm.lower()] = dict(name=nm, gene_id=r["gene_id"], bm=f(r["baseMean"]), lfc=f(r["log2FoldChange"]), se=f(r["lfcSE"]), padj=f(r["padj"]), bt=r.get("gene_biotype", ""))
if GENE.lower() not in D: sys.exit("%s is not in the DE table (removed by the pre-filter?)" % GENE)
universe = [d["name"] for d in D.values() if d["bt"] == "protein_coding" and d["bm"] >= MIN_BASEMEAN]; say("%d expressed protein-coding genes (baseMean >= %g) will be scanned" % (len(universe), MIN_BASEMEAN))
if GENE not in universe: universe.append(GENE)
tx = M.read_gtf(universe); fa = M.Fasta(M.GENOME); motifs, src = M.load_motifs(); pm = M.prep(motifs); say("%d motifs from %s" % (len(pm), src))

def promoter_seq(g):
    tl = tx.get(g.lower())
    if not tl: return None
    p = M.promoter_of(tl); gs, ge = (p["tss"] - M.UP, p["tss"] + M.DOWN - 1) if p["strand"] == "+" else (p["tss"] - M.DOWN + 1, p["tss"] + M.UP)
    if p["chrom"] not in fa.idx: return None
    s = fa.get(p["chrom"], gs, ge); s = s if p["strand"] == "+" else M.revcomp(s)
    return s if len(s) >= 200 else None
def tf_expr(name):                                                   # best baseMean among the TF names in a JASPAR name (Fos::Jun -> Fos, Jun)
    toks = [t.strip().lower() for t in re.split(r"::", re.sub(r"\(var\.\d+\)", "", name))]
    v = [D[t]["bm"] for t in toks if t in D]; return (max(v) if v else np.nan), toks

q = promoter_seq(GENE)
if q is None: sys.exit("could not cut the %s promoter (gene not in the GTF / genome?)" % GENE)
s_top = M.best_scores(q, pm); order = np.argsort(-s_top); keep = []
for k in order:
    if s_top[k] < REL: break
    bm, toks = tf_expr(pm[k][1])
    if TFS and not any(t in TFS for t in toks): continue
    if not TFS and not (bm >= MIN_TF_EXPR): continue
    keep.append(k)
    if len(keep) >= M_MAX: break
if not keep: sys.exit("no motif passes in the %s promoter with an expressed TF (try REL=0.85 or TFS=...)" % GENE)
pm2 = [pm[k] for k in keep]; say("%s promoter: %d motifs kept (score >= %g, TF expressed): %s ..." % (GENE, len(pm2), REL, ", ".join(m[1] for m in pm2[:12])))
names, S = [], []
for i, g in enumerate(universe):
    s = promoter_seq(g)
    if s is None: continue
    names.append(g); S.append(M.best_scores(s, pm2))
    if len(names) % 200 == 0: say("  scanned %d promoters" % len(names))
S = np.array(S); P = S >= REL; frac = P.mean(0); w = -np.log2(np.maximum(frac, 1.0 / len(names))); say("%d promoters scanned" % len(names))
score = (P * w).sum(1); nsh = P.sum(1)
rank = score.argsort().argsort(); pctile = 100.0 * rank / max(1, len(score) - 1); gi = {g: i for i, g in enumerate(names)}
rows = []
for g, i in gi.items():
    if g.lower() == GENE.lower(): continue
    sh = [pm2[k][1] for k in np.where(P[i])[0]]; d = D.get(g.lower(), {})
    rows.append(dict(gene=g, shared_score=score[i], n_shared=int(nsh[i]), percentile=pctile[i], baseMean=d.get("bm", np.nan), log2FC=d.get("lfc", np.nan), lfcSE=d.get("se", np.nan), padj=d.get("padj", np.nan), shared_TFs=";".join(sh)))
rows.sort(key=lambda r: -r["shared_score"]); keys = list(rows[0].keys())
with open(os.path.join(OUT, "top1_shared_tf_genes.tsv"), "w") as o:
    o.write("\t".join(keys) + "\n"); [o.write("\t".join(("%.4g" % r[k]) if isinstance(r[k], float) else str(r[k]) for k in keys) + "\n") for r in rows]
with open(os.path.join(OUT, "top1_motifs_used.tsv"), "w") as o:
    o.write("motif_id\tTF\ttop1_score\tpct_promoters_with_motif\tweight\tTF_baseMean\n")
    for j, k in enumerate(keep): o.write("%s\t%s\t%.3f\t%.1f\t%.2f\t%s\n" % (pm[k][0], pm[k][1], s_top[k], 100 * frac[j], w[j], "%.0f" % tf_expr(pm[k][1])[0] if np.isfinite(tf_expr(pm[k][1])[0]) else "NA"))
say("\ntop genes sharing %s promoter motifs:" % GENE)
for r in rows[:TOP_N]: say("  %-14s score %6.1f  shared %3d  log2FC %+5.2f  padj %-8s %s" % (r["gene"], r["shared_score"], r["n_shared"], r["log2FC"], "%.2g" % r["padj"] if np.isfinite(r["padj"]) else "NA", r["shared_TFs"][:60]))

def stars(p): return "" if not np.isfinite(p) else "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else "ns"
top = [r for r in rows[:TOP_N] if np.isfinite(r["log2FC"])]; n = len(top); x = np.arange(n)
fig, axs = plt.subplots(2, 1, figsize=(max(7, 0.45 * n + 3), 7.2), sharex=True, gridspec_kw=dict(height_ratios=[1, 2.2], hspace=0.08))
axs[0].bar(x, [r["shared_score"] for r in top], color="#9AA5B1", width=0.7); axs[0].set_ylabel("shared-motif score\n(rarity-weighted)", fontsize=8, color=INK)
for xi, r in zip(x, top): axs[0].text(xi, r["shared_score"], str(r["n_shared"]), ha="center", va="bottom", fontsize=6, color=INK)
axs[0].set_title("Genes whose promoters share transcription-factor motifs with %s: %s differential expression\n(top panel: score and number of shared motifs; bottom: shrunk log2FC +/- 1 SE, padj stars)" % (GENE, COMP.replace("_vs_", " vs ").replace("_", "+")), fontsize=9.5, color=INK, fontweight="bold", loc="left")
lf = np.array([r["log2FC"] for r in top]); se = np.array([r["lfcSE"] for r in top]); axs[1].bar(x, lf, color=FILL, width=0.7); axs[1].errorbar(x, lf, yerr=se, fmt="none", ecolor=INK, elinewidth=0.6, capsize=2); axs[1].axhline(0, color=INK, lw=0.6)
for xi, r in zip(x, top): sg = stars(r["padj"]); axs[1].text(xi, r["log2FC"] + (1 if r["log2FC"] >= 0 else -1) * (r["lfcSE"] + 0.03) if np.isfinite(r["lfcSE"]) else r["log2FC"], sg, ha="center", va="bottom" if r["log2FC"] >= 0 else "top", fontsize=7, color=INK)
axs[1].set_xticks(x); axs[1].set_xticklabels([r["gene"] + ("†" if r["baseMean"] < 50 else "") for r in top], rotation=60, ha="right", fontsize=8, style="italic"); axs[1].set_ylabel("log2 fold-change (%s)" % COMP.replace("_vs_", " vs ").replace("_", "+"), fontsize=8, color=INK)
for a in axs: a.tick_params(colors=INK, labelsize=8); [sp.set_color("#9AA5B1") for sp in a.spines.values()]
fig.text(0.01, 0.005, "dagger = baseMean < 50.  Shared motif = JASPAR matrix scoring >= %.2f in both promoters; weight = -log2(share of all promoters with the motif).  A shared motif is possible co-regulation, not evidence of it." % REL, fontsize=6.5, color="grey")
for e in ("png", "pdf"): fig.savefig(os.path.join(OUT, "bars_top1_shared_tf_%s.%s" % (COMP, e)), dpi=300, bbox_inches="tight")
plt.close(fig)
# motif map: top genes x the motifs most shared
cols = [j for j in np.argsort(-P[[gi[r["gene"]] for r in top]].sum(0)) if P[[gi[r["gene"]] for r in top], j].sum() > 0][:40]
if cols:
    Sm = np.array([[S[gi[r["gene"]], j] for j in cols] for r in top]); fig, ax = plt.subplots(figsize=(0.32 * len(cols) + 4, 0.26 * n + 2.2)); Sm2 = np.where(Sm >= REL, Sm, np.nan)
    im = ax.imshow(np.ma.masked_invalid(Sm2), aspect="auto", cmap="viridis", vmin=REL, vmax=1.0, interpolation="nearest"); ax.set_xticks(range(len(cols))); ax.set_xticklabels([pm2[j][1] for j in cols], rotation=75, ha="right", fontsize=7)
    ax.set_yticks(range(n)); ax.set_yticklabels([r["gene"] for r in top], fontsize=7.5, style="italic"); fig.colorbar(im, ax=ax, label="best relative score (blank = no site >= %.2f)" % REL, fraction=0.04)
    ax.set_title("Promoter motifs shared with %s (rows = top genes; columns = the most-shared %s-promoter motifs)" % (GENE, GENE), fontsize=8.5, color=INK, loc="left"); fig.tight_layout()
    for e in ("png", "pdf"): fig.savefig(os.path.join(OUT, "top1_shared_motif_map.%s" % e), dpi=300, bbox_inches="tight")
    plt.close(fig)
say("\nWrote results to %s" % OUT)
