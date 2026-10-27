#!/usr/bin/env python3
"""Heatmap of the EBOV amino-acid changes (from 117_ebov_aa_mutations_ifn.py): the untreated EBOV samples are GREY (their own frequency printed in the tile), the IFNB+EBOV and IFNG+EBOV samples are coloured by their CHANGE from the untreated EBOV mean
(red = more abundant than untreated EBOV, blue = less, GREY = no change - the same grey as the untreated tiles; change = frequency in the sample - mean frequency of the untreated EBOV samples, in percentage points; the number in the tile is that change).
One row per amino-acid change (gene + residue change), rows grouped by gene; blank/hatched = sample has < MIN_DEPTH reads there.  Colour scale is symmetric and compressive (symlog), so a +30 point change and a +1 point change can be seen together.
Reads construct_variants/results/ebov_aa/ebov_aa_mutations.tsv (run 117 first).   Writes ebov_aa_heatmap_vs_untreated_pN.png/.pdf (all genes, paged), ebov_aa_heatmap_<gene>.png/.pdf, and ebov_aa_delta_table.tsv in the same folder.
Options (environment): AVERAGE=1 (average the replicates of each condition first: three columns, EBOV | IFNB+EBOV | IFNG+EBOV; mean frequency of the covered replicates; files are named ebov_aa_heatmap_mean_*)  MIN_DELTA=0 (keep rows where some IFN sample differs from untreated by >= this many percentage points; e.g. MIN_DELTA=2)  MAX_ROWS=40  GENES=model_GP,model_NP (only these)  SORT=position|delta  PER_GENE=1  PROJ=...
Usage: python 118_ebov_aa_heatmap_vs_untreated.py     (numpy, pandas, matplotlib)
"""
import os, sys
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import SymLogNorm
from matplotlib.patches import Patch

PROJ = os.environ.get("PROJ", os.path.join(os.path.expanduser("~"), "Documents/aviti_rnaseq/2025/aviti/20251009-100725_Maury_25161_WivfFacICOTHFAUJwSPeHBEFfCfNtCRpUSgVGUkz/results/Maury_25161"))
OUT = os.path.join(PROJ, "construct_variants", "results", "ebov_aa")
MIN_DELTA = float(os.environ.get("MIN_DELTA", "0")); MAX_ROWS = int(os.environ.get("MAX_ROWS", "40")); SORT = os.environ.get("SORT", "position"); PER_GENE = os.environ.get("PER_GENE", "1") == "1"; AVERAGE = os.environ.get("AVERAGE", "0") == "1"; GENES = [g.strip() for g in os.environ.get("GENES", "").split(",") if g.strip()]
INK = "#1F2933"; GC = {"Untreated": "#6B7280", "IFNB": "#2B6CB0", "IFNG": "#D97706"}
tp = os.path.join(OUT, "ebov_aa_mutations.tsv")
if not os.path.exists(tp): sys.exit("Missing %s - run 117_ebov_aa_mutations_ifn.py first" % tp)
T = pd.read_csv(tp, sep="\t", dtype={"feature": str}); samples = [c[4:] for c in T.columns if c.startswith("vaf_")]
si = pd.read_csv(os.path.join(PROJ, "sample_info.csv"), sep=None, engine="python", dtype=str)
ifn = dict(zip(si["sample_id"], si["condition1"].fillna("None"))); gmap = {s: ("IFNB" if ifn.get(s) == "IFNB" else "IFNG" if ifn.get(s) in ("IFNY", "IFNG") else "Untreated") for s in samples}
samples = sorted(samples, key=lambda s: (["Untreated", "IFNB", "IFNG"].index(gmap[s]), s)); U = [s for s in samples if gmap[s] == "Untreated"]
if GENES: T = T[T.feature.isin(GENES)]
V = 100 * T[["vaf_" + s for s in samples]].values; mu = np.nanmean(100 * T[["vaf_" + s for s in U]].values, axis=1) if U else np.full(len(T), np.nan)
Dl = V - mu[:, None]
for j, s in enumerate(samples):
    if gmap[s] == "Untreated": Dl[:, j] = np.nan          # untreated columns are grey, not coloured
T = T.assign(mean_untreated_pct=mu, max_abs_delta=np.nanmax(np.abs(Dl), axis=1) if Dl.shape[1] else np.nan)
for j, s in enumerate(samples): T["delta_" + s] = Dl[:, j]
T["mean_delta_IFNB"] = [np.nanmean([Dl[i, j] for j, s in enumerate(samples) if gmap[s] == "IFNB"] or [np.nan]) for i in range(len(T))]; T["mean_delta_IFNG"] = [np.nanmean([Dl[i, j] for j, s in enumerate(samples) if gmap[s] == "IFNG"] or [np.nan]) for i in range(len(T))]
T = T[T.max_abs_delta.fillna(0) >= MIN_DELTA]
if not len(T): sys.exit("No row has a change >= %g percentage points" % MIN_DELTA)
T["gene"] = T.feature.str.replace("model_", "", regex=False)
T = T.sort_values(["feature", "residue"]) if SORT == "position" else T.assign(_s=T[["mean_delta_IFNB", "mean_delta_IFNG"]].mean(axis=1)).sort_values(["feature", "_s"], ascending=[True, False])
T.drop(columns=[c for c in ["_s"] if c in T]).to_csv(os.path.join(OUT, "ebov_aa_delta_table.tsv"), sep="\t", index=False, float_format="%.4g")
Dall = np.array([T["delta_" + s].values for s in samples]).T; lim = max(5.0, float(np.nanpercentile(np.abs(Dall[np.isfinite(Dall)]), 99))) if np.isfinite(Dall).any() else 5.0
norm = SymLogNorm(linthresh=1.0, linscale=0.8, vmin=-lim, vmax=lim, base=10); GREY = "#C9CDD3"      # same grey as the untreated tiles: no change from untreated = grey
cm = matplotlib.colors.LinearSegmentedColormap.from_list("blue_grey_red", ["#2166AC", "#92C5DE", GREY, "#F4A582", "#B2182B"]); cm.set_bad((0, 0, 0, 0))
NAME = {"Untreated": "EBOV", "IFNB": "IFNB+EBOV", "IFNG": "IFNG+EBOV"}
def columns(d):
    """-> labels, groups, Vm (frequency %), Dm (change from untreated mean, points; NaN for untreated columns)"""
    if not AVERAGE:
        return (["%s %s" % (NAME[gmap[s]], s) for s in samples], [gmap[s] for s in samples], np.array([100 * d["vaf_" + s].values for s in samples]).T, np.array([d["delta_" + s].values for s in samples]).T)
    grps = ["Untreated", "IFNB", "IFNG"]; n_of = {g: sum(gmap[s] == g for s in samples) for g in grps}
    Vm = 100 * np.array([d.mean_untreated.values, d.mean_IFNB.values, d.mean_IFNG.values]).T; Dm = np.array([np.full(len(d), np.nan), Vm[:, 1] - Vm[:, 0], Vm[:, 2] - Vm[:, 0]]).T
    return (["%s (mean of %d)" % (NAME[g], n_of[g]) for g in grps], grps, Vm, Dm)
if AVERAGE:
    Dall = np.array([T.mean_IFNB.values * 100 - T.mean_untreated.values * 100, T.mean_IFNG.values * 100 - T.mean_untreated.values * 100]); lim = max(5.0, float(np.nanpercentile(np.abs(Dall[np.isfinite(Dall)]), 99))) if np.isfinite(Dall).any() else 5.0
    norm = SymLogNorm(linthresh=1.0, linscale=0.8, vmin=-lim, vmax=lim, base=10)
def style(ax): ax.tick_params(colors=INK, labelsize=8); [sp.set_color("#9AA5B1") for sp in ax.spines.values()]
def draw(d, name, title):
    n = len(d); labels, groups, Vm, Dm = columns(d); nc = len(labels); edges = [j - 0.5 for j in range(1, nc) if groups[j] != groups[j - 1]]
    fig = plt.figure(figsize=(max(7, (1.1 if AVERAGE else 0.55) * nc + 6.5), 1.7 + 0.27 * n + 1.5)); gs = fig.add_gridspec(1, 3, width_ratios=[max(3.2, (1.1 if AVERAGE else 0.55) * nc), 2.6, 0.25], wspace=0.04)
    ax = fig.add_subplot(gs[0, 0]); tx = fig.add_subplot(gs[0, 1]); cax = fig.add_subplot(gs[0, 2])
    G = np.full(Dm.shape, np.nan); G[:, [j for j, g in enumerate(groups) if g == "Untreated"]] = 1.0
    ax.imshow(np.ma.masked_invalid(G), aspect="auto", cmap=plt.matplotlib.colors.ListedColormap([GREY]), vmin=0, vmax=1, interpolation="nearest")
    im = ax.imshow(np.ma.masked_invalid(Dm), aspect="auto", cmap=cm, norm=norm, interpolation="nearest")
    for i in range(n):
        for j in range(nc):
            if not np.isfinite(Vm[i, j]): ax.add_patch(plt.Rectangle((j - 0.5, i - 0.5), 1, 1, fill=False, hatch="////", edgecolor="#D1D5DB", lw=0)); continue
            if groups[j] == "Untreated": ax.text(j, i, "%.1f" % Vm[i, j], ha="center", va="center", fontsize=7 if AVERAGE else 6, color=INK)
            else: ax.text(j, i, "%+.1f" % Dm[i, j], ha="center", va="center", fontsize=7 if AVERAGE else 6, color="white" if abs(Dm[i, j]) > 6 else INK)
    ax.set_yticks(range(n)); ax.set_yticklabels(["%s %s" % (r.gene, r.label) for r in d.itertuples()], fontsize=7.5, color=INK)
    ax.set_xticks(range(nc)); ax.set_xticklabels(labels, rotation=75 if not AVERAGE else 30, ha="right", fontsize=7 if not AVERAGE else 8)
    for tl, g in zip(ax.get_xticklabels(), groups): tl.set_color(GC[g]); tl.set_fontweight("bold")
    for e in edges: ax.axvline(e, color=INK, lw=1.2)
    gl = list(d.feature); [ax.axhline(i - 0.5, color=INK, lw=1.0) for i in range(1, n) if gl[i] != gl[i - 1]]
    ax.set_xticks(np.arange(-0.5, nc, 1), minor=True); ax.set_yticks(np.arange(-0.5, n, 1), minor=True); ax.grid(which="minor", color="white", lw=0.5); ax.tick_params(which="minor", length=0); style(ax)
    ax.set_title(title, fontsize=10, color=INK, fontweight="bold", loc="left", pad=10)
    tx.set_xlim(0, 3); tx.set_ylim(n - 0.5, -0.5); tx.axis("off")
    for x, h in zip([0.45, 1.5, 2.55], ["untreated\nmean %", "IFNB mean\nchange (pts)", "IFNG mean\nchange (pts)"]): tx.text(x, -0.9, h, ha="center", va="bottom", fontsize=6.5, color=INK, fontweight="bold")
    for i, r in enumerate(d.itertuples()):
        for x, v, f in zip([0.45, 1.5, 2.55], [r.mean_untreated_pct, r.mean_delta_IFNB, r.mean_delta_IFNG], ["%.2f", "%+.2f", "%+.2f"]): tx.text(x, i, "n/a" if not np.isfinite(v) else f % v, ha="center", va="center", fontsize=6.5, color=INK)
    cb = fig.colorbar(im, cax=cax); cb.set_label("change from untreated EBOV mean (percentage points)", fontsize=7); cb.ax.tick_params(labelsize=6)
    ax.set_xlabel("grey tiles = untreated EBOV (%s frequency, %%); coloured tiles = change from the untreated mean (red up, blue down); hatched = below the depth threshold" % ("mean" if AVERAGE else "own"), fontsize=6.5, color="grey", labelpad=6)
    fig.savefig(os.path.join(OUT, name + ".png"), dpi=300, bbox_inches="tight"); fig.savefig(os.path.join(OUT, name + ".pdf"), bbox_inches="tight"); plt.close(fig)
PFX = "ebov_aa_heatmap_mean" if AVERAGE else "ebov_aa_heatmap"
pages = [T.iloc[i:i + MAX_ROWS] for i in range(0, len(T), MAX_ROWS)]
for pi, d in enumerate(pages): draw(d, PFX + "_vs_untreated_p%d" % (pi + 1), "EBOV amino-acid changes%s: change with IFN relative to untreated EBOV%s" % (" (replicates averaged)" if AVERAGE else "", "  (page %d of %d)" % (pi + 1, len(pages)) if len(pages) > 1 else ""))
if PER_GENE:
    for f, d in T.groupby("feature", sort=False):
        for pi, dd in enumerate([d.iloc[i:i + MAX_ROWS] for i in range(0, len(d), MAX_ROWS)]): draw(dd, "%s_%s%s" % (PFX, f, "" if len(d) <= MAX_ROWS else "_p%d" % (pi + 1)), "%s: amino-acid changes, change with IFN relative to untreated EBOV" % f)
print("%d amino-acid changes in %d genes; colour scale +/- %.1f points; wrote figures + ebov_aa_delta_table.tsv to %s" % (len(T), T.feature.nunique(), lim, OUT))
