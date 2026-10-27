#!/usr/bin/env python3
"""Residue-labelled view of the construct variants that are shared by at least MIN_SAMPLES samples, one figure per construct gene (reads results/variants.tsv written by 113_construct_variants.py).
Top panel   : lollipop along the protein (x = residue number); height = number of samples carrying the variant, colour = effect (missense / silent / nonsense / non-coding or indel), label = residue change (e.g. V34I = Val34 -> Ile).
Bottom panel: one row per shared variant, one column per sample (grouped and coloured by condition); tile colour = variant frequency in that sample, blank = not called (below the thresholds of 113, or too little coverage).
Features whose record does not start with ATG / is not a multiple of 3 are shown in nucleotide positions (c.123A>G).   Substitution labels use the reference FASTA reading frame.   Indels have no residue label (c.pos del / ins).
Options (environment): MIN_SAMPLES=2  MAX_ROWS=45 (rows per heatmap page; more variants -> extra pages)  PROJ=...      Writes to construct_variants/results/shared_variants/: shared_variants_<feature>.png/.pdf, shared_variants_all.tsv, shared_variants_overview.png/.pdf
Usage: python 114_construct_shared_variants.py     (numpy, pandas, matplotlib)
"""
import os, sys, re
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.patches import Patch

PROJ = os.environ.get("PROJ", os.path.join(os.path.expanduser("~"), "Documents/aviti_rnaseq/2025/aviti/20251009-100725_Maury_25161_WivfFacICOTHFAUJwSPeHBEFfCfNtCRpUSgVGUkz/results/Maury_25161"))
RES = os.path.join(PROJ, "construct_variants", "results"); OUT = os.path.join(RES, "shared_variants"); os.makedirs(OUT, exist_ok=True)
MIN_SAMPLES = int(os.environ.get("MIN_SAMPLES", "2")); MAX_ROWS = int(os.environ.get("MAX_ROWS", "45"))
INK = "#1F2933"
ORDER = ["Mock_Mock", "IFNB_Mock", "IFNY_Mock", "Mock_EBOV", "IFNB_EBOV", "IFNY_EBOV", "Mock_MARV", "IFNB_MARV", "IFNY_MARV", "Mock_EVLP", "Mock_MVLP"]
LAB = {"Mock_Mock": "Mock", "IFNB_Mock": "IFNB", "IFNY_Mock": "IFNG", "Mock_EBOV": "EBOV", "IFNB_EBOV": "IFNB+EBOV", "IFNY_EBOV": "IFNG+EBOV", "Mock_MARV": "MARV", "IFNB_MARV": "IFNB+MARV", "IFNY_MARV": "IFNG+MARV", "Mock_EVLP": "EVLP", "Mock_MVLP": "MVLP"}
EFF = {"missense": "#C0392B", "silent": "#7F8C8D", "nonsense": "#111111", "non-coding / indel": "#2B6CB0"}

vp = os.path.join(RES, "variants.tsv")
if not os.path.exists(vp) or os.path.getsize(vp) < 5: sys.exit("No variants.tsv - run 113_construct_variants.py first (or no variants were called)")
V = pd.read_csv(vp, sep="\t", dtype={"feature": str}); meta = ["feature", "pos", "ref", "alt", "n_samples", "mean_vaf", "max_vaf", "min_vaf", "n_covered", "frac_samples", "aa_change", "class"]
samples = [c for c in V.columns if c not in meta]
si = pd.read_csv(os.path.join(PROJ, "sample_info.csv"), sep=None, engine="python", dtype=str)
g = (si["condition1"].fillna("Mock").replace({"None": "Mock", "none": "Mock"}) + "_" + si["condition2"].fillna("Mock").replace({"None": "Mock", "none": "Mock"})); grp = dict(zip(si["sample_id"], g))
samples = sorted(samples, key=lambda s: (ORDER.index(grp[s]) if grp.get(s) in ORDER else 99, s)); cond = {s: LAB.get(grp.get(s, "?"), grp.get(s, "?")) for s in samples}
conds = list(dict.fromkeys(cond[s] for s in samples)); cpal = dict(zip(conds, plt.cm.tab20(np.linspace(0, 1, 20))[[0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 1, 3][:len(conds)]] if len(conds) <= 12 else plt.cm.tab20(np.linspace(0, 1, len(conds)))))
refs, n = {}, None
for line in open(os.path.join(PROJ, "ref", "recombinant_model.fa")):
    line = line.strip()
    if line.startswith(">"): n = line[1:].split()[0]; refs[n] = []
    elif n: refs[n].append(line.upper())
refs = {k: "".join(v) for k, v in refs.items()}; is_cds = {k: (len(s) % 3 == 0 and s.startswith("ATG")) for k, s in refs.items()}

def residue(r):
    """-> (x position, label, effect)"""
    m = re.match(r"([A-Z\*\?])(\d+)([A-Z\*\?]) \((\w+)\)", r.aa_change) if isinstance(r.aa_change, str) else None
    if m: return int(m.group(2)), "%s%s%s" % (m.group(1), m.group(2), m.group(3)), m.group(4)
    kind = r.alt if r.alt in ("del", "ins") else "%s>%s" % (r.ref, r.alt)
    x = (r.pos - 1) // 3 + 1 if is_cds[r.feature] else r.pos
    return x, "c.%d%s" % (r.pos, kind), "non-coding / indel"
Sh = V[V.n_samples >= MIN_SAMPLES].copy()
if not len(Sh): sys.exit("No variant is called in >= %d samples" % MIN_SAMPLES)
res = Sh.apply(residue, axis=1, result_type="expand"); Sh["x"], Sh["label"], Sh["effect"] = res[0], res[1], res[2]
Sh["sample_list"] = [", ".join("%s:%s" % (cond[s], s) for s in samples if pd.notna(r[s])) for _, r in Sh.iterrows()]
Sh.sort_values(["feature", "pos"])[["feature", "pos", "x", "label", "ref", "alt", "effect", "n_samples", "n_covered", "mean_vaf", "max_vaf", "class", "sample_list"]].to_csv(os.path.join(OUT, "shared_variants_all.tsv"), sep="\t", index=False, float_format="%.4f")
print("%d variants in >= %d samples across %d features:" % (len(Sh), MIN_SAMPLES, Sh.feature.nunique()))
for f in refs: print("  %-24s %d" % (f, (Sh.feature == f).sum()))
def style(ax): ax.tick_params(colors=INK, labelsize=8); [sp.set_color("#9AA5B1") for sp in ax.spines.values()]
def save(fig, name): fig.savefig(os.path.join(OUT, name + ".png"), dpi=300, bbox_inches="tight"); fig.savefig(os.path.join(OUT, name + ".pdf"), bbox_inches="tight"); plt.close(fig)
cm = plt.get_cmap("YlOrRd").copy(); cm.set_bad("#F3F4F6")

for f, seq in refs.items():
    d = Sh[Sh.feature == f].sort_values("pos")
    if not len(d): continue
    L = len(seq) // 3 if is_cds[f] else len(seq); unit = "residue" if is_cds[f] else "nucleotide"
    pages = [d.iloc[i:i + MAX_ROWS] for i in range(0, len(d), MAX_ROWS)]
    for pi, dp in enumerate(pages):
        fig = plt.figure(figsize=(max(9, 0.32 * len(samples) + 6), 3.4 + 0.24 * len(dp) + 1.2)); gs = fig.add_gridspec(2, 2, height_ratios=[2.3, 0.24 * len(dp) + 0.9], width_ratios=[40, 1], hspace=0.45, wspace=0.03)
        a1 = fig.add_subplot(gs[0, 0]); a2 = fig.add_subplot(gs[1, 0]); cax = fig.add_subplot(gs[1, 1])
        a1.hlines(0, 0, L + 1, color="#9AA5B1", lw=4, zorder=1)
        for r in dp.itertuples(): a1.vlines(r.x, 0, r.n_samples, color=EFF[r.effect], lw=1.0, zorder=2); a1.scatter(r.x, r.n_samples, s=18 + 45 * r.mean_vaf, color=EFF[r.effect], edgecolor="white", linewidth=0.5, zorder=3)
        # label with simple vertical staggering to limit overlaps
        last = {}
        for k, r in enumerate(dp.sort_values("x").itertuples()):
            lvl = k % 3; a1.annotate(r.label, (r.x, r.n_samples), xytext=(0, 6 + 7 * lvl), textcoords="offset points", ha="center", fontsize=6.5, color=INK, rotation=0, arrowprops=dict(arrowstyle="-", lw=0.3, color="#9AA5B1") if lvl else None)
        a1.set_xlim(0, L + 1); a1.set_ylim(0, max(len(samples) * 1.0, dp.n_samples.max()) * 1.18); a1.set_ylabel("samples carrying the variant", fontsize=8, color=INK); a1.set_xlabel("%s number in %s  (%d %ss)   [point size = mean variant frequency]" % (unit, f, L, unit), fontsize=8, color=INK); style(a1)
        a1.legend(handles=[Patch(color=c, label=e) for e, c in EFF.items() if (dp.effect == e).any()], fontsize=7, frameon=False, loc="upper right", ncol=2)
        a1.set_title("%s: variants shared by >= %d samples%s" % (f, MIN_SAMPLES, "  (page %d of %d)" % (pi + 1, len(pages)) if len(pages) > 1 else ""), fontsize=10, color=INK, fontweight="bold", loc="left")
        samples_p = [s for s in samples if dp[s].notna().any()]      # drop samples with no call on this page (e.g. Mock, no construct reads)
        M = np.array([[r[s] if pd.notna(r[s]) else np.nan for s in samples_p] for _, r in dp.iterrows()], dtype=float)
        im = a2.imshow(np.ma.masked_invalid(M), aspect="auto", cmap=cm, norm=LogNorm(vmin=0.02, vmax=1.0), interpolation="nearest")
        a2.set_yticks(range(len(dp))); a2.set_yticklabels(["%s  (%s)" % (r.label, r.effect.split(" ")[0]) for r in dp.itertuples()], fontsize=7, color=INK)
        a2.set_xticks(range(len(samples_p))); a2.set_xticklabels(["%s %s" % (cond[s], s) for s in samples_p], rotation=75, ha="right", fontsize=6.5, color=INK)
        for tl, s in zip(a2.get_xticklabels(), samples_p): tl.set_color(cpal[cond[s]]); tl.set_fontweight("bold")
        for i in range(len(dp)):
            for j in range(len(samples_p)):
                if np.isfinite(M[i, j]): a2.text(j, i, "%d" % round(100 * M[i, j]) if M[i, j] >= 0.1 else "%.0f" % (100 * M[i, j]) if M[i, j] >= 0.095 else "<10", ha="center", va="center", fontsize=5, color="white" if M[i, j] > 0.35 else INK)
        a2.set_xticks(np.arange(-0.5, len(samples_p), 1), minor=True); a2.set_yticks(np.arange(-0.5, len(dp), 1), minor=True); a2.grid(which="minor", color="white", lw=0.6); a2.tick_params(which="minor", length=0); style(a2)
        cb = fig.colorbar(im, cax=cax); cb.set_label("variant frequency", fontsize=7); cb.ax.tick_params(labelsize=6)
        a2.set_title("numbers in tiles = variant frequency (%); blank = not called in that sample", fontsize=7.5, color="grey", loc="left")
        fig.legend(handles=[Patch(color=cpal[c], label=c) for c in conds], fontsize=6, ncol=min(8, len(conds)), frameon=False, loc="lower center", bbox_to_anchor=(0.45, -0.07), title="condition (colour of the sample labels)", title_fontsize=6)
        save(fig, "shared_variants_%s%s" % (f, "" if len(pages) == 1 else "_p%d" % (pi + 1)))

# overview: all features, number of shared variants and how many are missense
cnt = Sh.groupby(["feature", "effect"]).size().unstack(fill_value=0).reindex(list(refs)).fillna(0)
for e in EFF:
    if e not in cnt: cnt[e] = 0
cnt = cnt[list(EFF)]; fig, ax = plt.subplots(figsize=(max(6, 0.6 * len(cnt) + 3), 4)); bot = np.zeros(len(cnt))
for e, c in EFF.items(): ax.bar(range(len(cnt)), cnt[e].values, bottom=bot, color=c, label=e); bot += cnt[e].values
ax.set_xticks(range(len(cnt))); ax.set_xticklabels(cnt.index, rotation=60, ha="right", fontsize=8, color=INK); ax.set_ylabel("variants shared by >= %d samples" % MIN_SAMPLES, fontsize=8, color=INK); ax.legend(fontsize=7, frameon=False); style(ax)
ax.set_title("Shared construct variants per gene", fontsize=10, color=INK, fontweight="bold", loc="left"); fig.tight_layout(); save(fig, "shared_variants_overview")
print("\nWrote figures + shared_variants_all.tsv to", OUT)
