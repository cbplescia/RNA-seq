#!/usr/bin/env python3
"""Are the construct sequences identical to the reference in every sample, and if not, where and how often do they differ?   (needs the pileups from 112_construct_variants_align.sh)

Per sample and per construct feature (each record of ref/recombinant_model.fa = one CDS / transcribed unit) the script counts, at every base, how many reads show the reference base, another base (A/C/G/T), a deletion or an insertion.
  identity      = 100 * (1 - non-reference bases / all bases) over positions with depth >= MIN_DEPTH   (this includes sequencing error, ~0.1-0.5 %, so 100.00 % is never expected)
  variant site  = a position with depth >= MIN_DEPTH, >= MIN_ALT reads and a fraction >= VAF_MIN for the same non-reference allele (and, with STRAND_FILTER=1, seen on both strands) - i.e. above the error floor
  shared        = the same variant is called in >= SHARED_FRAC of the samples that have coverage there; a shared variant at ~100 % means the reference file differs from the sequence that was really expressed (or transfected);
                  a variant in only some samples / replicates at an intermediate fraction is genuine heterogeneity (mixed population, RNA editing, mutation) or an artefact of low depth.
Caveats: RT / PCR / sequencing error and read-mapping artefacts all look like low-frequency variants; the thresholds are what separate them. Variants of a few % in a few reads are NOT reliable. Samples without construct (Mock) show up as 'no coverage'.
Outputs (PROJ/construct_variants/results/): feature_sample_summary.tsv, variants.tsv (every called variant, with VAF in each sample and the amino-acid change), per_base_counts_<feature>.tsv, and figures
  construct_identity_overview, variants_<feature> (depth track + heatmap of the non-reference fraction per sample + variant lollipops), construct_variant_spectrum, construct_variant_frequencies  (.png / .pdf).
Options (environment): MIN_DEPTH=30 MIN_ALT=5 VAF_MIN=0.02 STRAND_FILTER=1 SHARED_FRAC=0.5 N_LABEL=12 PROJ=...      Usage: python 113_construct_variants.py     (standard library + numpy, pandas, matplotlib)
"""
import os, sys, glob, re, gzip
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, ListedColormap
from matplotlib.patches import Patch

PROJ = os.environ.get("PROJ", os.path.join(os.path.expanduser("~"), "Documents/aviti_rnaseq/2025/aviti/20251009-100725_Maury_25161_WivfFacICOTHFAUJwSPeHBEFfCfNtCRpUSgVGUkz/results/Maury_25161"))
W = os.path.join(PROJ, "construct_variants"); OUT = os.path.join(W, "results"); os.makedirs(OUT, exist_ok=True)
MIN_DEPTH = int(os.environ.get("MIN_DEPTH", "30")); MIN_ALT = int(os.environ.get("MIN_ALT", "5")); VAF_MIN = float(os.environ.get("VAF_MIN", "0.02"))
STRAND_FILTER = os.environ.get("STRAND_FILTER", "1") == "1"; SHARED_FRAC = float(os.environ.get("SHARED_FRAC", "0.5")); N_LABEL = int(os.environ.get("N_LABEL", "12"))
INK = "#1F2933"
ORDER = ["Mock_Mock", "IFNB_Mock", "IFNY_Mock", "Mock_EBOV", "IFNB_EBOV", "IFNY_EBOV", "Mock_MARV", "IFNB_MARV", "IFNY_MARV", "Mock_EVLP", "Mock_MVLP"]
LAB = {"Mock_Mock": "Mock", "IFNB_Mock": "IFNB", "IFNY_Mock": "IFNG", "Mock_EBOV": "EBOV", "IFNB_EBOV": "IFNB+EBOV", "IFNY_EBOV": "IFNG+EBOV", "Mock_MARV": "MARV", "IFNB_MARV": "IFNB+MARV", "IFNY_MARV": "IFNG+MARV", "Mock_EVLP": "EVLP", "Mock_MVLP": "MVLP"}
BASES = "ACGT"; ALTS = ["A", "C", "G", "T", "del", "ins"]
CODON = {}
_aa = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
for i, a in enumerate("TCAG"):
    for j, b in enumerate("TCAG"):
        for k, c in enumerate("TCAG"): CODON[a + b + c] = _aa[16 * i + 4 * j + k]

def read_fasta(path):
    d, n = {}, None
    for line in open(path):
        line = line.strip()
        if line.startswith(">"): n = line[1:].split()[0]; d[n] = []
        elif n: d[n].append(line.upper())
    return {k: "".join(v) for k, v in d.items()}

def parse_pileup(path, refs):
    """-> {feature: (counts[L,6], fwd[L,6])}; columns A,C,G,T,del,ins"""
    C = {f: np.zeros((len(s), 6), dtype=np.int64) for f, s in refs.items()}; F = {f: np.zeros((len(s), 6), dtype=np.int64) for f, s in refs.items()}
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as fh:
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) < 5 or p[0] not in C: continue
            pos = int(p[1]) - 1; ref = p[2].upper(); bs = p[4]
            if int(p[3]) == 0 or ref not in BASES: continue
            ri = BASES.index(ref); c, f = C[p[0]][pos], F[p[0]][pos]; i, n = 0, len(bs)
            while i < n:
                ch = bs[i]
                if ch == "^": i += 2; continue
                if ch == "$": i += 1; continue
                if ch in ".,": c[ri] += 1; f[ri] += (ch == "."); i += 1; continue
                if ch in "ACGTacgt": k = BASES.index(ch.upper()); c[k] += 1; f[k] += ch.isupper(); i += 1; continue
                if ch in "*#": c[4] += 1; f[4] += (ch == "*"); i += 1; continue
                if ch in "+-":
                    j = i + 1
                    while j < n and bs[j].isdigit(): j += 1
                    ln = int(bs[i + 1:j])
                    if ch == "+": c[5] += 1; f[5] += bs[j:j + ln].isupper()
                    i = j + ln; continue
                i += 1     # N, n and anything else
    return C, F

# ---- samples ------------------------------------------------------------------------------------------------------
refs = read_fasta(os.path.join(PROJ, "ref", "recombinant_model.fa")); feats = list(refs)
pus = sorted(glob.glob(os.path.join(W, "pileup", "*.pileup.txt*")))
if not pus: sys.exit("No pileups in %s - run 112_construct_variants_align.sh first" % os.path.join(W, "pileup"))
si = pd.read_csv(os.path.join(PROJ, "sample_info.csv"), sep=None, engine="python", dtype=str)
si["ifn"] = si["condition1"].fillna("Mock").replace({"None": "Mock", "none": "Mock"}); si["virus"] = si["condition2"].fillna("Mock").replace({"None": "Mock", "none": "Mock"}); si["group"] = si["ifn"] + "_" + si["virus"]
grp = dict(zip(si["sample_id"], si["group"]))
samples = [os.path.basename(p).split(".pileup")[0] for p in pus]
samples = sorted(samples, key=lambda s: (ORDER.index(grp.get(s, "")) if grp.get(s, "") in ORDER else 99, s)); path_of = {os.path.basename(p).split(".pileup")[0]: p for p in pus}
cond = {s: LAB.get(grp.get(s, "?"), grp.get(s, "?")) for s in samples}
conds = list(dict.fromkeys(cond[s] for s in samples)); cpal = dict(zip(conds, plt.cm.tab20(np.linspace(0, 1, 20))[[0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 1, 3][:len(conds)]] if len(conds) <= 12 else plt.cm.tab20(np.linspace(0, 1, len(conds)))))
print("%d samples, %d features: %s" % (len(samples), len(feats), ", ".join("%s(%d nt)" % (f, len(refs[f])) for f in feats)))
DATA = {}
for s in samples:
    print("  reading", s, flush=True); DATA[s] = parse_pileup(path_of[s], refs)

def depth_of(c): return c[:, :5].sum(axis=1)          # A,C,G,T,del
def nonref(c, seq): ri = np.array([BASES.index(b) if b in BASES else 0 for b in seq]); d = c[:, :5].copy(); d[np.arange(len(seq)), ri] = 0; return d

# ---- per sample x feature summary + variant calling -------------------------------------------------------------------
rows, vrows = [], []
for f in feats:
    seq = refs[f]; L = len(seq)
    for s in samples:
        c, fw = DATA[s][0][f], DATA[s][1][f]; dp = depth_of(c); ok = dp >= MIN_DEPTH; nr = nonref(c, seq)
        vaf = np.where(dp[:, None] > 0, c[:, :5] / np.maximum(dp[:, None], 1), 0.0); nrf = nonref(c, seq)
        calls = np.zeros(L, dtype=bool)
        for a in range(5):
            rev = c[:, a] - fw[:, a]; sok = ((fw[:, a] >= 1) & (rev >= 1)) if STRAND_FILTER else True
            m = ok & (nrf[:, a] >= MIN_ALT) & (nrf[:, a] / np.maximum(dp, 1) >= VAF_MIN) & sok
            for pos in np.where(m)[0]: vrows.append((f, pos + 1, seq[pos], ALTS[a], s, int(nrf[pos, a]), int(dp[pos]), nrf[pos, a] / dp[pos])); calls[pos] = True
        tot = int(dp[ok].sum()); mis = int(nr[ok].sum()); nins = int(c[ok, 5].sum())
        rows.append(dict(feature=f, length=L, sample=s, condition=cond[s], median_depth=float(np.median(dp)), pct_positions_covered=100 * ok.mean(), positions_covered=int(ok.sum()),
                         identity_pct=(100 * (1 - mis / tot)) if tot else np.nan, mismatch_bases=mis, total_bases=tot, insertion_events=nins, variant_sites=int(calls.sum()),
                         pct_covered_positions_with_no_variant=(100 * (1 - calls[ok].sum() / ok.sum())) if ok.sum() else np.nan))
S = pd.DataFrame(rows); S.to_csv(os.path.join(OUT, "feature_sample_summary.tsv"), sep="\t", index=False, float_format="%.4f")
V = pd.DataFrame(vrows, columns=["feature", "pos", "ref", "alt", "sample", "alt_reads", "depth", "vaf"])
# amino-acid consequence for substitutions inside a CDS (record starts with ATG, length multiple of 3)
def aa_change(f, pos, ref, alt):
    seq = refs[f]
    if len(seq) % 3 or not seq.startswith("ATG") or alt not in BASES: return ""
    cod = (pos - 1) // 3; i = (pos - 1) % 3; c0 = seq[3 * cod:3 * cod + 3]; c1 = c0[:i] + alt + c0[i + 1:]
    a0, a1 = CODON.get(c0, "?"), CODON.get(c1, "?")
    return "%s%d%s (%s)" % (a0, cod + 1, a1, "silent" if a0 == a1 else "nonsense" if a1 == "*" else "missense")
cover = {f: {s: depth_of(DATA[s][0][f]) >= MIN_DEPTH for s in samples} for f in feats}
if len(V):
    key = ["feature", "pos", "ref", "alt"]
    agg = V.groupby(key).agg(n_samples=("sample", "nunique"), mean_vaf=("vaf", "mean"), max_vaf=("vaf", "max"), min_vaf=("vaf", "min")).reset_index()
    agg["n_covered"] = [sum(cover[r.feature][s][r.pos - 1] for s in samples) for r in agg.itertuples()]
    agg["frac_samples"] = agg.n_samples / agg.n_covered.clip(lower=1)
    agg["aa_change"] = [aa_change(r.feature, r.pos, r.ref, r.alt) for r in agg.itertuples()]
    agg["class"] = np.where((agg.frac_samples >= SHARED_FRAC) & (agg.mean_vaf >= 0.9), "reference differs (shared, ~100%)", np.where(agg.frac_samples >= SHARED_FRAC, "shared, partial", "sample-specific"))
    wide = V.pivot_table(index=key, columns="sample", values="vaf").reindex(columns=samples).reset_index()
    out = agg.merge(wide, on=key).sort_values(["feature", "pos"]); out.to_csv(os.path.join(OUT, "variants.tsv"), sep="\t", index=False, float_format="%.4f")
else:
    agg = pd.DataFrame(columns=["feature", "pos", "ref", "alt", "n_samples", "mean_vaf", "max_vaf", "n_covered", "frac_samples", "aa_change", "class"])
    pd.DataFrame().to_csv(os.path.join(OUT, "variants.tsv"), sep="\t", index=False)
print("\nPer feature (samples with any construct coverage):")
for f in feats:
    sub = S[(S.feature == f) & (S.positions_covered > 0)]
    print("  %-22s %2d samples covered, identity %.3f-%.3f %%, variant sites per sample median %g, distinct variants %d" % (f, len(sub), sub.identity_pct.min() if len(sub) else np.nan, sub.identity_pct.max() if len(sub) else np.nan, sub.variant_sites.median() if len(sub) else 0, int((agg.feature == f).sum())))
if len(agg):
    print("\nVariant classes:"); print(agg["class"].value_counts().to_string())

def save(fig, name): fig.savefig(os.path.join(OUT, name + ".png"), dpi=300, bbox_inches="tight"); fig.savefig(os.path.join(OUT, name + ".pdf"), bbox_inches="tight"); plt.close(fig)
def style(ax): ax.tick_params(colors=INK, labelsize=8); [sp.set_color("#9AA5B1") for sp in ax.spines.values()]

# ---- figure 1: identity overview (sample x feature) --------------------------------------------------------------------
nS, nF = len(samples), len(feats)
fig, axs = plt.subplots(1, 3, figsize=(max(9, 0.55 * nF * 3 + 6), max(5, 0.28 * nS + 2.5)), sharey=True)
panels = [("identity_pct", "Identity to reference (%)\nover covered positions", "viridis", None), ("variant_sites", "Called variant sites\n(above error floor)", "magma_r", None), ("median_depth", "Median depth\n(reads per base)", "Blues", "log")]
for ax, (col, ttl, cm, norm) in zip(axs, panels):
    M = np.full((nS, nF), np.nan)
    for i, s in enumerate(samples):
        for j, f in enumerate(feats):
            r = S[(S["sample"] == s) & (S.feature == f)].iloc[0]
            if r.positions_covered > 0 or col == "median_depth": M[i, j] = r[col]
    Mm = np.ma.masked_invalid(M); cmap = plt.get_cmap(cm).copy(); cmap.set_bad("#E5E7EB")
    if norm == "log": im = ax.imshow(np.ma.masked_less_equal(Mm, 0), aspect="auto", cmap=cmap, norm=LogNorm(vmin=1, vmax=max(10, np.nanmax(M[np.isfinite(M)]) if np.isfinite(M).any() else 10)))
    elif col == "identity_pct": im = ax.imshow(Mm, aspect="auto", cmap=cmap, vmin=min(99, np.nanmin(M) if np.isfinite(M).any() else 99), vmax=100)
    else: im = ax.imshow(Mm, aspect="auto", cmap=cmap, vmin=0)
    for i in range(nS):
        for j in range(nF):
            if np.isfinite(M[i, j]) and (col != "median_depth" or M[i, j] >= 1): ax.text(j, i, ("%.2f" % M[i, j]) if col == "identity_pct" else ("%d" % M[i, j]), ha="center", va="center", fontsize=6, color="white" if (col == "identity_pct" and im.norm(M[i, j]) < 0.55) or (col == "variant_sites" and im.norm(M[i, j]) > 0.5) or (col == "median_depth" and im.norm(M[i, j]) > 0.45) else INK)
    ax.set_xticks(range(nF)); ax.set_xticklabels(feats, rotation=60, ha="right", fontsize=8, color=INK); ax.set_title(ttl, fontsize=9, color=INK); style(ax)
    cb = fig.colorbar(im, ax=ax, fraction=0.04, pad=0.02); cb.ax.tick_params(labelsize=7)
axs[0].set_yticks(range(nS)); axs[0].set_yticklabels(["%s  %s" % (cond[s], s) for s in samples], fontsize=7, color=INK)
fig.suptitle("Are the construct sequences identical to the reference?  (grey = fewer than %d reads at every position, e.g. no construct)" % MIN_DEPTH, fontsize=10, color=INK, fontweight="bold")
fig.tight_layout(); save(fig, "construct_identity_overview")

# ---- figure 2: per feature landscape --------------------------------------------------------------------------------------
vcm = plt.get_cmap("YlOrRd").copy(); vcm.set_bad("#E5E7EB"); vcm.set_under("#F8F9FA")
for f in feats:
    seq = refs[f]; L = len(seq); cs = [s for s in samples if cover[f][s].any()]
    if not cs: print("  %s: no sample has >= %d reads at any position - skipping figure" % (f, MIN_DEPTH)); continue
    table = pd.DataFrame({"pos": np.arange(1, L + 1), "ref": list(seq)})
    for s in cs:
        c = DATA[s][0][f]; d = depth_of(c); table[s + "_depth"] = d; table[s + "_nonref_frac"] = np.where(d > 0, nonref(c, seq).sum(axis=1) / np.maximum(d, 1), np.nan)
    table.to_csv(os.path.join(OUT, "per_base_counts_%s.tsv" % f), sep="\t", index=False, float_format="%.5f")
    fig = plt.figure(figsize=(12, 3.0 + 0.2 * len(cs) + 2.2)); gs = fig.add_gridspec(3, 2, height_ratios=[1.3, 0.22 * len(cs) + 0.8, 1.5], width_ratios=[40, 1], hspace=0.12, wspace=0.02)
    a1 = fig.add_subplot(gs[0, 0]); a2 = fig.add_subplot(gs[1, 0], sharex=a1); a3 = fig.add_subplot(gs[2, 0], sharex=a1); cax = fig.add_subplot(gs[1, 1])
    x = np.arange(1, L + 1)
    for cd in dict.fromkeys(cond[s] for s in cs):
        ds = np.array([depth_of(DATA[s][0][f]) for s in cs if cond[s] == cd]); a1.plot(x, np.median(ds, axis=0) + 1, color=cpal[cd], lw=1.0, label=cd)
    a1.set_yscale("log"); a1.axhline(MIN_DEPTH, color="grey", lw=0.6, ls=":"); a1.set_ylabel("depth (median\nper condition)", fontsize=8, color=INK); a1.legend(fontsize=6, ncol=min(6, len(conds)), frameon=False, loc="upper center", bbox_to_anchor=(0.5, 1.5)); plt.setp(a1.get_xticklabels(), visible=False); style(a1)
    H = np.array([np.where(cover[f][s], np.maximum(nonref(DATA[s][0][f], seq).sum(axis=1) / np.maximum(depth_of(DATA[s][0][f]), 1), 1e-5), np.nan) for s in cs])
    im = a2.imshow(np.ma.masked_invalid(H), aspect="auto", cmap=vcm, norm=LogNorm(vmin=1e-3, vmax=1), extent=(0.5, L + 0.5, len(cs) - 0.5, -0.5), interpolation="nearest")
    a2.set_yticks(range(len(cs))); a2.set_yticklabels(["%s %s" % (cond[s], s) for s in cs], fontsize=6, color=INK); plt.setp(a2.get_xticklabels(), visible=False); style(a2)
    for i, s in enumerate(cs): a2.add_patch(plt.Rectangle((-0.0, i - 0.5), -L * 0.006, 1, color=cpal[cond[s]], clip_on=False, transform=a2.transData))
    cb = fig.colorbar(im, cax=cax); cb.set_label("non-reference fraction", fontsize=7); cb.ax.tick_params(labelsize=6)
    sub = agg[agg.feature == f].copy()
    if len(sub):
        cm_ = {"reference differs (shared, ~100%)": "#C0392B", "shared, partial": "#E08E0B", "sample-specific": "#2B6CB0"}
        for k, g in sub.groupby("class"): a3.vlines(g.pos, VAF_MIN * 50, g.mean_vaf * 100, color=cm_[k], lw=0.8); a3.scatter(g.pos, g.mean_vaf * 100, s=8 + 60 * g.frac_samples, color=cm_[k], edgecolor="white", linewidth=0.4, zorder=3, label=k)
        top = sub.sort_values(["n_samples", "mean_vaf"], ascending=False).head(N_LABEL)
        for r in top.itertuples(): a3.annotate("%s%d%s%s" % (r.ref, r.pos, r.alt, (" " + r.aa_change.split(" ")[0]) if r.aa_change else ""), (r.pos, r.mean_vaf * 100), xytext=(0, 5), textcoords="offset points", ha="center", fontsize=6, color=INK, rotation=60)
        a3.set_yscale("log"); a3.set_ylim(VAF_MIN * 50, 600); a3.legend(fontsize=6, frameon=False, loc="upper right", ncol=3)
    else: a3.text(0.5, 0.5, "no variant above the error floor (>= %d reads and >= %.1f%%)" % (MIN_ALT, 100 * VAF_MIN), ha="center", va="center", transform=a3.transAxes, fontsize=9, color=INK)
    a3.set_ylabel("mean frequency in the\nsamples with the variant (%)", fontsize=8, color=INK); a3.set_xlabel("position in %s (nt)   [point size = share of covered samples carrying the variant]" % f, fontsize=8, color=INK); style(a3)
    a1.set_xlim(0.5, L + 0.5); a1.set_title("%s (%d nt): coverage, non-reference fraction per sample, and called variants" % (f, L), fontsize=10, color=INK, fontweight="bold", loc="left", pad=22)
    fig.legend(handles=[Patch(color=cpal[c], label=c) for c in dict.fromkeys(cond[s] for s in cs)], fontsize=6, ncol=6, frameon=False, loc="lower center", bbox_to_anchor=(0.5, -0.04), title="row labels = condition + sample", title_fontsize=6)
    save(fig, "variants_%s" % f)

# ---- figure 3: substitution spectrum by condition (rate per reference base) -------------------------------------------------
types = [r + ">" + a for r in BASES for a in BASES if r != a]; R = pd.DataFrame(0.0, index=list(dict.fromkeys(cond[s] for s in samples)), columns=types); Dn = R.copy()
for s in samples:
    for f in feats:
        c = DATA[s][0][f]; seq = refs[f]; ok = depth_of(c) >= MIN_DEPTH
        for r in BASES:
            m = ok & np.array([b == r for b in seq]); tot = c[m, :4].sum()
            for k, a in enumerate(BASES):
                if a != r: R.loc[cond[s], r + ">" + a] += c[m, k].sum(); Dn.loc[cond[s], r + ">" + a] += tot
rate = (R / Dn.replace(0, np.nan)) * 100
if rate.notna().any().any():
    fig, ax = plt.subplots(figsize=(8.5, 0.45 * len(rate) + 2.2)); im = ax.imshow(np.ma.masked_invalid(rate.values), aspect="auto", cmap="magma_r", vmin=0)
    ax.set_xticks(range(12)); ax.set_xticklabels(types, fontsize=8, color=INK); ax.set_yticks(range(len(rate))); ax.set_yticklabels(rate.index, fontsize=8, color=INK)
    for i in range(rate.shape[0]):
        for j in range(12):
            if np.isfinite(rate.values[i, j]): ax.text(j, i, "%.3f" % rate.values[i, j], ha="center", va="center", fontsize=6, color="white" if im.norm(rate.values[i, j]) > 0.5 else INK)
    fig.colorbar(im, ax=ax, label="% of reads at that reference base", fraction=0.04, pad=0.02); ax.set_title("Mismatch spectrum of construct reads, by condition\n(all non-reference bases incl. sequencing error; e.g. an A>G or C>T excess would suggest RNA editing, uniform noise suggests sequencing error)", fontsize=8.5, color=INK); style(ax)
    fig.tight_layout(); save(fig, "construct_variant_spectrum"); rate.to_csv(os.path.join(OUT, "mismatch_spectrum_by_condition.tsv"), sep="\t", float_format="%.5f")

# ---- figure 4: frequency distribution of calls ------------------------------------------------------------------------------
if len(V):
    fig, axs = plt.subplots(1, 2, figsize=(10, 3.8)); bins = np.logspace(np.log10(VAF_MIN), 0, 25)
    for k, col in [("reference differs (shared, ~100%)", "#C0392B"), ("shared, partial", "#E08E0B"), ("sample-specific", "#2B6CB0")]:
        ks = agg[agg["class"] == k]; vv = V.merge(ks[["feature", "pos", "ref", "alt"]], on=["feature", "pos", "ref", "alt"])
        if len(vv): axs[0].hist(vv.vaf, bins=bins, color=col, alpha=0.8, label="%s (%d variants)" % (k, len(ks)))
    axs[0].set_xscale("log"); axs[0].set_xlabel("variant frequency in a sample", fontsize=8, color=INK); axs[0].set_ylabel("sample-level calls", fontsize=8, color=INK); axs[0].legend(fontsize=6, frameon=False); style(axs[0])
    perS = V.groupby(["sample"]).size().reindex(samples).fillna(0)
    axs[1].bar(range(nS), perS.values, color=[cpal[cond[s]] for s in samples]); axs[1].set_xticks(range(nS)); axs[1].set_xticklabels(["%s %s" % (cond[s], s) for s in samples], rotation=70, ha="right", fontsize=6); axs[1].set_ylabel("variant calls (all features)", fontsize=8, color=INK); style(axs[1])
    fig.suptitle("Called construct variants: frequency distribution and number per sample", fontsize=10, color=INK, fontweight="bold"); fig.tight_layout(); save(fig, "construct_variant_frequencies")
print("\nWrote results to", OUT)
