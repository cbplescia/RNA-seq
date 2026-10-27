#!/usr/bin/env python3
"""Amino-acid changes in the EBOV construct genes, and how abundant / enriched each one is with vs without interferon.   EBOV samples only (Mock_EBOV = untreated, IFNB_EBOV, IFNY_EBOV); MARV, VLP and mock-infected samples are ignored.
Works directly on the pileups from 112_construct_variants_align.sh (no threshold is applied to the abundance numbers).
  candidate  : a substitution that, in at least one EBOV sample, has depth >= MIN_DEPTH, >= MIN_ALT reads and frequency >= VAF_MIN for the same alternative base (and, with STRAND_FILTER=1, on both strands)
  abundance  : frequency (alt reads / depth) in every EBOV sample with depth >= MIN_DEPTH, even where it was too low to be 'called'; mean per group = untreated | IFNB | IFNG
  enrichment : log2[(mean frequency in IFN samples + EPS) / (mean frequency in untreated + EPS)]   (IFN = IFNB and IFNG pooled; EPS = 0.2 %, so a rise from 0 to 0.1 % does not look like a huge fold change)
  p-value    : exact permutation test (all ways of relabelling the covered samples) on the difference of the mean frequencies, BH-adjusted across variants. With ~3 untreated and ~6 IFN samples the smallest possible p is ~0.01, so treat p as a ranking, not as proof.
  amino-acid change: from the reading frame of the reference record (starts with ATG, length multiple of 3); changes are labelled W104* (nonsense), V34I (missense) ; silent changes only with SHOW_SILENT=1. Indels and non-coding records are not shown.
Outputs (construct_variants/results/ebov_aa/): aa_<gene>.png/.pdf (lollipop of enrichment along the protein + sample heatmap + table of group means), ebov_aa_overview.png/.pdf (untreated vs IFN abundance, enrichment vs p), ebov_aa_mutations.tsv.
Options (environment): MIN_DEPTH=30 MIN_ALT=5 VAF_MIN=0.02 STRAND_FILTER=1 EPS=0.002 SHOW_SILENT=0 MAX_ROWS=30 N_LABEL=12 VIRUS=EBOV PROJ=...      Usage: python 117_ebov_aa_mutations_ifn.py     (numpy, pandas, matplotlib)
"""
import os, sys, glob, itertools, gzip
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.patches import Patch

PROJ = os.environ.get("PROJ", os.path.join(os.path.expanduser("~"), "Documents/aviti_rnaseq/2025/aviti/20251009-100725_Maury_25161_WivfFacICOTHFAUJwSPeHBEFfCfNtCRpUSgVGUkz/results/Maury_25161"))
W = os.path.join(PROJ, "construct_variants"); OUT = os.path.join(W, "results", "ebov_aa"); os.makedirs(OUT, exist_ok=True)
MIN_DEPTH = int(os.environ.get("MIN_DEPTH", "30")); MIN_ALT = int(os.environ.get("MIN_ALT", "5")); VAF_MIN = float(os.environ.get("VAF_MIN", "0.02")); STRAND_FILTER = os.environ.get("STRAND_FILTER", "1") == "1"
EPS = float(os.environ.get("EPS", "0.002")); SHOW_SILENT = os.environ.get("SHOW_SILENT", "0") == "1"; MAX_ROWS = int(os.environ.get("MAX_ROWS", "30")); N_LABEL = int(os.environ.get("N_LABEL", "12")); VIRUS = os.environ.get("VIRUS", "EBOV")
INK = "#1F2933"; BASES = "ACGT"
GC = {"Untreated": "#6B7280", "IFNB": "#2B6CB0", "IFNG": "#D97706"}; EFF = {"missense": "#C0392B", "nonsense": "#111111", "silent": "#9AA5B1"}
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
                i += 1
    return C, F

refs = read_fasta(os.path.join(PROJ, "ref", "recombinant_model.fa")); is_cds = {k: (len(s) % 3 == 0 and s.startswith("ATG")) for k, s in refs.items()}
si = pd.read_csv(os.path.join(PROJ, "sample_info.csv"), sep=None, engine="python", dtype=str)
si["ifn"] = si["condition1"].fillna("Mock").replace({"None": "Mock", "none": "Mock"}); si["virus"] = si["condition2"].fillna("Mock").replace({"None": "Mock", "none": "Mock"})
sel = si[si["virus"] == VIRUS]; gmap = {r.sample_id: ("IFNB" if r.ifn == "IFNB" else "IFNG" if r.ifn in ("IFNY", "IFNG") else "Untreated") for r in sel.itertuples()}
files = {os.path.basename(p).split(".pileup")[0]: p for p in glob.glob(os.path.join(W, "pileup", "*.pileup.txt*"))}
samples = [s for s in sorted(gmap, key=lambda x: (["Untreated", "IFNB", "IFNG"].index(gmap[x]), x)) if s in files]
if not samples: sys.exit("No %s pileups found in %s" % (VIRUS, os.path.join(W, "pileup")))
print("%s samples: %s" % (VIRUS, ", ".join("%s(%s)" % (s, gmap[s]) for s in samples)))
DATA = {}
for s in samples: print("  reading", s, flush=True); DATA[s] = parse_pileup(files[s], refs)

# ---- candidates and abundance ---------------------------------------------------------------------------------------------
rows = []
for f, seq in refs.items():
    if not is_cds[f]: continue
    L = len(seq); depth = {s: DATA[s][0][f][:, :5].sum(axis=1) for s in samples}; cand = set()
    for s in samples:
        c, fw = DATA[s][0][f], DATA[s][1][f]; dp = depth[s]
        for a in range(4):
            isref = np.array([BASES[a] == b for b in seq]); alt = np.where(isref, 0, c[:, a]); rev = alt - np.where(isref, 0, fw[:, a])
            m = (dp >= MIN_DEPTH) & (alt >= MIN_ALT) & (alt / np.maximum(dp, 1) >= VAF_MIN) & ((np.where(isref, 0, fw[:, a]) >= 1) & (rev >= 1) if STRAND_FILTER else True)
            for pos in np.where(m)[0]: cand.add((pos, a))
    for pos, a in sorted(cand):
        cod = pos // 3; i = pos % 3; c0 = seq[3 * cod:3 * cod + 3]; c1 = c0[:i] + BASES[a] + c0[i + 1:]; a0, a1 = CODON.get(c0, "?"), CODON.get(c1, "?")
        eff = "silent" if a0 == a1 else "nonsense" if a1 == "*" else "missense"
        if eff == "silent" and not SHOW_SILENT: continue
        vaf = {s: (DATA[s][0][f][pos, a] / depth[s][pos] if depth[s][pos] >= MIN_DEPTH else np.nan) for s in samples}
        called = {s: (depth[s][pos] >= MIN_DEPTH and DATA[s][0][f][pos, a] >= MIN_ALT and vaf[s] >= VAF_MIN) for s in samples}
        rows.append(dict(feature=f, nt_pos=pos + 1, residue=cod + 1, ref_aa=a0, alt_aa=a1, label="%s%d%s" % (a0, cod + 1, a1), codon_change="%s>%s" % (c0, c1), nt_change="%s%d%s" % (seq[pos], pos + 1, BASES[a]), effect=eff, vaf=vaf, called=called))
if not rows: sys.exit("No amino-acid change passes the candidate filters (MIN_DEPTH=%d MIN_ALT=%d VAF_MIN=%g) in the %s samples - lower them (e.g. VAF_MIN=0.01)" % (MIN_DEPTH, MIN_ALT, VAF_MIN, VIRUS))
def perm_p(x, y):
    """exact two-sided permutation p for difference of means (x = untreated, y = IFN)"""
    v = np.array(list(x) + list(y)); n1 = len(x); obs = abs(np.mean(y) - np.mean(x)); cnt = tot = 0
    for idx in itertools.combinations(range(len(v)), n1):
        m = np.ones(len(v), bool); m[list(idx)] = False; d = abs(v[m].mean() - v[list(idx)].mean()); tot += 1; cnt += d >= obs - 1e-12
    return cnt / tot
recs = []
for r in rows:
    g = {k: [r["vaf"][s] for s in samples if gmap[s] == k and np.isfinite(r["vaf"][s])] for k in GC}
    mu = {k: (np.mean(v) if v else np.nan) for k, v in g.items()}; ifn_v = g["IFNB"] + g["IFNG"]; un = g["Untreated"]
    enr = np.log2((np.mean(ifn_v) + EPS) / (np.mean(un) + EPS)) if (ifn_v and un) else np.nan; p = perm_p(un, ifn_v) if (ifn_v and un) else np.nan
    ncall = {k: sum(r["called"][s] for s in samples if gmap[s] == k) for k in GC}; ncov = {k: len(g[k]) for k in GC}
    recs.append(dict(feature=r["feature"], residue=r["residue"], label=r["label"], effect=r["effect"], nt_change=r["nt_change"], codon_change=r["codon_change"], mean_untreated=mu["Untreated"], mean_IFNB=mu["IFNB"], mean_IFNG=mu["IFNG"],
                     mean_IFN=np.mean(ifn_v) if ifn_v else np.nan, log2_enrichment=enr, perm_p=p, called_untreated="%d/%d" % (ncall["Untreated"], ncov["Untreated"]), called_IFNB="%d/%d" % (ncall["IFNB"], ncov["IFNB"]), called_IFNG="%d/%d" % (ncall["IFNG"], ncov["IFNG"]),
                     **{"vaf_" + s: r["vaf"][s] for s in samples}))
T = pd.DataFrame(recs)
pv = T.perm_p.values; ok = np.isfinite(pv); q = np.full(len(T), np.nan); o = np.argsort(pv[ok]); ranks = np.empty(ok.sum()); adj = pv[ok][o] * ok.sum() / (np.arange(ok.sum()) + 1); adj = np.minimum.accumulate(adj[::-1])[::-1]; ranks[o] = np.minimum(adj, 1); q[ok] = ranks; T["FDR"] = q
T.sort_values(["feature", "residue"]).to_csv(os.path.join(OUT, "ebov_aa_mutations.tsv"), sep="\t", index=False, float_format="%.5g")
print("\n%d amino-acid changes (%s) in %d genes; enriched with IFN (log2 > 1): %d, depleted (< -1): %d" % (len(T), "incl. silent" if SHOW_SILENT else "missense + nonsense", T.feature.nunique(), (T.log2_enrichment > 1).sum(), (T.log2_enrichment < -1).sum()))
print(T.groupby("feature").agg(changes=("label", "size"), mean_vaf_untreated=("mean_untreated", "mean"), mean_vaf_IFN=("mean_IFN", "mean")).round(4).to_string())

def style(ax): ax.tick_params(colors=INK, labelsize=8); [sp.set_color("#9AA5B1") for sp in ax.spines.values()]
def save(fig, name): fig.savefig(os.path.join(OUT, name + ".png"), dpi=300, bbox_inches="tight"); fig.savefig(os.path.join(OUT, name + ".pdf"), bbox_inches="tight"); plt.close(fig)
cm = plt.get_cmap("YlOrRd").copy(); cm.set_bad("#E5E7EB")
edges = [j - 0.5 for j in range(1, len(samples)) if gmap[samples[j]] != gmap[samples[j - 1]]]
for f in refs:
    d = T[T.feature == f].sort_values("residue")
    if not len(d): continue
    L = len(refs[f]) // 3; pages = [d.iloc[i:i + MAX_ROWS] for i in range(0, len(d), MAX_ROWS)]
    for pi, dp in enumerate(pages):
        n = len(dp); fig = plt.figure(figsize=(max(11, 0.38 * len(samples) + 9), 3.6 + 0.26 * n + 1.3)); gs = fig.add_gridspec(2, 3, height_ratios=[2.3, 0.26 * n + 0.9], width_ratios=[30, 17, 1], hspace=0.7, wspace=0.04)
        a1 = fig.add_subplot(gs[0, 0]); a2 = fig.add_subplot(gs[1, 0]); a3 = fig.add_subplot(gs[1, 1]); cax = fig.add_subplot(gs[1, 2])
        a1.axhline(0, color="#9AA5B1", lw=3, zorder=1); dd = dp[np.isfinite(dp.log2_enrichment)]
        for r in dd.itertuples():
            a1.vlines(r.residue, 0, np.clip(r.log2_enrichment, -6, 6), color=EFF[r.effect], lw=1.0, zorder=2)
            a1.scatter(r.residue, np.clip(r.log2_enrichment, -6, 6), s=14 + 900 * np.nanmax([r.mean_untreated, r.mean_IFN]), color=EFF[r.effect], edgecolor="white", linewidth=0.5, zorder=3)
        for k, r in enumerate(dd.sort_values("residue").itertuples()): a1.annotate(r.label, (r.residue, np.clip(r.log2_enrichment, -6, 6)), xytext=(0, 6 + 7 * (k % 3) if r.log2_enrichment >= 0 else -(12 + 7 * (k % 3))), textcoords="offset points", ha="center", fontsize=6.5, color=INK)
        a1.set_xlim(0, L + 1); a1.set_ylim(-7.4, 7.8); a1.set_ylabel("log2 enrichment\n(IFN vs untreated)", fontsize=8, color=INK); a1.set_xlabel("residue number in %s (%d residues)   [point size = highest group mean frequency]" % (f, L), fontsize=8, color=INK); style(a1)
        a1.legend(handles=[Patch(color=c, label=e) for e, c in EFF.items() if (dp.effect == e).any()], fontsize=7, frameon=False, loc="upper right", ncol=3)
        a1.set_title("%s (%s samples): amino-acid changes, enrichment with interferon%s" % (f, VIRUS, "  (page %d of %d)" % (pi + 1, len(pages)) if len(pages) > 1 else ""), fontsize=10, color=INK, fontweight="bold", loc="left")
        M = np.array([[r["vaf_" + s] for s in samples] for _, r in dp.iterrows()], dtype=float)
        im = a2.imshow(np.ma.masked_invalid(np.maximum(M, 1e-4)), aspect="auto", cmap=cm, norm=LogNorm(vmin=0.002, vmax=1.0), interpolation="nearest")
        for i in range(n):
            for j in range(len(samples)):
                if np.isfinite(M[i, j]): a2.text(j, i, ("%.1f" % (100 * M[i, j])) if M[i, j] < 0.1 else "%d" % round(100 * M[i, j]), ha="center", va="center", fontsize=5, color="white" if M[i, j] > 0.2 else INK)
        a2.set_yticks(range(n)); a2.set_yticklabels(["%s (%s)" % (r.label, r.nt_change) for r in dp.itertuples()], fontsize=7, color=INK)
        a2.set_xticks(range(len(samples))); a2.set_xticklabels(["%s %s" % ({"Untreated": "EBOV", "IFNB": "IFNB+EBOV", "IFNG": "IFNG+EBOV"}[gmap[s]], s) for s in samples], rotation=75, ha="right", fontsize=6.5)
        for tl, s in zip(a2.get_xticklabels(), samples): tl.set_color(GC[gmap[s]]); tl.set_fontweight("bold")
        for e in edges: a2.axvline(e, color=INK, lw=1.2)
        a2.set_xticks(np.arange(-0.5, len(samples), 1), minor=True); a2.set_yticks(np.arange(-0.5, n, 1), minor=True); a2.grid(which="minor", color="white", lw=0.5); a2.tick_params(which="minor", length=0); style(a2)
        a2.set_title("frequency (%%) of the amino-acid-changing base in each sample; grey = < %d reads" % MIN_DEPTH, fontsize=7.5, color="grey", loc="left")
        cb = fig.colorbar(im, cax=cax); cb.set_label("frequency", fontsize=7); cb.ax.tick_params(labelsize=6)
        a3.set_xlim(0, 5); a3.set_ylim(n - 0.5, -0.5); a3.axis("off")
        for x, h in zip([0.35, 1.2, 2.05, 3.0, 4.0], ["untreated\nmean %", "IFNB\nmean %", "IFNG\nmean %", "log2 enr.\n(p)", "called in\nU / B / G"]): a3.text(x, -1.1, h, ha="center", va="bottom", fontsize=6.5, color=INK, fontweight="bold")
        for i, r in enumerate(dp.itertuples()):
            for x, v in zip([0.35, 1.2, 2.05], [r.mean_untreated, r.mean_IFNB, r.mean_IFNG]): a3.text(x, i, "n/a" if not np.isfinite(v) else "%.2f" % (100 * v), ha="center", va="center", fontsize=6.5, color=INK)
            a3.text(3.0, i, "n/a" if not np.isfinite(r.log2_enrichment) else "%+.1f (%.2g)" % (r.log2_enrichment, r.perm_p), ha="center", va="center", fontsize=6.5, color=INK, fontweight="bold" if (np.isfinite(r.perm_p) and r.perm_p < 0.05) else "normal")
            a3.text(4.0, i, "%s  %s  %s" % (r.called_untreated, r.called_IFNB, r.called_IFNG), ha="center", va="center", fontsize=6, color="grey")
        save(fig, "aa_%s%s" % (f, "" if len(pages) == 1 else "_p%d" % (pi + 1)))

# ---- overview ----------------------------------------------------------------------------------------------------------------
fig, axs = plt.subplots(1, 2, figsize=(12, 5.2)); fl = 0.0005
x = np.maximum(T.mean_untreated.fillna(0), fl); y = np.maximum(T.mean_IFN.fillna(0), fl)
for e, c in EFF.items():
    m = T.effect == e
    if m.any(): axs[0].scatter(100 * x[m], 100 * y[m], s=26, color=c, edgecolor="white", linewidth=0.4, label=e, alpha=0.9)
lim = [100 * fl * 0.8, 120]; axs[0].plot(lim, lim, color="#9AA5B1", lw=0.8, ls="--"); axs[0].set_xscale("log"); axs[0].set_yscale("log"); axs[0].set_xlim(lim); axs[0].set_ylim(lim)
axs[0].set_xlabel("mean frequency, untreated %s (%%)" % VIRUS, fontsize=8, color=INK); axs[0].set_ylabel("mean frequency, IFNB + IFNG %s (%%)" % VIRUS, fontsize=8, color=INK); axs[0].legend(fontsize=7, frameon=False); style(axs[0])
lab = T[np.isfinite(T.log2_enrichment)].assign(score=lambda d: d.log2_enrichment.abs() * np.maximum(d.mean_IFN.fillna(0), d.mean_untreated.fillna(0)) ** 0.5).sort_values("score", ascending=False).head(N_LABEL)
for r in lab.itertuples(): axs[0].annotate("%s %s" % (r.feature.replace("model_", ""), r.label), (100 * max(r.mean_untreated if np.isfinite(r.mean_untreated) else 0, fl), 100 * max(r.mean_IFN if np.isfinite(r.mean_IFN) else 0, fl)), xytext=(4, 3), textcoords="offset points", fontsize=6, color=INK)
axs[0].set_title("above the dashed line = more abundant with IFN", fontsize=9, color=INK, loc="left")
v = T[np.isfinite(T.perm_p)]
for e, c in EFF.items():
    m = v.effect == e
    if m.any(): axs[1].scatter(v.log2_enrichment[m], -np.log10(v.perm_p[m]), s=26, color=c, edgecolor="white", linewidth=0.4, alpha=0.9)
axs[1].axvline(0, color="#9AA5B1", lw=0.8); axs[1].axhline(-np.log10(0.05), color="#9AA5B1", lw=0.8, ls=":"); axs[1].set_xlabel("log2 enrichment (IFN vs untreated)", fontsize=8, color=INK); axs[1].set_ylabel("-log10 permutation p", fontsize=8, color=INK); style(axs[1])
for r in v.assign(sc=lambda d: -np.log10(d.perm_p) * d.log2_enrichment.abs()).sort_values("sc", ascending=False).head(N_LABEL).itertuples(): axs[1].annotate("%s %s" % (r.feature.replace("model_", ""), r.label), (r.log2_enrichment, -np.log10(r.perm_p)), xytext=(4, 3), textcoords="offset points", fontsize=6, color=INK)
axs[1].set_title("right of 0 = enriched with IFN; dotted line p = 0.05 (few samples: use as a ranking)", fontsize=9, color=INK, loc="left")
fig.suptitle("%s construct: amino-acid changes, abundance and interferon enrichment (%d changes, %d genes)" % (VIRUS, len(T), T.feature.nunique()), fontsize=10, color=INK, fontweight="bold"); fig.tight_layout(); save(fig, "ebov_aa_overview")
print("\nWrote figures + ebov_aa_mutations.tsv to", OUT)
