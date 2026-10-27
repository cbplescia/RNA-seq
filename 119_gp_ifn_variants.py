#!/usr/bin/env python3
"""Preliminary-data figure: which substitutions in the GP gene of the EBOV construct are associated with IFNG treatment?   (IFNG+EBOV vs untreated EBOV; IFNB, MARV, VLP and mock samples are left out.)
Reads the pileups from 112_construct_variants_align.sh. Nothing is pre-selected: EVERY possible single-base substitution along the gene is measured in every replicate, and the figure shows how the IFNG-associated ones stand out from all the others.
  panel A  Manhattan-style plot: x = residue, y = change in frequency (IFNG+EBOV mean - untreated EBOV mean, percentage points) for every substitution with coverage in both groups; grey = all, coloured = candidates (labelled with the amino-acid change)
  panel B  the candidates replicate by replicate (untreated vs IFNG+EBOV), so you can see the frequencies are reproducible across replicates and not driven by one sample
  candidate = every covered IFNG replicate is higher than every covered untreated replicate (complete separation), the mean difference is >= MIN_DELTA points, the IFNG mean is >= MIN_IFNG %, and >= MIN_ALT reads support it in >= 2 IFNG replicates;
              needs >= 2 covered replicates per group (depth >= MIN_DEPTH). Exact one-sided permutation p is reported; with 3 vs 3 replicates the smallest possible p is 0.05, so this is hypothesis-generating, not proof.
WHAT THIS CAN AND CANNOT SHOW: it shows variants that are more frequent in IFNG-treated EBOV samples than in untreated ones. It cannot by itself show that IFNG CAUSES them (RNA-seq cannot separate a new mutation from RT/sequencing error, RNA editing, selection, or a batch difference);
a variant that is the same in every sample (e.g. residue 104 at ~3-5 % in all EBOV samples) is NOT reported - it is background. If no variant passes, the figure says so.
Outputs (construct_variants/results/gp_ifng/): gp_ifng_variants.png/.pdf, gp_ifng_candidates.tsv (all candidates with replicate frequencies), gp_ifng_all_substitutions.tsv (every substitution tested).
Options (environment): IFN=IFNG (default) or IFNB (IFNB+EBOV vs untreated EBOV; outputs in gp_ifnb/)  GENE=model_GP  MIN_DEPTH=30 MIN_ALT=5 MIN_DELTA=1.0 MIN_IFNG=2.0 N_PANELS=8 N_LABEL=15  PROJ=...      Usage: python 119_gp_ifn_variants.py     (numpy, pandas, matplotlib)
"""
import os, sys, glob, gzip, itertools
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

PROJ = os.environ.get("PROJ", os.path.join(os.path.expanduser("~"), "Documents/aviti_rnaseq/2025/aviti/20251009-100725_Maury_25161_WivfFacICOTHFAUJwSPeHBEFfCfNtCRpUSgVGUkz/results/Maury_25161"))
W = os.path.join(PROJ, "construct_variants"); TREAT = os.environ.get("IFN", "IFNG").upper().replace("IFNY", "IFNG"); assert TREAT in ("IFNG", "IFNB"), "IFN must be IFNG or IFNB"; TAG = TREAT.lower()
OUT = os.path.join(W, "results", "gp_" + TAG); os.makedirs(OUT, exist_ok=True)
GENE = os.environ.get("GENE", "model_GP"); MIN_DEPTH = int(os.environ.get("MIN_DEPTH", "30")); MIN_ALT = int(os.environ.get("MIN_ALT", "5")); MIN_DELTA = float(os.environ.get("MIN_DELTA", "1.0")); MIN_IFNG = float(os.environ.get("MIN_IFNG", "2.0"))
N_PANELS = int(os.environ.get("N_PANELS", "8")); N_LABEL = int(os.environ.get("N_LABEL", "15")); INK = "#1F2933"; BASES = "ACGT"
C_UN = "#6B7280"; C_G = "#2B6CB0" if os.environ.get("IFN", "IFNG").upper() == "IFNB" else "#D97706"; EFF = {"missense": "#C0392B", "nonsense": "#111111", "silent": "#2B8A6E"}
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
refs = read_fasta(os.path.join(PROJ, "ref", "recombinant_model.fa"))
if GENE not in refs: sys.exit("%s not in the construct FASTA. Records: %s" % (GENE, ", ".join(refs)))
seq = refs[GENE]; L = len(seq); is_cds = (L % 3 == 0 and seq.startswith("ATG"))

def parse(path):
    C = np.zeros((L, 4), dtype=np.int64); D = np.zeros(L, dtype=np.int64); opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as fh:
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) < 5 or p[0] != GENE or int(p[3]) == 0: continue
            pos = int(p[1]) - 1; ref = p[2].upper(); bs = p[4]
            if ref not in BASES: continue
            ri = BASES.index(ref); i, n = 0, len(bs); tot = 0
            while i < n:
                ch = bs[i]
                if ch == "^": i += 2; continue
                if ch == "$": i += 1; continue
                if ch in ".,": C[pos, ri] += 1; tot += 1; i += 1; continue
                if ch in "ACGTacgt": C[pos, BASES.index(ch.upper())] += 1; tot += 1; i += 1; continue
                if ch in "*#": tot += 1; i += 1; continue
                if ch in "+-":
                    j = i + 1
                    while j < n and bs[j].isdigit(): j += 1
                    i = j + int(bs[i + 1:j]); continue
                i += 1
            D[pos] = tot
    return C, D
si = pd.read_csv(os.path.join(PROJ, "sample_info.csv"), sep=None, engine="python", dtype=str)
si["ifn"] = si["condition1"].fillna("Mock").replace({"None": "Mock", "none": "Mock"}); si["virus"] = si["condition2"].fillna("Mock").replace({"None": "Mock", "none": "Mock"})
files = {os.path.basename(p).split(".pileup")[0]: p for p in glob.glob(os.path.join(W, "pileup", "*.pileup.txt*"))}
UN = [r.sample_id for r in si.itertuples() if r.virus == "EBOV" and r.ifn == "Mock" and r.sample_id in files]; GG = [r.sample_id for r in si.itertuples() if r.virus == "EBOV" and r.ifn in (("IFNB",) if TREAT == "IFNB" else ("IFNY", "IFNG")) and r.sample_id in files]
if not UN or not GG: sys.exit("Need untreated EBOV and %s+EBOV pileups (found %d and %d)" % (TREAT, len(UN), len(GG)))
print("untreated EBOV: %s | %s+EBOV: %s | gene %s (%d nt)" % (", ".join(UN), TREAT, ", ".join(GG), GENE, L))
DATA = {s: parse(files[s]) for s in UN + GG}

def perm_p_onesided(x, y):
    v = np.array(list(x) + list(y)); n1 = len(x); obs = np.mean(y) - np.mean(x); cnt = tot = 0
    for idx in itertools.combinations(range(len(v)), n1):
        m = np.ones(len(v), bool); m[list(idx)] = False; d = v[m].mean() - v[list(idx)].mean(); tot += 1; cnt += d >= obs - 1e-12
    return cnt / tot
rows = []
for pos in range(L):
    for a in range(4):
        if BASES[a] == seq[pos]: continue
        fu = [100 * DATA[s][0][pos, a] / DATA[s][1][pos] if DATA[s][1][pos] >= MIN_DEPTH else np.nan for s in UN]; fg = [100 * DATA[s][0][pos, a] / DATA[s][1][pos] if DATA[s][1][pos] >= MIN_DEPTH else np.nan for s in GG]
        cu = [v for v in fu if np.isfinite(v)]; cg = [v for v in fg if np.isfinite(v)]
        if len(cu) < 2 or len(cg) < 2: continue
        nalt = sum(DATA[s][0][pos, a] >= MIN_ALT for s in GG if DATA[s][1][pos] >= MIN_DEPTH)
        eff = ""; lab = "c.%d%s>%s" % (pos + 1, seq[pos], BASES[a]); res = (pos // 3 + 1) if is_cds else pos + 1
        if is_cds:
            cod = pos // 3; i = pos % 3; c0 = seq[3 * cod:3 * cod + 3]; c1 = c0[:i] + BASES[a] + c0[i + 1:]; a0, a1 = CODON.get(c0, "?"), CODON.get(c1, "?"); eff = "silent" if a0 == a1 else "nonsense" if a1 == "*" else "missense"; lab = "%s%d%s" % (a0, cod + 1, a1)
        delta = np.mean(cg) - np.mean(cu); cand = (min(cg) > max(cu)) and delta >= MIN_DELTA and np.mean(cg) >= MIN_IFNG and nalt >= 2
        rows.append(dict(pos=pos + 1, residue=res, ref=seq[pos], alt=BASES[a], label=lab, effect=eff, mean_untreated=np.mean(cu), mean_treated=np.mean(cg), delta=delta, candidate=cand, p_onesided=perm_p_onesided(cu, cg) if cand else np.nan,
                         **{"untreated_" + s: v for s, v in zip(UN, fu)}, **{TREAT + "_" + s: v for s, v in zip(GG, fg)}))
T = pd.DataFrame(rows); C = T[T.candidate].sort_values("delta", ascending=False)
T.to_csv(os.path.join(OUT, "gp_%s_all_substitutions.tsv" % TAG), sep="\t", index=False, float_format="%.4g"); C.to_csv(os.path.join(OUT, "gp_%s_candidates.tsv" % TAG), sep="\t", index=False, float_format="%.4g")
print("%d substitutions tested; %d candidates (complete separation, delta >= %g points, %s mean >= %g %%)" % (len(T), len(C), MIN_DELTA, TREAT, MIN_IFNG))
if len(C): print(C[["label", "effect", "pos", "mean_untreated", "mean_treated", "delta", "p_onesided"]].round(3).to_string(index=False))

def style(ax): ax.tick_params(colors=INK, labelsize=8); [sp.set_color("#9AA5B1") for sp in ax.spines.values()]
nP = min(N_PANELS, len(C)); ncol = 4 if nP > 3 else max(nP, 1); nrow = max(1, int(np.ceil(nP / ncol)))
fig = plt.figure(figsize=(13, 4.6 + 2.5 * nrow)); gs = fig.add_gridspec(1 + nrow, ncol, height_ratios=[1.9] + [1.2] * nrow, hspace=0.75, wspace=0.45)
a1 = fig.add_subplot(gs[0, :]); nc = T[~T.candidate]
a1.scatter(nc.residue, nc.delta, s=5, color="#C9CDD3", linewidth=0, zorder=1, rasterized=True)
for e, c in EFF.items():
    m = C[C.effect == e] if is_cds else (C if e == "missense" else C.iloc[0:0])
    if len(m): a1.scatter(m.residue, m.delta, s=34, color=c, edgecolor="white", linewidth=0.5, zorder=3, label="candidate, %s" % e)
a1.axhline(0, color="#9AA5B1", lw=0.8); a1.axhline(MIN_DELTA, color="#9AA5B1", lw=0.8, ls=":")
for k, r in enumerate(C.head(N_LABEL).sort_values("residue").itertuples()): a1.annotate(r.label, (r.residue, r.delta), xytext=(0, 6 + 8 * (k % 3)), textcoords="offset points", ha="center", fontsize=7, color=INK, arrowprops=dict(arrowstyle="-", lw=0.3, color="#9AA5B1") if k % 3 else None)
a1.set_xlim(0, (L // 3 if is_cds else L) + 1); a1.set_xlabel("residue number in %s" % GENE if is_cds else "nucleotide position", fontsize=9, color=INK); a1.set_ylabel("%s+EBOV minus untreated EBOV\n(frequency, percentage points)" % TREAT, fontsize=9, color=INK); style(a1)
if len(C): a1.legend(fontsize=7.5, frameon=False, loc="upper right")
else: a1.text(0.5, 0.8, "no substitution passes the candidate criteria", transform=a1.transAxes, ha="center", fontsize=11, color=INK)
a1.set_title("A. %s: every substitution tested (%d, grey) - %s-associated candidates highlighted (%d)" % (GENE, len(T), TREAT, len(C)), fontsize=10.5, color=INK, fontweight="bold", loc="left")
for k in range(nP):
    r = C.iloc[k]; ax = fig.add_subplot(gs[1 + k // ncol, k % ncol])
    for gi, (names, col, xs) in enumerate([(UN, C_UN, 0), (GG, C_G, 1)]):
        vals = [r[("untreated_" if gi == 0 else TREAT + "_") + s] for s in names]; vv = [v for v in vals if np.isfinite(v)]
        ax.bar(xs, np.mean(vv), width=0.6, color=col, alpha=0.35, zorder=1); ax.scatter(np.full(len(vv), xs) + np.linspace(-0.12, 0.12, len(vv)), vv, color=col, s=26, edgecolor="white", linewidth=0.5, zorder=3)
    ax.set_xticks([0, 1]); ax.set_xticklabels(["EBOV\n(n=%d)" % sum(np.isfinite(r["untreated_" + s]) for s in UN), "%s+EBOV\n(n=%d)" % (TREAT, sum(np.isfinite(r[TREAT + "_" + s]) for s in GG))], fontsize=7.5); ax.set_xlim(-0.6, 1.6)
    ax.set_ylabel("frequency (%)", fontsize=8, color=INK); style(ax); ax.set_title("B%d. %s   +%.1f pts, p = %.2g" % (k + 1, r.label, r.delta, r.p_onesided), fontsize=8.5, color=INK, fontweight="bold")
fig.suptitle("Preliminary data: GP substitutions associated with %s treatment in EBOV-infected cells" % TREAT, fontsize=12, color=INK, fontweight="bold", y=0.995)
fig.text(0.5, 0.005, "Panel B (below A): top candidates replicate by replicate; each dot = one replicate, bar = mean. Frequencies are alt reads / depth; positions with < %d reads in a replicate are excluded. Variants identical in all samples are background and are not highlighted." % MIN_DEPTH, ha="center", fontsize=7, color="grey")
for e in ("png", "pdf"): fig.savefig(os.path.join(OUT, "gp_%s_variants." % TAG + e), dpi=300, bbox_inches="tight")
print("wrote", OUT)
