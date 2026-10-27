#!/usr/bin/env python3
"""Look up ONE position (or one codon) of ONE construct gene in every sample, with no thresholds: depth, reference count and the count / frequency of every other base, deletion and insertion (from the pileups of 112).
Use it when a mutation you care about does not show up in 113/114/115 (those only list calls above the error floor) - this shows what is really in the reads, sample by sample, and draws it.
Usage:  python 116_variant_lookup.py GENE RESIDUE            e.g.  python 116_variant_lookup.py model_GP 104        (codon 104 = nt 310-312, all three positions are shown; needs a CDS record)
        python 116_variant_lookup.py GENE nt:POS             e.g.  python 116_variant_lookup.py model_GP nt:310     (one nucleotide)
Writes construct_variants/results/lookup_<GENE>_<where>.tsv and .png/.pdf (stacked bars of the non-reference fraction per sample, grouped by condition; numbers = depth).   Needs numpy, pandas, matplotlib.
"""
import os, sys, glob, re
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

if len(sys.argv) != 3: sys.exit(__doc__)
GENE, WHERE = sys.argv[1], sys.argv[2]
PROJ = os.environ.get("PROJ", os.path.join(os.path.expanduser("~"), "Documents/aviti_rnaseq/2025/aviti/20251009-100725_Maury_25161_WivfFacICOTHFAUJwSPeHBEFfCfNtCRpUSgVGUkz/results/Maury_25161"))
W = os.path.join(PROJ, "construct_variants"); RES = os.path.join(W, "results"); os.makedirs(RES, exist_ok=True); INK = "#1F2933"
ORDER = ["Mock_Mock", "IFNB_Mock", "IFNY_Mock", "Mock_EBOV", "IFNB_EBOV", "IFNY_EBOV", "Mock_MARV", "IFNB_MARV", "IFNY_MARV", "Mock_EVLP", "Mock_MVLP"]
LAB = {"Mock_Mock": "Mock", "IFNB_Mock": "IFNB", "IFNY_Mock": "IFNG", "Mock_EBOV": "EBOV", "IFNB_EBOV": "IFNB+EBOV", "IFNY_EBOV": "IFNG+EBOV", "Mock_MARV": "MARV", "IFNB_MARV": "IFNB+MARV", "IFNY_MARV": "IFNG+MARV", "Mock_EVLP": "EVLP", "Mock_MVLP": "MVLP"}
CODON = {}
_aa = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
for i, a in enumerate("TCAG"):
    for j, b in enumerate("TCAG"):
        for k, c in enumerate("TCAG"): CODON[a + b + c] = _aa[16 * i + 4 * j + k]
seq, n = {}, None
for line in open(os.path.join(PROJ, "ref", "recombinant_model.fa")):
    line = line.strip()
    if line.startswith(">"): n = line[1:].split()[0]; seq[n] = []
    elif n: seq[n].append(line.upper())
seq = {k: "".join(v) for k, v in seq.items()}
if GENE not in seq: sys.exit("%s not in the construct FASTA. Records: %s" % (GENE, ", ".join(seq)))
s = seq[GENE]
if WHERE.startswith("nt:"): pos = [int(WHERE[3:])]; tag = "nt%d" % pos[0]; title = "%s nt %d (%s)" % (GENE, pos[0], s[pos[0] - 1])
else:
    r = int(WHERE); pos = [3 * (r - 1) + 1, 3 * (r - 1) + 2, 3 * (r - 1) + 3]; cod = s[pos[0] - 1:pos[2]]; tag = "res%d" % r
    title = "%s residue %d: codon %s = %s (nt %d-%d)" % (GENE, r, cod, CODON.get(cod, "?"), pos[0], pos[2])
if max(pos) > len(s): sys.exit("position beyond the %d nt of %s" % (len(s), GENE))
si = pd.read_csv(os.path.join(PROJ, "sample_info.csv"), sep=None, engine="python", dtype=str)
g = si["condition1"].fillna("Mock").replace({"None": "Mock", "none": "Mock"}) + "_" + si["condition2"].fillna("Mock").replace({"None": "Mock", "none": "Mock"}); grp = dict(zip(si["sample_id"], g))
files = {os.path.basename(p).split(".pileup")[0]: p for p in glob.glob(os.path.join(W, "pileup", "*.pileup.txt*"))}
if not files: sys.exit("No pileups in %s" % os.path.join(W, "pileup"))
samples = sorted(files, key=lambda x: (ORDER.index(grp[x]) if grp.get(x) in ORDER else 99, x)); BASES = "ACGT"; COLS = ["A", "C", "G", "T", "del", "ins"]
def counts(path):
    out = {p: np.zeros(6, dtype=int) for p in pos}; refc = {p: "N" for p in pos}
    for line in open(path):
        q = line.split("\t")
        if len(q) < 5 or q[0] != GENE: continue
        p = int(q[1])
        if p not in out: continue
        refc[p] = q[2].upper()
        if int(q[3]) == 0: continue          # no reads (samtools prints "*" as the base string)
        refc[p] = q[2].upper(); bs = q[4]; ri = BASES.find(refc[p]); i, m = 0, len(bs)
        while i < m:
            ch = bs[i]
            if ch == "^": i += 2; continue
            if ch == "$": i += 1; continue
            if ch in ".,":
                if ri >= 0: out[p][ri] += 1
                i += 1; continue
            if ch in "ACGTacgt": out[p][BASES.index(ch.upper())] += 1; i += 1; continue
            if ch in "*#": out[p][4] += 1; i += 1; continue
            if ch in "+-":
                j = i + 1
                while j < m and bs[j].isdigit(): j += 1
                if ch == "+": out[p][5] += 1
                i = j + int(bs[i + 1:j]); continue
            i += 1
    return out, refc
rows = []
for sm in samples:
    c, rc = counts(files[sm])
    for p in pos:
        dp = int(c[p][:5].sum()); r = {"sample": sm, "condition": LAB.get(grp.get(sm, "?"), grp.get(sm, "?")), "pos": p, "ref": rc[p], "depth": dp}
        for k, nme in enumerate(COLS): r["n_" + nme] = int(c[p][k]); r["f_" + nme] = (c[p][k] / dp if dp else np.nan)
        r["nonref_fraction"] = (1 - c[p][BASES.find(rc[p])] / dp) if (dp and rc[p] in BASES) else np.nan; rows.append(r)
T = pd.DataFrame(rows); out = os.path.join(RES, "lookup_%s_%s.tsv" % (GENE, tag)); T.to_csv(out, sep="\t", index=False, float_format="%.5f")
pd.set_option("display.width", 220); show = T[["sample", "condition", "pos", "ref", "depth", "n_A", "n_C", "n_G", "n_T", "n_del", "n_ins", "nonref_fraction"]].copy(); show["nonref_fraction"] = (100 * show["nonref_fraction"]).round(2)
print(title); print(show.to_string(index=False)); print("\nwrote", out)
cols = {"A": "#2E7D32", "C": "#1565C0", "G": "#F9A825", "T": "#C62828", "del": "#424242", "ins": "#8E24AA"}
fig, axs = plt.subplots(len(pos), 1, figsize=(max(7, 0.38 * len(samples) + 3), 2.6 * len(pos) + 1.4), sharex=True, squeeze=False)
for ax, p in zip(axs[:, 0], pos):
    d = T[T.pos == p].set_index("sample").reindex(samples); bot = np.zeros(len(samples)); ref = d["ref"].iloc[0]
    for nme in COLS:
        if nme == ref: continue
        v = (100 * d["f_" + nme].fillna(0)).values; ax.bar(range(len(samples)), v, bottom=bot, color=cols[nme], label=nme, width=0.8); bot += v
    for i, sm in enumerate(samples): ax.text(i, bot[i] + 0.02 * max(bot.max(), 1), "%d" % d["depth"].iloc[i], ha="center", va="bottom", fontsize=5, color="grey", rotation=90)
    ax.set_ylabel("non-reference %% of reads\nnt %d (ref %s)" % (p, ref), fontsize=8, color=INK); ax.set_ylim(0, max(bot.max() * 1.35, 1)); ax.legend(fontsize=6, ncol=6, frameon=False, loc="upper left"); ax.tick_params(colors=INK, labelsize=7)
axs[-1, 0].set_xticks(range(len(samples))); axs[-1, 0].set_xticklabels(["%s %s" % (LAB.get(grp.get(sm, "?"), grp.get(sm, "?")), sm) for sm in samples], rotation=75, ha="right", fontsize=7, color=INK)
fig.suptitle(title + "\n(grey numbers = read depth; no thresholds applied)", fontsize=10, color=INK, fontweight="bold"); fig.tight_layout()
for e in ("png", "pdf"): fig.savefig(os.path.join(RES, "lookup_%s_%s.%s" % (GENE, tag, e)), dpi=300, bbox_inches="tight")
print("wrote figure lookup_%s_%s.png/.pdf in %s" % (GENE, tag, RES))
