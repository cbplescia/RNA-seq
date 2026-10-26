#!/usr/bin/env python3
"""Figures for the output of 97_lysine_scan.py (needs matplotlib; run 97 first).

  <name>_lysine_map.png / .pdf              position map: every lysine along the sequence, the flagged ones coloured and labelled (SUMO consensus, histone-like, K/R cluster),
                                            plus the lysine density and the K+R (basic) density in a sliding window
  <name>_lysine_context_heatmap.png / .pdf  amino-acid enrichment around the lysines (positions -W..+W): what residues are over- / under-represented next to a lysine compared with the whole sequence
  <name>_lysine_flag_counts.png / .pdf      how many lysines carry each flag
Usage:  python 98_lysine_scan_plots.py eGP.FASTA     (reads eGP.FASTA and eGP_lysines.tsv from the current folder; the first record is plotted, or RECORD=2 for the second; DENSITY_WINDOW=21; LABEL_ALL=1 also labels the K/R-cluster sites)
"""
import os, sys, csv
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

AA = "ACDEFGHIKLMNPQRSTVWY"
COL = {"SUMO_consensus": "#C53030", "histone_like": "#2B6CB0", "K_cluster": "#E08214", "none": "#9AA0A6"}
DW = int(os.environ.get("DENSITY_WINDOW", "21"))


def read_fasta(path):
    recs, name, seq = [], None, []
    for line in open(path):
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if name is not None:
                recs.append((name, "".join(seq).upper().replace("*", "")))
            name, seq = line[1:].strip(), []
        else:
            seq.append("".join(c for c in line if c.isalpha() or c == "*"))
    if name is not None:
        recs.append((name, "".join(seq).upper().replace("*", "")))
    return recs


def sliding(seq, residues, w):
    h = w // 2
    out = []
    for i in range(len(seq)):
        seg = seq[max(0, i - h): i + h + 1]
        out.append(100.0 * sum(c in residues for c in seg) / len(seg))
    return out


def main():
    if len(sys.argv) < 2:
        sys.exit("Usage: python 98_lysine_scan_plots.py eGP.FASTA   (run 97_lysine_scan.py first)")
    fasta = sys.argv[1]
    base = os.path.splitext(os.path.basename(fasta))[0]
    tsv = f"{base}_lysines.tsv"
    if not os.path.exists(tsv):
        sys.exit(f"{tsv} not found - run:  python 97_lysine_scan.py {fasta}")
    recs = read_fasta(fasta)
    ridx = int(os.environ.get("RECORD", "1")) - 1
    name, seq = recs[ridx]
    short = name.split()[0]
    rows = [r for r in csv.DictReader(open(tsv), delimiter="\t") if r["record"] == short]
    if not rows:
        sys.exit(f"No rows for record '{short}' in {tsv}")
    L = len(seq)
    pos = [int(r["position"]) for r in rows]
    fl = {int(r["position"]): [x for x in r["flags"].split(",") if x] for r in rows}
    win = {int(r["position"]): r["window"] for r in rows}
    flagged = [p for p in pos if fl[p]]
    print(f"{short}: {L} aa, {len(pos)} lysines, {len(flagged)} with a flag")

    # ---- 1. position map ----
    fig, ax = plt.subplots(3, 1, figsize=(12, 6.5), sharex=True, gridspec_kw={"height_ratios": [2.2, 1.2, 1.2], "hspace": 0.12})
    a = ax[0]
    for p in pos:
        a.vlines(p, 0, 1, color=COL["none"], lw=0.8)
    for p in flagged:
        top = fl[p][0]
        a.vlines(p, 0, 1, color=COL[top], lw=1.8)
        a.plot(p, 1.0, "o", color=COL[top], ms=5)
    ymax = 1.0
    last = -999
    lvl = 0
    LABEL_ALL = os.environ.get("LABEL_ALL") == "1"            # by default only SUMO-consensus / histone-like sites are labelled; K/R-cluster sites are too many to label legibly
    for p in sorted(flagged):
        if not LABEL_ALL and fl[p] == ["K_cluster"]:
            continue
        lvl = (lvl + 1) % 3 if p - last < 18 else 0
        a.text(p, 1.08 + 0.12 * lvl, f"K{p}", rotation=90, ha="center", va="bottom", fontsize=6.5, color=COL[fl[p][0]])
        last = p
    a.set_ylim(0, 1.7)
    a.set_yticks([])
    a.set_title(f"Lysines in {short} ({L} aa; {len(pos)} K, {len(flagged)} flagged)", loc="left", fontsize=11, fontweight="bold")
    a.legend(handles=[Line2D([0], [0], color=COL["none"], lw=1.5, label="lysine (no flag)")] + [Line2D([0], [0], color=COL[k], lw=2, label=k) for k in ("SUMO_consensus", "histone_like", "K_cluster")],
             loc="upper right", fontsize=8, ncol=4, frameon=False)
    for s in ("top", "right", "left"):
        a.spines[s].set_visible(False)
    xs = list(range(1, L + 1))
    ax[1].fill_between(xs, sliding(seq, "K", DW), color="#2B6CB0", alpha=0.35, lw=0)
    ax[1].plot(xs, sliding(seq, "K", DW), color="#2B6CB0", lw=1)
    ax[1].set_ylabel(f"% K\n({DW}-aa window)", fontsize=8)
    ax[2].fill_between(xs, sliding(seq, "KR", DW), color="#E08214", alpha=0.35, lw=0)
    ax[2].plot(xs, sliding(seq, "KR", DW), color="#E08214", lw=1)
    ax[2].set_ylabel(f"% K+R\n({DW}-aa window)", fontsize=8)
    ax[2].set_xlabel("residue position")
    ax[2].set_xlim(1, L)
    for b in ax[1:]:
        for s in ("top", "right"):
            b.spines[s].set_visible(False)
    fig.savefig(f"{base}_lysine_map.png", dpi=300, bbox_inches="tight")
    fig.savefig(f"{base}_lysine_map.pdf", bbox_inches="tight")
    plt.close(fig)

    # ---- 2. context enrichment heatmap ----
    W = max(len(w) for w in win.values()) // 2
    bg = {c: seq.count(c) / L for c in AA}
    mat = []
    for c in AA:
        row = []
        for j in range(-W, W + 1):
            n = sum(1 for p in pos if win[p][W + j] == c)
            tot = sum(1 for p in pos if win[p][W + j] in AA)
            f = n / tot if tot else 0.0
            row.append(float("nan") if j == 0 else ((f + 0.5 / max(tot, 1)) / (bg[c] + 0.5 / L)))
        mat.append(row)
    import math
    lm = [[(math.log2(v) if v == v and v > 0 else float("nan")) for v in r] for r in mat]
    lim = max(0.5, min(3, max(abs(v) for r in lm for v in r if v == v)))
    fig, a = plt.subplots(figsize=(0.45 * (2 * W + 1) + 2, 6))
    im = a.imshow(lm, cmap="RdBu_r", vmin=-lim, vmax=lim, aspect="auto")
    a.set_xticks(range(2 * W + 1))
    a.set_xticklabels([f"{j:+d}" if j else "K" for j in range(-W, W + 1)], fontsize=8)
    a.set_yticks(range(len(AA)))
    a.set_yticklabels(list(AA), fontsize=9)
    a.set_xlabel("position relative to the lysine")
    a.set_title(f"Residues around the {len(pos)} lysines of {short}", loc="left", fontsize=11, fontweight="bold")
    cb = fig.colorbar(im, ax=a, fraction=0.03, pad=0.02)
    cb.set_label("log2 enrichment vs the sequence composition", fontsize=8)
    fig.text(0.01, -0.02, "Small numbers of lysines make this noisy: read it as a pattern check, not a statistic.", fontsize=7, color="grey")
    fig.savefig(f"{base}_lysine_context_heatmap.png", dpi=300, bbox_inches="tight")
    fig.savefig(f"{base}_lysine_context_heatmap.pdf", bbox_inches="tight")
    plt.close(fig)

    # ---- 3. flag counts ----
    counts = {k: sum(1 for p in pos if k in fl[p]) for k in ("SUMO_consensus", "histone_like", "K_cluster")}
    counts["no flag"] = sum(1 for p in pos if not fl[p])
    fig, a = plt.subplots(figsize=(5.5, 3.6))
    ks = list(counts)
    a.bar(ks, [counts[k] for k in ks], color=[COL.get(k, COL["none"]) for k in ks], width=0.65)
    for i, k in enumerate(ks):
        a.text(i, counts[k] + 0.1, str(counts[k]), ha="center", fontsize=9)
    a.set_ylabel("lysines")
    a.set_title(f"{short}: lysines by flag", loc="left", fontsize=11, fontweight="bold")
    a.tick_params(axis="x", labelsize=8)
    for s in ("top", "right"):
        a.spines[s].set_visible(False)
    fig.savefig(f"{base}_lysine_flag_counts.png", dpi=300, bbox_inches="tight")
    fig.savefig(f"{base}_lysine_flag_counts.pdf", bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {base}_lysine_map, {base}_lysine_context_heatmap and {base}_lysine_flag_counts (.png and .pdf)")


if __name__ == "__main__":
    main()
