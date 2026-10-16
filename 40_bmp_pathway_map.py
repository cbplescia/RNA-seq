#!/usr/bin/env python3
"""Pathway map of BMP (LBPA) biosynthesis and the acyl-chain supply behind it, with each enzyme's IFNB-vs-Mock
transcript log2FC (colour; stars = padj) and the matching pilot-lipidomics readouts.

Usage:
  python3 40_bmp_pathway_map.py --rna DE_v2_IFNB_vs_Mock_protein_coding.csv \
      --lipid 01_lipidomics_ranked_IFNB_vs_Untreated.csv --out bmp_pathway_map.png
Needs pandas and matplotlib (conda install matplotlib pandas).

The pathway is drawn from current understanding of BMP/LBPA synthesis (PA -> CDP-DAG -> PGP -> PG in ER/mitochondria,
PG delivered to the late endosome/lysosome and converted to BMP by PLA2G15 + CLN5 transacylation; PLD3/PLD4 are drawn
as a proposed alternative route). It is a simplified map, not a complete one.
"""
import argparse, re
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Rectangle
from matplotlib.colors import LinearSegmentedColormap, Normalize

ap = argparse.ArgumentParser()
ap.add_argument("--rna", required=True); ap.add_argument("--lipid", required=True)
ap.add_argument("--out", default="bmp_pathway_map.png")
args = ap.parse_args()

de = pd.read_csv(args.rna).set_index("gene_name")
cmap = LinearSegmentedColormap.from_list("div", ["#2C5F9E", "#F7F7F7", "#B2382C"])
norm = Normalize(-1.2, 1.2)
INK = "#1F2933"

def gene_info(g):
    if g not in de.index: return None
    r = de.loc[g]; p = r.padj
    stars = "***" if p < 1e-3 else "**" if p < 1e-2 else "*" if p < 0.05 else "ns"
    return r.log2FoldChange, stars

# ---------------- lipid readouts ----------------
L = pd.read_csv(args.lipid)
bad = L.Lipid_ID.str.contains(r"\+D\d|-D\d|\(s\)|COOH|OOH") | L.any_missing_Untreated | L.any_missing_IFNB
L = L[~bad].copy()
L["ch"] = L.Lipid_ID.apply(lambda s: re.findall(r"(?<![0-9])(\d{1,2}):(\d{1,2})", s))
bmp_all = L[L.Class.isin(["LBPA", "BMP"])]
bmp = bmp_all[bmp_all.ch.str.len() >= 2]
total_lfc = np.log2(bmp_all.B_Norm.sum() / bmp_all.Untreated_Norm.sum())
col = "log2FC_IFNB_vs_Untreated"
def feat(df, fn):
    m = df.ch.apply(fn); return df[col][m].mean(), df[col][~m].mean(), int(m.sum())
f181 = feat(bmp, lambda c: ("18", "1") in c)
fdha = feat(bmp, lambda c: ("22", "6") in c)
pct = lambda name: (100 * bmp_all.loc[bmp_all.Lipid_ID == name, "Untreated_Norm"].sum() / bmp_all.Untreated_Norm.sum(),
                    100 * bmp_all.loc[bmp_all.Lipid_ID == name, "B_Norm"].sum() / bmp_all.B_Norm.sum())
u181, i181 = pct("LBPA(18:1_18:1)")

# ---------------- drawing helpers ----------------
fig, ax = plt.subplots(figsize=(16, 11.5))
ax.set_xlim(0, 100); ax.set_ylim(0, 100); ax.axis("off")

def region(x0, y0, x1, y1, color, label=None, ly=None):
    ax.add_patch(Rectangle((x0, y0), x1 - x0, y1 - y0, fc=color, ec="none", zorder=0))
    if label: ax.text(x0 + 0.8, ly if ly else y1 - 1.2, label, fontsize=9, color="#555", style="italic", va="top", zorder=1)

def node(x, y, text, w=15, h=4.4, fc="white", ec="#555", bold=False, fs=9):
    ax.add_patch(FancyBboxPatch((x - w / 2, y - h / 2), w, h, boxstyle="round,pad=0.25,rounding_size=0.8",
                                fc=fc, ec=ec, lw=1.1, zorder=3))
    ax.text(x, y, text, ha="center", va="center", fontsize=fs, color=INK, fontweight="bold" if bold else "normal", zorder=4)

def arrow(p, q, ls="-", color="#555", rad=0.0, lw=1.4):
    ax.add_patch(FancyArrowPatch(p, q, arrowstyle="-|>", mutation_scale=13, lw=lw, ls=ls, color=color,
                                 connectionstyle=f"arc3,rad={rad}", zorder=2, shrinkA=0, shrinkB=0))

def badge(x, y, g, w=12.6, h=2.9, left=False):
    """gene badge centred at (x, y)"""
    info = gene_info(g)
    if info is None:
        ax.add_patch(FancyBboxPatch((x - w / 2, y - h / 2), w, h, boxstyle="round,pad=0.12,rounding_size=0.5",
                                    fc="white", ec="#999", lw=1, ls="--", zorder=5))
        ax.text(x, y, f"{g}  n.d.", ha="center", va="center", fontsize=8, color="#777", zorder=6, style="italic")
        return
    v, s = info
    fc = cmap(norm(v)); lum = 0.299 * fc[0] + 0.587 * fc[1] + 0.114 * fc[2]
    ax.add_patch(FancyBboxPatch((x - w / 2, y - h / 2), w, h, boxstyle="round,pad=0.12,rounding_size=0.5",
                                fc=fc, ec="#333" if s != "ns" else "#999", lw=1.3 if s != "ns" else 0.8, zorder=5))
    ax.text(x, y, f"{g}  {v:+.2f} {s if s != 'ns' else ''}".rstrip(), ha="center", va="center", fontsize=8.2,
            color="white" if lum < 0.45 else INK, fontweight="bold" if s != "ns" else "normal", zorder=6)

def badges(x, ymid, genes, dy=3.3, **kw):
    n = len(genes)
    for i, g in enumerate(genes):
        badge(x, ymid + (n - 1) * dy / 2 - i * dy, g, **kw)

# ---------------- regions ----------------
region(0.5, 0.5, 49.5, 99.5, "#F7F7F3", "Cytosol / ER: fatty acid supply and phospholipid acyl remodeling", 99.2)
region(51, 45, 99.5, 99.5, "#FFF8E5", "ER / mitochondria: PG synthesis", 99.2)
region(51, 0.5, 99.5, 45, "#E8F1FA", "Late endosome / lysosome", 44.7)
ax.text(25, 103.2, "ACYL-CHAIN SUPPLY", ha="center", fontsize=12, fontweight="bold", color=INK)
ax.text(75, 103.2, "BMP (LBPA) BIOSYNTHESIS", ha="center", fontsize=12, fontweight="bold", color=INK)

# ---------------- left: saturated / monounsaturated path ----------------
X1 = 11
node(X1, 91, "Acetyl-CoA"); node(X1, 77, "Palmitate 16:0"); node(X1, 63, "Stearate 18:0"); node(X1, 49, "Oleate 18:1", fc="#EAF2FB", bold=True)
for a, b in [(91, 77), (77, 63), (63, 49)]: arrow((X1, a - 2.4), (X1, b + 2.4))
badges(X1 + 8.7, 84, ["Acaca", "Fasn"]); badges(X1 + 8.7, 70, ["Elovl6"]); badges(X1 + 8.7, 56, ["Scd1", "Scd2"])
ax.text(19.8, 48.9, "(Scd1/2 also\nmake 16:1)", fontsize=6.8, color="#666", style="italic", va="center")

# ---------------- left: PUFA path ----------------
XP = 34
node(XP, 91, "Linoleate 18:2\n(diet / serum)", w=16, h=5.2); node(XP, 77, "18:3", w=16); node(XP, 63, "20:3", w=16)
node(XP, 49, "Arachidonate 20:4", w=16, fc="#FDEDE3")
node(36, 36, "DHA 22:6\n(mostly uptake)", w=11.5, h=5.6, fc="#FDEDE3", bold=True, fs=8.2)
for a, b in [(88.4, 79.4), (74.6, 65.4), (60.6, 51.4)]: arrow((XP, a), (XP, b))
arrow((36, 46.6), (36, 38.9), ls="--")
badges(XP + 9, 83, ["Fads2"]); badges(XP + 9, 70, ["Elovl5"]); badges(XP + 9, 56, ["Fads1"]); badges(XP + 9, 42.8, ["Elovl2"])

# ---------------- left: activation and remodeling ----------------
node(24, 27, "Acyl-CoA pool  (16:0, 18:1, 20:4, 22:6 ...)", w=40, h=4.6)
arrow((X1, 46.6), (X1, 29.5)); arrow((26.8, 46.6), (26.8, 29.5)); arrow((36, 33.1), (36, 29.5))
badges(X1 + 8.7, 38, ["Acsl1", "Acsl3", "Acsl4"], dy=3.0)
ax.add_patch(FancyBboxPatch((3, 3.2), 42, 17.2, boxstyle="round,pad=0.25,rounding_size=0.8", fc="white", ec="#555", lw=1.1, zorder=3))
ax.text(24, 18.4, "Acyl incorporation and remodeling (PA / PG / PC / PE)", ha="center", fontsize=9, fontweight="bold", color=INK, zorder=4)
for i, g in enumerate(["Agpat3", "Lpcat3", "Mboat1"]): badge(10.5 + i * 13.5, 13.2, g)
for i, g in enumerate(["Mboat7", "Pla2g4a", "Pla2g6"]): badge(10.5 + i * 13.5, 8.2, g)
arrow((24, 24.6), (24, 20.9))

# ---------------- right: PG synthesis ----------------
XR = 68
node(XR, 91, "Phosphatidic acid (PA)", w=19); node(XR, 77, "CDP-diacylglycerol", w=19); node(XR, 63, "PG-phosphate (PGP)", w=19)
node(XR, 50, "Phosphatidylglycerol (PG)", w=21, bold=True)
for a, b in [(91, 77), (77, 63), (63, 50)]: arrow((XR, a - 2.4), (XR, b + 2.4))
badges(XR + 9.2, 84, ["Cds1", "Cds2"]); badges(XR + 9.2, 70, ["Pgs1"]); badges(XR + 9.2, 57, ["Ptpmt1"])

# ---------------- right: lysosome ----------------
node(XR, 33, "Lyso-PG (LPG)", w=17); node(XR, 18, "BMP / LBPA", w=19, h=5, fc="#E9DDF3", bold=True, fs=10)
arrow((XR, 47.6), (XR, 35.4)); arrow((XR, 30.6), (XR, 20.8))
ax.text(XR - 1.5, 40.2, "PG delivered to\nthe lysosome", fontsize=7.5, color="#555", style="italic", ha="right", va="center")
badge(XR + 9.2, 41.4, "Pla2g15"); badge(XR + 9.2, 25.8, "Cln5")
ax.text(XR + 2, 41.4, "", fontsize=7)
ax.text(XR + 5.5, 22.7, "BMP synthase\n(transacylation)", fontsize=7, color="#555", style="italic", ha="left", va="center")
arrow((XR - 4.5, 20.6), (XR - 4.5, 30.8), ls="--", color="#8A6BB0")
ax.text(XR - 5.8, 22.6, "turnover:\nPla2g15", fontsize=7, color="#8A6BB0", ha="right", style="italic", va="center")
arrow((XR + 10.6, 50.2), (XR + 9.8, 19.6), ls="--", color="#2F7F6F", rad=-0.5)
badges(88.6, 34.0, ["Pld3", "Pld4"], dy=3.2, w=11.5)
ax.text(88.6, 27.6, "proposed alternative\nPG -> BMP route", fontsize=7, color="#2F7F6F", ha="center", style="italic")

# ---------------- chain-composition link ----------------
arrow((45.4, 11.5), (58.2, 17.4), ls="--", color="#8A6BB0", rad=0.15)
ax.text(53.8, 12.9, "acyl composition\nof PG / BMP", fontsize=7.2, color="#8A6BB0", ha="center", style="italic", va="center", zorder=6)

# ---------------- lipid readouts ----------------
def readout(x, y, w, h, title, body, tag, tag_color):
    ax.add_patch(FancyBboxPatch((x - w / 2, y - h / 2), w, h, boxstyle="round,pad=0.25,rounding_size=0.8",
                                fc="#F3ECF9", ec="#6A3D9A", lw=1.3, zorder=3))
    ax.text(x, y + h / 2 - 0.9, title, ha="center", va="top", fontsize=8.2, fontweight="bold", color="#4A2A70", zorder=4)
    ax.text(x, y + 0.3, body, ha="center", va="center", fontsize=7.6, color=INK, zorder=4)
    ax.text(x, y - h / 2 + 1.3, tag, ha="center", va="center", fontsize=7.4, color=tag_color, fontweight="bold", zorder=4)
ax.text(75.5, 13.2, "Pilot lipidomics, BMP / LBPA (IFNB vs untreated)", ha="center", fontsize=8.8, fontweight="bold", color="#4A2A70")
readout(59.5, 5.4, 15, 9.2, "Total BMP signal", f"{total_lfc:+.2f} log2FC", "synthesis / turnover genes mixed", "#666")
readout(75.5, 5.4, 15, 9.2, "BMP with 18:1", f"{f181[0]:+.2f} vs {f181[1]:+.2f} without (n={f181[2]})\n18:1/18:1: {u181:.0f}% -> {i181:.0f}% of pool", "matches Scd1/2 down", "#2F7F6F")
readout(91.5, 5.4, 15, 9.2, "BMP with 22:6 (DHA)", f"{fdha[0]:+.2f} vs {fdha[1]:+.2f} without (n={fdha[2]})", "opposes PUFA-supply genes down", "#B2382C")

# ---------------- colour bar & footer ----------------
cax = fig.add_axes([0.40, 0.085, 0.2, 0.012])
cb = fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap), cax=cax, orientation="horizontal")
cb.set_label("enzyme transcript log2FC, IFNB vs Mock  (bold = padj < 0.05; * < 0.05, ** < 0.01, *** < 0.001)", fontsize=8.5)
cb.ax.tick_params(labelsize=8)
fig.text(0.5, 0.02, "Simplified pathway. Dashed box = not detected in the RNA-seq. mRNA level is not enzyme activity; lipid values are one measurement per condition.",
         ha="center", fontsize=8, color="#555")
fig.savefig(args.out, dpi=200, bbox_inches="tight", facecolor="white")
print("saved", args.out)
