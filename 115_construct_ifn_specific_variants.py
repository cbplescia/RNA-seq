#!/usr/bin/env python3
"""Construct variants that belong to the interferon-treated conditions (IFNB and/or IFNG) and not to the matching untreated virus conditions   (reads results/ from 113_construct_variants.py).
Comparison set: samples that carry the construct (>= MIN_DEPTH reads at that position). IFN group = samples with IFNB or IFNG (e.g. IFNB+EBOV, IFNG+MARV); reference group = construct-bearing samples without IFN (EBOV, MARV, EVLP ...).
For every variant (positions where the SAME alternative base is called, per 113):
  IFN-specific : called in >= MIN_CARRIERS IFN samples and in NONE of the >= MIN_REF covered untreated construct samples
  IFN-enriched : not strictly specific, but carried by a larger share of IFN samples than of untreated samples (Fisher exact p < P_CUT)
  tag IFNB-only / IFNG-only / both IFNs tells which treatment the carriers come from.
Figures per construct gene (shared_ifn/ifn_variants_<gene>.png/.pdf): protein lollipop (x = residue, y = carriers among IFN samples, labelled with the residue change, e.g. V34I) and a heatmap variants x samples (columns grouped untreated | IFNB | IFNG;
tile = variant frequency, dot = covered but not called, grey = too little coverage).  Plus ifn_variants_overview (counts per gene) and ifn_specific_variants.tsv (all tested variants: carriers, coverage, frequencies, Fisher p).
CAUTION: a variant seen only in IFN samples can reflect IFN biology (e.g. ADAR/APOBEC editing, selection), different virus/construct batches, or simple chance with few samples; check the carrier counts and that the IFN carriers are not a single replicate set.
Options (environment): FOCUS=model_GP:104 (show just these residues - gene:residue, comma-separated - even if they fail every filter below, e.g. carried by one sample only; combine with VIRUS=EBOV MIN_CARRIERS=1)  VIRUS=EBOV (only the samples infected with that virus, e.g. EBOV: untreated EBOV vs IFNB+EBOV vs IFNG+EBOV; MARV etc. are ignored; output goes to shared_ifn_EBOV/)  ALL_SHARED=1 (also show variants shared among the IFN samples that are
  ALSO present in untreated samples, marked with open squares)  MIN_DEPTH=30 MIN_CARRIERS=2 MIN_REF=2 P_CUT=0.05 SHOW_ENRICHED=1 MAX_ROWS=40 PROJ=...        Usage: python 115_construct_ifn_specific_variants.py     (numpy, pandas, matplotlib)
"""
import os, sys, re, math
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.patches import Patch
from matplotlib.lines import Line2D

PROJ = os.environ.get("PROJ", os.path.join(os.path.expanduser("~"), "Documents/aviti_rnaseq/2025/aviti/20251009-100725_Maury_25161_WivfFacICOTHFAUJwSPeHBEFfCfNtCRpUSgVGUkz/results/Maury_25161"))
RES = os.path.join(PROJ, "construct_variants", "results"); VIRUS = os.environ.get("VIRUS", "").strip(); FOCUS = [t.strip() for t in os.environ.get("FOCUS", "").split(",") if t.strip()]     # FOCUS=model_GP:104,model_NP:50 -> show only these residues, whatever the filters say
ALL_SHARED = os.environ.get("ALL_SHARED", "0") == "1"
OUT = os.path.join(RES, "shared_ifn" + ("_" + VIRUS if VIRUS else "")); os.makedirs(OUT, exist_ok=True)
MIN_DEPTH = int(os.environ.get("MIN_DEPTH", "30")); MIN_CARRIERS = int(os.environ.get("MIN_CARRIERS", "2")); MIN_REF = int(os.environ.get("MIN_REF", "2")); P_CUT = float(os.environ.get("P_CUT", "0.05"))
SHOW_ENRICHED = os.environ.get("SHOW_ENRICHED", "1") == "1"; MAX_ROWS = int(os.environ.get("MAX_ROWS", "40")); INK = "#1F2933"
ORDER = ["Mock_Mock", "IFNB_Mock", "IFNY_Mock", "Mock_EBOV", "IFNB_EBOV", "IFNY_EBOV", "Mock_MARV", "IFNB_MARV", "IFNY_MARV", "Mock_EVLP", "Mock_MVLP"]
LAB = {"Mock_Mock": "Mock", "IFNB_Mock": "IFNB", "IFNY_Mock": "IFNG", "Mock_EBOV": "EBOV", "IFNB_EBOV": "IFNB+EBOV", "IFNY_EBOV": "IFNG+EBOV", "Mock_MARV": "MARV", "IFNB_MARV": "IFNB+MARV", "IFNY_MARV": "IFNG+MARV", "Mock_EVLP": "EVLP", "Mock_MVLP": "MVLP"}
TAGC = {"IFNB-only": "#2B6CB0", "IFNG-only": "#D97706", "IFNB + IFNG": "#7B3FA0"}

def fisher_p(a, b, c, d):
    """two-sided Fisher exact p for [[a,b],[c,d]] (standard library only)"""
    r1, r2, c1, n = a + b, c + d, a + c, a + b + c + d
    if n == 0: return 1.0
    def pr(x): return math.comb(r1, x) * math.comb(r2, c1 - x) / math.comb(n, c1)
    lo, hi = max(0, c1 - r2), min(r1, c1); p0 = pr(a)
    return min(1.0, sum(pr(x) for x in range(lo, hi + 1) if pr(x) <= p0 * (1 + 1e-9)))

vp = os.path.join(RES, "variants.tsv")
if not os.path.exists(vp) or os.path.getsize(vp) < 5: sys.exit("No variants.tsv - run 113_construct_variants.py first")
V = pd.read_csv(vp, sep="\t", dtype={"feature": str}); meta = ["feature", "pos", "ref", "alt", "n_samples", "mean_vaf", "max_vaf", "min_vaf", "n_covered", "frac_samples", "aa_change", "class"]
samples = [c for c in V.columns if c not in meta]
si = pd.read_csv(os.path.join(PROJ, "sample_info.csv"), sep=None, engine="python", dtype=str)
si["ifn"] = si["condition1"].fillna("Mock").replace({"None": "Mock", "none": "Mock"}); si["virus"] = si["condition2"].fillna("Mock").replace({"None": "Mock", "none": "Mock"}); si["group"] = si["ifn"] + "_" + si["virus"]
ifn = dict(zip(si["sample_id"], si["ifn"])); grp = dict(zip(si["sample_id"], si["group"])); vir = dict(zip(si["sample_id"], si["virus"]))
if VIRUS: samples = [s_ for s_ in samples if vir.get(s_) == VIRUS]; print("restricting to %s samples: %s" % (VIRUS, ", ".join(samples)))
def ifn_class(s): return "IFNB" if ifn.get(s) == "IFNB" else "IFNG" if ifn.get(s) in ("IFNY", "IFNG") else "none"
samples = sorted(samples, key=lambda s: ({"none": 0, "IFNB": 1, "IFNG": 2}[ifn_class(s)], ORDER.index(grp[s]) if grp.get(s) in ORDER else 99, s)); cond = {s: LAB.get(grp.get(s, "?"), grp.get(s, "?")) for s in samples}
refs, n = {}, None
for line in open(os.path.join(PROJ, "ref", "recombinant_model.fa")):
    line = line.strip()
    if line.startswith(">"): n = line[1:].split()[0]; refs[n] = []
    elif n: refs[n].append(line.upper())
refs = {k: "".join(v) for k, v in refs.items()}; is_cds = {k: (len(s) % 3 == 0 and s.startswith("ATG")) for k, s in refs.items()}
cov = {}
for f in refs:
    p = os.path.join(RES, "per_base_counts_%s.tsv" % f)
    if os.path.exists(p):
        t = pd.read_csv(p, sep="\t"); cov[f] = {s: (t[s + "_depth"].values >= MIN_DEPTH) for s in samples if (s + "_depth") in t.columns}
def residue(r):
    m = re.match(r"([A-Z\*\?])(\d+)([A-Z\*\?]) \((\w+)\)", r.aa_change) if isinstance(r.aa_change, str) else None
    if m: return int(m.group(2)), "%s%s%s" % (m.group(1), m.group(2), m.group(3)), m.group(4)
    kind = r.alt if r.alt in ("del", "ins") else "%s>%s" % (r.ref, r.alt)
    return ((r.pos - 1) // 3 + 1 if is_cds[r.feature] else r.pos), "c.%d%s" % (r.pos, kind), "non-coding / indel"

rows = []
for _, r in V.iterrows():
    f = r["feature"]
    if f not in cov: continue
    c = cov[f]; pos = int(r["pos"]) - 1
    ci = [s for s in samples if ifn_class(s) != "none" and c.get(s) is not None and c[s][pos]]; cr = [s for s in samples if ifn_class(s) == "none" and c.get(s) is not None and c[s][pos]]
    ki = [s for s in ci if pd.notna(r[s])]; kr = [s for s in cr if pd.notna(r[s])]
    if not ki: continue
    a, b, cc, dd = len(ki), len(ci) - len(ki), len(kr), len(cr) - len(kr); p = fisher_p(a, b, cc, dd)
    tags = {ifn_class(s) for s in ki}; tag = "IFNB + IFNG" if len(tags) == 2 else ("IFNB-only" if "IFNB" in tags else "IFNG-only")
    cls = "IFN-specific" if (len(ki) >= MIN_CARRIERS and len(kr) == 0 and len(cr) >= MIN_REF) else ("IFN-enriched" if (len(ci) and len(cr) and a / len(ci) > cc / len(cr) and p < P_CUT and len(ki) >= MIN_CARRIERS) else ("IFN-shared (also untreated)" if len(ki) >= MIN_CARRIERS else ""))
    x, lab, eff = residue(type("R", (), dict(feature=f, pos=int(r["pos"]), ref=r["ref"], alt=r["alt"], aa_change=r["aa_change"]))())
    rows.append(dict(feature=f, pos=int(r["pos"]), ref=r["ref"], alt=r["alt"], x=x, label=lab, effect=eff, IFN_carriers=len(ki), IFN_covered=len(ci), untreated_carriers=len(kr), untreated_covered=len(cr), fisher_p=p, ifn_tag=tag, ifn_class=cls,
                     mean_vaf_IFN=float(np.nanmean([r[s] for s in ki])), carriers_IFN=", ".join("%s:%s" % (cond[s], s) for s in ki), carriers_untreated=", ".join("%s:%s" % (cond[s], s) for s in kr)))
T = pd.DataFrame(rows)
if not len(T): sys.exit("No variants were called in any IFN sample with coverage")
T.sort_values(["feature", "pos"]).to_csv(os.path.join(OUT, "ifn_specific_variants.tsv"), sep="\t", index=False, float_format="%.4g")
focus_ok = lambda r: ("%s:%d" % (r.feature, r.x)) in FOCUS
S = T[T.apply(focus_ok, axis=1)] if FOCUS else T[T.ifn_class.isin(["IFN-specific"] + (["IFN-enriched"] if SHOW_ENRICHED else []) + (["IFN-shared (also untreated)"] if ALL_SHARED else []))]
print("tested %d variants called in IFN samples; IFN-specific: %d, IFN-enriched: %d, shared by >= %d IFN samples but also in untreated: %d" % (len(T), (T.ifn_class == "IFN-specific").sum(), (T.ifn_class == "IFN-enriched").sum(), MIN_CARRIERS, (T.ifn_class == "IFN-shared (also untreated)").sum()))
for f in refs: print("  %-24s specific %d  enriched %d" % (f, ((T.feature == f) & (T.ifn_class == "IFN-specific")).sum(), ((T.feature == f) & (T.ifn_class == "IFN-enriched")).sum()))
def style(ax): ax.tick_params(colors=INK, labelsize=8); [sp.set_color("#9AA5B1") for sp in ax.spines.values()]
def save(fig, name): fig.savefig(os.path.join(OUT, name + ".png"), dpi=300, bbox_inches="tight"); fig.savefig(os.path.join(OUT, name + ".pdf"), bbox_inches="tight"); plt.close(fig)
cm = plt.get_cmap("YlOrRd").copy(); cm.set_bad((0, 0, 0, 0))
Vi = V.set_index(["feature", "pos", "alt"])

for f, seq in refs.items():
    d = S[S.feature == f].sort_values("pos")
    if not len(d) or f not in cov: continue
    cs = [s for s in samples if cov[f].get(s) is not None and cov[f][s].any()]
    L = len(seq) // 3 if is_cds[f] else len(seq); unit = "residue" if is_cds[f] else "nucleotide"
    pages = [d.iloc[i:i + MAX_ROWS] for i in range(0, len(d), MAX_ROWS)]
    for pi, dp in enumerate(pages):
        fig = plt.figure(figsize=(max(9, 0.34 * len(cs) + 6), 3.6 + 0.24 * len(dp) + 1.3)); gs = fig.add_gridspec(2, 2, height_ratios=[2.3, 0.24 * len(dp) + 0.9], width_ratios=[40, 1], hspace=0.75, wspace=0.03)
        a1 = fig.add_subplot(gs[0, 0]); a2 = fig.add_subplot(gs[1, 0]); cax = fig.add_subplot(gs[1, 1])
        a1.hlines(0, 0, L + 1, color="#9AA5B1", lw=4, zorder=1)
        for r in dp.itertuples():
            filled = r.ifn_class == "IFN-specific"; mk = "o" if filled else ("D" if r.ifn_class == "IFN-enriched" else "s"); a1.vlines(r.x, 0, r.IFN_carriers, color=TAGC[r.ifn_tag], lw=1.0, zorder=2)
            a1.scatter(r.x, r.IFN_carriers, s=40, marker=mk, facecolor=TAGC[r.ifn_tag] if filled else "white", edgecolor=TAGC[r.ifn_tag], linewidth=1.2, zorder=3)
        for k, r in enumerate(dp.sort_values("x").itertuples()):
            lvl = k % 3; a1.annotate(r.label, (r.x, r.IFN_carriers), xytext=(0, 7 + 7 * lvl), textcoords="offset points", ha="center", fontsize=6.5, color=INK, arrowprops=dict(arrowstyle="-", lw=0.3, color="#9AA5B1") if lvl else None)
        nI = max(1, sum(1 for s in cs if ifn_class(s) != "none")); a1.set_xlim(0, L + 1); a1.set_ylim(0, max(nI, dp.IFN_carriers.max()) * 1.6); a1.set_ylabel("IFN samples carrying the variant", fontsize=8, color=INK)
        a1.set_xlabel("%s number in %s  (%d %ss)" % (unit, f, L, unit), fontsize=8, color=INK); style(a1)
        h = [Patch(color=c, label=t) for t, c in TAGC.items() if (dp.ifn_tag == t).any()] + [Line2D([], [], marker="o", color="grey", ls="", label="IFN-specific (absent from untreated)"), Line2D([], [], marker="D", markerfacecolor="white", color="grey", ls="", label="IFN-enriched (Fisher p < %.2g)" % P_CUT)] + ([Line2D([], [], marker="s", markerfacecolor="white", color="grey", ls="", label="shared by IFN samples, also in untreated")] if (dp.ifn_class == "IFN-shared (also untreated)").any() else []) + ([Line2D([], [], marker="s", markerfacecolor="white", color="grey", ls="", label="focus residue (fails the specificity filters)")] if (dp.ifn_class == "").any() else [])
        a1.legend(handles=h, fontsize=6.5, frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.2), ncol=3); a1.yaxis.set_major_locator(plt.matplotlib.ticker.MaxNLocator(integer=True))
        a1.set_title("%s%s: variants shared by IFNB / IFNG samples%s" % (f, " (%s only)" % VIRUS if VIRUS else "", "  (page %d of %d)" % (pi + 1, len(pages)) if len(pages) > 1 else ""), fontsize=10, color=INK, fontweight="bold", loc="left")
        M = np.full((len(dp), len(cs)), np.nan); C = np.zeros((len(dp), len(cs)), dtype=bool)
        for i, r in enumerate(dp.itertuples()):
            row = Vi.loc[(f, r.pos, r.alt)]; row = row.iloc[0] if isinstance(row, pd.DataFrame) else row
            for j, s in enumerate(cs):
                C[i, j] = bool(cov[f][s][r.pos - 1]); M[i, j] = row[s] if pd.notna(row[s]) else np.nan
        a2.imshow(np.ma.masked_where(C, np.ones(M.shape)), aspect="auto", cmap=plt.matplotlib.colors.ListedColormap(["#E5E7EB"]), interpolation="nearest")      # grey = too little coverage
        im = a2.imshow(np.ma.masked_invalid(M), aspect="auto", cmap=cm, norm=LogNorm(vmin=0.02, vmax=1.0), interpolation="nearest")
        for i in range(len(dp)):
            for j in range(len(cs)):
                if np.isfinite(M[i, j]): a2.text(j, i, "%d" % max(1, round(100 * M[i, j])), ha="center", va="center", fontsize=5, color="white" if M[i, j] > 0.35 else INK)
                elif C[i, j]: a2.text(j, i, "·", ha="center", va="center", fontsize=8, color="#9AA5B1")
        a2.set_yticks(range(len(dp))); a2.set_yticklabels(["%s  (%s; %d/%d IFN vs %d/%d untreated)" % (r.label, r.ifn_tag, r.IFN_carriers, r.IFN_covered, r.untreated_carriers, r.untreated_covered) for r in dp.itertuples()], fontsize=6.5, color=INK)
        a2.set_xticks(range(len(cs))); a2.set_xticklabels(["%s %s" % (cond[s], s) for s in cs], rotation=75, ha="right", fontsize=6.5)
        for tl, s in zip(a2.get_xticklabels(), cs): tl.set_color({"none": "#374151", "IFNB": TAGC["IFNB-only"], "IFNG": TAGC["IFNG-only"]}[ifn_class(s)]); tl.set_fontweight("bold")
        edges = [j - 0.5 for j in range(1, len(cs)) if ifn_class(cs[j]) != ifn_class(cs[j - 1])]
        for e in edges: a2.axvline(e, color=INK, lw=1.2)
        a2.set_xticks(np.arange(-0.5, len(cs), 1), minor=True); a2.set_yticks(np.arange(-0.5, len(dp), 1), minor=True); a2.grid(which="minor", color="white", lw=0.5); a2.tick_params(which="minor", length=0); style(a2)
        cb = fig.colorbar(im, cax=cax); cb.set_label("variant frequency", fontsize=7); cb.ax.tick_params(labelsize=6)
        a2.set_title("columns: untreated construct samples | IFNB | IFNG   (numbers = frequency %%, dot = covered but not called, grey = <%d reads)" % MIN_DEPTH, fontsize=7.5, color="grey", loc="left")
        save(fig, "ifn_variants_%s%s" % (f, "" if len(pages) == 1 else "_p%d" % (pi + 1)))

cnt = S.groupby(["feature", "ifn_tag"]).size().unstack(fill_value=0).reindex(list(refs)).fillna(0)
for t in TAGC:
    if t not in cnt: cnt[t] = 0
cnt = cnt[list(TAGC)]; fig, ax = plt.subplots(figsize=(max(6, 0.6 * len(cnt) + 3), 4)); bot = np.zeros(len(cnt))
for t, c in TAGC.items(): ax.bar(range(len(cnt)), cnt[t].values, bottom=bot, color=c, label=t); bot += cnt[t].values
ax.set_xticks(range(len(cnt))); ax.set_xticklabels(cnt.index, rotation=60, ha="right", fontsize=8, color=INK); ax.set_ylabel("IFN-specific%s variants" % (" / enriched" if SHOW_ENRICHED else ""), fontsize=8, color=INK); ax.legend(fontsize=7, frameon=False); style(ax)
ax.set_title("Construct variants shared by the interferon conditions%s, per gene" % (" (%s only)" % VIRUS if VIRUS else ""), fontsize=10, color=INK, fontweight="bold", loc="left"); fig.tight_layout(); save(fig, "ifn_variants_overview")
print("\nWrote figures + ifn_specific_variants.tsv to", OUT)
