#!/usr/bin/env python3
"""Which transcription-factor motifs mark IFNB-induced genes that EBOV does NOT suppress, compared with those it fully suppresses or partly restricts?

Input: requant_v2_outputs/figures_ifnb_ebov_retention/ifnb_induced_groups.csv from 107_ifnb_induced_gene_groups.R (groups: fully suppressed / partly restricted / not suppressed).
A BACKGROUND set (random protein-coding genes that IFNB does not induce, taken from transcriptome_blocks.csv) is scanned too, so generic motifs (GC-rich, SP-like) are not mistaken for IFN-specific ones.
For every gene: the promoter of its canonical transcript (UP bp upstream .. DOWN bp downstream of the TSS, default 1500 / 500) is cut from your GRCm39 reference (GTF + genome FASTA), scanned on both strands with JASPAR
vertebrate matrices (relative score >= REL, default 0.88; JASPAR download, JASPAR_FILE=..., or a small built-in IFN panel if neither works). A motif counts as "present" in a promoter if it has >= 1 site.
Comparison: for each motif, the share of promoters with a site in each group; Fisher exact tests (not suppressed vs fully suppressed = the main question; each group vs background), Benjamini-Hochberg FDR across motifs.
Outputs in requant_v2_outputs/promoter_tf_groups/:
  tf_group_comparison.tsv        every motif: % of promoters with a site per group, odds ratios, p values, FDR
  promoter_motif_best_scores.tsv gene x motif best relative score (the raw scan)       genes_scanned.tsv   the genes, groups and promoter coordinates
  tf_group_dotplot.png / .pdf    the most informative TFs x groups (colour = % promoters with a site), with the TF's own expression (IFNB, EBOV, IFNB+EBOV vs Mock) alongside
  tf_volcano_not_vs_full.png     effect size and significance of every TF motif, not suppressed vs fully suppressed
  tf_family_by_group.png         the same comparison at the level of TF families (IRF/STAT, NF-kB, AP-1/CREB, SP/KLF, C/EBP, ETS, other)
IMPORTANT: motif presence is a possible-binding-site signal, not evidence of binding; a promoter-only scan misses enhancers; with few genes per group, differences are suggestive. Groups are small, so look at
effect sizes and the dot plot together, not at the FDR alone. Real evidence needs ChIP-seq / ATAC data from mouse macrophages (public tracks are mostly mm10; your reference is GRCm39).
Usage:  python 108_promoter_tf_groups.py     (environment options: UP=1500 DOWN=500 REL=0.88 N_BG=300 SEED=1 JASPAR_FILE=path PROJ=path)     Needs python 3, numpy, scipy, matplotlib.  Runs several minutes (about 1 s per 3 promoters).
"""
import os, re, sys, csv, gzip, urllib.request, random
import numpy as np
from scipy.stats import fisher_exact

PROJ = os.environ.get("PROJ", os.path.join(os.path.expanduser("~"), "Documents/aviti_rnaseq/2025/aviti/20251009-100725_Maury_25161_WivfFacICOTHFAUJwSPeHBEFfCfNtCRpUSgVGUkz/results/Maury_25161"))
RNA = os.path.join(PROJ, "requant_v2_outputs")
GTF = os.environ.get("GTF", os.path.join(PROJ, "ref", "Mus_musculus.GRCm39.109.gtf")); GENOME = os.environ.get("GENOME", os.path.join(PROJ, "ref", "Mus_musculus.GRCm39.dna.primary_assembly.fa"))
UP, DOWN, REL = int(os.environ.get("UP", "1500")), int(os.environ.get("DOWN", "500")), float(os.environ.get("REL", "0.88")); N_BG, SEED = int(os.environ.get("N_BG", "300")), int(os.environ.get("SEED", "1")); MIN_COUNT = int(os.environ.get("MIN_COUNT", "3")); LOR = float(os.environ.get("LOR", "1")); N_LAB = int(os.environ.get("N_LAB", "10"))
OUT = os.path.join(RNA, "promoter_tf_groups"); os.makedirs(OUT, exist_ok=True)
GROUPS = ["fully suppressed", "partly restricted", "not suppressed", "background"]
GCOL = {"fully suppressed": "#762A83", "partly restricted": "#E08214", "not suppressed": "#1B7837", "background": "#9AA0A6"}
JASPAR_URLS = ["https://jaspar.elixir.no/download/data/2024/CORE/JASPAR2024_CORE_vertebrates_non-redundant_pfms_jaspar.txt", "https://jaspar.genereg.net/download/data/2024/CORE/JASPAR2024_CORE_vertebrates_non-redundant_pfms_jaspar.txt",
               "https://jaspar2022.genereg.net/download/data/2022/CORE/JASPAR2022_CORE_vertebrates_non-redundant_pfms_jaspar.txt"]
COMP = str.maketrans("ACGTacgtN", "TGCAtgcaN")
def revcomp(s): return s.translate(COMP)[::-1]
def say(*a): print(*a, flush=True)

# ---- genome / GTF (as in 106) ----
def build_fai(path):
    fai = path + ".fai"
    if os.path.exists(fai): return fai
    say(f"building the FASTA index {fai} (about a minute) ...")
    rows, name, length, lw, lb, start, pos = [], None, 0, None, None, 0, 0
    with open(path, "rb") as fh:
        for line in fh:
            n = len(line)
            if line.startswith(b">"):
                if name is not None: rows.append((name, length, start, lw, lb))
                name, length, start, lw, lb = line[1:].split()[0].decode(), 0, pos + n, None, None
            else:
                s = len(line.rstrip(b"\r\n"))
                if lw is None: lw, lb = s, n
                length += s
            pos += n
        if name is not None: rows.append((name, length, start, lw, lb))
    with open(fai, "w") as o:
        for r in rows: o.write("\t".join(map(str, r)) + "\n")
    return fai
class Fasta:
    def __init__(self, path):
        self.fh = open(path, "rb"); self.idx = {}
        for l in open(build_fai(path)):
            n, ln, st, lw, lb = l.split()[:5]; self.idx[n] = (int(ln), int(st), int(lw), int(lb))
    def get(self, chrom, start, end):
        ln, st, lw, lb = self.idx[chrom]; start, end = max(1, start), min(ln, end); a, b = start - 1, end
        off = st + (a // lw) * lb + a % lw; nbytes = (b - a) + ((b - 1) // lw - a // lw) * (lb - lw) + (lb - lw)
        self.fh.seek(off); return re.sub(r"\s", "", self.fh.read(nbytes).decode())[: b - a].upper()

def read_gtf(names):                                           # one pass: transcripts of the wanted gene symbols
    want = {n.lower() for n in names}; tx = {}
    op = gzip.open if GTF.endswith(".gz") else open
    with op(GTF, "rt") as fh:
        for line in fh:
            if line.startswith("#"): continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "transcript": continue
            m = re.search(r'gene_name "([^"]+)"', f[8])
            if not m or m.group(1).lower() not in want: continue
            bt = re.search(r'transcript_biotype "([^"]+)"', f[8])
            tx.setdefault(m.group(1).lower(), []).append(dict(chrom=f[0], strand=f[6], tss=int(f[3]) if f[6] == "+" else int(f[4]), biotype=bt.group(1) if bt else "", canon='tag "Ensembl_canonical"' in f[8]))
    return tx
def promoter_of(tl):
    pc = [t for t in tl if t["biotype"] in ("protein_coding", "")] or tl; cl = []
    for t in sorted(pc, key=lambda x: x["tss"]):
        for c in cl:
            if abs(c["tss"] - t["tss"]) <= 200: c["m"].append(t); break
        else: cl.append(dict(tss=t["tss"], m=[t]))
    for c in cl: c["canon"] = any(m["canon"] for m in c["m"]); c["tss"] = int(np.median([m["tss"] for m in c["m"]]))
    c = sorted(cl, key=lambda c: (-int(c["canon"]), -len(c["m"])))[0]; t0 = c["m"][0]
    return dict(chrom=t0["chrom"], strand=t0["strand"], tss=c["tss"])

# ---- motifs (as in 106) ----
def parse_jaspar(text):
    ms, cur = [], None
    for line in text.splitlines():
        line = line.strip()
        if line.startswith(">"): p = line[1:].split(None, 1); cur = dict(id=p[0], name=p[1] if len(p) > 1 else p[0], rows=[]); ms.append(cur)
        elif cur and line and line[0] in "ACGT": cur["rows"].append([float(x) for x in re.findall(r"[-\d.]+", line.split("[", 1)[1])])
    return [(m["id"], m["name"], np.array(m["rows"]).T) for m in ms if len(m["rows"]) == 4]
IUPAC = {"A": "A", "C": "C", "G": "G", "T": "T", "R": "AG", "Y": "CT", "S": "CG", "W": "AT", "K": "GT", "M": "AC", "N": "ACGT"}
BUILTIN = {"ISRE_IRF9/STAT1/STAT2": "GAAANNGAAACT", "IRF_ISRE-like": "GAAANNGAAA", "GAS_STAT1": "TTCNNNGAA", "NFKB": "GGGRNNYYCC", "SP1": "GGGCGG", "AP1": "TGASTCA", "CREB/ATF": "TGACGTCA", "CEBP": "TTGCGCAA", "ETS": "GGAAGT", "GATA": "WGATAR"}
def builtin_motifs():
    out = []
    for name, pat in BUILTIN.items():
        c = np.zeros((len(pat), 4))
        for i, ch in enumerate(pat):
            for b in IUPAC[ch]: c[i, "ACGT".index(b)] = 20.0 / len(IUPAC[ch])
        out.append((f"builtin:{name}", name, c + 0.1))
    return out
def load_motifs():
    jf = os.environ.get("JASPAR_FILE")
    if jf and os.path.exists(jf): return parse_jaspar(open(jf).read()), f"JASPAR file {jf}"
    cache = os.path.join(RNA, "promoter_scan", "jaspar_pfms.txt")
    for c in (cache, os.path.join(OUT, "jaspar_pfms.txt")):
        if os.path.exists(c) and os.path.getsize(c) > 100000: return parse_jaspar(open(c).read()), f"JASPAR (cached {c})"
    for u in JASPAR_URLS:
        try:
            txt = urllib.request.urlopen(u, timeout=60).read().decode(); ms = parse_jaspar(txt)
            if len(ms) > 100: open(os.path.join(OUT, "jaspar_pfms.txt"), "w").write(txt); return ms, f"JASPAR download ({u.split('/')[-1]})"
        except Exception as e: say(f"  could not get {u.split('/')[2]}: {type(e).__name__}")
    say("JASPAR could not be downloaded - using the small built-in panel (rerun with JASPAR_FILE=... for the full set)")
    return builtin_motifs(), "built-in IFN-relevant panel"
def prep(motifs):
    out = []
    for mid, name, counts in motifs:
        n = counts.sum(1, keepdims=True); W = np.log2(((counts + 0.2) / (n + 0.8)) / 0.25); mn, mx = W.min(1).sum(), W.max(1).sum()
        if mx > mn: out.append((mid, name, W, mn, mx))
    return out
def best_scores(seq, pm):
    enc = np.array([("ACGT".find(c)) for c in seq]); L = len(seq); res = np.full(len(pm), -1.0)
    for k, (mid, name, W, mn, mx) in enumerate(pm):
        w = W.shape[0]
        if L < w: continue
        for e in (enc, np.where(enc >= 0, 3 - enc, -1)[::-1]):
            s = np.zeros(L - w + 1); bad = np.zeros(L - w + 1, bool)
            for j in range(w):
                col = e[j:L - w + 1 + j]; bad |= col < 0; s += W[j, np.where(col < 0, 0, col)]
            r = (s - mn) / (mx - mn); r[bad] = -1; res[k] = max(res[k], r.max())
    return res

def bh(p):
    p = np.array(p, float); n = len(p); o = np.argsort(p); q = np.empty(n); prev = 1.0
    for rank, i in enumerate(o[::-1]): k = n - rank; prev = min(prev, p[i] * n / k); q[i] = prev
    return q

def main():
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    gf = os.path.join(RNA, "figures_ifnb_ebov_retention", "ifnb_induced_groups.csv")
    if not os.path.exists(gf): sys.exit(f"Run 107_ifnb_induced_gene_groups.R first (missing {gf})")
    for p in (GTF, GENOME):
        if not os.path.exists(p): sys.exit(f"Missing {p}")
    genes = {}
    for r in csv.DictReader(open(gf)):
        if r["group"] in GROUPS and r["gene_name"]: genes[r["gene_name"]] = r["group"]
    blocks = os.path.join(RNA, "figures_mock_ifnb_ebov", "transcriptome_blocks.csv")
    if not os.path.exists(blocks): sys.exit(f"Missing {blocks} (script 51 output)")
    pool = [r["gene_name"] for r in csv.DictReader(open(blocks)) if r.get("resp") == "None" and r.get("gene_name") and r["gene_name"] not in genes]
    random.seed(SEED); bg = random.sample(pool, min(N_BG, len(pool)))
    for g in bg: genes[g] = "background"
    say("genes per group: " + ", ".join(f"{g} {sum(1 for v in genes.values() if v == g)}" for g in GROUPS))
    tx = read_gtf(genes); fa = Fasta(GENOME); motifs, src = load_motifs(); pm = prep(motifs); say(f"{len(pm)} motifs from {src}; REL >= {REL}; promoter {UP} up / {DOWN} down")
    names, scores, info = [], [], []
    for i, (g, grp) in enumerate(genes.items()):
        tl = tx.get(g.lower())
        if not tl: continue
        p = promoter_of(tl)
        gs, ge = (p["tss"] - UP, p["tss"] + DOWN - 1) if p["strand"] == "+" else (p["tss"] - DOWN + 1, p["tss"] + UP)
        if p["chrom"] not in fa.idx: continue
        seq = fa.get(p["chrom"], gs, ge); seq = seq if p["strand"] == "+" else revcomp(seq)
        if len(seq) < 200: continue
        names.append(g); scores.append(best_scores(seq, pm)); info.append((g, grp, p["chrom"], p["strand"], p["tss"]))
        if (len(names)) % 50 == 0: say(f"  scanned {len(names)} promoters")
    S = np.array(scores); say(f"{len(names)} promoters scanned")
    with open(f"{OUT}/genes_scanned.tsv", "w") as o:
        o.write("gene\tgroup\tchrom\tstrand\tTSS\n"); [o.write("\t".join(map(str, r)) + "\n") for r in info]
    with open(f"{OUT}/promoter_motif_best_scores.tsv", "w") as o:
        o.write("gene\tgroup\t" + "\t".join(f"{m[0]}|{m[1]}" for m in pm) + "\n")
        for r, row in zip(info, S): o.write(f"{r[0]}\t{r[1]}\t" + "\t".join(f"{x:.3f}" for x in row) + "\n")
    grp_arr = np.array([r[1] for r in info]); P = S >= REL
    n_g = {g: int((grp_arr == g).sum()) for g in GROUPS}; say("promoters per group: " + str(n_g))
    rows = []
    for k, (mid, name, *_ ) in enumerate(pm):
        c = {g: int(P[grp_arr == g, k].sum()) for g in GROUPS}
        def fe(a, b):                                            # group a vs group b
            t = [[c[a], n_g[a] - c[a]], [c[b], n_g[b] - c[b]]]; o, p = fisher_exact(t)
            lo = np.log2(((c[a] + 0.5) * (n_g[b] - c[b] + 0.5)) / ((n_g[a] - c[a] + 0.5) * (c[b] + 0.5)))        # 0.5 added to every cell, so a zero count gives a finite value
            return (lo, p)
        l1, p1 = fe("not suppressed", "fully suppressed"); l2, p2 = fe("not suppressed", "partly restricted"); l3, p3 = fe("fully suppressed", "background"); l4, p4 = fe("not suppressed", "background"); l5, p5 = fe("partly restricted", "background")
        rows.append(dict(motif_id=mid, TF=name, n_full=c["fully suppressed"], n_partial=c["partly restricted"], n_not=c["not suppressed"], n_bg=c["background"], **{f"pct_{g.split()[0]}": 100.0 * c[g] / max(1, n_g[g]) for g in GROUPS}, log2OR_not_vs_full=l1, p_not_vs_full=p1, log2OR_not_vs_partial=l2, p_not_vs_partial=p2,
                         log2OR_full_vs_bg=l3, p_full_vs_bg=p3, log2OR_partial_vs_bg=l5, p_partial_vs_bg=p5, log2OR_not_vs_bg=l4, p_not_vs_bg=p4))
    for key in ("not_vs_full", "not_vs_partial", "full_vs_bg", "partial_vs_bg", "not_vs_bg"):
        q = bh([r[f"p_{key}"] for r in rows])
        for r, v in zip(rows, q): r[f"FDR_{key}"] = v
    for r in rows:
        if r["log2OR_not_vs_full"] >= LOR and r["p_not_vs_full"] < 0.05 and r["n_not"] >= MIN_COUNT: r["class"] = "NOT-suppressed-associated"
        elif r["log2OR_not_vs_full"] <= -LOR and r["p_not_vs_full"] < 0.05 and r["n_full"] >= MIN_COUNT: r["class"] = "fully-suppressed-associated"
        elif min(r["pct_fully"], r["pct_not"]) >= 20 and abs(r["log2OR_not_vs_full"]) < 0.5: r["class"] = "shared (common to both)"
        else: r["class"] = "no clear difference"
    rows.sort(key=lambda r: r["p_not_vs_full"]); keys = list(rows[0].keys())
    with open(f"{OUT}/tf_group_comparison.tsv", "w") as o:
        o.write("\t".join(keys) + "\n"); [o.write("\t".join(f"{r[k]:.4g}" if isinstance(r[k], float) else str(r[k]) for k in keys) + "\n") for r in rows]
    say("top motifs, not suppressed vs fully suppressed (positive log2OR = more common in NOT-suppressed genes):")
    for r in rows[:12]: say(f"   {r['TF']:<28s} {r['pct_fully']:5.1f}% (full) {r['pct_partly']:5.1f}% (partial) {r['pct_not']:5.1f}% (not) {r['pct_background']:5.1f}% (bg)  log2OR {r['log2OR_not_vs_full']:+.2f}  p {r['p_not_vs_full']:.3g}  FDR {r['FDR_not_vs_full']:.2g}")
    # ---- expression of the TFs ----
    expr = {}
    for nm, col in (("IFNB_vs_Mock", "IFNB"), ("EBOV_vs_Mock", "EBOV"), ("IFNB_EBOV_vs_Mock", "IFNB+EBOV")):
        f = os.path.join(RNA, f"DE_v2_{nm}_annotated.csv")
        if os.path.exists(f):
            for r in csv.DictReader(open(f)): expr.setdefault(r["gene_name"].lower(), {})[col] = float(r["log2FoldChange"]) if r["log2FoldChange"] not in ("", "NA") else np.nan
    def tfexpr(name):
        for tok in re.split(r"::", re.sub(r"\(var\.\d+\)", "", name)):
            if tok.strip().lower() in expr: return expr[tok.strip().lower()]
        return None
    # ---- dotplot of the most informative TFs (one motif per TF name) ----
    seen, pick = set(), []
    def add(lst):
        for r in lst:
            if r["TF"] not in seen: seen.add(r["TF"]); pick.append(r)
    add(rows[:24]); add(sorted(rows, key=lambda r: r["p_not_vs_bg"])[:8]); add(sorted(rows, key=lambda r: r["p_full_vs_bg"])[:8])
    pick = sorted(pick, key=lambda r: -r["log2OR_not_vs_full"])[:40]
    fig = plt.figure(figsize=(9.5, 0.27 * len(pick) + 2.6)); gs = fig.add_gridspec(1, 3, width_ratios=[4, 0.25, 2.2], wspace=0.05)
    ax = fig.add_subplot(gs[0]); ax2 = fig.add_subplot(gs[2], sharey=ax)
    keysg = ["pct_fully", "pct_partly", "pct_not", "pct_background"]
    for j, kk in enumerate(keysg):
        ax.scatter([j] * len(pick), range(len(pick)), c=[r[kk] for r in pick], cmap="viridis", vmin=0, vmax=max(40, max(r[kk] for r in pick)), s=120, edgecolors="none")
    for i, r in enumerate(pick):
        if r["FDR_not_vs_full"] < 0.05: ax.text(3.45, i, "*", va="center", fontsize=10)
        elif r["p_not_vs_full"] < 0.05: ax.text(3.45, i, "+", va="center", fontsize=9, color="grey")
    ax.set_xticks(range(4)); ax.set_xticklabels([f"{g} (n={n_g[g]})" for g in GROUPS], fontsize=8, rotation=25, ha="right"); ax.set_yticks(range(len(pick))); ax.set_yticklabels([r["TF"] for r in pick], fontsize=7.5); ax.invert_yaxis(); ax.set_xlim(-0.5, 3.8)
    ax.set_title("TF motifs in promoters: % with a site", loc="left", fontsize=9, fontweight="bold"); sm = plt.cm.ScalarMappable(cmap="viridis"); sm.set_clim(0, max(40, max(r[kk] for r in pick for kk in keysg)))
    fig.colorbar(sm, ax=ax, fraction=0.04, pad=0.12, label="% promoters with a site")
    M = np.array([[(tfexpr(r["TF"]) or {}).get(c, np.nan) for c in ("IFNB", "EBOV", "IFNB+EBOV")] for r in pick])
    im = ax2.imshow(M, aspect="auto", cmap="RdBu_r", vmin=-3, vmax=3); ax2.set_xticks(range(3)); ax2.set_xticklabels(["IFNB", "EBOV", "IFNB\n+EBOV"], fontsize=7); ax2.tick_params(axis="y", left=False, labelleft=False)
    ax2.set_title("TF gene log2FC\nvs Mock", fontsize=8); ax.set_xlabel("* FDR < 0.05, + p < 0.05 (not suppressed vs fully suppressed). Sorted: top = more common in NOT-suppressed promoters.", fontsize=6.5, color="grey", loc="left")
    for e in ("png", "pdf"): fig.savefig(f"{OUT}/tf_group_dotplot.{e}", dpi=300, bbox_inches="tight")
    plt.close(fig)
    # ---- volcano: x = log2 odds ratio (not suppressed vs fully suppressed); labels = the largest |log2 OR| among motifs that pass the p / count filter ----
    fig, ax = plt.subplots(figsize=(7.5, 6)); x = np.array([r["log2OR_not_vs_full"] for r in rows]); y = -np.log10([r["p_not_vs_full"] for r in rows])
    colr = ["#1B7837" if r["class"].startswith("NOT") else "#762A83" if r["class"].startswith("fully") else "#9AA0A6" for r in rows]
    ax.scatter(x, y, s=12, c=colr, alpha=0.8, linewidths=0); ax.axvline(0, color="black", lw=0.6); ax.axvline(LOR, color="grey", lw=0.5, ls="--"); ax.axvline(-LOR, color="grey", lw=0.5, ls="--"); ax.axhline(-np.log10(0.05), color="grey", lw=0.5, ls="--")
    def labelled(cls):
        sel = [(abs(r["log2OR_not_vs_full"]), r) for r in rows if r["class"] == cls]; seen_tf = set(); out = []
        for _, r in sorted(sel, key=lambda t: -t[0]):
            if r["TF"] not in seen_tf: seen_tf.add(r["TF"]); out.append(r)
            if len(out) >= N_LAB: break
        return out
    try:
        from adjustText import adjust_text
    except Exception: adjust_text = None
    texts = []
    for cls, colz in (("NOT-suppressed-associated", "#1B7837"), ("fully-suppressed-associated", "#762A83")):
        for r in labelled(cls): texts.append(ax.text(r["log2OR_not_vs_full"], -np.log10(r["p_not_vs_full"]), r["TF"], fontsize=6.5, color=colz))
    if adjust_text and texts: adjust_text(texts, ax=ax, arrowprops=dict(arrowstyle="-", color="grey", lw=0.4))
    ax.set_xlabel("log2 odds ratio  (positive = motif more common in NOT-suppressed promoters;  negative = more common in fully suppressed)", fontsize=7.5); ax.set_ylabel("-log10 p (Fisher exact)"); ax.spines[["top", "right"]].set_visible(False)
    ax.set_title(f"TF motifs: not suppressed ({n_g['not suppressed']} genes) vs fully suppressed ({n_g['fully suppressed']} genes) by EBOV", loc="left", fontsize=9, fontweight="bold")
    ax.text(0.0, -0.13, f"green = NOT-suppressed-associated  (log2OR >= {LOR:g}, p < 0.05, >= {MIN_COUNT} promoters)\npurple = fully-suppressed-associated  (log2OR <= -{LOR:g}, p < 0.05)\nlabels: the {N_LAB} largest |log2 OR| per side", transform=ax.transAxes, ha="left", va="top", fontsize=6.5, color="grey")
    for e in ("png", "pdf"): fig.savefig(f"{OUT}/tf_volcano_not_vs_full.{e}", dpi=300, bbox_inches="tight")
    plt.close(fig)
    cnt = {}
    for r in rows: cnt[r["class"]] = cnt.get(r["class"], 0) + 1
    say("motif classes (not suppressed vs fully suppressed): " + ", ".join(f"{k} {v}" for k, v in sorted(cnt.items())))
    # ---- family level ----
    FAM = [("IRF / STAT", r"IRF|STAT|ISRE"), ("NF-kB / REL", r"NFKB|REL|NF-?KB"), ("AP-1 / CREB / ATF", r"FOS|JUN|ATF|CREB|BATF|AP1|MAF"), ("SP / KLF", r"^SP\d|KLF|^SP1"), ("C/EBP", r"CEBP"), ("ETS", r"^ETS|ELF|ELK|SPI|FLI|ERG|GABP")]
    fam_of = lambda n: next((f for f, rg in FAM if re.search(rg, n, re.I)), "other")
    fams = [f for f, _ in FAM] + ["other"]; fidx = {f: [k for k, m in enumerate(pm) if fam_of(m[1]) == f] for f in fams}
    fig, ax = plt.subplots(figsize=(8.5, 4)); wd = 0.2
    for j, g in enumerate(GROUPS):
        v = [100.0 * (P[grp_arr == g][:, fidx[f]].any(1).mean() if len(fidx[f]) else 0) for f in fams]; ax.bar(np.arange(len(fams)) + j * wd, v, wd, color=GCOL[g], label=f"{g} (n={n_g[g]})")
    ax.set_xticks(np.arange(len(fams)) + 1.5 * wd); ax.set_xticklabels(fams, fontsize=8, rotation=20, ha="right"); ax.set_ylabel("% promoters with >= 1 site of the family"); ax.legend(frameon=False, fontsize=7); ax.spines[["top", "right"]].set_visible(False)
    ax.set_title("TF families in the promoters of IFNB-induced genes, by effect of EBOV", loc="left", fontsize=9, fontweight="bold")
    for e in ("png", "pdf"): fig.savefig(f"{OUT}/tf_family_by_group.{e}", dpi=300, bbox_inches="tight")
    plt.close(fig); say(f"done. Everything is in {OUT}")

if __name__ == "__main__":
    main()
