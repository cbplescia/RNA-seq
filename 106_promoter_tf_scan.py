#!/usr/bin/env python3
"""Find the promoters of selected mouse genes in YOUR reference genome and scan them for transcription-factor binding sites.

Genes (default): Ifitm3, Srms, Cert1, Axl.   Reference: the Ensembl GRCm39 files your RNA-seq was mapped to (<project>/ref/Mus_musculus.GRCm39.109.gtf and
Mus_musculus.GRCm39.dna.primary_assembly.fa; override with GTF=... and GENOME=...).

What it does
  1. reads the gene's transcripts from the GTF, finds the transcription start sites (TSS; clusters within 200 bp are merged, the canonical transcript's cluster is marked) and cuts the promoter
     region from the genome FASTA: UP bp upstream and DOWN bp downstream of each TSS (default 2000 / 500), reverse-complemented for genes on the minus strand so every sequence reads 5' -> 3' along the gene
  2. scans each promoter on both strands with JASPAR vertebrate position-weight matrices (the file is downloaded if the URL below is reachable, or give JASPAR_FILE=your.txt; if neither works, a small
     built-in panel of IFN-relevant motifs is used so you still get a result). A site is reported when its relative PWM score >= REL (default 0.88; 1.0 = the best possible match)
  3. joins the transcription factors with your RNA-seq (DE_v2_IFNB_vs_Mock / EBOV_vs_Mock annotated tables): is that TF's gene expressed in your cells, and is it changed by IFNB or EBOV?
IMPORTANT: a motif match is a POSSIBLE binding site, not evidence of binding - short motifs match by chance (hundreds of hits per kb at low thresholds) and many TFs share similar motifs. Real
binding / enhancer evidence needs ChIP-seq or ATAC / H3K27ac data from mouse macrophages (e.g. ChIP-Atlas, ENCODE, GEO). Public mouse tracks are mostly mm10 (GRCm38); your reference is GRCm39,
so coordinates would need a liftover. Distal enhancers cannot be found from sequence alone: widen the window (UP=10000) only as a rough guide.
Outputs in <project>/requant_v2_outputs/promoter_scan/:
  tss_table.tsv, promoters.bed (GRCm39), promoters.fasta, motif_hits.tsv (every site: TF, matrix, score, position relative to the TSS, genomic position, sequence),
  tf_summary.tsv (per gene and TF: sites, best score), tf_expression.tsv (TFs with sites + your expression data), tf_hits_heatmap.png / .pdf, promoter_map.png / .pdf
Usage:  python 106_promoter_tf_scan.py     (options as environment variables: GENES=Ifitm3,Srms,Cert1,Axl  UP=2000  DOWN=500  REL=0.88  JASPAR_FILE=path  PROJ=path)
Needs: python 3, numpy, matplotlib   (pip install numpy matplotlib)
"""
import os, re, sys, csv, gzip, urllib.request
import numpy as np

PROJ = os.environ.get("PROJ", os.path.join(os.path.expanduser("~"), "Documents/aviti_rnaseq/2025/aviti/20251009-100725_Maury_25161_WivfFacICOTHFAUJwSPeHBEFfCfNtCRpUSgVGUkz/results/Maury_25161"))
GTF = os.environ.get("GTF", os.path.join(PROJ, "ref", "Mus_musculus.GRCm39.109.gtf"))
GENOME = os.environ.get("GENOME", os.path.join(PROJ, "ref", "Mus_musculus.GRCm39.dna.primary_assembly.fa"))
GENES = [g.strip() for g in os.environ.get("GENES", "Ifitm3,Srms,Cert1,Axl").split(",") if g.strip()]
UP, DOWN, REL = int(os.environ.get("UP", "2000")), int(os.environ.get("DOWN", "500")), float(os.environ.get("REL", "0.88"))
OUT = os.path.join(PROJ, "requant_v2_outputs", "promoter_scan"); os.makedirs(OUT, exist_ok=True)
RNA = os.path.join(PROJ, "requant_v2_outputs")
JASPAR_URLS = ["https://jaspar.elixir.no/download/data/2024/CORE/JASPAR2024_CORE_vertebrates_non-redundant_pfms_jaspar.txt",
               "https://jaspar.genereg.net/download/data/2024/CORE/JASPAR2024_CORE_vertebrates_non-redundant_pfms_jaspar.txt",
               "https://jaspar2022.genereg.net/download/data/2022/CORE/JASPAR2022_CORE_vertebrates_non-redundant_pfms_jaspar.txt"]
COMP = str.maketrans("ACGTacgtN", "TGCAtgcaN")
def revcomp(s): return s.translate(COMP)[::-1]
def say(*a): print(*a, flush=True)

# ------------------------------------------------------------------ genome ----
def build_fai(path):
    fai = path + ".fai"
    if os.path.exists(fai): return fai
    say(f"building the FASTA index {fai} (one pass over the genome, about a minute) ...")
    rows, name, off, length, lw, lb = [], None, 0, 0, None, None
    with open(path, "rb") as fh:
        pos = 0
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
    def get(self, chrom, start, end):                        # 1-based inclusive
        ln, st, lw, lb = self.idx[chrom]; start, end = max(1, start), min(ln, end)
        a, b = start - 1, end
        off = st + (a // lw) * lb + a % lw; nbytes = (b - a) + ((b - 1) // lw - a // lw) * (lb - lw) + (lb - lw)
        self.fh.seek(off); return re.sub(r"\s", "", self.fh.read(nbytes).decode())[: b - a].upper()

# --------------------------------------------------------------------- GTF ----
def tss_for(gene):
    tx = []
    op = gzip.open if GTF.endswith(".gz") else open
    with op(GTF, "rt") as fh:
        for line in fh:
            if line.startswith("#"): continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "transcript": continue
            a = f[8]
            m = re.search(r'gene_name "([^"]+)"', a)
            if not m or m.group(1).lower() != gene.lower(): continue
            tid = re.search(r'transcript_id "([^"]+)"', a).group(1); gid = re.search(r'gene_id "([^"]+)"', a).group(1)
            bt = re.search(r'transcript_biotype "([^"]+)"', a); canon = 'tag "Ensembl_canonical"' in a
            tss = int(f[3]) if f[6] == "+" else int(f[4])
            tx.append(dict(gene=m.group(1), gene_id=gid, tid=tid, chrom=f[0], strand=f[6], tss=tss, biotype=bt.group(1) if bt else "", canonical=canon))
    return tx

def promoters(gene):
    tx = tss_for(gene)
    if not tx: return []
    pc = [t for t in tx if t["biotype"] in ("protein_coding", "")] or tx
    cl = []
    for t in sorted(pc, key=lambda x: x["tss"]):
        for c in cl:
            if abs(c["tss"] - t["tss"]) <= 200: c["members"].append(t); break
        else: cl.append(dict(tss=t["tss"], members=[t]))
    for c in cl:
        c["n"] = len(c["members"]); c["canonical"] = any(m["canonical"] for m in c["members"]); c["tss"] = int(np.median([m["tss"] for m in c["members"]]))
    cl.sort(key=lambda c: (-int(c["canonical"]), -c["n"])); out = []
    for k, c in enumerate(cl[:3]):
        t0 = c["members"][0]
        out.append(dict(gene=t0["gene"], gene_id=t0["gene_id"], chrom=t0["chrom"], strand=t0["strand"], tss=c["tss"], n_tx=c["n"], canonical=c["canonical"], name=f"P{k + 1}", tids=",".join(m["tid"] for m in c["members"][:4])))
    return out

# ------------------------------------------------------------------- motifs ----
def parse_jaspar(text):
    ms, cur = [], None
    for line in text.splitlines():
        line = line.strip()
        if line.startswith(">"):
            p = line[1:].split(None, 1); cur = dict(id=p[0], name=p[1] if len(p) > 1 else p[0], rows=[]); ms.append(cur)
        elif cur and line and line[0] in "ACGT":
            cur["rows"].append([float(x) for x in re.findall(r"[-\d.]+", line.split("[", 1)[1])])
    out = []
    for m in ms:
        if len(m["rows"]) == 4:
            c = np.array(m["rows"]).T; out.append((m["id"], m["name"], c))        # width x 4 counts (A,C,G,T)
    return out

IUPAC = {"A": "A", "C": "C", "G": "G", "T": "T", "R": "AG", "Y": "CT", "S": "CG", "W": "AT", "K": "GT", "M": "AC", "N": "ACGT"}
BUILTIN = {"ISRE_IRF9/STAT1/STAT2": "GAAANNGAAACT", "IRF_ISRE-like": "GAAANNGAAA", "GAS_STAT1": "TTCNNNGAA", "NFKB": "GGGRNNYYCC", "SP1": "GGGCGG", "AP1": "TGASTCA", "CREB/ATF": "TGACGTCA", "CEBP": "TTGCGCAA", "ETS": "GGAAGT", "GATA": "WGATAR"}
def builtin_motifs():
    out = []
    for name, pat in BUILTIN.items():
        c = np.full((len(pat), 4), 0.0)
        for i, ch in enumerate(pat):
            for b in IUPAC[ch]: c[i, "ACGT".index(b)] = 20.0 / len(IUPAC[ch])
        out.append((f"builtin:{name}", name, c + 0.1))
    return out

def load_motifs():
    jf = os.environ.get("JASPAR_FILE")
    if jf and os.path.exists(jf): return parse_jaspar(open(jf).read()), f"JASPAR file {jf}"
    cache = os.path.join(OUT, "jaspar_pfms.txt")
    if os.path.exists(cache) and os.path.getsize(cache) > 100000: return parse_jaspar(open(cache).read()), f"JASPAR (cached {cache})"
    for u in JASPAR_URLS:
        try:
            txt = urllib.request.urlopen(u, timeout=60).read().decode(); ms = parse_jaspar(txt)
            if len(ms) > 100: open(cache, "w").write(txt); return ms, f"JASPAR download ({u.split('/')[-1]})"
        except Exception as e: say(f"  could not get {u.split('/')[2]}: {type(e).__name__}")
    say("JASPAR could not be downloaded - using the small built-in motif panel (download the vertebrate non-redundant PFMs in JASPAR format from jaspar.elixir.no and rerun with JASPAR_FILE=...)")
    return builtin_motifs(), "built-in IFN-relevant panel"

def pwm(counts, pc=0.8):
    n = counts.sum(1, keepdims=True); p = (counts + pc * 0.25) / (n + pc); return np.log2(p / 0.25)

def scan(seq, motifs, rel):
    enc = np.array([("ACGT".find(c)) for c in seq]); L = len(seq); hits = []
    for mid, name, counts in motifs:
        w = counts.shape[0]
        if L < w: continue
        W = pwm(counts); mn, mx = W.min(1).sum(), W.max(1).sum()
        if mx <= mn: continue
        for strand, e in (("+", enc), ("-", np.where(enc >= 0, 3 - enc, -1)[::-1])):
            s = np.zeros(L - w + 1); bad = np.zeros(L - w + 1, bool)
            for j in range(w):
                col = e[j:L - w + 1 + j]; bad |= col < 0; s += W[j, np.where(col < 0, 0, col)]
            r = (s - mn) / (mx - mn); r[bad] = -1
            for i in np.flatnonzero(r >= rel):
                start = i if strand == "+" else L - (i + w)                       # position on the given (5'->3' of the gene) sequence
                hits.append((mid, name, float(r[i]), strand, int(start), int(start + w - 1)))
    return hits

def tfgenes(name):                                                              # JASPAR name -> candidate mouse gene symbols
    out = []
    for tok in re.split(r"::", re.sub(r"\(var\.\d+\)", "", name)): out.append(tok.strip().lower())
    return out

# --------------------------------------------------------------------- main ----
def main():
    for p in (GTF, GENOME):
        if not os.path.exists(p): sys.exit(f"Missing {p} - set GTF=... / GENOME=... or PROJ=...")
    fa = Fasta(GENOME); motifs, src = load_motifs(); say(f"{len(motifs)} motifs from {src}; relative score >= {REL}; promoter = {UP} bp upstream .. {DOWN} bp downstream of each TSS")
    proms = []
    for g in GENES:
        ps = promoters(g)
        for alt in {"cert1": ["Col4a3bp"], "gba1": ["Gba"]}.get(g.lower(), []):          # older Ensembl symbols
            if not ps: ps = promoters(alt)
        if not ps: say(f"  {g}: not found in the GTF (check the symbol)"); continue
        say(f"  {g}: {len(ps)} promoter(s) - " + "; ".join(f"{p['name']} chr{p['chrom']}:{p['tss']}({p['strand']}) {p['n_tx']} transcript(s){' canonical' if p['canonical'] else ''}" for p in ps)); proms += ps
    if not proms: sys.exit("No promoters found")
    hits = []
    with open(f"{OUT}/tss_table.tsv", "w") as t, open(f"{OUT}/promoters.bed", "w") as b, open(f"{OUT}/promoters.fasta", "w") as fo:
        t.write("gene\tpromoter\tchrom\tstrand\tTSS\ttranscripts_at_TSS\tcanonical\ttranscript_ids\n")
        for p in proms:
            if p["strand"] == "+": gs, ge = p["tss"] - UP, p["tss"] + DOWN - 1
            else: gs, ge = p["tss"] - DOWN + 1, p["tss"] + UP
            seq = fa.get(p["chrom"], gs, ge); seq = seq if p["strand"] == "+" else revcomp(seq); p.update(gs=gs, ge=ge, seq=seq)
            t.write(f"{p['gene']}\t{p['name']}\t{p['chrom']}\t{p['strand']}\t{p['tss']}\t{p['n_tx']}\t{p['canonical']}\t{p['tids']}\n")
            b.write(f"{p['chrom']}\t{gs - 1}\t{ge}\t{p['gene']}_{p['name']}\t0\t{p['strand']}\n"); fo.write(f">{p['gene']}_{p['name']} chr{p['chrom']}:{gs}-{ge}({p['strand']}) TSS at index {UP}\n{seq}\n")
            for mid, name, sc, st, a, z in scan(seq, motifs, REL):
                rel = lambda i: i - UP if i < UP else i - UP + 1
                gpos = (gs + a, gs + z) if p["strand"] == "+" else (ge - z, ge - a)
                hits.append(dict(gene=p["gene"], promoter=p["name"], chrom=p["chrom"], motif_id=mid, TF=name, score=sc, site_strand=st, rel_start=rel(a), rel_end=rel(z), g_start=gpos[0], g_end=gpos[1], seq=seq[a:z + 1]))
    cols = ["gene", "promoter", "chrom", "motif_id", "TF", "score", "site_strand", "rel_start", "rel_end", "g_start", "g_end", "seq"]
    with open(f"{OUT}/motif_hits.tsv", "w") as o:
        o.write("\t".join(cols) + "\n")
        for h in sorted(hits, key=lambda h: (h["gene"], h["promoter"], h["rel_start"])): o.write("\t".join(f"{h[c]:.3f}" if c == "score" else str(h[c]) for c in cols) + "\n")
    say(f"{len(hits)} motif sites in {len(proms)} promoter(s); table written")
    # --- summary per gene x TF
    S = {}
    for h in hits:
        k = (h["gene"], h["TF"]); d = S.setdefault(k, dict(n=0, best=0.0, pos=[])); d["n"] += 1; d["best"] = max(d["best"], h["score"]); d["pos"].append(h["rel_start"])
    # --- RNA-seq expression of the TFs
    expr = {}
    for nm, col in (("IFNB_vs_Mock", "IFNB"), ("EBOV_vs_Mock", "EBOV")):
        f = os.path.join(RNA, f"DE_v2_{nm}_annotated.csv")
        if os.path.exists(f):
            for r in csv.DictReader(open(f)):
                e = expr.setdefault(r["gene_name"].lower(), {}); e["baseMean"] = r.get("baseMean"); e[col] = (r.get("log2FoldChange"), r.get("padj"))
    def tfexpr(name):
        for g in tfgenes(name):
            if g in expr: return g, expr[g]
        return None, None
    rows = []
    for (gene, tf), d in sorted(S.items(), key=lambda kv: (kv[0][0], -kv[1]["best"])):
        g, e = tfexpr(tf); f = lambda c, i: (e.get(c) or ("", ""))[i] if e else ""
        rows.append([gene, tf, d["n"], f"{d['best']:.3f}", ",".join(map(str, sorted(d["pos"])[:12])), g or "", e["baseMean"] if e else "", f("IFNB", 0), f("IFNB", 1), f("EBOV", 0), f("EBOV", 1)])
    with open(f"{OUT}/tf_summary.tsv", "w") as o:
        o.write("gene\tTF_motif\tn_sites\tbest_score\tsite_positions_rel_TSS\tmouse_gene_in_your_data\tbaseMean\tlog2FC_IFNB\tpadj_IFNB\tlog2FC_EBOV\tpadj_EBOV\n")
        for r in rows: o.write("\t".join(map(str, r)) + "\n")
    with open(f"{OUT}/tf_expression.tsv", "w") as o:
        o.write("TF_motif\tgenes_with_sites\tmouse_gene_in_your_data\tbaseMean\tlog2FC_IFNB\tpadj_IFNB\tlog2FC_EBOV\tpadj_EBOV\n")
        seen = {}
        for r in rows: seen.setdefault(r[1], dict(genes=[], r=r))["genes"].append(r[0])
        for tf, d in sorted(seen.items(), key=lambda kv: -len(kv[1]["genes"])): o.write("\t".join(map(str, [tf, ",".join(sorted(set(d["genes"]))), d["r"][5], d["r"][6], d["r"][7], d["r"][8], d["r"][9], d["r"][10]])) + "\n")
    for gene in GENES:
        top = [r for r in rows if r[0] == gene][:8]
        if top: say(f"  {gene}: top motifs " + ", ".join(f"{r[1]} ({r[2]} sites, {r[3]})" for r in top))
    # --- figures
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    fam = lambda n: ("IRF/STAT (ISRE, GAS)" if re.search(r"IRF|STAT|ISRE|GAS", n, re.I) else "NF-kB / REL" if re.search(r"NFKB|REL|NF-?KB", n, re.I) else "AP-1 / CREB / ATF" if re.search(r"FOS|JUN|ATF|CREB|BATF|AP1|MAF", n, re.I) else "SP / KLF" if re.search(r"^SP\d|KLF|SP1", n, re.I) else "CEBP" if re.search(r"CEBP", n, re.I) else "ETS" if re.search(r"^ETS|ELF|ELK|SPI|FLI|ERG|GABP", n, re.I) else "other")
    FC = {"IRF/STAT (ISRE, GAS)": "#C53030", "NF-kB / REL": "#2B6CB0", "AP-1 / CREB / ATF": "#E08214", "SP / KLF": "#1B7837", "CEBP": "#762A83", "ETS": "#8C510A", "other": "#9AA0A6"}
    fig, axs = plt.subplots(len(proms), 1, figsize=(12, 1.6 * len(proms) + 1.2), sharex=True); axs = np.atleast_1d(axs)
    byp = {}
    for h in hits: byp.setdefault((h["gene"], h["promoter"]), []).append(h)
    for ax, p in zip(axs, proms):
        for h in byp.get((p["gene"], p["name"]), []): ax.vlines((h["rel_start"] + h["rel_end"]) / 2, 0, 1, color=FC[fam(h["TF"])], lw=0.8, alpha=0.8)
        ax.axvline(0, color="black", lw=1.2); ax.set_yticks([]); ax.set_ylabel(f"{p['gene']}\n{p['name']}", rotation=0, ha="right", va="center", fontsize=8); ax.set_xlim(-UP, DOWN)
        for s in ("top", "right", "left"): ax.spines[s].set_visible(False)
    axs[-1].set_xlabel("position relative to the TSS (bp); black line = TSS"); axs[0].legend(handles=[plt.Line2D([0], [0], color=c, lw=2, label=k) for k, c in FC.items()], ncol=4, fontsize=7, frameon=False, loc="upper left", bbox_to_anchor=(0, 1.9))
    for e in ("png", "pdf"): fig.savefig(f"{OUT}/promoter_map.{e}", dpi=300, bbox_inches="tight")
    plt.close(fig)
    best = {}
    for r in rows: best.setdefault(r[1], {})[r[0]] = float(r[3])
    ordtf = sorted(best, key=lambda t: (-len(best[t]), -max(best[t].values())))[:60]
    if ordtf:
        M = np.array([[best[t].get(g, np.nan) for g in GENES] for t in ordtf])
        fig, ax = plt.subplots(figsize=(0.9 * len(GENES) + 4, 0.2 * len(ordtf) + 2))
        im = ax.imshow(M, aspect="auto", cmap="viridis", vmin=REL, vmax=1); ax.set_xticks(range(len(GENES))); ax.set_xticklabels(GENES, rotation=30, ha="right"); ax.set_yticks(range(len(ordtf))); ax.set_yticklabels(ordtf, fontsize=6.5)
        fig.colorbar(im, ax=ax, fraction=0.04, label="best relative score"); ax.set_title(f"TF motifs with a site in the promoters (top {len(ordtf)}, shared first)", loc="left", fontsize=9, fontweight="bold")
        for e in ("png", "pdf"): fig.savefig(f"{OUT}/tf_hits_heatmap.{e}", dpi=300, bbox_inches="tight")
        plt.close(fig)
    say(f"done. Everything is in {OUT}")

if __name__ == "__main__":
    main()
