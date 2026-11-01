#!/usr/bin/env python3
"""Step 1 of the intron / exon check: write a FASTA of INTRON sequences for a subset of genes (your gene + N_CONTROL expressed protein-coding control genes), to be added to the Salmon index (130_intron_salmon.sh).
Why a subset: all mouse introns are ~800 Mb and would make the index too large for a laptop; a few thousand control genes (default 2,000, ~100-150 Mb of introns) are enough to learn how exon and intron signal move together in EBOV vs Mock.
Reads that come from the other genes' introns simply fall on the genome decoys of the index, so they do not distort the counts of the genes that are included.
Per gene: exons of all its transcripts are merged; the gaps inside the gene (introns) are cut out, minus anything that overlaps an exon of ANY other gene (so overlapping genes do not contaminate each other), pieces shorter than MIN_INTRON (100 nt) are dropped;
sequences are written in the gene's own orientation (reverse-complemented for minus-strand genes) so the stranded library (ISR) counts them as sense. Names are <gene_id>_I<k>.
Outputs in <PROJ>/ref/: introns_subset.fa, introns_tx2gene.tsv (name -> gene_id), intron_genes.tsv (the genes used), and <PROJ>/intron_samples.txt (the sample ids of the two groups, with their group).
Env: GENE (Top1), N_CONTROL (2000), MIN_BASEMEAN (10), MIN_INTRON (100), SEED (1), GROUPS (Mock_Mock,Mock_EBOV: the groups whose samples are listed in intron_samples.txt), PROJ.      Usage: python 129_build_intron_reference.py      (needs 108_promoter_tf_groups.py in the same folder for the FASTA reader)
"""
import os, re, sys, csv, random, bisect, importlib.util
from collections import defaultdict

here = os.path.dirname(os.path.abspath(__file__)); p108 = os.path.join(here, "108_promoter_tf_groups.py")
if not os.path.exists(p108): sys.exit("Put 108_promoter_tf_groups.py in the same folder as this script (its FASTA reader is reused)")
spec = importlib.util.spec_from_file_location("p108", p108); M = importlib.util.module_from_spec(spec); spec.loader.exec_module(M)
PROJ = M.PROJ; GTF = M.GTF; GENOME = M.GENOME; RNA = M.RNA; say = M.say
GENE = os.environ.get("GENE", "Top1"); N_CONTROL = int(os.environ.get("N_CONTROL", "2000")); MIN_BASEMEAN = float(os.environ.get("MIN_BASEMEAN", "10")); MIN_INTRON = int(os.environ.get("MIN_INTRON", "100")); SEED = int(os.environ.get("SEED", "1"))
GROUPS = [g.strip() for g in os.environ.get("GROUPS", "Mock_Mock,Mock_EBOV").split(",")]; REFD = os.path.join(PROJ, "ref"); os.makedirs(REFD, exist_ok=True)
def nov(x): return re.sub(r"\.\d+$", "", x)

de = os.path.join(RNA, "DE_v2_EBOV_vs_Mock_annotated.csv"); expr = {}; sym2id = {}
for r in csv.DictReader(open(de)):
    try: bm = float(r["baseMean"])
    except Exception: bm = 0.0
    gid = nov(r["gene_id"]); expr[gid] = (bm, r.get("gene_biotype", ""), r.get("gene_name", ""))
    if r.get("gene_name"): sym2id.setdefault(r["gene_name"].lower(), gid)
if GENE.lower() not in sym2id: sys.exit("%s is not in the DE table" % GENE)
gene_id = sym2id[GENE.lower()]
say("reading exons from the GTF ...")
exons = defaultdict(list); ginfo = {}                      # gene -> [(start,end)], gene -> (chrom, strand)
allex = defaultdict(list)                                  # chrom -> exons of every gene (for subtracting)
with open(GTF) as fh:
    for line in fh:
        if line.startswith("#"): continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 9 or f[2] != "exon": continue
        m = re.search(r'gene_id "([^"]+)"', f[8])
        if not m: continue
        g = nov(m.group(1)); s, e = int(f[3]), int(f[4]); exons[g].append((s, e)); ginfo[g] = (f[0], f[6]); allex[f[0]].append((s, e, g))
def merge(iv):
    out = []
    for s, e in sorted(iv):
        if out and s <= out[-1][1] + 1: out[-1][1] = max(out[-1][1], e)
        else: out.append([s, e])
    return out
merged_all = {c: merge([(s, e) for s, e, g in v]) for c, v in allex.items()}
starts = {c: [s for s, e in v] for c, v in merged_all.items()}
def overlapping_other(chrom, s, e):                        # merged exon blocks (of all genes) that overlap [s,e]
    v = merged_all.get(chrom, []); st = starts.get(chrom, []); i = max(0, bisect.bisect_right(st, e) - 1); out = []
    while i >= 0 and i < len(v):
        a, b = v[i]
        if b < s: break
        if a <= e: out.append((a, b))
        i -= 1
    return sorted(out)
def introns_of(g):
    chrom, strand = ginfo[g]; ex = merge(exons[g]); res = []
    for (a1, b1), (a2, b2) in zip(ex[:-1], ex[1:]):
        s, e = b1 + 1, a2 - 1
        if e - s + 1 < MIN_INTRON: continue
        cur = s
        for a, b in overlapping_other(chrom, s, e):          # remove sequence belonging to exons of other genes (the gene's own exons are outside the gap)
            if a > cur: res.append((cur, a - 1))
            cur = max(cur, b + 1)
        if cur <= e: res.append((cur, e))
    return [(s, e) for s, e in res if e - s + 1 >= MIN_INTRON]
cands = [g for g, (bm, bt, nm) in expr.items() if bt == "protein_coding" and bm >= MIN_BASEMEAN and g in exons and g != gene_id and len(merge(exons[g])) > 1]
random.seed(SEED); ctrl = random.sample(cands, min(N_CONTROL, len(cands))); use = [gene_id] + ctrl; say("%s + %d control genes" % (GENE, len(ctrl)))
fa = M.Fasta(GENOME); n_int = n_nt = 0
with open(os.path.join(REFD, "introns_subset.fa"), "w") as fo, open(os.path.join(REFD, "introns_tx2gene.tsv"), "w") as mp, open(os.path.join(REFD, "intron_genes.tsv"), "w") as gl:
    mp.write("TXNAME\tGENEID\n"); gl.write("gene_id\tgene_name\trole\tn_introns\tintron_nt\n")
    for g in use:
        if g not in ginfo or ginfo[g][0] not in fa.idx: continue
        chrom, strand = ginfo[g]; iv = introns_of(g); k = 0; tot = 0
        for s, e in iv:
            seq = fa.get(chrom, s, e)
            if strand == "-": seq = M.revcomp(seq)
            if len(seq) < MIN_INTRON or seq.count("N") > 0.2 * len(seq): continue
            k += 1; tot += len(seq); name = "%s_I%d" % (g, k); fo.write(">%s\n%s\n" % (name, "\n".join(seq[i:i + 80] for i in range(0, len(seq), 80)))); mp.write("%s\t%s\n" % (name, g))
        n_int += k; n_nt += tot; gl.write("%s\t%s\t%s\t%d\t%d\n" % (g, expr.get(g, (0, "", g))[2], "query" if g == gene_id else "control", k, tot))
say("wrote %d introns (%.1f Mb) to %s" % (n_int, n_nt / 1e6, os.path.join(REFD, "introns_subset.fa")))
si = os.path.join(PROJ, "sample_info.csv"); first = open(si).readline(); sep = "\t" if "\t" in first else ","
rows = list(csv.DictReader(open(si), delimiter=sep)); out = []
for r in rows:
    c1 = r.get("condition1") or "None"; c2 = r.get("condition2") or "None"; ifn = "Mock" if c1 in ("None", "none", "") else c1; vir = "Mock" if c2 in ("None", "none", "") else c2; grp = "%s_%s" % (ifn, vir)
    if grp in GROUPS: out.append((r["sample_id"], grp))
open(os.path.join(PROJ, "intron_samples.txt"), "w").write("".join("%s\t%s\n" % x for x in out)); say("samples for the intron quantification: " + ", ".join("%s(%s)" % x for x in out))
