#!/usr/bin/env python3
"""Which OTHER genes share sequence with GENE (default Top1) at the level Salmon sees - 31-mers - so reads could be split between them or mis-assigned?
Takes every 31-mer (forward strand, as in the transcripts) of the gene's transcripts in ref/Mus_musculus.GRCm39.transcripts.fa, then streams all other transcripts and counts, per other gene, how many of those 31-mers they also contain.
A gene that shares a sizeable fraction of the 31-mers (processed pseudogene / retrocopy, paralog, overlapping transcript) is a candidate source of mis-assigned reads. Sharing with other isoforms of the SAME gene is ignored.
Not covered: unannotated retrocopies (they are not in the transcript FASTA; your decoy-aware index uses the genome as decoys, which absorbs most such reads) - for those, BLAT / BLAST the gene's CDS against the GRCm39 genome.
Outputs (requant_v2_outputs/figures_<gene>_mechanism/): <gene>_kmer_sharing.tsv (every other gene that shares >= 1 31-mer: shared 31-mers, % of the gene's 31-mers, biotype, transcript).
Env: GENE (default Top1), K (31), TX_FA, PROJ.   Usage: python 127_kmer_paralog_check.py      (python 3 + numpy; a few minutes)
"""
import os, re, sys
import numpy as np

PROJ = os.environ.get("PROJ", os.path.join(os.path.expanduser("~"), "Documents/aviti_rnaseq/2025/aviti/20251009-100725_Maury_25161_WivfFacICOTHFAUJwSPeHBEFfCfNtCRpUSgVGUkz/results/Maury_25161"))
GENE = os.environ.get("GENE", "Top1"); K = int(os.environ.get("K", "31")); TX_FA = os.environ.get("TX_FA", os.path.join(PROJ, "ref", "Mus_musculus.GRCm39.transcripts.fa"))
OUT = os.path.join(PROJ, "requant_v2_outputs", "figures_%s_mechanism" % GENE.lower()); os.makedirs(OUT, exist_ok=True)
if not os.path.exists(TX_FA): sys.exit("Missing %s - set TX_FA=..." % TX_FA)
LUT = np.full(256, 4, dtype=np.int64)
for i, b in enumerate("ACGT"): LUT[ord(b)] = i; LUT[ord(b.lower())] = i
def kmers(seq):
    a = LUT[np.frombuffer(seq.encode(), dtype=np.uint8)]; n = len(a) - K + 1
    if n <= 0: return np.empty(0, dtype=np.int64)
    val = np.zeros(n, dtype=np.int64); bad = np.zeros(n, dtype=bool)
    for j in range(K): c = a[j:j + n]; bad |= c > 3; val = (val << 2) | np.minimum(c, 3)
    return val[~bad]
def records():
    hdr, seq = None, []
    with open(TX_FA) as fh:
        for line in fh:
            if line.startswith(">"):
                if hdr is not None: yield hdr, "".join(seq)
                hdr, seq = line[1:].strip(), []
            else: seq.append(line.strip())
    if hdr is not None: yield hdr, "".join(seq)
def meta(h):
    g = re.search(r"gene_symbol:(\S+)", h); gi = re.search(r"gene:(ENSMUSG\d+)", h); bt = re.search(r"gene_biotype:(\S+)", h)
    return (g.group(1) if g else (gi.group(1) if gi else "?")), (gi.group(1) if gi else "?"), (bt.group(1) if bt else "")
mine, own_ids = [], set()
for h, s in records():
    sym, gid, bt = meta(h)
    if sym.lower() == GENE.lower(): mine.append(kmers(s)); own_ids.add(gid)
if not mine: sys.exit("No transcripts with gene_symbol:%s found in %s (check the FASTA header format)" % (GENE, TX_FA))
top = np.unique(np.concatenate(mine)); print("%s: %d transcripts, %d distinct %d-mers" % (GENE, len(mine), len(top), K), flush=True)
hits = {}; n = 0
for h, s in records():
    sym, gid, bt = meta(h)
    if gid in own_ids or sym.lower() == GENE.lower(): continue
    n += 1
    km = kmers(s)
    if not len(km): continue
    sh = np.unique(km[np.isin(km, top, assume_unique=False)])
    if len(sh):
        d = hits.setdefault((sym, gid, bt), set()); d.update(sh.tolist())
    if n % 20000 == 0: print("  %d transcripts checked" % n, flush=True)
rows = sorted(((len(v), k) for k, v in hits.items()), reverse=True)
with open(os.path.join(OUT, "%s_kmer_sharing.tsv" % GENE.lower()), "w") as o:
    o.write("gene\tgene_id\tbiotype\tshared_%dmers\tpct_of_%s_%dmers\n" % (K, GENE, K))
    for c, (sym, gid, bt) in rows: o.write("%s\t%s\t%s\t%d\t%.2f\n" % (sym, gid, bt, c, 100.0 * c / len(top)))
print("\n%d other genes share at least one %d-mer with %s. Top 20:" % (len(rows), K, GENE))
for c, (sym, gid, bt) in rows[:20]: print("  %-14s %-18s %-26s %6d shared  (%.2f %% of %s's %d-mers)" % (sym, gid, bt, c, 100.0 * c / len(top), GENE, K))
print("\nA gene sharing more than a few percent (especially a processed pseudogene / Gm gene) is worth a closer look: reads from it can be counted as %s, or %s reads moved to it." % (GENE, GENE))
print("Wrote %s" % os.path.join(OUT, "%s_kmer_sharing.tsv" % GENE.lower()))
