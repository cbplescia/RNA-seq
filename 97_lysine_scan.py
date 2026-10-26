#!/usr/bin/env python3
"""Lysine inventory and rule-based flags for a protein FASTA (standard library only).

This is NOT a trained predictor. Peptidyl-lysine modifications (ubiquitination, acetylation, methylation, SUMOylation, hydroxylation ...) are
predicted properly by trained models such as MusiteDeep (web server musite.net or its offline version), GPS-Uber / GPS-SUMO or similar - see the notes printed at the end. What this script does:
  * lists every lysine with its position and a +/- WINDOW residue context (the 15-mer / 21-mer many web predictors want),
  * flags a few sequence rules that are known but weak:
      SUMO_consensus   psi-K-x-[DE]   (psi = I, L, V, M, F)  - the classic SUMOylation motif; real sites often do not match it, and matches are not always modified
      histone_like     [AG]-K-[ST]    - the A/G-K-S/T context of histone H3 / H4 lysine methylation and acetylation sites (weak outside histones)
      K_cluster        >= 3 K/R within +/- 3 residues - a basic patch, where acetylation / methylation / ubiquitination are often reported
  * writes eGP_lysines.tsv, a windows FASTA (one 2*WINDOW+1-mer per lysine) for pasting into predictors, and a per-record summary.
Usage:  python 97_lysine_scan.py eGP.FASTA   (default: the first file in the folder named egp.fasta / eGP.fasta, any case; WINDOW=10 python ... for a different window)
"""
import os, sys, glob, re

WINDOW = int(os.environ.get("WINDOW", "10"))


def read_fasta(path):
    recs, name, seq = [], None, []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    recs.append((name, "".join(seq).upper()))
                name, seq = line[1:].strip(), []
            else:
                seq.append(re.sub(r"[^A-Za-z*]", "", line))
        if name is not None:
            recs.append((name, "".join(seq).upper()))
    return [(n, s.replace("*", "")) for n, s in recs]


def find_default():
    for p in glob.glob("*"):
        if os.path.basename(p).lower() in ("egp.fasta", "egp.fa", "egp.faa", "egp.fasta.txt"):
            return p
    sys.exit("No FASTA given and no eGP.FASTA found in this folder. Usage: python 97_lysine_scan.py eGP.FASTA")


def flags(seq, i):
    f = []
    if i + 2 < len(seq) and i >= 1 and seq[i - 1] in "ILVMF" and seq[i + 2] in "DE":
        f.append("SUMO_consensus")
    if 1 <= i < len(seq) - 1 and seq[i - 1] in "AG" and seq[i + 1] in "ST":
        f.append("histone_like")
    lo, hi = max(0, i - 3), min(len(seq), i + 4)
    if sum(1 for c in seq[lo:hi] if c in "KR") >= 3:
        f.append("K_cluster")
    return f


def main():
    path = sys.argv[1] if len(sys.argv) > 1 else find_default()
    recs = read_fasta(path)
    if not recs:
        sys.exit(f"No sequences found in {path}")
    base = os.path.splitext(os.path.basename(path))[0]
    tsv, fa = f"{base}_lysines.tsv", f"{base}_lysine_windows.fasta"
    with open(tsv, "w") as t, open(fa, "w") as w:
        t.write("record\tposition\twindow\tflags\n")
        for name, seq in recs:
            ks = [i for i, c in enumerate(seq) if c == "K"]
            nflag = {"SUMO_consensus": 0, "histone_like": 0, "K_cluster": 0}
            for i in ks:
                a, b = i - WINDOW, i + WINDOW + 1
                win = ("-" * max(0, -a)) + seq[max(0, a):b] + ("-" * max(0, b - len(seq)))
                fl = flags(seq, i)
                for x in fl:
                    nflag[x] += 1
                t.write(f"{name.split()[0]}\t{i + 1}\t{win}\t{','.join(fl)}\n")
                w.write(f">{name.split()[0]}_K{i + 1}\n{win}\n")
            print(f"\n{name}\n  length {len(seq)} aa; lysines {len(ks)} ({100 * len(ks) / max(1, len(seq)):.1f}%)")
            print("  flags: " + ", ".join(f"{k} {v}" for k, v in nflag.items()))
            for i in ks:
                fl = flags(seq, i)
                if fl:
                    print(f"    K{i + 1:<5d} {seq[max(0, i - 4):i]}[K]{seq[i + 1:i + 5]}  {','.join(fl)}")
    print(f"\nwrote {tsv} and {fa}")
    print("""
Next step - an actual predictor (these are trained models; check each tool's current input rules and licence):
  * MusiteDeep: web server https://www.musite.net/ (upload the FASTA, choose the lysine PTM types such as ubiquitination, SUMOylation, acetyllysine, methyllysine, hydroxylysine);
    an offline version is also distributed by the authors for large jobs.
  * GPS-SUMO / GPS-Uber and related tools from the GPS suite (gps.biocuckoo.cn) for SUMOylation / ubiquitination.
  * Predictions are probabilities, not observations: confirm important sites by mass spectrometry or mutagenesis (K to R).""")


if __name__ == "__main__":
    main()
