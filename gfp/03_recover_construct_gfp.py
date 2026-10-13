#!/usr/bin/env python3
"""Reconstruct your construct's GFP nucleotide sequence from the reads.

The screen found reads that encode GFP protein but share no exact nucleotide
k-mer with the reference, which means the construct carries a different codon
variant. Those reads *are* the construct's sequence. This anchors each one to
the reference through its translation, piles them up, and calls a consensus.

The result is the reporter's nucleotide sequence as it exists in your virus,
derived from your own data rather than from a plasmid map.

  python3 -I 03_recover_construct_gfp.py ref/gfp.fa gfp_screen_aa/hits_*.txt

Options (environment variables):
  MIN_DEPTH   minimum reads to call a position, default 3
  PEPK        peptide anchor length, must match the screen's, default 9
  OUT         output FASTA path, default ref/construct_gfp_recovered.fa

Assumes no insertions or deletions relative to the reference protein. Each read
is placed independently, so a local indel costs that read rather than the run,
but a construct carrying a genuine indel will show as a run of disagreement -
which the per-position report will make visible.
"""

import os
import sys
from collections import Counter, defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from gfp_tools import (CODONS, find_chromophore, longest_orf,  # noqa: E402
                       read_fasta, revcomp, translate)

MIN_DEPTH = int(os.environ.get("MIN_DEPTH", 3))
PEPK = int(os.environ.get("PEPK", 9))
OUT = os.environ.get("OUT", "ref/construct_gfp_recovered.fa")
# Reads anchored inside the ORF often extend past the reference CDS. Piling up
# those overhangs recovers the construct's own start codon (which a reference
# with a different N-terminus places outside its coordinates) and whatever
# sequence flanks the cassette.
FLANK = int(os.environ.get("FLANK", 150))


def ref_cds(path):
    recs = read_fasta(path)
    if not recs:
        sys.exit(f"no records in {path}")
    name, seq = recs[0]
    if translate(seq)[:-1].count("*"):
        found = longest_orf(seq)
        if found:
            seq = found[0]
    return name, seq


def placements(read, aa_pos_of, pepk):
    """Offsets into the reference CDS implied by this read's peptide anchors.

    Returns (oriented_read, offset) or None. Several anchors in one read must
    agree, which rejects a read placed by a peptide that recurs in the protein.
    """
    for strand in (read, revcomp(read)):
        votes = Counter()
        for frame in range(3):
            s = strand[frame:]
            n = (len(s) // 3) * 3
            if n < pepk * 3:
                continue
            aa = "".join(CODONS.get(s[i:i + 3], "X") for i in range(0, n, 3))
            for i in range(len(aa) - pepk + 1):
                p = aa_pos_of.get(aa[i:i + pepk])
                if p is not None:
                    # Reference nt position that this read's position 0 sits at.
                    votes[p * 3 - (frame + i * 3)] += 1
        if votes:
            best, count = votes.most_common(1)[0]
            # An unambiguous placement: no competing offset with equal support.
            if sum(1 for v in votes.values() if v == count) == 1:
                return strand, best
    return None


def main():
    if len(sys.argv) < 3:
        sys.exit(__doc__)

    name, ref = ref_cds(sys.argv[1])
    ref_aa = translate(ref).rstrip("*")
    aa_pos_of = {ref_aa[i:i + PEPK]: i for i in range(len(ref_aa) - PEPK + 1)}

    pile = defaultdict(Counter)
    n_reads = n_placed = n_ambiguous = 0
    per_file = {}

    for path in sys.argv[2:]:
        try:
            reads = [l.strip().upper() for l in open(path) if l.strip()]
        except OSError as e:
            print(f"  {path}: {e}", file=sys.stderr)
            continue
        placed_here = 0
        for r in reads:
            n_reads += 1
            got = placements(r, aa_pos_of, PEPK)
            if got is None:
                n_ambiguous += 1
                continue
            oriented, off = got
            n_placed += 1
            placed_here += 1
            for j, base in enumerate(oriented):
                pos = off + j
                if -FLANK <= pos < len(ref) + FLANK and base in "ACGT":
                    pile[pos][base] += 1
        per_file[path.rsplit("/", 1)[-1]] = (len(reads), placed_here)

    print(f"reference      {name}  ({len(ref)} nt, {len(ref_aa)} aa)")
    print(f"reads read     {n_reads}")
    print(f"reads placed   {n_placed}"
          f"  ({n_ambiguous} unplaceable or ambiguous)")
    if not n_placed:
        sys.exit("\nNothing could be placed - check that these hit files came "
                 "from a MODE=aa run against this same reference.")

    def call(i):
        c = pile.get(i)
        d = sum(c.values()) if c else 0
        return (c.most_common(1)[0][0] if d >= MIN_DEPTH else "N"), d

    # Consensus over the reference CDS, which is what the identity stats use.
    cons, depths, low = [], [], 0
    for i in range(len(ref)):
        b, d = call(i)
        cons.append(b)
        depths.append(d)
        if b == "N":
            low += 1
    consensus = "".join(cons)

    called = len(ref) - low
    mean_depth = sum(depths) / len(depths)
    print(f"positions      {called}/{len(ref)} called at depth >= {MIN_DEPTH}"
          f"  ({called / len(ref):.0%})")
    print(f"mean depth     {mean_depth:.0f}x")

    # Flanks: how far the overhangs reach, contiguously, either side.
    up = []
    for i in range(-1, -FLANK - 1, -1):
        b, _ = call(i)
        if b == "N":
            break
        up.append(b)
    upstream = "".join(reversed(up))
    down = []
    for i in range(len(ref), len(ref) + FLANK):
        b, _ = call(i)
        if b == "N":
            break
        down.append(b)
    downstream = "".join(down)
    if FLANK:
        print(f"flanks         {len(upstream)} nt upstream, "
              f"{len(downstream)} nt downstream recovered")

    # Divergence from the reference, and what it means
    comparable = [(a, b) for a, b in zip(ref, consensus) if b != "N"]
    if comparable:
        nt_id = sum(a == b for a, b in comparable) / len(comparable)
        print(f"nt identity    {nt_id:.1%} to the reference over called positions")

    cons_aa = translate(consensus)
    aa_pairs = [(a, b) for a, b in zip(ref_aa, cons_aa) if b not in "X*"]
    if aa_pairs:
        aa_id = sum(a == b for a, b in aa_pairs) / len(aa_pairs)
        print(f"aa identity    {aa_id:.1%} to the reference protein")
        diffs = [(i + 1, a, b) for i, (a, b) in enumerate(zip(ref_aa, cons_aa))
                 if b not in "X*" and a != b]
        if diffs:
            shown = ", ".join(f"{a}{i}{b}" for i, a, b in diffs[:12])
            more = f" (+{len(diffs) - 12} more)" if len(diffs) > 12 else ""
            print(f"aa differences {shown}{more}")
        else:
            print("aa differences none - identical protein, different codons")

    chrom = find_chromophore(cons_aa)
    if chrom:
        print(f"chromophore    {chrom[1]} at residue {chrom[0]} - GFP confirmed")
    else:
        print("chromophore    NOT FOUND in the consensus - treat with suspicion")

    # Coverage gaps are what decide whether this is usable as a reference.
    gaps, start = [], None
    for i, d in enumerate(depths):
        if d < MIN_DEPTH and start is None:
            start = i
        elif d >= MIN_DEPTH and start is not None:
            gaps.append((start + 1, i))
            start = None
    if start is not None:
        gaps.append((start + 1, len(ref)))
    if gaps:
        print(f"gaps           {len(gaps)}: " +
              ", ".join(f"{a}-{b}" for a, b in gaps[:8]) +
              (" ..." if len(gaps) > 8 else ""))

    print("\nper file (reads / placed)")
    for f, (tot, pl) in sorted(per_file.items()):
        print(f"  {f:<28} {tot:>5} / {pl}")

    # The construct's own ORF, found in the flank-extended sequence. A
    # reference with a different N-terminus puts the real start codon outside
    # its own coordinates, so this is where an inserted residue shows up - and
    # where the record to actually use comes from.
    full = upstream + consensus + downstream
    if FLANK and "N" not in full:
        found = longest_orf(full)
        if found:
            cds, strand, start = found
            rel = start - len(upstream)
            caa = translate(cds)
            print(f"\nconstruct ORF  {len(cds)} nt ({len(caa) - 1} aa) on the "
                  f"{strand} strand, starting {rel:+d} nt relative to the "
                  "reference CDS")
            print(f"    starts     {caa[:14]!r}")
            c = find_chromophore(caa)
            if c:
                print(f"    chromophore {c[1]} at residue {c[0]}")
            if len(cds) != len(ref):
                print(f"    NOTE       {abs(len(cds) - len(ref))} nt "
                      f"{'longer' if len(cds) > len(ref) else 'shorter'} than "
                      "the reference - an insertion or deletion relative to it, "
                      "which is why a residue may read as substituted")
            consensus, cds_offset = cds, rel
        else:
            cds_offset = 0
    else:
        cds_offset = 0

    if upstream:
        print(f"\nupstream  ...{upstream[-60:]}")
    if downstream:
        print(f"downstream   {downstream[:60]}...")

    os.makedirs(os.path.dirname(OUT) or ".", exist_ok=True)
    label = "construct_GFP_recovered"
    with open(OUT, "w") as fh:
        fh.write(f">{label} offset={cds_offset:+d} depth={mean_depth:.0f}x\n")
        for i in range(0, len(consensus), 60):
            fh.write(consensus[i:i + 60] + "\n")
    print(f"\nwrote {OUT}  ({len(consensus)} nt)")

    if low:
        print(f"\n{low} position(s) are N. Raise depth by re-running the screen "
              "with KEEP_HITS=0 (keep every matching read) and, if needed, a "
              "larger MAX_READS, then pool more samples here.")
    else:
        print("\nFully called. Validate it, then use it as the reporter record "
              "in recombinant_model.fa:\n"
              f"    python3 -I gfp_tools.py validate {OUT}")


if __name__ == "__main__":
    main()
