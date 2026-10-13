#!/usr/bin/env bash
# Pull the GFP coding sequence out of a construct map into ref/gfp.fa.
#
# Usage:
#   ./00_make_gfp_fasta.sh my_construct.gb            # any CDS whose name looks like GFP
#   ./00_make_gfp_fasta.sh my_construct.gb eGFP       # a specific feature name
#
# Works on GenBank (.gb/.gbk) exports, which is what SnapGene, Benchling and
# Addgene all produce. If you already have a plain FASTA of the reporter, skip
# this and point GFP_FA at it.
set -euo pipefail
cd "$(dirname "$0")" && source ../config.sh

GB="${1:-}"
WANT="${2:-}"
: "${GFP_FA:=$PROJ/ref/gfp.fa}"

if [[ -z "$GB" || ! -s "$GB" ]]; then
    sed -n '2,12p' "$0" >&2
    exit 1
fi

# Write via a temp file so a failed extraction cannot truncate a good gfp.fa.
TMP_FA="$(mktemp)"
trap 'rm -f "$TMP_FA"' EXIT

python3 -I - "$GB" "$WANT" > "$TMP_FA" <<'PY'
import re, sys

path, want = sys.argv[1], sys.argv[2].lower()

text = open(path, errors="replace").read()

# ORIGIN block -> the plasmid sequence
m = re.search(r"\nORIGIN(.*?)\n//", text, re.S)
if not m:
    sys.exit(f"{path}: no ORIGIN block - is this really GenBank format?")
seq = re.sub(r"[^acgtnACGTN]", "", m.group(1)).upper()

# FEATURES block -> every CDS/gene with a location and a name
feats = []
fblock = text.split("FEATURES", 1)[-1].split("\nORIGIN", 1)[0]
cur = None
for line in fblock.splitlines():
    head = re.match(r"^ {5}(\S+)\s+(\S.*)$", line)
    if head:
        if cur:
            feats.append(cur)
        cur = {"type": head.group(1), "loc": head.group(2).strip(), "names": []}
    elif cur is not None:
        q = re.match(r'^\s+/(?:gene|label|product|note)="?([^"]*)"?', line)
        if q:
            cur["names"].append(q.group(1))
if cur:
    feats.append(cur)

GFPISH = re.compile(r"gfp|egfp|sfgfp|zsgreen|mgreen|green.?fluor", re.I)


def extract(loc):
    rc = "complement" in loc
    nums = [int(x) for x in re.findall(r"\d+", loc)]
    if len(nums) < 2:
        return None
    # join(...) spans are concatenated in the order written
    parts = [(nums[i], nums[i + 1]) for i in range(0, len(nums) - 1, 2)]
    s = "".join(seq[a - 1:b] for a, b in parts)
    if rc:
        s = s.translate(str.maketrans("ACGTN", "TGCAN"))[::-1]
    return s


hits = []
for f in feats:
    if f["type"] not in ("CDS", "gene", "misc_feature"):
        continue
    label = next((n for n in f["names"] if n), "")
    if want:
        ok = any(want in n.lower() for n in f["names"])
    else:
        ok = any(GFPISH.search(n) for n in f["names"])
    if ok:
        s = extract(f["loc"])
        if s:
            hits.append((label or "GFP", f["type"], s))

if not hits:
    names = sorted({n for f in feats for n in f["names"] if n})
    sys.exit("No GFP-like feature found. Features present:\n  " +
             "\n  ".join(names[:60]) +
             "\n\nRe-run with the exact name as the second argument.")

# Prefer a CDS over a bare annotation, and the longest if several match.
hits.sort(key=lambda h: (h[1] != "CDS", -len(h[2])))
name, ftype, s = hits[0]
clean = re.sub(r"[^A-Za-z0-9_.-]", "_", name)
print(f">{clean}")
for i in range(0, len(s), 60):
    print(s[i:i + 60])

print(f"took {ftype} '{name}', {len(s)} nt", file=sys.stderr)
if len(hits) > 1:
    print("other matches (ignored): " +
          ", ".join(f"{h[0]}({len(h[2])}nt)" for h in hits[1:]), file=sys.stderr)
PY

mkdir -p "$(dirname "$GFP_FA")"
mv "$TMP_FA" "$GFP_FA"

echo
echo "Wrote $GFP_FA"
echo
python3 -I gfp_tools.py validate "$GFP_FA"
