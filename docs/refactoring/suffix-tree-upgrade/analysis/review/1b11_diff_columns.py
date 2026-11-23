# Compares two generated COMPETITOR_MATRIX.md files (old, new): every cell outside the Seqeron column must be equal (1b.11). usage: python3 -I 1b11_diff_columns.py OLD.md NEW.md
import re, sys
def read(p): return open(p, encoding="utf-8").read().split("\n")
old, new = read(sys.argv[1]), read(sys.argv[2])
def cells(l): return [c.strip() for c in l.strip().strip("|").split("|")]
def sections(lines):
    s, cur = {}, None
    for l in lines:
        m = re.match(r"## (\d)\.", l)
        if m: cur = m.group(1); s[cur] = []
        elif cur: s[cur].append(l)
        else: s.setdefault("0", []).append(l)
    return s
so, sn = sections(old), sections(new)
print("sections:", sorted(so), sorted(sn))
# sections 0 (preamble), 1 (competitors), 6 (evidence) must be byte-identical
for k in ("0", "1", "6"):
    print("section", k, "identical:", so[k] == sn[k])
# section 3 and 4: per-row compare, Seqeron column index 3 (props) / 2 (funcs)
def rows(lines, prefix):
    return [cells(l) for l in lines if l.startswith("| " + prefix) or "weighted mean" in l or "group coverage" in l or "TOTAL coverage" in l]
changed = {"3": [], "4": []}; other = 0
for sec, col0 in (("3", 3), ("4", 2)):
    a_rows = [cells(l) for l in so[sec] if l.startswith("|") and not l.startswith("|---")]
    b_rows = [cells(l) for l in sn[sec] if l.startswith("|") and not l.startswith("|---")]
    assert len(a_rows) == len(b_rows), (sec, len(a_rows), len(b_rows))
    for a, b in zip(a_rows, b_rows):
        assert len(a) == len(b)
        for i, (x, y) in enumerate(zip(a, b)):
            if x != y:
                if i == col0: changed[sec].append((a[0] or a[1], x, y))
                else: other += 1; print("OTHER DIFF", sec, a[:2], i, x, y)
    # non-table lines identical?
    na = [l for l in so[sec] if not l.startswith("|")]; nb = [l for l in sn[sec] if not l.startswith("|")]
    print("section", sec, "non-table lines identical:", na == nb)
print("changed Seqeron cells: properties/summary rows", len(changed["3"]), "| functions/summary rows", len(changed["4"]), "| changed cells in other columns:", other)
# section 2: ranking
def rank(lines):
    r = {}
    for l in lines:
        m = re.match(r"\| (\d+) \| (.+?) \| ([\d.]+) \| ([\d.]+)% \|", l)
        if m: r[m.group(2).replace(" **(this project)**", "")] = (m.group(3), m.group(4), int(m.group(1)))
    return r
ro, rn = rank(so["2"]), rank(sn["2"])
for t in ro:
    if t != "Seqeron":
        if ro[t][:2] != rn[t][:2]: print("OTHER DIFF rank values", t, ro[t], rn[t])
        elif ro[t][2] != rn[t][2]: print("  rank number only:", t, ro[t][2], "->", rn[t][2])
print("Seqeron ranking line:", ro["Seqeron"], "->", rn["Seqeron"])
# section 2 text lines: show the differing ones
for a, b in zip([l for l in so["2"] if not l.startswith("|")], [l for l in sn["2"] if not l.startswith("|")]):
    if a != b: print("TEXT LINE DIFF (section 2):\n  old:", a[:260], "\n  new:", b[:260])
# section 5: gap list
ga = [cells(l) for l in so["5"] if l.startswith("| F")]; gb = [cells(l) for l in sn["5"] if l.startswith("| F")]
print("gap list entries old/new:", len(ga), len(gb))
# competitor columns in gap list: 'best' value can only change where Seqeron changed
fa = {g[0]: g for g in ga}; fb = {g[0]: g for g in gb}
for f in sorted(set(fa) | set(fb)):
    if f not in fb: print("  gap entry removed:", f, fa[f][1], "Seqeron", fa[f][2], "best", fa[f][3])
    elif f not in fa: print("  gap entry added:", f)
    elif fa[f] != fb[f]: print("  gap entry changed:", f, fa[f][2:], "->", fb[f][2:])
print("--- changed Seqeron cells")
for sec in ("3", "4"):
    for r in changed[sec]: print(sec, r)
