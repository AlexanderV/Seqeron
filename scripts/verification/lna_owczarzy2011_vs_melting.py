#!/usr/bin/env python3
"""Check the Owczarzy 2011 LNA tables of PrimerDesigner.cs against MELTING 5.2.0's data files, and recompute the Tm of the
counterexample of PrimerProbeProperties.LnaTm_AddingLnaCPosition_NeverLowersTm independently of the C# code.

Usage:  python3 -I scripts/verification/lna_owczarzy2011_vs_melting.py <data_dir> <repo_root>

<data_dir> must hold MELTING 5.2.0 data files (data only, nothing is executed), e.g. downloaded from
  https://raw.githubusercontent.com/mohakjain/TmCalculator/master/MELTING5.2.0/Data/<file>
  Owczarzy2011lockedmn.xml  (single LNA)   Owczarzy2011lockedTandemmn.xml  (two consecutive LNAs)
  Owczarzy2011lockedmmn.xml (LNA + mismatch)
and saves the first one as owc11_single.xml (or pass the original name via the OWC11_SINGLE env var).

Part 1: every entry of the three C# dictionaries equals the XML entry (XML enthalpy is cal/mol, the code uses kcal/mol).
Part 2: CCCCGCTCGCCCCGCC with LNA at {1,2,3,4,5,7,8,10,11,12,13,14} vs the same plus 9. Every step has a locked base, so
dH and dS are the sum of XML step values plus the SantaLucia 1998 initiation of two terminal G.C pairs (+0.1 kcal/mol,
-2.8 cal/mol/K each); Tm = 1000 dH / (dS + R ln(CT/4)) with R = 1.9872, CT = 0.5 uM, then the Owczarzy 2004 monovalent
salt correction for 50 mM Na+. Expected: dH -179.909 -> -195.320, dS -435.374 -> -475.501, Tm 100.29881 -> 100.20802.
"""
import math, os, re, sys, xml.etree.ElementTree as ET

data_dir, repo = sys.argv[1], sys.argv[2]
cs = open(os.path.join(repo, "src/Seqeron/Algorithms/Seqeron.Genomics.MolTools/PrimerDesigner.cs"), encoding="utf-8").read()


def cs_table(name):
    m = re.search(r"Dictionary<string, \(double DeltaH, double DeltaS\)> %s = new\(StringComparer\.Ordinal\)\s*\{(.*?)\n    \};" % name, cs, re.S)
    out = {}
    for k, h, s in re.findall(r'\["([^"]+)"\]\s*=\s*\(\s*([-0-9.]+)\s*,\s*([-0-9.]+)\s*\)', m.group(1)):
        out[k] = (float(h), float(s))          # a later duplicate wins, like a Dictionary indexer assignment
    return out


def xml_table(fn):
    root = ET.parse(os.path.join(data_dir, fn)).getroot()   # data file, parsed as data only
    return {e.get("sequence"): (float(e.find("enthalpy").text) / 1000.0, float(e.find("entropy").text)) for e in root.iter("modified")}


single_file = os.environ.get("OWC11_SINGLE", "owc11_single.xml")
ok = True
for cname, fn in (("Owczarzy2011LnaSingle", single_file), ("Owczarzy2011LnaConsecutive", "Owczarzy2011lockedTandemmn.xml"),
                  ("Owczarzy2011LnaMismatch", "Owczarzy2011lockedmmn.xml")):
    c, x = cs_table(cname), xml_table(fn)
    diffs = [k for k in c if k in x and (abs(c[k][0] - x[k][0]) > 1e-9 or abs(c[k][1] - x[k][1]) > 1e-9)]
    extra = sorted(set(c) ^ set(x))
    print(f"{cname}: code {len(c)}, xml {len(x)}, key mismatches {len(extra)}, value differences {len(diffs)}")
    ok &= not diffs and not extra

single, tandem = xml_table(single_file), xml_table("Owczarzy2011lockedTandemmn.xml")
seq = "CCCCGCTCGCCCCGCC"
comp = "".join({"A": "T", "T": "A", "C": "G", "G": "C"}[c] for c in seq)
base = {1, 2, 3, 4, 5, 7, 8, 10, 11, 12, 13, 14}
res = {}
for name, s in (("before", base), ("after", base | {9})):
    locked = [i in s for i in range(len(seq))]
    H = S = 0.0
    for i in range(len(seq) - 1):
        a, b = locked[i], locked[i + 1]
        key = seq[i] + ("L" if a else "") + seq[i + 1] + ("L" if b else "") + "/" + comp[i] + comp[i + 1]
        v = tandem[key] if (a and b) else single[key]       # every step of this example has at least one locked base
        H += v[0]; S += v[1]
    H += 2 * 0.1; S += 2 * -2.8
    R, CT, Na = 1.9872, 0.5e-6, 0.05
    tm1 = H * 1000.0 / (S + R * math.log(CT / 4.0))
    f = sum(c in "GC" for c in seq) / len(seq)
    ln = math.log(Na)
    tm = 1.0 / (1.0 / tm1 + (4.29 * f - 3.95) * 1e-5 * ln + 9.40e-6 * ln * ln) - 273.15
    res[name] = (H, S, tm)
    print(f"{name}: dH={H:.3f} dS={S:.3f} Tm={tm:.5f}")
print(f"dTm = {res['after'][2] - res['before'][2]:+.5f}")
sys.exit(0 if ok else 1)
