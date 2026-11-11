#!/usr/bin/env python3
"""Round 2 of the E4 audit: the survivors of round 1 that are not equivalent (B06, E06) against the tests added after round 1, the replacements of the mutants that did not compile (B12, B26, E07, E08, E09, P09)
and mutants of the exit-code helper `E4.Healthy` (the condition E05 was moved into a function so that it can be tested). usage: make_e4_round2_mutants.py <repo root>"""
import json, os, sys

ROOT = sys.argv[1] if len(sys.argv) > 1 else '/home/user/Seqeron'
APP = 'apps/SuffixTree.EmProto/'
OUT = os.path.join(ROOT, 'docs/refactoring/suffix-tree-upgrade/analysis/mutation')
M = {}


def add(file, id_, desc, old, new):
    M.setdefault(file, []).append({'id': id_, 'desc': desc, 'old': old, 'new': new})


b = APP + 'Ph4/LcpBuild.cs'
add(b, 'B06', 'B2: a chunk with too few tuples is accepted (round 1 survivor, re-run)', 'if (reader.Count != length) throw new InvalidOperationException($"text chunk {c} holds {reader.Count} Phi tuples, expected {length}: the suffix array is not a permutation");\n                int got;\n                while ((got = reader.Read(read.AsSpan(readEntries))) > 0)\n                {\n                    for (int i = 0; i < got; i++)\n                    {\n                        PhiTuple t',
    'if (reader.Count > length) throw new InvalidOperationException($"text chunk {c} holds {reader.Count} Phi tuples, expected {length}: the suffix array is not a permutation");\n                int got;\n                while ((got = reader.Read(read.AsSpan(readEntries))) > 0)\n                {\n                    for (int i = 0; i < got; i++)\n                    {\n                        PhiTuple t')
add(b, 'B06b', 'bounded B2: a chunk with too few tuples is accepted', 'if (reader.Count != length) throw new InvalidOperationException($"text chunk {c} holds {reader.Count} Phi tuples, expected {length}: the suffix array is not a permutation");\n                int got;\n                while ((got = reader.Read(read.AsSpan(readEntries))) > 0)\n                {\n                    for (int i = 0; i < got; i++)\n                    {\n                        PhiPair t',
    'if (reader.Count > length) throw new InvalidOperationException($"text chunk {c} holds {reader.Count} Phi tuples, expected {length}: the suffix array is not a permutation");\n                int got;\n                while ((got = reader.Read(read.AsSpan(readEntries))) > 0)\n                {\n                    for (int i = 0; i < got; i++)\n                    {\n                        PhiPair t')
add(b, 'B12b', 'B3: the consumer is told the wrong first rank', '_consumer?.Consume(firstRank, dense);', '_consumer?.Consume(firstRank + 1, dense);')
add(b, 'B26b', 'Merge: the restart cost of the last task replaces the sum', '_restart += t.Restart;', '_restart = t.Restart;')

p = APP + 'Ph4/PlcpScan.cs'
add(p, 'P09b', 'PLCP loop: every value >= 255 counted twice', '                e8++;\n', '                e8 += 2;\n')

e = APP + 'Experiments/E4.cs'
add(e, 'E06', 'check-lcp: structure scan accepts a missing pair (round 1 survivor, re-run)', 'if (width == 4 ? overflowCount != 0 : escapesInBase != overflowCount) problems++;', '')
add(e, 'E07b', 'check-lcp: only a missing pair counts as a sample mismatch', 'if (got != expected) mismatches++;', 'if (got != expected && got < 0) mismatches++;')
add(e, 'E08b', 'check-lcp: a pair ahead of its slot is accepted', 'if (pr != rank || pr <= previousRank || pv < escape) problems++;', 'if (pr < rank || pr <= previousRank || pv < escape) problems++;')
add(e, 'E09b', 'check-lcp: a pair value below the escape is accepted', 'if (pr != rank || pr <= previousRank || pv < escape) problems++;', 'if (pr != rank || pr <= previousRank || pv < 1) problems++;')
add(e, 'H01', 'Healthy: the cap itself is exceeded', '=> rssAnonMb <= capMb && scratchEntriesLeft == 0', '=> rssAnonMb < capMb && scratchEntriesLeft == 0')
add(e, 'H02', 'Healthy: a leftover scratch entry is accepted (E05)', 'rssAnonMb <= capMb && scratchEntriesLeft == 0 && scratchLiveBytes == 0', 'rssAnonMb <= capMb && scratchLiveBytes == 0')
add(e, 'H03', 'Healthy: tracked scratch bytes left are accepted', 'scratchLiveBytes == 0 && arenaLiveBytes == 0;', 'arenaLiveBytes == 0;')
add(e, 'H04', 'Healthy: native bytes still charged are accepted', 'scratchLiveBytes == 0 && arenaLiveBytes == 0;', 'scratchLiveBytes == 0;')
add(e, 'H05', 'Healthy: the RssAnon cap is not checked', '=> rssAnonMb <= capMb && scratchEntriesLeft == 0', '=> scratchEntriesLeft == 0')

names = {APP + 'Ph4/LcpBuild.cs': 'build', APP + 'Ph4/PlcpScan.cs': 'plcpscan', APP + 'Experiments/E4.cs': 'e4'}
for f, items in M.items():
    src = open(os.path.join(ROOT, f)).read()
    for it in items:
        c = src.count(it['old'])
        if c != 1:
            print(f"PATTERN {it['id']} occurs {c} times in {f}")
            sys.exit(1)
    json.dump(items, open(os.path.join(OUT, f'e4_round2_{names[f]}_mutants.json'), 'w'), indent=1)
print(sum(len(v) for v in M.values()), 'mutants in', len(M), 'files')
