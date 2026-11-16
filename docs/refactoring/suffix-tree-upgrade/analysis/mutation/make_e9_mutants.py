#!/usr/bin/env python3
"""Writes the hand-made single-site mutants of the new production code of E9 (apps/SuffixTree.EmProto/Ph6: PsiStream, BwtPass, IsaByDist and PsiFromIsa, PsiFiles, LinkSmoke) as one JSON list per file
   (e9_<file>_listed.json: every candidate written; e9_<file>_mutants.json: the ones handed to run_e9_mutants.py, here the same list). Each pattern must occur exactly once in its file (checked here).
   Round 1 is written against the sources of commit 8d79fc47 (the code the audit ran on); round 2 (mutants that did not compile or survived round 1, run after the tests were strengthened and LinkSmoke gained its two counters) against the working tree.
   usage: make_e9_mutants.py <repo root> [round1-source-revision, default 8d79fc47]"""
import json, os, sys

ROOT = sys.argv[1] if len(sys.argv) > 1 else '/home/user/Seqeron'
REV = sys.argv[2] if len(sys.argv) > 2 else '8d79fc47'
import subprocess
PH6 = 'apps/SuffixTree.EmProto/Ph6/'

M = {}
def m(file, id, desc, old, new):
    M.setdefault(file, []).append({'id': id, 'desc': desc, 'old': old, 'new': new})

# ---- PsiStream.cs: the sigma-way distribution, C[], the sentinel row, the multi-pass loop
f = 'PsiStream.cs'
m(f, 'P01', 'C[] base without the empty suffix: start[0] = 0', 'start[0] = 1;', 'start[0] = 0;')
m(f, 'P02', 'C[] ignores the counts: start[c+1] = start[c]', 'start[c + 1] = start[c] + counts[c];', 'start[c + 1] = start[c];')
m(f, 'P03', 'C[] adds one per symbol: start[c+1] = start[c] + counts[c] + 1', 'start[c + 1] = start[c] + counts[c];', 'start[c + 1] = start[c] + counts[c] + 1;')
m(f, 'P04', 'sentinel row: skip index off by one (primary - first + 1)', 'long skip = meta.PrimaryRow - first;', 'long skip = meta.PrimaryRow - first + 1;')
m(f, 'P05', 'sentinel row at the start of a block is not skipped (skip > 0)', 'if (skip >= 0 && skip < take)', 'if (skip > 0 && skip < take)')
m(f, 'P06', 'rows after the sentinel get the row number of the sentinel (firstRow of the second scan without +1)', 'cursors.Scan(b + skip + 1, take - (int)skip - 1, first + skip + 1);', 'cursors.Scan(b + skip + 1, take - (int)skip - 1, first + skip);')
m(f, 'P07', 'the first scan includes the sentinel row (skip + 1)', 'cursors.Scan(b, (int)skip, first);', 'cursors.Scan(b, (int)skip + 1, first);')
m(f, 'P08', 'the second scan starts at the sentinel row (no +1 in the pointer)', 'cursors.Scan(b + skip + 1, take - (int)skip - 1, first + skip + 1);', 'cursors.Scan(b + skip, take - (int)skip, first + skip + 1);')
m(f, 'P09', 'symbol of the pass off by one (bytes[j] - lo + 1)', 'uint s = (uint)(bytes[j] - lo);', 'uint s = (uint)(bytes[j] - lo + 1);')
m(f, 'P10', 'row written off by one (firstRow + j + 1)', '*p = (uint)(firstRow + j);', '*p = (uint)(firstRow + j + 1);')
m(f, 'P11', 'cursor advances by two entries', 'cur[s] = (nuint)(p + 1);', 'cur[s] = (nuint)(p + 2);')
m(f, 'P12', 'the cursor window is the larger of the window and the rest (Max)', 'int entries = (int)Math.Min(_entriesPerWindow, remaining);', 'int entries = (int)Math.Max(_entriesPerWindow, remaining);')
m(f, 'P13', 'the next window of a cursor starts one entry later', '_next[s] += entries;', '_next[s] += entries + 1;')
m(f, 'P14', 'a cursor that has used its range is not refused (remaining < 0)', 'if (remaining <= 0) throw new InvalidDataException', 'if (remaining < 0) throw new InvalidDataException')
m(f, 'P15', 'under-full cursors are not detected inside the last window', 'if (_next[s] != _limit[s] || _cur[s] != _end[s])', 'if (_next[s] != _limit[s])')
m(f, 'P16', 'windows not fully mapped are not detected', 'if (_next[s] != _limit[s] || _cur[s] != _end[s])', 'if (_cur[s] != _end[s])')
m(f, 'P17', 'mapped bound without the factor 4 (entries instead of bytes)', 'long bound = group * (long)entriesPerWindow * 4;', 'long bound = group * (long)entriesPerWindow;')
m(f, 'P18', 'passes by floor division: the last partial group is dropped', '(sigma + o.SymbolsPerPass - 1) / o.SymbolsPerPass;', 'sigma / o.SymbolsPerPass;')
m(f, 'P19', 'the last group runs past sigma (hi without Math.Min)', 'int hi = Math.Min(sigma, lo + o.SymbolsPerPass);', 'int hi = lo + o.SymbolsPerPass;')
m(f, 'P20', 'the group start uses the open cursors instead of the limit (pass * group)', 'int lo = pass * o.SymbolsPerPass;', 'int lo = pass * group;')
m(f, 'P21', 'convention swapped: entry 0 is the primary row unless cyclic', 'o.CyclicEntryZero ? (uint)meta.PrimaryRow : 0u', 'o.CyclicEntryZero ? 0u : (uint)meta.PrimaryRow')
m(f, 'P22', 'cyclic entry 0 off by one', 'o.CyclicEntryZero ? (uint)meta.PrimaryRow : 0u', 'o.CyclicEntryZero ? (uint)meta.PrimaryRow + 1 : 0u')
m(f, 'P23', 'primary row == n + 1 accepted', 'meta.PrimaryRow >= rows', 'meta.PrimaryRow > rows')
m(f, 'P24', 'counts that sum to more than n accepted', 'counts.Sum() != n', 'counts.Sum() < n')
m(f, 'P25', 'a BWT file longer than n + 1 accepted', 'if (bwt.Length != rows)', 'if (bwt.Length < rows)')
m(f, 'P26', 'rows = n', 'long rows = n + 1;', 'long rows = n;')
m(f, 'P27', 'zero symbols per pass accepted', 'ArgumentOutOfRangeException.ThrowIfLessThan(o.SymbolsPerPass, 1);', 'ArgumentOutOfRangeException.ThrowIfLessThan(o.SymbolsPerPass, 0);')
m(f, 'P28', 'window below four bytes gives zero entries', 'int entriesPerWindow = Math.Max(1, o.WindowBytes / 4);', 'int entriesPerWindow = o.WindowBytes / 4;')
m(f, 'P29', 'read buffer not bounded by the number of rows', 'int readLen = (int)Math.Min(o.ReadBytes, rows);', 'int readLen = o.ReadBytes;')
m(f, 'P30', 'the BWT is read from the wrong offset', 'bwt.ReadExactAt(block, first);', 'bwt.ReadExactAt(block, first + 1);')
m(f, 'P31', 'window count of the passes is not accumulated', 'windows += cursors.Windows;', 'windows = cursors.Windows;')
m(f, 'P32', 'BWT bytes read reported for one pass', 'passes * rows, 4 * rows);', 'rows, 4 * rows);')
m(f, 'P33', 'the window cap is not charged (reservation of a different size)', 'new ArenaBudget(o.WindowCapBytes).Reserve(bound, "psi.windows");', 'new ArenaBudget(o.WindowCapBytes).Reserve(bound / 2, "psi.windows");')
m(f, 'P34', 'the shared budget is not released', 'budget.Release(bound);', 'budget.Release(0);')
m(f, 'P35', 'symbols of other passes are counted as this pass (s > span)', 'if (s >= span) continue; // a symbol of another pass', 'if (s > span) continue; // a symbol of another pass')
m(f, 'P36', 'the first entry of the cursor range is entry start (without the +1 of the empty suffix) via counts offset', 'start[c + 1] = start[c] + counts[c];\n', 'start[c + 1] = start[c] + counts[c] - (c == 0 ? 1 : 0);\n')

# ---- BwtPass.cs
f = 'BwtPass.cs'
m(f, 'B01', 'the symbol before the suffix is taken at p instead of p - 1', 'uint s = text.Symbol(p - 1);', 'uint s = text.Symbol(p);')
m(f, 'B02', 'the sentinel row is the row of suffix 1', 'if (p == 0)\n', 'if (p == 1)\n')
m(f, 'B03', 'counts add two', 'counts[s]++;', 'counts[s] += 2;')
m(f, 'B04', 'row 0 is not checked to be the empty suffix', 'if (first + i == 0 && p != n)', 'if (first + i == 0 && p == n)')
m(f, 'B05', 'primary row off by one', 'primary = first + i;', 'primary = first + i + 1;')
m(f, 'B06', 'symbol equal to sigma accepted', 'if (s >= sigma)', 'if (s > sigma)')
m(f, 'B07', 'an SA entry equal to n + 1 accepted', 'if (p > n) throw', 'if (p > n + 1) throw')
m(f, 'B08', 'two rows with suffix 0 accepted', 'if (primary >= 0) throw', 'if (primary >= 1) throw')
m(f, 'B09', 'the dummy byte of the sentinel row is 1', 'output[i] = 0;\n                        continue;', 'output[i] = 1;\n                        continue;')
m(f, 'B10', 'SA file of the wrong length accepted (<)', 'if (sa.Length != 4 * rows)', 'if (sa.Length < 4 * rows)')

# ---- IsaByDist.cs
f = 'IsaByDist.cs'
m(f, 'D01', 'chunk of a position off by one', 'long chunk = p / chunkPositions;', 'long chunk = (p + 1) / chunkPositions;')
m(f, 'D02', 'rank in the tuple off by one', 'buf[filled[chunk]++] = ((ulong)p << 32) | (uint)(first + i);', 'buf[filled[chunk]++] = ((ulong)p << 32) | (uint)(first + i + 1);')
m(f, 'D03', 'buffer entries from 4-byte tuples', 'int bufEntries = o.IoBufferBytes / 8;', 'int bufEntries = o.IoBufferBytes / 4;')
m(f, 'D04', 'chunk count by floor division', 'long chunks = (rows + chunkPositions - 1) / chunkPositions;', 'long chunks = rows / chunkPositions;')
m(f, 'D05', 'the chunk array is not initialised to the missing marker', 'array.Fill(uint.MaxValue);', 'array.Clear();')
m(f, 'D06', 'the last chunk array is as long as a full chunk', 'int length = (int)Math.Min(chunkPositions, rows - startPos);', 'int length = (int)Math.Min(chunkPositions, rows);')
m(f, 'D07', 'duplicate positions are not detected when they occur', 'if (slot != uint.MaxValue) throw new InvalidDataException($"the suffix array is not a permutation: position {(long)(v >> 32)} occurs twice");', '')
m(f, 'D08', 'missing positions are not detected', 'if (array[k] == uint.MaxValue) throw new InvalidDataException', 'if (array[k] == uint.MaxValue - 1) throw new InvalidDataException')
m(f, 'D09', 'entries beyond n accepted', 'if (p > n) throw new InvalidDataException($"SA[{first + i}] = {p} lies beyond n = {n}");', '')
m(f, 'D10', 'the tuples of a flushed buffer are written at the start of the file', 'f.WriteAt(MemoryMarshal.AsBytes(buf[..filled[chunk]]), flushed[chunk] * 8);', 'f.WriteAt(MemoryMarshal.AsBytes(buf[..filled[chunk]]), 0);')
m(f, 'D11', 'the chunk files are read from the wrong offset', 'f.ReadExactAt(MemoryMarshal.AsBytes(part), at * 8);', 'f.ReadExactAt(MemoryMarshal.AsBytes(part), at * 8 + 8);')
m(f, 'D12', 'the in-memory remainder of a chunk is not scattered', 'Scatter(all.Slice(checked((int)(chunk * bufEntries)), filled[chunk]), array, startPos);', 'Scatter(all.Slice(checked((int)(chunk * bufEntries)), 0), array, startPos);')

# ---- PsiFromIsa (IsaByDist.cs) and PsiFiles.cs
f = 'IsaByDist.cs'
m(f, 'F01', 'Psi reads ISA[SA[i]] instead of ISA[SA[i] + 1]', 'output[i] = table[p + 1];', 'output[i] = table[p];')
m(f, 'F02', 'Psi entry 0 is 1', 'output[i] = 0;\n                        continue;', 'output[i] = 1;\n                        continue;')
m(f, 'F03', 'SA entry n at a row >= 1 accepted', 'if (p >= n) throw new InvalidDataException', 'if (p > n) throw new InvalidDataException')
f = 'PsiFiles.cs'
m(f, 'L01', 'LoadForP1 keeps entry 0', 'if (psi.Length > 0) psi[0] = -1;', 'if (psi.Length > 0) psi[0] = psi[0];')
m(f, 'L02', 'the oracle file stores P1\'s -1 in entry 0', 'if (psi.Length > 0) psi[0] = 0;', 'if (psi.Length > 0) psi[0] = psi[0];')
m(f, 'L03', 'a ragged file is accepted', 'if (fs.Length % 4 != 0) throw', 'if (fs.Length % 4 != 5) throw')

# ---- LinkSmoke.cs
f = 'LinkSmoke.cs'
m(f, 'S01', 'sample A of P1: half without the -4000', 'int half = Math.Min(sample, Math.Max(1, (drawn.Length - 4000) / 2));', 'int half = Math.Min(sample, Math.Max(1, drawn.Length / 2));')
m(f, 'S02', 'step counts are not compared', 'if (threw || okA != okB || stepsA != stepsB || (okA && a != b)) mismatches++;', 'if (threw || okA != okB || (okA && a != b)) mismatches++;')
m(f, 'S03', 'targets are not compared', 'if (threw || okA != okB || stepsA != stepsB || (okA && a != b)) mismatches++;', 'if (threw || okA != okB || stepsA != stepsB) mismatches++;')
m(f, 'S04', 'kernel exceptions are not counted', 'if (threw || okA != okB', 'if (okA != okB')
m(f, 'S05', 'the empty-node guard of the P1 sampling is dropped', 'if (nodes.Count > 0 && p1Sampling)', 'if (p1Sampling)')
m(f, 'S06', 'the definition is compared with the other run only', 'else if (expected is not null && a != expected[k]) defMismatches++;', 'else if (expected is not null && false) defMismatches++;')

def source(f, rev):
    if rev is None:
        return open(os.path.join(ROOT, PH6, f)).read()
    return subprocess.run(['git', '-C', ROOT, 'show', '%s:%s%s' % (rev, PH6, f)], capture_output=True, text=True, check=True).stdout

def write(table, prefix, rev):
    for f, lst in table.items():
        src = source(f, rev)
        for x in lst:
            c = src.count(x['old'])
            if c != 1:
                print('PATTERN-NOT-UNIQUE', prefix, f, x['id'], c)
                sys.exit(1)
            x['file'] = PH6 + f
        base = os.path.join(os.path.dirname(os.path.abspath(__file__)), prefix + f.replace('.cs', '').lower())
        json.dump(lst, open(base + '_listed.json', 'w'), indent=1)

write(M, 'e9_', REV)
tot = sum(len(v) for v in M.values())
print('round 1 mutants listed:', tot, {k: len(v) for k, v in M.items()})

# ---- round 2: replacements of the mutants that did not compile (analyzer errors S1125, S1481, S1656), the mutants that survived with a test that kills them now, and S05 (its pattern had changed before the round-1 list was made)
R2 = {}
def r2(file, id, desc, old, new):
    R2.setdefault(file, []).append({'id': id, 'desc': desc, 'old': old, 'new': new})
r2('PsiStream.cs', 'P14', 'a cursor that has used its range is not refused (remaining < 0) [message of the overflow is asserted now]', 'if (remaining <= 0) throw new InvalidDataException', 'if (remaining < 0) throw new InvalidDataException')
r2('PsiStream.cs', 'P16', 'windows not fully mapped are not detected [windows of one entry are tested now]', 'if (_next[s] != _limit[s] || _cur[s] != _end[s])', 'if (_cur[s] != _end[s])')
r2('PsiStream.cs', 'P24', 'counts that sum to more than n accepted [N - 1 is tested now]', 'counts.Sum() != n', 'counts.Sum() < n')
r2('BwtPass.cs', 'B10', 'SA file of the wrong length accepted (<) [a longer file is tested now]', 'if (sa.Length != 4 * rows)', 'if (sa.Length < 4 * rows)')
r2('PsiFiles.cs', 'L01', 'LoadForP1 keeps entry 0 as the file has it (0): the mutant maps to 0 instead of -1', 'if (psi.Length > 0) psi[0] = -1;', 'if (psi.Length > 0) psi[0] = 0;')
r2('PsiFiles.cs', 'L02', 'the oracle file stores P1\'s -1 in entry 0', 'if (psi.Length > 0) psi[0] = 0;', 'if (psi.Length > 0) psi[0] = -1;')
r2('LinkSmoke.cs', 'S02', 'step counts compared one way only (<=)', 'bool sameRun = okA == okB && stepsA == stepsB;', 'bool sameRun = okA == okB && stepsA <= stepsB;')
r2('LinkSmoke.cs', 'S03', 'targets are not compared', 'if (threw || !sameRun || (okA && a != b)) mismatches++;', 'if (threw || !sameRun) mismatches++;')
r2('LinkSmoke.cs', 'S04', 'kernel exceptions are not counted as mismatches', 'if (threw || !sameRun || (okA && a != b)) mismatches++;', 'if (!sameRun || (okA && a != b)) mismatches++;')
r2('LinkSmoke.cs', 'S05', 'the empty-node guard of the P1 sampling is dropped', 'if (nodes.Count > 0 && p1Sampling)', 'if (p1Sampling)')
r2('LinkSmoke.cs', 'S06', 'the definition is compared for equality instead of difference', 'else if (expected is not null && a != expected[k]) defMismatches++;', 'else if (expected is not null && a == expected[k]) defMismatches++;')
write(R2, 'e9_round2_', None)

# ---- round 3: P16 again after its test asserts which check refuses (the symbol with the lowered count overflows in a later pass with the same exception type, so the type alone cannot tell the checks apart)
R3 = {'PsiStream.cs': [dict(R2['PsiStream.cs'][1], desc='windows not fully mapped are not detected [the message of the end-of-pass check is asserted for windows of one and two entries]')]}
write(R3, 'e9_round3_', None)
print('round 2 mutants listed:', sum(len(v) for v in R2.values()))
