#!/usr/bin/env python3
"""Writes the hand-made mutant lists of the E4 audit (one JSON per production file; each mutant is a single-site edit: old text -> new text, `old` occurring exactly once in the file)
and checks that every pattern is unique in the current source. usage: make_e4_mutants.py <repo root>"""
import json, os, sys

ROOT = sys.argv[1] if len(sys.argv) > 1 else '/home/user/Seqeron'
APP = 'apps/SuffixTree.EmProto/'
OUT = os.path.join(ROOT, 'docs/refactoring/suffix-tree-upgrade/analysis/mutation')

M = {}


def add(file, id_, desc, old, new):
    M.setdefault(file, []).append({'id': id_, 'desc': desc, 'old': old, 'new': new})


# --- LcpGeometry.cs
g = APP + 'Ph4/LcpGeometry.cs'
add(g, 'G01', 'TextChunks: floor instead of ceiling', 'TextChunks = (n + chunkPositions - 1) / chunkPositions;', 'TextChunks = n / chunkPositions;')
add(g, 'G02', 'RankChunks: ceil(n / R) instead of ceil((n + 1) / R)', 'RankChunks = (n + chunkPositions) / chunkPositions;', 'RankChunks = (n + chunkPositions - 1) / chunkPositions;')
add(g, 'G03', 'TextChunkLength: the last chunk is as long as the others', 'Math.Min(R, N - (chunk * R));', 'R;')
add(g, 'G04', 'RankChunkLength: the empty suffix rank is not counted', 'Math.Min(R, N + 1 - (chunk * R));', 'Math.Min(R, N - (chunk * R));')
add(g, 'G05', 'the fast division is not used for R = 2', 'chunkPositions is >= 2 and <= uint.MaxValue;', 'chunkPositions is >= 3 and <= uint.MaxValue;')
add(g, 'G06', 'ChunkOf: R = 1 falls back to a single chunk', 'return R == 1 ? x : 0;', 'return 0;')
add(g, 'G07', 'a text of 2^32 - 1 symbols is accepted', 'if (n >= uint.MaxValue) throw', 'if (n > uint.MaxValue) throw')

# --- Spool.cs
s = APP + 'Ph4/Spool.cs'
add(s, 'S01', 'Append: the reserved offset is the end of the reservation', 'long at = Interlocked.Add(ref _entries[chunk], items.Length) - items.Length;', 'long at = Interlocked.Add(ref _entries[chunk], items.Length);')
add(s, 'S02', 'Append: the tracker is not told about the bytes', 'ScratchTracker.Add((long)items.Length * sizeof(T));', '')
add(s, 'S03', 'Delete: the tracker is not told about the deletion', 'ScratchTracker.Remove(count * sizeof(T));', '')
add(s, 'S04', 'Reader.Read: reads past the end of the chunk', 'int take = (int)Math.Min(buffer.Length, Count - _position);', 'int take = buffer.Length;')
add(s, 'S05', 'Reader.Read: the position is not advanced (loops for ever)', '_position += take;', '_position += 0;')
add(s, 'S06', 'Writer: a buffer holds zero tuples when bufferBytes < sizeof(T)', '_capacity = Math.Max(1, bufferBytes / sizeof(T));', '_capacity = bufferBytes / sizeof(T);')
add(s, 'S07', 'Writer.Add: flush one tuple too late (writes into the next buffer)', 'if (++filled == _capacity)', 'if (++filled > _capacity)')
add(s, 'S08', 'Writer.Flush: the first chunk of the writer is ignored (wrong file in later rounds)', '_spool.Append(_firstChunk + local,', '_spool.Append(local,')
add(s, 'S09', 'Spool.Dispose: the last chunk file is left', 'for (long c = 0; c < _entries.Length; c++) Delete(c);', 'for (long c = 0; c < _entries.Length - 1; c++) Delete(c);')
add(s, 'S10', 'Writer.FlushAll: the filled count is not reset (tuples written twice)', '            Flush(local, _filled[local]);\n            _filled[local] = 0;\n        }\n    }', '            Flush(local, _filled[local]);\n        }\n    }')
add(s, 'S11', 'OpenReader: a chunk with entries is treated as empty', 'return new SpoolReader<T>(count > 0 ? new CountedFile', 'return new SpoolReader<T>(count > 1 ? new CountedFile')

# --- PlcpScan.cs
p = APP + 'Ph4/PlcpScan.cs'
add(p, 'P01', 'limit from the smaller of the two positions', 'long limit = n - Math.Max(j, p);', 'long limit = n - Math.Min(j, p);')
add(p, 'P02', 'the carried length is not decremented', 'if (l > 0) l--;', 'if (l > 0) l -= 0;')
add(p, 'P03', 'escape >= 255 counted from 256', 'if (l >= 255)', 'if (l > 255)')
add(p, 'P04', 'escape >= 65535 counted from 65536', 'if (l >= 65535) e16++;', 'if (l > 65535) e16++;')
add(p, 'P05', 'FirstValue taken at the second position', 'if (k == 0) first = l;', 'if (k == 1) first = l;')
add(p, 'P06', 'the carry survives a position without a predecessor', '                l = 0; // j is the smallest nonempty suffix: its predecessor is the empty suffix\n', '')
add(p, 'P07', 'MaxPlcp keeps the smallest', 'if (l > max) max = l;', 'if (l < max) max = l;')
add(p, 'P08', 'a position without a predecessor keeps a stale value', '                phi[k] = 0;\n                continue;', '                continue;')
add(p, 'P09', 'the escapes are not counted when the value reaches 255', '                e8++;\n', '')

# --- LcpOutput.cs
o = APP + 'Ph4/LcpOutput.cs'
add(o, 'O01', 'Select: U8 only below n / 16', 'if (16 * escapes8 <= n) return LcpKind.U8;', 'if (16 * escapes8 < n) return LcpKind.U8;')
add(o, 'O02', 'Select: U16 only below n / 8', 'if (8 * escapes16 <= n) return LcpKind.U16;', 'if (8 * escapes16 < n) return LcpKind.U16;')
add(o, 'O03', 'Select: equal counts rejected', 'if (escapes16 > escapes8) throw', 'if (escapes16 >= escapes8) throw')
add(o, 'O04', 'Bytes: U16 entries have four bytes', 'LcpKind.U16 => 2, _ => 4 };\n\n    /// <summary>The escape', 'LcpKind.U16 => 4, _ => 4 };\n\n    /// <summary>The escape')
add(o, 'O05', 'Escape: U16 escape is 65534', 'LcpKind.U16 => ushort.MaxValue, _ => uint.MaxValue };', 'LcpKind.U16 => ushort.MaxValue - 1, _ => uint.MaxValue };')
add(o, 'O06', 'OverflowCount: U8 and U16 swapped', 'LcpKind.U8 => escapes8, LcpKind.U16 => escapes16, _ => 0 };', 'LcpKind.U8 => escapes16, LcpKind.U16 => escapes8, _ => 0 };')
add(o, 'O07', 'Write: ranks may arrive out of order', 'if (firstRank != _nextRank) throw', 'if (firstRank < _nextRank) throw')
add(o, 'O08', 'Write: the announced number of escapes is a bound that is reached', 'if (_overflowWritten + pairs > _overflowExpected) throw', 'if (_overflowWritten + pairs >= _overflowExpected) throw')
add(o, 'O09', 'Complete: the last rank is not required', 'if (_nextRank != n + 1) throw', 'if (_nextRank != n) throw')
add(o, 'O10', 'Convert U8: the value 255 is not an escape', '                    if (v >= _escape)\n                    {\n                        list[2 * pairs] = (uint)(firstRank + i);\n                        list[(2 * pairs) + 1] = v;\n                        pairs++;\n                        bytes[i] = (byte)_escape;',
    '                    if (v > _escape)\n                    {\n                        list[2 * pairs] = (uint)(firstRank + i);\n                        list[(2 * pairs) + 1] = v;\n                        pairs++;\n                        bytes[i] = (byte)_escape;')
add(o, 'O11', 'Convert: the rank of an escape is relative to the slice', 'list[2 * pairs] = (uint)(firstRank + i);\n                        list[(2 * pairs) + 1] = v;\n                        pairs++;\n                        u16[i]', 'list[2 * pairs] = (uint)i;\n                        list[(2 * pairs) + 1] = v;\n                        pairs++;\n                        u16[i]')
add(o, 'O12', 'Write: the last slice is a full slice', 'ReadOnlySpan<uint> part = values.Slice(at, Math.Min(_slice, values.Length - at));', 'ReadOnlySpan<uint> part = values.Slice(at, _slice);')
add(o, 'O13', 'the base file is one entry short', '_base = new MappedWindowWriter(basePath, (n + 1) * _width, windowBytes);', '_base = new MappedWindowWriter(basePath, n * _width, windowBytes);')
add(o, 'O14', 'the overflow file holds four bytes per entry', 'new MappedWindowWriter(overflowPath, overflowCount * 8, windowBytes);', 'new MappedWindowWriter(overflowPath, overflowCount * 4, windowBytes);')
add(o, 'O15', 'U32 values are not copied', 'part.CopyTo(MemoryMarshal.Cast<byte, uint>(bytes));', '')

# --- LcpBuild.cs
b = APP + 'Ph4/LcpBuild.cs'
add(b, 'B01', 'the file length check ignores the empty suffix entry', 'if (length != (n + 1) * 4) throw', 'if (length != n * 4) throw')
add(b, 'B02', 'sa[0] smaller than n accepted', 'if (first[0] != n) throw', 'if (first[0] > n) throw')
add(b, 'B03', 'B1: the predecessor of the first entry is 0', 'uint previous = (uint)n; // SA[0], the empty suffix\n        for (long first = 1; first <= n; first += readEntries)', 'uint previous = 0; // SA[0], the empty suffix\n        for (long first = 1; first <= n; first += readEntries)')
add(b, 'B04', 'B1: an entry equal to n accepted', 'if (s >= n) throw new ArgumentException($"sa[{first + i}] = {s} is outside [0, {n})");', 'if (s > n) throw new ArgumentException($"sa[{first + i}] = {s} is outside [0, {n})");')
add(b, 'B05', 'B1: the rank of the tuple is off by one', 'new PhiTuple(s, previous, (uint)(first + i))', 'new PhiTuple(s, previous, (uint)(first + i - 1))')
add(b, 'B06', 'B2: a chunk with too few tuples is accepted', 'if (reader.Count != length) throw new InvalidOperationException($"text chunk {c} holds {reader.Count} Phi tuples, expected {length}: the suffix array is not a permutation");\n                int got;\n                while ((got = reader.Read(read.AsSpan(readEntries))) > 0)\n                {\n                    for (int i = 0; i < got; i++)\n                    {\n                        PhiTuple t',
    'if (reader.Count > length) throw new InvalidOperationException($"text chunk {c} holds {reader.Count} Phi tuples, expected {length}: the suffix array is not a permutation");\n                int got;\n                while ((got = reader.Read(read.AsSpan(readEntries))) > 0)\n                {\n                    for (int i = 0; i < got; i++)\n                    {\n                        PhiTuple t')
add(b, 'B07', 'B2: the rank array is filled with the position', 'rank[(int)at] = t.Rank;', 'rank[(int)at] = t.Position;')
add(b, 'B08', 'B2: the Phi chunk file is not deleted', 'phiSpool.Delete(c); // the chunk is in memory now: its file is scratch no more', '')
add(b, 'B09', 'B2: the rank tuple carries the Phi value instead of the PLCP value', 'new RankTuple(r, phi[k])', 'new RankTuple(r, rank[k])')
add(b, 'B10', 'B3: rank chunk 0 expects a tuple for the empty suffix', 'long expected = c == 0 ? length - 1 : length;', 'long expected = length;')
add(b, 'B11', 'B3: LCP[0] is whatever the buffer holds', 'if (c == 0) dense[0] = 0; // LCP[0] = 0', '')
add(b, 'B12', 'B3: the consumer is not called', '_consumer?.Consume(firstRank, dense);', '')
add(b, 'B13', 'the round count of the bounded mode is floored', 'long rounds = _g.TextChunks == 0 ? 0 : (_g.TextChunks + perRound - 1) / perRound;', 'long rounds = _g.TextChunks == 0 ? 0 : _g.TextChunks / perRound;')
add(b, 'B14', 'bounded B1: the chunk range of a round includes the next chunk', 'if (chunk >= first && chunk < last) writer.Add', 'if (chunk >= first && chunk <= last) writer.Add')
add(b, 'B15', 'bounded B2: tasks of a later round start at chunk 0', 'long c = Volatile.Read(ref _roundFirst) + Interlocked.Increment(ref _nextChunk) - 1;', 'long c = Interlocked.Increment(ref _nextChunk) - 1;')
add(b, 'B16', 'bounded B2: the PLCP file bytes are not tracked as scratch', 'ScratchTracker.Add(length * 4L);', '')
add(b, 'B17', 'bounded: the PLCP file is left behind', 'File.Delete(plcpPath);\n            ScratchTracker.Remove(_plcpFileBytes);', 'ScratchTracker.Remove(_plcpFileBytes);')
add(b, 'B18', 'bounded: the tracker keeps the PLCP file', 'ScratchTracker.Remove(_plcpFileBytes);', '')
add(b, 'B19', 'bounded B3: the suffix array block is re-read from the start of the chunk', 'sa.ReadExactAt(MemoryMarshal.AsBytes(read.AsSpan(take)), (start + done) * 4);', 'sa.ReadExactAt(MemoryMarshal.AsBytes(read.AsSpan(take)), start * 4);')
add(b, 'B20', 'bounded B3: the gather reads the neighbour', 'dense[done + i] = plcp[s];', 'dense[done + i] = plcp[s == 0 ? 0 : s - 1];')
add(b, 'B21', 'the width is chosen from swapped counts', 'LcpWidth.Select(n, _escapes8, _escapes16)', 'LcpWidth.Select(n, _escapes16, _escapes8)')
add(b, 'B22', 'the task count may exceed the number of chunks by one', 'workers = (int)Math.Max(1, Math.Min(workers, chunks));', 'workers = (int)Math.Max(1, Math.Min(workers, chunks + 1));')
add(b, 'B23', 'the task count ignores the budget', 'long fit = (_budget.CapBytes - _budget.LiveBytes) / Math.Max(1, perTaskBytes);', 'long fit = int.MaxValue;')
add(b, 'B24', 'Merge: the escape counts of a task replace the totals', '_escapes8 += t.Escapes8;', '_escapes8 = t.Escapes8;')
add(b, 'B25', 'Merge: the largest PLCP of the last task wins', '_maxPlcp = Math.Max(_maxPlcp, t.Max);', '_maxPlcp = t.Max;')
add(b, 'B26', 'Merge: the restart cost is not summed', '_restart += t.Restart;', '')
add(b, 'B27', 'the base mode reports no round', '_stats.Rounds = 1;', '_stats.Rounds = 0;')
add(b, 'B28', 'WorkerBytes: Phi and rank chunk counted once', 'return (8 * chunk) + (Math.Max(1, o.ReadBufferBytes / sizeof(PhiTuple))', 'return (4 * chunk) + (Math.Max(1, o.ReadBufferBytes / sizeof(PhiTuple))')
add(b, 'B29', 'a failed worker does not stop the others (error only read at the end)', 'while (Volatile.Read(ref _error) is null)\n        {\n            long c = Interlocked.Increment(ref _nextChunk) - 1;\n            if (c >= _g.TextChunks) break;\n            int length = (int)_g.TextChunkLength(c);\n            long start = c * _g.R;\n            using (SpoolReader<PhiTuple>', 'while (true)\n        {\n            long c = Interlocked.Increment(ref _nextChunk) - 1;\n            if (c >= _g.TextChunks) break;\n            int length = (int)_g.TextChunkLength(c);\n            long start = c * _g.R;\n            using (SpoolReader<PhiTuple>')
add(b, 'B30', 'B2: the rank tuples are not flushed', 'if (Volatile.Read(ref _error) is null) writer.FlushAll();', '')
add(b, 'B31', 'B1: the Phi tuples are not flushed', '        writer.FlushAll();\n    }\n\n    /// <summary>B2 of one task: chunks in any order;', '    }\n\n    /// <summary>B2 of one task: chunks in any order;')
add(b, 'B32', 'the stats keep the escape counts of the last chunk only (no sum in the task)', 'Escapes8 += r.Escapes8;', 'Escapes8 = r.Escapes8;')
add(b, 'B33', 'the overflow count comes from the wrong tier', '_stats.OverflowCount = LcpWidth.OverflowCount(_stats.Kind, _escapes8, _escapes16);', '_stats.OverflowCount = LcpWidth.OverflowCount(LcpKind.U8, _escapes8, _escapes16);')

# --- ExternalLcp.cs (nothing but options and stats; the entry point)
# --- E4.cs
e = APP + 'Experiments/E4.cs'
add(e, 'E01', 'R scaled down rounds towards zero', 'Math.Max(3, (long)Math.Ceiling(ChunkRef * f))', 'Math.Max(3, (long)Math.Floor(ChunkRef * f))')
add(e, 'E02', 'rounds: floor of P / r chunks per round', '(chunks + long.Parse(rounds, CultureInfo.InvariantCulture) - 1) / long.Parse(rounds, CultureInfo.InvariantCulture)', 'chunks / long.Parse(rounds, CultureInfo.InvariantCulture)')
add(e, 'E03', 'B2 workers: twice what the reference geometry allows', 'budgetCap / Math.Max(1, perTask)', 'budgetCap / Math.Max(1, perTask / 2)')
add(e, 'E04', 'the exit code ignores the oracle', 'ok &= equal;', '')
add(e, 'E05', 'the exit code ignores leftover scratch', 'left.Length == 0 && scratchLive == 0', 'scratchLive == 0')
add(e, 'E06', 'check-lcp: structure scan accepts a missing pair', 'if (width == 4 ? overflowCount != 0 : escapesInBase != overflowCount) problems++;', '')
add(e, 'E07', 'check-lcp: samples are not compared', 'if (got != expected) mismatches++;', '')
add(e, 'E08', 'check-lcp: the rank of a pair is not compared', 'if (pr != rank || pr <= previousRank || pv < escape) problems++;', 'if (pr <= previousRank || pv < escape) problems++;')
add(e, 'E09', 'check-lcp: the value of a pair is not required to reach the escape', 'if (pr != rank || pr <= previousRank || pv < escape) problems++;', 'if (pr != rank || pr <= previousRank) problems++;')
add(e, 'E10', 'check-lcp: the last rank is not sampled', 'if (s == 1) rank = n;', '')
add(e, 'E11', 'lcp-oracle: the overflow ranks are written as values', 'pairs[2 * i] = (uint)lcp.OverflowRanks[i];', 'pairs[2 * i] = (uint)lcp.OverflowValues[i];')

for f, items in M.items():
    src = open(os.path.join(ROOT, f)).read()
    for it in items:
        c = src.count(it['old'])
        if c != 1:
            print(f"PATTERN {it['id']} occurs {c} times in {f}")
            sys.exit(1)
        if it['old'] == it['new']:
            print(f"PATTERN {it['id']} is not a change")
            sys.exit(1)
    name = {'LcpGeometry': 'geometry', 'Spool': 'spool', 'PlcpScan': 'plcpscan', 'LcpOutput': 'output', 'LcpBuild': 'build', 'E4': 'e4'}[os.path.basename(f)[:-3]]
    json.dump(items, open(os.path.join(OUT, f'e4_{name}_mutants.json'), 'w'), indent=1)
print(sum(len(v) for v in M.values()), 'mutants in', len(M), 'files')
