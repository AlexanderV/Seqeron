#!/usr/bin/env python3
"""Writes the hand-made mutant lists of the E3b audit (one JSON per production file; each mutant is a single-site edit: old text -> new text, `old` occurring exactly once in the file)
and checks that every pattern is unique in the current source. usage: make_e3b_mutants.py <repo root>"""
import json, os, sys

ROOT = sys.argv[1] if len(sys.argv) > 1 else '/home/user/Seqeron'
APP = 'apps/SuffixTree.EmProto/'
OUT = os.path.join(ROOT, 'docs/refactoring/suffix-tree-upgrade/analysis/mutation')

M = {}

def add(file, id_, desc, old, new):
    M.setdefault(file, []).append({'id': id_, 'desc': desc, 'old': old, 'new': new})

# --- SampleGeometry.cs
g = APP + 'Dcs/SampleGeometry.cs'
add(g, 'G01', 'Slots: one block too few', 'Slots = ((n / Period) + 1) * Classes;', 'Slots = (n / Period) * Classes;')
add(g, 'G02', 'ClassCount: a residue equal to n has no member', 'ClassCount[c] = r <= n ? ((n - r) / Period) + 1 : 0;', 'ClassCount[c] = r < n ? ((n - r) / Period) + 1 : 0;')
add(g, 'G03', 'RedStart: separator of the class not counted', 'RedStart[c + 1] = ClassStart[c + 1] + c + 1;', 'RedStart[c + 1] = ClassStart[c + 1] + c;')
add(g, 'G04', 'SlotOf: block stride is the period instead of |D|', 'public long SlotOf(long p) => (p / Period * Classes)'.replace('public ', 'internal ') , 'internal long SlotOf(long p) => (p / Period * Period)')
add(g, 'G05', 'TrySlotOfRed: the separator index k == count is taken as a member', 'if (k >= ClassCount[lo])', 'if (k > ClassCount[lo])')
add(g, 'G06', 'RedIndexOfSlot: block taken as slot modulo |D|', 'RedStart[(int)(slot % Classes)] + (slot / Classes);', 'RedStart[(int)(slot % Classes)] + (slot % Classes);')
add(g, 'G07', 'TrySlotOfRed: binary search keeps the class whose start is strictly below q', 'if (RedStart[mid] <= q) lo = mid;', 'if (RedStart[mid] < q) lo = mid;')
add(g, 'G08', 'PositionOfSlot: residue taken from the block index', 'Cover.Residues[(int)(slot % Classes)];', 'Cover.Residues[(int)(slot / Classes % Classes)];')

# --- MappedWindowWriter.cs
w = APP + 'Dcs/MappedWindowWriter.cs'
add(w, 'W01', 'WriteAt: a write at the end of the window does not open the next one (loops for ever)', 'offset < _viewStart || offset >= _viewEnd', 'offset < _viewStart || offset > _viewEnd')
add(w, 'W02', 'OpenWindow: window not clipped to the file end', 'long size = Math.Min(_window, Length - offset);', 'long size = _window;')
add(w, 'W03', 'Append: position not advanced', '_next += data.Length;', '_next += 0;')
add(w, 'W04', 'WriteAt: the last byte of the file cannot be written', 'offset + data.Length > Length)', 'offset + data.Length >= Length)')
add(w, 'W05', 'WriteAt: mapped write counter not updated', 'MappedIo.AddWrite(take);', '')
add(w, 'W06', 'OpenWindow: window counter not updated', 'MappedIo.AddWindow(size);', '')

# --- Runs.cs
r = APP + 'Dcs/Runs.cs'
add(r, 'R01', 'Create: a run of exactly the spill size goes to a file', 'if ((long)sorted.Length * 8 <= spillBytes)', 'if ((long)sorted.Length * 8 < spillBytes)')
add(r, 'R02', 'ReadAt: reads at least the buffer (Max instead of Min)', 'int take = (int)Math.Min(dest.Length, left);', 'int take = (int)Math.Max(dest.Length, left);')
add(r, 'R03', 'Drain: a record at the limit is scattered too', 'if (slot >= limit) break;', 'if (slot > limit) break;')
add(r, 'R04', 'Drain: the bucket offset is not added', 'window[(int)(slot - windowStart)] = offset + (uint)v;', 'window[(int)(slot - windowStart)] = (uint)v;')
add(r, 'R05', 'Refill: a drained run is refilled (loops for ever)', 'if (_consumed >= run.Count) return false;', 'if (_consumed > run.Count) return false;')
add(r, 'R06', 'Create: the tracker is not told about the file', 'ScratchTracker.Add((long)sorted.Length * 8);', '')
add(r, 'R07', 'Dispose: the tracker is not told about the deletion', 'ScratchTracker.Remove(Count * 8);', '')

# --- SampleNamer.cs
n = APP + 'Dcs/SampleNamer.cs'
add(n, 'N01', 'PrefixOrder: equal X-prefix needs l > X', 'if (l >= x) return 0;', 'if (l > x) return 0;')
add(n, 'N02', 'PrefixOrder: the shorter suffix is the larger', 'if (l == maxLen) return pa > pb ? -1 : 1;', 'if (l == maxLen) return pa > pb ? 1 : -1;')
add(n, 'N03', 'AssignNames: a group is continued over the boundary when the prefixes differ', '_prefix.Compare(prev.LastPos, leaf.FirstPos) == 0;', '_prefix.Compare(prev.LastPos, leaf.FirstPos) != 0;')
add(n, 'N04', 'AssignNames: continued bucket counts its first group again', 'names += leaf.Distinct - 1;', 'names += leaf.Distinct;')
add(n, 'N05', 'AssignNames: the offset of a continued bucket', 'leaf.Offset = (uint)(names - 1);', 'leaf.Offset = (uint)names;')
add(n, 'N06', 'AssignNames: singleton group at a boundary counted at 2', 'names += leaf.Distinct;\n                if (openSize == 1) singles++;', 'names += leaf.Distinct;\n                if (openSize == 2) singles++;')
add(n, 'N07', 'AssignNames: non-unique records are the singles', 'r.NonUnique = _g.Samples - singles;', 'r.NonUnique = singles;')
add(n, 'N08', 'AssignNames: the last group of a bucket is not discounted from the inner singles', 'if (leaf.LastLen == 1) inner--;', '')
add(n, 'N09', 'InRange: the lower bound record itself is excluded', 'if (key == lo.Key && Cmp(lo.Pos, p) > 0) return false;', 'if (key == lo.Key && Cmp(lo.Pos, p) >= 0) return false;')
add(n, 'N10', 'InRange: the upper bound record is included', 'if (key == hi.Key && Cmp(p, hi.Pos) >= 0) return false;', 'if (key == hi.Key && Cmp(p, hi.Pos) > 0) return false;')
add(n, 'N11', 'Cmp: no position tie-break', 'return c != 0 ? c : a.CompareTo(b);\n        }\n\n        /// <summary>One scan that stores', 'return c;\n        }\n\n        /// <summary>One scan that stores')
add(n, 'N12', 'Collect: overflow detected one record late', 'if (stored == cap)', 'if (stored > cap)')
add(n, 'N13', 'Split: fan-out from count instead of 2 x count (output unchanged)', '(2 * count + bucket - 1) / bucket', '(count + bucket - 1) / bucket')
add(n, 'N14', 'Split: equal splitter indices kept (empty children only)', 'if (i <= last) continue;', 'if (i < last) continue;')
add(n, 'N15', 'Split: reservoir replacement may hit the slot after the last', 'if (j < (ulong)cap) reservoir[(int)j] = (uint)p;', 'if (j <= (ulong)cap) reservoir[(int)j] = (uint)p;')
add(n, 'N16', 'Split: the last reservoir record is not sorted', 'Array.Sort(reservoir, 0, size, new RecordComparer(text, owner._prefix));', 'Array.Sort(reservoir, 0, size - 1, new RecordComparer(text, owner._prefix));')
add(n, 'N17', 'Finish: insertion sort unstable on equal keys (output unchanged)', 'while (j >= 0 && keys[j] > k)', 'while (j >= 0 && keys[j] >= k)')
add(n, 'N18', 'RadixSort: wrong digit in the histogram', 'hist[256 + (int)((k >> 8) & 0xFF)]++;', 'hist[256 + (int)((k >> 9) & 0xFF)]++;')
add(n, 'N19', 'Finish: a new group starts when the keys are equal', 'if (i == 0 || keys[i] != keys[i - 1] || owner._prefix.Compare(pos[i - 1], pos[i]) != 0)', 'if (i == 0 || keys[i] == keys[i - 1] || owner._prefix.Compare(pos[i - 1], pos[i]) != 0)')
add(n, 'N20', 'Finish: slot shifted by 31 bits', 'values[i] = ((ulong)owner._g.SlotOf(pos[i]) << 32) | (uint)distinct;', 'values[i] = ((ulong)owner._g.SlotOf(pos[i]) << 31) | (uint)distinct;')
add(n, 'N21', 'Finish: the last group of the bucket not counted as a singleton', 'if (lastLen == 1) singles++;', '')
add(n, 'N22', 'Finish: first group length of a one-group bucket not set (unused)', 'if (distinct == 1) firstLen = lastLen;', '')
add(n, 'N23', 'Merge: fair share multiplied instead of divided', 'long fair = _o.MergeBufferBudget / 8 / Math.Max(1, fileRuns);', 'long fair = _o.MergeBufferBudget / 8 * Math.Max(1, fileRuns);')
add(n, 'N24', 'Merge: the window check only fails on too many records (defensive)', 'if (got != end - start) throw', 'if (got > end - start) throw')
add(n, 'N25', 'Collect: buffers sized by the capacity, not by the sample', 'int cap = (int)Math.Min(owner._o.BucketCapacity, Math.Max(1, owner._g.Samples));', 'int cap = (int)owner._o.BucketCapacity;')
add(n, 'N26', 'Worker: file runs not counted', '_fileRuns += w.FileRuns;', '')
add(n, 'N27', 'Run: the root is split only when the sample exceeds twice the capacity', '_g.Samples > _o.BucketCapacity, 1)', '_g.Samples > 2 * _o.BucketCapacity, 1)')
add(n, 'N28', 'Merge: window size taken as the capacity of a bucket instead of the chunk', 'int windowSlots = (int)Math.Min(Math.Max(1, _o.ChunkSlots), s);', 'int windowSlots = (int)Math.Min(Math.Max(1, _o.BucketCapacity), s);')

# --- NameSinks.cs
s = APP + 'Dcs/NameSinks.cs'
add(s, 'S01', 'RankFileSink: slot offset in 8-byte units', 'public void Accept(long firstSlot, ReadOnlySpan<uint> names) => _writer.WriteAt(firstSlot * 4, MemoryMarshal.AsBytes(names));', 'public void Accept(long firstSlot, ReadOnlySpan<uint> names) => _writer.WriteAt(firstSlot * 8, MemoryMarshal.AsBytes(names));')
add(s, 'S02', 'ReducedFileSink: little endian', 'BinaryPrimitives.WriteUInt32BigEndian(buf[(_filled[c] * 4)..], names[i]);', 'BinaryPrimitives.WriteUInt32LittleEndian(buf[(_filled[c] * 4)..], names[i]);')
add(s, 'S03', 'ReducedFileSink: flush offset ignores what was flushed', '_file.WriteAt(buf[..(_filled[c] * 4)], (_g.RedStart[c] + _flushed[c]) * 4);', '_file.WriteAt(buf[..(_filled[c] * 4)], _g.RedStart[c] * 4);')
add(s, 'S04', 'ReducedFileSink: no padding after the symbols', 'long bytes = (geometry.RedLength * 4) + PackedText.PadBytes;', 'long bytes = geometry.RedLength * 4;')
add(s, 'S05', 'ReducedFileSink: a full cursor buffer is flushed one entry late', 'if (++_filled[c] == _bufEntries) Flush(c, buf);', 'if (++_filled[c] > _bufEntries) Flush(c, buf);')
add(s, 'S06', 'ReducedRamSink: class offset replaced by the class index', 'red[(int)(_g.RedStart[c] + block)] = (int)names[i];', 'red[(int)(c + block)] = (int)names[i];')
add(s, 'S07', 'ReducedRamSink: block not advanced at the end of a row', 'block++;\n            }', 'block += 0;\n            }')
add(s, 'S08', 'RankFileSink: the scratch output is not tracked', 'if (scratch) ScratchTracker.Add(slots * 4);', '')

# --- SaInverter.cs
i = APP + 'Dcs/SaInverter.cs'
add(i, 'I01', 'Invert: separators are treated as real positions', 'if (!g.TrySlotOfRed(entries[i], out long slot)) continue; // a separator', 'if (g.TrySlotOfRed(entries[i], out long slot)) continue; // a separator')
add(i, 'I02', 'Invert: the entry of rank 1 is skipped instead of rank 0', 'if (rank == 0) continue; // the empty suffix', 'if (rank == 1) continue; // the empty suffix')
add(i, 'I03', 'Invert: chunk by modulo', 'long chunk = slot / chunkSlots;', 'long chunk = slot % chunkSlots;')
add(i, 'I04', 'Invert: chunk array not cleared', 'array.Clear();', '')
add(i, 'I05', 'Invert: a full tuple buffer is flushed one entry late', 'if (filled[chunk] == bufEntries) Flush(chunk, buf, filled, flushed, o, tag);', 'if (filled[chunk] > bufEntries) Flush(chunk, buf, filled, flushed, o, tag);')
add(i, 'I06', 'Invert: tuple files not untracked', 'ScratchTracker.Remove(flushed[chunk] * 8);', '')
add(i, 'I07', 'Invert: tuple file not deleted', 'File.Delete(path);', '')
add(i, 'I08', 'Flush: tuple file written from the start every time', 'f.WriteAt(MemoryMarshal.AsBytes(buf[..filled[chunk]]), flushed[chunk] * 8);', 'f.WriteAt(MemoryMarshal.AsBytes(buf[..filled[chunk]]), 0);')

# --- ExternalDcs.cs
d = APP + 'Dcs/ExternalDcs.cs'
add(d, 'D01', 'BuildLevel: a reduced string of exactly N_base symbols recurses', 'if (reduced <= o.BaseCase || reduced >= text.Length)', 'if (reduced < o.BaseCase || reduced >= text.Length)')
add(d, 'D02', 'BuildLevel: a reduced string as long as the text recurses', 'if (reduced <= o.BaseCase || reduced >= text.Length)', 'if (reduced <= o.BaseCase || reduced > text.Length)')
add(d, 'D03', 'BuildLevel: level rank file not deleted', 'File.Delete(childRank);', '')
add(d, 'D04', 'BuildLevel: suffix array file not untracked', 'ScratchTracker.Remove((reduced + 1) * 4);\n        ls.BelowSeconds', 'ls.BelowSeconds')
add(d, 'D05', 'BaseCase: the inverse starts at rank 0', 'for (int k = 1; k <= n; k++) rank[saSpan[k]] = k;', 'for (int k = 0; k <= n; k++) rank[saSpan[k]] = k;')
add(d, 'D06', 'BaseCase: unused slots read the table', 'window[i] = slot < g.Samples ? (uint)rank[(int)g.RedIndexOfSlot(slot)] : 0;', 'window[i] = slot <= g.Samples ? (uint)rank[(int)g.RedIndexOfSlot(slot)] : 0;')
add(d, 'D07', 'Build: sigma of the base case one too small', 'BaseCase(sink.Red, g, ls.DistinctNames + 1, rankPath', 'BaseCase(sink.Red, g, ls.DistinctNames, rankPath')
add(d, 'D08', 'Validation: bucket capacity 1 accepted', 'ThrowIfLessThan(options.BucketCapacity, 2)', 'ThrowIfLessThan(options.BucketCapacity, 1)')
add(d, 'D09', 'Validation: reservoir of 2 accepted', 'ThrowIfLessThan(options.ReservoirSize, 3)', 'ThrowIfLessThan(options.ReservoirSize, 2)')
add(d, 'D10', 'BuildLevel: the unique branch does not return after the merge', 'ls.MergeSeconds = Seconds(t0);\n            return;\n        }\n\n        long reduced', 'ls.MergeSeconds = Seconds(t0);\n        }\n\n        long reduced')
add(d, 'D11', 'BuildLevel: level counter of the child not advanced', 'BuildLevel(red, cover, level + 1, childRank, true, o, tag, budget, stats);', 'BuildLevel(red, cover, level, childRank, true, o, tag, budget, stats);')
add(d, 'D12', 'BuildLevel: external suffix sort seeded and sized by the naming capacity', 'LeafCapacity = o.SortLeafCapacity,', 'LeafCapacity = o.BucketCapacity,')

# --- E3b.cs
e = APP + 'Experiments/E3b.cs'
add(e, 'C01', 'Scale: scaled values rounded down to the floor', 'return Math.Max(floor, (long)Math.Round(reference * f));', 'return Math.Min(floor, (long)Math.Round(reference * f));')
add(e, 'C02', 'CompareFiles: a length mismatch is not a difference', 'if (gx != gy) return (false, at + common);', 'if (gx != gy) return (true, -1);')
add(e, 'C03', 'PackRandom: unused bits of the last byte not cleared', 'buf[take - 1] &= (byte)(0xFF << (8 - (int)((n * 2) % 8)));', '')
add(e, 'C04', 'Build: the verify result does not change the exit code', 'ok &= equal;', '')
add(e, 'C05', 'CompareFiles: first difference reported one early', 'if (diff < common) return (false, at + diff);', 'if (diff < common) return (false, at + diff - 1);')

# --- ArenaBudget.cs
a = 'apps/SuffixTree.EmProto/Harness/ArenaBudget.cs'
add(a, 'T01', 'TouchPages: every second page', 'for (long i = 0; i < bytes; i += 4096)', 'for (long i = 0; i < bytes; i += 8192)')
add(a, 'T02', 'Allocate: pages touched when not asked', 'if (_touch) TouchPages', 'if (!_touch) TouchPages')

# --- round 2: the survivors of round 1 that are not equivalent (against the tests added after round 1) and two mutants that did not compile in round 1, in a compiling form
R2 = {}

def add2(file, id_, desc, old, new):
    R2.setdefault(file, []).append({'id': id_, 'desc': desc, 'old': old, 'new': new})

for f, lst in M.items():
    for m in lst:
        if m['id'] in ('N14', 'N24', 'N27', 'N28', 'D12', 'C04'):
            add2(f, m['id'], m['desc'], m['old'], m['new'])
add2(n, 'N26b', 'Worker: file runs counted one too many per worker', '_fileRuns += w.FileRuns;', '_fileRuns += w.FileRuns + 1;')
add2(e, 'C03b', 'PackRandom: unused bits of the last byte not cleared (whole statement removed)', 'if (written + take == dataBytes && (n * 2) % 8 != 0) buf[take - 1] &= (byte)(0xFF << (8 - (int)((n * 2) % 8)));', '')

bad = 0
for f, lst in R2.items():
    src = open(os.path.join(ROOT, f)).read()
    for m in lst:
        if src.count(m['old']) != 1:
            print('NOT UNIQUE (round 2)', f, m['id'], src.count(m['old']))
            bad += 1
    name = os.path.splitext(os.path.basename(f))[0].lower()
    json.dump(lst, open(os.path.join(OUT, f'e3b_round2_{name}_mutants.json'), 'w'), indent=1)
for f, lst in M.items():
    src = open(os.path.join(ROOT, f)).read()
    for m in lst:
        if src.count(m['old']) != 1:
            print('NOT UNIQUE', f, m['id'], src.count(m['old']))
            bad += 1
    name = os.path.splitext(os.path.basename(f))[0].lower()
    json.dump(lst, open(os.path.join(OUT, f'e3b_{name}_mutants.json'), 'w'), indent=1)
total = sum(len(v) for v in M.values())
print(total, 'mutants,', bad, 'patterns not unique')
