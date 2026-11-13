#!/usr/bin/env python3
"""Hand-made single-site mutants of the new production code of E5 (PH5). Writes e5_<name>_mutants.json (id, desc, old, new[, also]) next to this script and checks that every `old` pattern is unique in its file.
usage: make_e5_mutants.py [repo root]"""
import json, os, sys

ROOT = sys.argv[1] if len(sys.argv) > 1 else '/home/user/Seqeron'
HERE = os.path.dirname(os.path.abspath(__file__))
APP = 'apps/SuffixTree.EmProto/'

M = {}
M['spillstack'] = (APP + 'Ph5/SpillStack.cs', [
    ('S01', 'the stack advances one entry early (block of B - 1 entries)', 'if (_count == _block) Advance();', 'if (_count == _block - 1) Advance();'),
    ('S02', 'Pop retreats also at height 0', 'if (--_count == 0 && _height > 0) Retreat();', 'if (--_count == 0) Retreat();'),
    ('S03', 'every evicted block is written, also a clean one', 'if (_otherDirty) WriteBlock(_topBlock - 1, _other);', 'if (true) WriteBlock(_topBlock - 1, _other);'),
    ('S04', 'the block below is always dirty after an advance', '_otherDirty = _curDirty;\n        _curDirty = false;\n        _otherBelow = true;', '_otherDirty = true;\n        _curDirty = false;\n        _otherBelow = true;'),
    ('S05', 'a block that was read back is dirty', '            ReadBlock(_topBlock - 1, _other);\n            _otherDirty = false;', '            ReadBlock(_topBlock - 1, _other);\n            _otherDirty = true;'),
    ('S06', 'the new top block after an advance starts dirty (the push that follows sets it anyway)', '_curDirty = false;\n        _otherBelow = true;', '_curDirty = true;\n        _otherBelow = true;'),
    ('S07', 'a refilled block is one entry short', '        _topBlock--;\n        _count = _block;', '        _topBlock--;\n        _count = _block - 1;'),
    ('S08', 'Retreat reads when the block below is resident and not otherwise', 'if (!_otherBelow)\n        {\n            ReadBlock', 'if (_otherBelow)\n        {\n            ReadBlock'),
    ('S09', 'a block is read from the next block\'s offset', '_file!.ReadExactAt(new Span<byte>(buffer, _blockBytes), block * _blockBytes);', '_file!.ReadExactAt(new Span<byte>(buffer, _blockBytes), (block + 1) * _blockBytes);'),
    ('S10', 'the unflushed byte count does not grow (no fsync)', '_unflushed += _blockBytes;', '_unflushed += 0;'),
    ('S11', 'Dispose keeps the spill file', '        File.Delete(_path);\n        ScratchTracker.Remove(_fileBytes);', '        ScratchTracker.Remove(_fileBytes);'),
    ('S12', 'the peak height lags one push', 'if (++_height > _peak) _peak = _height;', 'if (_height++ > _peak) _peak = _height;'),
    ('S13', 'an update through TopForUpdate does not mark the block dirty', '        _curDirty = true;\n        return ref _cur[_count - 1];', '        return ref _cur[_count - 1];'),
    ('S14', 'a block size of 1 is accepted', 'ArgumentOutOfRangeException.ThrowIfLessThan(blockEntries, 2);', 'ArgumentOutOfRangeException.ThrowIfLessThan(blockEntries, 1);'),
    ('S15', 'Advance writes only when the other buffer is not below', 'if (_otherBelow)\n        {\n            // the other buffer holds', 'if (!_otherBelow)\n        {\n            // the other buffer holds'),
    ('S16', 'blocks read are not counted', '        _blocksRead++;\n', ''),
    ('S17', 'bytes written are not counted', '        _bytesWritten += _blockBytes;\n', ''),
    ('S18', 'a refilled block does not become the top (the buffers are not swapped)', '        T* top = _other;\n        _other = _cur;\n        _cur = top;\n        _curDirty = _otherDirty;', '        _curDirty = _otherDirty;'),
    ('S19', 'the file high-water mark is not tracked', '            _fileBytes = end;\n', ''),
    ('S20', 'the empty stack does not refuse Top', 'if (_height == 0) ThrowEmpty();\n            return ref _cur[_count - 1];', 'return ref _cur[_count - 1];'),
])
M['windowwriter'] = (APP + 'Ph5/AlignedWindowWriter.cs', [
    ('W01', 'a window that is not a power of two is accepted', 'windowBytes < 4 || !BitOperations.IsPow2(windowBytes)', 'windowBytes < 4'),
    ('W02', 'the bound of Write32 is one entry too wide', '(ulong)index >= (ulong)(Length >> 2)', '(ulong)index > (ulong)(Length >> 2)'),
    ('W03', 'the window that is released is the most recently used', '_slots[i].LastUse < oldest', '_slots[i].LastUse > oldest'),
    ('W04', 'a window is flushed also when flushing is off', 'if (_flush) FlushSlot(ref _slots[victim]);', 'FlushSlot(ref _slots[victim]);'),
    ('W05', 'a window mapped again is counted as new', 'if ((_seen[w >> 6] & bit) != 0)', 'if ((_seen[w >> 6] & bit) == 0)'),
    ('W06', 'Dispose does not flush', '            if (_flush) FlushSlot(ref _slots[i]);\n            Release(ref _slots[i]);', '            Release(ref _slots[i]);'),
    ('W07', 'Write does not cut the piece at the window end', 'int take = (int)Math.Min(data.Length, _slots[HotSlot()].Size - inside);', 'int take = data.Length;'),
    ('W08', 'the last window is mapped at the full window size', 'long size = Math.Min(_window, Length - start);', 'long size = _window;'),
    ('W09', 'evictions are not counted', '            _evictions++;\n', ''),
    ('W10', 'a hit does not make the window the most recent', '                _slots[i].LastUse = ++_clock;\n                _hotWindow = w;', '                _hotWindow = w;'),
    ('W11', 'Flush() flushes nothing', '            if (_slots[i].View is not null) FlushSlot(ref _slots[i]);\n        }\n    }', '            if (_slots[i].View is null) FlushSlot(ref _slots[i]);\n        }\n    }'),
    ('W12', 'a flush is not counted in bytes', '        _flushedBytes += slot.Size;', '        _flushedBytes += 0;'),
    ('W13', 'writes are not counted', '        _writes++;\n', ''),
    ('W14', 'distinct windows are not counted', '            _distinct++;\n', ''),
])
M['compactlcpfile'] = (APP + 'Ph5/CompactLcpFile.cs', [
    ('C01', 'a 4-byte base array is refused', 'Width is not (1 or 2 or 4)', 'Width is not (1 or 2)'),
    ('C02', 'a pair is accepted when only one of rank and value is wrong (forward)', 'if (rank != first + i || value < _escape) throw new InvalidDataException($"overflow pair {pair - 1}', 'if (rank != first + i && value < _escape) throw new InvalidDataException($"overflow pair {pair - 1}'),
    ('C03', 'extra pairs are not detected', 'if (pair != OverflowCount) throw', 'if (pair > OverflowCount) throw'),
    ('C04', 'the backward cursor starts one pair too far', 'long pair = OverflowCount - 1;', 'long pair = OverflowCount;'),
    ('C05', 'the backward end check expects a pair left', 'if (pair != -1) throw', 'if (pair != 0) throw'),
    ('C06', 'the last block is one entry short', 'int count = (int)Math.Min(_blockEntries, _n + 1 - first);', 'int count = (int)Math.Min(_blockEntries, _n - first);'),
    ('C07', 'the number of pairs is a quarter of the file', 'OverflowCount = overflowLength / 8;', 'OverflowCount = overflowLength / 4;'),
    ('C08', 'the overflow file is accepted when it holds half pairs', 'if (overflowLength % 8 != 0)', 'if (overflowLength % 4 != 0)'),
    ('C09', 'the overflow window going backwards starts at the index', 'long start = backward ? Math.Max(0, index - capacity + 1) : index;', 'long start = index;'),
    ('C10', 'a value below the escape is accepted in a pair (backward)', 'if (rank != first + i || value < _escape) throw new InvalidDataException($"overflow pair {pair + 1}', 'if (rank != first + i) throw new InvalidDataException($"overflow pair {pair + 1}'),
    ('C11', 'the base length check accepts any multiple', 'baseLength != (n + 1) * Width', 'baseLength < (n + 1) * Width'),
    ('C12', 'U16 entries are read as bytes', '                ReadOnlySpan<ushort> u16 = MemoryMarshal.Cast<byte, ushort>(bytes);\n                for (int i = 0; i < count; i++) _values[i] = u16[i];', '                for (int i = 0; i < count; i++) _values[i] = bytes[i];'),
    ('C13', 'a 4-byte base array with an overflow list is accepted', 'if (Width == 4 && OverflowCount != 0) throw', 'if (Width == 4 && OverflowCount < 0) throw'),
])
M['childscan'] = (APP + 'Ph5/ChildScan.cs', [
    ('A01', 'LCP[1] is not checked', 'if (rank <= 1)', 'if (rank == 0)'),
    ('A02', 'a value of n + 1 is accepted', 'else if (cur > _n)', 'else if (cur > _n + 1)'),
    ('A03', 'the last-child flag compares with <=', 'Close(e.Ell, e.Lb, i - 1, e.First, cur < top.Ell);', 'Close(e.Ell, e.Lb, i - 1, e.First, cur <= top.Ell);'),
    ('A04', 'the left boundary of a new interval is not the one of the last popped interval', '            left = e.Lb;\n', ''),
    ('A05', 'the PSV of a position equal to the root is the root boundary', 'psv = top.Ell == 0 ? ChildScan.None : top.Lb;', 'psv = top.Lb;'),
    ('A06', 'an equal value pushes a new interval', 'if (top.Ell < cur)\n        {\n            _stack.Push', 'if (top.Ell <= cur)\n        {\n            _stack.Push'),
    ('A07', 'the first l-index of the root is not recorded separately', 'else if (top.First == ChildScan.None)', 'else if (top.First == ChildScan.None && false)'),
    ('A08', 'the root counts as an internal node', '        if (ell > 0)\n        {\n            _internal++;', '        if (ell >= 0)\n        {\n            _internal++;'),
    ('A09', 'the deepest interval is the last of equal depth', 'if (ell > _lrsDepth)', 'if (ell >= _lrsDepth)'),
    ('A10', 'NODE_COUNT is one too high', 'NodeCount = _n == 0 ? 1 : 1 + _internal + _n + 1,', 'NodeCount = _n == 0 ? 1 : 2 + _internal + _n + 1,'),
    ('A11', 'the remaining intervals are closed as non-last children', 'Close(e.Ell, e.Lb, (uint)_n, e.First, true);', 'Close(e.Ell, e.Lb, (uint)_n, e.First, false);'),
    ('A12', 'chunks may arrive out of order', 'if (firstRank != _nextRank) throw', 'if (firstRank < 0) throw'),
    ('A13', 'the last chunk is refused', 'if (firstRank + lcp.Length > _n + 1) throw', 'if (firstRank + lcp.Length >= _n + 1) throw'),
    ('A14', 'Complete does not need every rank', 'if (_nextRank != _n + 1) throw', 'if (_nextRank > _n + 1) throw'),
    ('A15', 'a text of 2^32 - 2 symbols is accepted', 'if (n >= uint.MaxValue - 1) throw', 'if (n > uint.MaxValue - 1) throw'),
    ('A16', 'the tagged codec accepts one symbol too many', 'options.Codec == ChildCodec.Tagged && n > TaggedChild.MaxLength', 'options.Codec == ChildCodec.Tagged && n > TaggedChild.MaxLength + 1'),
    ('A17', 'an up field is written with the down tag', 'else writer.Write32(rb, UpTag | first);', 'else writer.Write32(rb, DownTag | first);'),
    ('A18', 'the packed codec writes every next link', 'if (previous != ChildScan.None && k != previous + 1) writer.Write32(previous, k);', 'if (previous != ChildScan.None) writer.Write32(previous, k);'),
    ('A19', 'rank 0 gets no PSV entry', '                    _sink.Position(0, ChildScan.None);\n', ''),
    ('A20', 'the maximum LCP is the smallest', 'if (cur > _maxLcp) _maxLcp = cur;', 'if (cur < _maxLcp) _maxLcp = cur;'),
    ('A21', 'the tagged sink is chosen for the packed codec', '(not null, null, ChildCodec.Tagged) => new ScanCore<TaggedChildSink>', '(not null, null, ChildCodec.Packed) => new ScanCore<TaggedChildSink>'),
    ('A22', 'the first l-index of an interval is its left boundary', 'new IntervalEntry { Ell = cur, Lb = left, First = i, Last = i }', 'new IntervalEntry { Ell = cur, Lb = left, First = left, Last = i }'),
    ('A23', 'a next link is written from the wrong l-index', '_sink.LIndex(t.Last, i);\n            t.Last = i;', '_sink.LIndex(t.First, i);\n            t.Last = i;'),
    ('A24', 'the PSV sink writes one entry late', 'public void Position(uint i, uint psv) => writer.Write32(i, psv);', 'public void Position(uint i, uint psv) => writer.Write32(i == 0 ? 0 : i - 1, psv);'),
    ('A25', 'the stack peak is not reported', 'StackPeak = _stack.PeakHeight,', 'StackPeak = 0,'),
    ('A26', 'the root is not closed with the last-child flag set when it is popped by Complete (first l-index lost)', 'Close(e.Ell, e.Lb, (uint)_n, e.First, true);\n            }', 'Close(e.Ell, e.Lb, (uint)_n, e.Last, true);\n            }'),
])
M['nsvscan'] = (APP + 'Ph5/NsvScan.cs', [
    ('N01', 'equal values are not popped', 'stack.Top.Value >= x) stack.Pop();', 'stack.Top.Value > x) stack.Pop();'),
    ('N02', 'the "no smaller value" is n', 'uint none = (uint)(n + 1);', 'uint none = (uint)n;'),
    ('N03', 'the pushed index is one too high', 'stack.Push(new ValueEntry { Value = x, Index = rank });', 'stack.Push(new ValueEntry { Value = x, Index = rank + 1 });'),
    ('N04', 'the NSV windows are not flushed', 'new AlignedWindowWriter(nsvPath, 4 * (n + 1), options.WindowBytes, 2, options.FlushWindows);', 'new AlignedWindowWriter(nsvPath, 4 * (n + 1), options.WindowBytes, 2, false);'),
    ('N05', 'the scan goes forwards through each block', 'for (int j = values.Length - 1; j >= 0; j--)', 'for (int j = 0; j < values.Length; j++)'),
    ('N06', 'the stack peak is not reported', 'StackPeak = stack.PeakHeight,', 'StackPeak = 0,'),
])
M['externalchild'] = (APP + 'Ph5/ExternalChild.cs', [
    ('X01', 'the parallel option is inverted', 'outputs.Nsv is not null && options.ParallelBackward ?', 'outputs.Nsv is not null && !options.ParallelBackward ?'),
    ('X02', 'the outputs of the forward scan are swapped', 'new ChildScanOutputs { ChildPath = outputs.Child, PsvPath = outputs.Psv }', 'new ChildScanOutputs { ChildPath = outputs.Psv, PsvPath = outputs.Child }'),
    ('X03', 'the volume of the LCP files counts the base array only', 'LcpBytes = new FileInfo(lcpBase).Length + new FileInfo(lcpOverflow).Length', 'LcpBytes = new FileInfo(lcpBase).Length'),
    ('X04', 'the NSV pass is skipped when it is not concurrent', 'else if (outputs.Nsv is not null) backward = Backward', 'else if (outputs.Nsv is null) backward = Backward'),
])
M['nextsym'] = (APP + 'Ph5/NextSymScan.cs', [
    ('M01', 'the end mark is stored one position late', 'destination[i] = p < n ? (byte)text.Symbol(p) : NoSymbol;', 'destination[i] = p <= n ? (byte)text.Symbol(p) : NoSymbol;'),
    ('M02', 'the LCP is not added', 'long p = (long)sa[i] + lcp[i];', 'long p = (long)sa[i];'),
    ('M03', 'a 9-bit text is accepted', 'if (text.Bits > 8) throw', 'if (text.Bits > 9) throw'),
    ('M04', 'a suffix array file one entry short is accepted', 'if (sa.Length != (n + 1) * 4) throw', 'if (sa.Length < n * 4) throw'),
    ('M05', 'the block is written one rank early', 'writer.Write(first + at, bytes.AsSpan(0, count));', 'writer.Write(first + at + 1, bytes.AsSpan(0, count - 1));'),
])
M['e5'] = (APP + 'Experiments/E5.cs', [
    ('E01', 'the RssAnon cap is exclusive', 'rssAnonMb <= capMb && scratchEntriesLeft == 0', 'rssAnonMb < capMb && scratchEntriesLeft == 0'),
    ('E02', 'Healthy ignores the arena', 'scratchLiveBytes == 0 && arenaLiveBytes == 0', 'scratchLiveBytes == 0'),
    ('E03', 'Healthy ignores the tracked scratch', 'scratchEntriesLeft == 0 && scratchLiveBytes == 0', 'scratchEntriesLeft == 0'),
    ('E04', 'a difference of one file does not fail the run', 'ok &= equal;', 'ok |= equal;'),
    ('E05', 'the layout "stats" is refused', 'layout is not ("a" or "f" or "both" or "stats")', 'layout is not ("a" or "f" or "both")'),
    ('E06', 'windows are flushed only with --no-flush', 'FlushWindows = !o.Flag("no-flush"),', 'FlushWindows = o.Flag("no-flush"),'),
    ('E07', 'the statistics comparison ignores NODE_COUNT', 'bool stats = expected["nodecount"] == f.NodeCount.ToString(CultureInfo.InvariantCulture)\n            && ', 'bool stats = true\n            && '),
    ('E08', 'the statistics comparison ignores LRS_RANK', '\n            && expected["lrsrank"] == f.LrsRank.ToString(CultureInfo.InvariantCulture);', ';'),
    ('E09', 'the packed codec is chosen for tagged', '"tagged" => ChildCodec.Tagged,', '"tagged" => ChildCodec.Packed,'),
    ('E10', 'the oracle writes the packed array as the tagged one', 'WriteRaw(outPrefix + ".cldT", (uint[])typeof(TaggedChild)', 'WriteRaw(outPrefix + ".cldT", (uint[])typeof(PackedChild)'),
    ('E11', 'the scan result of the NEXTSYM verification is ignored', 'ok = equal;\n        }\n        return ok;', 'ok = true;\n        }\n        return ok;'),
    ('E12', 'the fused scan is created at the second chunk', 'if (_scan is null)\n            {\n                sampler?.Mark();', 'if (_scan is null && firstRank > 0)\n            {\n                sampler?.Mark();'),
])
M['dirtysampler'] = (APP + 'Harness/DirtySampler.cs', [
    ('D01', 'the mapping filter is a suffix match', 'inside = slash >= 0 && line[slash..].StartsWith(prefix);', 'inside = slash >= 0 && line[slash..].EndsWith(prefix);'),
    ('D02', 'Shared_Dirty is not summed', 'else if (inside && (line.StartsWith(KeyPrivate) || line.StartsWith(KeyShared)))', 'else if (inside && line.StartsWith(KeyPrivate))'),
    ('D03', 'every mapping is summed', 'inside = slash >= 0 && line[slash..].StartsWith(prefix);', 'inside = slash >= 0;'),
    ('D04', 'Private_Dirty is not summed', 'else if (inside && (line.StartsWith(KeyPrivate) || line.StartsWith(KeyShared)))', 'else if (inside && line.StartsWith(KeyShared))'),
])

bad = 0
total = 0
for name, (path, mutants) in M.items():
    src = open(os.path.join(ROOT, path)).read()
    out = []
    for mid, desc, old, new in mutants:
        total += 1
        n = src.count(old)
        if n != 1:
            print('NOT UNIQUE (%d): %s %s' % (n, name, mid))
            bad += 1
        out.append({'id': mid, 'desc': desc, 'old': old, 'new': new})
    json.dump(out, open(os.path.join(HERE, 'e5_%s_mutants.json' % name), 'w'), indent=1)
print('%d mutants in %d files, %d patterns not unique' % (total, len(M), bad))
