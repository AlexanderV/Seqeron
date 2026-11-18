#!/usr/bin/env python3
"""Round 3 of the audit of 1b.7a: fresh hand-made mutants, written after rounds 1 and 2 had shown which parts of the code the tests pin hardest and which they pin least (1b7a_r3_<file>_listed.json).
   Each pattern must occur exactly once in its file (checked here). usage: make_1b7a_round3.py <repo root>"""
import json, os, sys

ROOT = sys.argv[1]
SRC = 'src/SuffixTree/Algorithms/SuffixTree.Esa/External/'


def m(i, desc, old, new):
    return {'id': i, 'desc': desc, 'old': old, 'new': new}


FILES = {
    'dist': ('Dist.cs', [
        m('J01', 'the last source read is as long as a full one', 'int take = (int)Math.Min(readEntries, sourceEntries - first);\n            Span<uint> entries', 'int take = (int)Math.Min(readEntries, sourceEntries);\n            Span<uint> entries'),
        m('J02', 'the start slot of a chunk is one short', 'long startSlot = chunk * chunkEntries;', 'long startSlot = (chunk * chunkEntries) - (chunk > 0 ? 1 : 0);'),
        m('J03', 'the last chunk array is as long as a full one', 'int length = (int)Math.Min(chunkEntries, slots - startSlot);', 'int length = (int)Math.Min(chunkEntries, slots);'),
        m('J04', 'the chunk array is one entry short', '(int)Math.Min(chunkEntries, slots), zero: false, "dist.chunk"', '(int)Math.Min(chunkEntries, slots) - (slots > 1 ? 1 : 0), zero: false, "dist.chunk"'),
        m('J05', 'the source is read in steps one entry longer than the buffer', 'for (long first = 0; first < sourceEntries; first += readEntries)', 'for (long first = 0; first < sourceEntries; first += readEntries + 1)'),
        m('J06', 'a scattered slot gets the key', '            slot = e.Value;\n', '            slot = default(TValue) is uint ? (TValue)(object)e.Key : e.Value;\n'),
        m('J07', 'the ISA key is flipped in its lowest bit', 'key = sourceValue;\n        value = (uint)index;', 'key = sourceValue ^ 1;\n        value = (uint)index;'),
        m('J08', 'the Phi mapper checks against n + 1', 'var mapper = new PhiMapper(n);', 'var mapper = new PhiMapper(n + 1);'),
        m('J09', 'the exactly-once check starts at the second slot', 'for (int k = 0; k < length; k++)\n                {\n                    if (array[k].Equals(target.Fill))', 'for (int k = 1; k < length; k++)\n                {\n                    if (array[k].Equals(target.Fill))'),
        m('J10', 'the written-twice test compares with the value', 'if (exactlyOnce && !slot.Equals(fill)) ThrowWrittenTwice(e.Key);', 'if (exactlyOnce && !slot.Equals(e.Value)) ThrowWrittenTwice(e.Key);'),
        m('J11', 'the tuple size is taken from the key only', 'int tupleBytes = sizeof(DistEntry<TValue>);', 'int tupleBytes = sizeof(uint) + 4;'),
        m('J12', 'the target is sized in tuples', 'new MappedWindowWriter(targetPath, slots * sizeof(TValue), options.WindowBytes, counters)', 'new MappedWindowWriter(targetPath, slots * sizeof(DistEntry<TValue>), options.WindowBytes, counters)'),
    ]),
    'spillstack': ('SpillStack.cs', [
        m('K01', 'a block is written one block too far', '_file.WriteAt(new ReadOnlySpan<byte>(buffer, _blockBytes), block * _blockBytes);', '_file.WriteAt(new ReadOnlySpan<byte>(buffer, _blockBytes), (block + 1) * _blockBytes);'),
        m('K02', 'the evicted block is the wrong one', '            if (_otherDirty) WriteBlock(_topBlock - 1, _other);', '            if (_otherDirty) WriteBlock(_topBlock, _other);'),
        m('K03', 'a block is read from the wrong offset', '_file!.ReadExactAt(new Span<byte>(buffer, _blockBytes), block * _blockBytes);', '_file!.ReadExactAt(new Span<byte>(buffer, _blockBytes), (block * _blockBytes) + 1);'),
        m('K04', 'blocks are sized as 16-byte entries', '_blockBytes = checked(blockEntries * sizeof(T));', '_blockBytes = checked(blockEntries * 16);'),
        m('K05', 'a retreat goes two blocks down', '        _topBlock--;\n        _count = _block;', '        _topBlock -= 2;\n        _count = _block;'),
        m('K06', 'the unflushed count is not reset after a flush', '        _flushes++;\n        _unflushed = 0;', '        _flushes++;'),
        m('K07', 'flushes are not counted', '        _flushes++;\n', ''),
        m('K08', 'a push does not advance the count', '        _cur[_count++] = value;', '        _cur[_count] = value;'),
        m('K09', 'the stack takes its second buffer from the first', '_other = arena.Allocate<T>(blockEntries, zero: false, "spill.b").Pointer;', '_other = _cur;'),
        m('K10', 'a read-back block is not marked as resident below nothing', '        _otherDirty = false;\n        _otherBelow = false;\n        _topBlock--;', '        _otherDirty = false;\n        _topBlock--;'),
    ]),
    'procio': ('ProcIo.cs', [
        m('L01', 'a newline may precede the number', "while (i < text.Length && (text[i] == (byte)' ' || text[i] == (byte)'\\t')) i++;", "while (i < text.Length && (text[i] == (byte)' ' || text[i] == (byte)'\\t' || text[i] == (byte)'\\n')) i++;"),
        m('L02', 'the search restarts at the same position after a mismatch', 'from = position + 1;', 'from = position + key.Length;'),
        m('L03', 'the digit loop stops at the first zero', "text[i] >= (byte)'0' && text[i] <= (byte)'9')", "text[i] > (byte)'0' && text[i] <= (byte)'9')"),
        m('L04', 'the proc snapshot buffer is 16 bytes', 'Span<byte> buffer = stackalloc byte[2048];', 'Span<byte> buffer = stackalloc byte[16];'),
        m('L05', 'the real source stops after one read', '            while (total < buffer.Length)\n            {\n                int got = RandomAccess.Read(handle, buffer[total..], total);', '            while (total < 1)\n            {\n                int got = RandomAccess.Read(handle, buffer[total..], total);'),
    ]),
    'rssreader': ('RssReader.cs', [
        m('L06', 'the status buffer is 64 bytes', 'Span<byte> buffer = stackalloc byte[8192];', 'Span<byte> buffer = stackalloc byte[64];'),
        m('L07', 'Stop takes no sample when the thread was never started', '        _started = false;\n        return new RssStats(_available,', '        _started = false;\n        return new RssStats(true,'),
        m('L08', 'the sampler reports available before it has a sample', '        _available = true;\n        Record(first);', '        _available = true;'),
        m('L09', 'the stats report the maximum RssAnon as the last', 'Interlocked.Read(ref _maxThreads), Interlocked.Read(ref _lastAnon));', 'Interlocked.Read(ref _maxThreads), Interlocked.Read(ref _maxAnon));'),
    ]),
    'arenabudget': ('ArenaBudget.cs', [
        m('P01', 'the peak is the largest single request', 'Atomic.RaiseTo(ref _peak, now);', 'Atomic.RaiseTo(ref _peak, bytes);'),
        m('P02', 'the arena charges elements instead of bytes', 'bytes = (long)length * sizeof(T);', 'bytes = length;'),
        m('P03', 'ReleaseTo keeps the sizes of the blocks it freed', '            freed += _sizes[i];\n            _sizes.RemoveAt(i);', '            freed += _sizes[i];'),
        m('P04', 'the reservation message omits the live bytes', "({bytes} bytes; {seen} live)", "({bytes} bytes)"),
    ]),
    'distscratch': ('DistScratch.cs', [
        m('N01', 'the chunk file name does not contain the chunk', 'Path.Combine(_dir, $"{_tag}-{chunk:D8}.bin")', 'Path.Combine(_dir, $"{_tag}.bin")'),
        m('N02', 'the writer capacity counts bytes', '_capacity = Math.Max(1, bufferBytes / sizeof(T));', '_capacity = Math.Max(1, bufferBytes);'),
        m('N03', 'the pending count is the capacity', 'internal int Pending(long chunk) => _filled[Local(chunk)];', 'internal int Pending(long chunk) => Local(chunk) >= 0 ? _capacity : 0;'),
    ]),
    'mappedwindowwriter': ('MappedWindowWriter.cs', [
        m('M01', 'the window is mapped read-only', 'CreateViewAccessor(offset, size, MemoryMappedFileAccess.Write)', 'CreateViewAccessor(offset, size, MemoryMappedFileAccess.Read)'),
        m('M02', 'the window ends one byte early', '        _viewEnd = offset + size;', '        _viewEnd = offset + size - 1;'),
    ]),
    'iocounters': ('IoCounters.cs', [
        m('O01', 'a snapshot reports the scratch peak as the live bytes', 'MappedWindows, MappedMaxWindowBytes, ScratchLiveBytes, ScratchPeakBytes);', 'MappedWindows, MappedMaxWindowBytes, ScratchPeakBytes, ScratchPeakBytes);'),
        m('O02', 'Since subtracts the mapped windows from the writes', '            MappedWriteBytes - earlier.MappedWriteBytes,\n            MappedWindows - earlier.MappedWindows,', '            MappedWriteBytes - earlier.MappedWriteBytes,\n            MappedWindows,'),
    ]),
}

total = 0
for base, (fname, muts) in FILES.items():
    path = os.path.join(ROOT, SRC + fname)
    text = open(path).read()
    out = []
    for x in muts:
        n = text.count(x['old'])
        if n != 1:
            print('PATTERN %s occurs %d times in %s: %s' % (x['id'], n, fname, x['old'][:70].replace('\n', '\\n')))
            continue
        out.append(dict(x, file=SRC + fname))
    json.dump(out, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), '1b7a_r3_%s_listed.json' % base), 'w'), indent=1)
    total += len(out)
    print('%-24s %3d mutants' % (fname, len(out)))
print('total', total)
