#!/usr/bin/env python3
"""Round 3 of the mutation audit of 1b.7c: fresh mutants of the code that the tests written after round 1 reach (the plumbing of every parameter into its stage, the final-length check, the ranges of the CHILD writer, the
   buffers of PH4, the window bookkeeping), appended to the lists that make_1b7c_round2.py wrote for the survivors of round 2 (1b7c_r3_<file>_listed.json; run make_1b7c_round2.py r3 r2 first).
   usage: make_1b7c_fresh.py <repo root>"""
import json, os, sys

ROOT = sys.argv[1]
SRC = 'src/SuffixTree/Algorithms/SuffixTree.Esa/External/'


def m(i, desc, old, new):
    return {'id': i, 'desc': desc, 'old': old, 'new': new}


FRESH = {
    'LcpChildBuilder.cs': ('lcpchildbuilder', [
        m('F01', 'PH4 does not get the read buffer of the rule', 'ReadBufferBytes = p.LcpReadBytes,', 'ReadBufferBytes = 1 << 20,'),
        m('F02', 'PH4 does not get the tuple buffer of the rule', 'TupleBufferBytes = p.LcpTupleBytes,', 'TupleBufferBytes = 32 * 1024,'),
        m('F03', 'PH4 does not get the window of the rule', 'WindowBytes = p.WindowBytes,', 'WindowBytes = 1 << 20,'),
        m('F04', 'PH4 does not get the chunk of the rule', 'ChunkPositions = p.LcpChunk,', 'ChunkPositions = 8_000_000,'),
        m('F06', 'PH5 does not get the block of the rule', 'BlockEntries = p.BlockEntries,', 'BlockEntries = 256 * 1024,'),
        m('F07', 'PH5 does not get the window of the rule', 'WindowBytes = p.ChildWindowBytes,', 'WindowBytes = 1 << 20,'),
        m('F08', 'PH5 does not get the slots of the rule', 'WindowSlots = p.ChildSlots,', 'WindowSlots = 8,'),
        m('F09', 'PH5 does not flush its windows', 'FlushWindows = true,', 'FlushWindows = false,'),
        m('F10', 'PH5 does not get the spill flush of the rule', 'SpillFlushBytes = p.SpillFlushBytes,\n                Codec', 'SpillFlushBytes = 0,\n                Codec'),
        m('F11', 'PH5 always uses the packed codec', 'Codec = options.Codec,\n                Budget = childBudget', 'Codec = ChildCodec.Packed,\n                Budget = childBudget'),
        m('F12', 'PH4 is not charged to its budget', 'Budget = lcpBudget,', 'Budget = null,'),
        m('F13', 'PH5 is not charged to its budget', 'Budget = childBudget,', 'Budget = null,'),
        m('F14', 'the PH4 peak is the PH5 peak', 'ArenaPeakBytes = lcpBudget.PeakBytes,', 'ArenaPeakBytes = childBudget.PeakBytes,'),
        m('F15', 'the PH5 peak is the PH4 peak', 'ChildArenaPeakBytes = childBudget.PeakBytes,', 'ChildArenaPeakBytes = lcpBudget.PeakBytes,'),
        m('F16', 'the result carries no counters', 'Io = counters.Snapshot(),', 'Io = default,'),
        m('F17', 'the result carries no time', 'Seconds = seconds,', 'Seconds = 0,'),
        m('F18', 'the result carries no reservation outcome', 'Reservation = lcp.Reservation,', 'Reservation = default,'),
        m('F19', 'the trace hook is not handed to PH4', '                Trace = options.Trace,\n', ''),
        m('F21', 'the RSS sampler is not stopped for the result', 'RssStats rss = sampler?.Stop() ?? default;', 'RssStats rss = default;'),
        m('F25', 'CheckLengths: the base LCP is one entry short', 'Expect(outputs.Lcp, LcpWidth.Bytes(lcp.Kind) * (n + 1), "the base LCP");', 'Expect(outputs.Lcp, LcpWidth.Bytes(lcp.Kind) * n, "the base LCP");'),
        m('F26', 'CheckLengths: the overflow pairs are 4 bytes', 'Expect(outputs.Overflow, 8 * lcp.OverflowCount, "the overflow list");', 'Expect(outputs.Overflow, 4 * lcp.OverflowCount, "the overflow list");'),
        m('F27', 'CheckLengths: the child table is one entry short', 'Expect(childPath, 4 * (n + 1), "the child table");', 'Expect(childPath, 4 * n, "the child table");'),
        m('F28', 'CheckLengths: a shorter file is accepted', 'if (got != bytes) throw new InvalidOperationException($"{what} has', 'if (got > bytes) throw new InvalidOperationException($"{what} has'),
        m('F29', 'CheckLengths: a longer file is accepted', 'if (got != bytes) throw new InvalidOperationException($"{what} has', 'if (got < bytes) throw new InvalidOperationException($"{what} has'),
        m('F26b', 'the lazy scan does not dispose its scan', '            _scan?.Dispose();\n            _file?.Dispose();', '            _file?.Dispose();'),
        m('F26c', 'the lazy scan does not dispose its file', '            _scan?.Dispose();\n            _file?.Dispose();', '            _scan?.Dispose();'),
    ]),
    'Lcp/LcpBuild.cs': ('lcpbuild', [
        m('F30', 'B1: the blocks overlap by one entry', 'for (long first = 1; first <= n; first += readEntries)', 'for (long first = 1; first <= n; first += readEntries - 1)'),
        m('F31', 'B1: the last block takes one entry too few', 'int take = (int)Math.Min(readEntries, n + 1 - first);\n            sa.ReadExactAt(MemoryMarshal.AsBytes(read.AsSpan(take)), first * 4);\n            for (int i = 0; i < take; i++)\n            {\n                uint s = read[i];\n                if (s >= n)', 'int take = (int)Math.Min(readEntries, n - first);\n            sa.ReadExactAt(MemoryMarshal.AsBytes(read.AsSpan(take)), first * 4);\n            for (int i = 0; i < take; i++)\n            {\n                uint s = read[i];\n                if (s >= n)'),
        m('F32', 'B2: the Phi value of a position is one too large', 'phi[(int)at] = t.Phi;', 'phi[(int)at] = t.Phi + 1;'),
        m('F33', 'B2: the rank of a position is one too large', 'rank[(int)at] = t.Rank;', 'rank[(int)at] = t.Rank + 1;'),
        m('F34', 'B3: the value of a rank is one too large', 'dense[(int)at] = t.Value;', 'dense[(int)at] = t.Value + 1;'),
        m('F35', 'B3: the block handed on is one entry short', 'Deliver(output, start, dense.AsSpan(length));', 'Deliver(output, start, dense.AsSpan(length - 1));'),
        m('F36', 'the workers are the threads, never the request', 'int workers = _o.Workers > 0 ? _o.Workers : _o.Threads;', 'int workers = _o.Threads;'),
        m('F37', 'a read buffer smaller than a tuple holds two tuples', 'private int ReadEntries(int tupleBytes) => Math.Max(1, _o.ReadBufferBytes / tupleBytes);', 'private int ReadEntries(int tupleBytes) => Math.Max(2, _o.ReadBufferBytes / tupleBytes);'),
        m('F39', 'the text length is reported one too large', '_stats.N = n;', '_stats.N = n + 1;'),
        m('F40', 'B3: the first rank of a chunk is one less', '                long start = c * _g.R;\n                long expected', '                long start = c * (_g.R - 1);\n                long expected'),
        m('F43', 'the header check reads SA[1]', 'sa.ReadExactAt(MemoryMarshal.AsBytes(first), 0);', 'sa.ReadExactAt(MemoryMarshal.AsBytes(first), 4);'),
        m('F44', 'WorkerBytes: the tuple buffers are charged for the text chunks', 'sizeof(RankTuple) * g.RankChunks;', 'sizeof(RankTuple) * g.TextChunks;'),
        m('F45', 'WorkerBytes: the read buffer is sized for rank tuples', '(Math.Max(1, o.ReadBufferBytes / sizeof(PhiTuple)) * (long)sizeof(PhiTuple))', '(Math.Max(1, o.ReadBufferBytes / sizeof(RankTuple)) * (long)sizeof(PhiTuple))'),
    ]),
    'Child/AlignedWindowWriter.cs': ('alignedwindowwriter', [
        m('F50', 'the bit set of the windows has no room for a partial last window', 'long windows = (length + windowBytes - 1) >> _shift;', 'long windows = length >> _shift;'),
        m('F51', 'the bit set of the windows is one word short', '_seen = new ulong[(windows + 63) / 64];', '_seen = new ulong[windows / 64];'),
        m('F52', 'a free slot is not preferred to a used one', '                victim = i;\n                oldest = long.MinValue;', '                victim = i;'),
        m('F53', 'a free slot found earlier does not stop the search for the oldest', 'else if (oldest != long.MinValue && _slots[i].LastUse < oldest)', 'else if (_slots[i].LastUse < oldest)'),
        m('F54', 'the flushes are not counted', '        slot.View!.Flush();\n        _flushes++;', '        slot.View!.Flush();\n        _flushes += 0;'),
        m('F55', 'the size of a window is not recorded', '        slot.Size = size;\n', ''),
        m('F56', 'a range beyond the end of the file is accepted', 'if (start > file.Length - length)', 'if (start > file.Length)'),
    ]),
    'Child/ChildScan.cs': ('childscan', [
        m('F60', 'the left boundary of a pushed interval is the rank itself', 'uint left = i - 1;', 'uint left = i;'),
        m('F61', 'a closing interval is closed with the stale top', '            _stack.Pop();\n            top = ref _stack.Top;\n            Close(e.Ell, e.Lb, i - 1, e.First, cur < top.Ell);', '            _stack.Pop();\n            Close(e.Ell, e.Lb, i - 1, e.First, cur < top.Ell);'),
        m('F62', 'the stack of the empty text is closed too', '        if (_n > 0)\n        {\n            while (_stack.Height > 0)', '        if (_n >= 0)\n        {\n            while (_stack.Height > 0)'),
        m('F64', 'a pushed interval has no last l-index', 'Lb = left, First = i, Last = i });', 'Lb = left, First = i, Last = 0 });'),
        m('F65', 'the root is recognised by a first l-index of 0', 'else if (top.First == ChildScan.None)', 'else if (top.First == 0)'),
        m('F67', 'the scan ignores the start of its range', 'new AlignedWindowWriter(file, start, 4 * (n + 1),', 'new AlignedWindowWriter(file, 0, 4 * (n + 1),'),
    ]),
    'Lcp/LcpOutputWriter.cs': ('lcpoutputwriter', [
        m('F70', 'only half of the pairs are appended', '_overflow.Append(MemoryMarshal.AsBytes(_pairs.AsSpan(pairs * 2)));', '_overflow.Append(MemoryMarshal.AsBytes(_pairs.AsSpan(pairs)));'),
        m('F71', 'a block is cut into slices of the larger size', 'Math.Min(_slice, values.Length - at)', 'Math.Max(_slice, values.Length - at)'),
        m('F72', 'the byte buffer holds one byte per slice entry whatever the width', '_bytes = arena.Allocate<byte>(_slice * _width, zero: false, "lcp.out");', '_bytes = arena.Allocate<byte>(_slice, zero: false, "lcp.out");'),
    ]),
}

HERE = os.path.dirname(os.path.abspath(__file__))
total = 0
for name in sorted(FRESH):
    base, lst = FRESH[name]
    text = open(os.path.join(ROOT, SRC + name)).read()
    for x in lst:
        if text.count(x['old']) != 1:
            raise SystemExit('%s %s: pattern occurs %d times' % (name, x['id'], text.count(x['old'])))
        x['file'] = SRC + name
    out = os.path.join(HERE, '1b7c_r3_%s_listed.json' % base)
    cur = json.load(open(out)) if os.path.exists(out) else []
    have = {x['id'] for x in cur}
    cur += [x for x in lst if x['id'] not in have]
    json.dump(cur, open(out, 'w'), indent=1)
    total += len(lst)
    print('%-24s +%d fresh, %d in the list' % (base, len(lst), len(cur)))
print('fresh mutants:', total)
