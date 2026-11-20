#!/usr/bin/env python3
"""Hand-made single-site mutants of the code that work package 1b.9b added or changed, and of the code whose tests it closed gaps for (prefix 1b9b_): the guard counter of the arena, the BUILD_INFO room check, the
   padding gaps, the PSI primary row and LRS_DEPTH comparisons of the Structure level, the bounds of the mapped child table, the 32-bit LCP reader, the equality case of the corrupt-leaf-depth guard (navigator and contract tree),
   the tie-break of the longest repeat (three scans), the reader count of the facade for in-memory and file-backed trees, and the in-memory pattern descent. Each mutant is written to one JSON list per source file
   (1b9b_<file>_listed.json); every `old` must occur exactly once in its file and must change it (checked here). The runner (run_1b9b_mutants.py) builds the mutant in a separate worktree and runs the unit tests of the new
   and the closed-gap tests, and for a survivor the budgeted harness of 1b.8.
   usage: make_1b9b_mutants.py <repo root>"""
import json, os, sys

ROOT = sys.argv[1]
HERE = os.path.dirname(os.path.abspath(__file__))
SRC = 'src/SuffixTree/Algorithms/SuffixTree.Esa/'


def m(i, desc, old, new):
    return {'id': '1b9b_' + i, 'desc': desc, 'old': old, 'new': new}


FILES = {
    SRC + 'Memory/NativeArena.cs': [
        m('A01', 'the guard counter never moves', "if (options.GuardBytes > 0) Interlocked.Increment(ref _guardedCreated);", "if (options.GuardBytes < 0) Interlocked.Increment(ref _guardedCreated);"),
        m('A02', 'the guard counter counts every arena', "if (options.GuardBytes > 0) Interlocked.Increment(ref _guardedCreated);", "if (options.GuardBytes >= 0) Interlocked.Increment(ref _guardedCreated);"),
        m('A03', 'the guard counter is not atomic', "if (options.GuardBytes > 0) Interlocked.Increment(ref _guardedCreated);", "if (options.GuardBytes > 0) _guardedCreated++;"),
    ],
    SRC + 'External/Sesa/SesaBuilder.cs': [
        m('B01', 'B102: the BWT pass gets no guards', "ArenaGuards = options.ArenaGuards,\n            Counters = counters,\n            Files = options.Files,", "ArenaGuards = null,\n            Counters = counters,\n            Files = options.Files,"),
        m('B02', 'the PSI pass gets no guards', "ArenaGuards = options.ArenaGuards,\n                Counters = counters,\n            });\n            double secPsi", "ArenaGuards = null,\n                Counters = counters,\n            });\n            double secPsi"),
        m('B03', 'B118: BUILD_INFO may overrun its room by 4 KiB', "if (length > capacity) throw", "if (length > capacity + 4096) throw"),
        m('B04', 'BUILD_INFO may overrun its room by one byte', "if (length > capacity) throw", "if (length > capacity + 1) throw"),
        m('B05', 'BUILD_INFO that fills its room exactly is refused', "if (length > capacity) throw", "if (length >= capacity) throw"),
        m('B06', 'the room check is not made by the builder', "CheckBuildInfoFits(infoBytes.Length, plan.BuildInfoCapacity);", "_ = plan.BuildInfoCapacity;"),
        m('B07', 'the SA pass gets no guards', "ArenaGuards = options.ArenaGuards,\n                VerifyRanks = options.VerifyRanks,", "ArenaGuards = null,\n                VerifyRanks = options.VerifyRanks,"),
        m('B08', 'the LCP and child pass gets no guards', "ArenaGuards = options.ArenaGuards,\n                MeasureRss = options.MeasureRss,\n                ProcSource = options.ProcSource,\n                RssAnonCapMb", "ArenaGuards = null,\n                MeasureRss = options.MeasureRss,\n                ProcSource = options.ProcSource,\n                RssAnonCapMb"),
    ],
    SRC + 'External/Sesa/SesaVerifier.cs': [
        m('V01', 'VE7: the first eight bytes of a padding gap are not checked', "if (s.Offset > end && new ReadOnlySpan<byte>(file.Data + end, checked((int)(s.Offset - end))).IndexOfAnyExcept((byte)0) >= 0)", "if (s.Offset > end + 8 && new ReadOnlySpan<byte>(file.Data + end + 8, checked((int)(s.Offset - end - 8))).IndexOfAnyExcept((byte)0) >= 0)"),
        m('V02', 'only the first eight bytes of a padding gap are checked', "new ReadOnlySpan<byte>(file.Data + end, checked((int)(s.Offset - end))).IndexOfAnyExcept((byte)0) >= 0)", "new ReadOnlySpan<byte>(file.Data + end, checked((int)Math.Min(8, s.Offset - end))).IndexOfAnyExcept((byte)0) >= 0)"),
        m('V03', 'VE9: a PSI primary row below the real one is accepted', "if (ps.Psi!.Param0 != r) throw", "if (ps.Psi!.Param0 > r) throw"),
        m('V04', 'a PSI primary row above the real one is accepted', "if (ps.Psi!.Param0 != r) throw", "if (ps.Psi!.Param0 < r) throw"),
        m('V05', 'VE12: LRS_DEPTH above the real one is accepted', "if (h.LrsDepth != sink.LrsDepth) throw", "if (h.LrsDepth < sink.LrsDepth) throw"),
        m('V06', 'LRS_DEPTH below the real one is accepted', "if (h.LrsDepth != sink.LrsDepth) throw", "if (h.LrsDepth > sink.LrsDepth) throw"),
        m('V07', 'the verifier keeps the last deepest interval, not the first', "            if (ell > LrsDepth)", "            if (ell >= LrsDepth)"),
    ],
    SRC + 'ChildLayouts.cs': [
        m('C01', 'the first l-index may be read one cell past the table', "if (_checked && (uint)at >= (uint)_count) ThrowOutside(at);", "if (_checked && (uint)at > (uint)_count) ThrowOutside(at);"),
        m('C02', 'the next l-index of a tagged table may be read one cell past it', "if (_checked && (uint)k >= (uint)_count) ThrowOutside(k);\n            uint cell", "if (_checked && (uint)k > (uint)_count) ThrowOutside(k);\n            uint cell"),
        m('C03', 'the next l-index of a packed table may be read one cell past it', "if (_checked && (uint)k >= (uint)_count) ThrowOutside(k);\n        int lk", "if (_checked && (uint)k > (uint)_count) ThrowOutside(k);\n        int lk"),
    ],
    SRC + 'MappedSearch.cs': [
        m('S01', 'S03: a 32-bit LCP value up to length + 5 is accepted', "if (TAccess.IsChecked && (uint)v > (uint)length) MappedSearch.ThrowDamaged();", "if (TAccess.IsChecked && (uint)v > (uint)length + 5) MappedSearch.ThrowDamaged();"),
        m('S02', 'a 32-bit LCP value equal to the length is refused', "if (TAccess.IsChecked && (uint)v > (uint)length) MappedSearch.ThrowDamaged();", "if (TAccess.IsChecked && (uint)v >= (uint)length) MappedSearch.ThrowDamaged();"),
        m('S04', 'the symbol read keeps the bits that the shift moved out of the 16', "return (int)(((w << (int)(bit & 7)) & 0xFFFF) >> (16 - bits));", "return (int)((w << (int)(bit & 7)) >> (16 - bits));"),
        m('S05', 'the symbol read takes one bit too many', "return (int)(((w << (int)(bit & 7)) & 0xFFFF) >> (16 - bits));", "return (int)(((w << (int)(bit & 7)) & 0xFFFF) >> (15 - bits));"),
        m('S06', 'the symbol read starts at the bit offset modulo 4', "return (int)(((w << (int)(bit & 7)) & 0xFFFF) >> (16 - bits));", "return (int)(((w << (int)(bit & 3)) & 0xFFFF) >> (16 - bits));"),
        m('S07', 'the symbol read does not turn the big-endian stream around', "uint w = BinaryPrimitives.ReverseEndianness(Unsafe.ReadUnaligned<ushort>(text + (bit >> 3)));", "uint w = Unsafe.ReadUnaligned<ushort>(text + (bit >> 3));"),
        m('S03', 'a 32-bit LCP value is never checked', "if (TAccess.IsChecked && (uint)v > (uint)length) MappedSearch.ThrowDamaged();", "if (TAccess.IsChecked && false) MappedSearch.ThrowDamaged();"),
    ],
    SRC + 'EsaNavigator.cs': [
        m('N01', 'the navigator accepts a leaf as deep as its parent', "if (_isMapped && _n - saA + 1 <= parentDepth) ThrowCorruptChildTable();", "if (_isMapped && _n - saA + 1 < parentDepth) ThrowCorruptChildTable();"),
    ],
    SRC + 'EsaSuffixIndex.Tree.cs': [
        m('T01', 'the contract tree accepts a leaf as deep as its parent', "if (_isMapped && _n - saA + 1 <= parentDepth) ThrowCorruptChildTable();", "if (_isMapped && _n - saA + 1 < parentDepth) ThrowCorruptChildTable();"),
    ],
    SRC + 'EsaSuffixTree.cs': [
        m('F01', 'an in-memory call does not look at the disposed flag', "ObjectDisposedException.ThrowIf(Volatile.Read(ref _disposed) != 0, this);\n            return new ReadScope(this, -1);", "_ = Volatile.Read(ref _disposed);\n            return new ReadScope(this, -1);"),
        m('F02', 'an in-memory call throws for a tree that is not disposed (the flag test is inverted)', "ObjectDisposedException.ThrowIf(Volatile.Read(ref _disposed) != 0, this);\n            return new ReadScope(this, -1);", "ObjectDisposedException.ThrowIf(Volatile.Read(ref _disposed) == 0, this);\n            return new ReadScope(this, -1);"),
        m('F03', 'a tree over a container takes the in-memory path of Enter (no reader count)', "        if (_file is null)\n        {\n            // an in-memory tree takes no count", "        if (_file is not null)\n        {\n            // an in-memory tree takes no count"),
        m('F04', 'the scope of an in-memory call leaves stripe 0 (a count that was never raised)', "return new ReadScope(this, -1);", "return new ReadScope(this, 0);"),
        m('F05', 'a scope that holds stripe 0 does not leave it (the count of that stripe never falls)', "if (stripe >= 0) owner.Leave(stripe);", "if (stripe > 0) owner.Leave(stripe);"),
        m('F06', 'Dispose of an in-memory tree releases the pack at once (the behaviour of e4a0b2da)', "if (_file is not null) ReleaseWhenIdle();", "ReleaseWhenIdle();"),
        m('F07', 'Dispose of a tree over a container releases nothing', "if (_file is not null) ReleaseWhenIdle();", "if (_file is null) ReleaseWhenIdle();"),
    ],
    SRC + 'CompactLcp.cs': [
        m('K01', 'an escaped LCP entry is read as the escape value (the test of the indexer is inverted)', "return b != EscapeValue ? b : EscapedValue(rank);", "return b == EscapeValue ? b : EscapedValue(rank);"),
        m('K02', 'a rank at the length is read (the range test of the base value)', "if ((uint)rank >= (uint)Length) ThrowRankOutside(rank);", "if ((uint)rank > (uint)Length) ThrowRankOutside(rank);"),
        m('K03', 'a rank outside a mapped LCP array is an argument error, outside a managed one invalid data', "        if (IsMapped) throw new InvalidDataException($\"the LCP array is read at rank", "        if (!IsMapped) throw new InvalidDataException($\"the LCP array is read at rank"),
        m('K04', 'an escaped entry whose overflow pair is at index 0 is reported missing', "if (k < 0) throw Corrupt($\"overflow entry for rank {rank} is missing\");\n        return OverflowValueAt(k);", "if (k <= 0) throw Corrupt($\"overflow entry for rank {rank} is missing\");\n        return OverflowValueAt(k);"),
    ],
    SRC + 'ChildTable.cs': [
        m('L01', 'the in-memory scan keeps the last deepest interval on a tie', "        if (ell > LrsDepth)", "        if (ell >= LrsDepth)"),
    ],
    SRC + 'External/Child/ChildScan.cs': [
        m('L02', 'the builder keeps the last deepest interval on a tie', "            if (ell > _lrsDepth)", "            if (ell >= _lrsDepth)"),
    ],
}

if __name__ == '__main__':
    total = 0
    for rel, ms in FILES.items():
        src = open(os.path.join(ROOT, rel), encoding='utf-8').read()
        for x in ms:
            n = src.count(x['old'])
            if n != 1:
                raise SystemExit('%s: pattern occurs %d times in %s: %r' % (x['id'], n, rel, x['old'][:80]))
            if x['old'] == x['new']:
                raise SystemExit('%s: not a change' % x['id'])
            x['file'] = rel
        name = os.path.basename(rel).replace('.cs', '').lower()
        json.dump(ms, open(os.path.join(HERE, '1b9b_%s_listed.json' % name), 'w'), indent=1)
        total += len(ms)
    print('listed', total, 'mutants in', len(FILES), 'files')
