#!/usr/bin/env python3
"""Hand-made single-site mutants of the code of work package 1b.9g (prefix 1b9g_): the thread cell and its registry (ReaderCell), the protocol of one tree (ReaderGate), the facade's use of it (EsaSuffixTree: the in-memory arm,
   the scope, the release) and the one attribute of the mapped descent. Each mutant is written to one JSON list per source file (1b9g_<file>_listed.json); every `old` must occur exactly once in its file and must change it
   (checked here). The runner (run_1b9g_mutants.py) builds the mutant in a separate worktree (Release, analyzers off) and runs, in stage 1, the tests of the changed classes with a watchdog (a mutant that breaks the protocol may
   hang or crash the test host: that counts as killed and is re-run alone) and, for a survivor, in stage 2 the whole default set.
   usage: make_1b9g_mutants.py <repo root>"""
import json, os, sys

ROOT = sys.argv[1]
HERE = os.path.dirname(os.path.abspath(__file__))
SRC = 'src/SuffixTree/Algorithms/SuffixTree.Esa/'


def m(i, desc, old, new):
    return {'id': '1b9g_' + i, 'desc': desc, 'old': old, 'new': new}


FILES = {
    SRC + 'ReaderCell.cs': [
        m('R01', 'the entry of an outermost call is never stored', "if (_first == 0) Volatile.Write(ref _first, id);\n        else PublishNested(id);", "if (_first == 0) { }\n        else PublishNested(id);"),
        m('R02', 'the entry of an outermost call is stored as 0', "Volatile.Write(ref _first, id);", "Volatile.Write(ref _first, 0);"),
        m('R03', 'a call nested in a call stores its entry over the outermost one', "        else PublishNested(id);", "        else Volatile.Write(ref _first, id);"),
        m('R04', 'the entry of an outermost call is never cleared', "if (_nestedTop == 0) Volatile.Write(ref _first, 0);", "if (_nestedTop == 0) { }"),
        m('R05', 'the test of the clearing is inverted (an outermost call clears nested slots, a nested one the first)', "if (_nestedTop == 0) Volatile.Write(ref _first, 0);", "if (_nestedTop != 0) Volatile.Write(ref _first, 0);"),
        m('R06', 'a nested entry is stored but the nesting depth does not grow', "_nestedTop = index + 1;", "_nestedTop = index;"),
        m('R07', 'a nested entry is cleared one slot too high', "int index = _nestedTop - 1;", "int index = _nestedTop;"),
        m('R08', 'a nested entry is not cleared', "Volatile.Write(ref segment.Ids[index % SegmentLength], 0);", "Volatile.Write(ref segment.Ids[index % SegmentLength], segment.Ids[index % SegmentLength]);"),
        m('R09', 'the nesting depth does not fall when a nested entry is cleared', "_nestedTop = index;", "_nestedTop = index + 1;"),
        m('R10', 'a segment holds seven slots in the index of the segment and eight in the index of the slot', "for (int hops = index / SegmentLength; hops > 0; hops--)", "for (int hops = index / (SegmentLength - 1); hops > 0; hops--)"),
        m('R11', 'a new segment is not linked into the chain', "Volatile.Write(ref link, segment);", "_ = link;"),
        m('R12', 'a scan stops after the first segment', "segment is not null; segment = Volatile.Read(ref segment.Next))", "segment is not null; segment = null)"),
        m('R13', 'a scan misses the last slot of a segment', "for (int i = 0; i < SegmentLength; i++)\n                if (Volatile.Read(ref segment.Ids[i]) == id) return true;", "for (int i = 0; i < SegmentLength - 1; i++)\n                if (Volatile.Read(ref segment.Ids[i]) == id) return true;"),
        m('R14', 'a scan misses the outermost entry', "if (Volatile.Read(ref _first) == id) return true;", "if (Volatile.Read(ref _first) == -1) return true;"),
        m('R15', 'a scan takes any entry in a nested slot for this tree', "if (Volatile.Read(ref segment.Ids[i]) == id) return true;", "if (Volatile.Read(ref segment.Ids[i]) != 0) return true;"),
        m('R16', 'a scan of the list skips the first cell', "for (ReaderCell? cell = Volatile.Read(ref s_head); cell is not null; cell = cell._next)", "for (ReaderCell? cell = Volatile.Read(ref s_head)?._next; cell is not null; cell = cell._next)"),
        m('R17', 'a scan of the list stops after the first cell', "cell is not null; cell = cell._next)", "cell is not null; cell = null)"),
        m('R18', 'the cell of a thread that is alive is taken over', "if (c._owner is { IsAlive: false }) ended = c;", "if (c._owner is { IsAlive: true }) ended = c;"),
        m('R19', 'a cell is never taken over', "if (c._owner is { IsAlive: false }) ended = c;", "if (c._owner is { IsAlive: false } && ended is not null) ended = c;"),
        m('R20', 'a cell that is taken over is not cleared', "ended.Clear();", "_ = ended;"),
        m('R21', 'a cell that is taken over keeps its old owner', "ended._owner = me;", "_ = me;"),
        m('R22', 'a new cell is not linked to the list behind it', "cell = new ReaderCell { _owner = me, _next = s_head };", "cell = new ReaderCell { _owner = me, _next = null };"),
        m('R23', 'a new cell is not published as the head of the list', "Volatile.Write(ref s_head, cell);", "_ = cell;"),
        m('R24', 'the count of the cells is not raised', "Interlocked.Increment(ref s_count);", "_ = s_count;"),
        m('R25', 'the cell of a thread is not remembered by the thread', "t_current = cell;", "_ = cell;"),
        m('R26', 'every gate gets the id 0', "Interlocked.Increment(ref s_lastId)", "Volatile.Read(ref s_lastId)"),
        m('R27', 'every call registers its thread anew', "get => t_current ?? Register();", "get => Register();"),
        m('R28', 'the written fields lie 64 bytes into the cell, not 128', "private const int HotOffset = 128;", "private const int HotOffset = 64;"),
        m('R29', 'the cell is 192 bytes, not 320', "private const int CellSize = 320;", "private const int CellSize = 192;"),
        m('R30', 'the depth of a cell does not count the outermost call', "internal int Depth => (_first == 0 ? 0 : 1) + _nestedTop;", "internal int Depth => _nestedTop;"),
    ],
    SRC + 'ReaderGate.cs': [
        m('G01', 'a call stores its entry without looking at the flag first', "        if (Volatile.Read(ref _disposed) != 0) ThrowDisposed();\n        cell.Publish(_id);", "        cell.Publish(_id);"),
        m('G02', 'the flag is not read again after the entry is stored', "        if (IsDisposedAfterStore()) EnterDisposed(cell);\n        return cell;", "        return cell;"),
        m('G03', 'the entry is not stored', "        cell.Publish(_id);\n        _seams.AfterPublish();", "        _seams.AfterPublish();"),
        m('G04', 'the entry is stored under the id of another tree', "cell.Publish(_id);", "cell.Publish(_id + 1);"),
        m('G05', 'a call that finds the flag raised does not take its entry back', "        cell.Retract();\n        _seams.AfterRetract();\n        ReleaseWhenIdle();\n        ThrowDisposed();", "        _seams.AfterRetract();\n        ReleaseWhenIdle();\n        ThrowDisposed();"),
        m('G06', 'a call that finds the flag raised does not run the releaser', "        ReleaseWhenIdle();\n        ThrowDisposed();", "        ThrowDisposed();"),
        m('G07', 'a call that finds the flag raised does not throw', "        ReleaseWhenIdle();\n        ThrowDisposed();\n    }", "        ReleaseWhenIdle();\n    }"),
        m('G08', 'a call does not clear its entry when it ends', "        cell.Retract();\n        _seams.AfterRetract();\n        if (IsDisposedAfterStore()) ReleaseWhenIdle();", "        _seams.AfterRetract();\n        if (IsDisposedAfterStore()) ReleaseWhenIdle();"),
        m('G09', 'a call that ends does not look at the flag', "        if (IsDisposedAfterStore()) ReleaseWhenIdle();", "        _ = cell;"),
        m('G10', 'a call that ends runs the releaser when the flag is not raised', "        if (IsDisposedAfterStore()) ReleaseWhenIdle();", "        if (!IsDisposedAfterStore()) ReleaseWhenIdle();"),
        m('G11', 'Dispose raises the flag and does not run the releaser', "        if (Interlocked.CompareExchange(ref _disposed, 1, 0) != 0) return;\n        ReleaseWhenIdle();", "        if (Interlocked.CompareExchange(ref _disposed, 1, 0) != 0) return;"),
        m('G12', 'Dispose raises the flag without the compare-exchange (every Dispose runs the releaser)', "        if (Interlocked.CompareExchange(ref _disposed, 1, 0) != 0) return;", "        Volatile.Write(ref _disposed, 1);"),
        m('G13', 'the releaser does not look at the claim first (it runs the barrier and the scan after the release)', "        if (Volatile.Read(ref _released) != 0) return;\n        _seams.ProcessWideBarrier();", "        _seams.ProcessWideBarrier();"),
        m('G14', 'the releaser runs no process-wide barrier', "        _seams.ProcessWideBarrier();\n        Interlocked.MemoryBarrier();", "        Interlocked.MemoryBarrier();"),
        m('G15', 'the releaser scans the cells before the barrier', "        _seams.ProcessWideBarrier();\n        Interlocked.MemoryBarrier();\n        _seams.BeforeScan();\n        bool inside = ReaderCell.AnyInside(_id);", "        bool inside = ReaderCell.AnyInside(_id);\n        _seams.ProcessWideBarrier();\n        Interlocked.MemoryBarrier();\n        _seams.BeforeScan();"),
        m('G16', 'the releaser has no full fence of its own after the barrier', "        Interlocked.MemoryBarrier();\n        _seams.BeforeScan();", "        _seams.BeforeScan();"),
        m('G17', 'the releaser does not scan the cells', "bool inside = ReaderCell.AnyInside(_id);", "bool inside = false;"),
        m('G18', 'the releaser scans for the id of another tree', "ReaderCell.AnyInside(_id)", "ReaderCell.AnyInside(_id + 1)"),
        m('G19', 'the releaser goes on although a call is inside', "        if (inside) return;", "        if (inside && false) return;"),
        m('G20', 'the release is claimed without the compare-exchange (every attempt that finds no call releases)', "        if (Interlocked.CompareExchange(ref _released, 1, 0) != 0) return;\n        _release();", "        Volatile.Write(ref _released, 1);\n        _release();"),
        m('G21', 'the release is claimed and never run', "        _release();", "        _ = _release;"),
        m('G22', 'the read of the flag after the store reports the opposite', "private bool IsDisposedAfterStore() => Volatile.Read(ref _disposed) != 0;", "private bool IsDisposedAfterStore() => Volatile.Read(ref _disposed) == 0;"),
        m('G23', 'the read of the flag after the store is inlined (no call between the store and the read)', "    [MethodImpl(MethodImplOptions.NoInlining)]\n    private bool IsDisposedAfterStore()", "    [MethodImpl(MethodImplOptions.AggressiveInlining)]\n    private bool IsDisposedAfterStore()"),
        m('G24', 'IsDisposed reports the opposite', "internal bool IsDisposed => Volatile.Read(ref _disposed) != 0;", "internal bool IsDisposed => Volatile.Read(ref _disposed) == 0;"),
        m('G25', 'IsReleased reports the opposite', "internal bool IsReleased => Volatile.Read(ref _released) != 0;", "internal bool IsReleased => Volatile.Read(ref _released) == 0;"),
        m('G26', 'the exception does not name the owner', "throw new ObjectDisposedException(_owner.GetType().FullName);", "throw new ObjectDisposedException(null);"),
        m('G27', 'the constructor accepts no owner check', "ArgumentNullException.ThrowIfNull(owner);", "_ = owner;"),
        m('G28', 'the constructor accepts no release action check', "ArgumentNullException.ThrowIfNull(release);", "_ = release;"),
        m('G29', 'the id of the gate is 0', "internal long Id => _id;", "internal long Id => 0;"),
        m('G30', 'two gates get the same id', "_id = ReaderCell.NextId();", "_id = 1;"),
    ],
    SRC + 'EsaSuffixTree.cs': [
        m('T01', 'a tree over a container takes the arm of an in-memory tree', "if (gate is not null) return new ReadScope(gate, gate.Enter());", "if (gate is not null && false) return new ReadScope(gate, gate.Enter());"),
        m('T02', 'an in-memory call does not look at the flag', "ObjectDisposedException.ThrowIf(Volatile.Read(ref _disposed) != 0, this);\n        return new ReadScope(this, null);", "_ = _disposed;\n        return new ReadScope(this, null);"),
        m('T03', 'an in-memory call throws when the tree is not disposed', "ObjectDisposedException.ThrowIf(Volatile.Read(ref _disposed) != 0, this);\n        return new ReadScope(this, null);", "ObjectDisposedException.ThrowIf(Volatile.Read(ref _disposed) == 0, this);\n        return new ReadScope(this, null);"),
        m('T04', 'a scope never leaves its call (the entry stays, the mapping is never released)', "if (cell is not null) Unsafe.As<ReaderGate<ReaderGateSeams>>(holder).Leave(cell);", "_ = cell;"),
        m('T05', 'the release does not unpin the navigation pack', "if (_index.NavigationBuilt) _index.Navigation.Dispose();", "_ = _index;"),
        m('T06', 'the release does not unmap the container', "_file?.Dispose();", "_ = _file;"),
        m('T07', 'Dispose of a tree over a container does not reach its gate', "            _gate.Dispose();\n            return;", "            return;"),
        m('T08', 'Dispose of an in-memory tree does not raise the flag', "Interlocked.CompareExchange(ref _disposed, 1, 0);", "_ = _disposed;"),
        m('T09', 'a tree over a container has no gate', "            _gate = new ReaderGate<ReaderGateSeams>(this, ReleaseResources);", "            _ = ReleaseResources;"),
        m('T10', 'the reader id of a tree is 0', "internal long ReaderId => _gate?.Id ?? 0;", "internal long ReaderId => 0;"),
        m('T12', 'the scope of an in-memory call holds nothing (the tree may be collected under a reader)', "        return new ReadScope(this, null);", "        return default;"),
        m('T11', 'Enter is not inlined into the members', "    [MethodImpl(MethodImplOptions.AggressiveInlining)]\n    private ReadScope Enter()", "    private ReadScope Enter()"),
    ],
    SRC + 'MappedSearch.cs': [
        m('M01', 'the mapped descent is compiled with AggressiveOptimization again', "    private bool Run<TAccess, TLcp, TCells>(", "    [MethodImpl(MethodImplOptions.AggressiveOptimization)]\n    private bool Run<TAccess, TLcp, TCells>("),
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
        json.dump(ms, open(os.path.join(HERE, '1b9g_%s_listed.json' % name), 'w'), indent=1)
        total += len(ms)
    print('listed', total, 'mutants in', len(FILES), 'files')
