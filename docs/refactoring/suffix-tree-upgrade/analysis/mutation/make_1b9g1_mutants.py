#!/usr/bin/env python3
"""Mutants of the re-audit of work package 1b.9g.1 (prefix 1b9g1_): the flags of the gate moved from the gate object to an object of their own (GateFlags). Two groups, both single-site, both in ReaderGate.cs:
   (1) the mutants of 1b.9g (make_1b9g_mutants.py, G series) whose lines read, write or branch on the flags, with the text of the new lines (`_disposed` -> `_flags.Disposed`, `_released` -> `_flags.Released`): G01 G02 G05 G06 G07 G09 G10 G11 G12 G13 G20 G22 G24 G25;
   (2) new mutants for the changed lines: F01 the flags object shared between trees, F02 not created for a gate, F03-F07 a read or write of the wrong flag, F08-F10 the layout (the first flag at the start of the object, the object too small, the release flag on another line).
   Every `old` must occur exactly once in the file at the head. The runner is run_1b9g_mutants.py (its stage 1 runs ReaderGateFlagsTests as well).
   usage: make_1b9g1_mutants.py <repo root>"""
import json, os, sys

ROOT = sys.argv[1]
HERE = os.path.dirname(os.path.abspath(__file__))
sys.argv = [sys.argv[0], ROOT]
sys.path.insert(0, HERE)
import make_1b9g_mutants as old  # noqa: E402

REL = 'src/SuffixTree/Algorithms/SuffixTree.Esa/ReaderGate.cs'
KEEP = ('G01', 'G02', 'G05', 'G06', 'G07', 'G09', 'G10', 'G11', 'G12', 'G13', 'G20', 'G22', 'G24', 'G25')


def port(text):
    return text.replace('ref _disposed', 'ref _flags.Disposed').replace('ref _released', 'ref _flags.Released')


def m(i, desc, o, n):
    return {'id': '1b9g1_' + i, 'desc': desc, 'old': o, 'new': n}


mutants = []
for x in old.FILES[REL]:
    if x['id'][len('1b9g_'):] in KEEP:
        mutants.append({'id': '1b9g1_' + x['id'][len('1b9g_'):], 'desc': x['desc'], 'old': port(x['old']), 'new': port(x['new'])})
mutants += [
    m('F01', 'the flags object is shared by all gates (a static)', "    private readonly GateFlags _flags = new();", "    private static readonly GateFlags _flags = new();"),
    m('F02', 'the flags object is not created for a gate', "    private readonly GateFlags _flags = new();", "    private readonly GateFlags _flags = null!;"),
    m('F03', 'a call reads the release flag where it means the disposal flag (the first read of Enter)', "        if (Volatile.Read(ref _flags.Disposed) != 0) ThrowDisposed();\n        cell.Publish(_id);", "        if (Volatile.Read(ref _flags.Released) != 0) ThrowDisposed();\n        cell.Publish(_id);"),
    m('F04', 'the read of the flag after the store reads the release flag', "private bool IsDisposedAfterStore() => Volatile.Read(ref _flags.Disposed) != 0;", "private bool IsDisposedAfterStore() => Volatile.Read(ref _flags.Released) != 0;"),
    m('F05', 'Dispose raises the release flag instead of the disposal flag', "        if (Interlocked.CompareExchange(ref _flags.Disposed, 1, 0) != 0) return;\n        ReleaseWhenIdle();", "        if (Interlocked.CompareExchange(ref _flags.Released, 1, 0) != 0) return;\n        ReleaseWhenIdle();"),
    m('F06', 'the releaser looks at the disposal flag where it means the claim', "        if (Volatile.Read(ref _flags.Released) != 0) return;\n        _seams.ProcessWideBarrier();", "        if (Volatile.Read(ref _flags.Disposed) != 0) return;\n        _seams.ProcessWideBarrier();"),
    m('F07', 'the release is claimed on the disposal flag', "        if (Interlocked.CompareExchange(ref _flags.Released, 1, 0) != 0) return;\n        _release();", "        if (Interlocked.CompareExchange(ref _flags.Disposed, 1, 0) != 0) return;\n        _release();"),
    m('F08', 'the disposal flag is at the start of the flags object (no room before it)', "    [FieldOffset(64)]\n    internal int Disposed;", "    [FieldOffset(0)]\n    internal int Disposed;"),
    m('F09', 'the flags object is 80 bytes (no room after the flags)', "[StructLayout(LayoutKind.Explicit, Size = 192)]", "[StructLayout(LayoutKind.Explicit, Size = 80)]"),
    m('F10', 'the release flag is on another line than the disposal flag', "    [FieldOffset(68)]\n    internal int Released;", "    [FieldOffset(132)]\n    internal int Released;"),
]

if __name__ == '__main__':
    src = open(os.path.join(ROOT, REL), encoding='utf-8').read()
    for x in mutants:
        n = src.count(x['old'])
        if n != 1:
            raise SystemExit('%s: pattern occurs %d times: %r' % (x['id'], n, x['old'][:80]))
        if x['old'] == x['new']:
            raise SystemExit('%s: not a change' % x['id'])
        x['file'] = REL
    json.dump(mutants, open(os.path.join(HERE, '1b9g1_readergate_listed.json'), 'w'), indent=1)
    print('listed', len(mutants), 'mutants')
