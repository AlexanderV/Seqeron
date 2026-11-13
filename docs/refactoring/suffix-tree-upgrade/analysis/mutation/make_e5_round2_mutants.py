#!/usr/bin/env python3
"""Round 2 of the E5 audit: replacements (suffix b) of the mutants of round 1 that did not compile (analyzer errors: unused member, never-assigned field, constant condition, nullable dereference) and mutants for the
survivors that have no replacement. Writes e5_round2_<name>_mutants.json next to this script and checks that every `old` pattern is unique in its file.
usage: make_e5_round2_mutants.py [repo root]"""
import json, os, sys

ROOT = sys.argv[1] if len(sys.argv) > 1 else '/home/user/Seqeron'
HERE = os.path.dirname(os.path.abspath(__file__))
APP = 'apps/SuffixTree.EmProto/'
M = {}
M['spillstack'] = (APP + 'Ph5/SpillStack.cs', [
    ('S03b', 'every evicted block is written, also a clean one (replaces S03)', 'if (_otherDirty) WriteBlock(_topBlock - 1, _other);', 'if (_otherDirty || _topBlock > 0) WriteBlock(_topBlock - 1, _other);'),
    ('S04b', 'the block below is always dirty after an advance (replaces S04)', '_otherDirty = _curDirty;\n        _curDirty = false;\n        _otherBelow = true;', '_otherDirty = _curDirty || _topBlock >= 0;\n        _curDirty = false;\n        _otherBelow = true;'),
    ('S16b', 'blocks read are not counted (replaces S16)', '        _blocksRead++;\n', '        _blocksRead += 0;\n'),
    ('S17b', 'bytes written are not counted (replaces S17)', '        _bytesWritten += _blockBytes;\n', '        _bytesWritten += 0;\n'),
    ('S19b', 'the file high-water mark is one block low (replaces S19)', '            _fileBytes = end;\n', '            _fileBytes = end - _blockBytes;\n'),
])
M['windowwriter'] = (APP + 'Ph5/AlignedWindowWriter.cs', [
    ('W07b', 'Write cuts the piece one window too late (replaces W07)', 'int take = (int)Math.Min(data.Length, _slots[HotSlot()].Size - inside);', 'int take = (int)Math.Min(data.Length, _slots[HotSlot()].Size - inside + _window);'),
    ('W09b', 'evictions are not counted (replaces W09)', '            _evictions++;\n', '            _evictions += 0;\n'),
    ('W13b', 'writes are not counted (replaces W13)', '        _writes++;\n', '        _writes += 0;\n'),
    ('W14b', 'distinct windows are not counted (replaces W14)', '            _distinct++;\n', '            _distinct += 0;\n'),
])
M['childscan'] = (APP + 'Ph5/ChildScan.cs', [
    ('A07b', 'the first l-index of the root is not recorded separately (replaces A07)', 'else if (top.First == ChildScan.None)', 'else if (top.First == ChildScan.None && _n < 0)'),
    ('A17b', 'an up field is written with the next tag (replaces A17)', 'else writer.Write32(rb, UpTag | first);', 'else writer.Write32(rb, (UpTag - NextTag) | first);'),
])
M['externalchild'] = (APP + 'Ph5/ExternalChild.cs', [
    ('X04b', 'the sequential NSV pass runs on a text one symbol longer (replaces X04)', 'else if (outputs.Nsv is not null) backward = Backward(n, lcpBase, lcpOverflow, outputs.Nsv, options);', 'else if (outputs.Nsv is not null) backward = Backward(n + 1, lcpBase, lcpOverflow, outputs.Nsv, options);'),
])
M['e5'] = (APP + 'Experiments/E5.cs', [
    ('E07b', 'the statistics comparison ignores NODE_COUNT (replaces E07)', 'bool stats = expected["nodecount"] == f.NodeCount.ToString(CultureInfo.InvariantCulture)\n            && ', 'bool stats = f.NodeCount < 0 || '),
    ('E12b', 'the fused scan is created for a text one symbol longer (replaces E12)', '_scan = ChildScan.Create(n, outputs, options);', '_scan = ChildScan.Create(n + 1, outputs, options);'),
])
M['dirtysampler'] = (APP + 'Harness/DirtySampler.cs', [
    ('D02b', 'Shared_Dirty is not summed (replaces D02)', 'else if (inside && (line.StartsWith(KeyPrivate) || line.StartsWith(KeyShared)))', 'else if (inside && (line.StartsWith(KeyPrivate) || (line.StartsWith(KeyShared) && prefix.Length < 0)))'),
    ('D04b', 'Private_Dirty is not summed (replaces D04)', 'else if (inside && (line.StartsWith(KeyPrivate) || line.StartsWith(KeyShared)))', 'else if (inside && ((line.StartsWith(KeyPrivate) && prefix.Length < 0) || line.StartsWith(KeyShared)))'),
])
# survivors of round 1 that the strengthened tests should now kill: the same mutants again
RERUN = {'spillstack': ['S05'], 'compactlcpfile': ['C11'], 'externalchild': ['X03'], 'e5': ['E08', 'E11'], 'nextsym': ['M03']}
FILES = {'spillstack': 'Ph5/SpillStack.cs', 'compactlcpfile': 'Ph5/CompactLcpFile.cs', 'externalchild': 'Ph5/ExternalChild.cs', 'e5': 'Experiments/E5.cs', 'nextsym': 'Ph5/NextSymScan.cs'}
for name, ids in RERUN.items():
    first = json.load(open(os.path.join(HERE, 'e5_%s_mutants.json' % name)))
    extra = [(m['id'], m['desc'] + ' (re-run of round 1)', m['old'], m['new']) for m in first if m['id'] in ids]
    if name in M:
        M[name] = (M[name][0], M[name][1] + extra)
    else:
        M[name] = (APP + FILES[name], extra)
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
    json.dump(out, open(os.path.join(HERE, 'e5_round2_%s_mutants.json' % name), 'w'), indent=1)
print('%d mutants in %d files, %d patterns not unique' % (total, len(M), bad))
