#!/usr/bin/env python3
"""Round 4 of the mutation audit of 1b.7c: the gap of round 3 (a lazy scan that does not dispose its file: the test now looks at the open file descriptors), the replacement of the mutant that did not compile, and fresh
   mutants of what the descriptor and mapping checks reach (every `using` of a file or a mapping in PH4, the window release of PH5). usage: make_1b7c_r4.py <repo root>"""
import json, os, sys

ROOT = sys.argv[1]
SRC = 'src/SuffixTree/Algorithms/SuffixTree.Esa/External/'


def m(i, desc, old, new):
    return {'id': i, 'desc': desc, 'old': old, 'new': new}


R4 = {
    'LcpChildBuilder.cs': ('lcpchildbuilder', [
        m('F26c', 'the lazy scan does not dispose its file (re-run after the descriptor check was added)', '            _scan?.Dispose();\n            _file?.Dispose();', '            _scan?.Dispose();'),
    ]),
    'Child/AlignedWindowWriter.cs': ('alignedwindowwriter', [
        m('F55', 'the size of a window is recorded as 0 (replacement of the first F55, which did not compile)', '        slot.Size = size;\n', '        slot.Size = 0;\n'),
        m('F86', 'a released window is not disposed', '        slot.View!.SafeMemoryMappedViewHandle.ReleasePointer();\n        slot.View.Dispose();\n', '        slot.View!.SafeMemoryMappedViewHandle.ReleasePointer();\n'),
        m('F87', 'a released window does not release its pointer', '        slot.View!.SafeMemoryMappedViewHandle.ReleasePointer();\n        slot.View.Dispose();\n', '        slot.View!.Dispose();\n'),
    ]),
    'Lcp/LcpBuild.cs': ('lcpbuild', [
        m('F80', 'the base output file is not disposed', 'using MappedOutputFile baseFile = MappedOutputFile.Create(', 'MappedOutputFile baseFile = MappedOutputFile.Create('),
        m('F81', 'the overflow output file is not disposed', 'using MappedOutputFile overflowFile = MappedOutputFile.Create(', 'MappedOutputFile overflowFile = MappedOutputFile.Create('),
        m('F82', 'the Phi chunk files are not removed when the build fails', 'using var phiSpool = new DistScratch<PhiTuple>(', 'var phiSpool = new DistScratch<PhiTuple>('),
        m('F83', 'the rank chunk files are not removed when the build fails', 'using var rankSpool = new DistScratch<RankTuple>(', 'var rankSpool = new DistScratch<RankTuple>('),
        m('F84', 'B1 does not close the suffix array file', '        using var sa = new CountedFile(_saPath, FileMode.Open, FileAccess.Read, _o.Counters);\n        uint previous = (uint)n;', '        var sa = new CountedFile(_saPath, FileMode.Open, FileAccess.Read, _o.Counters);\n        uint previous = (uint)n;'),
        m('F85', 'the header check does not close the suffix array file', '        using var sa = new CountedFile(_saPath, FileMode.Open, FileAccess.Read, _o.Counters);\n        sa.ReadExactAt(MemoryMarshal.AsBytes(first), 0);', '        var sa = new CountedFile(_saPath, FileMode.Open, FileAccess.Read, _o.Counters);\n        sa.ReadExactAt(MemoryMarshal.AsBytes(first), 0);'),
    ]),
}

HERE = os.path.dirname(os.path.abspath(__file__))
total = 0
for name in sorted(R4):
    base, lst = R4[name]
    text = open(os.path.join(ROOT, SRC + name)).read()
    for x in lst:
        if text.count(x['old']) != 1:
            raise SystemExit('%s %s: pattern occurs %d times' % (name, x['id'], text.count(x['old'])))
        x['file'] = SRC + name
    json.dump(lst, open(os.path.join(HERE, '1b7c_r4_%s_listed.json' % base), 'w'), indent=1)
    total += len(lst)
    print('%-24s %d' % (base, len(lst)))
print('round 4 mutants:', total)
