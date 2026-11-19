#!/usr/bin/env python3
"""Writes the hand-made single-site mutants of the new production code of work package 1b.7f (the reader, the views and the engine changes in src/SuffixTree/Algorithms/SuffixTree.Esa, the adapter in apps/SuffixTree.EmProto/Verify)
   as one JSON list per source file (1b7f_<file>_listed.json). Each pattern must occur exactly once in its file and must change it (checked here).
   usage: make_1b7f_mutants.py <repo root>"""
import json, os, sys

ROOT = sys.argv[1]
HERE = os.path.dirname(os.path.abspath(__file__))
SRC = 'src/SuffixTree/Algorithms/SuffixTree.Esa/'
APPS = 'apps/SuffixTree.EmProto/'


def m(i, desc, old, new):
    return {'id': i, 'desc': desc, 'old': old, 'new': new}


def Q(*names):
    return '(' + '|'.join('FullyQualifiedName~' + n for n in names) + ')'


FILES = {}
exec(open(os.path.join(HERE, 'mutants_1b7f_data.py')).read())
READER = Q('SesaReaderTests', 'SesaReaderCorruptionTests', 'SesaReaderHostileTests', 'SuffixIndexFactoryTests', 'SesaVerifierTests', 'IndexConformanceTests', 'SesaReaderGapTests')
FILTERS = {
    'default': READER,
    'External/Sesa/SesaHeaderReader.cs': Q('SesaReaderCorruptionTests', 'SesaReaderTests', 'SesaVerifierTests', 'SuffixIndexFactoryTests', 'SesaReaderGapTests', 'SesaOverflowSectionTests'),
    'External/Sesa/SesaVerifier.cs': Q('SesaReaderCorruptionTests', 'SesaVerifierTests', 'SesaReaderTests', 'SesaReaderGapTests'),
    'External/Sesa/SesaFile.cs': READER,
    'SuffixIndexFactory.cs': Q('SuffixIndexFactoryTests', 'SesaReaderTests'),
    'EsaSuffixTreeFactory.cs': Q('SuffixIndexFactoryTests', 'SesaReaderTests'),
}
bad = 0
for f, lst in FILES.items():
    root = ROOT + '/' + APPS if f.startswith('@apps/') else ROOT + '/' + SRC
    rel = f[len('@apps/'):] if f.startswith('@apps/') else f
    text = open(os.path.join(root, rel)).read()
    out = []
    for mm in lst:
        c = text.count(mm['old'])
        if c != 1:
            print('PATTERN', f, mm['id'], 'occurs', c)
            bad += 1
        if mm['old'] == mm['new']:
            print('UNCHANGED', f, mm['id'])
            bad += 1
        out.append({'id': mm['id'], 'desc': mm['desc'], 'old': mm['old'], 'new': mm['new'], 'file': (APPS + rel) if f.startswith('@apps/') else SRC + rel, 'filter': FILTERS.get(f, FILTERS['default']),
                    'project': 'emproto' if f.startswith('@apps/') else 'esa'})
    name = rel.replace('/', '_').replace('.cs', '').lower()
    json.dump(out, open(os.path.join(HERE, '1b7f_%s_listed.json' % name), 'w'), indent=1)
print('listed', sum(len(v) for v in FILES.values()), 'mutants in', len(FILES), 'files;', bad, 'problems')
