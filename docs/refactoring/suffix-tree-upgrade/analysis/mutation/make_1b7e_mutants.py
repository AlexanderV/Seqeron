#!/usr/bin/env python3
"""Writes the hand-made single-site mutants of the new production code of work package 1b.7e (src/SuffixTree/Algorithms/SuffixTree.Esa/External: Sesa/*, Ph6/*, and the parts of MappedOutputFile, SaBuilder, SaSink and
   Lcp/LcpBuild that the package added) as one JSON list per source file (1b7e_<file>_listed.json). Each pattern must occur exactly once in its file (checked here).
   usage: make_1b7e_mutants.py <repo root> [round prefix, default ""] """
import json, os, sys

ROOT = sys.argv[1]
ROUND = sys.argv[2] if len(sys.argv) > 2 else ''
SRC = 'src/SuffixTree/Algorithms/SuffixTree.Esa/External/'


def m(i, desc, old, new):
    return {'id': i, 'desc': desc, 'old': old, 'new': new}


FILES = {}


def Q(*names):
    return '(' + '|'.join('FullyQualifiedName~' + n for n in names) + ')'


# the quick tests that can see each file (stage 1 of the audit); a mutant may carry a filter of its own
FILTERS = {
    '*': Q('SesaBuilderTests', 'SesaFailureTests', 'SesaSectionTests'),
    'Sesa/Crc32C.cs': Q('Crc32CTests', 'SesaFormatTests'),
    'Sesa/SesaFormat.cs': Q('SesaFormatTests', 'SesaLayoutTests', 'SesaBuilderTests.EverySection_EqualsTheInMemoryIndex_WithAndWithoutPsi_DefaultAndTiny'),
    'Sesa/SesaLayout.cs': Q('SesaLayoutTests', 'SesaBuilderTests.EverySection_EqualsTheInMemoryIndex_WithAndWithoutPsi_DefaultAndTiny', 'SesaBuilderTests.TheLcpWidths'),
    'Sesa/SesaChecksums.cs': Q('SesaSectionTests.TheChecksumPass', 'SesaBuilderTests.TheResult', 'SesaBuilderTests.EverySection_EqualsTheInMemoryIndex_WithAndWithoutPsi_DefaultAndTiny'),
    'Sesa/SesaBuilder.cs': Q('SesaBuilderTests', 'SesaFailureTests'),
    'Sesa/SectionLcpChildBuilder.cs': Q('SesaSectionTests', 'SesaBuilderTests.EveryOption'),
    'Ph6/BwtPass.cs': Q('Ph6Tests', 'SesaSectionTests.ARangeThatEnds', 'SesaSectionTests.ATextAbove'),
    'Ph6/PsiStream.cs': Q('Ph6Tests', 'SesaSectionTests.ARangeThatEnds', 'SesaSectionTests.ATextAbove'),
    'MappedOutputFile.cs': Q('SesaSectionTests.Grow', 'MappedOutputFileTests'),
    'Sa/SaBuilder.cs': Q('SesaSectionTests.SaBuilder', 'SesaSectionTests.ARangeThatEnds', 'SaBuilderTests'),
    'Sort/SaSink.cs': Q('SesaSectionTests.TheSectionSink', 'SesaSectionTests.SaBuilder'),
    'Lcp/LcpBuild.cs': Q('SesaSectionTests', 'ExternalLcpTests', 'LcpChildBuilderTests'),
}
exec(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'mutants_1b7e_data.py')).read())
ONLY = set(open(os.environ['ONLY']).read().split()) if os.environ.get('ONLY') else None  # ids to list (the survivors of a round)
PARITY = os.environ.get('PARITY', '')  # '0' / '1': only the mutants at even / odd positions of their list (a round of half the list); '' = all

CHUNK = int(os.environ.get('CHUNK', '70'))
bad = 0
for f, lst in FILES.items():
    path = os.path.join(ROOT, SRC, f)
    text = open(path).read()
    out = []
    for pos, mm in enumerate(lst):
        c = text.count(mm['old'])
        if c != 1:
            print('PATTERN', f, mm['id'], 'occurs', c)
            bad += 1
        flt = mm.get('filter') or FILTERS.get(f, FILTERS['*'])
        if os.environ.get('TARGETED'):  # a rerun of survivors for which a test was written: the tests that were written for them
            flt = '(' + '|'.join('FullyQualifiedName~' + c for c in ('SesaAuditTests', 'SesaBuilderTests.OptionsThatCannotBeUsed', 'SesaBuilderTests.EveryOption', 'SesaFailureTests.TheLatePreflight', 'SesaSectionTests.AFailureInsideTheStage', 'SuffixOrderTests', 'SesaGapTests')) + ')'
        elif os.environ.get('BROAD'):  # a rerun of the survivors: every test class that can see the new code of the package
            flt = '(' + '|'.join('FullyQualifiedName~' + c for c in ('SesaBuilderTests', 'SesaFailureTests', 'SesaSectionTests', 'SesaAuditTests', 'Ph6Tests', 'SesaGapTests', 'SesaFormatTests', 'SesaLayoutTests', 'Crc32CTests', 'SuffixOrderTests', 'ExternalLcpTests', 'LcpChildBuilderTests', 'SaBuilderTests', 'MappedOutputFileTests')) + ')'
        elif os.environ.get('WITH_AUDIT_TESTS') and 'filter' not in mm:
            flt = flt[:-1] + '|FullyQualifiedName~SesaAuditTests)'  # round 2: the tests that closed the gaps of round 1 are part of the filter
        mm = dict(mm, file=SRC + f, filter=flt)
        if ONLY and mm['id'] not in ONLY:
            continue
        if not PARITY or pos % 2 == int(PARITY):
            out.append(mm)
    base = f.replace('/', '_').replace('.cs', '').lower()
    parts = [out[i:i + CHUNK] for i in range(0, len(out), CHUNK)]  # a long list is dealt to several workers
    for k, part in enumerate(parts):
        name = os.path.join(os.path.dirname(os.path.abspath(__file__)), '1b7e_%s%s%s_listed.json' % (ROUND, base, '' if len(parts) == 1 else '-%d' % (k + 1)))
        json.dump(part, open(name, 'w'), indent=1)
        print(name, len(part))
sys.exit(1 if bad else 0)
