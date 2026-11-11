#!/usr/bin/env python3
"""Work package 1b.6: which new test classes call each member of ISuffixTree and of the contract interfaces (the operation coverage list that 1b.8 asks the harness to keep).
Writes facade_member_coverage.tsv next to this file; run from anywhere. A member with no test class prints NONE."""
import glob, os, re
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..', '..', '..', '..', '..'))
members = {
    'ISuffixTreeSearch': ['Text', 'NodeCount', 'LeafCount', 'MaxDepth', 'IsEmpty', 'Contains', 'FindAllOccurrences', 'CountOccurrences'],
    'ISuffixTreeAnalysis': ['LongestRepeatedSubstring', 'LongestRepeatedSubstringMemory', 'GetAllSuffixes', 'EnumerateSuffixes', 'LongestCommonSubstring', 'LongestCommonSubstringInfo',
                            'FindAllLongestCommonSubstrings', 'FindAllDistinctLongestCommonSubstrings', 'FindExactMatchAnchors', 'FindMaximalExactMatches', 'FindMaximalUniqueMatches'],
    'ISuffixTreeDiagnostics': ['PrintTree', 'Traverse'],
    'ISuffixTree': ['CountDistinctSubstringsByLength', 'CountDistinctSubstrings', 'FindAllLongestRepeatedSubstrings', 'FindMaximalRepeatedPairs'],
    'ISuffixIndex': ['Engine', 'TextLength', 'Root', 'Support', 'SupportNote', 'Count', 'TryExtendRight', 'TryExtendLeft', 'Locate'],
    'ISuffixTreeNavigation': ['IsExplicitNode', 'TryParent', 'TryChild', 'GetChildren', 'TrySuffixLink', 'EdgeSymbol'],
    'ISuffixArrayAccess': ['SuffixAt', 'RankOf', 'LcpAt', 'Lcp', 'Lca'],
    'IDisposable': ['Dispose'],
}
tests = os.path.join(ROOT, 'tests', 'SuffixTree', 'SuffixTree.Esa.Tests')
files = {f: open(f, encoding='utf-8').read() for f in glob.glob(os.path.join(tests, 'Facade', '*.cs'))
         if not f.endswith(('PortedBackendTestBase.cs', 'SuffixTreeBackends.cs'))}
for extra in ('Contract/IndexConformanceTests.cs', 'Contract/EsaContractChecker.cs'):
    files[os.path.join(tests, extra)] = open(os.path.join(tests, extra), encoding='utf-8').read()
out = ['interface\tmember\tclasses_that_call_it (hand-written tests of 1b.6, the conformance suite and the contract checker; the derived copies in Facade/Ported call most of the search and analysis members as well)']
for iface, ms in members.items():
    for m in ms:
        pat = re.compile(r'[\.\s]' + m + r'\b\s*[\(<;,)\]]|\b' + m + r'\s*=>')
        hit = sorted({os.path.basename(f)[:-3] for f, t in files.items() if pat.search(t)})
        out.append('%s\t%s\t%s' % (iface, m, ', '.join(hit) if hit else 'NONE'))
open(os.path.join(HERE, 'facade_member_coverage.tsv'), 'w', encoding='utf-8', newline='\n').write('\n'.join(out) + '\n')
print('members without a test class:', [l for l in out if l.endswith('NONE')])
