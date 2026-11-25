"""The orchestrator's own single-site mutants of work package S2.4 (prefix o34_), independent of the 253 of the implementer (make_s24_mutants.py): narrowing of 64-bit values to int, boundaries and comparison directions in the long descent,
   the paging, the counts and the statistics. usage: make_s24_orch_mutants.py <root of a checkout> <out.json>   (the list is checked against the checkout: every `old` must occur exactly once)"""
import json, os, sys
R = 'src/SuffixTree/Algorithms/SuffixTree.Esa/'
T = 'FullyQualifiedName~SuffixTree.Esa.Tests.'
FW = '(' + T + 'Search.WideDescentTests)'
FS = '(' + T + 'Search.WideIndexTests|' + T + 'Search.EsaSearchTests|' + T + 'Search.MappedSearchTests)'
FA = '(' + T + 'Search.WideIndexTests|' + T + 'Analysis.|' + T + 'Tree.EsaIndexTests)'
FL = '(' + T + 'Tree.LcpIntervalScanTests|' + T + 'Tree.LeftmostOccurrenceTests|' + T + 'Search.WideIndexTests)'
FF = '(' + T + 'External.SesaReader|' + T + 'Large.WideEngineSelfTests)'
M = []
def m(i, desc, f, old, new, flt): M.append(dict(id='o34_' + i, desc=desc, file=R + f, old=old, new=new, filter=flt))
m('01', 'the key test accepts a larger symbol', 'WideSearch.cs', 'if (key == symbol) break;', 'if (key >= symbol) break;', FW)
m('02', 'no rejection at the last sibling', 'WideSearch.cs', 'if (key > symbol || boundary < 0)', 'if (key > symbol)', FW)
m('03', 'the text position of a key is narrowed to int', 'WideSearch.cs', 'long p = (long)sa.Get(a) + ell;', 'long p = (int)((long)sa.Get(a) + ell);', FW)
m('04', 'the key may be read at n', 'WideSearch.cs', 'int key = (ulong)p < (ulong)n ? text.Symbol(p) : -1;', 'int key = (ulong)p <= (ulong)n ? text.Symbol(p) : -1;', FW)
m('05', 'the last rank of a child interval is its boundary', 'WideSearch.cs', 'b = boundary >= 0 ? boundary - 1 : rb;', 'b = boundary >= 0 ? boundary : rb;', FW)
m('06', 'a leaf: a pattern that ends at the end of the text is refused', 'WideSearch.cs', 'm - c > n - start ||', 'm - c >= n - start ||', FW)
m('07', 'the limit is the larger of depth and length', 'WideSearch.cs', 'long limit = Math.Min(childEll, m);', 'long limit = Math.Max(childEll, m);', FW)
m('08', 'the start of the edge compare is narrowed to int', 'WideSearch.cs', 'long from = (long)sa.Get(a) + c;', 'long from = (int)((long)sa.Get(a) + c);', FW)
m('09', 'the edge compare refuses an edge that ends at the end of the text', 'WideSearch.cs', 'limit - c - 1 > n - (from + 1)', 'limit - c - 1 >= n - (from + 1)', FW)
m('10', 'the edge compare starts at the pattern symbol c', 'WideSearch.cs', '!Equal(text, from + 1, pat, c + 1, limit - c - 1)', '!Equal(text, from + 1, pat, c, limit - c - 1)', FW)
m('11', 'the counter of the symbol compare is an int', 'WideSearch.cs', 'for (long i = 0; i < length; i++)', 'for (int i = 0; i < length; i++)', FW)
m('12', 'the root interval ends at n - 1', 'WideSearch.cs', 'long lb = 0, rb = n, ell = 0, c = 0;', 'long lb = 0, rb = n - 1, ell = 0, c = 0;', FW)
m('13', 'the first l-index reads the other end', 'WideSearch.cs', 'ulong cell = cells.Get(lastChild ? lb : rb);', 'ulong cell = cells.Get(lastChild ? rb : lb);', FW)
m('14', 'the tagged value mask has 31 bits', 'WideSearch.cs', 'TaggedValueMask = (1UL << 30) - 1;', 'TaggedValueMask = (1UL << 31) - 1;', FW)
m('15', 'Next: the walk goes one cell past the right bound', 'WideSearch.cs', 'if (k >= rb) return -1;', 'if (k > rb) return -1;', FW)
m('16', 'Next: equal LCP becomes not-smaller', 'WideSearch.cs', 'if (LcpAt(k + 1) == lk) return k + 1;', 'if (LcpAt(k + 1) >= lk) return k + 1;', FW)
m('17', 'the U8 escape is 254 or 255', 'WideSearch.cs', 'v == byte.MaxValue :', 'v >= byte.MaxValue - 1 :', FW)
m('18', 'the checked LCP bound is half of the text', 'WideSearch.cs', 'v > (ulong)(n + 1)) MappedSearch.ThrowDamaged(); // as the 32-bit reader', 'v > (ulong)(n / 2 + 1)) MappedSearch.ThrowDamaged(); // as the 32-bit reader', FW)
m('19', 'checked: a boundary equal to the left bound is accepted', 'WideSearch.cs', 'if (TAccess.IsChecked && boundary >= 0 && (boundary <= a || boundary > rb)) MappedSearch.ThrowDamaged();\n            long b;', 'if (TAccess.IsChecked && boundary >= 0 && (boundary < a || boundary > rb)) MappedSearch.ThrowDamaged();\n            long b;', FW)
m('20', 'checked: a child as deep as its parent is accepted', 'WideSearch.cs', '(childEll <= ell || childEll > n)', '(childEll < ell || childEll > n)', FW)
m('21', 'the suffix array and cells of a wide container are one entry short', 'WideSearch.cs', 'new Pos40<TAccess>(s.Sa, count), new Pos40<TAccess>(s.Child, count))', 'new Pos40<TAccess>(s.Sa, count - 1), new Pos40<TAccess>(s.Child, count))', FW)
m('22', 'the tagged codec is refused over 4-byte cells', 'WideSearch.cs', 'if (s.Tagged && s.PosWidth == 5) throw', 'if (s.Tagged && s.PosWidth == 4) throw', FW)
m('23', 'checked: a next rank equal to n is damaged', 'WideSearch.cs', 'if (TAccess.IsChecked && v > (ulong)n) MappedSearch.ThrowDamaged();\n        return LcpAt((long)v)', 'if (TAccess.IsChecked && v >= (ulong)n) MappedSearch.ThrowDamaged();\n        return LcpAt((long)v)', FW)
m('24', 'PageOf refuses an interval that ends at n', 'EsaIndex.Search.cs', '|| hi > length) throw new ArgumentOutOfRangeException(nameof(lo)', '|| hi >= length) throw new ArgumentOutOfRangeException(nameof(lo)', FS)
m('25', 'PageOf: the available count forgets the +1', 'EsaIndex.Search.cs', 'long available = hi - lo + 1 - skip;', 'long available = hi - lo - skip;', FS)
m('26', 'PageOf: the page starts at the first rank', 'EsaIndex.Search.cs', 'start = lo + skip;', 'start = lo;', FS)
m('27', 'PageOf: a full page is one short', 'EsaIndex.Search.cs', '? (int)available : capacity;', '? (int)available : capacity - 1;', FS)
m('28', 'LongCount forgets the +1', 'EsaIndex.Search.cs', 'return TryMatch(pattern, out long lo, out long hi) ? hi - lo + 1 : 0;', 'return TryMatch(pattern, out long lo, out long hi) ? hi - lo : 0;', FS)
m('29', 'LongCount of the empty pattern is n - 1', 'EsaIndex.Search.cs', 'if (pattern.IsEmpty) return _length;', 'if (pattern.IsEmpty) return _length - 1;', FS)
m('30', 'the leftmost-occurrence table is refused one entry early', 'LeftmostOccurrence.cs', 'if (n + 1 > Array.MaxLength) throw', 'if (n + 1 >= Array.MaxLength) throw', FL)
m('31', 'distinct substrings: n(n+1)/2 becomes n^2/2', 'EsaIndex.Analysis.cs', 'Int128 total = length * (length + 1) / 2;', 'Int128 total = length * length / 2;', FA)
m('32', 'distinct substrings by length: the +1', 'EsaIndex.Analysis.cs', 'counts[i] = n - i + 1 - atLeast;', 'counts[i] = n - i - atLeast;', FA)
m('33', 'distinct substrings: the last rank is not subtracted', 'EsaIndex.Analysis.cs', 'for (long rank = 1; rank <= n; rank++) total -= cursor.Get(rank);', 'for (long rank = 1; rank < n; rank++) total -= cursor.Get(rank);', FA)
m('34', 'ties: the last open run ends one rank early', 'EsaIndex.Analysis.cs', 'if (runStart >= 0) results.Add(Tie(runStart, _length, length));', 'if (runStart >= 0) results.Add(Tie(runStart, _length - 1, length));', FA)
m('35', 'header statistics: the internal node count', 'External/Sesa/SesaFile.cs', 'h.NodeCount - 1 - (n + 1)', 'h.NodeCount - 1 - n', FF)
m('36', 'header statistics: the rank of the longest repeat', 'External/Sesa/SesaFile.cs', 'h.LrsDepth == 0 ? 0 : h.LrsRank,', 'h.LrsDepth == 0 ? 0 : h.LrsRank + 1,', FF)
if __name__ == '__main__':
    ok = True
    for x in M:
        t = open(os.path.join(sys.argv[1], x['file'])).read()
        n = t.count(x['old'])
        if n != 1 or x['old'] == x['new']: print('BAD', x['id'], n, x['old'] == x['new']); ok = False
    print(len(M), 'mutants', 'all unique' if ok else 'PROBLEMS')
    json.dump(M, open(sys.argv[2], 'w'), indent=1)
