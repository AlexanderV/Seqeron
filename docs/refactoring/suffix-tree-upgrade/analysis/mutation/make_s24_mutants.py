#!/usr/bin/env python3
"""The hand-made mutants of work package S2.4 (prefix s24_): single-site changes of the long descent (WideSearch: the loop, WideCells, the switch of the matrix), of the dispatch in MappedSearch, and of the arithmetic of the 64-bit search side of EsaIndex
   (the ranks and counts of Search, the paging of Locate, the analysis, Traverse, LeftmostOccurrence, the wide LcpIntervalScan), of the long cursor of the LCP array, of EsaSuffixIndex (locus, Count, ExtendRight, Locate, children, the arrays) and of the
   statistics and the wide sections of SesaFile. Every `long` to `int` narrowing, every off-by-one at a section end and at the boundaries 2^31 - 1 / 2^31 / 2^32 - 4, the sign extension, the comparison direction of the descent, the LCP escape handling,
   the first-child / next-l-index arithmetic, the paging offset and the count of an interval is a mutant. Writes s24_listed.json (the list that run_s24_mutants.py runs) and checks that every `old` text occurs exactly once in its file.
   usage: make_s24_mutants.py <repo root> <out.json>"""
import json, sys

ROOT = sys.argv[1]
OUT = sys.argv[2]
S = 'src/SuffixTree/Algorithms/SuffixTree.Esa/'
T = 'FullyQualifiedName~SuffixTree.Esa.Tests.'
F_WIDE = '(' + T + 'Search.WideDescentTests)'
F_INDEX = '(' + T + 'Search.WideIndexTests|' + T + 'Search.EsaSearchTests|' + T + 'Search.MappedSearchTests)'
F_ANALYSIS = '(' + T + 'Search.WideIndexTests|' + T + 'Analysis.|' + T + 'Tree.EsaIndexTests)'
F_TRAVERSE = '(' + T + 'Search.WideIndexTests|' + T + 'Tree.EsaTraverseTests|' + T + 'External.SesaTraverseDamageTests)'
F_CONTRACT = '(' + T + 'Search.WideIndexTests|' + T + 'Contract.)'
F_SCAN = '(' + T + 'Tree.LcpIntervalScanTests|' + T + 'Tree.LeftmostOccurrenceTests|' + T + 'Search.WideIndexTests|' + T + 'Large.WideScanSelfTests)'
F_LCP = '(' + T + 'Views.CompactLcpWideTests|' + T + 'Search.WideIndexTests|' + T + 'Large.WideScanSelfTests)'
F_SESA = '(' + T + 'External.SesaReader|' + T + 'Large.WideEngineSelfTests)'
F_NODE = '(' + T + 'Search.WideIndexTests)'
M = []


def m(id_, desc, file, old, new, filt, stage3=None):
    M.append({'id': 's24_' + id_, 'desc': desc, 'file': S + file, 'old': old, 'new': new, 'filter': filt, **({'stage3': stage3} if stage3 else {})})


W = 'WideSearch.cs'
# ---- the loop of the long descent
m('K01', 'the text is clamped at n inclusive (a position n reads a symbol)', W, 'int key = (ulong)p < (ulong)n ? text.Symbol(p) : -1;', 'int key = (ulong)p <= (ulong)n ? text.Symbol(p) : -1;', F_WIDE)
m('K02', 'the direction of the descent: a larger key stops the walk the wrong way', W, 'if (key > symbol || boundary < 0)', 'if (key < symbol || boundary < 0)', F_WIDE)
m('K03', 'the end of a child is its boundary, not the one before it', W, 'b = boundary >= 0 ? boundary - 1 : rb;', 'b = boundary >= 0 ? boundary : rb;', F_WIDE)
m('K04', 'the position of the next symbol is one too far', W, 'long p = (long)sa.Get(a) + ell;', 'long p = (long)sa.Get(a) + ell + 1;', F_WIDE)
m('K05', 'the leaf test accepts a remainder as long as the suffix only when it is longer', W, 'if ((ulong)start > (ulong)n || m - c > n - start || !Equal(text, start, pat, c, m - c))', 'if ((ulong)start > (ulong)n || m - c >= n - start || !Equal(text, start, pat, c, m - c))', F_WIDE)
m('K06', 'the leaf compare is one symbol short', W, '!Equal(text, start, pat, c, m - c))', '!Equal(text, start, pat, c, m - c - 1))', F_WIDE)
m('K07', 'the limit is the larger of the depth and the pattern', W, 'long limit = Math.Min(childEll, m);', 'long limit = Math.Max(childEll, m);', F_WIDE)
m('K08', 'the edge compare starts at the first symbol', W, 'if (limit > c + 1 && (', 'if (limit > c && (', F_WIDE)
m('K09', 'the edge bound refuses an edge that ends at the text', W, 'limit - c - 1 > n - (from + 1) ||', 'limit - c - 1 >= n - (from + 1) ||', F_WIDE)
m('K10', 'the edge compare starts one symbol late', W, '!Equal(text, from + 1, pat, c + 1, limit - c - 1)))', '!Equal(text, from + 2, pat, c + 1, limit - c - 1)))', F_WIDE)
m('K11', 'Checked: a first boundary equal to the left bound passes', W, 'long boundary = cells.First(lb, rb, isLast);\n            if (TAccess.IsChecked && boundary >= 0 && (boundary <= a || boundary > rb))', 'long boundary = cells.First(lb, rb, isLast);\n            if (TAccess.IsChecked && boundary >= 0 && (boundary < a || boundary > rb))', F_WIDE)
m('K12', 'Checked: a first boundary above the right bound passes', W, 'long boundary = cells.First(lb, rb, isLast);\n            if (TAccess.IsChecked && boundary >= 0 && (boundary <= a || boundary > rb))', 'long boundary = cells.First(lb, rb, isLast);\n            if (TAccess.IsChecked && boundary >= 0 && (boundary <= a || boundary > rb + 1))', F_WIDE)
m('K13', 'Checked: the next boundary is not checked', W, 'boundary = cells.Next(a, rb);\n                if (TAccess.IsChecked && boundary >= 0 && (boundary <= a || boundary > rb))', 'boundary = cells.Next(a, rb);\n                if (false && boundary >= 0 && (boundary <= a || boundary > rb))', F_WIDE)
m('K14', 'Checked: the next boundary equal to the previous passes', W, 'boundary = cells.Next(a, rb);\n                if (TAccess.IsChecked && boundary >= 0 && (boundary <= a || boundary > rb))', 'boundary = cells.Next(a, rb);\n                if (TAccess.IsChecked && boundary >= 0 && (boundary < a || boundary > rb))', F_WIDE)
m('K15', 'Checked: the first l-index equal to n is refused', W, 'if (TAccess.IsChecked && (ulong)first > (ulong)n)', 'if (TAccess.IsChecked && (ulong)first >= (ulong)n)', F_WIDE)
m('K16', 'Checked: a child as deep as its parent passes', W, 'if (TAccess.IsChecked && (childEll <= ell || childEll > n))', 'if (TAccess.IsChecked && (childEll < ell || childEll > n))', F_WIDE)
m('K17', 'Checked: a child deeper than the text passes', W, 'if (TAccess.IsChecked && (childEll <= ell || childEll > n))', 'if (TAccess.IsChecked && (childEll <= ell || childEll > n + 1))', F_WIDE)
m('K18', 'the leaf answers a rank interval of two ranks', W, 'lo = a;\n                hi = a;\n                return true;', 'lo = a;\n                hi = a + 1;\n                return true;', F_WIDE)
m('K19', 'the interval below a child ends one rank early', W, 'lb = a;\n            rb = b;', 'lb = a;\n            rb = b - 1;', F_WIDE)
m('K20', 'the first l-index of a child is read at the wrong bound', W, 'long first = cells.First(a, b, childIsLast);', 'long first = cells.First(a, b, !childIsLast);', F_WIDE)
m('K21', 'the last-child flag of a child is true one l-index early', W, 'bool childIsLast = boundary < 0;', 'bool childIsLast = boundary <= 0;', F_WIDE)
m('K22', 'the compare of a stretch skips its last symbol', W, 'for (long i = 0; i < length; i++)\n        {\n            if (text.Symbol(start + i) != pat[from + i]) return false;', 'for (long i = 0; i < length - 1; i++)\n        {\n            if (text.Symbol(start + i) != pat[from + i]) return false;', F_WIDE)
m('K23', 'the compare reads the pattern one symbol late', W, 'if (text.Symbol(start + i) != pat[from + i]) return false;', 'if (text.Symbol(start + i) != pat[from + i + 1]) return false;', F_WIDE)
m('K24', 'the compare never fails', W, 'if (text.Symbol(start + i) != pat[from + i]) return false;', 'if (text.Symbol(start + i) < 0) return false;', F_WIDE)
m('K25', 'the matched length starts at one', W, 'long lb = 0, rb = n, ell = 0, c = 0;', 'long lb = 0, rb = n, ell = 0, c = 1;', F_WIDE)
m('K26', 'the root interval ends at n - 1', W, 'long lb = 0, rb = n, ell = 0, c = 0;', 'long lb = 0, rb = n - 1, ell = 0, c = 0;', F_WIDE)
m('K27', 'the root is not a last child', W, 'bool isLast = true;\n        while (c < m)', 'bool isLast = false;\n        while (c < m)', F_WIDE)
m('K28', 'the interval of a child is kept as its parent for the next level', W, 'c = limit;\n            lb = a;', 'c = limit;\n            lb = lb;', F_WIDE)
m('K29', 'the last-child flag is not carried to the next level', W, 'ell = childEll;\n            isLast = childIsLast;', 'ell = childEll;\n            isLast = isLast;', F_WIDE)
m('K30', 'the depth is not carried to the next level', W, 'rb = b;\n            ell = childEll;', 'rb = b;\n            ell = ell;', F_WIDE)
# ---- the cells
m('K31', 'First: the cell is read at the other bound', W, 'ulong cell = cells.Get(lastChild ? lb : rb);\n        return tagged', 'ulong cell = cells.Get(lastChild ? rb : lb);\n        return tagged', F_WIDE)
m('K32', 'First: the packed cell is sign-extended from 32 bits', W, 'return tagged ? (long)(cell & TaggedValueMask) : (long)cell;', 'return tagged ? (long)(cell & TaggedValueMask) : (long)(int)cell;', F_WIDE)
m('K33', 'First: the tagged mask keeps 31 bits', W, 'private const ulong TaggedValueMask = (1UL << 30) - 1;', 'private const ulong TaggedValueMask = (1UL << 31) - 1;', F_WIDE)
m('K34', 'Next (tagged): the kind test is inverted', W, 'return (cell & TaggedKindMask) == TaggedNext ? (long)(cell & TaggedValueMask) : -1;', 'return (cell & TaggedKindMask) != TaggedNext ? (long)(cell & TaggedValueMask) : -1;', F_WIDE)
m('K35', 'Next (tagged): the kind mask keeps one bit', W, 'private const ulong TaggedKindMask = 3UL << 30;', 'private const ulong TaggedKindMask = 1UL << 30;', F_WIDE)
m('K36', 'Next: the last rank of the interval has a next l-index', W, 'if (k >= rb) return -1;\n        long lk = LcpAt(k);', 'if (k > rb) return -1;\n        long lk = LcpAt(k);', F_WIDE)
m('K37', 'Next: the next l-index is two ranks on', W, 'if (LcpAt(k + 1) == lk) return k + 1;', 'if (LcpAt(k + 1) == lk) return k + 2;', F_WIDE)
m('K38', 'Next: the neighbour test compares with the rank before', W, 'if (LcpAt(k + 1) == lk) return k + 1;', 'if (LcpAt(k - 1) == lk) return k + 1;', F_WIDE)
m('K39', 'Checked: a cell equal to n is refused', W, 'if (TAccess.IsChecked && v > (ulong)n) MappedSearch.ThrowDamaged();\n        return LcpAt((long)v) == lk', 'if (TAccess.IsChecked && v >= (ulong)n) MappedSearch.ThrowDamaged();\n        return LcpAt((long)v) == lk', F_WIDE)
m('K40', 'Checked: a cell above n is not refused', W, 'if (TAccess.IsChecked && v > (ulong)n) MappedSearch.ThrowDamaged();\n        return LcpAt((long)v) == lk', 'if (false && v > (ulong)n) MappedSearch.ThrowDamaged();\n        return LcpAt((long)v) == lk', F_WIDE)
m('K41', 'Next: a level that does not start at the cell value is accepted', W, 'return LcpAt((long)v) == lk ? (long)v : -1;', 'return LcpAt((long)v) != lk ? (long)v : -1;', F_WIDE)
m('K42', 'the 8-bit escape is 254', W, 'v == byte.MaxValue : TLcp.ElementBytes == 2', 'v == byte.MaxValue - 1 : TLcp.ElementBytes == 2', F_WIDE)
m('K43', 'the 16-bit escape is 65,534', W, 'TLcp.ElementBytes == 2 && v == ushort.MaxValue)', 'TLcp.ElementBytes == 2 && v == ushort.MaxValue - 1)', F_WIDE)
m('K44', 'the escape is looked up one rank on', W, 'return escapes[rank];', 'return escapes[rank + 1];', F_WIDE)
m('K45', 'the 32-bit value is not checked against the text', W, 'if (TAccess.IsChecked && TLcp.ElementBytes >= 4 && v > (ulong)(n + 1))', 'if (false && TLcp.ElementBytes >= 4 && v > (ulong)(n + 1))', F_WIDE)
m('K46', 'the 32-bit value check starts at 8 bytes', W, 'TLcp.ElementBytes >= 4 && v > (ulong)(n + 1)', 'TLcp.ElementBytes >= 8 && v > (ulong)(n + 1)', F_WIDE)
m('K47', 'the 32-bit value check refuses a value equal to the length of the array', W, 'TLcp.ElementBytes >= 4 && v > (ulong)(n + 1)', 'TLcp.ElementBytes >= 4 && v > (ulong)n', F_WIDE)
m('K48', 'the LCP value is sign-extended from 32 bits', W, 'return (long)v;\n    }\n}\n\n/// <summary>', 'return (long)(int)v;\n    }\n}\n\n/// <summary>', F_WIDE)
m('K49', 'Depth reads the LCP array one rank on', W, 'public long Depth(long first) => LcpAt(first);', 'public long Depth(long first) => LcpAt(first + 1);', F_WIDE)
# ---- the pattern and the text
m('K50', 'the pattern is readable at its length', W, 'if ((ulong)index >= (ulong)_symbols.Length) ThrowIndex(index);', 'if ((ulong)index > (ulong)_symbols.Length) ThrowIndex(index);', F_WIDE)
m('K51', 'the pattern is not bounded', W, 'if ((ulong)index >= (ulong)_symbols.Length) ThrowIndex(index);', 'if (false) ThrowIndex(index);', F_WIDE)
m('K52', 'the packed text is read one symbol on', W, 'public int Symbol(long position) => MappedSearch.Symbol(text, bits, position);', 'public int Symbol(long position) => MappedSearch.Symbol(text, bits, position + 1);', F_WIDE)
m('K53', 'the pattern says its length minus one', W, 'public long Length => _symbols.Length;', 'public long Length => _symbols.Length - 1;', F_WIDE)
# ---- the switch of the matrix
m('K54', 'the stored width 5 takes the 4-byte readers', W, 'return s.PosWidth == 5\n            ? SelectLcp<TAccess, Pos40<TAccess>>', 'return s.PosWidth != 5\n            ? SelectLcp<TAccess, Pos40<TAccess>>', F_WIDE)
m('K55', 'the 8-bit and the 16-bit LCP arrays are swapped', W, 'LcpKind.U8 => Make<TAccess, TPos, Pos8<TAccess>>(s, sa, child, new Pos8<TAccess>(s.LcpBase, count)),', 'LcpKind.U8 => Make<TAccess, TPos, Pos16<TAccess>>(s, sa, child, new Pos16<TAccess>(s.LcpBase, count)),', F_WIDE)
m('K56', 'the 32-bit and the 64-bit LCP arrays are swapped', W, 'LcpKind.U32 => Make<TAccess, TPos, Pos32<TAccess>>(s, sa, child, new Pos32<TAccess>(s.LcpBase, count)),', 'LcpKind.U32 => Make<TAccess, TPos, Pos64<TAccess>>(s, sa, child, new Pos64<TAccess>(s.LcpBase, count)),', F_WIDE)
m('K57', 'the readers have n elements', W, 'long count = s.N + 1;\n        return s.PosWidth == 5', 'long count = s.N;\n        return s.PosWidth == 5', F_WIDE)
m('K58', 'the LCP readers have n elements', W, 'long count = s.N + 1;\n        return s.LcpKind switch', 'long count = s.N;\n        return s.LcpKind switch', F_WIDE)
m('K59', 'Checked and Trusted are swapped', W, 'return s.IsChecked ? Select<CheckedAccess>(s) : Select<TrustedAccess>(s);', 'return s.IsChecked ? Select<TrustedAccess>(s) : Select<CheckedAccess>(s);', F_WIDE)
m('K60', 'a tagged table over 5-byte cells is accepted', W, 'if (s.Tagged && s.PosWidth == 5) throw', 'if (s.Tagged && s.PosWidth == 55) throw', F_WIDE)
m('K61', 'a width of 6 is accepted', W, 'if (s.PosWidth is not (4 or 5)) throw', 'if (s.PosWidth is not (4 or 5 or 6)) throw', F_WIDE)
m('K62', 'the tagged flag is dropped', W, 'new WideCells<TAccess, TPos, TLcp>(child, lcpBase, s.Lcp, s.Tagged, s.N)', 'new WideCells<TAccess, TPos, TLcp>(child, lcpBase, s.Lcp, false, s.N)', F_WIDE)
m('K63', 'the cells read the suffix array section', W, 'new WideCells<TAccess, TPos, TLcp>(child, lcpBase, s.Lcp, s.Tagged, s.N)', 'new WideCells<TAccess, TPos, TLcp>(sa, lcpBase, s.Lcp, s.Tagged, s.N)', F_WIDE)
m('K64', 'the text has one bit more per symbol', W, 'new PackedSymbols(s.Text, s.Bits), s.N);', 'new PackedSymbols(s.Text, s.Bits + 1), s.N);', F_WIDE)

MS = 'MappedSearch.cs'
m('M01', 'the narrow test of the container is inverted', MS, 'if (!file.IsNarrow)', 'if (file.IsNarrow)', F_WIDE + '|' + F_INDEX.strip('()'))
m('M02', 'the long descent is not selected by its variant', MS, 'private const int WideVariant = 12;', 'private const int WideVariant = 11;', F_WIDE + '|' + F_INDEX.strip('()'))
m('M03', 'the last 32-bit variant falls through to the long descent', MS, 'case 11: {', 'case 12: {', F_INDEX)
m('M04', 'the long overload serves the narrow container with the long descent', MS, 'if (_wide is not null) return _wide.TryMatch(pattern, out lo, out hi);', 'if (_wide is null) return _wide!.TryMatch(pattern, out lo, out hi);', F_INDEX)
m('M05', 'the widened interval ends at its left bound', MS, 'bool found = TryMatch(pattern, out int l, out int h);\n        lo = l;\n        hi = h;', 'bool found = TryMatch(pattern, out int l, out int h);\n        lo = l;\n        hi = l;', F_INDEX)
m('M06', 'the 32-bit shape of an absent pattern is not the empty interval', MS, 'if (!found)\n        {\n            lo = 0;\n            hi = -1;\n            return false;\n        }', 'if (!found)\n        {\n            lo = 0;\n            hi = 0;\n            return false;\n        }', F_INDEX)
m('M07', 'the 32-bit shape answers the left bound for the right', MS, 'hi = IntBridge.Length(h, "EsaIndex.TryMatch, CountOccurrences or FindAllOccurrences', 'hi = IntBridge.Length(l, "EsaIndex.TryMatch, CountOccurrences or FindAllOccurrences', F_INDEX)
m('M08', 'a rank above int.MaxValue wraps in the 32-bit shape', MS, 'lo = IntBridge.Length(l, "EsaIndex.TryMatch, CountOccurrences or FindAllOccurrences (the 64-bit members are TryMatch with long ranks, LongCount and Locate)");', 'lo = unchecked((int)l);', F_INDEX, stage3='(' + T + 'Large.WideEngineSparseLargeTests)')
m('M09', 'the constructor of the 32-bit loop narrows the text length with an offset', MS, '_n = (int)file.TextReader.Length; // IsNarrow', '_n = (int)file.TextReader.Length - 1; // IsNarrow', F_INDEX)
m('M10', 'a long descent over views of a test is not marked wide', MS, 'internal bool IsWide => _wide is not null;', 'internal bool IsWide => _wide is null;', F_WIDE + '|' + F_INDEX.strip('()'))

E = 'EsaIndex.cs'
m('E01', 'the in-memory length is one short', E, '_length = text.Length;\n        _mapped = false;', '_length = text.Length - 1;\n        _mapped = false;', F_INDEX)
m('E02', 'an index over views is not mapped', E, '_length = text.Length;\n        _mapped = true;', '_length = text.Length;\n        _mapped = false;', F_INDEX)
m('E03', 'the 32-bit length is an unchecked narrowing', E, 'internal int Length => IntBridge.Length(_length, "EsaIndex.Length (the 64-bit member is TextLength)");', 'internal int Length => unchecked((int)_length);', F_INDEX, stage3='(' + T + 'Large.WideEngineSparseLargeTests)')
m('E04', 'the 64-bit length is the 32-bit one', E, 'internal long TextLength => _length;', 'internal long TextLength => (int)_length;', F_INDEX, stage3='(' + T + 'Large.WideEngineSparseLargeTests)')

SR = 'EsaIndex.Search.cs'
m('S01', 'Contains: the empty pattern is no substring', SR, 'internal bool Contains(ReadOnlySpan<char> pattern) => pattern.IsEmpty || TryMatch(pattern, out long _, out long _);', 'internal bool Contains(ReadOnlySpan<char> pattern) => pattern.IsEmpty && TryMatch(pattern, out long _, out long _);', F_INDEX)
m('S02', 'CountOccurrences: the empty pattern counts nothing', SR, 'if (pattern.IsEmpty) return Length;\n        return TryMatch(pattern, out int lo, out int hi) ? hi - lo + 1 : 0;', 'if (pattern.IsEmpty) return 0;\n        return TryMatch(pattern, out int lo, out int hi) ? hi - lo + 1 : 0;', F_INDEX)
m('S03', 'CountOccurrences: the count of an interval is hi - lo', SR, 'return TryMatch(pattern, out int lo, out int hi) ? hi - lo + 1 : 0;', 'return TryMatch(pattern, out int lo, out int hi) ? hi - lo : 0;', F_INDEX)
m('S04', 'LongCount: the empty pattern counts nothing', SR, 'if (pattern.IsEmpty) return _length;', 'if (pattern.IsEmpty) return 0;', F_INDEX)
m('S05', 'LongCount: the count of an interval is hi - lo', SR, 'return TryMatch(pattern, out long lo, out long hi) ? hi - lo + 1 : 0;', 'return TryMatch(pattern, out long lo, out long hi) ? hi - lo : 0;', F_INDEX)
m('S06', 'LongCount: an absent pattern occurs once', SR, 'return TryMatch(pattern, out long lo, out long hi) ? hi - lo + 1 : 0;', 'return TryMatch(pattern, out long lo, out long hi) ? hi - lo + 1 : 1;', F_INDEX)
m('S07', 'LongCount: the count is narrowed to 32 bits', SR, 'return TryMatch(pattern, out long lo, out long hi) ? hi - lo + 1 : 0;', 'return TryMatch(pattern, out long lo, out long hi) ? (int)(hi - lo + 1) : 0;', F_INDEX, stage3='(' + T + 'Large.WideEngineSparseLargeTests)')
m('S08', 'FindAllOccurrences: the list is one position short', SR, 'return Sa.ToArray(lo, hi - lo + 1);', 'return Sa.ToArray(lo, hi - lo);', F_INDEX)
m('S09', 'PageOf: a negative skip is accepted', SR, 'ArgumentOutOfRangeException.ThrowIfNegative(skip);\n        if ((ulong)lo > (ulong)hi || hi > length)', 'if ((ulong)lo > (ulong)hi || hi > length)', F_INDEX)
m('S10', 'PageOf: an interval of one rank is refused', SR, 'if ((ulong)lo > (ulong)hi || hi > length) throw', 'if ((ulong)lo >= (ulong)hi || hi > length) throw', F_INDEX)
m('S11', 'PageOf: the last rank is refused', SR, 'if ((ulong)lo > (ulong)hi || hi > length) throw', 'if ((ulong)lo > (ulong)hi || hi >= length) throw', F_INDEX)
m('S12', 'PageOf: a negative left bound is accepted', SR, 'if ((ulong)lo > (ulong)hi || hi > length) throw', 'if (lo > hi || hi > length) throw', F_INDEX)
m('S13', 'PageOf: the available count forgets the first rank', SR, 'long available = hi - lo + 1 - skip;', 'long available = hi - lo - skip;', F_INDEX)
m('S14', 'PageOf: the page starts one rank late', SR, 'start = lo + skip;\n        return available', 'start = lo + skip + 1;\n        return available', F_INDEX)
m('S15', 'PageOf: the page starts at the skip alone', SR, 'start = lo + skip;\n        return available', 'start = skip;\n        return available', F_INDEX)
m('S16', 'PageOf: the page is as long as what is available', SR, 'return available < capacity ? (int)available : capacity; // (int): below the capacity, which is an int', 'return available > capacity ? (int)available : capacity; // (int): below the capacity, which is an int', F_INDEX)
m('S17', 'PageOf: the count is narrowed from the whole interval', SR, 'return available < capacity ? (int)available : capacity; // (int): below the capacity, which is an int', 'return (int)Math.Min(hi - lo + 1, capacity);', F_INDEX)
m('S18', 'PageOf: the skip is narrowed to 32 bits', SR, 'start = lo + skip;\n        return available', 'start = lo + (int)skip;\n        return available', F_INDEX, stage3='(' + T + 'Large.WideEngineLargeTests)')
m('S19', 'Locate (int positions): a page of one position is not copied', SR, 'int count = PageOf(_length, lo, hi, skip, positions.Length, out long start);\n        if (count > 0) Sa.CopyTo(start, positions[..count]);\n        return count;\n    }\n\n    /// <summary>The same for 64-bit positions', 'int count = PageOf(_length, lo, hi, skip, positions.Length, out long start);\n        if (count > 1) Sa.CopyTo(start, positions[..count]);\n        return count;\n    }\n\n    /// <summary>The same for 64-bit positions', F_INDEX)
m('S20', 'Locate (long positions): the page is copied from the interval start', SR, 'int count = PageOf(_length, lo, hi, skip, positions.Length, out long start);\n        if (count > 0) Sa.CopyTo(start, positions[..count]);\n        return count;\n    }\n\n    /// <summary>\n    /// Finds the suffix-array interval', 'int count = PageOf(_length, lo, hi, skip, positions.Length, out long start);\n        if (count > 0) Sa.CopyTo(lo, positions[..count]);\n        return count;\n    }\n\n    /// <summary>\n    /// Finds the suffix-array interval', F_INDEX)
m('S21', 'TryMatch (int): the empty text has a match', SR, 'if (_length == 0) return false;\n        if (IsMapped) return MatchMapped(pattern, out lo, out hi);\n        return Prepared', 'if (_length == 1) return false;\n        if (IsMapped) return MatchMapped(pattern, out lo, out hi);\n        return Prepared', F_INDEX)
m('S22', 'TryMatch (long): the empty text has a match', SR, 'if (_length == 0) return false;\n        if (IsMapped) return MatchMapped(pattern, out lo, out hi);\n        bool found', 'if (_length == 1) return false;\n        if (IsMapped) return MatchMapped(pattern, out lo, out hi);\n        bool found', F_INDEX)
m('S23', 'TryMatch (long): the mapped test is inverted', SR, 'if (IsMapped) return MatchMapped(pattern, out lo, out hi);\n        bool found', 'if (!IsMapped) return MatchMapped(pattern, out lo, out hi);\n        bool found', F_INDEX)
m('S24', 'TryMatch (long): the interval of the in-memory index ends at its left bound', SR, 'bool found = TryMatch(pattern, out int l, out int h);\n        lo = l;\n        hi = h;\n        return found;\n    }\n\n    private bool Match<TSym>', 'bool found = TryMatch(pattern, out int l, out int h);\n        lo = l;\n        hi = l;\n        return found;\n    }\n\n    private bool Match<TSym>', F_INDEX)
m('S25', 'Dense: a letter just above the alphabet maps to a symbol', SR, 'int s = c < map.Length ? map[c] : -1;\n            if (s < 0) return false;\n            dense[i] = (byte)s;', 'int s = c <= map.Length ? map[c] : -1;\n            if (s < 0) return false;\n            dense[i] = (byte)s;', F_INDEX)
m('S26', 'the long reference compares one symbol late at a leaf', SR, 'long start = sa[a] + c;\n                if ((ulong)start > (ulong)n || m - c > n - start || !EdgeEquals(text, start, pat[checked((int)c)..], edge))', 'long start = sa[a] + c;\n                if ((ulong)start > (ulong)n || m - c > n - start || !EdgeEquals(text, start + 1, pat[checked((int)c)..], edge))', F_INDEX)
m('S27', 'the long reference clamps the text at n inclusive', SR, 'int key = (ulong)p < (ulong)n ? (int)text.Symbol(p) : -1;\n                if (key == symbol) break;\n                if (key > symbol || boundary < 0)\n                {\n                    lo = 0;\n                    hi = -1;\n                    return false;\n                }\n                a = boundary;\n                boundary = child.NextLIndex(boundary, rb);\n                if (isChecked && boundary >= 0 && (boundary <= a || boundary > rb)) ThrowDamagedChildTable();\n            }\n\n            bool childIsLast = boundary < 0;\n            if (a == b)\n            {\n                long start', 'int key = (ulong)p <= (ulong)n ? (int)text.Symbol(p) : -1;\n                if (key == symbol) break;\n                if (key > symbol || boundary < 0)\n                {\n                    lo = 0;\n                    hi = -1;\n                    return false;\n                }\n                a = boundary;\n                boundary = child.NextLIndex(boundary, rb);\n                if (isChecked && boundary >= 0 && (boundary <= a || boundary > rb)) ThrowDamagedChildTable();\n            }\n\n            bool childIsLast = boundary < 0;\n            if (a == b)\n            {\n                long start', F_INDEX)
m('S28', 'the long reference takes the larger limit', SR, 'long limit = Math.Min(childEll, m);\n            long from = sa[a] + c;', 'long limit = Math.Max(childEll, m);\n            long from = sa[a] + c;', F_INDEX)

AN = 'EsaIndex.Analysis.cs'
m('A01', 'LrsPosition: no repeat gives position 0', AN, 'internal long LrsPosition => Statistics.LrsDepth == 0 ? -1 : Sa[Statistics.LrsRank];', 'internal long LrsPosition => Statistics.LrsDepth == 0 ? 0 : Sa[Statistics.LrsRank];', F_ANALYSIS)
m('A02', 'LrsPosition: the rank of the position is one on', AN, 'internal long LrsPosition => Statistics.LrsDepth == 0 ? -1 : Sa[Statistics.LrsRank];', 'internal long LrsPosition => Statistics.LrsDepth == 0 ? -1 : Sa[Statistics.LrsRank + 1];', F_ANALYSIS)
m('A03', 'LrsLength is the statistics depth minus one', AN, 'internal long LrsLength => Statistics.LrsDepth;', 'internal long LrsLength => Statistics.LrsDepth - 1;', F_ANALYSIS)
m('A04', 'the ties: a run opens at its own rank', AN, 'if (runStart < 0) runStart = rank - 1;', 'if (runStart < 0) runStart = rank;', F_ANALYSIS)
m('A05', 'the ties: a run closes at its own rank', AN, 'results.Add(Tie(runStart, rank - 1, length));', 'results.Add(Tie(runStart, rank, length));', F_ANALYSIS)
m('A06', 'the ties: the run at the end of the array stops one rank early', AN, 'if (runStart >= 0) results.Add(Tie(runStart, _length, length));', 'if (runStart >= 0) results.Add(Tie(runStart, _length - 1, length));', F_ANALYSIS)
m('A07', 'the ties: the scan starts at rank 2', AN, 'long runStart = -1;      // left boundary of the interval whose run of the maximum value is open\n        for (long rank = 1;', 'long runStart = -1;      // left boundary of the interval whose run of the maximum value is open\n        for (long rank = 2;', F_ANALYSIS)
m('A08', 'the ties: the positions of a tie lose their last rank', AN, 'int[] positions = Sa.ToArray(lb, IntBridge.Length(rb - lb + 1, "EsaIndex.FindAllLongestRepeatedSubstrings (the positions of one longest repeat)"));', 'int[] positions = Sa.ToArray(lb, IntBridge.Length(rb - lb, "EsaIndex.FindAllLongestRepeatedSubstrings (the positions of one longest repeat)"));', F_ANALYSIS)
m('A09', 'the ties: the repeat is read from the second position', AN, 'return (TextView.Substring(positions[0], IntBridge.Length(length,', 'return (TextView.Substring(positions[^1], IntBridge.Length(length,', F_ANALYSIS)
m('A10', 'the distinct substrings: n (n + 1) / 2 is n n / 2', AN, 'Int128 total = length * (length + 1) / 2;', 'Int128 total = length * length / 2;', F_ANALYSIS)
m('A11', 'the distinct substrings: the LCP values are added', AN, 'for (long rank = 1; rank <= n; rank++) total -= cursor.Get(rank);', 'for (long rank = 1; rank <= n; rank++) total += cursor.Get(rank);', F_ANALYSIS)
m('A12', 'the distinct substrings: the last rank is not read', AN, 'for (long rank = 1; rank <= n; rank++) total -= cursor.Get(rank);', 'for (long rank = 1; rank < n; rank++) total -= cursor.Get(rank);', F_ANALYSIS)
m('A13', 'the distinct substrings: the sum is taken in 64 bits (a^n above 3 x 10^9 wraps)', AN, 'Int128 length = n;\n        Int128 total = length * (length + 1) / 2;', 'long length = n;\n        Int128 total = length * (length + 1) / 2;', F_ANALYSIS, stage3='(' + T + 'Large.WideEngineLargeTests)')
m('A14', 'the distinct substrings: the 64-bit result is not checked', AN, 'if (total > long.MaxValue) throw new NotSupportedException($"the number of distinct substrings is {total}, more than 64 bits hold");\n        return (long)total;', 'return (long)total;', F_ANALYSIS)
m('A15', 'the histogram: the suffix count n - i + 1 loses its one', AN, 'counts[i] = n - i + 1 - atLeast;', 'counts[i] = n - i - atLeast;', F_ANALYSIS)
m('A16', 'the histogram: the length 1 is not counted', AN, 'for (int i = m; i >= 1; i--)', 'for (int i = m; i > 1; i--)', F_ANALYSIS)
m('A17', 'the histogram: the empty substring is not counted', AN, 'counts[0] = 1;', 'counts[0] = 0;', F_ANALYSIS)
m('A18', 'the histogram: a value is clamped one too late', AN, 'if (value > 0) counts[value < m ? (int)value : m]++;', 'if (value > 0) counts[value < m ? (int)value : m - 1]++;', F_ANALYSIS)
m('A19', 'the suffixes: the list starts at rank 0', AN, 'for (int rank = 1; rank <= suffixes.Length; rank++) suffixes[rank - 1] = TextView.Substring(sa[(long)rank]);', 'for (int rank = 0; rank < suffixes.Length; rank++) suffixes[rank] = TextView.Substring(sa[(long)rank]);', F_ANALYSIS)
m('A20', 'the suffixes: the lazy walk stops one suffix short', AN, 'for (long rank = 1; rank <= _length; rank++) yield return TextView.Substring(sa[rank]);', 'for (long rank = 1; rank < _length; rank++) yield return TextView.Substring(sa[rank]);', F_ANALYSIS)
m('A21', 'the repeat is read at the next rank', AN, 'TextView.Substring(Sa[Statistics.LrsRank], IntBridge.Length(Statistics.LrsDepth, "EsaIndex.LongestRepeatedSubstring', 'TextView.Substring(Sa[Statistics.LrsRank + 1], IntBridge.Length(Statistics.LrsDepth, "EsaIndex.LongestRepeatedSubstring', F_ANALYSIS)

TR = 'EsaIndex.Traverse.cs'
m('T01', 'the root event counts n leaves', TR, 'visitor.VisitNode(0, 0, n + 1, CountChildren(firstOfRoot, n), 0);', 'visitor.VisitNode(0, 0, n, CountChildren(firstOfRoot, n), 0);', F_TRAVERSE)
m('T02', 'the event of the empty text has a leaf', TR, 'visitor.VisitNode(0, 0, 0, 0, 0);\n            return;', 'visitor.VisitNode(0, 0, 1, 0, 0);\n            return;', F_TRAVERSE)
m('T03', 'the end of the last child is the right bound of its parent minus one', TR, 'long b = lastChild ? f.Rb : f.NextBoundary - 1;', 'long b = lastChild ? f.Rb - 1 : f.NextBoundary - 1;', F_TRAVERSE)
m('T04', 'the end of a child is its next boundary', TR, 'long b = lastChild ? f.Rb : f.NextBoundary - 1;', 'long b = lastChild ? f.Rb : f.NextBoundary;', F_TRAVERSE)
m('T05', 'the next child starts at the end of this one', TR, 'frames[top] = new TraverseFrame { A = b + 1,', 'frames[top] = new TraverseFrame { A = b,', F_TRAVERSE)
m('T06', 'the branch key is read one symbol on', TR, 'long p = sa[a] + f.Ell;\n            visitor.EnterBranch((ulong)p < (ulong)n ? text[p] : -1);', 'long p = sa[a] + f.Ell + 1;\n            visitor.EnterBranch((ulong)p < (ulong)n ? text[p] : -1);', F_TRAVERSE)
m('T07', 'the key at the end of the text is a symbol', TR, 'visitor.EnterBranch((ulong)p < (ulong)n ? text[p] : -1);', 'visitor.EnterBranch((ulong)p <= (ulong)n ? text[p] : -1);', F_TRAVERSE)
m('T08', 'the depth of a leaf is its own', TR, 'visitor.VisitNode(p, -1, 1, 0, f.Ell);', 'visitor.VisitNode(p, -1, 1, 0, f.Ell + 1);', F_TRAVERSE)
m('T09', 'a leaf has two leaves', TR, 'visitor.VisitNode(p, -1, 1, 0, f.Ell);', 'visitor.VisitNode(p, -1, 2, 0, f.Ell);', F_TRAVERSE)
m('T10', 'the start of an internal node ignores the parent depth', TR, 'visitor.VisitNode(occ + f.Ell, occ + ell, b - a + 1,', 'visitor.VisitNode(occ, occ + ell, b - a + 1,', F_TRAVERSE)
m('T11', 'the end of an internal node is at the parent depth', TR, 'visitor.VisitNode(occ + f.Ell, occ + ell, b - a + 1,', 'visitor.VisitNode(occ + f.Ell, occ + f.Ell, b - a + 1,', F_TRAVERSE)
m('T12', 'the leaf count of an internal node is hi - lo', TR, 'visitor.VisitNode(occ + f.Ell, occ + ell, b - a + 1,', 'visitor.VisitNode(occ + f.Ell, occ + ell, b - a,', F_TRAVERSE)
m('T13', 'the children of a node are counted from zero', TR, 'long count = 1;\n        for (long k = first; k >= 0;)', 'long count = 0;\n        for (long k = first; k >= 0;)', F_TRAVERSE)
m('T14', 'the first l-index equal to the left bound passes', TR, 'if (first <= lb || first > rb) throw DamagedChildTable();', 'if (first < lb || first > rb) throw DamagedChildTable();', F_TRAVERSE)
m('T15', 'a position equal to n is outside the text', TR, 'if ((ulong)position > (ulong)n) throw new InvalidDataException("the suffix array holds a position outside the text: the file is damaged");', 'if ((ulong)position >= (ulong)n) throw new InvalidDataException("the suffix array holds a position outside the text: the file is damaged");', F_TRAVERSE)
m('T16', 'the bit of a position is a 32-bit shift', TR, 'ulong bit = 1UL << (int)(position & 63);', 'ulong bit = 1UL << (int)(position & 31);', F_TRAVERSE)
m('T17', 'the bitmap has one word too few', TR, 'ulong[]? seen = mapped ? new ulong[(n >> 6) + 1] : null;', 'ulong[]? seen = mapped ? new ulong[n >> 6] : null;', F_TRAVERSE)
m('T18', 'the bitmap word of a position is the position', TR, 'ref ulong word = ref seen[position >> 6];', 'ref ulong word = ref seen[position >> 7];', F_TRAVERSE)
m('T19', 'the 32-bit visitor: a leaf end of 0 is boundless', TR, 'endIndex < 0 ? -1 : IntBridge.Length(endIndex, Member)', 'endIndex <= 0 ? -1 : IntBridge.Length(endIndex, Member)', F_TRAVERSE)
m('T20', 'the 32-bit visitor: the depth is narrowed without a check', TR, 'IntBridge.Length(leafCount, Member), IntBridge.Length(childCount, Member), IntBridge.Length(depth, Member));', 'IntBridge.Length(leafCount, Member), IntBridge.Length(childCount, Member), unchecked((int)depth));', F_TRAVERSE)
m('T21', 'the 32-bit visitor: the leaf count is narrowed without a check', TR, 'IntBridge.Length(startIndex, Member), endIndex < 0 ? -1 : IntBridge.Length(endIndex, Member), IntBridge.Length(leafCount, Member),', 'IntBridge.Length(startIndex, Member), endIndex < 0 ? -1 : IntBridge.Length(endIndex, Member), unchecked((int)leafCount),', F_TRAVERSE)
m('T22', 'the 32-bit visitor: the start is narrowed without a check', TR, 'inner.VisitNode(IntBridge.Length(startIndex, Member),', 'inner.VisitNode(unchecked((int)startIndex),', F_TRAVERSE)

LO = 'LeftmostOccurrence.cs'
m('L01', 'the table is refused one rank late', LO, 'if (n + 1 > Array.MaxLength) throw', 'if (n > Array.MaxLength) throw', F_SCAN)
m('L02', 'the LCP array of another length is accepted', LO, 'if (lcp.Length != n + 1) throw', 'if (lcp.Length != n) throw', F_SCAN)
m('L03', 'the occurrence of an interval is one too far', LO, 'occurrence[first[sp]] = IntBridge.Value(min[sp],', 'occurrence[first[sp]] = IntBridge.Value(min[sp] + 1,', F_SCAN)
m('L04', 'the occurrence is not carried to the parent', LO, 'carry = min[sp];\n                if (sp > 0 && carry < min[sp - 1]) min[sp - 1] = carry;', 'carry = leaf;\n                if (sp > 0 && carry < min[sp - 1]) min[sp - 1] = carry;', F_SCAN)
m('L05', 'the parent keeps the larger occurrence', LO, 'if (sp > 0 && carry < min[sp - 1]) min[sp - 1] = carry;', 'if (sp > 0 && carry > min[sp - 1]) min[sp - 1] = carry;', F_SCAN)
m('L06', 'an interval is closed when its value equals the next', LO, 'while (sp > 0 && ell[sp - 1] > cur)', 'while (sp > 0 && ell[sp - 1] >= cur)', F_SCAN)
m('L07', 'an interval is pushed when its value equals the top', LO, 'if (ell[top] < cur)\n            {\n                if (sp == capacity)', 'if (ell[top] <= cur)\n            {\n                if (sp == capacity)', F_SCAN)
m('L08', 'the first l-index of an interval is one rank on', LO, 'ell[sp] = cur;\n                first[sp] = i;', 'ell[sp] = cur;\n                first[sp] = i + 1;', F_SCAN)
m('L09', 'the last value is not read', LO, 'long cur = i <= n ? stream.Get(i) : -1;', 'long cur = i < n ? stream.Get(i) : -1;', F_SCAN)
m('L10', 'the leaf of the previous rank is the leaf of this one', LO, 'long leaf = sa[i - 1];', 'long leaf = sa[i];', F_SCAN)
m('L11', 'a damaged LCP array without a first l-index is accepted', LO, 'if (first[sp] < 0) throw DamagedLcp();', 'if (first[sp] < -1) throw DamagedLcp();', F_SCAN)
m('L12', 'the occurrence is narrowed without a check', LO, 'occurrence[first[sp]] = IntBridge.Value(min[sp], "a leftmost occurrence (a position of a text below 2^31 symbols)");', 'occurrence[first[sp]] = (int)min[sp];', F_SCAN)

SC = 'LcpIntervalScan.cs'
m('W01', 'the wide scan: a zero value is invalid', SC, 'if (i <= n && (cur < 0 || (i == 1 && cur != 0)))\n                throw new ArgumentException($"LCP[{i}] = {cur} is not a valid LCP value (LCP[1] must be 0 and no value may be negative)");\n\n            long left', 'if (i <= n && (cur <= 0 || (i == 1 && cur != 0)))\n                throw new ArgumentException($"LCP[{i}] = {cur} is not a valid LCP value (LCP[1] must be 0 and no value may be negative)");\n\n            long left', F_SCAN)
m('W02', 'the wide scan: an interval starts at its own rank', SC, 'long left = i - 1;', 'long left = i;', F_SCAN)
m('W03', 'the wide scan: an interval closes at its own rank', SC, 'long b = lb[sp];\n                sink.Close(ell[sp], b, i - 1, first[sp], sp == 0 || cur < ell[sp - 1]);', 'long b = lb[sp];\n                sink.Close(ell[sp], b, i, first[sp], sp == 0 || cur < ell[sp - 1]);', F_SCAN)
m('W04', 'the wide scan: the last-child flag of an equal parent', SC, 'long b = lb[sp];\n                sink.Close(ell[sp], b, i - 1, first[sp], sp == 0 || cur < ell[sp - 1]);', 'long b = lb[sp];\n                sink.Close(ell[sp], b, i - 1, first[sp], sp == 0 || cur <= ell[sp - 1]);', F_SCAN)
m('W05', 'the wide scan: an interval is pushed when its value equals the top', SC, 'int top = sp - 1;\n            if (ell[top] < cur)\n            {\n                if (sp == capacity)\n                {\n                    if (capacity > Array.MaxLength / 2)', 'int top = sp - 1;\n            if (ell[top] <= cur)\n            {\n                if (sp == capacity)\n                {\n                    if (capacity > Array.MaxLength / 2)', F_SCAN)
m('W06', 'the wide scan: the l-index link is to the first', SC, 'sink.LIndex(last[top], i);\n                last[top] = i;\n            }\n        }\n        return peak;\n    }\n\n    /// <summary>The header statistics', 'sink.LIndex(first[top], i);\n                last[top] = i;\n            }\n        }\n        return peak;\n    }\n\n    /// <summary>The header statistics', F_SCAN)
m('W07', 'the wide scan: the last l-index is not recorded', SC, 'sink.LIndex(last[top], i);\n                last[top] = i;\n            }\n        }\n        return peak;\n    }\n\n    /// <summary>The header statistics', 'sink.LIndex(last[top], i);\n            }\n        }\n        return peak;\n    }\n\n    /// <summary>The header statistics', F_SCAN)
m('W08', 'the wide scan: the peak starts at zero', SC, 'long peak = 1;\n        for (long i = 1; i <= n + 1; i++)', 'long peak = 0;\n        for (long i = 1; i <= n + 1; i++)', F_SCAN)
m('W09', 'the wide scan: the empty text has a peak', SC, 'if (n == 0) return 0;\n\n        int capacity = InitialCapacity;\n        long[]', 'if (n == 0) return 1;\n\n        int capacity = InitialCapacity;\n        long[]', F_SCAN)
m('W10', 'the wide statistics: the node count lacks the leaves', SC, ': new EsaStatistics(n, 1 + sink.Internal + (n + 1), sink.Internal, n, n, sink.LrsDepth, sink.LrsRank, peak);', ': new EsaStatistics(n, 1 + sink.Internal + n, sink.Internal, n, n, sink.LrsDepth, sink.LrsRank, peak);', F_SCAN)
m('W11', 'the wide statistics: the root counts as an internal node', SC, 'public void Close(long ell, long lb, long rb, long first, bool lastChild)\n    {\n        if (ell <= 0) return;\n        Internal++;', 'public void Close(long ell, long lb, long rb, long first, bool lastChild)\n    {\n        if (ell < 0) return;\n        Internal++;', F_SCAN)
m('W12', 'the wide statistics: the last of equal depths wins', SC, 'public void Close(long ell, long lb, long rb, long first, bool lastChild)\n    {\n        if (ell <= 0) return;\n        Internal++;\n        if (ell > LrsDepth)', 'public void Close(long ell, long lb, long rb, long first, bool lastChild)\n    {\n        if (ell <= 0) return;\n        Internal++;\n        if (ell >= LrsDepth)', F_SCAN)
m('W13', 'the wide statistics: the repeat starts at the right bound', SC, 'LrsDepth = ell;\n            LrsRank = lb;\n        }\n    }\n}\n\n/// <summary>The scan at 64 bits', 'LrsDepth = ell;\n            LrsRank = rb;\n        }\n    }\n}\n\n/// <summary>The scan at 64 bits', F_SCAN)
m('W14', 'the wide stream starts at rank 0', SC, 'private long _next;\n\n    internal WideCompactLcpStream(CompactLcp lcp)\n    {\n        _cursor = lcp.CreateLongCursor();\n        _next = 1;', 'private long _next;\n\n    internal WideCompactLcpStream(CompactLcp lcp)\n    {\n        _cursor = lcp.CreateLongCursor();\n        _next = 0;', F_SCAN)

CL = 'CompactLcp.cs'
m('C01', 'the long cursor reads a wide array through the 32-bit cursor', CL, '_wide = owner.IsWide || owner.Length > int.MaxValue;', '_wide = owner.Length > int.MaxValue;', F_LCP, stage3='(' + T + 'Large.WideEngineLargeTests)')
m('C02', 'the long cursor never reads through the store', CL, '_wide = owner.IsWide || owner.Length > int.MaxValue;', '_wide = false;', F_LCP, stage3='(' + T + 'Large.WideEngineLargeTests)')
m('C03', 'the long cursor reads the rank one on', CL, 'return _narrow.Get((int)rank); // (int): inside the array', 'return _narrow.Get((int)rank + 1); // (int): inside the array', F_LCP)
m('C04', 'the long cursor does not check its rank', CL, 'if ((ulong)rank >= (ulong)_owner.Length) _owner.ThrowRankOutside(rank);\n            return _narrow', 'return _narrow', F_LCP)
m('C05', 'the long cursor narrows a long rank of a wide array', CL, 'if (_wide) return _owner[rank];', 'if (_wide) return _owner[(int)rank];', F_LCP, stage3='(' + T + 'Large.WideEngineLargeTests)')

SI = 'EsaSuffixIndex.cs'
m('X01', 'Locate: the empty suffix is a position', SI, 'long from = Math.Max(locus.Lo, 1);', 'long from = locus.Lo;', F_CONTRACT)
m('X02', 'Locate: a locus of one rank has no position', SI, 'if (from > locus.Hi) return 0;\n        return _index.Locate(', 'if (from >= locus.Hi) return 0;\n        return _index.Locate(', F_CONTRACT)
m('X03', 'Locate: the interval ends one rank early', SI, 'return _index.Locate(from, locus.Hi, skip, positions);', 'return _index.Locate(from, locus.Hi - 1, skip, positions);', F_CONTRACT)
m('X04', 'Locate: the page is narrowed to 32-bit ranks', SI, 'return _index.Locate(from, locus.Hi, skip, positions);', 'return _index.Locate((int)from, locus.Hi, skip, positions);', F_CONTRACT, stage3='(' + T + 'Large.WideEngineSparseLargeTests)')
m('X05', 'ExtendRight: the symbol inside an edge is read one position on', SI, 'long at = node.Pos + length;', 'long at = node.Pos + length + 1;', F_CONTRACT)
m('X06', 'ExtendRight: the end of the text is a symbol', SI, 'if ((ulong)at >= (ulong)_n || text[at] != symbol) return false;', 'if ((ulong)at > (ulong)_n || text[at] != symbol) return false;', F_CONTRACT)
m('X07', 'ExtendRight: the symbol of an edge is not compared', SI, 'if ((ulong)at >= (ulong)_n || text[at] != symbol) return false;', 'if ((ulong)at >= (ulong)_n) return false;', F_CONTRACT)
m('X08', 'ExtendRight: the position inside an edge is narrowed', SI, 'long at = node.Pos + length;', 'long at = (int)node.Pos + length;', F_CONTRACT, stage3='(' + T + 'Large.WideEngineSparseLargeTests)')
m('X09', 'ExtendRight: the locus is one symbol short', SI, 'result = new IndexLocus(node.Lb, node.Rb, length, node.Depth);', 'result = new IndexLocus(node.Lb, node.Rb, length - 1, node.Depth);', F_CONTRACT)
m('X10', 'the locus check accepts a hi above n', SI, 'bool valid = locus.Lo >= 0 && locus.Lo <= locus.Hi && locus.Hi <= n && locus.Length >= 0', 'bool valid = locus.Lo >= 0 && locus.Lo <= locus.Hi && locus.Hi <= n + 1 && locus.Length >= 0', F_CONTRACT)
m('X11', 'Resolve: the root has the depth of the text', SI, 'if (depth == 0) return new EsaNode64(0, _n, 0, 0, _n, true);', 'if (depth == 0) return new EsaNode64(0, _n, 0, 0, 0, true);', F_CONTRACT)
m('X12', 'Resolve: the right neighbour of the last rank is read', SI, 'long right = rb < _n ? lcp[rb + 1] : -1;', 'long right = rb <= _n ? lcp[rb + 1] : -1;', F_CONTRACT)
m('X13', 'Resolve: the last-child test is not strict', SI, 'bool lastChild = left > right;', 'bool lastChild = left >= right;', F_CONTRACT)
m('X14', 'Resolve: the parent depth is the left LCP', SI, 'long parentDepth = Math.Max(left, right);', 'long parentDepth = left;', F_CONTRACT)
m('X15', 'Resolve: the depth of a leaf lacks the terminator', SI, 'if (lb == rb) valid = depth == _n - pos + 1;', 'if (lb == rb) valid = depth == _n - pos;', F_CONTRACT)
m('X16', 'Resolve: the first l-index may be the left bound', SI, 'valid = first > lb && first <= rb && lcp[first] == depth;', 'valid = first >= lb && first <= rb && lcp[first] == depth;', F_CONTRACT)
m('X17', 'Resolve: the depth of an interval is not compared', SI, 'valid = first > lb && first <= rb && lcp[first] == depth;', 'valid = first > lb && first <= rb;', F_CONTRACT)
m('X18', 'Resolve: the matched length may equal the parent depth', SI, 'bool valid = parentDepth < locus.Length;', 'bool valid = parentDepth <= locus.Length;', F_CONTRACT)
m('X19', 'CheckRank: the rank n is outside', SI, 'if ((ulong)rank > (ulong)_n) throw new ArgumentOutOfRangeException(parameterName, rank,', 'if ((ulong)rank >= (ulong)_n) throw new ArgumentOutOfRangeException(parameterName, rank,', F_CONTRACT)
m('X20', 'TryExtendLeft: a symbol at position 0 has a predecessor', SI, 'if (position > 0 && text[position - 1] == symbol)', 'if (position >= 0 && text[position - 1] == symbol)', F_CONTRACT)
m('X21', 'TryExtendLeft: the single rank locus has the depth of the suffix', SI, 'long aux = low == high ? _n - sa[low] + 1 : pack.Pyramid.Min(low + 1, high);', 'long aux = low == high ? _n - sa[low] : pack.Pyramid.Min(low + 1, high);', F_CONTRACT)
m('X22', 'TryExtendLeft: the new locus is as long as the old', SI, 'result = new IndexLocus(low, high, at.Length + 1, aux);', 'result = new IndexLocus(low, high, at.Length, aux);', F_CONTRACT)
m('X23', 'TryExtendLeft: the range minimum starts at the left rank', SI, 'pack.Pyramid.Min(low + 1, high);', 'pack.Pyramid.Min(low, high);', F_CONTRACT)
m('X24', 'TryExtendRight: the symbol 65535 is no symbol', SI, 'public bool TryExtendRight(in IndexLocus locus, int symbol, out IndexLocus result)\n    {\n        if (symbol is < 0 or > char.MaxValue)', 'public bool TryExtendRight(in IndexLocus locus, int symbol, out IndexLocus result)\n    {\n        if (symbol is < 0 or >= char.MaxValue)', F_CONTRACT)
m('X25', 'Count: the empty suffix is a position', SI, 'return locus.Hi - Math.Max(locus.Lo, 1) + 1;', 'return locus.Hi - locus.Lo + 1;', F_CONTRACT)

TRE = 'EsaSuffixIndex.Tree.cs'
m('Y01', 'FindChild: the symbol is read at the parent depth minus one', TRE, 'long at = saA + node.Depth;\n            int symbol', 'long at = saA + node.Depth - 1;\n            int symbol', F_CONTRACT)
m('Y02', 'FindChild: a larger symbol continues the walk', TRE, 'if (symbol > key || boundary < 0)', 'if (symbol < key || boundary < 0)', F_CONTRACT)
m('Y03', 'FindChild: the child ends at its boundary', TRE, 'long b = boundary >= 0 ? boundary - 1 : rb;\n                child = MakeChild', 'long b = boundary >= 0 ? boundary : rb;\n                child = MakeChild', F_CONTRACT)
m('Y04', 'FindChild: the position is clamped at n inclusive', TRE, 'int symbol = (ulong)at < (ulong)_n ? text[at] : -1;', 'int symbol = (ulong)at <= (ulong)_n ? text[at] : -1;', F_CONTRACT)
m('Y05', 'MakeChild: the depth of a leaf lacks the terminator', TRE, 'return new EsaNode64(a, a, _n - saA + 1, parentDepth, saA, childIsLast);', 'return new EsaNode64(a, a, _n - saA, parentDepth, saA, childIsLast);', F_CONTRACT)
m('Y06', 'MakeChild: a first l-index at the left bound passes', TRE, 'if (first <= a || first > b) ThrowCorruptChildTable();', 'if (first < a || first > b) ThrowCorruptChildTable();', F_CONTRACT)
m('Y07', 'MakeChild: a child as deep as its parent passes', TRE, 'if (_isMapped && depth <= parentDepth) ThrowCorruptChildTable();', 'if (_isMapped && depth < parentDepth) ThrowCorruptChildTable();', F_CONTRACT)
m('Y08', 'MakeChild: a leaf no longer than its parent passes', TRE, 'if (_isMapped && _n - saA + 1 <= parentDepth) ThrowCorruptChildTable();', 'if (_isMapped && _n - saA + 1 < parentDepth) ThrowCorruptChildTable();', F_CONTRACT)
m('Y09', 'CheckBoundary: a boundary at the previous rank passes', TRE, 'if (boundary >= 0 && (boundary <= previous || boundary > rb)) ThrowCorruptChildTable();', 'if (boundary >= 0 && (boundary < previous || boundary > rb)) ThrowCorruptChildTable();', F_CONTRACT)
m('Y10', 'GetChildren: a destination as long as the count is too short', TRE, 'if (count > children.Length) return count;', 'if (count >= children.Length) return count;', F_CONTRACT)
m('Y11', 'GetChildren: the children are counted from one', TRE, 'int count = 0;\n        long boundary = table.FirstLIndex(node.Lb, rb, node.LastChild);', 'int count = 1;\n        long boundary = table.FirstLIndex(node.Lb, rb, node.LastChild);', F_CONTRACT)
m('Y12', 'EdgeSymbol: an offset equal to n is inside', TRE, 'if (offset >= _n) return -1;', 'if (offset > _n) return -1;', F_CONTRACT)
m('Y13', 'EdgeSymbol: the end of the text is a symbol', TRE, 'return index < _n ? _index.TextView[index] : -1;', 'return index <= _n ? _index.TextView[index] : -1;', F_CONTRACT)
m('Y14', 'EdgeSymbol: the edge starts at the node depth', TRE, 'long index = resolved.Pos + resolved.ParentDepth + offset;', 'long index = resolved.Pos + resolved.Depth + offset;', F_CONTRACT)
m('Y15', 'TryParent: the root is the parent of every node', TRE, 'if (parentDepth == 0)\n        {\n            parent = Root;', 'if (parentDepth <= 1)\n        {\n            parent = Root;', F_CONTRACT)
m('Y16', 'TryParent: the left border is read at the right bound', TRE, 'long left = pyramid.PrevLess(node.Lb, parentDepth);', 'long left = pyramid.PrevLess(node.Rb, parentDepth);', F_CONTRACT)
m('Y17', 'TryParent: the right border is one rank short', TRE, 'parent = new IndexLocus(left, right - 1, parentDepth, parentDepth);', 'parent = new IndexLocus(left, right, parentDepth, parentDepth);', F_CONTRACT)

AR = 'EsaSuffixIndex.Arrays.cs'
m('R01', 'Lcp: the minimum starts at the lower rank', AR, 'return _index.Navigation.Pyramid.Min(Math.Min(rankA, rankB) + 1, Math.Max(rankA, rankB));', 'return _index.Navigation.Pyramid.Min(Math.Min(rankA, rankB), Math.Max(rankA, rankB));', F_CONTRACT)
m('R02', 'Lcp: equal ranks have the depth of the suffix without its end', AR, 'if (rankA == rankB) return _n - _index.Sa[rankA];', 'if (rankA == rankB) return _n - _index.Sa[rankA] + 1;', F_CONTRACT)
m('R03', 'Lca: a leaf has the length of its suffix', AR, 'long depth = _n - _index.Sa[rankA] + 1;', 'long depth = _n - _index.Sa[rankA];', F_CONTRACT)
m('R04', 'Lca: the right border is not one rank back', AR, 'pyramid.NextLess(hi + 1, value) - 1, value, value);', 'pyramid.NextLess(hi + 1, value), value, value);', F_CONTRACT)
m('R05', 'Lca: the left border is read at the right rank', AR, 'return new IndexLocus(pyramid.PrevLess(lo, value),', 'return new IndexLocus(pyramid.PrevLess(hi, value),', F_CONTRACT)
m('R06', 'Lca: the minimum starts at the lower rank', AR, 'long value = pyramid.Min(lo + 1, hi);', 'long value = pyramid.Min(lo, hi);', F_CONTRACT)
m('R07', 'Lca: the empty text is not the root', AR, 'if (_n == 0) return Root;', 'if (_n == 1) return Root;', F_CONTRACT)
m('R08', 'SuffixAt: the rank is read one on', AR, 'return _index.Sa[rank];', 'return _index.Sa[rank + 1];', F_CONTRACT)
m('R09', 'LcpAt: the rank is read one on', AR, 'return _index.Lcp[rank];', 'return _index.Lcp[rank + 1];', F_CONTRACT)
m('R10', 'SuffixAt: the rank is narrowed to 32 bits', AR, 'return _index.Sa[rank];', 'return _index.Sa[(int)rank];', F_CONTRACT, stage3='(' + T + 'Large.WideEngineSparseLargeTests)')
m('R11', 'LcpAt: the rank is narrowed to 32 bits', AR, 'return _index.Lcp[rank];', 'return _index.Lcp[(int)rank];', F_CONTRACT, stage3='(' + T + 'Large.WideEngineSparseLargeTests)')

SF = 'External/Sesa/SesaFile.cs'
m('F01', 'the narrow test admits n = int.MaxValue', SF, 'bool narrow = n < int.MaxValue && posWidth == 4 && parsed.LcpKind != LcpKind.U64;', 'bool narrow = n <= int.MaxValue && posWidth == 4 && parsed.LcpKind != LcpKind.U64;', F_SESA, stage3='(' + T + 'Large.WideEngineSparseLargeTests)')
m('F02', 'the narrow test admits n = int.MaxValue - 2 only', SF, 'bool narrow = n < int.MaxValue && posWidth == 4 && parsed.LcpKind != LcpKind.U64;', 'bool narrow = n < int.MaxValue - 1 && posWidth == 4 && parsed.LcpKind != LcpKind.U64;', F_SESA, stage3='(' + T + 'Large.WideEngineSparseLargeTests)')
m('F03', 'the statistics: the internal nodes are NODE_COUNT - n', SF, 'h.NodeCount - 1 - (n + 1)', 'h.NodeCount - 1 - n', F_SESA, stage3='(' + T + 'Large.WideEngineSparseLargeTests)')
m('F04', 'the statistics: a text without a repeat has the LRS rank of the header', SF, 'h.LrsDepth == 0 ? 0 : h.LrsRank', 'h.LrsRank', F_SESA)
m('F05', 'the statistics: the length is the node count', SF, 'return new EsaStatistics(n, h.NodeCount,', 'return new EsaStatistics(h.NodeCount, h.NodeCount,', F_SESA, stage3='(' + T + 'Large.WideEngineSparseLargeTests)')
m('F06', 'the statistics: the length is narrowed', SF, 'long n = h.TextLength;\n            if (n == 0) return new EsaStatistics(0, 1, 0, 0, 0, 0, 0, 0);', 'long n = (int)h.TextLength;\n            if (n == 0) return new EsaStatistics(0, 1, 0, 0, 0, 0, 0, 0);', F_SESA, stage3='(' + T + 'Large.WideEngineSparseLargeTests)')
m('F07', 'the wide sections: the tagged flag is inverted', SF, '_parsed.Child.Codec == SesaFormat.ChildCodecTagged, _data + _parsed.Sa.Offset', '_parsed.Child.Codec != SesaFormat.ChildCodecTagged, _data + _parsed.Sa.Offset', F_SESA, stage3='(' + T + 'Large.WideEngineSparseLargeTests)')
m('F08', 'the wide sections: the suffix array and the child table are swapped', SF, '_data + _parsed.Sa.Offset, _data + _parsed.Child.Offset, _data + _parsed.Lcp.Offset', '_data + _parsed.Child.Offset, _data + _parsed.Sa.Offset, _data + _parsed.Lcp.Offset', F_SESA, stage3='(' + T + 'Large.WideEngineSparseLargeTests)')
m('F09', 'the wide sections: the text offset is the LCP offset', SF, '_data + _parsed.Text.Offset, _parsed.Bits, _parsed.Header.TextLength);', '_data + _parsed.Lcp.Offset, _parsed.Bits, _parsed.Header.TextLength);', F_SESA, stage3='(' + T + 'Large.WideEngineSparseLargeTests)')
m('F10', 'the wide sections: the stored width is 4', SF, 'new(IsChecked, _parsed.Header.PosWidth, _parsed.LcpKind,', 'new(IsChecked, 4, _parsed.LcpKind,', F_SESA, stage3='(' + T + 'Large.WideEngineSparseLargeTests)')
m('F11', 'the wide sections: the access mode is inverted', SF, 'new(IsChecked, _parsed.Header.PosWidth', 'new(!IsChecked, _parsed.Header.PosWidth', F_SESA, stage3='(' + T + 'Large.WideEngineSparseLargeTests)')

NO = 'EsaNode.cs'
m('N01', 'EsaNode64: a leaf test on the right bound only', NO, 'internal bool IsLeaf => Lb == Rb;\n\n    /// <summary>True for the null node.</summary>\n    internal bool IsNull => Lb < 0;\n\n    /// <summary>Length of the incoming edge, <c>Depth - ParentDepth</c> (0 for the root).</summary>\n    internal long EdgeLength => Depth - ParentDepth;\n\n    /// <summary>Number of ranks in the interval (the leaf count; the root counts the empty suffix too).</summary>\n    internal long RankCount => Rb - Lb + 1;\n\n    /// <summary>The node as the 32-bit navigator', 'internal bool IsLeaf => Lb == Rb;\n\n    /// <summary>True for the null node.</summary>\n    internal bool IsNull => Lb < 0;\n\n    /// <summary>Length of the incoming edge, <c>Depth - ParentDepth</c> (0 for the root).</summary>\n    internal long EdgeLength => Depth - ParentDepth;\n\n    /// <summary>Number of ranks in the interval (the leaf count; the root counts the empty suffix too).</summary>\n    internal long RankCount => Rb - Lb;\n\n    /// <summary>The node as the 32-bit navigator', F_NODE)
m('N02', 'EsaNode64.ToNode: the depth and the parent depth are swapped', NO, 'new(IntBridge.Length(Lb, member), IntBridge.Length(Rb, member), IntBridge.Length(Depth, member), IntBridge.Length(ParentDepth, member), IntBridge.Length(Pos, member), LastChild);', 'new(IntBridge.Length(Lb, member), IntBridge.Length(Rb, member), IntBridge.Length(ParentDepth, member), IntBridge.Length(Depth, member), IntBridge.Length(Pos, member), LastChild);', F_NODE)
m('N03', 'EsaNode64.Of: the last-child flag is lost', NO, 'new(node.Lb, node.Rb, node.Depth, node.ParentDepth, node.Pos, node.LastChild);', 'new(node.Lb, node.Rb, node.Depth, node.ParentDepth, node.Pos, false);', F_NODE)
m('N04', 'EsaNode64.ToNode: the position is narrowed without a check', NO, 'IntBridge.Length(Pos, member), LastChild);', 'unchecked((int)Pos), LastChild);', F_NODE, stage3='(' + T + 'Large.WideEngineSparseLargeTests)')

# the check: every old text occurs exactly once
bad = []
for x in M:
    text = open(ROOT + '/' + x['file']).read()
    n = text.count(x['old'])
    if n != 1:
        bad.append((x['id'], n))
if bad:
    print('NOT UNIQUE:', bad)
    sys.exit(1)
ids = [x['id'] for x in M]
assert len(ids) == len(set(ids))
json.dump(M, open(OUT, 'w'), indent=1)
print(len(M), 'mutants written to', OUT)
