#!/usr/bin/env python3
"""Hand-made single-site mutants of the production engine for the harness audit of work package 1b.8 (prefix 1b8_): the search and analysis arms of EsaIndex, the traversal, the contract index, the navigator, the navigation
   pack, the mapped (file-backed) arms of 1b.7f, the validation of the container. Each mutant is written to one JSON list per source file (1b8_<file>_listed.json); every `old` must occur exactly once in its file and must change it
   (checked here). The runner (run_1b8_mutants.py) builds the mutant in a separate worktree and runs the budgeted harness against it (a short budget, 60-90 s).
   usage: make_1b8_mutants.py <repo root>
   environment: ONLY = a file with ids to list (the survivors of a round), ROUND = a prefix such as r2 for the names of the lists"""
import json, os, sys

ROOT = sys.argv[1]
ONLY = set(open(os.environ['ONLY']).read().split()) if os.environ.get('ONLY') else None
ROUND = os.environ.get('ROUND', '')
HERE = os.path.dirname(os.path.abspath(__file__))
SRC = 'src/SuffixTree/Algorithms/SuffixTree.Esa/'


def m(i, desc, old, new):
    return {'id': '1b8_' + i, 'desc': desc, 'old': old, 'new': new}


FILES = {
    'EsaIndex.Search.cs': [
        m('S01', 'the long overload of Locate ignores the skip', "        int start = (int)(lo + skip);\n        SaView sa = Sa;", "        int start = lo;\n        SaView sa = Sa;"),
        m('S02', 'the descent of the in-memory trees does not stop at the end of the keys of an interval',
          "int key = p < n ? ToInt(text[p]) : -1;\n                if (key == symbol) break;\n                if (key > symbol || boundary < 0)", "int key = p < n ? ToInt(text[p]) : -1;\n                if (key == symbol) break;\n                if (key > symbol)"),
        m('S03', 'the edge check of the descent leaves out one symbol', "if (limit > c + 1 && !text.Slice(from + 1, limit - c - 1).SequenceEqual(pat.Slice(c + 1, limit - c - 1)))", "if (limit > c + 2 && !text.Slice(from + 1, limit - c - 1).SequenceEqual(pat.Slice(c + 1, limit - c - 1)))"),
        m('S04', 'a pattern longer than the rest of the suffix is sliced', "if (m - c > n - start || !text.Slice(start, m - c).SequenceEqual(pat.Slice(c)))", "if (!text.Slice(start, m - c).SequenceEqual(pat.Slice(c)))"),
        m('S05', 'the mapped descent translates the pattern to a wrong symbol', "            pat[i] = (byte)s;", "            pat[i] = (byte)(s ^ 1);"),
        m('S06', 'the mapped descent accepts a child as deep as its parent', "if (isChecked && (childEll <= ell || childEll > n)) ThrowDamagedChildTable();", "if (isChecked && (childEll < ell || childEll > n)) ThrowDamagedChildTable();"),
        m('S07', 'TryMatch does not refuse the empty text', "        if (Length == 0) return false;\n        if (IsMapped) return MatchMapped", "        if (IsMapped) return MatchMapped"),
        m('S08', 'an edge that ends exactly at the end of the text is refused (mapped)', "if ((uint)start > (uint)text.Length || expected.Length > text.Length - start) return false;", "if ((uint)start > (uint)text.Length || expected.Length >= text.Length - start) return false;"),
        m('S09', 'the mapped descent does not stop at a leaf mismatch of the first symbol', "                if ((uint)start > (uint)n || m - c > n - start || !EdgeEquals(text, start, pat.Slice(c), edge))", "                if ((uint)start > (uint)n || !EdgeEquals(text, start, pat.Slice(c), edge))"),
    ],
    'EsaIndex.Analysis.cs': [
        m('A01', 'a longest repeat that starts at rank 1 is dropped', "        if (runStart >= 0) results.Add(Tie(runStart, Length, length));", "        if (runStart > 0) results.Add(Tie(runStart, Length, length));"),
        m('A02', 'the number of distinct substrings is n(n-1)/2 minus the LCP sum', "long total = n * (n + 1) / 2;", "long total = n * (n - 1) / 2;"),
        m('A03', 'distinct substrings by length skip the LCP value 1', "if (value > 0) counts[Math.Min(value, m)]++;", "if (value > 1) counts[Math.Min(value, m)]++;"),
        m('A04', 'distinct substrings by length: the count of length i is off by one', "counts[i] = n - i + 1 - atLeast;", "counts[i] = n - i - atLeast;"),
        m('A05', 'distinct substrings by length: the empty substring is not counted', "        counts[0] = 1;\n        return counts;", "        counts[0] = 0;\n        return counts;"),
        m('A06', 'GetAllSuffixes leaves the last suffix out', "for (int rank = 1; rank <= Length; rank++) suffixes[rank - 1] = TextView.Substring(sa[rank]);", "for (int rank = 1; rank < Length; rank++) suffixes[rank - 1] = TextView.Substring(sa[rank]);"),
        m('A07', 'the positions of a longest repeat are not sorted', "        Array.Sort(positions);\n        return (TextView.Substring(positions[0], length), positions);", "        return (TextView.Substring(positions[0], length), positions);"),
        m('A08', 'the longest repeats are not ordered by first occurrence', "        results.Sort(static (a, b) => a.Positions[0].CompareTo(b.Positions[0]));\n", ""),
        m('A09', 'EnumerateSuffixes leaves the last suffix out', "for (int rank = 1; rank <= Length; rank++) yield return TextView.Substring(sa[rank]);", "for (int rank = 1; rank < Length; rank++) yield return TextView.Substring(sa[rank]);"),
    ],
    'EsaIndex.Matching.cs': [
        m('M01', 'MEM accepts a minimum length below 1', "        ArgumentNullException.ThrowIfNull(query);\n        ArgumentOutOfRangeException.ThrowIfLessThan(minLength, 1);\n        EsaNavigator nav = CreateNavigator();\n        return SuffixTreeAlgorithms.FindMaximalExactMatches<",
          "        ArgumentNullException.ThrowIfNull(query);\n        EsaNavigator nav = CreateNavigator();\n        return SuffixTreeAlgorithms.FindMaximalExactMatches<"),
        m('M02', 'MUM accepts an undefined uniqueness', "        if (uniqueness is not (MumUniqueness.Both or MumUniqueness.Reference))\n            throw new ArgumentOutOfRangeException(nameof(uniqueness), uniqueness, \"Undefined MumUniqueness value.\");\n", ""),
        m('M03', 'the info of no common substring reports position 0', "if (results.PositionsInText.Count == 0) return (string.Empty, -1, -1);", "if (results.PositionsInText.Count == 0) return (string.Empty, 0, 0);"),
        m('M04', 'the anchors accept a null query (the reference check is left out)', "        ArgumentNullException.ThrowIfNull(query);\n        EsaNavigator nav = CreateNavigator();\n        return SuffixTreeAlgorithms.FindExactMatchAnchors<", "        EsaNavigator nav = CreateNavigator();\n        return SuffixTreeAlgorithms.FindExactMatchAnchors<"),
    ],
    'EsaIndex.Traverse.cs': [
        m('T01', 'the root reports n leaves', "visitor.VisitNode(0, 0, n + 1, CountChildren(firstOfRoot, n), 0);", "visitor.VisitNode(0, 0, n, CountChildren(firstOfRoot, n), 0);"),
        m('T02', 'a leaf reports a wrong parent depth', "visitor.VisitNode(p, -1, 1, 0, f.Ell);", "visitor.VisitNode(p, -1, 1, 0, f.Ell + 1);"),
        m('T03', 'the terminator edge has the key 0', "visitor.EnterBranch((uint)p < (uint)n ? text[p] : -1);", "visitor.EnterBranch((uint)p < (uint)n ? text[p] : 0);"),
        m('T04', 'an internal edge records the first suffix instead of the leftmost occurrence', "int occ = occurrence[first];", "int occ = sa[a];"),
        m('T05', 'the repeated-leaf check of the file-backed Traverse is removed (reverts the fix of 1b.8)', "                if (seen is not null) MarkLeaf(seen, sa[a], n);\n", ""),
        m('T06', 'the mapped walk does not check that a child interval is not empty', "if (mapped && (b < a || b > f.Rb || (next >= 0 && (next <= f.NextBoundary || next > f.Rb)))) throw DamagedChildTable();", "if (mapped && (b > f.Rb || (next >= 0 && (next <= f.NextBoundary || next > f.Rb)))) throw DamagedChildTable();"),
        m('T07', 'the mapped walk accepts a child as deep as its parent', "if (mapped && (ell <= f.Ell || ell > n)) throw DamagedChildTable();", "if (mapped && (ell < f.Ell || ell > n)) throw DamagedChildTable();"),
        m('T08', 'the first l-index at the left bound of an interval is accepted', "if (first <= lb || first > rb) throw DamagedChildTable();", "if (first < lb || first > rb) throw DamagedChildTable();"),
        m('T09', 'a suffix-array entry below zero reaches the bitmap of the repeated-leaf check', "if ((uint)position > (uint)n) throw new InvalidDataException(\"the suffix array holds a position outside the text", "if (position > n) throw new InvalidDataException(\"the suffix array holds a position outside the text"),
    ],
    'LeftmostOccurrence.cs': [
        m('O01', 'an LCP[1] other than 0 is not refused (reverts a fix of 1b.8)', "                if (first[sp] < 0) throw DamagedLcp(); // only an LCP[1] other than 0 leaves the base interval without a first l-index\n", ""),
        m('O02', 'a negative LCP value is not refused (reverts a fix of 1b.8)', "            if (sp == 0) throw DamagedLcp(); // a negative LCP value popped the base interval\n", ""),
        m('O03', 'the occurrence of an interval is the largest position, not the leftmost', "if (sa[i - 1] < min[sp - 1]) min[sp - 1] = sa[i - 1];", "if (sa[i - 1] > min[sp - 1]) min[sp - 1] = sa[i - 1];"),
    ],
    'EsaSuffixIndex.cs': [
        m('X02', 'Locate of a locus includes the empty suffix of rank 0', "int from = (int)Math.Max(locus.Lo, 1);", "int from = (int)locus.Lo;"),
        m('X03', 'the extension of a point inside an edge compares one symbol beyond the text', "if (at >= (uint)_n || text[(int)at] != symbol) return false;", "if (at > (uint)_n || text[(int)at] != symbol) return false;"),
        m('X04', 'a locus with an aux one beyond the text is refused', "locus.Aux >= locus.Length && locus.Aux <= n + 1;", "locus.Aux >= locus.Length && locus.Aux <= n;"),
        m('X05', 'the last-child flag of a node takes the equal side', "bool lastChild = left > right;", "bool lastChild = left >= right;"),
        m('X06', 'the rank n is refused', "if ((ulong)rank > (ulong)_n) throw new ArgumentOutOfRangeException(parameterName, rank,", "if ((ulong)rank >= (ulong)_n) throw new ArgumentOutOfRangeException(parameterName, rank,"),
        m('X07', 'a locus whose parent is as deep as the locus is accepted', "bool valid = parentDepth < locus.Length;", "bool valid = parentDepth <= locus.Length;"),
    ],
    'EsaSuffixIndex.Arrays.cs': [
        m('X01', 'Lcp of a rank with itself is one too long', "if (rankA == rankB) return _n - _index.Sa[(int)rankA];", "if (rankA == rankB) return _n - _index.Sa[(int)rankA] + 1;"),
        m('X08', 'Lca of two ranks of one text position has the wrong depth', "long depth = _n - _index.Sa[rank] + 1;", "long depth = _n - _index.Sa[rank];"),
        m('X09', 'Lca does not return the root for an LCP of 0', "        if (value == 0) return Root;\n        return new IndexLocus(pyramid.PrevLess", "        return new IndexLocus(pyramid.PrevLess"),
        m('X10', 'Lca includes one rank too many on the right', "pyramid.NextLess(hi + 1, value) - 1, value, value);", "pyramid.NextLess(hi + 1, value), value, value);"),
    ],
    'EsaNavigator.cs': [
        m('N01', 'CollectLeaves includes the rank of the empty suffix', "int from = Math.Max(node.Lb, 1);", "int from = node.Lb;"),
        m('N02', 'the child lookup does not stop at the end of the keys', "            if (symbol > key || boundary < 0)\n            {\n                child = EsaNode.Null;\n                return false;\n            }\n            a = boundary;\n            boundary = table.NextLIndex(boundary, rb);", "            if (symbol > key)\n            {\n                child = EsaNode.Null;\n                return false;\n            }\n            a = boundary;\n            boundary = table.NextLIndex(boundary, rb);"),
        m('N03', 'the edge symbol ignores the parent depth', "uint index = (uint)node.Pos + (uint)node.ParentDepth + (uint)offset;", "uint index = (uint)node.Pos + (uint)offset;"),
    ],
    'EsaNavigationPack.cs': [
        m('P01', 'the suffix link of a node of depth 2 is the root', "        if (t == 0) return Root;                                  // an internal node of depth 1", "        if (t <= 1) return Root;                                  // an internal node of depth 1"),
        m('P02', 'the suffix link reads the inverse suffix array one rank off', "int a = isa[node.Pos + 1];", "int a = isa[node.Pos];"),
        m('P03', 'the suffix link flags the last child on a tie', "return new EsaNode(left, right - 1, t, parentDepth, sa[left], leftValue > rightValue);", "return new EsaNode(left, right - 1, t, parentDepth, sa[left], leftValue >= rightValue);"),
        m('P04', 'the pack does not refuse a suffix array that repeats a position (mapped)', "if ((uint)position >= (uint)count || isaBuf[position] != -1)", "if ((uint)position >= (uint)count)"),
    ],
    'SaView.cs': [
        m('V01', 'a copy that ends exactly at the end of the suffix array is refused', "if (start < 0 || destination.Length > _count - start) ThrowOutside(start);", "if (start < 0 || destination.Length >= _count - start) ThrowOutside(start);"),
    ],
    'EsaText.cs': [
        m('E01', 'a substring that ends exactly at the end of the packed text is refused', "if (start < 0 || length < 0 || start > _length - length) ThrowOutside(start);", "if (start < 0 || length < 0 || start >= _length - length) ThrowOutside(start);"),
    ],
    'LcpPyramid.cs': [
        m('L01', 'the previous smaller value of position 0 is position 0', "        if (pos < 0) return -1;\n        ArgumentOutOfRangeException.ThrowIfGreaterThanOrEqual(pos, _length);", "        if (pos <= 0) return -1;\n        ArgumentOutOfRangeException.ThrowIfGreaterThanOrEqual(pos, _length);"),
        m('L02', 'the next smaller value of a position behind the end is scanned', "        if (pos >= _length) return _length;\n        ArgumentOutOfRangeException.ThrowIfNegative(pos);", "        if (pos > _length) return _length;\n        ArgumentOutOfRangeException.ThrowIfNegative(pos);"),
    ],
    'ChildLayouts.cs': [
        m('C01', 'the mapped child table: the next l-index is always the stored cell', "        int v = _cld[k];\n        return _lcp[v] == lk ? v : -1;\n    }\n\n    private void ThrowOutside(int cell) => throw new InvalidDataException($\"the child table is read at cell", "        int v = _cld[k];\n        return v;\n    }\n\n    private void ThrowOutside(int cell) => throw new InvalidDataException($\"the child table is read at cell"),
        m('C02', 'the mapped child table: the walk continues at the right bound', "        if (k >= rb) return -1;\n        if (_checked && (uint)k >= (uint)_count) ThrowOutside(k);", "        if (k > rb) return -1;\n        if (_checked && (uint)k >= (uint)_count) ThrowOutside(k);"),
    ],
    'External/Sesa/SesaFile.cs': [
        m('F01', 'the Structure level is never checked', "if (options.Validation >= EsaValidationLevel.Structure) SesaVerifier.CheckStructure(file);", "if (options.Validation > EsaValidationLevel.Structure) SesaVerifier.CheckStructure(file);"),
        m('F02', 'the Sections level is never checked at Sections', "if (options.Validation >= EsaValidationLevel.Sections) SesaVerifier.CheckSections(file, hardwareCrc);", "if (options.Validation > EsaValidationLevel.Sections) SesaVerifier.CheckSections(file, hardwareCrc);"),
    ],
    'External/Sesa/SesaHeaderReader.cs': [
        m('H01', 'a FILE_SIZE smaller than the file is accepted (trailing bytes)', "if (fileSize != length) throw new InvalidDataException", "if (fileSize > length) throw new InvalidDataException"),
        m('H02', 'the header checksum is not compared', "if (crc != stored) throw new InvalidDataException", "if (crc != stored && false) throw new InvalidDataException"),
        m('H03', 'the reserved header word may be anything', "if (U32(p, 120) != 0) throw new InvalidDataException(\"the reserved header word is not zero\");\n", ""),
        m('H04', 'COUNTS may sum to less than the text', "if (sum != (ulong)n) throw new InvalidDataException($\"COUNTS sums to {sum}, not to the {n} symbols of the text\");", "if (sum > (ulong)n) throw new InvalidDataException($\"COUNTS sums to {sum}, not to the {n} symbols of the text\");"),
    ],
    'External/Sesa/SesaVerifier.cs': [
        m('W01', 'an LCP value that is smaller than the true one is accepted', "if (stored != h) throw new InvalidDataException($\"LCP[{r}] is {stored}, but the suffixes", "if (stored > h) throw new InvalidDataException($\"LCP[{r}] is {stored}, but the suffixes"),
        m('W02', 'NODE_COUNT may be larger than the true one', "if (h.NodeCount != nodeCount) throw new InvalidDataException", "if (h.NodeCount < nodeCount) throw new InvalidDataException"),
        m('W03', 'LRS_DEPTH is not compared', "if (h.LrsDepth != sink.LrsDepth) throw new InvalidDataException", "if (h.LrsDepth != sink.LrsDepth && false) throw new InvalidDataException"),
        m('W04', 'the symbols of the text are not compared with SIGMA', "if (s >= sigma) throw new InvalidDataException", "if (s > sigma) throw new InvalidDataException"),
    ],
    'SuffixIndexFactory.cs': [
        m('Y01', 'the dispatcher passes no options to the container reader', "return EsaSuffixTreeFactory.Open(path, options);", "return EsaSuffixTreeFactory.Open(path, null);"),
    ],
    'EsaSuffixTree.cs': [
        m('Z01', 'the facade forwards the span overload of the count to the string search of the value', "        return _index.CountOccurrences(pattern);\n    }\n\n    // ------------------------------------------------------------------ ISuffixTreeAnalysis", "        return _index.CountOccurrences(pattern) + (pattern.Length == 7 ? 1 : 0);\n    }\n\n    // ------------------------------------------------------------------ ISuffixTreeAnalysis"),
        m('Z02', 'the facade drops the uniqueness of MUM', "return _index.FindMaximalUniqueMatches(query, minLength, uniqueness);", "return _index.FindMaximalUniqueMatches(query, minLength);"),
        m('Z03', 'the node count is one too large', "return CheckedNodeCount(_statistics.NodeCount);", "return CheckedNodeCount(_statistics.NodeCount + 1);"),
    ],
}

bad = 0
total = 0
for f, lst in FILES.items():
    text = open(os.path.join(ROOT, SRC + f)).read()
    out = []
    for mm in lst:
        if ONLY is not None and mm['id'] not in ONLY:
            continue
        c = text.count(mm['old'])
        if c != 1:
            print('PATTERN', f, mm['id'], 'occurs', c)
            bad += 1
        if mm['old'] == mm['new']:
            print('UNCHANGED', f, mm['id'])
            bad += 1
        out.append({'id': mm['id'], 'desc': mm['desc'], 'old': mm['old'], 'new': mm['new'], 'file': SRC + f})
    total += len(out)
    name = f.replace('/', '_').replace('.cs', '').lower()
    if ONLY is not None and not out:
        continue
    json.dump(out, open(os.path.join(HERE, '1b8_%s%s_listed.json' % (ROUND + '_' if ROUND else '', name)), 'w'), indent=1)
print('listed', total, 'mutants in', len(FILES), 'files;', bad, 'problems')
