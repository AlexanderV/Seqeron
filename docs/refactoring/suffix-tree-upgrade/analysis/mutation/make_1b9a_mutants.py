#!/usr/bin/env python3
"""Hand-made single-site mutants of the code that work package 1b.9a added or changed (prefix 1b9a_): the specialised pattern descent over a container (MappedSearch: the dispatch, the checks of the Checked mode, the
   64-bit window decode and the packed pattern, the LCP readers, the two child codecs), the tagged arm of the mapped child table, CHILD codec 2 in the header check, the builder, the section entry, BUILD_INFO, the
   Structure level and the e8-verify adapter. Each mutant is written to one JSON list per source file (1b9a_<file>_listed.json); every `old` must occur exactly once in its file and must change it (checked here).
   The runner (run_1b9a_mutants.py) builds the mutant in a separate worktree and runs the unit tests of the new code, and for a survivor the budgeted harness of 1b.8 (a short budget).
   usage: make_1b9a_mutants.py <repo root>"""
import json, os, sys

ROOT = sys.argv[1]
HERE = os.path.dirname(os.path.abspath(__file__))
SRC = 'src/SuffixTree/Algorithms/SuffixTree.Esa/'
APP = 'apps/SuffixTree.EmProto/'


def m(i, desc, old, new):
    return {'id': '1b9a_' + i, 'desc': desc, 'old': old, 'new': new}


FILES = {
    SRC + 'MappedSearch.cs': [
        m('M01', 'the access mode of the dispatch is swapped', "var lcp = new Lcp8Reader(_lcp.MappedBase, _lcp); return Run<CheckedAccess, Lcp8Reader, PackedCells<CheckedAccess, Lcp8Reader>>", "var lcp = new Lcp8Reader(_lcp.MappedBase, _lcp); return Run<TrustedAccess, Lcp8Reader, PackedCells<TrustedAccess, Lcp8Reader>>"),
        m('M02', 'the variant of a trusted file starts at 0 (the Checked loops run)', "(file.IsChecked ? 0 : 6)", "(file.IsChecked ? 0 : 0)"),
        m('M03', 'a 16-bit LCP file runs the 8-bit loop', "LcpKind.U16 => 1,", "LcpKind.U16 => 0,"),
        m('M04', 'a tagged child table runs the packed loop', "+ (file.Child.Mapped.IsTagged ? 1 : 0);", "+ 0;"),
        m('M05', 'the packed pattern has one word too few', "+ 63) / 64) + 2;", "+ 63) / 64) + 1;"),
        m('M06', 'the walk accepts a first l-index at the left bound (checked)', "if (TAccess.IsChecked && boundary >= 0 && (boundary <= a || boundary > rb)) ThrowDamaged();\n            int b;", "if (TAccess.IsChecked && boundary >= 0 && (boundary < a || boundary > rb)) ThrowDamaged();\n            int b;"),
        m('M07', 'the walk refuses an l-index at the right bound (checked)', "if (TAccess.IsChecked && boundary >= 0 && (boundary <= a || boundary > rb)) ThrowDamaged();\n            int b;", "if (TAccess.IsChecked && boundary >= 0 && (boundary <= a || boundary >= rb)) ThrowDamaged();\n            int b;"),
        m('M08', 'the l-indices that follow the first are not checked', "boundary = cells.Next(a, rb);\n                if (TAccess.IsChecked && boundary >= 0 && (boundary <= a || boundary > rb)) ThrowDamaged();", "boundary = cells.Next(a, rb);"),
        m('M09', 'the children are not left at the end of the keys of the interval', "if (key > symbol || boundary < 0)", "if (key > symbol)"),
        m('M10', 'a key that is larger or equal ends the search of the child', "if (key == symbol) break;", "if (key >= symbol) break;"),
        m('M11', 'a text position below zero is read (the clamp to the text is unsigned)', "int key = (uint)p < (uint)n ? Symbol(text, bits, p) : -1;", "int key = p < n ? Symbol(text, bits, p) : -1;"),
        m('M12', 'the leaf check does not compare the length of the rest with the suffix', "if ((uint)start > (uint)n || m - c > n - start || !Equal(", "if ((uint)start > (uint)n || !Equal("),
        m('M13', 'a leaf at the end of the text is refused', "if ((uint)start > (uint)n || m - c > n - start", "if ((uint)start >= (uint)n || m - c > n - start"),
        m('M14', 'the first l-index of a child may lie one past the table (checked)', "if (TAccess.IsChecked && (uint)first > (uint)n) ThrowDamaged();", "if (TAccess.IsChecked && (uint)first > (uint)n + 1) ThrowDamaged();"),
        m('M15', 'a child as deep as its parent is accepted (checked)', "(childEll <= ell || childEll > n)", "(childEll < ell || childEll > n)"),
        m('M16', 'a child deeper than the text is accepted (checked)', "(childEll <= ell || childEll > n)", "(childEll <= ell || childEll > n + 1)"),
        m('M17', 'the edge compare does not check that the stretch lies inside the text', "|| limit - c - 1 > n - (from + 1) || !Equal(", "|| !Equal("),
        m('M18', 'the edge is cut at the longer of the child value and the pattern', "int limit = Math.Min(childEll, m);", "int limit = Math.Max(childEll, m);"),
        m('M19', 'the edge compare leaves out one symbol', "if (limit > c + 1 && ((uint)from >= (uint)n", "if (limit > c + 2 && ((uint)from >= (uint)n"),
        m('M20', 'the edge compare starts at the first symbol of the edge, not behind it', "!Equal(text, bits, from + 1, pat, packed, c + 1, limit - c - 1, ref patternPacked)", "!Equal(text, bits, from, pat, packed, c, limit - c, ref patternPacked)"),
        m('M21', 'the interval of a match loses its last rank', "        lo = lb;\n        hi = rb;\n        return true;\n    }\n\n    /// <summary>The 64 stream bits", "        lo = lb;\n        hi = rb - 1;\n        return true;\n    }\n\n    /// <summary>The 64 stream bits"),
        m('M22', 'the window leaves out the carry byte', "return (w << s) | ((ulong)q[8] >> (8 - s));", "return w << s;"),
        m('M23', 'the symbol is cut one bit too far down', "(int)(Window(text, bits, i) >> (64 - bits));", "(int)(Window(text, bits, i) >> (63 - bits));"),
        m('M24', 'the stretch of a window is one symbol short', "int per = 64 / bits;\n        for (int done = 0; done < length; done += per)", "int per = 63 / bits;\n        for (int done = 0; done < length; done += per)"),
        m('M25', 'the compared window of the last stretch is one bit short', "if ((x >> (64 - (take * bits))) != 0) return false;", "if ((x >> (63 - (take * bits))) != 0) return false;"),
        m('M26', 'the second word of a pattern window is read for an aligned offset too', "return s == 0 ? high : high | (packed[w + 1] >> (64 - s));", "return high | (packed[w + 1] >> (64 - s));"),
        m('M27', 'a symbol that fits the word exactly is split', "if (room >= 0)\n            {\n                words[w] |= v << room;", "if (room > 0)\n            {\n                words[w] |= v << room;"),
        m('M28', 'the packed pattern is not cleared before it is filled', "        words.Clear();\n        long bit = 0;", "        long bit = 0;"),
        m('M29', 'the pattern is never packed', "if (!patternPacked)\n        {\n            PackPattern(pattern, packed, bits);", "if (patternPacked)\n        {\n            PackPattern(pattern, packed, bits);"),
        m('M30', 'the packed pattern lives on the heap from the first word on', "Span<ulong> packed = words <= 128 ? stackalloc ulong[words] : new ulong[words];", "Span<ulong> packed = words <= 1 ? stackalloc ulong[words] : new ulong[words];"),
        m('M31', 'the tagged cell is not masked to its rank', "public int First(int lb, int rb, bool lastChild) => (int)(cld[lastChild ? lb : rb] & ValueMask);", "public int First(int lb, int rb, bool lastChild) => (int)cld[lastChild ? lb : rb];"),
        m('M32', 'every tagged field with a rank is a next link', "return (cell & KindMask) == NextKind ? (int)(cell & ValueMask) : -1;", "return (cell & KindMask) != 0 ? (int)(cell & ValueMask) : -1;"),
        m('M33', 'the tagged kind is told by one bit', "return (cell & KindMask) == NextKind ?", "return (cell & NextKind) == NextKind ?"),
        m('M34', 'the packed walk does not stop at the right bound', "        if (k >= rb) return -1;\n        int lk = lcp[k];\n        if (lcp[k + 1] == lk) return k + 1;\n        int v = cld[k];\n        if (TAccess.IsChecked", "        int lk = lcp[k];\n        if (lcp[k + 1] == lk) return k + 1;\n        int v = cld[k];\n        if (TAccess.IsChecked"),
        m('M35', 'the packed walk takes a deeper l-index for the next one', "if (lcp[k + 1] == lk) return k + 1;", "if (lcp[k + 1] >= lk) return k + 1;"),
        m('M36', 'the packed next link is not validated against the table (checked)', "if (TAccess.IsChecked && (uint)v > (uint)n) MappedSearch.ThrowDamaged();\n        return lcp[v] == lk ? v : -1;", "return lcp[v] == lk ? v : -1;"),
        m('M37', 'the packed next link is taken without the LCP test', "return lcp[v] == lk ? v : -1;", "return v;"),
        m('M38', 'an escaped 8-bit LCP entry is returned as the escape', "return v != byte.MaxValue ? v : lcp[rank];", "return v;"),
        m('M39', 'an escaped 16-bit LCP entry is returned as the escape', "return v != ushort.MaxValue ? v : lcp[rank];", "return v;"),
        m('M40', 'a 32-bit LCP value above the table is accepted (checked)', "if (TAccess.IsChecked && (uint)v > (uint)length) MappedSearch.ThrowDamaged();", "if (TAccess.IsChecked && (uint)v > (uint)length + 1) MappedSearch.ThrowDamaged();"),
    ],
    SRC + 'ChildLayouts.cs': [
        m('C01', 'the tagged cell is not masked in FirstLIndex of the mapped table', "        return _tagged ? cell & TaggedValueMask : cell;", "        return cell;"),
        m('C02', 'every tagged field is a next link (mapped table)', "return (cell & TaggedKindMask) == TaggedNext ? (int)(cell & TaggedValueMask) : -1;", "return (cell & TaggedKindMask) != 0 ? (int)(cell & TaggedValueMask) : -1;"),
        m('C03', 'the tagged next link of the mapped table is not read as tagged', "        if (_tagged)\n        {\n            if (_checked && (uint)k >= (uint)_count) ThrowOutside(k);", "        if (!_tagged)\n        {\n            if (_checked && (uint)k >= (uint)_count) ThrowOutside(k);"),
        m('C04', 'a tagged cell outside the section is read in checked mode (the bound of NextLIndex)', "            if (_checked && (uint)k >= (uint)_count) ThrowOutside(k);\n            uint cell = (uint)_cld[k];", "            uint cell = (uint)_cld[k];"),
    ],
    SRC + 'ChildTable.cs': [
        m('T01', 'a mapped tagged table reports the packed layout', "mapped.IsTagged ? ChildLayoutKind.Tagged : ChildLayoutKind.Packed", "ChildLayoutKind.Packed"),
    ],
    SRC + 'External/Sesa/SesaFile.cs': [
        m('F01', 'the codec of the section is not given to the mapped table', "isChecked, parsed.Child.Codec == SesaFormat.ChildCodecTagged));", "isChecked, false));"),
    ],
    SRC + 'External/Sesa/SesaHeaderReader.cs': [
        m('H01', 'codec 2 is not supported by the header check', "if (codec is not (SesaFormat.ChildCodecPacked or SesaFormat.ChildCodecTagged)) throw", "if (codec is not SesaFormat.ChildCodecPacked) throw"),
        m('H02', 'the codec check lets an undefined codec 3 through', "if (codec is not (SesaFormat.ChildCodecPacked or SesaFormat.ChildCodecTagged)) throw", "if (codec is not (SesaFormat.ChildCodecPacked or SesaFormat.ChildCodecTagged or 3)) throw"),
        m('H03', 'the limit of the tagged ranks is one rank too small', "if (codec == SesaFormat.ChildCodecTagged && textLength > TaggedChild.MaxLength) throw", "if (codec == SesaFormat.ChildCodecTagged && textLength >= TaggedChild.MaxLength) throw"),
        m('H04', 'the limit of the tagged ranks is one rank too large', "if (codec == SesaFormat.ChildCodecTagged && textLength > TaggedChild.MaxLength) throw", "if (codec == SesaFormat.ChildCodecTagged && textLength > TaggedChild.MaxLength + 1L) throw"),
        m('H05', 'the entry of the child table is checked against codec 1 only', "CheckArray(child, \"CHILD\", rows, 4, child.Codec);", "CheckArray(child, \"CHILD\", rows, 4, 1);"),
    ],
    SRC + 'External/Sesa/SesaBuilder.cs': [
        m('B01', 'the default codec is the packed one', "ChildCodec codec = requested ?? (n <= TaggedChild.MaxLength ? ChildCodec.Tagged : ChildCodec.Packed);", "ChildCodec codec = requested ?? ChildCodec.Packed;"),
        m('B02', 'the default codec switches one text length too early', "requested ?? (n <= TaggedChild.MaxLength ? ChildCodec.Tagged : ChildCodec.Packed)", "requested ?? (n < TaggedChild.MaxLength ? ChildCodec.Tagged : ChildCodec.Packed)"),
        m('B03', 'a text that is too long for the tagged codec is accepted when it is requested', "if (codec == ChildCodec.Tagged && n > TaggedChild.MaxLength) throw new NotSupportedException($\"the tagged child table holds ranks in 30 bits (CHILD codec 2); a text of {n} symbols needs the packed codec (CHILD codec 1)\");\n        return codec;", "return codec;"),
        m('B04', 'the stage builds the packed codec whatever was chosen', "                Codec = childCodec,", "                Codec = ChildCodec.Packed,"),
        m('B05', 'the section entry says codec 1 whatever was built', "SesaFormat.KindChild => new SesaSection(kind, required, offset, length, rows, 4, childCodec == ChildCodec.Tagged ? SesaFormat.ChildCodecTagged : SesaFormat.ChildCodecPacked, 0, 0, crcs[i]),", "SesaFormat.KindChild => new SesaSection(kind, required, offset, length, rows, 4, SesaFormat.ChildCodecPacked, 0, 0, crcs[i]),"),
        m('B06', 'BUILD_INFO names the layout the other way around', "child=layout A, {(childCodec == ChildCodec.Tagged ? \"tagged\" : \"packed\")}", "child=layout A, {(childCodec == ChildCodec.Tagged ? \"packed\" : \"tagged\")}"),
        m('B07', 'BUILD_INFO is written with the default codec, not the built one', "string buildInfo = BuildInfoText(info, placed.Kind, placed.OverflowCount, options.IncludePsi, childCodec);", "string buildInfo = BuildInfoText(info, placed.Kind, placed.OverflowCount, options.IncludePsi, ChildCodecFor(n, null));"),
    ],
    SRC + 'External/Sesa/SectionLcpChildBuilder.cs': [
        m('S01', 'the stage accepts the tagged codec for a text beyond its ranks', "if (options.Codec == ChildCodec.Tagged && n > TaggedChild.MaxLength) throw new NotSupportedException(", "if (options.Codec == ChildCodec.Tagged && n > TaggedChild.MaxLength + (1L << 40)) throw new NotSupportedException("),
    ],
    SRC + 'External/Sesa/SesaVerifier.cs': [
        m('V01', 'the verifier expects the packed layout of a tagged file', "var sink = new ChildCompareSink(cld, n, written, ps.Child.Codec == SesaFormat.ChildCodecTagged);", "var sink = new ChildCompareSink(cld, n, written, false);"),
        m('V02', 'the verifier expects a next link without its tag', "if (_tagged) Expect(previous, NextKind | k);", "if (_tagged) Expect(previous, k);"),
        m('V03', 'the verifier expects the up kind for the down field', "int kind = lastChild ? DownKind : UpKind;", "int kind = lastChild ? UpKind : DownKind;"),
        m('V04', 'the verifier lets a tagged next link of the successor k + 1 go unchecked', "if (previous < 0) return;\n            if (_tagged) Expect(previous, NextKind | k);\n            else if (k != previous + 1) Expect(previous, k);", "if (previous < 0) return;\n            if (_tagged) { if (k != previous + 1) Expect(previous, NextKind | k); }\n            else if (k != previous + 1) Expect(previous, k);"),
    ],
    APP + 'Verify/OutputVerifier.cs': [
        m('A01', 'the adapter expects the packed layout of a tagged container', "source.ChildTagged ?? o.ChildTagged);", "o.ChildTagged);"),
        m('A02', 'the adapter expects a next link without its tag', "if (tagged) Expect(top.Last, NextTag | (uint)i);", "if (tagged) Expect(top.Last, (uint)i);"),
        m('A03', 'the adapter expects the up kind for the down field', "uint kind = lastChild ? DownTag : UpTag;", "uint kind = lastChild ? UpTag : DownTag;"),
    ],
    APP + 'Verify/VerifySources.cs': [
        m('A04', 'a container reports the packed codec whatever it holds', "public bool? ChildTagged => _file.Parsed.Child.Codec == SesaFormat.ChildCodecTagged;", "public bool? ChildTagged => false;"),
    ],
}

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
        if rel.startswith('apps/'):
            x['project'] = 'emproto'
    name = os.path.basename(rel).replace('.cs', '').lower()
    json.dump(ms, open(os.path.join(HERE, '1b9a_%s_listed.json' % name), 'w'), indent=1)
    total += len(ms)
print('listed', total, 'mutants in', len(FILES), 'files')
