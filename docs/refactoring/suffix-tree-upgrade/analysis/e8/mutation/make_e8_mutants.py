#!/usr/bin/env python3
"""Writes the hand-made single-site mutants of the verifier of E8 (apps/SuffixTree.EmProto/Verify and the E8 command glue), one JSON list per source file (e8_<file>_listed.json); every `old` pattern must occur exactly once in its file (checked here).
   usage: make_e8_mutants.py <repo root>"""
import json, os, sys

ROOT = sys.argv[1] if len(sys.argv) > 1 else '.'
V = 'apps/SuffixTree.EmProto/Verify/'
M = {}
def m(file, id, desc, old, new):
    M.setdefault(file, []).append({'id': id, 'desc': desc, 'old': old, 'new': new})

f = 'OutputVerifier.cs'
# ---- files and metadata
m(f, 'A01', 'the packed text is expected without its 64 padding bytes', 'Length("pk", PackedText.TotalBytes(n, bits));', 'Length("pk", PackedText.DataBytes(n, bits));')
m(f, 'A02', 'the suffix array file is expected to have n entries', 'Length("sa", 4 * (n + 1));', 'Length("sa", 4 * n);')
m(f, 'A03', 'the child table file is expected to have n entries', 'Length("cld", 4 * (n + 1));', 'Length("cld", 4 * n);')
m(f, 'A04', 'the LCP file is expected to have n entries', 'Length("lcp", width * (n + 1));', 'Length("lcp", width * n);')
m(f, 'A05', 'overflow pairs are expected to have 4 bytes', 'Length("ovf", 8 * overflow);', 'Length("ovf", 4 * overflow);')
m(f, 'A06', 'stats n greater than meta n is not reported by the files check', 'if (stats.N != n) filesError.Add', 'if (stats.N < n) filesError.Add')
m(f, 'A07', 'counts that add up to more than n pass the metadata check', 'if (m.Counts.Sum() != m.N)', 'if (m.Counts.Sum() < m.N)')
m(f, 'A08', 'a width above the one the alphabet needs passes the metadata check', 'if (m.Bits != Math.Max(1, BitsFor(m.Sigma)))', 'if (m.Bits < Math.Max(1, BitsFor(m.Sigma)))')
m(f, 'A09', 'equal neighbouring letters pass the metadata check', 'if (m.Letters[i] <= m.Letters[i - 1])', 'if (m.Letters[i] < m.Letters[i - 1])')
m(f, 'A10', 'a letter that never occurs passes the metadata check', 'if (m.Counts.Any(c => c < 1))', 'if (m.Counts.Any(c => c < 0))')
# ---- text
m(f, 'B01', 'the last padding byte is not checked', 'for (long i = dataBytes; i < file.Bytes; i++)', 'for (long i = dataBytes; i < file.Bytes - 1; i++)')
m(f, 'B02', 'the unused bits of the last data byte are not checked', 'if (usedBits != 0 && (file.Pointer[dataBytes - 1]', 'if (false && usedBits != 0 && (file.Pointer[dataBytes - 1]')
m(f, 'B03', 'the unused bits are counted modulo 4', 'int usedBits = (int)((n * bits) & 7);', 'int usedBits = (int)((n * bits) & 3);')
m(f, 'B04', 'a longer input passes the length check', 'if (length != text.Length)', 'if (length < text.Length)')
m(f, 'B05', 'the SHA-256 of the input is not compared with the metadata', 'if (sha != meta.Sha256)', 'if (sha == "")')
m(f, 'B06', 'a wrong symbol count is reported only when it is too low', 'if (total[s] != expected)', 'if (total[s] < expected)')
# ---- the rank pass
m(f, 'C01', 'a duplicate suffix array entry is not noticed', 'if ((old & mask) != 0) totals.Permutation.Add', 'if ((old & mask) == 7) totals.Permutation.Add')
m(f, 'C02', 'the bitmap uses 32 bits per word index but 64 bits per mask', 'ulong old = Interlocked.Or(ref bits[q >> 6], mask);', 'ulong old = Interlocked.Or(ref bits[q >> 5], mask);')
m(f, 'C03', 'the bit of an entry is q modulo 32', 'ulong mask = 1UL << (int)(q & 63);', 'ulong mask = 1UL << (int)(q & 31);')
m(f, 'C04', 'an entry above n is not reported', 'if (q > n) totals.Permutation.Add', 'if (q > n + 1) totals.Permutation.Add')
m(f, 'C05', 'a value of exactly 255 is not counted as an escape', 'if (stored >= 255) local.Escapes8++;', 'if (stored > 255) local.Escapes8++;')
m(f, 'C06', 'a value of exactly 65535 is not counted', 'if (stored >= 65535) local.Escapes16++;', 'if (stored > 65535) local.Escapes16++;')
m(f, 'C07', 'SA[0] is not compared with n', 'if (q != n) totals.Permutation.Add(0,', 'if (q > n) totals.Permutation.Add(0,')
m(f, 'C08', 'LCP[0] is not checked', 'if (stored != 0) totals.Lcp.Add(0,', 'if (stored > 1_000_000_000) totals.Lcp.Add(0,')
m(f, 'C09', 'the common prefix is bounded by the shorter start instead of the shorter suffix', 'long maxLen = n - Math.Max(p, q);', 'long maxLen = n - Math.Min(p, q);')
m(f, 'C10', 'a cap of 0 means a limit of 0', 'long lim = cap > 0 ? Math.Min(maxLen, cap) : maxLen;', 'long lim = cap >= 0 ? Math.Min(maxLen, cap) : maxLen;')
m(f, 'C11', 'a pair where one suffix ends is treated as a long pair', 'if (l < lim || lim == maxLen)', 'if (l < lim)')
m(f, 'C12', 'the suffix that starts earlier is the shorter one', 'bool ordered = l == maxLen ? p > q :', 'bool ordered = l == maxLen ? p < q :')
m(f, 'C13', 'equal symbols at the mismatch pass the order check', 'text.Symbol(p + l) < text.Symbol(q + l);', 'text.Symbol(p + l) <= text.Symbol(q + l);')
m(f, 'C14', 'a stored LCP that is too large is not reported (exact path)', 'if (stored != l) totals.Lcp.Add', 'if (stored < l) totals.Lcp.Add')
m(f, 'C15', 'a stored LCP above the shorter suffix is accepted for a long pair', 'if (stored < cap || stored > maxLen)', 'if (stored < cap)')
m(f, 'C16', 'for a long pair that ends at the shorter suffix the order is reversed', 'if (stored == maxLen) ok = p > q;', 'if (stored == maxLen) ok = p < q;')
m(f, 'C17', 'for a long pair the order of the mismatch symbols is reversed', 'ok = a < b;', 'ok = a > b;')
m(f, 'C18', 'no probes of the matching part of a long pair', 'for (int k = 0; k < o.Probes && span > 0; k++)', 'for (int k = 0; k < 0 && span > 0; k++)')
m(f, 'C19', 'the first probe is two symbols before the end instead of one', 'long off = k == 0 ? stored - 1 :', 'long off = k == 0 ? stored - 2 :')
m(f, 'C20', 'the probe compares for inequality as equality', 'if (text.Symbol(p + off) != text.Symbol(q + off))', 'if (text.Symbol(p + off) == text.Symbol(q + off))')
m(f, 'C21', 'a pair value equal to the escape is refused', 'lcp.ValueOfPair(pair) >= lcp.Escape)', 'lcp.ValueOfPair(pair) > lcp.Escape)')
m(f, 'C22', 'the first pair of a block is looked up from the start of the list', 'long pair = lcp.HasEscape && lcp.PairCount > 0 ? lcp.LowerBound(first) : 0;', 'long pair = 0;')
m(f, 'C23', 'a pair left over in a block is not reported', 'if (lcp.HasEscape && lcp.PairCount > 0 && pair != lcp.LowerBound(end)) totals.Lcp.Add(', 'if (false && lcp.HasEscape && lcp.PairCount > 0 && pair != lcp.LowerBound(end)) totals.Lcp.Add(')
m(f, 'C24', 'the width rule takes 17 times the escapes', 'if (16 * escapes8 <= n) return LcpKind.U8;', 'if (17 * escapes8 <= n) return LcpKind.U8;')
m(f, 'C25', 'the width rule for U16 takes 9 times the escapes', 'if (8 * escapes16 <= n) return LcpKind.U16;', 'if (9 * escapes16 <= n) return LcpKind.U16;')
# ---- the scan
m(f, 'D01', 'the left boundary of a new interval is the rank itself', 'uint left = (uint)(i - 1);', 'uint left = (uint)i;')
m(f, 'D02', 'a child is the last child when the next value equals the parent', 'Close(e.Ell, e.Lb, i - 1, e.First, cur < stack.Top.Ell);', 'Close(e.Ell, e.Lb, i - 1, e.First, cur <= stack.Top.Ell);')
m(f, 'D03', 'an interval of the same value is pushed again', 'if (top.Ell < cur)\n            {\n                if (!stack.TryPush', 'if (top.Ell <= cur)\n            {\n                if (!stack.TryPush')
m(f, 'D04', 'the next link is expected when the successor is two ranks away', 'if (i != top.Last + 1L) Expect(top.Last, (uint)i);', 'if (i != top.Last + 2L) Expect(top.Last, (uint)i);')
m(f, 'D05', 'up and down fields are swapped', 'Expect(lastChild ? lb : rb, first);', 'Expect(lastChild ? rb : lb, first);')
m(f, 'D06', 'the deepest interval is the last one among equals', 'if (ell > r.LrsDepth)', 'if (ell >= r.LrsDepth)')
m(f, 'D07', 'the LRS rank is the right boundary', 'r.LrsRank = lb;', 'r.LrsRank = (uint)rb;')
m(f, 'D08', 'the first maximum is the last maximum', 'if (cur > r.MaxLcp)', 'if (cur >= r.MaxLcp)')
m(f, 'D09', 'the first maximum rank is off by one', 'r.FirstMaxRank = i;', 'r.FirstMaxRank = i - 1;')
m(f, 'D10', 'cells holding 1 are not counted as set', 'if (cld[k] != 0) count++;', 'if (cld[k] > 1) count++;')
m(f, 'D11', 'the root is closed as a non-last child', 'Close(e.Ell, e.Lb, n, e.First, true);', 'Close(e.Ell, e.Lb, n, e.First, false);')
m(f, 'D12', 'the internal count includes the root', 'if (ell > 0)\n            {\n                r.Internal++;', 'if (ell >= 0)\n            {\n                r.Internal++;')
m(f, 'D13', 'the node count omits the empty suffix', 'var expected = new StatsInfo(1 + scan.Internal + n + 1,', 'var expected = new StatsInfo(1 + scan.Internal + n,')
m(f, 'D14', 'a text with the stack bound exceeded is reported as passed', 'r.Aborted = $"the stack of open intervals would exceed', 'r.Aborted = "" + $"the stack of open intervals would exceed')
m(f, 'D15', 'the hashes are of the wrong file extension list order (sha of pk for stats)', 'string sha = HashFile(o.Prefix + "." + e, o.IoBufferBytes);', 'string sha = HashFile(o.Prefix + "." + (e == "stats" ? "pk" : e), o.IoBufferBytes);')
m(f, 'D16', 'the hash reads only the first buffer of a file', 'at += got;\n        }\n        return Convert.ToHexStringLower', 'at += length;\n        }\n        return Convert.ToHexStringLower')

f = 'MetaFiles.cs'
m(f, 'E01', 'a statistics file without a final line feed passes', "if (!text.EndsWith('\\n'))", 'if (false)')
m(f, 'E02', 'a file with more lines than keys passes', 'if (lines.Length != keys.Length)', 'if (lines.Length < keys.Length)')
m(f, 'E03', 'a line that only contains its key passes', 'if (!lines[i].StartsWith(prefix, StringComparison.Ordinal))', 'if (!lines[i].Contains(keys[i], StringComparison.Ordinal))')
m(f, 'E04', 'numbers with a leading zero are accepted', '(s.Length > 1 && s[0] == \'0\') || ', '')
m(f, 'E05', 'letters_hex of the wrong length passes', 'v[3].Length != 2 * sigma', 'v[3].Length < 2 * sigma')
m(f, 'E06', 'a counts list of the wrong length passes', 'if (counts.Length != sigma ||', 'if (counts.Length < sigma ||')
m(f, 'E07', 'a SHA-256 of the wrong length passes', 'v[5].Length != 64', 'v[5].Length < 64')
m(f, 'E08', 'lower-case letters_hex passes', "c is >= '0' and <= '9' or >= 'A' and <= 'F'))\n        {\n            error = \"sigma, bits", "c is >= '0' and <= '9' or >= 'A' and <= 'F' or >= 'a' and <= 'f'))\n        {\n            error = \"sigma, bits")
m(f, 'E09', 'upper-case hash digits pass', "v[5].All(c => c is >= '0' and <= '9' or >= 'a' and <= 'f')", "v[5].All(c => c is >= '0' and <= '9' or >= 'a' and <= 'f' or >= 'A' and <= 'F')")
m(f, 'E10', 'a sigma of 257 passes', 'sigma is < 1 or > 256', 'sigma is < 1 or > 257')

f = 'LcpView.cs'
m(f, 'F01', 'the 2-byte base is read as 4 bytes', '2 => ((ushort*)basePointer)[rank],', '2 => ((uint*)basePointer)[rank],')
m(f, 'F02', 'the 1-byte escape is 65535', '1 => byte.MaxValue,', '1 => ushort.MaxValue,')
m(f, 'F03', 'the lower bound skips a pair whose rank equals the key', 'if (pairs[2 * mid] < rank) lo = mid + 1;', 'if (pairs[2 * mid] <= rank) lo = mid + 1;')
m(f, 'F04', 'the lower bound never sees the last pair', 'long lo = 0, hi = PairCount;', 'long lo = 0, hi = PairCount - 1;')
m(f, 'F05', 'the value of a pair is its rank', 'internal uint ValueOfPair(long k) => pairs[(2 * k) + 1];', 'internal uint ValueOfPair(long k) => pairs[2 * k];')
m(f, 'F06', 'the 4-byte base has an escape of 0', '2 => ushort.MaxValue, _ => uint.MaxValue };', '2 => ushort.MaxValue, _ => 0 };')

f = 'IntervalStack.cs'
m(f, 'G01', 'the capacity is one entry too large', 'if (Count >= Capacity) return false;', 'if (Count > Capacity) return false;')
m(f, 'G02', 'the peak is not tracked', 'if (Count > Peak) Peak = Count;', 'if (Count > Peak + 1000000) Peak = Count;')
m(f, 'G03', 'a pop removes two entries', 'internal void Pop() => Count--;', 'internal void Pop() => Count -= 2;')
m(f, 'G04', 'the block index of a push uses the wrong shift', 'int block = (int)(Count >> _blockBits);', 'int block = (int)(Count >> (_blockBits + 1));')

f = 'E8.cs'
m(f, 'H01', 'the exit code of a failed verification is 0', 'return res.AllPassed ? 0 : 1;', 'return 0;')
m(f, 'H02', 'the whole check is reported as pass when one is skipped', 'r.Row("E8", c, "verify_all", res.AllPassed ? "pass" : "FAIL"', 'r.Row("E8", c, "verify_all", res.Checks.All(x => x.Status != CheckStatus.Fail) ? "pass" : "FAIL"')
m(f, 'H03', 'a skipped check is reported as pass', 'CheckStatus.Fail => "FAIL", _ => "skipped" }', 'CheckStatus.Fail => "FAIL", _ => "pass" }')
m(f, 'H04', 'the default LCE cap is 0', 'LceCap = o.Long("lce-cap", 65536),', 'LceCap = o.Long("lce-cap", 0),')
m(f, 'H05', 'the stack cap is given in entries instead of Mi entries', 'StackCapEntries = o.Long("stack-cap-m", 64) << 20,', 'StackCapEntries = o.Long("stack-cap-m", 64),')
m(f, 'H06', 'no-child still checks and no-hash still hashes (flags ignored)', 'Hash = !o.Flag("no-hash"),', 'Hash = true,')
m(f, 'H07', 'the sums file lists the hashes in reverse order', 'new[] { "pk", "pk.meta", "sa", "lcp", "ovf", "cld", "stats" }.Where(res.Sha256.ContainsKey)', 'new[] { "pk", "pk.meta", "sa", "lcp", "ovf", "cld", "stats" }.Reverse().Where(res.Sha256.ContainsKey)')

total = 0
for file, items in M.items():
    path = os.path.join(ROOT, V + file) if file != 'E8.cs' else os.path.join(ROOT, 'apps/SuffixTree.EmProto/Experiments/E8.cs')
    src = open(path).read()
    for it in items:
        n = src.count(it['old'])
        if n != 1:
            print('PATTERN-%s %s %s' % ('MISSING' if n == 0 else 'NOT-UNIQUE', file, it['id']))
        it['file'] = os.path.relpath(path, ROOT)
    total += len(items)
    json.dump(items, open('e8_%s_listed.json' % file.replace('.cs', '').lower(), 'w'), indent=1)
print('%d mutants in %d files' % (total, len(M)))
