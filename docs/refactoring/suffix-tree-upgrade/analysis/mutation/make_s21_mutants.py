#!/usr/bin/env python3
"""The hand-made mutants of work package S2.1 (prefix s21_): single-site changes of the test-support code of tests/SuffixTree/SuffixTree.Esa.Tests/Large (the closed forms, the CRC combine, the sparse container
   writer, the 5-byte re-encoder and its reader). Writes s21_listed.json (the list that run_s21_mutants.py runs) and checks that every `old` text occurs exactly once in its file.
   usage: make_s21_mutants.py <repo root> <out.json> [round]"""
import json, sys

ROOT = sys.argv[1]
OUT = sys.argv[2]
ROUND = sys.argv[3] if len(sys.argv) > 3 else '1'
L = 'tests/SuffixTree/SuffixTree.Esa.Tests/Large/'
F_CF = 'FullyQualifiedName~Large.ClosedFormSelfTests'
F_CRC = 'FullyQualifiedName~Large.Crc32CCombineTests'
F_SP = 'FullyQualifiedName~Large.SparseContainerTests'
F_ENC = 'FullyQualifiedName~Large.Sesa5EncoderTests'
STAGE1 = {'PowerOfA.cs': F_CF, 'PeriodicAkB.cs': F_CF, 'DeBruijnRepeat.cs': F_CF, 'Crc32CCombine.cs': F_CRC, 'SparseContainer.cs': F_SP, 'Sesa5Encoder.cs': F_ENC, 'Sesa5Image.cs': F_ENC}

M = []


def m(id_, desc, file, old, new):
    M.append({'id': 's21_' + id_, 'desc': desc, 'file': L + file, 'old': old, 'new': new, 'filter': '(' + STAGE1[file] + ')'})


# ---- a^n
m('P01', 'a^n SuffixAt: off by one', 'PowerOfA.cs', 'return Length - rank;', 'return Length - rank - 1;')
m('P02', 'a^n RankOf: off by one', 'PowerOfA.cs', 'return Length - position;', 'return Length - position + 1;')
m('P03', 'a^n Lcp: r instead of r - 1', 'PowerOfA.cs', 'return rank <= 1 ? 0 : rank - 1;', 'return rank <= 1 ? 0 : rank;')
m('P04', 'a^n LcpSum: n narrowed to int', 'PowerOfA.cs', '(Int128)Length * (Length - 1) / 2', '(Int128)(int)Length * (Length - 1) / 2')
m('P05', 'a^n InternalCount: n - 2', 'PowerOfA.cs', 'Length == 0 ? 0 : Length - 1;', 'Length == 0 ? 0 : Length - 2;')
m('P06', 'a^n LrsRank: n - 2', 'PowerOfA.cs', 'internal override long LrsRank => Length >= 2 ? Length - 1 : -1;', 'internal override long LrsRank => Length >= 2 ? Length - 2 : -1;')
m('P07', 'a^n RangeOfRun: first rank L + 1', 'PowerOfA.cs', 'new RankRange(length, Length);', 'new RankRange(length + 1, Length);')
m('P08', 'a^n RangeOfRun: a run of n symbols is absent', 'PowerOfA.cs', 'return length > Length ? RankRange.Empty', 'return length >= Length ? RankRange.Empty')
m('P09', 'a^n Down: one past the end', 'PowerOfA.cs', 'return i < Length ? i + 1 : 0;', 'return i <= Length ? i + 1 : 0;')
m('P10', 'a^n SuffixAt: the rank narrowed to int', 'PowerOfA.cs', 'return Length - rank;', 'return Length - (long)(int)rank;')
m('P11', 'a^n SymbolCount: n - 1', 'PowerOfA.cs', 'code == 0 ? Length : 0;', 'code == 0 ? Length - 1 : 0;')

# ---- (a^k b)*
m('K01', '(a^k b)* SuffixAt: the position of the j-th element of a block is off by one period', 'PeriodicAkB.cs', 'return ((_q - 1 - j) * _p) + phase;', 'return ((_q - j) * _p) + phase;')
m('K02', '(a^k b)* SuffixAt: the rank narrowed to int', 'PeriodicAkB.cs', 'if (rank <= _t) return Length - rank; // the pure runs (all ranks when q = 0)\n        long x = rank - 1 - _t;', 'if (rank <= _t) return Length - rank; // the pure runs (all ranks when q = 0)\n        long x = (int)(rank - 1 - _t);')
m('K03', '(a^k b)* RankOf: the partial period starts one late', 'PeriodicAkB.cs', 'if (position >= _q * _p) return Length - position;', 'if (position > _q * _p) return Length - position;')
m('K04', '(a^k b)* RankOf: the index in the block is off by one', 'PeriodicAkB.cs', 'return Rank(phase, _q - 1 - c);', 'return Rank(phase, _q - c);')
m('K05', '(a^k b)* Lcp: the first element of block 0 has LCP t + 1', 'PeriodicAkB.cs', 'if (rank == _t + 1) return _t;', 'if (rank == _t + 1) return _t + 1;')
m('K06', '(a^k b)* Lcp: t + j p - phase + 1', 'PeriodicAkB.cs', 'return j == 0 ? _k - phase : _t + (j * _p) - phase;', 'return j == 0 ? _k - phase : _t + (j * _p) - phase + 1;')
m('K07', '(a^k b)* Lcp: j narrowed to int', 'PeriodicAkB.cs', 'return j == 0 ? _k - phase : _t + (j * _p) - phase;', 'return j == 0 ? _k - phase : _t + ((int)j * _p) - phase;')
m('K08', '(a^k b)* InternalCount: max(t, k) for max(t, k - 1)', 'PeriodicAkB.cs', 'Math.Max(_t, _k - 1);\n\n    internal override long LrsDepth', 'Math.Max(_t, _k);\n\n    internal override long LrsDepth')
m('K09', '(a^k b)* LrsRank: t + q', 'PeriodicAkB.cs', 'if (_q >= 2) return _t + _q - 1;', 'if (_q >= 2) return _t + _q;')
m('K10', '(a^k b)* RangeOfRun: the last rank misses one phase', 'PeriodicAkB.cs', 'long last = _t + ((_k - length + 1) * _q);', 'long last = _t + ((_k - length) * _q);')
m('K11', '(a^k b)* RangeOfRun: the pure run a^t is not a start of the range', 'PeriodicAkB.cs', 'long first = length <= _t ? length : _t + 1;', 'long first = length < _t ? length : _t + 1;')
m('K12', '(a^k b)* Occurrences: one period too many', 'PeriodicAkB.cs', 'long cMax = room / _p;', 'long cMax = (room + _p - 1) / _p;')
m('K13', '(a^k b)* Down: the chain extends one element too far', 'PeriodicAkB.cs', 'return j <= _q - 2 ? i + 1 : 0;', 'return j <= _q - 1 ? i + 1 : 0;')
m('K14', '(a^k b)* Down: the chain of block 0 is a last child also when t < k', 'PeriodicAkB.cs', '(phase >= 1 || _t == _k) ? i + 1 : 0;', '(phase >= 0 || _t == _k) ? i + 1 : 0;')
m('K15', '(a^k b)* Up: the chain (0, 1) ends one late', 'PeriodicAkB.cs', 'j == _t + _q ? _t + 2 : 0;', 'j == _t + _q + 1 ? _t + 2 : 0;')
m('K16', '(a^k b)* Up: an interval of value k beyond max(t, k - 1)', 'PeriodicAkB.cs', 'if (v >= 1 && v <= Math.Max(_t, _k - 1)) return FirstOfValue(v);', 'if (v >= 1 && v <= Math.Max(_t, _k)) return FirstOfValue(v);')
m('K17', '(a^k b)* Next: the second l-index of the root is Rank(k, 1)', 'PeriodicAkB.cs', 'if (k == 1) return Rank(_k, 0);', 'if (k == 1) return Rank(_k, 1);')
m('K18', '(a^k b)* Next: intervals of value up to min(t, k) have two l-indices', 'PeriodicAkB.cs', 'k - 1 <= Math.Min(_t, _k - 1) ?', 'k - 1 <= Math.Min(_t, _k) ?')
m('K19', '(a^k b)* LcpSum: (q - 1) k t becomes q k t', 'PeriodicAkB.cs', '((q - 1) * k * t)', '(q * k * t)')
m('K20', '(a^k b)* InfiniteCodeAt: the b is one symbol early', 'PeriodicAkB.cs', 'position % _p == _k ? 1 : 0;', 'position % _p == _k - 1 ? 1 : 0;')
m('K21', '(a^k b)* FirstOfValue: v < t for v <= t', 'PeriodicAkB.cs', 'v <= _t ? v + 1 : 1 + _t + ((_k - v) * _q);', 'v < _t ? v + 1 : 1 + _t + ((_k - v) * _q);')
m('K22', '(a^k b)* RangeOfFactor: the start narrowed to int', 'PeriodicAkB.cs', 'long phase = start % _p;', 'long phase = (int)start % _p;')

# ---- u^m, de Bruijn
m('D01', 'u^m: the rank of a short suffix lacks the empty suffix', 'DeBruijnRepeat.cs', '_specialRank[l - 1] = 1 + smallerLong + smallerShort;', '_specialRank[l - 1] = smallerLong + smallerShort;')
m('D02', 'u^m: the windows of the short suffixes are not taken out of the count', 'DeBruijnRepeat.cs', 'long smallerLong = (_m * below) - Enumerable.Range(1, specials).Count(other => _windowOfRotation[_period - other] < below);', 'long smallerLong = _m * below;')
m('D03', 'u^m: the window LCP stops one digit early', 'DeBruijnRepeat.cs', 'while (common < _order && (a / divisor)', 'while (common < _order - 1 && (a / divisor)')
m('D04', 'u^m: InternalCount + 1', 'DeBruijnRepeat.cs', 'return words + (Length - _period - _order + 1);', 'return words + (Length - _period - _order + 2);')
m('D05', 'u^m: RangeOf forgets the occurrences that wrap around the period', 'DeBruijnRepeat.cs', 'long count = (_m * span) - wrapped;', 'long count = _m * span;')
m('D06', 'u^m: InfiniteCodeAt narrows the position', 'DeBruijnRepeat.cs', '_u[position % _period];', '_u[(int)position % _period];')
m('D07', 'u^m: RankOf, the index in a block is off by one', 'DeBruijnRepeat.cs', 'long pure = (_m * _windowOfRotation[r]) + (_m - 1 - c);', 'long pure = (_m * _windowOfRotation[r]) + (_m - c);')
m('D08', 'u^m: SuffixAt skips a special index one late', 'DeBruijnRepeat.cs', 'foreach (long q in _sortedPure)\n            if (q <= pure) pure++;', 'foreach (long q in _sortedPure)\n            if (q < pure) pure++;')
m('D09', 'u^m: InRotation counts one position more', 'DeBruijnRepeat.cs', 'long cMax = room / _period;', 'long cMax = (room / _period) + 1;')
m('D10', 'u^m: LrsRank + 1', 'DeBruijnRepeat.cs', 'internal override long LrsRank => RankOf(_period);', 'internal override long LrsRank => RankOf(_period) + 1;')
m('D11', 'u^m: the short suffixes are those of length at most d', 'DeBruijnRepeat.cs', 'position < Length && Length - position < _order;', 'position < Length && Length - position <= _order;')

# ---- CRC combine
m('C01', 'CRC: one bit of the polynomial', 'Crc32CCombine.cs', 'private const uint ReflectedPolynomial = 0x82F63B78;', 'private const uint ReflectedPolynomial = 0x82F63B79;')
m('C02', 'CRC: Multiply stops after 31 coefficients', 'Crc32CCombine.cs', 'for (int i = 0; i < 32; i++)', 'for (int i = 0; i < 31; i++)')
m('C03', 'CRC: Multiply adds with OR, not XOR', 'Crc32CCombine.cs', 'if ((a & m) != 0) product ^= b;', 'if ((a & m) != 0) product |= b;')
m('C04', 'CRC: the base of the power is x^2', 'Crc32CCombine.cs', 'uint square = One >> 1; // x', 'uint square = One >> 2; // x')
m('C05', 'CRC: 7 bits per byte', 'Crc32CCombine.cs', 'checked(8UL * (ulong)zeroBytes)', 'checked(7UL * (ulong)zeroBytes)')
m('C06', 'CRC: OfZeros misses the final inversion', 'Crc32CCombine.cs', 'Shift(0xFFFFFFFF, length) ^ 0xFFFFFFFF;', 'Shift(0xFFFFFFFF, length);')
m('C07', 'CRC: wrong carry, the zero run before a patch is one byte short', 'Crc32CCombine.cs', 'if (offset > at) crc = Combine(crc, OfZeros(offset - at), offset - at);', 'if (offset > at) crc = Combine(crc, OfZeros(offset - at), offset - at - 1);')
m('C08', 'CRC: the length of a patch is one too many', 'Crc32CCombine.cs', 'crc = Combine(crc, Crc32C.Compute(data), data.Length);', 'crc = Combine(crc, Crc32C.Compute(data), data.Length + 1);')
m('C09', 'CRC: a tail of one zero byte is not combined', 'Crc32CCombine.cs', 'if (length > at) crc = Combine(crc, OfZeros(length - at), length - at);', 'if (length > at + 1) crc = Combine(crc, OfZeros(length - at), length - at);')
m('C10', 'CRC: the number of zero bytes narrowed to int', 'Crc32CCombine.cs', 'checked(8UL * (ulong)zeroBytes)', 'checked(8UL * (ulong)(int)zeroBytes)')
m('C11', 'CRC: Combine joins the wrong pieces', 'Crc32CCombine.cs', 'Shift(crcOfFirst, lengthOfSecond) ^ crcOfSecond;', 'Shift(crcOfSecond, lengthOfSecond) ^ crcOfFirst;')
m('C12', 'CRC: Multiply reduces nothing', 'Crc32CCombine.cs', 'b = (b & 1) != 0 ? (b >> 1) ^ ReflectedPolynomial : b >> 1;', 'b >>= 1;')

# ---- the sparse container writer
m('S01', 'sparse: a position section has 2 bytes of padding, not 3', 'SparseContainer.cs', '(width == 5 ? ZeroPadding5 : 0);', '(width == 5 ? ZeroPadding5 - 1 : 0);')
m('S02', 'sparse: TEXT has 63 bytes of padding', 'SparseContainer.cs', '/ 8 + TextPadding;', '/ 8 + TextPadding - 1;')
m('S03', 'sparse: a wrong offset in the section table (CHILD overlaps SA)', 'SparseContainer.cs', 'long childOffset = SesaFormat.AlignedOffset(saOffset + arrayLength, arrayLength);', 'long childOffset = SesaFormat.AlignedOffset(saOffset, arrayLength);')
m('S04', 'sparse: width 5 from 2^32 - 4', 'SparseContainer.cs', 'n <= MaxContainer4 ? 4 : 5;', 'n < MaxContainer4 ? 4 : 5;')
m('S05', 'sparse: U8 up to 255', 'SparseContainer.cs', 'maxLcp <= 254 ? 1', 'maxLcp <= 255 ? 1')
m('S06', 'sparse: a width-5 container says version 1.0', 'SparseContainer.cs', 'BinaryPrimitives.WriteUInt16LittleEndian(page.AsSpan(10), 1);', 'BinaryPrimitives.WriteUInt16LittleEndian(page.AsSpan(10), 0);')
m('S07', 'sparse: the header CRC is not recomputed after the version', 'SparseContainer.cs', 'BinaryPrimitives.WriteUInt16LittleEndian(page.AsSpan(10), 1);\n            ResealHeader(page);', 'BinaryPrimitives.WriteUInt16LittleEndian(page.AsSpan(10), 1);')
m('S08', 'sparse: big-endian entries', 'SparseContainer.cs', 'target[offset + i] = (byte)(value >> (8 * i));', 'target[offset + i] = (byte)(value >> (8 * (width - 1 - i)));')
m('S09', 'sparse: the section CRC ignores the patches', 'SparseContainer.cs', 'Crc = Crc32CCombine.OfSparse(s.Length, parts) }', 'Crc = Crc32CCombine.OfZeros(s.Length) }')
m('S10', 'sparse: the tail of the entries starts one rank late', 'SparseContainer.cs', 'long tailFrom = Math.Max(head, rows - entriesPerEnd);', 'long tailFrom = Math.Max(head, rows - entriesPerEnd + 1);')
m('S11', 'sparse: the CHILD cell is the down field only', 'SparseContainer.cs', 'PutPosition(child, at, w, f.PackedCell(r));', 'PutPosition(child, at, w, f.Down(r));')
m('S12', 'sparse: the dense limit is 128 MiB', 'SparseContainer.cs', 'if (dense > LargeDisk.MaxDenseBytes) throw', 'if (dense > 2 * LargeDisk.MaxDenseBytes) throw')
m('S13', 'sparse: a position of exactly 2^40 is accepted at width 5', 'SparseContainer.cs', '(width == 5 && value >= 1L << 40)', '(width == 5 && value > 1L << 40)')
m('S14', 'sparse: the text bits are least significant first', 'SparseContainer.cs', 'data[(bit + b) >> 3] |= (byte)(0x80 >> (int)((bit + b) & 7));', 'data[(bit + b) >> 3] |= (byte)(0x01 << (int)((bit + b) & 7));')
m('S15', 'sparse: PSI of the last suffix', 'SparseContainer.cs', 'return position + 1 >= f.Length ? 0 : f.RankOf(position + 1);', 'return position + 1 > f.Length ? 0 : f.RankOf(position + 1);')
m('S16', 'sparse: the LCP is aligned for the narrowest codec of the width', 'SparseContainer.cs', 'int widestLcp = width == 5 ? 8 : 4;', 'int widestLcp = width == 5 ? 4 : 2;')

# ---- the re-encoder and its reader
m('E01', 'encoder: 2 bytes of padding', 'Sesa5Encoder.cs', 'private const int Padding = 3;', 'private const int Padding = 2;')
m('E02', 'encoder: swapped endianness', 'Sesa5Encoder.cs', 'bytes[(5 * i) + j] = (byte)(values[i] >> (8 * j));', 'bytes[(5 * i) + j] = (byte)(values[i] >> (8 * (4 - j)));')
m('E03', 'encoder: the fifth byte is not stored', 'Sesa5Encoder.cs', 'for (int j = 0; j < 5; j++) bytes[(5 * i) + j]', 'for (int j = 0; j < 4; j++) bytes[(5 * i) + j]')
m('E04', 'encoder: a next link to k is dropped, not the one to k + 1', 'Sesa5Encoder.cs', '(kind == 1 && value == (uint)k + 1)', '(kind == 1 && value == (uint)k)')
m('E05', 'encoder: a next link to k + 1 is kept', 'Sesa5Encoder.cs', 'packed[k] = kind == 0 || (kind == 1 && value == (uint)k + 1) ? 0 : value;', 'packed[k] = kind == 0 ? 0 : value;')
m('E06', 'encoder: the version stays 1.0', 'Sesa5Encoder.cs', 'BinaryPrimitives.WriteUInt16LittleEndian(page.AsSpan(10), 1);', 'BinaryPrimitives.WriteUInt16LittleEndian(page.AsSpan(10), 0);')
m('E07', 'encoder: the header CRC is not recomputed after the version', 'Sesa5Encoder.cs', 'BinaryPrimitives.WriteUInt16LittleEndian(page.AsSpan(10), 1);\n        SparseContainer.ResealHeader(page);', 'BinaryPrimitives.WriteUInt16LittleEndian(page.AsSpan(10), 1);')
m('E08', 'encoder: a section is not aligned', 'Sesa5Encoder.cs', 'long at = SesaFormat.AlignedOffset(offset, payload.Length);', 'long at = offset;')
m('E09', 'encoder: the CRC of the 1.0 section is kept', 'Sesa5Encoder.cs', 's.HasCrc ? Crc32C.Compute(payload) : 0', 's.HasCrc ? s.Crc : 0')
m('E10', 'encoder: POS_WIDTH stays 4', 'Sesa5Encoder.cs', 'PosWidth = 5,', 'PosWidth = 4,')
m('E11', 'encoder: FILE_SIZE is 64 too large', 'Sesa5Encoder.cs', 'FileSize = offset,', 'FileSize = offset + 64,')
m('E12', 'encoder: a pair of the overflow list has 8 bytes', 'Sesa5Encoder.cs', 'return (PackPositions(pairs), 10, s.Codec);', 'return (PackPositions(pairs), 8, s.Codec);')
m('E13', 'encoder: ELEM_WIDTH 4 in the table', 'Sesa5Encoder.cs', 'return (PackPositions(cells), 5, codec);', 'return (PackPositions(cells), 4, codec);')
m('E14', 'encoder: BUILD_INFO still says tagged', 'Sesa5Encoder.cs', '.Replace("child=layout A, tagged", "child=layout A, packed", StringComparison.Ordinal);', ';')
m('E15', 'encoder: a tagged CHILD keeps its codec id', 'Sesa5Encoder.cs', 'codec = SesaFormat.ChildCodecPacked;', 'codec = s.Codec;')
m('E16', 'encoder: LRS_RANK is lost', 'Sesa5Encoder.cs', 'LrsRank = o.LrsRank,', 'LrsRank = -1,')
m('E17', 'encoder: TEXT_SHA256 is lost', 'Sesa5Encoder.cs', 'TextSha256 = o.TextSha256,', 'TextSha256 = new byte[32],')
m('E18', 'encoder: a container of width 5 is re-encoded', 'Sesa5Encoder.cs', 'if (o.PosWidth != 4) throw', 'if (o.PosWidth != 4 && false) throw')
m('E19', 'encoder: a damaged header CRC is not noticed', 'Sesa5Encoder.cs', 'if (o.HeaderCrc != o.ExpectedHeaderCrc()) throw', 'if (o.HeaderCrc != o.ExpectedHeaderCrc() && false) throw')
m('E20', 'encoder: a wrong FILE_SIZE is not noticed', 'Sesa5Encoder.cs', 'if (o.FileSize != container.Length) throw', 'if (o.FileSize != container.Length && false) throw')
m('E21', 'encoder: the primary row of PSI is lost', 'Sesa5Encoder.cs', 's.Count, elemWidth, codec, s.Param0, s.Param1,', 's.Count, elemWidth, codec, 0, s.Param1,')
m('E22', 'encoder: the sections are laid out in table order, not file order', 'Sesa5Encoder.cs', 'foreach (OracleSection s in o.Sections.OrderBy(x => x.Offset))', 'foreach (OracleSection s in o.Sections.OrderByDescending(x => x.Offset))')
m('I01', 'reader: a 4-byte mask', 'Sesa5Image.cs', 'return BinaryPrimitives.ReadUInt64LittleEndian(load) & Mask40;', 'return BinaryPrimitives.ReadUInt64LittleEndian(load) & 0xFFFF_FFFFUL;')
m('I02', 'reader: the value of a pair at +4', 'Sesa5Image.cs', '(10 * i) + 5), 8)) & Mask40;', '(10 * i) + 4), 8)) & Mask40;')
m('I03', 'reader: the padding is 2 bytes', 'Sesa5Image.cs', 'file.AsSpan((int)(s.Offset + (width * s.Count)), 3)', 'file.AsSpan((int)(s.Offset + (width * s.Count)), 2)')

if ROUND != '1':
    exec(open(sys.argv[4]).read())  # a second round: extra mutants defined in the file named by the fourth argument

for x in M:
    text = open(ROOT + '/' + x['file']).read()
    n = text.count(x['old'])
    if n != 1:
        raise SystemExit('%s: the old text occurs %d times in %s' % (x['id'], n, x['file']))
    if x['old'] == x['new']:
        raise SystemExit('%s: no change' % x['id'])
json.dump(M, open(OUT, 'w'), indent=1)
print(len(M), 'mutants written to', OUT)
