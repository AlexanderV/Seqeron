#!/usr/bin/env python3
"""The hand-made mutants of work package S2.3 (prefix s23_): single-site changes of the position readers, the views over them (SaView, CompactLcp + LcpStore, MappedChild + ChildCells, ChildTable, LcpPyramid, EsaText), the checked
   bridges (IntBridge) and the container reader's views (SesaFile). Every `long` to `int` narrowing of the new code, every constant of the 40-bit decode (mask, stride, padding), the overflow-pair width, the U64 branch, the bound of
   the Checked mode, the sign extension and an off-by-one at a section end is a mutant. Writes s23_listed.json (the list that run_s23_mutants.py runs) and checks that every `old` text occurs exactly once in its file.
   usage: make_s23_mutants.py <repo root> <out.json>"""
import json, sys

ROOT = sys.argv[1]
OUT = sys.argv[2]
S = 'src/SuffixTree/Algorithms/SuffixTree.Esa/'
V = 'FullyQualifiedName~SuffixTree.Esa.Tests.Views.'
F_READER = '(' + V + 'PosReaderTests)'
F_SA = '(' + V + 'SaViewTests)'
F_LCP = '(' + V + 'CompactLcpWideTests)'
F_CHILD = '(' + V + 'ChildCellsTests)'
F_TEXT = '(' + V + 'TextViewTests)'
F_PYR = '(' + V + 'LcpPyramidWideTests|FullyQualifiedName~SuffixTree.Esa.Tests.Navigation.LcpPyramidTests)'
F_BRIDGE = '(' + V + 'IntBridgeTests)'
F_SF = '(FullyQualifiedName~SuffixTree.Esa.Tests.External.SesaReader|FullyQualifiedName~SuffixTree.Esa.Tests.External.SesaFailure)'
FILTER = {'PosReaders.cs': F_READER, 'SaView.cs': F_SA, 'CompactLcp.cs': F_LCP, 'LcpStore.cs': F_LCP, 'ChildLayouts.cs': F_CHILD, 'ChildTable.cs': F_CHILD, 'LcpPyramid.cs': F_PYR,
          'EsaText.cs': F_TEXT, 'IntBridge.cs': F_BRIDGE, 'External/Sesa/SesaFile.cs': F_SF, 'External/Sesa/SesaVerifier.cs': F_SF}
M = []


def m(id_, desc, file, old, new):
    M.append({'id': 's23_' + id_, 'desc': desc, 'file': S + file, 'old': old, 'new': new, 'filter': FILTER[file]})


P = 'PosReaders.cs'
# ---- the readers
m('R01', 'Pos8: the bound is element > count', P, 'if (TAccess.IsChecked && (ulong)index >= (ulong)_count) PosReader.ThrowOutside(index, _count);\n        return _p[index];', 'if (TAccess.IsChecked && (ulong)index > (ulong)_count) PosReader.ThrowOutside(index, _count);\n        return _p[index];')
m('R02', 'Pos16: the bound is element > count', P, 'if (TAccess.IsChecked && (ulong)index >= (ulong)_count) PosReader.ThrowOutside(index, _count);\n        return Unsafe.ReadUnaligned<ushort>(', 'if (TAccess.IsChecked && (ulong)index > (ulong)_count) PosReader.ThrowOutside(index, _count);\n        return Unsafe.ReadUnaligned<ushort>(')
m('R03', 'Pos32: the bound is element > count', P, 'if (TAccess.IsChecked && (ulong)index >= (ulong)_count) PosReader.ThrowOutside(index, _count);\n        return Unsafe.ReadUnaligned<uint>(', 'if (TAccess.IsChecked && (ulong)index > (ulong)_count) PosReader.ThrowOutside(index, _count);\n        return Unsafe.ReadUnaligned<uint>(')
m('R04', 'Pos40: the bound is element > count', P, 'if (TAccess.IsChecked && (ulong)index >= (ulong)_count) PosReader.ThrowOutside(index, _count);\n        return Unsafe.ReadUnaligned<ulong>(_p + (5 * index))', 'if (TAccess.IsChecked && (ulong)index > (ulong)_count) PosReader.ThrowOutside(index, _count);\n        return Unsafe.ReadUnaligned<ulong>(_p + (5 * index))')
m('R05', 'Pos64: the bound is element > count', P, 'if (TAccess.IsChecked && (ulong)index >= (ulong)_count) PosReader.ThrowOutside(index, _count);\n        return Unsafe.ReadUnaligned<ulong>(_p + (8 * index))', 'if (TAccess.IsChecked && (ulong)index > (ulong)_count) PosReader.ThrowOutside(index, _count);\n        return Unsafe.ReadUnaligned<ulong>(_p + (8 * index))')
m('R06', 'Pos8: Checked mode does not check', P, 'if (TAccess.IsChecked && (ulong)index >= (ulong)_count) PosReader.ThrowOutside(index, _count);\n        return _p[index];', 'if (false && (ulong)index >= (ulong)_count) PosReader.ThrowOutside(index, _count);\n        return _p[index];')
m('R07', 'Pos16: Checked mode does not check', P, 'if (TAccess.IsChecked && (ulong)index >= (ulong)_count) PosReader.ThrowOutside(index, _count);\n        return Unsafe.ReadUnaligned<ushort>(', 'if (false && (ulong)index >= (ulong)_count) PosReader.ThrowOutside(index, _count);\n        return Unsafe.ReadUnaligned<ushort>(')
m('R08', 'Pos32: Checked mode does not check', P, 'if (TAccess.IsChecked && (ulong)index >= (ulong)_count) PosReader.ThrowOutside(index, _count);\n        return Unsafe.ReadUnaligned<uint>(', 'if (false && (ulong)index >= (ulong)_count) PosReader.ThrowOutside(index, _count);\n        return Unsafe.ReadUnaligned<uint>(')
m('R09', 'Pos40: Checked mode does not check', P, 'if (TAccess.IsChecked && (ulong)index >= (ulong)_count) PosReader.ThrowOutside(index, _count);\n        return Unsafe.ReadUnaligned<ulong>(_p + (5 * index))', 'if (false && (ulong)index >= (ulong)_count) PosReader.ThrowOutside(index, _count);\n        return Unsafe.ReadUnaligned<ulong>(_p + (5 * index))')
m('R10', 'Pos64: Checked mode does not check', P, 'if (TAccess.IsChecked && (ulong)index >= (ulong)_count) PosReader.ThrowOutside(index, _count);\n        return Unsafe.ReadUnaligned<ulong>(_p + (8 * index))', 'if (false && (ulong)index >= (ulong)_count) PosReader.ThrowOutside(index, _count);\n        return Unsafe.ReadUnaligned<ulong>(_p + (8 * index))')
m('R11', 'Pos40: the mask keeps 32 bits', P, 'return Unsafe.ReadUnaligned<ulong>(_p + (5 * index)) & PosReader.Mask40;', 'return Unsafe.ReadUnaligned<ulong>(_p + (5 * index)) & 0xFFFF_FFFFUL;')
m('R12', 'Pos40: the mask keeps 48 bits', P, 'return Unsafe.ReadUnaligned<ulong>(_p + (5 * index)) & PosReader.Mask40;', 'return Unsafe.ReadUnaligned<ulong>(_p + (5 * index)) & 0xFFFF_FFFF_FFFFUL;')
m('R13', 'Pos40: the stride is 4', P, 'return Unsafe.ReadUnaligned<ulong>(_p + (5 * index)) & PosReader.Mask40;', 'return Unsafe.ReadUnaligned<ulong>(_p + (4 * index)) & PosReader.Mask40;')
m('R14', 'Pos40: the index is narrowed to int in the address', P, 'return Unsafe.ReadUnaligned<ulong>(_p + (5 * index)) & PosReader.Mask40;', 'return Unsafe.ReadUnaligned<ulong>(_p + (5 * (int)index)) & PosReader.Mask40;')
m('R15', 'Pos32: the index is narrowed to int in the address', P, 'return Unsafe.ReadUnaligned<uint>(_p + (4 * index));', 'return Unsafe.ReadUnaligned<uint>(_p + (4 * (int)index));')
m('R16', 'Pos64: the index is narrowed to int in the address', P, 'return Unsafe.ReadUnaligned<ulong>(_p + (8 * index));', 'return Unsafe.ReadUnaligned<ulong>(_p + (8 * (int)index));')
m('R17', 'Pos16: the index is narrowed to int in the address', P, 'return Unsafe.ReadUnaligned<ushort>(_p + (2 * index));', 'return Unsafe.ReadUnaligned<ushort>(_p + (2 * (int)index));')
m('R18', 'Pos8: the index is narrowed to int', P, 'return _p[index];', 'return _p[(int)index];')
m('R19', 'Pos32: the word is sign-extended', P, 'return Unsafe.ReadUnaligned<uint>(_p + (4 * index));', 'return (ulong)Unsafe.ReadUnaligned<int>(_p + (4 * index));')
m('R20', 'Pos16: the word is sign-extended', P, 'return Unsafe.ReadUnaligned<ushort>(_p + (2 * index));', 'return (ulong)Unsafe.ReadUnaligned<short>(_p + (2 * index));')
m('R21', 'PosInt: the entry is sign-extended', P, 'public ulong Get(long index) => (uint)_values[index];', 'public ulong Get(long index) => (ulong)_values[index];')
m('R22', 'Mask40 loses its lowest bit', P, 'internal const ulong Mask40 = 0xFF_FFFF_FFFFUL;', 'internal const ulong Mask40 = 0xFF_FFFF_FFFEUL;')
m('R23', 'the padding is 4 bytes', P, 'internal const int Padding40 = 3;', 'internal const int Padding40 = 4;')
m('R24', 'Pos40 says 4 bytes per element', P, 'public static int ElementBytes => 5;', 'public static int ElementBytes => 4;')
m('R25', 'a negative count is accepted', P, 'internal static void CheckCount(long count) => ArgumentOutOfRangeException.ThrowIfNegative(count);', 'internal static void CheckCount(long count) { }')
m('R26', 'PosInt counts one element less', P, 'public long Count => _values.Length;', 'public long Count => _values.Length - 1;')
m('R27', 'Pos40 counts one element less', P, 'public static int ElementBytes => 5;\n\n    public long Count => _count;', 'public static int ElementBytes => 5;\n\n    public long Count => _count - 1;')
m('R28', 'Pos64 says 4 bytes per element', P, 'public static int ElementBytes => 8;', 'public static int ElementBytes => 4;')

SV = 'SaView.cs'
m('S01', 'SaView.FromMapped accepts a width of 8', SV, 'if (posWidth is not (4 or 5)) throw new ArgumentOutOfRangeException(nameof(posWidth), posWidth, "the stored width of a position is 4 or 5 bytes");', 'if (posWidth is not (4 or 5 or 8)) throw new ArgumentOutOfRangeException(nameof(posWidth), posWidth, "the stored width of a position is 4 or 5 bytes");')
m('S02', 'SaView.FromMapped sets the 5-byte flag for width 4', SV, 'if (posWidth == 5) flags |= FiveBytes;', 'if (posWidth == 4) flags |= FiveBytes;')
m('S03', 'SaView.Load40: the mask keeps 32 bits', SV, 'Unsafe.ReadUnaligned<ulong>((byte*)_pointer + (5 * rank)) & PosReader.Mask40;', 'Unsafe.ReadUnaligned<ulong>((byte*)_pointer + (5 * rank)) & 0xFFFF_FFFFUL;')
m('S04', 'SaView.Load40: the rank is narrowed to int', SV, 'Unsafe.ReadUnaligned<ulong>((byte*)_pointer + (5 * rank)) & PosReader.Mask40;', 'Unsafe.ReadUnaligned<ulong>((byte*)_pointer + (5 * (int)rank)) & PosReader.Mask40;')
m('S05', 'SaView.Load40: the stride is 4', SV, 'Unsafe.ReadUnaligned<ulong>((byte*)_pointer + (5 * rank)) & PosReader.Mask40;', 'Unsafe.ReadUnaligned<ulong>((byte*)_pointer + (4 * rank)) & PosReader.Mask40;')
m('S06', 'SaView[long]: the bound is rank > count', SV, 'if ((_flags & Checked) != 0 && (ulong)rank >= (ulong)_count) ThrowOutside(rank);\n            return (_flags & FiveBytes) != 0 ? (long)Load40(rank) : (uint)_pointer[rank];', 'if ((_flags & Checked) != 0 && (ulong)rank > (ulong)_count) ThrowOutside(rank);\n            return (_flags & FiveBytes) != 0 ? (long)Load40(rank) : (uint)_pointer[rank];')
m('S07', 'SaView[long]: a 4-byte entry is sign-extended', SV, 'return (_flags & FiveBytes) != 0 ? (long)Load40(rank) : (uint)_pointer[rank];\n        }\n    }\n\n    /// <summary>The same read, always bounded', 'return (_flags & FiveBytes) != 0 ? (long)Load40(rank) : _pointer[rank];\n        }\n    }\n\n    /// <summary>The same read, always bounded')
m('S08', 'SaView[int]: the bound is rank > count', SV, 'if ((_flags & Checked) != 0 && (ulong)rank >= (ulong)_count) ThrowOutside(rank);\n            if ((_flags & Bridged) != 0) return IntBridge.Value(', 'if ((_flags & Checked) != 0 && (ulong)rank > (ulong)_count) ThrowOutside(rank);\n            if ((_flags & Bridged) != 0) return IntBridge.Value(')
m('S09', 'SaView[int]: a 5-byte entry is narrowed without the check', SV, 'return IntBridge.Value(ReadUnchecked(rank), "a suffix array entry");\n            return _pointer[rank];', 'return (int)ReadUnchecked(rank);\n            return _pointer[rank];')
m('S10', 'SaView.At: the bound is rank > count', SV, 'if ((ulong)rank >= (ulong)_count) ThrowOutside(rank);\n        return this[rank];', 'if ((ulong)rank > (ulong)_count) ThrowOutside(rank);\n        return this[rank];')
m('S11', 'SaView.CopyTo(int): the range may end one past the section', SV, 'if (start < 0 || destination.Length > _count - start) ThrowOutside(start);\n        if (_array is not null) _array.AsSpan((int)start', 'if (start < 0 || destination.Length > _count - start + 1) ThrowOutside(start);\n        if (_array is not null) _array.AsSpan((int)start')
m('S12', 'SaView.CopyTo(long): the range may end one past the section', SV, 'if (start < 0 || destination.Length > _count - start) ThrowOutside(start);\n        for (int i = 0; i < destination.Length; i++) destination[i] = ReadUnchecked(start + i);', 'if (start < 0 || destination.Length > _count - start + 1) ThrowOutside(start);\n        for (int i = 0; i < destination.Length; i++) destination[i] = ReadUnchecked(start + i);')
m('S13', 'SaView.CopyTo(int): a 5-byte entry is narrowed without the check', SV, 'destination[i] = IntBridge.Value(ReadUnchecked(start + i), "a suffix array entry");', 'destination[i] = (int)ReadUnchecked(start + i);')
m('S14', 'SaView.CopyTo(int): the start is narrowed in the pointer', SV, 'new ReadOnlySpan<int>(_pointer + start, destination.Length)', 'new ReadOnlySpan<int>(_pointer + (int)start, destination.Length)')
m('S15', 'SaView.ReadUnchecked: the 5-byte test is inverted', SV, 'return (_flags & FiveBytes) != 0 ? (long)Load40(rank) : (uint)_pointer[rank];\n    }', 'return (_flags & FiveBytes) == 0 ? (long)Load40(rank) : (uint)_pointer[rank];\n    }')
m('S16', 'SaView.StoredWidth: 4 and 5 swapped', SV, 'internal int StoredWidth => (_flags & FiveBytes) != 0 ? 5 : 4;', 'internal int StoredWidth => (_flags & FiveBytes) != 0 ? 4 : 5;')
m('S17', 'SaView.Pointer: of a 5-byte section is returned', SV, 'internal int* Pointer => (_flags & FiveBytes) == 0 ? _pointer : throw', 'internal int* Pointer => (_flags & FiveBytes) != 0 ? _pointer : throw')
m('S18', 'SaView.ToArray: a negative count is not refused', SV, 'ArgumentOutOfRangeException.ThrowIfNegative(count);\n        var result = new int[count];', 'var result = new int[count];')
m('S19', 'SaView.ToLongArray: a negative count is not refused', SV, 'ArgumentOutOfRangeException.ThrowIfNegative(count);\n        var result = new long[count];\n        CopyTo(start, result);\n        return result;\n    }\n\n    /// <summary>\n    /// The pointer', 'var result = new long[count];\n        CopyTo(start, result);\n        return result;\n    }\n\n    /// <summary>\n    /// The pointer')
m('S20', 'SaView.FromMapped: a negative count is accepted', SV, 'ArgumentOutOfRangeException.ThrowIfNegative(count);\n        return new SaView((int*)section, count, FlagsOf(', 'return new SaView((int*)section, count, FlagsOf(')
m('S21', 'SaView<TPos>.At: the bound is rank > count', SV, 'if ((ulong)rank >= (ulong)_reader.Count) SaViewErrors.Outside(rank, _reader.Count);', 'if ((ulong)rank > (ulong)_reader.Count) SaViewErrors.Outside(rank, _reader.Count);')
m('S22', 'SaView<TPos>: the range may end one past the section', SV, 'if (start < 0 || length > _reader.Count - start) SaViewErrors.Outside(start, _reader.Count);', 'if (start < 0 || length > _reader.Count - start + 1) SaViewErrors.Outside(start, _reader.Count);')
m('S23', 'SaView<TPos>.CopyTo(int): the position is narrowed without the check', SV, 'destination[i] = IntBridge.Value((long)_reader.Get(start + i), "a suffix array entry");', 'destination[i] = (int)_reader.Get(start + i);')
m('S24', 'SaView<TPos>.CopyTo(long): the start is narrowed', SV, 'for (int i = 0; i < destination.Length; i++) destination[i] = (long)_reader.Get(start + i);', 'for (int i = 0; i < destination.Length; i++) destination[i] = (long)_reader.Get((int)start + i);')
m('S25', 'SaView<TPos>[long]: the entry is narrowed to int', SV, 'get => (long)_reader.Get(rank);', 'get => (int)_reader.Get(rank);')
m('S26', 'SaView<TPos>: a range starting at 0 is refused', SV, 'if (start < 0 || length > _reader.Count - start) SaViewErrors.Outside(start, _reader.Count);', 'if (start <= 0 || length > _reader.Count - start) SaViewErrors.Outside(start, _reader.Count);')
m('S27', 'SaView.FlagsOf: a 4-byte section of more than int.MaxValue entries is not marked', SV, 'else if (count > int.MaxValue) flags |= FourBytesWide;', 'else if (count > long.MaxValue) flags |= FourBytesWide;')
m('S28', 'SaView: the bridged forms are the 5-byte section only', SV, 'private const byte Bridged = FiveBytes | FourBytesWide;', 'private const byte Bridged = FiveBytes;')
m('S29', 'SaView.this[int]: the bridge is tested for the 5-byte flag only', SV, 'if ((_flags & Bridged) != 0) return IntBridge.Value(ReadUnchecked(rank), "a suffix array entry");', 'if ((_flags & FiveBytes) != 0) return IntBridge.Value(ReadUnchecked(rank), "a suffix array entry");')
m('S30', 'SaView.CopyTo(int): the bridge is tested for the 5-byte flag only', SV, 'else if ((_flags & Bridged) == 0) new ReadOnlySpan<int>(_pointer + start, destination.Length).CopyTo(destination);', 'else if ((_flags & FiveBytes) == 0) new ReadOnlySpan<int>(_pointer + start, destination.Length).CopyTo(destination);')
m('S31', 'SaView.CopyTo(long): the start is narrowed', SV, 'destination[i] = ReadUnchecked(start + i);', 'destination[i] = ReadUnchecked((int)start + i);')
m('S32', 'SaView.ReadUnchecked: the 4-byte index is narrowed', SV, ': (uint)_pointer[rank];\n    }\n\n    private ulong Load40', ': (uint)_pointer[(int)rank];\n    }\n\n    private ulong Load40')
m('S33', 'SaView[long]: the 4-byte index is narrowed', SV, '(uint)_pointer[rank];\n        }\n    }\n\n    /// <summary>The same read, always bounded in both modes', '(uint)_pointer[(int)rank];\n        }\n    }\n\n    /// <summary>The same read, always bounded in both modes')
m('S34', 'SaView.Load40: the stride is computed in int', SV, '(byte*)_pointer + (5 * rank)', '(byte*)_pointer + (5 * (int)rank)')

IB = 'IntBridge.cs'
m('B01', 'IntBridge.Length refuses int.MaxValue', IB, 'if ((ulong)length > int.MaxValue) ThrowLength(length, member);', 'if ((ulong)length >= int.MaxValue) ThrowLength(length, member);')
m('B02', 'IntBridge.Value refuses int.MaxValue', IB, 'if ((ulong)value > int.MaxValue) ThrowValue(value, what);', 'if ((ulong)value >= int.MaxValue) ThrowValue(value, what);')
m('B03', 'IntBridge.Length compares as uint', IB, 'if ((ulong)length > int.MaxValue) ThrowLength(length, member);', 'if ((uint)length > int.MaxValue) ThrowLength(length, member);')
m('B04', 'IntBridge.Value compares as uint', IB, 'if ((ulong)value > int.MaxValue) ThrowValue(value, what);', 'if ((uint)value > int.MaxValue) ThrowValue(value, what);')
m('B05', 'IntBridge.Length accepts a negative length', IB, 'if ((ulong)length > int.MaxValue) ThrowLength(length, member);', 'if (length > int.MaxValue) ThrowLength(length, member);')
m('B06', 'IntBridge.Value accepts a negative value', IB, 'if ((ulong)value > int.MaxValue) ThrowValue(value, what);', 'if (value > int.MaxValue) ThrowValue(value, what);')

CL = 'CompactLcp.cs'
m('C01', 'FromMappedWide: a U64 array may have an overflow list', CL, 'if (kind is LcpKind.U32 or LcpKind.U64 && pairCount != 0) throw new InvalidDataException("an LCP array of 32 or 64 bits per entry has no overflow list");\n        long fields', 'if (kind is LcpKind.U32 && pairCount != 0) throw new InvalidDataException("an LCP array of 32 or 64 bits per entry has no overflow list");\n        long fields')
m('C02', 'FromMappedWide: the pair fields of width 5 are read as 4 bytes', CL, 'posWidth == 5\n            ? FromReaders(kind, baseEntries, new Pos40<TrustedAccess>(pairs, fields), isChecked)', 'posWidth == 4\n            ? FromReaders(kind, baseEntries, new Pos40<TrustedAccess>(pairs, fields), isChecked)')
m('C03', 'FromMappedWide: one field per pair', CL, 'long fields = 2 * pairCount;', 'long fields = pairCount;')
m('C04', 'FromReaders: the width of the base reader is not checked', CL, 'if (TBase.ElementBytes != LcpStore.WidthOf(kind)) throw', 'if (false) throw')
m('C05', 'FromReaders: an odd number of fields is accepted', CL, 'if ((pairs.Count & 1) != 0) throw', 'if ((pairs.Count & 2) != 0) throw')
m('C06', 'FromReaders: a U32 array with an overflow list is accepted', CL, 'if (kind is LcpKind.U32 or LcpKind.U64 && pairs.Count != 0) throw', 'if (kind is LcpKind.U32 && pairs.Count != 0) throw')
m('C07', 'this[long] of the managed form: the rank is compared as uint', CL, 'if (_store is not null) return _store.Value(rank);\n            if ((ulong)rank >= (ulong)Length) ThrowRankOutside(rank);', 'if (_store is not null) return _store.Value(rank);\n            if ((uint)rank >= (uint)Length) ThrowRankOutside(rank);')
m('C08', 'PayloadBytes: a pair of the wide form is 8 bytes', CL, '((_store is null ? 8L : 2L * _store.PosWidth) * OverflowCount)', '((_store is null ? 8L : 8L) * OverflowCount)')
m('C09', 'OverflowCount of the wide form is 0', CL, 'if (_store is not null) return _store.PairCount;', 'if (_store is not null) return 0;')
m('C10', 'OverflowCountInt: narrowed without the check', CL, 'if (_store is not null) return IntBridge.Length(_store.PairCount, "the overflow list of the LCP array");', 'if (_store is not null) return (int)_store.PairCount;')
m('C11', 'FindOverflow of the wide form: narrowed without the check', CL, 'return found < 0 ? -1 : IntBridge.Length(found, "an index of the overflow list of the LCP array");', 'return found < 0 ? -1 : (int)found;')
m('C12', 'OverflowRankAt of the wide form: narrowed without the check', CL, 'if (_store is not null) return IntBridge.Value(_store.PairRank(index), "a rank of the overflow list of the LCP array");', 'if (_store is not null) return (int)_store.PairRank(index);')
m('C13', 'OverflowValueAt of the wide form: narrowed without the check', CL, 'if (_store is not null) return IntBridge.Value(_store.PairValue(index), "a value of the overflow list of the LCP array");', 'if (_store is not null) return (int)_store.PairValue(index);')
m('C14', 'ToArray: the length is narrowed without the check', CL, 'var result = new int[IntBridge.Length(Length, "the LCP array (ToArray)")];', 'var result = new int[(int)Length];')
m('C15', 'BaseValue of the wide form: the entry is narrowed without the check', CL, 'return _store is null ? MappedValue(rank) : IntBridge.Value(_store.BaseEntry(rank), "an LCP value");', 'return _store is null ? MappedValue(rank) : (int)_store.BaseEntry(rank);')
m('C16', 'FromMapped: a U64 array is read as U32', CL, 'if (kind == LcpKind.U64) throw new ArgumentOutOfRangeException(nameof(kind), kind, "the 32-bit form reads no U64 array: use FromMappedWide");', '')
m('C17', 'FromDense: the managed form accepts U64', CL, 'if (kind == LcpKind.U64) throw new ArgumentOutOfRangeException(nameof(kind), kind, "the managed form has no U64 kind (LCP codec 4 is read from a container)");\n        if (lcp.IsEmpty)', 'if (lcp.IsEmpty)')
m('C18', 'FromMappedWide: a width of 8 is accepted', CL, 'if (posWidth is not (4 or 5)) throw new ArgumentOutOfRangeException(nameof(posWidth), posWidth, "the stored width of a position is 4 or 5 bytes");\n        if (kind is', 'if (posWidth is not (4 or 5 or 8)) throw new ArgumentOutOfRangeException(nameof(posWidth), posWidth, "the stored width of a position is 4 or 5 bytes");\n        if (kind is')
m('C19', 'FromMappedWide: a length of 0 is accepted', CL, 'ArgumentOutOfRangeException.ThrowIfLessThan(length, 1);\n        ArgumentOutOfRangeException.ThrowIfNegative(pairCount);\n        if (posWidth', 'ArgumentOutOfRangeException.ThrowIfLessThan(length, 0);\n        ArgumentOutOfRangeException.ThrowIfNegative(pairCount);\n        if (posWidth')

LS = 'LcpStore.cs'
m('L01', 'LcpStore.BaseEntry: the bound is rank > length', LS, 'if ((ulong)rank >= (ulong)_length) throw new InvalidDataException($"the LCP array is read at rank', 'if ((ulong)rank > (ulong)_length) throw new InvalidDataException($"the LCP array is read at rank')
m('L02', 'LcpStore.BaseEntry: a value above long.MaxValue is not refused', LS, 'if (v > long.MaxValue || (_checked && v > (ulong)_length))', 'if (_checked && v > (ulong)_length)')
m('L03', 'LcpStore.BaseEntry: a value equal to the length is refused', LS, '(_checked && v > (ulong)_length)', '(_checked && v >= (ulong)_length)')
m('L04', 'LcpStore.BaseEntry: Trusted mode checks as well', LS, '(_checked && v > (ulong)_length)', '(v > (ulong)_length)')
m('L05', 'LcpStore: the escape of U8 is 65,535', LS, '_escape = kind == LcpKind.U8 ? byte.MaxValue : ushort.MaxValue;', '_escape = ushort.MaxValue;')
m('L06', 'LcpStore.Value: an escape is not resolved', LS, 'if (!_hasEscape || (ulong)b != _escape) return b;', 'if (!_hasEscape || (ulong)b == _escape) return b;')
m('L07', 'LcpStore.Find: the search moves right by 0', LS, 'if (v < (ulong)rank) lo = mid + 1;', 'if (v <= (ulong)rank) lo = mid + 1;')
m('L08', 'LcpStore.Find: the pair of rank r is read at element mid', LS, 'ulong v = _pairs.Get(2 * mid);', 'ulong v = _pairs.Get(mid);')
m('L09', 'LcpStore.PairRank: element index instead of 2 x index', LS, 'internal override long PairRank(long index) => Field(2 * index, "the overflow list");', 'internal override long PairRank(long index) => Field(index, "the overflow list");')
m('L10', 'LcpStore.PairValue: the value is read at 2 x index', LS, 'internal override long PairValue(long index) => Field((2 * index) + 1, "the overflow list");', 'internal override long PairValue(long index) => Field(2 * index, "the overflow list");')
m('L11', 'LcpStore: the pair count is a quarter of the fields', LS, '_pairCount = pairs.Count >> 1;', '_pairCount = pairs.Count >> 2;')
m('L12', 'LcpStore.Find: the middle is narrowed to int', LS, 'long mid = lo + ((hi - lo) >> 1);', 'long mid = (int)(lo + ((hi - lo) >> 1));')
m('L13', 'LcpStore.Field: a field above int.MaxValue is refused', LS, 'if (v > long.MaxValue) throw new InvalidDataException($"{what} holds {v}: the file is damaged");', 'if (v > int.MaxValue) throw new InvalidDataException($"{what} holds {v}: the file is damaged");')
m('L14', 'LcpStore.WidthOf: U32 is 8 bytes', LS, 'LcpKind.U32 => 4,', 'LcpKind.U32 => 8,')
m('L15', 'LcpStore.Value: the value of the escape is the rank field', LS, 'return PairValue(k);', 'return PairRank(k);')
m('L16', 'LcpStore.Find: the search upper bound is the count', LS, 'long lo = 0, hi = _pairCount - 1;', 'long lo = 0, hi = _pairCount;')
m('L17', 'LcpStore.BaseEntry: the value is narrowed to int', LS, 'return (long)v;\n    }\n\n    internal override long Value(long rank)', 'return (int)v;\n    }\n\n    internal override long Value(long rank)')
m('L18', 'LcpStore.Find: the sought rank is narrowed', LS, 'if (v == (ulong)rank) return mid;', 'if (v == (ulong)(int)rank) return mid;')

CC = 'ChildLayouts.cs'
m('H01', 'ChildCells: the cell of a last child is at the right bound', CC, 'ulong cell = _reader.Get(lastChild ? lb : rb);', 'ulong cell = _reader.Get(lastChild ? rb : lb);')
m('H02', 'ChildCells.NextLIndex: k = rb is not the last', CC, 'if (k >= rb) return -1;\n        long lk = _lcp[k];', 'if (k > rb) return -1;\n        long lk = _lcp[k];')
m('H03', 'ChildCells.NextLIndex: k + 1 is the next when the LCP differs', CC, 'if (_lcp[k + 1] == lk) return k + 1;\n        long v = (long)_reader.Get(k);', 'if (_lcp[k + 1] != lk) return k + 1;\n        long v = (long)_reader.Get(k);')
m('H04', 'ChildCells.NextLIndex: a cell with another LCP is the next', CC, 'long v = (long)_reader.Get(k);\n        return _lcp[v] == lk ? v : -1;', 'long v = (long)_reader.Get(k);\n        return _lcp[v] != lk ? v : -1;')
m('H05', 'ChildCells.NextLIndex (tagged): any kind but next is the next', CC, 'return (cell & TaggedKindMask) == TaggedNext ? (long)(cell & TaggedValueMask) : -1;', 'return (cell & TaggedKindMask) != TaggedNext ? (long)(cell & TaggedValueMask) : -1;')
m('H06', 'ChildCells: the tagged value has 29 bits', CC, 'private const ulong TaggedValueMask = (1UL << 30) - 1;', 'private const ulong TaggedValueMask = (1UL << 29) - 1;')
m('H07', 'ChildCells: the next kind is 2', CC, 'private const ulong TaggedNext = 1UL << 30;', 'private const ulong TaggedNext = 2UL << 30;')
m('H08', 'ChildCells: a tagged table over 5-byte cells is accepted', CC, 'if (tagged && TPos.ElementBytes != 4) throw', 'if (tagged && TPos.ElementBytes > 5) throw')
m('H09', 'MappedChild.FromMappedWide: a tagged table over 5-byte cells is accepted', CC, 'if (tagged && posWidth == 5) throw new ArgumentException', 'if (tagged && posWidth == 55) throw new ArgumentException')
m('H10', 'MappedChild.PayloadBytes: no padding', CC, 'internal long PayloadBytes => _five ? (5 * _count) + PosReader.Padding40 : 4 * _count;', 'internal long PayloadBytes => _five ? (5 * _count) : 4 * _count;')
m('H11', 'MappedChild.FirstLIndex(int) of 5-byte cells: narrowed without the check', CC, 'if (_bridged) return IntBridge.Value(FirstLIndex((long)lb, rb, lastChild), "a first l-index");', 'if (_bridged) return (int)FirstLIndex((long)lb, rb, lastChild);')
m('H12', 'MappedChild.NextLIndex(int) of 5-byte cells: narrowed without the check', CC, 'return next < 0 ? -1 : IntBridge.Value(next, "an l-index");', 'return (int)next;')
m('H13', 'MappedChild.FirstLIndex(long): Checked mode reads as Trusted', CC, 'if (_five) return _checked ? Cells(new Pos40<CheckedAccess>((byte*)_cld, _count)).FirstLIndex(lb, rb, lastChild) : Cells(new Pos40<TrustedAccess>((byte*)_cld, _count)).FirstLIndex(lb, rb, lastChild);\n        return _checked ? Cells(new Pos32<CheckedAccess>((byte*)_cld, _count)).FirstLIndex(lb, rb, lastChild)', 'if (_five) return !_checked ? Cells(new Pos40<CheckedAccess>((byte*)_cld, _count)).FirstLIndex(lb, rb, lastChild) : Cells(new Pos40<TrustedAccess>((byte*)_cld, _count)).FirstLIndex(lb, rb, lastChild);\n        return _checked ? Cells(new Pos32<CheckedAccess>((byte*)_cld, _count)).FirstLIndex(lb, rb, lastChild)')
m('H14', 'MappedChild.NextLIndex(long): the 4-byte table is read as 5-byte', CC, 'if (_five) return _checked ? Cells(new Pos40<CheckedAccess>((byte*)_cld, _count)).NextLIndex(k, rb) : Cells(new Pos40<TrustedAccess>((byte*)_cld, _count)).NextLIndex(k, rb);', 'if (!_five) return _checked ? Cells(new Pos40<CheckedAccess>((byte*)_cld, _count)).NextLIndex(k, rb) : Cells(new Pos40<TrustedAccess>((byte*)_cld, _count)).NextLIndex(k, rb);')
m('H15', 'MappedChild.Pointer: the pointer of a 5-byte table is returned', CC, 'internal int* Pointer => _five ? throw new NotSupportedException', 'internal int* Pointer => !_five ? throw new NotSupportedException')
m('H16', 'MappedChild.Cells: the table is never tagged', CC, 'where TPos : struct, IPosReader => new(reader, _lcp, _tagged);', 'where TPos : struct, IPosReader => new(reader, _lcp, false);')
m('H17', 'MappedChild: a negative count is accepted', CC, 'ArgumentOutOfRangeException.ThrowIfNegative(count);\n        _cld = cld;', '_cld = cld;')
m('H18', 'MappedChild.FromMappedWide: a width of 8 is accepted', CC, 'if (posWidth is not (4 or 5)) throw new ArgumentOutOfRangeException(nameof(posWidth), posWidth, "the stored width of a cell is 4 or 5 bytes");', 'if (posWidth is not (4 or 5 or 8)) throw new ArgumentOutOfRangeException(nameof(posWidth), posWidth, "the stored width of a cell is 4 or 5 bytes");')
m('H19', 'ChildCells.FirstLIndex: the tagged value keeps the kind bits', CC, 'return _tagged ? (long)(cell & TaggedValueMask) : (long)cell;', 'return _tagged ? (long)cell : (long)cell;')
m('H20', 'ChildCells.NextLIndex: the next of a packed cell is read at k + 1', CC, 'long v = (long)_reader.Get(k);\n        return _lcp[v] == lk', 'long v = (long)_reader.Get(k + 1);\n        return _lcp[v] == lk')
m('H21', 'MappedChild.NextLIndex(int) of 4-byte cells: the right bound is not tested', CC, 'if (k >= rb) return -1;\n        if (_checked && (ulong)k >= (ulong)_count) ThrowOutside(k);\n        int lk = _lcp[k];', 'if (k > rb) return -1;\n        if (_checked && (ulong)k >= (ulong)_count) ThrowOutside(k);\n        int lk = _lcp[k];')
m('H22', 'MappedChild: a table of more than int.MaxValue 4-byte cells is not bridged', CC, '_bridged = five || count > int.MaxValue;', '_bridged = five;')
m('H23', 'MappedChild.FirstLIndex(int): the bridge is tested for the 5-byte flag only', CC, 'public int FirstLIndex(int lb, int rb, bool lastChild)\n    {\n        if (_bridged) return', 'public int FirstLIndex(int lb, int rb, bool lastChild)\n    {\n        if (_five) return')
m('H24', 'MappedChild.NextLIndex(int): the bridge is tested for the 5-byte flag only', CC, 'public int NextLIndex(int k, int rb)\n    {\n        if (_bridged)\n        {', 'public int NextLIndex(int k, int rb)\n    {\n        if (_five)\n        {')

CT = 'ChildTable.cs'
m('T01', 'ChildTable.Rank accepts int.MaxValue + 1', CT, 'private static int Rank(long rank, string name) => (ulong)rank <= int.MaxValue ?', 'private static int Rank(long rank, string name) => (ulong)rank <= (ulong)int.MaxValue + 1 ?')
m('T02', 'ChildTable.Rank compares as uint', CT, 'private static int Rank(long rank, string name) => (ulong)rank <= int.MaxValue ?', 'private static int Rank(long rank, string name) => (uint)rank <= int.MaxValue ?')
m('T03', 'ChildTable.FirstLIndex(long): a mapped table is read as managed', CT, 'internal long FirstLIndex(long lb, long rb, bool lastChild) => IsMapped ? Mapped.FirstLIndex(lb, rb, lastChild) :', 'internal long FirstLIndex(long lb, long rb, bool lastChild) => !IsMapped ? Mapped.FirstLIndex(lb, rb, lastChild) :')
m('T04', 'ChildTable.NextLIndex(long): the right bound is not converted', CT, ': NextLIndex(Rank(k, nameof(k)), Rank(rb, nameof(rb)));', ': NextLIndex(Rank(k, nameof(k)), (int)rb);')

LP = 'LcpPyramid.cs'
m('P01', 'LcpPyramid.Bound: the bound is narrowed to int', LP, 'private static int Bound(long bound) => (int)Math.Clamp(bound, int.MinValue, int.MaxValue);', 'private static int Bound(long bound) => (int)bound;')
m('P02', 'LcpPyramid.Rank accepts int.MaxValue + 1', LP, 'private static int Rank(long rank, string name) => (ulong)rank <= int.MaxValue ?', 'private static int Rank(long rank, string name) => (ulong)rank <= (ulong)int.MaxValue + 1 ?')
m('P03', 'LcpPyramid.PrevLess(long): a position at the length is accepted', LP, 'ArgumentOutOfRangeException.ThrowIfGreaterThanOrEqual(pos, _longLength);\n        return _wideLcp is not null ? PrevLessWide', 'ArgumentOutOfRangeException.ThrowIfGreaterThan(pos, _longLength);\n        return _wideLcp is not null ? PrevLessWide')
m('P04', 'LcpPyramid.NextLess(long): a position at the length reads', LP, 'if (pos >= _longLength) return _longLength;\n        ArgumentOutOfRangeException.ThrowIfNegative(pos);\n        return _wideLcp', 'if (pos > _longLength) return _longLength;\n        ArgumentOutOfRangeException.ThrowIfNegative(pos);\n        return _wideLcp')
m('P05', 'LcpPyramid.ValueAtWide: the bound is rank > length', LP, 'if ((ulong)rank >= (ulong)_longLength) throw new ArgumentOutOfRangeException(nameof(rank));\n        return _wideLcp![rank];', 'if ((ulong)rank > (ulong)_longLength) throw new ArgumentOutOfRangeException(nameof(rank));\n        return _wideLcp![rank];')
m('P06', 'LcpPyramid.MinWide: the first partial block ends one entry early', LP, 'for (long j = lo, end = lo | (long)(f - 1); j <= end; j++) min = Math.Min(min, lcp[j]);', 'for (long j = lo, end = lo | (long)(f - 2); j <= end; j++) min = Math.Min(min, lcp[j]);')
m('P07', 'LcpPyramid.MinWide: the last partial block starts one entry late', LP, 'for (long j = hi & ~(long)(f - 1); j <= hi; j++) min = Math.Min(min, lcp[j]);', 'for (long j = (hi & ~(long)(f - 1)) + 1; j <= hi; j++) min = Math.Min(min, lcp[j]);')
m('P08', 'LcpPyramid.MinWide: the whole blocks start at the first partial one', LP, 'int l = (int)(loBlock + 1), r = (int)(hiBlock - 1);', 'int l = (int)loBlock, r = (int)(hiBlock - 1);')
m('P09', 'LcpPyramid.PrevLessWide: the block index is narrowed before the shift', LP, 'int idx = (int)(pos >> _logFanout);\n        for (int k = 1; k <= _levelCount; k++)\n        {\n            int off = _offset[k - 1];\n            for (int b = idx - 1', 'int idx = (int)pos >> _logFanout;\n        for (int k = 1; k <= _levelCount; k++)\n        {\n            int off = _offset[k - 1];\n            for (int b = idx - 1')
m('P10', 'LcpPyramid.NextLessWide: the block index is narrowed before the shift', LP, 'int idx = (int)(pos >> _logFanout);\n        for (int k = 1; k <= _levelCount; k++)\n        {\n            int off = _offset[k - 1];\n            for (int b = idx + 1', 'int idx = (int)pos >> _logFanout;\n        for (int k = 1; k <= _levelCount; k++)\n        {\n            int off = _offset[k - 1];\n            for (int b = idx + 1')
m('P11', 'LcpPyramid.RightmostBelowWide: the first position of a block is computed in int', LP, 'long first = (long)b * f;\n        for (long position = Math.Min(first + f - 1, _longLength - 1)', 'long first = b * f;\n        for (long position = Math.Min(first + f - 1, _longLength - 1)')
m('P12', 'LcpPyramid.LeftmostBelowWide: the end of a block is computed in int', LP, 'long end = Math.Min(((long)b * f) + f - 1, _longLength - 1);\n        for (long position = (long)b * f;', 'long end = Math.Min((b * f) + f - 1, _longLength - 1);\n        for (long position = (long)b * f;')
m('P13', 'LcpPyramid.BordersWide: a right border at the length reads', LP, 'rightValue = right < _longLength ? lcp[right] : -1;', 'rightValue = right <= _longLength ? lcp[right] : -1;')
m('P14', 'LcpPyramid.BordersWide: a left border of 0 has no value', LP, 'leftValue = left >= 0 ? lcp[left] : -1;', 'leftValue = left > 0 ? lcp[left] : -1;')
m('P15', 'LcpPyramid.Min(long): the right end may be at the length (non-wide form)', LP, 'if ((ulong)hi >= (ulong)_longLength) throw new ArgumentOutOfRangeException(nameof(hi)', 'if ((ulong)hi > (ulong)_longLength) throw new ArgumentOutOfRangeException(nameof(hi)')
m('P16', 'LcpPyramid.Min(long): lo = hi is refused (wide form)', LP, 'if ((ulong)lo > (ulong)hi || hi >= _longLength) throw new ArgumentOutOfRangeException(nameof(lo)', 'if ((ulong)lo >= (ulong)hi || hi >= _longLength) throw new ArgumentOutOfRangeException(nameof(lo)')
m('P17', 'LcpPyramid.Borders(long): b at the length is accepted', LP, 'if ((ulong)a >= (ulong)_longLength || b < a || b >= _longLength) throw new ArgumentOutOfRangeException(nameof(a));\n        Borders((int)a', 'if ((ulong)a >= (ulong)_longLength || b < a || b > _longLength) throw new ArgumentOutOfRangeException(nameof(a));\n        Borders((int)a')
m('P18', 'LcpPyramid.FillWide: a block minimum of int.MaxValue is refused', LP, 'if (min > int.MaxValue) throw new NotSupportedException($"the LCP array has a block whose minimum', 'if (min >= int.MaxValue) throw new NotSupportedException($"the LCP array has a block whose minimum')
m('P19', 'LcpPyramid.FillWide: the block start is computed in int', LP, 'long start = b * f;\n            long end = Math.Min(_longLength, start + f);', 'long start = (int)b * (int)f;\n            long end = Math.Min(_longLength, start + f);')
m('P20', 'LcpPyramid: the levels of a wide pyramid are not filled above the first', LP, 'FillUpperLevels();\n    }\n\n    private void FillUpperLevels()', 'if (_levelCount < 0) FillUpperLevels();\n    }\n\n    private void FillUpperLevels()')
m('P21', 'LcpPyramid.ValueAt(int): the wide guard is dropped', LP, 'if (_wideLcp is not null) return IntBridge.Value(ValueAtWide(rank), "an LCP value");', '')
m('P22', 'LcpPyramid.ctor: the sizes of the levels are rounded down', LP, 'size = (size + fanout - 1) / fanout;', 'size = size / fanout;')
m('P23', 'LcpPyramid.ctor: the block count is checked against int too late', LP, 'if (total + size > int.MaxValue) throw new NotSupportedException', 'if (total + size > long.MaxValue) throw new NotSupportedException')
m('P24', 'LcpPyramid.PrevLessWide: the bound is compared with <=', LP, 'for (long j = pos, start = pos & ~(long)(f - 1); j >= start; j--)\n            if (lcp[j] < bound) return j;', 'for (long j = pos, start = pos & ~(long)(f - 1); j >= start; j--)\n            if (lcp[j] <= bound) return j;')
m('P25', 'LcpPyramid.NextLessWide: the last block position is the length', LP, 'for (long j = pos, end = Math.Min(pos | (long)(f - 1), _longLength - 1); j <= end; j++)', 'for (long j = pos, end = Math.Min(pos | (long)(f - 1), _longLength); j <= end; j++)')
m('P26', 'LcpPyramid.PrevLess(int) of the wide form: -1 (no entry) is a damaged number', LP, 'return found < 0 ? -1 : IntBridge.Value(found, "a rank");', 'return IntBridge.Value(found, "a rank");')
m('P27', 'LcpPyramid.Borders(int) of the wide form: a left border of -1 is a damaged number', LP, 'left = l < 0 ? -1 : IntBridge.Value(l, "a rank");', 'left = IntBridge.Value(l, "a rank");')
m('P28', 'LcpPyramid.NextLess(int) of the wide form: a rank above int is a damaged file', LP, 'IntBridge.Length(NextLess((long)pos, bound), "a rank of the LCP array (NextLess)")', 'IntBridge.Value(NextLess((long)pos, bound), "a rank")')
m('P29', 'LcpPyramid.Borders(int) of the wide form: a right border above int is a damaged file', LP, 'right = IntBridge.Length(r, "a rank of the LCP array (Borders)");', 'right = IntBridge.Value(r, "a rank");')

ET = 'EsaText.cs'
m('E01', 'EsaText[long]: the bound is index > length', ET, 'if (_string is not null) return _string[checked((int)index)];\n            if ((ulong)index >= (ulong)_length) ThrowOutside(index);', 'if (_string is not null) return _string[checked((int)index)];\n            if ((ulong)index > (ulong)_length) ThrowOutside(index);')
m('E17', 'EsaText[int]: the bound is index > length', ET, 'if (_string is not null) return _string[index];\n            if ((ulong)index >= (ulong)_length) ThrowOutside(index);', 'if (_string is not null) return _string[index];\n            if ((ulong)index > (ulong)_length) ThrowOutside(index);')
m('E02', 'EsaText[long]: the index is compared as uint', ET, 'if (_string is not null) return _string[checked((int)index)];\n            if ((ulong)index >= (ulong)_length) ThrowOutside(index);', 'if (_string is not null) return _string[checked((int)index)];\n            if ((uint)index >= (uint)_length) ThrowOutside(index);')
m('E03', 'EsaText.CheckRange: a range may end one past the text', ET, 'if (start < 0 || length < 0 || start > _length - length) ThrowOutside(start);', 'if (start < 0 || length < 0 || start > _length - length + 1) ThrowOutside(start);')
m('E04', 'EsaText.CheckRange: a negative length is accepted', ET, 'if (start < 0 || length < 0 || start > _length - length) ThrowOutside(start);', 'if (start < 0 || start > _length - length) ThrowOutside(start);')
m('E05', 'EsaText.Substring(long): a start at the length is refused', ET, 'if ((ulong)start > (ulong)_length) ThrowOutside(start);\n        return Substring(start, IntBridge.Length', 'if ((ulong)start >= (ulong)_length) ThrowOutside(start);\n        return Substring(start, IntBridge.Length')
m('E06', 'EsaText.Substring(long): the rest is narrowed without the check', ET, 'return Substring(start, IntBridge.Length(_length - start, "the rest of the text (Substring)"));', 'return Substring(start, (int)(_length - start));')
m('E07', 'PackedTextSource.Length: narrowed without the check', ET, 'public int Length => IntBridge.Length(text.Length, "the length of the text (ITextSource.Length; the 64-bit member is IIndexText.Length)");', 'public int Length => (int)text.Length;')
m('E08', 'PackedTextSource[int]: the index is compared as uint', ET, 'if ((ulong)index >= (ulong)text.Length) throw new IndexOutOfRangeException', 'if ((uint)index >= (uint)text.Length) throw new IndexOutOfRangeException')
m('E09', 'PackedTextSource.CheckRange: a range may end one past the text', ET, 'if (start > text.Length - length) throw new ArgumentOutOfRangeException(nameof(length)', 'if (start > text.Length - length + 1) throw new ArgumentOutOfRangeException(nameof(length)')
m('E10', 'PackedTextSource.CheckRange: a negative start is accepted', ET, 'ArgumentOutOfRangeException.ThrowIfNegative(start);\n            ArgumentOutOfRangeException.ThrowIfNegative(length);\n            if (start > text.Length - length)', 'ArgumentOutOfRangeException.ThrowIfNegative(length);\n            if (start > text.Length - length)')
m('E11', 'PackedIndexText.SymbolAt: the terminator starts one late', ET, 'return position >= text.Length ? -1 : text[position];', 'return position > text.Length ? -1 : text[position];')
m('E12', 'PackedIndexText.CopyTo: a start beyond the length by one is accepted', ET, 'ArgumentOutOfRangeException.ThrowIfGreaterThan(start, text.Length);\n            int count', 'ArgumentOutOfRangeException.ThrowIfGreaterThan(start, text.Length + 1);\n            int count')
m('E13', 'PackedIndexText.CopyTo: the count is not clipped to the text', ET, 'int count = (int)Math.Min(destination.Length, text.Length - start);\n            text.CopyTo(start, destination[..count]);', 'int count = destination.Length;\n            text.CopyTo(start, destination[..count]);')
m('E14', 'EsaText.CopyChars: the start is narrowed', ET, '_reader!.ReadSymbols(start + done, buffer[..take]);', '_reader!.ReadSymbols((int)start + done, buffer[..take]);')
m('E15', 'EsaText[long]: the packed symbol is read at the index narrowed to int', ET, 'if (_string is not null) return _string[checked((int)index)];\n            if ((ulong)index >= (ulong)_length) ThrowOutside(index);\n            return Letter(_reader!.Symbol(index));', 'if (_string is not null) return _string[checked((int)index)];\n            if ((ulong)index >= (ulong)_length) ThrowOutside(index);\n            return Letter(_reader!.Symbol((int)index));')
m('E16', 'PackedIndexText.SymbolAt: a negative position is accepted', ET, 'ArgumentOutOfRangeException.ThrowIfNegative(position);\n            return position >= text.Length', 'return position >= text.Length')

SF = 'External/Sesa/SesaFile.cs'
m('F01', 'SesaFile: a container of int.MaxValue symbols gets the 32-bit views', SF, 'bool narrow = n < int.MaxValue && posWidth == 4', 'bool narrow = n <= int.MaxValue && posWidth == 4')
m('F02', 'SesaFile: the 32-bit views for every width', SF, 'bool narrow = n < int.MaxValue && posWidth == 4', 'bool narrow = n < int.MaxValue')
m('F03', 'SesaFile: the suffix array view has n entries', SF, '_sa = SaView.FromMapped(data + parsed.Sa.Offset, n + 1, posWidth, isChecked);', '_sa = SaView.FromMapped(data + parsed.Sa.Offset, n, posWidth, isChecked);')
m('F04', 'SesaFile.Statistics: the length is narrowed without the check', SF, 'int n = IntBridge.Length(h.TextLength, "the length of the text (EsaStatistics is 32-bit until S2.4)");', 'int n = (int)h.TextLength;')
m('F05', 'SesaFile: the wide LCP array takes the width 4 for width 5', SF, '_lcp = CompactLcp.FromMappedWide(parsed.LcpKind, n + 1, data + parsed.Lcp.Offset, pairs, parsed.Lcp.Param0, posWidth, isChecked);', '_lcp = CompactLcp.FromMappedWide(parsed.LcpKind, n + 1, data + parsed.Lcp.Offset, pairs, parsed.Lcp.Param0, 4, isChecked);')
m('F06', 'SesaFile.Payload: a section above int.MaxValue is narrowed', SF, 'if (section.Length > int.MaxValue) throw new NotSupportedException', 'if (section.Length > long.MaxValue) throw new NotSupportedException')
m('F07', 'SesaFile: the text view has n + 1 symbols', SF, '_text = EsaText.FromPacked(_reader, parsed.Letters, n);', '_text = EsaText.FromPacked(_reader, parsed.Letters, n + 1);')
m('F08', 'SesaFile: the wide child table is not tagged when the file says so', SF, 'bool tagged = parsed.Child.Codec == SesaFormat.ChildCodecTagged;', 'bool tagged = false;')

SZ = 'External/Sesa/SesaVerifier.cs'
m('Z01', 'SesaVerifier.IsZero: a gap of non-zero bytes is accepted', SZ, 'IndexOfAnyExcept((byte)0) >= 0) return false;', 'IndexOfAnyExcept((byte)0) >= 0) return true;')
m('Z02', 'SesaVerifier.IsZero: the last block of a gap is not looked at', SZ, 'for (long done = 0; done < length; done += ZeroScanBlock)', 'for (long done = 0; done + ZeroScanBlock < length; done += ZeroScanBlock)')

# check every pattern
import os
bad = 0
for x in M:
    text = open(os.path.join(ROOT, x['file'])).read()
    c = text.count(x['old'])
    if c != 1:
        print('PATTERN', x['id'], 'occurs', c, 'times in', x['file'])
        bad += 1
    if x['old'] == x['new']:
        print('NO-OP', x['id'])
        bad += 1
ids = [x['id'] for x in M]
assert len(ids) == len(set(ids)), 'duplicate ids'
json.dump(M, open(OUT, 'w'), indent=1)
print(len(M), 'mutants,', bad, 'problems')
