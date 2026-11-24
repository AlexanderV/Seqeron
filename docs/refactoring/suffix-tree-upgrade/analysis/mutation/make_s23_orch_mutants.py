import json, sys
R='src/SuffixTree/Algorithms/SuffixTree.Esa/'
F1='(FullyQualifiedName~SuffixTree.Esa.Tests.Views|FullyQualifiedName~Large.GeneratedViewSelfTests)'
M=[]
def m(i, desc, f, old, new): M.append(dict(id='o33_'+i, desc=desc, file=R+f, old=old, new=new, filter=F1))
m('01','Mask40 keeps 39 bits','PosReaders.cs','0xFF_FFFF_FFFFUL','0x7F_FFFF_FFFFUL')
m('02','Padding40 is 2','PosReaders.cs','internal const int Padding40 = 3;','internal const int Padding40 = 2;')
m('03','Pos40 stride 4','PosReaders.cs','_p + (5 * index)) & PosReader.Mask40','_p + (4 * index)) & PosReader.Mask40')
m('04','Pos40 bound is > count','PosReaders.cs','if (TAccess.IsChecked && (ulong)index >= (ulong)_count) PosReader.ThrowOutside(index, _count);\n        return Unsafe.ReadUnaligned<ulong>(_p + (5 * index))','if (TAccess.IsChecked && (ulong)index > (ulong)_count) PosReader.ThrowOutside(index, _count);\n        return Unsafe.ReadUnaligned<ulong>(_p + (5 * index))')
m('05','Pos16 stride 1','PosReaders.cs','Unsafe.ReadUnaligned<ushort>(_p + (2 * index))','Unsafe.ReadUnaligned<ushort>(_p + index)')
m('06','Pos32 offset narrowed to int','PosReaders.cs','Unsafe.ReadUnaligned<uint>(_p + (4 * index))','Unsafe.ReadUnaligned<uint>(_p + (int)(4 * index))')
m('07','Pos64 offset narrowed to int','PosReaders.cs','Unsafe.ReadUnaligned<ulong>(_p + (8 * index))','Unsafe.ReadUnaligned<ulong>(_p + (int)(8 * index))')
m('08','PosInt sign-extends','PosReaders.cs','public ulong Get(long index) => (uint)_values[index];','public ulong Get(long index) => (ulong)_values[index];')
m('09','SaView: a 4-byte section is wide only above uint.MaxValue','SaView.cs','else if (count > int.MaxValue) flags |= FourBytesWide;','else if (count > uint.MaxValue) flags |= FourBytesWide;')
m('10','SaView[long]: 4-byte entries sign-extended','SaView.cs','return (_flags & FiveBytes) != 0 ? (long)Load40(rank) : (uint)_pointer[rank];\n        }\n    }\n\n    /// <summary>The same read, always bounded in both modes','return (_flags & FiveBytes) != 0 ? (long)Load40(rank) : _pointer[rank];\n        }\n    }\n\n    /// <summary>The same read, always bounded in both modes')
m('11','SaView.Load40 without mask','SaView.cs','+ (5 * rank)) & PosReader.Mask40;','+ (5 * rank));')
m('12','SaView.CopyTo(long[]) range is >=','SaView.cs','Span<long> destination)\n    {\n        if (start < 0 || destination.Length > _count - start) ThrowOutside(start);','Span<long> destination)\n    {\n        if (start < 0 || destination.Length >= _count - start) ThrowOutside(start);')
m('13','SaView.At bound is >','SaView.cs','if ((ulong)rank >= (ulong)_count) ThrowOutside(rank);\n        return this[rank];','if ((ulong)rank > (ulong)_count) ThrowOutside(rank);\n        return this[rank];')
m('14','SaView<T>.CheckRange accepts a negative start','SaView.cs','if (start < 0 || length > _reader.Count - start) SaViewErrors.Outside(start, _reader.Count);','if (length > _reader.Count - start) SaViewErrors.Outside(start, _reader.Count);')
m('15','IntBridge.Length lets a negative length through','IntBridge.cs','if ((ulong)length > int.MaxValue) ThrowLength(length, member);','if (length > int.MaxValue) ThrowLength(length, member);')
m('16','IntBridge.Value lets a negative value through','IntBridge.cs','if ((ulong)value > int.MaxValue) ThrowValue(value, what);','if (value > int.MaxValue) ThrowValue(value, what);')
m('17','LcpStore: pair count is count/4','LcpStore.cs','_pairCount = pairs.Count >> 1;','_pairCount = pairs.Count >> 2;')
m('18','LcpStore: escape of U8 and U16 swapped','LcpStore.cs','_escape = kind == LcpKind.U8 ? byte.MaxValue : ushort.MaxValue;','_escape = kind == LcpKind.U8 ? ushort.MaxValue : byte.MaxValue;')
m('19','LcpStore.Find: the search goes the wrong way','LcpStore.cs','if (v < (ulong)rank) lo = mid + 1;\n            else hi = mid - 1;','if (v < (ulong)rank) hi = mid - 1;\n            else lo = mid + 1;')
m('20','LcpStore.PairValue reads the rank field','LcpStore.cs','Field((2 * index) + 1, "the overflow list")','Field((2 * index), "the overflow list")')
m('21','LcpStore: the Checked bound on U32/U64 values is always off','LcpStore.cs','(_checked && v > (ulong)_length)','(false && v > (ulong)_length)')
m('22','LcpStore.Value resolves non-escapes','LcpStore.cs','(ulong)b != _escape) return b;','(ulong)b == _escape) return b;')
m('23','CompactLcp wide: U64 may carry an overflow list','CompactLcp.cs','if (kind is LcpKind.U32 or LcpKind.U64 && pairCount != 0) throw new InvalidDataException("an LCP array of 32 or 64 bits per entry has no overflow list");\n        long fields','if (kind is LcpKind.U32 && pairCount != 0) throw new InvalidDataException("an LCP array of 32 or 64 bits per entry has no overflow list");\n        long fields')
m('24','CompactLcp wide: width-5 pairs read as 4 bytes','CompactLcp.cs','? FromReaders(kind, baseEntries, new Pos40<TrustedAccess>(pairs, fields), isChecked)','? FromReaders(kind, baseEntries, new Pos32<TrustedAccess>(pairs, fields), isChecked)')
m('25','SesaFile: the narrow form is used at n = int.MaxValue','External/Sesa/SesaFile.cs','bool narrow = n < int.MaxValue &&','bool narrow = n <= int.MaxValue &&')
m('26','SesaVerifier.IsZero: the gap is one narrowed span','External/Sesa/SesaVerifier.cs','(int)Math.Min(ZeroScanBlock, length - done)','(int)(length - done)')
m('27','LcpPyramid: the level total is checked against uint.MaxValue','LcpPyramid.cs','if (total + size > int.MaxValue) throw new NotSupportedException','if (total + size > uint.MaxValue) throw new NotSupportedException')
m('28','EsaText.CopyTo copies the whole destination','EsaText.cs','int count = (int)Math.Min(destination.Length, text.Length - start);','int count = destination.Length;')
m('29','EsaText.CopyTo start bound is >=','EsaText.cs','ArgumentOutOfRangeException.ThrowIfGreaterThan(start, text.Length);','ArgumentOutOfRangeException.ThrowIfGreaterThanOrEqual(start, text.Length);')
m('30','SaView: the 4-byte narrow read bridges nothing (Bridged is FiveBytes only)','SaView.cs','private const byte Bridged = FiveBytes | FourBytesWide;','private const byte Bridged = FiveBytes;')
import os
ok=True
for x in M:
    t=open(os.path.join(sys.argv[1],x['file'])).read()
    n=t.count(x['old'])
    if n!=1 or x['old']==x['new']: print('BAD',x['id'],n, x['old']==x['new']); ok=False
print(len(M),'mutants', 'all unique' if ok else 'PROBLEMS')
json.dump(M,open(sys.argv[2],'w'),indent=1)
