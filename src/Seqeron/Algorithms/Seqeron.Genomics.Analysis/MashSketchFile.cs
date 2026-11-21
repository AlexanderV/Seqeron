using System.Text;

namespace Seqeron.Genomics.Analysis;

/// <summary>One sketch of a Mash <c>.msh</c> file (<c>MinHash.capnp</c> <c>ReferenceList.Reference</c>).</summary>
/// <param name="Name">The sketch name (the file name, or the sequence id with <c>-i</c>).</param>
/// <param name="Comment">The comment.</param>
/// <param name="Length">The length of the sketched sequence (<c>length64</c>).</param>
/// <param name="Hashes">The bottom hashes, ascending (<c>hashes32</c> or <c>hashes64</c>).</param>
/// <param name="Counts">The multiplicity of each hash of a read-set sketch (<c>counts32</c>), else null.</param>
public sealed record MashSketchFileReference(string Name, string Comment, ulong Length, IReadOnlyList<ulong> Hashes, IReadOnlyList<uint>? Counts = null);

/// <summary>The content of a Mash <c>.msh</c> file: the sketch parameters and the sketches (<c>MinHash.capnp</c>, written by <c>Sketch::writeToCapnp</c> of marbl/Mash v2.3).</summary>
/// <param name="KmerSize">The k-mer size (<c>kmerSize</c>).</param>
/// <param name="MinHashesPerWindow">The sketch size s (<c>minHashesPerWindow</c>).</param>
/// <param name="HashSeed">The hash seed (<c>hashSeed</c>, stored XOR its default 42).</param>
/// <param name="Alphabet">The alphabet (<c>alphabet</c>: <c>ACGT</c> for nucleotides).</param>
/// <param name="Noncanonical"><c>noncanonical</c>: the strand is preserved (<c>-n</c>).</param>
/// <param name="PreserveCase"><c>preserveCase</c> (<c>-Z</c>).</param>
/// <param name="Concatenated"><c>concatenated</c>: one sketch per file (true unless <c>-i</c>).</param>
/// <param name="Use64">The hashes are 64-bit (<c>hashes64</c>); false: 32-bit (<c>hashes32</c>). Mash derives it from the alphabet size and k (<c>alphabetSize^k &gt; 2^32</c>).</param>
/// <param name="References">The sketches.</param>
/// <param name="WindowSize"><c>windowSize</c> (0).</param>
/// <param name="Error"><c>error</c> (0).</param>
public sealed record MashSketchFileData(int KmerSize, int MinHashesPerWindow, uint HashSeed, string Alphabet, bool Noncanonical, bool PreserveCase, bool Concatenated, bool Use64,
    IReadOnlyList<MashSketchFileReference> References, int WindowSize = 0, float Error = 0f);

/// <summary>
/// The Mash sketch file (<c>mash sketch -o x.msh</c>): a Cap'n Proto message of <c>MinHash.capnp</c>, unpacked, written by <c>capnp::MallocMessageBuilder</c> + <c>writeMessageToFd</c>. This class reproduces the builder byte for byte:
/// the objects are allocated in the order of <c>Sketch::writeToCapnp</c> (root pointer, root struct, reference list struct, list of references, per reference name / comment / hashes / counts, the empty locus list,
/// the alphabet), in segments of 1 024, 1 024, 2 048, 4 096, ... words (<c>MallocMessageBuilder</c>: <c>SUGGESTED_FIRST_SEGMENT_WORDS</c> and <c>GROW_HEURISTICALLY</c>); an object that does not fit into the segment
/// of its pointer is placed in the segment last created if it fits there (<c>BuilderArena::allocate</c>) or in a new one, behind a one-word landing pad that a far pointer references
/// (<c>WireHelpers::allocate</c>). A5-10b, B09 F193; the Cap'n Proto sources are capnproto v1.0.1 <c>message.c++</c>, <c>arena.c++</c>, <c>layout.c++</c> (the library the installed mash 2.3 links).
/// </summary>
public static class MashSketchFile
{
    private const int FirstSegmentWords = 1024;
    private const uint DefaultHashSeed = 42;

    // ---------------------------------------------------------------------------------------------------- writer

    /// <summary>The bytes of the <c>.msh</c> file of <paramref name="data"/>, identical to those of <c>mash sketch</c> / <c>mash paste</c>.</summary>
    public static byte[] Write(MashSketchFileData data)
    {
        ArgumentNullException.ThrowIfNull(data);
        var b = new Builder();
        // the root pointer and the root struct (3 data words, 4 pointers)
        var root = b.InitRoot(3, 4);
        // referenceListOld @4 (pointer 0) when the seed is the default 42, else referenceList @11 (pointer 3); locusList @5 (pointer 1); alphabet @8 (pointer 2)
        int listPointer = data.HashSeed == DefaultHashSeed ? 0 : 3;
        var refList = b.InitStruct(root, listPointer, 0, 1);
        var refs = b.InitStructList(refList, 0, data.References.Count, 2, 7);
        for (int i = 0; i < data.References.Count; i++)
        {
            var r = data.References[i];
            var e = refs.Element(i);
            b.SetText(e, 2, r.Name);
            b.SetText(e, 3, r.Comment);
            b.SetData64(e, 1, r.Length);
            if (r.Hashes.Count != 0)
            {
                if (data.Use64)
                    b.SetUInt64List(e, 5, r.Hashes);
                else
                    b.SetUInt32List(e, 4, r.Hashes.Select(h => (uint)h).ToArray());
                if (r.Counts is { Count: > 0 })
                    b.SetUInt32List(e, 6, r.Counts.ToArray());
            }
        }

        var loci = b.InitStruct(root, 1, 0, 1);
        b.InitStructList(loci, 0, 0, 3, 0);
        b.SetText(root, 2, data.Alphabet);
        b.SetData32(root, 0, (uint)data.KmerSize);
        b.SetData32(root, 1, (uint)data.WindowSize);
        b.SetData32(root, 2, (uint)data.MinHashesPerWindow);
        ulong flags = (data.Concatenated ? 1UL : 0) | (data.Noncanonical ? 2UL : 0) | (data.PreserveCase ? 4UL : 0);
        b.SetData32(root, 3, (uint)flags);
        b.SetData32(root, 4, BitConverter.SingleToUInt32Bits(data.Error));
        b.SetData32(root, 5, data.HashSeed ^ DefaultHashSeed);
        return b.ToBytes();
    }

    // a struct (or composite list element) in the arena: its segment, its first word, its data words
    private readonly record struct Obj(Builder B, int Segment, int Word, int DataWords);

    private sealed class Builder
    {
        private sealed class Segment(int capacity)
        {
            public readonly ulong[] Words = new ulong[capacity];
            public int Used;
        }

        private readonly List<Segment> _segments = [];
        private int _nextSize = FirstSegmentWords;
        private int _withSpace;

        // MallocMessageBuilder::allocateSegment
        private int NewSegment(int minimum)
        {
            int size = Math.Max(minimum, _nextSize);
            _segments.Add(new Segment(size));
            if (_segments.Count == 1)
                _nextSize = size;
            else
                _nextSize += size;
            return _segments.Count - 1;
        }

        private bool TryAllocate(int segment, int amount, out int word)
        {
            var s = _segments[segment];
            if (s.Used + amount > s.Words.Length)
            {
                word = -1;
                return false;
            }

            word = s.Used;
            s.Used += amount;
            return true;
        }

        // BuilderArena::allocate
        private (int Segment, int Word) ArenaAllocate(int amount)
        {
            if (_segments.Count == 0)
            {
                int first = NewSegment(amount);
                _withSpace = first;
                TryAllocate(first, amount, out int w0);
                return (first, w0);
            }

            if (TryAllocate(_withSpace, amount, out int w))
                return (_withSpace, w);
            int created = NewSegment(amount);
            _withSpace = created;
            TryAllocate(created, amount, out int w2);
            return (created, w2);
        }

        // WireHelpers::allocate: the pointer at (refSegment, refWord) is set to point at the new object of `amount` words; returns the object's segment and first word
        // and the location of the pointer word that has to receive the type information (the pointer itself, or the landing pad of a far pointer)
        private (int Segment, int Word, int PtrSegment, int PtrWord) Allocate(int refSegment, int refWord, int amount)
        {
            if (TryAllocate(refSegment, amount, out int word))
                return (refSegment, word, refSegment, refWord);
            var (segment, pad) = ArenaAllocate(amount + 1);
            // far pointer in the original location: landing pad at (segment, pad), the object follows the pad
            ulong far = ((ulong)(uint)pad << 3) | 2UL | ((ulong)(uint)segment << 32);
            _segments[refSegment].Words[refWord] = far;
            return (segment, pad + 1, segment, pad);
        }

        private void SetPointer(int ptrSegment, int ptrWord, int targetSegment, int targetWord, ulong kindBits, ulong upper)
        {
            long offset = targetWord - (ptrWord + 1);   // from the end of the pointer; the pad of a far pointer targets the next word (offset 0)
            _segments[ptrSegment].Words[ptrWord] = ((ulong)(uint)(int)offset << 2 & 0xFFFFFFFCUL) | kindBits | (upper << 32);
        }

        public Obj InitRoot(int dataWords, int pointers)
        {
            var (seg, word) = ArenaAllocate(1);     // the root pointer
            var (s2, w2, ps, pw) = Allocate(seg, word, dataWords + pointers);
            SetPointer(ps, pw, s2, w2, 0, (ulong)(uint)dataWords | ((ulong)(uint)pointers << 16));
            return new Obj(this, s2, w2, dataWords);
        }

        public Obj InitStruct(Obj parent, int pointerIndex, int dataWords, int pointers)
        {
            int refWord = parent.Word + parent.DataWords + pointerIndex;
            var (seg, word, ps, pw) = Allocate(parent.Segment, refWord, dataWords + pointers);
            SetPointer(ps, pw, seg, word, 0, (ulong)(uint)dataWords | ((ulong)(uint)pointers << 16));
            return new Obj(this, seg, word, dataWords);
        }

        // a list of structs (composite): the tag word, then count elements; the returned Obj stands for "element 0" and is advanced with Element(i)
        public ListObj InitStructList(Obj parent, int pointerIndex, int count, int dataWords, int pointers)
        {
            int refWord = parent.Word + parent.DataWords + pointerIndex;
            int per = dataWords + pointers;
            var (seg, word, ps, pw) = Allocate(parent.Segment, refWord, 1 + count * per);
            // list pointer, element size 7 (inline composite), word count without the tag
            SetPointer(ps, pw, seg, word, 1, 7UL | ((ulong)(uint)(count * per) << 3));
            _segments[seg].Words[word] = ((ulong)(uint)count << 2) | (ulong)(uint)dataWords << 32 | (ulong)(uint)pointers << 48;   // the tag: a struct pointer whose offset field is the element count
            return new ListObj(this, seg, word + 1, dataWords, per);
        }

        public void SetText(Obj o, int pointerIndex, string text)
        {
            byte[] bytes = Encoding.UTF8.GetBytes(text);
            WriteBytes(o, pointerIndex, bytes, 2, bytes.Length + 1);
        }

        public void SetUInt32List(Obj o, int pointerIndex, uint[] values)
        {
            var bytes = new byte[values.Length * 4];
            for (int i = 0; i < values.Length; i++)
                BitConverter.TryWriteBytes(bytes.AsSpan(i * 4), values[i]);
            WriteBytes(o, pointerIndex, bytes, 4, values.Length);
        }

        public void SetUInt64List(Obj o, int pointerIndex, IReadOnlyList<ulong> values)
        {
            var bytes = new byte[values.Count * 8];
            for (int i = 0; i < values.Count; i++)
                BitConverter.TryWriteBytes(bytes.AsSpan(i * 8), values[i]);
            WriteBytes(o, pointerIndex, bytes, 5, values.Count);
        }

        // elementSize codes: 2 = byte, 4 = four bytes, 5 = eight bytes
        private void WriteBytes(Obj o, int pointerIndex, byte[] content, int elementSizeCode, int count)
        {
            int refWord = o.Word + o.DataWords + pointerIndex;
            int byteSize = elementSizeCode == 2 ? count : content.Length;
            int words = (byteSize + 7) / 8;
            var (seg, word, ps, pw) = Allocate(o.Segment, refWord, words);
            SetPointer(ps, pw, seg, word, 1, (ulong)(uint)elementSizeCode | ((ulong)(uint)count << 3));
            for (int i = 0; i < content.Length; i++)
                _segments[seg].Words[word + i / 8] |= (ulong)content[i] << (8 * (i % 8));
        }

        public void SetData32(Obj o, int index, uint value)
        {
            ref ulong w = ref _segments[o.Segment].Words[o.Word + index / 2];
            w = (w & ~(0xFFFFFFFFUL << (32 * (index % 2)))) | ((ulong)value << (32 * (index % 2)));
        }

        public void SetData64(Obj o, int index, ulong value) => _segments[o.Segment].Words[o.Word + index] = value;

        public byte[] ToBytes()
        {
            int n = _segments.Count;
            int headerWords = (1 + n + 1) / 2;
            int total = headerWords + _segments.Sum(s => s.Used);
            var bytes = new byte[total * 8];
            BitConverter.TryWriteBytes(bytes.AsSpan(0), (uint)(n - 1));
            for (int i = 0; i < n; i++)
                BitConverter.TryWriteBytes(bytes.AsSpan(4 + 4 * i), (uint)_segments[i].Used);
            int pos = headerWords * 8;
            foreach (var s in _segments)
            {
                for (int i = 0; i < s.Used; i++)
                {
                    BitConverter.TryWriteBytes(bytes.AsSpan(pos), s.Words[i]);
                    pos += 8;
                }
            }

            return bytes;
        }
    }

    // ---------------------------------------------------------------------------------------------------- reader

    /// <summary>
    /// Reads a <c>.msh</c> file (a Cap'n Proto message of <c>MinHash.capnp</c>: any segmentation, single and double far pointers). <c>referenceList</c> (written when the seed is not 42) is preferred to
    /// <c>referenceListOld</c>; the hash width is derived from the alphabet and k as Mash derives it.
    /// </summary>
    /// <exception cref="InvalidDataException">The bytes are not a well-formed message of this schema.</exception>
    public static MashSketchFileData Read(ReadOnlySpan<byte> bytes)
    {
        var m = new Message(bytes);
        var root = m.Struct(0, 0);
        uint w0 = (uint)root.Data(0), w1 = (uint)(root.Data(0) >> 32);
        ulong flags = root.Data(1);
        uint minHashes = (uint)flags;
        bool concatenated = (flags >> 32 & 1) != 0, noncanonical = (flags >> 33 & 1) != 0, preserveCase = (flags >> 34 & 1) != 0;
        ulong third = root.Data(2);
        float error = BitConverter.UInt32BitsToSingle((uint)third);
        uint seed = (uint)(third >> 32) ^ DefaultHashSeed;
        string alphabet = root.Text(2);
        var list = root.Pointer(3);
        if (list is null || list.Value.IsNull)
            list = root.Pointer(0);
        var refs = new List<MashSketchFileReference>();
        if (list is { IsNull: false } l)
        {
            var listStruct = m.StructAt(l);
            var elements = listStruct.Pointer(0);
            if (elements is { IsNull: false } e)
            {
                foreach (var r in m.CompositeList(e))
                {
                    ulong length = r.Data(1);
                    if (length == 0)
                        length = (uint)r.Data(0);
                    var h32 = r.Pointer(4);
                    var h64 = r.Pointer(5);
                    var counts = r.Pointer(6);
                    var hashes = new List<ulong>();
                    if (h64 is { IsNull: false } p64)
                    {
                        hashes.AddRange(m.PrimitiveList(p64, 8));
                    }
                    else if (h32 is { IsNull: false } p32)
                    {
                        hashes.AddRange(m.PrimitiveList(p32, 4));
                    }

                    IReadOnlyList<uint>? countList = counts is { IsNull: false } pc ? [.. m.PrimitiveList(pc, 4).Select(x => (uint)x)] : null;
                    refs.Add(new MashSketchFileReference(r.Text(2), r.Text(3), length, hashes, countList));
                }
            }
        }

        string letters = AlphabetString(alphabet.Length == 0 ? "ACGT" : alphabet, preserveCase);
        return new MashSketchFileData((int)w0, (int)minHashes, seed, letters, noncanonical, preserveCase, concatenated, Math.Pow(letters.Length, (int)w0) > 4294967296.0, refs, (int)w1, error);
    }

    /// <summary><c>setAlphabetFromString</c> + <c>getAlphabetAsString</c>: the distinct characters in byte order, upper-cased unless the case is preserved.</summary>
    public static string AlphabetString(string characters, bool preserveCase)
    {
        ArgumentNullException.ThrowIfNull(characters);
        var set = new SortedSet<char>();
        foreach (char c in characters)
            set.Add(!preserveCase && c is > (char)96 and < (char)123 ? (char)(c - 32) : c);
        return string.Concat(set);
    }

    private readonly record struct Pointer(int Segment, int Word, ulong Value)
    {
        public bool IsNull => Value == 0;
    }

    private sealed class Message
    {
        private readonly ulong[][] _segments;

        public Message(ReadOnlySpan<byte> bytes)
        {
            if (bytes.Length < 8)
                throw new InvalidDataException("A Cap'n Proto message starts with its segment table.");
            uint n = BitConverter.ToUInt32(bytes) + 1;
            if (n > 512)
                throw new InvalidDataException("Too many segments.");
            int header = (int)(4 + 4 * n);
            if (header % 8 != 0)
                header += 4;
            if (bytes.Length < header)
                throw new InvalidDataException("Truncated segment table.");
            _segments = new ulong[n][];
            long pos = header;
            for (int i = 0; i < n; i++)
            {
                long words = BitConverter.ToUInt32(bytes[(4 + 4 * i)..]);
                if (pos + words * 8 > bytes.Length)
                    throw new InvalidDataException("Truncated segment.");
                var seg = new ulong[words];
                for (long j = 0; j < words; j++)
                    seg[j] = BitConverter.ToUInt64(bytes[(int)(pos + j * 8)..]);
                _segments[i] = seg;
                pos += words * 8;
            }
        }

        private ulong WordAt(int segment, int word)
        {
            if (segment < 0 || segment >= _segments.Length || word < 0 || word >= _segments[segment].Length)
                throw new InvalidDataException("Pointer outside the message.");
            return _segments[segment][word];
        }

        // follow far pointers: the pointer whose target is the object and the object's first word
        private (ulong Ptr, int Segment, int Target) Resolve(Pointer p)
        {
            ulong w = p.Value;
            if ((w & 3) != 2)
            {
                int offset = (int)(uint)(w & 0xFFFFFFFC) >> 2;
                return (w, p.Segment, p.Word + 1 + offset);
            }

            int seg = (int)(w >> 32), pad = (int)((w >> 3) & 0x1FFFFFFF);
            ulong landing = WordAt(seg, pad);
            if ((w & 4) == 0)
            {
                int offset = (int)(uint)(landing & 0xFFFFFFFC) >> 2;
                return (landing, seg, pad + 1 + offset);
            }

            // double far: the first pad word is a far pointer to the content, the second the tag
            int seg2 = (int)(landing >> 32), off2 = (int)((landing >> 3) & 0x1FFFFFFF);
            return (WordAt(seg, pad + 1), seg2, off2);
        }

        public StructView Struct(int segment, int word) => StructAt(new Pointer(segment, word, WordAt(segment, word)));

        public StructView StructAt(Pointer p)
        {
            var (ptr, seg, target) = Resolve(p);
            if ((ptr & 3) != 0)
                throw new InvalidDataException("A struct pointer was expected.");
            int data = (int)((ptr >> 32) & 0xFFFF), pointers = (int)(ptr >> 48);
            return new StructView(this, seg, target, data, pointers);
        }

        public IEnumerable<StructView> CompositeList(Pointer p)
        {
            var (ptr, seg, target) = Resolve(p);
            if ((ptr & 3) != 1 || (ptr >> 32 & 7) != 7)
                throw new InvalidDataException("A composite list pointer was expected.");
            ulong tag = WordAt(seg, target);
            int count = (int)(uint)(tag & 0xFFFFFFFC) >> 2;
            int data = (int)((tag >> 32) & 0xFFFF), pointers = (int)(tag >> 48);
            for (int i = 0; i < count; i++)
                yield return new StructView(this, seg, target + 1 + i * (data + pointers), data, pointers);
        }

        // element sizes: 2 = bytes (Text with the terminating NUL), 4 = 4 bytes, 5 = 8 bytes
        public IEnumerable<ulong> PrimitiveList(Pointer p, int bytesPerElement)
        {
            var (ptr, seg, target) = Resolve(p);
            if ((ptr & 3) != 1)
                throw new InvalidDataException("A list pointer was expected.");
            int count = (int)(ptr >> 35);
            for (int i = 0; i < count; i++)
            {
                int bit = i * bytesPerElement * 8;
                ulong w = WordAt(seg, target + bit / 64);
                yield return bytesPerElement == 8 ? w : w >> (bit % 64) & 0xFFFFFFFFUL;
            }
        }

        public string TextOf(Pointer p)
        {
            var (ptr, seg, target) = Resolve(p);
            if ((ptr & 3) != 1)
                throw new InvalidDataException("A text pointer was expected.");
            int count = (int)(ptr >> 35);
            var bytes = new byte[Math.Max(count - 1, 0)];
            for (int i = 0; i < bytes.Length; i++)
                bytes[i] = (byte)(WordAt(seg, target + i / 8) >> (8 * (i % 8)));
            return Encoding.UTF8.GetString(bytes);
        }

        public Pointer PointerAt(int segment, int word) => new(segment, word, WordAt(segment, word));
    }

    private readonly record struct StructView(Message M, int Segment, int Word, int DataWords, int PointerCount)
    {
        public ulong Data(int index) => index < DataWords ? M.PointerAt(Segment, Word + index).Value : 0;

        public Pointer? Pointer(int index) => index < PointerCount ? M.PointerAt(Segment, Word + DataWords + index) : null;

        public string Text(int index) => Pointer(index) is { IsNull: false } p ? M.TextOf(p) : "";
    }

    private readonly record struct ListObj(Builder B, int Segment, int FirstWord, int DataWords, int PerElement)
    {
        public Obj Element(int index) => new(B, Segment, FirstWord + index * PerElement, DataWords);
    }
}
