using System.Buffers.Binary;
using System.Runtime.InteropServices;
using SuffixTree.Esa.External;
using SuffixTree.Esa.Tests.TestSupport;
using SuffixTree.Indexing;

namespace SuffixTree.Esa.Tests.External;

/// <summary>
/// The views of the reader (work package 1b.7f) one by one, on memory that the test owns: the suffix array view, the text view over a packed section, the compact LCP and the child table in their mapped forms, and the cases of the
/// whole reader that the other classes reach only through a damaged header or a damaged body (every width of the LCP, a symbol above sigma, a container that is shorter than its page). Written to close the survivors of the mutation
/// audit of the reader and kept as the unit tests of the views.
/// </summary>
[TestFixture]
public unsafe class SesaReaderGapTests
{
    // ---------------------------------------------------------------- SaView

    [Test]
    public void SaView_OfAnArray_ReadsAndCopies()
    {
        int[] data = [5, 4, 3, 2, 1, 0];
        var view = new SaView(data);
        Assert.That(view.Length, Is.EqualTo(6));
        Assert.That(view[0], Is.EqualTo(5));
        Assert.That(view[5], Is.EqualTo(0));
        Assert.That(view.Array, Is.SameAs(data));
        Assert.That(view.ToArray(2, 3), Is.EqualTo(new[] { 3, 2, 1 }));
        Assert.That(view.At(4), Is.EqualTo(1));
        Assert.That(() => view.At(6), Throws.TypeOf<InvalidDataException>());
        Assert.That(() => view.At(-1), Throws.TypeOf<InvalidDataException>());
    }

    [Test]
    public void SaView_OfMappedMemory_ChecksEveryReadAgainstTheSection()
    {
        int[] data = [10, 11, 12, 13, 14, 15, 16, 17];
        fixed (int* p = data)
        {
            var section = new SaView(p, 5, isChecked: true);
            Assert.That(section.Length, Is.EqualTo(5));
            Assert.That(section[0], Is.EqualTo(10));
            Assert.That(section[4], Is.EqualTo(14));
            Assert.That(() => section[5], Throws.TypeOf<InvalidDataException>().With.Message.Contains("rank 5"), "the count itself is outside");
            Assert.That(() => section[-1], Throws.TypeOf<InvalidDataException>());
            Assert.That(() => section[int.MinValue], Throws.TypeOf<InvalidDataException>());
            Assert.That(section.At(4), Is.EqualTo(14));
            Assert.That(() => section.At(5), Throws.TypeOf<InvalidDataException>(), "the 64-bit read is checked at the count");
            Assert.That(() => section.At(long.MaxValue), Throws.TypeOf<InvalidDataException>());
            Assert.That(() => section.At(5L + (1L << 32)), Throws.TypeOf<InvalidDataException>(), "no truncation to 32 bits");
            Assert.That(section.ToArray(1, 4), Is.EqualTo(new[] { 11, 12, 13, 14 }));
            Assert.That(section.ToArray(5, 0), Is.Empty);
            var dest = new int[3];
            section.CopyTo(2, dest);
            Assert.That(dest, Is.EqualTo(new[] { 12, 13, 14 }));
            Assert.That(() => section.CopyTo(3, new int[3]), Throws.TypeOf<InvalidDataException>(), "a bulk read that ends one past the section");
            Assert.That(() => section.CopyTo(-1, new int[2]), Throws.TypeOf<InvalidDataException>(), "a bulk read that starts below the section");
            Assert.That(() => section.CopyTo(6, []), Throws.TypeOf<InvalidDataException>());
            Assert.That(() => section.ToArray(0, 6), Throws.TypeOf<InvalidDataException>());

            var trusted = new SaView(p, 5, isChecked: false);
            Assert.That(trusted[5], Is.EqualTo(15), "trusted mode does not check the single read");
            Assert.That(() => trusted.CopyTo(3, new int[3]), Throws.TypeOf<InvalidDataException>(), "a bulk read is always checked");
        }
    }

    // ---------------------------------------------------------------- EsaText over a packed section

    private static EsaText Packed(string text, out byte[] image, out GCHandle pin)
    {
        char[] letters = [.. text.Distinct().OrderBy(c => c)];
        var symbols = new byte[text.Length];
        for (int i = 0; i < text.Length; i++) symbols[i] = (byte)Array.IndexOf(letters, text[i]);
        int bits = TextOracle.BitsOf(Math.Max(1, letters.Length));
        image = TextOracle.Pack(symbols, bits);
        pin = GCHandle.Alloc(image, GCHandleType.Pinned);
        var reader = MappedTextReader.FromMemory((byte*)pin.AddrOfPinnedObject(), image.Length, text.Length, bits);
        return EsaText.FromPacked(reader, letters, text.Length);
    }

    [Test]
    public void EsaText_OverAPackedSection_IsTheStringItPacks()
    {
        string text = StringFamilies.Random(10_000, 4, 3);
        EsaText t = Packed(text, out _, out GCHandle pin);
        try
        {
            Assert.That(t.IsPacked, Is.True);
            Assert.That(t.Length, Is.EqualTo(10_000));
            Assert.That(t[0], Is.EqualTo(text[0]));
            Assert.That(t[9_999], Is.EqualTo(text[9_999]));
            Assert.That(() => t[10_000], Throws.TypeOf<InvalidDataException>().With.Message.Contains("position 10000"), "the length itself is outside");
            Assert.That(() => t[-1], Throws.TypeOf<InvalidDataException>());
            Assert.That(t.Substring(0, 10_000), Is.EqualTo(text), "more than one chunk of 4,096 symbols");
            Assert.That(t.Substring(4_095, 4_097), Is.EqualTo(text.Substring(4_095, 4_097)), "a range across the chunk boundaries");
            Assert.That(t.Substring(9_990), Is.EqualTo(text[9_990..]));
            Assert.That(t.Substring(10_000), Is.Empty);
            Assert.That(t.Substring(10_000, 0), Is.Empty);
            Assert.That(t.Substring(5, 0), Is.Empty);
            Assert.That(() => t.Substring(10_001), Throws.TypeOf<InvalidDataException>());
            Assert.That(() => t.Substring(9_999, 2), Throws.TypeOf<InvalidDataException>(), "a range that ends one past the text");
            Assert.That(() => t.Substring(-1, 2), Throws.TypeOf<InvalidDataException>());
            Assert.That(() => t.Substring(0, -1), Throws.TypeOf<InvalidDataException>(), "a negative length");
            Assert.That(t.AsMemory(100, 50).ToString(), Is.EqualTo(text.Substring(100, 50)));
            var dest = new char[300];
            t.CopyTo(9_700, dest);
            Assert.That(new string(dest), Is.EqualTo(text.Substring(9_700, 300)));
            Assert.That(() => t.CopyTo(9_701, new char[300]), Throws.TypeOf<InvalidDataException>());

            ITextSource source = t.AsTextSource();
            Assert.That(source.Length, Is.EqualTo(10_000));
            Assert.That(source[77], Is.EqualTo(text[77]));
            Assert.That(source.Substring(20, 30), Is.EqualTo(text.Substring(20, 30)));
            Assert.That(source.Slice(20, 30).ToString(), Is.EqualTo(text.Substring(20, 30)));
            Assert.That(new string([.. source]), Is.EqualTo(text), "the enumerator");
            Assert.That(source.ToString(), Is.EqualTo(text));

            IIndexText index = t.AsIndexText();
            Assert.That(index.Length, Is.EqualTo(10_000));
            Assert.That(index.SymbolAt(0), Is.EqualTo(text[0]));
            Assert.That(index.SymbolAt(9_999), Is.EqualTo(text[9_999]));
            Assert.That(index.SymbolAt(10_000), Is.EqualTo(-1), "the terminator at the length");
            Assert.That(index.SymbolAt(long.MaxValue), Is.EqualTo(-1));
            Assert.That(() => index.SymbolAt(-1), Throws.InstanceOf<ArgumentOutOfRangeException>());
            var part = new char[40];
            Assert.That(index.CopyTo(9_980, part), Is.EqualTo(20), "clipped to the text");
            Assert.That(new string(part, 0, 20), Is.EqualTo(text[9_980..]));
            Assert.That(index.CopyTo(10_000, part), Is.Zero);
            Assert.That(index.CopyTo(0, part), Is.EqualTo(40));
            Assert.That(() => index.CopyTo(10_001, part), Throws.InstanceOf<ArgumentOutOfRangeException>());
            Assert.That(() => index.CopyTo(-1, part), Throws.InstanceOf<ArgumentOutOfRangeException>());
        }
        finally
        {
            pin.Free();
        }
    }

    [Test]
    public void EsaText_APackedSymbolWithoutALetter_IsInvalidData()
    {
        // sigma 3 on 2 bits: the symbol 3 is not a letter; the section is damaged after the header was written
        string text = "abcabcabc";
        EsaText t = Packed(text, out byte[] image, out GCHandle pin);
        try
        {
            image[0] |= 0xC0; // the first symbol becomes 3
            Assert.That(() => t[0], Throws.TypeOf<InvalidDataException>().With.Message.Contains("symbol 3"));
            Assert.That(() => t.Substring(0, 4), Throws.TypeOf<InvalidDataException>());
            Assert.That(t[1], Is.EqualTo('b'));
        }
        finally
        {
            pin.Free();
        }
    }

    [Test]
    public void EsaText_OverAString_IsTheString()
    {
        EsaText t = EsaText.FromString("banana");
        Assert.That(t.IsPacked, Is.False);
        Assert.That((t.Length, t[0], t[5]), Is.EqualTo((6, 'b', 'a')));
        Assert.That(t.Substring(1, 3), Is.EqualTo("ana"));
        Assert.That(t.Substring(2), Is.EqualTo("nana"));
        Assert.That(t.AsMemory(0, 3).ToString(), Is.EqualTo("ban"));
        var dest = new char[4];
        t.CopyTo(2, dest);
        Assert.That(new string(dest), Is.EqualTo("nana"));
        Assert.That(t.AsTextSource(), Is.InstanceOf<StringTextSource>());
        Assert.That(() => t.Reader, Throws.InstanceOf<InvalidOperationException>());
        Assert.That(() => t.Alphabet, Throws.InstanceOf<InvalidOperationException>());
    }

    // ---------------------------------------------------------------- CompactLcp over mapped sections

    private static (byte[] Base, uint[] Pairs) LcpSections(LcpKind kind, int[] values)
    {
        int escape = kind switch { LcpKind.U8 => 255, LcpKind.U16 => 65535, _ => int.MaxValue };
        var pairs = new List<uint>();
        byte[] bytes;
        switch (kind)
        {
            case LcpKind.U8:
                bytes = new byte[values.Length];
                for (int i = 0; i < values.Length; i++)
                {
                    if (values[i] >= escape)
                    {
                        bytes[i] = 255;
                        pairs.Add((uint)i);
                        pairs.Add((uint)values[i]);
                    }
                    else
                    {
                        bytes[i] = (byte)values[i];
                    }
                }
                break;
            case LcpKind.U16:
                bytes = new byte[2 * values.Length];
                for (int i = 0; i < values.Length; i++)
                {
                    if (values[i] >= escape)
                    {
                        BinaryPrimitives.WriteUInt16LittleEndian(bytes.AsSpan(2 * i), 65535);
                        pairs.Add((uint)i);
                        pairs.Add((uint)values[i]);
                    }
                    else
                    {
                        BinaryPrimitives.WriteUInt16LittleEndian(bytes.AsSpan(2 * i), (ushort)values[i]);
                    }
                }
                break;
            default:
                bytes = new byte[4 * values.Length];
                for (int i = 0; i < values.Length; i++) BinaryPrimitives.WriteInt32LittleEndian(bytes.AsSpan(4 * i), values[i]);
                break;
        }
        return (bytes, [.. pairs]);
    }

    [TestCase("U8")]
    [TestCase("U16")]
    [TestCase("U32")]
    public void CompactLcp_Mapped_ReadsEveryKind_WithItsOverflowList(string kindName)
    {
        var kind = Enum.Parse<LcpKind>(kindName);
        int[] values = kind switch
        {
            LcpKind.U8 => [0, 0, 1, 5, 300, 2, 307, 254, 255, 256, 3, 0],
            LcpKind.U16 => [0, 0, 1, 5, 70_000, 2, 70_007, 254, 255, 65534, 65535, 70_001, 3, 0],
            _ => [0, 0, 1, 5, 13, 2, 12, 11, 10, 9, 8, 7, 3, 0], // 32-bit values are checked against n + 1 = 14
        };
        (byte[] bytes, uint[] pairs) = LcpSections(kind, values);
        fixed (byte* b = bytes)
        fixed (uint* pp = pairs)
        {
            CompactLcp lcp = CompactLcp.FromMapped(kind, values.Length, b, pairs.Length == 0 ? null : pp, pairs.Length / 2, isChecked: true);
            Assert.That(lcp.IsMapped, Is.True);
            Assert.That(lcp.Kind, Is.EqualTo(kind));
            Assert.That(lcp.Length, Is.EqualTo(values.Length));
            Assert.That(lcp.OverflowCount, Is.EqualTo(pairs.Length / 2));
            for (int i = 0; i < values.Length; i++) Assert.That(lcp[i], Is.EqualTo(values[i]), $"rank {i}");
            Assert.That(lcp.ToArray(), Is.EqualTo(values));
            CompactLcp.Cursor cursor = lcp.CreateCursor();
            for (int i = 0; i < values.Length; i++) Assert.That(cursor.Get(i), Is.EqualTo(values[i]), $"cursor, rank {i}");
            CompactLcp.Cursor back = lcp.CreateCursor();
            for (int i = values.Length - 1; i >= 0; i--) Assert.That(back.Get(i), Is.EqualTo(values[i]), $"cursor backwards, rank {i}");
            Assert.That(lcp.OverflowRanks, Is.EqualTo(pairs.Where((_, i) => i % 2 == 0).Select(v => (int)v)));
            Assert.That(lcp.OverflowValues, Is.EqualTo(pairs.Where((_, i) => i % 2 == 1).Select(v => (int)v)));
            Assert.That(() => lcp[values.Length], Throws.TypeOf<InvalidDataException>());
            Assert.That(() => lcp[-1], Throws.TypeOf<InvalidDataException>());
        }
    }

    [Test]
    public void CompactLcp_Mapped_AnEscapeWithoutItsPair_OrAWordAboveIntMax_IsInvalidData()
    {
        (byte[] bytes, uint[] pairs) = LcpSections(LcpKind.U8, [0, 0, 300, 400, 2]);
        fixed (byte* b = bytes)
        fixed (uint* pp = pairs)
        {
            // the pair of rank 3 is dropped: the slot holds the escape and nothing resolves it
            CompactLcp missing = CompactLcp.FromMapped(LcpKind.U8, 5, b, pp, 1, isChecked: true);
            Assert.That(missing[2], Is.EqualTo(300));
            Assert.That(() => missing[3], Throws.TypeOf<InvalidDataException>().With.Message.Contains("rank 3"));
            Assert.That(() => missing.CreateCursor().Get(3), Throws.TypeOf<InvalidDataException>());
            // a cursor at a later rank passes an earlier escape that has no pair
            CompactLcp.Cursor c = missing.CreateCursor();
            Assert.That(() => { _ = c.Get(2); _ = c.Get(3); }, Throws.TypeOf<InvalidDataException>());
        }
        uint[] huge = [2, 0x80000000u];
        fixed (byte* b = bytes)
        fixed (uint* pp = huge)
        {
            CompactLcp damaged = CompactLcp.FromMapped(LcpKind.U8, 5, b, pp, 1, isChecked: true);
            Assert.That(() => damaged[2], Throws.TypeOf<InvalidDataException>().With.Message.Contains("overflow list"));
        }
        Assert.That(() => CompactLcp.FromMapped(LcpKind.U32, 5, null, null, 1, true), Throws.TypeOf<InvalidDataException>(), "a 32-bit array has no overflow list");
    }

    [Test]
    public void CompactLcp_Mapped_U32_ChecksTheValuesAgainstTheTextInCheckedMode_Only()
    {
        (byte[] bytes, _) = LcpSections(LcpKind.U32, [0, 0, 3, 6, 2, 2_000_000_000]);
        fixed (byte* b = bytes)
        {
            CompactLcp strict = CompactLcp.FromMapped(LcpKind.U32, 6, b, null, 0, isChecked: true);
            Assert.That(strict[3], Is.EqualTo(6), "the length itself is a value (an LCP is at most n, the array has n + 1 entries)");
            Assert.That(() => strict[5], Throws.TypeOf<InvalidDataException>().With.Message.Contains("larger than the text"));
            CompactLcp lax = CompactLcp.FromMapped(LcpKind.U32, 6, b, null, 0, isChecked: false);
            Assert.That(lax[5], Is.EqualTo(2_000_000_000));
        }
    }

    // ---------------------------------------------------------------- MappedChild

    [Test]
    public void MappedChild_AnswersLikeTheTableOfTheInMemoryIndex()
    {
        string text = "mississippimississippi";
        LcpChildCase c = LcpChildCase.Of(text);
        uint[] cells = c.Packed;
        fixed (uint* cp = cells)
        {
            (byte[] bytes, uint[] pairs) = LcpSections(c.Lcp.Kind, c.Lcp.ToArray());
            fixed (byte* b = bytes)
            fixed (uint* pp = pairs)
            {
                CompactLcp lcp = CompactLcp.FromMapped(c.Lcp.Kind, c.N + 1, b, pairs.Length == 0 ? null : pp, pairs.Length / 2, true);
                var mapped = new MappedChild((int*)cp, c.N + 1, lcp, isChecked: true);
                ChildTable reference = c.Index.Child;
                Assert.That(mapped.Length, Is.EqualTo(c.N + 1));
                Assert.That(mapped.PayloadBytes, Is.EqualTo(4L * (c.N + 1)));
                ChildTable table = ChildTable.FromMapped(mapped);
                Assert.That(table.IsMapped, Is.True);
                Assert.That(table.PayloadBytes, Is.EqualTo(4L * (c.N + 1)));
                for (int lb = 0; lb < c.N; lb++)
                {
                    for (int rb = lb + 1; rb <= c.N; rb++)
                    {
                        foreach (bool last in new[] { true, false })
                        {
                            if (lb == 0 && rb == c.N && !last) continue;
                            Assert.That(mapped.FirstLIndex(lb, rb, last), Is.EqualTo(reference.FirstLIndex(lb, rb, last)), $"first [{lb},{rb}] {last}");
                            Assert.That(table.FirstLIndex(lb, rb, last), Is.EqualTo(reference.FirstLIndex(lb, rb, last)), $"table first [{lb},{rb}] {last}");
                        }
                    }
                }
                for (int k = 1; k <= c.N; k++)
                {
                    for (int rb = k; rb <= c.N; rb++)
                    {
                        // NextLIndex is only meaningful for an l-index k of an interval that ends at rb; compare where the reference has a defined answer
                        int expected = reference.NextLIndex(k, rb);
                        Assert.That(mapped.NextLIndex(k, rb), Is.EqualTo(expected), $"next {k} {rb}");
                        Assert.That(table.NextLIndex(k, rb), Is.EqualTo(expected), $"table next {k} {rb}");
                    }
                }
                Assert.That(() => mapped.FirstLIndex(c.N + 1, c.N + 1, true), Throws.TypeOf<InvalidDataException>());
                Assert.That(() => mapped.FirstLIndex(-1, 3, true), Throws.TypeOf<InvalidDataException>());
                Assert.That(() => mapped.NextLIndex(-1, 3), Throws.TypeOf<InvalidDataException>());
                var trusted = new MappedChild((int*)cp, c.N + 1, lcp, isChecked: false);
                Assert.That(trusted.FirstLIndex(c.N, c.N, true), Is.EqualTo(reference.FirstLIndex(c.N, c.N, true)));
            }
        }
    }

    // ---------------------------------------------------------------- whole containers

    [Test]
    public void EveryWidthOfTheLcp_IsReadBack_AtTheStructureLevel_AndEqualsTheInMemoryTree()
    {
        string random = StringFamilies.Random(3_700, 4, 17);
        foreach ((string label, string text) in new[] { ("U8 with overflow", random + random[..320]), ("U16 with overflow", new string('A', 70_000) + StringFamilies.RandomDna(70_000, 4)), ("U32", new string('A', 200_000)) })
        {
            LcpChildCase c = LcpChildCase.Of(text);
            Assert.That(c.Lcp.Kind, Is.EqualTo(label.StartsWith("U8", StringComparison.Ordinal) ? LcpKind.U8 : label.StartsWith("U16", StringComparison.Ordinal) ? LcpKind.U16 : LcpKind.U32), label);
            using var dir = new TempDir();
            using EsaSuffixTree tree = SesaTestFile.OpenImage(dir, SesaTestFile.Image(c), EsaValidationLevel.Structure, sha: true);
            using EsaSuffixTree memory = EsaSuffixTree.Build(text);
            Assert.That(tree.NodeCount, Is.EqualTo(memory.NodeCount), label);
            Assert.That(tree.LongestRepeatedSubstring(), Is.EqualTo(memory.LongestRepeatedSubstring()), label);
            Assert.That(tree.CountDistinctSubstrings(), Is.EqualTo(memory.CountDistinctSubstrings()), label);
            Assert.That(tree.CountDistinctSubstringsByLength(300), Is.EqualTo(memory.CountDistinctSubstringsByLength(300)), label);
            Assert.That(tree.FindAllLongestRepeatedSubstrings().Select(t => (t.Substring, t.Positions)), Is.EqualTo(memory.FindAllLongestRepeatedSubstrings().Select(t => (t.Substring, t.Positions))), label);
            var rng = new SplitMix64(9);
            for (int i = 0; i < 100; i++)
            {
                string p = text.Substring(rng.NextInt(text.Length - 40), 1 + rng.NextInt(40));
                Assert.That(tree.CountOccurrences(p), Is.EqualTo(memory.CountOccurrences(p)), $"{label}: '{p}'");
            }
            Assert.That(tree.FindMaximalExactMatches(text[5..300], 12), Is.EqualTo(memory.FindMaximalExactMatches(text[5..300], 12)), label);
        }
    }

    [Test]
    public void ADamagedU32Lcp_InCheckedMode_IsInvalidData_NotAHang()
    {
        string text = new string('A', 200_000);
        LcpChildCase c = LcpChildCase.Of(text);
        Assert.That(c.Lcp.Kind, Is.EqualTo(LcpKind.U32));
        byte[] image = SesaTestFile.Image(c);
        OracleSection lcp = SesaOracle.FromBytes(image).Get(SesaFormat.KindLcp);
        BinaryPrimitives.WriteUInt32LittleEndian(image.AsSpan((int)lcp.Offset + (4 * 1000)), 0x7FFFFFF0);
        BinaryPrimitives.WriteUInt32LittleEndian(image.AsSpan((int)lcp.Offset + (4 * 50_000)), 0xFFFFFFFF);
        SesaTestFile.Reseal(image);
        using var dir = new TempDir();
        using EsaSuffixTree tree = SesaTestFile.OpenImage(dir, image);
        Assert.That(() => tree.CountDistinctSubstrings(), Throws.TypeOf<InvalidDataException>());
        Assert.That(() => tree.Traverse(new Visitor()), Throws.TypeOf<InvalidDataException>());
        Assert.That(() => SesaTestFile.OpenImage(dir, image, EsaValidationLevel.Structure, name: "s.sesa").Dispose(), Throws.TypeOf<InvalidDataException>());
    }

    private sealed class Visitor : ISuffixTreeVisitor
    {
        public void VisitNode(int startIndex, int endIndex, int leafCount, int childCount, int depth)
        {
        }

        public void EnterBranch(int key)
        {
        }

        public void ExitBranch()
        {
        }
    }

    [Test]
    public void ATextSymbolAboveSigma_IsInvalidData_WhereItIsRead_AtTheHeaderLevel()
    {
        LcpChildCase c = LcpChildCase.Of("abcabcabcabcab");
        byte[] image = SesaTestFile.Image(c);
        OracleSection text = SesaOracle.FromBytes(image).Get(SesaFormat.KindText);
        image[text.Offset + 1] |= 0xC0; // the symbol at position 4 (2 bits each) becomes 3, which is not below sigma = 3
        using var dir = new TempDir();
        using EsaSuffixTree tree = SesaTestFile.OpenImage(dir, image);
        Assert.That(() => tree.GetAllSuffixes(), Throws.TypeOf<InvalidDataException>());
        Assert.That(() => tree.LongestRepeatedSubstring(), Throws.TypeOf<InvalidDataException>().Or.Null);
        Assert.That(() => SesaTestFile.OpenImage(dir, image, EsaValidationLevel.Sections, name: "t.sesa").Dispose(), Throws.TypeOf<InvalidDataException>(), "the checksum of the section sees it");
    }

    [Test]
    public void AContainerInMemory_ShorterThanItsHeaderPage_IsInvalidData()
    {
        byte[] image = SesaTestFile.Image(LcpChildCase.Of("banana"));
        fixed (byte* p = image)
        {
            nint address = (nint)p;
            foreach (int length in new[] { 0, 1, 8, 100, 4095 })
                Assert.That(() => SesaFile.FromMemory((byte*)address, length, new EsaOpenOptions()), Throws.TypeOf<InvalidDataException>().With.Message.Contains("truncated"), $"{length} bytes");
            using SesaFile ok = SesaFile.FromMemory(p, image.Length, new EsaOpenOptions());
            Assert.That(ok.Header.TextLength, Is.EqualTo(6));
            int total = image.Length;
            Assert.That(() => SesaFile.FromMemory((byte*)address, total, new EsaOpenOptions { Validation = (EsaValidationLevel)7 }), Throws.TypeOf<ArgumentOutOfRangeException>());
            Assert.That(() => SesaFile.FromMemory((byte*)address, total, new EsaOpenOptions { Access = (EsaAccessMode)(-1) }), Throws.TypeOf<ArgumentOutOfRangeException>());
        }
    }

    [Test]
    public void TheModeOfTheOptions_ReachesTheViews()
    {
        using var dir = new TempDir();
        string path = SesaTestFile.Write(dir, "banana");
        using SesaFile strict = SesaFile.Open(path, new EsaOpenOptions { Access = EsaAccessMode.Checked });
        Assert.That(strict.IsChecked, Is.True);
        Assert.That(() => strict.Sa[7], Throws.TypeOf<InvalidDataException>());
        using SesaFile lax = SesaFile.Open(path, new EsaOpenOptions { Access = EsaAccessMode.Trusted });
        Assert.That(lax.IsChecked, Is.False);
        Assert.That(strict.Sa.Length, Is.EqualTo(7));
        Assert.That(strict.Sa[0], Is.EqualTo(6));
        Assert.That(strict.Statistics, Is.EqualTo(LcpChildCase.Of("banana").Statistics with { StackPeak = 0 }));
        Assert.That(strict.PsiPointer is not null, Is.True);
        Assert.That(strict.PsiPointer[0], Is.Zero);
        SesaOracle o = SesaOracle.FromBytes(File.ReadAllBytes(path));
        Assert.That(strict.PsiPointer[3], Is.EqualTo(BinaryPrimitives.ReadUInt32LittleEndian(o.Payload(SesaFormat.KindPsi).AsSpan(12))));
    }

    // ---------------------------------------------------------------- gap V10 of the audit of 1b.7e (a test of the writer, found while reading the container back)

    [Test]
    public void AnOverflowWriterThatCannotBeCreated_LeavesNoOpenWriter_OfTheBaseArray()
    {
        // mutant V10: the failure of the second writer of the section stage left the first one (the base array of the LCP) open; the file could then not be completed
        using var dir = new TempDir();
        string random = StringFamilies.Random(1_500, 4, 21);
        string text = random + random[..300];
        LcpChildCase c = LcpChildCase.Of(text);
        Assert.That(c.Lcp.OverflowCount, Is.GreaterThan(0), "the text needs an overflow list: the placer is asked for its range");
        string textPath = TextOracle.WriteFile(dir, "t.pk", c.Symbols, c.Bits);
        long childOffset = 4096 + 64;
        long reserved = childOffset + (4L * (c.N + 1)) + 64;
        string path = dir.File("c.bin");
        using MappedOutputFile file = MappedOutputFile.Create(path, reserved, new MappedFileOptions { Reservation = ReservationMode.Sparse });
        using (MappedWindowWriter w = file.CreateWriter(128, 4L * (c.N + 1), 256)) w.Append(MemoryMarshal.AsBytes(c.Sa.AsSpan()));
        using MappedTextReader reader = MappedTextReader.Open(textPath, c.N, c.Bits);
        Func<LcpPlacementRequest, LcpPlacement> place = request =>
        {
            long lcp = file.Length;
            file.Grow(lcp + (LcpWidth.Bytes(request.Kind) * (c.N + 1L)) + 8 + (8 * request.OverflowCount));
            return new LcpPlacement(file, lcp, file.Length + 10_000); // the base array is where it should be, the overflow list is outside the file
        };
        Assert.That(() => SectionLcpChildBuilder.Build(path, 128, reader, file, childOffset, place, LcpChildRun.Options(dir.Sub("s"), LcpChildRun.TinyTuning(1))), Throws.TypeOf<ArgumentOutOfRangeException>());
        Assert.That(() => file.Complete(false), Throws.Nothing, "every writer that the stage opened was closed");
    }
}
