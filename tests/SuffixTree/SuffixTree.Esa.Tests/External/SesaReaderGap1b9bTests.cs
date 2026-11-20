using System.Buffers.Binary;
using SuffixTree.Esa.External;
using SuffixTree.Esa.Tests.TestSupport;
using SuffixTree.Indexing;

namespace SuffixTree.Esa.Tests.External;

/// <summary>
/// The test gaps of the orchestrator's verifications of 1b.8 and 1b.9a that work package 1b.9b closes: the <c>Checked</c> bound of a cell of the mapped child table (an entry that points at cell <c>count</c>, one past the section, which
/// is still memory of the process when another section follows); the equality case of the corrupt-leaf-depth guard (<c>n - SA[a] + 1 &lt;= parentDepth</c>) in the navigator and in the contract tree; and the bound of a 32-bit LCP value
/// in <c>Checked</c> mode (<c>S03</c>: a value from <c>length + 1</c> to <c>length + 5</c> is refused).
/// </summary>
[TestFixture]
public unsafe class SesaReaderGap1b9bTests
{
    // ---------------------------------------------------------------- the bound of a mapped child cell

    [TestCase(false, "mississippimississippi")]
    [TestCase(true, "mississippimississippi")]
    [TestCase(false, "abracadabra")]
    [TestCase(true, "abracadabra")]
    public void AnEntryThatPointsAtCellCount_IsRefused_InCheckedMode_ForBothCodecs(bool tagged, string text)
    {
        LcpChildCase c = LcpChildCase.Of(text);
        int count = c.N + 1;
        // four cells of "the next section" behind the table, so that a read at cell count (one past the section) is a read of memory that exists and that holds a plausible value
        var cells = new uint[count + 4];
        (tagged ? c.Tagged : c.Packed).CopyTo(cells, 0);
        for (int i = count; i < cells.Length; i++) cells[i] = (tagged ? 1u << 30 : 0u) | 2u;
        (byte[] bytes, uint[] pairs) = SesaReaderGapTests.LcpSections(c.Lcp.Kind, c.Lcp.ToArray());
        fixed (uint* cp = cells)
        fixed (byte* b = bytes)
        fixed (uint* pp = pairs)
        {
            CompactLcp lcp = CompactLcp.FromMapped(c.Lcp.Kind, count, b, pairs.Length == 0 ? null : pp, pairs.Length / 2, isChecked: true);
            var checkedTable = new MappedChild((int*)cp, count, lcp, isChecked: true, tagged: tagged);
            const string Refusal = "the child table is read at cell";
            // the first l-index of an interval that ends (or begins) at the cell behind the table
            Assert.That(() => checkedTable.FirstLIndex(0, count, false), Throws.TypeOf<InvalidDataException>().With.Message.Contains(Refusal), "FirstLIndex at rb = count");
            Assert.That(() => checkedTable.FirstLIndex(count, count + 1, true), Throws.TypeOf<InvalidDataException>().With.Message.Contains(Refusal), "FirstLIndex at lb = count");
            // the l-index that follows: an entry that names cell count
            Assert.That(() => checkedTable.NextLIndex(count, count + 2), Throws.TypeOf<InvalidDataException>().With.Message.Contains(Refusal), "NextLIndex at k = count");
            Assert.That(() => checkedTable.NextLIndex(count + 1, count + 3), Throws.TypeOf<InvalidDataException>().With.Message.Contains(Refusal), "NextLIndex at k = count + 1");
            // the last cell of the section is still read
            Assert.That(() => checkedTable.FirstLIndex(0, count - 1, false), Throws.Nothing, "FirstLIndex at the last cell");
            Assert.That(() => checkedTable.NextLIndex(count - 1, count - 1), Throws.Nothing, "NextLIndex at the last cell");
            // the trusted mode has no check (it may read what is there; here that memory exists)
            var trusted = new MappedChild((int*)cp, count, lcp, isChecked: false, tagged: tagged);
            Assert.That(() => trusted.FirstLIndex(0, count, false), Throws.Nothing);
            if (tagged) Assert.That(() => trusted.NextLIndex(count, count + 2), Throws.Nothing);
        }
    }

    // ---------------------------------------------------------------- the equality case of the corrupt-leaf-depth guard

    /// <summary>Every (internal node of depth at least 1, leaf child) of the clean tree of the text, from the in-memory navigator.</summary>
    private static List<(EsaNode Parent, int LeafRank)> LeafChildren(EsaIndex index)
    {
        EsaNavigator nav = index.CreateNavigator();
        var found = new List<(EsaNode, int)>();
        var stack = new Stack<EsaNode>();
        stack.Push(nav.Root);
        var children = new List<EsaNode>();
        while (stack.Count > 0)
        {
            EsaNode v = stack.Pop();
            nav.GetChildren(v, children);
            foreach (EsaNode child in children.ToArray())
            {
                if (child.IsLeaf)
                {
                    if (v.Depth >= 1) found.Add((v, child.Lb));
                }
                else
                {
                    stack.Push(child);
                }
            }
        }
        return found;
    }

    /// <summary>
    /// The image whose suffix array has the entries of rank <paramref name="leafRank"/> and of the rank of position <c>n + 1 - parentDepth</c> exchanged (still a permutation, so the inverse suffix array can be built): the leaf of that
    /// rank now has the depth <c>n - SA + 1</c> of exactly its parent, the equality case of the guard. Null when no such position exists.
    /// </summary>
    private static byte[]? ImageWithALeafAsDeepAsItsParent(LcpChildCase c, in EsaNode parent, int leafRank, ChildCodec codec)
    {
        int position = c.N + 1 - parent.Depth;
        if (position < 0 || position > c.N) return null;
        int other = Array.IndexOf(c.Sa, position);
        if (other < 0 || other == leafRank) return null;
        byte[] image = SesaTestFile.Image(c, codec: codec);
        OracleSection sa = SesaTestFile.Section(image, SesaFormat.KindSa);
        int at = (int)sa.Offset;
        int old = BinaryPrimitives.ReadInt32LittleEndian(image.AsSpan(at + (4 * leafRank)));
        BinaryPrimitives.WriteInt32LittleEndian(image.AsSpan(at + (4 * leafRank)), position);
        BinaryPrimitives.WriteInt32LittleEndian(image.AsSpan(at + (4 * other)), old);
        SesaTestFile.Reseal(image);
        return image;
    }

    [TestCase("abab")]
    [TestCase("abracadabra")]
    [TestCase("mississippi")]
    [TestCase("banana")]
    [TestCase("aaaaaaaaaa")]
    [TestCase("GATTACAGATTACAGGTA")]
    public void ALeafAsDeepAsItsParent_IsRefused_ByTheNavigator_ForBothCodecs(string text)
    {
        LcpChildCase c = LcpChildCase.Of(text);
        int tested = 0;
        foreach (ChildCodec codec in new[] { ChildCodec.Packed, ChildCodec.Tagged })
        {
            foreach ((EsaNode parent, int leafRank) in LeafChildren(c.Index))
            {
                byte[]? image = ImageWithALeafAsDeepAsItsParent(c, parent, leafRank, codec);
                if (image is null) continue;
                using var dir = new TempDir();
                using EsaSuffixTree tree = SesaTestFile.OpenImage(dir, image, EsaValidationLevel.Header, EsaAccessMode.Checked);
                EsaNavigator nav = tree.Index.CreateNavigator();
                EsaNode node = parent;
                Assert.That(() => nav.GetChildren(node, new List<EsaNode>()), Throws.TypeOf<InvalidDataException>().With.Message.Contains("child table"), $"'{text}' {codec}: leaf {leafRank} under the node of depth {parent.Depth} [{parent.Lb}..{parent.Rb}]");
                tested++;
            }
        }
        Assert.That(tested, Is.GreaterThan(0), $"'{text}' has a leaf under an internal node whose depth is reachable by a position");
    }

    [TestCase("abab")]
    [TestCase("abracadabra")]
    [TestCase("mississippi")]
    [TestCase("banana")]
    [TestCase("aaaaaaaaaa")]
    [TestCase("GATTACAGATTACAGGTA")]
    public void ALeafAsDeepAsItsParent_IsRefused_ByTheContractTree_ForBothCodecs(string text)
    {
        LcpChildCase c = LcpChildCase.Of(text);
        int tested = 0;
        foreach (ChildCodec codec in new[] { ChildCodec.Packed, ChildCodec.Tagged })
        {
            foreach ((EsaNode parent, int leafRank) in LeafChildren(c.Index))
            {
                byte[]? image = ImageWithALeafAsDeepAsItsParent(c, parent, leafRank, codec);
                if (image is null) continue;
                using var dir = new TempDir();
                using EsaSuffixTree tree = SesaTestFile.OpenImage(dir, image, EsaValidationLevel.Header, EsaAccessMode.Checked);
                var locus = new IndexLocus(parent.Lb, parent.Rb, parent.Depth, parent.Depth);
                var buffer = new IndexLocus[text.Length + 2];
                Assert.That(() => tree.GetChildren(locus, buffer), Throws.TypeOf<InvalidDataException>().With.Message.Contains("child table"), $"'{text}' {codec}: leaf {leafRank} under the node of depth {parent.Depth} [{parent.Lb}..{parent.Rb}]");
                tested++;
            }
        }
        Assert.That(tested, Is.GreaterThan(0));
    }

    [Test]
    public void TheCleanTree_OfTheSameTexts_HasNoLeafAsDeepAsItsParent_AndAnswersTheSameWalk()
    {
        // the control: the equality case is damage, not a property of a valid tree (every leaf is strictly deeper than its parent)
        foreach (string text in new[] { "abab", "abracadabra", "mississippi", "banana", "aaaaaaaaaa", "GATTACAGATTACAGGTA" })
        {
            LcpChildCase c = LcpChildCase.Of(text);
            foreach ((EsaNode parent, int leafRank) in LeafChildren(c.Index))
            {
                EsaNavigator nav = c.Index.CreateNavigator();
                var children = new List<EsaNode>();
                nav.GetChildren(parent, children);
                EsaNode leaf = children.Single(ch => ch.Lb == leafRank && ch.IsLeaf);
                Assert.That(leaf.Depth, Is.GreaterThan(parent.Depth), $"'{text}'");
            }
        }
    }

    // ---------------------------------------------------------------- S03: the bound of a 32-bit LCP value in Checked mode

    [Test]
    public void A32BitLcpValueAboveTheLength_IsRefusedInCheckedMode_FromLengthPlusOneOn()
    {
        const int length = 50;
        int[] values = new int[length + 12];
        int[] interesting = [0, 1, length - 1, length, length + 1, length + 2, length + 3, length + 4, length + 5, 1000, int.MaxValue, -1, int.MinValue];
        fixed (int* data = values)
        {
            var checkedReader = new Lcp32Reader<CheckedAccess>(data, length);
            var trustedReader = new Lcp32Reader<TrustedAccess>(data, length);
            foreach (int v in interesting)
            {
                values[3] = v;
                bool valid = (uint)v <= (uint)length;
                if (valid) Assert.That(checkedReader[3], Is.EqualTo(v), $"value {v} is within the length {length}");
                else Assert.That(() => checkedReader[3], Throws.TypeOf<InvalidDataException>(), $"value {v} is above the length {length}");
                Assert.That(trustedReader[3], Is.EqualTo(v), $"the trusted reader does not check ({v})");
            }
        }
    }
}
