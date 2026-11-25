using System.Buffers.Binary;
using SuffixTree.Esa.External;
using SuffixTree.Esa.Tests.External;
using SuffixTree.Esa.Tests.TestSupport;
using SuffixTree.Indexing;

namespace SuffixTree.Esa.Tests.Navigation;

/// <summary>
/// Work package S2.5, the navigation pack from the container (ADR-005 section 3.3, "The navigation pack from the container"; written before the engine: against the pack of S2.4, which always builds the inverse suffix array and the
/// pyramid in memory, the assertions about the heap and about the damaged PSI and RMQ sections are red). A container that has the PSI section and the RMQ section serves <c>Parent</c>, <c>SuffixLink</c>, <c>Lcp</c> and <c>Lca</c> with no
/// heap at all; the pack of a container equals the pack of the index that was built in memory, for every operation and every node of the small texts, in both access modes; a container without them builds what it lacks (the pyramid, the
/// inverse suffix array) as before; a damaged section is refused in Checked mode as <see cref="InvalidDataException"/> and never ends in anything else.
/// </summary>
[TestFixture]
public class ContainerPackTests
{
    private static readonly EsaAccessMode[] Modes = [EsaAccessMode.Checked, EsaAccessMode.Trusted];

    public static IEnumerable<string> Texts()
    {
        var rng = new SplitMix64(0x25B00C);
        yield return "mississippi";
        yield return "abracadabraabracadabra";
        yield return "aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa";
        yield return "abababababababababababababababababababab";
        yield return "banana";
        yield return "a";
        yield return RandomText(rng, 400, 4);
        yield return RandomText(rng, 300, 2);
        yield return RandomText(rng, 250, 20);
        string block = RandomText(rng, 60, 4);
        yield return RandomText(rng, 80, 4) + block + RandomText(rng, 80, 4) + block + new string('A', 300);
    }

    private static string RandomText(SplitMix64 rng, int n, int sigma)
    {
        var chars = new char[n];
        for (int i = 0; i < n; i++) chars[i] = (char)('A' + rng.NextInt(sigma));
        return new string(chars);
    }

    private static EsaOpenOptions Options(EsaAccessMode access, long inPlaceFrom) => new() { Validation = EsaValidationLevel.Header, Access = access, NavigationInPlaceFrom = inPlaceFrom };

    /// <summary>
    /// The index over the container <paramref name="image"/> (written to a file of its own, mapped), with the file to dispose. <paramref name="inPlaceFrom"/> is the text length from which the 32-bit navigator reads PSI in place (0: always, the
    /// path that the texts of the tests would otherwise never take; the default policy is <see cref="EsaNavigationPack.InPlaceFrom"/>).
    /// </summary>
    private static (SesaFile File, EsaIndex Index) Open(TempDir dir, byte[] image, EsaAccessMode access, long inPlaceFrom = 0)
    {
        string path = dir.File("c-" + Guid.NewGuid().ToString("N") + ".sesa");
        File.WriteAllBytes(path, image);
        SesaFile file = SesaFile.Open(path, Options(access, inPlaceFrom));
        return (file, EsaIndex.FromFile(file));
    }

    private static long LevelBytes(string text)
    {
        LcpChildCase c = LcpChildCase.Of(text);
        return SesaTestFile.RmqBytes(c.Lcp.ToArray(), SesaFormat.RmqFanout).Length;
    }

    // ------------------------------------------------------------------------------------------ the pack from the container equals the pack of the index in memory

    [Test]
    public void TheContainerWithPsiAndRmq_AnswersEveryNavigationOperationAsTheMemoryIndexDoes_AndHoldsNoHeap([ValueSource(nameof(Texts))] string text, [ValueSource(nameof(Modes))] EsaAccessMode access)
    {
        using var dir = new TempDir();
        byte[] image = SesaTestFile.Image(LcpChildCase.Of(text), psi: true, rmq: true);
        (SesaFile file, EsaIndex index) = Open(dir, image, access);
        using (file)
        {
            var w = new EsaSuffixIndex(index);
            var m = new EsaSuffixIndex(EsaIndex.Build(text));
            // the suffix link, the parent, Lcp and Lca need neither the inverse suffix array nor a pyramid in memory
            foreach (IndexLocus node in WideNavigationTests.Nodes(m))
            {
                Assert.That(w.TrySuffixLink(node, out IndexLocus link), Is.EqualTo(m.TrySuffixLink(node, out IndexLocus expected)), $"TrySuffixLink {node}");
                Assert.That(link, Is.EqualTo(expected), $"link of {node}");
                Assert.That(w.TryParent(node, out IndexLocus parent), Is.EqualTo(m.TryParent(node, out IndexLocus expectedParent)), $"TryParent {node}");
                Assert.That(parent, Is.EqualTo(expectedParent), $"parent of {node}");
            }
            int step = Math.Max(1, text.Length / 30);
            for (long a = 0; a <= text.Length; a += step)
            {
                for (long b = 0; b <= text.Length; b += step)
                {
                    Assert.That(w.Lcp(a, b), Is.EqualTo(m.Lcp(a, b)), $"Lcp({a}, {b})");
                    Assert.That(w.Lca(a, b), Is.EqualTo(m.Lca(a, b)), $"Lca({a}, {b})");
                }
            }
            Assert.That(index.Navigation.PayloadBytes, Is.EqualTo(0), "no heap ISA and no heap pyramid: PSI and RMQ are read where they lie");
            // what needs the inverse suffix array still works: it is built on the first call, below the limit of the heap
            WideNavigationTests.CompareAll("container with PSI and RMQ", m, w, text);
            Assert.That(index.Navigation.PayloadBytes, Is.EqualTo(4L * (text.Length + 1)), "the inverse suffix array is the only thing in memory, built by RankOf / TryExtendLeft");
        }
    }

    [Test]
    public void TheContainerWithPsiButNoRmq_BuildsOnlyThePyramid_AndTheOneWithoutPsiBuildsTheInverseSuffixArrayToo([ValueSource(nameof(Texts))] string text, [ValueSource(nameof(Modes))] EsaAccessMode access)
    {
        using var dir = new TempDir();
        long levels = LevelBytes(text);
        var m = new EsaSuffixIndex(EsaIndex.Build(text));

        (SesaFile file, EsaIndex index) = Open(dir, SesaTestFile.Image(LcpChildCase.Of(text), psi: true, rmq: false), access);
        using (file)
        {
            var w = new EsaSuffixIndex(index);
            foreach (IndexLocus node in WideNavigationTests.Nodes(m))
            {
                Assert.That(w.TrySuffixLink(node, out IndexLocus link), Is.EqualTo(m.TrySuffixLink(node, out IndexLocus expected)), $"TrySuffixLink {node}");
                Assert.That(link, Is.EqualTo(expected), $"link of {node}");
            }
            Assert.That(index.Navigation.PayloadBytes, Is.EqualTo(levels), "PSI in place, the pyramid built in memory (0.13 bytes per character)");
        }

        (SesaFile noPsi, EsaIndex indexNoPsi) = Open(dir, SesaTestFile.Image(LcpChildCase.Of(text), psi: false, rmq: true), access);
        using (noPsi)
        {
            var w = new EsaSuffixIndex(indexNoPsi);
            WideNavigationTests.CompareAll("container with RMQ and no PSI", m, w, text);
            Assert.That(indexNoPsi.Navigation.PayloadBytes, Is.EqualTo(4L * (text.Length + 1)), "RMQ in place, the inverse suffix array built (4 bytes per character: the S1 pack without its pyramid)");
        }

        (SesaFile neither, EsaIndex indexNeither) = Open(dir, SesaTestFile.Image(LcpChildCase.Of(text), psi: false, rmq: false), access);
        using (neither)
        {
            var w = new EsaSuffixIndex(indexNeither);
            WideNavigationTests.CompareAll("container with neither", m, w, text);
            Assert.That(indexNeither.Navigation.PayloadBytes, Is.EqualTo((4L * (text.Length + 1)) + levels), "what the pack of S1 holds, no regression for a container of S1");
        }
    }

    [Test]
    public void BelowTheInPlaceLimit_TheContainerWithPsiBuildsTheInverseSuffixArrayAsS1Does_AndTheLimitIsAtLeast([Values(false, true)] bool rmq)
    {
        // the policy of S2.5: PSI is read in place from InPlaceFrom symbols on; below it the pack is the pack of S1 (the walk along the text reads ISA in order, PSI by rank: measured 13 - 15 % slower)
        const string text = "abracadabraabracadabraabracadabra";
        using var dir = new TempDir();
        byte[] image = SesaTestFile.Image(LcpChildCase.Of(text), psi: true, rmq: rmq);
        var m = new EsaSuffixIndex(EsaIndex.Build(text));
        foreach ((long from, bool inPlace) in new[] { (EsaNavigationPack.InPlaceFrom, false), ((long)text.Length + 1, false), ((long)text.Length, true), (0L, true) })
        {
            (SesaFile file, EsaIndex index) = Open(dir, image, EsaAccessMode.Checked, from);
            using (file)
            {
                var w = new EsaSuffixIndex(index);
                foreach (IndexLocus node in WideNavigationTests.Nodes(m))
                    Assert.That(w.TrySuffixLink(node, out IndexLocus link), Is.EqualTo(m.TrySuffixLink(node, out IndexLocus expected)), $"{node}");
                Assert.That(index.Navigation.UsesPsi, Is.EqualTo(inPlace), $"InPlaceFrom {from}: PSI is read in place from n = {text.Length} symbols on");
                long ranges = rmq ? 0 : LevelBytes(text);
                Assert.That(index.Navigation.PayloadBytes, Is.EqualTo((inPlace ? 0 : 4L * (text.Length + 1)) + ranges), $"InPlaceFrom {from}");
            }
        }
    }

    [Test]
    public void ThePackOfAnIndexInMemory_HoldsTheInverseSuffixArrayAndThePyramid_AsBefore([ValueSource(nameof(Texts))] string text)
    {
        EsaIndex index = EsaIndex.Build(text);
        Assert.That(index.Navigation.PayloadBytes, Is.EqualTo((4L * (text.Length + 1)) + LevelBytes(text)));
    }

    // ------------------------------------------------------------------------------------------ the damaged sections

    /// <summary>The first internal node of depth at least 2 that has a suffix link, and the rank of its left and right bound.</summary>
    private static (IndexLocus Node, long Lb, long Rb) FirstInternalNode(EsaSuffixIndex idx)
    {
        foreach (IndexLocus node in WideNavigationTests.Nodes(idx))
        {
            if (node.Length >= 2 && node.Lo != node.Hi) return (node, node.Lo, node.Hi);
        }
        throw new InvalidOperationException("the text has no internal node of depth 2");
    }

    [Test]
    public void APsiEntryThatIsNoRank_IsRefusedAsDamage_NotReadOutside([ValueSource(nameof(Modes))] EsaAccessMode access)
    {
        const string text = "mississippimississippi";
        using var dir = new TempDir();
        byte[] image = SesaTestFile.Image(LcpChildCase.Of(text), psi: true, rmq: true);
        (IndexLocus node, long lb, _) = FirstInternalNode(new EsaSuffixIndex(EsaIndex.Build(text)));
        OracleSection psi = SesaTestFile.Section(image, SesaFormat.KindPsi);
        BinaryPrimitives.WriteUInt32LittleEndian(image.AsSpan((int)(psi.Offset + (4 * lb))), uint.MaxValue); // a rank that no array has
        (SesaFile file, EsaIndex index) = Open(dir, image, access);
        using (file)
        {
            var w = new EsaSuffixIndex(index);
            Assert.That(() => w.TrySuffixLink(node, out _), Throws.TypeOf<InvalidDataException>(), "the PSI entry of the left bound of the node is not a rank of the text");
        }
    }

    [Test]
    public void APsiPairThatIsNotAscending_IsRefusedAsDamage([ValueSource(nameof(Modes))] EsaAccessMode access)
    {
        const string text = "mississippimississippi";
        using var dir = new TempDir();
        byte[] image = SesaTestFile.Image(LcpChildCase.Of(text), psi: true, rmq: true);
        (IndexLocus node, long lb, long rb) = FirstInternalNode(new EsaSuffixIndex(EsaIndex.Build(text)));
        OracleSection psi = SesaTestFile.Section(image, SesaFormat.KindPsi);
        Span<byte> left = image.AsSpan((int)(psi.Offset + (4 * lb)), 4), right = image.AsSpan((int)(psi.Offset + (4 * rb)), 4);
        byte[] t = left.ToArray();
        right.CopyTo(left);
        t.CopyTo(right); // PSI is ascending inside the interval of a node: swapped, the right end is below the left
        (SesaFile file, EsaIndex index) = Open(dir, image, access);
        using (file)
        {
            var w = new EsaSuffixIndex(index);
            Assert.That(() => w.TrySuffixLink(node, out _), Throws.TypeOf<InvalidDataException>());
        }
    }

    [Test]
    public void AnRmqSectionWhoseEntryDoesNotDescribeThePyramidOfTheText_IsRefusedAtTheBuildOfThePack([ValueSource(nameof(Modes))] EsaAccessMode access)
    {
        const string text = "mississippimississippi";
        using var dir = new TempDir();
        byte[] good = SesaTestFile.Image(LcpChildCase.Of(text), psi: true, rmq: true);
        string[] names = ["COUNT", "ELEM_WIDTH", "CODEC", "PARAM0 (the fanout)", "PARAM1", "PARAM0 not a power of two"];
        for (int field = 0; field < names.Length; field++)
        {
            byte[] image = (byte[])good.Clone();
            SesaOracle oracle = SesaOracle.FromBytes(image);
            int entry = 128 + (64 * oracle.Get(SesaFormat.KindRmq).Index);
            switch (field)
            {
                case 0: BinaryPrimitives.WriteInt64LittleEndian(image.AsSpan(entry + 24), oracle.Get(SesaFormat.KindRmq).Count - 1); break;
                case 1: BinaryPrimitives.WriteInt32LittleEndian(image.AsSpan(entry + 32), 8); break;
                case 2: BinaryPrimitives.WriteInt32LittleEndian(image.AsSpan(entry + 36), 2); break;
                case 3: BinaryPrimitives.WriteInt64LittleEndian(image.AsSpan(entry + 40), 16); break;
                case 4: BinaryPrimitives.WriteInt64LittleEndian(image.AsSpan(entry + 48), 1); break;
                default: BinaryPrimitives.WriteInt64LittleEndian(image.AsSpan(entry + 40), 33); break;
            }
            SesaTestFile.ResealHeader(image);
            (SesaFile file, EsaIndex index) = Open(dir, image, access);
            using (file)
            {
                Assert.That(() => _ = index.Navigation, Throws.TypeOf<InvalidDataException>(), $"RMQ with a wrong {names[field]}");
            }
        }
    }

    [Test]
    public void ARmqPyramidWhoseLevelsAreAllZero_IsRefusedAsDamage_WhenALinkNeedsTheBlockOfAMinimumAboveZero([ValueSource(nameof(Modes))] EsaAccessMode access)
    {
        // u^60: the 60 suffixes of one rotation are adjacent and share long prefixes, so a block of the LCP array holds only large values, and a pyramid that says "there is an entry of 0 here" leads the walk into a block that has none
        string u = StringFamilies.RandomDna(12, 77);
        string text = string.Concat(Enumerable.Repeat(u, 60));
        using var dir = new TempDir();
        byte[] image = SesaTestFile.Image(LcpChildCase.Of(text), psi: true, rmq: true);
        OracleSection rmq = SesaTestFile.Section(image, SesaFormat.KindRmq);
        image.AsSpan((int)rmq.Offset, (int)rmq.Length).Clear();
        (SesaFile file, EsaIndex index) = Open(dir, image, access);
        using (file)
        {
            var w = new EsaSuffixIndex(index);
            var m = new EsaSuffixIndex(EsaIndex.Build(text));
            int refused = 0;
            foreach (IndexLocus node in WideNavigationTests.Nodes(m))
            {
                try
                {
                    Assert.That(w.TrySuffixLink(node, out IndexLocus link), Is.EqualTo(m.TrySuffixLink(node, out IndexLocus expected)));
                    Assert.That(link, Is.EqualTo(expected));
                }
                catch (InvalidDataException)
                {
                    refused++;
                }
            }
            Assert.That(refused, Is.GreaterThan(0), "a pyramid that lies about the blocks is seen by the walk that goes down into a block without an entry below the bound");
        }
    }

    /// <summary>One answer of every navigation operation over every node of <paramref name="idx"/>, or the type of the exception that ended it.</summary>
    private static string Outcome(EsaSuffixIndex idx, EsaIndex index, string text)
    {
        try
        {
            var sb = new System.Text.StringBuilder();
            foreach (IndexLocus node in WideNavigationTests.Nodes(idx))
            {
                sb.Append(idx.TrySuffixLink(node, out IndexLocus link)).Append(link).Append(idx.TryParent(node, out IndexLocus parent)).Append(parent).Append(';');
            }
            for (long a = 0; a <= text.Length; a += 7) sb.Append(idx.Lcp(a, text.Length - a)).Append(idx.Lca(a, text.Length - a)).Append(';');
            // the generic algorithms walk suffix links through the 32-bit navigator: a damaged PSI or RMQ must end in InvalidDataException, never in an index error of the algorithm
            string query = text[2..Math.Min(text.Length, 30)] + "q" + text[1..Math.Min(text.Length, 12)];
            sb.Append(index.FindMaximalExactMatches(query, 3).Count).Append(index.FindMaximalUniqueMatches(query, 3).Count).Append(index.FindExactMatchAnchors(query, 3).Count);
            return sb.ToString();
        }
        catch (InvalidDataException)
        {
            return "InvalidDataException";
        }
    }

    [TestCase(KindOfSection.Psi)]
    [TestCase(KindOfSection.Rmq)]
    public void EveryByteOfThePsiAndTheRmqSection_Damaged_GivesTheAnswerOfTheUndamagedContainerOrInvalidData_InCheckedMode(KindOfSection which)
    {
        const string text = "mississippimississippiabracadabra";
        using var dir = new TempDir();
        byte[] good = SesaTestFile.Image(LcpChildCase.Of(text), psi: true, rmq: true);
        int kind = which == KindOfSection.Psi ? SesaFormat.KindPsi : SesaFormat.KindRmq;
        OracleSection section = SesaTestFile.Section(good, kind);
        string baseline;
        {
            (SesaFile file, EsaIndex index) = Open(dir, good, EsaAccessMode.Checked);
            using (file) baseline = Outcome(new EsaSuffixIndex(index), index, text);
        }
        Assert.That(baseline, Is.Not.EqualTo("InvalidDataException"));
        int refused = 0;
        for (long at = 0; at < section.Length; at++)
        {
            foreach (byte mask in new byte[] { 0x01, 0x80, 0xFF })
            {
                byte[] image = (byte[])good.Clone();
                image[section.Offset + at] ^= mask;
                (SesaFile file, EsaIndex index) = Open(dir, image, EsaAccessMode.Checked);
                using (file)
                {
                    string outcome = Outcome(new EsaSuffixIndex(index), index, text);
                    if (outcome == "InvalidDataException") refused++;
                }
            }
        }
        Assert.That(refused, Is.GreaterThan(0), "some of the damage is seen: the pack reads the section and checks what it reads");
    }

    /// <summary>The two sections that the pack reads in place.</summary>
    public enum KindOfSection
    {
        Psi,
        Rmq,
    }
}
