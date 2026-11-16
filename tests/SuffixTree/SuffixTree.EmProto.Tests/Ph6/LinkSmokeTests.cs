using SuffixTree.EmProto.Ph6;
using SuffixTree.EmProto.Text;
using SuffixTree.Esa;
using SuffixTree.Esa.Bench.Links;
using SuffixTree.Esa.Tests.TestSupport;

namespace SuffixTree.EmProto.Tests.Ph6;

/// <summary>
/// Gate G3 of E9: the link kernels of prototype P1 (ADR-001 section 5) return the same suffix-link targets with the Psi written by the streaming code and the ISA from DIST as with the in-RAM tables, and both equal the
/// definition (the rank interval of the label minus its first symbol) on the small adversarial texts. The smoke test must also be able to fail: a Psi or an ISA with two swapped entries is detected.
/// </summary>
[TestFixture]
public unsafe class LinkSmokeTests
{
    private static readonly int[] Fanouts = [2, 16];

    /// <summary>Every table of the streaming code for one text: ISA by DIST, Psi from the BWT (BWT from the SA file and the packed text), Psi from the ISA.</summary>
    internal sealed record Streamed(int[] IsaDist, int[] PsiStream, int[] PsiFromIsa, int[] PsiStreamCyclic);

    internal static Streamed Stream(Ph6Case c)
    {
        using var dir = new Ph6Dir();
        string sa = c.WriteSa(dir);
        string pk = c.WritePacked(dir);
        string bwt = dir.File("x.bwt");
        BwtMeta meta;
        using (var mapped = new MappedFile(pk, writable: false))
        using (var text = new PackedText(mapped.Pointer, c.N, Math.Max(1, PackedText.BitsFor(c.Sigma)), mapped, ownsOwner: false))
            meta = BwtPass.Build(sa, text, c.Sigma, bwt, new BwtOptions { ReadBytes = 28, WindowBytes = 16 });
        PsiStream.Build(bwt, meta, dir.File("psi"), new PsiOptions { WindowBytes = 12, SymbolsPerPass = 7 });
        PsiStream.Build(bwt, meta, dir.File("psic"), new PsiOptions { WindowBytes = 64, CyclicEntryZero = true });
        IsaByDist.Build(sa, c.N, dir.File("isa"), new DistOptions { ScratchDir = dir.File("scratch"), ChunkPositions = 5 + (c.N % 7), IoBufferBytes = 24, ReadBytes = 36, WindowBytes = 20 });
        PsiFromIsa.Build(sa, dir.File("isa"), c.N, dir.File("psi2"), 20, 12);
        return new Streamed(PsiFiles.ReadInts(dir.File("isa")), PsiFiles.LoadForP1(dir.File("psi")), PsiFiles.LoadForP1(dir.File("psi2")), PsiFiles.LoadForP1(dir.File("psic")));
    }

    /// <summary>The definitional target of a node: the rank interval of the suffixes that start with the label minus its first symbol (as the tests of P1).</summary>
    internal static Func<LinkNode, LinkTarget> Definition(string text, int[] sa, int[] lcp) => node =>
    {
        int n = text.Length;
        if (node.Depth == 1) return new LinkTarget(0, n, -1, true);
        string y = text.Substring(sa[node.Lb] + 1, node.Depth - 1);
        int lo = -1, hi = -1;
        for (int r = 1; r <= n; r++)
        {
            if (!text.AsSpan(sa[r]).StartsWith(y, StringComparison.Ordinal)) continue;
            if (lo < 0) lo = r;
            else if (hi != r - 1) Assert.Fail($"'{text}': the suffixes starting with '{y}' are not contiguous");
            hi = r;
        }
        if (lo < 0) Assert.Fail($"'{text}': the label of node ({node.Lb}, {node.Rb}, {node.Depth}) minus its first symbol does not occur");
        int right = hi + 1 <= n ? lcp[hi + 1] : -1;
        return new LinkTarget(lo, hi, Math.Max(lcp[lo], right), lcp[lo] > right);
    };

    private static void CheckText(string text, string label)
    {
        Ph6Case c = Ph6Case.OfNaive(text);
        EsaIndex index = EsaIndex.Build(text, ChildLayoutKind.Tagged);
        Assert.That(index.SuffixArray, Is.EqualTo(c.Sa), $"{label}: the two suffix arrays");
        int[] lcp = NaiveLcp.Build(text, c.Sa);
        Streamed s = Stream(c);
        var def = Definition(text, c.Sa, lcp);
        long expectedNodes = BruteLcpTree.Intervals(lcp).Count(x => x.Ell >= 1);
        foreach (int fanout in Fanouts)
        {
            foreach (int[] psi in new[] { s.PsiStream, s.PsiFromIsa, s.PsiStreamCyclic })
            {
                SmokeResult r = LinkSmoke.Run(index, c.IsaRam, c.PsiRam, s.IsaDist, psi, int.MaxValue, 1, 3_000_000, fanout, p1Sampling: false, def);
                Assert.That(r.NodeCount, Is.EqualTo(expectedNodes), $"{label}: internal nodes");
                Assert.That(r.SampleCount, Is.EqualTo(expectedNodes), $"{label}: every node is checked");
                Assert.That(r.Cells, Has.Length.EqualTo(8), $"{label}: four kernels x two sources");
                Assert.That(r.Cells.All(x => x.Nodes == expectedNodes), Is.True, $"{label}: every cell covers every node");
                Assert.That(r.Mismatches, Is.EqualTo(0), $"{label}: streamed tables against the in-RAM tables");
                Assert.That(r.DefinitionMismatches, Is.EqualTo(0), $"{label}: against the definition");
                Assert.That(r.Ok, Is.True);
            }
        }
    }

    [Test]
    public void AdversarialFamilies_AllNodes_AllKernelsAllSources()
    {
        foreach (int n in new[] { 1, 2, 3, 5, 8, 17, 64, 100, 257 })
            foreach (var (name, text) in StringFamilies.Adversarial(n, 11))
            {
                if (text.Distinct().Count() > 256) continue;
                CheckText(text, $"{name} n {n}");
            }
    }

    [Test]
    public void ExhaustiveSmallStrings()
    {
        foreach (string s in ExhaustiveStrings.All(2, 9).Concat(ExhaustiveStrings.All(3, 5)).Concat(ExhaustiveStrings.All(4, 4)))
            if (s.Length > 0) CheckText(s, $"'{s}'");
    }

    [Test]
    public void RandomTexts_TablesAgree_OnASampleOfNodes()
    {
        int v = 0;
        foreach (int sigma in new[] { 2, 4, 20, 65 })
        {
            Ph6Case c = Ph6Case.Of(Ph6Case.RandomText(20000, sigma, (ulong)(40 + v++)));
            EsaIndex index = EsaIndex.Build(c.Text, ChildLayoutKind.Tagged);
            Streamed s = Stream(c);
            SmokeResult r = LinkSmoke.Run(index, c.IsaRam, c.PsiRam, s.IsaDist, s.PsiStream, 3000, 1, 3_000_000, 16, p1Sampling: false, null);
            Assert.That(r.Ok, Is.True, $"sigma {sigma}");
            Assert.That(r.SampleCount, Is.EqualTo(3000));
            Assert.That(r.Cells.Sum(x => x.GaveUp), Is.EqualTo(0));
            SmokeResult r2 = LinkSmoke.Run(index, c.IsaRam, c.PsiRam, s.IsaDist, s.PsiFromIsa, 3000, 1, 3_000_000, 16, p1Sampling: false, null);
            Assert.That(r2.Ok, Is.True, $"sigma {sigma} (Psi from ISA)");
        }
    }

    [Test]
    public void P1Sampling_IsTheSampleAOfTheBenchmark()
    {
        // LinksCommands: drawn = nodes.Sample(2 * sample + 4000, seed); half = min(sample, max(1, (drawn.Length - 4000) / 2)); sample A = drawn[..half]
        Ph6Case c = Ph6Case.Of(Ph6Case.RandomText(30000, 4, 9));
        EsaIndex index = EsaIndex.Build(c.Text, ChildLayoutKind.Tagged);
        var nodes = LinkNodeSet.Enumerate(index.Lcp);
        Streamed s = Stream(c);
        foreach (int sample in new[] { 100, 5000, 100000 })
        {
            SmokeResult r = LinkSmoke.Run(index, c.IsaRam, c.PsiRam, s.IsaDist, s.PsiStream, sample, 1, 3_000_000, 16, p1Sampling: true, null);
            int drawn = Math.Min(2 * sample + 4000, nodes.Count);
            int expected = Math.Min(sample, Math.Max(1, (drawn - 4000) / 2));
            Assert.That(r.SampleCount, Is.EqualTo(expected), $"sample {sample}");
            Assert.That(r.NodeCount, Is.EqualTo(nodes.Count));
        }
    }

    [Test]
    public void ASwappedPsiEntry_IsDetected()
    {
        Ph6Case c = Ph6Case.Of(Ph6Case.RandomText(3000, 4, 5));
        EsaIndex index = EsaIndex.Build(c.Text, ChildLayoutKind.Tagged);
        Streamed s = Stream(c);
        var bad = (int[])s.PsiStream.Clone();
        for (int i = 10; i < 2990; i += 7) (bad[i], bad[i + 1]) = (bad[i + 1], bad[i]);
        SmokeResult r = LinkSmoke.Run(index, c.IsaRam, c.PsiRam, s.IsaDist, bad, int.MaxValue, 1, 3_000_000, 16, p1Sampling: false, null);
        Assert.That(r.Ok, Is.False);
        Assert.That(r.Cells.Where(x => x.Source == EndpointSource.Psi).Sum(x => x.Mismatches), Is.GreaterThan(0));
        Assert.That(r.Cells.Where(x => x.Source == EndpointSource.Isa).Sum(x => x.Mismatches), Is.EqualTo(0), "the ISA kernels do not read Psi");

        var badIsa = (int[])s.IsaDist.Clone();
        for (int i = 10; i < 2990; i += 7) (badIsa[i], badIsa[i + 1]) = (badIsa[i + 1], badIsa[i]);
        SmokeResult r2 = LinkSmoke.Run(index, c.IsaRam, c.PsiRam, badIsa, s.PsiStream, int.MaxValue, 1, 3_000_000, 16, p1Sampling: false, null);
        Assert.That(r2.Cells.Where(x => x.Source == EndpointSource.Isa).Sum(x => x.Mismatches), Is.GreaterThan(0));
        Assert.That(r2.Cells.Where(x => x.Source == EndpointSource.Psi).Sum(x => x.Mismatches), Is.EqualTo(0), "the Psi kernels do not read the ISA");
    }

    [Test]
    public void ADefinitionThatDisagrees_IsDetected()
    {
        Ph6Case c = Ph6Case.OfNaive("abracadabra");
        EsaIndex index = EsaIndex.Build(c.Text, ChildLayoutKind.Tagged);
        Streamed s = Stream(c);
        SmokeResult r = LinkSmoke.Run(index, c.IsaRam, c.PsiRam, s.IsaDist, s.PsiStream, int.MaxValue, 1, 3_000_000, 16, p1Sampling: false, node => new LinkTarget(0, 0, 0, false));
        Assert.That(r.DefinitionMismatches, Is.GreaterThan(0));
        Assert.That(r.Ok, Is.False);
        Assert.That(r.Mismatches, Is.EqualTo(0), "the tables themselves agree");
    }

    [Test]
    public void TablesOfTheWrongLength_AreRefused()
    {
        Ph6Case c = Ph6Case.OfNaive("abracadabra");
        EsaIndex index = EsaIndex.Build(c.Text, ChildLayoutKind.Tagged);
        Streamed s = Stream(c);
        Assert.Throws<ArgumentException>(() => LinkSmoke.Run(index, c.IsaRam, c.PsiRam, s.IsaDist[..^1], s.PsiStream, 10, 1, 1000, 16, false, null));
        Assert.Throws<ArgumentException>(() => LinkSmoke.Run(index, c.IsaRam, c.PsiRam, s.IsaDist, s.PsiStream[..^1], 10, 1, 1000, 16, false, null));
    }
}
