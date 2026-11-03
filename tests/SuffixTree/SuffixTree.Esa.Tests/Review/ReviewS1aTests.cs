using SuffixTree.Esa.Tests.TestSupport;

namespace SuffixTree.Esa.Tests.Review;

/// <summary>
/// Tests added by the independent review of stage S1a (docs/refactoring/suffix-tree-upgrade/adr/REVIEW-S1a.md). Each of the first group
/// kills a mutant of the reference kernels that survived every earlier test (the mutant table is in the review); the second group pins
/// contracts that the later unsafe kernels (1o.1) must keep and that safe code only enforces through exceptions.
/// </summary>
[TestFixture]
public class ReviewS1aTests
{
    // ---- mutant killers ---------------------------------------------------------------------------------------------------------

    [Test]
    public void SuffixArrayBuffer_WithTheWrongLength_IsRejected_InBothDirections()
    {
        // mutant M40: `sa.Length < n + 1` instead of `!=` accepted an oversized buffer, and CompactNames/SolveReduced then derive n from sa.Length
        Assert.DoesNotThrow(() => SaisReference.Build(new IntSymbols(new[] { 0, 1, 0 }), 2, new int[4]));
        Assert.Throws<ArgumentException>(() => SaisReference.Build(new IntSymbols(new[] { 0, 1, 0 }), 2, new int[3]));
        Assert.Throws<ArgumentException>(() => SaisReference.Build(new IntSymbols(new[] { 0, 1, 0 }), 2, new int[5]));
        Assert.Throws<ArgumentException>(() => SaisReference.Build(new IntSymbols(new[] { 0, 1, 0 }), 2, new int[100]));
        Assert.Throws<ArgumentException>(() => SaisReference.Build(new IntSymbols(Array.Empty<int>()), 0, new int[2]));
    }

    [Test]
    public void ComputePlcp_RejectsASuffixArrayOfTheWrongLength_InBothDirections()
    {
        // mutant M41: `sa.Length < n + 1` instead of `!=`
        var reader = new IntSymbols(new[] { 1, 0, 1, 0 });
        int[] sa = SaisReference.Build(reader, 2);
        Assert.That(LcpReference.ComputePlcp(reader, sa), Has.Length.EqualTo(4));
        Assert.Throws<ArgumentException>(() => LcpReference.ComputePlcp(new IntSymbols(new[] { 1, 0, 1, 0 }), sa.AsSpan(0, 4)));
        // a valid array followed by one extra entry (sa[0] = n is right, so only the length check can reject it)
        int[] oversized = new int[sa.Length + 1];
        sa.CopyTo(oversized, 0);
        Assert.Throws<ArgumentException>(() => LcpReference.ComputePlcp(new IntSymbols(new[] { 1, 0, 1, 0 }), oversized));
    }

    [Test]
    public void AllLmsSubstringsDistinct_NeedNoRecursion()
    {
        // mutant M08: `names <= m` recursed although every name is unique; the suffix array stays correct, only the work and the level count change
        foreach (int n in new[] { 2, 3, 5, 10, 100, 1000, 5000 })
        {
            string text = StringFamilies.AllDistinct(n, 11);
            var trace = new SaisTrace();
            int[] sa = SaisReference.BuildSuffixArray(text, trace);
            Assert.That(SuffixArrayVerifier.Check(text, sa), Is.Null, $"n {n}");
            Assert.That(trace.Depth, Is.EqualTo(1), $"n {n}: all symbols are distinct, so all LMS substrings are distinct and the reduced problem is solved directly");
        }
        // the converse: repeated LMS substrings do recurse
        var repeated = new SaisTrace();
        SaisReference.BuildSuffixArray(StringFamilies.AbRepeat(1000), repeated);
        Assert.That(repeated.Depth, Is.GreaterThanOrEqualTo(2));
    }


    [Test]
    public void Build_WithAnExplicitBuffer_RejectsANegativeAlphabetSize()
    {
        // mutant M60: the argument check was tested only on the overload that allocates the array
        Assert.Throws<ArgumentOutOfRangeException>(() => SaisReference.Build(new IntSymbols(new[] { 0, 1 }), -1, new int[3]));
        Assert.Throws<ArgumentOutOfRangeException>(() => SaisReference.Build(new IntSymbols(Array.Empty<int>()), -1, new int[1]));
    }

    [Test]
    public void CompactLcpBuilder_EnforcesItsPreconditions()
    {
        // mutants M50 (fewer escapes written than announced), M51 (the same escaped rank twice), M58b (zero length)
        var announcedTwo = new CompactLcp.Builder(LcpKind.U8, 5, 2);
        announcedTwo.Set(1, 300);
        Assert.Throws<InvalidOperationException>(() => announcedTwo.Complete());

        var duplicate = new CompactLcp.Builder(LcpKind.U8, 5, 2);
        duplicate.Set(2, 300);
        Assert.Throws<InvalidOperationException>(() => duplicate.Set(2, 400));
        Assert.Throws<InvalidOperationException>(() => duplicate.Set(1, 400));

        var exact = new CompactLcp.Builder(LcpKind.U8, 5, 2);
        exact.Set(1, 300);
        exact.Set(3, 400);
        CompactLcp done = exact.Complete();
        Assert.That(done.ToArray(), Is.EqualTo(new[] { 0, 300, 0, 400, 0 }));

        Assert.Throws<ArgumentOutOfRangeException>(() => new CompactLcp.Builder(LcpKind.U8, 0, 0));
        Assert.Throws<ArgumentException>(() => new CompactLcp.Builder(LcpKind.U32, 3, 1));
    }

    // ---- contracts the optimized (pointer) kernels must keep ----------------------------------------------------------------------

    [Test]
    public void SuffixArrayBuild_DoesNotDependOnTheInitialContentOfTheBuffer()
    {
        // the kernels of 1o.1 may use GC.AllocateUninitializedArray (ADR-002 section 3): no entry may be read before it is written
        var rng = new SplitMix64(2026);
        foreach (var (name, text) in StringFamilies.Adversarial(700, 5))
        {
            PreparedText p = TextPreparation.Prepare(text);
            int[] expected = p.Width == SymbolWidth.Byte
                ? SaisReference.Build(p.ByteReader(), p.AlphabetSize)
                : SaisReference.Build(p.UShortReader(), p.AlphabetSize);
            foreach (int fill in new[] { 0, -1, int.MaxValue, int.MinValue, p.Length, p.Length + 1, 1 })
            {
                var sa = new int[p.Length + 1];
                Array.Fill(sa, fill);
                Fill(p, sa);
                if (!sa.AsSpan().SequenceEqual(expected)) Assert.Fail($"{name}: result depends on the initial buffer content {fill}");
            }
            var noisy = new int[p.Length + 1];
            for (int i = 0; i < noisy.Length; i++) noisy[i] = (int)(rng.Next() & 0x7FFFFFFF);
            Fill(p, noisy);
            if (!noisy.AsSpan().SequenceEqual(expected)) Assert.Fail($"{name}: result depends on random initial buffer content");
        }

        static void Fill(PreparedText p, int[] sa)
        {
            if (p.Width == SymbolWidth.Byte) SaisReference.Build(p.ByteReader(), p.AlphabetSize, sa);
            else SaisReference.Build(p.UShortReader(), p.AlphabetSize, sa);
        }
    }


    [Test]
    public void ComputePlcp_OnAnInvalidSuffixArray_NeverReadsOutsideTheText()
    {
        // mutant M25: the bound `j + l < n` is dead for a valid suffix array (the predecessor of j never has suffix j as a proper prefix), but it is the
        // only guard against reading text[n] when the array comes from an untrusted source (a corrupt container in S1b). The pointer kernels of 1o.1
        // have no bounds check, so they must keep the bound; garbage in must give garbage out (every PLCP[j] within [0, n - j]), never a read past the end.
        var rng = new SplitMix64(77);
        for (int i = 0; i < 4000; i++)
        {
            int n = 1 + rng.NextInt(40);
            string text = StringFamilies.Random(n, 1 + rng.NextInt(3), rng.Next());
            PreparedText p = TextPreparation.Prepare(text);
            var sa = new int[n + 1];
            sa[0] = n;
            for (int k = 1; k <= n; k++) sa[k] = k - 1;
            for (int k = n; k > 1; k--)
            {
                int r = 1 + rng.NextInt(k);
                (sa[k], sa[r]) = (sa[r], sa[k]);
            }
            int[] plcp = p.Width == SymbolWidth.Byte ? LcpReference.ComputePlcp(p.ByteReader(), sa) : LcpReference.ComputePlcp(p.UShortReader(), sa);
            for (int j = 0; j < n; j++)
                if (plcp[j] < 0 || plcp[j] > n - j) Assert.Fail($"'{text}' sa [{string.Join(',', sa)}]: PLCP[{j}] = {plcp[j]} is outside [0, {n - j}]");
        }
        // the reversed suffix array of a^n puts every suffix after the one that has it as a proper prefix
        string aaaa = StringFamilies.AllSame(50);
        PreparedText pa = TextPreparation.Prepare(aaaa);
        int[] reversed = Enumerable.Range(0, 51).Select(k => k == 0 ? 50 : k - 1).ToArray();
        Assert.DoesNotThrow(() => LcpReference.ComputePlcp(pa.ByteReader(), reversed));
    }

    [Test]
    public void ConcurrentBuilds_OverASharedPreparedText_AreIdenticalToTheSerialBuild()
    {
        // thread safety / determinism of the reference: no mutable static state, a prepared text is read-only
        var texts = new[]
        {
            StringFamilies.RandomDna(60_000, 3),
            StringFamilies.Fibonacci(40_000),
            StringFamilies.Random(30_000, 300, 4),
            StringFamilies.AbRepeat(30_001),
        };
        foreach (string text in texts)
        {
            var (saSerial, lcpSerial) = EsaPipeline.Build(text);
            string saHash = EsaPipeline.HashSa(saSerial), lcpHash = EsaPipeline.HashLcp(lcpSerial, text.Length);
            var failures = new System.Collections.Concurrent.ConcurrentBag<string>();
            Parallel.For(0, 12, new ParallelOptions { MaxDegreeOfParallelism = 4 }, _ =>
            {
                var (sa, lcp) = EsaPipeline.Build(text);
                if (EsaPipeline.HashSa(sa) != saHash || EsaPipeline.HashLcp(lcp, text.Length) != lcpHash) failures.Add("hash differs");
            });
            Assert.That(failures, Is.Empty, $"n = {text.Length}");
        }
    }

    [Test]
    public void OrderIsByUtf16CodeUnit_NotByCodePoint()
    {
        // U+E000 followed by U+10000 (D800 DC00): in code unit order the surrogate sorts BEFORE U+E000, in code point order it would sort after
        string s = "𐀀";
        Assert.That(SaisReference.BuildSuffixArray(s), Is.EqualTo(new[] { 3, 1, 2, 0 }));
        Assert.That(NaiveSuffixArray.Build(s), Is.EqualTo(new[] { 3, 1, 2, 0 }));
    }
}
