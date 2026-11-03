using SuffixTree.Esa.Tests.TestSupport;

namespace SuffixTree.Esa.Tests.Kernels;

/// <summary>
/// The kernel-vs-reference registry (work package 1o.0): every registered engine runs through <see cref="KernelHarness"/> on the default corpus, the
/// intermediate-state oracle is checked against the reference, and the harness is shown to catch each class of kernel defect with fault-injected
/// test engines (wrong SA, wrong trace, wrong bucket arrays / LMS order / names, wrong LCP, wrong PLCP, out-of-bounds write, read of an uninitialised cell).
/// </summary>
[TestFixture]
public class KernelRegistryTests
{
    private const string Banana = "banana";
    private static readonly string Fib = StringFamilies.Fibonacci(200);

    public static IEnumerable<TestCaseData> SaisEngines() => KernelRegistry.SaisEngines.Select(e => new TestCaseData(e.Name).SetName("Sais_" + e.Name));

    public static IEnumerable<TestCaseData> LcpEngines() => KernelRegistry.LcpEngines.Select(e => new TestCaseData(e.Name).SetName("Lcp_" + e.Name));

    [TestCaseSource(nameof(SaisEngines))]
    public void EverySaisEngine_EqualsTheReference_OnTheDefaultCorpus(string engineName)
    {
        SaisEngine engine = KernelRegistry.SaisEngines.Single(e => e.Name == engineName);
        int count = 0;
        foreach (var (name, text) in KernelCorpus.Default())
        {
            string? failure = KernelHarness.CheckSais(engine, text, allPoisons: text.Length <= 2000);
            if (failure is not null) Assert.Fail($"{engine.Name} on {name} (n = {text.Length}): {failure}");
            count++;
        }
        TestContext.Out.WriteLine($"{engine.Name}: {count} texts equal to SaisReference, the naive SA / verifier and the trace, canaries intact");
    }

    [TestCaseSource(nameof(LcpEngines))]
    public void EveryLcpEngine_EqualsTheReference_OnTheDefaultCorpus(string engineName)
    {
        LcpEngine engine = KernelRegistry.LcpEngines.Single(e => e.Name == engineName);
        int count = 0;
        foreach (var (name, text) in KernelCorpus.Default())
        {
            string? failure = KernelHarness.CheckLcp(engine, text, allPoisons: text.Length <= 2000);
            if (failure is not null) Assert.Fail($"{engine.Name} on {name} (n = {text.Length}): {failure}");
            count++;
        }
        TestContext.Out.WriteLine($"{engine.Name}: {count} texts equal to LcpReference, the naive LCP / Kasai, canaries intact");
    }

    // ---- the state oracle itself ----

    [Test]
    public void TheStateOracle_AgreesWithTheReferenceTrace_OnTheDefaultCorpus()
    {
        foreach (var (name, text) in KernelCorpus.Default())
        {
            PreparedText prepared = TextPreparation.Prepare(text);
            using var arena = new NativeArena();
            SaisRun reference = KernelRegistry.RunReferenceSais(prepared, arena);
            List<SaisLevelState> oracle = SaisStateOracle.Compute(KernelHarness.Symbols(prepared), prepared.AlphabetSize);
            Assert.That(oracle, Has.Count.EqualTo(reference.Levels.Count), name);
            for (int i = 0; i < oracle.Count; i++)
                Assert.That(new SaisLevel(oracle[i].Length, oracle[i].Sigma, oracle[i].LmsPositions.Length), Is.EqualTo(reference.Levels[i]), $"{name} level {i}");
        }
    }

    [Test]
    public void TheStateOracle_NamesAreOrderPreserving_TheReducedSuffixOrderIsTheOrderOfTheLmsSuffixes()
    {
        // the property that makes SA-IS correct: sorting the LMS suffixes = sorting the reduced string; checked against the naive SA of the text
        foreach (var (name, text) in KernelCorpus.Default().Where(c => c.Text.Length is > 1 and <= 2000).Take(3000))
        {
            PreparedText prepared = TextPreparation.Prepare(text);
            SaisLevelState level0 = SaisStateOracle.Compute(KernelHarness.Symbols(prepared), prepared.AlphabetSize)[0];
            if (level0.LmsPositions.Length < 2) continue;
            string reduced = new(level0.Reduced.Select(v => (char)v).ToArray());
            int[] reducedSa = NaiveSuffixArray.Build(reduced);
            int[] viaReduced = reducedSa.Skip(1).Select(r => level0.LmsPositions[r]).ToArray();
            int[] naive = NaiveSuffixArray.Build(text);
            var lms = new HashSet<int>(level0.LmsPositions);
            int[] viaText = naive.Where(lms.Contains).ToArray();
            Assert.That(viaReduced, Is.EqualTo(viaText), name);
        }
    }

    [Test]
    public void TheStateOracle_KnownExample_Banana()
    {
        PreparedText prepared = TextPreparation.Prepare(Banana); // a=0 b=1 n=2: b a n a n a; types L S L S L L; LMS positions 1, 3
        List<SaisLevelState> states = SaisStateOracle.Compute(KernelHarness.Symbols(prepared), prepared.AlphabetSize);
        SaisLevelState l0 = states[0];
        Assert.That(l0.Counts, Is.EqualTo(new[] { 3, 1, 2 }));
        Assert.That(l0.Heads, Is.EqualTo(new[] { 1, 4, 5 }));
        Assert.That(l0.Tails, Is.EqualTo(new[] { 4, 5, 7 }));
        Assert.That(l0.LmsPositions, Is.EqualTo(new[] { 1, 3 }));
        Assert.That(l0.NameCount, Is.EqualTo(2)); // "ana" (positions 1..3) and "ana$" (3..6, ends in the sentinel) differ
        Assert.That(states, Has.Count.EqualTo(1));
    }

    // ---- the harness catches each class of defect ----

    [Test]
    public void AWrongSuffixArray_IsReported()
        => Assert.That(KernelHarness.CheckSais(Faulty.Sais(sa => (sa[1], sa[2]) = (sa[2], sa[1])), Banana), Does.StartWith("poison").And.Contain("SA["));

    [Test]
    public void AWrongRecursionTrace_IsReported()
        => Assert.That(KernelHarness.CheckSais(Faulty.SaisTrace(l => l with { LmsCount = l.LmsCount + 1 }), Fib), Does.Contain("trace level"));

    [Test]
    public void AMissingRecursionLevel_IsReported()
        => Assert.That(KernelHarness.CheckSais(Faulty.SaisTrace(null), Fib), Does.Contain("recursion levels"));

    [Test]
    public void ACorrectStatesReport_IsAccepted_AndEachWrongStateIsNamed()
    {
        Assert.That(KernelHarness.CheckSais(Faulty.States(s => s), Fib), Is.Null);
        Assert.That(KernelHarness.CheckSais(Faulty.States(s => s with { Counts = [.. s.Counts.Select((c, i) => i == 0 ? c + 1 : c)] }), Fib), Does.Contain("symbol counts"));
        Assert.That(KernelHarness.CheckSais(Faulty.States(s => s with { Heads = [.. s.Heads.Select((c, i) => i == 1 ? c + 1 : c)] }), Fib), Does.Contain("bucket heads"));
        Assert.That(KernelHarness.CheckSais(Faulty.States(s => s with { Tails = [.. s.Tails.Select((c, i) => i == 1 ? c - 1 : c)] }), Fib), Does.Contain("bucket tails"));
        Assert.That(KernelHarness.CheckSais(Faulty.States(s => s with { LmsPositions = [.. s.LmsPositions.Reverse()] }), Fib), Does.Contain("LMS positions"));
        Assert.That(KernelHarness.CheckSais(Faulty.States(s => s with { Reduced = [.. s.Reduced.Select((v, i) => i == 0 ? v + 1 : v)] }), Fib), Does.Contain("reduced string"));
        Assert.That(KernelHarness.CheckSais(Faulty.States(s => s with { NameCount = s.NameCount + 1 }), Fib), Does.Contain("names"));
        Assert.That(KernelHarness.CheckSais(Faulty.States(s => s with { LmsSorted = [.. s.LmsSorted.Reverse()] }), Fib), Does.Contain("LMS order"));
        Assert.That(KernelHarness.CheckSais(Faulty.States(s => s with { LmsSorted = [.. s.LmsSorted.Skip(1), -1] }), Fib), Does.Contain("LMS order"));
    }

    [Test]
    public void AnEqualNameSwap_InTheLmsOrder_IsAccepted_ButANameChangingSwapIsNot()
    {
        // LMS substrings with the same name may legitimately come out of the induction in either order
        string text = "abababab" + "cabababab";
        PreparedText prepared = TextPreparation.Prepare(text);
        SaisLevelState l0 = SaisStateOracle.Compute(KernelHarness.Symbols(prepared), prepared.AlphabetSize)[0];
        int[] nameOf = SaisStateOracle.NameOfPosition(l0);
        int a = Array.FindIndex(l0.LmsSorted, p => nameOf[p] == nameOf[l0.LmsSorted[1]] && p != l0.LmsSorted[1]);
        Assert.That(a, Is.GreaterThanOrEqualTo(0), "the text must contain two LMS substrings with the same name");
        Assert.That(KernelHarness.CheckSais(Faulty.States(s => s with { LmsSorted = SwapIndices(s.LmsSorted, 1, a) }), text), Is.Null);
        Assert.That(KernelHarness.CheckSais(Faulty.States(s => s with { LmsSorted = SwapIndices(s.LmsSorted, 0, s.LmsSorted.Length - 1) }), text), Does.Contain("LMS order"));
    }

    [Test]
    public void AnOutOfBoundsWrite_IsReported_AsAMemoryViolation_InEveryBuildMode()
    {
        string expected = EsaBuild.BoundsChecked ? "OutOfBounds" : "CanaryOverwritten";
        foreach (int offset in new[] { -1, 0, 5 })
        {
            string? failure = KernelHarness.CheckSais(Faulty.SaisWritingOutsideItsWorkspace(offset), Banana);
            Assert.That(failure, Does.StartWith("memory violation").And.Contain(expected), $"offset {offset}");
        }
        string? lcp = KernelHarness.CheckLcp(Faulty.LcpWritingOutsideItsWorkspace(), Banana);
        Assert.That(lcp, Does.StartWith("memory violation").And.Contain(expected));
    }

    [Test]
    public void AReadOfAnUninitialisedCell_IsReported_BecauseThePoisonPatternChangesTheResult()
    {
        // correct only when the uninitialised cell happens to hold the first poison pattern (0xCD...), wrong with the others
        string? failure = KernelHarness.CheckSais(Faulty.SaisReadingUninitialisedMemory(), Fib);
        Assert.That(failure, Does.StartWith("poison 0x00").And.Contain("SA["));
        failure = KernelHarness.CheckLcp(Faulty.LcpReadingUninitialisedMemory(), Fib);
        Assert.That(failure, Does.StartWith("poison 0x00").And.Contain("LCP["));
    }

    [Test]
    public void AWrongLcp_AWrongKind_AWrongEscapeCount_AWrongPlcp_AreReported()
    {
        Assert.That(KernelHarness.CheckLcp(Faulty.Lcp((lcp, plcp) => (WithChange(lcp, 5, 99), plcp)), Fib), Does.Contain("LCP["));
        Assert.That(KernelHarness.CheckLcp(Faulty.Lcp((lcp, plcp) => (CompactLcp.FromDense(lcp.ToArray(), LcpKind.U32), plcp)), Fib), Does.Contain("kind"));
        Assert.That(KernelHarness.CheckLcp(Faulty.Lcp((lcp, plcp) => (lcp, WithChange(plcp, 3, 1))), Fib), Does.Contain("PLCP["));
        string many = StringFamilies.AllSame(600); // 345 values >= 255: the rule picks U16 (U8 + overflow would be larger)
        Assert.That(KernelHarness.CheckLcp(Faulty.Lcp((lcp, plcp) => (CompactLcp.FromDense(lcp.ToArray(), LcpKind.U8), plcp)), many), Does.Contain("kind"));
        Assert.That(KernelHarness.CheckLcp(Faulty.Lcp((lcp, plcp) => (lcp, plcp)), many), Is.Null);
    }

    [Test]
    public void AnEngineThatThrows_IsReported()
    {
        string? failure = KernelHarness.CheckSais(new SaisEngine("thrower", (_, _) => throw new InvalidOperationException("boom")), Banana);
        Assert.That(failure, Does.StartWith("exception InvalidOperationException").And.Contain("boom"));
    }

    [Test]
    public void ThePipelineAdapter_RunsEngines_ForTheFuzzHarnessAndTheVectorTests()
    {
        var (sa, lcp) = KernelRegistry.Pipeline(KernelRegistry.ReferenceSais, KernelRegistry.ReferenceLcp)(Banana);
        Assert.That(sa, Is.EqualTo(NaiveSuffixArray.Build(Banana)));
        Assert.That(lcp.ToArray(), Is.EqualTo(NaiveLcp.Build(Banana, sa)));
    }

    private static int[] SwapIndices(int[] values, int i, int j)
    {
        int[] copy = [.. values];
        (copy[i], copy[j]) = (copy[j], copy[i]);
        return copy;
    }

    private static CompactLcp WithChange(CompactLcp lcp, int index, int value)
    {
        int[] dense = lcp.ToArray();
        dense[index] = value;
        return CompactLcp.FromDense(dense);
    }

    private static int[] WithChange(int[]? plcp, int index, int value)
    {
        int[] copy = [.. plcp!];
        copy[index] = value;
        return copy;
    }

    /// <summary>Test-only engines with a deliberate defect each (never part of the registry).</summary>
    private static class Faulty
    {
        internal static SaisEngine Sais(Action<int[]> corrupt) => new("faulty-sa", (p, a) =>
        {
            SaisRun r = KernelRegistry.RunReferenceSais(p, a);
            corrupt(r.Sa);
            return r;
        });

        internal static SaisEngine SaisTrace(Func<SaisLevel, SaisLevel>? change) => new("faulty-trace", (p, a) =>
        {
            SaisRun r = KernelRegistry.RunReferenceSais(p, a);
            return change is null ? r with { Levels = r.Levels.Take(r.Levels.Count - 1).ToList() } : r with { Levels = r.Levels.Select(change).ToList() };
        });

        internal static SaisEngine States(Func<SaisLevelState, SaisLevelState> change) => new("faulty-states", (p, a) =>
        {
            SaisRun r = KernelRegistry.RunReferenceSais(p, a);
            List<SaisLevelState> states = SaisStateOracle.Compute(KernelHarness.Symbols(p), p.AlphabetSize);
            states[0] = change(states[0]);
            return r with { States = states };
        });

        internal static SaisEngine SaisWritingOutsideItsWorkspace(int offset) => new("overflow", (p, a) =>
        {
            SaisRun r = KernelRegistry.RunReferenceSais(p, a);
            Buf<int> work = a.Allocate<int>(p.Length + 1, zero: true, "workspace");
            work[offset < 0 ? offset : p.Length + 1 + offset] = 1; // one int before the block / `offset` ints after its end
            return r;
        });

        internal static LcpEngine LcpWritingOutsideItsWorkspace() => new("overflow-lcp", (p, sa, a) =>
        {
            LcpRun r = KernelRegistry.RunReferenceLcp(p, sa, a);
            Buf<int> phi = a.Allocate<int>(p.Length, zero: true, "phi");
            phi[p.Length] = 1; // the classic off-by-one past the end of Phi/PLCP
            return r;
        });

        internal static SaisEngine SaisReadingUninitialisedMemory() => new("uninitialised", (p, a) =>
        {
            SaisRun r = KernelRegistry.RunReferenceSais(p, a);
            Buf<int> scratch = a.Allocate<int>(4, zero: false, "scratch");
            if (scratch[0] != unchecked((int)0xCDCDCDCD)) (r.Sa[1], r.Sa[2]) = (r.Sa[2], r.Sa[1]);
            return r;
        });

        internal static LcpEngine LcpReadingUninitialisedMemory() => new("uninitialised-lcp", (p, sa, a) =>
        {
            LcpRun r = KernelRegistry.RunReferenceLcp(p, sa, a);
            Buf<int> scratch = a.Allocate<int>(4, zero: false, "scratch");
            if (scratch[0] == unchecked((int)0xCDCDCDCD)) return r;
            int[] dense = r.Lcp.ToArray();
            dense[4] += 1;
            return new LcpRun(CompactLcp.FromDense(dense), r.Plcp);
        });

        internal static LcpEngine Lcp(Func<CompactLcp, int[]?, (CompactLcp Lcp, int[]? Plcp)> change) => new("faulty-lcp", (p, sa, a) =>
        {
            LcpRun r = KernelRegistry.RunReferenceLcp(p, sa, a);
            var (lcp, plcp) = change(r.Lcp, r.Plcp);
            return new LcpRun(lcp, plcp);
        });
    }
}
