using SuffixTree.Esa.Tests.TestSupport;

namespace SuffixTree.Esa.Tests.Kernels;

/// <summary>
/// 1o.1b: the SA-IS pointer kernel (<see cref="SaisKernel"/>) beyond the generic registry run: tiny texts, argument checks, the selector and the public entry
/// points, the trace phases, "no recursion when all names are distinct" (review MIN-3, M08), the LMS-substring predicate of the kernel's naming against the
/// definition on the same enumeration as the reference (review T2: M10, M11), and the layout invariants at the limits of the name slots.
/// </summary>
[TestFixture]
public class SaisKernelTests
{
    private static SaisEngine Engine => KernelRegistry.PointerSais;

    [Test]
    public void TinyTexts_n0_n1_n2_n3_AndKnownExamples()
    {
        foreach (string text in new[] { "", "a", "aa", "ab", "ba", "aaa", "aba", "abc", "cab", "banana", "mississippi", "abracadabra", "aabaabaab" })
        {
            KernelHarness.AssertSais(Engine, text);
            PreparedText p = TextPreparation.Prepare(text);
            Assert.That(SaisKernel.Build(p), Is.EqualTo(NaiveSuffixArray.Build(text)), $"public entry on '{text}'");
        }
    }

    [Test]
    public void TheSelectorAndTheEntryPoints_GiveTheReferenceResult_AndTheTraceIsFilled()
    {
        string text = StringFamilies.RandomDna(20_000, 5);
        PreparedText p = TextPreparation.Prepare(text);
        int[] reference = SaisKernel.Build(SaisKernelKind.Reference, p);
        Assert.That(SaisKernel.Build(SaisKernelKind.Pointer, p), Is.EqualTo(reference));
        using var arena = new NativeArena(ArenaOptions.Guarded());
        Assert.That(SaisKernel.Build(p, arena), Is.EqualTo(reference));
        arena.VerifyCanaries();

        var referenceTrace = new SaisTrace();
        var kernelTrace = new SaisTrace();
        SaisKernel.Build(SaisKernelKind.Reference, p, referenceTrace);
        SaisKernel.Build(SaisKernelKind.Pointer, p, kernelTrace);
        Assert.That(kernelTrace.Levels, Is.EqualTo(referenceTrace.Levels));
        foreach (SaisPhase phase in Enum.GetValues<SaisPhase>())
            Assert.That(kernelTrace.TotalTicks(phase) > 0, Is.EqualTo(referenceTrace.TotalTicks(phase) > 0), $"phase {phase} is timed by both or by neither");
    }

    [Test]
    public void ArgumentChecks_AreTheOnesOfTheReference()
    {
        using var arena = new NativeArena(ArenaOptions.Guarded());
        Buf<byte> text = arena.Allocate<byte>(3, zero: true, "text");
        Buf<int> sa = arena.Allocate<int>(4, zero: false, "sa");
        text[0] = 0;
        text[1] = 2; // outside [0, 2)
        text[2] = 1;
        Assert.Throws<ArgumentException>(() => SaisKernel.Build(text, 3, 2, sa, arena, default(NoProbe)));
        Assert.Throws<ArgumentOutOfRangeException>(() => SaisKernel.Build(text, 3, -1, sa, arena, default(NoProbe)));
        Assert.Throws<ArgumentOutOfRangeException>(() => SaisKernel.Build(text, -1, 2, sa, arena, default(NoProbe)));
        Assert.Throws<NotSupportedException>(() => SaisKernel.Build(text, int.MaxValue, 2, sa, arena, default(NoProbe)));
        arena.VerifyCanaries();
    }

    [Test]
    public void AllLmsSubstringsDistinct_NeedNoRecursion()
    {
        // review MIN-3 (mutant M08): with names == m the reduced SA is placed directly; recursing anyway gives the same SA but one useless level
        string text = StringFamilies.AllDistinct(2000, 1);
        PreparedText p = TextPreparation.Prepare(text);
        var trace = new SaisTrace();
        SaisKernel.Build(SaisKernelKind.Pointer, p, trace);
        Assert.That(trace.Depth, Is.EqualTo(1));
    }

    [Test]
    public void TheKernelsLmsSubstringPredicate_MatchesTheDefinition_OnAllSmallStrings()
    {
        // review T2 (mutants M10, M11): the type comparison and the sentinel case are invisible to every SA-level test, so the predicate that the kernel's naming uses is tested on its own
        long pairs = 0, equalPairs = 0;
        foreach (var (sigma, maxLen) in new[] { (2, 14), (3, 9), (4, 7) })
            foreach (string s in ExhaustiveStrings.All(sigma, maxLen))
            {
                PreparedText p = TextPreparation.Prepare(s);
                ByteSymbols reader = p.ByteReader();
                bool[] isS = SaisReference.ClassifyTypes(reader);
                byte[] types = [.. isS.Select(b => b ? (byte)1 : (byte)0)];
                int n = s.Length;
                var lms = new List<int>();
                for (int i = 1; i < n; i++)
                    if (SaisReference.IsLms(isS, i)) lms.Add(i);
                lms.Add(n);
                for (int x = 0; x + 1 < lms.Count; x++)
                    for (int y = 0; y + 1 < lms.Count; y++)
                    {
                        if (x == y) continue;
                        bool expected = Definition(s, isS, lms, x, y);
                        pairs++;
                        if (expected) equalPairs++;
                        Assert.That(SaisSafePasses.LmsSubstringsEqual(reader, types, lms[x], lms[y]), Is.EqualTo(expected), $"'{s}': LMS substrings at {lms[x]} and {lms[y]}");
                    }
            }
        Assert.That(equalPairs, Is.GreaterThan(10_000));
        TestContext.Out.WriteLine($"{pairs} ordered LMS pairs checked against the definition, {equalPairs} equal");
    }

    [Test]
    public void TheKernelsTypesAndLmsCount_EqualTheReferences_OnTheDefaultCorpus()
    {
        // the classification is fused with the LMS count and the types are bytes: the states of the harness compare m per level; this pins the type array itself
        foreach (var (name, text) in KernelCorpus.Default().Where(c => c.Text.Length is > 0 and <= 2000))
        {
            PreparedText p = TextPreparation.Prepare(text);
            bool[] expected = p.Width == SymbolWidth.Byte ? SaisReference.ClassifyTypes(p.ByteReader()) : SaisReference.ClassifyTypes(p.UShortReader());
            using var arena = new NativeArena(ArenaOptions.Guarded());
            Buf<int> sa = arena.Allocate<int>(p.Length + 1, zero: false, "sa");
            var probe = new SaisStateRecorder();
            if (p.Width == SymbolWidth.Byte)
            {
                Buf<byte> t = arena.Allocate<byte>(p.Length, zero: false, "text");
                p.Bytes.CopyTo(t.AsSpan(p.Length));
                SaisKernel.Build(t, p.Length, p.AlphabetSize, sa, arena, new RecordingSaisProbe(probe));
            }
            else
            {
                Buf<ushort> t = arena.Allocate<ushort>(p.Length, zero: false, "text");
                p.UShorts.CopyTo(t.AsSpan(p.Length));
                SaisKernel.Build(t, p.Length, p.AlphabetSize, sa, arena, new RecordingSaisProbe(probe));
            }
            int lmsExpected = SaisReference.CountLms(expected);
            Assert.That(probe.ToStates()[0].LmsPositions, Has.Length.EqualTo(lmsExpected), name);
        }
    }

    private static bool Definition(string s, bool[] isS, List<int> lms, int x, int y)
    {
        int ax = lms[x], bx = lms[x + 1], ay = lms[y], by = lms[y + 1];
        if (bx == s.Length || by == s.Length) return false;
        if (bx - ax != by - ay) return false;
        for (int d = 0; d <= bx - ax; d++)
            if (s[ax + d] != s[ay + d] || isS[ax + d] != isS[ay + d]) return false;
        return true;
    }
}
