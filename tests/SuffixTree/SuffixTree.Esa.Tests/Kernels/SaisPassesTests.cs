using SuffixTree.Esa.Tests.TestSupport;

namespace SuffixTree.Esa.Tests.Kernels;

/// <summary>
/// 1o.1c: the induction passes (L and S, both rounds) and the naming of the SA-IS kernel as pointer loops, compared with their definitions in situ.
/// Three nets, strongest first:
/// <list type="number">
/// <item>every pass of every level of the real pipeline is replayed on the state that the kernel itself reported before it (the <see cref="ISaisProbe.SaState"/> snapshots)
/// with the verbatim safe copy <see cref="SaisSafePasses"/>, and the whole suffix array (all n + 1 entries, ties and empty slots included) must be bit-identical to the
/// state the pointer pass produced; this includes the recursion levels, whose texts are the reduced strings of the previous level;</item>
/// <item>the states of round 2 are compared with their DEFINITION (the true suffix array of the level: sorted LMS suffixes at the ends of their buckets; L-type suffixes
/// at their final indices after the L pass; the whole suffix array after the S pass);</item>
/// <item>the predicate of the naming against the definition of equal LMS substrings (kernel T2, in <see cref="SaisKernelTests"/>).</item>
/// </list>
/// The passes of the kernel write a dummy value to <c>SA[0]</c> when they do not induce (no branch), so the snapshots also pin that <c>SA[0] = n</c> is restored.
/// </summary>
[TestFixture]
public class SaisPassesTests
{
    [Test]
    public void EveryInductionAndNamingPass_EqualsTheSafeCopy_StateForState_OnTheDefaultCorpus()
    {
        int texts = 0;
        long passes = 0;
        foreach (var (name, text) in KernelCorpus.Default())
        {
            if (text.Length == 0) continue;
            texts++;
            passes += ReplayPasses(name, text);
        }
        Assert.That(texts, Is.GreaterThan(7000));
        TestContext.Out.WriteLine($"{texts} texts, {passes} pass replays (levels x 5 passes), all bit-identical");
    }

    [TestCase("banana")]
    [TestCase("mississippi")]
    [TestCase("abracadabra")]
    [TestCase("aaaaaaaaaaaaaaaa")]
    [TestCase("abababababababab")]
    [TestCase("abcabcabcabcabc")]
    public void KnownTexts_ReplayAndDefinition(string text)
    {
        Assert.That(ReplayPasses(text, text), Is.GreaterThanOrEqualTo(5));
    }

    [Test]
    public void TheStatesOfThePasses_AreTheirDefinition_OnTheDefaultCorpus()
    {
        int levelsChecked = 0;
        foreach (var (name, text) in KernelCorpus.Default())
        {
            if (text.Length == 0) continue;
            PreparedText prepared = TextPreparation.Prepare(text);
            SaisRun run = RunWithStates(prepared);
            int[] levelText = KernelHarness.Symbols(prepared);
            List<SaisLevelState> oracle = SaisStateOracle.Compute(levelText, prepared.AlphabetSize);
            for (int level = 0; level < run.States!.Count; level++)
            {
                SaisLevelState st = run.States[level];
                if (levelText.Length <= 2000 && CheckDefinition(name, level, levelText, st, oracle[level])) levelsChecked++;
                levelText = st.Reduced;
            }
        }
        Assert.That(levelsChecked, Is.GreaterThan(7000));
    }

    [Test]
    public void TheSnapshotsAreReportedForEveryPass()
    {
        PreparedText prepared = TextPreparation.Prepare("mississippi");
        SaisRun run = RunWithStates(prepared);
        SaisLevelState level0 = run.States![0];
        Assert.That(level0.SaStates, Is.Not.Null);
        foreach (SaisPhase phase in new[] { SaisPhase.Round1Place, SaisPhase.Round1InduceL, SaisPhase.Round1InduceS, SaisPhase.CompactLms, SaisPhase.Name, SaisPhase.Round2Place, SaisPhase.Round2InduceL, SaisPhase.Round2InduceS })
            Assert.That(level0.SaStates!.ContainsKey(phase), Is.True, phase.ToString());
        Assert.That(level0.SaStates![SaisPhase.Round2InduceS], Is.EqualTo(NaiveSuffixArray.Build("mississippi")));
        Assert.That(level0.SaStates[SaisPhase.Round1InduceL][0], Is.EqualTo(11), "SA[0] = n after the L pass (the dummy writes are undone)");
    }

    [Test]
    public void TheReplay_ReportsATamperedSnapshot()
    {
        // fault injection: a state that differs in one entry in any of the five passes (or in SA[0]) must be named
        foreach (SaisPhase phase in new[] { SaisPhase.Round1InduceL, SaisPhase.Round1InduceS, SaisPhase.Name, SaisPhase.Round2InduceL, SaisPhase.Round2InduceS })
            foreach (int index in new[] { 0, 3, 7 })
            {
                var ex = Assert.Throws<AssertionException>(() => ReplayPasses("tampered", "mississippi", snap =>
                {
                    int[] copy = (int[])snap[phase].Clone();
                    copy[index] ^= 1;
                    snap[phase] = copy;
                }), $"{phase}[{index}]");
                Assert.That(ex!.Message, Does.Contain($"SA[{index}]"));
            }
    }

    private static SaisRun RunWithStates(PreparedText prepared)
    {
        using var arena = new NativeArena(ArenaOptions.Guarded());
        SaisRun run = KernelRegistry.PointerSaisWithStates.Run(prepared, arena);
        arena.VerifyCanaries();
        return run;
    }

    /// <summary>Replays the five passes of every level on the kernel's own snapshots; returns the number of replays.</summary>
    private static int ReplayPasses(string name, string text, Action<Dictionary<SaisPhase, int[]>>? tamper = null)
    {
        PreparedText prepared = TextPreparation.Prepare(text);
        SaisRun run = RunWithStates(prepared);
        int[] levelText = KernelHarness.Symbols(prepared);
        int replays = 0;
        for (int level = 0; level < run.States!.Count; level++)
        {
            SaisLevelState st = run.States[level];
            IReadOnlyDictionary<SaisPhase, int[]> snap = st.SaStates ?? throw new InvalidOperationException("the probe reports no snapshots");
            if (tamper is not null && level == 0)
            {
                var copy = new Dictionary<SaisPhase, int[]>(snap);
                tamper(copy);
                snap = copy;
            }
            byte[] types = Types(levelText);
            string where = $"'{name}' level {level} (n = {levelText.Length})";

            AssertSame(where + " round 1 induce L", snap[SaisPhase.Round1InduceL], SafeInduceL(levelText, types, st.Heads, snap[SaisPhase.Round1Place]));
            AssertSame(where + " round 1 induce S", snap[SaisPhase.Round1InduceS], SafeInduceS(levelText, types, st.Tails, snap[SaisPhase.Round1InduceL]));
            var (afterName, names) = SafeName(levelText, types, snap[SaisPhase.CompactLms], st.LmsPositions.Length);
            AssertSame(where + " naming", snap[SaisPhase.Name], afterName);
            Assert.That(names, Is.EqualTo(st.NameCount), where + " number of names");
            AssertSame(where + " round 2 induce L", snap[SaisPhase.Round2InduceL], SafeInduceL(levelText, types, st.Heads, snap[SaisPhase.Round2Place]));
            AssertSame(where + " round 2 induce S", snap[SaisPhase.Round2InduceS], SafeInduceS(levelText, types, st.Tails, snap[SaisPhase.Round2InduceL]));
            replays += 5;
            levelText = st.Reduced;
        }
        return replays;
    }

    private static bool CheckDefinition(string name, int level, int[] text, SaisLevelState st, SaisLevelState oracle)
    {
        int n = text.Length;
        byte[] types = Types(text);
        int[] sa = TrueSuffixArray(text, st.Sigma);
        IReadOnlyDictionary<SaisPhase, int[]> snap = st.SaStates!;
        string where = $"'{name}' level {level} (n = {n})";

        // round 1, after placing the LMS positions: each at the end of its bucket, ascending by position inside the bucket, everything else empty, SA[0] = n
        var placed1 = new int[n + 1];
        Array.Fill(placed1, -1);
        placed1[0] = n;
        int[] fill1 = (int[])st.Tails.Clone();
        for (int i = n - 1; i >= 1; i--)
            if (types[i] == 1 && types[i - 1] == 0) placed1[--fill1[text[i]]] = i;
        AssertSame(where + " definition: round 1 after placing the LMS positions", snap[SaisPhase.Round1Place], placed1);

        // naming: SA[0] = n, the sorted LMS positions in SA[1..m], the name of the LMS position p at SA[1 + m + (p >> 1)], every other entry -1
        int m = st.LmsPositions.Length;
        int[] afterName = new int[n + 1];
        Array.Fill(afterName, -1);
        afterName[0] = n;
        int[] compacted = snap[SaisPhase.CompactLms];
        for (int i = 1; i <= m; i++) afterName[i] = compacted[i];
        int[] nameOf = SaisStateOracle.NameOfPosition(oracle);
        foreach (int p in oracle.LmsPositions) afterName[1 + m + (p >> 1)] = nameOf[p];
        AssertSame(where + " definition: after naming", snap[SaisPhase.Name], afterName);

        // the suffix array of the reduced string, as it stands after the recursion or the direct placement (n + 1 layout)
        AssertSame(where + " definition: suffix array of the reduced string", snap[SaisPhase.Reduced], TrueSuffixArray(st.Reduced, Math.Max(st.NameCount, 1)));

        // sorted LMS suffixes at the ends of their buckets, everything else empty, SA[0] = n
        var placed = new int[n + 1];
        Array.Fill(placed, -1);
        placed[0] = n;
        int[] fill = (int[])st.Tails.Clone();
        for (int i = n; i >= 1; i--)
        {
            int p = sa[i];
            if (types[p] == 1 && p > 0 && types[p - 1] == 0) placed[--fill[text[p]]] = p;
        }
        AssertSame(where + " definition: after placing the sorted LMS suffixes", snap[SaisPhase.Round2Place], placed);

        // the L-type suffixes sit at their final indices, the LMS suffixes still at the ends of the buckets
        int[] afterL = (int[])placed.Clone();
        for (int i = 1; i <= n; i++)
            if (types[sa[i]] == 0) afterL[i] = sa[i];
        AssertSame(where + " definition: after the L pass", snap[SaisPhase.Round2InduceL], afterL);

        AssertSame(where + " definition: after the S pass (the suffix array)", snap[SaisPhase.Round2InduceS], sa);
        return true;
    }

    private static int[] TrueSuffixArray(int[] text, int sigma)
    {
        if (text.Length <= 64)
        {
            // fully independent: sort the suffixes by comparing them
            var idx = Enumerable.Range(0, text.Length).ToList();
            idx.Sort((a, b) =>
            {
                int k = 0;
                while (a + k < text.Length && b + k < text.Length && text[a + k] == text[b + k]) k++;
                if (a + k == text.Length) return b + k == text.Length ? 0 : -1;
                if (b + k == text.Length) return 1;
                return text[a + k].CompareTo(text[b + k]);
            });
            return [text.Length, .. idx];
        }
        return SaisReference.Build(new IntSymbols(text), sigma);
    }

    /// <summary>The types by the definition: 1 = S-type, position n (the virtual sentinel) is S-type.</summary>
    private static byte[] Types(int[] t)
    {
        int n = t.Length;
        var isS = new byte[n + 1];
        isS[n] = 1;
        // position n - 1 is L-type (it is larger than the sentinel); below it the definition: smaller than the successor, or equal and the successor is S-type
        for (int i = n - 2; i >= 0; i--)
            isS[i] = t[i] < t[i + 1] || (t[i] == t[i + 1] && isS[i + 1] == 1) ? (byte)1 : (byte)0;
        return isS;
    }

    private static int[] SafeInduceL(int[] text, byte[] types, int[] heads, int[] before)
    {
        int[] sa = (int[])before.Clone();
        SaisSafePasses.InduceL(new IntSymbols(text), sa, types, (int[])heads.Clone());
        return sa;
    }

    private static int[] SafeInduceS(int[] text, byte[] types, int[] tails, int[] before)
    {
        int[] sa = (int[])before.Clone();
        SaisSafePasses.InduceS(new IntSymbols(text), sa, types, (int[])tails.Clone());
        return sa;
    }

    private static (int[] Sa, int Names) SafeName(int[] text, byte[] types, int[] before, int m)
    {
        int[] sa = (int[])before.Clone();
        int names = SaisSafePasses.NameLmsSubstrings(new IntSymbols(text), sa, types, m);
        return (sa, names);
    }

    private static void AssertSame(string what, int[] actual, int[] expected)
    {
        int k = Math.Min(actual.Length, expected.Length);
        for (int i = 0; i < k; i++)
            if (actual[i] != expected[i])
                Assert.Fail($"{what}: SA[{i}] = {actual[i]}, expected {expected[i]}");
        Assert.That(actual.Length, Is.EqualTo(expected.Length), what + ": length");
    }
}
