using SuffixTree.Esa.Tests.TestSupport;

namespace SuffixTree.Esa.Tests.Fuzz;

/// <summary>
/// Default-run tests of the differential fuzz harness itself (work package 1a.7): determinism of the generators, case ids, the minimizer, fault
/// detection (a deliberately broken engine must be found and shrunk), the time budget, the regression corpus, and a small smoke run of every
/// tier. The acceptance-size runs are <see cref="DifferentialFuzzTests"/> (category <c>Fuzz</c>, explicit).
/// </summary>
[TestFixture]
public class FuzzHarnessTests
{
    // ---- a deliberately broken engine: the real pipeline plus a fault triggered by a property of the text ----

    private static EsaBuilder SaFault(Func<string, bool> trigger) => text =>
    {
        var (sa, lcp) = EsaPipeline.Build(text);
        if (sa.Length >= 3 && trigger(text)) (sa[1], sa[2]) = (sa[2], sa[1]);
        return (sa, lcp);
    };

    private static EsaBuilder LcpFault(Func<string, bool> trigger) => text =>
    {
        var (sa, lcp) = EsaPipeline.Build(text);
        if (text.Length >= 1 && trigger(text))
        {
            int[] values = lcp.ToArray();
            values[text.Length]++;
            lcp = CompactLcp.FromDense(values);
        }
        return (sa, lcp);
    };

    private static bool ContainsOrdinal(string text, string part) => text.Contains(part, StringComparison.Ordinal);

    [Test]
    public void CaseId_RoundTrips_AndRejectsGarbage()
    {
        var c = new FuzzCase("block-copy", 0xDEADBEEFCAFEF00DUL, 1234, 27, 99);
        Assert.That(c.Id, Is.EqualTo("block-copy:deadbeefcafef00d:1234:27:99"));
        Assert.That(FuzzCase.Parse(c.Id), Is.EqualTo(c));
        foreach (string bad in new[] { "", "a:b", "x:zz:1:1:1", "x:1:-1:1:1", "x:1:1:1", ":1:1:1:1", "x:1:1:1:1:1" })
            Assert.Throws<FormatException>(() => FuzzCase.Parse(bad), bad);
    }

    [Test]
    public void EveryFamily_IsDeterministic_HasExactLength_AndStaysInsideItsAlphabet()
    {
        foreach (string family in FuzzTexts.Families)
        {
            foreach (int n in new[] { 0, 1, 2, 3, 7, 64, 257, 1000 })
            {
                foreach (int sigma in new[] { 1, 2, 5, 27, 256, 257, 1000 })
                {
                    int s = family == "unicode-mix" ? Math.Min(sigma, 10) : family == "permutation" ? n : sigma;
                    var c = new FuzzCase(family, 12345UL + (ulong)n, n, s, 7);
                    string a = c.Text(), b = c.Text();
                    Assert.That(a, Is.EqualTo(b), c.Id + " is not deterministic");
                    Assert.That(a, Has.Length.EqualTo(n), c.Id);
                    if (family != "permutation") Assert.That(a.Distinct().Count(), Is.LessThanOrEqualTo(s), c.Id + " exceeds its alphabet");
                }
            }
        }
        Assert.Throws<ArgumentException>(() => new FuzzCase("nope", 1, 5, 2, 0).Text());
    }

    [Test]
    public void Pick_CoversEveryFamily_BothSymbolWidths_AndEveryTierProducesTextsOfTheRequestedLength()
    {
        var rng = new SplitMix64(99);
        var families = new HashSet<string>();
        bool narrow = false, wide = false;
        foreach (FuzzTier tier in new[] { FuzzTier.Tiny, FuzzTier.Small, FuzzTier.Large })
        {
            for (long i = 0; i < 400; i++)
            {
                var c = FuzzRunner.CaseAt(99, tier, i);
                Assert.That(c.N, Is.InRange(tier.MinN, tier.MaxN), c.Id);
                if (tier == FuzzTier.Large && i > 25) continue; // generating 10^5 units 400 times is not needed to check the picker
                string t = c.Text();
                Assert.That(t, Has.Length.EqualTo(c.N), c.Id);
                families.Add(c.Family);
                if (c.Sigma > 256) wide = true; else narrow = true;
            }
        }
        for (int i = 0; i < 20000; i++) families.Add(FuzzTexts.Pick(rng, 100).Family);
        Assert.That(families, Is.EquivalentTo(FuzzTexts.Families));
        Assert.That(narrow && wide, Is.True);
    }

    [Test]
    public void CaseAt_IsAPureFunction_AndConsecutiveCasesAreNotShiftedCopiesOfOneStream()
    {
        var ids = new HashSet<string>();
        var seeds = new HashSet<ulong>();
        for (long i = 0; i < 2000; i++)
        {
            var c = FuzzRunner.CaseAt(7, FuzzTier.Small, i);
            Assert.That(FuzzRunner.CaseAt(7, FuzzTier.Small, i), Is.EqualTo(c));
            ids.Add(c.Id);
            seeds.Add(c.Seed);
        }
        Assert.That(ids, Has.Count.EqualTo(2000));
        Assert.That(seeds, Has.Count.EqualTo(2000));
        Assert.That(FuzzRunner.CaseAt(8, FuzzTier.Small, 0), Is.Not.EqualTo(FuzzRunner.CaseAt(7, FuzzTier.Small, 0)));
        Assert.That(FuzzRunner.CaseAt(7, FuzzTier.Tiny, 0), Is.Not.EqualTo(FuzzRunner.CaseAt(7, FuzzTier.Small, 0)));
    }

    [Test]
    public void Minimizer_ShrinksToTheSmallestFailingText()
    {
        // fails iff the text contains "aba": the unique minimal failing text is "aba"
        string big = new string('c', 300) + "xyz" + new string('b', 200) + "aba" + new string('q', 400);
        string min = FuzzMinimizer.Minimize(big, t => ContainsOrdinal(t, "aba"));
        Assert.That(min, Is.EqualTo("aba"));

        // symbols are merged and renamed: fails iff some symbol occurs three times in a row
        string runs = "qwertyuuuiopasdf";
        min = FuzzMinimizer.Minimize(runs, t => t.Length >= 3 && Enumerable.Range(0, t.Length - 2).Any(i => t[i] == t[i + 1] && t[i + 1] == t[i + 2]));
        Assert.That(min, Is.EqualTo("aaa"));

        // 1-minimality on a harder predicate: needs a 'b' before an 'a' and a 'c' somewhere after both
        string text = "dddbdddadddcddd";
        string m2 = FuzzMinimizer.Minimize(text, t => t.IndexOf('b') is var b && b >= 0 && t.IndexOf('a', b + 1) is var a && a > b && t.IndexOf('c', a + 1) > a);
        Assert.That(m2, Is.EqualTo("bac"), "the predicate depends on the literal letters b, a, c in this order, so only the d units can go");
    }

    [Test]
    public void Minimizer_ResultIsOneMinimal_ForAnInjectedEngineFault()
    {
        EsaBuilder faulty = SaFault(t => ContainsOrdinal(t, "abab"));
        string start = string.Concat(Enumerable.Repeat("abcabab", 40));
        Assert.That(FuzzOracle.Check(start, faulty), Is.Not.Null);
        string min = FuzzMinimizer.Minimize(start, t => FuzzOracle.Check(t, faulty) is { Op: "sa-naive" });
        Assert.That(min, Is.EqualTo("abab"));
        Assert.That(FuzzOracle.Check(min, faulty), Is.Not.Null);
        for (int i = 0; i < min.Length; i++)
            Assert.That(FuzzOracle.Check(min.Remove(i, 1), faulty), Is.Null, $"deleting unit {i} of '{min}' keeps the failure: not minimal");
    }

    [Test]
    public void Runner_FindsAnInjectedSaFault_RecordsTheCase_AndShrinksIt()
    {
        EsaBuilder faulty = SaFault(t => ContainsOrdinal(t, "aba"));
        var summary = FuzzRunner.Run(new FuzzOptions { Tier = FuzzTier.Tiny, Cases = 50_000, Candidate = faulty, MaxFailures = 3 });
        Assert.That(summary.Failures, Is.Not.Empty, summary.ToString());
        foreach (var f in summary.Failures)
        {
            Assert.That(f.Op, Is.EqualTo("sa-naive"));
            Assert.That(f.MinimizedText, Is.EqualTo("aba"));
            string original = f.Case.Text();
            Assert.That(ContainsOrdinal(original, "aba"), Is.True, f.Case.Id);
            Assert.That(FuzzOracle.Check(original, faulty), Is.Not.Null, "the recorded id must reproduce the failure");
            Assert.That(FuzzOracle.Check(original, FuzzOracle.ReferencePipeline), Is.Null, "the real pipeline is fine on the same text");
            Assert.That(FuzzRegressions.Parse(FuzzRegressions.Format(f))[0], Is.EqualTo(new FuzzRegression(f.Case, f.Op, "aba")));
        }
    }

    [Test]
    public void Runner_FindsAnInjectedLcpFault_InTheNaiveAndInTheVerifierRange()
    {
        var small = FuzzRunner.Run(new FuzzOptions { Tier = FuzzTier.Small, Cases = 2_000, Candidate = LcpFault(t => ContainsOrdinal(t, "aaa")), MaxFailures = 1 });
        Assert.That(small.Failures, Has.Count.EqualTo(1), small.ToString());
        Assert.That(small.Failures[0].Op, Is.EqualTo("lcp-naive"));
        Assert.That(small.Failures[0].MinimizedText, Is.EqualTo("aaa"));

        var large = FuzzRunner.Run(new FuzzOptions { Tier = FuzzTier.Large, Cases = 40, Candidate = LcpFault(t => t.Length > 0), MaxFailures = 1 });
        Assert.That(large.Failures, Has.Count.EqualTo(1), large.ToString());
        Assert.That(large.Failures[0].Op, Is.EqualTo("lcp-kasai"));
        Assert.That(large.Failures[0].MinimizedText, Has.Length.EqualTo(1), "shrinking crosses the naive limit: the failure class (lcp) is kept, not the oracle");

        var saLarge = FuzzRunner.Run(new FuzzOptions { Tier = FuzzTier.Large, Cases = 40, Candidate = SaFault(t => ContainsOrdinal(t, "ab")), MaxFailures = 1 });
        Assert.That(saLarge.Failures, Has.Count.EqualTo(1), saLarge.ToString());
        Assert.That(saLarge.Failures[0].Op, Is.EqualTo("sa-verify"));
        Assert.That(saLarge.Failures[0].MinimizedText, Is.EqualTo("ab"));
    }

    [Test]
    public void Runner_ReportsAnEngineExceptionAsAFailure_AndComparesWithAReference()
    {
        EsaBuilder throwing = text => text.Contains("ab", StringComparison.Ordinal) ? throw new InvalidOperationException("boom") : EsaPipeline.Build(text);
        var summary = FuzzRunner.Run(new FuzzOptions { Tier = FuzzTier.Tiny, Cases = 20_000, Candidate = throwing, MaxFailures = 1 });
        Assert.That(summary.Failures, Has.Count.EqualTo(1));
        Assert.That(summary.Failures[0].Op, Is.EqualTo("exception:InvalidOperationException"));
        Assert.That(summary.Failures[0].MinimizedText, Is.EqualTo("ab"));

        Assert.That(FuzzOracle.Check("mississippi", FuzzOracle.ReferencePipeline, FuzzOracle.ReferencePipeline), Is.Null);
        EsaBuilder swapped = SaFault(_ => true);
        Assert.That(FuzzOracle.Check("mississippi", swapped, FuzzOracle.ReferencePipeline), Is.Not.Null);
    }

    [Test]
    public void Runner_StopsAtTheTimeBudget_AndAMatchingRunIsReproducible()
    {
        var limited = FuzzRunner.Run(new FuzzOptions { Tier = FuzzTier.Tiny, Cases = 1_000_000_000, TimeBudget = TimeSpan.FromMilliseconds(150) });
        Assert.That(limited.BudgetExhausted, Is.True);
        Assert.That(limited.Cases, Is.GreaterThan(0).And.LessThan(1_000_000_000));
        Assert.That(limited.Failures, Is.Empty);
        Assert.That(limited.Elapsed, Is.LessThan(TimeSpan.FromSeconds(30)));

        // the same seed and case count give the same sequence of cases: first/last ids and the per-family counts agree
        var a = FuzzRunner.Run(new FuzzOptions { MasterSeed = 5, Tier = FuzzTier.Tiny, Cases = 500 });
        var b = FuzzRunner.Run(new FuzzOptions { MasterSeed = 5, Tier = FuzzTier.Tiny, Cases = 500 });
        Assert.That((a.FirstCaseId, a.LastCaseId, a.TotalUnits), Is.EqualTo((b.FirstCaseId, b.LastCaseId, b.TotalUnits)));
        Assert.That(a.PerFamily, Is.EquivalentTo(b.PerFamily));
    }

    [Test]
    public void RegressionFormat_RoundTrips_NulSurrogatesAndTheEmptyText()
    {
        var c = new FuzzCase("unicode-mix", 77, 4, 3, 0);
        foreach (string text in new[] { "", "a", "\0￿", "𐀀x", "\uDFFF\uD800" })
        {
            string line = FuzzRegressions.Format(new FuzzFailure(c, "sa-naive", "d", text));
            var parsed = FuzzRegressions.Parse("# comment\n\n" + line + "\r\n");
            Assert.That(parsed, Has.Count.EqualTo(1));
            Assert.That(parsed[0], Is.EqualTo(new FuzzRegression(c, "sa-naive", text)));
        }
        Assert.Throws<FormatException>(() => FuzzRegressions.Parse("only-one-field"));
        Assert.Throws<FormatException>(() => FuzzRegressions.Parse("bad-id\top\t0061"));
        Assert.Throws<FormatException>(() => FuzzRegressions.Parse($"{c.Id}\top\tzzzz"));
    }

    [Test]
    public void CommittedRegressionCorpus_PassesTheFullPipeline()
    {
        var corpus = FuzzRegressions.Load();
        TestContext.Out.WriteLine($"{corpus.Count} regression entr{(corpus.Count == 1 ? "y" : "ies")}");
        foreach (var r in corpus)
        {
            Assert.That(FuzzOracle.Check(r.MinimizedText, FuzzOracle.ReferencePipeline), Is.Null, $"minimized text of {r.Case.Id} ({r.Op})");
            Assert.That(FuzzOracle.Check(r.Case.Text(), FuzzOracle.ReferencePipeline), Is.Null, $"regenerated case {r.Case.Id} ({r.Op})");
        }
    }

    [Test]
    public void SmokeRun_EveryTier_FindsNothing()
    {
        foreach (var (tier, cases) in new[] { (FuzzTier.Tiny, 3000L), (FuzzTier.Small, 300L), (FuzzTier.Large, 4L) })
        {
            var s = FuzzRunner.Run(new FuzzOptions { MasterSeed = 20261103, Tier = tier, Cases = cases, Reference = FuzzOracle.ReferencePipeline });
            TestContext.Out.WriteLine(s.ToString());
            Assert.That(s.Failures, Is.Empty, string.Join("\n", s.Failures.Select(FuzzRegressions.Format)));
            Assert.That(s.Cases, Is.EqualTo(cases));
        }
    }
}
