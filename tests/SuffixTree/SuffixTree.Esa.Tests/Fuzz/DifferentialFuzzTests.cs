using System.Globalization;

namespace SuffixTree.Esa.Tests.Fuzz;

/// <summary>
/// Acceptance-size differential fuzz of work package 1a.7: <c>SaisReference</c> + <c>LcpReference</c> against the naive SA/LCP (n &lt;= 2000) and
/// against the independent verifier + Kasai (n &gt; 2000). Category <c>Fuzz</c>, <c>[Explicit]</c>. Run all three tiers with
/// <c>dotnet test tests/SuffixTree/SuffixTree.Esa.Tests/SuffixTree.Esa.Tests.csproj -c Release --filter Category=Fuzz</c>.
/// Environment parameters: <c>ESA_FUZZ_SEED</c> (master seed, default 1), <c>ESA_FUZZ_SCALE</c> (fraction of the acceptance case counts, default 1),
/// <c>ESA_FUZZ_BUDGET_S</c> (time budget in seconds per tier, default none; a budget-limited run is reported as such and is not the acceptance run).
/// A failure prints <c>family:seed:n:sigma:param</c>, the failed operation and the minimized text, writes the regression lines to
/// <c>fuzz-failures.txt</c> in the test work directory, and the lines are then appended to <c>Fuzz/regressions.txt</c> once the bug is fixed.
/// </summary>
[TestFixture]
[Category("Fuzz")]
[Explicit("Long differential fuzz (about 4 minutes at the acceptance sizes); run with --filter Category=Fuzz (see PROGRESS-S1.md, Log 1a.7)")]
public class DifferentialFuzzTests
{
    [Test]
    public void Tiny_1e6_Cases_UpTo64() => RunTier(FuzzTier.Tiny);

    [Test]
    public void Small_1e5_Cases_UpTo2000() => RunTier(FuzzTier.Small);

    [Test]
    public void Large_1e3_Cases_UpTo1e5() => RunTier(FuzzTier.Large);

    private static void RunTier(FuzzTier tier)
    {
        ulong seed = ulong.Parse(Environment.GetEnvironmentVariable("ESA_FUZZ_SEED") ?? "1", CultureInfo.InvariantCulture);
        double scale = double.Parse(Environment.GetEnvironmentVariable("ESA_FUZZ_SCALE") ?? "1", CultureInfo.InvariantCulture);
        string? budgetText = Environment.GetEnvironmentVariable("ESA_FUZZ_BUDGET_S");
        TimeSpan? budget = budgetText is null ? null : TimeSpan.FromSeconds(double.Parse(budgetText, CultureInfo.InvariantCulture));
        long cases = Math.Max(1, (long)(tier.Cases * scale));

        var summary = FuzzRunner.Run(new FuzzOptions
        {
            MasterSeed = seed,
            Tier = tier,
            Cases = cases,
            TimeBudget = budget,
            Log = TestContext.Progress.WriteLine,
        });
        TestContext.Out.WriteLine(summary.ToString());
        TestContext.Out.WriteLine(summary.Cases == tier.Cases && !summary.BudgetExhausted
            ? $"acceptance size reached: {summary.Cases} cases"
            : $"NOT the acceptance run: {summary.Cases} of {tier.Cases} cases (scale {scale}, budget {budgetText ?? "none"})");

        if (summary.Failures.Count > 0)
        {
            string lines = string.Join("\n", summary.Failures.Select(FuzzRegressions.Format));
            File.WriteAllText(Path.Combine(TestContext.CurrentContext.WorkDirectory, "fuzz-failures.txt"), lines + "\n");
            Assert.Fail($"{summary.Failures.Count} mismatch(es); regression lines:\n{lines}\nfirst: {summary.Failures[0].Op}: {summary.Failures[0].Detail}");
        }
    }
}
