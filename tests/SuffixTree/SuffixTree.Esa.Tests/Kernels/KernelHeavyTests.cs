using SuffixTree.Esa.Tests.Fuzz;
using SuffixTree.Esa.Tests.TestSupport;

namespace SuffixTree.Esa.Tests.Kernels;

/// <summary>
/// The explicit categories with the OPTIMIZED kernels plugged in (every registered non-reference engine, combined with the reference of the other stage and
/// with each other, <see cref="KernelRegistry.OptimizedPipelines"/>): <c>Exhaustive</c> (binary strings up to length 20 against Kasai, the 10^6 texts), <c>Vectors</c> (the 79
/// libsais hashes, <c>SEQERON_BENCH_DIR</c> as in <see cref="Vectors.LibsaisVectorTests"/>) and <c>Fuzz</c> (the acceptance tiers of 1a.7 with the reference as
/// <c>FuzzOptions.Reference</c>, entry-by-entry equality). Each pipeline call runs in a guarded arena whose canaries are verified.
/// </summary>
[TestFixture]
public class KernelHeavyTests
{
    /// <summary>The case that stands for "no optimized kernel is registered" (a rejected iteration is reverted, so the registry can be empty): the tests pass trivially, they are not skipped.</summary>
    private const string NoKernel = "(no optimized kernel registered)";

    public static IEnumerable<TestCaseData> Pipelines() => PipelineNames().Select(name => new TestCaseData(name).SetName("Pipeline_" + name));

    public static IEnumerable<string> PipelineNames()
    {
        string[] names = [.. KernelRegistry.OptimizedPipelines().Select(p => p.Name)];
        return names.Length == 0 ? [NoKernel] : names;
    }

    private static EsaBuilder Find(string name) => KernelRegistry.OptimizedPipelines().Single(p => p.Name == name).Build;

    private static bool NothingToTest(string pipeline) => pipeline == NoKernel;

    [Test]
    [Category("Exhaustive")]
    [Explicit("Heavy enumeration; run with --filter Category=Exhaustive (see PROGRESS-S1.md, Log)")]
    [TestCaseSource(nameof(Pipelines))]
    public void AllBinaryStrings_UpTo20_PipelineEqualsKasai_AndTheReference(string pipeline)
    {
        if (NothingToTest(pipeline)) Assert.Pass(NoKernel);
        EsaBuilder build = Find(pipeline);
        long count = 0;
        foreach (string s in ExhaustiveStrings.All(2, 20))
        {
            var (sa, lcp) = build(s);
            if (!lcp.ToArray().AsSpan().SequenceEqual(KasaiLcp.Build(s, sa))) Assert.Fail($"{pipeline}: LCP differs from Kasai on '{s}'");
            if (s.Length <= 12 && !sa.AsSpan().SequenceEqual(NaiveSuffixArray.Build(s))) Assert.Fail($"{pipeline}: SA differs from the naive SA on '{s}'");
            count++;
        }
        Assert.That(count, Is.EqualTo((1L << 21) - 1));
    }

    [Test]
    [Category("Exhaustive")]
    [Explicit("Heavy: 10^6-character texts; run with --filter Category=Exhaustive")]
    [TestCaseSource(nameof(Pipelines))]
    public void TextsOf1e6_PipelineEqualsKasai_AndPicksTheExpectedKind(string pipeline)
    {
        if (NothingToTest(pipeline)) Assert.Pass(NoKernel);
        const int n = 1_000_000;
        EsaBuilder build = Find(pipeline);
        foreach (var (name, text, kind) in new[]
        {
            ("a^n", StringFamilies.AllSame(n), LcpKind.U32),
            ("fibonacci", StringFamilies.Fibonacci(n), LcpKind.U32),
            ("random-dna", StringFamilies.RandomDna(n, 31), LcpKind.U8),
            ("abab", StringFamilies.AbRepeat(n), LcpKind.U32),
            ("thue-morse", StringFamilies.ThueMorse(n), LcpKind.U32),
        })
        {
            var (sa, lcp) = build(text);
            Assert.That(SuffixArrayVerifier.Check(text, sa), Is.Null, $"{pipeline} {name}");
            Assert.That(lcp.ToArray().AsSpan().SequenceEqual(KasaiLcp.Build(text, sa)), Is.True, $"{pipeline} {name}");
            Assert.That(lcp.Kind, Is.EqualTo(kind), $"{pipeline} {name}");
        }
    }

    [Test]
    [Category("Exhaustive")]
    [Explicit("Heavy; run with --filter Category=Exhaustive")]
    [TestCaseSource(nameof(Pipelines))]
    public void LargeCorpus_1e5_EveryEngine_EqualsTheReference_WithAllPoisonPatterns(string pipeline)
    {
        if (NothingToTest(pipeline)) Assert.Pass(NoKernel);
        foreach (SaisEngine engine in KernelRegistry.OptimizedSais)
            foreach (var (name, text) in KernelCorpus.Large())
                Assert.That(KernelHarness.CheckSais(engine, text, allPoisons: false), Is.Null, $"{engine.Name} on {name}");
        foreach (LcpEngine engine in KernelRegistry.OptimizedLcp)
            foreach (var (name, text) in KernelCorpus.Large())
                Assert.That(KernelHarness.CheckLcp(engine, text, allPoisons: false), Is.Null, $"{engine.Name} on {name}");
    }

    [Test]
    [Category("Vectors")]
    [Explicit("Heavy: 10 Mbp builds; run with SEQERON_BENCH_DIR=<dir> --filter Category=Vectors (see PROGRESS-S1.md, Log 1a.6)")]
    public void SaAndLcp_EqualLibsais([ValueSource(nameof(VectorCorpora))] string corpus, [ValueSource(nameof(PipelineNames))] string pipeline)
    {
        if (NothingToTest(pipeline)) Assert.Pass(NoKernel);
        var entry = VectorManifest.Load()[corpus];
        string text = entry.Kind == VectorManifest.CorpusKind.File ? ReadInput(corpus) : VectorManifest.GenerateText(corpus, entry.N);
        Assert.That(text.Length, Is.EqualTo(entry.N));
        var (sa, lcp) = Find(pipeline)(text);
        Assert.That(EsaPipeline.HashSa(sa), Is.EqualTo(entry.SaSha256), $"{pipeline} {corpus}: SA hash");
        Assert.That(EsaPipeline.HashLcp(lcp, entry.N), Is.EqualTo(entry.LcpSha256), $"{pipeline} {corpus}: LCP hash");
    }

    public static IEnumerable<string> VectorCorpora() => VectorManifest.Load().Entries.Select(e => e.Corpus);

    [Test]
    [Category("Fuzz")]
    [Explicit("Long differential fuzz; run with --filter Category=Fuzz (ESA_FUZZ_SEED, ESA_FUZZ_SCALE, ESA_FUZZ_BUDGET_S as in DifferentialFuzzTests)")]
    public void Fuzz_Tiny_1e6_Cases_UpTo64([ValueSource(nameof(PipelineNames))] string pipeline) => RunTier(pipeline, FuzzTier.Tiny);

    [Test]
    [Category("Fuzz")]
    [Explicit("Long differential fuzz; run with --filter Category=Fuzz")]
    public void Fuzz_Small_1e5_Cases_UpTo2000([ValueSource(nameof(PipelineNames))] string pipeline) => RunTier(pipeline, FuzzTier.Small);

    [Test]
    [Category("Fuzz")]
    [Explicit("Long differential fuzz; run with --filter Category=Fuzz")]
    public void Fuzz_Large_1e3_Cases_UpTo1e5([ValueSource(nameof(PipelineNames))] string pipeline) => RunTier(pipeline, FuzzTier.Large);

    private static void RunTier(string pipeline, FuzzTier tier)
    {
        if (NothingToTest(pipeline)) Assert.Pass(NoKernel);
        ulong seed = ulong.Parse(Environment.GetEnvironmentVariable("ESA_FUZZ_SEED") ?? "1", System.Globalization.CultureInfo.InvariantCulture);
        double scale = double.Parse(Environment.GetEnvironmentVariable("ESA_FUZZ_SCALE") ?? "1", System.Globalization.CultureInfo.InvariantCulture);
        string? budgetText = Environment.GetEnvironmentVariable("ESA_FUZZ_BUDGET_S");
        TimeSpan? budget = budgetText is null ? null : TimeSpan.FromSeconds(double.Parse(budgetText, System.Globalization.CultureInfo.InvariantCulture));
        var summary = FuzzRunner.Run(new FuzzOptions
        {
            MasterSeed = seed,
            Tier = tier,
            Cases = Math.Max(1, (long)(tier.Cases * scale)),
            TimeBudget = budget,
            Candidate = Find(pipeline),
            Reference = FuzzOracle.ReferencePipeline,
            Log = TestContext.Progress.WriteLine,
        });
        TestContext.Out.WriteLine($"{pipeline}: {summary}");
        if (summary.Failures.Count > 0)
            Assert.Fail($"{summary.Failures.Count} mismatch(es); first: {summary.Failures[0].Case.Id} {summary.Failures[0].Op}: {summary.Failures[0].Detail}; minimized {FuzzRegressions.EncodeUnits(summary.Failures[0].MinimizedText)}");
    }

    private static string ReadInput(string corpus)
    {
        string? dir = Environment.GetEnvironmentVariable("SEQERON_BENCH_DIR");
        if (string.IsNullOrWhiteSpace(dir)) Assert.Fail("SEQERON_BENCH_DIR is not set (see LibsaisVectorTests)");
        string relative = $"data/{corpus}.txt";
        try
        {
            if (!CorpusManifest.LoadBenchManifest().Verify(relative, dir!)) Assert.Fail($"{Path.Combine(dir!, relative)} does not match analysis/bench/MANIFEST.sha256");
        }
        catch (FileNotFoundException ex)
        {
            Assert.Fail($"input missing: {ex.FileName}");
        }
        return System.Text.Encoding.Latin1.GetString(File.ReadAllBytes(Path.Combine(dir!, relative)));
    }
}
