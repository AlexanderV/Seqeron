using System.Text;
using SuffixTree.EmProto.Tests.WholePipeline;
using SuffixTree.EmProto.Verify;
using SuffixTree.Esa;
using SuffixTree.Esa.Tests.TestSupport;
using Prod = SuffixTree.Esa.External;
using ProtoPh5 = SuffixTree.EmProto.Ph5;
using ProtoPipe = SuffixTree.EmProto.Pipeline;

namespace SuffixTree.EmProto.Tests.External;

/// <summary>Helpers of the differential tests of work package 1b.7c: the prototype <c>e6</c> makes the inputs (packed text, suffix array) and the reference outputs, the production stage runs on the prototype's own files.</summary>
internal static class LcpChildDiff
{
    internal static ProtoPipe.EndToEndResult RunProto(DiffDir dir, string inputPath, string name, int threads = 4, ProtoPipe.EndToEndTuning? tuning = null, bool exact = false, int period = 1024)
        => ProtoPipe.EndToEnd.Run(new ProtoPipe.EndToEndOptions
        {
            InputPath = inputPath,
            OutPrefix = dir.File(name),
            ScratchDir = dir.File(name + "-scratch"),
            Threads = threads,
            SyncOutputs = false,
            Tuning = tuning ?? new ProtoPipe.EndToEndTuning(),
            Period = period,
            ExactPeriod = exact,
        });

    /// <summary>The production stage on the packed text and the suffix array of a prototype run; the outputs go to <paramref name="outName"/> next to them.</summary>
    internal static (Prod.LcpChildResult Result, Prod.LcpChildOutputs Files) RunProd(DiffDir dir, ProtoPipe.EndToEndResult proto, string outName, Prod.LcpChildOptions options, Action<Prod.MappedTextReader>? configure = null)
    {
        var files = new Prod.LcpChildOutputs(dir.File(outName + ".lcp"), dir.File(outName + ".ovf"), dir.File(outName + ".cld"), dir.File(outName + ".stats"));
        using Prod.MappedTextReader text = Prod.MappedTextReader.Open(proto.Files.Packed, proto.Prepare.N, proto.Prepare.Bits);
        configure?.Invoke(text);
        Prod.LcpChildResult result = Prod.LcpChildBuilder.Build(proto.Files.Sa, text, files, options with { Sigma = proto.Prepare.Sigma });
        return (result, files);
    }

    internal static void AssertSameFiles(ProtoPipe.EndToEndResult proto, Prod.LcpChildOutputs prod, string label)
    {
        Assert.That(DiffInputs.SameBytes(proto.Files.Lcp, prod.Lcp), Is.True, $"{label}: .lcp differs from the prototype's");
        Assert.That(DiffInputs.SameBytes(proto.Files.Overflow, prod.Overflow), Is.True, $"{label}: .ovf differs from the prototype's");
        Assert.That(DiffInputs.SameBytes(proto.Files.Child, prod.Child), Is.True, $"{label}: .cld differs from the prototype's");
        Assert.That(File.ReadAllText(prod.Stats!), Is.EqualTo(File.ReadAllText(proto.Files.Stats)), $"{label}: .stats differs from the prototype's");
    }

    internal static void AssertSameCounts(ProtoPipe.EndToEndResult proto, Prod.LcpChildResult prod, string label)
    {
        Assert.That(prod.Lcp.Kind, Is.EqualTo(proto.Lcp.Kind), label);
        Assert.That((prod.Lcp.Escapes8, prod.Lcp.Escapes16, prod.Lcp.OverflowCount, prod.Lcp.MaxPlcp), Is.EqualTo((proto.Lcp.Escapes8, proto.Lcp.Escapes16, proto.Lcp.OverflowCount, proto.Lcp.MaxPlcp)), $"{label}: counts of the LCP build");
        Assert.That((prod.Scan.NodeCount, prod.Scan.InternalCount, prod.Scan.LrsDepth, prod.Scan.LrsRank, prod.Scan.MaxLcp), Is.EqualTo((proto.Scan.NodeCount, proto.Scan.InternalCount, proto.Scan.LrsDepth, proto.Scan.LrsRank, proto.Scan.MaxLcp)), $"{label}: statistics of the scan");
        Assert.That(prod.StatsText(proto.Prepare.Sigma), Is.EqualTo(ProtoPipe.EndToEnd.StatsText(proto)), $"{label}: statistics text");
    }

    internal static string TextFile(DiffDir dir, string name, string text)
    {
        string path = dir.File(name + ".txt");
        File.WriteAllBytes(path, Encoding.Latin1.GetBytes(text));
        return path;
    }

    internal static string DataFile(string name)
    {
        string path = Path.Combine(DiffInputs.DataDir, name + ".txt");
        if (!File.Exists(path)) Assert.Ignore($"{path} is not present (generate it with analysis/bench/gen/gen_inputs.py)");
        return path;
    }
}

/// <summary>
/// G2 of work package 1b.7c: the production stage (PH4 base mode fused with PH5, layout A, codec Packed) on the suffix array file of the prototype <c>e6</c> writes the <c>.lcp .ovf .cld .stats</c> files byte for byte as <c>e6</c> did, with the
/// same counts, the same counted I/O, and the verifier of E8 (<c>e8-verify</c>) passes all eleven checks on them. The prototype code is unchanged.
/// </summary>
[TestFixture]
public class LcpChildDifferentialTests
{
    private static Prod.LcpChildOptions Default => new() { ScratchDir = string.Empty };

    private static void Differential(DiffDir dir, string input, string label, Prod.LcpChildOptions? options = null, int threads = 4)
    {
        ProtoPipe.EndToEndResult proto = LcpChildDiff.RunProto(dir, input, "proto", threads);
        (Prod.LcpChildResult result, Prod.LcpChildOutputs files) = LcpChildDiff.RunProd(dir, proto, "prod", (options ?? Default) with { ScratchDir = dir.File("prod-scratch"), Threads = threads });
        LcpChildDiff.AssertSameFiles(proto, files, label);
        LcpChildDiff.AssertSameCounts(proto, result, label);
        Assert.That(Directory.Exists(dir.File("prod-scratch")), Is.False, $"{label}: scratch");
    }

    [Test]
    public void TheBenchmarkFiles_AreByteIdentical_ToTheE6Outputs()
    {
        foreach (string name in new[] { "prot_1m", "dna_10m", "rep_10m" })
        {
            string input = LcpChildDiff.DataFile(name);
            using var dir = new DiffDir();
            Differential(dir, input, name);
        }
    }

    [Test]
    public void TheFamilies_AtTenToTheSixth_AreByteIdentical_ToTheE6Outputs()
    {
        foreach ((string name, string text) in DiffInputs.Families(1_000_000))
        {
            using var dir = new DiffDir();
            string input = LcpChildDiff.TextFile(dir, "in", text);
            Differential(dir, input, name);
        }
    }

    [Test]
    public void RandomDnaAndRepeatRich_AtTenToTheSixth_WithTheScaledGeometry_AreByteIdentical()
    {
        foreach ((string name, string text) in new[] { ("dna", StringFamilies.RandomDna(1_000_000, 1)), ("repeat-rich", LcpRepeatRich(1_000_000)), ("n-runs", NRuns(1_000_000)) })
        {
            using var dir = new DiffDir();
            string input = LcpChildDiff.TextFile(dir, "in", text);
            Differential(dir, input, name, new Prod.LcpChildOptions { ScratchDir = string.Empty, ScaleN = 3_000_000_000 });
        }
    }

    private static string LcpRepeatRich(int n)
    {
        char[] t = StringFamilies.RandomDna(n, 19).ToCharArray();
        var rng = new SplitMix64(1019);
        for (int copy = 0; copy < 6; copy++)
        {
            int len = n / 12;
            int from = rng.NextInt(n - len), to = rng.NextInt(n - len);
            Array.Copy(t, from, t, to, len);
            for (int s = 0; s < len / 50; s++) t[to + rng.NextInt(len)] = "ACGT"[rng.NextInt(4)];
        }
        return new string(t);
    }

    private static string NRuns(int n)
    {
        char[] t = StringFamilies.RandomDna(n, 20).ToCharArray();
        var rng = new SplitMix64(1020);
        for (int run = 0; run < 3; run++) Array.Fill(t, 'N', rng.NextInt(n - (n / 15)), n / 15);
        return new string(t);
    }

    [Test]
    public void SmallTexts_WithRotatingTinyParameters_AreByteIdentical_ToTheE6Outputs()
    {
        var rng = new SplitMix64(1_818);
        string[] letters = ["a", "ab", "acgt", "abcdefgh"];
        for (int v = 0; v < 120; v++)
        {
            int n = 1 + rng.NextInt(v % 4 == 0 ? 600 : 90);
            string alphabet = letters[rng.NextInt(letters.Length)];
            var sb = new StringBuilder();
            for (int i = 0; i < n; i++) sb.Append(alphabet[rng.NextInt(alphabet.Length)]);
            string text = v % 5 == 0 ? new string('a', n) : v % 5 == 1 ? StringFamilies.Fibonacci(n) : sb.ToString();
            using var dir = new DiffDir();
            string input = LcpChildDiff.TextFile(dir, "in", text);
            ProtoPipe.EndToEndResult proto = LcpChildDiff.RunProto(dir, input, "proto", 1 + (v % 3), PipelineChecks.Tiny(v), exact: true, period: PipelineChecks.Period(v));
            var tuning = new Prod.LcpChildTuning
            {
                LcpChunk = 3 + ((v * 11) % 48),
                LcpWorkers = 1 + (v % 3),
                LcpReadBytes = 16 + (8 * (v % 8)),
                LcpTupleBytes = 48 + (4 * (v % 5)),
                WindowBytes = 64 + (32 * (v % 6)),
                BlockEntries = 2 + (v % 5),
                ChildWindowBytes = 16 << (v % 3),
                ChildSlots = 1 + (v % 3),
                SpillFlushBytes = 64,
            };
            (Prod.LcpChildResult result, Prod.LcpChildOutputs files) = LcpChildDiff.RunProd(dir, proto, "prod", new Prod.LcpChildOptions { ScratchDir = dir.File("prod-scratch"), Threads = 1 + (v % 4), Tuning = tuning, ArenaGuards = ArenaOptions.Guarded(0xCD) });
            LcpChildDiff.AssertSameFiles(proto, files, $"v={v} n={n}");
            LcpChildDiff.AssertSameCounts(proto, result, $"v={v} n={n}");
        }
    }

    [Test]
    public void TheCountedIo_IsTheSame_AsThePrototypesInThePhase()
    {
        foreach ((string name, string text) in new[] { ("dna", StringFamilies.RandomDna(400_000, 3)), ("a^n", new string('a', 300_000)) })
        {
            using var dir = new DiffDir();
            string input = LcpChildDiff.TextFile(dir, "in", text);
            ProtoPipe.EndToEndResult proto = LcpChildDiff.RunProto(dir, input, "proto");
            ProtoPipe.PhaseRecord phase = proto.Phase(ProtoPipe.EndToEnd.PhaseLcpChild);
            var counters = new Prod.IoCounters();
            (Prod.LcpChildResult result, _) = LcpChildDiff.RunProd(dir, proto, "prod", new Prod.LcpChildOptions { ScratchDir = dir.File("prod-scratch"), Counters = counters });
            Assert.That(result.Io.ReadBytes, Is.EqualTo(phase.ExplicitReadBytes), $"{name}: explicit reads");
            Assert.That(result.Io.WriteBytes, Is.EqualTo(phase.ExplicitWriteBytes), $"{name}: explicit writes");
            Assert.That(result.Io.MappedWriteBytes - (4 * result.Scan.ChildWindows.Writes), Is.EqualTo(phase.MappedWriteBytes), $"{name}: the mapped writes of the LCP files (the CHILD windows are counted in addition)");
            Assert.That(result.Scan.Spill, Is.EqualTo(ToProd(proto.Scan.Spill)), $"{name}: the spill, field by field");
            Assert.That(result.Scan.ChildWindows, Is.EqualTo(ToProd(proto.Scan.ChildWindows)), $"{name}: the CHILD windows");
            Assert.That(result.Scan.StackPeak, Is.EqualTo(proto.Scan.StackPeak), $"{name}: the stack peak");
        }
    }

    private static Prod.SpillStats ToProd(ProtoPh5.SpillStats s) => new(s.BlocksWritten, s.BlocksRead, s.BytesWritten, s.BytesRead, s.PeakFileBytes, s.Flushes, s.CleanEvictions);

    private static Prod.WindowStats ToProd(ProtoPh5.WindowStats s) => new(s.Opened, s.Reopened, s.Flushes, s.FlushedBytes, s.Writes, s.Evictions, s.DistinctWindows);

    [Test]
    public void ThePrototypeAndTheProduction_SpillTheSameBytes_OnATheN()
    {
        using var dir = new DiffDir();
        string input = LcpChildDiff.TextFile(dir, "in", new string('a', 1_000_000));
        ProtoPipe.EndToEndResult proto = LcpChildDiff.RunProto(dir, input, "proto");
        (Prod.LcpChildResult result, _) = LcpChildDiff.RunProd(dir, proto, "prod", new Prod.LcpChildOptions { ScratchDir = dir.File("prod-scratch") });
        Assert.That(result.Scan.Spill.BytesWritten, Is.EqualTo(proto.Scan.Spill.BytesWritten));
        Assert.That(result.Scan.Spill.BytesWritten, Is.GreaterThan(15_000_000), "16 bytes per entry beyond two blocks of 256 K entries");
        Assert.That(result.Scan.Spill.BytesRead, Is.EqualTo(proto.Scan.Spill.BytesRead));
    }

    [Test]
    public void TheVerifierOfE8_PassesAllElevenChecks_OnTheProductionOutputs()
    {
        foreach ((string name, string text) in new[] { ("dna", StringFamilies.RandomDna(300_000, 3)), ("rep", LcpRepeatRich(300_000)), ("fib", StringFamilies.Fibonacci(200_000)), ("a^n", new string('a', 50_000)) })
        {
            using var dir = new DiffDir();
            string input = LcpChildDiff.TextFile(dir, "in", text);
            ProtoPipe.EndToEndResult proto = LcpChildDiff.RunProto(dir, input, "proto");
            // the production outputs under a prefix of their own, with the prototype's text, metadata and suffix array (hard links: the verifier reads files by prefix)
            string prefix = dir.File("check");
            File.Copy(proto.Files.Packed, prefix + ".pk");
            File.Copy(proto.Files.PackedMeta, prefix + ".pk.meta");
            File.Copy(proto.Files.Sa, prefix + ".sa");
            var files = new Prod.LcpChildOutputs(prefix + ".lcp", prefix + ".ovf", prefix + ".cld", prefix + ".stats");
            using (Prod.MappedTextReader reader = Prod.MappedTextReader.Open(prefix + ".pk", proto.Prepare.N, proto.Prepare.Bits))
                Prod.LcpChildBuilder.Build(prefix + ".sa", reader, files, new Prod.LcpChildOptions { ScratchDir = dir.File("prod-scratch"), Sigma = proto.Prepare.Sigma });
            VerifyResult v = OutputVerifier.Run(new VerifyOptions { Prefix = prefix, InputPath = input, ScratchDir = dir.File("verify-scratch"), Threads = 2 });
            Assert.That(v.Checks, Has.Count.EqualTo(11), name);
            Assert.That(v.AllPassed, Is.True, $"{name}: {string.Join("; ", v.Checks.Where(c => c.Status != CheckStatus.Pass).Select(c => $"{c.Name} {c.Status} {c.Detail}"))}");
        }
    }

    [Test]
    public void TheScalarPath_GivesTheSameFiles_AsThePrototypesVectorPath()
    {
        using var dir = new DiffDir();
        string input = LcpChildDiff.TextFile(dir, "in", LcpRepeatRich(500_000));
        ProtoPipe.EndToEndResult proto = LcpChildDiff.RunProto(dir, input, "proto");
        (Prod.LcpChildResult result, Prod.LcpChildOutputs files) = LcpChildDiff.RunProd(dir, proto, "prod", new Prod.LcpChildOptions { ScratchDir = dir.File("prod-scratch") }, text => text.UseVectors = false);
        LcpChildDiff.AssertSameFiles(proto, files, "scalar");
        Assert.That(result.Lcp.CompareSteps, Is.GreaterThan(proto.Lcp.CompareSteps), "the scalar path takes more (smaller) steps for the same values");
    }
}
