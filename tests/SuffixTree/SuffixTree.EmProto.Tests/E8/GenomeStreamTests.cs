using System.Security.Cryptography;
using SuffixTree.EmProto.Experiments;
using SuffixTree.EmProto.Gen;
using SuffixTree.EmProto.Harness;
using SuffixTree.EmProto.Tests.WholePipeline;

namespace SuffixTree.EmProto.Tests.E8;

/// <summary>
/// The streamed genome-like generator of the E8 package (ADR-004 section 5, row E8): the same text as <c>Families.GenomeLike</c> (the family of E3 to E7) as ASCII bytes, for any buffer size, with a checksum, in bounded memory.
/// </summary>
[TestFixture]
public class GenomeStreamTests : TempDirFixture
{
    private static byte[] Stream(long n, ulong seed, int buffer)
    {
        var ms = new MemoryStream();
        GenomeStream.Generate(n, seed, block => ms.Write(block), buffer);
        return ms.ToArray();
    }

    [TestCase(1L, 20240002UL)]
    [TestCase(2L, 20240002UL)]
    [TestCase(7L, 1UL)]
    [TestCase(199L, 20240002UL)]
    [TestCase(1000L, 20240002UL)]
    [TestCase(65536L, 1UL)]
    [TestCase(100_000L, 20240002UL)]
    [TestCase(1_234_567L, 42UL)]
    [TestCase(3_000_000L, 20240002UL)]
    public void Equals_TheInRamFamily_ByteForByte(long n, ulong seed)
    {
        byte[] expected = Families.ToAscii(Families.GenomeLike(n, seed));
        Assert.That(Stream(n, seed, 1 << 20), Is.EqualTo(expected));
    }

    [Test]
    public void EveryLengthUpTo300_EqualsTheInRamFamily()
    {
        for (long n = 1; n <= 300; n++)
            Assert.That(Stream(n, 20240002, 64), Is.EqualTo(Families.ToAscii(Families.GenomeLike(n, 20240002))), $"n = {n}");
    }

    [Test]
    public void IdenticalForAnyBufferSize()
    {
        byte[] reference = Stream(150_000, 20240002, 1 << 20);
        foreach (int buffer in new[] { 1, 2, 7, 64, 4096, 4097, 100_000, 150_000, 1 << 24 })
            Assert.That(Stream(150_000, 20240002, buffer), Is.EqualTo(reference), $"buffer {buffer}");
    }

    [Test]
    public void Blocks_NeverExceedTheBuffer_AndSumToN()
    {
        long total = 0;
        int largest = 0;
        GenomeStream.Generate(500_000, 20240002, block => { total += block.Length; largest = Math.Max(largest, block.Length); }, 4096);
        Assert.That(total, Is.EqualTo(500_000));
        Assert.That(largest, Is.LessThanOrEqualTo(4096));
    }

    [Test]
    public void Result_HasTheSha256_TheCounts_AndTheLength()
    {
        byte[] bytes = Stream(250_000, 20240002, 1 << 20);
        GenomeStreamResult r = GenomeStream.Generate(250_000, 20240002, _ => { }, 1 << 16);
        Assert.That(r.N, Is.EqualTo(250_000));
        Assert.That(r.Seed, Is.EqualTo(20240002UL));
        Assert.That(r.Sha256, Is.EqualTo(Convert.ToHexStringLower(SHA256.HashData(bytes))));
        long[] counts = GenomeStream.Letters.Select(c => (long)bytes.Count(b => b == c)).ToArray();
        Assert.That(r.Counts, Is.EqualTo(counts));
        Assert.That(counts.Sum(), Is.EqualTo(250_000));
    }

    [Test]
    public void StatisticalCharacter_NRunsAndRepeatsAreThere()
    {
        // the family: about 0.85 % N in six runs, 40 % interspersed repeats with 10-15 % divergence, 4 % satellites; the counts of A, C, G, T of the unique background are equal
        GenomeStreamResult r = GenomeStream.Generate(4_000_000, 20240002, _ => { }, 1 << 20);
        double nShare = r.Counts[3] / 4e6;
        Assert.That(nShare, Is.InRange(0.0075, 0.0095));
        long acgt = r.Counts[0] + r.Counts[1] + r.Counts[2] + r.Counts[4];
        foreach (int c in new[] { 0, 1, 2, 4 }) Assert.That(r.Counts[c] / (double)acgt, Is.InRange(0.2, 0.3), "base " + GenomeStream.Letters[c]);
    }

    [Test]
    public void BoundedMemory_TheAllocationsDoNotGrowWithN()
    {
        long before = GC.GetAllocatedBytesForCurrentThread();
        GenomeStream.Generate(40_000_000, 20240002, _ => { }, 1 << 20);
        long allocated = GC.GetAllocatedBytesForCurrentThread() - before;
        Assert.That(allocated, Is.LessThan(24L << 20), "allocations for 40 M symbols");
    }

    [Test]
    public void SumsText_HasEveryField()
    {
        GenomeStreamResult r = GenomeStream.Generate(10_000, 7, _ => { }, 1 << 10);
        Dictionary<string, string> kv = GenomeStream.SumsText(r).Split('\n', StringSplitOptions.RemoveEmptyEntries).Select(l => l.Split('=', 2)).ToDictionary(p => p[0], p => p[1]);
        Assert.That(kv["n"], Is.EqualTo("10000"));
        Assert.That(kv["seed"], Is.EqualTo("7"));
        Assert.That(kv["sha256"], Is.EqualTo(r.Sha256));
        Assert.That(kv["alphabet"], Is.EqualTo("ACGNT"));
        Assert.That(kv["counts"], Is.EqualTo(string.Join(',', r.Counts)));
        Assert.That(kv.ContainsKey("format"), Is.True);
    }

    [Test]
    public void ZeroOrNegative_N_IsRefused()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => GenomeStream.Generate(0, 1, _ => { }));
        Assert.Throws<ArgumentOutOfRangeException>(() => GenomeStream.Generate(-5, 1, _ => { }));
        Assert.Throws<ArgumentOutOfRangeException>(() => GenomeStream.Generate(10, 1, _ => { }, 0));
    }

    [Test]
    public void Command_WritesTheFile_AndTheChecksumFile()
    {
        string file = Path.Combine(Run("gen"), "g.txt");
        string sums = file + ".sums";
        int rc = Cli.Run(["e8-gen", "--n", "123456", "--out", file, "--sums", sums, "--rows", Path.Combine(Dir, "rows.tsv")]);
        Assert.That(rc, Is.EqualTo(0));
        byte[] bytes = File.ReadAllBytes(file);
        Assert.That(bytes, Is.EqualTo(Families.ToAscii(Families.GenomeLike(123456, 20240002))));
        Dictionary<string, string> kv = File.ReadAllLines(sums).Select(l => l.Split('=', 2)).ToDictionary(p => p[0], p => p[1]);
        Assert.That(kv["sha256"], Is.EqualTo(Convert.ToHexStringLower(SHA256.HashData(bytes))));
        Assert.That(kv["n"], Is.EqualTo("123456"));
        Assert.That(File.ReadAllText(Path.Combine(Dir, "rows.tsv")), Does.Contain("\tE8\t").And.Contain("gen_sha256"));
    }

    [Test]
    public void Command_AcceptsDecimalSuffixes_AndASeed_AndRefusesGarbage()
    {
        string a = Path.Combine(Run("suf"), "a.txt");
        Assert.That(Cli.Run(["e8-gen", "--n", "2k", "--out", a, "--seed", "5"]), Is.EqualTo(0));
        Assert.That(File.ReadAllBytes(a), Is.EqualTo(Families.ToAscii(Families.GenomeLike(2000, 5))));
        Assert.That(Cli.Run(["e8-gen", "--n", "abc", "--out", a]), Is.EqualTo(2));
        Assert.That(Cli.Run(["e8-gen", "--out", a]), Is.EqualTo(2));
        Assert.That(E8Parsing.Count("1g"), Is.EqualTo(1_000_000_000L));
        Assert.That(E8Parsing.Count("2.2G"), Is.EqualTo(2_200_000_000L));
        Assert.That(E8Parsing.Count("3000m"), Is.EqualTo(3_000_000_000L));
        Assert.That(E8Parsing.Count("123"), Is.EqualTo(123L));
    }

    [Test]
    public void Command_Refuses_ToOverwriteWithoutForce_AndRemovesAPartialFileOnFailure()
    {
        string a = Path.Combine(Run("ow"), "a.txt");
        File.WriteAllText(a, "x");
        Assert.That(Cli.Run(["e8-gen", "--n", "100", "--out", a]), Is.EqualTo(2));
        Assert.That(File.ReadAllText(a), Is.EqualTo("x"));
        Assert.That(Cli.Run(["e8-gen", "--n", "100", "--out", a, "--force"]), Is.EqualTo(0));
        Assert.That(new FileInfo(a).Length, Is.EqualTo(100));
    }

    [Test]
    public void EnvCommand_ReportsTheArchitecture_TheCores_AndTheVectorPath()
    {
        string rows = Path.Combine(Dir, "env.tsv");
        Assert.That(Cli.Run(["e8-env", "--rows", rows]), Is.EqualTo(0));
        string[][] lines = File.ReadAllLines(rows).Select(l => l.Split('\t')).ToArray();
        string Value(string metric) => lines.First(l => l[3] == metric)[4];
        Assert.That(int.Parse(Value("env_cores")), Is.EqualTo(Environment.ProcessorCount));
        Assert.That(Value("env_process_arch"), Is.EqualTo(System.Runtime.InteropServices.RuntimeInformation.ProcessArchitecture.ToString()));
        Assert.That(Value("env_page_size"), Is.EqualTo(Environment.SystemPageSize.ToString()));
        Assert.That(Value("env_little_endian"), Is.EqualTo("yes"));
        Assert.That(Value("env_lce_path"), Is.EqualTo(System.Runtime.Intrinsics.X86.Avx2.IsSupported ? "avx2 256-bit vectors" : "scalar 64-bit words"));
        Assert.That(long.Parse(Value("env_gc_available_memory")), Is.GreaterThan(0));
        foreach (string m in new[] { "env_framework", "env_os", "env_os_arch", "env_runtime_id", "env_gc_hard_limit", "env_avx2", "env_advsimd" }) Assert.That(Value(m), Is.Not.Empty, m);
    }
}
