using System.Security.Cryptography;
using SuffixTree.EmProto.Tests.Ph5;
using SuffixTree.EmProto.Tests.WholePipeline;
using SuffixTree.EmProto.Verify;
using SuffixTree.Esa.Tests.TestSupport;
using static SuffixTree.EmProto.Tests.E8.VerifierHarness;

namespace SuffixTree.EmProto.Tests.E8;

/// <summary>
/// The verifier of the E8 package (ADR-004 section 5, row E8): clean outputs pass, and every kind of injected fault is detected by the check that owns it (swapped suffix-array entries, a changed LCP value, a corrupt child entry,
/// a truncated file, an off-by-one in n, ...); every single-byte flip of the outputs of small texts is detected; bounded memory options never turn into a pass.
/// </summary>
[TestFixture]
public class VerifierTests : TempDirFixture
{
    private static readonly string[] Small = ["a", "ab", "banana", "mississippi", "abracadabra", "aaaaaaaaab", "abababababab", "acgtacgtnnnnacgtacgt"];

    private string Reference(string text, string name = "ref") => FromReference(Run(name), text);

    // ---- clean outputs pass

    [Test]
    public void CleanPipelineOutputs_PassEveryCheck()
    {
        int v = 0;
        foreach (string s in Small)
        {
            string prefix = FromPipeline(Run($"p{v}"), s, v);
            VerifyResult r = Check(prefix, 1 + (v % 3), 3 + (v % 5), input: v % 2 == 0);
            Assert.That(r.AllPassed, Is.True, $"'{s}': {Describe(r)}");
            string[] expected = OutputVerifier.Names.Where(n => n != "input" || v % 2 == 0).ToArray();
            Assert.That(r.Checks.Select(c => c.Name), Is.EqualTo(expected), $"'{s}'");
            v++;
        }
    }

    [Test]
    public void CleanReferenceOutputs_AllWidthsOfTheLcp_Pass()
    {
        foreach ((string label, string text, string width, bool overflow) in new[]
        {
            ("u8 plain", StringFamilies.RandomDna(5000, 3), "1", false),
            ("u8 + escapes", ChildSamples.U8WithEscapes(), "1", true),
            ("u16 + overflow", new string('a', 70_000), "2", true),
            ("u32", new string('a', 200_000), "4", false),
        })
        {
            string prefix = Reference(text, label.Replace(' ', '_'));
            VerifyResult r = OutputVerifier.Run(Options(prefix, 3, 4096));
            Assert.That(r.AllPassed, Is.True, $"{label}: {Describe(r)}");
            Assert.That(r.Facts["lcpwidth"], Is.EqualTo(width), label);
            Assert.That(long.Parse(r.Facts["overflow"]) > 0, Is.EqualTo(overflow), label);
        }
    }

    [Test]
    public void Exhaustive_SmallStrings_AllPass()
    {
        int count = 0;
        foreach ((int sigma, int max) in new[] { (2, 10), (3, 5) })
        {
            string d = Run("x" + sigma);
            foreach (string s in ExhaustiveStrings.All(sigma, max))
            {
                if (s.Length == 0) continue;
                string prefix = FromReference(d, s);
                VerifyResult r = OutputVerifier.Run(Options(prefix, 1 + (count % 2), 2 + (count % 7), hash: false));
                if (!r.AllPassed) Assert.Fail($"'{s}': {Describe(r)}");
                count++;
            }
        }
        Assert.That(count, Is.EqualTo(2046 + 363));
    }

    [Test]
    public void GenomeLikeText_WithNRunsAndSatellites_Passes_WithTheDefaultCap_AndWithoutCap()
    {
        string text = System.Text.Encoding.ASCII.GetString(global::SuffixTree.EmProto.Gen.Families.ToAscii(global::SuffixTree.EmProto.Gen.Families.GenomeLike(30_000, 20240002)));
        string prefix = Reference(text, "gl");
        foreach (long cap in new long[] { 0, 64, 65536 })
        {
            VerifyResult r = OutputVerifier.Run(Options(prefix, 3, 997, cap));
            Assert.That(r.AllPassed, Is.True, $"cap {cap}: {Describe(r)}");
        }
    }

    [Test]
    public void ResultsAreTheSame_ForAnyThreadCount_AndBlockSize()
    {
        string prefix = Reference(SmallWithEscapes(), "tb");
        VerifyResult first = OutputVerifier.Run(Options(prefix, 1, 1 << 20, 16));
        foreach ((int t, int b) in new[] { (1, 1), (2, 7), (4, 64), (3, 1000), (4, 4097) })
        {
            VerifyResult r = OutputVerifier.Run(Options(prefix, t, b, 16));
            Assert.That(r.AllPassed, Is.True, $"{t} threads, block {b}: {Describe(r)}");
            Assert.That(r.Facts["long_pairs"], Is.EqualTo(first.Facts["long_pairs"]), $"{t} threads, block {b}");
            Assert.That(r.Sha256, Is.EqualTo(first.Sha256));
        }
    }

    // ---- injected faults

    [Test]
    public void SwappedSuffixArrayEntries_AreDetected_ByTheOrderCheck()
    {
        string clean = Reference("mississippi");
        foreach ((int i, int j) in new[] { (1, 2), (3, 4), (2, 9), (1, 11), (5, 6), (0, 1) })
        {
            string prefix = Clone(clean, Run($"s{i}_{j}"));
            uint a = ReadUInt(prefix + ".sa", i), b = ReadUInt(prefix + ".sa", j);
            WriteUInt(prefix + ".sa", i, b);
            WriteUInt(prefix + ".sa", j, a);
            VerifyResult r = Check(prefix);
            Assert.That(r.AllPassed, Is.False, $"swap {i},{j}");
            Assert.That(r.Failed("sa_order") || r.Failed("sa_permutation"), Is.True, $"swap {i},{j}: {Describe(r)}");
            Assert.That(r.Passed("sa_permutation") || i == 0 || j == 0, Is.True, $"a swap keeps the permutation: {Describe(r)}");
        }
    }

    [Test]
    public void SwapInALargerText_IsDetected_AtAnyDistance()
    {
        string clean = Reference(StringFamilies.RandomDna(3000, 9), "big");
        var rng = new SplitMix64(5);
        for (int k = 0; k < 25; k++)
        {
            int i = 1 + rng.NextInt(3000), j = 1 + rng.NextInt(3000);
            if (i == j) continue;
            string prefix = Clone(clean, Run($"w{k}"));
            uint a = ReadUInt(prefix + ".sa", i), b = ReadUInt(prefix + ".sa", j);
            WriteUInt(prefix + ".sa", i, b);
            WriteUInt(prefix + ".sa", j, a);
            VerifyResult r = Check(prefix, 3, 100);
            Assert.That(r.Failed("sa_order"), Is.True, $"swap {i},{j}: {Describe(r)}");
        }
    }

    [Test]
    public void DuplicatedOrOutOfRangeSuffixArrayEntry_IsDetected_ByThePermutationCheck()
    {
        string clean = Reference("banana");
        string d1 = Clone(clean, Run("dup"));
        WriteUInt(d1 + ".sa", 3, ReadUInt(d1 + ".sa", 4));
        Assert.That(Check(d1).Failed("sa_permutation"), Is.True);
        string d2 = Clone(clean, Run("range"));
        WriteUInt(d2 + ".sa", 3, 7);
        Assert.That(Check(d2).Failed("sa_permutation"), Is.True);
        string d3 = Clone(clean, Run("huge"));
        WriteUInt(d3 + ".sa", 3, uint.MaxValue);
        VerifyResult r3 = Check(d3);
        Assert.That(r3.Failed("sa_permutation"), Is.True, Describe(r3));
    }

    [Test]
    public void SaEntryZero_MustBeN()
    {
        string clean = Reference("banana");
        string d = Clone(clean, Run("sa0"));
        uint first = ReadUInt(d + ".sa", 0), other = ReadUInt(d + ".sa", 1);
        WriteUInt(d + ".sa", 0, other);
        WriteUInt(d + ".sa", 1, first);
        VerifyResult r = Check(d);
        Assert.That(r.AllPassed, Is.False);
        Assert.That(r.Failed("sa_permutation") || r.Failed("sa_order"), Is.True, Describe(r));
    }

    [Test]
    public void OneChangedLcpValue_IsDetected_ByTheLcpCheck_AndByNothingElseInTheSa()
    {
        string clean = Reference("mississippi");
        for (int rank = 0; rank <= 11; rank++)
        {
            foreach (int delta in new[] { 1, -1, 5 })
            {
                string prefix = Clone(clean, Run($"l{rank}_{delta}"));
                byte old = File.ReadAllBytes(prefix + ".lcp")[rank];
                int changed = old + delta;
                if (changed < 0 || changed > 250) continue;
                FlipByte(prefix + ".lcp", rank, (byte)(old ^ changed));
                VerifyResult r = Check(prefix);
                Assert.That(r.Failed("lcp_values"), Is.True, $"rank {rank} {old} -> {changed}: {Describe(r)}");
                Assert.That(r.Passed("sa_order") && r.Passed("sa_permutation"), Is.True, $"the SA itself is intact: {Describe(r)}");
            }
        }
    }

    [Test]
    public void ChangedOverflowEntries_AreDetected_WithAndWithoutTheCap()
    {
        string clean = Reference(SmallWithEscapes(), "ovf");
        long pairs = new FileInfo(clean + ".ovf").Length / 8;
        Assert.That(pairs, Is.GreaterThan(10));
        foreach (long cap in new long[] { 0, 16, 100 })
        {
            foreach (int pair in new[] { 0, 1, (int)pairs / 2, (int)pairs - 1 })
            {
                foreach (int field in new[] { 0, 1 })
                {
                    foreach (int delta in new[] { 1, -1 })
                    {
                        string prefix = Clone(clean, Run($"o{cap}_{pair}_{field}_{delta}"));
                        uint v = ReadUInt(prefix + ".ovf", (2L * pair) + field);
                        WriteUInt(prefix + ".ovf", (2L * pair) + field, (uint)(v + delta));
                        VerifyResult r = OutputVerifier.Run(Options(prefix, 2, 50, cap));
                        Assert.That(r.Failed("lcp_values"), Is.True, $"cap {cap}, pair {pair} field {field} {delta:+0;-0}: {Describe(r)}");
                    }
                }
            }
        }
    }

    [Test]
    public void CorruptChildEntries_AreDetected_ByTheChildCheck_AndOnlyByIt()
    {
        string clean = Reference("mississippi");
        long cells = new FileInfo(clean + ".cld").Length / 4;
        for (long k = 0; k < cells; k++)
        {
            foreach (uint value in new uint[] { ReadUInt(clean + ".cld", k) + 1, 0, 1, uint.MaxValue, 12 })
            {
                if (value == ReadUInt(clean + ".cld", k)) continue;
                string prefix = Clone(clean, Run("c"));
                WriteUInt(prefix + ".cld", k, value);
                VerifyResult r = Check(prefix);
                Assert.That(r.Failed("child"), Is.True, $"cell {k} = {value}: {Describe(r)}");
                Assert.That(r.Passed("sa_order") && r.Passed("lcp_values") && r.Passed("stats"), Is.True, Describe(r));
            }
        }
    }

    [Test]
    public void CorruptChildEntry_InALargeTable_IsDetected()
    {
        string clean = Reference(StringFamilies.RandomDna(4000, 13), "cb");
        var rng = new SplitMix64(7);
        for (int t = 0; t < 30; t++)
        {
            string prefix = Clone(clean, Run("cb" + t));
            long k = rng.NextInt(4001);
            WriteUInt(prefix + ".cld", k, ReadUInt(prefix + ".cld", k) ^ (1u << rng.NextInt(12)));
            Assert.That(Check(prefix, 2, 500).Failed("child"), Is.True, $"cell {k}");
        }
    }

    [Test]
    public void TruncatedOrMissingFile_IsDetected_WithoutACrash()
    {
        string clean = Reference(SmallWithEscapes(), "tr");
        foreach (string e in Extensions)
        {
            long length = new FileInfo(clean + "." + e).Length;
            foreach (string how in new[] { "minus1", "zero", "delete", "plus1" })
            {
                string prefix = Clone(clean, Run($"t_{e}_{how}"));
                string f = prefix + "." + e;
                if (how == "minus1") Truncate(f, length - 1);
                else if (how == "zero") Truncate(f, 0);
                else if (how == "delete") File.Delete(f);
                else File.AppendAllText(f, e is "pk.meta" or "stats" ? "x" : "\0");
                VerifyResult r = OutputVerifier.Run(Options(prefix, 2, 50));
                Assert.That(r.AllPassed, Is.False, $"{e} {how}");
                Assert.That(r.Checks.Any(c => c.Status == CheckStatus.Fail), Is.True, $"{e} {how}: {Describe(r)}");
                Assert.That(r.Failed("files") || r.Failed("meta") || r.Failed("stats"), Is.True, $"{e} {how}: {Describe(r)}");
            }
        }
    }

    [Test]
    public void OffByOneInN_IsDetected_InTheMetadata_TheStatistics_AndTheFileLengths()
    {
        string clean = Reference("mississippi");
        foreach (int delta in new[] { 1, -1 })
        {
            string m = Clone(clean, Run($"mn{delta}"));
            SetKey(m + ".pk.meta", "n", (11 + delta).ToString());
            Assert.That(Check(m).AllPassed, Is.False, $"meta n {delta:+0;-0}");
            string s = Clone(clean, Run($"sn{delta}"));
            SetKey(s + ".stats", "n", (11 + delta).ToString());
            Assert.That(Check(s).AllPassed, Is.False, $"stats n {delta:+0;-0}");
            string both = Clone(clean, Run($"bn{delta}"));
            SetKey(both + ".pk.meta", "n", (11 + delta).ToString());
            SetKey(both + ".stats", "n", (11 + delta).ToString());
            VerifyResult r = Check(both);
            Assert.That(r.Failed("files"), Is.True, $"meta and stats agree on a wrong n, the files do not: {Describe(r)}");
        }
    }

    [Test]
    public void FlippedTextSymbol_AndNonZeroPadding_AreDetected_ByTheTextCheck()
    {
        string clean = Reference("mississippi");
        string d = Clone(clean, Run("tx"));
        FlipByte(d + ".pk", 1, 0x20);
        VerifyResult r = Check(d);
        Assert.That(r.Failed("text"), Is.True, Describe(r));
        string p = Clone(clean, Run("pad"));
        long length = new FileInfo(p + ".pk").Length;
        FlipByte(p + ".pk", length - 3, 0x01);
        Assert.That(Check(p).Failed("text"), Is.True, "padding");
    }

    [Test]
    public void TamperedStatistics_AreDetected_KeyByKey()
    {
        string clean = Reference(SmallWithEscapes(), "st");
        foreach (string key in new[] { "nodecount", "internal", "lrsdepth", "lrsrank", "n", "sigma", "lcpwidth", "overflow", "maxlcp" })
        {
            foreach (int delta in new[] { 1, -1 })
            {
                string prefix = Clone(clean, Run($"k{key}{delta}"));
                SetKey(prefix + ".stats", key, (long.Parse(Key(prefix + ".stats", key)) + delta).ToString());
                VerifyResult r = OutputVerifier.Run(Options(prefix, 2, 100));
                Assert.That(r.AllPassed, Is.False, $"{key} {delta:+0;-0}");
                Assert.That(r.Failed("stats") || r.Failed("files") || r.Failed("meta") || r.Failed("lcp_width"), Is.True, $"{key} {delta:+0;-0}: {Describe(r)}");
            }
        }
    }

    [Test]
    public void ALcpWidthThatBreaksTheRule_IsDetected_EvenWhenEveryValueIsRight()
    {
        string text = StringFamilies.RandomDna(2000, 4);
        string clean = Reference(text, "wr");
        string d = Clone(clean, Run("w16"));
        byte[] u8 = File.ReadAllBytes(d + ".lcp");
        var u16 = new byte[u8.Length * 2];
        for (int i = 0; i < u8.Length; i++) u16[2 * i] = u8[i];
        File.WriteAllBytes(d + ".lcp", u16);
        SetKey(d + ".stats", "lcpwidth", "2");
        VerifyResult r = Check(d, 2, 100);
        Assert.That(r.Passed("lcp_values"), Is.True, Describe(r));
        Assert.That(r.Failed("lcp_width"), Is.True, Describe(r));
    }

    [Test]
    public void ChangedInput_IsDetected_ByTheInputCheck()
    {
        string clean = Reference("mississippi");
        string d = Clone(clean, Run("in1"));
        byte[] input = File.ReadAllBytes(d.Replace("out", "in") + ".txt");
        _ = input;
        string path = Path.Combine(Path.GetDirectoryName(d)!, "in.txt");
        byte[] bytes = File.ReadAllBytes(path);
        bytes[4] = (byte)'z';
        File.WriteAllBytes(path, bytes);
        Assert.That(Check(d).Failed("input"), Is.True);
        string e = Clone(clean, Run("in2"));
        File.AppendAllText(Path.Combine(Path.GetDirectoryName(e)!, "in.txt"), "a");
        Assert.That(Check(e).Failed("input"), Is.True, "a longer input");
        string f = Clone(clean, Run("in3"));
        Truncate(Path.Combine(Path.GetDirectoryName(f)!, "in.txt"), 10);
        Assert.That(Check(f).Failed("input"), Is.True, "a shorter input");
        string g = Clone(clean, Run("in4"));
        File.Delete(Path.Combine(Path.GetDirectoryName(g)!, "in.txt"));
        Assert.That(Check(g).Failed("input"), Is.True, "a missing input");
    }

    [Test]
    public void EverySingleByteFlip_OfTheOutputs_IsDetected()
    {
        int flips = 0;
        foreach (string text in new[] { "mississippi", "banana", "aaaaaaaa", "acgtnacgtacgtnnt", "a" })
        {
            string clean = Reference(text, "f" + text.Length + text[0]);
            string[] metaLines = File.ReadAllLines(clean + ".pk.meta");
            long nameStart = File.ReadAllBytes(clean + ".pk.meta").Length - (metaLines[^1].Length + 1); // the `name=` line is not checked
            foreach (string e in Extensions)
            {
                long length = new FileInfo(clean + "." + e).Length;
                for (long off = 0; off < length; off++)
                {
                    if (e == "pk.meta" && off >= nameStart) continue;
                    foreach (byte mask in new byte[] { 0x01, 0x80 })
                    {
                        string prefix = Clone(clean, Run("b"));
                        FlipByte(prefix + "." + e, off, mask);
                        VerifyResult r = OutputVerifier.Run(Options(prefix, 1, 3, hash: false));
                        Assert.That(r.AllPassed, Is.False, $"'{text}' {e} byte {off} ^ {mask:X2}");
                        flips++;
                        Directory.Delete(Path.GetDirectoryName(prefix)!, true);
                    }
                }
            }
        }
        Assert.That(flips, Is.GreaterThan(1500));
    }

    [Test]
    public void FlipsInALargerText_WithOverflow_AreDetected_ForTheOverflowFile_AndSampledElsewhere()
    {
        string clean = Reference(SmallWithEscapes(), "fl");
        var rng = new SplitMix64(11);
        int flips = 0;
        foreach (string e in new[] { "ovf", "lcp", "sa", "cld", "pk" })
        {
            long length = new FileInfo(clean + "." + e).Length;
            long stride = e == "ovf" ? 1 : Math.Max(1, length / 40);
            for (long off = rng.NextInt(3); off < length; off += stride)
            {
                string prefix = Clone(clean, Run("lf"));
                FlipByte(prefix + "." + e, off, (byte)(1 << rng.NextInt(8)));
                VerifyResult r = OutputVerifier.Run(Options(prefix, 2, 500, 32, hash: false));
                Assert.That(r.AllPassed, Is.False, $"{e} byte {off}");
                flips++;
                Directory.Delete(Path.GetDirectoryName(prefix)!, true);
            }
        }
        Assert.That(flips, Is.GreaterThan(300));
    }

    // ---- the cap, the stack bound, the optional checks, the hashes

    [Test]
    public void LongPairs_AreCounted_AndCheckedAtTheirEndsAndByProbes()
    {
        string prefix = Reference(SmallWithEscapes(), "lp");
        VerifyResult capped = OutputVerifier.Run(Options(prefix, 2, 100, 16));
        VerifyResult exact = OutputVerifier.Run(Options(prefix, 2, 100, 0));
        Assert.That(capped.AllPassed && exact.AllPassed, Is.True, Describe(capped) + " | " + Describe(exact));
        Assert.That(long.Parse(capped.Facts["long_pairs"]), Is.GreaterThan(100));
        Assert.That(long.Parse(exact.Facts["long_pairs"]), Is.EqualTo(0));
        // a text with no pair above the cap
        Assert.That(long.Parse(OutputVerifier.Run(Options(prefix, 2, 100, 100_000)).Facts["long_pairs"]), Is.EqualTo(0));
    }

    [Test]
    public void LongPairProbes_CatchACorruptionInTheMatchingPart_OfAnSaEntry()
    {
        // two copies of a 300-symbol block: the adjacent suffixes of the copies share 300 symbols; the SA entry of one copy is replaced by a suffix that agrees on the first 20 symbols only
        string block = StringFamilies.RandomDna(300, 31);
        string text = StringFamilies.RandomDna(1000, 32) + block + StringFamilies.RandomDna(50, 33) + block;
        string clean = Reference(text, "pr");
        int caught = 0, tried = 0;
        var rng = new SplitMix64(3);
        for (int t = 0; t < 40; t++)
        {
            string prefix = Clone(clean, Run("pr" + t));
            int i = 1 + rng.NextInt(text.Length), j = 1 + rng.NextInt(text.Length);
            if (i == j) continue;
            uint a = ReadUInt(prefix + ".sa", i), b = ReadUInt(prefix + ".sa", j);
            WriteUInt(prefix + ".sa", i, b);
            WriteUInt(prefix + ".sa", j, a);
            tried++;
            if (!OutputVerifier.Run(Options(prefix, 2, 64, 16)).AllPassed) caught++;
        }
        Assert.That(caught, Is.EqualTo(tried));
    }

    [Test]
    public void AStackBoundBelowTheLcp_SkipsTheScan_NeverPasses()
    {
        string prefix = Reference(new string('a', 2000), "sk");
        VerifyResult r = OutputVerifier.Run(Options(prefix, 2, 100, stackCap: 100));
        Assert.That(r.Find("child")!.Status, Is.EqualTo(CheckStatus.Skipped), Describe(r));
        Assert.That(r.AllPassed, Is.False);
        Assert.That(r.Passed("sa_order") && r.Passed("lcp_values"), Is.True, Describe(r));
        Assert.That(OutputVerifier.Run(Options(prefix, 2, 100, stackCap: 2001)).AllPassed, Is.True);
        Assert.That(OutputVerifier.Run(Options(prefix, 2, 100, stackCap: 2000)).Find("child")!.Status, Is.EqualTo(CheckStatus.Skipped), "the stack of a^n holds n + 1 entries");
    }

    [Test]
    public void OptionalChecks_AreOmittedWhenSwitchedOff()
    {
        string prefix = Reference("mississippi");
        VerifyResult noChild = OutputVerifier.Run(Options(prefix, 1, 4, child: false));
        Assert.That(noChild.Find("child"), Is.Null);
        Assert.That(noChild.Find("stats"), Is.Null);
        Assert.That(noChild.AllPassed, Is.True, Describe(noChild));
        VerifyResult noHash = OutputVerifier.Run(Options(prefix, 1, 4, hash: false));
        Assert.That(noHash.Find("sha256"), Is.Null);
        Assert.That(noHash.Sha256, Is.Empty);
        VerifyResult noInput = OutputVerifier.Run(Options(prefix, 1, 4, input: false));
        Assert.That(noInput.Find("input"), Is.Null);
        Assert.That(noInput.AllPassed, Is.True, Describe(noInput));
    }

    [Test]
    public void Sha256_IsTheHashOfEveryFile_AndChangesWithAFlip()
    {
        string prefix = Reference(SmallWithEscapes(), "sh");
        VerifyResult r = Check(prefix, 2, 500, 64);
        Assert.That(r.Sha256.Keys.OrderBy(k => k), Is.EqualTo(Extensions.OrderBy(k => k)));
        foreach (string e in Extensions) Assert.That(r.Sha256[e], Is.EqualTo(Convert.ToHexStringLower(SHA256.HashData(File.ReadAllBytes(prefix + "." + e)))), e);
        string d = Clone(prefix, Run("shf"));
        FlipByte(d + ".sa", 40, 0x04);
        VerifyResult r2 = Check(d, 2, 500, 64);
        Assert.That(r2.Sha256["sa"], Is.Not.EqualTo(r.Sha256["sa"]));
        Assert.That(r2.Sha256["cld"], Is.EqualTo(r.Sha256["cld"]));
    }

    [Test]
    public void ABlockLargerThanTheText_AndOneRankPerBlock_GiveTheSameVerdict()
    {
        string clean = Reference("abracadabra");
        foreach (int block in new[] { 1, 2, 12, 13, 1 << 20 })
            Assert.That(Check(clean, 3, block).AllPassed, Is.True, $"block {block}");
    }

    [Test]
    public void TheBitmapFile_IsRemoved_AndNothingElseIsLeftBehind()
    {
        string d = Run("clean");
        string prefix = FromReference(d, "mississippi");
        string[] before = Directory.GetFiles(d).Select(Path.GetFileName).OrderBy(x => x).ToArray()!;
        Assert.That(Check(prefix).AllPassed, Is.True);
        Assert.That(Directory.GetFiles(d).Select(Path.GetFileName).OrderBy(x => x), Is.EqualTo(before));
        string failing = Clone(prefix, Run("clean2"));
        WriteUInt(failing + ".sa", 3, 5);
        Check(failing);
        Assert.That(Directory.GetFiles(Path.GetDirectoryName(failing)!).Length, Is.EqualTo(before.Length));
    }

    // ---- the command

    [Test]
    public void Command_ExitCodes_Rows_AndSums()
    {
        string clean = Reference(SmallWithEscapes(), "cmd");
        string rows = Path.Combine(Dir, "rows.tsv");
        string sums = Path.Combine(Dir, "sums.txt");
        int ok = Cli.Run(["e8-verify", "--prefix", clean, "--input", Path.Combine(Path.GetDirectoryName(clean)!, "in.txt"), "--threads", "2", "--lce-cap", "100", "--rows", rows, "--sums", sums, "--case", "t"]);
        Assert.That(ok, Is.EqualTo(0));
        string text = File.ReadAllText(rows);
        foreach (string name in OutputVerifier.Names) Assert.That(text, Does.Contain($"\tverify_{name}\t"), name);
        Assert.That(text, Does.Contain("\tverify_all\tpass\t"));
        Assert.That(text, Does.Contain("\tsha256_sa\t"));
        Assert.That(File.ReadAllLines(sums).Select(l => l.Split('=', 2)[0]), Is.EqualTo(Extensions));
        string bad = Clone(clean, Run("cmdbad"));
        FlipByte(bad + ".cld", 8, 0x01);
        Assert.That(Cli.Run(["e8-verify", "--prefix", bad, "--rows", Path.Combine(Dir, "rows2.tsv")]), Is.EqualTo(1));
        Assert.That(File.ReadAllText(Path.Combine(Dir, "rows2.tsv")), Does.Contain("\tverify_all\tFAIL\t").And.Contain("\tverify_child\tFAIL\t"));
        Assert.That(Cli.Run(["e8-verify"]), Is.EqualTo(2));
        Assert.That(Cli.Run(["e8-verify", "--prefix", Path.Combine(Dir, "nothing")]), Is.EqualTo(1), "missing files are a failed verification, not a crash");
    }
}
