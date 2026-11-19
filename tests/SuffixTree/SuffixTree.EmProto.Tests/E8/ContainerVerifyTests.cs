using System.Buffers.Binary;
using SuffixTree.Esa.External;
using SuffixTree.Esa.Tests.TestSupport;
using SuffixTree.EmProto.Tests.WholePipeline;
using SuffixTree.EmProto.Verify;
using static SuffixTree.EmProto.Tests.E8.VerifierHarness;

namespace SuffixTree.EmProto.Tests.E8;

/// <summary>
/// The container adapter of <c>e8-verify</c> (work package 1b.7f, <c>--container</c>): the eleven checks run on the sections of a SESA container read in place from its mapping, with no extraction to files. A container written by
/// the end-to-end builder passes every check; each injected fault in a section is found by the check that owns it; the result equals that of the same sections extracted to files.
/// </summary>
[TestFixture]
public class ContainerVerifyTests : TempDirFixture
{
    private static readonly string[] Small = ["a", "ab", "banana", "mississippi", "abracadabra", "aaaaaaaaab", "abababababab", "acgtacgtnnnnacgtacgt", "GATTACAGATTACAGGTAGATTACA"];

    /// <summary>Builds the container of the text with the end-to-end builder; the input text is <c>in.txt</c> next to it.</summary>
    private string Build(string text, string name, bool psi = true)
    {
        string dir = Run(name);
        string input = Path.Combine(dir, "in.txt");
        File.WriteAllBytes(input, PipelineChecks.Latin1(text));
        string target = Path.Combine(dir, "out.sesa");
        SesaBuilder.Build(input, target, new SesaBuildOptions { Threads = 2, Fsync = false, MeasureRss = false, IncludePsi = psi });
        return target;
    }

    private static VerifyOptions Options(string container, bool input = true)
        => new()
        {
            Prefix = container,
            InputPath = input ? Path.Combine(Path.GetDirectoryName(container)!, "in.txt") : null,
            Threads = 2,
            BlockRanks = 5,
            ScanBlock = 7,
            IoBufferBytes = 13,
            ScratchDir = Path.GetDirectoryName(container),
        };

    private static VerifyResult Check(string container, bool input = true)
    {
        using var source = new ContainerVerifySource(container);
        return OutputVerifier.Run(Options(container, input), source);
    }

    [Test]
    public void ABuiltContainer_PassesAllElevenChecks()
    {
        int v = 0;
        foreach (string s in Small)
        {
            string container = Build(s, $"c{v}", psi: v % 2 == 0);
            VerifyResult r = Check(container, input: v % 3 != 0);
            Assert.That(r.AllPassed, Is.True, $"'{s}': {Describe(r)}");
            Assert.That(r.Checks.Count, Is.EqualTo(v % 3 != 0 ? 11 : 10), $"'{s}'");
            Assert.That(r.Checks.Select(c => c.Name), Is.EqualTo(OutputVerifier.Names.Where(n => n != "input" || v % 3 != 0)));
            v++;
        }
    }

    [Test]
    public void AContainerWithAnOverflowList_PassesAllElevenChecks_AndTheFactsSayWhat()
    {
        string container = Build(SmallWithEscapes(), "esc");
        VerifyResult r = Check(container);
        Assert.That(r.AllPassed, Is.True, Describe(r));
        Assert.That(long.Parse(r.Facts["overflow"]), Is.GreaterThan(0));
        Assert.That(r.Facts["lcpwidth"], Is.EqualTo("1"));
        Assert.That(r.Facts["n"], Is.EqualTo(SmallWithEscapes().Length.ToString()));
    }

    [Test]
    public void TheResult_EqualsThatOfTheSameSectionsExtractedToFiles()
    {
        string container = Build(SmallWithEscapes(), "x");
        VerifyResult inPlace = Check(container);
        // the same seven outputs as files, written from the sections: pk, sa, lcp, ovf, cld verbatim, pk.meta and stats as the adapter derives them
        string dir = Run("extracted");
        string prefix = Path.Combine(dir, "out");
        using (var source = new ContainerVerifySource(container))
        {
            foreach (string e in new[] { "pk", "sa", "lcp", "ovf", "cld" })
            {
                using var block = source.Open(e);
                byte[] bytes;
                unsafe
                {
                    bytes = new ReadOnlySpan<byte>(block.Pointer, (int)block.Bytes).ToArray();
                }
                File.WriteAllBytes($"{prefix}.{e}", bytes);
            }
            File.WriteAllText(prefix + ".pk.meta", source.ReadText("pk.meta"));
            File.WriteAllText(prefix + ".stats", source.ReadText("stats"));
        }
        File.Copy(Path.Combine(Path.GetDirectoryName(container)!, "in.txt"), Path.Combine(dir, "in.txt"));
        VerifyResult files = CheckFiles(prefix);
        Assert.That(files.AllPassed, Is.True, Describe(files));
        Assert.That(inPlace.Checks.Select(c => (c.Name, c.Status)), Is.EqualTo(files.Checks.Select(c => (c.Name, c.Status))));
        Assert.That(inPlace.Facts, Is.EqualTo(files.Facts));
        Assert.That(inPlace.Sha256, Is.EqualTo(files.Sha256), "the SHA-256 of every section equals that of the extracted file");
    }

    [Test]
    public void NothingIsExtracted_TheDirectoryHoldsTheContainerAndTheInputOnly()
    {
        string container = Build("mississippi", "d");
        string dir = Path.GetDirectoryName(container)!;
        string[] before = Directory.GetFileSystemEntries(dir).Order().ToArray();
        Assert.That(Check(container).AllPassed, Is.True);
        Assert.That(Directory.GetFileSystemEntries(dir).Order(), Is.EqualTo(before), "the permutation bitmap is a temporary file that is removed; no section is written out");
    }

    // ---- injected faults, in place

    private static (byte[] Image, SesaOracleLike O) Image(string container) => (File.ReadAllBytes(container), SesaOracleLike.Of(File.ReadAllBytes(container)));

    /// <summary>The section table of a container image, decoded for the tests of this class only.</summary>
    private sealed record SesaOracleLike(Dictionary<int, (long Offset, long Length, int Index)> Sections)
    {
        internal static SesaOracleLike Of(byte[] image)
        {
            int count = BinaryPrimitives.ReadInt32LittleEndian(image.AsSpan(48));
            var map = new Dictionary<int, (long, long, int)>();
            for (int i = 0; i < count; i++)
            {
                ReadOnlySpan<byte> e = image.AsSpan(128 + (64 * i), 64);
                map[BinaryPrimitives.ReadInt32LittleEndian(e)] = (BinaryPrimitives.ReadInt64LittleEndian(e[8..]), BinaryPrimitives.ReadInt64LittleEndian(e[16..]), i);
            }
            return new SesaOracleLike(map);
        }
    }

    /// <summary>The header CRC-32C again after an edit of the header or the table (the adapter opens the container at its header level).</summary>
    private static void ResealHeader(byte[] image)
    {
        int count = BinaryPrimitives.ReadInt32LittleEndian(image.AsSpan(48));
        uint crc = Crc32C.Append(Crc32C.Append(0, image.AsSpan(0, 124)), image.AsSpan(128, 64 * count));
        BinaryPrimitives.WriteUInt32LittleEndian(image.AsSpan(124), crc);
    }

    private VerifyResult Damaged(string text, string name, Action<byte[], SesaOracleLike> damage, bool input = true)
    {
        string container = Build(text, name);
        (byte[] image, SesaOracleLike o) = Image(container);
        damage(image, o);
        File.WriteAllBytes(container, image);
        return Check(container, input);
    }

    [Test]
    public void ASwappedPairOfSuffixArrayEntries_IsFoundByTheOrderCheck()
    {
        VerifyResult r = Damaged("GATTACAGATTACAGGTAGATTACA", "swap", (img, o) =>
        {
            long at = o.Sections[SesaFormat.KindSa].Offset + (4 * 5);
            Span<byte> a = img.AsSpan((int)at, 4), b = img.AsSpan((int)at + 4, 4);
            byte[] t = a.ToArray();
            b.CopyTo(a);
            t.CopyTo(b);
        });
        Assert.That(r.Failed("sa_order"), Is.True, Describe(r));
        Assert.That(r.Passed("sa_permutation"), Is.True, Describe(r));
    }

    [Test]
    public void ADuplicatedSuffixArrayEntry_IsFoundByThePermutationCheck()
    {
        VerifyResult r = Damaged("GATTACAGATTACAGGTAGATTACA", "dup", (img, o) =>
        {
            long at = o.Sections[SesaFormat.KindSa].Offset;
            img.AsSpan((int)at + (4 * 7), 4).CopyTo(img.AsSpan((int)at + (4 * 6), 4));
        });
        Assert.That(r.Failed("sa_permutation"), Is.True, Describe(r));
    }

    [Test]
    public void AWrongLcpValue_IsFoundByTheLcpCheck()
    {
        VerifyResult r = Damaged("GATTACAGATTACAGGTAGATTACA", "lcp", (img, o) => img[o.Sections[SesaFormat.KindLcp].Offset + 9] ^= 0x01);
        Assert.That(r.Failed("lcp_values"), Is.True, Describe(r));
    }

    [Test]
    public void AWrongChildCell_IsFoundByTheChildCheck()
    {
        VerifyResult r = Damaged("GATTACAGATTACAGGTAGATTACA", "child", (img, o) =>
        {
            long cells = o.Sections[SesaFormat.KindChild].Length / 4;
            long at = o.Sections[SesaFormat.KindChild].Offset;
            for (long k = 1; k < cells; k++)
            {
                int p = (int)(at + (4 * k));
                if (BinaryPrimitives.ReadUInt32LittleEndian(img.AsSpan(p)) != 0)
                {
                    BinaryPrimitives.WriteUInt32LittleEndian(img.AsSpan(p), BinaryPrimitives.ReadUInt32LittleEndian(img.AsSpan(p)) ^ 1);
                    return;
                }
            }
        });
        Assert.That(r.Failed("child"), Is.True, Describe(r));
    }

    [Test]
    public void AChangedTextSymbol_IsFoundByTheTextCheck_AndTheInputCheck()
    {
        VerifyResult r = Damaged("GATTACAGATTACAGGTAGATTACA", "text", (img, o) => img[o.Sections[SesaFormat.KindText].Offset] ^= 0x80);
        Assert.That(r.Failed("text"), Is.True, Describe(r));
        Assert.That(r.Failed("input"), Is.True, Describe(r));
    }

    [Test]
    public void AWrongNodeCountInTheHeader_IsFoundByTheStatsCheck()
    {
        VerifyResult r = Damaged("GATTACAGATTACAGGTAGATTACA", "stats", (img, _) =>
        {
            BinaryPrimitives.WriteInt64LittleEndian(img.AsSpan(24), BinaryPrimitives.ReadInt64LittleEndian(img.AsSpan(24)) + 1);
            ResealHeader(img);
        });
        Assert.That(r.Failed("stats"), Is.True, Describe(r));
    }

    [Test]
    public void AWrongSymbolCount_IsFoundByTheTextCheck()
    {
        VerifyResult r = Damaged("GATTACAGATTACAGGTAGATTACA", "counts", (img, o) =>
        {
            // the COUNTS section sits in the table; move one count to the next letter (the sum still equals n, so the open at the header level accepts it)
            int count = BinaryPrimitives.ReadInt32LittleEndian(img.AsSpan(48));
            for (int i = 0; i < count; i++)
            {
                int e = 128 + (64 * i);
                if (BinaryPrimitives.ReadInt32LittleEndian(img.AsSpan(e)) != SesaFormat.KindCounts) continue;
                long off = BinaryPrimitives.ReadInt64LittleEndian(img.AsSpan(e + 8));
                BinaryPrimitives.WriteInt64LittleEndian(img.AsSpan((int)off), BinaryPrimitives.ReadInt64LittleEndian(img.AsSpan((int)off)) + 1);
                BinaryPrimitives.WriteInt64LittleEndian(img.AsSpan((int)off + 8), BinaryPrimitives.ReadInt64LittleEndian(img.AsSpan((int)off + 8)) - 1);
            }
        });
        Assert.That(r.Failed("text"), Is.True, Describe(r));
    }

    [Test]
    public void TheTextHash_OfTheHeader_IsComparedWithTheInput()
    {
        VerifyResult r = Damaged("GATTACAGATTACAGGTAGATTACA", "sha", (img, _) =>
        {
            img[80] ^= 0xFF;
            ResealHeader(img);
        });
        Assert.That(r.Failed("input"), Is.True, Describe(r));
    }

    // ---- the command

    [Test]
    public void TheCommand_VerifiesAContainer_ExitCodeZero_AndWritesTheRows()
    {
        string container = Build(SmallWithEscapes(), "cli");
        string rows = Path.Combine(Path.GetDirectoryName(container)!, "rows.tsv");
        string sums = Path.Combine(Path.GetDirectoryName(container)!, "sums.txt");
        int rc = Cli.Run(["e8-verify", "--container", container, "--input", Path.Combine(Path.GetDirectoryName(container)!, "in.txt"), "--quiet", "--rows", rows, "--sums", sums, "--threads", "2"]);
        Assert.That(rc, Is.EqualTo(0));
        string text = File.ReadAllText(rows);
        Assert.That(text, Does.Contain("verify_all").And.Contain("11 of 11 checks passed"));
        Assert.That(text, Does.Contain("the payload of the section of out.sesa that plays '.sa'"));
        Assert.That(File.ReadAllLines(sums).Select(l => l.Split('=')[0]), Is.EqualTo(new[] { "pk", "pk.meta", "sa", "lcp", "ovf", "cld", "stats" }));
    }

    [Test]
    public void TheCommand_ExitsOne_OnADamagedContainer_AndTwo_OnAFileThatIsNotOne()
    {
        string container = Build("GATTACAGATTACAGGTAGATTACA", "bad");
        (byte[] image, SesaOracleLike o) = Image(container);
        image[o.Sections[SesaFormat.KindSa].Offset + 20] ^= 0x04;
        File.WriteAllBytes(container, image);
        Assert.That(Cli.Run(["e8-verify", "--container", container, "--quiet", "--threads", "2"]), Is.EqualTo(1));
        string junk = Path.Combine(Path.GetDirectoryName(container)!, "junk.bin");
        File.WriteAllBytes(junk, new byte[9000]);
        Assert.That(Cli.Run(["e8-verify", "--container", junk, "--quiet"]), Is.EqualTo(2));
        Assert.That(Cli.Run(["e8-verify", "--container", Path.Combine(Path.GetDirectoryName(container)!, "missing.sesa"), "--quiet"]), Is.EqualTo(2));
        Assert.That(Cli.Run(["e8-verify", "--quiet"]), Is.EqualTo(2), "neither --prefix nor --container");
    }

    [Test]
    public void ThePrefixMode_IsUnchanged()
    {
        string prefix = FromReference(Run("pre"), "mississippi");
        VerifyResult r = CheckFiles(prefix);
        Assert.That(r.AllPassed, Is.True, Describe(r));
        Assert.That(Cli.Run(["e8-verify", "--prefix", prefix, "--quiet", "--threads", "2"]), Is.EqualTo(0));
    }

    private static VerifyResult CheckFiles(string prefix) => OutputVerifier.Run(VerifierHarness.Options(prefix));
}
