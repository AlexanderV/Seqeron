using System.Text;
using SuffixTree.Esa.External;

namespace SuffixTree.Esa.Tests.External;

/// <summary>
/// The pieces of the resumable build that stand alone (work package 1b.9e): the format of <c>state.json</c> (a state is only what its CRC and its shape say; a damaged, cut or edited file is refused, never an
/// exception), the build lock (managed: <see cref="FileShare.None"/>), and the reopening of the container of an interrupted build.
/// </summary>
[TestFixture]
public sealed class SesaResumeFormatTests
{
    private static SesaResumeState Sample() => new()
    {
        N = 12345,
        Bits = 2,
        Sigma = 4,
        Sha256 = new string('a', 32) + new string('0', 32),
        Letters = [65, 67, 71, 84],
        Counts = [3000, 3000, 3345, 3000],
        Psi = true,
        ChildCodec = "tagged",
        ReservedLength = 987_654,
        Stages =
        [
            new SesaStageRecord(SesaStages.Packed, [new SesaFileCheck("text", true, 4096, 3100, 0xDEADBEEF), new SesaFileCheck("alphabet", true, 8192, 8, 1), new SesaFileCheck("counts", true, 8200, 32, 2)], new Dictionary<string, long>()),
            new SesaStageRecord(SesaStages.Sa, [new SesaFileCheck("sa", true, 16384, 49384, 0x0BADF00D)], new Dictionary<string, long>()),
            new SesaStageRecord(SesaStages.Lcp, [new SesaFileCheck("child", true, 70000, 49384, 7), new SesaFileCheck("lcp", true, 200000, 12346, 8)], new Dictionary<string, long> { ["kind"] = 0, ["overflow"] = 0, ["final_length"] = 987_654, ["node_count"] = 24_000 }),
            new SesaStageRecord(SesaStages.Bwt, [new SesaFileCheck("ph6/bwt.bin", false, 0, 12346, 9)], new Dictionary<string, long> { ["primary_row"] = 4711 }),
        ],
    };

    private static string Json(byte[] file) => Encoding.UTF8.GetString(file);

    [Test]
    public void ASerialisedState_ParsesBackToItself_AndIsValidJsonWithTheCrcLast()
    {
        SesaResumeState s = Sample();
        byte[] file = SesaResumeFormat.Serialize(s);
        Assert.That(SesaResumeFormat.TryParse(file, out SesaResumeState? back), Is.True);
        Assert.That(back!.N, Is.EqualTo(s.N));
        Assert.That(back.Sha256, Is.EqualTo(s.Sha256));
        Assert.That(back.Letters, Is.EqualTo(s.Letters));
        Assert.That(back.Counts, Is.EqualTo(s.Counts));
        Assert.That((back.Psi, back.ChildCodec, back.ReservedLength, back.Writer, back.Version), Is.EqualTo((true, "tagged", 987_654L, SesaResumeFormat.WriterId, 1)));
        Assert.That(back.Stages.Select(x => x.Name), Is.EqualTo(new[] { "packed", "sa", "lcp", "bwt" }));
        Assert.That(back.Stages[0].Checks[0], Is.EqualTo(s.Stages[0].Checks[0]));
        Assert.That(back.Stages[2].Value("final_length"), Is.EqualTo(987_654L));
        Assert.That(back.Stages[3].Checks[0].InContainer, Is.False);
        Assert.That(SesaResumeFormat.Serialize(back), Is.EqualTo(file), "serialising again gives the same bytes");
        using var doc = System.Text.Json.JsonDocument.Parse(file);
        string last = doc.RootElement.EnumerateObject().Last().Name;
        Assert.That(last, Is.EqualTo("crc32c"), "the CRC is the last member");
        Assert.That(Json(file), Does.EndWith("\"}\n"));
    }

    [Test]
    public void EveryCutAndEverySingleByteFlipOfAStateFile_IsRefused()
    {
        byte[] file = SesaResumeFormat.Serialize(Sample());
        for (int cut = 0; cut < file.Length - 1; cut++)
        {
            Assert.That(SesaResumeFormat.TryParse(file.AsSpan(0, cut), out SesaResumeState? s), Is.False, $"a file cut to {cut} bytes");
            Assert.That(s, Is.Null);
        }
        Assert.That(SesaResumeFormat.TryParse(file.AsSpan(0, file.Length - 1), out _), Is.True, "the newline at the end is not data");
        for (int at = 0; at < file.Length; at++)
        {
            for (int bit = 0; bit < 8; bit++)
            {
                byte[] bad = (byte[])file.Clone();
                bad[at] ^= (byte)(1 << bit);
                Assert.That(SesaResumeFormat.TryParse(bad, out _), Is.False, $"byte {at} bit {bit} flipped");
            }
        }
        Assert.That(SesaResumeFormat.TryParse(file.Concat(new byte[] { 0 }).ToArray(), out _), Is.False, "a byte appended");
        Assert.That(SesaResumeFormat.TryParse(new byte[SesaResumeFormat.MaxBytes + 1], out _), Is.False, "too long");
        Assert.That(SesaResumeFormat.TryParse([], out _), Is.False, "empty");
    }

    [Test]
    public void AFileThatIsNotWhatTheWriterProduces_IsRefusedEvenWithAValidCrc()
    {
        SesaResumeState s = Sample();
        string good = Json(SesaResumeFormat.Serialize(s));
        string body = good[..good.LastIndexOf(",\"crc32c\"", StringComparison.Ordinal)];
        byte[] Sealed(string b) => SesaResumeFormat.Seal(Encoding.UTF8.GetBytes(b));
        Assert.That(SesaResumeFormat.TryParse(Sealed(body), out _), Is.True, "the helper reproduces the writer's own file");

        var edits = new (string Name, string From, string To)[]
        {
            ("another version", "\"version\":1", "\"version\":2"),
            ("another format", "\"sesa-resume\"", "\"sesa-other\""),
            ("negative n", "\"n\":12345", "\"n\":-1"),
            ("a short hash", new string('a', 32), "aa"),
            ("an upper-case hash", new string('a', 32), new string('A', 32)),
            ("sigma and letters differ", "\"sigma\":4", "\"sigma\":5"),
            ("a letter above 255", "\"letters\":[65,", "\"letters\":[300,"),
            ("a negative count", "\"counts\":[3000,", "\"counts\":[-3,"),
            ("an unknown codec", "\"child_codec\":\"tagged\"", "\"child_codec\":\"zip\""),
            ("an unknown stage", "\"name\":\"sa\"", "\"name\":\"zz\""),
            ("stages out of order", "\"name\":\"sa\"", "\"name\":\"bwt\""),
            ("a repeated stage", "\"name\":\"lcp\"", "\"name\":\"sa\""),
            ("a negative offset", "\"offset\":4096", "\"offset\":-4096"),
            ("a negative length", "\"length\":3100", "\"length\":-3100"),
            ("a scratch file above the scratch directory", "\"role\":\"ph6/bwt.bin\"", "\"role\":\"../x\""),
            ("a scratch file at an absolute path", "\"role\":\"ph6/bwt.bin\"", "\"role\":\"/etc/x\""),
            ("a scratch file with an offset", "\"container\":false,\"offset\":0", "\"container\":false,\"offset\":5"),
            ("a crc that is not hex", "\"crc32c\":\"deadbeef\"", "\"crc32c\":\"zzzzzzzz\""),
            ("a missing member", "\"psi\":true,", string.Empty),
            ("a string where a number is", "\"bits\":2", "\"bits\":\"2\""),
        };
        foreach ((string name, string from, string to) in edits)
        {
            Assert.That(body, Does.Contain(from), name);
            Assert.That(SesaResumeFormat.TryParse(Sealed(body.Replace(from, to, StringComparison.Ordinal)), out SesaResumeState? parsed), Is.False, name);
            Assert.That(parsed, Is.Null, name);
        }
        Assert.That(SesaResumeFormat.TryParse(Sealed("[1,2"), out _), Is.False, "an array");
        Assert.That(SesaResumeFormat.TryParse(Sealed("{\"format\":\"sesa-resume\""), out _), Is.False, "an object without members");
        Assert.That(SesaResumeFormat.TryParse(Sealed("garbage"), out _), Is.False, "text that is not JSON");
    }

    [Test]
    public void TheStagesThatCannotBeTogether_AreRefused()
    {
        SesaResumeState s = Sample();
        SesaStageRecord Rec(string name) => new(name, [], new Dictionary<string, long>());
        Assert.That(SesaResumeFormat.TryParse(SesaResumeFormat.Serialize(s with { Stages = [Rec("packed"), Rec("ranks")] }), out _), Is.True);
        Assert.That(SesaResumeFormat.TryParse(SesaResumeFormat.Serialize(s with { Stages = [Rec("ranks")] }), out _), Is.False, "a stage without the first one");
        Assert.That(SesaResumeFormat.TryParse(SesaResumeFormat.Serialize(s with { Stages = [Rec("packed"), Rec("lcp")] }), out _), Is.False, "LCP without SA");
        Assert.That(SesaResumeFormat.TryParse(SesaResumeFormat.Serialize(s with { Stages = [Rec("packed"), Rec("sa"), Rec("bwt")] }), out _), Is.False, "BWT without LCP");
        Assert.That(SesaResumeFormat.TryParse(SesaResumeFormat.Serialize(s with { Stages = [Rec("packed"), Rec("sa"), Rec("lcp"), Rec("psi")], Psi = false }), out _), Is.False, "PSI in a build without PSI");
        Assert.That(SesaResumeFormat.TryParse(SesaResumeFormat.Serialize(s with { Stages = [] }), out _), Is.True, "no stage at all is a state that holds nothing");
    }

    [Test]
    public void TheWriterIdOfTheState_IsTheWriterOfTheBuildInfo()
    {
        PackedTextInfo info = new(4, 1, [97, 98], [2, 2], new string('0', 64));
        Assert.That(SesaBuilder.BuildInfoText(info, LcpKind.U8, 0, true, ChildCodec.Tagged), Does.Contain("writer=" + SesaResumeFormat.WriterId + "\n"));
    }

    [Test]
    public void WithAndTruncate_KeepTheOrder_AndASupersededStageGoes()
    {
        SesaStageRecord Rec(string name) => new(name, [], new Dictionary<string, long>());
        SesaResumeState s = Sample() with { Stages = [] };
        s = s.With(Rec("packed")).With(Rec("ranks"));
        Assert.That(s.Stages.Select(x => x.Name), Is.EqualTo(new[] { "packed", "ranks" }));
        s = s.With(Rec("sa"));
        Assert.That(s.Stages.Select(x => x.Name), Is.EqualTo(new[] { "packed", "sa" }), "SA supersedes the rank file");
        s = s.With(Rec("lcp")).With(Rec("bwt"));
        Assert.That(s.Stages.Select(x => x.Name), Is.EqualTo(new[] { "packed", "sa", "lcp", "bwt" }));
        s = s.With(Rec("psi"));
        Assert.That(s.Stages.Select(x => x.Name), Is.EqualTo(new[] { "packed", "sa", "lcp", "psi" }), "PSI supersedes the BWT file");
        Assert.That(s.Truncate(2).Stages.Select(x => x.Name), Is.EqualTo(new[] { "packed", "sa" }));
        Assert.That(s.With(Rec("sa")).Stages.Select(x => x.Name), Is.EqualTo(new[] { "packed", "sa", "lcp", "psi" }), "a stage that is recorded again replaces itself");
        Assert.That(s.Find("lcp"), Is.Not.Null);
        Assert.That(s.Find("bwt"), Is.Null);
    }

    // ------------------------------------------------------------------ the lock

    [Test]
    public void TheLock_ExcludesASecondHolder_UntilTheFirstIsDone_AndRemovesItsFile()
    {
        using var dir = new TempDir();
        string path = dir.File("x.sesa.lock");
        SesaBuildLock first = SesaBuildLock.Acquire(path);
        Assert.That(File.Exists(path), Is.True);
        IOException? e = Assert.Throws<IOException>(() => SesaBuildLock.Acquire(path));
        Assert.That(e!.Message, Does.Contain("lock"));
        Assert.That(File.Exists(path), Is.True, "the refused one did not remove the holder's file");
        first.Dispose();
        first.Dispose();
        Assert.That(File.Exists(path), Is.False, "the holder removes the file");
        using SesaBuildLock second = SesaBuildLock.Acquire(path);
        Assert.That(File.Exists(path), Is.True);
    }

    [Test]
    public void ALockFileThatAKilledProcessLeft_IsOnlyAName()
    {
        using var dir = new TempDir();
        string path = dir.File("x.sesa.lock");
        File.WriteAllText(path, "left by a build that was killed");
        using (SesaBuildLock.Acquire(path)) Assert.That(File.Exists(path), Is.True);
        Assert.That(File.Exists(path), Is.False);
    }

    [Test]
    public void TheLock_InADirectoryThatIsGone_IsAnIoExceptionWithoutAClaimThatSomeoneElseHoldsIt()
    {
        using var dir = new TempDir();
        Exception? e = Assert.Catch(() => SesaBuildLock.Acquire(Path.Combine(dir.File("missing"), "x.lock")));
        Assert.That(e, Is.InstanceOf<DirectoryNotFoundException>());
    }

    // ------------------------------------------------------------------ reopening the container

    [Test]
    public void OpenExisting_KeepsWhatTheFileHolds_ReportsReopened_AndWritesThrough()
    {
        using var dir = new TempDir();
        string path = dir.File("c.bin");
        File.WriteAllBytes(path, Enumerable.Range(0, 10_000).Select(i => (byte)(i % 251)).ToArray());
        using (MappedOutputFile file = MappedOutputFile.OpenExisting(path))
        {
            Assert.That((file.Length, file.Outcome, file.Path), Is.EqualTo((10_000L, ReservationOutcome.Reopened, path)));
            using MappedWindowWriter w = file.CreateWriter(100, 50, 16);
            w.Append(Enumerable.Repeat((byte)0xEE, 50).ToArray());
        }
        byte[] bytes = File.ReadAllBytes(path);
        Assert.That(bytes.Length, Is.EqualTo(10_000));
        Assert.That(bytes[99], Is.EqualTo((byte)(99 % 251)));
        Assert.That(bytes[100], Is.EqualTo((byte)0xEE));
        Assert.That(bytes[149], Is.EqualTo((byte)0xEE));
        Assert.That(bytes[150], Is.EqualTo((byte)(150 % 251)), "everything else is as it was");
    }

    [Test]
    public void OpenExisting_CutsALongerFileBack_AndLeavesAShorterOne()
    {
        using var dir = new TempDir();
        string path = dir.File("c.bin");
        File.WriteAllBytes(path, new byte[10_000]);
        using (MappedOutputFile file = MappedOutputFile.OpenExisting(path, 6_000)) Assert.That(file.Length, Is.EqualTo(6_000));
        Assert.That(new FileInfo(path).Length, Is.EqualTo(6_000));
        using (MappedOutputFile file = MappedOutputFile.OpenExisting(path, 8_000)) Assert.That(file.Length, Is.EqualTo(6_000), "a file shorter than the cut is not extended");
        using (MappedOutputFile file = MappedOutputFile.OpenExisting(path, 6_000)) Assert.That(file.Length, Is.EqualTo(6_000), "a file of exactly the length is not touched");
        Assert.Throws<ArgumentOutOfRangeException>(() => MappedOutputFile.OpenExisting(path, -1));
    }

    [Test]
    public void OpenExisting_OfAFileThatIsGone_IsFileNotFound_AndCreatesNothing()
    {
        using var dir = new TempDir();
        string path = dir.File("gone.bin");
        Assert.Throws<FileNotFoundException>(() => MappedOutputFile.OpenExisting(path));
        Assert.Throws<FileNotFoundException>(() => MappedOutputFile.OpenExisting(path, 100));
        Assert.That(File.Exists(path), Is.False);
    }

    [Test]
    public void AReopenedContainer_GrowsLikeAFreshOne_AndTheGrowthIsReservedOnce()
    {
        using var dir = new TempDir();
        string path = dir.File("c.bin");
        File.WriteAllBytes(path, new byte[8192]);
        var counters = new IoCounters();
        using MappedOutputFile file = MappedOutputFile.OpenExisting(path);
        Assert.That(file.Grow(8192 + 5000, null, counters), Is.EqualTo(ReservationOutcome.Filled));
        Assert.That((file.Length, new FileInfo(path).Length, counters.ReservedWriteBytes), Is.EqualTo((13192L, 13192L, 5000L)));
    }
}
