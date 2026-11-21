using System.Collections.Concurrent;
using System.Security.Cryptography;
using Microsoft.Win32.SafeHandles;
using SuffixTree.Esa.External;
using SuffixTree.Esa.Tests.TestSupport;

namespace SuffixTree.Esa.Tests.External;

/// <summary>
/// Work package 1b.9e.1: a resumed build trusts only what a kept, CRC-verified stage vouches for. Every byte of the temporary container that lies inside its length and outside a kept range (the header page, the padding between the sections, the room of
/// BUILD_INFO behind its text, the slots of CHILD that the scan never writes, the whole part that belongs to the stages that are not kept, the part behind the reserved length) must be exactly what a new file has there (zero) when the
/// build goes on, whatever another process, a failing disk or a power loss put there: the finished container is then the container of an uninterrupted build.
/// </summary>
[TestFixture]
public sealed class SesaResume1b9e1Tests
{
    private static readonly string s_text = StringFamilies.RandomDna(1500, 11) + new string('A', 300) + StringFamilies.RandomDna(900, 12);

    internal static readonly Dictionary<string, string> Texts = new()
    {
        ["dna"] = s_text,
        ["u16"] = StringFamilies.RandomDna(30_000, 5) + new string('A', 2_600),
        ["overflow"] = StringFamilies.RandomDna(50_000, 6) + new string('A', 700),
        ["u32"] = new string('a', 100_000),
    };

    /// <summary>The boundaries at which a resumable build has a state (the first three boundaries are before the first record).</summary>
    internal static readonly string[] Boundaries = [.. SesaHooks.All.Skip(3)];

    private sealed record Built(byte[] Bytes, SesaBuildResult Result);

    private static readonly ConcurrentDictionary<string, Built> s_built = new();

    private static byte[] Read(string path)
    {
        using var fs = new FileStream(path, FileMode.Open, FileAccess.Read, FileShare.ReadWrite | FileShare.Delete);
        var bytes = new byte[fs.Length];
        fs.ReadExactly(bytes);
        return bytes;
    }

    /// <summary>An uninterrupted build without Resume: the bytes and the result (the oracle of every test here).</summary>
    private static Built Reference(string text, bool psi = true) => s_built.GetOrAdd($"{text.GetHashCode(StringComparison.Ordinal)}-{text.Length}-{psi}", _ =>
    {
        using var c = new ResumeCase(text);
        SesaBuildResult r = SesaBuilder.Build(c.Input, c.Target, new SesaBuildOptions { Threads = 2, IncludePsi = psi, Fsync = false, MeasureRss = false });
        return new Built(Read(c.Target), r);
    });

    private static SesaBuildOptions Opts(bool psi = true, int threads = 2) => SesaRun.Tiny(3, threads, psi) with { Resume = true, MeasureRss = false };

    /// <summary>The tiny parameters for the short text (every multi-bucket and multi-chunk path), the defaults for the long ones (which need a wider LCP).</summary>
    private static SesaBuildOptions OptsFor(string kind, int threads = 2) => kind == "dna" ? Opts(threads: threads) : new SesaBuildOptions { Threads = threads, Fsync = false, MeasureRss = false, Resume = true };

    private static void Fails(ResumeCase c, SesaBuildOptions options, string point)
    {
        int seen = 0;
        Assert.Throws<InjectedFailureException>(() => SesaBuilder.Build(c.Input, c.Target, options with { Hook = p => { if (p == point && ++seen == 1) throw new InjectedFailureException(p); } }));
    }

    private static void AssertSame(ResumeCase c, byte[] reference, string label)
    {
        byte[] got = Read(c.Target);
        Assert.That(got.Length, Is.EqualTo(reference.Length), $"{label}: length");
        Assert.That(SHA256.HashData(got), Is.EqualTo(SHA256.HashData(reference)), $"{label}: the bytes of an uninterrupted build");
    }

    private static void AssertOnlyTarget(ResumeCase c, string label) => Assert.That(c.Entries, Is.EqualTo(new[] { "x.sesa" }), $"{label}: nothing but the target is left");

    /// <summary>The ranges of the container that the stages of the state vouch for, sorted.</summary>
    private static List<(long Offset, long Length)> KeptRanges(ResumeCase c)
    {
        Assert.That(SesaResumeFormat.TryParse(File.ReadAllBytes(c.State), out SesaResumeState? s), Is.True, "the state of the failed run");
        return [.. s!.Stages.SelectMany(x => x.Checks).Where(k => k.InContainer).Select(k => (k.Offset, k.Length)).Order()];
    }

    /// <summary>What is left of <c>[0, length)</c> when the ranges are taken away, as spans.</summary>
    private static List<(long Offset, long Length)> Complement(IReadOnlyList<(long Offset, long Length)> ranges, long length)
    {
        var gaps = new List<(long, long)>();
        long at = 0;
        foreach ((long offset, long len) in ranges.OrderBy(r => r.Offset))
        {
            if (offset > at) gaps.Add((at, Math.Min(offset, length) - at));
            at = Math.Max(at, offset + len);
            if (at >= length) break;
        }
        if (at < length) gaps.Add((at, length - at));
        return [.. gaps.Where(g => g.Item2 > 0)];
    }

    /// <summary>The bytes of the finished container that no section has: the header page, the padding between the sections, the room of BUILD_INFO behind its text.</summary>
    private static List<(long Offset, long Length)> Gaps(Built b) => Complement([.. b.Result.Sections.Select(s => (s.Offset, s.Length))], b.Bytes.Length);

    private static byte Garbage(long offset, int pattern) => pattern == 0 ? (byte)0xA5 : (byte)((((offset * 2654435761L) >> 11) & 0xFF) | 1);

    private static void Fill(string path, IEnumerable<(long Offset, long Length)> spans, int pattern)
    {
        using var fs = new FileStream(path, FileMode.Open, FileAccess.ReadWrite, FileShare.ReadWrite);
        long length = fs.Length;
        foreach ((long offset, long len) in spans)
        {
            long end = Math.Min(offset + len, length);
            if (end <= offset) continue;
            var block = new byte[end - offset];
            for (long i = 0; i < block.Length; i++) block[i] = Garbage(offset + i, pattern);
            fs.Position = offset;
            fs.Write(block);
        }
    }

    /// <summary>Changes one byte (never to the value it has).</summary>
    private static void Damage(string path, long offset)
    {
        using var fs = new FileStream(path, FileMode.Open, FileAccess.ReadWrite, FileShare.ReadWrite);
        if (offset >= fs.Length) return;
        fs.Position = offset;
        int b = fs.ReadByte();
        fs.Position = offset;
        fs.WriteByte(b == 0xA5 ? (byte)0x5A : (byte)0xA5);
    }

    private static string[] DoneAt(string point, bool psi) => SesaResume1b9eTests.DoneAt(point, psi);

    // ------------------------------------------------------------------ the two cases of the orchestrator's verification, as tests

    [Test]
    public void AFlippedLastByteOfTheRecordedChildRange_IsNotInTheContainer_TheLcpStageIsDoneAgain()
    {
        // the slot is one that the child writer never writes: the CRC of the recorded range catches the flip, the stage is redone, and the redone stage must not leave the flip behind
        Built reference = Reference(s_text);
        using var c = new ResumeCase(s_text);
        Fails(c, Opts(), SesaHooks.LcpChildBuilt);
        SesaPlan plan = SesaLayout.Plan(PackedTextWriter.Scan(c.Input), true);
        long last = plan.ChildOffset + plan.ChildLength - 1;
        Assert.That(reference.Bytes[last], Is.EqualTo(0), "the last byte of CHILD is a slot that the scan does not set");
        Damage(c.Tmp, last);
        SesaBuildResult r = SesaBuilder.Build(c.Input, c.Target, Opts());
        Assert.That(r.Resume!.Skipped, Is.EqualTo(new[] { "packed", "sa" }), "the damaged stage is not kept");
        AssertSame(c, reference.Bytes, "child last byte");
        AssertOnlyTarget(c, "child last byte");
    }

    [Test]
    public void GarbageFromTheEndOfTheLastRecordedRangeToTheEndOfTheFile_IsNotInTheContainer([Values(SesaHooks.Packed, SesaHooks.SaBuilt, SesaHooks.LcpChildBuilt, SesaHooks.PsiBuilt)] string point, [Values(0, 1)] int pattern)
    {
        Built reference = Reference(s_text);
        using var c = new ResumeCase(s_text);
        Fails(c, Opts(), point);
        List<(long Offset, long Length)> kept = KeptRanges(c);
        long end = kept.Max(k => k.Offset + k.Length);
        long length = new FileInfo(c.Tmp).Length;
        Fill(c.Tmp, [(end, length - end)], pattern);
        SesaBuildResult r = SesaBuilder.Build(c.Input, c.Target, Opts());
        Assert.That(r.Resume!.Skipped, Is.EqualTo(DoneAt(point, true)));
        AssertSame(c, reference.Bytes, $"{point} {pattern}");
        AssertOnlyTarget(c, point);
    }

    // ------------------------------------------------------------------ everything outside a kept range, at every boundary

    [Test]
    public void GarbageEverywhereOutsideTheKeptRanges_IsNotInTheContainer_AtEveryBoundary_OnEveryKindOfTail(
        [ValueSource(nameof(Boundaries))] string point, [Values(0, 1)] int pattern, [Values("dna", "u16", "overflow", "u32")] string kind)
    {
        string text = Texts[kind];
        Built reference = Reference(text);
        using var c = new ResumeCase(text);
        Fails(c, OptsFor(kind), point);
        List<(long Offset, long Length)> outside = Complement(KeptRanges(c), new FileInfo(c.Tmp).Length);
        Assert.That(outside, Is.Not.Empty, "there is something to damage");
        Fill(c.Tmp, outside, pattern);
        SesaBuildResult r = SesaBuilder.Build(c.Input, c.Target, OptsFor(kind));
        Assert.That(r.Resume!.Skipped, Is.EqualTo(DoneAt(point, true)), "garbage outside the kept ranges does not make a kept stage fail its CRC");
        AssertSame(c, reference.Bytes, $"{kind} {point} {pattern}");
        AssertOnlyTarget(c, point);
    }

    [Test]
    public void GarbageEverywhereOutsideTheKeptRanges_AlsoWithoutPsi([Values(SesaHooks.Packed, SesaHooks.SaBuilt, SesaHooks.LcpChildBuilt, SesaHooks.InfoWritten, SesaHooks.MagicWritten)] string point, [Values(0, 1)] int pattern)
    {
        Built reference = Reference(s_text, psi: false);
        using var c = new ResumeCase(s_text);
        Fails(c, Opts(psi: false), point);
        Fill(c.Tmp, Complement(KeptRanges(c), new FileInfo(c.Tmp).Length), pattern);
        SesaBuildResult r = SesaBuilder.Build(c.Input, c.Target, Opts(psi: false));
        Assert.That(r.Resume!.Skipped, Is.EqualTo(DoneAt(point, false)));
        AssertSame(c, reference.Bytes, $"no psi {point} {pattern}");
        AssertOnlyTarget(c, point);
    }

    [Test]
    public void GarbageInAContainerOfTheLargeAlignment_IsNotInTheContainer([Values(SesaHooks.Packed, SesaHooks.SaBuilt, SesaHooks.LcpChildBuilt, SesaHooks.MagicWritten)] string point)
    {
        // from 1 MiB a section starts at a multiple of 4096: the padding is up to 4,095 bytes (the 2,464 bytes between COUNTS and SA of the orchestrator's case at 10^7)
        string text = StringFamilies.RandomDna(300_000, 21);
        Built reference = Reference(text);
        Assert.That(Gaps(reference).Max(g => g.Length), Is.GreaterThan(64), "a gap of the large alignment exists");
        using var c = new ResumeCase(text);
        SesaBuildOptions o = new() { Threads = 2, Fsync = false, MeasureRss = false, Resume = true };
        Fails(c, o, point);
        Fill(c.Tmp, Complement(KeptRanges(c), new FileInfo(c.Tmp).Length), 0);
        SesaBuildResult r = SesaBuilder.Build(c.Input, c.Target, o);
        Assert.That(r.Resume!.Skipped, Is.EqualTo(DoneAt(point, true)));
        AssertSame(c, reference.Bytes, point);
        AssertOnlyTarget(c, point);
    }

    // ------------------------------------------------------------------ one bad byte in a gap, in the header page, in a slot that is never written

    public static IEnumerable<TestCaseData> GapCases()
    {
        Built reference = Reference(s_text);
        List<(long Offset, long Length)> gaps = Gaps(reference);
        foreach (string point in Boundaries)
        {
            for (int g = 0; g < gaps.Count; g++)
            {
                foreach (string where in new[] { "first", "last" })
                {
                    long at = where == "first" ? gaps[g].Offset : gaps[g].Offset + gaps[g].Length - 1;
                    yield return new TestCaseData(point, at).SetName($"OneBadByteInAGapOrTheHeaderPage({point}, gap {g} [{gaps[g].Offset}, +{gaps[g].Length}) {where})");
                }
            }
        }
    }

    [TestCaseSource(nameof(GapCases))]
    public void OneBadByteInAGapOfTheLayout_OrInTheHeaderPage_IsNotInTheContainer(string point, long at)
    {
        Built reference = Reference(s_text);
        using var c = new ResumeCase(s_text);
        Fails(c, Opts(), point);
        Assert.That(new FileInfo(c.Tmp).Length, Is.GreaterThan(0));
        Damage(c.Tmp, at);
        SesaBuildResult r = SesaBuilder.Build(c.Input, c.Target, Opts());
        Assert.That(r.Resume!.Skipped, Is.EqualTo(DoneAt(point, true)), "a byte outside every recorded range does not make a kept stage fail");
        AssertSame(c, reference.Bytes, $"{point} @{at}");
        AssertOnlyTarget(c, point);
    }

    [Test]
    public void ADamagedSlotOfChildThatTheScanNeverWrites_IsZeroAgain_AtEveryBoundary_OnEveryKindOfTail(
        [ValueSource(nameof(Boundaries))] string point, [Values("dna", "u16", "overflow")] string kind, [Values("last", "first")] string which)
    {
        string text = Texts[kind];
        Built reference = Reference(text);
        using var c = new ResumeCase(text);
        Fails(c, OptsFor(kind), point);
        SesaPlan plan = SesaLayout.Plan(PackedTextWriter.Scan(c.Input), true);
        // a word of CHILD that is zero in the finished container is a slot that the scan does not set
        IEnumerable<long> words = Enumerable.Range(0, (int)(plan.ChildLength / 4)).Select(i => plan.ChildOffset + (4L * i));
        long slot = (which == "last" ? words.Reverse() : words).First(w => reference.Bytes.AsSpan((int)w, 4).SequenceEqual(new byte[4]));
        Fill(c.Tmp, [(slot, 4)], 0);
        SesaBuildResult r = SesaBuilder.Build(c.Input, c.Target, OptsFor(kind));
        AssertSame(c, reference.Bytes, $"{kind} {point} slot {slot}");
        AssertOnlyTarget(c, point);
        string[] done = DoneAt(point, true);
        Assert.That(r.Resume!.Skipped, Is.EqualTo(done.Contains(SesaStages.Lcp) ? done.Where(s => s is SesaStages.Packed or SesaStages.Sa).ToArray() : done), "a damaged recorded CHILD fails its CRC, the later stages go with it");
    }

    // ------------------------------------------------------------------ a resume that is stopped again, then resumed

    [Test]
    public void AResumeThatIsStoppedRightAfterTheNormalisation_AndResumedAgain_IsByteIdentical([ValueSource(nameof(Boundaries))] string point, [Values(SesaHooks.Created, SesaHooks.SaBuilt)] string second)
    {
        // the second boundary is reached whether the stage ran or was found done: the damage is put right (or not yet) when the second run dies, and the third run must still finish with the bytes of an uninterrupted build
        Built reference = Reference(s_text);
        using var c = new ResumeCase(s_text);
        Fails(c, Opts(), point);
        Fill(c.Tmp, Complement(KeptRanges(c), new FileInfo(c.Tmp).Length), 0);
        Fails(c, Opts(threads: 1), second);
        Fill(c.Tmp, Complement(KeptRanges(c), new FileInfo(c.Tmp).Length), 1); // the second run is damaged again after its death
        SesaBuildResult r = SesaBuilder.Build(c.Input, c.Target, Opts(threads: 4));
        Assert.That(r.Resume!.Accepted, Is.True);
        AssertSame(c, reference.Bytes, $"{point} then {second}");
        AssertOnlyTarget(c, point);
    }

    // ------------------------------------------------------------------ what the build reports

    [Test]
    public void TheReport_SaysWhatWasRead_AndWhatWasWritten([Values(SesaHooks.Packed, SesaHooks.SaBuilt, SesaHooks.LcpChildBuilt)] string point, [Values(false, true)] bool damaged)
    {
        // at these boundaries a killed run has written nothing outside the ranges of its kept stages, so a clean file has nothing to clear and one damaged byte costs one unit
        Built reference = Reference(s_text);
        using var c = new ResumeCase(s_text);
        Fails(c, Opts(), point);
        long length = new FileInfo(c.Tmp).Length;
        List<(long Offset, long Length)> spans = Complement(KeptRanges(c), length);
        long outside = spans.Sum(g => g.Length);
        if (damaged)
        {
            Damage(c.Tmp, spans[0].Offset);
            Damage(c.Tmp, spans[^1].Offset + spans[^1].Length - 1);
        }
        SesaBuildResult r = SesaBuilder.Build(c.Input, c.Target, Opts());
        Assert.That(r.Resume!.ClearScannedBytes, Is.EqualTo(outside), "everything outside the kept ranges is read once");
        if (damaged) Assert.That(r.Resume.ClearedBytes, Is.InRange(2, 2 * ZeroFill.UnitBytes), "and the two units that held a damaged byte are written");
        else Assert.That(r.Resume.ClearedBytes, Is.EqualTo(0), "a file with nothing wrong in it is not written");
        Assert.That(r.Resume.ClearSeconds, Is.GreaterThan(0));
        AssertSame(c, reference.Bytes, point);
    }

    [Test]
    public void TheKeptRanges_AreThoseOfTheStagesThatVerified_AndOfTheContainerOnly([Values(SesaHooks.PsiBuilt, SesaHooks.BwtBuilt, SesaHooks.SaBuilt)] string point, [Values(false, true)] bool damageSa)
    {
        // the ranges that the normalisation keeps are the ranges of the container that the stages which passed their CRC vouch for: not those of a stage that the state names and that failed, and not the scratch files of a kept stage
        using var c = new ResumeCase(s_text);
        Fails(c, Opts(), point);
        List<(long Offset, long Length)> named = KeptRanges(c);
        PackedTextInfo info = PackedTextWriter.Scan(c.Input);
        SesaPlan plan = SesaLayout.Plan(info, true);
        if (damageSa) Damage(c.Tmp, plan.SaOffset + 100);
        SesaResumeRun run = SesaResumeRun.Load(Opts(), c.Tmp, c.Scratch, info, plan, SesaBuilder.ChildCodecFor(info.N, null), new IoCounters());
        if (!damageSa) Assert.That(run.KeptContainerRanges().Order(), Is.EqualTo(named), "every range of the container that the kept stages name, and nothing else");
        else
        {
            Assert.That(run.Skipped, Is.EqualTo(new[] { "packed" }));
            Assert.That(run.KeptContainerRanges().Order(), Is.EqualTo(new[] { (plan.TextOffset, plan.TextLength), (plan.AlphabetOffset, plan.AlphabetLength), (plan.CountsOffset, plan.CountsLength) }.Order()));
        }
    }

    [Test]
    public void AFreshBuildWithResume_ClearsNothing_AndTheDefaultBuildHasNoReport()
    {
        using var c = new ResumeCase(s_text);
        SesaBuildResult fresh = SesaBuilder.Build(c.Input, c.Target, Opts());
        Assert.That((fresh.Resume!.ClearScannedBytes, fresh.Resume.ClearedBytes, fresh.Resume.ClearSeconds), Is.EqualTo((0L, 0L, 0.0)));
        File.Delete(c.Target);
        SesaBuildResult plain = SesaBuilder.Build(c.Input, c.Target, Opts() with { Resume = false });
        Assert.That(plain.Resume, Is.Null);
    }

    [Test]
    public void TheNormalisation_WritesInPlace_TheReservedZerosAreOnlyThoseOfTheGrowth_AndEveryUnitOfGarbageIsWritten()
    {
        // the file keeps its length and its blocks: the normalisation writes in place (the reservation of the file is what the creation made it)
        // and it is not the reservation of a growth: IoCounters.ReservedWriteBytes counts what PH4 reserves, nothing else
        Built reference = Reference(s_text);
        using var c = new ResumeCase(s_text);
        Fails(c, Opts(), SesaHooks.SaBuilt);
        long length = new FileInfo(c.Tmp).Length;
        long kept = KeptRanges(c).Sum(k => k.Length);
        Fill(c.Tmp, Complement(KeptRanges(c), length), 0);
        var counters = new IoCounters();
        SesaBuildResult r = SesaBuilder.Build(c.Input, c.Target, Opts() with { Counters = counters });
        Assert.That(r.Resume!.ClearScannedBytes, Is.EqualTo(length - kept));
        Assert.That(r.Resume.ClearedBytes, Is.EqualTo(length - kept), "every unit was garbage");
        Assert.That(counters.ReservedWriteBytes, Is.EqualTo(r.GrownBytes), "the zeros of the reservation are the growth's own: the normalisation is not counted as a reservation");
        Assert.That(r.GrownBytes, Is.EqualTo(reference.Result.FileSize - length), "the growth is the one of PH4");
        AssertSame(c, reference.Bytes, "tail");
    }
}

/// <summary>The normalisation of a reopened container by itself (work package 1b.9e.1): <see cref="ZeroFill.ClearOutside"/> and the <c>keep</c> ranges of <see cref="MappedOutputFile.OpenExisting"/>, on files of a few kilobytes to a few megabytes.</summary>
[TestFixture]
public sealed class ZeroFillClearingTests
{
    /// <summary>A file of <paramref name="length"/> bytes of which about half the 4,096-byte units hold a non-zero byte, the rest are zero.</summary>
    private static byte[] RandomFile(string path, int length, int seed)
    {
        var rnd = new Random(seed);
        var bytes = new byte[length];
        for (int unit = 0; unit < length; unit += ZeroFill.UnitBytes)
        {
            if (rnd.Next(2) == 0) continue;
            int end = Math.Min(length, unit + ZeroFill.UnitBytes);
            for (int i = 0; i < 1 + rnd.Next(3); i++) bytes[rnd.Next(unit, end)] = (byte)(1 + rnd.Next(255));
        }
        File.WriteAllBytes(path, bytes);
        return bytes;
    }

    private static ClearResult Clear(string path, long length, IReadOnlyList<(long Offset, long Length)> keep)
    {
        using SafeFileHandle handle = File.OpenHandle(path, FileMode.Open, FileAccess.ReadWrite, FileShare.ReadWrite);
        return ZeroFill.ClearOutside(handle, length, keep);
    }

    [Test]
    public void ClearOutside_MakesEverythingOutsideTheKeptRangesZero_AndTouchesNothingInsideThem([Range(0, 29)] int seed)
    {
        using var dir = new TempDir();
        string path = dir.File("f.bin");
        var rnd = new Random(seed);
        int length = seed % 3 == 0 ? rnd.Next(1, 5000) : rnd.Next(ZeroFill.ChunkBytes - 3000, (2 * ZeroFill.ChunkBytes) + 90_000);
        byte[] before = RandomFile(path, length, seed);
        var keep = new List<(long Offset, long Length)>();
        for (int i = 0, count = rnd.Next(0, 6); i < count; i++)
        {
            long offset = rnd.Next(-10, length + 20); // ranges may start before the file, overlap one another, touch, nest, or reach beyond the end
            offset = Math.Max(0, offset);
            keep.Add((offset, rnd.Next(1, (length / 3) + 2)));
        }
        var inside = new bool[length];
        foreach ((long offset, long len) in keep)
        {
            for (long i = offset; i < Math.Min(length, offset + len); i++) inside[i] = true;
        }
        ClearResult result = Clear(path, length, keep);
        byte[] after = File.ReadAllBytes(path);
        Assert.That(after.Length, Is.EqualTo(length), "the length is not changed");
        long outsideBytes = 0, nonZeroOutside = 0, expectedCleared = 0;
        for (int i = 0; i < length; i++)
        {
            if (inside[i]) Assert.That(after[i], Is.EqualTo(before[i]), $"byte {i} is inside a kept range and is not touched");
            else
            {
                Assert.That(after[i], Is.EqualTo(0), $"byte {i} is outside every kept range and is zero");
                outsideBytes++;
                if (before[i] != 0) nonZeroOutside++;
            }
        }
        // the accounting, from the definition: a gap is a maximal run of bytes outside every range, it is read in chunks, a unit of 4,096 bytes from the start of its chunk is written when it held a non-zero byte
        for (int i = 0; i < length;)
        {
            if (inside[i]) { i++; continue; }
            int gapEnd = i;
            while (gapEnd < length && !inside[gapEnd]) gapEnd++;
            for (int chunk = i; chunk < gapEnd; chunk += ZeroFill.ChunkBytes)
            {
                int chunkEnd = Math.Min(gapEnd, chunk + ZeroFill.ChunkBytes);
                for (int unit = chunk; unit < chunkEnd; unit += ZeroFill.UnitBytes)
                {
                    int unitEnd = Math.Min(chunkEnd, unit + ZeroFill.UnitBytes);
                    if (before.AsSpan(unit, unitEnd - unit).IndexOfAnyExcept((byte)0) >= 0) expectedCleared += unitEnd - unit;
                }
            }
            i = gapEnd;
        }
        Assert.That(result.ScannedBytes, Is.EqualTo(outsideBytes), "every byte outside the kept ranges was read once");
        Assert.That(result.ClearedBytes, Is.EqualTo(expectedCleared), "a unit is written when it held a non-zero byte, and only then");
        if (nonZeroOutside == 0) Assert.That(result.ClearedBytes, Is.EqualTo(0), "nothing that is zero is written");
        Assert.That(Clear(path, length, keep), Is.EqualTo(new ClearResult(outsideBytes, 0)), "a second pass finds nothing to write");
    }

    [Test]
    public void ClearOutside_DoesNotTakeTheStaleTailOfItsBufferForData([Values(100, 4096, 4097, 9000)] int tail)
    {
        // a gap of one chunk and a short tail: the first chunk is all garbage, the tail is zero, and the buffer still holds garbage behind the tail
        using var dir = new TempDir();
        string path = dir.File("f.bin");
        var bytes = new byte[ZeroFill.ChunkBytes + tail];
        Array.Fill(bytes, (byte)0xA5, 0, ZeroFill.ChunkBytes);
        File.WriteAllBytes(path, bytes);
        Assert.That(Clear(path, bytes.Length, []), Is.EqualTo(new ClearResult(bytes.Length, ZeroFill.ChunkBytes)));
        Assert.That(File.ReadAllBytes(path), Is.All.EqualTo(0));
    }

    [Test]
    public void ClearOutside_WritesTheLastUnitOfAGapOnlyAsFarAsTheGapGoes([Values(1, 100, 4095, 4096, 4097, 5000, 8191)] int length)
    {
        using var dir = new TempDir();
        string path = dir.File("f.bin");
        var bytes = new byte[length + 30];
        bytes[length - 1] = 0xA5; // the last byte of the gap
        bytes[length] = 0x77;     // the first byte of the kept range behind it
        bytes[length + 29] = 0x66;
        File.WriteAllBytes(path, bytes);
        ClearResult r = Clear(path, bytes.Length, [(length, 30)]);
        Assert.That(r, Is.EqualTo(new ClearResult(length, length - ((length - 1) / ZeroFill.UnitBytes * ZeroFill.UnitBytes))), "the unit that holds the byte, as far as the gap reaches");
        byte[] after = File.ReadAllBytes(path);
        Assert.That((after[length - 1], after[length], after[length + 29]), Is.EqualTo(((byte)0, (byte)0x77, (byte)0x66)), "the gap is clear, the bytes behind it are not touched");
    }

    [Test]
    public void ClearOutside_FindsOneByteAtEveryEdge_OfAGapOfALengthThatIsNotAMultipleOfAUnit(
        [Values(0, 1, 4094, 4095, 4096, 4097, ZeroFill.ChunkBytes - 1, ZeroFill.ChunkBytes, ZeroFill.ChunkBytes + 1, ZeroFill.ChunkBytes + 4097)] int at, [Values("start", "middle", "end")] string gap)
    {
        // the gap is [100, 100 + 1 MiB + 9000), with a kept range before and behind it
        const int start = 100, gapLength = ZeroFill.ChunkBytes + 9000;
        int length = start + gapLength + 50;
        using var dir = new TempDir();
        string path = dir.File("f.bin");
        File.WriteAllBytes(path, new byte[length]);
        long position = gap switch { "start" => start + at, "middle" => start + (gapLength / 2) + (at % 7), _ => start + gapLength - 1 - at };
        if (position >= start + gapLength || position < start) Assert.Ignore("the position is outside the gap");
        using (var fs = new FileStream(path, FileMode.Open, FileAccess.Write))
        {
            fs.Position = position;
            fs.WriteByte(0xA5);
            fs.Position = 10;
            fs.WriteByte(0x11); // in the kept range before
            fs.Position = start + gapLength + 7;
            fs.WriteByte(0x22); // in the kept range behind
        }
        ClearResult result = Clear(path, length, [(0, start), (start + gapLength, 50)]);
        byte[] after = File.ReadAllBytes(path);
        Assert.That(after.Count(b => b != 0), Is.EqualTo(2), "only the two bytes inside the kept ranges are left");
        Assert.That((after[10], after[start + gapLength + 7]), Is.EqualTo(((byte)0x11, (byte)0x22)));
        Assert.That(result.ScannedBytes, Is.EqualTo(gapLength));
        Assert.That(result.ClearedBytes, Is.InRange(1, ZeroFill.UnitBytes));
    }

    [Test]
    public void ClearOutside_EmptyKeepClearsTheWholeFile_AFullKeepReadsNothing_AnEmptyFileIsNothing()
    {
        using var dir = new TempDir();
        string path = dir.File("f.bin");
        byte[] before = RandomFile(path, 3 * ZeroFill.UnitBytes, 5);
        Assert.That(Clear(path, before.Length, [(0, before.Length)]), Is.EqualTo(new ClearResult(0, 0)));
        Assert.That(File.ReadAllBytes(path), Is.EqualTo(before));
        Assert.That(Clear(path, before.Length, [(0, 100), (50, 100), (150, before.Length)]).ScannedBytes, Is.EqualTo(0), "overlapping and over-long ranges cover the file");
        Assert.That(Clear(path, before.Length, []).ScannedBytes, Is.EqualTo(before.Length));
        Assert.That(File.ReadAllBytes(path), Is.All.EqualTo(0));
        File.WriteAllBytes(path, []);
        Assert.That(Clear(path, 0, []), Is.EqualTo(new ClearResult(0, 0)));
    }

    [Test]
    public void ClearOutside_RefusesWhatItCannotUse()
    {
        using var dir = new TempDir();
        string path = dir.File("f.bin");
        File.WriteAllBytes(path, new byte[10]);
        using SafeFileHandle handle = File.OpenHandle(path, FileMode.Open, FileAccess.ReadWrite, FileShare.ReadWrite);
        Assert.Throws<ArgumentNullException>(() => ZeroFill.ClearOutside(null!, 10, []));
        Assert.Throws<ArgumentNullException>(() => ZeroFill.ClearOutside(handle, 10, null!));
        Assert.Throws<ArgumentOutOfRangeException>(() => ZeroFill.ClearOutside(handle, -1, []));
        Assert.Throws<ArgumentOutOfRangeException>(() => MappedOutputFile.OpenExisting(path, -1, []));
    }

    [Test]
    public void ClearOutside_AFileShorterThanItsLength_IsAnIoException()
    {
        using var dir = new TempDir();
        string path = dir.File("f.bin");
        File.WriteAllBytes(path, new byte[1000]);
        Assert.Throws<IOException>(() => Clear(path, 5000, [(0, 10)]));
    }

    [Test]
    public void OpenExisting_WithKeep_CutsThenClears_AndReportsIt_WithoutKeepItLeavesEveryByte()
    {
        using var dir = new TempDir();
        string path = dir.File("f.bin");
        byte[] before = RandomFile(path, 20_000, 9);
        using (MappedOutputFile untouched = MappedOutputFile.OpenExisting(path, 15_000))
        {
            Assert.That(untouched.Length, Is.EqualTo(15_000));
            Assert.That(untouched.Cleared, Is.EqualTo(default(ClearResult)), "not asked, not done");
        }
        Assert.That(File.ReadAllBytes(path), Is.EqualTo(before.Take(15_000)), "cut, nothing else changed");
        using (MappedOutputFile cleared = MappedOutputFile.OpenExisting(path, 12_000, [(0, 4096), (8192, 100)]))
        {
            Assert.That(cleared.Length, Is.EqualTo(12_000));
            Assert.That(cleared.Cleared.ScannedBytes, Is.EqualTo(12_000 - 4096 - 100));
            Assert.That(cleared.Outcome, Is.EqualTo(ReservationOutcome.Reopened));
        }
        byte[] after = File.ReadAllBytes(path);
        Assert.That(after.Length, Is.EqualTo(12_000));
        Assert.That(after.Take(4096), Is.EqualTo(before.Take(4096)));
        Assert.That(after.Skip(4096).Take(4096), Is.All.EqualTo(0));
        Assert.That(after.Skip(8192).Take(100), Is.EqualTo(before.Skip(8192).Take(100)));
        Assert.That(after.Skip(8292), Is.All.EqualTo(0));
    }
}
