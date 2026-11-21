using System.Collections.Concurrent;
using System.Security.Cryptography;
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
}
