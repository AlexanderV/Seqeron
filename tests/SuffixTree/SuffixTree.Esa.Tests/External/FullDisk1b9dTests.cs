using System.Globalization;
using SuffixTree.Esa.External;
using SuffixTree.Esa.Tests.TestSupport;

namespace SuffixTree.Esa.Tests.External;

/// <summary>
/// A full disk at any point of the container build (work package 1b.9d): the late pre-flight and the early one are advisory (another process can fill the volume after them), so every page that a mapped writer touches
/// must have its block already, whichever phase writes it. The growth of the container behind its reserved length (the LCP above one byte per entry, and the overflow list) was sparse until this package: a page of it that
/// the volume could not give was a signal that ended the process. The tests here fail the reservation of that range by injection, fill a real tmpfs at every boundary of the build and give the build a real volume of
/// every size between a third of what it needs and more than it needs, with a pre-flight that cannot tell the free space: the build is a success or an <see cref="IOException"/>, nothing is left, and the process lives.
/// </summary>
[TestFixture]
public sealed class FullDisk1b9dTests
{
    /// <summary>An LCP of two bytes per entry and no overflow list: 32,600 symbols, a run of 2,600 equal letters.</summary>
    private static string TextU16 => StringFamilies.RandomDna(30_000, 5) + new string('A', 2_600);

    /// <summary>An LCP of four bytes per entry (the run is longer than 65,535).</summary>
    private static string TextU32 => new('a', 100_000);

    /// <summary>An LCP of one byte per entry with an overflow list (445 pairs, 8 bytes each).</summary>
    private static string TextOverflow => StringFamilies.RandomDna(50_000, 6) + new string('A', 700);

    /// <summary>An LCP of one byte per entry and no overflow list: the file does not grow.</summary>
    private static string TextU8 => StringFamilies.RandomDna(20_000, 7);

    private static SesaBuildOptions Options(int threads = 2, bool psi = true) => new() { Threads = threads, IncludePsi = psi, Fsync = false, MeasureRss = false };

    private static void AssertNothingLeft(string label, string input, string target, string outDir)
    {
        string tmp = target + ".tmp";
        Assert.That(File.Exists(target), Is.False, $"{label}: a target was published");
        Assert.That(File.Exists(tmp), Is.False, $"{label}: the temporary file is left");
        Assert.That(Directory.Exists(tmp + ".d"), Is.False, $"{label}: the scratch directory is left");
        Assert.That(Directory.GetFileSystemEntries(outDir), Is.Empty, $"{label}: something is left in the output directory");
        ProcFds.AssertClosed(label, [input, target, tmp], [tmp + ".d"]);
        Assert.That(ProcMaps.Contains(tmp) ?? false, Is.False, $"{label}: the temporary file is still mapped");
    }

    private static Exception? Run(string input, string target, SesaBuildOptions options)
    {
        try
        {
            SesaBuilder.Build(input, target, options);
            return null;
        }
        catch (Exception e)
        {
            return e;
        }
    }

    private static (string Input, string Target, string OutDir) Setup(TempDir dir, string text, string root = "out")
    {
        string input = SesaRun.WriteInput(dir, text);
        string outDir = dir.Sub(root);
        Directory.CreateDirectory(outDir);
        return (input, Path.Combine(outDir, "x.sesa"), outDir);
    }

    // ------------------------------------------------------------------ the kinds of tail

    [Test]
    public void TheTextsOfTheTests_HaveTheKindsOfTailThatTheyAreMeantFor()
    {
        using var dir = new TempDir();
        foreach ((string name, string text, LcpKind kind, bool grows) in new[]
        {
            ("u8", TextU8, LcpKind.U8, false),
            ("overflow", TextOverflow, LcpKind.U8, true),
            ("u16", TextU16, LcpKind.U16, true),
            ("u32", TextU32, LcpKind.U32, true),
        })
        {
            (string input, string target, string outDir) = Setup(dir, text, "k" + name);
            SesaBuildResult r = SesaBuilder.Build(input, target, Options());
            Assert.That(r.LcpKind, Is.EqualTo(kind), name);
            Assert.That(r.GrownBytes > 0, Is.EqualTo(grows), name);
            Assert.That(r.OverflowCount > 0, Is.EqualTo(name == "overflow"), name);
            Assert.That(outDir, Is.Not.Null);
        }
    }

    // ------------------------------------------------------------------ the growth is reserved before the LCP is written

    [Test]
    public void TheGrowthOfTheContainer_IsReserved_BeforeAnyByteOfTheLcpIsWritten_AndTheResultSaysSo()
    {
        foreach ((string name, string text) in new[] { ("overflow", TextOverflow), ("u16", TextU16), ("u32", TextU32) })
        {
            using var dir = new TempDir();
            (string input, string target, _) = Setup(dir, text);
            var source = new TailSource { Can = true, Free = -1 };
            var counters = new IoCounters();
            List<(long Offset, long Length)>? atPlacement = null;
            SesaBuildOptions options = Options() with
            {
                Files = new MappedFileOptions { Source = source },
                Counters = counters,
                Hook = point =>
                {
                    if (point == SesaHooks.LcpPlaced) atPlacement = [.. source.Reserved];
                },
            };
            SesaBuildResult r = SesaBuilder.Build(input, target, options);
            PackedTextInfo info = PackedTextWriter.Scan(input);
            SesaPlan plan = SesaLayout.Plan(info, true);
            Assert.That(r.GrownBytes, Is.GreaterThan(0), name);
            Assert.That(atPlacement, Is.EqualTo(new[] { (plan.ReservedLength, r.FileSize - plan.ReservedLength) }), $"{name}: the range that was added is reserved when the LCP is placed, before PH4 writes it");
            Assert.That(r.GrownBytes, Is.EqualTo(r.FileSize - plan.ReservedLength), name);
            Assert.That(r.GrowthReservation, Is.EqualTo(ReservationOutcome.Filled), name);
            Assert.That(counters.ReservedWriteBytes, Is.EqualTo(r.GrownBytes), $"{name}: every reserved byte is counted");
            Assert.That(r.Io.WriteBytes, Is.Not.Zero);
            SesaOracle.Read(target).AssertConforms(name, true, SesaFormat.ChildCodecTagged);
        }
    }

    [Test]
    public void ABuildWhoseLcpDoesNotGrow_ReservesNothingMore_AndCountsNothing()
    {
        using var dir = new TempDir();
        (string input, string target, _) = Setup(dir, TextU8);
        var source = new TailSource { Can = true, Free = -1 };
        var counters = new IoCounters();
        SesaBuildResult r = SesaBuilder.Build(input, target, Options() with { Files = new MappedFileOptions { Source = source }, Counters = counters });
        Assert.That(source.Reserved, Is.Empty);
        Assert.That(r.GrowthReservation, Is.Null);
        Assert.That(r.GrownBytes, Is.Zero);
        Assert.That(counters.ReservedWriteBytes, Is.Zero);
    }

    [Test]
    public void TheSparseMode_StillGrowsTheContainerSparse_AndReservesNothing()
    {
        using var dir = new TempDir();
        (string input, string target, _) = Setup(dir, TextU16);
        var source = new TailSource { Can = true, Free = 1L << 40 };
        var counters = new IoCounters();
        SesaBuildResult r = SesaBuilder.Build(input, target, Options() with { Files = new MappedFileOptions { Source = source, Reservation = ReservationMode.Sparse }, Counters = counters });
        Assert.That(source.Reserved, Is.Empty);
        Assert.That(r.Reservation, Is.EqualTo(ReservationOutcome.Sparse));
        Assert.That(r.GrowthReservation, Is.EqualTo(ReservationOutcome.FreeSpaceChecked));
        Assert.That(counters.ReservedWriteBytes, Is.Zero);
    }

    // ------------------------------------------------------------------ a disk that fills inside the reservation

    [Test]
    public void ADiskThatFillsInTheReservationOfTheTail_IsAnIoException_AndLeavesNothing()
    {
        foreach ((string name, string text) in new[] { ("overflow", TextOverflow), ("u16", TextU16), ("u32", TextU32) })
        {
            using var dir = new TempDir();
            (string input, string target, string outDir) = Setup(dir, text);
            SesaBuildResult baseline = SesaBuilder.Build(input, target, Options());
            File.Delete(target);
            long grown = baseline.GrownBytes;
            foreach (long after in new[] { 0L, 1L, grown / 2, grown - 1 })
            {
                foreach ((int threads, bool psi) in new[] { (1, true), (3, false) })
                {
                    var source = new TailSource { Can = true, Free = -1, FailAfterBytes = after };
                    var counters = new IoCounters();
                    Exception? e = Run(input, target, Options(threads, psi) with { Files = new MappedFileOptions { Source = source }, Counters = counters });
                    Assert.That(e, Is.InstanceOf<IOException>(), $"{name} after {after} bytes, threads {threads}, psi {psi}: {e}");
                    Assert.That(e!.Message, Does.Contain("disk full"));
                    Assert.That(source.Reserved, Has.Count.EqualTo(1), "the reservation is asked for once");
                    Assert.That(counters.ReservedWriteBytes, Is.EqualTo(after));
                    Assert.That(counters.ScratchLiveBytes, Is.Zero, "the scratch is off the books");
                    AssertNothingLeft($"{name} after {after}", input, target, outDir);
                }
            }
        }
    }

    // ------------------------------------------------------------------ a real volume

    /// <summary>The bytes of a volume that holds the build of the text: the container, the scratch at its peak and a little room.</summary>
    private static long VolumeBytesFor(TempDir dir, string text, bool psi, out SesaBuildResult baseline)
    {
        (string input, string target, _) = Setup(dir, text, "base");
        baseline = SesaBuilder.Build(input, target, Options(2, psi));
        long bytes = baseline.FileSize + baseline.ScratchPeakBytes + (256 << 10);
        File.Delete(target);
        return bytes;
    }

    [Test]
    public void TheVolumeFillsRightAfterTheGrowth_TheLcpIsStillWritten_AndTheBuildEndsWithoutASignal()
    {
        foreach ((string name, string text) in new[] { ("overflow", TextOverflow), ("u16", TextU16), ("u32", TextU32) })
        {
            using var dir = new TempDir();
            long bytes = VolumeBytesFor(dir, text, false, out SesaBuildResult baseline);
            Assert.That(baseline.GrownBytes, Is.GreaterThan(0));
            using SmallVolume volume = SmallVolume.MountOrIgnore(dir, "mnt", bytes);
            string outDir = volume.File("out");
            Directory.CreateDirectory(outDir);
            string input = SesaRun.WriteInput(dir, text, "in2.txt");
            string target = Path.Combine(outDir, "x.sesa");
            long freeAtFill = -1;
            SesaBuildOptions options = Options(2, false) with
            {
                Files = new MappedFileOptions { Source = new BlindSource() },
                Hook = point =>
                {
                    if (point != SesaHooks.LcpPlaced) return;
                    volume.FillToTheBrim(); // after the late pre-flight, after the growth: the disk is gone
                    freeAtFill = volume.Free;
                },
            };
            SesaBuildResult r = SesaBuilder.Build(input, target, options);
            Assert.That(freeAtFill, Is.Zero, $"{name}: the volume really was full while PH4 wrote the LCP");
            Assert.That(r.GrownBytes, Is.EqualTo(baseline.GrownBytes));
            volume.RemoveBallast();
            SesaOracle.Read(target).AssertConforms(name, false, SesaFormat.ChildCodecTagged);
            Assert.That(File.ReadAllBytes(target), Is.EqualTo(ReadBaseline(dir, text, false)), $"{name}: the same bytes as a build on a volume with room");
        }
    }

    private static byte[] ReadBaseline(TempDir dir, string text, bool psi)
    {
        (string input, string target, _) = Setup(dir, text, "ref" + (psi ? "p" : "n"));
        SesaBuilder.Build(input, target, Options(2, psi));
        return File.ReadAllBytes(target);
    }

    [Test]
    public void AFullDisk_AtEveryBoundaryOfTheBuild_IsASuccessOrAnIoException_NeverASignal_AndNothingIsLeft()
    {
        foreach ((string name, string text) in new[] { ("u16", TextU16), ("overflow", TextOverflow) })
        {
            foreach (bool psi in new[] { true, false })
            {
                using var dir = new TempDir();
                long bytes = VolumeBytesFor(dir, text, psi, out SesaBuildResult baseline);
                string[] points = psi ? SesaHooks.All : [.. SesaHooks.All.Where(h => h is not (SesaHooks.BwtBuilt or SesaHooks.PsiBuilt))];
                string input = SesaRun.WriteInput(dir, text, "in2.txt");
                foreach (string point in points)
                {
                    using SmallVolume volume = SmallVolume.MountOrIgnore(dir, "mnt-" + point, bytes);
                    string outDir = volume.File("out");
                    Directory.CreateDirectory(outDir);
                    string target = Path.Combine(outDir, "x.sesa");
                    SesaBuildOptions options = Options(2, psi) with
                    {
                        Files = new MappedFileOptions { Source = new BlindSource() },
                        Hook = p =>
                        {
                            if (p == point) volume.FillToTheBrim();
                        },
                    };
                    Exception? e = Run(input, target, options);
                    string label = $"{name} psi {psi} full at {point}";
                    volume.RemoveBallast();
                    if (e is null)
                    {
                        SesaOracle.Read(target).AssertConforms(label, psi, SesaFormat.ChildCodecTagged);
                        Assert.That(new FileInfo(target).Length, Is.EqualTo(baseline.FileSize), label);
                        continue;
                    }
                    Assert.That(e, Is.InstanceOf<IOException>(), $"{label}: {e}");
                    AssertNothingLeft(label, input, target, outDir);
                }
            }
        }
    }

    [Test]
    public void AFullDisk_AtEveryPointThatThePhasesReport_IsASuccessOrAnIoException_NeverASignal_AndNothingIsLeft()
    {
        foreach ((string name, string text) in new[] { ("u16", TextU16), ("u8", TextU8) })
        {
            using var dir = new TempDir();
            var seen = new SortedSet<string>(StringComparer.Ordinal);
            string input = SesaRun.WriteInput(dir, text, "in2.txt");
            (_, string baseTarget, _) = Setup(dir, text, "base");
            SesaBuildResult baseline = SesaBuilder.Build(input, baseTarget, Options(3) with
            {
                SaTrace = point => { lock (seen) seen.Add("sa:" + point); },
                LcpTrace = point => { lock (seen) seen.Add("lcp:" + point); },
            });
            File.Delete(baseTarget);
            long bytes = baseline.FileSize + baseline.ScratchPeakBytes + (256 << 10);
            Assert.That(seen, Has.Count.GreaterThanOrEqualTo(8), "the phases report their points");
            foreach (string point in seen)
            {
                using SmallVolume volume = SmallVolume.MountOrIgnore(dir, "mnt-" + point.Replace(':', '-'), bytes);
                string outDir = volume.File("out");
                Directory.CreateDirectory(outDir);
                string target = Path.Combine(outDir, "x.sesa");
                SesaBuildOptions options = Options(3) with
                {
                    Files = new MappedFileOptions { Source = new BlindSource() },
                    SaTrace = p => { if ("sa:" + p == point) volume.FillToTheBrim(); },
                    LcpTrace = p => { if ("lcp:" + p == point) volume.FillToTheBrim(); },
                };
                Exception? e = Run(input, target, options);
                string label = $"{name} full at {point}";
                volume.RemoveBallast();
                if (e is null)
                {
                    SesaOracle.Read(target).AssertConforms(label, true, SesaFormat.ChildCodecTagged);
                    continue;
                }
                Assert.That(e, Is.InstanceOf<IOException>(), $"{label}: {e}");
                AssertNothingLeft(label, input, target, outDir);
            }
        }
    }

    [Test]
    public void AVolumeOfAnySize_GivesASuccessOrAnIoException_NeverASignal_WhateverThePreflightKnows()
    {
        int builds = 0, failures = 0;
        foreach ((string name, string text) in new[] { ("u16", TextU16), ("u32", TextU32), ("overflow", TextOverflow) })
        {
            using var dir = new TempDir();
            long full = VolumeBytesFor(dir, text, true, out SesaBuildResult baseline);
            string input = SesaRun.WriteInput(dir, text, "in2.txt");
            long lowest = baseline.FileSize / 3;
            long step = Math.Max(4096, (full - lowest) / 20 / 4096 * 4096);
            for (long bytes = lowest; bytes <= full + step; bytes += step)
            {
                using SmallVolume volume = SmallVolume.MountOrIgnore(dir, $"mnt-{bytes}", bytes);
                string outDir = volume.File("out");
                Directory.CreateDirectory(outDir);
                string target = Path.Combine(outDir, "x.sesa");
                Exception? e = Run(input, target, Options() with { Files = new MappedFileOptions { Source = new BlindSource() } });
                string label = $"{name} on {bytes.ToString(CultureInfo.InvariantCulture)} bytes";
                builds++;
                if (e is null)
                {
                    SesaOracle.Read(target).AssertConforms(label, true, SesaFormat.ChildCodecTagged);
                    continue;
                }
                failures++;
                Assert.That(e, Is.InstanceOf<IOException>(), $"{label}: {e}");
                AssertNothingLeft(label, input, target, outDir);
            }
        }
        Assert.That(builds, Is.GreaterThan(50));
        Assert.That(failures, Is.GreaterThan(builds / 4), "the sizes reach below what the build needs, so the failure paths ran");
        Assert.That(builds - failures, Is.GreaterThan(5), "and above it, so the success path ran");
    }
}
