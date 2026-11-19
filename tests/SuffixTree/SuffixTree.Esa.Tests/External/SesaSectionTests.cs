using System.Runtime.InteropServices;
using SuffixTree.Esa.External;
using SuffixTree.Esa.Tests.TestSupport;

namespace SuffixTree.Esa.Tests.External;

/// <summary>The stages of the builder that were given an overload or an option for the container (work package 1b.7e): <see cref="MappedOutputFile.Grow"/>, the suffix array as a section, PH4 over a range with a placed LCP, and PH4 + PH5 into sections.</summary>
[TestFixture]
public sealed class SesaSectionTests
{
    private static long Align64(long x) => (x + 63) / 64 * 64;

    // ------------------------------------------------------------------ MappedOutputFile.Grow

    [Test]
    public void Grow_KeepsTheContent_MapsTheNewTail_AndUpdatesTheLength()
    {
        using var dir = new TempDir();
        string path = dir.File("g.bin");
        using MappedOutputFile file = MappedOutputFile.Create(path, 100_000);
        using (MappedWindowWriter w = file.CreateWriter(0, 100_000, 4096)) w.Append(Enumerable.Range(0, 100_000).Select(i => (byte)(i * 7)).ToArray());
        file.Complete(false);
        ReservationOutcome outcome = file.Grow(250_000);
        Assert.That(outcome, Is.EqualTo(ReservationOutcome.FreeSpaceChecked).Or.EqualTo(ReservationOutcome.Unchecked));
        Assert.That(file.Length, Is.EqualTo(250_000));
        Assert.That(new FileInfo(path).Length, Is.EqualTo(250_000));
        using (MappedWindowWriter w = file.CreateWriter(100_000, 150_000, 4096)) w.Append(Enumerable.Range(0, 150_000).Select(i => (byte)(i * 11 + 1)).ToArray());
        using (MappedWindowWriter w = file.CreateWriter(99_990, 20, 16)) w.Append(Enumerable.Range(0, 20).Select(i => (byte)200).ToArray());
        file.Complete(true);
        byte[] bytes = File.ReadAllBytes(path);
        for (int i = 0; i < 99_990; i++) if (bytes[i] != (byte)(i * 7)) Assert.Fail($"old byte {i} changed");
        for (int i = 99_990; i < 100_010; i++) if (bytes[i] != 200) Assert.Fail($"straddling write at {i}");
        for (int i = 100_010; i < 250_000; i++) if (bytes[i] != (byte)((i - 100_000) * 11 + 1)) Assert.Fail($"new byte {i}");
    }

    [Test]
    public void Grow_RefusesWhileAWriterIsOpen_ToShrink_AndAfterDisposal_AndKeepsTheFileAsItWas()
    {
        using var dir = new TempDir();
        string path = dir.File("g.bin");
        MappedOutputFile file = MappedOutputFile.Create(path, 1000);
        MappedWindowWriter w = file.CreateWriter(0, 1000, 64);
        Assert.Throws<InvalidOperationException>(() => file.Grow(2000));
        Assert.That(new FileInfo(path).Length, Is.EqualTo(1000));
        w.Dispose();
        Assert.Throws<ArgumentOutOfRangeException>(() => file.Grow(999), "a file does not shrink");
        Assert.That(file.Grow(1000), Is.EqualTo(ReservationOutcome.Sparse), "no growth, no check");
        Assert.That(file.Length, Is.EqualTo(1000));
        Assert.That(file.Grow(3000), Is.Not.EqualTo(ReservationOutcome.Sparse));
        Assert.That(file.Grow(3000), Is.EqualTo(ReservationOutcome.Sparse), "a second call to the same length does nothing");
        Assert.That(file.Grow(4000), Is.Not.EqualTo(ReservationOutcome.Sparse), "and it can grow again");
        file.Dispose();
        Assert.Throws<ObjectDisposedException>(() => file.Grow(5000));
        Assert.That(ProcMaps.Contains(path) ?? false, Is.False);
    }

    [Test]
    public void Grow_ChecksTheFreeSpaceOfTheVolume_AndAFailureLeavesTheFileUsable()
    {
        using var dir = new TempDir();
        string path = dir.File("g.bin");
        using MappedOutputFile file = MappedOutputFile.Create(path, 1000);
        var tight = new FakeReservationSource { CanPreallocate = true, Free = 1999 };
        var e = Assert.Throws<IOException>(() => file.Grow(3000, new MappedFileOptions { Source = tight }));
        Assert.That(e!.Message, Does.Contain("free space"));
        Assert.That((file.Length, new FileInfo(path).Length), Is.EqualTo((1000L, 1000L)));
        Assert.That(tight.LastDirectory, Is.EqualTo(Path.GetDirectoryName(Path.GetFullPath(path))));
        using (MappedWindowWriter w = file.CreateWriter(0, 1000, 64)) w.Append(new byte[1000]);
        var enough = new FakeReservationSource { CanPreallocate = true, Free = 2000 };
        Assert.That(file.Grow(3000, new MappedFileOptions { Source = enough }), Is.EqualTo(ReservationOutcome.FreeSpaceChecked), "exactly the growth is enough");
        var unknown = new FakeReservationSource { CanPreallocate = false, Free = -1 };
        Assert.That(file.Grow(3500, new MappedFileOptions { Source = unknown }), Is.EqualTo(ReservationOutcome.Unchecked));
        Assert.That(unknown.FreeCalls, Is.EqualTo(1));
        Assert.That(unknown.PreallocateCalls, Is.Zero, "the BCL reserves only for a file that it creates");
    }

    [Test]
    public void Grow_ToMoreThan2To32Bytes_IsSparse_AndTheTailIsWritable()
    {
        using var dir = new TempDir();
        string path = dir.File("g.bin");
        using MappedOutputFile file = MappedOutputFile.Create(path, 4096);
        long target = 5_000_000_000;
        file.Grow(target, new MappedFileOptions { Source = new FakeReservationSource { CanPreallocate = false, Free = -1 } });
        Assert.That(file.Length, Is.EqualTo(target));
        using (MappedWindowWriter w = file.CreateWriter(target - 100, 100, 64)) w.Append([.. Enumerable.Range(0, 100).Select(i => (byte)(i + 1))]);
        file.Complete(false);
        var read = new byte[100];
        using (var h = File.OpenHandle(path)) RandomAccess.Read(h, read, target - 100);
        Assert.That(read, Is.EqualTo(Enumerable.Range(0, 100).Select(i => (byte)(i + 1)).ToArray()));
    }

    // ------------------------------------------------------------------ the suffix array as a section

    private static (MappedOutputFile File, string Path) Container(TempDir dir, long length, string name = "c.bin")
    {
        string path = dir.File(name);
        MappedOutputFile file = MappedOutputFile.Create(path, length, new MappedFileOptions { Reservation = ReservationMode.Sparse });
        using (MappedWindowWriter w = file.CreateWriter(512)) w.Append(Enumerable.Repeat((byte)0xA5, (int)length).ToArray());
        return (file, path);
    }

    [Test]
    public void SaBuilder_IntoASection_WritesTheOracleAtTheOffset_AndNothingElse()
    {
        using var dir = new TempDir();
        int v = 0;
        foreach (string text in new[] { "", "a", "ab", "mississippi", "abracadabra", "aaaaaaaaaaaaaaaa", StringFamilies.RandomDna(900, 4), StringFamilies.Fibonacci(700), StringFamilies.RandomProtein(400, 5) })
        {
            v++;
            SaCase c = SaCase.Of(text);
            string textPath = SaRun.WriteText(dir, c, $"-{v}");
            long saOffset = new long[] { 64, 4096, 12 }[v % 3];
            long size = saOffset + (4L * (c.N + 1)) + 77;
            (MappedOutputFile file, string path) = Container(dir, size, $"c{v}.bin");
            using (file)
            {
                using MappedTextReader reader = MappedTextReader.Open(textPath, c.N, c.Bits);
                string scratch = dir.Sub($"scratch-{v}");
                SaBuildResult r = SaBuilder.Build(reader, file, saOffset, SaTiny.Build(scratch, v, 1 + (v % 3)));
                Assert.That(Directory.Exists(scratch), Is.False, "scratch is removed");
                Assert.That(r.Reservation, Is.EqualTo(file.Outcome));
                file.Complete(false); // every writer of the stage is closed; completing the file is the owner's
                byte[] bytes = File.ReadAllBytes(path);
                uint[] sa = MemoryMarshal.Cast<byte, uint>(bytes.AsSpan((int)saOffset, 4 * (c.N + 1))).ToArray();
                Assert.That(sa, Is.EqualTo(c.Sa), $"'{text}'");
                for (int i = 0; i < bytes.Length; i++)
                {
                    if ((i < saOffset || i >= saOffset + (4L * (c.N + 1))) && bytes[i] != 0xA5) Assert.Fail($"'{text}': byte {i} outside the range changed");
                }
            }
        }
    }

    [Test]
    public void SaBuilder_IntoASection_RefusesARangeOutsideTheFile_BeforeAnythingIsDone()
    {
        using var dir = new TempDir();
        SaCase c = SaCase.Of("mississippi");
        string textPath = SaRun.WriteText(dir, c);
        (MappedOutputFile file, _) = Container(dir, 100);
        using (file)
        {
            using MappedTextReader reader = MappedTextReader.Open(textPath, c.N, c.Bits);
            string scratch = dir.Sub("s");
            Assert.Throws<ArgumentOutOfRangeException>(() => SaBuilder.Build(reader, file, 60, SaTiny.Build(scratch, 1)));
            Assert.Throws<ArgumentOutOfRangeException>(() => SaBuilder.Build(reader, file, -1, SaTiny.Build(scratch, 1)));
            Assert.That(Directory.Exists(scratch), Is.False, "no scratch directory was made");
        }
    }

    [Test]
    public void SaBuilder_IntoASection_OnAFailure_LeavesTheFileToItsOwner_AndNothingElse()
    {
        using var dir = new TempDir();
        string text = StringFamilies.RandomDna(2000, 3);
        SaCase c = SaCase.Of(text);
        string textPath = SaRun.WriteText(dir, c);
        int expectedCreations;
        {
            var counting = new FailAtSource(int.MaxValue);
            (MappedOutputFile f, _) = Container(dir, 64 + (4L * 2001), "count.bin");
            using (f)
            {
                using MappedTextReader reader = MappedTextReader.Open(textPath, c.N, c.Bits);
                SaBuilder.Build(reader, f, 64, SaTiny.Build(dir.Sub("s0"), 1) with { Files = new MappedFileOptions { Source = counting } });
                expectedCreations = counting.Calls;
            }
        }
        Assert.That(expectedCreations, Is.GreaterThanOrEqualTo(1), "the rank file is created with a reservation");
        for (int k = 1; k <= expectedCreations; k++)
        {
            (MappedOutputFile file, string path) = Container(dir, 64 + (4L * 2001), $"f{k}.bin");
            using (file)
            {
                using MappedTextReader reader = MappedTextReader.Open(textPath, c.N, c.Bits);
                var counters = new IoCounters();
                string scratch = dir.Sub($"s{k}");
                var e = Assert.Throws<IOException>(() => SaBuilder.Build(reader, file, 64, SaTiny.Build(scratch, 1) with { Counters = counters, Files = new MappedFileOptions { Source = new FailAtSource(k) } }));
                Assert.That(e!.Message, Does.Contain("disk full"));
                Assert.That(File.Exists(path), Is.True, "the file belongs to the caller");
                Assert.That(Directory.Exists(scratch), Is.False);
                Assert.That(counters.ScratchLiveBytes, Is.Zero);
                file.Complete(false);
            }
        }
        (MappedOutputFile capped, _) = Container(dir, 64 + (4L * 2001), "cap.bin");
        using (capped)
        {
            using MappedTextReader reader = MappedTextReader.Open(textPath, c.N, c.Bits);
            Assert.Throws<ArenaBudgetExceededException>(() => SaBuilder.Build(reader, capped, 64, SaTiny.Build(dir.Sub("sc"), 1) with { ArenaCapBytes = 1 }));
            capped.Complete(false);
        }
    }

    [Test]
    public void TheSectionSink_WritesEntriesRelativeToItsRange_AndRefusesWhatIsOutsideIt()
    {
        using var dir = new TempDir();
        (MappedOutputFile file, string path) = Container(dir, 1000);
        using (file)
        {
            var counters = new IoCounters();
            var sink = new MappedSaSectionSink(file, 100, 200, 64, counters);
            using (ISaSliceWriter w = sink.CreateWriter())
            {
                w.Write(0, [1u, 2u, 3u]);
                w.Write(47, [9u, 10u, 11u]);
                Assert.Throws<ArgumentOutOfRangeException>(() => w.Write(49, [1u, 2u]));
            }
            file.Complete(false);
            byte[] b = File.ReadAllBytes(path);
            Assert.That(MemoryMarshal.Cast<byte, uint>(b.AsSpan(100, 12)).ToArray(), Is.EqualTo(new[] { 1u, 2u, 3u }));
            Assert.That(MemoryMarshal.Cast<byte, uint>(b.AsSpan(100 + (47 * 4), 12)).ToArray(), Is.EqualTo(new[] { 9u, 10u, 11u }));
            Assert.That(b[99], Is.EqualTo(0xA5));
            Assert.That(b[300], Is.EqualTo(0xA5));
            Assert.That(counters.MappedWriteBytes, Is.EqualTo(24));
        }
    }

    // ------------------------------------------------------------------ PH4 over a range with a placed LCP, and PH4 + PH5 into sections

    /// <summary>A container in which the suffix array sits at <c>saOffset</c> and the child table is placed after it; the LCP goes behind the child table when the placer is asked, growing the file.</summary>
    private sealed class SectionRun
    {
        internal required LcpChildCase Case { get; init; }

        internal required string Path { get; init; }

        internal required string ScratchDir { get; init; }

        internal required LcpChildResult Result { get; init; }

        internal long SaOffset { get; init; }

        internal long ChildOffset { get; init; }

        internal long LcpOffset { get; set; }

        internal long OverflowOffset { get; set; }

        internal List<LcpPlacementRequest> Requests { get; } = [];

        internal long LengthAtPlacement { get; set; }

        internal byte[]? ChildAtPlacement { get; set; }
    }

    private static SectionRun RunSections(TempDir dir, string text, LcpChildOptions options, string tag = "", int[]? saOverride = null, long saOffset = 128)
    {
        LcpChildCase c = LcpChildCase.Of(text);
        int[] sa = saOverride ?? c.Sa;
        string textPath = TextOracle.WriteFile(dir, $"t{tag}.pk", c.Symbols, c.Bits);
        long childOffset = Align64(saOffset + (4L * (c.N + 1))) + 64;
        long reserved = childOffset + (4L * (c.N + 1)) + 64;
        string path = dir.File($"c{tag}.bin");
        using MappedOutputFile file = MappedOutputFile.Create(path, reserved, new MappedFileOptions { Reservation = ReservationMode.Sparse });
        using (MappedWindowWriter w = file.CreateWriter(saOffset, 4L * (c.N + 1), 256)) w.Append(MemoryMarshal.AsBytes(sa.AsSpan()));
        file.Complete(false);
        var run = new SectionRun { Case = c, Path = path, ScratchDir = options.ScratchDir, Result = null!, SaOffset = saOffset, ChildOffset = childOffset };
        using MappedTextReader reader = MappedTextReader.Open(textPath, c.N, c.Bits);
        Func<LcpPlacementRequest, LcpPlacement> place = request =>
        {
            run.Requests.Add(request);
            run.LengthAtPlacement = file.Length;
            using (var fs = new FileStream(path, FileMode.Open, FileAccess.Read, FileShare.ReadWrite)) { run.ChildAtPlacement = new byte[4 * (c.N + 1)]; fs.Seek(childOffset, SeekOrigin.Begin); fs.ReadExactly(run.ChildAtPlacement); }
            int w = LcpWidth.Bytes(request.Kind);
            long lcp = Align64(reserved);
            long overflow = Align64(lcp + (w * (c.N + 1L)));
            file.Grow(overflow + (8 * request.OverflowCount));
            run.LcpOffset = lcp;
            run.OverflowOffset = overflow;
            return new LcpPlacement(file, lcp, overflow);
        };
        LcpChildResult result = SectionLcpChildBuilder.Build(path, saOffset, reader, file, childOffset, place, options);
        file.Complete(false);
        var final = new SectionRun
        {
            Case = c, Path = path, ScratchDir = options.ScratchDir, Result = result, SaOffset = saOffset, ChildOffset = childOffset,
            LcpOffset = run.LcpOffset, OverflowOffset = run.OverflowOffset, LengthAtPlacement = run.LengthAtPlacement, ChildAtPlacement = run.ChildAtPlacement,
        };
        final.Requests.AddRange(run.Requests);
        return final;
    }

    private static void CompareSections(SectionRun r, string label)
    {
        LcpChildCase c = r.Case;
        byte[] bytes = File.ReadAllBytes(r.Path);
        int n = c.N;
        Assert.That(r.Result.Lcp.Kind, Is.EqualTo(c.Lcp.Kind), $"{label}: width");
        Assert.That(bytes.AsSpan((int)r.LcpOffset, LcpChildRun.Width(c.Lcp.Kind) * (n + 1)).SequenceEqual(c.LcpBaseBytes()), Is.True, $"{label}: LCP");
        Assert.That(bytes.AsSpan((int)r.OverflowOffset, 8 * c.Lcp.OverflowCount).SequenceEqual(c.OverflowBytes()), Is.True, $"{label}: overflow");
        Assert.That(MemoryMarshal.Cast<byte, uint>(bytes.AsSpan((int)r.ChildOffset, 4 * (n + 1))).ToArray(), Is.EqualTo(c.Packed), $"{label}: CHILD");
        Assert.That((r.Result.Scan.NodeCount, r.Result.Scan.InternalCount, (int)r.Result.Scan.LrsDepth, (int)r.Result.Scan.LrsRank, r.Result.Scan.MaxLcp), Is.EqualTo((c.Statistics.NodeCount, (long)c.Statistics.InternalCount, c.Statistics.LrsDepth, c.Statistics.LrsRank, (long)c.MaxLcp)), $"{label}: statistics");
        Assert.That(r.Requests, Has.Count.EqualTo(1), $"{label}: the placer is asked once");
        Assert.That(r.Requests[0], Is.EqualTo(new LcpPlacementRequest(c.Lcp.Kind, c.Lcp.OverflowCount, c.MaxLcp)), $"{label}: what the placer is told");
        Assert.That(Directory.Exists(Path.Combine(r.ScratchDir, "ph4")), Is.False, $"{label}: PH4 scratch");
        Assert.That(Directory.Exists(Path.Combine(r.ScratchDir, "ph5")), Is.False, $"{label}: PH5 scratch");
        Assert.That(Directory.Exists(r.ScratchDir), Is.False, $"{label}: the scratch root that the stage made");
        ProcFds.AssertClosed(label, [r.Path], [r.ScratchDir]);
    }

    [Test]
    public void TheSectionStage_WritesLcpOverflowAndChildAsTheOracleHasThem_ForRotatingTinyParameters()
    {
        using var dir = new TempDir();
        int v = 0;
        foreach (string text in new[] { "", "a", "ab", "aa", "mississippi", "abracadabra", new string('a', 300), StringFamilies.RandomDna(1500, 3), StringFamilies.Fibonacci(800), StringFamilies.ThueMorse(800), LcpChildRun.NRuns(900, 6) + new string('Q', 400) })
        {
            v++;
            LcpChildOptions options = LcpChildRun.Options(dir.Sub($"scr{v}"), LcpChildRun.TinyTuning(v), 1 + (v % 3), sigma: null);
            SectionRun r = RunSections(dir, text, options, $"-{v}", saOffset: new long[] { 128, 4096, 12 }[v % 3]);
            CompareSections(r, $"'{(text.Length > 20 ? text[..20] : text)}' variant {v}");
        }
    }

    [Test]
    public void TheSectionStage_EveryStringOfSmallAlphabets_EqualsTheOracle()
    {
        using var dir = new TempDir();
        int v = 0;
        foreach ((int sigma, int maxLength) in new[] { (1, 6), (2, 7), (3, 4), (4, 3) })
        {
            foreach (byte[] s in TextOracle.AllTexts(sigma, maxLength))
            {
                v++;
                string text = new([.. s.Select(b => (char)('a' + b))]);
                LcpChildOptions options = LcpChildRun.Options(dir.Sub("scr"), LcpChildRun.TinyTuning(v), 1 + (v % 3));
                SectionRun r = RunSections(dir, text, options, "-x");
                CompareSections(r, $"'{text}' variant {v}");
            }
        }
    }

    [Test]
    public void ThePlacer_IsAskedOnce_AfterB2_BeforeAnyEntryOfTheChildTableIsWritten()
    {
        using var dir = new TempDir();
        LcpChildOptions options = LcpChildRun.Options(dir.Sub("scr"), LcpChildRun.TinyTuning(2));
        string text = StringFamilies.RandomDna(2000, 1) + new string('G', 500);
        SectionRun r = RunSections(dir, text, options);
        CompareSections(r, "placer");
        Assert.That(r.ChildAtPlacement!.All(b => b == 0), Is.True, "the child table is still empty when the width is known (it is written in B3)");
        Assert.That(r.LengthAtPlacement, Is.LessThan(new FileInfo(r.Path).Length), "the file grows in the placer");
    }

    [Test]
    public void ASuffixArrayRangeThatIsWrong_IsRefusedWithInvalidData_AndTheScratchIsRemoved_ButTheContainerIsNotTouched()
    {
        using var dir = new TempDir();
        string text = StringFamilies.RandomDna(600, 9);
        LcpChildCase c = LcpChildCase.Of(text);
        (string Label, Action<int[]> Damage)[] damages =
        [
            ("SA[0] is not n", sa => sa[0] = 3),
            ("an entry out of range", sa => sa[100] = c.N + 5),
            ("a repeated entry", sa => sa[200] = sa[201]),
        ];
        foreach ((string label, Action<int[]> damage) in damages)
        {
            int[] sa = (int[])c.Sa.Clone();
            damage(sa);
            string scratch = dir.Sub("scr");
            Assert.Throws<InvalidDataException>(() => RunSections(dir, text, LcpChildRun.Options(scratch, LcpChildRun.TinyTuning(3)), "-bad", sa), label);
            Assert.That(Directory.Exists(scratch), Is.False, $"{label}: scratch");
            Assert.That(File.Exists(dir.File("c-bad.bin")), Is.True, $"{label}: the container is the owner's, not removed by the stage");
            ProcFds.AssertClosed(label, [dir.File("c-bad.bin")], [scratch]);
        }
    }

    [Test]
    public void ASuffixArrayInTheWrongOrderThatIsStillAPermutation_IsNotDetectedByTheStages_AKnownBoundary()
    {
        // The builder trusts the order of the suffix array that PH3 writes: PH4 refuses only what is not a permutation (1b.7c). The verifier's sa_order check (e8-verify) and the Structure level of the reader (1b.7f) catch this.
        using var dir = new TempDir();
        string text = StringFamilies.RandomDna(600, 9);
        LcpChildCase c = LcpChildCase.Of(text);
        int[] sa = (int[])c.Sa.Clone();
        (sa[100], sa[101]) = (sa[101], sa[100]);
        SectionRun wrong = RunSections(dir, text, LcpChildRun.Options(dir.Sub("scr"), LcpChildRun.TinyTuning(3)), "-wrong", sa);
        byte[] bytes = File.ReadAllBytes(wrong.Path);
        Assert.That(MemoryMarshal.Cast<byte, uint>(bytes.AsSpan((int)wrong.ChildOffset, 4 * (c.N + 1))).ToArray(), Is.Not.EqualTo(c.Packed), "the child table of the wrong order is wrong, and nothing said so");
        Assert.That(bytes.AsSpan((int)wrong.LcpOffset, c.N + 1).SequenceEqual(c.LcpBaseBytes()), Is.False, "so is the LCP");
        // the same for PH6: BWT and PSI of a wrong order are built without a complaint
        string textPath = TextOracle.WriteFile(dir, "w.pk", c.Symbols, c.Bits);
        using MappedTextReader reader = MappedTextReader.Open(textPath, c.N, c.Bits);
        BwtResult b = BwtPass.Build(wrong.Path, wrong.SaOffset, reader, c.Sigma, dir.File("w.bwt"));
        Assert.That(b.PrimaryRow, Is.EqualTo(Array.IndexOf(sa, 0)));
    }

    [Test]
    public void APlacerThatThrows_OrReturnsARangeOutsideTheFile_FailsTheStage_AndKeepsTheType()
    {
        using var dir = new TempDir();
        string text = StringFamilies.RandomDna(400, 2);
        LcpChildCase c = LcpChildCase.Of(text);
        string textPath = TextOracle.WriteFile(dir, "t.pk", c.Symbols, c.Bits);
        long childOffset = 4096 + 64;
        long reserved = childOffset + (4L * (c.N + 1)) + 64;
        foreach (bool outside in new[] { false, true })
        {
            string path = dir.File($"c{outside}.bin");
            using MappedOutputFile file = MappedOutputFile.Create(path, reserved, new MappedFileOptions { Reservation = ReservationMode.Sparse });
            using (MappedWindowWriter w = file.CreateWriter(128, 4L * (c.N + 1), 256)) w.Append(MemoryMarshal.AsBytes(c.Sa.AsSpan()));
            using MappedTextReader reader = MappedTextReader.Open(textPath, c.N, c.Bits);
            string scratch = dir.Sub($"s{outside}");
            Func<LcpPlacementRequest, LcpPlacement> place = outside ? _ => new LcpPlacement(file, reserved, reserved + 10_000) : _ => throw new InjectedFailureException("placer");
            Type expected = outside ? typeof(ArgumentOutOfRangeException) : typeof(InjectedFailureException);
            Exception? e = null;
            try
            {
                SectionLcpChildBuilder.Build(path, 128, reader, file, childOffset, place, LcpChildRun.Options(scratch, LcpChildRun.TinyTuning(1)));
            }
            catch (Exception ex)
            {
                e = ex;
            }
            Assert.That(e, Is.InstanceOf(expected), $"outside={outside}: {e}");
            Assert.That(Directory.Exists(scratch), Is.False);
            file.Complete(false); // every writer the stage opened was closed
        }
    }

    [Test]
    public void TheSectionStage_RefusesAChildRangeThatIsNotThere_AndANullPlacer()
    {
        using var dir = new TempDir();
        string textPath = TextOracle.WriteFile(dir, "t.pk", [0, 1, 0], 1);
        (MappedOutputFile file, string path) = Container(dir, 200);
        using (file)
        {
            using MappedTextReader reader = MappedTextReader.Open(textPath, 3, 1);
            LcpChildOptions o = LcpChildRun.Options(dir.Sub("s"), LcpChildRun.TinyTuning(0));
            Assert.Throws<ArgumentNullException>(() => SectionLcpChildBuilder.Build(path, 64, reader, file, 128, null!, o));
            Assert.Throws<ArgumentOutOfRangeException>(() => SectionLcpChildBuilder.Build(path, 64, reader, file, 190, _ => default, o), "the child range ends outside the file");
            Assert.Throws<ArgumentOutOfRangeException>(() => SectionLcpChildBuilder.Build(path, -1, reader, file, 128, _ => default, o));
        }
    }

}
