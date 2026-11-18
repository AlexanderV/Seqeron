using System.Security.Cryptography;
using SuffixTree.Esa.External;
using SuffixTree.Esa.Tests.TestSupport;

namespace SuffixTree.Esa.Tests.External;

/// <summary>The streaming writer of the packed text (PH0 of ADR-004 section 4.2 as far as the codec needs it): round trip, determinism, atomic finish, failure modes, sections inside a larger file.</summary>
[TestFixture]
public class PackedTextWriterTests
{
    private static byte[] RandomInput(int n, int sigma, ulong seed, int offset = 'a')
    {
        var rng = new SplitMix64(seed);
        var b = new byte[n];
        for (int i = 0; i < n; i++) b[i] = (byte)(offset + rng.NextInt(sigma));
        return b;
    }

    private static string Input(TempDir dir, byte[] data, string name = "in.txt")
    {
        string p = dir.File(name);
        File.WriteAllBytes(p, data);
        return p;
    }

    private static void CheckRoundTrip(TempDir dir, byte[] input, string label)
    {
        string txt = Input(dir, input);
        string pk = dir.File("out.pk");
        PackedTextInfo info = PackedTextWriter.WriteFile(txt, pk);
        (byte[] dense, byte[] letters) = TextOracle.Dense(input);
        int bits = TextOracle.BitsOf(letters.Length);
        Assert.That(info.N, Is.EqualTo(input.Length), label);
        Assert.That(info.Bits, Is.EqualTo(bits), label);
        Assert.That(info.Letters, Is.EqualTo(letters), $"{label}: increasing byte order");
        Assert.That(info.Counts, Is.EqualTo(letters.Select(l => (long)input.Count(b => b == l)).ToArray()), label);
        Assert.That(info.Sha256, Is.EqualTo(Convert.ToHexStringLower(SHA256.HashData(input))), label);
        Assert.That(File.ReadAllBytes(pk), Is.EqualTo(TextOracle.Pack(dense, bits)), $"{label}: the file");
        using MappedTextReader r = MappedTextReader.Open(pk, info.N, info.Bits);
        r.Verify(info);
        var back = new byte[input.Length];
        r.ReadSymbols(0, back);
        Assert.That(back.Select(s => letters[s]).ToArray(), Is.EqualTo(input), $"{label}: reads back as the input");
    }

    [Test]
    public void RoundTrip_ForEverySigma_AndSeveralLengths()
    {
        using var dir = new TempDir();
        foreach (int sigma in new[] { 1, 2, 3, 4, 5, 7, 8, 9, 16, 17, 100, 255, 256 })
        {
            foreach (int n in new[] { 1, 7, 8, 9, 1000, 100_003 })
            {
                byte[] input = RandomInput(n, sigma, (ulong)((sigma * 1000) + n), sigma == 256 ? 0 : 3);
                for (int i = 0; i < Math.Min(sigma, n); i++) input[i] = (byte)(i + (sigma == 256 ? 0 : 3)); // every letter occurs
                CheckRoundTrip(dir, input, $"sigma={sigma} n={n}");
            }
        }
    }

    [Test]
    public void TheWidth_IsTheSmallestBThatHoldsSigma_AndAtLeastOne()
    {
        using var dir = new TempDir();
        foreach ((int sigma, int bits) in new[] { (1, 1), (2, 1), (3, 2), (4, 2), (5, 3), (8, 3), (9, 4), (16, 4), (17, 5), (129, 8), (256, 8) })
        {
            var data = new byte[3 * sigma];
            for (int i = 0; i < data.Length; i++) data[i] = (byte)(((i * 7) % sigma) + (sigma == 256 ? 0 : 3));
            for (int i = 0; i < sigma; i++) data[i] = (byte)(i + (sigma == 256 ? 0 : 3));
            PackedTextInfo info = PackedTextWriter.WriteFile(Input(dir, data), dir.File("w.pk"));
            Assert.That((info.Sigma, info.Bits), Is.EqualTo((sigma, bits)));
        }
    }

    [Test]
    public void NContainingDna_And_AllByteValues_AreOrdinaryAlphabets()
    {
        using var dir = new TempDir();
        var rng = new SplitMix64(3);
        byte[] dna = new byte[50_000];
        byte[] alphabet = "ACGTN"u8.ToArray();
        for (int i = 0; i < dna.Length; i++) dna[i] = alphabet[rng.NextInt(i % 1000 < 100 ? 5 : 4)];
        Array.Fill(dna, (byte)'N', 20_000, 5_000);
        CheckRoundTrip(dir, dna, "ACGTN");
        PackedTextInfo info = PackedTextWriter.WriteFile(Input(dir, dna), dir.File("n.pk"));
        Assert.That((info.Sigma, info.Bits), Is.EqualTo((5, 3)), "three bits: 5 letters; ADR-004 6a Q3 writes 4 bits for A, C, G, T, N and IUPAC, the layout is b = ceil(log2 sigma) as in the prototype");
        var all = new byte[2000];
        for (int i = 0; i < all.Length; i++) all[i] = (byte)(i % 256);
        CheckRoundTrip(dir, all, "all 256 byte values");
        Assert.That(File.ReadAllBytes(dir.File("out.pk")).Take(2000), Is.EqualTo(all), "8 bits per symbol: the file is the text");
    }

    [Test]
    public void MinBits_RaisesTheWidth_WithoutChangingTheLayoutRules()
    {
        using var dir = new TempDir();
        byte[] input = RandomInput(1003, 4, 9, 'A');
        PackedTextInfo info = PackedTextWriter.WriteFile(Input(dir, input), dir.File("m.pk"), new PackedTextWriterOptions { MinBits = 4 });
        Assert.That(info.Bits, Is.EqualTo(4));
        (byte[] dense, _) = TextOracle.Dense(input);
        Assert.That(File.ReadAllBytes(dir.File("m.pk")), Is.EqualTo(TextOracle.Pack(dense, 4)));
        using MappedTextReader r = MappedTextReader.Open(dir.File("m.pk"), info.N, info.Bits);
        r.Verify(info);
        PackedTextInfo wide = PackedTextWriter.WriteFile(Input(dir, input), dir.File("m2.pk"), new PackedTextWriterOptions { MinBits = 2 });
        Assert.That(wide.Bits, Is.EqualTo(2), "a minimum below the needed width changes nothing");
        Assert.Throws<ArgumentOutOfRangeException>(() => PackedTextWriter.WriteFile(Input(dir, input), dir.File("m3.pk"), new PackedTextWriterOptions { MinBits = 0 }));
        Assert.Throws<ArgumentOutOfRangeException>(() => PackedTextWriter.WriteFile(Input(dir, input), dir.File("m4.pk"), new PackedTextWriterOptions { MinBits = 9 }));
    }

    [Test]
    public void TheFile_DoesNotDependOnTheBlockSizeOrTheWindow_AndTwoWritesAreByteIdentical()
    {
        using var dir = new TempDir();
        byte[] input = RandomInput(10_007, 3, 11);
        string txt = Input(dir, input);
        byte[]? first = null;
        PackedTextInfo? firstInfo = null;
        foreach (int block in new[] { 8, 16, 24, 40, 1000, 1 << 20 })
        {
            foreach (int window in new[] { 1, 7, 4096, 1 << 20 })
            {
                string pk = dir.File($"p{block}_{window}.pk");
                PackedTextInfo info = PackedTextWriter.WriteFile(txt, pk, new PackedTextWriterOptions { BlockSymbols = block, WindowBytes = window });
                byte[] bytes = File.ReadAllBytes(pk);
                first ??= bytes;
                firstInfo ??= info;
                Assert.That(bytes, Is.EqualTo(first), $"block {block} window {window}");
                Assert.That(info.Sha256, Is.EqualTo(firstInfo.Sha256));
                Assert.That(info.Counts, Is.EqualTo(firstInfo.Counts));
            }
        }
        PackedTextWriter.WriteFile(txt, dir.File("again.pk"));
        Assert.That(File.ReadAllBytes(dir.File("again.pk")), Is.EqualTo(first), "a second write of the same input");
    }

    [TestCase(0)]
    [TestCase(7)]
    [TestCase(9)]
    [TestCase(100)]
    [TestCase(-8)]
    public void ABlockThatIsNotAPositiveMultipleOfEight_IsRefused(int block)
    {
        using var dir = new TempDir();
        string txt = Input(dir, "abc"u8.ToArray());
        Assert.Throws<ArgumentOutOfRangeException>(() => PackedTextWriter.WriteFile(txt, dir.File("p.pk"), new PackedTextWriterOptions { BlockSymbols = block }));
        Assert.Throws<ArgumentOutOfRangeException>(() => PackedTextWriter.Scan(txt, new PackedTextWriterOptions { BlockSymbols = block }));
        Assert.That(Directory.GetFiles(dir.Path).Select(Path.GetFileName), Is.EquivalentTo(new[] { "in.txt" }), "nothing was created");
    }

    [Test]
    public void SingleSymbolTexts_AndTheEmptyText()
    {
        using var dir = new TempDir();
        PackedTextInfo one = PackedTextWriter.WriteFile(Input(dir, "aaaaaaaaaaaa"u8.ToArray()), dir.File("a.pk"));
        Assert.That((one.Sigma, one.Bits, one.N), Is.EqualTo((1, 1, 12L)));
        Assert.That(new FileInfo(dir.File("a.pk")).Length, Is.EqualTo(2 + 64));
        Assert.That(File.ReadAllBytes(dir.File("a.pk")), Is.All.Zero, "the only symbol is 0");
        PackedTextInfo z = PackedTextWriter.WriteFile(Input(dir, "z"u8.ToArray()), dir.File("z.pk"));
        Assert.That((z.Sigma, z.Bits, z.N), Is.EqualTo((1, 1, 1L)));
        Assert.That(new FileInfo(dir.File("z.pk")).Length, Is.EqualTo(1 + 64));

        PackedTextInfo empty = PackedTextWriter.WriteFile(Input(dir, [], "empty.txt"), dir.File("e.pk"));
        Assert.That((empty.N, empty.Sigma, empty.Bits), Is.EqualTo((0L, 0, 1)));
        Assert.That(empty.Letters, Is.Empty);
        Assert.That(empty.Counts, Is.Empty);
        Assert.That(empty.Sha256, Is.EqualTo(Convert.ToHexStringLower(SHA256.HashData([]))));
        Assert.That(File.ReadAllBytes(dir.File("e.pk")), Is.EqualTo(new byte[64]));
        using MappedTextReader r = MappedTextReader.Open(dir.File("e.pk"), 0, 1);
        r.Verify(empty);
    }

    [Test]
    public void ScanAlone_FindsWhatTheWriteReports_AndCreatesNothing()
    {
        using var dir = new TempDir();
        byte[] input = RandomInput(5000, 5, 7, 'A');
        string txt = Input(dir, input);
        PackedTextInfo scan = PackedTextWriter.Scan(txt);
        PackedTextInfo written = PackedTextWriter.WriteFile(txt, dir.File("s.pk"));
        Assert.That((scan.N, scan.Bits, scan.Sha256), Is.EqualTo((written.N, written.Bits, written.Sha256)));
        Assert.That(scan.Letters, Is.EqualTo(written.Letters));
        Assert.That(scan.Counts, Is.EqualTo(written.Counts));
        Assert.That(Directory.GetFiles(dir.Path).Select(Path.GetFileName), Is.EquivalentTo(new[] { "in.txt", "s.pk" }), "the scan created nothing");
    }

    [Test]
    public void TheInput_IsReadTwice_TheOutputIsWrittenThroughMappedWindows_NothingIsWrittenExplicitly()
    {
        using var dir = new TempDir();
        byte[] input = RandomInput(100_000, 4, 5);
        var c = new IoCounters();
        PackedTextWriter.WriteFile(Input(dir, input), dir.File("p.pk"), new PackedTextWriterOptions { BlockSymbols = 4096, WindowBytes = 1000, Counters = c });
        Assert.That(c.ReadBytes, Is.EqualTo(200_000), "two passes");
        Assert.That(c.WriteBytes, Is.Zero);
        Assert.That(c.MappedWriteBytes, Is.EqualTo(25_000 + 64), "2 bits per symbol and the padding");
        Assert.That(c.MappedMaxWindowBytes, Is.LessThanOrEqualTo(1000));
        Assert.That(c.MappedWindows, Is.EqualTo((25_064 + 999) / 1000));
    }

    [Test]
    public void ASymbolicLinkToTheText_IsTheText_NotTheLink()
    {
        using var dir = new TempDir();
        byte[] input = RandomInput(7001, 4, 8);
        string target = Input(dir, input);
        string link = dir.File("link.txt");
        File.CreateSymbolicLink(link, target);
        Assert.That(new FileInfo(link).Length, Is.Not.EqualTo(input.Length), "FileInfo.Length of a link is the length of its target path: the reason for the handle");
        Assert.That(PackedTextWriter.InputLength(link), Is.EqualTo(7001));
        PackedTextInfo viaLink = PackedTextWriter.WriteFile(link, dir.File("l.pk"));
        PackedTextWriter.WriteFile(target, dir.File("d.pk"));
        Assert.That(viaLink.N, Is.EqualTo(7001));
        Assert.That(File.ReadAllBytes(dir.File("l.pk")), Is.EqualTo(File.ReadAllBytes(dir.File("d.pk"))));
    }

    // ---- limits

    [Test]
    public void ALengthAboveTheLimit_IsRefusedBeforeAnythingIsWritten()
    {
        using var dir = new TempDir();
        string txt = Input(dir, RandomInput(101, 4, 1));
        var o = new PackedTextWriterOptions { MaxSymbols = 100 };
        NotSupportedException? e = Assert.Throws<NotSupportedException>(() => PackedTextWriter.WriteFile(txt, dir.File("p.pk"), o));
        Assert.That(e!.Message, Does.Contain("101").And.Contain("100"));
        Assert.Throws<NotSupportedException>(() => PackedTextWriter.Scan(txt, o));
        Assert.That(Directory.GetFiles(dir.Path).Select(Path.GetFileName), Is.EquivalentTo(new[] { "in.txt" }));
        PackedTextWriter.WriteFile(Input(dir, RandomInput(100, 4, 1), "ok.txt"), dir.File("ok.pk"), o);
        Assert.That(PackedTextFormat.MaxSymbols, Is.EqualTo(uint.MaxValue - 3L), "the default is the design limit of the 32-bit positions");
    }

    // ---- atomic finish and failure modes

    [Test]
    public void TheTarget_AppearsOnlyWhenTheFileIsComplete_AndNoTemporaryFileIsLeft()
    {
        using var dir = new TempDir();
        byte[] input = RandomInput(50_000, 4, 3);
        string txt = Input(dir, input);
        string pk = dir.File("p.pk");
        var seen = new List<string[]>();
        PackedTextWriter.WriteFile(txt, pk, new PackedTextWriterOptions { BetweenPasses = () => seen.Add(Directory.GetFiles(dir.Path).Select(Path.GetFileName).OfType<string>().Order().ToArray()) });
        Assert.That(seen.Single(), Is.EquivalentTo(new[] { "in.txt" }), "between the passes the target does not exist yet");
        Assert.That(Directory.GetFiles(dir.Path).Select(Path.GetFileName), Is.EquivalentTo(new[] { "in.txt", "p.pk" }), "after: the target only");
    }

    [Test]
    public void AFailure_LeavesAnExistingTargetUntouched_AndNoTemporaryFile()
    {
        using var dir = new TempDir();
        string pk = dir.File("p.pk");
        File.WriteAllBytes(pk, [1, 2, 3]);
        byte[] data = RandomInput(50_000, 4, 3);
        string txt = Input(dir, data);

        // the input is replaced after the alphabet scan (same length): detected at the end of the second pass
        InvalidOperationException? e = Assert.Throws<InvalidOperationException>(() => PackedTextWriter.WriteFile(txt, pk, new PackedTextWriterOptions
        {
            BlockSymbols = 1024,
            WindowBytes = 4096,
            BetweenPasses = () => File.WriteAllBytes(txt, RandomInput(50_000, 4, 99)),
        }));
        Assert.That(e!.Message, Does.Contain("changed between the two passes"));
        Assert.That(File.ReadAllBytes(pk), Is.EqualTo(new byte[] { 1, 2, 3 }));
        Assert.That(Directory.GetFiles(dir.Path).Select(Path.GetFileName), Is.EquivalentTo(new[] { "in.txt", "p.pk" }));

        // the input is truncated after the scan: refused when the second pass starts, after the temporary file exists
        File.WriteAllBytes(txt, data);
        Assert.Throws<InvalidOperationException>(() => PackedTextWriter.WriteFile(txt, pk, new PackedTextWriterOptions { BetweenPasses = () => File.WriteAllBytes(txt, RandomInput(4_000, 4, 3)) }));
        Assert.That(File.ReadAllBytes(pk), Is.EqualTo(new byte[] { 1, 2, 3 }));
        Assert.That(Directory.GetFiles(dir.Path).Select(Path.GetFileName), Is.EquivalentTo(new[] { "in.txt", "p.pk" }));

        // grown by a few bytes: the alphabet may be the same, the length differs
        File.WriteAllBytes(txt, data);
        Assert.Throws<InvalidOperationException>(() => PackedTextWriter.WriteFile(txt, pk, new PackedTextWriterOptions { BetweenPasses = () => File.WriteAllBytes(txt, [.. data, (byte)'a']) }));
        Assert.That(File.ReadAllBytes(pk), Is.EqualTo(new byte[] { 1, 2, 3 }));
        Assert.That(Directory.GetFiles(dir.Path).Select(Path.GetFileName), Is.EquivalentTo(new[] { "in.txt", "p.pk" }));
    }

    [Test]
    public void ASuccessfulWrite_ReplacesAnExistingTarget()
    {
        using var dir = new TempDir();
        string pk = dir.File("p.pk");
        File.WriteAllBytes(pk, new byte[10_000]);
        byte[] data = RandomInput(1_000, 4, 3);
        PackedTextWriter.WriteFile(Input(dir, data), pk);
        (byte[] dense, _) = TextOracle.Dense(data);
        Assert.That(File.ReadAllBytes(pk), Is.EqualTo(TextOracle.Pack(dense, 2)));
    }

    [Test]
    public void AMissingDirectory_AndAMissingInput_AreSurfaced_AndNothingIsLeftBehind()
    {
        using var dir = new TempDir();
        string txt = Input(dir, RandomInput(100, 4, 1));
        Assert.Throws<DirectoryNotFoundException>(() => PackedTextWriter.WriteFile(txt, dir.File("no/such/p.pk")));
        Assert.Throws<FileNotFoundException>(() => PackedTextWriter.WriteFile(dir.File("missing.txt"), dir.File("p.pk")));
        Assert.Throws<FileNotFoundException>(() => PackedTextWriter.Scan(dir.File("missing.txt")));
        Assert.Throws<FileNotFoundException>(() => PackedTextWriter.InputLength(dir.File("missing.txt")));
        Assert.That(Directory.GetFiles(dir.Path).Select(Path.GetFileName), Is.EquivalentTo(new[] { "in.txt" }));
        Assert.Throws<ArgumentNullException>(() => PackedTextWriter.WriteFile(null!, dir.File("p.pk")));
        Assert.Throws<ArgumentNullException>(() => PackedTextWriter.WriteFile(txt, null!));
    }

    [Test]
    public void ADiskThatCannotTakeTheFile_IsSurfaced_AsAnIoException_AndTheTargetIsNotCreated()
    {
        if (!File.Exists("/dev/full") || !OperatingSystem.IsLinux()) Assert.Ignore("no /dev/full on this platform");
        foreach (ReservationMode mode in new[] { ReservationMode.Sparse, ReservationMode.Preallocate })
        {
            using var dir = new TempDir();
            string txt = Input(dir, RandomInput(10_000, 4, 1));
            string pk = dir.File("p.pk");
            File.CreateSymbolicLink(pk + ".tmp", "/dev/full"); // a link, so that nothing here can delete the device: the temporary file is the device
            Assert.Throws<IOException>(() => PackedTextWriter.WriteFile(txt, pk, new PackedTextWriterOptions { File = new MappedFileOptions { Reservation = mode } }), mode.ToString());
            Assert.That(File.Exists(pk), Is.False, $"{mode}: no target");
            Assert.That(File.Exists("/dev/full"), Is.True);
        }
    }

    [Test]
    public void AFullDisk_IsRaisedUpFront_ByThePreallocation_AndByTheFreeSpaceCheck_NotAsAFaultOnAMappedPage()
    {
        using var dir = new TempDir();
        string txt = Input(dir, RandomInput(10_000, 4, 1));
        // a source that behaves like a platform without fallocate and with 100 bytes of free space
        var small = new FakeReservationSource { CanPreallocate = false, Free = 100 };
        IOException? e = Assert.Throws<IOException>(() => PackedTextWriter.WriteFile(txt, dir.File("p.pk"), new PackedTextWriterOptions { File = new MappedFileOptions { Reservation = ReservationMode.Preallocate, Source = small } }));
        Assert.That(e!.Message, Does.Contain("space").IgnoreCase);
        Assert.That(Directory.GetFiles(dir.Path).Select(Path.GetFileName), Is.EquivalentTo(new[] { "in.txt" }));
        // a source whose fallocate reports a full disk
        var full = new FakeReservationSource { CanPreallocate = true, PreallocateThrows = new IOException("disk full") };
        Assert.Throws<IOException>(() => PackedTextWriter.WriteFile(txt, dir.File("p.pk"), new PackedTextWriterOptions { File = new MappedFileOptions { Reservation = ReservationMode.Preallocate, Source = full } }));
        Assert.That(Directory.GetFiles(dir.Path).Select(Path.GetFileName), Is.EquivalentTo(new[] { "in.txt" }));
    }

    // ---- a section inside a larger file (what the container writer does), above 2^32 bytes, on a simulated platform without fallocate

    [Test]
    public void PackIntoASectionOfALargerFile_WritesTheSectionOnly_AndTheReaderFindsIt()
    {
        using var dir = new TempDir();
        byte[] input = RandomInput(30_001, 5, 21, 'A');
        string txt = Input(dir, input);
        PackedTextInfo info = PackedTextWriter.Scan(txt);
        const long front = 4096;
        long total = info.TotalBytes;
        string path = dir.File("sesa.tmp");
        using (MappedOutputFile file = MappedOutputFile.Create(path, front + total + 8192))
        {
            using (MappedWindowWriter w = file.CreateWriter(front, total, 1 << 12))
            {
                // the region is dirty (an earlier build under the same name was not deleted): the writer must write every byte of the section
                w.WriteAt(0, Enumerable.Repeat((byte)0xFF, (int)total).ToArray());
                PackedTextWriter.Pack(txt, info, w, new PackedTextWriterOptions { BlockSymbols = 4096 });
                Assert.That(w.Position, Is.EqualTo(total));
            }
            file.Complete();
        }
        byte[] bytes = File.ReadAllBytes(path);
        (byte[] dense, _) = TextOracle.Dense(input);
        Assert.That(bytes[(int)front..(int)(front + total)], Is.EqualTo(TextOracle.Pack(dense, info.Bits)), "the section");
        Assert.That(bytes[..(int)front], Is.All.Zero);
        Assert.That(bytes[(int)(front + total)..], Is.All.Zero);
        using MappedTextReader r = MappedTextReader.OpenRegion(path, front, total, info.N, info.Bits);
        r.Verify(info);
    }

    [Test]
    public void Pack_RefusesATargetThatIsNotTheSizeOfTheSection_OrAnInputThatChanged()
    {
        using var dir = new TempDir();
        byte[] input = RandomInput(1_000, 4, 3);
        string txt = Input(dir, input);
        PackedTextInfo info = PackedTextWriter.Scan(txt);
        using MappedOutputFile file = MappedOutputFile.Create(dir.File("x.bin"), 10_000);
        using (MappedWindowWriter small = file.CreateWriter(0, info.TotalBytes - 1)) Assert.Throws<ArgumentException>(() => PackedTextWriter.Pack(txt, info, small));
        using (MappedWindowWriter big = file.CreateWriter(0, info.TotalBytes + 1)) Assert.Throws<ArgumentException>(() => PackedTextWriter.Pack(txt, info, big));
        using (MappedWindowWriter used = file.CreateWriter(0, info.TotalBytes))
        {
            used.Append([1]);
            Assert.Throws<ArgumentException>(() => PackedTextWriter.Pack(txt, info, used), "the target must be at its start");
        }
        File.WriteAllBytes(txt, RandomInput(1_001, 4, 3));
        using (MappedWindowWriter ok = file.CreateWriter(0, info.TotalBytes)) Assert.Throws<InvalidOperationException>(() => PackedTextWriter.Pack(txt, info, ok), "the length changed since the scan");
        File.WriteAllBytes(txt, RandomInput(1_000, 4, 4));
        using (MappedWindowWriter ok = file.CreateWriter(0, info.TotalBytes))
        {
            InvalidOperationException? e = Assert.Throws<InvalidOperationException>(() => PackedTextWriter.Pack(txt, info, ok));
            Assert.That(e!.Message, Does.Contain("changed"));
        }
    }

    [Test]
    public void OnASimulatedPlatformWithoutFallocateAndWithoutKnownFreeSpace_WriterAndReaderWorkWithoutThrowing()
    {
        using var dir = new TempDir();
        byte[] input = RandomInput(20_000, 4, 9);
        string txt = Input(dir, input);
        foreach (FakeReservationSource source in new[]
        {
            new FakeReservationSource { CanPreallocate = false, Free = -1 },
            new FakeReservationSource { CanPreallocate = false, Free = long.MaxValue },
            new FakeReservationSource { CanPreallocate = false, Free = 20_000 },
        })
        {
            string pk = dir.File("p.pk");
            PackedTextInfo info = PackedTextWriter.WriteFile(txt, pk, new PackedTextWriterOptions { File = new MappedFileOptions { Reservation = ReservationMode.Preallocate, Source = source } });
            using MappedTextReader r = MappedTextReader.Open(pk, info.N, info.Bits);
            r.Verify(info);
            Assert.That(source.PreallocateCalls, Is.Zero, "no fallocate where there is none");
        }
    }

    [Test]
    public void ASectionAtAFileOffsetAboveTwoToThe32_IsWrittenAndRead()
    {
        using var dir = new TempDir();
        byte[] input = RandomInput(9_999, 7, 5, 'a');
        string txt = Input(dir, input);
        PackedTextInfo info = PackedTextWriter.Scan(txt);
        const long front = 5_000_000_000L;
        string path = dir.File("far.tmp");
        using (MappedOutputFile file = MappedOutputFile.Create(path, front + info.TotalBytes, new MappedFileOptions { Reservation = ReservationMode.Sparse }))
        {
            using MappedWindowWriter w = file.CreateWriter(front, info.TotalBytes);
            PackedTextWriter.Pack(txt, info, w);
        }
        using MappedTextReader r = MappedTextReader.OpenRegion(path, front, info.TotalBytes, info.N, info.Bits);
        r.Verify(info);
        File.Delete(path);
    }
}
