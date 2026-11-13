using SuffixTree.EmProto.Harness;
using SuffixTree.EmProto.Ph4;
using SuffixTree.EmProto.Ph5;
using SuffixTree.EmProto.Tests.Build;
using SuffixTree.EmProto.Tests.Ph4;
using SuffixTree.EmProto.Text;
using SuffixTree.Esa;
using SuffixTree.Esa.Tests.TestSupport;

namespace SuffixTree.EmProto.Tests.Ph5;

/// <summary>
/// PH5 of ADR-004 (the lcp-interval scan with the external stack; layout A and layout F) against the in-RAM structures of the product (gate G1 of the prototype E5): the single-array child table entry by entry in both
/// codecs and as the three plain arrays, PSV and NSV against the LCP pyramid, NODE_COUNT and LRS against the header statistics, the child enumeration through every layout and the logical hash against the Ukkonen tree.
/// </summary>
[TestFixture]
public class ChildScanTests
{
    private string _dir = "";

    [SetUp]
    public void SetUp()
    {
        _dir = Path.Combine(Path.GetTempPath(), "e5-" + Guid.NewGuid().ToString("N"));
        Directory.CreateDirectory(_dir);
    }

    [TearDown]
    public void TearDown() => Directory.Delete(_dir, true);

    private string Scratch => Path.Combine(_dir, "scratch");

    [Test]
    public void Banana_And_Friends_AreEqualToTheInRamStructures()
    {
        foreach (string s in new[] { "a", "aa", "ab", "aba", "banana", "mississippi", "abracadabra", "aaaaaaaaab", "abababababab" })
        {
            ChildReference r = ChildReference.Of(s);
            foreach (ChildCodec codec in Enum.GetValues<ChildCodec>())
                ChildChecks.RunAndCompare(r, ChildChecks.Tiny(Scratch, codec: codec), _dir, $"'{s}' {codec}");
        }
    }

    [Test]
    public void Exhaustive_Binary1To12_Ternary1To6_AllOutputsEqual_TinyBlocksAndWindows()
    {
        int count = 0;
        foreach ((int sigma, int max) in new[] { (2, 12), (3, 6) })
        {
            foreach (string s in ExhaustiveStrings.All(sigma, max))
            {
                if (s.Length == 0) continue;
                ChildReference r = ChildReference.Of(s);
                int k = count++;
                // block 2..5, window 16..64 bytes, 1..3 slots rotate through the strings so that every spill, refill and eviction path is met
                ChildChecks.RunAndCompare(r, ChildChecks.Tiny(Scratch, blockEntries: 2 + (k % 4), windowBytes: 16 << (k % 3), slots: 1 + (k % 3), codec: ChildCodec.Tagged), _dir, $"'{s}' tagged");
                ChildChecks.RunAndCompare(r, ChildChecks.Tiny(Scratch, blockEntries: 2 + ((k + 1) % 4), windowBytes: 16 << ((k + 1) % 3), slots: 1 + ((k + 2) % 3), codec: ChildCodec.Packed), _dir, $"'{s}' packed", both: false);
            }
        }
        Assert.That(count, Is.EqualTo(8190 + 1092), "binary 1-12 and ternary 1-6");
    }

    [Test]
    public void Traversal_ThroughEveryLayout_EqualsTheInRamIndex_AndTheUkkonenHash_OnTinyTexts()
    {
        int count = 0;
        var texts = new List<string> { "banana", "mississippi", "abracadabra", "a", "aaaa", "abababab", "aabaabaab" };
        foreach ((int sigma, int max) in new[] { (2, 8), (3, 5) })
            texts.AddRange(ExhaustiveStrings.All(sigma, max).Where(s => s.Length > 0));
        foreach (string s in texts)
        {
            ChildReference r = ChildReference.Of(s);
            string t = Path.Combine(_dir, "t.bin"), p = Path.Combine(_dir, "p.bin");
            ChildBuildStats st = ChildChecks.RunAndCompare(r, ChildChecks.Tiny(Scratch, blockEntries: 2 + (count % 3)), _dir, $"'{s}'");
            File.Copy(Path.Combine(_dir, "child.bin"), t, true);
            ChildChecks.RunAndCompare(r, ChildChecks.Tiny(Scratch, codec: ChildCodec.Packed), _dir, $"'{s}' packed");
            File.Copy(Path.Combine(_dir, "child.bin"), p, true);
            ChildChecks.CompareTraversals(r, t, p, Path.Combine(_dir, "nsv.bin"), $"'{s}'", ukkonen: true);
            Assert.That(st.Backward, Is.Not.Null);
            count++;
        }
        Assert.That(count, Is.GreaterThan(500));
    }

    [Test]
    public void AdversarialFamilies_Differential_TinyBlocksAndWindows()
    {
        var rng = new Random(20241105);
        int cases = 0;
        foreach (int n in new[] { 200, 1000, 5000 })
        {
            foreach (ulong seed in new ulong[] { 1, 2 })
            {
                foreach ((string name, string text) in StringFamilies.Adversarial(n, seed))
                {
                    if (text.Length == 0) continue;
                    ChildReference r = ChildReference.Of(text);
                    ChildCodec codec = rng.Next(2) == 0 ? ChildCodec.Tagged : ChildCodec.Packed;
                    ChildChecks.RunAndCompare(r, ChildChecks.Tiny(Scratch, blockEntries: 2 + rng.Next(49), windowBytes: 16 << rng.Next(7), slots: 1 + rng.Next(4), codec: codec), _dir, $"{name}-n{n} {codec}", both: codec == ChildCodec.Tagged);
                    cases++;
                }
            }
        }
        Assert.That(cases, Is.GreaterThan(150));
    }

    [Test]
    public void BothCodecsAndLayoutF_OnTheFamiliesAt100000()
    {
        foreach ((string name, string text) in StringFamilies.Adversarial(100_000, 7))
        {
            if (text.Length == 0) continue;
            ChildReference r = ChildReference.Of(text);
            ChildChecks.RunAndCompare(r, ChildChecks.Tiny(Scratch, blockEntries: 4096, windowBytes: 1 << 16, slots: 4, codec: ChildCodec.Tagged), _dir, $"{name} tagged");
            ChildChecks.RunAndCompare(r, ChildChecks.Tiny(Scratch, blockEntries: 4096, windowBytes: 1 << 16, slots: 4, codec: ChildCodec.Packed), _dir, $"{name} packed", both: false);
        }
    }

    [Test]
    public void TheDeepFamiliesAt1000000_AreEqual()
    {
        foreach ((string name, string text) in new[] { ("all-same", StringFamilies.AllSame(1_000_000)), ("fibonacci", StringFamilies.Fibonacci(1_000_000)), ("thue-morse", StringFamilies.ThueMorse(1_000_000)) })
        {
            ChildReference r = ChildReference.Of(text);
            ChildChecks.RunAndCompare(r, ChildChecks.Tiny(Scratch, blockEntries: 8192, windowBytes: 1 << 16, slots: 4), _dir, name);
        }
    }

    [Test]
    public void SpillVolume_OfTheAllSameText_IsTheBlockGranularityOfSixteenBytesTimesMaxLcpPlusOne()
    {
        foreach ((int n, int block) in new[] { (1000, 7), (5000, 64), (20_000, 100), (4096, 4096) })
        {
            ChildReference r = ChildReference.Of(StringFamilies.AllSame(n));
            ChildBuildStats st = ChildChecks.RunAndCompare(r, ChildChecks.Tiny(Scratch, blockEntries: block, windowBytes: 4096, slots: 4), _dir, $"a^{n}", both: false);
            ChildScanResult f = st.Forward;
            Assert.That(f.MaxLcp, Is.EqualTo(n - 1), "max LCP of a^n");
            Assert.That(f.StackPeak, Is.EqualTo(n), "height = max LCP + 1 (every entry pushes, nothing pops before the end)");
            long blocks = Math.Max(0, ((n - 1) / block) - 1);
            Assert.That(f.Spill.BlocksWritten, Is.EqualTo(blocks));
            Assert.That(f.Spill.BytesWritten, Is.EqualTo(blocks * block * 16));
            Assert.That(f.Spill.PeakFileBytes, Is.EqualTo(blocks * block * 16));
            Assert.That(f.Spill.BlocksRead, Is.EqualTo(blocks), "every spilled block is read back once");
            long expected = 16L * (f.MaxLcp + 1);
            Assert.That(expected - f.Spill.BytesWritten, Is.InRange(16L * (block + 1), 16L * 2 * block), "spill = 16 B x (max LCP + 1) up to two blocks");
        }
    }

    [Test]
    public void NoSpill_WhenTheStackStaysInItsTwoBlocks()
    {
        ChildReference r = ChildReference.Of(StringFamilies.RandomDna(5000, 3));
        ChildBuildStats st = ChildChecks.RunAndCompare(r, ChildChecks.Tiny(Scratch, blockEntries: 1024), _dir, "random", both: false);
        Assert.That(st.Forward.Spill.BlocksWritten, Is.Zero);
        Assert.That(st.Forward.StackPeak, Is.LessThan(40));
    }

    [Test]
    public void TheHook_FedInRandomChunks_GivesTheFileModeOutputs()
    {
        var rng = new Random(3);
        foreach (string text in new[] { StringFamilies.Fibonacci(3000), StringFamilies.RandomDna(2000, 4), StringFamilies.AllSame(500), "mississippi" })
        {
            ChildReference r = ChildReference.Of(text);
            uint[] dense = ChildFiles.Dense(r.Lcp);
            string child = Path.Combine(_dir, "hook.child"), psv = Path.Combine(_dir, "hook.psv");
            using (ChildScan scan = ChildScan.Create(r.N, new ChildScanOutputs { ChildPath = child, PsvPath = psv }, ChildChecks.Tiny(Scratch, blockEntries: 5)))
            {
                int at = 0;
                while (at < dense.Length)
                {
                    int take = Math.Min(dense.Length - at, 1 + rng.Next(40));
                    scan.Consume(at, dense.AsSpan(at, take));
                    at += take;
                }
                ChildScanResult res = scan.Complete();
                ChildChecks.CompareStats(r, res, "hook");
            }
            ChildChecks.CompareChild(r, child, ChildCodec.Tagged, "hook");
            ChildChecks.CompareArrays(r.Psv, ChildFiles.ReadSmaller(psv), "PSV", "hook");
        }
    }

    [Test]
    public void TheHook_OfPh4_FusedBuild_GivesTheFileModeOutputs()
    {
        foreach (string text in new[] { StringFamilies.Fibonacci(2000), StringFamilies.RandomDna(3000, 8), StringFamilies.AllSame(700) })
        {
            ChildReference r = ChildReference.Of(text);
            DenseCase c = BuilderTestSupport.Dense("fused", text);
            uint[] sa = r.Index.SuffixArray.Select(x => (uint)x).ToArray();
            using PackedText packed = PackedText.Pack(c.Symbols, c.MinBits);
            string child = Path.Combine(_dir, "fused.child"), psv = Path.Combine(_dir, "fused.psv");
            ExtLcpOptions lo = LcpHarness.Options(Path.Combine(_dir, "lcpscratch"), r: 61, threads: 2);
            string saPath = LcpHarness.WriteSa(sa, _dir);
            using (ChildScan scan = ChildScan.Create(r.N, new ChildScanOutputs { ChildPath = child, PsvPath = psv }, ChildChecks.Tiny(Scratch, blockEntries: 3)))
            {
                _ = ExternalLcp.Build(saPath, packed, Path.Combine(_dir, "f.lcp"), Path.Combine(_dir, "f.ovf"), lo, scan);
                ChildChecks.CompareStats(r, scan.Complete(), "fused");
            }
            ChildChecks.CompareChild(r, child, ChildCodec.Tagged, "fused");
            ChildChecks.CompareArrays(r.Psv, ChildFiles.ReadSmaller(psv), "PSV", "fused");
            // the NSV pass reads the LCP file PH4 wrote
            string nsv = Path.Combine(_dir, "fused.nsv");
            using (var file = new CompactLcpFile(Path.Combine(_dir, "f.lcp"), Path.Combine(_dir, "f.ovf"), r.N))
                _ = NsvScan.Run(r.N, file, nsv, ChildChecks.Tiny(Scratch, blockEntries: 3));
            ChildChecks.CompareArrays(r.Nsv, ChildFiles.ReadSmaller(nsv), "NSV", "fused");
        }
    }

    [Test]
    public void Consume_RefusesChunksOutOfOrder_AndCompleteNeedsEveryRank()
    {
        ChildReference r = ChildReference.Of("banana");
        uint[] dense = ChildFiles.Dense(r.Lcp);
        using ChildScan scan = ChildScan.Create(r.N, new ChildScanOutputs { ChildPath = Path.Combine(_dir, "o.child") }, ChildChecks.Tiny(Scratch));
        Assert.Throws<InvalidOperationException>(() => scan.Consume(1, dense.AsSpan(1, 3)), "the first chunk starts at rank 0");
        scan.Consume(0, dense.AsSpan(0, 3));
        Assert.Throws<InvalidOperationException>(() => scan.Consume(4, dense.AsSpan(4, 2)), "a gap");
        Assert.Throws<InvalidOperationException>(() => scan.Consume(2, dense.AsSpan(2, 2)), "an overlap");
        Assert.Throws<InvalidOperationException>(() => scan.Complete(), "ranks 3.. are missing");
        scan.Consume(3, dense.AsSpan(3, dense.Length - 3));
        Assert.Throws<InvalidOperationException>(() => scan.Consume(dense.Length, dense.AsSpan(0, 1)), "more ranks than n + 1");
        _ = scan.Complete();
    }

    [Test]
    public void ACorruptLcp_IsRefused()
    {
        using ChildScan scan = ChildScan.Create(3, new ChildScanOutputs { ChildPath = Path.Combine(_dir, "c.child") }, ChildChecks.Tiny(Scratch));
        Assert.Throws<InvalidDataException>(() => scan.Consume(0, new uint[] { 0, 1, 0, 0 }), "LCP[1] must be 0");
    }

    [Test]
    public void Scratch_IsEmptyAfterTheScan_AndWindowsStayWithinTheirSlots()
    {
        ChildReference r = ChildReference.Of(StringFamilies.Fibonacci(20_000));
        ChildBuildStats st = ChildChecks.RunAndCompare(r, ChildChecks.Tiny(Scratch, blockEntries: 16, windowBytes: 256, slots: 3), _dir, "fib");
        Assert.That(st.Forward.ChildWindows.Opened, Is.GreaterThan(st.Forward.ChildWindows.DistinctWindows - 1));
        Assert.That(st.Forward.ChildWindows.Evictions, Is.GreaterThan(0), "the small window set was recycled");
        Assert.That(st.Forward.Spill.BlocksWritten, Is.GreaterThan(0));
    }
}
