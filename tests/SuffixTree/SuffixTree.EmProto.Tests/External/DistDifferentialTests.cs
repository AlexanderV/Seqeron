using System.Runtime.InteropServices;
using SuffixTree.EmProto.Dcs;
using SuffixTree.Esa;
using SuffixTree.Esa.Tests.TestSupport;
using Prod = SuffixTree.Esa.External;
using Proto = SuffixTree.EmProto.Ph6;
using ProtoBuild = SuffixTree.EmProto.Build;
using ProtoHarness = SuffixTree.EmProto.Harness;

namespace SuffixTree.EmProto.Tests.External;

/// <summary>
/// G2 of work package 1b.7a, part 1: the production <c>Dist</c> (<c>src/SuffixTree/Algorithms/SuffixTree.Esa/External</c>) writes byte-identical ISA files to the prototype <c>IsaByDist</c> and the rank files of
/// <c>SaInverter</c>, with the same statistics, the same explicit and mapped I/O and the same peak scratch. The prototype has no Phi file (PH4 holds Phi per chunk in RAM), so the Phi file is compared with the
/// in-RAM scatter here and with the in-RAM LCP in the Esa tests.
/// </summary>
[TestFixture]
public class DistDifferentialTests
{
    private static void CompareIsa(DiffDir dir, string saPath, long n, long r, int io, int read, int window, string label)
    {
        string scratchProto = dir.File("scratch-proto");
        string scratchProd = dir.File("scratch-prod");
        string isaProto = dir.File("isa.proto");
        string isaProd = dir.File("isa.prod");

        ProtoHarness.IoCounters.Reset();
        ProtoBuild.ScratchTracker.Reset();
        MappedIo.Reset();
        Proto.DistStats a = Proto.IsaByDist.Build(saPath, n, isaProto, new Proto.DistOptions { ScratchDir = scratchProto, ChunkPositions = r, IoBufferBytes = io, ReadBytes = read, WindowBytes = window });
        long protoRead = ProtoHarness.IoCounters.ReadBytes;
        long protoWrite = ProtoHarness.IoCounters.WriteBytes;
        long protoMapped = MappedIo.WriteBytes;
        long protoScratchPeak = ProtoBuild.ScratchTracker.Peak;

        var counters = new Prod.IoCounters();
        Prod.DistStats b = Prod.Dist.BuildIsa(saPath, n, isaProd, new Prod.DistOptions { ScratchDir = scratchProd, ChunkEntries = r, IoBufferBytes = io, ReadBytes = read, WindowBytes = window, Counters = counters });

        Assert.That(DiffInputs.SameBytes(isaProto, isaProd), Is.True, $"{label}: ISA files differ");
        Assert.That((b.Chunks, b.TupleBytesWritten, b.TupleBytesRead), Is.EqualTo((a.Chunks, a.TupleBytesWritten, a.TupleBytesRead)), $"{label}: stats");
        Assert.That(counters.ReadBytes, Is.EqualTo(protoRead), $"{label}: explicit read bytes");
        Assert.That(counters.WriteBytes, Is.EqualTo(protoWrite), $"{label}: explicit write bytes");
        Assert.That(counters.MappedWriteBytes, Is.EqualTo(protoMapped), $"{label}: mapped write bytes");
        Assert.That(counters.ScratchPeakBytes, Is.EqualTo(protoScratchPeak), $"{label}: peak scratch");
        Assert.That(Directory.GetFileSystemEntries(scratchProd), Is.Empty);
        File.Delete(isaProto);
        File.Delete(isaProd);
    }

    [Test]
    public void Isa_SmallTexts_EveryOptionCombination_EqualsThePrototype()
    {
        var rng = new SplitMix64(11);
        int cases = 0;
        foreach (string text in new[] { "", "a", "banana", "mississippi", StringFamilies.AllSame(300), StringFamilies.Fibonacci(500), StringFamilies.ThueMorse(500), StringFamilies.RandomDna(800, 3), StringFamilies.Random(400, 20, 4) })
        {
            using var dir = new DiffDir();
            string sa = DiffInputs.WriteSa(dir, text);
            foreach (long r in new long[] { 1, 3, 17, 100, 1_000_000 })
            {
                foreach (int io in new[] { 8, 24, 64, 4096 })
                {
                    int read = 4 * (1 + rng.NextInt(300));
                    int window = new[] { 1, 7, 64, 1 << 20 }[rng.NextInt(4)];
                    CompareIsa(dir, sa, text.Length, r, io, read, window, $"'{(text.Length > 12 ? text[..12] : text)}'({text.Length}) R={r} io={io} read={read} window={window}");
                    cases++;
                }
            }
        }
        Assert.That(cases, Is.EqualTo(9 * 5 * 4));
    }

    [Test]
    public void Isa_AdversarialFamilies_OneMillion_EqualsThePrototype()
    {
        foreach ((string name, string text) in DiffInputs.Families(1_000_000))
        {
            using var dir = new DiffDir();
            string sa = DiffInputs.WriteSa(dir, text);
            CompareIsa(dir, sa, text.Length, 266_667, 64 * 1024, 1 << 20, 1 << 20, name + " scale-down R");
            CompareIsa(dir, sa, text.Length, 8_000_000, 64 * 1024, 1 << 20, 1 << 20, name + " real R");
        }
    }

    [TestCase("dna_10m")]
    [TestCase("prot_1m")]
    [TestCase("rep_10m")]
    public void Isa_DataFiles_EqualsThePrototype(string name)
    {
        string text = DiffInputs.DataText(name);
        using var dir = new DiffDir();
        string sa = DiffInputs.WriteSa(dir, text);
        long scaled = Math.Max(1, text.Length * 8_000_000L / 3_000_000_000L * 1); // the 3x10^9 ratio: P = 375
        CompareIsa(dir, sa, text.Length, Math.Max(1, (text.Length + 1 + 374) / 375), Math.Max(8, 64 * 1024 * text.Length / 100_000_000), 1 << 20, 1 << 20, name + " P = 375");
        CompareIsa(dir, sa, text.Length, 8_000_000, 64 * 1024, 1 << 20, 1 << 20, name + " real R");
        Assert.That(scaled, Is.GreaterThan(0));
    }

    // ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    // the rank file of PH1 (SaInverter): the same primitive with a skipping mapper and a fill of zero

    /// <summary>The mapper of the rank file: the entry <c>i</c> of the suffix array of the reduced string, <c>SA[i]</c> a reduced position; real positions go to their slot with the rank <c>i</c>, separators and entry 0 are skipped.</summary>
    private readonly struct RankMapper(SampleGeometry g) : Prod.IDistMapper<uint>
    {
        public bool TryMap(long index, uint sourceValue, out uint key, out uint value)
        {
            key = 0;
            value = 0;
            if (index == 0 || !g.TrySlotOfRed(sourceValue, out long slot)) return false;
            key = (uint)slot;
            value = (uint)index;
            return true;
        }
    }

    [Test]
    public void RankFile_EqualsSaInverter_ForRandomSuffixArraysOfTheReducedString()
    {
        var rng = new SplitMix64(21);
        int cases = 0;
        foreach (int x in new[] { 7, 13 })
        {
            DifferenceCover cover = DifferenceCover.Create(x);
            foreach (long n in new long[] { 0, 1, 20, 93, 400 })
            {
                var g = new SampleGeometry(n, cover);
                int reduced = (int)g.RedLength;
                int[] perm = Enumerable.Range(0, reduced).OrderBy(_ => rng.Next()).ToArray();
                var sa = new uint[reduced + 1];
                sa[0] = (uint)reduced;
                for (int i = 0; i < reduced; i++) sa[i + 1] = (uint)perm[i];
                using var dir = new DiffDir();
                string saPath = dir.File("sa.bin");
                File.WriteAllBytes(saPath, MemoryMarshal.AsBytes(sa.AsSpan()).ToArray());
                foreach ((long chunk, int ioBytes) in new[] { (1L, 16), (5L, 24), (g.Slots + 3, 4096), (33L, 100) })
                {
                    string scratch = dir.File("scratch");
                    Directory.CreateDirectory(scratch);
                    ProtoBuild.ScratchTracker.Reset();
                    ExtDcsOptions o = new() { ScratchDir = scratch, ChunkSlots = chunk, IoBufferBytes = ioBytes, WindowBytes = 64 };
                    SaInverter.Invert(saPath, reduced, g, dir.File("rank.proto"), true, o, "inv", new ProtoHarness.ArenaBudget(1 << 24));

                    var mapper = new RankMapper(g);
                    var counters = new Prod.IoCounters();
                    // SaInverter reads IoBufferBytes / 4 source entries at a time and keeps IoBufferBytes / 8 tuples per chunk
                    Prod.Dist.Scatter(saPath, reduced + 1L, ref mapper, dir.File("rank.prod"), new Prod.DistTarget<uint>(g.Slots, 0u, false),
                        new Prod.DistOptions { ScratchDir = scratch, ChunkEntries = chunk, IoBufferBytes = ioBytes, ReadBytes = ioBytes, WindowBytes = 64, Counters = counters, Arena = ArenaOptions.Guarded(0xCD) });
                    Assert.That(DiffInputs.SameBytes(dir.File("rank.proto"), dir.File("rank.prod")), Is.True, $"X={x} n={n} chunk={chunk} io={ioBytes}");
                    Assert.That(Directory.GetFileSystemEntries(scratch), Is.Empty);
                    Assert.That(counters.ScratchLiveBytes, Is.Zero);
                    cases++;
                }
            }
        }
        Assert.That(cases, Is.EqualTo(40));
    }

    // ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    // Phi: no file in the prototype, so the in-RAM scatter is the oracle

    [Test]
    public void Phi_OfTheFamiliesAndOfTheDataFiles_EqualsTheInRamScatter()
    {
        var texts = new List<(string Name, string Text)>(DiffInputs.Families(1_000_000));
        foreach (string name in new[] { "dna_10m", "prot_1m", "rep_10m" })
        {
            string path = Path.Combine(DiffInputs.DataDir, name + ".txt");
            if (File.Exists(path)) texts.Add((name, File.ReadAllText(path).TrimEnd('\n', '\r')));
        }
        foreach ((string name, string text) in texts)
        {
            using var dir = new DiffDir();
            string saPath = DiffInputs.WriteSa(dir, text);
            uint[] sa = MemoryMarshal.Cast<byte, uint>(File.ReadAllBytes(saPath)).ToArray();
            var phi = new uint[sa.Length - 1];
            for (int i = 1; i < sa.Length; i++) phi[sa[i]] = sa[i - 1];
            Prod.Dist.BuildPhi(saPath, text.Length, dir.File("phi.bin"), new Prod.DistOptions { ScratchDir = dir.File("scratch"), ChunkEntries = Math.Max(1, (text.Length + 374) / 375), IoBufferBytes = 64 * 1024, Arena = ArenaOptions.Guarded(0xCD) });
            Assert.That(File.ReadAllBytes(dir.File("phi.bin")), Is.EqualTo(MemoryMarshal.AsBytes(phi.AsSpan()).ToArray()), name);
        }
    }
}
