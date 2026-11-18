using System.Globalization;
using SuffixTree.Esa.Tests.TestSupport;
using Prod = SuffixTree.Esa.External;
using Proto = SuffixTree.EmProto.Text;
using ProtoDcs = SuffixTree.EmProto.Dcs;
using ProtoPipe = SuffixTree.EmProto.Pipeline;

namespace SuffixTree.EmProto.Tests.External;

/// <summary>
/// G2 of work package 1b.7b: the production writer of the packed text (TEXT codec 4) produces byte-identical files to the prototype PH0 (<c>Prepare</c>), with the same alphabet, counts, width and SHA-256; the production reader
/// returns the same keys, symbols, windows, common prefixes and comparisons as the prototype <c>PackedText</c> on the same files, with the vector and the scalar path; the production mapped writer writes the same bytes as the
/// prototype's, also at file offsets above 2^32. The prototype code is unchanged.
/// </summary>
[TestFixture]
public class TextDifferentialTests
{
    private static byte[] Random(int n, int sigma, ulong seed, int offset = 'a')
    {
        var rng = new SplitMix64(seed);
        var b = new byte[n];
        for (int i = 0; i < n; i++) b[i] = (byte)(offset + rng.NextInt(sigma));
        return b;
    }

    private static byte[] DnaWithN(int n, ulong seed)
    {
        var rng = new SplitMix64(seed);
        byte[] alphabet = "ACGTN"u8.ToArray();
        byte[] iupac = "RYKMSWBDHV"u8.ToArray();
        var b = new byte[n];
        for (int i = 0; i < n; i++) b[i] = alphabet[rng.NextInt(4)];
        for (int run = 0; run < 8; run++)
        {
            int at = rng.NextInt(n - 5000);
            Array.Fill(b, (byte)'N', at, 1 + rng.NextInt(4000));
        }
        for (int k = 0; k < 200; k++) b[rng.NextInt(n)] = iupac[rng.NextInt(iupac.Length)];
        return b;
    }

    private static string Meta(string pk, string key) => File.ReadAllLines(pk + ".meta").First(l => l.StartsWith(key + "=", StringComparison.Ordinal))[(key.Length + 1)..];

    private static void SameAsPrototype(DiffDir dir, string input, string label, int block = 1 << 20, int window = 1 << 20)
    {
        string protoPk = dir.File("proto.pk");
        string prodPk = dir.File("prod.pk");
        ProtoPipe.PrepareResult proto = ProtoPipe.Prepare.Run(input, protoPk, block, window);
        Prod.PackedTextInfo prod = Prod.PackedTextWriter.WriteFile(input, prodPk, new Prod.PackedTextWriterOptions { BlockSymbols = block, WindowBytes = window });
        Assert.That(DiffInputs.SameBytes(protoPk, prodPk), Is.True, $"{label}: the packed files differ");
        Assert.That((prod.N, prod.Bits, prod.Sigma, prod.Sha256), Is.EqualTo((proto.N, proto.Bits, proto.Sigma, proto.Sha256)), label);
        Assert.That(prod.Letters, Is.EqualTo(proto.Letters), label);
        Assert.That(prod.Counts, Is.EqualTo(proto.Counts), label);
        File.Delete(protoPk);
        File.Delete(prodPk);
        File.Delete(protoPk + ".meta");
    }

    // ---- the writer

    [Test]
    public void Writer_EqualsPrototypePrepare_OnSmallInputsOfEverySigma_AndEveryLengthTo40()
    {
        using var dir = new DiffDir();
        string input = dir.File("in.txt");
        foreach (int sigma in new[] { 1, 2, 3, 4, 5, 8, 9, 16, 17, 100, 256 })
        {
            for (int n = 0; n <= 40; n++)
            {
                File.WriteAllBytes(input, Random(n, sigma, (ulong)((sigma * 100) + n), sigma == 256 ? 0 : 5));
                SameAsPrototype(dir, input, $"sigma={sigma} n={n}", 8, 64);
            }
        }
    }

    [Test]
    public void Writer_EqualsPrototypePrepare_OnTheGeneratedTexts_AnNContainingText_AndA256SymbolText()
    {
        using var dir = new DiffDir();
        string input = dir.File("in.txt");
        var cases = new List<(string Name, byte[] Data)>
        {
            ("dna", Random(2_000_003, 4, 1, 'A')),
            ("dna with N runs and IUPAC letters", DnaWithN(2_000_001, 2)),
            ("256 symbols", Random(1_000_000, 256, 3, 0)),
            ("protein", Random(1_000_000, 20, 4, 'A')),
            ("binary", Random(777_777, 2, 5, '0')),
        };
        foreach ((string name, byte[] data) in cases)
        {
            File.WriteAllBytes(input, data);
            SameAsPrototype(dir, input, name);
            SameAsPrototype(dir, input, name + " (small blocks and windows)", 4096, 1000);
        }
    }

    [TestCase("dna_10m")]
    [TestCase("prot_1m")]
    [TestCase("rep_10m")]
    public void Writer_EqualsPrototypePrepare_OnTheBenchmarkData(string name)
    {
        string path = Path.Combine(DiffInputs.DataDir, name + ".txt");
        if (!File.Exists(path)) Assert.Ignore($"{path} is not present (generate it with analysis/bench/gen/gen_inputs.py)");
        using var dir = new DiffDir();
        SameAsPrototype(dir, path, name);
    }

    [Test]
    public void Writer_EqualsTheE9SetupPackedFile_OfDna100m()
    {
        string root = Environment.GetEnvironmentVariable("SEQERON_BENCH_DIR") ?? "/tmp/stbench";
        string text = Path.Combine(root, "e9", "dna_100m.txt");
        string pk = Path.Combine(root, "e9", "dna_100m.pk");
        if (!File.Exists(text) || !File.Exists(pk)) Assert.Ignore("the E9 setup (analysis/em/run_e9_setup.sh) is not present");
        using var dir = new DiffDir();
        Prod.PackedTextInfo prod = Prod.PackedTextWriter.WriteFile(text, dir.File("prod.pk"));
        Assert.That(DiffInputs.SameBytes(pk, dir.File("prod.pk")), Is.True, "the file of the E6 pipeline of the E9 setup");
        Assert.That(prod.N.ToString(CultureInfo.InvariantCulture), Is.EqualTo(Meta(pk, "n")));
        Assert.That(prod.Bits.ToString(CultureInfo.InvariantCulture), Is.EqualTo(Meta(pk, "bits")));
        Assert.That(string.Join(',', prod.Counts), Is.EqualTo(Meta(pk, "counts")));
        Assert.That(prod.Sha256, Is.EqualTo(Meta(pk, "sha256")));
        Assert.That(Convert.ToHexString(prod.Letters), Is.EqualTo(Meta(pk, "letters_hex")));
    }

    // ---- the reader

    private static void SameReads(string pk, long n, int bits, int positions, string label)
    {
        using var mapped = new Proto.MappedFile(pk, writable: false);
        unsafe
        {
            using var proto = new Proto.PackedText(mapped.Pointer, n, bits, mapped, ownsOwner: false);
            using Prod.MappedTextReader prod = Prod.MappedTextReader.Open(pk, n, bits);
            var rng = new SplitMix64(11);
            for (int round = 0; round < positions; round++)
            {
                long i = (long)(rng.Next() % (ulong)(n + 1));
                if (proto.Key(i) != prod.Key(i)) Assert.Fail($"{label}: Key({i}) differs");
                if (proto.Symbol(i) != prod.Symbol(i)) Assert.Fail($"{label}: Symbol({i}) differs");
                if (round % 8 == 0 && proto.Window(i) != prod.Window(i)) Assert.Fail($"{label}: Window({i}) differs");
            }
            // the same on the scalar path of both and on the vector path of both (the keys do not depend on it; the common prefix does)
            foreach (bool vectors in new[] { false, true })
            {
                proto.UseVectors = vectors && System.Runtime.Intrinsics.X86.Avx2.IsSupported;
                prod.UseVectors = vectors;
                Assert.That(prod.UseVectors, Is.EqualTo(proto.UseVectors), $"{label}: the flag is guarded by the hardware in both");
                for (int round = 0; round < Math.Max(1000, positions / 20); round++)
                {
                    long i = (long)(rng.Next() % (ulong)n);
                    long j = round % 3 == 0 ? Math.Min(n - 1, i + 1 + (long)(rng.Next() % 5000UL)) : (long)(rng.Next() % (ulong)n);
                    long max = n - Math.Max(i, j);
                    long cap = round % 2 == 0 ? max : Math.Min(max, 1 + (long)(rng.Next() % 400UL));
                    long s1 = 0, s2 = 0;
                    long a = proto.Lce(i, j, cap, ref s1);
                    long b = prod.Lce(i, j, cap, ref s2);
                    if (a != b || s1 != s2) Assert.Fail($"{label}: Lce({i},{j},{cap}) = {b} in {s2} steps, prototype {a} in {s1} steps (vectors {vectors})");
                    if (i != j)
                    {
                        int c1 = proto.CompareSuffixes(i, j, 0, ref s1, out long l1);
                        int c2 = prod.CompareSuffixes(i, j, 0, ref s2, out long l2);
                        if (Math.Sign(c1) != Math.Sign(c2) || l1 != l2) Assert.Fail($"{label}: CompareSuffixes({i},{j}) = {c2}/{l2}, prototype {c1}/{l1} (vectors {vectors})");
                    }
                }
            }
        }
    }

    [TestCase(1)]
    [TestCase(2)]
    [TestCase(3)]
    [TestCase(4)]
    [TestCase(5)]
    [TestCase(8)]
    public void Reader_EqualsPrototypePackedText_At10To6RandomPositionsPerWidth_VectorAndScalar(int bits)
    {
        using var dir = new DiffDir();
        string input = dir.File("in.txt");
        File.WriteAllBytes(input, Random(3_000_000, 1 << bits, (ulong)bits, bits == 8 ? 0 : 'A'));
        // a repeat in the text: the common prefixes are long there
        byte[] data = File.ReadAllBytes(input);
        Array.Copy(data, 100_000, data, 2_000_000, 400_000);
        File.WriteAllBytes(input, data);
        Prod.PackedTextInfo info = Prod.PackedTextWriter.WriteFile(input, dir.File("t.pk"));
        Assert.That(info.Bits, Is.EqualTo(bits));
        SameReads(dir.File("t.pk"), info.N, info.Bits, 1_000_000, $"b={bits}");
    }

    [TestCase("dna_10m")]
    [TestCase("rep_10m")]
    [TestCase("prot_1m")]
    public void Reader_EqualsPrototypePackedText_OnTheBenchmarkData(string name)
    {
        string path = Path.Combine(DiffInputs.DataDir, name + ".txt");
        if (!File.Exists(path)) Assert.Ignore($"{path} is not present");
        using var dir = new DiffDir();
        Prod.PackedTextInfo info = Prod.PackedTextWriter.WriteFile(path, dir.File("t.pk"));
        SameReads(dir.File("t.pk"), info.N, info.Bits, 1_000_000, name);
        using Prod.MappedTextReader r = Prod.MappedTextReader.Open(dir.File("t.pk"), info.N, info.Bits);
        r.Verify(info);
    }

    [Test]
    public void Reader_EqualsPrototypePackedText_OnTheE9PackedFileOfDna100m()
    {
        string root = Environment.GetEnvironmentVariable("SEQERON_BENCH_DIR") ?? "/tmp/stbench";
        string pk = Path.Combine(root, "e9", "dna_100m.pk");
        if (!File.Exists(pk)) Assert.Ignore("the E9 setup is not present");
        long n = long.Parse(Meta(pk, "n"), CultureInfo.InvariantCulture);
        int bits = int.Parse(Meta(pk, "bits"), CultureInfo.InvariantCulture);
        SameReads(pk, n, bits, 1_000_000, "dna_100m");
    }

    // ---- the mapped writer

    [Test]
    public void MappedWriter_WritesTheSameBytesAsThePrototype_ForTheSameSequenceOfWrites()
    {
        using var dir = new DiffDir();
        var rng = new SplitMix64(4);
        foreach (int window in new[] { 61, 4096, 1 << 20 })
        {
            const int length = 300_000;
            var script = new List<(long At, byte[] Data, bool Append)>();
            long next = 0;
            for (int k = 0; k < 400; k++)
            {
                bool append = k % 3 != 0 && next < length - 2000;
                int take = 1 + rng.NextInt(1500);
                long at = append ? next : rng.NextInt(length - take);
                var data = new byte[take];
                for (int i = 0; i < take; i++) data[i] = (byte)rng.NextInt(256);
                script.Add((at, data, append));
                if (append) next += take;
            }
            string a = dir.File($"a{window}.bin"), b = dir.File($"b{window}.bin");
            using (var proto = new ProtoDcs.MappedWindowWriter(a, length, window))
            using (var prod = new Prod.MappedWindowWriter(b, length, window))
            {
                foreach ((long at, byte[] data, bool append) in script)
                {
                    if (append)
                    {
                        proto.Append(data);
                        prod.Append(data);
                    }
                    else
                    {
                        proto.WriteAt(at, data);
                        prod.WriteAt(at, data);
                    }
                    Assert.That(prod.Position, Is.EqualTo(proto.Position));
                }
            }
            Assert.That(DiffInputs.SameBytes(a, b), Is.True, $"window {window}");
        }
    }

    [Test]
    public void MappedWriter_WritesTheSameBytesAsThePrototype_AtFileOffsetsAboveTwoToThe32_OnSparseFiles()
    {
        using var dir = new DiffDir();
        const long length = 5_000_000_000L;
        string a = dir.File("a.bin"), b = dir.File("b.bin");
        (long At, byte[] Data)[] writes = [(4_294_967_290L, [5, 5, 5, 5, 5, 5, 5, 5, 5, 5]), (4_900_000_000L, [1, 2, 3, 4]), (4_294_967_295L, [9, 8]), (length - 2, [7, 7])];
        using (var proto = new ProtoDcs.MappedWindowWriter(a, length, 1 << 20))
        using (var prod = new Prod.MappedWindowWriter(b, length, 1 << 20))
        {
            foreach ((long at, byte[] data) in writes)
            {
                proto.WriteAt(at, data);
                prod.WriteAt(at, data);
            }
        }
        foreach (string p in new[] { a, b })
        {
            using var fs = new FileStream(p, FileMode.Open, FileAccess.Read);
            Assert.That(fs.Length, Is.EqualTo(length));
            var buf = new byte[10];
            fs.Position = 4_294_967_290L;
            fs.ReadExactly(buf);
            Assert.That(buf, Is.EqualTo(new byte[] { 5, 5, 5, 5, 5, 9, 8, 5, 5, 5 }), p);
            fs.Position = 4_900_000_000L;
            fs.ReadExactly(buf.AsSpan(0, 4));
            Assert.That(buf[..4], Is.EqualTo(new byte[] { 1, 2, 3, 4 }));
            fs.Position = length - 2;
            fs.ReadExactly(buf.AsSpan(0, 2));
            Assert.That(buf[..2], Is.EqualTo(new byte[] { 7, 7 }));
        }
        File.Delete(a);
        File.Delete(b);
    }
}
