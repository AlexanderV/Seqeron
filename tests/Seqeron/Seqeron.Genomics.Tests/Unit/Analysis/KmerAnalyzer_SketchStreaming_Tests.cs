// KMER-DIST-001 (audit round 6, A6-1; B09 for B06) — the sketch builders hash a rolling 2-bit window instead of building a string per
// distinct k-mer. The previous implementation (CountKmers -> HashSet<string> -> hash of the k ASCII bytes) is kept here as the test oracle:
// the new builders must return EQUAL sketches (k, size, canonical, seed, length, hashes, abundances) on random inputs of every kind the old
// code accepted (case, N / IUPAC / U, short and null records, canonical or forward, k <= 32 and > 32), and must not allocate per k-mer.
// Reference: Mash 2.3 Sketch.cpp / MinHashHeap.cpp (only the s smallest hashes are kept), sourmash 4.9.4 (only hashes <= max_hash are kept).

namespace Seqeron.Genomics.Tests.Unit.Analysis;

[TestFixture]
public class KmerAnalyzer_SketchStreaming_Tests
{
    // ---- the old implementation, verbatim (B09 F95 / F138), as oracle -------------------------------------------------------------

    private static ulong OldH1(string kmer, uint seed)
    {
        var bytes = new byte[kmer.Length];
        for (int i = 0; i < bytes.Length; i++)
            bytes[i] = (byte)kmer[i];
        return KmerAnalyzer.MurmurHash3X64_128(bytes, seed).H1;
    }

    private static KmerCountingOptions Options(bool canonical)
        => canonical ? new KmerCountingOptions(Canonical: true) : new KmerCountingOptions(AcgtOnly: true);

    private static MinHashSketch OldMinHash(IEnumerable<string> records, int k, int sketchSize, bool canonical, uint seed)
    {
        bool use64 = k > 16;
        var kmers = new HashSet<string>(StringComparer.Ordinal);
        long length = 0;
        foreach (var record in records)
        {
            var seq = record ?? string.Empty;
            if (seq.Length < k)
                continue;
            length += seq.Length;
            kmers.UnionWith(KmerAnalyzer.CountKmers(seq, k, Options(canonical)).Keys);
        }

        var hashes = kmers.Select(m => OldH1(m, seed)).Select(h => use64 ? h : (uint)h).ToArray();
        Array.Sort(hashes);
        var bottom = new List<ulong>();
        foreach (var h in hashes)
        {
            if (bottom.Count == sketchSize)
                break;
            if (bottom.Count == 0 || bottom[^1] != h)
                bottom.Add(h);
        }

        return new MinHashSketch(k, sketchSize, canonical, seed, use64, length, bottom.ToArray());
    }

    private static FracMinHashSketch OldFrac(string sequence, int k, long scaled, bool canonical, uint seed, bool trackAbundance)
    {
        ulong maxHash = KmerAnalyzer.FracMinHashMaxHash(scaled);
        var counts = KmerAnalyzer.CountKmers(sequence ?? string.Empty, k, Options(canonical));
        var kept = new Dictionary<ulong, long>();
        foreach (var (kmer, count) in counts)
        {
            ulong h = OldH1(kmer, seed);
            if (h <= maxHash && !kept.TryAdd(h, count))
                kept[h] += count;
        }

        var hashes = kept.Keys.ToArray();
        Array.Sort(hashes);
        long[]? abundances = trackAbundance ? hashes.Select(h => kept[h]).ToArray() : null;
        return new FracMinHashSketch(k, scaled, maxHash, canonical, seed, hashes, abundances);
    }

    // ---- random inputs ----------------------------------------------------------------------------------------------------------

    private static string RandomSequence(Random rng, int length)
    {
        string alphabet = rng.Next(6) switch
        {
            0 => "ACGT",
            1 => "ACGTacgt",
            2 => "ACGTN",
            3 => "ACGTUacgtn",
            4 => "ACGTRYKM-*x",
            _ => "ACGT",
        };
        bool repeats = rng.Next(3) == 0;                    // low-complexity / periodic inputs: many repeated and palindromic k-mers
        var chars = new char[length];
        for (int i = 0; i < length; i++)
            chars[i] = repeats && i >= 7 && rng.Next(10) != 0 ? chars[i - 7] : alphabet[rng.Next(alphabet.Length)];
        return new string(chars);
    }

    private static void AssertSame(MinHashSketch actual, MinHashSketch expected, string context)
    {
        Assert.That(actual.K, Is.EqualTo(expected.K), context);
        Assert.That(actual.SketchSize, Is.EqualTo(expected.SketchSize), context);
        Assert.That(actual.Canonical, Is.EqualTo(expected.Canonical), context);
        Assert.That(actual.Seed, Is.EqualTo(expected.Seed), context);
        Assert.That(actual.Use64, Is.EqualTo(expected.Use64), context);
        Assert.That(actual.Length, Is.EqualTo(expected.Length), context);
        Assert.That(actual.Hashes.ToArray(), Is.EqualTo(expected.Hashes.ToArray()), context);
    }

    // ---- equality with the old implementation -----------------------------------------------------------------------------------

    [Test]
    public void CreateMinHashSketch_Equals_TheMaterialisingImplementation_On_1200_RandomRecordSets()
    {
        var rng = new Random(6001);
        int nonEmpty = 0;
        for (int t = 0; t < 1200; t++)
        {
            int k = rng.Next(5) == 0 ? rng.Next(28, 33) : rng.Next(1, 33);
            int s = new[] { 1, 2, 3, 7, 20, 100, 1000, 100000 }[rng.Next(8)];
            bool canonical = rng.Next(2) == 0;
            uint seed = rng.Next(3) == 0 ? KmerAnalyzer.DefaultMashSeed : (uint)rng.NextInt64(0, 1L << 32);
            var records = new List<string?>();
            int count = rng.Next(0, 5);
            for (int r = 0; r < count; r++)
                records.Add(rng.Next(8) == 0 ? null : RandomSequence(rng, rng.Next(0, 5) == 0 ? rng.Next(0, k + 2) : rng.Next(k, 400)));

            var oracle = OldMinHash(records!, k, s, canonical, seed);
            var actual = KmerAnalyzer.CreateMinHashSketch(records!, k, s, canonical, seed);
            AssertSame(actual, oracle, $"t={t} k={k} s={s} canonical={canonical} seed={seed}");
            if (oracle.Hashes.Count > 0)
                nonEmpty++;
        }

        Assert.That(nonEmpty, Is.GreaterThan(600), "many random sets must hold at least one k-mer, else the comparison is vacuous");
    }

    [Test]
    public void CreateMinHashSketch_SingleSequenceOverload_Equals_TheRecordSetOverload()
    {
        var rng = new Random(6002);
        for (int t = 0; t < 200; t++)
        {
            int k = rng.Next(1, 33);
            string seq = RandomSequence(rng, rng.Next(k, 300));
            AssertSame(KmerAnalyzer.CreateMinHashSketch(seq, k, 50), OldMinHash([seq], k, 50, true, KmerAnalyzer.DefaultMashSeed), $"t={t} k={k}");
        }
    }

    [Test]
    public void CreateFracMinHashSketch_Equals_TheMaterialisingImplementation_On_1200_RandomSequences_K_UpTo_60()
    {
        var rng = new Random(6003);
        int nonEmpty = 0, longK = 0;
        for (int t = 0; t < 1200; t++)
        {
            int k = rng.Next(3) == 0 ? rng.Next(33, 61) : rng.Next(1, 33);
            long scaled = new long[] { 1, 1, 2, 3, 5, 20, 1000 }[rng.Next(7)];
            bool canonical = rng.Next(2) == 0;
            bool abundance = rng.Next(2) == 0;
            uint seed = rng.Next(3) == 0 ? KmerAnalyzer.DefaultSourmashSeed : (uint)rng.NextInt64(0, 1L << 32);
            string? seq = rng.Next(20) == 0 ? null : RandomSequence(rng, rng.Next(0, 6) == 0 ? rng.Next(0, k + 2) : rng.Next(k, 500));

            var oracle = OldFrac(seq!, k, scaled, canonical, seed, abundance);
            var actual = KmerAnalyzer.CreateFracMinHashSketch(seq!, k, scaled, canonical, seed, abundance);
            string context = $"t={t} k={k} scaled={scaled} canonical={canonical} seed={seed} abundance={abundance}";
            Assert.That(actual.Hashes.ToArray(), Is.EqualTo(oracle.Hashes.ToArray()), context);
            Assert.That(actual.Abundances?.ToArray(), Is.EqualTo(oracle.Abundances?.ToArray()), context);
            Assert.That(actual.MaxHash, Is.EqualTo(oracle.MaxHash), context);
            Assert.That((actual.K, actual.Scaled, actual.Canonical, actual.Seed), Is.EqualTo((oracle.K, oracle.Scaled, oracle.Canonical, oracle.Seed)), context);
            if (oracle.Hashes.Count > 0)
            {
                nonEmpty++;
                if (k > 32)
                    longK++;
            }
        }

        Assert.That(nonEmpty, Is.GreaterThan(450));
        Assert.That(longK, Is.GreaterThan(100), "the k > 32 path (byte buffers) must be exercised");
    }

    [Test]
    public void AcgtCode_Accepts_Exactly_TheLetters_ToUpperInvariantMapsToACGT()
    {
        // the streaming reader reads the symbols itself; CountKmers upper-cases first: the two must agree on every UTF-16 code unit
        for (int c = 0; c <= char.MaxValue; c++)
        {
            string seq = new string((char)c, 1) + "ACGTAC";
            bool oldAccepts = KmerAnalyzer.CountKmers(seq, 7, new KmerCountingOptions(AcgtOnly: true)).Count > 0;
            // 7-mers of "<c>ACGTAC" exist only when c is an ACGT letter (the window of the whole string)
            bool newAccepts = KmerAnalyzer.CreateFracMinHashSketch(seq, 7, 1, canonical: false).Hashes.Count > 0;
            Assert.That(newAccepts, Is.EqualTo(oldAccepts), $"U+{c:X4}");
        }
    }

    [Test]
    public void Sketch_Of_PalindromicAndHomopolymerInputs_Equals_TheOracle()
    {
        foreach (string seq in new[] { new string('A', 200), new string('a', 90) + new string('T', 90), "ACGT" + "ACGT" + "ACGT" + "ACGTACGTACGTACGT", "ACGTACGTAC" })
        {
            foreach (int k in new[] { 1, 2, 4, 5, 16, 17, 31, 32 })
                foreach (bool canonical in new[] { true, false })
                    AssertSame(KmerAnalyzer.CreateMinHashSketch(seq, k, 1000, canonical), OldMinHash([seq], k, 1000, canonical, 42), $"k={k} canonical={canonical}");
        }
    }

    // ---- resources and cancellation ---------------------------------------------------------------------------------------------

    private static string RandomGenome(int length)
    {
        ulong x = 12345;
        return string.Create(length, 0, (span, _) =>
        {
            for (int i = 0; i < span.Length; i++)
            {
                x = unchecked(x * 6364136223846793005UL + 1442695040888963407UL);
                span[i] = "ACGT"[(int)(x >> 62)];
            }
        });
    }

    [Test]
    public void CreateMinHashSketch_Of_A_4Mbp_Genome_Allocates_No_String_Per_Kmer()
    {
        string genome = RandomGenome(4_000_000);
        KmerAnalyzer.CreateMinHashSketch(genome.Substring(0, 5000), 21);                 // JIT / static initialisation outside the measurement
        long before = GC.GetAllocatedBytesForCurrentThread();
        var sketch = KmerAnalyzer.CreateMinHashSketch(genome, 21, 1000);
        long allocated = GC.GetAllocatedBytesForCurrentThread() - before;

        Assert.That(sketch.Hashes.Count, Is.EqualTo(1000));
        Assert.That(sketch.Length, Is.EqualTo(4_000_000));
        // the old implementation allocated about 250 bytes per distinct k-mer (~1 GB here); the bounded set allocates O(s log(n / s)) nodes
        Assert.That(allocated, Is.LessThan(16L * 1024 * 1024), $"allocated {allocated:N0} bytes");
    }

    [Test]
    public void CreateFracMinHashSketch_Of_A_4Mbp_Genome_Stores_Only_The_Retained_Hashes()
    {
        string genome = RandomGenome(4_000_000);
        KmerAnalyzer.CreateFracMinHashSketch(genome.Substring(0, 5000), 31, 1000);
        long before = GC.GetAllocatedBytesForCurrentThread();
        var sketch = KmerAnalyzer.CreateFracMinHashSketch(genome, 31, 1000);
        long allocated = GC.GetAllocatedBytesForCurrentThread() - before;

        Assert.That(sketch.Hashes.Count, Is.InRange(3500, 4500));                         // 4e6 / 1000 = 4000 retained hashes
        Assert.That(allocated, Is.LessThan(16L * 1024 * 1024), $"allocated {allocated:N0} bytes");
    }

    [Test]
    public void CreateMinHashSketch_ObservesTheToken_BetweenRecordsAndInsideALongRecord()
    {
        using var already = new CancellationTokenSource();
        already.Cancel();
        Assert.Throws<OperationCanceledException>(() => KmerAnalyzer.CreateMinHashSketch("ACGTACGTACGTACGTACGTACGT", 5, 10, true, 42, already.Token));
        Assert.Throws<OperationCanceledException>(() => KmerAnalyzer.CreateFracMinHashSketch("ACGTACGTACGTACGTACGTACGT", 5, 1, true, 42, false, already.Token));

        // cancelled while the SECOND record is being read: the poll at the start of the next record (and the one inside a long record) throws
        using var cts = new CancellationTokenSource();
        IEnumerable<string> Records()
        {
            yield return RandomGenome(10_000);
            yield return RandomGenome(10_000);
            cts.Cancel();
            yield return RandomGenome(10_000);
            Assert.Fail("a record after the cancellation was requested must not be reached");
        }

        Assert.Throws<OperationCanceledException>(() => KmerAnalyzer.CreateMinHashSketch(Records(), 21, 100, true, 42, cts.Token));

        // inside one long record: the sink of a counting probe cannot be reached, so the cancelled token must stop the scan at the first poll of the record
        Assert.Throws<OperationCanceledException>(() => KmerAnalyzer.CreateMinHashSketch(RandomGenome(100_000), 21, 100, true, 42, already.Token));
    }
}
