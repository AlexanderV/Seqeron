// SEQ-COMPLEX-001 / SEQ-COMPLEX-WINDOW-001 — completeness audit WP5 (B04 F37–F39).
// Evidence: docs/Evidence/SEQ-COMPLEX-WINDOW-001-Evidence.md, docs/Evidence/SEQ-COMPLEX-001-Evidence.md
// TestSpecs: tests/TestSpecs/SEQ-COMPLEX-001.md, tests/TestSpecs/SEQ-COMPLEX-WINDOW-001.md
// Sources: Troyanskaya et al. (2002) Bioinformatics 18:679 (LC over an alphabet of size a); Rosalind LING (a = 4,
//          sample ATTTGGATT → 0.875); Gabrielian & Bolshoy (1999) Comput Chem 23:263 (word length bounded by W);
//          BBTools BBMap 40.02 jgi/BBDuk.java maskLowEntropy + tracker/EntropyTracker.java (entropyk=5, entropywindow=50,
//          normalisation 1/ln(window k-mers), windows with an undefined base skipped).
// Expected values: exact-Fraction Python brute force (LC), scipy/brute force (windows) and the output of the real
// `bbduk.sh ... entropymask=t` (BBMap 40.02), NOT the implementation.

namespace Seqeron.Genomics.Tests.Unit.Analysis;

[TestFixture]
public class SequenceComplexity_AuditWp5_Tests
{
    #region F37 — LC with a fixed alphabet size

    // Exact fractions from the brute-force definition LC = Σ V_i / Σ min(a^i, N − i + 1).
    [TestCase("ATTTGGATT", 100, 4, 7.0 / 8.0, TestName = "LC_Fixed_RosalindLingSample")]
    [TestCase("AACCAACC", 8, 2, 8.0 / 9.0, TestName = "LC_Fixed_BinaryAlphabet")]
    [TestCase("AACCAACC", 8, 4, 3.0 / 4.0, TestName = "LC_Fixed_SameStringDnaAlphabet")]
    [TestCase("MKVLAAGIVGLLLAA", 15, 20, 9.0 / 10.0, TestName = "LC_Fixed_ProteinAlphabet20")]
    [TestCase("ACGTTGCAACGTTGCAAGT", 20, 4, 138.0 / 173.0, TestName = "LC_Fixed_SuffixTreePath_a4")]
    [TestCase("ACGTTGCAACGTTGCAAGT", 20, 6, 46.0 / 59.0, TestName = "LC_Fixed_SuffixTreePath_a6")]
    [TestCase("GATTACAGATTACAGATTACA", 21, 4, 41.0 / 70.0, TestName = "LC_Fixed_Repeat")]
    [TestCase("AAAAAAAA", 8, 1, 1.0, TestName = "LC_Fixed_UnaryAlphabet")]
    public void CalculateLinguisticComplexity_FixedAlphabet_MatchesBruteForce(string seq, int m, int a, double expected)
    {
        Assert.Multiple(() =>
        {
            Assert.That(SequenceComplexity.CalculateLinguisticComplexity(seq, m, a), Is.EqualTo(expected).Within(1e-15));
            Assert.That(SequenceComplexity.CalculateLinguisticComplexity(seq.ToLowerInvariant(), m, a), Is.EqualTo(expected).Within(1e-15));
            if (seq.All(c => "ACGT".Contains(c)))
                Assert.That(SequenceComplexity.CalculateLinguisticComplexity(new DnaSequence(seq), m, a),
                    Is.EqualTo(expected).Within(1e-15), "DnaSequence (suffix-tree) path");
        });
    }

    [Test]
    public void CalculateLinguisticComplexity_FixedAlphabet4_EqualsInferredForDna()
    {
        // Pure DNA infers a = 4, so the fixed-alphabet overload with a = 4 must agree.
        const string seq = "ATGCTAGCATGCAATGGATTACA";
        foreach (int m in new[] { 1, 3, 6, 10, 13, 100 })
            Assert.That(SequenceComplexity.CalculateLinguisticComplexity(seq, m, 4),
                Is.EqualTo(SequenceComplexity.CalculateLinguisticComplexity(seq, m)), $"m={m}");
    }

    [Test]
    public void CalculateLinguisticComplexity_FixedAlphabet_InvalidArguments_Throw()
    {
        Assert.Multiple(() =>
        {
            Assert.Throws<ArgumentException>(() => SequenceComplexity.CalculateLinguisticComplexity("ACGTN", 5, 4),
                "5 distinct symbols > a = 4");
            Assert.Throws<ArgumentException>(() => SequenceComplexity.CalculateLinguisticComplexity(new DnaSequence("ACGT"), 5, 3));
            Assert.Throws<ArgumentOutOfRangeException>(() => SequenceComplexity.CalculateLinguisticComplexity("ACGT", 5, 0));
            Assert.Throws<ArgumentOutOfRangeException>(() => SequenceComplexity.CalculateLinguisticComplexity("ACGT", 0, 4));
            Assert.Throws<ArgumentNullException>(() => SequenceComplexity.CalculateLinguisticComplexity((DnaSequence)null!, 5, 4));
            Assert.That(SequenceComplexity.CalculateLinguisticComplexity("", 5, 4), Is.EqualTo(0));
            Assert.That(SequenceComplexity.CalculateLinguisticComplexity((string)null!, 5, 4), Is.EqualTo(0));
        });
    }

    [Test]
    public void CalculateLinguisticComplexity_HugeAlphabet_NoOverflow()
    {
        // a^i saturates at N: with a = int.MaxValue every V_max,i = N − i + 1, so LC(ACGT, m = 4) = 10/10.
        Assert.That(SequenceComplexity.CalculateLinguisticComplexity("ACGT", 4, int.MaxValue), Is.EqualTo(1.0));
    }

    #endregion

    #region F38 — windowed LC word-length cap, string overloads, threshold validation

    [Test]
    public void CalculateWindowedComplexity_LcMaxWordLength_ChangesOnlyLc()
    {
        // ACGTACGT (w = 8): all-length LC (m = 8): V = 4,4,4,4,4,3,2,1 (26) / maxima 4,7,6,5,4,3,2,1 (32) = 13/16.
        // Default m = 6 keeps the historical 23/29.
        var seq = new DnaSequence("ACGTACGT");
        var def = SequenceComplexity.CalculateWindowedComplexity(seq, 8, 8).Single();
        var full = SequenceComplexity.CalculateWindowedComplexity(seq, 8, 8, lcMaxWordLength: 8).Single();
        var m1 = SequenceComplexity.CalculateWindowedComplexity(seq, 8, 8, lcMaxWordLength: 1).Single();
        Assert.Multiple(() =>
        {
            Assert.That(def.LinguisticComplexity, Is.EqualTo(23.0 / 29.0).Within(1e-15));
            Assert.That(full.LinguisticComplexity, Is.EqualTo(13.0 / 16.0).Within(1e-15));
            Assert.That(m1.LinguisticComplexity, Is.EqualTo(1.0).Within(1e-15));
            Assert.That(full.ShannonEntropy, Is.EqualTo(def.ShannonEntropy));
        });
    }

    [Test]
    public void CalculateWindowedComplexity_String_SkipsWindowsWithUndefinedBases()
    {
        // w = 8, s = 4 over 24 bp with an N at index 10: windows starting at 4 and 8 contain it and are skipped
        // (BBDuk scores only windows with ns() == 0); lower case and U are defined bases.
        const string seq = "acgtacgtaaNaaaaaacguacgt";
        var pts = SequenceComplexity.CalculateWindowedComplexity(seq, 8, 4).ToList();
        Assert.Multiple(() =>
        {
            Assert.That(pts.Select(p => p.WindowStart), Is.EqualTo(new[] { 0, 12, 16 }));
            Assert.That(pts[0].ShannonEntropy, Is.EqualTo(2.0).Within(1e-12));
            Assert.That(pts[0].LinguisticComplexity, Is.EqualTo(23.0 / 29.0).Within(1e-15));
            Assert.That(pts[2].ShannonEntropy, Is.EqualTo(2.0).Within(1e-12), "ACGUACGT: U counts as T");
        });
    }

    [Test]
    public void CalculateWindowedComplexity_String_EqualsDnaSequenceOverloadOnAcgt()
    {
        const string seq = "ATGCATGCAAAAAAAAAAAAGGCCTTAGGCATCGATCGATTTTTTTACGACGACGACG";
        var a = SequenceComplexity.CalculateWindowedComplexity(new DnaSequence(seq), 12, 5, 7).ToList();
        var b = SequenceComplexity.CalculateWindowedComplexity(seq.ToLowerInvariant(), 12, 5, 7).ToList();
        Assert.That(b, Is.EqualTo(a));
    }

    [Test]
    public void CalculateWindowedComplexity_InvalidArguments_ThrowEagerly()
    {
        Assert.Multiple(() =>
        {
            Assert.Throws<ArgumentOutOfRangeException>(() => SequenceComplexity.CalculateWindowedComplexity(new DnaSequence("ACGT"), 2, 1, 0));
            Assert.Throws<ArgumentOutOfRangeException>(() => SequenceComplexity.CalculateWindowedComplexity("ACGT", 0));
            Assert.Throws<ArgumentOutOfRangeException>(() => SequenceComplexity.CalculateWindowedComplexity("ACGT", 2, 0));
            Assert.Throws<ArgumentNullException>(() => SequenceComplexity.CalculateWindowedComplexity((string)null!));
        });
    }

    [TestCase(double.NaN)]
    [TestCase(double.PositiveInfinity)]
    [TestCase(double.NegativeInfinity)]
    [TestCase(-0.1)]
    public void FindLowComplexityRegions_InvalidThreshold_ThrowsEagerly(double threshold)
    {
        Assert.Multiple(() =>
        {
            Assert.Throws<ArgumentOutOfRangeException>(() => SequenceComplexity.FindLowComplexityRegions(new DnaSequence("ACGTACGT"), 4, threshold));
            Assert.Throws<ArgumentOutOfRangeException>(() => SequenceComplexity.FindLowComplexityRegions("ACGTACGT", 4, threshold));
        });
    }

    [Test]
    public void FindLowComplexityRegions_String_NWindowsNeverFlagged()
    {
        // ATGC×5 + A×20 + N + A×20 + ATGC×5, w = 10, t = 0.5: the 10 windows covering the N are skipped, so the
        // poly-A tract splits into 19..39 and 41..62 (flagged windows 19..30 and 41..53 — all-A windows plus the
        // transition windows "CAAAAAAAAA" / "AAAAAAAAAT" with H = 0.469; scipy/brute force).
        string seq = string.Concat(Enumerable.Repeat("ATGC", 5)) + new string('A', 20) + "N" + new string('A', 20)
                     + string.Concat(Enumerable.Repeat("ATGC", 5));
        var regions = SequenceComplexity.FindLowComplexityRegions(seq, 10, 0.5).ToList();
        Assert.That(regions.Select(r => (r.Start, r.End)), Is.EqualTo(new[] { (19, 39), (41, 62) }));
        Assert.That(regions.Select(r => r.MinEntropy), Is.All.EqualTo(0.0));
    }

    [Test]
    public void FindLowComplexityRegions_String_EqualsDnaSequenceOverloadOnAcgt()
    {
        string seq = string.Concat(Enumerable.Repeat("ATGC", 20)) + new string('A', 64) + string.Concat(Enumerable.Repeat("ATGC", 20));
        var a = SequenceComplexity.FindLowComplexityRegions(new DnaSequence(seq), 20, 0.5).ToList();
        var b = SequenceComplexity.FindLowComplexityRegions(seq.ToLowerInvariant(), 20, 0.5).ToList();
        Assert.That(b, Is.EqualTo(a));
        Assert.That(b.Select(r => (r.Start, r.End)), Is.EqualTo(new[] { (79, 145) }));
    }

    #endregion

    #region F39 — BBDuk maskLowEntropy (k-mer entropy) — values from bbduk.sh 40.02

    private const string S1 = "CGGAGCCTGTTCCTGTACCATTATCTCTTCAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAATACCCTGAAGAGGATCTACAGATGCAAAGC";
    private const string S4 = "CGGCTCACAAGGATGATGGGCCATAGTTAGCTTCGCCAAAACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACCGGAGCCTGTTCCTGTACCATTATCTCTTC";

    // Masked intervals (inclusive) from `bbduk.sh entropy=c entropymask=t entropywindow=w entropyk=k` on each sequence.
    [TestCase(S1, 50, 5, 0.5, 11, 88, TestName = "Bbduk_PolyA_Defaults")]
    [TestCase("cggagccuguuccuguaccauuaucucuucaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaauacccugaagaggaucuacagaugcaaagc",
        50, 5, 0.5, 11, 88, TestName = "Bbduk_LowerCaseRna_SameAsDna")]
    [TestCase(S4, 50, 5, 0.7, 17, 143, TestName = "Bbduk_AcRepeat_Cutoff07")]
    [TestCase(S4, 20, 3, 0.6, 34, 125, TestName = "Bbduk_AcRepeat_w20_k3")]
    [TestCase("CGGCTCACAAGGATGATGGGCCATAGTTAGCTTCGCCAAAAAGAAGAAGAAGAAGAAGAAGAAGAAGAAGTACCCTGAAGAGGATCTACAGATGCAAAGC",
        25, 2, 0.55, 32, 74, TestName = "Bbduk_AagRepeat_w25_k2")]
    public void FindLowEntropyRegionsBbduk_MatchesBbduk(string seq, int w, int k, double cutoff, int start, int end)
    {
        var regions = SequenceComplexity.FindLowEntropyRegionsBbduk(seq, cutoff, w, k);
        Assert.That(regions.Select(r => (r.Start, r.End)), Is.EqualTo(new[] { (start, end) }));
        Assert.That(regions[0].Sequence, Is.EqualTo(seq.Substring(start, end - start + 1)));
    }

    [Test]
    public void FindLowEntropyRegionsBbduk_MinEntropy_IsNormalisedKmerEntropy()
    {
        // Lowest failing-window value −Σ p ln p / ln 46 (w = 50, k = 5; exact double 0.2674965368023891 for S1,
        // 0.1810425967800402 for S4), reported in single precision as BBDuk compares it.
        Assert.Multiple(() =>
        {
            Assert.That(SequenceComplexity.FindLowEntropyRegionsBbduk(S1, 0.5).Single().MinEntropy,
                Is.EqualTo((double)(float)0.2674965368023891));
            Assert.That(SequenceComplexity.FindLowEntropyRegionsBbduk(S4, 0.7).Single().MinEntropy,
                Is.EqualTo((double)(float)0.1810425967800402));
        });
    }

    [Test]
    public void FindLowEntropyRegionsBbduk_UndefinedBaseOrShortRead_NothingMasked()
    {
        // bbduk.sh (defaults, entropy=0.5): an N inside the 50-bp poly-A leaves no N-free failing window → unmasked;
        // a 49-bp read is shorter than the window → unmasked even with entropy=1.
        const string withN = "CGGAGCCTGTTCCTGTACCATTATCTCTTCAAAAAAAAAAAAAAAAAAAAAAAAANAAAAAAAAAAAAAAAAAAAAAAAATACCCTGAAGAGGATCTACAGATGCAAAGC";
        Assert.Multiple(() =>
        {
            Assert.That(SequenceComplexity.FindLowEntropyRegionsBbduk(withN, 0.5), Is.Empty);
            Assert.That(SequenceComplexity.FindLowEntropyRegionsBbduk("GTCCATACAGAAGTGATACGGACTTTTGGGCAGGCCTCATCGGGTGAGT", 1.0), Is.Empty);
            Assert.That(SequenceComplexity.FindLowEntropyRegionsBbduk("", 0.5), Is.Empty);
        });
    }

    [Test]
    public void FindLowEntropyRegionsBbduk_K1_EqualsShannonScanWithRescaledCutoff()
    {
        // BBDuk's 1-mer value is H_bits / log2(w), so FindLowComplexityRegions(w, t) ≡ BBDuk(entropyk = 1,
        // entropy = t / log2 w) (1 600 random cases vs bbduk.sh: 0 mismatches).
        string seq = string.Concat(Enumerable.Repeat("ATGC", 12)) + "AAAAAAAAAAAACAAAAAAAAAAAGAAAAAA" + new string('T', 9)
                     + string.Concat(Enumerable.Repeat("GATC", 12)) + "ACACACACACACAC";
        foreach (int w in new[] { 8, 16, 20 })
            foreach (double t in new[] { 0.5, 1.0, 1.3 })
            {
                var shannon = SequenceComplexity.FindLowComplexityRegions(seq, w, t).Select(r => (r.Start, r.End));
                var bbduk = SequenceComplexity.FindLowEntropyRegionsBbduk(seq, t / Math.Log2(w), w, 1).Select(r => (r.Start, r.End));
                Assert.That(bbduk, Is.EqualTo(shannon), $"w={w} t={t}");
            }
    }

    [Test]
    public void FindLowEntropyRegionsBbduk_InvalidArguments_Throw()
    {
        Assert.Multiple(() =>
        {
            Assert.Throws<ArgumentNullException>(() => SequenceComplexity.FindLowEntropyRegionsBbduk(null!, 0.5));
            Assert.Throws<ArgumentOutOfRangeException>(() => SequenceComplexity.FindLowEntropyRegionsBbduk(S1, -0.01));
            Assert.Throws<ArgumentOutOfRangeException>(() => SequenceComplexity.FindLowEntropyRegionsBbduk(S1, 1.01));
            Assert.Throws<ArgumentOutOfRangeException>(() => SequenceComplexity.FindLowEntropyRegionsBbduk(S1, double.NaN));
            Assert.Throws<ArgumentOutOfRangeException>(() => SequenceComplexity.FindLowEntropyRegionsBbduk(S1, 0.5, 50, 0));
            Assert.Throws<ArgumentOutOfRangeException>(() => SequenceComplexity.FindLowEntropyRegionsBbduk(S1, 0.5, 50, 16));
            Assert.Throws<ArgumentOutOfRangeException>(() => SequenceComplexity.FindLowEntropyRegionsBbduk(S1, 0.5, 5, 5));
        });
    }

    #endregion
}
