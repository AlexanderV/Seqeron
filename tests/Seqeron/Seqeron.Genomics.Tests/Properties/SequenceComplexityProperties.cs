using FsCheck;
using FsCheck.Fluent;
using FsCheck.NUnit;

namespace Seqeron.Genomics.Tests.Properties;

/// <summary>
/// Property-based tests for sequence-complexity measures (SequenceComplexity): compression-based
/// (Lempel-Ziv) complexity, DUST score, k-mer entropy, and the windowed complexity profile.
///
/// Test Units: SEQ-COMPLEX-COMPRESS-001, SEQ-COMPLEX-DUST-001, SEQ-COMPLEX-KMER-001, SEQ-COMPLEX-WINDOW-001
/// </summary>
[TestFixture]
[Category("Property")]
[Category("Complexity")]
public class SequenceComplexityProperties
{
    private static Arbitrary<string> DnaArbitrary(int minLen) =>
        Gen.Elements('A', 'C', 'G', 'T').ArrayOf().Where(a => a.Length >= minLen)
            .Select(a => new string(a)).ToArbitrary();

    // A long, base-cycling sequence used as the high-complexity reference in monotonicity controls.
    private const string Diverse = "ACGTACGTACGTAGCTAGCTTGCAACGTTACGGATCAGTCAGGCATTACGGTACGTAGCT";

    #region SEQ-COMPLEX-COMPRESS-001: R: ratio ≥ 0; M: repetitive → lower ratio; D: deterministic

    // EstimateCompressionRatio is the normalized Lempel-Ziv complexity; lower values mean more
    // repetitive (more compressible) sequences. (Normalized LZ is non-negative; the checklist's
    // (0,1] band is the asymptotic range.)

    /// <summary>INV-1 (R): the compression complexity is finite and non-negative.</summary>
    [FsCheck.NUnit.Property]
    public Property Compression_IsNonNegative()
    {
        return Prop.ForAll(DnaArbitrary(10), seq =>
        {
            double r = SequenceComplexity.EstimateCompressionRatio(seq);
            return (double.IsFinite(r) && r >= 0).Label($"compression ratio {r} must be ≥ 0");
        });
    }

    /// <summary>
    /// INV-2 (M): a homopolymer (maximally repetitive) has a lower compression complexity than a
    /// base-cycling diverse sequence of the same length.
    /// </summary>
    [Test]
    [Category("Property")]
    public void Compression_Repetitive_LowerThanDiverse()
    {
        double repetitive = SequenceComplexity.EstimateCompressionRatio(new string('A', Diverse.Length));
        double diverse = SequenceComplexity.EstimateCompressionRatio(Diverse);
        Assert.That(repetitive, Is.LessThan(diverse), "repetitive sequence must compress more (lower ratio)");
    }

    /// <summary>INV-3 (D): compression complexity is deterministic.</summary>
    [FsCheck.NUnit.Property]
    public Property Compression_IsDeterministic()
    {
        return Prop.ForAll(DnaArbitrary(10), seq =>
            (SequenceComplexity.EstimateCompressionRatio(seq) == SequenceComplexity.EstimateCompressionRatio(seq))
                .Label("EstimateCompressionRatio must be deterministic"));
    }

    #endregion

    #region SEQ-COMPLEX-DUST-001: R: DUST score ≥ 0; M: low-complexity → higher score; D: deterministic

    // CalculateDustScore (Morgulis et al. 2006): higher score = lower complexity.

    /// <summary>INV-1 (R): the DUST score is non-negative.</summary>
    [FsCheck.NUnit.Property]
    public Property Dust_IsNonNegative()
    {
        return Prop.ForAll(DnaArbitrary(10), seq =>
            (SequenceComplexity.CalculateDustScore(seq) >= 0).Label("DUST score must be ≥ 0"));
    }

    /// <summary>INV-2 (M): a homopolymer scores higher (lower complexity) than a diverse sequence.</summary>
    [Test]
    [Category("Property")]
    public void Dust_LowComplexity_HigherScore()
    {
        double homo = SequenceComplexity.CalculateDustScore(new string('A', Diverse.Length));
        double diverse = SequenceComplexity.CalculateDustScore(Diverse);
        Assert.That(homo, Is.GreaterThan(diverse), "low-complexity sequence must have a higher DUST score");
    }

    /// <summary>INV-3 (D): DUST score is deterministic.</summary>
    [FsCheck.NUnit.Property]
    public Property Dust_IsDeterministic()
    {
        return Prop.ForAll(DnaArbitrary(10), seq =>
            (SequenceComplexity.CalculateDustScore(seq) == SequenceComplexity.CalculateDustScore(seq))
                .Label("CalculateDustScore must be deterministic"));
    }

    #endregion

    #region SEQ-COMPLEX-KMER-001: R: entropy ≥ 0; M: more distinct k-mers → higher; P: homopolymer → 0; D: deterministic

    /// <summary>INV-1 (R): k-mer entropy is non-negative.</summary>
    [FsCheck.NUnit.Property]
    public Property KmerEntropy_IsNonNegative()
    {
        return Prop.ForAll(DnaArbitrary(10), seq =>
            (SequenceComplexity.CalculateKmerEntropy(seq, 2) >= -1e-12).Label("k-mer entropy must be ≥ 0"));
    }

    /// <summary>
    /// INV-2 (P + M): a homopolymer has zero k-mer entropy (one distinct k-mer); a diverse sequence
    /// with many distinct k-mers has strictly higher entropy.
    /// </summary>
    [Test]
    [Category("Property")]
    public void KmerEntropy_HomopolymerZero_DiverseHigher()
    {
        double homo = SequenceComplexity.CalculateKmerEntropy(new string('A', 30), 2);
        double diverse = SequenceComplexity.CalculateKmerEntropy(Diverse, 2);
        Assert.Multiple(() =>
        {
            Assert.That(homo, Is.EqualTo(0.0).Within(1e-9), "homopolymer → single k-mer → zero entropy");
            Assert.That(diverse, Is.GreaterThan(homo), "more distinct k-mers → higher entropy");
        });
    }

    /// <summary>INV-3 (D): k-mer entropy is deterministic.</summary>
    [FsCheck.NUnit.Property]
    public Property KmerEntropy_IsDeterministic()
    {
        return Prop.ForAll(DnaArbitrary(10), seq =>
            (SequenceComplexity.CalculateKmerEntropy(seq, 2) == SequenceComplexity.CalculateKmerEntropy(seq, 2))
                .Label("CalculateKmerEntropy must be deterministic"));
    }

    #endregion

    #region SEQ-COMPLEX-WINDOW-001: R: each window score ∈ [0,1]; P: window count = len−w+1; D: deterministic

    private const int WindowSize = 10;

    /// <summary>
    /// INV-1 (R + P): with step 1 there are len−w+1 windows, each with linguistic complexity in [0,1],
    /// non-negative Shannon entropy, and valid coordinates.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property WindowedComplexity_CountAndRange()
    {
        return Prop.ForAll(DnaArbitrary(20), seq =>
        {
            var points = SequenceComplexity.CalculateWindowedComplexity(new DnaSequence(seq), WindowSize, stepSize: 1).ToList();
            bool countOk = points.Count == seq.Length - WindowSize + 1;
            bool rangeOk = points.All(p =>
                p.LinguisticComplexity is >= 0.0 and <= 1.0 + 1e-9 &&
                p.ShannonEntropy >= -1e-12 &&
                p.WindowStart >= 0 && p.WindowEnd < seq.Length && p.WindowStart <= p.WindowEnd);
            return (countOk && rangeOk).Label($"count={points.Count}, expected {seq.Length - WindowSize + 1}");
        });
    }

    /// <summary>INV-2 (D): the windowed complexity profile is deterministic.</summary>
    [FsCheck.NUnit.Property]
    public Property WindowedComplexity_IsDeterministic()
    {
        return Prop.ForAll(DnaArbitrary(20), seq =>
        {
            var a = SequenceComplexity.CalculateWindowedComplexity(new DnaSequence(seq), WindowSize, 1).Select(p => p.LinguisticComplexity).ToList();
            var b = SequenceComplexity.CalculateWindowedComplexity(new DnaSequence(seq), WindowSize, 1).Select(p => p.LinguisticComplexity).ToList();
            return a.SequenceEqual(b).Label("CalculateWindowedComplexity must be deterministic");
        });
    }

    #endregion
}
