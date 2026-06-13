// SEQ-DINUC-001 — Dinucleotide Analysis (frequencies, O/E relative abundance, codon frequencies)
// Evidence: docs/Evidence/SEQ-DINUC-001-Evidence.md
// TestSpec: tests/TestSpecs/SEQ-DINUC-001.md
// Source: Karlin S. (1998). Pervasive properties of the genomic signature. PMC126251.
//         Karlin S, Burge C (1995). Trends Genet 11(7):283-290.
//         Nakamura Y et al. Kazusa Codon Usage Database (CUTG) readme.

using Seqeron.Genomics.Analysis;

namespace Seqeron.Genomics.Tests;

[TestFixture]
public class SequenceStatistics_CalculateDinucleotide_Tests
{
    // Expected values are exact rationals derived from the Karlin odds-ratio
    // definition and count/total frequency definitions in the Evidence, computed
    // independently of the implementation.
    private const double Tolerance = 1e-10;

    #region CalculateDinucleotideRatios

    // M1 — rho_XY = f_XY/(f_X f_Y) for ATGCGCGT (A=1,T=2,G=3,C=2; N=8; 7 dinucleotide positions).
    // Evidence: Karlin PMC126251. rho_GC=rho_CG=(2/7)/((3/8)(2/8))=64/21; rho_AT=(1/7)/((1/8)(2/8))=32/7;
    //           rho_TG=rho_GT=(1/7)/((2/8)(3/8))=(1/7)/((3/8)(2/8))=32/21.
    [Test]
    public void CalculateDinucleotideRatios_AtgcgcgtSequence_ReturnsExactOddsRatios()
    {
        var ratios = SequenceStatistics.CalculateDinucleotideRatios("ATGCGCGT");

        Assert.Multiple(() =>
        {
            Assert.That(ratios["GC"], Is.EqualTo(64.0 / 21.0).Within(Tolerance),
                "rho_GC = (2/7)/((3/8)(2/8)) = 64/21 per Karlin odds ratio (INV-03)");
            Assert.That(ratios["CG"], Is.EqualTo(64.0 / 21.0).Within(Tolerance),
                "rho_CG = (2/7)/((2/8)(3/8)) = 64/21");
            Assert.That(ratios["AT"], Is.EqualTo(32.0 / 7.0).Within(Tolerance),
                "rho_AT = (1/7)/((1/8)(2/8)) = 32/7");
            Assert.That(ratios["TG"], Is.EqualTo(32.0 / 21.0).Within(Tolerance),
                "rho_TG = (1/7)/((2/8)(3/8)) = 32/21");
            Assert.That(ratios["GT"], Is.EqualTo(32.0 / 21.0).Within(Tolerance),
                "rho_GT = (1/7)/((3/8)(2/8)) = 32/21");
        });
    }

    // M7 — no-bias baseline: homopolymer AAAA has f_AA=1 and f_A=1, so rho_AA = 1/(1*1) = 1.0 exactly.
    // Evidence: Karlin PMC126251 (r=1 means no bias, INV-03).
    [Test]
    public void CalculateDinucleotideRatios_Homopolymer_ReturnsOddsRatioOne()
    {
        var ratios = SequenceStatistics.CalculateDinucleotideRatios("AAAA");

        Assert.That(ratios["AA"], Is.EqualTo(1.0).Within(Tolerance),
            "rho_AA = f_AA/(f_A*f_A) = 1/(1*1) = 1.0, the no-bias baseline (INV-03)");
    }

    // S4 — division-by-zero guard: in GATTACA there is no C-before-G context but base C is present;
    //      use a sequence missing a base so an expected product is 0. AATT has no G or C, so any
    //      hypothetical GC would have expected 0; we verify the present dinucleotides and that no
    //      ratio is non-finite. Construct ACGT-only-missing-G: ATAT -> bases A,T only.
    // Evidence: expected=0 guard returns ratio 0 (Evidence corner case).
    [Test]
    public void CalculateDinucleotideRatios_AbsentBaseContext_DoesNotProduceInfinity()
    {
        // "ATAT": A=2,T=2,G=0,C=0. Dinucs: AT,TA,AT. f_AT=2/3, f_A=2/4, f_T=2/4 -> rho_AT=(2/3)/(0.25)=8/3.
        var ratios = SequenceStatistics.CalculateDinucleotideRatios("ATAT");

        Assert.Multiple(() =>
        {
            Assert.That(ratios["AT"], Is.EqualTo(8.0 / 3.0).Within(Tolerance),
                "rho_AT = (2/3)/((2/4)(2/4)) = 8/3");
            Assert.That(ratios.Values, Has.All.Matches<double>(double.IsFinite),
                "no ratio is infinite/NaN even though G and C are absent (expected=0 guard)");
        });
    }

    // S1 — guards: null, empty, and length < 2 all return an empty dictionary.
    [Test]
    public void CalculateDinucleotideRatios_NullEmptyOrTooShort_ReturnsEmpty()
    {
        Assert.Multiple(() =>
        {
            Assert.That(SequenceStatistics.CalculateDinucleotideRatios(null!), Is.Empty,
                "null input yields empty dictionary (input guard)");
            Assert.That(SequenceStatistics.CalculateDinucleotideRatios(string.Empty), Is.Empty,
                "empty input yields empty dictionary (input guard)");
            Assert.That(SequenceStatistics.CalculateDinucleotideRatios("A"), Is.Empty,
                "length < 2 yields empty dictionary (no dinucleotide positions)");
        });
    }

    // C1 — case-insensitivity: lowercase input produces the same ratios as uppercase (M1 values).
    [Test]
    public void CalculateDinucleotideRatios_LowercaseInput_MatchesUppercase()
    {
        var ratios = SequenceStatistics.CalculateDinucleotideRatios("atgcgcgt");

        Assert.That(ratios["GC"], Is.EqualTo(64.0 / 21.0).Within(Tolerance),
            "lowercase is normalized via ToUpperInvariant; rho_GC equals the uppercase value 64/21");
    }

    #endregion

    #region CalculateDinucleotideFrequencies

    // M2 — normalized dinucleotide frequencies for ATGCGCGT (7 positions): count/(N-1).
    // Evidence: Karlin PMC126251 (f_XY normalized). GC=CG=2/7; AT=TG=GT=1/7.
    [Test]
    public void CalculateDinucleotideFrequencies_AtgcgcgtSequence_ReturnsExactFrequencies()
    {
        var freq = SequenceStatistics.CalculateDinucleotideFrequencies("ATGCGCGT");

        Assert.Multiple(() =>
        {
            Assert.That(freq["GC"], Is.EqualTo(2.0 / 7.0).Within(Tolerance),
                "f_GC = 2/7 (2 of 7 dinucleotide positions)");
            Assert.That(freq["CG"], Is.EqualTo(2.0 / 7.0).Within(Tolerance),
                "f_CG = 2/7");
            Assert.That(freq["AT"], Is.EqualTo(1.0 / 7.0).Within(Tolerance),
                "f_AT = 1/7");
            Assert.That(freq["TG"], Is.EqualTo(1.0 / 7.0).Within(Tolerance),
                "f_TG = 1/7");
            Assert.That(freq["GT"], Is.EqualTo(1.0 / 7.0).Within(Tolerance),
                "f_GT = 1/7");
        });
    }

    // M3 — INV-01: dinucleotide frequencies sum to 1.0.
    [Test]
    public void CalculateDinucleotideFrequencies_AtgcgcgtSequence_SumsToOne()
    {
        var freq = SequenceStatistics.CalculateDinucleotideFrequencies("ATGCGCGT");

        Assert.That(freq.Values.Sum(), Is.EqualTo(1.0).Within(Tolerance),
            "normalized frequencies over all valid dinucleotide positions sum to 1.0 (INV-01)");
    }

    // S2 — guards: null, empty, length < 2 return empty.
    [Test]
    public void CalculateDinucleotideFrequencies_NullEmptyOrTooShort_ReturnsEmpty()
    {
        Assert.Multiple(() =>
        {
            Assert.That(SequenceStatistics.CalculateDinucleotideFrequencies(null!), Is.Empty,
                "null input yields empty dictionary");
            Assert.That(SequenceStatistics.CalculateDinucleotideFrequencies(string.Empty), Is.Empty,
                "empty input yields empty dictionary");
            Assert.That(SequenceStatistics.CalculateDinucleotideFrequencies("A"), Is.Empty,
                "length < 2 yields empty dictionary");
        });
    }

    // C2 — RNA U handling: U is a valid base; AUGCGC includes the AU dinucleotide.
    [Test]
    public void CalculateDinucleotideFrequencies_RnaInput_IncludesUracilDinucleotides()
    {
        // AUGCGC -> dinucs AU,UG,GC,CG,GC (5). AU=1/5, GC=2/5.
        var freq = SequenceStatistics.CalculateDinucleotideFrequencies("AUGCGC");

        Assert.Multiple(() =>
        {
            Assert.That(freq.ContainsKey("AU"), Is.True,
                "U is in the dinucleotide alphabet, so AU is counted (RNA support)");
            Assert.That(freq["AU"], Is.EqualTo(1.0 / 5.0).Within(Tolerance),
                "f_AU = 1/5 (1 of 5 dinucleotide positions)");
            Assert.That(freq["GC"], Is.EqualTo(2.0 / 5.0).Within(Tolerance),
                "f_GC = 2/5");
        });
    }

    #endregion

    #region CalculateCodonFrequencies

    // M4 — frame 0 of ATGATGAAA: codons ATG, ATG, AAA -> ATG=2/3, AAA=1/3.
    // Evidence: Kazusa CUTG (count/total over non-overlapping triplets).
    [Test]
    public void CalculateCodonFrequencies_Frame0_ReturnsExactFrequencies()
    {
        var freq = SequenceStatistics.CalculateCodonFrequencies("ATGATGAAA", readingFrame: 0);

        Assert.Multiple(() =>
        {
            Assert.That(freq["ATG"], Is.EqualTo(2.0 / 3.0).Within(Tolerance),
                "ATG occurs 2 of 3 codons in frame 0");
            Assert.That(freq["AAA"], Is.EqualTo(1.0 / 3.0).Within(Tolerance),
                "AAA occurs 1 of 3 codons in frame 0");
            Assert.That(freq.Keys, Is.EquivalentTo(new[] { "ATG", "AAA" }),
                "only ATG and AAA are counted in frame 0");
        });
    }

    // M5 — frame 1 of ATGATGAAA: from index 1 -> TGA, TGA (then 'AA' leftover) -> TGA=1.0.
    // Evidence: Kazusa CUTG (frame offset, non-overlapping, trailing bases ignored).
    [Test]
    public void CalculateCodonFrequencies_Frame1_ReturnsExactFrequencies()
    {
        var freq = SequenceStatistics.CalculateCodonFrequencies("ATGATGAAA", readingFrame: 1);

        Assert.Multiple(() =>
        {
            Assert.That(freq["TGA"], Is.EqualTo(1.0).Within(Tolerance),
                "TGA is the only codon in frame 1 (both triplets), so frequency 1.0");
            Assert.That(freq.Keys, Is.EquivalentTo(new[] { "TGA" }),
                "trailing 'AA' is ignored; only TGA is counted");
        });
    }

    // M6 — INV-02: codon frequencies sum to 1.0.
    [Test]
    public void CalculateCodonFrequencies_Frame0_SumsToOne()
    {
        var freq = SequenceStatistics.CalculateCodonFrequencies("ATGATGAAA", readingFrame: 0);

        Assert.That(freq.Values.Sum(), Is.EqualTo(1.0).Within(Tolerance),
            "codon frequencies over all counted codons sum to 1.0 (INV-02)");
    }

    // S5 — non-ACGT triplets are excluded: ATGNNNAAA frame 0 -> ATG, NNN(excluded), AAA -> each 0.5.
    // Evidence: Kazusa CUTG (ambiguous codons excluded).
    [Test]
    public void CalculateCodonFrequencies_NonAcgtTriplet_ExcludedFromCounts()
    {
        var freq = SequenceStatistics.CalculateCodonFrequencies("ATGNNNAAA", readingFrame: 0);

        Assert.Multiple(() =>
        {
            Assert.That(freq.ContainsKey("NNN"), Is.False,
                "NNN contains non-ACGT bases and is excluded");
            Assert.That(freq["ATG"], Is.EqualTo(0.5).Within(Tolerance),
                "ATG is 1 of 2 valid codons -> 0.5");
            Assert.That(freq["AAA"], Is.EqualTo(0.5).Within(Tolerance),
                "AAA is 1 of 2 valid codons -> 0.5");
        });
    }

    // S6 — trailing 1-2 bases are ignored: ATGAA frame 0 -> only ATG; ATG=1.0.
    // Evidence: Kazusa CUTG (non-overlapping triplets; remainder ignored).
    [Test]
    public void CalculateCodonFrequencies_TrailingBases_Ignored()
    {
        var freq = SequenceStatistics.CalculateCodonFrequencies("ATGAA", readingFrame: 0);

        Assert.Multiple(() =>
        {
            Assert.That(freq["ATG"], Is.EqualTo(1.0).Within(Tolerance),
                "only the full ATG triplet is counted; trailing 'AA' is ignored");
            Assert.That(freq.Keys, Is.EquivalentTo(new[] { "ATG" }),
                "no partial codon is added");
        });
    }

    // S3 — guards: null, empty, length < 3 return empty.
    [Test]
    public void CalculateCodonFrequencies_NullEmptyOrTooShort_ReturnsEmpty()
    {
        Assert.Multiple(() =>
        {
            Assert.That(SequenceStatistics.CalculateCodonFrequencies(null!), Is.Empty,
                "null input yields empty dictionary");
            Assert.That(SequenceStatistics.CalculateCodonFrequencies(string.Empty), Is.Empty,
                "empty input yields empty dictionary");
            Assert.That(SequenceStatistics.CalculateCodonFrequencies("AT"), Is.Empty,
                "length < 3 yields empty dictionary (no full codon)");
        });
    }

    #endregion
}
