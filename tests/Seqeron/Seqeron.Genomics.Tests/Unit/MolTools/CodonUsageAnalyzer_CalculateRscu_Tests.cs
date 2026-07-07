// CODON-RSCU-001 — Relative Synonymous Codon Usage (RSCU) and codon counting
// Evidence: docs/Evidence/CODON-RSCU-001-Evidence.md
// TestSpec: tests/TestSpecs/CODON-RSCU-001.md
// Source: Sharp PM, Tuohy TMF, Mosurski KR (1986). Nucleic Acids Res. 14(13):5125-5143.
//         Formula RSCU = (n_i * x_ij) / sum_k x_ik per GenomicSig (CRAN) / LIRMM.

namespace Seqeron.Genomics.Tests.Unit.MolTools;

[TestFixture]
public class CodonUsageAnalyzer_CalculateRscu_Tests
{
    // Leu (6-fold) synonymous family per the standard genetic code.
    private static readonly string[] LeuCodons = { "TTA", "TTG", "CTT", "CTC", "CTA", "CTG" };

    #region CalculateRscu

    // M1 — Leu 6-fold family: CTG x3, CTA x1; n_i=6, total=4.
    // RSCU(CTG)=6*3/4=4.5, RSCU(CTA)=6*1/4=1.5, others=0. [GenomicSig formula; Evidence dataset]
    [Test]
    public void CalculateRscu_SixfoldLeuFamily_ComputesExactFormulaValues()
    {
        var rscu = CodonUsageAnalyzer.CalculateRscu(new DnaSequence("CTGCTGCTGCTA"));

        Assert.Multiple(() =>
        {
            Assert.That(rscu["CTG"], Is.EqualTo(4.5).Within(1e-10),
                "RSCU(CTG)=6*3/4=4.5 for n_i=6, x=3, total=4");
            Assert.That(rscu["CTA"], Is.EqualTo(1.5).Within(1e-10),
                "RSCU(CTA)=6*1/4=1.5 for n_i=6, x=1, total=4");
            Assert.That(rscu["TTA"], Is.EqualTo(0.0).Within(1e-10), "TTA unused in Leu family -> RSCU 0");
            Assert.That(rscu["TTG"], Is.EqualTo(0.0).Within(1e-10), "TTG unused in Leu family -> RSCU 0");
            Assert.That(rscu["CTT"], Is.EqualTo(0.0).Within(1e-10), "CTT unused in Leu family -> RSCU 0");
            Assert.That(rscu["CTC"], Is.EqualTo(0.0).Within(1e-10), "CTC unused in Leu family -> RSCU 0");
        });
    }

    // M2 — Phe 2-fold biased: TTT x2, TTC x1; n_i=2, total=3.
    // RSCU(TTT)=2*2/3=4/3, RSCU(TTC)=2*1/3=2/3. [GenomicSig formula]
    [Test]
    public void CalculateRscu_TwofoldBiasedFamily_ComputesExactFractions()
    {
        var rscu = CodonUsageAnalyzer.CalculateRscu(new DnaSequence("TTTTTTTTC"));

        Assert.Multiple(() =>
        {
            Assert.That(rscu["TTT"], Is.EqualTo(4.0 / 3.0).Within(1e-10),
                "RSCU(TTT)=2*2/3=4/3 for n_i=2, x=2, total=3");
            Assert.That(rscu["TTC"], Is.EqualTo(2.0 / 3.0).Within(1e-10),
                "RSCU(TTC)=2*1/3=2/3 for n_i=2, x=1, total=3");
        });
    }

    // M3 — Unbiased: equal usage within Phe family -> RSCU 1 (no bias). [seqinr/LIRMM no-bias value 1]
    [Test]
    public void CalculateRscu_EqualUsageWithinFamily_ReturnsOne()
    {
        var rscu = CodonUsageAnalyzer.CalculateRscu(new DnaSequence("TTTTTC"));

        Assert.Multiple(() =>
        {
            Assert.That(rscu["TTT"], Is.EqualTo(1.0).Within(1e-10),
                "Equal usage of synonymous codons -> RSCU 1 (no bias)");
            Assert.That(rscu["TTC"], Is.EqualTo(1.0).Within(1e-10),
                "Equal usage of synonymous codons -> RSCU 1 (no bias)");
        });
    }

    // M4 — Single-codon family (Met=ATG, n_i=1) -> RSCU 1 when present. [GenomicSig bounds; cubar]
    [Test]
    public void CalculateRscu_SingleCodonFamily_ReturnsOne()
    {
        var rscu = CodonUsageAnalyzer.CalculateRscu(new DnaSequence("ATGATG"));

        Assert.That(rscu["ATG"], Is.EqualTo(1.0).Within(1e-10),
            "Met has one codon (n_i=1) -> RSCU=1 whenever present");
    }

    // M5 — INV-04: RSCU values within a present family sum to n_i (=6 for Leu). [derivation from formula]
    [Test]
    public void CalculateRscu_PresentFamily_SumsToDegeneracy()
    {
        var rscu = CodonUsageAnalyzer.CalculateRscu(new DnaSequence("CTGCTGCTGCTA"));

        double familySum = LeuCodons.Sum(c => rscu[c]);

        Assert.That(familySum, Is.EqualTo(6.0).Within(1e-10),
            "Sum of RSCU over a present family equals its degeneracy n_i (6 for Leu)");
    }

    // C1 — INV-03: every RSCU value is within [0, n_i] (0..6 for Leu). [GenomicSig bounds]
    [Test]
    public void CalculateRscu_PresentFamily_ValuesWithinBounds()
    {
        var rscu = CodonUsageAnalyzer.CalculateRscu(new DnaSequence("CTGCTGCTGCTA"));

        Assert.Multiple(() =>
        {
            foreach (var c in LeuCodons)
            {
                Assert.That(rscu[c], Is.InRange(0.0, 6.0),
                    $"RSCU({c}) must lie in [0, n_i]=[0,6] for the Leu family");
            }
        });
    }

    // C2 — Absent synonymous family (0/0): the repository convention is to return 0 for
    // every codon of a family that never appears (no pseudocount), since canonical RSCU is
    // undefined for 0/0. Input `ATGATG` contains only Met, so the entire Leu family is absent.
    // [TestSpec Assumption #1 / Evidence corner case "Absent family (0/0) → repository returns 0";
    //  cubar est_rscu documents this as implementation-defined.]
    [Test]
    public void CalculateRscu_AbsentFamily_ReturnsZeroForEveryCodon()
    {
        var rscu = CodonUsageAnalyzer.CalculateRscu(new DnaSequence("ATGATG"));

        Assert.Multiple(() =>
        {
            foreach (var c in LeuCodons)
            {
                Assert.That(rscu[c], Is.EqualTo(0.0).Within(1e-10),
                    $"Leu family is absent in 'ATGATG' (0/0) -> repository returns RSCU({c})=0");
            }
            Assert.That(rscu["TTT"], Is.EqualTo(0.0).Within(1e-10),
                "Phe family is absent (0/0) -> RSCU(TTT)=0");
        });
    }

    // C3 — Stop codons treated as one 3-fold synonymous family (TAA/TAG/TGA), degeneracy 3.
    // Input `TAATAGTGA` uses each stop once: total=3, RSCU = 3*1/3 = 1.0 for each.
    // [TestSpec Assumption #2 / Evidence corner case: repository groups the three standard
    //  stop codons as a synonymous family of size 3 and computes RSCU like any family. The
    //  family-ratio formula RSCU=(n_i*x)/Σx (Source 3 / LIRMM) gives 1.0 for equal usage.]
    [Test]
    public void CalculateRscu_StopCodonFamily_TreatedAsThreeFoldFamily()
    {
        var rscu = CodonUsageAnalyzer.CalculateRscu(new DnaSequence("TAATAGTGA"));

        Assert.Multiple(() =>
        {
            Assert.That(rscu["TAA"], Is.EqualTo(1.0).Within(1e-10),
                "Stop family size 3, each used once -> RSCU(TAA)=3*1/3=1.0");
            Assert.That(rscu["TAG"], Is.EqualTo(1.0).Within(1e-10),
                "Stop family size 3, each used once -> RSCU(TAG)=3*1/3=1.0");
            Assert.That(rscu["TGA"], Is.EqualTo(1.0).Within(1e-10),
                "Stop family size 3, each used once -> RSCU(TGA)=3*1/3=1.0");
        });
    }

    // C4 — Stop family biased: `TAATAATGA` (TAA x2, TGA x1, TAG x0); total=3, n_i=3.
    // RSCU(TAA)=3*2/3=2.0, RSCU(TGA)=3*1/3=1.0, RSCU(TAG)=3*0/3=0.0; sum over family = 3 = n_i.
    // [Source 3 / LIRMM family-ratio formula applied to the 3-fold stop family.]
    [Test]
    public void CalculateRscu_StopFamilyBiased_ComputesFamilyRatioAndSumsToDegeneracy()
    {
        var rscu = CodonUsageAnalyzer.CalculateRscu(new DnaSequence("TAATAATGA"));

        Assert.Multiple(() =>
        {
            Assert.That(rscu["TAA"], Is.EqualTo(2.0).Within(1e-10),
                "RSCU(TAA)=3*2/3=2.0 for n_i=3, x=2, total=3");
            Assert.That(rscu["TGA"], Is.EqualTo(1.0).Within(1e-10),
                "RSCU(TGA)=3*1/3=1.0 for n_i=3, x=1, total=3");
            Assert.That(rscu["TAG"], Is.EqualTo(0.0).Within(1e-10),
                "RSCU(TAG)=3*0/3=0.0 (unused in present stop family)");
            Assert.That(rscu["TAA"] + rscu["TAG"] + rscu["TGA"], Is.EqualTo(3.0).Within(1e-10),
                "Sum of RSCU over the present stop family equals its degeneracy n_i=3");
        });
    }

    // S1 — null DnaSequence throws.
    [Test]
    public void CalculateRscu_NullDnaSequence_ThrowsArgumentNullException()
    {
        Assert.Throws<System.ArgumentNullException>(
            () => CodonUsageAnalyzer.CalculateRscu((DnaSequence)null!),
            "Null DnaSequence is an invalid input and must throw ArgumentNullException");
    }

    // S2 — empty string returns empty dictionary.
    [Test]
    public void CalculateRscu_EmptyString_ReturnsEmpty()
    {
        var rscu = CodonUsageAnalyzer.CalculateRscu("");

        Assert.That(rscu, Is.Empty, "Empty input yields no codons and therefore an empty RSCU map");
    }

    // S5 — string overload delegates to the same computation (exact M1 values).
    [Test]
    public void CalculateRscu_StringOverload_ProducesSameValuesAsDnaSequence()
    {
        var rscu = CodonUsageAnalyzer.CalculateRscu("CTGCTGCTGCTA");

        Assert.Multiple(() =>
        {
            Assert.That(rscu["CTG"], Is.EqualTo(4.5).Within(1e-10), "String overload delegates: RSCU(CTG)=4.5");
            Assert.That(rscu["CTA"], Is.EqualTo(1.5).Within(1e-10), "String overload delegates: RSCU(CTA)=1.5");
        });
    }

    #endregion

    #region CountCodons

    // M6 — non-overlapping triplet counting from offset 0. [repository contract / Kazusa convention]
    [Test]
    public void CountCodons_CodingSequence_CountsNonOverlappingTriplets()
    {
        var counts = CodonUsageAnalyzer.CountCodons(new DnaSequence("ATGAAATGA"));

        Assert.Multiple(() =>
        {
            Assert.That(counts["ATG"], Is.EqualTo(1), "ATG occurs once as a frame-0 triplet");
            Assert.That(counts["AAA"], Is.EqualTo(1), "AAA occurs once as a frame-0 triplet");
            Assert.That(counts["TGA"], Is.EqualTo(1), "TGA occurs once as a frame-0 triplet");
            Assert.That(counts.Values.Sum(), Is.EqualTo(3), "Three full triplets counted");
        });
    }

    // M7 — repeated codon counted with full multiplicity.
    [Test]
    public void CountCodons_RepeatedCodon_CountsMultiplicity()
    {
        var counts = CodonUsageAnalyzer.CountCodons(new DnaSequence("ATGATGATG"));

        Assert.That(counts["ATG"], Is.EqualTo(3), "Three consecutive ATG triplets -> count 3");
    }

    // M8 — trailing 1-2 bases are not a codon and are ignored.
    [Test]
    public void CountCodons_TrailingPartialCodon_IsIgnored()
    {
        var counts = CodonUsageAnalyzer.CountCodons(new DnaSequence("ATGAA"));

        Assert.Multiple(() =>
        {
            Assert.That(counts["ATG"], Is.EqualTo(1), "Only the full ATG triplet is counted");
            Assert.That(counts.Values.Sum(), Is.EqualTo(1), "Trailing 'AA' is not a triplet and is ignored");
        });
    }

    // M9 — triplets with non-ACGT characters are excluded (string overload).
    [Test]
    public void CountCodons_NonAcgtTriplet_IsExcluded()
    {
        var counts = CodonUsageAnalyzer.CountCodons("ATGNNNAAA");

        Assert.Multiple(() =>
        {
            Assert.That(counts["ATG"], Is.EqualTo(1), "Valid ATG triplet counted");
            Assert.That(counts["AAA"], Is.EqualTo(1), "Valid AAA triplet counted");
            Assert.That(counts.ContainsKey("NNN"), Is.False, "Triplet with non-ACGT base is excluded");
            Assert.That(counts.Values.Sum(), Is.EqualTo(2), "Only the two ACGT triplets are counted");
        });
    }

    // S3 — null DnaSequence throws.
    [Test]
    public void CountCodons_NullDnaSequence_ThrowsArgumentNullException()
    {
        Assert.Throws<System.ArgumentNullException>(
            () => CodonUsageAnalyzer.CountCodons((DnaSequence)null!),
            "Null DnaSequence is an invalid input and must throw ArgumentNullException");
    }

    // S4 — empty string returns empty dictionary.
    [Test]
    public void CountCodons_EmptyString_ReturnsEmpty()
    {
        var counts = CodonUsageAnalyzer.CountCodons("");

        Assert.That(counts, Is.Empty, "Empty input yields no triplets and therefore an empty count map");
    }

    // S6 — string overload delegates and is case-insensitive (uppercases first).
    [Test]
    public void CountCodons_StringOverloadLowercase_IsCaseInsensitive()
    {
        var counts = CodonUsageAnalyzer.CountCodons("atgaaatga");

        Assert.Multiple(() =>
        {
            Assert.That(counts["ATG"], Is.EqualTo(1), "Lowercase input is uppercased before counting");
            Assert.That(counts["AAA"], Is.EqualTo(1), "Lowercase input is uppercased before counting");
            Assert.That(counts["TGA"], Is.EqualTo(1), "Lowercase input is uppercased before counting");
        });
    }

    #endregion
}
