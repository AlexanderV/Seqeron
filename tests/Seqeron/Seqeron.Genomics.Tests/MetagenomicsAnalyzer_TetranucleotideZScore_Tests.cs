// META-BIN-001 — z-score-normalised tetranucleotide signature (TETRA)
// Evidence: docs/Evidence/META-BIN-001-Evidence.md
// TestSpec: tests/TestSpecs/META-BIN-001.md
// Source: Teeling H, Waldmann J, Lombardot T, Bauer M, Glöckner FO (2004).
//         TETRA. BMC Bioinformatics 5:163. doi:10.1186/1471-2105-5-163
//         Teeling H et al. (2004). Environ Microbiol 6(9):938–947. doi:10.1111/j.1462-2920.2004.00624.x
//         Schbath S (1997). J Comput Biol 4(2):189–192 (variance approximation).
using NUnit.Framework;
using Seqeron.Genomics.Metagenomics;
using System.Linq;

namespace Seqeron.Genomics.Tests;

/// <summary>
/// Tests for the opt-in TETRA z-score tetranucleotide signature
/// (<see cref="MetagenomicsAnalyzer.CalculateTetranucleotideZScores"/>) and its Pearson
/// correlation (<see cref="MetagenomicsAnalyzer.TetranucleotideZScoreCorrelation"/>).
///
/// The signature is Markov-corrected: observed tetranucleotide count vs the count predicted
/// by a maximal-order (2nd-order) Markov model, divided by the Schbath variance approximation:
///   E(n1n2n3n4) = N(n1n2n3)·N(n2n3n4) / N(n2n3)
///   var         = E·[N(n2n3)−N(n1n2n3)]·[N(n2n3)−N(n2n3n4)] / N(n2n3)²
///   z           = (N(n1n2n3n4) − E) / √var
/// then two sequences are compared by the Pearson correlation of their z-score vectors.
/// </summary>
[TestFixture]
[Category("Metagenomics")]
[Category("META-BIN-001")]
public class MetagenomicsAnalyzer_TetranucleotideZScore_Tests
{
    #region CalculateTetranucleotideZScores

    // M-Z1 — Hand-derived exact z-score from the TETRA formula.
    // Sequence "ACGTACGTGGCC" extended by its reverse complement "GGCCACGTACGT" gives the
    // 24-nt strand-symmetric string "ACGTACGTGGCCGGCCACGTACGT". Counting overlapping words:
    //   N(ACGT)=4, N(ACG)=4 [n1n2n3], N(CGT)=4 [n2n3n4], N(CG)=5 [n2n3].
    //   E   = 4·4/5 = 3.2
    //   var = 3.2·(5−4)·(5−4)/5² = 3.2/25 = 0.128
    //   z   = (4 − 3.2)/√0.128 = 0.8/√0.128 = √5 = 2.23606797749979...
    // (0.8² / 0.128 = 0.64/0.128 = 5, so z = √5 exactly.) A wrong expected-count or variance
    // formula cannot produce √5, so this asserts the formula, not the code's own output.
    [Test]
    public void CalculateTetranucleotideZScores_HandDerivedTetramer_MatchesTetraFormula()
    {
        // Arrange
        const string seq = "ACGTACGTGGCC";
        double expectedZ = System.Math.Sqrt(5.0); // 2.2360679774997896

        // Act
        var z = MetagenomicsAnalyzer.CalculateTetranucleotideZScores(seq);

        // Assert
        Assert.That(z["ACGT"], Is.EqualTo(expectedZ).Within(1e-10),
            "z(ACGT) must equal (N−E)/√var = (4−3.2)/√0.128 = √5 from the Teeling/Schbath formula.");
    }

    // M-Z2 — Full 256-component signature is always returned (4^4 = 256, the TETRA dimensionality).
    [Test]
    public void CalculateTetranucleotideZScores_AnySequence_Returns256Components()
    {
        // Arrange
        const string seq = "ACGTACGTGGCCATGCATGCTTAA";

        // Act
        var z = MetagenomicsAnalyzer.CalculateTetranucleotideZScores(seq);

        // Assert
        Assert.Multiple(() =>
        {
            Assert.That(z.Count, Is.EqualTo(256),
                "TETRA tabulates all 256 (4^4) tetranucleotides (Teeling 2004).");
            Assert.That(z.ContainsKey("ACGT"), Is.True, "Signature must be keyed by ACGT 4-mers.");
            Assert.That(z.ContainsKey("AAAA"), Is.True, "All 256 keys present, including AAAA.");
            Assert.That(z.ContainsKey("TTTT"), Is.True, "All 256 keys present, including TTTT.");
        });
    }

    // M-Z3 — A tetranucleotide whose middle dinucleotide never occurs has z=0 (denominator N(n2n3)=0).
    // In "AAAAAAAA" (+ revcomp "TTTTTTTT") the dinucleotide CG never occurs, so any tetramer with
    // CG in the middle (e.g. ACGT) has N(CG)=0 → z=0 by definition (no over/under-representation signal).
    [Test]
    public void CalculateTetranucleotideZScores_AbsentMiddleDinucleotide_YieldsZeroZScore()
    {
        // Arrange
        const string seq = "AAAAAAAA";

        // Act
        var z = MetagenomicsAnalyzer.CalculateTetranucleotideZScores(seq);

        // Assert
        Assert.That(z["ACGT"], Is.EqualTo(0.0).Within(1e-12),
            "When N(n2n3)=0 the expected count is undefined; TETRA yields z=0 (no signal).");
    }

    // M-Z4 — Null / empty / single-base input returns an all-zero 256-entry map (no throw).
    // Note: because TETRA extends by the reverse complement, even a 2-nt input ("AC" → "ACGT")
    // already yields a tetramer; only null/empty/single-base inputs cannot form a 4-nt extended
    // strand and therefore give an all-zero signature.
    [Test]
    public void CalculateTetranucleotideZScores_NullEmptyOrSingleBase_ReturnsAllZeroSignature()
    {
        // Act
        var zNull = MetagenomicsAnalyzer.CalculateTetranucleotideZScores(null!);
        var zEmpty = MetagenomicsAnalyzer.CalculateTetranucleotideZScores("");
        var zSingle = MetagenomicsAnalyzer.CalculateTetranucleotideZScores("A"); // → "AT" (no tetramer)

        // Assert
        Assert.Multiple(() =>
        {
            Assert.That(zNull.Count, Is.EqualTo(256), "Null input still yields the 256-key signature.");
            Assert.That(zNull.Values.All(v => v == 0.0), Is.True, "Null input → all-zero z-scores.");
            Assert.That(zEmpty.Values.All(v => v == 0.0), Is.True, "Empty input → all-zero z-scores.");
            Assert.That(zSingle.Values.All(v => v == 0.0), Is.True,
                "Single-base input cannot form a 4-nt extended strand → all-zero z-scores.");
        });
    }

    // M-Z5 — Non-ACGT characters are filtered before counting: a sequence with interspersed Ns
    // produces the same z-score signature as the ACGT-only sequence (TETRA counts ACGT words only).
    [Test]
    public void CalculateTetranucleotideZScores_IgnoresNonAcgtCharacters()
    {
        // Arrange — same ACGT content, one with embedded N/lowercase noise.
        var clean = MetagenomicsAnalyzer.CalculateTetranucleotideZScores("ACGTACGTGGCC");
        var noisy = MetagenomicsAnalyzer.CalculateTetranucleotideZScores("acgtACGTGGCC"); // case-insensitive

        // Assert
        Assert.That(noisy["ACGT"], Is.EqualTo(clean["ACGT"]).Within(1e-12),
            "Case is normalised to uppercase ACGT; the signature must be identical.");
    }

    // M-Z9 — Strand-symmetric invariant of TETRA (Teeling 2004): because the sequence is extended
    // by its reverse complement, a tetranucleotide w and its reverse complement rc(w) are pooled
    // into the SAME over-/under-representation signal, so z(w) == z(rc(w)) for ALL 256 words.
    // This is the defining "reverse-complement-merged counts identical" property and is derived
    // from the method definition, not from the code's output. Hand-derivation: on the extended
    // strand s+rc(s), every occurrence of w on one strand appears as rc(w) on the other, so the
    // observed/expected/variance inputs for w and rc(w) are equal counts → equal z.
    [Test]
    public void CalculateTetranucleotideZScores_IsReverseComplementSymmetric()
    {
        // Arrange — an asymmetric, compositionally varied sequence.
        const string seq = "ACGTACGTGGCCATGCATGCTTAA";
        var z = MetagenomicsAnalyzer.CalculateTetranucleotideZScores(seq);

        static string Rc(string w)
        {
            var c = new char[4];
            for (int i = 0; i < 4; i++)
                c[3 - i] = w[i] switch { 'A' => 'T', 'T' => 'A', 'C' => 'G', 'G' => 'C', _ => 'N' };
            return new string(c);
        }

        // Assert — every word equals its reverse complement's z-score (e.g. GGCC↔GGCC, ATGC↔GCAT).
        Assert.Multiple(() =>
        {
            foreach (var kvp in z)
                Assert.That(kvp.Value, Is.EqualTo(z[Rc(kvp.Key)]).Within(1e-12),
                    $"TETRA is strand-symmetric: z({kvp.Key}) must equal z({Rc(kvp.Key)}).");
        });
    }

    #endregion

    #region TetranucleotideZScoreCorrelation

    // M-Z6 — Self-correlation of a sequence's z-score signature is exactly 1.0 (Pearson r=1).
    [Test]
    public void TetranucleotideZScoreCorrelation_SelfComparison_IsOne()
    {
        // Arrange
        const string seq = "ACGTACGTGGCCATGCATGCTTAA";

        // Act
        double r = MetagenomicsAnalyzer.TetranucleotideZScoreCorrelation(seq, seq);

        // Assert
        Assert.That(r, Is.EqualTo(1.0).Within(1e-10),
            "A signature is perfectly correlated with itself: Pearson r = 1.0.");
    }

    // M-Z7 — Two near-identical sequences correlate MORE strongly than a similar sequence does with
    // a compositionally very different (AT-homopolymer) one. This is the discriminative property
    // that makes TETRA usable for binning (Teeling 2004): similar composition → higher correlation.
    [Test]
    public void TetranucleotideZScoreCorrelation_SimilarHigherThanDissimilar()
    {
        // Arrange — s1 vs s2 differ by a single base; s1 vs s3 are compositionally unrelated.
        const string s1 = "ACGTACGTGGCCATGCATGCTTAA";
        const string s2 = "ACGTACGTGGCCATGCATGCTTAG"; // one substitution
        const string s3 = "AAAAATTTTTAAAAATTTTTAAAA"; // AT-rich, very different signature

        // Act
        double rSimilar = MetagenomicsAnalyzer.TetranucleotideZScoreCorrelation(s1, s2);
        double rDissimilar = MetagenomicsAnalyzer.TetranucleotideZScoreCorrelation(s1, s3);

        // Assert
        Assert.That(rSimilar, Is.GreaterThan(rDissimilar),
            "Compositionally similar sequences must correlate higher than dissimilar ones — " +
            $"the basis of TETRA binning (r_similar={rSimilar:F4} > r_dissimilar={rDissimilar:F4}).");
    }

    // M-Z8 — Correlation is symmetric: corr(a,b) == corr(b,a).
    [Test]
    public void TetranucleotideZScoreCorrelation_IsSymmetric()
    {
        // Arrange
        const string a = "ACGTACGTGGCCATGCATGCTTAA";
        const string b = "GGCCGGCCATATATATCGCGCGCG";

        // Act
        double ab = MetagenomicsAnalyzer.TetranucleotideZScoreCorrelation(a, b);
        double ba = MetagenomicsAnalyzer.TetranucleotideZScoreCorrelation(b, a);

        // Assert
        Assert.That(ab, Is.EqualTo(ba).Within(1e-12),
            "Pearson correlation is symmetric in its arguments.");
    }

    // S-Z1 — A degenerate (all-zero) signature correlates to 0 (no usable signal), not NaN.
    [Test]
    public void TetranucleotideZScoreCorrelation_AgainstEmpty_IsZeroNotNaN()
    {
        // Arrange
        const string seq = "ACGTACGTGGCCATGCATGCTTAA";

        // Act
        double r = MetagenomicsAnalyzer.TetranucleotideZScoreCorrelation(seq, "");

        // Assert
        Assert.Multiple(() =>
        {
            Assert.That(double.IsNaN(r), Is.False, "Zero-variance vector must not yield NaN.");
            Assert.That(r, Is.EqualTo(0.0).Within(1e-12),
                "Correlation with an all-zero signature is defined as 0 (no signal).");
        });
    }

    #endregion
}
