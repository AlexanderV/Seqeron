// SEQ-PI-001 — Isoelectric Point (pI) Calculation
// Evidence: docs/Evidence/SEQ-PI-001-Evidence.md
// TestSpec: tests/TestSpecs/SEQ-PI-001.md
// Source: EMBOSS iep (Epk.dat pKa scale), https://emboss.sourceforge.net/emboss/apps/iep.html;
//         Peptides charge_pI.cpp (Osorio et al. 2015) Henderson-Hasselbalch net-charge model.

using Seqeron.Genomics.Analysis;

namespace Seqeron.Genomics.Tests;

[TestFixture]
public class SequenceStatistics_CalculateIsoelectricPoint_Tests
{
    // pI is returned rounded to 2 decimals (bisection precision 0.01); expected values are
    // derived from the EMBOSS pKa scale, whose charge function was confirmed against the
    // Peptides EMBOSS worked example (charge 3.037398/2.914112/0.7184524 at pH 5/7/9).
    private const double Tolerance = 0.01;

    #region CalculateIsoelectricPoint

    // M1 — basic reference peptide "FLPVLAGLTPSIVPKLVCLLTKKC".
    // Evidence: Peptides EMBOSS net charge is +0.7184524 at pH 9, so the zero-charge pH lies
    // above pH 9; bisection on the EMBOSS scale yields pI = 9.67.
    [Test]
    public void CalculateIsoelectricPoint_BasicReferencePeptide_Returns967()
    {
        double pi = SequenceStatistics.CalculateIsoelectricPoint("FLPVLAGLTPSIVPKLVCLLTKKC");

        Assert.That(pi, Is.EqualTo(9.67).Within(Tolerance),
            "EMBOSS-scale net charge is still +0.72 at pH 9, so pI is basic (9.67); validates the charge formula + pKa set");
    }

    // M2 / INV-01 — pI bounds 0..14 for any input, anchored to exact sourced values.
    // Evidence: bisection is confined to [0, 14] (EMBOSS iep). The exact pI values
    // (DDDDDDDD = 2.96, RRRRRRRR = 13.35) were reproduced by an independent reference
    // implementation built from the EMBOSS Epk.dat pKa scale + Henderson-Hasselbalch
    // charge formula (the same formula that reproduces the Peptides worked example).
    // Asserting the exact values — not just the bounds — means this test fails against a
    // deliberately-wrong implementation that merely returns something inside [0,14].
    [Test]
    public void CalculateIsoelectricPoint_AnySequence_StaysWithinPhBoundsWithExactValues()
    {
        double acidic = SequenceStatistics.CalculateIsoelectricPoint("DDDDDDDD");
        double basic = SequenceStatistics.CalculateIsoelectricPoint("RRRRRRRR");

        Assert.Multiple(() =>
        {
            Assert.That(acidic, Is.EqualTo(2.96).Within(Tolerance),
                "EMBOSS-scale pI of an 8-Asp peptide is 2.96 (acidic-dominated, within [0,14])");
            Assert.That(basic, Is.EqualTo(13.35).Within(Tolerance),
                "EMBOSS-scale pI of an 8-Arg peptide is 13.35 (highly basic, within [0,14])");
            Assert.That(acidic, Is.GreaterThanOrEqualTo(0.0).And.LessThanOrEqualTo(14.0),
                "INV-01: pI must lie in [0,14] for an acidic peptide");
            Assert.That(basic, Is.GreaterThanOrEqualTo(0.0).And.LessThanOrEqualTo(14.0),
                "INV-01: pI must lie in [0,14] for a basic peptide");
        });
    }

    // M3 / INV-04 — termini-only "A": pI = midpoint of N-term (8.6) and C-term (3.6) pKa = 6.10.
    // Evidence: EMBOSS pKa; with no ionizable side chains the two terminal terms cancel at the midpoint.
    [Test]
    public void CalculateIsoelectricPoint_TerminiOnlyAlanine_ReturnsPkaMidpoint610()
    {
        double pi = SequenceStatistics.CalculateIsoelectricPoint("A");

        Assert.That(pi, Is.EqualTo(6.10).Within(Tolerance),
            "INV-04: no side chains, so pI = (N-term 8.6 + C-term 3.6)/2 = 6.10");
    }

    // M4 — acidic-only "DDDD": four Asp pull pI down to 3.23 (EMBOSS scale).
    // Evidence: derived from EMBOSS pKa (Asp 3.9, C-term 3.6).
    [Test]
    public void CalculateIsoelectricPoint_AcidicTetraAspartate_Returns323()
    {
        double pi = SequenceStatistics.CalculateIsoelectricPoint("DDDD");

        Assert.That(pi, Is.EqualTo(3.23).Within(Tolerance),
            "Four acidic Asp residues drive pI into the acidic range (3.23) on the EMBOSS scale");
    }

    // M5 — basic-only "KKKK": four Lys push pI up to 11.27 (EMBOSS scale).
    // Evidence: derived from EMBOSS pKa (Lys 10.8, N-term 8.6).
    [Test]
    public void CalculateIsoelectricPoint_BasicTetraLysine_Returns1127()
    {
        double pi = SequenceStatistics.CalculateIsoelectricPoint("KKKK");

        Assert.That(pi, Is.EqualTo(11.27).Within(Tolerance),
            "Four basic Lys residues drive pI into the basic range (11.27) on the EMBOSS scale");
    }

    // M6 — all-20 residues "ACDEFGHIKLMNPQRSTVWY": EMBOSS-scale pI = 7.36.
    // Evidence: derived from EMBOSS pKa (note: the Bjellqvist/seqinr value 6.78454 is a different scale).
    [Test]
    public void CalculateIsoelectricPoint_AllTwentyResidues_Returns736OnEmbossScale()
    {
        double pi = SequenceStatistics.CalculateIsoelectricPoint("ACDEFGHIKLMNPQRSTVWY");

        Assert.That(pi, Is.EqualTo(7.36).Within(Tolerance),
            "One of each residue gives pI = 7.36 on the EMBOSS scale (distinct from Bjellqvist 6.78454)");
    }

    // M7 — empty string → neutral sentinel 7.0 (pI undefined for zero-length protein).
    // Evidence: ASSUMPTION input-guard convention.
    [Test]
    public void CalculateIsoelectricPoint_EmptyString_ReturnsNeutralSeven()
    {
        double pi = SequenceStatistics.CalculateIsoelectricPoint("");

        Assert.That(pi, Is.EqualTo(7.0),
            "Empty input has no defined pI; the documented input-guard sentinel is 7.0");
    }

    // M8 — null → neutral sentinel 7.0.
    // Evidence: ASSUMPTION input-guard convention.
    [Test]
    public void CalculateIsoelectricPoint_Null_ReturnsNeutralSeven()
    {
        double pi = SequenceStatistics.CalculateIsoelectricPoint(null!);

        Assert.That(pi, Is.EqualTo(7.0),
            "Null input has no defined pI; the documented input-guard sentinel is 7.0");
    }

    // S1 — single Asp "D": pI = 3.75 (one acidic side chain + termini).
    // Evidence: derived from EMBOSS pKa.
    [Test]
    public void CalculateIsoelectricPoint_SingleAspartate_Returns375()
    {
        double pi = SequenceStatistics.CalculateIsoelectricPoint("D");

        Assert.That(pi, Is.EqualTo(3.75).Within(Tolerance),
            "One acidic Asp plus termini gives pI = 3.75 on the EMBOSS scale");
    }

    // S2 — single Lys "K": pI = 9.70 (one basic side chain + termini).
    // Evidence: derived from EMBOSS pKa.
    [Test]
    public void CalculateIsoelectricPoint_SingleLysine_Returns970()
    {
        double pi = SequenceStatistics.CalculateIsoelectricPoint("K");

        Assert.That(pi, Is.EqualTo(9.70).Within(Tolerance),
            "One basic Lys plus termini gives pI = 9.70 on the EMBOSS scale");
    }

    // S3 — case-insensitivity: lowercase input yields the same pI as uppercase.
    // Evidence: input is normalized via ToUpperInvariant.
    [Test]
    public void CalculateIsoelectricPoint_LowercaseInput_MatchesUppercase()
    {
        double lower = SequenceStatistics.CalculateIsoelectricPoint("dddd");
        double upper = SequenceStatistics.CalculateIsoelectricPoint("DDDD");

        Assert.That(lower, Is.EqualTo(upper).Within(Tolerance),
            "pI is case-insensitive: lowercase 'dddd' must equal uppercase 'DDDD'");
    }

    // C1 / INV-02 — order-independence: pI is composition-only, so permutations are equal.
    // Evidence: EMBOSS "no electrostatic interactions" — charge summed over counts, not positions.
    [Test]
    public void CalculateIsoelectricPoint_PermutedSequence_HasIdenticalPi()
    {
        double dk = SequenceStatistics.CalculateIsoelectricPoint("DKDK");
        double kd = SequenceStatistics.CalculateIsoelectricPoint("KDDK");

        Assert.That(dk, Is.EqualTo(kd).Within(Tolerance),
            "INV-02: pI depends only on composition, so reordering the same residues gives the same pI");
    }

    // Coverage — non-ionizable characters (gaps, whitespace, punctuation, non-standard
    // residues) are ignored, never throw, and leave the result equal to the termini-only pI.
    // Evidence: doc §3.3 "non-standard residues, gaps, or whitespace do not throw"; the
    // composition model counts only the nine ionizable groups, so a string with no ionizable
    // side chains yields the termini-only midpoint 6.10 (verified by the independent reference).
    [Test]
    public void CalculateIsoelectricPoint_NonIonizableCharacters_IgnoredEqualsTerminiOnly()
    {
        double withNoise = SequenceStatistics.CalculateIsoelectricPoint("A B!G");
        double nonStandard = SequenceStatistics.CalculateIsoelectricPoint("XZ");

        Assert.Multiple(() =>
        {
            Assert.That(withNoise, Is.EqualTo(6.10).Within(Tolerance),
                "Whitespace/punctuation are ignored, so 'A B!G' has only termini → pI 6.10");
            Assert.That(nonStandard, Is.EqualTo(6.10).Within(Tolerance),
                "Non-ionizable residues X/Z contribute no charge → termini-only pI 6.10");
        });
    }

    #endregion
}
