using NUnit.Framework;
using Seqeron.Genomics.Annotation;
using static Seqeron.Genomics.Annotation.SpliceSitePredictor;

namespace Seqeron.Genomics.Tests;

/// <summary>
/// SPLICE-PREDICT-001: Canonical test file for Gene Structure Prediction.
/// Tests <see cref="SpliceSitePredictor.PredictGeneStructure"/>,
/// <see cref="SpliceSitePredictor.PredictIntrons"/>, and internal helpers
/// (SelectNonOverlappingIntrons, DeriveExons, GenerateSplicedSequence).
/// Evidence: Breathnach &amp; Chambon (1981), Burge et al. (1999),
/// Gilbert (1978), Shapiro &amp; Senapathy (1987), Alberts et al. (2002).
/// </summary>
[TestFixture]
public class SpliceSitePredictor_GeneStructure_Tests
{
    // ── Shared test sequences ──────────────────────────────────────────
    // Two-exon gene: exon1(35) + GU-intron(83) + exon2(35) = 153 nt
    // Intron = GUAAGU(6) + 60×A + UUUUUUUUUUUUUU(14) + CAG(3) = 83 nt

    private const string Exon1 = "AUGCCCAAAGGGCCCUUUAAAGGGCCCUUUAAAGC"; // 35 nt
    private const string Donor = "GUAAGU";
    private const string IntronBody = "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"; // 60 A
    private const string Ppt = "UUUUUUUUUUUUUU"; // 14 nt PPT
    private const string Acceptor = "CAG";
    private const string Exon2 = "GCCUUUAAAGGGCCCUUUAAAGGGCCCUUUAAAGC"; // 35 nt

    private static readonly string TwoExonSequence =
        Exon1 + Donor + IntronBody + Ppt + Acceptor + Exon2;

    // IntronPart used for reference calculations
    private static readonly int IntronPartLength = (Donor + IntronBody + Ppt + Acceptor).Length;

    // Single-exon gene: 50 nt, no GT/GU dinucleotide
    private const string SingleExonSequence =
        "AACCAACCAACCAACCAACCAACCAACCAACCAACCAACCAACCAACCAA";

    #region M1: Empty/Null Input Returns Empty Structure — Trivially Correct

    [Test]
    public void PredictGeneStructure_EmptyString_ReturnsEmptyStructure()
    {
        var structure = PredictGeneStructure("");

        Assert.Multiple(() =>
        {
            Assert.That(structure.Exons, Is.Empty,
                "Empty input must produce no exons");
            Assert.That(structure.Introns, Is.Empty,
                "Empty input must produce no introns");
            Assert.That(structure.SplicedSequence, Is.Empty,
                "Empty input must produce empty spliced sequence");
            Assert.That(structure.OverallScore, Is.EqualTo(0),
                "Empty input must produce zero overall score");
        });
    }

    [Test]
    public void PredictGeneStructure_NullInput_ReturnsEmptyStructure()
    {
        var structure = PredictGeneStructure(null!);

        Assert.Multiple(() =>
        {
            Assert.That(structure.Exons, Is.Empty,
                "Null input must produce no exons");
            Assert.That(structure.Introns, Is.Empty,
                "Null input must produce no introns");
        });
    }

    #endregion

    #region M2: Single-Exon Gene (No Introns) — Gilbert (1978)

    [Test]
    public void PredictGeneStructure_NoSpliceSites_SingleExon()
    {
        // Sequence without GU dinucleotide: cannot form any intron
        var structure = PredictGeneStructure(SingleExonSequence, minExonLength: 10, minScore: 0.5);

        Assert.Multiple(() =>
        {
            Assert.That(structure.Introns, Has.Count.EqualTo(0),
                "Sequence without splice sites must have zero introns");
            Assert.That(structure.Exons, Has.Count.EqualTo(1),
                "Single-exon gene must have exactly one exon");
            Assert.That(structure.Exons[0].Type, Is.EqualTo(ExonType.Single),
                "Sole exon must be typed as Single per Gilbert (1978)");
            Assert.That(structure.Exons[0].Start, Is.EqualTo(0),
                "Single exon must start at position 0");
            Assert.That(structure.Exons[0].End, Is.EqualTo(SingleExonSequence.Length - 1),
                "Single exon must end at last position");
        });
    }

    #endregion

    #region M3: Two-Exon Gene with GT-AG Intron — Breathnach & Chambon (1981)

    [Test]
    public void PredictGeneStructure_ClearGTAG_FindsOneIntronTwoExons()
    {
        var structure = PredictGeneStructure(TwoExonSequence, minExonLength: 5, minScore: 0.2);

        Assert.Multiple(() =>
        {
            // Exact counts — designed sequence has one clear GT-AG intron
            Assert.That(structure.Introns, Has.Count.EqualTo(1),
                "Designed two-exon sequence must yield exactly 1 intron");
            Assert.That(structure.Exons, Has.Count.EqualTo(2),
                "One intron splits the sequence into exactly 2 exons");

            var intron = structure.Introns[0];

            // Intron boundaries per GT-AG rule — Breathnach & Chambon (1981)
            Assert.That(intron.Start, Is.EqualTo(Exon1.Length),
                "Intron must start at the GU donor (position after Exon1)");
            Assert.That(intron.End, Is.EqualTo(Exon1.Length + IntronPartLength - 1),
                "Intron must end at the G of AG acceptor");
            Assert.That(intron.Length, Is.EqualTo(IntronPartLength),
                $"Intron length = donor(6) + body(60) + PPT(14) + acceptor(3) = {IntronPartLength} nt");

            // Intron sequence boundaries: starts with GU, ends with AG — S2
            Assert.That(intron.Sequence, Does.StartWith("GU"),
                "Intron must start with GU dinucleotide per GT-AG rule");
            Assert.That(intron.Sequence, Does.EndWith("AG"),
                "Intron must end with AG dinucleotide per GT-AG rule");

            // INV-3: exon + intron positions cover entire sequence
            int coveredPositions = structure.Exons.Sum(e => e.Length) + structure.Introns.Sum(i => i.Length);
            Assert.That(coveredPositions, Is.EqualTo(TwoExonSequence.Length),
                "Sum of exon and intron lengths must equal total sequence length (INV-3)");
        });
    }

    #endregion

    #region M4: Spliced Sequence Excludes Intron — Gilbert (1978)

    [Test]
    public void PredictGeneStructure_SplicedSequence_ExactContent()
    {
        var structure = PredictGeneStructure(TwoExonSequence, minExonLength: 5, minScore: 0.2);

        Assert.That(structure.Introns, Has.Count.GreaterThan(0),
            "Prerequisite: introns must be found for this test to be meaningful");

        Assert.Multiple(() =>
        {
            // INV-5: length = total − intron length
            int totalIntronLength = structure.Introns.Sum(i => i.Length);
            Assert.That(structure.SplicedSequence.Length,
                Is.EqualTo(TwoExonSequence.Length - totalIntronLength),
                "Spliced length = input length − total intron length (INV-5)");

            // INV-4: spliced sequence = concatenation of exon sequences
            string exonConcat = string.Join("", structure.Exons.Select(e => e.Sequence));
            Assert.That(structure.SplicedSequence, Is.EqualTo(exonConcat),
                "Spliced sequence must equal concatenation of exon sequences (INV-4)");

            // Exact expected value: Exon1 + Exon2 — Gilbert (1978)
            Assert.That(structure.SplicedSequence, Is.EqualTo(Exon1 + Exon2),
                "Spliced sequence = Exon1 + Exon2 (intron excised)");
            Assert.That(structure.SplicedSequence.Length, Is.EqualTo(Exon1.Length + Exon2.Length),
                $"Spliced length = {Exon1.Length} + {Exon2.Length} = {Exon1.Length + Exon2.Length} nt");
        });
    }

    #endregion

    #region M5: Intron MinLength Respected — Algorithm Parameter

    [Test]
    public void PredictIntrons_MinLength_FiltersShortCandidates()
    {
        // TwoExonSequence has intron candidates ≥79 nt — found when minIntronLength=60
        var intronsFound = PredictIntrons(TwoExonSequence, minIntronLength: 60, minScore: 0.2).ToList();
        Assert.That(intronsFound, Is.Not.Empty,
            "Introns must be found when minIntronLength (60) ≤ candidate lengths");
        Assert.That(intronsFound.All(i => i.Length >= 60), Is.True,
            "Every returned intron must satisfy Length ≥ minIntronLength (INV-8)");

        // Raise minIntronLength above all candidate lengths → must filter all out
        var intronsFiltered = PredictIntrons(TwoExonSequence, minIntronLength: 200, minScore: 0.2).ToList();
        Assert.That(intronsFiltered, Is.Empty,
            "No introns when minIntronLength (200) > all candidate lengths");
    }

    #endregion

    #region M6: Intron MaxLength Respected — Algorithm Parameter

    [Test]
    public void PredictIntrons_MaxLength_FiltersLongCandidates()
    {
        // TwoExonSequence intron candidates are ≥79 nt — found with default maxIntronLength
        var intronsFound = PredictIntrons(TwoExonSequence, minIntronLength: 60, minScore: 0.2).ToList();
        Assert.That(intronsFound, Is.Not.Empty,
            "Introns must be found with default maxIntronLength (100000)");
        Assert.That(intronsFound.All(i => i.Length <= 100000), Is.True,
            "Every returned intron must satisfy Length ≤ maxIntronLength (INV-9)");

        // Lower maxIntronLength below all candidate lengths → must filter all out
        var intronsFiltered = PredictIntrons(TwoExonSequence, minIntronLength: 60, maxIntronLength: 70, minScore: 0.2).ToList();
        Assert.That(intronsFiltered, Is.Empty,
            "No introns when maxIntronLength (70) < all candidate lengths (≥79)");
    }

    #endregion

    #region M7: Intron Type U2 for GT-AG — Burge et al. (1999)

    [Test]
    public void PredictIntrons_GTAG_ClassifiedAsU2()
    {
        var introns = PredictIntrons(TwoExonSequence, minIntronLength: 50, minScore: 0.2).ToList();

        // Filter to GU-donor introns by checking the intron sequence start
        var guIntrons = introns.Where(i => i.Sequence.StartsWith("GU")).ToList();

        Assert.That(guIntrons, Is.Not.Empty,
            "Should find at least one GU-donor intron in the test sequence");

        // Per Burge et al. (1999): GT-AG introns are spliced by the U2 (major) spliceosome
        foreach (var intron in guIntrons)
        {
            Assert.That(intron.Type, Is.EqualTo(IntronType.U2),
                "GT-AG intron must be classified as U2 per Burge et al. (1999)");
        }
    }

    #endregion

    #region M8: Exon Types Assigned Correctly — Gilbert (1978)

    [Test]
    public void PredictGeneStructure_MultiExon_CorrectExonTypes()
    {
        var structure = PredictGeneStructure(TwoExonSequence, minExonLength: 5, minScore: 0.2);

        Assert.That(structure.Exons, Has.Count.EqualTo(2),
            "Prerequisite: two exons from the designed two-exon sequence");

        Assert.Multiple(() =>
        {
            Assert.That(structure.Exons[0].Type, Is.EqualTo(ExonType.Initial),
                "First exon must be typed as Initial per Gilbert (1978)");
            Assert.That(structure.Exons[^1].Type, Is.EqualTo(ExonType.Terminal),
                "Last exon must be typed as Terminal per Gilbert (1978)");
        });
    }

    #endregion

    #region M9: Exon Phase Tracks Reading Frame — Alberts et al. (2002)

    [Test]
    public void PredictGeneStructure_ExonPhase_CumulativeLengthMod3()
    {
        var structure = PredictGeneStructure(TwoExonSequence, minExonLength: 5, minScore: 0.2);

        Assert.That(structure.Exons, Has.Count.EqualTo(2),
            "Prerequisite: two exons from the designed two-exon sequence");

        Assert.Multiple(() =>
        {
            // First exon phase must be 0 — trivially, no preceding exons
            Assert.That(structure.Exons[0].Phase, Is.EqualTo(0),
                "First exon phase must be 0 per reading frame convention — Alberts et al. (2002)");

            // Second exon phase = Exon1.Length mod 3 = 35 mod 3 = 2
            Assert.That(structure.Exons[1].Phase, Is.EqualTo(Exon1.Length % 3),
                $"Second exon phase = ({Exon1.Length} mod 3) = {Exon1.Length % 3} — cumulative length rule (INV-6)");

            // General formula verification
            int cumulative = 0;
            for (int i = 0; i < structure.Exons.Count; i++)
            {
                Assert.That(structure.Exons[i].Phase, Is.EqualTo(cumulative % 3),
                    $"Exon {i} phase must be (sum of preceding exon lengths) mod 3 = {cumulative % 3} (INV-6)");
                cumulative += structure.Exons[i].Length;
            }
        });
    }

    #endregion

    #region M10: Score Range [0, 1] — Normalization Invariant

    [Test]
    public void PredictIntrons_Scores_InValidRange()
    {
        var introns = PredictIntrons(TwoExonSequence, minIntronLength: 50, minScore: 0.1).ToList();

        Assert.That(introns, Is.Not.Empty,
            "Prerequisite: at least one intron must be found for score validation");

        foreach (var intron in introns)
        {
            Assert.Multiple(() =>
            {
                Assert.That(intron.Score, Is.GreaterThanOrEqualTo(0).And.LessThanOrEqualTo(1),
                    $"Intron score {intron.Score} must be in [0,1] (INV-10)");
                Assert.That(intron.DonorSite.Score, Is.GreaterThanOrEqualTo(0).And.LessThanOrEqualTo(1),
                    $"Donor score must be in [0,1]");
                Assert.That(intron.AcceptorSite.Score, Is.GreaterThanOrEqualTo(0).And.LessThanOrEqualTo(1),
                    $"Acceptor score must be in [0,1]");
            });
        }
    }

    [Test]
    public void PredictGeneStructure_OverallScore_InValidRange()
    {
        var structure = PredictGeneStructure(TwoExonSequence, minExonLength: 5, minScore: 0.2);

        Assert.That(structure.OverallScore, Is.GreaterThanOrEqualTo(0).And.LessThanOrEqualTo(1),
            "Overall gene structure score must be in [0,1]");
    }

    #endregion

    #region S1: Non-Overlapping Intron Selection — Algorithm Invariant

    [Test]
    public void PredictGeneStructure_Introns_DoNotOverlap()
    {
        var structure = PredictGeneStructure(TwoExonSequence, minExonLength: 5, minScore: 0.2);

        var sortedIntrons = structure.Introns.OrderBy(i => i.Start).ToList();
        for (int i = 1; i < sortedIntrons.Count; i++)
        {
            Assert.That(sortedIntrons[i].Start, Is.GreaterThan(sortedIntrons[i - 1].End),
                $"Intron {i} (start={sortedIntrons[i].Start}) must start after intron {i - 1} "
                + $"(end={sortedIntrons[i - 1].End}) — non-overlapping invariant (INV-2)");
        }
    }

    #endregion

    #region S2: DNA T-Equivalence — Implementation Converts T→U

    [Test]
    public void PredictGeneStructure_DNA_SameAsRNA()
    {
        string dnaSequence = TwoExonSequence.Replace('U', 'T');

        var rnaStructure = PredictGeneStructure(TwoExonSequence, minExonLength: 5, minScore: 0.2);
        var dnaStructure = PredictGeneStructure(dnaSequence, minExonLength: 5, minScore: 0.2);

        Assert.Multiple(() =>
        {
            Assert.That(dnaStructure.Introns.Count, Is.EqualTo(rnaStructure.Introns.Count),
                "DNA (T) and RNA (U) input must produce same number of introns");
            Assert.That(dnaStructure.Exons.Count, Is.EqualTo(rnaStructure.Exons.Count),
                "DNA (T) and RNA (U) input must produce same number of exons");
            Assert.That(dnaStructure.SplicedSequence, Is.EqualTo(rnaStructure.SplicedSequence),
                "DNA and RNA input must produce identical spliced sequence (T→U normalization)");
            Assert.That(dnaStructure.OverallScore, Is.EqualTo(rnaStructure.OverallScore).Within(1e-10),
                "DNA and RNA input must produce identical overall score");
        });
    }

    #endregion

    #region S3: Intron Sequence Matches Substring — Definition

    [Test]
    public void PredictIntrons_IntronSequence_MatchesInputSubstring()
    {
        var introns = PredictIntrons(TwoExonSequence, minIntronLength: 50, minScore: 0.2).ToList();
        string upper = TwoExonSequence.ToUpperInvariant().Replace('T', 'U');

        foreach (var intron in introns)
        {
            string expected = upper.Substring(intron.Start, intron.Length);
            Assert.That(intron.Sequence, Is.EqualTo(expected),
                $"Intron.Sequence must equal input[{intron.Start}..{intron.Start + intron.Length}]");
        }
    }

    #endregion

    #region S4: Higher Threshold Filters More — Parameter Semantics

    [Test]
    public void PredictIntrons_HigherThreshold_FewerOrEqualIntrons()
    {
        var intronsLow = PredictIntrons(TwoExonSequence, minScore: 0.2).ToList();
        var intronsHigh = PredictIntrons(TwoExonSequence, minScore: 0.6).ToList();

        Assert.That(intronsHigh.Count, Is.LessThanOrEqualTo(intronsLow.Count),
            "Higher minScore threshold must produce fewer or equal introns than lower threshold");
    }

    #endregion

    #region S5: Case Insensitivity — Implementation Uses ToUpperInvariant

    [Test]
    public void PredictGeneStructure_LowercaseInput_SameResult()
    {
        string lower = TwoExonSequence.ToLowerInvariant();

        var upperResult = PredictGeneStructure(TwoExonSequence, minExonLength: 5, minScore: 0.2);
        var lowerResult = PredictGeneStructure(lower, minExonLength: 5, minScore: 0.2);

        Assert.Multiple(() =>
        {
            Assert.That(lowerResult.Introns.Count, Is.EqualTo(upperResult.Introns.Count),
                "Lowercase input must produce same intron count as uppercase");
            Assert.That(lowerResult.Exons.Count, Is.EqualTo(upperResult.Exons.Count),
                "Lowercase input must produce same exon count as uppercase");
            Assert.That(lowerResult.SplicedSequence, Is.EqualTo(upperResult.SplicedSequence),
                "Lowercase input must produce identical spliced sequence");
            Assert.That(lowerResult.OverallScore, Is.EqualTo(upperResult.OverallScore).Within(1e-10),
                "Lowercase input must produce identical overall score");
        });
    }

    #endregion

    #region C1: Overall Score = Mean of Intron Scores — Implementation Invariant

    [Test]
    public void PredictGeneStructure_OverallScore_EqualsMeanOfIntronScores()
    {
        var structure = PredictGeneStructure(TwoExonSequence, minExonLength: 5, minScore: 0.2);

        Assert.That(structure.Introns, Has.Count.GreaterThan(0),
            "Prerequisite: introns must be found for mean score verification");

        double expectedMean = structure.Introns.Average(i => i.Score);
        Assert.That(structure.OverallScore, Is.EqualTo(expectedMean).Within(1e-10),
            "OverallScore must equal arithmetic mean of intron scores (INV-7)");
    }

    #endregion

    #region C2: No Introns → Spliced Equals Original — Identity

    [Test]
    public void PredictGeneStructure_NoIntrons_SplicedEqualsInput()
    {
        var structure = PredictGeneStructure(SingleExonSequence, minExonLength: 10, minScore: 0.5);

        Assert.That(structure.Introns, Has.Count.EqualTo(0),
            "SingleExonSequence (no GU dinucleotide) must yield zero introns");

        string expectedSpliced = SingleExonSequence.ToUpperInvariant().Replace('T', 'U');
        Assert.That(structure.SplicedSequence, Is.EqualTo(expectedSpliced),
            "When no introns exist, spliced sequence must equal the (uppercased RNA) input");
    }

    #endregion

    #region PredictIntrons: Empty/Null → Empty — Guard Clause

    [Test]
    public void PredictIntrons_EmptyInput_ReturnsEmpty()
    {
        var introns = PredictIntrons("").ToList();
        Assert.That(introns, Is.Empty, "Empty input must yield no introns");
    }

    [Test]
    public void PredictIntrons_NullInput_ReturnsEmpty()
    {
        var introns = PredictIntrons(null!).ToList();
        Assert.That(introns, Is.Empty, "Null input must yield no introns");
    }

    #endregion
}
