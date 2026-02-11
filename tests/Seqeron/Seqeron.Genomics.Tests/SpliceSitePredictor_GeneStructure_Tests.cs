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
    public void PredictGeneStructure_ClearGTAG_FindsOneIntron()
    {
        var structure = PredictGeneStructure(TwoExonSequence, minExonLength: 5, minScore: 0.2);

        Assert.That(structure.Introns, Has.Count.GreaterThanOrEqualTo(1),
            "Sequence with clear GT-AG intron (GUAAGU...PPT...CAG) must detect at least one intron "
            + "per Breathnach & Chambon (1981)");
    }

    [Test]
    public void PredictGeneStructure_ClearGTAG_ProducesTwoExons()
    {
        var structure = PredictGeneStructure(TwoExonSequence, minExonLength: 5, minScore: 0.2);

        if (structure.Introns.Count > 0)
        {
            Assert.That(structure.Exons, Has.Count.GreaterThanOrEqualTo(2),
                "Gene with at least one intron must have at least two exons");
        }
    }

    #endregion

    #region M4: Spliced Sequence Excludes Intron — Gilbert (1978)

    [Test]
    public void PredictGeneStructure_SplicedSequence_ShorterThanOriginal()
    {
        var structure = PredictGeneStructure(TwoExonSequence, minExonLength: 5, minScore: 0.2);

        if (structure.Introns.Count > 0)
        {
            // Intron.Length = End - Start, but GenerateSplicedSequence removes
            // positions Start..End inclusive = (End - Start + 1) characters.
            int totalCharsRemoved = structure.Introns.Sum(i => i.End - i.Start + 1);
            string upperInput = TwoExonSequence.ToUpperInvariant().Replace('T', 'U');

            Assert.Multiple(() =>
            {
                Assert.That(structure.SplicedSequence.Length,
                    Is.EqualTo(upperInput.Length - totalCharsRemoved),
                    "Spliced sequence length must equal input length minus total removed characters (INV-5)");
                Assert.That(structure.SplicedSequence.Length, Is.LessThan(upperInput.Length),
                    "Spliced sequence must be shorter than original when introns are removed");
            });
        }
    }

    #endregion

    #region M5: Intron MinLength Respected — Algorithm Parameter

    [Test]
    public void PredictIntrons_MinLength_FiltersShortCandidates()
    {
        // Short intron-like pattern: GT...AG with only ~10 nt between
        const string shortIntron = "CAGGUAAGUAAAAAAAAAAAACAGGGCCCC";
        var introns = PredictIntrons(shortIntron, minIntronLength: 50).ToList();

        Assert.That(introns.All(i => i.Length >= 50), Is.True,
            "No intron shorter than minIntronLength must be returned (INV-8)");
    }

    [Test]
    public void PredictIntrons_LargerMinLength_EmptyForShortSequence()
    {
        // Full intron is ~83 nt, set minIntronLength above that
        var introns = PredictIntrons(TwoExonSequence, minIntronLength: 200, minScore: 0.2).ToList();

        Assert.That(introns, Is.Empty,
            "When minIntronLength exceeds available intron length, no introns should be returned");
    }

    #endregion

    #region M6: Intron MaxLength Respected — Algorithm Parameter

    [Test]
    public void PredictIntrons_MaxLength_FiltersLongCandidates()
    {
        string longBody = new string('A', 600);
        string sequence = $"CAGGUAAGU{longBody}UUUUUUUUUUUUUUCAGG";

        var introns = PredictIntrons(sequence, maxIntronLength: 500, minScore: 0.1).ToList();

        Assert.That(introns.All(i => i.Length <= 500), Is.True,
            "No intron longer than maxIntronLength must be returned (INV-9)");
    }

    #endregion

    #region M7: Intron Type U2 for GT-AG — Burge et al. (1999)

    [Test]
    public void PredictIntrons_GTAG_NotClassifiedAsU12()
    {
        var introns = PredictIntrons(TwoExonSequence, minIntronLength: 50, minScore: 0.2).ToList();

        // Filter to introns with canonical GU donor (SpliceSiteType.Donor, not U12Donor)
        var canonicalIntrons = introns
            .Where(i => i.DonorSite.Type == SpliceSiteType.Donor
                     && i.AcceptorSite.Type == SpliceSiteType.Acceptor)
            .ToList();

        Assert.That(canonicalIntrons, Is.Not.Empty,
            "Should find at least one intron with canonical GU donor and AG acceptor");

        // Per Burge et al. (1999): GT-AG introns belong to the U2 (major) spliceosome.
        // Implementation note: DetermineIntronType classifies by donor/acceptor Type fields;
        // canonical GU-AG introns with Type=Donor are never classified as U12.
        foreach (var intron in canonicalIntrons)
        {
            Assert.That(intron.Type, Is.Not.EqualTo(IntronType.U12),
                "Canonical GU-AG intron must not be classified as U12 (minor spliceosome)");
        }
    }

    #endregion

    #region M8: Exon Types Assigned Correctly — Gilbert (1978)

    [Test]
    public void PredictGeneStructure_MultiExon_CorrectExonTypes()
    {
        var structure = PredictGeneStructure(TwoExonSequence, minExonLength: 5, minScore: 0.2);

        if (structure.Exons.Count >= 2)
        {
            Assert.Multiple(() =>
            {
                Assert.That(structure.Exons[0].Type, Is.EqualTo(ExonType.Initial),
                    "First exon must be typed as Initial per Gilbert (1978)");
                Assert.That(structure.Exons[^1].Type, Is.EqualTo(ExonType.Terminal),
                    "Last exon must be typed as Terminal per Gilbert (1978)");
            });
        }
    }

    #endregion

    #region M9: Exon Phase Tracks Reading Frame — Alberts et al. (2002)

    [Test]
    public void PredictGeneStructure_ExonPhase_CumulativeLengthMod3()
    {
        var structure = PredictGeneStructure(TwoExonSequence, minExonLength: 5, minScore: 0.2);

        if (structure.Exons.Count >= 2)
        {
            // First exon phase must be 0 (start of reading frame)
            Assert.That(structure.Exons[0].Phase, Is.EqualTo(0),
                "First exon phase must be 0 per reading frame convention — Alberts et al. (2002)");

            // Subsequent exon phases must equal cumulative preceding exon length mod 3
            int cumulative = 0;
            for (int i = 0; i < structure.Exons.Count; i++)
            {
                Assert.That(structure.Exons[i].Phase, Is.EqualTo(cumulative % 3),
                    $"Exon {i} phase must be (sum of preceding exon lengths) mod 3 = {cumulative % 3} (INV-6)");
                cumulative += structure.Exons[i].Length;
            }
        }
    }

    #endregion

    #region M10: Score Range [0, 1] — Normalization Invariant

    [Test]
    public void PredictIntrons_Scores_InValidRange()
    {
        var introns = PredictIntrons(TwoExonSequence, minIntronLength: 50, minScore: 0.1).ToList();

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
        });
    }

    #endregion

    #region C1: Overall Score = Mean of Intron Scores — Implementation Invariant

    [Test]
    public void PredictGeneStructure_OverallScore_EqualsMeanOfIntronScores()
    {
        var structure = PredictGeneStructure(TwoExonSequence, minExonLength: 5, minScore: 0.2);

        if (structure.Introns.Count > 0)
        {
            double expectedMean = structure.Introns.Average(i => i.Score);
            Assert.That(structure.OverallScore, Is.EqualTo(expectedMean).Within(1e-10),
                "OverallScore must equal arithmetic mean of intron scores (INV-7)");
        }
        else
        {
            Assert.That(structure.OverallScore, Is.EqualTo(0),
                "OverallScore must be 0 when no introns are found");
        }
    }

    #endregion

    #region C2: No Introns → Spliced Equals Original — Identity

    [Test]
    public void PredictGeneStructure_NoIntrons_SplicedEqualsInput()
    {
        var structure = PredictGeneStructure(SingleExonSequence, minExonLength: 10, minScore: 0.5);

        if (structure.Introns.Count == 0)
        {
            string expectedSpliced = SingleExonSequence.ToUpperInvariant().Replace('T', 'U');
            Assert.That(structure.SplicedSequence, Is.EqualTo(expectedSpliced),
                "When no introns exist, spliced sequence must equal the (uppercased RNA) input");
        }
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
