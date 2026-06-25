using System;
using System.Collections.Generic;
using System.Linq;
using NUnit.Framework;
using static Seqeron.Genomics.Annotation.VariantAnnotator;
using VariantType = Seqeron.Genomics.Annotation.VariantAnnotator.VariantType;
using Variant = Seqeron.Genomics.Annotation.VariantAnnotator.Variant;

namespace Seqeron.Genomics.Tests;

/// <summary>
/// VARIANT-ANNOT-001 mutation killers (batch 6): the closed-form numeric helpers the canonical
/// fixture only smoke-tested — variant-type classification / left-normalization boundaries, the
/// ACMG-style evidence-point accumulation, classification thresholds and confidence formula
/// (Richards et al. 2015, ACMG/AMP), conserved-element segmentation (PhastCons), regulatory
/// overlap, and PWM/IUPAC transcription-factor motif scoring.
/// </summary>
[TestFixture]
public class VariantAnnotator_MutationKillers6_Tests
{
    private const double Tol = 1e-9;

    #region ClassifyVariant / NormalizeVariant boundaries

    [Test]
    public void ClassifyVariant_RefLongerButNotPureDeletion_IsIndel()
        // ref longer, alt length 1 but alt is NOT a prefix of ref ⇒ Indel, not Deletion.
        => Assert.That(ClassifyVariant("ATG", "C"), Is.EqualTo(VariantType.Indel));

    [Test]
    public void ClassifyVariant_AltLongerButNotPureInsertion_IsIndel()
        // alt longer, ref length 1 but ref is NOT a prefix of alt ⇒ Indel, not Insertion.
        => Assert.That(ClassifyVariant("C", "ATG"), Is.EqualTo(VariantType.Indel));

    [Test]
    public void NormalizeVariant_SuffixTrimStopsWhenAlternateReachesLengthOne()
    {
        // AAG/AG: trim suffix G ⇒ AA/A; trimming must STOP (alt now length 1) ⇒ deletion AA→A at 100.
        var v = NormalizeVariant("chr1", 100, "AAG", "AG");
        Assert.That(v.Position, Is.EqualTo(100));
        Assert.That(v.Reference, Is.EqualTo("AA"));
        Assert.That(v.Alternate, Is.EqualTo("A"));
        Assert.That(v.Type, Is.EqualTo(VariantType.Deletion));
    }

    [Test]
    public void NormalizeVariant_SuffixTrimStopsWhenReferenceReachesLengthOne()
    {
        // GCC/TGCC: trim suffix C,C ⇒ G/TG; ref now length 1 ⇒ STOP ⇒ G→TG (not a clean insertion prefix) = Indel.
        var v = NormalizeVariant("chr1", 100, "GCC", "TGCC");
        Assert.That(v.Position, Is.EqualTo(100));
        Assert.That(v.Reference, Is.EqualTo("G"));
        Assert.That(v.Alternate, Is.EqualTo("TG"));
        Assert.That(v.Type, Is.EqualTo(VariantType.Indel));
    }

    [Test]
    public void NormalizeVariant_PrefixTrimStopsWhenReferenceReachesLengthOne()
    {
        // AC/ACG: trim prefix A (pos→101) ⇒ C/CG; ref now length 1 ⇒ STOP ⇒ insertion at 101.
        var v = NormalizeVariant("chr1", 100, "AC", "ACG");
        Assert.That(v.Position, Is.EqualTo(101));
        Assert.That(v.Reference, Is.EqualTo("C"));
        Assert.That(v.Alternate, Is.EqualTo("CG"));
        Assert.That(v.Type, Is.EqualTo(VariantType.Insertion));
    }

    [Test]
    public void NormalizeVariant_PrefixTrimStopsWhenAlternateReachesLengthOne()
    {
        // CGA/CG: trim prefix C (pos→101) ⇒ GA/G; alt now length 1 ⇒ STOP ⇒ deletion at 101.
        var v = NormalizeVariant("chr1", 100, "CGA", "CG");
        Assert.That(v.Position, Is.EqualTo(101));
        Assert.That(v.Reference, Is.EqualTo("GA"));
        Assert.That(v.Alternate, Is.EqualTo("G"));
        Assert.That(v.Type, Is.EqualTo(VariantType.Deletion));
    }

    #endregion

    #region PredictPathogenicity — evidence points

    private static VariantAnnotation Ann(ConsequenceType c, ImpactLevel imp, double? sift = null, double? poly = null) =>
        new(new Variant("chr1", 100, "A", "G", VariantType.SNV),
            "ENST", "ENSG", "GENE", c, imp, null, null, null, null, sift, poly, null, null, null);

    // Net evidence points (= pathogenicPoints − benignPoints, surfaced as ClinicalSignificance)
    // driven purely by population frequency on a missense variant with no other evidence.
    // PM2 extreme < 0.01% ⇒ +2; PM2 rare < 0.1% ⇒ +1; BS1 > 1% ⇒ −2; BA1 > 5% ⇒ −4 (Richards 2015).
    [TestCase(0.00005, 2)]   // < 0.0001 ⇒ PM2 strong
    [TestCase(0.0001, 1)]    // boundary: NOT < 0.0001 ⇒ falls to < 0.001 ⇒ +1
    [TestCase(0.0005, 1)]    // < 0.001 ⇒ PM2
    [TestCase(0.001, 0)]     // boundary: NOT < 0.001, NOT > 0.05/0.01 ⇒ 0
    [TestCase(0.01, 0)]      // boundary: NOT > 0.01 ⇒ 0
    [TestCase(0.02, -2)]     // > 0.01 ⇒ BS1
    [TestCase(0.05, -2)]     // boundary: NOT > 0.05 ⇒ falls to > 0.01 ⇒ BS1
    [TestCase(0.1, -4)]      // > 0.05 ⇒ BA1
    public void PredictPathogenicity_PopulationFrequencyPoints(double freq, int expectedNet)
    {
        var p = PredictPathogenicity(Ann(ConsequenceType.MissenseVariant, ImpactLevel.Moderate), populationFrequency: freq);
        Assert.That(p.ClinicalSignificance!.Value, Is.EqualTo(expectedNet).Within(Tol));
    }

    [TestCase(0.04, 1)]   // SIFT < 0.05 ⇒ PP3 deleterious ⇒ +1
    [TestCase(0.05, 0)]   // boundary: NOT < 0.05 ⇒ 0
    public void PredictPathogenicity_SiftDeleteriousPoint(double sift, int expectedNet)
    {
        var p = PredictPathogenicity(Ann(ConsequenceType.MissenseVariant, ImpactLevel.Moderate, sift: sift));
        Assert.That(p.ClinicalSignificance!.Value, Is.EqualTo(expectedNet).Within(Tol));
    }

    [TestCase(0.95, 1)]    // PolyPhen > 0.908 ⇒ PP3 probably damaging ⇒ +1
    [TestCase(0.908, 0)]   // boundary: NOT > 0.908 ⇒ 0
    public void PredictPathogenicity_PolyphenDamagingPoint(double poly, int expectedNet)
    {
        var p = PredictPathogenicity(Ann(ConsequenceType.MissenseVariant, ImpactLevel.Moderate, poly: poly));
        Assert.That(p.ClinicalSignificance!.Value, Is.EqualTo(expectedNet).Within(Tol));
    }

    [TestCase(0.96, 0.05, -1)]  // BP4: SIFT > 0.95 AND PolyPhen < 0.1 ⇒ benign +1
    [TestCase(0.95, 0.05, 0)]   // boundary: SIFT NOT > 0.95 ⇒ no BP4
    [TestCase(0.96, 0.1, 0)]    // boundary: PolyPhen NOT < 0.1 ⇒ no BP4
    public void PredictPathogenicity_Bp4ComputationalBenign(double sift, double poly, int expectedNet)
    {
        var p = PredictPathogenicity(Ann(ConsequenceType.MissenseVariant, ImpactLevel.Moderate, sift: sift, poly: poly));
        Assert.That(p.ClinicalSignificance!.Value, Is.EqualTo(expectedNet).Within(Tol));
    }

    [Test]
    public void PredictPathogenicity_Bp7SynonymousBenignPoint()
    {
        var p = PredictPathogenicity(Ann(ConsequenceType.SynonymousVariant, ImpactLevel.Low));
        Assert.That(p.ClinicalSignificance!.Value, Is.EqualTo(-1).Within(Tol)); // BP7 ⇒ benign +1
    }

    [TestCase(5.0, 1)]    // conservation > 4.0 ⇒ PP3 conserved ⇒ +1
    [TestCase(4.0, 0)]    // boundary: NOT > 4.0
    [TestCase(-3.0, -1)]  // conservation < -2.0 ⇒ BP4 non-conserved ⇒ benign +1
    [TestCase(-2.0, 0)]   // boundary: NOT < -2.0
    public void PredictPathogenicity_ConservationPoints(double cons, int expectedNet)
    {
        var p = PredictPathogenicity(Ann(ConsequenceType.MissenseVariant, ImpactLevel.Moderate), conservationScore: cons);
        Assert.That(p.ClinicalSignificance!.Value, Is.EqualTo(expectedNet).Within(Tol));
    }

    [Test]
    public void PredictPathogenicity_ClinVarRequiresBothInClinvarAndSignificance()
    {
        // significance present but inClinvar=false ⇒ NO points (the guard is a conjunction).
        var noFlag = PredictPathogenicity(Ann(ConsequenceType.MissenseVariant, ImpactLevel.Moderate),
            inClinvar: false, clinvarSignificance: "Pathogenic");
        Assert.That(noFlag.ClinicalSignificance!.Value, Is.EqualTo(0).Within(Tol));

        var pathogenic = PredictPathogenicity(Ann(ConsequenceType.MissenseVariant, ImpactLevel.Moderate),
            inClinvar: true, clinvarSignificance: "Pathogenic");
        Assert.That(pathogenic.ClinicalSignificance!.Value, Is.EqualTo(4).Within(Tol)); // PP5 +4

        var benign = PredictPathogenicity(Ann(ConsequenceType.MissenseVariant, ImpactLevel.Moderate),
            inClinvar: true, clinvarSignificance: "Benign");
        Assert.That(benign.ClinicalSignificance!.Value, Is.EqualTo(-4).Within(Tol)); // BP6 −4
    }

    [Test]
    public void PredictPathogenicity_FunctionalLofPoint()
    {
        var p = PredictPathogenicity(Ann(ConsequenceType.MissenseVariant, ImpactLevel.Moderate),
            functionalEvidence: new[] { "LOF functional assay" });
        Assert.That(p.ClinicalSignificance!.Value, Is.EqualTo(4).Within(Tol)); // PS3 +4
    }

    #endregion

    #region ClassifyByPoints thresholds + CalculateConfidence

    [Test]
    public void Classify_NetExactly10_IsPathogenic()
    {
        // PVS1 high impact(8) + PM2 extreme(2) ⇒ net 10 ⇒ Pathogenic (inclusive ≥ 10).
        var p = PredictPathogenicity(Ann(ConsequenceType.StopGained, ImpactLevel.High), populationFrequency: 0.00005);
        Assert.That(p.ClinicalSignificance!.Value, Is.EqualTo(10).Within(Tol));
        Assert.That(p.Classification, Is.EqualTo(PathogenicityClass.Pathogenic));
    }

    [Test]
    public void Classify_PathogenicPointsAtLeast10EvenWhenNetBelow10_IsPathogenic()
    {
        // High(8)+PM2(2)=10 pathogenic points; conservation −3 ⇒ benign 1 ⇒ net 9, but pathogenic ≥ 10 ⇒ Pathogenic.
        var p = PredictPathogenicity(Ann(ConsequenceType.StopGained, ImpactLevel.High),
            populationFrequency: 0.00005, conservationScore: -3.0);
        Assert.That(p.ClinicalSignificance!.Value, Is.EqualTo(9).Within(Tol));
        Assert.That(p.Classification, Is.EqualTo(PathogenicityClass.Pathogenic));
    }

    [Test]
    public void Classify_NetExactly6_IsLikelyPathogenic_AndConfidenceFormula()
    {
        // High(8) pathogenic + BS1(2) benign ⇒ net 6 ⇒ Likely Pathogenic.
        // confidence = min(0.99, 0.5 + |8−2|/(8+2) * 0.5) = 0.5 + 0.6*0.5 = 0.8.
        var p = PredictPathogenicity(Ann(ConsequenceType.StopGained, ImpactLevel.High), populationFrequency: 0.02);
        Assert.That(p.ClinicalSignificance!.Value, Is.EqualTo(6).Within(Tol));
        Assert.That(p.Classification, Is.EqualTo(PathogenicityClass.LikelyPathogenic));
        Assert.That(p.ConfidenceScore, Is.EqualTo(0.8).Within(Tol));
        Assert.That(p.IsActionable, Is.True);
    }

    [Test]
    public void Classify_NetExactlyMinus6_IsBenign_AndConfidenceClamped()
    {
        // BA1(4) + BP7 synonymous(1) + non-conserved(1) ⇒ benign 6, net −6 ⇒ Benign.
        // confidence = min(0.99, 0.5 + 6/6*0.5) = min(0.99, 1.0) = 0.99 (clamp).
        var p = PredictPathogenicity(Ann(ConsequenceType.SynonymousVariant, ImpactLevel.Low),
            populationFrequency: 0.1, conservationScore: -3.0);
        Assert.That(p.ClinicalSignificance!.Value, Is.EqualTo(-6).Within(Tol));
        Assert.That(p.Classification, Is.EqualTo(PathogenicityClass.Benign));
        Assert.That(p.ConfidenceScore, Is.EqualTo(0.99).Within(Tol));
        Assert.That(p.IsActionable, Is.False);
    }

    #endregion

    #region FindConservedElements (PhastCons segmentation)

    private static IEnumerable<ConservationScore> Run(string chrom, int start, int end, double phastCons)
    {
        for (int p = start; p <= end; p++)
            yield return new ConservationScore(chrom, p, 0, phastCons, 0, 0);
    }

    [Test]
    public void FindConservedElements_TwoBlocksSplitByGapOver10bp()
    {
        // Two ≥20-bp conserved runs separated by a 20-bp gap (> 10) ⇒ two distinct elements.
        var scores = Run("chr1", 1, 20, 0.9).Concat(Run("chr1", 40, 59, 0.9)).ToList();
        var elems = FindConservedElements(scores, threshold: 0.8, minLength: 20).ToList();

        Assert.That(elems, Has.Count.EqualTo(2));
        Assert.That(elems[0].Start, Is.EqualTo(1));
        Assert.That(elems[0].End, Is.EqualTo(20));
        Assert.That(elems[0].Score, Is.EqualTo(0.9).Within(Tol));
        Assert.That(elems[1].Start, Is.EqualTo(40));
        Assert.That(elems[1].End, Is.EqualTo(59));
    }

    [Test]
    public void FindConservedElements_ThresholdIsInclusive()
    {
        // PhastCons exactly == threshold counts as conserved (≥, not >).
        var scores = Run("chr1", 1, 20, 0.9).ToList();
        var elems = FindConservedElements(scores, threshold: 0.9, minLength: 20).ToList();
        Assert.That(elems, Has.Count.EqualTo(1));
        Assert.That(elems[0].Start, Is.EqualTo(1));
        Assert.That(elems[0].End, Is.EqualTo(20));
    }

    [Test]
    public void FindConservedElements_MinLengthIsInclusiveBoundary()
    {
        var twenty = Run("chr1", 1, 20, 0.9).ToList(); // length exactly 20
        Assert.That(FindConservedElements(twenty, threshold: 0.8, minLength: 20).Count(), Is.EqualTo(1)); // 20 ≥ 20 ⇒ kept
        Assert.That(FindConservedElements(twenty, threshold: 0.8, minLength: 21).Count(), Is.EqualTo(0)); // 20 < 21 ⇒ dropped

        var nineteen = Run("chr1", 1, 19, 0.9).ToList(); // length 19
        Assert.That(FindConservedElements(nineteen, threshold: 0.8, minLength: 20).Count(), Is.EqualTo(0));
    }

    #endregion

    #region AnnotateRegulatoryElements (overlap)

    [Test]
    public void AnnotateRegulatoryElements_ReportsOnlyOverlappingRegionsOnTheSameChromosome()
    {
        // 2-bp variant at 100..101 (ref length 2). Overlap iff varEnd(101) ≥ region.Start AND varStart(100) ≤ region.End.
        var v = new Variant("chr1", 100, "AT", "GC", VariantType.MNV);
        var regions = new (string, int, int, string, string?, double?, IReadOnlyList<string>)[]
        {
            ("chr1", 101, 200, "enhancer", null, null, new List<string>()), // touches varEnd ⇒ overlap
            ("chr1", 102, 200, "enhancer", null, null, new List<string>()), // starts past varEnd ⇒ no overlap
            ("chr1",  50, 100, "promoter", null, null, new List<string>()), // ends at varStart ⇒ overlap
            ("chr1",  50,  99, "promoter", null, null, new List<string>()), // ends before varStart ⇒ no overlap
            ("chr2", 100, 200, "silencer", null, null, new List<string>()), // wrong chromosome ⇒ skipped
        };

        var hits = AnnotateRegulatoryElements(v, regions).ToList();
        Assert.That(hits.Select(h => h.Start).OrderBy(s => s), Is.EqualTo(new[] { 50, 101 }));
    }

    #endregion

    #region PredictTfBindingChange / ScoreMotif (PWM + IUPAC)

    [Test]
    public void PredictTfBindingChange_PlainMotif_ExactRefAndAltScores()
    {
        // Reference contains a perfect ATG at offset 4..6; the SNV at offset 5 (T→C) breaks the middle base.
        var v = new Variant("chr1", 100, "T", "C", VariantType.SNV);
        var motifs = new[] { ("TF1", "ATG", 0.5) };
        var result = PredictTfBindingChange(v, motifs, "AAACATGCAAA", contextOffset: 5).Single();

        Assert.That(result.RefScore, Is.EqualTo(1.0).Within(Tol));       // perfect 3/3 match
        Assert.That(result.AltScore, Is.EqualTo(2.0 / 3.0).Within(Tol)); // 2/3 after disruption
        Assert.That(result.ScoreDifference, Is.EqualTo(1.0 / 3.0).Within(Tol));
    }

    [Test]
    public void PredictTfBindingChange_IupacMotif_AllDegenerateCodesMatch()
    {
        // Motif RYWSN scores a perfect 5/5 against ACAGT (R=A, Y=C, W=A, S=G, N=T); the SNV (G→A)
        // breaks only the S position ⇒ 4/5.
        var v = new Variant("chr1", 100, "G", "A", VariantType.SNV);
        var motifs = new[] { ("TF2", "RYWSN", 0.5) };
        var result = PredictTfBindingChange(v, motifs, "TACAGTTTT", contextOffset: 4).Single();

        Assert.That(result.RefScore, Is.EqualTo(1.0).Within(Tol));
        Assert.That(result.AltScore, Is.EqualTo(0.8).Within(Tol));
    }

    [Test]
    public void PredictTfBindingChange_OffsetAtOrPastContextEnd_YieldsNothing()
    {
        var v = new Variant("chr1", 100, "A", "C", VariantType.SNV);
        var motifs = new[] { ("TF3", "ATG", 0.5) };
        Assert.That(PredictTfBindingChange(v, motifs, "AAAAA", contextOffset: 5), Is.Empty);
    }

    [Test]
    public void PredictTfBindingChange_NonSnvYieldsNothing()
    {
        // A non-SNV must short-circuit before any scoring, even when its alternate base would disrupt the motif.
        var v = new Variant("chr1", 100, "AG", "A", VariantType.Deletion);
        var motifs = new[] { ("TF4", "RYWSN", 0.5) };
        Assert.That(PredictTfBindingChange(v, motifs, "TACAGTTTT", contextOffset: 4), Is.Empty);
    }

    #endregion
}
