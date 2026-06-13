// VARIANT-ANNOT-001 — Variant Annotation (functional impact / SO consequences)
// Evidence: docs/Evidence/VARIANT-ANNOT-001-Evidence.md
// TestSpec: tests/TestSpecs/VARIANT-ANNOT-001.md
// Source: McLaren W, et al. (2016) The Ensembl Variant Effect Predictor, Genome Biology 17:122,
//         https://doi.org/10.1186/s13059-016-0974-4; Ensembl ensembl-variation release/110
//         Utils/Constants.pm (impact/rank) and Utils/VariationEffect.pm (predicates);
//         NCBI gc.prt Standard code (transl_table 1).
using NUnit.Framework;
using Seqeron.Genomics;
using System;
using System.Collections.Generic;
using System.Globalization;
using System.Linq;
using System.Threading;

namespace Seqeron.Genomics.Tests;

[TestFixture]
public class VariantAnnotator_FunctionalImpact_Tests
{
    // Single-exon transcript on '+' strand. Exon (90..130) wraps a CDS (100..120),
    // so coding variants sit inside the exon (no splice-site trigger) and the CDS
    // begins at genomic 100. referenceSequence[0] == genomic 100 (sequenceStart=100).
    // Reference CDS layout (codon : genomic) used by the coding tests:
    //   codon1 ATG (100-102)  codon2 CAA (103-105)  codon3 GAA (106-108)
    //   codon4 TTA (109-111)  codon5 TAA (112-114)  codon6 NAA (115-117)  codon7 GGG (118-120)
    private const string RefWindow = "ATGCAAGAATTATAANAAGGG"; // 21 bases, genomic 100..120
    private const int RefStart = 100;

    private static VariantAnnotator.Transcript CodingTranscript()
    {
        var exons = new List<(int, int)> { (90, 130) };
        var codingExons = new List<(int, int)> { (100, 120) };
        return new VariantAnnotator.Transcript(
            "ENST_TEST", "ENSG_TEST", "GENE_TEST", "chr1",
            90, 130, '+', exons, codingExons, 100, 120);
    }

    private static VariantAnnotator.FunctionalImpact Predict(int position, string refAllele, string altAllele,
        VariantAnnotator.VariantType type)
    {
        var variant = new VariantAnnotator.Variant("chr1", position, refAllele, altAllele, type);
        return VariantAnnotator.PredictFunctionalImpact(variant, CodingTranscript(), RefWindow, RefStart);
    }

    #region PredictFunctionalImpact — coding substitutions

    // M1 — GAA (Glu/E, codon3 @106-108) base2 A>T -> GTA (Val/V): missense, MODERATE. Src: VariationEffect.pm + gc.prt
    [Test]
    public void PredictFunctionalImpact_MissenseSnv_ReturnsMissenseModerate()
    {
        var fi = Predict(107, "A", "T", VariantAnnotator.VariantType.SNV);

        Assert.Multiple(() =>
        {
            Assert.That(fi.Consequence, Is.EqualTo(VariantAnnotator.ConsequenceType.MissenseVariant),
                "GAA(Glu)->GTA(Val) changes the amino acid with preserved length => missense_variant");
            Assert.That(fi.Impact, Is.EqualTo(VariantAnnotator.ImpactLevel.Moderate),
                "missense_variant impact is MODERATE in Constants.pm");
            Assert.That(fi.AminoAcidChange, Is.EqualTo("p.E3V"),
                "ref Glu (E) at protein position 3 becomes Val (V)");
        });
    }

    // M2 — TTA (Leu/L, codon4 @109-111) base3 A>G -> TTG (Leu/L): synonymous, LOW. Src: VariationEffect.pm + gc.prt
    [Test]
    public void PredictFunctionalImpact_SynonymousSnv_ReturnsSynonymousLow()
    {
        var fi = Predict(111, "A", "G", VariantAnnotator.VariantType.SNV);

        Assert.Multiple(() =>
        {
            Assert.That(fi.Consequence, Is.EqualTo(VariantAnnotator.ConsequenceType.SynonymousVariant),
                "TTA(Leu)->TTG(Leu) leaves the amino acid unchanged => synonymous_variant");
            Assert.That(fi.Impact, Is.EqualTo(VariantAnnotator.ImpactLevel.Low),
                "synonymous_variant impact is LOW in Constants.pm");
            Assert.That(fi.AminoAcidChange, Is.EqualTo("p.L4="),
                "synonymous change at protein position 4 (Leu)");
        });
    }

    // M3 — CAA (Gln/Q, codon2 @103-105) base1 C>T -> TAA (Stop): stop_gained, HIGH. Src: VariationEffect.pm + gc.prt
    [Test]
    public void PredictFunctionalImpact_PrematureStop_ReturnsStopGainedHigh()
    {
        var fi = Predict(103, "C", "T", VariantAnnotator.VariantType.SNV);

        Assert.Multiple(() =>
        {
            Assert.That(fi.Consequence, Is.EqualTo(VariantAnnotator.ConsequenceType.StopGained),
                "CAA(Gln)->TAA(Stop) introduces a premature stop absent from the reference => stop_gained");
            Assert.That(fi.Impact, Is.EqualTo(VariantAnnotator.ImpactLevel.High),
                "stop_gained impact is HIGH in Constants.pm");
            Assert.That(fi.AminoAcidChange, Is.EqualTo("p.Q2*"),
                "Gln (Q) at protein position 2 becomes a stop (*)");
        });
    }

    // M4 — TAA (Stop, codon5 @112-114) base1 T>C -> CAA (Gln/Q): stop_lost, HIGH. Src: VariationEffect.pm + gc.prt
    [Test]
    public void PredictFunctionalImpact_StopLoss_ReturnsStopLostHigh()
    {
        var fi = Predict(112, "T", "C", VariantAnnotator.VariantType.SNV);

        Assert.Multiple(() =>
        {
            Assert.That(fi.Consequence, Is.EqualTo(VariantAnnotator.ConsequenceType.StopLost),
                "TAA(Stop)->CAA(Gln) removes the reference stop => stop_lost");
            Assert.That(fi.Impact, Is.EqualTo(VariantAnnotator.ImpactLevel.High),
                "stop_lost impact is HIGH in Constants.pm");
        });
    }

    // M5 — ATG (Met/M, start codon1 @100-102) base3 G>C -> ATC (Ile, not a start): start_lost, HIGH.
    // Src: VariationEffect.pm start_lost precedence; gc.prt starts {ATG,TTG,CTG}.
    [Test]
    public void PredictFunctionalImpact_StartCodonChange_ReturnsStartLostHigh()
    {
        var fi = Predict(102, "G", "C", VariantAnnotator.VariantType.SNV);

        Assert.Multiple(() =>
        {
            Assert.That(fi.Consequence, Is.EqualTo(VariantAnnotator.ConsequenceType.StartLost),
                "ATG->ATC at the canonical start codon disrupts initiation => start_lost (overrides missense)");
            Assert.That(fi.Impact, Is.EqualTo(VariantAnnotator.ImpactLevel.High),
                "start_lost impact is HIGH in Constants.pm");
        });
    }

    // S1 — NAA (codon6 @115-117) ambiguous reference codon: not synonymous even if alt equals ref translation.
    // Src: VariationEffect.pm synonymous predicate excludes peptides containing X.
    [Test]
    public void PredictFunctionalImpact_AmbiguousCodon_NotSynonymous()
    {
        // base3 A>A would be a no-op; use base3 A>G on an N-containing codon: NAG still translates to X.
        var fi = Predict(117, "A", "G", VariantAnnotator.VariantType.SNV);

        Assert.That(fi.Consequence, Is.Not.EqualTo(VariantAnnotator.ConsequenceType.SynonymousVariant),
            "a codon containing N translates to X and is excluded from synonymous_variant");
    }

    #endregion

    #region PredictFunctionalImpact — coding indels (length rule)

    // M6 — coding deletion ref "AC" alt "A" (Δ-1, not a multiple of 3): frameshift, HIGH. Src: VariationEffect.pm
    [Test]
    public void PredictFunctionalImpact_CodingDeletion_ReturnsFrameshiftHigh()
    {
        var variant = new VariantAnnotator.Variant("chr1", 106, "AC", "A", VariantAnnotator.VariantType.Deletion);
        var fi = VariantAnnotator.PredictFunctionalImpact(variant, CodingTranscript(), RefWindow, RefStart);

        Assert.Multiple(() =>
        {
            Assert.That(fi.Consequence, Is.EqualTo(VariantAnnotator.ConsequenceType.FrameshiftVariant),
                "a 1-base coding deletion shifts the reading frame => frameshift_variant");
            Assert.That(fi.Impact, Is.EqualTo(VariantAnnotator.ImpactLevel.High),
                "frameshift_variant impact is HIGH in Constants.pm");
        });
    }

    // M7 — coding insertion ref "A" alt "ATTT" (Δ+3, multiple of 3): inframe_insertion, MODERATE. Src: VariationEffect.pm
    [Test]
    public void PredictFunctionalImpact_InframeInsertion_ReturnsInframeInsertionModerate()
    {
        var variant = new VariantAnnotator.Variant("chr1", 106, "A", "ATTT", VariantAnnotator.VariantType.Insertion);
        var fi = VariantAnnotator.PredictFunctionalImpact(variant, CodingTranscript(), RefWindow, RefStart);

        Assert.Multiple(() =>
        {
            Assert.That(fi.Consequence, Is.EqualTo(VariantAnnotator.ConsequenceType.InframeInsertion),
                "a 3-base coding insertion preserves the frame and lengthens the protein => inframe_insertion");
            Assert.That(fi.Impact, Is.EqualTo(VariantAnnotator.ImpactLevel.Moderate),
                "inframe_insertion impact is MODERATE in Constants.pm");
        });
    }

    // M8 — coding deletion ref "ATTT" alt "A" (Δ-3, multiple of 3): inframe_deletion, MODERATE. Src: VariationEffect.pm
    [Test]
    public void PredictFunctionalImpact_InframeDeletion_ReturnsInframeDeletionModerate()
    {
        var variant = new VariantAnnotator.Variant("chr1", 106, "ATTT", "A", VariantAnnotator.VariantType.Deletion);
        var fi = VariantAnnotator.PredictFunctionalImpact(variant, CodingTranscript(), RefWindow, RefStart);

        Assert.Multiple(() =>
        {
            Assert.That(fi.Consequence, Is.EqualTo(VariantAnnotator.ConsequenceType.InframeDeletion),
                "a 3-base coding deletion preserves the frame and shortens the protein => inframe_deletion");
            Assert.That(fi.Impact, Is.EqualTo(VariantAnnotator.ImpactLevel.Moderate),
                "inframe_deletion impact is MODERATE in Constants.pm");
        });
    }

    #endregion

    #region GetImpactLevel — Constants.pm impact mapping

    // M9 — HIGH-impact terms per Constants.pm.
    [Test]
    public void GetImpactLevel_HighTerms_ReturnsHigh()
    {
        var high = new[]
        {
            VariantAnnotator.ConsequenceType.StopGained,
            VariantAnnotator.ConsequenceType.FrameshiftVariant,
            VariantAnnotator.ConsequenceType.StopLost,
            VariantAnnotator.ConsequenceType.StartLost,
            VariantAnnotator.ConsequenceType.SpliceDonorVariant,
            VariantAnnotator.ConsequenceType.SpliceAcceptorVariant,
        };

        Assert.Multiple(() =>
        {
            foreach (var c in high)
                Assert.That(VariantAnnotator.GetImpactLevel(c), Is.EqualTo(VariantAnnotator.ImpactLevel.High),
                    $"{c} is HIGH impact in Constants.pm");
        });
    }

    // M10 — MODERATE-impact terms per Constants.pm.
    [Test]
    public void GetImpactLevel_ModerateTerms_ReturnsModerate()
    {
        var moderate = new[]
        {
            VariantAnnotator.ConsequenceType.MissenseVariant,
            VariantAnnotator.ConsequenceType.InframeInsertion,
            VariantAnnotator.ConsequenceType.InframeDeletion,
        };

        Assert.Multiple(() =>
        {
            foreach (var c in moderate)
                Assert.That(VariantAnnotator.GetImpactLevel(c), Is.EqualTo(VariantAnnotator.ImpactLevel.Moderate),
                    $"{c} is MODERATE impact in Constants.pm");
        });
    }

    // M11 — LOW-impact terms per Constants.pm.
    [Test]
    public void GetImpactLevel_LowTerms_ReturnsLow()
    {
        Assert.Multiple(() =>
        {
            Assert.That(VariantAnnotator.GetImpactLevel(VariantAnnotator.ConsequenceType.SynonymousVariant),
                Is.EqualTo(VariantAnnotator.ImpactLevel.Low), "synonymous_variant is LOW in Constants.pm");
            Assert.That(VariantAnnotator.GetImpactLevel(VariantAnnotator.ConsequenceType.SpliceRegionVariant),
                Is.EqualTo(VariantAnnotator.ImpactLevel.Low), "splice_region_variant is LOW in Constants.pm");
        });
    }

    // M12 — MODIFIER-impact terms per Constants.pm.
    [Test]
    public void GetImpactLevel_ModifierTerms_ReturnsModifier()
    {
        var modifier = new[]
        {
            VariantAnnotator.ConsequenceType.IntronVariant,
            VariantAnnotator.ConsequenceType.FivePrimeUtrVariant,
            VariantAnnotator.ConsequenceType.ThreePrimeUtrVariant,
            VariantAnnotator.ConsequenceType.UpstreamGeneVariant,
            VariantAnnotator.ConsequenceType.DownstreamGeneVariant,
            VariantAnnotator.ConsequenceType.IntergenicVariant,
        };

        Assert.Multiple(() =>
        {
            foreach (var c in modifier)
                Assert.That(VariantAnnotator.GetImpactLevel(c), Is.EqualTo(VariantAnnotator.ImpactLevel.Modifier),
                    $"{c} is MODIFIER impact in Constants.pm");
        });
    }

    #endregion

    #region Annotate — most-severe selection

    // M13 — variant overlapping two transcripts: a coding missense (rank 13) is more severe
    // than an intron_variant (rank 28); Annotate returns the missense. The same genomic
    // position is exonic/coding in transcript A and intronic in transcript B.
    // Src: McLaren 2016 most-severe reporting; Constants.pm ranks (missense=13, intron=28).
    [Test]
    public void Annotate_TwoTranscripts_ReturnsMostSevere()
    {
        // Transcript A: single exon (90..130) wrapping CDS (100..120). Variant @107 is the
        // 2nd base of codon3 (GAA @106-108); A>T -> GTA (Glu->Val) = missense (rank 13).
        var transcriptA = new VariantAnnotator.Transcript(
            "ENST_A", "ENSG", "GENE", "chr1", 90, 130, '+',
            new List<(int, int)> { (90, 130) }, new List<(int, int)> { (100, 120) }, 100, 120);

        // Transcript B: two exons (90..95) and (120..130) with a wide intron 96..119. Genomic
        // 107 lies deep in the intron (>8 bp from both exon boundaries, clear of splice regions)
        // => intron_variant (rank 28).
        var transcriptB = new VariantAnnotator.Transcript(
            "ENST_B", "ENSG", "GENE", "chr1", 90, 130, '+',
            new List<(int, int)> { (90, 95), (120, 130) },
            new List<(int, int)> { (120, 125) }, 120, 125);

        // referenceSequence genomic 100..120 = "ATGCAAGAATTATAANAAGGG"; codon3 @106-108 = "GAA".
        var variant = new VariantAnnotator.Variant("chr1", 107, "A", "T", VariantAnnotator.VariantType.SNV);

        var result = VariantAnnotator.Annotate(
            new[] { variant }, new[] { transcriptA, transcriptB }, RefWindow, RefStart).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(result, Has.Count.EqualTo(1), "Annotate returns one annotation per variant");
            Assert.That(result[0].Consequence, Is.EqualTo(VariantAnnotator.ConsequenceType.MissenseVariant),
                "missense (rank 13) is more severe than intron_variant (rank 28); Annotate returns missense");
            Assert.That(result[0].TranscriptId, Is.EqualTo("ENST_A"),
                "the most severe consequence comes from the coding transcript A");
        });
    }

    // S2 — variant far from any transcript: intergenic, MODIFIER. Src: Constants.pm intergenic_variant.
    [Test]
    public void Annotate_VariantFarFromTranscripts_ReturnsIntergenic()
    {
        var variant = new VariantAnnotator.Variant("chr1", 5_000_000, "A", "G", VariantAnnotator.VariantType.SNV);
        var transcript = CodingTranscript();

        var result = VariantAnnotator.Annotate(new[] { variant }, new[] { transcript }).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(result, Has.Count.EqualTo(1), "one annotation per variant");
            Assert.That(result[0].Consequence, Is.EqualTo(VariantAnnotator.ConsequenceType.IntergenicVariant),
                "a variant with no nearby transcript is intergenic_variant");
            Assert.That(result[0].Impact, Is.EqualTo(VariantAnnotator.ImpactLevel.Modifier),
                "intergenic_variant impact is MODIFIER in Constants.pm");
        });
    }

    #endregion

    #region VCF formatting — culture invariance

    // S3 — VCF INFO numeric fields use a '.' decimal separator regardless of current culture.
    // Src: VCF 4.x fields are locale-independent ASCII.
    [Test]
    public void FormatAsVcfInfo_NonInvariantCulture_UsesDotDecimal()
    {
        var original = Thread.CurrentThread.CurrentCulture;
        try
        {
            // de-DE uses ',' as the decimal separator; the output must still use '.'.
            Thread.CurrentThread.CurrentCulture = CultureInfo.GetCultureInfo("de-DE");

            var variant = new VariantAnnotator.Variant("chr1", 1000, "A", "G", VariantAnnotator.VariantType.SNV);
            var annotation = new VariantAnnotator.VariantAnnotation(
                variant, "ENST00001", "ENSG00001", "BRCA1",
                VariantAnnotator.ConsequenceType.MissenseVariant,
                VariantAnnotator.ImpactLevel.Moderate,
                "c.100A>G", "p.R100W", 100, 100,
                0.01, 0.95, null, null, null);

            string info = VariantAnnotator.FormatAsVcfInfo(annotation);

            Assert.Multiple(() =>
            {
                Assert.That(info, Does.Contain("SIFT=0.010"),
                    "SIFT must format with an invariant '.' decimal separator even under de-DE");
                Assert.That(info, Does.Contain("POLYPHEN=0.950"),
                    "POLYPHEN must format with an invariant '.' decimal separator even under de-DE");
            });
        }
        finally
        {
            Thread.CurrentThread.CurrentCulture = original;
        }
    }

    #endregion

    #region Validation / edge cases

    // C1 — Annotate(null variants) throws ArgumentNullException.
    [Test]
    public void Annotate_NullVariants_Throws()
    {
        Assert.Throws<ArgumentNullException>(
            () => VariantAnnotator.Annotate(null!, new List<VariantAnnotator.Transcript>()).ToList());
    }

    // C1b — Annotate(null annotations) throws ArgumentNullException.
    [Test]
    public void Annotate_NullAnnotations_Throws()
    {
        Assert.Throws<ArgumentNullException>(
            () => VariantAnnotator.Annotate(new List<VariantAnnotator.Variant>(), null!).ToList());
    }

    // C2 — PredictFunctionalImpact with empty reference sequence throws ArgumentException.
    [Test]
    public void PredictFunctionalImpact_EmptyReferenceSequence_Throws()
    {
        var variant = new VariantAnnotator.Variant("chr1", 107, "A", "T", VariantAnnotator.VariantType.SNV);
        Assert.Throws<ArgumentException>(
            () => VariantAnnotator.PredictFunctionalImpact(variant, CodingTranscript(), "", RefStart));
    }

    #endregion
}
