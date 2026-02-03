using System;
using System.Collections.Generic;
using System.Linq;

namespace Seqeron.Genomics.Annotation;

/// <summary>
/// Provides variant annotation and effect prediction algorithms.
/// Includes VEP-like annotation, pathogenicity prediction, and functional impact assessment.
/// </summary>
public static class VariantAnnotator
{
    #region Enums

    /// <summary>
    /// Types of genetic variants.
    /// </summary>
    public enum VariantType
    {
        SNV,           // Single nucleotide variant
        Insertion,
        Deletion,
        MNV,           // Multi-nucleotide variant
        Complex,
        Indel
    }

    /// <summary>
    /// Variant consequence types (VEP-like).
    /// </summary>
    public enum ConsequenceType
    {
        // High impact
        TranscriptAblation,
        SpliceAcceptorVariant,
        SpliceDonorVariant,
        StopGained,
        FrameshiftVariant,
        StopLost,
        StartLost,
        TranscriptAmplification,

        // Moderate impact
        InframeInsertion,
        InframeDeletion,
        MissenseVariant,
        ProteinAlteringVariant,

        // Low impact
        SpliceRegionVariant,
        IncompleteTerminalCodonVariant,
        StartRetained,
        StopRetained,
        SynonymousVariant,

        // Modifier
        CodingSequenceVariant,
        MatureMirnaVariant,
        FivePrimeUtrVariant,
        ThreePrimeUtrVariant,
        NonCodingTranscriptExonVariant,
        IntronVariant,
        NmdTranscriptVariant,
        NonCodingTranscriptVariant,
        UpstreamGeneVariant,
        DownstreamGeneVariant,
        TfbsAblation,
        TfbsAmplification,
        TfBindingSiteVariant,
        RegulatoryRegionAblation,
        RegulatoryRegionAmplification,
        FeatureElongation,
        RegulatoryRegionVariant,
        FeatureTruncation,
        IntergenicVariant
    }

    /// <summary>
    /// Impact severity levels.
    /// </summary>
    public enum ImpactLevel
    {
        High,
        Moderate,
        Low,
        Modifier
    }

    /// <summary>
    /// Pathogenicity classification.
    /// </summary>
    public enum PathogenicityClass
    {
        Pathogenic,
        LikelyPathogenic,
        UncertainSignificance,
        LikelyBenign,
        Benign
    }

    #endregion

    #region Records

    /// <summary>
    /// Represents a genetic variant.
    /// </summary>
    public readonly record struct Variant(
        string Chromosome,
        int Position,
        string Reference,
        string Alternate,
        VariantType Type,
        double? Quality = null,
        string? Id = null);

    /// <summary>
    /// Transcript information for annotation.
    /// </summary>
    public readonly record struct Transcript(
        string TranscriptId,
        string GeneId,
        string GeneName,
        string Chromosome,
        int Start,
        int End,
        char Strand,
        IReadOnlyList<(int Start, int End)> Exons,
        IReadOnlyList<(int Start, int End)> CodingExons,
        int? CdsStart,
        int? CdsEnd);

    /// <summary>
    /// Variant annotation result.
    /// </summary>
    public readonly record struct VariantAnnotation(
        Variant Variant,
        string TranscriptId,
        string GeneId,
        string GeneName,
        ConsequenceType Consequence,
        ImpactLevel Impact,
        string? CodonChange,
        string? AminoAcidChange,
        int? ProteinPosition,
        int? CdsPosition,
        double? SiftScore,
        double? PolyphenScore,
        double? CaddScore,
        string? ExistingVariation,
        IReadOnlyDictionary<string, double>? PopulationFrequencies);

    /// <summary>
    /// Pathogenicity prediction result.
    /// </summary>
    public readonly record struct PathogenicityPrediction(
        Variant Variant,
        PathogenicityClass Classification,
        double ConfidenceScore,
        IReadOnlyList<string> EvidenceCriteria,
        double? ClinicalSignificance,
        bool IsActionable);

    /// <summary>
    /// Conservation score result.
    /// </summary>
    public readonly record struct ConservationScore(
        string Chromosome,
        int Position,
        double PhyloP,
        double PhastCons,
        double Gerp,
        int ConservedSpeciesCount);

    /// <summary>
    /// Regulatory annotation.
    /// </summary>
    public readonly record struct RegulatoryAnnotation(
        string Chromosome,
        int Start,
        int End,
        string FeatureType,
        string? CellType,
        double? Score,
        IReadOnlyList<string> TranscriptionFactors);

    #endregion

    #region Variant Classification

    /// <summary>
    /// Classifies variant type from reference and alternate alleles.
    /// </summary>
    public static VariantType ClassifyVariant(string reference, string alternate)
    {
        if (string.IsNullOrEmpty(reference) || string.IsNullOrEmpty(alternate))
            return VariantType.Complex;

        reference = reference.ToUpperInvariant();
        alternate = alternate.ToUpperInvariant();

        if (reference.Length == 1 && alternate.Length == 1)
            return VariantType.SNV;

        if (reference.Length == alternate.Length)
            return reference.Length > 1 ? VariantType.MNV : VariantType.SNV;

        if (reference.Length > alternate.Length)
        {
            // Check if it's a pure deletion
            if (alternate.Length == 1 && reference.StartsWith(alternate))
                return VariantType.Deletion;
            return VariantType.Indel;
        }

        // alternate.Length > reference.Length
        if (reference.Length == 1 && alternate.StartsWith(reference))
            return VariantType.Insertion;

        return VariantType.Indel;
    }

    /// <summary>
    /// Normalizes variant representation (left-aligns and trims).
    /// </summary>
    public static Variant NormalizeVariant(
        string chromosome,
        int position,
        string reference,
        string alternate,
        string? referenceSequence = null)
    {
        reference = reference.ToUpperInvariant();
        alternate = alternate.ToUpperInvariant();

        // Trim common suffix
        while (reference.Length > 1 && alternate.Length > 1 &&
               reference[^1] == alternate[^1])
        {
            reference = reference[..^1];
            alternate = alternate[..^1];
        }

        // Trim common prefix
        while (reference.Length > 1 && alternate.Length > 1 &&
               reference[0] == alternate[0])
        {
            reference = reference[1..];
            alternate = alternate[1..];
            position++;
        }

        var type = ClassifyVariant(reference, alternate);
        return new Variant(chromosome, position, reference, alternate, type);
    }

    #endregion

    #region Consequence Prediction

    /// <summary>
    /// Predicts variant consequences on transcripts.
    /// </summary>
    public static IEnumerable<VariantAnnotation> AnnotateVariant(
        Variant variant,
        IEnumerable<Transcript> transcripts,
        string? referenceSequence = null,
        IReadOnlyDictionary<string, double>? populationFrequencies = null)
    {
        var relevantTranscripts = transcripts
            .Where(t => t.Chromosome == variant.Chromosome &&
                       OverlapsOrNear(variant.Position, variant.Position + variant.Reference.Length - 1,
                                     t.Start, t.End, upstream: 5000, downstream: 500))
            .ToList();

        if (relevantTranscripts.Count == 0)
        {
            // Intergenic variant
            yield return new VariantAnnotation(
                variant, "", "", "",
                ConsequenceType.IntergenicVariant,
                ImpactLevel.Modifier,
                null, null, null, null,
                null, null, null, null, populationFrequencies);
            yield break;
        }

        foreach (var transcript in relevantTranscripts)
        {
            var consequence = DetermineConsequence(variant, transcript, referenceSequence);
            var impact = GetImpactLevel(consequence);

            string? codonChange = null;
            string? aaChange = null;
            int? proteinPos = null;
            int? cdsPos = null;

            if (IsCodingConsequence(consequence) && transcript.CdsStart.HasValue)
            {
                cdsPos = CalculateCdsPosition(variant.Position, transcript);
                if (cdsPos.HasValue)
                {
                    proteinPos = (cdsPos.Value - 1) / 3 + 1;
                    (codonChange, aaChange) = PredictAminoAcidChange(
                        variant, transcript, cdsPos.Value, referenceSequence);
                }
            }

            // Calculate prediction scores
            double? sift = null;
            double? polyphen = null;

            if (consequence == ConsequenceType.MissenseVariant && aaChange != null)
            {
                sift = PredictSiftScore(aaChange);
                polyphen = PredictPolyphenScore(aaChange);
            }

            yield return new VariantAnnotation(
                variant,
                transcript.TranscriptId,
                transcript.GeneId,
                transcript.GeneName,
                consequence,
                impact,
                codonChange,
                aaChange,
                proteinPos,
                cdsPos,
                sift,
                polyphen,
                null,
                null,
                populationFrequencies);
        }
    }

    /// <summary>
    /// Determines the most severe consequence for a variant on a transcript.
    /// </summary>
    private static ConsequenceType DetermineConsequence(
        Variant variant,
        Transcript transcript,
        string? referenceSequence)
    {
        int varStart = variant.Position;
        int varEnd = variant.Position + variant.Reference.Length - 1;

        // Check if upstream
        if (transcript.Strand == '+' && varEnd < transcript.Start)
        {
            return transcript.Start - varEnd <= 5000 ?
                ConsequenceType.UpstreamGeneVariant :
                ConsequenceType.IntergenicVariant;
        }
        if (transcript.Strand == '-' && varStart > transcript.End)
        {
            return varStart - transcript.End <= 5000 ?
                ConsequenceType.UpstreamGeneVariant :
                ConsequenceType.IntergenicVariant;
        }

        // Check if downstream
        if (transcript.Strand == '+' && varStart > transcript.End)
        {
            return varStart - transcript.End <= 500 ?
                ConsequenceType.DownstreamGeneVariant :
                ConsequenceType.IntergenicVariant;
        }
        if (transcript.Strand == '-' && varEnd < transcript.Start)
        {
            return transcript.Start - varEnd <= 500 ?
                ConsequenceType.DownstreamGeneVariant :
                ConsequenceType.IntergenicVariant;
        }

        // Check splice sites
        var spliceConsequence = CheckSpliceSite(variant, transcript);
        if (spliceConsequence.HasValue)
            return spliceConsequence.Value;

        // Check if in exon
        bool inExon = transcript.Exons.Any(e => varStart <= e.End && varEnd >= e.Start);

        if (!inExon)
            return ConsequenceType.IntronVariant;

        // Check if in coding region
        if (!transcript.CdsStart.HasValue || !transcript.CdsEnd.HasValue)
            return ConsequenceType.NonCodingTranscriptExonVariant;

        bool inCoding = transcript.CodingExons.Any(e => varStart <= e.End && varEnd >= e.Start);

        if (!inCoding)
        {
            // In UTR
            if (transcript.Strand == '+')
            {
                return varEnd < transcript.CdsStart.Value ?
                    ConsequenceType.FivePrimeUtrVariant :
                    ConsequenceType.ThreePrimeUtrVariant;
            }
            else
            {
                return varStart > transcript.CdsEnd.Value ?
                    ConsequenceType.FivePrimeUtrVariant :
                    ConsequenceType.ThreePrimeUtrVariant;
            }
        }

        // Coding variant - determine specific consequence
        return DetermineCodingConsequence(variant, transcript, referenceSequence);
    }

    /// <summary>
    /// Checks for splice site variants.
    /// </summary>
    private static ConsequenceType? CheckSpliceSite(Variant variant, Transcript transcript)
    {
        int varStart = variant.Position;
        int varEnd = variant.Position + variant.Reference.Length - 1;

        foreach (var exon in transcript.Exons)
        {
            // Splice acceptor (2bp at start of exon)
            if (transcript.Strand == '+')
            {
                if (varStart <= exon.Start && varEnd >= exon.Start - 1)
                    return ConsequenceType.SpliceAcceptorVariant;
                if (varStart <= exon.End + 2 && varEnd >= exon.End + 1)
                    return ConsequenceType.SpliceDonorVariant;
            }
            else
            {
                if (varStart <= exon.End + 2 && varEnd >= exon.End + 1)
                    return ConsequenceType.SpliceAcceptorVariant;
                if (varStart <= exon.Start && varEnd >= exon.Start - 1)
                    return ConsequenceType.SpliceDonorVariant;
            }

            // Splice region (3-8 bp into intron)
            int donorStart = transcript.Strand == '+' ? exon.End + 3 : exon.Start - 8;
            int donorEnd = transcript.Strand == '+' ? exon.End + 8 : exon.Start - 3;
            int acceptorStart = transcript.Strand == '+' ? exon.Start - 8 : exon.End + 3;
            int acceptorEnd = transcript.Strand == '+' ? exon.Start - 3 : exon.End + 8;

            if ((varStart <= donorEnd && varEnd >= donorStart) ||
                (varStart <= acceptorEnd && varEnd >= acceptorStart))
            {
                return ConsequenceType.SpliceRegionVariant;
            }
        }

        return null;
    }

    /// <summary>
    /// Determines consequence for coding variants.
    /// </summary>
    private static ConsequenceType DetermineCodingConsequence(
        Variant variant,
        Transcript transcript,
        string? referenceSequence)
    {
        // Frameshift check for indels
        if (variant.Type == VariantType.Insertion || variant.Type == VariantType.Deletion ||
            variant.Type == VariantType.Indel)
        {
            int lengthDiff = variant.Alternate.Length - variant.Reference.Length;
            if (lengthDiff % 3 != 0)
                return ConsequenceType.FrameshiftVariant;

            return lengthDiff > 0 ?
                ConsequenceType.InframeInsertion :
                ConsequenceType.InframeDeletion;
        }

        // SNV/MNV - check for stop codon, start codon, synonymous, missense
        // Simplified: would need actual codon lookup
        if (IsStartCodonPosition(variant.Position, transcript))
            return ConsequenceType.StartLost;

        // Default to missense for coding SNVs
        return ConsequenceType.MissenseVariant;
    }

    /// <summary>
    /// Checks if position is at start codon.
    /// </summary>
    private static bool IsStartCodonPosition(int position, Transcript transcript)
    {
        if (!transcript.CdsStart.HasValue) return false;

        if (transcript.Strand == '+')
            return position >= transcript.CdsStart.Value && position <= transcript.CdsStart.Value + 2;
        else
            return position >= transcript.CdsEnd!.Value - 2 && position <= transcript.CdsEnd.Value;
    }

    #endregion

    #region Impact and Scores

    /// <summary>
    /// Gets impact level for a consequence type.
    /// </summary>
    public static ImpactLevel GetImpactLevel(ConsequenceType consequence)
    {
        return consequence switch
        {
            ConsequenceType.TranscriptAblation => ImpactLevel.High,
            ConsequenceType.SpliceAcceptorVariant => ImpactLevel.High,
            ConsequenceType.SpliceDonorVariant => ImpactLevel.High,
            ConsequenceType.StopGained => ImpactLevel.High,
            ConsequenceType.FrameshiftVariant => ImpactLevel.High,
            ConsequenceType.StopLost => ImpactLevel.High,
            ConsequenceType.StartLost => ImpactLevel.High,
            ConsequenceType.TranscriptAmplification => ImpactLevel.High,

            ConsequenceType.InframeInsertion => ImpactLevel.Moderate,
            ConsequenceType.InframeDeletion => ImpactLevel.Moderate,
            ConsequenceType.MissenseVariant => ImpactLevel.Moderate,
            ConsequenceType.ProteinAlteringVariant => ImpactLevel.Moderate,

            ConsequenceType.SpliceRegionVariant => ImpactLevel.Low,
            ConsequenceType.IncompleteTerminalCodonVariant => ImpactLevel.Low,
            ConsequenceType.StartRetained => ImpactLevel.Low,
            ConsequenceType.StopRetained => ImpactLevel.Low,
            ConsequenceType.SynonymousVariant => ImpactLevel.Low,

            _ => ImpactLevel.Modifier
        };
    }

    /// <summary>
    /// Checks if consequence affects coding sequence.
    /// </summary>
    private static bool IsCodingConsequence(ConsequenceType consequence)
    {
        return consequence switch
        {
            ConsequenceType.StopGained => true,
            ConsequenceType.FrameshiftVariant => true,
            ConsequenceType.StopLost => true,
            ConsequenceType.StartLost => true,
            ConsequenceType.InframeInsertion => true,
            ConsequenceType.InframeDeletion => true,
            ConsequenceType.MissenseVariant => true,
            ConsequenceType.SynonymousVariant => true,
            ConsequenceType.ProteinAlteringVariant => true,
            _ => false
        };
    }

    /// <summary>
    /// Calculates CDS position for a variant.
    /// </summary>
    private static int? CalculateCdsPosition(int genomicPosition, Transcript transcript)
    {
        if (!transcript.CdsStart.HasValue || !transcript.CdsEnd.HasValue)
            return null;

        int cdsPos = 0;
        var codingExons = transcript.CodingExons.OrderBy(e =>
            transcript.Strand == '+' ? e.Start : -e.End).ToList();

        foreach (var exon in codingExons)
        {
            if (genomicPosition >= exon.Start && genomicPosition <= exon.End)
            {
                if (transcript.Strand == '+')
                    return cdsPos + (genomicPosition - exon.Start) + 1;
                else
                    return cdsPos + (exon.End - genomicPosition) + 1;
            }
            cdsPos += exon.End - exon.Start + 1;
        }

        return null;
    }

    /// <summary>
    /// Predicts amino acid change from variant.
    /// </summary>
    private static (string? CodonChange, string? AaChange) PredictAminoAcidChange(
        Variant variant,
        Transcript transcript,
        int cdsPosition,
        string? referenceSequence)
    {
        // Simplified amino acid prediction
        // In real implementation, would use actual codon table and sequence lookup

        if (variant.Type != VariantType.SNV)
            return (null, null);

        // Get codon position
        int codonNumber = (cdsPosition - 1) / 3 + 1;
        int positionInCodon = (cdsPosition - 1) % 3;

        // Placeholder for actual codon lookup
        // Would need reference sequence to get actual codons
        string refAa = "X";
        string altAa = "Y";

        if (refAa == altAa)
            return ($"c.{cdsPosition}{variant.Reference}>{variant.Alternate}",
                    $"p.{refAa}{codonNumber}=");

        return ($"c.{cdsPosition}{variant.Reference}>{variant.Alternate}",
                $"p.{refAa}{codonNumber}{altAa}");
    }

    /// <summary>
    /// Predicts SIFT-like score (deleterious < 0.05).
    /// </summary>
    private static double PredictSiftScore(string aminoAcidChange)
    {
        if (string.IsNullOrEmpty(aminoAcidChange))
            return 1.0;

        // Simplified scoring based on amino acid properties
        // Real SIFT uses sequence alignment and conservation

        // Parse amino acid change
        if (!TryParseAaChange(aminoAcidChange, out char refAa, out char altAa))
            return 0.5;

        // Check for stop gained
        if (altAa == '*')
            return 0.0;

        // Check for synonymous
        if (refAa == altAa)
            return 1.0;

        // Score based on biochemical similarity
        double score = GetAminoAcidSimilarity(refAa, altAa);
        return score;
    }

    /// <summary>
    /// Predicts PolyPhen-like score (probably damaging > 0.908).
    /// </summary>
    private static double PredictPolyphenScore(string aminoAcidChange)
    {
        if (string.IsNullOrEmpty(aminoAcidChange))
            return 0.0;

        if (!TryParseAaChange(aminoAcidChange, out char refAa, out char altAa))
            return 0.5;

        if (altAa == '*')
            return 1.0;

        if (refAa == altAa)
            return 0.0;

        // Inverse of SIFT-like score
        double similarity = GetAminoAcidSimilarity(refAa, altAa);
        return 1.0 - similarity;
    }

    /// <summary>
    /// Parses amino acid change notation.
    /// </summary>
    private static bool TryParseAaChange(string aaChange, out char refAa, out char altAa)
    {
        refAa = 'X';
        altAa = 'X';

        if (aaChange.StartsWith("p.") && aaChange.Length >= 5)
        {
            refAa = aaChange[2];
            altAa = aaChange[^1];
            return true;
        }

        return false;
    }

    /// <summary>
    /// Calculates biochemical similarity between amino acids.
    /// </summary>
    private static double GetAminoAcidSimilarity(char aa1, char aa2)
    {
        // BLOSUM62-inspired scoring
        var aromatic = new HashSet<char> { 'F', 'W', 'Y' };
        var aliphatic = new HashSet<char> { 'A', 'V', 'L', 'I', 'M' };
        var polar = new HashSet<char> { 'S', 'T', 'N', 'Q' };
        var positive = new HashSet<char> { 'K', 'R', 'H' };
        var negative = new HashSet<char> { 'D', 'E' };
        var special = new HashSet<char> { 'G', 'P', 'C' };

        if (aa1 == aa2) return 1.0;

        // Same biochemical group
        if ((aromatic.Contains(aa1) && aromatic.Contains(aa2)) ||
            (aliphatic.Contains(aa1) && aliphatic.Contains(aa2)) ||
            (polar.Contains(aa1) && polar.Contains(aa2)) ||
            (positive.Contains(aa1) && positive.Contains(aa2)) ||
            (negative.Contains(aa1) && negative.Contains(aa2)))
        {
            return 0.7;
        }

        // Glycine and proline are special
        if (special.Contains(aa1) || special.Contains(aa2))
            return 0.1;

        // Charge change
        if ((positive.Contains(aa1) && negative.Contains(aa2)) ||
            (negative.Contains(aa1) && positive.Contains(aa2)))
        {
            return 0.05;
        }

        return 0.3;
    }

    #endregion

    #region Pathogenicity Prediction

    /// <summary>
    /// Predicts variant pathogenicity using ACMG-like criteria.
    /// </summary>
    public static PathogenicityPrediction PredictPathogenicity(
        VariantAnnotation annotation,
        double? populationFrequency = null,
        double? conservationScore = null,
        bool inClinvar = false,
        string? clinvarSignificance = null,
        IEnumerable<string>? functionalEvidence = null)
    {
        var evidence = new List<string>();
        int pathogenicPoints = 0;
        int benignPoints = 0;

        // PVS1: Null variant in gene where LOF is known mechanism
        if (annotation.Impact == ImpactLevel.High)
        {
            pathogenicPoints += 8;
            evidence.Add("PVS1: Null variant (high impact)");
        }

        // PS1/PM5: Same/similar amino acid change as established pathogenic
        // (Simplified - would need database lookup)

        // PM1: Located in mutational hot spot or critical domain
        // (Would need domain annotation)

        // PM2: Absent from controls
        if (populationFrequency.HasValue)
        {
            if (populationFrequency < 0.0001)
            {
                pathogenicPoints += 2;
                evidence.Add("PM2: Extremely rare (AF < 0.01%)");
            }
            else if (populationFrequency < 0.001)
            {
                pathogenicPoints += 1;
                evidence.Add("PM2: Rare (AF < 0.1%)");
            }
            else if (populationFrequency > 0.05)
            {
                benignPoints += 4;
                evidence.Add("BA1: MAF > 5%");
            }
            else if (populationFrequency > 0.01)
            {
                benignPoints += 2;
                evidence.Add("BS1: MAF > 1%");
            }
        }

        // PP2/BP1: Missense in gene with constraint
        if (annotation.Consequence == ConsequenceType.MissenseVariant)
        {
            if (annotation.SiftScore.HasValue && annotation.SiftScore < 0.05)
            {
                pathogenicPoints += 1;
                evidence.Add("PP3: SIFT deleterious");
            }
            if (annotation.PolyphenScore.HasValue && annotation.PolyphenScore > 0.908)
            {
                pathogenicPoints += 1;
                evidence.Add("PP3: PolyPhen probably damaging");
            }
        }

        // BP4: Computational evidence suggests no impact
        if (annotation.SiftScore.HasValue && annotation.SiftScore > 0.95 &&
            annotation.PolyphenScore.HasValue && annotation.PolyphenScore < 0.1)
        {
            benignPoints += 1;
            evidence.Add("BP4: Computational evidence benign");
        }

        // BP7: Synonymous with no splicing impact
        if (annotation.Consequence == ConsequenceType.SynonymousVariant)
        {
            benignPoints += 1;
            evidence.Add("BP7: Synonymous variant");
        }

        // Conservation score
        if (conservationScore.HasValue)
        {
            if (conservationScore > 4.0)
            {
                pathogenicPoints += 1;
                evidence.Add("PP3: Highly conserved position");
            }
            else if (conservationScore < -2.0)
            {
                benignPoints += 1;
                evidence.Add("BP4: Non-conserved position");
            }
        }

        // ClinVar evidence
        if (inClinvar && clinvarSignificance != null)
        {
            if (clinvarSignificance.Contains("Pathogenic", StringComparison.OrdinalIgnoreCase))
            {
                pathogenicPoints += 4;
                evidence.Add("PP5: ClinVar pathogenic");
            }
            else if (clinvarSignificance.Contains("Benign", StringComparison.OrdinalIgnoreCase))
            {
                benignPoints += 4;
                evidence.Add("BP6: ClinVar benign");
            }
        }

        // Functional evidence
        if (functionalEvidence != null)
        {
            foreach (var func in functionalEvidence)
            {
                if (func.Contains("LOF", StringComparison.OrdinalIgnoreCase))
                {
                    pathogenicPoints += 4;
                    evidence.Add($"PS3: {func}");
                }
            }
        }

        // Classify based on points
        var classification = ClassifyByPoints(pathogenicPoints, benignPoints);
        double confidence = CalculateConfidence(pathogenicPoints, benignPoints);
        bool actionable = classification == PathogenicityClass.Pathogenic ||
                         classification == PathogenicityClass.LikelyPathogenic;

        return new PathogenicityPrediction(
            annotation.Variant,
            classification,
            confidence,
            evidence,
            pathogenicPoints - benignPoints,
            actionable);
    }

    /// <summary>
    /// Classifies pathogenicity based on evidence points.
    /// </summary>
    private static PathogenicityClass ClassifyByPoints(int pathogenic, int benign)
    {
        int net = pathogenic - benign;

        if (net >= 10 || pathogenic >= 10)
            return PathogenicityClass.Pathogenic;
        if (net >= 6)
            return PathogenicityClass.LikelyPathogenic;
        if (net <= -6)
            return PathogenicityClass.Benign;
        if (net <= -2)
            return PathogenicityClass.LikelyBenign;

        return PathogenicityClass.UncertainSignificance;
    }

    /// <summary>
    /// Calculates confidence score for classification.
    /// </summary>
    private static double CalculateConfidence(int pathogenic, int benign)
    {
        int total = pathogenic + benign;
        if (total == 0) return 0.5;

        double ratio = Math.Abs(pathogenic - benign) / (double)total;
        return Math.Min(0.99, 0.5 + ratio * 0.5);
    }

    #endregion

    #region Conservation Analysis

    /// <summary>
    /// Calculates conservation scores for positions.
    /// </summary>
    public static IEnumerable<ConservationScore> CalculateConservation(
        IEnumerable<(string Chromosome, int Position, IReadOnlyList<char> SpeciesAlleles)> alignedPositions)
    {
        foreach (var (chromosome, position, alleles) in alignedPositions)
        {
            if (alleles.Count == 0)
            {
                yield return new ConservationScore(chromosome, position, 0, 0, 0, 0);
                continue;
            }

            // Calculate conservation metrics
            char reference = alleles[0];
            int conserved = alleles.Count(a => a == reference);
            int total = alleles.Count;

            double conservationFraction = conserved / (double)total;

            // PhyloP-like score (-14 to +6 scale)
            double phyloP = (conservationFraction - 0.5) * 12;

            // PhastCons-like score (0 to 1)
            double phastCons = conservationFraction;

            // GERP-like score
            double expectedConserved = total * 0.25; // Random expectation
            double gerp = (conserved - expectedConserved) / Math.Max(1, total - expectedConserved);
            gerp = Math.Max(-12.36, Math.Min(6.18, gerp * 6));

            yield return new ConservationScore(
                chromosome,
                position,
                phyloP,
                phastCons,
                gerp,
                conserved);
        }
    }

    /// <summary>
    /// Identifies conserved elements from conservation scores.
    /// </summary>
    public static IEnumerable<(string Chromosome, int Start, int End, double Score)> FindConservedElements(
        IEnumerable<ConservationScore> scores,
        double threshold = 0.8,
        int minLength = 20)
    {
        var byChromosome = scores.GroupBy(s => s.Chromosome);

        foreach (var chromGroup in byChromosome)
        {
            var sortedScores = chromGroup.OrderBy(s => s.Position).ToList();
            int? elementStart = null;
            int lastPosition = -1;
            double scoreSum = 0;
            int count = 0;

            foreach (var score in sortedScores)
            {
                bool isConserved = score.PhastCons >= threshold;

                if (isConserved)
                {
                    if (!elementStart.HasValue || score.Position > lastPosition + 10)
                    {
                        // End previous element if exists
                        if (elementStart.HasValue && lastPosition - elementStart.Value + 1 >= minLength)
                        {
                            yield return (chromGroup.Key, elementStart.Value, lastPosition, scoreSum / count);
                        }

                        // Start new element
                        elementStart = score.Position;
                        scoreSum = 0;
                        count = 0;
                    }

                    scoreSum += score.PhastCons;
                    count++;
                    lastPosition = score.Position;
                }
                else if (elementStart.HasValue && score.Position > lastPosition + 10)
                {
                    // End element
                    if (lastPosition - elementStart.Value + 1 >= minLength)
                    {
                        yield return (chromGroup.Key, elementStart.Value, lastPosition, scoreSum / count);
                    }
                    elementStart = null;
                }
            }

            // Handle last element
            if (elementStart.HasValue && lastPosition - elementStart.Value + 1 >= minLength)
            {
                yield return (chromGroup.Key, elementStart.Value, lastPosition, scoreSum / count);
            }
        }
    }

    #endregion

    #region Regulatory Annotation

    /// <summary>
    /// Annotates variants with regulatory elements.
    /// </summary>
    public static IEnumerable<RegulatoryAnnotation> AnnotateRegulatoryElements(
        Variant variant,
        IEnumerable<(string Chromosome, int Start, int End, string Type, string? CellType, double? Score, IReadOnlyList<string> TFs)> regulatoryRegions)
    {
        int varStart = variant.Position;
        int varEnd = variant.Position + variant.Reference.Length - 1;

        foreach (var region in regulatoryRegions)
        {
            if (region.Chromosome != variant.Chromosome)
                continue;

            if (varEnd >= region.Start && varStart <= region.End)
            {
                yield return new RegulatoryAnnotation(
                    region.Chromosome,
                    region.Start,
                    region.End,
                    region.Type,
                    region.CellType,
                    region.Score,
                    region.TFs);
            }
        }
    }

    /// <summary>
    /// Predicts transcription factor binding disruption.
    /// </summary>
    public static IEnumerable<(string TfName, double RefScore, double AltScore, double ScoreDifference)> PredictTfBindingChange(
        Variant variant,
        IEnumerable<(string TfName, string Motif, double Threshold)> motifs,
        string referenceContext,
        int contextOffset = 20)
    {
        if (variant.Type != VariantType.SNV)
            yield break;

        foreach (var (tfName, motif, threshold) in motifs)
        {
            // Get reference and alternate contexts
            int motifStart = Math.Max(0, contextOffset - motif.Length + 1);
            int motifEnd = Math.Min(referenceContext.Length, contextOffset + motif.Length);

            // Score reference
            double refScore = ScoreMotif(referenceContext, motif, motifStart, motifEnd);

            // Create alternate sequence
            if (contextOffset < referenceContext.Length)
            {
                char[] altContext = referenceContext.ToCharArray();
                altContext[contextOffset] = variant.Alternate[0];

                double altScore = ScoreMotif(new string(altContext), motif, motifStart, motifEnd);

                if (Math.Abs(refScore - altScore) > 0.1)
                {
                    yield return (tfName, refScore, altScore, refScore - altScore);
                }
            }
        }
    }

    /// <summary>
    /// Scores a motif match.
    /// </summary>
    private static double ScoreMotif(string sequence, string motif, int start, int end)
    {
        double bestScore = 0;

        for (int i = start; i <= end - motif.Length; i++)
        {
            int matches = 0;
            for (int j = 0; j < motif.Length; j++)
            {
                char m = char.ToUpperInvariant(motif[j]);
                char s = char.ToUpperInvariant(sequence[i + j]);

                if (m == s || m == 'N')
                    matches++;
                else if (m == 'R' && (s == 'A' || s == 'G'))
                    matches++;
                else if (m == 'Y' && (s == 'C' || s == 'T'))
                    matches++;
                else if (m == 'W' && (s == 'A' || s == 'T'))
                    matches++;
                else if (m == 'S' && (s == 'G' || s == 'C'))
                    matches++;
            }

            double score = matches / (double)motif.Length;
            bestScore = Math.Max(bestScore, score);
        }

        return bestScore;
    }

    #endregion

    #region Batch Processing

    /// <summary>
    /// Annotates multiple variants efficiently.
    /// </summary>
    public static IEnumerable<IGrouping<Variant, VariantAnnotation>> AnnotateVariants(
        IEnumerable<Variant> variants,
        IEnumerable<Transcript> transcripts,
        IReadOnlyDictionary<string, double>? populationFrequencies = null)
    {
        var transcriptList = transcripts.ToList();
        var transcriptIndex = BuildTranscriptIndex(transcriptList);

        foreach (var variant in variants)
        {
            var relevantTranscripts = GetRelevantTranscripts(variant, transcriptIndex);
            var annotations = AnnotateVariant(variant, relevantTranscripts, null, populationFrequencies);

            yield return new VariantAnnotationGroup(variant, annotations.ToList());
        }
    }

    /// <summary>
    /// Builds spatial index for transcripts.
    /// </summary>
    private static Dictionary<string, List<Transcript>> BuildTranscriptIndex(List<Transcript> transcripts)
    {
        return transcripts.GroupBy(t => t.Chromosome)
                         .ToDictionary(g => g.Key, g => g.OrderBy(t => t.Start).ToList());
    }

    /// <summary>
    /// Gets transcripts near a variant.
    /// </summary>
    private static IEnumerable<Transcript> GetRelevantTranscripts(
        Variant variant,
        Dictionary<string, List<Transcript>> index)
    {
        if (!index.TryGetValue(variant.Chromosome, out var transcripts))
            yield break;

        int varPos = variant.Position;
        int searchStart = varPos - 5000;
        int searchEnd = varPos + variant.Reference.Length + 500;

        foreach (var transcript in transcripts)
        {
            if (transcript.End < searchStart)
                continue;
            if (transcript.Start > searchEnd)
                break;

            yield return transcript;
        }
    }

    /// <summary>
    /// Helper class for grouping annotations.
    /// </summary>
    private class VariantAnnotationGroup : IGrouping<Variant, VariantAnnotation>
    {
        private readonly Variant _variant;
        private readonly List<VariantAnnotation> _annotations;

        public VariantAnnotationGroup(Variant variant, List<VariantAnnotation> annotations)
        {
            _variant = variant;
            _annotations = annotations;
        }

        public Variant Key => _variant;

        public IEnumerator<VariantAnnotation> GetEnumerator() => _annotations.GetEnumerator();
        System.Collections.IEnumerator System.Collections.IEnumerable.GetEnumerator() => GetEnumerator();
    }

    #endregion

    #region VCF Integration

    /// <summary>
    /// Parses variant from VCF fields.
    /// </summary>
    public static Variant ParseVcfVariant(
        string chromosome,
        int position,
        string id,
        string reference,
        string alternate,
        double? quality = null)
    {
        var type = ClassifyVariant(reference, alternate);
        return new Variant(chromosome, position, reference, alternate, type, quality, id);
    }

    /// <summary>
    /// Formats annotation as VCF INFO field.
    /// </summary>
    public static string FormatAsVcfInfo(VariantAnnotation annotation)
    {
        var parts = new List<string>
        {
            $"GENE={annotation.GeneName}",
            $"TRANSCRIPT={annotation.TranscriptId}",
            $"CONSEQUENCE={annotation.Consequence}",
            $"IMPACT={annotation.Impact}"
        };

        if (annotation.AminoAcidChange != null)
            parts.Add($"HGVSP={annotation.AminoAcidChange}");

        if (annotation.CodonChange != null)
            parts.Add($"HGVSC={annotation.CodonChange}");

        if (annotation.SiftScore.HasValue)
            parts.Add($"SIFT={annotation.SiftScore:F3}");

        if (annotation.PolyphenScore.HasValue)
            parts.Add($"POLYPHEN={annotation.PolyphenScore:F3}");

        return string.Join(";", parts);
    }

    #endregion

    #region Helpers

    /// <summary>
    /// Checks if regions overlap with optional flanking.
    /// </summary>
    private static bool OverlapsOrNear(int start1, int end1, int start2, int end2,
                                        int upstream = 0, int downstream = 0)
    {
        return end1 >= (start2 - upstream) && start1 <= (end2 + downstream);
    }

    #endregion
}
