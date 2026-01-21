using System;
using System.Collections.Generic;
using System.Linq;

namespace SuffixTree.Genomics;

/// <summary>
/// Identifies variants (SNPs, insertions, deletions) between sequences.
/// Compares a query sequence against a reference to detect mutations.
/// </summary>
public static class VariantCaller
{
    #region Variant Detection

    /// <summary>
    /// Detects all variants between a query and reference sequence.
    /// </summary>
    /// <param name="reference">Reference DNA sequence.</param>
    /// <param name="query">Query DNA sequence to compare.</param>
    /// <returns>Collection of detected variants.</returns>
    public static IEnumerable<Variant> CallVariants(DnaSequence reference, DnaSequence query)
    {
        ArgumentNullException.ThrowIfNull(reference);
        ArgumentNullException.ThrowIfNull(query);

        return CallVariantsCore(reference.Sequence, query.Sequence);
    }

    /// <summary>
    /// Detects variants from aligned sequences.
    /// </summary>
    /// <param name="alignedReference">Aligned reference sequence (may contain gaps).</param>
    /// <param name="alignedQuery">Aligned query sequence (may contain gaps).</param>
    /// <returns>Collection of detected variants.</returns>
    public static IEnumerable<Variant> CallVariantsFromAlignment(
        string alignedReference,
        string alignedQuery)
    {
        if (string.IsNullOrEmpty(alignedReference) || string.IsNullOrEmpty(alignedQuery))
            yield break;

        if (alignedReference.Length != alignedQuery.Length)
            throw new ArgumentException("Aligned sequences must have the same length.");

        int refPos = 0;
        int queryPos = 0;

        for (int i = 0; i < alignedReference.Length; i++)
        {
            char refBase = alignedReference[i];
            char queryBase = alignedQuery[i];

            if (refBase == '-' && queryBase != '-')
            {
                // Insertion in query
                yield return new Variant(
                    Position: refPos,
                    ReferenceAllele: "-",
                    AlternateAllele: queryBase.ToString(),
                    Type: VariantType.Insertion,
                    QueryPosition: queryPos);
                queryPos++;
            }
            else if (refBase != '-' && queryBase == '-')
            {
                // Deletion in query
                yield return new Variant(
                    Position: refPos,
                    ReferenceAllele: refBase.ToString(),
                    AlternateAllele: "-",
                    Type: VariantType.Deletion,
                    QueryPosition: queryPos);
                refPos++;
            }
            else if (refBase != '-' && queryBase != '-' && refBase != queryBase)
            {
                // SNP
                yield return new Variant(
                    Position: refPos,
                    ReferenceAllele: refBase.ToString(),
                    AlternateAllele: queryBase.ToString(),
                    Type: VariantType.SNP,
                    QueryPosition: queryPos);
                refPos++;
                queryPos++;
            }
            else
            {
                // Match
                if (refBase != '-') refPos++;
                if (queryBase != '-') queryPos++;
            }
        }
    }

    private static IEnumerable<Variant> CallVariantsCore(string reference, string query)
    {
        // Align sequences first
        var alignment = SequenceAligner.GlobalAlign(reference, query);

        return CallVariantsFromAlignment(alignment.AlignedSequence1, alignment.AlignedSequence2);
    }

    #endregion

    #region SNP Detection

    /// <summary>
    /// Detects only SNPs (Single Nucleotide Polymorphisms).
    /// </summary>
    public static IEnumerable<Variant> FindSnps(DnaSequence reference, DnaSequence query)
    {
        return CallVariants(reference, query).Where(v => v.Type == VariantType.SNP);
    }

    /// <summary>
    /// Detects SNPs from aligned sequences (faster, no alignment needed).
    /// </summary>
    public static IEnumerable<Variant> FindSnpsDirect(string reference, string query)
    {
        if (string.IsNullOrEmpty(reference) || string.IsNullOrEmpty(query))
            yield break;

        int minLen = Math.Min(reference.Length, query.Length);

        for (int i = 0; i < minLen; i++)
        {
            if (reference[i] != query[i])
            {
                yield return new Variant(
                    Position: i,
                    ReferenceAllele: reference[i].ToString(),
                    AlternateAllele: query[i].ToString(),
                    Type: VariantType.SNP,
                    QueryPosition: i);
            }
        }
    }

    #endregion

    #region Indel Detection

    /// <summary>
    /// Detects only insertions.
    /// </summary>
    public static IEnumerable<Variant> FindInsertions(DnaSequence reference, DnaSequence query)
    {
        return CallVariants(reference, query).Where(v => v.Type == VariantType.Insertion);
    }

    /// <summary>
    /// Detects only deletions.
    /// </summary>
    public static IEnumerable<Variant> FindDeletions(DnaSequence reference, DnaSequence query)
    {
        return CallVariants(reference, query).Where(v => v.Type == VariantType.Deletion);
    }

    /// <summary>
    /// Detects insertions and deletions (indels).
    /// </summary>
    public static IEnumerable<Variant> FindIndels(DnaSequence reference, DnaSequence query)
    {
        return CallVariants(reference, query)
            .Where(v => v.Type == VariantType.Insertion || v.Type == VariantType.Deletion);
    }

    #endregion

    #region Mutation Classification

    /// <summary>
    /// Classifies a SNP as transition or transversion.
    /// </summary>
    /// <param name="variant">The SNP variant to classify.</param>
    /// <returns>Mutation type classification.</returns>
    public static MutationType ClassifyMutation(Variant variant)
    {
        if (variant.Type != VariantType.SNP)
            return MutationType.Other;

        char refBase = char.ToUpperInvariant(variant.ReferenceAllele[0]);
        char altBase = char.ToUpperInvariant(variant.AlternateAllele[0]);

        bool refPurine = refBase is 'A' or 'G';
        bool altPurine = altBase is 'A' or 'G';

        // Transition: purine <-> purine or pyrimidine <-> pyrimidine
        // Transversion: purine <-> pyrimidine
        return refPurine == altPurine ? MutationType.Transition : MutationType.Transversion;
    }

    /// <summary>
    /// Calculates the transition/transversion ratio (Ti/Tv).
    /// </summary>
    public static double CalculateTiTvRatio(IEnumerable<Variant> variants)
    {
        ArgumentNullException.ThrowIfNull(variants);

        int transitions = 0;
        int transversions = 0;

        foreach (var variant in variants.Where(v => v.Type == VariantType.SNP))
        {
            var mutationType = ClassifyMutation(variant);
            if (mutationType == MutationType.Transition)
                transitions++;
            else if (mutationType == MutationType.Transversion)
                transversions++;
        }

        return transversions > 0 ? (double)transitions / transversions : 0;
    }

    #endregion

    #region Variant Statistics

    /// <summary>
    /// Calculates comprehensive variant statistics.
    /// </summary>
    public static VariantStatistics CalculateStatistics(
        DnaSequence reference,
        DnaSequence query)
    {
        ArgumentNullException.ThrowIfNull(reference);
        ArgumentNullException.ThrowIfNull(query);

        var variants = CallVariants(reference, query).ToList();

        int snps = variants.Count(v => v.Type == VariantType.SNP);
        int insertions = variants.Count(v => v.Type == VariantType.Insertion);
        int deletions = variants.Count(v => v.Type == VariantType.Deletion);

        double titvRatio = CalculateTiTvRatio(variants);
        double variantDensity = reference.Length > 0
            ? (double)variants.Count / reference.Length * 1000
            : 0;

        return new VariantStatistics(
            TotalVariants: variants.Count,
            Snps: snps,
            Insertions: insertions,
            Deletions: deletions,
            TiTvRatio: titvRatio,
            VariantDensity: variantDensity,
            ReferenceLength: reference.Length,
            QueryLength: query.Length);
    }

    #endregion

    #region Variant Effect Prediction

    /// <summary>
    /// Predicts the effect of a variant in a coding sequence.
    /// </summary>
    /// <param name="variant">The variant to analyze.</param>
    /// <param name="codingSequence">The complete coding sequence.</param>
    /// <param name="variantPosition">Position of variant in the coding sequence.</param>
    /// <returns>Predicted variant effect.</returns>
    public static VariantEffect PredictEffect(
        Variant variant,
        DnaSequence codingSequence,
        int variantPosition)
    {
        ArgumentNullException.ThrowIfNull(codingSequence);

        if (variant.Type != VariantType.SNP)
        {
            // Frameshift for indels in coding regions
            return variant.Type == VariantType.Insertion || variant.Type == VariantType.Deletion
                ? VariantEffect.Frameshift
                : VariantEffect.Unknown;
        }

        string seq = codingSequence.Sequence;
        if (variantPosition < 0 || variantPosition >= seq.Length)
            return VariantEffect.Unknown;

        // Find codon position
        int codonStart = (variantPosition / 3) * 3;
        if (codonStart + 3 > seq.Length)
            return VariantEffect.Unknown;

        // Get reference and mutant codons
        string refCodon = seq.Substring(codonStart, 3);
        char[] mutCodonChars = refCodon.ToCharArray();
        mutCodonChars[variantPosition - codonStart] = variant.AlternateAllele[0];
        string mutCodon = new string(mutCodonChars);

        // Translate codons
        char refAa = GeneticCode.Standard.Translate(refCodon);
        char mutAa = GeneticCode.Standard.Translate(mutCodon);

        if (refAa == mutAa)
            return VariantEffect.Synonymous;

        if (mutAa == '*')
            return VariantEffect.Nonsense;

        if (refAa == '*')
            return VariantEffect.StopLoss;

        return VariantEffect.Missense;
    }

    /// <summary>
    /// Annotates all variants in a coding sequence.
    /// </summary>
    public static IEnumerable<AnnotatedVariant> AnnotateVariants(
        DnaSequence reference,
        DnaSequence query,
        bool isCodingSequence = false)
    {
        ArgumentNullException.ThrowIfNull(reference);
        ArgumentNullException.ThrowIfNull(query);

        foreach (var variant in CallVariants(reference, query))
        {
            var effect = isCodingSequence
                ? PredictEffect(variant, reference, variant.Position)
                : VariantEffect.Unknown;

            var mutationType = variant.Type == VariantType.SNP
                ? ClassifyMutation(variant)
                : MutationType.Other;

            yield return new AnnotatedVariant(
                Variant: variant,
                Effect: effect,
                MutationType: mutationType);
        }
    }

    #endregion

    #region VCF Output

    /// <summary>
    /// Formats variants as VCF (Variant Call Format) lines.
    /// </summary>
    public static IEnumerable<string> ToVcfLines(
        IEnumerable<Variant> variants,
        string chromosome = "chr1",
        string sampleName = "SAMPLE")
    {
        ArgumentNullException.ThrowIfNull(variants);

        // VCF header
        yield return "##fileformat=VCFv4.2";
        yield return $"##source=SuffixTree.Genomics";
        yield return $"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sampleName}";

        foreach (var variant in variants)
        {
            string refAllele = variant.ReferenceAllele == "-" ? "." : variant.ReferenceAllele;
            string altAllele = variant.AlternateAllele == "-" ? "." : variant.AlternateAllele;

            yield return $"{chromosome}\t{variant.Position + 1}\t.\t{refAllele}\t{altAllele}\t.\tPASS\t.\tGT\t0/1";
        }
    }

    #endregion
}

/// <summary>
/// Type of genetic variant.
/// </summary>
public enum VariantType
{
    SNP,
    Insertion,
    Deletion,
    MNP,  // Multi-nucleotide polymorphism
    Complex
}

/// <summary>
/// Type of point mutation.
/// </summary>
public enum MutationType
{
    Transition,    // Purine <-> Purine or Pyrimidine <-> Pyrimidine
    Transversion,  // Purine <-> Pyrimidine
    Other
}

/// <summary>
/// Effect of a variant on protein.
/// </summary>
public enum VariantEffect
{
    Synonymous,    // No amino acid change
    Missense,      // Amino acid change
    Nonsense,      // Creates stop codon
    StopLoss,      // Removes stop codon
    Frameshift,    // Insertion/deletion causes frameshift
    Unknown
}

/// <summary>
/// A detected genetic variant.
/// </summary>
public readonly record struct Variant(
    int Position,
    string ReferenceAllele,
    string AlternateAllele,
    VariantType Type,
    int QueryPosition);

/// <summary>
/// A variant with functional annotation.
/// </summary>
public readonly record struct AnnotatedVariant(
    Variant Variant,
    VariantEffect Effect,
    MutationType MutationType);

/// <summary>
/// Statistics about variants between two sequences.
/// </summary>
public readonly record struct VariantStatistics(
    int TotalVariants,
    int Snps,
    int Insertions,
    int Deletions,
    double TiTvRatio,
    double VariantDensity,
    int ReferenceLength,
    int QueryLength);
