namespace Seqeron.Genomics.Annotation;

/// <summary>
/// Identifies variants (SNPs, insertions, deletions) between sequences.
/// Compares a query sequence against a reference to detect mutations.
/// </summary>
public static class VariantCaller
{
    // Alignment gap character. A reference-side gap marks an insertion in the query and a
    // query-side gap marks a deletion (Danecek et al. 2011, Bioinformatics 27(15):2156-2158,
    // doi:10.1093/bioinformatics/btr330 — VCF stores SNPs, insertions and deletions).
    private const char GapChar = '-';
    private const string GapAllele = "-";

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
            return Enumerable.Empty<Variant>();
        if (alignedReference.Length != alignedQuery.Length)
            throw new ArgumentException("Aligned sequences must have the same length.");
        return CallVariantsFromAlignmentCore(alignedReference, alignedQuery);
    }

    private static IEnumerable<Variant> CallVariantsFromAlignmentCore(
        string alignedReference,
        string alignedQuery)
    {
        int refPos = 0;
        int queryPos = 0;

        for (int i = 0; i < alignedReference.Length; i++)
        {
            char refBase = alignedReference[i];
            char queryBase = alignedQuery[i];

            if (refBase == GapChar && queryBase != GapChar)
            {
                // Insertion in query
                yield return new Variant(
                    Position: refPos,
                    ReferenceAllele: GapAllele,
                    AlternateAllele: queryBase.ToString(),
                    Type: VariantType.Insertion,
                    QueryPosition: queryPos);
                queryPos++;
            }
            else if (refBase != GapChar && queryBase == GapChar)
            {
                // Deletion in query
                yield return new Variant(
                    Position: refPos,
                    ReferenceAllele: refBase.ToString(),
                    AlternateAllele: GapAllele,
                    Type: VariantType.Deletion,
                    QueryPosition: queryPos);
                refPos++;
            }
            else if (refBase != GapChar && queryBase != GapChar && refBase != queryBase)
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
                // Match (or gap-gap column, which advances neither coordinate)
                if (refBase != GapChar) refPos++;
                if (queryBase != GapChar) queryPos++;
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
        return AnnotateVariantsCore(reference, query, isCodingSequence);
    }

    private static IEnumerable<AnnotatedVariant> AnnotateVariantsCore(DnaSequence reference, DnaSequence query, bool isCodingSequence)
    {
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
        return ToVcfLinesCore(variants, chromosome, sampleName);
    }

    private static IEnumerable<string> ToVcfLinesCore(IEnumerable<Variant> variants, string chromosome, string sampleName)
    {
        // VCF header
        yield return "##fileformat=VCFv4.2";
        yield return $"##source=Seqeron.Genomics";
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
/// A detected genetic variant. <see cref="Position"/> is 0-based (internal/array convention).
/// </summary>
/// <remarks>
/// The VCF specification uses 1-based POS coordinates (VCFv4.2/4.3 §1.4.1, "the reference
/// position, with the 1st base having position 1"). Use <see cref="VcfPosition"/> for the
/// VCF coordinate, or <see cref="VariantCaller.ToVcfLines"/> which already emits Position+1.
/// Source: VCF v4.3 specification, https://samtools.github.io/hts-specs/VCFv4.3.pdf .
/// </remarks>
public readonly record struct Variant(
    int Position,
    string ReferenceAllele,
    string AlternateAllele,
    VariantType Type,
    int QueryPosition)
{
    /// <summary>
    /// The variant position in VCF 1-based POS coordinates (= <see cref="Position"/> + 1).
    /// Opt-in accessor for callers that consume the record directly and want the VCF
    /// convention without re-deriving the +1 offset. VCF v4.3 §1.4.1.
    /// </summary>
    public int VcfPosition => Position + 1;
}

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
