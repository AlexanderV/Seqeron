namespace Seqeron.Genomics.Oncology;

public static partial class OncologyAnalyzer
{
    #region HLA nomenclature parsing and allele-specific HLA LOH (ONCO-HLA-001)

    /// <summary>
    /// Minimum number of colon-separated numeric fields in a valid HLA allele name. Source: WHO HLA
    /// Nomenclature (hla.alleles.org "Naming Alleles"); Marsh et al. (2010), Tissue Antigens 75(4):291–455 —
    /// "All alleles receive at least a four digit name, which corresponds to the first two sets of digits".
    /// </summary>
    private const int HlaMinFieldCount = 2;

    /// <summary>
    /// Maximum number of colon-separated numeric fields in a valid HLA allele name (type : protein :
    /// synonymous coding : non-coding). Source: WHO HLA Nomenclature — "up to four sets of digits separated
    /// by colons" (hla.alleles.org "Naming Alleles"; Marsh et al. 2010).
    /// </summary>
    private const int HlaMaxFieldCount = 4;

    /// <summary>
    /// Allele copy-number threshold below which an HLA allele is classified as lost. Source (verbatim):
    /// McGranahan et al. (2017), Cell 171(6):1259–1271 (PMC5720478) — "A copy number &lt; 0.5, is classified
    /// as subject to loss, and thereby indicative of LOH." The comparison is strict (CN must be &lt; 0.5).
    /// </summary>
    public const double HlaLohCopyNumberThreshold = 0.5;

    /// <summary>
    /// Allelic-imbalance significance threshold required before HLA LOH may be called. Source (verbatim):
    /// McGranahan et al. (2017), Cell 171(6):1259–1271 — "To avoid over-calling LOH, we also calculate a p
    /// value relating to allelic imbalance for each HLA gene. Allelic imbalance is determined if p &lt; 0.01
    /// using the paired Student's t-Test." The comparison is strict (p must be &lt; 0.01). Corroborated by the
    /// paired <c>t.test(..., paired=TRUE)</c> in mskcc/lohhla <c>LOHHLAscript.R</c>.
    /// </summary>
    public const double HlaLohAllelicImbalancePValueThreshold = 0.01;

    /// <summary>HLA expression-status suffix per WHO HLA Nomenclature ("Naming Alleles", hla.alleles.org).</summary>
    public enum HlaExpressionSuffix
    {
        /// <summary>No suffix — normally expressed allele.</summary>
        None,

        /// <summary><c>N</c> — Null allele (not expressed).</summary>
        Null,

        /// <summary><c>L</c> — Low cell-surface expression.</summary>
        Low,

        /// <summary><c>S</c> — Secreted molecule only (not on the cell surface).</summary>
        Secreted,

        /// <summary><c>C</c> — present in the Cytoplasm but not on the cell surface.</summary>
        Cytoplasm,

        /// <summary><c>A</c> — Aberrant expression (uncertain whether expressed).</summary>
        Aberrant,

        /// <summary><c>Q</c> — Questionable expression.</summary>
        Questionable,
    }

    /// <summary>
    /// Which of the two homologous HLA alleles at a locus was lost (had copy number below the loss threshold).
    /// </summary>
    public enum HlaLostAllele
    {
        /// <summary>No allele was lost (locus retained, no allele-specific LOH).</summary>
        None,

        /// <summary>The first allele (allele 1) is the lost allele.</summary>
        Allele1,

        /// <summary>The second allele (allele 2) is the lost allele.</summary>
        Allele2,

        /// <summary>Both alleles fell below the loss threshold (homozygous loss; not allele-specific LOH).</summary>
        Both,
    }

    /// <summary>
    /// A parsed HLA allele name per the WHO HLA Nomenclature. Format: <c>HLA-&lt;Gene&gt;*F1:F2[:F3[:F4]][suffix]</c>.
    /// Source: hla.alleles.org "Naming Alleles"; Marsh et al. (2010), Tissue Antigens 75(4):291–455.
    /// </summary>
    /// <param name="Gene">Gene name (e.g. "A", "B", "C", "DRB1"), upper-cased.</param>
    /// <param name="Fields">The 2–4 numeric fields, each as its original digit string (e.g. ["02","01"]).
    /// Field 1 = type/allele group; Field 2 = specific HLA protein; Field 3 = synonymous coding substitutions;
    /// Field 4 = non-coding differences.</param>
    /// <param name="Suffix">Optional expression-status suffix (None when absent).</param>
    public readonly record struct HlaAllele(string Gene, IReadOnlyList<string> Fields, HlaExpressionSuffix Suffix)
    {
        /// <summary>The allele-group (first) field — the serological "type" digits.</summary>
        public string AlleleGroup => Fields[0];

        /// <summary>The specific-HLA-protein (second) field.</summary>
        public string Protein => Fields[1];

        /// <summary>The canonical normalized allele name, e.g. <c>HLA-A*02:01:01:02L</c>.</summary>
        public string Name
        {
            get
            {
                string body = "HLA-" + Gene + "*" + string.Join(":", Fields);
                return Suffix == HlaExpressionSuffix.None ? body : body + SuffixLetter(Suffix);
            }
        }
    }

    /// <summary>
    /// Caller-supplied allele-specific copy-number evidence at one HLA gene, as produced by an HLA copy-number
    /// caller (e.g. LOHHLA): the estimated copy number of each of the two homologous alleles plus the
    /// allelic-imbalance p value. Source: McGranahan et al. (2017) — LOHHLA reports per-allele copy number
    /// (<c>HLA_type1copyNum</c>, <c>HLA_type2copyNum</c>) and a paired-t-test allelic-imbalance p value.
    /// </summary>
    /// <param name="Allele1">Allele 1 name (informational).</param>
    /// <param name="Allele1CopyNumber">Estimated copy number of allele 1 (≥ 0).</param>
    /// <param name="Allele2">Allele 2 name (informational).</param>
    /// <param name="Allele2CopyNumber">Estimated copy number of allele 2 (≥ 0).</param>
    /// <param name="AllelicImbalancePValue">Paired Student's t-test p value for allelic imbalance, in [0, 1].</param>
    public readonly record struct HlaAlleleCopyNumber(
        string Allele1,
        double Allele1CopyNumber,
        string Allele2,
        double Allele2CopyNumber,
        double AllelicImbalancePValue);

    /// <summary>Result of an allele-specific HLA LOH determination (LOHHLA classification).</summary>
    /// <param name="IsLoh">True iff allele-specific HLA LOH was called (one allele lost with significant imbalance).</param>
    /// <param name="LostAllele">Which allele was lost (<see cref="HlaLostAllele.None"/> when no LOH).</param>
    /// <param name="AllelicImbalanceSignificant">True iff the allelic-imbalance p value is below the threshold.</param>
    public readonly record struct HlaLohResult(bool IsLoh, HlaLostAllele LostAllele, bool AllelicImbalanceSignificant);

    private static char SuffixLetter(HlaExpressionSuffix suffix) => suffix switch
    {
        HlaExpressionSuffix.Null => 'N',
        HlaExpressionSuffix.Low => 'L',
        HlaExpressionSuffix.Secreted => 'S',
        HlaExpressionSuffix.Cytoplasm => 'C',
        HlaExpressionSuffix.Aberrant => 'A',
        HlaExpressionSuffix.Questionable => 'Q',
        _ => '\0',
    };

    private static bool TryMapSuffix(char letter, out HlaExpressionSuffix suffix)
    {
        // WHO HLA Nomenclature expression-status suffixes (hla.alleles.org "Naming Alleles").
        switch (letter)
        {
            case 'N': suffix = HlaExpressionSuffix.Null; return true;
            case 'L': suffix = HlaExpressionSuffix.Low; return true;
            case 'S': suffix = HlaExpressionSuffix.Secreted; return true;
            case 'C': suffix = HlaExpressionSuffix.Cytoplasm; return true;
            case 'A': suffix = HlaExpressionSuffix.Aberrant; return true;
            case 'Q': suffix = HlaExpressionSuffix.Questionable; return true;
            default: suffix = HlaExpressionSuffix.None; return false;
        }
    }

    /// <summary>
    /// Parses and validates an HLA allele name per the WHO HLA Nomenclature
    /// (<c>HLA-&lt;Gene&gt;*F1:F2[:F3[:F4]][suffix]</c>). The input is trimmed and the gene name upper-cased;
    /// the <c>HLA-</c> prefix is mandatory, the gene is separated from the fields by <c>*</c>, fields are
    /// colon-separated digit groups, and an optional trailing expression-status letter (N/L/S/C/A/Q) may
    /// follow the last field. Source: hla.alleles.org "Naming Alleles"; Marsh et al. (2010), Tissue Antigens
    /// 75(4):291–455.
    /// </summary>
    /// <param name="alleleName">The HLA allele name to parse.</param>
    /// <returns>The parsed <see cref="HlaAllele"/>.</returns>
    /// <exception cref="ArgumentNullException"><paramref name="alleleName"/> is null.</exception>
    /// <exception cref="ArgumentException"><paramref name="alleleName"/> is empty or whitespace.</exception>
    /// <exception cref="FormatException">The name does not conform to WHO HLA nomenclature (missing
    /// <c>HLA-</c> prefix or <c>*</c>; fewer than 2 or more than 4 fields; a non-numeric field; an invalid
    /// trailing suffix).</exception>
    public static HlaAllele ParseHlaAllele(string alleleName)
    {
        ArgumentNullException.ThrowIfNull(alleleName);
        if (string.IsNullOrWhiteSpace(alleleName))
        {
            throw new ArgumentException("HLA allele name must not be empty or whitespace.", nameof(alleleName));
        }

        string text = alleleName.Trim();

        // The "HLA-" prefix is mandatory (case-insensitive); the gene is separated from the fields by '*'.
        const string Prefix = "HLA-";
        if (!text.StartsWith(Prefix, StringComparison.OrdinalIgnoreCase))
        {
            throw new FormatException($"HLA allele name must start with '{Prefix}': '{alleleName}'.");
        }

        string afterPrefix = text.Substring(Prefix.Length);
        int star = afterPrefix.IndexOf('*');
        if (star <= 0 || star == afterPrefix.Length - 1)
        {
            throw new FormatException($"HLA allele name must have a gene and field block separated by '*': '{alleleName}'.");
        }

        string gene = afterPrefix.Substring(0, star).ToUpperInvariant();
        string fieldBlock = afterPrefix.Substring(star + 1);

        // Strip an optional single trailing expression-status suffix (only when the last char is a letter).
        var suffix = HlaExpressionSuffix.None;
        char last = fieldBlock[fieldBlock.Length - 1];
        if (char.IsLetter(last))
        {
            if (!TryMapSuffix(char.ToUpperInvariant(last), out suffix))
            {
                throw new FormatException($"Invalid HLA expression-status suffix '{last}' (allowed: N, L, S, C, A, Q): '{alleleName}'.");
            }

            fieldBlock = fieldBlock.Substring(0, fieldBlock.Length - 1);
        }

        string[] fields = fieldBlock.Split(':');
        if (fields.Length < HlaMinFieldCount || fields.Length > HlaMaxFieldCount)
        {
            throw new FormatException(
                $"HLA allele name must have between {HlaMinFieldCount} and {HlaMaxFieldCount} colon-separated fields, found {fields.Length}: '{alleleName}'.");
        }

        foreach (string field in fields)
        {
            // WHO HLA Nomenclature fields are ASCII decimal-digit groups ('0'–'9'); char.IsDigit also accepts
            // non-ASCII Unicode decimal digits (e.g. fullwidth '０', Arabic-Indic '٢'), which are not valid
            // nomenclature and must be rejected as malformed.
            if (field.Length == 0 || !field.All(static c => c is >= '0' and <= '9'))
            {
                throw new FormatException($"HLA allele field '{field}' must be a non-empty digit group: '{alleleName}'.");
            }
        }

        return new HlaAllele(gene, fields, suffix);
    }

    /// <summary>
    /// Non-throwing variant of <see cref="ParseHlaAllele(string)"/>. Returns true and the parsed allele on
    /// success; false and the default value when the name is null/empty or does not conform to WHO HLA
    /// nomenclature.
    /// </summary>
    /// <param name="alleleName">The HLA allele name to parse.</param>
    /// <param name="allele">The parsed allele on success; default otherwise.</param>
    /// <returns>True if parsing succeeded; otherwise false.</returns>
    public static bool TryParseHlaAllele(string? alleleName, out HlaAllele allele)
    {
        if (alleleName is null)
        {
            allele = default;
            return false;
        }

        try
        {
            allele = ParseHlaAllele(alleleName);
            return true;
        }
        catch (Exception ex) when (ex is FormatException or ArgumentException)
        {
            allele = default;
            return false;
        }
    }

    /// <summary>
    /// Classifies allele-specific HLA loss of heterozygosity from caller-supplied per-allele copy number and
    /// the allelic-imbalance p value, per the LOHHLA decision rule. HLA LOH is called iff <b>exactly one</b>
    /// of the two homologous alleles has copy number strictly below
    /// <see cref="HlaLohCopyNumberThreshold"/> (0.5) <b>and</b> the allelic-imbalance p value is strictly
    /// below <see cref="HlaLohAllelicImbalancePValueThreshold"/> (0.01). If both alleles fall below the loss
    /// threshold the locus is reported as homozygous loss (<see cref="HlaLostAllele.Both"/>), which is not
    /// allele-specific LOH. Source: McGranahan et al. (2017), Cell 171(6):1259–1271 (PMC5720478) — "A copy
    /// number &lt; 0.5, is classified as subject to loss, and thereby indicative of LOH" and "Allelic
    /// imbalance is determined if p &lt; 0.01 using the paired Student's t-Test".
    /// </summary>
    /// <param name="alleleCopyNumber">Caller-supplied per-allele copy number and allelic-imbalance p value.</param>
    /// <returns>The HLA LOH classification.</returns>
    /// <exception cref="ArgumentException">A copy number is negative, or the p value is outside [0, 1].</exception>
    public static HlaLohResult DetectHlaLoh(HlaAlleleCopyNumber alleleCopyNumber)
    {
        // Copy numbers must be finite and non-negative. NaN must be rejected explicitly: `NaN < 0` is false, so
        // a NaN copy number would otherwise slip past the bound and leak into the loss-threshold comparison.
        if (!double.IsFinite(alleleCopyNumber.Allele1CopyNumber) || !double.IsFinite(alleleCopyNumber.Allele2CopyNumber)
            || alleleCopyNumber.Allele1CopyNumber < 0 || alleleCopyNumber.Allele2CopyNumber < 0)
        {
            throw new ArgumentException("HLA allele copy numbers must be finite and non-negative.", nameof(alleleCopyNumber));
        }

        // The p value must be finite and in [0, 1]. NaN must be rejected explicitly: both `NaN < 0.0` and
        // `NaN > 1.0` are false, so a NaN p value would otherwise pass this check and silently be treated as
        // not-significant in the `< 0.01` test rather than reported as malformed input.
        if (!double.IsFinite(alleleCopyNumber.AllelicImbalancePValue)
            || alleleCopyNumber.AllelicImbalancePValue is < 0.0 or > 1.0)
        {
            throw new ArgumentException("Allelic-imbalance p value must be in [0, 1].", nameof(alleleCopyNumber));
        }

        // LOHHLA over-calling guard: significant allelic imbalance (paired t-test p < 0.01) is required.
        bool imbalanceSignificant = alleleCopyNumber.AllelicImbalancePValue < HlaLohAllelicImbalancePValueThreshold;

        bool allele1Lost = alleleCopyNumber.Allele1CopyNumber < HlaLohCopyNumberThreshold;
        bool allele2Lost = alleleCopyNumber.Allele2CopyNumber < HlaLohCopyNumberThreshold;

        if (allele1Lost && allele2Lost)
        {
            // Both homologs below 0.5 → homozygous loss, not allele-specific LOH (Evidence assumption).
            return new HlaLohResult(false, HlaLostAllele.Both, imbalanceSignificant);
        }

        if (imbalanceSignificant && allele1Lost)
        {
            return new HlaLohResult(true, HlaLostAllele.Allele1, true);
        }

        if (imbalanceSignificant && allele2Lost)
        {
            return new HlaLohResult(true, HlaLostAllele.Allele2, true);
        }

        return new HlaLohResult(false, HlaLostAllele.None, imbalanceSignificant);
    }

    #endregion


    #region Clinical Actionability (OncoKB Therapeutic Levels of Evidence)

    /// <summary>
    /// OncoKB therapeutic level of evidence assigned to a biomarker–drug association, indicating how
    /// strongly an alteration is predictive of sensitivity (Levels 1, 2, 3A, 3B, 4) or resistance
    /// (Levels R1, R2) to a therapy. Source: Chakravarty D et al. (2017), JCO Precis Oncol 2017:1–16;
    /// definitions verbatim from the OncoKB Therapeutic Levels of Evidence (V2) document.
    /// </summary>
    public enum OncoKbLevel
    {
        /// <summary>No leveled therapeutic association (variant is not actionable on this axis).</summary>
        None,

        /// <summary>
        /// Level R2 — Investigational Resistance: "Compelling clinical evidence supports the biomarker as
        /// being predictive of resistance to a drug." Lowest in the combined order. Source: OncoKB Levels V2.
        /// </summary>
        R2,

        /// <summary>
        /// Level 4 — Hypothetical: "Compelling biological evidence supports the biomarker as being
        /// predictive of response to a drug." Source: OncoKB Levels V2.
        /// </summary>
        Level4,

        /// <summary>
        /// Level 3B — Investigational: "Standard care or investigational biomarker predictive of response
        /// to an FDA-approved or investigational drug in another indication." Source: OncoKB Levels V2.
        /// </summary>
        Level3B,

        /// <summary>
        /// Level 3A — Investigational: "Compelling clinical evidence supports the biomarker as being
        /// predictive of response to a drug in this indication." Source: OncoKB Levels V2.
        /// </summary>
        Level3A,

        /// <summary>
        /// Level 2 — Standard Care: "Standard care biomarker recommended by the NCCN or other professional
        /// guidelines predictive of response to an FDA-approved drug in this indication." Source: OncoKB Levels V2.
        /// </summary>
        Level2,

        /// <summary>
        /// Level 1 — Standard Care: "FDA-recognized biomarker predictive of response to an FDA-approved drug
        /// in this indication." Source: OncoKB Levels V2.
        /// </summary>
        Level1,

        /// <summary>
        /// Level R1 — Standard Care Resistance: "Standard care biomarker predictive of resistance to an
        /// FDA-approved drug in this indication." Highest in the combined order. Source: OncoKB Levels V2.
        /// </summary>
        R1
    }

    /// <summary>Levels that denote sensitivity (response) to a therapy. Source: OncoKB Levels V2 (1/2/3A/3B/4).</summary>
    private static readonly HashSet<OncoKbLevel> SensitivityLevels = new()
    {
        OncoKbLevel.Level1, OncoKbLevel.Level2, OncoKbLevel.Level3A, OncoKbLevel.Level3B, OncoKbLevel.Level4
    };

    /// <summary>Levels that denote resistance to a therapy. Source: OncoKB Levels V2 (R1/R2).</summary>
    private static readonly HashSet<OncoKbLevel> ResistanceLevels = new()
    {
        OncoKbLevel.R1, OncoKbLevel.R2
    };

    /// <summary>
    /// Levels categorized as "standard care" by OncoKB (as opposed to investigational/hypothetical):
    /// Levels 1, 2 (sensitivity) and R1 (resistance). Source: OncoKB Curation SOP v3 — "The highest levels
    /// of evidence, Levels 1 and 2, refer to the standard implications... Level R1 refers to the standard
    /// implications for resistance"; Levels 3A/3B/4/R2 are investigational/hypothetical.
    /// </summary>
    private static readonly HashSet<OncoKbLevel> StandardCareLevels = new()
    {
        OncoKbLevel.Level1, OncoKbLevel.Level2, OncoKbLevel.R1
    };

    /// <summary>
    /// One caller-supplied biomarker–drug therapeutic association from a precision-oncology knowledgebase
    /// (e.g. an OncoKB export). The library does not embed the OncoKB curated content (3,000+ alterations
    /// across 418 genes, Chakravarty 2017); the caller performs the lookup and supplies the relevant rows.
    /// </summary>
    /// <param name="Drug">Therapy name the association refers to.</param>
    /// <param name="Level">OncoKB therapeutic level of evidence for this drug–variant association.</param>
    public readonly record struct TherapyAssociation(string Drug, OncoKbLevel Level);

    /// <summary>
    /// Caller-supplied evidence for one variant's clinical actionability: the gene/protein change (for
    /// reporting) and the set of leveled drug associations curated for it. Mirrors the framework boundary of
    /// <see cref="CancerVariantAnnotationInput"/> — actionability comes from a caller-supplied knowledgebase.
    /// </summary>
    /// <param name="Gene">Gene symbol the variant falls in.</param>
    /// <param name="ProteinChange">HGVS protein change (e.g. p.V600E); informational.</param>
    /// <param name="Associations">Leveled drug associations from the knowledgebase (may be empty, never null).</param>
    public readonly record struct VariantActionabilityInput
    {
        /// <summary>Gene symbol the variant falls in.</summary>
        public string Gene { get; }

        /// <summary>HGVS protein change (e.g. p.V600E); informational.</summary>
        public string ProteinChange { get; }

        /// <summary>Leveled drug associations from the caller-supplied knowledgebase.</summary>
        public IReadOnlyList<TherapyAssociation> Associations { get; }

        /// <summary>Creates an actionability input. <paramref name="associations"/> must not be null.</summary>
        /// <exception cref="ArgumentNullException"><paramref name="associations"/> is null.</exception>
        public VariantActionabilityInput(
            string gene, string proteinChange, IReadOnlyList<TherapyAssociation> associations)
        {
            ArgumentNullException.ThrowIfNull(associations);
            Gene = gene;
            ProteinChange = proteinChange;
            Associations = associations;
        }
    }

    /// <summary>
    /// Per-variant clinical actionability assessment under the OncoKB therapeutic levels system.
    /// </summary>
    /// <param name="Variant">The variant evidence that was assessed.</param>
    /// <param name="HighestSensitiveLevel">Highest sensitivity level (1 &gt; 2 &gt; 3A &gt; 3B &gt; 4), or None.</param>
    /// <param name="HighestResistanceLevel">Highest resistance level (R1 &gt; R2), or None.</param>
    /// <param name="HighestCombinedLevel">Highest level over both axes (R1 &gt; 1 &gt; 2 &gt; 3A &gt; 3B &gt; 4 &gt; R2), or None.</param>
    public readonly record struct ActionabilityAssessment(
        VariantActionabilityInput Variant,
        OncoKbLevel HighestSensitiveLevel,
        OncoKbLevel HighestResistanceLevel,
        OncoKbLevel HighestCombinedLevel)
    {
        /// <summary>True when the variant has at least one leveled therapeutic association.</summary>
        public bool IsActionable => HighestCombinedLevel != OncoKbLevel.None;
    }

    /// <summary>
    /// Compares two OncoKB levels by the combined actionability order R1 &gt; 1 &gt; 2 &gt; 3A &gt; 3B &gt; 4 &gt; R2.
    /// Returns a positive number when <paramref name="a"/> is more actionable than <paramref name="b"/>,
    /// negative when less, zero when equal. <see cref="OncoKbLevel.None"/> ranks below every leveled value.
    /// Source: oncokb-annotator README HIGHEST_LEVEL order.
    /// </summary>
    /// <param name="a">First level.</param>
    /// <param name="b">Second level.</param>
    /// <returns>Sign indicates which level is higher in the combined order.</returns>
    public static int CompareLevels(OncoKbLevel a, OncoKbLevel b)
        // The enum is declared in ascending actionability (None lowest, R1 highest), so the underlying
        // integer value already encodes the combined order; comparing values is the documented precedence.
        => ((int)a).CompareTo((int)b);

    /// <summary>
    /// Returns the highest OncoKB level over a set of levels using the combined order R1 &gt; 1 &gt; 2 &gt;
    /// 3A &gt; 3B &gt; 4 &gt; R2, restricted to the levels in <paramref name="allowed"/> (used to compute the
    /// sensitivity-only and resistance-only maxima). Returns <see cref="OncoKbLevel.None"/> when no allowed
    /// level is present. Source: oncokb-annotator README HIGHEST_*_LEVEL orders.
    /// </summary>
    private static OncoKbLevel HighestLevel(
        IReadOnlyList<TherapyAssociation> associations, HashSet<OncoKbLevel>? allowed)
    {
        OncoKbLevel best = OncoKbLevel.None;
        foreach (var association in associations)
        {
            if (allowed is not null && !allowed.Contains(association.Level))
            {
                continue;
            }

            if (CompareLevels(association.Level, best) > 0)
            {
                best = association.Level;
            }
        }

        return best;
    }

    /// <summary>
    /// Classifies a single variant's clinical actionability to the highest OncoKB therapeutic level over all
    /// its caller-supplied drug associations, under the combined order R1 &gt; 1 &gt; 2 &gt; 3A &gt; 3B &gt; 4
    /// &gt; R2. Returns <see cref="OncoKbLevel.None"/> when the variant has no leveled association (not
    /// actionable). Source: Chakravarty D et al. (2017); oncokb-annotator README HIGHEST_LEVEL.
    /// </summary>
    /// <param name="variant">Caller-supplied variant actionability evidence.</param>
    /// <returns>The highest combined OncoKB level, or <see cref="OncoKbLevel.None"/>.</returns>
    public static OncoKbLevel ClassifyActionabilityLevel(VariantActionabilityInput variant)
    {
        ArgumentNullException.ThrowIfNull(variant.Associations);
        return HighestLevel(variant.Associations, allowed: null);
    }

    /// <summary>
    /// Assesses the clinical actionability of each variant under the OncoKB therapeutic levels system,
    /// computing the highest sensitivity level (1 &gt; 2 &gt; 3A &gt; 3B &gt; 4), highest resistance level
    /// (R1 &gt; R2), and highest combined level (R1 &gt; 1 &gt; 2 &gt; 3A &gt; 3B &gt; 4 &gt; R2) from the
    /// caller-supplied drug associations. Output preserves input order, one entry per variant. Source:
    /// Chakravarty D et al. (2017); oncokb-annotator README HIGHEST_LEVEL / HIGHEST_SENSITIVE_LEVEL /
    /// HIGHEST_RESISTANCE_LEVEL.
    /// </summary>
    /// <param name="variants">Caller-supplied variant actionability evidence records.</param>
    /// <returns>One <see cref="ActionabilityAssessment"/> per input variant, in input order.</returns>
    /// <exception cref="ArgumentNullException"><paramref name="variants"/> is null.</exception>
    public static IReadOnlyList<ActionabilityAssessment> AssessActionability(
        IEnumerable<VariantActionabilityInput> variants)
    {
        ArgumentNullException.ThrowIfNull(variants);

        var assessments = new List<ActionabilityAssessment>();
        foreach (var variant in variants)
        {
            ArgumentNullException.ThrowIfNull(variant.Associations);

            OncoKbLevel sensitive = HighestLevel(variant.Associations, SensitivityLevels);
            OncoKbLevel resistance = HighestLevel(variant.Associations, ResistanceLevels);
            OncoKbLevel combined = HighestLevel(variant.Associations, allowed: null);

            assessments.Add(new ActionabilityAssessment(variant, sensitive, resistance, combined));
        }

        return assessments;
    }

    /// <summary>
    /// Returns the caller-supplied therapy associations for a variant ordered by descending OncoKB level
    /// (most actionable first) under the combined order R1 &gt; 1 &gt; 2 &gt; 3A &gt; 3B &gt; 4 &gt; R2. Thin
    /// presentation wrapper over the knowledgebase rows; returns an empty list (never null) when there are no
    /// associations. Source: oncokb-annotator README HIGHEST_LEVEL order.
    /// </summary>
    /// <param name="variant">Caller-supplied variant actionability evidence.</param>
    /// <returns>The associations ordered most-actionable first.</returns>
    public static IReadOnlyList<TherapyAssociation> GetTherapyRecommendations(VariantActionabilityInput variant)
    {
        ArgumentNullException.ThrowIfNull(variant.Associations);

        var ordered = new List<TherapyAssociation>(variant.Associations);
        // Descending by combined order: higher enum value = more actionable, so negate the comparison.
        ordered.Sort((x, y) => CompareLevels(y.Level, x.Level));
        return ordered;
    }

    /// <summary>
    /// True when the level is one of OncoKB's "standard care" levels (1, 2, R1) as opposed to the
    /// investigational/hypothetical levels (3A, 3B, 4, R2). Source: OncoKB Curation SOP v3.
    /// </summary>
    /// <param name="level">An OncoKB therapeutic level.</param>
    /// <returns>True for standard-care levels (1, 2, R1).</returns>
    public static bool IsStandardCare(OncoKbLevel level) => StandardCareLevels.Contains(level);

    #endregion


    #region Complex Somatic Rearrangement (Chromothripsis) Classification

    /// <summary>
    /// Minimum number of oscillating copy-number changes (per-segment CN state transitions) required by
    /// the first-pass chromothripsis screen. Source: Magrangeas et al. (2011), Blood 118(3):675–678 — the
    /// lowest first-pass operational cutoff of "10, 20, or 50 oscillating copy number changes" cited by the
    /// Korbel &amp; Campbell (2013) framework (Cell 152:1226–1236; review PMC3861665). Default = 10.
    /// </summary>
    public const int MinOscillatingCopyNumberChanges = 10;

    /// <summary>
    /// Maximum number of distinct copy-number states permitted for a chromothripsis call. Korbel &amp;
    /// Campbell (2013) define the hallmark profile as oscillation between (canonically) two — and at most
    /// two-or-three — copy-number states, in contrast with progressive amplification (many ascending
    /// states). Source: Korbel &amp; Campbell (2013), Cell 152:1226–1236 (criterion B). Default = 3.
    /// </summary>
    public const int MaxChromothripsisCopyNumberStates = 3;

    /// <summary>
    /// Minimum number of clustered intrachromosomal structural variants for an event to be eligible for a
    /// chromothripsis call. Source: Cortés-Ciriano et al. (2020), Nat. Genet. 52:331–341 — focal events
    /// "comprising fewer than six SVs" are excluded. Default = 6.
    /// </summary>
    public const int MinChromothripsisSvBurden = 6;

    /// <summary>
    /// Number of adjacent oscillating copy-number segments at or above which a chromothripsis call is
    /// high-confidence. Source: Cortés-Ciriano et al. (2020) — high-confidence calls "display oscillations
    /// between two states in at least seven adjacent segments". Default = 7.
    /// </summary>
    public const int HighConfidenceOscillatingSegments = 7;

    /// <summary>
    /// Minimum number of adjacent oscillating copy-number segments for a low-confidence chromothripsis
    /// signal. Source: Cortés-Ciriano et al. (2020) — low-confidence calls "involve between four and six
    /// segments". Default = 4 (the band [4, 6] is low-confidence; ≥ 7 is high-confidence).
    /// </summary>
    public const int LowConfidenceOscillatingSegments = 4;

    /// <summary>
    /// Coefficient of variation of inter-breakpoint distances expected under the random-breakpoint null.
    /// Source: Korbel &amp; Campbell (2013) — "the null hypothesis of random breakpoints predicts that the
    /// distance between breakpoints should be distributed exponentially"; the exponential distribution has
    /// CV = 1, so over-dispersion toward many short gaps with few long gaps (clustering) gives CV &gt; 1.
    /// </summary>
    private const double ExponentialNullCoefficientOfVariation = 1.0;

    /// <summary>
    /// Classification of a chromosome's somatic structural-rearrangement profile.
    /// </summary>
    public enum ComplexRearrangementType
    {
        /// <summary>Does not meet the chromothripsis hallmark criteria.</summary>
        NotComplex,

        /// <summary>
        /// Chromothripsis: clustered breakpoints with oscillation between ≤ 3 (canonically 2) copy-number
        /// states and sufficient oscillation/SV burden (Korbel &amp; Campbell 2013; Cortés-Ciriano 2020).
        /// </summary>
        Chromothripsis
    }

    /// <summary>
    /// Confidence tier for a chromothripsis signal, from the number of adjacent oscillating segments.
    /// </summary>
    public enum ChromothripsisConfidence
    {
        /// <summary>Fewer than four adjacent oscillating segments — no chromothripsis signal.</summary>
        None,

        /// <summary>Four to six adjacent oscillating segments — low-confidence (Cortés-Ciriano 2020).</summary>
        Low,

        /// <summary>Seven or more adjacent oscillating segments — high-confidence (Cortés-Ciriano 2020).</summary>
        High
    }

    /// <summary>
    /// Input for <see cref="ClassifyComplexRearrangement"/>: the per-segment copy-number states along one
    /// chromosomal region and the number of clustered intrachromosomal structural variants supporting it.
    /// </summary>
    /// <param name="SegmentCopyNumbers">Per-segment integer copy numbers in genomic order along the region.</param>
    /// <param name="StructuralVariantCount">Number of clustered intrachromosomal SVs in the region.</param>
    public readonly record struct ComplexRearrangementInput(
        IReadOnlyList<int> SegmentCopyNumbers,
        int StructuralVariantCount);

    /// <summary>
    /// Result of complex-rearrangement classification.
    /// </summary>
    /// <param name="Type">The classification (chromothripsis or not).</param>
    /// <param name="Confidence">The confidence tier derived from adjacent oscillating-segment count.</param>
    /// <param name="OscillationCount">Number of per-segment copy-number state transitions.</param>
    /// <param name="OscillatingSegmentCount">Number of segments participating in the oscillation.</param>
    /// <param name="DistinctStateCount">Number of distinct copy-number states in the profile.</param>
    /// <param name="StructuralVariantCount">Clustered intrachromosomal SV burden of the region.</param>
    public readonly record struct ComplexRearrangementResult(
        ComplexRearrangementType Type,
        ChromothripsisConfidence Confidence,
        int OscillationCount,
        int OscillatingSegmentCount,
        int DistinctStateCount,
        int StructuralVariantCount);

    /// <summary>
    /// Summary of a breakpoint-clustering test against the random-breakpoint exponential null.
    /// </summary>
    /// <param name="BreakpointCount">Number of breakpoints provided.</param>
    /// <param name="MeanGap">Mean inter-breakpoint distance.</param>
    /// <param name="CoefficientOfVariation">Standard deviation / mean of inter-breakpoint distances.</param>
    /// <param name="IsClustered">True when CV &gt; 1 (over-dispersed relative to the exponential null).</param>
    public readonly record struct BreakpointClusteringResult(
        int BreakpointCount,
        double MeanGap,
        double CoefficientOfVariation,
        bool IsClustered);

    /// <summary>
    /// Counts the number of oscillating copy-number changes along a region: the number of adjacent segments
    /// whose copy-number state differs from the immediately preceding segment. This is the "oscillating
    /// copy number changes" quantity used by the first-pass chromothripsis screen.
    /// Source: Magrangeas et al. (2011), Blood 118(3):675–678; Korbel &amp; Campbell (2013), Cell 152:1226–1236.
    /// For <c>n</c> segments the count is in [0, n−1]; fewer than two segments yields 0.
    /// </summary>
    /// <param name="segmentCopyNumbers">Per-segment integer copy numbers in genomic order.</param>
    /// <returns>The number of adjacent copy-number state transitions.</returns>
    /// <exception cref="ArgumentNullException"><paramref name="segmentCopyNumbers"/> is null.</exception>
    public static int CountCopyNumberStateOscillations(IReadOnlyList<int> segmentCopyNumbers)
    {
        ArgumentNullException.ThrowIfNull(segmentCopyNumbers);

        int transitions = 0;
        for (int i = 1; i < segmentCopyNumbers.Count; i++)
        {
            if (segmentCopyNumbers[i] != segmentCopyNumbers[i - 1])
            {
                transitions++;
            }
        }

        return transitions;
    }

    /// <summary>
    /// Tests a set of genomic breakpoint positions for clustering against the random-breakpoint null.
    /// Under the null of uniformly random breakpoints the inter-breakpoint distances are exponentially
    /// distributed, which has a coefficient of variation (CV = sd/mean) of 1; over-dispersion toward many
    /// short gaps with a few long gaps (a tight cluster plus outliers) gives CV &gt; 1, which flags
    /// clustering. Source: Korbel &amp; Campbell (2013), Cell 152:1226–1236 (criterion A, exponential null).
    /// </summary>
    /// <param name="breakpointPositions">Genomic breakpoint coordinates (any order); sorted internally.</param>
    /// <returns>The clustering summary; with fewer than three breakpoints clustering cannot be assessed
    /// (fewer than two gaps), so <see cref="BreakpointClusteringResult.IsClustered"/> is false.</returns>
    /// <exception cref="ArgumentNullException"><paramref name="breakpointPositions"/> is null.</exception>
    public static BreakpointClusteringResult TestBreakpointClustering(IReadOnlyList<long> breakpointPositions)
    {
        ArgumentNullException.ThrowIfNull(breakpointPositions);

        // Need at least two gaps (three breakpoints) to define a CV; otherwise clustering is undefined.
        if (breakpointPositions.Count < 3)
        {
            return new BreakpointClusteringResult(breakpointPositions.Count, 0.0, 0.0, false);
        }

        var sorted = breakpointPositions.OrderBy(p => p).ToArray();
        int gapCount = sorted.Length - 1;
        var gaps = new double[gapCount];
        for (int i = 0; i < gapCount; i++)
        {
            gaps[i] = sorted[i + 1] - sorted[i];
        }

        double mean = gaps.Average();
        if (mean <= 0.0)
        {
            // All breakpoints coincide: degenerate, treat as not assessable.
            return new BreakpointClusteringResult(breakpointPositions.Count, 0.0, 0.0, false);
        }

        double variance = gaps.Sum(g => (g - mean) * (g - mean)) / gapCount;
        double cv = Math.Sqrt(variance) / mean;
        bool clustered = cv > ExponentialNullCoefficientOfVariation;

        return new BreakpointClusteringResult(breakpointPositions.Count, mean, cv, clustered);
    }

    /// <summary>
    /// Classifies a chromosomal region's somatic structural-rearrangement profile as chromothripsis or not,
    /// applying the Korbel &amp; Campbell (2013) hallmark criteria together with the Cortés-Ciriano (2020)
    /// operational thresholds. A region is called <see cref="ComplexRearrangementType.Chromothripsis"/> when
    /// ALL of the following hold: (i) the copy-number profile oscillates between at most
    /// <see cref="MaxChromothripsisCopyNumberStates"/> (canonically 2) distinct states — the two-state
    /// hallmark, excluding progressive amplification; (ii) it has at least
    /// <see cref="MinOscillatingCopyNumberChanges"/> oscillating copy-number changes (first-pass screen);
    /// and (iii) the clustered intrachromosomal SV burden is at least <see cref="MinChromothripsisSvBurden"/>.
    /// The confidence tier is derived independently from the number of adjacent oscillating segments
    /// (≥ <see cref="HighConfidenceOscillatingSegments"/> → High; [<see cref="LowConfidenceOscillatingSegments"/>, 6] → Low; else None).
    /// Source: Korbel &amp; Campbell (2013), Cell 152:1226–1236; Cortés-Ciriano et al. (2020), Nat. Genet. 52:331–341;
    /// Magrangeas et al. (2011), Blood 118(3):675–678.
    /// </summary>
    /// <param name="input">Per-segment copy numbers and the clustered SV burden of the region.</param>
    /// <returns>The classification result.</returns>
    /// <exception cref="ArgumentNullException"><paramref name="input"/>.SegmentCopyNumbers is null.</exception>
    public static ComplexRearrangementResult ClassifyComplexRearrangement(ComplexRearrangementInput input)
    {
        ArgumentNullException.ThrowIfNull(input.SegmentCopyNumbers);

        var states = input.SegmentCopyNumbers;
        int oscillations = CountCopyNumberStateOscillations(states);

        // Segments participating in an oscillation: a run of k transitions spans k+1 segments.
        int oscillatingSegments = oscillations > 0 ? oscillations + 1 : 0;

        int distinctStates = states.Count == 0 ? 0 : states.Distinct().Count();

        // Confidence tier from adjacent oscillating-segment count (Cortés-Ciriano 2020).
        ChromothripsisConfidence confidence;
        if (oscillatingSegments >= HighConfidenceOscillatingSegments)
        {
            confidence = ChromothripsisConfidence.High;
        }
        else if (oscillatingSegments >= LowConfidenceOscillatingSegments)
        {
            confidence = ChromothripsisConfidence.Low;
        }
        else
        {
            confidence = ChromothripsisConfidence.None;
        }

        // Chromothripsis hallmark gate: two-state oscillation + ≥10 oscillations + ≥6 clustered SVs.
        bool isChromothripsis =
            distinctStates >= 2 &&
            distinctStates <= MaxChromothripsisCopyNumberStates &&
            oscillations >= MinOscillatingCopyNumberChanges &&
            input.StructuralVariantCount >= MinChromothripsisSvBurden;

        var type = isChromothripsis
            ? ComplexRearrangementType.Chromothripsis
            : ComplexRearrangementType.NotComplex;

        return new ComplexRearrangementResult(
            type,
            confidence,
            oscillations,
            oscillatingSegments,
            distinctStates,
            input.StructuralVariantCount);
    }

    #endregion

}
