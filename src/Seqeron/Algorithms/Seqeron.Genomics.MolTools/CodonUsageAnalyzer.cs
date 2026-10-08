namespace Seqeron.Genomics.MolTools;

/// <summary>
/// Analyzes codon usage patterns in coding sequences.
/// Useful for studying codon bias, gene expression optimization, and evolutionary analysis.
/// </summary>
public static class CodonUsageAnalyzer
{
    #region Codon Usage Tables

    /// <summary>
    /// Counts codon occurrences in a coding sequence.
    /// </summary>
    /// <param name="sequence">Coding DNA sequence (must be multiple of 3).</param>
    /// <returns>Dictionary of codon counts.</returns>
    public static Dictionary<string, int> CountCodons(DnaSequence sequence)
    {
        ArgumentNullException.ThrowIfNull(sequence);
        return CountCodonsCore(sequence.Sequence);
    }

    /// <summary>
    /// Counts codon occurrences in a raw sequence string.
    /// </summary>
    /// <remarks>
    /// Input is case-insensitive and may be DNA (T) or RNA (U): U is read as T, as in
    /// CodonW <c>ident_codon</c> (Peden 1999, codon_us.c: 'T','t','U','u' → the same base).
    /// Codons are reported in DNA spelling. Triplets containing any other symbol are
    /// skipped without shifting the frame; a trailing partial triplet is ignored.
    /// Delegates to the canonical Core counter <see cref="Translator.CountCodons(string?, int)"/>.
    /// </remarks>
    public static Dictionary<string, int> CountCodons(string sequence) => CountCodonsCore(sequence);

    // Frame-0 codon counts (canonical Core counter: upper-case, U read as T, ambiguous triplets
    // skipped frame-preservingly, trailing partial triplet dropped).
    private static Dictionary<string, int> CountCodonsCore(string? sequence) =>
        Translator.CountCodons(sequence);

    /// <summary>
    /// Frame-0 complete triplets of a coding sequence (DNA spelling, upper case; RNA U read as T);
    /// an ambiguous triplet is <c>null</c> (skipped by callers without shifting the frame).
    /// Delegates to the canonical Core splitter <see cref="Translator.SplitInFrameCodons(string?, int)"/>;
    /// shared by <see cref="CountCodons(string)"/>, CAI and the per-position codon screens of
    /// <see cref="CodonOptimizer"/>.
    /// </summary>
    internal static string?[] SplitInFrameCodons(string? sequence) =>
        Translator.SplitInFrameCodons(sequence);

    #endregion

    #region RSCU (Relative Synonymous Codon Usage)

    /// <summary>
    /// Calculates Relative Synonymous Codon Usage (RSCU) under the Standard genetic code
    /// (NCBI table 1). See <see cref="CalculateRscu(IReadOnlyDictionary{string, int}, GeneticCode)"/>.
    /// </summary>
    public static Dictionary<string, double> CalculateRscu(DnaSequence sequence) =>
        CalculateRscu(sequence, GeneticCode.Standard);

    /// <summary>
    /// Calculates RSCU of a coding sequence under the given genetic code.
    /// See <see cref="CalculateRscu(IReadOnlyDictionary{string, int}, GeneticCode)"/>.
    /// </summary>
    /// <exception cref="ArgumentNullException"><paramref name="sequence"/> or <paramref name="code"/> is null.</exception>
    public static Dictionary<string, double> CalculateRscu(DnaSequence sequence, GeneticCode code)
    {
        ArgumentNullException.ThrowIfNull(sequence);
        ArgumentNullException.ThrowIfNull(code);
        return CalculateRscu(CountCodonsCore(sequence.Sequence), code);
    }

    /// <summary>
    /// Calculates RSCU from a raw sequence string (DNA or RNA, case-insensitive) under the
    /// Standard genetic code (NCBI table 1). Null/empty returns an empty dictionary.
    /// </summary>
    public static Dictionary<string, double> CalculateRscu(string sequence) =>
        CalculateRscu(sequence, GeneticCode.Standard);

    /// <summary>
    /// Calculates RSCU from a raw sequence string (DNA or RNA, case-insensitive; counted as
    /// by <see cref="CountCodons(string)"/>) under the given genetic code.
    /// Null/empty returns an empty dictionary.
    /// </summary>
    /// <exception cref="ArgumentNullException"><paramref name="code"/> is null.</exception>
    public static Dictionary<string, double> CalculateRscu(string sequence, GeneticCode code)
    {
        ArgumentNullException.ThrowIfNull(code);
        if (string.IsNullOrEmpty(sequence))
            return new Dictionary<string, double>();

        return CalculateRscu(CountCodonsCore(sequence), code);
    }

    /// <summary>
    /// Calculates Relative Synonymous Codon Usage (RSCU) from codon counts.
    /// For codon j of an amino acid with n synonymous codons and observed counts x,
    /// RSCU = x_j / ((1/n) * sum_k x_k) = (n * x_j) / sum_k x_k
    /// (Sharp, Tuohy &amp; Mosurski 1986, Nucleic Acids Res. 14(13):5125-5143;
    /// Sharp &amp; Li 1987, Nucleic Acids Res. 15(3):1281-1295).
    /// RSCU = 1 means no bias, &gt; 1 over-represented, &lt; 1 under-represented.
    /// </summary>
    /// <remarks>
    /// Follows the CodonW reference implementation (Peden 1999, codon_us.c
    /// <c>rscu_usage_out</c>): synonymous families are taken from <paramref name="code"/>
    /// ("RSCU values are genetic code dependent"); all 64 codons are reported, including the
    /// termination codons, which form one synonymous family ('*'); single-codon families
    /// (e.g. Met, Trp in table 1) are 1 when present; every codon of a family that does not
    /// occur (0/0) is reported as 0. Codons that are context-dependent stops in NCBI tables
    /// 27/28/31 belong to the family of the amino acid they encode
    /// (<see cref="GeneticCode.CodonTable"/>; Biopython <c>forward_table</c>).
    /// </remarks>
    /// <param name="codonCounts">Codon counts keyed by uppercase DNA codon (as returned by
    /// <see cref="CountCodons(string)"/>); other keys are ignored.</param>
    /// <param name="code">Genetic code defining the synonymous families.</param>
    /// <returns>RSCU for each of the 64 codons (uppercase DNA spelling).</returns>
    /// <exception cref="ArgumentNullException">An argument is null.</exception>
    public static Dictionary<string, double> CalculateRscu(
        IReadOnlyDictionary<string, int> codonCounts, GeneticCode code)
    {
        ArgumentNullException.ThrowIfNull(codonCounts);
        ArgumentNullException.ThrowIfNull(code);

        var rscu = new Dictionary<string, double>(64);

        foreach (var family in SynonymousFamilies(code))
        {
            int familySize = family.Count();
            long familyTotal = 0;
            foreach (var codon in family)
                familyTotal += codonCounts.GetValueOrDefault(codon, 0);

            foreach (var codon in family)
            {
                // RSCU = n·x / Σx; an absent family (Σx = 0) is 0 for every member (CodonW).
                rscu[codon] = familyTotal > 0
                    ? (double)familySize * codonCounts.GetValueOrDefault(codon, 0) / familyTotal
                    : 0.0;
            }
        }

        return rscu;
    }

    // Synonymous codon families (DNA spelling) of a genetic code, keyed by the encoded
    // amino acid ('*' = termination), in NCBI codon order. Single grouping used by RSCU, CAI,
    // ENC and GC3s. Context-dependent stops of tables 27/28/31 fall in their amino-acid family.
    private static IEnumerable<IGrouping<char, string>> SynonymousFamilies(GeneticCode code) =>
        code.CodonTable.GroupBy(kv => kv.Value, kv => kv.Key.Replace('U', 'T'));

    #endregion

    #region CAI (Codon Adaptation Index)

    // CodonW 1.4.4 cai_out (Peden 1999, codon_us.c): "these codons have fitness of zero
    // (<.0001) are adjusted to 0.01" — the Bulmer (1988) substitution also used by seqinr
    // cai(zero.threshold = 0.0001, zero.to = 0.01).
    private const double EffectivelyZeroAdaptiveness = 0.0001;
    private const double ZeroAdaptivenessSubstitute = 0.01;

    /// <summary>
    /// Calculates the Codon Adaptation Index (CAI) under the Standard genetic code (NCBI table 1).
    /// See <see cref="CalculateCai(string, IReadOnlyDictionary{string, double}, GeneticCode)"/>.
    /// </summary>
    /// <param name="sequence">Coding sequence to analyze.</param>
    /// <param name="referenceRscu">Reference RSCU or w values (e.g. from highly expressed genes).</param>
    public static double CalculateCai(DnaSequence sequence, Dictionary<string, double> referenceRscu) =>
        CalculateCai(sequence, referenceRscu, GeneticCode.Standard);

    /// <summary>
    /// Calculates the CAI of a coding sequence under the given genetic code.
    /// See <see cref="CalculateCai(string, IReadOnlyDictionary{string, double}, GeneticCode)"/>.
    /// </summary>
    /// <exception cref="ArgumentNullException">An argument is null.</exception>
    public static double CalculateCai(
        DnaSequence sequence, IReadOnlyDictionary<string, double> referenceRscu, GeneticCode code)
    {
        ArgumentNullException.ThrowIfNull(sequence);
        ArgumentNullException.ThrowIfNull(referenceRscu);
        ArgumentNullException.ThrowIfNull(code);

        return CalculateCaiCore(sequence.Sequence, referenceRscu, code, excludeSingleCodonFamilies: true);
    }

    /// <summary>
    /// Calculates CAI from a raw sequence string under the Standard genetic code (NCBI table 1).
    /// See <see cref="CalculateCai(string, IReadOnlyDictionary{string, double}, GeneticCode)"/>.
    /// </summary>
    public static double CalculateCai(string sequence, Dictionary<string, double> referenceRscu) =>
        CalculateCai(sequence, referenceRscu, GeneticCode.Standard);

    /// <summary>
    /// Calculates the Codon Adaptation Index of Sharp &amp; Li (1987, Nucleic Acids Res.
    /// 15(3):1281-1295): <c>CAI = exp((1/L) Σ ln w_k)</c>, the geometric mean of the relative
    /// adaptiveness <c>w_ij = RSCU_ij / RSCU_imax = X_ij / X_imax</c> of the gene's L scored codons.
    /// </summary>
    /// <remarks>
    /// Follows the CodonW reference implementation (Peden 1999, codon_us.c <c>cai_out</c>; seqinr
    /// <c>cai</c> reproduces it):
    /// <list type="bullet">
    /// <item>Termination codons and single-codon ("non-synonymous") families — Met and Trp in
    /// table 1 — are excluded; both are genetic-code dependent (<paramref name="code"/>). Sharp &amp; Li
    /// (1987) state single-codon families should be excluded (quoted by Xia 2007, Evol. Bioinform.
    /// 3:53-58), since their w is always 1.</item>
    /// <item>A relative adaptiveness below 0.0001 (a codon absent from the reference) is replaced by
    /// 0.01 (CodonW; Bulmer 1988), so such a codon lowers CAI instead of forcing it to 0 or being
    /// silently dropped.</item>
    /// <item>Input is case-insensitive DNA or RNA (U read as T); triplets containing any other
    /// symbol are skipped without shifting the frame; a trailing partial triplet is ignored.</item>
    /// </list>
    /// <paramref name="referenceRscu"/> may hold RSCU values or w values (uppercase DNA keys):
    /// each value is divided by the maximum of its synonymous family, so a w table (family
    /// maximum 1) is used unchanged. A missing key counts as 0. A family whose values are all 0
    /// carries no reference information and its codons are not scored. Codons that are
    /// context-dependent stops in NCBI tables 27/28/31 are scored in the family of the amino acid
    /// they encode (<see cref="GeneticCode.CodonTable"/>, as for RSCU).
    /// Returns 0 when no codon is scored (empty input, only stops/Met/Trp).
    /// </remarks>
    /// <param name="sequence">Coding sequence (frame 0); null/empty returns 0.</param>
    /// <param name="referenceRscu">Reference RSCU or w values keyed by DNA codon.</param>
    /// <param name="code">Genetic code defining stops and synonymous families.</param>
    /// <exception cref="ArgumentNullException"><paramref name="referenceRscu"/> or <paramref name="code"/> is null.</exception>
    /// <exception cref="ArgumentOutOfRangeException">A reference value used is negative or not finite.</exception>
    public static double CalculateCai(
        string sequence, IReadOnlyDictionary<string, double> referenceRscu, GeneticCode code)
    {
        return CalculateCai(sequence, referenceRscu, code, excludeSingleCodonFamilies: true);
    }

    /// <summary>
    /// Canonical CAI core shared with <c>CodonOptimizer.CalculateCAI</c>.
    /// <paramref name="excludeSingleCodonFamilies"/> = false scores single-codon families with
    /// w = 1 (EMBOSS <c>ajCodCalcCaiSeq</c> convention), a non-Sharp &amp; Li opt-in.
    /// </summary>
    internal static double CalculateCai(
        string sequence, IReadOnlyDictionary<string, double> reference, GeneticCode code,
        bool excludeSingleCodonFamilies)
    {
        ArgumentNullException.ThrowIfNull(reference);
        ArgumentNullException.ThrowIfNull(code);
        if (string.IsNullOrEmpty(sequence))
            return 0;

        return CalculateCaiCore(sequence, reference, code, excludeSingleCodonFamilies);
    }

    private static double CalculateCaiCore(
        string seq, IReadOnlyDictionary<string, double> reference, GeneticCode code,
        bool excludeSingleCodonFamilies)
    {
        var w = RelativeAdaptiveness(reference, code, excludeSingleCodonFamilies);

        // CAI = geometric mean of w over the L scored codons, computed as
        // exp((1/L) Σ ln w_k) (Sharp & Li 1987; CodonW natural-log summation).
        double logSum = 0;
        long scored = 0;

        // Canonical in-frame codons (ambiguous triplets are null → not scored, frame kept).
        foreach (string? codon in SplitInFrameCodons(seq))
        {
            if (codon is not null && w.TryGetValue(codon, out double wk))
            {
                logSum += Math.Log(wk);
                scored++;
            }
        }

        return scored > 0 ? Math.Exp(logSum / scored) : 0;
    }

    // w_ij = X_ij / X_imax per synonymous family of the code (Sharp & Li 1987), with the CodonW
    // 0.0001 → 0.01 substitution. Only scored codons are keys (DNA spelling).
    private static Dictionary<string, double> RelativeAdaptiveness(
        IReadOnlyDictionary<string, double> reference, GeneticCode code, bool excludeSingleCodonFamilies)
    {
        var w = new Dictionary<string, double>(64);

        foreach (var family in SynonymousFamilies(code))
        {
            if (family.Key == '*') continue;
            var codons = family.ToList();
            if (excludeSingleCodonFamilies && codons.Count == 1) continue;

            double max = 0;
            foreach (var codon in codons)
                max = Math.Max(max, ReferenceValue(reference, codon));
            if (max <= 0) continue; // no reference information for this amino acid

            foreach (var codon in codons)
            {
                double value = ReferenceValue(reference, codon) / max;
                w[codon] = value < EffectivelyZeroAdaptiveness ? ZeroAdaptivenessSubstitute : value;
            }
        }

        return w;
    }

    private static double ReferenceValue(IReadOnlyDictionary<string, double> reference, string codon)
    {
        double value = reference.GetValueOrDefault(codon, 0);
        if (value < 0 || !double.IsFinite(value))
            throw new ArgumentOutOfRangeException(nameof(reference), value,
                $"Reference value for codon {codon} must be a finite non-negative number.");
        return value;
    }

    /// <summary>
    /// Reference relative-adaptiveness (w) table for E. coli very highly expressed genes,
    /// from Sharp &amp; Li (1987), Nucleic Acids Res. 15(13):1281-1295, as reproduced in
    /// Biopython's <c>SharpEcoliIndex</c> (Bio.SeqUtils.CodonUsageIndices, v1.79).
    /// Values are w = f_codon / f_max-synonym in [0,1]; the most-used codon of each amino
    /// acid is 1.0. Stop codons are not part of CAI and are listed as 0.0.
    /// Suitable as the <c>referenceRscu</c> argument of <see cref="CalculateCai(string, Dictionary{string, double})"/>:
    /// since CAI rescales each value by its family maximum, passing w (max 1.0) reproduces w.
    /// </summary>
    public static Dictionary<string, double> EColiOptimalCodons => new()
    {
        ["TTT"] = 0.296, ["TTC"] = 1.000, ["TTA"] = 0.020, ["TTG"] = 0.020,
        ["CTT"] = 0.042, ["CTC"] = 0.037, ["CTA"] = 0.007, ["CTG"] = 1.000,
        ["ATT"] = 0.185, ["ATC"] = 1.000, ["ATA"] = 0.003, ["ATG"] = 1.000,
        ["GTT"] = 1.000, ["GTC"] = 0.066, ["GTA"] = 0.495, ["GTG"] = 0.221,
        ["TCT"] = 1.000, ["TCC"] = 0.744, ["TCA"] = 0.077, ["TCG"] = 0.017,
        ["CCT"] = 0.070, ["CCC"] = 0.012, ["CCA"] = 0.135, ["CCG"] = 1.000,
        ["ACT"] = 0.965, ["ACC"] = 1.000, ["ACA"] = 0.076, ["ACG"] = 0.099,
        ["GCT"] = 1.000, ["GCC"] = 0.122, ["GCA"] = 0.586, ["GCG"] = 0.424,
        ["TAT"] = 0.239, ["TAC"] = 1.000, ["TAA"] = 0.000, ["TAG"] = 0.000,
        ["CAT"] = 0.291, ["CAC"] = 1.000, ["CAA"] = 0.124, ["CAG"] = 1.000,
        ["AAT"] = 0.051, ["AAC"] = 1.000, ["AAA"] = 1.000, ["AAG"] = 0.253,
        ["GAT"] = 0.434, ["GAC"] = 1.000, ["GAA"] = 1.000, ["GAG"] = 0.259,
        ["TGT"] = 0.500, ["TGC"] = 1.000, ["TGA"] = 0.000, ["TGG"] = 1.000,
        ["CGT"] = 1.000, ["CGC"] = 0.356, ["CGA"] = 0.004, ["CGG"] = 0.004,
        ["AGT"] = 0.085, ["AGC"] = 0.410, ["AGA"] = 0.004, ["AGG"] = 0.002,
        ["GGT"] = 1.000, ["GGC"] = 0.724, ["GGA"] = 0.010, ["GGG"] = 0.019
    };

    /// <summary>
    /// Reference RSCU table for <i>Homo sapiens</i>, derived from the Kazusa codon-usage
    /// database (Nakamura, Gojobori &amp; Ikemura 2000, Nucleic Acids Res. 28(1):292; species
    /// Homo sapiens [gbpri], 93,487 CDS / 40,662,582 codons, accessed 2026-06-13).
    /// RSCU_j = n·x_j / Σ_k x_k over the n synonymous codons of each amino acid, computed
    /// from the published per-thousand frequencies (Sharp, Tuohy &amp; Mosurski 1986).
    /// Single-codon families (Met, Trp) are 1.0; the '*' column holds the RSCU of the three
    /// stop codons treated as one family (not used by CAI).
    /// </summary>
    public static Dictionary<string, double> HumanOptimalCodons => new()
    {
        ["TTT"] = 0.9288, ["TTC"] = 1.0712, ["TTA"] = 0.4611, ["TTG"] = 0.7725,
        ["CTT"] = 0.7904, ["CTC"] = 1.1737, ["CTA"] = 0.4311, ["CTG"] = 2.3713,
        ["ATT"] = 1.0835, ["ATC"] = 1.4086, ["ATA"] = 0.5079, ["ATG"] = 1.0000,
        ["GTT"] = 0.7249, ["GTC"] = 0.9555, ["GTA"] = 0.4679, ["GTG"] = 1.8517,
        ["TCT"] = 1.1245, ["TCC"] = 1.3095, ["TCA"] = 0.9026, ["TCG"] = 0.3255,
        ["CCT"] = 1.1457, ["CCC"] = 1.2962, ["CCA"] = 1.1064, ["CCG"] = 0.4517,
        ["ACT"] = 0.9850, ["ACC"] = 1.4211, ["ACA"] = 1.1353, ["ACG"] = 0.4586,
        ["GCT"] = 1.0620, ["GCC"] = 1.5988, ["GCA"] = 0.9120, ["GCG"] = 0.4271,
        ["TAT"] = 0.8873, ["TAC"] = 1.1127, ["TAA"] = 0.8824, ["TAG"] = 0.7059,
        ["CAT"] = 0.8385, ["CAC"] = 1.1615, ["CAA"] = 0.5290, ["CAG"] = 1.4710,
        ["AAT"] = 0.9418, ["AAC"] = 1.0582, ["AAA"] = 0.8668, ["AAG"] = 1.1332,
        ["GAT"] = 0.9296, ["GAC"] = 1.0704, ["GAA"] = 0.8455, ["GAG"] = 1.1545,
        ["TGT"] = 0.9138, ["TGC"] = 1.0862, ["TGA"] = 1.4118, ["TGG"] = 1.0000,
        ["CGT"] = 0.4762, ["CGC"] = 1.1005, ["CGA"] = 0.6561, ["CGG"] = 1.2063,
        ["AGT"] = 0.8952, ["AGC"] = 1.4427, ["AGA"] = 1.2910, ["AGG"] = 1.2698,
        ["GGT"] = 0.6545, ["GGC"] = 1.3455, ["GGA"] = 1.0000, ["GGG"] = 1.0000
    };

    #endregion

    #region Effective Number of Codons (ENC)

    // CodonW 1.4.4 enc_out (Peden 1999, codon_us.c): an amino acid enters its class average
    // only when its homozygosity estimate exceeds this threshold ("if (bb > 0.0000001)"), i.e.
    // F̂ = 0 (every observed codon used once) is treated as not estimable, like n ≤ 1.
    private const double MinEstimableHomozygosity = 0.0000001;

    /// <summary>
    /// Calculates the Effective Number of Codons (ENC / Nc) under the Standard genetic code
    /// (NCBI table 1). See <see cref="CalculateEnc(string, GeneticCode)"/>.
    /// </summary>
    /// <exception cref="ArgumentNullException"><paramref name="sequence"/> is null.</exception>
    public static double CalculateEnc(DnaSequence sequence) =>
        CalculateEnc(sequence, GeneticCode.Standard);

    /// <summary>
    /// Calculates the Effective Number of Codons under the given genetic code.
    /// See <see cref="CalculateEnc(string, GeneticCode)"/>.
    /// </summary>
    /// <exception cref="ArgumentNullException">An argument is null.</exception>
    public static double CalculateEnc(DnaSequence sequence, GeneticCode code)
    {
        ArgumentNullException.ThrowIfNull(sequence);
        ArgumentNullException.ThrowIfNull(code);
        return CalculateEncCore(CountCodonsCore(sequence.Sequence), code);
    }

    /// <summary>
    /// Calculates ENC from a raw sequence string under the Standard genetic code (NCBI table 1).
    /// See <see cref="CalculateEnc(string, GeneticCode)"/>.
    /// </summary>
    public static double CalculateEnc(string sequence) =>
        CalculateEnc(sequence, GeneticCode.Standard);

    /// <summary>
    /// Calculates the Effective Number of Codons of Wright (1990, Gene 87:23-29):
    /// <c>Nc = K₁ + Σ_z K_z / F̄_z</c>, where K_z is the number of amino acids with z
    /// synonymous codons in <paramref name="code"/> (table 1: Nc = 2 + 9/F̄₂ + 1/F̄₃ + 5/F̄₄ + 3/F̄₆),
    /// F̄_z is the mean over the class of the codon homozygosity
    /// <c>F̂ = (n·Σ p_i² − 1)/(n − 1)</c> (p_i = n_i/n, n codons of that amino acid).
    /// Nc ranges from 20 (one codon per amino acid) to the number of sense codons (61 in table 1).
    /// </summary>
    /// <remarks>
    /// Follows the CodonW reference implementation (Peden 1999, codon_us.c <c>enc_out</c>):
    /// <list type="bullet">
    /// <item>Synonymous classes are taken from <paramref name="code"/> (CodonW <c>-enc</c> honours
    /// <c>-code</c>); stop codons are excluded. Codons that are context-dependent stops in NCBI
    /// tables 27/28/31 belong to the amino acid they encode (<see cref="GeneticCode.CodonTable"/>,
    /// as for RSCU/CAI).</item>
    /// <item>An amino acid with n ≤ 1, or with F̂ = 0 (every observed codon used once), is not
    /// estimable and is left out of its class average (Wright 1990 Eq. 4).</item>
    /// <item>If the single 3-fold amino acid (Ile) is not estimable, F̄₃ = (F̄₂ + F̄₄)/2
    /// (Wright 1990).</item>
    /// <item>If any other synonymous class has no estimable amino acid, Nc is not calculated
    /// ("the gene is either too short or has extremely skewed amino acid usage", Wright 1990 /
    /// CodonW) and this method returns 0 — never a value in the valid range.</item>
    /// <item>Values above the number of sense codons of the code (61 for table 1, CodonW's cap)
    /// are re-adjusted down to it (Wright 1990).</item>
    /// <item>Input is case-insensitive DNA or RNA (U read as T); triplets containing any other
    /// symbol are skipped without shifting the frame; a trailing partial triplet is ignored.</item>
    /// </list>
    /// </remarks>
    /// <param name="sequence">Coding sequence (frame 0); null/empty returns 0.</param>
    /// <param name="code">Genetic code defining the synonymous classes.</param>
    /// <returns>Nc, or 0 when Nc cannot be calculated.</returns>
    /// <exception cref="ArgumentNullException"><paramref name="code"/> is null.</exception>
    public static double CalculateEnc(string sequence, GeneticCode code)
    {
        ArgumentNullException.ThrowIfNull(code);
        if (string.IsNullOrEmpty(sequence))
            return 0;

        return CalculateEncCore(CountCodonsCore(sequence), code);
    }

    private static double CalculateEncCore(IReadOnlyDictionary<string, int> counts, GeneticCode code)
    {
        // Per degeneracy class z: K_z (amino acids in the class) and Σ F̂ / number of estimable
        // amino acids (CodonW fold[z], totb[z], numaa[z]).
        var aminoAcidsInClass = new SortedDictionary<int, int>();
        var sumF = new Dictionary<int, double>();
        var estimable = new Dictionary<int, int>();
        int senseCodons = 0;

        foreach (var family in SynonymousFamilies(code))
        {
            if (family.Key == '*') continue; // termination codons are not amino acids

            int z = family.Count();
            senseCodons += z;
            aminoAcidsInClass[z] = aminoAcidsInClass.GetValueOrDefault(z) + 1;

            long n = 0;
            foreach (var codon in family)
                n += counts.GetValueOrDefault(codon, 0);
            if (n <= 1) continue; // F̂ undefined (denominator n − 1)

            // Wright Eq. (1): F̂ = (n·Σ p_i² − 1)/(n − 1).
            double sumPSquared = 0;
            foreach (var codon in family)
            {
                double p = (double)counts.GetValueOrDefault(codon, 0) / n;
                sumPSquared += p * p;
            }
            double f = (n * sumPSquared - 1) / (n - 1);
            if (f <= MinEstimableHomozygosity) continue;

            sumF[z] = sumF.GetValueOrDefault(z) + f;
            estimable[z] = estimable.GetValueOrDefault(z) + 1;
        }

        double ClassMean(int z) => sumF[z] / estimable[z];

        double enc = aminoAcidsInClass.GetValueOrDefault(1); // single-codon amino acids (F = 1)
        foreach (var (z, aminoAcids) in aminoAcidsInClass)
        {
            if (z == 1) continue;

            double meanF;
            if (estimable.ContainsKey(z))
                meanF = ClassMean(z);
            else if (z == 3 && aminoAcids == 1 && estimable.ContainsKey(2) && estimable.ContainsKey(4))
                meanF = (ClassMean(2) + ClassMean(4)) / 2.0; // Ile absent: F̄₃ = (F̄₂ + F̄₄)/2
            else
                return 0; // empty synonymous class: Nc not calculated (CodonW "*****")

            enc += aminoAcids / meanF;
        }

        // F̂ ≤ 1 bounds Nc below by the number of amino acids; clamp only guards rounding.
        int aminoAcidCount = aminoAcidsInClass.Values.Sum();
        return Math.Min(senseCodons, Math.Max(aminoAcidCount, enc));
    }

    #endregion

    #region Codon Usage Statistics

    /// <summary>
    /// Gets codon usage statistics under the Standard genetic code (NCBI table 1).
    /// See <see cref="GetStatistics(string, GeneticCode)"/>.
    /// </summary>
    /// <exception cref="ArgumentNullException"><paramref name="sequence"/> is null.</exception>
    public static CodonUsageStatistics GetStatistics(DnaSequence sequence) =>
        GetStatistics(sequence, GeneticCode.Standard);

    /// <summary>
    /// Gets codon usage statistics under the given genetic code.
    /// See <see cref="GetStatistics(string, GeneticCode)"/>.
    /// </summary>
    /// <exception cref="ArgumentNullException">An argument is null.</exception>
    public static CodonUsageStatistics GetStatistics(DnaSequence sequence, GeneticCode code)
    {
        ArgumentNullException.ThrowIfNull(sequence);
        ArgumentNullException.ThrowIfNull(code);
        return GetStatisticsCore(CountCodonsCore(sequence.Sequence), code);
    }

    /// <summary>
    /// Gets codon usage statistics under the Standard genetic code (NCBI table 1).
    /// See <see cref="GetStatistics(string, GeneticCode)"/>.
    /// </summary>
    public static CodonUsageStatistics GetStatistics(string sequence) =>
        GetStatistics(sequence, GeneticCode.Standard);

    /// <summary>
    /// Gets codon usage statistics of a coding sequence (frame 0): codon counts, RSCU, ENC,
    /// total codons, GC at codon positions 1/2/3 and GC3s.
    /// </summary>
    /// <remarks>
    /// <list type="bullet">
    /// <item><c>TotalCodons</c>, <c>Gc1</c>/<c>Gc2</c>/<c>Gc3</c> (percent) are taken over every
    /// valid in-frame codon, termination codons included — EMBOSS <c>cusp</c>
    /// (ajcod.c <c>ajCodWrite</c>: "1st/2nd/3rd letter GC" over all 64 codons / CodonCount);
    /// <c>OverallGc</c> = (GC1+GC2+GC3)/3 is cusp's "Coding GC". (CodonW <c>-gc</c>/GC1-3
    /// exclude stop codons, so they differ from these values exactly when the gene contains
    /// stop codons.)</item>
    /// <item><c>Gc3s</c> (percent) = G+C at the third position of synonymous codons, i.e. codons
    /// whose amino acid has more than one codon in <paramref name="code"/>, termination codons
    /// excluded ("excluding Met, Trp and termination codons" in the Standard code; Peden 1999
    /// §1.8.2.1.3; CodonW 1.4.4 <c>gc_out</c>, genetic-code dependent via <c>how_synon</c>).
    /// CodonW reports the same quantity as a fraction. 0 when there is no synonymous codon.</item>
    /// <item><c>Rscu</c> and <c>Enc</c> use <paramref name="code"/>; see
    /// <see cref="CalculateRscu(IReadOnlyDictionary{string, int}, GeneticCode)"/> and
    /// <see cref="CalculateEnc(string, GeneticCode)"/>.</item>
    /// <item>Input is case-insensitive DNA or RNA (U read as T); triplets containing any other
    /// symbol are skipped without shifting the frame; a trailing partial triplet is ignored.
    /// Null/empty input returns all-zero statistics with empty tables.</item>
    /// </list>
    /// </remarks>
    /// <exception cref="ArgumentNullException"><paramref name="code"/> is null.</exception>
    public static CodonUsageStatistics GetStatistics(string sequence, GeneticCode code)
    {
        ArgumentNullException.ThrowIfNull(code);
        if (string.IsNullOrEmpty(sequence))
            return new CodonUsageStatistics(
                new Dictionary<string, int>(),
                new Dictionary<string, double>(),
                0, 0, 0, 0, 0, 0);

        return GetStatisticsCore(CountCodonsCore(sequence), code);
    }

    private static CodonUsageStatistics GetStatisticsCore(Dictionary<string, int> counts, GeneticCode code)
    {
        var rscu = CalculateRscu(counts, code);
        double enc = CalculateEncCore(counts, code);

        // Synonymous codons of this code: members of a sense family with more than one codon
        // (CodonW how_synon: ds[codon] > 1 and not a stop).
        var synonymous = new HashSet<string>(
            SynonymousFamilies(code)
                .Where(family => family.Key != '*' && family.Count() > 1)
                .SelectMany(family => family));

        long totalCodons = 0, gc1 = 0, gc2 = 0, gc3 = 0, synonymousCodons = 0, gc3s = 0;
        foreach (var (codon, n) in counts)
        {
            totalCodons += n;
            if (IsGC(codon[0])) gc1 += n;
            if (IsGC(codon[1])) gc2 += n;
            if (IsGC(codon[2])) gc3 += n;

            if (synonymous.Contains(codon))
            {
                synonymousCodons += n;
                if (IsGC(codon[2])) gc3s += n;
            }
        }

        static double Percent(long part, long whole) => whole > 0 ? 100.0 * part / whole : 0;

        return new CodonUsageStatistics(
            CodonCounts: counts,
            Rscu: rscu,
            Enc: enc,
            TotalCodons: (int)totalCodons,
            Gc1: Percent(gc1, totalCodons),
            Gc2: Percent(gc2, totalCodons),
            Gc3: Percent(gc3, totalCodons),
            Gc3s: Percent(gc3s, synonymousCodons));
    }

    private static bool IsGC(char c) => c is 'G' or 'C';

    #endregion
}

/// <summary>
/// Comprehensive codon usage statistics.
/// </summary>
public readonly record struct CodonUsageStatistics(
    IReadOnlyDictionary<string, int> CodonCounts,
    IReadOnlyDictionary<string, double> Rscu,
    double Enc,
    int TotalCodons,
    double Gc1,
    double Gc2,
    double Gc3,
    double Gc3s)
{
    /// <summary>
    /// Gets the overall GC content of the coding sequence.
    /// </summary>
    public double OverallGc => (Gc1 + Gc2 + Gc3) / 3;
}
