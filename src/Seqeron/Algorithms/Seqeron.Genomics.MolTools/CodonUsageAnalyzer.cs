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

    /// <summary>
    /// Calculates RSCU of a coding sequence under the given genetic code with a pseudocount
    /// added to every codon count. See
    /// <see cref="CalculateRscu(IReadOnlyDictionary{string, int}, GeneticCode, double)"/>.
    /// </summary>
    /// <exception cref="ArgumentNullException"><paramref name="sequence"/> or <paramref name="code"/> is null.</exception>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="pseudocount"/> is negative, NaN or infinite.</exception>
    public static Dictionary<string, double> CalculateRscu(DnaSequence sequence, GeneticCode code, double pseudocount)
    {
        ArgumentNullException.ThrowIfNull(sequence);
        ArgumentNullException.ThrowIfNull(code);
        ValidatePseudocount(pseudocount);
        return CalculateRscu(CountCodonsCore(sequence.Sequence), code, pseudocount);
    }

    /// <summary>
    /// Calculates RSCU from a raw sequence string (DNA or RNA, case-insensitive; counted as by
    /// <see cref="CountCodons(string)"/>) under the given genetic code with a pseudocount added to
    /// every codon count. Null/empty returns an empty dictionary. See
    /// <see cref="CalculateRscu(IReadOnlyDictionary{string, int}, GeneticCode, double)"/>.
    /// </summary>
    /// <exception cref="ArgumentNullException"><paramref name="code"/> is null.</exception>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="pseudocount"/> is negative, NaN or infinite.</exception>
    public static Dictionary<string, double> CalculateRscu(string sequence, GeneticCode code, double pseudocount)
    {
        ArgumentNullException.ThrowIfNull(code);
        ValidatePseudocount(pseudocount);
        if (string.IsNullOrEmpty(sequence))
            return new Dictionary<string, double>();

        return CalculateRscu(CountCodonsCore(sequence), code, pseudocount);
    }

    /// <summary>
    /// Calculates RSCU from codon counts with a pseudocount <c>c</c> added to every codon count
    /// before the within-family normalisation:
    /// RSCU_j = n·(x_j + c) / Σ_k (x_k + c) = n·(x_j + c) / (Σ_k x_k + n·c).
    /// </summary>
    /// <remarks>
    /// The pseudocount is applied exactly as in codonbias 0.5.0
    /// (<c>scores.RelativeSynonymousCodonUsage(pseudocount=c)</c>: <c>counts + pseudocount</c>
    /// normalised per amino acid, divided by the uniform 1/n) and in the R package cubar
    /// (<c>est_rscu(pseudo_cnt = c, level = "amino_acid")</c>: <c>prop = (cts + c)/sum(cts + c)</c>,
    /// <c>rscu = prop/mean(prop)</c>). Families are whole amino acids plus the termination family
    /// ('*'), as in <see cref="CalculateRscu(IReadOnlyDictionary{string, int}, GeneticCode)"/>
    /// (codonbias <c>ignore_stop=False</c>, cubar <c>incl_stop = TRUE</c>). With <c>c &gt; 0</c>
    /// a family absent from the gene is uniform (every member 1) instead of CodonW's 0.
    /// <c>c = 0</c> (the default overload) is the CodonW behaviour, bit for bit.
    /// </remarks>
    /// <param name="codonCounts">Codon counts keyed by uppercase DNA codon; other keys are ignored.</param>
    /// <param name="code">Genetic code defining the synonymous families.</param>
    /// <param name="pseudocount">Non-negative pseudocount added to every codon count (0 = none).</param>
    /// <returns>RSCU for each of the 64 codons (uppercase DNA spelling).</returns>
    /// <exception cref="ArgumentNullException">An argument is null.</exception>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="pseudocount"/> is negative, NaN or infinite.</exception>
    public static Dictionary<string, double> CalculateRscu(
        IReadOnlyDictionary<string, int> codonCounts, GeneticCode code, double pseudocount) =>
        CalculateRscu(codonCounts, code, pseudocount, RscuFamilyLevel.AminoAcid);

    /// <summary>
    /// Calculates RSCU of a coding sequence with a pseudocount, grouping synonymous codons at the
    /// chosen <paramref name="level"/>. See
    /// <see cref="CalculateRscu(IReadOnlyDictionary{string, int}, GeneticCode, double, RscuFamilyLevel)"/>.
    /// </summary>
    /// <exception cref="ArgumentNullException"><paramref name="sequence"/> or <paramref name="code"/> is null.</exception>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="pseudocount"/> is negative, NaN or infinite,
    /// or <paramref name="level"/> is undefined.</exception>
    public static Dictionary<string, double> CalculateRscu(DnaSequence sequence, GeneticCode code, double pseudocount,
        RscuFamilyLevel level)
    {
        ArgumentNullException.ThrowIfNull(sequence);
        ArgumentNullException.ThrowIfNull(code);
        ValidatePseudocount(pseudocount);
        ValidateFamilyLevel(level);
        return CalculateRscu(CountCodonsCore(sequence.Sequence), code, pseudocount, level);
    }

    /// <summary>
    /// Calculates RSCU from a raw sequence string (DNA or RNA, case-insensitive; counted as by
    /// <see cref="CountCodons(string)"/>) with a pseudocount, grouping synonymous codons at the chosen
    /// <paramref name="level"/>. Null/empty returns an empty dictionary. See
    /// <see cref="CalculateRscu(IReadOnlyDictionary{string, int}, GeneticCode, double, RscuFamilyLevel)"/>.
    /// </summary>
    /// <exception cref="ArgumentNullException"><paramref name="code"/> is null.</exception>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="pseudocount"/> is negative, NaN or infinite,
    /// or <paramref name="level"/> is undefined.</exception>
    public static Dictionary<string, double> CalculateRscu(string sequence, GeneticCode code, double pseudocount,
        RscuFamilyLevel level)
    {
        ArgumentNullException.ThrowIfNull(code);
        ValidatePseudocount(pseudocount);
        ValidateFamilyLevel(level);
        if (string.IsNullOrEmpty(sequence))
            return new Dictionary<string, double>();

        return CalculateRscu(CountCodonsCore(sequence), code, pseudocount, level);
    }

    /// <summary>
    /// Calculates RSCU from codon counts with a pseudocount <c>c</c>, grouping synonymous codons either by
    /// amino acid (<see cref="RscuFamilyLevel.AminoAcid"/>, the CodonW / Sharp &amp; Li family, as
    /// <see cref="CalculateRscu(IReadOnlyDictionary{string, int}, GeneticCode, double)"/>) or by codon
    /// sub-family (<see cref="RscuFamilyLevel.Subfamily"/>): RSCU_j = n·(x_j + c) / Σ_k (x_k + c) over the
    /// n codons of j's group.
    /// </summary>
    /// <remarks>
    /// <see cref="RscuFamilyLevel.Subfamily"/> reproduces the R package cubar (Liu et al.,
    /// <c>est_rscu(level = "subfam")</c>, its default): <c>get_codon_table()</c> labels every codon
    /// <c>subfam = paste(amino_acid, substr(codon, 1, 2), sep = "_")</c>, so an amino acid whose codons
    /// differ in the first two bases splits into sub-families (table 1: Leu → Leu_TT {TTA, TTG} and
    /// Leu_CT {CTN}; Ser → Ser_TC {TCN} and Ser_AG {AGT, AGC}; Arg → Arg_CG {CGN} and Arg_AG {AGA, AGG};
    /// the termination codons → *_TA {TAA, TAG} and *_TG {TGA}); within each group
    /// <c>prop = (cts + c)/sum(cts + c)</c> and <c>rscu = prop/mean(prop)</c>. Every amino acid of every
    /// genetic code is split by the same rule (e.g. table 2: Ser_TC, Ser_AG; Arg_CG; *_TA, *_AG). The
    /// termination sub-families are reported, as cubar <c>incl_stop = TRUE</c> (all 64 codons, matching
    /// the amino-acid level). With <c>c = 0</c> a group absent from the gene is reported as 0 for every
    /// member (the CodonW convention of this class; cubar returns NaN = 0/0).
    /// </remarks>
    /// <param name="codonCounts">Codon counts keyed by uppercase DNA codon; other keys are ignored.</param>
    /// <param name="code">Genetic code defining the synonymous families.</param>
    /// <param name="pseudocount">Non-negative pseudocount added to every codon count (0 = none).</param>
    /// <param name="level">Synonymous grouping (amino acid or cubar sub-family).</param>
    /// <returns>RSCU for each of the 64 codons (uppercase DNA spelling).</returns>
    /// <exception cref="ArgumentNullException">An argument is null.</exception>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="pseudocount"/> is negative, NaN or infinite,
    /// or <paramref name="level"/> is undefined.</exception>
    public static Dictionary<string, double> CalculateRscu(
        IReadOnlyDictionary<string, int> codonCounts, GeneticCode code, double pseudocount, RscuFamilyLevel level)
    {
        ArgumentNullException.ThrowIfNull(codonCounts);
        ArgumentNullException.ThrowIfNull(code);
        ValidatePseudocount(pseudocount);
        ValidateFamilyLevel(level);
        if (pseudocount == 0 && level == RscuFamilyLevel.AminoAcid)
            return CalculateRscu(codonCounts, code);

        var families = level == RscuFamilyLevel.AminoAcid
            ? SynonymousFamilies(code).Select(g => g.ToList())
            : SynonymousSubfamilies(code);
        var rscu = new Dictionary<string, double>(64);

        foreach (var family in families)
        {
            int familySize = family.Count;
            long familyTotal = 0;
            foreach (var codon in family)
                familyTotal += codonCounts.GetValueOrDefault(codon, 0);

            // Σ(x + c) = Σx + n·c; with c = 0 an absent group is 0 for every member (CodonW).
            double denominator = familyTotal + familySize * pseudocount;
            foreach (var codon in family)
                rscu[codon] = denominator > 0
                    ? familySize * (codonCounts.GetValueOrDefault(codon, 0) + pseudocount) / denominator
                    : 0.0;
        }

        return rscu;
    }

    // cubar get_codon_table(): subfam = amino acid + first two codon bases (DNA spelling), NCBI codon order.
    private static IEnumerable<List<string>> SynonymousSubfamilies(GeneticCode code) =>
        code.CodonTable
            .GroupBy(kv => (kv.Value, kv.Key.Substring(0, 2)), kv => kv.Key.Replace('U', 'T'))
            .Select(g => g.ToList());

    private static void ValidateFamilyLevel(RscuFamilyLevel level)
    {
        if (level is not (RscuFamilyLevel.AminoAcid or RscuFamilyLevel.Subfamily))
            throw new ArgumentOutOfRangeException(nameof(level), level, "Unknown RSCU family level.");
    }

    private static void ValidatePseudocount(double pseudocount)
    {
        if (!double.IsFinite(pseudocount) || pseudocount < 0)
            throw new ArgumentOutOfRangeException(nameof(pseudocount), pseudocount,
                "Pseudocount must be a finite non-negative number.");
    }

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
    /// <returns>Nc, or 0 when Nc cannot be calculated (use
    /// <see cref="TryCalculateEnc(string, GeneticCode, out double)"/> to get an explicit
    /// not-calculable result).</returns>
    /// <exception cref="ArgumentNullException"><paramref name="code"/> is null.</exception>
    public static double CalculateEnc(string sequence, GeneticCode code)
    {
        TryCalculateEnc(sequence, code, out double enc);
        return enc;
    }

    /// <summary>
    /// Tries to calculate the Effective Number of Codons under the Standard genetic code
    /// (NCBI table 1). See <see cref="TryCalculateEnc(string, GeneticCode, out double)"/>.
    /// </summary>
    /// <exception cref="ArgumentNullException"><paramref name="sequence"/> is null.</exception>
    public static bool TryCalculateEnc(DnaSequence sequence, out double enc) =>
        TryCalculateEnc(sequence, GeneticCode.Standard, out enc);

    /// <summary>
    /// Tries to calculate the Effective Number of Codons under the given genetic code.
    /// See <see cref="TryCalculateEnc(string, GeneticCode, out double)"/>.
    /// </summary>
    /// <exception cref="ArgumentNullException">An argument is null.</exception>
    public static bool TryCalculateEnc(DnaSequence sequence, GeneticCode code, out double enc)
    {
        ArgumentNullException.ThrowIfNull(sequence);
        ArgumentNullException.ThrowIfNull(code);
        return TryCalculateEncCore(CountCodonsCore(sequence.Sequence), code, out enc);
    }

    /// <summary>
    /// Tries to calculate ENC from a raw sequence string under the Standard genetic code
    /// (NCBI table 1). See <see cref="TryCalculateEnc(string, GeneticCode, out double)"/>.
    /// </summary>
    public static bool TryCalculateEnc(string sequence, out double enc) =>
        TryCalculateEnc(sequence, GeneticCode.Standard, out enc);

    /// <summary>
    /// Tries to calculate the Effective Number of Codons (Wright 1990; CodonW 1.4.4
    /// <c>enc_out</c>) with an explicit not-calculable result instead of the 0 sentinel of
    /// <see cref="CalculateEnc(string, GeneticCode)"/>. Same definition and conventions as
    /// <see cref="CalculateEnc(string, GeneticCode)"/>.
    /// </summary>
    /// <remarks>
    /// Returns <c>false</c> exactly where CodonW prints <c>*****</c> ("Nc was not calculated":
    /// a synonymous class other than the single-member Ile class has no estimable amino acid —
    /// the gene is too short or its amino-acid usage too skewed, Wright 1990) and for
    /// null/empty input or input without any valid codon.
    /// </remarks>
    /// <param name="sequence">Coding sequence (frame 0).</param>
    /// <param name="code">Genetic code defining the synonymous classes.</param>
    /// <param name="enc">Nc when calculable; otherwise 0.</param>
    /// <returns><c>true</c> when Nc could be calculated.</returns>
    /// <exception cref="ArgumentNullException"><paramref name="code"/> is null.</exception>
    public static bool TryCalculateEnc(string sequence, GeneticCode code, out double enc)
    {
        ArgumentNullException.ThrowIfNull(code);
        if (string.IsNullOrEmpty(sequence))
        {
            enc = 0;
            return false;
        }

        return TryCalculateEncCore(CountCodonsCore(sequence), code, out enc);
    }

    /// <summary>
    /// Calculates the Effective Number of Codons with an explicit estimator and options.
    /// See <see cref="TryCalculateEnc(string, GeneticCode, EncOptions, out double)"/>.
    /// </summary>
    /// <returns>Nc, or 0 when Nc cannot be calculated.</returns>
    /// <exception cref="ArgumentNullException"><paramref name="code"/> or <paramref name="options"/> is null.</exception>
    /// <exception cref="ArgumentException"><paramref name="options"/> is inconsistent (see <see cref="EncOptions"/>).</exception>
    public static double CalculateEnc(string sequence, GeneticCode code, EncOptions options)
    {
        TryCalculateEnc(sequence, code, options, out double enc);
        return enc;
    }

    /// <summary>
    /// Calculates the Effective Number of Codons with an explicit estimator and options.
    /// See <see cref="TryCalculateEnc(string, GeneticCode, EncOptions, out double)"/>.
    /// </summary>
    /// <exception cref="ArgumentNullException">An argument is null.</exception>
    /// <exception cref="ArgumentException"><paramref name="options"/> is inconsistent (see <see cref="EncOptions"/>).</exception>
    public static double CalculateEnc(DnaSequence sequence, GeneticCode code, EncOptions options)
    {
        ArgumentNullException.ThrowIfNull(sequence);
        return CalculateEnc(sequence.Sequence, code, options);
    }

    /// <summary>
    /// Tries to calculate the Effective Number of Codons with an explicit estimator and options.
    /// See <see cref="TryCalculateEnc(string, GeneticCode, EncOptions, out double)"/>.
    /// </summary>
    /// <exception cref="ArgumentNullException">An argument is null.</exception>
    /// <exception cref="ArgumentException"><paramref name="options"/> is inconsistent (see <see cref="EncOptions"/>).</exception>
    public static bool TryCalculateEnc(DnaSequence sequence, GeneticCode code, EncOptions options, out double enc)
    {
        ArgumentNullException.ThrowIfNull(sequence);
        return TryCalculateEnc(sequence.Sequence, code, options, out enc);
    }

    /// <summary>
    /// Tries to calculate the Effective Number of Codons with the estimator selected by
    /// <see cref="EncOptions.Method"/>:
    /// <list type="bullet">
    /// <item><see cref="EncMethod.Wright1990"/> (default): identical to
    /// <see cref="TryCalculateEnc(string, GeneticCode, out double)"/> (CodonW 1.4.4 <c>enc_out</c>);
    /// <see cref="EncOptions.CodonWCap"/> re-adjusts to CodonW's fixed 61 instead of the code's
    /// sense-codon count.</item>
    /// <item><see cref="EncMethod.Novembre2002"/>: Novembre (2002, Mol Biol Evol 19:1390) Nc′ —
    /// per amino acid χ² = Σ (p_i − e_i)²/e_i of the observed codon frequencies against the
    /// frequencies e_i expected from the background nucleotide composition, and
    /// <c>F = (N·χ² + N − k)/(k(N − 1))</c>; unweighted class means; no pseudocount.</item>
    /// <item><see cref="EncMethod.SunYangXia2013"/>: Sun, Yang &amp; Xia (2013, Mol Biol Evol
    /// 30:191) — codon counts with a pseudocount (default 1), <c>F = (χ² + 1)/k</c> (= Σ p_i²
    /// without background correction) and class means weighted by the amino-acid counts.</item>
    /// <item><see cref="EncMethod.EncPrimeNc"/> / <see cref="EncMethod.EncPrimeNcPrime"/>: the
    /// <c>Nc</c> / <c>Ncp</c> columns of ENCprime (Novembre's reference program,
    /// github.com/jnovembre/ENCprime, <c>src/ENCprime.c</c>, 2006 version) — F (Wright Eq. 1) or
    /// F′ = (χ² + n − k)/(k(n − 1)) only for amino acids observed n &gt; 5 times; e_i from the
    /// background base frequencies over the amino acid's variable codon positions, renormalised within
    /// the amino acid, no pseudocount (background = the gene, as <c>SeqCount -n</c> on the same FASTA,
    /// or <see cref="EncOptions.Background"/>); an amino acid with an e_i &lt; 1e-6 uses χ² = −1 as the
    /// program does; unweighted class means; a missing class is imputed with 1/k, a missing 3-fold
    /// class with (F̄₂ + F̄₄)/2 when both exist; every class mean is floored at 1/k; Nc = Σ K_z/F̄_z
    /// without any cap (so Nc′ may exceed the sense-codon count).</item>
    /// </list>
    /// </summary>
    /// <remarks>
    /// The Novembre / Sun–Yang–Xia estimators follow the reference implementation codonbias 0.5.0
    /// <c>EffectiveNumberOfCodons(k_mer=1)</c> (<c>scores.py</c> <c>_calc_score</c>/<c>_calc_F</c>):
    /// <list type="bullet">
    /// <item>Counts n_i of the sense codons of <paramref name="code"/> plus the pseudocount;
    /// N = Σ n_i per amino acid; p_i = n_i/N. Six-fold families are not split.</item>
    /// <item>Background codon frequencies e_i: product of the background base frequencies
    /// (A/C/G/T counts over the whole sequence — or <see cref="EncOptions.Background"/> — each + 1,
    /// all positions pooled), renormalised within each amino acid; uniform 1/k without
    /// background correction.</item>
    /// <item>An amino acid enters its class mean only if N &gt; 1 and 1e-6 &lt; F &lt; ∞; weighted
    /// means use N.</item>
    /// <item>A class without an estimable amino acid is imputed with F̄ = 1/z (no bias); a missing
    /// 3-fold class is then replaced by (F̄₂ + F̄₄)/2 — so these estimators are always
    /// calculable for a gene with at least one sense codon.</item>
    /// <item>Nc = Σ_z K_z/F̄_z re-adjusted down to the code's sense-codon count (no lower clamp:
    /// with background correction Nc may fall below the number of amino acids).</item>
    /// </list>
    /// Returns <c>false</c> (enc = 0) for null/empty input, for input without any sense codon of
    /// <paramref name="code"/> (codonbias would report the uninformative maximum) and — Wright
    /// only — where CodonW prints <c>*****</c>.
    /// </remarks>
    /// <param name="sequence">Coding sequence (frame 0), DNA or RNA, case-insensitive.</param>
    /// <param name="code">Genetic code defining the synonymous classes.</param>
    /// <param name="options">Estimator and its parameters.</param>
    /// <param name="enc">Nc when calculable; otherwise 0.</param>
    /// <returns><c>true</c> when Nc could be calculated.</returns>
    /// <exception cref="ArgumentNullException"><paramref name="code"/> or <paramref name="options"/> is null.</exception>
    /// <exception cref="ArgumentException"><paramref name="options"/> is inconsistent (see <see cref="EncOptions"/>).</exception>
    public static bool TryCalculateEnc(string sequence, GeneticCode code, EncOptions options, out double enc)
    {
        ArgumentNullException.ThrowIfNull(code);
        ArgumentNullException.ThrowIfNull(options);
        ValidateEncOptions(options);
        enc = 0;
        if (string.IsNullOrEmpty(sequence))
            return false;

        var counts = CountCodonsCore(sequence);
        if (options.Method == EncMethod.Wright1990)
            return TryCalculateEncCore(counts, code, out enc, options.CodonWCap ? CodonWEncCap : null);
        if (options.Method is EncMethod.EncPrimeNc or EncMethod.EncPrimeNcPrime)
            return TryCalculateEncPrimeCore(sequence, counts, code, options, out enc);

        return TryCalculateChiSquareEncCore(sequence, counts, code, options, out enc);
    }

    // ENCprime (Novembre 2002 reference program, github.com/jnovembre/ENCprime, src/ENCprime.c,
    // JN 2/28/06 version): calc_F estimates F and F′ only for amino acids observed more than 5 times.
    private const int EncPrimeMinAminoAcidCount = 5;

    // ENCprime test_null: an expected codon frequency below eps² = 1e-6 marks the amino acid's χ² as -1.
    private const double EncPrimeMinExpected = 1e-6;

    // ENCprime test_null: χ² sentinel of a single-codon amino acid (redun ≤ 1).
    private const double EncPrimeSingleCodonChiSquare = -2.0;

    // ENCprime test_null: χ² sentinel set when an expected codon frequency is below 1e-6.
    private const double EncPrimeBadExpectedChiSquare = -1.0;

    // ENCprime / SeqCount codon order (T, C, A, G at each position: TTT, TTC, TTA, TTG, TCT, …).
    private static int EncPrimeCodonOrder(string codon)
    {
        static int Rank(char b) => b switch { 'T' => 0, 'C' => 1, 'A' => 2, _ => 3 };
        return Rank(codon[0]) * 16 + Rank(codon[1]) * 4 + Rank(codon[2]);
    }

    // ENCprime Nc (Wright 1990) and Nc′ (Novembre 2002), reproducing src/ENCprime.c
    // calc_cod_data / test_null / calc_F / calc_Nc for one sequence (see EncMethod.EncPrimeNc).
    private static bool TryCalculateEncPrimeCore(
        string sequence, IReadOnlyDictionary<string, int> counts, GeneticCode code, EncOptions options, out double enc)
    {
        enc = 0;
        bool prime = options.Method == EncMethod.EncPrimeNcPrime;
        var families = SynonymousFamilies(code)
            .Where(f => f.Key != '*')
            .Select(f => f.OrderBy(EncPrimeCodonOrder).ToArray())
            .ToList();

        // No sense-codon guard: like the program, a gene without an amino acid observed more
        // than 5 times reports the uninformative Σ K_z·z (every class imputed with 1/z).
        // Background base frequencies f_A, f_C, f_G, f_T (SeqCount -n: base counts / length, no
        // pseudocount). Only their ratios matter: e_i is renormalised within the amino acid.
        double[]? baseFrequencies = prime ? EncPrimeBaseFrequencies(options.Background ?? sequence) : null;

        var aminoAcidsInClass = new SortedDictionary<int, int>(); // n_fold_o
        var sumF = new Dictionary<int, double>();                 // F_tot / F_tot_p
        var estimable = new Dictionary<int, int>();               // n_fold / n_fold_p

        foreach (var codons in families)
        {
            int k = codons.Length;
            aminoAcidsInClass[k] = aminoAcidsInClass.GetValueOrDefault(k) + 1;

            long n = 0;
            foreach (var codon in codons)
                n += counts.GetValueOrDefault(codon, 0);
            if (n <= EncPrimeMinAminoAcidCount)
                continue; // calc_F: F = F′ = -1 (not estimated)

            double f;
            if (!prime)
            {
                // Wright (1990) Eq. 1: F = (n·Σp² − 1)/(n − 1).
                double sumPSquared = 0;
                foreach (var codon in codons)
                {
                    double p = (double)counts.GetValueOrDefault(codon, 0) / n;
                    sumPSquared += p * p;
                }
                f = (n * sumPSquared - 1.0) / (n - 1.0);
            }
            else
            {
                double chiSquare = EncPrimeChiSquare(codons, counts, n, baseFrequencies!);
                f = (chiSquare + n - k) / (k * (n - 1.0)); // Novembre (2002) Eq. 6
            }

            sumF[k] = sumF.GetValueOrDefault(k) + f;
            estimable[k] = estimable.GetValueOrDefault(k) + 1;
        }

        // calc_Nc: class mean over estimable amino acids; a missing 3-fold class takes the
        // arithmetic mean of the 2- and 4-fold means (JN 2/28/06 fix), any other missing class 1/k;
        // the mean is then floored at 1/k. No cap and no lower clamp.
        double ClassMean(int z) => sumF[z] / estimable[z];
        double nc = 0;
        foreach (var (z, aminoAcids) in aminoAcidsInClass)
        {
            double meanF;
            if (estimable.ContainsKey(z))
                meanF = ClassMean(z);
            else if (z == 3 && estimable.ContainsKey(2) && estimable.ContainsKey(4))
                meanF = (ClassMean(2) + ClassMean(4)) * 0.5;
            else
                meanF = 1.0 / z;

            if (meanF < 1.0 / z)
                meanF = 1.0 / z;
            nc += aminoAcids / meanF;
        }

        enc = nc;
        return true;
    }

    // ENCprime test_null χ² = Σ n·(p_i − e_i)²/e_i over the codons of one amino acid in ENCprime's
    // codon order, with e_i = Π f(base) over the amino acid's variable codon positions,
    // renormalised within the amino acid (calc_cod_data). Sentinels as in the program: -2 for a
    // single-codon amino acid; an e_i < 1e-6 resets χ² to -1 and accumulation continues.
    private static double EncPrimeChiSquare(
        string[] codons, IReadOnlyDictionary<string, int> counts, long n, double[] baseFrequencies)
    {
        int k = codons.Length;
        if (k <= 1)
            return EncPrimeSingleCodonChiSquare;

        var variable = new bool[3];
        for (int pos = 0; pos < 3; pos++)
            variable[pos] = codons.Any(c => c[pos] != codons[0][pos]);

        var numerator = new double[k];
        double denominator = 0;
        for (int i = 0; i < k; i++)
        {
            double e = 1;
            for (int pos = 0; pos < 3; pos++)
                if (variable[pos])
                    e *= baseFrequencies[BaseIndex(codons[i][pos])];
            numerator[i] = e;
            denominator += e;
        }
        if (denominator == 0.0)
            denominator = 1; // calc_cod_data: "if (ex_denom == 0.0) ex_denom = 1"

        double chiSquare = 0;
        for (int i = 0; i < k; i++)
        {
            double expected = numerator[i] / denominator;
            if (expected < EncPrimeMinExpected)
            {
                chiSquare = EncPrimeBadExpectedChiSquare;
                continue;
            }
            double bias = (double)counts.GetValueOrDefault(codons[i], 0) / n - expected;
            chiSquare += bias * bias / expected * n;
        }
        return chiSquare;
    }

    // SeqCount -n: A/C/G/T counts over the whole sequence (case-insensitive; U read as T, the
    // library's RNA convention), divided by the total; other symbols ignored. All zero when the
    // sequence has no A/C/G/T (every e_i is then 0 and χ² = -1, as in ENCprime).
    private static double[] EncPrimeBaseFrequencies(string sequence)
    {
        var counts = new double[4];
        foreach (char ch in sequence)
        {
            int i = BaseIndex(char.ToUpperInvariant(ch));
            if (i >= 0) counts[i]++;
        }
        double total = counts[0] + counts[1] + counts[2] + counts[3];
        if (total > 0)
            for (int i = 0; i < 4; i++) counts[i] /= total;
        return counts;
    }

    // CodonW 1.4.4 enc_out prints every Nc above 61 as "61.00", whatever the genetic code
    // ("else if (enc_tot <= 61) … else fprintf(foutput, \"61.00%c\", sp)").
    private const int CodonWEncCap = 61;

    // codonbias 0.5.0 EffectiveNumberOfCodons._calc_score: an amino acid enters its class mean
    // only when "(N > 1) & (F > 1e-6) & np.isfinite(F)".
    private const double MinChiSquareEstimableHomozygosity = 1e-6;

    // Novembre (2002) Nc′ and Sun, Yang & Xia (2013) Nc, as implemented by codonbias 0.5.0
    // (EffectiveNumberOfCodons._calc_F / _calc_BNC / _calc_BCC / _calc_score, k_mer = 1).
    private static bool TryCalculateChiSquareEncCore(
        string sequence, IReadOnlyDictionary<string, int> counts, GeneticCode code, EncOptions options, out double enc)
    {
        enc = 0;
        bool robust = options.Method == EncMethod.SunYangXia2013;
        double pseudocount = options.EffectivePseudocount;
        bool weighted = options.EffectiveWeightedAverage;
        double[]? baseFrequencies = options.EffectiveBackgroundCorrection
            ? BackgroundBaseFrequencies(options.Background ?? sequence)
            : null;

        var families = SynonymousFamilies(code).Where(f => f.Key != '*').ToList();
        long observedSense = 0;
        foreach (var family in families)
            foreach (var codon in family)
                observedSense += counts.GetValueOrDefault(codon, 0);
        if (observedSense == 0)
            return false; // no sense codon: nothing to estimate

        var aminoAcidsInClass = new SortedDictionary<int, int>();
        var sumF = new Dictionary<int, double>();
        var sumWeight = new Dictionary<int, double>();
        int senseCodons = 0;

        foreach (var family in families)
        {
            string[] codons = family.ToArray();
            int k = codons.Length;
            senseCodons += k;
            aminoAcidsInClass[k] = aminoAcidsInClass.GetValueOrDefault(k) + 1;

            var n = new double[k];
            double total = 0;
            for (int i = 0; i < k; i++)
            {
                n[i] = counts.GetValueOrDefault(codons[i], 0) + pseudocount;
                total += n[i];
            }

            // Expected (background) codon frequencies within the family.
            var expected = new double[k];
            if (baseFrequencies is null)
            {
                for (int i = 0; i < k; i++) expected[i] = 1.0 / k;
            }
            else
            {
                double familySum = 0;
                for (int i = 0; i < k; i++)
                {
                    string c = codons[i];
                    expected[i] = baseFrequencies[BaseIndex(c[0])] * baseFrequencies[BaseIndex(c[1])] * baseFrequencies[BaseIndex(c[2])];
                    familySum += expected[i];
                }
                for (int i = 0; i < k; i++) expected[i] /= familySum;
            }

            // χ² of the observed against the expected codon frequencies (0 when N = 0).
            double chiSquare = 0;
            if (total > 0)
            {
                for (int i = 0; i < k; i++)
                {
                    double d = n[i] / total - expected[i];
                    chiSquare += d * d / expected[i];
                }
            }

            double f = robust
                ? (chiSquare + 1) / k                                // Sun, Yang & Xia 2013
                : (total * chiSquare + total - k) / (total - 1) / k; // Novembre 2002
            if (!(total > 1) || !(f > MinChiSquareEstimableHomozygosity) || !double.IsFinite(f))
                continue;

            double weight = weighted ? total : 1.0;
            sumF[k] = sumF.GetValueOrDefault(k) + f * weight;
            sumWeight[k] = sumWeight.GetValueOrDefault(k) + weight;
        }

        // Class means; a class without an estimable amino acid is imputed with 1/z, and a
        // missing 3-fold class with (F̄₂ + F̄₄)/2 (codonbias "Wright/Alexaki imputation").
        var meanF = new SortedDictionary<int, double>();
        foreach (int z in aminoAcidsInClass.Keys)
            meanF[z] = sumWeight.TryGetValue(z, out double w) && w > 0 ? sumF[z] / w : 1.0 / z;
        if (aminoAcidsInClass.ContainsKey(3) && !sumWeight.ContainsKey(3)
            && meanF.TryGetValue(2, out double f2) && meanF.TryGetValue(4, out double f4))
            meanF[3] = 0.5 * (f2 + f4);

        double nc = 0;
        foreach (var (z, aminoAcids) in aminoAcidsInClass)
            nc += aminoAcids / meanF[z];

        enc = Math.Min(senseCodons, nc);
        return true;
    }

    // Background base frequencies (A, C, G, T): counts of A/C/G/T(U) over the whole sequence,
    // case-insensitive, other symbols ignored, each + 1, normalised (codonbias _calc_BNC).
    private static double[] BackgroundBaseFrequencies(string sequence)
    {
        var counts = new double[] { 1, 1, 1, 1 };
        foreach (char ch in sequence)
        {
            int i = BaseIndex(char.ToUpperInvariant(ch));
            if (i >= 0) counts[i]++;
        }
        double total = counts[0] + counts[1] + counts[2] + counts[3];
        for (int i = 0; i < 4; i++) counts[i] /= total;
        return counts;
    }

    private static int BaseIndex(char b) => b switch
    {
        'A' => 0,
        'C' => 1,
        'G' => 2,
        'T' or 'U' => 3,
        _ => -1,
    };

    private static void ValidateEncOptions(EncOptions options)
    {
        switch (options.Method)
        {
            case EncMethod.Wright1990:
                if (options.Pseudocount is not null || options.BackgroundCorrection is not null || options.WeightedAverage is not null || options.Background is not null)
                    throw new ArgumentException(
                        "Pseudocount, BackgroundCorrection, WeightedAverage and Background apply to the Novembre2002 and SunYangXia2013 estimators only.",
                        nameof(options));
                break;
            case EncMethod.Novembre2002:
            case EncMethod.SunYangXia2013:
                if (options.CodonWCap)
                    throw new ArgumentException("CodonWCap applies to the Wright1990 estimator only (CodonW computes Wright's Nc).", nameof(options));
                if (options.Pseudocount is double pc && (!double.IsFinite(pc) || pc < 0))
                    throw new ArgumentException($"Pseudocount must be a finite non-negative number (got {pc}).", nameof(options));
                if (options.Background is not null && !options.EffectiveBackgroundCorrection)
                    throw new ArgumentException("Background requires background correction.", nameof(options));
                break;
            case EncMethod.EncPrimeNc:
            case EncMethod.EncPrimeNcPrime:
                if (options.CodonWCap || options.Pseudocount is not null || options.BackgroundCorrection is not null || options.WeightedAverage is not null)
                    throw new ArgumentException(
                        "CodonWCap, Pseudocount, BackgroundCorrection and WeightedAverage do not apply to the ENCprime estimators (the program has fixed rules).",
                        nameof(options));
                if (options.Background is not null && options.Method == EncMethod.EncPrimeNc)
                    throw new ArgumentException("Background applies to EncPrimeNcPrime only (ENCprime's Nc does not use the nucleotide composition).", nameof(options));
                break;
            default:
                throw new ArgumentException($"Unknown ENC method: {options.Method}.", nameof(options));
        }
    }

    private static double CalculateEncCore(IReadOnlyDictionary<string, int> counts, GeneticCode code)
    {
        TryCalculateEncCore(counts, code, out double enc);
        return enc;
    }

    private static bool TryCalculateEncCore(
        IReadOnlyDictionary<string, int> counts, GeneticCode code, out double enc, int? cap = null)
    {
        enc = 0;
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

        double nc = aminoAcidsInClass.GetValueOrDefault(1); // single-codon amino acids (F = 1)
        foreach (var (z, aminoAcids) in aminoAcidsInClass)
        {
            if (z == 1) continue;

            double meanF;
            if (estimable.ContainsKey(z))
                meanF = ClassMean(z);
            else if (z == 3 && aminoAcids == 1 && estimable.ContainsKey(2) && estimable.ContainsKey(4))
                meanF = (ClassMean(2) + ClassMean(4)) / 2.0; // Ile absent: F̄₃ = (F̄₂ + F̄₄)/2
            else
                return false; // empty synonymous class: Nc not calculated (CodonW "*****")

            nc += aminoAcids / meanF;
        }

        // F̂ ≤ 1 bounds Nc below by the number of amino acids; clamp only guards rounding.
        int aminoAcidCount = aminoAcidsInClass.Values.Sum();
        enc = Math.Min(cap ?? senseCodons, Math.Max(aminoAcidCount, nc));
        return true;
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
        return GetStatisticsCore(CountCodonsCore(sequence.Sequence), code, CodonGcConvention.EmbossCusp);
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
    /// stop codons; use <see cref="GetStatistics(string, GeneticCode, CodonGcConvention)"/> with
    /// <see cref="CodonGcConvention.CodonW"/> for CodonW's values.)</item>
    /// <item><c>Gc3s</c> (percent) = G+C at the third position of synonymous codons, i.e. codons
    /// whose amino acid has more than one codon in <paramref name="code"/>, termination codons
    /// excluded ("excluding Met, Trp and termination codons" in the Standard code; Peden 1999
    /// §1.8.2.1.3; CodonW 1.4.4 <c>gc_out</c>, genetic-code dependent via <c>how_synon</c>).
    /// CodonW reports the same quantity as a fraction. 0 when there is no synonymous codon.</item>
    /// <item><c>Rscu</c> and <c>Enc</c> use <paramref name="code"/> (<c>Enc</c> = 0 and
    /// <c>IsEncCalculable</c> = false when Nc is not calculable); see
    /// <see cref="CalculateRscu(IReadOnlyDictionary{string, int}, GeneticCode)"/> and
    /// <see cref="CalculateEnc(string, GeneticCode)"/>.</item>
    /// <item>Input is case-insensitive DNA or RNA (U read as T); triplets containing any other
    /// symbol are skipped without shifting the frame; a trailing partial triplet is ignored.
    /// Null/empty input returns all-zero statistics with empty tables.</item>
    /// </list>
    /// </remarks>
    /// <exception cref="ArgumentNullException"><paramref name="code"/> is null.</exception>
    public static CodonUsageStatistics GetStatistics(string sequence, GeneticCode code) =>
        GetStatistics(sequence, code, CodonGcConvention.EmbossCusp);

    /// <summary>
    /// Gets codon usage statistics under the given genetic code with the chosen convention for
    /// the per-position GC values. See <see cref="GetStatistics(string, GeneticCode, CodonGcConvention)"/>.
    /// </summary>
    /// <exception cref="ArgumentNullException">An argument is null.</exception>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="gcConvention"/> is not defined.</exception>
    public static CodonUsageStatistics GetStatistics(DnaSequence sequence, GeneticCode code, CodonGcConvention gcConvention)
    {
        ArgumentNullException.ThrowIfNull(sequence);
        ArgumentNullException.ThrowIfNull(code);
        ValidateGcConvention(gcConvention);
        return GetStatisticsCore(CountCodonsCore(sequence.Sequence), code, gcConvention);
    }

    /// <summary>
    /// Gets codon usage statistics with the chosen convention for <c>Gc1</c>/<c>Gc2</c>/<c>Gc3</c>
    /// (and hence <c>OverallGc</c>). Everything else is as in
    /// <see cref="GetStatistics(string, GeneticCode)"/>.
    /// </summary>
    /// <remarks>
    /// <list type="bullet">
    /// <item><see cref="CodonGcConvention.EmbossCusp"/> (default of the other overloads):
    /// GC1/GC2/GC3 over every valid in-frame codon, termination codons included (EMBOSS
    /// <c>cusp</c>).</item>
    /// <item><see cref="CodonGcConvention.CodonW"/>: GC1/GC2/GC3 over the sense codons only —
    /// codons that are stops in <paramref name="code"/> are left out of numerator and
    /// denominator (CodonW 1.4.4 <c>codon_us.c gc_out</c>: <c>if (pcu-&gt;ca[id] == 11) continue;</c>,
    /// denominator <c>Len_aa</c>). <c>OverallGc</c> = (GC1+GC2+GC3)/3 then equals CodonW's
    /// <c>GC</c> index (<c>-gc</c>). CodonW reports fractions; the values here are percent.
    /// A gene without sense codons gives 0 (CodonW writes no value).</item>
    /// <item><c>TotalCodons</c>, <c>Gc3s</c> (already stop-free, CodonW <c>-gc3s</c>), counts,
    /// RSCU and ENC do not depend on the convention.</item>
    /// </list>
    /// </remarks>
    /// <exception cref="ArgumentNullException"><paramref name="code"/> is null.</exception>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="gcConvention"/> is not defined.</exception>
    public static CodonUsageStatistics GetStatistics(string sequence, GeneticCode code, CodonGcConvention gcConvention)
    {
        ArgumentNullException.ThrowIfNull(code);
        ValidateGcConvention(gcConvention);
        if (string.IsNullOrEmpty(sequence))
            return new CodonUsageStatistics(
                new Dictionary<string, int>(),
                new Dictionary<string, double>(),
                0, 0, 0, 0, 0, 0);

        return GetStatisticsCore(CountCodonsCore(sequence), code, gcConvention);
    }

    private static void ValidateGcConvention(CodonGcConvention gcConvention)
    {
        if (gcConvention is not (CodonGcConvention.EmbossCusp or CodonGcConvention.CodonW))
            throw new ArgumentOutOfRangeException(nameof(gcConvention), gcConvention, "Unknown GC convention.");
    }

    private static CodonUsageStatistics GetStatisticsCore(
        Dictionary<string, int> counts, GeneticCode code, CodonGcConvention gcConvention)
    {
        var rscu = CalculateRscu(counts, code);
        double enc = CalculateEncCore(counts, code);

        // Synonymous codons of this code: members of a sense family with more than one codon
        // (CodonW how_synon: ds[codon] > 1 and not a stop).
        var synonymous = new HashSet<string>(
            SynonymousFamilies(code)
                .Where(family => family.Key != '*' && family.Count() > 1)
                .SelectMany(family => family));

        // CodonW gc_out skips termination codons of the selected code for GC1-3 (ca == 11).
        var stops = gcConvention == CodonGcConvention.CodonW
            ? new HashSet<string>(SynonymousFamilies(code).Where(f => f.Key == '*').SelectMany(f => f))
            : null;

        long totalCodons = 0, gcCodons = 0, gc1 = 0, gc2 = 0, gc3 = 0, synonymousCodons = 0, gc3s = 0;
        foreach (var (codon, n) in counts)
        {
            totalCodons += n;

            if (synonymous.Contains(codon))
            {
                synonymousCodons += n;
                if (IsGC(codon[2])) gc3s += n;
            }

            if (stops is not null && stops.Contains(codon)) continue;
            gcCodons += n;
            if (IsGC(codon[0])) gc1 += n;
            if (IsGC(codon[1])) gc2 += n;
            if (IsGC(codon[2])) gc3 += n;
        }

        static double Percent(long part, long whole) => whole > 0 ? 100.0 * part / whole : 0;

        return new CodonUsageStatistics(
            CodonCounts: counts,
            Rscu: rscu,
            Enc: enc,
            TotalCodons: (int)totalCodons,
            Gc1: Percent(gc1, gcCodons),
            Gc2: Percent(gc2, gcCodons),
            Gc3: Percent(gc3, gcCodons),
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

    /// <summary>
    /// Whether <see cref="Enc"/> is a calculated Nc. <c>false</c> when Nc could not be
    /// calculated (CodonW prints <c>*****</c>) — <see cref="Enc"/> is then 0; a calculated Nc is
    /// never below the number of amino acids of the code. See
    /// <see cref="CodonUsageAnalyzer.TryCalculateEnc(string, GeneticCode, out double)"/>.
    /// </summary>
    public bool IsEncCalculable => Enc > 0;
}

/// <summary>
/// Convention for the per-codon-position GC values of <see cref="CodonUsageStatistics"/>.
/// </summary>
public enum CodonGcConvention
{
    /// <summary>EMBOSS <c>cusp</c>: GC1/GC2/GC3 over all valid codons, termination codons included (default).</summary>
    EmbossCusp = 0,

    /// <summary>CodonW 1.4.4 <c>gc_out</c>: GC1/GC2/GC3 (and GC) over sense codons only, termination codons of the genetic code excluded.</summary>
    CodonW = 1,
}

/// <summary>
/// Estimator of the Effective Number of Codons (see
/// <see cref="CodonUsageAnalyzer.TryCalculateEnc(string, GeneticCode, EncOptions, out double)"/>).
/// </summary>
public enum EncMethod
{
    /// <summary>Wright (1990) Nc with the CodonW 1.4.4 <c>enc_out</c> conventions (default).</summary>
    Wright1990 = 0,

    /// <summary>Novembre (2002) background-corrected Nc′ (χ² against background codon frequencies; codonbias 0.5.0 <c>bg_correction=True, robust=False, pseudocount=0, mean='unweighted'</c>).</summary>
    Novembre2002 = 1,

    /// <summary>Sun, Yang &amp; Xia (2013) Nc: pseudocount, F = (χ² + 1)/k, count-weighted class means (codonbias 0.5.0 defaults <c>robust=True, pseudocount=1, mean='weighted'</c>).</summary>
    SunYangXia2013 = 2,

    /// <summary>
    /// Wright's Nc exactly as printed by ENCprime (Novembre 2002 reference program, column <c>Nc</c>):
    /// F = (nΣp² − 1)/(n − 1) only for amino acids with n &gt; 5; a missing class is imputed with 1/k
    /// (3-fold: (F̄₂ + F̄₄)/2); class means floored at 1/k; no cap.
    /// </summary>
    EncPrimeNc = 3,

    /// <summary>
    /// Novembre's Nc′ exactly as printed by ENCprime (column <c>Ncp</c>): F′ = (χ² + n − k)/(k(n − 1)) with
    /// χ² = Σ n(p_i − e_i)²/e_i against the background base composition (no pseudocount; the gene itself or
    /// <see cref="EncOptions.Background"/>); same n &gt; 5 rule, imputation, 1/k floor and no cap as <see cref="EncPrimeNc"/>.
    /// </summary>
    EncPrimeNcPrime = 4,
}

/// <summary>
/// Options of <see cref="CodonUsageAnalyzer.TryCalculateEnc(string, GeneticCode, EncOptions, out double)"/>.
/// The default instance is Wright (1990) / CodonW with the code's sense-codon cap, identical to
/// <see cref="CodonUsageAnalyzer.TryCalculateEnc(string, GeneticCode, out double)"/>.
/// </summary>
/// <remarks>
/// <see cref="Pseudocount"/>, <see cref="BackgroundCorrection"/>, <see cref="WeightedAverage"/>
/// and <see cref="Background"/> apply to <see cref="EncMethod.Novembre2002"/> and
/// <see cref="EncMethod.SunYangXia2013"/> only (null = the method's published default);
/// <see cref="CodonWCap"/> applies to <see cref="EncMethod.Wright1990"/> only; the ENCprime
/// estimators take only <see cref="Background"/> (Nc′ only). Setting an option
/// that does not apply, a negative or non-finite pseudocount, or a <see cref="Background"/>
/// without background correction throws <see cref="ArgumentException"/>.
/// </remarks>
public sealed record EncOptions
{
    /// <summary>Default options (Wright 1990, sense-codon cap).</summary>
    public static EncOptions Default { get; } = new();

    /// <summary>Estimator (default <see cref="EncMethod.Wright1990"/>).</summary>
    public EncMethod Method { get; init; } = EncMethod.Wright1990;

    /// <summary>Pseudocount added to every sense-codon count (codonbias <c>pseudocount</c>); null = 0 for Novembre 2002, 1 for Sun–Yang–Xia 2013.</summary>
    public double? Pseudocount { get; init; }

    /// <summary>Novembre (2002) background correction of the expected codon frequencies (codonbias <c>bg_correction</c>); null = true for Novembre 2002, false for Sun–Yang–Xia 2013.</summary>
    public bool? BackgroundCorrection { get; init; }

    /// <summary>Class means weighted by amino-acid counts (codonbias <c>mean='weighted'</c>); null = false for Novembre 2002, true for Sun–Yang–Xia 2013.</summary>
    public bool? WeightedAverage { get; init; }

    /// <summary>Sequence whose base composition defines the background (codonbias <c>get_score(background=…)</c>; ENCprime: the <c>SeqCount -n</c> nucleotide-composition file); null = the gene itself. Requires background correction (Novembre2002 / SunYangXia2013) or <see cref="EncMethod.EncPrimeNcPrime"/>.</summary>
    public string? Background { get; init; }

    /// <summary>Wright 1990 only: re-adjust Nc to CodonW's fixed 61 (codon_us.c <c>enc_out</c>) instead of the genetic code's sense-codon count (differs for codes with ≠ 61 sense codons, e.g. 60 for table 2, 62 for tables 3 and 9).</summary>
    public bool CodonWCap { get; init; }

    internal double EffectivePseudocount => Pseudocount ?? (Method == EncMethod.SunYangXia2013 ? 1.0 : 0.0);

    internal bool EffectiveBackgroundCorrection => BackgroundCorrection ?? Method == EncMethod.Novembre2002;

    internal bool EffectiveWeightedAverage => WeightedAverage ?? Method == EncMethod.SunYangXia2013;

}

/// <summary>
/// Synonymous-codon grouping of
/// <see cref="CodonUsageAnalyzer.CalculateRscu(IReadOnlyDictionary{string, int}, GeneticCode, double, RscuFamilyLevel)"/>.
/// </summary>
public enum RscuFamilyLevel
{
    /// <summary>One family per amino acid (and one termination family) — Sharp &amp; Li 1987, CodonW, cubar <c>level = "amino_acid"</c> (default).</summary>
    AminoAcid = 0,

    /// <summary>
    /// cubar <c>level = "subfam"</c> (cubar's default): codons grouped by amino acid and their first two bases
    /// (<c>subfam = paste(amino_acid, substr(codon, 1, 2), sep = "_")</c>), so Leu, Ser and Arg (table 1) split into a
    /// four-fold and a two-fold sub-family.
    /// </summary>
    Subfamily = 1,
}
