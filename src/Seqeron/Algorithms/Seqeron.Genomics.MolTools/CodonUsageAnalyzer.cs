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
    /// </remarks>
    public static Dictionary<string, int> CountCodons(string sequence)
    {
        if (string.IsNullOrEmpty(sequence))
            return new Dictionary<string, int>();

        return CountCodonsCore(NormalizeCodingSequence(sequence));
    }

    // Upper-cases and reads RNA U as T (CodonW ident_codon treats T/t/U/u identically;
    // EMBOSS ajBaseAlphaToBin maps U to the T bit), so DNA and RNA spellings of the same
    // coding sequence give identical codon counts and indices.
    private static string NormalizeCodingSequence(string sequence) =>
        sequence.ToUpperInvariant().Replace('U', 'T');

    private static Dictionary<string, int> CountCodonsCore(string seq)
    {
        var counts = new Dictionary<string, int>();

        for (int i = 0; i + 3 <= seq.Length; i += 3)
        {
            string codon = seq.Substring(i, 3);
            if (IsValidCodon(codon))
            {
                counts.TryGetValue(codon, out int count);
                counts[codon] = count + 1;
            }
        }

        return counts;
    }

    private static bool IsValidCodon(string codon)
    {
        return codon.Length == 3 && codon.All(c => c is 'A' or 'C' or 'G' or 'T');
    }

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

        return CalculateRscu(CountCodonsCore(NormalizeCodingSequence(sequence)), code);
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
            long familyTotal = 0;
            foreach (var codon in family)
                familyTotal += codonCounts.GetValueOrDefault(codon, 0);

            foreach (var codon in family)
            {
                // RSCU = n·x / Σx; an absent family (Σx = 0) is 0 for every member (CodonW).
                rscu[codon] = familyTotal > 0
                    ? (double)family.Count * codonCounts.GetValueOrDefault(codon, 0) / familyTotal
                    : 0.0;
            }
        }

        return rscu;
    }

    // Synonymous codon families (DNA spelling) of a genetic code, keyed by the encoded
    // amino acid ('*' = termination), in NCBI codon order.
    private static IEnumerable<List<string>> SynonymousFamilies(GeneticCode code) =>
        code.CodonTable
            .GroupBy(kv => kv.Value, kv => kv.Key.Replace('U', 'T'))
            .Select(g => g.ToList());

    #endregion

    #region CAI (Codon Adaptation Index)

    /// <summary>
    /// Calculates Codon Adaptation Index (CAI) using a reference codon table.
    /// CAI measures how well codon usage matches highly expressed genes.
    /// Range: 0-1, where 1 means optimal codon usage.
    /// </summary>
    /// <param name="sequence">Coding sequence to analyze.</param>
    /// <param name="referenceRscu">RSCU values from reference set (e.g., highly expressed genes).</param>
    public static double CalculateCai(DnaSequence sequence, Dictionary<string, double> referenceRscu)
    {
        ArgumentNullException.ThrowIfNull(sequence);
        ArgumentNullException.ThrowIfNull(referenceRscu);

        return CalculateCaiCore(sequence.Sequence, referenceRscu);
    }

    /// <summary>
    /// Calculates CAI from a raw sequence string.
    /// </summary>
    public static double CalculateCai(string sequence, Dictionary<string, double> referenceRscu)
    {
        if (string.IsNullOrEmpty(sequence))
            return 0;

        return CalculateCaiCore(NormalizeCodingSequence(sequence), referenceRscu);
    }

    private static double CalculateCaiCore(string seq, Dictionary<string, double> referenceRscu)
    {
        // Relative adaptiveness w_i = f_i / max(f_j) over the synonymous family of the
        // codon's amino acid (Sharp & Li 1987, Nucleic Acids Res. 15:1281-1295).
        // Non-synonymous codons (single-codon amino acids Met/Trp) and termination codons
        // are excluded from CAI (Sharp & Li 1987; CodonW codon-usage indices; EMBOSS cai).
        var relativeAdaptiveness = new Dictionary<string, double>();

        foreach (var aaGroup in CodonToAminoAcid.GroupBy(kv => kv.Value))
        {
            // Exclude termination codons ('*') and single-codon amino acids (Met, Trp):
            // they carry no synonymous bias and are not counted in CAI.
            if (aaGroup.Key == '*') continue;

            var synonymousCodons = aaGroup.Select(kv => kv.Key).ToList();
            if (synonymousCodons.Count == 1) continue;

            double maxRscu = synonymousCodons.Max(c => referenceRscu.GetValueOrDefault(c, 0));

            foreach (var codon in synonymousCodons)
            {
                double rscu = referenceRscu.GetValueOrDefault(codon, 0);
                relativeAdaptiveness[codon] = maxRscu > 0 ? rscu / maxRscu : 0;
            }
        }

        // CAI = geometric mean of w over the L scored codons, computed as
        // exp((1/L) Σ ln w_i) for numerical stability (Sharp & Li 1987, Eq. 2).
        double logSum = 0;
        int codonCount = 0;

        for (int i = 0; i + 3 <= seq.Length; i += 3)
        {
            string codon = seq.Substring(i, 3);
            if (IsValidCodon(codon) && relativeAdaptiveness.TryGetValue(codon, out double w) && w > 0)
            {
                logSum += Math.Log(w);
                codonCount++;
            }
        }

        return codonCount > 0 ? Math.Exp(logSum / codonCount) : 0;
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

    // --- Constants per Wright F. (1990) Gene 87(1):23–29, as reproduced verbatim in
    //     Fuglsang A. (2004) BBRC 317:957–964 (Eqs. 1–5a). ---

    // Number of synonymous-codon amino acids in each degeneracy class of the
    // standard (NCBI table 1) genetic code, used as the numerators of Wright Eq. (3).
    private const int TwoFoldAminoAcidCount = 9;   // 9 doublets (His, Gln, …)
    private const int ThreeFoldAminoAcidCount = 1; // 1 triplet  (Ile)
    private const int FourFoldAminoAcidCount = 5;  // 5 quartets (Ala, Gly, Pro, Thr, Val)
    private const int SixFoldAminoAcidCount = 3;   // 3 sextets  (Leu, Ser, Arg)

    // The two single-codon amino acids Met (ATG) and Trp (TGG) each contribute exactly
    // one effective codon; this is the constant "2" in Wright Eq. (3).
    private const double SingleCodonAminoAcidContribution = 2.0;

    // Wright Eq. (3): if Nc exceeds 61 it is re-adjusted down to 61 (the maximum number
    // of sense codons in the standard genetic code).
    private const double MaxEffectiveCodons = 61.0;

    // Structural lower bound: every degeneracy class collapsed to one codon gives Nc = 20
    // (the extreme-bias limit stated by Wright/Fuglsang).
    private const double MinEffectiveCodons = 20.0;

    /// <summary>
    /// Calculates the Effective Number of Codons (ENC / Nc) per Wright (1990).
    /// Nc ranges from 20 (extreme bias — one codon per amino acid) to 61 (no bias).
    /// </summary>
    public static double CalculateEnc(DnaSequence sequence)
    {
        ArgumentNullException.ThrowIfNull(sequence);
        return CalculateEncCore(sequence.Sequence);
    }

    /// <summary>
    /// Calculates ENC from a raw sequence string. Null/empty returns 0.
    /// </summary>
    public static double CalculateEnc(string sequence)
    {
        if (string.IsNullOrEmpty(sequence))
            return 0;

        return CalculateEncCore(NormalizeCodingSequence(sequence));
    }

    private static double CalculateEncCore(string seq)
    {
        var counts = CountCodonsCore(seq);

        // Wright Eq. (1): per-amino-acid codon homozygosity, grouped by degeneracy class.
        var fByDegeneracy = new Dictionary<int, List<double>>();

        foreach (var aaGroup in CodonToAminoAcid.GroupBy(kv => kv.Value))
        {
            // Exclude stop codons ('*'); they are not amino acids and not counted in Nc.
            if (aaGroup.Key == '*') continue;

            var synonymousCodons = aaGroup.Select(kv => kv.Key).ToList();
            int degeneracy = synonymousCodons.Count;

            if (degeneracy == 1) continue; // Met / Trp handled by SingleCodonAminoAcidContribution.

            int n = synonymousCodons.Sum(c => counts.GetValueOrDefault(c, 0));
            if (n <= 1) continue; // F̂ undefined for n ≤ 1 (denominator n − 1); Fuglsang 2004.

            // Wright Eq. (1): F̂ = (n·Σ p_i² − 1)/(n − 1), p_i = n_i/n.
            double sumPSquared = 0;
            foreach (var codon in synonymousCodons)
            {
                double p = (double)counts.GetValueOrDefault(codon, 0) / n;
                sumPSquared += p * p;
            }
            double f = (n * sumPSquared - 1) / (n - 1);

            if (!fByDegeneracy.TryGetValue(degeneracy, out var list))
            {
                list = new List<double>();
                fByDegeneracy[degeneracy] = list;
            }
            list.Add(f);
        }

        // Wright Eq. (4): the class average F̂ substitutes for any amino acid that cannot
        // be estimated within the same degeneracy class.
        double? f2 = AverageOrNull(fByDegeneracy, 2);
        double? f3 = AverageOrNull(fByDegeneracy, 3);
        double? f4 = AverageOrNull(fByDegeneracy, 4);
        double? f6 = AverageOrNull(fByDegeneracy, 6);

        // Wright Eq. (5a): when isoleucine (the only 3-fold amino acid) cannot be
        // estimated, F̂₃ = (F̂₂ + F̂₄)/2.
        if (f3 is null && f2 is not null && f4 is not null)
            f3 = (f2.Value + f4.Value) / 2.0;

        // Wright Eq. (3): Nc = 2 + 9/F̂₂ + 1/F̂₃ + 5/F̂₄ + 3/F̂₆.
        // A class with no estimable F̂ (and, for Ile, no Eq. 5a fallback) contributes its
        // full codon count, i.e. all its codons are assumed effectively present.
        double enc = SingleCodonAminoAcidContribution
            + ClassContribution(TwoFoldAminoAcidCount, f2)
            + ClassContribution(ThreeFoldAminoAcidCount, f3)
            + ClassContribution(FourFoldAminoAcidCount, f4)
            + ClassContribution(SixFoldAminoAcidCount, f6);

        return Math.Min(MaxEffectiveCodons, Math.Max(MinEffectiveCodons, enc));
    }

    private static double? AverageOrNull(Dictionary<int, List<double>> fByDegeneracy, int degeneracy)
        => fByDegeneracy.TryGetValue(degeneracy, out var list) && list.Count > 0
            ? list.Average()
            : null;

    private static double ClassContribution(int aminoAcidCount, double? averageF)
        // No estimable homozygosity for the class ⇒ assume all codons of every amino acid
        // in the class are effectively in use (contribution equals the codon count).
        => averageF is double f && f > 0 ? aminoAcidCount / f : aminoAcidCount;

    #endregion

    #region Codon Usage Statistics

    /// <summary>
    /// Gets comprehensive codon usage statistics.
    /// </summary>
    public static CodonUsageStatistics GetStatistics(DnaSequence sequence)
    {
        ArgumentNullException.ThrowIfNull(sequence);
        return GetStatisticsCore(sequence.Sequence);
    }

    /// <summary>
    /// Gets codon usage statistics from a raw sequence string.
    /// </summary>
    public static CodonUsageStatistics GetStatistics(string sequence)
    {
        if (string.IsNullOrEmpty(sequence))
            return new CodonUsageStatistics(
                new Dictionary<string, int>(),
                new Dictionary<string, double>(),
                0, 0, 0, 0, 0, 0);

        return GetStatisticsCore(NormalizeCodingSequence(sequence));
    }

    private static CodonUsageStatistics GetStatisticsCore(string seq)
    {
        var counts = CountCodonsCore(seq);
        var rscu = CalculateRscu(counts, GeneticCode.Standard);
        double enc = CalculateEncCore(seq);

        int totalCodons = counts.Values.Sum();

        // GC content at codon positions 1/2/3 over all valid codons (EMBOSS cusp:
        // "1st/2nd/3rd letter GC"). Reported as a percentage of valid codons.
        int gc1 = 0, gc2 = 0, gc3 = 0;
        int positionCount = 0;

        // GC3s: frequency of G/C at the THIRD position of *synonymous* codons, i.e.
        // excluding Met, Trp and termination codons (Peden 1999, CodonW thesis §1.8.2.1.3:
        // "the frequency of G or C nucleotides present at the third position of synonymous
        // codons (i.e. excluding Met, Trp and termination codons)").
        int gc3sCount = 0;          // numerator: synonymous codons with G/C at position 3
        int synonymousCodonCount = 0; // denominator: codons at synonymous third positions

        for (int i = 0; i + 3 <= seq.Length; i += 3)
        {
            string codon = seq.Substring(i, 3);
            if (IsValidCodon(codon))
            {
                gc1 += IsGC(codon[0]) ? 1 : 0;
                gc2 += IsGC(codon[1]) ? 1 : 0;
                gc3 += IsGC(codon[2]) ? 1 : 0;
                positionCount++;

                if (IsSynonymousAtThirdPosition(codon))
                {
                    synonymousCodonCount++;
                    gc3sCount += IsGC(codon[2]) ? 1 : 0;
                }
            }
        }

        double gc1Percent = positionCount > 0 ? (double)gc1 / positionCount * 100 : 0;
        double gc2Percent = positionCount > 0 ? (double)gc2 / positionCount * 100 : 0;
        double gc3Percent = positionCount > 0 ? (double)gc3 / positionCount * 100 : 0;
        // GC3s expressed as a percentage for consistency with GC1/GC2/GC3 above
        // (CodonW reports it as a fraction in [0,1]).
        double gc3s = synonymousCodonCount > 0 ? (double)gc3sCount / synonymousCodonCount * 100 : 0;

        return new CodonUsageStatistics(
            CodonCounts: counts,
            Rscu: rscu,
            Enc: enc,
            TotalCodons: totalCodons,
            Gc1: gc1Percent,
            Gc2: gc2Percent,
            Gc3: gc3Percent,
            Gc3s: gc3s);
    }

    private static bool IsGC(char c) => c is 'G' or 'C';

    // A codon is "synonymous at the third position" iff its amino acid has more than one
    // codon (degeneracy > 1). This excludes Met (ATG), Trp (TGG) and the three stop codons,
    // exactly the set CodonW omits from GC3s (Peden 1999, §1.8.2.1.3).
    private static bool IsSynonymousAtThirdPosition(string codon)
    {
        if (!CodonToAminoAcid.TryGetValue(codon, out char aa) || aa == '*')
            return false;
        return CodonToAminoAcid.Count(kv => kv.Value == aa) > 1;
    }

    #endregion

    #region Codon Table

    // Standard genetic code (NCBI table 1) in DNA spelling, derived from the canonical
    // GeneticCode (no private copy of the code).
    private static readonly IReadOnlyDictionary<string, char> CodonToAminoAcid =
        GeneticCode.Standard.CodonTable.ToDictionary(kv => kv.Key.Replace('U', 'T'), kv => kv.Value);

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
