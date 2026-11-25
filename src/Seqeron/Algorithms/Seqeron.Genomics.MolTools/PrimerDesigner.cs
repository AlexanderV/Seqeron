using System.Text;

namespace Seqeron.Genomics.MolTools;

/// <summary>
/// Designs PCR primers for DNA sequences with various quality criteria.
/// </summary>
public static partial class PrimerDesigner
{

    /// <summary>
    /// Former fixed search flank (bp). <see cref="DesignPrimers"/> now uses Primer3's search region
    /// (included region + <see cref="PrimerPairOptions.ProductSizeRanges"/>); kept for source compatibility.
    /// </summary>
    [Obsolete("DesignPrimers now uses Primer3's search region: PrimerPairOptions.IncludedRegion and ProductSizeRanges.")]
    public const int PrimerSearchFlank = 200;

    /// <summary>
    /// Library default maximum |Tm_forward − Tm_reverse| (°C) for a primer pair (Addgene "within 5°C", protocols/primer-design, re-verified 2026-10-14;
    /// stricter than Primer3's PRIMER_PAIR_MAX_DIFF_TM default of 100) — the default of
    /// <see cref="PrimerPairOptions.MaxTmDifference"/>.
    /// </summary>
    public const double MaxPairTmDifference = 5.0;

    /// <summary>Primer3's default PRIMER_PAIR_MAX_DIFF_TM (°C), <c>pr_set_default_global_args_1</c>.</summary>
    public const double Primer3MaxPairTmDifference = 100.0;

    /// <summary>
    /// Primer3 product melting temperature (<c>oligotm.c</c> <c>long_seq_tm</c>, PRIMER_PAIR_k_PRODUCT_TM):
    /// Tm = 81.5 + 16.6·log10([Mon]_eq/1000) + 41·(G+C)/N − 600/N, with [Mon]_eq = [Mon] +
    /// 120·√([Mg²⁺] − [dNTP]) (mM) and N the product length (G/C counted case-insensitively; every
    /// character counts towards N). Defaults are Primer3's primer conditions (50 mM, 1.5 mM Mg²⁺, 0.6 mM dNTP).
    /// </summary>
    /// <param name="product">Product (amplicon) sequence.</param>
    /// <param name="monovalentMillimolar">Monovalent cation concentration, mM (≥ 0).</param>
    /// <param name="divalentMillimolar">Mg²⁺ concentration, mM (≥ 0).</param>
    /// <param name="dntpMillimolar">dNTP concentration, mM (≥ 0).</param>
    /// <returns>Product Tm in °C.</returns>
    /// <exception cref="ArgumentException">Empty product.</exception>
    /// <exception cref="ArgumentOutOfRangeException">Negative concentrations or zero total cation.</exception>
    public static double CalculateProductMeltingTemperaturePrimer3(
        string product,
        double monovalentMillimolar = Primer3MonovalentMillimolar,
        double divalentMillimolar = Primer3DivalentMillimolar,
        double dntpMillimolar = Primer3DntpMillimolar)
    {
        if (string.IsNullOrEmpty(product))
            throw new ArgumentException("Product cannot be null or empty.", nameof(product));
        if (!(monovalentMillimolar >= 0) || !(divalentMillimolar >= 0) || !(dntpMillimolar >= 0))
            throw new ArgumentOutOfRangeException(nameof(monovalentMillimolar), "Concentrations must be ≥ 0 mM.");
        double monovalentEq = Primer3MonovalentEquivalent(monovalentMillimolar, divalentMillimolar, dntpMillimolar);
        if (!(monovalentEq > 0))
            throw new ArgumentOutOfRangeException(nameof(monovalentMillimolar), "Total monovalent-equivalent cation concentration must be > 0 mM.");
        // long_seq_tm: G+C (canonical counter, case-insensitive) over the full product length.
        int gc = product.AsSpan().CountGcAndValidNucleotides().GcCount;
        return LongSeqTm(gc, product.Length, monovalentEq);
    }

    // long_seq_tm (oligotm.c) with DMSO = formamide = 0; the canonical salt-adjusted GC formula.
    [System.ComponentModel.EditorBrowsable(System.ComponentModel.EditorBrowsableState.Never)] // shared with the separately licensed Seqeron.Genomics.OligoThermodynamics assembly; not part of the API
    public static double LongSeqTm(int gc, int length, double monovalentEqMillimolar) =>
        ThermoConstants.CalculateSaltAdjustedTm((double)gc / length, length, monovalentEqMillimolar / 1000.0);

    // Primer3 is_lowercase_masked: only a / c / g / t count as masked.
    [System.ComponentModel.EditorBrowsable(System.ComponentModel.EditorBrowsableState.Never)] // shared with the separately licensed Seqeron.Genomics.OligoThermodynamics assembly; not part of the API
    public static bool IsLowercaseMaskedBase(char c) => c is 'a' or 'c' or 'g' or 't';

    /// <summary>
    /// Calculates the "basic" melting temperature for DNA primers (OligoCalc basic Tm,
    /// Kibbe 2007, NAR 35:W43): Wallace rule Tm = 2(A+T) + 4(G+C) (Thein &amp; Wallace 1986) for
    /// &lt; 14 counted bases, and Tm = 64.9 + 41·(G+C − 16.4)/N (customarily attributed to
    /// Marmur &amp; Doty 1962) for ≥ 14 counted bases. Both formulas assume fixed standard conditions
    /// (50 nM primer, 50 mM Na+, pH 7.0); use <see cref="CalculateMeltingTemperatureWithSalt"/>
    /// for another [Na+], or <see cref="CalculateMeltingTemperatureNN"/> for a nearest-neighbor Tm.
    /// A, C, G, T and U are counted (case-insensitive; U read as T, as Biopython <c>MeltingTemp._check</c>
    /// back-transcribes RNA for <c>Tm_Wallace</c>/<c>Tm_GC</c>); all other characters are ignored.
    /// Delegates to the canonical <see cref="ThermoConstants.CalculateBasicTm"/>.
    /// </summary>
    public static double CalculateMeltingTemperature(string primer) => ThermoConstants.CalculateBasicTm(primer);

    /// <summary>
    /// Calculates the OligoCalc "salt adjusted" melting temperature (Kibbe 2007, NAR 35:W43),
    /// the [Na+]-aware counterpart of <see cref="CalculateMeltingTemperature(string)"/>:
    /// <list type="bullet">
    /// <item>&lt; 14 valid bases: Tm = 2(A+T) + 4(G+C) − 16.6·log10(0.050) + 16.6·log10([Na+]) —
    /// the Wallace rule (defined at 50 mM Na+) shifted by the Schildkraut–Lifson relative correction;</item>
    /// <item>≥ 14 valid bases: Tm = 100.5 + 41·(G+C)/N − 820/N + 16.6·log10([Na+]).</item>
    /// </list>
    /// [Na+] enters in mol/L (converted from the mM argument). The basic formulas already assume
    /// 50 mM Na+, so the salt term is never added on top of them. Result rounded to one decimal.
    /// A/C/G/T/U are counted as in <see cref="CalculateMeltingTemperature(string)"/> (U read as T);
    /// returns 0 for null/empty input or no counted bases.
    /// </summary>
    /// <param name="primer">Primer sequence.</param>
    /// <param name="naConcentration">Na+ concentration in mM (default: 50). Must be &gt; 0.</param>
    /// <returns>Salt-adjusted melting temperature in °C.</returns>
    /// <exception cref="ArgumentOutOfRangeException">When <paramref name="naConcentration"/> is not a positive finite number.</exception>
    public static double CalculateMeltingTemperatureWithSalt(string primer, double naConcentration = 50)
    {
        if (!(naConcentration > 0) || double.IsInfinity(naConcentration))
            throw new ArgumentOutOfRangeException(nameof(naConcentration), naConcentration,
                "Na+ concentration must be a positive finite value in mM.");

        if (string.IsNullOrEmpty(primer))
            return 0;

        var (at, gc) = ThermoConstants.CountBasicTmBases(primer);
        if (at + gc == 0)
            return 0;

        return Math.Round(
            ThermoConstants.CalculateOligoCalcSaltAdjustedTm(at, gc, naConcentration / 1000.0), 1);
    }

    // ---- Primer3 default primer Tm (PRIMER-DESIGN-001) -------------------------------------
    // Reproduces Primer3's seqtm()/oligotm() with PRIMER_TM_FORMULA = 1 (SantaLucia 1998) and
    // PRIMER_SALT_CORRECTIONS = 1 (SantaLucia 1998), the Primer3 ≥ 2.0 defaults, including the
    // divalent→monovalent equivalence of von Ahsen et al. (2001) and the long_seq_tm fallback
    // for oligos longer than MAX_NN_TM_LENGTH (36). Source: primer3 src/oligotm.c (oligotm,
    // seqtm, long_seq_tm, divalent_to_monovalent, symmetry; SantaLucia_1998_dH/dS tables);
    // SantaLucia (1998) PNAS 95:1460 Table 2 (NN + terminal-initiation terms); von Ahsen et al.
    // (2001) Clin Chem 47:1956. Cross-checked against primer3-py 2.3.1 primer3.calc_tm(max_nn_length=36)
    // (calc_tm's own default is max_nn_length=60; see the maxNearestNeighborLength overload).

    /// <summary>Primer3 MAX_NN_TM_LENGTH: longer oligos use the long_seq_tm GC/length formula.</summary>
    public const int Primer3MaxNnTmLength = 36;
    private const double Primer3DivalentFactor = 120.0;

    /// <summary>
    /// Calculates GC content as a percentage.
    /// </summary>
    public static double CalculateGcContent(string sequence) =>
        string.IsNullOrEmpty(sequence) ? 0 : sequence.CalculateGcContentFast();

    /// <summary>
    /// Finds the longest homopolymer run (consecutive identical nucleotides, case-insensitive) — the
    /// quantity Primer3 limits with <c>PRIMER_MAX_POLY_X</c> (a primer fails when the run exceeds
    /// it; <c>libprimer3.cc</c> <c>_pr_violates_poly_x</c>). As in Primer3, N is a wildcard that takes
    /// the worst case: a forward scan assigns each N to the preceding non-N base, a reverse scan to
    /// the following one, and the longer run is reported (e.g. ANA → 3, GNGNG → 5, ANGNG → 4).
    /// </summary>
    /// <returns>0 for null/empty input, otherwise ≥ 1.</returns>
    public static int FindLongestHomopolymer(string sequence)
    {
        if (string.IsNullOrEmpty(sequence))
            return 0;

        string seq = sequence.ToUpperInvariant();
        int len = seq.Length;

        // Forward scan (N counts as the last seen non-N base; leading Ns as the first non-N base).
        char lastNonN = seq[0];
        if (lastNonN == 'N')
        {
            for (int i = 1; i < len; i++)
                if (seq[i] != 'N') { lastNonN = seq[i]; break; }
        }
        int run = 1, maxRun = 1;
        bool hasN = false;
        for (int i = 1; i < len; i++)
        {
            if (seq[i] == 'N') { hasN = true; run++; }
            else if (seq[i] == lastNonN) run++;
            else { run = 1; lastNonN = seq[i]; }
            if (run > maxRun) maxRun = run;
        }
        if (!hasN)
            return maxRun;

        // Reverse scan (N counts as the next non-N base).
        lastNonN = seq[len - 1];
        if (lastNonN == 'N')
        {
            for (int i = len - 2; i >= 0; i--)
                if (seq[i] != 'N') { lastNonN = seq[i]; break; }
        }
        run = 1;
        for (int i = len - 2; i >= 0; i--)
        {
            if (seq[i] == 'N' || seq[i] == lastNonN) run++;
            else { run = 1; lastNonN = seq[i]; }
            if (run > maxRun) maxRun = run;
        }
        return maxRun;
    }

    /// <summary>
    /// Library screen (unsourced heuristic, no Primer3 counterpart): the number of consecutive copies of the
    /// most-repeated 2-mer (e.g. ATATAT = 3 units of AT), case-insensitive; inputs shorter than 4 nt return 0.
    /// Any 2-mer counts, including homo-dinucleotides (AAAA = 2 units of AA), and no base is treated as a
    /// wildcard. <see cref="PrimerParameters.MaxDinucleotideRepeats"/> limits this count in
    /// <see cref="EvaluatePrimer"/> (library default 4 in <see cref="DefaultParameters"/>). No authoritative
    /// published definition of this screen or of the limit 4 was found (audit round 3, A3-8): Primer3 has no
    /// dinucleotide-repeat setting (<c>primer3_manual.htm</c> lists only PRIMER_MAX_POLY_X, the mononucleotide run =
    /// <see cref="FindLongestHomopolymer"/>, and repeat-library mispriming PRIMER_MISPRIMING_LIBRARY), so
    /// <see cref="Primer3DefaultParameters"/> disables it (<see cref="int.MaxValue"/>). For Primer3's sourced
    /// alternatives use PRIMER_MAX_POLY_X (<see cref="PrimerParameters.MaxHomopolymer"/>) and the mispriming
    /// library PRIMER_MISPRIMING_LIBRARY (<see cref="PrimerParameters.MisprimingLibrary"/>).
    /// </summary>
    /// <param name="sequence">Nucleotide sequence.</param>
    /// <returns>Unit count of the longest dinucleotide tandem repeat (≥ 1 for inputs of ≥ 4 nt), 0 below 4 nt.</returns>
    public static int FindLongestDinucleotideRepeat(string sequence)
    {
        if (string.IsNullOrEmpty(sequence) || sequence.Length < 4)
            return 0;

        var seq = sequence.ToUpperInvariant();
        int maxRepeats = 0;

        for (int i = 0; i < seq.Length - 3; i++)
        {
            string dinuc = seq.Substring(i, 2);
            int repeats = 1;
            int j = i + 2;

            while (j + 1 < seq.Length && seq.Substring(j, 2) == dinuc)
            {
                repeats++;
                j += 2;
            }

            maxRepeats = Math.Max(maxRepeats, repeats);
        }

        return maxRepeats;
    }

    /// <summary>
    /// Library screen (sequence-only heuristic, no Primer3 counterpart): <c>true</c> when the sequence contains
    /// two non-overlapping segments of <paramref name="minStemLength"/> bases that are exact Watson–Crick reverse
    /// complements of each other (an antiparallel stem) separated by at least
    /// <paramref name="minLoopLength"/> unpaired bases. Case-insensitive; only A·T and G·C pairs (a stem window with
    /// N or another non-ACGT symbol never pairs, at any length); no G·T wobble, mismatches or energies.
    /// The default minimum loop of 3 nt is Primer3's hairpin minimum (<c>thal.c</c> <c>min_hrpn_loop = 3</c>); the
    /// default minimum stem of 4 bp is an unsourced library threshold (no published definition of a
    /// "≥ 4-bp stem + ≥ 3-nt loop" rule was found, audit round 3, A3-8). It is not equivalent to Primer3's
    /// hairpin screen: e.g. <c>AAAACCCTTTT</c> is flagged although Primer3 <c>calc_hairpin</c> finds no structure,
    /// and <c>CAGTAAAACCCTTTTGCAGC</c> is flagged although its hairpin Tm is 37.65 °C (primer3-py 2.3.1,
    /// 50 mM / 1.5 mM / 0.6 mM / 50 nM), below PRIMER_MAX_HAIRPIN_TH = 47 °C. The sourced Primer3 hairpin screen
    /// (PRIMER_HAIRPIN_TH) is <see cref="CalculatePrimer3OligoStructure"/> /
    /// <see cref="CalculateHairpinThermodynamicsNtthal(string, double)"/>, which
    /// <see cref="EvaluatePrimer"/> uses by default; primer evaluation uses this screen only with
    /// <see cref="PrimerStructureScreen.Heuristic"/>. It is also the fallback hairpin screen of the probe designers
    /// (<see cref="ProbeDesigner.ProbeParameters"/>: probes longer than the THAL_MAX_ALIGN of the thermodynamic screen,
    /// non-ACGT probes, <see cref="ProbeDesigner.ProbeStructureScreen.Heuristic"/>; audit round 5, A5-2). Uses an O(n²) scan below 100 nt and a suffix tree
    /// at ≥ 100 nt (identical results).
    /// </summary>
    /// <param name="sequence">DNA sequence to check.</param>
    /// <param name="minStemLength">Minimum stem length (default 4; ≥ 1 — a stem has at least one base pair).</param>
    /// <param name="minLoopLength">Minimum loop length (default 3; ≥ 0 — 0 allows the two arms to abut, below
    /// Primer3's 3-nt hairpin minimum).</param>
    /// <returns>True if hairpin potential detected; <c>false</c> for a null or empty sequence.</returns>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="minStemLength"/> &lt; 1 or
    /// <paramref name="minLoopLength"/> &lt; 0.</exception>
    public static bool HasHairpinPotential(string sequence, int minStemLength = 4, int minLoopLength = 3)
    {
        // Audit round 7, A7-3: a 0-bp stem is no stem (it returned true for any long-enough sequence) and a negative
        // stem / loop crashed in Substring.
        ArgumentOutOfRangeException.ThrowIfLessThan(minStemLength, 1);
        ArgumentOutOfRangeException.ThrowIfNegative(minLoopLength);
        if (string.IsNullOrEmpty(sequence) || sequence.Length < 2L * minStemLength + minLoopLength)
            return false;

        var seq = sequence.ToUpperInvariant();

        // For short sequences (typical primers), use simple O(n²) approach
        // Break-even point is ~100bp based on suffix tree construction overhead
        if (seq.Length < 100)
        {
            return HasHairpinPotentialSimple(seq, minStemLength, minLoopLength);
        }

        // For longer sequences, use suffix tree for O(n) lookup
        return HasHairpinPotentialWithSuffixTree(seq, minStemLength, minLoopLength);
    }

    /// <summary>
    /// Simple O(n²) hairpin detection for short sequences.
    /// </summary>
    private static bool HasHairpinPotentialSimple(string seq, int minStemLength, int minLoopLength)
    {
        // Check for self-complementary regions
        for (int i = 0; i <= seq.Length - minStemLength; i++)
        {
            string fragment = seq.Substring(i, minStemLength);
            // Look for complementary sequence at least minLoopLength positions away
            for (int j = i + minStemLength + minLoopLength; j <= seq.Length - minStemLength; j++)
            {
                string target = seq.Substring(j, minStemLength);
                if (AreComplementary(fragment, OligoSymmetry.Reverse(target)))
                    return true;
            }
        }

        return false;
    }

    private static readonly System.Buffers.SearchValues<char> AcgtBases = System.Buffers.SearchValues.Create("ACGT");

    /// <summary>
    /// Suffix tree-based O(n) hairpin detection for long sequences.
    /// 
    /// Algorithm: A hairpin forms when a substring S at position i is complementary
    /// to a substring at position j (in reverse). This is equivalent to:
    /// - seq[i..i+k] being present in revComp at some position p
    /// - The positions must satisfy: j = n - p - k, and j >= i + k + minLoopLength
    /// 
    /// We build a suffix tree on seq and search for all substrings of revComp,
    /// checking if any match satisfies the loop constraint.
    /// </summary>
    private static bool HasHairpinPotentialWithSuffixTree(string seq, int minStemLength, int minLoopLength)
    {
        var revComp = DnaSequence.GetReverseComplementString(seq);
        var tree = global::SuffixTree.SuffixTree.Build(seq);

        // For each position in revComp, find matches in seq via suffix tree
        // and check if they form valid hairpin (sufficient loop distance)
        int n = seq.Length;

        // Slide through revComp looking for stems
        for (int p = 0; p <= n - minStemLength; p++)
        {
            var pattern = revComp.AsSpan(p, minStemLength);
            // Only Watson–Crick A·T / G·C pairs form a stem (IsComplementary, as in the < 100-nt scan): the IUPAC
            // reverse complement maps N→N, S→S, W→W, so a stem window with a non-ACGT base would otherwise match itself;
            // and it maps U→A, so the source window seq[n-p-k..n-p) must itself be ACGT (audit round 6, A6-1).
            if (pattern.ContainsAnyExcept(AcgtBases)
                || seq.AsSpan(n - p - minStemLength, minStemLength).ContainsAnyExcept(AcgtBases))
                continue;
            var matches = tree.FindAllOccurrences(pattern);

            foreach (int i in matches)
            {
                // Position in revComp p corresponds to position (n - p - minStemLength) in seq
                // when we reverse complement back
                int j = n - p - minStemLength;

                // Check if positions form valid hairpin: j >= i + minStemLength + minLoopLength
                // Also check i and j don't overlap with the stem itself
                if (j >= i + minStemLength + minLoopLength && j + minStemLength <= n)
                {
                    return true;
                }
                // Also check the reverse case: i is the 3' stem, j is the 5' stem
                if (i >= j + minStemLength + minLoopLength && i + minStemLength <= n)
                {
                    return true;
                }
            }
        }

        return false;
    }

    /// <summary>
    /// Checks whether two primers can form a 3′-end primer-dimer using Primer3's alignment-based
    /// (non-thermodynamic, <c>PRIMER_THERMODYNAMIC_OLIGO_ALIGNMENT=0</c>) pair 3′-complementarity:
    /// returns <c>true</c> when <see cref="CalculatePrimerDimerEndComplementarity"/> ≥
    /// <paramref name="minComplementarity"/>. For ACGT primers the score is an integer, so the
    /// default threshold 4 is exactly Primer3's default <c>PRIMER_PAIR_MAX_COMPL_END = 3.00</c>
    /// (a pair fails when <c>compl_end &gt; 3</c>). The thermodynamic counterpart (Primer3's default
    /// mode) is <see cref="CalculatePrimer3PairComplementarity"/>.
    /// </summary>
    /// <param name="primer1">First primer (5′→3′).</param>
    /// <param name="primer2">Second primer (5′→3′).</param>
    /// <param name="minComplementarity">Minimum 3′-anchored complementarity score that flags a dimer (default 4).</param>
    /// <returns><c>false</c> for null/empty primers.</returns>
    public static bool HasPrimerDimer(string primer1, string primer2, int minComplementarity = 4)
    {
        if (string.IsNullOrEmpty(primer1) || string.IsNullOrEmpty(primer2))
            return false;
        return CalculatePrimerDimerEndComplementarity(primer1, primer2) >= minComplementarity;
    }

    /// <summary>
    /// Primer3 alignment-mode pair 3′-complementarity <c>compl_end</c> (Primer3 <c>libprimer3.cc</c>
    /// <c>characterize_pair</c>, <c>PRIMER_PAIR_COMPL_END</c> with
    /// <c>PRIMER_THERMODYNAMIC_OLIGO_ALIGNMENT=0</c>; Rozen &amp; Skaletsky 2000; Untergasser et al. 2012):
    /// the maximum of <c>align(p1, revcomp(p2))</c> and <c>align(p2, revcomp(p1))</c>, where
    /// <c>align</c> is the <c>dpal</c> end-anchored (<c>DPAL_GLOBAL_END</c>) alignment — it must end
    /// at the 3′-terminal base of the first sequence — scored +1 per complementary base pair, −1 per
    /// mismatch, −0.25 against N (any non-ACGT character is scored as N), −2 per single-base gap
    /// (max gap 1), and floored at 0. Two primers whose 3′-terminal k bases are reverse complements
    /// score k; a homopolymer against itself (e.g. A₈/A₈) scores 0 because it cannot pair.
    /// </summary>
    /// <param name="primer1">First primer (5′→3′), case-insensitive.</param>
    /// <param name="primer2">Second primer (5′→3′), case-insensitive.</param>
    /// <returns>The Primer3 <c>compl_end</c> score (≥ 0); 0 for null/empty input.</returns>
    public static double CalculatePrimerDimerEndComplementarity(string primer1, string primer2)
    {
        if (string.IsNullOrEmpty(primer1) || string.IsNullOrEmpty(primer2))
            return 0;
        string p1 = primer1.ToUpperInvariant();
        string p2 = primer2.ToUpperInvariant();
        // Primer3 compares s1 with s2 taken from the same (top) strand, i.e. the right primer's
        // template-strand copy = revcomp(right primer); then also s2_rev against s1_rev.
        double a = DpalGlobalEndScore(p1, ReverseComplementPrimer3(p2));
        double b = DpalGlobalEndScore(p2, ReverseComplementPrimer3(p1));
        return Math.Max(a, b);
    }

    /// <summary>
    /// Primer3 alignment-mode self 3′-complementarity <c>self_end</c>
    /// (<c>oligo_compl</c>: <c>align(oligo, revcomp(oligo), DPAL_GLOBAL_END)</c>, Primer3
    /// <c>PRIMER_LEFT/RIGHT_SELF_END</c>; default limit <c>PRIMER_MAX_SELF_END = 3.00</c>).
    /// </summary>
    /// <param name="primer">Primer (5′→3′), case-insensitive.</param>
    /// <returns>The Primer3 <c>self_end</c> score (≥ 0); 0 for null/empty input.</returns>
    public static double CalculatePrimerSelfEndComplementarity(string primer)
    {
        if (string.IsNullOrEmpty(primer))
            return 0;
        string p = primer.ToUpperInvariant();
        return DpalGlobalEndScore(p, ReverseComplementPrimer3(p));
    }

    /// <summary>Primer3 default PRIMER_MAX_SELF_ANY (alignment mode, <c>pr_set_default_global_args_1</c>).</summary>
    public const double Primer3MaxSelfAny = 8.0;

    /// <summary>Primer3 default PRIMER_MAX_SELF_END (alignment mode).</summary>
    public const double Primer3MaxSelfEnd = 3.0;

    /// <summary>Primer3 default PRIMER_PAIR_MAX_COMPL_ANY (alignment mode).</summary>
    public const double Primer3MaxPairComplAny = 8.0;

    /// <summary>Primer3 default PRIMER_PAIR_MAX_COMPL_END (alignment mode).</summary>
    public const double Primer3MaxPairComplEnd = 3.0;

    /// <summary>Primer3 default PRIMER_INTERNAL_MAX_SELF_ANY and PRIMER_INTERNAL_MAX_SELF_END (alignment mode, 12.00).</summary>
    public const double Primer3InternalMaxSelfComplementarity = 12.0;

    /// <summary>Reverse complement as Primer3's <c>p3_reverse_complement</c> writes it: ACGT complemented, every other character becomes N.</summary>
    /// <exception cref="ArgumentNullException"><paramref name="seq"/> is null.</exception>
    [System.ComponentModel.EditorBrowsable(System.ComponentModel.EditorBrowsableState.Never)] // shared with the separately licensed Seqeron.Genomics.OligoThermodynamics assembly; not part of the API
    public static string ReverseComplementPrimer3(string seq)
    {
        ArgumentNullException.ThrowIfNull(seq);
        var chars = new char[seq.Length];
        for (int i = 0; i < seq.Length; i++)
        {
            chars[seq.Length - 1 - i] = seq[i] switch
            {
                'A' => 'T',
                'T' => 'A',
                'G' => 'C',
                'C' => 'G',
                _ => 'N',
            };
        }
        return new string(chars);
    }

    // dpal default primer-picking scoring (dpal.c set_dpal_args): identity matrix ×100.
    public const int DpalMatch = 100, DpalMismatch = -100, DpalN = -25, DpalGap = -200;

    [System.ComponentModel.EditorBrowsable(System.ComponentModel.EditorBrowsableState.Never)] // shared with the separately licensed Seqeron.Genomics.OligoThermodynamics assembly; not part of the API
    public static int DpalSsm(char x, char y)
    {
        bool xn = x is not ('A' or 'C' or 'G' or 'T');
        bool yn = y is not ('A' or 'C' or 'G' or 'T');
        if (xn || yn) return DpalN;
        return x == y ? DpalMatch : DpalMismatch;
    }

    /// <summary>
    /// Primer3 <c>align(X, Y, DPAL_GLOBAL_END)</c> score / 100 floored at 0 (libprimer3.cc
    /// <c>align</c>), computed by a line-by-line port of dpal.c
    /// <c>_dpal_long_nopath_maxgap1_global_end</c> (the routine Primer3 runs in DPM_FAST mode with
    /// max_gap = 1); inputs too short for that routine (|X| ≤ 3 or |Y| = 1, where the C code reads
    /// past the sequence end) use dpal.c's <c>_dpal_generic</c> recurrence for GLOBAL_END.
    /// </summary>
    /// <exception cref="ArgumentNullException">An argument is null.</exception>
    [System.ComponentModel.EditorBrowsable(System.ComponentModel.EditorBrowsableState.Never)] // shared with the separately licensed Seqeron.Genomics.OligoThermodynamics assembly; not part of the API
    public static double DpalGlobalEndScore(string x, string y)
    {
        ArgumentNullException.ThrowIfNull(x);
        ArgumentNullException.ThrowIfNull(y);
        if (x.Length == 0 || y.Length == 0)
            return 0.0; // nothing to align
        int xlen = x.Length, ylen = y.Length;
        // The fast routine reads Y[ylen] (C's NUL terminator) when |X| ≤ 3 and Y[1] when |Y| = 1.
        int smax = xlen < 4 || ylen < 2 ? DpalGlobalEndGeneric(x, y) : DpalGlobalEndFast(x, y);
        return smax < 0 ? 0.0 : smax / 100.0;
    }

    private static int DpalGlobalEndFast(string x, string y)
    {
        int xlen = x.Length, ylen = y.Length;
        const int gap = DpalGap;
        var s0 = new int[xlen];
        var s1 = new int[xlen];
        var s2 = new int[xlen];
        int score, a;

        int smax = DpalSsm(x[xlen - 1], y[0]);
        for (int j = 0; j < xlen; j++) s0[j] = DpalSsm(x[j], y[0]);

        s1[0] = DpalSsm(x[0], y[1]);
        for (int j = 1; j < xlen; j++)
        {
            score = s0[j - 1];
            if (j > 1 && (a = s0[j - 2] + gap) > score) score = a;
            score += DpalSsm(x[j], y[1]);
            if (score > smax && j == xlen - 1) smax = score;
            s1[j] = score;
        }

        int k = ylen - xlen / 2 + 1;
        if (k < 1) k = 1;

        // Rectangular part.
        for (int j = 2; j < k + 1; j++)
        {
            s2[0] = DpalSsm(x[0], y[j]);
            score = s1[0];
            if ((a = s0[0] + gap) > score) score = a;
            score += DpalSsm(x[1], y[j]);
            s2[1] = score;
            for (int i = 2; i < xlen - 1; i++)
            {
                score = s1[i - 2];
                if ((a = s0[i - 1]) > score) score = a;
                score += gap;
                if ((a = s1[i - 1]) > score) score = a;
                score += DpalSsm(x[i], y[j]);
                s2[i] = score;
            }
            score = s1[xlen - 3];
            if ((a = s0[xlen - 2]) > score) score = a;
            score += gap;
            if ((a = s1[xlen - 2]) > score) score = a;
            score += DpalSsm(x[xlen - 1], y[j]);
            s2[xlen - 1] = score;
            if (score > smax) smax = score;
            (s0, s1, s2) = (s1, s2, s0);
        }

        // Triangular part (cells left of the band are not recomputed, exactly as dpal.c).
        int t = 2;
        for (int j = k + 1; j < ylen; j++)
        {
            for (int i = t; i < xlen - 1; i++)
            {
                score = s1[i - 2];
                if ((a = s0[i - 1]) > score) score = a;
                score += gap;
                if ((a = s1[i - 1]) > score) score = a;
                score += DpalSsm(x[i], y[j]);
                s2[i] = score;
            }
            t += 2;
            score = s1[xlen - 3];
            if ((a = s0[xlen - 2]) > score) score = a;
            score += gap;
            if ((a = s1[xlen - 2]) > score) score = a;
            score += DpalSsm(x[xlen - 1], y[j]);
            s2[xlen - 1] = score;
            if (score > smax) smax = score;
            (s0, s1, s2) = (s1, s2, s0);
        }
        return smax;
    }

    // dpal.c _dpal_generic with flag DPAL_GLOBAL_END and max_gap = 1 (score only).
    private static int DpalGlobalEndGeneric(string x, string y)
    {
        int xlen = x.Length, ylen = y.Length;
        var sm = new int[xlen, ylen];
        for (int i = 0; i < xlen; i++) sm[i, 0] = DpalSsm(x[i], y[0]);
        int smax = sm[xlen - 1, 0];
        for (int j = 0; j < ylen; j++) sm[0, j] = DpalSsm(x[0], y[j]);
        for (int i = 1; i < xlen; i++)
        {
            for (int j = 1; j < ylen; j++)
            {
                long a = sm[i - 1, j - 1];
                long b = i > 1 ? (long)sm[i - 2, j - 1] + DpalGap : long.MinValue;
                long c = j > 1 ? (long)sm[i - 1, j - 2] + DpalGap : long.MinValue;
                long best;
                if (a >= b && a >= c) best = a;
                else if (b > a && b >= c) best = b;
                else best = c;
                int score = (int)(best + DpalSsm(x[i], y[j]));
                if (score >= smax && i == xlen - 1) smax = score;
                sm[i, j] = score;
            }
        }
        return smax;
    }

    /// <summary>
    /// 3′-end stability: the duplex ΔG°37 (kcal/mol, 1 M NaCl) of the last five bases (the whole
    /// primer when it is shorter than five), computed exactly as Primer3's
    /// <c>end_oligodg(seq, 5, santalucia)</c> (<c>oligotm.c</c> <c>oligodg</c>; the value Primer3 reports
    /// as <c>PRIMER_{LEFT,RIGHT}_n_END_STABILITY</c> and limits with <c>PRIMER_MAX_END_STABILITY</c>):
    /// SantaLucia (1998) Table 1 unified NN ΔG°37 values summed over the steps, plus initiation
    /// +1.96, +0.05 per terminal A·T base pair and +0.43 for a self-complementary (even-length)
    /// sequence. For a 5-mer this is identical to SantaLucia (1998)'s "initiation with terminal G·C
    /// +0.98 / terminal A·T +1.03" form (0.98+0.98 = 1.96, 1.03−0.98 = 0.05). N is accepted with
    /// Primer3's N-row values (NA 0.58, NC 1.30, NG 1.28, NT 0.88, NN 0.58 negated, etc.) and no A·T
    /// penalty.
    /// <para><b>Sign convention:</b> this method returns the physical ΔG (negative = stable); Primer3's
    /// END_STABILITY is the same quantity with the opposite sign (GCGCG → −6.86 here, 6.86 in
    /// Primer3; TATAT → −0.86 / 0.86).</para>
    /// </summary>
    /// <param name="sequence">Primer sequence (5′→3′), case-insensitive.</param>
    /// <returns>ΔG°37 in kcal/mol; 0 for null/empty input; <c>double.NaN</c> when the 3′ window
    /// contains a character other than A, C, G, T, N (Primer3 <c>OLIGOTM_ERROR</c>).</returns>
    public static double Calculate3PrimeStability(string sequence)
    {
        if (string.IsNullOrEmpty(sequence))
            return 0;

        string window = (sequence.Length > 5 ? sequence[^5..] : sequence).ToUpperInvariant();
        var idx = new int[window.Length];
        for (int i = 0; i < window.Length; i++)
        {
            int b = window[i] switch { 'A' => 0, 'C' => 1, 'G' => 2, 'T' => 3, 'N' => 4, _ => -1 };
            if (b < 0) return double.NaN;
            idx[i] = b;
        }

        // oligodg (santalucia): dg = −1960 [−430 if symmetric] [−50 per terminal A/T] + Σ table, in cal/mol
        // of −ΔG; the method returns ΔG = −dg/1000.
        int dg = -1960;
        if (OligoSymmetry.IsSymmetric(window)) dg += -430; // oligotm.c symmetry() = thal.c symmetry_thermo on upper case
        if (window[0] is 'A' or 'T') dg += -50;
        for (int i = 0; i + 1 < idx.Length; i++)
            dg += Primer3SantaLucia1998Dg[idx[i], idx[i + 1]];
        if (window[^1] is 'A' or 'T') dg += -50;

        return -dg / 1000.0;
    }

    // oligotm.c SantaLucia_1998_dG (−ΔG°37, cal/mol), rows/cols A, C, G, T, N.
    private static readonly int[,] Primer3SantaLucia1998Dg =
    {
        { 1000, 1440, 1280,  880,  880 },
        { 1450, 1840, 2170, 1280, 1450 },
        { 1300, 2240, 1840, 1440, 1300 },
        {  580, 1300, 1450, 1000,  580 },
        {  580, 1300, 1280,  880,  580 },
    };

    // ---- Nearest-neighbour salt-corrected Tm (PRIMER-NNTM-001, opt-in) ----------
    // SantaLucia & Hicks (2004) Watson-Crick NN ΔH°/ΔS° (1 M NaCl; Biopython DNA_NN4 — NOT the
    // SantaLucia 1998 / Allawi & SantaLucia 1997 set: AA/TT is −7.6/−21.3 here, −7.9/−22.2 in 1998,
    // and 2004 uses one duplex-initiation term plus a terminal A·T penalty instead of per-end
    // initiation) with the bimolecular Tm equation and published salt corrections. This is an OPT-IN
    // design Tm: the default CalculateMeltingTemperature (Wallace / Marmur-Doty) is unchanged.
    // All arithmetic is delegated to the canonical NN core ThermoConstants.CalculateNearestNeighborDuplex
    // (a line-by-line port of Biopython Tm_NN); this class only maps its API onto it:
    //   parameter set = SantaLuciaHicks2004 (DNA_NN4); R = 1.9872 (SantaLucia & Hicks 2004 Eq. 3;
    //   Biopython uses 1.987); C_T total → dnac1 = dnac2 = C_T/2 (k = C_T/4), self-complementary →
    //   dnac1 = C_T (k = C_T); self-complementarity detected from the sequence; salt mode → Biopython
    //   salt_correction method 0 / 5 / 6 / 7 with [Na⁺] (and, for method 7 only, [Mg²⁺], [dNTPs]).
    // Sources: SantaLucia J, Hicks D (2004) Annu Rev Biophys Biomol Struct 33:415, Table 1 + Eq. 3/5;
    //   Owczarzy R et al. (2004) Biochemistry 43:3537 (monovalent); Owczarzy R et al. (2008)
    //   Biochemistry 47:5336 (Mg²⁺/dNTP); Biopython 1.88 Bio.SeqUtils.MeltingTemp (reference).

    /// <summary>R of SantaLucia &amp; Hicks (2004) Eq. 3, used by the NN Tm methods of this class.</summary>
    private const double GasConstant = ThermoConstants.NnGasConstantSantaLuciaHicks2004;

    // Strand-concentration divisor x in Tm = ΔH°/(ΔS° + R·ln(C_T/x)) (SantaLucia & Hicks 2004 Eq. 3):
    // x = 4 for non-self-complementary, x = 1 for self-complementary duplexes (used by the dimer Tm).
    private const double NonSelfComplementaryFactor = 4.0;
    private const double SelfComplementaryFactor = 1.0;

    // Kelvin-to-Celsius offset.
    public const double KelvinOffset = 273.15;

    // Default total strand concentration C_T = 0.5 µM (a common PCR primer working
    // concentration). Exposed as a parameter; the caller may override.
    private const double DefaultStrandConcentrationMolar = 0.5e-6;

    // Default total strand concentration for the intermolecular dimer Tm. Primer3 / ntthal
    // uses dna_conc = 50 nM (thal.c default a->dna_conc = 50, lines 829/844); the opt-in dimer
    // methods adopt the same convention so they reproduce the ntthal reference out of the box.
    public const double DefaultDimerStrandConcentrationMolar = 50e-9;

    // The NN parameter set of this class's NN Tm, hairpin and dimer helpers (SantaLucia & Hicks 2004).
    private const NnParameterSet DesignNnParameterSet = NnParameterSet.SantaLuciaHicks2004;

    // SantaLucia & Hicks (2004) Table 1 initiation-type terms of DesignNnParameterSet
    // (duplex initiation +0.2/−5.7, terminal A·T +2.2/+6.9, symmetry 0/−1.4).
    private static readonly NnInitiationTerms DesignNnInitiation =
        ThermoConstants.GetNearestNeighborInitiation(DesignNnParameterSet);

    // Watson–Crick stack of the top-strand dinucleotide (upper-case ACGT) from DesignNnParameterSet.
    private static bool TryGetDesignStack(string dinucleotide, out (double DeltaH, double DeltaS) p) =>
        ThermoConstants.TryGetNearestNeighborStack(DesignNnParameterSet, dinucleotide, out p);

    /// <summary>Salt-correction mode for <see cref="CalculateMeltingTemperatureNN"/>.</summary>
    public enum SaltCorrectionMode
    {
        /// <summary>No correction — Tm at the SantaLucia 1 M NaCl reference state (Biopython saltcorr 0).</summary>
        None,

        /// <summary>
        /// SantaLucia (1998) / SantaLucia &amp; Hicks (2004) Eq. 5 entropy correction (Biopython saltcorr 5):
        /// ΔS°[Na] = ΔS°[1 M] + 0.368·(N − 1)·ln[Na⁺], N = oligo length (N − 1 = half the duplex phosphates);
        /// applied to ΔS° before the Tm equation.
        /// </summary>
        SantaLuciaEntropy,

        /// <summary>
        /// Owczarzy et al. (2004) monovalent 1/Tm correction (Biochemistry 43:3537; Biopython saltcorr 6):
        /// 1/Tm = 1/Tm(1 M) + (4.29·f(GC) − 3.95)·10⁻⁵·ln[Na⁺] + 9.40·10⁻⁶·ln²[Na⁺]. [Mg²⁺] and dNTPs are ignored.
        /// </summary>
        Owczarzy2004Monovalent,

        /// <summary>
        /// Owczarzy et al. (2008) divalent Mg²⁺ (and dNTP-adjusted) correction (Biochemistry 47:5336; Biopython
        /// saltcorr 7); reduces to the 2004 monovalent form when √[Mg²⁺]/[Na⁺] &lt; 0.22.
        /// </summary>
        Owczarzy2008Divalent
    }

    // SantaLucia (1998) / SantaLucia & Hicks (2004) Eq. 5 entropy salt-correction coefficient (0.368),
    // used by the Primer3 seqtm port and the dimer Tm.
    public const double SantaLuciaEntropySaltCoefficient = 0.368;

    /// <summary>
    /// Computes the duplex ΔH° (kcal/mol) and ΔS° (cal/(K·mol), 1 M NaCl) of a DNA oligonucleotide and its
    /// perfect complement with the SantaLucia &amp; Hicks (2004) nearest-neighbour parameters (Biopython
    /// <c>DNA_NN4</c>): duplex initiation, a terminal A·T penalty per A·T-closed end, the stacks and (for a
    /// self-complementary sequence) the symmetry correction. Equals Biopython
    /// <c>Tm_NN(seq, nn_table=DNA_NN4, selfcomp=…)</c>'s ΔH°/ΔS°; computed by
    /// <see cref="ThermoConstants.CalculateNearestNeighborThermodynamics"/>.
    /// </summary>
    /// <param name="sequence">DNA sequence (one strand, 5'→3'; case-insensitive).</param>
    /// <returns>(ΔH°, ΔS°, IsSelfComplementary) or <c>null</c> if the sequence is empty,
    /// shorter than 2 bases, or contains a non-ACGT character.</returns>
    public static (double DeltaH, double DeltaS, bool IsSelfComplementary)? CalculateNearestNeighborThermodynamics(string sequence)
    {
        if (string.IsNullOrEmpty(sequence) || sequence.Length < 2)
            return null;

        string seq = sequence.ToUpperInvariant();
        if (!IsAcgtOnly(seq))
            return null; // non-ACGT base present

        bool selfComp = OligoSymmetry.IsSymmetric(seq);
        var (dH, dS) = ThermoConstants.CalculateNearestNeighborThermodynamics(
            seq, parameterSet: DesignNnParameterSet, selfComplementary: selfComp, check: false);
        return (dH, dS, selfComp);
    }

    /// <summary>
    /// Computes the design melting temperature (°C) of a primer/oligonucleotide using the
    /// SantaLucia &amp; Hicks (2004) nearest-neighbour thermodynamics (Biopython <c>DNA_NN4</c>) and the
    /// bimolecular Tm equation, with an optional published salt correction. <b>Opt-in</b>: the default
    /// <see cref="CalculateMeltingTemperature(string)"/> (Wallace / Marmur-Doty) is unchanged.
    /// <para>
    /// Tm = ΔH°·1000 / (ΔS° + R·ln(C_T / x)) − 273.15, with R = 1.9872 cal/(K·mol),
    /// x = 4 for a non-self-complementary duplex and x = 1 for a self-complementary one
    /// (SantaLucia &amp; Hicks 2004, Eq. 3); salt corrections per <paramref name="saltMode"/>.
    /// Identical to Biopython <c>Tm_NN(seq, nn_table=DNA_NN4, dnac1=dnac2=C_T/2 (self-complementary:
    /// dnac1=C_T), selfcomp, Na, Mg, dNTPs, saltcorr=0/5/6/7)</c> with Biopython's R = 1.987 replaced by
    /// 1.9872 (≈ +0.005 °C); computed by the canonical <see cref="ThermoConstants.CalculateNearestNeighborDuplex"/>.
    /// </para>
    /// </summary>
    /// <param name="primer">DNA primer sequence (5'→3'). Must be ≥ 2 ACGT bases.</param>
    /// <param name="strandConcentrationMolar">Total strand concentration C_T in mol/L
    /// (default 0.5 µM).</param>
    /// <param name="sodiumMolar">Monovalent cation ([Na⁺]+[K⁺]+[Tris]/2) concentration in
    /// mol/L (default 0.05 M = 50 mM).</param>
    /// <param name="magnesiumMolar">[Mg²⁺] in mol/L (default 0; only used by the
    /// <see cref="SaltCorrectionMode.Owczarzy2008Divalent"/> mode).</param>
    /// <param name="dntpMolar">Total dNTP concentration in mol/L (default 0; sequesters Mg²⁺
    /// in the divalent mode).</param>
    /// <param name="saltMode">Which salt correction to apply (default
    /// <see cref="SaltCorrectionMode.Owczarzy2004Monovalent"/>).</param>
    /// <returns>The nearest-neighbour Tm in °C, or <c>double.NaN</c> if the sequence is
    /// empty, shorter than 2 bases, or contains a non-ACGT character.</returns>
    public static double CalculateMeltingTemperatureNN(
        string primer,
        double strandConcentrationMolar = DefaultStrandConcentrationMolar,
        double sodiumMolar = ThermoConstants.DefaultNaConcentration,
        double magnesiumMolar = 0.0,
        double dntpMolar = 0.0,
        SaltCorrectionMode saltMode = SaltCorrectionMode.Owczarzy2004Monovalent)
    {
        ValidateNnConditions(strandConcentrationMolar, sodiumMolar, magnesiumMolar, dntpMolar);

        var thermo = CalculateNearestNeighborThermodynamics(primer);
        if (thermo is null)
            return double.NaN;

        var (dH, dS, selfComp) = thermo.Value;
        return NnTm(dH, dS, primer.ToUpperInvariant(), selfComp,
            strandConcentrationMolar, sodiumMolar, magnesiumMolar, dntpMolar, saltMode);
    }

    // Parameter-domain guards: the Tm equation takes R·ln(C_T/x) and the salt corrections take
    // ln([Na⁺]) / ln([Mg²⁺]); a zero or negative concentration would make ln undefined. Reject these
    // out-of-domain inputs explicitly (Biopython raises ValueError for them too).
    private static void ValidateNnConditions(
        double strandConcentrationMolar, double sodiumMolar, double magnesiumMolar, double dntpMolar)
    {
        if (!(strandConcentrationMolar > 0))
            throw new ArgumentOutOfRangeException(nameof(strandConcentrationMolar),
                strandConcentrationMolar, "Strand concentration C_T must be > 0 mol/L.");
        if (!(sodiumMolar > 0))
            throw new ArgumentOutOfRangeException(nameof(sodiumMolar),
                sodiumMolar, "Monovalent cation concentration [Na⁺] must be > 0 mol/L.");
        if (!(magnesiumMolar >= 0))
            throw new ArgumentOutOfRangeException(nameof(magnesiumMolar),
                magnesiumMolar, "[Mg²⁺] must be ≥ 0 mol/L.");
        if (!(dntpMolar >= 0))
            throw new ArgumentOutOfRangeException(nameof(dntpMolar),
                dntpMolar, "Total dNTP concentration must be ≥ 0 mol/L.");
    }

    // Maps this class's NN Tm API onto the canonical Biopython Tm_NN core (see the region comment):
    // C_T → dnac1 = dnac2 = C_T/2 (k = C_T/4) or dnac1 = C_T for a self-complementary duplex (k = C_T);
    // mol/L → mM; salt mode → method 0/5/6/7 ([Mg²⁺]/[dNTPs] only for method 7); R = 1.9872.
    private static double NnTm(
        double deltaH, double deltaS, string saltSequence, bool selfComp,
        double strandConcentrationMolar, double sodiumMolar, double magnesiumMolar, double dntpMolar,
        SaltCorrectionMode saltMode, double gasConstant = GasConstant)
    {
        double dnac = selfComp ? strandConcentrationMolar * 1e9 : strandConcentrationMolar * 1e9 / 2.0;
        var method = saltMode switch
        {
            SaltCorrectionMode.None => NnSaltCorrection.None,
            SaltCorrectionMode.SantaLuciaEntropy => NnSaltCorrection.SantaLucia1998Entropy,
            SaltCorrectionMode.Owczarzy2004Monovalent => NnSaltCorrection.Owczarzy2004,
            SaltCorrectionMode.Owczarzy2008Divalent => NnSaltCorrection.Owczarzy2008,
            _ => throw new ArgumentOutOfRangeException(nameof(saltMode), saltMode, "Unknown salt-correction mode.")
        };
        bool divalent = method == NnSaltCorrection.Owczarzy2008;
        try
        {
            return ThermoConstants.CalculateNearestNeighborTmFromThermodynamics(
                deltaH, deltaS, saltSequence, dnac1: dnac, dnac2: dnac, selfComplementary: selfComp,
                sodium: sodiumMolar * 1000.0,
                magnesium: divalent ? magnesiumMolar * 1000.0 : 0,
                dntps: divalent ? dntpMolar * 1000.0 : 0,
                saltCorrection: method, gasConstant: gasConstant).MeltingTemperature;
        }
        catch (ArgumentException)
        {
            return double.NaN; // degenerate duplex the reference rejects (Tm_NN ZeroDivisionError)
        }
    }

    // ---- NN internal-mismatch + dangling-end Tm (PRIMER-NNTM-001, opt-in extension) -------
    // Extends the perfect-match NN model to a probe–target duplex with internal mismatches,
    // terminal mismatches and/or a single unpaired dangling base at either end, exactly as
    // Biopython Tm_NN(seq, c_seq, shift, nn_table=DNA_NN4, tmm_table=DNA_TMM1, imm_table=DNA_IMM1,
    // de_table=DNA_DE1) does (the core ThermoConstants.CalculateNearestNeighborDuplex is a port of it).
    //
    // Convention of this API: the top strand is 5'→3'; the bottom strand is supplied 3'→5' (the
    // complement of the top read in the SAME left-to-right order, NOT the reverse complement), so
    // column i pairs top[i] with bottom[i]. A '.' as the first or last character of a strand marks
    // the missing partner of a single dangling base on the other strand. It maps onto Tm_NN as
    // seq = top without '.', c_seq = bottom without '.', shift = (leading '.' of top) − (leading '.'
    // of bottom). Tm_NN then scores, in this order: dangling ends (Bommarito et al. 2000), terminal
    // mismatches (SantaLucia & Peyret 2001, Biopython DNA_TMM1), the initiation terms — the terminal
    // A·T penalty is taken from the first and last base of the TOP strand without '.' (Tm_NN's
    // `ends = seq[0] + seq[-1]`) — and every remaining stack from the internal-mismatch table
    // (Allawi & SantaLucia 1997/1998, Peyret et al. 1999) or the Watson-Crick table.

    /// <summary>
    /// Computes the duplex ΔH° (kcal/mol) and ΔS° (cal/(K·mol), 1 M NaCl) of a probe–target DNA duplex
    /// that may contain internal single mismatches (Allawi &amp; SantaLucia 1997/1998; Peyret et al. 1999),
    /// terminal mismatches (SantaLucia &amp; Peyret 2001) and/or a single dangling end at either end
    /// (Bommarito et al. 2000), with the SantaLucia &amp; Hicks (2004) Watson–Crick parameters. Identical to
    /// Biopython <c>Tm_NN(seq, c_seq, shift, nn_table=DNA_NN4)</c>'s ΔH°/ΔS° (see the mapping in the
    /// remarks of <see cref="CalculateMeltingTemperatureNNMismatch"/>).
    /// </summary>
    /// <param name="topStrand">Top strand 5'→3'. May start/end with a single '.' marking a
    /// dangling base of the bottom strand.</param>
    /// <param name="bottomStrand">Bottom strand written 3'→5' (the complement of the top
    /// read left-to-right, NOT the reverse complement), so base i pairs with top base i.
    /// May start/end with a single '.' marking a dangling base of the top strand.</param>
    /// <returns>(ΔH°, ΔS°, IsSelfComplementary) or <c>null</c> if the strands are null,
    /// unequal length, shorter than two columns, contain a character other than A/C/G/T (case-insensitive)
    /// or a '.' that is not a single terminal marker facing a base, or contain a stack with no parameter
    /// (e.g. adjacent mismatches other than the tandem G·T motifs).</returns>
    public static (double DeltaH, double DeltaS, bool IsSelfComplementary)? CalculateNearestNeighborThermodynamicsMismatch(
        string topStrand, string bottomStrand)
    {
        var mapped = MapMismatchDuplex(topStrand, bottomStrand);
        if (mapped is null)
            return null;

        var (seq, cSeq, shift, selfComp) = mapped.Value;
        try
        {
            var (dH, dS) = ThermoConstants.CalculateNearestNeighborThermodynamics(
                seq, cSeq, shift, DesignNnParameterSet, selfComp, check: false, strict: true);
            return (dH, dS, selfComp);
        }
        catch (ArgumentException)
        {
            return null; // a neighbour pair with no thermodynamic parameter
        }
    }

    // Validates the column-aligned (top, bottom) pair and maps it to Tm_NN's (seq, c_seq, shift).
    private static (string Seq, string CSeq, int Shift, bool SelfComp)? MapMismatchDuplex(string top, string bottom)
    {
        if (top is null || bottom is null)
            return null;
        string t = top.ToUpperInvariant();
        string b = bottom.ToUpperInvariant();
        int n = t.Length;
        if (n != b.Length || n < 2)
            return null;

        for (int i = 0; i < n; i++)
        {
            char x = t[i], y = b[i];
            bool xDot = x == '.', yDot = y == '.';
            if ((xDot || yDot) && i != 0 && i != n - 1)
                return null; // a dangling marker must be terminal
            if (xDot && yDot)
                return null; // a column with no base at all
            if ((!xDot && x is not ('A' or 'C' or 'G' or 'T')) || (!yDot && y is not ('A' or 'C' or 'G' or 'T')))
                return null;
        }

        string seq = t.Replace(".", string.Empty, StringComparison.Ordinal);
        string cSeq = b.Replace(".", string.Empty, StringComparison.Ordinal);
        if (seq.Length < 1 || cSeq.Length < 1)
            return null;
        int shift = (t[0] == '.' ? 1 : 0) - (b[0] == '.' ? 1 : 0);

        // Symmetry term only for a fully paired, self-complementary duplex (Tm_NN's selfcomp flag).
        bool hasDangling = t.Contains('.') || b.Contains('.');
        bool selfComp = !hasDangling && OligoSymmetry.IsSymmetric(t)
                        && string.Equals(b, Complement(t), StringComparison.Ordinal);
        return (seq, cSeq, shift, selfComp);
    }

    /// <summary>
    /// Computes the design melting temperature (°C) for a probe–target DNA duplex that may contain
    /// internal and terminal mismatches and/or a single dangling end at either end, with the
    /// SantaLucia &amp; Hicks (2004) Watson–Crick parameters plus the internal-mismatch, terminal-mismatch and
    /// Bommarito (2000) dangling-end terms (<see cref="CalculateNearestNeighborThermodynamicsMismatch"/>), the
    /// same bimolecular Tm equation and salt corrections as <see cref="CalculateMeltingTemperatureNN"/>.
    /// <b>Opt-in extension</b>: a fully paired duplex through this path equals
    /// <see cref="CalculateMeltingTemperatureNN"/>.
    /// </summary>
    /// <remarks>
    /// Equals Biopython <c>Tm_NN(seq, c_seq=…, shift=…, nn_table=DNA_NN4, dnac1=dnac2=C_T/2, Na, Mg, dNTPs,
    /// saltcorr=0/5/6/7)</c> with R = 1.9872, where seq/c_seq are the strands without '.', shift = (leading '.'
    /// of the top) − (leading '.' of the bottom). The salt correction uses the top strand without '.'
    /// (its length and GC fraction), as Tm_NN does.
    /// </remarks>
    /// <param name="topStrand">Top strand 5'→3' (may carry a leading/trailing '.' dangling-end marker).</param>
    /// <param name="bottomStrand">Bottom strand 3'→5', aligned base-for-base under the top
    /// (complement direction, NOT reverse complement; may carry a '.' dangling-end marker).</param>
    /// <param name="strandConcentrationMolar">Total strand concentration C_T in mol/L (default 0.5 µM).</param>
    /// <param name="sodiumMolar">Monovalent cation concentration in mol/L (default 50 mM).</param>
    /// <param name="magnesiumMolar">[Mg²⁺] in mol/L (default 0; only used by the divalent mode).</param>
    /// <param name="dntpMolar">Total dNTP concentration in mol/L (default 0).</param>
    /// <param name="saltMode">Salt correction to apply (default Owczarzy2004Monovalent).</param>
    /// <returns>The NN Tm in °C, or <c>double.NaN</c> if the duplex is not computable
    /// (see <see cref="CalculateNearestNeighborThermodynamicsMismatch"/>).</returns>
    /// <exception cref="ArgumentOutOfRangeException">A non-positive C_T or [Na⁺], or a negative [Mg²⁺]/[dNTPs].</exception>
    public static double CalculateMeltingTemperatureNNMismatch(
        string topStrand,
        string bottomStrand,
        double strandConcentrationMolar = DefaultStrandConcentrationMolar,
        double sodiumMolar = ThermoConstants.DefaultNaConcentration,
        double magnesiumMolar = 0.0,
        double dntpMolar = 0.0,
        SaltCorrectionMode saltMode = SaltCorrectionMode.Owczarzy2004Monovalent)
    {
        ValidateNnConditions(strandConcentrationMolar, sodiumMolar, magnesiumMolar, dntpMolar);

        var thermo = CalculateNearestNeighborThermodynamicsMismatch(topStrand, bottomStrand);
        if (thermo is null)
            return double.NaN;

        var (dH, dS, selfComp) = thermo.Value;
        string topPaired = topStrand.ToUpperInvariant().Replace(".", string.Empty, StringComparison.Ordinal);
        return NnTm(dH, dS, topPaired, selfComp,
            strandConcentrationMolar, sodiumMolar, magnesiumMolar, dntpMolar, saltMode);
    }

    // ---- LNA (locked nucleic acid)-modified NN thermodynamics + Tm (PROBE-LNATM-001) ----------
    // A DNA oligonucleotide with LNA monomers at given INTERNAL positions, hybridised to a DNA strand
    // (perfect complement or a supplied 3'→5' target), modelled by one of two published LNA·DNA
    // nearest-neighbour models, exactly as the MELTING 5 reference implementation (Dumousseau et al.
    // 2012, BMC Bioinformatics 13:101; melting5.jar 5.2.0 as shipped in Bioconductor rmelting:
    // LockedAcidNNMethod / McTigue04LockedAcid / Owczarzy11LockedAcid / Owczarzy11TandemLockedAcid /
    // Owczarzy11SingleMismatchLockedAcid, NearestNeighborMode.getAppropriatePatternModel) realises them:
    //
    //  * Base DNA model (both LNA models): SantaLucia (1998) unified = Allawi & SantaLucia (1997) stacks,
    //    initiation per terminal A·T (+2.3/+4.1) or G·C (+0.1/−2.8) pair (Biopython DNA_NN3, MELTING
    //    "all97", the DNA/DNA default of MELTING). McTigue et al. (2004) derived their LNA increments on
    //    this unified set; Owczarzy et al. (2011)'s single-LNA ΔH° equal those increments added to it
    //    (e.g. TTL/AA −5.574 = −7.9 + 2.326). Steps without LNA: Watson–Crick or internal-mismatch
    //    (DNA_IMM1) value, looked up as Biopython Tm_NN (ThermoConstants.TryGetNearestNeighborDuplexStep).
    //  * McTigue, Peterson & Kahn (2004) Biochemistry 43:5388 (MELTING "mct04"): for an ISOLATED LNA, each
    //    of its two flanking steps = DNA stack + ΔΔH°/ΔΔS° increment (32 increments, McTigue2004lockedmn.xml).
    //    McTigue parameterised single internal LNAs only; a run of ≥ 2 consecutive LNAs is handled — as in
    //    MELTING (getAppropriatePatternModel → tandem/single-mismatch LNA model whatever -lck says) — by the
    //    Owczarzy (2011) tables below.
    //  * Owczarzy, You, Groth & Tataurov (2011) Biochemistry 50:9352 (MELTING "owc11", its default): every
    //    step containing an LNA takes a COMPLETE NN value (no DNA base added): one LNA in the step → 32
    //    single-LNA parameters (Owczarzy2011lockedmn.xml); two LNAs, both pairs Watson–Crick → 16
    //    consecutive-LNA parameters (Owczarzy2011lockedTandemmn.xml); two LNAs, a mismatch opposite one of
    //    them → LNA-mismatch parameters (Owczarzy2011lockedmmn.xml). A mismatched LNA is therefore
    //    computable only when both of its neighbours are LNAs (the paper's +X+M+Y triplets).
    //  * Terminal LNAs and terminal mismatches are not parameterised by either paper → not computable.
    //    MELTING 5.2.0 rejects terminal mismatches ("No method for terminal mismatches"), but it does evaluate a
    //    terminal LNA: its isApplicable guard (warning "not established for terminal locked nucleic acids", then
    //    false) compares the end base pair with the literal "L"/"-" and never fires, so it applies the internal-LNA
    //    doublet parameter at the duplex end (e.g. CLCATTGCTACC, 0.1 mM, 1 M Na⁺ → 66.6565 °C). That extrapolation
    //    outside the parameters' validated (internal) range is deliberately not reproduced (B07 F66).
    //  * Tm = 1000·ΔH° / (ΔS° + R·ln(C_T/x)) − 273.15 (x = 4, 1 if self-complementary) with the salt
    //    correction of saltMode (MELTING default for Na⁺ only = Owczarzy et al. 2004 Eq. 22 = Biopython
    //    method 6). MELTING's code uses R = 1.99 (NearestNeighborMode.computesMeltingTemperature); this
    //    class uses R = 1.9872 by default — pass gasConstant: 1.99 for bit-exact MELTING parity.
    //
    // Data anomalies kept exactly as in the reference implementation (the primary tables are behind
    // the ACS paywall; see docs/Validation/review-2026-09/B07.md F26): the mismatch file has no
    // "GLAL/CA" entry but a double-mismatch key "CLAL/CA", and lists "GLTL/TA" twice (MELTING's
    // HashMap keeps the last: −14.213/−40.041); an independent transcription (iCarrin/Bio_dpt
    // lna_tm.py) has ΔS° −21.535 for consecutive +T+C where MELTING has −21.735 (MELTING kept).
    // Reference cross-check (melting5.jar, Na = 1 M): CCATT(L)GCTACC, 1e-4 M — mct04 63.61426 °C,
    // owc11 63.48299 °C; GA(L)C(L)C 12.94323 °C (rmelting test-method.locked.R) — reproduced bit-exactly.

    /// <summary>Published LNA·DNA nearest-neighbour models for <see cref="CalculateMeltingTemperatureNNLna(string, IReadOnlyCollection{int}, LnaNearestNeighborModel, string?, double, double, double, double, SaltCorrectionMode, double)"/>.</summary>
    public enum LnaNearestNeighborModel
    {
        /// <summary>
        /// Owczarzy, You, Groth &amp; Tataurov (2011) Biochemistry 50:9352 — complete NN parameters for steps with
        /// one LNA, two consecutive LNAs, and an LNA·DNA mismatch inside an LNA triplet (MELTING "owc11", its default).
        /// </summary>
        Owczarzy2011,

        /// <summary>
        /// McTigue, Peterson &amp; Kahn (2004) Biochemistry 43:5388 — ΔΔH°/ΔΔS° increments for an isolated internal
        /// LNA added to the SantaLucia (1998) unified DNA stacks (MELTING "mct04"); runs of consecutive LNAs use
        /// the <see cref="Owczarzy2011"/> parameters, as MELTING does.
        /// </summary>
        McTigue2004
    }

    // The base DNA NN set of both LNA models: SantaLucia (1998) unified / Allawi & SantaLucia (1997).
    private const NnParameterSet LnaBaseParameterSet = NnParameterSet.AllawiSantaLucia1997;

    // LNA parameter tables (ΔH° kcal/mol, ΔS° cal/(K·mol), 1 M Na⁺), transcribed from the MELTING 5.2.0
    // data files (cal/mol ÷ 1000). Key = MELTING notation: top strand 5'→3' with "L" after each locked
    // base / the two opposite bottom-strand bases 3'→5' (e.g. "TTL/AA": step TT, 3' base locked).
    // McTigue (2004): increments ΔΔH°/ΔΔS° (McTigue2004lockedmn.xml).
    private static readonly Dictionary<string, (double DeltaH, double DeltaS)> McTigue2004LnaIncrements = new(StringComparer.Ordinal)
    {
        ["ALA/TT"] = (0.707, 2.5), ["ALT/TA"] = (2.282, 7.5),
        ["ALG/TC"] = (0.264, 2.6), ["ALC/TG"] = (1.131, 4.1),
        ["TLA/AT"] = (-0.046, 1.6), ["TLT/AA"] = (1.528, 5.3),
        ["TLG/AC"] = (-1.540, -3), ["TLC/AG"] = (1.893, 6.7),
        ["GLA/CT"] = (3.162, 10.5), ["GLT/CA"] = (-0.212, 0.1),
        ["GLG/CC"] = (-2.844, -6.7), ["GLC/CG"] = (-0.360, -0.3),
        ["CLA/GT"] = (1.049, 4.3), ["CLT/GA"] = (0.708, 4.2),
        ["CLG/GC"] = (0.785, 3.7), ["CLC/GG"] = (2.096, 8),
        ["AAL/TT"] = (0.992, 4.1), ["ATL/TA"] = (1.816, 6.9),
        ["AGL/TC"] = (-1.200, -1.8), ["ACL/TG"] = (2.890, 10.6),
        ["TAL/AT"] = (1.591, 5.3), ["TTL/AA"] = (2.326, 8.1),
        ["TGL/AC"] = (2.165, 7.2), ["TCL/AG"] = (0.609, 3.2),
        ["GAL/CT"] = (0.444, 2.9), ["GTL/CA"] = (-0.635, -0.3),
        ["GGL/CC"] = (-0.943, -0.9), ["GCL/CG"] = (-0.925, -1.1),
        ["CAL/GT"] = (1.358, 4.4), ["CTL/GA"] = (-1.671, -4.1),
        ["CGL/GC"] = (-0.276, -0.7), ["CCL/GG"] = (2.063, 7.6)
    };

    // Owczarzy (2011): complete NN parameters of a step with one LNA (Owczarzy2011lockedmn.xml).
    private static readonly Dictionary<string, (double DeltaH, double DeltaS)> Owczarzy2011LnaSingle = new(StringComparer.Ordinal)
    {
        ["ALA/TT"] = (-7.193, -19.723), ["ALC/TG"] = (-7.269, -18.336),
        ["ALG/TC"] = (-7.536, -18.387), ["ALT/TA"] = (-4.918, -12.943),
        ["CLA/GT"] = (-7.451, -18.38), ["CLC/GG"] = (-5.904, -11.904),
        ["CLG/GC"] = (-9.815, -23.491), ["CLT/GA"] = (-7.092, -16.825),
        ["GLA/CT"] = (-5.038, -11.656), ["GLC/CG"] = (-10.160, -24.651),
        ["GLG/CC"] = (-10.844, -26.58), ["GLT/CA"] = (-8.612, -22.327),
        ["TLA/AT"] = (-7.246, -19.738), ["TLC/AG"] = (-6.307, -15.515),
        ["TLG/AC"] = (-10.040, -25.744), ["TLT/AA"] = (-6.372, -16.902),
        ["AAL/TT"] = (-6.908, -18.135), ["ACL/TG"] = (-5.510, -11.824),
        ["AGL/TC"] = (-9.000, -22.826), ["ATL/TA"] = (-5.384, -13.537),
        ["CAL/GT"] = (-7.142, -18.333), ["CCL/GG"] = (-5.937, -12.335),
        ["CGL/GC"] = (-10.876, -27.918), ["CTL/GA"] = (-9.471, -25.07),
        ["GAL/CT"] = (-7.756, -19.302), ["GCL/CG"] = (-10.725, -25.511),
        ["GGL/CC"] = (-8.943, -20.833), ["GTL/CA"] = (-9.035, -22.742),
        ["TAL/AT"] = (-5.609, -16.019), ["TCL/AG"] = (-7.591, -19.031),
        ["TGL/AC"] = (-6.335, -15.537), ["TTL/AA"] = (-5.574, -14.149)
    };

    // Owczarzy (2011): two consecutive LNAs, Watson–Crick (Owczarzy2011lockedTandemmn.xml).
    private static readonly Dictionary<string, (double DeltaH, double DeltaS)> Owczarzy2011LnaConsecutive = new(StringComparer.Ordinal)
    {
        ["ALAL/TT"] = (-9.991, -27.175), ["ALCL/TG"] = (-11.389, -28.963),
        ["ALGL/TC"] = (-12.793, -31.607), ["ALTL/TA"] = (-14.703, -40.75),
        ["CLAL/GT"] = (-14.177, -35.498), ["CLCL/GG"] = (-15.399, -36.375),
        ["CLGL/GC"] = (-14.558, -35.239), ["CLTL/GA"] = (-15.737, -41.218),
        ["GLAL/CT"] = (-13.959, -35.097), ["GLCL/CG"] = (-16.109, -40.738),
        ["GLGL/CC"] = (-13.022, -29.673), ["GLTL/CA"] = (-17.361, -45.858),
        ["TLAL/AT"] = (-10.318, -26.108), ["TLCL/AG"] = (-9.166, -21.735),
        ["TLGL/AC"] = (-10.046, -22.591), ["TLTL/AA"] = (-10.419, -27.683)
    };

    // Owczarzy (2011): two consecutive LNAs with a mismatch opposite one of them (Owczarzy2011lockedmmn.xml;
    // 98 entries, "GLTL/TA" listed twice — the later value is kept, as MELTING's HashMap does).
    private static readonly Dictionary<string, (double DeltaH, double DeltaS)> Owczarzy2011LnaMismatch = new(StringComparer.Ordinal)
    {
        ["ALAL/AT"] = (-3.826, -13.109), ["ALCL/AG"] = (-2.367, -7.322),
        ["ALGL/AC"] = (-4.849, -13.007), ["ALTL/AA"] = (-5.049, -17.514),
        ["ALAL/TA"] = (-4.229, -15.16), ["CLAL/GA"] = (-5.878, -17.663),
        ["CLAL/CA"] = (-8.558, -23.976), ["TLAL/AA"] = (2.074, 3.446),
        ["CLAL/CT"] = (2.218, 4.75), ["CLCL/CG"] = (1.127, 1.826),
        ["CLGL/CC"] = (-10.903, -32.025), ["CLTL/CA"] = (-2.053, -10.517),
        ["ALCL/TC"] = (1.065, -1.403), ["CLCL/GC"] = (-9.522, -27.024),
        ["GLCL/CC"] = (-4.767, -14.897), ["TLCL/AC"] = (4.114, 9.258),
        ["GLAL/GT"] = (-2.920, -9.387), ["GLCL/GG"] = (-8.139, -21.784),
        ["GLGL/GC"] = (-5.149, -12.508), ["GLTL/GA"] = (-8.991, -27.311),
        ["ALGL/TG"] = (-4.980, -15.426), ["CLGL/GG"] = (-4.441, -12.158),
        ["GLGL/CG"] = (-13.505, -36.021), ["TLGL/AG"] = (-2.775, -9.286),
        ["TLAL/TT"] = (-3.744, -12.149), ["TLCL/TG"] = (-4.387, -13.52),
        ["TLGL/TC"] = (-6.346, -16.629), ["TLTL/TA"] = (-7.697, -25.049),
        ["ALTL/TT"] = (-4.207, -14.307), ["CLTL/GT"] = (-8.176, -22.962),
        ["GLTL/CT"] = (-7.241, -20.622), ["TLTL/AT"] = (-2.051, -7.055),
        ["ALAL/CT"] = (-1.362, -5.551), ["ALCL/CG"] = (-1.759, -6.511),
        ["ALGL/CC"] = (-6.549, -18.073), ["ALTL/CA"] = (-3.563, -14.105),
        ["ALAL/TC"] = (-2.078, -10.088), ["CLAL/GC"] = (-5.868, -16.952),
        ["GLAL/CC"] = (-8.477, -24.565), ["TLAL/AC"] = (2.690, 4.965),
        ["CLAL/AT"] = (-9.844, -29.673), ["CLCL/AG"] = (-3.761, -11.204),
        ["CLGL/AC"] = (-9.845, -27.316), ["CLTL/AA"] = (-3.389, -12.517),
        ["ALCL/TA"] = (0.753, -0.503), ["CLCL/GA"] = (-12.714, -35.555),
        ["GLCL/CA"] = (-12.658, -35.729), ["TLCL/AA"] = (-1.719, -7.023),
        ["ALAL/GT"] = (2.193, 4.374), ["ALCL/GG"] = (-8.453, -22.672),
        ["ALGL/GC"] = (-1.164, -2.532), ["ALTL/GA"] = (-7.418, -24.066),
        ["ALAL/TG"] = (-1.963, -9.013), ["CLAL/GG"] = (-8.712, -23.779),
        ["GLAL/CG"] = (-7.875, -21.661), ["TLAL/AG"] = (3.207, 7.156),
        ["GLAL/AT"] = (-2.914, -9.402), ["GLCL/AG"] = (-9.131, -25.347),
        ["GLGL/AC"] = (-2.154, -3.871), ["GLTL/AA"] = (-8.515, -26.313),
        ["ALGL/TA"] = (-6.691, -21.148), ["CLGL/GA"] = (-3.960, -10.588),
        ["GLGL/CA"] = (-12.898, -34.656), ["TLGL/AA"] = (0.334, -0.44),
        ["CLAL/TT"] = (0.382, -0.579), ["CLCL/TG"] = (-2.716, -8),
        ["CLGL/TC"] = (-10.363, -29.315), ["CLTL/TA"] = (-5.783, -20.173),
        ["ALCL/TT"] = (-0.692, -5.278), ["CLCL/GT"] = (-10.288, -28.503),
        ["GLCL/CT"] = (-9.062, -26.356), ["TLCL/AT"] = (2.073, 3.968),
        ["TLAL/CT"] = (-5.485, -17.347), ["TLCL/CG"] = (1.451, 1.556),
        ["TLGL/CC"] = (-7.213, -20.128), ["TLTL/CA"] = (-2.397, -11.371),
        ["ALTL/TC"] = (-0.633, -5.801), ["CLTL/GC"] = (-6.868, -21),
        ["GLTL/CC"] = (-5.853, -16.643), ["TLTL/AC"] = (0.211, -1.446),
        ["GLAL/TT"] = (-5.551, -15.398), ["GLCL/TG"] = (-14.943, -40.148),
        ["GLGL/TC"] = (-8.110, -18.349), ["GLTL/TA"] = (-14.213, -40.041),
        ["ALGL/TT"] = (-7.130, -20.786), ["CLGL/GT"] = (-14.862, -39.43),
        ["GLGL/CT"] = (-14.622, -37.51), ["TLGL/AT"] = (-6.703, -18.111),
        ["TLAL/GT"] = (-4.612, -14.039), ["TLCL/GG"] = (-9.798, -26.406),
        ["TLGL/GC"] = (-4.519, -11.065), ["TLTL/GA"] = (-4.523, -15.693),
        ["ALTL/TG"] = (-2.364, -8.834), ["CLTL/GG"] = (-11.396, -30.732),
        ["GLTL/CG"] = (-6.233, -15.933), ["TLTL/AG"] = (-2.960, -9.305)
    };

    /// <summary>
    /// Duplex ΔH° (kcal/mol) and ΔS° (cal/(K·mol), 1 M Na⁺) of a DNA oligonucleotide with <b>internal</b> LNA
    /// monomers at <paramref name="lnaPositions"/> paired with its perfect DNA complement, by the default
    /// <see cref="LnaNearestNeighborModel.Owczarzy2011"/> LNA·DNA nearest-neighbour model on the SantaLucia (1998)
    /// unified DNA parameters (MELTING 5 default). See
    /// <see cref="CalculateNearestNeighborThermodynamicsLna(string, IReadOnlyCollection{int}, LnaNearestNeighborModel, string?)"/>.
    /// </summary>
    /// <param name="sequence">DNA sequence (5'→3'; case-insensitive), ≥ 2 ACGT bases.</param>
    /// <param name="lnaPositions">Zero-based LNA positions (order/duplicates tolerated); terminal or out-of-range → not computable.</param>
    /// <returns>(ΔH°, ΔS°, IsSelfComplementary) or <c>null</c> when not computable.</returns>
    public static (double DeltaH, double DeltaS, bool IsSelfComplementary)? CalculateNearestNeighborThermodynamicsLna(
        string sequence,
        IReadOnlyCollection<int> lnaPositions) =>
        CalculateNearestNeighborThermodynamicsLna(sequence, lnaPositions, LnaNearestNeighborModel.Owczarzy2011);

    /// <summary>
    /// Duplex ΔH° (kcal/mol) and ΔS° (cal/(K·mol), 1 M Na⁺) of a DNA oligonucleotide carrying <b>internal</b> LNA
    /// (locked nucleic acid) monomers, hybridised to a DNA strand, by a published LNA·DNA nearest-neighbour model
    /// (see <see cref="LnaNearestNeighborModel"/>) on the SantaLucia (1998) unified DNA parameters (Allawi &amp;
    /// SantaLucia 1997; initiation per terminal A·T / G·C pair) — the MELTING 5 implementation, reproduced
    /// bit-exactly. Steps without an LNA use the Watson–Crick or (for an internal DNA mismatch) the Allawi /
    /// SantaLucia / Peyret internal-mismatch parameters.
    /// </summary>
    /// <param name="sequence">LNA-modified oligonucleotide (DNA letters, 5'→3'; case-insensitive), ≥ 2 ACGT bases.</param>
    /// <param name="lnaPositions">Zero-based positions of the LNA monomers in <paramref name="sequence"/>; order and
    /// duplicates are tolerated. A terminal (0 or length − 1) or out-of-range position is not parameterised (neither
    /// paper covers terminal LNAs; MELTING 5.2.0 extrapolates them only because its terminal-LNA guard never fires).</param>
    /// <param name="model">LNA nearest-neighbour model.</param>
    /// <param name="target">The DNA strand opposite <paramref name="sequence"/>, written 3'→5' (base i pairs with
    /// sequence[i]; same length); <c>null</c> = the perfect complement. Internal mismatches are allowed opposite a DNA
    /// base (DNA_IMM1) or opposite the central LNA of three consecutive LNAs (Owczarzy 2011); the two terminal pairs
    /// must be Watson–Crick.</param>
    /// <returns>(ΔH°, ΔS°, IsSelfComplementary) or <c>null</c> when not computable: empty / &lt; 2 bases / non-ACGT
    /// strand, target of another length, terminal or out-of-range LNA, terminal mismatch, or a step with no
    /// published parameter (e.g. a mismatch opposite an isolated LNA, any mismatch with
    /// <see cref="LnaNearestNeighborModel.McTigue2004"/> at an isolated LNA).</returns>
    /// <exception cref="ArgumentNullException"><paramref name="lnaPositions"/> is null.</exception>
    public static (double DeltaH, double DeltaS, bool IsSelfComplementary)? CalculateNearestNeighborThermodynamicsLna(
        string sequence,
        IReadOnlyCollection<int> lnaPositions,
        LnaNearestNeighborModel model,
        string? target = null)
    {
        ArgumentNullException.ThrowIfNull(lnaPositions);
        if (model is not (LnaNearestNeighborModel.Owczarzy2011 or LnaNearestNeighborModel.McTigue2004))
            throw new ArgumentOutOfRangeException(nameof(model), model, "Unknown LNA nearest-neighbour model.");
        if (string.IsNullOrEmpty(sequence) || sequence.Length < 2)
            return null;

        string seq = sequence.ToUpperInvariant();
        if (!IsAcgtOnly(seq))
            return null;
        string bottom = target is null ? Complement(seq) : target.ToUpperInvariant();
        if (bottom.Length != seq.Length || !IsAcgtOnly(bottom))
            return null;

        int n = seq.Length;
        var locked = new bool[n];
        foreach (int pos in lnaPositions)
        {
            // Terminal LNAs are not parameterised by either model (MELTING's intended isApplicable rejection).
            if (pos <= 0 || pos >= n - 1)
                return null;
            locked[pos] = true;
        }

        // Terminal pairs must be Watson–Crick (no terminal-mismatch LNA parameters).
        if (!IsWatsonCrickPair(seq[0], bottom[0]) || !IsWatsonCrickPair(seq[n - 1], bottom[n - 1]))
            return null;

        // An LNA-modified strand paired with an unmodified DNA strand is never a symmetric duplex (the two
        // strands differ chemically), so only an LNA-free self-complementary sequence with its own complement
        // is treated as self-complementary (symmetry term, x = 1); MELTING likewise rejects -self with LNAs.
        bool selfComp = target is null && !locked.Contains(true) && OligoSymmetry.IsSymmetric(seq);
        var init = ThermoConstants.GetNearestNeighborInitiation(LnaBaseParameterSet);
        double dH = init.Initiation.DeltaH, dS = init.Initiation.DeltaS;
        var oneOrAll = seq.Any(c => c is 'G' or 'C') ? init.OneGC : init.AllAT;
        dH += oneOrAll.DeltaH; dS += oneOrAll.DeltaS;
        foreach (char end in new[] { seq[0], seq[n - 1] })
        {
            var t = end is 'A' or 'T' ? init.TerminalAT : init.TerminalGC;
            dH += t.DeltaH; dS += t.DeltaS;
        }
        if (seq[0] == 'T') { dH += init.FiveTerminalTA.DeltaH; dS += init.FiveTerminalTA.DeltaS; }
        if (seq[n - 1] == 'A') { dH += init.FiveTerminalTA.DeltaH; dS += init.FiveTerminalTA.DeltaS; }

        for (int i = 0; i < n - 1; i++)
        {
            string top = seq.Substring(i, 2);
            string bot = bottom.Substring(i, 2);
            int lockedCount = (locked[i] ? 1 : 0) + (locked[i + 1] ? 1 : 0);

            if (lockedCount == 0)
            {
                if (!ThermoConstants.TryGetNearestNeighborDuplexStep(LnaBaseParameterSet, top, bot, out var p))
                    return null;
                dH += p.DeltaH; dS += p.DeltaS;
                continue;
            }

            bool isolated = (locked[i] && IsIsolatedLna(locked, i)) || (locked[i + 1] && IsIsolatedLna(locked, i + 1));
            if (model == LnaNearestNeighborModel.McTigue2004 && isolated)
            {
                // McTigue (2004): DNA stack + increment (Watson–Crick steps only).
                if (!IsWatsonCrickPair(top[0], bot[0]) || !IsWatsonCrickPair(top[1], bot[1])
                    || !ThermoConstants.TryGetNearestNeighborStack(LnaBaseParameterSet, top, out var stack)
                    || !McTigue2004LnaIncrements.TryGetValue(LnaKey(top, bot, locked[i], locked[i + 1]), out var inc))
                    return null;
                dH += stack.DeltaH + inc.DeltaH;
                dS += stack.DeltaS + inc.DeltaS;
                continue;
            }

            // Owczarzy (2011): complete parameter of the LNA-containing step.
            bool watsonCrick = IsWatsonCrickPair(top[0], bot[0]) && IsWatsonCrickPair(top[1], bot[1]);
            var table = Owczarzy2011LnaSingle;
            if (lockedCount == 2)
                table = watsonCrick ? Owczarzy2011LnaConsecutive : Owczarzy2011LnaMismatch;
            if (!table.TryGetValue(LnaKey(top, bot, locked[i], locked[i + 1]), out var v))
                return null;
            dH += v.DeltaH; dS += v.DeltaS;
        }

        if (selfComp)
        {
            dH += init.Symmetry.DeltaH; dS += init.Symmetry.DeltaS;
        }
        return (dH, dS, selfComp);

        static bool IsIsolatedLna(bool[] l, int p) =>
            (p == 0 || !l[p - 1]) && (p == l.Length - 1 || !l[p + 1]);

        static string LnaKey(string top, string bot, bool lock0, bool lock1) =>
            string.Concat(top[0].ToString(), lock0 ? "L" : "", top[1].ToString(), lock1 ? "L" : "", "/", bot);
    }

    /// <summary>
    /// LNA-adjusted nearest-neighbour Tm (°C) of a DNA oligonucleotide with <b>internal</b> LNA monomers paired with its
    /// perfect DNA complement, by the default <see cref="LnaNearestNeighborModel.Owczarzy2011"/> model (MELTING 5
    /// default) at the stated conditions. See
    /// <see cref="CalculateMeltingTemperatureNNLna(string, IReadOnlyCollection{int}, LnaNearestNeighborModel, string?, double, double, double, double, SaltCorrectionMode, double)"/>.
    /// </summary>
    /// <param name="sequence">DNA sequence (5'→3'), ≥ 2 ACGT bases.</param>
    /// <param name="lnaPositions">Zero-based internal LNA positions.</param>
    /// <param name="strandConcentrationMolar">Total strand concentration C_T in mol/L (default 0.5 µM).</param>
    /// <param name="sodiumMolar">Monovalent cation concentration in mol/L (default 50 mM).</param>
    /// <param name="magnesiumMolar">[Mg²⁺] in mol/L (default 0; only used by the divalent mode).</param>
    /// <param name="dntpMolar">Total dNTP concentration in mol/L (default 0).</param>
    /// <param name="saltMode">Salt correction (default Owczarzy2004Monovalent, as MELTING for Na⁺ only).</param>
    /// <returns>Tm in °C, or <c>double.NaN</c> when not computable.</returns>
    /// <exception cref="ArgumentOutOfRangeException">C_T ≤ 0, [Na⁺] ≤ 0, [Mg²⁺] &lt; 0 or [dNTP] &lt; 0.</exception>
    public static double CalculateMeltingTemperatureNNLna(
        string sequence,
        IReadOnlyCollection<int> lnaPositions,
        double strandConcentrationMolar = DefaultStrandConcentrationMolar,
        double sodiumMolar = ThermoConstants.DefaultNaConcentration,
        double magnesiumMolar = 0.0,
        double dntpMolar = 0.0,
        SaltCorrectionMode saltMode = SaltCorrectionMode.Owczarzy2004Monovalent) =>
        CalculateMeltingTemperatureNNLna(sequence, lnaPositions, LnaNearestNeighborModel.Owczarzy2011, null,
            strandConcentrationMolar, sodiumMolar, magnesiumMolar, dntpMolar, saltMode);

    /// <summary>
    /// LNA-adjusted nearest-neighbour melting temperature (°C) of an LNA-modified DNA oligonucleotide hybridised to a
    /// DNA strand: ΔH°/ΔS° from
    /// <see cref="CalculateNearestNeighborThermodynamicsLna(string, IReadOnlyCollection{int}, LnaNearestNeighborModel, string?)"/>,
    /// Tm = 1000·ΔH° / (ΔS° + R·ln(C_T/x)) − 273.15 (x = 4; 1 for a self-complementary duplex) with the salt correction
    /// of <paramref name="saltMode"/> evaluated on the DNA letters of <paramref name="sequence"/> (the canonical
    /// <see cref="ThermoConstants.CalculateNearestNeighborTmFromThermodynamics"/>). With
    /// <c>gasConstant: 1.99</c> and <see cref="SaltCorrectionMode.Owczarzy2004Monovalent"/> this equals MELTING 5
    /// (<c>-H dnadna -P C_T -E Na=…</c>, <c>-lck mct04|owc11</c>) exactly.
    /// </summary>
    /// <param name="sequence">LNA-modified oligonucleotide (DNA letters, 5'→3'), ≥ 2 ACGT bases.</param>
    /// <param name="lnaPositions">Zero-based internal LNA positions.</param>
    /// <param name="model">LNA nearest-neighbour model.</param>
    /// <param name="target">Opposite DNA strand 3'→5' (same length), or <c>null</c> for the perfect complement.</param>
    /// <param name="strandConcentrationMolar">Total strand concentration C_T in mol/L (default 0.5 µM).</param>
    /// <param name="sodiumMolar">Monovalent cation concentration in mol/L (default 50 mM).</param>
    /// <param name="magnesiumMolar">[Mg²⁺] in mol/L (default 0; only used by the divalent mode).</param>
    /// <param name="dntpMolar">Total dNTP concentration in mol/L (default 0).</param>
    /// <param name="saltMode">Salt correction (default Owczarzy2004Monovalent).</param>
    /// <param name="gasConstant">R in cal/(K·mol) (default 1.9872, SantaLucia &amp; Hicks 2004; MELTING uses 1.99).</param>
    /// <returns>Tm in °C, or <c>double.NaN</c> when not computable.</returns>
    /// <exception cref="ArgumentOutOfRangeException">C_T ≤ 0, [Na⁺] ≤ 0, [Mg²⁺] &lt; 0, [dNTP] &lt; 0, R ≤ 0 or an unknown model.</exception>
    /// <exception cref="ArgumentNullException"><paramref name="lnaPositions"/> is null.</exception>
    public static double CalculateMeltingTemperatureNNLna(
        string sequence,
        IReadOnlyCollection<int> lnaPositions,
        LnaNearestNeighborModel model,
        string? target = null,
        double strandConcentrationMolar = DefaultStrandConcentrationMolar,
        double sodiumMolar = ThermoConstants.DefaultNaConcentration,
        double magnesiumMolar = 0.0,
        double dntpMolar = 0.0,
        SaltCorrectionMode saltMode = SaltCorrectionMode.Owczarzy2004Monovalent,
        double gasConstant = GasConstant)
    {
        ValidateNnConditions(strandConcentrationMolar, sodiumMolar, magnesiumMolar, dntpMolar);
        if (!(gasConstant > 0) || double.IsInfinity(gasConstant))
            throw new ArgumentOutOfRangeException(nameof(gasConstant), gasConstant, "R must be a positive finite number.");

        var thermo = CalculateNearestNeighborThermodynamicsLna(sequence, lnaPositions, model, target);
        if (thermo is null)
            return double.NaN;

        var (dH, dS, selfComp) = thermo.Value;
        return NnTm(dH, dS, sequence.ToUpperInvariant(), selfComp,
            strandConcentrationMolar, sodiumMolar, magnesiumMolar, dntpMolar, saltMode, gasConstant);
    }

    // ---- DNA hairpin folding + secondary-structure (hairpin) Tm (PRIMER-TM-001, opt-in) ----
    // Finds the most stable intramolecular hairpin (a single stem closing one hairpin loop)
    // of a DNA oligo and computes its ΔH°/ΔS°/ΔG°37 and unimolecular melting temperature.
    // The perfect-match CalculateMeltingTemperatureNN and the default Wallace/Marmur-Doty Tm
    // are UNCHANGED; this is a new opt-in capability.
    //
    // Model (SantaLucia & Hicks 2004, Annu Rev Biophys 33:415, "Hairpin Loops", Eqs 8–11):
    //   ΔG°37(hairpin) = Σ stem NN stacks (Table 1)  +  ΔG°37(hairpin loop of N) (Table 4)
    //   "To compute the stability of a complete hairpin + stem, one simply adds the
    //    salt-corrected base pair NN contributions (Table 1; Equation 3) to the loop energy
    //    from Equations 8–10." (paper, p.428).
    //   The stem contributes only its nearest-neighbour STACKS (Table 1 propagation terms) —
    //   the bimolecular duplex-initiation term (+0.2/−5.7) is a TWO-strand nucleation cost and
    //   does NOT apply to a unimolecular hairpin (the loop-initiation term is the nucleation
    //   cost instead). Loop ΔH° = 0 for every loop size, and the loop ΔS° increment is
    //   ΔS° = −ΔG°37 × 1000 / 310.15 (Table 4 note: "ΔS° = ΔG°37 × 1000/310.15"; the loop is
    //   destabilising so ΔG°37 > 0 → ΔS° < 0).
    //   Two-state hairpin Tm (Eq 11) is UNIMOLECULAR/concentration-independent:
    //       Tm = ΔH° × 1000 / ΔS° − 273.15        (NO R·ln(C_T/x) strand-concentration term).
    //
    // NOT bundled (honest residual): the supplementary triloop/tetraloop bonus tables (length-3
    // and length-4 special loops) and the terminal-mismatch increment (the first mismatch stack
    // closing loops of length ≥4) are separate Annual-Reviews supplementary tables not embedded
    // here. They are exposed as an OPT-IN caller-supplied additive ΔG°37/ΔH° adjustment
    // (default 0) so a caller who has those tables can supply the increment; without it the
    // result is the stem-stack + loop-initiation core, which is exact and fully sourced.
    // The complete model (thal.c: terminal mismatches, special loops, bulges/internal loops,
    // dangles, Mg²⁺/dNTP) is selected with StructureModel.Ntthal on the model overloads below,
    // which delegate to the canonical ntthal engine (audit round 3, A3-13).

    /// <summary>Hairpin loop ΔG°37 increment (kcal/mol, 1 M NaCl) by loop size (number of
    /// unpaired loop nucleotides). SantaLucia &amp; Hicks (2004) Table 4 "Hairpin loops" column;
    /// sizes 3–30 are tabulated. ΔH° = 0 for all sizes; ΔS° = −ΔG°37·1000/310.15.</summary>
    private static readonly Dictionary<int, double> HairpinLoopInitiationDeltaG = new()
    {
        [3] = 3.5, [4] = 3.5, [5] = 3.3, [6] = 4.0, [7] = 4.2, [8] = 4.3, [9] = 4.5,
        [10] = 4.6, [12] = 5.0, [14] = 5.1, [16] = 5.3, [18] = 5.5, [20] = 5.7,
        [25] = 6.1, [30] = 6.3
    };

    // SantaLucia & Hicks (2004): minimum sterically allowed hairpin loop size is 3 nt
    // ("Hairpin loops with lengths shorter than 3 are sterically prohibited.").
    private const int MinHairpinLoopSize = 3;

    // Reference temperature for the ΔG°37 ↔ ΔS° conversion (310.15 K = 37 °C).
    // SantaLucia & Hicks (2004) Table 4 note: ΔS° = ΔG°37 × 1000/310.15.
    private const double ReferenceTemperatureKelvin = 310.15;

    // Jacobson-Stockmayer entropic extrapolation coefficient for loop sizes beyond the
    // tabulated lengths. SantaLucia & Hicks (2004) Eq. 7: ΔG°37(loop-n) =
    // ΔG°37(loop-x) + 2.44·R·310.15·ln(n/x); the 2.44 coefficient is from recent DNA
    // kinetics measurements (ref 22), preferred over the older 1.75.
    private const double JacobsonStockmayerCoefficient = 2.44;

    /// <summary>
    /// Hairpin loop ΔG°37 (kcal/mol) for a loop of <paramref name="loopSize"/> unpaired
    /// nucleotides. Tabulated sizes return Table 4 directly; non-tabulated sizes are filled by
    /// the Jacobson-Stockmayer extrapolation from the largest tabulated size ≤ n
    /// (SantaLucia &amp; Hicks 2004 Eq. 7). Loop sizes &lt; 3 are sterically prohibited.
    /// </summary>
    private static double HairpinLoopDeltaG(int loopSize)
    {
        if (HairpinLoopInitiationDeltaG.TryGetValue(loopSize, out double dg))
            return dg;

        // Jacobson-Stockmayer from the largest tabulated x ≤ loopSize.
        int x = 0;
        foreach (int size in HairpinLoopInitiationDeltaG.Keys)
            if (size <= loopSize && size > x) x = size;

        return HairpinLoopInitiationDeltaG[x]
               + JacobsonStockmayerCoefficient * GasConstant * ReferenceTemperatureKelvin
                 * 1e-3 * Math.Log((double)loopSize / x);
    }

    private static bool IsWatsonCrickPair(char a, char b) =>
        (a == 'A' && b == 'T') || (a == 'T' && b == 'A') ||
        (a == 'G' && b == 'C') || (a == 'C' && b == 'G');

    /// <summary>
    /// The most stable intramolecular DNA hairpin found in an oligo: the closing stem span,
    /// stem length (base pairs), loop size, and the hairpin ΔH° (kcal/mol), ΔS° (cal/(K·mol)),
    /// and ΔG°37 (kcal/mol).
    /// </summary>
    /// <param name="StemStart">5'-most index (0-based) of the stem on the input strand.</param>
    /// <param name="StemEnd">3'-most index (0-based) of the stem on the input strand.</param>
    /// <param name="StemLength">Number of base pairs in the stem.</param>
    /// <param name="LoopSize">Number of unpaired loop nucleotides closed by the stem.</param>
    /// <param name="DeltaH">Hairpin ΔH° in kcal/mol.</param>
    /// <param name="DeltaS">Hairpin ΔS° in cal/(K·mol).</param>
    /// <param name="DeltaG37">Hairpin ΔG°37 in kcal/mol (negative = stable).</param>
    public readonly record struct HairpinResult(
        int StemStart, int StemEnd, int StemLength, int LoopSize,
        double DeltaH, double DeltaS, double DeltaG37);

    /// <summary>
    /// Finds the most stable (minimum ΔG°37) intramolecular DNA <b>hairpin</b> — a single
    /// Watson-Crick stem closing one hairpin loop — in <paramref name="sequence"/>, using the
    /// SantaLucia &amp; Hicks (2004) Table 1 nearest-neighbour stem stacks and their
    /// Table 4 hairpin-loop initiation increments. <b>Opt-in</b>: the duplex Tm methods
    /// are unchanged. Every stem length from <paramref name="minStemLength"/> up to the maximal
    /// Watson–Crick extension of each closing pair is a candidate as long as it closes a loop of
    /// ≥ 3 nt (so GGGGCCCC folds into a 2-bp stem + 3-nt loop). Returns <c>null</c> when the sequence is empty, contains a non-ACGT
    /// character, or admits no hairpin at all (no stem of ≥ 2 bp can close a loop of ≥ 3 nt,
    /// e.g. a homopolymer such as poly-A). The full thal.c hairpin model is the
    /// <see cref="FindMostStableHairpin(string, StructureModel, double?, double?, double?, double?, int?)"/> overload with
    /// <see cref="StructureModel.Ntthal"/>.
    /// <para>
    /// Model: ΔG°37 = Σ stem NN stacks (Table 1) + ΔG°37(loop of N) (Table 4); the bimolecular
    /// duplex-initiation term is intentionally excluded for this unimolecular structure. Loop
    /// ΔH° = 0; loop ΔS° = −ΔG°37·1000/310.15. The supplementary triloop/tetraloop and
    /// terminal-mismatch increments are not bundled (see <paramref name="loopBonusDeltaG37"/>).
    /// </para>
    /// </summary>
    /// <param name="sequence">DNA oligo (5'→3').</param>
    /// <param name="minStemLength">Minimum stem length in base pairs (default 2 → at least one
    /// NN stack). Must be ≥ 2.</param>
    /// <param name="loopBonusDeltaG37">Optional caller-supplied additive ΔG°37 increment
    /// (kcal/mol) for the terminal-mismatch / special triloop-tetraloop bonus that is NOT
    /// bundled (default 0). Added to the loop free energy; its ΔS° contribution follows the
    /// same −ΔG·1000/310.15 rule (ΔH° contribution 0), consistent with the Table 4 loop model.</param>
    /// <returns>The most stable hairpin, or <c>null</c> if none exists / invalid input.</returns>
    public static HairpinResult? FindMostStableHairpin(
        string sequence,
        int minStemLength = 2,
        double loopBonusDeltaG37 = 0.0)
    {
        if (string.IsNullOrEmpty(sequence) || minStemLength < 2)
            return null;

        string seq = sequence.ToUpperInvariant();
        int n = seq.Length;
        foreach (char c in seq)
            if (c is not ('A' or 'C' or 'G' or 'T'))
                return null; // non-ACGT base present

        HairpinResult? best = null;

        // For every candidate outermost closing pair (i, j), extend the stem inward as far as
        // Watson-Crick pairing allows, then close the remaining inner bases as a hairpin loop.
        for (int i = 0; i < n; i++)
        {
            for (int j = i + 1; j < n; j++)
            {
                if (!IsWatsonCrickPair(seq[i], seq[j]))
                    continue;

                double dH = 0.0, dS = 0.0;
                // Extend the stem from the outermost pair (i, j) inward. EVERY stem length that closes a
                // loop of ≥ 3 nt is a candidate (not only the maximal Watson–Crick extension): a stem that
                // could pair further but would then close a < 3-nt loop (e.g. GGGGCCCC, GCGCGCGCGC) still
                // forms the shorter stem + loop, as the loop bases are simply unpaired (SantaLucia & Hicks
                // 2004 Eqs 8–10; INV-01/INV-02).
                int a = i, b = j;
                int stemPairs = 0;
                while (a < b && IsWatsonCrickPair(seq[a], seq[b]))
                {
                    if (stemPairs > 0)
                    {
                        // NN stack between pair (a-1, b+1) and (a, b): key = 5'-strand dinucleotide seq[a-1..a].
                        string step = seq.Substring(a - 1, 2);
                        if (!TryGetDesignStack(step, out var p))
                            break;
                        dH += p.DeltaH;
                        dS += p.DeltaS;
                    }
                    stemPairs++;

                    // Innermost pair is (a, b); the loop is the bases strictly between them.
                    int loopSize = b - a - 1;
                    if (loopSize < MinHairpinLoopSize)
                        break; // a longer stem only shrinks the loop further

                    if (stemPairs >= minStemLength)
                    {
                        double loopDg = HairpinLoopDeltaG(loopSize) + loopBonusDeltaG37;
                        // Loop ΔH° = 0; loop ΔS° = −ΔG°37·1000/310.15 (destabilising loop).
                        double loopDs = -loopDg * 1000.0 / ReferenceTemperatureKelvin;

                        double totalDh = dH;                       // loop ΔH° contribution is 0
                        double totalDs = dS + loopDs;
                        double dG37 = totalDh - ReferenceTemperatureKelvin * totalDs / 1000.0;

                        if (best is null || dG37 < best.Value.DeltaG37)
                            best = new HairpinResult(i, j, stemPairs, loopSize, totalDh, totalDs, dG37);
                    }

                    a++;
                    b--;
                }
            }
        }

        return best;
    }

    /// <summary>
    /// Computes the secondary-structure (hairpin) melting temperature (°C) of a DNA oligo:
    /// finds its most stable intramolecular hairpin (<see cref="FindMostStableHairpin"/>) and
    /// returns the <b>unimolecular</b> two-state Tm = ΔH°·1000/ΔS° − 273.15
    /// (SantaLucia &amp; Hicks 2004, Eq. 11). A hairpin is intramolecular, so the Tm is
    /// concentration-independent: there is <b>no</b> R·ln(C_T/x) strand-concentration term.
    /// <b>Opt-in</b>: the duplex Tm methods (<see cref="CalculateMeltingTemperatureNN"/>) and the
    /// default <see cref="CalculateMeltingTemperature(string)"/> are unchanged.
    /// </summary>
    /// <param name="sequence">DNA oligo (5'→3').</param>
    /// <param name="minStemLength">Minimum stem length in base pairs (default 2).</param>
    /// <param name="loopBonusDeltaG37">Optional caller-supplied terminal-mismatch / special-loop
    /// ΔG°37 increment (default 0; not bundled — see <see cref="FindMostStableHairpin"/>).</param>
    /// <returns>The hairpin Tm in °C, or <c>double.NaN</c> if no hairpin exists / invalid input.</returns>
    public static double CalculateHairpinMeltingTemperature(
        string sequence,
        int minStemLength = 2,
        double loopBonusDeltaG37 = 0.0)
    {
        var hairpin = FindMostStableHairpin(sequence, minStemLength, loopBonusDeltaG37);
        if (hairpin is null)
            return double.NaN;

        var h = hairpin.Value;
        // Unimolecular: NO concentration term (Eq. 11).
        return (h.DeltaH * 1000.0) / h.DeltaS - KelvinOffset;
    }

    // ---- Self-dimer / hetero-dimer (intermolecular) Tm via thermodynamic alignment ----
    // PRIMER-TM-001, opt-in. Finds the most stable INTERMOLECULAR antiparallel duplex
    // between two oligonucleotides (self-dimer = an oligo against a second copy of itself;
    // hetero/cross-dimer = two different oligos) and returns its NN ΔH°/ΔS° and the
    // bimolecular Tm. Reuses the SantaLucia & Hicks (2004) Table 1 NN stacking table
    // (ThermoConstants NnParameterSet.SantaLuciaHicks2004), the terminal-A·T penalty, the duplex-initiation term and the
    // 0.368 entropy salt coefficient already used by CalculateMeltingTemperatureNN.
    // The duplex Tm / hairpin Tm / default Tm methods and their defaults are UNCHANGED.
    //
    // Model (gapless contiguous-WC scorer on the ntthal duplex terms — SantaLucia & Hicks 2004 unified NN;
    // the full ntthal DP with mismatches/loops/bulges/overhangs is NtthalDimer via CalculateDimerThermodynamicsNtthal):
    //   For each gapless antiparallel offset of strand2 (read 3'→5') under strand1 (5'→3'),
    //   each maximal contiguous run of Watson-Crick pairs (≥ 1 NN stack) is a candidate
    //   duplex with
    //     ΔH° = ΔH°_init + Σ stacks ΔH° + ΔH°_AT-penalty(per A·T-closed end),
    //     ΔS° = ΔS°_init + Σ stacks ΔS° + ΔS°_AT-penalty(per A·T-closed end)
    //           + 0.368·N_stacks·ln[Na⁺]   (salt correction baked into ΔS°, ntthal saltCorrectS),
    //   and bimolecular
    //     Tm = ΔH°·1000 / (ΔS° + R·ln(C_T / x)) − 273.15,
    //   with x = 1 when BOTH oligos are reverse-complement palindromes (ntthal symmetry_thermo),
    //   else x = 4. ntthal keeps the candidate with the highest Tm — so does this method.
    //
    // Sources (retrieved & extracted this session, 2026-06-25):
    //   SantaLucia J, Hicks D (2004) Annu Rev Biophys 33:415-440 — unified NN parameters
    //     (Table 1, ThermoConstants NnParameterSet.SantaLuciaHicks2004) + the bimolecular Tm Eq. 3 + Eq. 5 entropy
    //     salt correction (0.368 coefficient).
    //   Untergasser A et al. (2012) Nucleic Acids Res 40:e115 (Primer3 2.0) — the ntthal
    //     thermodynamic-alignment engine for oligo dimers.
    //   Primer3 `thal.c` (primer3-py vendored libprimer3, retrieved
    //     https://raw.githubusercontent.com/libnano/primer3-py/master/primer3/src/libprimer3/thal.c):
    //     dplx_init_H=200 cal, dplx_init_S=−5.7 (lines 588-589); AT_H=2200, AT_S=6.9
    //     (lines 128-129); saltCorrectS = 0.368·ln((mv+120·√max(0,dv−dntp))/1000) per stack
    //     (lines 623-624, 1042); RC = R·ln(dna_conc/1e9) when both strands symmetric else
    //     R·ln(dna_conc/4e9) (lines 590-593); symmetry_thermo = reverse-complement palindrome
    //     (line 2771). dna_conc is in nM, so /1e9 → mol/L with x=1 and /4e9 → x=4.
    //   Cross-checked against primer3-py 2.3.0 calc_homodimer / calc_heterodimer
    //     (mv=50, dv=0, dntp=0, dna_conc=50 nM): this method reproduces ntthal's ΔH°, ΔS°
    //     and Tm to machine precision for every case whose optimal structure is a contiguous
    //     Watson-Crick duplex (e.g. GCGCGCGC, ACGTACGTACGT, ATCGATCGATCG/CGATCGATCGAT,
    //     CGATCGATCG self-dimer, GCATGC, GGGGCCCC). By design this record keeps the contiguous-WC
    //     optimum (DNA_Dimer_Tm.md §5.3/§5.4): ntthal's mismatch / loop / bulge / terminal-overhang
    //     (tstack2 / dangle) terms are in the full port NtthalDimer — CalculateDimerThermodynamicsNtthal,
    //     CalculateSelfDimerMeltingTemperature and CalculateDimerMeltingTemperature use it (bit-exact to
    //     primer3-py 2.3.1, incl. dv/dntp, temp_c, max_loop; PRIMER-DIMER-001).

    /// <summary>
    /// The most stable intermolecular DNA duplex (self- or hetero-dimer) found between two
    /// oligonucleotides: the aligned spans on each strand, the number of base pairs, and the
    /// dimer ΔH° (kcal/mol), ΔS° (cal/(K·mol)) and ΔG°37 (kcal/mol).
    /// </summary>
    /// <param name="Strand1Start">5'-most aligned index (0-based) on strand 1.</param>
    /// <param name="Strand2Start">5'-most aligned index (0-based) on strand 2.</param>
    /// <param name="BasePairs">Number of contiguous Watson-Crick base pairs in the duplex.</param>
    /// <param name="DeltaH">Dimer ΔH° in kcal/mol (salt-independent).</param>
    /// <param name="DeltaS">Dimer ΔS° in cal/(K·mol), including the 0.368 salt correction.</param>
    /// <param name="DeltaG37">Dimer ΔG°37 = ΔH° − 310.15·ΔS°/1000 in kcal/mol (negative = stable).</param>
    public readonly record struct DimerResult(
        int Strand1Start, int Strand2Start, int BasePairs,
        double DeltaH, double DeltaS, double DeltaG37);

    /// <summary>
    /// Finds the most stable (highest-Tm) gapless, contiguous Watson–Crick intermolecular DNA duplex between two
    /// oligonucleotides, scored with the <c>ntthal</c> duplex terms over the SantaLucia &amp; Hicks (2004) unified
    /// nearest-neighbour model (the full <c>ntthal</c> alignment with mismatches, loops, bulges and terminal
    /// overhangs is <see cref="CalculateDimerThermodynamicsNtthal(string, string, NtthalAlignmentMode, double, double, double, double)"/>,
    /// and the <see cref="FindMostStableDimer(string, string, StructureModel, double?, double?, double?, double?, double?, int?)"/>
    /// overload with <see cref="StructureModel.Ntthal"/> maps that structure onto this record). A <b>self-dimer</b> is obtained by passing the same
    /// sequence as both strands; a <b>hetero/cross-dimer</b> by passing two different sequences.
    /// <b>Opt-in</b>: the duplex (<see cref="CalculateMeltingTemperatureNN"/>) and hairpin Tm
    /// methods, and the default <see cref="CalculateMeltingTemperature(string)"/>, are unchanged.
    /// <para>
    /// The two strands are aligned antiparallel (strand 2 read 3'→5' under strand 1 5'→3') over
    /// every gapless offset; each maximal contiguous Watson-Crick run of ≥ 2 bp is scored as a
    /// duplex with ΔH° = init + Σ stacks + terminal-A·T penalty per A·T-closed end, and ΔS° the
    /// same plus the 0.368·N<sub>stacks</sub>·ln[Na⁺] salt correction. The candidate with the
    /// highest bimolecular Tm is returned (the ΔS° in the result already includes the salt term
    /// for the supplied <paramref name="sodiumMolar"/>).
    /// </para>
    /// </summary>
    /// <param name="strand1">First DNA oligo (5'→3'). Must contain ≥ 2 ACGT bases.</param>
    /// <param name="strand2">Second DNA oligo (5'→3'); the same string as
    /// <paramref name="strand1"/> for a self-dimer. Must contain ≥ 2 ACGT bases.</param>
    /// <param name="sodiumMolar">Monovalent cation concentration in mol/L (default 50 mM); only
    /// the entropy salt correction depends on it.</param>
    /// <param name="strandConcentrationMolar">Total strand concentration C_T in mol/L for the
    /// bimolecular Tm used to rank candidates (default 50 nM, the Primer3/ntthal convention).</param>
    /// <returns>The most stable dimer, or <c>null</c> if either strand is null/&lt; 2 bases/contains
    /// a non-ACGT character, or no duplex of ≥ 2 contiguous base pairs exists between them.</returns>
    public static DimerResult? FindMostStableDimer(
        string strand1,
        string strand2,
        double sodiumMolar = ThermoConstants.DefaultNaConcentration,
        double strandConcentrationMolar = DefaultDimerStrandConcentrationMolar)
    {
        if (string.IsNullOrEmpty(strand1) || string.IsNullOrEmpty(strand2))
            return null;

        string s1 = strand1.ToUpperInvariant();
        string s2 = strand2.ToUpperInvariant();
        if (s1.Length < 2 || s2.Length < 2)
            return null;
        foreach (char c in s1)
            if (c is not ('A' or 'C' or 'G' or 'T')) return null;
        foreach (char c in s2)
            if (c is not ('A' or 'C' or 'G' or 'T')) return null;

        // x = 1 only when both strands are reverse-complement palindromes (ntthal symmetry_thermo).
        bool symmetric = OligoSymmetry.IsSymmetric(s1) && OligoSymmetry.IsSymmetric(s2);
        double x = symmetric ? SelfComplementaryFactor : NonSelfComplementaryFactor;
        // Strand-concentration term R·ln(C_T / x) for the bimolecular Tm (constant over candidates).
        double rcTerm = GasConstant * Math.Log(strandConcentrationMolar / x);
        double saltPerStack = SantaLuciaEntropySaltCoefficient * Math.Log(sodiumMolar);

        // Strand 2 read 3'→5' so its index i pairs base-for-base under strand 1 read 5'→3'.
        string s2Rev = OligoSymmetry.Reverse(s2);
        int n = s1.Length, m = s2.Length;

        DimerResult? best = null;
        double bestTm = double.NegativeInfinity;

        // Slide strand 2 across strand 1 over every gapless antiparallel offset.
        for (int offset = -(m - 1); offset < n; offset++)
        {
            // Within this offset, walk the overlap and split it into maximal contiguous
            // Watson-Crick runs; each run of ≥ 2 bp is a candidate duplex.
            int runStart = -1; // strand-1 index where the current WC run started
            for (int i = 0; i <= n; i++)
            {
                int j = i - offset; // index into s2Rev paired with s1[i]
                bool paired = i < n && j >= 0 && j < m && IsWatsonCrickPair(s1[i], s2Rev[j]);

                if (paired && runStart < 0)
                    runStart = i;

                if (!paired && runStart >= 0)
                {
                    int runEnd = i - 1; // inclusive strand-1 index of the run end
                    EvaluateRun(s1, s2Rev, offset, runStart, runEnd, saltPerStack, rcTerm,
                                ref best, ref bestTm);
                    runStart = -1;
                }
            }
        }

        return best;

        // Scores one contiguous Watson-Crick run [runStart..runEnd] on strand 1 and keeps it
        // if its bimolecular Tm exceeds the best found so far.
        static void EvaluateRun(
            string s1, string s2Rev, int offset, int runStart, int runEnd,
            double saltPerStack, double rcTerm, ref DimerResult? best, ref double bestTm)
        {
            int basePairs = runEnd - runStart + 1;
            if (basePairs < MinDimerBasePairs)
                return;

            double dH = DesignNnInitiation.Initiation.DeltaH;
            double dS = DesignNnInitiation.Initiation.DeltaS;
            for (int k = runStart; k < runEnd; k++)
            {
                string step = s1.Substring(k, 2);
                // A contiguous WC run pairs every column; the stack is the perfect-match NN
                // keyed by the strand-1 dinucleotide (its complement is the strand-2 stack).
                if (!TryGetDesignStack(step, out var p))
                    return;
                dH += p.DeltaH;
                dS += p.DeltaS;
            }

            // Terminal A·T penalty per duplex end closing with an A·T pair.
            if (s1[runStart] is 'A' or 'T') { dH += DesignNnInitiation.TerminalAT.DeltaH; dS += DesignNnInitiation.TerminalAT.DeltaS; }
            if (s1[runEnd] is 'A' or 'T') { dH += DesignNnInitiation.TerminalAT.DeltaH; dS += DesignNnInitiation.TerminalAT.DeltaS; }

            int stacks = basePairs - 1;
            dS += stacks * saltPerStack; // 0.368·N_stacks·ln[Na⁺] (ntthal saltCorrectS)

            double tmKelvin = (dH * 1000.0) / (dS + rcTerm);
            double tmCelsius = tmKelvin - KelvinOffset;

            if (tmCelsius > bestTm)
            {
                bestTm = tmCelsius;
                int strand2Start5 = s2Rev.Length - 1 - (runEnd - offset); // → 5'→3' index on strand 2
                double dG37 = dH - ReferenceTemperatureKelvin * dS / 1000.0;
                best = new DimerResult(runStart, strand2Start5, basePairs, dH, dS, dG37);
            }
        }
    }

    // Minimum base pairs for a dimer duplex (at least one NN stack). ntthal requires a paired
    // region; a single base pair has no stacking energy and is not a duplex.
    private const int MinDimerBasePairs = 2;

    /// <summary>
    /// Full <c>ntthal</c> dimer thermodynamics of the most stable intermolecular duplex.
    /// </summary>
    /// <param name="DeltaH">Dimer ΔH° in kcal/mol (salt-independent).</param>
    /// <param name="DeltaS">Dimer ΔS° in cal/(K·mol), including the N·saltCorrection term.</param>
    /// <param name="DeltaG37">Dimer ΔG° = ΔH° − T·ΔS°/1000 in kcal/mol (negative = stable) at the
    /// analysis temperature T (310.15 K = 37 °C unless an overload with <c>temperatureCelsius</c> is used;
    /// ntthal / primer3-py <c>temp_c</c>). ntthal may report a positive ΔG for its optimal structure.</param>
    /// <param name="TmCelsius">Bimolecular melting temperature in °C.</param>
    /// <param name="BasePairs">Number of paired bases in the optimal structure.</param>
    public readonly record struct DimerThermodynamics(
        double DeltaH, double DeltaS, double DeltaG37, double TmCelsius, int BasePairs);

    /// <summary>
    /// Full <c>ntthal</c> intramolecular-hairpin thermodynamics of the most stable hairpin
    /// (reproduces primer3-py <c>calc_hairpin</c>; special tri/tetraloop bonuses applied).
    /// </summary>
    /// <param name="DeltaH">Hairpin ΔH° in kcal/mol (salt-independent).</param>
    /// <param name="DeltaS">Hairpin ΔS° in cal/(K·mol), including the (N/2−1)·saltCorrection term.</param>
    /// <param name="DeltaG37">Hairpin ΔG = ΔH° − T·ΔS°/1000 in kcal/mol (negative = stable) at T = 310.15 K
    /// (37 °C), or at <c>temperatureCelsius</c> + 273.15 for the overloads taking it (primer3-py <c>temp_c</c>).</param>
    /// <param name="TmCelsius">Unimolecular melting temperature in °C (no strand-concentration term).</param>
    /// <param name="BasePairs">ntthal N/2: half the number of paired positions among bases 1..len−1 of the
    /// optimal structure (the count thal.c uses in the (N/2 − 1)·saltCorrection term).</param>
    public readonly record struct HairpinThermodynamics(
        double DeltaH, double DeltaS, double DeltaG37, double TmCelsius, int BasePairs);

    // Non-empty and A/C/G/T only (case-insensitive): the canonical Core predicate SequenceExtensions.IsValidDna.
    [System.ComponentModel.EditorBrowsable(System.ComponentModel.EditorBrowsableState.Never)] // shared with the separately licensed Seqeron.Genomics.OligoThermodynamics assembly; not part of the API
    public static bool IsAcgtOnly(string? s) => !string.IsNullOrEmpty(s) && s.AsSpan().IsValidDna();

    /// <summary>
    /// Primer3's default (thermodynamic) limit, in °C, on every secondary-structure Tm:
    /// <c>PRIMER_MAX_SELF_ANY_TH</c>, <c>PRIMER_MAX_SELF_END_TH</c>, <c>PRIMER_MAX_HAIRPIN_TH</c>,
    /// <c>PRIMER_PAIR_MAX_COMPL_ANY_TH</c>, <c>PRIMER_PAIR_MAX_COMPL_END_TH</c> = 47.0
    /// (<c>libprimer3.cc</c> <c>pr_set_default_global_args_2</c>). A value is a violation when it is
    /// strictly greater than the limit.
    /// </summary>
    public const double Primer3MaxStructureTm = 47.0;

    /// <summary>Primer3 primer default PRIMER_SALT_MONOVALENT = 50 mM (<c>libprimer3.c</c> <c>p_args.salt_conc</c>).</summary>
    public const double Primer3MonovalentMillimolar = 50.0;

    /// <summary>Primer3 primer default PRIMER_SALT_DIVALENT = 1.5 mM (<c>pr_set_default_global_args_2</c> <c>p_args.divalent_conc</c>).</summary>
    public const double Primer3DivalentMillimolar = 1.5;

    /// <summary>Primer3 primer default PRIMER_DNTP_CONC = 0.6 mM (<c>pr_set_default_global_args_2</c> <c>p_args.dntp_conc</c>).</summary>
    public const double Primer3DntpMillimolar = 0.6;

    /// <summary>Primer3 primer default PRIMER_DNA_CONC = 50 nM (<c>p_args.dna_conc</c>).</summary>
    public const double Primer3DnaConcentrationNanomolar = 50.0;

    // Primer3 _pr_data_control messages for a GC weight without PRIMER_[INTERNAL_]OPT_GC_PERCENT.
    public const string PrimerGcOptimumUndefinedMessage =
        "Primer GC content is part of objective function while optimum gc_content is not defined (Primer3 _pr_data_control).";
    public const string ProbeGcOptimumUndefinedMessage =
        "Hyb probe GC content is part of objective function while optimum gc_content is not defined (Primer3 _pr_data_control).";

    /// <summary>
    /// The Primer3 manual's documented PRIMER_OPT_GC_PERCENT (50 %). It is <b>not</b> applied by default: Primer3's code keeps
    /// the optimum undefined (<c>DEFAULT_OPT_GC_PERCENT</c> = <c>PR_UNDEFINED_INT_OPT</c>, <c>pr_set_default_global_args_1</c>)
    /// and its <c>_pr_data_control</c> rejects a non-zero PRIMER_(INTERNAL_)WT_GC_PERCENT_GT/_LT without an explicit optimum,
    /// as <see cref="PrimerParameters.OptimalGcPercent"/> / <see cref="ProbeDesigner.Primer3ProbeSettings.OptGcPercent"/> /
    /// <see cref="Primer3Optima.OptGcPercent"/> = null now do.
    /// </summary>
    [Obsolete("Primer3 leaves PRIMER_OPT_GC_PERCENT undefined by default; set PrimerParameters.OptimalGcPercent explicitly.")]
    public const double Primer3DefaultOptGcPercent = 50.0;

    /// <summary>
    /// Primer3 default PRIMER_MAX_END_STABILITY = 100 kcal/mol (<c>pr_set_default_global_args_1</c>
    /// <c>max_end_stability</c>): a left/right primer fails when its <c>end_stability</c> =
    /// <c>end_oligodg(seq, 5)</c> = −<see cref="Calculate3PrimeStability"/> is strictly greater. The largest
    /// value an ACGT 3′ pentamer reaches is 6.86 (GCGCG), so the default never rejects a primer.
    /// </summary>
    public const double Primer3MaxEndStability = 100.0;

    /// <summary>Primer3 default PRIMER_GC_CLAMP = 0 (no G/C required at the 3′ end; <c>gc_clamp</c>).</summary>
    public const int Primer3GcClamp = 0;

    /// <summary>
    /// Primer3 default PRIMER_MAX_END_GC = 5 (<c>max_end_gc</c>): the maximum number of G/C among the five 3′-most
    /// bases of a left/right primer; Primer3 only checks it when it is &lt; 5. Legal range 0–5.
    /// </summary>
    public const int Primer3MaxEndGc = 5;

    // Primer3 data control of one oligo-condition set (_pr_data_control: "Illegal value for primer salt or dna
    // concentration" / "… divalent salt or dNTP concentration"): salt and DNA concentration > 0, divalent ≥ 0;
    // a negative dNTP is also rejected here (seqtm has no meaning for it). NaN / ∞ are rejected.
    [System.ComponentModel.EditorBrowsable(System.ComponentModel.EditorBrowsableState.Never)] // shared with the separately licensed Seqeron.Genomics.OligoThermodynamics assembly; not part of the API
    public static void ValidatePrimer3Conditions(
        double monovalentMillimolar, double divalentMillimolar, double dntpMillimolar, double dnaConcentrationNanomolar,
        string paramName)
    {
        if (!(monovalentMillimolar > 0) || double.IsInfinity(monovalentMillimolar)
            || !(dnaConcentrationNanomolar > 0) || double.IsInfinity(dnaConcentrationNanomolar))
            throw new ArgumentOutOfRangeException(paramName,
                "Illegal value for salt or DNA concentration (Primer3: PRIMER_SALT_MONOVALENT and PRIMER_DNA_CONC must be > 0).");
        if (!(divalentMillimolar >= 0) || double.IsInfinity(divalentMillimolar)
            || !(dntpMillimolar >= 0) || double.IsInfinity(dntpMillimolar))
            throw new ArgumentOutOfRangeException(paramName,
                "Illegal value for divalent salt or dNTP concentration (Primer3: PRIMER_SALT_DIVALENT and PRIMER_DNTP_CONC must be ≥ 0).");
    }

    // oligotm.c divalent_to_monovalent added to the monovalent salt: [Mon] + 120·√([Mg²⁺] − [dNTP]) (mM); no divalent ⇒
    // dNTP ignored, Mg ≤ dNTP ⇒ no contribution.
    [System.ComponentModel.EditorBrowsable(System.ComponentModel.EditorBrowsableState.Never)] // shared with the separately licensed Seqeron.Genomics.OligoThermodynamics assembly; not part of the API
    public static double Primer3MonovalentEquivalent(double monovalentMillimolar, double divalentMillimolar, double dntpMillimolar)
    {
        double freeDivalent = divalentMillimolar == 0 ? 0 : Math.Max(0, divalentMillimolar - dntpMillimolar);
        return monovalentMillimolar + Primer3DivalentFactor * Math.Sqrt(freeDivalent);
    }

    /// <summary>Primer3 hybridization-probe (internal-oligo) default PRIMER_INTERNAL_DNA_CONC = 50 nM (<c>libprimer3.cc</c> <c>o_args.dna_conc</c>).</summary>
    public const double Primer3InternalDnaConcentrationNanomolar = 50.0;

    /// <summary>Primer3 hybridization-probe default PRIMER_INTERNAL_SALT_MONOVALENT = 50 mM (<c>o_args.salt_conc</c>).</summary>
    public const double Primer3InternalMonovalentMillimolar = 50.0;

    /// <summary>Primer3 hybridization-probe default PRIMER_INTERNAL_SALT_DIVALENT = 0 mM (<c>o_args.divalent_conc</c>).</summary>
    public const double Primer3InternalDivalentMillimolar = 0.0;

    /// <summary>Primer3 hybridization-probe default PRIMER_INTERNAL_DNTP_CONC = 0 mM (<c>o_args.dntp_conc</c>).</summary>
    public const double Primer3InternalDntpMillimolar = 0.0;

    /// <summary>
    /// Primer3 thermodynamic secondary-structure values of one primer (Tm in °C; 0 when ntthal finds
    /// no structure or the Tm is below 0 °C, as <c>align_thermod</c> reports): <c>PRIMER_*_SELF_ANY_TH</c>, <c>PRIMER_*_SELF_END_TH</c>, <c>PRIMER_*_HAIRPIN_TH</c>.
    /// </summary>
    /// <param name="SelfAnyTh">Self-dimer Tm, ntthal ANY of (primer, primer).</param>
    /// <param name="SelfEndTh">3′-anchored self-dimer Tm, ntthal END1 of (primer, primer).</param>
    /// <param name="HairpinTh">Hairpin Tm, ntthal HAIRPIN of the primer.</param>
    public readonly record struct Primer3OligoStructure(double SelfAnyTh, double SelfEndTh, double HairpinTh)
    {
        /// <summary>True when any value exceeds <paramref name="maxTm"/> (Primer3 rejects the primer).</summary>
        public bool Exceeds(double maxTm = Primer3MaxStructureTm) =>
            SelfAnyTh > maxTm || SelfEndTh > maxTm || HairpinTh > maxTm;
    }

    /// <summary>
    /// Primer3 thermodynamic pair complementarity (Tm in °C; 0 when no structure or Tm &lt; 0 °C):
    /// <c>PRIMER_PAIR_COMPL_ANY_TH</c> and <c>PRIMER_PAIR_COMPL_END_TH</c>.
    /// </summary>
    /// <param name="ComplAnyTh">Hetero-dimer Tm, ntthal ANY of (left, right).</param>
    /// <param name="ComplEndTh">Max of ntthal END1/END2 of (left, right) and of (rc(right), rc(left)).</param>
    public readonly record struct Primer3PairComplementarity(double ComplAnyTh, double ComplEndTh)
    {
        /// <summary>True when either value exceeds <paramref name="maxTm"/> (Primer3 rejects the pair).</summary>
        public bool Exceeds(double maxTm = Primer3MaxStructureTm) => ComplAnyTh > maxTm || ComplEndTh > maxTm;
    }

    // Watson-Crick complement (same left-to-right order) via the canonical Core per-base complement
    // (SequenceExtensions.TryGetComplement); identical to the former local ACGT switch on the
    // validated ACGT duplexes it is applied to.
    private static string Complement(string seq) =>
        string.Create(seq.Length, seq, static (dest, src) => src.AsSpan().TryGetComplement(dest));

    /// <summary>
    /// Primer3's fixed <c>weights.temp_cutoff</c> (= 5 °C, <c>pr_set_default_global_args_1</c>; not a
    /// user-settable tag). In thermodynamic mode a secondary-structure Tm within this many degrees of
    /// the primer Tm is penalised linearly, otherwise by a reciprocal term.
    /// </summary>
    public const double Primer3TempCutoff = 5.0;

    // p_obj_fn thermodynamic secondary-structure term (libprimer3.cc):
    //   if (Tm − temp_cutoff) ≤ s : w · (s − (Tm − temp_cutoff − 1))
    //   else                       : w · 1 / (Tm − temp_cutoff + 1 − s)
    [System.ComponentModel.EditorBrowsable(System.ComponentModel.EditorBrowsableState.Never)] // shared with the separately licensed Seqeron.Genomics.OligoThermodynamics assembly; not part of the API
    public static double ThermodynamicStructurePenalty(double weight, double primerTm, double structureTm)
    {
        if (weight == 0)
            return 0;
        double threshold = primerTm - Primer3TempCutoff;
        return threshold <= structureTm
            ? weight * (structureTm - (threshold - 1.0))
            : weight * (1.0 / (threshold + 1.0 - structureTm));
    }

    private static bool AreComplementary(string s1, string s2)
    {
        if (s1.Length != s2.Length) return false;
        for (int i = 0; i < s1.Length; i++)
        {
            if (!IsComplementary(s1[i], s2[i]))
                return false;
        }
        return true;
    }

    private static bool IsComplementary(char c1, char c2) =>
        (c1 == 'A' && c2 == 'T') || (c1 == 'T' && c2 == 'A') ||
        (c1 == 'G' && c2 == 'C') || (c1 == 'C' && c2 == 'G');
}

/// <summary>
/// Secondary-structure screen applied by <see cref="PrimerDesigner.EvaluatePrimer"/> and
/// <see cref="PrimerDesigner.DesignPrimers"/>.
/// </summary>
public enum PrimerStructureScreen
{
    /// <summary>
    /// Primer3's default thermodynamic screen (<c>PRIMER_THERMODYNAMIC_OLIGO_ALIGNMENT=1</c>): a
    /// primer is rejected when its ntthal self-dimer, 3′ self-dimer or hairpin Tm exceeds 47 °C
    /// (<see cref="PrimerDesigner.CalculatePrimer3OligoStructure"/>), a pair when its ntthal
    /// hetero-dimer or 3′ hetero-dimer Tm exceeds 47 °C
    /// (<see cref="PrimerDesigner.CalculatePrimer3PairComplementarity"/>).
    /// </summary>
    Primer3Thermodynamic = 0,

    /// <summary>
    /// Sequence-only library screen: <see cref="PrimerDesigner.HasHairpinPotential"/> (a ≥ 4-bp
    /// Watson–Crick stem closing a ≥ 3-nt loop; unsourced library heuristic, not Primer3's hairpin screen) per
    /// primer and <see cref="PrimerDesigner.HasPrimerDimer"/>
    /// (Primer3 alignment-mode pair 3′ complementarity ≥ 4) per pair.
    /// </summary>
    Heuristic = 1,

    /// <summary>
    /// Primer3's alignment-score screen (<c>PRIMER_THERMODYNAMIC_OLIGO_ALIGNMENT=0</c>, <c>libprimer3.cc</c>
    /// <c>oligo_compl</c> / <c>characterize_pair</c>): a primer is rejected when its dpal self-complementarity
    /// <c>self_any</c> (<see cref="PrimerDesigner.CalculatePrimerSelfAnyComplementarity"/>) exceeds
    /// <see cref="PrimerParameters.MaxSelfAny"/> (PRIMER_MAX_SELF_ANY, default 8.00) or its 3′ self-complementarity
    /// <c>self_end</c> (<see cref="PrimerDesigner.CalculatePrimerSelfEndComplementarity"/>) exceeds
    /// <see cref="PrimerParameters.MaxSelfEnd"/> (PRIMER_MAX_SELF_END, 3.00); a pair when <c>compl_any</c>
    /// (<see cref="PrimerDesigner.CalculatePrimerDimerAnyComplementarity"/>) exceeds
    /// <see cref="PrimerPairOptions.MaxComplAny"/> (PRIMER_PAIR_MAX_COMPL_ANY, 8.00) or <c>align(left, rc(right),
    /// GLOBAL_END)</c> exceeds <see cref="PrimerPairOptions.MaxComplEnd"/> (PRIMER_PAIR_MAX_COMPL_END, 3.00) or the
    /// reverse orientation <c>align(right, rc(left), GLOBAL_END)</c>, when larger, exceeds PRIMER_MAX_SELF_END
    /// (Primer3 compares that one with the per-primer limit). No hairpin value exists in this mode. An internal
    /// oligo (<see cref="PrimerPairOptions.PickInternalOligo"/>) is screened with PRIMER_INTERNAL_MAX_SELF_ANY /
    /// _SELF_END (12.00). Verified against primer3-py 2.3.1 <c>design_primers</c> with
    /// PRIMER_THERMODYNAMIC_OLIGO_ALIGNMENT = 0.
    /// </summary>
    Primer3Alignment = 2,
}

/// <summary>
/// A primer candidate with quality metrics. <see cref="MeltingTemperature"/> is the Primer3-default
/// Tm rounded to 0.1 °C; <see cref="Score"/> is an informational heuristic (higher is better);
/// <see cref="Penalty"/> is the unrounded Primer3 per-primer penalty (lower is better) that
/// <see cref="PrimerDesigner.DesignPrimers"/> ranks by. Under the default
/// <see cref="PrimerStructureScreen.Primer3Thermodynamic"/> screen <see cref="SelfAnyTh"/>,
/// <see cref="SelfEndTh"/> and <see cref="HairpinTh"/> carry Primer3's PRIMER_*_SELF_ANY_TH /
/// _SELF_END_TH / _HAIRPIN_TH (°C) and <see cref="HasHairpin"/> means HairpinTh &gt; 47 °C; under the
/// heuristic screen they are <c>null</c> and <see cref="HasHairpin"/> is
/// <see cref="PrimerDesigner.HasHairpinPotential"/>.
/// </summary>
public sealed record PrimerCandidate(
    string Sequence,
    int Position,
    bool IsForward,
    int Length,
    double GcContent,
    double MeltingTemperature,
    int HomopolymerLength,
    bool HasHairpin,
    double Stability3Prime,
    bool IsValid,
    IReadOnlyList<string> Issues,
    double Score,
    double Penalty = 0.0,
    double? SelfAnyTh = null,
    double? SelfEndTh = null,
    double? HairpinTh = null)
{
    /// <summary>
    /// Primer3 alignment-mode PRIMER_*_SELF_ANY (<see cref="PrimerDesigner.CalculatePrimerSelfAnyComplementarity"/>);
    /// set under <see cref="PrimerStructureScreen.Primer3Alignment"/>, otherwise <c>null</c>.
    /// </summary>
    public double? SelfAny { get; init; }

    /// <summary>
    /// Primer3 alignment-mode PRIMER_*_SELF_END (<see cref="PrimerDesigner.CalculatePrimerSelfEndComplementarity"/>);
    /// set under <see cref="PrimerStructureScreen.Primer3Alignment"/>, otherwise <c>null</c>.
    /// </summary>
    public double? SelfEnd { get; init; }

    /// <summary>
    /// Primer3 PRIMER_LEFT/RIGHT_n_LIBRARY_MISPRIMING score (<see cref="PrimerDesigner.CalculateLibraryMispriming"/>) when
    /// <see cref="PrimerParameters.MisprimingLibrary"/> is set, otherwise <c>null</c>.
    /// </summary>
    public double? LibraryMispriming { get; init; }

    /// <summary>The library entry named in PRIMER_LEFT/RIGHT_n_LIBRARY_MISPRIMING (with <see cref="LibraryMispriming"/>).</summary>
    public string? LibraryMisprimingName { get; init; }

    /// <summary>
    /// Primer3 PRIMER_LEFT/RIGHT_n_TEMPLATE_MISPRIMING (alignment mode) or _TEMPLATE_MISPRIMING_TH (thermodynamic mode,
    /// °C) — <see cref="TemplateMisprimingScore.Max"/> of <see cref="PrimerDesigner.CalculateTemplateMispriming"/> —
    /// for a primer of a pair returned by <see cref="PrimerDesigner.DesignPrimers"/> / <see cref="PrimerDesigner.DesignPrimerPairs"/>
    /// when Primer3 computes it (a template-mispriming limit or weight is set); otherwise <c>null</c>.
    /// </summary>
    public double? TemplateMispriming { get; init; }

    /// <summary>
    /// Primer3 PRIMER_LEFT/RIGHT_n_BOUND: the fraction (%) bound at <see cref="PrimerParameters.AnnealingTemperature"/>
    /// (<see cref="PrimerDesigner.CalculateFractionBoundPrimer3"/>) when that is &gt; 0 and the value is defined
    /// (≤ 36 ACGT bases); otherwise <c>null</c>.
    /// </summary>
    public double? Bound { get; init; }

    /// <summary>
    /// Primer3 PRIMER_LEFT/RIGHT_n_POSITION_PENALTY (<see cref="PrimerDesigner.CalculatePositionPenaltyPrimer3"/>) of a
    /// designed primer when <see cref="PrimerPairOptions.InsidePenalty"/> / <see cref="PrimerPairOptions.OutsidePenalty"/>
    /// are not Primer3's defaults (−1 / 0); otherwise <c>null</c>.
    /// </summary>
    public double? PositionPenalty { get; init; }

    /// <summary>
    /// Primer3 PRIMER_LEFT/RIGHT_n_MIN_SEQ_QUALITY: the primer's minimum base quality
    /// (<see cref="PrimerDesigner.CalculateSequenceQualityPrimer3"/>) when SEQUENCE_QUALITY
    /// (<see cref="PrimerPairOptions.SequenceQuality"/>) is given; otherwise <c>null</c>.
    /// </summary>
    public int? MinSequenceQuality { get; init; }

    /// <summary>
    /// The primer's predicted PCR failure rate (<see cref="PrimerDesigner.CalculateMaskFailureRatePrimer3"/>, the value Primer3
    /// weights with PRIMER_WT_MASK_FAILURE_RATE) when <see cref="PrimerPairOptions.MaskTemplate"/> is set; otherwise <c>null</c>.
    /// </summary>
    public double? MaskFailureRate { get; init; }
}

/// <summary>
/// Result of primer pair design. For a valid pair the optional values are Primer3's
/// PRIMER_PAIR_k_PENALTY (<see cref="PairPenalty"/>), PRIMER_PAIR_k_PRODUCT_TM (<see cref="ProductTm"/>),
/// PRIMER_PAIR_k_COMPL_ANY_TH / _COMPL_END_TH (<see cref="ComplAnyTh"/>/<see cref="ComplEndTh"/>, thermodynamic
/// screen only) and PRIMER_INTERNAL_k_* (<see cref="InternalOligo"/>, when an internal oligo was requested).
/// </summary>
/// <param name="Forward">Forward (left) primer.</param>
/// <param name="Reverse">Reverse (right) primer.</param>
/// <param name="IsValid">True when the pair satisfies every constraint.</param>
/// <param name="Message">Human-readable status.</param>
/// <param name="ProductSize">Product size (bp).</param>
/// <param name="PairPenalty">Primer3 pair penalty (<c>obj_fn</c>) of a valid pair.</param>
/// <param name="ProductTm">Primer3 product Tm (<c>long_seq_tm</c>, °C) of a valid pair.</param>
/// <param name="ComplAnyTh">Pair hetero-dimer Tm (ntthal ANY, °C).</param>
/// <param name="ComplEndTh">Pair 3′ hetero-dimer Tm (ntthal END1/END2, °C).</param>
/// <param name="InternalOligo">Internal hybridization oligo picked for the pair (PRIMER_PICK_INTERNAL_OLIGO).</param>
public sealed record PrimerPairResult(
    PrimerCandidate? Forward,
    PrimerCandidate? Reverse,
    bool IsValid,
    string Message,
    int ProductSize,
    double? PairPenalty = null,
    double? ProductTm = null,
    double? ComplAnyTh = null,
    double? ComplEndTh = null,
    ProbeDesigner.Primer3Probe? InternalOligo = null)
{
    /// <summary>
    /// Primer3 alignment-mode PRIMER_PAIR_k_COMPL_ANY (<see cref="PrimerDesigner.CalculatePrimerDimerAnyComplementarity"/>)
    /// of a valid pair under <see cref="PrimerStructureScreen.Primer3Alignment"/>; otherwise <c>null</c>.
    /// </summary>
    public double? ComplAny { get; init; }

    /// <summary>
    /// Primer3 alignment-mode PRIMER_PAIR_k_COMPL_END (<see cref="PrimerDesigner.CalculatePrimerDimerEndComplementarity"/>)
    /// of a valid pair under <see cref="PrimerStructureScreen.Primer3Alignment"/>; otherwise <c>null</c>.
    /// </summary>
    public double? ComplEnd { get; init; }

    /// <summary>
    /// Primer3 PRIMER_PAIR_k_LIBRARY_MISPRIMING score of a valid pair (<c>pair_repeat_sim</c>: the maximum over library
    /// entries of the integer part of the left + right primer scores) when <see cref="PrimerParameters.MisprimingLibrary"/>
    /// is set; otherwise <c>null</c>.
    /// </summary>
    public double? LibraryMispriming { get; init; }

    /// <summary>The library entry named in PRIMER_PAIR_k_LIBRARY_MISPRIMING (with <see cref="LibraryMispriming"/>).</summary>
    public string? LibraryMisprimingName { get; init; }

    /// <summary>
    /// Primer3 PRIMER_PAIR_k_TEMPLATE_MISPRIMING (alignment mode) or _TEMPLATE_MISPRIMING_TH (thermodynamic mode, °C) of a
    /// valid pair — max(left same-strand + right other-strand, left other-strand + right same-strand)
    /// (<see cref="TemplateMisprimingScore"/>) — when a pair template-mispriming limit or weight is set; otherwise <c>null</c>.
    /// </summary>
    public double? TemplateMispriming { get; init; }
}

/// <summary>A Primer3 product-size range (PRIMER_PRODUCT_SIZE_RANGE element), inclusive, in bp.</summary>
/// <param name="Min">Smallest product size.</param>
/// <param name="Max">Largest product size.</param>
public readonly record struct ProductSizeRange(int Min, int Max);

/// <summary>
/// Weights of Primer3's pair objective function <c>obj_fn</c> (<c>libprimer3.cc</c>), the
/// <c>PRIMER_PAIR_WT_*</c> tags; defaults are Primer3's (<c>pr_set_default_global_args_1</c>):
/// PRIMER_PAIR_WT_PR_PENALTY = 1, every other weight 0. The thermodynamic complementarity terms use
/// Primer3's fixed <c>temp_cutoff</c> = 5 °C relative to the lower primer Tm (<see cref="PrimerDesigner.Primer3TempCutoff"/>).
/// </summary>
/// <param name="PrimerPenalty">PRIMER_PAIR_WT_PR_PENALTY (× sum of the two primer penalties).</param>
/// <param name="InternalOligoPenalty">PRIMER_PAIR_WT_IO_PENALTY (× internal-oligo penalty).</param>
/// <param name="DiffTm">PRIMER_PAIR_WT_DIFF_TM (× |Tm_left − Tm_right|).</param>
/// <param name="ComplAnyTh">PRIMER_PAIR_WT_COMPL_ANY_TH (thermodynamic screen only).</param>
/// <param name="ComplEndTh">PRIMER_PAIR_WT_COMPL_END_TH (thermodynamic screen only).</param>
/// <param name="ProductTmLt">PRIMER_PAIR_WT_PRODUCT_TM_LT (needs <see cref="PrimerPairOptions.ProductOptTm"/>).</param>
/// <param name="ProductTmGt">PRIMER_PAIR_WT_PRODUCT_TM_GT (needs <see cref="PrimerPairOptions.ProductOptTm"/>).</param>
/// <param name="ProductSizeLt">PRIMER_PAIR_WT_PRODUCT_SIZE_LT (needs <see cref="PrimerPairOptions.ProductOptSize"/>).</param>
/// <param name="ProductSizeGt">PRIMER_PAIR_WT_PRODUCT_SIZE_GT (needs <see cref="PrimerPairOptions.ProductOptSize"/>).</param>
public sealed record Primer3PairWeights(
    double PrimerPenalty = 1.0,
    double InternalOligoPenalty = 0.0,
    double DiffTm = 0.0,
    double ComplAnyTh = 0.0,
    double ComplEndTh = 0.0,
    double ProductTmLt = 0.0,
    double ProductTmGt = 0.0,
    double ProductSizeLt = 0.0,
    double ProductSizeGt = 0.0)
{
    /// <summary>PRIMER_PAIR_WT_COMPL_ANY (× PRIMER_PAIR_COMPL_ANY; <see cref="PrimerStructureScreen.Primer3Alignment"/> only; Primer3 default 0).</summary>
    public double ComplAny { get; init; }

    /// <summary>PRIMER_PAIR_WT_COMPL_END (× PRIMER_PAIR_COMPL_END; <see cref="PrimerStructureScreen.Primer3Alignment"/> only; Primer3 default 0).</summary>
    public double ComplEnd { get; init; }

    /// <summary>
    /// PRIMER_PAIR_WT_LIBRARY_MISPRIMING (× PRIMER_PAIR_LIBRARY_MISPRIMING, <see cref="PrimerPairResult.LibraryMispriming"/>;
    /// Primer3 default 0). Non-zero requires <see cref="PrimerParameters.MisprimingLibrary"/>.
    /// </summary>
    public double LibraryMispriming { get; init; }

    /// <summary>
    /// PRIMER_PAIR_WT_TEMPLATE_MISPRIMING (× PRIMER_PAIR_k_TEMPLATE_MISPRIMING, <see cref="PrimerPairResult.TemplateMispriming"/>;
    /// alignment mode, <see cref="PrimerParameters.ThermodynamicTemplateAlignment"/> false; Primer3 default 0).
    /// </summary>
    public double TemplateMispriming { get; init; }

    /// <summary>
    /// PRIMER_PAIR_WT_TEMPLATE_MISPRIMING_TH (thermodynamic mode; Primer3 default 0): the 5 °C <c>temp_cutoff</c> rule
    /// relative to the lower primer Tm, as for <see cref="ComplAnyTh"/>.
    /// </summary>
    public double TemplateMisprimingTh { get; init; }
}
