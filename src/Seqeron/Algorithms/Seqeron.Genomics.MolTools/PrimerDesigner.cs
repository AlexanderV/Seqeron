using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.MolTools;

/// <summary>
/// Designs PCR primers for DNA sequences with various quality criteria.
/// </summary>
public static class PrimerDesigner
{
    /// <summary>
    /// Default primer design parameters.
    /// </summary>
    public static readonly PrimerParameters DefaultParameters = new(
        MinLength: 18,
        MaxLength: 25,
        OptimalLength: 20,
        MinGcContent: 40,
        MaxGcContent: 60,
        MinTm: 57,
        MaxTm: 63,
        OptimalTm: 60,
        MaxHomopolymer: 4,
        MaxDinucleotideRepeats: 4,
        Avoid3PrimeGC: false,
        Check3PrimeStability: true
    );

    /// <summary>
    /// Designs forward and reverse primers for a target region.
    /// </summary>
    /// <param name="template">The DNA template sequence.</param>
    /// <param name="targetStart">Start position of target region.</param>
    /// <param name="targetEnd">End position of target region.</param>
    /// <param name="parameters">Primer design parameters (optional).</param>
    /// <returns>Primer pair result.</returns>
    public static PrimerPairResult DesignPrimers(
        DnaSequence template,
        int targetStart,
        int targetEnd,
        PrimerParameters? parameters = null)
    {
        var param = parameters ?? DefaultParameters;

        if (targetStart < 0 || targetEnd >= template.Length || targetStart >= targetEnd)
            throw new ArgumentException("Invalid target region.");

        // Design forward primer (upstream of target)
        var forwardCandidates = new List<PrimerCandidate>();
        int forwardSearchStart = Math.Max(0, targetStart - 200);
        int forwardSearchEnd = targetStart;

        for (int start = forwardSearchStart; start < forwardSearchEnd; start++)
        {
            for (int len = param.MinLength; len <= param.MaxLength && start + len <= targetStart; len++)
            {
                var candidate = EvaluatePrimer(template.Sequence.Substring(start, len), start, true, param);
                if (candidate.IsValid)
                    forwardCandidates.Add(candidate);
            }
        }

        // Design reverse primer (downstream of target, on reverse complement)
        var reverseCandidates = new List<PrimerCandidate>();
        int reverseSearchStart = targetEnd;
        int reverseSearchEnd = Math.Min(template.Length, targetEnd + 200);

        for (int end = reverseSearchStart + param.MinLength; end <= reverseSearchEnd; end++)
        {
            for (int len = param.MinLength; len <= param.MaxLength && end - len >= targetEnd; len++)
            {
                int start = end - len;
                var seq = template.Sequence.Substring(start, len);
                var revComp = new DnaSequence(seq).ReverseComplement().Sequence;
                var candidate = EvaluatePrimer(revComp, start, false, param);
                if (candidate.IsValid)
                    reverseCandidates.Add(candidate);
            }
        }

        // Select best pair
        var bestForward = forwardCandidates
            .OrderByDescending(c => c.Score)
            .FirstOrDefault();

        var bestReverse = reverseCandidates
            .OrderByDescending(c => c.Score)
            .FirstOrDefault();

        if (bestForward == null || bestReverse == null)
        {
            return new PrimerPairResult(
                null, null, false,
                "Could not find valid primers for the target region.",
                0
            );
        }

        // Check primer pair compatibility
        double tmDiff = Math.Abs(bestForward.MeltingTemperature - bestReverse.MeltingTemperature);
        bool isCompatible = tmDiff <= 5.0 && !HasPrimerDimer(bestForward.Sequence, bestReverse.Sequence);

        int productSize = bestReverse.Position + bestReverse.Sequence.Length - bestForward.Position;

        return new PrimerPairResult(
            Forward: bestForward,
            Reverse: bestReverse,
            IsValid: isCompatible,
            Message: isCompatible ? "Valid primer pair found." : $"Primer Tm difference: {tmDiff:F1}°C",
            ProductSize: productSize
        );
    }

    /// <summary>
    /// Evaluates a single primer candidate.
    /// </summary>
    public static PrimerCandidate EvaluatePrimer(
        string sequence,
        int position,
        bool isForward,
        PrimerParameters? parameters = null)
    {
        var param = parameters ?? DefaultParameters;
        var seq = sequence.ToUpperInvariant();

        double gcContent = CalculateGcContent(seq);
        double tm = CalculateMeltingTemperature(seq);
        int homopolymer = FindLongestHomopolymer(seq);
        int dinucRepeat = FindLongestDinucleotideRepeat(seq);
        bool hasHairpin = HasHairpinPotential(seq);
        double stability3Prime = Calculate3PrimeStability(seq);

        var issues = new List<string>();

        // Validate against parameters
        if (seq.Length < param.MinLength || seq.Length > param.MaxLength)
            issues.Add($"Length {seq.Length} outside range [{param.MinLength}-{param.MaxLength}]");

        if (gcContent < param.MinGcContent || gcContent > param.MaxGcContent)
            issues.Add($"GC content {gcContent:F1}% outside range [{param.MinGcContent}-{param.MaxGcContent}]%");

        if (tm < param.MinTm || tm > param.MaxTm)
            issues.Add($"Tm {tm:F1}°C outside range [{param.MinTm}-{param.MaxTm}]°C");

        if (homopolymer > param.MaxHomopolymer)
            issues.Add($"Homopolymer run of {homopolymer} exceeds max {param.MaxHomopolymer}");

        if (dinucRepeat > param.MaxDinucleotideRepeats)
            issues.Add($"Dinucleotide repeat of {dinucRepeat} exceeds max {param.MaxDinucleotideRepeats}");

        if (hasHairpin)
            issues.Add("Potential hairpin structure detected");

        if (param.Check3PrimeStability && stability3Prime < -9)
            issues.Add($"3' end too stable (ΔG = {stability3Prime:F1} kcal/mol)");

        // Check 3' end for GC clamp
        if (param.Avoid3PrimeGC && seq.Length >= 2)
        {
            string last2 = seq.Substring(seq.Length - 2);
            int gcCount = last2.Count(c => c == 'G' || c == 'C');
            if (gcCount == 0)
                issues.Add("No GC clamp at 3' end");
        }

        bool isValid = issues.Count == 0;

        // Calculate score
        double score = CalculatePrimerScore(seq, gcContent, tm, homopolymer, param);

        return new PrimerCandidate(
            Sequence: seq,
            Position: position,
            IsForward: isForward,
            Length: seq.Length,
            GcContent: Math.Round(gcContent, 1),
            MeltingTemperature: Math.Round(tm, 1),
            HomopolymerLength: homopolymer,
            HasHairpin: hasHairpin,
            Stability3Prime: Math.Round(stability3Prime, 1),
            IsValid: isValid,
            Issues: issues.AsReadOnly(),
            Score: Math.Round(score, 2)
        );
    }

    /// <summary>
    /// Calculates the melting temperature for DNA primers.
    /// Uses Wallace rule for short primers (&lt; 14 valid bases)
    /// and Marmur-Doty formula for longer primers (≥ 14 valid bases).
    /// Only standard DNA bases (A, C, G, T) are recognized;
    /// all other characters are ignored.
    /// </summary>
    public static double CalculateMeltingTemperature(string primer)
    {
        if (string.IsNullOrEmpty(primer))
            return 0;

        var seq = primer.ToUpperInvariant();

        int at = seq.Count(c => c == 'A' || c == 'T');
        int gc = seq.Count(c => c == 'G' || c == 'C');
        int validLength = at + gc;

        if (validLength == 0)
            return 0;

        // For short primers (< 14 valid bases), use Wallace rule
        if (validLength < ThermoConstants.WallaceMaxLength)
        {
            return ThermoConstants.CalculateWallaceTm(at, gc);
        }

        // For longer primers, use Marmur-Doty formula
        return Math.Max(0, ThermoConstants.CalculateMarmurDotyTm(gc, validLength));
    }

    /// <summary>
    /// Calculates the melting temperature with salt correction.
    /// </summary>
    /// <param name="primer">Primer sequence.</param>
    /// <param name="naConcentration">Na+ concentration in mM (default: 50).</param>
    /// <returns>Corrected melting temperature in °C.</returns>
    public static double CalculateMeltingTemperatureWithSalt(string primer, double naConcentration = 50)
    {
        if (string.IsNullOrEmpty(primer))
            return 0;

        double baseTm = CalculateMeltingTemperature(primer);
        double saltCorrection = ThermoConstants.CalculateSaltCorrection(naConcentration);
        return Math.Round(baseTm + saltCorrection, 1);
    }

    /// <summary>
    /// Calculates GC content as a percentage.
    /// </summary>
    public static double CalculateGcContent(string sequence) =>
        string.IsNullOrEmpty(sequence) ? 0 : sequence.CalculateGcContentFast();

    /// <summary>
    /// Finds the longest homopolymer run (consecutive identical nucleotides).
    /// </summary>
    public static int FindLongestHomopolymer(string sequence)
    {
        if (string.IsNullOrEmpty(sequence))
            return 0;

        int maxRun = 1;
        int currentRun = 1;

        for (int i = 1; i < sequence.Length; i++)
        {
            if (char.ToUpperInvariant(sequence[i]) == char.ToUpperInvariant(sequence[i - 1]))
            {
                currentRun++;
                maxRun = Math.Max(maxRun, currentRun);
            }
            else
            {
                currentRun = 1;
            }
        }

        return maxRun;
    }

    /// <summary>
    /// Finds the longest dinucleotide repeat (e.g., ATATAT).
    /// </summary>
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
    /// Checks if primer has potential to form hairpin structure.
    /// Uses O(n²) algorithm for short sequences, suffix tree O(n) for long sequences.
    /// </summary>
    /// <param name="sequence">DNA sequence to check.</param>
    /// <param name="minStemLength">Minimum stem length (default 4).</param>
    /// <param name="minLoopLength">Minimum loop length (default 3).</param>
    /// <returns>True if hairpin potential detected.</returns>
    public static bool HasHairpinPotential(string sequence, int minStemLength = 4, int minLoopLength = 3)
    {
        if (string.IsNullOrEmpty(sequence) || sequence.Length < minStemLength * 2 + minLoopLength)
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
                if (AreComplementary(fragment, Reverse(target)))
                    return true;
            }
        }

        return false;
    }

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
    /// Checks if two primers can form primer-dimer.
    /// </summary>
    public static bool HasPrimerDimer(string primer1, string primer2, int minComplementarity = 4)
    {
        if (string.IsNullOrEmpty(primer1) || string.IsNullOrEmpty(primer2))
            return false;

        var seq1 = primer1.ToUpperInvariant();
        var seq2 = DnaSequence.GetReverseComplementString(primer2.ToUpperInvariant());

        // Check 3' end complementarity (most problematic for extension)
        int checkLength = Math.Min(8, Math.Min(seq1.Length, seq2.Length));
        string end1 = seq1.Substring(seq1.Length - checkLength);
        string end2 = seq2.Substring(0, checkLength);

        int complementary = 0;
        for (int i = 0; i < checkLength; i++)
        {
            if (IsComplementary(end1[i], end2[i]))
                complementary++;
        }

        return complementary >= minComplementarity;
    }

    /// <summary>
    /// Calculates the stability of the 3' end (last 5 bases) as a duplex ΔG°37.
    /// Uses SantaLucia (1998) unified nearest-neighbor parameters with initiation.
    /// More negative = more stable = potentially problematic.
    /// Matches Primer3 PRIMER_MAX_END_STABILITY calculation.
    /// </summary>
    public static double Calculate3PrimeStability(string sequence)
    {
        if (string.IsNullOrEmpty(sequence) || sequence.Length < 5)
            return 0;

        var seq = sequence.ToUpperInvariant();
        string last5 = seq.Substring(seq.Length - 5);

        // Nearest-neighbor ΔG°37 values in kcal/mol
        // Source: SantaLucia (1998) PNAS 95:1460-65, Table 1, unified parameters (1 M NaCl)
        var deltaG = new Dictionary<string, double>
        {
            ["AA"] = -1.0,
            ["TT"] = -1.0,
            ["AT"] = -0.88,
            ["TA"] = -0.58,
            ["CA"] = -1.45,
            ["TG"] = -1.45,
            ["GT"] = -1.44,
            ["AC"] = -1.44,
            ["CT"] = -1.28,
            ["AG"] = -1.28,
            ["GA"] = -1.30,
            ["TC"] = -1.30,
            ["CG"] = -2.17,
            ["GC"] = -2.24,
            ["GG"] = -1.84,
            ["CC"] = -1.84
        };

        double totalDeltaG = 0;
        for (int i = 0; i < last5.Length - 1; i++)
        {
            string dinuc = last5.Substring(i, 2);
            if (deltaG.TryGetValue(dinuc, out double dg))
                totalDeltaG += dg;
        }

        // Initiation parameters per SantaLucia (1998) Table 1:
        // Init w/terminal G·C: +0.98 kcal/mol
        // Init w/terminal A·T: +1.03 kcal/mol
        // Primer3 PRIMER_MAX_END_STABILITY includes these (GCGCG = -6.86, TATAT = -0.86).
        totalDeltaG += IsGC(last5[0]) ? 0.98 : 1.03;
        totalDeltaG += IsGC(last5[^1]) ? 0.98 : 1.03;

        return totalDeltaG;

        static bool IsGC(char c) => c is 'G' or 'C';
    }

    // ---- Nearest-neighbour salt-corrected Tm (PRIMER-TM-001, opt-in) ----------
    // SantaLucia (1998) unified Watson-Crick NN ΔH°/ΔS° (1 M NaCl) with the
    // bimolecular Tm equation, plus published monovalent (Owczarzy 2004) and
    // divalent (Owczarzy 2008) salt corrections. This is an OPT-IN design Tm: the
    // default CalculateMeltingTemperature (Wallace / Marmur-Doty) is unchanged.
    // Sources:
    //   SantaLucia J (1998) PNAS 95(4):1460-65, Table 1 (unified NN parameters),
    //     Eq. 3 Tm = ΔH°·1000 / (ΔS° + R·ln(C_T/x)) − 273.15;
    //   SantaLucia J, Hicks D (2004) Annu Rev Biophys 33:415-440, Table 1 + Eq. 5
    //     (cross-check of the unified parameters and the entropy salt correction);
    //   Owczarzy R et al. (2004) Biochemistry 43:3537-54 (monovalent correction);
    //   Owczarzy R et al. (2008) Biochemistry 47:5336-53 (divalent Mg²⁺ correction);
    //   Biopython Bio.SeqUtils.MeltingTemp (DNA_NN4 table, salt_correction methods
    //     6 and 7 — reference implementation, cross-checked verbatim).

    /// <summary>
    /// SantaLucia (1998) unified Watson-Crick nearest-neighbour parameters at 1 M NaCl,
    /// as (ΔH° in kcal/mol, ΔS° in cal/(K·mol)). 5'→3' dinucleotide keys; the reverse
    /// strand is implied by Watson-Crick pairing (e.g. AC pairs with the GT NN, hence
    /// AC and GT share parameters). Source: SantaLucia &amp; Hicks (2004) Table 1
    /// (identical to SantaLucia 1998); cross-checked against Biopython DNA_NN4.
    /// </summary>
    private static readonly Dictionary<string, (double DeltaH, double DeltaS)> NnUnifiedParams = new()
    {
        ["AA"] = (-7.6, -21.3), ["TT"] = (-7.6, -21.3),
        ["AT"] = (-7.2, -20.4),
        ["TA"] = (-7.2, -21.3),
        ["CA"] = (-8.5, -22.7), ["TG"] = (-8.5, -22.7),
        ["GT"] = (-8.4, -22.4), ["AC"] = (-8.4, -22.4),
        ["CT"] = (-7.8, -21.0), ["AG"] = (-7.8, -21.0),
        ["GA"] = (-8.2, -22.2), ["TC"] = (-8.2, -22.2),
        ["CG"] = (-10.6, -27.2),
        ["GC"] = (-9.8, -24.4),
        ["GG"] = (-8.0, -19.9), ["CC"] = (-8.0, -19.9)
    };

    // Duplex-initiation term (per duplex). SantaLucia & Hicks (2004) Table 1: ΔH°=+0.2, ΔS°=−5.7.
    private const double NnInitDeltaH = 0.2;    // kcal/mol
    private const double NnInitDeltaS = -5.7;   // cal/(K·mol)

    // Terminal A·T penalty, applied once per duplex end that closes with an A·T pair.
    // SantaLucia & Hicks (2004) Table 1: ΔH°=+2.2, ΔS°=+6.9.
    private const double NnTerminalAtDeltaH = 2.2;   // kcal/mol
    private const double NnTerminalAtDeltaS = 6.9;   // cal/(K·mol)

    // Symmetry correction, applied once for a self-complementary duplex.
    // SantaLucia & Hicks (2004) Table 1: ΔH°=0.0, ΔS°=−1.4.
    private const double NnSymmetryDeltaS = -1.4;    // cal/(K·mol)

    // Gas constant R in cal/(K·mol) for the Tm equation. SantaLucia & Hicks (2004) Eq. 3: R = 1.9872.
    private const double GasConstant = 1.9872;

    // Strand-concentration divisor x in Tm = ΔH°/(ΔS° + R·ln(C_T/x)).
    // SantaLucia & Hicks (2004) Eq. 3: x = 4 for non-self-complementary, x = 1 for self-complementary.
    private const double NonSelfComplementaryFactor = 4.0;
    private const double SelfComplementaryFactor = 1.0;

    // Kelvin-to-Celsius offset.
    private const double KelvinOffset = 273.15;

    // Default total strand concentration C_T = 0.5 µM (a common PCR primer working
    // concentration). Exposed as a parameter; the caller may override.
    private const double DefaultStrandConcentrationMolar = 0.5e-6;

    // Owczarzy (2004) monovalent (Na⁺) correction, 1/Tm form:
    //   1/Tm[Na] = 1/Tm[1M] + (4.29·f(GC) − 3.95)·1e-5·ln[Na⁺] + 9.40e-6·(ln[Na⁺])²
    // Source: Owczarzy et al. (2004) Biochemistry 43:3537-54; coefficients per the
    // Biopython salt_correction method 6 (cross-checked).
    private const double Owczarzy2004GcCoefficient = 4.29e-5;
    private const double Owczarzy2004Constant = 3.95e-5;
    private const double Owczarzy2004QuadraticCoefficient = 9.40e-6;

    /// <summary>Salt-correction mode for <see cref="CalculateMeltingTemperatureNN"/>.</summary>
    public enum SaltCorrectionMode
    {
        /// <summary>No correction — Tm at the SantaLucia 1 M NaCl reference state.</summary>
        None,

        /// <summary>
        /// SantaLucia &amp; Hicks (2004) Eq. 5 entropy correction:
        /// ΔS°[Na] = ΔS°[1 M] + 0.368·(N/2)·ln[Na⁺], N = total phosphates = 2·(length−1).
        /// Fully primary-sourced; applied to ΔS° before the Tm equation.
        /// </summary>
        SantaLuciaEntropy,

        /// <summary>
        /// Owczarzy et al. (2004) monovalent quadratic 1/Tm correction (Biochemistry 43:3537).
        /// </summary>
        Owczarzy2004Monovalent,

        /// <summary>
        /// Owczarzy et al. (2008) divalent Mg²⁺ (and dNTP-adjusted) correction (Biochemistry 47:5336);
        /// reduces to the 2004 monovalent form when the divalent ratio is negligible.
        /// </summary>
        Owczarzy2008Divalent
    }

    // SantaLucia & Hicks (2004) Eq. 5 entropy salt-correction coefficient (0.368).
    private const double SantaLuciaEntropySaltCoefficient = 0.368;

    // Owczarzy (2008) divalent correction coefficients (Biopython salt_correction method 7;
    // Owczarzy et al. 2008 Biochemistry 47:5336). The base coefficients a..g and the
    // regime-dependent reparameterisations of a, d, g below are taken verbatim.
    private const double Owc2008A = 3.92e-5;
    private const double Owc2008B = -0.911e-5;
    private const double Owc2008C = 6.26e-5;
    private const double Owc2008D = 1.42e-5;
    private const double Owc2008E = -48.2e-5;
    private const double Owc2008F = 52.5e-5;
    private const double Owc2008G = 8.31e-5;
    private const double Owc2008MonovalentRatioLow = 0.22;   // R = √[Mg²⁺]/[Mon] threshold below which monovalent dominates
    private const double Owc2008MonovalentRatioHigh = 6.0;   // R threshold above which divalent dominates
    private const double DntpMgAssociationConstant = 3.0e4;  // dNTP·Mg²⁺ association constant Ka

    /// <summary>
    /// Computes the duplex ΔH° (kcal/mol) and ΔS° (cal/(K·mol)) of a DNA oligonucleotide
    /// using the SantaLucia (1998) unified nearest-neighbour parameters, including the
    /// duplex-initiation term, a terminal A·T penalty per A·T-closed end, and (for a
    /// self-complementary sequence) the symmetry correction. Only A/C/G/T are summed;
    /// any other character makes the NN lookup fail and the result is reported as not
    /// computable (returns <c>null</c>). Source: SantaLucia &amp; Hicks (2004) Eq. 1 + Table 1.
    /// </summary>
    /// <param name="sequence">DNA sequence (one strand, 5'→3').</param>
    /// <returns>(ΔH°, ΔS°, IsSelfComplementary) or <c>null</c> if the sequence is empty,
    /// shorter than 2 bases, or contains a non-ACGT character.</returns>
    public static (double DeltaH, double DeltaS, bool IsSelfComplementary)? CalculateNearestNeighborThermodynamics(string sequence)
    {
        if (string.IsNullOrEmpty(sequence))
            return null;

        string seq = sequence.ToUpperInvariant();
        if (seq.Length < 2)
            return null;

        double dH = NnInitDeltaH;
        double dS = NnInitDeltaS;

        for (int i = 0; i < seq.Length - 1; i++)
        {
            string dinuc = seq.Substring(i, 2);
            if (!NnUnifiedParams.TryGetValue(dinuc, out var p))
                return null; // non-ACGT base present
            dH += p.DeltaH;
            dS += p.DeltaS;
        }

        // Terminal A·T penalty per end that closes with an A·T pair (A or T terminus).
        if (seq[0] is 'A' or 'T') { dH += NnTerminalAtDeltaH; dS += NnTerminalAtDeltaS; }
        if (seq[^1] is 'A' or 'T') { dH += NnTerminalAtDeltaH; dS += NnTerminalAtDeltaS; }

        bool selfComp = IsSelfComplementary(seq);
        if (selfComp)
            dS += NnSymmetryDeltaS; // symmetry correction (ΔH° contribution is 0)

        return (dH, dS, selfComp);
    }

    /// <summary>
    /// Computes the design melting temperature (°C) of a primer/oligonucleotide using the
    /// SantaLucia (1998) unified nearest-neighbour thermodynamics and the bimolecular Tm
    /// equation, with an optional published salt correction. <b>Opt-in</b>: the default
    /// <see cref="CalculateMeltingTemperature(string)"/> (Wallace / Marmur-Doty) is unchanged.
    /// <para>
    /// Tm = ΔH°·1000 / (ΔS° + R·ln(C_T / x)) − 273.15, with R = 1.9872 cal/(K·mol),
    /// x = 4 for a non-self-complementary duplex and x = 1 for a self-complementary one
    /// (SantaLucia &amp; Hicks 2004, Eq. 3). Salt corrections per
    /// <paramref name="saltMode"/>.
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
        var thermo = CalculateNearestNeighborThermodynamics(primer);
        if (thermo is null)
            return double.NaN;

        var (dH, dS, selfComp) = thermo.Value;
        int length = primer.Length;
        double x = selfComp ? SelfComplementaryFactor : NonSelfComplementaryFactor;

        // SantaLucia & Hicks (2004) Eq. 5: salt-correct ΔS° before the Tm equation.
        double dSeff = dS;
        if (saltMode == SaltCorrectionMode.SantaLuciaEntropy)
        {
            // N = total phosphates in the duplex = 2·(length − 1) (paper's 6-bp duplex → 10).
            double phosphates = 2.0 * (length - 1);
            dSeff += SantaLuciaEntropySaltCoefficient * (phosphates / 2.0) * Math.Log(sodiumMolar);
        }

        // Tm in Kelvin (ΔH° converted kcal → cal via ·1000).
        double tmKelvin = (dH * 1000.0) / (dSeff + GasConstant * Math.Log(strandConcentrationMolar / x));

        switch (saltMode)
        {
            case SaltCorrectionMode.Owczarzy2004Monovalent:
                tmKelvin = ApplyOwczarzy2004(tmKelvin, primer, sodiumMolar);
                break;
            case SaltCorrectionMode.Owczarzy2008Divalent:
                tmKelvin = ApplyOwczarzy2008(tmKelvin, primer, sodiumMolar, magnesiumMolar, dntpMolar);
                break;
            case SaltCorrectionMode.None:
            case SaltCorrectionMode.SantaLuciaEntropy:
            default:
                break;
        }

        return tmKelvin - KelvinOffset;
    }

    // ---- NN internal-mismatch + dangling-end Tm (PRIMER-TM-001, opt-in extension) -------
    // Extends the perfect-match NN model with published internal single-mismatch and
    // single-dangling-end ΔH°/ΔS° terms so the NN Tm can be computed for a probe–target
    // duplex that contains an internal mismatch and/or an unpaired dangling end. The
    // perfect-match CalculateMeltingTemperatureNN above is UNCHANGED; this is opt-in.
    //
    // Convention (mirrors Biopython Bio.SeqUtils.MeltingTemp.Tm_NN with imm_table=DNA_IMM,
    // de_table=DNA_DE): the top strand is 5'→3'; the bottom strand is supplied 3'→5'
    // (i.e. the complement of the top read in the SAME left-to-right order, NOT the
    // reverse complement), so position i of the bottom is the base paired under position i
    // of the top. A '.' in either strand marks the single unpaired base of a dangling end.
    // A nearest-neighbour key is "topPair/bottomPair" (two bases of each strand, slash-
    // separated); an internal mismatch is looked up in the mismatch table, trying the
    // forward key then its character-reverse, exactly as Tm_NN does.
    //
    // Sources (retrieved & cross-checked this session):
    //   Internal single mismatches — Allawi & SantaLucia (1997) Biochemistry 36:10581
    //     (G·T); Allawi & SantaLucia (1998) Biochemistry 37:9435 (G·A), 37:2170 (C·T),
    //     Nucleic Acids Res 26:2694 (A·C, C·C variants); Peyret et al. (1999) Biochemistry
    //     38:3468 (A·A, C·C, G·G, T·T). Values transcribed verbatim from Biopython DNA_IMM
    //     and cross-checked against the SantaLucia & Hicks (2004) Table 2 worked example
    //     (5'-GGACTGACG-3'/3'-CCTGGCTGC-5' → ΔG°37 ≈ −8.3 kcal/mol).
    //   Single dangling ends — Bommarito, Peyret & SantaLucia (2000) Nucleic Acids Res
    //     28:1929. Values transcribed from Biopython DNA_DE and cross-checked term-by-term
    //     against SantaLucia & Hicks (2004) Table 3 ΔH° (all 32 entries reproduce exactly).

    /// <summary>
    /// Allawi/SantaLucia/Peyret internal single-mismatch nearest-neighbour parameters
    /// (ΔH° kcal/mol, ΔS° cal/(K·mol)) at 1 M NaCl. Key = "topPair/bottomPair" with the
    /// bottom strand written 3'→5' (complement direction); the mismatched base pair is one
    /// of the two columns. Transcribed verbatim from Biopython <c>DNA_IMM</c> (the Watson-
    /// Crick / inosine entries are excluded — only A/C/G/T single mismatches are kept).
    /// Source: Allawi &amp; SantaLucia (1997/1998); Peyret et al. (1999); cross-checked
    /// against SantaLucia &amp; Hicks (2004) Table 2.
    /// </summary>
    private static readonly Dictionary<string, (double DeltaH, double DeltaS)> NnInternalMismatch = new()
    {
        ["AG/TT"] = (1.0, 0.9), ["AT/TG"] = (-2.5, -8.3), ["CG/GT"] = (-4.1, -11.7),
        ["CT/GG"] = (-2.8, -8.0), ["GG/CT"] = (3.3, 10.4), ["GG/TT"] = (5.8, 16.3),
        ["GT/CG"] = (-4.4, -12.3), ["GT/TG"] = (4.1, 9.5), ["TG/AT"] = (-0.1, -1.7),
        ["TG/GT"] = (-1.4, -6.2), ["TT/AG"] = (-1.3, -5.3), ["AA/TG"] = (-0.6, -2.3),
        ["AG/TA"] = (-0.7, -2.3), ["CA/GG"] = (-0.7, -2.3), ["CG/GA"] = (-4.0, -13.2),
        ["GA/CG"] = (-0.6, -1.0), ["GG/CA"] = (0.5, 3.2), ["TA/AG"] = (0.7, 0.7),
        ["TG/AA"] = (3.0, 7.4), ["AC/TT"] = (0.7, 0.2), ["AT/TC"] = (-1.2, -6.2),
        ["CC/GT"] = (-0.8, -4.5), ["CT/GC"] = (-1.5, -6.1), ["GC/CT"] = (2.3, 5.4),
        ["GT/CC"] = (5.2, 13.5), ["TC/AT"] = (1.2, 0.7), ["TT/AC"] = (1.0, 0.7),
        ["AA/TC"] = (2.3, 4.6), ["AC/TA"] = (5.3, 14.6), ["CA/GC"] = (1.9, 3.7),
        ["CC/GA"] = (0.6, -0.6), ["GA/CC"] = (5.2, 14.2), ["GC/CA"] = (-0.7, -3.8),
        ["TA/AC"] = (3.4, 8.0), ["TC/AA"] = (7.6, 20.2), ["AA/TA"] = (1.2, 1.7),
        ["CA/GA"] = (-0.9, -4.2), ["GA/CA"] = (-2.9, -9.8), ["TA/AA"] = (4.7, 12.9),
        ["AC/TC"] = (0.0, -4.4), ["CC/GC"] = (-1.5, -7.2), ["GC/CC"] = (3.6, 8.9),
        ["TC/AC"] = (6.1, 16.4), ["AG/TG"] = (-3.1, -9.5), ["CG/GG"] = (-4.9, -15.3),
        ["GG/CG"] = (-6.0, -15.8), ["TG/AG"] = (1.6, 3.6), ["AT/TT"] = (-2.7, -10.8),
        ["CT/GT"] = (-5.0, -15.8), ["GT/CT"] = (-2.2, -8.4), ["TT/AT"] = (0.2, -1.5)
    };

    /// <summary>
    /// Bommarito et al. (2000) single dangling-end nearest-neighbour parameters
    /// (ΔH° kcal/mol, ΔS° cal/(K·mol)) at 1 M NaCl. The '.' marks the unpaired (dangling)
    /// base. Key form for a 5'-side (left) dangling end is "topPair/bottomPair" of the
    /// first two columns; for a 3'-side (right) dangling end the reversed last two columns
    /// of each strand are used (per Tm_NN). Transcribed verbatim from Biopython
    /// <c>DNA_DE</c>; cross-checked term-by-term against SantaLucia &amp; Hicks (2004)
    /// Table 3 ΔH°. Source: Bommarito, Peyret &amp; SantaLucia (2000) NAR 28:1929.
    /// </summary>
    private static readonly Dictionary<string, (double DeltaH, double DeltaS)> NnDanglingEnd = new()
    {
        ["AA/.T"] = (0.2, 2.3), ["AC/.G"] = (-6.3, -17.1), ["AG/.C"] = (-3.7, -10.0), ["AT/.A"] = (-2.9, -7.6),
        ["CA/.T"] = (0.6, 3.3), ["CC/.G"] = (-4.4, -12.6), ["CG/.C"] = (-4.0, -11.9), ["CT/.A"] = (-4.1, -13.0),
        ["GA/.T"] = (-1.1, -1.6), ["GC/.G"] = (-5.1, -14.0), ["GG/.C"] = (-3.9, -10.9), ["GT/.A"] = (-4.2, -15.0),
        ["TA/.T"] = (-6.9, -20.0), ["TC/.G"] = (-4.0, -10.9), ["TG/.C"] = (-4.9, -13.8), ["TT/.A"] = (-0.2, -0.5),
        [".A/AT"] = (-0.7, -0.8), [".C/AG"] = (-2.1, -3.9), [".G/AC"] = (-5.9, -16.5), [".T/AA"] = (-0.5, -1.1),
        [".A/CT"] = (4.4, 14.9), [".C/CG"] = (-0.2, -0.1), [".G/CC"] = (-2.6, -7.4), [".T/CA"] = (4.7, 14.2),
        [".A/GT"] = (-1.6, -3.6), [".C/GG"] = (-3.9, -11.2), [".G/GC"] = (-3.2, -10.4), [".T/GA"] = (-4.1, -13.1),
        [".A/TT"] = (2.9, 10.4), [".C/TG"] = (-4.4, -13.1), [".G/TC"] = (-5.2, -15.0), [".T/TA"] = (-3.8, -12.6)
    };

    /// <summary>
    /// Computes the duplex ΔH° (kcal/mol) and ΔS° (cal/(K·mol)) for a probe–target DNA
    /// duplex that may contain a single internal mismatch and/or a single dangling end,
    /// using the SantaLucia (1998) Watson-Crick NN parameters together with the Allawi/
    /// SantaLucia/Peyret internal-mismatch and Bommarito (2000) dangling-end NN terms.
    /// Mirrors Biopython <c>Tm_NN(..., imm_table=DNA_IMM, de_table=DNA_DE)</c>.
    /// </summary>
    /// <param name="topStrand">Top strand 5'→3'. May start/end with a single '.' marking a
    /// dangling end on the bottom strand.</param>
    /// <param name="bottomStrand">Bottom strand written 3'→5' (the complement of the top
    /// read left-to-right, NOT the reverse complement), so base i pairs with top base i.
    /// May start/end with a single '.' marking a dangling end on the top strand.</param>
    /// <returns>(ΔH°, ΔS°, IsSelfComplementary) or <c>null</c> if the strands are null,
    /// unequal length, shorter than two columns, or contain a stack with no NN parameter
    /// (e.g. two adjacent mismatches, a tandem mismatch, or a non-ACGT character).</returns>
    public static (double DeltaH, double DeltaS, bool IsSelfComplementary)? CalculateNearestNeighborThermodynamicsMismatch(
        string topStrand, string bottomStrand)
    {
        if (topStrand is null || bottomStrand is null)
            return null;

        string top = topStrand.ToUpperInvariant();
        string bot = bottomStrand.ToUpperInvariant();
        if (top.Length != bot.Length || top.Length < 2)
            return null;

        double dH = NnInitDeltaH;
        double dS = NnInitDeltaS;

        // Terminal A·T penalty per end that closes with an A·T pair, using the (un-dotted)
        // top-strand termini exactly as Tm_NN computes `ends = seq[0] + seq[-1]`.
        if (top[0] is 'A' or 'T') { dH += NnTerminalAtDeltaH; dS += NnTerminalAtDeltaS; }
        if (top[^1] is 'A' or 'T') { dH += NnTerminalAtDeltaH; dS += NnTerminalAtDeltaS; }

        // Symmetry correction only for a fully paired self-complementary duplex.
        bool hasDangling = top.Contains('.') || bot.Contains('.');
        bool selfComp = !hasDangling && IsSelfComplementary(top)
                        && string.Equals(bot, Complement(top), StringComparison.Ordinal);
        if (selfComp)
            dS += NnSymmetryDeltaS;

        string ts = top, tc = bot;

        // Left (5'-side) dangling end.
        if (ts[0] == '.' || tc[0] == '.')
        {
            string leftDe = ts.Substring(0, 2) + "/" + tc.Substring(0, 2);
            if (!NnDanglingEnd.TryGetValue(leftDe, out var ld))
                return null;
            dH += ld.DeltaH; dS += ld.DeltaS;
            ts = ts.Substring(1); tc = tc.Substring(1);
        }

        // Right (3'-side) dangling end: reversed last two columns of each strand.
        if (ts[^1] == '.' || tc[^1] == '.')
        {
            string rightDe = Reverse(tc.Substring(tc.Length - 2)) + "/" + Reverse(ts.Substring(ts.Length - 2));
            if (!NnDanglingEnd.TryGetValue(rightDe, out var rd))
                return null;
            dH += rd.DeltaH; dS += rd.DeltaS;
            ts = ts.Substring(0, ts.Length - 1); tc = tc.Substring(0, tc.Length - 1);
        }

        // Nearest-neighbour stack over the paired region (Watson-Crick or internal mismatch).
        for (int i = 0; i < ts.Length - 1; i++)
        {
            string key = ts.Substring(i, 2) + "/" + tc.Substring(i, 2);
            if (!TryNnOrMismatch(key, out var p))
                return null;
            dH += p.DeltaH; dS += p.DeltaS;
        }

        return (dH, dS, selfComp);

        static bool TryNnOrMismatch(string key, out (double DeltaH, double DeltaS) p)
        {
            // key = "t0t1/b0b1" (bottom written 3'→5' aligned). The stack is a perfect
            // Watson-Crick step ONLY when both columns are complementary pairs; in that
            // case the perfect-match table (keyed by the top dinucleotide) applies. Any
            // non-WC column makes it an internal mismatch → use the mismatch table.
            char t0 = key[0], t1 = key[1], b0 = key[3], b1 = key[4];
            bool col0Wc = IsWatsonCrick(t0, b0);
            bool col1Wc = IsWatsonCrick(t1, b1);
            if (col0Wc && col1Wc)
            {
                string top2 = key.Substring(0, 2);
                if (NnUnifiedParams.TryGetValue(top2, out p)) return true;
                p = default;
                return false;
            }

            string rev = Reverse(key);
            if (NnInternalMismatch.TryGetValue(key, out p)) return true;
            if (NnInternalMismatch.TryGetValue(rev, out p)) return true;
            p = default;
            return false;
        }

        static bool IsWatsonCrick(char a, char b) =>
            (a == 'A' && b == 'T') || (a == 'T' && b == 'A') ||
            (a == 'G' && b == 'C') || (a == 'C' && b == 'G');
    }

    /// <summary>
    /// Computes the design melting temperature (°C) for a probe–target DNA duplex that may
    /// contain a single internal mismatch and/or a single dangling end, using the
    /// SantaLucia (1998) NN parameters plus the Allawi/SantaLucia/Peyret internal-mismatch
    /// and Bommarito (2000) dangling-end terms, the same bimolecular Tm equation and salt
    /// corrections as <see cref="CalculateMeltingTemperatureNN"/>. <b>Opt-in extension</b>:
    /// the perfect-match <see cref="CalculateMeltingTemperatureNN"/> is unchanged, and a
    /// fully paired duplex through this path equals it.
    /// </summary>
    /// <param name="topStrand">Top strand 5'→3' (may carry a leading/trailing '.' dangling-end marker).</param>
    /// <param name="bottomStrand">Bottom strand 3'→5', aligned base-for-base under the top
    /// (complement direction, NOT reverse complement; may carry a '.' dangling-end marker).</param>
    /// <param name="strandConcentrationMolar">Total strand concentration C_T in mol/L (default 0.5 µM).</param>
    /// <param name="sodiumMolar">Monovalent cation concentration in mol/L (default 50 mM).</param>
    /// <param name="magnesiumMolar">[Mg²⁺] in mol/L (default 0; only used by the divalent mode).</param>
    /// <param name="dntpMolar">Total dNTP concentration in mol/L (default 0).</param>
    /// <param name="saltMode">Salt correction to apply (default Owczarzy2004Monovalent).</param>
    /// <returns>The NN Tm in °C, or <c>double.NaN</c> if the duplex is not computable
    /// (null/unequal-length strands, &lt; 2 columns, or a stack with no NN parameter).</returns>
    public static double CalculateMeltingTemperatureNNMismatch(
        string topStrand,
        string bottomStrand,
        double strandConcentrationMolar = DefaultStrandConcentrationMolar,
        double sodiumMolar = ThermoConstants.DefaultNaConcentration,
        double magnesiumMolar = 0.0,
        double dntpMolar = 0.0,
        SaltCorrectionMode saltMode = SaltCorrectionMode.Owczarzy2004Monovalent)
    {
        var thermo = CalculateNearestNeighborThermodynamicsMismatch(topStrand, bottomStrand);
        if (thermo is null)
            return double.NaN;

        var (dH, dS, selfComp) = thermo.Value;

        // Use the paired-base count (excluding any dangling '.') for salt-correction length
        // and GC fraction, consistent with the duplex actually formed.
        string topPaired = topStrand.ToUpperInvariant().Replace(".", string.Empty);
        int length = topPaired.Length;
        double x = selfComp ? SelfComplementaryFactor : NonSelfComplementaryFactor;

        double dSeff = dS;
        if (saltMode == SaltCorrectionMode.SantaLuciaEntropy)
        {
            double phosphates = 2.0 * (length - 1);
            dSeff += SantaLuciaEntropySaltCoefficient * (phosphates / 2.0) * Math.Log(sodiumMolar);
        }

        double tmKelvin = (dH * 1000.0) / (dSeff + GasConstant * Math.Log(strandConcentrationMolar / x));

        switch (saltMode)
        {
            case SaltCorrectionMode.Owczarzy2004Monovalent:
                tmKelvin = ApplyOwczarzy2004(tmKelvin, topPaired, sodiumMolar);
                break;
            case SaltCorrectionMode.Owczarzy2008Divalent:
                tmKelvin = ApplyOwczarzy2008(tmKelvin, topPaired, sodiumMolar, magnesiumMolar, dntpMolar);
                break;
            case SaltCorrectionMode.None:
            case SaltCorrectionMode.SantaLuciaEntropy:
            default:
                break;
        }

        return tmKelvin - KelvinOffset;
    }

    // ---- LNA (locked nucleic acid)-adjusted NN Tm (PROBE-DESIGN-001, opt-in extension) ----
    // Extends the perfect-match DNA NN model with the McTigue, Peterson & Kahn (2004) LNA-DNA
    // nearest-neighbour increments so the NN Tm can be computed for a DNA oligo carrying one or
    // more INTERNAL LNA substitutions on one strand. The LNA value is an ADDITIVE increment
    // (ΔΔH°, ΔΔS°) added to the underlying DNA NN stack for each step containing the LNA base,
    // exactly as the MELTING 5 reference implementation realises it (McTigue04LockedAcid.java:
    // DNA NN sum, then `enthalpy += lockedAcidValue.getEnthalpy()`). The perfect-match
    // CalculateMeltingTemperatureNN above is UNCHANGED; this is opt-in.
    //
    // Convention: the increment for a step is keyed by the DNA dinucleotide step (e.g. "TT")
    // and which of the two bases of that step is locked (0 = the 5' base, 1 = the 3' base),
    // matching the paper's MX_L / X_L N notation (XML keys e.g. "TTL/AA" = step TT, 3' base
    // locked; "TLG/AC" = step TG, 5' base locked). Terminal LNA positions are NOT parameterised
    // by McTigue (2004) and are rejected (return not-computable), per MELTING isApplicable.
    //
    // Source (retrieved & cross-checked this session, 2026-06-24):
    //   McTigue PM, Peterson RJ, Kahn JD (2004) Biochemistry 43:5388-5405, DOI 10.1021/bi035976d
    //   — ΔΔH°/ΔΔS° for all 32 LNA+DNA:DNA nearest neighbours.
    //   Parameter values transcribed VERBATIM from the MELTING 5 data file
    //   "McTigue2004lockedmn.xml" (Dumousseau et al. 2012, BMC Bioinformatics 13:101; mirrored in
    //   aravind-j/rmelting). MELTING stores ΔΔH°/ΔΔS° in cal/mol and cal/(mol·K); the kcal/mol
    //   values below are XML_value / 1000. Worked example reproduced: CCATT(L)GCTACC at C=1e-4,
    //   Na=1 → ΔH°=-80.014 kcal/mol, ΔS°=-216.6 cal/(mol·K), Tm=63.528 °C (MELTING mct04: 63.614).

    /// <summary>Which base of a nearest-neighbour dinucleotide step is the LNA monomer.</summary>
    private enum LnaStepPosition
    {
        /// <summary>The 5' (first) base of the step is locked — McTigue X_L N (key e.g. "TLG/AC").</summary>
        FivePrime = 0,

        /// <summary>The 3' (second) base of the step is locked — McTigue MX_L (key e.g. "TTL/AA").</summary>
        ThreePrime = 1
    }

    /// <summary>
    /// McTigue, Peterson &amp; Kahn (2004) LNA-DNA nearest-neighbour increments
    /// (ΔΔH° in kcal/mol, ΔΔS° in cal/(K·mol)), added to the base DNA NN stack for the step
    /// containing the LNA base. Key = (DNA dinucleotide step 5'→3', which base is locked).
    /// All 32 nearest neighbours are present (16 with the 5' base locked, 16 with the 3' base
    /// locked). Values transcribed verbatim from MELTING 5 <c>McTigue2004lockedmn.xml</c>
    /// (cal/mol ÷ 1000 = kcal/mol); the XML <c>sequence</c> key is shown in the comment.
    /// Source: McTigue et al. (2004) Biochemistry 43:5388-5405 (DOI 10.1021/bi035976d).
    /// </summary>
    private static readonly Dictionary<(string Step, LnaStepPosition Locked), (double DeltaH, double DeltaS)> McTigueLnaIncrements = new()
    {
        // 5'-base locked (X_L N): XML key "XLY/comp".
        [("AA", LnaStepPosition.FivePrime)] = (0.707, 2.5),    // ALA/TT
        [("AT", LnaStepPosition.FivePrime)] = (2.282, 7.5),    // ALT/TA
        [("AG", LnaStepPosition.FivePrime)] = (0.264, 2.6),    // ALG/TC
        [("AC", LnaStepPosition.FivePrime)] = (1.131, 4.1),    // ALC/TG
        [("TA", LnaStepPosition.FivePrime)] = (-0.046, 1.6),   // TLA/AT
        [("TT", LnaStepPosition.FivePrime)] = (1.528, 5.3),    // TLT/AA
        [("TG", LnaStepPosition.FivePrime)] = (-1.540, -3.0),  // TLG/AC
        [("TC", LnaStepPosition.FivePrime)] = (1.893, 6.7),    // TLC/AG
        [("GA", LnaStepPosition.FivePrime)] = (3.162, 10.5),   // GLA/CT
        [("GT", LnaStepPosition.FivePrime)] = (-0.212, 0.1),   // GLT/CA
        [("GG", LnaStepPosition.FivePrime)] = (-2.844, -6.7),  // GLG/CC
        [("GC", LnaStepPosition.FivePrime)] = (-0.360, -0.3),  // GLC/CG
        [("CA", LnaStepPosition.FivePrime)] = (1.049, 4.3),    // CLA/GT
        [("CT", LnaStepPosition.FivePrime)] = (0.708, 4.2),    // CLT/GA
        [("CG", LnaStepPosition.FivePrime)] = (0.785, 3.7),    // CLG/GC
        [("CC", LnaStepPosition.FivePrime)] = (2.096, 8.0),    // CLC/GG

        // 3'-base locked (MX_L): XML key "MXL/comp".
        [("AA", LnaStepPosition.ThreePrime)] = (0.992, 4.1),   // AAL/TT
        [("AT", LnaStepPosition.ThreePrime)] = (1.816, 6.9),   // ATL/TA
        [("AG", LnaStepPosition.ThreePrime)] = (-1.200, -1.8), // AGL/TC
        [("AC", LnaStepPosition.ThreePrime)] = (2.890, 10.6),  // ACL/TG
        [("TA", LnaStepPosition.ThreePrime)] = (1.591, 5.3),   // TAL/AT
        [("TT", LnaStepPosition.ThreePrime)] = (2.326, 8.1),   // TTL/AA
        [("TG", LnaStepPosition.ThreePrime)] = (2.165, 7.2),   // TGL/AC
        [("TC", LnaStepPosition.ThreePrime)] = (0.609, 3.2),   // TCL/AG
        [("GA", LnaStepPosition.ThreePrime)] = (0.444, 2.9),   // GAL/CT
        [("GT", LnaStepPosition.ThreePrime)] = (-0.635, -0.3), // GTL/CA
        [("GG", LnaStepPosition.ThreePrime)] = (-0.943, -0.9), // GGL/CC
        [("GC", LnaStepPosition.ThreePrime)] = (-0.925, -1.1), // GCL/CG
        [("CA", LnaStepPosition.ThreePrime)] = (1.358, 4.4),   // CAL/GT
        [("CT", LnaStepPosition.ThreePrime)] = (-1.671, -4.1), // CTL/GA
        [("CG", LnaStepPosition.ThreePrime)] = (-0.276, -0.7), // CGL/GC
        [("CC", LnaStepPosition.ThreePrime)] = (2.063, 7.6)    // CCL/GG
    };

    /// <summary>
    /// Computes the duplex ΔH° (kcal/mol) and ΔS° (cal/(K·mol)) of a DNA oligonucleotide that
    /// carries one or more <b>internal</b> LNA (locked nucleic acid) substitutions, by adding the
    /// McTigue, Peterson &amp; Kahn (2004) LNA-DNA nearest-neighbour increments to the SantaLucia
    /// (1998) DNA NN stack (initiation + terminal-A·T + symmetry are computed on the underlying
    /// DNA sequence, unchanged). <b>Opt-in</b>: the perfect-match
    /// <see cref="CalculateNearestNeighborThermodynamics"/> is unchanged, and an empty
    /// <paramref name="lnaPositions"/> reproduces it exactly.
    /// </summary>
    /// <param name="sequence">DNA sequence (one strand, 5'→3'); the LNA monomers are at the given
    /// positions of this sequence. Must be ≥ 2 ACGT bases.</param>
    /// <param name="lnaPositions">Zero-based positions of the LNA monomers within
    /// <paramref name="sequence"/>. Order and duplicates are tolerated. A <b>terminal</b> position
    /// (0 or length−1) is not parameterised by McTigue (2004) and makes the result not computable.</param>
    /// <returns>(ΔH°, ΔS°, IsSelfComplementary) of the LNA-substituted duplex, or <c>null</c> if the
    /// sequence is empty/&lt; 2 bases/contains a non-ACGT base, or any LNA position is out of range
    /// or terminal.</returns>
    public static (double DeltaH, double DeltaS, bool IsSelfComplementary)? CalculateNearestNeighborThermodynamicsLna(
        string sequence,
        IReadOnlyCollection<int> lnaPositions)
    {
        ArgumentNullException.ThrowIfNull(lnaPositions);

        var dna = CalculateNearestNeighborThermodynamics(sequence);
        if (dna is null)
            return null;

        string seq = sequence.ToUpperInvariant();
        var locked = new HashSet<int>();
        foreach (int pos in lnaPositions)
        {
            // McTigue (2004) parameters are for internal LNA only — reject terminal/out-of-range.
            if (pos <= 0 || pos >= seq.Length - 1)
                return null;
            locked.Add(pos);
        }

        var (dH, dS, selfComp) = dna.Value;

        // Add the McTigue increment to each NN step (i, i+1) that contains an LNA base.
        for (int i = 0; i < seq.Length - 1; i++)
        {
            string step = seq.Substring(i, 2);
            if (locked.Contains(i)
                && McTigueLnaIncrements.TryGetValue((step, LnaStepPosition.FivePrime), out var inc5))
            {
                dH += inc5.DeltaH; dS += inc5.DeltaS;
            }
            if (locked.Contains(i + 1)
                && McTigueLnaIncrements.TryGetValue((step, LnaStepPosition.ThreePrime), out var inc3))
            {
                dH += inc3.DeltaH; dS += inc3.DeltaS;
            }
        }

        return (dH, dS, selfComp);
    }

    /// <summary>
    /// Computes the design melting temperature (°C) of a DNA oligonucleotide carrying one or more
    /// <b>internal</b> LNA substitutions, using the McTigue (2004) LNA-DNA nearest-neighbour
    /// increments on top of the SantaLucia (1998) DNA NN model, with the same bimolecular Tm
    /// equation and optional salt corrections as <see cref="CalculateMeltingTemperatureNN"/>.
    /// <b>Opt-in</b>: the perfect-match <see cref="CalculateMeltingTemperatureNN"/> is unchanged,
    /// and an empty <paramref name="lnaPositions"/> equals it exactly.
    /// </summary>
    /// <param name="sequence">DNA sequence (5'→3'). Must be ≥ 2 ACGT bases.</param>
    /// <param name="lnaPositions">Zero-based positions of the internal LNA monomers (see
    /// <see cref="CalculateNearestNeighborThermodynamicsLna"/>).</param>
    /// <param name="strandConcentrationMolar">Total strand concentration C_T in mol/L (default 0.5 µM).</param>
    /// <param name="sodiumMolar">Monovalent cation concentration in mol/L (default 50 mM).</param>
    /// <param name="magnesiumMolar">[Mg²⁺] in mol/L (default 0; only used by the divalent mode).</param>
    /// <param name="dntpMolar">Total dNTP concentration in mol/L (default 0).</param>
    /// <param name="saltMode">Salt correction to apply (default Owczarzy2004Monovalent).</param>
    /// <returns>The LNA-adjusted NN Tm in °C, or <c>double.NaN</c> if the duplex is not computable
    /// (empty/&lt; 2 bases/non-ACGT, or an out-of-range/terminal LNA position).</returns>
    public static double CalculateMeltingTemperatureNNLna(
        string sequence,
        IReadOnlyCollection<int> lnaPositions,
        double strandConcentrationMolar = DefaultStrandConcentrationMolar,
        double sodiumMolar = ThermoConstants.DefaultNaConcentration,
        double magnesiumMolar = 0.0,
        double dntpMolar = 0.0,
        SaltCorrectionMode saltMode = SaltCorrectionMode.Owczarzy2004Monovalent)
    {
        var thermo = CalculateNearestNeighborThermodynamicsLna(sequence, lnaPositions);
        if (thermo is null)
            return double.NaN;

        var (dH, dS, selfComp) = thermo.Value;
        int length = sequence.Length;
        double x = selfComp ? SelfComplementaryFactor : NonSelfComplementaryFactor;

        double dSeff = dS;
        if (saltMode == SaltCorrectionMode.SantaLuciaEntropy)
        {
            double phosphates = 2.0 * (length - 1);
            dSeff += SantaLuciaEntropySaltCoefficient * (phosphates / 2.0) * Math.Log(sodiumMolar);
        }

        double tmKelvin = (dH * 1000.0) / (dSeff + GasConstant * Math.Log(strandConcentrationMolar / x));

        switch (saltMode)
        {
            case SaltCorrectionMode.Owczarzy2004Monovalent:
                tmKelvin = ApplyOwczarzy2004(tmKelvin, sequence.ToUpperInvariant(), sodiumMolar);
                break;
            case SaltCorrectionMode.Owczarzy2008Divalent:
                tmKelvin = ApplyOwczarzy2008(tmKelvin, sequence.ToUpperInvariant(), sodiumMolar, magnesiumMolar, dntpMolar);
                break;
            case SaltCorrectionMode.None:
            case SaltCorrectionMode.SantaLuciaEntropy:
            default:
                break;
        }

        return tmKelvin - KelvinOffset;
    }

    /// <summary>Watson-Crick complement of an ACGT string (same 5'→3'/left-to-right order).</summary>
    private static string Complement(string seq)
    {
        var sb = new StringBuilder(seq.Length);
        foreach (char c in seq)
        {
            sb.Append(c switch
            {
                'A' => 'T',
                'T' => 'A',
                'G' => 'C',
                'C' => 'G',
                _ => c
            });
        }
        return sb.ToString();
    }

    /// <summary>GC fraction over A/C/G/T bases only (denominator excludes non-ACGT).</summary>
    private static double GcFraction(string seq)
    {
        int gc = 0, valid = 0;
        foreach (char c in seq)
        {
            if (c is 'G' or 'C') { gc++; valid++; }
            else if (c is 'A' or 'T') valid++;
        }
        return valid == 0 ? 0.0 : (double)gc / valid;
    }

    // Owczarzy (2004) monovalent correction in 1/Tm form (Kelvin).
    private static double ApplyOwczarzy2004(double tmKelvin, string seq, double sodiumMolar)
    {
        double lnNa = Math.Log(sodiumMolar);
        double fgc = GcFraction(seq);
        double corr = (Owczarzy2004GcCoefficient * fgc - Owczarzy2004Constant) * lnNa
                      + Owczarzy2004QuadraticCoefficient * lnNa * lnNa;
        return 1.0 / (1.0 / tmKelvin + corr);
    }

    // Owczarzy (2008) divalent (Mg²⁺/dNTP) correction in 1/Tm form (Kelvin),
    // reproducing Biopython salt_correction method 7.
    private static double ApplyOwczarzy2008(
        double tmKelvin, string seq, double sodiumMolar, double magnesiumMolar, double dntpMolar)
    {
        double mon = sodiumMolar;
        double mg = magnesiumMolar;

        // Free Mg²⁺ after dNTP chelation (quadratic solution with Ka).
        if (dntpMolar > 0)
        {
            double ka = DntpMgAssociationConstant;
            mg = (-(ka * dntpMolar - ka * mg + 1.0)
                  + Math.Sqrt((ka * dntpMolar - ka * mg + 1.0) * (ka * dntpMolar - ka * mg + 1.0)
                              + 4.0 * ka * mg)) / (2.0 * ka);
        }

        double fgc = GcFraction(seq);
        double corr;

        // If essentially no divalent ion, fall back to the monovalent 2004 form.
        if (mg <= 0 && mon > 0)
        {
            double ln = Math.Log(mon);
            corr = (Owczarzy2004GcCoefficient * fgc - Owczarzy2004Constant) * ln
                   + Owczarzy2004QuadraticCoefficient * ln * ln;
            return 1.0 / (1.0 / tmKelvin + corr);
        }

        double a = Owc2008A, b = Owc2008B, c = Owc2008C, d = Owc2008D, e = Owc2008E, f = Owc2008F, g = Owc2008G;
        double r = mon > 0 ? Math.Sqrt(mg) / mon : double.PositiveInfinity;

        if (mon > 0 && r < Owc2008MonovalentRatioLow)
        {
            // Monovalent dominates → 2004 form.
            double ln = Math.Log(mon);
            corr = (Owczarzy2004GcCoefficient * fgc - Owczarzy2004Constant) * ln
                   + Owczarzy2004QuadraticCoefficient * ln * ln;
            return 1.0 / (1.0 / tmKelvin + corr);
        }

        if (mon > 0 && r < Owc2008MonovalentRatioHigh)
        {
            // Mixed regime: reparameterise a, d, g (Owczarzy 2008 / Biopython method 7).
            double lnMon = Math.Log(mon);
            a = 3.92e-5 * (0.843 - 0.352 * Math.Sqrt(mon) * lnMon);
            d = 1.42e-5 * (1.279 - 4.03e-3 * lnMon - 8.03e-3 * lnMon * lnMon);
            g = 8.31e-5 * (0.486 - 0.258 * lnMon + 5.25e-3 * lnMon * lnMon * lnMon);
        }

        double lnMg = Math.Log(mg);
        corr = a + b * lnMg + fgc * (c + d * lnMg)
               + (1.0 / (2.0 * (seq.Length - 1))) * (e + f * lnMg + g * lnMg * lnMg);
        return 1.0 / (1.0 / tmKelvin + corr);
    }

    /// <summary>True if the sequence equals its own reverse complement (self-complementary).</summary>
    private static bool IsSelfComplementary(string seq)
    {
        int n = seq.Length;
        if (n % 2 != 0) return false; // odd-length cannot be self-complementary
        for (int i = 0; i < n; i++)
        {
            char a = seq[i];
            char b = seq[n - 1 - i];
            bool pair = (a == 'A' && b == 'T') || (a == 'T' && b == 'A')
                     || (a == 'G' && b == 'C') || (a == 'C' && b == 'G');
            if (!pair) return false;
        }
        return true;
    }

    /// <summary>
    /// Generates all possible primers for a region.
    /// </summary>
    public static IEnumerable<PrimerCandidate> GeneratePrimerCandidates(
        DnaSequence template,
        int regionStart,
        int regionEnd,
        bool forward = true,
        PrimerParameters? parameters = null)
    {
        var param = parameters ?? DefaultParameters;

        for (int start = regionStart; start + param.MinLength <= regionEnd; start++)
        {
            for (int len = param.MinLength; len <= param.MaxLength && start + len <= regionEnd; len++)
            {
                var seq = template.Sequence.Substring(start, len);
                if (!forward)
                    seq = new DnaSequence(seq).ReverseComplement().Sequence;

                yield return EvaluatePrimer(seq, start, forward, param);
            }
        }
    }

    // ---- Primer3 weighted penalty objective (PRIMER-TM-001) -------------------
    // Reproduces the per-primer objective function `p_obj_fn` (left/right primer
    // branch) from Primer3's reference source `libprimer3.cc`, with the documented
    // default weights and optima. Lower penalty = better primer, exactly as Primer3.
    // Source: Primer3 source libprimer3.cc p_obj_fn / pr_set_default_global_args_2;
    //         Primer3 manual §19 "HOW PRIMER3 CALCULATES THE PENALTY VALUE";
    //         Untergasser et al. (2012) NAR 40(15):e115; Koressaar & Remm (2007).

    /// <summary>
    /// Primer3 default per-primer objective weights, taken verbatim from
    /// <c>pr_set_default_global_args_2</c> in Primer3's <c>libprimer3.cc</c>:
    /// length/Tm weights = 1; GC, self-complementarity and N weights = 0.
    /// Source: Primer3 source (branch main), function pr_set_default_global_args_2.
    /// </summary>
    public static readonly Primer3PenaltyWeights DefaultPrimer3Weights = new(
        TmGt: 1.0,           // PRIMER_WT_TM_GT  (libprimer3.cc: weights.temp_gt = 1)
        TmLt: 1.0,           // PRIMER_WT_TM_LT  (weights.temp_lt = 1)
        SizeGt: 1.0,         // PRIMER_WT_SIZE_GT (weights.length_gt = 1)
        SizeLt: 1.0,         // PRIMER_WT_SIZE_LT (weights.length_lt = 1)
        GcGt: 0.0,           // PRIMER_WT_GC_PERCENT_GT (weights.gc_content_gt = 0)
        GcLt: 0.0,           // PRIMER_WT_GC_PERCENT_LT (weights.gc_content_lt = 0)
        SelfAny: 0.0,        // PRIMER_WT_SELF_ANY (weights.compl_any = 0)
        SelfEnd: 0.0,        // PRIMER_WT_SELF_END (weights.compl_end = 0)
        NumNs: 0.0           // PRIMER_WT_NUM_NS  (weights.num_ns = 0)
    );

    /// <summary>
    /// Primer3 default per-primer optima. OPT_TM = 60 °C and OPT_SIZE = 20 bases are
    /// from <c>libprimer3.cc</c> (opt_tm = 60.0, opt_size = 20); OPT_GC_PERCENT = 50.0 %
    /// is the published default in the Primer3 manual (PRIMER_OPT_GC_PERCENT, default 50.0).
    /// </summary>
    public static readonly Primer3Optima DefaultPrimer3Optima = new(
        OptTm: 60.0,         // PRIMER_OPT_TM (libprimer3.cc: opt_tm = 60.0) °C
        OptSize: 20,         // PRIMER_OPT_SIZE (libprimer3.cc: opt_size = 20) bases
        OptGcPercent: 50.0   // PRIMER_OPT_GC_PERCENT (manual default 50.0) %
    );

    /// <summary>
    /// Computes the Primer3 per-primer penalty (objective function value) for a single
    /// primer, faithfully reproducing the left/right-primer branch of Primer3's
    /// <c>p_obj_fn</c>. The penalty is the weighted sum of one-sided deviations of Tm,
    /// length and GC% from their optima, plus the weighted self-/3'-complementarity and
    /// number-of-Ns terms. Each term is added only when its weight is non-zero and the
    /// deviation has the matching sign, so the result is always ≥ 0; <b>lower is better</b>,
    /// exactly as Primer3 sorts candidates.
    /// </summary>
    /// <param name="inputs">Measured primer properties (Tm in °C, length in bases, GC in
    /// percent 0–100, self/3' local-alignment scores, count of N bases).</param>
    /// <param name="weights">Objective weights; defaults to <see cref="DefaultPrimer3Weights"/>.</param>
    /// <param name="optima">Parameter optima; defaults to <see cref="DefaultPrimer3Optima"/>.</param>
    /// <returns>The Primer3 objective-function value (penalty); 0 means every term is at its optimum.</returns>
    public static double CalculatePrimer3Penalty(
        Primer3PenaltyInputs inputs,
        Primer3PenaltyWeights? weights = null,
        Primer3Optima? optima = null)
    {
        var w = weights ?? DefaultPrimer3Weights;
        var o = optima ?? DefaultPrimer3Optima;

        double sum = 0.0;

        // Tm term: one-sided, separate _gt / _lt weights (p_obj_fn temp_gt / temp_lt).
        if (w.TmGt != 0 && inputs.Tm > o.OptTm)
            sum += w.TmGt * (inputs.Tm - o.OptTm);
        if (w.TmLt != 0 && inputs.Tm < o.OptTm)
            sum += w.TmLt * (o.OptTm - inputs.Tm);

        // GC% term (gc_content is a percentage 0–100 in libprimer3.cc, line 3856).
        if (w.GcGt != 0 && inputs.GcPercent > o.OptGcPercent)
            sum += w.GcGt * (inputs.GcPercent - o.OptGcPercent);
        if (w.GcLt != 0 && inputs.GcPercent < o.OptGcPercent)
            sum += w.GcLt * (o.OptGcPercent - inputs.GcPercent);

        // Length/size term (p_obj_fn length_lt / length_gt).
        if (w.SizeLt != 0 && inputs.Length < o.OptSize)
            sum += w.SizeLt * (o.OptSize - inputs.Length);
        if (w.SizeGt != 0 && inputs.Length > o.OptSize)
            sum += w.SizeGt * (inputs.Length - o.OptSize);

        // Self-complementarity terms (non-thermodynamic branch: compl_any, compl_end).
        if (w.SelfAny != 0)
            sum += w.SelfAny * inputs.SelfAny;
        if (w.SelfEnd != 0)
            sum += w.SelfEnd * inputs.SelfEnd;

        // Number-of-Ns term (num_ns).
        if (w.NumNs != 0)
            sum += w.NumNs * inputs.NumNs;

        return sum;
    }

    private static double CalculatePrimerScore(string seq, double gc, double tm, int homopolymer, PrimerParameters param)
    {
        double score = 100;

        // Penalize for deviation from optimal length
        score -= Math.Abs(seq.Length - param.OptimalLength) * 2;

        // Penalize for deviation from optimal Tm
        score -= Math.Abs(tm - param.OptimalTm) * 2;

        // Penalize for deviation from 50% GC
        score -= Math.Abs(gc - 50) * 0.5;

        // Penalize for homopolymers
        score -= homopolymer * 5;

        // Bonus for GC clamp at 3' end
        if (seq.Length >= 2)
        {
            char last = seq[^1];
            if (last == 'G' || last == 'C')
                score += 5;
        }

        return Math.Max(0, score);
    }

    private static string Reverse(string s)
    {
        var chars = s.ToCharArray();
        Array.Reverse(chars);
        return new string(chars);
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
/// Parameters for primer design.
/// </summary>
public readonly record struct PrimerParameters(
    int MinLength,
    int MaxLength,
    int OptimalLength,
    double MinGcContent,
    double MaxGcContent,
    double MinTm,
    double MaxTm,
    double OptimalTm,
    int MaxHomopolymer,
    int MaxDinucleotideRepeats,
    bool Avoid3PrimeGC,
    bool Check3PrimeStability);

/// <summary>
/// A primer candidate with quality metrics.
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
    double Score);

/// <summary>
/// Measured properties of a single primer used as input to the Primer3 penalty
/// objective (<see cref="PrimerDesigner.CalculatePrimer3Penalty"/>). Units mirror
/// Primer3's <c>p_obj_fn</c>: Tm in °C, length in bases, GC in percent (0–100),
/// self/3'-complementarity as local-alignment scores, and the count of N bases.
/// </summary>
/// <param name="Tm">Primer melting temperature in °C.</param>
/// <param name="Length">Primer length in bases.</param>
/// <param name="GcPercent">GC content as a percentage in [0, 100].</param>
/// <param name="SelfAny">Self-complementarity local-alignment score (PRIMER_SELF_ANY); 0 if unused.</param>
/// <param name="SelfEnd">3'-self-complementarity local-alignment score (PRIMER_SELF_END); 0 if unused.</param>
/// <param name="NumNs">Number of ambiguous N bases in the primer.</param>
public readonly record struct Primer3PenaltyInputs(
    double Tm,
    int Length,
    double GcPercent,
    double SelfAny = 0.0,
    double SelfEnd = 0.0,
    int NumNs = 0);

/// <summary>
/// Weights for the Primer3 per-primer penalty objective (the <c>PRIMER_WT_*</c>
/// parameters). Tm, GC and size each have separate "greater-than" (_gt) and
/// "less-than" (_lt) weights, applied one-sidedly relative to the optimum.
/// Defaults are <see cref="PrimerDesigner.DefaultPrimer3Weights"/>.
/// </summary>
public readonly record struct Primer3PenaltyWeights(
    double TmGt,
    double TmLt,
    double SizeGt,
    double SizeLt,
    double GcGt,
    double GcLt,
    double SelfAny,
    double SelfEnd,
    double NumNs);

/// <summary>
/// Optimal parameter values for the Primer3 penalty objective
/// (<c>PRIMER_OPT_TM</c>, <c>PRIMER_OPT_SIZE</c>, <c>PRIMER_OPT_GC_PERCENT</c>).
/// Defaults are <see cref="PrimerDesigner.DefaultPrimer3Optima"/>.
/// </summary>
public readonly record struct Primer3Optima(
    double OptTm,
    int OptSize,
    double OptGcPercent);

/// <summary>
/// Result of primer pair design.
/// </summary>
public sealed record PrimerPairResult(
    PrimerCandidate? Forward,
    PrimerCandidate? Reverse,
    bool IsValid,
    string Message,
    int ProductSize);
