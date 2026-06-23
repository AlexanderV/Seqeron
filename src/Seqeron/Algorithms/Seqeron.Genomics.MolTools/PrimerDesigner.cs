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
