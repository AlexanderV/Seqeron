using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace SuffixTree.Genomics
{
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
            MinTm: 55,
            MaxTm: 65,
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
        /// Calculates the melting temperature using the nearest-neighbor method.
        /// Simplified version using basic formula for primers.
        /// </summary>
        public static double CalculateMeltingTemperature(string primer)
        {
            if (string.IsNullOrEmpty(primer))
                return 0;

            var seq = primer.ToUpperInvariant();

            // For short primers (< 14 bp), use Wallace rule
            if (seq.Length < 14)
            {
                int at = seq.Count(c => c == 'A' || c == 'T');
                int gc = seq.Count(c => c == 'G' || c == 'C');
                return 2 * at + 4 * gc;
            }

            // For longer primers, use basic nearest-neighbor approximation
            // Tm = 64.9 + 41 * (G+C-16.4) / N
            int gcCount = seq.Count(c => c == 'G' || c == 'C');
            double tm = 64.9 + 41.0 * (gcCount - 16.4) / seq.Length;

            return Math.Max(0, tm);
        }

        /// <summary>
        /// Calculates the melting temperature with salt correction.
        /// </summary>
        /// <param name="primer">Primer sequence.</param>
        /// <param name="naConcentration">Na+ concentration in mM (default: 50).</param>
        /// <returns>Corrected melting temperature in °C.</returns>
        public static double CalculateMeltingTemperatureWithSalt(string primer, double naConcentration = 50)
        {
            double baseTm = CalculateMeltingTemperature(primer);
            
            // Salt correction: Tm = Tm_basic + 16.6 * log10([Na+])
            double saltCorrection = 16.6 * Math.Log10(naConcentration / 1000.0);
            
            return Math.Round(baseTm + saltCorrection, 1);
        }

        /// <summary>
        /// Calculates GC content as a percentage.
        /// </summary>
        public static double CalculateGcContent(string sequence)
        {
            if (string.IsNullOrEmpty(sequence))
                return 0;

            var seq = sequence.ToUpperInvariant();
            int gcCount = seq.Count(c => c == 'G' || c == 'C');
            return (double)gcCount / seq.Length * 100;
        }

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
        /// </summary>
        public static bool HasHairpinPotential(string sequence, int minStemLength = 4)
        {
            if (string.IsNullOrEmpty(sequence) || sequence.Length < minStemLength * 2 + 3)
                return false;

            var seq = sequence.ToUpperInvariant();
            var revComp = GetReverseComplement(seq);

            // Check for self-complementary regions
            for (int i = 0; i <= seq.Length - minStemLength; i++)
            {
                string fragment = seq.Substring(i, minStemLength);
                // Look for complementary sequence at least 3 positions away
                for (int j = i + minStemLength + 3; j <= seq.Length - minStemLength; j++)
                {
                    string target = seq.Substring(j, minStemLength);
                    if (AreComplementary(fragment, Reverse(target)))
                        return true;
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
            var seq2 = GetReverseComplement(primer2.ToUpperInvariant());

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
        /// Calculates the stability of the 3' end (last 5 bases).
        /// More negative = more stable = potentially problematic.
        /// </summary>
        public static double Calculate3PrimeStability(string sequence)
        {
            if (string.IsNullOrEmpty(sequence) || sequence.Length < 5)
                return 0;

            var seq = sequence.ToUpperInvariant();
            string last5 = seq.Substring(seq.Length - 5);

            // Simplified ΔG calculation using nearest-neighbor values
            // Values are ΔG in kcal/mol for each dinucleotide
            var deltaG = new Dictionary<string, double>
            {
                ["AA"] = -1.0, ["TT"] = -1.0,
                ["AT"] = -0.88, ["TA"] = -0.58,
                ["CA"] = -1.45, ["TG"] = -1.45,
                ["GT"] = -1.44, ["AC"] = -1.44,
                ["CT"] = -1.28, ["AG"] = -1.28,
                ["GA"] = -1.30, ["TC"] = -1.30,
                ["CG"] = -2.17, ["GC"] = -2.24,
                ["GG"] = -1.84, ["CC"] = -1.84
            };

            double totalDeltaG = 0;
            for (int i = 0; i < last5.Length - 1; i++)
            {
                string dinuc = last5.Substring(i, 2);
                if (deltaG.TryGetValue(dinuc, out double dg))
                    totalDeltaG += dg;
            }

            return totalDeltaG;
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

        private static string GetReverseComplement(string sequence)
        {
            var sb = new StringBuilder(sequence.Length);
            for (int i = sequence.Length - 1; i >= 0; i--)
            {
                sb.Append(sequence[i] switch
                {
                    'A' => 'T', 'T' => 'A',
                    'C' => 'G', 'G' => 'C',
                    _ => sequence[i]
                });
            }
            return sb.ToString();
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
    /// Result of primer pair design.
    /// </summary>
    public sealed record PrimerPairResult(
        PrimerCandidate? Forward,
        PrimerCandidate? Reverse,
        bool IsValid,
        string Message,
        int ProductSize);
}
