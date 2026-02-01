using System;
using System.Collections.Generic;
using System.Linq;
using System.Threading;

namespace Seqeron.Genomics
{
    /// <summary>
    /// Performs approximate pattern matching with support for mismatches, insertions, and deletions.
    /// </summary>
    public static class ApproximateMatcher
    {
        /// <summary>
        /// Finds all approximate matches of a pattern in a sequence with at most k mismatches.
        /// Uses Hamming distance (substitutions only).
        /// </summary>
        /// <param name="sequence">The sequence to search in.</param>
        /// <param name="pattern">The pattern to find.</param>
        /// <param name="maxMismatches">Maximum number of allowed mismatches.</param>
        /// <returns>Enumerable of match results.</returns>
        public static IEnumerable<ApproximateMatchResult> FindWithMismatches(
            string sequence, string pattern, int maxMismatches)
        {
            if (string.IsNullOrEmpty(sequence) || string.IsNullOrEmpty(pattern))
                yield break;

            if (maxMismatches < 0)
                throw new ArgumentOutOfRangeException(nameof(maxMismatches), "Cannot be negative.");

            var seq = sequence.ToUpperInvariant();
            var pat = pattern.ToUpperInvariant();

            if (pat.Length > seq.Length)
                yield break;

            for (int i = 0; i <= seq.Length - pat.Length; i++)
            {
                int mismatches = 0;
                var positions = new List<int>();

                for (int j = 0; j < pat.Length && mismatches <= maxMismatches; j++)
                {
                    if (seq[i + j] != pat[j])
                    {
                        mismatches++;
                        positions.Add(j);
                    }
                }

                if (mismatches <= maxMismatches)
                {
                    yield return new ApproximateMatchResult(
                        i,
                        seq.Substring(i, pat.Length),
                        mismatches,
                        positions.AsReadOnly(),
                        MismatchType.Substitution
                    );
                }
            }
        }

        /// <summary>
        /// Finds all approximate matches with cancellation support.
        /// </summary>
        /// <param name="sequence">The sequence to search in.</param>
        /// <param name="pattern">The pattern to find.</param>
        /// <param name="maxMismatches">Maximum number of allowed mismatches.</param>
        /// <param name="cancellationToken">Cancellation token for long-running operations.</param>
        /// <returns>Enumerable of match results.</returns>
        public static IEnumerable<ApproximateMatchResult> FindWithMismatches(
            string sequence,
            string pattern,
            int maxMismatches,
            CancellationToken cancellationToken)
        {
            return CancellableOperations.FindWithMismatches(sequence, pattern, maxMismatches, cancellationToken);
        }

        /// <summary>
        /// Finds all approximate matches of a pattern in a DNA sequence with at most k mismatches.
        /// </summary>
        public static IEnumerable<ApproximateMatchResult> FindWithMismatches(
            DnaSequence sequence, string pattern, int maxMismatches)
        {
            return FindWithMismatches(sequence.Sequence, pattern, maxMismatches);
        }

        /// <summary>
        /// Finds all approximate matches in a DNA sequence with cancellation support.
        /// </summary>
        public static IEnumerable<ApproximateMatchResult> FindWithMismatches(
            DnaSequence sequence,
            string pattern,
            int maxMismatches,
            CancellationToken cancellationToken)
        {
            return FindWithMismatches(sequence.Sequence, pattern, maxMismatches, cancellationToken);
        }

        /// <summary>
        /// Finds all approximate matches using edit distance (Levenshtein distance).
        /// Allows substitutions, insertions, and deletions.
        /// </summary>
        /// <param name="sequence">The sequence to search in.</param>
        /// <param name="pattern">The pattern to find.</param>
        /// <param name="maxEdits">Maximum edit distance allowed.</param>
        /// <returns>Enumerable of match results.</returns>
        public static IEnumerable<ApproximateMatchResult> FindWithEdits(
            string sequence, string pattern, int maxEdits)
        {
            if (string.IsNullOrEmpty(sequence) || string.IsNullOrEmpty(pattern))
                yield break;

            if (maxEdits < 0)
                throw new ArgumentOutOfRangeException(nameof(maxEdits), "Cannot be negative.");

            var seq = sequence.ToUpperInvariant();
            var pat = pattern.ToUpperInvariant();

            // Use sliding window with variable length (pattern ± maxEdits)
            int minLen = Math.Max(1, pat.Length - maxEdits);
            int maxLen = pat.Length + maxEdits;

            for (int i = 0; i <= seq.Length - minLen; i++)
            {
                for (int len = minLen; len <= maxLen && i + len <= seq.Length; len++)
                {
                    string window = seq.Substring(i, len);
                    int distance = EditDistance(pat, window);

                    if (distance <= maxEdits)
                    {
                        yield return new ApproximateMatchResult(
                            i,
                            window,
                            distance,
                            Array.Empty<int>().ToList().AsReadOnly(),
                            (window.Length == pat.Length && distance == HammingDistanceFast(pat, window))
                                ? MismatchType.Substitution
                                : MismatchType.Edit
                        );
                    }
                }
            }
        }

        /// <summary>
        /// Finds all approximate matches using edit distance in a DNA sequence.
        /// </summary>
        public static IEnumerable<ApproximateMatchResult> FindWithEdits(
            DnaSequence sequence, string pattern, int maxEdits)
        {
            return FindWithEdits(sequence.Sequence, pattern, maxEdits);
        }

        /// <summary>
        /// Calculates the Hamming distance between two strings of equal length.
        /// </summary>
        /// <param name="s1">First string.</param>
        /// <param name="s2">Second string.</param>
        /// <returns>Number of positions with different characters.</returns>
        public static int HammingDistance(string s1, string s2)
        {
            if (s1 == null || s2 == null)
                throw new ArgumentNullException(s1 == null ? nameof(s1) : nameof(s2));

            if (s1.Length != s2.Length)
                throw new ArgumentException("Strings must have equal length for Hamming distance.");

            int distance = 0;
            for (int i = 0; i < s1.Length; i++)
            {
                if (char.ToUpperInvariant(s1[i]) != char.ToUpperInvariant(s2[i]))
                    distance++;
            }
            return distance;
        }

        /// <summary>
        /// Calculates the edit distance (Levenshtein distance) between two strings.
        /// </summary>
        /// <param name="s1">First string.</param>
        /// <param name="s2">Second string.</param>
        /// <returns>Minimum number of edits (insertions, deletions, substitutions) needed.</returns>
        public static int EditDistance(string s1, string s2)
        {
            if (s1 == null || s2 == null)
                throw new ArgumentNullException(s1 == null ? nameof(s1) : nameof(s2));

            s1 = s1.ToUpperInvariant();
            s2 = s2.ToUpperInvariant();

            int m = s1.Length;
            int n = s2.Length;

            // Optimize for edge cases
            if (m == 0) return n;
            if (n == 0) return m;

            // Use two rows instead of full matrix for memory efficiency
            var prev = new int[n + 1];
            var curr = new int[n + 1];

            // Initialize first row
            for (int j = 0; j <= n; j++)
                prev[j] = j;

            for (int i = 1; i <= m; i++)
            {
                curr[0] = i;

                for (int j = 1; j <= n; j++)
                {
                    int cost = s1[i - 1] == s2[j - 1] ? 0 : 1;
                    curr[j] = Math.Min(
                        Math.Min(
                            prev[j] + 1,      // deletion
                            curr[j - 1] + 1   // insertion
                        ),
                        prev[j - 1] + cost    // substitution
                    );
                }

                // Swap rows
                (prev, curr) = (curr, prev);
            }

            return prev[n];
        }

        /// <summary>
        /// Finds the best approximate match (minimum distance) of a pattern in a sequence.
        /// </summary>
        /// <param name="sequence">The sequence to search in.</param>
        /// <param name="pattern">The pattern to find.</param>
        /// <returns>The best match result, or null if sequence is too short.</returns>
        public static ApproximateMatchResult? FindBestMatch(string sequence, string pattern)
        {
            if (string.IsNullOrEmpty(sequence) || string.IsNullOrEmpty(pattern))
                return null;

            var seq = sequence.ToUpperInvariant();
            var pat = pattern.ToUpperInvariant();

            if (pat.Length > seq.Length)
                return null;

            ApproximateMatchResult? best = null;
            int bestDistance = int.MaxValue;

            // Check exact length windows
            for (int i = 0; i <= seq.Length - pat.Length; i++)
            {
                string window = seq.Substring(i, pat.Length);
                int distance = HammingDistance(pat, window);

                if (distance < bestDistance)
                {
                    bestDistance = distance;
                    var positions = new List<int>();
                    for (int j = 0; j < pat.Length; j++)
                    {
                        if (window[j] != pat[j])
                            positions.Add(j);
                    }

                    best = new ApproximateMatchResult(
                        i, window, distance, positions.AsReadOnly(), MismatchType.Substitution
                    );

                    if (distance == 0)
                        return best; // Perfect match found
                }
            }

            return best;
        }

        /// <summary>
        /// Counts the number of approximate occurrences of a pattern in a sequence.
        /// </summary>
        public static int CountApproximateOccurrences(string sequence, string pattern, int maxMismatches)
        {
            return FindWithMismatches(sequence, pattern, maxMismatches).Count();
        }

        /// <summary>
        /// Finds the most frequent approximate k-mers (with up to d mismatches).
        /// </summary>
        /// <param name="sequence">The sequence to analyze.</param>
        /// <param name="k">K-mer length.</param>
        /// <param name="d">Maximum mismatches.</param>
        /// <returns>List of most frequent k-mers with their counts.</returns>
        public static IEnumerable<(string Kmer, int Count)> FindFrequentKmersWithMismatches(
            string sequence, int k, int d)
        {
            if (string.IsNullOrEmpty(sequence))
                yield break;

            if (k <= 0)
                throw new ArgumentOutOfRangeException(nameof(k), "K must be positive.");

            if (d < 0)
                throw new ArgumentOutOfRangeException(nameof(d), "D cannot be negative.");

            var seq = sequence.ToUpperInvariant();
            var counts = new Dictionary<string, int>();

            // For each k-mer in sequence, count all patterns that match it with ≤d mismatches
            for (int i = 0; i <= seq.Length - k; i++)
            {
                string kmer = seq.Substring(i, k);

                // Generate all patterns within d mismatches
                foreach (string neighbor in GenerateNeighbors(kmer, d))
                {
                    if (!counts.TryAdd(neighbor, 1))
                        counts[neighbor]++;
                }
            }

            // Find maximum count
            int maxCount = counts.Values.Max();

            // Return all k-mers with maximum count
            foreach (var kvp in counts.Where(kvp => kvp.Value == maxCount))
            {
                yield return (kvp.Key, kvp.Value);
            }
        }

        /// <summary>
        /// Generates all DNA neighbors within d mismatches.
        /// </summary>
        private static IEnumerable<string> GenerateNeighbors(string pattern, int d)
        {
            if (d == 0)
            {
                yield return pattern;
                yield break;
            }

            if (pattern.Length == 1)
            {
                foreach (char c in "ACGT")
                    yield return c.ToString();
                yield break;
            }

            char first = pattern[0];
            string suffix = pattern.Substring(1);

            foreach (string neighborSuffix in GenerateNeighbors(suffix, d))
            {
                if (HammingDistanceFast(suffix, neighborSuffix) < d)
                {
                    // Can change first character
                    foreach (char c in "ACGT")
                        yield return c + neighborSuffix;
                }
                else
                {
                    // Must keep first character
                    yield return first + neighborSuffix;
                }
            }
        }

        private static int HammingDistanceFast(string s1, string s2)
        {
            int distance = 0;
            int minLen = Math.Min(s1.Length, s2.Length);
            for (int i = 0; i < minLen; i++)
            {
                if (s1[i] != s2[i])
                    distance++;
            }
            return distance + Math.Abs(s1.Length - s2.Length);
        }
    }

    /// <summary>
    /// Type of mismatch in approximate matching.
    /// </summary>
    public enum MismatchType
    {
        /// <summary>Only substitutions (Hamming distance).</summary>
        Substitution,

        /// <summary>Insertions, deletions, and substitutions (edit distance).</summary>
        Edit
    }

    /// <summary>
    /// Result of an approximate pattern match.
    /// </summary>
    public readonly record struct ApproximateMatchResult(
        int Position,
        string MatchedSequence,
        int Distance,
        IReadOnlyList<int> MismatchPositions,
        MismatchType MismatchType)
    {
        /// <summary>
        /// Gets whether this is an exact match (distance = 0).
        /// </summary>
        public bool IsExact => Distance == 0;
    }
}
