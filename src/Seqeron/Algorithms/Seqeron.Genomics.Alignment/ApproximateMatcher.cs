namespace Seqeron.Genomics.Alignment
{
    /// <summary>
    /// Performs approximate pattern matching with support for mismatches, insertions, and deletions.
    /// </summary>
    public static class ApproximateMatcher
    {
        // DNA alphabet over which the d-neighborhood is generated.
        // Per Compeau & Pevzner, Bioinformatics Algorithms ch.1 (ROSALIND BA1N),
        // Neighbors(Pattern, d) enumerates substitutions over {A, C, G, T}.
        private const string DnaAlphabet = "ACGT";

        /// <summary>
        /// Finds all approximate matches of a pattern in a sequence with at most k mismatches.
        /// Uses Hamming distance (substitutions only, no indels): every 0-based start
        /// i ∈ [0, n − m] whose length-m window has Hamming distance ≤ k from the pattern is
        /// reported, overlapping occurrences included (Navarro 2001 "k mismatches";
        /// Compeau &amp; Pevzner ch.1, ROSALIND BA1H). Comparison is case-insensitive and
        /// MatchedSequence is returned upper-cased; MismatchPositions are pattern-relative.
        /// A null/empty sequence or pattern, or a pattern longer than the sequence, yields no matches.
        /// </summary>
        /// <param name="sequence">The sequence to search in.</param>
        /// <param name="pattern">The pattern to find.</param>
        /// <param name="maxMismatches">Maximum number of allowed mismatches.</param>
        /// <returns>Enumerable of match results.</returns>
        public static IEnumerable<ApproximateMatchResult> FindWithMismatches(
            string sequence, string pattern, int maxMismatches)
        {
            return FindWithMismatches(sequence, pattern, maxMismatches, CancellationToken.None);
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
            if (maxMismatches < 0)
                throw new ArgumentOutOfRangeException(nameof(maxMismatches), "Cannot be negative.");
            return FindWithMismatchesCore(sequence, pattern, maxMismatches, cancellationToken);
        }

        private static IEnumerable<ApproximateMatchResult> FindWithMismatchesCore(
            string sequence,
            string pattern,
            int maxMismatches,
            CancellationToken cancellationToken)
        {
            if (string.IsNullOrEmpty(sequence) || string.IsNullOrEmpty(pattern))
                yield break;

            var seq = sequence.ToUpperInvariant();
            var pat = pattern.ToUpperInvariant();

            if (pat.Length > seq.Length)
                yield break;

            const int checkInterval = 1000;
            for (int i = 0; i <= seq.Length - pat.Length; i++)
            {
                if (i % checkInterval == 0)
                    cancellationToken.ThrowIfCancellationRequested();

                var positions = CollectMismatchPositions(seq, i, pat, maxMismatches);
                int mismatches = positions.Count;

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
        /// Finds all approximate matches using edit distance (Levenshtein distance; unit-cost
        /// substitutions, insertions and deletions). Every substring T[i..i+len) of the
        /// upper-cased sequence with ed(pattern, T[i..i+len)) ≤ <paramref name="maxEdits"/> is
        /// reported, ordered by start position then window length (len ∈ [max(1, m − k), m + k];
        /// empty windows are excluded). The set of window end positions i+len−1 equals the
        /// Sellers (1980) end-position set returned by <see cref="FindEditEndPositions(string, string, int)"/>.
        /// Implementation: for each start i one start-anchored Wagner–Fischer/Sellers column DP
        /// (Navarro 2001 §5.1) yields the distance to every window length in a single pass,
        /// with Ukkonen's (1985) cut-off once the column minimum exceeds k (column minima never
        /// decrease), i.e. O(n·m·(m+k)) instead of re-running a full DP per window.
        /// MismatchType is Substitution when the window has the pattern's length and the edit
        /// distance equals the Hamming distance, otherwise Edit; MismatchPositions is empty.
        /// </summary>
        /// <param name="sequence">The sequence to search in.</param>
        /// <param name="pattern">The pattern to find.</param>
        /// <param name="maxEdits">Maximum edit distance allowed.</param>
        /// <returns>Enumerable of match results.</returns>
        public static IEnumerable<ApproximateMatchResult> FindWithEdits(
            string sequence, string pattern, int maxEdits)
        {
            if (maxEdits < 0)
                throw new ArgumentOutOfRangeException(nameof(maxEdits), "Cannot be negative.");
            return FindWithEditsCore(sequence, pattern, maxEdits);
        }

        private static IEnumerable<ApproximateMatchResult> FindWithEditsCore(
            string sequence, string pattern, int maxEdits)
        {
            if (string.IsNullOrEmpty(sequence) || string.IsNullOrEmpty(pattern))
                yield break;

            var seq = sequence.ToUpperInvariant();
            var pat = pattern.ToUpperInvariant();
            int m = pat.Length;

            // Window lengths considered: pattern ± maxEdits (no longer/shorter window can be
            // within maxEdits, since ed ≥ |length difference|); empty windows excluded.
            int minLen = Math.Max(1, m - maxEdits);
            int maxLen = m + maxEdits;

            var prev = new int[m + 1];
            var curr = new int[m + 1];

            for (int i = 0; i <= seq.Length - minLen; i++)
            {
                // Column 0 of the DP anchored at start i: ed(pat[0..r), "") = r.
                for (int r = 0; r <= m; r++)
                    prev[r] = r;

                for (int len = 1; len <= maxLen && i + len <= seq.Length; len++)
                {
                    int columnMin = AdvanceColumn(pat, seq[i + len - 1], prev, curr, len);
                    (prev, curr) = (curr, prev);

                    int distance = prev[m]; // = ed(pat, seq[i..i+len))
                    if (len >= minLen && distance <= maxEdits)
                    {
                        string window = seq.Substring(i, len);
                        yield return new ApproximateMatchResult(
                            i,
                            window,
                            distance,
                            Array.Empty<int>(),
                            (len == m && distance == pat.AsSpan().HammingDistance(window.AsSpan()))
                                ? MismatchType.Substitution
                                : MismatchType.Edit
                        );
                    }

                    // Ukkonen cut-off: column minima are non-decreasing, so no longer window
                    // from this start can come back within maxEdits.
                    if (columnMin > maxEdits)
                        break;
                }
            }
        }

        /// <summary>
        /// Sellers (1980) approximate string matching ("k differences" problem; Navarro 2001 §5.1):
        /// reports every 0-based end position j of the (upper-cased) sequence at which some
        /// substring ending at j is within edit distance <paramref name="maxEdits"/> of the
        /// pattern, together with that minimum distance C[m, j] = min_i ed(pattern, T[i..j]).
        /// The DP is the Wagner–Fischer recurrence with a free start in the text (C[0, j] = 0).
        /// A null/empty sequence or pattern yields no matches. O(n·m) time, O(m) space.
        /// </summary>
        /// <param name="sequence">The text to search in.</param>
        /// <param name="pattern">The pattern to find.</param>
        /// <param name="maxEdits">Maximum edit distance allowed (k ≥ 0).</param>
        /// <returns>(EndPosition, Distance) pairs in increasing end position.</returns>
        public static IEnumerable<(int EndPosition, int Distance)> FindEditEndPositions(
            string sequence, string pattern, int maxEdits)
        {
            if (maxEdits < 0)
                throw new ArgumentOutOfRangeException(nameof(maxEdits), "Cannot be negative.");
            return FindEditEndPositionsCore(sequence, pattern, maxEdits);
        }

        private static IEnumerable<(int EndPosition, int Distance)> FindEditEndPositionsCore(
            string sequence, string pattern, int maxEdits)
        {
            if (string.IsNullOrEmpty(sequence) || string.IsNullOrEmpty(pattern))
                yield break;

            var seq = sequence.ToUpperInvariant();
            var pat = pattern.ToUpperInvariant();
            int m = pat.Length;

            var prev = new int[m + 1];
            var curr = new int[m + 1];
            for (int r = 0; r <= m; r++)
                prev[r] = r;

            for (int j = 0; j < seq.Length; j++)
            {
                // Free start in the text: C[0, j] = 0 (Sellers 1980).
                AdvanceColumn(pat, seq[j], prev, curr, 0);
                (prev, curr) = (curr, prev);

                if (prev[m] <= maxEdits)
                    yield return (j, prev[m]);
            }
        }

        /// <summary>
        /// One column step of the unit-cost edit-distance DP (Wagner &amp; Fischer 1974):
        /// given column <paramref name="prev"/> (pattern prefixes vs text up to the previous
        /// character) computes column <paramref name="curr"/> for text character
        /// <paramref name="c"/> with top cell <paramref name="top"/>
        /// (= column index for global distance, 0 for Sellers' free text start).
        /// Returns the column minimum.
        /// </summary>
        private static int AdvanceColumn(string pat, char c, int[] prev, int[] curr, int top)
        {
            curr[0] = top;
            int min = top;
            for (int r = 1; r <= pat.Length; r++)
            {
                int cost = pat[r - 1] == c ? 0 : 1;
                int v = Math.Min(
                    Math.Min(
                        prev[r] + 1,      // text character unmatched (insertion into pattern)
                        curr[r - 1] + 1   // pattern character unmatched (deletion from pattern)
                    ),
                    prev[r - 1] + cost    // match / substitution
                );
                curr[r] = v;
                if (v < min)
                    min = v;
            }
            return min;
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
        /// Calculates the Hamming distance between two strings of equal length
        /// (number of positions whose symbols differ; Hamming 1950; ROSALIND HAMM).
        /// Comparison is case-insensitive (ordinal, invariant upper-casing).
        /// </summary>
        /// <exception cref="ArgumentNullException">Either string is null.</exception>
        /// <exception cref="ArgumentException">The strings differ in length (the distance is undefined).</exception>
        /// <param name="s1">First string.</param>
        /// <param name="s2">Second string.</param>
        /// <returns>Number of positions with different characters.</returns>
        public static int HammingDistance(string s1, string s2)
        {
            if (s1 == null || s2 == null)
                throw new ArgumentNullException(s1 == null ? nameof(s1) : nameof(s2));

            if (s1.Length != s2.Length)
                throw new ArgumentException("Strings must have equal length for Hamming distance.");

            // Delegate to the canonical case-insensitive span implementation
            // (SequenceExtensions.HammingDistance) rather than duplicating the count.
            return s1.AsSpan().HammingDistance(s2.AsSpan());
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

            int m = s1.Length;
            int n = s2.Length;

            // Optimize for edge cases
            if (m == 0) return n;
            if (n == 0) return m;

            // Two-column Wagner–Fischer DP (s1 = rows, s2 = columns) sharing the
            // column kernel used by FindWithEdits / FindEditEndPositions.
            var prev = new int[m + 1];
            var curr = new int[m + 1];
            for (int r = 0; r <= m; r++)
                prev[r] = r;

            for (int j = 1; j <= n; j++)
            {
                AdvanceColumn(s1, s2[j - 1], prev, curr, j);
                (prev, curr) = (curr, prev);
            }

            return prev[m];
        }

        /// <summary>
        /// Finds the best approximate match (minimum Hamming distance) of a pattern over all
        /// equal-length windows of the sequence. When several windows tie for the minimum
        /// distance, the leftmost one is returned (an exact match short-circuits the scan).
        /// Distance follows the Hamming definition of an approximate occurrence in
        /// Compeau &amp; Pevzner, Bioinformatics Algorithms ch.1 (ROSALIND BA1H).
        /// </summary>
        /// <param name="sequence">The sequence to search in.</param>
        /// <param name="pattern">The pattern to find.</param>
        /// <returns>The best match result, or null if either input is empty or the pattern is longer than the sequence.</returns>
        public static ApproximateMatchResult? FindBestMatch(string sequence, string pattern)
        {
            if (string.IsNullOrEmpty(sequence) || string.IsNullOrEmpty(pattern))
                return null;

            var seq = sequence.ToUpperInvariant();
            var pat = pattern.ToUpperInvariant();

            if (pat.Length > seq.Length)
                return null;

            int m = pat.Length;
            int bestPosition = -1;
            int bestDistance = int.MaxValue;

            // Leftmost window with the strictly smallest Hamming distance (canonical span
            // Hamming, no per-window allocation); an exact match cannot be beaten.
            for (int i = 0; i <= seq.Length - m && bestDistance > 0; i++)
            {
                int distance = seq.AsSpan(i, m).HammingDistance(pat.AsSpan());
                if (distance < bestDistance)
                {
                    bestDistance = distance;
                    bestPosition = i;
                }
            }

            string window = seq.Substring(bestPosition, m);
            var positions = CollectMismatchPositions(seq, bestPosition, pat, m);

            return new ApproximateMatchResult(
                bestPosition, window, bestDistance, positions.AsReadOnly(), MismatchType.Substitution);
        }

        /// <summary>
        /// Pattern-relative indices j at which <paramref name="seq"/>[start + j] ≠ <paramref name="pat"/>[j]
        /// (both already upper-cased), in ascending order. Scanning stops as soon as more than
        /// <paramref name="limit"/> mismatches have been seen (the list then holds limit + 1 entries),
        /// so a caller can reject the window by <c>Count &gt; limit</c>. Shared by
        /// <see cref="FindWithMismatches(string, string, int, CancellationToken)"/> and <see cref="FindBestMatch"/>.
        /// </summary>
        private static List<int> CollectMismatchPositions(string seq, int start, string pat, int limit)
        {
            var positions = new List<int>();
            for (int j = 0; j < pat.Length && positions.Count <= limit; j++)
            {
                if (seq[start + j] != pat[j])
                    positions.Add(j);
            }
            return positions;
        }

        /// <summary>
        /// Counts the number of approximate occurrences of a pattern in a sequence, i.e.
        /// Count_d(Text, Pattern) — the number of start positions where the equal-length window
        /// has Hamming distance at most <paramref name="maxMismatches"/> from the pattern
        /// (Compeau &amp; Pevzner, Bioinformatics Algorithms ch.1; ROSALIND BA1H/BA1I).
        /// </summary>
        public static int CountApproximateOccurrences(string sequence, string pattern, int maxMismatches)
        {
            return FindWithMismatches(sequence, pattern, maxMismatches).Count();
        }

        /// <summary>
        /// Finds the most frequent k-mers with up to d mismatches (the Frequent Words with
        /// Mismatches Problem; Compeau &amp; Pevzner, Bioinformatics Algorithms ch.1, ROSALIND BA1I).
        /// Each window's full d-neighborhood (Hamming ball, which includes the window itself) is
        /// tallied; ALL k-mers achieving the maximum Count_d are returned. The result may include
        /// k-mers that do not occur exactly in the sequence. Input is case-insensitive; reported
        /// k-mers are DNA k-mers over {A, C, G, T} only (BA1N), so a window's non-ACGT symbol always
        /// costs one mismatch. Windows are first counted with the canonical k-mer counter
        /// (SequenceExtensions.CountKmersSpan) and each distinct window's neighborhood is weighted
        /// by its multiplicity. Returns nothing when no DNA k-mer has a positive Count_d.
        /// Result order is unspecified.
        /// </summary>
        /// <param name="sequence">The sequence to analyze.</param>
        /// <param name="k">K-mer length.</param>
        /// <param name="d">Maximum mismatches.</param>
        /// <returns>All most-frequent k-mers (ties included) with their Count_d.</returns>
        public static IEnumerable<(string Kmer, int Count)> FindFrequentKmersWithMismatches(
            string sequence, int k, int d)
        {
            if (k <= 0)
                throw new ArgumentOutOfRangeException(nameof(k), "K must be positive.");
            if (d < 0)
                throw new ArgumentOutOfRangeException(nameof(d), "D cannot be negative.");
            return FindFrequentKmersWithMismatchesCore(sequence, k, d);
        }

        private static IEnumerable<(string Kmer, int Count)> FindFrequentKmersWithMismatchesCore(
            string sequence, int k, int d)
        {
            if (string.IsNullOrEmpty(sequence))
                yield break;

            // Distinct windows with multiplicities from the canonical (case-insensitive) k-mer
            // counter, so each distinct window's d-neighborhood is generated once and weighted by
            // its multiplicity: Count_d(P) = Σ_w mult(w)·[HD(P, w) ≤ d] (BA1I tally).
            var windowCounts = sequence.AsSpan().CountKmersSpan(k);
            var counts = new Dictionary<string, int>();
            var buffer = new char[k];

            foreach (var (window, multiplicity) in windowCounts)
            {
                foreach (string neighbor in GenerateNeighbors(window, d, buffer))
                {
                    counts[neighbor] = counts.TryGetValue(neighbor, out int c) ? c + multiplicity : multiplicity;
                }
            }

            // No length-k window could be cut (e.g. k > sequence length), or no DNA k-mer lies
            // within d of any window (all windows carry ≥ d+1 non-ACGT symbols): every Count_d is
            // 0, so there is no most-frequent k-mer and the result is empty.
            if (counts.Count == 0)
                yield break;

            int maxCount = counts.Values.Max();

            // Return all k-mers with maximum count
            foreach (var kvp in counts.Where(kvp => kvp.Value == maxCount))
            {
                yield return (kvp.Key, kvp.Value);
            }
        }

        /// <summary>
        /// Neighbors(Pattern, d) — the d-neighborhood of Pattern: every k-mer over {A, C, G, T}
        /// whose Hamming distance from Pattern does not exceed d (Compeau &amp; Pevzner ch.1,
        /// ROSALIND BA1N), each exactly once. Positions are filled left to right, spending one
        /// unit of the mismatch budget per substituted base. For an ACGT pattern this is the same
        /// set as the textbook recursive Neighbors (identity included); unlike that recursion,
        /// which assumes an ACGT pattern and would emit non-DNA strings such as ANA for ANG, a
        /// non-ACGT symbol can never be kept and always costs one mismatch.
        /// </summary>
        private static IEnumerable<string> GenerateNeighbors(string pattern, int d, char[] buffer)
        {
            var result = new List<string>();
            CollectNeighbors(pattern, 0, d, buffer, result);
            return result;
        }

        private static void CollectNeighbors(string pattern, int position, int budget, char[] buffer, List<string> result)
        {
            if (position == pattern.Length)
            {
                result.Add(new string(buffer, 0, pattern.Length));
                return;
            }

            char original = pattern[position];
            foreach (char c in DnaAlphabet)
            {
                int cost = c == original ? 0 : 1;
                if (cost > budget)
                    continue;
                buffer[position] = c;
                CollectNeighbors(pattern, position + 1, budget - cost, buffer, result);
            }
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

