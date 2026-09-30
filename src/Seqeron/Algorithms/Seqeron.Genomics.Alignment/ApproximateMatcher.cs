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
        /// Each hit carries its optimal alignment (<see cref="ApproximateMatchResult.Alignment"/>:
        /// pattern = query, window = target; edlib extended CIGAR '=', 'X', 'I', 'D') recovered by
        /// traceback over the same DP columns with the deterministic tie-break of
        /// <see cref="GetEditAlignment(string, string)"/> (diagonal, then 'I', then 'D', walking back
        /// from the end). MismatchPositions lists the pattern-relative indices of the substituted
        /// ('X') pattern characters of that alignment. MismatchType is Substitution when the
        /// alignment has no indel — which, under this tie-break, holds exactly when the window has
        /// the pattern's length and the edit distance equals the Hamming distance (then the
        /// MismatchPositions are the Hamming mismatch indices) — otherwise Edit.
        /// Space: O(m·(m+k)) for the start-anchored DP columns kept for the traceback.
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
            // (long arithmetic: maxEdits may be as large as int.MaxValue.)
            int minLen = Math.Max(1, m - maxEdits);
            int maxLen = (int)Math.Min((long)m + maxEdits, seq.Length);

            // Start-anchored DP columns kept for the traceback: cols[len] = column after
            // T[i..i+len); cols[0] is the fixed boundary column ed(pat[0..r), "") = r.
            var cols = new int[maxLen + 1][];
            for (int len = 0; len <= maxLen; len++)
                cols[len] = new int[m + 1];
            for (int r = 0; r <= m; r++)
                cols[0][r] = r;

            for (int i = 0; i <= seq.Length - minLen; i++)
            {
                for (int len = 1; len <= maxLen && i + len <= seq.Length; len++)
                {
                    int columnMin = AdvanceColumn(pat, seq[i + len - 1], cols[len - 1], cols[len], len);

                    int distance = cols[len][m]; // = ed(pat, seq[i..i+len))
                    if (len >= minLen && distance <= maxEdits)
                    {
                        var alignment = Traceback(pat, seq, i, len, cols);
                        yield return new ApproximateMatchResult(
                            i,
                            seq.Substring(i, len),
                            distance,
                            alignment.SubstitutionPositions,
                            alignment.HasIndels ? MismatchType.Edit : MismatchType.Substitution
                        )
                        { Alignment = alignment };
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
        /// A null/empty sequence or pattern yields no matches. Engine: Myers (1999) bit-parallel
        /// column scan (⌈m/64⌉ words, O(⌈m/64⌉·n) time), identical to the DP reference
        /// <see cref="FindEditEndPositionsDp"/>.
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

            // Myers (1999) bit-parallel engine with the Sellers boundary (horizontal delta 0 in
            // row 0 ⇔ C[0, j] = 0); column scores are identical to the DP reference
            // FindEditEndPositionsDp (locked by exhaustive/random equality tests).
            var myers = new MyersBitVector(pat, global: false);
            for (int j = 0; j < seq.Length; j++)
            {
                int score = myers.Advance(seq[j]);
                if (score <= maxEdits)
                    yield return (j, score);
            }
        }

        /// <summary>
        /// Reference (non-bit-parallel) Sellers DP over the shared column kernel
        /// <see cref="AdvanceColumn"/>: the oracle against which the Myers engine of
        /// <see cref="FindEditEndPositions(string, string, int)"/> is tested. Same contract.
        /// </summary>
        internal static IEnumerable<(int EndPosition, int Distance)> FindEditEndPositionsDp(
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
        /// Calculates the edit distance (Levenshtein distance) between two strings
        /// (unit-cost substitutions, insertions, deletions; case-sensitive ordinal comparison).
        /// Engine: Myers (1999) bit-parallel algorithm in the global (Needleman–Wunsch) form of
        /// Hyyrö (2003) — horizontal delta +1 entering row 0 of every column — with ⌈m/64⌉ 64-bit
        /// blocks (the shorter string is the bit-vector "pattern"), as in edlib's
        /// <c>calculateBlock</c>. O(⌈min(m,n)/64⌉·max(m,n)) time. Results are identical to the
        /// Wagner–Fischer reference <see cref="EditDistanceDp"/> (exhaustive/random equality tests).
        /// </summary>
        /// <param name="s1">First string.</param>
        /// <param name="s2">Second string.</param>
        /// <returns>Minimum number of edits (insertions, deletions, substitutions) needed.</returns>
        public static int EditDistance(string s1, string s2)
        {
            if (s1 == null || s2 == null)
                throw new ArgumentNullException(s1 == null ? nameof(s1) : nameof(s2));

            // Distance is symmetric: the shorter string becomes the bit-vector pattern.
            var (pattern, text) = s1.Length <= s2.Length ? (s1, s2) : (s2, s1);
            if (pattern.Length == 0)
                return text.Length;

            var myers = new MyersBitVector(pattern, global: true);
            int score = pattern.Length;
            foreach (char c in text)
                score = myers.Advance(c);
            return score;
        }

        /// <summary>
        /// Reference two-column Wagner–Fischer DP (s1 = rows, s2 = columns) over the shared column
        /// kernel <see cref="AdvanceColumn"/>: the oracle against which the Myers engine of
        /// <see cref="EditDistance(string, string)"/> is tested. Same contract.
        /// </summary>
        internal static int EditDistanceDp(string s1, string s2)
        {
            if (s1 == null || s2 == null)
                throw new ArgumentNullException(s1 == null ? nameof(s1) : nameof(s2));

            int m = s1.Length;
            int n = s2.Length;
            if (m == 0) return n;
            if (n == 0) return m;

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
        /// Optimal global Levenshtein alignment of <paramref name="query"/> against
        /// <paramref name="target"/> (edlib NW mode; query = rows, target = columns). The
        /// operations follow edlib's convention: '=' match, 'X' mismatch, 'I' query character
        /// absent from the target (insertion to the query / deletion from the target), 'D' target
        /// character absent from the query; <see cref="EditAlignment.Cigar"/> is the edlib
        /// EXTENDED CIGAR and <see cref="EditAlignment.StandardCigar"/> the STANDARD one ('M' for
        /// '=' and 'X'). Among co-optimal paths the traceback — walking back from (m, n) — takes
        /// the diagonal move whenever it is optimal, then 'I', then 'D'. (edlib's traceback tries
        /// 'I', then 'D', then the diagonal, so it may return a different co-optimal path with the
        /// same cost; the diagonal-first rule is what makes a substitution-only optimum come back
        /// as the Hamming path.) Case-sensitive, like <see cref="EditDistance"/>; the distance
        /// equals <see cref="EditDistance"/>. Full Wagner–Fischer matrix: O(m·n) time and space.
        /// </summary>
        /// <exception cref="ArgumentNullException">Either string is null.</exception>
        public static EditAlignment GetEditAlignment(string query, string target)
        {
            if (query == null || target == null)
                throw new ArgumentNullException(query == null ? nameof(query) : nameof(target));

            int m = query.Length;
            int n = target.Length;
            var cols = new int[n + 1][];
            cols[0] = new int[m + 1];
            for (int r = 0; r <= m; r++)
                cols[0][r] = r;
            for (int j = 1; j <= n; j++)
            {
                cols[j] = new int[m + 1];
                AdvanceColumn(query, target[j - 1], cols[j - 1], cols[j], j);
            }

            return Traceback(query, target, 0, n, cols);
        }

        /// <summary>
        /// Traceback over DP columns <paramref name="cols"/>[0..len] of <paramref name="query"/>
        /// against <paramref name="text"/>[start..start+len) whose row 0 and column 0 are the
        /// global boundaries (C[0, j] = j, C[r, 0] = r). From (m, len) back to (0, 0) it takes
        /// the diagonal ('=' / 'X') when C[r−1, j−1] + cost = C[r, j], else 'I' when
        /// C[r−1, j] + 1 = C[r, j], else 'D' (C[r, j−1] + 1 = C[r, j] then holds by the recurrence).
        /// </summary>
        private static EditAlignment Traceback(string query, string text, int start, int len, int[][] cols)
        {
            int r = query.Length;
            int j = len;
            var ops = new char[r + len];
            int p = ops.Length;
            while (r > 0 || j > 0)
            {
                int v = cols[j][r];
                if (r > 0 && j > 0)
                {
                    bool same = query[r - 1] == text[start + j - 1];
                    if (cols[j - 1][r - 1] + (same ? 0 : 1) == v)
                    {
                        ops[--p] = same ? '=' : 'X';
                        r--;
                        j--;
                        continue;
                    }
                }

                if (r > 0 && cols[j][r - 1] + 1 == v)
                {
                    ops[--p] = 'I';
                    r--;
                }
                else
                {
                    ops[--p] = 'D';
                    j--;
                }
            }

            return new EditAlignment(
                cols[len][query.Length],
                new string(ops, p, ops.Length - p),
                query,
                text.Substring(start, len));
        }

        /// <summary>
        /// Myers (1999, J. ACM 46(3):395–415) bit-vector edit-distance column engine, block-based
        /// (⌈m/64⌉ words, bit r of block b = pattern row 64b + r + 1), transcribed from edlib's
        /// <c>calculateBlock</c> (Šošić &amp; Šikić 2017). Vertical deltas start at +1 (column 0:
        /// C[r, 0] = r). The horizontal delta entering row 0 is +1 for global distance
        /// (Hyyrö 2003: C[0, j] = j) and 0 for Sellers search (C[0, j] = 0); each block's
        /// outgoing delta feeds the next block. The score C[m, j] is tracked through the
        /// horizontal delta at row m (bit (m − 1) mod 64 of the last block); rows below m in the
        /// last word never influence rows above (carries and shifts only move upward).
        /// </summary>
        private sealed class MyersBitVector
        {
            private readonly ulong[]?[] _peqAscii = new ulong[]?[128];
            private readonly Dictionary<char, ulong[]> _peqOther = new();
            private readonly ulong[] _zero;
            private readonly ulong[] _pv;
            private readonly ulong[] _mv;
            private readonly int _lastBit;
            private readonly int _topDelta;
            private int _score;

            public MyersBitVector(string pattern, bool global)
            {
                int m = pattern.Length;
                int blocks = (m + 63) / 64;
                _zero = new ulong[blocks];
                _pv = new ulong[blocks];
                _mv = new ulong[blocks];
                Array.Fill(_pv, ulong.MaxValue);
                _lastBit = (m - 1) % 64;
                _topDelta = global ? 1 : 0;
                _score = m;

                // Peq[c] bit r ⇔ pattern[r] == c (Myers 1999 §3).
                for (int r = 0; r < m; r++)
                {
                    char c = pattern[r];
                    ulong[] eq;
                    if (c < 128)
                        eq = _peqAscii[c] ??= new ulong[blocks];
                    else if (!_peqOther.TryGetValue(c, out eq!))
                        _peqOther[c] = eq = new ulong[blocks];
                    eq[r / 64] |= 1UL << (r % 64);
                }
            }

            /// <summary>Advances one text column; returns C[m, j].</summary>
            public int Advance(char c)
            {
                ulong[]? eq;
                if (c < 128)
                    eq = _peqAscii[c];
                else
                    _peqOther.TryGetValue(c, out eq);
                eq ??= _zero;

                int h = _topDelta;
                int last = _pv.Length - 1;
                for (int b = 0; b <= last; b++)
                    h = AdvanceBlock(ref _pv[b], ref _mv[b], eq[b], h, b == last ? _lastBit : 63);

                _score += h;
                return _score;
            }

            /// <summary>
            /// edlib <c>calculateBlock</c> (Myers' Advance_Block with Hyyrö's horizontal input):
            /// <paramref name="hin"/> ∈ {−1, 0, +1} enters at the block's top row; returns the
            /// horizontal delta at bit <paramref name="outBit"/>.
            /// </summary>
            private static int AdvanceBlock(ref ulong pv, ref ulong mv, ulong eq, int hin, int outBit)
            {
                ulong hinIsNeg = hin < 0 ? 1UL : 0UL;
                ulong xv = eq | mv;
                eq |= hinIsNeg;
                ulong xh = (((eq & pv) + pv) ^ pv) | eq;
                ulong ph = mv | ~(xh | pv);
                ulong mh = pv & xh;

                int hout = (int)((ph >> outBit) & 1UL) - (int)((mh >> outBit) & 1UL);

                ph <<= 1;
                mh <<= 1;
                mh |= hinIsNeg;
                if (hin > 0)
                    ph |= 1UL;

                pv = mh | ~(xv | ph);
                mv = ph & xv;
                return hout;
            }
        }

        /// <summary>
        /// Optimal string alignment distance (restricted Damerau–Levenshtein; Damerau 1964,
        /// Boytsov 2011): unit-cost insertion, deletion, substitution and transposition of two
        /// adjacent characters, with the restriction that no substring is edited more than once
        /// (so it violates the triangle inequality: OSA(CA, ABC) = 3). Recurrence:
        /// d[i, j] = min(d[i−1, j] + 1, d[i, j−1] + 1, d[i−1, j−1] + [a_i ≠ b_j],
        /// d[i−2, j−2] + 1 if a_i = b_{j−1} and a_{i−1} = b_j). Case-sensitive (ordinal).
        /// O(m·n) time, O(n) space. Reference: rapidfuzz.distance.OSA.
        /// </summary>
        /// <exception cref="ArgumentNullException">Either string is null.</exception>
        public static int OptimalStringAlignmentDistance(string s1, string s2)
        {
            if (s1 == null || s2 == null)
                throw new ArgumentNullException(s1 == null ? nameof(s1) : nameof(s2));

            int m = s1.Length;
            int n = s2.Length;
            if (m == 0) return n;
            if (n == 0) return m;

            // Three rolling rows over s2: d[i−2], d[i−1], d[i].
            var prev2 = new int[n + 1];
            var prev = new int[n + 1];
            var curr = new int[n + 1];
            for (int j = 0; j <= n; j++)
                prev[j] = j;

            for (int i = 1; i <= m; i++)
            {
                curr[0] = i;
                for (int j = 1; j <= n; j++)
                {
                    int cost = s1[i - 1] == s2[j - 1] ? 0 : 1;
                    int v = Math.Min(Math.Min(prev[j] + 1, curr[j - 1] + 1), prev[j - 1] + cost);
                    if (i > 1 && j > 1 && s1[i - 1] == s2[j - 2] && s1[i - 2] == s2[j - 1])
                        v = Math.Min(v, prev2[j - 2] + 1);
                    curr[j] = v;
                }
                (prev2, prev, curr) = (prev, curr, prev2);
            }

            return prev[n];
        }

        /// <summary>
        /// True (unrestricted) Damerau–Levenshtein distance — Lowrance &amp; Wagner (1975,
        /// J. ACM 22(2):177–183) with unit costs: insertions, deletions, substitutions and
        /// transpositions of adjacent characters, where the transposed characters may be further
        /// separated by insertions/deletions (DL(CA, ABC) = 2 via CA → AC → ABC). A metric.
        /// Algorithm: d[i, j] = min(d[i−1, j−1] + [a_i ≠ b_j], d[i, j−1] + 1, d[i−1, j] + 1,
        /// d[k−1, l−1] + (i−k−1) + 1 + (j−l−1)), k = last row whose character equals b_j,
        /// l = last column ≤ j−1 whose character equals a_i (the "da" table over the alphabet).
        /// Case-sensitive (ordinal). O(m·n) time and space. Reference:
        /// rapidfuzz.distance.DamerauLevenshtein, jellyfish.damerau_levenshtein_distance.
        /// </summary>
        /// <exception cref="ArgumentNullException">Either string is null.</exception>
        public static int DamerauLevenshteinDistance(string s1, string s2)
        {
            if (s1 == null || s2 == null)
                throw new ArgumentNullException(s1 == null ? nameof(s1) : nameof(s2));

            int m = s1.Length;
            int n = s2.Length;
            if (m == 0) return n;
            if (n == 0) return m;

            // d is offset by one so that index 0 holds the "max distance" sentinel row/column
            // (Lowrance–Wagner's d[−1, ·] / d[·, −1]).
            int maxDist = m + n;
            var d = new int[m + 2, n + 2];
            d[0, 0] = maxDist;
            for (int i = 0; i <= m; i++)
            {
                d[i + 1, 0] = maxDist;
                d[i + 1, 1] = i;
            }
            for (int j = 0; j <= n; j++)
            {
                d[0, j + 1] = maxDist;
                d[1, j + 1] = j;
            }

            var da = new Dictionary<char, int>(); // last row (1-based) where each character occurred in s1
            for (int i = 1; i <= m; i++)
            {
                int db = 0; // last column (1-based) in this row where s2[j] == s1[i]
                for (int j = 1; j <= n; j++)
                {
                    int k = da.TryGetValue(s2[j - 1], out int row) ? row : 0;
                    int l = db;
                    int cost;
                    if (s1[i - 1] == s2[j - 1])
                    {
                        cost = 0;
                        db = j;
                    }
                    else
                    {
                        cost = 1;
                    }

                    d[i + 1, j + 1] = Math.Min(
                        Math.Min(d[i, j] + cost, d[i + 1, j] + 1),
                        Math.Min(d[i, j + 1] + 1, d[k, l] + (i - k - 1) + 1 + (j - l - 1)));
                }
                da[s1[i - 1]] = i;
            }

            return d[m + 1, n + 1];
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
            // No length-k window could be cut (e.g. k > sequence length), or no DNA k-mer lies
            // within d of any window (all windows carry ≥ d+1 non-ACGT symbols): every Count_d is
            // 0, so there is no most-frequent k-mer and the result is empty.
            foreach (var entry in MaxCountKmers(TallyNeighborhoodCounts(sequence, k, d)))
                yield return entry;
        }

        /// <summary>
        /// Count_d(Text, P) for every DNA k-mer P with a positive count (BA1I tally): distinct
        /// windows with multiplicities from the canonical (case-insensitive) k-mer counter, so each
        /// distinct window's d-neighborhood is generated once and weighted by its multiplicity:
        /// Count_d(P) = Σ_w mult(w)·[HD(P, w) ≤ d]. Empty for a null/empty sequence.
        /// </summary>
        private static Dictionary<string, int> TallyNeighborhoodCounts(string sequence, int k, int d)
        {
            var counts = new Dictionary<string, int>();
            if (string.IsNullOrEmpty(sequence))
                return counts;

            var windowCounts = sequence.AsSpan().CountKmersSpan(k);
            var buffer = new char[k];

            foreach (var (window, multiplicity) in windowCounts)
            {
                foreach (string neighbor in GenerateNeighbors(window, d, buffer))
                {
                    counts[neighbor] = counts.TryGetValue(neighbor, out int c) ? c + multiplicity : multiplicity;
                }
            }

            return counts;
        }

        /// <summary>All entries achieving the maximum count (nothing for an empty tally).</summary>
        private static IEnumerable<(string Kmer, int Count)> MaxCountKmers(Dictionary<string, int> counts)
        {
            if (counts.Count == 0)
                yield break;

            int maxCount = counts.Values.Max();
            foreach (var kvp in counts.Where(kvp => kvp.Value == maxCount))
                yield return (kvp.Key, kvp.Value);
        }

        /// <summary>
        /// Frequent Words with Mismatches and Reverse Complements (Compeau &amp; Pevzner,
        /// Bioinformatics Algorithms ch.1; ROSALIND BA1J): all DNA k-mers P maximising
        /// Count_d(Text, P) + Count_d(Text, ReverseComplement(P)). Reuses the BA1I tally
        /// (<see cref="FindFrequentKmersWithMismatches"/>: same case handling, same ACGT-only
        /// neighbourhoods) and the canonical reverse complement
        /// (<see cref="DnaSequence.GetReverseComplementString"/>). A reverse-complement palindrome
        /// P = rc(P) scores 2·Count_d(Text, P), as in the problem's definition. The result set is
        /// closed under reverse complement. Result order is unspecified.
        /// Sample: ACGTTGCATGTCGCATGATGCATGAGAGCT, k = 4, d = 1 → {ATGT, ACAT} (score 9).
        /// </summary>
        /// <param name="sequence">The sequence to analyze.</param>
        /// <param name="k">K-mer length.</param>
        /// <param name="d">Maximum mismatches.</param>
        /// <returns>All maximising k-mers (ties included) with Count_d(P) + Count_d(rc(P)).</returns>
        public static IEnumerable<(string Kmer, int Count)> FindFrequentKmersWithMismatchesAndReverseComplements(
            string sequence, int k, int d)
        {
            if (k <= 0)
                throw new ArgumentOutOfRangeException(nameof(k), "K must be positive.");
            if (d < 0)
                throw new ArgumentOutOfRangeException(nameof(d), "D cannot be negative.");
            return FindFrequentKmersWithMismatchesAndReverseComplementsCore(sequence, k, d);
        }

        private static IEnumerable<(string Kmer, int Count)> FindFrequentKmersWithMismatchesAndReverseComplementsCore(
            string sequence, int k, int d)
        {
            var counts = TallyNeighborhoodCounts(sequence, k, d);

            // A k-mer has a positive combined score iff it or its reverse complement has a positive
            // Count_d, so scoring every tallied k-mer and its reverse complement covers all candidates.
            var combined = new Dictionary<string, int>();
            foreach (string kmer in counts.Keys)
            {
                string rc = DnaSequence.GetReverseComplementString(kmer);
                int score = counts[kmer] + counts.GetValueOrDefault(rc);
                combined[kmer] = score;
                combined[rc] = score;
            }

            foreach (var entry in MaxCountKmers(combined))
                yield return entry;
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

        /// <summary>
        /// Optimal alignment of the pattern (query) against <see cref="MatchedSequence"/> (target)
        /// for edit-distance results (<see cref="ApproximateMatcher.FindWithEdits(string, string, int)"/>);
        /// null for Hamming results, whose alignment is the ungapped diagonal.
        /// </summary>
        public EditAlignment? Alignment { get; init; }
    }

    /// <summary>
    /// A unit-cost (Levenshtein) alignment path between a query and a target, in edlib's
    /// convention (Šošić &amp; Šikić 2017, edlib.h): '=' match, 'X' mismatch, 'I' insertion to the
    /// query (a query character with no target counterpart, EDLIB_EDOP_INSERT), 'D' deletion from
    /// the query (a target character with no query counterpart, EDLIB_EDOP_DELETE).
    /// </summary>
    public sealed record EditAlignment
    {
        internal EditAlignment(int distance, string operations, string query, string target)
        {
            Distance = distance;
            Operations = operations;

            var alignedQuery = new System.Text.StringBuilder(operations.Length);
            var alignedTarget = new System.Text.StringBuilder(operations.Length);
            var substitutions = new List<int>();
            int q = 0, t = 0;
            foreach (char op in operations)
            {
                switch (op)
                {
                    case 'I':
                        alignedQuery.Append(query[q++]);
                        alignedTarget.Append('-');
                        HasIndels = true;
                        break;
                    case 'D':
                        alignedQuery.Append('-');
                        alignedTarget.Append(target[t++]);
                        HasIndels = true;
                        break;
                    default:
                        if (op == 'X')
                            substitutions.Add(q);
                        alignedQuery.Append(query[q++]);
                        alignedTarget.Append(target[t++]);
                        break;
                }
            }

            AlignedQuery = alignedQuery.ToString();
            AlignedTarget = alignedTarget.ToString();
            SubstitutionPositions = substitutions.AsReadOnly();
            Cigar = RunLength(operations, standard: false);
            StandardCigar = RunLength(operations, standard: true);
        }

        /// <summary>Number of non-'=' operations (the edit distance).</summary>
        public int Distance { get; }

        /// <summary>One operation per alignment column: '=', 'X', 'I' or 'D'.</summary>
        public string Operations { get; }

        /// <summary>edlib EXTENDED CIGAR (run-length '=', 'X', 'I', 'D'), e.g. "2=1D2=".</summary>
        public string Cigar { get; }

        /// <summary>edlib STANDARD CIGAR ('M' for both '=' and 'X'), e.g. "2M1D2M".</summary>
        public string StandardCigar { get; }

        /// <summary>Query with '-' in 'D' columns.</summary>
        public string AlignedQuery { get; }

        /// <summary>Target with '-' in 'I' columns.</summary>
        public string AlignedTarget { get; }

        /// <summary>Query-relative (0-based) indices of the substituted ('X') query characters.</summary>
        public IReadOnlyList<int> SubstitutionPositions { get; }

        /// <summary>True when the path contains an insertion or deletion.</summary>
        public bool HasIndels { get; }

        /// <summary>Value equality over the path and the aligned strings.</summary>
        public bool Equals(EditAlignment? other) =>
            other is not null && Distance == other.Distance && Operations == other.Operations
            && AlignedQuery == other.AlignedQuery && AlignedTarget == other.AlignedTarget;

        /// <inheritdoc />
        public override int GetHashCode() => HashCode.Combine(Distance, Operations, AlignedQuery, AlignedTarget);

        /// <summary>edlibAlignmentToCigar: counts of consecutive identical move characters.</summary>
        private static string RunLength(string operations, bool standard)
        {
            var sb = new System.Text.StringBuilder();
            int i = 0;
            while (i < operations.Length)
            {
                char c = Map(operations[i]);
                int run = 1;
                while (i + run < operations.Length && Map(operations[i + run]) == c)
                    run++;
                sb.Append(run).Append(c);
                i += run;
            }
            return sb.ToString();

            char Map(char op) => standard && (op == '=' || op == 'X') ? 'M' : op;
        }
    }
}

