using System.Numerics;

namespace Seqeron.Genomics.Alignment
{
    /// <summary>
    /// Weighted Damerau distances (Lowrance &amp; Wagner 1975) and transposition-aware tracebacks.
    /// </summary>
    public static partial class ApproximateMatcher
    {
        /// <summary>
        /// Weighted optimal string alignment (restricted Damerau–Levenshtein) distance:
        /// d[i, 0] = i·Deletion, d[0, j] = j·Insertion,
        /// d[i, j] = min(d[i−1, j] + Deletion, d[i, j−1] + Insertion,
        /// d[i−1, j−1] + [a_i ≠ b_j]·Substitution,
        /// d[i−2, j−2] + Transposition if a_i = b_{j−1} and a_{i−1} = b_j).
        /// Costs follow <see cref="EditCosts"/> (s1 → s2: an insertion adds a character of s2, a deletion
        /// removes a character of s1); no relation between the costs is required (the value is defined
        /// by the recurrence). Equal to R stringdist <c>method = "osa"</c> with
        /// <c>weight = c(Insertion, Deletion, Substitution, Transposition)</c> — stringdist's
        /// "d" weight is charged for characters of its second argument, i.e. our insertion.
        /// Unit costs return exactly <see cref="OptimalStringAlignmentDistance(string, string)"/>;
        /// uniform costs c return c times it. Otherwise O(m·n) time, O(n) space. Case-sensitive.
        /// </summary>
        /// <exception cref="ArgumentNullException">Either string is null.</exception>
        /// <exception cref="ArgumentOutOfRangeException">A cost is negative.</exception>
        /// <exception cref="OverflowException">The distance exceeds <see cref="int.MaxValue"/>.</exception>
        public static int OptimalStringAlignmentDistance(string s1, string s2, DamerauCosts costs)
        {
            if (s1 == null || s2 == null)
                throw new ArgumentNullException(s1 == null ? nameof(s1) : nameof(s2));
            ValidateDamerauCosts(costs, unrestricted: false);

            if (costs.IsUniform)
                return checked(costs.Insertion * OptimalStringAlignmentDistance(s1, s2));

            return checked((int)OsaDistanceRows<long>(
                s1, s2, costs.Insertion, costs.Deletion, costs.Substitution, costs.Transposition));
        }

        /// <summary>
        /// The single OSA (restricted Damerau–Levenshtein) DP, three rolling rows over <paramref name="s2"/>
        /// (d[i−2], d[i−1], d[i]); unit costs in 32-bit cells (<see cref="OptimalStringAlignmentDistance(string, string)"/>),
        /// weighted in 64-bit cells. d[i, j] = min(min(d[i−1, j] + del, d[i, j−1] + ins), d[i−1, j−1] + [a_i ≠ b_j]·sub),
        /// then min with d[i−2, j−2] + tr when a_i = b_{j−1} and a_{i−1} = b_j.
        /// </summary>
        private static T OsaDistanceRows<T>(string s1, string s2, T ins, T del, T sub, T tr)
            where T : struct, INumber<T>
        {
            int m = s1.Length;
            int n = s2.Length;
            var prev2 = new T[n + 1];
            var prev = new T[n + 1];
            var curr = new T[n + 1];
            for (int j = 0; j <= n; j++)
                prev[j] = T.CreateTruncating(j) * ins;

            for (int i = 1; i <= m; i++)
            {
                curr[0] = T.CreateTruncating(i) * del;
                for (int j = 1; j <= n; j++)
                {
                    T v = T.Min(T.Min(prev[j] + del, curr[j - 1] + ins),
                        prev[j - 1] + (s1[i - 1] == s2[j - 1] ? T.Zero : sub));
                    if (i > 1 && j > 1 && s1[i - 1] == s2[j - 2] && s1[i - 2] == s2[j - 1])
                        v = T.Min(v, prev2[j - 2] + tr);
                    curr[j] = v;
                }
                (prev2, prev, curr) = (prev, curr, prev2);
            }

            return prev[n];
        }

        /// <summary>
        /// Weighted optimal string alignment distance with the transposition cost given separately;
        /// see <see cref="OptimalStringAlignmentDistance(string, string, DamerauCosts)"/>.
        /// </summary>
        /// <exception cref="ArgumentNullException">Either string is null.</exception>
        /// <exception cref="ArgumentOutOfRangeException">A cost is negative.</exception>
        /// <exception cref="OverflowException">The distance exceeds <see cref="int.MaxValue"/>.</exception>
        public static int OptimalStringAlignmentDistance(string s1, string s2, EditCosts costs, int transpositionCost) =>
            OptimalStringAlignmentDistance(s1, s2, new DamerauCosts(costs, transpositionCost));

        /// <summary>
        /// Weighted (unrestricted) Damerau–Levenshtein distance — the extended string-to-string
        /// correction problem of Lowrance &amp; Wagner (1975, J. ACM 22(2):177–183, "An extension of the
        /// string-to-string correction problem"): the minimum total cost of insertions (W_I),
        /// deletions (W_D), substitutions (W_C) and swaps of two adjacent characters (W_S) turning s1
        /// into s2, with no restriction on editing a substring twice. Recurrence (Lowrance–Wagner):
        /// H[i, 0] = i·W_D, H[0, j] = j·W_I,
        /// H[i, j] = min(H[i−1, j−1] + [a_i ≠ b_j]·W_C, H[i, j−1] + W_I, H[i−1, j] + W_D,
        /// H[k−1, l−1] + (i−k−1)·W_D + W_S + (j−l−1)·W_I), where k is the last row &lt; i with
        /// a_k = b_j and l the last column &lt; j with b_l = a_i (the term is absent when either does not
        /// exist). Lowrance &amp; Wagner prove the recurrence is exact when 2·W_S ≥ W_I + W_D; for
        /// cheaper swaps the problem is NP-complete in general (Wagner 1975) and the recurrence is no
        /// longer the minimum (e.g. W_I = W_D = 2, W_S = 1: ABC → BCA costs 2 by two swaps, the
        /// recurrence gives 4), so such costs are rejected.
        /// </summary>
        /// <remarks>
        /// Unit costs return exactly <see cref="DamerauLevenshteinDistance(string, string)"/> (linear-space
        /// Zhao–Sahni engine); uniform costs c return c times it. Other costs use the Lowrance–Wagner full
        /// matrix: O(m·n) time and space. Verified against an exhaustive Dijkstra search over all edit
        /// sequences on small strings. R stringdist <c>method = "dl"</c> agrees only when the deletion,
        /// insertion and transposition weights are equal: it charges the characters deleted / inserted
        /// between the transposed pair at the transposition weight instead of W_D / W_I (and can then
        /// report less than the true minimum, e.g. CAAAACA → CA with W_I = 1, W_D = 7, W_C = W_S = 5:
        /// stringdist 32, true minimum 35 = five deletions). Case-sensitive.
        /// </remarks>
        /// <exception cref="ArgumentNullException">Either string is null.</exception>
        /// <exception cref="ArgumentOutOfRangeException">A cost is negative.</exception>
        /// <exception cref="ArgumentException">2·Transposition &lt; Insertion + Deletion (Lowrance–Wagner condition violated).</exception>
        /// <exception cref="OverflowException">The distance exceeds <see cref="int.MaxValue"/>.</exception>
        public static int DamerauLevenshteinDistance(string s1, string s2, DamerauCosts costs)
        {
            if (s1 == null || s2 == null)
                throw new ArgumentNullException(s1 == null ? nameof(s1) : nameof(s2));
            ValidateDamerauCosts(costs, unrestricted: true);

            if (costs.IsUniform)
                return checked(costs.Insertion * DamerauLevenshteinDistance(s1, s2));

            var h = FillDamerauMatrix(s1, s2, costs, restricted: false);
            return checked((int)h[s1.Length, s2.Length]);
        }

        /// <summary>
        /// Weighted unrestricted Damerau–Levenshtein distance with the transposition cost given
        /// separately; see <see cref="DamerauLevenshteinDistance(string, string, DamerauCosts)"/>.
        /// </summary>
        /// <exception cref="ArgumentNullException">Either string is null.</exception>
        /// <exception cref="ArgumentOutOfRangeException">A cost is negative.</exception>
        /// <exception cref="ArgumentException">2·transpositionCost &lt; Insertion + Deletion.</exception>
        /// <exception cref="OverflowException">The distance exceeds <see cref="int.MaxValue"/>.</exception>
        public static int DamerauLevenshteinDistance(string s1, string s2, EditCosts costs, int transpositionCost) =>
            DamerauLevenshteinDistance(s1, s2, new DamerauCosts(costs, transpositionCost));

        /// <summary>
        /// Optimal unit-cost Damerau–Levenshtein edit script (Lowrance–Wagner trace) turning s1 into s2;
        /// see <see cref="GetDamerauLevenshteinAlignment(string, string, DamerauCosts)"/>.
        /// <see cref="DamerauAlignment.Distance"/> equals <see cref="DamerauLevenshteinDistance(string, string)"/>.
        /// </summary>
        /// <exception cref="ArgumentNullException">Either string is null.</exception>
        public static DamerauAlignment GetDamerauLevenshteinAlignment(string s1, string s2) =>
            GetDamerauLevenshteinAlignment(s1, s2, DamerauCosts.Unit);

        /// <summary>
        /// Optimal weighted Damerau–Levenshtein edit script turning s1 into s2: the Lowrance–Wagner (1975)
        /// matrix of <see cref="DamerauLevenshteinDistance(string, string, DamerauCosts)"/> followed by a
        /// traceback from (m, n). Each step is a match, a substitution, an insertion (a character of s2),
        /// a deletion (a character of s1) or a Lowrance–Wagner transposition block: a_k … a_i with
        /// a_k = b_j and a_i = b_l becomes b_l … b_j by deleting the i−k−1 characters between the pair,
        /// swapping it, and inserting the j−l−1 characters b_{l+1..j−1} between them, at cost
        /// (i−k−1)·Deletion + Transposition + (j−l−1)·Insertion. Among co-optimal steps the traceback
        /// prefers the diagonal (match / substitution), then the transposition, then the deletion, then
        /// the insertion. <see cref="DamerauAlignment.Distance"/> is the summed cost of the operations
        /// and equals the distance; applying the operations to s1 yields s2.
        /// O(m·n) time and space. Same cost validation (2·Transposition ≥ Insertion + Deletion).
        /// </summary>
        /// <exception cref="ArgumentNullException">Either string is null.</exception>
        /// <exception cref="ArgumentOutOfRangeException">A cost is negative.</exception>
        /// <exception cref="ArgumentException">2·Transposition &lt; Insertion + Deletion.</exception>
        /// <exception cref="OverflowException">The distance exceeds <see cref="int.MaxValue"/>.</exception>
        public static DamerauAlignment GetDamerauLevenshteinAlignment(string s1, string s2, DamerauCosts costs)
        {
            if (s1 == null || s2 == null)
                throw new ArgumentNullException(s1 == null ? nameof(s1) : nameof(s2));
            ValidateDamerauCosts(costs, unrestricted: true);

            var h = FillDamerauMatrix(s1, s2, costs, restricted: false);
            return TracebackDamerau(s1, s2, costs, h, restricted: false);
        }

        /// <summary>
        /// Optimal unit-cost optimal-string-alignment edit script; see
        /// <see cref="GetOptimalStringAlignment(string, string, DamerauCosts)"/>.
        /// <see cref="DamerauAlignment.Distance"/> equals <see cref="OptimalStringAlignmentDistance(string, string)"/>.
        /// </summary>
        /// <exception cref="ArgumentNullException">Either string is null.</exception>
        public static DamerauAlignment GetOptimalStringAlignment(string s1, string s2) =>
            GetOptimalStringAlignment(s1, s2, DamerauCosts.Unit);

        /// <summary>
        /// Optimal edit script under the (weighted) optimal string alignment recurrence of
        /// <see cref="OptimalStringAlignmentDistance(string, string, DamerauCosts)"/>: matches,
        /// substitutions, insertions, deletions and swaps of two adjacent characters
        /// (a_{i−1} a_i = b_j b_{j−1}; a transposition block of length 2 on both sides, nothing
        /// inserted or deleted between). Same tie-break and result contract as
        /// <see cref="GetDamerauLevenshteinAlignment(string, string, DamerauCosts)"/>. O(m·n) time and space.
        /// </summary>
        /// <exception cref="ArgumentNullException">Either string is null.</exception>
        /// <exception cref="ArgumentOutOfRangeException">A cost is negative.</exception>
        /// <exception cref="OverflowException">The distance exceeds <see cref="int.MaxValue"/>.</exception>
        public static DamerauAlignment GetOptimalStringAlignment(string s1, string s2, DamerauCosts costs)
        {
            if (s1 == null || s2 == null)
                throw new ArgumentNullException(s1 == null ? nameof(s1) : nameof(s2));
            ValidateDamerauCosts(costs, unrestricted: false);

            var h = FillDamerauMatrix(s1, s2, costs, restricted: true);
            return TracebackDamerau(s1, s2, costs, h, restricted: true);
        }

        private static void ValidateDamerauCosts(DamerauCosts costs, bool unrestricted)
        {
            ValidateCosts(new EditCosts(costs.Insertion, costs.Deletion, costs.Substitution));
            if (costs.Transposition < 0)
                throw new ArgumentOutOfRangeException(nameof(costs), "Transposition cost cannot be negative.");
            if (unrestricted && 2L * costs.Transposition < (long)costs.Insertion + costs.Deletion)
                throw new ArgumentException(
                    "Unrestricted Damerau–Levenshtein requires 2·Transposition ≥ Insertion + Deletion " +
                    "(Lowrance & Wagner 1975); use the optimal string alignment distance for cheaper transpositions.",
                    nameof(costs));
        }

        /// <summary>
        /// Full (m + 1) × (n + 1) Lowrance–Wagner matrix (restricted = false) or optimal string
        /// alignment matrix (restricted = true) under weighted costs.
        /// </summary>
        private static long[,] FillDamerauMatrix(string s1, string s2, DamerauCosts costs, bool restricted)
        {
            int m = s1.Length;
            int n = s2.Length;
            long ins = costs.Insertion, del = costs.Deletion, sub = costs.Substitution, tr = costs.Transposition;
            var h = new long[m + 1, n + 1];
            for (int i = 0; i <= m; i++)
                h[i, 0] = i * del;
            for (int j = 0; j <= n; j++)
                h[0, j] = j * ins;

            var da = new Dictionary<char, int>(); // last row (1-based) where each character occurred in s1
            for (int i = 1; i <= m; i++)
            {
                char a = s1[i - 1];
                int db = 0; // last column (1-based) < j in this row with b_l = a_i
                for (int j = 1; j <= n; j++)
                {
                    char b = s2[j - 1];
                    long v = Math.Min(Math.Min(h[i - 1, j - 1] + (a == b ? 0 : sub), h[i, j - 1] + ins), h[i - 1, j] + del);
                    if (restricted)
                    {
                        if (i > 1 && j > 1 && a == s2[j - 2] && s1[i - 2] == b)
                            v = Math.Min(v, h[i - 2, j - 2] + tr);
                    }
                    else
                    {
                        int k = da.TryGetValue(b, out int row) ? row : 0;
                        int l = db;
                        if (k >= 1 && l >= 1)
                            v = Math.Min(v, h[k - 1, l - 1] + (i - k - 1) * del + tr + (j - l - 1) * ins);
                    }
                    if (a == b)
                        db = j;
                    h[i, j] = v;
                }
                da[a] = i;
            }

            return h;
        }

        private static DamerauAlignment TracebackDamerau(string s1, string s2, DamerauCosts costs, long[,] h, bool restricted)
        {
            long ins = costs.Insertion, del = costs.Deletion, sub = costs.Substitution, tr = costs.Transposition;
            Dictionary<char, List<int>>? rows1 = null, cols2 = null;
            if (!restricted)
            {
                rows1 = Positions(s1);
                cols2 = Positions(s2);
            }

            var ops = new List<DamerauEditOperation>();
            int i = s1.Length;
            int j = s2.Length;
            while (i > 0 || j > 0)
            {
                long v = h[i, j];
                if (i > 0 && j > 0)
                {
                    bool same = s1[i - 1] == s2[j - 1];
                    long c = same ? 0 : sub;
                    if (h[i - 1, j - 1] + c == v)
                    {
                        ops.Add(new DamerauEditOperation(
                            same ? DamerauOperationKind.Match : DamerauOperationKind.Substitution,
                            i - 1, 1, j - 1, 1, checked((int)c)));
                        i--;
                        j--;
                        continue;
                    }

                    int k = 0, l = 0;
                    if (restricted)
                    {
                        if (i > 1 && j > 1 && s1[i - 1] == s2[j - 2] && s1[i - 2] == s2[j - 1])
                        {
                            k = i - 1;
                            l = j - 1;
                        }
                    }
                    else
                    {
                        k = LastBefore(rows1!, s2[j - 1], i);
                        l = LastBefore(cols2!, s1[i - 1], j);
                    }

                    if (k >= 1 && l >= 1)
                    {
                        long tc = (i - k - 1) * del + tr + (j - l - 1) * ins;
                        if (h[k - 1, l - 1] + tc == v)
                        {
                            ops.Add(new DamerauEditOperation(
                                DamerauOperationKind.Transposition, k - 1, i - k + 1, l - 1, j - l + 1, checked((int)tc)));
                            i = k - 1;
                            j = l - 1;
                            continue;
                        }
                    }
                }

                if (i > 0 && h[i - 1, j] + del == v)
                {
                    ops.Add(new DamerauEditOperation(DamerauOperationKind.Deletion, i - 1, 1, j, 0, checked((int)del)));
                    i--;
                }
                else
                {
                    ops.Add(new DamerauEditOperation(DamerauOperationKind.Insertion, i, 0, j - 1, 1, checked((int)ins)));
                    j--;
                }
            }

            ops.Reverse();
            return new DamerauAlignment(checked((int)h[s1.Length, s2.Length]), ops);

            static Dictionary<char, List<int>> Positions(string s)
            {
                var map = new Dictionary<char, List<int>>();
                for (int p = 0; p < s.Length; p++)
                {
                    if (!map.TryGetValue(s[p], out var list))
                        map[s[p]] = list = new List<int>();
                    list.Add(p + 1);
                }
                return map;
            }

            // Largest 1-based position < limit holding c, or 0.
            static int LastBefore(Dictionary<char, List<int>> map, char c, int limit)
            {
                if (!map.TryGetValue(c, out var list))
                    return 0;
                int idx = list.BinarySearch(limit);
                idx = idx >= 0 ? idx - 1 : ~idx - 1;
                return idx >= 0 ? list[idx] : 0;
            }
        }
    }

    /// <summary>
    /// Non-negative additive costs for the weighted Damerau distances of <see cref="ApproximateMatcher"/>
    /// (Lowrance &amp; Wagner 1975: W_I, W_D, W_C, W_S), in the <see cref="EditCosts"/> convention
    /// (s1 → s2: an insertion adds a character of s2, a deletion removes a character of s1) plus the
    /// cost of swapping two adjacent characters.
    /// </summary>
    /// <param name="Insertion">Cost of inserting one character of s2 (W_I).</param>
    /// <param name="Deletion">Cost of deleting one character of s1 (W_D).</param>
    /// <param name="Substitution">Cost of replacing one character (W_C).</param>
    /// <param name="Transposition">Cost of swapping two adjacent characters (W_S).</param>
    public readonly record struct DamerauCosts(int Insertion, int Deletion, int Substitution, int Transposition)
    {
        /// <summary>Combines Levenshtein costs with a transposition cost.</summary>
        public DamerauCosts(EditCosts costs, int transposition)
            : this(costs.Insertion, costs.Deletion, costs.Substitution, transposition)
        {
        }

        /// <summary>Unit costs (1, 1, 1, 1): plain Damerau–Levenshtein / optimal string alignment distance.</summary>
        public static DamerauCosts Unit => new(1, 1, 1, 1);

        /// <summary>True when all four costs are equal.</summary>
        internal bool IsUniform =>
            Insertion == Deletion && Insertion == Substitution && Insertion == Transposition;
    }

    /// <summary>Kind of a <see cref="DamerauEditOperation"/>.</summary>
    public enum DamerauOperationKind
    {
        /// <summary>a_i = b_j kept (cost 0).</summary>
        Match,

        /// <summary>a_i replaced by b_j (substitution cost).</summary>
        Substitution,

        /// <summary>b_j inserted (insertion cost); SourceLength 0.</summary>
        Insertion,

        /// <summary>a_i deleted (deletion cost); TargetLength 0.</summary>
        Deletion,

        /// <summary>
        /// Lowrance–Wagner transposition block: source a_k … a_i (a_k = b_j, a_i = b_l) becomes target
        /// b_l … b_j — the SourceLength − 2 characters between the pair are deleted, the pair is
        /// swapped, and the TargetLength − 2 characters b_{l+1..j−1} are inserted between them.
        /// The plain adjacent swap has SourceLength = TargetLength = 2.
        /// </summary>
        Transposition,
    }

    /// <summary>
    /// One step of a <see cref="DamerauAlignment"/>: it consumes source (s1) characters
    /// [SourcePosition, SourcePosition + SourceLength) and produces target (s2) characters
    /// [TargetPosition, TargetPosition + TargetLength) (0-based).
    /// </summary>
    /// <param name="Kind">Operation kind.</param>
    /// <param name="SourcePosition">0-based start in s1 (for an insertion: the s1 index before which it occurs).</param>
    /// <param name="SourceLength">Number of s1 characters consumed.</param>
    /// <param name="TargetPosition">0-based start in s2 (for a deletion: the s2 index before which it occurs).</param>
    /// <param name="TargetLength">Number of s2 characters produced.</param>
    /// <param name="Cost">Cost of the step under the costs used.</param>
    public readonly record struct DamerauEditOperation(
        DamerauOperationKind Kind, int SourcePosition, int SourceLength, int TargetPosition, int TargetLength, int Cost);

    /// <summary>
    /// A transposition-aware edit script from s1 to s2 (Lowrance–Wagner trace), returned by
    /// <see cref="ApproximateMatcher.GetDamerauLevenshteinAlignment(string, string, DamerauCosts)"/> and
    /// <see cref="ApproximateMatcher.GetOptimalStringAlignment(string, string, DamerauCosts)"/>.
    /// </summary>
    public sealed record DamerauAlignment
    {
        internal DamerauAlignment(int distance, List<DamerauEditOperation> operations)
        {
            Distance = distance;
            Operations = operations.AsReadOnly();

            var sb = new System.Text.StringBuilder();
            foreach (var op in operations)
            {
                switch (op.Kind)
                {
                    case DamerauOperationKind.Match: sb.Append('='); break;
                    case DamerauOperationKind.Substitution: sb.Append('X'); break;
                    case DamerauOperationKind.Deletion: sb.Append('I'); break;
                    case DamerauOperationKind.Insertion: sb.Append('D'); break;
                    default:
                        TranspositionCount++;
                        sb.Append('T').Append('i', op.SourceLength - 2).Append('d', op.TargetLength - 2);
                        break;
                }
            }
            Script = sb.ToString();
        }

        /// <summary>Summed cost of <see cref="Operations"/> (the distance under the costs used).</summary>
        public int Distance { get; }

        /// <summary>The edit steps, left to right.</summary>
        public IReadOnlyList<DamerauEditOperation> Operations { get; }

        /// <summary>
        /// Compact script in the <see cref="EditAlignment.Operations"/> (edlib) letters: '=' match,
        /// 'X' substitution, 'I' an s1 character absent from s2 (deletion), 'D' an s2 character absent
        /// from s1 (insertion); a transposition block is 'T' followed by one 'i' per character deleted
        /// and one 'd' per character inserted between the swapped pair (an adjacent swap is "T").
        /// </summary>
        public string Script { get; }

        /// <summary>Number of transposition blocks.</summary>
        public int TranspositionCount { get; }

        /// <summary>Value equality over the distance and the operations.</summary>
        public bool Equals(DamerauAlignment? other) =>
            other is not null && Distance == other.Distance && Operations.SequenceEqual(other.Operations);

        /// <inheritdoc />
        public override int GetHashCode() => HashCode.Combine(Distance, Script);
    }
}
