namespace Seqeron.Genomics.Alignment;

/// <summary>Options of <see cref="NucmerClusters.Cluster"/>; the defaults are the ones <c>nucmer</c> hands to <c>mgaps</c> (<c>-g 90 -c 65 -D 5 -d 0.12</c>).</summary>
public sealed record NucmerClusterOptions
{
    /// <summary>Largest separation between two matches of a cluster (<c>mgaps -s</c> / <c>nucmer -g</c>, default 90).</summary>
    public long MaxSeparation { get; init; } = 90;

    /// <summary>Smallest cluster length that is reported (<c>mgaps -l</c> / <c>nucmer -c</c>, default 65).</summary>
    public long MinClusterLength { get; init; } = 65;

    /// <summary>Fixed diagonal difference that joins two matches (<c>mgaps -d</c> / <c>nucmer -D</c>, default 5).</summary>
    public int DiagonalDifference { get; init; } = 5;

    /// <summary>Fraction of the separation that is added to the allowed diagonal difference (<c>mgaps -f</c> / <c>nucmer -d</c>, default 0.12).</summary>
    public double DiagonalFraction { get; init; } = 0.12;

    /// <summary>Cluster length = last end - first start in the first sequence instead of the sum of the match lengths (<c>mgaps -e</c>).</summary>
    public bool UseExtent { get; init; }
}

/// <summary>One match of a <see cref="NucmerCluster"/> as <c>mgaps</c> prints it: the start in both sequences and the length after the overlap with the previous match has been cut off.</summary>
/// <param name="Start1">Start in the first sequence (the coordinate system of the input).</param>
/// <param name="Start2">Start in the second sequence.</param>
/// <param name="Length">Length left after the cut.</param>
/// <param name="Overlap">How far the match was shortened at its start because it overlaps the previous match of the cluster (0 = untouched).</param>
public readonly record struct NucmerClusterMatch(long Start1, long Start2, long Length, long Overlap);

/// <summary>A cluster of exact matches: a chain of matches that lie within the allowed diagonal difference of each other and are ordered along both sequences.</summary>
/// <param name="Matches">The matches in order, adjusted as the program prints them.</param>
public sealed record NucmerCluster(IReadOnlyList<NucmerClusterMatch> Matches);

/// <summary>
/// Clustering of exact matches into longer inexact matches by diagonal and separation: the <c>mgaps</c> program of MUMmer 3 (Delcher et al. 2002, Nucleic Acids Res 30:2478; the program that <c>nucmer</c> runs
/// between <c>mummer</c> and <c>postnuc</c>). Port of <c>mgaps</c> (<c>src/tigr/mgaps.cc</c>, <c>include/mummer/mgaps.hh</c> of the MUMmer sources, the algorithm unchanged since 3.x): matches are
/// sorted by their start in the second sequence, matches inside a repeat are removed and overlapping matches of one diagonal are merged (<c>Filter_Matches</c>), neighbours whose
/// diagonals differ by at most <c>max(d, f · separation)</c> within the largest separation are joined by union-find, and every component is cut into chains: repeatedly the chain of
/// the largest score (match length minus overlap minus diagonal difference, a dynamic program over the component) is taken out and reported when its length reaches the minimum.
/// </summary>
public static class NucmerClusters
{
    private sealed class MatchT
    {
        public long Start1, Start2, Len;
        public long SimpleScore, SimpleAdj;
        public int SimpleFrom;
        public int ClusterId;
        public bool Good = true;
        public bool Tentative;
    }

    private sealed class UnionFind
    {
        private int[] _uf = [];

        public void Reset(int size)
        {
            _uf = new int[size + 1];
            Array.Fill(_uf, -1);
        }

        public int Find(int a)
        {
            if (_uf[a] < 0)
                return a;
            int i = a;
            while (_uf[i] > 0)
                i = _uf[i];
            for (int j = a, k; _uf[j] != i; j = k)
            {
                k = _uf[j];
                _uf[j] = i;
            }

            return i;
        }

        public void Union(int a, int b)
        {
            if (a == b)
                return;
            if (_uf[a] < _uf[b])
            {
                _uf[a] += _uf[b];
                _uf[b] = a;
            }
            else
            {
                _uf[b] += _uf[a];
                _uf[a] = b;
            }
        }
    }

    /// <summary>
    /// Clusters the matches (<c>mgaps</c> for one query sequence). The matches may be given in any order; a coordinate origin of 0 or 1 gives the same clusters. Two matches with the same
    /// start in both sequences are not supported (the program's sort leaves their order unspecified).
    /// </summary>
    /// <param name="matches">The exact matches: start in the first (reference) and second (query) sequence and the length.</param>
    /// <param name="options">Clustering options (null = nucmer defaults).</param>
    /// <returns>The reported clusters in the order the program prints them.</returns>
    /// <exception cref="ArgumentNullException"><paramref name="matches"/> is null.</exception>
    /// <exception cref="ArgumentException">A match has a non-positive length or two matches coincide.</exception>
    public static IReadOnlyList<NucmerCluster> Cluster(IEnumerable<(long Start1, long Start2, long Length)> matches, NucmerClusterOptions? options = null)
    {
        ArgumentNullException.ThrowIfNull(matches);
        var o = options ?? new NucmerClusterOptions();
        var raw = matches.ToList();
        if (raw.Any(m => m.Length <= 0))
            throw new ArgumentException("A match must have a positive length.", nameof(matches));
        if (raw.Select(m => (m.Start1, m.Start2)).Distinct().Count() != raw.Count)
            throw new ArgumentException("Two matches have the same start in both sequences.", nameof(matches));

        var clusters = new List<NucmerCluster>();
        int n = raw.Count;
        if (n == 0)
            return clusters;

        // 1-based array as in the program
        var a = new MatchT[n + 1];
        for (int i = 0; i < n; i++)
            a[i + 1] = new MatchT { Start1 = raw[i].Start1, Start2 = raw[i].Start2, Len = raw[i].Length };
        var uf = new UnionFind();
        uf.Reset(n);

        Array.Sort(a, 1, n, Comparer<MatchT>.Create(ByStart2));
        n = FilterMatches(a, 1, n);

        for (int i = 1; i < n; i++)
        {
            long iEnd = a[i].Start2 + a[i].Len;
            long iDiag = a[i].Start2 - a[i].Start1;
            for (int j = i + 1; j <= n; j++)
            {
                long sep = a[j].Start2 - iEnd;
                if (sep > o.MaxSeparation)
                    break;
                long diagDiff = Math.Abs(a[j].Start2 - a[j].Start1 - iDiag);
                if (diagDiff <= Math.Max(o.DiagonalDifference, (int)(o.DiagonalFraction * sep)))
                    uf.Union(uf.Find(i), uf.Find(j));
            }
        }

        for (int i = 1; i <= n; i++)
        {
            a[i].ClusterId = uf.Find(i);
            a[i].Good = false;
        }

        Array.Sort(a, 1, n, Comparer<MatchT>.Create((x, y) => x.ClusterId != y.ClusterId ? x.ClusterId.CompareTo(y.ClusterId) : ByStart2(x, y)));

        for (int i = 1, size; i <= n; i += size)
        {
            int j = i + 1;
            while (j <= n && a[i].ClusterId == a[j].ClusterId)
                j++;
            size = j - i;
            ProcessCluster(a, i, size, o, clusters);
        }

        return clusters;
    }

    private static int ByStart2(MatchT x, MatchT y) => x.Start2 != y.Start2 ? x.Start2.CompareTo(y.Start2) : x.Start1.CompareTo(y.Start1);

    // Filter_Matches on A[first .. first + count - 1]; returns the new count
    private static int FilterMatches(MatchT[] arr, int first, int count)
    {
        for (int i = 0; i < count - 1; i++)
        {
            var ai = arr[first + i];
            if (!ai.Good)
                continue;
            long iDiag = ai.Start2 - ai.Start1;
            long iEnd = ai.Start2 + ai.Len;
            for (int j = i + 1; j < count && arr[first + j].Start2 <= iEnd; j++)
            {
                var aj = arr[first + j];
                if (!aj.Good)
                    continue;
                long jDiag = aj.Start2 - aj.Start1;
                if (iDiag == jDiag)
                {
                    long jExtent = aj.Len + aj.Start2 - ai.Start2;
                    if (jExtent > ai.Len)
                    {
                        ai.Len = jExtent;
                        iEnd = ai.Start2 + jExtent;
                    }

                    aj.Good = false;
                }
                else if (ai.Start1 == aj.Start1)
                {
                    long olap = ai.Start2 + ai.Len - aj.Start2;
                    if (ResolveRepeat(ai, aj, olap))
                        break;
                }
                else if (ai.Start2 == aj.Start2)
                {
                    long olap = ai.Start1 + ai.Len - aj.Start1;
                    if (ResolveRepeat(ai, aj, olap))
                        break;
                }
            }
        }

        // std::remove_if(!Good): the kept matches keep their order
        int kept = 0;
        for (int i = 0; i < count; i++)
        {
            if (arr[first + i].Good)
                arr[first + kept++] = arr[first + i];
        }

        return kept;
    }

    // the two repeat branches of Filter_Matches (same start in one sequence); true = leave the inner loop (ai is gone)
    private static bool ResolveRepeat(MatchT ai, MatchT aj, long olap)
    {
        if (ai.Len < aj.Len)
        {
            if (olap >= ai.Len / 2)
            {
                ai.Good = false;
                return true;
            }
        }
        else if (aj.Len < ai.Len)
        {
            if (olap >= aj.Len / 2)
                aj.Good = false;
        }
        else if (olap >= ai.Len / 2)
        {
            aj.Tentative = true;
            if (ai.Tentative)
            {
                ai.Good = false;
                return true;
            }
        }

        return false;
    }

    // Process_Cluster on A[first .. first + count - 1]
    private static void ProcessCluster(MatchT[] arr, int first, int count, NucmerClusterOptions o, List<NucmerCluster> output)
    {
        var a = arr.AsSpan(first, count).ToArray();
        int n = count;
        while (n > 0)
        {
            for (int i = 0; i < n; i++)
            {
                a[i].SimpleScore = a[i].Len;
                a[i].SimpleAdj = 0;
                a[i].SimpleFrom = -1;
                for (int j = 0; j < i; j++)
                {
                    long olap1 = a[j].Start1 + a[j].Len - a[i].Start1;
                    long olap2 = a[j].Start2 + a[j].Len - a[i].Start2;
                    long olap = Math.Max(Math.Max(0, olap1), olap2);
                    long pen = olap + Math.Abs(a[i].Start2 - a[i].Start1 - (a[j].Start2 - a[j].Start1));
                    if (a[j].SimpleScore + a[i].Len - pen > a[i].SimpleScore)
                    {
                        a[i].SimpleFrom = j;
                        a[i].SimpleScore = a[j].SimpleScore + a[i].Len - pen;
                        a[i].SimpleAdj = olap;
                    }
                }
            }

            int best = 0;
            for (int i = 1; i < n; i++)
                if (a[i].SimpleScore > a[best].SimpleScore)
                    best = i;
            long total = 0, hi = long.MinValue, lo = long.MaxValue;
            for (int i = best; i >= 0; i = a[i].SimpleFrom)
            {
                a[i].Good = true;
                total += a[i].Len;
                hi = Math.Max(hi, a[i].Start1 + a[i].Len);
                lo = Math.Min(lo, a[i].Start1);
            }

            long score = o.UseExtent ? hi - lo : total;
            if (score >= o.MinClusterLength)
                output.Add(Print(a.Take(n).Where(m => m.Good).ToList()));

            // compact: drop the chain
            int kept = 0;
            for (int i = 0; i < n; i++)
                if (!a[i].Good)
                    a[kept++] = a[i];
            n = kept;
        }
    }

    // Print_Cluster: the first match as it is, every later one shortened by its overlap
    private static NucmerCluster Print(List<MatchT> cluster)
    {
        var matches = new List<NucmerClusterMatch> { new(cluster[0].Start1, cluster[0].Start2, cluster[0].Len, 0) };
        for (int i = 1; i < cluster.Count; i++)
        {
            long adj = cluster[i].SimpleAdj;
            matches.Add(new NucmerClusterMatch(cluster[i].Start1 + adj, cluster[i].Start2 + adj, cluster[i].Len - adj, adj));
        }

        return new NucmerCluster(matches);
    }
}
