using System.Text;
using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Analysis;

/// <summary>A named sequence record of a <see cref="MashInput"/> (FASTA: the first word of the header is the id, the rest the comment).</summary>
public readonly record struct MashRecord(string Id, string Comment, string Sequence);

/// <summary>One input file of a Mash command: its name (what the program prints, the path as given) and its records.</summary>
public sealed record MashInput(string Name, IReadOnlyList<MashRecord> Records);

/// <summary>
/// The text commands of Mash 2.3 on top of the sketching and comparison of <see cref="KmerAnalyzer"/> (audit A5-10, B09 F174): <c>mash bounds</c>, <c>mash dist</c> (with <c>-t</c>, <c>-v</c>, <c>-d</c>, <c>-i</c>) and
/// <c>mash triangle</c> (Phylip matrix, edge list <c>-E</c>, <c>-v</c>, <c>-d</c>), written as the program writes them (<c>CommandBounds.cpp</c>, <c>CommandDistance.cpp</c>, <c>CommandTriangle.cpp</c>, marbl/Mash v2.3 at raw.githubusercontent.com).
/// </summary>
public static class MashCommands
{
    // C++ `cout << double`: %g with six significant digits ("inf" / "nan" as glibc prints them)
    private static string Num(double v)
    {
        if (double.IsNaN(v))
            return "nan";
        if (double.IsInfinity(v))
            return v > 0 ? "inf" : "-inf";
        return CFormat.General(v, 6);
    }

    /// <summary>
    /// <c>mash bounds -k K -p P</c>: the table of error bounds of the Mash distance and of the screen distance for nine sketch sizes and eight distances (<c>CommandBounds.cpp</c>): for the binomial model of the number of shared
    /// hashes the smallest <c>x</c> whose cumulative probability exceeds <c>(1 - p) / 2</c> gives the lower Jaccard bound, converted back to a distance and printed as the error (bound minus distance, <c>inf</c> when
    /// no hash is certain to be shared).
    /// </summary>
    /// <param name="k">K-mer size (1 - 32, default 21).</param>
    /// <param name="probability">Confidence of the bounds (0 - 1, default 0.99).</param>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="k"/> or <paramref name="probability"/> outside the program's range.</exception>
    public static string Bounds(int k = 21, double probability = 0.99)
    {
        if (k is < 1 or > 32)
            throw new ArgumentOutOfRangeException(nameof(k), k, "k must be in 1..32 (mash bounds -k).");
        if (!(probability > 0 && probability < 1))
            throw new ArgumentOutOfRangeException(nameof(probability), probability, "p must be in (0, 1) (mash bounds -p).");
        int[] sketchSizes = [100, 500, 1000, 5000, 10000, 50000, 100000, 500000, 1000000];
        double[] dists = [0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4];
        double q2 = (1.0 - probability) / 2.0;
        var sb = new StringBuilder();
        sb.Append('\n').Append("Parameters (run with -h for details):\n").Append("   k:   ").Append(k).Append('\n').Append("   p:   ").Append(Num(probability)).Append("\n\n");
        for (int screen = 0; screen < 2; screen++)
        {
            sb.Append(screen == 1 ? "\tScreen distance\n" : "\tMash distance\n").Append("Sketch");
            foreach (double d in dists)
                sb.Append('\t').Append(Num(d));
            sb.Append('\n');
            foreach (int s in sketchSizes)
            {
                sb.Append(s);
                foreach (double dist in dists)
                {
                    double m2j = screen == 1 ? Math.Pow(1.0 - dist, k) : 1.0 / (2.0 * Math.Exp(k * dist) - 1.0);
                    // the smallest x < s whose cumulative probability exceeds q2 (s when none does): the CDF is monotone, so bisect instead of the program's linear scan
                    int lo = 0, hi = s;
                    while (lo < hi)
                    {
                        int mid = lo + (hi - lo) / 2;
                        if (StatisticsHelper.BinomialCdf(mid, s, m2j) > q2)
                            hi = mid;
                        else
                            lo = mid + 1;
                    }

                    int x = lo;
                    double je = (double)x / s;
                    double j2m = screen == 1 ? 1.0 - Math.Pow(je, 1.0 / k) : -1.0 / k * Math.Log(2.0 * je / (1.0 + je));
                    sb.Append('\t').Append(Num(j2m - dist));
                }

                sb.Append('\n');
            }

            sb.Append('\n');
        }

        return sb.ToString();
    }

    private sealed record Named(string Name, MinHashSketch Sketch);

    private static List<Named> Sketches(IReadOnlyList<MashInput> inputs, int k, int sketchSize, uint seed, bool canonical, bool individual)
    {
        var list = new List<Named>();
        foreach (var input in inputs)
        {
            if (individual)
            {
                foreach (var r in input.Records.Where(r => r.Sequence.Length >= k))              // Sketch.cpp: a sequence shorter than k is not sketched
                    list.Add(new Named(r.Id, KmerAnalyzer.CreateMinHashSketch(r.Sequence, k, sketchSize, canonical, seed)));
            }
            else
            {
                list.Add(new Named(input.Name, KmerAnalyzer.CreateMinHashSketch(input.Records.Select(r => r.Sequence), k, sketchSize, canonical, seed)));
            }
        }

        return list;
    }

    private static bool Passes(MashComparison c, double maxPValue, double maxDistance) => c.Distance <= maxDistance && c.PValue <= maxPValue;

    /// <summary>
    /// <c>mash dist [-t] [-v P] [-d D] [-i] reference... query...</c>: for every query (in input order) and every reference, one line <c>reference TAB query TAB distance TAB p-value TAB shared/total</c> of the pairs whose p-value is at most
    /// <paramref name="maxPValue"/> and distance at most <paramref name="maxDistance"/>; with <paramref name="table"/> a header <c>#query TAB reference...</c> and a row per query with the distances (a failing pair leaves its field blank).
    /// </summary>
    /// <param name="references">The reference files.</param>
    /// <param name="queries">The query files.</param>
    /// <param name="k">K-mer size (<c>-k</c>, 1 - 32).</param>
    /// <param name="sketchSize">Sketch size (<c>-s</c>).</param>
    /// <param name="seed">Hash seed (<c>-S</c>).</param>
    /// <param name="canonical">False for <c>-n</c> (strand preserved).</param>
    /// <param name="individual"><c>-i</c>: one sketch per sequence (named by its id) instead of per file.</param>
    /// <param name="table"><c>-t</c>: table output.</param>
    /// <param name="maxPValue"><c>-v</c> (0 - 1).</param>
    /// <param name="maxDistance"><c>-d</c> (0 - 1).</param>
    public static string Dist(IReadOnlyList<MashInput> references, IReadOnlyList<MashInput> queries, int k = 21, int sketchSize = KmerAnalyzer.DefaultMashSketchSize, uint seed = KmerAnalyzer.DefaultMashSeed,
        bool canonical = true, bool individual = false, bool table = false, double maxPValue = 1.0, double maxDistance = 1.0)
    {
        ArgumentNullException.ThrowIfNull(references);
        ArgumentNullException.ThrowIfNull(queries);
        RequireFraction(maxPValue, nameof(maxPValue), "-v");
        RequireFraction(maxDistance, nameof(maxDistance), "-d");
        var refs = Sketches(references, k, sketchSize, seed, canonical, individual);
        var qs = Sketches(queries, k, sketchSize, seed, canonical, individual);
        var sb = new StringBuilder();
        if (table)
        {
            sb.Append("#query");
            foreach (var r in refs)
                sb.Append('\t').Append(r.Name);
            sb.Append('\n');
        }

        foreach (var q in qs)
        {
            if (table)
                sb.Append(q.Name);
            foreach (var r in refs)
            {
                var c = KmerAnalyzer.CompareMinHashSketches(r.Sketch, q.Sketch);
                bool pass = Passes(c, maxPValue, maxDistance);
                if (table)
                {
                    sb.Append('\t');
                    if (pass)
                        sb.Append(Num(c.Distance));
                }
                else if (pass)
                {
                    sb.Append(r.Name).Append('\t').Append(q.Name).Append('\t').Append(Num(c.Distance)).Append('\t').Append(Num(c.PValue)).Append('\t').Append(c.SharedHashes).Append('/').Append(c.Denominator).Append('\n');
                }
            }

            if (table)
                sb.Append('\n');
        }

        return sb.ToString();
    }

    /// <summary>
    /// <c>mash triangle [-E] [-v P] [-d D] [-i] sequence...</c>: the lower-triangular distance matrix in relaxed Phylip format (the number of sequences, then one row per sequence: its name and the distances to the sequences before it), or with
    /// <paramref name="edgeList"/> (implied by a p-value or distance filter in the program) the edge list <c>name TAB name TAB distance TAB p-value TAB shared/total</c> of the pairs that pass the filters.
    /// </summary>
    public static string Triangle(IReadOnlyList<MashInput> inputs, int k = 21, int sketchSize = KmerAnalyzer.DefaultMashSketchSize, uint seed = KmerAnalyzer.DefaultMashSeed,
        bool canonical = true, bool individual = false, bool edgeList = false, double maxPValue = 1.0, double maxDistance = 1.0)
    {
        ArgumentNullException.ThrowIfNull(inputs);
        RequireFraction(maxPValue, nameof(maxPValue), "-v");
        RequireFraction(maxDistance, nameof(maxDistance), "-d");
        edgeList |= maxPValue < 1.0 || maxDistance < 1.0;
        var sketches = Sketches(inputs, k, sketchSize, seed, canonical, individual);
        var sb = new StringBuilder();
        if (!edgeList)
            sb.Append('\t').Append(sketches.Count).Append('\n');
        for (int index = 0; index < sketches.Count; index++)
        {
            if (!edgeList)
                sb.Append(sketches[index].Name);
            for (int i = 0; i < index; i++)
            {
                var c = KmerAnalyzer.CompareMinHashSketches(sketches[i].Sketch, sketches[index].Sketch);
                if (edgeList)
                {
                    if (Passes(c, maxPValue, maxDistance))
                        sb.Append(sketches[index].Name).Append('\t').Append(sketches[i].Name).Append('\t').Append(Num(c.Distance)).Append('\t').Append(Num(c.PValue)).Append('\t').Append(c.SharedHashes).Append('/').Append(c.Denominator).Append('\n');
                }
                else
                {
                    sb.Append('\t').Append(Num(c.Distance));
                }
            }

            if (!edgeList)
                sb.Append('\n');
        }

        return sb.ToString();
    }

    private static void RequireFraction(double value, string name, string option)
    {
        if (!(value >= 0 && value <= 1))
            throw new ArgumentOutOfRangeException(name, value, $"{option} must be in 0..1 (mash).");
    }
}
