// RNA-PKPREDICT-001 — exhaustive reference for the single-knot search space (test-only, deliberately table-free).
//
// Enumerates EVERY nested pair set of the flanks and of the three knot loops and EVERY canonical knot
// (a·u·b·v·a'·w·b', Reeder & Giegerich 2004 rules 1-3) explicitly, scoring nested parts with the shared
// structure evaluator (`EvaluateStructureEnergy`, = ViennaRNA eval_structure + NNDB terms) and the knot with the
// documented formula. It shares no code with the dynamic programme (no maxPairs / helix tables, no span order).

using Seqeron.Genomics.Analysis;

namespace Seqeron.Genomics.Tests.Unit.Analysis;

internal sealed class PknotsRgBruteForce
{
    private const double PkInit = 9.0, Npp = 0.3;
    private readonly string _seq;
    private readonly int _n;
    private readonly byte[] _codes;
    private readonly Dictionary<(int, int), double> _flank = new();
    private readonly Dictionary<(int, int), double> _segment = new();

    internal PknotsRgBruteForce(string seq)
    {
        _seq = seq;
        _n = seq.Length;
        _codes = TurnerModel.Encode(seq);
    }

    /// <summary>Optimum of "at most one canonical H-type knot in the exterior loop, loops folded pseudoknot-free".</summary>
    internal double SingleKnotOptimum()
    {
        double best = FlankMin(0, _n - 1);
        var m = TurnerModel.At37;
        for (int i = 0; i < _n; i++)
            for (int j = i + 10; j < _n; j++)
                for (int k = i + 1; k < j; k++)
                {
                    int alpha = Chain(i, k);
                    if (alpha < 2) continue;
                    for (int l = i + alpha + 1; l < j; l++)
                    {
                        int beta0 = Chain(l, j);
                        if (beta0 < 2) continue;
                        int beta = Math.Min(beta0, Math.Min(k - alpha - l + 1, j - k - 2));
                        if (beta < 2) continue;
                        double e = PkInit + 3 * Npp;
                        for (int t = 0; t < alpha - 1; t++) e += m.InteriorTerm(_codes, i + t, k - t, i + t + 1, k - t - 1);
                        for (int t = 0; t < beta - 1; t++) e += m.InteriorTerm(_codes, l + t, j - t, l + t + 1, j - t - 1);
                        e += m.ExteriorStem(_codes, i, k) + m.ExteriorStem(_codes, l, j)
                           + m.InnerStemEnd(_codes, i + alpha - 1, k - alpha + 1)
                           + m.InnerStemEnd(_codes, l + beta - 1, j - beta + 1);
                        e += Segment(i + alpha + 1, l - 1) + Segment(l + beta, k - alpha) + Segment(k + 1, j - beta - 2);
                        e += FlankMin(0, i - 1) + FlankMin(j + 1, _n - 1);
                        if (e < best) best = e;
                    }
                }
        return best;
    }

    // number of consecutive pairs (p,q),(p+1,q-1),… (strands must not meet)
    private int Chain(int p, int q)
    {
        int t = 0;
        while (p + t < q - t && TurnerModel.CanPair(_codes[p + t], _codes[q - t])) t++;
        return t;
    }

    // min over all nested pair sets inside [lo,hi] of the structure energy (exterior loop, unpaired nucleotides free)
    private double FlankMin(int lo, int hi)
    {
        if (hi < lo) return 0.0;
        if (_flank.TryGetValue((lo, hi), out double v)) return v;
        double best = double.MaxValue;
        foreach (var set in AllNested(lo, hi))
            best = Math.Min(best, Evaluate(set));
        return _flank[(lo, hi)] = best;
    }

    // knot-loop segment: min over nested sets of structure energy + 0.3 per unpaired nucleotide
    private double Segment(int lo, int hi)
    {
        if (hi < lo) return 0.0;
        if (_segment.TryGetValue((lo, hi), out double v)) return v;
        double best = double.MaxValue;
        foreach (var set in AllNested(lo, hi))
            best = Math.Min(best, Evaluate(set) + Npp * (hi - lo + 1 - 2 * set.Count));
        return _segment[(lo, hi)] = best;
    }

    private double Evaluate(List<(int, int)> set)
    {
        var db = new char[_n];
        Array.Fill(db, '.');
        foreach (var (a, b) in set) { db[a] = '('; db[b] = ')'; }
        return RnaSecondaryStructure.EvaluateStructureEnergy(_seq, new string(db));
    }

    private IEnumerable<List<(int, int)>> AllNested(int lo, int hi)
    {
        if (hi < lo) { yield return new List<(int, int)>(); yield break; }
        foreach (var rest in AllNested(lo + 1, hi)) yield return rest;           // lo unpaired
        for (int q = lo + 4; q <= hi; q++)                                       // lo pairs with q (hairpin >= 3)
        {
            if (!TurnerModel.CanPair(_codes[lo], _codes[q])) continue;
            foreach (var inner in AllNested(lo + 1, q - 1))
                foreach (var after in AllNested(q + 1, hi))
                {
                    var set = new List<(int, int)>(inner.Count + after.Count + 1) { (lo, q) };
                    set.AddRange(inner);
                    set.AddRange(after);
                    yield return set;
                }
        }
    }
}
