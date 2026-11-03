using SuffixTree.Esa.Tests.TestSupport;

namespace SuffixTree.Esa.Tests.Fuzz;

/// <summary>A detected mismatch: <see cref="Op"/> names the failed comparison, <see cref="Detail"/> describes the first difference.</summary>
internal readonly record struct FuzzMismatch(string Op, string Detail)
{
    /// <summary>
    /// The class of the failed comparison: <c>sa</c> for <c>sa-*</c>, <c>lcp</c> for <c>lcp-*</c>, otherwise the operation itself. The minimizer keeps
    /// a failure of the same class, because shrinking a text below <see cref="FuzzOracle.NaiveLimit"/> switches the oracle (verifier to naive) and with it the
    /// name of the failed comparison, not the defect.
    /// </summary>
    public string Class => Op.StartsWith("sa-", StringComparison.Ordinal) ? "sa" : Op.StartsWith("lcp-", StringComparison.Ordinal) ? "lcp" : Op;
}

/// <summary>
/// The oracles of the differential fuzz harness v1 (SA/LCP level). Up to <see cref="NaiveLimit"/> code units the candidate is compared with the
/// definitional naive suffix array and LCP; above it with the linear independent checks (<see cref="SuffixArrayVerifier"/> and Kasai's LCP from
/// the candidate's own SA, which is valid once the verifier accepted the SA). An optional reference engine is compared entry by entry
/// (optimized kernels against <see cref="SaisReference"/>, later work packages). Operations: <c>sa-naive</c>, <c>lcp-naive</c>, <c>sa-verify</c>,
/// <c>lcp-kasai</c>, <c>sa-reference</c>, <c>lcp-reference</c>, <c>shape</c>, <c>exception:&lt;Type&gt;</c>.
/// </summary>
internal static class FuzzOracle
{
    /// <summary>Largest text length compared with the naive O(n^2 log n) oracle.</summary>
    internal const int NaiveLimit = 2000;

    /// <summary>The reference pipeline: text preparation, <see cref="SaisReference"/>, <see cref="LcpReference"/>.</summary>
    internal static readonly EsaBuilder ReferencePipeline = EsaPipeline.Build;

    /// <summary>Returns null when the candidate agrees with the oracles, otherwise the first mismatch.</summary>
    internal static FuzzMismatch? Check(string text, EsaBuilder candidate, EsaBuilder? reference = null)
    {
        int n = text.Length;
        int[] sa, lcp;
        try
        {
            var (s, l) = candidate(text);
            sa = s;
            if (l.Length != n + 1) return new FuzzMismatch("shape", $"LCP has {l.Length} entries, expected {n + 1}");
            lcp = l.ToArray();
        }
        catch (Exception ex)
        {
            return new FuzzMismatch("exception:" + ex.GetType().Name, ex.Message);
        }
        if (sa.Length != n + 1) return new FuzzMismatch("shape", $"SA has {sa.Length} entries, expected {n + 1}");

        if (n <= NaiveLimit)
        {
            int[] naiveSa = NaiveSuffixArray.Build(text);
            int d = FirstDifference(sa, naiveSa);
            if (d >= 0) return new FuzzMismatch("sa-naive", $"SA[{d}] = {sa[d]}, naive {naiveSa[d]}");
            int[] naiveLcp = NaiveLcp.Build(text, naiveSa);
            d = FirstDifference(lcp, naiveLcp);
            if (d >= 0) return new FuzzMismatch("lcp-naive", $"LCP[{d}] = {lcp[d]}, naive {naiveLcp[d]}");
        }
        else
        {
            string? verdict = SuffixArrayVerifier.Check(text, sa);
            if (verdict is not null) return new FuzzMismatch("sa-verify", verdict);
            int[] kasai = KasaiLcp.Build(text, sa);
            int d = FirstDifference(lcp, kasai);
            if (d >= 0) return new FuzzMismatch("lcp-kasai", $"LCP[{d}] = {lcp[d]}, Kasai {kasai[d]}");
        }

        if (reference is not null)
        {
            var (rsa, rlcp) = reference(text);
            int d = FirstDifference(sa, rsa);
            if (d >= 0) return new FuzzMismatch("sa-reference", $"SA[{d}] = {sa[d]}, reference {rsa[d]}");
            d = FirstDifference(lcp, rlcp.ToArray());
            if (d >= 0) return new FuzzMismatch("lcp-reference", $"LCP[{d}] = {lcp[d]}, reference {rlcp[d]}");
        }
        return null;
    }

    private static int FirstDifference(int[] a, int[] b)
    {
        int k = Math.Min(a.Length, b.Length);
        for (int i = 0; i < k; i++)
            if (a[i] != b[i]) return i;
        return a.Length == b.Length ? -1 : k;
    }
}
