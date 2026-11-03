using System.Numerics;

namespace SuffixTree.Esa;

/// <summary>
/// Reference LCP construction by Phi / PLCP (Karkkainen, Manzini, Puglisi: "Permuted longest-common-prefix array", CPM 2009;
/// ADR-002 section 4), in safe code, generic over the element type of the suffix array like <see cref="SaisReference{TIdx}"/>
/// (the <c>long</c> instantiation is an oracle only). Three passes, linear time:
/// <list type="number">
/// <item><c>Phi[SA[i]] = SA[i-1]</c> for <c>i = 1..n</c> (for <c>i = 1</c> the predecessor is the empty suffix <c>n</c>);</item>
/// <item>in text order <c>PLCP[j] = lcp(j, Phi[j])</c> with the Kasai invariant <c>PLCP[j+1] &gt;= PLCP[j] - 1</c> (at most 2n symbol
///   comparisons), written over <c>Phi</c> in place;</item>
/// <item><c>LCP[r] = PLCP[SA[r]]</c> for <c>r = 1..n</c>, stored directly in the smallest <see cref="CompactLcp"/> representation
///   (the number of escapes is counted on PLCP first, because the multiset of LCP values equals that of PLCP).</item>
/// </list>
/// </summary>
/// <typeparam name="TIdx">Signed integer type of the suffix array (<c>int</c> or <c>long</c>).</typeparam>
internal static class LcpReference<TIdx>
    where TIdx : unmanaged, IBinaryInteger<TIdx>, ISignedNumber<TIdx>
{
    private static int ToInt(TIdx value) => int.CreateChecked(value);

    private static TIdx FromInt(int value) => TIdx.CreateChecked(value);

    /// <summary>Builds the compact LCP array of <paramref name="text"/> from its suffix array (n + 1 layout, <c>sa[0] = n</c>).</summary>
    internal static CompactLcp Build<TReader>(TReader text, ReadOnlySpan<TIdx> sa, LcpTrace? trace = null)
        where TReader : ISymbolReader, allows ref struct
    {
        int n = text.Length;
        TIdx[] plcp = ComputePlcp(text, sa, trace);

        long tick = trace is null ? 0 : System.Diagnostics.Stopwatch.GetTimestamp();
        long escapes8 = 0, escapes16 = 0;
        foreach (TIdx value in plcp)
        {
            int v = ToInt(value);
            if (v >= 255) escapes8++;
            if (v >= 65535) escapes16++;
        }
        LcpKind kind = CompactLcp.SelectKind(n, escapes8, escapes16);
        int escapes = kind switch { LcpKind.U8 => (int)escapes8, LcpKind.U16 => (int)escapes16, _ => 0 };

        if (trace is not null)
        {
            long now = System.Diagnostics.Stopwatch.GetTimestamp();
            trace.CountTicks += now - tick;
            tick = now;
        }
        var builder = new CompactLcp.Builder(kind, n + 1, escapes);
        for (int r = 1; r <= n; r++) builder.Set(r, ToInt(plcp[ToInt(sa[r])]));
        CompactLcp result = builder.Complete();
        if (trace is not null) trace.OutputTicks += System.Diagnostics.Stopwatch.GetTimestamp() - tick;
        return result;
    }

    /// <summary>The permuted LCP array: <c>PLCP[j]</c> is the LCP of suffix j and its predecessor in suffix-array order (0 when that is the empty suffix). Length n.</summary>
    internal static TIdx[] ComputePlcp<TReader>(TReader text, ReadOnlySpan<TIdx> sa, LcpTrace? trace = null)
        where TReader : ISymbolReader, allows ref struct
    {
        int n = text.Length;
        if (sa.Length != n + 1)
            throw new ArgumentException($"the suffix array needs n + 1 = {n + 1} entries, got {sa.Length}", nameof(sa));
        if (ToInt(sa[0]) != n)
            throw new ArgumentException("sa[0] must be n, the empty suffix", nameof(sa));

        long tick = trace is null ? 0 : System.Diagnostics.Stopwatch.GetTimestamp();
        var phi = new TIdx[n];
        for (int i = 1; i <= n; i++) phi[ToInt(sa[i])] = sa[i - 1];
        if (trace is not null)
        {
            long now = System.Diagnostics.Stopwatch.GetTimestamp();
            trace.PhiTicks += now - tick;
            tick = now;
        }

        int l = 0;
        for (int j = 0; j < n; j++)
        {
            int p = ToInt(phi[j]);
            if (p == n)
            {
                l = 0;
                phi[j] = TIdx.Zero;
                continue;
            }
            while (j + l < n && p + l < n && text[j + l] == text[p + l]) l++;
            phi[j] = FromInt(l);
            if (l > 0) l--;
        }
        if (trace is not null) trace.PlcpTicks += System.Diagnostics.Stopwatch.GetTimestamp() - tick;
        return phi;
    }
}
