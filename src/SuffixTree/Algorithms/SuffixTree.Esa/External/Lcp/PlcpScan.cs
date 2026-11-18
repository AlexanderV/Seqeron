namespace SuffixTree.Esa.External;

/// <summary>What the PLCP loop of one chunk found.</summary>
internal readonly record struct PlcpChunkResult(long Escapes8, long Escapes16, long CompareSteps, long MaxPlcp, long FirstValue);

/// <summary>
/// The PLCP loop of <c>LcpKernel.ComputePlcp</c> on one chunk of Phi in text order: <c>PLCP[j] = lcp(j, Phi[j])</c> with the Kasai invariant <c>PLCP[j+1] &gt;= PLCP[j] - 1</c>, written over Phi in place, the carried
/// length <c>l</c> compared with the vector common prefix of <see cref="MappedTextReader"/>, the values at or above 255 and 65,535 counted in the loop. The carry across chunks does not exist (the chunks are processed in
/// any order by any task), so a chunk starts with <c>l = 0</c>: its first position compares <see cref="PlcpChunkResult.FirstValue"/> symbols from scratch, the whole cost of the restart (at most the maximum PLCP of the
/// chunk; n / 2 symbols per chunk on average for a^n). Every read of the text is clamped by the reader, so a Phi that is not the one of the suffix array cannot make the loop read outside the text.
/// </summary>
internal static class PlcpScan
{
    /// <summary>Iterations ahead at which the text line of a later suffix is requested (the plateau of 8-64 in the prototype's prefetch experiment).</summary>
    internal const int PrefetchDistance = 16;

    /// <summary>Phi to PLCP for the positions <c>firstIndex .. firstIndex + count - 1</c>; <paramref name="phi"/> holds <c>Phi[j]</c> (the predecessor in rank order, <c>n</c> for the empty suffix) and is overwritten by <c>PLCP[j]</c>.</summary>
    [System.Runtime.CompilerServices.MethodImpl(System.Runtime.CompilerServices.MethodImplOptions.AggressiveOptimization)]
    internal static PlcpChunkResult Run(MappedTextReader text, long firstIndex, Buf<uint> phi, int count)
    {
        long n = text.Length;
        long steps = 0, e8 = 0, e16 = 0, max = 0, first = 0;
        long l = 0;
        for (int k = 0; k < count; k++)
        {
            if (k + PrefetchDistance < count) text.Prefetch(phi[k + PrefetchDistance]);
            long j = firstIndex + k;
            long p = phi[k];
            if (p == n)
            {
                l = 0; // j is the smallest nonempty suffix: its predecessor is the empty suffix
                phi[k] = 0;
                continue;
            }
            long limit = n - Math.Max(j, p); // symbols that both suffixes still have
            if (l < limit) l += text.Lce(j + l, p + l, limit - l, ref steps);
            phi[k] = (uint)l;
            if (k == 0) first = l;
            if (l >= 255)
            {
                e8++;
                if (l >= 65535) e16++;
            }
            if (l > max) max = l;
            if (l > 0) l--;
        }
        return new PlcpChunkResult(e8, e16, steps, max, first);
    }
}
