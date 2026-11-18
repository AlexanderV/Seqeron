namespace SuffixTree.Esa.External;

/// <summary>What the PLCP loop of one chunk found.</summary>
internal readonly record struct PlcpChunkResult(long Escapes8, long Escapes16, long CompareSteps, long MaxPlcp, long FirstValue);

/// <summary>
/// The PLCP loop of <c>LcpKernel.ComputePlcp</c> on one chunk of Phi in text order: <c>PLCP[j] = lcp(j, Phi[j])</c> with the Kasai invariant <c>PLCP[j+1] &gt;= PLCP[j] - 1</c>, written over Phi in place, the carried
/// length <c>l</c> compared with the vector common prefix of <see cref="MappedTextReader"/>, the values at or above 255 and 65,535 counted in the loop. A chunk starts with <c>l = 0</c> (chunks are processed in
/// any order by any task); <see cref="PlcpChunkResult.FirstValue"/> is the whole cost of that restart.
/// </summary>
internal static class PlcpScan
{
    /// <summary>Iterations ahead at which the text line of a later suffix is requested.</summary>
    internal const int PrefetchDistance = 16;

    /// <summary>Phi to PLCP for the positions <c>firstIndex .. firstIndex + count - 1</c>; <paramref name="phi"/> holds <c>Phi[j]</c> (the predecessor in rank order, <c>n</c> for the empty suffix) and is overwritten by <c>PLCP[j]</c>.</summary>
    internal static PlcpChunkResult Run(MappedTextReader text, long firstIndex, Buf<uint> phi, int count) => throw new NotImplementedException();
}
