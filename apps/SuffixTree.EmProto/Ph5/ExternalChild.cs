using System.Diagnostics;

namespace SuffixTree.EmProto.Ph5;

/// <summary>The files of a PH5 run (null = off): CHILD of layout A, PSV and NSV of layout F.</summary>
internal sealed record ChildBuildOutputs
{
    internal string? Child { get; init; }

    internal string? Psv { get; init; }

    internal string? Nsv { get; init; }
}

/// <summary>What a PH5 run found and did (the forward result, plus the backward pass when NSV was asked for).</summary>
internal sealed class ChildBuildStats
{
    internal required ChildScanResult Forward { get; init; }

    internal NsvResult? Backward { get; init; }

    /// <summary>Wall seconds of the forward pass including the reading of the LCP file.</summary>
    internal double ForwardWallSeconds { get; init; }

    /// <summary>Bytes of the LCP files read (base and overflow) by the forward pass, and by the backward pass.</summary>
    internal long LcpBytes { get; init; }
}

/// <summary>PH5 of ADR-004 on a given compact LCP file (the E4 output form), for any prepared SA / LCP.</summary>
internal static class ExternalChild
{
    /// <summary>
    /// Runs the forward scan (layout A to <see cref="ChildBuildOutputs.Child"/> and / or the PSV to <see cref="ChildBuildOutputs.Psv"/>, always the statistics) over the LCP file, then, when
    /// <see cref="ChildBuildOutputs.Nsv"/> is set, the backward scan that reads the file backwards.
    /// </summary>
    internal static ChildBuildStats BuildFromFiles(long n, string lcpBase, string lcpOverflow, ChildBuildOutputs outputs, ChildBuildOptions options)
    {
        ArgumentNullException.ThrowIfNull(outputs);
        ArgumentNullException.ThrowIfNull(options);
        ChildScanResult forward;
        Task<NsvResult>? parallel = outputs.Nsv is not null && options.ParallelBackward ? Task.Run(() => Backward(n, lcpBase, lcpOverflow, outputs.Nsv, options)) : null;
        double wall;
        try
        {
            long t0 = Stopwatch.GetTimestamp();
            using (var file = new CompactLcpFile(lcpBase, lcpOverflow, n, options.ReadBufferBytes))
            using (ChildScan scan = ChildScan.Create(n, new ChildScanOutputs { ChildPath = outputs.Child, PsvPath = outputs.Psv }, options))
            {
                file.ForEachForward(scan.Consume);
                forward = scan.Complete();
            }
            wall = (Stopwatch.GetTimestamp() - t0) / (double)Stopwatch.Frequency;
        }
        catch
        {
            if (parallel is not null) ((IAsyncResult)parallel).AsyncWaitHandle.WaitOne();
            throw;
        }
        NsvResult? backward = null;
        if (parallel is not null) backward = parallel.GetAwaiter().GetResult();
        else if (outputs.Nsv is not null) backward = Backward(n, lcpBase, lcpOverflow, outputs.Nsv, options);
        return new ChildBuildStats { Forward = forward, Backward = backward, ForwardWallSeconds = wall, LcpBytes = new FileInfo(lcpBase).Length + new FileInfo(lcpOverflow).Length };
    }

    private static NsvResult Backward(long n, string lcpBase, string lcpOverflow, string nsvPath, ChildBuildOptions options)
    {
        using var file = new CompactLcpFile(lcpBase, lcpOverflow, n, options.ReadBufferBytes);
        return NsvScan.Run(n, file, nsvPath, options);
    }
}
