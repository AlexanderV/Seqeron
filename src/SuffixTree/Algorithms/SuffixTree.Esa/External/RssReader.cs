#pragma warning disable S2325, CA1822, S1144, IDE0051, IDE0060, S1172 // red stub of work package 1b.7a: every member throws NotImplementedException; the implementation replaces this file

namespace SuffixTree.Esa.External;

/// <summary>One reading of <c>/proc/self/status</c> in KiB as the kernel prints it; <see cref="Available"/> is false where the file or its <c>RssAnon</c> line does not exist.</summary>
internal readonly record struct RssSnapshot(bool Available, long RssAnonKb, long RssFileKb, long VmHwmKb, long Threads)
{
    /// <summary>The reading of a platform without the information.</summary>
    internal static RssSnapshot Unavailable => throw new NotImplementedException();

    /// <summary>Anonymous resident memory in MB (10^6 bytes).</summary>
    internal double RssAnonMb => throw new NotImplementedException();
}

/// <summary>Reads <c>RssAnon</c> (anonymous resident pages: the quantity of the working-set cap), <c>RssFile</c>, <c>VmHWM</c> and <c>Threads</c>; never throws.</summary>
internal sealed class RssReader
{
    internal RssReader(IProcSource? source = null) => throw new NotImplementedException();

    internal RssSnapshot Read() => throw new NotImplementedException();
}

/// <summary>Maxima observed by an <see cref="RssSampler"/>.</summary>
internal readonly record struct RssStats(bool Available, long Samples, long MaxRssAnonKb, long MaxRssFileKb, long MaxVmHwmKb, long MaxThreads, long LastRssAnonKb);

/// <summary>Background sampler of an <see cref="RssReader"/>; where the reader is unavailable no thread is started and <see cref="RssStats.Available"/> is false.</summary>
internal sealed class RssSampler : IDisposable
{
    internal RssSampler(RssReader reader, int intervalMs = 10) => throw new NotImplementedException();

    internal void Start() => throw new NotImplementedException();

    internal RssStats Stop() => throw new NotImplementedException();

    public void Dispose() => throw new NotImplementedException();
}
