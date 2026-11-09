namespace SuffixTree.Indexing;

/// <summary>How well an engine answers an <see cref="IndexOperation"/>.</summary>
public enum IndexSupport
{
    /// <summary>The engine cannot do it; calling the member throws <see cref="System.NotSupportedException"/> and never fakes an answer.</summary>
    Unsupported = 0,

    /// <summary>Correct, but slower than the engine class norm or needing an auxiliary structure that is built lazily; <see cref="ISuffixIndex.SupportNote"/> states the cost.</summary>
    Degraded = 1,

    /// <summary>Answered at the speed that is normal for the engine class.</summary>
    Native = 2,
}
