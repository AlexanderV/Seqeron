#pragma warning disable S2325, S4144, S1144, S3453, CA1822 // stub of the red run (work package 1b.4): every member throws NotImplementedException; the implementation replaces this file
using SuffixTree.Indexing;

namespace SuffixTree.Esa;

internal sealed partial class EsaIndex
{
    /// <summary>
    /// The navigation pack (ISA and the LCP pyramid), built on first use: once, thread-safe, heap-resident. Core operations (search, count, locate, <c>Traverse</c>, the native analysis of
    /// 1b.3) never call it. Stub of the red run.
    /// </summary>
    internal EsaNavigationPack Navigation => throw new NotImplementedException();

    /// <summary>True once <see cref="Navigation"/> has been built (the indicator that work package 1b.5 maps to the capability flag).</summary>
    internal bool NavigationBuilt => throw new NotImplementedException();

    /// <summary>How many times the pack was built (1 after the first use, however many threads asked).</summary>
    internal int NavigationBuildCount => throw new NotImplementedException();

    /// <summary><see cref="IndexSupport.Degraded"/> until the pack exists (the first suffix-link user pays O(n)), <see cref="IndexSupport.Native"/> afterwards.</summary>
    internal IndexSupport SuffixLinkSupport => throw new NotImplementedException();

    /// <summary>A navigator over this index (builds the pack if it does not exist yet).</summary>
    internal EsaNavigator CreateNavigator() => throw new NotImplementedException();
}
