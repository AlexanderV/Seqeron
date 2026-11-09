using SuffixTree.Indexing;

namespace SuffixTree.Esa;

internal sealed partial class EsaIndex
{
    private readonly object _navigationLock = new();
    private EsaNavigationPack? _navigation;
    private int _navigationBuilds;

    /// <summary>
    /// The navigation pack (ISA and the LCP pyramid), built on first use: once, thread-safe, heap-resident (ADR-001 R-SL). Core operations (search, count, locate, <c>Traverse</c>, the native
    /// analysis of 1b.3) never call it; the first suffix-link user pays O(n) (0.1 s at 10 Mbp).
    /// </summary>
    internal EsaNavigationPack Navigation => Volatile.Read(ref _navigation) ?? BuildNavigation();

    /// <summary>True once <see cref="Navigation"/> has been built (the indicator that work package 1b.5 maps to the capability flag).</summary>
    internal bool NavigationBuilt => Volatile.Read(ref _navigation) is not null;

    /// <summary>How many times the pack was built (1 after the first use, however many threads asked).</summary>
    internal int NavigationBuildCount => Volatile.Read(ref _navigationBuilds);

    /// <summary><see cref="IndexSupport.Degraded"/> until the pack exists (the first suffix-link user pays O(n)), <see cref="IndexSupport.Native"/> afterwards.</summary>
    internal IndexSupport SuffixLinkSupport => NavigationBuilt ? IndexSupport.Native : IndexSupport.Degraded;

    /// <summary>A navigator over this index (builds the pack if it does not exist yet).</summary>
    internal EsaNavigator CreateNavigator() => new(this);

    private EsaNavigationPack BuildNavigation()
    {
        lock (_navigationLock)
        {
            EsaNavigationPack? pack = _navigation;
            if (pack is null)
            {
                pack = EsaNavigationPack.Build(this);
                _navigationBuilds++;
                Volatile.Write(ref _navigation, pack);
            }
            return pack;
        }
    }
}
