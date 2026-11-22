using System.Reflection;

namespace SuffixTree.Esa.Tests.Search;

/// <summary>
/// The decision of work package 1b.9g about the compilation of the mapped descent (ADR-004 section 5a, "1b.9g"): <c>MappedSearch.Run</c> carries no <c>AggressiveOptimization</c>. That attribute compiles the loop once, fully optimised and
/// without a profile, and keeps it out of tiered compilation and dynamic PGO; the descent of a short pattern, a few loads per step in the cache, ran 25 % slower than the same loop of the in-memory tree, which is profiled
/// (locate of an 8-mer on 10 Mbp of DNA, index level: 288 ns against 200 ns; without the attribute 198 ns against 200 ns). The price is the warm-up of the tiered code, about 10-40 ms in the first 10,000-20,000 queries of a
/// process. The attribute is performance only, so no behaviour test can see it; this test pins the choice so that a change of it is a change of this file and of the ADR.
/// </summary>
[TestFixture]
public class MappedSearchTiering1b9gTests
{
    [Test]
    public void TheMappedDescent_IsLeftToTieredCompilation_SoThatItGetsTheProfileOfDynamicPgo()
    {
        MethodInfo run = typeof(MappedSearch).GetMethod("Run", BindingFlags.Instance | BindingFlags.NonPublic)!;
        Assert.That(run.IsGenericMethodDefinition, Is.True, "the descent is one generic method, instantiated for the access mode, the LCP width and the child codec");
        Assert.That(run.MethodImplementationFlags.HasFlag(MethodImplAttributes.AggressiveOptimization), Is.False);
    }
}
