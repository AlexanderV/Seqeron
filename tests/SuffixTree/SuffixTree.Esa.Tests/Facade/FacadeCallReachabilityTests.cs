using System.Reflection;
using SuffixTree.Esa.Tests.External;
using SuffixTree.Esa.Tests.TestSupport;

namespace SuffixTree.Esa.Tests.Facade;

/// <summary>
/// The scope of a call keeps its tree reachable for the length of the call (work package 1b.9g restores what 1b.9b stated and no test held: "the scope still holds the tree"). The navigation pack of a tree pins the suffix array
/// with GC handles that the finalizer of its <c>PinnedArray</c> frees, a reader walks suffix links through raw pointers, and an optimizing JIT regards <c>this</c> as dead after its last use: when a member of the facade has handed the
/// index to the callee, the tree may be unreachable while the callee still reads, unless the scope of the call holds it. That is a fact of optimized code that a test cannot force (tier-0 code keeps every local alive, tier-1 code
/// appears after the call counts of the tiered JIT; a collection inside a call showed no difference in either), so the test reads the scope itself: the value that <c>Enter</c> returns for an in-memory tree holds the tree, and for a
/// tree over a container it holds the gate, which holds the tree.
/// </summary>
[TestFixture]
public class FacadeCallReachabilityTests
{
    private static object Enter(EsaSuffixTree tree) => typeof(EsaSuffixTree).GetMethod("Enter", BindingFlags.Instance | BindingFlags.NonPublic)!.Invoke(tree, null)!;

    private static List<object?> FieldValues(object o) => [.. o.GetType().GetFields(BindingFlags.Instance | BindingFlags.NonPublic | BindingFlags.Public).Select(f => f.GetValue(o))];

    [Test]
    public void TheScopeOfACallOfAnInMemoryTree_HoldsTheTree()
    {
        using EsaSuffixTree tree = EsaSuffixTree.Build(StringFamilies.RandomDna(500, 31));
        object scope = Enter(tree);
        try
        {
            Assert.That(FieldValues(scope).Any(v => ReferenceEquals(v, tree)), Is.True, "the scope holds the tree: the collector may not finalize the pins of its pack under a reader");
        }
        finally
        {
            ((IDisposable)scope).Dispose();
        }
    }

    [Test]
    public void TheScopeOfACallOfATreeOverAContainer_HoldsTheGate_WhichHoldsTheTree_AndTheCellOfTheThread()
    {
        using var dir = new TempDir();
        using EsaSuffixTree tree = SesaTestFile.Open(SesaTestFile.Write(dir, StringFamilies.RandomDna(500, 32)));
        object scope = Enter(tree);
        try
        {
            List<object?> values = FieldValues(scope);
            object? gate = values.FirstOrDefault(v => v is not null && v.GetType().Name.StartsWith("ReaderGate", StringComparison.Ordinal));
            Assert.That(gate, Is.Not.Null, "the scope holds the gate");
            object? owner = gate!.GetType().GetField("_owner", BindingFlags.Instance | BindingFlags.NonPublic)!.GetValue(gate);
            Assert.That(owner, Is.SameAs(tree), "and the gate holds the tree");
            Assert.That(values.Contains(ReaderCell.Current), Is.True, "and the scope holds the cell of the calling thread, which it clears when the call ends");
            Assert.That(ReaderCell.AnyInside(tree.ReaderId), Is.True);
        }
        finally
        {
            ((IDisposable)scope).Dispose();
        }

        Assert.That(ReaderCell.AnyInside(tree.ReaderId), Is.False);
    }
}
