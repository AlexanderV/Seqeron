using System.Reflection;

namespace SuffixTree.Esa.Tests.Navigation;

/// <summary>
/// The navigator's walks over the child table end in an <see cref="InvalidOperationException"/> on a corrupt table instead of looping or reading outside the interval (work package 1b.4 keeps every
/// loop bounded: the l-indices of an interval must strictly increase and stay inside it; the <c>Checked</c> mode of 1b.7 builds on this). The tables are corrupted through reflection on the
/// private arrays of the layouts; the index is rebuilt for every case.
/// </summary>
[TestFixture]
public class CorruptTableTests
{
    private const string Text = "abracadabra";

    private static T Field<T>(object owner, string name)
        => (T)owner.GetType().GetField(name, BindingFlags.NonPublic | BindingFlags.Instance)!.GetValue(owner)!;

    private static void AssertWalksRefuse(EsaNavigator nav, string what)
    {
        var list = new List<EsaNode>();
        Assert.Throws<InvalidOperationException>(() => nav.GetChildren(nav.Root, list), $"{what}: GetChildren");
        Assert.Throws<InvalidOperationException>(() => nav.TryGetChild(nav.Root, 'z', out _), $"{what}: TryGetChild of a key behind every child");
        Assert.Throws<InvalidOperationException>(() => nav.TryGetChild(nav.Root, 'b', out _), $"{what}: TryGetChild of a key in the middle");
    }

    [Test]
    public void Reference3_SelfLoopAndBoundaryBeyondTheInterval_AreRefused()
    {
        EsaIndex index = EsaIndex.Build(Text, ChildLayoutKind.Reference3);
        EsaNavigator nav = index.CreateNavigator();
        object table = index.Child.Reference3;
        int[] next = Field<int[]>(table, "_next");
        int[] down = Field<int[]>(table, "_down");
        int first = index.Child.FirstLIndex(0, index.Length, true);
        int second = next[first];

        next[first] = first;                       // the l-index after the first points at itself
        AssertWalksRefuse(nav, "self loop");
        next[first] = index.Length + 3;            // ... or beyond the interval
        AssertWalksRefuse(nav, "boundary beyond the interval");
        next[first] = second;

        int saved = down[0];
        down[0] = 0;                               // the first l-index of the root is its own left boundary
        AssertWalksRefuse(nav, "first boundary at the left border");
        down[0] = saved;

        var list = new List<EsaNode>();
        nav.GetChildren(nav.Root, list);
        Assert.That(list, Has.Count.EqualTo(6), "the repaired table works again");
    }

    [Test]
    public void Tagged_SelfLoopAndFirstBoundaryAtTheLeftBorder_AreRefused()
    {
        EsaIndex index = EsaIndex.Build(Text, ChildLayoutKind.Tagged);
        EsaNavigator nav = index.CreateNavigator();
        object table = index.Child.Tagged;
        uint[] cld = Field<uint[]>(table, "_cld");
        int first = index.Child.FirstLIndex(0, index.Length, true);
        uint savedNext = cld[first];
        uint savedRoot = cld[0];

        cld[first] = TaggedChild.PackNext(first);
        AssertWalksRefuse(nav, "self loop");
        cld[first] = savedNext;

        cld[0] = TaggedChild.PackDown(0);
        AssertWalksRefuse(nav, "first boundary at the left border");
        cld[0] = savedRoot;

        var list = new List<EsaNode>();
        nav.GetChildren(nav.Root, list);
        Assert.That(list, Has.Count.EqualTo(6), "the repaired table works again");
    }
}
