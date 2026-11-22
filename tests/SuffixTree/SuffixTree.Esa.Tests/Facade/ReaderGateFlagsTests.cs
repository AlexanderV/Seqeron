using System.Reflection;
using System.Runtime.InteropServices;

namespace SuffixTree.Esa.Tests.Facade;

/// <summary>
/// Where the flags of a gate live (work package 1b.9g.1). Every call of a tree over a container reads the flag of its gate, twice; the gate object is small (about 56 bytes) and sits among whatever the allocator placed next to it, and in some processes
/// that neighbour is written often: the line of the flag then bounces between cores at four threads (<c>ReadBench mt:</c> on DNA, 0.95 of the index, 1b.9g; the same build with the flags in an object of their own: 1.00). The flags are therefore in a
/// separate object of explicit layout whose fields have at least a cache line of room on each side, one object per gate. The layout is pinned here by reflection (a machine cannot show the neighbour that is written, so the room is the property to test);
/// the behaviour that must not change (a gate disposes, releases and refuses on its own flags) is tested by <see cref="ReaderGateTests"/>.
/// </summary>
[TestFixture]
public class ReaderGateFlagsTests
{
    private const BindingFlags Instance = BindingFlags.Instance | BindingFlags.NonPublic | BindingFlags.Public;

    private static ReaderGate<ReaderGateSeams> NewGate(Action? release = null) => new(new object(), release ?? (() => { }));

    private static FieldInfo FlagsField() => typeof(ReaderGate<ReaderGateSeams>).GetField("_flags", Instance)
        ?? throw new AssertionException("the gate has no field `_flags`: its flags are not in an object of their own");

    private static object FlagsOf(ReaderGate<ReaderGateSeams> gate) => FlagsField().GetValue(gate) ?? throw new AssertionException("the flags object of the gate is null");

    [Test]
    public void TheFlagsEveryCallReads_LiveInAnObjectOfTheirOwn_ThatIsNotTheGate()
    {
        FieldInfo field = FlagsField();
        Assert.Multiple(() =>
        {
            Assert.That(field.FieldType.IsClass && !field.FieldType.IsValueType, Is.True, "a reference to a separate object, not a struct inside the gate");
            Assert.That(field.IsInitOnly, Is.True, "set once, in the constructor: a call reads the reference and then the flag");
            FieldInfo[] own = typeof(ReaderGate<ReaderGateSeams>).GetFields(Instance);
            Assert.That(own.Count(f => f.FieldType == typeof(int)), Is.EqualTo(0), "no flag is left among the fields of the gate");
        });
    }

    [Test]
    public void TheFlagsObject_HasAtLeastACacheLineOfRoomOnBothSidesOfTheFlags()
    {
        Type flags = FlagsField().FieldType;
        StructLayoutAttribute layout = flags.StructLayoutAttribute!;
        FieldInfo[] fields = flags.GetFields(Instance);
        Assert.That(layout.Value, Is.EqualTo(LayoutKind.Explicit));
        Assert.That(fields.Select(f => f.Name), Is.EquivalentTo(new[] { "Disposed", "Released" }));
        int[] offsets = [.. fields.Select(f => f.GetCustomAttribute<FieldOffsetAttribute>()!.Value)];
        Assert.Multiple(() =>
        {
            Assert.That(offsets.Min(), Is.GreaterThanOrEqualTo(64), "a cache line before the first flag (inside the object: the header and the table pointer come on top)");
            Assert.That(layout.Size - (offsets.Max() + sizeof(int)), Is.GreaterThanOrEqualTo(64), "and a cache line after the last one");
            Assert.That(offsets.Distinct().Count(), Is.EqualTo(2), "the two flags are two words, not one");
            Assert.That(offsets.Max() - offsets.Min(), Is.LessThan(64), "the flag a call reads and the one the release reads share one line: nothing else is on it");
        });
    }

    [Test]
    public void EveryGateHasItsOwnFlagsObject_AndTheFlagsOfTwoGatesAreIndependent()
    {
        int releasedA = 0, releasedB = 0;
        var a = NewGate(() => releasedA++);
        var b = NewGate(() => releasedB++);
        Assert.That(FlagsOf(a), Is.Not.SameAs(FlagsOf(b)), "one flags object per gate");

        a.Dispose();

        Assert.Multiple(() =>
        {
            Assert.That(a.IsDisposed, Is.True);
            Assert.That(a.IsReleased, Is.True);
            Assert.That(releasedA, Is.EqualTo(1));
            Assert.That(b.IsDisposed, Is.False, "the flag of another gate is not raised");
            Assert.That(b.IsReleased, Is.False);
            Assert.That(releasedB, Is.EqualTo(0));
        });
        Assert.DoesNotThrow(() => b.Leave(b.Enter()));
        Assert.Throws<ObjectDisposedException>(() => a.Enter());
        b.Dispose();
        Assert.That(releasedB, Is.EqualTo(1));
    }

    [Test]
    public void TheFlagsObject_IsCreatedWithTheGate_AndStartsLowered()
    {
        var gate = NewGate();
        object flags = FlagsOf(gate);
        Type type = flags.GetType();
        Assert.Multiple(() =>
        {
            Assert.That(type.GetField("Disposed", Instance)!.GetValue(flags), Is.EqualTo(0));
            Assert.That(type.GetField("Released", Instance)!.GetValue(flags), Is.EqualTo(0));
            Assert.That(gate.IsDisposed, Is.False);
            Assert.That(gate.IsReleased, Is.False);
        });
        gate.Dispose();
        Assert.Multiple(() =>
        {
            Assert.That(type.GetField("Disposed", Instance)!.GetValue(flags), Is.EqualTo(1), "Dispose raises the flag in the flags object");
            Assert.That(type.GetField("Released", Instance)!.GetValue(flags), Is.EqualTo(1), "and the release claim is in it too");
        });
    }
}
