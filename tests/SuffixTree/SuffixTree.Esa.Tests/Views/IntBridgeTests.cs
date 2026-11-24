namespace SuffixTree.Esa.Tests.Views;

/// <summary>The checked bridges to the 32-bit engine (work package S2.3): the boundary at <see cref="int.MaxValue"/> and the numbers that a narrowing cast would wrap (2^32 + 1 is 1 as a 32-bit number).</summary>
[TestFixture]
public class IntBridgeTests
{
    [Test]
    public void Length_TakesEveryIntegerFromZeroToIntMaxValue_AndNothingAbove()
    {
        Assert.That(IntBridge.Length(0, "x"), Is.EqualTo(0));
        Assert.That(IntBridge.Length(1, "x"), Is.EqualTo(1));
        Assert.That(IntBridge.Length(int.MaxValue - 1, "x"), Is.EqualTo(int.MaxValue - 1));
        Assert.That(IntBridge.Length(int.MaxValue, "x"), Is.EqualTo(int.MaxValue));
        foreach (long length in new[] { int.MaxValue + 1L, (1L << 32) - 1, 1L << 32, (1L << 32) + 1, (1L << 33) + 5, long.MaxValue, -1, -2, int.MinValue, long.MinValue })
        {
            Assert.That(() => IntBridge.Length(length, "the thing"), Throws.TypeOf<NotSupportedException>().With.Message.Contains("the thing").And.Message.Contains(length.ToString(System.Globalization.CultureInfo.InvariantCulture)), $"length {length}");
        }
    }

    [Test]
    public void Value_TakesEveryIntegerFromZeroToIntMaxValue_AndNothingAbove()
    {
        Assert.That(IntBridge.Value(0, "x"), Is.EqualTo(0));
        Assert.That(IntBridge.Value(int.MaxValue - 1, "x"), Is.EqualTo(int.MaxValue - 1));
        Assert.That(IntBridge.Value(int.MaxValue, "x"), Is.EqualTo(int.MaxValue));
        foreach (long value in new[] { int.MaxValue + 1L, (1L << 32) - 1, 1L << 32, (1L << 32) + 1, (1L << 33) + 5, long.MaxValue, -1, -2, int.MinValue, long.MinValue })
        {
            Assert.That(() => IntBridge.Value(value, "the entry"), Throws.TypeOf<InvalidDataException>().With.Message.Contains("the entry").And.Message.Contains("damaged"), $"value {value}");
        }
    }
}
