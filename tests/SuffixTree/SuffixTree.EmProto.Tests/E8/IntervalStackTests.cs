using SuffixTree.EmProto.Verify;

namespace SuffixTree.EmProto.Tests.E8;

/// <summary>The bounded stack of the child scan: blocks allocated on demand, a hard capacity, the peak.</summary>
[TestFixture]
public class IntervalStackTests
{
    [Test]
    public void PushPopTop_AcrossBlocks_KeepsEveryEntry()
    {
        foreach (int blockBits in new[] { 0, 1, 2, 3, 20 })
        {
            var stack = new IntervalStack(100, blockBits);
            for (uint i = 0; i < 70; i++)
            {
                Assert.That(stack.TryPush(new OpenInterval { Ell = i, Lb = i + 1, First = i + 2, Last = i + 3 }), Is.True);
                Assert.That(stack.Top.Ell, Is.EqualTo(i));
                Assert.That(stack.Count, Is.EqualTo(i + 1));
            }
            for (int i = 69; i >= 0; i--)
            {
                OpenInterval e = stack.Top;
                Assert.That((e.Ell, e.Lb, e.First, e.Last), Is.EqualTo(((uint)i, (uint)i + 1, (uint)i + 2, (uint)i + 3)), $"block bits {blockBits}, entry {i}");
                stack.Pop();
            }
            Assert.That(stack.Count, Is.EqualTo(0));
            Assert.That(stack.Peak, Is.EqualTo(70));
        }
    }

    [Test]
    public void TheTopCanBeUpdatedInPlace()
    {
        var stack = new IntervalStack(10, 1);
        for (uint i = 0; i < 5; i++) stack.TryPush(new OpenInterval { Ell = i, First = IntervalStack.None, Last = IntervalStack.None });
        stack.Top.First = 42;
        stack.Top.Last = 43;
        Assert.That((stack.Top.First, stack.Top.Last), Is.EqualTo((42u, 43u)));
        stack.Pop();
        Assert.That(stack.Top.First, Is.EqualTo(IntervalStack.None));
    }

    [Test]
    public void TheCapacityIsHard_AndATooFullPushChangesNothing()
    {
        var stack = new IntervalStack(5, 2);
        for (uint i = 0; i < 5; i++) Assert.That(stack.TryPush(new OpenInterval { Ell = i }), Is.True);
        Assert.That(stack.TryPush(new OpenInterval { Ell = 99 }), Is.False);
        Assert.That(stack.Count, Is.EqualTo(5));
        Assert.That(stack.Top.Ell, Is.EqualTo(4u));
        stack.Pop();
        Assert.That(stack.TryPush(new OpenInterval { Ell = 7 }), Is.True);
        Assert.That(stack.Peak, Is.EqualTo(5));
        Assert.That(stack.Capacity, Is.EqualTo(5));
    }

    [Test]
    public void AnEmptyStack_HasNoBlocksAndAZeroPeak()
    {
        var stack = new IntervalStack(1L << 40);
        Assert.That(stack.Count, Is.EqualTo(0));
        Assert.That(stack.Peak, Is.EqualTo(0));
    }
}
