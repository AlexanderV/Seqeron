using System.Runtime;
using System.Runtime.CompilerServices;

namespace SuffixTree.Esa.Tests.Navigation;

/// <summary>The array whose address stays fixed (work package 1b.4): the navigation pack reads the suffix array, the LCP base array, ISA and the pyramid through pointers between calls.</summary>
[TestFixture]
public unsafe class PinnedArrayTests
{
    [Test]
    public void Allocate_ZeroedOrNot_LengthAndName()
    {
        using var zeroed = PinnedArray<int>.Allocate(1000, zero: true, "zeroed");
        Assert.That(zeroed.Length, Is.EqualTo(1000));
        Assert.That(zeroed.Name, Is.EqualTo("zeroed"));
        Assert.That(zeroed.Array.All(v => v == 0));
        Assert.That(zeroed.HoldsHandle, Is.False, "an array on the pinned object heap needs no handle");
        using var raw = PinnedArray<byte>.Allocate(0, zero: false, "empty");
        Assert.That(raw.Length, Is.EqualTo(0));
        Assert.Throws<ArgumentOutOfRangeException>(() => PinnedArray<int>.Allocate(-1, zero: true, "bad"));
    }

    [Test]
    public void Pin_ExistingArray_BufIsTheArray()
    {
        var array = new int[500];
        using var pinned = PinnedArray<int>.Pin(array, "pinned");
        Assert.That(pinned.HoldsHandle);
        Assert.That(pinned.Array, Is.SameAs(array));
        for (int i = 0; i < array.Length; i++) pinned.Buf[i] = i * 3;
        for (int i = 0; i < array.Length; i++) Assert.That(array[i], Is.EqualTo(i * 3), "a write through the Buf is a write to the array");
        array[17] = -5;
        Assert.That(pinned.Buf[17], Is.EqualTo(-5), "and the other way round");
        Assert.Throws<ArgumentNullException>(() => PinnedArray<int>.Pin(null!, "null"));
    }

    [Test]
    public void Pin_SurvivesCompactingCollections()
    {
        var array = new int[200]; // small: lives in generation 0 and would move without the pin
        for (int i = 0; i < array.Length; i++) array[i] = i;
        using var pinned = PinnedArray<int>.Pin(array, "small");
        int* before = (int*)Unsafe.AsPointer(ref pinned.Buf[0]);
        for (int round = 0; round < 3; round++)
        {
            GCSettings.LargeObjectHeapCompactionMode = GCLargeObjectHeapCompactionMode.CompactOnce;
            GC.Collect(2, GCCollectionMode.Forced, blocking: true, compacting: true);
        }
        Assert.That((int*)Unsafe.AsPointer(ref pinned.Buf[0]) == before, "the address did not change");
        for (int i = 0; i < array.Length; i++) Assert.That(pinned.Buf[i], Is.EqualTo(i));
    }

    [Test]
    public void Allocate_SurvivesCompactingCollections()
    {
        using var heap = PinnedArray<int>.Allocate(300, zero: true, "poh");
        for (int i = 0; i < 300; i++) heap.Buf[i] = i ^ 0x55;
        GC.Collect(2, GCCollectionMode.Forced, blocking: true, compacting: true);
        for (int i = 0; i < 300; i++) Assert.That(heap.Buf[i], Is.EqualTo(i ^ 0x55));
    }

    [Test]
    public void Dispose_FreesTheHandle_AndIsIdempotent()
    {
        var pinned = PinnedArray<ushort>.Pin(new ushort[10], "x");
        Assert.That(pinned.HoldsHandle);
        pinned.Dispose();
        Assert.That(pinned.HoldsHandle, Is.False);
        Assert.DoesNotThrow(pinned.Dispose);
        var heap = PinnedArray<ushort>.Allocate(10, zero: true, "y");
        Assert.DoesNotThrow(heap.Dispose);
    }

    [MethodImpl(MethodImplOptions.NoInlining)]
    private static WeakReference PinThenDispose()
    {
        var array = new int[1000];
        var pinned = PinnedArray<int>.Pin(array, "to be disposed");
        var weak = new WeakReference(array);
        pinned.Dispose();
        return weak;
    }

    [MethodImpl(MethodImplOptions.NoInlining)]
    private static WeakReference PinThenForget()
    {
        var array = new int[1000];
        var weak = new WeakReference(array);
        _ = PinnedArray<int>.Pin(array, "to be finalized");
        return weak;
    }

    [Test]
    public void Dispose_ReleasesThePin_SoTheArrayCanBeCollected()
    {
        WeakReference weak = PinThenDispose();
        GC.Collect();
        GC.WaitForPendingFinalizers();
        GC.Collect();
        Assert.That(weak.IsAlive, Is.False, "a pin that Dispose did not free keeps the array alive for ever");
    }

    [Test]
    public void Finalizer_ReleasesThePin_OfAnObjectThatWasNeverDisposed()
    {
        WeakReference weak = PinThenForget();
        for (int i = 0; i < 3; i++)
        {
            GC.Collect();
            GC.WaitForPendingFinalizers();
        }
        Assert.That(weak.IsAlive, Is.False, "the finalizer of the owner frees the pin");
    }

    [Test]
    public void CheckedBuild_OutOfBoundsAccessIsReported()
    {
        if (!EsaBuild.BoundsChecked) Assert.Ignore("only the checked build (-p:EsaBounds=true) bounds-checks a Buf");
        using var pinned = PinnedArray<int>.Pin(new int[4], "four");
        Assert.Throws<EsaMemoryViolationException>(() => _ = pinned.Buf[4]);
        Assert.Throws<EsaMemoryViolationException>(() => _ = pinned.Buf[-1]);
    }
}
