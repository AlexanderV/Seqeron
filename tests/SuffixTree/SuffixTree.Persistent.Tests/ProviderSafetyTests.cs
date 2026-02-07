using System;
using System.IO;
using NUnit.Framework;
using SuffixTree.Persistent;

namespace SuffixTree.Persistent.Tests;

/// <summary>
/// Tests for P5 (ObjectDisposedException instead of NRE after Dispose),
/// P9 (Allocate(negative) must not corrupt state).
/// P2 (EnsureCapacity/TrimToSize exception safety) is fixed in GREEN phase.
/// Written RED-first: these tests expose existing bugs.
/// </summary>
[TestFixture]
public class ProviderSafetyTests
{
    private string _tempFile = string.Empty;

    [SetUp]
    public void SetUp()
    {
        _tempFile = Path.GetTempFileName();
    }

    [TearDown]
    public void TearDown()
    {
        try { if (File.Exists(_tempFile)) File.Delete(_tempFile); } catch { }
    }

    // ─── P5: HeapStorageProvider must throw ObjectDisposedException ──

    [Test]
    public void HeapProvider_ReadInt32_AfterDispose_ThrowsODE()
    {
        var p = new HeapStorageProvider();
        p.Allocate(4);
        p.WriteInt32(0, 1);
        p.Dispose();
        Assert.Throws<ObjectDisposedException>(() => p.ReadInt32(0));
    }

    [Test]
    public void HeapProvider_WriteInt32_AfterDispose_ThrowsODE()
    {
        var p = new HeapStorageProvider();
        p.Allocate(4);
        p.Dispose();
        Assert.Throws<ObjectDisposedException>(() => p.WriteInt32(0, 1));
    }

    [Test]
    public void HeapProvider_ReadBytes_AfterDispose_ThrowsODE()
    {
        var p = new HeapStorageProvider();
        p.Allocate(4);
        p.Dispose();
        Assert.Throws<ObjectDisposedException>(() => p.ReadBytes(0, new byte[4], 0, 4));
    }

    [Test]
    public void HeapProvider_WriteBytes_AfterDispose_ThrowsODE()
    {
        var p = new HeapStorageProvider();
        p.Allocate(4);
        p.Dispose();
        Assert.Throws<ObjectDisposedException>(() => p.WriteBytes(0, new byte[4], 0, 4));
    }

    [Test]
    public void HeapProvider_Allocate_AfterDispose_ThrowsODE()
    {
        var p = new HeapStorageProvider();
        p.Dispose();
        Assert.Throws<ObjectDisposedException>(() => p.Allocate(4));
    }

    // ─── P5: MappedFileStorageProvider must throw ObjectDisposedException ──

    [Test]
    public void MappedProvider_ReadInt32_AfterDispose_ThrowsODE()
    {
        var p = new MappedFileStorageProvider(_tempFile);
        p.Allocate(4);
        p.WriteInt32(0, 1);
        p.Dispose();
        Assert.Throws<ObjectDisposedException>(() => p.ReadInt32(0));
    }

    [Test]
    public void MappedProvider_WriteInt32_AfterDispose_ThrowsODE()
    {
        var p = new MappedFileStorageProvider(_tempFile);
        p.Allocate(4);
        p.Dispose();
        Assert.Throws<ObjectDisposedException>(() => p.WriteInt32(0, 1));
    }

    [Test]
    public void MappedProvider_ReadBytes_AfterDispose_ThrowsODE()
    {
        var p = new MappedFileStorageProvider(_tempFile);
        p.Allocate(4);
        p.Dispose();
        Assert.Throws<ObjectDisposedException>(() => p.ReadBytes(0, new byte[4], 0, 4));
    }

    [Test]
    public void MappedProvider_WriteBytes_AfterDispose_ThrowsODE()
    {
        var p = new MappedFileStorageProvider(_tempFile);
        p.Allocate(4);
        p.Dispose();
        Assert.Throws<ObjectDisposedException>(() => p.WriteBytes(0, new byte[4], 0, 4));
    }

    // ─── P9: Allocate(negative) must reject and preserve state ──────

    [Test]
    public void HeapProvider_AllocateNegative_Throws()
    {
        var p = new HeapStorageProvider();
        p.Allocate(8);
        Assert.Throws<ArgumentOutOfRangeException>(() => p.Allocate(-1));
        // Position must remain 8 — not corrupted
        Assert.That(p.Size, Is.EqualTo(8));
        p.Dispose();
    }

    [Test]
    public void MappedProvider_AllocateNegative_Throws()
    {
        var p = new MappedFileStorageProvider(_tempFile);
        p.Allocate(8);
        Assert.Throws<ArgumentOutOfRangeException>(() => p.Allocate(-1));
        Assert.That(p.Size, Is.EqualTo(8));
        p.Dispose();
    }

    // ─── P2: EnsureCapacity / TrimToSize preserve data (sanity) ─────

    [Test]
    public void MappedProvider_EnsureCapacity_PreservesData()
    {
        var p = new MappedFileStorageProvider(_tempFile, initialCapacity: 128);
        p.Allocate(4);
        p.WriteInt32(0, 0xDEAD);

        // Force expansion beyond initial capacity
        p.EnsureCapacity(1024);

        Assert.That(p.ReadInt32(0), Is.EqualTo(0xDEAD));
        p.Dispose();
    }

    [Test]
    public void MappedProvider_TrimToSize_PreservesData()
    {
        var p = new MappedFileStorageProvider(_tempFile, initialCapacity: 65536);
        p.Allocate(4);
        p.WriteInt32(0, 0xBEEF);

        p.TrimToSize();

        Assert.That(p.ReadInt32(0), Is.EqualTo(0xBEEF));
        p.Dispose();
    }

    // ─── S15: Allocate must not corrupt _position if EnsureCapacity throws ───

    [Test]
    public void HeapProvider_Allocate_FailedExpansion_DoesNotCorruptPosition()
    {
        // HeapStorageProvider caps at int.MaxValue.
        // Allocate a small block, then try to allocate a block that would push
        // _position beyond int.MaxValue.  The allocation must fail AND Size must
        // remain at the value it had before the failed call.
        var p = new HeapStorageProvider(initialCapacity: 128);
        p.Allocate(64);
        long sizeBefore = p.Size;
        Assert.That(sizeBefore, Is.EqualTo(64));

        // Request that would push _position beyond 2 GB
        Assert.Throws<InvalidOperationException>(() => p.Allocate(int.MaxValue));

        // Critical: Size must be unchanged after failed allocation
        Assert.That(p.Size, Is.EqualTo(sizeBefore),
            "S15: Failed Allocate must not corrupt _position");
        p.Dispose();
    }
}
