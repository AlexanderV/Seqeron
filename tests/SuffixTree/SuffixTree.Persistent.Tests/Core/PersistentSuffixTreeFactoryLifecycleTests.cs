using System;
using System.IO;
using NUnit.Framework;
using SuffixTree;
using SuffixTree.Persistent;

namespace SuffixTree.Persistent.Tests.Core;

[TestFixture]
[Category("Core")]
public class PersistentSuffixTreeFactoryLifecycleTests
{
    [Test]
    public void CreatePersistent_UsesConcreteDisposableType()
    {
        using var tree = PersistentSuffixTreeFactory.CreatePersistent(new StringTextSource("banana"));

        Assert.Multiple(() =>
        {
            Assert.That(tree, Is.TypeOf<PersistentSuffixTree>());
            Assert.That(tree.Contains("ana"), Is.True);
            Assert.That(tree.CountOccurrences("na"), Is.EqualTo(2));
        });
    }

    [Test]
    public void LoadPersistent_UsesConcreteDisposableType()
    {
        string filePath = Path.Combine(Path.GetTempPath(), $"ST_{Guid.NewGuid():N}.tree");
        try
        {
            using (var created = PersistentSuffixTreeFactory.CreatePersistent(new StringTextSource("abracadabra"), filePath))
            {
                Assert.That(created.Contains("abra"), Is.True);
            }

            using var loaded = PersistentSuffixTreeFactory.LoadPersistent(filePath);

            Assert.Multiple(() =>
            {
                Assert.That(loaded, Is.TypeOf<PersistentSuffixTree>());
                Assert.That(loaded.Contains("cada"), Is.True);
                Assert.That(loaded.CountOccurrences("a"), Is.EqualTo(5));
            });
        }
        finally
        {
            try { File.Delete(filePath); } catch { }
        }
    }

    [Test]
    public void CreatePersistent_NullText_ThrowsArgumentNullException()
    {
        Assert.Throws<ArgumentNullException>(() => PersistentSuffixTreeFactory.CreatePersistent(null!));
    }

    [Test]
    public void LoadPersistent_EmptyPath_ThrowsArgumentException()
    {
        Assert.Throws<ArgumentException>(() => PersistentSuffixTreeFactory.LoadPersistent(string.Empty));
    }
}
