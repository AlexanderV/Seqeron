using NUnit.Framework;
using SuffixTree;

namespace SuffixTree.Persistent.Tests;

/// <summary>
/// Runs the full <see cref="SuffixTreeTestBase"/> suite against the Large (v3, 64-bit offset)
/// storage format to verify functional parity with the default Compact (v4) format.
/// </summary>
[TestFixture]
public class LargeFormatSuffixTreeTests : SuffixTreeTestBase
{
    protected override ISuffixTree CreateTree(string text)
    {
        var storage = new HeapStorageProvider();
        var builder = new PersistentSuffixTreeBuilder(storage, NodeLayout.Large);
        var textSource = new StringTextSource(text);
        long root = builder.Build(textSource);
        return new PersistentSuffixTree(storage, root, textSource, NodeLayout.Large);
    }
}
