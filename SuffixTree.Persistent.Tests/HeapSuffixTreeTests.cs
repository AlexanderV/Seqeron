using System;
using NUnit.Framework;
using SuffixTree;

namespace SuffixTree.Persistent.Tests
{
    [TestFixture]
    public class HeapSuffixTreeTests : SuffixTreeTestBase
    {
        protected override ISuffixTree CreateTree(string text)
        {
            return PersistentSuffixTreeFactory.Create(text);
        }
    }
}
