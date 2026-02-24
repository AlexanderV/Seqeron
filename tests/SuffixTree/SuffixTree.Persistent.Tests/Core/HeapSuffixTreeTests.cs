using System;
using NUnit.Framework;
using SuffixTree;

namespace SuffixTree.Persistent.Tests.Core
{
    [TestFixture]
    [Category("Core")]
    public class HeapSuffixTreeTests : SuffixTreeTestBase
    {
        protected override ISuffixTree CreateTree(string text)
        {
            return PersistentSuffixTreeFactory.Create(new StringTextSource(text));
        }
    }
}

