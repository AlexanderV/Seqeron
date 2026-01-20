using System;
using NUnit.Framework;

namespace SuffixTree.Tests.Core
{
    /// <summary>
    /// Tests for diagnostic functionality (ToString, PrintTree, etc.)
    /// </summary>
    [TestFixture]
    public class DiagnosticsTests
    {
        #region ToString

        [Test]
        public void ToString_EmptyTree_DoesNotThrow()
        {
            var st = SuffixTree.Build("");

            Assert.DoesNotThrow(() => _ = st.ToString());
        }

        [Test]
        public void ToString_NonEmptyTree_ReturnsNonEmpty()
        {
            var st = SuffixTree.Build("hello");

            var result = st.ToString();

            Assert.That(result, Is.Not.Null.And.Not.Empty);
        }

        #endregion

        #region PrintTree

        [Test]
        public void PrintTree_EmptyTree_DoesNotThrow()
        {
            var st = SuffixTree.Build("");

            Assert.DoesNotThrow(() => _ = st.PrintTree());
        }

        [Test]
        public void PrintTree_NonEmptyTree_ReturnsNonEmpty()
        {
            var st = SuffixTree.Build("hello");

            var result = st.PrintTree();

            Assert.That(result, Is.Not.Null.And.Not.Empty);
        }

        #endregion

        #region Comparison

        [Test]
        public void PrintTree_DifferentFromToString()
        {
            var st = SuffixTree.Build("hello");

            var printResult = st.PrintTree();
            var toStringResult = st.ToString();

            // PrintTree is typically more detailed than ToString
            Assert.That(printResult.Length, Is.GreaterThanOrEqualTo(toStringResult.Length));
        }

        #endregion
    }
}
