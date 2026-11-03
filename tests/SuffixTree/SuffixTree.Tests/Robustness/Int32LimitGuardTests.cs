using System;
using System.Collections;
using System.Collections.Generic;
using NUnit.Framework;

namespace SuffixTree.Tests.Robustness
{
    /// <summary>
    /// ST-S0 WP 0.2: a text whose worst-case node count (2n + 1) does not fit in Int32 must be rejected
    /// with a clear exception before any work is done. The oversized texts are fake <see cref="ITextSource"/>
    /// instances that report a huge length and refuse to serve content, so nothing large is allocated.
    /// </summary>
    [TestFixture]
    [Category("Robustness")]
    public class Int32LimitGuardTests
    {
        private const int SmallestRejectedLength = 1_073_741_824; // 2n + 1 = 2^31 + 1 > int.MaxValue

        private sealed class HugeLengthTextSource : ITextSource
        {
            public HugeLengthTextSource(int length) { Length = length; }
            public int Length { get; }
            public char this[int index] => throw new NotSupportedException("content must not be read");
            public string Substring(int start, int length) => throw new NotSupportedException("content must not be read");
            public ReadOnlySpan<char> Slice(int start, int length) => throw new NotSupportedException("content must not be read");
            public IEnumerator<char> GetEnumerator() => throw new NotSupportedException("content must not be read");
            IEnumerator IEnumerable.GetEnumerator() => throw new NotSupportedException("content must not be read");
        }

        [Test]
        public void Build_TextOneBeyondLimit_ThrowsArgumentException_BeforeReadingContent()
        {
            var ex = Assert.Throws<ArgumentException>(() => SuffixTree.Build(new HugeLengthTextSource(SmallestRejectedLength)));

            Assert.That(ex!.ParamName, Is.EqualTo("value"));
            Assert.That(ex.Message, Does.Contain("1073741824").And.Contain("1073741823"));
        }

        [Test]
        public void Build_TextExactlyAtLimit_PassesTheLengthGuard()
        {
            // The largest legal length must not be rejected: construction starts and only then trips over
            // the fake source's refusal to serve content.
            Assert.Throws<NotSupportedException>(() => SuffixTree.Build(new HugeLengthTextSource(SuffixTreeLimits.MaxTextLength)));
        }

        [Test]
        public void MaxTextLength_IsTheLargestLengthWhoseWorstCaseNodeCountFitsInt32()
        {
            Assert.Multiple(() =>
            {
                Assert.That(SuffixTreeLimits.MaxTextLength, Is.EqualTo(1_073_741_823));
                Assert.That(2L * SuffixTreeLimits.MaxTextLength + 1, Is.LessThanOrEqualTo(int.MaxValue));
                Assert.That(2L * (SuffixTreeLimits.MaxTextLength + 1L) + 1, Is.GreaterThan(int.MaxValue));
            });
        }

        [Test]
        public void Build_TextOfIntMaxValue_ThrowsArgumentException()
        {
            Assert.Throws<ArgumentException>(() => SuffixTree.Build(new HugeLengthTextSource(int.MaxValue)));
        }

        [Test]
        public void TryBuild_TextBeyondLimit_ReturnsFalse_DoesNotThrow()
        {
            bool ok = false;
            SuffixTree? tree = null;

            Assert.DoesNotThrow(() => ok = SuffixTree.TryBuild(new HugeLengthTextSource(SmallestRejectedLength), out tree));

            Assert.That(ok, Is.False);
            Assert.That(tree, Is.Null);
        }
    }
}
