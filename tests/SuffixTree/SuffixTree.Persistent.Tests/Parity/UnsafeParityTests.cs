using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using NUnit.Framework;
using SuffixTree;

namespace SuffixTree.Persistent.Tests.Parity
{
    [TestFixture]
    [Category("Parity")]
    public class UnsafeParityTests
    {
        private static readonly string[] TestStrings =
        {
            "banana",
            "mississippi",
            "abracadabra",
            "aaaaaaaaaa",
            "abcdefghij",
            "a",
            "",
            "random123!@#",
            "🧬αβγ🧪$",
            "repetitive-repetitive-repetitive"
        };

        private static void AssertQueryParity(ISuffixTree safeTree, ISuffixTree unsafeTree, string pattern, string context)
        {
            Assert.That(unsafeTree.Contains(pattern), Is.EqualTo(safeTree.Contains(pattern)),
                $"Contains mismatch for pattern '{pattern}' ({context})");
            Assert.That(unsafeTree.CountOccurrences(pattern), Is.EqualTo(safeTree.CountOccurrences(pattern)),
                $"CountOccurrences mismatch for pattern '{pattern}' ({context})");

            var safeOccs = safeTree.FindAllOccurrences(pattern).OrderBy(x => x).ToList();
            var unsafeOccs = unsafeTree.FindAllOccurrences(pattern).OrderBy(x => x).ToList();
            Assert.That(unsafeOccs, Is.EqualTo(safeOccs),
                $"FindAllOccurrences mismatch for pattern '{pattern}' ({context})");
        }

        private static string BuildMissingPattern(string text)
        {
            string[] candidates = { "\u0000", "\uFFFF", "§§§", text + "\u0000", text + "\uFFFF" };
            foreach (var candidate in candidates)
            {
                if (!text.Contains(candidate, StringComparison.Ordinal))
                    return candidate;
            }

            return text + "|missing|";
        }

        private static void CleanupTempFiles(string tempFile)
        {
            if (File.Exists(tempFile))
            {
                try { File.Delete(tempFile); } catch { }
            }

            var tmpChild = tempFile + ".children.tmp";
            if (File.Exists(tmpChild)) try { File.Delete(tmpChild); } catch { }

            var tmpDepth = tempFile + ".depth.tmp";
            if (File.Exists(tmpDepth)) try { File.Delete(tmpDepth); } catch { }
        }

        [Test]
        public void Parity_SafeVsUnsafe_FastPaths_Substrings([ValueSource(nameof(TestStrings))] string text)
        {
            string tempFile = Path.GetTempFileName();
            try
            {
                using var safeTreeDisposable = PersistentSuffixTreeFactory.Create(new StringTextSource(text)) as IDisposable;
                var safeTree = (ISuffixTree)safeTreeDisposable!;

                using var unsafeTreeDisposable = PersistentSuffixTreeFactory.Create(new StringTextSource(text), tempFile) as IDisposable;
                var unsafeTree = (ISuffixTree)unsafeTreeDisposable!;

                Assert.Multiple(() =>
                {
                    Assert.That(unsafeTree.LeafCount, Is.EqualTo(safeTree.LeafCount), "LeafCount mismatch");
                    Assert.That(unsafeTree.NodeCount, Is.EqualTo(safeTree.NodeCount), "NodeCount mismatch");
                    Assert.That(unsafeTree.MaxDepth, Is.EqualTo(safeTree.MaxDepth), "MaxDepth mismatch");
                    Assert.That(unsafeTree.IsEmpty, Is.EqualTo(safeTree.IsEmpty), "IsEmpty mismatch");

                    AssertQueryParity(safeTree, unsafeTree, string.Empty, "empty pattern contract");
                    AssertQueryParity(safeTree, unsafeTree, BuildMissingPattern(text), "missing pattern");

                    if (!string.IsNullOrEmpty(text))
                    {
                        AssertQueryParity(safeTree, unsafeTree, text, "full text");
                        AssertQueryParity(safeTree, unsafeTree, text.Substring(0, 1), "first char");
                        AssertQueryParity(safeTree, unsafeTree, text.Substring(text.Length - 1, 1), "last char");
                        AssertQueryParity(safeTree, unsafeTree, text.Substring(0, Math.Min(8, text.Length)), "prefix");
                        AssertQueryParity(safeTree, unsafeTree, text.Substring(Math.Max(0, text.Length - Math.Min(8, text.Length))), "suffix");

                        var random = new Random(42);
                        for (int i = 0; i < 30; i++)
                        {
                            int start = random.Next(text.Length);
                            int len = random.Next(1, text.Length - start + 1);
                            string sub = text.Substring(start, len);
                            AssertQueryParity(safeTree, unsafeTree, sub, $"random substring #{i}");
                        }
                    }
                });
            }
            finally
            {
                CleanupTempFiles(tempFile);
            }
        }

        [Test]
        public void Parity_SafeVsUnsafe_ControlAndUnicodePatterns()
        {
            const string text = "A\u0000B\uFFFFC🧬αβγ🧪$XYZXYZXYZ";
            string tempFile = Path.GetTempFileName();

            try
            {
                using var safeTreeDisposable = PersistentSuffixTreeFactory.Create(new StringTextSource(text)) as IDisposable;
                var safeTree = (ISuffixTree)safeTreeDisposable!;

                using var unsafeTreeDisposable = PersistentSuffixTreeFactory.Create(new StringTextSource(text), tempFile) as IDisposable;
                var unsafeTree = (ISuffixTree)unsafeTreeDisposable!;

                var patterns = new List<string>
                {
                    "\u0000",
                    "\uFFFF",
                    "🧬αβγ",
                    "XYZXYZ",
                    "A\u0000B",
                    "βγ🧪$",
                    "not-present"
                };

                foreach (string pattern in patterns)
                {
                    AssertQueryParity(safeTree, unsafeTree, pattern, "control/unicode");
                }
            }
            finally
            {
                CleanupTempFiles(tempFile);
            }
        }
    }
}
