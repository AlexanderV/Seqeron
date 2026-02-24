using System;
using System.IO;
using System.Linq;
using NUnit.Framework;
using SuffixTree;

namespace SuffixTree.Persistent.Tests.Parity
{
    [TestFixture]
    public class UnsafeParityTests
    {
        private static readonly string[] TestStrings = new[]
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

        [Test]
        public void Parity_SafeVsUnsafe_FastPaths_Substrings([ValueSource(nameof(TestStrings))] string text)
        {
            string tempFile = Path.GetTempFileName();
            try
            {
                // Create tree in Heap (safe paths)
                using var safeTreeDisposable = PersistentSuffixTreeFactory.Create(new StringTextSource(text)) as IDisposable;
                var safeTree = (ISuffixTree)safeTreeDisposable!;

                // Create tree in MappedFile (unsafe/fast paths with raw pointers)
                using var unsafeTreeDisposable = PersistentSuffixTreeFactory.Create(new StringTextSource(text), tempFile) as IDisposable;
                var unsafeTree = (ISuffixTree)unsafeTreeDisposable!;

                Assert.Multiple(() =>
                {
                    Assert.That(unsafeTree.LeafCount, Is.EqualTo(safeTree.LeafCount), "LeafCount mismatch");
                    
                    if (text.Length > 0)
                    {
                        var random = new Random(42);
                        for (int i = 0; i < 20; i++)
                        {
                            int start = random.Next(text.Length);
                            int len = random.Next(1, text.Length - start + 1);
                            string sub = text.Substring(start, len);

                            // Test TryGetChild vs TryGetChildFast (via Contains/CountOccurrences)
                            Assert.That(unsafeTree.Contains(sub), Is.EqualTo(safeTree.Contains(sub)), $"Contains mismatch for '{sub}'");
                            Assert.That(unsafeTree.CountOccurrences(sub), Is.EqualTo(safeTree.CountOccurrences(sub)), $"CountOccurrences mismatch for '{sub}'");

                            // Test CollectLeavesSequential vs CollectLeavesUnsafe (via FindAllOccurrences)
                            var safeOccs = safeTree.FindAllOccurrences(sub).ToList();
                            var unsafeOccs = unsafeTree.FindAllOccurrences(sub).ToList();
                            
                            safeOccs.Sort();
                            unsafeOccs.Sort();
                            
                            Assert.That(unsafeOccs, Is.EqualTo(safeOccs), $"FindAllOccurrences mismatch for '{sub}'");
                        }
                        
                        // Also test a pattern that does not exist
                        string missing = text + "X!";
                        Assert.That(unsafeTree.Contains(missing), Is.EqualTo(safeTree.Contains(missing)), "Contains mismatch for missing pattern");
                        Assert.That(unsafeTree.FindAllOccurrences(missing), Is.EqualTo(safeTree.FindAllOccurrences(missing)), "FindAllOccurrences mismatch for missing pattern");
                    }
                });
            }
            finally
            {
                if (File.Exists(tempFile))
                {
                    try { File.Delete(tempFile); } catch { }
                }
                // Cleanup temporary files created by PersistentSuffixTreeFactory that use original file path prefix
                var tmpChild = tempFile + ".children.tmp";
                if (File.Exists(tmpChild)) try { File.Delete(tmpChild); } catch { }
                var tmpDepth = tempFile + ".depth.tmp";
                if (File.Exists(tmpDepth)) try { File.Delete(tmpDepth); } catch { }
            }
        }
        
        [Test]
        public void Parity_SafeVsUnsafe_RandomStrings()
        {
            var random = new Random(1337);
            for (int i = 0; i < 10; i++)
            {
                int len = random.Next(50, 500);
                string text = new string(Enumerable.Range(0, len).Select(_ => (char)random.Next('a', 'e')).ToArray());

                string tempFile = Path.GetTempFileName();
                try
                {
                    using var safeTreeDisposable = PersistentSuffixTreeFactory.Create(new StringTextSource(text)) as IDisposable;
                    var safeTree = (ISuffixTree)safeTreeDisposable!;

                    using var unsafeTreeDisposable = PersistentSuffixTreeFactory.Create(new StringTextSource(text), tempFile) as IDisposable;
                    var unsafeTree = (ISuffixTree)unsafeTreeDisposable!;

                    string pattern = text.Substring(random.Next(len / 2), 5);
                    Assert.That(unsafeTree.Contains(pattern), Is.EqualTo(safeTree.Contains(pattern)), $"Contains mismatch for '{pattern}' on random string {i}");
                    
                    var safeOccs = safeTree.FindAllOccurrences(pattern).ToList();
                    var unsafeOccs = unsafeTree.FindAllOccurrences(pattern).ToList();
                    safeOccs.Sort();
                    unsafeOccs.Sort();
                    Assert.That(unsafeOccs, Is.EqualTo(safeOccs), $"FindAllOccurrences mismatch for '{pattern}' on random string {i}");
                }
                finally
                {
                    if (File.Exists(tempFile)) try { File.Delete(tempFile); } catch { }
                    var tmpChild = tempFile + ".children.tmp";
                    if (File.Exists(tmpChild)) try { File.Delete(tmpChild); } catch { }
                    var tmpDepth = tempFile + ".depth.tmp";
                    if (File.Exists(tmpDepth)) try { File.Delete(tmpDepth); } catch { }
                }
            }
        }
    }
}
