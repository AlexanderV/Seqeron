using System;
using System.Collections.Generic;
using System.Linq;
using NUnit.Framework;
using SuffixTree;

namespace SuffixTree.Persistent.Tests
{
    [TestFixture]
    public class ParityTests
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
        public void Parity_WithReferenceImplementation_AllMethods([ValueSource(nameof(TestStrings))] string text)
        {
            var reference = global::SuffixTree.SuffixTree.Build(text);
            using (var persistentDisposable = PersistentSuffixTreeFactory.Create(text) as IDisposable)
            {
                var persistent = (ISuffixTree)persistentDisposable!;

                Assert.Multiple(() =>
                {
                    // 1. Text parity
                    Assert.That(persistent.Text, Is.EqualTo(reference.Text), "Text mismatch");

                    // 2. Node/Leaf count parity
                    // Note: Node counts might slightly differ if builders use different splitting strategies 
                    // (though Ukkonen should be consistent), but LeafCount must be identical.
                    Assert.That(persistent.LeafCount, Is.EqualTo(reference.LeafCount), "Total LeafCount mismatch");

                    // 3. Substring existence (random samples from text)
                    if (text.Length > 0)
                    {
                        var random = new Random(42);
                        for (int i = 0; i < 20; i++)
                        {
                            int start = random.Next(text.Length);
                            int len = random.Next(1, text.Length - start + 1);
                            string sub = text.Substring(start, len);

                            Assert.That(persistent.Contains(sub), Is.True, $"Contains failed for valid substring '{sub}'");
                            Assert.That(persistent.CountOccurrences(sub), Is.EqualTo(reference.CountOccurrences(sub)), $"CountOccurrences mismatch for '{sub}'");
                            
                            var refOccs = reference.FindAllOccurrences(sub).ToList();
                            var perOccs = persistent.FindAllOccurrences(sub).ToList();
                            refOccs.Sort();
                            perOccs.Sort();
                            Assert.That(perOccs, Is.EqualTo(refOccs), $"FindAllOccurrences mismatch for '{sub}'");
                        }
                    }

                    // 4. Longest Repeated Substring
                    Assert.That(persistent.LongestRepeatedSubstring(), Is.EqualTo(reference.LongestRepeatedSubstring()), "LRS mismatch");

                    // 5. Longest Common Substring with another string
                    string other = "some_other_shared_stuff_ana_missi";
                    Assert.That(persistent.LongestCommonSubstring(other), Is.EqualTo(reference.LongestCommonSubstring(other)), "LCS mismatch");
                });
            }
        }

        [Test]
        public void Parity_RandomStrings_Comprehensive()
        {
            var random = new Random(1337);
            for (int i = 0; i < 10; i++)
            {
                int len = random.Next(50, 500);
                string text = new string(Enumerable.Range(0, len).Select(_ => (char)random.Next('a', 'e')).ToArray());
                
                var reference = global::SuffixTree.SuffixTree.Build(text);
                using (var persistentDisposable = PersistentSuffixTreeFactory.Create(text) as IDisposable)
                {
                    var persistent = (ISuffixTree)persistentDisposable!;

                    Assert.That(persistent.LeafCount, Is.EqualTo(reference.LeafCount), $"LeafCount mismatch for random string index {i}");
                    
                    // Test some patterns
                    string pattern = text.Substring(random.Next(len/2), 5);
                    Assert.That(persistent.CountOccurrences(pattern), Is.EqualTo(reference.CountOccurrences(pattern)));
                    
                    Assert.That(persistent.LongestRepeatedSubstring(), Is.EqualTo(reference.LongestRepeatedSubstring()));
                }
            }
        }
    }
}
