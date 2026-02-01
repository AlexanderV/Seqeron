using System;

namespace SuffixTree.Console
{
    class Program
    {
        static void Main(string[] args)
        {
            // Extensive stress test
            var r = new Random(12345);
            int totalTests = 0;
            int totalSubstrings = 0;

            System.Console.WriteLine("Running exhaustive stress tests...");

            // Test 1: Small strings with small alphabet (exhaustive)
            System.Console.WriteLine("\n[Phase 1] Small strings, small alphabet (exhaustive)");
            for (int length = 1; length <= 50; length++)
            {
                for (int trial = 0; trial < 5; trial++)
                {
                    var chars = new char[length];
                    for (int i = 0; i < length; i++)
                        chars[i] = (char)('a' + r.Next(4));

                    var s = new string(chars);
                    var tree = SuffixTree.Build(s);

                    // Verify ALL substrings
                    for (int i = 0; i < s.Length; i++)
                    {
                        for (int len = 1; len <= s.Length - i; len++)
                        {
                            var substr = s.Substring(i, len);
                            if (!tree.Contains(substr))
                            {
                                System.Console.WriteLine($"FAIL! String '{s}': Substring '{substr}' not found!");
                                return;
                            }
                            totalSubstrings++;
                        }
                    }
                    totalTests++;
                }
            }
            System.Console.WriteLine($"  Phase 1 complete: {totalTests} tests, {totalSubstrings:N0} substrings");

            // Test 2: Large strings with large alphabet (spot check)
            System.Console.WriteLine("\n[Phase 2] Large strings, large alphabet (spot check)");
            const string FULL_ALPHABET = "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789";
            for (int length = 100; length <= 5000; length += 100)
            {
                var chars = new char[length];
                for (int i = 0; i < length; i++)
                    chars[i] = FULL_ALPHABET[r.Next(FULL_ALPHABET.Length)];

                var s = new string(chars);
                var tree = SuffixTree.Build(s);

                // Spot check: verify 500 random substrings
                for (int check = 0; check < 500; check++)
                {
                    int start = r.Next(s.Length);
                    int len = r.Next(1, Math.Min(50, s.Length - start + 1));
                    var substr = s.Substring(start, len);
                    if (!tree.Contains(substr))
                    {
                        System.Console.WriteLine($"FAIL! Length {length}: Substring '{substr}' not found!");
                        return;
                    }
                    totalSubstrings++;
                }

                // Verify all suffixes
                for (int i = 0; i < s.Length; i++)
                {
                    var suffix = s.Substring(i);
                    if (!tree.Contains(suffix))
                    {
                        System.Console.WriteLine($"FAIL! Length {length}: Suffix at {i} not found!");
                        return;
                    }
                }

                totalTests++;
                if (length % 1000 == 0)
                    System.Console.WriteLine($"  Length {length} complete...");
            }

            // Test 3: Edge cases
            System.Console.WriteLine("\n[Phase 3] Edge cases");
            string[] edgeCases = {
                "a",
                "aa",
                "ab",
                "aaa",
                "aba",
                "abab",
                "abcabc",
                new string('a', 1000),
                new string('a', 1000) + "b",
                "abcabxabcd", // Classic Ukkonen test
                "mississippi",
                "abracadabra",
                "aaaabaaaabaaaabaaa"
            };

            foreach (var s in edgeCases)
            {
                var tree = SuffixTree.Build(s);
                var display = s.Length > 30 ? s.Substring(0, 27) + "..." : s;

                // Verify ALL substrings
                for (int i = 0; i < s.Length; i++)
                {
                    for (int len = 1; len <= s.Length - i; len++)
                    {
                        var substr = s.Substring(i, len);
                        if (!tree.Contains(substr))
                        {
                            System.Console.WriteLine($"FAIL! Edge case '{display}': Substring '{substr}' not found!");
                            return;
                        }
                        totalSubstrings++;
                    }
                }
                totalTests++;
                System.Console.WriteLine($"  '{display}' OK");
            }

            System.Console.WriteLine();
            System.Console.WriteLine($"=== ALL TESTS PASSED ===");
            System.Console.WriteLine($"Total tests: {totalTests}");
            System.Console.WriteLine($"Total substrings verified: {totalSubstrings:N0}");
        }
    }
}