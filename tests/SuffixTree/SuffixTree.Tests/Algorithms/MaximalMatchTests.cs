using System;
using System.Collections.Generic;
using System.Linq;
using NUnit.Framework;

namespace SuffixTree.Tests.Algorithms
{
    /// <summary>
    /// FindMaximalExactMatches / FindMaximalUniqueMatches — the MUMmer 3 match sets
    /// (Kurtz et al. 2004, Genome Biol 5:R12) on the forward strand.
    ///
    /// Expected values marked "MUMmer 3.23" were produced by the MUMmer 3.23 binary
    /// (Ubuntu package mummer 3.23+dfsg-8) with the command in the comment, converted from
    /// MUMmer's 1-based (ref, query, length) to 0-based (PositionInText, PositionInQuery, Length)
    /// and ordered by (query, text) position. MUMmer prints -maxmatch/-mumreference by query
    /// position and -mum by reference position; the library orders every set by (query, text).
    /// </summary>
    [TestFixture]
    [Category("Algorithms")]
    public class MaximalMatchTests
    {
        private static (int, int, int)[] Mems(string text, string query, int minLength)
            => SuffixTree.Build(text).FindMaximalExactMatches(query, minLength)
                .Select(m => (m.PositionInText, m.PositionInQuery, m.Length)).ToArray();

        private static (int, int, int)[] Mums(string text, string query, int minLength, MumUniqueness mode)
            => SuffixTree.Build(text).FindMaximalUniqueMatches(query, minLength, mode)
                .Select(m => (m.PositionInText, m.PositionInQuery, m.Length)).ToArray();

        #region MUMmer 3.23 literal outputs

        // mummer -maxmatch -l 3 ref.fa qry.fa   (ref GATTACAGATTACA, qry ATTACAGGATTACAT)
        //   2 1 7 / 9 1 6 / 1 8 7 / 8 8 7
        [Test]
        public void Maxmatch_GattacaRepeat_MatchesMummer()
        {
            Assert.That(Mems("GATTACAGATTACA", "ATTACAGGATTACAT", 3),
                Is.EqualTo(new[] { (1, 0, 7), (8, 0, 6), (0, 7, 7), (7, 7, 7) }));
        }

        // mummer -mum -l 3 and mummer -mumreference -l 3 (same inputs): 2 1 7
        [Test]
        public void Mum_GattacaRepeat_MatchesMummer()
        {
            Assert.That(Mums("GATTACAGATTACA", "ATTACAGGATTACAT", 3, MumUniqueness.Both),
                Is.EqualTo(new[] { (1, 0, 7) }));
            Assert.That(Mums("GATTACAGATTACA", "ATTACAGGATTACAT", 3, MumUniqueness.Reference),
                Is.EqualTo(new[] { (1, 0, 7) }));
        }

        // mummer -maxmatch -l 3 (ref ACGTACGTTACGT, qry TACGTACGTT): 4 1 5 / 9 1 5 / 1 2 9 / 9 5 5 / 1 6 4
        // mummer -mum -l 3 and -mumreference -l 3: 1 2 9
        [Test]
        public void TandemRepeat_MatchesMummer()
        {
            Assert.That(Mems("ACGTACGTTACGT", "TACGTACGTT", 3),
                Is.EqualTo(new[] { (3, 0, 5), (8, 0, 5), (0, 1, 9), (8, 4, 5), (0, 5, 4) }));
            Assert.That(Mums("ACGTACGTTACGT", "TACGTACGTT", 3, MumUniqueness.Both), Is.EqualTo(new[] { (0, 1, 9) }));
            Assert.That(Mums("ACGTACGTTACGT", "TACGTACGTT", 3, MumUniqueness.Reference), Is.EqualTo(new[] { (0, 1, 9) }));
        }

        // ref TTTCCTCATGCAATTCAAAACCATGTCCGTAATGTAGGCG, qry TTTCCTCATGCAAATAGCCTCATGCATCCGTAATGTAG
        // mummer -maxmatch -l 4:     1 1 13 / 22 7 4 / 16 11 4 / 4 18 9 / 22 21 4 / 26 27 12 / 23 33 4
        // mummer -mumreference -l 4: 1 1 13 / 16 11 4 / 4 18 9 / 26 27 12
        // mummer -mum -l 4:          1 1 13 / 16 11 4 / 26 27 12
        [Test]
        public void ThreeSetsDiffer_MatchesMummer()
        {
            const string reference = "TTTCCTCATGCAATTCAAAACCATGTCCGTAATGTAGGCG";
            const string query = "TTTCCTCATGCAAATAGCCTCATGCATCCGTAATGTAG";
            Assert.That(Mems(reference, query, 4), Is.EqualTo(new[]
                { (0, 0, 13), (21, 6, 4), (15, 10, 4), (3, 17, 9), (21, 20, 4), (25, 26, 12), (22, 32, 4) }));
            Assert.That(Mums(reference, query, 4, MumUniqueness.Reference), Is.EqualTo(new[]
                { (0, 0, 13), (15, 10, 4), (3, 17, 9), (25, 26, 12) }));
            Assert.That(Mums(reference, query, 4, MumUniqueness.Both), Is.EqualTo(new[]
                { (0, 0, 13), (15, 10, 4), (25, 26, 12) }));
        }

        // ref ACGT, qry TTAG, -l 1.
        // mummer -maxmatch -l 1 and -mumreference -l 1: 4 1 1 / 4 2 1 / 1 3 1 / 3 4 1
        // mummer -mum -l 1: 3 4 1 only. MUMmer's mumuniqueinquery (cleanMUMcand.c) starts its sweep
        // with dbright = 0, so a length-1 candidate at reference position 0 (right end 0) is dropped;
        // "A" (ref 0, qry 2) occurs once in each string and is a MUM by definition — kept here.
        [Test]
        public void LengthOne_MumAtReferenceStart_KeptDespiteMummerSweepArtefact()
        {
            Assert.That(Mems("ACGT", "TTAG", 1), Is.EqualTo(new[] { (3, 0, 1), (3, 1, 1), (0, 2, 1), (2, 3, 1) }));
            Assert.That(Mums("ACGT", "TTAG", 1, MumUniqueness.Reference),
                Is.EqualTo(new[] { (3, 0, 1), (3, 1, 1), (0, 2, 1), (2, 3, 1) }));
            Assert.That(Mums("ACGT", "TTAG", 1, MumUniqueness.Both), Is.EqualTo(new[] { (0, 2, 1), (2, 3, 1) }));
        }

        #endregion

        #region Definition cases

        [Test]
        public void Mems_ReportEveryReferenceOccurrence()
        {
            // "ACG" occurs at 0, 4 and 8; each occurrence is left-maximal (start / preceded by T)
            // and right-maximal (followed by T vs query end).
            Assert.That(Mems("ACGTACGTACGT", "ACG", 3), Is.EqualTo(new[] { (0, 0, 3), (4, 0, 3), (8, 0, 3) }));
            Assert.That(Mums("ACGTACGTACGT", "ACG", 3, MumUniqueness.Reference), Is.Empty);
        }

        [Test]
        public void Mems_HomopolymerOnlyBoundaryAnchoredMatches()
        {
            // Left-maximal needs r = 0 or q = 0; right-maximal needs an end of either string.
            // text AAAA, query AA, minLength 1: (r, 0, min(2, 4 - r)) for r = 0..3, and (0, 1, 1).
            Assert.That(Mems("AAAA", "AA", 1),
                Is.EqualTo(new[] { (0, 0, 2), (1, 0, 2), (2, 0, 2), (3, 0, 1), (0, 1, 1) }));
        }

        [Test]
        public void Mums_RepeatInQuery_ExcludedOnlyForBoth()
        {
            // "GATTACA" is unique in the text but occurs twice in the query.
            Assert.That(Mums("CCGATTACACC", "GATTACATTGATTACA", 5, MumUniqueness.Reference),
                Is.EqualTo(new[] { (2, 0, 7), (2, 9, 7) }));
            Assert.That(Mums("CCGATTACACC", "GATTACATTGATTACA", 5, MumUniqueness.Both), Is.Empty);
        }

        [Test]
        public void Mums_DefaultModeIsBoth()
        {
            var tree = SuffixTree.Build("CCGATTACACC");
            Assert.That(tree.FindMaximalUniqueMatches("GATTACATTGATTACA", 5), Is.Empty);
        }

        [Test]
        public void FindExactMatchAnchors_IsSubsetOfFullMemSet()
        {
            var tree = SuffixTree.Build("ACGTACGTTACGT");
            var all = tree.FindMaximalExactMatches("TACGTACGTT", 3).ToHashSet();
            Assert.That(tree.FindExactMatchAnchors("TACGTACGTT", 3), Is.SubsetOf(all));
        }

        #endregion

        #region Guards

        [Test]
        public void NullQuery_Throws()
        {
            var tree = SuffixTree.Build("banana");
            Assert.Throws<ArgumentNullException>(() => tree.FindMaximalExactMatches(null!, 3));
            Assert.Throws<ArgumentNullException>(() => tree.FindMaximalUniqueMatches(null!, 3));
        }

        [TestCase(0)]
        [TestCase(-1)]
        public void MinLengthBelowOne_Throws(int minLength)
        {
            var tree = SuffixTree.Build("banana");
            Assert.Throws<ArgumentOutOfRangeException>(() => tree.FindMaximalExactMatches("ana", minLength));
            Assert.Throws<ArgumentOutOfRangeException>(() => tree.FindMaximalUniqueMatches("ana", minLength));
        }

        [Test]
        public void UndefinedUniqueness_Throws()
        {
            var tree = SuffixTree.Build("banana");
            Assert.Throws<ArgumentOutOfRangeException>(() => tree.FindMaximalUniqueMatches("ana", 2, (MumUniqueness)7));
        }

        [Test]
        public void EmptyInputsOrMinLengthAboveLengths_ReturnEmpty()
        {
            Assert.That(SuffixTree.Build("banana").FindMaximalExactMatches("", 1), Is.Empty);
            Assert.That(SuffixTree.Build("").FindMaximalExactMatches("banana", 1), Is.Empty);
            Assert.That(SuffixTree.Build("").FindMaximalUniqueMatches("banana", 1, MumUniqueness.Reference), Is.Empty);
            Assert.That(SuffixTree.Build("banana").FindMaximalExactMatches("banana", 7), Is.Empty);
            Assert.That(SuffixTree.Build("banana").FindMaximalUniqueMatches("banana", 7), Is.Empty);
        }

        #endregion

        #region Seeded property test vs brute force

        [Test]
        public void RandomPairs_EqualBruteForceDefinition()
        {
            var rng = new Random(20260930);
            for (int iter = 0; iter < 150; iter++)
            {
                string alphabet = iter % 3 == 0 ? "AB" : "ACGT";
                string text = RandomString(rng, rng.Next(0, 60), alphabet);
                string query = rng.Next(2) == 0 || text.Length == 0
                    ? RandomString(rng, rng.Next(0, 60), alphabet)
                    : Mutate(rng, text.Substring(rng.Next(text.Length)), alphabet);
                int minLength = rng.Next(1, 6);

                var tree = SuffixTree.Build(text);
                var mems = BruteForceMems(text, query, minLength);
                string ctx = $"iter {iter}: text={text} query={query} L={minLength}";
                Assert.That(tree.FindMaximalExactMatches(query, minLength), Is.EqualTo(mems), ctx);
                Assert.That(tree.FindMaximalUniqueMatches(query, minLength, MumUniqueness.Reference),
                    Is.EqualTo(mems.Where(m => Count(text, query.Substring(m.Item2, m.Item3)) == 1).ToList()), ctx);
                Assert.That(tree.FindMaximalUniqueMatches(query, minLength, MumUniqueness.Both),
                    Is.EqualTo(mems.Where(m => Count(text, query.Substring(m.Item2, m.Item3)) == 1
                                            && Count(query, query.Substring(m.Item2, m.Item3)) == 1).ToList()), ctx);
            }
        }

        private static List<(int, int, int)> BruteForceMems(string text, string query, int minLength)
        {
            var result = new List<(int, int, int)>();
            for (int q = 0; q < query.Length; q++)
                for (int r = 0; r < text.Length; r++)
                {
                    if (q > 0 && r > 0 && text[r - 1] == query[q - 1]) continue;
                    int k = 0;
                    while (r + k < text.Length && q + k < query.Length && text[r + k] == query[q + k]) k++;
                    if (k >= minLength) result.Add((r, q, k));
                }
            return result; // ascending by (q, r)
        }

        private static int Count(string s, string p)
        {
            int c = 0;
            for (int i = 0; i + p.Length <= s.Length; i++)
                if (string.CompareOrdinal(s, i, p, 0, p.Length) == 0) c++;
            return c;
        }

        private static string RandomString(Random rng, int length, string alphabet)
        {
            var c = new char[length];
            for (int i = 0; i < length; i++) c[i] = alphabet[rng.Next(alphabet.Length)];
            return new string(c);
        }

        private static string Mutate(Random rng, string s, string alphabet)
        {
            var c = s.ToCharArray();
            for (int i = 0; i < c.Length; i++)
                if (rng.Next(10) == 0) c[i] = alphabet[rng.Next(alphabet.Length)];
            return new string(c);
        }

        #endregion
    }
}
