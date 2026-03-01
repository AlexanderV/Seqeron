using NUnit.Framework;

namespace Seqeron.Genomics.Tests
{
    /// <summary>
    /// Test suite for PAT-APPROX-002: Approximate Matching (Edit Distance).
    /// 
    /// Tests EditDistance (Levenshtein distance) and FindWithEdits methods.
    /// Evidence sources: Wikipedia, Rosetta Code, Navarro (2001).
    /// 
    /// Canonical test vectors from:
    /// - Wikipedia: "kitten" → "sitting" = 3, "flaw" → "lawn" = 2
    /// - Rosetta Code: "rosettacode" → "raisethysword" = 8
    /// </summary>
    [TestFixture]
    public class ApproximateMatcher_EditDistance_Tests
    {
        #region EditDistance - MUST Tests (Evidence-Backed)

        [Test]
        [Description("M01: Identity property - identical strings have distance 0")]
        public void EditDistance_IdenticalStrings_ReturnsZero()
        {
            Assert.Multiple(() =>
            {
                Assert.That(ApproximateMatcher.EditDistance("ACGT", "ACGT"), Is.EqualTo(0));
                Assert.That(ApproximateMatcher.EditDistance("kitten", "kitten"), Is.EqualTo(0));
                Assert.That(ApproximateMatcher.EditDistance("a", "a"), Is.EqualTo(0));
            });
        }

        [Test]
        [Description("M02: Empty string - distance equals length of non-empty string (Wikipedia definition)")]
        public void EditDistance_EmptyAndNonEmpty_ReturnsLength()
        {
            Assert.Multiple(() =>
            {
                Assert.That(ApproximateMatcher.EditDistance("", "ACGT"), Is.EqualTo(4));
                Assert.That(ApproximateMatcher.EditDistance("ACGT", ""), Is.EqualTo(4));
                Assert.That(ApproximateMatcher.EditDistance("", "abc"), Is.EqualTo(3));
                Assert.That(ApproximateMatcher.EditDistance("abc", ""), Is.EqualTo(3));
            });
        }

        [Test]
        [Description("M03: Canonical example - 'kitten' to 'sitting' = 3 (Wikipedia, Rosetta Code)")]
        public void EditDistance_KittenSitting_ReturnsThree()
        {
            // Wikipedia: k→s (sub), e→i (sub), +g (insert) = 3 edits
            int distance = ApproximateMatcher.EditDistance("kitten", "sitting");
            Assert.That(distance, Is.EqualTo(3));
        }

        [Test]
        [Description("M04: Canonical example - 'rosettacode' to 'raisethysword' = 8 (Rosetta Code)")]
        public void EditDistance_RosettacodeRaisethysword_ReturnsEight()
        {
            int distance = ApproximateMatcher.EditDistance("rosettacode", "raisethysword");
            Assert.That(distance, Is.EqualTo(8));
        }

        [Test]
        [Description("M05: Symmetry property - EditDistance(a,b) = EditDistance(b,a) (Metric property from Wikipedia)")]
        public void EditDistance_Symmetry_CommutativeProperty()
        {
            Assert.Multiple(() =>
            {
                Assert.That(
                    ApproximateMatcher.EditDistance("kitten", "sitting"),
                    Is.EqualTo(ApproximateMatcher.EditDistance("sitting", "kitten")));
                Assert.That(
                    ApproximateMatcher.EditDistance("rosettacode", "raisethysword"),
                    Is.EqualTo(ApproximateMatcher.EditDistance("raisethysword", "rosettacode")));
                Assert.That(
                    ApproximateMatcher.EditDistance("flaw", "lawn"),
                    Is.EqualTo(ApproximateMatcher.EditDistance("lawn", "flaw")));
            });
        }

        [Test]
        [Description("M06: Single substitution operation returns 1")]
        public void EditDistance_SingleSubstitution_ReturnsOne()
        {
            Assert.Multiple(() =>
            {
                Assert.That(ApproximateMatcher.EditDistance("ACGT", "ACGG"), Is.EqualTo(1)); // T→G
                Assert.That(ApproximateMatcher.EditDistance("cat", "bat"), Is.EqualTo(1)); // c→b
            });
        }

        [Test]
        [Description("M07: Single insertion operation returns 1")]
        public void EditDistance_SingleInsertion_ReturnsOne()
        {
            Assert.Multiple(() =>
            {
                Assert.That(ApproximateMatcher.EditDistance("ACGT", "ACGGT"), Is.EqualTo(1)); // insert G
                Assert.That(ApproximateMatcher.EditDistance("cat", "cats"), Is.EqualTo(1)); // insert s
            });
        }

        [Test]
        [Description("M08: Single deletion operation returns 1")]
        public void EditDistance_SingleDeletion_ReturnsOne()
        {
            Assert.Multiple(() =>
            {
                Assert.That(ApproximateMatcher.EditDistance("ACGT", "ACT"), Is.EqualTo(1)); // delete G
                Assert.That(ApproximateMatcher.EditDistance("cats", "cat"), Is.EqualTo(1)); // delete s
            });
        }

        [Test]
        [Description("M09: Null input throws ArgumentNullException")]
        public void EditDistance_NullInput_ThrowsArgumentNullException()
        {
            Assert.Multiple(() =>
            {
                Assert.Throws<ArgumentNullException>(() => ApproximateMatcher.EditDistance(null!, "test"));
                Assert.Throws<ArgumentNullException>(() => ApproximateMatcher.EditDistance("test", null!));
            });
        }

        [Test]
        [Description("M10: Case-sensitive comparison per standard Levenshtein definition (Wikipedia: head(a) = head(b))")]
        public void EditDistance_CaseSensitive_DistinguishesCase()
        {
            // Standard Levenshtein compares characters as distinct symbols.
            // Wikipedia pseudocode: 'if s[i] = t[j]' — no case normalization.
            Assert.Multiple(() =>
            {
                Assert.That(ApproximateMatcher.EditDistance("A", "a"), Is.EqualTo(1),
                    "'A' and 'a' are distinct characters: 1 substitution");
                Assert.That(ApproximateMatcher.EditDistance("Kitten", "kitten"), Is.EqualTo(1),
                    "Only 'K' vs 'k' differs: 1 substitution");
                Assert.That(ApproximateMatcher.EditDistance("ABC", "abc"), Is.EqualTo(3),
                    "All 3 characters differ in case: 3 substitutions");
            });
        }

        [Test]
        [Description("M13: 'flaw' to 'lawn' = 2, demonstrating Levenshtein < Hamming (Wikipedia example)")]
        public void EditDistance_FlawLawn_ReturnsTwo()
        {
            // Wikipedia: This demonstrates Levenshtein (2) can be less than Hamming (4)
            // Edits: delete 'f', insert 'n' = 2 edits
            int distance = ApproximateMatcher.EditDistance("flaw", "lawn");
            Assert.That(distance, Is.EqualTo(2));
        }

        #endregion

        #region EditDistance - SHOULD Tests (Good Practice)

        [Test]
        [Description("S01: 'saturday' to 'sunday' = 3 (Rosetta Code canonical example)")]
        public void EditDistance_SaturdaySunday_ReturnsThree()
        {
            int distance = ApproximateMatcher.EditDistance("saturday", "sunday");
            Assert.That(distance, Is.EqualTo(3));
        }

        [Test]
        [Description("S02: 'stop' to 'tops' = 2 (Rosetta Code)")]
        public void EditDistance_StopTops_ReturnsTwo()
        {
            int distance = ApproximateMatcher.EditDistance("stop", "tops");
            Assert.That(distance, Is.EqualTo(2));
        }

        [Test]
        [Description("S03: Both empty strings returns 0")]
        public void EditDistance_BothEmpty_ReturnsZero()
        {
            int distance = ApproximateMatcher.EditDistance("", "");
            Assert.That(distance, Is.EqualTo(0));
        }

        [Test]
        [Description("S07: Triangle inequality holds - d(a,c) ≤ d(a,b) + d(b,c) (Metric property)")]
        public void EditDistance_TriangleInequality_Holds()
        {
            // Wikipedia decomposition: kitten → sitten (1 sub) → sitting (2 edits)
            int ab = ApproximateMatcher.EditDistance("kitten", "sitten");
            int bc = ApproximateMatcher.EditDistance("sitten", "sitting");
            int ac = ApproximateMatcher.EditDistance("kitten", "sitting");

            Assert.Multiple(() =>
            {
                Assert.That(ab, Is.EqualTo(1), "kitten→sitten: 1 substitution (k→s)");
                Assert.That(bc, Is.EqualTo(2), "sitten→sitting: 2 edits (e→i, +g)");
                Assert.That(ac, Is.EqualTo(3), "kitten→sitting: 3 edits (canonical)");
                Assert.That(ac, Is.LessThanOrEqualTo(ab + bc),
                    "Triangle inequality: d(a,c) ≤ d(a,b) + d(b,c)");
            });
        }

        [Test]
        [Description("S08: Distance bounds - |len(a) - len(b)| ≤ distance ≤ max(len(a), len(b)) (Wikipedia)")]
        public void EditDistance_Bounds_WithinExpectedRange()
        {
            Assert.Multiple(() =>
            {
                // Same-length: "flaw"(4) vs "lawn"(4) → distance=2, bounds: 0 ≤ 2 ≤ 4
                Assert.That(ApproximateMatcher.EditDistance("flaw", "lawn"), Is.EqualTo(2));
                Assert.That(2, Is.GreaterThanOrEqualTo(Math.Abs(4 - 4)));
                Assert.That(2, Is.LessThanOrEqualTo(Math.Max(4, 4)));

                // Different-length: "saturday"(8) vs "sunday"(6) → distance=3, bounds: 2 ≤ 3 ≤ 8
                Assert.That(ApproximateMatcher.EditDistance("saturday", "sunday"), Is.EqualTo(3));
                Assert.That(3, Is.GreaterThanOrEqualTo(Math.Abs(8 - 6)));
                Assert.That(3, Is.LessThanOrEqualTo(Math.Max(8, 6)));
            });
        }

        #endregion

        #region FindWithEdits - MUST Tests

        [Test]
        [Description("M11: FindWithEdits with maxEdits=0 behaves like exact matching")]
        public void FindWithEdits_ExactMatch_Found()
        {
            var matches = ApproximateMatcher.FindWithEdits("ACGTACGT", "ACGT", 0).ToList();

            Assert.Multiple(() =>
            {
                Assert.That(matches, Has.Count.EqualTo(2));
                Assert.That(matches.All(m => m.Distance == 0), Is.True);
                Assert.That(matches[0].Position, Is.EqualTo(0));
                Assert.That(matches[1].Position, Is.EqualTo(4));
            });
        }

        [Test]
        [Description("M12: FindWithEdits with negative maxEdits throws ArgumentOutOfRangeException")]
        public void FindWithEdits_NegativeMaxEdits_ThrowsException()
        {
            Assert.Throws<ArgumentOutOfRangeException>(() =>
                ApproximateMatcher.FindWithEdits("ACGT", "AC", -1).ToList());
        }

        #endregion

        #region FindWithEdits - SHOULD Tests

        [Test]
        [Description("S04: FindWithEdits finds matches with substitution")]
        public void FindWithEdits_WithSubstitution_Found()
        {
            // "ACT" in "ACGT" with maxEdits=1 produces 3 windows:
            //   pos=0 "AC"   dist=1 (Edit)   — deletion of T from pattern
            //   pos=0 "ACG"  dist=1 (Sub)    — substitution G→T
            //   pos=0 "ACGT" dist=1 (Edit)   — insertion of T in pattern
            var matches = ApproximateMatcher.FindWithEdits("ACGT", "ACT", 1).ToList();

            Assert.Multiple(() =>
            {
                Assert.That(matches, Has.Count.EqualTo(3));

                var subMatch = matches.Single(m => m.MismatchType == MismatchType.Substitution);
                Assert.That(subMatch.Position, Is.EqualTo(0));
                Assert.That(subMatch.MatchedSequence, Is.EqualTo("ACG"));
                Assert.That(subMatch.Distance, Is.EqualTo(1));
            });
        }

        [Test]
        [Description("S05: FindWithEdits finds matches with insertion in pattern")]
        public void FindWithEdits_WithInsertion_Found()
        {
            // Pattern "ACGGT"(5) in text "ACGT"(4) with maxEdits=1
            // Only window: pos=0 "ACGT" dist=1 (Edit) — extra G deleted from pattern
            var matches = ApproximateMatcher.FindWithEdits("ACGT", "ACGGT", 1).ToList();

            Assert.Multiple(() =>
            {
                Assert.That(matches, Has.Count.EqualTo(1));
                Assert.That(matches[0].Position, Is.EqualTo(0));
                Assert.That(matches[0].MatchedSequence, Is.EqualTo("ACGT"));
                Assert.That(matches[0].Distance, Is.EqualTo(1));
                Assert.That(matches[0].MismatchType, Is.EqualTo(MismatchType.Edit));
            });
        }

        [Test]
        [Description("S06: Empty inputs return empty results")]
        public void FindWithEdits_EmptyInputs_ReturnsEmpty()
        {
            Assert.Multiple(() =>
            {
                Assert.That(ApproximateMatcher.FindWithEdits("", "ACGT", 1).ToList(), Is.Empty);
                Assert.That(ApproximateMatcher.FindWithEdits("ACGT", "", 1).ToList(), Is.Empty);
            });
        }

        [Test]
        [Description("S09: FindWithEdits finds matches requiring deletion from pattern")]
        public void FindWithEdits_WithDeletion_Found()
        {
            // Pattern "ACG"(3) in text "AC"(2) with maxEdits=1
            // Only window: pos=0 "AC" dist=1 (Edit) — G deleted from pattern to match text
            var matches = ApproximateMatcher.FindWithEdits("AC", "ACG", 1).ToList();

            Assert.Multiple(() =>
            {
                Assert.That(matches, Has.Count.EqualTo(1));
                Assert.That(matches[0].Position, Is.EqualTo(0));
                Assert.That(matches[0].MatchedSequence, Is.EqualTo("AC"));
                Assert.That(matches[0].Distance, Is.EqualTo(1));
                Assert.That(matches[0].MismatchType, Is.EqualTo(MismatchType.Edit));
            });
        }

        #endregion

        #region FindWithEdits - COULD Tests (Wrapper Verification)

        [Test]
        [Description("C01: DnaSequence overload delegates to string version")]
        public void FindWithEdits_DnaSequenceOverload_DelegatesToStringVersion()
        {
            var dnaSeq = new DnaSequence("ACGTACGT");
            var fromDna = ApproximateMatcher.FindWithEdits(dnaSeq, "ACGT", 0).ToList();
            var fromString = ApproximateMatcher.FindWithEdits("ACGTACGT", "ACGT", 0).ToList();

            Assert.Multiple(() =>
            {
                Assert.That(fromDna.Count, Is.EqualTo(fromString.Count));
                for (int i = 0; i < fromDna.Count; i++)
                {
                    Assert.That(fromDna[i].Position, Is.EqualTo(fromString[i].Position));
                    Assert.That(fromDna[i].Distance, Is.EqualTo(fromString[i].Distance));
                }
            });
        }

        [Test]
        [Description("S10: 'sleep' to 'fleeting' = 5 (Rosetta Code)")]
        public void EditDistance_SleepFleeting_ReturnsFive()
        {
            int distance = ApproximateMatcher.EditDistance("sleep", "fleeting");
            Assert.That(distance, Is.EqualTo(5));
        }

        #endregion
    }
}
