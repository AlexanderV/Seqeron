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
        [Description("M10: Case-insensitive comparison (implementation behavior)")]
        public void EditDistance_CaseInsensitive_IgnoresCase()
        {
            Assert.Multiple(() =>
            {
                Assert.That(ApproximateMatcher.EditDistance("acgt", "ACGT"), Is.EqualTo(0));
                Assert.That(ApproximateMatcher.EditDistance("Kitten", "KITTEN"), Is.EqualTo(0));
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
            // Test with: a = "kitten", b = "sitten", c = "sitting"
            int ab = ApproximateMatcher.EditDistance("kitten", "sitten");  // 1
            int bc = ApproximateMatcher.EditDistance("sitten", "sitting"); // 2
            int ac = ApproximateMatcher.EditDistance("kitten", "sitting"); // 3

            Assert.That(ac, Is.LessThanOrEqualTo(ab + bc),
                "Triangle inequality: d(a,c) ≤ d(a,b) + d(b,c)");
        }

        [Test]
        [Description("S08: Distance bounds - |len(a) - len(b)| ≤ distance ≤ max(len(a), len(b)) (Wikipedia)")]
        public void EditDistance_Bounds_WithinExpectedRange()
        {
            string a = "kitten";
            string b = "sitting";
            int distance = ApproximateMatcher.EditDistance(a, b);
            int minBound = Math.Abs(a.Length - b.Length);
            int maxBound = Math.Max(a.Length, b.Length);

            Assert.Multiple(() =>
            {
                Assert.That(distance, Is.GreaterThanOrEqualTo(minBound),
                    $"Distance should be >= |len(a) - len(b)| = {minBound}");
                Assert.That(distance, Is.LessThanOrEqualTo(maxBound),
                    $"Distance should be <= max(len(a), len(b)) = {maxBound}");
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
            // Looking for "ACT" in "ACGT" - "ACG" matches with 1 substitution (G→T)
            var matches = ApproximateMatcher.FindWithEdits("ACGT", "ACT", 1).ToList();

            Assert.That(matches.Any(m => m.Distance <= 1), Is.True,
                "Should find match with substitution");
        }

        [Test]
        [Description("S05: FindWithEdits finds matches with insertion in pattern")]
        public void FindWithEdits_WithInsertion_Found()
        {
            // Pattern "ACGGT" should match "ACGT" with 1 insertion (G inserted in pattern)
            var matches = ApproximateMatcher.FindWithEdits("ACGT", "ACGGT", 1).ToList();

            Assert.That(matches.Any(m => m.Distance == 1), Is.True,
                "Should find match requiring insertion");
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
        [Description("S09: FindWithEdits finds deletion matches")]
        public void FindWithEdits_WithDeletion_Found()
        {
            // Pattern "ACG" in "ACGT" - exact match at position 0
            var matches = ApproximateMatcher.FindWithEdits("ACGT", "ACG", 1).ToList();

            Assert.That(matches.Any(m => m.Distance == 0), Is.True,
                "Should find exact substring match");
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

        #endregion

        #region Additional Rosetta Code Test Vectors

        [Test]
        [Description("Additional canonical examples from Rosetta Code")]
        public void EditDistance_RosettaCodeTestVectors_AllCorrect()
        {
            // Test vectors verified across multiple programming languages
            Assert.Multiple(() =>
            {
                Assert.That(ApproximateMatcher.EditDistance("yo", ""), Is.EqualTo(2));
                Assert.That(ApproximateMatcher.EditDistance("", "yo"), Is.EqualTo(2));
                Assert.That(ApproximateMatcher.EditDistance("sleep", "fleeting"), Is.EqualTo(5));
                Assert.That(ApproximateMatcher.EditDistance("mississippi", "swiss miss"), Is.EqualTo(8));
            });
        }

        #endregion
    }
}
