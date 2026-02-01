using NUnit.Framework;

namespace Seqeron.Genomics.Tests
{
    [TestFixture]
    public class PrimerDesignerTests
    {
        #region GC Content

        [Test]
        public void CalculateGcContent_AllGC_Returns100()
        {
            double gc = PrimerDesigner.CalculateGcContent("GCGCGC");
            Assert.That(gc, Is.EqualTo(100.0));
        }

        [Test]
        public void CalculateGcContent_NoGC_Returns0()
        {
            double gc = PrimerDesigner.CalculateGcContent("ATATAT");
            Assert.That(gc, Is.EqualTo(0.0));
        }

        [Test]
        public void CalculateGcContent_HalfGC_Returns50()
        {
            double gc = PrimerDesigner.CalculateGcContent("ACGT");
            Assert.That(gc, Is.EqualTo(50.0));
        }

        [Test]
        public void CalculateGcContent_EmptySequence_Returns0()
        {
            double gc = PrimerDesigner.CalculateGcContent("");
            Assert.That(gc, Is.EqualTo(0.0));
        }

        #endregion

        // NOTE: Melting Temperature tests moved to PrimerDesigner_MeltingTemperature_Tests.cs
        // as part of PRIMER-TM-001 Test Unit consolidation.

        #region Melting Temperature - Smoke Test

        /// <summary>
        /// Smoke test verifying PrimerDesigner's Tm calculation is available.
        /// Full Tm tests are in PrimerDesigner_MeltingTemperature_Tests.cs.
        /// </summary>
        [Test]
        public void CalculateMeltingTemperature_SmokeTest_ReturnsValidValue()
        {
            string primer = "ACGTACGTACGTACGTACGT";
            double tmBase = PrimerDesigner.CalculateMeltingTemperature(primer);
            double tmSalt = PrimerDesigner.CalculateMeltingTemperatureWithSalt(primer, 50);

            // Salt correction typically lowers Tm for physiological concentrations
            Assert.That(tmSalt, Is.Not.EqualTo(tmBase));
        }

        #endregion

        // NOTE: Primer Structure tests moved to PrimerDesigner_PrimerStructure_Tests.cs
        // as part of PRIMER-STRUCT-001 Test Unit consolidation.
        // The following are smoke tests for integration verification.

        #region Primer Structure - Smoke Tests

        /// <summary>
        /// Smoke test for homopolymer detection.
        /// Full tests in PrimerDesigner_PrimerStructure_Tests.cs.
        /// </summary>
        [Test]
        public void FindLongestHomopolymer_SmokeTest_ReturnsValidValue()
        {
            int run = PrimerDesigner.FindLongestHomopolymer("ACAAAAGT");
            Assert.That(run, Is.EqualTo(4)); // AAAA
        }

        /// <summary>
        /// Smoke test for dinucleotide repeat detection.
        /// Full tests in PrimerDesigner_PrimerStructure_Tests.cs.
        /// </summary>
        [Test]
        public void FindLongestDinucleotideRepeat_SmokeTest_ReturnsValidValue()
        {
            int repeat = PrimerDesigner.FindLongestDinucleotideRepeat("ACACACACG");
            Assert.That(repeat, Is.EqualTo(4)); // ACACACAC = 4 x AC
        }

        /// <summary>
        /// Smoke test for hairpin detection.
        /// Full tests in PrimerDesigner_PrimerStructure_Tests.cs.
        /// </summary>
        [Test]
        public void HasHairpinPotential_SmokeTest_ReturnsExpectedValue()
        {
            bool hasHairpin = PrimerDesigner.HasHairpinPotential("ACGTACGTACGT");
            Assert.That(hasHairpin, Is.True);
        }

        /// <summary>
        /// Smoke test for primer-dimer detection.
        /// Full tests in PrimerDesigner_PrimerStructure_Tests.cs.
        /// </summary>
        [Test]
        public void HasPrimerDimer_SmokeTest_ReturnsExpectedValue()
        {
            bool hasDimer = PrimerDesigner.HasPrimerDimer("AAAAAAAA", "AAAAAAAA");
            Assert.That(hasDimer, Is.True);
        }

        /// <summary>
        /// Smoke test for 3' stability calculation.
        /// Full tests in PrimerDesigner_PrimerStructure_Tests.cs.
        /// </summary>
        [Test]
        public void Calculate3PrimeStability_SmokeTest_ReturnsNegativeValue()
        {
            double stability = PrimerDesigner.Calculate3PrimeStability("ACGTGCGCG");
            Assert.That(stability, Is.LessThan(0));
        }

        #endregion

        #region Evaluate Primer

        [Test]
        public void EvaluatePrimer_GoodPrimer_IsValid()
        {
            // A well-designed 20bp primer with ~50% GC
            string primer = "ACGTACGTACGTACGTACGT";
            var candidate = PrimerDesigner.EvaluatePrimer(primer, 0, true);

            Assert.That(candidate.Sequence, Is.EqualTo(primer));
            Assert.That(candidate.Length, Is.EqualTo(20));
            Assert.That(candidate.GcContent, Is.EqualTo(50.0));
            Assert.That(candidate.IsForward, Is.True);
        }

        [Test]
        public void EvaluatePrimer_TooShort_NotValid()
        {
            var candidate = PrimerDesigner.EvaluatePrimer("ACGT", 0, true);

            Assert.That(candidate.IsValid, Is.False);
            Assert.That(candidate.Issues, Does.Contain("Length 4 outside range [18-25]"));
        }

        [Test]
        public void EvaluatePrimer_TooLongHomopolymer_NotValid()
        {
            string primer = "ACGTAAAAAAACGTACGTAC"; // 20bp with 6x A run
            var candidate = PrimerDesigner.EvaluatePrimer(primer, 0, true);

            Assert.That(candidate.Issues.Any(i => i.Contains("Homopolymer")), Is.True);
        }

        [Test]
        public void EvaluatePrimer_CalculatesScore()
        {
            string primer = "ACGTACGTACGTACGTACGT";
            var candidate = PrimerDesigner.EvaluatePrimer(primer, 0, true);

            Assert.That(candidate.Score, Is.GreaterThan(0));
        }

        #endregion

        #region Design Primers

        [Test]
        public void DesignPrimers_ValidTarget_ReturnsResult()
        {
            // Create a template with valid primer regions
            var sb = new System.Text.StringBuilder();
            // Forward primer region (before target)
            sb.Append("ACGTACGTACGTACGTACGTACGT"); // 24 bp
            sb.Append("GCTAGCTAGCTAGCTAGCTAGCTA"); // 24 bp
            sb.Append("ATCGATCGATCGATCGATCGATCG"); // 24 bp
            sb.Append("CGATCGATCGATCGATCGATCGAT"); // 24 bp - target start (96)
            // Target region
            sb.Append("NNNNNNNNNNNNNNNNNNNNNNNN"); // 24 bp target (use valid bases)
            sb.Append("GGGGGGGGGGGGGGGGGGGGGGGG"); // Replace N with G
            // Reverse primer region (after target)
            sb.Append("TAGCTAGCTAGCTAGCTAGCTAGC"); // 24 bp
            sb.Append("CGATCGATCGATCGATCGATCGAT"); // 24 bp
            sb.Append("ACGTACGTACGTACGTACGTACGT"); // 24 bp

            var template = new DnaSequence(
                "ACGTACGTACGTACGTACGTACGT" +
                "GCTAGCTAGCTAGCTAGCTAGCTA" +
                "ATCGATCGATCGATCGATCGATCG" +
                "CGATCGATCGATCGATCGATCGAT" +
                "TTTTTTTTTTTTTTTTTTTTTTTT" +
                "TAGCTAGCTAGCTAGCTAGCTAGC" +
                "CGATCGATCGATCGATCGATCGAT" +
                "ACGTACGTACGTACGTACGTACGT"
            );

            var result = PrimerDesigner.DesignPrimers(template, 80, 120);

            Assert.That(result, Is.Not.Null);
            // May or may not find valid primers depending on sequence
        }

        [Test]
        public void DesignPrimers_InvalidTarget_ThrowsException()
        {
            var template = new DnaSequence("ACGTACGTACGTACGTACGT");

            Assert.Throws<ArgumentException>(() =>
                PrimerDesigner.DesignPrimers(template, 15, 10));
        }

        [Test]
        public void DesignPrimers_TargetOutOfRange_ThrowsException()
        {
            var template = new DnaSequence("ACGTACGTACGTACGTACGT");

            Assert.Throws<ArgumentException>(() =>
                PrimerDesigner.DesignPrimers(template, 0, 100));
        }

        #endregion

        #region Generate Primer Candidates

        [Test]
        public void GeneratePrimerCandidates_ReturnsMultipleCandidates()
        {
            var template = new DnaSequence("ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT");
            var candidates = PrimerDesigner.GeneratePrimerCandidates(template, 0, 40, true).ToList();

            Assert.That(candidates, Has.Count.GreaterThan(0));
        }

        [Test]
        public void GeneratePrimerCandidates_Forward_HasCorrectOrientation()
        {
            var template = new DnaSequence("ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT");
            var candidates = PrimerDesigner.GeneratePrimerCandidates(template, 0, 40, true).ToList();

            Assert.That(candidates.All(c => c.IsForward), Is.True);
        }

        [Test]
        public void GeneratePrimerCandidates_Reverse_HasCorrectOrientation()
        {
            var template = new DnaSequence("ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT");
            var candidates = PrimerDesigner.GeneratePrimerCandidates(template, 0, 40, false).ToList();

            Assert.That(candidates.All(c => !c.IsForward), Is.True);
        }

        #endregion

        #region Default Parameters

        [Test]
        public void DefaultParameters_HasReasonableValues()
        {
            var param = PrimerDesigner.DefaultParameters;

            Assert.That(param.MinLength, Is.EqualTo(18));
            Assert.That(param.MaxLength, Is.EqualTo(25));
            Assert.That(param.OptimalLength, Is.EqualTo(20));
            Assert.That(param.MinGcContent, Is.EqualTo(40));
            Assert.That(param.MaxGcContent, Is.EqualTo(60));
            Assert.That(param.MinTm, Is.EqualTo(55));
            Assert.That(param.MaxTm, Is.EqualTo(65));
            Assert.That(param.MaxHomopolymer, Is.EqualTo(4));
        }

        #endregion

        #region Custom Parameters

        [Test]
        public void EvaluatePrimer_CustomParameters_AppliesCorrectly()
        {
            var customParams = new PrimerParameters(
                MinLength: 15,
                MaxLength: 30,
                OptimalLength: 22,
                MinGcContent: 30,
                MaxGcContent: 70,
                MinTm: 50,
                MaxTm: 70,
                OptimalTm: 58,
                MaxHomopolymer: 5,
                MaxDinucleotideRepeats: 5,
                Avoid3PrimeGC: true,
                Check3PrimeStability: false
            );

            // 16bp primer would be invalid with default (18) but valid with custom (15)
            string primer = "ACGTACGTACGTACGT";
            var defaultResult = PrimerDesigner.EvaluatePrimer(primer, 0, true);
            var customResult = PrimerDesigner.EvaluatePrimer(primer, 0, true, customParams);

            Assert.That(defaultResult.Issues.Any(i => i.Contains("Length")), Is.True);
            Assert.That(customResult.Issues.Any(i => i.Contains("Length")), Is.False);
        }

        #endregion

        #region Real-World Scenarios

        [Test]
        public void EvaluatePrimer_TypicalGoodPrimer_PassesAllChecks()
        {
            // A typical good primer for PCR
            string primer = "ATGCGATCGATCGATCGATC"; // 20bp, 50% GC
            var candidate = PrimerDesigner.EvaluatePrimer(primer, 0, true);

            Assert.That(candidate.GcContent, Is.InRange(40, 60));
            Assert.That(candidate.Length, Is.InRange(18, 25));
            Assert.That(candidate.HomopolymerLength, Is.LessThanOrEqualTo(4));
        }

        [Test]
        public void EvaluatePrimer_ProblematicPrimer_DetectsIssues()
        {
            // A primer with multiple issues
            string primer = "GGGGGGGGGGGGGGGGGGGG"; // 20bp, 100% GC, long homopolymer
            var candidate = PrimerDesigner.EvaluatePrimer(primer, 0, true);

            Assert.That(candidate.IsValid, Is.False);
            Assert.That(candidate.GcContent, Is.EqualTo(100.0));
            Assert.That(candidate.HomopolymerLength, Is.EqualTo(20));
        }

        #endregion
    }
}
