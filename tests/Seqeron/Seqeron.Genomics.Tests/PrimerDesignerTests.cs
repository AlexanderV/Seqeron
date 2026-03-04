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
            Assert.That(param.MinTm, Is.EqualTo(57));
            Assert.That(param.MaxTm, Is.EqualTo(63));
            Assert.That(param.MaxHomopolymer, Is.EqualTo(4));
        }

        #endregion

        #region Mutation-Killing Tests — EvaluatePrimer Tm Boundary Conditions

        /// <summary>
        /// Kills survived ||→&amp;&amp; mutation on EvaluatePrimer line 145:
        /// <c>if (tm &lt; param.MinTm || tm &gt; param.MaxTm)</c>
        /// Tests case where ONLY tm &lt; MinTm is true (Tm is NOT above MaxTm).
        /// With ||: issue flagged. With &amp;&amp;: no issue (both must be true).
        /// </summary>
        [Test]
        public void EvaluatePrimer_TmOnlyBelowMin_FlagsTmIssue()
        {
            // 20bp, 50% GC → Tm ≈ 51.8°C (Marmur-Doty formula)
            // Primer3 MinTm=57: tm(51.8) < 57 is TRUE, tm(51.8) > 63 is FALSE
            string primer = "ATGCGATCGATCGATCGATC"; // 20bp, 10 GC
            var param = new PrimerParameters(
                MinLength: 18, MaxLength: 25, OptimalLength: 20,
                MinGcContent: 0, MaxGcContent: 100,
                MinTm: 57, MaxTm: 63,
                OptimalTm: 60, MaxHomopolymer: 100, MaxDinucleotideRepeats: 100,
                Avoid3PrimeGC: false, Check3PrimeStability: false);

            var candidate = PrimerDesigner.EvaluatePrimer(primer, 0, true, param);

            Assert.That(candidate.Issues, Has.Some.Contains("Tm"),
                "Tm below MinTm should flag issue even when Tm is not above MaxTm");
        }

        /// <summary>
        /// Kills survived ||→&amp;&amp; mutation (reverse clause): tests case where
        /// ONLY tm &gt; MaxTm is true (Tm is NOT below MinTm).
        /// </summary>
        [Test]
        public void EvaluatePrimer_TmOnlyAboveMax_FlagsTmIssue()
        {
            // 20bp, 90% GC → Tm ≈ 68.2°C (Marmur-Doty)
            // Primer3 MaxTm=63: tm(68.2) > 63 is TRUE, tm(68.2) < 57 is FALSE
            string primer = "GCGCGCGCGCGCGCGCGCAT"; // 20bp, 18 GC
            var param = new PrimerParameters(
                MinLength: 18, MaxLength: 25, OptimalLength: 20,
                MinGcContent: 0, MaxGcContent: 100,
                MinTm: 57, MaxTm: 63,
                OptimalTm: 60, MaxHomopolymer: 100, MaxDinucleotideRepeats: 100,
                Avoid3PrimeGC: false, Check3PrimeStability: false);

            var candidate = PrimerDesigner.EvaluatePrimer(primer, 0, true, param);

            Assert.That(candidate.Issues, Has.Some.Contains("Tm"),
                "Tm above MaxTm should flag issue even when Tm is not below MinTm");
        }

        /// <summary>
        /// Kills survived &lt;→&lt;= boundary mutation on MinTm comparison.
        /// When Tm == MinTm, the original &lt; returns false (no issue).
        /// With &lt;= mutation, it returns true (issue added).
        /// </summary>
        [Test]
        public void EvaluatePrimer_TmExactlyAtMinTm_NoTmIssue()
        {
            string primer = "ATGCGATCGATCGATCGATC"; // 20bp
            double exactTm = PrimerDesigner.CalculateMeltingTemperature(primer);

            var param = new PrimerParameters(
                MinLength: 18, MaxLength: 25, OptimalLength: 20,
                MinGcContent: 0, MaxGcContent: 100,
                MinTm: exactTm, MaxTm: 100,
                OptimalTm: 60, MaxHomopolymer: 100, MaxDinucleotideRepeats: 100,
                Avoid3PrimeGC: false, Check3PrimeStability: false);

            var candidate = PrimerDesigner.EvaluatePrimer(primer, 0, true, param);

            Assert.That(candidate.Issues, Has.None.Contains("Tm"),
                "Tm exactly at MinTm boundary should not flag issue");
        }

        /// <summary>
        /// Kills survived &gt;→&gt;= boundary mutation on MaxTm comparison.
        /// When Tm == MaxTm, the original &gt; returns false (no issue).
        /// With &gt;= mutation, it returns true (issue added).
        /// </summary>
        [Test]
        public void EvaluatePrimer_TmExactlyAtMaxTm_NoTmIssue()
        {
            string primer = "ATGCGATCGATCGATCGATC"; // 20bp
            double exactTm = PrimerDesigner.CalculateMeltingTemperature(primer);

            var param = new PrimerParameters(
                MinLength: 18, MaxLength: 25, OptimalLength: 20,
                MinGcContent: 0, MaxGcContent: 100,
                MinTm: 0, MaxTm: exactTm,
                OptimalTm: 60, MaxHomopolymer: 100, MaxDinucleotideRepeats: 100,
                Avoid3PrimeGC: false, Check3PrimeStability: false);

            var candidate = PrimerDesigner.EvaluatePrimer(primer, 0, true, param);

            Assert.That(candidate.Issues, Has.None.Contains("Tm"),
                "Tm exactly at MaxTm boundary should not flag issue");
        }

        #endregion

        #region Mutation-Killing Tests — HasHairpinPotential Guard Clause

        /// <summary>
        /// Kills survived ||→&amp;&amp; mutation on HasHairpinPotential line 299:
        /// <c>if (string.IsNullOrEmpty(sequence) || sequence.Length &lt; minStemLength * 2 + minLoopLength)</c>
        /// With null input: || short-circuits to return false. With &amp;&amp;: evaluates
        /// null.Length → NullReferenceException.
        /// </summary>
        [Test]
        public void HasHairpinPotential_NullSequence_ReturnsFalse()
        {
            bool result = PrimerDesigner.HasHairpinPotential(null!);
            Assert.That(result, Is.False, "Null sequence should return false without exception");
        }

        [Test]
        public void HasHairpinPotential_EmptySequence_ReturnsFalse()
        {
            bool result = PrimerDesigner.HasHairpinPotential("");
            Assert.That(result, Is.False, "Empty sequence should return false without exception");
        }

        /// <summary>
        /// Kills survived &lt;→&lt;= boundary mutation on length check.
        /// With default params (stem=4, loop=3): threshold = 4*2+3 = 11.
        /// A sequence of exactly 11 chars should proceed past the guard (original &lt;).
        /// With &lt;= mutation, it would return false prematurely.
        /// </summary>
        [Test]
        public void HasHairpinPotential_SequenceExactlyAtThreshold_DoesNotReturnEarly()
        {
            // 11 chars = 4*2+3, exactly at threshold
            // With <: 11 < 11 = false → proceeds to hairpin check
            // With <=: 11 <= 11 = true → returns false prematurely
            // Use a self-complementary sequence that forms a hairpin: stem GCGC, loop AAA, stem GCGC
            string sequence = "GCGCAAAGCGC"; // 11 chars: stem(4) + loop(3) + stem(4)

            // The method should NOT return false from the guard clause
            // Whether it detects a hairpin depends on complementarity check
            // The important thing: with <= mutation, it returns false regardless of sequence content
            bool result = PrimerDesigner.HasHairpinPotential(sequence, minStemLength: 4, minLoopLength: 3);

            // GCGC complement is CGCG; reverse of GCGC is CGCG; complement of reverse is GCGC
            // AreComplementary checks if fragment matches complement of reversed target
            // The hairpin detection should at least run (not short-circuit)
            // We verify by testing a known hairpin-forming sequence at exact threshold
            Assert.That(result, Is.True,
                "Self-complementary sequence at exact threshold length should detect hairpin");
        }

        #endregion
    }
}
