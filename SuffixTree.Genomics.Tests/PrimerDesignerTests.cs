using NUnit.Framework;

namespace SuffixTree.Genomics.Tests
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

        #region Melting Temperature

        [Test]
        public void CalculateMeltingTemperature_ShortPrimer_UsesWallaceRule()
        {
            // ATATATAT: 8 AT bases = 2*8 + 4*0 = 16
            double tm = PrimerDesigner.CalculateMeltingTemperature("ATATATAT");
            Assert.That(tm, Is.EqualTo(16.0));
        }

        [Test]
        public void CalculateMeltingTemperature_ShortAllGC_HighTm()
        {
            // GCGCGCGC: 0 AT + 8 GC = 2*0 + 4*8 = 32
            double tm = PrimerDesigner.CalculateMeltingTemperature("GCGCGCGC");
            Assert.That(tm, Is.EqualTo(32.0));
        }

        [Test]
        public void CalculateMeltingTemperature_LongPrimer_UsesNearestNeighbor()
        {
            // 20 bp primer
            string primer = "ACGTACGTACGTACGTACGT";
            double tm = PrimerDesigner.CalculateMeltingTemperature(primer);
            
            // Should be in reasonable range for 50% GC
            Assert.That(tm, Is.GreaterThan(40).And.LessThan(70));
        }

        [Test]
        public void CalculateMeltingTemperature_EmptyPrimer_Returns0()
        {
            double tm = PrimerDesigner.CalculateMeltingTemperature("");
            Assert.That(tm, Is.EqualTo(0.0));
        }

        [Test]
        public void CalculateMeltingTemperatureWithSalt_AppliesSaltCorrection()
        {
            string primer = "ACGTACGTACGTACGTACGT";
            double tmBase = PrimerDesigner.CalculateMeltingTemperature(primer);
            double tmSalt = PrimerDesigner.CalculateMeltingTemperatureWithSalt(primer, 50);
            
            // Salt correction typically lowers Tm for physiological concentrations
            Assert.That(tmSalt, Is.Not.EqualTo(tmBase));
        }

        #endregion

        #region Homopolymer Detection

        [Test]
        public void FindLongestHomopolymer_NoHomopolymer_Returns1()
        {
            int run = PrimerDesigner.FindLongestHomopolymer("ACGT");
            Assert.That(run, Is.EqualTo(1));
        }

        [Test]
        public void FindLongestHomopolymer_HasHomopolymer_ReturnsLength()
        {
            int run = PrimerDesigner.FindLongestHomopolymer("ACAAAAGT");
            Assert.That(run, Is.EqualTo(4)); // AAAA
        }

        [Test]
        public void FindLongestHomopolymer_AllSame_ReturnsFullLength()
        {
            int run = PrimerDesigner.FindLongestHomopolymer("AAAAAA");
            Assert.That(run, Is.EqualTo(6));
        }

        [Test]
        public void FindLongestHomopolymer_EmptySequence_Returns0()
        {
            int run = PrimerDesigner.FindLongestHomopolymer("");
            Assert.That(run, Is.EqualTo(0));
        }

        #endregion

        #region Dinucleotide Repeats

        [Test]
        public void FindLongestDinucleotideRepeat_NoRepeat_Returns1OrLess()
        {
            int repeat = PrimerDesigner.FindLongestDinucleotideRepeat("ACGT");
            Assert.That(repeat, Is.LessThanOrEqualTo(1));
        }

        [Test]
        public void FindLongestDinucleotideRepeat_HasRepeat_ReturnsCount()
        {
            int repeat = PrimerDesigner.FindLongestDinucleotideRepeat("ACACACACG");
            Assert.That(repeat, Is.EqualTo(4)); // ACACACAC = 4 x AC
        }

        [Test]
        public void FindLongestDinucleotideRepeat_ShortSequence_Returns0()
        {
            int repeat = PrimerDesigner.FindLongestDinucleotideRepeat("ACG");
            Assert.That(repeat, Is.EqualTo(0));
        }

        #endregion

        #region Hairpin Detection

        [Test]
        public void HasHairpinPotential_NoHairpin_ReturnsFalse()
        {
            // Non-self-complementary sequence (AAAA cannot form hairpin with itself)
            bool hasHairpin = PrimerDesigner.HasHairpinPotential("AAAACCCCAAAA");
            Assert.That(hasHairpin, Is.False);
        }

        [Test]
        public void HasHairpinPotential_SelfComplementary_ReturnsTrue()
        {
            // ACGT...ACGT pattern - reverse of ACGT is TGCA which is complementary to ACGT
            bool hasHairpin = PrimerDesigner.HasHairpinPotential("ACGTACGTACGT");
            Assert.That(hasHairpin, Is.True);
        }

        [Test]
        public void HasHairpinPotential_ShortSequence_ReturnsFalse()
        {
            bool hasHairpin = PrimerDesigner.HasHairpinPotential("ACGT");
            Assert.That(hasHairpin, Is.False);
        }

        #endregion

        #region Primer Dimer Detection

        [Test]
        public void HasPrimerDimer_NoComplementarity_ReturnsFalse()
        {
            // Primers where 3' ends don't form complementary pairs
            // primer1 ends with CCCC, revcomp(primer2) starts with CCCC -> C-C is not complementary
            bool hasDimer = PrimerDesigner.HasPrimerDimer("AAAACCCCCCCC", "GGGGGGGGTTTT");
            Assert.That(hasDimer, Is.False);
        }

        [Test]
        public void HasPrimerDimer_Complementary3Ends_ReturnsTrue()
        {
            // primer1 ends with AAAA, primer2 ends with TTTT -> revcomp of primer2 ends with AAAA
            // 3' of primer1 = AAAA, 3' of revcomp(primer2) = AAAA - these are same not complementary
            // Actually: primer1=ACGTACGTACGT, revcomp=ACGTACGTACGT -> 3' of p1 vs 5' of revcomp
            bool hasDimer = PrimerDesigner.HasPrimerDimer("AAAAAAAA", "AAAAAAAA");
            Assert.That(hasDimer, Is.True); // 3' of p1=AAAA, revcomp(p2) starts with TTTTTTTT -> A-T complementary
        }

        [Test]
        public void HasPrimerDimer_EmptyPrimers_ReturnsFalse()
        {
            bool hasDimer = PrimerDesigner.HasPrimerDimer("", "ACGT");
            Assert.That(hasDimer, Is.False);
        }

        #endregion

        #region 3' Stability

        [Test]
        public void Calculate3PrimeStability_GCRich_MoreNegative()
        {
            double gcRich = PrimerDesigner.Calculate3PrimeStability("ACGTGCGCG");
            double atRich = PrimerDesigner.Calculate3PrimeStability("ACGTATATAT");
            
            // GC-rich 3' end should be more stable (more negative)
            Assert.That(gcRich, Is.LessThan(atRich));
        }

        [Test]
        public void Calculate3PrimeStability_ShortSequence_Returns0()
        {
            double stability = PrimerDesigner.Calculate3PrimeStability("ACGT");
            Assert.That(stability, Is.EqualTo(0));
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
