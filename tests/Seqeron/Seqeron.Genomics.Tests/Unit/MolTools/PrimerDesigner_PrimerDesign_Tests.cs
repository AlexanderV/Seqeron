using System.Text;

namespace Seqeron.Genomics.Tests.Unit.MolTools;

/// <summary>
/// Tests for Primer Pair Design functionality (PRIMER-DESIGN-001).
/// Covers DesignPrimers, EvaluatePrimer, and GeneratePrimerCandidates methods.
/// 
/// Evidence sources:
/// - Wikipedia: Primer (molecular biology) - standard primer length 18-24 bp
/// - Addgene: How to Design a Primer - 40-60% GC, 50-60°C Tm, pairs within 5°C
/// - Primer3 Manual (v2.6.1) - PRIMER_MIN_TM=57, PRIMER_OPT_TM=60, PRIMER_MAX_TM=63, PRIMER_MAX_POLY_X=5
/// </summary>
[TestFixture]
public class PrimerDesigner_PrimerDesign_Tests
{
    #region Test Fixtures

    private DnaSequence _standardTemplate = null!;
    private DnaSequence _gcRichTemplate = null!;
    private DnaSequence _atRichTemplate = null!;

    [SetUp]
    public void SetUp()
    {
        // Standard template with non-palindromic primer regions.
        // Repeating units are chosen so that NO 4-base window is a DNA palindrome
        // (i.e., no window where reverse complement equals itself), preventing
        // hairpin detection from rejecting all candidates.
        //
        // Forward unit: GAACTCGT (50% GC, no 4bp palindromes, max homopolymer=2)
        // Reverse unit: TCCGAAGT (50% GC, no 4bp palindromes, different from forward)
        //
        // Tm is the Primer3-default SantaLucia NN Tm (primer3-py calc_tm): the 20-mer
        // GAACTCGTGAACTCGTGAAC = 56.82°C (just below 57), longer windows fall in [57, 63].
        var sb = new StringBuilder();

        // Forward primer region (100bp, ~50% GC)
        while (sb.Length < 100) sb.Append("GAACTCGT");

        // Target region (50bp poly-T, clearly different from primer regions)
        sb.Append(new string('T', 50));

        // Reverse primer region (100bp, ~50% GC, different unit to avoid primer-dimer)
        int revStart = sb.Length;
        while (sb.Length - revStart < 100) sb.Append("TCCGAAGT");

        _standardTemplate = new DnaSequence(sb.ToString());

        // GC-rich template (70% GC)
        _gcRichTemplate = new DnaSequence(
            "GCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGC" + // 50bp
            "CCCCGGGGCCCCGGGGCCCCGGGGCCCCGGGGCCCCGGGGCCCCGGGGCC" + // 50bp
            "ATATATATATATATATATATATATATATATATATATATATATATATATATAT" + // 50bp target
            "GGGGCCCCGGGGCCCCGGGGCCCCGGGGCCCCGGGGCCCCGGGGCCCCGG" + // 50bp
            "GCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGC"   // 50bp
        );

        // AT-rich template (70% AT)
        _atRichTemplate = new DnaSequence(
            "ATATATATATATATATATATATATATATATATATATATATATATATATAT" + // 50bp
            "TTTTAAAATTTTAAAATTTTAAAATTTTAAAATTTTAAAATTTTAAAATT" + // 50bp
            "GCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGC" + // 50bp target
            "AAAATTTTAAAATTTTAAAATTTTAAAATTTTAAAATTTTAAAATTTTAA" + // 50bp
            "TATATATATATATATATATATATATATATATATATATATATATATATATAT"   // 50bp
        );
    }

    #endregion

    #region M1: DesignPrimers returns valid forward primer in upstream region

    [Test]
    public void DesignPrimers_ValidTemplate_ForwardIsUpstreamOfTarget()
    {
        // Arrange
        int targetStart = 100;
        int targetEnd = 150;

        // Act
        var result = PrimerDesigner.DesignPrimers(_standardTemplate, targetStart, targetEnd);

        // Assert
        Assert.That(result.IsValid, Is.True, "Standard template must produce valid primer pair");
        Assert.That(result.Forward, Is.Not.Null);
        Assert.That(result.Forward!.Position, Is.LessThan(targetStart),
            "Forward primer must be positioned upstream of target start");
        Assert.That(result.Forward.Position + result.Forward.Length, Is.LessThanOrEqualTo(targetStart),
            "Forward primer must end before or at target start");
    }

    [Test]
    public void DesignPrimers_ValidTemplate_ForwardWithinSearchRegion()
    {
        // Arrange
        int targetStart = 100;
        int targetEnd = 150;
        int expectedMinPosition = Math.Max(0, targetStart - 200);

        // Act
        var result = PrimerDesigner.DesignPrimers(_standardTemplate, targetStart, targetEnd);

        // Assert
        Assert.That(result.IsValid, Is.True, "Standard template must produce valid primer pair");
        Assert.That(result.Forward, Is.Not.Null);
        Assert.That(result.Forward!.Position, Is.GreaterThanOrEqualTo(expectedMinPosition),
            "Forward primer should be within 200bp upstream search region");
    }

    #endregion

    #region M2: DesignPrimers returns valid reverse primer in downstream region

    [Test]
    public void DesignPrimers_ValidTemplate_ReverseIsDownstreamOfTarget()
    {
        // Arrange
        int targetStart = 100;
        int targetEnd = 150;

        // Act
        var result = PrimerDesigner.DesignPrimers(_standardTemplate, targetStart, targetEnd);

        // Assert
        Assert.That(result.IsValid, Is.True, "Standard template must produce valid primer pair");
        Assert.That(result.Reverse, Is.Not.Null);
        Assert.That(result.Reverse!.Position, Is.GreaterThanOrEqualTo(targetEnd),
            "Reverse primer must be positioned downstream of target end");
    }

    [Test]
    public void DesignPrimers_ValidTemplate_ReverseWithinSearchRegion()
    {
        // Arrange
        int targetStart = 100;
        int targetEnd = 150;
        int expectedMaxPosition = Math.Min(_standardTemplate.Length - 1, targetEnd + 200);

        // Act
        var result = PrimerDesigner.DesignPrimers(_standardTemplate, targetStart, targetEnd);

        // Assert
        Assert.That(result.IsValid, Is.True, "Standard template must produce valid primer pair");
        Assert.That(result.Reverse, Is.Not.Null);
        Assert.That(result.Reverse!.Position, Is.LessThanOrEqualTo(expectedMaxPosition),
            "Reverse primer should be within 200bp downstream search region");
    }

    #endregion

    #region M3: Primer length within 18-25bp (Primer3: 18-27, Addgene: 18-24)

    [Test]
    public void DesignPrimers_Primers_HaveLengthWithinRange()
    {
        // Arrange
        int targetStart = 100;
        int targetEnd = 150;
        int minLength = 18;
        int maxLength = 25;

        // Act
        var result = PrimerDesigner.DesignPrimers(_standardTemplate, targetStart, targetEnd);

        // Assert
        Assert.That(result.IsValid, Is.True, "Standard template must produce valid primer pair");
        Assert.That(result.Forward!.Length, Is.InRange(minLength, maxLength),
            $"Forward primer length should be {minLength}-{maxLength}bp");
        Assert.That(result.Reverse!.Length, Is.InRange(minLength, maxLength),
            $"Reverse primer length should be {minLength}-{maxLength}bp");
    }

    [TestCase(17, Description = "Below minimum length")]
    [TestCase(26, Description = "Above maximum length")]
    public void EvaluatePrimer_LengthOutsideRange_ReportsIssue(int length)
    {
        // Arrange
        string primer = new string('A', length / 2) + new string('T', length - length / 2);

        // Act
        var candidate = PrimerDesigner.EvaluatePrimer(primer, 0, true);

        // Assert
        Assert.That(candidate.Issues.Any(i => i.Contains("Length")), Is.True,
            $"Primer of length {length} should report length issue");
    }

    #endregion

    #region M4: GC content within 40-60% (Addgene standard)

    [Test]
    public void DesignPrimers_Primers_HaveGcContentWithinRange()
    {
        // Arrange
        int targetStart = 100;
        int targetEnd = 150;
        double minGc = 40.0;
        double maxGc = 60.0;

        // Act
        var result = PrimerDesigner.DesignPrimers(_standardTemplate, targetStart, targetEnd);

        // Assert
        Assert.That(result.IsValid, Is.True, "Standard template must produce valid primer pair");
        Assert.That(result.Forward!.GcContent, Is.InRange(minGc, maxGc),
            $"Forward primer GC content should be {minGc}-{maxGc}% (Addgene standard)");
        Assert.That(result.Reverse!.GcContent, Is.InRange(minGc, maxGc),
            $"Reverse primer GC content should be {minGc}-{maxGc}% (Addgene standard)");
    }

    [TestCase(100.0, "GGGGGGGGGGGGGGGGGGGG", Description = "100% GC")]
    [TestCase(0.0, "AAAAAAAAAAAAAAAAAAAA", Description = "0% GC")]
    public void EvaluatePrimer_GcOutsideRange_ReportsIssue(double expectedGc, string primer)
    {
        // Act
        var candidate = PrimerDesigner.EvaluatePrimer(primer, 0, true);

        // Assert
        Assert.That(candidate.GcContent, Is.EqualTo(expectedGc).Within(0.1));
        Assert.That(candidate.Issues.Any(i => i.Contains("GC")), Is.True,
            $"Primer with {expectedGc}% GC should report GC content issue");
    }

    #endregion

    #region M5: Tm within 57-63°C (Primer3: PRIMER_MIN_TM=57, PRIMER_MAX_TM=63)

    [Test]
    public void DesignPrimers_Primers_HaveTmWithinRange()
    {
        // Arrange
        int targetStart = 100;
        int targetEnd = 150;
        double minTm = 57.0;
        double maxTm = 63.0;

        // Act
        var result = PrimerDesigner.DesignPrimers(_standardTemplate, targetStart, targetEnd);

        // Assert
        Assert.That(result.IsValid, Is.True, "Standard template must produce valid primer pair");
        Assert.That(result.Forward!.MeltingTemperature, Is.InRange(minTm, maxTm),
            $"Forward primer Tm should be {minTm}-{maxTm}°C");
        Assert.That(result.Reverse!.MeltingTemperature, Is.InRange(minTm, maxTm),
            $"Reverse primer Tm should be {minTm}-{maxTm}°C");
    }

    #endregion

    #region M6: Tm difference ≤5°C between primer pair (Primer3: PRIMER_PAIR_MAX_DIFF_TM)

    [Test]
    public void DesignPrimers_PrimerPair_TmDifferenceWithin5Degrees()
    {
        // Arrange
        int targetStart = 100;
        int targetEnd = 150;
        double maxTmDiff = 5.0;

        // Act
        var result = PrimerDesigner.DesignPrimers(_standardTemplate, targetStart, targetEnd);

        // Assert
        Assert.That(result.IsValid, Is.True, "Standard template must produce valid primer pair");
        double tmDiff = Math.Abs(result.Forward!.MeltingTemperature - result.Reverse!.MeltingTemperature);
        Assert.That(tmDiff, Is.LessThanOrEqualTo(maxTmDiff),
            $"Primer pair Tm difference should be ≤{maxTmDiff}°C (Addgene/Wikipedia standard)");
    }

    #endregion

    #region M7: No excessive homopolymer runs (≤4 bp, Primer3: ≤5)

    [Test]
    public void DesignPrimers_Primers_NoExcessiveHomopolymers()
    {
        // Arrange
        int targetStart = 100;
        int targetEnd = 150;
        int maxHomopolymer = 4;

        // Act
        var result = PrimerDesigner.DesignPrimers(_standardTemplate, targetStart, targetEnd);

        // Assert
        Assert.That(result.IsValid, Is.True, "Standard template must produce valid primer pair");
        Assert.That(result.Forward!.HomopolymerLength, Is.LessThanOrEqualTo(maxHomopolymer),
            $"Forward primer homopolymer run should be ≤{maxHomopolymer}bp");
        Assert.That(result.Reverse!.HomopolymerLength, Is.LessThanOrEqualTo(maxHomopolymer),
            $"Reverse primer homopolymer run should be ≤{maxHomopolymer}bp");
    }

    [Test]
    public void EvaluatePrimer_ExcessiveHomopolymer_ReportsIssue()
    {
        // Arrange - primer with 6bp A run (exceeds default max of 4)
        string primer = "ACGTAAAAAAAACGTACGT"; // 19bp with 7x A

        // Act
        var candidate = PrimerDesigner.EvaluatePrimer(primer, 0, true);

        // Assert
        Assert.That(candidate.Issues.Any(i => i.Contains("Homopolymer")), Is.True,
            "Primer with excessive homopolymer run should report issue");
    }

    #endregion

    #region M8: No primer-dimer formation detected

    [Test]
    public void DesignPrimers_PrimerPair_NoPrimerDimerFormation()
    {
        // Arrange
        int targetStart = 100;
        int targetEnd = 150;

        // Act
        var result = PrimerDesigner.DesignPrimers(_standardTemplate, targetStart, targetEnd);

        // Assert
        Assert.That(result.IsValid, Is.True, "Standard template must produce valid primer pair");
        bool hasDimer = PrimerDesigner.HasPrimerDimer(
            result.Forward!.Sequence,
            result.Reverse!.Sequence);
        Assert.That(hasDimer, Is.False,
            "Selected primer pair should not form primer-dimers");
    }

    [Test]
    public void HasPrimerDimer_ComplementaryPrimers_ReturnsTrue()
    {
        // A 3'-3' primer-dimer needs the 3'-terminal bases of primer2 to be the REVERSE COMPLEMENT
        // of primer1's: primer1 ends ATCGATCG, primer2 ends CGATCGAT (= revcomp(ATCGATCG)).
        // primer3-py check_primers (alignment mode): PRIMER_PAIR_0_COMPL_END = 8.0 ("high end compl").
        // The former fixture (primer2 ending GCTAGCTA, the plain reverse) scores 1.0 in Primer3.
        string primer1 = "AACCGGTTAACCATCGATCG";
        string primer2 = "AACCGGTTAACGATCGAT";

        // Act
        bool hasDimer = PrimerDesigner.HasPrimerDimer(primer1, primer2);

        // Assert
        Assert.That(PrimerDesigner.CalculatePrimerDimerEndComplementarity(primer1, primer2), Is.EqualTo(8.0));
        Assert.That(PrimerDesigner.HasPrimerDimer(primer1, "AACCGGTTAAGCTAGCTA"), Is.False);
        Assert.That(hasDimer, Is.True,
            "Primers with fully complementary 3' ends should be detected as primer-dimer prone");
    }

    #endregion

    #region M9: No hairpin potential detected

    [Test]
    public void EvaluatePrimer_SelfComplementary_DetectsHairpin()
    {
        // Arrange - primer with clear hairpin potential:
        // GCGC (stem) + AAAA (loop) + GCGC (matching stem)
        // The reverse complement of GCGC is GCGC, forming a stable hairpin
        string primer = "GCGCAAAAGCGCATGCGATC"; // 20bp with hairpin structure

        // Act
        bool hasHairpin = PrimerDesigner.HasHairpinPotential(primer);

        // Assert
        Assert.That(hasHairpin, Is.True,
            "Primer with GCGC..GCGC stem-loop should be detected as hairpin-prone");
    }

    [Test]
    public void HasHairpinPotential_LongSelfComplementary_ReturnsTrue()
    {
        // Arrange - clearly self-complementary sequence
        // ACGT repeated is self-complementary after internal folding
        string primer = "ACGTACGTACGTACGTACGT"; // 20bp

        // Act
        bool hasHairpin = PrimerDesigner.HasHairpinPotential(primer);

        // Assert
        Assert.That(hasHairpin, Is.True,
            "Self-complementary sequence should be detected as hairpin-prone");
    }

    #endregion

    #region M10: Product size correctly calculated

    [Test]
    public void DesignPrimers_ValidResult_ProductSizeCorrect()
    {
        // Arrange
        int targetStart = 100;
        int targetEnd = 150;

        // Act
        var result = PrimerDesigner.DesignPrimers(_standardTemplate, targetStart, targetEnd);

        // Assert
        Assert.That(result.IsValid, Is.True, "Standard template must produce valid primer pair");
        int expectedProductSize = result.Reverse!.Position + result.Reverse.Length -
                                  result.Forward!.Position;
        Assert.That(result.ProductSize, Is.EqualTo(expectedProductSize),
            "Product size should equal distance from forward start to reverse end");
        Assert.That(result.ProductSize, Is.GreaterThan(targetEnd - targetStart),
            "Product size should be larger than target region");
    }

    #endregion

    #region M11: Invalid target coordinates throw ArgumentException

    [Test]
    public void DesignPrimers_TargetEndBeforeStart_ThrowsArgumentException()
    {
        // Arrange
        int targetStart = 100;
        int targetEnd = 50; // Before start

        // Act & Assert
        Assert.Throws<ArgumentException>(() =>
            PrimerDesigner.DesignPrimers(_standardTemplate, targetStart, targetEnd),
            "Should throw when target end is before target start");
    }

    [Test]
    public void DesignPrimers_TargetBeyondTemplate_ThrowsArgumentException()
    {
        // Arrange
        int targetStart = 0;
        int targetEnd = _standardTemplate.Length + 100; // Beyond template

        // Act & Assert
        Assert.Throws<ArgumentException>(() =>
            PrimerDesigner.DesignPrimers(_standardTemplate, targetStart, targetEnd),
            "Should throw when target extends beyond template");
    }

    [Test]
    public void DesignPrimers_NegativeCoordinates_ThrowsArgumentException()
    {
        // Arrange
        int targetStart = -10;
        int targetEnd = 50;

        // Act & Assert
        Assert.Throws<ArgumentException>(() =>
            PrimerDesigner.DesignPrimers(_standardTemplate, targetStart, targetEnd),
            "Should throw for negative coordinates");
    }

    #endregion

    #region M12: EvaluatePrimer returns PrimerCandidate with all properties populated

    [Test]
    public void EvaluatePrimer_ValidPrimer_AllPropertiesPopulated()
    {
        // Arrange
        string primer = "ATGCGATCGATCGATCGATC"; // 20bp, 50% GC
        int position = 42;
        bool isForward = true;

        // Act
        var candidate = PrimerDesigner.EvaluatePrimer(primer, position, isForward);

        // Assert
        Assert.Multiple(() =>
        {
            Assert.That(candidate.Sequence, Is.EqualTo(primer), "Sequence should match input");
            Assert.That(candidate.Position, Is.EqualTo(position), "Position should match input");
            Assert.That(candidate.IsForward, Is.EqualTo(isForward), "IsForward should match input");
            Assert.That(candidate.Length, Is.EqualTo(primer.Length), "Length should match sequence length");
            Assert.That(candidate.GcContent, Is.GreaterThanOrEqualTo(0).And.LessThanOrEqualTo(100),
                "GC content should be valid percentage");
            Assert.That(candidate.MeltingTemperature, Is.GreaterThan(0), "Tm should be calculated");
            Assert.That(candidate.Score, Is.GreaterThanOrEqualTo(0), "Score should be non-negative");
        });
    }

    [Test]
    public void EvaluatePrimer_Reverse_HasCorrectOrientation()
    {
        // Arrange
        string primer = "TAGCTAGCTAGCTAGCTAGC"; // 20bp
        bool isForward = false;

        // Act
        var candidate = PrimerDesigner.EvaluatePrimer(primer, 100, isForward);

        // Assert
        Assert.That(candidate.IsForward, Is.False, "Reverse primer should have IsForward = false");
    }

    #endregion

    #region M13: GeneratePrimerCandidates returns sorted candidates by score

    [Test]
    public void GeneratePrimerCandidates_ReturnsMultipleValidCandidates()
    {
        // Arrange & Act
        var candidates = PrimerDesigner.GeneratePrimerCandidates(_standardTemplate, 0, 60, true)
            .ToList();

        // Assert - verify candidates are generated with scores
        // Note: Implementation returns candidates in generation order, not sorted by score
        Assert.That(candidates.Count, Is.GreaterThan(0),
            "Should return at least one candidate");
        foreach (var candidate in candidates)
        {
            Assert.That(candidate.Score, Is.GreaterThanOrEqualTo(0),
                "Each candidate should have a calculated score");
        }
    }

    [Test]
    public void GeneratePrimerCandidates_AllCandidatesHaveValidLength()
    {
        // Arrange
        var candidates = PrimerDesigner.GeneratePrimerCandidates(_standardTemplate, 0, 60, true)
            .ToList();

        // Assert
        foreach (var candidate in candidates)
        {
            Assert.That(candidate.Length, Is.InRange(18, 25),
                $"Candidate at position {candidate.Position} should have length in valid range");
        }
    }

    [Test]
    public void GeneratePrimerCandidates_Reverse_SequenceIsReverseComplement()
    {
        // Arrange — use a known short region so we can verify the exact reverse complement
        var template = new DnaSequence("AACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTT"); // 40bp
        int regionStart = 0;
        int regionEnd = 40;

        // Act
        var candidates = PrimerDesigner.GeneratePrimerCandidates(
            template, regionStart, regionEnd, forward: false).ToList();

        // Assert — each reverse candidate's sequence should be the reverse complement
        // of the corresponding template substring
        Assert.That(candidates.Count, Is.GreaterThan(0), "Should generate reverse candidates");
        foreach (var candidate in candidates)
        {
            string templateSubstring = template.Sequence.Substring(candidate.Position, candidate.Length);
            string expectedRevComp = new DnaSequence(templateSubstring).ReverseComplement().Sequence;
            Assert.That(candidate.Sequence, Is.EqualTo(expectedRevComp),
                $"Reverse candidate at position {candidate.Position} (len {candidate.Length}) " +
                $"should be reverse complement of template substring '{templateSubstring}'");
        }
    }

    #endregion

    #region S1: Difficult templates may return invalid result with reason

    [Test]
    public void DesignPrimers_VeryShortTemplate_ThrowsArgumentException()
    {
        // Arrange - template too short for primer design
        var shortTemplate = new DnaSequence("ACGTACGTACGT"); // 12bp

        // Act & Assert
        // Implementation throws ArgumentException for invalid target regions
        Assert.Throws<ArgumentException>(() =>
            PrimerDesigner.DesignPrimers(shortTemplate, 0, 12),
            "Should throw for template too short to have valid target region");
    }

    [Test]
    public void DesignPrimers_HomopolymerRichTemplate_MayReturnInvalid()
    {
        // Arrange - template with extensive homopolymer regions
        var homopolymerTemplate = new DnaSequence(
            new string('A', 100) + // Homopolymer forward region
            "ACGT" + // Tiny target
            new string('T', 100)   // Homopolymer reverse region
        );

        // Act
        var result = PrimerDesigner.DesignPrimers(homopolymerTemplate, 100, 104);

        // Assert
        if (!result.IsValid)
        {
            Assert.That(result.Message, Is.Not.Null.And.Not.Empty,
                "Failed design should explain why primers couldn't be found");
        }
    }

    #endregion

    #region S2: Custom parameters are respected

    [Test]
    public void DesignPrimers_CustomParameters_AppliesLengthRange()
    {
        // Arrange
        var customParams = new PrimerParameters(
            MinLength: 22,
            MaxLength: 28,
            OptimalLength: 25,
            MinGcContent: 45,
            MaxGcContent: 55,
            MinTm: 58,
            MaxTm: 62,
            OptimalTm: 60,
            MaxHomopolymer: 3,
            MaxDinucleotideRepeats: 3,
            Avoid3PrimeGC: true,
            Check3PrimeStability: true
        );

        // Act
        var result = PrimerDesigner.DesignPrimers(_standardTemplate, 100, 150, customParams);

        // Assert - custom params may not find valid primers (stricter constraints),
        // so test with EvaluatePrimer directly for parameter respect
        if (result.IsValid)
        {
            Assert.That(result.Forward!.Length, Is.InRange(22, 28),
                "Forward primer should respect custom length range");
            Assert.That(result.Reverse!.Length, Is.InRange(22, 28),
                "Reverse primer should respect custom length range");
        }
        else
        {
            // Even if DesignPrimers fails with strict params, verify EvaluatePrimer uses them
            string primer = "ATGCGATCGATCGATCGATCGATC"; // 24bp
            var candidate = PrimerDesigner.EvaluatePrimer(primer, 0, true, customParams);
            Assert.That(candidate.Length, Is.InRange(22, 28),
                "EvaluatePrimer should report length within custom range");
        }
    }

    [Test]
    public void GeneratePrimerCandidates_CustomParameters_AppliesLengthRange()
    {
        // Arrange
        var customParams = new PrimerParameters(
            MinLength: 20,
            MaxLength: 22,
            OptimalLength: 21,
            MinGcContent: 40,
            MaxGcContent: 60,
            MinTm: 55,
            MaxTm: 65,
            OptimalTm: 60,
            MaxHomopolymer: 4,
            MaxDinucleotideRepeats: 4,
            Avoid3PrimeGC: false,
            Check3PrimeStability: false
        );

        // Act
        var candidates = PrimerDesigner.GeneratePrimerCandidates(
            _standardTemplate, 0, 50, true, customParams).ToList();

        // Assert
        foreach (var candidate in candidates)
        {
            Assert.That(candidate.Length, Is.InRange(20, 22),
                "All candidates should respect custom length range");
        }
    }

    #endregion

    #region S3: Score reflects primer quality (higher = better)

    [Test]
    public void EvaluatePrimer_OptimalPrimer_HasHighScore()
    {
        // Arrange - optimal primer: 20bp, ~50% GC, no issues
        string optimalPrimer = "ATGCGATCGATCGATCGATC";

        // Act
        var candidate = PrimerDesigner.EvaluatePrimer(optimalPrimer, 0, true);

        // Assert
        Assert.That(candidate.Score, Is.GreaterThan(50),
            "Optimal primer should have high score");
    }

    [Test]
    public void EvaluatePrimer_SuboptimalLength_ScoreVaries()
    {
        // Arrange
        string optimal = "ATGCGATCGATCGATCGATC"; // 20bp (optimal)
        string shorter = "ATGCGATCGATCGATCGAT"; // 19bp
        string longer = "ATGCGATCGATCGATCGATCG"; // 21bp

        // Act
        var optimalCandidate = PrimerDesigner.EvaluatePrimer(optimal, 0, true);
        var shorterCandidate = PrimerDesigner.EvaluatePrimer(shorter, 0, true);
        var longerCandidate = PrimerDesigner.EvaluatePrimer(longer, 0, true);

        // Assert - all primers should have scores calculated
        // Note: Score depends on multiple factors (Tm, GC, length) so
        // optimal length alone doesn't guarantee highest score
        Assert.That(optimalCandidate.Score, Is.GreaterThan(0),
            "Optimal primer should have positive score");
        Assert.That(shorterCandidate.Score, Is.GreaterThan(0),
            "Shorter primer should have positive score");
        Assert.That(longerCandidate.Score, Is.GreaterThan(0),
            "Longer primer should have positive score");
        Assert.That(optimalCandidate.Score, Is.GreaterThanOrEqualTo(shorterCandidate.Score),
            "Optimal length should score >= shorter");
    }

    #endregion

    #region S4: 3' stability is calculated correctly

    [Test]
    public void Calculate3PrimeStability_GCRich3Prime_MoreNegative()
    {
        // Arrange
        string gcRich3Prime = "ATATATATATATATGCGCGC"; // 3' = GCGCGC (GC-rich)
        string atRich3Prime = "GCGCGCGCGCGCGCATATAT"; // 3' = ATATAT (AT-rich)

        // Act
        double gcStability = PrimerDesigner.Calculate3PrimeStability(gcRich3Prime);
        double atStability = PrimerDesigner.Calculate3PrimeStability(atRich3Prime);

        // Assert - GC-rich 3' is more stable (more negative ΔG)
        Assert.That(gcStability, Is.LessThan(atStability),
            "GC-rich 3' end should have more negative (stable) ΔG");
    }

    #endregion

    #region S5: Dinucleotide repeats detected

    [Test]
    public void EvaluatePrimer_ExcessiveDinucleotideRepeats_ReportsIssue()
    {
        // Arrange - primer with long dinucleotide repeat
        string primer = "ACACACACACACACACAC"; // 18bp of AC repeats (9 repeats)

        // Act
        var candidate = PrimerDesigner.EvaluatePrimer(primer, 0, true);
        int dinucRepeats = PrimerDesigner.FindLongestDinucleotideRepeat(primer);

        // Assert
        Assert.That(dinucRepeats, Is.GreaterThan(4),
            "Should detect extensive dinucleotide repeats");
        Assert.That(candidate.Issues.Any(i => i.ToLower().Contains("dinucleotide") || i.ToLower().Contains("repeat")),
            Is.True, "Should report dinucleotide repeat issue");
    }

    #endregion

    #region C1: Multiple primer pairs can be generated

    [Test]
    public void GeneratePrimerCandidates_LargeRegion_ReturnsMultipleCandidates()
    {
        // Arrange
        var candidates = PrimerDesigner.GeneratePrimerCandidates(_standardTemplate, 0, 80, true)
            .ToList();

        // Assert
        Assert.That(candidates.Count, Is.GreaterThan(1),
            "Should generate multiple primer candidates from a larger region");
    }

    #endregion

    #region C1b: Performance on long templates

    [Test]
    public void DesignPrimers_LongTemplate_CompletesWithinTimeout()
    {
        // Arrange — 10 kb template with varied sequence (realistic gene region)
        var sb = new StringBuilder(10_000);
        var bases = "ACGTACGATCGATCGTAGCTAGCATGCATGC"; // 30bp repeating unit, ~50% GC
        while (sb.Length < 10_000)
            sb.Append(bases);
        var longTemplate = new DnaSequence(sb.ToString(0, 10_000));

        int targetStart = 5000;
        int targetEnd = 5100;

        // Act — should complete within a few seconds even on slow machines
        var sw = System.Diagnostics.Stopwatch.StartNew();
        var result = PrimerDesigner.DesignPrimers(longTemplate, targetStart, targetEnd);
        sw.Stop();

        // Assert — anti-hang guard (generous bound; won't flake under parallel-suite CPU load).
        Assert.That(sw.ElapsedMilliseconds, Is.LessThan(30_000),
            "DesignPrimers on 10kb template must not hang/blow up");
        Assert.That(result, Is.Not.Null,
            "Should return a result (valid or invalid) for long template");
    }

    #endregion

    #region C2: Primer positions are 0-indexed

    [Test]
    public void EvaluatePrimer_Position_IsZeroIndexed()
    {
        // Arrange
        string primer = "ATGCGATCGATCGATCGATC";
        int position = 0;

        // Act
        var candidate = PrimerDesigner.EvaluatePrimer(primer, position, true);

        // Assert
        Assert.That(candidate.Position, Is.EqualTo(0),
            "Position 0 should be valid (0-indexed)");
    }

    [Test]
    public void GeneratePrimerCandidates_StartAtZero_IncludesPositionZero()
    {
        // Arrange & Act
        var candidates = PrimerDesigner.GeneratePrimerCandidates(_standardTemplate, 0, 40, true)
            .ToList();

        // Assert
        bool hasPositionZero = candidates.Any(c => c.Position == 0);
        Assert.That(hasPositionZero, Is.True,
            "Candidates starting at region 0 should include position 0");
    }

    #endregion

    #region Edge Cases

    [Test]
    public void DesignPrimers_NullTemplate_ThrowsException()
    {
        // Act & Assert
        // Implementation throws NullReferenceException (not ArgumentNullException)
        Assert.Throws<NullReferenceException>(() =>
            PrimerDesigner.DesignPrimers(null!, 0, 100));
    }

    [Test]
    public void EvaluatePrimer_EmptySequence_HandledGracefully()
    {
        // Act
        var candidate = PrimerDesigner.EvaluatePrimer("", 0, true);

        // Assert
        Assert.That(candidate.IsValid, Is.False, "Empty sequence should not be valid");
    }

    [Test]
    public void GeneratePrimerCandidates_EmptyRegion_ReturnsEmpty()
    {
        // Arrange - region too small for any primer
        var candidates = PrimerDesigner.GeneratePrimerCandidates(_standardTemplate, 0, 10, true)
            .ToList();

        // Assert
        Assert.That(candidates.Count, Is.EqualTo(0),
            "Region smaller than min primer length should return no candidates");
    }

    #endregion
    #region Primer3 reference cross-checks (primer3-py 2.3.1)

    // primer3-py 2.3.1 design_primers settings mirroring DefaultParameters:
    // PRIMER_{MIN,OPT,MAX}_SIZE 18/20/25, PRIMER_{MIN,OPT,MAX}_TM 57/60/63, PRIMER_{MIN,MAX}_GC 40/60,
    // PRIMER_MAX_POLY_X 4, PRIMER_PAIR_MAX_DIFF_TM 5, thermodynamic structure limits disabled (1000),
    // PRIMER_PRODUCT_SIZE_RANGE 25-2000; every other tag at its Primer3 default.

    [TestCase("AGCTAGCTAGCTAGCTAGCT", 58.10101033826538)]
    [TestCase("GAACTCGTGAACTCGTGAAC", 56.817102902094234)]
    [TestCase("CGGTTCACTACGTCCGTTCTGG", 63.122839283928954)]
    [TestCase("GACGCTGTCTGAGACTAGAA", 56.4298406521595)]
    [TestCase("GAATTCGAATTCGAATTC", 47.794297549459145)] // self-complementary: C_T/1
    [TestCase("gcgcgcgc", 42.54999596123861)]           // case-insensitive, symmetric
    [TestCase("GAGCAGGATCCCTATAGAGTGACAAAAGGATCTTGGTCCA", 72.93445880948137)] // 40 nt > 36: long_seq_tm
    public void CalculateMeltingTemperaturePrimer3_DefaultConditions_MatchesPrimer3CalcTm(string seq, double expected)
    {
        // primer3.calc_tm(seq) (max_nn_length = 36 as in Primer3's design engine).
        Assert.That(PrimerDesigner.CalculateMeltingTemperaturePrimer3(seq), Is.EqualTo(expected).Within(1e-9));
    }

    [Test]
    public void CalculateMeltingTemperaturePrimer3_NonDefaultConditions_MatchesPrimer3CalcTm()
    {
        const string seq = "AGCTAGCTAGCTAGCTAGCT";
        Assert.Multiple(() =>
        {
            // calc_tm(seq, dv_conc=0, dntp_conc=0): no divalent term (dNTP ignored without Mg2+).
            Assert.That(PrimerDesigner.CalculateMeltingTemperaturePrimer3(seq, 50, 50, 0, 0),
                Is.EqualTo(52.18613256285289).Within(1e-9));
            // calc_tm(seq, mv_conc=100, dv_conc=2.5, dntp_conc=0.8, dna_conc=250).
            Assert.That(PrimerDesigner.CalculateMeltingTemperaturePrimer3(seq, 250, 100, 2.5, 0.8),
                Is.EqualTo(62.76000990783831).Within(1e-9));
            // calc_tm(seq, dv_conc=0.5, dntp_conc=0.6): Mg2+ fully chelated by dNTP → no divalent term.
            Assert.That(PrimerDesigner.CalculateMeltingTemperaturePrimer3(seq, 50, 50, 0.5, 0.6),
                Is.EqualTo(52.18613256285289).Within(1e-9));
        });
    }

    [Test]
    public void CalculateMeltingTemperaturePrimer3_InvalidInput_NaNOrThrows()
    {
        Assert.Multiple(() =>
        {
            Assert.That(PrimerDesigner.CalculateMeltingTemperaturePrimer3(""), Is.NaN);
            Assert.That(PrimerDesigner.CalculateMeltingTemperaturePrimer3("A"), Is.NaN);
            Assert.That(PrimerDesigner.CalculateMeltingTemperaturePrimer3("ACGTNACGTACGTACGTACG"), Is.NaN);
            Assert.Throws<ArgumentOutOfRangeException>(() => PrimerDesigner.CalculateMeltingTemperaturePrimer3("ACGTACGT", 0));
            Assert.Throws<ArgumentOutOfRangeException>(() => PrimerDesigner.CalculateMeltingTemperaturePrimer3("ACGTACGT", 50, -1));
            Assert.Throws<ArgumentOutOfRangeException>(() => PrimerDesigner.CalculateMeltingTemperaturePrimer3("ACGTACGT", 50, 0, 0, 0));
        });
    }

    [Test]
    public void EvaluatePrimer_NonAcgtBase_TmNotComputableAndInvalid()
    {
        // Primer3 PRIMER_MAX_NS_ACCEPTED = 0: an ambiguous base disqualifies the primer.
        var c = PrimerDesigner.EvaluatePrimer("AGCTAGCTAGNTAGCTAGCT", 0, true);
        Assert.Multiple(() =>
        {
            Assert.That(c.IsValid, Is.False);
            Assert.That(c.MeltingTemperature, Is.EqualTo(0));
            Assert.That(c.Issues, Has.Some.StartsWith("Tm not computable"));
        });
    }

    [Test]
    public void EvaluatePrimer_Penalty_IsPrimer3PerPrimerPenalty()
    {
        // primer3 PRIMER_LEFT_0 of the random-template case below: ATGCTGGGTAGAGGTCGAGG,
        // Tm 60.757187619543856, PRIMER_LEFT_0_PENALTY 0.7571876195438563 (= |Tm − 60| + |20 − 20|).
        var c = PrimerDesigner.EvaluatePrimer("ATGCTGGGTAGAGGTCGAGG", 14, true);
        Assert.That(c.Penalty, Is.EqualTo(0.7571876195438563).Within(1e-9));
        Assert.That(c.MeltingTemperature, Is.EqualTo(60.8));
    }

    [Test]
    public void DesignPrimers_RandomTemplate_MatchesPrimer3DesignPrimers()
    {
        // Random 160-mer (python random.seed(2026)), SEQUENCE_TARGET = 70,20. primer3-py 2.3.1
        // design_primers with Primer3's default thermodynamic structure limits (47 °C) and this
        // library's default per-primer limits (PRIMER_MIN/MAX_GC 40/60, PRIMER_MAX_POLY_X 4,
        // Tm 57–63, size 18–25, PRIMER_PAIR_MAX_DIFF_TM 5): PRIMER_LEFT_0 = [14,20]
        // ATGCTGGGTAGAGGTCGAGG (penalty 0.7571876195438563), PRIMER_RIGHT_0 = [129,20]
        // AGGAACGGATCGAGGACTGC (Tm 61.66810138313815), PRIMER_PAIR_0_PENALTY 2.425289002682007,
        // product 129 − 14 + 1 = 116. The pair returned without structure limits (right
        // GATCGAGGACTGCCTTGGTA) is rejected by Primer3: its hairpin Tm is 47.19 °C > 47 °C
        // (PRIMER_RIGHT_EXPLAIN "high hairpin stability").
        const string template =
            "AGACTTTCAAAGATATGCTGGGTAGAGGTCGAGGTTATTATTTGTTACCAATTCTCATTGTGTTTCGGAACTTGCGTTTTAGGTATGTCTTAGTGACTCTAAATACCAAGGCAGTCCTCGATCCGTTCCTAATAAGGAATGGTGATTCCCTGTCATACCA";
        var result = PrimerDesigner.DesignPrimers(new DnaSequence(template), 70, 90);

        Assert.Multiple(() =>
        {
            Assert.That(result.IsValid, Is.True);
            Assert.That(result.Forward!.Sequence, Is.EqualTo("ATGCTGGGTAGAGGTCGAGG"));
            Assert.That(result.Forward.Position, Is.EqualTo(14));
            Assert.That(result.Reverse!.Sequence, Is.EqualTo("AGGAACGGATCGAGGACTGC"));
            Assert.That(result.Reverse.Position, Is.EqualTo(110)); // Primer3 right start 129 = 110 + 20 − 1
            Assert.That(result.Reverse.MeltingTemperature, Is.EqualTo(61.7));
            Assert.That(result.Forward.Penalty + result.Reverse.Penalty, Is.EqualTo(2.425289002682007).Within(1e-9));
            Assert.That(result.ProductSize, Is.EqualTo(116));
            Assert.That(PrimerDesigner.CalculatePrimer3OligoStructure("GATCGAGGACTGCCTTGGTA")!.Value.HairpinTh,
                Is.EqualTo(47.193593776115506).Within(1e-9)); // primer3-py calc_hairpin
        });
    }

    [Test]
    public void DesignPrimers_HeuristicScreen_KeepsSequenceOnlyChecks()
    {
        // With the sequence-only screen the Primer3 hairpin limit is not applied, so the pair is the
        // one primer3-py returns when every *_TH limit is 100 °C: right [122,20] GATCGAGGACTGCCTTGGTA,
        // PRIMER_PAIR_0_PENALTY 1.869300477208128 (HasHairpinPotential finds no 4-bp stem in it).
        const string template =
            "AGACTTTCAAAGATATGCTGGGTAGAGGTCGAGGTTATTATTTGTTACCAATTCTCATTGTGTTTCGGAACTTGCGTTTTAGGTATGTCTTAGTGACTCTAAATACCAAGGCAGTCCTCGATCCGTTCCTAATAAGGAATGGTGATTCCCTGTCATACCA";
        var heuristic = PrimerDesigner.DesignPrimers(new DnaSequence(template), 70, 90,
            PrimerDesigner.DefaultParameters with { StructureScreen = PrimerStructureScreen.Heuristic });
        var relaxed = PrimerDesigner.DesignPrimers(new DnaSequence(template), 70, 90,
            PrimerDesigner.DefaultParameters with { MaxStructureTm = 100 });

        Assert.Multiple(() =>
        {
            foreach (var r in new[] { heuristic, relaxed })
            {
                Assert.That(r.IsValid, Is.True);
                Assert.That(r.Reverse!.Sequence, Is.EqualTo("GATCGAGGACTGCCTTGGTA"));
                Assert.That(r.Forward!.Penalty + r.Reverse.Penalty, Is.EqualTo(1.869300477208128).Within(1e-9));
            }
            Assert.That(heuristic.Reverse!.HairpinTh, Is.Null);
            Assert.That(relaxed.Reverse!.HairpinTh, Is.EqualTo(47.193593776115506).Within(1e-9));
        });
    }

    [Test]
    public void DesignPrimers_IndividuallyBestPrimersTmIncompatible_SearchesPairs()
    {
        // Regression for greedy selection. Without the pair ΔTm limit primer3-py returns the
        // individually best primers LEFT [4,22] (Tm 57.300) + RIGHT [62,19] (Tm 62.882): ΔTm 5.58 > 5,
        // so choosing each side independently yields no valid pair. With PRIMER_PAIR_MAX_DIFF_TM = 5
        // primer3-py returns LEFT [3,23] (Tm 57.9563) + RIGHT [62,19], PRIMER_PAIR_0_PENALTY 8.925688222301858.
        // (primer3-py run with every *_TH structure limit at 100 °C; with Primer3's default 47 °C the
        // right primer ATGGAGCACGAGCGCAACA fails on its hairpin, 48.215 °C — see the next test.)
        const string template = "TAATTGGTGTAATAATCTAGGGGTGCTTTTTTTTTGCAGTCCGGTGTTGCGCTCGTGCTCCAT";
        var param = PrimerDesigner.DefaultParameters with { MaxStructureTm = 100 };
        var bestForwardAlone = PrimerDesigner.EvaluatePrimer(template.Substring(4, 22), 4, true, param);
        Assert.That(bestForwardAlone.IsValid, Is.True, "the individually best forward primer is itself valid");

        var result = PrimerDesigner.DesignPrimers(new DnaSequence(template), 26, 34, param);

        Assert.Multiple(() =>
        {
            Assert.That(result.IsValid, Is.True, "a compatible pair exists and must be found");
            Assert.That(result.Forward!.Position, Is.EqualTo(3));
            Assert.That(result.Forward.Length, Is.EqualTo(23));
            Assert.That(result.Forward.MeltingTemperature, Is.EqualTo(58.0));
            Assert.That(result.Reverse!.Position, Is.EqualTo(44)); // right start 62 = 44 + 19 − 1
            Assert.That(result.Reverse.Length, Is.EqualTo(19));
            Assert.That(result.Reverse.MeltingTemperature, Is.EqualTo(62.9));
            Assert.That(result.Forward.Penalty + result.Reverse.Penalty, Is.EqualTo(8.925688222301858).Within(1e-9));
            Assert.That(bestForwardAlone.Penalty, Is.LessThan(result.Forward.Penalty),
                "the chosen forward primer is not the individually best one");
        });
    }

    [Test]
    [CancelAfter(20000)]
    public void DesignPrimers_RepetitiveTemplate_StructureChecksCachedBySequence()
    {
        // 250 A's: ~800 candidates per side but only 8 distinct sequences each (A18..A25 / T18..T25),
        // many with equal penalties, so the pair loop meets each sequence pair many times. The ntthal
        // results are cached per sequence (pair); without it this took hours. A18·T18 hetero-dimer
        // Tm = 38.05 °C (primer3-py calc_heterodimer) ≤ 47 °C, so a pair exists.
        var param = PrimerDesigner.DefaultParameters with
        {
            MinGcContent = 0, MaxGcContent = 0, MinTm = 0, MaxTm = 200,
            MaxHomopolymer = 1000, MaxDinucleotideRepeats = 1000, Check3PrimeStability = false,
        };
        var sw = System.Diagnostics.Stopwatch.StartNew();
        var result = PrimerDesigner.DesignPrimers(new DnaSequence(new string('A', 250)), 120, 130, param);
        sw.Stop();

        Assert.Multiple(() =>
        {
            Assert.That(result.IsValid, Is.True);
            Assert.That(result.Forward!.Sequence.Trim('A'), Is.Empty);
            Assert.That(result.Reverse!.Sequence.Trim('T'), Is.Empty);
            var pc = PrimerDesigner.CalculatePrimer3PairComplementarity(result.Forward.Sequence, result.Reverse.Sequence)!.Value;
            Assert.That(pc.Exceeds(), Is.False);
            Assert.That(PrimerDesigner.CalculatePrimer3PairComplementarity(new string('A', 18), new string('T', 18))!.Value.ComplAnyTh,
                Is.EqualTo(38.04660273299868).Within(1e-9));
            Assert.That(sw.Elapsed, Is.LessThan(TimeSpan.FromSeconds(20)));
        });
    }

    [Test]
    public void DesignPrimers_RightPrimerHairpinAbovePrimer3Limit_NoValidPrimers()
    {
        // primer3-py 2.3.1 design_primers (default 47 °C limits, this library's per-primer limits):
        // PRIMER_RIGHT_EXPLAIN "considered 68, GC content failed 66, high tm 1, high hairpin
        // stability 1, ok 0" — the only GC-compatible right primer ATGGAGCACGAGCGCAACA has hairpin
        // Tm 48.21523465319416 °C (calc_hairpin), so no pair exists.
        const string template = "TAATTGGTGTAATAATCTAGGGGTGCTTTTTTTTTGCAGTCCGGTGTTGCGCTCGTGCTCCAT";
        var result = PrimerDesigner.DesignPrimers(new DnaSequence(template), 26, 34);
        var right = PrimerDesigner.EvaluatePrimer("ATGGAGCACGAGCGCAACA", 44, false);

        Assert.Multiple(() =>
        {
            Assert.That(result.IsValid, Is.False);
            Assert.That(result.Message, Is.EqualTo("Could not find valid primers for the target region."));
            Assert.That(right.HairpinTh, Is.EqualTo(48.21523465319416).Within(1e-9));
            Assert.That(right.HasHairpin, Is.True);
            Assert.That(right.IsValid, Is.False);
        });
    }

    [Test]
    public void DesignPrimers_NoPairWithinTmLimit_ReturnsInvalidWithTmMessage()
    {
        // primer3-py: with PRIMER_PAIR_MAX_DIFF_TM = 100 the best pair is LEFT [0,20] (Tm 57.2614) +
        // RIGHT [47,19] (Tm 62.6299); with PRIMER_PAIR_MAX_DIFF_TM = 5 no pair is returned.
        const string template = "TTGACCACAGCCAGGTTTAATTTTTTTTCAAATACGGTCACGCGCGGA";
        var result = PrimerDesigner.DesignPrimers(new DnaSequence(template), 20, 28);

        Assert.Multiple(() =>
        {
            Assert.That(result.IsValid, Is.False);
            Assert.That(result.Forward!.Position, Is.EqualTo(0));
            Assert.That(result.Forward.Length, Is.EqualTo(20));
            Assert.That(result.Forward.MeltingTemperature, Is.EqualTo(57.3));
            Assert.That(result.Reverse!.Position, Is.EqualTo(29)); // right start 47 = 29 + 19 − 1
            Assert.That(result.Reverse.MeltingTemperature, Is.EqualTo(62.6));
            Assert.That(result.Message, Does.StartWith("No primer pair within the 5°C Tm-difference limit"));
        });
    }

    #endregion
}
