namespace Seqeron.Genomics.Tests.Unit.Analysis;

/// <summary>
/// Evidence-based tests for RepeatFinder.FindMicrosatellites (REP-STR-001).
/// 
/// Sources:
/// - Wikipedia: Microsatellite (https://en.wikipedia.org/wiki/Microsatellite)
/// - Wikipedia: Trinucleotide repeat disorder
/// - Richard GF et al. (2008) MMBR
/// - Tóth G et al. (2000) Genome Research
/// </summary>
[TestFixture]
public class RepeatFinder_Microsatellite_Tests
{
    #region Repeat Type Detection (Evidence: Wikipedia - microsatellite classification)

    /// <summary>
    /// Evidence: Wikipedia states mononucleotide repeats are 1 bp units.
    /// Example: Poly-A tracts are common in genomes.
    /// </summary>
    [Test]
    public void FindMicrosatellites_MononucleotideRepeat_DetectsPolyATract()
    {
        var sequence = new DnaSequence("ACGTAAAAAACGT");

        var results = RepeatFinder.FindMicrosatellites(sequence, 1, 6, 3).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(results, Has.Count.EqualTo(1));
            Assert.That(results[0].RepeatUnit, Is.EqualTo("A"));
            Assert.That(results[0].RepeatCount, Is.EqualTo(6));
            Assert.That(results[0].Position, Is.EqualTo(4));
            Assert.That(results[0].RepeatType, Is.EqualTo(RepeatType.Mononucleotide));
            Assert.That(results[0].TotalLength, Is.EqualTo(6));
        });
    }

    /// <summary>
    /// Evidence: Wikipedia - "TATATATATA is a dinucleotide microsatellite"
    /// CA/AC repeats are among the most common in eukaryotic genomes.
    /// </summary>
    [Test]
    public void FindMicrosatellites_DinucleotideRepeat_DetectsCaRepeat()
    {
        var sequence = new DnaSequence("AAACACACACACAAA");

        var results = RepeatFinder.FindMicrosatellites(sequence, 2, 6, 3).ToList();

        Assert.Multiple(() =>
        {
            // One maximal period-2 run S[2..13) = ACACACACACA (11 bp): reported ONCE at its left end with
            // ⌊11/2⌋ = 5 complete copies; the rotation CA×5 at 3 is the same locus and is not re-reported
            // (MISA leftmost match / pytrf run start; independent brute-force maximal-run reference agrees).
            Assert.That(results, Has.Count.EqualTo(1));
            Assert.That(results[0].RepeatUnit, Is.EqualTo("AC"));
            Assert.That(results[0].RepeatCount, Is.EqualTo(5));
            Assert.That(results[0].Position, Is.EqualTo(2));
            Assert.That(results[0].TotalLength, Is.EqualTo(10));
            Assert.That(results[0].RepeatType, Is.EqualTo(RepeatType.Dinucleotide));
        });
    }

    /// <summary>
    /// Evidence: Wikipedia - exact quote: "TATATATATA is a dinucleotide microsatellite".
    /// Verifies detection using Wikipedia's canonical dinucleotide example.
    /// Source: https://en.wikipedia.org/wiki/Microsatellite (Structures, locations, and functions)
    /// </summary>
    [Test]
    public void FindMicrosatellites_WikipediaTataExample_DetectsDinucleotide()
    {
        // Exact Wikipedia example: "TATATATATA is a dinucleotide microsatellite"
        // TA×5 = 10 bp
        var sequence = new DnaSequence("TATATATATA");

        var results = RepeatFinder.FindMicrosatellites(sequence, 2, 2, 3).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(results, Has.Count.EqualTo(1));
            Assert.That(results[0].RepeatUnit, Is.EqualTo("TA"));
            Assert.That(results[0].RepeatCount, Is.EqualTo(5));
            Assert.That(results[0].RepeatType, Is.EqualTo(RepeatType.Dinucleotide));
            Assert.That(results[0].TotalLength, Is.EqualTo(10));
            Assert.That(results[0].Position, Is.EqualTo(0));
        });
    }

    /// <summary>
    /// Evidence: Wikipedia - exact quote: "GTCGTCGTCGTCGTC is a trinucleotide microsatellite".
    /// Verifies detection using Wikipedia's canonical trinucleotide example.
    /// Source: https://en.wikipedia.org/wiki/Microsatellite (Structures, locations, and functions)
    /// </summary>
    [Test]
    public void FindMicrosatellites_WikipediaGtcExample_DetectsTrinucleotide()
    {
        // Exact Wikipedia example: "GTCGTCGTCGTCGTC is a trinucleotide microsatellite"
        // GTC×5 = 15 bp
        var sequence = new DnaSequence("GTCGTCGTCGTCGTC");

        var results = RepeatFinder.FindMicrosatellites(sequence, 3, 3, 3).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(results, Has.Count.EqualTo(1));
            Assert.That(results[0].RepeatUnit, Is.EqualTo("GTC"));
            Assert.That(results[0].RepeatCount, Is.EqualTo(5));
            Assert.That(results[0].RepeatType, Is.EqualTo(RepeatType.Trinucleotide));
            Assert.That(results[0].TotalLength, Is.EqualTo(15));
            Assert.That(results[0].Position, Is.EqualTo(0));
        });
    }

    /// <summary>
    /// Evidence: Wikipedia Trinucleotide repeat disorder - CAG repeats cause Huntington's disease.
    /// HD: normal 6-35 repeats, pathogenic 36-250 repeats.
    /// </summary>
    [Test]
    public void FindMicrosatellites_TrinucleotideRepeat_DetectsCagExpansion()
    {
        // Simulate a pathogenic-range CAG expansion (5 repeats for test brevity)
        var sequence = new DnaSequence("ATGCAGCAGCAGCAGCAGTGA");

        var results = RepeatFinder.FindMicrosatellites(sequence, 3, 3, 3).ToList();

        Assert.Multiple(() =>
        {
            // The period-3 run is S[2..18) = G + (CAG)×5 (16 bp, starts at the G because S[2] = S[5]);
            // reported once at its left end as GCA×5 (MISA / pytrf report the run-start phase). The
            // rotation CAG×5 at 3 is the same locus and is not re-reported.
            Assert.That(results, Has.Count.EqualTo(1));
            Assert.That(results[0].RepeatUnit, Is.EqualTo("GCA"));
            Assert.That(results[0].RepeatCount, Is.EqualTo(5));
            Assert.That(results[0].Position, Is.EqualTo(2));
            Assert.That(results[0].TotalLength, Is.EqualTo(15));
            Assert.That(results[0].RepeatType, Is.EqualTo(RepeatType.Trinucleotide));
        });
    }

    /// <summary>
    /// Evidence: Wikipedia - forensic markers use tetra- and pentanucleotide repeats
    /// for higher accuracy and reduced PCR stutter.
    /// </summary>
    [Test]
    public void FindMicrosatellites_TetranucleotideRepeat_DetectsGataMarker()
    {
        // GATA is a common forensic marker pattern
        var sequence = new DnaSequence("AAGATAGATAGATAGATAAA");

        var results = RepeatFinder.FindMicrosatellites(sequence, 4, 4, 3).ToList();

        Assert.Multiple(() =>
        {
            // Period-4 run S[1..18) = A(GATA)×4 (17 bp, S[0]=A ≠ S[4]=T) → AGAT×4 at 1 (⌊17/4⌋ = 4),
            // reported once; the GATA×4 rotation at 2 is the same locus.
            Assert.That(results, Has.Count.EqualTo(1));
            Assert.That(results[0].RepeatUnit, Is.EqualTo("AGAT"));
            Assert.That(results[0].RepeatCount, Is.EqualTo(4));
            Assert.That(results[0].Position, Is.EqualTo(1));
            Assert.That(results[0].TotalLength, Is.EqualTo(16));
            Assert.That(results[0].RepeatType, Is.EqualTo(RepeatType.Tetranucleotide));
        });
    }

    /// <summary>
    /// Evidence: Wikipedia History - the first microsatellite was characterized in 1984
    /// by Weller, Jeffreys et al. as a polymorphic GGAT repeat in the human myoglobin gene.
    /// Source: https://en.wikipedia.org/wiki/Microsatellite (History)
    /// </summary>
    [Test]
    public void FindMicrosatellites_GgatMyoglobinRepeat_FirstEverCharacterizedMicrosatellite()
    {
        // GGAT is the historically first microsatellite ever described (Weller et al. 1984)
        var sequence = new DnaSequence("AAGGATGGATGGATGGATAA");

        var results = RepeatFinder.FindMicrosatellites(sequence, 4, 4, 3).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(results, Has.Count.EqualTo(1));
            Assert.That(results[0].RepeatUnit, Is.EqualTo("GGAT"));
            Assert.That(results[0].RepeatCount, Is.EqualTo(4));
            Assert.That(results[0].Position, Is.EqualTo(2));
            Assert.That(results[0].TotalLength, Is.EqualTo(16));
            Assert.That(results[0].RepeatType, Is.EqualTo(RepeatType.Tetranucleotide));
        });
    }

    /// <summary>
    /// Evidence: Wikipedia - microsatellites can be up to 6 bp (hexanucleotide).
    /// </summary>
    [Test]
    public void FindMicrosatellites_HexanucleotideRepeat_DetectsSixBpUnit()
    {
        // GAATTC is the EcoRI recognition site, used as hexanucleotide example
        var sequence = new DnaSequence("AAAGAATTCGAATTCGAATTCAAA");

        var results = RepeatFinder.FindMicrosatellites(sequence, 6, 6, 3).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(results, Has.Count.EqualTo(1));
            Assert.That(results[0].RepeatUnit, Is.EqualTo("GAATTC"));
            Assert.That(results[0].RepeatCount, Is.EqualTo(3));
            Assert.That(results[0].RepeatType, Is.EqualTo(RepeatType.Hexanucleotide));
        });
    }

    /// <summary>
    /// Verify RepeatType enum correctly maps unit lengths 1-6.
    /// Uses non-redundant repeat units (units that aren't just repetitions of smaller patterns).
    /// Tests for presence of correct RepeatType rather than exact count.
    /// </summary>
    [TestCase(1, "A", RepeatType.Mononucleotide)]
    [TestCase(2, "CA", RepeatType.Dinucleotide)]
    [TestCase(3, "CAG", RepeatType.Trinucleotide)]
    [TestCase(4, "GATA", RepeatType.Tetranucleotide)]
    [TestCase(5, "GATAC", RepeatType.Pentanucleotide)]
    [TestCase(6, "GAATTC", RepeatType.Hexanucleotide)]
    public void FindMicrosatellites_RepeatTypeClassification_MatchesUnitLength(int unitLength, string unit, RepeatType expectedType)
    {
        // Create a sequence with 5 repeats of the given unit
        var sequence = new DnaSequence(string.Concat(Enumerable.Repeat(unit, 5)));

        var results = RepeatFinder.FindMicrosatellites(sequence, unitLength, unitLength, 3).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(results, Has.Count.EqualTo(1), $"Should find exactly one repeat of {unit}");
            Assert.That(results[0].RepeatUnit, Is.EqualTo(unit), $"Unit should be {unit}");
            Assert.That(results[0].RepeatCount, Is.EqualTo(5), $"Count should be 5");
            Assert.That(results[0].Position, Is.EqualTo(0), $"Position should be 0");
            Assert.That(results[0].RepeatType, Is.EqualTo(expectedType), $"Type should be {expectedType}");
            Assert.That(results[0].TotalLength, Is.EqualTo(unitLength * 5), $"TotalLength should be {unitLength * 5}");
        });
    }

    #endregion

    #region Invariant Tests (Contract Validation)

    /// <summary>
    /// Invariant: TotalLength == RepeatUnit.Length × RepeatCount
    /// </summary>
    [Test]
    public void FindMicrosatellites_TotalLengthInvariant_AlwaysCorrect()
    {
        var sequence = new DnaSequence("AAAAAACGTCGTCGTACACACACATGATGATG");

        var results = RepeatFinder.FindMicrosatellites(sequence, 1, 6, 3).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(results, Has.Count.GreaterThanOrEqualTo(1), "Should find at least one repeat");
            foreach (var result in results)
            {
                Assert.That(result.TotalLength, Is.EqualTo(result.RepeatUnit.Length * result.RepeatCount),
                    $"TotalLength invariant violated for {result.RepeatUnit}");
            }
        });
    }

    /// <summary>
    /// Invariant: FullSequence == RepeatUnit repeated RepeatCount times
    /// </summary>
    [Test]
    public void FindMicrosatellites_FullSequenceInvariant_AlwaysCorrect()
    {
        var sequence = new DnaSequence("CAGCAGCAGCAGCAG");

        var results = RepeatFinder.FindMicrosatellites(sequence, 3, 3, 3).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(results, Has.Count.GreaterThanOrEqualTo(1));
            foreach (var result in results)
            {
                var expected = string.Concat(Enumerable.Repeat(result.RepeatUnit, result.RepeatCount));
                Assert.That(result.FullSequence, Is.EqualTo(expected),
                    $"FullSequence invariant violated for {result.RepeatUnit}");
            }
        });
    }

    /// <summary>
    /// Invariant: Position is within valid range [0, sequence.Length - TotalLength]
    /// </summary>
    [Test]
    public void FindMicrosatellites_PositionInvariant_WithinValidRange()
    {
        var sequence = new DnaSequence("ACGTAAAAAACGT");

        var results = RepeatFinder.FindMicrosatellites(sequence, 1, 6, 3).ToList();

        Assert.Multiple(() =>
        {
            foreach (var result in results)
            {
                Assert.That(result.Position, Is.GreaterThanOrEqualTo(0),
                    "Position should not be negative");
                Assert.That(result.Position, Is.LessThanOrEqualTo(sequence.Length - result.TotalLength),
                    "Position + TotalLength should not exceed sequence length");
            }
        });
    }

    /// <summary>
    /// Invariant: The sequence at the reported position matches FullSequence
    /// </summary>
    [Test]
    public void FindMicrosatellites_SequenceAtPosition_MatchesFullSequence()
    {
        var sequenceStr = "ACGTAAAAAACGTCGTCGTCGT";
        var sequence = new DnaSequence(sequenceStr);

        var results = RepeatFinder.FindMicrosatellites(sequence, 1, 6, 3).ToList();

        Assert.Multiple(() =>
        {
            foreach (var result in results)
            {
                var actualSequence = sequenceStr.Substring(result.Position, result.TotalLength);
                Assert.That(actualSequence, Is.EqualTo(result.FullSequence),
                    $"Sequence at position {result.Position} should match FullSequence");
            }
        });
    }

    /// <summary>
    /// Invariant: All results have RepeatCount >= minRepeats
    /// </summary>
    [TestCase(3)]
    [TestCase(4)]
    [TestCase(5)]
    public void FindMicrosatellites_MinRepeatsInvariant_AlwaysRespected(int minRepeats)
    {
        var sequence = new DnaSequence("AAAAAAAAAA"); // 10 A's

        var results = RepeatFinder.FindMicrosatellites(sequence, 1, 1, minRepeats).ToList();

        Assert.Multiple(() =>
        {
            foreach (var result in results)
            {
                Assert.That(result.RepeatCount, Is.GreaterThanOrEqualTo(minRepeats),
                    $"RepeatCount should be >= {minRepeats}");
            }
        });
    }

    #endregion

    #region Edge Cases

    /// <summary>
    /// Edge case: Empty sequence should return empty results.
    /// </summary>
    [Test]
    public void FindMicrosatellites_EmptySequence_ReturnsEmpty()
    {
        var results = RepeatFinder.FindMicrosatellites("", 1, 6, 3).ToList();

        Assert.That(results, Is.Empty);
    }

    /// <summary>
    /// Edge case: Sequence too short for any repeat.
    /// </summary>
    [Test]
    public void FindMicrosatellites_SequenceTooShort_ReturnsEmpty()
    {
        var results = RepeatFinder.FindMicrosatellites("AT", 2, 2, 3).ToList();

        Assert.That(results, Is.Empty);
    }

    /// <summary>
    /// Edge case: Exactly minRepeats should be detected.
    /// </summary>
    [Test]
    public void FindMicrosatellites_ExactlyMinRepeats_IsDetected()
    {
        var sequence = new DnaSequence("ATATAT"); // AT x 3 exactly

        var results = RepeatFinder.FindMicrosatellites(sequence, 2, 2, 3).ToList();

        Assert.That(results, Has.Count.EqualTo(1));
        Assert.That(results[0].RepeatCount, Is.EqualTo(3));
    }

    /// <summary>
    /// Edge case: Below minRepeats threshold should not be detected.
    /// </summary>
    [Test]
    public void FindMicrosatellites_BelowMinRepeats_NotDetected()
    {
        var sequence = new DnaSequence("ATAT"); // AT x 2, below minRepeats=3

        var results = RepeatFinder.FindMicrosatellites(sequence, 2, 2, 3).ToList();

        Assert.That(results, Is.Empty);
    }

    /// <summary>
    /// Edge case: Entire sequence is one repeat.
    /// </summary>
    [Test]
    public void FindMicrosatellites_EntireSequenceIsRepeat_CorrectCount()
    {
        var sequence = new DnaSequence("CAGCAGCAGCAGCAGCAGCAGCAGCAGCAG"); // CAG x 10

        var results = RepeatFinder.FindMicrosatellites(sequence, 3, 3, 3).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(results, Has.Count.EqualTo(1));
            Assert.That(results[0].RepeatCount, Is.EqualTo(10));
            Assert.That(results[0].Position, Is.EqualTo(0));
            Assert.That(results[0].TotalLength, Is.EqualTo(30));
        });
    }

    /// <summary>
    /// Edge case: Multiple different repeats in sequence.
    /// </summary>
    [Test]
    public void FindMicrosatellites_MultipleDifferentRepeats_FindsAll()
    {
        var sequence = new DnaSequence("AAAAAACGTCGTCGTACACACAC");

        var results = RepeatFinder.FindMicrosatellites(sequence, 1, 6, 3).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(results, Has.Count.EqualTo(3));
            Assert.That(results[0].RepeatUnit, Is.EqualTo("A"));
            Assert.That(results[0].RepeatCount, Is.EqualTo(6));
            Assert.That(results[0].Position, Is.EqualTo(0));
            Assert.That(results[1].RepeatUnit, Is.EqualTo("AC"));
            Assert.That(results[1].RepeatCount, Is.EqualTo(4));
            Assert.That(results[1].Position, Is.EqualTo(15));
            Assert.That(results[2].RepeatUnit, Is.EqualTo("CGT"));
            Assert.That(results[2].RepeatCount, Is.EqualTo(3));
            Assert.That(results[2].Position, Is.EqualTo(6));
        });
    }

    /// <summary>
    /// Case insensitivity: lowercase input should be handled.
    /// </summary>
    [Test]
    public void FindMicrosatellites_LowercaseInput_HandledCorrectly()
    {
        var results = RepeatFinder.FindMicrosatellites("cagcagcagcag", 3, 3, 3).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(results, Has.Count.EqualTo(1));
            Assert.That(results[0].RepeatUnit, Is.EqualTo("CAG"));
            Assert.That(results[0].RepeatCount, Is.EqualTo(4));
        });
    }

    /// <summary>
    /// M12: Redundant unit filtering — "ATAT" (4bp) should be reduced to "AT" (2bp) repeated more times.
    /// Implementation skips units that are repetitions of smaller patterns (e.g., ATAT = AT×2).
    /// </summary>
    [Test]
    public void FindMicrosatellites_RedundantUnitFiltering_ReportsSmallestUnit()
    {
        // ATATATAT = AT×4; searching unit range 2-4 should NOT report ATAT×2
        var sequence = new DnaSequence("ATATATAT");

        var results = RepeatFinder.FindMicrosatellites(sequence, 2, 4, 2).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(results, Has.Count.EqualTo(1), "Should find exactly one repeat (AT, not ATAT)");
            Assert.That(results[0].RepeatUnit, Is.EqualTo("AT"));
            Assert.That(results[0].RepeatCount, Is.EqualTo(4));
            Assert.That(results[0].Position, Is.EqualTo(0));
            Assert.That(results[0].TotalLength, Is.EqualTo(8));
            Assert.That(results[0].RepeatType, Is.EqualTo(RepeatType.Dinucleotide));
        });
    }

    /// <summary>
    /// S05: Non-standard characters (N) — DnaSequence rejects N (only ACGT valid); the string overload
    /// accepts N but never reports a unit containing a non-ACGT symbol: MISA searches <c>[acgt]{p}</c> motifs
    /// only and pytrf skips N, so a run of N (assembly gap) is not a microsatellite.
    /// </summary>
    [Test]
    public void FindMicrosatellites_NonStandardCharacterN_DnaSequenceRejectsStringOverloadSkipsN()
    {
        // DnaSequence constructor rejects N
        Assert.Throws<ArgumentException>(() => new DnaSequence("AAANNNAAACGT"));

        var results = RepeatFinder.FindMicrosatellites("AAANNNAAACGT", 1, 6, 3).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(results, Has.Count.EqualTo(2));
            Assert.That(results[0].RepeatUnit, Is.EqualTo("A"));
            Assert.That(results[0].RepeatCount, Is.EqualTo(3));
            Assert.That(results[0].Position, Is.EqualTo(0));
            Assert.That(results[1].RepeatUnit, Is.EqualTo("A"));
            Assert.That(results[1].RepeatCount, Is.EqualTo(3));
            Assert.That(results[1].Position, Is.EqualTo(6));
            Assert.That(RepeatFinder.FindMicrosatellites("NNNNNNNN", 1, 6, 3), Is.Empty);
            Assert.That(RepeatFinder.FindMicrosatellites("ANANANAN", 1, 6, 3), Is.Empty);
        });
    }

    /// <summary>
    /// Maximal-run semantics (each locus reported once per unit length). Expected values produced by an
    /// independent brute-force maximal-repetition reference (Kolpakov &amp; Kucherov 1999) and agreeing with
    /// MISA (leftmost regex <c>([ACGT]{p})\2{k-1,}</c>) and pytrf 1.5.0 <c>STRFinder</c> run-start output.
    /// </summary>
    [TestCase("ATATATA", 2, 2, 3, "0,AT,3")]              // partial trailing A is not a copy; TA×3 rotation at 1 not re-reported
    [TestCase("ATATATAT", 2, 2, 2, "0,AT,4")]             // TA×3 at 1 is a suffix of the same run
    [TestCase("CAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCA", 3, 3, 3, "0,CAG,10")]
    [TestCase("AAAAAACACACAC", 1, 2, 3, "0,A,6;5,AC,4")]  // AC run starts at 5 (S[4]=A ≠ S[6]=C)
    [TestCase("ACACACGCGCGC", 2, 2, 3, "0,AC,3;5,CG,3")]  // two distinct period-2 runs overlapping by 1 base
    [TestCase("AAGATAGATAGATAGATAAA", 1, 6, 3, "17,A,3;1,AGAT,4")]
    public void FindMicrosatellites_MaximalRuns_EachLocusReportedOnce(
        string sequence, int minUnit, int maxUnit, int minRepeats, string expected)
    {
        var results = RepeatFinder.FindMicrosatellites(sequence, minUnit, maxUnit, minRepeats)
            .Select(r => $"{r.Position},{r.RepeatUnit},{r.RepeatCount}");

        Assert.That(string.Join(";", results), Is.EqualTo(expected));
    }

    /// <summary>
    /// The cancellable DnaSequence overload validates parameters like the others (minUnitLength = 0
    /// previously reached the scan with an empty unit, which never terminates).
    /// </summary>
    [Test]
    public void FindMicrosatellites_CancellableDnaOverload_InvalidParameters_Throw()
    {
        var dna = new DnaSequence("ACACAC");
        Assert.Multiple(() =>
        {
            Assert.Throws<ArgumentOutOfRangeException>(() =>
                RepeatFinder.FindMicrosatellites(dna, 0, 6, 3, CancellationToken.None).ToList());
            Assert.Throws<ArgumentOutOfRangeException>(() =>
                RepeatFinder.FindMicrosatellites(dna, 3, 2, 3, CancellationToken.None).ToList());
            Assert.Throws<ArgumentOutOfRangeException>(() =>
                RepeatFinder.FindMicrosatellites(dna, 1, 6, 1, CancellationToken.None).ToList());
        });
    }

    /// <summary>
    /// S06: Adjacent different repeat types are detected independently.
    /// A poly-A tract followed immediately by CAG repeats should yield both.
    /// </summary>
    [Test]
    public void FindMicrosatellites_AdjacentDifferentRepeatTypes_DetectsBoth()
    {
        // A×5 immediately followed by CAG×3
        var sequence = new DnaSequence("AAAAACAGCAGCAG");

        var results = RepeatFinder.FindMicrosatellites(sequence, 1, 6, 3).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(results, Has.Count.EqualTo(2));
            Assert.That(results[0].RepeatUnit, Is.EqualTo("A"));
            Assert.That(results[0].RepeatCount, Is.EqualTo(5));
            Assert.That(results[0].Position, Is.EqualTo(0));
            Assert.That(results[0].RepeatType, Is.EqualTo(RepeatType.Mononucleotide));
            Assert.That(results[1].RepeatUnit, Is.EqualTo("CAG"));
            Assert.That(results[1].RepeatCount, Is.EqualTo(3));
            Assert.That(results[1].Position, Is.EqualTo(5));
            Assert.That(results[1].RepeatType, Is.EqualTo(RepeatType.Trinucleotide));
        });
    }

    #endregion

    #region API Overload Tests

    /// <summary>
    /// String overload should produce same results as DnaSequence overload.
    /// </summary>
    [Test]
    public void FindMicrosatellites_StringOverload_ProducesSameResults()
    {
        const string sequenceStr = "CAGCAGCAGCAG";
        var dnaSequence = new DnaSequence(sequenceStr);

        var stringResults = RepeatFinder.FindMicrosatellites(sequenceStr, 3, 3, 3).ToList();
        var dnaResults = RepeatFinder.FindMicrosatellites(dnaSequence, 3, 3, 3).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(stringResults, Has.Count.EqualTo(dnaResults.Count));
            for (int i = 0; i < stringResults.Count; i++)
            {
                Assert.That(stringResults[i].RepeatUnit, Is.EqualTo(dnaResults[i].RepeatUnit));
                Assert.That(stringResults[i].RepeatCount, Is.EqualTo(dnaResults[i].RepeatCount));
                Assert.That(stringResults[i].Position, Is.EqualTo(dnaResults[i].Position));
            }
        });
    }

    #endregion

    #region Error Handling

    /// <summary>
    /// Null DnaSequence should throw ArgumentNullException.
    /// </summary>
    [Test]
    public void FindMicrosatellites_NullDnaSequence_ThrowsArgumentNullException()
    {
        Assert.Throws<ArgumentNullException>(() =>
            RepeatFinder.FindMicrosatellites((DnaSequence)null!, 1, 6, 3).ToList());
    }

    /// <summary>
    /// minUnitLength < 1 should throw ArgumentOutOfRangeException.
    /// </summary>
    [Test]
    public void FindMicrosatellites_MinUnitLengthZero_ThrowsArgumentOutOfRangeException()
    {
        var sequence = new DnaSequence("ACGT");

        Assert.Throws<ArgumentOutOfRangeException>(() =>
            RepeatFinder.FindMicrosatellites(sequence, 0, 6, 3).ToList());
    }

    /// <summary>
    /// maxUnitLength < minUnitLength should throw ArgumentOutOfRangeException.
    /// </summary>
    [Test]
    public void FindMicrosatellites_MaxLessThanMin_ThrowsArgumentOutOfRangeException()
    {
        var sequence = new DnaSequence("ACGT");

        Assert.Throws<ArgumentOutOfRangeException>(() =>
            RepeatFinder.FindMicrosatellites(sequence, 5, 4, 3).ToList());
    }

    /// <summary>
    /// minRepeats < 2 should throw ArgumentOutOfRangeException.
    /// (A "repeat" requires at least 2 occurrences by definition)
    /// </summary>
    [Test]
    public void FindMicrosatellites_MinRepeatsOne_ThrowsArgumentOutOfRangeException()
    {
        var sequence = new DnaSequence("ACGT");

        Assert.Throws<ArgumentOutOfRangeException>(() =>
            RepeatFinder.FindMicrosatellites(sequence, 1, 6, 1).ToList());
    }

    /// <summary>
    /// String overload: minUnitLength < 1 should throw ArgumentOutOfRangeException.
    /// Ensures API parity with DnaSequence overload.
    /// </summary>
    [Test]
    public void FindMicrosatellites_StringOverload_MinUnitLengthZero_Throws()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() =>
            RepeatFinder.FindMicrosatellites("ACGT", 0, 6, 3).ToList());
    }

    /// <summary>
    /// String overload: maxUnitLength < minUnitLength should throw ArgumentOutOfRangeException.
    /// Ensures API parity with DnaSequence overload.
    /// </summary>
    [Test]
    public void FindMicrosatellites_StringOverload_MaxLessThanMin_Throws()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() =>
            RepeatFinder.FindMicrosatellites("ACGT", 5, 4, 3).ToList());
    }

    /// <summary>
    /// String overload: minRepeats < 2 should throw ArgumentOutOfRangeException.
    /// Ensures API parity with DnaSequence overload.
    /// </summary>
    [Test]
    public void FindMicrosatellites_StringOverload_MinRepeatsOne_Throws()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() =>
            RepeatFinder.FindMicrosatellites("ACGT", 1, 6, 1).ToList());
    }

    #endregion

    // NOTE: GetTandemRepeatSummary tests are in RepeatFinderTests.cs (canonical location)
    // Duplicate tests removed as part of REP-TANDEM-001 consolidation.

    #region Cancellation Smoke Test

    /// <summary>
    /// Cancellation overload should complete normally when not cancelled.
    /// (Deep cancellation testing is in PerformanceExtensionsTests)
    /// </summary>
    [Test]
    public void FindMicrosatellites_WithCancellationToken_CompletesNormally()
    {
        var sequence = new DnaSequence("ATATATATATAT");
        using var cts = new CancellationTokenSource();

        var results = RepeatFinder.FindMicrosatellites(sequence, 1, 6, 3, cts.Token).ToList();

        Assert.That(results, Is.Not.Empty);
    }

    #endregion
}
