using NUnit.Framework;
using Seqeron.Genomics;
using System.Linq;
using System.Threading;

namespace Seqeron.Genomics.Tests;

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
            Assert.That(results, Has.Count.GreaterThanOrEqualTo(1));
            var caRepeat = results.First(r => r.RepeatUnit == "CA" || r.RepeatUnit == "AC");
            Assert.That(caRepeat.RepeatCount, Is.GreaterThanOrEqualTo(3));
            Assert.That(caRepeat.RepeatType, Is.EqualTo(RepeatType.Dinucleotide));
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
            Assert.That(results, Has.Count.GreaterThanOrEqualTo(1));
            var cagRepeat = results.First(r => r.RepeatUnit == "CAG");
            Assert.That(cagRepeat.RepeatCount, Is.EqualTo(5));
            Assert.That(cagRepeat.RepeatType, Is.EqualTo(RepeatType.Trinucleotide));
            Assert.That(cagRepeat.Position, Is.EqualTo(3));
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
            Assert.That(results, Has.Count.GreaterThanOrEqualTo(1));
            // GATA family includes rotations: GATA, ATAG, TAGA, AGAT
            var gataRepeat = results.First(r =>
                r.RepeatUnit == "GATA" || r.RepeatUnit == "ATAG" ||
                r.RepeatUnit == "TAGA" || r.RepeatUnit == "AGAT");
            Assert.That(gataRepeat.RepeatType, Is.EqualTo(RepeatType.Tetranucleotide));
            Assert.That(gataRepeat.RepeatCount, Is.GreaterThanOrEqualTo(3));
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

        Assert.That(results, Has.Count.GreaterThanOrEqualTo(1), $"Should find at least one repeat of {unit}");
        Assert.That(results.Any(r => r.RepeatType == expectedType),
            Is.True, $"Should have a {expectedType} repeat");
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
            Assert.That(results, Has.Count.GreaterThanOrEqualTo(2),
                "Should find at least mononucleotide A and other repeats");
            Assert.That(results.Any(r => r.RepeatUnit == "A" && r.RepeatCount >= 5),
                Is.True, "Should find poly-A tract");
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
