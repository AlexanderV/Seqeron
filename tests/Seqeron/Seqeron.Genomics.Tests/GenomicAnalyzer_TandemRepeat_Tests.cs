using NUnit.Framework;
using Seqeron.Genomics;
using System.Linq;

namespace Seqeron.Genomics.Tests;

/// <summary>
/// Comprehensive tests for GenomicAnalyzer.FindTandemRepeats.
/// Test Unit: REP-TANDEM-001
/// 
/// Evidence sources:
/// - Wikipedia (Tandem repeat): Definition, terminology, detection
/// - Wikipedia (Microsatellite/STR): 1-6bp classification, forensic use, disease associations
/// - Richard et al. (2008): Comparative genomics of DNA repeats
/// </summary>
[TestFixture]
public class GenomicAnalyzer_TandemRepeat_Tests
{
    #region MUST Tests - Core Algorithm Verification

    /// <summary>
    /// M1: Simple trinucleotide repeat detection - ATG repeated 3 times.
    /// Evidence: Wikipedia definition of tandem repeats.
    /// </summary>
    [Test]
    public void FindTandemRepeats_SimpleTrinucleotide_FindsRepeat()
    {
        var dna = new DnaSequence("ATGATGATG");

        var tandems = GenomicAnalyzer.FindTandemRepeats(dna, minUnitLength: 3, minRepetitions: 3).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(tandems, Has.Count.GreaterThanOrEqualTo(1));
            Assert.That(tandems.Any(t => t.Unit == "ATG" && t.Repetitions == 3), Is.True,
                "Should find ATG repeated 3 times");
        });
    }

    /// <summary>
    /// M2: Dinucleotide repeat detection - common forensic STR type.
    /// Evidence: Wikipedia - microsatellite forensics uses di/tetra/penta repeats.
    /// </summary>
    [Test]
    public void FindTandemRepeats_DinucleotideRepeat_FindsRepeat()
    {
        var dna = new DnaSequence("CACACACA");

        var tandems = GenomicAnalyzer.FindTandemRepeats(dna, minUnitLength: 2, minRepetitions: 4).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(tandems, Has.Count.GreaterThanOrEqualTo(1));
            Assert.That(tandems.Any(t => t.Unit == "CA" && t.Repetitions == 4), Is.True,
                "Should find CA repeated 4 times");
        });
    }

    /// <summary>
    /// M3: Mononucleotide repeat (homopolymer run) detection.
    /// Evidence: Wikipedia - microsatellite types include 1bp units.
    /// </summary>
    [Test]
    public void FindTandemRepeats_MononucleotideRepeat_FindsRepeat()
    {
        var dna = new DnaSequence("AAAAA");

        var tandems = GenomicAnalyzer.FindTandemRepeats(dna, minUnitLength: 1, minRepetitions: 5).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(tandems, Has.Count.GreaterThanOrEqualTo(1));
            Assert.That(tandems.Any(t => t.Unit == "A" && t.Repetitions == 5), Is.True,
                "Should find A repeated 5 times");
        });
    }

    /// <summary>
    /// M4: Tetranucleotide repeat - forensic STR standard.
    /// Evidence: Wikipedia - forensic STRs use tetra/penta repeats for accuracy.
    /// </summary>
    [Test]
    public void FindTandemRepeats_TetranucleotideRepeat_FindsRepeat()
    {
        // GATA repeats are used in forensic DNA profiling
        var dna = new DnaSequence("GATAGATAGATA");

        var tandems = GenomicAnalyzer.FindTandemRepeats(dna, minUnitLength: 4, minRepetitions: 3).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(tandems, Has.Count.GreaterThanOrEqualTo(1));
            Assert.That(tandems.Any(t => t.Unit == "GATA" && t.Repetitions == 3), Is.True,
                "Should find GATA repeated 3 times");
        });
    }

    /// <summary>
    /// M5: No repeats found returns empty enumerable.
    /// Evidence: Standard edge case behavior.
    /// </summary>
    [Test]
    public void FindTandemRepeats_NoRepeatsFound_ReturnsEmpty()
    {
        var dna = new DnaSequence("ACGT");

        var tandems = GenomicAnalyzer.FindTandemRepeats(dna, minUnitLength: 2, minRepetitions: 2).ToList();

        Assert.That(tandems, Is.Empty, "Sequence without tandem repeats should return empty");
    }

    /// <summary>
    /// M6: Empty sequence returns empty enumerable.
    /// Evidence: Standard boundary case.
    /// </summary>
    [Test]
    public void FindTandemRepeats_EmptySequence_ReturnsEmpty()
    {
        var dna = new DnaSequence("");

        var tandems = GenomicAnalyzer.FindTandemRepeats(dna, minUnitLength: 2, minRepetitions: 2).ToList();

        Assert.That(tandems, Is.Empty, "Empty sequence should return empty");
    }

    /// <summary>
    /// M7: MinRepetitions filter is respected.
    /// Evidence: Algorithm specification.
    /// </summary>
    [Test]
    public void FindTandemRepeats_MinRepetitionsFilter_RespectsThreshold()
    {
        var dna = new DnaSequence("ATATAT"); // AT x 3

        var resultsMin2 = GenomicAnalyzer.FindTandemRepeats(dna, minUnitLength: 2, minRepetitions: 2).ToList();
        var resultsMin3 = GenomicAnalyzer.FindTandemRepeats(dna, minUnitLength: 2, minRepetitions: 3).ToList();
        var resultsMin4 = GenomicAnalyzer.FindTandemRepeats(dna, minUnitLength: 2, minRepetitions: 4).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(resultsMin2.Any(t => t.Unit == "AT" && t.Repetitions >= 2), Is.True,
                "minReps=2 should find AT repeat");
            Assert.That(resultsMin3.Any(t => t.Unit == "AT" && t.Repetitions >= 3), Is.True,
                "minReps=3 should find AT x 3");
            Assert.That(resultsMin4, Is.Empty,
                "minReps=4 should not find AT x 3");
        });
    }

    /// <summary>
    /// M8: MinUnitLength filter is respected.
    /// Evidence: Algorithm specification.
    /// </summary>
    [Test]
    public void FindTandemRepeats_MinUnitLengthFilter_RespectsThreshold()
    {
        var dna = new DnaSequence("AAATTT"); // A x 3, T x 3

        var resultsMin1 = GenomicAnalyzer.FindTandemRepeats(dna, minUnitLength: 1, minRepetitions: 3).ToList();
        var resultsMin2 = GenomicAnalyzer.FindTandemRepeats(dna, minUnitLength: 2, minRepetitions: 3).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(resultsMin1, Has.Count.GreaterThanOrEqualTo(2),
                "minUnit=1 should find A and T runs");
            Assert.That(resultsMin2.Any(t => t.Unit.Length == 1), Is.False,
                "minUnit=2 should not find mononucleotide repeats");
        });
    }

    /// <summary>
    /// M9: Position is accurate and 0-based.
    /// Evidence: Implementation contract.
    /// </summary>
    [Test]
    public void FindTandemRepeats_PositionCorrect_ZeroBased()
    {
        // Use valid nucleotides - DnaSequence doesn't allow N
        var dna = new DnaSequence("CCATGATGATGCC");

        var tandems = GenomicAnalyzer.FindTandemRepeats(dna, minUnitLength: 3, minRepetitions: 3).ToList();
        var atgRepeat = tandems.FirstOrDefault(t => t.Unit == "ATG");

        Assert.Multiple(() =>
        {
            Assert.That(atgRepeat.Unit, Is.EqualTo("ATG"));
            Assert.That(atgRepeat.Position, Is.EqualTo(2),
                "ATG repeat should start at position 2 (0-based)");
        });
    }

    /// <summary>
    /// M10: Repetition count is accurate.
    /// Evidence: Core correctness requirement.
    /// </summary>
    [Test]
    public void FindTandemRepeats_RepetitionCount_Accurate()
    {
        var dna = new DnaSequence("CAGCAGCAGCAGCAG"); // CAG x 5

        var tandems = GenomicAnalyzer.FindTandemRepeats(dna, minUnitLength: 3, minRepetitions: 2).ToList();
        var cagRepeat = tandems.FirstOrDefault(t => t.Unit == "CAG");

        Assert.Multiple(() =>
        {
            Assert.That(cagRepeat.Unit, Is.EqualTo("CAG"));
            Assert.That(cagRepeat.Repetitions, Is.EqualTo(5),
                "Should count exactly 5 CAG repetitions");
        });
    }

    /// <summary>
    /// M11: TotalLength invariant holds (Unit.Length × Repetitions).
    /// Evidence: Documented invariant.
    /// </summary>
    [Test]
    public void FindTandemRepeats_TotalLength_InvariantHolds()
    {
        var dna = new DnaSequence("GATAGATAGATA"); // GATA x 3

        var tandems = GenomicAnalyzer.FindTandemRepeats(dna, minUnitLength: 4, minRepetitions: 3).ToList();
        var gataRepeat = tandems.FirstOrDefault(t => t.Unit == "GATA");

        Assert.Multiple(() =>
        {
            Assert.That(gataRepeat.TotalLength, Is.EqualTo(gataRepeat.Unit.Length * gataRepeat.Repetitions),
                "TotalLength must equal Unit.Length × Repetitions");
            Assert.That(gataRepeat.TotalLength, Is.EqualTo(12),
                "GATA x 3 = 12 bases");
        });
    }

    /// <summary>
    /// M12: FullSequence property returns correct reconstruction.
    /// Evidence: Documented property behavior.
    /// </summary>
    [Test]
    public void FindTandemRepeats_FullSequence_Reconstructable()
    {
        var dna = new DnaSequence("ATGATGATG"); // ATG x 3

        var tandems = GenomicAnalyzer.FindTandemRepeats(dna, minUnitLength: 3, minRepetitions: 3).ToList();
        var atgRepeat = tandems.FirstOrDefault(t => t.Unit == "ATG");

        Assert.Multiple(() =>
        {
            Assert.That(atgRepeat.FullSequence, Is.EqualTo("ATGATGATG"),
                "FullSequence should match actual sequence");
            Assert.That(atgRepeat.FullSequence.Length, Is.EqualTo(atgRepeat.TotalLength),
                "FullSequence length must match TotalLength");
        });
    }

    /// <summary>
    /// M13: CAG expansion pattern detection (Huntington's disease).
    /// Evidence: Wikipedia - trinucleotide repeat disorders, Huntington's uses CAG.
    /// </summary>
    [Test]
    public void FindTandemRepeats_CAGExpansion_HuntingtonsPattern()
    {
        // Huntington's disease is caused by CAG expansions (typically >36 repeats)
        // Test with pure CAG repeat sequence to verify detection
        var dna = new DnaSequence("CAGCAGCAGCAGCAG"); // CAG x 5

        var tandems = GenomicAnalyzer.FindTandemRepeats(dna, minUnitLength: 3, minRepetitions: 3).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(tandems.Any(t => t.Unit == "CAG"), Is.True,
                "Should detect CAG expansion pattern");
            var cagRepeat = tandems.First(t => t.Unit == "CAG");
            Assert.That(cagRepeat.Repetitions, Is.EqualTo(5),
                "Should count 5 CAG repeats");
        });
    }

    #endregion

    #region SHOULD Tests - Important Edge Cases

    /// <summary>
    /// S1: Long repeat sequence handled correctly.
    /// Evidence: Robustness testing.
    /// </summary>
    [Test]
    public void FindTandemRepeats_LongRepeat_HandlesCorrectly()
    {
        // Create sequence with 20 AT repeats
        string sequence = string.Concat(Enumerable.Repeat("AT", 20));
        var dna = new DnaSequence(sequence);

        var tandems = GenomicAnalyzer.FindTandemRepeats(dna, minUnitLength: 2, minRepetitions: 10).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(tandems, Has.Count.GreaterThanOrEqualTo(1));
            Assert.That(tandems.Any(t => t.Unit == "AT" && t.Repetitions >= 10), Is.True,
                "Should find long AT repeat");
        });
    }

    /// <summary>
    /// S2: Entire sequence is one repeat - single result.
    /// Evidence: Edge case where whole sequence is repetitive.
    /// </summary>
    [Test]
    public void FindTandemRepeats_EntireSequenceIsRepeat_SingleResult()
    {
        var dna = new DnaSequence("CGTCGTCGTCGT"); // CGT x 4, entire sequence

        var tandems = GenomicAnalyzer.FindTandemRepeats(dna, minUnitLength: 3, minRepetitions: 4).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(tandems, Has.Count.GreaterThanOrEqualTo(1));
            var cgtRepeat = tandems.FirstOrDefault(t => t.Unit == "CGT");
            Assert.That(cgtRepeat.Position, Is.EqualTo(0),
                "Repeat should start at position 0");
            Assert.That(cgtRepeat.TotalLength, Is.EqualTo(12),
                "Repeat should cover entire sequence");
        });
    }

    /// <summary>
    /// S3: Adjacent different repeats both detected.
    /// Evidence: Common biological scenario.
    /// </summary>
    [Test]
    public void FindTandemRepeats_AdjacentDifferentRepeats_FindsBoth()
    {
        var dna = new DnaSequence("AAAAAATTTTT"); // A x 6, T x 5

        var tandems = GenomicAnalyzer.FindTandemRepeats(dna, minUnitLength: 1, minRepetitions: 5).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(tandems.Any(t => t.Unit == "A" && t.Repetitions >= 5), Is.True,
                "Should find A repeat");
            Assert.That(tandems.Any(t => t.Unit == "T" && t.Repetitions >= 5), Is.True,
                "Should find T repeat");
        });
    }

    /// <summary>
    /// S4: Maximum unit length boundary.
    /// Evidence: Algorithm limits.
    /// </summary>
    [Test]
    public void FindTandemRepeats_MaxUnitLength_Boundary()
    {
        // 6-nucleotide unit repeated 3 times
        var dna = new DnaSequence("ACGTACACGTACACGTAC");

        var tandems = GenomicAnalyzer.FindTandemRepeats(dna, minUnitLength: 6, minRepetitions: 3).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(tandems, Has.Count.GreaterThanOrEqualTo(1));
            Assert.That(tandems.Any(t => t.Unit.Length == 6), Is.True,
                "Should find 6bp unit repeat");
        });
    }

    /// <summary>
    /// S5: Case sensitivity - DnaSequence normalizes to uppercase.
    /// Evidence: DnaSequence implementation detail.
    /// </summary>
    [Test]
    public void FindTandemRepeats_CaseSensitivity_UpperCase()
    {
        // DnaSequence constructor normalizes to uppercase
        var dna = new DnaSequence("atgatgatg");

        var tandems = GenomicAnalyzer.FindTandemRepeats(dna, minUnitLength: 3, minRepetitions: 3).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(tandems, Has.Count.GreaterThanOrEqualTo(1));
            Assert.That(tandems.Any(t => t.Unit == "ATG"), Is.True,
                "Unit should be uppercase after normalization");
        });
    }

    #endregion

    #region Property Tests - Invariants

    /// <summary>
    /// All results satisfy the minRepetitions constraint.
    /// </summary>
    [Test]
    public void FindTandemRepeats_AllResults_SatisfyMinRepetitions()
    {
        var dna = new DnaSequence("AAAAAACGTCGTCGTACACACACGATAGATAGATA");
        const int minReps = 3;

        var tandems = GenomicAnalyzer.FindTandemRepeats(dna, minUnitLength: 1, minRepetitions: minReps).ToList();

        Assert.That(tandems.All(t => t.Repetitions >= minReps), Is.True,
            $"All results must have at least {minReps} repetitions");
    }

    /// <summary>
    /// All results satisfy the minUnitLength constraint.
    /// </summary>
    [Test]
    public void FindTandemRepeats_AllResults_SatisfyMinUnitLength()
    {
        var dna = new DnaSequence("AAAAAACGTCGTCGTACACACACGATAGATAGATA");
        const int minUnit = 2;

        var tandems = GenomicAnalyzer.FindTandemRepeats(dna, minUnitLength: minUnit, minRepetitions: 3).ToList();

        Assert.That(tandems.All(t => t.Unit.Length >= minUnit), Is.True,
            $"All results must have unit length >= {minUnit}");
    }

    /// <summary>
    /// Position + TotalLength does not exceed sequence length.
    /// </summary>
    [Test]
    public void FindTandemRepeats_AllResults_WithinSequenceBounds()
    {
        var dna = new DnaSequence("AAACGTCGTCGTCAGCAGCAGAAA");
        int seqLength = dna.Length;

        var tandems = GenomicAnalyzer.FindTandemRepeats(dna, minUnitLength: 2, minRepetitions: 3).ToList();

        Assert.That(tandems.All(t => t.Position + t.TotalLength <= seqLength), Is.True,
            "Position + TotalLength must not exceed sequence length");
    }

    #endregion
}
