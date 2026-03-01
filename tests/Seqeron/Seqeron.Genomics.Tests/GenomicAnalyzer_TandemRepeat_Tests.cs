using NUnit.Framework;
using Seqeron.Genomics;
using System.Diagnostics;
using System.Linq;

namespace Seqeron.Genomics.Tests;

/// <summary>
/// Comprehensive tests for GenomicAnalyzer.FindTandemRepeats.
/// Test Unit: REP-TANDEM-001
/// 
/// Evidence sources:
/// - Wikipedia (Tandem repeat): Definition, terminology, detection
/// - Wikipedia (Microsatellite/STR): 1–6 bp classification, forensic use, disease associations
/// - Richard et al. (2008): Comparative genomics of DNA repeats
/// </summary>
[TestFixture]
public class GenomicAnalyzer_TandemRepeat_Tests
{
    #region MUST Tests - Core Algorithm Verification

    /// <summary>
    /// M1: Simple trinucleotide repeat detection — ATG repeated 3 times.
    /// Evidence: Wikipedia (Tandem repeat): definition — adjacent repeating pattern.
    /// </summary>
    [Test]
    public void FindTandemRepeats_SimpleTrinucleotide_FindsRepeat()
    {
        var dna = new DnaSequence("ATGATGATG");

        var tandems = GenomicAnalyzer.FindTandemRepeats(dna, minUnitLength: 3, minRepetitions: 3).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(tandems, Has.Count.EqualTo(1), "Exactly one repeat at unitLen=3");
            Assert.That(tandems[0].Unit, Is.EqualTo("ATG"));
            Assert.That(tandems[0].Position, Is.EqualTo(0));
            Assert.That(tandems[0].Repetitions, Is.EqualTo(3));
        });
    }

    /// <summary>
    /// M2: Dinucleotide repeat detection — most common microsatellite type in human genome.
    /// Evidence: Wikipedia (Microsatellite) — 50,000–100,000 dinucleotide loci in human genome.
    /// </summary>
    [Test]
    public void FindTandemRepeats_DinucleotideRepeat_FindsRepeat()
    {
        var dna = new DnaSequence("CACACACA");

        var tandems = GenomicAnalyzer.FindTandemRepeats(dna, minUnitLength: 2, minRepetitions: 4).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(tandems, Has.Count.EqualTo(1), "Exactly one repeat at unitLen=2, minReps=4");
            Assert.That(tandems[0].Unit, Is.EqualTo("CA"));
            Assert.That(tandems[0].Position, Is.EqualTo(0));
            Assert.That(tandems[0].Repetitions, Is.EqualTo(4));
        });
    }

    /// <summary>
    /// M3: Mononucleotide repeat (homopolymer run) detection.
    /// Evidence: Wikipedia (Microsatellite) — 1–6 bp classification includes 1 bp units.
    /// </summary>
    [Test]
    public void FindTandemRepeats_MononucleotideRepeat_FindsRepeat()
    {
        var dna = new DnaSequence("AAAAA");

        var tandems = GenomicAnalyzer.FindTandemRepeats(dna, minUnitLength: 1, minRepetitions: 5).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(tandems, Has.Count.EqualTo(1), "Exactly one repeat at unitLen=1, minReps=5");
            Assert.That(tandems[0].Unit, Is.EqualTo("A"));
            Assert.That(tandems[0].Position, Is.EqualTo(0));
            Assert.That(tandems[0].Repetitions, Is.EqualTo(5));
        });
    }

    /// <summary>
    /// M4: Tetranucleotide repeat — forensic STR standard.
    /// Evidence: Wikipedia (Microsatellite) — forensic STRs are tetra-/pentanucleotide repeats.
    /// </summary>
    [Test]
    public void FindTandemRepeats_TetranucleotideRepeat_FindsRepeat()
    {
        // GATA repeats are used in forensic DNA profiling
        var dna = new DnaSequence("GATAGATAGATA");

        var tandems = GenomicAnalyzer.FindTandemRepeats(dna, minUnitLength: 4, minRepetitions: 3).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(tandems, Has.Count.EqualTo(1), "Exactly one repeat at unitLen=4, minReps=3");
            Assert.That(tandems[0].Unit, Is.EqualTo("GATA"));
            Assert.That(tandems[0].Position, Is.EqualTo(0));
            Assert.That(tandems[0].Repetitions, Is.EqualTo(3));
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
    /// M7: MinRepetitions filter is respected — only repeats meeting threshold are returned.
    /// Evidence: Algorithm specification.
    /// </summary>
    [Test]
    public void FindTandemRepeats_MinRepetitionsFilter_RespectsThreshold()
    {
        var dna = new DnaSequence("ATATAT"); // AT × 3

        var resultsMin2 = GenomicAnalyzer.FindTandemRepeats(dna, minUnitLength: 2, minRepetitions: 2).ToList();
        var resultsMin3 = GenomicAnalyzer.FindTandemRepeats(dna, minUnitLength: 2, minRepetitions: 3).ToList();
        var resultsMin4 = GenomicAnalyzer.FindTandemRepeats(dna, minUnitLength: 2, minRepetitions: 4).ToList();

        Assert.Multiple(() =>
        {
            // minReps=2: AT detected with exactly 3 repetitions
            var atMin2 = resultsMin2.First(t => t.Unit == "AT");
            Assert.That(atMin2.Repetitions, Is.EqualTo(3), "AT × 3 with minReps=2");

            // minReps=3: still found
            var atMin3 = resultsMin3.First(t => t.Unit == "AT");
            Assert.That(atMin3.Repetitions, Is.EqualTo(3), "AT × 3 with minReps=3");

            // minReps=4: not found
            Assert.That(resultsMin4, Is.Empty, "minReps=4 should not find AT × 3");
        });
    }

    /// <summary>
    /// M8: MinUnitLength filter is respected — units below threshold are excluded.
    /// Evidence: Algorithm specification.
    /// </summary>
    [Test]
    public void FindTandemRepeats_MinUnitLengthFilter_RespectsThreshold()
    {
        var dna = new DnaSequence("AAATTT"); // A × 3, T × 3

        var resultsMin1 = GenomicAnalyzer.FindTandemRepeats(dna, minUnitLength: 1, minRepetitions: 3).ToList();
        var resultsMin2 = GenomicAnalyzer.FindTandemRepeats(dna, minUnitLength: 2, minRepetitions: 3).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(resultsMin1, Has.Count.EqualTo(2), "minUnit=1 should find A × 3 and T × 3");
            Assert.That(resultsMin1[0].Unit, Is.EqualTo("A"));
            Assert.That(resultsMin1[1].Unit, Is.EqualTo("T"));
            Assert.That(resultsMin2, Is.Empty, "minUnit=2 should find no di+ repeats in AAATTT");
        });
    }

    /// <summary>
    /// M9: Position is accurate and 0-based.
    /// Evidence: Implementation contract.
    /// </summary>
    [Test]
    public void FindTandemRepeats_PositionCorrect_ZeroBased()
    {
        var dna = new DnaSequence("CCATGATGATGCC");

        var tandems = GenomicAnalyzer.FindTandemRepeats(dna, minUnitLength: 3, minRepetitions: 3).ToList();
        var atgRepeat = tandems.First(t => t.Unit == "ATG");

        Assert.Multiple(() =>
        {
            Assert.That(atgRepeat.Unit, Is.EqualTo("ATG"));
            Assert.That(atgRepeat.Position, Is.EqualTo(2),
                "ATG repeat should start at position 2 (0-based)");
        });
    }

    /// <summary>
    /// M10: Repetition count is accurate for CAG × 5 (Huntington's-relevant trinucleotide).
    /// Evidence: Wikipedia (Microsatellite) — trinucleotide repeat disorders;
    ///           Huntington's disease caused by CAG expansions (>36 repeats pathogenic).
    /// </summary>
    [Test]
    public void FindTandemRepeats_RepetitionCount_Accurate()
    {
        var dna = new DnaSequence("CAGCAGCAGCAGCAG"); // CAG × 5

        var tandems = GenomicAnalyzer.FindTandemRepeats(dna, minUnitLength: 3, minRepetitions: 2).ToList();
        var cagRepeat = tandems.First(t => t.Unit == "CAG");

        Assert.Multiple(() =>
        {
            Assert.That(cagRepeat.Unit, Is.EqualTo("CAG"));
            Assert.That(cagRepeat.Position, Is.EqualTo(0));
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
        var dna = new DnaSequence("GATAGATAGATA"); // GATA × 3

        var tandems = GenomicAnalyzer.FindTandemRepeats(dna, minUnitLength: 4, minRepetitions: 3).ToList();
        var gataRepeat = tandems.First(t => t.Unit == "GATA");

        Assert.Multiple(() =>
        {
            Assert.That(gataRepeat.TotalLength, Is.EqualTo(gataRepeat.Unit.Length * gataRepeat.Repetitions),
                "TotalLength must equal Unit.Length × Repetitions");
            Assert.That(gataRepeat.TotalLength, Is.EqualTo(12),
                "GATA × 3 = 12 bases");
        });
    }

    /// <summary>
    /// M12: FullSequence property returns correct reconstruction.
    /// Evidence: Documented property behavior.
    /// </summary>
    [Test]
    public void FindTandemRepeats_FullSequence_Reconstructable()
    {
        var dna = new DnaSequence("ATGATGATG"); // ATG × 3

        var tandems = GenomicAnalyzer.FindTandemRepeats(dna, minUnitLength: 3, minRepetitions: 3).ToList();
        var atgRepeat = tandems.First(t => t.Unit == "ATG");

        Assert.Multiple(() =>
        {
            Assert.That(atgRepeat.FullSequence, Is.EqualTo("ATGATGATG"),
                "FullSequence should match actual sequence");
            Assert.That(atgRepeat.FullSequence.Length, Is.EqualTo(atgRepeat.TotalLength),
                "FullSequence length must match TotalLength");
        });
    }

    /// <summary>
    /// M13: Pentanucleotide repeat detection — forensic STR standard.
    /// Evidence: Wikipedia (Microsatellite) — forensic STRs are tetra-/pentanucleotide repeats.
    /// Note: HUMTH01 locus uses TCAT; D18S51 uses AGAA; Penta E uses AAAGA.
    /// </summary>
    [Test]
    public void FindTandemRepeats_PentanucleotideRepeat_ForensicSTR()
    {
        // Pentanucleotide repeat — Penta E forensic locus uses AAAGA
        var dna = new DnaSequence("AAAGAAAAGAAAAGA"); // AAAGA × 3

        var tandems = GenomicAnalyzer.FindTandemRepeats(dna, minUnitLength: 5, minRepetitions: 3).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(tandems, Has.Count.EqualTo(1), "Exactly one pentanucleotide repeat");
            Assert.That(tandems[0].Unit, Is.EqualTo("AAAGA"));
            Assert.That(tandems[0].Position, Is.EqualTo(0));
            Assert.That(tandems[0].Repetitions, Is.EqualTo(3));
            Assert.That(tandems[0].TotalLength, Is.EqualTo(15));
        });
    }

    #endregion

    #region SHOULD Tests - Important Edge Cases

    /// <summary>
    /// S1: Long repeat sequence — verify exact count for AT × 20.
    /// Evidence: Robustness testing.
    /// </summary>
    [Test]
    public void FindTandemRepeats_LongRepeat_HandlesCorrectly()
    {
        string sequence = string.Concat(Enumerable.Repeat("AT", 20));
        var dna = new DnaSequence(sequence);

        var tandems = GenomicAnalyzer.FindTandemRepeats(dna, minUnitLength: 2, minRepetitions: 20).ToList();
        var atRepeat = tandems.First(t => t.Unit == "AT");

        Assert.Multiple(() =>
        {
            Assert.That(atRepeat.Position, Is.EqualTo(0));
            Assert.That(atRepeat.Repetitions, Is.EqualTo(20), "Should count exactly 20 AT repeats");
            Assert.That(atRepeat.TotalLength, Is.EqualTo(40));
        });
    }

    /// <summary>
    /// S2: Entire sequence is one repeat — covers full length.
    /// Evidence: Edge case where whole sequence is repetitive.
    /// </summary>
    [Test]
    public void FindTandemRepeats_EntireSequenceIsRepeat_CoversFullLength()
    {
        var dna = new DnaSequence("CGTCGTCGTCGT"); // CGT × 4, entire sequence

        var tandems = GenomicAnalyzer.FindTandemRepeats(dna, minUnitLength: 3, minRepetitions: 4).ToList();
        var cgtRepeat = tandems.First(t => t.Unit == "CGT");

        Assert.Multiple(() =>
        {
            Assert.That(cgtRepeat.Position, Is.EqualTo(0));
            Assert.That(cgtRepeat.Repetitions, Is.EqualTo(4));
            Assert.That(cgtRepeat.TotalLength, Is.EqualTo(12),
                "Repeat should cover entire 12-base sequence");
        });
    }

    /// <summary>
    /// S3: Adjacent different mononucleotide repeats both detected with exact counts.
    /// Evidence: Common biological scenario.
    /// </summary>
    [Test]
    public void FindTandemRepeats_AdjacentDifferentRepeats_FindsBoth()
    {
        var dna = new DnaSequence("AAAAAATTTTT"); // A × 6, T × 5

        var tandems = GenomicAnalyzer.FindTandemRepeats(dna, minUnitLength: 1, minRepetitions: 5).ToList();

        Assert.Multiple(() =>
        {
            var aRepeat = tandems.First(t => t.Unit == "A");
            Assert.That(aRepeat.Position, Is.EqualTo(0));
            Assert.That(aRepeat.Repetitions, Is.EqualTo(6), "A × 6");

            var tRepeat = tandems.First(t => t.Unit == "T");
            Assert.That(tRepeat.Position, Is.EqualTo(6));
            Assert.That(tRepeat.Repetitions, Is.EqualTo(5), "T × 5");
        });
    }

    /// <summary>
    /// S4: Hexanucleotide repeat — upper boundary of microsatellite range.
    /// Evidence: Wikipedia (Microsatellite) — 1–6 bp classification.
    /// </summary>
    [Test]
    public void FindTandemRepeats_HexanucleotideRepeat_Boundary()
    {
        // ACGTAC × 3 = 18 bases
        var dna = new DnaSequence("ACGTACACGTACACGTAC");

        var tandems = GenomicAnalyzer.FindTandemRepeats(dna, minUnitLength: 6, minRepetitions: 3).ToList();
        var hexRepeat = tandems.First(t => t.Unit.Length == 6);

        Assert.Multiple(() =>
        {
            Assert.That(hexRepeat.Unit, Is.EqualTo("ACGTAC"));
            Assert.That(hexRepeat.Position, Is.EqualTo(0));
            Assert.That(hexRepeat.Repetitions, Is.EqualTo(3));
            Assert.That(hexRepeat.TotalLength, Is.EqualTo(18));
        });
    }

    /// <summary>
    /// S5: Case sensitivity — DnaSequence normalizes to uppercase.
    /// Evidence: DnaSequence normalizes input to uppercase.
    /// </summary>
    [Test]
    public void FindTandemRepeats_CaseSensitivity_UpperCase()
    {
        var dna = new DnaSequence("atgatgatg");

        var tandems = GenomicAnalyzer.FindTandemRepeats(dna, minUnitLength: 3, minRepetitions: 3).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(tandems, Has.Count.EqualTo(1));
            Assert.That(tandems[0].Unit, Is.EqualTo("ATG"),
                "Unit should be uppercase after normalization");
            Assert.That(tandems[0].Repetitions, Is.EqualTo(3));
        });
    }

    #endregion

    #region COULD Tests

    /// <summary>
    /// C1: Performance baseline — O(n²) brute-force should complete for a 2 kb sequence.
    /// Evidence: Algorithm performance contract.
    /// </summary>
    [Test]
    public void FindTandemRepeats_PerformanceBaseline_MediumSequence()
    {
        // 2,000 bp sequence: ACGT repeated 500 times
        string sequence = string.Concat(Enumerable.Repeat("ACGT", 500));
        var dna = new DnaSequence(sequence);

        var sw = Stopwatch.StartNew();
        var tandems = GenomicAnalyzer.FindTandemRepeats(dna, minUnitLength: 2, minRepetitions: 3).ToList();
        sw.Stop();

        Assert.That(sw.ElapsedMilliseconds, Is.LessThan(30_000),
            "2 kb sequence should complete within 30 seconds (debug-safe)");
    }

    /// <summary>
    /// C2: TTAGGG telomere repeat — hexanucleotide repeat in vertebrate telomeres.
    /// Evidence: Wikipedia (Microsatellite) — vertebrate telomeres have TTAGGG repeat motif.
    /// </summary>
    [Test]
    public void FindTandemRepeats_TelomereRepeat_TTAGGG()
    {
        // Vertebrate telomere repeat: TTAGGG × 4
        var dna = new DnaSequence("TTAGGGTTAGGGTTAGGGTTAGGG");

        var tandems = GenomicAnalyzer.FindTandemRepeats(dna, minUnitLength: 6, minRepetitions: 4).ToList();
        var telomere = tandems.First(t => t.Unit == "TTAGGG");

        Assert.Multiple(() =>
        {
            Assert.That(telomere.Position, Is.EqualTo(0));
            Assert.That(telomere.Repetitions, Is.EqualTo(4));
            Assert.That(telomere.TotalLength, Is.EqualTo(24));
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
