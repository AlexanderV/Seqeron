namespace Seqeron.Genomics.Tests.Snapshots;

/// <summary>
/// Snapshot tests for RNA secondary structure prediction.
/// Verifies structure prediction, stem-loop finding, and energy calculations.
///
/// Test Units: RNA-STRUCT-001, RNA-STEMLOOP-001, RNA-ENERGY-001
/// </summary>
[TestFixture]
[Category("Snapshot")]
[Category("Analysis")]
public class RnaSnapshotTests
{
    // A sequence designed to form a stem-loop: GGGAAACCC pairs G-C with AAA loop
    private const string TestRna = "GGGAAACCC";

    [Test]
    public Task PredictStructure_KnownHairpin_MatchesSnapshot()
    {
        var structure = RnaSecondaryStructure.PredictStructure(TestRna);
        return Verify(new
        {
            structure.Sequence,
            structure.DotBracket,
            BasePairCount = structure.BasePairs.Count,
            StemLoopCount = structure.StemLoops.Count,
            MFE = Math.Round(structure.MinimumFreeEnergy, 4)
        });
    }

    [Test]
    public Task FindStemLoops_KnownSequence_MatchesSnapshot()
    {
        string rna = "GGGGAAAACCCCUUUUUUGGGAAACCC";
        var stemLoops = RnaSecondaryStructure.FindStemLoops(rna, minStemLength: 3)
            .Select(sl => new
            {
                sl.Start,
                sl.End,
                StemLength = sl.Stem.Length,
                LoopSize = sl.Loop.Size,
                Energy = Math.Round(sl.TotalFreeEnergy, 4)
            })
            .ToList();

        return Verify(new { StemLoopCount = stemLoops.Count, StemLoops = stemLoops });
    }

    [Test]
    public Task MinimumFreeEnergy_KnownSequences_MatchesSnapshot()
    {
        var sequences = new Dictionary<string, string>
        {
            ["Simple_Hairpin"] = "GGGAAACCC",
            ["Longer_Hairpin"] = "GGGGAAAACCCC",
            ["PolyA_NoStructure"] = "AAAAAAAAAA"
        };

        var energies = sequences.ToDictionary(
            kv => kv.Key,
            kv => Math.Round(RnaSecondaryStructure.CalculateMinimumFreeEnergy(kv.Value), 4));

        return Verify(new { Energies = energies });
    }
}
