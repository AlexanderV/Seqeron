namespace Seqeron.Genomics.Tests.Snapshots;

/// <summary>
/// Snapshot tests for codon optimization and usage analysis.
/// Verifies optimization results, CAI calculation, and rare codon detection.
///
/// Test Units: CODON-OPT-001, CODON-CAI-001, CODON-USAGE-001, CODON-RARE-001
/// </summary>
[TestFixture]
[Category("Snapshot")]
[Category("MolTools")]
public class CodonSnapshotTests
{
    // A synthetic coding sequence (multiples of 3)
    private const string TestCodingSequence = "ATGAAAGCGTTCAAGCGTACTGCGTGA";

    [Test]
    public Task OptimizeSequence_EColiTarget_MatchesSnapshot()
    {
        var result = CodonOptimizer.OptimizeSequence(TestCodingSequence, CodonOptimizer.EColiK12);
        return Verify(new
        {
            result.OriginalSequence,
            result.OptimizedSequence,
            result.ProteinSequence,
            OriginalCAI = Math.Round(result.OriginalCAI, 4),
            OptimizedCAI = Math.Round(result.OptimizedCAI, 4),
            GcOriginal = Math.Round(result.GcContentOriginal, 4),
            GcOptimized = Math.Round(result.GcContentOptimized, 4),
            result.ChangedCodons
        });
    }

    [Test]
    public Task CalculateCAI_KnownSequence_MatchesSnapshot()
    {
        double caiEcoli = CodonOptimizer.CalculateCAI(TestCodingSequence, CodonOptimizer.EColiK12);
        double caiHuman = CodonOptimizer.CalculateCAI(TestCodingSequence, CodonOptimizer.Human);
        double caiYeast = CodonOptimizer.CalculateCAI(TestCodingSequence, CodonOptimizer.Yeast);

        return Verify(new
        {
            CAI_EColi = Math.Round(caiEcoli, 4),
            CAI_Human = Math.Round(caiHuman, 4),
            CAI_Yeast = Math.Round(caiYeast, 4)
        });
    }

    [Test]
    public Task CodonUsage_KnownSequence_MatchesSnapshot()
    {
        var usage = CodonOptimizer.CalculateCodonUsage(TestCodingSequence);
        var sorted = usage.OrderBy(kv => kv.Key)
            .ToDictionary(kv => kv.Key, kv => kv.Value);
        return Verify(new { CodonCounts = sorted });
    }

    [Test]
    public Task FindRareCodons_KnownSequence_MatchesSnapshot()
    {
        var rareCodons = CodonOptimizer.FindRareCodons(TestCodingSequence, CodonOptimizer.EColiK12)
            .Select(r => new { r.Position, r.Codon, r.AminoAcid, Frequency = Math.Round(r.Frequency, 4) })
            .ToList();
        return Verify(new { RareCodonCount = rareCodons.Count, RareCodons = rareCodons });
    }
}
