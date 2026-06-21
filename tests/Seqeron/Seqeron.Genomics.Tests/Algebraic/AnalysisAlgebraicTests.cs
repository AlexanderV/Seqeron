using FsCheck;
using FsCheck.Fluent;
using FsCheck.NUnit;

namespace Seqeron.Genomics.Tests.Algebraic;

/// <summary>
/// Algebraic-law tests for the Analysis area (genomic similarity).
///
/// Algebraic testing pins the reflexive maximum and symmetry of the k-mer
/// Jaccard similarity.
/// — docs/checklists/06_ALGEBRAIC_TESTING.md §Description, row 179.
/// </summary>
[TestFixture]
[Category("Algebraic")]
[Category("Analysis")]
public class AnalysisAlgebraicTests
{
    private static Arbitrary<string> DnaArbitrary(int minLen) =>
        Gen.Elements('A', 'C', 'G', 'T').ArrayOf().Where(a => a.Length >= minLen)
            .Select(a => new string(a)).ToArbitrary();

    // ═══════════════════════════════════════════════════════════════════════
    // Unit: GENOMIC-SIMILARITY-001 — k-mer Jaccard similarity (Analysis), row 179.
    //
    // Model: similarity = |A∩B| / |A∪B| × 100 over the two sequences' k-mer sets —
    //        a symmetric Jaccard index whose self-comparison is the maximum (100).
    //   — docs/algorithms/Genomic_Analysis; GenomicAnalyzer.CalculateSimilarity.
    //
    // Laws (row 179): ID — sim(x, x) = 100 (max; Jaccard 1 on a sequence's own
    //                 k-mer set).  COMM — sim(a, b) = sim(b, a).
    // ═══════════════════════════════════════════════════════════════════════

    [FsCheck.NUnit.Property]
    public Property Similarity_Identity_SelfIsMax()
    {
        return Prop.ForAll(DnaArbitrary(minLen: 6), seq =>
        {
            double sim = GenomicAnalyzer.CalculateSimilarity(new DnaSequence(seq), new DnaSequence(seq));
            return (System.Math.Abs(sim - 100.0) < 1e-9).Label($"sim(x,x)={sim} for \"{seq}\"");
        });
    }

    [FsCheck.NUnit.Property]
    public Property Similarity_Commutative_Symmetric()
    {
        return Prop.ForAll(DnaArbitrary(minLen: 6), DnaArbitrary(minLen: 6), (a, b) =>
        {
            double ab = GenomicAnalyzer.CalculateSimilarity(new DnaSequence(a), new DnaSequence(b));
            double ba = GenomicAnalyzer.CalculateSimilarity(new DnaSequence(b), new DnaSequence(a));
            return (System.Math.Abs(ab - ba) < 1e-9).Label($"sim(a,b)={ab} != sim(b,a)={ba}");
        });
    }
}
