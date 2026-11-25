namespace Seqeron.Genomics.Tests.Mutation;

/// <summary>
/// Targeted mutation-killing tests for ProbeDesigner.cs (checklist 04 rows 24/25:
/// PROBE-DESIGN-001, PROBE-VALID-001).
///
/// The canonical suite asserts invariant RANGES (score∈[0,1], GC∈[0,1]); it leaves the
/// exact formulas un-pinned, so Stryker arithmetic/relational/logical mutants survived
/// (37.6% baseline). These tests transcribe each documented formula/rule as ground truth
/// and assert the production value EQUALS it, so a mutated operator diverges and is killed.
/// </summary>
[TestFixture]
public class ProbeDesignerMutationTests
{
    // ── Oligo analysis: molecular weight = Biopython molecular_weight (single-stranded) ──

    [Test]
    [TestCase("ACGT", 1253.8027)]   // Biopython molecular_weight("ACGT", "DNA")
    [TestCase("A", 331.2218)]       // single base, no phosphodiester water
    [TestCase("AA", 644.4283)]      // 2·331.2218 − 18.01528
    [TestCase("AAA", 957.6348)]     // (n−1)·water term
    [TestCase("U", 324.1813)]       // RNA inferred (U, no T): UMP 324.1813 (was dUMP-like 308.2)
    public void CalculateMolecularWeight_MatchesBiopython(string seq, double expected)
    {
        ProbeDesigner.CalculateMolecularWeight(seq).Should().BeApproximately(expected, 1e-6);
    }

    [Test]
    public void CalculateMolecularWeight_ExplicitType_DelegatesToCanonical()
    {
        // Biopython molecular_weight("ACGU", "RNA") = 1303.7737; DNA table skips U.
        ProbeDesigner.CalculateMolecularWeight("acgu", isDna: false).Should().BeApproximately(1303.7737, 1e-6);
        ProbeDesigner.CalculateMolecularWeight("ACGT", isDna: true)
            .Should().Be(SequenceStatistics.CalculateNucleotideMolecularWeight("ACGT", true));
    }

    // ── Oligo analysis: extinction coefficient (Σ per-base ε260) ─────────────────────

    [Test]
    [TestCase("ACGT", 43000.0)] // 15400+7400+11500+8700
    [TestCase("AA", 30800.0)]   // 2·15400
    [TestCase("U", 9900.0)]     // uracil branch
    public void CalculateExtinctionCoefficient_MatchesPerBaseSum(string seq, double expected)
    {
        ProbeDesigner.CalculateExtinctionCoefficient(seq).Should().BeApproximately(expected, 1e-6);
    }

    // ── Oligo analysis: concentration (Beer–Lambert c = A/(ε·l)·1e6) ─────────────────

    [Test]
    [TestCase(1.0, 10000.0, 1.0, 100.0)]  // 1/10000·1e6
    [TestCase(3.0, 10000.0, 1.0, 300.0)]  // absorbance scales numerator
    [TestCase(1.0, 10000.0, 2.0, 50.0)]   // path length doubles denominator
    [TestCase(1.0, 20000.0, 1.0, 50.0)]   // extinction doubles denominator
    public void CalculateConcentration_MatchesBeerLambert(double a260, double eps, double path, double expected)
    {
        ProbeDesigner.CalculateConcentration(a260, eps, path).Should().BeApproximately(expected, 1e-6);
    }

    // ── ValidateProbe: specificity score (0 hits→0, 1→1, N→1/N) + off-target rule ────

    private const string Probe8 = "ACGTTGCA";
}
