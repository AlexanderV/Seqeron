namespace Seqeron.Genomics.Tests.Algebraic;

using BasePair = RnaSecondaryStructure.BasePair;
using BasePairType = RnaSecondaryStructure.BasePairType;

/// <summary>
/// Algebraic-law tests for the RnaStructure area (free energy, dot-bracket I/O,
/// base pairing).
///
/// Algebraic testing pins the additive nearest-neighbour structure of the free
/// energy, the parse∘format round-trip of dot-bracket notation, and the symmetry
/// of the base-pairing relation.
/// — docs/checklists/06_ALGEBRAIC_TESTING.md §Description, rows 73, 149, 153.
/// </summary>
[TestFixture]
[Category("Algebraic")]
[Category("RnaStructure")]
public class RnaStructureAlgebraicTests
{
    private const string RnaBases = "ACGU";

    // ═══════════════════════════════════════════════════════════════════════
    // Unit: RNA-ENERGY-001 — RNA free energy (RnaStructure)
    // Checklist: docs/checklists/06_ALGEBRAIC_TESTING.md, row 73.
    //
    // Model: the Turner nearest-neighbour free energy of a helix is the sum of its
    //        independent stacking terms (plus terminal penalties); an unstructured
    //        sequence (no base pairs) has ΔG = 0.
    //   — docs/algorithms/RNA_Structure; RnaSecondaryStructure.CalculateMinimumFreeEnergy /
    //     CalculateStemEnergy.
    //
    // Laws under test (checklist row 73):
    //   • ID   — no structure → ΔG = 0 (a homopolymer cannot pair; empty stem = 0).
    //   • DIST — the helix energy is additive over independent stacks: a uniform
    //            G-C stem of k pairs has energy (k−1) × (single-stack energy).
    // ═══════════════════════════════════════════════════════════════════════

    /// <summary>ID: an unstructured sequence has zero free energy.</summary>
    [Test]
    public void Energy_Identity_NoStructureIsZero()
    {
        // Poly-A cannot form any base pair (A pairs only with U), so ΔG = 0.
        RnaSecondaryStructure.CalculateMinimumFreeEnergy("AAAAAAAAAAAA").Should().Be(0.0);
        // An empty base-pair set has zero stem energy.
        RnaSecondaryStructure.CalculateStemEnergy("", new List<BasePair>()).Should().Be(0.0);
    }

    private static List<BasePair> UniformGcStem(int k)
    {
        var pairs = new List<BasePair>();
        for (int i = 0; i < k; i++)
            pairs.Add(new BasePair(i, 2 * k - 1 - i, 'G', 'C', BasePairType.WatsonCrick));
        return pairs;
    }

    /// <summary>
    /// DIST: stacking energy is additive — a uniform G-C stem of k pairs has
    /// (k−1) identical, independent stacking contributions, so its energy is
    /// exactly (k−1) times the single-stack energy.
    /// </summary>
    [Test]
    public void Energy_Distributive_AdditiveOverIndependentStacks()
    {
        double singleStack = RnaSecondaryStructure.CalculateStemEnergy("", UniformGcStem(2)); // 1 stack
        singleStack.Should().BeLessThan(0); // a G-C stack is stabilizing

        for (int k = 2; k <= 8; k++)
        {
            double energy = RnaSecondaryStructure.CalculateStemEnergy("", UniformGcStem(k));
            energy.Should().BeApproximately((k - 1) * singleStack, 1e-9,
                $"a {k}-pair G-C stem should have (k-1) independent stacks");
        }
    }

    // ═══════════════════════════════════════════════════════════════════════
    // Unit: RNA-DOTBRACKET-001 — Dot-bracket notation (RnaStructure)
    // Checklist: docs/checklists/06_ALGEBRAIC_TESTING.md, row 149.
    //
    // Model: dot-bracket notation encodes a non-crossing secondary structure as
    //        matched ()/. symbols; for a nested structure the notation is unique,
    //        so reconstructing it from the parsed pairs reproduces the input string.
    //   — docs/algorithms/RnaStructure/Dot_Bracket_Notation.md;
    //     RnaSecondaryStructure.ParseDotBracket / ValidateDotBracket.
    //
    // Laws under test (checklist row 149):
    //   • RT — format(parse(s)) = s for canonical nested notation (parse inverts
    //          the unique encoding).
    //   • ID — empty / null input → no pairs (a valid, pair-free structure).
    // ═══════════════════════════════════════════════════════════════════════

    /// <summary>Reconstructs canonical nested dot-bracket from a pair set (the
    /// unique inverse for non-crossing structures).</summary>
    private static string Format(IEnumerable<(int P1, int P2)> pairs, int length)
    {
        var s = new char[length];
        System.Array.Fill(s, '.');
        foreach (var (p1, p2) in pairs)
        {
            s[p1] = '(';
            s[p2] = ')';
        }
        return new string(s);
    }

    [Test]
    public void DotBracket_RoundTrip_FormatOfParseIsIdentity()
    {
        foreach (var s in new[] { "((..))", "(())", "((((....))))", "(.(.).)", ".((...)).", "...", "()()" })
        {
            RnaSecondaryStructure.ValidateDotBracket(s).Should().BeTrue($"\"{s}\" is well-formed");
            var pairs = RnaSecondaryStructure.ParseDotBracket(s).ToList();
            Format(pairs, s.Length).Should().Be(s, $"round-trip of \"{s}\"");
        }
    }

    [Test]
    public void DotBracket_Identity_EmptyHasNoPairs()
    {
        RnaSecondaryStructure.ParseDotBracket("").Should().BeEmpty();
        RnaSecondaryStructure.ParseDotBracket(null!).Should().BeEmpty();
        RnaSecondaryStructure.ValidateDotBracket("").Should().BeTrue();
    }

    // ═══════════════════════════════════════════════════════════════════════
    // Unit: RNA-PAIR-001 — Base pairing relation (RnaStructure)
    // Checklist: docs/checklists/06_ALGEBRAIC_TESTING.md, row 153.
    //
    // Model: the can-pair relation is symmetric (Watson–Crick A-U, G-C and wobble
    //        G-U pair regardless of order) and is a deterministic pure function.
    //   — docs/algorithms/RnaStructure; RnaSecondaryStructure.CanPair / GetBasePairType.
    //
    // Laws under test (checklist row 153):
    //   • COMM  — canPair(a, b) = canPair(b, a); the pair type is order-independent.
    //   • IDEMP — deterministic: repeated calls give identical results.
    // ═══════════════════════════════════════════════════════════════════════

    [Test]
    public void CanPair_Commutative_Symmetric()
    {
        foreach (char a in RnaBases)
            foreach (char b in RnaBases)
            {
                RnaSecondaryStructure.CanPair(a, b)
                    .Should().Be(RnaSecondaryStructure.CanPair(b, a), $"canPair('{a}','{b}')");
                RnaSecondaryStructure.GetBasePairType(a, b)
                    .Should().Be(RnaSecondaryStructure.GetBasePairType(b, a), $"type('{a}','{b}')");
            }
    }

    [Test]
    public void CanPair_Idempotent_Deterministic()
    {
        foreach (char a in RnaBases)
            foreach (char b in RnaBases)
                RnaSecondaryStructure.CanPair(a, b).Should().Be(RnaSecondaryStructure.CanPair(a, b));
    }

    /// <summary>Witness: the three canonical pairings hold and a non-pair (A-G) does not.</summary>
    [Test]
    public void CanPair_WorkedExamples()
    {
        RnaSecondaryStructure.CanPair('A', 'U').Should().BeTrue();
        RnaSecondaryStructure.CanPair('G', 'C').Should().BeTrue();
        RnaSecondaryStructure.CanPair('G', 'U').Should().BeTrue(); // wobble
        RnaSecondaryStructure.CanPair('A', 'G').Should().BeFalse();
    }

    // ═══════════════════════════════════════════════════════════════════════
    // Unit: RNA-ACCESS-001 — McCaskill unpaired (accessibility) probabilities (RnaStructure)
    // Checklist: docs/checklists/06_ALGEBRAIC_TESTING.md, row 238.
    //
    // Model: the equilibrium probability that a window is wholly unpaired is the exact
    //        McCaskill ensemble ratio Z_open/Z (Bernhart et al. 2006; McCaskill 1990).
    //        The atomic (length-1) region accessibility equals the per-base unpaired
    //        probability, and a region's joint unpaired probability factorises into the
    //        per-base product ONLY in the independence floor (no pair can form) — there
    //        the whole region, like every base, is unpaired with probability 1.
    //   — docs/Validation/reports/RNA-ACCESS-001.md; RnaSecondaryStructure.
    //     CalculateUnpairedProbabilities / CalculateRegionUnpairedProbability.
    //
    // Laws under test (checklist row 238):
    //   • ID    — under the independence floor (no pairable bases) the full-length region
    //             unpaired probability equals the product over bases (= 1).
    //   • DIST  — region(L=1 @ i) = per-base p_unpaired(i): the atomic accessibility
    //             consistency identity (both = Z_forbid(i)/Z).
    //   • IDEMP — deterministic: repeated evaluation gives the identical probabilities.
    // ═══════════════════════════════════════════════════════════════════════

    /// <summary>
    /// ID: a poly-A sequence has no Watson–Crick/GU pair, so Z = 1 and every base is
    /// unpaired with probability 1 (the independence floor). There the full-length
    /// region's joint unpaired probability equals the product of the per-base ones (= 1).
    /// </summary>
    [Test]
    public void Accessibility_Identity_FloorRegionEqualsProductOverBases()
    {
        const string polyA = "AAAAAAAAAAAAAAAAAAAA"; // 20 nt, no pairable bases
        int n = polyA.Length;

        var perBase = RnaSecondaryStructure.CalculateUnpairedProbabilities(polyA).UnpairedProbabilities;
        double product = perBase.Aggregate(1.0, (acc, p) => acc * p);

        double fullRegion = RnaSecondaryStructure.CalculateRegionUnpairedProbability(
            polyA, windowEnd: n - 1, windowLength: n);

        product.Should().BeApproximately(1.0, 1e-12, "every base is unpaired with probability 1 under the floor");
        fullRegion.Should().BeApproximately(product, 1e-12,
            "the full region's joint unpaired probability factorises into the per-base product under independence");
    }

    /// <summary>
    /// DIST: the length-1 region accessibility ending at i is, by definition, the per-base
    /// unpaired probability of i — both are the exact ratio Z_forbid(i)/Z. Verified on a
    /// genuinely structured sequence (a G-C hairpin) where the probabilities are non-trivial.
    /// </summary>
    [Test]
    public void Accessibility_Distributive_UnitWindowEqualsPerBase()
    {
        const string hairpin = "GGGGGAAAACCCCC";
        var perBase = RnaSecondaryStructure.CalculateUnpairedProbabilities(hairpin).UnpairedProbabilities;

        // A genuinely folded sequence: at least one base must be substantially paired.
        perBase.Min().Should().BeLessThan(0.99, "the hairpin must actually fold");

        for (int i = 0; i < hairpin.Length; i++)
        {
            double unit = RnaSecondaryStructure.CalculateRegionUnpairedProbability(
                hairpin, windowEnd: i, windowLength: 1);
            unit.Should().BeApproximately(perBase[i], 1e-9, $"region(L=1 @ {i}) must equal p_unpaired({i})");
        }
    }

    [Test]
    public void Accessibility_Idempotent_Deterministic()
    {
        const string seq = "GGGACAUGUCCCAAGG";

        var a = RnaSecondaryStructure.CalculateUnpairedProbabilities(seq).UnpairedProbabilities;
        var b = RnaSecondaryStructure.CalculateUnpairedProbabilities(seq).UnpairedProbabilities;
        a.Should().Equal(b);

        double r1 = RnaSecondaryStructure.CalculateRegionUnpairedProbability(seq, windowEnd: 10, windowLength: 8);
        double r2 = RnaSecondaryStructure.CalculateRegionUnpairedProbability(seq, windowEnd: 10, windowLength: 8);
        r2.Should().Be(r1);
    }
}
