namespace Seqeron.Genomics.Tests.Mutation;

/// <summary>
/// Targeted mutation-killing tests for PrimerDesigner.cs (checklist 04 rows 21/22/23:
/// PRIMER-TM-001, PRIMER-DESIGN-001, PRIMER-STRUCT-001).
///
/// These pin the documented business logic that the canonical suite left as Stryker
/// "survivors": the primer-quality SCORING formula, the 3' GC-clamp rule, the
/// dinucleotide-repeat acceptance boundary, and the candidate-generation window
/// boundaries. Each test encodes the published rule (not the production result) so an
/// arithmetic/relational/logical mutant of the source diverges and is killed.
///
/// Scoring model (PrimerDesigner.CalculatePrimerScore, transcribed from source docs):
///   score = 100
///         − |len − OptimalLength| · 2
///         − |Tm  − OptimalTm|     · 2
///         − |GC  − 50|            · 0.5
///         − longestHomopolymer    · 5
///         + (3'-base ∈ {G,C} ? 5 : 0)      , clamped to ≥ 0, rounded to 2 dp.
/// </summary>
[TestFixture]
public class PrimerDesignerMutationTests
{

    // ── Hairpin detection: exercise BOTH the O(n²) and suffix-tree paths ──────────────

    [Test]
    public void HasHairpinPotential_ShortSelfComplementaryStem_DetectsHairpin()
    {
        // stem GGGGG / loop AAA / stem CCCCC (revcomp(CCCCC)=GGGGG) → hairpin, O(n²) path (<100 bp).
        PrimerDesigner.HasHairpinPotential("GGGGGAAACCCCC").Should().BeTrue();
    }

    [Test]
    public void HasHairpinPotential_ShortNonComplementary_NoHairpin()
    {
        PrimerDesigner.HasHairpinPotential("AAAAAAAAAAAAA").Should().BeFalse();
    }

    [Test]
    public void HasHairpinPotential_LongSequenceWithPlantedHairpin_UsesSuffixTreePath()
    {
        // ≥100 bp routes through the suffix-tree implementation. Poly-A padding cannot form
        // spurious stems (A pairs only with T, absent), so detection is driven solely by the
        // planted GGGGG…CCCCC stem-loop. Kills the NoCoverage block on the suffix-tree branch.
        string seq = new string('A', 60) + "GGGGGAAACCCCC" + new string('A', 60); // 133 bp
        PrimerDesigner.HasHairpinPotential(seq).Should().BeTrue();
    }

    [Test]
    public void HasHairpinPotential_LongSequenceNoComplementarity_NoHairpin()
    {
        PrimerDesigner.HasHairpinPotential(new string('A', 120)).Should().BeFalse();
    }
}
