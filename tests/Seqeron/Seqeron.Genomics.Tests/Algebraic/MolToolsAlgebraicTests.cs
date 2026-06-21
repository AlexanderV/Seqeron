using System.Linq;
using FsCheck;
using FsCheck.Fluent;
using FsCheck.NUnit;

namespace Seqeron.Genomics.Tests.Algebraic;

/// <summary>
/// Algebraic-law tests for the MolTools area (CRISPR off-target scoring,
/// restriction digestion).
///
/// Algebraic testing pins the formal equations the scoring/transform functions
/// must obey: the identity value on a perfect/empty input and the monotonic /
/// conservation structure as the input degrades or is partitioned.
/// — docs/checklists/06_ALGEBRAIC_TESTING.md §Description, rows 20, 27.
/// </summary>
[TestFixture]
[Category("Algebraic")]
[Category("MolTools")]
public class MolToolsAlgebraicTests
{
    private static Arbitrary<string> Guide20Arbitrary() =>
        Gen.Elements('A', 'C', 'G', 'T')
            .ArrayOf(20)
            .Select(a => new string(a))
            .ToArbitrary();

    private static char NextBase(char c) => c switch
    {
        'A' => 'C', 'C' => 'G', 'G' => 'T', 'T' => 'A',
        _ => 'A'
    };

    // ═══════════════════════════════════════════════════════════════════════
    // Unit: CRISPR-OFF-001 — Off-target scoring (MolTools)
    // Checklist: docs/checklists/06_ALGEBRAIC_TESTING.md, row 20.
    //
    // Model: the Doench 2016 CFD off-target score is the product over the 20
    //        protospacer positions of a per-position mismatch penalty (∈ [0,1],
    //        exactly 1.0 at a matched position) times the PAM-activity score. A
    //        perfect 20/20 match against a canonical NGG PAM scores exactly 1.0.
    //        The MIT/Hsu 2013 single-hit score scores a perfect match at 100.
    //   — docs/algorithms/MolTools/Off_Target_Analysis.md; Doench 2016 / Hsu 2013;
    //     CrisprDesigner.CalculateCfdScore / CalculateMitHitScore.
    //
    // Laws under test (checklist row 20):
    //   • ID   — 0 mismatches → max score: CFD(g, g, NGG) = 1.0 and the MIT hit
    //            score of a self-comparison is 100.
    //   • DIST — score decreases monotonically as mismatches accumulate: because
    //            CFD multiplies one extra penalty factor ∈ [0,1] per additional
    //            mismatch, growing the mismatch set never increases the score.
    // ═══════════════════════════════════════════════════════════════════════

    /// <summary>
    /// ID: a guide scored against itself is a perfect match — CFD = 1.0 (canonical
    /// NGG PAM) and MIT hit score = 100, the maxima of the two scales.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property OffTarget_Identity_PerfectMatchIsMaxScore()
    {
        return Prop.ForAll(Guide20Arbitrary(), guide =>
        {
            double cfd = CrisprDesigner.CalculateCfdScore(guide, guide, "AGG");
            double mit = CrisprDesigner.CalculateMitHitScore(guide, guide);
            return (Math.Abs(cfd - 1.0) < 1e-12 && Math.Abs(mit - 100.0) < 1e-9)
                .Label($"perfect match not max: CFD={cfd}, MIT={mit} for \"{guide}\"");
        });
    }

    /// <summary>
    /// DIST: CFD is monotonically non-increasing as mismatches accumulate. Starting
    /// from the perfect off-target (= guide), flip one more position to a real
    /// mismatch at each step and assert the score never rises along the chain.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property OffTarget_Distributive_MonotoneNonIncreasingInMismatches()
    {
        return Prop.ForAll(Guide20Arbitrary(), guide =>
        {
            var current = guide.ToCharArray();
            double prev = CrisprDesigner.CalculateCfdScore(guide, new string(current), "AGG");
            for (int i = 0; i < 20; i++)
            {
                current[i] = NextBase(guide[i]); // introduce one more guaranteed mismatch
                double next = CrisprDesigner.CalculateCfdScore(guide, new string(current), "AGG");
                if (next > prev + 1e-12)
                    return false.Label($"CFD rose from {prev} to {next} at mismatch #{i + 1} for \"{guide}\"");
                prev = next;
            }
            return true.ToProperty();
        });
    }

    /// <summary>
    /// DIST witness: one mismatch strictly lowers CFD below the perfect 1.0, so the
    /// monotonicity law is exercised non-vacuously against a concrete drop.
    /// </summary>
    [Test]
    public void OffTarget_Distributive_OneMismatchLowersScore()
    {
        const string guide = "GACGTTGCAACGTTGCAACG";
        double perfect = CrisprDesigner.CalculateCfdScore(guide, guide, "AGG");
        var oneOff = guide.ToCharArray();
        oneOff[10] = NextBase(guide[10]);
        double mutated = CrisprDesigner.CalculateCfdScore(guide, new string(oneOff), "AGG");

        perfect.Should().Be(1.0);
        mutated.Should().BeLessThan(perfect);
        mutated.Should().BeGreaterThanOrEqualTo(0.0);
    }

    // ═══════════════════════════════════════════════════════════════════════
    // Unit: RESTR-DIGEST-001 — Restriction digest simulation (MolTools)
    // Checklist: docs/checklists/06_ALGEBRAIC_TESTING.md, row 27.
    //
    // Model: a linear digest with k distinct forward-strand cut positions splits
    //        the sequence into k+1 contiguous, non-overlapping fragments that
    //        tile the molecule end to end. With zero cuts the whole sequence is a
    //        single fragment.
    //   — docs/algorithms/MolTools/Restriction_Digest_Simulation.md;
    //     RestrictionAnalyzer.Digest.
    //
    // Laws under test (checklist row 27):
    //   • ID   — 0 cut sites → exactly 1 fragment equal to the full sequence.
    //   • DIST — conservation of length: Σ fragment.Length = sequence.Length, and
    //            (stronger) the ordered concatenation of fragments reconstructs
    //            the original sequence — no base is created or destroyed by the cut.
    // ═══════════════════════════════════════════════════════════════════════

    private static Arbitrary<string> DnaArbitrary(int minLen) =>
        Gen.Elements('A', 'C', 'G', 'T')
            .ArrayOf()
            .Where(a => a.Length >= minLen)
            .Select(a => new string(a))
            .ToArbitrary();

    /// <summary>
    /// ID: a sequence with no EcoRI site (no GAATTC) digests into exactly one
    /// fragment equal to the whole input.
    /// </summary>
    [Test]
    public void Digest_Identity_NoSiteYieldsWholeSequence()
    {
        var seq = new DnaSequence("AAAACCCCGGGGTTTTAAAACCCC");
        var fragments = RestrictionAnalyzer.Digest(seq, "EcoRI").ToList();
        fragments.Should().HaveCount(1);
        fragments[0].Sequence.Should().Be(seq.Sequence);
        fragments[0].Length.Should().Be(seq.Length);
    }

    /// <summary>
    /// DIST: length conservation — the EcoRI fragments always sum to the input
    /// length and concatenate back to the input, for any DNA (cut or not).
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Digest_Distributive_LengthIsConserved()
    {
        return Prop.ForAll(DnaArbitrary(minLen: 1), seq =>
        {
            var dna = new DnaSequence(seq);
            var fragments = RestrictionAnalyzer.Digest(dna, "EcoRI").ToList();
            int sum = fragments.Sum(f => f.Length);
            string reconstructed = string.Concat(fragments.Select(f => f.Sequence));
            return (sum == dna.Length && reconstructed == dna.Sequence)
                .Label($"sum={sum}, len={dna.Length}, reconstructed==input:{reconstructed == dna.Sequence}");
        });
    }

    /// <summary>
    /// DIST witness: EcoRI cuts AAAGAATTCAAA once, giving two fragments whose
    /// lengths sum to 12 and which concatenate to the original.
    /// </summary>
    [Test]
    public void Digest_Distributive_WorkedCut()
    {
        var seq = new DnaSequence("AAAGAATTCAAA");
        var fragments = RestrictionAnalyzer.Digest(seq, "EcoRI").ToList();
        fragments.Count.Should().BeGreaterThan(1);
        fragments.Sum(f => f.Length).Should().Be(seq.Length);
        string.Concat(fragments.Select(f => f.Sequence)).Should().Be(seq.Sequence);
    }
}
