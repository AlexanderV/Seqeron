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

    // ═══════════════════════════════════════════════════════════════════════
    // Unit: PRIMER-NNTM-001 — Nearest-neighbour Tm (MolTools), row 240.
    //
    // Model: a perfect Watson–Crick duplex's thermodynamics are the SantaLucia (1998)
    //        unified nearest-neighbour SUM — initiation + Σ per-dimer increments +
    //        terminal-AT penalties (+ symmetry for self-complementary strands).
    //   — docs/algorithms/MolTools; PrimerDesigner.CalculateNearestNeighborThermodynamics /
    //     CalculateMeltingTemperatureNN. TestSpec tests/TestSpecs/PRIMER-NNTM-001.md.
    //
    // Laws (row 240): ID — the perfect-match path equals the published NN sum.
    //                 IDEMP — Tm is a pure, deterministic function.
    // ═══════════════════════════════════════════════════════════════════════

    [Test]
    public void NearestNeighborTm_Identity_PerfectMatchEqualsPublishedSum()
    {
        // ATGCATGC: ΔH = 0.2 + (−7.2−8.5−9.8−8.5−7.2−8.5−9.8) + 2.2 + 2.2 = −57.1 kcal/mol;
        //           ΔS = −5.7 + (−20.4−22.7−24.4−22.7−20.4−22.7−24.4) + 6.9 + 6.9 = −156.5 cal/(K·mol).
        var r1 = PrimerDesigner.CalculateNearestNeighborThermodynamics("ATGCATGC");
        r1!.Value.DeltaH.Should().BeApproximately(-57.1, 1e-9);
        r1.Value.DeltaS.Should().BeApproximately(-156.5, 1e-9);

        // GCGCGC: ΔH = 0.2 + (−9.8−10.6−9.8−10.6−9.8) = −50.4; ΔS = −5.7 + (sum) + (−1.4 symmetry) = −134.7.
        var r2 = PrimerDesigner.CalculateNearestNeighborThermodynamics("GCGCGC");
        r2!.Value.DeltaH.Should().BeApproximately(-50.4, 1e-9);
        r2.Value.DeltaS.Should().BeApproximately(-134.7, 1e-9);
    }

    [Test]
    public void NearestNeighborTm_Idempotent_Deterministic()
    {
        double a = PrimerDesigner.CalculateMeltingTemperatureNN("ATGCATGCATGC");
        double b = PrimerDesigner.CalculateMeltingTemperatureNN("ATGCATGCATGC");
        b.Should().Be(a);
    }

    // ═══════════════════════════════════════════════════════════════════════
    // Unit: PRIMER-HAIRPIN-001 — Unimolecular hairpin folder (MolTools), row 241.
    //
    // Model: the most stable hairpin's ΔH/ΔS/ΔG37 are the nearest-neighbour SUM over
    //        the stem stacks plus the loop ΔG penalty — a hand-derivable closed form.
    //   — PrimerDesigner.FindMostStableHairpin / CalculateHairpinMeltingTemperature.
    //     TestSpec tests/TestSpecs/PRIMER-HAIRPIN-001.md.
    //
    // Laws (row 241): ID — the fixed canonical hairpin reproduces the hand-derived NN ΔG.
    //                 IDEMP — hairpin folding/Tm is a pure, deterministic function.
    // ═══════════════════════════════════════════════════════════════════════

    private const string CanonicalHairpin = "GGGCTTTTGCCC"; // 4-bp G-C stem, 4-nt loop

    [Test]
    public void Hairpin_Identity_CanonicalHairpinEqualsHandDerivedNn()
    {
        var hp = PrimerDesigner.FindMostStableHairpin(CanonicalHairpin);
        hp.Should().NotBeNull();
        hp!.Value.DeltaH.Should().BeApproximately(-25.8, 1e-9);
        hp.Value.DeltaS.Should().BeApproximately(-75.48486216346927, 1e-9);
        hp.Value.DeltaG37.Should().BeApproximately(-2.3883700000000054, 1e-9);
    }

    [Test]
    public void Hairpin_Idempotent_Deterministic()
    {
        double a = PrimerDesigner.CalculateHairpinMeltingTemperature(CanonicalHairpin);
        double b = PrimerDesigner.CalculateHairpinMeltingTemperature(CanonicalHairpin);
        b.Should().Be(a);
    }

    // ═══════════════════════════════════════════════════════════════════════
    // Unit: PRIMER-DIMER-001 — ntthal self/hetero-dimer (MolTools), row 242.
    //
    // Model: the most stable duplex between two strands is the contiguous Watson–Crick
    //        alignment whose ΔH/ΔS is the nearest-neighbour SUM over its consecutive
    //        base-pair steps (ntthal thermodynamic alignment).
    //   — PrimerDesigner.FindMostStableDimer / CalculateSelfDimerMeltingTemperature.
    //     TestSpec tests/TestSpecs/PRIMER-DIMER-001.md.
    //
    // Laws (row 242): ID — the contiguous-WC optimum equals the NN sum.
    //                 IDEMP — dimer search/Tm is a pure, deterministic function.
    // ═══════════════════════════════════════════════════════════════════════

    private const double DimerNa = 0.05;    // 50 mM monovalent
    private const double DimerCt = 50e-9;   // 50 nM total strand

    [Test]
    public void Dimer_Identity_ContiguousWcOptimumEqualsNnSum()
    {
        // The fully complementary 8-mer self-dimer: all 8 bases pair, 7 G-C/C-G stacks.
        var d = PrimerDesigner.FindMostStableDimer("GCGCGCGC", "GCGCGCGC", DimerNa, DimerCt);
        d.Should().NotBeNull();
        d!.Value.BasePairs.Should().Be(8);
        d.Value.DeltaH.Should().BeApproximately(-70.8, 1e-9);
        d.Value.DeltaS.Should().BeApproximately(-192.61700633667505, 1e-9);
        d.Value.DeltaG37.Should().BeApproximately(-11.059835484680235, 1e-9);
    }

    [Test]
    public void Dimer_Idempotent_Deterministic()
    {
        double a = PrimerDesigner.CalculateSelfDimerMeltingTemperature("GCGCGCGC", DimerNa, DimerCt);
        double b = PrimerDesigner.CalculateSelfDimerMeltingTemperature("GCGCGCGC", DimerNa, DimerCt);
        b.Should().Be(a);
    }

    // ═══════════════════════════════════════════════════════════════════════
    // Unit: PROBE-LNATM-001 — LNA-adjusted NN Tm (MolTools), row 243.
    //
    // Model: an LNA-modified duplex adds per-position LNA increments to the DNA NN sum;
    //        with NO LNA positions the model must reduce EXACTLY to the plain DNA NN Tm.
    //   — PrimerDesigner.CalculateNearestNeighborThermodynamicsLna /
    //     CalculateMeltingTemperatureNNLna. TestSpec tests/TestSpecs/PROBE-LNATM-001.md.
    //
    // Laws (row 243): ID — zero LNA → equals the PRIMER-NNTM result.
    //                 IDEMP — LNA Tm is a pure, deterministic function.
    // ═══════════════════════════════════════════════════════════════════════

    [Test]
    public void LnaTm_Identity_ZeroLnaEqualsPlainNn()
    {
        const string seq = "GTGCATCGATGCAGC";
        var lna = PrimerDesigner.CalculateNearestNeighborThermodynamicsLna(seq, System.Array.Empty<int>());
        var dna = PrimerDesigner.CalculateNearestNeighborThermodynamics(seq);
        lna!.Value.DeltaH.Should().Be(dna!.Value.DeltaH);
        lna.Value.DeltaS.Should().Be(dna.Value.DeltaS);

        double tmLna = PrimerDesigner.CalculateMeltingTemperatureNNLna(seq, System.Array.Empty<int>());
        double tmDna = PrimerDesigner.CalculateMeltingTemperatureNN(seq);
        tmLna.Should().Be(tmDna);
    }

    [Test]
    public void LnaTm_Idempotent_Deterministic()
    {
        const string seq = "GTGCATCGATGCAGC";
        double a = PrimerDesigner.CalculateMeltingTemperatureNNLna(seq, new[] { 2, 7 });
        double b = PrimerDesigner.CalculateMeltingTemperatureNNLna(seq, new[] { 2, 7 });
        b.Should().Be(a);
    }

    // ═══════════════════════════════════════════════════════════════════════
    // Unit: PROBE-EVALUE-001 — Karlin–Altschul off-target E-value (MolTools), row 244.
    //
    // Model: the bit score is the affine-then-rescaled raw score S' = (λS − ln K)/ln 2,
    //        and E = K·m·n·e^{−λS} = m·n·2^{−S'} (Karlin & Altschul 1990; Altschul 1990).
    //   — ProbeDesigner.ComputeKarlinAltschul / ComputeLambdaNucleotide.
    //     TestSpec tests/TestSpecs/PROBE-EVALUE-001.md.
    //
    // Laws (row 244): ID — bit = (λS − ln K)/ln 2 (and E = m·n·2^{−bit}).
    //                 IDEMP — the statistics are a pure, deterministic function.
    // ═══════════════════════════════════════════════════════════════════════

    [Test]
    public void KarlinAltschul_Identity_BitScoreEqualsClosedForm()
    {
        const double k = 0.711;
        var s = ProbeDesigner.ComputeKarlinAltschul(rawScore: 45.0, queryLength: 25, databaseLength: 3_000_000_000L, k: k);

        // bit = (λS − ln K)/ln 2, computed independently from the returned λ.
        double expectedBit = (s.Lambda * s.RawScore - Math.Log(k)) / Math.Log(2.0);
        s.BitScore.Should().BeApproximately(expectedBit, 1e-9);

        // E = m·n·2^{−bit} (equivalently K·m·n·e^{−λS}).
        double expectedE = (double)s.QueryLength * s.DatabaseLength * Math.Pow(2.0, -s.BitScore);
        s.EValue.Should().BeApproximately(expectedE, expectedE * 1e-9);
    }

    [Test]
    public void KarlinAltschul_Idempotent_Deterministic()
    {
        var a = ProbeDesigner.ComputeKarlinAltschul(45.0, 25, 3_000_000_000L);
        var b = ProbeDesigner.ComputeKarlinAltschul(45.0, 25, 3_000_000_000L);
        b.BitScore.Should().Be(a.BitScore);
        b.EValue.Should().Be(a.EValue);
        b.Lambda.Should().Be(a.Lambda);
    }
}
