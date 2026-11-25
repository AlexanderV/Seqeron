// PRIMER-TM-001 — Self-dimer / hetero-dimer (intermolecular) Tm via thermodynamic alignment (opt-in)
// Evidence: docs/Evidence/PRIMER-TM-001-DIMER-Evidence.md
// TestSpec: tests/TestSpecs/PRIMER-TM-001-DIMER.md
// Sources:
//   SantaLucia J, Hicks D (2004). Annu Rev Biophys 33:415-440 — unified NN ΔH°/ΔS° (Table 1),
//     bimolecular Tm Eq. 3 (x=1 self-comp / x=4 non), Eq. 5 entropy salt correction (0.368).
//   Untergasser A et al. (2012). Nucleic Acids Res 40:e115 — Primer3/ntthal thermodynamic alignment.
//   primer3-py 2.3.0 (calc_homodimer / calc_heterodimer; mv=50, dv=0, dntp=0, dna_conc=50 nM)
//     — reference ΔH/ΔS/Tm captured this session (ntthal engine).

namespace Seqeron.Genomics.Tests.Unit.MolTools;

/// <summary>
/// Tests for the opt-in self-/hetero-dimer Tm added under PRIMER-TM-001. The MUST cases use
/// values hand-derived from the SantaLucia &amp; Hicks (2004) Table 1 NN parameters
/// (independent of the implementation) and/or captured from primer3-py 2.3.0 (the ntthal
/// reference). Hand-derived ΔH°/ΔS°/Tm are asserted to 1e-9; primer3 parity to 1e-3.
/// </summary>
[TestFixture]
public class PrimerDesigner_DimerTm_Tests
{
    private const double ExactTol = 1e-9;   // hand-derived (matches both derivation and primer3)
    private const double ParityTol = 1e-3;  // primer3-py reference rounded to 4 dp

    // primer3-py default dimer conditions captured this session.
    private const double Na = 0.05;     // 50 mM monovalent
    private const double Ct = 50e-9;    // 50 nM total strand (ntthal dna_conc)


    // M1 — GCGCGCGC self-dimer, full 8-bp Watson-Crick duplex (palindrome → x=1).
    // Hand-derived from SantaLucia & Hicks (2004) Table 1:
    //   stacks 4·GC(-9.8,-24.4) + 3·CG(-10.6,-27.2); init(+0.2,-5.7); no A·T end;
    //   ΔH° = -70.8 kcal/mol; ΔS° = -184.9 + 7·0.368·ln(0.05) = -192.61700633667505;
    //   ΔG°37 = ΔH° - 310.15·ΔS°/1000 = -11.059835484680235.
    [Test]
    public void FindMostStableDimer_GcgcgcgcSelfDimer_MatchesHandDerivedThermodynamics()
    {
        var result = PrimerDesigner.FindMostStableDimer("GCGCGCGC", "GCGCGCGC", Na, Ct);

        Assert.That(result, Is.Not.Null, "GCGCGCGC forms a fully complementary self-dimer.");
        var d = result!.Value;
        Assert.Multiple(() =>
        {
            Assert.That(d.BasePairs, Is.EqualTo(8),
                "The most stable self-dimer is the full 8-bp Watson-Crick duplex.");
            Assert.That(d.DeltaH, Is.EqualTo(-70.8).Within(ExactTol),
                "ΔH° = init(+0.2) + 4·GC(-9.8) + 3·CG(-10.6) = -70.8 kcal/mol (salt-independent).");
            Assert.That(d.DeltaS, Is.EqualTo(-192.61700633667505).Within(ExactTol),
                "ΔS° = -5.7 + 4·(-24.4) + 3·(-27.2) + 7·0.368·ln(0.05) (Eq. 5 salt correction).");
            Assert.That(d.DeltaG37, Is.EqualTo(-11.059835484680235).Within(ExactTol),
                "ΔG°37 = ΔH° - 310.15·ΔS°/1000.");
        });
    }

    // M12 — alignment spans: the optimal self-dimer starts at index 0 on both strands.
    [Test]
    public void FindMostStableDimer_GcgcgcgcSelfDimer_ReportsFullAlignmentSpans()
    {
        var d = PrimerDesigner.FindMostStableDimer("GCGCGCGC", "GCGCGCGC", Na, Ct)!.Value;
        Assert.Multiple(() =>
        {
            Assert.That(d.Strand1Start, Is.EqualTo(0), "Duplex covers strand 1 from index 0.");
            Assert.That(d.Strand2Start, Is.EqualTo(0), "Duplex covers strand 2 from index 0.");
        });
    }

    // M3 — TGCATGCATG / CATGCATGCA: non-palindromic hetero-dimer (x=4).
    // Hand-derived: stacks TG,GC,CA,AT,TG,GC,CA,AT,TG; init(+0.2,-5.7); one A·T end (5'-T).
    //   ΔH° = -74.1 kcal/mol; ΔS° = -211.8218652900108; Tm(x=4) = 25.659587124835923 °C.
    [Test]
    public void FindMostStableDimer_NonPalindromicHeteroDimer_MatchesHandDerivedThermodynamics()
    {
        var result = PrimerDesigner.FindMostStableDimer("TGCATGCATG", "CATGCATGCA", Na, Ct);

        Assert.That(result, Is.Not.Null, "These strands form a 10-bp Watson-Crick hetero-dimer.");
        var d = result!.Value;
        Assert.Multiple(() =>
        {
            Assert.That(d.BasePairs, Is.EqualTo(10), "Full 10-bp duplex.");
            Assert.That(d.DeltaH, Is.EqualTo(-74.1).Within(ExactTol),
                "ΔH° = init + Σ stacks + one terminal-A·T penalty (5'-T end).");
            Assert.That(d.DeltaS, Is.EqualTo(-211.8218652900108).Within(ExactTol),
                "ΔS° includes one A·T-penalty (+6.9) and 9·0.368·ln(0.05) salt correction.");
        });
    }

    // C1 — most-stable selection: a longer complementary core is preferred over a short one.
    // 'GCGCGCGCAA' / 'GCGCGCGCAA' can form a short 2-bp AA-edge pairing or the strong 8-bp
    // GCGCGCGC core; the 8-bp core has the higher Tm and must be the one returned.
    [Test]
    public void FindMostStableDimer_PrefersStrongerLongerDuplex()
    {
        var d = PrimerDesigner.FindMostStableDimer("GCGCGCGCAA", "GCGCGCGCAA", Na, Ct)!.Value;
        Assert.That(d.BasePairs, Is.EqualTo(8),
            "The 8-bp GCGCGCGC core is the highest-Tm duplex; the AA tail cannot self-pair.");
    }



    // Reference values: primer3-py 2.3.1 calc_heterodimer / calc_homodimer / calc_end_stability
    // (thal.c vendored at libnano/primer3-py tag v2.3.1), captured 2026-10-01; ΔG in cal/mol.
    // Before the fix (LSH/RSH tstack2-vs-A·T selection, thal.c keeps T1 = −∞ unless a dangling-end
    // branch is taken) these optima were off by e.g. ΔG +13.8 / −208.5 / +150 / −1014 cal/mol.
    private const double P3Mv = 0.050, P3Dv = 0.0015, P3Dntp = 0.0006;

}
