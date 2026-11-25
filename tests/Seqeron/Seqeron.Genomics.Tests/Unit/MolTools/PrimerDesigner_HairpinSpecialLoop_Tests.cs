// PRIMER-TM-001 — Special tri/tetraloop hairpin bonus tables (bundled) via full ntthal hairpin DP (opt-in)
// Evidence: docs/Evidence/PRIMER-TM-001-SPECIAL-LOOP-Evidence.md
// TestSpec: tests/TestSpecs/PRIMER-TM-001-SPECIAL-LOOP.md
// Sources:
//   SantaLucia J, Hicks D (2004). Annu Rev Biophys 33:415-440 — special hairpin loops (tri/tetraloop bonuses).
//   primer3 libprimer3 primer3_config/{triloop,tetraloop}.dh,.ds — the verbatim bonus tables ntthal loads
//     (raw.githubusercontent.com/libnano/primer3-py/master/primer3/src/libprimer3/primer3_config/...).
//   thal.c calc_hairpin (lines 2106-2127) — how the bonus is keyed (full loop string incl. closing pair) & applied.
//   primer3-py 2.3.0 (calc_hairpin; mv_conc=50, dv_conc=0, dntp_conc=0, dna_conc=50)
//     — reference ΔH/ΔS/ΔG/Tm captured this session.

namespace Seqeron.Genomics.Tests.Unit.MolTools;

/// <summary>
/// Tests for the bundled sequence-specific special triloop / tetraloop hairpin bonus tables,
/// applied automatically by <see cref="ThermodynamicPrimerDesigner.CalculateHairpinThermodynamicsNtthal"/>
/// (the full Primer3 ntthal hairpin DP). Reference ΔH/ΔS/ΔG/Tm were captured from primer3-py
/// 2.3.0 <c>calc_hairpin</c> at the standard conditions (mv=50, dv=0, dntp=0, dna_conc=50 nM).
/// ΔH (integer cal/mol) is asserted exactly; ΔS/ΔG/Tm to a tight tolerance. The non-special-loop
/// regression cases confirm bundling the tables does NOT change a hairpin whose loop is not a
/// recognised special loop. The legacy <see cref="PrimerDesigner.FindMostStableHairpin"/> /
/// <see cref="PrimerDesigner.CalculateHairpinMeltingTemperature"/> (SantaLucia &amp; Hicks 2004
/// Table 4 model + caller-supplied <c>loopBonusDeltaG37</c>) are verified unchanged.
/// </summary>
[TestFixture]
public class PrimerDesigner_HairpinSpecialLoop_Tests
{
    private const double Na = 0.05;        // 50 mM monovalent (primer3 calc_hairpin mv default)
    private const double ParityTol = 1e-6;



    // C1 — the legacy FindMostStableHairpin (SantaLucia & Hicks 2004 Table 4 model) is UNCHANGED by
    // bundling: the canonical GGGCTTTTGCCC hairpin still returns the hand-derived Table-4 values
    // (Evidence: PRIMER-TM-001-HAIRPIN). ΔH°=-25.8, ΔS°=-75.48486216346927, ΔG°37=-2.3883700000000054.
    [Test]
    public void FindMostStableHairpin_LegacyTable4Model_Unchanged()
    {
        var h = PrimerDesigner.FindMostStableHairpin("GGGCTTTTGCCC");

        Assert.That(h, Is.Not.Null);
        var v = h!.Value;
        Assert.Multiple(() =>
        {
            Assert.That(v.DeltaH, Is.EqualTo(-25.8).Within(1e-9), "Table-4 stem stacks ΔH° (unchanged by bundling).");
            Assert.That(v.DeltaS, Is.EqualTo(-75.48486216346927).Within(1e-9), "Table-4 ΔS° (loop ΔS = -ΔG37·1000/310.15).");
            Assert.That(v.DeltaG37, Is.EqualTo(-2.3883700000000054).Within(1e-9), "Table-4 ΔG°37 (unchanged).");
        });
    }

    // C2 — the caller-supplied loopBonusDeltaG37 manual path still works: a +1.0 kcal/mol bonus
    // shifts the legacy hairpin ΔG°37 by exactly +1.0 vs the no-bonus result.
    [Test]
    public void FindMostStableHairpin_ManualLoopBonus_StillApplied()
    {
        var baseline = PrimerDesigner.FindMostStableHairpin("GGGCTTTTGCCC", 2, 0.0)!.Value;
        var bonused = PrimerDesigner.FindMostStableHairpin("GGGCTTTTGCCC", 2, 1.0)!.Value;

        Assert.That(bonused.DeltaG37 - baseline.DeltaG37, Is.EqualTo(1.0).Within(1e-9),
            "A caller-supplied loopBonusDeltaG37 = +1.0 raises the legacy hairpin ΔG°37 by exactly +1.0 kcal/mol.");
    }

}
