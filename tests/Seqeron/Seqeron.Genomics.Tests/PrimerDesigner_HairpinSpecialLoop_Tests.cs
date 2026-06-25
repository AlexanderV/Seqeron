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

using NUnit.Framework;
using Seqeron.Genomics.MolTools;

namespace Seqeron.Genomics.Tests;

/// <summary>
/// Tests for the bundled sequence-specific special triloop / tetraloop hairpin bonus tables,
/// applied automatically by <see cref="PrimerDesigner.CalculateHairpinThermodynamicsNtthal"/>
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
    private const double ParityTol = 1e-6; // primer3-py reference; ntthal port matches to machine precision

    #region CalculateHairpinThermodynamicsNtthal — special tetraloop (bundled bonus)

    // M1 — recognised GNRA tetraloop CGAAAG (closing C·G, loop GAAA). primer3-py calc_hairpin:
    //   dH=-40900 cal/mol, dS=-114.1872884299936, dG=-5484.812493437487, Tm=85.03347700825856 °C.
    // The tetraloop bonus (tetraloop.dh CGAAAG = -1100; tetraloop.ds = 0) is the difference from a
    // generic 4-nt loop — without the bundled table dH would be -39800.
    [Test]
    public void Ntthal_RecognisedTetraloop_MatchesPrimer3WithBundledBonus()
    {
        var r = PrimerDesigner.CalculateHairpinThermodynamicsNtthal("GGGGCGAAAGCCCC", Na);

        Assert.That(r, Is.Not.Null, "GGGGCGAAAGCCCC folds into a hairpin (5-bp stem, GAAA tetraloop).");
        var h = r!.Value;
        Assert.Multiple(() =>
        {
            Assert.That(h.DeltaH * 1000.0, Is.EqualTo(-40900.0).Within(ParityTol),
                "ΔH° = stem stacks + tstack2 terminal mismatch + tetraloop bonus (-1100 cal/mol); " +
                "matches primer3 calc_hairpin. Absent the bundled table this would be -39800.");
            Assert.That(h.DeltaS, Is.EqualTo(-114.1872884299936).Within(ParityTol),
                "ΔS° (incl. (N/2-1)·saltCorrection) matches primer3 calc_hairpin exactly.");
            Assert.That(h.DeltaG37 * 1000.0, Is.EqualTo(-5484.812493437487).Within(1e-3),
                "ΔG°37 = ΔH° - 310.15·ΔS° matches primer3 calc_hairpin.");
            Assert.That(h.TmCelsius, Is.EqualTo(85.03347700825856).Within(ParityTol),
                "Unimolecular Tm = ΔH°/ΔS° - 273.15 matches primer3 calc_hairpin.");
        });
    }

    // M2 — recognised tetraloop GGGGAC (closing G·C, loop GGGA). primer3-py:
    //   dH=-34000, dS=-94.1872884299836, dG=-4787.81249344059, Tm=87.8328944728006.
    [Test]
    public void Ntthal_SecondRecognisedTetraloop_MatchesPrimer3()
    {
        var r = PrimerDesigner.CalculateHairpinThermodynamicsNtthal("GGGGGGGACCCCC", Na);

        Assert.That(r, Is.Not.Null);
        var h = r!.Value;
        Assert.Multiple(() =>
        {
            Assert.That(h.DeltaH * 1000.0, Is.EqualTo(-34000.0).Within(ParityTol),
                "ΔH° matches primer3 calc_hairpin (tetraloop GGGGAC bonus applied).");
            Assert.That(h.DeltaS, Is.EqualTo(-94.1872884299836).Within(ParityTol), "ΔS° matches primer3.");
            Assert.That(h.TmCelsius, Is.EqualTo(87.8328944728006).Within(ParityTol), "Tm matches primer3.");
        });
    }

    #endregion

    #region CalculateHairpinThermodynamicsNtthal — special triloop (bundled bonus)

    // M3 — recognised triloop CGAAG (closing C·G, loop GAA). primer3-py calc_hairpin:
    //   dH=-27800, dS=-77.68485895331574, dG=-3706.040995629125, Tm=84.7060915802943.
    // The triloop bonus (triloop.dh CGAAG = -2000; triloop.ds = 0) plus the closing-A·T penalty
    // gate is the difference from a generic 3-nt loop.
    [Test]
    public void Ntthal_RecognisedTriloop_MatchesPrimer3WithBundledBonus()
    {
        var r = PrimerDesigner.CalculateHairpinThermodynamicsNtthal("GGGCGAAGCCC", Na);

        Assert.That(r, Is.Not.Null, "GGGCGAAGCCC folds into a hairpin (4-bp stem, GAA triloop).");
        var h = r!.Value;
        Assert.Multiple(() =>
        {
            Assert.That(h.DeltaH * 1000.0, Is.EqualTo(-27800.0).Within(ParityTol),
                "ΔH° = stem stacks + closing-A·T penalty + triloop bonus (-2000 cal/mol); matches primer3.");
            Assert.That(h.DeltaS, Is.EqualTo(-77.68485895331574).Within(ParityTol), "ΔS° matches primer3.");
            Assert.That(h.DeltaG37 * 1000.0, Is.EqualTo(-3706.040995629125).Within(1e-3), "ΔG°37 matches primer3.");
            Assert.That(h.TmCelsius, Is.EqualTo(84.7060915802943).Within(ParityTol), "Tm matches primer3.");
        });
    }

    // M4 — recognised triloop GGAAC (closing G·C, loop GAA). primer3-py:
    //   dH=-26000, dS=-73.18485895331571, dG=-3301.715995629134, Tm=82.11474153055735.
    [Test]
    public void Ntthal_SecondRecognisedTriloop_MatchesPrimer3()
    {
        var r = PrimerDesigner.CalculateHairpinThermodynamicsNtthal("GGGGGAACCCC", Na);

        Assert.That(r, Is.Not.Null);
        var h = r!.Value;
        Assert.Multiple(() =>
        {
            Assert.That(h.DeltaH * 1000.0, Is.EqualTo(-26000.0).Within(ParityTol), "ΔH° matches primer3 (triloop GGAAC bonus).");
            Assert.That(h.DeltaS, Is.EqualTo(-73.18485895331571).Within(ParityTol), "ΔS° matches primer3.");
            Assert.That(h.TmCelsius, Is.EqualTo(82.11474153055735).Within(ParityTol), "Tm matches primer3.");
        });
    }

    #endregion

    #region Non-special-loop regression (bundling must NOT change these)

    // S1 — 4-nt loop TTTT (closing C·G → key CTTTG, NOT in tetraloop.dh) is unchanged.
    //   primer3-py: dH=-32400, dS=-94.58485895332572, dG=-3064.5059956260297, Tm=69.39954078842845.
    [Test]
    public void Ntthal_NonSpecialTetraLoop_UnchangedByBundledTables()
    {
        var r = PrimerDesigner.CalculateHairpinThermodynamicsNtthal("GGGCTTTTGCCC", Na);

        Assert.That(r, Is.Not.Null);
        var h = r!.Value;
        Assert.Multiple(() =>
        {
            Assert.That(h.DeltaH * 1000.0, Is.EqualTo(-32400.0).Within(ParityTol),
                "A non-special 4-nt loop gets NO special bonus; ΔH° matches primer3 (generic loop + tstack2).");
            Assert.That(h.DeltaS, Is.EqualTo(-94.58485895332572).Within(ParityTol), "ΔS° unchanged, matches primer3.");
            Assert.That(h.TmCelsius, Is.EqualTo(69.39954078842845).Within(ParityTol), "Tm unchanged, matches primer3.");
        });
    }

    // S2 — 5-nt loop (never a special loop; only 3/4-nt loops have bonus tables) is unchanged.
    //   primer3-py: dH=-30100, dS=-87.74485895332572, dG=-2885.9319956260306, Tm=69.89004085311882.
    [Test]
    public void Ntthal_NonSpecialFiveLoop_UnchangedByBundledTables()
    {
        var r = PrimerDesigner.CalculateHairpinThermodynamicsNtthal("GGGCAAAAAGCCC", Na);

        Assert.That(r, Is.Not.Null);
        var h = r!.Value;
        Assert.Multiple(() =>
        {
            Assert.That(h.DeltaH * 1000.0, Is.EqualTo(-30100.0).Within(ParityTol), "5-nt loop: no special bonus; matches primer3.");
            Assert.That(h.DeltaS, Is.EqualTo(-87.74485895332572).Within(ParityTol), "ΔS° matches primer3.");
            Assert.That(h.TmCelsius, Is.EqualTo(69.89004085311882).Within(ParityTol), "Tm matches primer3.");
        });
    }

    #endregion

    #region Edge cases / invalid input

    // S3 — a homopolymer forms no hairpin (ntthal no_structure) → null, matching primer3 structure_found=False.
    [Test]
    public void Ntthal_Homopolymer_ReturnsNull()
    {
        Assert.That(PrimerDesigner.CalculateHairpinThermodynamicsNtthal("AAAAAAAAAAAA", Na), Is.Null,
            "poly-A admits no Watson-Crick stem → no hairpin (primer3 structure_found=False).");
    }

    // S4 — too short to close a 3-nt loop with a stem → null.
    [Test]
    public void Ntthal_TooShort_ReturnsNull()
    {
        Assert.That(PrimerDesigner.CalculateHairpinThermodynamicsNtthal("GCGC", Na), Is.Null,
            "GCGC cannot close a ≥3-nt hairpin loop → null (primer3 structure_found=False).");
    }

    // S5 — null / empty / non-ACGT input → null.
    [Test]
    public void Ntthal_InvalidInput_ReturnsNull()
    {
        Assert.Multiple(() =>
        {
            Assert.That(PrimerDesigner.CalculateHairpinThermodynamicsNtthal(null!, Na), Is.Null, "null → null.");
            Assert.That(PrimerDesigner.CalculateHairpinThermodynamicsNtthal("", Na), Is.Null, "empty → null.");
            Assert.That(PrimerDesigner.CalculateHairpinThermodynamicsNtthal("GGGCNAAGCCC", Na), Is.Null, "non-ACGT (N) → null.");
        });
    }

    #endregion

    #region Legacy Table-4 hairpin model unchanged (loopBonusDeltaG37 path still works)

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

    #endregion
}
