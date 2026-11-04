// PRIMER-HAIRPIN-001 / PRIMER-DIMER-001 — StructureModel opt-in (audit round 3, A3-13; F57 hairpin, F58 dimer)
// TestSpecs: tests/TestSpecs/PRIMER-HAIRPIN-001.md, tests/TestSpecs/PRIMER-DIMER-001.md
// Reference: primer3-py 2.3.1 calc_hairpin / calc_heterodimer / calc_homodimer (output_structure=True);
//   values below were captured from it (Tm °C; dh, dg cal/mol → /1000 kcal/mol; ds cal/(K·mol)); the
//   span / stem / loop fields are read off its ascii_structure_lines.

namespace Seqeron.Genomics.Tests.Unit.MolTools;

/// <summary>
/// <see cref="PrimerDesigner.StructureModel.Ntthal"/> makes <see cref="PrimerDesigner.FindMostStableHairpin(string, PrimerDesigner.StructureModel, double?, double?, double?, double?, int?)"/>,
/// <see cref="PrimerDesigner.CalculateHairpinMeltingTemperature(string, PrimerDesigner.StructureModel, double?, double?, double?, double?, int?)"/>
/// and <see cref="PrimerDesigner.FindMostStableDimer(string, string, PrimerDesigner.StructureModel, double?, double?, double?, double?, double?, int?)"/>
/// report the full thal.c structure (terminal mismatches, special loops, bulges, internal loops, dangles, Mg²⁺/dNTP);
/// <see cref="PrimerDesigner.StructureModel.SingleHelix"/> is the unchanged default.
/// </summary>
[TestFixture]
public class PrimerDesigner_StructureModel_Tests
{
    private const double Tol = 1e-9;
    private const PrimerDesigner.StructureModel Ntthal = PrimerDesigner.StructureModel.Ntthal;
    private const PrimerDesigner.StructureModel SingleHelix = PrimerDesigner.StructureModel.SingleHelix;

    #region Hairpin

    // calc_hairpin defaults (mv 50 mM, dv 1.5 mM, dNTP 0.6 mM, 37 °C, max_loop 30).
    //  GGGCTTTTGCCC: SEQ ////----\\\\ — 4-bp stem + TTTT loop, with the thal terminal mismatch the single-stem
    //    core lacks (its Tm is 68.64 °C).
    //  CCCTGAGTCCGAGGAGAGGGT: SEQ ////---//----\\-\\\\- — stem (0,19)(1,18)(2,17)(3,16) + internal loop/bulge +
    //    (7,14)(8,13) closing a 4-nt loop; ΔG at 37 °C is positive.
    //  TTTGCCACTAATAATATGATCAACCGGAGGGTCTCCATT (F13 former discrepancy): pairs (21,31)(23,29)(24,28) — a bulge
    //    and a 3-nt loop.
    [TestCase("GGGCTTTTGCCC", 0, 11, 4, 4, -32.4, -93.71131562841765, -3.335435457846266, 72.59266493570397)]
    [TestCase("CCCTGAGTCCGAGGAGAGGGT", 0, 19, 6, 4, -45.6, -148.34283687051104, 0.40853085538899177, 34.24603584502296)]
    [TestCase("TTTGCCACTAATAATATGATCAACCGGAGGGTCTCCATT", 21, 31, 3, 3, -3.5, -9.811315628417656, -0.457020457846264, 83.58095561848427)]
    public void FindMostStableHairpin_Ntthal_MatchesPrimer3PyCalcHairpin(
        string seq, int stemStart, int stemEnd, int stemLength, int loopSize,
        double dH, double dS, double dG, double tm)
    {
        var h = PrimerDesigner.FindMostStableHairpin(seq, Ntthal)!.Value;
        Assert.Multiple(() =>
        {
            Assert.That((h.StemStart, h.StemEnd, h.StemLength, h.LoopSize),
                Is.EqualTo((stemStart, stemEnd, stemLength, loopSize)), "structure from ascii_structure_lines");
            Assert.That(h.DeltaH, Is.EqualTo(dH).Within(Tol), "dh/1000");
            Assert.That(h.DeltaS, Is.EqualTo(dS).Within(Tol), "ds");
            Assert.That(h.DeltaG37, Is.EqualTo(dG).Within(Tol), "dg/1000");
            Assert.That(PrimerDesigner.CalculateHairpinMeltingTemperature(seq, Ntthal), Is.EqualTo(tm).Within(Tol), "tm");
        });
    }

    // calc_hairpin(..., mv 100, dv 3, dntp 0.8, temp_c 55, max_loop 10) and (mv 50, dv 0, dntp 0, temp_c 60,
    // max_loop 2) → structure_found False.
    [Test]
    public void FindMostStableHairpin_Ntthal_Conditions_MatchPrimer3Py()
    {
        Assert.Multiple(() =>
        {
            Assert.That(PrimerDesigner.FindMostStableHairpin("CCCTGAGTCCGAGGAGAGGGT", Ntthal, 0.05, 0.0, 0.0, 60.0, 2), Is.Null);
            Assert.That(double.IsNaN(PrimerDesigner.CalculateHairpinMeltingTemperature(
                "CCCTGAGTCCGAGGAGAGGGT", Ntthal, 0.05, 0.0, 0.0, 60.0, 2)), Is.True);
            // Delegation: the opt-in is the canonical ntthal engine, not a second implementation.
            var t = PrimerDesigner.CalculateHairpinThermodynamicsNtthal("GGGAGACAGTAGTCGCCCAT", 0.1, 0.003, 0.0008, 55.0, 10)!.Value;
            var h = PrimerDesigner.FindMostStableHairpin("gggagacagtagtcgcccat", Ntthal, 0.1, 0.003, 0.0008, 55.0, 10)!.Value;
            Assert.That((h.DeltaH, h.DeltaS, h.DeltaG37), Is.EqualTo((t.DeltaH, t.DeltaS, t.DeltaG37)));
            Assert.That(PrimerDesigner.CalculateHairpinMeltingTemperature("GGGAGACAGTAGTCGCCCAT", Ntthal, 0.1, 0.003, 0.0008, 55.0, 10),
                Is.EqualTo(t.TmCelsius));
        });
    }

    [Test]
    public void FindMostStableHairpin_SingleHelix_IsTheDefaultModel()
    {
        foreach (var seq in new[] { "GGGCTTTTGCCC", "GCGCGCGCGC", "GCCAAAGGC", "AAAAAAAAAAAA" })
        {
            Assert.That(PrimerDesigner.FindMostStableHairpin(seq, SingleHelix), Is.EqualTo(PrimerDesigner.FindMostStableHairpin(seq)));
            Assert.That(PrimerDesigner.CalculateHairpinMeltingTemperature(seq, SingleHelix),
                Is.EqualTo(PrimerDesigner.CalculateHairpinMeltingTemperature(seq)));
        }
    }

    [Test]
    public void FindMostStableHairpin_Model_GuardsAndLimits()
    {
        Assert.Multiple(() =>
        {
            Assert.Throws<ArgumentException>(() => PrimerDesigner.FindMostStableHairpin("GGGCTTTTGCCC", SingleHelix, sodiumMolar: 0.05),
                "SingleHelix has fixed conditions (1 M NaCl).");
            Assert.Throws<ArgumentException>(() => PrimerDesigner.CalculateHairpinMeltingTemperature("GGGCTTTTGCCC", SingleHelix, maxLoop: 30));
            Assert.Throws<ArgumentOutOfRangeException>(() => PrimerDesigner.FindMostStableHairpin("GGGCTTTTGCCC", (PrimerDesigner.StructureModel)7));
            Assert.Throws<ArgumentException>(() => PrimerDesigner.FindMostStableHairpin(new string('A', 30) + new string('T', 31), Ntthal),
                "thal.c THAL_MAX_ALIGN: primer3-py raises for a 61-nt hairpin.");
            Assert.Throws<ArgumentOutOfRangeException>(() => PrimerDesigner.FindMostStableHairpin("GGGCTTTTGCCC", Ntthal, maxLoop: 31));
            Assert.That(PrimerDesigner.FindMostStableHairpin("GGGCNNNNGCCC", Ntthal), Is.Null);
            Assert.That(PrimerDesigner.FindMostStableHairpin("", Ntthal), Is.Null);
            Assert.That(PrimerDesigner.FindMostStableHairpin("AAAAAAAAAAAA", Ntthal), Is.Null, "ntthal no_structure");
        });
    }

    #endregion

    #region Dimer

    // calc_heterodimer / calc_homodimer at the primer3-py defaults (mv 50, dv 1.5, dntp 0.6, dna 50 nM, 37 °C, 30):
    //  TCAGGTCAGCTAGGCATC / GATGCCTAGATGACCTGA: 17 bp with one internal G·A mismatch (the contiguous scorer can
    //    only see an 8- or 9-bp run).
    //  GCGCGCAAAA / AAAAGCGCGC: 6 bp, the AAAA tails are terminal overhangs (strand-2 5'-most paired base = 4).
    //  ACGTTGCAAGCTTGCAACGT self-dimer (calc_homodimer): 20 bp.
    [TestCase("TCAGGTCAGCTAGGCATC", "GATGCCTAGATGACCTGA", 0, 0, 17, -121.0, -338.75052502741124, -15.936524662748416)]
    [TestCase("GCGCGCAAAA", "AAAAGCGCGC", 0, 4, 6, -60.0, -163.12828907104912, -9.405761144614116)]
    [TestCase("ACGTTGCAAGCTTGCAACGT", "ACGTTGCAAGCTTGCAACGT", 0, 0, 20, -160.8, -443.34749847006265, -23.29577334951007)]
    public void FindMostStableDimer_Ntthal_MatchesPrimer3Py(
        string a, string b, int start1, int start2, int basePairs, double dH, double dS, double dG)
    {
        var d = PrimerDesigner.FindMostStableDimer(a, b, Ntthal)!.Value;
        Assert.Multiple(() =>
        {
            Assert.That((d.Strand1Start, d.Strand2Start, d.BasePairs), Is.EqualTo((start1, start2, basePairs)),
                "from ascii_structure_lines");
            Assert.That(d.DeltaH, Is.EqualTo(dH).Within(Tol), "dh/1000");
            Assert.That(d.DeltaS, Is.EqualTo(dS).Within(Tol), "ds");
            Assert.That(d.DeltaG37, Is.EqualTo(dG).Within(Tol), "dg/1000");
        });
    }

    // calc_heterodimer(..., mv 100, dv 3, dntp 0.8, dna_conc 250, temp_c 55, max_loop 10): ds -335.6376679534587,
    // dg -10860.499261072546, tm 55.11775596157162 (the Tm is CalculateDimerThermodynamicsNtthal's).
    [Test]
    public void FindMostStableDimer_Ntthal_Conditions_MatchPrimer3Py()
    {
        var d = PrimerDesigner.FindMostStableDimer("TCAGGTCAGCTAGGCATC", "GATGCCTAGATGACCTGA", Ntthal,
            0.1, 0.003, 0.0008, 250e-9, 55.0, 10)!.Value;
        var t = PrimerDesigner.CalculateDimerThermodynamicsNtthal("TCAGGTCAGCTAGGCATC", "GATGCCTAGATGACCTGA",
            PrimerDesigner.NtthalAlignmentMode.Any, 0.1, 0.003, 0.0008, 250e-9, 55.0, 10)!.Value;
        Assert.Multiple(() =>
        {
            Assert.That((d.Strand1Start, d.Strand2Start, d.BasePairs), Is.EqualTo((0, 0, 17)));
            Assert.That(d.DeltaH, Is.EqualTo(-121.0).Within(Tol));
            Assert.That(d.DeltaS, Is.EqualTo(-335.6376679534587).Within(Tol));
            Assert.That(d.DeltaG37, Is.EqualTo(-10.860499261072546).Within(Tol));
            Assert.That(t.TmCelsius, Is.EqualTo(55.11775596157162).Within(Tol));
            Assert.That((d.DeltaH, d.DeltaS, d.DeltaG37, d.BasePairs), Is.EqualTo((t.DeltaH, t.DeltaS, t.DeltaG37, t.BasePairs)),
                "delegates to the canonical ntthal dimer engine");
        });
    }

    [Test]
    public void FindMostStableDimer_SingleHelix_IsTheDefaultModel()
    {
        Assert.Multiple(() =>
        {
            Assert.That(PrimerDesigner.FindMostStableDimer("GCGCGCGC", "GCGCGCGC", SingleHelix),
                Is.EqualTo(PrimerDesigner.FindMostStableDimer("GCGCGCGC", "GCGCGCGC")));
            Assert.That(PrimerDesigner.FindMostStableDimer("TCAGGTCAGCTAGGCATC", "GATGCCTAGATGACCTGA", SingleHelix, 0.2, strandConcentrationMolar: 1e-6),
                Is.EqualTo(PrimerDesigner.FindMostStableDimer("TCAGGTCAGCTAGGCATC", "GATGCCTAGATGACCTGA", 0.2, 1e-6)));
            Assert.That(PrimerDesigner.FindMostStableDimer("TCAGGTCAGCTAGGCATC", "GATGCCTAGATGACCTGA", SingleHelix)!.Value.BasePairs,
                Is.LessThan(17), "the contiguous scorer cannot bridge the internal mismatch");
        });
    }

    [Test]
    public void FindMostStableDimer_Model_GuardsAndLimits()
    {
        Assert.Multiple(() =>
        {
            Assert.Throws<ArgumentException>(() => PrimerDesigner.FindMostStableDimer("GCGCGCGC", "GCGCGCGC", SingleHelix, divalentMolar: 0.0015),
                "SingleHelix has no divalent-cation model.");
            Assert.Throws<ArgumentException>(() => PrimerDesigner.FindMostStableDimer("GCGCGCGC", "GCGCGCGC", SingleHelix, temperatureCelsius: 37));
            Assert.Throws<ArgumentOutOfRangeException>(() => PrimerDesigner.FindMostStableDimer("GCGC", "GCGC", (PrimerDesigner.StructureModel)(-1)));
            string longA = new string('A', 61), longT = new string('T', 61);
            Assert.Throws<ArgumentException>(() => PrimerDesigner.FindMostStableDimer(longA, longT, Ntthal),
                "thal.c THAL_MAX_ALIGN: both strands > 60 nt.");
            Assert.Throws<ArgumentOutOfRangeException>(() => PrimerDesigner.FindMostStableDimer("GCGC", "GCGC", Ntthal, maxLoop: -1));
            Assert.That(PrimerDesigner.FindMostStableDimer("GCGCN", "GCGC", Ntthal), Is.Null);
            Assert.That(PrimerDesigner.FindMostStableDimer(null!, "GCGC", Ntthal), Is.Null);
        });
    }

    #endregion
}
