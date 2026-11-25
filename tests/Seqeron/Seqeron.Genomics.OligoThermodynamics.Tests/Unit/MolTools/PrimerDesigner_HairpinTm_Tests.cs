// SPDX-License-Identifier: GPL-2.0-or-later
// Tests of the port of Primer3 ntthal (thal.c) and of the Primer3-style primer / probe design; part of Seqeron.Genomics.OligoThermodynamics.Tests, licensed under the GNU General Public License, version 2 or (at your option) any later version (see LICENSE in the production project).
// PRIMER-TM-001 — DNA hairpin folding + secondary-structure (hairpin) Tm (opt-in)
// Evidence: docs/Evidence/PRIMER-TM-001-HAIRPIN-Evidence.md
// TestSpec: tests/TestSpecs/PRIMER-TM-001-HAIRPIN.md
// Sources:
//   SantaLucia J, Hicks D (2004). Annu Rev Biophys 33:415-440 — Table 1 NN stacks,
//     Table 4 hairpin-loop ΔG°37 by size, Eqs 7-11 (hairpin model + unimolecular Tm).
//   SantaLucia J (1998). PNAS 95(4):1460-65 — unified NN ΔH°/ΔS° (stem stacks reused).
//   Vallone PM, Benight AS (1999). Biochemistry — hairpin Tm concentration-independence.

namespace Seqeron.Genomics.Tests.Unit.MolTools;

/// <summary>
/// Tests for the opt-in DNA hairpin folder + unimolecular hairpin Tm added under
/// PRIMER-TM-001. Expected ΔH°/ΔS°/ΔG°37/Tm are hand-derived from the SantaLucia
/// (1998) Table 1 NN stem stacks and the SantaLucia &amp; Hicks (2004) Table 4
/// hairpin-loop-initiation increments, independent of the implementation.
/// </summary>
[TestFixture]
public class PrimerDesigner_HairpinTm_Tests
{

    /// <summary>
    /// Former discrepancies of the pre-PRIMER-HAIRPIN-001 port (calc_hairpin defaults mv 50 mM,
    /// dv 1.5 mM, dNTP 0.6 mM, 37 °C, max loop 30): the port evaluated thal.c RSH with the hairpin
    /// neighbours s[j−1] instead of s[j+1] and recomputed its −∞ running Tm, used G2 = −1 instead of
    /// thal.c's global G2 = 0 in calc_terminal_bp, a different max5 tie/NaN rule, a 1e−6 / ∞==∞
    /// equal(), and an incomplete traceback (old values: 69.31 / 44.72 / 38.27 / 61.80 / 139.59 °C, and a
    /// structure for the last oligo).
    /// </summary>
    [TestCase("GGGAGACAGTAGTCGCCCAT", 64.43690682436392, -3852.3691446110024, -47400.0, -140.40828907105916)]
    [TestCase("GGGCGGGGGGGC", 41.3455300128216, -272.2040000031011, -19700.0, -62.63999999999)]
    [TestCase("TGTTGAATATCAGCG", 29.597157947010032, 530.6133132321556, -21700.0, -71.67697344263149)]
    [TestCase("TACCTCATTCTGATGGGGGT", 54.81257620875647, -1705.422915689429, -31400.0, -95.74263125684531)]
    [TestCase("TTTGCCACTAATAATATGATCAACCGGAGGGTCTCCATT", 83.58095561848427, -457.020457846264, -3500.0, -9.811315628417656)]
    public void CalculateHairpinThermodynamicsNtthal_FormerDiscrepancies_MatchPrimer3Py(
        string seq, double tm, double dgCal, double dhCal, double ds)
    {
        var r = ThermodynamicPrimerDesigner.CalculateHairpinThermodynamicsNtthal(seq, 0.050, 0.0015, 0.0006);
        Assert.That(r, Is.Not.Null);
        var v = r!.Value;
        Assert.Multiple(() =>
        {
            Assert.That(v.TmCelsius, Is.EqualTo(tm).Within(NtTmTol), "Tm");
            Assert.That(v.DeltaG37 * 1000.0, Is.EqualTo(dgCal).Within(NtDgTol), "ΔG (cal/mol)");
            Assert.That(v.DeltaH * 1000.0, Is.EqualTo(dhCal).Within(NtDgTol), "ΔH (cal/mol)");
            Assert.That(v.DeltaS, Is.EqualTo(ds).Within(NtDgTol), "ΔS (cal/(K·mol))");
        });
    }

    [Test]
    public void CalculateHairpinThermodynamicsNtthal_FormerFalseStructure_IsNoStructure()
    {
        // primer3-py: structure_found=False; the old port reported Tm 32.64 °C.
        Assert.That(ThermodynamicPrimerDesigner.CalculateHairpinThermodynamicsNtthal(
            "AACCCGATAAAAAAGTTACACTCACTAAGAACAAGGGGGCTGCAAAAACTTTC", 0.050, 0.0015, 0.0006), Is.Null);
    }

    /// <summary>
    /// temp_c enters thal.c calc_terminal_bp (exterior candidate accepted only when ΔH − T·ΔS &lt; G2 = 0),
    /// so it can change the selected structure and the Tm, besides ΔG; max_loop bounds the bulge /
    /// internal-loop search. Values: primer3-py 2.3.1 calc_hairpin at the listed conditions.
    /// </summary>
    [TestCase("GCGGAGCTTCGTGGGAACCAGAGACA", 50, 1.5, 0.6, 70, 30, 72.01060100981272, -101.35704210848417, -17400.0, -50.411315628417654)]
    [TestCase("GGTTAGACCTGACCACGCGTGATCCGCGCCTC", 100, 3, 0.8, 70, 30, 78.61203029226152, -910.7507328347856, -37200.0, -105.75331274126539)]
    [TestCase("GGTTAGACCTGACCACGCGTGATCCGCGCCTC", 10, 0, 0, 25, 30, 48.0973646573546, -4069.483467983191, -56600.0, -176.18821577064168)]
    [TestCase("CCCTGAGTCCGAGGAGAGGGT", 50, 1.5, 0.6, 37, 2, 36.192654595551744, 76.46931323215904, -29300.0, -94.7169734426315)]
    [TestCase("AGGCACTTTGAATACGAATCATGCCACGTTACATGG", 200, 10, 2, 55, 2, 53.26062797107494, 172.652631096822, -32400.0, -99.26147381105234)]
    [TestCase("TGCGGGGGCATCTGCAGGCCGTCGTCCGTCCTACG", 50, 0, 0, 37, 0, 47.74656152259223, -853.973995626031, -25500.0, -79.46485895332572)]
    public void CalculateHairpinThermodynamicsNtthal_ConditionsTemperatureAndMaxLoop_MatchPrimer3Py(
        string seq, double mvMm, double dvMm, double dntpMm, double tempC, int maxLoop,
        double tm, double dgCal, double dhCal, double ds)
    {
        var r = ThermodynamicPrimerDesigner.CalculateHairpinThermodynamicsNtthal(
            seq, mvMm / 1000.0, dvMm / 1000.0, dntpMm / 1000.0, tempC, maxLoop);
        Assert.That(r, Is.Not.Null);
        var v = r!.Value;
        Assert.Multiple(() =>
        {
            Assert.That(v.TmCelsius, Is.EqualTo(tm).Within(NtTmTol), "Tm");
            Assert.That(v.DeltaG37 * 1000.0, Is.EqualTo(dgCal).Within(NtDgTol), "ΔG at temp_c (cal/mol)");
            Assert.That(v.DeltaH * 1000.0, Is.EqualTo(dhCal).Within(NtDgTol), "ΔH (cal/mol)");
            Assert.That(v.DeltaS, Is.EqualTo(ds).Within(NtDgTol), "ΔS (cal/(K·mol))");
        });
    }

    /// <summary>thal.c drawHairpin / primer3-py ascii_structure_lines ('/' = 5′ partner, '\' = 3′ partner).</summary>
    [TestCase("GGGAGACAGTAGTCGCCCAT", 50, 1.5, 0.6, 37, 30, "SEQ\t///-///----\\\\\\-\\\\\\--")]
    [TestCase("GGTTAGACCTGACCACGCGTGATCCGCGCCTC", 10, 0, 0, 25, 30, "SEQ\t///---\\\\\\------////-----\\\\\\\\----")]
    [TestCase("TTTGCCACTAATAATATGATCAACCGGAGGGTCTCCATT", 50, 1.5, 0.6, 37, 30, "SEQ\t---------------------/-//---\\\\-\\-------")]
    [TestCase("ggtcgaaaccatt", 50, 1.5, 0.6, 37, 30, "SEQ\t///----\\\\\\---")]
    [TestCase("aaaaaaaaaaaa", 50, 1.5, 0.6, 37, 30, null)]
    public void CalculateHairpinStructureNtthal_AsciiStructure_MatchesPrimer3Py(
        string seq, double mvMm, double dvMm, double dntpMm, double tempC, int maxLoop, string? seqLine)
    {
        var r = ThermodynamicPrimerDesigner.CalculateHairpinStructureNtthal(
            seq, mvMm / 1000.0, dvMm / 1000.0, dntpMm / 1000.0, tempC, maxLoop);
        if (seqLine is null)
        {
            Assert.That(r, Is.Null, "primer3-py: structure_found=False, ascii_structure_lines None.");
            return;
        }
        Assert.That(r, Is.Not.Null);
        Assert.That(r!.AsciiStructureLines, Is.EqualTo(new[] { seqLine, "STR\t" + seq.ToUpperInvariant() }));
        var t = ThermodynamicPrimerDesigner.CalculateHairpinThermodynamicsNtthal(
            seq, mvMm / 1000.0, dvMm / 1000.0, dntpMm / 1000.0, tempC, maxLoop);
        Assert.That(r.Thermodynamics, Is.EqualTo(t!.Value), "Structure overload carries the same thermodynamics.");
    }

    [Test]
    public void CalculateHairpinThermodynamicsNtthal_ThalLimits()
    {
        // thal.c CHECK_ERROR: both "sequences" of a hairpin are the oligo, so > THAL_MAX_ALIGN (60) nt
        // raises in primer3-py ("At least one sequence must be equal to or shorter than 60bp …").
        Assert.Multiple(() =>
        {
            Assert.That(() => ThermodynamicPrimerDesigner.CalculateHairpinThermodynamicsNtthal(new string('A', 61), 0.05, 0.0015, 0.0006),
                NUnit.Framework.Throws.ArgumentException);
            // 60-mer accepted: primer3-py calc_hairpin("GC"*30) → tm 106.34790650000872, poly-A 60-mer → no structure.
            Assert.That(ThermodynamicPrimerDesigner.CalculateHairpinThermodynamicsNtthal(new string('A', 60), 0.05, 0.0015, 0.0006), Is.Null);
            Assert.That(ThermodynamicPrimerDesigner.CalculateHairpinThermodynamicsNtthal(
                string.Concat(Enumerable.Repeat("GC", 30)), 0.05, 0.0015, 0.0006)!.Value.TmCelsius,
                Is.EqualTo(106.34790650000872).Within(NtTmTol));
            // primer3-py max_loop setter: ValueError above 30.
            Assert.That(() => ThermodynamicPrimerDesigner.CalculateHairpinThermodynamicsNtthal("GGGCTTTTGCCC", 0.05, 0.0015, 0.0006, 37.0, 31),
                NUnit.Framework.Throws.TypeOf<ArgumentOutOfRangeException>());
        });
    }




    // Reference values: primer3-py 2.3.1 primer3.calc_hairpin(seq, mv_conc, dv_conc, dntp_conc,
    // temp_c, max_loop, output_structure=True) — tm (°C), dg (cal/mol), dh (cal/mol), ds (cal/(K·mol))
    // printed with repr(). Engine source: primer3-py v2.3.1 primer3/src/libprimer3/thal.c (type 4).
    private const double NtTmTol = 1e-9, NtDgTol = 1e-6;
}
