// SPDX-License-Identifier: GPL-2.0-or-later
// Tests of the port of Primer3 ntthal (thal.c) and of the Primer3-style primer / probe design; part of Seqeron.Genomics.OligoThermodynamics.Tests, licensed under the GNU General Public License, version 2 or (at your option) any later version (see LICENSE in the production project).
namespace Seqeron.Genomics.Tests.Unit.MolTools;

/// <summary>
/// Canonical tests for PRIMER-STRUCT-001: Primer Structure Analysis.
/// Tests for secondary structure detection in PCR primers.
/// 
/// Methods under test:
/// - PrimerDesigner.FindLongestHomopolymer(seq)
/// - PrimerDesigner.FindLongestDinucleotideRepeat(seq)
/// - PrimerDesigner.HasHairpinPotential(seq, minStemLength)
/// - PrimerDesigner.HasPrimerDimer(primer1, primer2, minComplementarity)
/// - PrimerDesigner.Calculate3PrimeStability(seq)
/// 
/// Evidence sources:
/// - Wikipedia (Primer, Primer dimer, Stem-loop, Nucleic acid thermodynamics)
/// - Primer3 Manual (primer3.org)
/// - SantaLucia (1998) PNAS 95:1460-65
/// </summary>
[TestFixture]
public class PrimerDesigner_PrimerStructure_Tests
{

    /// <summary>
    /// M10b — <c>HasHairpinPotential</c> is a sequence-only library screen (default stem 4 unsourced, loop 3 =
    /// Primer3 thal.c min_hrpn_loop), not Primer3's hairpin screen: it flags both sequences although Primer3's
    /// ntthal hairpin (PRIMER_HAIRPIN_TH, limit 47 °C) accepts them.
    /// Source: primer3-py 2.3.1 <c>calc_hairpin</c> at 50 mM Na⁺ / 1.5 mM Mg²⁺ / 0.6 mM dNTP / 50 nM:
    /// AAAACCCTTTT → structure_found False (Tm 0); CAGTAAAACCCTTTTGCAGC → Tm 37.65 °C.
    /// </summary>
    [TestCase("AAAACCCTTTT", 0.0)]
    [TestCase("CAGTAAAACCCTTTTGCAGC", 37.65)]
    public void HasHairpinPotential_LibraryScreen_DiffersFromPrimer3NtthalHairpin(string sequence, double primer3HairpinTm)
    {
        var primer3 = ThermodynamicPrimerDesigner.CalculatePrimer3OligoStructure(sequence);
        Assert.Multiple(() =>
        {
            Assert.That(PrimerDesigner.HasHairpinPotential(sequence), Is.True, "library screen flags the exact 4-bp stem");
            Assert.That(primer3, Is.Not.Null);
            Assert.That(primer3!.Value.HairpinTh, Is.EqualTo(primer3HairpinTm).Within(0.005), "Primer3 ntthal hairpin Tm");
            Assert.That(primer3.Value.HairpinTh, Is.LessThanOrEqualTo(PrimerDesigner.Primer3MaxStructureTm),
                "Primer3 PRIMER_MAX_HAIRPIN_TH accepts it");
        });
    }

    /// <summary>
    /// Regression (DUP_MAP §10): identical poly-A primers cannot pair (A·A is not a base pair), so
    /// they are not a primer-dimer. The former implementation compared primer1's 3' window with
    /// revcomp(primer2) by complementarity (parallel pairing) and flagged A₈/A₈.
    /// Source: Primer3 compl_end = align(A₈, revcomp(A₈)) = 0 (dpal.c compiled from the primer3
    /// source; PRIMER_LEFT_0_SELF_END of AAAAAAAA = 0.0 in primer3-py check_primers);
    /// primer3-py calc_heterodimer(A₂₀, A₂₀) → structure_found = False.
    /// </summary>
    [Test]
    public void HasPrimerDimer_IdenticalPolyA_IsNotADimer()
    {
        Assert.Multiple(() =>
        {
            Assert.That(PrimerDesigner.CalculatePrimerDimerEndComplementarity("AAAAAAAA", "AAAAAAAA"), Is.EqualTo(0.0));
            Assert.That(PrimerDesigner.HasPrimerDimer("AAAAAAAA", "AAAAAAAA"), Is.False);
            Assert.That(ThermodynamicPrimerDesigner.CalculateDimerThermodynamicsNtthal(
                new string('A', 20), new string('A', 20)), Is.Null);
        });
    }

    [TestCaseSource(nameof(Primer3StructureCases))]
    public void Primer3ThermodynamicStructure_MatchesPrimer3CheckPrimers(string left, string right, double[] expected)
    {
        var l = ThermodynamicPrimerDesigner.CalculatePrimer3OligoStructure(left)!.Value;
        var r = ThermodynamicPrimerDesigner.CalculatePrimer3OligoStructure(right)!.Value;
        var pair = ThermodynamicPrimerDesigner.CalculatePrimer3PairComplementarity(left, right)!.Value;
        double[] actual = { l.SelfAnyTh, l.SelfEndTh, l.HairpinTh, r.SelfAnyTh, r.SelfEndTh, r.HairpinTh, pair.ComplAnyTh, pair.ComplEndTh };
        Assert.That(actual, Is.EqualTo(expected).Within(1e-9));
    }

    /// <summary>
    /// ntthal END1/END2 alignment types and the divalent/dNTP salt term.
    /// Source: primer3-py 2.3.1 calc_end_stability(a, b) (END1; END2(a, b) = END1(b, a)) and
    /// calc_heterodimer at mv 50, dv 1.5, dntp 0.6, dna 50 nM.
    /// </summary>
    [TestCase("TTCAGTCAGTCAGTGGCC", "ACTGACTGACTGAGGCC", ThermodynamicPrimerDesigner.NtthalAlignmentMode.End1, 10.942030550739162, -5.051857560300428)]
    [TestCase("TTCAGTCAGTCAGTGGCC", "ACTGACTGACTGAGGCC", ThermodynamicPrimerDesigner.NtthalAlignmentMode.End2, 13.481280192208146, -5.668962560300446)]
    [TestCase("TTCAGTCAGTCAGTGGCC", "ACTGACTGACTGAGGCC", ThermodynamicPrimerDesigner.NtthalAlignmentMode.Any, 41.167106798577606, -12.54145474705899)]
    [TestCase("TTCAGTCAGTCAGTAAAAGGGG", "ACTGACTGACTGACCCCTTTT", ThermodynamicPrimerDesigner.NtthalAlignmentMode.End1, 22.803405505733053, -8.035343602454173)]
    public void CalculateDimerThermodynamicsNtthal_AlignmentModes_MatchPrimer3Py(
        string a, string b, ThermodynamicPrimerDesigner.NtthalAlignmentMode mode, double tm, double dg)
    {
        var d = ThermodynamicPrimerDesigner.CalculateDimerThermodynamicsNtthal(a, b, mode, 0.050, 0.0015, 0.0006, 50e-9)!.Value;
        Assert.Multiple(() =>
        {
            Assert.That(d.TmCelsius, Is.EqualTo(tm).Within(1e-9));
            Assert.That(d.DeltaG37, Is.EqualTo(dg).Within(1e-9));
        });
    }

    /// <summary>
    /// Divalent cations enter ntthal hairpins only through saltCorrectS.
    /// Source: primer3-py 2.3.1 calc_hairpin('GGGGAAAACCCCATATGCAG') = 58.09337111624046 °C at
    /// dv 1.5/dntp 0.6 and 54.28830868675732 °C at dv = dntp = 0.
    /// </summary>
    [Test]
    public void CalculateHairpinThermodynamicsNtthal_Divalent_MatchesPrimer3Py()
    {
        Assert.Multiple(() =>
        {
            Assert.That(ThermodynamicPrimerDesigner.CalculateHairpinThermodynamicsNtthal("GGGGAAAACCCCATATGCAG", 0.050, 0.0015, 0.0006)!.Value.TmCelsius,
                Is.EqualTo(58.09337111624046).Within(1e-9));
            Assert.That(ThermodynamicPrimerDesigner.CalculateHairpinThermodynamicsNtthal("GGGGAAAACCCCATATGCAG", 0.050, 0.0, 0.0)!.Value.TmCelsius,
                Is.EqualTo(54.28830868675732).Within(1e-9));
        });
    }

    /// <summary>
    /// Poly-A cannot form any structure: every Primer3 structure value is 0 (align_thermod).
    /// Non-ACGT input has no ntthal structure (null).
    /// </summary>
    [Test]
    public void Primer3ThermodynamicStructure_NoStructureAndInvalidInput()
    {
        var a = ThermodynamicPrimerDesigner.CalculatePrimer3OligoStructure(new string('A', 20))!.Value;
        var pair = ThermodynamicPrimerDesigner.CalculatePrimer3PairComplementarity(new string('A', 20), new string('A', 20))!.Value;
        Assert.Multiple(() =>
        {
            Assert.That(a, Is.EqualTo(new PrimerDesigner.Primer3OligoStructure(0, 0, 0)));
            Assert.That(pair, Is.EqualTo(new PrimerDesigner.Primer3PairComplementarity(0, 0)));
            Assert.That(ThermodynamicPrimerDesigner.CalculatePrimer3OligoStructure("ACGTNACGT"), Is.Null);
            Assert.That(ThermodynamicPrimerDesigner.CalculatePrimer3PairComplementarity("ACGT", ""), Is.Null);
        });
    }

    /// <summary>
    /// EvaluatePrimer applies Primer3's default limits (47 °C): ACGTACGTACGTACGTACGT has
    /// self-any 60.49, self-end 60.49 and hairpin 71.76 °C (primer3-py check_primers), so it is
    /// rejected for all three; the heuristic screen reports only the stem-loop and no Tm values.
    /// </summary>
    [Test]
    public void EvaluatePrimer_ThermodynamicScreen_ReportsPrimer3Values()
    {
        var param = ThermodynamicPrimerDesigner.DefaultParameters with { MinGcContent = 0, MaxGcContent = 100, MinTm = 0, MaxTm = 100 };
        var thermo = ThermodynamicPrimerDesigner.EvaluatePrimer("ACGTACGTACGTACGTACGT", 0, true, param);
        var heur = ThermodynamicPrimerDesigner.EvaluatePrimer("ACGTACGTACGTACGTACGT", 0, true,
            param with { StructureScreen = PrimerStructureScreen.Heuristic });
        var relaxed = ThermodynamicPrimerDesigner.EvaluatePrimer("ACGTACGTACGTACGTACGT", 0, true, param with { MaxStructureTm = 100 });
        Assert.Multiple(() =>
        {
            Assert.That(thermo.SelfAnyTh, Is.EqualTo(60.49048111491186).Within(1e-9));
            Assert.That(thermo.SelfEndTh, Is.EqualTo(60.49048111491186).Within(1e-9));
            Assert.That(thermo.HairpinTh, Is.EqualTo(71.76211301142945).Within(1e-9));
            Assert.That(thermo.HasHairpin, Is.True);
            Assert.That(thermo.IsValid, Is.False);
            Assert.That(thermo.Issues.Count(i => i.Contains("Primer3 PRIMER_MAX_")), Is.EqualTo(3));
            Assert.That(heur.SelfAnyTh, Is.Null);
            Assert.That(heur.HasHairpin, Is.True);
            Assert.That(heur.Issues, Has.Some.EqualTo("Potential hairpin structure detected"));
            Assert.That(relaxed.HasHairpin, Is.False);
            Assert.That(relaxed.Issues.Any(i => i.Contains("Primer3 PRIMER_MAX_")), Is.False);
        });
    }

    [Test]
    public void EvaluatePrimer_Null_ThrowsArgumentNull()
    {
        // A7-4: was NullReferenceException in EvaluatePrimerCore (DesignPrimers & co. throw ArgumentNullException).
        var ex = Assert.Throws<ArgumentNullException>(() => ThermodynamicPrimerDesigner.EvaluatePrimer(null!, 0, true));
        Assert.That(ex!.ParamName, Is.EqualTo("sequence"));
    }

    [Test]
    public void GeneratePrimerCandidates_NullTemplate_ThrowsArgumentNull()
    {
        Assert.Throws<ArgumentNullException>(() => ThermodynamicPrimerDesigner.GeneratePrimerCandidates(null!, 0, 30));
    }

    [TestCase(-1, 30, "regionStart")]
    [TestCase(int.MinValue, 30, "regionStart")]
    [TestCase(int.MaxValue, 40, "regionStart")]
    [TestCase(0, 41, "regionEnd")]
    [TestCase(0, int.MaxValue, "regionEnd")]
    public void GeneratePrimerCandidates_RegionOutsideTemplate_ThrowsArgumentOutOfRangeEagerly(int start, int end, string param)
    {
        // Sweep: these crashed with a raw Substring AOORE (or start + MinLength overflow) during enumeration.
        var template = new DnaSequence("ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"); // 40 nt
        var ex = Assert.Throws<ArgumentOutOfRangeException>(() => ThermodynamicPrimerDesigner.GeneratePrimerCandidates(template, start, end));
        Assert.That(ex!.ParamName, Is.EqualTo(param));
    }

    [Test]
    public void GeneratePrimerCandidates_FullAndEmptyRegion_Unchanged()
    {
        var template = new DnaSequence("ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"); // 40 nt
        var all = ThermodynamicPrimerDesigner.GeneratePrimerCandidates(template, 0, 40).ToList();
        int min = ThermodynamicPrimerDesigner.DefaultParameters.MinLength, max = ThermodynamicPrimerDesigner.DefaultParameters.MaxLength;
        int expected = 0;
        for (int s = 0; s + min <= 40; s++)
            expected += Math.Min(max, 40 - s) - min + 1;
        Assert.That(all, Has.Count.EqualTo(expected));
        Assert.That(ThermodynamicPrimerDesigner.GeneratePrimerCandidates(template, 40, 40), Is.Empty);
        Assert.That(ThermodynamicPrimerDesigner.GeneratePrimerCandidates(template, 30, 10), Is.Empty);
    }

    [Test]
    public void CalculateTemplateMispriming_PositionPlusLengthOverflow_ThrowsArgumentOutOfRange()
    {
        // Sweep: position + length overflowed past the site check and crashed in Substring ('startIndex').
        var template = new DnaSequence("ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT");
        var ex = Assert.Throws<ArgumentOutOfRangeException>(
            () => ThermodynamicPrimerDesigner.CalculateTemplateMispriming(template, int.MaxValue, 20, true));
        Assert.That(ex!.ParamName, Is.EqualTo("position"));
    }




    /// <summary>
    /// Primer3's default (PRIMER_THERMODYNAMIC_OLIGO_ALIGNMENT=1) per-primer and pair structure values.
    /// Source: primer3-py 2.3.1 design_primers, PRIMER_TASK=check_primers, default conditions
    /// (50 mM monovalent, 1.5 mM Mg²⁺, 0.6 mM dNTP, 50 nM): PRIMER_{LEFT,RIGHT}_0_SELF_ANY_TH,
    /// _SELF_END_TH, _HAIRPIN_TH, PRIMER_PAIR_0_COMPL_ANY_TH, PRIMER_PAIR_0_COMPL_END_TH.
    /// </summary>
    private static readonly object[] Primer3StructureCases =
    {
        new object[] { "GGGGAAAACCCCATATGCAG", "CTGCATATGGGGTTTTCCCA",
            new[] { 18.30105741642683, 0.0, 58.09337111624046, 6.200192585916625, 6.200192585916625, 66.14453842604473, 57.51407534362289, 57.51407534362289 } },
        new object[] { "TTCAGTCAGTCAGTAAAAGGGG", "ACTGACTGACTGACCCCTTTT",
            new[] { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 42.22316686789662, 43.12474243211591 } },
        new object[] { "ACGTACGTACGTACGTACGT", "TGCAGCATGCATGCAATGCA",
            new[] { 60.49048111491186, 60.49048111491186, 71.76211301142945, 41.76336394643283, 41.76336394643283, 76.36745306969078, 0.0, 0.0 } },
        new object[] { "TGTCGGAAAACAGCCAGGCTAACG", "CGGTTGATCGTGAGTCAAATTCAGC",
            new[] { 15.644106552541643, 12.3535987107486, 36.400151579667295, 0.0, 0.35225774955779343, 37.336259822139596, 0.0, 0.0 } },
    };
}
