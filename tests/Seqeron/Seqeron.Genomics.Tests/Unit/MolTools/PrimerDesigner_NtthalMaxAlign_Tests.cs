// PRIMER-DIMER-001 / PRIMER-HAIRPIN-001 / PROBE-DESIGN-001 — opt-in THAL_MAX_ALIGN override (audit round 3, A3-9, F55)
// Evidence: docs/Evidence/PRIMER-TM-001-HAIRPIN-Evidence.md, docs/Validation/review-2026-09/B07.md (F55)
// TestSpec: tests/TestSpecs/PROBE-DESIGN-001.md
// Source: primer3-py 2.3.1 primer3/src/libprimer3/thal.h (#ifndef THAL_MAX_ALIGN / #define THAL_MAX_ALIGN 60) and
//   thal.c thal_check_errors (the only use of THAL_MAX_ALIGN; DP tables are allocated from the actual lengths).
// Reference: thal.c + thal_parameters.c of primer3-py 2.3.1 compiled with -DTHAL_MAX_ALIGN=10000, driven through
//   thal() with THL_GENERAL / print_output 0 (= primer3-py calc_hairpin / calc_homodimer / calc_end_stability) and
//   the bundled primer3_config tables. The 60-default build raises "At least one sequence must be equal to or
//   shorter than 60bp …" for every case below. 1680 random cases (560 oligos of 61–120 nt × hairpin/ANY/END1,
//   four condition sets) agree to |ΔTm| = 0, |ΔG| ≤ 1.5e-11 cal/mol.

namespace Seqeron.Genomics.Tests.Unit.MolTools;

[TestFixture]
public class PrimerDesigner_NtthalMaxAlign_Tests
{
    private const double TmTol = 1e-9;
    private const double KcalTol = 1e-9;

    private const string Random61 = "AGACTTTCAAAGATATGCTGGGTAGAGGTCGAGGTTATTATTTGTTACCAATTCTCATTGT";
    private const string Random120 =
        "GTTTCGGAACTTGCGTTTTAGGTATGTCTTAGTGACTCTAAATACCAAGGCAGTCCTCGATCCGTTCCTAATAAGGAATGGTGATTCCCTGTCATACCAATCTACCCCCTGTTATGCGCG";
    // 20-bp stem + GAAA loop with random flanks (80 nt).
    private const string StemLoop80 =
        "AGCGCAGCGGCAGATTTTGTCGTTAGACCAATGTCGAAAGACATTGGTCTAACGACAAACAAGCAGGAGGCGGAATGTAA";
    // 12-bp stem + TTTT loop with random flanks (82 nt).
    private const string StemLoop82 =
        "CTTAGGTGGATAGGGAGTGAGCAACAAACGACAGAAGGTATGTTTTCATACCTTCTGTGATCGTTTCTCCCATGCCAAGTTG";
    private static readonly string AcgtRepeat64 = string.Concat(Enumerable.Repeat("ACGT", 16));

    // mv / dv / dntp mM, dna nM, temp_c, max_loop, sequence; expected thal.c (THAL_MAX_ALIGN 10000) Tm °C,
    // ΔG cal/mol at temp_c, ΔH cal/mol, ΔS cal/(K·mol).
    private static IEnumerable<TestCaseData> HairpinCases()
    {
        yield return new TestCaseData(50.0, 0.0, 0.0, 37.0, 30, Random61,
            33.980529935122263, 268.39250656251897, -27300.0, -88.887288429993617).SetName("Hairpin_Random61_ProbeConditions");
        yield return new TestCaseData(50.0, 1.5, 0.6, 37.0, 30, Random120,
            49.326287932858634, -5454.5445784378098, -142700.0, -442.51315628425664).SetName("Hairpin_Random120_PrimerConditions");
        yield return new TestCaseData(100.0, 3.0, 0.8, 55.0, 20, StemLoop80,
            82.407778125350774, -13420.305967659981, -174100.0, -489.65318918890756).SetName("Hairpin_StemLoop80_NonDefault");
        yield return new TestCaseData(50.0, 0.0, 0.0, 37.0, 30, StemLoop82,
            68.102992186803021, -8631.2894759260962, -94700.0, -277.50672424334647).SetName("Hairpin_StemLoop82_ProbeConditions");
    }

    private static IEnumerable<TestCaseData> DimerCases()
    {
        yield return new TestCaseData(PrimerDesigner.NtthalAlignmentMode.Any, 50.0, 0.0, 0.0, 50.0, 37.0, 30, Random61,
            -0.96817310228368569, -4017.7154846864287, -51600.0, -153.41700633665508).SetName("SelfAny_Random61");
        yield return new TestCaseData(PrimerDesigner.NtthalAlignmentMode.End1, 50.0, 0.0, 0.0, 50.0, 37.0, 30, Random61,
            -79.597058062067646, -854.3414978145629, -17200.0, -52.702429476657869).SetName("SelfEnd_Random61");
        yield return new TestCaseData(PrimerDesigner.NtthalAlignmentMode.Any, 50.0, 1.5, 0.6, 50.0, 37.0, 30, Random120,
            29.127204238873617, -4011.6553857972031, -276600.0, -878.8919703827271).SetName("SelfAny_Random120");
        yield return new TestCaseData(PrimerDesigner.NtthalAlignmentMode.End1, 50.0, 1.5, 0.6, 50.0, 37.0, 30, Random120,
            -1.0124740646909913, -5488.7636867709516, -41000.0, -114.49697344262148).SetName("SelfEnd_Random120");
        yield return new TestCaseData(PrimerDesigner.NtthalAlignmentMode.Any, 100.0, 3.0, 0.8, 250.0, 55.0, 20, StemLoop80,
            63.277378655777511, -20471.609859263641, -392400.0, -1133.4096911191114).SetName("SelfAny_StemLoop80_NonDefault");
        yield return new TestCaseData(PrimerDesigner.NtthalAlignmentMode.End1, 100.0, 3.0, 0.8, 250.0, 55.0, 20, StemLoop80,
            60.562591020483637, -16726.186576626846, -354500.0, -1029.3274826249374).SetName("SelfEnd_StemLoop80_NonDefault");
        yield return new TestCaseData(PrimerDesigner.NtthalAlignmentMode.Any, 50.0, 0.0, 0.0, 50.0, 37.0, 30, AcgtRepeat64,
            72.769499884640766, -66385.62436212186, -541800.0, -1532.8530570300763).SetName("SelfAny_AcgtRepeat64");
    }

    [TestCaseSource(nameof(HairpinCases))]
    public void CalculateHairpinThermodynamicsNtthal_RaisedMaxAlign_MatchesThalCompiledWithLargerMaxAlign(
        double mv, double dv, double dntp, double tempC, int maxLoop, string seq,
        double tm, double dgCal, double dhCal, double ds)
    {
        var h = PrimerDesigner.CalculateHairpinThermodynamicsNtthal(seq, mv / 1000, dv / 1000, dntp / 1000, tempC, maxLoop, 120)!.Value;
        var st = PrimerDesigner.CalculateHairpinStructureNtthal(seq, mv / 1000, dv / 1000, dntp / 1000, tempC, maxLoop, 120)!;
        Assert.Multiple(() =>
        {
            Assert.That(h.TmCelsius, Is.EqualTo(tm).Within(TmTol));
            Assert.That(h.DeltaG37, Is.EqualTo(dgCal / 1000).Within(KcalTol));
            Assert.That(h.DeltaH, Is.EqualTo(dhCal / 1000).Within(KcalTol));
            Assert.That(h.DeltaS, Is.EqualTo(ds).Within(1e-9));
            Assert.That(st.Thermodynamics, Is.EqualTo(h));
            Assert.That(st.AsciiStructureLines, Has.Count.EqualTo(2));
            Assert.That(st.AsciiStructureLines[1], Is.EqualTo("STR\t" + seq));
            // Primer3's own limit is unchanged: the default overloads refuse > 60 nt as primer3-py does.
            Assert.That(() => PrimerDesigner.CalculateHairpinThermodynamicsNtthal(seq, mv / 1000, dv / 1000, dntp / 1000, tempC, maxLoop),
                NUnit.Framework.Throws.ArgumentException.With.Message.Contains("shorter than 60bp"));
        });
    }

    [TestCaseSource(nameof(DimerCases))]
    public void CalculateDimerThermodynamicsNtthal_RaisedMaxAlign_MatchesThalCompiledWithLargerMaxAlign(
        PrimerDesigner.NtthalAlignmentMode mode, double mv, double dv, double dntp, double dnaNm, double tempC, int maxLoop,
        string seq, double tm, double dgCal, double dhCal, double ds)
    {
        var d = PrimerDesigner.CalculateDimerThermodynamicsNtthal(
            seq, seq, mode, mv / 1000, dv / 1000, dntp / 1000, dnaNm * 1e-9, tempC, maxLoop, 120)!.Value;
        var st = PrimerDesigner.CalculateDimerStructureNtthal(
            seq, seq, mode, mv / 1000, dv / 1000, dntp / 1000, dnaNm * 1e-9, tempC, maxLoop, 120)!;
        Assert.Multiple(() =>
        {
            Assert.That(d.TmCelsius, Is.EqualTo(tm).Within(TmTol));
            Assert.That(d.DeltaG37, Is.EqualTo(dgCal / 1000).Within(KcalTol));
            Assert.That(d.DeltaH, Is.EqualTo(dhCal / 1000).Within(KcalTol));
            Assert.That(d.DeltaS, Is.EqualTo(ds).Within(1e-9));
            Assert.That(st.Thermodynamics, Is.EqualTo(d));
            Assert.That(st.AsciiStructureLines, Has.Count.EqualTo(4));
            Assert.That(() => PrimerDesigner.CalculateDimerThermodynamicsNtthal(
                    seq, seq, mode, mv / 1000, dv / 1000, dntp / 1000, dnaNm * 1e-9, tempC, maxLoop),
                NUnit.Framework.Throws.ArgumentException.With.Message.Contains("shorter than 60bp"));
        });
    }

    [Test]
    public void CalculatePrimer3OligoStructure_RaisedMaxAlign_ScreensLongOligo()
    {
        // align_thermod: Tm < 0 → 0 (Random61 self_any −0.968, self_end −79.6).
        var r61 = PrimerDesigner.CalculatePrimer3OligoStructure(Random61, 50, 0, 0, 50, 61)!.Value;
        var r82 = PrimerDesigner.CalculatePrimer3OligoStructure(StemLoop82, 50, 0, 0, 50, 100)!.Value;
        Assert.Multiple(() =>
        {
            Assert.That(r61.SelfAnyTh, Is.EqualTo(0.0));
            Assert.That(r61.SelfEndTh, Is.EqualTo(0.0));
            Assert.That(r61.HairpinTh, Is.EqualTo(33.980529935122263).Within(TmTol));
            Assert.That(r82.SelfAnyTh, Is.EqualTo(45.91673334496744).Within(TmTol));
            Assert.That(r82.SelfEndTh, Is.EqualTo(43.010478511327449).Within(TmTol));
            Assert.That(r82.HairpinTh, Is.EqualTo(68.102992186803021).Within(TmTol));
            Assert.That(() => PrimerDesigner.CalculatePrimer3OligoStructure(Random61, 50, 0, 0, 50),
                NUnit.Framework.Throws.ArgumentException);
        });
    }

    [Test]
    public void RaisedMaxAlign_LeavesOligosWithinSixtyUnchanged()
    {
        const string s = "GCGGAGCTTCGTGGGAACCAGAGACA"; // primer3-py calc_hairpin (50/1.5/0.6): Tm 62.713271406705132
        Assert.Multiple(() =>
        {
            Assert.That(PrimerDesigner.CalculateHairpinThermodynamicsNtthal(s, 0.05, 0.0015, 0.0006, 37, 30, 10000),
                Is.EqualTo(PrimerDesigner.CalculateHairpinThermodynamicsNtthal(s, 0.05, 0.0015, 0.0006, 37, 30)));
            Assert.That(PrimerDesigner.CalculateHairpinThermodynamicsNtthal(s, 0.05, 0.0015, 0.0006, 37, 30, 10000)!.Value.TmCelsius,
                Is.EqualTo(62.713271406705132).Within(TmTol));
            Assert.That(PrimerDesigner.CalculatePrimer3OligoStructure(s, 50, 1.5, 0.6, 50, 10000),
                Is.EqualTo(PrimerDesigner.CalculatePrimer3OligoStructure(s, 50, 1.5, 0.6, 50)));
        });
    }

    [Test]
    public void RaisedMaxAlign_Limits()
    {
        Assert.Multiple(() =>
        {
            Assert.That(PrimerDesigner.NtthalMaxAlignLength, Is.EqualTo(60));
            Assert.That(PrimerDesigner.NtthalMaxSequenceLength, Is.EqualTo(10000));
            Assert.That(() => PrimerDesigner.CalculateHairpinThermodynamicsNtthal(Random61, 0.05, 0, 0, 37, 30, 59),
                NUnit.Framework.Throws.TypeOf<ArgumentOutOfRangeException>());
            Assert.That(() => PrimerDesigner.CalculateHairpinThermodynamicsNtthal(Random61, 0.05, 0, 0, 37, 30, 10001),
                NUnit.Framework.Throws.TypeOf<ArgumentOutOfRangeException>());
            // The override is the THAL_MAX_ALIGN value: a 61-mer needs ≥ 61 (thal.c message carries the value).
            Assert.That(() => PrimerDesigner.CalculateHairpinThermodynamicsNtthal(StemLoop80, 0.05, 0, 0, 37, 30, 79),
                NUnit.Framework.Throws.ArgumentException.With.Message.Contains("shorter than 79bp"));
            Assert.That(PrimerDesigner.CalculateHairpinThermodynamicsNtthal(StemLoop80, 0.05, 0, 0, 37, 30, 80), Is.Not.Null);
            Assert.That(() => PrimerDesigner.CalculateDimerThermodynamicsNtthal(StemLoop80, StemLoop80,
                    PrimerDesigner.NtthalAlignmentMode.Any, 0.05, 0, 0, 5e-8, 37, 30, 79),
                NUnit.Framework.Throws.ArgumentException);
            Assert.That(() => PrimerDesigner.CalculatePrimer3OligoStructure(Random61, 50, 0, 0, 50, 0),
                NUnit.Framework.Throws.TypeOf<ArgumentOutOfRangeException>());
        });
    }

    [Test]
    public void ValidateProbe_ThermodynamicScreenMaxLength_ScreensLongProbeWithNtthal()
    {
        var conditions = new ProbeDesigner.ProbeParameters(20, 120, -1000, 1000, 0, 1, 100, true, 0.3);
        var fallback = ProbeDesigner.ValidateProbe(AcgtRepeat64, Enumerable.Empty<string>(), conditions: conditions);
        var thermo = ProbeDesigner.ValidateProbe(AcgtRepeat64, Enumerable.Empty<string>(),
            conditions: conditions with { ThermodynamicScreenMaxLength = 64 });
        var stemLoop = ProbeDesigner.ValidateProbe(StemLoop82, Enumerable.Empty<string>(),
            conditions: conditions with { ThermodynamicScreenMaxLength = 100 });
        Assert.Multiple(() =>
        {
            // Default THAL_MAX_ALIGN 60: unchanged fallback (Primer3 alignment-mode self_any 64.00).
            Assert.That(fallback.ThermodynamicScreen, Is.False);
            Assert.That(fallback.HairpinTm, Is.Null);
            Assert.That(fallback.Issues, Has.Some.EqualTo("Self-complementarity: Primer3 self_any 64.00 exceeds 12.00"));
            // Opt-in: ntthal at 50 nM / 50 mM / 0 / 0 (thal.c THAL_MAX_ALIGN 10000: hairpin 77.098, ANY = END1 72.769).
            Assert.That(thermo.ThermodynamicScreen, Is.True);
            Assert.That(thermo.HairpinTm, Is.EqualTo(77.098245155727511).Within(TmTol));
            Assert.That(thermo.SelfDimerTm, Is.EqualTo(72.769499884640766).Within(TmTol));
            Assert.That(thermo.SelfEndDimerTm, Is.EqualTo(72.769499884640766).Within(TmTol));
            Assert.That(thermo.Issues, Has.Some.EqualTo("Self-complementarity: ntthal self-dimer Tm 72.8°C exceeds 47°C"));
            Assert.That(thermo.Issues, Has.Some.EqualTo("Potential secondary structure formation: ntthal hairpin Tm 77.1°C exceeds 47°C"));
            Assert.That(thermo.SelfAny, Is.EqualTo(64.0), "alignment-mode values are still reported");
            // 82-nt stem-loop: hairpin 68.10 > 47 flagged; self-dimer 45.92 / 43.01 ≤ 47 not flagged.
            Assert.That(stemLoop.ThermodynamicScreen, Is.True);
            Assert.That(stemLoop.HasSecondaryStructure, Is.True);
            Assert.That(stemLoop.HairpinTm, Is.EqualTo(68.102992186803021).Within(TmTol));
            Assert.That(stemLoop.SelfDimerTm, Is.EqualTo(45.91673334496744).Within(TmTol));
            Assert.That(stemLoop.Issues, Has.None.Contain("Self-complementarity"));
            Assert.That(() => ProbeDesigner.ValidateProbe(StemLoop82, Enumerable.Empty<string>(),
                    conditions: conditions with { ThermodynamicScreenMaxLength = 59 }),
                NUnit.Framework.Throws.TypeOf<ArgumentOutOfRangeException>());
        });
    }

    [Test]
    public void DesignProbes_ThermodynamicScreenMaxLength_ScreensLongProbeWithNtthal()
    {
        var p = new ProbeDesigner.ProbeParameters(82, 82, -1000, 1000, 0, 1, 100, true, 1.0);
        var byDefault = ProbeDesigner.DesignProbes(StemLoop82, p).Single();
        var optIn = ProbeDesigner.DesignProbes(StemLoop82, p with { ThermodynamicScreenMaxLength = 82 }).Single();
        Assert.Multiple(() =>
        {
            Assert.That(byDefault.Warnings, Has.None.Contain("(ntthal)"));
            Assert.That(optIn.Warnings, Has.Some.EqualTo("Hairpin Tm 68.1°C exceeds 47°C (ntthal)"));
            Assert.That(optIn.Warnings, Has.None.Contain("Self-dimer Tm"));
            Assert.That(() => ProbeDesigner.DesignProbes(StemLoop82, p with { ThermodynamicScreenMaxLength = 10001 }),
                NUnit.Framework.Throws.TypeOf<ArgumentOutOfRangeException>());
        });
    }

    // --- A3-27 (F56): AssessCrossHybridization site duplex Tm and DesignMolecularBeacon stem-loop Tm ---

    // 75-nt probe = Random120[0..75); non-target A = flanks + the probe with substitutions at 10, 35, 60 (site 75 nt);
    // non-target B = flanks + probe[20..60) (site 40 nt).
    private static readonly string Probe75 = Random120.Substring(0, 75);
    private const string CrossNonTargetA =
        "GATTACAGATTACAGTTTCGGAACATGCGTTTTAGGTATGTCTTAGTGAATCTAAATACCAAGGCAGTCCTCGAACCGTTCCTAATAAGCCTTGGAACC";
    private static readonly string CrossNonTargetB = string.Concat("TTTTTGGGGG", Probe75.AsSpan(20, 40), "AAAAACCCCC");
    private static ProbeDesigner.ProbeParameters ProbeConditions => new(20, 120, -1000, 1000, 0, 1, 100, true, 0.3);

    [Test]
    public void AssessCrossHybridization_MaxAlignOptIn_SiteDuplexTm_MatchesThal()
    {
        var nonTargets = new[] { CrossNonTargetA, CrossNonTargetB };
        var byDefault = ProbeDesigner.AssessCrossHybridization(Probe75, nonTargets);
        var optIn = ProbeDesigner.AssessCrossHybridization(Probe75, nonTargets,
            conditions: ProbeConditions with { ThermodynamicScreenMaxLength = 120 });
        Assert.Multiple(() =>
        {
            Assert.That((optIn[0].SiteStart, optIn[0].SiteEnd), Is.EqualTo((14, 88)));
            // Probe 75 nt and site 75 nt both > 60: thal.c THAL_MAX_ALIGN 60 refuses (primer3-py raises) → null.
            Assert.That(byDefault[0].DuplexTm, Is.Null);
            // thal.c compiled with -DTHAL_MAX_ALIGN=10000, THAL_ANY probe vs revcomp(site), 50 mM / 0 / 0 / 50 nM, 37 °C.
            Assert.That(optIn[0].DuplexTm, Is.EqualTo(66.457703046655695).Within(TmTol));
            // Sites ≤ 60 nt (7, 40, 13 nt): computed by default too (thal_check_errors needs one strand ≤ 60) and
            // equal to primer3-py calc_heterodimer(probe, revcomp(site)); the opt-in leaves them unchanged.
            foreach (var r in new[] { byDefault, optIn })
            {
                Assert.That(r[1].DuplexTm, Is.EqualTo(19.05924515571178).Within(TmTol));
                Assert.That(r[2].DuplexTm, Is.EqualTo(63.99555300270714).Within(TmTol));
                Assert.That(r[3].DuplexTm, Is.EqualTo(-33.229679404433625).Within(TmTol));
            }
            Assert.That(() => ProbeDesigner.AssessCrossHybridization(Probe75, nonTargets,
                    conditions: ProbeConditions with { ThermodynamicScreenMaxLength = 59 }),
                NUnit.Framework.Throws.TypeOf<ArgumentOutOfRangeException>());
        });
    }

    [Test]
    public void ValidateProbe_ThermodynamicScreenMaxLength_AppliesToNonTargetDuplexTm()
    {
        var v = ProbeDesigner.ValidateProbe(Probe75, Enumerable.Empty<string>(),
            conditions: ProbeConditions with { ThermodynamicScreenMaxLength = 120 },
            nonTargetSequences: new[] { CrossNonTargetA }, maxDuplexTm: 60);
        Assert.Multiple(() =>
        {
            // thal.c (THAL_MAX_ALIGN 10000) at 50/0/0/50 nM: hairpin 36.318725375690178, ANY 23.325175622075108,
            // END1 17.579447109879538; duplex with the site 66.457703046655695.
            Assert.That(v.ThermodynamicScreen, Is.True);
            Assert.That(v.HairpinTm, Is.EqualTo(36.318725375690178).Within(TmTol));
            Assert.That(v.SelfDimerTm, Is.EqualTo(23.325175622075108).Within(TmTol));
            Assert.That(v.SelfEndDimerTm, Is.EqualTo(17.579447109879538).Within(TmTol));
            Assert.That(v.CrossHybridization[0].DuplexTm, Is.EqualTo(66.457703046655695).Within(TmTol));
            Assert.That(v.CrossHybridization[0].ExceedsDuplexTmThreshold, Is.True);
            Assert.That(v.Issues, Has.Some.EqualTo(
                "Cross-hybridization risk with non-target 0: identity 96%, longest contiguous match 24 nt (Kane 2000), site duplex Tm 66.5°C > 60°C"));
        });
    }

    [Test]
    public void DesignMolecularBeacon_MaxAlignOptIn_ReportsStemLoopTmOfLongBeacon()
    {
        // Loop 60 + stem 7 → 74-nt beacon GGGCCCC + Random120[0..60) + GGGGCCC. thal.c -DTHAL_MAX_ALIGN=10000 hairpin at
        // 50 mM / 0 / 0, 37 °C: 45.449128344157543 °C (the 60-nt build and primer3-py calc_hairpin raise).
        const string beacon = "GGGCCCCGTTTCGGAACTTGCGTTTTAGGTATGTCTTAGTGACTCTAAATACCAAGGCAGTCCTCGAGGGGCCC";
        var byDefault = ProbeDesigner.DesignMolecularBeacon(Random120, 60, 7)!.Value;
        var optIn = ProbeDesigner.DesignMolecularBeacon(Random120, 60, 7, maxAlignLength: 100)!.Value;
        var detect = ProbeDesigner.DesignMolecularBeacon(Random120, 60, 7, detectionTemperatureCelsius: 60, maxAlignLength: 100)!.Value;
        var hp = PrimerDesigner.CalculateHairpinThermodynamicsNtthal(beacon, 0.05, 0, 0, 37, 30, 100)!.Value;
        Assert.Multiple(() =>
        {
            Assert.That(optIn.Sequence, Is.EqualTo(beacon));
            Assert.That(hp.TmCelsius, Is.EqualTo(45.449128344157543).Within(TmTol));
            Assert.That(byDefault.Warnings, Is.EqualTo(new[] { "Stem: 7bp, Loop: 60bp" }));
            Assert.That(optIn.Warnings, Is.EqualTo(new[] { "Stem: 7bp, Loop: 60bp", "Stem-loop (hairpin) Tm 45.4°C (ntthal)" }));
            Assert.That(detect.Warnings, Has.Some.EqualTo(
                "Stem-loop Tm 45.4°C is less than 7 °C above the detection temperature 60°C"));
            Assert.That((optIn.Start, optIn.End, optIn.Tm), Is.EqualTo((byDefault.Start, byDefault.End, byDefault.Tm)));
            Assert.That(() => ProbeDesigner.DesignMolecularBeacon(Random120, 60, 7, maxAlignLength: 59),
                NUnit.Framework.Throws.TypeOf<ArgumentOutOfRangeException>());
            Assert.That(() => ProbeDesigner.DesignMolecularBeacon(Random120, 60, 7, maxAlignLength: 10001),
                NUnit.Framework.Throws.TypeOf<ArgumentOutOfRangeException>());
        });
    }
}
