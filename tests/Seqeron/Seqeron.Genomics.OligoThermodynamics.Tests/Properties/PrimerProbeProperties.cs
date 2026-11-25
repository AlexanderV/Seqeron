// SPDX-License-Identifier: GPL-2.0-or-later
// Tests of the port of Primer3 ntthal (thal.c) and of the Primer3-style primer / probe design; part of Seqeron.Genomics.OligoThermodynamics.Tests, licensed under the GNU General Public License, version 2 or (at your option) any later version (see LICENSE in the production project).
using System.Text;
using FsCheck;
using FsCheck.Fluent;

namespace Seqeron.Genomics.Tests.Properties;

/// <summary>
/// Property-based tests for primer/probe design: melting temperature, hairpin, dimer, probe validation.
///
/// Test Units: PRIMER-TM-001, PRIMER-DESIGN-001, PRIMER-STRUCT-001, PROBE-DESIGN-001, PROBE-VALID-001
/// </summary>
[TestFixture]
[Category("Property")]
[Category("MolTools")]
public class PrimerProbeProperties
{


    /// <summary>
    /// Generates a randomized parameter window (the design "config knobs") together with a
    /// template long enough to host both flanking search regions and a target between them.
    /// The length/GC/Tm windows are deliberately varied so the invariants are exercised across
    /// configurations rather than a single hard-coded setting (cf. ADVANCED_TESTING_CHECKLIST
    /// "multiple configuration knobs": length range + GC range + Tm range).
    /// </summary>
    private static Arbitrary<(string template, int targetStart, int targetEnd, PrimerParameters param)>
        DesignScenarioArbitrary() =>
        (from minLen in Gen.Choose(15, 20)
         from extra in Gen.Choose(0, 8)
         let maxLen = minLen + extra
         from minGc in Gen.Choose(20, 45)
         from gcSpan in Gen.Choose(15, 50)
         let maxGc = Math.Min(100, minGc + gcSpan)
         from minTm in Gen.Choose(40, 58)
         from tmSpan in Gen.Choose(8, 30)
         let maxTm = minTm + tmSpan
         // Template must hold: forward flank (>= maxLen) + target (>= 1) + reverse flank (>= maxLen).
         from fwdFlank in Gen.Choose(maxLen, maxLen + 40)
         from targetSpan in Gen.Choose(1, 30)
         from revFlank in Gen.Choose(maxLen, maxLen + 40)
         let total = fwdFlank + targetSpan + revFlank
         from chars in Gen.Elements('A', 'C', 'G', 'T').ArrayOf(total)
         let param = new PrimerParameters(
             MinLength: minLen, MaxLength: maxLen, OptimalLength: (minLen + maxLen) / 2,
             MinGcContent: minGc, MaxGcContent: maxGc,
             MinTm: minTm, MaxTm: maxTm, OptimalTm: (minTm + maxTm) / 2.0,
             MaxHomopolymer: 4, MaxDinucleotideRepeats: 4,
             Avoid3PrimeGC: false, Check3PrimeStability: true)
         select (new string(chars), fwdFlank, fwdFlank + targetSpan, param)).ToArbitrary();

    /// <summary>
    /// Generates a candidate-enumeration scenario: a template, a region inside it, an orientation,
    /// and a randomized parameter window. Used to test <see cref="ThermodynamicPrimerDesigner.GeneratePrimerCandidates"/>.
    /// </summary>
    private static Arbitrary<(string template, int regionStart, int regionEnd, bool forward, PrimerParameters param)>
        CandidateScenarioArbitrary() =>
        (from minLen in Gen.Choose(15, 20)
         from extra in Gen.Choose(0, 8)
         let maxLen = minLen + extra
         from minGc in Gen.Choose(20, 45)
         from gcSpan in Gen.Choose(15, 50)
         let maxGc = Math.Min(100, minGc + gcSpan)
         from minTm in Gen.Choose(40, 58)
         from tmSpan in Gen.Choose(8, 30)
         let maxTm = minTm + tmSpan
         from regionLen in Gen.Choose(maxLen, maxLen + 60)
         from forward in Gen.Elements(true, false)
         from chars in Gen.Elements('A', 'C', 'G', 'T').ArrayOf(regionLen)
         let param = new PrimerParameters(
             MinLength: minLen, MaxLength: maxLen, OptimalLength: (minLen + maxLen) / 2,
             MinGcContent: minGc, MaxGcContent: maxGc,
             MinTm: minTm, MaxTm: maxTm, OptimalTm: (minTm + maxTm) / 2.0,
             MaxHomopolymer: 4, MaxDinucleotideRepeats: 4,
             Avoid3PrimeGC: false, Check3PrimeStability: true)
         select (new string(chars), 0, regionLen, forward, param)).ToArbitrary();

    /// <summary>
    /// Independent re-implementation of the pair-validity predicate (INV-02), derived from the
    /// doc contract rather than reused production code: a pair is valid iff |Tm_f − Tm_r| ≤ 5
    /// AND the two primers do not form a primer-dimer.
    /// </summary>
    // Primer3 compares unrounded Tm values (PrimerCandidate.MeltingTemperature is rounded to 0.1 °C).
    /// <summary>
    /// Independent pair-validity predicate of Primer_Design.md §2.4 INV-02 at the library defaults
    /// (<see cref="PrimerPairOptions.Default"/>): product size inside a PRIMER_PRODUCT_SIZE_RANGE range
    /// (default 100-300), |Tm_f − Tm_r| ≤ PRIMER_PAIR_MAX_DIFF_TM (5 °C, unrounded Primer3 Tm) and no
    /// primer-dimer by Primer3's default thermodynamic screen (PRIMER_PAIR_MAX_COMPL_ANY_TH /
    /// PRIMER_PAIR_MAX_COMPL_END_TH: compl_any_th / compl_end_th ≤ 47 °C at the default primer conditions).
    /// </summary>
    private static bool ExpectedPairValid(PrimerCandidate fwd, PrimerCandidate rev)
    {
        var opt = PrimerPairOptions.Default;
        int product = rev.Position + rev.Sequence.Length - fwd.Position;
        bool sizeOk = opt.ProductSizeRanges.Any(r => product >= r.Min && product <= r.Max);
        bool tmOk = Math.Abs(ThermodynamicPrimerDesigner.CalculateMeltingTemperaturePrimer3(fwd.Sequence)
                             - ThermodynamicPrimerDesigner.CalculateMeltingTemperaturePrimer3(rev.Sequence)) <= opt.MaxTmDifference;
        var compl = ThermodynamicPrimerDesigner.CalculatePrimer3PairComplementarity(fwd.Sequence, rev.Sequence);
        bool dimerOk = compl is { } c
                       && c.ComplAnyTh <= PrimerDesigner.Primer3MaxStructureTm
                       && c.ComplEndTh <= PrimerDesigner.Primer3MaxStructureTm;
        return sizeOk && tmOk && dimerOk;
    }

    /// <summary>
    /// INV-R-len (negation): an invalid-length primer is never marked valid by
    /// <see cref="ThermodynamicPrimerDesigner.EvaluatePrimer"/>. Builds primers one base shorter than MinLength
    /// and one longer than MaxLength and asserts IsValid is false in both cases.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property EvaluatePrimer_LengthOutsideWindow_IsNeverValid()
    {
        var arb = (from minLen in Gen.Choose(16, 22)
                   from extra in Gen.Choose(2, 8)
                   let maxLen = minLen + extra
                   from chars in Gen.Elements('A', 'C', 'G', 'T').ArrayOf(maxLen + 6)
                   select (new string(chars), minLen, maxLen)).ToArbitrary();

        return Prop.ForAll(arb, t =>
        {
            var (seq, minLen, maxLen) = t;
            var param = new PrimerParameters(minLen, maxLen, (minLen + maxLen) / 2,
                0, 100, 0, 200, 60, 99, 99, false, false);

            string tooShort = seq.Substring(0, minLen - 1);
            string tooLong = seq.Substring(0, maxLen + 1);

            bool shortInvalid = !ThermodynamicPrimerDesigner.EvaluatePrimer(tooShort, 0, true, param).IsValid;
            bool longInvalid = !ThermodynamicPrimerDesigner.EvaluatePrimer(tooLong, 0, true, param).IsValid;

            return (shortInvalid && longInvalid)
                .Label($"len-window violation accepted: short={shortInvalid}, long={longInvalid} for [{minLen},{maxLen}]");
        });
    }



    /// <summary>
    /// INV-03 (product size): when <see cref="ThermodynamicPrimerDesigner.DesignPrimers"/> returns both primers,
    /// ProductSize equals reverse.Position + reverse.Sequence.Length − forward.Position, recomputed
    /// independently from the returned candidates. Source: Primer_Design.md §2.4 INV-03.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property DesignPrimers_ProductSize_MatchesIndependentFormula()
    {
        return Prop.ForAll(DesignScenarioArbitrary(), s =>
        {
            var dna = new DnaSequence(s.template);
            var result = ThermodynamicPrimerDesigner.DesignPrimers(dna, s.targetStart, s.targetEnd, s.param);

            if (result.Forward is null || result.Reverse is null)
                return true.ToProperty(); // INV-03 only applies when both primers exist.

            int expected = result.Reverse.Position + result.Reverse.Sequence.Length - result.Forward.Position;
            return (result.ProductSize == expected)
                .Label($"ProductSize {result.ProductSize} != {expected}");
        });
    }

    /// <summary>
    /// INV-02 (pair validity): result.IsValid ⟺ (product size in range AND |Tm_f − Tm_r| ≤ 5 AND no
    /// thermodynamic primer-dimer), with the predicate recomputed independently (ExpectedPairValid)
    /// from the returned candidates (on failure DesignPrimers returns the individually best primers,
    /// which then must violate a pair constraint — otherwise the exhaustive search would have paired them).
    /// Source: Primer_Design.md §2.4 INV-02.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property DesignPrimers_PairValidity_MatchesIndependentPredicate()
    {
        return Prop.ForAll(DesignScenarioArbitrary(), s =>
        {
            var dna = new DnaSequence(s.template);
            var result = ThermodynamicPrimerDesigner.DesignPrimers(dna, s.targetStart, s.targetEnd, s.param);

            if (result.Forward is null || result.Reverse is null)
                return !result.IsValid ? true.ToProperty()  // INV-01: missing side ⇒ invalid.
                    : false.ToProperty().Label("missing primer but IsValid==true");

            bool expected = ExpectedPairValid(result.Forward, result.Reverse);
            return (result.IsValid == expected)
                .Label($"IsValid {result.IsValid} != predicate {expected} " +
                       $"(Tm {result.Forward.MeltingTemperature}/{result.Reverse.MeltingTemperature})");
        });
    }

    /// <summary>
    /// INV-R/P/R on the chosen pair: whenever DesignPrimers selects primers, EACH selected primer
    /// satisfies the configured length, GC, and Tm windows (selection only ever draws from valid
    /// candidates). Ties the three filter invariants to the public design entry point.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property DesignPrimers_SelectedPrimers_SatisfyAllWindows()
    {
        return Prop.ForAll(DesignScenarioArbitrary(), s =>
        {
            var dna = new DnaSequence(s.template);
            var result = ThermodynamicPrimerDesigner.DesignPrimers(dna, s.targetStart, s.targetEnd, s.param);

            foreach (var c in new[] { result.Forward, result.Reverse })
            {
                if (c is null) continue;
                bool ok = c.Length >= s.param.MinLength && c.Length <= s.param.MaxLength
                          && c.GcContent >= s.param.MinGcContent && c.GcContent <= s.param.MaxGcContent
                          && c.MeltingTemperature >= s.param.MinTm && c.MeltingTemperature <= s.param.MaxTm;
                if (!ok)
                    return false.Label(
                        $"selected primer len={c.Length} gc={c.GcContent} tm={c.MeltingTemperature} " +
                        $"violates windows len[{s.param.MinLength},{s.param.MaxLength}] " +
                        $"gc[{s.param.MinGcContent},{s.param.MaxGcContent}] tm[{s.param.MinTm},{s.param.MaxTm}]");
            }

            return true.ToProperty();
        });
    }



    /// <summary>
    /// INV-D (DesignPrimers): identical inputs yield identical pair results — same validity,
    /// product size, message, and selected primer sequences/positions.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property DesignPrimers_IsDeterministic()
    {
        return Prop.ForAll(DesignScenarioArbitrary(), s =>
        {
            var dna = new DnaSequence(s.template);
            var a = ThermodynamicPrimerDesigner.DesignPrimers(dna, s.targetStart, s.targetEnd, s.param);
            var b = ThermodynamicPrimerDesigner.DesignPrimers(dna, s.targetStart, s.targetEnd, s.param);

            bool same = a.IsValid == b.IsValid
                        && a.ProductSize == b.ProductSize
                        && a.Message == b.Message
                        && a.Forward?.Sequence == b.Forward?.Sequence
                        && a.Forward?.Position == b.Forward?.Position
                        && a.Reverse?.Sequence == b.Reverse?.Sequence
                        && a.Reverse?.Position == b.Reverse?.Position;

            return same.Label("DesignPrimers not deterministic for identical inputs");
        });
    }

    /// <summary>
    /// INV-D (GeneratePrimerCandidates): the candidate list is identical across repeated calls —
    /// same count, sequences, positions, validity flags, and scores.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property GeneratePrimerCandidates_IsDeterministic()
    {
        return Prop.ForAll(CandidateScenarioArbitrary(), s =>
        {
            var dna = new DnaSequence(s.template);
            var a = ThermodynamicPrimerDesigner.GeneratePrimerCandidates(dna, s.regionStart, s.regionEnd, s.forward, s.param).ToList();
            var b = ThermodynamicPrimerDesigner.GeneratePrimerCandidates(dna, s.regionStart, s.regionEnd, s.forward, s.param).ToList();

            if (a.Count != b.Count)
                return false.Label($"candidate count differs: {a.Count} vs {b.Count}");

            for (int i = 0; i < a.Count; i++)
            {
                if (a[i].Sequence != b[i].Sequence || a[i].Position != b[i].Position
                    || a[i].IsValid != b[i].IsValid || a[i].Score != b[i].Score)
                    return false.Label($"candidate {i} differs between runs");
            }

            return true.ToProperty();
        });
    }



    /// <summary>
    /// Edge (§6.1): an invalid target region throws ArgumentException — targetStart &lt; 0,
    /// targetEnd ≥ template.Length, or targetStart ≥ targetEnd. Source: explicit guard in DesignPrimers.
    /// </summary>
    [TestCase(-1, 50, TestName = "targetStart < 0 throws")]
    [TestCase(0, 500, TestName = "targetEnd >= length throws")]
    [TestCase(60, 60, TestName = "targetStart == targetEnd throws")]
    [TestCase(80, 40, TestName = "targetStart > targetEnd throws")]
    [Category("Property")]
    public void DesignPrimers_InvalidTargetRegion_ThrowsArgumentException(int start, int end)
    {
        var dna = new DnaSequence(string.Concat(Enumerable.Repeat("ACGT", 75))); // length 300
        Assert.That(() => ThermodynamicPrimerDesigner.DesignPrimers(dna, start, end),
            NUnit.Framework.Throws.ArgumentException);
    }

    /// <summary>
    /// Evidence anchor (GcContent + Tm): "CGGTTCACTACGTCCGTTCTGG" has 13 G/C of 22 bases ⇒
    /// GC% = 1300/22 = 59.09 → 59.1 (rounded). EvaluatePrimer's Tm is the Primer3-default Tm
    /// (SantaLucia 1998 NN, 50 mM Na+, 1.5 mM Mg2+, 0.6 mM dNTP, 50 nM): primer3-py 2.3.1
    /// calc_tm = 63.1228 °C → 63.1, just above PRIMER_MAX_TM = 63 ⇒ the primer is rejected.
    /// </summary>
    [Test]
    [Category("Property")]
    public void EvaluatePrimer_KnownMer_HasComputedGcAndTm()
    {
        var c = ThermodynamicPrimerDesigner.EvaluatePrimer("CGGTTCACTACGTCCGTTCTGG", 0, true);
        Assert.Multiple(() =>
        {
            Assert.That(c.Length, Is.EqualTo(22));
            Assert.That(c.GcContent, Is.EqualTo(59.1).Within(0.05),
                "GC% = 100×13/22 = 59.09 → 59.1");
            Assert.That(c.MeltingTemperature, Is.EqualTo(63.1).Within(0.05),
                "primer3.calc_tm = 63.1228 → 63.1");
            Assert.That(c.IsValid, Is.False, "Tm 63.12 > PRIMER_MAX_TM 63");
            Assert.That(c.Issues, Has.Some.StartsWith("Tm"));
        });
    }

    /// <summary>
    /// Evidence anchor (valid pair + product size): 44 bp template = forward 20-mer
    /// "GATTCGAAGGGGATAGCGCA" (0-19) + target "AAAA" (20-23) + reverse complement of the reverse
    /// 20-mer "ATGGGCGTGGGCATAATACC" (24-43). primer3-py 2.3.1 design_primers (sizes 18/20/25,
    /// Tm 57/60/63, GC 40-60, poly-X 4, pair ΔTm ≤ 5, SEQUENCE_TARGET 20,4,
    /// PRIMER_PRODUCT_SIZE_RANGE 40-60 — the 44-bp template is below Primer3's default 100-300) returns
    /// PRIMER_LEFT_0 = [0,20] (Tm 59.9665, penalty 0.033504), PRIMER_RIGHT_0 = [43,20]
    /// (Tm 60.2512, penalty 0.251165), PRIMER_PAIR_0_PENALTY = 0.284669, product 44.
    /// ProductSize = reverse.Position(24) + reverse.Length(20) − forward.Position(0) = 44.
    /// </summary>
    [Test]
    [Category("Property")]
    public void DesignPrimers_KnownTemplate_ProductSizeEqualsSpan()
    {
        const string forwardPrimer = "GATTCGAAGGGGATAGCGCA";
        const string reversePrimer = "ATGGGCGTGGGCATAATACC";
        string reverseFlank = DnaSequence.GetReverseComplementString(reversePrimer);
        string template = forwardPrimer + "AAAA" + reverseFlank; // length 44

        var dna = new DnaSequence(template);
        var pairOptions = new PrimerPairOptions { ProductSizeRanges = [new ProductSizeRange(40, 60)] };
        var result = ThermodynamicPrimerDesigner.DesignPrimers(dna, targetStart: 20, targetEnd: 24, pairOptions: pairOptions);

        Assert.That(result.Forward, Is.Not.Null, "expected a forward primer on this template");
        Assert.That(result.Reverse, Is.Not.Null, "expected a reverse primer on this template");

        int expected = result.Reverse!.Position + result.Reverse.Sequence.Length - result.Forward!.Position;
        Assert.Multiple(() =>
        {
            Assert.That(result.IsValid, Is.True, "engineered pair must be compatible");
            Assert.That(result.Forward!.Sequence, Is.EqualTo(forwardPrimer));
            Assert.That(result.Reverse!.Sequence, Is.EqualTo(reversePrimer));
            Assert.That(result.Reverse.Position, Is.EqualTo(24));
            Assert.That(result.Forward.Penalty, Is.EqualTo(0.033504).Within(1e-6));
            Assert.That(result.Reverse.Penalty, Is.EqualTo(0.251165).Within(1e-6));
            Assert.That(result.Forward.Penalty + result.Reverse.Penalty, Is.EqualTo(0.284669).Within(1e-6));
            Assert.That(result.ProductSize, Is.EqualTo(expected));
            Assert.That(result.ProductSize, Is.EqualTo(44), "ProductSize = 24 + 20 − 0 = 44");
        });
    }

    /// <summary>
    /// Generates a randomized probe-design scenario: a target long enough to host probe windows,
    /// plus a randomized parameter window (the design "config knobs" — length, GC and Tm bands).
    /// The windows are deliberately varied so the invariants are exercised across configurations
    /// rather than a single hard-coded preset (cf. ADVANCED_TESTING_CHECKLIST multi-knob probe
    /// mention: length range + GC range + Tm range). Lengths are kept short (probe windows
    /// ~12-26 bp) so hundreds of cases run quickly while still spanning the length filter.
    /// </summary>
    private static Arbitrary<(string target, ThermodynamicProbeDesigner.ProbeParameters param)>
        ProbeScenarioArbitrary() =>
        (from minLen in Gen.Choose(12, 18)
         from extra in Gen.Choose(0, 8)
         let maxLen = minLen + extra
         from minGcPct in Gen.Choose(25, 50)
         from gcSpanPct in Gen.Choose(10, 40)
         let maxGcPct = Math.Min(95, minGcPct + gcSpanPct)
         from minTm in Gen.Choose(30, 60)
         from tmSpan in Gen.Choose(5, 35)
         let maxTm = minTm + tmSpan
         // Target must be at least maxLen long so windows exist; add slack for multiple windows.
         from slack in Gen.Choose(0, 30)
         from chars in Gen.Elements('A', 'C', 'G', 'T').ArrayOf(maxLen + slack)
         let param = new ThermodynamicProbeDesigner.ProbeParameters(
             MinLength: minLen, MaxLength: maxLen,
             MinTm: minTm, MaxTm: maxTm,
             MinGc: minGcPct / 100.0, MaxGc: maxGcPct / 100.0,
             MaxHomopolymer: 6,
             AvoidSecondaryStructure: false,
             MaxSelfComplementarity: 0.5)
         select (new string(chars), param)).ToArbitrary();

    /// <summary>
    /// Generates a target plus two NESTED Tm windows (B ⊆ A) over the same length/GC parameters,
    /// for the subset-monotonicity property. Both windows share identical length and GC bands so
    /// the only difference is the Tm acceptance band — window B is strictly no wider than window A.
    /// </summary>
    private static Arbitrary<(string target, ThermodynamicProbeDesigner.ProbeParameters wide, ThermodynamicProbeDesigner.ProbeParameters narrow)>
        NestedTmScenarioArbitrary() =>
        (from minLen in Gen.Choose(12, 18)
         from extra in Gen.Choose(0, 6)
         let maxLen = minLen + extra
         // Wide Tm window [aMin, aMax]; narrow window [bMin, bMax] nested inside it.
         from aMin in Gen.Choose(20, 50)
         from aSpan in Gen.Choose(30, 70)
         let aMax = aMin + aSpan
         from bLeftPad in Gen.Choose(0, aSpan / 2)
         from bRightPad in Gen.Choose(0, aSpan / 2)
         let bMin = aMin + bLeftPad
         let bMax = aMax - bRightPad
         from slack in Gen.Choose(10, 40)
         from chars in Gen.Elements('A', 'C', 'G', 'T').ArrayOf(maxLen + slack)
         let wide = new ThermodynamicProbeDesigner.ProbeParameters(
             MinLength: minLen, MaxLength: maxLen,
             MinTm: aMin, MaxTm: aMax,
             MinGc: 0.30, MaxGc: 0.70,
             MaxHomopolymer: 6, AvoidSecondaryStructure: false, MaxSelfComplementarity: 0.5)
         let narrow = wide with { MinTm = bMin, MaxTm = bMax }
         select (new string(chars), wide, narrow)).ToArbitrary();



    /// <summary>
    /// INV-R-len (strict, by construction): EVERY returned <see cref="ProbeDesigner.Probe"/> has
    /// Sequence.Length within [MinLength, MaxLength], with coordinates consistent
    /// (End − Start + 1 == Sequence.Length, Start ≥ 0, End &lt; target.Length). This is a HARD
    /// guarantee: the candidate loop only emits windows whose length lies in the configured range
    /// (ProbeDesigner.cs §DesignProbesOptimized). Exercised across randomized targets AND randomized
    /// length windows so the bound is a true config-driven invariant. Source: §4.1 step 3; the
    /// length filter is structural, never a soft penalty.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property DesignProbes_Length_StrictlyInConfiguredWindow()
    {
        return Prop.ForAll(ProbeScenarioArbitrary(), s =>
        {
            var (target, param) = s;
            var probes = ThermodynamicProbeDesigner.DesignProbes(target, param, maxProbes: 100000).ToList();

            foreach (var p in probes)
            {
                if (p.Sequence.Length < param.MinLength || p.Sequence.Length > param.MaxLength)
                    return false.Label(
                        $"length {p.Sequence.Length} outside [{param.MinLength},{param.MaxLength}]");
                if (p.End - p.Start + 1 != p.Sequence.Length)
                    return false.Label($"coords inconsistent: {p.Start}..{p.End} vs len {p.Sequence.Length}");
                if (p.Start < 0 || p.End >= target.Length)
                    return false.Label($"coords out of target bounds: {p.Start}..{p.End}, target len {target.Length}");
            }

            return true.ToProperty();
        });
    }



    /// <summary>
    /// INV-P-GC (RELAXED band, NOT strict): GC is only a −0.3 soft penalty in
    /// <c>EvaluateProbeWithGc</c>; the sole HARD GC bound is the early-rejection filter
    /// <c>gc &lt; MinGc − 0.1 || gc &gt; MaxGc + 0.1</c>. Therefore the guaranteed invariant is
    /// GcContent ∈ [MinGc − 0.1, MaxGc + 0.1] (a fraction) — NOT the strict [MinGc, MaxGc].
    /// The production GcContent is also recomputed independently (ExpectedGcFraction) to confirm
    /// the stored fraction is itself correct (INV-02: gcCount/length). Source: ProbeDesigner.cs
    /// §DesignProbesOptimized early-rejection; §2.2 penalty table; doc §6 "outside range" is a warning.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property DesignProbes_GcContent_WithinRelaxedBand_AndMatchesOracle()
    {
        return Prop.ForAll(ProbeScenarioArbitrary(), s =>
        {
            var (target, param) = s;
            var probes = ThermodynamicProbeDesigner.DesignProbes(target, param, maxProbes: 100000).ToList();

            foreach (var p in probes)
            {
                double oracleGc = ExpectedGcFraction(p.Sequence);
                if (Math.Abs(p.GcContent - oracleGc) > 1e-9)
                    return false.Label($"GcContent {p.GcContent} != oracle {oracleGc} for '{p.Sequence}'");

                // Relaxed hard band: [MinGc − 0.1, MaxGc + 0.1]. A small epsilon absorbs fp jitter.
                if (p.GcContent < param.MinGc - 0.1 - 1e-9 || p.GcContent > param.MaxGc + 0.1 + 1e-9)
                    return false.Label(
                        $"GcContent {p.GcContent} outside relaxed band " +
                        $"[{param.MinGc - 0.1},{param.MaxGc + 0.1}]");
            }

            return true.ToProperty();
        });
    }

    /// <summary>
    /// INV-P-GC-warning (equivalence): a returned probe carries a GC warning IFF its GcContent is
    /// strictly outside the configured [MinGc, MaxGc] band. The GC penalty branch in
    /// <c>EvaluateProbeWithGc</c> adds exactly the "GC content … outside range" warning precisely
    /// when <c>gc &lt; MinGc || gc &gt; MaxGc</c>. We recompute the GC fraction independently and
    /// verify the warning⇔outside-strict-range biconditional. This is the honest replacement for a
    /// (WRONG) strict-band assertion. Source: ProbeDesigner.cs §EvaluateProbeWithGc.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property DesignProbes_GcWarning_IffOutsideStrictBand()
    {
        return Prop.ForAll(ProbeScenarioArbitrary(), s =>
        {
            var (target, param) = s;
            var probes = ThermodynamicProbeDesigner.DesignProbes(target, param, maxProbes: 100000).ToList();

            foreach (var p in probes)
            {
                double gc = ExpectedGcFraction(p.Sequence);
                bool outsideStrict = gc < param.MinGc || gc > param.MaxGc;
                bool hasGcWarning = p.Warnings.Any(w => w.StartsWith("GC content"));

                if (outsideStrict != hasGcWarning)
                    return false.Label(
                        $"GC warning {hasGcWarning} != outsideStrict {outsideStrict} " +
                        $"(gc {gc}, band [{param.MinGc},{param.MaxGc}]) for '{p.Sequence}'");
            }

            return true.ToProperty();
        });
    }



    /// <summary>
    /// INV-M-Tm (subset monotonicity): Tm is only a −0.3 soft penalty (never a hard filter);
    /// probes are dropped solely when total Score ≤ 0, then the top <c>maxProbes</c> by score are
    /// returned. So for the SAME target and parameters identical except the Tm band, if window
    /// B ⊆ window A (B narrower) then count(B) ≤ count(A): each probe's score under B is ≤ its
    /// score under A (the only difference is possibly an extra −0.3 Tm penalty), so the
    /// positive-score set under B is a subset of that under A. A LARGE maxProbes (100000) is used
    /// so neither result saturates at the cap — otherwise both saturate and the relation is vacuous.
    /// This is the HONEST monotonic claim, NOT the naive "narrower window ⇒ fewer windows scanned".
    /// Source: ProbeDesigner.cs §EvaluateProbeWithGc Tm penalty + §DesignProbesOptimized score gate.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property DesignProbes_NarrowerTmWindow_NeverMoreProbes()
    {
        return Prop.ForAll(NestedTmScenarioArbitrary(), s =>
        {
            var (target, wide, narrow) = s;
            const int cap = 100000;
            int countWide = ThermodynamicProbeDesigner.DesignProbes(target, wide, maxProbes: cap).Count();
            int countNarrow = ThermodynamicProbeDesigner.DesignProbes(target, narrow, maxProbes: cap).Count();

            return (countNarrow <= countWide)
                .Label($"narrower Tm window produced MORE probes: " +
                       $"narrow [{narrow.MinTm},{narrow.MaxTm}]={countNarrow} > " +
                       $"wide [{wide.MinTm},{wide.MaxTm}]={countWide}");
        });
    }



    /// <summary>
    /// INV-D (Determinism): identical (target, params, maxProbes) ⇒ identical probe list —
    /// same sequences, starts, Tm, GcContent, Score, AND order. Two runs are compared element-wise.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property DesignProbes_IsDeterministic()
    {
        return Prop.ForAll(ProbeScenarioArbitrary(), s =>
        {
            var (target, param) = s;
            var a = ThermodynamicProbeDesigner.DesignProbes(target, param, maxProbes: 50).ToList();
            var b = ThermodynamicProbeDesigner.DesignProbes(target, param, maxProbes: 50).ToList();

            if (a.Count != b.Count)
                return false.Label($"count differs: {a.Count} vs {b.Count}");

            for (int i = 0; i < a.Count; i++)
            {
                if (a[i].Sequence != b[i].Sequence || a[i].Start != b[i].Start ||
                    !a[i].Tm.Equals(b[i].Tm) || !a[i].GcContent.Equals(b[i].GcContent) ||
                    !a[i].Score.Equals(b[i].Score))
                    return false.Label($"probe #{i} differs between runs");
            }

            return true.ToProperty();
        });
    }

    /// <summary>
    /// INV-Edge: null and empty targets yield no probes; a target shorter than MinLength yields
    /// no probes (the explicit early return in <see cref="ThermodynamicProbeDesigner.DesignProbes"/>).
    /// Source: ProbeDesigner.cs §DesignProbes guard; doc §6.1 edge cases.
    /// </summary>
    [Test]
    [Category("Property")]
    public void DesignProbes_NullEmptyOrTooShort_ReturnsNoProbes()
    {
        var param = ThermodynamicProbeDesigner.Defaults.qPCR; // MinLength 20
        Assert.Multiple(() =>
        {
            Assert.That(ThermodynamicProbeDesigner.DesignProbes(null!, param).ToList(), Is.Empty, "null target");
            Assert.That(ThermodynamicProbeDesigner.DesignProbes("", param).ToList(), Is.Empty, "empty target");
            Assert.That(ThermodynamicProbeDesigner.DesignProbes(new string('A', 19), param).ToList(), Is.Empty,
                "target shorter than MinLength (20)");
        });
    }



    /// <summary>
    /// Engineered length + GC anchor with hand-computed expectations. A 24-mer with exactly 50%
    /// GC, in a window forcing MinLength == MaxLength == 24 and a GC band that contains 0.50, must
    /// return a single probe whose Sequence is the whole target, whose length is 24, whose
    /// coordinates are 0..23, and whose GcContent equals the independently computed 12/24 = 0.5.
    /// Because GC is inside the strict band, NO GC warning is present.
    /// </summary>
    [Test]
    [Category("Property")]
    public void DesignProbes_EngineeredFiftyPercentGc_FixedLengthWindow()
    {
        // 24 bases: 12 G/C (ACGT repeated 6× → 6 G + 6 C = 12) ⇒ GC fraction = 12/24 = 0.5.
        const string target = "ACGTACGTACGTACGTACGTACGT";
        var param = new ThermodynamicProbeDesigner.ProbeParameters(
            MinLength: 24, MaxLength: 24,
            MinTm: 0, MaxTm: 200,        // wide Tm so no Tm penalty drops the probe
            MinGc: 0.40, MaxGc: 0.60,    // 0.50 is strictly inside
            MaxHomopolymer: 6, AvoidSecondaryStructure: false, MaxSelfComplementarity: 1.0);

        var probes = ThermodynamicProbeDesigner.DesignProbes(target, param, maxProbes: 10).ToList();

        Assert.That(probes, Has.Count.EqualTo(1), "exactly one length-24 window exists");
        var p = probes[0];
        Assert.Multiple(() =>
        {
            Assert.That(p.Sequence, Is.EqualTo(target));
            Assert.That(p.Sequence.Length, Is.EqualTo(24));
            Assert.That(p.Start, Is.EqualTo(0));
            Assert.That(p.End, Is.EqualTo(23));
            Assert.That(p.GcContent, Is.EqualTo(0.5).Within(1e-9), "12 G/C of 24 bases = 0.5");
            Assert.That(p.Warnings.Any(w => w.StartsWith("GC content")), Is.False,
                "GC 0.5 is inside [0.40,0.60] ⇒ no GC warning");
        });
    }

    /// <summary>
    /// Engineered warning anchor: a low-GC target (20% GC) whose fraction (0.20) is still inside
    /// the RELAXED band [MinGc − 0.1, MaxGc + 0.1] = [0.20, 0.70] but strictly BELOW the configured
    /// MinGc = 0.30 ⇒ the returned probe survives the early-rejection filter yet MUST carry a GC
    /// warning. Demonstrates the relaxed-band + warning-equivalence contract on a concrete case.
    /// </summary>
    [Test]
    [Category("Property")]
    public void DesignProbes_GcInsideRelaxedBandButBelowStrict_CarriesWarning()
    {
        // 20 bases, 4 of them G/C (G,C,G,C) ⇒ GC fraction = 4/20 = 0.20.
        const string target20 = "AAGCAAGCAAAAAAAAAAAA";
        var param = new ThermodynamicProbeDesigner.ProbeParameters(
            MinLength: 20, MaxLength: 20,
            MinTm: 0, MaxTm: 200,
            MinGc: 0.30, MaxGc: 0.60,   // strict band; relaxed band = [0.20, 0.70]
            MaxHomopolymer: 20, AvoidSecondaryStructure: false, MaxSelfComplementarity: 1.0);

        var probes = ThermodynamicProbeDesigner.DesignProbes(target20, param, maxProbes: 10).ToList();

        Assert.That(probes, Has.Count.EqualTo(1), "0.20 is the lower edge of the relaxed band ⇒ retained");
        var p = probes[0];
        Assert.Multiple(() =>
        {
            Assert.That(p.GcContent, Is.EqualTo(0.20).Within(1e-9), "4 G/C of 20 bases = 0.20");
            Assert.That(p.GcContent, Is.LessThan(param.MinGc), "0.20 < strict MinGc 0.30");
            Assert.That(p.Warnings.Any(w => w.StartsWith("GC content")), Is.True,
                "GC below strict MinGc ⇒ GC warning present");
        });
    }



    /// <summary>
    /// INV-01 (R): <c>SpecificityScore ∈ [0.0, 1.0]</c> for EVERY input. The score is mapped from
    /// a non-negative hit count to 0, 1, or 1/hits (all ≤ 1), so it can never leave the unit
    /// interval. A small epsilon absorbs fp jitter. Source: doc §2.4 INV-01; ProbeDesigner.cs
    /// §ValidateProbe specificity mapping.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property ValidateProbe_Specificity_InUnitInterval()
    {
        return Prop.ForAll(ProbeValidationScenarioArbitrary(), s =>
        {
            var (probe, refs, maxMismatches) = s;
            var v = ThermodynamicProbeDesigner.ValidateProbe(probe, refs, maxMismatches);
            return (v.SpecificityScore >= -1e-12 && v.SpecificityScore <= 1.0 + 1e-12)
                .Label($"SpecificityScore {v.SpecificityScore} outside [0,1] " +
                       $"(hits={v.OffTargetHits}, probe='{probe}')");
        });
    }

    /// <summary>
    /// INV-02 (R): <c>SelfComplementarity ∈ [0.0, 1.0]</c> AND equals the independent oracle
    /// (fraction of positions where probe[i] == revComp(probe)[i]) within 1e-9. The oracle uses its
    /// OWN complement table, so a wrong production complement map or counting bug is caught.
    /// Source: doc §2.4 INV-02; ProbeDesigner.cs §CalculateSelfComplementarity.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property ValidateProbe_SelfComplementarity_MatchesOracle_AndInRange()
    {
        return Prop.ForAll(ProbeValidationScenarioArbitrary(), s =>
        {
            var (probe, refs, maxMismatches) = s;
            var v = ThermodynamicProbeDesigner.ValidateProbe(probe, refs, maxMismatches);
            double oracle = OracleSelfComplementarity(probe);

            if (v.SelfComplementarity < -1e-12 || v.SelfComplementarity > 1.0 + 1e-12)
                return false.Label($"SelfComplementarity {v.SelfComplementarity} outside [0,1]");
            if (Math.Abs(v.SelfComplementarity - oracle) > 1e-9)
                return false.Label($"SelfComplementarity {v.SelfComplementarity} != oracle {oracle} " +
                                   $"for probe '{probe}'");
            return true.ToProperty();
        });
    }

    /// <summary>
    /// INV-03 (R, core business logic): <c>OffTargetHits ≥ 0</c> AND equals the independent
    /// substitution-only fixed-length sliding-window count summed across all references (the
    /// <see cref="OracleOffTargetHits"/> re-implementation from the doc, NOT production internals).
    /// This is the heart of the contract — the production count must match the oracle EXACTLY.
    /// Source: doc §2.4 INV-03, §4.1; ProbeDesigner.cs §FindApproximateMatches.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property ValidateProbe_OffTargetHits_MatchOracle_AndNonNegative()
    {
        return Prop.ForAll(ProbeValidationScenarioArbitrary(), s =>
        {
            var (probe, refs, maxMismatches) = s;
            var v = ThermodynamicProbeDesigner.ValidateProbe(probe, refs, maxMismatches);
            int oracle = OracleOffTargetHits(probe, refs, maxMismatches);

            if (v.OffTargetHits < 0)
                return false.Label($"OffTargetHits {v.OffTargetHits} is negative");
            return (v.OffTargetHits == oracle)
                .Label($"OffTargetHits {v.OffTargetHits} != oracle {oracle} " +
                       $"(probe='{probe}', maxMm={maxMismatches}, refs=[{string.Join(",", refs)}])");
        });
    }



    /// <summary>
    /// Specificity piecewise map (P): re-derived from the INDEPENDENT hit count, production's
    /// <c>SpecificityScore</c> must equal <c>0 hits ⇒ 0.0</c>, <c>1 hit ⇒ 1.0</c> (INV-04),
    /// <c>N&gt;1 hits ⇒ 1.0/N</c> (within 1e-9). This couples the score to the oracle hit count,
    /// not to production's own hits, so a mapping OR a counting error is caught. Source: doc §2.2
    /// piecewise definition, §2.4 INV-04.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property ValidateProbe_Specificity_FollowsPiecewiseMap()
    {
        return Prop.ForAll(ProbeValidationScenarioArbitrary(), s =>
        {
            var (probe, refs, maxMismatches) = s;
            var v = ThermodynamicProbeDesigner.ValidateProbe(probe, refs, maxMismatches);
            int oracleHits = OracleOffTargetHits(probe, refs, maxMismatches);
            double expected = OracleSpecificity(oracleHits);

            return (Math.Abs(v.SpecificityScore - expected) > 1e-9
                    ? false.Label($"Specificity {v.SpecificityScore} != map({oracleHits})={expected}")
                    : true.ToProperty());
        });
    }



    /// <summary>
    /// IsValid rule (R, pass/fail): <c>IsValid ⟺ Issues.Count == 0</c> (Probe_Validation.md §3.2 / §6:
    /// "IsValid is true exactly when no issue was recorded"; the former lenient rule "≤ 1 hit and
    /// fold-back fraction ≤ 0.4 ⇒ valid despite issues" had no source and was removed in B07 F22).
    /// Issues are recorded for &gt;1 off-target hits, a self-structure Tm above MaxStructureTm (or the
    /// fallback alignment limits) and secondary structure. We feed production's own Issues, so this
    /// isolates the BOOLEAN decision rule itself.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property ValidateProbe_IsValid_FollowsDecisionRule()
    {
        return Prop.ForAll(ProbeValidationScenarioArbitrary(), s =>
        {
            var (probe, refs, maxMismatches) = s;
            var v = ThermodynamicProbeDesigner.ValidateProbe(probe, refs, maxMismatches);
            bool expected = v.Issues.Count == 0;
            return (v.IsValid == expected)
                .Label($"IsValid {v.IsValid} != rule {expected} " +
                       $"(issues={v.Issues.Count}, hits={v.OffTargetHits}, selfComp={v.SelfComplementarity})");
        });
    }

    /// <summary>
    /// INV-D (Determinism): identical inputs ⇒ identical <see cref="ProbeDesigner.ProbeValidation"/> —
    /// every field (IsValid, SpecificityScore, OffTargetHits, SelfComplementarity,
    /// HasSecondaryStructure, and the Issues list element-wise) compared across two runs.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property ValidateProbe_IsDeterministic()
    {
        return Prop.ForAll(ProbeValidationScenarioArbitrary(), s =>
        {
            var (probe, refs, maxMismatches) = s;
            var a = ThermodynamicProbeDesigner.ValidateProbe(probe, refs, maxMismatches);
            var b = ThermodynamicProbeDesigner.ValidateProbe(probe, refs, maxMismatches);

            bool same = a.IsValid == b.IsValid
                        && a.SpecificityScore.Equals(b.SpecificityScore)
                        && a.OffTargetHits == b.OffTargetHits
                        && a.SelfComplementarity.Equals(b.SelfComplementarity)
                        && a.HasSecondaryStructure == b.HasSecondaryStructure
                        && a.Issues.SequenceEqual(b.Issues);
            return same.Label($"non-deterministic ProbeValidation for probe '{probe}'");
        });
    }



    /// <summary>
    /// Contract §3: a null probe sequence throws <see cref="ArgumentNullException"/>.
    /// </summary>
    [Test]
    [Category("Property")]
    public void ValidateProbe_NullProbe_Throws()
    {
        Assert.Throws<ArgumentNullException>(
            () => ThermodynamicProbeDesigner.ValidateProbe(null!, new[] { "ACGT" }));
    }

    /// <summary>
    /// Contract §3: a null reference collection throws <see cref="ArgumentNullException"/>.
    /// </summary>
    [Test]
    [Category("Property")]
    public void ValidateProbe_NullReferences_Throws()
    {
        Assert.Throws<ArgumentNullException>(
            () => ThermodynamicProbeDesigner.ValidateProbe("ACGT", null!));
    }

    /// <summary>
    /// Contract §3.3 / §6.1: an empty probe yields a structured invalid result —
    /// IsValid == false, SpecificityScore == 0.0, OffTargetHits == 0, and Issues contains the exact
    /// "Empty probe sequence" entry. This is a defined special case, NOT an exception.
    /// </summary>
    [Test]
    [Category("Property")]
    public void ValidateProbe_EmptyProbe_StructuredInvalidResult()
    {
        var v = ThermodynamicProbeDesigner.ValidateProbe("", new[] { "ACGTACGT", "TTTTTTTT" });
        Assert.Multiple(() =>
        {
            Assert.That(v.IsValid, Is.False, "empty probe must be invalid");
            Assert.That(v.SpecificityScore, Is.EqualTo(0.0).Within(1e-12));
            Assert.That(v.OffTargetHits, Is.EqualTo(0));
            Assert.That(v.SelfComplementarity, Is.EqualTo(0.0).Within(1e-12));
            Assert.That(v.HasSecondaryStructure, Is.False);
            Assert.That(v.Issues, Does.Contain("Empty probe sequence"));
        });
    }



    /// <summary>
    /// Engineered exact-occurrence anchor: an 8-mer probe present EXACTLY TWICE across the references
    /// (once in each reference, no near-matches in tolerance 0) ⇒ OffTargetHits hand-counted as 2 and
    /// Specificity = 1/2 = 0.5. The flanks are chosen so no additional ≤0-mismatch window exists.
    /// </summary>
    [Test]
    [Category("Property")]
    public void ValidateProbe_ProbeTwiceAcrossRefs_HitsTwo_SpecificityHalf()
    {
        const string probe = "ACGTACGT"; // 8-mer
        // Two references, each containing the probe exactly once; flanks avoid extra exact windows.
        var refs = new[]
        {
            "TTTACGTACGTTTT", // probe at index 3, once
            "GGGACGTACGTGGG"  // probe at index 3, once
        };

        var v = ThermodynamicProbeDesigner.ValidateProbe(probe, refs, maxMismatches: 0);

        // Independent hand count: exact (0-mismatch) windows of "ACGTACGT".
        int oracleHits = OracleOffTargetHits(probe, refs, 0);
        Assert.Multiple(() =>
        {
            Assert.That(oracleHits, Is.EqualTo(2), "oracle: probe occurs exactly twice");
            Assert.That(v.OffTargetHits, Is.EqualTo(2), "production matches hand count");
            Assert.That(v.SpecificityScore, Is.EqualTo(0.5).Within(1e-9), "1/2 = 0.5");
        });
    }

    /// <summary>
    /// Engineered single-mismatch anchor: a unique reference window that differs from the probe by
    /// exactly ONE substitution is counted as a hit when maxMismatches ≥ 1 (⇒ 1 hit, Specificity
    /// 1.0) but NOT counted when maxMismatches == 0 (⇒ 0 hits, Specificity 0.0). Demonstrates the
    /// substitution-tolerant window contract on a concrete case.
    /// </summary>
    [Test]
    [Category("Property")]
    public void ValidateProbe_SingleMismatchWithinTolerance_CountedAsHit()
    {
        const string probe = "ACGTACGT";
        // Reference contains the probe with one substitution at position 0 (A→T), unique window.
        var refs = new[] { "GG" + "TCGTACGT" + "GG" };

        var withTolerance = ThermodynamicProbeDesigner.ValidateProbe(probe, refs, maxMismatches: 1);
        var noTolerance = ThermodynamicProbeDesigner.ValidateProbe(probe, refs, maxMismatches: 0);

        Assert.Multiple(() =>
        {
            Assert.That(OracleOffTargetHits(probe, refs, 1), Is.EqualTo(1), "oracle: 1 hit at mm=1");
            Assert.That(withTolerance.OffTargetHits, Is.EqualTo(1), "1-mismatch window counted at mm=1");
            Assert.That(withTolerance.SpecificityScore, Is.EqualTo(1.0).Within(1e-9), "unique ⇒ 1.0");
            Assert.That(OracleOffTargetHits(probe, refs, 0), Is.EqualTo(0), "oracle: 0 hits at mm=0");
            Assert.That(noTolerance.OffTargetHits, Is.EqualTo(0), "1-mismatch window NOT counted at mm=0");
            Assert.That(noTolerance.SpecificityScore, Is.EqualTo(0.0).Within(1e-9), "no hit ⇒ 0.0");
        });
    }




    /// <summary>
    /// Independent GC fraction oracle: (#G + #C) / length on the uppercased sequence, computed
    /// from raw base counts — NOT routed through any production helper. Returns 0 for empty input.
    /// (The production code stores GcContent as a fraction in [0,1], cf. INV-02.)
    /// </summary>
    private static double ExpectedGcFraction(string seq)
    {
        if (seq.Length == 0) return 0;
        var (_, gc) = CountAtGc(seq);
        return (double)gc / seq.Length;
    }

    /// <summary>
    /// Independent self-complementarity oracle (INV-02): fraction of positions i where
    /// <c>seq[i] == reverseComplement(seq)[i]</c>, over the uppercased probe. Mirrors the
    /// definition in ProbeDesigner §CalculateSelfComplementarity but with its OWN revcomp.
    /// </summary>
    private static double OracleSelfComplementarity(string seq)
    {
        string up = seq.ToUpperInvariant();
        string rc = OracleReverseComplement(up);
        int matches = 0;
        for (int i = 0; i < up.Length; i++)
            if (up[i] == rc[i]) matches++;
        return up.Length == 0 ? 0.0 : matches / (double)up.Length;
    }

    /// <summary>
    /// Independent off-target counter (INV-03): substitution-only, fixed-length, ungapped sliding
    /// window over each reference (both probe and reference uppercased), counting every start
    /// position whose Hamming distance to the probe over the window is ≤ maxMismatches, summed
    /// across all references. Re-implemented from the doc §4.1 / §5.2 — NOT from production
    /// internals — so it is a true oracle for <c>FindApproximateMatches</c>.
    /// </summary>
    private static int OracleOffTargetHits(string probe, IEnumerable<string> references, int maxMismatches)
    {
        string p = probe.ToUpperInvariant();
        int total = 0;
        foreach (var reference in references)
        {
            string text = reference.ToUpperInvariant();
            for (int i = 0; i + p.Length <= text.Length; i++)
            {
                int mismatches = 0;
                for (int j = 0; j < p.Length; j++)
                {
                    if (text[i + j] != p[j]) mismatches++;
                    if (mismatches > maxMismatches) break;
                }
                if (mismatches <= maxMismatches) total++;
            }
        }
        return total;
    }

    /// <summary>
    /// Independent specificity map (INV-01 / INV-04 / §2.2): 0 hits ⇒ 0.0, 1 hit ⇒ 1.0,
    /// N&gt;1 hits ⇒ 1.0/N. Re-derived purely from a hit count, not from production.
    /// </summary>
    private static double OracleSpecificity(int offTargetHits) =>
        offTargetHits == 0 ? 0.0 : offTargetHits == 1 ? 1.0 : 1.0 / offTargetHits;


    /// <summary>
    /// Counts the bases the basic Tm uses (case-insensitive), partitioned into A/T and G/C; U counts
    /// as T (Biopython <c>MeltingTemp._check</c> back-transcription, ThermoConstants.CountBasicTmBases).
    /// </summary>
    private static (int at, int gc) CountAtGc(string s)
    {
        int at = 0, gc = 0;
        foreach (char ch in s.ToUpperInvariant())
        {
            if (ch is 'A' or 'T' or 'U') at++;
            else if (ch is 'G' or 'C') gc++;
        }
        return (at, gc);
    }



    /// <summary>
    /// Independent reverse-complement oracle over an uppercased ACGT string, built char-by-char
    /// from the complement map A↔T, C↔G and then reversed — deliberately NOT routed through any
    /// production helper, so a wrong complement table is caught. Non-ACGT chars are mapped to 'N'.
    /// </summary>
    private static string OracleReverseComplement(string seq)
    {
        var chars = new char[seq.Length];
        for (int i = 0; i < seq.Length; i++)
        {
            char c = char.ToUpperInvariant(seq[i]);
            chars[seq.Length - 1 - i] = c switch
            {
                'A' => 'T',
                'T' => 'A',
                'C' => 'G',
                'G' => 'C',
                _ => 'N'
            };
        }
        return new string(chars);
    }




    /// <summary>
    /// INV-R-len: Every VALID candidate returned by <see cref="ThermodynamicPrimerDesigner.GeneratePrimerCandidates"/>
    /// has Length within the configured [MinLength, MaxLength] window, and an out-of-window length is
    /// never accepted. Exercised across randomized templates AND randomized length windows so the
    /// bound is a true config-driven invariant, not a single-template coincidence.
    /// Source: Primer_Design.md §3 contract (length 18-25 default; window from parameters).
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property GeneratePrimerCandidates_ValidLength_InConfiguredWindow()
    {
        return Prop.ForAll(CandidateScenarioArbitrary(), s =>
        {
            var dna = new DnaSequence(s.template);
            var candidates = ThermodynamicPrimerDesigner
                .GeneratePrimerCandidates(dna, s.regionStart, s.regionEnd, s.forward, s.param)
                .ToList();

            bool allValidInWindow = candidates
                .Where(c => c.IsValid)
                .All(c => c.Length >= s.param.MinLength && c.Length <= s.param.MaxLength);

            return allValidInWindow.Label(
                $"valid candidate length outside [{s.param.MinLength},{s.param.MaxLength}] " +
                $"in region len {s.template.Length}");
        });
    }

    /// <summary>
    /// INV-P-GC: Every VALID candidate has GcContent (a percentage 0-100 in source) within the
    /// configured [MinGcContent, MaxGcContent] window; a candidate whose GC% falls outside the window
    /// is never marked valid. GcContent is recomputed independently from base counts (ExpectedGcPercent)
    /// to confirm the production value used for filtering is itself correct.
    /// Source: Primer_Design.md §3 / §4.2 (GC 40-60 default; window from parameters).
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property GeneratePrimerCandidates_ValidGcContent_InConfiguredWindow()
    {
        return Prop.ForAll(CandidateScenarioArbitrary(), s =>
        {
            var dna = new DnaSequence(s.template);
            var candidates = ThermodynamicPrimerDesigner
                .GeneratePrimerCandidates(dna, s.regionStart, s.regionEnd, s.forward, s.param)
                .ToList();

            foreach (var c in candidates)
            {
                // GcContent is rounded to 1 decimal in source; oracle within that tolerance.
                if (Math.Abs(c.GcContent - ExpectedGcPercent(c.Sequence)) > 0.06)
                    return false.Label($"GcContent {c.GcContent} != oracle {ExpectedGcPercent(c.Sequence)} for {c.Sequence}");

                if (c.IsValid &&
                    (c.GcContent < s.param.MinGcContent || c.GcContent > s.param.MaxGcContent))
                    return false.Label(
                        $"valid candidate GC {c.GcContent}% outside [{s.param.MinGcContent},{s.param.MaxGcContent}]");
            }

            return true.ToProperty();
        });
    }

    /// <summary>
    /// INV-R-Tm: Every VALID candidate has MeltingTemperature within the configured [MinTm, MaxTm]
    /// window; a candidate whose Tm is outside the window is never marked valid.
    /// Source: Primer_Design.md §3 (Tm 57-63 default; window from parameters).
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property GeneratePrimerCandidates_ValidTm_InConfiguredWindow()
    {
        return Prop.ForAll(CandidateScenarioArbitrary(), s =>
        {
            var dna = new DnaSequence(s.template);
            var candidates = ThermodynamicPrimerDesigner
                .GeneratePrimerCandidates(dna, s.regionStart, s.regionEnd, s.forward, s.param)
                .ToList();

            bool allValidInWindow = candidates
                .Where(c => c.IsValid)
                .All(c => c.MeltingTemperature >= s.param.MinTm && c.MeltingTemperature <= s.param.MaxTm);

            return allValidInWindow.Label(
                $"valid candidate Tm outside [{s.param.MinTm},{s.param.MaxTm}]");
        });
    }


    /// <summary>Independent GC% oracle: 100 × (#G + #C) / length, computed from raw base counts.</summary>
    private static double ExpectedGcPercent(string seq)
    {
        if (seq.Length == 0) return 0;
        var (_, gc) = CountAtGc(seq);
        return 100.0 * gc / seq.Length;
    }

    /// <summary>
    /// Generates a probe-validation scenario: a random ACGT probe (8-20 bp), a set of 0-4 reference
    /// sequences, and a mismatch tolerance. References are deliberately MIXED in kind so the
    /// off-target counter is exercised across the interesting regimes:
    /// (a) references that EMBED the probe verbatim 0/1/several times (exact hits),
    /// (b) references that embed a NEAR-match (the probe with 1-3 substitutions — within tolerance),
    /// (c) unrelated filler (no hits). Lengths are kept short so hundreds of cases run quickly.
    /// </summary>
    private static Arbitrary<(string probe, string[] refs, int maxMismatches)>
        ProbeValidationScenarioArbitrary() =>
        (from plen in Gen.Choose(8, 20)
         from pchars in Gen.Elements('A', 'C', 'G', 'T').ArrayOf(plen)
         let probe = new string(pchars)
         from maxMismatches in Gen.Choose(0, 4)
         from refCount in Gen.Choose(0, 4)
         from refs in BuildReferenceGen(probe, maxMismatches).ArrayOf(refCount)
         select (probe, refs, maxMismatches)).ToArbitrary();


    /// <summary>
    /// Builds a single reference sequence around the probe: random ACGT flanks plus, with varying
    /// probability, an embedded exact copy of the probe or a perturbed copy carrying a controlled
    /// number of substitutions (0..probeLen). The perturbed copy may or may not be within tolerance,
    /// so both "counted" and "not counted" near-matches are produced.
    /// </summary>
    private static Gen<string> BuildReferenceGen(string probe, int maxMismatches)
    {
        return
            from leftLen in Gen.Choose(0, 8)
            from rightLen in Gen.Choose(0, 8)
            from left in Gen.Elements('A', 'C', 'G', 'T').ArrayOf(leftLen)
            from right in Gen.Elements('A', 'C', 'G', 'T').ArrayOf(rightLen)
            from kind in Gen.Choose(0, 3) // 0=no insert, 1=exact, 2=near-match, 3=double exact
            from subs in Gen.Choose(0, probe.Length)
            let middle = kind switch
            {
                1 => probe,
                2 => PerturbProbe(probe, subs),
                3 => probe + probe,
                _ => string.Empty
            }
            select new string(left) + middle + new string(right);
    }


    /// <summary>
    /// Returns a copy of the probe with the first <paramref name="count"/> positions substituted to
    /// a DIFFERENT base (deterministic 'A'→'C', else→'A'), producing exactly <paramref name="count"/>
    /// guaranteed mismatches at known positions.
    /// </summary>
    private static string PerturbProbe(string probe, int count)
    {
        var chars = probe.ToCharArray();
        for (int i = 0; i < count && i < chars.Length; i++)
            chars[i] = chars[i] == 'A' ? 'C' : 'A';
        return new string(chars);
    }
}
