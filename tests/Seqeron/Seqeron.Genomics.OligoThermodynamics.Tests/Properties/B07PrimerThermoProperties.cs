// SPDX-License-Identifier: GPL-2.0-or-later
// Tests of the port of Primer3 ntthal (thal.c) and of the Primer3-style primer / probe design; part of Seqeron.Genomics.OligoThermodynamics.Tests, licensed under the GNU General Public License, version 2 or (at your option) any later version (see LICENSE in the production project).
using FsCheck;
using FsCheck.Fluent;
using static Seqeron.Genomics.MolTools.PrimerDesigner;

namespace Seqeron.Genomics.Tests.Properties;

/// <summary>
/// Heavy-tier property tests for the PRIMER / ntthal / ThermoConstants behaviour introduced by review batch B07
/// (docs/Validation/review-2026-09/B07.md): F55 (THAL_MAX_ALIGN opt-in), F57/F58 (StructureModel.Ntthal overloads and
/// the single-stem core's stem enumeration), F61/F64 (HasHairpinPotential on both sides of the 100-nt suffix-tree
/// switch), F62 (Primer3 Tm max_nn_length), F68a (NN shift beyond the overlap) and F72 (self-complementarity).
/// Oracles are independent re-computations (reverse complement, stem existence, substring restriction) or the
/// documented equalities between public overloads — never copies of the implementation.
///
/// Test Units: PRIMER-HAIRPIN-001, PRIMER-DIMER-001, PRIMER-STRUCT-001, PRIMER-DESIGN-001, PRIMER-NNTM-001.
/// </summary>
[TestFixture]
[Category("Property")]
[Category("MolTools")]
public class B07PrimerThermoProperties
{



    /// <summary>
    /// F57: <c>FindMostStableHairpin(seq, Ntthal, …)</c> ΔH/ΔS/ΔG and <c>CalculateHairpinMeltingTemperature(seq, Ntthal, …)</c>
    /// equal <c>CalculateHairpinThermodynamicsNtthal(seq, mv ?? 0.05, dv ?? 0.0015, dntp ?? 0.0006, temp ?? 37, maxLoop ?? 30)</c>
    /// (documented delegation; null ⇔ no structure ⇔ NaN Tm), and the mapped hairpin is a well-formed structure:
    /// 0 ≤ StemStart &lt; StemEnd &lt; n, StemLength ≥ 1, LoopSize ≥ 3 (thal.c <c>min_hrpn_loop</c>) and the span holds
    /// the stem's 2·StemLength paired bases plus the loop.
    /// </summary>
    [FsCheck.NUnit.Property(MaxTest = 150)]
    public Property FindMostStableHairpin_NtthalModel_EqualsHairpinThermodynamicsNtthal()
    {
        var gen = from s in HairpinOligoGen() from c in ConditionsGen() select (s, c);
        return Prop.ForAll(gen.ToArbitrary(), x =>
        {
            var (s, c) = x;
            var r = FindMostStableHairpin(s, StructureModel.Ntthal, c.Mv, c.Dv, c.Dntp, c.T, c.MaxLoop);
            double tm = CalculateHairpinMeltingTemperature(s, StructureModel.Ntthal, c.Mv, c.Dv, c.Dntp, c.T, c.MaxLoop);
            var t = CalculateHairpinThermodynamicsNtthal(s, c.Mv ?? 0.05, c.Dv ?? 0.0015, c.Dntp ?? 0.0006,
                c.T ?? NtthalDefaultTemperatureCelsius, c.MaxLoop ?? NtthalDefaultMaxLoop);
            if (t is null)
                return (r is null && double.IsNaN(tm)).Label($"no structure but r={r}, tm={tm} for {s}");
            if (r is null)
                return false.Label($"structure {t} but FindMostStableHairpin null for {s}");
            var h = r.Value;
            var v = t.Value;
            bool thermo = Close(h.DeltaH, v.DeltaH) && Close(h.DeltaS, v.DeltaS) && Close(h.DeltaG37, v.DeltaG37)
                          && Close(tm, v.TmCelsius);
            bool shape = h.StemStart >= 0 && h.StemStart < h.StemEnd && h.StemEnd < s.Length && h.StemLength >= 1
                         && h.LoopSize >= 3 && h.StemEnd - h.StemStart + 1 >= 2 * h.StemLength + h.LoopSize
                         && Wc(s[h.StemStart], s[h.StemEnd]);
            return (thermo && shape).Label($"{s} {c}: hairpin {h} tm {tm:R} vs {v}");
        });
    }

    /// <summary>
    /// F58: <c>FindMostStableDimer(a, b, Ntthal, …)</c> ΔH/ΔS/ΔG/BasePairs equal
    /// <c>CalculateDimerThermodynamicsNtthal(a, b, Any, mv ?? 0.05, dv ?? 0.0015, dntp ?? 0.0006, C ?? 50 nM, temp ?? 37,
    /// maxLoop ?? 30)</c>; the reported 5′-most paired bases lie inside their strands and pair Watson–Crick.
    /// </summary>
    [FsCheck.NUnit.Property(MaxTest = 150)]
    public Property FindMostStableDimer_NtthalModel_EqualsDimerThermodynamicsNtthal()
    {
        var mutatedRc =
            from a in AcgtGen(6, 40)
            from k in Gen.Choose(0, 3)
            from positions in Gen.Choose(0, 1000).ArrayOf(k)
            from bases in Gen.Elements('A', 'C', 'G', 'T').ArrayOf(k)
            select (a, Mutate(RevComp(a), positions, bases));
        var pairGen = Gen.Frequency(
            (1, from a in AcgtGen(5, 40) select (a, a)),
            (2, mutatedRc),
            (2, from a in AcgtGen(5, 40) from b in AcgtGen(5, 40) select (a, b)));
        var gen = from p in pairGen
                  from c in ConditionsGen()
                  from conc in Gen.Elements<double?>(null, 1e-9, 50e-9, 250e-9, 2e-6)
                  select (p.Item1, p.Item2, c, conc);
        return Prop.ForAll(gen.ToArbitrary(), x =>
        {
            var (a, b, c, conc) = x;
            var r = FindMostStableDimer(a, b, StructureModel.Ntthal, c.Mv, c.Dv, c.Dntp, conc, c.T, c.MaxLoop);
            var t = CalculateDimerThermodynamicsNtthal(a, b, NtthalAlignmentMode.Any, c.Mv ?? 0.05, c.Dv ?? 0.0015,
                c.Dntp ?? 0.0006, conc ?? 50e-9, c.T ?? NtthalDefaultTemperatureCelsius, c.MaxLoop ?? NtthalDefaultMaxLoop);
            if (t is null)
                return (r is null).Label($"no structure but r={r} for {a}/{b}");
            if (r is null)
                return false.Label($"structure {t} but FindMostStableDimer null for {a}/{b}");
            var d = r.Value;
            var v = t.Value;
            bool thermo = Close(d.DeltaH, v.DeltaH) && Close(d.DeltaS, v.DeltaS) && Close(d.DeltaG37, v.DeltaG37)
                          && d.BasePairs == v.BasePairs;
            bool shape = d.Strand1Start >= 0 && d.Strand1Start < a.Length && d.Strand2Start >= 0
                         && d.Strand2Start < b.Length && d.BasePairs >= 1;
            return (thermo && shape).Label($"{a}/{b} {c} C={conc}: {d} vs {v}");
        });
    }



    /// <summary>
    /// F55/F56: whenever at least one strand is ≤ 60 nt (the default THAL_MAX_ALIGN accepts it), every
    /// <c>maxAlignLength</c> overload (60–10000) of the ntthal hairpin / dimer thermodynamics and structure methods equals
    /// the default overload exactly (thermodynamics, no-structure and ASCII structure lines).
    /// </summary>
    [FsCheck.NUnit.Property(MaxTest = 100)]
    public Property NtthalMaxAlignOverloads_EqualDefault_ForStrandsUpTo60()
    {
        var gen = from s1 in Gen.Frequency((1, AcgtGen(5, 60)), (1, StemLoopGen()))
                  from s2 in Gen.Frequency((2, AcgtGen(5, 60)), (1, AcgtGen(61, 90)),
                      (2, from k in Gen.Choose(0, 2) select RevComp(s1).Substring(0, Math.Max(5, s1.Length - k))))
                  from mode in Gen.Elements(NtthalAlignmentMode.Any, NtthalAlignmentMode.End1, NtthalAlignmentMode.End2)
                  from maxAlign in Gen.Elements(60, 61, 75, 200, 1000, 10000)
                  from c in ConditionsGen()
                  select (s1, s2, mode, maxAlign, c);
        return Prop.ForAll(gen.ToArbitrary(), x =>
        {
            var (s1, s2, mode, maxAlign, c) = x;
            double mv = c.Mv ?? 0.05, dv = c.Dv ?? 0.0015, dntp = c.Dntp ?? 0.0006;
            double t = c.T ?? NtthalDefaultTemperatureCelsius;
            int ml = c.MaxLoop ?? NtthalDefaultMaxLoop;

            var hDefault = CalculateHairpinThermodynamicsNtthal(s1, mv, dv, dntp, t, ml);
            var hOpt = CalculateHairpinThermodynamicsNtthal(s1, mv, dv, dntp, t, ml, maxAlign);
            var hsDefault = CalculateHairpinStructureNtthal(s1, mv, dv, dntp, t, ml);
            var hsOpt = CalculateHairpinStructureNtthal(s1, mv, dv, dntp, t, ml, maxAlign);

            var dDefault = CalculateDimerThermodynamicsNtthal(s1, s2, mode, mv, dv, dntp, 50e-9, t, ml);
            var dOpt = CalculateDimerThermodynamicsNtthal(s1, s2, mode, mv, dv, dntp, 50e-9, t, ml, maxAlign);
            var dsDefault = CalculateDimerStructureNtthal(s1, s2, mode, mv, dv, dntp, 50e-9, t, ml);
            var dsOpt = CalculateDimerStructureNtthal(s1, s2, mode, mv, dv, dntp, 50e-9, t, ml, maxAlign);

            bool ok = Equals(hDefault, hOpt) && Equals(dDefault, dOpt)
                      && SameStructure(hsDefault?.Thermodynamics, hsDefault?.AsciiStructureLines,
                          hsOpt?.Thermodynamics, hsOpt?.AsciiStructureLines)
                      && SameStructure(dsDefault?.Thermodynamics, dsDefault?.AsciiStructureLines,
                          dsOpt?.Thermodynamics, dsOpt?.AsciiStructureLines);
            return ok.Label($"{s1}/{s2} {mode} maxAlign {maxAlign} {c}: hairpin {hDefault} vs {hOpt}; dimer {dDefault} vs {dOpt}");
        });
    }


    private static Gen<string> AcgtGen(int minLen, int maxLen) =>
        from n in Gen.Choose(minLen, maxLen)
        from chars in Gen.Elements('A', 'C', 'G', 'T').ArrayOf(n)
        select new string(chars);

    /// <summary>Independent Watson–Crick reverse complement (A↔T, C↔G; upper-case ACGT only).</summary>
    private static string RevComp(string s)
    {
        var r = new char[s.Length];
        for (int i = 0; i < s.Length; i++)
            r[s.Length - 1 - i] = s[i] switch { 'A' => 'T', 'T' => 'A', 'C' => 'G', 'G' => 'C', _ => 'N' };
        return new string(r);
    }

    private static bool Wc(char a, char b) =>
        (a, b) is ('A', 'T') or ('T', 'A') or ('C', 'G') or ('G', 'C');

    /// <summary>Designed stem-loop: flank + arm + loop + revcomp(arm) + flank, optionally one point mutation; ≤ 60 nt.</summary>
    private static Gen<string> StemLoopGen() =>
        from f5 in AcgtGen(0, 8)
        from arm in AcgtGen(3, 10)
        from loop in AcgtGen(3, 8)
        from f3 in AcgtGen(0, 8)
        from mutate in Gen.Elements(false, false, true)
        from pos in Gen.Choose(0, 1000)
        from b in Gen.Elements('A', 'C', 'G', 'T')
        let s = f5 + arm + loop + RevComp(arm) + f3
        select mutate ? s.Remove(pos % s.Length, 1).Insert(pos % s.Length, b.ToString()) : s;

    private static Gen<string> HairpinOligoGen() =>
        Gen.Frequency((2, AcgtGen(5, 60)), (3, StemLoopGen()));

    /// <summary>Ntthal conditions: all null (defaults) or random mv/dv/dNTP (mol/L), temp_c and max_loop.</summary>
    private static Gen<(double? Mv, double? Dv, double? Dntp, double? T, int? MaxLoop)> ConditionsGen() =>
        Gen.Frequency(
            (1, Gen.Constant<(double?, double?, double?, double?, int?)>((null, null, null, null, null))),
            (1, from mv in Gen.Choose(5, 500)
                from dv in Gen.Choose(0, 100)
                from dntp in Gen.Choose(0, 40)
                from t in Gen.Choose(20, 80)
                from ml in Gen.Choose(0, 30)
                select ((double?)(mv / 1000.0), (double?)(dv / 10000.0), (double?)(dntp / 10000.0), (double?)t, (int?)ml)));

    private static bool Close(double a, double b) =>
        a == b || Math.Abs(a - b) <= 1e-9 * Math.Max(1.0, Math.Max(Math.Abs(a), Math.Abs(b)));

    private static string Mutate(string s, int[] positions, char[] bases)
    {
        var a = s.ToCharArray();
        for (int i = 0; i < positions.Length; i++)
            a[positions[i] % a.Length] = bases[i];
        return new string(a);
    }

    private static bool SameStructure<T>(T? t1, IReadOnlyList<string>? l1, T? t2, IReadOnlyList<string>? l2)
        where T : struct =>
        Equals(t1, t2) && (l1 is null ? l2 is null : l2 is not null && l1.SequenceEqual(l2));




    /// <summary>
    /// F62: the five-argument <c>CalculateMeltingTemperaturePrimer3</c> is primer3-py <c>calc_tm(…, max_nn_length=36)</c>,
    /// so it equals the <c>maxNearestNeighborLength</c> overload with 36 for every oligo and condition; for oligos of
    /// ≤ 36 nt the NN path is taken under both limits, so the overload with calc_tm's default 60 gives the same value.
    /// </summary>
    [FsCheck.NUnit.Property(MaxTest = 200)]
    public Property MeltingTemperaturePrimer3_MaxNn36_EqualsFiveArgOverload()
    {
        var gen = from s in AcgtGen(2, 70)
                  from dna in Gen.Choose(1, 2000)
                  from mv in Gen.Choose(1, 500)
                  from dv in Gen.Choose(0, 100)
                  from dntp in Gen.Choose(0, 40)
                  select (s, (double)dna, (double)mv, dv / 10.0, dntp / 10.0);
        return Prop.ForAll(gen.ToArbitrary(), x =>
        {
            var (s, dna, mv, dv, dntp) = x;
            double five = CalculateMeltingTemperaturePrimer3(s, dna, mv, dv, dntp);
            double nn36 = CalculateMeltingTemperaturePrimer3(s, dna, mv, dv, dntp, 36);
            double nn60 = CalculateMeltingTemperaturePrimer3(s, dna, mv, dv, dntp, 60);
            bool ok = five.Equals(nn36) && (s.Length > 36 || five.Equals(nn60));
            return ok.Label($"{s} ({s.Length} nt): 5-arg {five:R}, nn36 {nn36:R}, nn60 {nn60:R}");
        });
    }
}
