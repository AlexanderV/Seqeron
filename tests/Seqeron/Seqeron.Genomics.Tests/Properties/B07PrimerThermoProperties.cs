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
    #region Generators and independent helpers

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

    #endregion

    #region F57 / F58 — StructureModel.Ntthal overloads = the ntthal thermodynamics methods

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

    private static string Mutate(string s, int[] positions, char[] bases)
    {
        var a = s.ToCharArray();
        for (int i = 0; i < positions.Length; i++)
            a[positions[i] % a.Length] = bases[i];
        return new string(a);
    }

    #endregion

    #region F57 — single-stem (SingleHelix) core enumerates every stem length

    /// <summary>Independent existence oracle: some (i, j) with (i, j) and (i+1, j−1) Watson–Crick and ≥ 3 bases between them.</summary>
    private static bool HasTwoBpStemClosingThreeNtLoop(string s)
    {
        for (int i = 0; i < s.Length; i++)
            for (int j = i + 6; j < s.Length; j++) // loop = (j−1) − (i+1) − 1 ≥ 3 ⇔ j ≥ i + 6
                if (Wc(s[i], s[j]) && Wc(s[i + 1], s[j - 1]))
                    return true;
        return false;
    }

    private static Gen<string> SingleHelixGen() =>
        Gen.Frequency(
            (3, AcgtGen(4, 30)),
            (2, StemLoopGen()),
            // Palindromic G^kC^k / (GC)^k: the maximal stem closes no loop, shorter stems do (F57 repro family).
            (1, from k in Gen.Choose(2, 8) select new string('G', k) + new string('C', k)),
            (1, from k in Gen.Choose(2, 8) select string.Concat(Enumerable.Repeat("GC", k))));

    /// <summary>
    /// F57 INV-02 both ways: the single-stem core returns a hairpin ⇔ some ≥ 2-bp Watson–Crick stem closes a ≥ 3-nt
    /// loop (before F57 GGGGCCCC / GCGCGCGCGC returned null), and the returned hairpin is geometrically exact:
    /// StemLength ≥ 2 consecutive WC pairs from (StemStart, StemEnd) inward closing LoopSize ≥ 3 unpaired bases.
    /// </summary>
    [FsCheck.NUnit.Property(MaxTest = 300)]
    public Property FindMostStableHairpin_SingleHelix_NonNullIffStemClosesLoop()
    {
        return Prop.ForAll(SingleHelixGen().ToArbitrary(), s =>
        {
            var r = FindMostStableHairpin(s);
            bool exists = HasTwoBpStemClosingThreeNtLoop(s);
            if (r is null)
                return (!exists).Label($"null but a 2-bp stem closes a ≥ 3-nt loop in {s}");
            var h = r.Value;
            bool geometry = h.StemLength >= 2 && h.LoopSize >= 3
                            && h.StemEnd - h.StemStart + 1 == 2 * h.StemLength + h.LoopSize
                            && Enumerable.Range(0, h.StemLength).All(t => Wc(s[h.StemStart + t], s[h.StemEnd - t]));
            return (exists && geometry).Label($"{s}: {h} (exists={exists})");
        });
    }

    /// <summary>
    /// F57 "never beaten by a shorter stem": the candidate set of <c>minStemLength = k</c> shrinks as k grows, so
    /// ΔG(k) = ΔG(2) for every k ≤ the optimum's StemLength (the optimum stays a candidate) and ΔG(k) ≥ ΔG(2) (or null)
    /// beyond it; restricting the sequence to the optimum's closing span keeps every shorter stem of that closing
    /// pair as a candidate and reproduces the optimum exactly, and no substring is ever more stable.
    /// </summary>
    [FsCheck.NUnit.Property(MaxTest = 200)]
    public Property FindMostStableHairpin_SingleHelix_ShorterStemsNeverBeatTheOptimum()
    {
        var gen = from s in Gen.Frequency((1, AcgtGen(8, 30)), (3, StemLoopGen()))
                  from a in Gen.Choose(0, 1000)
                  from b in Gen.Choose(0, 1000)
                  select (s, a, b);
        return Prop.ForAll(gen.ToArbitrary(), x =>
        {
            var (s, ra, rb) = x;
            var r = FindMostStableHairpin(s);
            if (r is null)
                return true.ToProperty();
            var h = r.Value;
            for (int k = 2; k <= h.StemLength; k++)
            {
                var rk = FindMostStableHairpin(s, k);
                if (rk is null || rk.Value.DeltaG37 != h.DeltaG37)
                    return false.Label($"minStem {k} ≤ {h.StemLength}: {rk} vs {h} for {s}");
            }
            var beyond = FindMostStableHairpin(s, h.StemLength + 1);
            if (beyond is { } bh && bh.DeltaG37 < h.DeltaG37)
                return false.Label($"minStem {h.StemLength + 1} more stable: {bh} vs {h} for {s}");

            var span = FindMostStableHairpin(s.Substring(h.StemStart, h.StemEnd - h.StemStart + 1));
            if (span is null || span.Value.DeltaG37 != h.DeltaG37 || span.Value.StemLength != h.StemLength)
                return false.Label($"closing span {s.Substring(h.StemStart, h.StemEnd - h.StemStart + 1)}: {span} vs {h}");

            int i = ra % s.Length, j = rb % s.Length;
            if (i > j) (i, j) = (j, i);
            var sub = FindMostStableHairpin(s.Substring(i, j - i + 1));
            return (sub is null || sub.Value.DeltaG37 >= h.DeltaG37)
                .Label($"substring [{i},{j}] more stable: {sub} vs {h} for {s}");
        });
    }

    #endregion

    #region F61 / F64 — HasHairpinPotential on both sides of the 100-nt suffix-tree switch

    /// <summary>
    /// F61/F64: for a &lt; 100-nt sequence over A/C/G/T/U/N/S in both cases (with engineered stems), the answer equals
    /// that of the same sequence padded with N to ≥ 100 nt (N never pairs, so padding moves the call onto the
    /// suffix-tree path without adding or removing a stem), and both paths are case-insensitive.
    /// </summary>
    [FsCheck.NUnit.Property(MaxTest = 200)]
    public Property HasHairpinPotential_ShortScanEqualsSuffixTreePath_WithNonAcgtAndLowercase()
    {
        var alphabet = Gen.Elements("ACGTUNSacgtuns".ToCharArray());
        var gen = from n in Gen.Choose(0, 99)
                  from chars in alphabet.ArrayOf(n)
                  from engineered in Gen.Elements(false, true)
                  from arm in AcgtGen(2, 7)
                  from loop in AcgtGen(0, 6)
                  from at in Gen.Choose(0, 1000)
                  from stem in Gen.Choose(1, 6)
                  from minLoop in Gen.Choose(0, 5)
                  from padTotal in Gen.Choose(0, 40)
                  from padLeft in Gen.Choose(0, 1000)
                  select (Build(new string(chars), engineered, arm, loop, at), stem, minLoop, padTotal, padLeft);
        return Prop.ForAll(gen.ToArbitrary(), x =>
        {
            var (core, stem, minLoop, padExtra, padLeftSeed) = x;
            int pad = Math.Max(0, 100 - core.Length) + padExtra;
            int left = padLeftSeed % (pad + 1);
            string padded = new string('N', left) + core + new string('N', pad - left);
            bool shortPath = HasHairpinPotential(core, stem, minLoop);
            bool suffixPath = HasHairpinPotential(padded, stem, minLoop);
            bool shortLower = HasHairpinPotential(core.ToLowerInvariant(), stem, minLoop);
            bool shortUpper = HasHairpinPotential(core.ToUpperInvariant(), stem, minLoop);
            bool suffixLower = HasHairpinPotential(padded.ToLowerInvariant(), stem, minLoop);
            return (shortPath == suffixPath && shortLower == shortPath && shortUpper == shortPath
                    && suffixLower == suffixPath)
                .Label($"core '{core}' (stem {stem}, loop {minLoop}): short {shortPath}, suffix {suffixPath}, " +
                       $"lower {shortLower}/{suffixLower}, upper {shortUpper}");
        });
    }

    private static string Build(string background, bool engineered, string arm, string loop, int at)
    {
        if (!engineered)
            return background;
        string insert = arm + loop + RevComp(arm);
        string s = background.Insert(at % (background.Length + 1), insert);
        return s.Length < 100 ? s : s.Substring(0, 99);
    }

    #endregion

    #region F55 — THAL_MAX_ALIGN overloads equal the default for ≤ 60-nt strands

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

    private static bool SameStructure<T>(T? t1, IReadOnlyList<string>? l1, T? t2, IReadOnlyList<string>? l2)
        where T : struct =>
        Equals(t1, t2) && (l1 is null ? l2 is null : l2 is not null && l1.SequenceEqual(l2));

    #endregion

    #region F62 — CalculateMeltingTemperaturePrimer3 max_nn_length

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

    #endregion

    #region F68a — NN shift beyond the overlap ≡ shift ±(len+1)

    /// <summary>
    /// F68a: every Tm_NN shift past the other strand leaves the same strands after Biopython drops the over-dangling
    /// ends, so the outcome (value, or exception type and message) for any shift ≥ len(c)+1 equals that for len(c)+1,
    /// and for any shift ≤ −(len+1) that for −(len+1) — including int.MinValue / int.MaxValue — and a failure is
    /// always an <see cref="ArgumentException"/> (never a raw allocation / overflow error).
    /// </summary>
    [FsCheck.NUnit.Property(MaxTest = 150)]
    public Property NearestNeighborTm_ShiftBeyondOverlap_EqualsBoundaryShift()
    {
        var gen = from s in AcgtGen(2, 20)
                  from c in Gen.Frequency((1, Gen.Constant<string?>(null)), (2, AcgtGen(2, 20).Select(v => (string?)v)))
                  from positive in Gen.Elements(true, false)
                  from extra in Gen.Elements(0, 1, 2, 5, 40, 1000, int.MaxValue)
                  from strict in Gen.Elements(true, false)
                  select (s, c, positive, extra, strict);
        return Prop.ForAll(gen.ToArbitrary(), x =>
        {
            var (s, c, positive, extra, strict) = x;
            int cLen = string.IsNullOrEmpty(c) ? s.Length : c.Length;
            long boundary = positive ? cLen + 1L : -(s.Length + 1L);
            long far = positive ? Math.Min(int.MaxValue, boundary + extra) : Math.Max(int.MinValue, boundary - extra);
            var atBoundary = Outcome(() => ThermoConstants.CalculateNearestNeighborTm(s, c, (int)boundary, strict: strict));
            var atFar = Outcome(() => ThermoConstants.CalculateNearestNeighborTm(s, c, (int)far, strict: strict));
            bool argOnly = (atBoundary.Error is null || atBoundary.Error is ArgumentException)
                           && (atFar.Error is null || atFar.Error is ArgumentException);
            bool same = atBoundary.Error is null
                ? atFar.Error is null && atFar.Value.Equals(atBoundary.Value)
                : atFar.Error is not null && atFar.Error.GetType() == atBoundary.Error.GetType()
                  && atFar.Error.Message == atBoundary.Error.Message;
            return (argOnly && same).Label($"{s}/{c} shift {boundary} → {atBoundary}, shift {far} → {atFar} (strict {strict})");
        });
    }

    private static (double Value, Exception? Error) Outcome(Func<double> f)
    {
        try { return (f(), null); }
        catch (Exception e) { return (double.NaN, e); }
    }

    #endregion

    #region F72 — self-complementarity

    /// <summary>
    /// F72: the NN thermodynamics' self-complementarity flag (now <c>NtthalDimer.IsSymmetric</c>, Primer3 <c>symmetry</c> /
    /// Biopython <c>selfcomp</c>) is exactly "the oligo equals its reverse complement", recomputed independently, on
    /// random oligos and on engineered palindromes s + revcomp(s).
    /// </summary>
    [FsCheck.NUnit.Property(MaxTest = 300)]
    public Property NearestNeighborThermodynamics_SelfComplementaryFlag_EqualsReverseComplementIdentity()
    {
        var gen = Gen.Frequency(
            (1, AcgtGen(2, 30)),
            (1, from h in AcgtGen(1, 15) select h + RevComp(h)));
        return Prop.ForAll(gen.ToArbitrary(), s =>
        {
            var r = CalculateNearestNeighborThermodynamics(s);
            bool expected = s == RevComp(s);
            return (r is { } v && v.IsSelfComplementary == expected)
                .Label($"{s}: flag {r?.IsSelfComplementary}, expected {expected}");
        });
    }

    #endregion
}
