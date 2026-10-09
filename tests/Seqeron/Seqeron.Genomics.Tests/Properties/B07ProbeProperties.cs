using FsCheck;
using FsCheck.Fluent;
using static Seqeron.Genomics.MolTools.ProbeDesigner;

namespace Seqeron.Genomics.Tests.Properties;

/// <summary>
/// Heavy-tier property tests for the ProbeDesigner behaviour introduced by review batch B07
/// (docs/Validation/review-2026-09/B07.md): F48 (G+C over the non-N bases), F50 (preset Tm scales), F51 (Primer3
/// internal-oligo penalty ranking), F52 (Kane-criteria hit count), F54 (BLAST+ K indexing), F55/F56 (THAL_MAX_ALIGN
/// opt-in: lazy screen and site duplex Tm), F60 (both-strand reference scan), F61/F63 (fallback stem-loop screen
/// reported, not deciding validity), F65 (no length adjustment for untabulated ungapped schemes) and F67 (tiling
/// with the end-anchored window). Oracles are independent re-computations (Hamming scans, G+C counts, interval
/// unions, |Tm − opt| + |len − opt|, reverse complements) or the documented equalities between public methods —
/// never copies of the implementation. Targets are short (≤ 120 nt) and the ntthal screen is used only where the
/// property concerns it (otherwise <see cref="ProbeStructureScreen.Heuristic"/>), so each case is fast.
///
/// Test Units: PROBE-DESIGN-001, PROBE-VALID-001, PROBE-EVALUE-001.
/// </summary>
[TestFixture]
[Category("Property")]
[Category("MolTools")]
public class B07ProbeProperties
{
    #region Generators and independent helpers

    private static Gen<string> SeqGen(int minLen, int maxLen, params char[] alphabet) =>
        from n in Gen.Choose(minLen, maxLen)
        from chars in Gen.Elements(alphabet.Length == 0 ? ['A', 'C', 'G', 'T'] : alphabet).ArrayOf(n)
        select new string(chars);

    private static Gen<string> AcgtGen(int minLen, int maxLen) => SeqGen(minLen, maxLen, 'A', 'C', 'G', 'T');

    /// <summary>Independent Watson–Crick reverse complement (upper-case ACGT; anything else → N).</summary>
    private static string RevComp(string s)
    {
        var r = new char[s.Length];
        for (int i = 0; i < s.Length; i++)
            r[s.Length - 1 - i] = s[i] switch { 'A' => 'T', 'T' => 'A', 'C' => 'G', 'G' => 'C', _ => 'N' };
        return new string(r);
    }

    private static string Mutate(string s, int[] positions, char[] bases)
    {
        var c = s.ToCharArray();
        for (int i = 0; i < positions.Length && c.Length > 0; i++)
            c[positions[i] % c.Length] = bases[i];
        return new string(c);
    }

    private static bool Close(double a, double b, double rel = 1e-9) =>
        a == b || Math.Abs(a - b) <= rel * Math.Max(1.0, Math.Max(Math.Abs(a), Math.Abs(b)));

    /// <summary>Base parameters with a wide Tm window and the cheap fallback screen (properties not about the screen).</summary>
    private static ProbeParameters Cheap(int minLen, int maxLen, double minTm = 0, double maxTm = 200) =>
        new ProbeParameters(minLen, maxLen, minTm, maxTm, 0.30, 0.70, 5, false, 0.5)
        {
            StructureScreen = ProbeStructureScreen.Heuristic,
        };

    /// <summary>Independent ungapped site oracle: (position, mismatch positions) of every window with ≤ k mismatches.</summary>
    private static List<(int Pos, List<int> Mm)> HammingSites(string reference, string pattern, int k)
    {
        var sites = new List<(int, List<int>)>();
        for (int p = 0; p + pattern.Length <= reference.Length; p++)
        {
            var mm = new List<int>();
            for (int i = 0; i < pattern.Length && mm.Count <= k; i++)
                if (reference[p + i] != pattern[i])
                    mm.Add(i);
            if (mm.Count <= k)
                sites.Add((p, mm));
        }
        return sites;
    }

    /// <summary>Kane et al. (2000) on an ungapped diagonal: identity (L − d)/L &gt; thr, or an identical run &gt; contig.</summary>
    private static bool Kane(List<int> mm, int length, double thr, int contig)
    {
        if ((double)(length - mm.Count) / length > thr)
            return true;
        int run = 0, best = 0, next = 0;
        for (int i = 0; i < length; i++)
        {
            if (next < mm.Count && mm[next] == i) { next++; run = 0; }
            else best = Math.Max(best, ++run);
        }
        return best > contig;
    }

    /// <summary>Independent (sites, Kane sites) count of a probe over references, optionally on both strands (union of positions).</summary>
    private static (int Sites, int Kane) SiteOracle(string probe, IEnumerable<string> refs, int k, bool both, double thr = 0.75, int contig = 15)
    {
        int sites = 0, kane = 0;
        foreach (var r in refs)
        {
            var map = new Dictionary<int, bool>();
            var patterns = both && RevComp(probe) != probe ? new[] { probe, RevComp(probe) } : new[] { probe };
            foreach (var pat in patterns)
                foreach (var (pos, mm) in HammingSites(r, pat, k))
                    map[pos] = (map.TryGetValue(pos, out bool v) && v) || Kane(mm, probe.Length, thr, contig);
            sites += map.Count;
            kane += map.Values.Count(v => v);
        }
        return (sites, kane);
    }

    /// <summary>A probe plus references holding planted copies (0–5 substitutions, either strand) in random flanks.</summary>
    private static Gen<(string Probe, string[] Refs, int K)> ProbeWithReferencesGen(int minLen, int maxLen) =>
        from probe in AcgtGen(minLen, maxLen)
        from nRefs in Gen.Choose(1, 3)
        from refs in (
            from copies in Gen.Choose(0, 3)
            from planted in (
                from d in Gen.Choose(0, 5)
                from pos in Gen.Choose(0, 1000).ArrayOf(d)
                from bases in Gen.Elements('A', 'C', 'G', 'T').ArrayOf(d)
                from rc in Gen.Elements(false, false, true)
                let m = Mutate(probe, pos, bases)
                select rc ? RevComp(m) : m).ArrayOf(copies)
            from flanks in AcgtGen(0, 25).ArrayOf(copies + 1)
            select string.Concat(flanks.Zip(planted.Append(""), (f, p) => f + p))).ArrayOf(nRefs)
        from k in Gen.Choose(0, 4)
        select (probe, refs, k);

    #endregion

    #region F48 — G+C fraction over the non-N bases on ACGTN targets

    /// <summary>
    /// F48: every designed probe's <c>GcContent</c> equals G+C over the A/C/G/T bases of the window (Primer3
    /// <c>gc_and_n_content</c>, N excluded; 0 for an all-N window) and <c>CalculateGcFractionFast(probe.Sequence)</c>,
    /// on random ACGTN targets (the prefix-sum path), both screens.
    /// </summary>
    [FsCheck.NUnit.Property(MaxTest = 100)]
    public Property DesignProbes_GcContent_EqualsGcOverNonNBases_OnAcgtnTargets()
    {
        var gen = from target in SeqGen(20, 70, 'A', 'C', 'G', 'T', 'N', 'a', 'c', 'g', 't', 'n')
                  from minLen in Gen.Choose(10, 18)
                  from extra in Gen.Choose(0, 4)
                  from avoid in Gen.Elements(false, true)
                  select (target, Cheap(minLen, minLen + extra) with { AvoidSecondaryStructure = avoid });
        return Prop.ForAll(gen.ToArbitrary(), x =>
        {
            var (target, param) = x;
            foreach (var p in DesignProbes(target, param, 100_000))
            {
                int gc = p.Sequence.Count(c => c is 'G' or 'C');
                int valid = p.Sequence.Count(c => c is 'A' or 'C' or 'G' or 'T');
                double oracle = valid == 0 ? 0 : (double)gc / valid;
                if (p.Sequence != target.ToUpperInvariant().Substring(p.Start, p.Sequence.Length))
                    return false.Label($"window mismatch at {p.Start}");
                if (!Close(p.GcContent, oracle) || !Close(p.GcContent, p.Sequence.CalculateGcFractionFast()))
                    return false.Label($"{p.Sequence}: GcContent {p.GcContent:R} vs oracle {oracle:R} / fast {p.Sequence.CalculateGcFractionFast():R}");
            }
            return true.ToProperty();
        });
    }

    #endregion

    #region F50 — preset Tm scales

    /// <summary>
    /// F50: for the long-probe presets (FISH / NorthernBlot / SouthernBlot) the Tm of every probe in the preset's length
    /// × G+C window lies in the documented attainable <c>long_seq_tm</c> range at the 50 mM Primer3 probe conditions
    /// (it depends only on length and G+C count), the window overlaps the Tm window, and <c>DesignProbes</c> reports
    /// exactly <c>CalculateMeltingTemperaturePrimer3</c> at the preset's conditions; any permutation of the bases keeps
    /// the Tm (long_seq_tm is composition-only).
    /// </summary>
    [FsCheck.NUnit.Property(MaxTest = 40)]
    public Property LongProbePresets_Tm_InDocumentedAttainableRange_AndPermutationInvariant()
    {
        var presets = new (string Name, ProbeParameters P, double Min, double Max)[]
        {
            ("FISH", Defaults.FISH, 71.2529020719779, 85.35290207197791),
            ("NorthernBlot", Defaults.NorthernBlot, 70.30290207197791, 82.5029020719779),
            ("SouthernBlot", Defaults.SouthernBlot, 70.32012061502427, 85.35290207197791),
        };
        var gen = from i in Gen.Choose(0, presets.Length - 1)
                  let pr = presets[i]
                  from len in Gen.Choose(pr.P.MinLength, Math.Min(pr.P.MaxLength, 260))
                  from gcFrac in Gen.Choose((int)Math.Ceiling(pr.P.MinGc * 1000), (int)(pr.P.MaxGc * 1000))
                  from seed in Gen.Choose(0, int.MaxValue)
                  select (pr, len, gcFrac / 1000.0, seed);
        return Prop.ForAll(gen.ToArbitrary(), x =>
        {
            var (pr, len, gcTarget, seed) = x;
            var rng = new Random(seed);
            int gcCount = Math.Clamp((int)Math.Ceiling(gcTarget * len), (int)Math.Ceiling(pr.P.MinGc * len), (int)Math.Floor(pr.P.MaxGc * len));
            var chars = Enumerable.Range(0, len).Select(j => j < gcCount ? (rng.Next(2) == 0 ? 'G' : 'C') : (rng.Next(2) == 0 ? 'A' : 'T')).ToArray();
            string seq = new(chars.OrderBy(_ => rng.Next()).ToArray());
            string perm = new(seq.OrderBy(_ => rng.Next()).ToArray());
            var p = pr.P;
            double tm = PrimerDesigner.CalculateMeltingTemperaturePrimer3(seq, p.DnaConcentrationNanomolar, p.MonovalentMillimolar,
                p.DivalentMillimolar, p.DntpMillimolar, p.MaxNearestNeighborLength);
            double tmPerm = PrimerDesigner.CalculateMeltingTemperaturePrimer3(perm, p.DnaConcentrationNanomolar, p.MonovalentMillimolar,
                p.DivalentMillimolar, p.DntpMillimolar, p.MaxNearestNeighborLength);
            var probe = DesignProbes(seq, p with { MinLength = len, MaxLength = len, StructureScreen = ProbeStructureScreen.Heuristic }, 1).Single();
            return (tm >= pr.Min - 1e-9 && tm <= pr.Max + 1e-9 && Close(tm, tmPerm, 1e-12) && Close(probe.Tm, tm, 1e-12)
                    && pr.Min <= p.MaxTm && pr.Max >= p.MinTm)
                .Label($"{pr.Name} len {len} gc {gcCount}: tm {tm:R} perm {tmPerm:R} design {probe.Tm:R} vs [{pr.Min}, {pr.Max}]");
        });
    }

    /// <summary>
    /// F50: the Microarray preset's Tm is the nearest-neighbour Tm over the whole oligo at OligoArray's 1 M Na⁺ / 1 µM
    /// (max_nn_length 60): <c>DesignProbes</c> reports <c>CalculateMeltingTemperaturePrimer3(seq, 1000, 1000, 0, 0, 60)</c>
    /// for random 50–60-mers in the G+C window, that value exceeds the 50 mM / 50 nM Primer3 value (more salt and more
    /// oligo stabilise the duplex), and a witness in the 82–90 °C window exists among the designs of the window.
    /// </summary>
    [FsCheck.NUnit.Property(MaxTest = 40)]
    public Property MicroarrayPreset_Tm_IsNearestNeighbourAtOligoArrayConditions()
    {
        var gen = from len in Gen.Choose(50, 60)
                  from gc in Gen.Choose(20, 30)
                  from seed in Gen.Choose(0, int.MaxValue)
                  select (len, gc, seed);
        return Prop.ForAll(gen.ToArbitrary(), x =>
        {
            var (len, gcCountRaw, seed) = x;
            var rng = new Random(seed);
            int gcCount = Math.Clamp(gcCountRaw, (int)Math.Ceiling(0.4 * len), (int)Math.Floor(0.6 * len));
            string seq = new(Enumerable.Range(0, len)
                .Select(j => j < gcCount ? (rng.Next(2) == 0 ? 'G' : 'C') : (rng.Next(2) == 0 ? 'A' : 'T'))
                .OrderBy(_ => rng.Next()).ToArray());
            var p = Defaults.Microarray;
            double tm = PrimerDesigner.CalculateMeltingTemperaturePrimer3(seq, 1000, 1000, 0, 0, 60);
            double tm50 = PrimerDesigner.CalculateMeltingTemperaturePrimer3(seq, 50, 50, 0, 0, 60);
            var probe = DesignProbes(seq, p with { MinLength = len, MaxLength = len, StructureScreen = ProbeStructureScreen.Heuristic }, 1).Single();
            bool tmWarning = probe.Warnings.Any(w => w.StartsWith("Tm ", StringComparison.Ordinal));
            bool outside = tm < p.MinTm || tm > p.MaxTm;
            return (Close(probe.Tm, tm, 1e-12) && tm > tm50 && tmWarning == outside)
                .Label($"{seq}: design {probe.Tm:R} vs {tm:R} (50 mM {tm50:R}), warning {tmWarning}");
        });
    }

    #endregion

    #region F51 — Primer3 internal-oligo penalty ranking

    /// <summary>
    /// F51: with <see cref="ProbeRanking.Primer3Penalty"/> the output is a reordering of the additive candidate set
    /// (same windows, same additive <c>Score</c>, same warnings); every penalty is |Tm − OptTm| + |len − OptLength|
    /// (Primer3 default internal-oligo weights) and null with the additive ranking; the order is Primer3
    /// <c>primer_rec_comp</c> (penalty ascending, start descending, length ascending); and the top-k output is the
    /// prefix of the full ranking.
    /// </summary>
    [FsCheck.NUnit.Property(MaxTest = 60)]
    public Property DesignProbes_Primer3Ranking_IsReorderingOfAdditiveSet_PrefixClosed()
    {
        var gen = from target in AcgtGen(30, 80)
                  from minLen in Gen.Choose(14, 20)
                  from extra in Gen.Choose(0, 5)
                  from minTm in Gen.Choose(30, 60)
                  from span in Gen.Choose(0, 25)
                  from optTmOff in Gen.Choose(0, 100)
                  from optLenOff in Gen.Choose(0, 5)
                  from k in Gen.Choose(1, 12)
                  from thermo in Gen.Elements(false, false, false, true)
                  let p = Cheap(minLen, minLen + extra, minTm, minTm + span) with
                  {
                      StructureScreen = thermo ? ProbeStructureScreen.Thermodynamic : ProbeStructureScreen.Heuristic,
                      AvoidSecondaryStructure = thermo,
                      OptTm = minTm + span * optTmOff / 100.0,
                      OptLength = minLen + Math.Min(optLenOff, extra),
                  }
                  select (target, p, k);
        return Prop.ForAll(gen.ToArbitrary(), x =>
        {
            var (target, p, k) = x;
            var additive = DesignProbes(target, p, 100_000).ToList();
            var pp = p with { Ranking = ProbeRanking.Primer3Penalty };
            var ranked = DesignProbes(target, pp, 100_000).ToList();
            var top = DesignProbes(target, pp, k).ToList();

            var key = (Probe q) => (q.Start, q.Sequence);
            var addById = additive.ToDictionary(key);
            if (ranked.Count != additive.Count || ranked.Any(q => !addById.ContainsKey(key(q))))
                return false.Label($"candidate sets differ: {additive.Count} vs {ranked.Count}");
            foreach (var q in ranked)
            {
                var a = addById[key(q)];
                if (q.Score != a.Score || q.Tm != a.Tm || !q.Warnings.SequenceEqual(a.Warnings) || a.Primer3Penalty is not null)
                    return false.Label($"record differs for {q.Start}/{q.Sequence}");
                double oracle = Math.Abs(q.Tm - p.OptTm) + Math.Abs(q.Sequence.Length - p.OptLength);
                if (q.Primer3Penalty is not double pen || !Close(pen, oracle))
                    return false.Label($"penalty {q.Primer3Penalty} vs {oracle} for {q.Sequence}");
            }
            for (int i = 1; i < ranked.Count; i++)
            {
                var (a, b) = (ranked[i - 1], ranked[i]);
                int c = a.Primer3Penalty!.Value.CompareTo(b.Primer3Penalty!.Value);
                if (c == 0) c = b.Start.CompareTo(a.Start);
                if (c == 0) c = a.Sequence.Length.CompareTo(b.Sequence.Length);
                if (c > 0)
                    return false.Label($"order broken at {i}: {a.Primer3Penalty}/{a.Start}/{a.Sequence.Length} before {b.Primer3Penalty}/{b.Start}/{b.Sequence.Length}");
            }
            bool prefix = top.Select(key).SequenceEqual(ranked.Take(k).Select(key));
            return prefix.Label($"top-{k} is not the prefix of the full ranking");
        });
    }

    #endregion

    #region F52 / F60 — Kane-criteria hit count, both strands

    /// <summary>
    /// F52: <c>OffTargetHits</c> and <c>CrossHybridizingHits</c> equal an independent Hamming-scan oracle (sites = windows
    /// with ≤ k mismatches; Kane sites = identity &gt; 0.75 or an identical run &gt; 15); CrossHybridizingHits ≤
    /// OffTargetHits, equal when L ≥ 13 and k ≤ 3 (every such site is &gt; 75 % identical); the off-target issue is raised
    /// iff CrossHybridizingHits &gt; 1; raising the identity threshold never increases the count.
    /// F60: with <c>bothStrands</c> both counts are the oracle's union over both orientations, ≥ the single-strand counts;
    /// reverse-complementing every reference leaves the both-strand counts unchanged; a reverse-palindromic probe gives
    /// equal counts in both modes.
    /// </summary>
    [FsCheck.NUnit.Property(MaxTest = 120)]
    public Property ValidateProbe_KaneHitCounts_MatchOracle_BothStrandsConsistent()
    {
        var gen = from x in ProbeWithReferencesGen(8, 30)
                  from palindrome in Gen.Elements(false, false, false, true)
                  let probe = palindrome ? x.Probe[..(x.Probe.Length / 2)] + RevComp(x.Probe[..(x.Probe.Length / 2)]) : x.Probe
                  select (probe, x.Refs, x.K);
        var cond = Defaults.qPCR with { StructureScreen = ProbeStructureScreen.Heuristic };
        return Prop.ForAll(gen.ToArbitrary(), x =>
        {
            var (probe, refs, k) = x;
            var single = ValidateProbe(probe, refs, k, conditions: cond);
            var both = ValidateProbe(probe, refs, k, conditions: cond, bothStrands: true);
            var bothRc = ValidateProbe(probe, refs.Select(RevComp).ToArray(), k, conditions: cond, bothStrands: true);
            var o1 = SiteOracle(probe, refs, k, false);
            var o2 = SiteOracle(probe, refs, k, true);
            string ctx = $"{probe} k={k} refs=[{string.Join(",", refs)}]";

            if (single.OffTargetHits != o1.Sites || single.CrossHybridizingHits != o1.Kane)
                return false.Label($"single {single.OffTargetHits}/{single.CrossHybridizingHits} vs oracle {o1} {ctx}");
            if (both.OffTargetHits != o2.Sites || both.CrossHybridizingHits != o2.Kane)
                return false.Label($"both {both.OffTargetHits}/{both.CrossHybridizingHits} vs oracle {o2} {ctx}");
            if (single.CrossHybridizingHits > single.OffTargetHits || both.OffTargetHits < single.OffTargetHits
                || both.CrossHybridizingHits < single.CrossHybridizingHits)
                return false.Label($"ordering violated {ctx}");
            if (probe.Length >= 13 && k <= 3 && single.CrossHybridizingHits != single.OffTargetHits)
                return false.Label($"L ≥ 13, k ≤ 3 but {single.CrossHybridizingHits} != {single.OffTargetHits} {ctx}");
            if (bothRc.OffTargetHits != both.OffTargetHits || bothRc.CrossHybridizingHits != both.CrossHybridizingHits)
                return false.Label($"rc(references) changed both-strand counts {bothRc.OffTargetHits}/{bothRc.CrossHybridizingHits} {ctx}");
            if (RevComp(probe) == probe && (both.OffTargetHits != single.OffTargetHits || both.CrossHybridizingHits != single.CrossHybridizingHits))
                return false.Label($"palindrome counts differ {ctx}");
            foreach (var v in new[] { single, both })
            {
                bool issue = v.Issues.Any(i => i.Contains("potential off-target sites", StringComparison.Ordinal));
                if (issue != v.CrossHybridizingHits > 1)
                    return false.Label($"issue {issue} but CrossHybridizingHits {v.CrossHybridizingHits} {ctx}");
            }
            int previous = int.MaxValue;
            foreach (double thr in new[] { 0.0, 0.5, 0.75, 0.85, 0.95, 1.0 })
            {
                var v = ValidateProbe(probe, refs, k, conditions: cond, maxNonTargetIdentity: thr);
                if (v.CrossHybridizingHits > previous || v.CrossHybridizingHits != SiteOracle(probe, refs, k, false, thr).Kane)
                    return false.Label($"threshold {thr}: {v.CrossHybridizingHits} (previous {previous}) {ctx}");
                previous = v.CrossHybridizingHits;
            }
            return true.ToProperty();
        });
    }

    #endregion

    #region F61 / F63 — fallback stem-loop screen reported, not deciding validity

    /// <summary>
    /// F61/F63: on the fallback path (Heuristic screen, or a &gt; 60-nt probe with the default THAL_MAX_ALIGN)
    /// <c>ThermodynamicScreen</c> is false, <c>HasSecondaryStructure == PrimerDesigner.HasHairpinPotential(probe)</c>,
    /// <c>Warnings</c> holds "Potential secondary structure formation" iff HasSecondaryStructure, no issue mentions
    /// secondary structure, and <c>IsValid</c> is the sourced rule alone: ≤ 1 Kane site and Primer3 alignment-mode
    /// self_any / self_end ≤ 12. <c>DesignProbes</c> on the fallback path carries "Potential secondary structure" iff
    /// AvoidSecondaryStructure and HasHairpinPotential(window).
    /// </summary>
    [FsCheck.NUnit.Property(MaxTest = 100)]
    public Property FallbackScreen_StemLoopFlag_EqualsHasHairpinPotential_AndDoesNotDecideValidity()
    {
        var stemLoop = from f5 in AcgtGen(0, 10)
                       from arm in AcgtGen(4, 9)
                       from loop in AcgtGen(3, 8)
                       from f3 in AcgtGen(0, 10)
                       select f5 + arm + loop + RevComp(arm) + f3;
        var gen = from longProbe in Gen.Elements(false, true)
                  from probe in longProbe ? AcgtGen(61, 90) : Gen.Frequency((1, AcgtGen(8, 50)), (1, stemLoop))
                  from lower in Gen.Elements(false, true)
                  from refCopies in Gen.Choose(0, 2)
                  from flank in AcgtGen(0, 20)
                  select (lower ? probe.ToLowerInvariant() : probe, longProbe,
                          new[] { flank + string.Concat(Enumerable.Repeat(probe + "TTTT", refCopies)) });
        return Prop.ForAll(gen.ToArbitrary(), x =>
        {
            var (probe, longProbe, refs) = x;
            var cond = longProbe ? (ProbeParameters?)null : Defaults.qPCR with { StructureScreen = ProbeStructureScreen.Heuristic };
            var v = ValidateProbe(probe, refs, conditions: cond);
            string up = probe.ToUpperInvariant();
            bool hairpin = PrimerDesigner.HasHairpinPotential(up);
            bool warn = v.Warnings.Contains("Potential secondary structure formation");
            bool structureIssue = v.Issues.Any(i => i.Contains("secondary structure", StringComparison.OrdinalIgnoreCase));
            bool expectedValid = v.CrossHybridizingHits <= 1 && v.SelfAny <= 12.0 && v.SelfEnd <= 12.0;
            bool ok = !v.ThermodynamicScreen && v.HairpinTm is null && v.HasSecondaryStructure == hairpin && warn == hairpin
                      && !structureIssue && v.IsValid == expectedValid && v.IsValid == (v.Issues.Count == 0)
                      && v.SelfAny == PrimerDesigner.CalculatePrimerSelfAnyComplementarity(up)
                      && v.SelfEnd == PrimerDesigner.CalculatePrimerSelfEndComplementarity(up);
            if (!ok)
                return false.Label($"{probe}: thermo {v.ThermodynamicScreen} hss {v.HasSecondaryStructure} hairpin {hairpin} warn {warn} " +
                                   $"valid {v.IsValid} expected {expectedValid} issues [{string.Join("; ", v.Issues)}]");

            if (longProbe)
                return true.ToProperty();
            var dp = Cheap(Math.Min(up.Length, 12), Math.Min(up.Length, 14)) with { AvoidSecondaryStructure = true };
            foreach (var p in DesignProbes(up, dp, 100_000))
            {
                bool flag = p.Warnings.Contains("Potential secondary structure");
                if (flag != PrimerDesigner.HasHairpinPotential(p.Sequence))
                    return false.Label($"DesignProbes {p.Sequence}: flag {flag}");
            }
            return true.ToProperty();
        });
    }

    #endregion

    #region F55 / F56 — THAL_MAX_ALIGN opt-in

    /// <summary>
    /// F55: with <c>ThermodynamicScreenMaxLength</c> raised to 100 on a 61–63-nt preset, the lazily screened top-k is the
    /// prefix of the full (exhaustive) ranking, the full ranking is score-descending with ties by (length, start), and each
    /// probe's self-dimer / hairpin warnings are the ntthal decisions of
    /// <c>CalculatePrimer3OligoStructure(seq, …, 100)</c> at the preset conditions (Tm &gt; MaxStructureTm) — the
    /// thermodynamic screen, not the fallback, ran on the &gt; 60-nt windows. Kept small (≤ 16 windows).
    /// </summary>
    [FsCheck.NUnit.Property(MaxTest = 6)]
    public Property DesignProbes_OptInScreenAbove60_LazyEqualsExhaustive_AndUsesNtthal()
    {
        var gen = from target in AcgtGen(64, 67)
                  from k in Gen.Choose(1, 4)
                  select (target, k);
        var p = new ProbeParameters(61, 63, 0, 200, 0.2, 0.8, 8, true, 0.5) { ThermodynamicScreenMaxLength = 100 };
        return Prop.ForAll(gen.ToArbitrary(), x =>
        {
            var (target, k) = x;
            var full = DesignProbes(target, p, 100_000).ToList();
            var top = DesignProbes(target, p, k).ToList();
            if (!top.SequenceEqual(full.Take(k), ProbeComparer.Instance))
                return false.Label("top-k differs from the prefix of the exhaustive ranking");
            for (int i = 1; i < full.Count; i++)
            {
                var (a, b) = (full[i - 1], full[i]);
                bool ordered = a.Score > b.Score || (a.Score == b.Score
                    && (a.Sequence.Length, a.Start).CompareTo((b.Sequence.Length, b.Start)) < 0);
                if (!ordered)
                    return false.Label($"order broken at {i}");
            }
            foreach (var q in full)
            {
                var s = PrimerDesigner.CalculatePrimer3OligoStructure(q.Sequence, p.MonovalentMillimolar, p.DivalentMillimolar,
                    p.DntpMillimolar, p.DnaConcentrationNanomolar, 100);
                if (s is not { } v)
                    return false.Label($"no ntthal structure for {q.Sequence}");
                bool selfDimer = v.SelfAnyTh > p.MaxStructureTm || v.SelfEndTh > p.MaxStructureTm;
                bool hairpin = v.HairpinTh > p.MaxStructureTm;
                if (q.Warnings.Any(w => w.StartsWith("Self-dimer Tm", StringComparison.Ordinal)) != selfDimer
                    || q.Warnings.Any(w => w.StartsWith("Hairpin Tm", StringComparison.Ordinal)) != hairpin
                    || q.Warnings.Any(w => w.StartsWith("Self-complementarity: Primer3", StringComparison.Ordinal)
                                           || w == "Potential secondary structure"))
                    return false.Label($"{q.Sequence}: warnings [{string.Join("; ", q.Warnings)}] vs ntthal {v}");
            }
            return true.ToProperty();
        });
    }

    private sealed class ProbeComparer : IEqualityComparer<Probe>
    {
        public static readonly ProbeComparer Instance = new();
        public bool Equals(Probe a, Probe b) =>
            a.Sequence == b.Sequence && a.Start == b.Start && a.Score == b.Score && a.Tm == b.Tm
            && a.Warnings.SequenceEqual(b.Warnings) && a.Primer3Penalty == b.Primer3Penalty;
        public int GetHashCode(Probe p) => HashCode.Combine(p.Sequence, p.Start);
    }

    /// <summary>
    /// F56: for probes ≤ 40 nt (so min(probe, site) ≤ 60) the site duplex Tm of <c>AssessCrossHybridization</c> equals
    /// the library's calc_heterodimer equivalent <c>CalculateDimerThermodynamicsNtthal(probe, revcomp(site), Any,
    /// mv, dv, dNTP, C, 37 °C, max loop 30)</c> at the conditions (0 when no duplex forms), and raising
    /// <c>ThermodynamicScreenMaxLength</c> (opt-in) changes nothing there.
    /// </summary>
    [FsCheck.NUnit.Property(MaxTest = 60)]
    public Property AssessCrossHybridization_DuplexTm_EqualsHeterodimer_OptInUnchangedUpTo60()
    {
        var gen = from x in ProbeWithReferencesGen(12, 40)
                  from maxAlign in Gen.Elements(61, 100, 500, 10_000)
                  from mv in Gen.Choose(10, 1000)
                  select (x.Probe, x.Refs, maxAlign, (double)mv);
        return Prop.ForAll(gen.ToArbitrary(), x =>
        {
            var (probe, refs, maxAlign, mv) = x;
            var cond = Defaults.qPCR with { MonovalentMillimolar = mv };
            var def = AssessCrossHybridization(probe, refs, conditions: cond);
            var opt = AssessCrossHybridization(probe, refs, conditions: cond with { ThermodynamicScreenMaxLength = maxAlign });
            if (def.Count != opt.Count)
                return false.Label("assessment count differs");
            for (int i = 0; i < def.Count; i++)
            {
                var a = def[i];
                if (a != opt[i])
                    return false.Label($"opt-in {maxAlign} changed assessment {i}: {a} vs {opt[i]}");
                if (a.SiteStart < 0)
                {
                    if (a.DuplexTm is not null)
                        return false.Label("DuplexTm without a site");
                    continue;
                }
                string strand = a.ReverseComplementStrand ? RevComp(refs[a.NonTargetIndex]) : refs[a.NonTargetIndex];
                string site = strand.Substring(a.SiteStart, a.SiteEnd - a.SiteStart + 1);
                var d = PrimerDesigner.CalculateDimerThermodynamicsNtthal(probe, RevComp(site), PrimerDesigner.NtthalAlignmentMode.Any,
                    cond.MonovalentMillimolar / 1000.0, cond.DivalentMillimolar / 1000.0, cond.DntpMillimolar / 1000.0,
                    cond.DnaConcentrationNanomolar * 1e-9, PrimerDesigner.NtthalDefaultTemperatureCelsius, PrimerDesigner.NtthalDefaultMaxLoop);
                double expected = d?.TmCelsius ?? 0.0;
                if (a.DuplexTm is not double tm || !Close(tm, expected, 1e-12))
                    return false.Label($"{probe} site {site}: DuplexTm {a.DuplexTm} vs heterodimer {expected:R}");
            }
            return true.ToProperty();
        });
    }

    #endregion

    #region F54 / F65 — BLAST+ K indexing and ungapped length adjustment

    private static readonly HashSet<(int, int)> BlastnTables =
        [(1, -5), (1, -4), (2, -7), (1, -3), (2, -5), (1, -2), (2, -3), (3, -4), (1, -1), (3, -2), (4, -5), (5, -4)];

    private static int Gcd(int a, int b) => b == 0 ? Math.Abs(a) : Gcd(b, a % b);

    /// <summary>Schemes with a defined λ under the uniform composition: match &gt; 0, mismatch &lt; 0, match + 3·mismatch &lt; 0.</summary>
    private static Gen<(int Match, int Mismatch)> SchemeGen(int maxMatch, int maxMismatch) =>
        (from m in Gen.Choose(1, maxMatch)
         from mm in Gen.Choose(1, maxMismatch)
         select (m, -mm)).Where(s => s.Item1 + 3 * s.Item2 < 0);

    /// <summary>
    /// F54: for every gcd = 1 scheme <see cref="KarlinKMethod.NcbiBlast"/> and <see cref="KarlinKMethod.ReducedLattice"/>
    /// give the same λ / K / H; for the closed-form schemes δ·(a, b) with reduced lowest score −1 or highest +1 they agree
    /// too, and <see cref="KarlinKMethod.ReducedLattice"/> is scale-invariant (K(δa, δb) = K(a, b), λ(δa, δb) = λ(a, b)/δ).
    /// Uniform and non-uniform compositions.
    /// </summary>
    [FsCheck.NUnit.Property(MaxTest = 200)]
    public Property UngappedK_NcbiBlastEqualsReducedLattice_ForCoprimeAndClosedFormSchemes()
    {
        var closedForm = from a in Gen.Choose(1, 6)
                         from b in Gen.Choose(1, 6)
                         from which in Gen.Elements(true, false)
                         let s = which ? (Match: 1, Mismatch: -b) : (Match: a, Mismatch: -1)
                         where s.Match + 3 * s.Mismatch < 0
                         from d in Gen.Choose(2, 6)
                         select (s.Match * d, s.Mismatch * d, d);
        var coprime = from s in SchemeGen(12, 20) where Gcd(s.Match, -s.Mismatch) == 1 select (s.Match, s.Mismatch, 1);
        var gen = from s in Gen.Frequency((1, closedForm), (1, coprime))
                  from pct in Gen.Frequency((2, Gen.Constant(25)), (1, Gen.Choose(15, 35)))
                  let freq = pct / 100.0
                  let pMatch = 4 * freq * freq
                  where pMatch * s.Item1 + (1 - pMatch) * s.Item2 < 0
                  select (s, freq);
        return Prop.ForAll(gen.ToArbitrary(), x =>
        {
            var ((match, mismatch, d), f) = x;
            var r = ComputeUngappedKarlinParameters(match, mismatch, f, KarlinKMethod.ReducedLattice);
            var n = ComputeUngappedKarlinParameters(match, mismatch, f, KarlinKMethod.NcbiBlast);
            bool same = Close(r.Lambda, n.Lambda, 1e-12) && Close(r.K, n.K, 1e-12) && Close(r.H, n.H, 1e-12);
            if (d > 1)
            {
                var red = ComputeUngappedKarlinParameters(match / d, mismatch / d, f, KarlinKMethod.ReducedLattice);
                same &= Close(r.K, red.K, 1e-9) && Close(r.Lambda * d, red.Lambda, 1e-9);
            }
            return same.Label($"{match}/{mismatch} δ {d} p {f}: reduced {r} ncbi {n}");
        });
    }

    /// <summary>
    /// F65: ungapped <c>ComputeBlastnStatistics</c> of a scheme whose gcd-reduced form has no BLAST+ <c>blastn_values_*</c>
    /// table applies no length adjustment (α = β = 0, ℓ = 0, search space m·n); tabulated schemes keep α = λ/H and
    /// ℓ = <c>ComputeLengthAdjustment(K, α/λ, β, m, n, N)</c> with search space (m − ℓ)·max(1, n − N·ℓ). λ / K / H are the
    /// NcbiBlast ungapped values either way.
    /// </summary>
    [FsCheck.NUnit.Property(MaxTest = 200)]
    public Property UngappedBlastnStatistics_LengthAdjustmentOnlyForTabulatedSchemes()
    {
        var gen = from s in SchemeGen(8, 12)
                  from m in Gen.Choose(10, 500)
                  from n in Gen.Choose(10, 2_000_000)
                  from bigN in Gen.Choose(1, 20)
                  from raw in Gen.Choose(0, 400)
                  select (s.Match, s.Mismatch, m, (long)n, bigN, raw);
        return Prop.ForAll(gen.ToArbitrary(), x =>
        {
            var (match, mismatch, m, n, bigN, raw) = x;
            int g = Gcd(match, -mismatch);
            bool tabulated = BlastnTables.Contains((match / g, mismatch / g));
            var st = ComputeBlastnStatistics(raw, m, n, bigN, new ScoringMatrix(match, mismatch, -5, -2), gapped: false);
            var u = ComputeUngappedKarlinParameters(match, mismatch, 0.25, KarlinKMethod.NcbiBlast);
            var p = st.Parameters;
            bool lambdaKh = p.Lambda == u.Lambda && p.K == u.K && p.H == u.H && !p.Gapped;
            if (!tabulated)
                return (lambdaKh && p.Alpha == 0 && p.Beta == 0 && st.LengthAdjustment == 0
                        && st.EffectiveSearchSpace == (double)m * n)
                    .Label($"untabulated {match}/{mismatch}: α {p.Alpha} β {p.Beta} ℓ {st.LengthAdjustment} space {st.EffectiveSearchSpace}");
            int ell = ComputeLengthAdjustment(p.K, p.Alpha / p.Lambda, p.Beta, m, n, bigN);
            double space = (m - ell) * Math.Max(1.0, n - (double)bigN * ell);
            return (lambdaKh && Close(p.Alpha, u.Lambda / u.H, 1e-12) && p.Beta == u.Beta && st.LengthAdjustment == ell
                    && Close(st.EffectiveSearchSpace, space, 1e-12))
                .Label($"tabulated {match}/{mismatch}: α {p.Alpha} vs {u.Lambda / u.H}, ℓ {st.LengthAdjustment} vs {ell}");
        });
    }

    #endregion

    #region F67 — tiling with the end-anchored window

    /// <summary>
    /// F67: <c>DesignTilingProbes</c> windows are the target's substrings at strictly increasing starts 0, step, 2·step, …
    /// (step = P − overlap) while they fit, plus one end-anchored window at L − P exactly when (L − P) mod step ≠ 0; the
    /// last start is always L − P; <c>Coverage</c> equals an independent interval-union count, = L whenever overlap ≥ 0;
    /// every window is P nt with End = Start + P − 1 and type Tiling; MeanTm / TmRange summarise the windows' Tm.
    /// </summary>
    [FsCheck.NUnit.Property(MaxTest = 150)]
    public Property DesignTilingProbes_GridPlusEndAnchor_CoverageEqualsIntervalUnion()
    {
        var gen = from target in AcgtGen(1, 260)
                  from pRaw in Gen.Choose(1, 80)
                  let probeLength = Math.Min(pRaw, target.Length)
                  from overlap in Gen.Frequency((3, Gen.Choose(0, probeLength - 1)), (1, Gen.Choose(-probeLength, -1)))
                  select (target, probeLength, overlap);
        return Prop.ForAll(gen.ToArbitrary(), x =>
        {
            var (target, P, overlap) = x;
            var param = Cheap(P, P);
            var set = DesignTilingProbes(target, P, overlap, param);
            int L = target.Length, step = P - overlap;
            var starts = set.Probes.Select(p => p.Start).ToList();
            var expected = new List<int>();
            for (long s = 0; s <= L - P; s += step)
                expected.Add((int)s);
            if ((L - P) % step != 0)
                expected.Add(L - P);
            var covered = new bool[L];
            foreach (var p in set.Probes)
                for (int i = p.Start; i <= p.End; i++)
                    covered[i] = true;
            int coverage = covered.Count(c => c);
            bool windows = set.Probes.All(p => p.Sequence == target.Substring(p.Start, P) && p.End == p.Start + P - 1
                                               && p.Type == ProbeType.Tiling);
            bool distinctSorted = starts.Zip(starts.Skip(1), (a, b) => a < b).All(t => t);
            bool spacing = Enumerable.Range(1, Math.Max(0, starts.Count - 2)).All(i => starts[i] - starts[i - 1] == step)
                           && (starts.Count < 2 || starts[^1] - starts[^2] == step || (L - P) % step != 0);
            bool summary = Close(set.MeanTm, set.Probes.Average(p => p.Tm))
                           && Close(set.TmRange, set.Probes.Max(p => p.Tm) - set.Probes.Min(p => p.Tm));
            bool ok = starts.SequenceEqual(expected) && starts[^1] == L - P && distinctSorted && spacing && windows && summary
                      && set.Coverage == coverage && (overlap < 0 || set.Coverage == L);
            return ok.Label($"L {L} P {P} overlap {overlap}: starts [{string.Join(",", starts)}] expected [{string.Join(",", expected)}] " +
                            $"coverage {set.Coverage} vs {coverage}");
        });
    }

    #endregion
}
