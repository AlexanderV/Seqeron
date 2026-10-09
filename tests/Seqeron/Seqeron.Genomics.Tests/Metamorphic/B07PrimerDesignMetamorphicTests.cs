using FsCheck;
using FsCheck.Fluent;
using static Seqeron.Genomics.MolTools.PrimerDesigner;

namespace Seqeron.Genomics.Tests.Metamorphic;

/// <summary>
/// Heavy-tier metamorphic / property tests for the Primer3 design features added by review batch B07
/// (docs/Validation/review-2026-09/B07.md F32–F47): search region and product-size ranges, PRIMER_NUM_RETURN,
/// PRIMER_PICK_INTERNAL_OLIGO, alignment mode, PRIMER_WT_END_STABILITY, reaction conditions, the 3′-end checks
/// (GC clamp / max end GC / max end stability / min 3′ distance), mispriming library, template mispriming,
/// fraction bound / position penalty, sequence quality, template masking and lowercase masking.
/// Every relation is either a containment (a stricter limit keeps a subset; a limit at the observed maximum changes
/// nothing), a symmetry (reverse-complementing the template swaps left and right) or an independently recomputed
/// contract of the returned pairs. Templates are random 150–210-nt ACGT strings, so each case is fast.
///
/// Test Units: PRIMER-DESIGN-001, PRIMER-TM-001, PRIMER-STRUCT-001.
/// </summary>
[TestFixture]
[Category("Metamorphic")]
[Category("MolTools")]
public class B07PrimerDesignMetamorphicTests
{
    #region Scenario generator and independent helpers

    private sealed record Scenario(string Template, int TargetStart, int TargetEnd, int Seed)
    {
        public override string ToString() => $"seed {Seed}: target [{TargetStart},{TargetEnd}) template {Template}";
    }

    private static string RandomAcgt(Random rng, int n)
    {
        var c = new char[n];
        for (int i = 0; i < n; i++)
            c[i] = "ACGT"[rng.Next(4)];
        return new string(c);
    }

    private static Arbitrary<Scenario> ScenarioArb() =>
        Gen.Choose(0, int.MaxValue - 1).Select(seed =>
        {
            var rng = new Random(seed);
            string t = RandomAcgt(rng, rng.Next(150, 211));
            int ts = t.Length / 2 - rng.Next(3, 10);
            return new Scenario(t, ts, ts + rng.Next(5, 20), seed);
        }).ToArbitrary();

    private static readonly PrimerParameters P3 = Primer3DefaultParameters;

    private static PrimerPairOptions Opts(int numReturn = 5) => PrimerPairOptions.Primer3Defaults with
    {
        ProductSizeRanges = [new ProductSizeRange(60, 200)],
        NumReturn = numReturn,
    };

    private static IReadOnlyList<PrimerPairResult> Design(Scenario s, PrimerParameters p, PrimerPairOptions o) =>
        DesignPrimerPairs(new DnaSequence(s.Template), s.TargetStart, s.TargetEnd, p, o);

    private static string Key(PrimerPairResult r) =>
        $"{r.Forward!.Position}+{r.Forward.Length}/{r.Reverse!.Position}+{r.Reverse.Length}:{r.PairPenalty:R}";

    private static string Keys(IEnumerable<PrimerPairResult> rs) => string.Join(" | ", rs.Select(Key));

    private static int Left3(PrimerCandidate f) => f.Position + f.Length - 1;

    private static string RevComp(string s)
    {
        var r = new char[s.Length];
        for (int i = 0; i < s.Length; i++)
            r[s.Length - 1 - i] = char.ToUpperInvariant(s[i]) switch { 'A' => 'T', 'T' => 'A', 'C' => 'G', 'G' => 'C', _ => 'N' };
        return new string(r);
    }

    #endregion

    #region F32 — PRIMER_NUM_RETURN, product size, included region, pair contract

    /// <summary>
    /// F32: Primer3's <c>choose_pair_or_triple</c> picks pairs one at a time, so the pairs for PRIMER_NUM_RETURN = k are
    /// the first k of those for any larger value, and <c>DesignPrimers</c> returns the first one.
    /// </summary>
    [FsCheck.NUnit.Property(MaxTest = 15)]
    public Property NumReturn_SmallerValue_IsPrefixOfLarger()
    {
        return Prop.ForAll(ScenarioArb(), s =>
        {
            var many = Design(s, P3, Opts(6));
            var few = Design(s, P3, Opts(2));
            var single = DesignPrimers(new DnaSequence(s.Template), s.TargetStart, s.TargetEnd, P3, Opts(6));
            bool prefix = Keys(few) == Keys(many.Take(2));
            bool first = many.Count == 0 ? !single.IsValid : single.IsValid && Key(single) == Key(many[0]);
            return (prefix && first).Label($"{s}: 6 → {Keys(many)}; 2 → {Keys(few)}; single {single.Message}")
                .Classify(many.Count == 0, "no pair").Classify(many.Count >= 2, "≥ 2 pairs");
        });
    }

    /// <summary>
    /// F32 contract of every returned pair, recomputed independently from the template: valid; product size inside one of
    /// the PRIMER_PRODUCT_SIZE_RANGE ranges and = right end − left start; both primers inside SEQUENCE_INCLUDED_REGION and
    /// outside the target; forward = template substring, reverse = reverse complement of its site; |ΔTm| (unrounded
    /// Primer3 Tm) ≤ PRIMER_PAIR_MAX_DIFF_TM; pair compl_any_th / compl_end_th ≤ 47 °C; pairs distinct.
    /// </summary>
    [FsCheck.NUnit.Property(MaxTest = 20)]
    public Property ReturnedPairs_SatisfyRangeRegionAndPairContract()
    {
        var gen = from s in ScenarioArb().Generator
                  from lo in Gen.Choose(50, 90)
                  from span in Gen.Choose(10, 60)
                  from incStart in Gen.Choose(0, 30)
                  from incTrim in Gen.Choose(0, 30)
                  from diffTm in Gen.Elements(2.0, 5.0, 100.0)
                  select (s, lo, span, incStart, incTrim, diffTm);
        return Prop.ForAll(gen.ToArbitrary(), x =>
        {
            var (s, lo, span, incStart, incTrim, diffTm) = x;
            var ranges = new[] { new ProductSizeRange(lo, lo + span), new ProductSizeRange(lo + span + 1, 200) };
            int incLen = s.Template.Length - incStart - incTrim;
            var opts = Opts(5) with { ProductSizeRanges = ranges, IncludedRegion = (incStart, incLen), MaxTmDifference = diffTm };
            var pairs = Design(s, P3, opts);
            foreach (var r in pairs)
            {
                var f = r.Forward!;
                var v = r.Reverse!;
                int product = v.Position + v.Length - f.Position;
                double dTm = Math.Abs(CalculateMeltingTemperaturePrimer3(f.Sequence) - CalculateMeltingTemperaturePrimer3(v.Sequence));
                bool ok = r.IsValid && r.ProductSize == product
                          && ranges.Any(g => product >= g.Min && product <= g.Max)
                          && f.Position >= incStart && v.Position + v.Length <= incStart + incLen
                          && Left3(f) < s.TargetStart && v.Position >= s.TargetEnd
                          && f.Sequence == s.Template.Substring(f.Position, f.Length)
                          && v.Sequence == RevComp(s.Template.Substring(v.Position, v.Length))
                          && dTm <= diffTm + 1e-9
                          && r.ComplAnyTh <= Primer3MaxStructureTm && r.ComplEndTh <= Primer3MaxStructureTm;
                if (!ok)
                    return false.Label($"{s}: ranges {lo}-{lo + span}, inc ({incStart},{incLen}), ΔTm≤{diffTm}: bad pair {Key(r)} dTm {dTm}");
            }
            bool distinct = pairs.Select(Key).Distinct().Count() == pairs.Count;
            return distinct.Label($"{s}: duplicate pairs {Keys(pairs)}");
        });
    }

    #endregion

    #region F40 — 3′-end checks and PRIMER_MIN_*_THREE_PRIME_DISTANCE

    /// <summary>
    /// F40: a stricter PRIMER_GC_CLAMP / PRIMER_MAX_END_GC / PRIMER_MAX_END_STABILITY keeps a subset of the valid
    /// candidates, and every valid candidate satisfies the check recomputed independently (last c bases G/C; ≤ e G/C
    /// among the last five; end stability −ΔG(3′ pentamer) = −<c>Calculate3PrimeStability</c> ≤ s).
    /// </summary>
    [FsCheck.NUnit.Property(MaxTest = 8)]
    public Property ThreePrimeEndChecks_StricterLimit_KeepsSubsetOfValidCandidates()
    {
        var gen = from s in ScenarioArb().Generator
                  from forward in Gen.Elements(true, false)
                  from clamp in Gen.Choose(0, 3)
                  from endGc in Gen.Choose(1, 4)
                  from stab in Gen.Choose(30, 90)
                  select (s, forward, clamp, endGc, stab / 10.0);
        return Prop.ForAll(gen.ToArbitrary(), x =>
        {
            var (s, forward, clamp, endGc, stab) = x;
            var dna = new DnaSequence(s.Template);
            int start = forward ? 0 : s.Template.Length - 45;
            HashSet<string> Valid(PrimerParameters p) =>
                GeneratePrimerCandidates(dna, start, start + 45, forward, p)
                    .Where(c => c.IsValid).Select(c => $"{c.Position}+{c.Length}").ToHashSet();
            var all = GeneratePrimerCandidates(dna, start, start + 45, forward, P3)
                .Where(c => c.IsValid).ToDictionary(c => $"{c.Position}+{c.Length}", c => c.Sequence);

            var loose = Valid(P3);
            var clampLoose = Valid(P3 with { GcClamp = clamp });
            var clampStrict = Valid(P3 with { GcClamp = clamp + 1 });
            var endLoose = Valid(P3 with { MaxEndGc = endGc });
            var endStrict = Valid(P3 with { MaxEndGc = endGc - 1 });
            var stabLoose = Valid(P3 with { MaxEndStability = stab });
            var stabStrict = Valid(P3 with { MaxEndStability = stab - 1.0 });

            bool subsets = clampStrict.IsSubsetOf(clampLoose) && clampLoose.IsSubsetOf(loose)
                           && endStrict.IsSubsetOf(endLoose) && endLoose.IsSubsetOf(loose)
                           && stabStrict.IsSubsetOf(stabLoose) && stabLoose.IsSubsetOf(loose);
            bool predicates =
                clampLoose.All(k => all[k][^clamp..].All(b => b is 'G' or 'C'))
                && endLoose.All(k => all[k][^5..].Count(b => b is 'G' or 'C') <= endGc)
                && stabLoose.All(k => -Calculate3PrimeStability(all[k]) <= stab + 1e-9)
                // Each check removes exactly the loose candidates violating it (nothing else changes).
                && loose.Where(k => all[k][^clamp..].All(b => b is 'G' or 'C')).ToHashSet().SetEquals(clampLoose)
                && loose.Where(k => all[k][^5..].Count(b => b is 'G' or 'C') <= endGc).ToHashSet().SetEquals(endLoose);
            return (subsets && predicates).Label(
                $"{s} fwd={forward} clamp {clamp} endGc {endGc} stab {stab}: |loose| {loose.Count}, clamp {clampLoose.Count}/{clampStrict.Count}, " +
                $"endGc {endLoose.Count}/{endStrict.Count}, stab {stabLoose.Count}/{stabStrict.Count}");
        });
    }

    /// <summary>
    /// F40: with PRIMER_MIN_THREE_PRIME_DISTANCE = d ≥ 1 any two returned pairs have left 3′ ends and right 3′ ends at
    /// least d bases apart (<c>left/right_oligo_in_pair_overlaps_used_oligo</c>); with d = 0 no primer is reused.
    /// </summary>
    [FsCheck.NUnit.Property(MaxTest = 15)]
    public Property MinThreePrimeDistance_ReturnedPairsKeepDistance()
    {
        var gen = from s in ScenarioArb().Generator from d in Gen.Choose(0, 8) select (s, d);
        return Prop.ForAll(gen.ToArbitrary(), x =>
        {
            var (s, d) = x;
            var pairs = Design(s, P3, Opts(5) with { MinThreePrimeDistance = d });
            for (int i = 0; i < pairs.Count; i++)
                for (int j = i + 1; j < pairs.Count; j++)
                {
                    var (a, b) = (pairs[i], pairs[j]);
                    bool ok = d == 0
                        ? a.Forward!.Sequence + a.Forward.Position != b.Forward!.Sequence + b.Forward.Position
                          && a.Reverse!.Sequence + a.Reverse.Position != b.Reverse!.Sequence + b.Reverse.Position
                        : Math.Abs(Left3(a.Forward!) - Left3(b.Forward!)) >= d
                          && Math.Abs(a.Reverse!.Position - b.Reverse!.Position) >= d;
                    if (!ok)
                        return false.Label($"{s} d={d}: pairs {Key(a)} and {Key(b)}");
                }
            return true.ToProperty();
        });
    }

    #endregion

    #region F46 — PRIMER_LOWERCASE_MASKING

    /// <summary>
    /// F46: the case-preserving string overload with an all-upper-case template (with or without LowercaseMasking), and
    /// a mixed-case template without LowercaseMasking, give exactly the <see cref="DnaSequence"/> result; with
    /// LowercaseMasking no returned left primer ends on, and no right primer starts on, a lower-case base, and the
    /// masked result is a valid design of the upper-case template under the same options.
    /// </summary>
    [FsCheck.NUnit.Property(MaxTest = 15)]
    public Property LowercaseMasking_UppercaseTemplateUnchanged_MaskedThreePrimeEndsRejected()
    {
        return Prop.ForAll(ScenarioArb(), s =>
        {
            var rng = new Random(s.Seed ^ 0x5eed);
            var mixed = s.Template.ToCharArray();
            for (int run = 0; run < 6; run++)
            {
                int at = rng.Next(mixed.Length), len = rng.Next(3, 15);
                for (int i = at; i < Math.Min(mixed.Length, at + len); i++)
                    mixed[i] = char.ToLowerInvariant(mixed[i]);
            }
            string lower = new(mixed);
            var opts = Opts(5);
            var reference = Keys(Design(s, P3, opts));
            var upperMasked = Keys(DesignPrimerPairs(s.Template, s.TargetStart, s.TargetEnd, P3, opts with { LowercaseMasking = true }));
            var mixedUnmasked = Keys(DesignPrimerPairs(lower, s.TargetStart, s.TargetEnd, P3, opts));
            var masked = DesignPrimerPairs(lower, s.TargetStart, s.TargetEnd, P3, opts with { LowercaseMasking = true });
            bool ends = masked.All(r => char.IsUpper(lower[Left3(r.Forward!)]) && char.IsUpper(lower[r.Reverse!.Position]));
            return (reference == upperMasked && reference == mixedUnmasked && ends)
                .Label($"{s} lower {lower}: ref {reference}; upper+mask {upperMasked}; mixed {mixedUnmasked}; masked {Keys(masked)}");
        });
    }

    #endregion

    #region F36 — Primer3 alignment mode

    /// <summary>
    /// F36: under <see cref="PrimerStructureScreen.Primer3Alignment"/> every returned primer has self_any ≤ 8 and
    /// self_end ≤ 3 (recomputed with the public dpal methods) and every pair compl_any ≤ 8 / compl_end ≤ 3 (Primer3
    /// PRIMER_[PAIR_]MAX_SELF/COMPL_ANY/END defaults).
    /// </summary>
    [FsCheck.NUnit.Property(MaxTest = 15)]
    public Property AlignmentMode_ReturnedPairsRespectDpalLimits()
    {
        var p = P3 with { StructureScreen = PrimerStructureScreen.Primer3Alignment };
        return Prop.ForAll(ScenarioArb(), s =>
        {
            foreach (var r in Design(s, p, Opts(4)))
            {
                var f = r.Forward!;
                var v = r.Reverse!;
                bool ok = CalculatePrimerSelfAnyComplementarity(f.Sequence) <= Primer3MaxSelfAny
                          && CalculatePrimerSelfAnyComplementarity(v.Sequence) <= Primer3MaxSelfAny
                          && CalculatePrimerSelfEndComplementarity(f.Sequence) <= Primer3MaxSelfEnd
                          && CalculatePrimerSelfEndComplementarity(v.Sequence) <= Primer3MaxSelfEnd
                          && CalculatePrimerDimerAnyComplementarity(f.Sequence, v.Sequence) <= Primer3MaxPairComplAny
                          && CalculatePrimerDimerEndComplementarity(f.Sequence, v.Sequence) <= Primer3MaxPairComplEnd
                          && r.ComplAny == CalculatePrimerDimerAnyComplementarity(f.Sequence, v.Sequence)
                          && f.SelfAny == CalculatePrimerSelfAnyComplementarity(f.Sequence);
                if (!ok)
                    return false.Label($"{s}: pair {Key(r)} selfAny {f.SelfAny} complAny {r.ComplAny} complEnd {r.ComplEnd}");
            }
            return true.ToProperty();
        });
    }

    #endregion

    #region F34 — PRIMER_PICK_INTERNAL_OLIGO

    /// <summary>
    /// F34: with PRIMER_PICK_INTERNAL_OLIGO every returned pair carries an internal oligo lying strictly between the
    /// primers (start after the left 3′ end, end before the right primer), equal to its template substring, with a Tm
    /// inside the PRIMER_INTERNAL_MIN/MAX_TM window; the primer pairs are a subset of the pairs designable without it.
    /// </summary>
    [FsCheck.NUnit.Property(MaxTest = 12)]
    public Property InternalOligo_LiesBetweenPrimers()
    {
        return Prop.ForAll(ScenarioArb(), s =>
        {
            var settings = new ProbeDesigner.Primer3ProbeSettings();
            var pairs = Design(s, P3, Opts(3) with { PickInternalOligo = true });
            foreach (var r in pairs)
            {
                if (r.InternalOligo is not { } h)
                    return false.Label($"{s}: pair {Key(r)} without internal oligo");
                bool ok = h.Start > Left3(r.Forward!) && h.Start + h.Length - 1 < r.Reverse!.Position
                          && string.Equals(h.Sequence, s.Template.Substring(h.Start, h.Length), StringComparison.OrdinalIgnoreCase)
                          && h.Tm >= settings.MinTm - 1e-9 && h.Tm <= settings.MaxTm + 1e-9;
                if (!ok)
                    return false.Label($"{s}: pair {Key(r)} internal {h}");
            }
            return true.ToProperty();
        });
    }

    #endregion

    #region F43 — template mispriming

    /// <summary>
    /// F43 (WP3-6 note): reverse-complementing the template turns a left primer at p into a right primer at
    /// L − p − len with the same oligo, and its template-mispriming scores (same strand, other strand) are unchanged —
    /// in alignment mode and in thermodynamic (ntthal) mode.
    /// </summary>
    [FsCheck.NUnit.Property(MaxTest = 40)]
    public Property TemplateMispriming_ReverseComplementTemplate_SwapsLeftAndRight()
    {
        var gen = from s in ScenarioArb().Generator
                  from len in Gen.Choose(15, 25)
                  from p in Gen.Choose(0, 1000)
                  from thermo in Gen.Elements(false, true)
                  select (s, len, p, thermo);
        return Prop.ForAll(gen.ToArbitrary(), x =>
        {
            var (s, len, pSeed, thermo) = x;
            string t = s.Template, rc = RevComp(t);
            int p = pSeed % (t.Length - len + 1);
            var left = CalculateTemplateMispriming(new DnaSequence(t), p, len, true, thermo);
            var right = CalculateTemplateMispriming(new DnaSequence(rc), t.Length - p - len, len, false, thermo);
            var right2 = CalculateTemplateMispriming(new DnaSequence(t), p, len, false, thermo);
            var left2 = CalculateTemplateMispriming(new DnaSequence(rc), t.Length - p - len, len, true, thermo);
            return (left == right && right2 == left2).Label($"{s} p={p} len={len} thermo={thermo}: {left} vs {right}; {right2} vs {left2}");
        });
    }

    /// <summary>
    /// F43 (WP3-6 note): a stricter PRIMER_MAX_TEMPLATE_MISPRIMING[_TH] keeps a subset of the valid primers, so (a) with
    /// the limit set to the largest score among the returned primers of a loose run nothing changes, and (b) with any
    /// lower limit every returned primer's score respects it.
    /// </summary>
    [FsCheck.NUnit.Property(MaxTest = 12)]
    public Property TemplateMispriming_StricterLimit_SubsetOfValidCandidates()
    {
        var gen = from s in ScenarioArb().Generator from thermo in Gen.Elements(false, true) from cut in Gen.Choose(1, 5) select (s, thermo, cut);
        return Prop.ForAll(gen.ToArbitrary(), x =>
        {
            var (s, thermo, cut) = x;
            PrimerParameters WithLimit(double v) => thermo
                ? P3 with { ThermodynamicTemplateAlignment = true, MaxTemplateMisprimingTh = v }
                : P3 with { MaxTemplateMispriming = v };
            var loose = Design(s, WithLimit(thermo ? 1000 : 32767), Opts(4));
            if (loose.Count == 0)
                return true.ToProperty();
            double max = loose.SelectMany(r => new[] { r.Forward!.TemplateMispriming!.Value, r.Reverse!.TemplateMispriming!.Value }).Max();
            var atMax = Design(s, WithLimit(max), Opts(4));
            double lower = Math.Max(0.5, max - cut);
            var strict = Design(s, WithLimit(lower), Opts(4));
            bool sound = strict.All(r => r.Forward!.TemplateMispriming <= lower && r.Reverse!.TemplateMispriming <= lower);
            return (Keys(atMax) == Keys(loose) && sound)
                .Label($"{s} thermo={thermo}: loose {Keys(loose)} max {max}; at max {Keys(atMax)}; strict {lower}: {Keys(strict)}");
        });
    }

    #endregion

    #region F41 — mispriming library

    private static PrimerMisprimingLibrary Library(IEnumerable<(string Name, string Seq)> entries) =>
        new(entries.Select(e => new KeyValuePair<string, string>(e.Name, e.Seq)));

    /// <summary>
    /// F41: the library mispriming score is a maximum over entries, so it is invariant to the entry order and
    /// non-decreasing when one entry's weight is raised.
    /// </summary>
    [FsCheck.NUnit.Property(MaxTest = 60)]
    public Property LibraryMispriming_OrderInvariant_MonotoneInWeight()
    {
        var gen = from seed in Gen.Choose(0, int.MaxValue - 1)
                  from forward in Gen.Elements(true, false)
                  select (seed, forward);
        return Prop.ForAll(gen.ToArbitrary(), x =>
        {
            var rng = new Random(x.seed);
            string primer = RandomAcgt(rng, rng.Next(18, 26));
            int n = rng.Next(1, 6);
            var entries = Enumerable.Range(0, n).Select(i =>
            {
                // Entries containing a mutated primer fragment so scores are non-trivial.
                string frag = rng.Next(2) == 0 ? primer.Substring(rng.Next(0, 8)) : RevComp(primer).Substring(rng.Next(0, 8));
                var c = (RandomAcgt(rng, rng.Next(0, 15)) + frag + RandomAcgt(rng, rng.Next(0, 15))).ToCharArray();
                c[rng.Next(c.Length)] = "ACGT"[rng.Next(4)];
                return (Seq: new string(c), W: 1 + rng.Next(0, 30) / 10.0);
            }).ToList();
            var baseLib = Library(entries.Select((e, i) => (string.Create(System.Globalization.CultureInfo.InvariantCulture, $"e{i}*{e.W:0.0}"), e.Seq)));
            var shuffled = Library(entries.Select((e, i) => (string.Create(System.Globalization.CultureInfo.InvariantCulture, $"e{i}*{e.W:0.0}"), e.Seq)).OrderBy(_ => rng.Next()).ToList());
            int bump = rng.Next(n);
            var heavier = Library(entries.Select((e, i) => (string.Create(System.Globalization.CultureInfo.InvariantCulture, $"e{i}*{(i == bump ? e.W + 2 : e.W):0.0}"), e.Seq)));
            double s0 = CalculateLibraryMispriming(primer, x.forward, baseLib).Score;
            double s1 = CalculateLibraryMispriming(primer, x.forward, shuffled).Score;
            double s2 = CalculateLibraryMispriming(primer, x.forward, heavier).Score;
            return (s0 == s1 && s2 >= s0).Label($"{primer} fwd={x.forward}: base {s0}, shuffled {s1}, heavier {s2}");
        });
    }

    /// <summary>
    /// F41: a stricter PRIMER_MAX_LIBRARY_MISPRIMING keeps a subset of the valid primers: at the largest score of the
    /// loose run's returned primers nothing changes, and any lower limit is respected by every returned primer.
    /// </summary>
    [FsCheck.NUnit.Property(MaxTest = 12)]
    public Property LibraryMispriming_StricterLimit_SubsetOfValidCandidates()
    {
        var gen = from s in ScenarioArb().Generator from cut in Gen.Choose(1, 4) select (s, cut);
        return Prop.ForAll(gen.ToArbitrary(), x =>
        {
            var (s, cut) = x;
            var rng = new Random(s.Seed);
            var lib = Library(Enumerable.Range(0, 4).Select(i =>
            {
                int at = rng.Next(0, s.Template.Length - 40);
                string frag = s.Template.Substring(at, 40);
                return ($"frag{i}", i % 2 == 0 ? frag : RevComp(frag));
            }));
            var opts = Opts(4) with { MaxLibraryMispriming = 1000 };
            PrimerParameters WithLimit(double v) => P3 with { MisprimingLibrary = lib, MaxLibraryMispriming = v };
            var loose = Design(s, WithLimit(1000), opts);
            if (loose.Count == 0)
                return true.ToProperty();
            double max = loose.SelectMany(r => new[] { r.Forward!.LibraryMispriming!.Value, r.Reverse!.LibraryMispriming!.Value }).Max();
            // Primer3 compares with the limit truncated to a C short, so the "no change" limit is ⌈max⌉.
            var atMax = Design(s, WithLimit(Math.Ceiling(max)), opts);
            double lower = Math.Max(1, Math.Floor(max) - cut);
            var strict = Design(s, WithLimit(lower), opts);
            bool sound = strict.All(r => r.Forward!.LibraryMispriming <= lower && r.Reverse!.LibraryMispriming <= lower);
            return (Keys(atMax) == Keys(loose) && sound)
                .Label($"{s}: loose {Keys(loose)} max {max}; at max {Keys(atMax)}; strict {lower}: {Keys(strict)}");
        });
    }

    #endregion

    #region F44 — fraction bound and position penalty

    /// <summary>
    /// F44: Primer3's fraction bound lies in [0, 100] and never increases with the annealing temperature
    /// (K = exp(−ΔG/RT) falls with T for a duplex with ΔH &lt; 0); NaN beyond 36 nt.
    /// </summary>
    [FsCheck.NUnit.Property(MaxTest = 200)]
    public Property FractionBound_InRange_NonIncreasingInAnnealingTemperature()
    {
        var gen = from seed in Gen.Choose(0, int.MaxValue - 1)
                  from len in Gen.Choose(8, 45)
                  from t1 in Gen.Choose(1, 1000)
                  from t2 in Gen.Choose(1, 1000)
                  select (seed, len, Math.Min(t1, t2) / 10.0, Math.Max(t1, t2) / 10.0);
        return Prop.ForAll(gen.ToArbitrary(), x =>
        {
            var (seed, len, lo, hi) = x;
            string p = RandomAcgt(new Random(seed), len);
            double bLo = CalculateFractionBoundPrimer3(p, lo);
            double bHi = CalculateFractionBoundPrimer3(p, hi);
            if (len > Primer3MaxNnTmLength)
                return (double.IsNaN(bLo) && double.IsNaN(bHi)).Label($"{p}: {bLo} {bHi}");
            return (bLo >= 0 && bLo <= 100 && bHi >= 0 && bHi <= 100 && bHi <= bLo + 1e-9)
                .Label($"{p}: bound({lo}) = {bLo:R}, bound({hi}) = {bHi:R}");
        });
    }

    /// <summary>
    /// F44: with non-default PRIMER_INSIDE/OUTSIDE_PENALTY every returned pair spans the target (left 3′ end before the
    /// right 3′ end, left 3′ ≤ target last base, right 3′ ≥ target first base) and each primer's position penalty is
    /// <c>CalculatePositionPenaltyPrimer3</c> of its own coordinates.
    /// </summary>
    [FsCheck.NUnit.Property(MaxTest = 15)]
    public Property PositionPenalties_ReturnedPairsSpanTarget()
    {
        var gen = from s in ScenarioArb().Generator
                  from inside in Gen.Elements(0.0, 0.5, 2.0)
                  from outside in Gen.Elements(0.1, 0.5, 1.0)
                  select (s, inside, outside);
        return Prop.ForAll(gen.ToArbitrary(), x =>
        {
            var (s, inside, outside) = x;
            var pairs = Design(s, P3, Opts(4) with { InsidePenalty = inside, OutsidePenalty = outside });
            foreach (var r in pairs)
            {
                var f = r.Forward!;
                var v = r.Reverse!;
                bool ok = Left3(f) < v.Position && Left3(f) <= s.TargetEnd - 1 && v.Position >= s.TargetStart
                          && f.PositionPenalty == CalculatePositionPenaltyPrimer3(f.Position, f.Length, true, s.TargetStart, s.TargetEnd, inside, outside)
                          && v.PositionPenalty == CalculatePositionPenaltyPrimer3(v.Position, v.Length, false, s.TargetStart, s.TargetEnd, inside, outside);
                if (!ok)
                    return false.Label($"{s} inside {inside} outside {outside}: pair {Key(r)} pos {f.PositionPenalty}/{v.PositionPenalty}");
            }
            return true.ToProperty();
        });
    }

    #endregion

    #region F45 — sequence quality

    /// <summary>
    /// F45: with SEQUENCE_QUALITY every returned primer's minimum base quality (recomputed from the array) is ≥
    /// PRIMER_MIN_QUALITY and equals <see cref="PrimerCandidate.MinSequenceQuality"/>, and the minimum over its five
    /// 3′-most bases (last five of a left primer, first five of a right primer) is ≥ PRIMER_MIN_END_QUALITY.
    /// </summary>
    [FsCheck.NUnit.Property(MaxTest = 12)]
    public Property SequenceQuality_ReturnedPrimersMeetMinimumQualities()
    {
        var gen = from s in ScenarioArb().Generator from minQ in Gen.Choose(10, 30) from minEnd in Gen.Choose(10, 40) select (s, minQ, minEnd);
        return Prop.ForAll(gen.ToArbitrary(), x =>
        {
            var (s, minQ, minEnd) = x;
            var rng = new Random(s.Seed ^ 0x9a1);
            var q = Enumerable.Range(0, s.Template.Length).Select(_ => rng.Next(10) == 0 ? rng.Next(0, 40) : rng.Next(25, 61)).ToArray();
            var pairs = Design(s, P3 with { MinQuality = minQ, MinEndQuality = minEnd }, Opts(4) with { SequenceQuality = q });
            foreach (var r in pairs)
            {
                var f = r.Forward!;
                var v = r.Reverse!;
                int fMin = q.Skip(f.Position).Take(f.Length).Min(), vMin = q.Skip(v.Position).Take(v.Length).Min();
                int fEnd = q.Skip(f.Position + f.Length - 5).Take(5).Min(), vEnd = q.Skip(v.Position).Take(5).Min();
                bool ok = fMin >= minQ && vMin >= minQ && fEnd >= minEnd && vEnd >= minEnd
                          && f.MinSequenceQuality == fMin && v.MinSequenceQuality == vMin;
                if (!ok)
                    return false.Label($"{s} minQ {minQ} minEnd {minEnd}: pair {Key(r)} min {fMin}/{vMin} end {fEnd}/{vEnd}");
            }
            return true.ToProperty();
        });
    }

    #endregion

    #region F38 / F39 — end-stability weight and reaction conditions

    /// <summary>
    /// F38: the PRIMER_WT_END_STABILITY term of <c>p_obj_fn</c> is weight × end_stability, so the penalty is affine in
    /// the weight with slope = the end stability (non-decreasing for end_stability ≥ 0).
    /// </summary>
    [FsCheck.NUnit.Property(MaxTest = 200)]
    public Property Penalty_EndStabilityWeight_AffineWithSlopeEndStability()
    {
        var gen = from tm in Gen.Choose(500, 700)
                  from len in Gen.Choose(15, 30)
                  from gc in Gen.Choose(20, 80)
                  from es in Gen.Choose(0, 120)
                  from w1 in Gen.Choose(0, 50)
                  from w2 in Gen.Choose(0, 50)
                  select (tm / 10.0, len, (double)gc, es / 10.0, w1 / 10.0, w2 / 10.0);
        return Prop.ForAll(gen.ToArbitrary(), x =>
        {
            var (tm, len, gc, es, w1, w2) = x;
            var inputs = new Primer3PenaltyInputs(tm, len, gc, EndStability: es);
            double p0 = CalculatePrimer3Penalty(inputs, DefaultPrimer3Weights with { EndStability = 0 });
            double p1 = CalculatePrimer3Penalty(inputs, DefaultPrimer3Weights with { EndStability = w1 });
            double p2 = CalculatePrimer3Penalty(inputs, DefaultPrimer3Weights with { EndStability = w2 });
            bool ok = Math.Abs(p1 - p0 - w1 * es) <= 1e-9 && Math.Abs(p2 - p0 - w2 * es) <= 1e-9
                      && (w1 <= w2 ? p1 <= p2 + 1e-12 : p2 <= p1 + 1e-12);
            return ok.Label($"Tm {tm} len {len} gc {gc} es {es}: p0 {p0:R} p({w1}) {p1:R} p({w2}) {p2:R}");
        });
    }

    /// <summary>
    /// F39: the Primer3 primer Tm (<c>seqtm</c>) never decreases with PRIMER_SALT_MONOVALENT, PRIMER_SALT_DIVALENT or
    /// PRIMER_DNA_CONC and never increases with PRIMER_DNTP_CONC ([Mon]_eq = Mon + 120·√(Mg − dNTP) enters only through
    /// the SantaLucia entropy salt term / <c>long_seq_tm</c>'s 16.6·log10; C_T through R·ln(C_T/x)).
    /// </summary>
    [FsCheck.NUnit.Property(MaxTest = 200)]
    public Property Primer3Tm_MonotoneInReactionConditions()
    {
        var gen = from seed in Gen.Choose(0, int.MaxValue - 1)
                  from len in Gen.Choose(10, 50)
                  from mv in Gen.Choose(1, 400)
                  from dv in Gen.Choose(0, 60)
                  from dntp in Gen.Choose(0, 60)
                  from dna in Gen.Choose(1, 1000)
                  from bump in Gen.Choose(1, 100)
                  select (seed, len, (double)mv, dv / 10.0, dntp / 10.0, (double)dna, (double)bump);
        return Prop.ForAll(gen.ToArbitrary(), x =>
        {
            var (seed, len, mv, dv, dntp, dna, bump) = x;
            string p = RandomAcgt(new Random(seed), len);
            double t = CalculateMeltingTemperaturePrimer3(p, dna, mv, dv, dntp);
            double tMv = CalculateMeltingTemperaturePrimer3(p, dna, mv + bump, dv, dntp);
            double tDv = CalculateMeltingTemperaturePrimer3(p, dna, mv, dv + bump / 10.0, dntp);
            double tDna = CalculateMeltingTemperaturePrimer3(p, dna + bump, mv, dv, dntp);
            double tDntp = CalculateMeltingTemperaturePrimer3(p, dna, mv, dv, dntp + bump / 10.0);
            const double eps = 1e-9;
            bool ok = tMv >= t - eps && tDv >= t - eps && tDna >= t - eps && tDntp <= t + eps;
            return ok.Label($"{p}: base {t:R}, +mv {tMv:R}, +dv {tDv:R}, +dna {tDna:R}, +dntp {tDntp:R}");
        });
    }

    #endregion

    #region F46 (masking part) — PRIMER_MASK_TEMPLATE

    private static PrimerMaskingKmerLists KmerLists(string template, Random rng)
    {
        Dictionary<string, int> Counts(int k)
        {
            var d = new Dictionary<string, int>();
            for (int i = 0; i + k <= template.Length; i += rng.Next(1, 4))
                d[template.Substring(i, k)] = rng.Next(0, 3) == 0 ? rng.Next(0, 10) : rng.Next(100, 100000);
            d.TryAdd(new string('A', k), 1);
            return d;
        }
        return new PrimerMaskingKmerLists(Counts(11), Counts(16));
    }

    /// <summary>
    /// F46 masking: the predicted failure rate lies in [0, 1) and is 0 below 16 nt; with PRIMER_MASK_TEMPLATE no returned
    /// left primer ends on a base masked (lower case) on <c>MaskTemplatePrimer3(…).Forward</c>, no right primer's 3′ base
    /// is masked on <c>.Reverse</c>, and each primer reports its failure rate.
    /// </summary>
    [FsCheck.NUnit.Property(MaxTest = 12)]
    public Property MaskTemplate_ReturnedPrimersAvoidMaskedThreePrimeEnds()
    {
        return Prop.ForAll(ScenarioArb(), s =>
        {
            var rng = new Random(s.Seed ^ 0x77);
            var lists = KmerLists(s.Template, rng);
            for (int i = 0; i < 20; i++)
            {
                string probe = RandomAcgt(rng, rng.Next(5, 30));
                double fr = CalculateMaskFailureRatePrimer3(probe, lists);
                if (!(fr >= 0 && fr < 1) || (probe.Length < 16 && fr != 0))
                    return false.Label($"failure rate {fr} for {probe}");
            }
            var (fwd, rev) = MaskTemplatePrimer3(s.Template, lists);
            var pairs = DesignPrimerPairs(s.Template, s.TargetStart, s.TargetEnd, P3,
                Opts(4) with { MaskTemplate = true, MaskKmerLists = lists });
            foreach (var r in pairs)
            {
                var f = r.Forward!;
                var v = r.Reverse!;
                bool ok = !char.IsLower(fwd[Left3(f)]) && !char.IsLower(rev[v.Position])
                          && f.MaskFailureRate == CalculateMaskFailureRatePrimer3(f.Sequence, lists)
                          && v.MaskFailureRate == CalculateMaskFailureRatePrimer3(v.Sequence, lists);
                if (!ok)
                    return false.Label($"{s}: pair {Key(r)} fwdMask {fwd} revMask {rev}");
            }
            return true.ToProperty();
        });
    }

    #endregion
}
