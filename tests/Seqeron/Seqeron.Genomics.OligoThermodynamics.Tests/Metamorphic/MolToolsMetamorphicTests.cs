// SPDX-License-Identifier: GPL-2.0-or-later
// Tests of the port of Primer3 ntthal (thal.c) and of the Primer3-style primer / probe design; part of Seqeron.Genomics.OligoThermodynamics.Tests, licensed under the GNU General Public License, version 2 or (at your option) any later version (see LICENSE in the production project).
namespace Seqeron.Genomics.Tests.Metamorphic;

/// <summary>
/// Metamorphic tests for the MolTools area.
///
/// Each test encodes a metamorphic relation (MR) — a property relating the outputs
/// of multiple runs under an input transformation, with no hardcoded oracle. The
/// relations are derived from the PAM-site-detection *definition*, not from the
/// current implementation's output: a PAM site exists wherever the strand-local
/// PAM motif occurs and the guide-length target fits within sequence bounds.
///
/// ───────────────────────────────────────────────────────────────────────────
/// Units: CRISPR-PAM-001 (PAM-site finding), CRISPR-GUIDE-001 (guide-RNA design),
///        CRISPR-OFF-001 (off-target scoring), PRIMER-TM-001 (primer melting temperature),
///        PRIMER-DESIGN-001 (PCR primer design),
///        PRIMER-STRUCT-001 (primer secondary structure: self-dimer / hairpin).
/// ───────────────────────────────────────────────────────────────────────────
/// Unit: CRISPR-PAM-001 — CRISPR PAM-site finding (MolTools).
/// Checklist: docs/checklists/02_METAMORPHIC_TESTING.md, row 18.
/// Relations:
///   • MON   — appending bases cannot remove existing PAM sites ⇒ count is
///             non-decreasing (the original region's sites persist).
///   • INV   — appending a region that creates NO new PAM motif (and no
///             junction-spanning motif) leaves the PAM-site count unchanged and
///             preserves the original sites exactly.
///   • SHIFT — prepending a non-PAM-creating flank F shifts every PAM-site
///             forward-strand Position by |F| and preserves the site set.
///
/// Source (semantics, pinned from spec + source, NOT from observed output):
///   docs/algorithms/MolTools/PAM_Site_Detection.md §2.2, §4 (both strands
///   scanned; forward NGG, reverse strand = revComp NGG = forward-strand CCN);
///   src/Seqeron/Algorithms/Seqeron.Genomics.MolTools/CrisprDesigner.cs
///   (FindPamSites / FindPamSitesCore). PamSite.Position is ALWAYS a forward-strand
///   coordinate (record doc-comment, CrisprDesigner.cs); 0-based.
///
/// Why the relations hold (SpCas9, PAM = NGG, 20-nt guide, PAM 3' of target):
///   A forward-strand PAM is an "NGG" (a GG dinucleotide) at index i, reported only
///   when its 20-nt target fits UPSTREAM (i ≥ 20). A reverse-strand PAM is an "NGG"
///   on the reverse complement, i.e. a "CCN" (a CC dinucleotide) on the forward
///   strand; its 20-nt target lies on the revComp, which maps to forward positions
///   DOWNSTREAM of the CC (forward [p+3, p+22] for a CC at forward index p), so a
///   reverse site needs room toward the 3' end.
///   • A "non-PAM-creating" flank/append must introduce neither a GG nor a CC
///     dinucleotide — internally OR across any junction. We use the alternating A/T
///     region "ATAT…" (no G, no C at all), which starts and ends in A/T, so every
///     junction with it is "<base>A" / "T<base>" → never GG or CC.
///   • MON (append, any region): appending only ADDS bases downstream. No upstream
///     motif is removed and no forward target (upstream) is shortened, so every
///     original site persists ⇒ count is non-decreasing. (This holds for ANY append,
///     even a PAM-rich one — hence MON is tested on arbitrary suffixes.)
///   • INV / SHIFT require EXACT preservation, so we additionally guard the relevant
///     end: appending downstream could only *add* a reverse (CC) site that was
///     previously starved of 3' room, and prepending upstream could only *add* a
///     forward (GG) site previously starved of 5' room. Padding every test body with
///     ≥22 A/T bases on BOTH ends means every in-body motif already has its full
///     20-nt target room, so neither prepend nor append can create a new in-bounds
///     site ⇒ the site set is preserved exactly.
///   • SHIFT positions: Position is reported on the FORWARD strand for both strands,
///     so a length-|F| prepend advances every site's Position by exactly |F|.
/// ───────────────────────────────────────────────────────────────────────────
/// </summary>
[TestFixture]
[Category("Metamorphic")]
public class MolToolsMetamorphicTests
{

    /// <summary>
    /// The set of VALID primer candidates a parameter set admits over [regionStart, regionEnd],
    /// keyed by identity (Position, Sequence, IsForward) so the comparison ignores run order.
    /// </summary>
    private static HashSet<(int, string, bool)> ValidPrimerIds(
        string template, int regionStart, int regionEnd, bool forward, PrimerParameters param)
        => ThermodynamicPrimerDesigner
            .GeneratePrimerCandidates(new DnaSequence(template), regionStart, regionEnd, forward, param)
            .Where(c => c.IsValid)
            .Select(PrimerId)
            .ToHashSet();

    /// <summary>Count of VALID primer candidates over the whole template, both orientations.</summary>
    private static int ValidPrimerCount(string template, PrimerParameters param)
        => ValidPrimerIds(template, 0, template.Length, forward: true, param).Count
         + ValidPrimerIds(template, 0, template.Length, forward: false, param).Count;


    [Test]
    [Description("MON: along a chain that widens [MinTm, MaxTm] (all other params fixed) the valid-primer set grows monotonically — each narrower set is a subset of the wider one.")]
    public void GeneratePrimerCandidates_WideningTmRange_YieldsSuperset_CountNonDecreasing()
    {
        // Increasingly wide Tm windows centred on the default optimum (60 °C). Each ⊇ the prior.
        (double Min, double Max)[] tmChain =
        {
            (59.5, 60.5),
            (58, 62),
            (55, 65),
            (50, 70),
            (0, 200),   // accept any Tm — the Tm filter is effectively disabled
        };

        foreach (var template in DesignTemplates())
        {
            foreach (var forward in new[] { true, false })
            {
                HashSet<(int, string, bool)>? previousSet = null;
                int previousCount = -1;

                foreach (var (min, max) in tmChain)
                {
                    var param = ThermodynamicPrimerDesigner.DefaultParameters with { MinTm = min, MaxTm = max };
                    var currentSet = ValidPrimerIds(template, 0, template.Length, forward, param);

                    if (previousSet is not null)
                    {
                        currentSet.IsSupersetOf(previousSet).Should().BeTrue(
                            because: $"widening the Tm window to [{min},{max}] keeps every Tm that already passed " +
                                     "passing (all other filters fixed), so the valid set is a superset of the narrower one");
                        currentSet.Count.Should().BeGreaterThanOrEqualTo(previousCount,
                            because: $"a wider Tm window [{min},{max}] admits-or-keeps each candidate, never drops one — count is non-decreasing");
                    }

                    previousSet = currentSet;
                    previousCount = currentSet.Count;
                }
            }
        }
    }

    [Test]
    [Description("MON: disabling the Tm filter (accept any Tm) yields at least as many valid primers as the default narrow window, holding all else fixed.")]
    public void GeneratePrimerCandidates_DisablingTmFilter_NeverFewerThanNarrowDefault()
    {
        foreach (var template in DesignTemplates())
        {
            var narrow = ThermodynamicPrimerDesigner.DefaultParameters; // 57–63 °C
            var wide = ThermodynamicPrimerDesigner.DefaultParameters with { MinTm = 0, MaxTm = 200 };

            int narrowCount = ValidPrimerCount(template, narrow);
            int wideCount = ValidPrimerCount(template, wide);

            wideCount.Should().BeGreaterThanOrEqualTo(narrowCount,
                because: "the only changed parameter is the Tm window, widened to accept everything; " +
                         "every primer valid under the narrow window stays valid, so the count cannot drop");
        }
    }



    [Test]
    [Description("SUB: along a chain that TIGHTENS [MinGcContent, MaxGcContent] (all other params fixed) the valid-primer set shrinks monotonically — each stricter set is a subset of the looser one.")]
    public void GeneratePrimerCandidates_TighteningGcBounds_YieldsSubset_CountNonIncreasing()
    {
        // Increasingly STRICT GC windows, each ⊆ the previous one (chain from loose to strict).
        (double Min, double Max)[] gcChain =
        {
            (0, 100),     // accept any GC% — filter disabled
            (30, 70),
            (40, 60),     // the default band
            (45, 55),
            (48, 52),
        };

        foreach (var template in DesignTemplates())
        {
            foreach (var forward in new[] { true, false })
            {
                HashSet<(int, string, bool)>? previousSet = null;
                int previousCount = int.MaxValue;

                foreach (var (min, max) in gcChain)
                {
                    var param = ThermodynamicPrimerDesigner.DefaultParameters with { MinGcContent = min, MaxGcContent = max };
                    var currentSet = ValidPrimerIds(template, 0, template.Length, forward, param);

                    if (previousSet is not null)
                    {
                        previousSet.IsSupersetOf(currentSet).Should().BeTrue(
                            because: $"tightening the GC window to [{min},{max}] can only DROP candidates whose fixed GC% " +
                                     "falls outside the narrower band, so the stricter set is a subset of the looser one");
                        currentSet.Count.Should().BeLessThanOrEqualTo(previousCount,
                            because: $"a stricter GC window [{min},{max}] removes-or-keeps each candidate, never adds one — count is non-increasing");
                    }

                    previousSet = currentSet;
                    previousCount = currentSet.Count;
                }
            }
        }
    }

    [Test]
    [Description("SUB: every primer valid under a stricter GC band is also valid under a looser band that contains it (membership is preserved when loosening).")]
    public void GeneratePrimerCandidates_StrictGcBand_SubsetOfLooseBand()
    {
        foreach (var template in DesignTemplates())
        {
            foreach (var forward in new[] { true, false })
            {
                var strict = ValidPrimerIds(template, 0, template.Length, forward,
                    ThermodynamicPrimerDesigner.DefaultParameters with { MinGcContent = 45, MaxGcContent = 55 });
                var loose = ValidPrimerIds(template, 0, template.Length, forward,
                    ThermodynamicPrimerDesigner.DefaultParameters with { MinGcContent = 40, MaxGcContent = 60 });

                loose.IsSupersetOf(strict).Should().BeTrue(
                    because: "[45,55] ⊆ [40,60]; a GC% inside the tighter band is inside the wider band too, " +
                             "and no other filter changed, so the strict valid set is a subset of the loose one");
            }
        }
    }



    [Test]
    [Description("MON: extending template T → T+X and the design region to cover it preserves every valid primer anchored in the original region EXACTLY (subset), so the count is non-decreasing.")]
    public void GeneratePrimerCandidates_LongerTemplate_OriginalCandidatesPreserved_CountNonDecreasing()
    {
        var param = ThermodynamicPrimerDesigner.DefaultParameters;

        foreach (var body in DesignTemplates())
        {
            foreach (var forward in new[] { true, false })
            {
                var baseIds = ValidPrimerIds(body, 0, body.Length, forward, param);

                // Append several distinct extensions, including GC-rich and fixed-seed random ones.
                foreach (var ext in new[] { "ACGTACGTACGT", "GGGGCCCCAAAA", RandomDna(30), RandomDna(60) })
                {
                    string extended = body + ext;
                    var extendedIds = ValidPrimerIds(extended, 0, extended.Length, forward, param);

                    extendedIds.IsSupersetOf(baseIds).Should().BeTrue(
                        because: $"the first {body.Length} bases are unchanged, so every window anchored in the original " +
                                 "region yields the IDENTICAL candidate; the extended region only ADDS windows, never removes one");
                    extendedIds.Count.Should().BeGreaterThanOrEqualTo(baseIds.Count,
                        because: "extending the template/region can only add candidate windows downstream — the original valid set survives, so the count is non-decreasing");
                }
            }
        }
    }

    [Test]
    [Description("MON: appending bases after the design region's end leaves the in-region valid set IDENTICAL — the appended tail lies outside [regionStart, regionEnd].")]
    public void GeneratePrimerCandidates_AppendOutsideFixedRegion_InRegionSetUnchanged()
    {
        var param = ThermodynamicPrimerDesigner.DefaultParameters;

        foreach (var body in DesignTemplates())
        {
            foreach (var forward in new[] { true, false })
            {
                int regionEnd = body.Length; // pin the design region to the original body
                var baseIds = ValidPrimerIds(body, 0, regionEnd, forward, param);

                foreach (var ext in new[] { "ACGTACGTACGT", RandomDna(40) })
                {
                    // Hold the region FIXED at [0, regionEnd]; the appended tail is outside it.
                    var sameRegionIds = ValidPrimerIds(body + ext, 0, regionEnd, forward, param);

                    sameRegionIds.SetEquals(baseIds).Should().BeTrue(
                        because: "GeneratePrimerCandidates only enumerates windows inside [regionStart, regionEnd]; " +
                                 "bases appended beyond regionEnd are never read, so the in-region valid set is exactly preserved");
                }
            }
        }
    }



    // ─────────────────────────────────────────────────────────────────────────
    // Unit: PROBE-DESIGN-001 — hybridization probe design (MolTools).
    // Checklist: docs/checklists/02_METAMORPHIC_TESTING.md, row 24.
    //
    // API under test (ProbeDesigner.cs):
    //   DesignProbes(target, ProbeParameters?, maxProbes)
    //     1. reject if target shorter than MinLength; uppercase target.
    //     2. enumerate EVERY (start, len) window with MinLength ≤ len ≤ MaxLength.
    //     3. EvaluateProbeWithGc: raw score starts at 1.0 and SOFT penalties subtract
    //        (doc §2.2): GC out of [MinGc,MaxGc] −0.3, Tm out of [MinTm,MaxTm] −0.3,
    //        homopolymer>Max −0.2, selfComp>Max −0.2, structure −0.15, repeats −0.1,
    //        5'/3' G·C −0.02 each. A window is a CANDIDATE iff raw score > 0
    //        (EvaluateProbeWithGc returns null at score ≤ 0) — except GC beyond the
    //        early-reject band [MinGc−0.1, MaxGc+0.1] is dropped BEFORE scoring.
    //     4. return the top-`maxProbes` candidates by descending score.
    //   DesignProbes(target, ISuffixTree genomeIndex, params, maxProbes, requireUnique)
    //     walks the SAME ranked candidate stream (every candidate, lazily — B07 F59) and
    //     (requireUnique) SKIPS any probe with CheckSpecificity < 1.0, i.e. any probe
    //     occurring more than once in the genome index (INV-03: spec = 0 / 1 / 1·hits⁻¹);
    //     with requireUnique = false it instead scales Score by that specificity and
    //     re-ranks on the scaled score (the candidate SET is unchanged either way).
    //
    // A probe's run-order-independent identity is (Start, Sequence); its full record
    // (Tm, GcContent, Score, Warnings) is a PURE function of that one window's substring
    // — independent of every OTHER candidate and of the rest of the target.
    //
    // Relation DIRECTIONS are derived from this filter/scoring structure (definition),
    // NOT from observed output. To isolate each metamorphic relation we lift the
    // top-K cap by designing with a very large `maxProbes` (AllProbes) so the returned
    // set is exactly {candidates : raw score > 0} (no truncation churn), then vary ONE
    // thing and hold the rest fixed.
    //
    //   • MON (wider Tm → ≥ probes): Tm is one independent SOFT penalty. Widening
    //     [MinTm,MaxTm] removes the −0.3 Tm penalty from every window whose Tm enters
    //     the wider band and changes NOTHING else, so each window's raw score is
    //     NON-DECREASING ⇒ {raw score > 0} grows ⇒ probe count is non-decreasing and the
    //     valid set is a SUPERSET. The exact, non-vacuous mechanism is asserted directly:
    //     for a probe shared by the narrow- and wide-Tm designs, every field but Score is
    //     identical and Score rises by EXACTLY the documented 0.3 iff the probe's Tm lay
    //     outside the narrow window (the OligoArray-condition nearest-neighbour Tm — 1 M Na⁺,
    //     1 µM — of balanced 50–60mers spans ~83–96 °C, partly above the 82–90 °C Microarray
    //     window, so the toggle is real).
    //
    //   • SUB (stricter uniqueness → ⊆ results): `requireUnique=true` applies an EXTRA
    //     conjunct (skip specificity < 1) over the SAME candidate stream that
    //     `requireUnique=false` yields (which only rescales and re-ranks it). With the cap lifted it therefore yields a SUBSET:
    //     ids(requireUnique) ⊆ ids(¬requireUnique), count non-increasing. A genome that
    //     fully duplicates the target makes EVERY probe non-unique, so the strict filter
    //     empties the result (the SUB endpoint), while the lenient design stays non-empty.
    //
    //   • INV (unrelated region append → same probes): a probe's whole record depends
    //     only on its own window. Appending a region downstream leaves every window that
    //     lies entirely in the original target BYTE-IDENTICAL, so each original probe is
    //     preserved EXACTLY (same Start/Sequence/Tm/GcContent/Score/Warnings) and the
    //     extended valid set is a SUPERSET. For EXACT set equality we append after a
    //     poly-G tail of length ≥ MaxLength using a poly-G extension: then every window
    //     that reaches the appended bases starts inside the tail (|tail| ≥ MaxLength), so
    //     it is all-G ⇒ GC = 1.0 > MaxGc+0.1 ⇒ early-rejected ⇒ the append creates NO new
    //     candidate and design(T) == design(T+X) exactly.
    //
    // Source (semantics pinned from spec, NOT from observed output):
    //   docs/algorithms/MolTools/Hybridization_Probe_Design.md §2.2 (penalty table),
    //   §2.4 (INV-01 raw-score-positive shortlist, INV-03 specificity 0/1/1·hits⁻¹),
    //   §3.1/§3.3 (requireUnique filters specificity < 1.0; top-maxProbes after ranking),
    //   §4.2 (Microarray default: length 50–60, Tm 82–90 °C at 1 M Na⁺ / 1 µM, GC 0.40–0.60).
    // ─────────────────────────────────────────────────────────────────────────

    /// <summary>Microarray default probe parameters (length 50–60, Tm 82–90 at 1 M Na⁺ / 1 µM, GC 0.40–0.60).</summary>
    private static ThermodynamicProbeDesigner.ProbeParameters Microarray => ThermodynamicProbeDesigner.Defaults.Microarray;

    /// <summary>All probes a parameter set admits over the target (cap lifted ⇒ the full score-positive set).</summary>
    private static List<ProbeDesigner.Probe> DesignAll(string target, ThermodynamicProbeDesigner.ProbeParameters param)
        => ThermodynamicProbeDesigner.DesignProbes(target, param, AllProbes).ToList();

    private static Dictionary<(int, string), ProbeDesigner.Probe> DesignById(string target, ThermodynamicProbeDesigner.ProbeParameters param)
        => DesignAll(target, param).ToDictionary(ProbeId);



    [Test]
    [Description("SUB: requireUnique=true applies an extra 'specificity = 1' filter over the same candidate stream, so its probe set is a subset of the requireUnique=false design (count non-increasing).")]
    public void DesignProbes_RequireUnique_YieldsSubsetOfLenientDesign()
    {
        foreach (var target in ProbeTargets())
        {
            string upper = target.ToUpperInvariant();
            // Genome where the FIRST part of the target is duplicated: probes drawn from that
            // region occur twice (non-unique), the rest occur once (unique). 'N' separators
            // carry no probe base, so they break any cross-junction match.
            int dupLen = Math.Min(upper.Length, Microarray.MinLength + 15);
            var genome = global::SuffixTree.SuffixTree.Build(upper + "NNNNNNNN" + upper[..dupLen]);

            var unique = ThermodynamicProbeDesigner.DesignProbes(target, genome, Microarray, AllProbes, requireUnique: true)
                .Select(ProbeId).ToHashSet();
            var lenient = ThermodynamicProbeDesigner.DesignProbes(target, genome, Microarray, AllProbes, requireUnique: false)
                .Select(ProbeId).ToHashSet();

            lenient.IsSupersetOf(unique).Should().BeTrue(
                because: "requiring genome uniqueness only SKIPS candidates with specificity < 1 from the same candidate stream, " +
                         "so the unique-only probe set is a subset of the lenient one");
            unique.Count.Should().BeLessThanOrEqualTo(lenient.Count,
                because: "a stricter uniqueness requirement removes-or-keeps each candidate, never adds one — count is non-increasing");
        }
    }

    [Test]
    [Description("SUB endpoint: against a genome that fully duplicates the target every probe is non-unique, so requireUnique empties the result while the lenient design stays non-empty.")]
    public void DesignProbes_FullyDuplicatedGenome_RequireUnique_EmptiesWhileLenientKeepsProbes()
    {
        foreach (var target in ProbeTargets())
        {
            string upper = target.ToUpperInvariant();
            // Two full copies of the target ⇒ every probe (a substring of the target) occurs ≥ 2×.
            var genome = global::SuffixTree.SuffixTree.Build(upper + "NNNNNNNN" + upper);

            int uniqueCount = ThermodynamicProbeDesigner.DesignProbes(target, genome, Microarray, AllProbes, requireUnique: true).Count();
            int lenientCount = ThermodynamicProbeDesigner.DesignProbes(target, genome, Microarray, AllProbes, requireUnique: false).Count();

            lenientCount.Should().BeGreaterThan(0,
                because: "the target yields probes, and the lenient design keeps non-unique candidates");
            uniqueCount.Should().Be(0,
                because: "every probe occurs in both target copies (specificity = 0.5 < 1), so the strict uniqueness filter removes them all");
        }
    }

    private static ProbeDesigner.ProbeValidation Validate(
        string probe, IEnumerable<string> references, int maxMismatches, double selfThr = SelfCompThreshold)
        => ThermodynamicProbeDesigner.ValidateProbe(probe, references, maxMismatches, selfThr);

    [Test]
    [Description("INV: a self-dimer is the hetero-dimer of S with a copy of itself, and the dimer is symmetric in strand order — dimer(A,B) equals dimer(B,A).")]
    public void Dimer_SelfDimer_EqualsHeteroDimerOfSWithItself()
    {
        foreach (string s in SelfDimerOligos)
        {
            string copy = new string(s.ToCharArray()); // distinct reference, same content

            double selfTm = ThermodynamicPrimerDesigner.CalculateSelfDimerMeltingTemperature(s);
            double heteroTm = ThermodynamicPrimerDesigner.CalculateDimerMeltingTemperature(s, copy);

            // Non-vacuity: S genuinely forms a self-dimer (a real duplex, not NaN).
            selfTm.Should().NotBe(double.NaN, because: $"'{s}' is partially self-complementary and forms a self-dimer");
            heteroTm.Should().BeApproximately(selfTm, 1e-9,
                because: $"the self-dimer of '{s}' is exactly the hetero-dimer of '{s}' with a copy of itself");

            var self = PrimerDesigner.FindMostStableDimer(s, s);
            var hetero = PrimerDesigner.FindMostStableDimer(s, copy);
            self.Should().NotBeNull();
            hetero.Should().Be(self, because: $"FindMostStableDimer('{s}','{s}') equals the hetero-dimer with a copy");
        }
    }



    [Test]
    [Description("MON: along a chain that widens [MinTm,MaxTm] (all else fixed) the valid-probe set grows monotonically — each narrower set is a subset of the wider one, count non-decreasing.")]
    public void DesignProbes_WideningTmWindow_YieldsSuperset_CountNonDecreasing()
    {
        // Increasingly wide Tm windows, each ⊇ the prior (Microarray default 82–90 outward).
        (double Min, double Max)[] tmChain =
        {
            (82, 90),   // Microarray default
            (75, 95),
            (50, 100),
            (0, 200),   // Tm filter effectively disabled
        };

        foreach (var target in ProbeTargets())
        {
            HashSet<(int, string)>? previousSet = null;
            int previousCount = -1;

            foreach (var (min, max) in tmChain)
            {
                var ids = DesignAll(target, Microarray with { MinTm = min, MaxTm = max })
                    .Select(ProbeId).ToHashSet();

                if (previousSet is not null)
                {
                    ids.IsSupersetOf(previousSet).Should().BeTrue(
                        because: $"widening the Tm window to [{min},{max}] only removes the −{ProbeTmPenalty} Tm penalty " +
                                 "(all other penalties fixed), so every window that already scored > 0 still does — the valid set is a superset");
                    ids.Count.Should().BeGreaterThanOrEqualTo(previousCount,
                        because: $"a wider Tm window [{min},{max}] can only raise raw scores, never lower one past 0 — probe count is non-decreasing");
                }

                previousSet = ids;
                previousCount = ids.Count;
            }
        }
    }

    [Test]
    [Description("MON mechanism: a probe shared by the narrow- and wide-Tm designs keeps every field but Score; its Score rises by EXACTLY the documented 0.3 iff its Tm was outside the narrow window.")]
    public void DesignProbes_WideningTmWindow_RaisesSharedProbeScoreByExactlyTheTmPenalty()
    {
        bool sawTmToggle = false;

        foreach (var target in ProbeTargets())
        {
            // narrow = Microarray default (82–90 °C); wide = Tm filter disabled.
            var narrow = DesignById(target, Microarray with { MinTm = 82, MaxTm = 90 });
            var wide = DesignById(target, Microarray with { MinTm = 0, MaxTm = 200 });

            foreach (var (id, probe) in narrow)
            {
                // A probe valid under the narrow Tm window stays valid when the window widens
                // (its raw score only rises), so the wide design must contain its twin.
                wide.Should().ContainKey(id,
                    because: "removing the Tm penalty cannot drop a probe that already scored > 0 under the narrow window");
                var twin = wide[id];

                twin.Sequence.Should().Be(probe.Sequence);
                twin.Tm.Should().Be(probe.Tm, because: "Tm is a pure function of the probe window, independent of the Tm parameter range");
                twin.GcContent.Should().Be(probe.GcContent, because: "GC content depends only on the window, not on the Tm range");

                // The ONLY score difference is the Tm penalty: present in narrow iff the probe's
                // Tm is outside [82,90]; never present in the wide (0–200) window.
                bool tmOutsideNarrow = probe.Tm < 82 || probe.Tm > 90;
                double expectedDelta = tmOutsideNarrow ? ProbeTmPenalty : 0.0;
                if (tmOutsideNarrow) sawTmToggle = true;

                twin.Score.Should().BeApproximately(probe.Score + expectedDelta, 1e-9,
                    because: tmOutsideNarrow
                        ? "the probe's Tm is outside the narrow window, so widening it removes exactly the documented 0.3 Tm penalty"
                        : "the probe's Tm is already inside the narrow window, so widening the window leaves its score unchanged");
            }
        }

        sawTmToggle.Should().BeTrue(
            because: "the 1 M / 1 µM nearest-neighbour Tm of balanced 50–60mers reaches above 90 °C, so at least one probe's Tm toggles — the relation is exercised, not vacuous");
    }



    [Test]
    [Description("INV: appending any unrelated downstream region preserves every original probe EXACTLY (same Start/Sequence/Tm/GC/Score/Warnings); the extended valid set is a superset.")]
    public void DesignProbes_AppendUnrelatedRegion_PreservesOriginalProbesExactly()
    {
        // Microarray preset (1 M Na⁺ / 1 µM, Tm 82–90, G+C 0.40–0.60, ntthal screen) with the length window narrowed
        // to 50–51 nt (B07 heavy tier, HEAVY-2): with the cap lifted every window gets the three ntthal recursions
        // (≈ 30 ms per 50–60-mer in Debug), so the full 50–60 window over the 25 designs below took 9–11 min (B06
        // "hang"). The relation — a window's record depends only on the window — is independent of the length range;
        // two lengths keep the multi-length enumeration / tie order in play. Same targets, same appends, same checks.
        var param = Microarray with { MaxLength = Microarray.MinLength + 1 };
        foreach (var target in ProbeTargets())
        {
            var baseProbes = DesignAll(target, param);
            baseProbes.Should().NotBeEmpty(because: "each target is constructed to yield Microarray probes");

            foreach (var ext in new[] { RandomDna(30), RandomDna(60), "GGGGCCCCAAAATTTT", NonPamRegion(40) })
            {
                var extById = DesignById(target + ext, param);

                foreach (var bp in baseProbes)
                {
                    extById.Should().ContainKey(ProbeId(bp),
                        because: "windows lying entirely in the original target are byte-identical after a downstream append, " +
                                 "so every original probe still appears");
                    var twin = extById[ProbeId(bp)];

                    twin.Sequence.Should().Be(bp.Sequence);
                    twin.Tm.Should().Be(bp.Tm, because: "Tm depends only on the probe window, which the downstream append does not touch");
                    twin.GcContent.Should().Be(bp.GcContent);
                    twin.Score.Should().Be(bp.Score, because: "the heuristic score reads only the probe's own window, so a distant append cannot move it");
                    twin.Warnings.Should().Equal(bp.Warnings, because: "warnings are recorded from the probe window alone");
                }
            }
        }
    }

    [Test]
    [Description("INV exact: with a poly-G tail ≥ MaxLength, a poly-G append creates no new candidate (every append-touching window is all-G ⇒ GC-rejected), so the designed probe set is exactly unchanged.")]
    public void DesignProbes_PolyGTailThenPolyGAppend_DesignSetExactlyUnchanged()
    {
        // good region (~50% GC, yields probes) + poly-G tail of length ≥ MaxLength.
        const string goodRegion = "ACGTGACTGACTGGATCAGTCAGTACGATCGATGCATGCATCGTAGCATGCATGCATGCAACGTGACTGACTGGATCAGT";
        string tail = new string('G', Microarray.MaxLength + 4); // ≥ MaxLength ⇒ append windows start inside the tail
        string baseSeq = goodRegion + tail;

        var baseById = DesignById(baseSeq, Microarray);
        baseById.Should().NotBeEmpty(because: "the balanced good region yields Microarray probes; the poly-G tail itself is GC-rejected");

        foreach (int appendLen in new[] { 1, 10, 30 })
        {
            string extended = baseSeq + new string('G', appendLen);
            var extById = DesignById(extended, Microarray);

            extById.Keys.Should().BeEquivalentTo(baseById.Keys,
                because: $"a poly-G×{appendLen} append after a poly-G tail of length ≥ MaxLength adds only all-G windows " +
                         "(GC = 100% > MaxGc+0.1), which are early-rejected — so no probe is added or removed");

            foreach (var (id, bp) in baseById)
            {
                var twin = extById[id];
                twin.Sequence.Should().Be(bp.Sequence);
                twin.Tm.Should().Be(bp.Tm);
                twin.GcContent.Should().Be(bp.GcContent);
                twin.Score.Should().Be(bp.Score, because: "the surviving probes are byte-identical windows, so every metric is preserved exactly");
            }
        }
    }


    [Test]
    [Description("MON: along an increasing maxMismatches chain the off-target hit count is non-decreasing and IsValid is non-increasing — lowering the threshold can only ADD passing probes, never remove one.")]
    public void ValidateProbe_RaisingMismatchTolerance_HitsNonDecreasing_ValidityNonIncreasing()
    {
        int[] toleranceChain = { 0, 1, 2, 3, 4, 5, 6 };

        foreach (var probe in CleanProbes())
        {
            var references = MismatchLadder(probe);

            // Precondition: the probe is "clean" (R = true), so IsValid ⇔ offTargetHits ≤ 1.
            var clean = Validate(probe, references, maxMismatches: 0);
            clean.SelfComplementarity.Should().BeLessThanOrEqualTo(SelfCompThreshold,
                because: "the construction requires a non-self-complementary probe so validity tracks only the off-target count");
            clean.HasSecondaryStructure.Should().BeFalse(
                because: "the construction requires a probe with no hairpin so validity tracks only the off-target count");

            int previousHits = -1;
            bool? previousValid = null;
            bool sawValid = false, sawInvalid = false;

            foreach (int k in toleranceChain)
            {
                var result = Validate(probe, references, maxMismatches: k);

                result.OffTargetHits.Should().BeGreaterThanOrEqualTo(previousHits,
                    because: $"raising the tolerance to {k} can only ADD approximate matches (matches(k) ⊆ matches(k+1)), so the hit count is non-decreasing");

                // IsValid ⇔ offTargetHits ≤ 1 for a clean probe — the algebra, checked exactly.
                result.IsValid.Should().Be(result.OffTargetHits <= 1,
                    because: "for a clean probe (R = true) the validity rule reduces to offTargetHits ≤ 1");

                if (previousValid is not null)
                    (result.IsValid && !previousValid.Value).Should().BeFalse(
                        because: $"a higher mismatch tolerance ({k}) never turns an invalid probe valid — validity is non-increasing in the threshold");

                sawValid |= result.IsValid;
                sawInvalid |= !result.IsValid;
                previousHits = result.OffTargetHits;
                previousValid = result.IsValid;
            }

            sawValid.Should().BeTrue(because: "at low tolerance only the exact on-target is found (1 hit) ⇒ the probe passes");
            sawInvalid.Should().BeTrue(because: "at high tolerance the near-duplicates are also found (>1 hit) ⇒ the probe fails — the relation is exercised, not vacuous");
        }
    }

    [Test]
    [Description("MON (set): the set of probes that PASS validation is downward-closed in maxMismatches — lowering the threshold yields a superset of passing probes (count non-decreasing).")]
    public void ValidateProbe_LoweringThreshold_PassingSetIsSuperset()
    {
        // Each probe carries its own mismatch ladder; pass/fail is decided per probe.
        var cases = CleanProbes().Select(p => (Probe: p, Refs: MismatchLadder(p))).ToList();
        int[] descendingThresholds = { 6, 4, 3, 2, 1, 0 };

        HashSet<string>? higherSet = null;
        foreach (int k in descendingThresholds)
        {
            var passing = cases.Where(c => Validate(c.Probe, c.Refs, k).IsValid).Select(c => c.Probe).ToHashSet();

            if (higherSet is not null)
                passing.IsSupersetOf(higherSet).Should().BeTrue(
                    because: $"lowering the threshold to {k} keeps every probe that already passed and may add more — the passing set is a superset");

            higherSet = passing;
        }
    }

    [Test]
    [Description("MON dual: the specificity score is non-increasing as the off-target hit count grows (1 hit → 1.0, N hits → 1/N).")]
    public void ValidateProbe_MoreOffTargetHits_SpecificityScoreNonIncreasing()
    {
        foreach (var probe in CleanProbes())
        {
            var references = MismatchLadder(probe);
            double previousScore = double.PositiveInfinity;

            foreach (int k in new[] { 0, 1, 2, 3, 4, 5, 6 })
            {
                var result = Validate(probe, references, k);

                result.SpecificityScore.Should().BeLessThanOrEqualTo(previousScore,
                    because: $"as the tolerance rises to {k} the hit count grows and specificity (1/hits for hits>1) is non-increasing");
                if (result.OffTargetHits >= 1)
                    result.SpecificityScore.Should().BeApproximately(1.0 / result.OffTargetHits, 1e-12,
                        because: "the spec maps hits>0 to 1/hits (and 1 hit to 1.0), independent of which references produced them");

                previousScore = result.SpecificityScore;
            }
        }
    }



    [Test]
    [Description("INV: ValidateProbe is deterministic — repeated calls on identical input return an identical validation record.")]
    public void ValidateProbe_SameInput_IsDeterministic()
    {
        foreach (var probe in CleanProbes())
        {
            var references = MismatchLadder(probe);
            foreach (int k in new[] { 0, 2, 4 })
            {
                var first = Validate(probe, references, k);
                for (int i = 0; i < 4; i++)
                    AssertSameValidation(Validate(probe, references, k), first,
                        "ValidateProbe is a pure function of (probe, references, params), so repeated calls must agree exactly");
            }
        }
    }

    [Test]
    [Description("INV: permuting the reference order leaves the validation record unchanged — offTargetHits is a sum over references, which is order-independent.")]
    public void ValidateProbe_ReferenceOrderPermutation_ResultUnchanged()
    {
        foreach (var probe in CleanProbes())
        {
            var references = MismatchLadder(probe);
            var reversed = Enumerable.Reverse(references).ToList();
            var rotated = references.Skip(1).Concat(references.Take(1)).ToList();

            foreach (int k in new[] { 0, 2, 4, 6 })
            {
                var baseline = Validate(probe, references, k);
                AssertSameValidation(Validate(probe, reversed, k), baseline,
                    "the total hit count sums over references regardless of order, and every other field reads the probe alone");
                AssertSameValidation(Validate(probe, rotated, k), baseline,
                    "rotating the reference list cannot change the summed hit count or the probe-only metrics");
            }
        }
    }

    [Test]
    [Description("INV: ValidateProbe is case-insensitive — lower/mixed casing the probe and references gives the same record, since both sides are uppercased before analysis.")]
    public void ValidateProbe_CaseFolding_ResultUnchanged()
    {
        foreach (var probe in CleanProbes())
        {
            var references = MismatchLadder(probe);
            foreach (int k in new[] { 0, 2, 4 })
            {
                var baseline = Validate(probe, references, k);

                AssertSameValidation(
                    Validate(probe.ToLowerInvariant(), references.Select(r => r.ToLowerInvariant()).ToList(), k),
                    baseline,
                    "the probe and each reference are uppercased before matching, so casing carries no information");

                var mixedRefs = references
                    .Select(r => new string(r.Select((ch, i) => i % 2 == 0 ? char.ToLowerInvariant(ch) : ch).ToArray()))
                    .ToList();
                AssertSameValidation(Validate(probe, mixedRefs, k), baseline,
                    "case folding is total, so any case pattern of the same letters yields the same validation");
            }
        }
    }


    /// <summary>
    /// A flank/append region guaranteed to create NO PAM motif on either strand:
    /// it contains only A and T, so it has neither a GG (forward NGG) nor a CC
    /// (reverse-strand CCN) dinucleotide, and—starting and ending in A/T—it forms
    /// none across any junction either.
    /// </summary>
    private static string NonPamRegion(int length)
    {
        var chars = new char[length];
        for (int i = 0; i < length; i++)
            chars[i] = (i % 2 == 0) ? 'A' : 'T';
        return new string(chars);
    }

    /// <summary>Generates a random DNA string of the given length over {A,C,G,T}.</summary>
    private static string RandomDna(int length)
    {
        const string bases = "ACGT";
        var chars = new char[length];
        for (int i = 0; i < length; i++)
            chars[i] = bases[Rng.Next(bases.Length)];
        return new string(chars);
    } // == ThermoConstants.WallaceMaxLength

    private static double Tm(string primer) => PrimerDesigner.CalculateMeltingTemperature(primer);




    // ─────────────────────────────────────────────────────────────────────────
    // Unit: PRIMER-DESIGN-001 — PCR primer design (MolTools).
    // Checklist: docs/checklists/02_METAMORPHIC_TESTING.md, row 22.
    //
    // API under test (PrimerDesigner.cs):
    //   GeneratePrimerCandidates(DnaSequence template, int regionStart, int regionEnd,
    //                            bool forward, PrimerParameters? parameters)
    //     → enumerates EVERY (start, len) window in [regionStart, regionEnd] with
    //       MinLength ≤ len ≤ MaxLength, evaluates each via EvaluatePrimer, and yields
    //       the resulting PrimerCandidate (valid AND invalid alike — it does NOT pre-filter).
    //   EvaluatePrimer(string seq, int position, bool isForward, PrimerParameters?)
    //     → sets IsValid = true iff EVERY screen passes; the relevant filters here are:
    //          MinLength ≤ Length ≤ MaxLength
    //          MinGcContent ≤ GcContent ≤ MaxGcContent
    //          MinTm ≤ MeltingTemperature ≤ MaxTm
    //          homopolymer ≤ MaxHomopolymer, dinuc-repeat ≤ MaxDinucleotideRepeats,
    //          !HasHairpin, (Check3PrimeStability ⇒ ΔG ≥ −9), (Avoid3PrimeGC ⇒ GC clamp).
    //
    // The MR phrase "primers" / "results" denotes the set of VALID candidates a parameter
    // set admits over a fixed region — i.e. GeneratePrimerCandidates(...).Where(IsValid).
    // A candidate's run-order-independent identity is (Position, Sequence, IsForward); its
    // metrics (GcContent, Tm, Length, IsValid) are a PURE function of that local window,
    // independent of every other candidate and of the rest of the template.
    //
    // Relation DIRECTIONS are derived from the FILTER STRUCTURE (definition), NOT from
    // observed output. Every filter is an independent conjunct, so widening exactly ONE
    // window while holding ALL other parameters fixed can only flip a candidate
    // invalid→valid (never valid→invalid):
    //
    //   • MON (wider Tm range → ≥ primers): EvaluatePrimer accepts a candidate's Tm iff
    //     MinTm ≤ Tm ≤ MaxTm. Widening the window to [MinTm', MaxTm'] ⊇ [MinTm, MaxTm]
    //     (MinTm' ≤ MinTm, MaxTm' ≥ MaxTm) keeps every previously-passing Tm passing and
    //     may admit more. With all other params fixed, the valid set is a SUPERSET and the
    //     count is non-decreasing: Valid(narrow Tm) ⊆ Valid(wide Tm).
    //
    //   • SUB (stricter GC% → ⊆ results): acceptance requires MinGc ≤ GC ≤ MaxGc.
    //     TIGHTENING the window to [MinGc', MaxGc'] ⊆ [MinGc, MaxGc] can only DROP
    //     candidates whose fixed GC% falls outside the narrower band; none are added.
    //     Along a tightening chain the valid set shrinks monotonically:
    //     Valid(strict GC) ⊆ Valid(loose GC); count non-increasing.
    //
    //   • MON (longer template → ≥ candidates): GeneratePrimerCandidates enumerates every
    //     window fitting inside [regionStart, regionEnd]. Extend the template T → T + X and
    //     the design region [0, |T|−1] → [0, |T+X|−1]. Every window that fit in the original
    //     region still fits in the extended one and—because the first |T| bases are
    //     byte-for-byte unchanged—evaluates to the IDENTICAL candidate (same Position,
    //     Sequence, IsForward, IsValid). So the original valid set is preserved EXACTLY as a
    //     subset and the count is non-decreasing: Valid(T) ⊆ Valid(T+X).
    //
    // For each relation only ONE parameter (or the region) is varied; ALL others are held
    // fixed — the essence of the metamorphic relation. Templates include hand-built
    // GC-balanced sequences that actually yield primers plus fixed-seed random ones.
    // ─────────────────────────────────────────────────────────────────────────

    /// <summary>
    /// Design templates rich enough to yield several valid primers under DefaultParameters:
    /// 40–60% GC, length ≥ 60, no long homopolymer/dinucleotide runs, plus fixed-seed random.
    /// </summary>
    private static IEnumerable<string> DesignTemplates()
    {
        // ~50% GC, varied composition, 60 nt.
        yield return "ACGTGACTGACTGGATCAGTCAGTACGATCGATGCATGCATCGTAGCATGCATGCATGCA";
        // Another balanced 72-nt template.
        yield return "TGCATGCAGTCAGTACGTACGATCGATCGTAGCTAGCATGCATGCATCGATCGATCAGTCAGTACGTACGTA";
        // GC-leaning but still mostly in-range, 64 nt.
        yield return "GCGATCGATGCATCGATCGTAGCTAGCATCGATCGATGCATGCATCGATCGATGCATCGATCGT";
        // Fixed-seed random templates (relations hold for arbitrary input too).
        yield return RandomDna(80);
        yield return RandomDna(120);
        yield return RandomDna(160);
    }

    /// <summary>Run-order-independent identity of a primer candidate.</summary>
    private static (int Position, string Sequence, bool Fwd) PrimerId(PrimerCandidate c)
        => (c.Position, c.Sequence, c.IsForward);

    /// <summary>Documented Tm-window penalty (doc §2.2): a probe whose Tm is out of [MinTm,MaxTm] loses 0.3.</summary>
    private const double ProbeTmPenalty = 0.3;

    /// <summary>
    /// Large probe cap that lifts the top-K truncation for these fixtures, so a design
    /// returns exactly {candidate windows : raw score &gt; 0} — making the set relations
    /// (superset / subset / equality) reflect the SCORING/FILTER definition, not the cut.
    /// </summary>
    private const int AllProbes = 100_000;

    /// <summary>Run-order-independent identity of a designed probe.</summary>
    private static (int Start, string Sequence) ProbeId(ProbeDesigner.Probe p) => (p.Start, p.Sequence);

    /// <summary>
    /// Targets ≥ MinLength with ~50% GC and varied composition, so several Microarray
    /// probes (50–60 nt) fit, plus fixed-seed random targets — relations must hold for
    /// arbitrary input too.
    /// </summary>
    private static IEnumerable<string> ProbeTargets()
    {
        yield return "ACGTGACTGACTGGATCAGTCAGTACGATCGATGCATGCATCGTAGCATGCATGCATGCAACGTGACTGACTGGATCAGT";
        yield return "TGCATGCAGTCAGTACGTACGATCGATCGTAGCTAGCATGCATGCATCGATCGATCAGTCAGTACGTACGTAGCATCGAT";
        yield return "GCGATCGATGCATCGATCGTAGCTAGCATCGATCGATGCATGCATCGATCGATGCATCGATCGTACGATCGTAGCTAGCA";
        yield return RandomDna(90);
        yield return RandomDna(120);
    }




    // ─────────────────────────────────────────────────────────────────────────
    // Unit: PROBE-VALID-001 — hybridization probe validation (MolTools).
    // Checklist: docs/checklists/02_METAMORPHIC_TESTING.md, row 25.
    //
    // API under test (ThermodynamicProbeDesigner.ValidateProbe):
    //   ValidateProbe(probe, references, maxMismatches = 3, selfComplementarityThreshold = 0.3)
    //     → ProbeValidation(IsValid, SpecificityScore, OffTargetHits, SelfComplementarity,
    //                       HasSecondaryStructure, Issues).
    //   offTargetHits = Σ_refs |{ i : Hamming(ref[i..i+|probe|], probe) ≤ maxMismatches }|
    //                   (ungapped fixed-length substitution scan, both sides uppercased).
    //   IsValid = issues.Count == 0  ∨  (offTargetHits ≤ 1 ∧ selfComp ≤ 0.4), where issues
    //   holds an off-target item iff offTargetHits > 1, a self-comp item iff
    //   selfComp > selfComplementarityThreshold, and a structure item iff hasStructure.
    //
    //   Algebra (definition, NOT observed output): substituting the issue conditions,
    //     IsValid = (offTargetHits ≤ 1) ∧ R,   R = ((selfComp ≤ thr ∧ ¬hasStructure) ∨ selfComp ≤ 0.4)
    //   and R is INDEPENDENT of maxMismatches (selfComp / hasStructure read the probe alone).
    //
    // Relation DIRECTIONS (derived from that algebra and the spec, NOT from output):
    //   Source: docs/algorithms/MolTools/Probe_Validation.md §2.2 (specificity map),
    //           §2.4 (INV-01..04), §4.1 step 5 (the IsValid rule), §5.2.
    //
    //   • MON (lower specificity threshold → more pass): `maxMismatches` is the off-target
    //     detection stringency — the mismatch tolerance below which a window counts as a
    //     cross-hybridization site. FindApproximateMatches is monotone in it:
    //     matches(k) ⊆ matches(k+1) (a window within k mismatches is within k+1), so
    //     offTargetHits(k) is NON-DECREASING in k. Since IsValid = (offTargetHits ≤ 1) ∧ R
    //     with R fixed, IsValid(k) is NON-INCREASING: lowering the threshold can only ADD
    //     passing probes (the passing set is downward-closed in k). The dual SpecificityScore
    //     (1 / hits for hits > 1) is likewise non-increasing as hits grow.
    //
    //   • INV (same input → same result): ValidateProbe is a pure function whose result
    //     depends only on the MULTISET of references (offTargetHits is a SUM over them) and
    //     on case-folded sequences (probe and each reference are uppercased). Hence the full
    //     record is invariant under (a) repeated calls, (b) permuting the reference order,
    //     and (c) changing the case of the probe and/or references.
    //
    // Construction: a "clean" probe has selfComp ≤ threshold AND no secondary structure, so
    // R = true and IsValid reduces to (offTargetHits ≤ 1) — isolating the maxMismatches MR.
    // References are exact-length copies of the probe mutated at a fixed number of positions
    // (via the same A↔T / C↔G flip used above), so each reference contributes one window at
    // a KNOWN Hamming distance and the hit count crosses 1 as the tolerance rises. The clean
    // probes' R-precondition is asserted from the returned fields, so the setup is self-checking.
    // ─────────────────────────────────────────────────────────────────────────

    /// <summary>Default self-complementarity issue threshold used by ValidateProbe.</summary>
    private const double SelfCompThreshold = 0.3;

    /// <summary>
    /// 20-nt probes that are NOT self-complementary and form no hairpin, so the validity
    /// rule reduces to (offTargetHits ≤ 1) — letting the maxMismatches relation be tested
    /// in isolation. (Each probe's R-precondition is verified in-test from the result.)
    /// </summary>
    private static IEnumerable<string> CleanProbes()
    {
        yield return "AAAAAAAAAAGGGGGGGGGG";
        yield return "AAAAAAAAAACCCCCCCCCC";
        yield return "TTTTTTTTTTCCCCCCCCCC";
    }

    /// <summary>
    /// References for a clean probe: an EXACT copy (0 mismatches) plus copies at 2, 4 and 6
    /// mismatches. Each is the probe's length, so it contributes exactly one window whose
    /// Hamming distance to the probe is the stated value; the hit count therefore steps up as
    /// the mismatch tolerance rises, crossing the validity boundary (offTargetHits ≤ 1).
    /// </summary>
    private static List<string> MismatchLadder(string probe) => new()
    {
        probe,                                  // 0 mismatches (exact on-target)
        MutateAt(probe, 0, 1),                  // 2 mismatches
        MutateAt(probe, 0, 1, 2, 3),            // 4 mismatches
        MutateAt(probe, 0, 1, 2, 3, 4, 5),      // 6 mismatches
    };

    private static void AssertSameValidation(ProbeDesigner.ProbeValidation a, ProbeDesigner.ProbeValidation b, string why)
    {
        a.IsValid.Should().Be(b.IsValid, because: why);
        a.SpecificityScore.Should().Be(b.SpecificityScore, because: why);
        a.OffTargetHits.Should().Be(b.OffTargetHits, because: why);
        a.SelfComplementarity.Should().Be(b.SelfComplementarity, because: why);
        a.HasSecondaryStructure.Should().Be(b.HasSecondaryStructure, because: why);
        a.Issues.Should().Equal(b.Issues, because: why);
    }


    /// <summary>Deterministic RNG — seed fixed so random inputs are reproducible.</summary>
    private static readonly Random Rng = new(20260619);


    // ═══════════════════════════════════════════════════════════════════
    //  PRIMER-DIMER-001 — most-stable primer dimer ΔG (MolTools)
    // ═══════════════════════════════════════════════════════════════════
    //
    // Theory (SantaLucia & Hicks 2004 unified NN; Primer3/ntthal;
    //   docs/algorithms/MolTools/DNA_Dimer_Tm.md):
    //   FindMostStableDimer aligns two oligos antiparallel over every gapless offset and scores the
    //   most stable contiguous Watson-Crick run; a SELF-dimer is the same method with both strands
    //   equal. Two metamorphic relations (checklist row 242):
    //
    //   • INV (self-dimer of S equals hetero-dimer(S,S)): the self-dimer is, by definition, the
    //     hetero-dimer of S with a second copy of S — and the dimer of two strands is symmetric in
    //     strand order (the same physical duplex), so dimer(A,B) and dimer(B,A) agree.
    //   • MON (extending the WC alignment lowers ΔG): each additional Watson-Crick base pair in the
    //     aligned run adds one more stabilising nearest-neighbour stack, so a longer complementary
    //     alignment gives a strictly more negative dimer ΔG°37.
    //
    // API under test: PrimerDesigner.FindMostStableDimer / .CalculateSelfDimerMeltingTemperature /
    //   .CalculateDimerMeltingTemperature (DimerResult).


    // Deterministic, partially self-complementary oligos that form a real self-dimer.
    private static readonly string[] SelfDimerOligos =
    {
        "CGATCGATCG",
        "GCGCGCGC",
        "ACGTACGTACGT",
        "GGGGCCCC",
    };

    /// <summary>
    /// Returns a copy of <paramref name="target"/> with a mismatch (vs itself) introduced
    /// at each position in <paramref name="positions"/> via <see cref="FlipBase"/>.
    /// </summary>
    private static string MutateAt(string target, params int[] positions)
    {
        var chars = target.ToCharArray();
        foreach (int p in positions)
            chars[p] = FlipBase(chars[p]);
        return new string(chars);
    }


    /// <summary>
    /// Flip a base to a guaranteed mismatch WITHOUT creating a new G/C run: A↔T and
    /// C↔G keep the base within its weak/strong pair, so neither a GG nor a CC
    /// dinucleotide can be introduced by the flip.
    /// </summary>
    private static char FlipBase(char b) => b switch
    {
        'A' => 'T',
        'T' => 'A',
        'C' => 'G',
        'G' => 'C',
        _ => throw new ArgumentException($"Unexpected base '{b}'.", nameof(b)),
    };
}
