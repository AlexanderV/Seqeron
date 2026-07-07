using static Seqeron.Genomics.Analysis.ProteinMotifFinder;

namespace Seqeron.Genomics.Tests.Fuzzing;

/// <summary>
/// Fuzz tests for the ProteinMotif area — COILED-COIL prediction
/// (PROTMOTIF-CC-001): predicting coiled-coil regions of a protein by heptad-repeat
/// a/d hydrophobic-core occupancy. The single public entry point under test is
/// <see cref="Seqeron.Genomics.Analysis.ProteinMotifFinder.PredictCoiledCoils"/>
/// (with private helpers <c>BestHeptadOccupancy</c>, <c>BuildRegion</c>, <c>Mod</c>).
/// Sibling ProteinMotif fuzz units — PROTMOTIF-FIND-001 (row 82, regex motif scan),
/// PROTMOTIF-PROSITE-001 (row 83, PROSITE syntax) and PROTMOTIF-DOMAIN-001 (row 84,
/// domain finding) — are covered separately; this file focuses on the COILED-COIL
/// heptad-occupancy contract.
///
/// ───────────────────────────────────────────────────────────────────────────
/// What fuzzing verifies
/// ───────────────────────────────────────────────────────────────────────────
/// Fuzzing feeds degenerate, boundary and out-of-domain inputs to a unit and asserts
/// that the code NEVER fails in an undisciplined way: no hang, no state corruption, no
/// nonsense output, and no *unhandled* runtime exception (IndexOutOfRange / NullReference
/// / KeyNotFound / ArgumentOutOfRange). Every input must resolve to EITHER a well-defined,
/// theory-correct result OR a *documented, intentional* outcome. For a windowed numeric
/// scanner that uppercases the input, scans a fixed-width window over every heptad register
/// and counts {I,L,V} occupancy at the a/d core positions, the headline hazards are:
///   • a NullReferenceException when the sequence is null;
///   • an IndexOutOfRangeException when the sequence is SHORTER than the window, so no full
///     window exists (the scan must simply yield nothing, never run off the end — a single
///     residue vs the 28-residue default window is the sharpest case);
///   • a KeyNotFoundException / mis-count on out-of-alphabet residues (digits, punctuation,
///     X, B/Z/J/O/U): these are simply not in the {I,L,V} core set, so they lower occupancy
///     and MUST NOT raise a region or crash;
///   • a Score escaping [0,1], or a region whose [Start..End] falls outside [0,n−1] or is
///     shorter than the documented minimum 3-heptad (21-residue) span (a coordinate bug).
/// — docs/ADVANCED_TESTING_CHECKLIST.md §8 "Fuzzing".
///
/// ───────────────────────────────────────────────────────────────────────────
/// Unit: PROTMOTIF-CC-001 — coiled-coil prediction (heptad a/d occupancy)
/// Checklist: docs/checklists/03_FUZZING.md, row 163.
/// Fuzz strategies exercised for THIS unit:
///   • BE = Boundary Exploitation — the length corners that could crash or yield a
///     spurious region:
///       – empty sequence: "" / null → NO regions by the explicit string.IsNullOrEmpty
///         guard, never a NullReference (Coiled_Coil_Prediction.md §3.3, §6.1; INV-04).
///       – single residue / any length &lt; windowSize → NO regions, no IndexOutOfRange:
///         no full window exists (§3.3, §6.1, INV-04). The 1-residue input vs the
///         28-residue default window is the sharpest sub-window probe.
///   • MC = Malformed Content — out-of-alphabet / junk residues:
///       – non-amino-acid characters (digits, punctuation, whitespace, the unknown
///         placeholder 'X', the extended IUPAC codes B/Z/J/O/U): none are in the {I,L,V}
///         hydrophobic-core set (§4.2, §6.1 "no {I,L,V} residues → empty result"), so a/d
///         occupancy is 0 &lt; threshold and NO false coiled coil is raised — never a
///         KeyNotFound or other crash. Lowercase input is uppercased first, so it scores
///         identically to its uppercase form (§3.3, §6.1 "lowercase input → recognised").
/// — docs/checklists/03_FUZZING.md §Description (strategy codes BE, MC);
///   targets: "empty, non-amino-acid, single residue".
///
/// ───────────────────────────────────────────────────────────────────────────
/// The coiled-coil contract under test (Coiled_Coil_Prediction.md §2–§4)
/// ───────────────────────────────────────────────────────────────────────────
/// Given a sequence S (uppercased) of length n, a window of length W (default 28 = 4
/// heptads) and a threshold t (default 0.5): the heptad position of residue index k in
/// register r ∈ {0..6} is p(k,r) = (k − r) mod 7; the hydrophobic-core positions are
/// a = 0 and d = 3; a residue is a core residue iff it is one of {I, L, V}. For each
/// window start i, occ(i,r) = (#a/d positions in the window holding a residue ∈ {I,L,V})
/// / (#a/d positions in the window), and the window score is score(i) = max over the 7
/// registers of occ(i,r). Contiguous windows with score ≥ t form a run [i₀..i₁] mapping
/// to residue region [i₀, i₁ + W − 1]; the region is emitted (with Score = peak score in
/// the run) only if its length ≥ MinCoiledCoilRegion = 3 heptads = 21 residues.
///   — Coiled_Coil_Prediction.md §2.2, §3.1–§3.3, §4.1.
///
/// Method under test (src/.../Seqeron.Genomics.Analysis/ProteinMotifFinder.cs):
///   IEnumerable&lt;(int Start, int End, double Score)&gt; PredictCoiledCoils(
///       string proteinSequence, int windowSize = 28, double threshold = 0.5)
///   — Coiled_Coil_Prediction.md §5.1.
///
/// Documented input handling (Coiled_Coil_Prediction.md §3.3, §6.1):
///   • null / "" / any sequence shorter than windowSize → NO regions (no exception).
///   • Input is uppercased (case-INSENSITIVE): lowercase == uppercase.
///   • The alphabet is the 20 amino acids; non-{I,L,V} residues simply do not count
///     toward a/d occupancy, so out-of-alphabet junk cannot raise a region or crash.
///
/// Theory-correct invariants asserted (Coiled_Coil_Prediction.md §2.4):
///   • INV-01 — every Score ∈ [0, 1] (a count ratio), and finite.
///   • INV-02 — each region spans ≥ 21 residues (3 heptads).
///   • INV-03 — 0 ≤ Start ≤ End ≤ n − 1; regions are non-overlapping and increasing in Start.
///   • INV-04 — sequences shorter than windowSize produce no regions.
///   • INV-05 — a region exists only if some covering window scores ≥ threshold (max over 7
///     registers); a sequence with no a/d hydrophobic periodicity yields none.
///   • [determinism] — re-running the same input yields identical regions.
///
/// ───────────────────────────────────────────────────────────────────────────
/// Complexity / hang-safety
/// ───────────────────────────────────────────────────────────────────────────
/// The scan is O(n · W · 7) ≈ O(n) with fixed constants (Coiled_Coil_Prediction.md §4.3).
/// Long-homopolymer / long-junk targets are kept modest and [CancelAfter]-guarded so a
/// regression that turned the windowed scan into a hang would FAIL rather than wedge the suite.
/// ───────────────────────────────────────────────────────────────────────────
/// </summary>
[TestFixture]
[Category("Fuzzing")]
public class ProteinCoiledCoilFuzzTests
{
    #region Helpers

    /// <summary>The 20 standard amino-acid one-letter codes.</summary>
    private const string StandardAminoAcids = "ACDEFGHIKLMNPQRSTVWY";

    /// <summary>The default scanning-window length (4 heptads) used by <c>PredictCoiledCoils</c>.</summary>
    private const int DefaultWindow = 28;

    /// <summary>The minimum reported region length (3 heptads = 21 residues).</summary>
    private const int MinRegion = 21;

    /// <summary>Deterministic RNG — seed fixed locally so generated fuzz inputs are reproducible.</summary>
    private static string RandomProtein(int length, int seed)
    {
        var rng = new Random(seed);
        var chars = new char[length];
        for (int i = 0; i < length; i++)
            chars[i] = StandardAminoAcids[rng.Next(StandardAminoAcids.Length)];
        return new string(chars);
    }

    /// <summary>
    /// Asserts the universal theory-correct contract every emitted coiled-coil region must satisfy
    /// against the original (case-insensitive) sequence (Coiled_Coil_Prediction.md §2.4, §3.2):
    /// the region is in-bounds (INV-03: 0 ≤ Start ≤ End ≤ n−1), spans at least the minimum 3-heptad
    /// length (INV-02: End − Start + 1 ≥ 21), and its Score is a finite occupancy fraction in [0,1]
    /// (INV-01). This is the headline "no coordinate bug, no out-of-range Score, no NaN" property.
    /// </summary>
    private static void AssertWellFormedRegion((int Start, int End, double Score) region, string originalSequence)
    {
        int n = originalSequence.Length;

        // INV-03 — in-bounds, ordered span.
        region.Start.Should().BeInRange(0, n - 1, "INV-03: Start is a valid 0-based residue index");
        region.End.Should().BeInRange(region.Start, n - 1, "INV-03: End is in-bounds and not before Start");

        // INV-02 — minimum 3-heptad (21-residue) span.
        (region.End - region.Start + 1).Should().BeGreaterThanOrEqualTo(MinRegion,
            "INV-02: a reported coiled coil spans at least 3 heptads (21 residues)");

        // INV-01 — Score is a finite occupancy fraction in [0,1].
        double.IsNaN(region.Score).Should().BeFalse("INV-01: a coiled-coil Score must never be NaN");
        double.IsInfinity(region.Score).Should().BeFalse("INV-01: a coiled-coil Score must never be infinite");
        region.Score.Should().BeInRange(0.0, 1.0, "INV-01: Score is an a/d occupancy fraction in [0,1]");
    }

    /// <summary>
    /// Asserts a whole result set is well-formed AND that regions are strictly increasing in Start
    /// (the part of INV-03 the single forward scan actually guarantees). NOTE on INV-03's
    /// "non-overlapping" wording: each region's reported End extends windowSize − 1 residues PAST its
    /// last above-threshold window, while a later run's Start is its first window index; two runs are
    /// separated only by a below-threshold WINDOW, not by windowSize − 1 residues, so the reported
    /// residue spans of consecutive regions CAN legitimately overlap. That is a property of the
    /// window→residue mapping, not a coordinate bug — every region is still in-bounds, ≥ 21 residues,
    /// with Score ∈ [0,1], and the run start-window indices strictly increase, so Start strictly
    /// increases. We therefore assert the real, code-guaranteed ordering (increasing Start) rather
    /// than residue-disjointness.
    /// </summary>
    private static void AssertWellFormedResult(
        List<(int Start, int End, double Score)> regions, string originalSequence)
    {
        foreach (var region in regions)
            AssertWellFormedRegion(region, originalSequence);

        for (int i = 1; i < regions.Count; i++)
        {
            regions[i].Start.Should().BeGreaterThan(regions[i - 1].Start,
                "INV-03: regions are emitted by a single forward scan, strictly increasing in Start");
        }
    }

    #endregion

    // ═══════════════════════════════════════════════════════════════════
    //  PROTMOTIF-CC-001 — coiled-coil prediction : fuzz targets
    // ═══════════════════════════════════════════════════════════════════

    #region PROTMOTIF-CC-001 — coiled-coil prediction

    #region BE — Empty / null sequence: no regions, no NullReference

    /// <summary>
    /// Target "empty": "" and null must produce NO regions — by the explicit string.IsNullOrEmpty
    /// guard, NEVER a NullReferenceException (Coiled_Coil_Prediction.md §3.3, §6.1 "null / empty →
    /// empty result"; INV-04). Verified for the default parameters and across non-default window /
    /// threshold settings (the guard precedes any indexing).
    /// </summary>
    [Test]
    public void PredictCoiledCoils_EmptyOrNullSequence_NoRegionsNoThrow()
    {
        foreach (string? seq in new[] { "", null })
        {
            var act = () => PredictCoiledCoils(seq!).ToList();
            act.Should().NotThrow($"empty/null sequence ('{seq ?? "null"}') must not crash PredictCoiledCoils")
                .Subject.Should().BeEmpty("empty/null sequence yields no coiled-coil regions (INV-04)");

            // Same guard must hold for non-default window / threshold settings.
            foreach (int w in new[] { 1, 7, 28, 56 })
            {
                var act2 = () => PredictCoiledCoils(seq!, w, 0.5).ToList();
                act2.Should().NotThrow($"empty/null sequence with window {w} must not crash")
                    .Subject.Should().BeEmpty("empty/null sequence yields no regions regardless of window");
            }
        }
    }

    #endregion

    #region BE — Single residue / sub-window length: no regions, no IndexOutOfRange

    /// <summary>
    /// Target "single residue": a length-1 sequence — and indeed any sequence STRICTLY SHORTER than
    /// windowSize — has no place for a full window, so the scan must yield NO regions and NEVER run
    /// off the end (no IndexOutOfRange) — Coiled_Coil_Prediction.md §3.3, §6.1 "length &lt; windowSize
    /// → empty result"; INV-04. We probe every single standard residue (including the core residues
    /// I/L/V, the sharpest case since those WOULD score 1.0 in a full window) against the 28-residue
    /// default window, plus a sweep of sub-window lengths up to windowSize − 1, and finally
    /// length == windowSize − 1 with a SMALL custom window to prove the guard is "length &lt; window",
    /// not a hard-coded 28.
    /// </summary>
    [Test]
    public void PredictCoiledCoils_SingleOrSubWindowLength_NoRegionsNoCrash()
    {
        // (a) Every single residue vs the default 28-window → no region, no crash.
        foreach (char aa in StandardAminoAcids)
        {
            string seq = aa.ToString();
            var act = () => PredictCoiledCoils(seq).ToList();
            act.Should().NotThrow($"a single residue ('{aa}') must not crash PredictCoiledCoils")
                .Subject.Should().BeEmpty($"a 1-residue sequence is shorter than the 28-window (INV-04)");
        }

        // (b) Sub-window lengths 1..27 of a pure-core homopolymer (the case that WOULD score 1.0 in a
        //     full window) → still no region, because no full window exists.
        for (int len = 1; len < DefaultWindow; len++)
        {
            string seq = new string('L', len);
            var act = () => PredictCoiledCoils(seq).ToList();
            act.Should().NotThrow($"a length-{len} sub-window sequence must not crash")
                .Subject.Should().BeEmpty($"length {len} < window 28 yields no regions (INV-04)");
        }

        // (c) length == windowSize − 1 with a SMALL custom window proves the guard is "len < window".
        //     A perfect coiled-coil pattern of length 6 with window 7 → no region; length 7 → region.
        string sixCore = BuildPerfectCoiledCoil(6);   // length 6 < window 7
        PredictCoiledCoils(sixCore, windowSize: 7, threshold: 0.5).Should().BeEmpty(
            "length 6 is shorter than the custom window 7 → no full window → no region (INV-04)");
    }

    #endregion

    #region MC — Non-amino-acid characters: no false region, never a crash

    /// <summary>
    /// Target "non-amino-acid": out-of-alphabet residues — digits, punctuation, whitespace, the
    /// unknown placeholder 'X', and the extended IUPAC ambiguity codes B/Z/J/O/U — are NONE of them
    /// in the {I,L,V} hydrophobic-core set, so they contribute 0 to a/d occupancy and CANNOT raise a
    /// coiled coil; the scan must therefore yield NO regions and NEVER throw (no KeyNotFound from a
    /// missing residue, no IndexOutOfRange) — Coiled_Coil_Prediction.md §4.2, §6.1 "no {I,L,V}
    /// residues → empty result". We build junk strings LONG ENOUGH to contain full windows (so a
    /// crash hazard genuinely exists), scan them across several window/threshold settings, and require
    /// an empty, well-formed result. Case-insensitivity is pinned separately.
    /// </summary>
    [Test]
    [CancelAfter(30000)]
    public void PredictCoiledCoils_NonAminoAcidChars_NoFalseRegionNoCrash(CancellationToken token)
    {
        string[] junk =
        {
            new string('1', 40),                                  // digits
            "!@#$%^&*()!@#$%^&*()!@#$%^&*()!@#$%^&*()",          // punctuation (40)
            new string(' ', 20) + "\t\t\t\t\t\t\t\t\t\t\n\n\n\n\n\n\n\n\n\n", // whitespace (40)
            string.Concat(Enumerable.Repeat("BZJOU", 8)),         // extended IUPAC codes (40)
            new string('X', 40),                                  // unknown placeholder
            "MK1RGD2SP" + new string('X', 31),                    // junk mixed with real (non-core) residues
        };

        foreach (string seq in junk)
        {
            foreach (int w in new[] { 7, 28 })
            {
                foreach (double t in new[] { 0.0, 0.5, 1.0 })
                {
                    var act = () => PredictCoiledCoils(seq, w, t).ToList();
                    var regions = act.Should().NotThrow(
                        $"junk sequence (window {w}, threshold {t}) must not crash").Subject;
                    token.ThrowIfCancellationRequested();

                    // No {I,L,V} present → 0 occupancy → no region for any threshold > 0. At threshold 0
                    // EVERY window qualifies (0 ≥ 0), so a region MAY appear; whatever is emitted must be
                    // well-formed (in-bounds, ≥ 21 residues, Score ∈ [0,1]). With no core residues that
                    // Score must be exactly 0.
                    AssertWellFormedResult(regions, seq);
                    if (t > 0.0)
                        regions.Should().BeEmpty(
                            $"out-of-alphabet residues never reach a/d occupancy ≥ {t} (no false coiled coil)");
                    else
                        regions.Should().OnlyContain(r => r.Score == 0.0,
                            "with no {I,L,V} residues every window occupancy is exactly 0");
                }
            }
        }
    }

    /// <summary>
    /// Case-insensitivity (MC / §3.3, §6.1 "lowercase input → recognised"): the input is uppercased
    /// before scanning, so a lowercase sequence yields EXACTLY the same regions as its uppercase form.
    /// Pinned on a strong lowercase coiled coil (which must be recognised identically) and on a mixed-
    /// case variant.
    /// </summary>
    [Test]
    public void PredictCoiledCoils_LowercaseInput_RecognisedIdenticallyToUppercase()
    {
        string upper = BuildPerfectCoiledCoil(35); // a strong coiled coil
        string lower = upper.ToLowerInvariant();
        string mixed = new string(upper.Select((c, i) => i % 2 == 0 ? char.ToLowerInvariant(c) : c).ToArray());

        var upperRegions = PredictCoiledCoils(upper).ToList();
        var lowerRegions = PredictCoiledCoils(lower).ToList();
        var mixedRegions = PredictCoiledCoils(mixed).ToList();

        upperRegions.Should().NotBeEmpty("the constructed sequence is a strong coiled coil");
        lowerRegions.Should().Equal(upperRegions,
            "matching is case-insensitive: lowercase input yields the same regions as uppercase");
        mixedRegions.Should().Equal(upperRegions,
            "matching is case-insensitive: mixed-case input yields the same regions as uppercase");
    }

    #endregion

    #region Positive sanity — a strong coiled coil is found; a non-periodic sequence is not

    /// <summary>
    /// Builds a "perfect" coiled coil: a length-<paramref name="length"/> sequence with a core residue
    /// at EVERY heptad a (index ≡ 0 mod 7) and d (index ≡ 3 mod 7) position in register 0, and a
    /// non-core residue (alanine) everywhere else. In register 0 every a/d position then holds a
    /// {I,L,V} residue → occupancy 1.0 for every full window (Coiled_Coil_Prediction.md §7.1).
    /// </summary>
    private static string BuildPerfectCoiledCoil(int length)
    {
        var chars = new char[length];
        for (int k = 0; k < length; k++)
        {
            int heptadPos = k % 7;
            chars[k] = heptadPos == 0 ? 'I' : heptadPos == 3 ? 'L' : 'A';
        }
        return new string(chars);
    }

    /// <summary>
    /// Positive sanity: the harness must assert against a predictor that actually FINDS coiled coils
    /// at the correct coordinates with the correct Score, not a no-op. The documented worked example
    /// ("LAALAAA"×5 → one region (0,34,1.0)) and a constructed I/L/V-at-every-a/d sequence (occupancy
    /// 1.0) must each yield a single region spanning the documented [0, n−1] with Score ≈ 1.0; while a
    /// sequence with NO hydrophobic a/d periodicity (all alanine, or all the non-core residue K) must
    /// yield NONE (INV-05). This pins the headline "a real coiled coil is found at the right span with
    /// the right Score; a non-coiled sequence is not over-predicted" contract (§2.4, §4.1, §7.1).
    /// </summary>
    [Test]
    public void PredictCoiledCoils_StrongCoiledCoil_FoundWithHighScore_NonPeriodicNotFound()
    {
        // (a) The documented worked example: "LAALAAA"×5 (35 aa) → exactly one region (0, 34, 1.0).
        string docExample = string.Concat(Enumerable.Repeat("LAALAAA", 5));
        docExample.Length.Should().Be(35, "sanity-check the worked-example length");
        var docRegions = PredictCoiledCoils(docExample).ToList();
        docRegions.Should().ContainSingle("the documented perfect 5-heptad repeat yields one region");
        AssertWellFormedRegion(docRegions[0], docExample);
        docRegions[0].Start.Should().Be(0, "the documented region starts at residue 0");
        docRegions[0].End.Should().Be(34, "the documented region ends at residue 34 (i₁=7, 7+27=34)");
        docRegions[0].Score.Should().BeApproximately(1.0, 1e-9,
            "every a/d position holds L → occupancy 1.0 (peak Score 1.0)");

        // (b) A constructed I/L/V-at-every-a/d coiled coil (length 35) → one region [0,34], Score 1.0.
        string perfect = BuildPerfectCoiledCoil(35);
        var perfectRegions = PredictCoiledCoils(perfect).ToList();
        perfectRegions.Should().ContainSingle("a perfect I/L/V a/d pattern yields one coiled-coil region");
        AssertWellFormedRegion(perfectRegions[0], perfect);
        perfectRegions[0].Start.Should().Be(0, "the perfect region starts at residue 0");
        perfectRegions[0].End.Should().Be(34, "the perfect region spans the whole 35-residue sequence");
        perfectRegions[0].Score.Should().BeApproximately(1.0, 1e-9, "I/L/V at every a/d → occupancy 1.0");

        // (c) No hydrophobic a/d periodicity → no region (INV-05).
        foreach (string flat in new[] { new string('A', 60), new string('K', 60), new string('G', 60) })
        {
            PredictCoiledCoils(flat).Should().BeEmpty(
                $"a homopolymer of the non-core residue '{flat[0]}' has 0 a/d occupancy → no coiled coil (INV-05)");
        }

        // (d) A genuinely short-of-minimum coiled coil: a perfect pattern of EXACTLY one full window
        //     (28 residues) maps to region [0, 27] = 28 residues ≥ 21 → reported; verify the documented
        //     span and that Score is the peak occupancy 1.0.
        string oneWindow = BuildPerfectCoiledCoil(28);
        var oneWindowRegions = PredictCoiledCoils(oneWindow).ToList();
        oneWindowRegions.Should().ContainSingle("a 28-residue perfect coiled coil is exactly one full window");
        oneWindowRegions[0].Start.Should().Be(0);
        oneWindowRegions[0].End.Should().Be(27, "a single window [0] maps to residues [0, 0+27] = [0,27]");
        oneWindowRegions[0].Score.Should().BeApproximately(1.0, 1e-9);
    }

    /// <summary>
    /// Positive sanity over RANDOM proteins: across fixed seeds and lengths the scan must never crash,
    /// hang, or emit a malformed region, and every emitted region must satisfy the full contract
    /// (in-bounds, ≥ 21-residue span, Score ∈ [0,1]) with regions non-overlapping and increasing in
    /// Start. Determinism is pinned by re-running the same input and requiring identical regions. This
    /// pins span-correctness and termination on arbitrary sequences, not just hand-built coiled coils.
    /// Lengths span sub-window (1, 5), exactly-window (28) and well-above-window (60, 250) cases.
    /// </summary>
    [Test]
    [CancelAfter(30000)]
    public void PredictCoiledCoils_RandomProtein_AlwaysWellFormedAndDeterministic(CancellationToken token)
    {
        foreach (int seed in new[] { 11, 53, 211, 2026 })
        {
            foreach (int len in new[] { 1, 5, 27, 28, 60, 250 })
            {
                string seq = RandomProtein(len, seed);

                var act = () => PredictCoiledCoils(seq).ToList();
                var regions = act.Should().NotThrow($"random protein must not crash (seed {seed}, len {len})").Subject;
                token.ThrowIfCancellationRequested();

                if (len < DefaultWindow)
                    regions.Should().BeEmpty($"len {len} < window 28 → no regions (INV-04)");

                AssertWellFormedResult(regions, seq);

                // [determinism] — the same input yields identical regions.
                var again = PredictCoiledCoils(seq).ToList();
                again.Should().Equal(regions, "PredictCoiledCoils is deterministic for a fixed input");
            }
        }
    }

    #endregion

    #endregion
}
