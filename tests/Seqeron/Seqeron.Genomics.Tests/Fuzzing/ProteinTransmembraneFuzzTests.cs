using static Seqeron.Genomics.Analysis.ProteinMotifFinder;

namespace Seqeron.Genomics.Tests.Fuzzing;

/// <summary>
/// Fuzz tests for the ProteinMotif area — TRANSMEMBRANE HELIX prediction
/// (PROTMOTIF-TM-001): predicting membrane-spanning α-helices in a protein
/// sequence with the Kyte &amp; Doolittle (1982) hydropathy method. The single
/// public entry point under test is
/// <see cref="ProteinMotifFinder.PredictTransmembraneHelices(string, int, double)"/>,
/// which computes the arithmetic-mean hydropathy over a sliding window
/// (private <c>CalculateHydropathyProfile</c> + <c>HydropathyScale</c>) and reports
/// each maximal run of above-threshold windows as a candidate segment. Sibling
/// ProteinMotif units — motif finding (PROTMOTIF-FIND-001, row 82), PROSITE
/// (row 83), domains (row 84), signal peptide / coiled-coil / low-complexity /
/// disorder (rows 163–167) — are covered separately; this file focuses on the
/// TRANSMEMBRANE hydropathy contract.
///
/// ───────────────────────────────────────────────────────────────────────────
/// What fuzzing verifies
/// ───────────────────────────────────────────────────────────────────────────
/// Fuzzing feeds degenerate, boundary and out-of-domain inputs to a unit and
/// asserts that the code NEVER fails in an undisciplined way: no hang, no state
/// corruption, no nonsense output, and no *unhandled* runtime exception
/// (IndexOutOfRange / NullReference / ArgumentOutOfRange / NaN). Every input must
/// resolve to EITHER a well-defined, theory-correct result OR a *documented,
/// intentional* outcome. For a sliding-window hydropathy filter that uppercases
/// the input, builds an arithmetic-mean profile and walks threshold-crossing runs,
/// the headline hazards are:
///   • a NullReferenceException when the sequence is null;
///   • an IndexOutOfRangeException at the SHORT / window-edge boundary, where a
///     sequence shorter than the window admits no window at all (must yield
///     nothing, never run off the end);
///   • a reported region whose End exceeds the last residue index (a coordinate
///     bug — the End must be clamped to length−1);
///   • a NaN Score on a window composed entirely of NON-STANDARD residues (the
///     mean divides by the count of scale-bearing residues; a 0/0 must NOT leak);
///   • a FALSE positive region on an all-hydrophilic sequence whose every window
///     mean is far below the threshold.
/// — docs/ADVANCED_TESTING_CHECKLIST.md §8 "Fuzzing".
///
/// ───────────────────────────────────────────────────────────────────────────
/// Unit: PROTMOTIF-TM-001 — transmembrane helix prediction (Kyte-Doolittle)
/// Checklist: docs/checklists/03_FUZZING.md, row 168.
/// Fuzz strategy exercised for THIS unit:
///   • BE = Boundary Exploitation — the composition / length corners that could
///     crash, hang, or produce a false or out-of-bounds region:
///       – all-hydrophilic (Arg/Lys/Asp homopolymers, h ≈ −4.5 / −3.9 / −3.5,
///         far below 1.6): every window mean is well below threshold → NO
///         transmembrane region, no false positive (Transmembrane_Helix_Prediction.md
///         §2.4 INV-03, §6.1 "all-hydrophilic → no segments").
///       – all-hydrophobic (Ile/Leu/Val homopolymers, h = 4.5 / 3.8 / 4.2, far
///         above 1.6 over a long stretch): ONE transmembrane region spanning the
///         whole run with Score = the residue's hydropathy value (the headline
///         positive). For a uniform run of length n ≥ w the region is (0, n−1, h)
///         (INV-02, INV-03; §7.1 worked example).
///       – short: a protein SHORTER than windowSize (19) admits no window → NO
///         regions per the guard (windowSize ≤ 0 OR length < windowSize → empty),
///         no IndexOutOfRange at the window edge (§3.3, §6.1, INV-04).
/// — docs/checklists/03_FUZZING.md §Description (strategy code BE);
///   targets: "all-hydrophilic, all-hydrophobic, short".
///
/// ───────────────────────────────────────────────────────────────────────────
/// The transmembrane-prediction contract under test
/// ───────────────────────────────────────────────────────────────────────────
/// Given a protein sequence S, window width w (default 19) and threshold T
/// (default 1.6), the Kyte-Doolittle profile P(i) is the arithmetic mean of the
/// per-residue hydropathy values over the window [i, i+w−1]:
///     P(i) = (1/w) · Σ_{j=i}^{i+w−1} h(S[j])           (§2.2, §2.3)
/// Non-standard residues (X/B/Z/* and any non-AA char) carry no scale value and
/// are EXCLUDED from the mean (the mean divides by the count of scale-bearing
/// residues; a window of only non-standard residues yields 0, never NaN — §3.3,
/// §5.2, §6.1). A maximal run of windows with P(i) ≥ T is reported as ONE segment
/// iff its residue span ≥ the minimum helix length (= w = 19); for the canonical
/// window any run with at least one passing window qualifies (§4.2). The reported
/// region is:
///   • Start = first above-threshold window's first residue (0-based);
///   • End   = last covered residue = lastPassingProfileIndex + w − 1, CLAMPED to
///     length−1 (0-based inclusive);
///   • Score = peak (maximum) window mean within the run (≥ T by INV-01).
/// Null / empty / windowSize ≤ 0 / length &lt; windowSize → empty, no exception.
///   — docs/algorithms/ProteinMotif/Transmembrane_Helix_Prediction.md §2.2–§2.4,
///     §3.1–§3.3, §4.1–§4.2, §5.4, §6.1.
///
/// Method under test (src/.../Seqeron.Genomics.Analysis/ProteinMotifFinder.cs):
///   IEnumerable&lt;(int Start, int End, double Score)&gt;
///       PredictTransmembraneHelices(string proteinSequence,
///           int windowSize = 19, double threshold = 1.6)
///   — Transmembrane_Helix_Prediction.md §5.1.
///
/// Theory-correct invariants asserted (Transmembrane_Helix_Prediction.md §2.4):
///   • INV-01 — every reported segment's peak Score ≥ T.
///   • INV-02 — 0 ≤ Start ≤ End ≤ length−1 (End clamped; no run-off-the-end).
///   • INV-03 — a uniform run of one residue (length ≥ w) gives P(i) = h(r)
///     everywhere, so the all-hydrophobic homopolymer yields (0, n−1, h(r)).
///   • INV-04 — null / empty / shorter-than-window / windowSize ≤ 0 → no segments.
///   • [finite-score] — Score is always finite (never NaN / ±∞), even on a window
///     of all non-standard residues.
///
/// ───────────────────────────────────────────────────────────────────────────
/// Complexity / hang-safety
/// ───────────────────────────────────────────────────────────────────────────
/// Profile + scan is O(n·w) with an O(n) linear scan over the profile
/// (§4.3); there is no recursion or backtracking. The long-homopolymer and
/// random-sequence targets are kept modest and [CancelAfter]-guarded so a
/// regression that turned the scan into a hang would FAIL rather than wedge the suite.
/// ───────────────────────────────────────────────────────────────────────────
/// </summary>
[TestFixture]
[Category("Fuzzing")]
public class ProteinTransmembraneFuzzTests
{
    #region Helpers

    /// <summary>Default Kyte-Doolittle window width (Kyte &amp; Doolittle 1982).</summary>
    private const int DefaultWindow = 19;

    /// <summary>Default mean-hydropathy threshold (Kyte &amp; Doolittle 1982).</summary>
    private const double DefaultThreshold = 1.6;

    /// <summary>The 20 standard amino-acid one-letter codes.</summary>
    private const string StandardAminoAcids = "ACDEFGHIKLMNPQRSTVWY";

    /// <summary>
    /// Kyte-Doolittle hydropathy values (Transmembrane_Helix_Prediction.md §2.2), mirrored here so
    /// the tests assert against the DOCUMENTED scale rather than reaching into the private source map.
    /// </summary>
    private static readonly IReadOnlyDictionary<char, double> Hydropathy = new Dictionary<char, double>
    {
        ['I'] = 4.5, ['V'] = 4.2, ['L'] = 3.8, ['F'] = 2.8, ['C'] = 2.5,
        ['M'] = 1.9, ['A'] = 1.8, ['G'] = -0.4, ['T'] = -0.7, ['S'] = -0.8,
        ['W'] = -0.9, ['Y'] = -1.3, ['P'] = -1.6, ['H'] = -3.2, ['E'] = -3.5,
        ['Q'] = -3.5, ['D'] = -3.5, ['N'] = -3.5, ['K'] = -3.9, ['R'] = -4.5,
    };

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
    /// Asserts the universal theory-correct contract every emitted segment must satisfy against the
    /// (case-insensitive) sequence of length <paramref name="n"/> and threshold <paramref name="threshold"/>
    /// (Transmembrane_Helix_Prediction.md §2.4, §3.2): in-bounds inclusive coordinates with the End clamped
    /// (INV-02: 0 ≤ Start ≤ End ≤ n−1), a residue span at least the minimum helix length (= the window
    /// width), and a finite peak Score that is ≥ the threshold (INV-01). This is the headline "no
    /// coordinate bug, no run-off-the-end, no NaN, no below-threshold region" property.
    /// </summary>
    private static void AssertWellFormedRegion(
        (int Start, int End, double Score) region, int n, int windowSize, double threshold)
    {
        // INV-02 — in-bounds, ordered, clamped coordinates.
        region.Start.Should().BeInRange(0, n - 1, "a segment Start is a valid 0-based residue index");
        region.End.Should().BeInRange(region.Start, n - 1,
            "a segment End is in-bounds (clamped to length−1) and not before its Start");

        // Minimum helix span: a reported region spans at least the window width (= MinTransmembraneHelixLength).
        (region.End - region.Start + 1).Should().BeGreaterThanOrEqualTo(windowSize,
            "a reported transmembrane segment spans at least the minimum helix length (the window width)");

        // [finite-score] + INV-01 — finite peak ≥ threshold.
        double.IsNaN(region.Score).Should().BeFalse("a segment Score must never be NaN");
        double.IsInfinity(region.Score).Should().BeFalse("a segment Score must never be infinite");
        region.Score.Should().BeGreaterThanOrEqualTo(threshold,
            "INV-01: a segment's peak window mean is at or above the detection threshold");
    }

    #endregion

    // ═══════════════════════════════════════════════════════════════════
    //  PROTMOTIF-TM-001 — transmembrane helix prediction : fuzz targets
    // ═══════════════════════════════════════════════════════════════════

    #region PROTMOTIF-TM-001 — transmembrane helix prediction

    #region BE — all-hydrophilic: no false-positive region

    /// <summary>
    /// Target "all-hydrophilic": a homopolymer of a strongly hydrophilic residue (Arg h = −4.5,
    /// Lys −3.9, Asp/Asn/Glu/Gln −3.5) has every window mean equal to that residue's (very negative)
    /// hydropathy (INV-03), which is far below the 1.6 threshold — so the prediction must yield NO
    /// transmembrane region, never a false positive, and never an IndexOutOfRange on the (long)
    /// profile (Transmembrane_Helix_Prediction.md §2.4 INV-03, §6.1 "all-hydrophilic → no segments").
    /// We probe every hydrophilic residue (h &lt; 0) at lengths well above the window so the guard is
    /// not what suppresses the region — it is genuinely the sub-threshold mean.
    /// </summary>
    [Test]
    public void Tm_AllHydrophilic_NoRegion()
    {
        // Every residue whose Kyte-Doolittle value is below the threshold, as a long homopolymer.
        foreach (char aa in StandardAminoAcids.Where(c => Hydropathy[c] < DefaultThreshold))
        {
            foreach (int len in new[] { DefaultWindow, 30, 100 })
            {
                string seq = new string(aa, len);
                var regions = PredictTransmembraneHelices(seq).ToList();

                // A residue below threshold can never produce an above-threshold window mean (INV-03).
                regions.Should().BeEmpty(
                    $"a length-{len} homopolymer of '{aa}' (h = {Hydropathy[aa]}) is below the 1.6 threshold "
                    + "→ no transmembrane region (no false positive)");
            }
        }

        // A canonical strongly-hydrophilic mixture (all residues h < 0) likewise yields nothing.
        var charged = PredictTransmembraneHelices(new string('R', 15) + new string('D', 15) + new string('K', 15))
            .ToList();
        charged.Should().BeEmpty("a fully charged/hydrophilic protein contains no membrane-spanning helix");
    }

    #endregion

    #region BE — all-hydrophobic: exactly one region spanning the whole run

    /// <summary>
    /// Target "all-hydrophobic" (the headline POSITIVE): a homopolymer of a strongly hydrophobic
    /// residue (Ile h = 4.5, Leu 3.8, Val 4.2, Phe 2.8, Cys 2.5, Met 1.9, Ala 1.8 — all &gt; 1.6) of
    /// length n ≥ w has every window mean equal to that residue's hydropathy (INV-03), all above the
    /// threshold, so the single maximal run maps to EXACTLY ONE region spanning residues (0, n−1)
    /// with Score = h(residue) (INV-01, INV-02, §7.1). We assert the documented span AND score for
    /// every above-threshold residue, at several lengths. A residue whose hydropathy is at/below the
    /// threshold (the rest of the alphabet) must NOT produce a region — the complementary check.
    /// </summary>
    [Test]
    public void Tm_AllHydrophobic_SingleRegionSpanningRun()
    {
        foreach (char aa in StandardAminoAcids.Where(c => Hydropathy[c] > DefaultThreshold))
        {
            foreach (int len in new[] { DefaultWindow, 25, 50 })
            {
                string seq = new string(aa, len);
                var regions = PredictTransmembraneHelices(seq).ToList();

                regions.Should().ContainSingle(
                    $"a length-{len} hydrophobic homopolymer of '{aa}' (h = {Hydropathy[aa]} > 1.6) yields "
                    + "exactly one transmembrane region");

                var region = regions[0];
                AssertWellFormedRegion(region, len, DefaultWindow, DefaultThreshold);

                // INV-03 / §7.1: a uniform run gives P(i) = h(r) everywhere; the single run spans the
                // whole sequence (Start 0 .. End n−1, clamped) with peak Score = the residue hydropathy.
                region.Start.Should().Be(0, "the uniform hydrophobic run opens at residue 0");
                region.End.Should().Be(len - 1, "the uniform hydrophobic run is clamped to the last residue index");
                region.Score.Should().BeApproximately(Hydropathy[aa], 1e-9,
                    $"INV-03: every window mean of a '{aa}' homopolymer equals h('{aa}') = {Hydropathy[aa]}");
            }
        }

        // Complement: a residue at/below the threshold (e.g. Gly, Ser, Pro, charged) yields no region
        // even as a long homopolymer — confirms the test is not trivially "any long homopolymer passes".
        foreach (char aa in StandardAminoAcids.Where(c => Hydropathy[c] <= DefaultThreshold))
        {
            PredictTransmembraneHelices(new string(aa, 50)).Should().BeEmpty(
                $"a homopolymer of '{aa}' (h = {Hydropathy[aa]} ≤ 1.6) is not hydrophobic enough to span a membrane");
        }
    }

    #endregion

    #region BE — short: shorter than the window → no region, no crash

    /// <summary>
    /// Target "short": a protein STRICTLY SHORTER than the window (default 19) admits no full window,
    /// so the prediction must yield NO region by the explicit guard
    /// (<c>windowSize ≤ 0 || length &lt; windowSize</c>) and NEVER an IndexOutOfRange at the window edge
    /// (Transmembrane_Helix_Prediction.md §3.3, §6.1, INV-04). We probe lengths 0..windowSize−1 of the
    /// MOST hydrophobic residue (Ile, h = 4.5) — a sequence that WOULD comfortably exceed the threshold
    /// if a window could be formed — to prove it is the length guard, not a low score, that suppresses
    /// the region. The boundary length == windowSize is the first length that CAN yield a region. We
    /// also exercise the <c>windowSize ≤ 0</c> guard and an oversized window.
    /// </summary>
    [Test]
    public void Tm_ShorterThanWindow_NoRegionNoCrash()
    {
        // (a) Every length below the window, of the most hydrophobic residue → no region, no crash.
        for (int len = 0; len < DefaultWindow; len++)
        {
            string seq = new string('I', len); // Ile h = 4.5, would pass the threshold if a window existed
            var act = () => PredictTransmembraneHelices(seq).ToList();
            act.Should().NotThrow($"a length-{len} sequence (< window {DefaultWindow}) must not crash")
                .Subject.Should().BeEmpty(
                    $"a length-{len} sequence is shorter than the window → no transmembrane region (INV-04)");
        }

        // (b) The boundary: exactly windowSize residues of Ile IS long enough for one window → one region.
        var atBoundary = PredictTransmembraneHelices(new string('I', DefaultWindow)).ToList();
        atBoundary.Should().ContainSingle("a length-19 hydrophobic run is exactly one window wide → one region");
        AssertWellFormedRegion(atBoundary[0], DefaultWindow, DefaultWindow, DefaultThreshold);
        atBoundary[0].Should().Be((0, DefaultWindow - 1, 4.5),
            "the single boundary-width window spans residues (0,18) with peak Score = h(Ile) = 4.5");

        // (c) windowSize ≤ 0 guard: even on a long, strongly hydrophobic sequence → no region, no crash
        //     (no divide-by-zero / no profile built).
        foreach (int badWindow in new[] { 0, -1, -19, int.MinValue })
        {
            var act = () => PredictTransmembraneHelices(new string('I', 60), badWindow).ToList();
            act.Should().NotThrow($"windowSize {badWindow} (≤ 0) must be guarded, not crash")
                .Subject.Should().BeEmpty($"windowSize {badWindow} ≤ 0 yields no regions per the guard");
        }

        // (d) Oversized window (> length) → no window can be formed → no region.
        PredictTransmembraneHelices(new string('I', 30), windowSize: 31).Should().BeEmpty(
            "a window wider than the sequence admits no full window → no region");

        // (e) null / empty → no region, no NullReference (INV-04).
        ((Func<List<(int, int, double)>>)(() => PredictTransmembraneHelices(null!).ToList()))
            .Should().NotThrow("a null sequence must be guarded").Subject.Should().BeEmpty();
        PredictTransmembraneHelices("").Should().BeEmpty("an empty sequence yields no region");
    }

    #endregion

    #region BE — non-standard residues: excluded from the mean, never NaN

    /// <summary>
    /// Defends the documented non-standard-residue handling (Transmembrane_Helix_Prediction.md §3.3,
    /// §5.2, §6.1): residues with no scale value (X/B/Z/* and any non-AA char) are EXCLUDED from the
    /// window mean, which divides by the count of scale-bearing residues. The pathological case is a
    /// window composed ENTIRELY of non-standard residues, where a naïve 0/0 would yield NaN — the
    /// source must instead emit 0 for that window and NEVER a NaN Score, and a NaN must never reach a
    /// reported region. We verify (a) an all-'X' sequence yields no region and never a NaN; (b) a
    /// hydrophobic run interrupted by an all-X stretch still produces only finite, in-bounds regions.
    /// </summary>
    [Test]
    public void Tm_AllNonStandardResidues_NoNaNNoRegion()
    {
        // (a) A long all-non-standard sequence: the only residues are off-scale, so every window mean is
        //     the documented 0 (count == 0 → 0, not NaN) → no above-threshold window → no region.
        foreach (string junkResidue in new[] { "X", "B", "Z", "U", "O", "J", "*", "1", "@" })
        {
            string seq = string.Concat(Enumerable.Repeat(junkResidue, 40));
            var regions = ((Func<List<(int Start, int End, double Score)>>)(
                    () => PredictTransmembraneHelices(seq).ToList()))
                .Should().NotThrow($"an all-'{junkResidue}' sequence must not crash").Subject;

            regions.Should().BeEmpty(
                $"an all-'{junkResidue}' window has mean 0 (not NaN, not above threshold) → no region");
            foreach (var r in regions)
                double.IsNaN(r.Score).Should().BeFalse("a window of only non-standard residues must not produce a NaN score");
        }

        // (b) Non-standard residues embedded in a hydrophobic run: any region must remain finite and
        //     in-bounds (the off-scale residues lower/skip the mean but never inject a NaN).
        string mixed = new string('I', 25) + new string('X', 10) + new string('L', 25);
        var mixedRegions = PredictTransmembraneHelices(mixed).ToList();
        foreach (var r in mixedRegions)
            AssertWellFormedRegion(r, mixed.Length, DefaultWindow, DefaultThreshold);
    }

    #endregion

    #region Positive sanity — a flanked hydrophobic stretch yields exactly one documented region

    /// <summary>
    /// Positive sanity (must assert against a predictor that actually FINDS the helix, not a no-op):
    /// the documented worked example — 10 hydrophilic (Asp) residues, a 20-residue hydrophobic (Leu)
    /// stretch, then 10 hydrophilic (Asp) residues — must yield EXACTLY ONE transmembrane region at
    /// the documented span (5, 34) with peak Score = h(Leu) = 3.8 (> 1.6)
    /// (Transmembrane_Helix_Prediction.md §7.1). The numerical walk-through: window 19, threshold 1.6;
    /// the profile has 22 points; means rise above 1.6 first at profile index 5 and last at index 16;
    /// any all-Leu window has mean 3.8; the last passing window (start 16) covers residues 16..34, so
    /// the single run maps to residues (5, 34). We additionally assert the all-hydrophilic flank yields
    /// none and a sub-window length yields none — pinning the three BE outcomes in one realistic shape.
    /// </summary>
    [Test]
    public void Tm_FlankedHydrophobicStretch_SingleDocumentedRegion()
    {
        // Exact §7.1 worked example: D×10 L×20 D×10.
        const string seq = "DDDDDDDDDD" + "LLLLLLLLLLLLLLLLLLLL" + "DDDDDDDDDD";
        seq.Length.Should().Be(40, "sanity-check the hand-built example length");

        var regions = PredictTransmembraneHelices(seq).ToList();
        regions.Should().ContainSingle("a single hydrophobic stretch flanked by hydrophilic residues yields one helix");

        var region = regions[0];
        AssertWellFormedRegion(region, seq.Length, DefaultWindow, DefaultThreshold);
        region.Start.Should().Be(5, "the documented worked example opens at residue 5 (§7.1)");
        region.End.Should().Be(34, "the documented worked example closes at residue 34 (last covered, clamped) (§7.1)");
        region.Score.Should().BeApproximately(3.8, 1e-9,
            "the documented worked example has peak Score = h(Leu) = 3.8 (§7.1)");

        // The hydrophilic flank alone yields nothing.
        PredictTransmembraneHelices(new string('D', 30)).Should().BeEmpty(
            "an all-hydrophilic sequence yields no transmembrane region");

        // A < 19-residue fragment of the hydrophobic core yields nothing (too short for a window).
        PredictTransmembraneHelices(new string('L', 18)).Should().BeEmpty(
            "an 18-residue sequence is shorter than the 19-residue window → no region");
    }

    /// <summary>
    /// Positive sanity over RANDOM proteins: across fixed seeds and lengths the prediction must never
    /// crash, hang, or emit a malformed region, and every emitted segment must satisfy the full
    /// contract (in-bounds clamped coordinates, span ≥ window width, finite peak Score ≥ threshold).
    /// Determinism is pinned by re-running the same input and requiring identical regions. This pins
    /// span-correctness, no-NaN, and termination on arbitrary sequences and a range of window sizes —
    /// not just hand-built homopolymers.
    /// </summary>
    [Test]
    [CancelAfter(30000)]
    public void Tm_RandomProtein_AlwaysWellFormedAndDeterministic(CancellationToken token)
    {
        foreach (int seed in new[] { 11, 53, 211, 2026 })
        {
            foreach (int len in new[] { 1, 18, 19, 60, 250 })
            {
                string seq = RandomProtein(len, seed);

                foreach (int window in new[] { DefaultWindow, 7, 25 })
                {
                    var act = () => PredictTransmembraneHelices(seq, window).ToList();
                    var regions = act.Should().NotThrow(
                        $"random protein must not crash (seed {seed}, len {len}, window {window})").Subject;
                    token.ThrowIfCancellationRequested();

                    foreach (var r in regions)
                        AssertWellFormedRegion(r, len, window, DefaultThreshold);

                    // Determinism: the same input yields identical regions.
                    var again = PredictTransmembraneHelices(seq, window).ToList();
                    again.Should().Equal(regions,
                        $"prediction is deterministic for a fixed input (seed {seed}, len {len}, window {window})");
                }
            }
        }
    }

    #endregion

    #endregion
}
