using static Seqeron.Genomics.Analysis.ProteinMotifFinder;

namespace Seqeron.Genomics.Tests.Fuzzing;

/// <summary>
/// Fuzz tests for the ProteinMotif area — LOW-COMPLEXITY REGION detection
/// (PROTMOTIF-LC-001). The single public entry point under test is
/// <see cref="Seqeron.Genomics.Analysis.ProteinMotifFinder.FindLowComplexityRegions"/>,
/// the SEG algorithm of Wootton &amp; Federhen (1993): a fixed-length window is slid
/// along the (upper-cased) sequence, each window's local complexity is measured as the
/// Shannon entropy of its residue composition in bits per residue
/// (K = −Σ pᵢ·log₂ pᵢ, max log₂(20) ≈ 4.322), and a two-pass trigger/extension rule
/// marks low-complexity segments. Sibling ProteinMotif fuzz units are covered separately:
/// PROTMOTIF-FIND-001 (motif finding), PROTMOTIF-COILEDCOIL-001, PROTMOTIF-COMMON-001;
/// this file focuses exclusively on the SEG low-complexity contract.
///
/// ───────────────────────────────────────────────────────────────────────────
/// What fuzzing verifies
/// ───────────────────────────────────────────────────────────────────────────
/// Fuzzing feeds degenerate, boundary and out-of-domain inputs to a unit and asserts
/// that the code NEVER fails in an undisciplined way: no hang or infinite loop, no
/// state corruption, no nonsense output, and no *unhandled* runtime exception
/// (IndexOutOfRange / NullReference). Every input must resolve to EITHER a well-defined,
/// theory-correct result OR a *documented, intentional* outcome (here the ONLY documented
/// throw is <c>ArgumentOutOfRangeException</c> for windowSize ≤ 0). For a windowed Shannon-
/// entropy scanner the headline hazards are:
///   • a NullReferenceException on a null sequence (must be guarded → empty result);
///   • an IndexOutOfRangeException when the sequence is SHORTER than the window, so no
///     complete window exists (must yield nothing, never run off the end);
///   • a NaN entropy from log₂(0): a residue with count 0 must never enter the −p·log₂ p
///     sum (the implementation visits each DISTINCT residue once), so K stays finite;
///   • a coordinate bug: a reported region [Start..End] must lie inside the sequence,
///     0 ≤ Start ≤ End ≤ n−1, with a finite Complexity in [0, log₂20].
/// — docs/ADVANCED_TESTING_CHECKLIST.md §8 "Fuzzing".
///
/// ───────────────────────────────────────────────────────────────────────────
/// Unit: PROTMOTIF-LC-001 — low-complexity region detection (SEG)
/// Checklist: docs/checklists/03_FUZZING.md, row 165.
/// Fuzz strategy exercised for THIS unit:
///   • BE = Boundary Exploitation — the composition / length / parameter corners that
///     could crash, hang or produce a false positive/negative:
///       – HOMOPOLYMER ("AAAA…"): every window has K = 0 ≤ K1, so the WHOLE tract is a
///         single low-complexity region with Complexity exactly 0 (the headline positive
///         outcome; INV-02, §6.1 "homopolymer tract → single region, complexity 0").
///       – HIGH-COMPLEXITY (a maximally diverse, all-20-residue mosaic where every window
///         holds 12 DISTINCT residues, K = log₂12 ≈ 3.585 > K2 = 2.5): NO region — no false
///         positive (§6.1 "fully diverse sequence → empty result").
///       – EMPTY / null → no regions (guard, no NullReference); a sequence SHORTER than the
///         window → no regions (no complete trigger window exists; no IndexOutOfRange).
///       – windowSize ≤ 0 (0, −1, int.MinValue) → ArgumentOutOfRangeException (§3.3, §6.1).
/// — docs/checklists/03_FUZZING.md §Description (strategy code BE);
///   targets: "homopolymer, high-complexity, empty".
///
/// ───────────────────────────────────────────────────────────────────────────
/// The SEG low-complexity contract under test
/// ───────────────────────────────────────────────────────────────────────────
/// Given a sequence S and parameters (W = windowSize, K1 = triggerComplexity,
/// K2 = extensionComplexity), for every window-start position compute
/// K = −Σᵢ pᵢ·log₂ pᵢ over the window's residue composition (pᵢ = countᵢ / W).
/// Group maximal runs of consecutive windows with K ≤ K2; emit a run as a region only if
/// at least one of its windows has K ≤ K1 (triggered). A region is reported as
/// (Start, End, Complexity) where Start = firstWindowStart, End = lastWindowStart + W − 1
/// (0-based inclusive), and Complexity = the MINIMUM window complexity inside the run.
///   — docs/algorithms/ProteinMotif/Low_Complexity_Region_Detection.md §2.2, §3, §4.
///
/// Method under test (src/.../Seqeron.Genomics.Analysis/ProteinMotifFinder.cs):
///   IEnumerable&lt;(int Start,int End,double Complexity)&gt; FindLowComplexityRegions(
///       string proteinSequence, int windowSize = 12,
///       double triggerComplexity = 2.2, double extensionComplexity = 2.5)
///   — §5.1. Default parameters W=12, K1=2.2, K2=2.5 (NCBI SEG man page / blast_seg.c).
///
/// Documented input handling (Low_Complexity_Region_Detection.md §3.3, §6.1):
///   • null / "" / length &lt; windowSize → EMPTY result (no complete window), no throw.
///   • windowSize ≤ 0 → ArgumentOutOfRangeException.
///   • case-INSENSITIVE: the sequence is upper-cased before counting; residue identity is
///     taken literally, so any non-standard char is simply an extra symbol in the alphabet.
///
/// Theory-correct invariants asserted (Low_Complexity_Region_Detection.md §2.4):
///   • INV-01 — 0 ≤ K ≤ log₂(20) for amino-acid windows (Shannon-entropy bound).
///   • INV-02 — a homopolymer window has K = 0.
///   • INV-03 — a reported region contains ≥1 window with K ≤ K1 and every window K ≤ K2;
///     therefore the reported Complexity (the run minimum) is ≤ K2.
///   • INV-04 — region boundaries are 0-based inclusive and lie within the sequence.
///   • INV-05 — output is deterministic, ordered by Start (single left-to-right scan).
///
/// ───────────────────────────────────────────────────────────────────────────
/// Complexity / hang-safety
/// ───────────────────────────────────────────────────────────────────────────
/// The scan is O(n·W): n−W+1 windows, each O(W) to tally composition (§4.3). The
/// homopolymer and long-mosaic targets are kept modest and [CancelAfter]-guarded so a
/// regression turning the linear scan into a hang would FAIL rather than wedge the suite.
/// ───────────────────────────────────────────────────────────────────────────
/// </summary>
[TestFixture]
[Category("Fuzzing")]
public class ProteinLowComplexityFuzzTests
{
    #region Helpers

    /// <summary>The 20 standard amino-acid one-letter codes.</summary>
    private const string StandardAminoAcids = "ACDEFGHIKLMNPQRSTVWY";

    /// <summary>Default SEG parameters under test (NCBI SEG man page / blast_seg.c).</summary>
    private const int DefaultWindow = 12;
    private const double DefaultK1 = 2.2;
    private const double DefaultK2 = 2.5;

    /// <summary>Maximum possible complexity for the 20-letter amino-acid alphabet: log₂(20).</summary>
    private static readonly double MaxAaComplexity = Math.Log2(20);

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
    /// A maximally diverse mosaic: residue i is the (i mod 20)-th standard amino acid. Every window
    /// of width ≥ 20 holds all 20 residues; every window of width 12 holds 12 DISTINCT residues, so
    /// K = log₂(12) ≈ 3.585 &gt; K2 = 2.5 everywhere — a high-complexity sequence with NO LCR.
    /// </summary>
    private static string DiverseProtein(int length)
    {
        var chars = new char[length];
        for (int i = 0; i < length; i++)
            chars[i] = StandardAminoAcids[i % StandardAminoAcids.Length];
        return new string(chars);
    }

    /// <summary>
    /// Asserts the universal theory-correct contract every emitted region must satisfy against a
    /// sequence of length <paramref name="n"/> with extension cutoff <paramref name="k2"/>
    /// (Low_Complexity_Region_Detection.md §2.4, §3.2):
    ///   • INV-04 — in-bounds, ordered span: 0 ≤ Start ≤ End ≤ n−1.
    ///   • INV-01 — the reported Complexity (a window minimum) is finite and in [0, log₂20].
    ///   • INV-03 — the reported Complexity is ≤ the extension cutoff K2 (every window in a reported
    ///     run satisfies K ≤ K2, so the run minimum does too) — the headline "no NaN, no out-of-range,
    ///     and the score reflects an actually-low-complexity run" property.
    /// </summary>
    private static void AssertWellFormedRegion((int Start, int End, double Complexity) region, int n, double k2)
    {
        region.Start.Should().BeInRange(0, n - 1, "INV-04: a region Start is a valid 0-based residue index");
        region.End.Should().BeInRange(region.Start, n - 1, "INV-04: a region End is in-bounds and not before its Start");

        double.IsNaN(region.Complexity).Should().BeFalse("INV-01: a region Complexity must never be NaN (no log₂(0))");
        double.IsInfinity(region.Complexity).Should().BeFalse("INV-01: a region Complexity must never be infinite");
        region.Complexity.Should().BeGreaterThanOrEqualTo(0.0, "INV-01: Shannon entropy is non-negative");
        region.Complexity.Should().BeLessThanOrEqualTo(MaxAaComplexity + 1e-9,
            "INV-01: amino-acid window entropy is bounded by log₂(20)");
        region.Complexity.Should().BeLessThanOrEqualTo(k2 + 1e-9,
            "INV-03: a reported region's minimum window complexity is ≤ the extension cutoff K2");
    }

    #endregion

    // ═══════════════════════════════════════════════════════════════════
    //  PROTMOTIF-LC-001 — low-complexity region detection (SEG) : fuzz targets
    // ═══════════════════════════════════════════════════════════════════

    #region PROTMOTIF-LC-001 — low-complexity region detection

    #region BE — Empty / null / shorter-than-window: no regions, no crash

    /// <summary>
    /// Targets "empty" and the short-input boundary: a null or empty sequence, and any sequence
    /// STRICTLY SHORTER than the window, must yield NO regions — by the explicit guard, never a
    /// NullReferenceException and never an IndexOutOfRangeException on the absent window
    /// (Low_Complexity_Region_Detection.md §3.3, §6.1 "null/empty → empty result";
    /// "length &lt; windowSize → empty result, no complete trigger window"). We probe null, "",
    /// every length 1..W−1 of a homopolymer (which WOULD be low-complexity if it were long enough),
    /// and the exact W−1 boundary, all of which must come back empty without throwing.
    /// </summary>
    [Test]
    public void FindLcr_EmptyNullOrShorterThanWindow_NoRegionsNoThrow()
    {
        foreach (string? seq in new[] { null, "" })
        {
            var act = () => FindLowComplexityRegions(seq!).ToList();
            act.Should().NotThrow($"null/empty sequence ('{seq ?? "null"}') must not crash")
                .Subject.Should().BeEmpty("a null/empty sequence has no windows, so no regions");
        }

        // Every length below the window — even a perfect homopolymer — has no complete window → empty.
        for (int len = 1; len < DefaultWindow; len++)
        {
            string homo = new string('A', len);
            var act = () => FindLowComplexityRegions(homo).ToList();
            act.Should().NotThrow($"a length-{len} sequence (< window {DefaultWindow}) must not crash")
                .Subject.Should().BeEmpty($"a length-{len} sequence is shorter than the window, so no complete window exists");
        }

        // Exactly W−1 of a homopolymer: still one short of a complete window → empty.
        FindLowComplexityRegions(new string('Q', DefaultWindow - 1)).Should().BeEmpty(
            "a sequence one residue shorter than the window yields no regions");

        // Exactly W of a homopolymer is the smallest sequence that CAN report (covered positively below);
        // here we only assert the W−1/W boundary does not throw.
        var boundary = () => FindLowComplexityRegions(new string('Q', DefaultWindow)).ToList();
        boundary.Should().NotThrow("a sequence of exactly the window length must not crash");
    }

    #endregion

    #region BE — windowSize ≤ 0: ArgumentOutOfRangeException

    /// <summary>
    /// Target boundary parameter: a non-positive window (0, −1, int.MinValue) is invalid and MUST
    /// throw <c>ArgumentOutOfRangeException</c> — the ONLY documented throw of this unit
    /// (Low_Complexity_Region_Detection.md §3.3, §6.1 "windowSize ≤ 0 → ArgumentOutOfRangeException").
    /// Because the method is an iterator (deferred execution), the throw must surface when the
    /// sequence is enumerated; we force enumeration with ToList(). The validation must fire even for
    /// a null/empty sequence (the window check precedes the empty guard in the contract).
    /// </summary>
    [Test]
    public void FindLcr_NonPositiveWindow_Throws()
    {
        foreach (int badWindow in new[] { 0, -1, -12, int.MinValue })
        {
            var act = () => FindLowComplexityRegions("ACDEFGHIKLMNPQRSTVWY", badWindow).ToList();
            act.Should().Throw<ArgumentOutOfRangeException>(
                    $"a non-positive window ({badWindow}) is invalid")
                .Which.ParamName.Should().Be("windowSize", "the offending parameter is identified");

            // Validation fires before the empty-sequence guard, so even null/empty must throw.
            var actEmpty = () => FindLowComplexityRegions("", badWindow).ToList();
            actEmpty.Should().Throw<ArgumentOutOfRangeException>(
                "windowSize validation precedes the empty-sequence guard");
        }
    }

    #endregion

    #region BE — Homopolymer: whole tract is one region, Complexity exactly 0

    /// <summary>
    /// Target "homopolymer" — the headline positive outcome. A homopolymer "AAAA…" of length n ≥ W
    /// has K = 0 in EVERY window (a single symbol has p = 1, −1·log₂1 = 0; INV-02), so every window
    /// triggers (0 ≤ K1) and the entire run merges into ONE region spanning the whole tract
    /// [0, n−1] with Complexity exactly 0 (the run minimum). We assert this for every standard
    /// residue and several lengths, and pin the documented poly-Q worked example
    /// (Low_Complexity_Region_Detection.md §6.1, §7.1).
    /// </summary>
    [Test]
    public void FindLcr_Homopolymer_SingleRegionWholeTractComplexityZero()
    {
        foreach (char aa in StandardAminoAcids)
        {
            foreach (int n in new[] { DefaultWindow, DefaultWindow + 1, 20, 50 })
            {
                string homo = new string(aa, n);
                var regions = FindLowComplexityRegions(homo).ToList();

                regions.Should().ContainSingle($"a length-{n} homopolymer of '{aa}' is one merged low-complexity region");
                var r = regions[0];
                AssertWellFormedRegion(r, n, DefaultK2);
                r.Start.Should().Be(0, "the homopolymer region starts at residue 0");
                r.End.Should().Be(n - 1, "the homopolymer region spans to the last residue (run = [0, lastStart + W − 1] = [0, n−1])");
                r.Complexity.Should().BeApproximately(0.0, 1e-12, "INV-02: a homopolymer window has Shannon entropy 0");
            }
        }
    }

    /// <summary>
    /// Documented poly-Q worked example (Low_Complexity_Region_Detection.md §7.1): a 20-residue
    /// glutamine tract flanked by diverse 8-residue segments yields ONE region whose core is the
    /// homopolymer (Complexity ≈ 0). The reported span must cover the Q tract and reproduce a
    /// well-formed region; the surrounding "MKLPRDST" flanks are diverse and must not themselves
    /// trigger a separate spurious LCR.
    /// </summary>
    [Test]
    public void FindLcr_PolyQTractInDiverseFlanks_ReportsHomopolymerCoreComplexityZero()
    {
        const string flank = "MKLPRDST";
        string seq = flank + new string('Q', 20) + flank;
        int qStart = flank.Length;          // 8
        int qEnd = qStart + 20 - 1;          // 27

        var regions = FindLowComplexityRegions(seq).ToList();

        regions.Should().NotBeEmpty("the poly-Q tract is a low-complexity region");
        foreach (var r in regions)
            AssertWellFormedRegion(r, seq.Length, DefaultK2);

        // Exactly one region whose minimum complexity is the homopolymer core (≈ 0) and which
        // covers the Q tract.
        regions.Should().ContainSingle("a single LCR spans the poly-Q tract amid diverse flanks");
        var region = regions[0];
        region.Complexity.Should().BeApproximately(0.0, 1e-12,
            "the homopolymer core contributes the minimum window complexity 0");
        region.Start.Should().BeLessThanOrEqualTo(qStart, "the region begins no later than the Q tract start");
        region.End.Should().BeGreaterThanOrEqualTo(qEnd, "the region extends at least to the Q tract end");
    }

    #endregion

    #region BE — High-complexity: maximally diverse sequence yields NO region (no false positive)

    /// <summary>
    /// Target "high-complexity" — the headline NEGATIVE outcome. A maximally diverse mosaic where
    /// residue i = (i mod 20)-th amino acid puts 12 DISTINCT residues in every width-12 window, so
    /// K = log₂(12) ≈ 3.585 &gt; K2 = 2.5 everywhere: NO window can trigger or extend, hence NO region
    /// (no false positive) — Low_Complexity_Region_Detection.md §6.1 "fully diverse sequence →
    /// empty result". We assert emptiness across several lengths AND confirm the underlying premise
    /// (the per-window complexity genuinely exceeds K2) is consistent, so the empty result is a true
    /// negative rather than a vacuous one.
    /// </summary>
    [Test]
    public void FindLcr_HighComplexityDiverseSequence_NoRegions()
    {
        foreach (int n in new[] { DefaultWindow, 20, 40, 100, 200 })
        {
            string diverse = DiverseProtein(n);
            var regions = FindLowComplexityRegions(diverse).ToList();
            regions.Should().BeEmpty(
                $"a maximally diverse length-{n} sequence has every window complexity ≈ log₂(12) > K2, so no LCR");
        }

        // Premise check: a single diverse window's complexity is log₂(12) ≈ 3.585 > K2 = 2.5, and
        // raising K2 above that threshold makes the WHOLE diverse sequence one trivially-extended
        // (but still untriggered, since K > K1) run → still no region. Confirms the negative is real:
        // even with K2 = 4 (above log₂12), nothing triggers because K > K1 = 2.2.
        string diverse40 = DiverseProtein(40);
        FindLowComplexityRegions(diverse40, DefaultWindow, triggerComplexity: DefaultK1, extensionComplexity: 4.0)
            .Should().BeEmpty("even extending over every diverse window, none triggers (K ≈ 3.585 > K1 = 2.2)");

        // And if BOTH cutoffs are pushed above log₂(12), the whole diverse sequence DOES become one
        // region — proving the prior emptiness was due to the cutoffs, not a no-op scanner.
        var forced = FindLowComplexityRegions(diverse40, DefaultWindow, triggerComplexity: 4.0, extensionComplexity: 4.0).ToList();
        forced.Should().ContainSingle("with both cutoffs above log₂(12) the diverse sequence is one big 'low-complexity' run");
        AssertWellFormedRegion(forced[0], diverse40.Length, 4.0);
        forced[0].Start.Should().Be(0);
        forced[0].End.Should().Be(diverse40.Length - 1, "the forced run spans the whole sequence");
    }

    #endregion

    #region BE — NaN / log₂(0) guard and case-insensitivity on arbitrary input

    /// <summary>
    /// Headline no-NaN property over arbitrary and adversarial composition: the per-window entropy
    /// must NEVER be NaN (a count-0 residue must never enter the −p·log₂ p sum — the implementation
    /// visits each DISTINCT residue once). We feed random proteins, biased mixtures, and out-of-
    /// alphabet junk (digits, punctuation, X) across seeds/lengths; every emitted region must be
    /// well-formed (finite Complexity in [0, log₂20], ≤ K2), the scan must terminate, and re-running
    /// must give the IDENTICAL ordered result (INV-05 determinism). Lowercase input must produce the
    /// SAME regions as its uppercase form (case-insensitive, §3.3).
    /// </summary>
    [Test]
    [CancelAfter(30000)]
    public void FindLcr_ArbitraryInput_AlwaysFiniteWellFormedDeterministic(CancellationToken token)
    {
        var inputs = new List<string>();
        foreach (int seed in new[] { 7, 31, 137, 2026 })
            foreach (int len in new[] { 12, 25, 60, 150 })
                inputs.Add(RandomProtein(len, seed));

        // Biased / adversarial compositions (low-complexity runs, junk, mixed alphabets).
        inputs.Add(new string('A', 30) + RandomProtein(30, 99));            // homopolymer + random
        inputs.Add("MK1RGD2SP" + new string('S', 20) + "!@#$%^&*()");        // junk + poly-S
        inputs.Add("XXXXXXXXXXXXXXXX");                                       // out-of-alphabet placeholder run
        inputs.Add("ABABABABABABABABABAB");                                   // period-2 repeat (high entropy by design)
        inputs.Add("   \t  \n  whitespace and text mixed in here  ");         // whitespace + literals

        foreach (string seq in inputs)
        {
            var act = () => FindLowComplexityRegions(seq).ToList();
            var regions = act.Should().NotThrow($"arbitrary input must not crash: '{seq[..Math.Min(seq.Length, 20)]}'").Subject;
            token.ThrowIfCancellationRequested();

            foreach (var r in regions)
                AssertWellFormedRegion(r, seq.Length, DefaultK2);

            // INV-05 — deterministic, ordered by Start.
            regions.Select(r => r.Start).Should().BeInAscendingOrder("INV-05: regions are ordered by Start");
            var again = FindLowComplexityRegions(seq).ToList();
            again.Should().Equal(regions, "INV-05: the scan is deterministic for a fixed input");
        }

        // Case-insensitivity: lowercase must give identical regions to uppercase.
        string mixed = "mklprdst" + new string('q', 20) + "MKLPRDST";
        var lower = FindLowComplexityRegions(mixed).ToList();
        var upper = FindLowComplexityRegions(mixed.ToUpperInvariant()).ToList();
        lower.Should().Equal(upper, "matching is case-insensitive: lowercase input yields the same regions as uppercase");
    }

    #endregion

    #region Positive sanity — a genuinely low-complexity run is found at the right span

    /// <summary>
    /// Positive sanity: the harness must assert against a scanner that actually FINDS low-complexity
    /// regions at the correct span, not a no-op. A long high-complexity sequence with an embedded
    /// homopolymer run long enough to fill the window must yield exactly one region whose minimum
    /// Complexity is ≈ 0 and whose span covers the run; and the documented 11-A/1-B numerical
    /// walk-through (K ≈ 0.4138 &lt; K1 = 2.2) must trigger when that biased window is present.
    /// (Low_Complexity_Region_Detection.md §2.4 INV-03/04, §7.1.)
    /// </summary>
    [Test]
    public void FindLcr_LowComplexityRunInDiverseFlanks_FoundAtCorrectSpan()
    {
        // A diverse 30-residue prefix, then a 16-A homopolymer, then a diverse 30-residue suffix.
        string diversePrefix = DiverseProtein(30);
        string diverseSuffix = DiverseProtein(30);
        const int runLen = 16;
        int runStart = diversePrefix.Length;          // 30
        int runEnd = runStart + runLen - 1;            // 45
        string seq = diversePrefix + new string('A', runLen) + diverseSuffix;

        var regions = FindLowComplexityRegions(seq).ToList();
        regions.Should().ContainSingle("exactly one homopolymer LCR sits amid diverse flanks");
        var region = regions[0];
        AssertWellFormedRegion(region, seq.Length, DefaultK2);
        region.Complexity.Should().BeApproximately(0.0, 1e-12, "the homopolymer core minimum complexity is 0");
        region.Start.Should().BeLessThanOrEqualTo(runStart, "the region begins no later than the homopolymer run");
        region.End.Should().BeGreaterThanOrEqualTo(runEnd, "the region extends at least to the homopolymer run end");

        // Documented numerical walk-through: an 11-A/1-B biased window has K ≈ 0.4138 < K1, so a
        // sequence whose only sub-window-complexity dip is such a biased composition still triggers.
        // Build a window of exactly 11 A + 1 'B'-like residue (use 'C') as a standalone W-length seq.
        string biasedWindow = new string('A', 11) + "C";   // length 12
        var biasedRegions = FindLowComplexityRegions(biasedWindow).ToList();
        biasedRegions.Should().ContainSingle("an 11:1 biased window (K ≈ 0.4138 < K1 = 2.2) triggers a low-complexity region");
        biasedRegions[0].Complexity.Should().BeApproximately(0.4138168503, 1e-9,
            "the documented worked-example complexity for an 11-A/1-other window");
        biasedRegions[0].Start.Should().Be(0);
        biasedRegions[0].End.Should().Be(biasedWindow.Length - 1, "the single window spans the whole 12-residue sequence");
    }

    #endregion

    #endregion
}
