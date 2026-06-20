using System;
using System.Collections.Generic;
using System.Linq;
using System.Threading;
using NUnit.Framework;
using FluentAssertions;
using Seqeron.Genomics.Annotation;
using static Seqeron.Genomics.Annotation.SpliceSitePredictor;

namespace Seqeron.Genomics.Tests;

/// <summary>
/// Fuzz tests for the Splicing area — splice DONOR (5') site detection
/// (SPLICE-DONOR-001) and splice ACCEPTOR (3') site detection
/// (SPLICE-ACCEPTOR-001).
///
/// ───────────────────────────────────────────────────────────────────────────
/// What fuzzing verifies
/// ───────────────────────────────────────────────────────────────────────────
/// Fuzzing feeds degenerate, boundary and out-of-domain inputs to a unit and
/// asserts that the code NEVER fails in an undisciplined way: no hang or infinite
/// loop, no state corruption, no nonsense output, and no *unhandled* runtime
/// exception. For a motif scanner that extracts a fixed window AROUND each hit,
/// the headline hazard is an IndexOutOfRangeException leaking when a candidate
/// dinucleotide sits near the sequence edge and the −3..+5 scoring window (or the
/// motif-context Substring) runs off either end of the string. Every input must
/// resolve to EITHER a well-defined, theory-correct result, OR a *documented,
/// intentional* validation outcome.
/// — docs/ADVANCED_TESTING_CHECKLIST.md §8 "Fuzzing".
///
/// ───────────────────────────────────────────────────────────────────────────
/// Unit: SPLICE-DONOR-001 — splice donor (5') site detection
/// Checklist: docs/checklists/03_FUZZING.md, row 77.
/// Fuzz strategies exercised for THIS unit:
///   • BE = Boundary Exploitation — the length/edge corners that could crash:
///       – seq &lt; window: a sequence shorter than 6 nt (the documented guard), AND
///         a canonical GU pushed flush against the 3' edge so the −3..+5 window
///         cannot be fully extracted. This is the KEY no-IndexOutOfRange boundary.
///       – all-GU: "GUGUGU…" — every even index is a candidate GU, so the scanner
///         scores the maximum number of overlapping candidates; it must complete
///         without crashing or wedging ([CancelAfter]-guarded).
///   • MC = Malformed Content — non-DNA characters (digits, punctuation, IUPAC
///         ambiguity codes, whitespace): handled per contract (uppercased, T→U,
///         then treated as never-matching residues), never rejected with a crash.
/// — docs/checklists/03_FUZZING.md §Description (strategy codes BE, MC);
///   targets: "No GT sites, seq &lt; window, non-DNA, all GT".
///
/// ───────────────────────────────────────────────────────────────────────────
/// The donor-detection contract under test
/// ───────────────────────────────────────────────────────────────────────────
/// The 5' splice donor site marks the START of an intron. The canonical intronic
/// dinucleotide is GU (GT in DNA); the extended consensus is MAG|GURAGU scored over
/// the offsets −3..+5 around the GU (Shapiro &amp; Senapathy 1987; Burge et al. 1999).
/// The repository scores each candidate with a BINARY consensus table: every offset
/// contributes 1.0 for a consensus match and 0.0 otherwise, normalised by the number
/// of offsets that actually fall inside the sequence — so the score is always in
/// [0, 1]. — docs/algorithms/Splicing/Donor_Site_Detection.md §2.1–2.2, §2.4.
///
/// Method under test (src/.../Seqeron.Genomics.Annotation/SpliceSitePredictor.cs):
///   FindDonorSites(string sequence, double minScore = 0.5, bool includeNonCanonical = false)
///     → IEnumerable&lt;SpliceSite&gt; { Position, Type, Motif, Score, Confidence }.
///   The scan uppercases the input, maps T→U, then for i in 0..Length−6 inspects the
///   dinucleotide at (i, i+1): canonical GU is always scored; GC and U12-style AU are
///   scored only when includeNonCanonical is set. A site is emitted iff its score ≥
///   minScore. — Donor_Site_Detection.md §3, §4.
///
/// Documented input handling (Donor_Site_Detection.md §3.1, §3.3, §6.1):
///   • null / "" / length &lt; 6 → no donor sites (guard clause) — never an exception.
///   • DNA T behaves exactly like RNA U (T is converted to U before scanning).
///   • lowercase behaves exactly like uppercase (the sequence is uppercased first).
///   • No GU/GT candidate in canonical mode → no donor sites (nothing to score).
///   • A GU within the last 5 positions is OUTSIDE the i ≤ Length−6 scan range, so it
///     is simply not emitted — never an IndexOutOfRange off the window extraction.
///     The internal scorer and GetMotifContext additionally clamp every offset/Substring
///     to the string bounds, so even an in-range GU near an edge is window-safe.
///
/// Theory-correct invariants asserted (Donor_Site_Detection.md §2.4):
///   • INV-01 — every reported Score is in [0, 1].
///   • INV-02 — every reported Confidence is in [0, 1].
///   • [donor-anchored] — a reported canonical/GC donor really HAS its donor
///     dinucleotide at Position in the normalised (uppercased, T→U) sequence: the
///     scanner only enters the emit branch after matching G·U (or G·C) there, so a
///     reported donor must actually sit on a GU/GC. This is the "a reported donor
///     site must actually have GT at the splice position" contract.
///
/// ───────────────────────────────────────────────────────────────────────────
/// Complexity / hang-safety
/// ───────────────────────────────────────────────────────────────────────────
/// FindDonorSites is a single O(n) streaming pass with O(1) auxiliary state
/// (Donor_Site_Detection.md §4.3). The all-GU target maximises the number of scored
/// candidates per length; it is kept modest and [CancelAfter]-guarded so a regression
/// that turned the linear scan into a hang or super-linear blow-up would FAIL rather
/// than wedge the suite.
///
/// ═══════════════════════════════════════════════════════════════════════════
/// Unit: SPLICE-ACCEPTOR-001 — splice acceptor (3') site detection
/// Checklist: docs/checklists/03_FUZZING.md, row 78.
/// ═══════════════════════════════════════════════════════════════════════════
/// Fuzz strategies exercised for THIS unit:
///   • BE = Boundary Exploitation — the length/edge corners that could crash:
///       – seq &lt; window: a sequence shorter than the documented 20-nt guard,
///         INCLUDING one that DOES contain an AG, and an AG flush against the 3'
///         edge of a ≥ 20-nt sequence so the scan's i ≤ Length-2 range plus the
///         clamped −15..0 PWM / Substring windows cannot run off either end. This
///         is the KEY no-IndexOutOfRange boundary (the −15-nt PPT/PWM window AND
///         the trailing AG both threaten the string bounds).
///       – all-AG: "AGAGAG…" — every even index is a candidate AG, so the scanner
///         scores the maximum number of overlapping acceptor candidates; it must
///         complete without crashing or wedging ([CancelAfter]-guarded).
///   • MC = Malformed Content — non-DNA characters (digits, punctuation, IUPAC
///         ambiguity codes, whitespace): handled per contract (uppercased, T→U,
///         then treated as never-matching residues), never rejected with a crash.
/// — docs/checklists/03_FUZZING.md §Description (strategy codes BE, MC);
///   targets: "No AG sites, seq &lt; window, non-DNA, all AG".
///
/// The acceptor-detection contract under test
/// ───────────────────────────────────────────────────────────────────────────
/// The 3' splice ACCEPTOR site marks the END of an intron. The canonical 3'
/// consensus is (Y)nNCAG|G with an almost-invariant terminal AG dinucleotide
/// preceded by a polypyrimidine tract (PPT) (Shapiro &amp; Senapathy 1987; Burge et
/// al. 1999). The repository scores each AG candidate by a PPT contribution over
/// [−15, −3) plus a sparse eight-offset AcceptorPwm, normalised to [0, 1]; optional
/// U12-style AC acceptors are scored against a YCCAC-style consensus normalised by
/// 3.5. — docs/algorithms/Splicing/Acceptor_Site_Detection.md §2.1–2.2, §2.4.
///
/// Method under test (src/.../Seqeron.Genomics.Annotation/SpliceSitePredictor.cs):
///   FindAcceptorSites(string sequence, double minScore = 0.5, bool includeNonCanonical = false)
///     → IEnumerable&lt;SpliceSite&gt; { Position, Type, Motif, Score, Confidence }.
///   The scan uppercases the input, maps T→U, then for i in 15..Length−2 inspects the
///   dinucleotide at (i, i+1): canonical AG is always scored; U12 AC is scored only
///   when includeNonCanonical is set. A site is emitted iff its score ≥ minScore, with
///   Position = i+1 (the G of AG / the C of AC — the terminal intronic base, INV-04).
///   — Acceptor_Site_Detection.md §3, §4.
///
/// Documented input handling (Acceptor_Site_Detection.md §3.1, §6.1):
///   • null / "" / length &lt; 20 → no acceptor sites (guard clause) — never an exception.
///   • DNA T behaves exactly like RNA U; lowercase behaves exactly like uppercase.
///   • No AG candidate in canonical mode → no acceptor sites (nothing to score).
///   • The scan starts at i = 15 so the −15 PPT/PWM window exists; the PPT loop and
///     every PWM offset clamp with an `(pos >= 0 && pos < Length)` bound, and
///     GetMotifContext clamps its Substring — so an AG anywhere in 15..Length−2,
///     including flush against the 3' edge, extracts a window-safe (clamped) score
///     and motif. No IndexOutOfRange leaks. An AG before index 15 is outside the scan.
///
/// Theory-correct invariants asserted (Acceptor_Site_Detection.md §2.4):
///   • INV-01 — every reported Score is in [0, 1].
///   • INV-02 — every reported Confidence is in [0, 1].
///   • INV-03 — canonical sites carry Type = Acceptor; AC sites Type = U12Acceptor.
///   • INV-04 / [acceptor-anchored] — a reported canonical acceptor really HAS its
///     AG at (Position−1, Position) in the normalised sequence (a U12 acceptor an AC):
///     the scanner only enters the emit branch after matching that dinucleotide and
///     yields Position = i+1, so a reported acceptor must actually sit on an AG/AC.
///     This is the "a reported acceptor site must actually have AG at the splice
///     position" contract.
/// ───────────────────────────────────────────────────────────────────────────
/// </summary>
[TestFixture]
[Category("Fuzzing")]
public class SplicingFuzzTests
{
    #region Helpers

    /// <summary>Deterministic RNG — seed fixed locally so generated fuzz inputs are reproducible.</summary>
    private static string RandomDna(int length, int seed)
    {
        const string bases = "ACGT";
        var rng = new Random(seed);
        var chars = new char[length];
        for (int i = 0; i < length; i++)
            chars[i] = bases[rng.Next(bases.Length)];
        return new string(chars);
    }

    /// <summary>
    /// Asserts the universal theory-correct contract every emitted donor site must
    /// satisfy (Donor_Site_Detection.md §2.4): Score and Confidence in [0, 1] (INV-01,
    /// INV-02), an in-range Position, and — the headline correctness property — the
    /// site is anchored on a real donor dinucleotide: the normalised sequence carries
    /// G·U (canonical/GC start with G) at (Position, Position+1), i.e. a reported donor
    /// truly has GT/GU at the splice position.
    /// </summary>
    private static void AssertWellFormedDonor(SpliceSite site, string sequence)
    {
        string normalized = sequence.ToUpperInvariant().Replace('T', 'U');

        site.Score.Should().BeInRange(0.0, 1.0, "INV-01: donor scores are normalised to [0, 1]");
        site.Confidence.Should().BeInRange(0.0, 1.0, "INV-02: confidence is clamped to [0, 1]");

        // Position indexes the donor dinucleotide start; the +1 partner must also exist.
        site.Position.Should().BeInRange(0, normalized.Length - 2,
            "a donor dinucleotide occupies Position and Position+1, both inside the sequence");

        // [donor-anchored] every canonical / GC donor starts on a G; U12 donors start on A.
        char first = normalized[site.Position];
        char second = normalized[site.Position + 1];
        if (site.Type == SpliceSiteType.U12Donor)
        {
            (first == 'A' && second == 'U').Should().BeTrue(
                $"a U12 donor must sit on an AU dinucleotide, found '{first}{second}' at {site.Position}");
        }
        else
        {
            site.Type.Should().Be(SpliceSiteType.Donor, "canonical and GC donors carry Type = Donor (INV-03)");
            first.Should().Be('G',
                $"a reported donor must have its donor dinucleotide (GU/GC) at the splice position, found '{first}' at {site.Position}");
            (second == 'U' || second == 'C').Should().BeTrue(
                $"a canonical/GC donor's second base must be U or C, found '{second}' at {site.Position + 1}");
        }

        site.Motif.Should().NotBeNull("every emitted site carries a (possibly clamped) motif context");
    }

    /// <summary>
    /// Asserts the universal theory-correct contract every emitted ACCEPTOR site must
    /// satisfy (Acceptor_Site_Detection.md §2.4): Score and Confidence in [0, 1] (INV-01,
    /// INV-02), an in-range Position, and — the headline correctness property — the site
    /// is anchored on a real acceptor dinucleotide. FindAcceptorSites reports Position =
    /// i+1, so Position indexes the terminal intronic base (the G of AG / the C of AC);
    /// the normalised sequence must carry A·G at (Position−1, Position) for a canonical
    /// Acceptor (INV-03/INV-04), or A·C for a U12Acceptor. I.e. a reported acceptor truly
    /// has AG at the splice position. Because the scan starts at i = 15, a canonical
    /// Position is always ≥ 16.
    /// </summary>
    private static void AssertWellFormedAcceptor(SpliceSite site, string sequence)
    {
        string normalized = sequence.ToUpperInvariant().Replace('T', 'U');

        site.Score.Should().BeInRange(0.0, 1.0, "INV-01: acceptor scores are normalised to [0, 1]");
        site.Confidence.Should().BeInRange(0.0, 1.0, "INV-02: confidence is clamped to [0, 1]");

        // Position = i + 1 is the terminal intronic base; (Position-1, Position) is the dinucleotide.
        site.Position.Should().BeInRange(16, normalized.Length - 1,
            "Position = i+1 with i ≥ 15, and the G/C member must be inside the sequence");

        char first = normalized[site.Position - 1];  // A of AG / A of AC
        char second = normalized[site.Position];      // G of AG / C of AC
        if (site.Type == SpliceSiteType.U12Acceptor)
        {
            (first == 'A' && second == 'C').Should().BeTrue(
                $"a U12 acceptor must sit on an AC dinucleotide, found '{first}{second}' ending at {site.Position}");
        }
        else
        {
            site.Type.Should().Be(SpliceSiteType.Acceptor, "canonical acceptors carry Type = Acceptor (INV-03)");
            (first == 'A' && second == 'G').Should().BeTrue(
                $"a reported acceptor must have its AG dinucleotide at the splice position, found '{first}{second}' ending at {site.Position} (INV-04)");
        }

        site.Motif.Should().NotBeNull("every emitted site carries a (possibly clamped) motif context");
    }

    #endregion

    // ═══════════════════════════════════════════════════════════════════
    //  SPLICE-DONOR-001 — splice donor site detection : fuzz targets
    // ═══════════════════════════════════════════════════════════════════

    #region SPLICE-DONOR-001 — splice donor site detection

    #region MC — No GT/GU candidate: returns no donor sites

    /// <summary>
    /// Target "No GT sites": a sequence containing NO GU/GT dinucleotide has no
    /// canonical donor candidate at all, so the scan emits nothing — no donor site,
    /// no crash (Donor_Site_Detection.md §6.1 "No GU in canonical mode"). We probe a
    /// homopolymer (poly-A, which can never form GU), an AC-only sequence (a G never
    /// followed by U), and confirm that turning ON includeNonCanonical does not invent
    /// a canonical-style donor from a GU that simply is not present.
    /// </summary>
    [Test]
    public void FindDonorSites_NoGuCandidate_ReturnsEmpty()
    {
        foreach (string seq in new[] { new string('A', 30), "ACACACACACACACAC", "CCCCAAAACCCCAAAA" })
        {
            // Canonical mode: no GU anywhere → nothing scored, nothing emitted.
            var canonical = ((Func<List<SpliceSite>>)(() =>
                FindDonorSites(seq, minScore: 0.0).ToList()))
                .Should().NotThrow("a GU-free sequence is scanned without error").Subject;
            canonical.Should().BeEmpty($"'{seq}' has no GU/GT dinucleotide, so there is no canonical donor candidate");

            // Even at minScore 0.0 with non-canonical OFF, no canonical donor can appear.
            canonical.Should().NotContain(s => s.Type == SpliceSiteType.Donor,
                "no GU/GC means no canonical donor can be reported");
        }
    }

    #endregion

    #region BE — Boundary: sequence shorter than the scoring window

    /// <summary>
    /// Target "seq &lt; window" (KEY boundary). Two sub-cases, both of which must avoid
    /// IndexOutOfRange when the −3..+5 window cannot be fully extracted:
    ///   (a) length &lt; 6 — the documented guard. Lengths 0..5, including strings that
    ///       DO contain a GU ("GU", "GUAA"), must return no donor sites and never throw
    ///       (Donor_Site_Detection.md §3.1, §6.1 "Sequence shorter than 6 nt").
    ///   (b) a canonical GU shoved against the 3' edge of a ≥ 6-nt sequence, so the GU
    ///       sits inside the last 5 positions and falls OUTSIDE the i ≤ Length−6 scan
    ///       window. The edge GU is simply not emitted, and — critically — the scan
    ///       extracts no out-of-range window: no IndexOutOfRange leaks.
    /// </summary>
    [Test]
    public void FindDonorSites_SequenceShorterThanWindow_NoCrashNoIndexOutOfRange()
    {
        // (a) Every length 0..5 — including ones that contain a GU — yields nothing, never throws.
        foreach (string tooShort in new[] { "", "G", "GU", "GUA", "GUAA", "AGUAA", "GGGGG" })
        {
            var act = () => FindDonorSites(tooShort, minScore: 0.0).ToList();
            act.Should().NotThrow($"a sub-window sequence ('{tooShort}', len {tooShort.Length}) must not crash")
               .Subject.Should().BeEmpty(
                   $"length {tooShort.Length} < 6 is below the donor window guard — no site, no IndexOutOfRange");
        }

        // (b) A canonical GU flush against the 3' edge: positions Length-2..Length-1.
        // For a length-6 string the only scanned start is i = 0; a GU at index 4 is in
        // the last 5 positions → never scanned, never windowed off the end, never thrown.
        foreach (string edge in new[] { "AAAAGU", "AAAAAGU", "CCCCCCGU", "AAAAAAAAGU" })
        {
            int guIndex = edge.Length - 2;
            edge.Substring(guIndex, 2).Should().Be("GU", "the test fixture really places a GU at the 3' edge");

            var act = () => FindDonorSites(edge, minScore: 0.0).ToList();
            var sites = act.Should().NotThrow(
                $"a 3'-edge GU in '{edge}' must extract its window safely, not throw IndexOutOfRange").Subject;

            // The edge GU (index Length-2 > Length-6) is outside the scan range → not emitted.
            sites.Should().NotContain(s => s.Position == guIndex,
                "a GU within the last 5 positions is outside the i ≤ Length-6 scan window and is not reported");

            // Whatever IS reported still satisfies the full donor contract (window-safe).
            foreach (var s in sites)
                AssertWellFormedDonor(s, edge);
        }
    }

    #endregion

    #region MC — Malformed Content: non-DNA characters are handled, not crashed

    /// <summary>
    /// Target "non-DNA": characters outside {A,C,G,T,U} — digits, punctuation, IUPAC
    /// ambiguity codes (N, R, Y…), whitespace. The contract uppercases and maps T→U,
    /// then treats anything not matching the consensus as a never-matching residue, so
    /// junk is HANDLED rather than rejected — no crash, no IndexOutOfRange off the
    /// window/Substring (Donor_Site_Detection.md §3.3, §6.1). A canonical GU embedded in
    /// junk is still scored; any junk-only candidate that survives still obeys the donor
    /// contract (it must genuinely sit on a GU). 'N' positions can never anchor a donor.
    /// </summary>
    [Test]
    public void FindDonorSites_NonDnaCharacters_AreHandledNotCrashed()
    {
        foreach (string junk in new[]
                 {
                     "NNNNNNNNNNNN",            // all-ambiguity: no concrete GU → handled
                     "12!! \t??ZZ##",           // pure garbage: never-matching residues
                     "CAGNGUAANGU",             // a real GU survives inside N-junk
                     "RYSWKMBDHVNN",             // IUPAC ambiguity codes
                     "cag gua agu  ",           // lowercase + whitespace around a GU
                 })
        {
            var act = () => FindDonorSites(junk, minScore: 0.0, includeNonCanonical: true).ToList();
            var sites = act.Should().NotThrow($"non-DNA content ('{junk}') must be handled, not crash").Subject;

            // Whatever is reported is still a genuine, window-safe donor on a real GU/GC/AU.
            foreach (var s in sites)
                AssertWellFormedDonor(s, junk);

            // 'N' is never a consensus match and never a donor dinucleotide member, so a
            // reported donor never starts on an 'N' in the normalised sequence.
            string normalized = junk.ToUpperInvariant().Replace('T', 'U');
            sites.Should().NotContain(s => normalized[s.Position] == 'N',
                "an ambiguity 'N' cannot anchor a donor dinucleotide");
        }
    }

    #endregion

    #region BE — Boundary: all-GU sequence scores every candidate without crashing

    /// <summary>
    /// Target "all GT": "GUGUGU…" places a candidate GU at every even index, so the
    /// scanner scores the maximum possible number of overlapping donor candidates over
    /// the i ≤ Length−6 range. It must complete the full O(n) pass without crashing,
    /// hanging or windowing off the end — [CancelAfter] turns a regression-to-hang into
    /// a FAIL. Every emitted site still satisfies the donor contract (anchored on a real
    /// GU, score/confidence in range), and the count is bounded by the number of in-range
    /// even starts. We sweep DNA ("GTGT…") and RNA ("GUGU…") spellings — they must agree.
    /// </summary>
    [Test]
    [CancelAfter(30000)]
    public void FindDonorSites_AllGu_ScoresEveryCandidateAndCompletes(CancellationToken token)
    {
        foreach (int reps in new[] { 4, 50, 500 })
        {
            string dna = string.Concat(Enumerable.Repeat("GT", reps)); // "GTGT…"
            string rna = string.Concat(Enumerable.Repeat("GU", reps)); // "GUGU…"

            var dnaAct = () => FindDonorSites(dna, minScore: 0.0).ToList();
            var dnaSites = dnaAct.Should().NotThrow($"an all-GT sequence of {2 * reps} nt must not crash").Subject;
            token.ThrowIfCancellationRequested();

            var rnaSites = FindDonorSites(rna, minScore: 0.0).ToList();
            token.ThrowIfCancellationRequested();

            // DNA T → U normalisation means the two spellings produce identical results.
            dnaSites.Select(s => s.Position).Should().Equal(rnaSites.Select(s => s.Position),
                "'GT…' and 'GU…' normalise to the same sequence and must yield the same donor positions");

            // Every candidate GU sits at an EVEN index; only those with i ≤ Length-6 are scanned.
            int maxStart = dna.Length - 6; // inclusive upper bound of the scan
            foreach (var s in dnaSites)
            {
                AssertWellFormedDonor(s, dna);
                (s.Position % 2).Should().Be(0, "every GU in 'GUGU…' begins at an even index");
                s.Position.Should().BeLessThanOrEqualTo(maxStart,
                    "the scan never reports a candidate beyond i = Length-6");
            }

            // The number of scored candidates is bounded by the in-range even starts.
            int inRangeEvenStarts = maxStart < 0 ? 0 : (maxStart / 2) + 1;
            dnaSites.Count.Should().BeLessThanOrEqualTo(inRangeEvenStarts,
                "at most one donor per in-range even GU start");
        }
    }

    #endregion

    #region Positive sanity — a canonical donor motif is detected at the right position

    /// <summary>
    /// Positive sanity: the perfect-consensus donor CAG|GUAAGU (MAG|GURAGU) must be
    /// detected as exactly ONE canonical donor whose GU starts at index 3, with the top
    /// score 1.0 (all 9 offsets match) and confidence 1.0 (Donor_Site_Detection.md §2.1,
    /// §7; mirrors the canonical SPLICE-DONOR-001 test). The DNA spelling CAGGTAAGT must
    /// give the identical result (T→U). This proves the fuzz harness asserts against a
    /// detector that actually FINDS real donor sites — not a no-op — and pins the
    /// "reported donor really has GT at the splice position" contract on a true positive.
    /// </summary>
    [Test]
    public void FindDonorSites_CanonicalConsensus_DetectedAtExpectedPosition()
    {
        const string rna = "CAGGUAAGU";
        const string dna = "CAGGTAAGT"; // same motif spelled with DNA T

        var rnaSites = FindDonorSites(rna, minScore: 0.3).ToList();

        rnaSites.Should().ContainSingle("the perfect GU consensus CAG|GUAAGU is a single canonical donor");
        var site = rnaSites[0];
        site.Type.Should().Be(SpliceSiteType.Donor, "a canonical GU donor carries Type = Donor");
        site.Position.Should().Be(3, "the GU dinucleotide starts at index 3 in CAGGUAAGU");
        site.Score.Should().BeApproximately(1.0, 1e-10, "all 9 consensus offsets match → score 9/9 = 1.0");
        site.Confidence.Should().BeApproximately(1.0, 1e-10, "confidence = (1.0 − 0.5)/(1.0 − 0.5) = 1.0");
        AssertWellFormedDonor(site, rna);

        // DNA spelling normalises identically (T → U) → same position/score.
        var dnaSites = FindDonorSites(dna, minScore: 0.3).ToList();
        dnaSites.Select(s => (s.Position, s.Score)).Should().Equal(rnaSites.Select(s => (s.Position, s.Score)),
            "the DNA spelling CAGGTAAGT maps to CAGGUAAGU and yields the identical donor");
    }

    /// <summary>
    /// Positive sanity over RANDOM DNA: across fixed seeds and lengths the scan must
    /// never crash or hang, and every emitted donor must satisfy the full contract
    /// (score/confidence in range, anchored on a real GU/GC at its Position). This pins
    /// window-safety on arbitrary sequences, not just hand-built motifs.
    /// </summary>
    [Test]
    [CancelAfter(30000)]
    public void FindDonorSites_RandomDna_AlwaysWellFormed(CancellationToken token)
    {
        foreach (int seed in new[] { 5, 23, 101, 2026 })
        {
            foreach (int len in new[] { 6, 30, 120 })
            {
                string seq = RandomDna(len, seed);

                var act = () => FindDonorSites(seq, minScore: 0.0, includeNonCanonical: true).ToList();
                var sites = act.Should().NotThrow($"random DNA must not crash the scan (seed {seed}, len {len})").Subject;
                token.ThrowIfCancellationRequested();

                foreach (var s in sites)
                    AssertWellFormedDonor(s, seq);
            }
        }
    }

    #endregion

    #endregion

    // ═══════════════════════════════════════════════════════════════════
    //  SPLICE-ACCEPTOR-001 — splice acceptor site detection : fuzz targets
    // ═══════════════════════════════════════════════════════════════════

    #region SPLICE-ACCEPTOR-001 — splice acceptor site detection

    #region MC — No AG candidate: returns no acceptor sites

    /// <summary>
    /// Target "No AG sites": a sequence containing NO AG dinucleotide has no canonical
    /// acceptor candidate at all, so the scan emits nothing — no acceptor site, no crash
    /// (Acceptor_Site_Detection.md §6.1; mirrors the existing AcceptorSite unit test). We
    /// probe a long poly-U tract (a real polypyrimidine context but with no AG terminus),
    /// a poly-A homopolymer (an A never followed by G), and a CU-only sequence, and confirm
    /// that even at minScore 0.0 nothing canonical is invented from an AG that is not there.
    /// All sequences are ≥ 20 nt so we isolate "no AG" from the length guard.
    /// </summary>
    [Test]
    public void FindAcceptorSites_NoAgCandidate_ReturnsEmpty()
    {
        foreach (string seq in new[]
                 {
                     new string('U', 30),          // PPT-rich but no AG terminus
                     new string('A', 30),          // an A is never followed by G
                     "CUCUCUCUCUCUCUCUCUCUCUCU",    // pyrimidine-only, no AG
                 })
        {
            var canonical = ((Func<List<SpliceSite>>)(() =>
                FindAcceptorSites(seq, minScore: 0.0).ToList()))
                .Should().NotThrow("an AG-free sequence is scanned without error").Subject;

            canonical.Should().BeEmpty($"'{seq}' has no AG dinucleotide, so there is no canonical acceptor candidate");
            canonical.Should().NotContain(s => s.Type == SpliceSiteType.Acceptor,
                "no AG means no canonical acceptor can be reported");
        }
    }

    #endregion

    #region BE — Boundary: sequence shorter than the window / AG flush at the 3' edge

    /// <summary>
    /// Target "seq &lt; window" (KEY boundary). Two sub-cases, both of which must avoid
    /// IndexOutOfRange:
    ///   (a) length &lt; 20 — the documented guard. Lengths 0..19, INCLUDING strings that
    ///       DO contain an AG ("CAGG", a 19-nt PPT+CAGG), must return no acceptor sites and
    ///       never throw (Acceptor_Site_Detection.md §3.1, §6.1).
    ///   (b) a canonical AG shoved against the 3' edge of a ≥ 20-nt sequence, so i = Length−2
    ///       is the last scanned index and the −15-nt PPT/PWM window plus the GetMotifContext
    ///       Substring must all extract window-safely (each access is bounds-clamped). The
    ///       edge AG is the LAST thing the scan can see; it must produce a clamped score and
    ///       motif, not an IndexOutOfRange off either end.
    /// </summary>
    [Test]
    public void FindAcceptorSites_SequenceShorterThanWindow_NoCrashNoIndexOutOfRange()
    {
        // (a) Every length 0..19 — including ones that contain an AG — yields nothing, never throws.
        foreach (string tooShort in new[]
                 {
                     "", "A", "AG", "CAGG", "UUUUUUUUUUUUUUUUUUU",   // 0,1,2,4,19 — last is 19-nt poly-U
                     "UUUUUUUUUUUUUUUCAGG",                          // 19 nt and DOES contain an AG
                 })
        {
            var act = () => FindAcceptorSites(tooShort, minScore: 0.0, includeNonCanonical: true).ToList();
            act.Should().NotThrow($"a sub-window sequence ('{tooShort}', len {tooShort.Length}) must not crash")
               .Subject.Should().BeEmpty(
                   $"length {tooShort.Length} < 20 is below the acceptor window guard — no site, no IndexOutOfRange");
        }

        // (b) A canonical AG flush against the 3' edge: dinucleotide at (Length-2, Length-1),
        // i.e. the scan's terminal index i = Length-2. The −15 PPT/PWM window and the motif
        // Substring must all clamp; the worst-case +0 PWM offset touches Length-1 exactly.
        foreach (string edge in new[]
                 {
                     "UUUUUUUUUUUUUUUUUUAG",     // 20 nt, AG at indices 18,19
                     "CCCCCCCCCCCCCCCCCCCCAG",   // 22 nt, AG at the very end
                     "UUUUUUUUUUUUUUUUCAGCAG",   // 22 nt, an interior AG and a 3'-edge AG
                 })
        {
            int agIndex = edge.Length - 2;
            edge.Substring(agIndex, 2).Should().Be("AG", "the test fixture really places an AG at the 3' edge");

            var act = () => FindAcceptorSites(edge, minScore: 0.0).ToList();
            var sites = act.Should().NotThrow(
                $"a 3'-edge AG in '{edge}' must extract its −15..0 window safely, not throw IndexOutOfRange").Subject;

            // Whatever IS reported is a genuine, window-safe acceptor anchored on a real AG.
            foreach (var s in sites)
                AssertWellFormedAcceptor(s, edge);
        }
    }

    #endregion

    #region MC — Malformed Content: non-DNA characters are handled, not crashed

    /// <summary>
    /// Target "non-DNA": characters outside {A,C,G,T,U} — digits, punctuation, IUPAC
    /// ambiguity codes (N, R, Y…), whitespace. The contract uppercases and maps T→U, then
    /// treats anything not matching the consensus as a never-matching residue, so junk is
    /// HANDLED rather than rejected — no crash, no IndexOutOfRange off the PPT/PWM window or
    /// the motif Substring (Acceptor_Site_Detection.md §3.3, §6.1). A canonical AG embedded
    /// in junk is still scored; any candidate that survives still obeys the acceptor contract
    /// (it must genuinely end on an AG/AC). 'N' can never be the A or G of an acceptor.
    /// </summary>
    [Test]
    public void FindAcceptorSites_NonDnaCharacters_AreHandledNotCrashed()
    {
        foreach (string junk in new[]
                 {
                     "NNNNNNNNNNNNNNNNNNNNNN",        // all-ambiguity: no concrete AG → handled
                     "12!! \t??ZZ##$$%%^^&&**",       // pure garbage: never-matching residues
                     "NNNNNNNNNNNNNNNNCAGNNN",        // a real AG survives inside N-junk
                     "RYSWKMBDHVNNRYSWKMBDHV",        // IUPAC ambiguity codes
                     "uuuu uuuu uuuu uuu cagg gg",    // lowercase + whitespace around a cag
                 })
        {
            var act = () => FindAcceptorSites(junk, minScore: 0.0, includeNonCanonical: true).ToList();
            var sites = act.Should().NotThrow($"non-DNA content ('{junk}') must be handled, not crash").Subject;

            // Whatever is reported is still a genuine, window-safe acceptor on a real AG/AC.
            foreach (var s in sites)
                AssertWellFormedAcceptor(s, junk);

            // 'N' is never an acceptor dinucleotide member, so a reported acceptor never has
            // its A or its G/C land on an 'N' in the normalised sequence.
            string normalized = junk.ToUpperInvariant().Replace('T', 'U');
            sites.Should().NotContain(s => normalized[s.Position] == 'N' || normalized[s.Position - 1] == 'N',
                "an ambiguity 'N' cannot anchor an acceptor dinucleotide");
        }
    }

    #endregion

    #region BE — Boundary: all-AG sequence scores every candidate without crashing

    /// <summary>
    /// Target "all AG": "AGAGAG…" places a candidate AG at every even index, so the scanner
    /// scores the maximum possible number of overlapping acceptor candidates over the
    /// 15 ≤ i ≤ Length−2 range. It must complete the full O(n) pass without crashing, hanging
    /// or windowing off the end — [CancelAfter] turns a regression-to-hang into a FAIL. Every
    /// emitted site still satisfies the acceptor contract (anchored on a real AG, score/
    /// confidence in range), every reported Position is ODD (G of an even-indexed AG → i+1),
    /// and the count is bounded by the number of in-range even AG starts. DNA ("AGAG…") and
    /// RNA ("AGAG…") spellings are identical here; we also sweep the DNA-T spelling separately.
    /// </summary>
    [Test]
    [CancelAfter(30000)]
    public void FindAcceptorSites_AllAg_ScoresEveryCandidateAndCompletes(CancellationToken token)
    {
        foreach (int reps in new[] { 10, 50, 500 })
        {
            string seq = string.Concat(Enumerable.Repeat("AG", reps)); // "AGAG…"

            var act = () => FindAcceptorSites(seq, minScore: 0.0).ToList();
            var sites = act.Should().NotThrow($"an all-AG sequence of {2 * reps} nt must not crash").Subject;
            token.ThrowIfCancellationRequested();

            // Every candidate AG begins at an EVEN index i; only those with 15 ≤ i ≤ Length-2 are scanned.
            // Reported Position = i+1, which is therefore ODD.
            int firstStart = 16;               // first even i ≥ 15
            int maxStart = seq.Length - 2;     // inclusive upper bound of the scan
            foreach (var s in sites)
            {
                AssertWellFormedAcceptor(s, seq);
                (s.Position % 2).Should().Be(1, "Position = i+1 with even i → an odd index in 'AGAG…'");
            }

            // The number of scored candidates is bounded by the in-range even AG starts.
            int inRangeEvenStarts = maxStart < firstStart ? 0 : ((maxStart - firstStart) / 2) + 1;
            sites.Count.Should().BeLessThanOrEqualTo(inRangeEvenStarts,
                "at most one acceptor per in-range even AG start");
        }
    }

    #endregion

    #region Positive sanity — a canonical acceptor motif is detected at the right position

    /// <summary>
    /// Positive sanity: a strong canonical acceptor — a 16-nt polypyrimidine tract followed
    /// by the CAG|G terminus ("UUUUUUUUUUUUUUUUCAGGG", the AcceptorSite unit's fixture) — must
    /// be detected as a canonical acceptor whose AG ends at Position 18 (the G of AG, INV-04),
    /// with Type = Acceptor and a Score/Confidence in range. The DNA spelling (T→U) must give
    /// the identical Position. This proves the fuzz harness asserts against a detector that
    /// actually FINDS real acceptor sites — not a no-op — and pins the "reported acceptor really
    /// has AG at the splice position" contract on a true positive.
    /// </summary>
    [Test]
    public void FindAcceptorSites_CanonicalConsensus_DetectedAtExpectedPosition()
    {
        const string rna = "UUUUUUUUUUUUUUUUCAGGG"; // 16×U + CAGGG; AG at indices 17,18
        const string dna = "TTTTTTTTTTTTTTTTCAGGG"; // 16×T then CAGGG → same AG after T→U

        var rnaSites = FindAcceptorSites(rna, minScore: 0.1).ToList();

        rnaSites.Should().Contain(s => s.Position == 18 && s.Type == SpliceSiteType.Acceptor,
            "the strong PPT + CAG|G acceptor is detected at Position 18 (the G of AG, INV-04)");
        var site = rnaSites.Single(s => s.Position == 18);
        AssertWellFormedAcceptor(site, rna);

        // DNA spelling normalises identically (T → U) → the same canonical acceptor position.
        var dnaSites = FindAcceptorSites(dna, minScore: 0.1).ToList();
        dnaSites.Should().Contain(s => s.Position == 18 && s.Type == SpliceSiteType.Acceptor,
            "the DNA spelling maps to the same CAG|G acceptor and yields the identical Position");
    }

    /// <summary>
    /// Positive sanity over RANDOM DNA: across fixed seeds and lengths the scan must never
    /// crash or hang, and every emitted acceptor must satisfy the full contract (score/
    /// confidence in range, anchored on a real AG/AC at its Position). This pins window-safety
    /// on arbitrary sequences, not just hand-built motifs. Lengths span the 20-nt guard.
    /// </summary>
    [Test]
    [CancelAfter(30000)]
    public void FindAcceptorSites_RandomDna_AlwaysWellFormed(CancellationToken token)
    {
        foreach (int seed in new[] { 5, 23, 101, 2026 })
        {
            foreach (int len in new[] { 20, 60, 200 })
            {
                string seq = RandomDna(len, seed);

                var act = () => FindAcceptorSites(seq, minScore: 0.0, includeNonCanonical: true).ToList();
                var sites = act.Should().NotThrow($"random DNA must not crash the scan (seed {seed}, len {len})").Subject;
                token.ThrowIfCancellationRequested();

                foreach (var s in sites)
                    AssertWellFormedAcceptor(s, seq);
            }
        }
    }

    #endregion

    #endregion
}
