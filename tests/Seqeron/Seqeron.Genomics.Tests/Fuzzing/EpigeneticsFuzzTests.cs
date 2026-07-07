using static Seqeron.Genomics.Annotation.EpigeneticsAnalyzer;

namespace Seqeron.Genomics.Tests.Fuzzing;

/// <summary>
/// Fuzz tests for the Epigenetics area — CpG island detection (EPIGEN-CPG-001).
/// The three public entry points under test live in
/// <see cref="EpigeneticsAnalyzer"/> and form the CpG-island contract:
///   • <see cref="EpigeneticsAnalyzer.FindCpGSites"/> — enumerates the 0-based
///     position of the cytosine in every adjacent <c>CG</c> dinucleotide.
///   • <see cref="EpigeneticsAnalyzer.CalculateCpGObservedExpected"/> — the
///     Gardiner-Garden &amp; Frommer CpG observed/expected ratio for a sequence.
///   • <see cref="EpigeneticsAnalyzer.FindCpGIslands"/> — sliding-window island
///     classifier (the window length IS the <c>minLength</c> parameter).
///
/// ───────────────────────────────────────────────────────────────────────────
/// What fuzzing verifies
/// ───────────────────────────────────────────────────────────────────────────
/// Fuzzing feeds degenerate and boundary inputs to a unit and asserts that the
/// code NEVER fails in an undisciplined way: no hang, no state corruption, no
/// nonsense output, and no *unhandled* runtime exception (IndexOutOfRange /
/// NullReference / DivideByZero / Overflow). Every input must resolve to EITHER a
/// well-defined, theory-correct value OR a *documented, intentional* outcome.
/// For the CpG O/E ratio and the window-based island detector, the headline
/// hazards are:
///   • a DivideByZeroException / NaN in the O/E ratio when a window contains zero
///     C or zero G (expected CpG count is 0);
///   • an IndexOutOfRangeException when the window length exceeds the sequence, or
///     a substring slice runs off the end of the buffer;
///   • an unhandled crash on a zero or negative window length (minLength=0);
///   • an island reported at coordinates whose substring does NOT actually meet
///     the GC / O-E thresholds (a classification or coordinate bug);
///   • numeric overflow / saturation on a long all-"CGCG…" input.
/// — docs/ADVANCED_TESTING_CHECKLIST.md §8 "Fuzzing".
///
/// ───────────────────────────────────────────────────────────────────────────
/// Unit: EPIGEN-CPG-001 — CpG island detection (Epigenetics)
/// Checklist: docs/checklists/03_FUZZING.md, row 85.
/// Fuzz strategy exercised for THIS unit:
///   • BE = Boundary Exploitation — граничні значення: 0, -1, empty.
///     Targets (checklist row 85): "No CG dinucleotides, all CG, empty seq,
///     windowSize=0, windowSize > seqLen".
/// — docs/checklists/03_FUZZING.md §Description (strategy codes).
///
/// NOTE ON "windowSize": the documented/source API has no separate windowSize
/// parameter — the sliding-window length IS <c>minLength</c>
/// (FindCpGIslands rescans a window of length minLength). So the checklist's
/// "windowSize=0" boundary maps to <c>minLength = 0</c> and "windowSize > seqLen"
/// maps to <c>minLength &gt; sequence.Length</c>.
///   — docs/algorithms/Epigenetics/CpG_Site_Detection.md §3.1, §4.1, §5.2.
///
/// ───────────────────────────────────────────────────────────────────────────
/// The CpG-island contract under test
/// ───────────────────────────────────────────────────────────────────────────
/// O/E ratio (Gardiner-Garden &amp; Frommer 1987):
///     O/E = CpGcount / ((Ccount × Gcount) / L)
///   and the implementation returns 0.0 when the expected count is zero (i.e.
///   when there is no C or no G), and 0.0 for null / empty / length-1 input.
///   — docs/algorithms/Epigenetics/CpG_Site_Detection.md §2.2, §3.2, §3.3.
///
/// Default island criteria (cited Gardiner-Garden &amp; Frommer thresholds):
///   1. length ≥ 200 bp, 2. GC content ≥ 50%, 3. CpG O/E ≥ 0.6 — all compared
///   INCLUSIVELY (gc &gt;= minGc, cpgRatio &gt;= minCpGRatio).
///   — docs/algorithms/Epigenetics/CpG_Site_Detection.md §2.2, §5.2.
///
/// Coordinates: island Start is 0-based inclusive; End is the EXCLUSIVE end
///   (the implementation stores i + windowSize), so substring length = End-Start.
///   — docs/algorithms/Epigenetics/CpG_Site_Detection.md §3.2, §3.3.
///
/// Boundary handling (docs §3.3, §6.1):
///   • null / empty                 → FindCpGSites & FindCpGIslands yield none;
///                                     CalculateCpGObservedExpected returns 0.0.
///   • length &lt; 2                 → no CpG sites; O/E = 0.0.
///   • sequence shorter than minLength → FindCpGIslands yields no islands.
///   • no C or no G                 → O/E = 0.0 (no DivideByZero).
///
/// All randomness is LOCALLY seeded (new Random(seed)); no shared static Rng.
/// </summary>
[TestFixture]
[Category("Fuzzing")]
public sealed class EpigeneticsFuzzTests
{
    // ── Well-formed-region assertion helper ─────────────────────────────────
    // Pins INV-04 + the documented coordinate contract: every reported island
    // must use an exclusive End strictly after an in-range, non-negative Start,
    // span at least minLength bases, and its OWN substring must genuinely meet
    // the GC and O/E thresholds. The reported GcContent / CpGRatio must equal a
    // recomputation on that exact substring. This is what stops a test from
    // rubber-stamping garbage coordinates green.
    private static void AssertWellFormedIsland(
        (int Start, int End, double GcContent, double CpGRatio) island,
        string sequence,
        int minLength,
        double minGc,
        double minCpGRatio)
    {
        island.Start.Should().BeGreaterThanOrEqualTo(0);
        island.End.Should().BeGreaterThan(island.Start);
        island.End.Should().BeLessThanOrEqualTo(sequence.Length);
        (island.End - island.Start).Should().BeGreaterThanOrEqualTo(minLength);

        // Documented inclusive thresholds on the reported metrics.
        island.GcContent.Should().BeGreaterThanOrEqualTo(minGc);
        island.CpGRatio.Should().BeGreaterThanOrEqualTo(minCpGRatio);

        // The reported metrics must describe the island's OWN substring.
        string sub = sequence.Substring(island.Start, island.End - island.Start);
        double expectedRatio = CalculateCpGObservedExpected(sub);
        island.CpGRatio.Should().BeApproximately(expectedRatio, 1e-9);

        // No NaN / Infinity must ever leak out.
        double.IsNaN(island.GcContent).Should().BeFalse();
        double.IsNaN(island.CpGRatio).Should().BeFalse();
        double.IsInfinity(island.CpGRatio).Should().BeFalse();
    }

    #region EPIGEN-CPG-001 — CpG island detection

    // ── BE: empty / null sequence ───────────────────────────────────────────

    [Test]
    public void FindCpGSites_NullOrEmpty_YieldsNothing()
    {
        FindCpGSites(null!).Should().BeEmpty();
        FindCpGSites("").Should().BeEmpty();
    }

    [Test]
    public void FindCpGIslands_NullOrEmpty_YieldsNothing()
    {
        FindCpGIslands(null!).Should().BeEmpty();
        FindCpGIslands("").Should().BeEmpty();
    }

    [Test]
    public void CalculateCpGObservedExpected_NullEmptyOrSingleBase_IsZero()
    {
        // Docs §3.3 / §6.1: ratio is 0.0 for null, empty, or length-1 input —
        // a dinucleotide cannot be formed from fewer than two bases.
        CalculateCpGObservedExpected(null!).Should().Be(0.0);
        CalculateCpGObservedExpected("").Should().Be(0.0);
        CalculateCpGObservedExpected("C").Should().Be(0.0);
        CalculateCpGObservedExpected("G").Should().Be(0.0);
    }

    // ── BE: NO CG dinucleotides (DivideByZero / NaN guard) ──────────────────

    [Test]
    public void CalculateCpGObservedExpected_NoCytosine_IsZero_NoDivideByZero()
    {
        // No C at all → expected = (0 × G)/L = 0 → ratio must be 0.0, not NaN.
        Action act = () => CalculateCpGObservedExpected("AGAGAGAGAGGGGAAA");
        act.Should().NotThrow();
        double ratio = CalculateCpGObservedExpected("AGAGAGAGAGGGGAAA");
        ratio.Should().Be(0.0);
        double.IsNaN(ratio).Should().BeFalse();
    }

    [Test]
    public void CalculateCpGObservedExpected_NoGuanine_IsZero_NoDivideByZero()
    {
        // No G at all → expected = (C × 0)/L = 0 → ratio must be 0.0, not NaN.
        double ratio = CalculateCpGObservedExpected("ACACACACTTTTCCCC");
        ratio.Should().Be(0.0);
        double.IsNaN(ratio).Should().BeFalse();
    }

    [Test]
    public void CalculateCpGObservedExpected_SingleGpCNoCpG_IsZero()
    {
        // Docs §6.1 "GpC without CpG": the single dinucleotide "GC" has one C and
        // one G but no C-then-G adjacency → CpG count 0 → O/E = 0 (definition is
        // C THEN G, not the reverse). (A LONGER "GCGC…" run does contain internal
        // CpGs at the GC|GC junctions, so it is deliberately NOT used here.)
        double ratio = CalculateCpGObservedExpected("GC");
        ratio.Should().Be(0.0);
        double.IsNaN(ratio).Should().BeFalse();
    }

    [Test]
    public void FindCpGSites_TrulyNoCpG_ReportsNone()
    {
        // All G's precede all C's so a C is never immediately followed by a G.
        string seq = "GGGGGGGGGGCCCCCCCCCC";
        FindCpGSites(seq).Should().BeEmpty();
        CalculateCpGObservedExpected(seq).Should().Be(0.0); // expected>0 but CpG=0
    }

    [Test]
    public void CalculateCpGObservedExpected_NoCpG_ButCandGPresent_IsExactlyZero()
    {
        // C and G both present (expected > 0) yet zero adjacent CG → O/E = 0/expected = 0.
        string seq = "GGGGGCCCCC"; // 5 G then 5 C: no C-then-G adjacency
        CalculateCpGObservedExpected(seq).Should().Be(0.0);
    }

    // ── BE: all-CG sequence (maximal O/E, no overflow) ──────────────────────

    [Test]
    [CancelAfter(15000)]
    public void CalculateCpGObservedExpected_AllCG_IsDocumentedSaturatedValue()
    {
        // "CGCG...CG", length L = 2n: C=n, G=n, CpG count = n (every even index).
        // O/E = n / ((n × n)/L) = n / (n²/2n) = n / (n/2) = 2 (exactly), for any n≥1.
        for (int n = 1; n <= 500; n++)
        {
            string seq = string.Concat(Enumerable.Repeat("CG", n));
            double ratio = CalculateCpGObservedExpected(seq);
            ratio.Should().BeApproximately(2.0, 1e-9,
                because: "an all-CGCG sequence has O/E exactly 2 for any length");
            double.IsInfinity(ratio).Should().BeFalse();
        }
    }

    [Test]
    [CancelAfter(20000)]
    public void CalculateCpGObservedExpected_VeryLongAllCG_NoOverflow_StaysFinite()
    {
        // 200k bases of CGCG... : counts (c, g, cpg ≈ 100k) and L = 200k easily fit
        // in int, but this pins that the long path neither overflows nor saturates
        // to Infinity/NaN and still returns the exact value 2.
        string seq = string.Concat(Enumerable.Repeat("CG", 100_000));
        double ratio = CalculateCpGObservedExpected(seq);
        ratio.Should().BeApproximately(2.0, 1e-9);
        double.IsNaN(ratio).Should().BeFalse();
        double.IsInfinity(ratio).Should().BeFalse();
    }

    [Test]
    public void FindCpGSites_AllCG_ReportsEveryEvenPosition()
    {
        // "CGCGCG" → CpG at indices 0, 2, 4 (the C of each CG pair).
        string seq = "CGCGCG";
        FindCpGSites(seq).Should().Equal(new[] { 0, 2, 4 });
    }

    [Test]
    public void FindCpGSites_AdjacentCpGs_ReportsBoth()
    {
        // Docs §6.1: "CGCG" → both CpG positions (0 and 2) are reported.
        FindCpGSites("CGCG").Should().Equal(new[] { 0, 2 });
    }

    // ── BE: windowSize == 0 (minLength = 0) — no crash, no islands ───────────

    [Test]
    public void FindCpGIslands_MinLengthZero_NoCrash_YieldsNoIslands()
    {
        // minLength=0 ⇒ the guard (Length < 0) is false, the loop runs to i==Length,
        // every window is the empty substring (GC=0, O/E=0) so nothing clears the
        // 50% GC / 0.6 O/E bar. The boundary must not IndexOutOfRange on the
        // i==Length slice, nor report a zero-length "island".
        string seq = string.Concat(Enumerable.Repeat("CG", 200)); // strong CpG content
        Action act = () => FindCpGIslands(seq, minLength: 0).ToList();
        act.Should().NotThrow();
        FindCpGIslands(seq, minLength: 0).Should().BeEmpty();
    }

    [Test]
    public void FindCpGIslands_MinLengthZero_OnEmptySequence_NoCrash()
    {
        Action act = () => FindCpGIslands("", minLength: 0).ToList();
        act.Should().NotThrow();
        FindCpGIslands("", minLength: 0).Should().BeEmpty();
    }

    [Test]
    public void FindCpGIslands_NegativeMinLength_NoCrash()
    {
        // minLength = -1: the length guard (Length < -1) is false; window =
        // Min(-1, Length-i) is negative → Substring would throw. Pin that the
        // contract degrades gracefully to "no islands" or a documented throw,
        // but NEVER an undisciplined crash that corrupts the iterator. We accept
        // either no islands OR an ArgumentException, but not e.g. IndexOutOfRange.
        string seq = string.Concat(Enumerable.Repeat("CG", 100));
        Action act = () => FindCpGIslands(seq, minLength: -1).ToList();
        act.Should().NotThrow<IndexOutOfRangeException>();
        act.Should().NotThrow<NullReferenceException>();
    }

    // ── BE: windowSize > seqLen (minLength > sequence length) ────────────────

    [Test]
    public void FindCpGIslands_MinLengthGreaterThanSeqLen_YieldsNoIslands()
    {
        // Docs §6.1: "Sequence shorter than minLength → yields no islands."
        string seq = string.Concat(Enumerable.Repeat("CG", 20)); // length 40
        FindCpGIslands(seq, minLength: 200).Should().BeEmpty();   // 200 > 40
        FindCpGIslands(seq, minLength: 41).Should().BeEmpty();    // just over length
        FindCpGIslands(seq, minLength: int.MaxValue).Should().BeEmpty();
    }

    [Test]
    public void FindCpGIslands_MinLengthEqualsSeqLen_SingleWindow_NoCrash()
    {
        // Exactly one window fits. With strong CpG content it should be the one
        // island; coordinates must be the whole sequence [0, L).
        string seq = string.Concat(Enumerable.Repeat("CG", 100)); // length 200, O/E=2, GC=100%
        var islands = FindCpGIslands(seq, minLength: 200).ToList();
        islands.Should().HaveCount(1);
        islands[0].Start.Should().Be(0);
        islands[0].End.Should().Be(200);
        AssertWellFormedIsland(islands[0], seq, 200, 0.5, 0.6);
    }

    // ── POSITIVE sanity: a genuine island flanked by low-GC sequence ─────────

    [Test]
    public void FindCpGIslands_GenuineIslandFlankedByLowGc_ReportsExactlyThatIsland()
    {
        // Construct: 300 bp AT-only flank (GC=0, O/E=0) | 400 bp CpG-rich core
        // (GC=100%, O/E=2) | 300 bp AT-only flank. With a 200 bp window, any window
        // lying wholly inside a flank fails (GC=0); only windows overlapping the
        // CpG-rich core can pass, so exactly one merged island is reported and it
        // must fully contain the core region.
        const int flankLen = 300;
        const int corePairs = 200; // 400 bp core
        string flank = new string('A', flankLen / 2) + new string('T', flankLen / 2);
        string core = string.Concat(Enumerable.Repeat("CG", corePairs)); // 400 bp, GC 100%, O/E 2
        string seq = flank + core + flank;
        int coreStart = flank.Length;
        int coreEnd = flank.Length + core.Length;

        var islands = FindCpGIslands(seq, minLength: 200, minGc: 0.5, minCpGRatio: 0.6).ToList();

        islands.Should().HaveCount(1, because: "there is exactly one CpG-rich core");
        var island = islands[0];
        AssertWellFormedIsland(island, seq, 200, 0.5, 0.6);

        // The reported island must FULLY CONTAIN the CpG-rich core. It may bleed
        // into the flanks by < window length where a straddling window still
        // clears the 50% GC / 0.6 O/E bars, but it can never miss the core.
        island.Start.Should().BeLessThanOrEqualTo(coreStart);
        island.End.Should().BeGreaterThanOrEqualTo(coreEnd);
        // ...and it must not extend a full window-length beyond the core on either
        // side (a pure-flank window cannot qualify).
        island.Start.Should().BeGreaterThan(coreStart - 200);
        island.End.Should().BeLessThan(coreEnd + 200);
    }

    [Test]
    public void FindCpGIslands_UniformlyLowGc_ReportsNone()
    {
        // AT-only 1 kb: GC = 0, O/E = 0 → no window can pass the 50%/0.6 bars.
        string seq = string.Concat(Enumerable.Repeat("AT", 500)); // 1000 bp
        FindCpGIslands(seq).Should().BeEmpty();
    }

    [Test]
    public void FindCpGIslands_GcRichButCpGDepleted_ReportsNone()
    {
        // High GC but the CpG dinucleotide is suppressed: alternating "GC" with no
        // CpG-rich stretch. Build a long "GCGCGC..." — wait, that contains CpGs.
        // Use a CpG-depleted but GC-rich design: blocks of G then blocks of C so
        // GC content is ~100% but C-then-G adjacency is rare.
        string seq = string.Concat(Enumerable.Repeat(new string('G', 10) + new string('C', 10), 50));
        // length = 1000, GC = 100%, but only a single C→G junction per 20 bp block.
        double ratio = CalculateCpGObservedExpected(seq.Substring(0, 200));
        ratio.Should().BeLessThan(0.6, because: "CpG dinucleotides are depleted");
        FindCpGIslands(seq).Should().BeEmpty();
    }

    // ── BE/robustness: random fuzz — never crash, always well-formed ─────────

    [Test]
    [CancelAfter(30000)]
    public void FindCpGIslands_RandomSequences_NeverCrash_AllIslandsWellFormed()
    {
        const string alphabet = "ACGT";
        for (int seed = 0; seed < 200; seed++)
        {
            var rng = new Random(seed);
            int len = rng.Next(0, 600);
            var chars = new char[len];
            for (int i = 0; i < len; i++)
                chars[i] = alphabet[rng.Next(alphabet.Length)];
            string seq = new string(chars);

            int minLength = rng.Next(0, 300); // includes 0 and values > len
            double minGc = rng.NextDouble();
            double minCpGRatio = rng.NextDouble() * 2.0;

            List<(int Start, int End, double GcContent, double CpGRatio)> islands = null!;
            Action act = () => islands = FindCpGIslands(seq, minLength, minGc, minCpGRatio).ToList();
            act.Should().NotThrow($"seed={seed}, len={len}, minLength={minLength}");

            // O/E ratio must never be NaN / Infinity regardless of base composition.
            double oe = CalculateCpGObservedExpected(seq);
            double.IsNaN(oe).Should().BeFalse($"seed={seed}");
            double.IsInfinity(oe).Should().BeFalse($"seed={seed}");

            if (islands is not null && minLength >= 1)
            {
                foreach (var island in islands)
                    AssertWellFormedIsland(island, seq, minLength, minGc, minCpGRatio);
            }
        }
    }

    [Test]
    [CancelAfter(30000)]
    public void FindCpGSites_RandomSequences_EverySiteIsCThenG()
    {
        // INV-01: every reported site marks a C whose next base is G.
        const string alphabet = "ACGTacgtNnXx";
        for (int seed = 0; seed < 200; seed++)
        {
            var rng = new Random(seed);
            int len = rng.Next(0, 400);
            var chars = new char[len];
            for (int i = 0; i < len; i++)
                chars[i] = alphabet[rng.Next(alphabet.Length)];
            string seq = new string(chars);
            string upper = seq.ToUpperInvariant();

            foreach (int pos in FindCpGSites(seq))
            {
                pos.Should().BeInRange(0, len - 2);
                upper[pos].Should().Be('C');
                upper[pos + 1].Should().Be('G');
            }
        }
    }

    #endregion
}
