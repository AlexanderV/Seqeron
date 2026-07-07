using static Seqeron.Genomics.Annotation.EpigeneticsAnalyzer;

namespace Seqeron.Genomics.Tests.Fuzzing;

/// <summary>
/// Fuzz tests for the Epigenetics area — per-CpG methylation level (β-value)
/// calling (EPIGEN-METHYL-001). The single public entry point under test lives in
/// <see cref="EpigeneticsAnalyzer"/>:
///   • <see cref="EpigeneticsAnalyzer.CalculateMethylationFromBisulfite(string, IEnumerable{ValueTuple{string, int}})"/>
///     — calls a per-CpG methylation level from aligned bisulfite reads against a
///     reference. At each reference CpG cytosine, a read base of <c>C</c> is a
///     methylated call and a read base of <c>T</c> is an unmethylated call; any
///     other read base is not a valid bisulfite call and is ignored. The reported
///     <c>MethylationLevel</c> is the Bismark β-fraction
///     <c>methylated / (methylated + unmethylated)</c> and <c>Coverage</c> is the
///     number of valid C/T calls at that site.
///
/// ───────────────────────────────────────────────────────────────────────────
/// What fuzzing verifies
/// ───────────────────────────────────────────────────────────────────────────
/// Fuzzing feeds degenerate and boundary inputs to a unit and asserts that the
/// code NEVER fails in an undisciplined way: no hang, no state corruption, no
/// nonsense output, and no *unhandled* runtime exception (DivideByZero / NaN /
/// IndexOutOfRange / NullReference). Every input must resolve to EITHER a
/// well-defined, theory-correct value OR a *documented, intentional* outcome.
/// For methylation-level calling the headline hazards are:
///   • a DivideByZeroException / NaN / Infinity when total reads = 0 at a site
///     (the β = methylated/total division) — the documented guard is to EXCLUDE
///     a zero-coverage CpG from the result, so no division ever runs on total 0;
///   • a β outside the closed interval [0, 1] (INV-04 — definition of a fraction);
///   • a NEGATIVE β (impossible from non-negative counts; would signal corruption);
///   • a NaN/Infinity β leaking from any code path.
/// — docs/ADVANCED_TESTING_CHECKLIST.md §8 "Fuzzing".
///
/// ───────────────────────────────────────────────────────────────────────────
/// Unit: EPIGEN-METHYL-001 — Methylation Analysis (Epigenetics)
/// Checklist: docs/checklists/03_FUZZING.md, row 185.
/// Fuzz strategy exercised for THIS unit:
///   • BE = Boundary Exploitation — граничні значення: 0, -1, MaxInt, empty.
///     Targets (checklist row 185): "no reads, all-methylated, zero coverage".
/// — docs/checklists/03_FUZZING.md §Description (strategy codes).
///
/// ───────────────────────────────────────────────────────────────────────────
/// The methylation-level contract under test
/// ───────────────────────────────────────────────────────────────────────────
/// Methylation call (Krueger &amp; Andrews 2011, Bismark) [2]: at a reference CpG
/// cytosine, a read base C = methylated, a read base T = unmethylated; the
/// per-CpG methylation level (β-value / Bismark fraction) is
///   level = methylated / (methylated + unmethylated)  ∈ [0, 1]   [3]
///   — docs/algorithms/Epigenetics/Methylation_Analysis.md §2.2;
///   — docs/algorithms/Epigenetics/Bisulfite_Sequencing_Analysis.md §2.2.
///
/// Documented invariant exercised here:
///   • INV-04: each per-site methylation level ∈ [0,1] (definition of a read
///     fraction). — docs/algorithms/Epigenetics/Methylation_Analysis.md §2.4.
///
/// Boundary / degenerate handling (Methylation_Analysis.md §"Degenerate";
/// Bisulfite_Sequencing_Analysis.md §"sites with no coverage are omitted"):
///   • null / empty reference                → empty enumeration (no crash).
///   • no reads at all                       → every CpG has total 0 → ALL sites
///     excluded → empty enumeration (the zero-coverage guard; NO DivideByZero).
///   • a CpG with zero coverage              → EXCLUDED from the result (its β is
///     undefined; the guard `if (total == 0) continue;` prevents meth/0).
///   • all-methylated (every valid call C)   → β = 1.0 exactly (upper bound).
///   • all-unmethylated (every valid call T) → β = 0.0 exactly (lower bound).
///   • read base ∉ {C, T} at a CpG C         → not a valid call, ignored
///     (does NOT inflate coverage, cannot push β out of [0,1]).
///   • read bases past the reference end     → ignored.
///
/// All randomness is LOCALLY seeded (new Random(seed)); no shared static Rng.
/// </summary>
[TestFixture]
[Category("Fuzzing")]
public sealed class EpigeneticMethylationLevelFuzzTests
{
    // ── Well-formed-result assertion helper ─────────────────────────────────
    // Pins the documented β contract against EVERY returned site, independent of
    // the inputs that produced it:
    //   • MethylationLevel must be a FINITE number (no NaN, no ±Infinity) — the
    //     direct evidence that the meth/total division never ran on total 0;
    //   • MethylationLevel ∈ [0, 1] (INV-04, definition of a read fraction);
    //   • Coverage > 0 — a returned site is, by the zero-coverage guard, covered;
    //   • Coverage is internally consistent: meth = round(level·coverage) is a
    //     non-negative integer in [0, coverage], i.e. β really is m/(m+u).
    // This is what stops a test from rubber-stamping a corrupted call green.
    private static void AssertWellFormedSite(MethylationSite site, string because)
    {
        double level = site.MethylationLevel;

        double.IsNaN(level).Should().BeFalse(
            $"β must never be NaN (no division by zero coverage): {because}");
        double.IsInfinity(level).Should().BeFalse(
            $"β must never be ±Infinity (no division by zero coverage): {because}");

        level.Should().BeGreaterThanOrEqualTo(0.0,
            $"INV-04: β = methylated/total ≥ 0 for non-negative counts: {because}");
        level.Should().BeLessThanOrEqualTo(1.0,
            $"INV-04: β = methylated/total ≤ 1: {because}");

        site.Coverage.Should().BeGreaterThan(0,
            $"a returned site is covered (zero-coverage sites are excluded): {because}");

        // Recover the integer methylated-read count and confirm β = m / coverage.
        long meth = (long)Math.Round(level * site.Coverage);
        meth.Should().BeGreaterThanOrEqualTo(0,
            $"methylated read count must be ≥ 0: {because}");
        meth.Should().BeLessThanOrEqualTo(site.Coverage,
            $"methylated read count must be ≤ coverage: {because}");
        ((double)meth / site.Coverage).Should().BeApproximately(level, 1e-9,
            $"β must equal methylated/total exactly: {because}");
    }

    private static void AssertAllWellFormed(IEnumerable<MethylationSite> sites, string because)
    {
        foreach (var site in sites)
            AssertWellFormedSite(site, because);
    }

    #region EPIGEN-METHYL-001 — Methylation level / β-value (CalculateMethylationFromBisulfite)

    // ── BE: empty / null reference → documented empty enumeration ────────────

    [Test]
    public void Methylation_NullReference_ReturnsEmpty_NoCrash()
    {
        // Docs §Degenerate: "Null/empty sequence → empty enumeration".
        Action act = () => CalculateMethylationFromBisulfite(null!, new[] { ("C", 0) }).ToList();
        act.Should().NotThrow();
        CalculateMethylationFromBisulfite(null!, new[] { ("C", 0) }).Should().BeEmpty();
    }

    [Test]
    public void Methylation_EmptyReference_ReturnsEmpty()
    {
        CalculateMethylationFromBisulfite("", new[] { ("C", 0) }).Should().BeEmpty();
    }

    [Test]
    public void Methylation_ReferenceWithNoCpG_ReturnsEmpty()
    {
        // No CpG dinucleotide → no callable site, regardless of the reads.
        var sites = CalculateMethylationFromBisulfite(
            "ATATATAT", new[] { ("C", 0), ("T", 3) }).ToList();
        sites.Should().BeEmpty(because: "there is no CpG cytosine to call");
    }

    // ── BE: "no reads" → all CpG sites zero-coverage → excluded ──────────────

    [Test]
    public void Methylation_NoReads_ReturnsEmpty_NoDivideByZero()
    {
        // THE zero-coverage hazard at scale: a reference full of CpGs but ZERO
        // reads. Every site has total 0 → the guard excludes every site. No
        // division by zero ever runs, and the enumeration is empty.
        string reference = string.Concat(Enumerable.Repeat("CG", 50)); // 50 CpGs
        Action act = () => CalculateMethylationFromBisulfite(
            reference, Array.Empty<(string, int)>()).ToList();
        act.Should().NotThrow();

        var sites = CalculateMethylationFromBisulfite(
            reference, Array.Empty<(string, int)>()).ToList();
        sites.Should().BeEmpty(
            because: "with no reads every CpG has zero coverage and is excluded");
    }

    [Test]
    public void Methylation_EmptyReadCollection_ReturnsEmpty()
    {
        CalculateMethylationFromBisulfite("ACGTACGT", new List<(string, int)>())
            .Should().BeEmpty();
    }

    // ── BE: a single CpG with zero coverage is excluded (guard direct hit) ────

    [Test]
    public void Methylation_UncoveredCpG_IsExcluded_NotNaN()
    {
        // Reference ACGTACGT has CpG cytosines at index 1 and 5. Cover only @1.
        // Site @5 has total 0 → it must NOT appear (β undefined), and must NOT
        // leak a NaN/Infinity placeholder site.
        var sites = CalculateMethylationFromBisulfite(
            "ACGTACGT", new[] { ("C", 1) }).ToList();

        sites.Should().NotContain(s => s.Position == 5,
            because: "a CpG with zero coverage has undefined β and is excluded");
        sites.Should().OnlyContain(s => !double.IsNaN(s.MethylationLevel)
                                     && !double.IsInfinity(s.MethylationLevel));
        AssertAllWellFormed(sites, "single covered CpG");
    }

    [Test]
    public void Methylation_OnlyNonCTCallsAtCpG_SiteHasZeroCoverage_AndIsExcluded()
    {
        // Reads that are neither C nor T at the CpG C are not valid calls, so the
        // site accrues ZERO coverage — exactly the zero-coverage path — and must
        // be excluded with no division by zero.
        var sites = CalculateMethylationFromBisulfite(
            "ACGTAA", new[] { ("A", 1), ("G", 1), ("N", 1) }).ToList();

        sites.Should().BeEmpty(
            because: "no valid C/T call → coverage 0 → site excluded, no meth/0");
    }

    // ── BE: all-methylated → β = 1.0 exactly (the upper bound) ────────────────

    [Test]
    public void Methylation_AllMethylated_BetaIsExactlyOne()
    {
        // Every covering read calls C at the CpG cytosine → meth == total → β = 1.0.
        var reads = Enumerable.Range(0, 12).Select(_ => ("C", 1)).ToArray();
        var sites = CalculateMethylationFromBisulfite("ACGTAA", reads).ToList();

        var site = sites.Single(s => s.Position == 1);
        site.MethylationLevel.Should().Be(1.0,
            because: "all 12 valid calls are C → 12/(12+0) = 1.0 (upper bound)");
        site.Coverage.Should().Be(12);
        AssertWellFormedSite(site, "all-methylated CpG");
    }

    [Test]
    public void Methylation_AllMethylated_AcrossManyCpGs_AllBetaOne()
    {
        // A run of CpGs, each fully covered by C calls. Every β must be exactly 1.0.
        string reference = string.Concat(Enumerable.Repeat("CG", 20)) + "A"; // CpGs at 0,2,4,...
        var reads = new List<(string, int)>();
        for (int pos = 0; pos < reference.Length - 1; pos += 2)
        {
            reads.Add(("C", pos)); // C at each CpG cytosine
            reads.Add(("C", pos));
        }

        var sites = CalculateMethylationFromBisulfite(reference, reads).ToList();
        sites.Should().NotBeEmpty();
        sites.Should().OnlyContain(s => s.MethylationLevel == 1.0,
            because: "every CpG is covered only by methylated (C) calls");
        AssertAllWellFormed(sites, "all-methylated multi-CpG");
    }

    // ── BE: all-unmethylated → β = 0.0 exactly (the lower bound) ──────────────

    [Test]
    public void Methylation_AllUnmethylated_BetaIsExactlyZero()
    {
        // Every covering read calls T at the CpG cytosine → meth == 0 → β = 0.0.
        var reads = Enumerable.Range(0, 8).Select(_ => ("T", 1)).ToArray();
        var sites = CalculateMethylationFromBisulfite("ACGTAA", reads).ToList();

        var site = sites.Single(s => s.Position == 1);
        site.MethylationLevel.Should().Be(0.0,
            because: "all 8 valid calls are T → 0/(0+8) = 0.0 (lower bound)");
        site.Coverage.Should().Be(8);
        AssertWellFormedSite(site, "all-unmethylated CpG");
    }

    // ── POSITIVE sanity: hand-computed β-values ──────────────────────────────

    [Test]
    public void Methylation_SevenMethylatedThreeUnmethylated_BetaIsZeroPointSeven()
    {
        // 7 C calls + 3 T calls at one CpG → β = 7/(7+3) = 0.7, coverage 10.
        var reads = Enumerable.Repeat(("C", 1), 7)
            .Concat(Enumerable.Repeat(("T", 1), 3))
            .ToArray();
        var sites = CalculateMethylationFromBisulfite("ACGTAA", reads).ToList();

        var site = sites.Single(s => s.Position == 1);
        site.MethylationLevel.Should().BeApproximately(0.7, 1e-12,
            because: "7 methylated / 10 total = 0.7 (Bismark β)");
        site.Coverage.Should().Be(10);
        AssertWellFormedSite(site, "7/3 mixed CpG");
    }

    [Test]
    public void Methylation_HalfMethylated_BetaIsExactlyHalf()
    {
        // One C and one T at a CpG → β = 1/2 = 0.5, coverage 2.
        var sites = CalculateMethylationFromBisulfite(
            "ACGTACGT", new[] { ("C", 1), ("T", 1) }).ToList();

        var site = sites.Single(s => s.Position == 1);
        site.MethylationLevel.Should().BeApproximately(0.5, 1e-12);
        site.Coverage.Should().Be(2);
        AssertWellFormedSite(site, "half-methylated CpG");
    }

    [Test]
    public void Methylation_NonCTCallIgnored_DoesNotInflateCoverageOrBeta()
    {
        // An 'A' read at the CpG C is not a valid call; only the T counts.
        var sites = CalculateMethylationFromBisulfite(
            "ACGTAA", new[] { ("A", 1), ("T", 1) }).ToList();

        var site = sites.Single(s => s.Position == 1);
        site.Coverage.Should().Be(1, because: "only the T call is a valid bisulfite call");
        site.MethylationLevel.Should().Be(0.0);
        AssertWellFormedSite(site, "non-CT ignored");
    }

    // ── BE: read coordinate boundaries (0, -1, MaxInt) ───────────────────────

    [Test]
    public void Methylation_ReadStartingBeforeReference_NegativeOffset_NoCrash()
    {
        // A read aligned to a negative start position must not throw; bases that
        // map to refPos < 0 simply do not match any CpG key.
        Action act = () => CalculateMethylationFromBisulfite(
            "ACGTACGT", new[] { ("CCCC", -2) }).ToList();
        act.Should().NotThrow();
        var sites = CalculateMethylationFromBisulfite(
            "ACGTACGT", new[] { ("CCCC", -2) }).ToList();
        AssertAllWellFormed(sites, "negative start offset");
    }

    [Test]
    public void Methylation_ReadStartingAtIntMaxValue_NoCrash()
    {
        // An absurd start offset must not overflow into a crash; it simply covers
        // nothing (the loop condition fails immediately).
        Action act = () => CalculateMethylationFromBisulfite(
            "ACGTACGT", new[] { ("CCCC", int.MaxValue) }).ToList();
        act.Should().NotThrow();
        CalculateMethylationFromBisulfite("ACGTACGT", new[] { ("CCCC", int.MaxValue) })
            .Should().BeEmpty(because: "no in-reference position is covered");
    }

    [Test]
    public void Methylation_ReadPastReferenceEnd_ExtraBasesIgnored()
    {
        // Reference length 6 (ACGTAA), CpG@1. Read 'CGGGGG' from pos 1 covers the
        // CpG@1 with a C; bases beyond the reference are ignored.
        var sites = CalculateMethylationFromBisulfite(
            "ACGTAA", new[] { ("CGGGGG", 1) }).ToList();

        var site = sites.Single(s => s.Position == 1);
        site.Coverage.Should().Be(1, because: "only the in-reference C call counts");
        site.MethylationLevel.Should().Be(1.0);
        AssertWellFormedSite(site, "read past reference end");
    }

    [Test]
    public void Methylation_EmptyReadString_NoCrash_CoversNothing()
    {
        Action act = () => CalculateMethylationFromBisulfite(
            "ACGTAA", new[] { ("", 1) }).ToList();
        act.Should().NotThrow();
        CalculateMethylationFromBisulfite("ACGTAA", new[] { ("", 1) })
            .Should().BeEmpty(because: "an empty read covers no position");
    }

    // ── BE/robustness: random fuzz — never crash, β always well-formed ───────

    [Test]
    [CancelAfter(30000)]
    public void Methylation_RandomReadsAndReferences_NeverCrash_BetaAlwaysInRange()
    {
        // Random references over an alphabet that produces CpGs, and random reads
        // with random (including out-of-range / negative) start offsets and read
        // bases drawn from C/T plus invalid bases. Every returned site must
        // satisfy INV-04 (β finite, ∈ [0,1]) and the coverage/β consistency check.
        const string refAlphabet = "ACGTacgt";
        const string readAlphabet = "CTNAGctna";
        for (int seed = 0; seed < 400; seed++)
        {
            var rng = new Random(seed);

            int refLen = rng.Next(0, 200);
            var refChars = new char[refLen];
            for (int i = 0; i < refLen; i++)
                refChars[i] = refAlphabet[rng.Next(refAlphabet.Length)];
            string reference = new string(refChars);

            int readCount = rng.Next(0, 40);
            var reads = new List<(string, int)>(readCount);
            for (int r = 0; r < readCount; r++)
            {
                int readLen = rng.Next(0, 12);
                var rc = new char[readLen];
                for (int i = 0; i < readLen; i++)
                    rc[i] = readAlphabet[rng.Next(readAlphabet.Length)];
                // Start offsets include negatives, in-range, and far past the end.
                int start = rng.Next(-5, refLen + 6);
                if (rng.Next(50) == 0) start = int.MaxValue; // extreme boundary
                reads.Add((new string(rc), start));
            }

            List<MethylationSite> sites = null!;
            Action act = () => sites = CalculateMethylationFromBisulfite(reference, reads).ToList();
            act.Should().NotThrow($"seed={seed}, refLen={refLen}, reads={readCount}");

            AssertAllWellFormed(sites, $"seed={seed}");
        }
    }

    [Test]
    [CancelAfter(30000)]
    public void Methylation_RandomPureCounts_BetaMatchesClosedFormFraction()
    {
        // For a SINGLE CpG, place m methylated (C) and u unmethylated (T) calls
        // and confirm β == m/(m+u) exactly, with m,u spanning 0..N including the
        // boundary 0 on each side (one-sided), but never both 0 (that is the
        // excluded zero-coverage case, covered separately).
        for (int seed = 0; seed < 300; seed++)
        {
            var rng = new Random(seed);
            int m = rng.Next(0, 60);
            int u = rng.Next(0, 60);
            if (m == 0 && u == 0) u = 1; // ensure non-zero coverage

            var reads = Enumerable.Repeat(("C", 1), m)
                .Concat(Enumerable.Repeat(("T", 1), u))
                .ToArray();

            var sites = CalculateMethylationFromBisulfite("ACGTAA", reads).ToList();
            var site = sites.Single(s => s.Position == 1);

            double expected = (double)m / (m + u);
            site.MethylationLevel.Should().BeApproximately(expected, 1e-12,
                $"seed={seed}: β = {m}/({m}+{u})");
            site.Coverage.Should().Be(m + u, $"seed={seed}: coverage = m+u");
            AssertWellFormedSite(site, $"seed={seed} pure counts");
        }
    }

    [Test]
    [CancelAfter(30000)]
    public void Methylation_AllMethylatedFuzz_BetaAlwaysOne()
    {
        // Fuzz the count and position: however many C calls land on a CpG, with
        // NO T calls, β must be exactly 1.0 (the upper bound) every time.
        for (int seed = 0; seed < 200; seed++)
        {
            var rng = new Random(seed);
            int n = rng.Next(1, 80);
            var reads = Enumerable.Repeat(("C", 1), n).ToArray();

            var site = CalculateMethylationFromBisulfite("ACGTAA", reads)
                .Single(s => s.Position == 1);

            site.MethylationLevel.Should().Be(1.0, $"seed={seed}: {n} C calls → β=1.0");
            site.Coverage.Should().Be(n, $"seed={seed}");
            AssertWellFormedSite(site, $"seed={seed} all-methylated");
        }
    }

    #endregion
}
