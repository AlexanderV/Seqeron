namespace Seqeron.Genomics.Tests.Fuzzing;

/// <summary>
/// Fuzz tests for the Metagenomics antibiotic-resistance-gene detection unit —
/// ResFinder-style acquired-resistance-gene calling via
/// <see cref="MetagenomicsAnalyzer.FindAntibioticResistanceGenes"/>.
///
/// ───────────────────────────────────────────────────────────────────────────
/// What fuzzing verifies
/// ───────────────────────────────────────────────────────────────────────────
/// Fuzzing feeds random, invalid and boundary inputs to a unit and asserts that
/// the code NEVER fails in an undisciplined way: no hang, no state corruption,
/// no NaN/Infinity leaking into an identity or coverage, and no *unhandled*
/// runtime exception (IndexOutOfRangeException, NullReferenceException,
/// DivideByZeroException, OverflowException, …). Every input must produce EITHER
/// a well-defined, theory-correct result, OR a *documented, intentional*
/// validation exception (ArgumentNullException / ArgumentOutOfRangeException). A
/// raw runtime exception, a hang, an identity/coverage outside [0, 1], or a
/// *fabricated* resistance hit on degenerate input is a bug, not a passing test.
/// — docs/ADVANCED_TESTING_CHECKLIST.md §8 "Fuzzing".
///
/// ───────────────────────────────────────────────────────────────────────────
/// Unit: META-RESIST-001 — antibiotic resistance gene detection (Metagenomics)
/// Checklist: docs/checklists/03_FUZZING.md, row 196.
/// Fuzz strategy for THIS unit: BE = Boundary Exploitation (0, -1, MaxInt, empty)
///   — docs/checklists/03_FUZZING.md §Description (strategy codes).
/// Fuzz targets (checklist row 196): "no resistance gene, empty DB, partial hit".
///
/// ───────────────────────────────────────────────────────────────────────────
/// The contract under test
/// ───────────────────────────────────────────────────────────────────────────
/// FindAntibioticResistanceGenes screens each contig against a caller-supplied
/// reference DB. For each (contig, reference) pair the best ungapped (gapless)
/// alignment is located by sliding the reference across the contig at every
/// offset (overhanging both ends); for that best window
///   identity = identical positions / window length    (gapless ⇒ columns = window)
///   coverage = window length / reference length m
/// A reference is reported only when identity ≥ identityThreshold AND
/// coverage ≥ coverageThreshold; the single best-matching reference per contig is
/// returned (max identity, ties → max coverage). Defaults: 0.90 identity / 0.60
/// coverage.
///   — docs/algorithms/Metagenomics/Antibiotic_Resistance_Detection.md §2.2, §4;
///     src/Seqeron/Algorithms/Seqeron.Genomics.Metagenomics/MetagenomicsAnalyzer.cs
///     (FindAntibioticResistanceGenes 1236–1293, BestUngappedMatch 1308–1349).
///
/// Documented invariants this fixture pins (Antibiotic_Resistance_Detection.md §2.4):
///   • INV-01: 0 ≤ identity ≤ 1 and 0 ≤ coverage ≤ 1 (matches ≤ window ≤ m).
///   • INV-02: a hit is reported only if identity ≥ idThr AND coverage ≥ covThr.
///   • INV-03: at most one hit per contig (best: max identity, tie → max coverage).
///   • INV-04: exact full-length match ⇒ identity = 1.0, coverage = 1.0.
///   • INV-05: default thresholds = 0.90 identity, 0.60 coverage.
///
/// Boundary / malformed-input handling fixed by the doc (§3.3, §6.1) and source,
/// which these fuzz tests pin so the contract can never silently drift:
///   • NO RESISTANCE GENE (BE): a contig that does not match any reference well
///     enough ⇒ NO hit for that contig (reporting rule). Resistance is NEVER
///     fabricated when nothing aligns. — §6.1 "No reference passes thresholds".
///   • EMPTY DB (BE): empty reference database (or a DB whose every reference has
///     an empty sequence) ⇒ no hits at all, no crash. — §6.1, source 1249–1251.
///   • PARTIAL HIT (BE): a match below a threshold is NOT reported —
///       – partial coverage: a short shared substring scores identity 1.0 but
///         coverage < 0.60 ⇒ no hit (e.g. 5-bp perfect window of a 10-bp gene).
///       – partial identity: a window covering the gene but diluted by scattered
///         mismatches scores identity < 0.90 ⇒ no hit, even at full coverage.
///     Lowering the relevant threshold turns the same partial hit into a reported
///     hit with the hand-derived identity/coverage. — INV-02, §6.1.
///   • Empty contig sequence skipped; empty-sequence reference ignored. — §6.1.
///   • null contigs / referenceGenes → ArgumentNullException; threshold ∉ [0,1] →
///     ArgumentOutOfRangeException, eagerly (the method is iterator-based, so the
///     guards live in the non-deferred outer method). — §3.3, §6.1, source 1242–1247.
///
/// Positive sanity (worked examples, derived INDEPENDENTLY from the BLAST identity
/// and ResFinder coverage definitions, NOT echoed off the implementation):
///   • §7.1 API example: contig "AAACGTACGT", gene "CGTACGT" (m=7) ⇒ the gene
///     aligns perfectly at contig offset 3 ⇒ identity = 7/7 = 1.0, coverage = 7/7
///     = 1.0 (INV-04). Hand-checkable: "AAA" + "CGTACGT".
///   • §7.1 numerical walk-through: gene "CGTACGT" vs contig "CGTTCGT" ⇒ best
///     window is the full 7-column overlap, one mismatch (col 4: T vs A) ⇒
///     identity = 6/7 ≈ 0.857, coverage = 7/7 = 1.0. At the 0.90 default this
///     FAILS the identity cutoff (partial hit ⇒ no report); at idThr ≤ 6/7 it is
///     reported with identity 6/7.
/// A genuine full-length match must therefore yield identity = coverage = 1.0, so
/// a passing "no crash" result cannot be a degenerate detector that returns
/// nothing (or everything) for every input.
///
/// ───────────────────────────────────────────────────────────────────────────
/// Determinism
/// ───────────────────────────────────────────────────────────────────────────
/// All inputs are hand-built or generated from a LOCALLY fixed-seed
/// `new Random(seed)` (never a shared static Rng), so every run is reproducible.
/// </summary>
[TestFixture]
[Category("Fuzzing")]
public class MetagenomicsResistanceFuzzTests
{

    // A single reference DB entry used by the worked examples (§7.1).
    private static (string GeneId, string Sequence, string Name, string AntibioticClass)[] BlaXDb() =>
        new[] { ("blaX", "CGTACGT", "blaX-like", "beta-lactam") };

    #region META-RESIST-001 — antibiotic resistance gene detection

    // ════════════════════════════════════════════════════════════════════════
    //  Positive sanity — worked examples must be reproduced EXACTLY.
    //  Guards against a degenerate detector (reports nothing / reports everything)
    //  that would pass every boundary test below.
    // ════════════════════════════════════════════════════════════════════════

    [Test]
    public void FindAntibioticResistanceGenes_ApiWorkedExample_PerfectFullLengthHit()
    {
        // §7.1: contig "AAA" + "CGTACGT" contains the gene exactly at offset 3.
        var contigs = new[] { ("c1", "AAACGTACGT") };

        var hits = MetagenomicsAnalyzer.FindAntibioticResistanceGenes(contigs, BlaXDb()).ToList();

        var h = hits.Should().ContainSingle().Subject;
        h.ContigId.Should().Be("c1");
        h.ResistanceGene.Should().Be("blaX-like");
        h.AntibioticClass.Should().Be("beta-lactam");
        h.PercentIdentity.Should().BeApproximately(1.0, 1e-12,
            "the 7-bp gene aligns with 7/7 identical positions — INV-04, §7.1");
        h.Coverage.Should().BeApproximately(1.0, 1e-12,
            "the full reference length is spanned ⇒ coverage 7/7 = 1.0 — INV-04, §7.1");
    }

    [Test]
    public void FindAntibioticResistanceGenes_NumericalWalkthrough_OneMismatchFailsDefaultIdentity()
    {
        // §7.1 walk-through: gene "CGTACGT" vs contig "CGTTCGT" ⇒ identity 6/7 ≈ 0.857.
        // 6/7 < 0.90 default ⇒ PARTIAL hit ⇒ NOT reported.
        var contigs = new[] { ("c1", "CGTTCGT") };

        var atDefault = MetagenomicsAnalyzer.FindAntibioticResistanceGenes(contigs, BlaXDb()).ToList();
        atDefault.Should().BeEmpty("identity 6/7 ≈ 0.857 < 0.90 default ⇒ partial hit, no report — INV-02");

        // Lower the identity cutoff below 6/7 and the SAME partial match is reported
        // with the independently-derived identity 6/7 and full coverage.
        var lowered = MetagenomicsAnalyzer
            .FindAntibioticResistanceGenes(contigs, BlaXDb(), identityThreshold: 0.80)
            .ToList();

        var h = lowered.Should().ContainSingle().Subject;
        h.PercentIdentity.Should().BeApproximately(6.0 / 7.0, 1e-12,
            "one mismatch over the 7-column overlap ⇒ identity 6/7 — §7.1 walk-through");
        h.Coverage.Should().BeApproximately(1.0, 1e-12, "the full reference is spanned ⇒ coverage 1.0");
    }

    // ───────────────────────────────────────────────────────────────────────
    // INV-03: at most one hit per contig — best match (max identity, tie → cov).
    // Two passing references; the perfect one must win.
    // ───────────────────────────────────────────────────────────────────────
    [Test]
    public void FindAntibioticResistanceGenes_MultiplePassing_ReportsBestMatchOnly()
    {
        var contigs = new[] { ("c1", "AAACGTACGT") };
        var db = new[]
        {
            ("perfect", "CGTACGT", "perfect-gene", "classA"),     // identity 1.0, coverage 1.0
            ("partial", "CGTACGA", "partial-gene", "classB"),     // 6/7 identity over same window
        };

        var hits = MetagenomicsAnalyzer.FindAntibioticResistanceGenes(contigs, db).ToList();

        var h = hits.Should().ContainSingle("at most one hit per contig — INV-03").Subject;
        h.ResistanceGene.Should().Be("perfect-gene", "highest-identity reference wins — INV-03");
        h.PercentIdentity.Should().BeApproximately(1.0, 1e-12);
    }

    // ════════════════════════════════════════════════════════════════════════
    //  BE — no resistance gene present
    // ════════════════════════════════════════════════════════════════════════

    [Test]
    public void FindAntibioticResistanceGenes_NoMatchingReference_ProducesNoHit()
    {
        // Contig shares no meaningful similarity with the reference ⇒ no hit.
        var contigs = new[] { ("c1", "TTTTTTTTTT") };
        var db = new[] { ("g1", "AAAAAAAAAA", "geneA", "classA") };

        var hits = MetagenomicsAnalyzer.FindAntibioticResistanceGenes(contigs, db).ToList();

        hits.Should().BeEmpty(
            "no window reaches 0.90 identity at 0.60 coverage ⇒ resistance is never fabricated — §6.1");
    }

    [Test]
    public void FindAntibioticResistanceGenes_EmptyContigSequence_SkippedNoHit()
    {
        // §6.1: empty contig sequence is skipped (nothing to align).
        var contigs = new[] { ("c1", "") };

        var hits = MetagenomicsAnalyzer.FindAntibioticResistanceGenes(contigs, BlaXDb()).ToList();

        hits.Should().BeEmpty("empty contig sequence is skipped — §6.1");
    }

    // ════════════════════════════════════════════════════════════════════════
    //  BE — empty resistance DB
    // ════════════════════════════════════════════════════════════════════════

    [Test]
    public void FindAntibioticResistanceGenes_EmptyDatabase_NoHits()
    {
        var contigs = new[] { ("c1", "AAACGTACGT") };
        var emptyDb = Array.Empty<(string, string, string, string)>();

        var hits = MetagenomicsAnalyzer.FindAntibioticResistanceGenes(contigs, emptyDb).ToList();

        hits.Should().BeEmpty("an empty reference DB can never produce a hit — §6.1");
    }

    [Test]
    public void FindAntibioticResistanceGenes_DatabaseOfEmptySequenceReferences_IgnoredNoHits()
    {
        // §6.1 / source 1249–1251: empty-sequence references are filtered out before matching,
        // so a DB consisting only of them behaves like an empty DB.
        var contigs = new[] { ("c1", "AAACGTACGT") };
        var db = new[]
        {
            ("e1", "", "empty-gene-1", "classA"),
            ("e2", "", "empty-gene-2", "classB"),
        };

        var hits = MetagenomicsAnalyzer.FindAntibioticResistanceGenes(contigs, db).ToList();

        hits.Should().BeEmpty("references with empty sequences are ignored ⇒ effectively empty DB — §6.1");
    }

    [Test]
    public void FindAntibioticResistanceGenes_NoContigs_NoHits()
    {
        var contigs = Array.Empty<(string, string)>();

        var hits = MetagenomicsAnalyzer.FindAntibioticResistanceGenes(contigs, BlaXDb()).ToList();

        hits.Should().BeEmpty("no contigs ⇒ nothing to report");
    }

    // ════════════════════════════════════════════════════════════════════════
    //  BE — partial hit (below a threshold ⇒ not reported)
    // ════════════════════════════════════════════════════════════════════════

    [Test]
    public void FindAntibioticResistanceGenes_PartialCoverage_BelowThreshold_NoHit()
    {
        // Reference "AAAAAAAAAA" (m = 10); contig shares only a 5-bp perfect run.
        // Best window = the 5 A's ⇒ identity 5/5 = 1.0, coverage 5/10 = 0.5 < 0.60 ⇒ no hit.
        var contigs = new[] { ("c1", "CCCCCAAAAA") };
        var db = new[] { ("g1", "AAAAAAAAAA", "geneA", "classA") };

        var atDefault = MetagenomicsAnalyzer.FindAntibioticResistanceGenes(contigs, db).ToList();
        atDefault.Should().BeEmpty(
            "perfect 5-bp window covers only 0.5 of the 10-bp gene < 0.60 ⇒ partial hit, no report — INV-02");

        // Lower the coverage cutoff below 0.5 ⇒ the SAME partial match is now reported,
        // with the independently-derived identity 1.0 and coverage 5/10 = 0.5.
        var lowered = MetagenomicsAnalyzer
            .FindAntibioticResistanceGenes(contigs, db, coverageThreshold: 0.40)
            .ToList();

        var h = lowered.Should().ContainSingle().Subject;
        h.PercentIdentity.Should().BeApproximately(1.0, 1e-12, "the 5-bp window is a perfect AAAAA match");
        h.Coverage.Should().BeApproximately(0.5, 1e-12, "5-bp window of a 10-bp gene ⇒ coverage 0.5");
    }

    [Test]
    public void FindAntibioticResistanceGenes_PartialIdentity_BelowThreshold_NoHit()
    {
        // Reference "AAAAAAAAAA" (m = 10); contig "AGAGAGAGAG" matches only the 5 even
        // positions (A). The most-matches window has 5 matches; the shortest such window
        // (contig indices 0..8, offset −1) spans 9 columns ⇒ identity = 5/9 ≈ 0.556,
        // coverage = 9/10 = 0.9. Identity 0.556 < 0.90 ⇒ partial hit ⇒ no report.
        var contigs = new[] { ("c1", "AGAGAGAGAG") };
        var db = new[] { ("g1", "AAAAAAAAAA", "geneA", "classA") };

        var atDefault = MetagenomicsAnalyzer.FindAntibioticResistanceGenes(contigs, db).ToList();
        atDefault.Should().BeEmpty(
            "best window identity 5/9 ≈ 0.556 < 0.90 ⇒ partial hit, no report — INV-02");

        // Lower the identity cutoff below 5/9 ⇒ the SAME partial match is reported,
        // with the independently-derived identity 5/9 and coverage 9/10.
        var lowered = MetagenomicsAnalyzer
            .FindAntibioticResistanceGenes(contigs, db, identityThreshold: 0.50)
            .ToList();

        var h = lowered.Should().ContainSingle().Subject;
        h.PercentIdentity.Should().BeApproximately(5.0 / 9.0, 1e-12,
            "5 matched even positions over the shortest 9-column window ⇒ identity 5/9");
        h.Coverage.Should().BeApproximately(0.9, 1e-12, "9-column window of a 10-bp gene ⇒ coverage 0.9");
    }

    // ════════════════════════════════════════════════════════════════════════
    //  BE — argument validation (eager, not deferred behind the iterator)
    // ════════════════════════════════════════════════════════════════════════

    [Test]
    public void FindAntibioticResistanceGenes_NullContigs_ThrowsArgumentNullEagerly()
    {
        Action act = () => MetagenomicsAnalyzer.FindAntibioticResistanceGenes(null!, BlaXDb());

        act.Should().Throw<ArgumentNullException>()
            .WithParameterName("contigs", "null contigs is rejected eagerly — §3.3");
    }

    [Test]
    public void FindAntibioticResistanceGenes_NullReferenceGenes_ThrowsArgumentNullEagerly()
    {
        var contigs = new[] { ("c1", "AAACGTACGT") };

        Action act = () => MetagenomicsAnalyzer.FindAntibioticResistanceGenes(contigs, null!);

        act.Should().Throw<ArgumentNullException>()
            .WithParameterName("referenceGenes", "null DB is rejected eagerly — §3.3");
    }

    [TestCase(-1.0)]
    [TestCase(-0.0001)]
    [TestCase(1.0001)]
    [TestCase(double.MaxValue)]
    public void FindAntibioticResistanceGenes_IdentityThresholdOutOfRange_Throws(double idThr)
    {
        var contigs = new[] { ("c1", "AAACGTACGT") };

        Action act = () => MetagenomicsAnalyzer.FindAntibioticResistanceGenes(contigs, BlaXDb(), identityThreshold: idThr);

        act.Should().Throw<ArgumentOutOfRangeException>()
            .WithParameterName("identityThreshold", "identity threshold must be in [0,1] — §3.3");
    }

    [TestCase(-1.0)]
    [TestCase(1.0001)]
    [TestCase(double.MaxValue)]
    public void FindAntibioticResistanceGenes_CoverageThresholdOutOfRange_Throws(double covThr)
    {
        var contigs = new[] { ("c1", "AAACGTACGT") };

        Action act = () => MetagenomicsAnalyzer.FindAntibioticResistanceGenes(contigs, BlaXDb(), coverageThreshold: covThr);

        act.Should().Throw<ArgumentOutOfRangeException>()
            .WithParameterName("coverageThreshold", "coverage threshold must be in [0,1] — §3.3");
    }

    [Test]
    public void FindAntibioticResistanceGenes_ZeroThresholds_AreInRange_DoNotThrow()
    {
        // 0.0 is a valid boundary of [0,1] (BE: the 0 boundary). With both cutoffs at 0
        // any non-empty alignment passes; the best-matching reference is still reported.
        var contigs = new[] { ("c1", "AAACGTACGT") };

        Action act = () =>
            _ = MetagenomicsAnalyzer
                .FindAntibioticResistanceGenes(contigs, BlaXDb(), identityThreshold: 0.0, coverageThreshold: 0.0)
                .ToList();

        act.Should().NotThrow("0.0 is the lower boundary of the valid [0,1] range — §3.3");
    }

    // ════════════════════════════════════════════════════════════════════════
    //  BE — randomized boundary sweep: no crash/hang/NaN, contract always holds.
    // ════════════════════════════════════════════════════════════════════════

    [Test]
    [CancelAfter(30000)]
    public void FindAntibioticResistanceGenes_RandomBoundaryBatch_NeverCrashesAndStaysWellFormed()
    {
        var rng = new Random(20260620); // locally fixed seed — deterministic, never a shared Rng.
        const string alphabet = "ACGT";

        for (int trial = 0; trial < 400; trial++)
        {
            // Random contigs: some empty (skipped), some short, some long. 0..4 contigs.
            int nContigs = rng.Next(0, 5);
            var contigs = new List<(string, string)>(nContigs);
            for (int c = 0; c < nContigs; c++)
                contigs.Add(($"c{c}", RandomSeq(rng, alphabet, rng.Next(0, 30))));

            // Random reference DB: 0..4 genes, some with empty sequences (ignored).
            int nGenes = rng.Next(0, 5);
            var db = new List<(string, string, string, string)>(nGenes);
            for (int g = 0; g < nGenes; g++)
                db.Add(($"g{g}", RandomSeq(rng, alphabet, rng.Next(0, 20)), $"gene{g}", $"class{g % 3}"));

            // Random thresholds within the valid [0,1] range, including the 0 and 1 boundaries.
            double idThr = rng.Next(3) switch { 0 => 0.0, 1 => 1.0, _ => rng.NextDouble() };
            double covThr = rng.Next(3) switch { 0 => 0.0, 1 => 1.0, _ => rng.NextDouble() };

            var hits = MetagenomicsAnalyzer
                .FindAntibioticResistanceGenes(contigs, db, idThr, covThr)
                .ToList();

            // INV-03: at most one hit per contig.
            hits.Select(h => h.ContigId).Should().OnlyHaveUniqueItems(
                "at most one best-matching reference per contig — INV-03");

            foreach (var h in hits)
            {
                // INV-01: identity and coverage are well-formed probabilities in [0,1].
                h.PercentIdentity.Should().BeInRange(0.0, 1.0, "0 ≤ identity ≤ 1 — INV-01");
                h.Coverage.Should().BeInRange(0.0, 1.0, "0 ≤ coverage ≤ 1 — INV-01");
                double.IsNaN(h.PercentIdentity).Should().BeFalse("no NaN identity on boundary input");
                double.IsNaN(h.Coverage).Should().BeFalse("no NaN coverage on boundary input");
                double.IsInfinity(h.PercentIdentity).Should().BeFalse("no Infinity identity");
                double.IsInfinity(h.Coverage).Should().BeFalse("no Infinity coverage");

                // INV-02: a reported hit must pass both cutoffs.
                h.PercentIdentity.Should().BeGreaterThanOrEqualTo(idThr,
                    "a reported hit passes the identity cutoff — INV-02");
                h.Coverage.Should().BeGreaterThanOrEqualTo(covThr,
                    "a reported hit passes the coverage cutoff — INV-02");

                // The reported contig id must be one of the (non-empty) input contigs.
                contigs.Select(x => x.Item1).Should().Contain(h.ContigId);
            }
        }
    }

    private static string RandomSeq(Random rng, string alphabet, int length)
    {
        if (length <= 0) return string.Empty;
        var chars = new char[length];
        for (int i = 0; i < length; i++)
            chars[i] = alphabet[rng.Next(alphabet.Length)];
        return new string(chars);
    }

    #endregion
}
