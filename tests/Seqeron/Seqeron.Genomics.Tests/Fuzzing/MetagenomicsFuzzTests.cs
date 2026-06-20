using System;
using System.Collections.Generic;
using System.Linq;
using NUnit.Framework;
using FluentAssertions;
using Seqeron.Genomics.Metagenomics;

namespace Seqeron.Genomics.Tests;

/// <summary>
/// Fuzz tests for the Metagenomics area — Kraken-style taxonomic read classification.
///
/// ───────────────────────────────────────────────────────────────────────────
/// What fuzzing verifies
/// ───────────────────────────────────────────────────────────────────────────
/// Fuzzing feeds random, invalid and boundary inputs to a unit and asserts that
/// the code NEVER fails in an undisciplined way: no hang, no state corruption,
/// and no *unhandled* runtime exception (IndexOutOfRangeException,
/// NullReferenceException, DivideByZeroException, OverflowException, …). Every
/// input must result in EITHER a well-defined, theory-correct value, OR a
/// *documented, intentional* validation exception (ArgumentNullException /
/// ArgumentOutOfRangeException). A raw runtime exception, a hang, or a *false
/// taxon assignment* on garbage input is a bug, not a passing test.
/// — docs/ADVANCED_TESTING_CHECKLIST.md §8 "Fuzzing".
///
/// ───────────────────────────────────────────────────────────────────────────
/// Unit: META-CLASS-001 — metagenomic read classification (Metagenomics)
/// Checklist: docs/checklists/03_FUZZING.md, row 53.
/// Fuzz strategies exercised for THIS unit:
///   • MC = Malformed Content — non-DNA reads (digits, IUPAC ambiguity codes,
///          gaps, unicode), reads that hit nothing in the database.
///   • BE = Boundary Exploitation — empty reference database, empty read,
///          extremely short read (shorter than k → no k-mers can be formed).
/// — docs/checklists/03_FUZZING.md §Description (strategy codes).
///
/// ───────────────────────────────────────────────────────────────────────────
/// The classification contract under test
/// ───────────────────────────────────────────────────────────────────────────
/// MetagenomicsAnalyzer.ClassifyReads implements the Kraken (Wood &amp; Salzberg
/// 2014, Genome Biology 15:R46) k-mer / lowest-common-ancestor classifier. Per
/// read it collects canonical k-mer hits against a canonical-k-mer → taxon-id
/// database, builds the classification tree over the hit taxa and their
/// ancestors weighted by k-mer count, finds the maximum-scoring root-to-leaf
/// (RTL) path, and assigns that path's leaf (LCA-of-leaves on ties).
///   — docs/algorithms/Metagenomics/Taxonomic_Classification.md §2–3;
///     src/Seqeron/Algorithms/Seqeron.Genomics.Metagenomics/MetagenomicsAnalyzer.cs
///     (ClassifyReads / ClassifyRead, lines 191–292).
///
/// Boundary / malformed-input handling fixed by the algorithm doc and source,
/// which these fuzz tests pin so they can never silently drift:
///   • Empty reference database → no canonical k-mer can ever match → the read's
///     classification tree is empty → it is assigned the root / Unclassified
///     taxon (TaxonomyTree.RootId). NOT a crash and — KEY — NOT a false taxon
///     assignment. (ClassifyRead, hitCounts.Count == 0 branch, line 241.)
///     — Taxonomic_Classification.md INV-04, §"Edge cases".
///   • Read with no database match → identical to the empty-database case:
///     Unclassified (root). Q (TotalKmers) is still counted for the
///     non-ambiguous k-mers that were queried.
///     — Taxonomic_Classification.md §"Edge cases" ("No database matches").
///   • Non-DNA read → ambiguous k-mers (those containing any non-A/C/G/T symbol
///     after upper-casing) are SKIPPED: not counted in Q, never matched. A read
///     made entirely of such symbols therefore queries zero k-mers and is
///     Unclassified (root). Handled, not rejected, not crashed.
///     (ClassifyRead, IsAcgtOnly guard, lines 228–229.)
///     — Taxonomic_Classification.md §"Edge cases" (ambiguous k-mers).
///   • Extremely short read (shorter than k) → no k-mer of length k can be
///     formed → the `sequence.Length >= k` guard (line 222) short-circuits the
///     whole k-mer loop → Q = 0, no hits, Unclassified (root). KEY div-by-zero /
///     empty-window boundary: Confidence is 0 when Q = 0 (BuildResult, line 298,
///     `totalKmers > 0 ? … : 0.0` — never divides by zero).
///     — Taxonomic_Classification.md §"Edge cases" (empty / shorter-than-k read).
///   • Argument validation: null reads / database / taxonomy → ArgumentNullException;
///     k ≤ 0 → ArgumentOutOfRangeException. These are *documented, intentional*
///     guards (lines 197–200). NB: ClassifyReads is a deferred iterator, so the
///     guards fire on the FIRST enumeration, not at the call itself.
///
/// Positive sanity: a read that exactly contains a database taxon's k-mers IS
/// classified to that taxon (or its correct ancestor), with Q &gt; 0 and a
/// positive RTL score — so a passing "no crash" result cannot be a classifier
/// that simply returns Unclassified for everything.
///
/// ───────────────────────────────────────────────────────────────────────────
/// Unit: META-ALPHA-001 — alpha diversity (Metagenomics)
/// Checklist: docs/checklists/03_FUZZING.md, row 55.
/// Fuzz strategy for THIS unit: BE = Boundary Exploitation.
/// Fuzz targets (checklist row 55): 0 species, 1 species, all equal abundance,
/// single sample with 0 (all-zero counts).
///
/// MetagenomicsAnalyzer.CalculateAlphaDiversity takes ONE taxon→abundance map
/// and reports six within-sample metrics. The boundary contract these fuzz tests
/// pin (Alpha_Diversity.md §2.2, §2.4, §6.1; CalculateAlphaDiversity / Calculate
/// ShannonIndex / CalculateSimpsonIndex, lines 480–536):
///   • 0 species / null / empty map → an all-zero AlphaDiversity record. KEY: no
///     DivideByZero and no log(0) → NaN — the `Count == 0` guard short-circuits
///     before any normalization (line 482). — §3.3, §6.1 (Empty input).
///   • single sample with 0 (every abundance ≤ 0, e.g. all-zero counts) → the
///     `Where(v => v > 0)` filter empties the value list ⇒ observedSpecies == 0 ⇒
///     the second guard returns the all-zero record (lines 485–489). KEY: the
///     p·ln(p) term is NEVER evaluated for p = 0, so no log(0) NaN ever appears.
///   • 1 species (one positive abundance) → ShannonIndex = 0 (entropy of a one-
///     category distribution, INV-02), SimpsonIndex = 1 (Σpᵢ² with a single pᵢ=1,
///     INV-03), InverseSimpson = 1, ObservedSpecies = 1, PielouEvenness = 0
///     (undefined for S ≤ 1, INV-05). — §6.1 (Single species).
///   • all-equal abundance over S species → ShannonIndex is MAXIMAL = ln(S)
///     (uniform maximizes Shannon at fixed richness — KEY), SimpsonIndex = 1/S
///     (NOT 1 − 1/S: this implementation reports the Simpson CONCENTRATION λ =
///     Σpᵢ², §2.2), InverseSimpson = S, PielouEvenness = 1. — §6.1 (equal-
///     abundance rows). Counts or proportions both work: Shannon/Simpson
///     normalize internally by the total positive abundance (§3.3).
///
/// NB — Simpson convention: unlike the 1 − Σpᵢ² "diversity" form, the SimpsonIndex
/// field here is Simpson's concentration λ = Σpᵢ² (Simpson, 1949; §2.2). Hence
/// single-species → 1 and uniform-over-S → 1/S, the reciprocals of the classic
/// diversity values. These tests assert the implementation's documented contract.
///
/// Positive sanity: a known asymmetric two-species sample with counts (3, 1)
/// ⇒ p = (0.75, 0.25), pinned EXACTLY against the natural-log formula:
/// Shannon H = −(0.75 ln 0.75 + 0.25 ln 0.25) ≈ 0.5623351446188083, Simpson
/// λ = 0.75² + 0.25² = 0.625 — so a passing "no crash" result cannot be a
/// degenerate diversity function that returns 0 (or 1) for everything.
///   — docs/algorithms/Metagenomics/Alpha_Diversity.md §2.2, §2.4, §6.1.
///
/// ───────────────────────────────────────────────────────────────────────────
/// Unit: META-BETA-001 — beta diversity (Metagenomics)
/// Checklist: docs/checklists/03_FUZZING.md, row 56.
/// Fuzz strategy for THIS unit: BE = Boundary Exploitation.
/// Fuzz targets (checklist row 56): identical samples, empty samples,
/// single-species samples.
///
/// MetagenomicsAnalyzer.CalculateBetaDiversity takes TWO taxon→abundance maps
/// and reports BETWEEN-sample dissimilarity as Bray-Curtis (abundance-sensitive)
/// and Jaccard (presence/absence) distances, plus shared / sample-specific taxon
/// counts. A taxon is "present" iff its abundance is strictly &gt; 0; missing keys
/// default to 0. The boundary contract these fuzz tests pin (Beta_Diversity.md
/// §2, §3.3, §6.1; CalculateBetaDiversity / CalculateBrayCurtis /
/// CalculateJaccardDistance, lines 572–631):
///   • Identical samples → d(a,a) = 0 for BOTH metrics (KEY): shared abundance
///     equals half the total ⇒ Bray-Curtis = 1 − 2·(S/2)/S = 0 (INV-BRAY-02),
///     and |A∩B| = |A∪B| ⇒ Jaccard = 0 (INV-JACCARD-02). — §6.1 (Identical).
///   • Empty samples / all-zero abundances → the union has no positively-present
///     taxon and the summed abundance is 0. KEY div-by-zero boundary: Bray-Curtis
///     guards `sumTotal &gt; 0 ? … : 0` (line 624) and Jaccard guards `total &gt; 0
///     ? … : 0` (line 630), so BOTH return 0 with NO DivideByZero on the empty
///     union. — §3.3, §6.1 (Empty dictionaries or all-zero abundances).
///   • Single-species samples →
///       – same single species, equal abundance → identical ⇒ both 0;
///       – same single species, different abundance → Jaccard = 0 (same set) but
///         Bray-Curtis = |a−b|/(a+b) &gt; 0 (abundance-sensitive, §2.A);
///       – two DIFFERENT single species → disjoint non-empty sets ⇒ no shared
///         positive taxon ⇒ Bray-Curtis = 1 (INV-BRAY-03) and Jaccard = 1
///         (INV-JACCARD-03), the maximal dissimilarity. — §6.1 (No shared taxa).
///   • Symmetry: d(a,b) = d(b,a) for both metrics (the formulas are symmetric in
///     the two samples — §2). All dissimilarities lie in [0, 1] (INV-BRAY-01,
///     INV-JACCARD-01). UniFracDistance is a hard-coded 0 placeholder (§5.4).
///
/// NB — null guards: CalculateBetaDiversity does NOT null-check its sample maps
/// (it dereferences `sample1.Keys` directly, line 578); a null map would throw a
/// raw NullReferenceException. That is outside the documented BE contract for this
/// row (empty/all-zero is the boundary, not null), so these tests pin the
/// documented empty-map behavior and do not assert on null inputs.
///
/// Positive sanity: a known asymmetric pair pins EXACT Bray-Curtis and Jaccard
/// against the documented formulas — sample1 {A:3, B:1} vs sample2 {A:1, C:3}:
///   Bray-Curtis = 1 − 2·min-sum/total-sum = 1 − 2·1/8 = 0.75
///     (Σmin = min(3,1)+min(1,0)+min(0,3) = 1; Σtotal = 4 + 4 = 8),
///   Jaccard = 1 − shared/(shared+u1+u2) = 1 − 1/3 ≈ 0.6666… (A shared, B & C
///     each unique) — so a passing "no crash" result cannot be a degenerate
///     dissimilarity that returns 0 (or 1) for everything. Symmetry and
///     identical→0 are checked alongside. — Beta_Diversity.md §2, §6.1.
///
/// ───────────────────────────────────────────────────────────────────────────
/// Determinism
/// ───────────────────────────────────────────────────────────────────────────
/// All inputs are either hand-built or generated from a LOCALLY fixed-seed
/// `new Random(seed)` (no shared static Rng). k = 4 throughout so every k-mer is
/// short enough to reason about by hand and the canonical mapping is controlled.
/// The taxonomy and database below are the same hand-built fixtures used by the
/// META-CLASS-001 unit / metamorphic suites, so the contract is pinned against a
/// known-correct reference.
/// </summary>
[TestFixture]
[Category("Fuzzing")]
public class MetagenomicsFuzzTests
{
    private const int K = 4;

    // ------------------------------------------------------------------
    // Hand-built taxonomy (NCBI-style ids), mirroring the META-CLASS-001 unit suite:
    //
    //                         root(1)
    //              ┌────────────┴───────────────┐
    //         Archaea(5)                    Bacteria(2) [domain]
    //          [domain]                         │
    //                                  Proteobacteria(3) [phylum]
    //                                           │
    //                              Gammaproteobacteria(4) [class]
    //                                           │
    //                            Enterobacteriaceae(10) [family]
    //                              ┌────────────┴───────────┐
    //                       Escherichia(20)[genus]   Salmonella(21)[genus]
    //                              │                        │
    //                       E.coli(100)[species]   S.enterica(200)[species]
    // ------------------------------------------------------------------
    private static TaxonomyTree BuildTaxonomy() => new(new[]
    {
        new TaxonNode(1,   "root",                "root",    1),
        new TaxonNode(5,   "Archaea",             "domain",  1),
        new TaxonNode(2,   "Bacteria",            "domain",  1),
        new TaxonNode(3,   "Proteobacteria",      "phylum",  2),
        new TaxonNode(4,   "Gammaproteobacteria", "class",   3),
        new TaxonNode(10,  "Enterobacteriaceae",  "family",  4),
        new TaxonNode(20,  "Escherichia",         "genus",   10),
        new TaxonNode(21,  "Salmonella",          "genus",   10),
        new TaxonNode(100, "Escherichia coli",    "species", 20),
        new TaxonNode(200, "Salmonella enterica", "species", 21),
    });

    // A reference sequence whose k-mers are wholly owned by E.coli(100). Built from a
    // simple repeat so its canonical 4-mers are easy to materialize and unique to one taxon.
    private const string EcoliReference =
        "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";

    private static Dictionary<string, int> BuildEcoliDatabase(TaxonomyTree taxonomy) =>
        MetagenomicsAnalyzer.BuildKmerDatabase(
            new[] { (100, EcoliReference) }, taxonomy, K);

    #region META-CLASS-001 — read classification

    // ───────────────────────────────────────────────────────────────────────
    // Fuzz target: EMPTY DATABASE.
    // No canonical k-mer can ever resolve to a taxon ⇒ empty classification tree
    // ⇒ root / Unclassified. KEY: must NOT be a crash and must NOT be a false
    // taxon assignment. — Taxonomic_Classification.md INV-04.
    // ───────────────────────────────────────────────────────────────────────
    [Test]
    public void ClassifyReads_EmptyDatabase_AssignsRootNotAFalseTaxon()
    {
        var taxonomy = BuildTaxonomy();
        var emptyDb = new Dictionary<string, int>();

        // A perfectly valid DNA read that WOULD match in a populated database.
        var reads = new[] { ("read-1", EcoliReference) };

        var result = MetagenomicsAnalyzer
            .ClassifyReads(reads, emptyDb, taxonomy, K)
            .ToList();

        result.Should().HaveCount(1, "one classification per input read, in input order");
        var c = result[0];

        c.TaxonId.Should().Be(TaxonomyTree.RootId,
            "with no database entries no k-mer can hit, so the read is unclassified (root) — INV-04");
        c.RtlScore.Should().Be(0, "an empty classification tree has no scored RTL path");
        c.MatchedKmers.Should().Be(0, "C = 0 when nothing in the clade was hit");
        c.Confidence.Should().Be(0.0, "C/Q with C = 0 is 0 (and never divides by zero)");
        new[] { 100, 200, 20, 21, 10 }.Should().NotContain(c.TaxonId,
            "an empty database must never fabricate a real-taxon assignment");
    }

    // ───────────────────────────────────────────────────────────────────────
    // Fuzz target: READ WITH NO MATCH (MC).
    // A valid DNA read whose k-mers are absent from a populated database ⇒
    // Unclassified (root), but Q (TotalKmers) is still counted.
    // — Taxonomic_Classification.md §"Edge cases" (No database matches).
    // ───────────────────────────────────────────────────────────────────────
    [Test]
    public void ClassifyReads_ValidReadNoMatch_AssignsRootButStillCountsQ()
    {
        var taxonomy = BuildTaxonomy();
        var db = BuildEcoliDatabase(taxonomy);

        // A run of 'G'/'T' whose canonical 4-mers (GGGG, GGGT, …) are not present
        // in the ACGT-repeat E.coli reference — valid DNA, but no hits.
        var reads = new[] { ("no-match", "GGGGGGGGTTTTTTTTGGGGGGGGTTTTTTTT") };

        var c = MetagenomicsAnalyzer
            .ClassifyReads(reads, db, taxonomy, K)
            .Single();

        c.TaxonId.Should().Be(TaxonomyTree.RootId,
            "a read with zero database hits is unclassified (root)");
        c.MatchedKmers.Should().Be(0, "no k-mer mapped to any taxon");
        c.TotalKmers.Should().BeGreaterThan(0,
            "Q still counts the non-ambiguous k-mers that were queried, even though none hit");
        c.Confidence.Should().Be(0.0, "C/Q = 0/Q = 0");
    }

    // ───────────────────────────────────────────────────────────────────────
    // Fuzz target: NON-DNA READ (MC).
    // Every k-mer contains a non-A/C/G/T symbol ⇒ all k-mers are skipped (not
    // counted in Q, never matched) ⇒ Q = 0, Unclassified (root). Handled, not
    // crashed, not rejected. — Taxonomic_Classification.md §"Edge cases".
    // ───────────────────────────────────────────────────────────────────────
    [Test]
    public void ClassifyReads_NonDnaRead_IsHandledAndUnclassified()
    {
        var taxonomy = BuildTaxonomy();
        var db = BuildEcoliDatabase(taxonomy);

        // Mixed garbage long enough to form k-mers, but every length-4 window
        // contains at least one non-ACGT symbol (IUPAC N/R, digits, gap, unicode).
        var reads = new[]
        {
            ("iupac",   "NNNNRRRRYYYYNNNN"),
            ("digits",  "1234567890123456"),
            ("gaps",    "----------------"),
            ("unicode", "ΑΒΓΔαβγδ★☢  ΑΒΓΔ"),
            ("lower-n", "nnnnnnnnnnnnnnnn"),
        };

        var results = MetagenomicsAnalyzer
            .ClassifyReads(reads, db, taxonomy, K)
            .ToList();

        results.Should().HaveCount(5);
        results.Should().OnlyContain(c => c.TaxonId == TaxonomyTree.RootId,
            "a read with no non-ambiguous k-mers cannot hit anything → unclassified (root)");
        results.Should().OnlyContain(c => c.TotalKmers == 0,
            "ambiguous k-mers are skipped and never counted in Q");
        results.Should().OnlyContain(c => c.Confidence == 0.0,
            "Q = 0 ⇒ Confidence is forced to 0.0 with no division by zero");
    }

    // ───────────────────────────────────────────────────────────────────────
    // Fuzz target: EXTREMELY SHORT READ (BE) — shorter than k.
    // No length-k window exists ⇒ the `sequence.Length >= k` guard short-circuits
    // the entire k-mer loop ⇒ Q = 0, Unclassified (root). KEY div-by-zero /
    // empty-window boundary. — Taxonomic_Classification.md §"Edge cases".
    // ───────────────────────────────────────────────────────────────────────
    [Test]
    public void ClassifyReads_ReadShorterThanK_NeverCrashesAndIsUnclassified()
    {
        var taxonomy = BuildTaxonomy();
        var db = BuildEcoliDatabase(taxonomy);

        // Empty read and reads of length 1..k-1 — none can form a single k-mer.
        var reads = new[]
        {
            ("empty",  ""),
            ("len-1",  "A"),
            ("len-2",  "AC"),
            ("len-3",  "ACG"),       // k-1 = 3, still no 4-mer
        };

        var results = MetagenomicsAnalyzer
            .ClassifyReads(reads, db, taxonomy, K)
            .ToList();

        results.Should().HaveCount(4);
        results.Should().OnlyContain(c => c.TaxonId == TaxonomyTree.RootId,
            "a read shorter than k forms no k-mers → unclassified (root)");
        results.Should().OnlyContain(c => c.TotalKmers == 0,
            "no k-mer window ⇒ Q = 0");
        results.Should().OnlyContain(c => c.Confidence == 0.0,
            "Q = 0 ⇒ Confidence 0.0, never a divide-by-zero on the empty boundary");
        results.Should().OnlyContain(c => c.RtlScore == 0,
            "no scored RTL path on an empty classification tree");
    }

    // ───────────────────────────────────────────────────────────────────────
    // Fuzz target: argument validation (BE) — deferred-iterator guards.
    // null reads / database / taxonomy → ArgumentNullException; k ≤ 0 →
    // ArgumentOutOfRangeException. These are *documented, intentional* guards,
    // not crashes. Because ClassifyReads is a deferred iterator, the guard fires
    // on the first enumeration.
    // ───────────────────────────────────────────────────────────────────────
    [Test]
    public void ClassifyReads_NullOrInvalidArguments_ThrowDocumentedExceptions()
    {
        var taxonomy = BuildTaxonomy();
        var db = BuildEcoliDatabase(taxonomy);
        var reads = new[] { ("r", EcoliReference) };

        // Deferred iterator → materialize with ToList() to trigger the guards.
        Action nullReads = () => MetagenomicsAnalyzer
            .ClassifyReads(null!, db, taxonomy, K).ToList();
        nullReads.Should().Throw<ArgumentNullException>()
            .WithParameterName("reads");

        Action nullDb = () => MetagenomicsAnalyzer
            .ClassifyReads(reads, null!, taxonomy, K).ToList();
        nullDb.Should().Throw<ArgumentNullException>()
            .WithParameterName("kmerDatabase");

        Action nullTaxonomy = () => MetagenomicsAnalyzer
            .ClassifyReads(reads, db, null!, K).ToList();
        nullTaxonomy.Should().Throw<ArgumentNullException>()
            .WithParameterName("taxonomy");

        foreach (int badK in new[] { 0, -1, int.MinValue })
        {
            Action invalidK = () => MetagenomicsAnalyzer
                .ClassifyReads(reads, db, taxonomy, badK).ToList();
            invalidK.Should().Throw<ArgumentOutOfRangeException>()
                .WithParameterName("k", $"k = {badK} is not positive");
        }
    }

    // ───────────────────────────────────────────────────────────────────────
    // Fuzz target: random malformed batch (MC, BE) under a time budget.
    // A deterministic, locally-seeded generator produces a mixed batch of empty,
    // short, non-DNA, and partially-valid reads. The classifier must process the
    // whole batch without crashing or hanging, and every result must be a
    // well-formed, in-tree taxon (root or a real node) — never an out-of-tree id
    // and never a divide-by-zero confidence.
    // ───────────────────────────────────────────────────────────────────────
    [Test]
    [CancelAfter(30000)]
    public void ClassifyReads_RandomMalformedBatch_NeverCrashesAndStaysWellFormed()
    {
        var taxonomy = BuildTaxonomy();
        var db = BuildEcoliDatabase(taxonomy);

        var rng = new Random(20260620); // locally fixed seed — deterministic
        const string alphabet = "ACGTacgtNRYWnX-09 \t★Α ";

        var reads = new List<(string, string)>();
        for (int i = 0; i < 200; i++)
        {
            int len = rng.Next(0, 24); // includes 0, < k, and ≥ k lengths
            var chars = new char[len];
            for (int j = 0; j < len; j++)
                chars[j] = alphabet[rng.Next(alphabet.Length)];
            reads.Add(($"read-{i}", new string(chars)));
        }

        var results = MetagenomicsAnalyzer
            .ClassifyReads(reads, db, taxonomy, K)
            .ToList();

        results.Should().HaveCount(200, "one classification per input read");
        results.Should().OnlyContain(c => taxonomy.Contains(c.TaxonId),
            "every assigned taxon must exist in the taxonomy tree (root or a real node)");
        results.Should().OnlyContain(c => c.TotalKmers >= 0 && c.MatchedKmers >= 0,
            "Q and C are non-negative counts");
        results.Should().OnlyContain(c => c.MatchedKmers <= c.TotalKmers,
            "C (clade hits) can never exceed Q (non-ambiguous k-mers queried)");
        results.Should().OnlyContain(c => c.Confidence >= 0.0 && c.Confidence <= 1.0,
            "C/Q confidence stays in [0, 1] with no divide-by-zero");
        // Reads with no successful query must be root, never a fabricated taxon.
        results.Where(c => c.TotalKmers == 0)
            .Should().OnlyContain(c => c.TaxonId == TaxonomyTree.RootId,
                "a read that queried zero k-mers cannot be assigned a real taxon");
    }

    // ───────────────────────────────────────────────────────────────────────
    // Positive sanity: a read that genuinely contains a taxon's k-mers IS
    // classified to that taxon (not Unclassified). Guards against a degenerate
    // "always root" classifier that would pass every boundary test above.
    // ───────────────────────────────────────────────────────────────────────
    [Test]
    public void ClassifyReads_ReadMatchingDatabaseTaxon_IsClassifiedToThatTaxon()
    {
        var taxonomy = BuildTaxonomy();
        var db = BuildEcoliDatabase(taxonomy);

        // The E.coli reference itself: all its canonical k-mers are owned by 100.
        var reads = new[] { ("ecoli-read", EcoliReference) };

        var c = MetagenomicsAnalyzer
            .ClassifyReads(reads, db, taxonomy, K)
            .Single();

        c.TaxonId.Should().Be(100,
            "a read built entirely from E.coli(100) k-mers must be classified to E.coli, not root");
        c.TaxonName.Should().Be("Escherichia coli");
        c.Species.Should().Be("Escherichia coli", "the species rank is read off the assigned lineage");
        c.TotalKmers.Should().BeGreaterThan(0, "Q counts the non-ambiguous k-mers queried");
        c.MatchedKmers.Should().BeGreaterThan(0, "C counts the k-mers supporting the assigned clade");
        c.RtlScore.Should().BeGreaterThan(0, "the winning root-to-leaf path has positive weight");
        c.Confidence.Should().BeInRange(0.0, 1.0).And.BeGreaterThan(0.0,
            "a genuine match yields positive C/Q confidence");
    }

    #endregion

    // ════════════════════════════════════════════════════════════════════════
    //
    //  META-PROF-001 — community / taxonomic profiling (Metagenomics).
    //  Checklist: docs/checklists/03_FUZZING.md, row 54.
    //
    //  What this unit profiles
    //  ───────────────────────
    //  MetagenomicsAnalyzer.GenerateTaxonomicProfile aggregates per-read
    //  TaxonomicClassification records into a community-level abundance summary.
    //  For a taxon i with count c_i, relative abundance is c_i / Σ_j c_j, where
    //  the denominator is the number of CLASSIFIED reads (kingdom not
    //  "Unclassified"/empty). It then computes species-level Shannon (H = −Σ pᵢ
    //  ln pᵢ) and Simpson (D = Σ pᵢ²) diversity.
    //    — docs/algorithms/Metagenomics/Taxonomic_Profile.md §2.2, §4.
    //
    //  Fuzz strategy for THIS unit: BE = Boundary Exploitation.
    //  Fuzz targets (checklist row 54): 0 reads, single read, all same taxon.
    //
    //  Boundary contract pinned by these tests (Taxonomic_Profile.md §3.3, §6.1):
    //    • 0 reads → TotalReads = 0, ClassifiedReads = 0, all abundance maps
    //      EMPTY, Shannon = Simpson = 0. KEY: no DivideByZero — abundance uses
    //      `total = classifiedReads > 0 ? classifiedReads : 1`, and the count
    //      dictionaries are empty so no map entry is even produced (§4.2).
    //    • single read → its kingdom/phylum/genus/species each appear at
    //      abundance 1.0 (a 1-category distribution); KingdomAbundance sums to 1
    //      (INV-03); Shannon = 0, Simpson = 1 (§6.1, single species).
    //    • all same taxon → exactly one entry per stored rank, all at 1.0; the
    //      KingdomAbundance map sums to 1; Shannon = 0, Simpson = 1.
    //    • KEY invariant INV-03: KingdomAbundance sums to 1 whenever
    //      ClassifiedReads > 0. Lower ranks sum to 1 only when every retained read
    //      has that rank populated (empty rank strings are skipped, §6.1).
    //    • all-unclassified reads → identical to 0 reads for the maps:
    //      ClassifiedReads = 0, empty maps, zero diversity (§6.1).
    //    • abundances are always in [0, 1].
    //
    //  Positive sanity: a known mixed set of classified reads gives the expected
    //  abundances (e.g. 3 E.coli + 1 S.enterica ⇒ species {E.coli:0.75,
    //  S.enterica:0.25}) summing to 1 — so a passing "no crash" result cannot be
    //  a profiler that just returns empty maps for everything.
    //
    //  Determinism: all classifications are hand-built or generated from a
    //  LOCALLY fixed-seed `new Random(seed)`. No shared static Rng.
    // ════════════════════════════════════════════════════════════════════════

    #region META-PROF-001 — community profiling

    // A minimal TaxonomicClassification carrying only the fields the profiler
    // reads (kingdom/phylum/genus/species + the unclassified filter on kingdom).
    // The taxonomy lineage strings mirror the META-CLASS-001 hand-built tree.
    private static MetagenomicsAnalyzer.TaxonomicClassification ClassifiedRead(
        string readId,
        string kingdom,
        string phylum,
        string genus,
        string species) =>
        new(
            ReadId: readId,
            TaxonId: 0,
            TaxonName: species,
            Rank: "species",
            RtlScore: 1,
            Confidence: 1.0,
            MatchedKmers: 1,
            TotalKmers: 1,
            Kingdom: kingdom,
            Phylum: phylum,
            Class: string.Empty,
            Order: string.Empty,
            Family: string.Empty,
            Genus: genus,
            Species: species);

    private static MetagenomicsAnalyzer.TaxonomicClassification Ecoli(string readId) =>
        ClassifiedRead(readId, "Bacteria", "Proteobacteria", "Escherichia", "Escherichia coli");

    private static MetagenomicsAnalyzer.TaxonomicClassification Senterica(string readId) =>
        ClassifiedRead(readId, "Bacteria", "Proteobacteria", "Salmonella", "Salmonella enterica");

    private static MetagenomicsAnalyzer.TaxonomicClassification Unclassified(string readId) =>
        ClassifiedRead(readId, "Unclassified", string.Empty, string.Empty, string.Empty);

    // ───────────────────────────────────────────────────────────────────────
    // Fuzz target: 0 READS (BE) — the div-by-zero boundary.
    // No classifications ⇒ TotalReads = 0, ClassifiedReads = 0, every abundance
    // map empty, zero diversity, and — KEY — no DivideByZeroException normalizing
    // abundance over 0 reads. — Taxonomic_Profile.md §3.3, §4.2, §6.1.
    // ───────────────────────────────────────────────────────────────────────
    [Test]
    public void GenerateTaxonomicProfile_ZeroReads_EmptyProfileNoDivideByZero()
    {
        var empty = Array.Empty<MetagenomicsAnalyzer.TaxonomicClassification>();

        MetagenomicsAnalyzer.TaxonomicProfile profile = default;
        Action act = () => profile = MetagenomicsAnalyzer.GenerateTaxonomicProfile(empty);

        act.Should().NotThrow(
            "normalizing abundance over 0 reads must never divide by zero — §4.2");

        profile.TotalReads.Should().Be(0, "no input classification records");
        profile.ClassifiedReads.Should().Be(0, "nothing survives the unclassified filter");
        profile.KingdomAbundance.Should().BeEmpty("no counts ⇒ no abundance entries");
        profile.PhylumAbundance.Should().BeEmpty();
        profile.GenusAbundance.Should().BeEmpty();
        profile.SpeciesAbundance.Should().BeEmpty();
        profile.ShannonDiversity.Should().Be(0.0, "an empty species distribution has zero diversity");
        profile.SimpsonDiversity.Should().Be(0.0);
    }

    // ───────────────────────────────────────────────────────────────────────
    // Fuzz target: SINGLE READ (BE).
    // One classified read ⇒ its kingdom/phylum/genus/species each appear at
    // abundance 1.0 (a single-category distribution). KingdomAbundance sums to 1
    // (INV-03). Diversity of one species: Shannon = 0, Simpson = 1. — §6.1.
    // ───────────────────────────────────────────────────────────────────────
    [Test]
    public void GenerateTaxonomicProfile_SingleRead_ThatTaxonAtAbundanceOne()
    {
        var reads = new[] { Ecoli("read-1") };

        var profile = MetagenomicsAnalyzer.GenerateTaxonomicProfile(reads);

        profile.TotalReads.Should().Be(1);
        profile.ClassifiedReads.Should().Be(1, "the single read is classified (kingdom = Bacteria)");

        profile.KingdomAbundance.Should().ContainKey("Bacteria")
            .WhoseValue.Should().Be(1.0, "the only read defines 100% of the kingdom abundance");
        profile.PhylumAbundance["Proteobacteria"].Should().Be(1.0);
        profile.GenusAbundance["Escherichia"].Should().Be(1.0);
        profile.SpeciesAbundance["Escherichia coli"].Should().Be(1.0);

        profile.KingdomAbundance.Values.Sum().Should().BeApproximately(1.0, 1e-12,
            "KingdomAbundance sums to 1 when ClassifiedReads > 0 — INV-03");
        profile.SpeciesAbundance.Values.Sum().Should().BeApproximately(1.0, 1e-12);

        profile.ShannonDiversity.Should().Be(0.0, "a single species has zero Shannon diversity");
        profile.SimpsonDiversity.Should().Be(1.0, "a single species has Simpson concentration 1");
    }

    // ───────────────────────────────────────────────────────────────────────
    // Fuzz target: ALL SAME TAXON (BE).
    // Many reads, all the same species ⇒ exactly one entry per stored rank, each
    // at abundance 1.0; maps sum to 1; Shannon = 0, Simpson = 1. — §6.1.
    // ───────────────────────────────────────────────────────────────────────
    [Test]
    public void GenerateTaxonomicProfile_AllSameTaxon_OneTaxonAtOneAllOthersAbsent()
    {
        var reads = Enumerable.Range(0, 50).Select(i => Ecoli($"read-{i}")).ToArray();

        var profile = MetagenomicsAnalyzer.GenerateTaxonomicProfile(reads);

        profile.TotalReads.Should().Be(50);
        profile.ClassifiedReads.Should().Be(50);

        profile.KingdomAbundance.Should().HaveCount(1).And.ContainKey("Bacteria");
        profile.SpeciesAbundance.Should().HaveCount(1, "all reads collapse to one species");
        profile.SpeciesAbundance["Escherichia coli"].Should().Be(1.0);
        profile.GenusAbundance["Escherichia"].Should().Be(1.0);

        // No other real taxon was fabricated.
        profile.SpeciesAbundance.Should().NotContainKey("Salmonella enterica");

        profile.KingdomAbundance.Values.Sum().Should().BeApproximately(1.0, 1e-12,
            "all mass concentrated on one taxon still sums to 1 — INV-03");

        profile.ShannonDiversity.Should().Be(0.0);
        profile.SimpsonDiversity.Should().Be(1.0);
    }

    // ───────────────────────────────────────────────────────────────────────
    // Fuzz target: ALL UNCLASSIFIED (BE).
    // Every read has kingdom "Unclassified" ⇒ all removed by the filter ⇒
    // ClassifiedReads = 0, empty maps, zero diversity — but TotalReads still
    // counts the inputs. No divide-by-zero. — Taxonomic_Profile.md §6.1.
    // ───────────────────────────────────────────────────────────────────────
    [Test]
    public void GenerateTaxonomicProfile_AllUnclassified_EmptyMapsButCountsTotalReads()
    {
        var reads = Enumerable.Range(0, 7).Select(i => Unclassified($"u-{i}")).ToArray();

        MetagenomicsAnalyzer.TaxonomicProfile profile = default;
        Action act = () => profile = MetagenomicsAnalyzer.GenerateTaxonomicProfile(reads);
        act.Should().NotThrow("0 classified reads must not divide by zero");

        profile.TotalReads.Should().Be(7, "TotalReads counts inputs before the unclassified filter");
        profile.ClassifiedReads.Should().Be(0, "every read is filtered out as unclassified");
        profile.KingdomAbundance.Should().BeEmpty();
        profile.SpeciesAbundance.Should().BeEmpty();
        profile.ShannonDiversity.Should().Be(0.0);
        profile.SimpsonDiversity.Should().Be(0.0);
    }

    // ───────────────────────────────────────────────────────────────────────
    // Fuzz target: random mixed batch (BE) — the KEY sum-to-1 invariant under a
    // deterministic, locally-seeded mix of two real taxa plus unclassified
    // noise. The profile must always be well-formed: every abundance in [0, 1],
    // KingdomAbundance summing to 1 whenever any read is classified, and never a
    // divide-by-zero. — Taxonomic_Profile.md INV-03.
    // ───────────────────────────────────────────────────────────────────────
    [Test]
    public void GenerateTaxonomicProfile_RandomMixedBatch_AbundancesInRangeAndSumToOne()
    {
        var rng = new Random(20260620); // locally fixed seed — deterministic

        var reads = new List<MetagenomicsAnalyzer.TaxonomicClassification>();
        int classifiedCount = 0;
        for (int i = 0; i < 300; i++)
        {
            int roll = rng.Next(3);
            if (roll == 0) { reads.Add(Ecoli($"r-{i}")); classifiedCount++; }
            else if (roll == 1) { reads.Add(Senterica($"r-{i}")); classifiedCount++; }
            else reads.Add(Unclassified($"r-{i}"));
        }

        var profile = MetagenomicsAnalyzer.GenerateTaxonomicProfile(reads);

        profile.TotalReads.Should().Be(300);
        profile.ClassifiedReads.Should().Be(classifiedCount);

        profile.KingdomAbundance.Values.Should().OnlyContain(v => v >= 0.0 && v <= 1.0,
            "every relative abundance lies in [0, 1]");
        profile.SpeciesAbundance.Values.Should().OnlyContain(v => v >= 0.0 && v <= 1.0);

        if (classifiedCount > 0)
        {
            profile.KingdomAbundance.Values.Sum().Should().BeApproximately(1.0, 1e-9,
                "KingdomAbundance sums to 1 whenever ClassifiedReads > 0 — INV-03 (KEY)");
            // Every retained read has phylum + genus + species populated here, so
            // those rank maps are also complete and sum to 1.
            profile.SpeciesAbundance.Values.Sum().Should().BeApproximately(1.0, 1e-9);
        }
    }

    // ───────────────────────────────────────────────────────────────────────
    // Positive sanity: a known mixed set yields the EXACT expected abundances,
    // in [0, 1] and summing to 1. Guards against a degenerate profiler that
    // returns empty maps (which would pass every boundary test above).
    //   3 × E.coli + 1 × S.enterica  ⇒  species {E.coli: 0.75, S.enterica: 0.25}.
    // ───────────────────────────────────────────────────────────────────────
    [Test]
    public void GenerateTaxonomicProfile_KnownMix_ExactExpectedAbundancesSummingToOne()
    {
        var reads = new[]
        {
            Ecoli("e1"), Ecoli("e2"), Ecoli("e3"),
            Senterica("s1"),
        };

        var profile = MetagenomicsAnalyzer.GenerateTaxonomicProfile(reads);

        profile.TotalReads.Should().Be(4);
        profile.ClassifiedReads.Should().Be(4);

        profile.SpeciesAbundance["Escherichia coli"].Should().BeApproximately(0.75, 1e-12,
            "3 of 4 classified reads are E.coli");
        profile.SpeciesAbundance["Salmonella enterica"].Should().BeApproximately(0.25, 1e-12,
            "1 of 4 classified reads is S.enterica");
        profile.SpeciesAbundance.Values.Sum().Should().BeApproximately(1.0, 1e-12,
            "the two species partition all classified reads");

        // Both share the same kingdom/phylum ⇒ those collapse to a single 1.0 entry.
        profile.KingdomAbundance["Bacteria"].Should().Be(1.0);
        profile.PhylumAbundance["Proteobacteria"].Should().Be(1.0);
        profile.GenusAbundance["Escherichia"].Should().BeApproximately(0.75, 1e-12);
        profile.GenusAbundance["Salmonella"].Should().BeApproximately(0.25, 1e-12);

        // A genuine two-species mix has positive Shannon diversity.
        profile.ShannonDiversity.Should().BeGreaterThan(0.0,
            "a two-species community is not a degenerate single-category distribution");
        profile.SimpsonDiversity.Should().BeInRange(0.0, 1.0)
            .And.BeLessThan(1.0, "two species ⇒ Simpson concentration below 1");
    }

    #endregion

    // ════════════════════════════════════════════════════════════════════════
    //
    //  META-ALPHA-001 — alpha diversity (Metagenomics).
    //  Checklist: docs/checklists/03_FUZZING.md, row 55. Strategy: BE.
    //  Entry point: MetagenomicsAnalyzer.CalculateAlphaDiversity
    //               (IReadOnlyDictionary<string, double> abundances).
    //  Contract pinned below: Alpha_Diversity.md §2.2, §2.4, §6.1.
    //
    //  Determinism: every input is hand-built or generated from a LOCALLY
    //  fixed-seed `new Random(seed)`. No shared static Rng.
    // ════════════════════════════════════════════════════════════════════════

    #region META-ALPHA-001 — alpha diversity

    // ───────────────────────────────────────────────────────────────────────
    // Fuzz target: 0 SPECIES (BE) — empty / null map.
    // No taxa ⇒ an all-zero AlphaDiversity record. KEY: the `Count == 0` guard
    // short-circuits before any normalization, so there is NO DivideByZero and
    // NO log(0) → NaN. — Alpha_Diversity.md §3.3, §6.1 (Empty input).
    // ───────────────────────────────────────────────────────────────────────
    [Test]
    public void CalculateAlphaDiversity_ZeroSpecies_AllZeroNoDivByZeroNoNaN()
    {
        foreach (var abundances in new IReadOnlyDictionary<string, double>?[]
                 {
                     new Dictionary<string, double>(),      // empty
                     null,                                   // null map
                 })
        {
            MetagenomicsAnalyzer.AlphaDiversity result = default;
            Action act = () => result = MetagenomicsAnalyzer.CalculateAlphaDiversity(abundances!);

            act.Should().NotThrow("0 species must not divide by zero or evaluate log(0)");

            result.ObservedSpecies.Should().Be(0, "no taxa ⇒ zero observed species");
            result.ShannonIndex.Should().Be(0.0, "empty distribution has zero Shannon entropy");
            result.ShannonIndex.Should().NotBe(double.NaN, "no p·ln(p) term is ever evaluated");
            result.SimpsonIndex.Should().Be(0.0);
            result.InverseSimpson.Should().Be(0.0, "no reciprocal of a zero Simpson concentration");
            result.Chao1Estimate.Should().Be(0.0);
            result.PielouEvenness.Should().Be(0.0);
        }
    }

    // ───────────────────────────────────────────────────────────────────────
    // Fuzz target: SINGLE SAMPLE WITH 0 (BE) — every abundance ≤ 0.
    // All-zero (and negative) counts ⇒ the `Where(v => v > 0)` filter empties
    // the value list ⇒ observedSpecies == 0 ⇒ all-zero record. KEY: the p·ln(p)
    // term is never reached for p = 0, so no log(0) NaN appears. — §3.3, §6.1.
    // ───────────────────────────────────────────────────────────────────────
    [Test]
    public void CalculateAlphaDiversity_AllNonPositiveCounts_DefinedAllZeroNoNaN()
    {
        var allZeroOrNegative = new Dictionary<string, double>
        {
            ["taxonA"] = 0.0,
            ["taxonB"] = 0.0,
            ["taxonC"] = -3.0,   // negative is also filtered out, never log()'d
            ["taxonD"] = 0.0,
        };

        MetagenomicsAnalyzer.AlphaDiversity result = default;
        Action act = () => result = MetagenomicsAnalyzer.CalculateAlphaDiversity(allZeroOrNegative);

        act.Should().NotThrow("an all-zero sample is a defined boundary, not a crash");

        result.ObservedSpecies.Should().Be(0,
            "no taxon has strictly positive abundance ⇒ zero observed species");
        result.ShannonIndex.Should().Be(0.0).And.NotBe(double.NaN,
            "p = 0 never reaches p·ln(p): no log(0) NaN");
        result.SimpsonIndex.Should().Be(0.0).And.NotBe(double.NaN);
        result.InverseSimpson.Should().Be(0.0);
        result.Chao1Estimate.Should().Be(0.0);
        result.PielouEvenness.Should().Be(0.0);
    }

    // ───────────────────────────────────────────────────────────────────────
    // Fuzz target: 1 SPECIES (BE).
    // One positive abundance ⇒ Shannon = 0 (INV-02), Simpson concentration = 1
    // (INV-03), InverseSimpson = 1, ObservedSpecies = 1, Pielou = 0 (undefined
    // for S ≤ 1, INV-05). The lone abundance value is irrelevant after internal
    // normalization. — Alpha_Diversity.md §6.1 (Single species).
    // ───────────────────────────────────────────────────────────────────────
    [Test]
    public void CalculateAlphaDiversity_SingleSpecies_ShannonZeroSimpsonOne()
    {
        // Try several positive magnitudes: normalization makes p = 1 regardless.
        foreach (double abundance in new[] { 1.0, 7.0, 1234.5, 0.001 })
        {
            var single = new Dictionary<string, double> { ["only-taxon"] = abundance };

            var result = MetagenomicsAnalyzer.CalculateAlphaDiversity(single);

            result.ObservedSpecies.Should().Be(1, "exactly one positive taxon");
            result.ShannonIndex.Should().Be(0.0,
                $"Shannon of a single species is 0 regardless of abundance ({abundance}) — INV-02");
            result.SimpsonIndex.Should().BeApproximately(1.0, 1e-12,
                "Σpᵢ² with a single pᵢ = 1 is 1 — INV-03 (concentration, not 1−Σpᵢ²)");
            result.InverseSimpson.Should().BeApproximately(1.0, 1e-12,
                "1 / Simpson = 1 / 1 = 1 — INV-04");
            result.PielouEvenness.Should().Be(0.0,
                "Pielou's evenness is 0 when S ≤ 1 (ln(1) = 0 boundary) — INV-05");
        }
    }

    // ───────────────────────────────────────────────────────────────────────
    // Fuzz target: ALL EQUAL ABUNDANCE over S species (BE) — the maximal-Shannon
    // boundary. Uniform abundances MAXIMIZE Shannon at fixed richness:
    //   Shannon = ln(S)  (KEY),
    //   Simpson concentration λ = Σ(1/S)² = S·(1/S²) = 1/S  (NOT 1 − 1/S),
    //   InverseSimpson = S, PielouEvenness = 1.
    // Pinned for S = 2, 3, 4. — Alpha_Diversity.md §2.2, §6.1.
    // ───────────────────────────────────────────────────────────────────────
    [Test]
    public void CalculateAlphaDiversity_AllEqualAbundance_ShannonIsLnS_Maximal()
    {
        foreach (int s in new[] { 2, 3, 4 })
        {
            // Equal counts of 5 each — normalization gives p_i = 1/S for every taxon.
            var uniform = Enumerable.Range(0, s)
                .ToDictionary(i => $"taxon-{i}", _ => 5.0);

            var result = MetagenomicsAnalyzer.CalculateAlphaDiversity(uniform);

            result.ObservedSpecies.Should().Be(s);
            result.ShannonIndex.Should().BeApproximately(Math.Log(s), 1e-12,
                $"uniform abundance MAXIMIZES Shannon at ln(S) = ln({s}) — KEY");
            result.SimpsonIndex.Should().BeApproximately(1.0 / s, 1e-12,
                $"Σ(1/S)² = 1/S = 1/{s} (Simpson concentration, §2.2)");
            result.InverseSimpson.Should().BeApproximately(s, 1e-12,
                $"1 / (1/S) = S = {s} effective species — INV-04");
            result.PielouEvenness.Should().BeApproximately(1.0, 1e-12,
                "H / ln(S) = ln(S) / ln(S) = 1 at maximal evenness");
        }
    }

    // ───────────────────────────────────────────────────────────────────────
    // Fuzz target: random non-negative batch (BE) under a time budget.
    // A deterministic, locally-seeded generator builds maps mixing zeros,
    // negatives and positives over a varying number of taxa. The metric must
    // process every map without crashing or hanging, and every result must be
    // well-formed: finite (no NaN/∞), Shannon ≥ 0, Simpson in (0, 1] when any
    // taxon is positive (else all-zero), and Shannon ≤ ln(S_obs).
    // ───────────────────────────────────────────────────────────────────────
    [Test]
    [CancelAfter(30000)]
    public void CalculateAlphaDiversity_RandomNonNegativeBatch_AlwaysFiniteAndWellFormed()
    {
        var rng = new Random(20260620); // locally fixed seed — deterministic

        for (int iter = 0; iter < 500; iter++)
        {
            int taxa = rng.Next(0, 12); // includes 0 taxa
            var map = new Dictionary<string, double>();
            for (int t = 0; t < taxa; t++)
            {
                // A mix of zeros, negatives and positive magnitudes.
                double v = rng.Next(4) switch
                {
                    0 => 0.0,
                    1 => -rng.NextDouble() * 10.0,
                    _ => rng.NextDouble() * 1000.0,
                };
                map[$"t-{t}"] = v;
            }

            MetagenomicsAnalyzer.AlphaDiversity r = default;
            Action act = () => r = MetagenomicsAnalyzer.CalculateAlphaDiversity(map);
            act.Should().NotThrow($"iteration {iter} must not crash");

            int positives = map.Values.Count(v => v > 0);

            r.ObservedSpecies.Should().Be(positives,
                "observed species = count of strictly positive abundances — INV-01");

            double.IsNaN(r.ShannonIndex).Should().BeFalse("Shannon is never NaN (no log(0))");
            double.IsNaN(r.SimpsonIndex).Should().BeFalse();
            double.IsInfinity(r.ShannonIndex).Should().BeFalse();
            double.IsInfinity(r.InverseSimpson).Should().BeFalse(
                "InverseSimpson is 0 (not ∞) when Simpson is 0");

            r.ShannonIndex.Should().BeGreaterThanOrEqualTo(0.0, "Shannon entropy is non-negative");
            r.SimpsonIndex.Should().BeGreaterThanOrEqualTo(0.0);

            if (positives == 0)
            {
                r.ShannonIndex.Should().Be(0.0);
                r.SimpsonIndex.Should().Be(0.0);
                r.InverseSimpson.Should().Be(0.0);
                r.PielouEvenness.Should().Be(0.0);
            }
            else
            {
                r.SimpsonIndex.Should().BeInRange(0.0, 1.0 + 1e-12,
                    "Simpson concentration lies in (0, 1]");
                // Shannon is bounded above by ln(S_obs) (max at uniform abundance).
                double maxH = positives > 1 ? Math.Log(positives) : 0.0;
                r.ShannonIndex.Should().BeLessThanOrEqualTo(maxH + 1e-9,
                    "Shannon never exceeds ln(S_obs), its uniform-abundance maximum");
            }
        }
    }

    // ───────────────────────────────────────────────────────────────────────
    // Positive sanity: a known asymmetric two-species sample pins EXACT Shannon
    // and Simpson values against the natural-log formula, guarding against a
    // degenerate diversity function that returns a constant (0 or 1) everywhere.
    //   counts (3, 1) ⇒ p = (0.75, 0.25)
    //   Shannon = −(0.75 ln 0.75 + 0.25 ln 0.25) ≈ 0.5623351446188083
    //   Simpson = 0.75² + 0.25² = 0.625, InverseSimpson = 1/0.625 = 1.6.
    // ───────────────────────────────────────────────────────────────────────
    [Test]
    public void CalculateAlphaDiversity_KnownAsymmetricMix_ExactShannonAndSimpson()
    {
        var sample = new Dictionary<string, double>
        {
            ["Escherichia coli"]    = 3.0,
            ["Salmonella enterica"] = 1.0,
        };

        var result = MetagenomicsAnalyzer.CalculateAlphaDiversity(sample);

        result.ObservedSpecies.Should().Be(2);

        const double expectedShannon = 0.5623351446188083; // −Σ pᵢ ln pᵢ for (0.75, 0.25)
        result.ShannonIndex.Should().BeApproximately(expectedShannon, 1e-12,
            "exact natural-log Shannon for p = (0.75, 0.25)");
        result.ShannonIndex.Should().BeGreaterThan(0.0,
            "a genuine two-species mix is not a degenerate single-category distribution");

        result.SimpsonIndex.Should().BeApproximately(0.625, 1e-12,
            "0.75² + 0.25² = 0.625 (Simpson concentration)");
        result.InverseSimpson.Should().BeApproximately(1.6, 1e-12,
            "1 / 0.625 = 1.6 effective species");

        // Uniform-maximum sanity: H is strictly below the ln(2) achievable at even split.
        result.ShannonIndex.Should().BeLessThan(Math.Log(2),
            "an uneven 3:1 split has lower Shannon than the even-split maximum ln(2)");
    }

    #endregion

    // ════════════════════════════════════════════════════════════════════════
    //
    //  META-BETA-001 — beta diversity (Metagenomics).
    //  Checklist: docs/checklists/03_FUZZING.md, row 56. Strategy: BE.
    //  Entry point: MetagenomicsAnalyzer.CalculateBetaDiversity
    //      (string s1Name, IReadOnlyDictionary<string,double> s1,
    //       string s2Name, IReadOnlyDictionary<string,double> s2).
    //  Contract pinned below: Beta_Diversity.md §2, §3.3, §6.1.
    //
    //  Determinism: every input is hand-built or generated from a LOCALLY
    //  fixed-seed `new Random(seed)`. No shared static Rng.
    // ════════════════════════════════════════════════════════════════════════

    #region META-BETA-001 — beta diversity

    // ───────────────────────────────────────────────────────────────────────
    // Fuzz target: IDENTICAL SAMPLES (BE) — the d(a,a) = 0 boundary (KEY).
    // A sample compared with itself has zero dissimilarity in BOTH metrics:
    // shared abundance is half the total ⇒ Bray-Curtis = 0 (INV-BRAY-02), and
    // |A∩B| = |A∪B| ⇒ Jaccard = 0 (INV-JACCARD-02). Pinned for several profiles.
    // — Beta_Diversity.md §6.1 (Identical samples).
    // ───────────────────────────────────────────────────────────────────────
    [Test]
    public void CalculateBetaDiversity_IdenticalSamples_ZeroDissimilarity()
    {
        var profiles = new IReadOnlyDictionary<string, double>[]
        {
            new Dictionary<string, double> { ["E.coli"] = 10.0 },
            new Dictionary<string, double> { ["E.coli"] = 3.0, ["S.enterica"] = 7.0 },
            new Dictionary<string, double> { ["a"] = 1.0, ["b"] = 2.0, ["c"] = 3.0, ["d"] = 4.0 },
            new Dictionary<string, double> { ["x"] = 0.001, ["y"] = 999.0 },
        };

        foreach (var p in profiles)
        {
            var r = MetagenomicsAnalyzer.CalculateBetaDiversity("s", p, "s-copy", p);

            r.BrayCurtis.Should().BeApproximately(0.0, 1e-12,
                "Bray-Curtis of a sample with itself is 0 — INV-BRAY-02 (d(a,a) = 0)");
            r.JaccardDistance.Should().BeApproximately(0.0, 1e-12,
                "Jaccard of identical species sets is 0 — INV-JACCARD-02");
            r.SharedSpecies.Should().Be(p.Count(kv => kv.Value > 0),
                "every positive-abundance taxon is shared with the identical copy");
            r.UniqueToSample1.Should().Be(0, "nothing is unique when comparing a sample with itself");
            r.UniqueToSample2.Should().Be(0);
        }
    }

    // ───────────────────────────────────────────────────────────────────────
    // Fuzz target: EMPTY / ALL-ZERO SAMPLES (BE) — the KEY div-by-zero boundary.
    // An empty union (or all-zero abundances) has no positively-present taxon and
    // a total abundance of 0. Both metrics GUARD their denominators and return 0
    // with NO DivideByZeroException on the 0 denominator. — Beta_Diversity.md §3.3.
    // ───────────────────────────────────────────────────────────────────────
    [Test]
    public void CalculateBetaDiversity_EmptyAndAllZeroSamples_DefinedZeroNoDivByZero()
    {
        var empty = new Dictionary<string, double>();
        var allZero = new Dictionary<string, double> { ["a"] = 0.0, ["b"] = 0.0, ["c"] = -2.0 };
        var nonEmpty = new Dictionary<string, double> { ["a"] = 5.0, ["b"] = 1.0 };

        // empty vs empty, all-zero vs all-zero, empty vs all-zero — every pairing
        // of "no informative comparison" must return 0/0 without dividing by zero.
        var noInfoPairs = new[]
        {
            (empty, empty),
            (allZero, allZero),
            (empty, allZero),
            (allZero, empty),
        };

        foreach (var (a, b) in noInfoPairs)
        {
            MetagenomicsAnalyzer.BetaDiversity r = default;
            Action act = () => r = MetagenomicsAnalyzer.CalculateBetaDiversity("a", a, "b", b);

            act.Should().NotThrow(
                "an empty / all-zero union must guard both denominators — no DivideByZero (§3.3)");

            r.BrayCurtis.Should().Be(0.0,
                "Σtotal = 0 ⇒ Bray-Curtis guard returns 0, never NaN or a division by zero");
            r.JaccardDistance.Should().Be(0.0,
                "the union has no positive taxon ⇒ Jaccard guard returns 0");
            r.SharedSpecies.Should().Be(0);
            r.UniqueToSample1.Should().Be(0);
            r.UniqueToSample2.Should().Be(0);
        }

        // A populated sample vs an empty one: no shared positive taxon, but the
        // populated side has positive total abundance ⇒ MAXIMAL dissimilarity 1,
        // still with no divide-by-zero. — INV-BRAY-03 / INV-JACCARD-03.
        var rMax = MetagenomicsAnalyzer.CalculateBetaDiversity("full", nonEmpty, "empty", empty);
        rMax.BrayCurtis.Should().BeApproximately(1.0, 1e-12,
            "no shared abundance but positive total ⇒ Bray-Curtis = 1 — INV-BRAY-03");
        rMax.JaccardDistance.Should().BeApproximately(1.0, 1e-12,
            "disjoint non-empty sets ⇒ Jaccard = 1 — INV-JACCARD-03");
        rMax.UniqueToSample1.Should().Be(2, "both taxa are unique to the populated sample");
        rMax.UniqueToSample2.Should().Be(0);
        rMax.SharedSpecies.Should().Be(0);
    }

    // ───────────────────────────────────────────────────────────────────────
    // Fuzz target: SINGLE-SPECIES SAMPLES (BE).
    //   • same species, equal abundance  → identical ⇒ both metrics 0;
    //   • same species, different abundance → Jaccard = 0 (same set) but
    //     Bray-Curtis = |a−b|/(a+b) > 0 (abundance-sensitive, §2.A);
    //   • two DIFFERENT single species → disjoint non-empty ⇒ both metrics 1.
    // — Beta_Diversity.md §6.1.
    // ───────────────────────────────────────────────────────────────────────
    [Test]
    public void CalculateBetaDiversity_SingleSpeciesSamples_DefinedAcrossTheBoundary()
    {
        // (a) Same single species, equal abundance → identical → 0 / 0.
        var sameEqualA = new Dictionary<string, double> { ["sp"] = 4.0 };
        var sameEqualB = new Dictionary<string, double> { ["sp"] = 4.0 };
        var rEqual = MetagenomicsAnalyzer.CalculateBetaDiversity("a", sameEqualA, "b", sameEqualB);
        rEqual.BrayCurtis.Should().BeApproximately(0.0, 1e-12, "identical single-species samples ⇒ 0");
        rEqual.JaccardDistance.Should().BeApproximately(0.0, 1e-12);
        rEqual.SharedSpecies.Should().Be(1);

        // (b) Same single species, different abundance → same SET (Jaccard 0) but
        //     Bray-Curtis = |3−1|/(3+1) = 0.5 (abundance-sensitive).
        var sameDiffA = new Dictionary<string, double> { ["sp"] = 3.0 };
        var sameDiffB = new Dictionary<string, double> { ["sp"] = 1.0 };
        var rDiff = MetagenomicsAnalyzer.CalculateBetaDiversity("a", sameDiffA, "b", sameDiffB);
        rDiff.JaccardDistance.Should().BeApproximately(0.0, 1e-12,
            "the presence/absence set is identical ⇒ Jaccard 0, even at different abundance");
        rDiff.BrayCurtis.Should().BeApproximately(0.5, 1e-12,
            "Bray-Curtis = 1 − 2·min(3,1)/(3+1) = 1 − 2·1/4 = 0.5 (abundance-sensitive, §2.A)");
        rDiff.SharedSpecies.Should().Be(1, "the single species is present in both");

        // (c) Two DIFFERENT single species → disjoint non-empty ⇒ maximal 1 / 1.
        var spX = new Dictionary<string, double> { ["X"] = 5.0 };
        var spY = new Dictionary<string, double> { ["Y"] = 5.0 };
        var rDisjoint = MetagenomicsAnalyzer.CalculateBetaDiversity("a", spX, "b", spY);
        rDisjoint.BrayCurtis.Should().BeApproximately(1.0, 1e-12,
            "no shared taxon, positive total ⇒ Bray-Curtis 1 — INV-BRAY-03");
        rDisjoint.JaccardDistance.Should().BeApproximately(1.0, 1e-12,
            "disjoint non-empty single-species sets ⇒ Jaccard 1 — INV-JACCARD-03");
        rDisjoint.SharedSpecies.Should().Be(0);
        rDisjoint.UniqueToSample1.Should().Be(1);
        rDisjoint.UniqueToSample2.Should().Be(1);
    }

    // ───────────────────────────────────────────────────────────────────────
    // Fuzz target: random non-negative pairs (BE) under a time budget.
    // A deterministic, locally-seeded generator builds pairs of maps mixing
    // zeros and positive abundances over overlapping taxon universes. Every
    // result must be well-formed: both metrics finite and in [0, 1], SYMMETRIC
    // (d(a,b) = d(b,a)), zero on the diagonal (d(a,a) = 0), and never a
    // divide-by-zero — regardless of how sparse / empty the inputs happen to be.
    // ───────────────────────────────────────────────────────────────────────
    [Test]
    [CancelAfter(30000)]
    public void CalculateBetaDiversity_RandomPairs_SymmetricInRangeZeroDiagonal()
    {
        var rng = new Random(20260620); // locally fixed seed — deterministic
        var taxa = new[] { "t0", "t1", "t2", "t3", "t4" };

        Dictionary<string, double> RandomSample()
        {
            var map = new Dictionary<string, double>();
            foreach (var t in taxa)
            {
                // ~40% absent, else a mix of zeros (absent) and positive magnitudes.
                if (rng.Next(5) == 0) continue;
                map[t] = rng.Next(3) switch
                {
                    0 => 0.0,                         // present key, but absent (≤ 0)
                    _ => rng.NextDouble() * 100.0,    // positive abundance
                };
            }
            return map;
        }

        for (int iter = 0; iter < 500; iter++)
        {
            var a = RandomSample();
            var b = RandomSample();

            MetagenomicsAnalyzer.BetaDiversity rab = default, rba = default, raa = default;
            Action act = () =>
            {
                rab = MetagenomicsAnalyzer.CalculateBetaDiversity("a", a, "b", b);
                rba = MetagenomicsAnalyzer.CalculateBetaDiversity("b", b, "a", a);
                raa = MetagenomicsAnalyzer.CalculateBetaDiversity("a", a, "a", a);
            };
            act.Should().NotThrow($"iteration {iter} must not crash or divide by zero");

            foreach (var d in new[] { rab.BrayCurtis, rab.JaccardDistance })
            {
                double.IsNaN(d).Should().BeFalse("dissimilarity is never NaN");
                double.IsInfinity(d).Should().BeFalse();
                d.Should().BeInRange(0.0, 1.0, "every dissimilarity lies in [0, 1]");
            }

            rab.BrayCurtis.Should().BeApproximately(rba.BrayCurtis, 1e-12,
                "Bray-Curtis is symmetric: d(a,b) = d(b,a)");
            rab.JaccardDistance.Should().BeApproximately(rba.JaccardDistance, 1e-12,
                "Jaccard is symmetric: d(a,b) = d(b,a)");

            raa.BrayCurtis.Should().BeApproximately(0.0, 1e-12, "d(a,a) = 0 — Bray-Curtis");
            raa.JaccardDistance.Should().BeApproximately(0.0, 1e-12, "d(a,a) = 0 — Jaccard");

            // SharedSpecies counts only taxa positive in BOTH samples.
            int expectedShared = taxa.Count(t =>
                a.GetValueOrDefault(t) > 0 && b.GetValueOrDefault(t) > 0);
            rab.SharedSpecies.Should().Be(expectedShared,
                "SharedSpecies = taxa with strictly positive abundance in both samples");
        }
    }

    // ───────────────────────────────────────────────────────────────────────
    // Positive sanity: a known asymmetric pair pins EXACT Bray-Curtis and Jaccard
    // against the documented formulas, guarding against a degenerate dissimilarity
    // that returns a constant (0 or 1) everywhere.
    //   s1 {A:3, B:1}  vs  s2 {A:1, C:3}
    //   Bray-Curtis = 1 − 2·(min(3,1)+min(1,0)+min(0,3)) / (4+4) = 1 − 2·1/8 = 0.75
    //   Jaccard     = 1 − shared/(shared+u1+u2) = 1 − 1/3 ≈ 0.6666…  (A shared)
    // ───────────────────────────────────────────────────────────────────────
    [Test]
    public void CalculateBetaDiversity_KnownAsymmetricPair_ExactBrayCurtisAndJaccard()
    {
        var s1 = new Dictionary<string, double> { ["A"] = 3.0, ["B"] = 1.0 };
        var s2 = new Dictionary<string, double> { ["A"] = 1.0, ["C"] = 3.0 };

        var r = MetagenomicsAnalyzer.CalculateBetaDiversity("s1", s1, "s2", s2);

        r.BrayCurtis.Should().BeApproximately(0.75, 1e-12,
            "1 − 2·Σmin/Σtotal = 1 − 2·1/8 = 0.75 (Σmin = 1, Σtotal = 8)");
        r.JaccardDistance.Should().BeApproximately(1.0 - 1.0 / 3.0, 1e-12,
            "1 − shared/(shared+u1+u2) = 1 − 1/3 (A shared; B, C each unique)");

        r.SharedSpecies.Should().Be(1, "only A is positive in both samples");
        r.UniqueToSample1.Should().Be(1, "B is unique to sample 1");
        r.UniqueToSample2.Should().Be(1, "C is unique to sample 2");

        // Non-degenerate: strictly between the identical (0) and disjoint (1) extremes.
        r.BrayCurtis.Should().BeInRange(0.0, 1.0).And.BeGreaterThan(0.0).And.BeLessThan(1.0);
        r.JaccardDistance.Should().BeInRange(0.0, 1.0).And.BeGreaterThan(0.0).And.BeLessThan(1.0);

        // Symmetry on the pinned pair: d(s1,s2) = d(s2,s1).
        var swapped = MetagenomicsAnalyzer.CalculateBetaDiversity("s2", s2, "s1", s1);
        swapped.BrayCurtis.Should().BeApproximately(r.BrayCurtis, 1e-12);
        swapped.JaccardDistance.Should().BeApproximately(r.JaccardDistance, 1e-12);

        r.UniFracDistance.Should().Be(0.0, "UniFrac is a hard-coded 0 placeholder — §5.4");
    }

    #endregion
}
