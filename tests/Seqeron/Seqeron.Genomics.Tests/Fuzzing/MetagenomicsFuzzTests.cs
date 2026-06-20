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
}
