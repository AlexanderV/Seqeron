using NUnit.Framework;
using Seqeron.Genomics.Metagenomics;
using System.Collections.Generic;
using System.Linq;

namespace Seqeron.Genomics.Tests.Unit.Metagenomics;

/// <summary>
/// Test suite for META-CLASS-001: Kraken-style taxonomic classification
/// (taxonomy tree + k-mer LCA database + per-read RTL classification).
///
/// Source: Wood DE, Salzberg SL. "Kraken: ultrafast metagenomic sequence classification
/// using exact alignments." Genome Biology 2014, 15:R46. doi:10.1186/gb-2014-15-3-r46.
/// Verbatim rules used here:
/// - DB build: a record is "a k-mer and the LCA of all organisms whose genomes contain that
///   k-mer"; when a k-mer's LCA is already set, store "the LCA of the stored value and the
///   current sequence's taxon".
/// - Classification: hit taxa and their ancestors form the classification tree; "Each node ...
///   is weighted with the number of k-mers ... that mapped to the taxon associated with that
///   node"; "each root-to-leaf (RTL) path ... is scored by ... the sum of all node weights along
///   the path. The maximum scoring RTL path ... is the classification path"; ties → "the LCA of
///   all those paths' leaves is selected". No hits → unclassified (root here).
///
/// All expected taxa below are hand-derived from these rules on a fixed hand-built taxonomy, and
/// every k-mer used is its own canonical form (kmer == min(kmer, revcomp)) at k = 4, so the
/// k-mer→taxon mapping is fully controlled and Q/C are countable by hand.
/// </summary>
[TestFixture]
public class MetagenomicsAnalyzer_TaxonomicClassification_Tests
{
    // ------------------------------------------------------------------
    // Hand-built taxonomy (NCBI-style ids):
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
    //                        ┌────────┴────────┐            │
    //                  E.coli(100)     E.fergusonii(101)  S.enterica(200)
    //                   [species]        [species]         [species]
    //
    // Hand-derived LCAs:
    //   Lca(100,101) = 20  (siblings → parent genus)
    //   Lca(100,200) = 10  (different genera → family)
    //   Lca(20,100)  = 20  (ancestor/descendant → ancestor)
    //   Lca(100,100) = 100 (self)
    //   Lca(100,5)   = 1   (disjoint top branches → root)
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
        new TaxonNode(101, "Escherichia fergusonii", "species", 20),
        new TaxonNode(200, "Salmonella enterica", "species", 21),
    });

    #region TaxonomyTree.Lca unit tests (hand-derived)

    [Test]
    [Description("LCA: siblings → their shared parent")]
    public void Lca_Siblings_ReturnsParent()
    {
        var t = BuildTaxonomy();
        // E.coli(100) and E.fergusonii(101) are siblings under Escherichia(20).
        Assert.That(t.Lca(100, 101), Is.EqualTo(20));
        // Order independent.
        Assert.That(t.Lca(101, 100), Is.EqualTo(20));
    }

    [Test]
    [Description("LCA: ancestor/descendant → the ancestor")]
    public void Lca_AncestorDescendant_ReturnsAncestor()
    {
        var t = BuildTaxonomy();
        Assert.That(t.Lca(20, 100), Is.EqualTo(20));
        Assert.That(t.Lca(100, 20), Is.EqualTo(20));
    }

    [Test]
    [Description("LCA: a node with itself → itself")]
    public void Lca_SameNode_ReturnsItself()
    {
        var t = BuildTaxonomy();
        Assert.That(t.Lca(100, 100), Is.EqualTo(100));
    }

    [Test]
    [Description("LCA: disjoint top-level branches → root")]
    public void Lca_DisjointBranches_ReturnsRoot()
    {
        var t = BuildTaxonomy();
        // Archaea(5) and E.coli(100) share only the root.
        Assert.That(t.Lca(5, 100), Is.EqualTo(1));
    }

    [Test]
    [Description("LCA: different genera resolve to their common family")]
    public void Lca_DifferentGenera_ReturnsFamily()
    {
        var t = BuildTaxonomy();
        // E.coli(100, genus Escherichia) vs S.enterica(200, genus Salmonella) → family(10).
        Assert.That(t.Lca(100, 200), Is.EqualTo(10));
    }

    [Test]
    [Description("LCA over a set folds pairwise: {100,101,200} → family(10)")]
    public void Lca_OfSet_FoldsPairwise()
    {
        var t = BuildTaxonomy();
        // Lca(100,101)=20, then Lca(20,200)=10.
        Assert.That(t.Lca(new[] { 100, 101, 200 }), Is.EqualTo(10));
    }

    [Test]
    [Description("Tree helpers: path-to-root, depth, ancestry")]
    public void Tree_PathAndAncestry_AreCorrect()
    {
        var t = BuildTaxonomy();
        Assert.Multiple(() =>
        {
            Assert.That(t.GetPathToRoot(100), Is.EqualTo(new[] { 100, 20, 10, 4, 3, 2, 1 }));
            Assert.That(t.GetDepth(1), Is.EqualTo(0));
            Assert.That(t.GetDepth(100), Is.EqualTo(6));
            Assert.That(t.IsAncestorOf(20, 100), Is.True);
            Assert.That(t.IsAncestorOf(100, 20), Is.False);
            Assert.That(t.IsAncestorOf(100, 100), Is.True, "a taxon is its own ancestor");
        });
    }

    #endregion

    #region BuildKmerDatabase tests

    [Test]
    [Description("Empty reference set → empty database")]
    public void BuildKmerDatabase_EmptyInput_ReturnsEmpty()
    {
        var t = BuildTaxonomy();
        var db = MetagenomicsAnalyzer.BuildKmerDatabase(
            new List<(int, string)>(), t, k: 4);
        Assert.That(db, Is.Empty);
    }

    [Test]
    [Description("Reference shorter than k → no k-mers")]
    public void BuildKmerDatabase_SequenceShorterThanK_ReturnsEmpty()
    {
        var t = BuildTaxonomy();
        var db = MetagenomicsAnalyzer.BuildKmerDatabase(
            new[] { (100, "AGC") }, t, k: 4);
        Assert.That(db, Is.Empty);
    }

    [Test]
    [Description("Single reference: each canonical k-mer maps to that taxon")]
    public void BuildKmerDatabase_SingleReference_MapsToTaxon()
    {
        var t = BuildTaxonomy();
        // "AAAACAA" (k=4): windows AAAA, AAAC, AACA, ACAA — all distinct, all self-canonical.
        var db = MetagenomicsAnalyzer.BuildKmerDatabase(
            new[] { (100, "AAAACAA") }, t, k: 4);

        Assert.Multiple(() =>
        {
            Assert.That(db.Count, Is.EqualTo(4), "7 - 4 + 1 = 4 distinct canonical k-mers");
            Assert.That(db.Values.All(v => v == 100), Is.True);
            Assert.That(db.ContainsKey("AAAA"), Is.True);
        });
    }

    [Test]
    [Description("Kraken DB-build LCA: a k-mer shared by two species collapses to their LCA (genus)")]
    public void BuildKmerDatabase_SharedKmer_CollapsesToLca()
    {
        var t = BuildTaxonomy();
        // Both references contain the (palindromic, self-canonical) k-mer AGCT:
        //   E.coli(100)        ref "AGCTAAAA" → AGCT, GCTA, CTAA, TAAA
        //   E.fergusonii(101)  ref "AGCTCCCC" → AGCT(shared), GCTC(canon GAGC), CTCC, TCCC(canon GGGA)
        // AGCT is owned by 100 and 101 ⇒ stored value = Lca(100,101) = Escherichia(20).
        var db = MetagenomicsAnalyzer.BuildKmerDatabase(
            new[] { (100, "AGCTAAAA"), (101, "AGCTCCCC") }, t, k: 4);

        Assert.Multiple(() =>
        {
            Assert.That(db["AGCT"], Is.EqualTo(20), "shared k-mer → LCA = genus Escherichia(20)");
            Assert.That(db["GCTA"], Is.EqualTo(100), "E.coli-only k-mer stays at species 100");
            Assert.That(db["GAGC"], Is.EqualTo(101), "E.fergusonii-only k-mer (canon of GCTC) stays at species 101");
        });
    }

    [Test]
    [Description("DB build skips ambiguous (non-ACGT) k-mers and uppercases input")]
    public void BuildKmerDatabase_AmbiguousSkipped_MixedCaseHandled()
    {
        var t = BuildTaxonomy();
        // lowercase + an N: only the all-ACGT windows are indexed; all keys uppercase ACGT.
        var db = MetagenomicsAnalyzer.BuildKmerDatabase(
            new[] { (100, "aaaaNcaa") }, t, k: 4);

        Assert.Multiple(() =>
        {
            Assert.That(db.Count, Is.GreaterThan(0));
            foreach (var key in db.Keys)
                Assert.That(key, Does.Match("^[ACGT]+$"), $"k-mer '{key}' must be uppercase ACGT");
        });
    }

    [Test]
    [Description("DB build rejects a reference taxon not present in the taxonomy tree")]
    public void BuildKmerDatabase_UnknownTaxon_Throws()
    {
        var t = BuildTaxonomy();
        Assert.Throws<KeyNotFoundException>(() =>
            MetagenomicsAnalyzer.BuildKmerDatabase(new[] { (999, "AAAACAA") }, t, k: 4));
    }

    #endregion

    #region ClassifyReads tests (RTL / LCA, hand-derived)

    private static MetagenomicsAnalyzer.TaxonomicClassification ClassifyOne(
        string sequence, IReadOnlyDictionary<string, int> db, TaxonomyTree tree, int k = 4)
        => MetagenomicsAnalyzer.ClassifyReads(new[] { ("r", sequence) }, db, tree, k).Single();

    [Test]
    [Description("Read hitting one species → that species (RTL path = single leaf)")]
    public void ClassifyReads_SingleSpecies_AssignsThatSpecies()
    {
        var t = BuildTaxonomy();
        // All 4 windows of "AAAACAA" map to E.coli(100).
        var db = new Dictionary<string, int>
        {
            ["AAAA"] = 100, ["AAAC"] = 100, ["AACA"] = 100, ["ACAA"] = 100,
        };

        var r = ClassifyOne("AAAACAA", db, t);

        Assert.Multiple(() =>
        {
            Assert.That(r.TaxonId, Is.EqualTo(100), "single leaf 100 → assigned E.coli");
            Assert.That(r.TaxonName, Is.EqualTo("Escherichia coli"));
            Assert.That(r.Rank, Is.EqualTo("species"));
            Assert.That(r.RtlScore, Is.EqualTo(4), "RTL score = sum of node weights on root→100 = 4");
            Assert.That(r.TotalKmers, Is.EqualTo(4), "Q = 4 non-ambiguous k-mers");
            Assert.That(r.MatchedKmers, Is.EqualTo(4), "C = clade(100) k-mers = 4");
            Assert.That(r.Confidence, Is.EqualTo(1.0).Within(1e-9));
            Assert.That(r.Genus, Is.EqualTo("Escherichia"));
            Assert.That(r.Species, Is.EqualTo("Escherichia coli"));
        });
    }

    [Test]
    [Description("Read splitting equally between two species of one genus → the genus (LCA of tied leaves)")]
    public void ClassifyReads_SplitWithinGenus_AssignsGenusLca()
    {
        var t = BuildTaxonomy();
        // 2 windows → E.coli(100), 2 windows → E.fergusonii(101).
        // Leaves {100,101} each score 2 (own weight only; sibling paths are disjoint below genus).
        // Tie → Lca(100,101) = Escherichia(20).
        var db = new Dictionary<string, int>
        {
            ["AAAA"] = 100, ["AAAC"] = 100, ["AACA"] = 101, ["ACAA"] = 101,
        };

        var r = ClassifyOne("AAAACAA", db, t);

        Assert.Multiple(() =>
        {
            Assert.That(r.TaxonId, Is.EqualTo(20), "tied species leaves → LCA = genus Escherichia(20)");
            Assert.That(r.TaxonName, Is.EqualTo("Escherichia"));
            Assert.That(r.Rank, Is.EqualTo("genus"));
            Assert.That(r.RtlScore, Is.EqualTo(2), "each tied RTL path scores 2");
            Assert.That(r.TotalKmers, Is.EqualTo(4));
            Assert.That(r.MatchedKmers, Is.EqualTo(4), "C = clade(20) = both species' 4 k-mers");
            Assert.That(r.Confidence, Is.EqualTo(1.0).Within(1e-9));
        });
    }

    [Test]
    [Description("Read splitting equally across two genera → their family (LCA of tied leaves)")]
    public void ClassifyReads_SplitAcrossGenera_AssignsFamilyLca()
    {
        var t = BuildTaxonomy();
        // 2 windows → E.coli(100), 2 windows → S.enterica(200). Tie → Lca(100,200)=family(10).
        var db = new Dictionary<string, int>
        {
            ["AAAA"] = 100, ["AAAC"] = 100, ["AACA"] = 200, ["ACAA"] = 200,
        };

        var r = ClassifyOne("AAAACAA", db, t);

        Assert.Multiple(() =>
        {
            Assert.That(r.TaxonId, Is.EqualTo(10), "tied across genera → LCA = Enterobacteriaceae(10)");
            Assert.That(r.Rank, Is.EqualTo("family"));
            Assert.That(r.RtlScore, Is.EqualTo(2));
            Assert.That(r.MatchedKmers, Is.EqualTo(4), "C = clade(10) = all 4 k-mers");
            Assert.That(r.Confidence, Is.EqualTo(1.0).Within(1e-9));
        });
    }

    [Test]
    [Description("RTL with ancestor weight: a clear single max path wins (no tie)")]
    public void ClassifyReads_RtlAncestorWeight_SingleWinner()
    {
        var t = BuildTaxonomy();
        // Hits: E.coli(100)×1, E.fergusonii(101)×2, genus Escherichia(20)×1 (a collapsed shared k-mer).
        // Leaves = {100,101} (20 is their ancestor → internal node, contributes weight to both paths).
        //   RTL(100) = w(100)+w(20) = 1+1 = 2
        //   RTL(101) = w(101)+w(20) = 2+1 = 3   ← unique maximum
        // Assigned = 101 (E.fergusonii). C = clade(101) = 2; Q = 4.
        var db = new Dictionary<string, int>
        {
            ["AAAA"] = 100, ["AAAC"] = 101, ["AACA"] = 101, ["ACAA"] = 20,
        };

        var r = ClassifyOne("AAAACAA", db, t);

        Assert.Multiple(() =>
        {
            Assert.That(r.TaxonId, Is.EqualTo(101), "max RTL path leaf = E.fergusonii(101)");
            Assert.That(r.RtlScore, Is.EqualTo(3), "RTL(101) = w(101)+w(20) = 2+1 = 3");
            Assert.That(r.TotalKmers, Is.EqualTo(4));
            Assert.That(r.MatchedKmers, Is.EqualTo(2), "C = clade(101) k-mers = 2 (genus hit is outside the clade)");
            Assert.That(r.Confidence, Is.EqualTo(0.5).Within(1e-9), "C/Q = 2/4");
        });
    }

    [Test]
    [Description("Read with no k-mer hits → unclassified (root, taxon 1), confidence 0")]
    public void ClassifyReads_NoHits_Unclassified()
    {
        var t = BuildTaxonomy();
        var db = new Dictionary<string, int> { ["GGGG"] = 100 }; // nothing the read produces
        var r = ClassifyOne("AAAACAA", db, t);

        Assert.Multiple(() =>
        {
            Assert.That(r.TaxonId, Is.EqualTo(TaxonomyTree.RootId), "no hits → root/unclassified");
            Assert.That(r.TaxonName, Is.EqualTo("root"));
            Assert.That(r.RtlScore, Is.EqualTo(0));
            Assert.That(r.MatchedKmers, Is.EqualTo(0));
            Assert.That(r.TotalKmers, Is.EqualTo(4), "Q still counts the 4 queried k-mers");
            Assert.That(r.Confidence, Is.EqualTo(0.0));
        });
    }

    [Test]
    [Description("Empty / short reads → unclassified with Q = 0")]
    public void ClassifyReads_EmptyOrShort_Unclassified()
    {
        var t = BuildTaxonomy();
        var db = new Dictionary<string, int> { ["AAAA"] = 100 };

        var empty = ClassifyOne("", db, t);
        var shortRead = ClassifyOne("AAA", db, t); // shorter than k=4

        Assert.Multiple(() =>
        {
            Assert.That(empty.TaxonId, Is.EqualTo(TaxonomyTree.RootId));
            Assert.That(empty.TotalKmers, Is.EqualTo(0));
            Assert.That(shortRead.TaxonId, Is.EqualTo(TaxonomyTree.RootId));
            Assert.That(shortRead.TotalKmers, Is.EqualTo(0));
        });
    }

    [Test]
    [Description("Ambiguous k-mers are excluded from Q (Kraken: only non-ambiguous k-mers are queried)")]
    public void ClassifyReads_AmbiguousKmers_ExcludedFromQ()
    {
        var t = BuildTaxonomy();
        var db = new Dictionary<string, int> { ["AAAA"] = 100 };
        // "NNNNNNN": every window contains N → Q = 0 → unclassified.
        var r = ClassifyOne("NNNNNNN", db, t);

        Assert.Multiple(() =>
        {
            Assert.That(r.TotalKmers, Is.EqualTo(0), "no non-ambiguous k-mers queried");
            Assert.That(r.TaxonId, Is.EqualTo(TaxonomyTree.RootId));
        });
    }

    [Test]
    [Description("Canonical lookup: a read window's reverse complement matches a self-canonical DB key")]
    public void ClassifyReads_CanonicalLookup_MatchesReverseComplement()
    {
        var t = BuildTaxonomy();
        // DB key AACC (self-canonical; its reverse complement is GGTT) → E.coli(100).
        // Read "AGGTT" windows: AGGT (canon ACCT, not in DB), GGTT (canon AACC → 100).
        // The read contains no literal "AACC" — the match is via canonicalization of GGTT. Q=2.
        var db = new Dictionary<string, int> { ["AACC"] = 100 };
        var r = ClassifyOne("AGGTT", db, t);

        Assert.Multiple(() =>
        {
            Assert.That(r.TaxonId, Is.EqualTo(100), "RC window GGTT canonicalizes to DB key AACC");
            Assert.That(r.TotalKmers, Is.EqualTo(2));
            Assert.That(r.MatchedKmers, Is.EqualTo(1));
        });
    }

    [Test]
    [Description("Output count equals input read count, in order")]
    public void ClassifyReads_OutputCountAndOrder_Preserved()
    {
        var t = BuildTaxonomy();
        var db = new Dictionary<string, int> { ["AAAA"] = 100 };
        var reads = new[]
        {
            ("read1", "AAAACAA"),
            ("read2", ""),
            ("read3", "GGGGGGGG"),
        };

        var results = MetagenomicsAnalyzer.ClassifyReads(reads, db, t, k: 4).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(results.Count, Is.EqualTo(3));
            Assert.That(results.Select(r => r.ReadId), Is.EqualTo(new[] { "read1", "read2", "read3" }));
        });
    }

    #endregion

    #region Invariants

    [Test]
    [Description("Invariant: confidence ∈ [0,1] for every read; MatchedKmers ≤ TotalKmers")]
    public void ClassifyReads_Invariants_Hold()
    {
        var t = BuildTaxonomy();
        var db = new Dictionary<string, int>
        {
            ["AAAA"] = 100, ["AAAC"] = 101, ["AACA"] = 200, ["ACAA"] = 20,
        };
        var reads = new[]
        {
            ("a", "AAAACAA"),
            ("b", "NNNNNNN"),
            ("c", ""),
            ("d", "GGGGGGGG"),
        };

        var results = MetagenomicsAnalyzer.ClassifyReads(reads, db, t, k: 4).ToList();

        foreach (var r in results)
        {
            Assert.That(r.Confidence, Is.InRange(0.0, 1.0), $"confidence for {r.ReadId}");
            Assert.That(r.MatchedKmers, Is.LessThanOrEqualTo(r.TotalKmers), $"C ≤ Q for {r.ReadId}");
        }
    }

    [Test]
    [Description("Invariant: unclassified reads have MatchedKmers = 0 and taxon = root")]
    public void ClassifyReads_Unclassified_ZeroMatched()
    {
        var t = BuildTaxonomy();
        var db = new Dictionary<string, int> { ["AAAA"] = 100 };
        var r = ClassifyOne("GGGGGGGG", db, t);

        Assert.Multiple(() =>
        {
            Assert.That(r.TaxonId, Is.EqualTo(TaxonomyTree.RootId));
            Assert.That(r.MatchedKmers, Is.EqualTo(0));
        });
    }

    [Test]
    [Description("Validation: null/invalid arguments are rejected")]
    public void ClassifyReads_InvalidArguments_Throw()
    {
        var t = BuildTaxonomy();
        var db = new Dictionary<string, int>();
        var reads = new[] { ("r", "AAAACAA") };

        Assert.Multiple(() =>
        {
            Assert.Throws<System.ArgumentNullException>(() =>
                MetagenomicsAnalyzer.ClassifyReads(null!, db, t, 4).ToList());
            Assert.Throws<System.ArgumentNullException>(() =>
                MetagenomicsAnalyzer.ClassifyReads(reads, null!, t, 4).ToList());
            Assert.Throws<System.ArgumentNullException>(() =>
                MetagenomicsAnalyzer.ClassifyReads(reads, db, null!, 4).ToList());
            Assert.Throws<System.ArgumentOutOfRangeException>(() =>
                MetagenomicsAnalyzer.ClassifyReads(reads, db, t, 0).ToList());
        });
    }

    [Test]
    [Description("TaxonomyTree construction rejects malformed trees")]
    public void TaxonomyTree_Construction_ValidatesShape()
    {
        Assert.Multiple(() =>
        {
            // No root.
            Assert.Throws<System.ArgumentException>(() => new TaxonomyTree(new[]
            {
                new TaxonNode(2, "a", "x", 1), // parent 1 missing & not self-parented
            }));
            // Two roots.
            Assert.Throws<System.ArgumentException>(() => new TaxonomyTree(new[]
            {
                new TaxonNode(1, "r1", "root", 1),
                new TaxonNode(2, "r2", "root", 2),
            }));
            // Duplicate id.
            Assert.Throws<System.ArgumentException>(() => new TaxonomyTree(new[]
            {
                new TaxonNode(1, "r", "root", 1),
                new TaxonNode(1, "dup", "x", 1),
            }));
        });
    }

    #endregion
}
