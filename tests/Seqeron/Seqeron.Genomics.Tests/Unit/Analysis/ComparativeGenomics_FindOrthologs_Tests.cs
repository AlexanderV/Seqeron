// COMPGEN-ORTHO-001 — Ortholog Identification (Reciprocal Best Hits) and Paralog Identification
// Evidence: docs/Evidence/COMPGEN-ORTHO-001-Evidence.md
// TestSpec: tests/TestSpecs/COMPGEN-ORTHO-001.md
// Source: Tatusov RL, Koonin EV, Lipman DJ (1997). Science 278:631-637 (symmetrical best hits).
//         Moreno-Hagelsieb G, Latimer K (2008). Bioinformatics 24(3):319-324 (RBH definition).
//         Fitch WM (1970). Syst Zool 19:99-106 (orthology/paralogy). Remm et al. (2001). JMB 314:1041-1052.

using NUnit.Framework;
using Seqeron.Genomics.Analysis;
using System;
using System.Collections.Generic;
using System.Linq;

namespace Seqeron.Genomics.Tests.Unit.Analysis;

[TestFixture]
public class ComparativeGenomics_FindOrthologs_Tests
{
    #region Helpers

    private static ComparativeGenomics.Gene Gene(string id, string genome, string seq)
        => new(id, genome, 0, seq.Length, '+', seq);

    // k-mer sets are computed over these sequences (k = 5):
    //   AcgtRepeat14 / AcgtRepeat14b are identical -> Jaccard 1.0
    //   AcgtPlusAA shares all of AcgtRepeat14's k-mers but adds 2 more -> Jaccard 0.667
    //   GcBlock / TtBlock are unrelated to the ACGT sequences (Jaccard 0.0 vs ACGT).
    //   NOTE: TtBlock and GcBlock are cyclic rotations of one another and share 8 of 12 5-mers
    //   (Jaccard 0.5) — but each is identical to its own genome-2 counterpart (Jaccard 1.0),
    //   so the 1.0 self-match always wins the best-hit ranking and the RBH matching is unaffected.
    private const string AcgtRepeat14 = "ACGTACGTACGTAC";    // k-mers {ACGTA,CGTAC,GTACG,TACGT}
    private const string AcgtPlusAA = "ACGTACGTACGTACAA";  // adds {GTACA,TACAA} -> vs AcgtRepeat14 Jaccard 0.667
    private const string TtBlock = "TTTTGGGGCCCCAAAA";      // 12 5-mers; vs ACGT Jaccard 0.0
    private const string GcBlock = "GGGGCCCCAAAATTTT";      // 12 5-mers; vs ACGT Jaccard 0.0; vs TtBlock Jaccard 0.5

    #endregion

    #region FindOrthologs — MUST Tests

    // M1 — Two mutual best-hit pairs are both returned (RBH).
    // Source: Moreno-Hagelsieb & Latimer (2008): orthologs find each other as best hit.
    [Test]
    public void FindOrthologs_MutualBestHits_ReturnsBothPairs()
    {
        // Arrange: a1==b1 (Jaccard 1.0), a2==b2 (Jaccard 1.0); cross-similarity 0.
        var g1 = new[] { Gene("a1", "G1", AcgtRepeat14), Gene("a2", "G1", TtBlock) };
        var g2 = new[] { Gene("b1", "G2", AcgtRepeat14), Gene("b2", "G2", TtBlock) };

        // Act
        var pairs = ComparativeGenomics.FindOrthologs(g1, g2).ToList();

        // Assert: exactly the two reciprocal pairs, with identity 1.0 each.
        Assert.Multiple(() =>
        {
            Assert.That(pairs, Has.Count.EqualTo(2), "two independent mutual best hits => two ortholog pairs");
            var byG1 = pairs.ToDictionary(p => p.Gene1Id, p => p.Gene2Id);
            Assert.That(byG1["a1"], Is.EqualTo("b1"), "a1 and b1 are mutual best hits (identical sequence)");
            Assert.That(byG1["a2"], Is.EqualTo("b2"), "a2 and b2 are mutual best hits (identical sequence)");
            Assert.That(pairs.Select(p => p.Identity), Is.All.EqualTo(1.0).Within(1e-10),
                "identical sequences => k-mer Jaccard identity 1.0");
            Assert.That(pairs.Select(p => p.Coverage), Is.All.EqualTo(1.0).Within(1e-10),
                "identical sequences => all shared 5-mers => coverage 1.0");
        });
    }

    // M2 — A one-directional best hit is NOT returned as an ortholog (reciprocity required).
    // Source: Tatusov et al. (1997) symmetrical best hits; Moreno-Hagelsieb & Latimer (2008).
    [Test]
    public void FindOrthologs_OneDirectionalBestHit_NotReturnedAsOrtholog()
    {
        // Arrange: a1==b1 (Jaccard 1.0); b2 shares all of a1's k-mers (Jaccard 0.667) so b2's best hit
        // into G1 is a1, but a1's best hit into G2 is b1 (1.0 > 0.667). b2 is therefore non-reciprocal.
        var g1 = new[] { Gene("a1", "G1", AcgtRepeat14) };
        var g2 = new[] { Gene("b1", "G2", AcgtRepeat14), Gene("b2", "G2", AcgtPlusAA) };

        // Act
        var pairs = ComparativeGenomics.FindOrthologs(g1, g2).ToList();

        // Assert: only the reciprocal pair a1<->b1; b2 excluded.
        Assert.Multiple(() =>
        {
            Assert.That(pairs, Has.Count.EqualTo(1), "only the reciprocal best-hit pair qualifies");
            Assert.That(pairs[0].Gene1Id, Is.EqualTo("a1"), "a1 is the genome-1 member");
            Assert.That(pairs[0].Gene2Id, Is.EqualTo("b1"), "a1's reciprocal best hit is b1, not b2");
            Assert.That(pairs.Any(p => p.Gene2Id == "b2"), Is.False,
                "b2's best hit is a1 but a1's best hit is b1 => non-reciprocal => excluded");
            Assert.That(pairs[0].Identity, Is.EqualTo(1.0).Within(1e-10),
                "the kept pair is the Jaccard-1.0 a1<->b1, not the 0.667 a1-b2 hit");
        });
    }

    // M3 — Pairs below the similarity/coverage threshold are excluded.
    // Source: Moreno-Hagelsieb & Latimer (2008): >= 50% coverage required to qualify.
    [Test]
    public void FindOrthologs_BelowThreshold_ReturnsNoPairs()
    {
        // Arrange: unrelated sequences (Jaccard 0.0, shared k-mers 0).
        var g1 = new[] { Gene("a1", "G1", AcgtRepeat14) };
        var g2 = new[] { Gene("b1", "G2", TtBlock) };

        // Act
        var pairs = ComparativeGenomics.FindOrthologs(g1, g2).ToList();

        // Assert
        Assert.That(pairs, Is.Empty, "no shared k-mers => below minIdentity/minCoverage => no orthologs");
    }

    // M4 — Empty genome yields no orthologs.
    // Source: a pair requires one gene from each genome (RBH definition).
    [Test]
    public void FindOrthologs_EmptyGenome_ReturnsEmpty()
    {
        // Arrange
        var g1 = new[] { Gene("a1", "G1", AcgtRepeat14) };
        var g2 = Array.Empty<ComparativeGenomics.Gene>();

        // Act
        var pairs = ComparativeGenomics.FindOrthologs(g1, g2).ToList();

        // Assert
        Assert.That(pairs, Is.Empty, "no genes in genome 2 => no ortholog pairs possible");
    }

    // M5 — Null inputs throw ArgumentNullException.
    [Test]
    public void FindOrthologs_NullInput_Throws()
    {
        var g = new[] { Gene("a1", "G1", AcgtRepeat14) };
        Assert.Multiple(() =>
        {
            Assert.Throws<ArgumentNullException>(
                () => ComparativeGenomics.FindOrthologs(null!, g).ToList(),
                "null genome1 must throw");
            Assert.Throws<ArgumentNullException>(
                () => ComparativeGenomics.FindOrthologs(g, null!).ToList(),
                "null genome2 must throw");
        });
    }

    #endregion

    #region FindOrthologs — SHOULD Tests

    // S1 — RBH is a matching: no genome-1 or genome-2 gene appears in two pairs (INV-2).
    [Test]
    public void FindOrthologs_MultipleGenes_YieldsMatchingNoGeneTwice()
    {
        // Arrange: three independent identical pairs.
        var g1 = new[] { Gene("a1", "G1", AcgtRepeat14), Gene("a2", "G1", TtBlock), Gene("a3", "G1", GcBlock) };
        var g2 = new[] { Gene("b1", "G2", AcgtRepeat14), Gene("b2", "G2", TtBlock), Gene("b3", "G2", GcBlock) };

        // Act
        var pairs = ComparativeGenomics.FindOrthologs(g1, g2).ToList();

        // Assert: each side's gene ids are distinct across all pairs.
        Assert.Multiple(() =>
        {
            Assert.That(pairs.Select(p => p.Gene1Id).Distinct().Count(), Is.EqualTo(pairs.Count),
                "no genome-1 gene appears in two pairs");
            Assert.That(pairs.Select(p => p.Gene2Id).Distinct().Count(), Is.EqualTo(pairs.Count),
                "no genome-2 gene appears in two pairs");
        });
    }

    // S2 — Every returned pair is reciprocal (INV-1): re-derive best hits and confirm symmetry.
    [Test]
    public void FindOrthologs_AllPairs_AreReciprocal()
    {
        // Arrange
        var g1 = new[] { Gene("a1", "G1", AcgtRepeat14), Gene("a2", "G1", TtBlock) };
        var g2 = new[] { Gene("b1", "G2", AcgtRepeat14), Gene("b2", "G2", TtBlock), Gene("b3", "G2", AcgtPlusAA) };

        // Act
        var pairs = ComparativeGenomics.FindOrthologs(g1, g2).ToList();
        // Independent re-check: best hit of g2 member back into g1 must be the g1 member.
        var backHits = ComparativeGenomics.FindOrthologs(g2, g1)
            .ToDictionary(p => p.Gene1Id, p => p.Gene2Id);

        // Assert
        Assert.Multiple(() =>
        {
            foreach (var p in pairs)
                Assert.That(backHits.TryGetValue(p.Gene2Id, out var back) && back == p.Gene1Id, Is.True,
                    $"pair {p.Gene1Id}<->{p.Gene2Id} must be reciprocal in both call directions");
        });
    }

    // S3 — A gene with no sequence is skipped (never paired).
    [Test]
    public void FindOrthologs_GeneWithoutSequence_IsSkipped()
    {
        // Arrange: a2 has no sequence.
        var g1 = new[]
        {
            Gene("a1", "G1", AcgtRepeat14),
            new ComparativeGenomics.Gene("a2", "G1", 0, 0, '+', Sequence: null),
        };
        var g2 = new[] { Gene("b1", "G2", AcgtRepeat14) };

        // Act
        var pairs = ComparativeGenomics.FindOrthologs(g1, g2).ToList();

        // Assert
        Assert.Multiple(() =>
        {
            Assert.That(pairs, Has.Count.EqualTo(1), "only the sequenced gene a1 can be an ortholog");
            Assert.That(pairs.Any(p => p.Gene1Id == "a2"), Is.False, "a2 has no sequence => skipped");
        });
    }

    // S5 — The public FindReciprocalBestHits entry point yields the same matching as FindOrthologs.
    // FindOrthologs delegates to FindReciprocalBestHits; assert the dedicated RBH entry point directly
    // (it is public and otherwise only exercised indirectly). Source: Moreno-Hagelsieb & Latimer (2008).
    [Test]
    public void FindReciprocalBestHits_SameAsFindOrthologs()
    {
        // Arrange
        var g1 = new[] { Gene("a1", "G1", AcgtRepeat14), Gene("a2", "G1", TtBlock) };
        var g2 = new[] { Gene("b1", "G2", AcgtRepeat14), Gene("b2", "G2", TtBlock) };

        // Act
        var rbh = ComparativeGenomics.FindReciprocalBestHits(g1, g2)
            .Select(p => (p.Gene1Id, p.Gene2Id)).OrderBy(p => p.Item1).ToList();
        var ortho = ComparativeGenomics.FindOrthologs(g1, g2)
            .Select(p => (p.Gene1Id, p.Gene2Id)).OrderBy(p => p.Item1).ToList();

        // Assert
        Assert.Multiple(() =>
        {
            Assert.That(rbh, Is.EqualTo(ortho), "FindOrthologs is the RBH criterion under another name");
            Assert.That(rbh, Is.EqualTo(new[] { ("a1", "b1"), ("a2", "b2") }),
                "the two reciprocal best-hit pairs");
        });
    }

    // S6 — FindReciprocalBestHits also validates null inputs (its own ArgumentNullException guards).
    [Test]
    public void FindReciprocalBestHits_NullInput_Throws()
    {
        var g = new[] { Gene("a1", "G1", AcgtRepeat14) };
        Assert.Multiple(() =>
        {
            Assert.Throws<ArgumentNullException>(
                () => ComparativeGenomics.FindReciprocalBestHits(null!, g).ToList(), "null genome1 must throw");
            Assert.Throws<ArgumentNullException>(
                () => ComparativeGenomics.FindReciprocalBestHits(g, null!).ToList(), "null genome2 must throw");
        });
    }

    #endregion

    #region FindOrthologs — COULD Tests

    // C1 — Deterministic: same input yields identical pair set across runs (INV-4).
    [Test]
    public void FindOrthologs_RunTwice_IsDeterministic()
    {
        // Arrange
        var g1 = new[] { Gene("a1", "G1", AcgtRepeat14), Gene("a2", "G1", TtBlock) };
        var g2 = new[] { Gene("b1", "G2", AcgtRepeat14), Gene("b2", "G2", TtBlock) };

        // Act
        var first = ComparativeGenomics.FindOrthologs(g1, g2).Select(p => (p.Gene1Id, p.Gene2Id)).ToList();
        var second = ComparativeGenomics.FindOrthologs(g1, g2).Select(p => (p.Gene1Id, p.Gene2Id)).ToList();

        // Assert
        Assert.That(second, Is.EqualTo(first), "RBH is deterministic for fixed input");
    }

    #endregion

    #region FindParalogs — MUST Tests

    // M6 — A duplicated within-genome gene pair is returned as a paralog pair; unrelated gene excluded.
    // Source: Fitch (1970) paralogy = within-genome duplication; Remm et al. (2001) in-paralog rule.
    [Test]
    public void FindParalogs_DuplicateGene_ReturnsParalogPair()
    {
        // Arrange: p1==p2 (Jaccard 1.0); q1 unrelated (Jaccard 0.0 to both).
        var genes = new[]
        {
            Gene("p1", "G1", GcBlock),
            Gene("p2", "G1", GcBlock),
            Gene("q1", "G1", AcgtRepeat14),
        };

        // Act
        var pairs = ComparativeGenomics.FindParalogs(genes).ToList();

        // Assert: exactly the {p1,p2} mutual best-hit pair.
        Assert.Multiple(() =>
        {
            Assert.That(pairs, Has.Count.EqualTo(1), "one within-genome mutual best-hit pair");
            var ids = new HashSet<string> { pairs[0].Gene1Id, pairs[0].Gene2Id };
            Assert.That(ids, Is.EquivalentTo(new[] { "p1", "p2" }), "p1 and p2 are mutual best hits");
            Assert.That(ids.Contains("q1"), Is.False, "unrelated q1 is not a paralog of either");
            Assert.That(pairs[0].Identity, Is.EqualTo(1.0).Within(1e-10), "duplicate => Jaccard identity 1.0");
            Assert.That(pairs[0].Coverage, Is.EqualTo(1.0).Within(1e-10), "duplicate => coverage 1.0");
        });
    }

    // M7 — A single-gene genome yields no paralogs.
    // Source: a paralog pair requires two within-genome genes.
    [Test]
    public void FindParalogs_SingleGene_ReturnsEmpty()
    {
        // Arrange
        var genes = new[] { Gene("p1", "G1", GcBlock) };

        // Act
        var pairs = ComparativeGenomics.FindParalogs(genes).ToList();

        // Assert
        Assert.That(pairs, Is.Empty, "fewer than two sequenced genes => no paralog pair possible");
    }

    // M8 — Null input throws ArgumentNullException.
    [Test]
    public void FindParalogs_NullInput_Throws()
    {
        Assert.Throws<ArgumentNullException>(
            () => ComparativeGenomics.FindParalogs(null!).ToList(),
            "null gene list must throw");
    }

    #endregion

    #region FindParalogs — SHOULD Tests

    // S4 — Paralog pair members are distinct genes (no gene paired with itself) (INV-3).
    [Test]
    public void FindParalogs_PairGenesAreDistinct()
    {
        // Arrange
        var genes = new[] { Gene("p1", "G1", GcBlock), Gene("p2", "G1", GcBlock) };

        // Act
        var pairs = ComparativeGenomics.FindParalogs(genes).ToList();

        // Assert
        Assert.Multiple(() =>
        {
            Assert.That(pairs, Has.Count.EqualTo(1), "one paralog pair");
            Assert.That(pairs[0].Gene1Id, Is.Not.EqualTo(pairs[0].Gene2Id),
                "a paralog pair is two distinct genes, never a gene with itself");
        });
    }

    #endregion
}
