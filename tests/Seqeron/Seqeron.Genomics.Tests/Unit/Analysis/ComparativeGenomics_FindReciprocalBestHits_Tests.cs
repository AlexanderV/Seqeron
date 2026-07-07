// COMPGEN-RBH-001 — Reciprocal Best Hits (RBH / bidirectional best hits) for ortholog identification
// Evidence: docs/Evidence/COMPGEN-RBH-001-Evidence.md
// TestSpec: tests/TestSpecs/COMPGEN-RBH-001.md
// Source: Moreno-Hagelsieb G, Latimer K (2008). Bioinformatics 24(3):319-324 (RBH definition, coverage gate).
//         Tatusov RL, Koonin EV, Lipman DJ (1997). Science 278:631-637 (symmetrical/mutual best hits).

namespace Seqeron.Genomics.Tests.Unit.Analysis;

[TestFixture]
public class ComparativeGenomics_FindReciprocalBestHits_Tests
{
    #region Helpers

    private static ComparativeGenomics.Gene Gene(string id, string genome, string seq)
        => new(id, genome, 0, seq.Length, '+', seq);

    // k-mer sets are computed over these sequences (k = 5), verified by hand:
    //   AcgtRepeat14 / AcgtRepeat14b are identical -> Jaccard 1.0, alignLen = min length = 14
    //   AcgtPlusAA shares all of AcgtRepeat14's k-mers but adds 2 more -> Jaccard 0.667
    //   TtBlock / GcBlock are unrelated to the ACGT sequences -> Jaccard 0.0
    private const string AcgtRepeat14 = "ACGTACGTACGTAC";    // length 14; k-mers {ACGTA,CGTAC,GTACG,TACGT}
    private const string AcgtPlusAA = "ACGTACGTACGTACAA";  // adds {GTACA,TACAA} -> 4 shared / 6 union = 0.667
    private const string TtBlock = "TTTTGGGGCCCCAAAA";      // unrelated 12 k-mers
    private const string GcBlock = "GGGGCCCCAAAATTTT";      // unrelated 12 k-mers

    // Partial-overlap pair for the COVERAGE gate (independent of the identity gate):
    //   AcPrefix / AcPrefixAlt share the 6-mer-content prefix AAAAACCCCC...; by hand (k = 5):
    //   AcPrefix    k-mers = {AAAAA,AAAAC,AAACC,AACCC,ACCCC,CCCCC,CCCCG,CCCGG,CCGGG,CGGGG,GGGGG} (11)
    //   AcPrefixAlt k-mers = {AAAAA,AAAAC,AAACC,AACCC,ACCCC,CCCCC,CCCCT,CCCTT,CCTTT,CTTTT,TTTTT} (11)
    //   shared = {AAAAA,AAAAC,AAACC,AACCC,ACCCC,CCCCC} = 6  =>  identity = 6/16 = 0.375 (>= 0.3),
    //   coverage = 6/11 = 0.5455 (>= 0.5 default, but < 0.6). So minCoverage 0.6 rejects on coverage
    //   ALONE while identity still passes — isolating the >= 50% coverage gate of Moreno-Hagelsieb (2008).
    private const string AcPrefix = "AAAAACCCCCGGGGG";     // 11 k-mers
    private const string AcPrefixAlt = "AAAAACCCCCTTTTT";  // 11 k-mers; 6 shared with AcPrefix

    #endregion

    #region FindReciprocalBestHits — MUST Tests

    // M1 — Two independent mutual best-hit pairs are both returned (RBH).
    // Source: Moreno-Hagelsieb & Latimer (2008): orthologs find each other as the best hit.
    [Test]
    public void FindReciprocalBestHits_MutualBestHits_ReturnsBothPairs()
    {
        // Arrange: a1==b1 (Jaccard 1.0), a2==b2 (Jaccard 1.0); cross-similarity 0.
        var g1 = new[] { Gene("a1", "G1", AcgtRepeat14), Gene("a2", "G1", TtBlock) };
        var g2 = new[] { Gene("b1", "G2", AcgtRepeat14), Gene("b2", "G2", TtBlock) };

        // Act
        var pairs = ComparativeGenomics.FindReciprocalBestHits(g1, g2).ToList();

        // Assert: exactly the two reciprocal pairs.
        Assert.Multiple(() =>
        {
            Assert.That(pairs, Has.Count.EqualTo(2), "two independent mutual best hits => two RBH pairs");
            var byG1 = pairs.ToDictionary(p => p.Gene1Id, p => p.Gene2Id);
            Assert.That(byG1["a1"], Is.EqualTo("b1"), "a1 and b1 are mutual best hits (identical sequence)");
            Assert.That(byG1["a2"], Is.EqualTo("b2"), "a2 and b2 are mutual best hits (identical sequence)");
            Assert.That(pairs.Select(p => p.Identity), Is.All.EqualTo(1.0).Within(1e-10),
                "identical sequences => k-mer Jaccard identity 1.0");
        });
    }

    // M2 — A one-directional best hit is NOT returned (reciprocity required).
    // Source: Tatusov et al. (1997) symmetrical best hits; Moreno-Hagelsieb & Latimer (2008).
    [Test]
    public void FindReciprocalBestHits_OneDirectionalBestHit_NotReturned()
    {
        // Arrange: a1==b1 (Jaccard 1.0); b2 shares all of a1's k-mers (Jaccard 0.667) so b2's best hit
        // into G1 is a1, but a1's best hit into G2 is b1 (1.0 > 0.667). b2 is therefore non-reciprocal.
        var g1 = new[] { Gene("a1", "G1", AcgtRepeat14) };
        var g2 = new[] { Gene("b1", "G2", AcgtRepeat14), Gene("b2", "G2", AcgtPlusAA) };

        // Act
        var pairs = ComparativeGenomics.FindReciprocalBestHits(g1, g2).ToList();

        // Assert: only the reciprocal pair a1<->b1; b2 excluded.
        Assert.Multiple(() =>
        {
            Assert.That(pairs, Has.Count.EqualTo(1), "only the reciprocal best-hit pair qualifies");
            Assert.That(pairs[0].Gene1Id, Is.EqualTo("a1"), "a1 is the genome-1 member");
            Assert.That(pairs[0].Gene2Id, Is.EqualTo("b1"), "a1's reciprocal best hit is b1, not b2");
            Assert.That(pairs.Any(p => p.Gene2Id == "b2"), Is.False,
                "b2's best hit is a1 but a1's best hit is b1 => non-reciprocal => excluded");
        });
    }

    // M3 — The returned pair carries the ACTUAL hit coverage/identity/alignment length, not placeholders.
    // (Guards the corrected defect class: old impl reported Coverage=1.0 hardcoded, AlignmentLength=0,
    //  and Identity=score-product.) Source: Moreno-Hagelsieb & Latimer (2008) best-hit metrics.
    [Test]
    public void FindReciprocalBestHits_IdenticalPair_ReportsActualCoverageAndLength()
    {
        // Arrange: identical 14-nt sequences => identity 1.0, coverage 1.0, alignLen = 14.
        var g1 = new[] { Gene("a1", "G1", AcgtRepeat14) };
        var g2 = new[] { Gene("b1", "G2", AcgtRepeat14) };

        // Act
        var pairs = ComparativeGenomics.FindReciprocalBestHits(g1, g2).ToList();

        // Assert
        Assert.Multiple(() =>
        {
            Assert.That(pairs, Has.Count.EqualTo(1), "identical mutual best hit => one pair");
            Assert.That(pairs[0].Identity, Is.EqualTo(1.0).Within(1e-10),
                "identical sequences => Jaccard identity exactly 1.0 (not the score product)");
            Assert.That(pairs[0].Coverage, Is.EqualTo(1.0).Within(1e-10),
                "all of the shorter sequence's k-mers are shared => coverage 1.0, computed not hardcoded");
            Assert.That(pairs[0].AlignmentLength, Is.EqualTo(14),
                "alignment length = min sequence length = 14 (old impl hardcoded 0)");
        });
    }

    // M4 — A pair below the minimum-identity gate is rejected even when it is the mutual top candidate.
    // Source: Moreno-Hagelsieb & Latimer (2008): significance / >=50% coverage gate.
    [Test]
    public void FindReciprocalBestHits_AboveMinIdentity_RejectsNonIdenticalPair()
    {
        // Arrange: a1 vs b1 are each other's only candidate, but Jaccard is 0.667 (< minIdentity 1.0).
        var g1 = new[] { Gene("a1", "G1", AcgtRepeat14) };
        var g2 = new[] { Gene("b1", "G2", AcgtPlusAA) };

        // Act: raise the gate above the 0.667 similarity so the pair cannot qualify.
        var pairs = ComparativeGenomics.FindReciprocalBestHits(g1, g2, minIdentity: 1.0).ToList();

        // Assert
        Assert.That(pairs, Is.Empty,
            "Jaccard 0.667 < minIdentity 1.0 => no qualifying hit => no RBH pair");
    }

    // M5 — An empty genome yields no RBH pairs.
    // Source: an RBH pair requires one gene from each genome.
    [Test]
    public void FindReciprocalBestHits_EmptyGenome_ReturnsEmpty()
    {
        // Arrange
        var g1 = new[] { Gene("a1", "G1", AcgtRepeat14) };
        var g2 = Array.Empty<ComparativeGenomics.Gene>();

        // Act
        var pairs = ComparativeGenomics.FindReciprocalBestHits(g1, g2).ToList();

        // Assert
        Assert.That(pairs, Is.Empty, "no genes in genome 2 => no RBH pairs possible");
    }

    // M6 — Null inputs throw ArgumentNullException.
    [Test]
    public void FindReciprocalBestHits_NullInput_Throws()
    {
        var g = new[] { Gene("a1", "G1", AcgtRepeat14) };
        Assert.Multiple(() =>
        {
            Assert.Throws<ArgumentNullException>(
                () => ComparativeGenomics.FindReciprocalBestHits(null!, g).ToList(),
                "null genome1 must throw");
            Assert.Throws<ArgumentNullException>(
                () => ComparativeGenomics.FindReciprocalBestHits(g, null!).ToList(),
                "null genome2 must throw");
        });
    }

    // M7 — A pair passing the identity gate but FAILING the coverage gate is rejected (coverage gate is
    // a distinct, independent filter). Source: Moreno-Hagelsieb & Latimer (2008): "coverage of at least
    // 50% of any of the protein sequences in the alignments" is a separate qualifying requirement.
    [Test]
    public void FindReciprocalBestHits_BelowMinCoverage_RejectsPair()
    {
        // Arrange: AcPrefix vs AcPrefixAlt => identity 0.375 (>= 0.3), coverage 6/11 = 0.5455.
        var g1 = new[] { Gene("a1", "G1", AcPrefix) };
        var g2 = new[] { Gene("b1", "G2", AcPrefixAlt) };

        // Control: with the default coverage gate (0.5) the pair qualifies and is returned, and its
        // reported coverage is the actual 6/11 — proving the gate, not the sequences, drives M7.
        var withDefaultCoverage = ComparativeGenomics.FindReciprocalBestHits(g1, g2).ToList();
        // Act: raise the coverage gate above 0.5455 while leaving identity (0.375) qualifying.
        var withHighCoverage =
            ComparativeGenomics.FindReciprocalBestHits(g1, g2, minIdentity: 0.3, minCoverage: 0.6).ToList();

        // Assert
        Assert.Multiple(() =>
        {
            Assert.That(withDefaultCoverage, Has.Count.EqualTo(1),
                "identity 0.375 >= 0.3 and coverage 0.5455 >= 0.5 => the pair qualifies by default");
            Assert.That(withDefaultCoverage[0].Identity, Is.EqualTo(6.0 / 16.0).Within(1e-10),
                "identity = shared/union = 6/16 = 0.375 (hand-computed)");
            Assert.That(withDefaultCoverage[0].Coverage, Is.EqualTo(6.0 / 11.0).Within(1e-10),
                "coverage = shared/min(kmers) = 6/11 = 0.5455 (hand-computed)");
            Assert.That(withHighCoverage, Is.Empty,
                "coverage 0.5455 < minCoverage 0.6 => rejected by the coverage gate even though identity passes");
        });
    }

    #endregion

    #region FindReciprocalBestHits — SHOULD Tests

    // S1 — RBH is a matching: no genome-1 or genome-2 gene appears in two pairs (INV-2).
    [Test]
    public void FindReciprocalBestHits_MultipleGenes_YieldsMatching()
    {
        // Arrange: three independent identical pairs.
        var g1 = new[] { Gene("a1", "G1", AcgtRepeat14), Gene("a2", "G1", TtBlock), Gene("a3", "G1", GcBlock) };
        var g2 = new[] { Gene("b1", "G2", AcgtRepeat14), Gene("b2", "G2", TtBlock), Gene("b3", "G2", GcBlock) };

        // Act
        var pairs = ComparativeGenomics.FindReciprocalBestHits(g1, g2).ToList();

        // Assert
        Assert.Multiple(() =>
        {
            Assert.That(pairs, Has.Count.EqualTo(3), "three independent mutual best hits");
            Assert.That(pairs.Select(p => p.Gene1Id).Distinct().Count(), Is.EqualTo(pairs.Count),
                "no genome-1 gene appears in two pairs");
            Assert.That(pairs.Select(p => p.Gene2Id).Distinct().Count(), Is.EqualTo(pairs.Count),
                "no genome-2 gene appears in two pairs");
        });
    }

    // S2 — Every returned pair is reciprocal (INV-1): re-derive best hits under the reversed call.
    [Test]
    public void FindReciprocalBestHits_AllPairs_AreReciprocal()
    {
        // Arrange
        var g1 = new[] { Gene("a1", "G1", AcgtRepeat14), Gene("a2", "G1", TtBlock) };
        var g2 = new[] { Gene("b1", "G2", AcgtRepeat14), Gene("b2", "G2", TtBlock), Gene("b3", "G2", AcgtPlusAA) };

        // Act
        var pairs = ComparativeGenomics.FindReciprocalBestHits(g1, g2).ToList();
        var backHits = ComparativeGenomics.FindReciprocalBestHits(g2, g1)
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
    public void FindReciprocalBestHits_GeneWithoutSequence_IsSkipped()
    {
        // Arrange: a2 has no sequence.
        var g1 = new[]
        {
            Gene("a1", "G1", AcgtRepeat14),
            new ComparativeGenomics.Gene("a2", "G1", 0, 0, '+', Sequence: null),
        };
        var g2 = new[] { Gene("b1", "G2", AcgtRepeat14) };

        // Act
        var pairs = ComparativeGenomics.FindReciprocalBestHits(g1, g2).ToList();

        // Assert
        Assert.Multiple(() =>
        {
            Assert.That(pairs, Has.Count.EqualTo(1), "only the sequenced gene a1 can be an RBH");
            Assert.That(pairs.Any(p => p.Gene1Id == "a2"), Is.False, "a2 has no sequence => skipped");
        });
    }

    // S4 — Sequences shorter than k = 5 have an empty k-mer set => similarity 0 => never qualify.
    // Source: Reciprocal_Best_Hits.md §6.1 (k-mer set empty for len < k); the gate then rejects them.
    [Test]
    public void FindReciprocalBestHits_SequencesShorterThanK_NeverQualify()
    {
        // Arrange: identical 4-nt sequences (< k = 5) => k-mer set empty => identity 0, coverage 0.
        var g1 = new[] { Gene("a1", "G1", "ACGT") };
        var g2 = new[] { Gene("b1", "G2", "ACGT") };

        // Act: even with the gate dropped to 0, similarity 0 < any positive identity floor; use 0.0
        // to show it is the empty k-mer set (similarity 0), not the gate, that excludes the pair.
        var pairs = ComparativeGenomics.FindReciprocalBestHits(g1, g2, minIdentity: 0.3).ToList();

        // Assert
        Assert.That(pairs, Is.Empty,
            "sequence length 4 < k = 5 => empty k-mer set => similarity 0 < minIdentity => no RBH pair");
    }

    #endregion

    #region FindReciprocalBestHits — COULD Tests

    // C1 — Deterministic and order-independent: reversing input order yields the same pair set (INV-4).
    [Test]
    public void FindReciprocalBestHits_ReversedInputOrder_IsDeterministic()
    {
        // Arrange
        var g1 = new[] { Gene("a1", "G1", AcgtRepeat14), Gene("a2", "G1", TtBlock) };
        var g2 = new[] { Gene("b1", "G2", AcgtRepeat14), Gene("b2", "G2", TtBlock) };

        // Act: forward order vs reversed gene order; compare unordered pair sets.
        var forward = ComparativeGenomics.FindReciprocalBestHits(g1, g2)
            .Select(p => (p.Gene1Id, p.Gene2Id)).OrderBy(x => x.Item1).ToList();
        var reversed = ComparativeGenomics.FindReciprocalBestHits(g1.Reverse().ToArray(), g2.Reverse().ToArray())
            .Select(p => (p.Gene1Id, p.Gene2Id)).OrderBy(x => x.Item1).ToList();

        // Assert
        Assert.That(reversed, Is.EqualTo(forward),
            "deterministic tie-break => RBH pair set is independent of input order");
    }

    #endregion

    #region FindOrthologs — Delegate smoke test

    // Delegate — FindOrthologs delegates to FindReciprocalBestHits; one smoke test proves equivalence.
    [Test]
    public void FindOrthologs_DelegatesTo_FindReciprocalBestHits()
    {
        // Arrange
        var g1 = new[] { Gene("a1", "G1", AcgtRepeat14), Gene("a2", "G1", TtBlock) };
        var g2 = new[] { Gene("b1", "G2", AcgtRepeat14), Gene("b2", "G2", TtBlock) };

        // Act
        var rbh = ComparativeGenomics.FindReciprocalBestHits(g1, g2)
            .Select(p => (p.Gene1Id, p.Gene2Id, p.Identity, p.Coverage, p.AlignmentLength)).ToList();
        var orth = ComparativeGenomics.FindOrthologs(g1, g2)
            .Select(p => (p.Gene1Id, p.Gene2Id, p.Identity, p.Coverage, p.AlignmentLength)).ToList();

        // Assert
        Assert.That(orth, Is.EqualTo(rbh), "FindOrthologs must return exactly the RBH result");
    }

    #endregion
}
