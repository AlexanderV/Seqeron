// ASSEMBLY-SCAFFOLD-001 — Scaffolding (joining ordered contigs with N-gaps)
// Evidence: docs/Evidence/ASSEMBLY-SCAFFOLD-001-Evidence.md
// TestSpec: tests/TestSpecs/ASSEMBLY-SCAFFOLD-001.md
// Source: Jackman SD et al. (2017). ABySS 2.0. Genome Research 27:768-777
//         (scaffold construction); NCBI AGP Specification v2.1 (unknown-gap = 100 N).

namespace Seqeron.Genomics.Tests.Unit.Alignment;

[TestFixture]
public class SequenceAssembler_Scaffold_Tests
{
    #region Scaffold

    // M1 — Two positive-gap links chain three contigs into one scaffold, each pair separated by a
    // run of N of length = the distance estimate. Source: Jackman et al. (2017) scaffold construction.
    [Test]
    public void Scaffold_TwoPositiveGapLinks_ChainsWithExactNRuns()
    {
        var contigs = new[] { "ACGT", "TTGG", "CCAA" };
        var links = new[] { (0, 1, 3), (1, 2, 2) };

        var scaffolds = SequenceAssembler.Scaffold(contigs, links);

        Assert.Multiple(() =>
        {
            Assert.That(scaffolds.Count, Is.EqualTo(1),
                "Both links are followed into a single path, so exactly one scaffold is produced (INV-04).");
            Assert.That(scaffolds[0], Is.EqualTo("ACGTNNNTTGGNNCCAA"),
                "Contigs concatenated interspersed with runs of N of length 3 then 2 (Jackman et al. 2017).");
        });
    }

    // M2 — A positive gap of size g emits exactly g gap characters between the two contigs.
    // Source: Jackman et al. (2017) — run length corresponds to the distance estimate (INV-01).
    [Test]
    public void Scaffold_PositiveGap_EmitsExactlyThatManyGapChars()
    {
        var contigs = new[] { "AA", "TT" };
        var links = new[] { (0, 1, 5) };

        var scaffolds = SequenceAssembler.Scaffold(contigs, links);

        Assert.Multiple(() =>
        {
            Assert.That(scaffolds[0], Is.EqualTo("AANNNNNTT"),
                "A gap estimate of 5 emits exactly 5 N between AA and TT (INV-01, Jackman et al. 2017).");
            Assert.That(scaffolds[0].Count(c => c == 'N'), Is.EqualTo(5),
                "Exactly 5 gap characters, no more and no fewer (run length = distance estimate).");
        });
    }

    // M3 — A negative estimate is a gap of unknown size: the GenBank/EMBL/DDBJ standard of 100 N is
    // emitted. Source: NCBI AGP Specification v2.1 (INV-02); Jackman et al. (2017) (negative=overlap).
    [Test]
    public void Scaffold_NegativeGap_EmitsHundredGapChars()
    {
        var contigs = new[] { "AAAA", "TTTT" };
        var links = new[] { (0, 1, -5) };

        var scaffolds = SequenceAssembler.Scaffold(contigs, links);

        Assert.Multiple(() =>
        {
            Assert.That(scaffolds[0].Count(c => c == 'N'), Is.EqualTo(100),
                "A negative estimate maps to the AGP unknown-gap length of 100 N (NCBI AGP Spec v2.1).");
            Assert.That(scaffolds[0], Is.EqualTo("AAAA" + new string('N', 100) + "TTTT"),
                "Scaffold = AAAA + 100 N + TTTT for an unresolved negative (overlap) estimate.");
            Assert.That(scaffolds[0].Length, Is.EqualTo(108),
                "Length = 4 + 100 + 4 = 108 (INV-03).");
        });
    }

    // M4 — A zero estimate is invalid as a 0-length gap, so it falls into the unknown-gap default of
    // 100 N. Source: NCBI AGP Specification v2.1 ("gap lines with zero length are not valid").
    [Test]
    public void Scaffold_ZeroGap_EmitsHundredGapChars()
    {
        var contigs = new[] { "AAAA", "TTTT" };
        var links = new[] { (0, 1, 0) };

        var scaffolds = SequenceAssembler.Scaffold(contigs, links);

        Assert.That(scaffolds[0].Count(c => c == 'N'), Is.EqualTo(100),
            "A zero-length gap is invalid (AGP) and maps to the unknown-gap default of 100 N (INV-02).");
    }

    // M5 — The gap-fill character is parameterized: a custom character is used verbatim in place of
    // N. Source: Jackman et al. (2017) gap is "a run of the character" (here parameterized).
    [Test]
    public void Scaffold_CustomGapCharacter_UsesItVerbatim()
    {
        var contigs = new[] { "AA", "TT" };
        var links = new[] { (0, 1, 2) };

        var scaffolds = SequenceAssembler.Scaffold(contigs, links, gapCharacter: 'X');

        Assert.That(scaffolds[0], Is.EqualTo("AAXXTT"),
            "The custom gap character 'X' fills the 2-character gap instead of N.");
    }

    // M6 — With no links, every contig is a length-1 path and becomes its own scaffold.
    // Source: Pop et al. (2004) — scaffold is a path of contigs (here trivial paths).
    [Test]
    public void Scaffold_NoLinks_OneScaffoldPerContig()
    {
        var contigs = new[] { "AAA", "CCC" };
        var links = Array.Empty<(int, int, int)>();

        var scaffolds = SequenceAssembler.Scaffold(contigs, links);

        Assert.Multiple(() =>
        {
            Assert.That(scaffolds.Count, Is.EqualTo(2),
                "No links means two single-contig scaffolds (length-1 paths).");
            Assert.That(scaffolds, Is.EqualTo(new[] { "AAA", "CCC" }),
                "Each contig is emitted verbatim, in ascending index order (INV-05).");
        });
    }

    // M7 — Length invariant over a followed path: total length = sum of contig lengths + sum of gaps.
    // Source: Jackman et al. (2017) concatenation + gap runs (INV-03).
    [Test]
    public void Scaffold_FollowedPath_LengthEqualsContigsPlusGaps()
    {
        var contigs = new[] { "ACGT", "TTGG", "CCAA" };
        var links = new[] { (0, 1, 3), (1, 2, 2) };

        var scaffolds = SequenceAssembler.Scaffold(contigs, links);

        const int contigChars = 4 + 4 + 4;
        const int gapChars = 3 + 2;
        Assert.That(scaffolds[0].Length, Is.EqualTo(contigChars + gapChars),
            "Scaffold length = sum of contig lengths (12) + sum of gaps (5) = 17 (INV-03).");
    }

    // S1 — A link to an already-placed contig is skipped: a contig is never used twice (INV-04).
    [Test]
    public void Scaffold_LinkToAlreadyPlacedContig_ContigNotReused()
    {
        var contigs = new[] { "AA", "TT" };
        // (0,1,1) places contig 1 after 0; (1,0,1) would re-place contig 0, which is already used.
        var links = new[] { (0, 1, 1), (1, 0, 1) };

        var scaffolds = SequenceAssembler.Scaffold(contigs, links);

        Assert.Multiple(() =>
        {
            Assert.That(scaffolds.Count, Is.EqualTo(1),
                "Contig 0 is already placed, so the back-link is skipped; one scaffold results (INV-04).");
            Assert.That(scaffolds[0], Is.EqualTo("AANTT"),
                "AA + 1 N + TT; contig 0 is not appended a second time.");
        });
    }

    // S2 — A link whose endpoint is out of range is ignored; the contigs remain unjoined.
    [Test]
    public void Scaffold_OutOfRangeLink_Ignored()
    {
        var contigs = new[] { "AA", "TT" };
        var links = new[] { (0, 5, 2) };

        var scaffolds = SequenceAssembler.Scaffold(contigs, links);

        Assert.That(scaffolds, Is.EqualTo(new[] { "AA", "TT" }),
            "An out-of-range contig2 index (5) is ignored, leaving two single-contig scaffolds.");
    }

    // S3 — A self-link (contig1 == contig2) is ignored.
    [Test]
    public void Scaffold_SelfLink_Ignored()
    {
        var contigs = new[] { "AA", "TT" };
        var links = new[] { (0, 0, 2) };

        var scaffolds = SequenceAssembler.Scaffold(contigs, links);

        Assert.That(scaffolds, Is.EqualTo(new[] { "AA", "TT" }),
            "A self-link cannot order a contig after itself and is ignored.");
    }

    // S4 — Null contigs must throw, matching sibling input-validation conventions (MergeContigs).
    [Test]
    public void Scaffold_NullContigs_ThrowsArgumentNullException()
    {
        var links = Array.Empty<(int, int, int)>();

        Assert.That(() => SequenceAssembler.Scaffold(null!, links),
            NUnit.Framework.Throws.TypeOf<ArgumentNullException>(),
            "Null contigs is invalid input and must throw ArgumentNullException.");
    }

    // S5 — Null links must throw, matching sibling input-validation conventions.
    [Test]
    public void Scaffold_NullLinks_ThrowsArgumentNullException()
    {
        var contigs = new[] { "AA" };

        Assert.That(() => SequenceAssembler.Scaffold(contigs, (IReadOnlyList<(int, int, int)>)null!),
            NUnit.Framework.Throws.TypeOf<ArgumentNullException>(),
            "Null links is invalid input and must throw ArgumentNullException.");
    }

    // S6 — A successor that is already placed elsewhere in the path is skipped, so a contig never
    // appears twice and the substrings stay verbatim in link order (INV-04/INV-05). Path
    // 0->1->2 is followed; the back-link 2->1 targets an already-placed contig and is dropped.
    // Source: Pop et al. (2004) path-of-distinct-contigs; Jackman et al. (2017) verbatim concatenation.
    [Test]
    public void Scaffold_SuccessorAlreadyPlaced_SkippedKeepingContigsDistinct()
    {
        var contigs = new[] { "AA", "CC", "GG" };
        var links = new[] { (0, 1, 2), (1, 2, 3), (2, 1, 4) };

        var scaffolds = SequenceAssembler.Scaffold(contigs, links);

        Assert.Multiple(() =>
        {
            Assert.That(scaffolds.Count, Is.EqualTo(1),
                "0->1->2 forms one path; the 2->1 back-link hits an already-placed contig (INV-04).");
            Assert.That(scaffolds[0], Is.EqualTo("AANNCCNNNGG"),
                "AA + 2 N + CC + 3 N + GG, each contig verbatim and used once (INV-05).");
            Assert.That(scaffolds[0].Length, Is.EqualTo(11),
                "Length = contigs (6) + gaps (2+3) = 11 (INV-03).");
        });
    }

    // S7 — When a contig has several forward links, the first declared link is followed (the unit
    // consumes pre-ordered links; input order is the documented tie-break). Following 0->1 ends the
    // path because contig 1 has no outgoing link, so the later 0->2 link is never taken and contig 2
    // stays its own length-1 scaffold. Source: Pop et al. (2004) path model; documented tie-break.
    [Test]
    public void Scaffold_MultipleForwardLinks_FollowsFirstDeclared()
    {
        var contigs = new[] { "AA", "CC", "GG" };
        var links = new[] { (0, 1, 2), (0, 2, 3) };

        var scaffolds = SequenceAssembler.Scaffold(contigs, links);

        Assert.Multiple(() =>
        {
            Assert.That(scaffolds, Is.EqualTo(new[] { "AANNCC", "GG" }),
                "First declared link 0->1 is followed (AA+2N+CC); contig 2 stays a length-1 scaffold.");
            Assert.That(scaffolds.Count, Is.EqualTo(2),
                "Two scaffolds: the 0->1 path and the unreached contig 2 (INV-04).");
        });
    }

    // C1 — An empty contig list yields an empty scaffold list (trivial identity).
    [Test]
    public void Scaffold_EmptyContigs_ReturnsEmpty()
    {
        var contigs = Array.Empty<string>();
        var links = Array.Empty<(int, int, int)>();

        var scaffolds = SequenceAssembler.Scaffold(contigs, links);

        Assert.That(scaffolds, Is.Empty,
            "No contigs means no scaffolds.");
    }

    #endregion

    #region Greedy path selection, orientation and overlaps (B10 review, F16-F18)

    // Expected values below were produced by the independent Python oracle recorded in
    // docs/Evidence/ASSEMBLY-SCAFFOLD-001-Evidence.md (Huson, Reinert & Myers 2002 Algorithm 5.1
    // greedy over contig end nodes; Biopython reverse complement; ABySS MergeContigs overlap rule;
    // NCBI AGP v2.1 unknown gap = 100).

    private static string Gap(int n) => new('N', n);

    // F16 — a link INTO a lower-index contig was dropped because scaffolds were started at every
    // unused index in ascending order (old: ["AA","TT"]).
    [Test]
    public void Scaffold_LinkIntoLowerIndexContig_IsJoined()
    {
        var scaffolds = SequenceAssembler.Scaffold(new[] { "AA", "TT" }, new[] { (1, 0, 2) });
        Assert.That(scaffolds, Is.EqualTo(new[] { "TTNNAA" }));
    }

    [Test]
    public void Scaffold_PredecessorOfFirstContig_StartsTheScaffold()
    {
        var scaffolds = SequenceAssembler.Scaffold(new[] { "AA", "CC", "GG" }, new[] { (2, 0, 1), (0, 1, 1) });
        Assert.That(scaffolds, Is.EqualTo(new[] { "GGNAANCC" }), "path 2→0→1 (old: AANCC, GG)");
    }

    // A second successor at an already-joined end is a conflict and is not used (old code
    // skipped the placed target and took the next link: CCNGG).
    [Test]
    public void Scaffold_SecondSuccessorOfJoinedEnd_IsRejected()
    {
        var scaffolds = SequenceAssembler.Scaffold(new[] { "AA", "CC", "GG" }, new[] { (1, 0, 1), (1, 2, 1) });
        Assert.That(scaffolds, Is.EqualTo(new[] { "CCNAA", "GG" }));
    }

    [Test]
    public void Scaffold_HeavierLinkWins_RegardlessOfInputOrder()
    {
        var contigs = new[] { "AAAA", "CCCC", "GGGG" };
        var links = new[]
        {
            new SequenceAssembler.ScaffoldLink(0, 1, 2, Weight: 1),
            new SequenceAssembler.ScaffoldLink(0, 2, 3, Weight: 5),
        };
        Assert.That(SequenceAssembler.Scaffold(contigs, links), Is.EqualTo(new[] { "AAAANNNGGGG", "CCCC" }),
            "mate-edges are taken in descending weight (Huson et al. 2002 Alg. 5.1)");
    }

    [Test]
    public void Scaffold_CycleClosingLink_IsRejected()
    {
        var scaffolds = SequenceAssembler.Scaffold(new[] { "AA", "CC", "GG" }, new[] { (2, 0, 1), (0, 1, 1), (1, 2, 1) });
        Assert.That(scaffolds, Is.EqualTo(new[] { "GGNAANCC" }), "1→2 would close the cycle 2→0→1→2");
    }

    [TestCase(false, true, "AACGNNNGCAA")]
    [TestCase(true, false, "GCAANNNAACG")]
    [TestCase(true, true, "TTGCNNNAACG")]
    public void Scaffold_OrientedLink_ReverseComplementsContigs(bool rev1, bool rev2, string expected)
    {
        var links = new[] { new SequenceAssembler.ScaffoldLink(0, 1, 3, rev1, rev2) };
        Assert.That(SequenceAssembler.Scaffold(new[] { "AACG", "TTGC" }, links), Is.EqualTo(new[] { expected }));
    }

    [Test]
    public void Scaffold_OrientedChain_FollowsBidirectedPath()
    {
        var contigs = new[] { "ACCA", "GGTT", "CATG" };
        var consistent = new[]
        {
            new SequenceAssembler.ScaffoldLink(0, 1, 2, Contig2Reversed: true),
            new SequenceAssembler.ScaffoldLink(1, 2, 4, Contig1Reversed: true),
        };
        var inconsistent = new[]
        {
            new SequenceAssembler.ScaffoldLink(0, 1, 2, Contig2Reversed: true),
            new SequenceAssembler.ScaffoldLink(1, 2, 4),
        };
        Assert.Multiple(() =>
        {
            Assert.That(SequenceAssembler.Scaffold(contigs, consistent), Is.EqualTo(new[] { "ACCANNAACCNNNNCATG" }));
            Assert.That(SequenceAssembler.Scaffold(contigs, inconsistent), Is.EqualTo(new[] { "ACCANNAACC", "CATG" }),
                "1+→2+ needs contig 1's 3′ end, already joined to contig 0");
        });
    }

    [TestCase("ACGTAC", "TACGGA", "ACGTACGGA")]
    [TestCase("ACGTNC", "TACGGA", "ACGTACGGA")]
    public void Scaffold_NegativeGapWithAgreeingOverlap_IsCollapsed(string a, string b, string expected)
    {
        Assert.That(SequenceAssembler.Scaffold(new[] { a, b }, new[] { (0, 1, -3) }), Is.EqualTo(new[] { expected }),
            "an overlap of −d agreeing bases is merged (Jackman et al. 2017; ABySS MergeContigs)");
    }

    [Test]
    public void Scaffold_NegativeGapWithMismatchedOverlap_EmitsUnknownGap()
    {
        Assert.That(SequenceAssembler.Scaffold(new[] { "ACGTAC", "TTCGGA" }, new[] { (0, 1, -3) }),
            Is.EqualTo(new[] { "ACGTAC" + Gap(100) + "TTCGGA" }));
    }

    [Test]
    public void Scaffold_NegativeGapOverlapOnReverseComplement_IsCollapsed()
    {
        var links = new[] { new SequenceAssembler.ScaffoldLink(0, 1, -2, Contig2Reversed: true) };
        Assert.That(SequenceAssembler.Scaffold(new[] { "GGATCC", "TTTGG" }, links), Is.EqualTo(new[] { "GGATCCAAA" }));
    }

    [Test]
    public void Scaffold_OverlapThenGap_Chain()
    {
        var scaffolds = SequenceAssembler.Scaffold(new[] { "AAACG", "CGTTT", "GGGG" }, new[] { (0, 1, -2), (1, 2, 3) });
        Assert.That(scaffolds, Is.EqualTo(new[] { "AAACGTTTNNNGGGG" }));
    }

    [Test]
    public void Scaffold_ScaffoldsOrderedByLowestContigIndex()
    {
        Assert.That(SequenceAssembler.Scaffold(new[] { "AAA", "CCC", "GGG" }, new[] { (2, 1, 4) }),
            Is.EqualTo(new[] { "AAA", "GGGNNNNCCC" }));
    }

    [Test]
    public void Scaffold_LinkWeightBelowOne_Throws()
    {
        var links = new[] { new SequenceAssembler.ScaffoldLink(0, 1, 3, Weight: 0) };
        Assert.Throws<ArgumentOutOfRangeException>(() => SequenceAssembler.Scaffold(new[] { "AA", "CC" }, links));
    }

    #endregion
}
