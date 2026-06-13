// ASSEMBLY-SCAFFOLD-001 — Scaffolding (joining ordered contigs with N-gaps)
// Evidence: docs/Evidence/ASSEMBLY-SCAFFOLD-001-Evidence.md
// TestSpec: tests/TestSpecs/ASSEMBLY-SCAFFOLD-001.md
// Source: Jackman SD et al. (2017). ABySS 2.0. Genome Research 27:768-777
//         (scaffold construction); NCBI AGP Specification v2.1 (unknown-gap = 100 N).

using System;
using System.Linq;
using NUnit.Framework;

using Seqeron.Genomics.Alignment;

namespace Seqeron.Genomics.Tests;

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

        Assert.That(() => SequenceAssembler.Scaffold(contigs, null!),
            NUnit.Framework.Throws.TypeOf<ArgumentNullException>(),
            "Null links is invalid input and must throw ArgumentNullException.");
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
}
