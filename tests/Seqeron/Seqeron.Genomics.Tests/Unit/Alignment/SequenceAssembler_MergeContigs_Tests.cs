// ASSEMBLY-MERGE-001 — Contig Merging (suffix–prefix overlap collapse)
// Evidence: docs/Evidence/ASSEMBLY-MERGE-001-Evidence.md
// TestSpec: tests/TestSpecs/ASSEMBLY-MERGE-001.md
// Source: Langmead B., "Assembly & shortest common superstring" / "Overlap Layout Consensus
//         assembly" (JHU lecture notes); MIT 7.91J Lecture 6 (MIT OCW, 2014).

using System;
using NUnit.Framework;

using Seqeron.Genomics.Alignment;

namespace Seqeron.Genomics.Tests.Unit.Alignment;

[TestFixture]
public class SequenceAssembler_MergeContigs_Tests
{
    #region MergeContigs

    // M1 — Published greedy-SCS merge: BAA + AAB, suffix "AA" of BAA == prefix "AA" of AAB
    // (overlap length 2) collapses to BAAB. Source: Langmead SCS notes greedy trace.
    [Test]
    public void MergeContigs_BaaAabOverlapTwo_ReturnsBaab()
    {
        string merged = SequenceAssembler.MergeContigs("BAA", "AAB", 2);

        Assert.Multiple(() =>
        {
            Assert.That(merged, Is.EqualTo("BAAB"),
                "BAA + AAB at overlap 2 keeps one copy of the shared 'AA' -> BAAB (Langmead SCS greedy trace).");
            Assert.That(merged.Length, Is.EqualTo(4),
                "Merged length = |BAA| + |AAB| - 2 = 3 + 3 - 2 = 4 (INV-01).");
        });
    }

    // M2 — Chaining the published overlaps of {AAA, AAB, ABB, BBB, BBA} (each overlap 2) yields
    // the shortest common superstring AAABBBA. Source: Langmead SCS notes worked example.
    [Test]
    public void MergeContigs_ChainAaaAabAbbBbbBba_ReturnsAaabbba()
    {
        // Build the superstring by successively collapsing the next string at overlap 2.
        string s = SequenceAssembler.MergeContigs("AAA", "AAB", 2);   // AAAB
        s = SequenceAssembler.MergeContigs(s, "ABB", 2);              // AAABB
        s = SequenceAssembler.MergeContigs(s, "BBB", 2);              // AAABBB
        s = SequenceAssembler.MergeContigs(s, "BBA", 2);              // AAABBBA

        Assert.Multiple(() =>
        {
            Assert.That(s, Is.EqualTo("AAABBBA"),
                "Chaining {AAA,AAB,ABB,BBB,BBA} at overlap 2 reconstructs the SCS AAABBBA (Langmead SCS notes).");
            Assert.That(s.Length, Is.EqualTo(7),
                "SCS length is 7 per Langmead SCS worked example.");
        });
    }

    // M3 — Overlap length 0 means no usable overlap, so the contigs are simply concatenated.
    // Source: Langmead SCS notes ("without requirement of 'shortest' ... just concatenate them").
    [Test]
    public void MergeContigs_OverlapZero_ConcatenatesContigs()
    {
        string merged = SequenceAssembler.MergeContigs("BAA", "AAB", 0);

        Assert.Multiple(() =>
        {
            Assert.That(merged, Is.EqualTo("BAAAAB"),
                "Overlap 0 = no overlap, so BAA and AAB are concatenated verbatim (Langmead SCS).");
            Assert.That(merged.Length, Is.EqualTo(6),
                "No overlap removed: length = |BAA| + |AAB| = 6 (INV-02).");
        });
    }

    // M4 — Boundary: overlap equals the shorter string's length (here 3 = |ACGT|? no, = |CGT...|).
    // ACGT suffix "CGT" == CGTAA prefix "CGT", overlap 3 -> ACGTAA. Source: overlap bounded by
    // min length and collapsed to one copy (Langmead SCS suffixPrefixMatch + merge).
    [Test]
    public void MergeContigs_OverlapEqualsMinLength_CollapsesSharedPrefix()
    {
        // contig2 prefix entirely consumed by the overlap (overlap 3 == |"CGT" shared region|).
        string merged = SequenceAssembler.MergeContigs("ACGT", "CGTAA", 3);

        Assert.Multiple(() =>
        {
            Assert.That(merged, Is.EqualTo("ACGTAA"),
                "ACGT + CGTAA at overlap 3 keeps one 'CGT' -> ACGTAA (Langmead SCS merge).");
            Assert.That(merged.Length, Is.EqualTo(6),
                "Length = 4 + 5 - 3 = 6 (INV-01).");
        });
    }

    // M5 — Length invariant for a valid overlap: GATTACA + ACATGAA, suffix "ACA" == prefix "ACA",
    // overlap 3 -> GATTACATGAA, length 7 + 7 - 3 = 11. Source: INV-01 (Langmead SCS).
    [Test]
    public void MergeContigs_ValidOverlap_LengthEqualsSumMinusOverlap()
    {
        string merged = SequenceAssembler.MergeContigs("GATTACA", "ACATGAA", 3);

        Assert.Multiple(() =>
        {
            Assert.That(merged, Is.EqualTo("GATTACATGAA"),
                "GATTACA + ACATGAA at overlap 3 keeps one 'ACA' -> GATTACATGAA (Langmead SCS merge).");
            Assert.That(merged.Length, Is.EqualTo(11),
                "Length = |GATTACA| + |ACATGAA| - 3 = 7 + 7 - 3 = 11 (INV-01).");
        });
    }

    // S1 — An overlap larger than the shorter contig is not a valid suffix/prefix overlap, so the
    // contigs are concatenated. Source: overlap bounded by min(|x|,|y|) (Langmead SCS).
    [Test]
    public void MergeContigs_OverlapExceedsMinLength_ConcatenatesContigs()
    {
        // min(|"AC"|, |"GTAA"|) = 2; overlap 3 > 2 is invalid.
        string merged = SequenceAssembler.MergeContigs("AC", "GTAA", 3);

        Assert.That(merged, Is.EqualTo("ACGTAA"),
            "Overlap 3 exceeds min length 2, so it is not a valid overlap -> concatenation ACGTAA (INV-03).");
    }

    // S2 — A negative overlap is non-positive, i.e. no overlap, so the contigs are concatenated.
    // Source: overlap length 0 = concatenation; non-positive is likewise no overlap (Langmead SCS).
    [Test]
    public void MergeContigs_NegativeOverlap_ConcatenatesContigs()
    {
        string merged = SequenceAssembler.MergeContigs("BAA", "AAB", -2);

        Assert.That(merged, Is.EqualTo("BAAAAB"),
            "Negative overlap is not a valid overlap, so BAA and AAB are concatenated (INV-03).");
    }

    // S3 — Null contig1 must throw, matching sibling input-validation conventions.
    [Test]
    public void MergeContigs_NullContig1_ThrowsArgumentNullException()
    {
        Assert.That(() => SequenceAssembler.MergeContigs(null!, "AAB", 1),
            NUnit.Framework.Throws.TypeOf<ArgumentNullException>(),
            "Null contig1 is invalid input and must throw ArgumentNullException.");
    }

    // S4 — Null contig2 must throw, matching sibling input-validation conventions.
    [Test]
    public void MergeContigs_NullContig2_ThrowsArgumentNullException()
    {
        Assert.That(() => SequenceAssembler.MergeContigs("BAA", null!, 1),
            NUnit.Framework.Throws.TypeOf<ArgumentNullException>(),
            "Null contig2 is invalid input and must throw ArgumentNullException.");
    }

    // C1 — Empty left contig with overlap 0 is the identity: result is the right contig.
    [Test]
    public void MergeContigs_EmptyContig1OverlapZero_ReturnsContig2()
    {
        string merged = SequenceAssembler.MergeContigs("", "AAB", 0);

        Assert.That(merged, Is.EqualTo("AAB"),
            "Empty + AAB with no overlap is the identity -> AAB.");
    }

    // C2 — Empty right contig with overlap 0 is the identity: result is the left contig.
    [Test]
    public void MergeContigs_EmptyContig2OverlapZero_ReturnsContig1()
    {
        string merged = SequenceAssembler.MergeContigs("BAA", "", 0);

        Assert.That(merged, Is.EqualTo("BAA"),
            "BAA + empty with no overlap is the identity -> BAA.");
    }

    // C3 — Containment property (INV-04): for a valid overlap the result starts with contig1 and
    // ends with the non-overlapped tail of contig2.
    [Test]
    public void MergeContigs_ValidOverlap_StartsWithContig1AndEndsWithContig2Tail()
    {
        const string c1 = "GATTACA";
        const string c2 = "ACATGAA";
        const int overlap = 3;

        string merged = SequenceAssembler.MergeContigs(c1, c2, overlap);

        Assert.Multiple(() =>
        {
            Assert.That(merged.StartsWith(c1, StringComparison.Ordinal), Is.True,
                "INV-04: contig1 is preserved verbatim as a prefix of the merged superstring.");
            Assert.That(merged.EndsWith(c2.Substring(overlap), StringComparison.Ordinal), Is.True,
                "INV-04: the non-overlapped tail of contig2 is the suffix of the merged superstring.");
        });
    }

    #endregion
}
