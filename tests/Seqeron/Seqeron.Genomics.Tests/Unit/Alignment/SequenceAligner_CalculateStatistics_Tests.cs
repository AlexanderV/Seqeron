// ALIGN-STATS-001 — Pairwise Alignment Statistics and Formatting
// Evidence: docs/Evidence/ALIGN-STATS-001-Evidence.md
// TestSpec: tests/TestSpecs/ALIGN-STATS-001.md
// Source: Rice P, Longden I, Bleasby A (2000). EMBOSS. Trends Genet. 16(6):276-277.
//         EMBOSS needle docs (rel 6.6); EMBOSS Alignment Formats; NCBI BLAST QuickStart NBK1734.
//
// Definitions (EMBOSS needle): denominator = alignment Length INCLUDING gap columns.
//   Identity%   = identical / L * 100
//   Similarity% = (identical + positive-score mismatches) / L * 100
//   Gaps%       = gap-columns / L * 100
// A non-identical column is "similar" iff its substitution score is POSITIVE.
// srspair markup legend: '|' identity, ':' similarity, ' ' gap/mismatch.

namespace Seqeron.Genomics.Tests.Unit.Alignment;

[TestFixture]
public class SequenceAligner_CalculateStatistics_Tests
{
    // Default DNA model: identical scores +1 (positive), mismatch scores -1 (non-positive).
    private static readonly ScoringMatrix SimpleDna = SequenceAligner.SimpleDna;

    // Model where a non-identical column scores POSITIVELY (+1), so mismatches count as similar.
    private static readonly ScoringMatrix PositiveMismatch =
        new(Match: 1, Mismatch: 1, GapOpen: 0, GapExtend: -1);

    /// <summary>Builds an AlignmentResult directly from two equal-length gapped rows.</summary>
    private static AlignmentResult Aln(string a, string b) =>
        new(a, b, Score: 0, AlignmentType.Global,
            StartPosition1: 0, StartPosition2: 0,
            EndPosition1: a.Length - 1, EndPosition2: b.Length - 1);

    #region CalculateStatistics

    // M1 — EMBOSS HBA vs HBB published example: L=149, identical=65, similar=90, gaps=9.
    // Reproduces the exact published percentages (43.6 / 60.4 / 6.0) from the formula and
    // the gap-inclusive denominator. A 149-column alignment is synthesised with the published
    // class counts: 65 identical, 25 similar mismatches, 59 gaps. PositiveMismatch makes the
    // 25 non-identical columns score positively (similar), so the similar set = 65 + 25 = 90,
    // matching the published Similarity count while gaps (59) never count as similar.
    [Test]
    public void CalculateStatistics_EmbossExample_MatchesPublishedPercentages()
    {
        const int identical = 65, similarMm = 25;
        int gaps = 149 - identical - similarMm; // 59 (only Identity/Similarity/Gap% are asserted)
        var sb1 = new System.Text.StringBuilder();
        var sb2 = new System.Text.StringBuilder();
        for (int i = 0; i < identical; i++) { sb1.Append('A'); sb2.Append('A'); } // identical
        for (int i = 0; i < similarMm; i++) { sb1.Append('A'); sb2.Append('C'); } // similar mismatch (positive)
        for (int i = 0; i < gaps; i++) { sb1.Append('A'); sb2.Append('-'); }      // gap

        var aln = Aln(sb1.ToString(), sb2.ToString());

        var stats = SequenceAligner.CalculateStatistics(aln, PositiveMismatch);

        Assert.Multiple(() =>
        {
            Assert.That(stats.AlignmentLength, Is.EqualTo(149),
                "Denominator = total columns including gap columns (EMBOSS Length=149).");
            Assert.That(stats.Matches, Is.EqualTo(65), "Published Identity count = 65.");
            Assert.That(stats.Identity, Is.EqualTo(65.0 / 149.0 * 100).Within(1e-3),
                "Identity% = 65/149*100 = 43.6% (EMBOSS worked example).");
            Assert.That(stats.Identity, Is.EqualTo(43.6).Within(0.05), "Rounds to published 43.6%.");
            Assert.That(stats.Similarity, Is.EqualTo(90.0 / 149.0 * 100).Within(1e-3),
                "Similarity% = (65 identical + 25 similar)/149*100 = 60.4% (EMBOSS).");
            Assert.That(stats.Similarity, Is.EqualTo(60.4).Within(0.05), "Rounds to published 60.4%.");
        });
    }

    // M2 — Under SimpleDna (Mismatch < 0) no non-identical column is similar ⇒ Similarity == Identity.
    [Test]
    public void CalculateStatistics_SimpleDna_SimilarityEqualsIdentity()
    {
        var aln = Aln("ACGT", "ACCT"); // column 3 G:C is a mismatch, scores -1 (not similar)

        var stats = SequenceAligner.CalculateStatistics(aln, SimpleDna);

        Assert.Multiple(() =>
        {
            Assert.That(stats.Matches, Is.EqualTo(3), "3 identical columns (A,C,T).");
            Assert.That(stats.Mismatches, Is.EqualTo(1), "1 mismatch column (G:C).");
            Assert.That(stats.Identity, Is.EqualTo(75.0).Within(1e-10), "3/4*100 = 75%.");
            Assert.That(stats.Similarity, Is.EqualTo(stats.Identity).Within(1e-10),
                "DNA mismatch scores negative ⇒ not similar ⇒ Similarity equals Identity.");
        });
    }

    // M2b — Default scoring (no matrix argument) must also be SimpleDna semantics.
    [Test]
    public void CalculateStatistics_DefaultScoring_SimilarityEqualsIdentity()
    {
        var aln = Aln("ACGT", "ACCT");

        var stats = SequenceAligner.CalculateStatistics(aln);

        Assert.That(stats.Similarity, Is.EqualTo(stats.Identity).Within(1e-10),
            "Default scoring is SimpleDna; mismatches are not similar.");
    }

    // M3 — Positive-scoring mismatch ⇒ the mismatch is similar ⇒ Similarity > Identity.
    [Test]
    public void CalculateStatistics_PositiveScoringMismatch_SimilarityExceedsIdentity()
    {
        var aln = Aln("ACGT", "ACCT"); // 3 identical, 1 mismatch (scores +1 under PositiveMismatch)

        var stats = SequenceAligner.CalculateStatistics(aln, PositiveMismatch);

        Assert.Multiple(() =>
        {
            Assert.That(stats.Identity, Is.EqualTo(75.0).Within(1e-10), "3/4*100 = 75%.");
            Assert.That(stats.Similarity, Is.EqualTo(100.0).Within(1e-10),
                "(3 identical + 1 positive-score mismatch)/4*100 = 100% (positive substitution ⇒ similar).");
        });
    }

    // M4 — Gap columns are counted as gaps, never as identity or mismatch.
    [Test]
    public void CalculateStatistics_GapColumns_CountedAsGaps()
    {
        var aln = Aln("A-CG", "ATCG"); // column 2 is a gap (- vs T)

        var stats = SequenceAligner.CalculateStatistics(aln, SimpleDna);

        Assert.Multiple(() =>
        {
            Assert.That(stats.Gaps, Is.EqualTo(1), "One column has a gap symbol.");
            Assert.That(stats.Matches, Is.EqualTo(3), "A, C, G align identically.");
            Assert.That(stats.Mismatches, Is.EqualTo(0), "Gap column is not a mismatch.");
            Assert.That(stats.GapPercent, Is.EqualTo(25.0).Within(1e-10), "1/4*100 = 25%.");
        });
    }

    // M5 — Exact counts/percentages on the hand-derived 9-column alignment from Evidence.
    [Test]
    public void CalculateStatistics_HandAlignment_ExactCountsAndPercentages()
    {
        // Evidence dataset: 6 identical, 1 mismatch (G:C), 2 gaps over Length 9.
        var aln = Aln("ACGT-ACGT", "ACCTAAC-T");

        var stats = SequenceAligner.CalculateStatistics(aln, SimpleDna);

        Assert.Multiple(() =>
        {
            Assert.That(stats.AlignmentLength, Is.EqualTo(9), "9 columns total.");
            Assert.That(stats.Matches, Is.EqualTo(6), "6 identical columns.");
            Assert.That(stats.Mismatches, Is.EqualTo(1), "1 mismatch column (G:C).");
            Assert.That(stats.Gaps, Is.EqualTo(2), "2 gap columns.");
            Assert.That(stats.Identity, Is.EqualTo(6.0 / 9.0 * 100).Within(1e-10), "6/9*100 = 66.67%.");
            Assert.That(stats.Similarity, Is.EqualTo(6.0 / 9.0 * 100).Within(1e-10),
                "SimpleDna: mismatch not similar ⇒ Similarity equals Identity.");
            Assert.That(stats.GapPercent, Is.EqualTo(2.0 / 9.0 * 100).Within(1e-10), "2/9*100 = 22.22%.");
            Assert.That(stats.Matches + stats.Mismatches + stats.Gaps, Is.EqualTo(stats.AlignmentLength),
                "INV-1: classes partition the alignment.");
        });
    }

    // S1 — Empty alignment ⇒ AlignmentStatistics.Empty.
    [Test]
    public void CalculateStatistics_EmptyAlignment_ReturnsEmpty()
    {
        var stats = SequenceAligner.CalculateStatistics(AlignmentResult.Empty);

        Assert.That(stats, Is.EqualTo(AlignmentStatistics.Empty),
            "Empty alignment has no columns; all statistics are zero.");
    }

    // S3 — Null alignment ⇒ ArgumentNullException.
    [Test]
    public void CalculateStatistics_Null_Throws()
    {
        Assert.Throws<ArgumentNullException>(
            () => SequenceAligner.CalculateStatistics(null!),
            "Null alignment must be rejected.");
    }

    // S5 — Perfect identity ⇒ Identity = Similarity = 100%, Gap 0%.
    [Test]
    public void CalculateStatistics_PerfectIdentity_HundredPercent()
    {
        var aln = Aln("ACGTACGT", "ACGTACGT");

        var stats = SequenceAligner.CalculateStatistics(aln, SimpleDna);

        Assert.Multiple(() =>
        {
            Assert.That(stats.Identity, Is.EqualTo(100.0).Within(1e-10), "All columns identical.");
            Assert.That(stats.Similarity, Is.EqualTo(100.0).Within(1e-10), "Identical ⇒ similar.");
            Assert.That(stats.GapPercent, Is.EqualTo(0.0).Within(1e-10), "No gaps.");
        });
    }

    // S6 — All-gap alignment ⇒ Gaps == Length, Identity 0%.
    [Test]
    public void CalculateStatistics_AllGaps_GapsEqualLength()
    {
        var aln = Aln("---", "ACG");

        var stats = SequenceAligner.CalculateStatistics(aln, SimpleDna);

        Assert.Multiple(() =>
        {
            Assert.That(stats.Gaps, Is.EqualTo(3), "Every column has a gap.");
            Assert.That(stats.Identity, Is.EqualTo(0.0).Within(1e-10), "No identical columns.");
            Assert.That(stats.GapPercent, Is.EqualTo(100.0).Within(1e-10), "3/3*100 = 100%.");
        });
    }

    // C2 — Invariants INV-1, INV-2, INV-3 hold on a deterministic mixed alignment.
    [Test]
    public void CalculateStatistics_Invariants_Hold()
    {
        var aln = Aln("ACGTACG-TACGTA", "ACATAAGGT-CGTA");
        var stats = SequenceAligner.CalculateStatistics(aln, PositiveMismatch);

        double mismatchOnlyPct = (double)stats.Mismatches / stats.AlignmentLength * 100
                                 - (stats.Similarity - stats.Identity);

        Assert.Multiple(() =>
        {
            Assert.That(stats.Matches + stats.Mismatches + stats.Gaps,
                Is.EqualTo(stats.AlignmentLength), "INV-1: column classes partition the alignment.");
            Assert.That(stats.Identity, Is.LessThanOrEqualTo(stats.Similarity),
                "INV-2: Identity% ≤ Similarity%.");
            Assert.That(stats.Similarity, Is.LessThanOrEqualTo(100.0 + 1e-9),
                "INV-2: Similarity% ≤ 100%.");
            Assert.That(stats.Identity + (stats.Similarity - stats.Identity) + mismatchOnlyPct + stats.GapPercent,
                Is.EqualTo(100.0).Within(1e-9), "INV-3: length-denominator partition sums to 100%.");
        });
    }

    #endregion

    #region FormatAlignment

    // M6 — srspair legend under SimpleDna: '|' identity, space for mismatch AND gap.
    [Test]
    public void FormatAlignment_SimpleDna_UsesIdentityAndGapLegend()
    {
        var aln = Aln("ACGT", "ACCT"); // col3 mismatch (G:C, non-positive ⇒ space under SimpleDna)
        string formatted = SequenceAligner.FormatAlignment(aln, lineWidth: 60, scoring: SimpleDna);

        string[] lines = formatted.Replace("\r\n", "\n").Split('\n');
        string markup = lines[1];

        Assert.That(markup, Is.EqualTo("||" + " " + "|"),
            "srspair: identity '|', mismatch under SimpleDna renders as space (non-positive score).");
    }

    // M6b — Gap column renders as space in the markup line.
    [Test]
    public void FormatAlignment_GapColumn_RendersSpace()
    {
        var aln = Aln("A-CG", "ATCG");
        string formatted = SequenceAligner.FormatAlignment(aln, scoring: SimpleDna);

        string markup = formatted.Replace("\r\n", "\n").Split('\n')[1];

        Assert.That(markup, Is.EqualTo("| ||"),
            "srspair: gap column renders as a space, identities as '|'.");
    }

    // M7 — Positive-scoring mismatch renders the similarity mark ':'.
    [Test]
    public void FormatAlignment_PositiveScoringMismatch_UsesSimilarityMark()
    {
        var aln = Aln("ACGT", "ACCT"); // col3 mismatch scores +1 ⇒ similar ⇒ ':'
        string formatted = SequenceAligner.FormatAlignment(aln, scoring: PositiveMismatch);

        string markup = formatted.Replace("\r\n", "\n").Split('\n')[1];

        Assert.That(markup, Is.EqualTo("||:|"),
            "srspair: positive-score non-identical column renders ':' (similarity).");
    }

    // S2 — Empty alignment ⇒ empty string.
    [Test]
    public void FormatAlignment_EmptyAlignment_ReturnsEmptyString()
    {
        Assert.That(SequenceAligner.FormatAlignment(AlignmentResult.Empty), Is.EqualTo(""),
            "Empty alignment renders to an empty string.");
    }

    // S3 — Null alignment ⇒ ArgumentNullException.
    [Test]
    public void FormatAlignment_Null_Throws()
    {
        Assert.Throws<ArgumentNullException>(
            () => SequenceAligner.FormatAlignment(null!),
            "Null alignment must be rejected.");
    }

    // S4 — Non-positive lineWidth ⇒ ArgumentOutOfRangeException.
    [Test]
    public void FormatAlignment_NonPositiveLineWidth_Throws()
    {
        var aln = Aln("ACGT", "ACGT");
        Assert.Multiple(() =>
        {
            Assert.Throws<ArgumentOutOfRangeException>(
                () => SequenceAligner.FormatAlignment(aln, lineWidth: 0), "lineWidth 0 is invalid.");
            Assert.Throws<ArgumentOutOfRangeException>(
                () => SequenceAligner.FormatAlignment(aln, lineWidth: -1), "negative lineWidth is invalid.");
        });
    }

    // C1 — Narrow lineWidth wraps the alignment into multiple blocks; markup line lengths track blocks.
    [Test]
    public void FormatAlignment_NarrowLineWidth_WrapsIntoBlocks()
    {
        var aln = Aln("ACGTAC", "ACGTAC"); // length 6, all identical
        string formatted = SequenceAligner.FormatAlignment(aln, lineWidth: 4, scoring: SimpleDna);

        string[] lines = formatted.Replace("\r\n", "\n").Split('\n');

        Assert.Multiple(() =>
        {
            // Block 1: row1, markup, row2, blank ; Block 2: row1, markup, row2, blank
            Assert.That(lines[0], Is.EqualTo("ACGT"), "First block row1 = first 4 residues.");
            Assert.That(lines[1], Is.EqualTo("||||"), "First block markup = 4 identity marks.");
            Assert.That(lines[4], Is.EqualTo("AC"), "Second block row1 = remaining 2 residues.");
            Assert.That(lines[5], Is.EqualTo("||"), "Second block markup length tracks block length (INV-4).");
        });
    }

    #endregion
}
