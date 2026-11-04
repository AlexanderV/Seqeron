// ALIGN-STATS-001 — Pairwise Alignment Statistics and Formatting
// Evidence: docs/Evidence/ALIGN-STATS-001-Evidence.md
// TestSpec: tests/TestSpecs/ALIGN-STATS-001.md
// Source: Rice P, Longden I, Bleasby A (2000). EMBOSS. Trends Genet. 16(6):276-277.
//         EMBOSS needle docs (rel 6.6); EMBOSS Alignment Formats; NCBI BLAST QuickStart NBK1734;
//         EMBOSS 6.6.0 source ajax/core/ajalign.c (alignConsStats, alignSim) and the needle/water 6.6.0 binary.
//
// Definitions (EMBOSS needle): denominator = alignment Length INCLUDING gap columns.
//   Identity%   = identical / L * 100
//   Similarity% = columns with a POSITIVE substitution score / L * 100 (BLAST "positives")
//   Gaps%       = columns with a gap in either row / L * 100
// srspair markup (needle 6.6.0 output): '|' identity, ':' similar (score > 0), '.' mismatch, ' ' gap column.
// Oracle table TestData/Alignment/needle_water_stats.tsv: EMBOSS needle/water 6.6.0 header + markup and
// Biopython 1.88 Alignment.counts(matrix) for 301 alignments (scripts/alignment/generate_alignment_stats_oracle.py).

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

    // M1 (formula) — the denominator/percent arithmetic of the EMBOSS HBA vs HBB example (L=149,
    // identical=65, similar=90). A 149-column alignment is synthesised with 65 identical, 25 similar
    // mismatches and 59 gap columns (NOT the real gap count: the literal published alignment is M1b).
    [Test]
    public void CalculateStatistics_EmbossFormula_MatchesPublishedPercentages()
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

    // M6 — srspair legend under SimpleDna: '|' identity, '.' mismatch, space for a gap column (needle 6.6.0).
    [Test]
    public void FormatAlignment_SimpleDna_UsesIdentityAndGapLegend()
    {
        var aln = Aln("ACGT", "ACCT"); // col3 mismatch (G:C, non-positive ⇒ '.' under SimpleDna)
        string formatted = SequenceAligner.FormatAlignment(aln, lineWidth: 60, scoring: SimpleDna);

        string[] lines = formatted.Replace("\r\n", "\n").Split('\n');
        string markup = lines[1];

        Assert.That(markup, Is.EqualTo("||.|"),
            "needle 6.6.0 srspair: identity '|', non-positive mismatch '.' (e.g. A/T in EDNAFULL), gap ' '.");
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

    #region EMBOSS / Biopython oracle (ALIGN-STATS-001 B09 F11-F14)

    private const string HbaRow =
        "MV-LSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHF-DLS-----HGSAQVKGHGKKVADALTNAVAHVDDMPNALSALSDLHAHKLRVDPVNFKLLSHCLLVTLAAHLPAEFTPAVHASLDKFLASVSTVLTSKYR";
    private const string HbbRow =
        "MVHLTPEEKSAVTALWGKV--NVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAVMGNPKVKAHGKKVLGAFSDGLAHLDNLKGTFATLSELHCDKLHVDPENFRLLGNVLVCVLAHHFGKEFTPPVQAAYQKVVAGVANALAHKYH";

    // M1b — the literal published EMBOSS needle example (HBA_HUMAN vs HBB_HUMAN, EBLOSUM62, gap open 10,
    // extend 0.5): Length 149, Identity 65/149 (43.6%), Similarity 90/149 (60.4%), Gaps 9/149 (6.0%),
    // reproduced with needle 6.6.0 (rows above are its output). Needs the real BLOSUM62 matrix.
    [Test]
    public void CalculateStatistics_EmbossHbaHbb_Blosum62_ReproducesPublishedExample()
    {
        var stats = SequenceAligner.CalculateStatistics(Aln(HbaRow, HbbRow), SubstitutionMatrix.Blosum62);

        Assert.Multiple(() =>
        {
            Assert.That(stats.AlignmentLength, Is.EqualTo(149));
            Assert.That(stats.Matches, Is.EqualTo(65), "Identity 65/149");
            Assert.That(stats.Positives, Is.EqualTo(90), "Similarity 90/149");
            Assert.That(stats.Gaps, Is.EqualTo(9), "Gaps 9/149");
            Assert.That(stats.Mismatches, Is.EqualTo(149 - 65 - 9));
            Assert.That(stats.Identity, Is.EqualTo(43.6).Within(0.05));
            Assert.That(stats.Similarity, Is.EqualTo(60.4).Within(0.05));
            Assert.That(stats.GapPercent, Is.EqualTo(6.0).Within(0.05));
            // Biopython 1.88 counts(BLOSUM62): open_gaps = 4 (HBA '-' at 3, 49, 55-59; HBB '--' at 20-21).
            Assert.That(stats.GapOpens, Is.EqualTo(4));
        });

        // The scalar DNA model cannot express BLOSUM62: only identities are positive under it.
        Assert.That(SequenceAligner.CalculateStatistics(Aln(HbaRow, HbbRow)).Positives, Is.EqualTo(65));
    }

    // M7b — the same alignment rendered like needle's srspair block (first 50 columns), BLOSUM62.
    [Test]
    public void FormatAlignment_EmbossHbaHbb_Blosum62_MatchesNeedleMarkupLine()
    {
        string formatted = SequenceAligner.FormatAlignment(Aln(HbaRow, HbbRow), 50, SubstitutionMatrix.Blosum62);

        string[] lines = formatted.Split('\n');
        Assert.Multiple(() =>
        {
            Assert.That(lines[0], Is.EqualTo("MV-LSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHF-D"));
            Assert.That(lines[1], Is.EqualTo("|| |:|.:|:.|.|.||||  :..|.|.|||.|:.:.:|.|:.:|..| |"), "needle 6.6.0 markup, block 1");
            Assert.That(lines[2], Is.EqualTo("MVHLTPEEKSAVTALWGKV--NVDEVGGEALGRLLVVYPWTQRFFESFGD"));
            Assert.That(lines[5], Is.EqualTo("||     .|:.:||.|||||..|.::.:||:|::....:.||:||..||."), "needle 6.6.0 markup, block 2");
            Assert.That(lines[9], Is.EqualTo("|||.||:||.:.|:..||.|...||||.|.|:..|.:|.|:..|..||."), "needle 6.6.0 markup, block 3");
        });
    }

    private sealed record OracleRow(
        int Id, string Program, string Matrix, string Row1, string Row2, string Markup,
        int Length, int Identity, int Similarity, int Gaps, int BpMismatches, int BpOpenGaps);

    private static List<OracleRow> ReadOracle()
    {
        const string name = "Seqeron.Genomics.Tests.TestData.Alignment.needle_water_stats.tsv";
        using var stream = System.Reflection.Assembly.GetExecutingAssembly().GetManifestResourceStream(name)
            ?? throw new InvalidOperationException($"Embedded oracle '{name}' not found.");
        using var reader = new StreamReader(stream);
        var rows = new List<OracleRow>();
        string? line;
        while ((line = reader.ReadLine()) is not null)
        {
            if (line.Length == 0 || line[0] == '#') continue;
            var f = line.Split('\t');
            rows.Add(new OracleRow(int.Parse(f[0]), f[1], f[2], f[3], f[4], f[5],
                int.Parse(f[6]), int.Parse(f[7]), int.Parse(f[8]), int.Parse(f[9]), int.Parse(f[10]), int.Parse(f[11])));
        }
        return rows;
    }

    private static SubstitutionMatrix MatrixOf(string emboss) => emboss switch
    {
        "EBLOSUM62" => SubstitutionMatrix.Blosum62,
        "EPAM250" => SubstitutionMatrix.Pam250,
        "EDNAFULL" => SubstitutionMatrix.EdnaFull,
        _ => throw new ArgumentOutOfRangeException(nameof(emboss))
    };

    // Every one of the 301 oracle alignments (needle + water 6.6.0; DNA, IUPAC DNA with EDNAFULL positives,
    // protein with BLOSUM62/PAM250 incl. B/Z/X): Length, Identity, Similarity, Gaps equal the EMBOSS header,
    // Mismatches and GapOpens equal Biopython 1.88 Alignment.counts, the markup equals needle's markup line.
    [Test]
    public void CalculateStatistics_And_FormatAlignment_MatchEmbossAndBiopython_OnOracleTable()
    {
        var rows = ReadOracle();
        Assert.That(rows, Has.Count.GreaterThanOrEqualTo(290));

        int iupacSimilarGtIdentity = 0;
        foreach (var r in rows)
        {
            var aln = Aln(r.Row1, r.Row2);
            var st = SequenceAligner.CalculateStatistics(aln, MatrixOf(r.Matrix));
            string ctx = $"#{r.Id} {r.Program} {r.Matrix}";

            Assert.That(st.AlignmentLength, Is.EqualTo(r.Length), ctx);
            Assert.That(st.Matches, Is.EqualTo(r.Identity), ctx + " identity");
            Assert.That(st.Positives, Is.EqualTo(r.Similarity), ctx + " similarity");
            Assert.That(st.Gaps, Is.EqualTo(r.Gaps), ctx + " gaps");
            Assert.That(st.Mismatches, Is.EqualTo(r.BpMismatches), ctx + " mismatches (Biopython)");
            Assert.That(st.GapOpens, Is.EqualTo(r.BpOpenGaps), ctx + " open_gaps (Biopython)");
            Assert.That(st.Matches + st.Mismatches + st.Gaps, Is.EqualTo(st.AlignmentLength), ctx);
            Assert.That(st.Identity, Is.EqualTo(100.0 * r.Identity / r.Length).Within(1e-9), ctx);
            Assert.That(st.Similarity, Is.EqualTo(100.0 * r.Similarity / r.Length).Within(1e-9), ctx);
            Assert.That(st.GapPercent, Is.EqualTo(100.0 * r.Gaps / r.Length).Within(1e-9), ctx);

            string formatted = SequenceAligner.FormatAlignment(aln, r.Length, MatrixOf(r.Matrix));
            Assert.That(formatted.Split('\n')[1], Is.EqualTo(r.Markup), ctx + " markup line");

            // Lower-case input is the same alignment (EMBOSS and the aligners are case-insensitive).
            var lower = SequenceAligner.CalculateStatistics(Aln(r.Row1.ToLowerInvariant(), r.Row2.ToLowerInvariant()), MatrixOf(r.Matrix));
            Assert.That(lower, Is.EqualTo(st), ctx + " lower-case");

            if (r.Matrix == "EDNAFULL" && r.Similarity > r.Identity) iupacSimilarGtIdentity++;
        }

        Assert.That(iupacSimilarGtIdentity, Is.GreaterThan(5), "oracle covers IUPAC DNA where similarity > identity");
    }

    // F11 — Residues are compared case-insensitively (EMBOSS; the aligners upper-case their input).
    [Test]
    public void CalculateStatistics_LowerCaseRows_AreCaseInsensitive()
    {
        var st = SequenceAligner.CalculateStatistics(Aln("acgt-a", "ACGTGa"));
        Assert.Multiple(() =>
        {
            Assert.That(st.Matches, Is.EqualTo(5));
            Assert.That(st.Mismatches, Is.EqualTo(0));
            Assert.That(st.Gaps, Is.EqualTo(1));
        });
        Assert.That(SequenceAligner.FormatAlignment(Aln("acgt", "ACGA")).Split('\n')[1], Is.EqualTo("|||."));
    }

    // F12 — rows of unequal length are rejected (was IndexOutOfRange / silently ignored tail).
    [Test]
    public void CalculateStatistics_And_FormatAlignment_UnequalRows_Throw()
    {
        var bad = Aln("ACGT", "AC");
        Assert.Throws<ArgumentException>(() => SequenceAligner.CalculateStatistics(bad));
        Assert.Throws<ArgumentException>(() => SequenceAligner.FormatAlignment(bad));
        var longer = Aln("AC", "ACGT");
        Assert.Throws<ArgumentException>(() => SequenceAligner.CalculateStatistics(longer));
    }

    // F13 — gap openings (BLAST tabular gapopen, Biopython open_gaps): maximal '-' runs per row.
    [TestCase("ACGT", "ACGT", 0)]
    [TestCase("AC--GT", "ACGGGT", 1)]
    [TestCase("--CGT", "ACGGT", 1)]                // leading gap run counts (Biopython open_left_gaps)
    [TestCase("ACG--", "ACGGT", 1)]                // trailing gap run counts
    [TestCase("A-C-G", "ATCTG", 2)]                // two separate runs
    [TestCase("A-CG", "AT-G", 2)]                  // a run in each row: two openings
    [TestCase("AC-GT", "ACT-T", 2)]                // insertion run directly followed by deletion run: two openings
    [TestCase("AC--GT", "AC--GT", 2)]              // gap-gap columns: one run per row (not a valid pairwise alignment)
    public void CalculateStatistics_GapOpens_CountsGapRuns(string a, string b, int expected)
    {
        Assert.That(SequenceAligner.CalculateStatistics(Aln(a, b)).GapOpens, Is.EqualTo(expected));
    }

    // F13b — a pairwise alignment computed by the library: GapOpens equals the number of gap runs it contains,
    // and the affine aligner never opens more gaps than the linear one scores (sanity) — here the BlastDna
    // example from ALIGN-GLOBAL-001 has exactly one 3-nt gap.
    [Test]
    public void CalculateStatistics_OnAffineGlobalAlignment_SingleLongGap_OneOpening()
    {
        var aln = SequenceAligner.GlobalAlignAffine("ACGTTTTACGT", "ACGTACGT", SequenceAligner.BlastDna);
        var st = SequenceAligner.CalculateStatistics(aln, SequenceAligner.BlastDna);
        Assert.Multiple(() =>
        {
            Assert.That(st.Gaps, Is.EqualTo(3));
            Assert.That(st.GapOpens, Is.EqualTo(1));
            Assert.That(st.Matches, Is.EqualTo(8));
            Assert.That(st.AlignmentLength, Is.EqualTo(11));
        });
    }

    // F14 — a column's similarity is its score > 0, identical columns included (EMBOSS alignConsStats,
    // Biopython positives): X/X in BLOSUM62 is identical but scores -1, so it is NOT similar; a scalar model with
    // Match <= 0 likewise has no similar identical columns.
    [Test]
    public void CalculateStatistics_IdenticalPairWithNonPositiveScore_IsNotSimilar()
    {
        var st = SequenceAligner.CalculateStatistics(Aln("AXK", "AXR"), SubstitutionMatrix.Blosum62);
        Assert.Multiple(() =>
        {
            Assert.That(st.Matches, Is.EqualTo(2), "A/A and X/X are identical");
            Assert.That(st.Positives, Is.EqualTo(2), "A/A (+4) and K/R (+2) are positive; X/X (-1) is not");
        });

        var zeroMatch = new ScoringMatrix(Match: 0, Mismatch: -1, GapOpen: 0, GapExtend: -1);
        var st2 = SequenceAligner.CalculateStatistics(Aln("ACGT", "ACGA"), zeroMatch);
        Assert.Multiple(() =>
        {
            Assert.That(st2.Matches, Is.EqualTo(3));
            Assert.That(st2.Positives, Is.EqualTo(0));
            Assert.That(st2.Similarity, Is.EqualTo(0.0));
        });
    }

    // F14b — markup of an identical pair with a non-positive own score is blank (needle 6.6.0: X/X in BLOSUM62).
    [Test]
    public void FormatAlignment_IdenticalPairWithNonPositiveScore_IsBlank()
    {
        string markup = SequenceAligner.FormatAlignment(Aln("AXK", "AXR"), 60, SubstitutionMatrix.Blosum62).Split('\n')[1];
        Assert.That(markup, Is.EqualTo("| :"), "A/A '|', X/X (-1) blank, K/R (+2) ':'");
    }

    #endregion

    #region SubstitutionMatrix presets vs published tables

    // Spot values from the published tables (Henikoff & Henikoff 1992 BLOSUM62; Dayhoff PAM250; NCBI NUC.4.4 /
    // EMBOSS EDNAFULL) as printed by Biopython 1.88 substitution_matrices and EMBOSS 6.6 data files; the full
    // 24x24 / 24x24 / 16x16 tables are exercised through the 301-alignment oracle above.
    [TestCase('A', 'A', 4)]
    [TestCase('W', 'W', 11)]
    [TestCase('C', 'C', 9)]
    [TestCase('I', 'L', 2)]
    [TestCase('A', 'S', 1)]
    [TestCase('A', 'T', 0)]
    [TestCase('X', 'X', -1)]
    [TestCase('W', 'D', -4)]
    [TestCase('*', '*', 1)]
    [TestCase('a', 'a', 4)]
    public void Blosum62_PublishedValues(char a, char b, int expected)
    {
        Assert.That(SubstitutionMatrix.Blosum62.GetScore(a, b), Is.EqualTo(expected));
        Assert.That(SubstitutionMatrix.Blosum62.GetScore(b, a), Is.EqualTo(expected), "symmetric");
    }

    [TestCase('W', 'W', 17)]
    [TestCase('C', 'C', 12)]
    [TestCase('A', 'A', 2)]
    [TestCase('W', 'C', -8)]
    [TestCase('Y', 'F', 7)]
    public void Pam250_PublishedValues(char a, char b, int expected) =>
        Assert.That(SubstitutionMatrix.Pam250.GetScore(a, b), Is.EqualTo(expected));

    [TestCase('A', 'A', 5)]
    [TestCase('A', 'C', -4)]
    [TestCase('A', 'R', 1)]    // R = A/G
    [TestCase('A', 'N', -2)]
    [TestCase('N', 'N', -1)]
    [TestCase('S', 'S', -1)]
    [TestCase('T', 'U', 5)]
    [TestCase('A', 'V', -1)]
    public void EdnaFull_PublishedValues(char a, char b, int expected) =>
        Assert.That(SubstitutionMatrix.EdnaFull.GetScore(a, b), Is.EqualTo(expected));

    [Test]
    public void SubstitutionMatrix_Tables_AreSymmetric_AndRejectUnknownSymbols()
    {
        foreach (var m in new[] { SubstitutionMatrix.Blosum62, SubstitutionMatrix.Pam250, SubstitutionMatrix.EdnaFull })
        {
            foreach (char a in m.Alphabet)
                foreach (char b in m.Alphabet)
                    Assert.That(m.GetScore(a, b), Is.EqualTo(m.GetScore(b, a)), $"{m.Name} {a}{b}");
            Assert.That(m.TryGetScore('J', 'A', out _), Is.False);
            Assert.Throws<ArgumentException>(() => m.GetScore('J', 'A'));
        }

        Assert.That(SubstitutionMatrix.Blosum62.Alphabet, Is.EqualTo("ARNDCQEGHILKMFPSTWYVBZX*"));
        Assert.That(SubstitutionMatrix.EdnaFull.Alphabet, Is.EqualTo("ATGCSWRYKMBVHDNU"));
        Assert.Throws<ArgumentException>(() => new SubstitutionMatrix("x", "AA", new int[2, 2]));
        Assert.Throws<ArgumentException>(() => new SubstitutionMatrix("x", "AB", new int[3, 3]));
        Assert.Throws<ArgumentNullException>(() => SequenceAligner.CalculateStatistics(Aln("A", "A"), (SubstitutionMatrix)null!));
    }

    // Scalar presets vs published schemes: BlastDna = blastn defaults (reward 2, penalty -3, gap existence 5,
    // extension 2; BLAST+ manual); HighIdentityDna match/mismatch = EDNAFULL / NUC.4.4 on A,C,G,T (+5/-4);
    // SimpleDna = textbook +1/-1 (Needleman-Wunsch examples). Gap parameters of HighIdentityDna/SimpleDna are
    // library choices, not those of a named tool.
    [Test]
    public void ScalarPresets_MatchPublishedMatchMismatch()
    {
        var edna = SubstitutionMatrix.EdnaFull;
        foreach (char a in "ACGT")
            foreach (char b in "ACGT")
            {
                int expected = a == b ? SequenceAligner.HighIdentityDna.Match : SequenceAligner.HighIdentityDna.Mismatch;
                Assert.That(edna.GetScore(a, b), Is.EqualTo(expected), $"{a}{b}");
            }

        Assert.That(SequenceAligner.BlastDna, Is.EqualTo(new ScoringMatrix(2, -3, -5, -2)));
    }

    #endregion
}
