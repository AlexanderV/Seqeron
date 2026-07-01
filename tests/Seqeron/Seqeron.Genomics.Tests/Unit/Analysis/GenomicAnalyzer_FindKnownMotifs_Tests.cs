// GENOMIC-MOTIFS-001 — Known Motif Search (multi-pattern exact substring matching)
// Evidence: docs/Evidence/GENOMIC-MOTIFS-001-Evidence.md
// TestSpec: tests/TestSpecs/GENOMIC-MOTIFS-001.md
// Source: Gusfield D. (1997). Algorithms on Strings, Trees and Sequences. Cambridge Univ. Press, ISBN 0-521-58519-8
//         (exact-matching definition, Tufts COMP 150GEN exact.html); Biopython Bio.Seq (search/count_overlap);
//         Wikipedia "Restriction site" (EcoRI=GAATTC) — all accessed 2026-06-13.

using System;
using System.Collections.Generic;
using System.Linq;
using NUnit.Framework;

namespace Seqeron.Genomics.Tests.Unit.Analysis;

[TestFixture]
public class GenomicAnalyzer_FindKnownMotifs_Tests
{
    #region FindKnownMotifs

    // M1 — Overlapping occurrences must all be reported. Gusfield/Tufts exact.html: P=aaa, T=aaaaa
    // has "three (overlapping) occurrences" -> 0-based starts {0,1,2}.
    [Test]
    public void FindKnownMotifs_OverlappingMotif_ReportsAllStarts()
    {
        var seq = new DnaSequence("AAAAA");

        var result = GenomicAnalyzer.FindKnownMotifs(seq, new[] { "AAA" });

        Assert.Multiple(() =>
        {
            Assert.That(result.ContainsKey("AAA"), Is.True, "AAA occurs in AAAAA so it must be a result key.");
            Assert.That(result["AAA"], Is.EqualTo(new[] { 0, 1, 2 }),
                "Exact matching reports all overlapping occurrences: AAA aligns at starts 0,1,2 in AAAAA (Gusfield/Tufts).");
        });
    }

    // M2 — EcoRI biological motif GAATTC (Wikipedia "Restriction site"). T=GAATTCAAAGAATTC: GAATTC starts
    // at index 0 (GAATTC) and index 9 (chars 9..14 = GAATTC), derived from the exact-matching definition.
    [Test]
    public void FindKnownMotifs_EcoRIMotif_FindsAllSites()
    {
        var seq = new DnaSequence("GAATTCAAAGAATTC");

        var result = GenomicAnalyzer.FindKnownMotifs(seq, new[] { "GAATTC" });

        Assert.Multiple(() =>
        {
            Assert.That(result["GAATTC"], Is.EqualTo(new[] { 0, 9 }),
                "EcoRI site GAATTC occurs at 0-based starts 0 and 9 in GAATTCAAAGAATTC (exact-matching set).");
            Assert.That(result.Count, Is.EqualTo(1), "Only the single queried motif appears in the result.");
        });
    }

    // M3 — Multi-motif set (Biopython Seq.search: per-motif hits). T=ACGTACGTAA:
    // ACGT at {0,4}; AA at {8} (chars 8..9); TTT absent and therefore omitted.
    [Test]
    public void FindKnownMotifs_MotifSet_ReturnsPerMotifPositions()
    {
        var seq = new DnaSequence("ACGTACGTAA");

        var result = GenomicAnalyzer.FindKnownMotifs(seq, new[] { "ACGT", "AA", "TTT" });

        Assert.Multiple(() =>
        {
            Assert.That(result["ACGT"], Is.EqualTo(new[] { 0, 4 }), "ACGT occurs at starts 0 and 4 in ACGTACGTAA.");
            Assert.That(result["AA"], Is.EqualTo(new[] { 8 }), "AA occurs only at start 8 (chars 8..9) in ACGTACGTAA.");
            Assert.That(result.ContainsKey("TTT"), Is.False, "TTT does not occur, so it is omitted from the result (INV-04).");
            Assert.That(result.Count, Is.EqualTo(2), "Exactly two of the three queried motifs occur.");
        });
    }

    // M4 — Positions are returned sorted ascending (deterministic contract, INV-03). The suffix tree
    // enumerates in DFS order; the implementation must sort. ACGT in ACGTACGTAA -> [0,4] ascending.
    [Test]
    public void FindKnownMotifs_Positions_AreSortedAscending()
    {
        var seq = new DnaSequence("ACGTACGTAA");

        var result = GenomicAnalyzer.FindKnownMotifs(seq, new[] { "ACGT" });
        var positions = result["ACGT"];

        Assert.That(positions, Is.Ordered.Ascending,
            "Each motif's positions must be sorted ascending for a stable, deterministic contract (INV-03).");
        Assert.That(positions, Is.EqualTo(new[] { 0, 4 }), "Sorted occurrence set of ACGT is [0,4].");
    }

    // M5 — Absent motif yields an empty result dictionary (empty occurrence set, INV-04).
    [Test]
    public void FindKnownMotifs_AbsentMotif_OmittedFromResult()
    {
        var seq = new DnaSequence("ACGT");

        var result = GenomicAnalyzer.FindKnownMotifs(seq, new[] { "TTTT" });

        Assert.That(result, Is.Empty, "A motif that does not occur produces no entry, leaving an empty result.");
    }

    // S1 — Empty motif set -> empty dictionary (degenerate input).
    [Test]
    public void FindKnownMotifs_EmptyMotifSet_ReturnsEmpty()
    {
        var seq = new DnaSequence("ACGTACGT");

        var result = GenomicAnalyzer.FindKnownMotifs(seq, Array.Empty<string>());

        Assert.That(result, Is.Empty, "No motifs to search -> empty result dictionary.");
    }

    // S2 — Lower-case motif is normalized and matched; result key is upper-cased (Biopython uppercase, INV-05).
    [Test]
    public void FindKnownMotifs_LowercaseMotif_NormalizedAndMatched()
    {
        var seq = new DnaSequence("GAATTC");

        var result = GenomicAnalyzer.FindKnownMotifs(seq, new[] { "gaattc" });

        Assert.Multiple(() =>
        {
            Assert.That(result.ContainsKey("GAATTC"), Is.True, "Motif is upper-cased; result is keyed by GAATTC (INV-05).");
            Assert.That(result.ContainsKey("gaattc"), Is.False, "The lower-case form is not a key after normalization.");
            Assert.That(result["GAATTC"], Is.EqualTo(new[] { 0 }), "GAATTC occurs once, at start 0.");
        });
    }

    // S3 — Empty/whitespace motif is skipped (empty string is not a motif; mirrors FindMotif). ASSUMPTION (Evidence §1).
    [Test]
    public void FindKnownMotifs_EmptyMotif_Skipped()
    {
        var seq = new DnaSequence("ACGT");

        var result = GenomicAnalyzer.FindKnownMotifs(seq, new[] { "", "   " });

        Assert.That(result, Is.Empty,
            "Empty/whitespace motifs are not motifs and are skipped (a suffix-tree search for \"\" would match every position).");
    }

    // S4 — Single-character motif counts every position, including consecutive (overlapping) ones.
    // T=AAAAA, motif A -> {0,1,2,3,4} (overlap-rule boundary, Gusfield/Tufts).
    [Test]
    public void FindKnownMotifs_SingleCharMotif_ReportsEveryPosition()
    {
        var seq = new DnaSequence("AAAAA");

        var result = GenomicAnalyzer.FindKnownMotifs(seq, new[] { "A" });

        Assert.That(result["A"], Is.EqualTo(new[] { 0, 1, 2, 3, 4 }),
            "A occurs at every position of AAAAA; all are reported (INV-02).");
    }

    // C1 — Duplicate motifs that normalize to the same upper-cased key collapse to one entry.
    [Test]
    public void FindKnownMotifs_DuplicateMotifKeys_Deduplicated()
    {
        var seq = new DnaSequence("ACAA");

        var result = GenomicAnalyzer.FindKnownMotifs(seq, new[] { "AA", "aa" });

        Assert.Multiple(() =>
        {
            Assert.That(result.Count, Is.EqualTo(1), "AA and aa normalize to the same key -> a single entry.");
            Assert.That(result["AA"], Is.EqualTo(new[] { 2 }), "AA occurs only at start 2 (chars 2..3) in ACAA.");
        });
    }

    // Edge — null motifs throws ArgumentNullException (input validation).
    [Test]
    public void FindKnownMotifs_NullMotifs_Throws()
    {
        var seq = new DnaSequence("ACGT");

        Assert.That(() => GenomicAnalyzer.FindKnownMotifs(seq, null!),
            NUnit.Framework.Throws.ArgumentNullException, "Null motif collection is invalid input.");
    }

    #endregion
}
