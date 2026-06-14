// MOTIF-GENERATE-001 — IUPAC-Degenerate Consensus Generation
// Evidence: docs/Evidence/MOTIF-GENERATE-001-Evidence.md
// TestSpec: tests/TestSpecs/MOTIF-GENERATE-001.md
// Source: Cornish-Bowden A. (1985). Nomenclature for incompletely specified bases in nucleic
//         acid sequences: recommendations 1984. Nucleic Acids Research 13(9):3021. DOI 10.1093/nar/13.9.3021.
//         UCSC IUPAC ambiguity codes; Wikipedia "Nucleic acid notation" (NC-IUB 1984 table);
//         DECIPHER ConsensusSequence (threshold-consensus mechanism).

using System;
using NUnit.Framework;
using Seqeron.Genomics.Analysis;

namespace Seqeron.Genomics.Tests;

/// <summary>
/// Canonical test class for MOTIF-GENERATE-001: IUPAC-degenerate consensus generation. Verifies
/// <see cref="MotifFinder.GenerateConsensus(System.Collections.Generic.IEnumerable{string})"/>
/// against the NC-IUB 1984 set→symbol mapping and the documented 25% inclusion threshold.
/// Expected symbols are derived from the authoritative IUPAC table, not from the implementation.
/// </summary>
[TestFixture]
public class MotifFinder_GenerateConsensus_Tests
{
    #region GenerateConsensus — MUST (NC-IUB set→symbol mapping)

    // M1 — {A,G} both above threshold → R (purine). Source: NC-IUB 1984 / UCSC.
    // n=2, threshold=0.5; A=1>0.5, G=1>0.5 → set {A,G} → R.
    [Test]
    public void GenerateConsensus_ColumnAG_ReturnsR()
    {
        var consensus = MotifFinder.GenerateConsensus(new[] { "A", "G" });
        Assert.That(consensus, Is.EqualTo("R"),
            "NC-IUB: the set {A,G} (purine) is encoded by the IUPAC symbol R.");
    }

    // M2 — {C,T} → Y (pyrimidine). Source: NC-IUB 1984 / UCSC.
    [Test]
    public void GenerateConsensus_ColumnCT_ReturnsY()
    {
        var consensus = MotifFinder.GenerateConsensus(new[] { "C", "T" });
        Assert.That(consensus, Is.EqualTo("Y"),
            "NC-IUB: the set {C,T} (pyrimidine) is encoded by the IUPAC symbol Y.");
    }

    // M3 — {C,G} → S (strong). Source: NC-IUB 1984 / UCSC.
    [Test]
    public void GenerateConsensus_ColumnCG_ReturnsS()
    {
        var consensus = MotifFinder.GenerateConsensus(new[] { "C", "G" });
        Assert.That(consensus, Is.EqualTo("S"),
            "NC-IUB: the set {C,G} (strong) is encoded by the IUPAC symbol S.");
    }

    // M4 — {A,T} → W (weak). Source: NC-IUB 1984 / UCSC.
    [Test]
    public void GenerateConsensus_ColumnAT_ReturnsW()
    {
        var consensus = MotifFinder.GenerateConsensus(new[] { "A", "T" });
        Assert.That(consensus, Is.EqualTo("W"),
            "NC-IUB: the set {A,T} (weak) is encoded by the IUPAC symbol W.");
    }

    // M5 — {G,T} → K (keto). Source: NC-IUB 1984 / UCSC.
    [Test]
    public void GenerateConsensus_ColumnGT_ReturnsK()
    {
        var consensus = MotifFinder.GenerateConsensus(new[] { "G", "T" });
        Assert.That(consensus, Is.EqualTo("K"),
            "NC-IUB: the set {G,T} (keto) is encoded by the IUPAC symbol K.");
    }

    // M6 — {A,C} → M (amino). Source: NC-IUB 1984 / UCSC.
    [Test]
    public void GenerateConsensus_ColumnAC_ReturnsM()
    {
        var consensus = MotifFinder.GenerateConsensus(new[] { "A", "C" });
        Assert.That(consensus, Is.EqualTo("M"),
            "NC-IUB: the set {A,C} (amino) is encoded by the IUPAC symbol M.");
    }

    // M7 — {C,G,T} → B (not-A). n=3, threshold=0.75; each count=1>0.75. Source: NC-IUB 1984 / UCSC.
    [Test]
    public void GenerateConsensus_ColumnCGT_ReturnsB()
    {
        var consensus = MotifFinder.GenerateConsensus(new[] { "C", "G", "T" });
        Assert.That(consensus, Is.EqualTo("B"),
            "NC-IUB: the set {C,G,T} (not-A) is encoded by the IUPAC symbol B.");
    }

    // M8 — {A,G,T} → D (not-C). Source: NC-IUB 1984 / UCSC.
    [Test]
    public void GenerateConsensus_ColumnAGT_ReturnsD()
    {
        var consensus = MotifFinder.GenerateConsensus(new[] { "A", "G", "T" });
        Assert.That(consensus, Is.EqualTo("D"),
            "NC-IUB: the set {A,G,T} (not-C) is encoded by the IUPAC symbol D.");
    }

    // M9 — {A,C,T} → H (not-G). Source: NC-IUB 1984 / UCSC.
    [Test]
    public void GenerateConsensus_ColumnACT_ReturnsH()
    {
        var consensus = MotifFinder.GenerateConsensus(new[] { "A", "C", "T" });
        Assert.That(consensus, Is.EqualTo("H"),
            "NC-IUB: the set {A,C,T} (not-G) is encoded by the IUPAC symbol H.");
    }

    // M10 — {A,C,G} → V (not-T). Source: NC-IUB 1984 / UCSC.
    [Test]
    public void GenerateConsensus_ColumnACG_ReturnsV()
    {
        var consensus = MotifFinder.GenerateConsensus(new[] { "A", "C", "G" });
        Assert.That(consensus, Is.EqualTo("V"),
            "NC-IUB: the set {A,C,G} (not-T) is encoded by the IUPAC symbol V.");
    }

    // M11 — Unanimous columns reproduce the input (INV-02): singleton set → standard base.
    [Test]
    public void GenerateConsensus_IdenticalSequences_ReturnsThatSequence()
    {
        var consensus = MotifFinder.GenerateConsensus(new[] { "ATGC", "ATGC", "ATGC" });
        Assert.That(consensus, Is.EqualTo("ATGC"),
            "Every column is unanimous; a singleton base set maps to its standard base (NC-IUB).");
    }

    // M12 — Multi-column mixed: col0 {A,G}→R, cols 1-3 unanimous → "RTGC". Source: NC-IUB 1984.
    [Test]
    public void GenerateConsensus_MixedFirstColumn_ReturnsRtgc()
    {
        var consensus = MotifFinder.GenerateConsensus(new[] { "ATGC", "GTGC" });
        Assert.That(consensus, Is.EqualTo("RTGC"),
            "Column 0 holds {A,G}→R; remaining columns are unanimous (T,G,C).");
    }

    // M13 — Strict 25% boundary (INV-05): a base at exactly the threshold is excluded.
    // n=4, threshold=1.0. col0/col1: A=4 → 'A'. col2: A=C=G=T=1, none >1.0 → fallback most-frequent
    // (alphabetical tie → 'A'). col3: A=1,G=1,T=2 → only T(2)>1.0 → singleton {T} → 'T'.
    [Test]
    public void GenerateConsensus_ExactlyQuarterBoundary_ExcludesBaseAtThreshold()
    {
        var consensus = MotifFinder.GenerateConsensus(new[] { "AAAA", "AAGT", "AACT", "AATT" });
        Assert.Multiple(() =>
        {
            Assert.That(consensus, Has.Length.EqualTo(4),
                "One IUPAC symbol per column (INV-01).");
            Assert.That(consensus[0], Is.EqualTo('A'), "Column 0 is all A.");
            Assert.That(consensus[1], Is.EqualTo('A'), "Column 1 is all A.");
            Assert.That(consensus[2], Is.EqualTo('A'),
                "Column 2 has four bases each at exactly 25% (count 1 = threshold 1.0); strict '>' " +
                "excludes all, so the fallback picks the most-frequent base (alphabetical tie → A).");
            Assert.That(consensus[3], Is.EqualTo('T'),
                "Column 3: only T (count 2 > threshold 1.0) passes; A and G at exactly 25% are dropped.");
        });
    }

    // M14 — Minority base below threshold is dropped before IUPAC encoding (DECIPHER threshold rule).
    // n=5, threshold=1.25. A=2>1.25, G=2>1.25, C=1≤1.25 dropped → set {A,G} → R (not a 3-base code).
    [Test]
    public void GenerateConsensus_MinorityBelowThreshold_DroppedYieldsR()
    {
        var consensus = MotifFinder.GenerateConsensus(new[] { "A", "A", "G", "G", "C" });
        Assert.That(consensus, Is.EqualTo("R"),
            "C (count 1 ≤ threshold 1.25) is below 25% and dropped; surviving {A,G} → R, not B/V.");
    }

    // M15 — No base passes the threshold → fallback to single most-frequent base (alphabetical tie).
    // n=4, threshold=1.0; each base count=1, none >1.0 → fallback → 'A' (not 'N').
    [Test]
    public void GenerateConsensus_FourEqualBases_FallbackToA_NotN()
    {
        var consensus = MotifFinder.GenerateConsensus(new[] { "A", "C", "G", "T" });
        Assert.That(consensus, Is.EqualTo("A"),
            "Four bases each at exactly 25%; strict '>' lets none pass, so the most-frequent " +
            "fallback (alphabetical tie) returns 'A', never the four-base symbol N.");
    }

    #endregion

    #region GenerateConsensus — SHOULD (edge / invariants)

    // S1 — Case-insensitive: lowercase input gives the same result as upper-case (INV / normalisation).
    [Test]
    public void GenerateConsensus_LowercaseInput_SameAsUppercase()
    {
        var lower = MotifFinder.GenerateConsensus(new[] { "atgc", "gtgc" });
        Assert.That(lower, Is.EqualTo("RTGC"),
            "Input is upper-cased before counting; lowercase yields the same RTGC consensus.");
    }

    // S2 — Empty collection → "" (INV-06).
    [Test]
    public void GenerateConsensus_EmptyCollection_ReturnsEmpty()
    {
        var consensus = MotifFinder.GenerateConsensus(Array.Empty<string>());
        Assert.That(consensus, Is.Empty,
            "An empty collection has no columns to summarise → empty string.");
    }

    // S3 — Output length equals the first sequence's length (INV-01).
    [Test]
    public void GenerateConsensus_OutputLength_MatchesFirstSequence()
    {
        var consensus = MotifFinder.GenerateConsensus(new[] { "ACGTACG", "ACGTACG" });
        Assert.That(consensus, Has.Length.EqualTo(7),
            "One symbol per column; the column count equals the first sequence's length.");
    }

    #endregion

    #region GenerateConsensus — COULD (guards)

    // C1 — Null collection throws ArgumentNullException (documented guard).
    [Test]
    public void GenerateConsensus_Null_ThrowsArgumentNullException()
    {
        Assert.That(() => MotifFinder.GenerateConsensus(null!),
            NUnit.Framework.Throws.TypeOf<ArgumentNullException>(),
            "A null sequence collection is rejected with ArgumentNullException.");
    }

    #endregion
}
