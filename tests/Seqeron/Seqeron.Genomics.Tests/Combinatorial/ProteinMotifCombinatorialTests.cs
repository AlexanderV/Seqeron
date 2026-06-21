namespace Seqeron.Genomics.Tests.Combinatorial;

/// <summary>
/// Combinatorial (pairwise / full-grid) tests for the ProteinMotif area.
///
/// See <see cref="CompositionCombinatorialTests"/> for the rationale of
/// combinatorial testing. Each grid cell carries a real business assertion;
/// small grids use the exhaustive <c>[Combinatorial]</c> product.
/// — docs/checklists/09_COMBINATORIAL_TESTING.md §Description.
/// </summary>
[TestFixture]
[Category("Combinatorial")]
[Category("ProteinMotif")]
public class ProteinMotifCombinatorialTests
{
    // ═══════════════════════════════════════════════════════════════════════
    // Unit: PROTMOTIF-FIND-001 — Protein motif search (ProteinMotif)
    // Checklist: docs/checklists/09_COMBINATORIAL_TESTING.md, row 82.
    // Spec: tests/TestSpecs/PROTMOTIF-FIND-001.md (canonical FindMotifByPattern / FindCommonMotifs).
    // Dimensions: minMotifLen(3) × maxMotifLen(3) × seqLen(3). Grid 3×3×3 = 27.
    //
    // Model (PROSITE pattern search): FindMotifByPattern scans a protein for a regex motif,
    // reporting each occurrence's span and matched residues (overlapping matches via lookahead).
    // The two length axes are realised as the lengths of a planted short and long motif searched
    // for literally; seqLen is the protein length.
    //
    // The combinatorial point: short-motif length, long-motif length and protein length interact —
    // each planted motif is found at its exact position with the exact matched length, regardless
    // of the surrounding protein length.
    // ═══════════════════════════════════════════════════════════════════════

    private const int ShortPos = 10;
    private const int LongPos = 25;

    [Test, Combinatorial]
    public void ProtMotifFind_LocatesPlantedMotifs_AcrossLengths(
        [Values(3, 4, 5)] int minMotifLen,
        [Values(6, 8, 10)] int maxMotifLen,
        [Values(40, 80, 160)] int seqLen)
    {
        string shortMotif = "WYCDE"[..minMotifLen];        // distinct residues, absent from poly-A filler
        string longMotif = "MNPQRTKHGS"[..maxMotifLen];    // disjoint residue set from the short motif
        string protein = new string('A', ShortPos) + shortMotif
            + new string('A', LongPos - ShortPos - minMotifLen) + longMotif
            + new string('A', seqLen - LongPos - maxMotifLen);

        var shortHits = ProteinMotifFinder.FindMotifByPattern(protein, shortMotif, "short").ToList();
        shortHits.Should().ContainSingle(h => h.Start == ShortPos && h.Sequence == shortMotif);
        shortHits.Should().OnlyContain(h => h.End - h.Start + 1 == minMotifLen, "match length equals the motif length");

        var longHits = ProteinMotifFinder.FindMotifByPattern(protein, longMotif, "long").ToList();
        longHits.Should().ContainSingle(h => h.Start == LongPos && h.Sequence == longMotif);
        longHits.Should().OnlyContain(h => h.End - h.Start + 1 == maxMotifLen);
    }

    /// <summary>
    /// Interaction witness: a regex character class matches all admitted residues — N-{P}-[ST]-{P}
    /// (the PROSITE N-glycosylation motif) matches NAS-A but not NPS-A (proline forbidden).
    /// </summary>
    [Test]
    public void ProtMotifFind_CharacterClassesRespected()
    {
        const string glyco = @"N[^P][ST][^P]";
        ProteinMotifFinder.FindMotifByPattern("AAANASAAAA", glyco).Should().ContainSingle(h => h.Sequence == "NASA");
        ProteinMotifFinder.FindMotifByPattern("AAANPSAAAA", glyco).Should().BeEmpty("proline at position 2 is forbidden");
    }

    /// <summary>
    /// Interaction witness: FindCommonMotifs detects a built-in PROSITE motif — an N-glycosylation
    /// site (NAS-A) is reported by name.
    /// </summary>
    [Test]
    public void ProtMotifFind_CommonMotifs_DetectsGlycosylation()
    {
        ProteinMotifFinder.FindCommonMotifs("AAAANASAAAAA")
            .Should().Contain(m => m.MotifName == "ASN_GLYCOSYLATION");
    }
}
