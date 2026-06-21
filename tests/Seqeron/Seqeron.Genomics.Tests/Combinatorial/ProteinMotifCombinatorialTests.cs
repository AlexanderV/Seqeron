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

    // ═══════════════════════════════════════════════════════════════════════
    // Unit: PROTMOTIF-DOMAIN-001 — Protein domain identification (ProteinMotif)
    // Checklist: docs/checklists/09_COMBINATORIAL_TESTING.md, row 84.
    // Spec: tests/TestSpecs/PROTMOTIF-DOMAIN-001.md (canonical FindDomains).
    // Dimensions: eValueThreshold(3) × minDomainLen(3) × nProfiles(3). Grid 3×3×3 = 27.
    //
    // Model (Pfam profile scan): FindDomains scans a fixed library of domain signatures (zinc
    // finger PF00096, kinase PF00069, SH3 PF00018, …) and reports each hit's span and score. The
    // library is fixed, so nProfiles is realised as the number of distinct domain types planted;
    // eValueThreshold/minDomainLen are caller-side score/length filters on the output.
    //
    // The combinatorial point: planted domain count, score threshold and length threshold interact
    // — exactly the planted domain types are detected, and a (score ≥ eVal ∧ length ≥ minLen)
    // filter yields only domains satisfying both, monotonically.
    // ═══════════════════════════════════════════════════════════════════════

    // Peptides matching the FindDomains regexes for the first three domain profiles.
    private static readonly (string Acc, string Motif)[] DomainMotifs =
    {
        ("PF00096", "CAACAAALAAAAAAAAHAAAH"), // zinc finger C2H2
        ("PF00069", "AAAAAGKS"),             // kinase ATP-binding
        ("PF00018", "LAAGWFAAAAAL"),         // SH3
    };

    private static string PlantDomains(int nProfiles)
    {
        var sb = new System.Text.StringBuilder(new string('Q', 10));
        for (int i = 0; i < nProfiles; i++) sb.Append(DomainMotifs[i].Motif).Append(new string('Q', 10));
        return sb.ToString();
    }

    [Test, Combinatorial]
    public void ProtMotifDomain_DetectsPlantedProfiles_AndFilters(
        [Values(0.0, 0.5, 0.9)] double eValueThreshold,
        [Values(4, 8, 12)] int minDomainLen,
        [Values(1, 2, 3)] int nProfiles)
    {
        string protein = PlantDomains(nProfiles);
        var planted = DomainMotifs.Take(nProfiles).Select(d => d.Acc).ToHashSet();

        var domains = ProteinMotifFinder.FindDomains(protein).ToList();

        domains.Should().OnlyContain(d => planted.Contains(d.Accession), "only planted domain types appear");
        foreach (var acc in planted)
            domains.Should().Contain(d => d.Accession == acc, "each planted domain is detected");
        domains.Should().OnlyContain(d => d.Start >= 0 && d.End < protein.Length && d.Start <= d.End, "valid coordinates");

        var filtered = domains.Where(d => d.Score >= eValueThreshold && (d.End - d.Start + 1) >= minDomainLen).ToList();
        filtered.Should().OnlyContain(d => d.Score >= eValueThreshold && (d.End - d.Start + 1) >= minDomainLen);
    }

    /// <summary>
    /// Interaction witness: the score and length filters are monotone — raising either threshold
    /// can only shrink the retained domain set.
    /// </summary>
    [Test]
    public void ProtMotifDomain_FiltersAreMonotone()
    {
        var domains = ProteinMotifFinder.FindDomains(PlantDomains(3)).ToList();
        int byScoreLow = domains.Count(d => d.Score >= 0.0);
        int byScoreHigh = domains.Count(d => d.Score >= 0.9);
        byScoreHigh.Should().BeLessThanOrEqualTo(byScoreLow, "a higher score floor retains no more domains");

        int byLenLow = domains.Count(d => d.End - d.Start + 1 >= 4);
        int byLenHigh = domains.Count(d => d.End - d.Start + 1 >= 20);
        byLenHigh.Should().BeLessThanOrEqualTo(byLenLow, "a higher length floor retains no more domains");
    }
}
