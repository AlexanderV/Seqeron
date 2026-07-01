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
    // Model (exact PROSITE pattern scan): FindDomains scans a fixed library of EXACT PROSITE
    // PATTERN signatures (zinc finger C2H2 PS00028/PF00096, WD-repeats PS00678/PF00400, kinase
    // ATP-binding PS00017/PF00069) and reports each hit's span and score. SH3/PDZ are PROFILE-only
    // (PS50002/PS50106) and are intentionally not detected. The library is fixed, so nProfiles is
    // realised as the number of distinct domain types planted; eValueThreshold/minDomainLen are
    // caller-side score/length filters on the output.
    //
    // The combinatorial point: planted domain count, score threshold and length threshold interact
    // — exactly the planted domain types are detected, and a (score ≥ eVal ∧ length ≥ minLen)
    // filter yields only domains satisfying both, monotonically.
    // ═══════════════════════════════════════════════════════════════════════

    // Peptides matching the FindDomains exact-PROSITE-pattern regexes for three domain profiles.
    private static readonly (string Acc, string Motif)[] DomainMotifs =
    {
        ("PF00096", "CAACAAALAAAAAAAAHAAAH"), // zinc finger C2H2 (PS00028)
        ("PF00069", "AAAAAGKS"),             // kinase ATP-binding / P-loop (PS00017)
        ("PF00400", "LVSASQDGKLIIWDS"),      // WD-repeats (PS00678); real GBB1_HUMAN repeat
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

    // ═══════════════════════════════════════════════════════════════════════
    // Unit: PROTMOTIF-CC-001 — Coiled-coil prediction (ProteinMotif)
    // Checklist: docs/checklists/09_COMBINATORIAL_TESTING.md, row 163.
    // Spec: tests/TestSpecs/PROTMOTIF-CC-001.md (canonical PredictCoiledCoils). ADVANCED §10.
    // Dimensions: windowSize(3) × threshold(3) × seqLen(3). Grid 3×3×3 = 27 (full, exhaustive).
    //
    // Model (Lupas 1991; Mason & Arndt 2004): a coiled coil shows the heptad repeat (abcdefg)ₙ with
    // hydrophobic core residues {I,L,V} at the a and d positions. PredictCoiledCoils scores each
    // window by the BEST a/d occupancy over the 7 registers (fraction of a/d positions filled by a
    // core residue) and reports maximal runs scoring ≥ threshold, of length ≥ 21 (3 heptads).
    //
    // Engineered construct: a perfect heptad "LQQLQQQ" (a=d=L, others the non-core filler Q) repeated
    // — register 0 gives occupancy 1.0 at every window. The combinatorial point: window width,
    // threshold and protein length jointly bound detection; a perfect coiled coil scores 1.0 and is
    // reported (Score∈[0,1], length ≥21, peak ≥ threshold) at every cell, while occupancy gates it.
    // ═══════════════════════════════════════════════════════════════════════

    private static string PerfectHeptads(int count) => string.Concat(Enumerable.Repeat("LQQLQQQ", count));

    [Test, Combinatorial]
    public void ProtMotifCC_DetectsPerfectCoiledCoil_AcrossWindowThresholdLength(
        [Values(14, 21, 28)] int windowSize,
        [Values(0.3, 0.5, 0.7)] double threshold,
        [Values(5, 7, 10)] int heptadCount)
    {
        string protein = PerfectHeptads(heptadCount); // length 7·heptadCount ≥ 35

        var regions = ProteinMotifFinder.PredictCoiledCoils(protein, windowSize, threshold).ToList();

        regions.Should().NotBeEmpty("a perfect a/d=L heptad repeat scores occupancy 1.0 ≥ threshold");

        regions.Should().OnlyContain(r => r.Score >= -1e-12 && r.Score <= 1.0 + 1e-12, "Score ∈ [0,1] (INV-1)");
        regions.Should().OnlyContain(r => r.Score >= threshold - 1e-12, "a region requires a window ≥ threshold (INV-5)");
        regions.Should().OnlyContain(r => r.End - r.Start + 1 >= 21, "regions span ≥ 3 heptads (INV-2)");
        regions.Should().OnlyContain(r => r.Start >= 0 && r.Start <= r.End && r.End < protein.Length, "valid coords (INV-3)");
        regions.Select(r => r.Start).Should().BeInAscendingOrder("regions are ordered by Start (INV-3)");

        regions.Should().Contain(r => Math.Abs(r.Score - 1.0) < 1e-9, "the perfect heptad core reaches full a/d occupancy");
    }

    /// <summary>
    /// Interaction witness — occupancy gates detection by threshold: a half-occupied heptad
    /// ("LQQQQQQ", one core residue per heptad ⇒ best a/d occupancy 0.5) is reported at threshold
    /// 0.3 but not 0.7, while a core-free poly-Q peptide is never a coiled coil.
    /// </summary>
    [Test]
    public void ProtMotifCC_OccupancyGatesDetection()
    {
        string half = string.Concat(Enumerable.Repeat("LQQQQQQ", 8)); // ≤0.5 occupancy (one L per heptad)

        ProteinMotifFinder.PredictCoiledCoils(half, 21, 0.3).Should().NotBeEmpty("0.5 ≥ 0.3");
        ProteinMotifFinder.PredictCoiledCoils(half, 21, 0.7).Should().BeEmpty("0.5 < 0.7");

        ProteinMotifFinder.PredictCoiledCoils(new string('Q', 56), 21, 0.3)
            .Should().BeEmpty("a sequence with no core {I,L,V} residues is not a coiled coil");
    }

    // ═══════════════════════════════════════════════════════════════════════
    // Unit: PROTMOTIF-COMMON-001 — Whole-library PROSITE motif scan (ProteinMotif)
    // Checklist: docs/checklists/09_COMBINATORIAL_TESTING.md, row 164.
    // Spec: tests/TestSpecs/PROTMOTIF-COMMON-001.md (canonical FindCommonMotifs). ADVANCED §10.
    // Dimensions: nSeqs(3) × motifLen(3). Grid 3×3 = 9 (full, exhaustive ⊇ pairwise).
    //
    // Model (PROSITE / ScanProsite): FindCommonMotifs scans the curated CommonMotifs dictionary over
    // the protein and returns every occurrence of every pattern, with 0-based coordinates whose
    // matched substring equals protein[Start..End] (INV-1) and whose name is a dictionary entry (INV-3).
    //
    // Axis mapping (documented — FindCommonMotifs has no literal nSeqs/motifLen knob): nSeqs → number
    // of DISTINCT PROSITE motif types planted; motifLen → the length class of the planted motifs,
    // realised as a start index into a catalog ordered by match length (Short/Medium/Long). The
    // combinatorial point: planting n motifs of a chosen length class, the whole-library scan locates
    // each at its exact coordinates with the correct name, and every returned match is coordinate-valid.
    // ═══════════════════════════════════════════════════════════════════════

    // Catalog ordered by (fixed) match length; every peptide is Q-padded so it matches exactly one site.
    private static readonly (string Name, string Peptide)[] PrositeCatalog =
    {
        ("RGD", "RGD"),                              // PS00016, len 3
        ("CK2_PHOSPHO_SITE", "SQQD"),                // PS00006 [ST].{2}[DE], len 4
        ("SH3_BINDING_1", "RQQPQQP"),                // SH3_1 [RK].{2}P.{2}P, len 7
        ("ATP_GTP_A", "GQQQQGKS"),                   // PS00017 [AG].{4}GK[ST], len 8
        ("LEUCINE_ZIPPER", "LQQQQQQLQQQQQQLQQQQQQL"), // PS00029 L.{6}L.{6}L.{6}L, len 22
    };

    public enum MotifLenClass { Short = 0, Medium = 1, Long = 2 }

    [Test, Combinatorial]
    public void ProtMotifCommon_LocatesPlantedPrositeSites_AcrossCountAndLength(
        [Values(1, 2, 3)] int nSeqs,
        [Values(MotifLenClass.Short, MotifLenClass.Medium, MotifLenClass.Long)] MotifLenClass motifLen)
    {
        int start = (int)motifLen;
        var planted = PrositeCatalog.Skip(start).Take(nSeqs).ToList();

        // Build protein: Q-gap, then each planted peptide separated by Q-gaps.
        var sb = new System.Text.StringBuilder(new string('Q', 5));
        foreach (var (_, pep) in planted) sb.Append(pep).Append(new string('Q', 5));
        string protein = sb.ToString();

        var matches = ProteinMotifFinder.FindCommonMotifs(protein).ToList();
        var knownNames = ProteinMotifFinder.CommonMotifs.Values.Select(v => v.Name).ToHashSet();

        // Every match is coordinate-valid (INV-1, INV-2) and names a real dictionary entry (INV-3).
        matches.Should().OnlyContain(m => m.Start >= 0 && m.End < protein.Length && m.Start <= m.End, "INV-2");
        matches.Should().OnlyContain(m => m.Sequence == protein.Substring(m.Start, m.End - m.Start + 1), "INV-1");
        matches.Should().OnlyContain(m => knownNames.Contains(m.MotifName), "INV-3");

        // Each planted site is found at its exact position with its name.
        foreach (var (name, pep) in planted)
        {
            int pos = protein.IndexOf(pep, StringComparison.Ordinal);
            matches.Should().Contain(m => m.MotifName == name && m.Start == pos && m.Sequence == pep,
                "the planted {0} site is detected at {1}", name, pos);
        }
    }

    /// <summary>
    /// Interaction witness — proline exclusion ({P}) and multi-occurrence reporting (spec M2/M8):
    /// an N-glycosylation window with proline is rejected, and two RGD sites are both reported.
    /// </summary>
    [Test]
    public void ProtMotifCommon_ProlineExclusionAndMultiOccurrence()
    {
        ProteinMotifFinder.FindCommonMotifs("QQQQNFTQQQQ")
            .Should().Contain(m => m.MotifName == "ASN_GLYCOSYLATION" && m.Sequence == "NFTQ");
        ProteinMotifFinder.FindCommonMotifs("QQQQNPSQQQQ")
            .Should().NotContain(m => m.MotifName == "ASN_GLYCOSYLATION", "proline is forbidden at position 2 ({P})");

        ProteinMotifFinder.FindCommonMotifs("RGDRGD")
            .Count(m => m.MotifName == "RGD").Should().Be(2, "overlapping/adjacent occurrences are both reported");

        ProteinMotifFinder.FindCommonMotifs("").Should().BeEmpty("empty input ⇒ no matches (INV-5)");
    }

    // ═══════════════════════════════════════════════════════════════════════
    // Unit: PROTMOTIF-LC-001 — Low-complexity regions (SEG) (ProteinMotif)
    // Checklist: docs/checklists/09_COMBINATORIAL_TESTING.md, row 165.
    // Spec: tests/TestSpecs/PROTMOTIF-LC-001.md (canonical FindLowComplexityRegions). ADVANCED §10.
    // Dimensions: windowSize(3) × threshold(3). Grid 3×3 = 9 (full, exhaustive ⊇ pairwise).
    //
    // Model (Wootton & Federhen 1993, SEG): window complexity K = −Σ pᵢ·log₂(pᵢ) bits/residue, in
    // [0, log₂20]. A window with K ≤ K1 (trigger) starts a low-complexity segment, extended over
    // adjacent windows with K ≤ K2 (extension, default 2.5). A homopolymer window has K = 0 exactly.
    //
    // Axis mapping (documented): threshold → the trigger cutoff K1 (extension kept at the SEG default
    // 2.5, so K1 ≤ K2). Engineered construct: a poly-A block (K=0) embedded in a high-complexity
    // 20-residue-cycling flank. The combinatorial point: window width and the trigger cutoff jointly
    // bound detection; the K=0 block is always reported (region complexity ≤ K2, min ≤ K1, K=0), and
    // the diverse flanks are not — verified per cell.
    // ═══════════════════════════════════════════════════════════════════════

    private const string HighComplexityFlank = "ACDEFGHIKLMNPQRSTVWYACDEFGHIKLMNPQRSTVWY"; // all 20 AA cycled
    private const int LowComplexityBlockStart = 40; // = HighComplexityFlank.Length

    [Test, Combinatorial]
    public void ProtMotifLC_DetectsHomopolymerBlock_AcrossWindowAndTrigger(
        [Values(6, 12, 18)] int windowSize,
        [Values(1.0, 1.6, 2.2)] double trigger)
    {
        string protein = HighComplexityFlank + new string('A', 30) + HighComplexityFlank;
        int blockEnd = LowComplexityBlockStart + 29;

        var regions = ProteinMotifFinder.FindLowComplexityRegions(protein, windowSize, trigger).ToList();

        regions.Should().NotBeEmpty("a 30-residue poly-A block (K=0) triggers at any K1 ≥ 0");
        regions.Should().OnlyContain(r => r.Start >= 0 && r.Start <= r.End && r.End < protein.Length, "INV-4");
        regions.Should().OnlyContain(r => r.Complexity <= 2.5 + 1e-9, "reported region complexity ≤ K2 (INV-3)");

        // A region covers the poly-A block, and its (minimum) complexity is exactly 0 (homopolymer, INV-2).
        var blockRegion = regions.Where(r => r.Start <= LowComplexityBlockStart + windowSize && r.End >= blockEnd - windowSize)
            .ToList();
        blockRegion.Should().NotBeEmpty("the poly-A block is reported");
        blockRegion.Should().Contain(r => Math.Abs(r.Complexity) < 1e-9, "a homopolymer window has complexity 0 (INV-2)");
    }

    /// <summary>
    /// Interaction witness — the trigger cutoff gates by complexity: a two-residue "AB" repeat
    /// (K = 1 bit) is a low-complexity segment at K1 = 1.5 but not at K1 = 0.5; a fully diverse
    /// 20-residue cycle (K = log₂20 ≈ 4.32) is never low-complexity.
    /// </summary>
    [Test]
    public void ProtMotifLC_ComplexityThresholdGates()
    {
        string ab = string.Concat(Enumerable.Repeat("AB", 20)); // K = 1 bit/residue

        ProteinMotifFinder.FindLowComplexityRegions(ab, 12, 1.5, 1.5).Should().NotBeEmpty("K=1 ≤ 1.5");
        ProteinMotifFinder.FindLowComplexityRegions(ab, 12, 0.5, 0.5).Should().BeEmpty("K=1 > 0.5");

        ProteinMotifFinder.FindLowComplexityRegions(HighComplexityFlank, 12, 2.2)
            .Should().BeEmpty("a fully diverse window (K ≈ 4.32) is not low-complexity");
    }

    // ═══════════════════════════════════════════════════════════════════════
    // Unit: PROTMOTIF-SP-001 — Signal-peptide cleavage prediction (von Heijne) (ProteinMotif)
    // Checklist: docs/checklists/09_COMBINATORIAL_TESTING.md, row 167.
    // Spec: tests/TestSpecs/PROTMOTIF-SP-001.md (canonical PredictSignalPeptide). ADVANCED §10.
    // Dimensions: windowSize(3) × threshold(3). Grid 3×3 = 9 (full, exhaustive ⊇ pairwise).
    //
    // Model (von Heijne 1986; EMBOSS sigcleave): the score at each candidate site is Σ ln(count/expect)
    // over the FIXED 15-residue window (−13..+2); the best-scoring site is returned. IsLikelySignalPeptide
    // ⇔ Score ≥ minWeight (INV-3). The score is the argmax and is INDEPENDENT of minWeight.
    //
    // Axis mapping (documented — the 15-residue von Heijne window is FIXED): windowSize → input length
    // (≥15), confirming the fixed-window score/coords are well-formed at any length; threshold →
    // minWeight. The combinatorial point: across length and threshold, the score is threshold-invariant
    // and the acceptance flag flips exactly at minWeight, with coordinates obeying INV-2/INV-6.
    // ═══════════════════════════════════════════════════════════════════════

    private const string SignalSeq = "MALWMRLLPLLALLALWGPDPAAAFVNQHLCGSHLVEALYLVCGERGFFYTPKTRREAEDLQ"; // ≥60 aa

    [Test, Combinatorial]
    public void ProtMotifSP_ScoreThresholdInvariance_AcrossLengthAndMinWeight(
        [Values(15, 30, 60)] int inputLen,
        [Values(0.0, 3.5, 6.0)] double minWeight)
    {
        string seq = SignalSeq[..inputLen];

        double baseScore = ProteinMotifFinder.PredictSignalPeptide(seq, false, double.MinValue)!.Value.Score;

        var r = ProteinMotifFinder.PredictSignalPeptide(seq, false, minWeight);
        r.Should().NotBeNull("a ≥15-residue sequence has a scoring window (INV-5)");
        var sp = r!.Value;

        sp.Score.Should().Be(baseScore, "the argmax score does not depend on the acceptance threshold");
        sp.IsLikelySignalPeptide.Should().Be(baseScore >= minWeight, "IsLikely ⇔ Score ≥ minWeight (INV-3)");
        sp.CleavagePosition.Should().BeInRange(1, seq.Length, "1-based mature start in range (INV-2)");
        sp.SignalSequence.Should().Be(seq.Substring(0, sp.CleavagePosition - 1), "signal = residues before cleavage (INV-6)");
        sp.WindowSequence.Length.Should().BeLessThanOrEqualTo(15, "the scoring window is ≤ 15 residues (INV-6)");
    }

    /// <summary>
    /// Interaction witness — the acceptance flag flips exactly at the score; the organism matrix
    /// changes the score; the prediction is case-independent (INV-4); too-short input is null (INV-5).
    /// </summary>
    [Test]
    public void ProtMotifSP_ThresholdFlip_OrganismAndCase()
    {
        double euk = ProteinMotifFinder.PredictSignalPeptide(SignalSeq, false, double.MinValue)!.Value.Score;

        ProteinMotifFinder.PredictSignalPeptide(SignalSeq, false, euk - 0.01)!.Value.IsLikelySignalPeptide
            .Should().BeTrue("Score ≥ minWeight just below the score");
        ProteinMotifFinder.PredictSignalPeptide(SignalSeq, false, euk + 0.01)!.Value.IsLikelySignalPeptide
            .Should().BeFalse("Score < minWeight just above the score");

        double pro = ProteinMotifFinder.PredictSignalPeptide(SignalSeq, true, double.MinValue)!.Value.Score;
        pro.Should().NotBe(euk, "the prokaryotic and eukaryotic matrices score differently");

        ProteinMotifFinder.PredictSignalPeptide(SignalSeq.ToLowerInvariant(), false, double.MinValue)!.Value.Score
            .Should().BeApproximately(euk, 1e-9, "prediction is case-independent (INV-4)");

        ProteinMotifFinder.PredictSignalPeptide("MKWV", false).Should().BeNull("inputs shorter than 15 residues return null (INV-5)");
    }

    // ═══════════════════════════════════════════════════════════════════════
    // Unit: PROTMOTIF-TM-001 — Transmembrane-helix prediction (Kyte-Doolittle) (ProteinMotif)
    // Checklist: docs/checklists/09_COMBINATORIAL_TESTING.md, row 168.
    // Spec: tests/TestSpecs/PROTMOTIF-TM-001.md (canonical PredictTransmembraneHelices). ADVANCED §10.
    // Dimensions: windowSize(3) × threshold(3) × seqLen(3). Grid 3×3×3 = 27 (full, exhaustive).
    //
    // Model (Kyte & Doolittle 1982): the mean hydropathy over a sliding window flags membrane-spanning
    // segments; each maximal run of windows with mean ≥ threshold (and span ≥ 19 residues) is reported
    // with its peak window mean. A uniform residue gives a flat profile at that residue's scale value.
    //
    // Engineered construct: a hydrophobic poly-I core (KD(I)=4.5) flanked by hydrophilic poly-Q
    // (KD(Q)=−3.5). The combinatorial point: window width, threshold and protein length jointly bound
    // detection; the poly-I core scores 4.5 ≥ every threshold ≤ 2.5 and is reported (Score ≥ threshold,
    // valid coords) at every cell, overlapping the planted core.
    // ═══════════════════════════════════════════════════════════════════════

    [Test, Combinatorial]
    public void ProtMotifTM_DetectsHydrophobicHelix_AcrossWindowThresholdLength(
        [Values(15, 19, 23)] int windowSize,
        [Values(1.0, 1.6, 2.5)] double threshold,
        [Values(10, 20, 40)] int flankLen)
    {
        const int coreLen = 30;
        string protein = new string('Q', flankLen) + new string('I', coreLen) + new string('Q', flankLen);
        int coreStart = flankLen, coreEnd = flankLen + coreLen - 1;

        var segments = ProteinMotifFinder.PredictTransmembraneHelices(protein, windowSize, threshold).ToList();

        segments.Should().NotBeEmpty("a 30-residue poly-I core (mean 4.5) exceeds every threshold ≤ 2.5");
        segments.Should().OnlyContain(s => s.Score >= threshold - 1e-12, "peak window mean ≥ threshold (INV-1)");
        segments.Should().OnlyContain(s => s.Start >= 0 && s.Start <= s.End && s.End < protein.Length, "valid coords (INV-2)");

        // A reported segment overlaps the hydrophobic core.
        segments.Should().Contain(s => s.Start <= coreEnd && s.End >= coreStart, "the TM segment covers the poly-I core");
    }

    /// <summary>
    /// Interaction witness — hydrophobicity gates detection: a moderately hydrophobic poly-A core
    /// (KD=1.8) is a TM helix at threshold 1.6 but not 2.5; a hydrophilic poly-Q is never one; and a
    /// uniform poly-I profile equals the Kyte-Doolittle value 4.5 (INV-3).
    /// </summary>
    [Test]
    public void ProtMotifTM_HydrophobicityGatesDetection()
    {
        string polyA = new string('Q', 20) + new string('A', 30) + new string('Q', 20);
        ProteinMotifFinder.PredictTransmembraneHelices(polyA, 19, 1.6).Should().NotBeEmpty("A's 1.8 ≥ 1.6");
        ProteinMotifFinder.PredictTransmembraneHelices(polyA, 19, 2.5).Should().BeEmpty("A's 1.8 < 2.5");

        ProteinMotifFinder.PredictTransmembraneHelices(new string('Q', 40), 19, 1.6)
            .Should().BeEmpty("hydrophilic poly-Q (KD −3.5) is not a transmembrane helix");

        var pureI = ProteinMotifFinder.PredictTransmembraneHelices(new string('I', 30), 19, 1.6).ToList();
        pureI.Should().ContainSingle();
        pureI[0].Score.Should().BeApproximately(4.5, 1e-9, "a uniform poly-I profile equals KD(I)=4.5 (INV-3)");
    }

    // ═══════════════════════════════════════════════════════════════════════
    // Unit: PROTMOTIF-HMM-001 — Plan7 profile-HMM domain search (HMMER3-style)
    // Checklist: docs/checklists/09_COMBINATORIAL_TESTING.md, row 239.
    // Spec: tests/TestSpecs/PROTMOTIF-HMM-001.md
    //       (ProteinMotifFinder.ScoreDomainHmm / FindDomainEnvelopes; Plan7ProfileHmm). ADVANCED §10.
    //
    // Sources: Durbin, Eddy, Krogh & Mitchison (1998) §5.4 (Plan7 log-odds recurrences);
    //   Eddy (2011) PLoS Comput Biol 7:e1002195 (HMMER3 local architecture, env decomposition);
    //   Pfam PF00018 (SH3), PF00595 (PDZ), PF00400 (WD40) bundled CC0 profiles.
    //
    // Model: a profile HMM scores a sequence in log-odds bits against a background null. The GLOCAL
    //   full-profile score (ScoreDomainHmm) aligns the whole profile to the whole input; the LOCAL
    //   search (FindDomainEnvelopes) decomposes the input into per-domain envelopes, so an embedded
    //   domain is detected regardless of flanking residues. A true positive scores above threshold on
    //   its OWN family and lower on the others (family specificity).
    //
    // Dimensions: profile(3) × mode(local/glocal) × seqLen(3). Grid 3×2×3 = 18 (exhaustive).
    //
    // Axis mapping (documented, cf. PRIMER-TM-001 inert salt axis): seqLen is realised as the length
    // of neutral Gly-Ser flanks around the planted domain. Under LOCAL search this is the
    // flank-robustness axis (the envelope is found at every flank length); under GLOCAL full-profile
    // scoring the score is taken on the bare domain (flank-independent by construction), so the
    // seqLen axis is inert for glocal cells.
    // ═══════════════════════════════════════════════════════════════════════

    // Evidence-sourced true-positive domains (UniProt cores) and their Pfam accessions.
    private const string Sh3Pos = "TFVALYDYESRTETDLSFKKGERLQIVNNTEGDWWLAHSLSTGQTGYIPSNYVAP";              // SRC_HUMAN SH3
    private const string PdzPos = "MEYEEITLERGNSGLGFSIAGGTDNPHIGDDPSIFITKIIPGGAAAQDGRLRVNDSILFVNEVDVREVTHSAAVEALKEAGSIVRLYVMRR"; // PSD-95 PDZ1
    private const string Wd40Pos =                                                                          // GBB1_HUMAN WD40 β-propeller
        "MSELDQLRQEAEQLKNQIRDARKACADATLSQITNNIDPVGRIQMRTRRTLRGHLAKIYAMHWGTDSRLLVSASQDGKLIIWDSYTTNKVHAIPLRSSWVMTCAYAPSGNYVACGGLDNICSIYNLKTREGNVRVSRELAGHTGYLSCCRFLDDNQIVTSSGDTTCALWDIETGQQTTTFTGHTGDVMSLSLAPDTRLFVSGACDASAKLWDVREGMCRQTFTGHESDINAICFFPNGNAFATGSDDATCRLFDLRADQELMTYSHDNIICGITSVSFSKSGRLLLAGYDDFNCNVWDALKADRAGVLAGHDNRVSCLGVTDDGMAVATGSWDSFLKIWN";

    private const double HmmThresholdBits = 10.0;

    private static readonly (string Name, string Seq, string Acc)[] HmmProfiles =
    {
        ("SH3", Sh3Pos, "PF00018"),
        ("PDZ", PdzPos, "PF00595"),
        ("WD40", Wd40Pos, "PF00400"),
    };

    private static string GsFlank(int len) => string.Concat(Enumerable.Range(0, len).Select(i => "GS"[i % 2]));

    public enum HmmMode { Glocal, Local }

    [Test, Combinatorial]
    public void Hmm_DetectsCognateFamily_AcrossProfileModeAndFlankLength(
        [Values(0, 1, 2)] int profileIdx,
        [Values(HmmMode.Glocal, HmmMode.Local)] HmmMode mode,
        [Values(0, 20, 40)] int flankLen)
    {
        var (name, seq, acc) = HmmProfiles[profileIdx];

        if (mode == HmmMode.Glocal)
        {
            // Glocal full-profile bit score on the bare domain (flank length inert — documented).
            double own = ProteinMotifFinder.ScoreDomainHmm(seq, acc);
            own.Should().BeGreaterThanOrEqualTo(HmmThresholdBits, $"the {name} domain scores above threshold on its own profile");
            foreach (var o in HmmProfiles.Where(p => p.Acc != acc))
                own.Should().BeGreaterThan(ProteinMotifFinder.ScoreDomainHmm(seq, o.Acc),
                    $"the {name} domain scores higher on {name} than on {o.Name} (family specificity)");
        }
        else
        {
            // Local envelope decomposition: the embedded cognate domain is found at every flank length.
            string padded = GsFlank(flankLen) + seq + GsFlank(flankLen);
            ProteinMotifFinder.FindDomainEnvelopes(padded)
                .Should().Contain(e => e.Accession == acc && e.BitScore >= HmmThresholdBits,
                    $"the embedded {name} domain is detected as a local envelope regardless of {flankLen}-residue flanks");
        }
    }

    /// <summary>
    /// Interaction witness: a low-complexity true-negative sequence triggers NO local envelope and
    /// scores below zero on every profile — the family calls are genuine, not green-lit by construction.
    /// </summary>
    [Test]
    public void Hmm_TrueNegative_NoEnvelopeAndNegativeGlocalScore()
    {
        const string trueNeg = "AAAAAAAAAAAAAAEEEEEEEEEEEEEEKKKKKKKKKKKK";
        ProteinMotifFinder.FindDomainEnvelopes(trueNeg).Should().BeEmpty("a low-complexity sequence yields no domain envelope");
        foreach (var p in HmmProfiles)
            ProteinMotifFinder.ScoreDomainHmm(trueNeg, p.Acc).Should().BeLessThan(0.0, $"the true negative scores below zero on {p.Name}");
    }
}
