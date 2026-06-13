// META-RESIST-001 — Antibiotic Resistance Gene Detection (ResFinder-style)
// Evidence: docs/Evidence/META-RESIST-001-Evidence.md
// TestSpec: tests/TestSpecs/META-RESIST-001.md
// Source: Zankari et al. (2012), J Antimicrob Chemother 67(11):2640-2644 (ResFinder);
//         Li H (2018) BLAST identity = matches/alignment-columns;
//         ResFinder GitHub README + Sci Rep (2023) thresholds 0.90 ID / 0.60 coverage.

using NUnit.Framework;
using Seqeron.Genomics.Metagenomics;
using System;
using System.Collections.Generic;
using System.Linq;

namespace Seqeron.Genomics.Tests;

[TestFixture]
public class MetagenomicsAnalyzer_FindAntibioticResistanceGenes_Tests
{
    private static IEnumerable<(string, string, string, string)> Db(
        params (string Id, string Seq, string Name, string Class)[] e) => e;

    private static IEnumerable<(string, string)> Contig(string id, string seq)
        => new[] { (id, seq) };

    #region FindAntibioticResistanceGenes — identity & coverage (M1, M2, M4)

    // M1 — Exact full-length match: contig contains the 7-base reference verbatim.
    // identity = 7/7 = 1.0, coverage = 7/7 = 1.0. Evidence: CARD "Perfect" = 100% full length; BLAST identity = matches/columns.
    [Test]
    public void FindAntibioticResistanceGenes_ExactFullLengthMatch_IdentityAndCoverageOne()
    {
        var contigs = Contig("c1", "AAACGTACGT");          // contains "CGTACGT"
        var db = Db(("blaX", "CGTACGT", "blaX-like", "beta-lactam"));

        var hit = MetagenomicsAnalyzer.FindAntibioticResistanceGenes(contigs, db).Single();

        Assert.Multiple(() =>
        {
            Assert.That(hit.ResistanceGene, Is.EqualTo("blaX-like"), "Best-matching reference gene name reported.");
            Assert.That(hit.AntibioticClass, Is.EqualTo("beta-lactam"), "Antibiotic class passed through from the DB entry.");
            Assert.That(hit.PercentIdentity, Is.EqualTo(1.0).Within(1e-10), "Exact match: 7/7 identical positions = 1.0.");
            Assert.That(hit.Coverage, Is.EqualTo(1.0).Within(1e-10), "Full-length match: 7/7 reference covered = 1.0.");
        });
    }

    // M2 — One substitution over a full-length 7-base alignment: identity = 6/7, coverage = 1.0.
    // Evidence: Li (2018) identity = matching bases / alignment columns (gapless ⇒ 7 columns).
    [Test]
    public void FindAntibioticResistanceGenes_SingleMismatchFullLength_IdentitySixSevenths()
    {
        var contigs = Contig("c1", "CGTTCGT");             // ref CGTACGT, position 4 A->T
        var db = Db(("blaX", "CGTACGT", "blaX-like", "beta-lactam"));

        // Threshold below 6/7 so the near-perfect hit is reported and we can read its identity.
        var hit = MetagenomicsAnalyzer
            .FindAntibioticResistanceGenes(contigs, db, identityThreshold: 0.80, coverageThreshold: 0.60)
            .Single();

        Assert.Multiple(() =>
        {
            Assert.That(hit.PercentIdentity, Is.EqualTo(6.0 / 7.0).Within(1e-10), "6 of 7 aligned positions identical = 6/7.");
            Assert.That(hit.Coverage, Is.EqualTo(1.0).Within(1e-10), "Window spans full reference length: 7/7 = 1.0.");
        });
    }

    // M4 — Partial hit above the coverage floor: 5 of 7 reference bases present, exact.
    // coverage = 5/7 ≈ 0.714 ≥ 0.60; identity over the matched window = 5/5 = 1.0.
    // Evidence: Zankari (2012) coverage vs reference length; Sci Rep (2023) edge genes retained ≥ floor.
    [Test]
    public void FindAntibioticResistanceGenes_PartialHitAboveCoverageFloor_Reported()
    {
        var contigs = Contig("c1", "TTCGTAC");             // matches first 5 of CGTACGT ("CGTAC") at the right end
        var db = Db(("blaX", "CGTACGT", "blaX-like", "beta-lactam"));

        var hit = MetagenomicsAnalyzer
            .FindAntibioticResistanceGenes(contigs, db, identityThreshold: 0.90, coverageThreshold: 0.60)
            .Single();

        Assert.Multiple(() =>
        {
            Assert.That(hit.PercentIdentity, Is.EqualTo(1.0).Within(1e-10), "All 5 aligned positions identical = 1.0.");
            Assert.That(hit.Coverage, Is.EqualTo(5.0 / 7.0).Within(1e-10), "5 of 7 reference bases covered = 5/7.");
        });
    }

    #endregion

    #region FindAntibioticResistanceGenes — reporting thresholds (M3, M5)

    // M3 — Contig-edge partial below the coverage floor: only 4 of 7 bases → coverage 4/7 ≈ 0.571 < 0.60.
    // Not reported. Evidence: Zankari (2012) "cover at least 2/5 of the reference"; default 0.60 floor.
    [Test]
    public void FindAntibioticResistanceGenes_PartialBelowCoverageFloor_NotReported()
    {
        var contigs = Contig("c1", "TTTCGTA");             // matches first 4 of CGTACGT ("CGTA")
        var db = Db(("blaX", "CGTACGT", "blaX-like", "beta-lactam"));

        var hits = MetagenomicsAnalyzer
            .FindAntibioticResistanceGenes(contigs, db)    // default 0.90 / 0.60
            .ToList();

        Assert.That(hits, Is.Empty, "Coverage 4/7 ≈ 0.571 is below the 0.60 floor → gene fragment not reported.");
    }

    // M5 — Below the identity threshold: 5/7 ≈ 0.714 identity over a full-length window, threshold 0.90.
    // Not reported. Evidence: Zankari (2012) sub-threshold matches are fragments/noise.
    [Test]
    public void FindAntibioticResistanceGenes_BelowIdentityThreshold_NotReported()
    {
        var contigs = Contig("c1", "CGAACTT");             // ref CGTACGT: positions 3,5,7 differ → 4/7? verify below
        var db = Db(("blaX", "CGTACGT", "blaX-like", "beta-lactam"));

        // CGTACGT vs CGAACTT aligned: C=C,G=G,T!=A,A=A,C!=C? compare pairwise:
        // ref C G T A C G T  /  qry C G A A C T T  -> matches at 1,2,4,5,7 = 5; identity 5/7 ≈ 0.714 < 0.90.
        var hits = MetagenomicsAnalyzer
            .FindAntibioticResistanceGenes(contigs, db)    // default 0.90 identity
            .ToList();

        Assert.That(hits, Is.Empty, "Best identity 5/7 ≈ 0.714 < 0.90 threshold → not reported.");
    }

    #endregion

    #region FindAntibioticResistanceGenes — best-match selection (M6, C1)

    // M6 — Two reference genes match the same contig; only the higher-identity one is reported.
    // Evidence: Zankari (2012) reports the "best-matching gene"; CARD RGI best hit.
    [Test]
    public void FindAntibioticResistanceGenes_TwoMatches_OnlyBestIdentityReported()
    {
        var contigs = Contig("c1", "CGTACGT");             // exact for geneA, 6/7 for geneB
        var db = Db(
            ("a", "CGTACGT", "geneA", "beta-lactam"),      // identity 1.0
            ("b", "CGTTCGT", "geneB", "aminoglycoside"));  // identity 6/7 vs contig

        var hits = MetagenomicsAnalyzer
            .FindAntibioticResistanceGenes(contigs, db, identityThreshold: 0.80, coverageThreshold: 0.60)
            .ToList();

        Assert.Multiple(() =>
        {
            Assert.That(hits, Has.Count.EqualTo(1), "Exactly one best-matching gene reported per contig.");
            Assert.That(hits[0].ResistanceGene, Is.EqualTo("geneA"), "Higher-identity gene (1.0 > 6/7) wins.");
            Assert.That(hits[0].PercentIdentity, Is.EqualTo(1.0).Within(1e-10), "Reported identity is the winner's.");
        });
    }

    // C1 — Identity tie broken by greater coverage.
    // Both genes are exact (identity 1.0); geneFull covers 7/7, genePart covers 5/7. Higher coverage wins.
    [Test]
    public void FindAntibioticResistanceGenes_IdentityTie_HigherCoverageWins()
    {
        var contigs = Contig("c1", "CGTACGT");
        // Both references are an exact gapless match within the contig (identity 1.0), so the
        // winner is decided by coverage: geneFull spans 7/7 of its reference; genePart's longer
        // 9-base reference is only 7/9 present (the trailing "GG" is absent), coverage 7/9 < 1.0.
        var db = Db(
            ("p", "CGTACGTGG", "genePart", "beta-lactam"), // identity 7/7=1.0 over window, coverage 7/9
            ("f", "CGTACGT", "geneFull", "beta-lactam"));   // identity 1.0, coverage 7/7 = 1.0

        var hits = MetagenomicsAnalyzer
            .FindAntibioticResistanceGenes(contigs, db, identityThreshold: 0.90, coverageThreshold: 0.60)
            .ToList();

        Assert.Multiple(() =>
        {
            Assert.That(hits, Has.Count.EqualTo(1), "One best match per contig.");
            Assert.That(hits[0].ResistanceGene, Is.EqualTo("geneFull"), "Identity tie (both 1.0) broken by higher coverage (7/7 > 7/9).");
            Assert.That(hits[0].Coverage, Is.EqualTo(1.0).Within(1e-10), "Winner's coverage is the full-length 1.0.");
        });
    }

    #endregion

    #region FindAntibioticResistanceGenes — defaults (M7)

    // M7 — Default thresholds equal the ResFinder operating values 0.90 ID / 0.60 coverage.
    // Evidence: ResFinder web service / GitHub README; Sci Rep (2023); JAC (2016).
    [Test]
    public void FindAntibioticResistanceGenes_DefaultThresholds_MatchResFinderValues()
    {
        Assert.Multiple(() =>
        {
            Assert.That(MetagenomicsAnalyzer.DefaultResistanceIdentityThreshold, Is.EqualTo(0.90).Within(1e-12),
                "Default identity threshold is the ResFinder web-service value 0.90.");
            Assert.That(MetagenomicsAnalyzer.DefaultResistanceCoverageThreshold, Is.EqualTo(0.60).Within(1e-12),
                "Default coverage threshold is the ResFinder value 0.60.");
        });
    }

    #endregion

    #region FindAntibioticResistanceGenes — validation & empty (S1–S5)

    // S1 — Null contigs throws ArgumentNullException.
    [Test]
    public void FindAntibioticResistanceGenes_NullContigs_Throws()
    {
        var db = Db(("a", "CGTACGT", "geneA", "beta-lactam"));
        Assert.Throws<ArgumentNullException>(
            () => MetagenomicsAnalyzer.FindAntibioticResistanceGenes(null!, db).ToList(),
            "Null contigs is a contract violation.");
    }

    // S2 — Null referenceGenes throws ArgumentNullException.
    [Test]
    public void FindAntibioticResistanceGenes_NullReferenceGenes_Throws()
    {
        var contigs = Contig("c1", "CGTACGT");
        Assert.Throws<ArgumentNullException>(
            () => MetagenomicsAnalyzer.FindAntibioticResistanceGenes(contigs, null!).ToList(),
            "Null reference database is a contract violation.");
    }

    // S3 — Identity threshold above 1 throws ArgumentOutOfRangeException.
    [Test]
    public void FindAntibioticResistanceGenes_IdentityThresholdOutOfRange_Throws()
    {
        var contigs = Contig("c1", "CGTACGT");
        var db = Db(("a", "CGTACGT", "geneA", "beta-lactam"));
        Assert.Throws<ArgumentOutOfRangeException>(
            () => MetagenomicsAnalyzer.FindAntibioticResistanceGenes(contigs, db, identityThreshold: 1.5).ToList(),
            "Identity threshold must be within [0, 1].");
    }

    // S4 — Coverage threshold below 0 throws ArgumentOutOfRangeException.
    [Test]
    public void FindAntibioticResistanceGenes_CoverageThresholdOutOfRange_Throws()
    {
        var contigs = Contig("c1", "CGTACGT");
        var db = Db(("a", "CGTACGT", "geneA", "beta-lactam"));
        Assert.Throws<ArgumentOutOfRangeException>(
            () => MetagenomicsAnalyzer.FindAntibioticResistanceGenes(contigs, db, coverageThreshold: -0.1).ToList(),
            "Coverage threshold must be within [0, 1].");
    }

    // S5 — Empty contig sequence and a non-matching contig both yield no hits.
    [Test]
    public void FindAntibioticResistanceGenes_EmptyOrNonMatching_ReturnsEmpty()
    {
        var db = Db(("a", "CGTACGT", "geneA", "beta-lactam"));

        var emptyHits = MetagenomicsAnalyzer.FindAntibioticResistanceGenes(Contig("c1", ""), db).ToList();
        var noMatchHits = MetagenomicsAnalyzer.FindAntibioticResistanceGenes(Contig("c2", "AAAAAAAAAA"), db).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(emptyHits, Is.Empty, "Empty contig sequence is skipped.");
            Assert.That(noMatchHits, Is.Empty, "No reference passes thresholds against an unrelated sequence.");
        });
    }

    #endregion
}
