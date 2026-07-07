// CODON-STATS-001 — Codon Usage Statistics (GetStatistics, CalculateCai, reference tables)
// Evidence: docs/Evidence/CODON-STATS-001-Evidence.md
// TestSpec: tests/TestSpecs/CODON-STATS-001.md
// Source: Sharp PM, Li W-H (1987). Nucleic Acids Res. 15(3):1281-1295 (CAI);
//         Peden JF (1999). CodonW thesis §1.8.2.1.3 (GC3s);
//         Biopython v1.79 SharpEcoliIndex (E. coli w); Kazusa H. sapiens (human RSCU).

namespace Seqeron.Genomics.Tests.Unit.MolTools;

[TestFixture]
public class CodonUsageAnalyzer_GetStatistics_Tests
{
    private const double Tol = 1e-10;

    #region CalculateCai

    // M1 — Sharp & Li 1987: CAI of an all-optimal E. coli sequence (each codon = its family's
    // w-max, w=1) is exactly 1.0. Sequence: Leu(CTG) Ile(ATC) Val(GTT) Ala(GCT) Arg(CGT) Lys(AAA).
    [Test]
    public void CalculateCai_AllOptimalEColiCodons_ReturnsOne()
    {
        double cai = CodonUsageAnalyzer.CalculateCai("CTGATCGTTGCTCGTAAA", CodonUsageAnalyzer.EColiOptimalCodons);
        Assert.That(cai, Is.EqualTo(1.0).Within(Tol),
            "All six codons are the most-adaptive synonym in E. coli (w=1) so the geometric mean is 1.0.");
    }

    // M2 — Sharp & Li 1987 geometric-mean derivation: GCTGCC = Ala GCT (w=1) + Ala GCC (w=0.122);
    // CAI = sqrt(1 * 0.122) = 0.34928498393146.
    [Test]
    public void CalculateCai_TwoAlaCodons_EqualsGeometricMeanOfW()
    {
        double cai = CodonUsageAnalyzer.CalculateCai("GCTGCC", CodonUsageAnalyzer.EColiOptimalCodons);
        Assert.That(cai, Is.EqualTo(Math.Sqrt(1.0 * 0.122)).Within(Tol),
            "CAI is the geometric mean of w(GCT)=1 and w(GCC)=0.122 per Sharp & Li 1987.");
    }

    // M3 — Sharp & Li 1987: suboptimal codons CTA(0.007) ATA(0.003) GTC(0.066);
    // CAI = cube root(0.007*0.003*0.066) = 0.01114947479545.
    [Test]
    public void CalculateCai_RareEColiCodons_EqualsCubeRootProduct()
    {
        double cai = CodonUsageAnalyzer.CalculateCai("CTAATAGTC", CodonUsageAnalyzer.EColiOptimalCodons);
        double expected = Math.Cbrt(0.007 * 0.003 * 0.066);
        Assert.That(cai, Is.EqualTo(expected).Within(Tol),
            "CAI is the geometric mean of the three rare codons' w values (Sharp & Li 1987).");
    }

    // M4 — seqinr / CodonW: CAI excludes single-codon amino acids (Met, Trp) and stop codons.
    // ATG(Met) TGG(Trp) TAA(stop) leaves zero scorable codons -> CAI = 0.
    [Test]
    public void CalculateCai_OnlyMetTrpStop_ReturnsZero()
    {
        double cai = CodonUsageAnalyzer.CalculateCai("ATGTGGTAA", CodonUsageAnalyzer.EColiOptimalCodons);
        Assert.That(cai, Is.EqualTo(0.0).Within(Tol),
            "Met, Trp and stop codons are excluded from CAI; no scorable codon remains so CAI is 0.");
    }

    // S6 / INV-1 — CAI is bounded in [0,1] for any sequence, and equals the exact geometric
    // mean of the scorable codons' E. coli w values (ATG=Met excluded). Sequence
    // ATG AAA TTT GGG CTG GTT AAA CGT -> scorable {AAA=1, TTT=0.296, GGG=0.019, CTG=1,
    // GTT=1, AAA=1, CGT=1}; CAI = (1*0.296*0.019*1*1*1*1)^(1/7) = 0.47706538020472955.
    // w values from Biopython SharpEcoliIndex (Sharp & Li 1987).
    [Test]
    public void CalculateCai_ArbitrarySequence_EqualsExactGeometricMeanAndIsBetweenZeroAndOne()
    {
        double cai = CodonUsageAnalyzer.CalculateCai("ATGAAATTTGGGCTGGTTAAACGT", CodonUsageAnalyzer.EColiOptimalCodons);
        double expected = Math.Pow(1.0 * 0.296 * 0.019 * 1.0 * 1.0 * 1.0 * 1.0, 1.0 / 7.0);
        Assert.Multiple(() =>
        {
            Assert.That(cai, Is.EqualTo(expected).Within(Tol),
                "CAI is the geometric mean over the 7 scorable codons (Met ATG excluded); Sharp & Li 1987.");
            Assert.That(cai, Is.GreaterThanOrEqualTo(0.0), "CAI is a geometric mean of non-negative weights.");
            Assert.That(cai, Is.LessThanOrEqualTo(1.0), "Each w <= 1, so the geometric mean <= 1 (Sharp & Li 1987).");
        });
    }

    // S2 — empty string returns CAI 0 (input contract).
    [Test]
    public void CalculateCai_EmptyString_ReturnsZero()
    {
        double cai = CodonUsageAnalyzer.CalculateCai("", CodonUsageAnalyzer.EColiOptimalCodons);
        Assert.That(cai, Is.EqualTo(0.0).Within(Tol), "Empty input has no codons; CAI is defined as 0.");
    }

    // S4 — null DnaSequence throws.
    [Test]
    public void CalculateCai_NullDnaSequence_Throws()
    {
        Assert.Throws<ArgumentNullException>(
            () => CodonUsageAnalyzer.CalculateCai((DnaSequence)null!, CodonUsageAnalyzer.EColiOptimalCodons),
            "A null DnaSequence must raise ArgumentNullException.");
    }

    // S5 — null reference table throws.
    [Test]
    public void CalculateCai_NullReference_Throws()
    {
        Assert.Throws<ArgumentNullException>(
            () => CodonUsageAnalyzer.CalculateCai(new DnaSequence("ATGAAA"), null!),
            "A null reference table must raise ArgumentNullException.");
    }

    #endregion

    #region GetStatistics — GC positions and GC3s

    // M5 — Peden §1.8.2.1.3: GC3s excludes Met/Trp/stop. ATG(Met, excluded) GCA(Ala, 3rd=A).
    // GC3s denominator = {GCA}, numerator 0 -> GC3s = 0; while GC3 over ALL positions = 50%
    // (G of ATG, A of GCA -> 1/2). This is the GC3s != GC3 discriminating case.
    [Test]
    public void GetStatistics_MetPlusAla_Gc3sExcludesMetButGc3DoesNot()
    {
        var stats = CodonUsageAnalyzer.GetStatistics("ATGGCA");
        Assert.Multiple(() =>
        {
            Assert.That(stats.Gc3s, Is.EqualTo(0.0).Within(Tol),
                "Only Ala GCA is synonymous; its 3rd base A is not G/C, so GC3s = 0 (Peden §1.8.2.1.3).");
            Assert.That(stats.Gc3, Is.EqualTo(50.0).Within(Tol),
                "GC3 over all positions counts ATG's G and GCA's A -> 1/2 = 50%, unlike GC3s.");
        });
    }

    // M6 — Peden §1.8.2.1.3: GC3s = fraction of synonymous codons with G/C at the 3rd position.
    // GCC(Ala, 3rd=C -> GC) and GCA(Ala, 3rd=A -> not) -> GC3s = 1/2 = 50%.
    [Test]
    public void GetStatistics_TwoAlaCodons_Gc3sIsSynonymousFraction()
    {
        var stats = CodonUsageAnalyzer.GetStatistics("GCCGCA");
        Assert.That(stats.Gc3s, Is.EqualTo(50.0).Within(Tol),
            "Two synonymous Ala codons; one has G/C at position 3, so GC3s = 50%.");
    }

    // M7 — EMBOSS cusp: GC1/GC2/GC3 are per-position G/C content. CTG GTT AAA:
    // pos1 = C,G,A -> 2/3; pos2 = T,T,A -> 0/3; pos3 = G,T,A -> 1/3.
    [Test]
    public void GetStatistics_ThreeCodons_GcPerPositionExact()
    {
        var stats = CodonUsageAnalyzer.GetStatistics("CTGGTTAAA");
        Assert.Multiple(() =>
        {
            Assert.That(stats.Gc1, Is.EqualTo(200.0 / 3.0).Within(Tol), "Position 1: C,G,A -> 2/3 are G/C = 66.6667%.");
            Assert.That(stats.Gc2, Is.EqualTo(0.0).Within(Tol), "Position 2: T,T,A -> none are G/C = 0%.");
            Assert.That(stats.Gc3, Is.EqualTo(100.0 / 3.0).Within(Tol), "Position 3: G,T,A -> 1/3 are G/C = 33.3333%.");
        });
    }

    // C1 / INV-6 — OverallGc = (Gc1+Gc2+Gc3)/3.
    [Test]
    public void GetStatistics_OverallGc_IsAverageOfThreePositions()
    {
        var stats = CodonUsageAnalyzer.GetStatistics("CTGGTTAAA");
        Assert.That(stats.OverallGc, Is.EqualTo((stats.Gc1 + stats.Gc2 + stats.Gc3) / 3.0).Within(Tol),
            "OverallGc is defined as the mean of the three positional GC percentages.");
    }

    #endregion

    #region GetStatistics — counts, RSCU, totals, edge cases

    // M8 — EMBOSS cusp "Number": TotalCodons = valid in-frame codons; per-codon counts.
    // CTGCTGGTTAAA -> CTG x2, GTT x1, AAA x1 -> 4 codons.
    [Test]
    public void GetStatistics_RepeatedCodons_CountsAndTotalExact()
    {
        var stats = CodonUsageAnalyzer.GetStatistics("CTGCTGGTTAAA");
        Assert.Multiple(() =>
        {
            Assert.That(stats.TotalCodons, Is.EqualTo(4), "Four valid in-frame codons present.");
            Assert.That(stats.CodonCounts["CTG"], Is.EqualTo(2), "CTG occurs twice.");
            Assert.That(stats.CodonCounts["GTT"], Is.EqualTo(1), "GTT occurs once.");
            Assert.That(stats.CodonCounts["AAA"], Is.EqualTo(1), "AAA occurs once.");
        });
    }

    // M11 — Sharp et al. 1986 RSCU: with only Leu codon CTG present (Leu is 6-fold),
    // RSCU(CTG) = n*x / sum = 6*2/2 = 6.0.
    [Test]
    public void GetStatistics_OnlyOneLeuCodon_RscuEqualsFamilySize()
    {
        var stats = CodonUsageAnalyzer.GetStatistics("CTGCTG");
        Assert.That(stats.Rscu["CTG"], Is.EqualTo(6.0).Within(Tol),
            "Leu is 6-fold; if only CTG is used, RSCU = 6 * x / sum = 6 (Sharp et al. 1986).");
    }

    // S7 — trailing partial codon (< 3 nt) is ignored. GCCGCAG (7 nt) -> 2 codons.
    [Test]
    public void GetStatistics_TrailingPartialCodon_Ignored()
    {
        var stats = CodonUsageAnalyzer.GetStatistics("GCCGCAG");
        Assert.That(stats.TotalCodons, Is.EqualTo(2),
            "The 7th nucleotide forms no complete codon and is not counted.");
    }

    // C2 — input is case-insensitive (normalized to upper-case).
    [Test]
    public void GetStatistics_Lowercase_MatchesUppercase()
    {
        var lower = CodonUsageAnalyzer.GetStatistics("ctgctg");
        var upper = CodonUsageAnalyzer.GetStatistics("CTGCTG");
        Assert.Multiple(() =>
        {
            Assert.That(lower.TotalCodons, Is.EqualTo(upper.TotalCodons), "Codon count is case-insensitive.");
            Assert.That(lower.CodonCounts["CTG"], Is.EqualTo(2), "Lowercase 'ctg' is normalized to 'CTG'.");
        });
    }

    // INV-5 / doc §6.1 — a codon containing a non-ACGT character is skipped entirely (not
    // counted in TotalCodons, counts, or GC positions). EMBOSS cusp counts only valid codons.
    // CTG NNN GTT -> the middle codon is invalid; only CTG and GTT remain (2 codons).
    [Test]
    public void GetStatistics_NonAcgtCodon_IsSkipped()
    {
        var stats = CodonUsageAnalyzer.GetStatistics("CTGNNNGTT");
        Assert.Multiple(() =>
        {
            Assert.That(stats.TotalCodons, Is.EqualTo(2), "The NNN codon is non-ACGT and is skipped.");
            Assert.That(stats.CodonCounts.ContainsKey("NNN"), Is.False, "Invalid codons are not recorded.");
            Assert.That(stats.CodonCounts["CTG"], Is.EqualTo(1), "CTG is the only valid Leu codon counted.");
            Assert.That(stats.CodonCounts["GTT"], Is.EqualTo(1), "GTT is counted.");
        });
    }

    // S1 — empty string returns a fully zeroed statistics record.
    [Test]
    public void GetStatistics_EmptyString_ReturnsZeroedStatistics()
    {
        var stats = CodonUsageAnalyzer.GetStatistics("");
        Assert.Multiple(() =>
        {
            Assert.That(stats.TotalCodons, Is.EqualTo(0), "No codons in empty input.");
            Assert.That(stats.CodonCounts, Is.Empty, "No codon counts in empty input.");
            Assert.That(stats.Rscu, Is.Empty, "No RSCU values in empty input.");
            Assert.That(stats.Enc, Is.EqualTo(0.0).Within(Tol), "ENC is 0 for empty input.");
            Assert.That(stats.Gc1, Is.EqualTo(0.0).Within(Tol), "GC1 is 0 for empty input.");
            Assert.That(stats.Gc3s, Is.EqualTo(0.0).Within(Tol), "GC3s is 0 for empty input.");
        });
    }

    // S3 — null DnaSequence throws.
    [Test]
    public void GetStatistics_NullDnaSequence_Throws()
    {
        Assert.Throws<ArgumentNullException>(
            () => CodonUsageAnalyzer.GetStatistics((DnaSequence)null!),
            "A null DnaSequence must raise ArgumentNullException.");
    }

    #endregion

    #region Reference Tables (E. coli / human)

    // M9 — Biopython SharpEcoliIndex (Sharp & Li 1987) exact w values.
    [Test]
    public void EColiOptimalCodons_MatchSharpLi1987WValues()
    {
        var t = CodonUsageAnalyzer.EColiOptimalCodons;
        Assert.Multiple(() =>
        {
            Assert.That(t["CTG"], Is.EqualTo(1.0).Within(Tol), "Leu CTG is the optimal codon (w=1).");
            Assert.That(t["GCC"], Is.EqualTo(0.122).Within(Tol), "Ala GCC w=0.122 (Sharp & Li 1987).");
            Assert.That(t["GCT"], Is.EqualTo(1.0).Within(Tol), "Ala GCT is the optimal codon (w=1).");
            Assert.That(t["CGT"], Is.EqualTo(1.0).Within(Tol), "Arg CGT is the optimal codon (w=1).");
            Assert.That(t["AGG"], Is.EqualTo(0.002).Within(Tol), "Arg AGG w=0.002 (Sharp & Li 1987).");
            Assert.That(t["TTT"], Is.EqualTo(0.296).Within(Tol), "Phe TTT w=0.296 (Sharp & Li 1987).");
            Assert.That(t.Count, Is.EqualTo(64), "All 64 codons are present (stops as 0.0).");
        });
    }

    // M10 — Kazusa-derived human RSCU exact values.
    [Test]
    public void HumanOptimalCodons_MatchKazusaDerivedRscu()
    {
        var t = CodonUsageAnalyzer.HumanOptimalCodons;
        Assert.Multiple(() =>
        {
            Assert.That(t["CTG"], Is.EqualTo(2.3713).Within(Tol), "Human Leu CTG RSCU from Kazusa frequencies.");
            Assert.That(t["GCC"], Is.EqualTo(1.5988).Within(Tol), "Human Ala GCC RSCU from Kazusa frequencies.");
            Assert.That(t["GTG"], Is.EqualTo(1.8517).Within(Tol), "Human Val GTG RSCU from Kazusa frequencies.");
            Assert.That(t["ATG"], Is.EqualTo(1.0).Within(Tol), "Met is single-codon: RSCU = 1.");
            Assert.That(t["TGG"], Is.EqualTo(1.0).Within(Tol), "Trp is single-codon: RSCU = 1.");
            Assert.That(t.Count, Is.EqualTo(64), "All 64 codons are present.");
        });
    }

    #endregion
}
