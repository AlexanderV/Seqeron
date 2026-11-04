// CODON-RARE-001 — review 2026-10 (FIN-B02 WP7): Sherlocc threshold comparison (A1-1) and
// per-thousand organism presets for %MinMax (A1-3).
// Evidence: docs/Evidence/CODON-RARE-001-Evidence.md; TestSpec: tests/TestSpecs/CODON-RARE-001.md
// Sources:
//   Sherlocc 1.0 (Chartier & Najmanovich; github.com/mtthchrtr/sherlocc, commit b3b4d15):
//     sherlocc.pl  `if($avrg <= $avrg_lim){ #We have a pause`; README "positions with a codon usage
//     frequency average equal or lower than the specified threshold ... considered as 'slow' positions".
//   Kazusa species=316407 per thousand: https://www.kazusa.or.jp/codon/cgi-bin/showcodon.cgi?species=316407
//     ("Escherichia coli W3110 [gbbct]: 4332 CDS's (1372057 codons)"), fetched 2026-10-14 (SRC RO1).
//   Kazusa species=4932 per thousand: CHARMING ScerCUB.txt (github.com/wrightgs/CHARMING, 2c75d1b).
//   Kazusa species=9606 counts (40 662 582 codons): PyPI codon-optimizer 0.2.0
//     reference_seq/Homo_sapiens_codon_frequency.csv.
//   %MinMax expected values: Clark-lab calculateMinMax (CHARMING.py; prints 2 dp) and an independent
//     Python port of the Clarke & Clark (2008) formula over the same tables (unrounded values).

namespace Seqeron.Genomics.Tests.Unit.MolTools;

[TestFixture]
public class CodonOptimizer_SherloccThreshold_PerThousandPresets_Tests
{
    // E. coli K12 preset (Kazusa 316407): UCU, UCC, UCG, AGU, GUA, GGG are exactly 0.15; CGC 0.40.
    private const string SixAtThresholdPlusCommon = "UCUUCCUCGAGUGUAGGGCGC";

    #region A1-1 — Sherlocc threshold comparison

    [Test]
    public void FindRareCodonClusters_DefaultIsSherloccInclusive_CodonsAtThresholdArePauses()
    {
        var clusters = CodonOptimizer.FindRareCodonClusters(SixAtThresholdPlusCommon, CodonOptimizer.EColiK12);
        Assert.That(clusters, Is.EqualTo(new[] { new CodonOptimizer.RareCodonCluster(0, 6, 6) }));
    }

    [Test]
    public void FindRareCodonClusters_StrictOption_ReproducesPreviousBehaviour()
    {
        var strict = CodonOptimizer.FindRareCodonClusters(
            SixAtThresholdPlusCommon, CodonOptimizer.EColiK12, 0.15, 7, 4, inclusiveThreshold: false);
        Assert.That(strict, Is.Empty);

        // Codons strictly below the threshold are pauses under both comparisons.
        string agaRun = string.Concat(Enumerable.Repeat("AGA", 7));
        Assert.That(CodonOptimizer.FindRareCodonClusters(agaRun, CodonOptimizer.EColiK12, 0.15, 7, 4, inclusiveThreshold: false),
            Is.EqualTo(CodonOptimizer.FindRareCodonClusters(agaRun, CodonOptimizer.EColiK12)));
    }

    [Test]
    public void FindRareCodonClusters_ExplicitInclusive_EqualsDefault()
    {
        Assert.That(
            CodonOptimizer.FindRareCodonClusters(SixAtThresholdPlusCommon, CodonOptimizer.EColiK12, 0.15, 7, 4, inclusiveThreshold: true),
            Is.EqualTo(CodonOptimizer.FindRareCodonClusters(SixAtThresholdPlusCommon, CodonOptimizer.EColiK12)));
    }

    [Test]
    public void FindRareCodons_DefaultStaysStrict_InclusiveOptionAddsCodonsAtThreshold()
    {
        var strict = CodonOptimizer.FindRareCodons("UCUAGA", CodonOptimizer.EColiK12).ToList();
        Assert.That(strict, Is.EqualTo(new[] { (3, "AGA", "R", 0.04) }));

        var inclusive = CodonOptimizer.FindRareCodons("UCUAGA", CodonOptimizer.EColiK12, 0.15, inclusiveThreshold: true).ToList();
        Assert.That(inclusive, Is.EqualTo(new[] { (0, "UCU", "S", 0.15), (3, "AGA", "R", 0.04) }));

        Assert.That(CodonOptimizer.FindRareCodons("UCUAGA", CodonOptimizer.EColiK12, 0.15, inclusiveThreshold: false),
            Is.EqualTo(strict));
    }

    #endregion

    #region A1-3 — per-thousand presets

    [Test]
    public void Presets_PerThousand_SourcedEntries()
    {
        // SRC RO1 (2026-10-14): was Is.Null (Kazusa 316407 unreachable in review 2026-10); the
        // column is now the published Kazusa 316407 one (exact 64-value check: EColiPerThousand_*).
        Assert.That(CodonOptimizer.EColiK12.PerThousand, Is.Not.Null);

        var yeast = CodonOptimizer.Yeast.PerThousand!;
        Assert.Multiple(() =>
        {
            Assert.That(yeast, Has.Count.EqualTo(64));
            Assert.That(yeast["UUU"], Is.EqualTo(26.1));
            Assert.That(yeast["GAA"], Is.EqualTo(45.6));
            Assert.That(yeast["UAG"], Is.EqualTo(0.5));
            Assert.That(yeast.Values.Sum(), Is.EqualTo(1000.5).Within(1e-9), "published one-decimal column");
        });

        var human = CodonOptimizer.Human.PerThousand!;
        Assert.Multiple(() =>
        {
            Assert.That(human, Has.Count.EqualTo(64));
            Assert.That(human["UUU"], Is.EqualTo(17.56646934028931).Within(1e-12)); // 1000*714298/40662582
            Assert.That(human["CUG"], Is.EqualTo(39.63843220777274).Within(1e-12)); // 1000*1611801/40662582
            Assert.That(human.Values.Sum(), Is.EqualTo(1000.0).Within(1e-9));
        });
    }

    // The per-thousand columns come from the SAME Kazusa entries as the stored fractions: their
    // per-amino-acid ratios round (2 dp) to the preset fractions — Human at all 64 codons; Yeast at
    // every non-stop codon, the stop family within 0.013 (its per-thousand values have one decimal).
    [Test]
    public void Presets_PerThousandRatios_MatchStoredFractions()
    {
        AssertRatiosMatch(CodonOptimizer.Human, stopTolerance: 0.005);
        AssertRatiosMatch(CodonOptimizer.Yeast, stopTolerance: 0.013);
    }

    private static void AssertRatiosMatch(CodonOptimizer.CodonUsageTable table, double stopTolerance)
    {
        var perThousand = table.PerThousand!;
        Assert.Multiple(() =>
        {
            foreach (var family in GeneticCode.Standard.CodonTable.GroupBy(kv => kv.Value))
            {
                double total = family.Sum(kv => perThousand[kv.Key]);
                double tolerance = family.Key == '*' ? stopTolerance : 0.005;
                foreach (var (codon, _) in family)
                    Assert.That(perThousand[codon] / total, Is.EqualTo(table.CodonFrequencies[codon]).Within(tolerance + 1e-12),
                        $"{table.OrganismName} {codon}");
            }
        });
    }

    private const string Gene20 = "AUGCUGAGAAAAGAAGCUCUGAGGUUUCCCGGAAUAUCUCGAUGGAAGACCCUAGUUCAG";

    // Yeast preset default = Clark-lab reference with Kazusa 4932 per thousand (CHARMING: 33.99;
    // 12.69, 17.03, 15.99) — identical to the per-thousand table test of the 2026-09 review.
    [Test]
    public void CalculateMinMaxProfile_YeastPreset_DefaultMatchesClarkLabReference()
    {
        Assert.That(CodonOptimizer.CalculateMinMaxProfile("CUGAGA", CodonOptimizer.Yeast, 2).Select(w => w.PercentMinMax),
            Is.EqualTo(new[] { 33.99209486166008 }).Within(1e-9));
        Assert.That(CodonOptimizer.CalculateMinMaxProfile(Gene20, CodonOptimizer.Yeast).Select(w => w.PercentMinMax),
            Is.EqualTo(new[] { 12.68791340950089, 17.028571428571425, 15.991062879029688 }).Within(1e-9));
    }

    // Human preset default = reference with Kazusa 9606 per thousand (CHARMING: 28.8, 23.26, 16.2;
    // UUA·CGU: −96.4); CodonFrequencies basis gives the equal-weight fraction values.
    [Test]
    public void CalculateMinMaxProfile_HumanPreset_DefaultUsesPerThousand_FractionBasisOnRequest()
    {
        Assert.That(CodonOptimizer.CalculateMinMaxProfile(Gene20, CodonOptimizer.Human).Select(w => w.PercentMinMax),
            Is.EqualTo(new[] { 28.796581754809075, 23.264123403923747, 16.195449171123403 }).Within(1e-9));
        Assert.That(CodonOptimizer.CalculateMinMaxProfile("UUACGU", CodonOptimizer.Human, 2)[0].PercentMinMax,
            Is.EqualTo(-96.40404591634571).Within(1e-9));

        var fractions = CodonOptimizer.CalculateMinMaxProfile(Gene20, CodonOptimizer.Human, 18,
            CodonOptimizer.MinMaxUsageBasis.CodonFrequencies);
        Assert.That(fractions.Select(w => w.PercentMinMax),
            Is.EqualTo(new[] { 19.650655021834133, 13.878080415045439, 13.803548247511891 }).Within(1e-9));
        Assert.That(CodonOptimizer.CalculateMinMaxProfile("UUACGU", CodonOptimizer.Human, 2,
            CodonOptimizer.MinMaxUsageBasis.CodonFrequencies)[0].PercentMinMax,
            Is.EqualTo(-94.54545454545456).Within(1e-9));
    }

    // SRC RO1 (2026-10-14): the E. coli preset now carries the Kazusa 316407 per-thousand column, so
    // the default basis gives the Clark-lab reference value (CHARMING calculateMinMax, run isolated on
    // the Kazusa 316407 /1000 column: 58.1; independent Python port 58.10113519091846). Was
    // 36.470588235294116 under both bases (no per-thousand column in review 2026-10); the
    // CodonFrequencies basis keeps that fraction value (CHARMING on the fractions: 36.47).
    [Test]
    public void CalculateMinMaxProfile_EColiPreset_DefaultPerThousand_FractionBasisUnchanged()
    {
        Assert.That(CodonOptimizer.CalculateMinMaxProfile("CUGAGA", CodonOptimizer.EColiK12, 2)[0].PercentMinMax,
            Is.EqualTo(58.10113519091846).Within(1e-9));
        Assert.That(CodonOptimizer.CalculateMinMaxProfile("CUGAGA", CodonOptimizer.EColiK12, 2,
            CodonOptimizer.MinMaxUsageBasis.CodonFrequencies)[0].PercentMinMax,
            Is.EqualTo(36.470588235294116).Within(1e-10));
    }

    // 20-codon gene, window 18 (3 windows): CHARMING 7.42, 7.37, -6.57 (port 7.417300184918832,
    // 7.371763071733311, -6.565186442534101); fraction basis: CHARMING -1.36, -0.9, -8.18 (port
    // -1.358411703239252, -0.9027081243731403, -8.181517058628705).
    [Test]
    public void CalculateMinMaxProfile_EColiPreset_Gene20_MatchesClarkLabReference()
    {
        Assert.That(CodonOptimizer.CalculateMinMaxProfile(Gene20, CodonOptimizer.EColiK12).Select(w => w.PercentMinMax),
            Is.EqualTo(new[] { 7.417300184918832, 7.371763071733311, -6.565186442534101 }).Within(1e-9));
        Assert.That(CodonOptimizer.CalculateMinMaxProfile(Gene20, CodonOptimizer.EColiK12, 18,
                CodonOptimizer.MinMaxUsageBasis.CodonFrequencies).Select(w => w.PercentMinMax),
            Is.EqualTo(new[] { -1.358411703239252, -0.9027081243731403, -8.181517058628705 }).Within(1e-9));
        Assert.That(CodonOptimizer.CalculateMinMaxProfile("UUACGU", CodonOptimizer.EColiK12, 2)[0].PercentMinMax,
            Is.EqualTo(16.202270381836932).Within(1e-9)); // CHARMING 16.2; fraction basis 31.18
    }

    // The 64 values equal the Kazusa 316407 "/1000" column verbatim (one decimal as published).
    [Test]
    public void EColiPerThousand_EqualsKazusa316407Verbatim()
    {
        var expected = new Dictionary<string, double>
        {
            ["UUU"] = 22.2, ["UUC"] = 16.5, ["UUA"] = 13.8, ["UUG"] = 13.6, ["UCU"] = 8.4, ["UCC"] = 8.6, ["UCA"] = 7.0, ["UCG"] = 8.9,
            ["UAU"] = 16.1, ["UAC"] = 12.2, ["UAA"] = 2.0, ["UAG"] = 0.2, ["UGU"] = 5.1, ["UGC"] = 6.4, ["UGA"] = 0.9, ["UGG"] = 15.2,
            ["CUU"] = 11.0, ["CUC"] = 11.1, ["CUA"] = 3.8, ["CUG"] = 53.1, ["CCU"] = 7.0, ["CCC"] = 5.5, ["CCA"] = 8.4, ["CCG"] = 23.4,
            ["CAU"] = 13.0, ["CAC"] = 9.8, ["CAA"] = 15.4, ["CAG"] = 29.0, ["CGU"] = 21.0, ["CGC"] = 22.3, ["CGA"] = 3.5, ["CGG"] = 5.4,
            ["AUU"] = 30.4, ["AUC"] = 25.2, ["AUA"] = 4.2, ["AUG"] = 27.8, ["ACU"] = 8.8, ["ACC"] = 23.5, ["ACA"] = 6.9, ["ACG"] = 14.4,
            ["AAU"] = 17.6, ["AAC"] = 21.6, ["AAA"] = 33.6, ["AAG"] = 10.3, ["AGU"] = 8.7, ["AGC"] = 16.1, ["AGA"] = 2.0, ["AGG"] = 1.1,
            ["GUU"] = 18.2, ["GUC"] = 15.3, ["GUA"] = 10.9, ["GUG"] = 26.3, ["GCU"] = 15.2, ["GCC"] = 25.7, ["GCA"] = 20.1, ["GCG"] = 33.9,
            ["GAU"] = 32.2, ["GAC"] = 19.1, ["GAA"] = 39.7, ["GAG"] = 18.0, ["GGU"] = 24.7, ["GGC"] = 29.8, ["GGA"] = 7.9, ["GGG"] = 11.0,
        };
        var actual = CodonOptimizer.EColiK12.PerThousand!;
        Assert.Multiple(() =>
        {
            Assert.That(actual, Has.Count.EqualTo(64));
            Assert.That(actual.Keys, Is.EquivalentTo(CodonOptimizer.EColiK12.CodonFrequencies.Keys));
            foreach (var (codon, value) in expected)
                Assert.That(actual[codon], Is.EqualTo(value), codon);
            Assert.That(actual.Values.Sum(), Is.EqualTo(1000.0).Within(1e-9), "published one-decimal column");
        });
    }

    // Same Kazusa entry as the stored fractions: the per-amino-acid ratios of the one-decimal
    // values are within 0.0055 of the fractions at all 64 codons (max 0.00548 at UAG, 0.00538 at
    // AAA/AAG — the 0.1 rounding of the /1000 column; the raw counts give ≤ 0.0049).
    [Test]
    public void EColiPerThousandRatios_MatchStoredFractions()
    {
        var perThousand = CodonOptimizer.EColiK12.PerThousand!;
        Assert.Multiple(() =>
        {
            foreach (var family in GeneticCode.Standard.CodonTable.GroupBy(kv => kv.Value))
            {
                double total = family.Sum(kv => perThousand[kv.Key]);
                foreach (var (codon, _) in family)
                    Assert.That(perThousand[codon] / total, Is.EqualTo(CodonOptimizer.EColiK12.CodonFrequencies[codon]).Within(0.0055),
                        codon);
            }
        });
    }

    // A caller table can carry its own per-thousand column; the default overload then uses it.
    [Test]
    public void CalculateMinMaxProfile_CustomTableWithPerThousand_UsesIt()
    {
        var custom = CodonOptimizer.CreateCodonUsageTable("yeast fractions", CodonOptimizer.Yeast.CodonFrequencies)
            with { PerThousand = CodonOptimizer.Yeast.PerThousand };
        Assert.That(CodonOptimizer.CalculateMinMaxProfile("CUGAGA", custom, 2)[0].PercentMinMax,
            Is.EqualTo(33.99209486166008).Within(1e-9));
    }

    [Test]
    public void CalculateMinMaxProfile_UndefinedBasis_Throws()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() =>
            CodonOptimizer.CalculateMinMaxProfile("CUGAGA", CodonOptimizer.Yeast, 2, (CodonOptimizer.MinMaxUsageBasis)7));
    }

    #endregion
}
