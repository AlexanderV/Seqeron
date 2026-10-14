// CODON-STATS-001 (audit A3-13, FIN-B02 WP23) — CodonW 1.4.4 Fop (Ikemura 1981) and CBI (Bennetzen & Hall 1982):
//   CalculateFop/TryCalculateFop, CalculateCbi/TryCalculateCbi, CodonWOptimalCodonSets,
//   GetStatistics(sequence, code, gcConvention, fopCbiOptimalCodons) → CodonUsageStatistics.Fop/Cbi.
// Evidence: docs/Evidence/CODON-STATS-001-Evidence.md; TestSpec: tests/TestSpecs/CODON-STATS-001.md (R17–R24)
// Reference: CodonW 1.4.4 (Ubuntu codonw_1.4.4.orig.tar.gz, compiled locally), codon_us.c fop_out / cbi_out:
//   Fop = opt/total (modified: (opt − rare)/total) and CBI = (opt − Σ naa/da)/(tot − Σ naa/da), both over the
//   amino acids with > 1 codon (not stop) having ≥ 1 optimal codon; 0.000 printed when the denominator is 0.
//   Locked TSV: TestData/Codon/codonw_fop_cbi.tsv — 70 genes × (8 built-in -f_type sets + 3 user sets) ×
//   -code 0..7 (NCBI 1,2,3,4,5,6,10,9) + modified Fop runs: 7 280 rows. CodonW prints %5.3f, so each value is
//   locked to ±0.0005 (half a unit of the last printed digit).

using System.Globalization;
using System.Reflection;

namespace Seqeron.Genomics.Tests.Unit.MolTools;

[TestFixture]
public class CodonUsageAnalyzer_FopCbi_Tests
{
    private const string Resource = "Seqeron.Genomics.Tests.TestData.Codon.codonw_fop_cbi.tsv";
    private const double PrintTol = 0.0005 + 1e-9;
    private const double Tol = 1e-12;

    private sealed record Row(string Set, int Table, string Gene, double Fop, double? Cbi, bool Modified);

    private sealed record Oracle(
        Dictionary<string, string> Genes,
        Dictionary<string, (IReadOnlySet<string> Optimal, IReadOnlySet<string> Rare)> UserSets,
        List<Row> Rows);

    private static readonly Lazy<Oracle> Data = new(Load);

    private static Oracle Load()
    {
        var asm = Assembly.GetExecutingAssembly();
        using var stream = asm.GetManifestResourceStream(Resource)
            ?? throw new InvalidOperationException($"Embedded oracle '{Resource}' not found.");
        using var reader = new StreamReader(stream);
        var genes = new Dictionary<string, string>();
        var users = new Dictionary<string, (IReadOnlySet<string>, IReadOnlySet<string>)>();
        var rows = new List<Row>();
        static IReadOnlySet<string> Set(string csv) =>
            csv.Split(',', StringSplitOptions.RemoveEmptyEntries).ToHashSet();
        string? line;
        while ((line = reader.ReadLine()) is not null)
        {
            if (line.Length == 0 || line[0] == '#') continue;
            var f = line.Split('\t');
            switch (f[0])
            {
                case "G": genes[f[1]] = f[2]; break;
                case "U": users[f[1]] = (Set(f[2]), Set(f[3])); break;
                default:
                    rows.Add(new Row(f[1], int.Parse(f[2], CultureInfo.InvariantCulture), f[3],
                        double.Parse(f[4], CultureInfo.InvariantCulture),
                        f[5] == "NA" ? null : double.Parse(f[5], CultureInfo.InvariantCulture),
                        f[6] == "1"));
                    break;
            }
        }
        return new Oracle(genes, users, rows);
    }

    private static (IReadOnlySet<string> Optimal, IReadOnlySet<string> Rare) ResolveSet(Oracle o, string name)
    {
        if (o.UserSets.TryGetValue(name, out var user)) return user;
        Assert.That(CodonWOptimalCodonSets.TryGetByName(name, out var set), Is.True, name);
        return (set!.OptimalCodons, set.RareCodons);
    }

    #region CodonW binary lock (R17, R22)

    [Test]
    public void FopCbi_MatchCodonWBinary_AllGenesSetsAndCodes()
    {
        var o = Data.Value;
        Assert.That(o.Genes, Has.Count.EqualTo(70));
        Assert.That(o.Rows, Has.Count.EqualTo(7280));

        int checkedFop = 0, checkedCbi = 0;
        Assert.Multiple(() =>
        {
            foreach (var row in o.Rows)
            {
                var (optimal, rare) = ResolveSet(o, row.Set);
                var code = GeneticCode.GetByTableNumber(row.Table);
                string seq = o.Genes[row.Gene];
                string id = $"{row.Set}/NCBI {row.Table}/{row.Gene}{(row.Modified ? "/modified" : "")}";

                bool fopOk = CodonUsageAnalyzer.TryCalculateFop(seq, out double fop, optimal, code,
                    row.Modified ? rare : null);
                Assert.That(fop, Is.EqualTo(row.Fop).Within(PrintTol), $"Fop {id}");
                if (!fopOk) Assert.That(row.Fop, Is.EqualTo(0.0), $"Fop not calculable ⇒ CodonW 0.000 {id}");
                checkedFop++;

                if (row.Cbi is double expectedCbi)
                {
                    bool cbiOk = CodonUsageAnalyzer.TryCalculateCbi(seq, out double cbi, optimal, code);
                    Assert.That(cbi, Is.EqualTo(expectedCbi).Within(PrintTol), $"CBI {id}");
                    if (!cbiOk) Assert.That(expectedCbi, Is.EqualTo(0.0), $"CBI not calculable ⇒ CodonW 0.000 {id}");
                    checkedCbi++;
                }
            }
        });
        Assert.That(checkedFop, Is.EqualTo(7280));
        Assert.That(checkedCbi, Is.EqualTo(6160));
    }

    [Test]
    public void ModifiedFop_BacillusSubtilis_DiffersFromOriginalAndCanBeNegative()
    {
        // CodonW's modified index (opt − rare)/total is not clamped (codon_us.c fop_out), so it goes
        // negative for genes rich in B. subtilis rare codons; the oracle contains such values.
        var o = Data.Value;
        var negatives = o.Rows.Where(r => r.Modified && r.Fop < 0).ToList();
        Assert.That(negatives, Is.Not.Empty);
        var r0 = negatives[0];
        var (optimal, rare) = ResolveSet(o, r0.Set);
        double modified = CodonUsageAnalyzer.CalculateFop(o.Genes[r0.Gene], optimal,
            GeneticCode.GetByTableNumber(r0.Table), rare);
        double original = CodonUsageAnalyzer.CalculateFop(o.Genes[r0.Gene], optimal,
            GeneticCode.GetByTableNumber(r0.Table));
        Assert.That(modified, Is.LessThan(0));
        Assert.That(original, Is.GreaterThanOrEqualTo(0));
    }

    #endregion

    #region Hand-computed values (R18–R19)

    // ATG CTG×4 AAA AAG TTT TTC GGC GGT GGA TAA under the default E. coli set:
    // optimal CTG×4, AAA, TTC, GGC, GGT = 8 of 11 scored codons (Met and the stop excluded).
    // N_ran = 4·1/6 + 2·1/2 + 2·1/2 + 3·2/4 = 25/6 ⇒ CBI = (8 − 25/6)/(11 − 25/6) = 23/41.
    // CodonW 1.4.4: Fop 0.727, CBI 0.561.
    private const string HandGene = "ATGCTGCTGCTGCTGAAAAAGTTTTTCGGCGGTGGATAA";

    [Test]
    public void Fop_HandGene_DefaultEColi_EightOfEleven()
    {
        Assert.That(CodonUsageAnalyzer.TryCalculateFop(HandGene, out double fop), Is.True);
        Assert.That(fop, Is.EqualTo(8.0 / 11.0).Within(Tol));
    }

    [Test]
    public void Cbi_HandGene_DefaultEColi_TwentyThreeOverFortyOne()
    {
        Assert.That(CodonUsageAnalyzer.TryCalculateCbi(HandGene, out double cbi), Is.True);
        Assert.That(cbi, Is.EqualTo(23.0 / 41.0).Within(Tol));
    }

    [Test]
    public void FopCbi_HandGene_VertebrateMito_FamiliesFollowCode()
    {
        // NCBI 2 (CodonW -code 1): Leu loses nothing here, but CodonW prints Fop 0.273 / CBI −0.055 with
        // the B. subtilis set (oracle row); under the E. coli set families are unchanged ⇒ same as table 1.
        var code = GeneticCode.GetByTableNumber(2);
        Assert.That(CodonUsageAnalyzer.CalculateFop(HandGene, null, code), Is.EqualTo(8.0 / 11.0).Within(Tol));
        var bsu = CodonWOptimalCodonSets.BacillusSubtilis.OptimalCodons;
        Assert.That(CodonUsageAnalyzer.CalculateFop(HandGene, bsu, code), Is.EqualTo(3.0 / 11.0).Within(Tol));
        Assert.That(CodonUsageAnalyzer.CalculateCbi(HandGene, bsu, code), Is.EqualTo(-0.055).Within(PrintTol));
    }

    #endregion

    #region Not calculable (R20–R21)

    [TestCase("ATGTGGATGTGGTAA", TestName = "FopCbi_MetTrpStopOnly_NotCalculable")]
    [TestCase("TAATAGTGA", TestName = "FopCbi_StopsOnly_NotCalculable")]
    [TestCase("ATGNNN", TestName = "FopCbi_NoScoredCodon_NotCalculable")]
    [TestCase("", TestName = "FopCbi_Empty_NotCalculable")]
    public void FopCbi_NoScoredCodon_TryFalseAndCodonWZero(string seq)
    {
        Assert.Multiple(() =>
        {
            Assert.That(CodonUsageAnalyzer.TryCalculateFop(seq, out double fop), Is.False);
            Assert.That(fop, Is.EqualTo(0.0));
            Assert.That(CodonUsageAnalyzer.TryCalculateCbi(seq, out double cbi), Is.False);
            Assert.That(cbi, Is.EqualTo(0.0));
            Assert.That(CodonUsageAnalyzer.CalculateFop(seq), Is.EqualTo(0.0));
            Assert.That(CodonUsageAnalyzer.CalculateCbi(seq), Is.EqualTo(0.0));
        });
    }

    [Test]
    public void Cbi_AllScoredCodonsOptimal_NotCalculable_FopOne()
    {
        // Phe = {TTT, TTC} both optimal: N_opt = N_tot = N_ran ⇒ 0/0, CodonW prints CBI 0.000 (oracle
        // user3 under NCBI 2: Met {ATA, ATG} and Trp {TGA, TGG} all optimal).
        var set = new HashSet<string> { "TTT", "TTC" };
        Assert.That(CodonUsageAnalyzer.TryCalculateCbi("TTTTTCTTT", out double cbi, set), Is.False);
        Assert.That(cbi, Is.EqualTo(0.0));
        Assert.That(CodonUsageAnalyzer.CalculateFop("TTTTTCTTT", set), Is.EqualTo(1.0));
    }

    #endregion

    #region Codon sets and validation (R23)

    [Test]
    public void BuiltInSets_TranscribeCodonW()
    {
        int[] optimalCounts = [24, 19, 15, 19, 22, 20, 20, 22];
        var all = CodonWOptimalCodonSets.All;
        Assert.That(all, Has.Count.EqualTo(8));
        for (int i = 0; i < all.Count; i++)
        {
            Assert.That(all[i].CodonWFType, Is.EqualTo(i));
            Assert.That(all[i].OptimalCodons, Has.Count.EqualTo(optimalCounts[i]), all[i].Organism);
            Assert.That(all[i].RareCodons, Has.Count.EqualTo(i == 1 ? 10 : 0), all[i].Organism);
            Assert.That(all[i].OptimalCodons.Overlaps(all[i].RareCodons), Is.False);
        }
        Assert.That(CodonWOptimalCodonSets.EscherichiaColi.OptimalCodons, Does.Contain("CTG").And.Not.Contain("CTA"));
    }

    [TestCase("EscherichiaColi", 0)]
    [TestCase("Escherichia coli", 0)]
    [TestCase("e. coli", 0)]
    [TestCase("B_subtilis", 1)]
    [TestCase("saccharomyces cerevisiae", 4)]
    [TestCase("NeurosporaCrassa", 7)]
    public void TryGetByName_Resolves(string name, int fType)
    {
        Assert.That(CodonWOptimalCodonSets.TryGetByName(name, out var set), Is.True);
        Assert.That(set!.CodonWFType, Is.EqualTo(fType));
    }

    [TestCase("Homo sapiens")]
    [TestCase("")]
    [TestCase(null)]
    public void TryGetByName_Unknown_False(string? name) =>
        Assert.That(CodonWOptimalCodonSets.TryGetByName(name, out _), Is.False);

    [Test]
    public void UserSet_RnaLowerCase_EqualsDna()
    {
        var dna = new HashSet<string> { "CTG", "AAA", "GGC" };
        var rna = new HashSet<string> { "cug", "aaa", "Ggc" };
        Assert.That(CodonUsageAnalyzer.CalculateFop(HandGene, rna), Is.EqualTo(CodonUsageAnalyzer.CalculateFop(HandGene, dna)));
        Assert.That(CodonUsageAnalyzer.CalculateCbi(HandGene, rna), Is.EqualTo(CodonUsageAnalyzer.CalculateCbi(HandGene, dna)));
    }

    [Test]
    public void InvalidCodonSets_Throw()
    {
        Assert.Throws<ArgumentException>(() => CodonUsageAnalyzer.CalculateFop(HandGene, new HashSet<string> { "CTGA" }));
        Assert.Throws<ArgumentException>(() => CodonUsageAnalyzer.CalculateCbi(HandGene, new HashSet<string> { "CTN" }));
        Assert.Throws<ArgumentException>(() => CodonUsageAnalyzer.CalculateFop(HandGene,
            new HashSet<string> { "CTG" }, null, new HashSet<string> { "CUG" }));
    }

    #endregion

    #region Overloads and statistics (R24)

    [Test]
    public void DnaSequenceOverloads_EqualString()
    {
        var dna = new DnaSequence(HandGene);
        Assert.That(CodonUsageAnalyzer.CalculateFop(dna), Is.EqualTo(CodonUsageAnalyzer.CalculateFop(HandGene)));
        Assert.That(CodonUsageAnalyzer.CalculateCbi(dna), Is.EqualTo(CodonUsageAnalyzer.CalculateCbi(HandGene)));
        Assert.That(CodonUsageAnalyzer.TryCalculateFop(dna, out double f), Is.True);
        Assert.That(CodonUsageAnalyzer.TryCalculateCbi(dna, out double c), Is.True);
        Assert.That((f, c), Is.EqualTo((8.0 / 11.0, 23.0 / 41.0)));
    }

    [Test]
    public void GetStatistics_FopCbiOptIn()
    {
        var plain = CodonUsageAnalyzer.GetStatistics(HandGene, GeneticCode.Standard, CodonGcConvention.EmbossCusp);
        var with = CodonUsageAnalyzer.GetStatistics(HandGene, GeneticCode.Standard, CodonGcConvention.EmbossCusp, null);
        Assert.Multiple(() =>
        {
            Assert.That(plain.Fop, Is.Null);
            Assert.That(plain.Cbi, Is.Null);
            Assert.That(with.Fop, Is.EqualTo(8.0 / 11.0).Within(Tol));
            Assert.That(with.Cbi, Is.EqualTo(23.0 / 41.0).Within(Tol));
            Assert.That(with.Enc, Is.EqualTo(plain.Enc));
            Assert.That(with.Gc3s, Is.EqualTo(plain.Gc3s));
            Assert.That(with.TotalCodons, Is.EqualTo(plain.TotalCodons));
            var none = CodonUsageAnalyzer.GetStatistics("ATGTGG", GeneticCode.Standard, CodonGcConvention.CodonW, null);
            Assert.That(none.Fop, Is.Null, "not calculable ⇒ null");
            Assert.That(none.Cbi, Is.Null);
        });
    }

    #endregion
}
