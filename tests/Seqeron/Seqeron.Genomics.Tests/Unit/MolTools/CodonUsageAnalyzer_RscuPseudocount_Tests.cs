// CODON-RSCU-001 (audit A2-9) — optional RSCU pseudocount:
//   CalculateRscu(counts | string | DnaSequence, GeneticCode, double pseudocount).
// Evidence: docs/Evidence/CODON-RSCU-001-Evidence.md; TestSpec: tests/TestSpecs/CODON-RSCU-001.md (M17–M21, S8)
// References:
//   codonbias 0.5.0 (PyPI codon-bias, scores.py RelativeSynonymousCodonUsage._weights_and_counts):
//     counts_pc = counts + pseudocount; P = counts_pc / Σ_aa counts_pc; RSCU = P / (1/n).
//   cubar (github.com/mt1022/cubar, R/codon_weight.R est_rscu): prop = (cts + pseudo_cnt)/sum(cts + pseudo_cnt),
//     rscu = prop/mean(prop) — the same formula at level = "amino_acid", incl_stop = TRUE.
//   Locked TSV: TestData/Codon/codonbias_rscu_pseudocount.tsv =
//     RelativeSynonymousCodonUsage(genetic_code=t, ignore_stop=False, pseudocount=c).get_weights(gene),
//     7 genes × tables 1/2 × c = 0.5/1/2.5 × 64 codons (2688 values).

using System.Globalization;
using System.Reflection;

namespace Seqeron.Genomics.Tests.Unit.MolTools;

[TestFixture]
public class CodonUsageAnalyzer_RscuPseudocount_Tests
{
    private const double Tol = 1e-12;
    private const string Resource = "Seqeron.Genomics.Tests.TestData.Codon.codonbias_rscu_pseudocount.tsv";

    private sealed record Row(int Table, double Pseudocount, string Gene, string Codon, double Rscu);

    private static (Dictionary<string, string> Genes, List<Row> Rows) Load()
    {
        var asm = Assembly.GetExecutingAssembly();
        using var stream = asm.GetManifestResourceStream(Resource)
            ?? throw new InvalidOperationException($"Embedded oracle '{Resource}' not found.");
        using var reader = new StreamReader(stream);
        var genes = new Dictionary<string, string>();
        var rows = new List<Row>();
        string? line;
        while ((line = reader.ReadLine()) is not null)
        {
            if (line.Length == 0 || line[0] == '#') continue;
            var f = line.Split('\t');
            if (f[0] == "G")
                genes[f[1]] = f[2];
            else
                rows.Add(new Row(int.Parse(f[1], CultureInfo.InvariantCulture),
                    double.Parse(f[2], CultureInfo.InvariantCulture), f[3], f[4],
                    double.Parse(f[5], CultureInfo.InvariantCulture)));
        }
        return (genes, rows);
    }

    // M17 — codonbias 0.5.0 lock: 7 genes × NCBI tables 1/2 × pseudocount 0.5/1/2.5 × 64 codons.
    [Test]
    public void CalculateRscu_Pseudocount_MatchesCodonbias()
    {
        var (genes, rows) = Load();
        Assert.That(genes, Has.Count.EqualTo(7), "7 locked genes");
        Assert.That(rows, Has.Count.EqualTo(7 * 2 * 3 * 64), "2688 locked RSCU values");

        var cache = new Dictionary<(int, double, string), Dictionary<string, double>>();
        Assert.Multiple(() =>
        {
            foreach (var r in rows)
            {
                var key = (r.Table, r.Pseudocount, r.Gene);
                if (!cache.TryGetValue(key, out var rscu))
                    cache[key] = rscu = CodonUsageAnalyzer.CalculateRscu(
                        genes[r.Gene], GeneticCode.GetByTableNumber(r.Table), r.Pseudocount);
                Assert.That(rscu[r.Codon], Is.EqualTo(r.Rscu).Within(Tol),
                    $"table {r.Table}, pseudocount {r.Pseudocount}, {r.Gene}, {r.Codon}");
            }
        });
    }

    // M18 — hand derivation: Phe TTT×2 TTC×1, c = 1 → TTT = 2·3/(3+2) = 1.2, TTC = 2·2/5 = 0.8;
    // absent families are uniform: Met (n=1) = 1, Leu (n=6) = 6·1/6 = 1, stops (n=3) = 1.
    [Test]
    public void CalculateRscu_PseudocountOne_HandDerivedValues()
    {
        var rscu = CodonUsageAnalyzer.CalculateRscu("TTTTTTTTC", GeneticCode.Standard, 1.0);

        Assert.Multiple(() =>
        {
            Assert.That(rscu, Has.Count.EqualTo(64), "all 64 codons reported");
            Assert.That(rscu["TTT"], Is.EqualTo(1.2).Within(Tol), "2·(2+1)/(3+2·1)");
            Assert.That(rscu["TTC"], Is.EqualTo(0.8).Within(Tol), "2·(1+1)/(3+2·1)");
            Assert.That(rscu["ATG"], Is.EqualTo(1.0).Within(Tol), "absent single-codon family → 1");
            Assert.That(rscu["CTG"], Is.EqualTo(1.0).Within(Tol), "absent Leu family → uniform 1");
            Assert.That(rscu["TGA"], Is.EqualTo(1.0).Within(Tol), "absent stop family → uniform 1");
        });
    }

    // M19 — table 2 families with c = 0.5: Met {ATA×1, ATG×2}: ATA = 2·1.5/4 = 0.75, ATG = 2·2.5/4 = 1.25;
    // stops {TAA×1, AGA×1, TAG, AGG}: TAA = AGA = 4·1.5/4 = 1.5, TAG = AGG = 4·0.5/4 = 0.5.
    [Test]
    public void CalculateRscu_Pseudocount_Table2Families()
    {
        var rscu = CodonUsageAnalyzer.CalculateRscu(
            "ATAATGATGTGATGGTGGAGATAA", GeneticCode.GetByTableNumber(2), 0.5);

        Assert.Multiple(() =>
        {
            Assert.That(rscu["ATA"], Is.EqualTo(0.75).Within(Tol), "Met(ATA,ATG): 2·1.5/(3+1)");
            Assert.That(rscu["ATG"], Is.EqualTo(1.25).Within(Tol), "Met(ATA,ATG): 2·2.5/(3+1)");
            Assert.That(rscu["TAA"], Is.EqualTo(1.5).Within(Tol), "stops(TAA,TAG,AGA,AGG): 4·1.5/(2+2)");
            Assert.That(rscu["AGA"], Is.EqualTo(1.5).Within(Tol), "AGA is a stop in table 2");
            Assert.That(rscu["TAG"], Is.EqualTo(0.5).Within(Tol), "4·0.5/(2+2)");
            Assert.That(rscu["AGG"], Is.EqualTo(0.5).Within(Tol), "4·0.5/(2+2)");
        });
    }

    // M20 — pseudocount 0 is the CodonW default, bit-identical to the existing overloads (absent family 0).
    [TestCase("TTTTTTTTC", 1)]
    [TestCase("ATAATGATGTGATGGTGGAGATAA", 2)]
    [TestCase("AGAAGGTAAATAATGTGATGG", 1)]
    [TestCase("auguuuNNNcugcugcuaUAA", 11)]
    public void CalculateRscu_PseudocountZero_IsBitIdenticalToDefault(string seq, int table)
    {
        var code = GeneticCode.GetByTableNumber(table);
        var expected = CodonUsageAnalyzer.CalculateRscu(seq, code);

        Assert.Multiple(() =>
        {
            Assert.That(CodonUsageAnalyzer.CalculateRscu(seq, code, 0.0), Is.EqualTo(expected), "string overload");
            Assert.That(CodonUsageAnalyzer.CalculateRscu(CodonUsageAnalyzer.CountCodons(seq), code, 0.0),
                Is.EqualTo(expected), "counts overload");
            Assert.That(expected["GGG"], Is.EqualTo(0.0), "absent family stays 0 without a pseudocount");
        });
    }

    // M21 — the three entry points agree for c > 0; family sums stay n (INV-4) with the pseudocount.
    [Test]
    public void CalculateRscu_Pseudocount_OverloadsAgree_AndFamilySumIsN()
    {
        const string seq = "ATGGCTGCAGCGAAAAAGTTTGGCGGTGGGCGTAGACGAAGCTCATCTTGA";
        var code = GeneticCode.Standard;
        var fromString = CodonUsageAnalyzer.CalculateRscu(seq, code, 1.0);
        var fromDna = CodonUsageAnalyzer.CalculateRscu(new DnaSequence(seq), code, 1.0);
        var fromCounts = CodonUsageAnalyzer.CalculateRscu(CodonUsageAnalyzer.CountCodons(seq), code, 1.0);

        Assert.Multiple(() =>
        {
            Assert.That(fromDna, Is.EqualTo(fromString), "DnaSequence overload");
            Assert.That(fromCounts, Is.EqualTo(fromString), "counts overload");
            foreach (var family in code.CodonTable.GroupBy(kv => kv.Value, kv => kv.Key.Replace('U', 'T')))
                Assert.That(family.Sum(c => fromString[c]), Is.EqualTo(family.Count()).Within(1e-12),
                    $"family '{family.Key}' sums to n");
        });
    }

    // S8 — invalid pseudocount and null arguments.
    [TestCase(-0.5)]
    [TestCase(double.NaN)]
    [TestCase(double.PositiveInfinity)]
    public void CalculateRscu_InvalidPseudocount_Throws(double pc)
    {
        Assert.Multiple(() =>
        {
            Assert.That(() => CodonUsageAnalyzer.CalculateRscu("TTT", GeneticCode.Standard, pc),
                NUnit.Framework.Throws.TypeOf<ArgumentOutOfRangeException>());
            Assert.That(() => CodonUsageAnalyzer.CalculateRscu(new DnaSequence("TTT"), GeneticCode.Standard, pc),
                NUnit.Framework.Throws.TypeOf<ArgumentOutOfRangeException>());
            Assert.That(() => CodonUsageAnalyzer.CalculateRscu(new Dictionary<string, int>(), GeneticCode.Standard, pc),
                NUnit.Framework.Throws.TypeOf<ArgumentOutOfRangeException>());
        });
    }

    [Test]
    public void CalculateRscu_Pseudocount_NullArguments_Throw()
    {
        Assert.Multiple(() =>
        {
            Assert.That(() => CodonUsageAnalyzer.CalculateRscu((DnaSequence)null!, GeneticCode.Standard, 1.0),
                NUnit.Framework.Throws.ArgumentNullException);
            Assert.That(() => CodonUsageAnalyzer.CalculateRscu("TTT", null!, 1.0), NUnit.Framework.Throws.ArgumentNullException);
            Assert.That(() => CodonUsageAnalyzer.CalculateRscu((IReadOnlyDictionary<string, int>)null!, GeneticCode.Standard, 1.0),
                NUnit.Framework.Throws.ArgumentNullException);
            Assert.That(CodonUsageAnalyzer.CalculateRscu("", GeneticCode.Standard, 1.0), Is.Empty, "empty → empty");
        });
    }
}
