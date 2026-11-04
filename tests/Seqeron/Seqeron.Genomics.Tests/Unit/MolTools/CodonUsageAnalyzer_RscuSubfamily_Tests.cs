// CODON-RSCU-001 (audit A3-4) — cubar sub-family RSCU:
//   CalculateRscu(counts | string | DnaSequence, GeneticCode, double pseudocount, RscuFamilyLevel).
// Evidence: docs/Evidence/CODON-RSCU-001-Evidence.md; TestSpec: tests/TestSpecs/CODON-RSCU-001.md (M22–M26, S9)
// References:
//   cubar (github.com/mt1022/cubar @350b295) R/codon_info.R get_codon_table:
//     subfam = paste(amino_acid, substr(codon, 1, 2), sep = '_');
//   R/codon_weight.R est_rscu(level = 'subfam'): prop = (cts + pseudo_cnt)/sum(cts + pseudo_cnt) by subfam,
//     rscu = prop/mean(prop) by subfam.
//   Locked TSV: TestData/Codon/cubar_rscu_subfam.tsv = independent Python port of est_rscu(level = "subfam",
//     incl_stop = TRUE), 7 genes × tables 1/2 × pseudo_cnt 0/0.5/1 × 64 codons (2688 values; NaN = 0/0).

using System.Globalization;
using System.Reflection;

namespace Seqeron.Genomics.Tests.Unit.MolTools;

[TestFixture]
public class CodonUsageAnalyzer_RscuSubfamily_Tests
{
    private const double Tol = 1e-12;
    private const string Resource = "Seqeron.Genomics.Tests.TestData.Codon.cubar_rscu_subfam.tsv";

    private sealed record Row(int Table, double Pseudocount, string Gene, string Codon, string Subfam, double Rscu);

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
                    double.Parse(f[2], CultureInfo.InvariantCulture), f[3], f[4], f[5],
                    f[6] == "nan" ? double.NaN : double.Parse(f[6], CultureInfo.InvariantCulture)));
        }
        return (genes, rows);
    }

    // M22 — cubar lock: 7 genes × NCBI tables 1/2 × pseudo_cnt 0/0.5/1 × 64 codons. cubar's NaN (0/0, an absent
    // sub-family with pseudo_cnt = 0) is this class's CodonW convention 0.
    [Test]
    public void CalculateRscu_Subfamily_MatchesCubar()
    {
        var (genes, rows) = Load();
        Assert.That(genes, Has.Count.EqualTo(7));
        Assert.That(rows, Has.Count.EqualTo(7 * 2 * 3 * 64));

        var cache = new Dictionary<(int, double, string), Dictionary<string, double>>();
        Assert.Multiple(() =>
        {
            foreach (var r in rows)
            {
                var key = (r.Table, r.Pseudocount, r.Gene);
                if (!cache.TryGetValue(key, out var rscu))
                    cache[key] = rscu = CodonUsageAnalyzer.CalculateRscu(genes[r.Gene],
                        GeneticCode.GetByTableNumber(r.Table), r.Pseudocount, RscuFamilyLevel.Subfamily);
                double expected = double.IsNaN(r.Rscu) ? 0.0 : r.Rscu;
                Assert.That(rscu[r.Codon], Is.EqualTo(expected).Within(Tol),
                    $"table {r.Table}, pseudo_cnt {r.Pseudocount}, {r.Gene}, {r.Codon} ({r.Subfam})");
            }
        });
    }

    // M23 — hand derivation, table 1, c = 0: CTGCTGCTGCTATTA → Leu_CT {CTT, CTC, CTA, CTG} = 0, 0, 1, 3 → 4·x/4:
    // CTA 1, CTG 3; Leu_TT {TTA 1, TTG 0} → 2·1/1 = 2, 0. Amino-acid level (6-fold): CTG 6·3/5 = 3.6, TTA 1.2.
    [Test]
    public void CalculateRscu_Subfamily_LeucineSplitsByFirstTwoBases()
    {
        const string gene = "CTGCTGCTGCTATTA";
        var sub = CodonUsageAnalyzer.CalculateRscu(gene, GeneticCode.Standard, 0, RscuFamilyLevel.Subfamily);
        var aa = CodonUsageAnalyzer.CalculateRscu(gene, GeneticCode.Standard, 0, RscuFamilyLevel.AminoAcid);

        Assert.Multiple(() =>
        {
            Assert.That(sub["CTG"], Is.EqualTo(3.0).Within(Tol));
            Assert.That(sub["CTA"], Is.EqualTo(1.0).Within(Tol));
            Assert.That(sub["CTT"], Is.EqualTo(0.0));
            Assert.That(sub["TTA"], Is.EqualTo(2.0).Within(Tol));
            Assert.That(sub["TTG"], Is.EqualTo(0.0));
            Assert.That(sub["AGT"], Is.EqualTo(0.0), "absent Ser_AG, c = 0 → 0 (cubar NaN)");
            Assert.That(aa["CTG"], Is.EqualTo(3.6).Within(Tol));
            Assert.That(aa["TTA"], Is.EqualTo(1.2).Within(Tol));
        });
    }

    // M24 — sub-family sizes: Σ RSCU over each sub-family equals its size; table 1 Ser_TC 4 / Ser_AG 2,
    // Arg_CG 4 / Arg_AG 2, *_TA 2 / *_TG 1, Ile_AT 3 apart from Met_AT 1. Table 2: Arg_AG does not exist
    // (AGA/AGG are stops → *_AG with *_TA {TAA, TAG}).
    [Test]
    public void CalculateRscu_Subfamily_GroupSumsEqualSubfamilySizes()
    {
        var t1 = CodonUsageAnalyzer.CalculateRscu("TCAAGTCGAAGATAATGAATAATG", GeneticCode.Standard, 1.0,
            RscuFamilyLevel.Subfamily);
        var t2 = CodonUsageAnalyzer.CalculateRscu("AGATAAAGGTAG", GeneticCode.GetByTableNumber(2), 0,
            RscuFamilyLevel.Subfamily);
        double Sum(Dictionary<string, double> r, params string[] c) => c.Sum(x => r[x]);

        Assert.Multiple(() =>
        {
            Assert.That(Sum(t1, "TCT", "TCC", "TCA", "TCG"), Is.EqualTo(4.0).Within(Tol));
            Assert.That(Sum(t1, "AGT", "AGC"), Is.EqualTo(2.0).Within(Tol));
            Assert.That(Sum(t1, "CGT", "CGC", "CGA", "CGG"), Is.EqualTo(4.0).Within(Tol));
            Assert.That(Sum(t1, "AGA", "AGG"), Is.EqualTo(2.0).Within(Tol));
            Assert.That(Sum(t1, "TAA", "TAG"), Is.EqualTo(2.0).Within(Tol));
            Assert.That(t1["TGA"], Is.EqualTo(1.0).Within(Tol));
            Assert.That(Sum(t1, "ATT", "ATC", "ATA"), Is.EqualTo(3.0).Within(Tol));
            Assert.That(t1["ATG"], Is.EqualTo(1.0).Within(Tol));
            Assert.That(t1, Has.Count.EqualTo(64));
            // table 2, c = 0: *_AG {AGA 1, AGG 1} → 1, 1; *_TA {TAA 1, TAG 1} → 1, 1.
            Assert.That(new[] { t2["AGA"], t2["AGG"], t2["TAA"], t2["TAG"] }, Is.All.EqualTo(1.0).Within(Tol));
        });
    }

    // M25 — AminoAcid level is the existing pseudocount overload; c = 0 is the CodonW overload.
    [TestCase("ATGGCTGCAGCGAAAAAGTTTGGCGGTGGGCGTAGACGAAGCTCATCTTGA", 1, 0.0)]
    [TestCase("ATGGCTGCAGCGAAAAAGTTTGGCGGTGGGCGTAGACGAAGCTCATCTTGA", 2, 1.0)]
    [TestCase("CTTCTCCTACTGTTATTGCTGCTGCTG", 11, 0.5)]
    public void CalculateRscu_AminoAcidLevel_EqualsExistingOverloads(string gene, int table, double c)
    {
        var code = GeneticCode.GetByTableNumber(table);
        var level = CodonUsageAnalyzer.CalculateRscu(gene, code, c, RscuFamilyLevel.AminoAcid);
        Assert.That(level, Is.EquivalentTo(CodonUsageAnalyzer.CalculateRscu(gene, code, c)));
        if (c == 0)
            Assert.That(level, Is.EquivalentTo(CodonUsageAnalyzer.CalculateRscu(gene, code)));
    }

    // M26 — string, DnaSequence and counts overloads agree at the sub-family level.
    [Test]
    public void CalculateRscu_Subfamily_OverloadsAgree()
    {
        const string gene = "ATGCTGTTAAGCTCGCGAAGGTAA";
        var code = GeneticCode.Standard;
        var s = CodonUsageAnalyzer.CalculateRscu(gene.ToLowerInvariant().Replace('T', 'U'), code, 1.0,
            RscuFamilyLevel.Subfamily);
        var d = CodonUsageAnalyzer.CalculateRscu(new DnaSequence(gene), code, 1.0, RscuFamilyLevel.Subfamily);
        var k = CodonUsageAnalyzer.CalculateRscu(CodonUsageAnalyzer.CountCodons(gene), code, 1.0,
            RscuFamilyLevel.Subfamily);
        Assert.Multiple(() =>
        {
            Assert.That(s, Is.EquivalentTo(d));
            Assert.That(k, Is.EquivalentTo(d));
        });
    }

    // S9 — argument checks.
    [Test]
    public void CalculateRscu_Subfamily_InvalidArguments()
    {
        Assert.Multiple(() =>
        {
            Assert.That(() => CodonUsageAnalyzer.CalculateRscu("ATG", GeneticCode.Standard, 0, (RscuFamilyLevel)5),
                NUnit.Framework.Throws.TypeOf<ArgumentOutOfRangeException>());
            Assert.That(() => CodonUsageAnalyzer.CalculateRscu("ATG", GeneticCode.Standard, -1, RscuFamilyLevel.Subfamily),
                NUnit.Framework.Throws.TypeOf<ArgumentOutOfRangeException>());
            Assert.That(() => CodonUsageAnalyzer.CalculateRscu("ATG", null!, 0, RscuFamilyLevel.Subfamily),
                NUnit.Framework.Throws.ArgumentNullException);
            Assert.That(() => CodonUsageAnalyzer.CalculateRscu((DnaSequence)null!, GeneticCode.Standard, 0,
                RscuFamilyLevel.Subfamily), NUnit.Framework.Throws.ArgumentNullException);
            Assert.That(() => CodonUsageAnalyzer.CalculateRscu((IReadOnlyDictionary<string, int>)null!,
                GeneticCode.Standard, 0, RscuFamilyLevel.Subfamily), NUnit.Framework.Throws.ArgumentNullException);
            Assert.That(CodonUsageAnalyzer.CalculateRscu("", GeneticCode.Standard, 0, RscuFamilyLevel.Subfamily),
                Is.Empty);
        });
    }
}
