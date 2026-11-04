// CODON-OPT-001 (B08 R14) — RemoveRestrictionSites / OptimizeSequence(avoidSites): linear cost, the site-length limit and eager validation.
//
// R14 (QUEUE.md): the pass re-scanned the whole sequence after every replaced codon and rebuilt the full candidate string for every candidate
// (O(n) each), so the cost was quadratic: the MCP tool remove_restriction_sites needed 8.2 s for 48 000 nt with the site N (10.2 s for 96 000 nt,
// 20.4 s for 96 000 nt of tandem GAATTC). The rewrite keeps one mutable buffer, restarts the scan only where a changed codon can have created
// or exposed an occurrence, and reads a candidate's effect off the buffer; the results are byte-identical:
//   * the golden TSV (TestData/Codon/remove_restriction_sites_golden.tsv) holds the outputs of the quadratic code (commit ea937485) for 464
//     random cases, replayed here; the new code also reproduced 3 x 3 105 further random cases of the same generator (seeds 11, 22, 33);
//   * the timing tests below fail on the old code by two orders of magnitude (the limits are generous: the new code needs 20-200 ms).
// B02's earlier finding stays true: with a GC constraint (BalancedOptimization / MinimizeSecondary through OptimizeSequence) every candidate is
// still scored against the whole sequence (the GC-breach score is a whole-sequence floating-point sum whose order is part of the result).
// O9-1 / O9-2: the sites are validated whatever the sequence is, under the caller's parameter name.

using System.Diagnostics;
using System.Globalization;
using System.Reflection;
using System.Text;

namespace Seqeron.Genomics.Tests.Unit.MolTools;

[TestFixture]
[Category("Codon")]
[Category("CODON-OPT-001")]
public class CodonOptimizer_RemoveRestrictionSites_Linear_Tests
{
    private const string Resource = "Seqeron.Genomics.Tests.TestData.Codon.remove_restriction_sites_golden.tsv";

    // New code: 20-200 ms for the sequences below; the quadratic code needed 10-20 s. Generous so that a slow runner does not fail it.
    private static readonly TimeSpan Limit = TimeSpan.FromSeconds(5);

    #region Golden replay

    private sealed record Row(int Id, string Mode, string Table, int Code, string Strategy, double GcMin, double GcMax, int? Window,
        string[] Sites, string Sequence, string Expected);

    private static readonly Lazy<List<Row>> Golden = new(LoadGolden);

    private static string Decode(string cell) => cell == "<empty>" ? "" : cell;

    private static List<Row> LoadGolden()
    {
        using var stream = Assembly.GetExecutingAssembly().GetManifestResourceStream(Resource)
            ?? throw new InvalidOperationException($"Embedded golden '{Resource}' not found.");
        using var reader = new StreamReader(stream);
        var rows = new List<Row>();
        string? line;
        while ((line = reader.ReadLine()) is not null)
        {
            if (line.Length == 0 || line[0] == '#')
                continue;
            var f = line.Split('\t');
            rows.Add(new Row(
                int.Parse(f[0], CultureInfo.InvariantCulture), f[1], f[2], int.Parse(f[3], CultureInfo.InvariantCulture), f[4],
                f[5] == "-" ? 0 : double.Parse(f[5], CultureInfo.InvariantCulture),
                f[6] == "-" ? 0 : double.Parse(f[6], CultureInfo.InvariantCulture),
                f[7] == "-" ? null : int.Parse(f[7], CultureInfo.InvariantCulture),
                f[8].Split(',').Select(Decode).ToArray(), Decode(f[9]), Decode(f[10])));
        }

        return rows;
    }

    private static CodonOptimizer.CodonUsageTable TableByName(string name) => name switch
    {
        "EColiK12" => CodonOptimizer.EColiK12,
        "Yeast" => CodonOptimizer.Yeast,
        "Human" => CodonOptimizer.Human,
        _ => throw new InvalidOperationException(name),
    };

    private static string Run(Row row)
    {
        try
        {
            var table = TableByName(row.Table);
            var code = GeneticCode.GetByTableNumber(row.Code);
            if (row.Mode == "R")
                return CodonOptimizer.RemoveRestrictionSites(row.Sequence, row.Sites, table, code);

            var strategy = Enum.Parse<CodonOptimizer.OptimizationStrategy>(row.Strategy);
            var r = CodonOptimizer.OptimizeSequence(row.Sequence, table, strategy, row.GcMin, row.GcMax, 0.15,
                CodonOptimizer.SecondaryStructureScoring.PairCount, row.Window, code, row.Sites);
            return r.OptimizedSequence + "|" + r.ChangedCodons;
        }
        catch (Exception e)
        {
            return "!" + e.GetType().Name;
        }
    }

    [Test]
    [Description("The 464 outputs of the quadratic implementation (commit ea937485) are reproduced byte for byte by the linear one")]
    public void Golden_EveryCase_IsReproducedByteForByte()
    {
        var rows = Golden.Value;
        Assert.That(rows, Has.Count.EqualTo(464));
        Assert.That(rows.Count(r => r.Mode == "R"), Is.EqualTo(324));
        Assert.That(rows.Count(r => r.Mode == "O"), Is.EqualTo(140));

        var mismatches = new List<string>();
        foreach (var row in rows)
        {
            string actual = Run(row);
            if (actual != row.Expected)
                mismatches.Add($"#{row.Id} {row.Mode} {row.Table} code {row.Code} sites [{string.Join(",", row.Sites)}] on {row.Sequence}: expected {row.Expected}, got {actual}");
        }

        Assert.That(mismatches, Is.Empty, string.Join(Environment.NewLine, mismatches.Take(10)));
    }

    [Test]
    public void Golden_Covers_ChangedAndUnchangedResults_AndTheRefusals()
    {
        var rows = Golden.Value;
        static string Complete(string s)
        {
            string rna = s.ToUpperInvariant().Replace('T', 'U');
            return rna[..(rna.Length - rna.Length % 3)];
        }

        int changed = rows.Count(r => r.Mode == "R" && !r.Expected.StartsWith('!') && r.Expected != Complete(r.Sequence));
        Assert.Multiple(() =>
        {
            Assert.That(changed, Is.GreaterThan(100), "cases in which codons were replaced");
            Assert.That(rows.Count(r => r.Expected.StartsWith('!')), Is.GreaterThanOrEqualTo(4), "refusals");
            Assert.That(rows.Where(r => r.Mode == "O").Select(r => r.Strategy).Distinct(), Is.EquivalentTo(new[] { "MaximizeCAI", "BalancedOptimization", "MinimizeSecondary" }));
        });
    }

    #endregion

    #region Linear cost (R14)

    private static string RandomDna(int codons, int seed)
    {
        var rnd = new Random(seed);
        var sb = new StringBuilder(codons * 3);
        for (int i = 0; i < codons; i++)
        {
            // Sense codons only, so that the stop codon cannot interrupt the protein comparison.
            string codon;
            do
            {
                codon = new string(new[] { "ACGT"[rnd.Next(4)], "ACGT"[rnd.Next(4)], "ACGT"[rnd.Next(4)] });
            } while (codon is "TAA" or "TAG" or "TGA");
            sb.Append(codon);
        }

        return sb.ToString();
    }

    private static string Protein(string rna)
    {
        var sb = new StringBuilder();
        for (int i = 0; i + 3 <= rna.Length; i += 3)
            sb.Append(GeneticCode.Standard.Translate(rna.Substring(i, 3)));
        return sb.ToString();
    }

    [Test]
    [Description("R14: the site N (every position matches, nothing can be removed) on 96 000 nt took 10.2 s")]
    public void SiteN_96000nt_FinishesQuickly_AndLeavesTheSequenceAsItWas()
    {
        string dna = RandomDna(32_000, 14);

        var sw = Stopwatch.StartNew();
        string result = CodonOptimizer.RemoveRestrictionSites(dna, ["N"], CodonOptimizer.EColiK12);
        sw.Stop();

        Assert.Multiple(() =>
        {
            Assert.That(sw.Elapsed, Is.LessThan(Limit), "elapsed");
            Assert.That(result, Is.EqualTo(dna.Replace('T', 'U')));
        });
    }

    [Test]
    [Description("R14: 16 000 tandem EcoRI sites (96 000 nt) took 20.4 s")]
    public void TandemEcoRI_96000nt_FinishesQuickly_AndEverySiteIsRemoved()
    {
        string dna = string.Concat(Enumerable.Repeat("GAATTC", 16_000));

        var sw = Stopwatch.StartNew();
        string result = CodonOptimizer.RemoveRestrictionSites(dna, ["GAATTC"], CodonOptimizer.EColiK12);
        sw.Stop();

        Assert.Multiple(() =>
        {
            Assert.That(sw.Elapsed, Is.LessThan(Limit), "elapsed");
            Assert.That(result, Has.Length.EqualTo(dna.Length));
            Assert.That(result, Does.Not.Contain("GAAUUC"), "no site is left");
            Assert.That(Protein(result), Is.EqualTo(Protein(dna.Replace('T', 'U'))), "the protein is preserved");
        });
    }

    [Test]
    [Description("R14: the longest permitted site (64 x N) on 96 000 nt took 10.3 s")]
    public void SiteOfTheLongestPermittedLength_96000nt_FinishesQuickly()
    {
        string dna = RandomDna(32_000, 64);
        string site = new('N', RestrictionAnalyzer.MaxRecognitionLength);

        var sw = Stopwatch.StartNew();
        string result = CodonOptimizer.RemoveRestrictionSites(dna, [site], CodonOptimizer.EColiK12);
        sw.Stop();

        Assert.Multiple(() =>
        {
            Assert.That(sw.Elapsed, Is.LessThan(Limit), "elapsed");
            Assert.That(result, Is.EqualTo(dna.Replace('T', 'U')));
        });
    }

    [Test]
    [Description("R14: ~23 scattered EcoRI sites in 96 000 random nt: all removed, protein kept, linear time")]
    public void ScatteredEcoRI_96000nt_AreAllRemoved_InLinearTime()
    {
        string dna = RandomDna(32_000, 99);
        string rna = dna.Replace('T', 'U');

        var sw = Stopwatch.StartNew();
        string result = CodonOptimizer.RemoveRestrictionSites(dna, ["GAATTC"], CodonOptimizer.EColiK12);
        sw.Stop();

        Assert.Multiple(() =>
        {
            Assert.That(sw.Elapsed, Is.LessThan(Limit), "elapsed");
            Assert.That(CountOccurrences(rna, "GAAUUC"), Is.GreaterThan(5), "the input does hold sites (test sanity)");
            Assert.That(result, Does.Not.Contain("GAAUUC"));
            Assert.That(Protein(result), Is.EqualTo(Protein(rna)), "the protein is preserved");
        });
    }

    [Test]
    [Description("R14: four enzymes (both strands, 7 patterns) in 96 000 random nt: protein kept, fewer sites, linear time")]
    public void SeveralEnzymes_96000nt_ReduceTheSites_InLinearTime()
    {
        string dna = RandomDna(32_000, 7);
        string[] sites = ["GAATTC", "GGATCC", "AAGCTT", "GGTCTC"];
        string[] patterns = sites.Concat(sites.Select(DnaSequence.GetReverseComplementString)).Distinct().ToArray();
        string rna = dna.Replace('T', 'U');

        var sw = Stopwatch.StartNew();
        string result = CodonOptimizer.RemoveRestrictionSites(dna, sites, CodonOptimizer.EColiK12);
        sw.Stop();

        int before = patterns.Sum(p => CountOccurrences(dna, p));
        int after = patterns.Sum(p => CountOccurrences(result.Replace('U', 'T'), p));
        Assert.Multiple(() =>
        {
            Assert.That(sw.Elapsed, Is.LessThan(Limit), "elapsed");
            Assert.That(before, Is.GreaterThan(50), "the input does hold sites (test sanity)");
            Assert.That(after, Is.LessThan(before / 10), $"{before} occurrences before, {after} after");
            Assert.That(Protein(result), Is.EqualTo(Protein(rna)), "the protein is preserved");
        });
    }

    private static int CountOccurrences(string text, string pattern)
    {
        int count = 0;
        for (int i = text.IndexOf(pattern, StringComparison.Ordinal); i >= 0; i = text.IndexOf(pattern, i + 1, StringComparison.Ordinal))
            count++;
        return count;
    }

    #endregion

    #region The site-length limit and eager validation (O9-1, O9-2)

    private static readonly string LongestSite = new('N', RestrictionAnalyzer.MaxRecognitionLength);
    private static readonly string TooLongSite = new('N', RestrictionAnalyzer.MaxRecognitionLength + 1);

    public static IEnumerable<TestCaseData> SequencesOfEveryShape()
    {
        yield return new TestCaseData("").SetName("empty");
        yield return new TestCaseData("AU").SetName("two_nucleotides");
        yield return new TestCaseData("AUG").SetName("one_codon");
        yield return new TestCaseData("AUGGCUAAAGAAUUCGGAUCC").SetName("seven_codons");
    }

    [TestCaseSource(nameof(SequencesOfEveryShape))]
    public void RemoveRestrictionSites_ASiteOverTheLimit_IsRefused_WhateverTheSequence(string sequence)
    {
        var e = Assert.Throws<ArgumentException>(() => CodonOptimizer.RemoveRestrictionSites(sequence, [TooLongSite], CodonOptimizer.EColiK12));

        Assert.Multiple(() =>
        {
            Assert.That(e!.ParamName, Is.EqualTo("restrictionSites"));
            Assert.That(e.Message, Does.Contain("64"));
            Assert.That(e.Message, Does.Contain("65"));
        });
    }

    [TestCaseSource(nameof(SequencesOfEveryShape))]
    public void RemoveRestrictionSites_ANonIupacSite_IsRefused_WhateverTheSequence(string sequence)
    {
        var e = Assert.Throws<ArgumentException>(() => CodonOptimizer.RemoveRestrictionSites(sequence, ["GAAXTC"], CodonOptimizer.EColiK12));

        Assert.Multiple(() =>
        {
            Assert.That(e!.ParamName, Is.EqualTo("restrictionSites"));
            Assert.That(e.Message, Does.Contain("GAAXTC"));
        });
    }

    [TestCaseSource(nameof(SequencesOfEveryShape))]
    public void RemoveRestrictionSites_ANullSiteList_IsRefused_WhateverTheSequence(string sequence)
    {
        var e = Assert.Throws<ArgumentNullException>(() => CodonOptimizer.RemoveRestrictionSites(sequence, null!, CodonOptimizer.EColiK12));

        Assert.That(e!.ParamName, Is.EqualTo("restrictionSites"));
    }

    [TestCaseSource(nameof(SequencesOfEveryShape))]
    public void RemoveRestrictionSites_TheLongestPermittedSite_AndRnaSpelling_AreAccepted_WhateverTheSequence(string sequence)
    {
        string result = CodonOptimizer.RemoveRestrictionSites(sequence, [LongestSite, "GAAUUC", ""], CodonOptimizer.EColiK12);

        Assert.That(result, Has.Length.EqualTo(sequence.Length - sequence.Length % 3));
    }

    [TestCaseSource(nameof(SequencesOfEveryShape))]
    public void OptimizeSequence_AnInvalidAvoidSite_IsRefused_WhateverTheSequence_UnderItsOwnParameterName(string sequence)
    {
        var nonIupac = Assert.Throws<ArgumentException>(() => CodonOptimizer.OptimizeSequence(sequence, CodonOptimizer.EColiK12, avoidSites: ["GAAXTC"]));
        var tooLong = Assert.Throws<ArgumentException>(() => CodonOptimizer.OptimizeSequence(sequence, CodonOptimizer.EColiK12, avoidSites: [TooLongSite]));

        Assert.Multiple(() =>
        {
            Assert.That(nonIupac!.ParamName, Is.EqualTo("avoidSites"));
            Assert.That(tooLong!.ParamName, Is.EqualTo("avoidSites"));
        });
    }

    [TestCaseSource(nameof(SequencesOfEveryShape))]
    public void OptimizeSequence_ValidAvoidSites_StillOptimize_WhateverTheSequence(string sequence)
    {
        var result = CodonOptimizer.OptimizeSequence(sequence, CodonOptimizer.EColiK12, avoidSites: [LongestSite, "GAATTC", ""]);

        Assert.That(result.OptimizedSequence, Has.Length.EqualTo(sequence.Length - sequence.Length % 3));
    }

    #endregion
}
