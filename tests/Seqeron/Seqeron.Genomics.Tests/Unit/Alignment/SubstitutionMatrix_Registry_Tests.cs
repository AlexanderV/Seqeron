using System.Globalization;
using NUnit.Framework;
using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Tests.Unit.Alignment;

/// <summary>
/// ALIGN-STATS-001 (audit A2-15, B09 F90): the named built-in matrices - every matrix of the EMBOSS 6.6.0 data directory (<c>Matrices.protein</c> / <c>Matrices.nucleotide</c> + <c>ENUC.4.4</c>), cell for cell
/// against the files (<c>emboss_matrices_oracle.tsv</c>, <c>scripts/alignment/generate_emboss_matrices.py</c>) and, for the NCBI aliases, against Biopython 1.88 <c>substitution_matrices</c> (hard-coded probe cells).
/// </summary>
[TestFixture]
[Category("Alignment")]
[Category("ALIGN-STATS-001")]
public class SubstitutionMatrix_Registry_Tests
{
    private sealed record Row(string Name, string Alphabet, int[] Scores);

    private static readonly Lazy<List<Row>> Oracle = new(() =>
    {
        const string name = "Seqeron.Genomics.Tests.TestData.Alignment.emboss_matrices_oracle.tsv";
        using var stream = System.Reflection.Assembly.GetExecutingAssembly().GetManifestResourceStream(name) ?? throw new InvalidOperationException(name);
        using var reader = new StreamReader(stream);
        reader.ReadLine();
        var rows = new List<Row>();
        string? line;
        while ((line = reader.ReadLine()) is not null)
        {
            if (line.Length == 0)
                continue;
            var f = line.Split('\t');
            rows.Add(new Row(f[0], f[1], f[2].Split(' ').Select(x => int.Parse(x, CultureInfo.InvariantCulture)).ToArray()));
        }

        return rows;
    });

    public static IEnumerable<TestCaseData> Names() => Oracle.Value.Select(r => new TestCaseData(r.Name).SetArgDisplayNames(r.Name));

    [Test]
    public void Oracle_Has70Matrices_EveryBuiltInNameIsCovered()
    {
        Assert.Multiple(() =>
        {
            Assert.That(Oracle.Value, Has.Count.EqualTo(70));
            Assert.That(SubstitutionMatrix.BuiltInNames, Is.EquivalentTo(Oracle.Value.Select(r => r.Name)));
            Assert.That(Oracle.Value.Count(r => r.Name.StartsWith("EBLOSUM")), Is.EqualTo(16), "EBLOSUM30 ... 90, EBLOSUM62-12, EBLOSUMN");
            Assert.That(Oracle.Value.Count(r => r.Name.StartsWith("EPAM")), Is.EqualTo(50), "EPAM10 ... EPAM500");
        });
    }

    [TestCaseSource(nameof(Names))]
    public void Get_EqualsTheEmbossFile_CellForCell(string name)
    {
        var r = Oracle.Value.First(x => x.Name == name);
        var m = SubstitutionMatrix.Get(name);
        int n = r.Alphabet.Length;
        Assert.That(m.Alphabet, Is.EqualTo(r.Alphabet));
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++)
                Assert.That(m.GetScore(r.Alphabet[i], r.Alphabet[j]), Is.EqualTo(r.Scores[i * n + j]), $"{name} {r.Alphabet[i]}/{r.Alphabet[j]}");
    }

    [Test]
    public void Get_IsCaseInsensitive_AcceptsNcbiNames_AndReturnsTheStaticInstances()
    {
        Assert.Multiple(() =>
        {
            Assert.That(SubstitutionMatrix.Get("ebloSUM62"), Is.SameAs(SubstitutionMatrix.Blosum62));
            Assert.That(SubstitutionMatrix.Get("BLOSUM62"), Is.SameAs(SubstitutionMatrix.Blosum62));
            Assert.That(SubstitutionMatrix.Get("pam250"), Is.SameAs(SubstitutionMatrix.Pam250));
            Assert.That(SubstitutionMatrix.Get("EDNAFULL"), Is.SameAs(SubstitutionMatrix.EdnaFull));
            Assert.That(SubstitutionMatrix.Get("NUC.4.4"), Is.SameAs(SubstitutionMatrix.EdnaFull), "the MCP's historic name for EDNAFULL");
            Assert.That(SubstitutionMatrix.Get("BLOSUM45"), Is.SameAs(SubstitutionMatrix.Get("EBLOSUM45")));
            Assert.That(SubstitutionMatrix.Get("PAM30").Name, Is.EqualTo("EPAM30"));
        });
    }

    /// <summary>Biopython 1.88 substitution_matrices.load(name)[a, b]; every probe equals the EMBOSS file as well.</summary>
    [TestCase("BLOSUM45", 'A', 'A', 5)]
    [TestCase("BLOSUM45", 'W', 'W', 15)]
    [TestCase("BLOSUM45", 'C', 'C', 12)]
    [TestCase("BLOSUM50", 'W', 'W', 15)]
    [TestCase("BLOSUM50", 'A', 'R', -2)]
    [TestCase("BLOSUM80", 'C', 'C', 13)]
    [TestCase("BLOSUM80", 'W', 'Y', 3)]
    [TestCase("BLOSUM90", 'W', 'W', 11)]
    [TestCase("PAM30", 'A', 'A', 6)]
    [TestCase("PAM30", 'W', 'W', 13)]
    [TestCase("PAM30", 'C', 'C', 10)]
    [TestCase("PAM70", 'A', 'A', 5)]
    [TestCase("PAM70", 'W', 'W', 13)]
    [TestCase("PAM250", 'W', 'W', 17)]
    [TestCase("NUC.4.4", 'A', 'A', 5)]
    [TestCase("NUC.4.4", 'A', 'T', -4)]
    public void NcbiAlias_ProbeCells_EqualBiopython(string name, char a, char b, int expected) =>
        Assert.That(SubstitutionMatrix.Get(name).GetScore(a, b), Is.EqualTo(expected));

    [Test]
    public void IsNucleotide_SeparatesTheNucleotideMatrices()
    {
        var nuc = SubstitutionMatrix.BuiltInNames.Where(n => SubstitutionMatrix.Get(n).IsNucleotide).ToArray();
        Assert.That(nuc, Is.EquivalentTo(new[] { "EDNAFULL", "EDNAMAT", "EDNASIMPLE", "ENUC.4.4" }));
        Assert.That(SubstitutionMatrix.Blosum62.IsNucleotide, Is.False);
    }

    [Test]
    public void Unknown_NullAndEmpty_AreRejected()
    {
        Assert.Multiple(() =>
        {
            Assert.That(SubstitutionMatrix.TryGet("BLOSUM999", out _), Is.False);
            Assert.That(SubstitutionMatrix.TryGet(null, out _), Is.False);
            Assert.That(SubstitutionMatrix.TryGet(" ", out _), Is.False);
            Assert.That(SubstitutionMatrix.TryGet("ENUC.4.2", out _), Is.False, "rectangular 15 x 4 table, not a pairwise matrix");
            Assert.Throws<ArgumentException>(() => SubstitutionMatrix.Get("nope"));
        });
    }
}
