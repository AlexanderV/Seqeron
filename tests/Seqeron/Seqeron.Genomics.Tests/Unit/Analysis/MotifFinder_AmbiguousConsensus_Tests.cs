// MOTIF-CONS-001 (B09 audit A3-13d, F147) — EMBOSS consambig: MotifFinder.GenerateEmbossAmbiguousConsensus.
// Oracle: scripts/analysis/generate_consambig_oracle.py — 1400 random alignments (nucleotide with all IUPAC letters, U, X, ?, gaps, lower case; protein with B Z J X U O * and gaps;
// ragged rows; type auto / -snucleotide / -sprotein), consambig 6.6.0 output. The NUL character written for a column without a scored amino acid is '@' in the file.
// Evidence: docs/Evidence/MOTIF-CONS-001-Evidence.md.

namespace Seqeron.Genomics.Tests.Unit.Analysis;

[TestFixture]
public class MotifFinder_AmbiguousConsensus_Tests
{
    private sealed record Row(string Type, string[] Rows, string Consensus);

    private static readonly Lazy<List<Row>> OracleRows = new(ReadOracleRows);

    private static List<Row> ReadOracleRows()
    {
        const string name = "Seqeron.Genomics.Tests.TestData.Analysis.consambig_oracle.tsv";
        using var stream = System.Reflection.Assembly.GetExecutingAssembly().GetManifestResourceStream(name)
            ?? throw new InvalidOperationException($"Embedded oracle '{name}' not found.");
        using var reader = new StreamReader(stream);
        var rows = new List<Row>();
        string? line;
        while ((line = reader.ReadLine()) is not null)
        {
            if (line.Length == 0) continue;
            var f = line.Split('\t');
            rows.Add(new Row(f[0], f[1].Split(','), f[2].Replace('@', '\0')));
        }
        return rows;
    }

    private static ConsensusResidueType TypeOf(string t) => t switch { "nuc" => ConsensusResidueType.Nucleotide, "prot" => ConsensusResidueType.Protein, _ => ConsensusResidueType.Auto };

    [Test]
    public void GenerateEmbossAmbiguousConsensus_EqualsConsambig_AllRows()
    {
        var rows = OracleRows.Value;
        Assert.That(rows, Has.Count.EqualTo(1400));
        int nuc = 0, protein = 0, ragged = 0, gapColumns = 0, nulColumns = 0, ambiguous = 0, forced = 0;
        foreach (var r in rows)
        {
            string actual = MotifFinder.GenerateEmbossAmbiguousConsensus(r.Rows, TypeOf(r.Type));
            Assert.That(actual, Is.EqualTo(r.Consensus), $"{r.Type} {string.Join(",", r.Rows)}");
            if (r.Rows[0].All(c => "ACGTUBDHKMNRSVWXY?.~-".Contains(char.ToUpperInvariant(c)))) nuc++; else protein++;
            if (r.Rows.Select(x => x.Length).Distinct().Count() > 1) ragged++;
            if (r.Consensus.Any(char.IsLower)) gapColumns++;
            if (r.Consensus.Contains('\0')) nulColumns++;
            if (r.Consensus.Any(c => "RYSWKMBDHV".Contains(char.ToUpperInvariant(c)))) ambiguous++;
            if (r.Type != "auto") forced++;
        }

        Assert.That(nuc, Is.GreaterThan(700));
        Assert.That(protein, Is.GreaterThan(400));
        Assert.That(ragged, Is.GreaterThan(100));
        Assert.That(gapColumns, Is.GreaterThan(500));
        Assert.That(nulColumns, Is.GreaterThan(10));
        Assert.That(ambiguous, Is.GreaterThan(400));
        Assert.That(forced, Is.GreaterThan(300));
    }

    [Test]
    public void WorkedValues_RunOnTheBinary()
    {
        // consambig on AC-GTN / AC-GTA / ARGGT-  ->  AVgGTn (V = A|C|G; the gap lower-cases; N = any)
        Assert.That(MotifFinder.GenerateEmbossAmbiguousConsensus(["AC-GTN", "AC-GTA", "ARGGT-"]), Is.EqualTo("AVgGTn"));
        // protein: D|N = B, I|L = J, A|D has no letter = X, an all-gap column is the NUL character
        Assert.That(MotifFinder.GenerateEmbossAmbiguousConsensus(["MK-LA", "MK-LD", "MN-LE", "DK-LJ"]), Is.EqualTo("XX\0LX"));
        Assert.That(MotifFinder.GenerateEmbossAmbiguousConsensus(["DI", "NL"], ConsensusResidueType.Protein), Is.EqualTo("BJ"));
        Assert.That(MotifFinder.GenerateEmbossAmbiguousConsensus(["ACGT", "AC", "ACGT"]), Is.EqualTo("ACgt"));   // the short row is padded with gaps
        // an all-gap nucleotide column: x, written n
        Assert.That(MotifFinder.GenerateEmbossAmbiguousConsensus(["AC-GT", "AC-GA"]), Is.EqualTo("ACnGW"));
    }

    [Test]
    public void Validation_AndGapOnlyDifferenceFromTheVotingConsensus()
    {
        Assert.Throws<ArgumentException>(() => MotifFinder.GenerateEmbossAmbiguousConsensus(["ACGT"]));
        Assert.Throws<ArgumentException>(() => MotifFinder.GenerateEmbossAmbiguousConsensus(["ACGT", "AC"], ConsensusResidueType.Auto, padRaggedRows: false));
        Assert.Throws<ArgumentNullException>(() => MotifFinder.GenerateEmbossAmbiguousConsensus(null!));
        // unlike cons there is no vote: one deviating base widens the code
        Assert.That(MotifFinder.GenerateEmbossAmbiguousConsensus(["AAAA", "AAAA", "AAAA", "AAAC"]), Is.EqualTo("AAAM"));
    }
}
