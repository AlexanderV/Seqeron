using NUnit.Framework;
using Seqeron.Genomics.Alignment;
using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Tests.Unit.Alignment;

/// <summary>
/// ALIGN-MULTI-001 (audit A2-14d, B09 F105): the alignment files of ClustalW 2.1 - <see cref="ClustalWAlignmentWriter"/> in the seven formats (<c>-OUTPUT=CLUSTAL | GCG | GDE | PHYLIP | PIR | NEXUS | FASTA</c>) with
/// <c>-OUTORDER</c>, <c>-CASE</c>, <c>-SEQNOS</c>, <c>-SEQNO_RANGE</c> and <c>-RANGE</c>, driven by the command line (<see cref="ClustalWCommandLine"/>), against the installed <c>clustalw</c> 2.1:
/// TestData/Alignment/clustalw_output_oracle.tsv (scripts/alignment/generate_clustalw_output_oracle.py), 2 310 files (330 per format; names of 1 - 40 characters, names with /start-end, every alignment option variant);
/// every output file must be equal byte for byte.
/// </summary>
[TestFixture]
[Category("Alignment")]
[Category("ALIGN-MULTI-001")]
public class SequenceAligner_ClustalWOutput_Tests
{
    private sealed record Row(string Id, string Type, string Format, string[] Names, string[] Seqs, string[] Args, string File);

    private static string Unescape(string s)
    {
        var sb = new System.Text.StringBuilder();
        bool escaped = false;
        foreach (char c in s)
        {
            if (escaped)
            {
                sb.Append(c switch { 'n' => '\n', 't' => '\t', _ => c });
                escaped = false;
            }
            else if (c == '\\')
            {
                escaped = true;
            }
            else
            {
                sb.Append(c);
            }
        }

        return sb.ToString();
    }

    private static readonly Lazy<List<Row>> Oracle = new(() =>
    {
        const string name = "Seqeron.Genomics.Tests.TestData.Alignment.clustalw_output_oracle.tsv";
        using var stream = TestDataResource.Open(name)
            ?? throw new InvalidOperationException($"Embedded oracle '{name}' not found.");
        using var reader = new StreamReader(stream);
        reader.ReadLine();
        var list = new List<Row>();
        string? line;
        while ((line = reader.ReadLine()) is not null)
        {
            if (line.Length == 0)
                continue;
            var f = line.Split('\t');
            list.Add(new Row(f[0], f[1], f[2], f[3].Split(';'), f[4].Split(';'), f[5].Split(' ', StringSplitOptions.RemoveEmptyEntries), Unescape(f[6])));
        }

        return list;
    });

    private static string Produce(Row r)
    {
        if (r.Args.Contains("-CONVERT"))
        {
            // clustalw -CONVERT: the unaligned input in the chosen format (rows of different lengths, weights 100, input order)
            string fasta = string.Concat(r.Names.Zip(r.Seqs, (n, q) => $">{n}\n{q}\n"));
            var run = ClustalWProgram.Run(["-INFILE=x.fa", .. r.Args], new Dictionary<string, string> { ["x.fa"] = fasta });
            return run.Files.Values.Single();
        }

        var cmd = ClustalWCommandLine.Parse(["-INFILE=x.fa", .. r.Args]);
        var seqs = r.Names.Zip(r.Seqs, (n, s) => new ClustalWSequence(n, s)).ToList();
        var result = r.Type == "dna" ? SequenceAligner.RunClustalW(seqs, cmd.Dna, cmd.Run) : SequenceAligner.RunClustalW(seqs, cmd.Protein, cmd.Run);
        return ClustalWAlignmentWriter.Write(result.Alignment!, cmd.Run.Output, cmd.Run.Structure);
    }

    [TestCase("CLUSTAL")]
    [TestCase("GCG")]
    [TestCase("GDE")]
    [TestCase("PHYLIP")]
    [TestCase("PIR")]
    [TestCase("NEXUS")]
    [TestCase("FASTA")]
    public void AlignmentFile_EqualsClustalW21ByteForByte_OnAllOracleSets(string format)
    {
        var rows = Oracle.Value.Where(r => r.Format == format).ToList();
        Assert.That(rows.Count(r => !r.Args.Contains("-CONVERT")), Is.GreaterThanOrEqualTo(300));
        var failures = new List<string>();
        foreach (var r in rows)
        {
            string got = Produce(r);
            if (got != r.File)
                failures.Add($"set {r.Id} ({r.Type}, {string.Join(' ', r.Args)})");
        }

        Assert.That(failures, Is.Empty, $"{failures.Count} / {rows.Count} files differ: {string.Join("; ", failures.Take(3))}");
    }

    [Test]
    public void Convert_UnalignedInput_IsWrittenInEveryFormat_AsTheProgramDoes()
    {
        var rows = Oracle.Value.Where(r => r.Args.Contains("-CONVERT")).ToList();
        Assert.That(rows.Count, Is.GreaterThanOrEqualTo(300));
        var failures = rows.Where(r => Produce(r) != r.File).Select(r => $"set {r.Id} ({r.Format}, {string.Join(' ', r.Args)})").ToList();
        Assert.That(failures, Is.Empty, $"{failures.Count} / {rows.Count} converted files differ: {string.Join("; ", failures.Take(3))}");
    }

    [Test]
    public void Oracle_CoversTheOptionsOfEveryFormat()
    {
        var rows = Oracle.Value;
        Assert.Multiple(() =>
        {
            Assert.That(rows.Count(r => r.Args.Contains("-SEQNO_RANGE=ON")), Is.GreaterThan(100));
            Assert.That(rows.Count(r => r.Args.Any(a => a.StartsWith("-RANGE=", StringComparison.Ordinal))), Is.GreaterThan(300));
            Assert.That(rows.Count(r => r.Args.Contains("-SEQNOS=ON")), Is.GreaterThan(100));
            Assert.That(rows.Count(r => r.Args.Contains("-CASE=UPPER")), Is.GreaterThan(100));
            Assert.That(rows.Count(r => r.Args.Contains("-OUTORDER=INPUT")), Is.GreaterThan(300));
            Assert.That(rows.Count(r => r.Type == "dna"), Is.GreaterThan(1000));
            Assert.That(rows.Count(r => r.Type == "protein"), Is.GreaterThan(1000));
        });
    }

    private static readonly ClustalWAlignment Example = new(
        ["human", "chimp", "mouse", "rat"],
        ["ACGTACGTTTGACCAGTAGGCATCGATGGATTACA-", "ACGTACGATTGACCAGTAGCCATCGATGGATTACA-", "ACGTAGGATTGACC-GTAGCCATCGATGG-TTACAC", "ACGTAGGATTGACC-GTAGCCATCGAAGG-TTACAC"],
        [0, 1, 2, 3],
        false,
        [357, 1, 214, 428]);

    [Test]
    public void Gcg_OfAFourSequenceAlignment_IsTheProgramsMsfFile()
    {
        // clustalw -INFILE=ex.fa -TYPE=DNA -ALIGN -OUTPUT=GCG -OUTORDER=INPUT (installed 2.1+lgpl-7): checksums, weights 35.7 / 0.1 / 21.4 / 42.8 and '.' for gaps
        string expected = "PileUp\n\n\n\n   MSF:   36  Type: N    Check:  6005   .. \n\n Name: human oo  Len:   36  Check:  6841  Weight:  35.7\n Name: chimp oo  Len:   36  Check:  6609  Weight:  0.1\n"
            + " Name: mouse oo  Len:   36  Check:  6534  Weight:  21.4\n Name: rat oo  Len:   36  Check:  6021  Weight:  42.8\n\n//\n\n\n\n"
            + "human           ACGTACGTTT GACCAGTAGG CATCGATGGA TTACA.\nchimp           ACGTACGATT GACCAGTAGC CATCGATGGA TTACA.\nmouse           ACGTAGGATT GACC.GTAGC CATCGATGG. TTACAC\n"
            + "rat             ACGTAGGATT GACC.GTAGC CATCGAAGG. TTACAC\n\n";
        string got = ClustalWAlignmentWriter.Write(Example, new ClustalWOutputOptions { Format = ClustalWOutputFormat.Gcg });
        Assert.That(got, Is.EqualTo(expected));
    }

    [Test]
    public void Clustal_WithSequenceNumbers_ListsTheRowsInAlignedOrderWithTheConservationLine()
    {
        // clustalw -INFILE=ex.fa -TYPE=DNA -ALIGN -OUTPUT=CLUSTAL -SEQNOS=ON: mouse, rat, chimp, human (the order of the progressive alignment)
        var seqs = new[]
        {
            new ClustalWSequence("human", "ACGTACGTTTGACCAGTAGGCATCGATGGATTACA"),
            new ClustalWSequence("chimp", "ACGTACGATTGACCAGTAGCCATCGATGGATTACA"),
            new ClustalWSequence("mouse", "ACGTAGGATTGACCGTAGCCATCGATGGTTACAC"),
            new ClustalWSequence("rat", "ACGTAGGATTGACCGTAGCCATCGAAGGTTACAC"),
        };
        var r = SequenceAligner.RunClustalW(seqs);
        string got = ClustalWAlignmentWriter.Write(r.Alignment!, new ClustalWOutputOptions { SequenceNumbers = true });
        string expected = "CLUSTAL 2.1 multiple sequence alignment\n\n\nmouse           ACGTAGGATTGACC-GTAGCCATCGATGG-TTACAC 34\nrat             ACGTAGGATTGACC-GTAGCCATCGAAGG-TTACAC 34\n"
            + "chimp           ACGTACGATTGACCAGTAGCCATCGATGGATTACA- 35\nhuman           ACGTACGTTTGACCAGTAGGCATCGATGGATTACA- 35\n                ***** * ****** **** ****** ** ***** \n";
        Assert.That(got, Is.EqualTo(expected));
    }

    [Test]
    public void Fasta_Pir_Phylip_Nexus_Gde_OfTheExample_FollowTheProgramsLayout()
    {
        string fasta = ClustalWAlignmentWriter.Write(Example, new ClustalWOutputOptions { Format = ClustalWOutputFormat.Fasta });
        string pir = ClustalWAlignmentWriter.Write(Example, new ClustalWOutputOptions { Format = ClustalWOutputFormat.Pir });
        string phylip = ClustalWAlignmentWriter.Write(Example, new ClustalWOutputOptions { Format = ClustalWOutputFormat.Phylip });
        string nexus = ClustalWAlignmentWriter.Write(Example, new ClustalWOutputOptions { Format = ClustalWOutputFormat.Nexus });
        string gde = ClustalWAlignmentWriter.Write(Example, new ClustalWOutputOptions { Format = ClustalWOutputFormat.Gde, Lowercase = false });
        Assert.Multiple(() =>
        {
            Assert.That(fasta, Does.StartWith(">human\nACGTACGTTTGACCAGTAGGCATCGATGGATTACA-\n>chimp\n"));
            Assert.That(pir, Does.StartWith(">DL;human\n\nACGTACGTTTGACCAGTAGGCATCGATGGATTACA-\n*\n>DL;chimp\n"));
            Assert.That(phylip, Does.StartWith("     4     36\nhuman      ACGTACGTTT GACCAGTAGG CATCGATGGA TTACA-"));
            Assert.That(nexus, Does.StartWith("#NEXUS\nBEGIN DATA;\ndimensions ntax=4 nchar=36;\nformat missing=?\ninterleave datatype=DNA gap= -;\n\nmatrix\nhuman       ACGTACGTTTGACCAGTAGGCATCGATGGATTACA-\n"));
            Assert.That(nexus, Does.EndWith("\n;\nend;\n"));
            Assert.That(gde, Does.StartWith("#human\nACGTACGTTTGACCAGTAGGCATCGATGGATTACA-\n#chimp\n"));
        });
    }

    [Test]
    public void Range_WritesTheRequestedColumns_AndSeqnoRangeNamesTheResidues()
    {
        var o = new ClustalWOutputOptions { Format = ClustalWOutputFormat.Fasta, Range = (5, 14), SequenceNumberRange = true };
        string fasta = ClustalWAlignmentWriter.Write(Example, o);
        // the program computes the last residue number as start + lastRes - gaps - 1 with lastRes = the LAST column (14), not the number of columns: 5 + 14 - 1 = 18
        Assert.That(fasta, Does.StartWith(">human/5-18\nACGTTTGACC\n>chimp/5-18\nACGATTGACC\n"));
    }

    [Test]
    public void InvalidAlignment_OrRange_IsRefused()
    {
        Assert.Throws<ArgumentException>(() => ClustalWAlignmentWriter.Write(new ClustalWAlignment(["a"], ["AC", "AC"], [0], false, [1000])));
        Assert.Throws<ArgumentException>(() => ClustalWAlignmentWriter.Write(Example, new ClustalWOutputOptions { Range = (0, 5) }));
        Assert.Throws<ArgumentException>(() => ClustalWAlignmentWriter.Write(Example, new ClustalWOutputOptions { Range = (100, 120) }));
        Assert.Throws<ArgumentNullException>(() => ClustalWAlignmentWriter.Write(null!));
    }
}
