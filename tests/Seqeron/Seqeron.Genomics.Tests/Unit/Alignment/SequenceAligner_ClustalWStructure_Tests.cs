using NUnit.Framework;
using Seqeron.Genomics.Alignment;
using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Tests.Unit.Alignment;

/// <summary>
/// ALIGN-MULTI-001 (audit A2-14c, B09 F106): the secondary-structure gap-penalty masks of ClustalW 2.1's profile commands. Profile files in CLUSTAL format (<c>!SS_</c> / <c>!GM_</c> lines) or GDE format
/// (<c>"SS_</c> / <c>"GM_</c> entries) are read by <see cref="ClustalWProfileFileReader"/>; <see cref="SequenceAligner.RunClustalWProfile(IReadOnlyList{ClustalWSequence}, IReadOnlyList{ClustalWSequence}, ClustalWOptions?, ClustalWRunOptions?)"/> /
/// <see cref="SequenceAligner.RunClustalWSequences(IReadOnlyList{ClustalWSequence}, IReadOnlyList{ClustalWSequence}, ClustalWOptions?, ClustalWRunOptions?)"/> use the masks, <see cref="ClustalWAlignmentWriter"/> writes them.
/// Oracle: TestData/Alignment/clustalw_structure_oracle.tsv (scripts/alignment/generate_clustalw_structure_oracle.py), 1 200 profile pairs with -NOSECSTR1/2, -SECSTROUT and the eight gap / end parameters; the whole output file
/// (rows and structure lines, CLUSTAL / GDE / FASTA / PIR) must be equal.
/// </summary>
[TestFixture]
[Category("Alignment")]
[Category("ALIGN-MULTI-001")]
public class SequenceAligner_ClustalWStructure_Tests
{
    private sealed record Row(string Id, string Command, string Type, string InputFormat, string File1, string File2, string[] Args, string Format, string Output);

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
        const string name = "Seqeron.Genomics.Tests.TestData.Alignment.clustalw_structure_oracle.tsv";
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
            list.Add(new Row(f[0], f[1], f[2], f[3], Unescape(f[4]), Unescape(f[5]), f[6].Split(' ', StringSplitOptions.RemoveEmptyEntries), f[7], Unescape(f[8])));
        }

        return list;
    });

    private static (string Output, bool RowsChangedByStructure) Produce(Row r)
    {
        bool dna = r.Type == "dna";
        var cmd = ClustalWCommandLine.Parse(["-TYPE=" + r.Type.ToUpperInvariant(), "-OUTPUT=" + r.Format, "-PROFILE1=a.in", "-PROFILE2=b.in", r.Command == "profile" ? "-PROFILE" : "-SEQUENCES", .. r.Args], new Dictionary<string, string>());
        var p1 = r.InputFormat == "clustal" ? ClustalWProfileFileReader.ParseClustal(r.File1) : ClustalWProfileFileReader.ParseGde(r.File1, dna);
        var p2 = r.InputFormat == "clustal" ? ClustalWProfileFileReader.ParseClustal(r.File2) : ClustalWProfileFileReader.ParseGde(r.File2, dna);
        var run = cmd.Run with { ProfileStructure1 = p1.Structure, ProfileStructure2 = p2.Structure };
        ClustalWRunResult With(ClustalWRunOptions o)
        {
            var a = p1.Sequences.ToList();
            var b = p2.Sequences.ToList();
            return dna
                ? (r.Command == "profile" ? SequenceAligner.RunClustalWProfile(a, b, cmd.Dna, o) : SequenceAligner.RunClustalWSequences(a, b, cmd.Dna, o))
                : (r.Command == "profile" ? SequenceAligner.RunClustalWProfile(a, b, cmd.Protein, o) : SequenceAligner.RunClustalWSequences(a, b, cmd.Protein, o));
        }

        var res = With(run);
        var plain = With(run with { ProfileStructure1 = null, ProfileStructure2 = null });
        return (ClustalWAlignmentWriter.Write(res.Alignment!, cmd.Run.Output, cmd.Run.Structure), !plain.Alignment!.Rows.SequenceEqual(res.Alignment!.Rows));
    }

    [TestCase("profile", "dna")]
    [TestCase("profile", "protein")]
    [TestCase("sequences", "dna")]
    [TestCase("sequences", "protein")]
    public void OutputFile_EqualsClustalW21_OnAllOracleSets(string command, string type)
    {
        var rows = Oracle.Value.Where(r => r.Command == command && r.Type == type).ToList();
        Assert.That(rows.Count, Is.GreaterThanOrEqualTo(300));
        var failures = new List<string>();
        int changed = 0;
        foreach (var r in rows)
        {
            var (output, structureChangesRows) = Produce(r);
            if (structureChangesRows)
                changed++;
            if (output != r.Output)
                failures.Add($"set {r.Id} ({r.InputFormat}, {r.Format}, {string.Join(' ', r.Args)})");
        }

        Assert.Multiple(() =>
        {
            Assert.That(failures, Is.Empty, $"{failures.Count} / {rows.Count} files differ: {string.Join("; ", failures.Take(3))}");
            Assert.That(changed, Is.GreaterThan(80), "the masks must change the rows of a good share of the sets (otherwise the oracle would not test them)");
        });
    }

    [Test]
    public void Oracle_UsesBothFileFormats_AllFourOutputsAndTheStructureOptions()
    {
        var rows = Oracle.Value;
        Assert.Multiple(() =>
        {
            Assert.That(rows.Count(r => r.InputFormat == "gde"), Is.GreaterThan(200));
            Assert.That(rows.Count(r => r.InputFormat == "clustal"), Is.GreaterThan(900));
            foreach (string f in new[] { "CLUSTAL", "GDE", "FASTA", "PIR" })
                Assert.That(rows.Any(r => r.Format == f), Is.True, f);
            foreach (string a in new[] { "-NOSECSTR1", "-NOSECSTR2", "-SECSTROUT=MASK", "-SECSTROUT=BOTH", "-SECSTROUT=NONE", "-HELIXGAP", "-STRANDGAP", "-LOOPGAP", "-TERMINALGAP", "-HELIXENDIN", "-HELIXENDOUT", "-STRANDENDIN", "-STRANDENDOUT" })
                Assert.That(rows.Any(r => r.Args.Any(x => x.StartsWith(a, StringComparison.Ordinal))), Is.True, a);
        });
    }

    private const string ProfileFile = "CLUSTAL 2.1 multiple sequence alignment\n\n\n!SS_helix_strand  ..AAAAAAAA$....BBBBBB%.....\nseq1            GATTACAGATTACAGATTACAGATTACAG\n"
        + "seq2            GATTACAGATTACAGATTTCAGATTACAG\nseq3            GATTACAGATTACAGATTACAGATTAC-G\n                ***************** ***** *** \n\n";

    [Test]
    public void HelixAndStrandProfile_GapPenaltiesAndStructureLines_AreTheProgramsOutput()
    {
        // clustalw -PROFILE1=sa.aln -PROFILE2=sb.fa -PROFILE -TYPE=DNA -OUTPUT=CLUSTAL -SECSTROUT=BOTH (installed 2.1+lgpl-7)
        var p1 = ClustalWProfileFileReader.ParseClustal(ProfileFile);
        var p2 = new[] { new ClustalWSequence("seq4", "GATTACAGGATTACAGATTACAGATTACAG"), new ClustalWSequence("seq5", "GATACAGATTACAGATTACAGATTACAGG") };
        var run = new ClustalWRunOptions { ProfileStructure1 = p1.Structure, Output = new ClustalWOutputOptions { Structure = ClustalWStructureOutput.Both } };
        var r = SequenceAligner.RunClustalWProfile(p1.Sequences.ToList(), p2, null, run);
        string expected = "CLUSTAL 2.1 multiple sequence alignment\n\n\n!SS_helix_strand  ..aaaAAAaaa....bBBBBBb.......--\n!GM_helix_strand  11222444222111224444422111111--\n"
            + "seq1            GATTACAGATTACAGATTACAGATTACAG--\nseq2            GATTACAGATTACAGATTTCAGATTACAG--\nseq3            GATTACAGATTACAGATTACAGATTAC-G--\n"
            + "seq4            -GATTACAGGATTACAGATTACAGATTACAG\nseq5            -GATACAGATTACAGATTACAGATTACAGG-\n                   *         * *    * *        \n";
        Assert.Multiple(() =>
        {
            Assert.That(p1.Structure, Is.EqualTo(new ClustalWProfileStructure("helix_strand", "..AAAAAAAA$....BBBBBB%.....", null)));
            Assert.That(ClustalWAlignmentWriter.Write(r.Alignment!, run.Output, run.Structure), Is.EqualTo(expected));
        });
        // with -HELIXGAP=9 the GDE file shows the stronger helix penalty (clustalw ... -OUTPUT=GDE -SECSTROUT=BOTH -HELIXGAP=9)
        var strong = SequenceAligner.RunClustalWProfile(p1.Sequences.ToList(), p2, null, run with { Structure = new ClustalWStructureOptions { HelixPenalty = 9 }, Output = new ClustalWOutputOptions { Format = ClustalWOutputFormat.Gde, Structure = ClustalWStructureOutput.Both } });
        string gde = ClustalWAlignmentWriter.Write(strong.Alignment!, new ClustalWOutputOptions { Format = ClustalWOutputFormat.Gde, Structure = ClustalWStructureOutput.Both });
        Assert.That(gde, Does.EndWith("\"SS_helix_strand\n..aaaAAAaaa....bBBBBBb.......--\n\"GM_helix_strand\n11222999222111224444422111111--\n"));
    }

    [Test]
    public void ProfileFileReader_Clustal_TakesTheFirstStructureLine_AndReportsAMissingLine()
    {
        var both = ClustalWProfileFileReader.ParseClustal("CLUSTAL 2.1 multiple sequence alignment\n\n\n!GM_m1  1234\n!SS_s1  AAAA\na               ACGT\nb               ACGA\n                ***.\n\n");
        Assert.Multiple(() =>
        {
            Assert.That(both.Sequences.Select(s => s.Name), Is.EqualTo(new[] { "a", "b" }));
            Assert.That(both.Sequences.Select(s => s.Residues), Is.EqualTo(new[] { "ACGT", "ACGA" }));
            Assert.That(both.Structure, Is.EqualTo(new ClustalWProfileStructure("m1", null, "1234")), "the first structure line of the file wins");
        });

        var broken = ClustalWProfileFileReader.ParseClustal("CLUSTAL 2.1 multiple sequence alignment\n\n\n!SS_s1  AAAA\na               ACGT\nb               ACGA\n                ***.\n\na               ACGT\nb               ACGA\n                ****\n\n");
        Assert.Multiple(() =>
        {
            Assert.That(broken.Structure, Is.Null);
            Assert.That(broken.StructureProblem, Is.EqualTo("bad secondary structure format"));
            Assert.That(broken.Sequences[0].Residues, Is.EqualTo("ACGTACGT"));
        });
    }

    [Test]
    public void ProfileFileReader_Gde_ReadsSequencesAndEntries()
    {
        var f = ClustalWProfileFileReader.ParseGde("#a\nacgt\n#b\nacga\n\"SS_helix\n.AAA\n", dna: true);
        Assert.Multiple(() =>
        {
            Assert.That(f.Sequences.Select(s => s.Name + ":" + s.Residues), Is.EqualTo(new[] { "a:ACGT", "b:ACGA" }));
            Assert.That(f.Structure, Is.EqualTo(new ClustalWProfileStructure("helix", ".AAA", null)));
        });
        Assert.That(ClustalWProfileFileReader.ParseGde("%p\nMKV\n", dna: true).Sequences, Is.Empty, "a protein file has no '#' entries");
    }

    [Test]
    public void NoSecStr_TurnsTheMaskOffForTheAlignment_AndForItsOutput()
    {
        var p1 = ClustalWProfileFileReader.ParseClustal(ProfileFile);
        var p2 = new[] { new ClustalWSequence("seq4", "GATTACAGGATTACAGATTACAGATTACAG"), new ClustalWSequence("seq5", "GATACAGATTACAGATTACAGATTACAGG") };
        var run = new ClustalWRunOptions { ProfileStructure1 = p1.Structure, Structure = new ClustalWStructureOptions { UseStructure1 = false } };
        var off = SequenceAligner.RunClustalWProfile(p1.Sequences.ToList(), p2, null, run);
        var none = SequenceAligner.RunClustalWProfile(p1.Sequences.ToList(), p2);
        Assert.Multiple(() =>
        {
            Assert.That(off.Alignment!.Rows, Is.EqualTo(none.Alignment!.Rows));
            Assert.That(ClustalWAlignmentWriter.Write(off.Alignment, run.Output, run.Structure), Does.Not.Contain("!SS_"));
        });
    }
}
