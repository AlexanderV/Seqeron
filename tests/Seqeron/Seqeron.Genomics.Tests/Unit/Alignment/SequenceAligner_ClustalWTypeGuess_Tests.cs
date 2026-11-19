using NUnit.Framework;
using Seqeron.Genomics.Alignment;
using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Tests.Unit.Alignment;

/// <summary>
/// ALIGN-MULTI-001 (audit A2-14, B09 F92): <see cref="SequenceAligner.GuessClustalWSequenceType"/> = ClustalW 2.1's residue-type guess when <c>-TYPE</c> is omitted
/// (<c>FileReader::readFile</c>: the flag is set from the SECOND sequence only; <c>Sequence::checkDNAFlag</c>: DNA when >= 85 % of its non-gap characters are A C G T U N).
/// Oracle: the unit (<c>bp</c> / <c>aa</c>) of the "Sequence n" lines printed by the installed <c>clustalw</c> 2.1 on 500 random sets
/// (<c>scripts/alignment/generate_clustalw_type_oracle.py</c>, <c>clustalw_type_oracle.tsv</c>).
/// </summary>
[TestFixture]
[Category("Alignment")]
[Category("ALIGN-MULTI-001")]
public class SequenceAligner_ClustalWTypeGuess_Tests
{
    private static List<(string[] Sequences, ClustalWSequenceType Expected)> ReadOracle()
    {
        const string name = "Seqeron.Genomics.Tests.TestData.Alignment.clustalw_type_oracle.tsv";
        using var stream = System.Reflection.Assembly.GetExecutingAssembly().GetManifestResourceStream(name) ?? throw new InvalidOperationException(name);
        using var reader = new StreamReader(stream);
        reader.ReadLine();
        var rows = new List<(string[], ClustalWSequenceType)>();
        string? line;
        while ((line = reader.ReadLine()) is not null)
        {
            if (line.Length == 0)
                continue;
            var f = line.Split('\t');
            rows.Add((f[0].Split(';'), f[1] == "dna" ? ClustalWSequenceType.Dna : ClustalWSequenceType.Protein));
        }

        return rows;
    }

    [Test]
    public void Guess_EqualsTheTypeTheInstalledClustalWPrints_On500Sets()
    {
        var rows = ReadOracle();
        var failures = rows.Where(r => SequenceAligner.GuessClustalWSequenceType(r.Sequences) != r.Expected)
            .Select(r => $"{string.Join(';', r.Sequences)} expected {r.Expected}").ToList();
        Assert.Multiple(() =>
        {
            Assert.That(rows, Has.Count.GreaterThanOrEqualTo(480));
            Assert.That(rows.Count(r => r.Expected == ClustalWSequenceType.Dna), Is.GreaterThan(100));
            Assert.That(rows.Count(r => r.Expected == ClustalWSequenceType.Protein), Is.GreaterThan(100));
            Assert.That(rows.Count(r => r.Sequences.Length == 1), Is.GreaterThan(20), "a single sequence is protein (the flag is only set at index 1)");
            Assert.That(failures, Is.Empty, string.Join("\n", failures.Take(3)));
        });
    }

    [Test]
    public void Guess_IsDecidedBySecondSequenceOnly_AndSingleSequenceIsProtein()
    {
        Assert.Multiple(() =>
        {
            Assert.That(SequenceAligner.GuessClustalWSequenceType(["MKVLAAGIV", "ACGTACGTAC", "MKVLAAGIV"]), Is.EqualTo(ClustalWSequenceType.Dna), "ClustalW printed 'bp' for all three (checked)");
            Assert.That(SequenceAligner.GuessClustalWSequenceType(["ACGTACGTAC", "ACGTRYGTAC"]), Is.EqualTo(ClustalWSequenceType.Protein), "8 of 10 = 80 % < 85 %");
            Assert.That(SequenceAligner.GuessClustalWSequenceType(["ACGTACGTAC"]), Is.EqualTo(ClustalWSequenceType.Protein));
            Assert.That(SequenceAligner.GuessClustalWSequenceType(["ACGTACGTAC", "acgtacgtac"]), Is.EqualTo(ClustalWSequenceType.Dna));
        });
    }
}
