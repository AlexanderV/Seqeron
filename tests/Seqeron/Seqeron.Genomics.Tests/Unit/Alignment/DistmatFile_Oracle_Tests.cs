using System.IO.Compression;
using System.Text.Json;
using NUnit.Framework;
using Seqeron.Genomics.Analysis;

namespace Seqeron.Genomics.Tests.Unit.Alignment;

/// <summary>
/// A5-8b (B09 F184): the whole <c>-outfile</c> text of EMBOSS 6.6.0 <c>distmat</c> (title, method line, <c>Using ambiguity codes</c>, <c>Using base positions N in the codon</c>, <c>Gap weighting is</c>, the numbered matrix with the names) written by
/// <see cref="MultipleAlignmentStatistics.DistanceMatrixFile"/> on 300 random alignments (every -nucmethod / -protmethod / -gapweight / -ambiguous / -position; TestData/Alignment/distmat_file_oracle.jsonl.gz,
/// scripts/alignment/generate_distmat_file_oracle.py). The sign of a <c>nan</c> (<c>-nan</c> of the program's x86 arithmetic) is not part of the comparison: the distances of the port come from .NET's math, whose NaN sign differs (as in the matrix test).
/// </summary>
[TestFixture]
public class DistmatFile_Oracle_Tests
{
    private static string NoNanSign(string text) => text.Replace("-nan", "nan", StringComparison.Ordinal);

    [Test]
    public void DistmatFile_EqualsEmboss()
    {
        const string name = "Seqeron.Genomics.Tests.TestData.Alignment.distmat_file_oracle.jsonl.gz";
        using var stream = System.Reflection.Assembly.GetExecutingAssembly().GetManifestResourceStream(name) ?? throw new InvalidOperationException($"Embedded oracle '{name}' not found.");
        using var gz = new GZipStream(stream, CompressionMode.Decompress);
        using var reader = new StreamReader(gz);
        int rows = 0;
        var failures = new List<string>();
        string? line;
        while ((line = reader.ReadLine()) is not null)
        {
            var j = JsonDocument.Parse(line).RootElement;
            var seqs = j.GetProperty("rows").EnumerateArray().Select(x => x.GetString()!).ToList();
            var names = j.GetProperty("names").EnumerateArray().Select(x => x.GetString()!).ToList();
            // EMBOSS decides the type from the first row (the generator gave both menus the method, protein methods are 0 - 2)
            bool nucleotide = Seqeron.Genomics.Infrastructure.EmbossSequenceReader.IsNucleotide(Seqeron.Genomics.Infrastructure.EmbossSequenceReader.CleanGapped(seqs[0]));
            int method = j.GetProperty("method").GetInt32();
            if (!nucleotide && method > 2)
                method = 0;
            var options = new MsaDistanceOptions
            {
                Method = (MsaDistanceMethod)method,
                GapWeight = j.GetProperty("gapweight").GetDouble(),
                Ambiguous = j.GetProperty("ambiguous").GetBoolean(),
                Position = j.GetProperty("position").GetInt32(),
                CalculateA = j.GetProperty("calculatea").GetBoolean(),
                ParameterA = j.GetProperty("parametera").GetDouble(),
            };
            string expected = j.GetProperty("report").GetString()!;
            string got = MultipleAlignmentStatistics.DistanceMatrixFile(names, MultipleAlignmentStatistics.CalculateDistanceMatrix(seqs, options), method, options, nucleotide);
            if (NoNanSign(got) != NoNanSign(expected))
                failures.Add($"row {rows}\n--- expected\n{expected}--- actual\n{got}");
            rows++;
        }

        Assert.That(rows, Is.GreaterThanOrEqualTo(300));
        Assert.That(failures, Is.Empty, $"{failures.Count} of {rows}:\n" + string.Join("\n", failures.Take(2)));
    }

    [Test]
    public void DistmatFile_RefusesAPositionTheProgramRefuses()
    {
        Assert.Throws<ArgumentException>(() => MultipleAlignmentStatistics.DistanceMatrixFile(["a", "b"], new double[2, 2], 0, new MsaDistanceOptions { Position = 23 }, true));
        Assert.Throws<ArgumentException>(() => MultipleAlignmentStatistics.DistanceMatrixFile(["a", "b", "c"], new double[2, 2], 0, new MsaDistanceOptions(), true));
    }
}
