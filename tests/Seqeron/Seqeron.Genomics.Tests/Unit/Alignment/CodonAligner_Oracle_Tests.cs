using System.IO.Compression;
using System.Reflection;
using System.Text.Json;
using NUnit.Framework;
using Seqeron.Genomics.Alignment;

namespace Seqeron.Genomics.Tests.Unit.Alignment;

/// <summary>
/// B09 WP67 (F307, R14-10b): <see cref="CodonAligner"/> equals Bio.Align.CodonAligner 1.88 (<c>_codonaligner.c</c>) on the rows of TestData/Alignment/codon_aligner_oracle.jsonl.gz (scripts/alignment/generate_codon_aligner_oracle.py): the
/// score, the score of the alignments, and every alignment (in the order of the iterator, as coordinates with the frame-shift vertices). Rows on which the reference reads memory its traceback never wrote (an optimal alignment that skips leading
/// residues: the iterator then lists one-vertex alignments for ever; a heap overflow in <c>PathGenerator_create_path</c>) have no defined result: the port refuses them with <see cref="InvalidOperationException"/>. A codon with a letter
/// outside IUPAC is the <c>TranslationError</c> of the reference (<see cref="ArgumentException"/> here).
/// </summary>
[TestFixture]
public class CodonAligner_Oracle_Tests
{
    private static string Join(IEnumerable<IEnumerable<long>> coordinates) => string.Join("|", coordinates.Select(r => string.Join(",", r)));

    private static bool Throws<T>(Action action) where T : Exception
    {
        try
        {
            action();
            return false;
        }
        catch (T)
        {
            return true;
        }
    }

    [Test]
    public void Aligner_EqualsBiopython_AllRows()
    {
        const string resource = "Seqeron.Genomics.Tests.TestData.Alignment.codon_aligner_oracle.jsonl.gz";
        using var stream = Assembly.GetExecutingAssembly().GetManifestResourceStream(resource) ?? throw new InvalidOperationException($"Embedded oracle '{resource}' not found.");
        using var gz = new GZipStream(stream, CompressionMode.Decompress);
        using var reader = new StreamReader(gz);
        var failures = new List<string>();
        int rows = 0, refused = 0, errors = 0, multi = 0;
        string? line;
        while ((line = reader.ReadLine()) is not null)
        {
            var j = JsonDocument.Parse(line).RootElement;
            rows++;
            var aligner = new CodonAligner(j.GetProperty("table").GetInt32())
            {
                MatchScore = j.GetProperty("match").GetDouble(),
                MismatchScore = j.GetProperty("mismatch").GetDouble(),
                Wildcard = j.GetProperty("wildcard").GetString()![0],
            };
            var fs = j.GetProperty("frameshift").EnumerateArray().Select(x => x.GetDouble()).ToArray();
            aligner.FrameshiftMinusTwoScore = fs[0];
            aligner.FrameshiftMinusOneScore = fs[1];
            aligner.FrameshiftPlusOneScore = fs[2];
            aligner.FrameshiftPlusTwoScore = fs[3];
            string protein = j.GetProperty("protein").GetString()!, dna = j.GetProperty("dna").GetString()!;
            string id = $"row {j.GetProperty("id")} {protein} {dna} table {j.GetProperty("table")}";
            if (j.TryGetProperty("error", out _))
            {
                errors++;
                if (!Throws<ArgumentException>(() => aligner.Score(protein, dna)))
                    failures.Add($"{id}: no error");
                continue;
            }

            double score = aligner.Score(protein, dna);
            if (Math.Abs(score - j.GetProperty("score").GetDouble()) > 1e-9)
            {
                failures.Add($"{id}: score {score} expected {j.GetProperty("score")}");
                continue;
            }

            if (j.GetProperty("undefined").GetBoolean())
            {
                refused++;
                if (!Throws<InvalidOperationException>(() => aligner.Align(protein, dna)))
                    failures.Add($"{id}: not refused");
                continue;
            }

            var alignments = aligner.Align(protein, dna);
            if (Math.Abs(alignments.Score - j.GetProperty("aligned_score").GetDouble()) > 1e-9)
                failures.Add($"{id}: aligned score {alignments.Score} expected {j.GetProperty("aligned_score")}");
            // len() of the reference is a count of paths; its iterator can list fewer (no residues) or more (it cycles): the list is the first min(total, listed) alignments of the iterator
            var want = j.GetProperty("alignments").EnumerateArray().Select(x => Join(x.EnumerateArray().Select(r => r.EnumerateArray().Select(v => v.GetInt64())))).ToList();
            int expectedCount = Math.Min(Math.Min(j.GetProperty("total").GetInt32(), want.Count), 1000);
            if (alignments.Count != expectedCount)
            {
                failures.Add($"{id}: {alignments.Count} alignments expected {expectedCount}");
                continue;
            }

            if (alignments.Count > 1)
                multi++;
            for (int k = 0; k < Math.Min(expectedCount, want.Count); k++)
            {
                string got = Join(alignments[k].Coordinates.Select(r => r.Select(v => v)));
                if (got != want[k])
                {
                    failures.Add($"{id}: alignment {k} {got} expected {want[k]}");
                    break;
                }
            }
        }

        Assert.That(failures, Is.Empty, string.Join(Environment.NewLine, failures.Take(10)));
        Assert.Multiple(() =>
        {
            Assert.That(rows, Is.GreaterThanOrEqualTo(1500));
            Assert.That(multi, Is.GreaterThan(300), "rows with several alignments");
            Assert.That(refused, Is.InRange(1, rows / 10), "rows without a defined result in the reference");
            Assert.That(errors, Is.EqualTo(1));
        });
    }

    [Test]
    public void Properties_AggregateFrameshiftScores_AsBiopython()
    {
        var aligner = new CodonAligner();
        Assert.Multiple(() =>
        {
            Assert.That(aligner.MatchScore, Is.EqualTo(1.0));
            Assert.That(aligner.MismatchScore, Is.EqualTo(0.0));
            Assert.That(aligner.Wildcard, Is.EqualTo('X'));
            Assert.That(aligner.Epsilon, Is.EqualTo(1e-6));
            Assert.That(aligner.FrameshiftScore, Is.EqualTo(-3.0));
        });
        aligner.FrameshiftPlusScore = -1.0;
        Assert.Multiple(() =>
        {
            Assert.That(aligner.FrameshiftPlusOneScore, Is.EqualTo(-1.0));
            Assert.That(aligner.FrameshiftPlusTwoScore, Is.EqualTo(-1.0));
            Assert.That(aligner.FrameshiftMinusScore, Is.EqualTo(-3.0));
            Assert.Throws<InvalidOperationException>(() => _ = aligner.FrameshiftScore);
            Assert.Throws<InvalidOperationException>(() => _ = aligner.FrameshiftOneScore);
        });
        aligner.FrameshiftScore = -2.0;
        Assert.That(aligner.FrameshiftTwoScore, Is.EqualTo(-2.0));
    }

    [Test]
    public void Align_FrameshiftExample_FromBiopythonDocumentation()
    {
        // Bio.Align.CodonAligner doctest: a frame shift between codon 7 and 8
        var aligner = new CodonAligner();
        var alignments = aligner.Align("MLGSNESRVCPKLSPS", "ATGCTGGGCTCGAACGAGTCCGTGTATGCCCTAAGCTGAGCCCGTCG");
        Assert.Multiple(() =>
        {
            Assert.That(alignments.Score, Is.EqualTo(13.0));
            Assert.That(alignments, Has.Count.EqualTo(1));
            Assert.That(Join(alignments[0].Coordinates.Select(r => r.Select(v => v))), Is.EqualTo("0,7,7,16|0,21,20,47"));
        });
    }
}
