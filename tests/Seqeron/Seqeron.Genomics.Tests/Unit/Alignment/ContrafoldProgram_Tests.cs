using System.Globalization;
using NUnit.Framework;
using Seqeron.Genomics.Alignment;

namespace Seqeron.Genomics.Tests.Unit.Alignment;

/// <summary>
/// ALIGN-MULTI-001 (audit A2-19c-rna, B09 WP86, F377): the port of CONTRAfold 2.02 (<c>contrafold predict --posteriors</c>), the program that <c>mafft --contrafold</c> runs. Its output is checked against the real program in two ways:
/// the 16 complete runs of <c>mafft --contrafold</c> of the Docker package (<c>SequenceAligner_MafftRna_Tests</c>; they equal byte for byte, but with the default structure weight 0 the alignment does not depend on the probabilities) and, as soon as run 4 of the package
/// is returned (<c>out/probes.json</c>: 111 sequences run directly through <c>contrafold predict --posteriors 0.01</c>), the byte-for-byte comparison of the posteriors. Here: the format of the output file and the properties that every table of
/// posterior probabilities of the model has (they would break if the inside / outside recursions, the parameter mapping or the float polynomials of <c>LogSpace.hpp</c> were wrong).
/// </summary>
[TestFixture]
[Category("Alignment")]
[Category("ALIGN-MULTI-001")]
public class ContrafoldProgram_Tests
{
    private static List<(int Left, int Right, double P)> Pairs(string text)
    {
        var list = new List<(int, int, double)>();
        foreach (var line in text.Split('\n', StringSplitOptions.RemoveEmptyEntries))
        {
            var f = line.Split(' ');
            int i = int.Parse(f[0], CultureInfo.InvariantCulture);
            foreach (var token in f.Skip(2))
            {
                var kv = token.Split(':');
                list.Add((i, int.Parse(kv[0], CultureInfo.InvariantCulture), double.Parse(kv[1], CultureInfo.InvariantCulture)));
            }
        }
        return list;
    }

    [Test]
    public void Output_HasOneLinePerNucleotide_WithTheLetterAndOnlyLaterPartners()
    {
        const string seq = "GGGGAAAAGCCCCUUUUGGGGAAAAAGCCCC";
        string text = ContrafoldProgram.PredictPosteriors(seq);
        var lines = text.Split('\n', StringSplitOptions.RemoveEmptyEntries);
        Assert.Multiple(() =>
        {
            Assert.That(text, Does.EndWith("\n"));
            Assert.That(lines, Has.Length.EqualTo(seq.Length));
            for (int i = 0; i < lines.Length; i++)
                Assert.That(lines[i], Does.StartWith((i + 1).ToString(CultureInfo.InvariantCulture) + " " + seq[i]));
            Assert.That(Pairs(text), Is.All.Matches<(int Left, int Right, double P)>(p => p.Right > p.Left && p.Right <= seq.Length && p.P >= 0.01 && p.P <= 1.0));
            Assert.That(Pairs(text), Is.Not.Empty);
        });
    }

    [Test]
    public void Posteriors_OfOnePosition_SumToAtMostOne()
    {
        string seq = "GGGCGCAAAGCGCCCUUUGGGAAACCCAAAGGGGAAAUCCCC";
        var pairs = Pairs(ContrafoldProgram.PredictPosteriors(seq, 0.0001));
        var total = new double[seq.Length + 2];
        foreach (var (l, r, p) in pairs)
        {
            total[l] += p;
            total[r] += p;
        }
        Assert.That(total.Max(), Is.LessThanOrEqualTo(1.0 + 1e-3));
    }

    /// <summary>The inside pass and the outside pass compute the same partition coefficient (the log of the sum over all structures of the weights of the model) in two different orders: they agree only when both recursions are right.</summary>
    [TestCase(40, 1)]
    [TestCase(65, 2)]
    [TestCase(90, 3)]
    [TestCase(130, 4)]
    [TestCase(25, 5)]
    public void InsideAndOutside_GiveTheSamePartitionCoefficient(int length, int seed)
    {
        var rng = new Random(seed);
        var sb = new System.Text.StringBuilder();
        for (int i = 0; i < length; i++)
            sb.Append(rng.Next(25) == 0 ? 'N' : "ACGU"[rng.Next(4)]);
        var engine = new ContrafoldProgram.Engine(sb.ToString());
        engine.ComputeInside();
        engine.ComputeOutside();
        engine.ComputePosterior();
        Assert.Multiple(() =>
        {
            Assert.That(engine.OutsideLogPartition, Is.EqualTo(engine.InsideLogPartition).Within(1e-3 * Math.Max(1.0, Math.Abs(engine.InsideLogPartition))));
            for (int i = 1; i <= length; i++)
            {
                float total = 0f;
                for (int j = 1; j <= length; j++)
                    total += i < j ? engine.Posterior(i, j) : j < i ? engine.Posterior(j, i) : 0f;
                Assert.That(total, Is.LessThanOrEqualTo(1.0f + 2e-2f), "position " + i);
            }
        });
    }

    [Test]
    public void ComplementaryStem_IsTheMostProbableStructure()
    {
        // a perfect hairpin: the outermost pair of the stem is among the pairs and its probability is high
        const string seq = "GGGGGCAAAAGCCCCC";
        var pairs = Pairs(ContrafoldProgram.PredictPosteriors(seq));
        var best = pairs.OrderByDescending(p => p.P).First();
        Assert.Multiple(() =>
        {
            Assert.That(best.P, Is.GreaterThan(0.5));
            Assert.That(pairs.Any(p => p.Left == 1 && p.Right == seq.Length), Is.True);
            // only the six allowed pairs (A-U, G-C, G-U): no pair of two equal letters
            Assert.That(pairs, Is.All.Matches<(int Left, int Right, double P)>(p => "AU UA GC CG GU UG".Contains(string.Concat(seq[p.Left - 1], seq[p.Right - 1]), StringComparison.Ordinal)));
        });
    }

    [Test]
    public void UnknownLetters_PairWithNothing_AndTIsUnknown()
    {
        var pairs = Pairs(ContrafoldProgram.PredictPosteriors("GGGGANNNNCCCCCNNNN"));
        Assert.That(pairs, Is.All.Matches<(int Left, int Right, double P)>(p => p.Left is not (6 or 7 or 8 or 9) && p.Right is not (6 or 7 or 8 or 9 or 15 or 16 or 17 or 18)));
        Assert.That(Pairs(ContrafoldProgram.PredictPosteriors("GGGGAAAACCCC")), Is.Not.Empty);
        Assert.That(Pairs(ContrafoldProgram.PredictPosteriors("GGGGAAAATTTT")), Is.Empty, "T is not a letter of the alphabet of the program (ACGU)");
    }

    [Test]
    public void Cutoff_FiltersTheTable()
    {
        const string seq = "GGGCGCAAAGCGCCCUUUGGGAAACCCAAAGGGGAAAUCCCC";
        var loose = Pairs(ContrafoldProgram.PredictPosteriors(seq, 0.01));
        var tight = Pairs(ContrafoldProgram.PredictPosteriors(seq, 0.5));
        Assert.Multiple(() =>
        {
            Assert.That(tight.Count, Is.LessThan(loose.Count));
            Assert.That(tight, Is.All.Matches<(int Left, int Right, double P)>(p => p.P >= 0.5));
        });
    }

    [Test]
    public void Deterministic_AndLowerCaseEqualsUpperCase()
    {
        const string seq = "GGGCGCAAAGCGCCCUUUGGGAAACCC";
        Assert.That(Pairs(ContrafoldProgram.PredictPosteriors(seq.ToLowerInvariant())), Is.EqualTo(Pairs(ContrafoldProgram.PredictPosteriors(seq))));
        string first = ContrafoldProgram.PredictPosteriors(seq);
        string second = ContrafoldProgram.PredictPosteriors(seq);
        Assert.That(second, Is.EqualTo(first));
    }

    [TestCase("")]
    public void ZeroLength_IsRefused(string seq) => Assert.Throws<ArgumentException>(() => ContrafoldProgram.PredictPosteriors(seq));

    [Test]
    public void SingleNucleotide_GivesOneLineWithoutPairs() => Assert.That(ContrafoldProgram.PredictPosteriors("G"), Is.EqualTo("1 G\n"));
}
