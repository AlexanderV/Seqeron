namespace Seqeron.Genomics.Tests.Fuzzing;

/// <summary>
/// Fuzz tests for the B05 follow-up consensus algorithms: arbitrary printable / non-ASCII characters,
/// extreme weights and thresholds, empty and long alignments. Every call must either return a string of
/// the alignment length or throw a documented ArgumentException subtype — never another exception.
/// </summary>
[TestFixture]
[Category("Fuzzing")]
public class AlignmentConsensusFuzzTests
{
    [Test]
    public void AllConsensusMethods_RandomInput_NoUndocumentedException()
    {
        var rng = new Random(7701);
        const string pool = "ACGTUNRYacgtu-.~*?XBZJOxbz 1#é中\ud83d";
        for (int trial = 0; trial < 1500; trial++)
        {
            int n = rng.Next(0, 8), length = rng.Next(0, 20);
            var rows = Enumerable.Range(0, n).Select(_ => new string(Enumerable.Range(0, length)
                .Select(_ => pool[rng.Next(pool.Length)]).ToArray())).ToArray();
            float[]? weights = rng.Next(3) == 0 ? rows.Select(_ => (float)(rng.NextDouble() * 1e6)).ToArray() : null;
            float? plurality = rng.Next(4) switch { 0 => null, 1 => float.PositiveInfinity, 2 => -1e30f, _ => (float)rng.NextDouble() * 10 };

            Check(() => MotifFinder.GenerateEmbossConsensus(rows, (ConsensusResidueType)rng.Next(0, 2), plurality, rng.Next(0, 10), null, weights), length);
            Check(() => MotifFinder.GenerateDumbConsensus(rows, rng.NextDouble() * 1.5 - 0.25, '#', rng.Next(2) == 0), n == 0 ? 0 : length);
            Check(() => MotifFinder.GenerateConsensus(rows, rng.NextDouble()), n == 0 ? 0 : length);
        }
    }

    [Test]
    public void Emboss_LongAlignment_Completes()
    {
        var rng = new Random(7702);
        var rows = Enumerable.Range(0, 40).Select(_ => new string(Enumerable.Range(0, 5000)
            .Select(_ => "ACGT-"[rng.Next(5)]).ToArray())).ToArray();
        Assert.That(MotifFinder.GenerateEmbossConsensus(rows).Length, Is.EqualTo(5000));
    }

    private static void Check(Func<string> call, int expectedLength)
    {
        try
        {
            Assert.That(call().Length, Is.EqualTo(expectedLength));
        }
        catch (ArgumentException)
        {
            // Documented: invalid characters, < 2 rows (EMBOSS), weight count mismatch, NaN thresholds.
        }
    }
}
