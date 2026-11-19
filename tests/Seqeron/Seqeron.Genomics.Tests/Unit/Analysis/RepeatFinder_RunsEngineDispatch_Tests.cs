namespace Seqeron.Genomics.Tests.Unit.Analysis;

/// <summary>
/// R15 (B04 request registered by B09 GENOMIC-TANDEM-001 F36, done as "B09 for B04"): <see cref="RepeatFinder.FindMicrosatellites(DnaSequence,int,int,int)"/>
/// (and the string / threshold-map / cancellable overloads) answer a wide unit-length range (more than
/// <see cref="RepeatFinder.RunsEngineMinUnitLengths"/> admissible unit lengths) with the run enumeration (Bannai et al. 2017) instead of one
/// linear scan per unit length. The result and its order must be identical to the scan.
/// </summary>
[TestFixture]
public class RepeatFinder_RunsEngineDispatch_Tests
{
    private sealed class SyncProgress(List<double> sink) : IProgress<double>
    {
        public void Report(double value) => sink.Add(value);
    }

    private static string RandomDna(Random rnd, int n, string alphabet = "ACGT") =>
        new(Enumerable.Range(0, n).Select(_ => alphabet[rnd.Next(alphabet.Length)]).ToArray());

    private static List<(string, int, int, int, RepeatType)> Key(IEnumerable<MicrosatelliteResult> hits) =>
        hits.Select(h => (h.RepeatUnit, h.Position, h.RepeatCount, h.TotalLength, h.RepeatType)).ToList();

    private static string Planted(Random rnd, int length)
    {
        var sb = new System.Text.StringBuilder();
        while (sb.Length < length)
        {
            if (rnd.Next(3) == 0)
                sb.Append(RandomDna(rnd, rnd.Next(1, 30)));
            else
            {
                string unit = RandomDna(rnd, rnd.Next(1, 60), rnd.Next(2) == 0 ? "AC" : "ACGT");
                for (int c = rnd.Next(2, 6); c > 0; c--) sb.Append(unit);
            }
        }
        return sb.ToString();
    }

    [Test]
    public void WideRange_AutoEqualsTheScan_IncludingMetadata_OnPlantedAndRandomSequences()
    {
        var rnd = new Random(1511);
        for (int iter = 0; iter < 30; iter++)
        {
            string seq = Planted(rnd, rnd.Next(1600, 2400));
            int minReps = rnd.Next(2, 4);
            int minUnit = rnd.Next(1, 4);
            // n / minReps > 256 + minUnit unit lengths, so Auto selects the run enumeration.
            Assert.That(seq.Length / minReps - minUnit + 1, Is.GreaterThan(RepeatFinder.RunsEngineMinUnitLengths), $"iteration {iter}");

            var auto = Key(RepeatFinder.FindMicrosatellites(seq, minUnit, int.MaxValue, minReps));
            var scan = Key(RepeatFinder.FindMicrosatellitesWithEngine(seq, minUnit, int.MaxValue, minReps, RepeatFinder.MicrosatelliteEngine.Scan));
            var runs = Key(RepeatFinder.FindMicrosatellitesWithEngine(seq, minUnit, int.MaxValue, minReps, RepeatFinder.MicrosatelliteEngine.Runs));
            Assert.That(auto, Is.EqualTo(scan), $"iteration {iter} auto vs scan");
            Assert.That(runs, Is.EqualTo(scan), $"iteration {iter} runs vs scan");
            Assert.That(auto, Is.Not.Empty);
        }
    }

    [Test]
    public void WideRange_NonAcgtSymbolsAndLowerCase_AreHandledLikeTheScan()
    {
        var rnd = new Random(8);
        string seq = Planted(rnd, 900).ToLowerInvariant();
        seq = seq[..200] + "NNNNNNNNNNNNNNNNNNNN" + seq[200..400] + "RYACGTACGTACGT" + seq[400..];
        var auto = Key(RepeatFinder.FindMicrosatellites(seq, 1, int.MaxValue, 2));
        var scan = Key(RepeatFinder.FindMicrosatellitesWithEngine(seq, 1, int.MaxValue, 2, RepeatFinder.MicrosatelliteEngine.Scan));
        Assert.That(auto, Is.EqualTo(scan));
        Assert.That(auto.Any(h => h.Item1.Contains('N')), Is.False);
    }

    [Test]
    public void WideThresholdMap_WithPerUnitThresholds_EqualsTheScan()
    {
        var rnd = new Random(42);
        string seq = Planted(rnd, 1400);
        // 300 unit lengths, per-unit thresholds 2..4: more than RunsEngineMinUnitLengths entries -> run enumeration.
        var map = Enumerable.Range(1, 300).ToDictionary(u => u, u => 2 + u % 3);
        var auto = Key(RepeatFinder.FindMicrosatellites(new DnaSequence(seq), map, CancellationToken.None, null));

        // the same map restricted to 256 entries uses the scan; chunk the map into two scans and merge by (unit, start)
        var low = map.Where(kv => kv.Key <= 200).ToDictionary(kv => kv.Key, kv => kv.Value);
        var high = map.Where(kv => kv.Key > 200).ToDictionary(kv => kv.Key, kv => kv.Value);
        var parts = Key(RepeatFinder.FindMicrosatellites(new DnaSequence(seq), low, CancellationToken.None, null))
            .Concat(Key(RepeatFinder.FindMicrosatellites(new DnaSequence(seq), high, CancellationToken.None, null))).ToList();
        Assert.That(auto, Is.EqualTo(parts));
        Assert.That(auto, Is.Not.Empty);
    }

    [Test]
    public void WideRange_HonoursCancellation_AndReportsFinalProgress()
    {
        var rnd = new Random(3);
        string seq = Planted(rnd, 1200);
        using var cts = new CancellationTokenSource();
        cts.Cancel();
        Assert.Throws<OperationCanceledException>(() =>
            RepeatFinder.FindMicrosatellites(new DnaSequence(seq), 1, int.MaxValue, 2, cts.Token).ToList());

        var values = new List<double>();
        var hits = RepeatFinder.FindMicrosatellites(new DnaSequence(seq), 1, int.MaxValue, 2, CancellationToken.None, new SyncProgress(values)).ToList();
        Assert.That(hits, Is.Not.Empty);
        Assert.That(values, Is.Not.Empty);
        Assert.That(values.Max(), Is.EqualTo(1.0));
        Assert.That(values, Is.Ordered);
    }

    [Test]
    public void MisaScanMode_IsNotRedirectedToTheRunsEngine()
    {
        // MisaRegex is a different scan (consumed non-primitive matches); a wide map must still run it.
        string seq = string.Concat(Enumerable.Repeat("ACGTTGCAAG", 80));
        var map = Enumerable.Range(1, 300).ToDictionary(u => u, _ => 2);
        var misa = RepeatFinder.FindMicrosatellites(new DnaSequence(seq), map, MicrosatelliteScanMode.MisaRegex, CancellationToken.None, null).ToList();
        var runs = RepeatFinder.FindMicrosatellites(new DnaSequence(seq), map, MicrosatelliteScanMode.MaximalRuns, CancellationToken.None, null).ToList();
        Assert.That(misa, Is.Not.Empty);
        Assert.That(runs, Is.Not.Empty);
        // MaximalRuns reports the primitive run once at its left end (period 10, 80 copies); MISA greedy match consumes the same run
        Assert.That(runs.Single(h => h.RepeatUnit.Length == 10).RepeatCount, Is.EqualTo(80));
    }

    [Test]
    public void GenomicAnalyzerFindTandemRepeats_IsAPureFacade_OverTheCanonicalEngineChoice()
    {
        var rnd = new Random(99);
        string seq = Planted(rnd, 1300);
        var viaFacade = GenomicAnalyzer.FindTandemRepeats(new DnaSequence(seq), 1, 2).Select(t => (t.Unit, t.Position, t.Repetitions)).ToList();
        var viaCanonical = RepeatFinder.FindMicrosatellites(new DnaSequence(seq), 1, int.MaxValue, 2).Select(m => (m.RepeatUnit, m.Position, m.RepeatCount)).ToList();
        Assert.That(viaFacade, Is.EqualTo(viaCanonical));
        Assert.That(GenomicAnalyzer.FindTandemRepeats(seq, 1, 2).Select(t => (t.Unit, t.Position, t.Repetitions)).ToList(), Is.EqualTo(viaCanonical));
    }
}
