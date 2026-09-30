using FsCheck;
using FsCheck.Fluent;

namespace Seqeron.Genomics.Tests.Properties;

/// <summary>
/// Heavy-tier properties for the REP-STR-001 / REP-TANDEM-001 MISA additions (review B04, audit WP3): per-unit-size
/// thresholds equal a brute force of maximal primitive runs; compound assembly equals an independent transcription of
/// misa.pl's loop (1-based coordinates as in the Perl source); MISA classes / Krait standard motifs equal brute forces of
/// their definitions; the progress contract holds on random input.
/// Test Units: REP-STR-001, REP-TANDEM-001.
/// </summary>
[TestFixture]
[Category("Property")]
[Category("Repeats")]
public class RepStrMisaProperties
{
    private static Gen<string> Dna(string alphabets, int minLen, int maxLen) =>
        from alpha in Gen.Elements(alphabets.Split('|'))
        from n in Gen.Choose(minLen, maxLen)
        from chars in Gen.Elements(alpha.ToCharArray()).ArrayOf(n)
        select new string(chars);

    private static Gen<Dictionary<int, int>> Thresholds() =>
        from count in Gen.Choose(1, 4)
        from sizes in Gen.Choose(1, 7).ArrayOf(count)
        from k in Gen.Choose(2, 4)
        from bump in Gen.Choose(0, 3)
        select sizes.Distinct().ToDictionary(p => p, p => p == 1 ? k + bump : k);

    private static bool Primitive(string u) =>
        !Enumerable.Range(1, u.Length - 1).Any(j => u.Length % j == 0 && string.Concat(Enumerable.Repeat(u[..j], u.Length / j)) == u);

    private static List<(int, string, int)> BruteRuns(string seq, IReadOnlyDictionary<int, int> th)
    {
        string s = seq.ToUpperInvariant();
        var o = new List<(int, string, int)>();
        foreach (var (p, k) in th.OrderBy(kv => kv.Key))
            for (int i = 0; i + p <= s.Length; i++)
            {
                string u = s.Substring(i, p);
                if (u.Any(c => "ACGT".IndexOf(c) < 0) || !Primitive(u)) continue;
                if (i > 0 && i - 1 + p < s.Length && s[i - 1] == s[i - 1 + p]) continue;
                int e = i + p;
                while (e < s.Length && s[e] == s[e - p]) e++;
                int c = (e - i) / p;
                if (c >= k) o.Add((i, u, c));
            }
        return o;
    }

    /// <summary>Transcription of misa.pl's compound loop (1-based start/end, $space = amb + 1).</summary>
    private static List<(int Start1, int End1, string Type, string Ssr)> MisaLoop(
        string seq, List<MicrosatelliteResult> ssrs, int amb)
    {
        var order = ssrs.OrderBy(m => m.Position).ToList(); // stable
        int nr = order.Count;
        int St(int i) => order[i].Position + 1;
        int En(int i) => order[i].Position + order[i].TotalLength;
        string M(int i) => $"({order[i].RepeatUnit}){order[i].RepeatCount}";
        var rows = new List<(int, int, string, string)>();
        int space = amb + 1;
        int x = 0;
        while (x < nr)
        {
            if (x + 1 >= nr || St(x + 1) - En(x) > space) { x++; continue; }
            int start = St(x);
            string type, ssr;
            if (St(x + 1) - En(x) < 1) { type = "c*"; ssr = M(x) + M(x + 1) + "*"; }
            else { type = "c"; ssr = M(x) + seq.Substring(En(x), St(x + 1) - En(x) - 1).ToLowerInvariant() + M(x + 1); }
            int end = En(x + 1);
            x++;
            while (x + 1 < nr && St(x + 1) - En(x) <= space)
            {
                if (St(x + 1) - En(x) < 1) { ssr += M(x + 1) + "*"; type = "c*"; }
                else ssr += seq.Substring(En(x), St(x + 1) - En(x) - 1).ToLowerInvariant() + M(x + 1);
                end = En(x + 1);
                x++;
            }
            x++;
            rows.Add((start, end, type, ssr));
        }
        return rows;
    }

    [FsCheck.NUnit.Property(MaxTest = 300)]
    public Property PerUnitSizeThresholds_EqualBruteForceMaximalPrimitiveRuns()
    {
        var gen = from s in Dna("ACGT|AT|AC|A|ACGTN|acgt", 0, 90)
                  from th in Thresholds()
                  select (s, th);
        return Prop.ForAll(gen.ToArbitrary(), c =>
        {
            var actual = RepeatFinder.FindMicrosatellites(c.s, c.th).Select(m => (m.Position, m.RepeatUnit, m.RepeatCount)).ToList();
            return actual.SequenceEqual(BruteRuns(c.s, c.th)).Label($"{c.s} {string.Join(",", c.th)}");
        });
    }

    [FsCheck.NUnit.Property(MaxTest = 300)]
    public Property CompoundAssembly_EqualsMisaLoopTranscription()
    {
        var gen = from s in Dna("ACGT|AT|AC|ACGTN|acgt", 0, 120)
                  from k in Gen.Choose(2, 3)
                  from amb in Gen.Elements(0, 1, 3, 10, 100)
                  select (s, k, amb);
        return Prop.ForAll(gen.ToArbitrary(), c =>
        {
            var ssrs = RepeatFinder.FindMicrosatellites(c.s, 1, 6, c.k).ToList();
            var expected = MisaLoop(c.s, ssrs, c.amb);
            var actual = RepeatFinder.FindCompoundMicrosatellites(c.s, c.k, c.amb)
                .Select(x => (x.Start + 1, x.End, x.MisaType, x.Notation)).ToList();
            bool shape = RepeatFinder.FindCompoundMicrosatellites(c.s, c.k, c.amb).All(x =>
                x.Components.Count >= 2 && x.Interruptions.Count == x.Components.Count - 1 && x.Length == x.End - x.Start);
            return (actual.SequenceEqual(expected) && shape).Label($"{c}");
        });
    }

    private static int KraitCompare(string a, string b)
    {
        const string order = "ATCG";
        for (int i = 0; i < a.Length; i++)
            if (a[i] != b[i]) return order.IndexOf(a[i]) - order.IndexOf(b[i]);
        return 0;
    }

    private static IEnumerable<string> Rotations(string m) => Enumerable.Range(0, m.Length).Select(r => m[r..] + m[..r]);

    [FsCheck.NUnit.Property(MaxTest = 500)]
    public Property MotifClasses_EqualDefinitions()
    {
        var gen = Dna("ACGT|AT|CG|ACG", 1, 8);
        return Prop.ForAll(gen.ToArbitrary(), m =>
        {
            string rc = new(m.Reverse().Select(ch => "ACGT"["TGCA".IndexOf(ch)]).ToArray());
            string comp = new(m.Select(ch => "ACGT"["TGCA".IndexOf(ch)]).ToArray());
            string rev = new(m.Reverse().ToArray());
            string f = Rotations(m).Min(StringComparer.Ordinal)!, r = Rotations(rc).Min(StringComparer.Ordinal)!;
            string misa = string.CompareOrdinal(f, r) < 0 ? $"{f}/{r}" : $"{r}/{f}";
            var sets = new[] { new[] { m }, Rotations(m).ToArray(), Rotations(rc).ToArray(), Rotations(comp).ToArray(), Rotations(rev).ToArray() };
            bool krait = Enumerable.Range(0, 5).All(level =>
            {
                var members = level == 0 ? sets[0] : sets.Skip(1).Take(level).SelectMany(x => x).ToArray();
                string best = members.Aggregate((a, b) => KraitCompare(b, a) < 0 ? b : a);
                return RepeatFinder.GetStandardMotif(m, level) == best;
            });
            bool invariant = Rotations(m).Concat(Rotations(rc)).All(x => RepeatFinder.GetCanonicalMotifClass(x) == misa);
            return (RepeatFinder.GetCanonicalMotifClass(m) == misa && krait && invariant).Label(m);
        });
    }

    [FsCheck.NUnit.Property(MaxTest = 100)]
    public Property Progress_MonotoneEndingAtOne_ResultsUnchanged()
    {
        var gen = from s in Dna("ACGT|AT|ACGTN", 0, 6000)
                  from k in Gen.Choose(2, 4)
                  select (s, k);
        return Prop.ForAll(gen.ToArbitrary(), c =>
        {
            var values = new List<double>();
            var progress = new SyncProgress(values.Add);
            var withProgress = RepeatFinder.FindMicrosatellites(c.s, 1, 6, c.k, CancellationToken.None, progress).ToList();
            bool ok = withProgress.SequenceEqual(RepeatFinder.FindMicrosatellites(c.s, 1, 6, c.k))
                      && values.Count >= 1 && values[^1] == 1.0
                      && values.Zip(values.Skip(1)).All(p => p.First <= p.Second)
                      && values.Take(values.Count - 1).All(v => v >= 0 && v < 1);
            return ok.Label($"len {c.s.Length} k {c.k}");
        });
    }

    private sealed class SyncProgress(Action<double> report) : IProgress<double>
    {
        public void Report(double value) => report(value);
    }
}
