using System.Reflection;
using Seqeron.Genomics.Annotation;

namespace Seqeron.Genomics.Tests.Unit.Annotation;

/// <summary>
/// SPLICE-PREDICT-001 (B11 A3, F45): <see cref="SpliceSitePredictor.PredictGeneStructureDp"/> - splice-site
/// log-odds + open-reading-frame dynamic program (GENSCAN / GeneSplicer style). Oracle 1: exhaustive search
/// over all intron sets (independent implementation, spliced-sequence based). Oracle 2: Ensembl GRCh38
/// annotated coding introns of HBB, TP53 and CFTR exons 2-5.
/// </summary>
[TestFixture]
public class SpliceSitePredictor_GeneStructureDp_Tests
{
    // ---------- Oracle 1: exhaustive search ----------

    private static readonly bool[,,] Stops = BuildStops();

    private static bool[,,] BuildStops()
    {
        var t = new bool[128, 128, 128];
        t['U', 'A', 'A'] = t['U', 'A', 'G'] = t['U', 'G', 'A'] = true;
        return t;
    }

    private static bool IsStop(string s, int i) => Stops[s[i], s[i + 1], s[i + 2]];

    private static double? CodingScorer(string exon, int offset)
    {
        // Deterministic pseudo score of (exon, offset) in [-3, 3).
        long h = offset * 7919L;
        foreach (char c in exon) h = (h * 31 + c) % 1_000_003;
        return (h % 600) / 100.0 - 3.0;
    }

    private sealed record Brute(double Score, List<(int D, int A)> Introns, double Second);

    private static Brute? BruteForce(string u, int start, List<(int, double)> donors, List<(int, double)> acceptors,
        int minExonLength, int minIntronLength, int maxIntronLength, bool scorer)
    {
        int minE = Math.Max(2, minExonLength), minI = Math.Max(4, minIntronLength);
        var pairs = new List<(int D, int A, double S)>();
        foreach (var (d, ds) in donors)
            foreach (var (a, ac) in acceptors)
                if (a - d + 1 >= minI && a - d + 1 <= maxIntronLength)
                    pairs.Add((d, a, ds + ac));
        pairs.Sort((x, y) => x.D != y.D ? x.D.CompareTo(y.D) : x.A.CompareTo(y.A));

        double best = double.NegativeInfinity, second = double.NegativeInfinity;
        List<(int D, int A)>? bestSet = null;
        var current = new List<(int D, int A, double S)>();

        void Consider(double score)
        {
            if (score > best + 1e-12) { second = best; best = score; bestSet = current.Select(p => (p.D, p.A)).ToList(); }
            else if (score > second + 1e-12 && score < best - 1e-12) second = score;
            else if (Math.Abs(score - best) <= 1e-12) second = best; // tie: optimum not unique
        }

        void Evaluate()
        {
            int n = current.Count;
            if (n == 0) { if (u.Length - start < 3) return; }
            else
            {
                if (current[0].D - start < Math.Max(3, minE)) return;
                if (u.Length - 1 - current[n - 1].A < minE) return;
            }

            var sb = new System.Text.StringBuilder();
            int pos = 0;
            var junctions = new List<int>();
            foreach (var p in current)
            {
                sb.Append(u, pos, p.D - pos);
                junctions.Add(sb.Length);
                pos = p.A + 1;
            }
            sb.Append(u, pos, u.Length - pos);
            string sp = sb.ToString();
            int t = -1;
            for (int i = start; i + 3 <= sp.Length; i += 3)
                if (IsStop(sp, i)) { t = i; break; }
            if (t < 0) return;
            if (junctions.Any(j => j > t + 2)) return;

            double score = current.Sum(p => p.S);
            if (scorer)
            {
                int coding = 0;
                int from = start;
                for (int k = 0; k <= n; k++)
                {
                    int toExcl = k < n
                        ? current[k].D
                        : (n == 0 ? t + 3 : current[n - 1].A + 1 + (t + 2 - junctions[n - 1]) + 1);
                    int offset = (3 - coding % 3) % 3;
                    score += CodingScorer(u.Substring(from, toExcl - from).Replace('U', 'T'), offset)!.Value;
                    coding += toExcl - from;
                    if (k < n) from = current[k].A + 1;
                }
            }
            Consider(score);
        }

        void Recurse(int next)
        {
            Evaluate();
            for (int i = next; i < pairs.Count; i++)
            {
                var p = pairs[i];
                if (current.Count > 0)
                {
                    var last = current[^1];
                    if (p.D <= last.A || p.D - last.A - 1 < minE) continue;
                }
                current.Add(p);
                Recurse(i + 1);
                current.RemoveAt(current.Count - 1);
            }
        }

        Recurse(0);
        return bestSet is null ? null : new Brute(best, bestSet, second);
    }

    private static MethodInfo Solve => typeof(SpliceSitePredictor).GetMethod(
        "SolveSpliceDp", BindingFlags.NonPublic | BindingFlags.Static)!;

    private static (double Score, List<(int, int)> Introns)? RunSolver(string u, int start,
        List<(int, double)> donors, List<(int, double)> acceptors, int minE, int minI, int maxI, bool scorer)
    {
        Func<string, int, double>? f = scorer ? (s, o) => CodingScorer(s, o)!.Value : null;
        var r = Solve.Invoke(null, new object?[] { u, start,
            donors.Select(d => (d.Item1, d.Item2)).ToList(), acceptors.Select(a => (a.Item1, a.Item2)).ToList(), minE, minI, maxI, f });
        if (r is null) return null;
        var t = r.GetType();
        double score = (double)t.GetProperty("Score")!.GetValue(r)!;
        var list = (List<(int Donor, int Acceptor)>)t.GetProperty("Introns")!.GetValue(r)!;
        return (score, list.Select(p => (p.Donor, p.Acceptor)).ToList());
    }

    [Test]
    public void Dp_EqualsExhaustiveSearch_RandomSequences()
    {
        var rng = new Random(20261125);
        int cases = 0, withIntrons = 0, nullCases = 0, uniqueChecked = 0;
        for (int iter = 0; iter < 4000 && cases < 1200; iter++)
        {
            int len = rng.Next(40, 121);
            var chars = new char[len];
            for (int i = 0; i < len; i++)
            {
                double r = rng.NextDouble();
                chars[i] = r < 0.34 ? 'A' : r < 0.68 ? 'C' : r < 0.84 ? 'G' : 'U';
            }
            // Plant a terminal stop region so that open frames are reachable.
            int tail = rng.Next(0, 3) == 0 ? -1 : len - 3 - rng.Next(0, 6);
            if (tail > 0) { chars[tail] = 'U'; chars[tail + 1] = rng.Next(2) == 0 ? 'A' : 'G'; chars[tail + 2] = 'A'; }
            string u = new(chars);
            int start = rng.Next(0, 8);
            int nd = rng.Next(1, 6), na = rng.Next(1, 6);
            var donors = Enumerable.Range(0, nd).Select(_ => rng.Next(start + 3, len - 4)).Distinct().OrderBy(x => x)
                .Select(p => (p, Math.Round(rng.NextDouble() * 12 - 3, 3))).ToList();
            var acceptors = Enumerable.Range(0, na).Select(_ => rng.Next(start + 6, len - 2)).Distinct().OrderBy(x => x)
                .Select(p => (p, Math.Round(rng.NextDouble() * 14 - 4, 3))).ToList();
            int minE = rng.Next(1, 9), minI = rng.Next(1, 20), maxI = rng.Next(minI, 90);
            bool scorer = rng.Next(3) == 0;

            var brute = BruteForce(u, start, donors, acceptors, minE, minI, maxI, scorer);
            var dp = RunSolver(u, start, donors, acceptors, minE, minI, maxI, scorer);
            cases++;
            if (brute is null)
            {
                Assert.That(dp, Is.Null, $"case {iter}: DP found a structure but exhaustive search none");
                nullCases++;
                continue;
            }
            Assert.That(dp, Is.Not.Null, $"case {iter}: exhaustive search found a structure, DP none");
            Assert.That(dp!.Value.Score, Is.EqualTo(brute.Score).Within(1e-9), $"case {iter}: score");
            if (brute.Introns.Count > 0) withIntrons++;
            if (brute.Score - brute.Second > 1e-9)
            {
                Assert.That(dp.Value.Introns, Is.EqualTo(brute.Introns), $"case {iter}: structure");
                uniqueChecked++;
            }
        }
        TestContext.Out.WriteLine($"cases={cases} optimum-with-introns={withIntrons} no-ORF={nullCases} structure-compared={uniqueChecked}");
        Assert.That(cases, Is.GreaterThanOrEqualTo(1000));
        Assert.That(withIntrons, Is.GreaterThan(200));
    }

    // ---------- Oracle 2: Ensembl annotated genes ----------

    private sealed record Gene(string Name, int Start, List<(int S, int E)> Introns, string Seq);

    private static List<Gene> LoadGenes()
    {
        var asm = Assembly.GetExecutingAssembly();
        string res = asm.GetManifestResourceNames().Single(n => n.EndsWith("ensembl_genes.tsv"));
        using var reader = new StreamReader(asm.GetManifestResourceStream(res)!);
        var genes = new List<Gene>();
        string? line;
        while ((line = reader.ReadLine()) is not null)
        {
            if (line.StartsWith('#') || line.Length == 0) continue;
            var f = line.Split('\t');
            var introns = f[2].Split(',', StringSplitOptions.RemoveEmptyEntries)
                .Select(x => x.Split('-')).Select(x => (int.Parse(x[0]), int.Parse(x[1]))).ToList();
            genes.Add(new Gene(f[0], int.Parse(f[1]), introns, f[3]));
        }
        return genes;
    }

    private static readonly Lazy<(IReadOnlyDictionary<string, double> Coding, IReadOnlyDictionary<string, double> Noncoding)> Hexamers = new(() =>
    {
        var asm = Assembly.GetExecutingAssembly();
        string res = asm.GetManifestResourceNames().Single(n => n.EndsWith("hexamer_train12.tsv"));
        using var reader = new StreamReader(asm.GetManifestResourceStream(res)!);
        var lines = new List<string>();
        string? line;
        while ((line = reader.ReadLine()) is not null) lines.Add(line);
        return GenomeAnnotator.CpatModel.ParseHexamerTable(lines);
    });

    // Sum of the CPAT in-frame hexamer log-likelihood ratios of the exon, in bits.
    private static double HexamerExonScore(string exon, int offset)
    {
        string e = exon.Substring(offset);
        if (e.Length < 6) return 0;
        int n = (e.Length - 6) / 3 + 1;
        return n * GenomeAnnotator.CalculateCodingPotential(e, Hexamers.Value.Coding, Hexamers.Value.Noncoding) / Math.Log(2);
    }

    private static (int Hit, int Predicted, int DonorOk, int AcceptorOk) Evaluate(Gene g, Func<string, int, double>? scorer)
    {
        var r = SpliceSitePredictor.PredictGeneStructureDp(g.Seq, g.Start, exonCodingScore: scorer);
        Assert.That(r, Is.Not.Null, g.Name);
        var pred = r!.Value.Structure.Introns.Select(i => (i.Start, i.End)).ToHashSet();
        return (g.Introns.Count(i => pred.Contains(i)), pred.Count,
            g.Introns.Count(i => pred.Any(q => q.Item1 == i.S)), g.Introns.Count(i => pred.Any(q => q.Item2 == i.E)));
    }

    [Test]
    public void Dp_EnsemblGenes_SpliceAndOrfOnly_Recall()
    {
        var res = LoadGenes().ToDictionary(g => g.Name, g => Evaluate(g, null));
        foreach (var kv in res) TestContext.Out.WriteLine($"{kv.Key} (ORF only): hit={kv.Value.Hit} predicted={kv.Value.Predicted} donorOk={kv.Value.DonorOk} acceptorOk={kv.Value.AcceptorOk}");
        // Greedy PredictGeneStructure recovers 0 of these 14 introns (B11 F32-F38 measurement).
        Assert.That(res["HBB"].Hit, Is.EqualTo(2));
        Assert.That(res["TP53"].Hit, Is.GreaterThanOrEqualTo(5));
        Assert.That(res.Values.Sum(v => v.Hit), Is.GreaterThanOrEqualTo(8));
    }

    [Test]
    public void Dp_EnsemblGenes_WithHexamerCodingScore_Recall()
    {
        var res = LoadGenes().ToDictionary(g => g.Name, g => Evaluate(g, HexamerExonScore));
        foreach (var kv in res) TestContext.Out.WriteLine($"{kv.Key} (+hexamer): hit={kv.Value.Hit} predicted={kv.Value.Predicted} donorOk={kv.Value.DonorOk} acceptorOk={kv.Value.AcceptorOk}");
        Assert.That(res["HBB"].Hit, Is.EqualTo(2));
        Assert.That(res["TP53"].Hit, Is.EqualTo(8));
        Assert.That(res.Values.Sum(v => v.Hit), Is.GreaterThanOrEqualTo(10));
    }

    [Test]
    public void Dp_Hbb_CodingSequenceIsTheHemoglobinBetaOrf()
    {
        var g = LoadGenes().Single(x => x.Name == "HBB");
        var r = SpliceSitePredictor.PredictGeneStructureDp(g.Seq)!.Value; // first AUG, no start given
        Assert.That(r.CdsStart, Is.EqualTo(g.Start));
        Assert.That(r.Structure.Introns.Select(i => (i.Start, i.End)), Is.EquivalentTo(g.Introns));
        Assert.That(r.CodingSequence.Length, Is.EqualTo(3 * 147 + 3)); // 147 aa + stop (UniProt P68871)
        string protein = Seqeron.Genomics.Core.Translator.Translate(new Seqeron.Genomics.Core.DnaSequence(r.CodingSequence.Replace('U', 'T'))).Sequence;
        Assert.That(protein, Does.StartWith("MVHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAVMGNPKVKAHGKKVLGAFSDGLAHLDNLKGTFATLSELHCDKLHVDPENFR"));
        Assert.That(protein, Does.EndWith("AHKYH*"));
        Assert.That(r.Structure.Exons.Sum(e => e.Length) + r.Structure.Introns.Sum(i => i.Length), Is.EqualTo(g.Seq.Length));
        Assert.That(r.Structure.SplicedSequence.Length, Is.EqualTo(g.Seq.Length - g.Introns.Sum(i => i.E - i.S + 1)));
        Assert.That(r.Structure.Exons.Select(e => e.Phase), Is.EqualTo(new int?[] { 0, (g.Introns[0].S - g.Start) % 3, (g.Introns[0].S - g.Start + g.Introns[1].S - g.Introns[0].E - 1) % 3 }));
    }

    [Test]
    public void Dp_TP53_GeneScale_FastAndBounded()
    {
        var g = LoadGenes().Single(x => x.Name == "TP53");
        var sw = System.Diagnostics.Stopwatch.StartNew();
        var r = SpliceSitePredictor.PredictGeneStructureDp(g.Seq, g.Start)!.Value;
        sw.Stop();
        TestContext.Out.WriteLine($"TP53 {g.Seq.Length} nt: {sw.ElapsedMilliseconds} ms, {r.Structure.Introns.Count} introns");
        Assert.That(sw.Elapsed.TotalSeconds, Is.LessThan(15));
    }

    // ---------- constructed cases and API contract ----------

    [Test]
    public void Dp_NoStartCodon_ReturnsNull() =>
        Assert.That(SpliceSitePredictor.PredictGeneStructureDp("CCCCCCCCCCCCCCCCCCCC"), Is.Null);

    [Test]
    public void Dp_EmptyInput_ReturnsNull() =>
        Assert.That(SpliceSitePredictor.PredictGeneStructureDp(""), Is.Null);

    [Test]
    public void Dp_SingleExonOrf_NoIntrons()
    {
        string seq = "CCAUG" + string.Concat(Enumerable.Repeat("GCC", 40)) + "UAAGG";
        var r = SpliceSitePredictor.PredictGeneStructureDp(seq)!.Value;
        Assert.That(r.Structure.Introns, Is.Empty);
        Assert.That(r.Structure.Exons.Single().Type, Is.EqualTo(SpliceSitePredictor.ExonType.Single));
        Assert.That(r.CdsStart, Is.EqualTo(2));
        Assert.That(r.CodingSequence, Is.EqualTo("AUG" + string.Concat(Enumerable.Repeat("GCC", 40)) + "UAA"));
        Assert.That(r.Structure.OverallScore, Is.EqualTo(0));
    }

    [Test]
    public void Dp_NoOpenFrame_ReturnsNull() =>
        Assert.That(SpliceSitePredictor.PredictGeneStructureDp("CCAUGGCCGCCGCCGCC"), Is.Null);

    [Test]
    public void Dp_StartPosition_OverridesFirstAug()
    {
        string seq = "AUGCCUAAGG" + "CCAUGGCCGCCUAAGG";
        var r = SpliceSitePredictor.PredictGeneStructureDp(seq, startPosition: 12)!.Value;
        Assert.That(r.CdsStart, Is.EqualTo(12));
        Assert.That(r.CodingSequence, Is.EqualTo("AUGGCCGCCUAA"));
    }

    [Test]
    public void Dp_InvalidArguments_Throw()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => SpliceSitePredictor.PredictGeneStructureDp("AUGAAAUAA", minExonLength: 0));
        Assert.Throws<ArgumentOutOfRangeException>(() => SpliceSitePredictor.PredictGeneStructureDp("AUGAAAUAA", minIntronLength: 0));
        Assert.Throws<ArgumentOutOfRangeException>(() => SpliceSitePredictor.PredictGeneStructureDp("AUGAAAUAA", minIntronLength: 60, maxIntronLength: 10));
        Assert.Throws<ArgumentOutOfRangeException>(() => SpliceSitePredictor.PredictGeneStructureDp("AUGAAAUAA", startPosition: 8));
    }

    [Test]
    public void Dp_IntronRemovalIsNeededForTheOrf_HbbIntron1RemovedKeepsFrame()
    {
        // HBB intron 1 (130 nt) is not a multiple of 3 and contains in-frame stops: the unspliced frame closes early.
        var g = LoadGenes().Single(x => x.Name == "HBB");
        string upper = g.Seq.Replace('T', 'U');
        int firstStop = -1;
        for (int i = g.Start; i + 3 <= upper.Length; i += 3)
            if (upper.Substring(i, 3) is "UAA" or "UAG" or "UGA") { firstStop = i; break; }
        Assert.That(firstStop, Is.LessThan(g.Introns[0].E), "unspliced ORF must be closed inside/before intron 1");
        var r = SpliceSitePredictor.PredictGeneStructureDp(g.Seq, g.Start)!.Value;
        Assert.That(r.Structure.Introns.Count, Is.GreaterThanOrEqualTo(1));
        Assert.That(r.CodingSequence.Length % 3, Is.EqualTo(0));
    }

    [Test]
    public void Greedy_PredictGeneStructure_Unchanged_ByTheDpOption()
    {
        var g = LoadGenes().Single(x => x.Name == "HBB");
        var greedy = SpliceSitePredictor.PredictGeneStructure(g.Seq);
        var dp = SpliceSitePredictor.PredictGeneStructureDp(g.Seq, g.Start)!.Value.Structure;
        Assert.That(greedy.Introns.Select(i => (i.Start, i.End)), Is.Not.EqualTo(dp.Introns.Select(i => (i.Start, i.End))));
        Assert.That(SpliceSitePredictor.PredictGeneStructure(g.Seq).OverallScore, Is.EqualTo(greedy.OverallScore));
    }
}
