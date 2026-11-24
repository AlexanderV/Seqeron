using System.Diagnostics;
using Seqeron.Genomics;
using Seqeron.Genomics.Alignment;
using Seqeron.Genomics.Analysis;
using Seqeron.Genomics.Core;
using Seqeron.Genomics.Infrastructure;

// B09 WP78 benchmark / golden harness (Release): `dotnet run -c Release -- <motif|etandem|similarity|orf|msa|msa2|golden>`; -p:Tree=<repository root> builds against another checkout (build the "before" with -p:DefineConstants=BASE, which leaves out msa2; the "before" numbers of docs/Validation/review-2026-09/B09.md F360 - F362 were taken on commit 12059b65a).
// `golden write <file.tsv>` writes the alignments of the guide-tree multiple aligners on 117 fixed random families; `golden check <file.tsv>` compares them with the file (the file in tests/.../TestData/Alignment was written by commit 12059b65a, before the exact speed-ups of F362).
var what = args.Length > 0 ? args[0] : "all";
var rng = new Random(78);
string Rand(int n, string a = "ACGT") { var c = new char[n]; for (int i = 0; i < n; i++) c[i] = a[rng.Next(a.Length)]; return new string(c); }
string Mut(string s, double r) { var c = s.ToCharArray(); for (int i = 0; i < c.Length; i++) if (rng.NextDouble() < r) c[i] = "ACGT"[rng.Next(4)]; return new string(c); }
void T(string name, Func<object?> f, int reps = 1)
{
    GC.Collect(); long a0 = GC.GetAllocatedBytesForCurrentThread(); var sw = Stopwatch.StartNew(); object? r = null;
    for (int i = 0; i < reps; i++) r = f();
    sw.Stop(); long a1 = GC.GetAllocatedBytesForCurrentThread();
    Console.WriteLine($"{name,-55} {sw.Elapsed.TotalSeconds / reps,9:F3} s  {(a1 - a0) / reps / 1048576.0,10:F1} MB  {r}");
}
if (what is "all" or "motif")
{
    var seq = new DnaSequence(Rand(2_000_000));
    var motifs = Enumerable.Range(0, 50).Select(_ => Rand(10)).ToList();
    var iupac = Enumerable.Range(0, 50).Select(_ => { var m = Rand(10).ToCharArray(); m[2] = 'N'; m[5] = 'R'; return new string(m); }).ToList();
    T("FindKnownMotifHits 50x10bp 2Mbp both mm=2", () => GenomicAnalyzer.FindKnownMotifHits(seq, motifs, true, false, 2).Sum(kv => kv.Value.Count));
    T("FindKnownMotifHits 50x10bp 2Mbp both mm=1", () => GenomicAnalyzer.FindKnownMotifHits(seq, motifs, true, false, 1).Sum(kv => kv.Value.Count));
    T("FindKnownMotifHits 50x10bp 2Mbp both mm=0", () => GenomicAnalyzer.FindKnownMotifHits(seq, motifs, true, false, 0).Sum(kv => kv.Value.Count));
    T("FindKnownMotifHits IUPAC mm=0", () => GenomicAnalyzer.FindKnownMotifHits(seq, iupac, true, true, 0).Sum(kv => kv.Value.Count));
    T("FindKnownMotifHits IUPAC mm=2", () => GenomicAnalyzer.FindKnownMotifHits(seq, iupac, true, true, 2).Sum(kv => kv.Value.Count));
}
if (what is "all" or "etandem")
{
    string s = Rand(2_000_000);
    T("FindTandemRepeatsEtandem 2Mbp n=10", () => GenomicAnalyzer.FindTandemRepeatsEtandem(s, 10, 10).Count);
    T("FindTandemRepeatsEtandem 2Mbp n=4..8", () => GenomicAnalyzer.FindTandemRepeatsEtandem(s, 4, 8).Count);
}
if (what is "all" or "similarity")
{
    var a = new DnaSequence(Rand(2_000_000)); var b = new DnaSequence(Mut(a.Sequence, 0.05));
    T("CalculateSimilarity 2Mbp k=21", () => GenomicAnalyzer.CalculateSimilarity(a, b, 21));
    T("CalculateSimilarity 2Mbp k=21 canonical", () => GenomicAnalyzer.CalculateSimilarity(a, b, 21, true));
    T("CalculateSimilarity 2Mbp k=12", () => GenomicAnalyzer.CalculateSimilarity(a, b, 12));
    T("CalculateSimilarity strings 2Mbp k=21", () => GenomicAnalyzer.CalculateSimilarity(a.Sequence, b.Sequence, 21));
}
if (what is "all" or "orf")
{
    var nested = new DnaSequence(string.Concat(Enumerable.Repeat("ATG", 8000)) + "TAA");
    T("FindOpenReadingFrames ATGx8000+TAA minLength 3 (count)", () => GenomicAnalyzer.FindOpenReadingFrames(nested, 3).Count());
    T("... + total Sequence chars", () => { long n = 0; foreach (var o in GenomicAnalyzer.FindOpenReadingFrames(nested, 3)) n += o.Sequence.Length; return n; });
}
if (what is "all" or "msa")
{
    string baseSeq = Rand(600);
    var fam = Enumerable.Range(0, 12).Select(_ => new DnaSequence(Mut(baseSeq, 0.1))).ToList();
    T("MultipleAlign 12x600", () => SequenceAligner.MultipleAlign(fam).TotalScore);
    T("MultipleAlignProgressive 12x600", () => SequenceAligner.MultipleAlignProgressive(fam).TotalScore);
    T("MultipleAlignIterative 12x600", () => SequenceAligner.MultipleAlignIterative(fam).TotalScore);
    if (args.Length > 1 && args[1] == "slow")
    {
        T("MultipleAlignConsistency 12x600", () => SequenceAligner.MultipleAlignConsistency(fam).TotalScore);
        T("MultipleAlignProbCons 12x600", () => SequenceAligner.MultipleAlignProbCons(fam.Select(x => x.Sequence)).TotalScore);
    }
    string b2 = Rand(1500);
    var big = Enumerable.Range(0, 30).Select(_ => new DnaSequence(Mut(b2, 0.1))).ToList();
    T("MultipleAlignProgressive 30x1500", () => SequenceAligner.MultipleAlignProgressive(big).TotalScore);
}

#if !BASE
if (what == "msa2")
{
    string b2 = Rand(1500);
    var big = Enumerable.Range(0, 30).Select(_ => new DnaSequence(Mut(b2, 0.1))).ToList();
    T("Progressive 30x1500 exact distances (FastDistanceCells=0)", () => SequenceAligner.MultipleAlignProgressive(big, null, new MsaOptions { FastDistanceCells = 0 }).TotalScore);
    T("Progressive 30x1500 default (k-tuple)", () => SequenceAligner.MultipleAlignProgressive(big, null, new MsaOptions()).TotalScore);
    T("Progressive 30x1500 k-tuple, ClustalWQuick for comparison", () => SequenceAligner.MultipleAlignProgressive(big, null, new MsaOptions { Distance = MsaDistanceMode.ClustalWQuick }).TotalScore);
    string b3 = Rand(600);
    var fam = Enumerable.Range(0, 12).Select(_ => new DnaSequence(Mut(b3, 0.1))).ToList();
    T("Iterative 12x600 exact", () => SequenceAligner.MultipleAlignIterative(fam, null, new MsaOptions { FastDistanceCells = 0 }).TotalScore);
    T("Iterative 12x600 k-tuple (FastDistanceCells=1)", () => SequenceAligner.MultipleAlignIterative(fam, null, new MsaOptions { FastDistanceCells = 1 }).TotalScore);
}

#endif

#if !BASE
if (what == "quality")
{
    // sum-of-pairs score of the progressive alignment with the exact distances against the k-tuple estimate (FastDistanceCells = 1), 8 families per row
    Console.WriteLine("kind        divergence  exact(avg SP)  ktuple(avg SP)  ratio   worst ratio");
    foreach (bool dna in new[] { true, false })
    foreach (double div in new[] { 0.05, 0.10, 0.20, 0.30, 0.45 })
    {
        double se = 0, sk = 0, worst = 1;
        for (int f = 0; f < 8; f++)
        {
            string alpha = dna ? "ACGT" : "ACDEFGHIKLMNPQRSTVWY";
            string bs = Rand(600, dna ? "ACGT" : alpha);
            var seqs = Enumerable.Range(0, 12).Select(_ => { var c = bs.ToCharArray(); for (int i = 0; i < c.Length; i++) if (rng.NextDouble() < div) c[i] = alpha[rng.Next(alpha.Length)]; return new string(c); }).ToList();
            double e, k;
            if (dna)
            {
                var d = seqs.Select(x => new DnaSequence(x)).ToList();
                e = SequenceAligner.MultipleAlignProgressive(d, null, new MsaOptions { FastDistanceCells = 0 }).TotalScore;
                k = SequenceAligner.MultipleAlignProgressive(d, null, new MsaOptions { FastDistanceCells = 1 }).TotalScore;
            }
            else
            {
                var sc = MatrixScoring.Emboss(SubstitutionMatrix.Blosum62);
                e = SequenceAligner.MultipleAlignProgressive(seqs, sc, new MsaOptions { FastDistanceCells = 0 }).TotalScore;
                k = SequenceAligner.MultipleAlignProgressive(seqs, sc, new MsaOptions { FastDistanceCells = 1 }).TotalScore;
            }
            se += e; sk += k; worst = Math.Min(worst, e == 0 ? 1 : (k / e));
        }
        Console.WriteLine($"{(dna ? "DNA" : "protein"),-10} {div,8:F2} {se / 8,14:F0} {sk / 8,15:F0} {sk / se,8:F4} {worst,9:F4}");
    }
}
#endif
if (what == "golden")
{
    var gr = new Random(780062);
    string GRand(int n, string a) { var c = new char[n]; for (int i = 0; i < n; i++) c[i] = a[gr.Next(a.Length)]; return new string(c); }
    string GMut(string s0, string a, double r, double ind) { var o = new System.Text.StringBuilder(); foreach (char c in s0) { double x = gr.NextDouble(); if (x < ind / 2) continue; if (x < ind) { o.Append(c).Append(a[gr.Next(a.Length)]); continue; } o.Append(gr.NextDouble() < r ? a[gr.Next(a.Length)] : c); } return o.Length == 0 ? s0[..1] : o.ToString(); }
    var lines = new List<string>();
    var protScoring = MatrixScoring.Emboss(SubstitutionMatrix.Blosum62);
    for (int k = 0; k < 117; k++)
    {
        bool dna = k < 60;
        string alpha = dna ? "ACGT" : "ACDEFGHIKLMNPQRSTVWY";
        int n = gr.Next(3, 9);
        int len = gr.Next(30, dna ? 160 : 110);
        string bs = GRand(len, alpha);
        double r = new[] { 0.05, 0.15, 0.3 }[gr.Next(3)], ind = new[] { 0.02, 0.06, 0.12 }[gr.Next(3)];
        var seqs = Enumerable.Range(0, n).Select(_ => GMut(bs, alpha, r, ind)).ToList();
        string kind, opts;
        MultipleAlignmentResult res;
        if (dna)
        {
            var dnaSeqs = seqs.Select(x => new DnaSequence(x)).ToList();
            switch (k % 6)
            {
                case 0: kind = "progressive"; opts = "linear"; res = SequenceAligner.MultipleAlignProgressive(dnaSeqs); break;
                case 1: kind = "iterative"; opts = "linear"; res = SequenceAligner.MultipleAlignIterative(dnaSeqs); break;
                case 2: kind = "progressive"; opts = "affine"; res = SequenceAligner.MultipleAlignProgressive(dnaSeqs, null, new MsaOptions { GapModel = MsaGapModel.Affine }); break;
                case 3: kind = "progressive"; opts = "nj"; res = SequenceAligner.MultipleAlignProgressive(dnaSeqs, null, new MsaOptions { GuideTree = MsaGuideTree.NeighborJoining }); break;
                case 4: kind = "consistency"; opts = "linear"; res = SequenceAligner.MultipleAlignConsistency(dnaSeqs.Take(5)); break;
                default: kind = "iterative"; opts = "affine"; res = SequenceAligner.MultipleAlignIterative(dnaSeqs, null, new MsaOptions { GapModel = MsaGapModel.Affine }); break;
            }
        }
        else
        {
            switch (k % 4)
            {
                case 0: kind = "progressive-protein"; opts = "linear"; res = SequenceAligner.MultipleAlignProgressive(seqs, protScoring); break;
                case 1: kind = "iterative-protein"; opts = "linear"; res = SequenceAligner.MultipleAlignIterative(seqs, protScoring); break;
                case 2: kind = "progressive-protein"; opts = "affine"; res = SequenceAligner.MultipleAlignProgressive(seqs, protScoring, new MsaOptions { GapModel = MsaGapModel.Affine }); break;
                default: kind = "progressive-protein"; opts = "nj"; res = SequenceAligner.MultipleAlignProgressive(seqs, protScoring, new MsaOptions { GuideTree = MsaGuideTree.NeighborJoining }); break;
            }
            if (kind == "progressive-protein" || kind == "iterative-protein") seqs = seqs.ToList();
        }
        if (kind == "consistency") seqs = seqs.Take(5).ToList();
        lines.Add(string.Join('\t', kind, opts, string.Join('|', seqs), string.Join('|', res.AlignedSequences), res.TotalScore.ToString("R", System.Globalization.CultureInfo.InvariantCulture)));
    }
    if (args.Length > 2 && args[1] == "write")
    {
        using var gz = new StreamWriter(new System.IO.Compression.GZipStream(File.Create(args[2]), System.IO.Compression.CompressionLevel.SmallestSize), new System.Text.UTF8Encoding(false)) { NewLine = "\n" };
        foreach (var l in lines) gz.WriteLine(l);
        Console.WriteLine("written " + lines.Count);
    }
    else if (args.Length > 2 && args[1] == "check")
    {
        var expected = new List<string>();
        using (var rd = new StreamReader(new System.IO.Compression.GZipStream(File.OpenRead(args[2]), System.IO.Compression.CompressionMode.Decompress))) { string? l; while ((l = rd.ReadLine()) != null) expected.Add(l); }
        int bad = 0;
        for (int i = 0; i < lines.Count; i++) if (i >= expected.Count || lines[i] != expected[i]) { bad++; Console.WriteLine("DIFF case " + i + " " + lines[i].Split('\t')[0] + " " + lines[i].Split('\t')[1]); }
        Console.WriteLine(lines.Count + " cases, " + bad + " differ");
    }
}
