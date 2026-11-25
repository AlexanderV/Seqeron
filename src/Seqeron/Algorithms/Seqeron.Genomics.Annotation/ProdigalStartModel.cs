namespace Seqeron.Genomics.Annotation;

/// <summary>How <see cref="GenomeAnnotator.FindRibosomeBindingSites"/> scores a Shine-Dalgarno hit.</summary>
public enum RbsScoreMode
{
    /// <summary>Legacy default: motif length / 6 (a motif-length heuristic, no spacer or genome model).</summary>
    MotifLength = 0,

    /// <summary>
    /// Prodigal log-odds weight of the hit's (motif, spacer) bin (Hyatt et al. 2010, BMC Bioinformatics 11:119):
    /// <c>ln(freq of the bin among gene starts / freq among all candidate starts)</c>, clamped to ±4, trained on the
    /// input genome itself (<see cref="RbsModel"/>).
    /// </summary>
    LogOdds = 1,
}

/// <summary>How <see cref="GenomeAnnotator.PredictGenes"/> chooses the start codon of a gene.</summary>
public enum StartSelection
{
    /// <summary>Legacy default: the most upstream in-frame ATG/GTG/TTG (longest ORF) per stop.</summary>
    LongestOrf = 0,

    /// <summary>
    /// Per stop, the start with the highest Prodigal start score (coding score + RBS log-odds + start-type weight +
    /// upstream composition) under a model trained on the input genome (<see cref="RbsModel"/>).
    /// </summary>
    ProdigalScore = 1,

    /// <summary>
    /// Prodigal's complete gene model: dynamic programming over all start/stop nodes with the trained scores followed
    /// by Prodigal's final start adjustment. The gene set is Prodigal's (minimum 90 nt incl. stop), not the ORF set.
    /// </summary>
    ProdigalDynamicProgramming = 2,
}

/// <summary>
/// Genome-trained Shine-Dalgarno / start-type model of Prodigal (Hyatt et al. 2010, BMC Bioinformatics 11:119), obtained
/// by Prodigal's self-supervised training on the input sequence (<see cref="GenomeAnnotator.TrainRbsModel(string)"/>).
/// </summary>
/// <remarks>
/// <b>28 RBS bins.</b> Bin 0 = no motif. Exact-match bins (a contiguous sub-motif of AGGAGG of length 3-6, scored
/// 2 for A and 3 for G minus 2, with spacer class 0 = 5-10 nt, 1 = 11-12 nt short / &lt;5 nt long, ... see
/// <see cref="GenomeAnnotator.ClassifyRbsStart"/>) are bins 1,2,3,6,10-13,15,16,20-27; the 5/6-mers with one mismatch
/// are bins 2-5,7-9,14,17-19. <see cref="BinWeights"/>[b] is the clamped natural-log odds of bin b.
/// </remarks>
public sealed class RbsModel
{
    internal RbsModel(double gc, double[] bias, double[] typeWt, double[] rbsWt, double[,] upsComp, double[]? geneDc,
        bool usesSd, double[] rbsReal, double[] rbsBackground, double[] typeReal)
    {
        GcContent = gc;
        Bias = bias;
        TypeWt = typeWt;
        RbsWt = rbsWt;
        UpsComp = upsComp;
        GeneDc = geneDc;
        UsesShineDalgarno = usesSd;
        BinFrequencyInGenes = rbsReal;
        BinFrequencyBackground = rbsBackground;
        TypeFrequencyInGenes = typeReal;
    }

    /// <summary>Prodigal start weight (<c>st_wt</c> = 4.35).</summary>
    public const double StartWeight = 4.35;

    internal double GcContent { get; }
    internal double[] Bias { get; }
    internal double[] TypeWt { get; }
    internal double[] RbsWt { get; }
    internal double[,] UpsComp { get; }
    internal double[]? GeneDc { get; }

    /// <summary>Natural-log odds weight of each of the 28 RBS bins (clamped to ±4).</summary>
    public IReadOnlyList<double> BinWeights => RbsWt;

    /// <summary>Log-odds weight of the start codon: ATG, GTG, TTG (clamped to ±4).</summary>
    public IReadOnlyList<double> StartTypeWeights => TypeWt;

    /// <summary>Prodigal's decision whether the genome uses Shine-Dalgarno motifs (<c>determine_sd_usage</c>).</summary>
    public bool UsesShineDalgarno { get; }

    /// <summary>Frequency of each RBS bin among the training genes' starts (last iteration; empty for gene-set training with no genes).</summary>
    public IReadOnlyList<double> BinFrequencyInGenes { get; }

    /// <summary>Background frequency of each RBS bin among all candidate starts (last iteration).</summary>
    public IReadOnlyList<double> BinFrequencyBackground { get; }

    /// <summary>Frequency of ATG, GTG, TTG among the training genes' starts (last iteration).</summary>
    public IReadOnlyList<double> TypeFrequencyInGenes { get; }
}

/// <summary>Best RBS bins found upstream of one start codon (Prodigal <c>rbs_score</c>).</summary>
/// <param name="ExactBin">Best exact-motif bin (0 = none).</param>
/// <param name="MismatchBin">Best one-mismatch bin (0 = none).</param>
/// <param name="Bin">Bin chosen by Prodigal's weight rule (exact vs mismatch).</param>
/// <param name="Weight">Log-odds weight of <paramref name="Bin"/> (0 when no model was given).</param>
public readonly record struct RbsStartClass(int ExactBin, int MismatchBin, int Bin, double Weight);

public static partial class GenomeAnnotator
{
    /// <summary>
    /// Trains Prodigal's start model on the input genome by Prodigal's own self-supervised procedure (single-genome
    /// training mode, Hyatt et al. 2010): all start/stop nodes, GC-frame-bias dynamic programming for an initial gene
    /// set, dicodon coding model, then ten iterations of <c>train_starts_sd</c> re-choosing each ORF's start by
    /// coding + RBS + start-type score and re-estimating the 28 RBS-bin and 3 start-type log-odds weights.
    /// </summary>
    /// <param name="genomeSequence">One contiguous genome (≥ 20 000 nt, as Prodigal requires); case-insensitive, non-ACGT read as C (Prodigal).</param>
    /// <exception cref="ArgumentException">Sequence shorter than 20 000 nt.</exception>
    /// <exception cref="NotSupportedException">The genome does not use Shine-Dalgarno motifs (Prodigal would switch to its non-SD motif trainer, which is not ported).</exception>
    public static RbsModel TrainRbsModel(string genomeSequence)
    {
        ArgumentNullException.ThrowIfNull(genomeSequence);
        if (genomeSequence.Length < 20000)
            throw new ArgumentException("Prodigal training needs at least 20000 nt (use a gene-set overload for short input).", nameof(genomeSequence));
        return ProdigalEngine.Train(genomeSequence);
    }

    /// <summary>
    /// Trains the RBS-bin and start-type log-odds weights from a caller-supplied gene set instead of Prodigal's
    /// self-supervised gene set: the same <c>train_starts_sd</c> estimator, with the given gene starts as the "real"
    /// starts and every candidate ATG/GTG/TTG start of the genome as the background (the genes' coordinates are fixed,
    /// so there is no coding-based start re-choice). The coding model is not trained (<see cref="StartSelection"/> modes
    /// that need it call <see cref="TrainRbsModel(string)"/>).
    /// </summary>
    /// <param name="genomeSequence">Genome sequence.</param>
    /// <param name="genes">Genes as returned by <see cref="PredictGenes"/>: 0-based <c>Start</c>, exclusive <c>End</c>, forward coordinates, <c>'+'</c>/<c>'-'</c> strand.</param>
    public static RbsModel TrainRbsModel(string genomeSequence, IEnumerable<GeneAnnotation> genes)
    {
        ArgumentNullException.ThrowIfNull(genomeSequence);
        ArgumentNullException.ThrowIfNull(genes);
        return ProdigalEngine.TrainFromGenes(genomeSequence, genes);
    }

    /// <summary>
    /// Classifies the Shine-Dalgarno site upstream of one start codon exactly as Prodigal's <c>rbs_score</c>: the best
    /// exact and one-mismatch AGGAGG-like motif in the window 6-20 nt upstream of the start, mapped to its bin.
    /// </summary>
    /// <param name="sequence">Genome sequence.</param>
    /// <param name="startIndex">0-based index of the first base of the start codon in forward coordinates; for a reverse-strand gene, the highest coordinate of the start codon (its first base on the reverse strand).</param>
    /// <param name="strand"><c>'+'</c> or <c>'-'</c>.</param>
    /// <param name="model">Optional model supplying the weights (also used to choose between the exact and mismatch bin).</param>
    public static RbsStartClass ClassifyRbsStart(string sequence, int startIndex, char strand = '+', RbsModel? model = null)
    {
        ArgumentNullException.ThrowIfNull(sequence);
        if (strand != '+' && strand != '-') throw new ArgumentException("strand must be '+' or '-'.", nameof(strand));
        if (startIndex < 0 || startIndex >= sequence.Length) throw new ArgumentOutOfRangeException(nameof(startIndex));
        return ProdigalEngine.ClassifyStart(sequence, startIndex, strand == '+' ? 1 : -1, model);
    }
}

/// <summary>
/// Clean-room port of the single-genome Prodigal pipeline (Hyatt et al. 2010, BMC Bioinformatics 11:119) from the
/// published algorithm; the formulas and constants (bin definitions, scoring, weights) are the specification.
/// Sequence encoding follows Prodigal: A=0, G=1, C=2, T=3 (complement = 3 - x), non-ACGT = C.
/// </summary>
internal static class ProdigalEngine
{
    private const int STOP = 3;
    private const int MinGene = 90, MinEdgeGene = 60, MaxSamOvlp = 60, OperDist = 60;
    private const double EdgeBonus = 0.74, EdgeUps = -1.0;
    private const int MaxNodeDist = 500, MaxOppOvlp = 200, GcWindow = 120;

    internal sealed class Node
    {
        public int Type, Edge, Ndx, Strand, StopVal, GcBias;
        public readonly int[] StarPtr = new int[3];
        public readonly double[] GcScore = new double[3];
        public double Cscore, Uscore, Tscore, Rscore, Sscore, Score;
        public readonly int[] Rbs = new int[2];
        public int Traceb = -1, Tracef = -1, OvMark = -1, Elim;
    }

    internal sealed class Training
    {
        public double Gc, StWt = 4.35;
        public double[] Bias = new double[3];
        public double[] TypeWt = new double[3];
        public double[] RbsWt = new double[28];
        public double[,] UpsComp = new double[32, 4];
        public double[] GeneDc = new double[4096];
        public bool UsesSd;
    }

    internal sealed class Seq
    {
        public required byte[] F { get; init; }
        public required byte[] R { get; init; }
        public int Len => F.Length;
        public static bool Stop(byte[] s, int n) => B(s, n) == 3 && ((B(s, n + 1) == 0 && B(s, n + 2) == 1) || (B(s, n + 1) == 1 && B(s, n + 2) == 0) || (B(s, n + 1) == 0 && B(s, n + 2) == 0));
        public static bool Start(byte[] s, int n) => (B(s, n) == 0 || B(s, n) == 1 || B(s, n) == 3) && B(s, n + 1) == 3 && B(s, n + 2) == 1;
        public static byte B(byte[] s, int n) => (uint)n < (uint)s.Length ? s[n] : (byte)0;
        public static int StartType(byte[] s, int n)
        {
            byte b = B(s, n);
            if (b == 0) return 0;
            return b == 1 ? 1 : 2;
        }
    }

    private static bool IsGc(byte[] s, int n) { byte b = Seq.B(s, n); return b == 1 || b == 2; }

    internal static Seq Encode(string dna)
    {
        var f = new byte[dna.Length];
        for (int i = 0; i < dna.Length; i++)
            f[i] = char.ToUpperInvariant(dna[i]) switch { 'A' => (byte)0, 'G' => (byte)1, 'T' or 'U' => (byte)3, _ => (byte)2 };
        var r = new byte[f.Length];
        for (int i = 0; i < f.Length; i++) r[f.Length - 1 - i] = (byte)(3 - f[i]);
        return new Seq { F = f, R = r };
    }

    private static int MerNdx(int len, byte[] s, int pos)
    {
        int ndx = 0;
        for (int k = 0; k < len; k++) ndx |= Seq.B(s, pos + k) << (2 * k);
        return ndx;
    }

    // ------------------------------------------------------------------ nodes
    private static List<Node> AddNodes(Seq q, bool closed)
    {
        var nodes = new List<Node>();
        int slen = q.Len, slmod = slen % 3;
        var last = new int[3]; var sawStart = new bool[3]; var minDist = new int[3];
        for (int strandPass = 0; strandPass < 2; strandPass++)
        {
            byte[] s = strandPass == 0 ? q.F : q.R;
            int strand = strandPass == 0 ? 1 : -1;
            for (int i = 0; i < 3; i++)
            {
                last[(i + slmod) % 3] = slen + i; sawStart[i % 3] = false; minDist[i % 3] = MinEdgeGene;
                if (!closed) while (last[(i + slmod) % 3] + 2 > slen - 1) last[(i + slmod) % 3] -= 3;
            }
            int Nd(int p) => strand == 1 ? p : slen - p - 1;
            for (int i = slen - 3; i >= 0; i--)
            {
                int m = i % 3;
                if (Seq.Stop(s, i))
                {
                    if (sawStart[m])
                    {
                        var n = new Node { Type = STOP, Strand = strand, Ndx = Nd(last[m]), StopVal = strand == 1 ? i : slen - i - 1 };
                        if (!Seq.Stop(s, last[m])) n.Edge = 1;
                        nodes.Add(n);
                    }
                    minDist[m] = MinGene; last[m] = i; sawStart[m] = false;
                    continue;
                }
                if (last[m] >= slen) continue;
                int sv = strand == 1 ? last[m] : slen - last[m] - 1;
                if (Seq.Start(s, i) && last[m] - i + 3 >= minDist[m])
                {
                    nodes.Add(new Node { Ndx = Nd(i), Type = Seq.StartType(s, i), Strand = strand, StopVal = sv });
                    sawStart[m] = true;
                }
                else if (i <= 2 && !closed && last[m] - i > MinEdgeGene)
                {
                    nodes.Add(new Node { Ndx = Nd(i), Type = 0, Strand = strand, StopVal = sv, Edge = 1 });
                    sawStart[m] = true;
                }
            }
            for (int i = 0; i < 3; i++)
            {
                if (!sawStart[i % 3]) continue;
                var n = new Node { Type = STOP, Strand = strand, Ndx = Nd(last[i % 3]), StopVal = strand == 1 ? i - 6 : slen - i + 5 };
                if (!Seq.Stop(s, last[i % 3])) n.Edge = 1;
                nodes.Add(n);
            }
        }
        // Prodigal's qsort(compare_nodes): ndx ascending, then strand 1 before -1 (stable here).
        return nodes.Select((n, k) => (n, k)).OrderBy(t => t.n.Ndx).ThenByDescending(t => t.n.Strand).ThenBy(t => t.k).Select(t => t.n).ToList();
    }

    private static int MaxFr(int a, int b, int c)
    {
        if (a > b) return a > c ? 0 : 2;
        return b > c ? 1 : 2;
    }

    private static int[] GcFramePlot(byte[] s)
    {
        int slen = s.Length;
        var fwd = new int[slen]; var bwd = new int[slen]; var tot = new int[slen]; var gp = new int[slen];
        Array.Fill(gp, -1);
        for (int j = 0; j < slen; j++)
        {
            fwd[j] = (j < 3 ? 0 : fwd[j - 3]) + (IsGc(s, j) ? 1 : 0);
            int k = slen - j - 1;
            bwd[k] = (j < 3 ? 0 : bwd[slen - j + 2]) + (IsGc(s, k) ? 1 : 0);
        }
        for (int i = 0; i < slen; i++)
        {
            tot[i] = fwd[i] + bwd[i] - (IsGc(s, i) ? 1 : 0);
            if (i - GcWindow / 2 >= 0) tot[i] -= fwd[i - GcWindow / 2];
            if (i + GcWindow / 2 < slen) tot[i] -= bwd[i + GcWindow / 2];
        }
        for (int i = 0; i < slen - 2; i += 3)
        {
            int win = MaxFr(tot[i], tot[i + 1], tot[i + 2]);
            gp[i] = gp[i + 1] = gp[i + 2] = win;
        }
        return gp;
    }

    private static void RecordGcBias(int[] gc, List<Node> nod, Training t)
    {
        int nn = nod.Count; if (nn == 0) return;
        var ctr = new int[3, 3]; var last = new int[3];
        for (int i = nn - 1; i >= 0; i--)
        {
            var n = nod[i]; int fr = n.Ndx % 3, frmod = 3 - fr;
            if (n.Strand == 1 && n.Type == STOP)
            {
                for (int j = 0; j < 3; j++) ctr[fr, j] = 0;
                last[fr] = n.Ndx;
                ctr[fr, (gc[n.Ndx] + frmod) % 3] = 1;
            }
            else if (n.Strand == 1)
            {
                for (int j = last[fr] - 3; j >= n.Ndx; j -= 3) ctr[fr, (gc[j] + frmod) % 3]++;
                n.GcBias = MaxFr(ctr[fr, 0], ctr[fr, 1], ctr[fr, 2]);
                for (int j = 0; j < 3; j++) n.GcScore[j] = 3.0 * ctr[fr, j] / (1.0 * (n.StopVal - n.Ndx + 3));
                last[fr] = n.Ndx;
            }
        }
        for (int i = 0; i < nn; i++)
        {
            var n = nod[i]; int fr = n.Ndx % 3, frmod = fr;
            if (n.Strand == -1 && n.Type == STOP)
            {
                for (int j = 0; j < 3; j++) ctr[fr, j] = 0;
                last[fr] = n.Ndx;
                ctr[fr, (3 - gc[n.Ndx] + frmod) % 3] = 1;
            }
            else if (n.Strand == -1)
            {
                for (int j = last[fr] + 3; j <= n.Ndx; j += 3) ctr[fr, (3 - gc[j] + frmod) % 3]++;
                n.GcBias = MaxFr(ctr[fr, 0], ctr[fr, 1], ctr[fr, 2]);
                for (int j = 0; j < 3; j++) n.GcScore[j] = 3.0 * ctr[fr, j] / (1.0 * (n.Ndx - n.StopVal + 3));
                last[fr] = n.Ndx;
            }
        }
        Array.Clear(t.Bias);
        foreach (var n in nod)
        {
            if (n.Type == STOP) continue;
            int len = Math.Abs(n.StopVal - n.Ndx) + 1;
            t.Bias[n.GcBias] += n.GcScore[n.GcBias] * len / 1000.0;
        }
        double tot = t.Bias[0] + t.Bias[1] + t.Bias[2];
        for (int i = 0; i < 3; i++) t.Bias[i] *= 3.0 / tot;
    }

    private static double IntergenicMod(Node n1, Node n2, Training t)
    {
        double rval = 0.0; bool ovlp = false;
        if ((n1.Strand == 1 && n2.Strand == 1 && (n1.Ndx + 2 == n2.Ndx || n1.Ndx - 1 == n2.Ndx)) ||
            (n1.Strand == -1 && n2.Strand == -1 && (n1.Ndx + 2 == n2.Ndx || n1.Ndx - 1 == n2.Ndx)))
        {
            if (n1.Strand == 1 && n2.Rscore < 0) rval -= n2.Rscore;
            if (n1.Strand == -1 && n1.Rscore < 0) rval -= n1.Rscore;
            if (n1.Strand == 1 && n2.Uscore < 0) rval -= n2.Uscore;
            if (n1.Strand == -1 && n1.Uscore < 0) rval -= n1.Uscore;
        }
        int dist = Math.Abs(n1.Ndx - n2.Ndx);
        if (n1.Strand == 1 && n2.Strand == 1 && n1.Ndx + 2 >= n2.Ndx) ovlp = true;
        else if (n1.Strand == -1 && n2.Strand == -1 && n1.Ndx >= n2.Ndx + 2) ovlp = true;
        if (dist > 3 * OperDist || n1.Strand != n2.Strand) rval -= 0.15 * t.StWt;
        else if ((dist <= OperDist && !ovlp) || dist < 0.25 * OperDist) rval += (2.0 - (double)dist / OperDist) * 0.15 * t.StWt;
        return rval;
    }

    private static void RecordOverlappingStarts(List<Node> nod, Training t, bool flag)
    {
        int nn = nod.Count;
        for (int i = 0; i < nn; i++)
        {
            var ni = nod[i];
            ni.StarPtr[0] = ni.StarPtr[1] = ni.StarPtr[2] = -1;
            if (ni.Type != STOP || ni.Edge == 1) continue;
            double maxSc = -100;
            if (ni.Strand == 1)
            {
                for (int j = i + 3; j >= 0; j--)
                {
                    if (j >= nn || nod[j].Ndx > ni.Ndx + 2) continue;
                    if (nod[j].Ndx + MaxSamOvlp < ni.Ndx) break;
                    var nj = nod[j];
                    if (nj.Strand == 1 && nj.Type != STOP)
                    {
                        if (nj.StopVal <= ni.Ndx) continue;
                        if (!flag && ni.StarPtr[nj.Ndx % 3] == -1) ni.StarPtr[nj.Ndx % 3] = j;
                        else if (flag && nj.Cscore + nj.Sscore + IntergenicMod(ni, nj, t) > maxSc)
                        {
                            ni.StarPtr[nj.Ndx % 3] = j;
                            maxSc = nj.Cscore + nj.Sscore + IntergenicMod(ni, nj, t);
                        }
                    }
                }
            }
            else
            {
                for (int j = i - 3; j < nn; j++)
                {
                    if (j < 0 || nod[j].Ndx < ni.Ndx - 2) continue;
                    if (nod[j].Ndx - MaxSamOvlp > ni.Ndx) break;
                    var nj = nod[j];
                    if (nj.Strand == -1 && nj.Type != STOP)
                    {
                        if (nj.StopVal >= ni.Ndx) continue;
                        if (!flag && ni.StarPtr[nj.Ndx % 3] == -1) ni.StarPtr[nj.Ndx % 3] = j;
                        else if (flag && nj.Cscore + nj.Sscore + IntergenicMod(nj, ni, t) > maxSc)
                        {
                            ni.StarPtr[nj.Ndx % 3] = j;
                            maxSc = nj.Cscore + nj.Sscore + IntergenicMod(nj, ni, t);
                        }
                    }
                }
            }
        }
    }

    // ------------------------------------------------------------------ dynamic programming
    private static int DProg(List<Node> nod, Training t, bool flag)
    {
        int nn = nod.Count; if (nn == 0) return -1;
        foreach (var n in nod) { n.Score = 0; n.Traceb = -1; n.Tracef = -1; }
        for (int i = 0; i < nn; i++)
        {
            int min = i < MaxNodeDist ? 0 : i - MaxNodeDist;
            var ni = nod[i];
            if (ni.Strand == -1 && ni.Type != STOP && nod[min].Ndx >= ni.StopVal)
                while (min >= 0 && ni.Ndx != ni.StopVal) min--;
            if (ni.Strand == 1 && ni.Type == STOP && nod[Math.Max(min, 0)].Ndx >= ni.StopVal)
                while (min >= 0 && ni.Ndx != ni.StopVal) min--;
            min = min < MaxNodeDist ? 0 : min - MaxNodeDist;
            for (int j = min; j < i; j++) ScoreConnection(nod, j, i, t, flag);
        }
        double maxSc = -1.0; int maxNdx = -1;
        for (int i = nn - 1; i >= 0; i--)
        {
            if (nod[i].Strand == 1 && nod[i].Type != STOP) continue;
            if (nod[i].Strand == -1 && nod[i].Type == STOP) continue;
            if (nod[i].Score > maxSc) { maxSc = nod[i].Score; maxNdx = i; }
        }
        if (maxNdx < 0) return -1;
        int path = maxNdx;
        while (nod[path].Traceb != -1)
        {
            int nxt = nod[path].Traceb;
            if (nod[path].Strand == -1 && nod[path].Type == STOP && nod[nxt].Strand == 1 && nod[nxt].Type == STOP &&
                nod[path].OvMark != -1 && nod[path].Ndx > nod[nxt].Ndx)
            {
                int tmp = nod[path].StarPtr[nod[path].OvMark];
                int i = tmp; while (nod[i].Ndx != nod[tmp].StopVal) i--;
                nod[path].Traceb = tmp; nod[tmp].Traceb = i; nod[i].OvMark = -1; nod[i].Traceb = nxt;
            }
            path = nod[path].Traceb;
        }
        path = maxNdx;
        while (nod[path].Traceb != -1)
        {
            int nxt = nod[path].Traceb;
            var p = nod[path]; var x = nod[nxt];
            if (p.Strand == -1 && p.Type != STOP && x.Strand == 1 && x.Type == STOP)
            {
                int i = path; while (nod[i].Ndx != p.StopVal) i--;
                p.Traceb = i; nod[i].Traceb = nxt;
            }
            if (p.Strand == 1 && p.Type == STOP && x.Strand == 1 && x.Type == STOP)
            {
                p.Traceb = x.StarPtr[p.Ndx % 3]; nod[p.Traceb].Traceb = nxt;
            }
            if (p.Strand == -1 && p.Type == STOP && x.Strand == -1 && x.Type == STOP)
            {
                p.Traceb = p.StarPtr[x.Ndx % 3]; nod[p.Traceb].Traceb = nxt;
            }
            path = nod[path].Traceb;
        }
        path = maxNdx;
        while (nod[path].Traceb != -1) { nod[nod[path].Traceb].Tracef = path; path = nod[path].Traceb; }
        return nod[maxNdx].Traceb == -1 ? -1 : maxNdx;
    }

    private static double GcMod(Node n, Training t) => t.Bias[0] * n.GcScore[0] + t.Bias[1] * n.GcScore[1] + t.Bias[2] * n.GcScore[2];

    private static void ScoreConnection(List<Node> nod, int p1, int p2, Training t, bool flag)
    {
        Node n1 = nod[p1], n2 = nod[p2], n3;
        int left = n1.Ndx, right = n2.Ndx, ovlp = 0, maxfr = -1;
        double score = 0.0, scrMod = 0.0;

        if (n1.Type != STOP && n2.Type != STOP && n1.Strand == n2.Strand) return;
        else if (n1.Strand == 1 && n1.Type != STOP && n2.Strand == -1) return;
        else if (n1.Strand == -1 && n1.Type == STOP && n2.Strand == 1) return;
        else if (n1.Strand == -1 && n1.Type != STOP && n2.Strand == 1 && n2.Type == STOP) return;

        if (n1.Traceb == -1 && n1.Strand == 1 && n1.Type == STOP) return;
        if (n1.Traceb == -1 && n1.Strand == -1 && n1.Type != STOP) return;

        if (n1.Strand == n2.Strand && n1.Strand == 1 && n1.Type != STOP && n2.Type == STOP)
        {
            if (n2.StopVal >= n1.Ndx) return;
            if (n1.Ndx % 3 != n2.Ndx % 3) return;
            right += 2;
            if (!flag) scrMod = GcMod(n1, t); else score = n1.Cscore + n1.Sscore;
        }
        else if (n1.Strand == n2.Strand && n1.Strand == -1 && n1.Type == STOP && n2.Type != STOP)
        {
            if (n1.StopVal <= n2.Ndx) return;
            if (n1.Ndx % 3 != n2.Ndx % 3) return;
            left -= 2;
            if (!flag) scrMod = GcMod(n2, t); else score = n2.Cscore + n2.Sscore;
        }
        else if (n1.Strand == 1 && n1.Type == STOP && n2.Strand == 1 && n2.Type != STOP)
        {
            left += 2;
            if (left >= right) return;
            if (flag) score = IntergenicMod(n1, n2, t);
        }
        else if (n1.Strand == 1 && n1.Type == STOP && n2.Strand == -1 && n2.Type == STOP)
        {
            left += 2; right -= 2;
            if (left >= right) return;
            double maxval = 0.0; maxfr = -1;
            for (int i = 0; i < 3; i++)
            {
                if (n2.StarPtr[i] == -1) continue;
                n3 = nod[n2.StarPtr[i]];
                ovlp = left - n3.StopVal + 3;
                if (ovlp <= 0 || ovlp >= MaxOppOvlp) continue;
                if (ovlp >= n3.Ndx - left) continue;
                if (n1.Traceb == -1) continue;
                if (ovlp >= n3.StopVal - nod[n1.Traceb].Ndx - 2) continue;
                if ((flag && n3.Cscore + n3.Sscore + IntergenicMod(n3, n2, t) > maxval) || (!flag && GcMod(n3, t) > maxval))
                {
                    maxfr = i;
                    maxval = n3.Cscore + n3.Sscore + IntergenicMod(n3, n2, t);
                }
            }
            if (maxfr != -1)
            {
                n3 = nod[n2.StarPtr[maxfr]];
                if (!flag) scrMod = GcMod(n3, t); else score = n3.Cscore + n3.Sscore + IntergenicMod(n3, n2, t);
            }
            else if (flag) score = IntergenicMod(n1, n2, t);
        }
        else if (n1.Strand == -1 && n1.Type != STOP && n2.Strand == -1 && n2.Type == STOP)
        {
            right -= 2;
            if (left >= right) return;
            if (flag) score = IntergenicMod(n1, n2, t);
        }
        else if (n1.Strand == -1 && n1.Type != STOP && n2.Strand == 1 && n2.Type != STOP)
        {
            if (left >= right) return;
            if (flag) score = IntergenicMod(n1, n2, t);
        }
        else if (n1.Strand == 1 && n2.Strand == 1 && n1.Type == STOP && n2.Type == STOP)
        {
            if (n2.StopVal >= n1.Ndx) return;
            if (n1.StarPtr[n2.Ndx % 3] == -1) return;
            n3 = nod[n1.StarPtr[n2.Ndx % 3]];
            left = n3.Ndx; right += 2;
            if (!flag) scrMod = GcMod(n3, t); else score = n3.Cscore + n3.Sscore + IntergenicMod(n1, n3, t);
        }
        else if (n1.Strand == -1 && n1.Type == STOP && n2.Strand == -1 && n2.Type == STOP)
        {
            if (n1.StopVal <= n2.Ndx) return;
            if (n2.StarPtr[n1.Ndx % 3] == -1) return;
            n3 = nod[n2.StarPtr[n1.Ndx % 3]];
            left -= 2; right = n3.Ndx;
            if (!flag) scrMod = GcMod(n3, t); else score = n3.Cscore + n3.Sscore + IntergenicMod(n3, n2, t);
        }
        else if (n1.Strand == 1 && n1.Type == STOP && n2.Strand == -1 && n2.Type != STOP)
        {
            if (n2.StopVal - 2 >= n1.Ndx + 2) return;
            ovlp = (n1.Ndx + 2) - (n2.StopVal - 2) + 1;
            if (ovlp >= MaxOppOvlp) return;
            if ((n1.Ndx + 2 - n2.StopVal - 2 + 1) >= (n2.Ndx - n1.Ndx + 3 + 1)) return;
            int bnd = n1.Traceb == -1 ? 0 : nod[n1.Traceb].Ndx;
            if ((n1.Ndx + 2 - n2.StopVal - 2 + 1) >= (n2.StopVal - 3 - bnd + 1)) return;
            left = n2.StopVal - 2;
            if (!flag) scrMod = GcMod(n2, t); else score = n2.Cscore + n2.Sscore - 0.15 * t.StWt;
        }

        if (!flag) score = (right - left + 1 - ovlp * 2) * scrMod;

        if (n1.Score + score >= n2.Score)
        {
            n2.Score = n1.Score + score;
            n2.Traceb = p1;
            n2.OvMark = maxfr;
        }
    }

    // ------------------------------------------------------------------ coding model
    private static void CalcDicodonGene(Training t, Seq q, List<Node> nod, int dbeg)
    {
        int slen = q.Len;
        var counts = new int[4096]; int glob = 0;
        var bg = new double[4096];
        int bgTotal = 0;
        var bgCounts = new int[4096];
        for (int i = 0; i < slen - 6 + 1; i++)
        {
            bgCounts[MerNdx(6, q.F, i)]++; bgCounts[MerNdx(6, q.R, i)]++; bgTotal += 2;
        }
        for (int i = 0; i < 4096; i++) bg[i] = bgCounts[i] * 1.0 / (bgTotal * 1.0);
        int left = -1, right = -1, inGene = 0, path = dbeg;
        while (path != -1)
        {
            var n = nod[path];
            if (n.Strand == -1 && n.Type != STOP) { inGene = -1; left = slen - n.Ndx - 1; }
            if (n.Strand == 1 && n.Type == STOP) { inGene = 1; right = n.Ndx + 2; }
            if (inGene == -1 && n.Strand == -1 && n.Type == STOP)
            {
                right = slen - n.Ndx + 1;
                for (int i = left; i < right - 5; i += 3) { counts[MerNdx(6, q.R, i)]++; glob++; }
                inGene = 0;
            }
            if (inGene == 1 && n.Strand == 1 && n.Type != STOP)
            {
                left = n.Ndx;
                for (int i = left; i < right - 5; i += 3) { counts[MerNdx(6, q.F, i)]++; glob++; }
                inGene = 0;
            }
            path = n.Traceb;
        }
        for (int i = 0; i < 4096; i++)
        {
            double prob = counts[i] * 1.0 / (glob * 1.0);
            double v;
            if (prob == 0 && bg[i] != 0) v = -5.0;
            else if (bg[i] == 0) v = 0.0;
            else v = Math.Log(prob / bg[i]);
            t.GeneDc[i] = Math.Clamp(v, -5.0, 5.0);
        }
    }

    private static void RawCodingScore(Seq q, List<Node> nod, Training t)
    {
        int nn = nod.Count, slen = q.Len;
        double noStop = (1 - t.Gc) * (1 - t.Gc) * t.Gc / 4.0 + (1 - t.Gc) * (1 - t.Gc) * (1 - t.Gc) / 8.0;
        noStop = 1 - noStop;
        var last = new int[3]; var score = new double[3];
        for (int i = 0; i < 3; i++) score[i] = 0.0;
        for (int i = nn - 1; i >= 0; i--)
        {
            var n = nod[i]; int fr = n.Ndx % 3;
            if (n.Strand == 1 && n.Type == STOP) { last[fr] = n.Ndx; score[fr] = 0.0; }
            else if (n.Strand == 1)
            {
                for (int j = last[fr] - 3; j >= n.Ndx; j -= 3) score[fr] += t.GeneDc[MerNdx(6, q.F, j)];
                n.Cscore = score[fr]; last[fr] = n.Ndx;
            }
        }
        for (int i = 0; i < 3; i++) score[i] = 0.0;
        for (int i = 0; i < nn; i++)
        {
            var n = nod[i]; int fr = n.Ndx % 3;
            if (n.Strand == -1 && n.Type == STOP) { last[fr] = n.Ndx; score[fr] = 0.0; }
            else if (n.Strand == -1)
            {
                for (int j = last[fr] + 3; j <= n.Ndx; j += 3) score[fr] += t.GeneDc[MerNdx(6, q.R, slen - j - 1)];
                n.Cscore = score[fr]; last[fr] = n.Ndx;
            }
        }
        for (int i = 0; i < 3; i++) score[i] = -10000.0;
        for (int i = 0; i < nn; i++)
        {
            var n = nod[i]; int fr = n.Ndx % 3;
            if (n.Strand == 1 && n.Type == STOP) score[fr] = -10000.0;
            else if (n.Strand == 1)
            {
                if (n.Cscore > score[fr]) score[fr] = n.Cscore; else n.Cscore -= score[fr] - n.Cscore;
            }
        }
        for (int i = 0; i < 3; i++) score[i] = -10000.0;
        for (int i = nn - 1; i >= 0; i--)
        {
            var n = nod[i]; int fr = n.Ndx % 3;
            if (n.Strand == -1 && n.Type == STOP) score[fr] = -10000.0;
            else if (n.Strand == -1)
            {
                if (n.Cscore > score[fr]) score[fr] = n.Cscore; else n.Cscore -= score[fr] - n.Cscore;
            }
        }
        double Lfac(double gsize)
        {
            double lfac;
            if (gsize > 1000.0)
            {
                lfac = Math.Log((1 - Math.Pow(noStop, 1000.0)) / Math.Pow(noStop, 1000.0));
                lfac -= Math.Log((1 - Math.Pow(noStop, 80)) / Math.Pow(noStop, 80));
                lfac *= (gsize - 80) / 920.0;
            }
            else
            {
                lfac = Math.Log((1 - Math.Pow(noStop, gsize)) / Math.Pow(noStop, gsize));
                lfac -= Math.Log((1 - Math.Pow(noStop, 80)) / Math.Pow(noStop, 80));
            }
            return lfac;
        }
        for (int i = 0; i < nn; i++)
        {
            var n = nod[i]; int fr = n.Ndx % 3;
            if (n.Strand == 1 && n.Type == STOP) score[fr] = -10000.0;
            else if (n.Strand == 1) ApplyLength(n, fr);
        }
        for (int i = nn - 1; i >= 0; i--)
        {
            var n = nod[i]; int fr = n.Ndx % 3;
            if (n.Strand == -1 && n.Type == STOP) score[fr] = -10000.0;
            else if (n.Strand == -1) ApplyLength(n, fr);
        }
        void ApplyLength(Node n, int fr)
        {
            double gsize = (double)(float)(Math.Abs(n.StopVal - n.Ndx) + 3.0) / 3.0;
            double lfac = Lfac(gsize);
            if (lfac > score[fr]) score[fr] = lfac;
            else lfac -= Math.Max(Math.Min(score[fr] - lfac, lfac), 0);
            if (lfac > 3.0 && n.Cscore < 0.5 * lfac) n.Cscore = 0.5 * lfac;
            n.Cscore += lfac;
        }
    }

    // ------------------------------------------------------------------ RBS
    private static double[] ExactMatchVector(byte[] s, int pos, int limit)
    {
        var match = new double[6];
        for (int i = 0; i < 6; i++) match[i] = -10.0;
        for (int i = 0; i < limit; i++)
        {
            if (pos + i < 0) continue;
            byte b = Seq.B(s, pos + i);
            if (i % 3 == 0 && b == 0) match[i] = 2.0;
            else if (i % 3 != 0 && b == 1) match[i] = 3.0;
        }
        return match;
    }

    /// <summary>Maps (exact motif score, spacer class) to the exact-match bin; -1 when no row applies (Prodigal keeps its previous value).</summary>
    internal static int ExactBinOf(double ctr, int disFlag) => (ctr, disFlag) switch
    {
        (6.0, 2) => 1, (6.0, 3) => 2, (8.0, 3) => 3, (9.0, 3) => 3, (6.0, 1) => 6,
        (11.0, 3) => 10, (12.0, 3) => 10, (14.0, 3) => 10,
        (8.0, 2) => 11, (9.0, 2) => 11, (8.0, 1) => 12, (9.0, 1) => 12,
        (6.0, 0) => 13, (8.0, 0) => 15, (9.0, 0) => 16,
        (11.0, 2) => 20, (11.0, 1) => 21, (11.0, 0) => 22,
        (12.0, 2) => 20, (12.0, 1) => 23, (12.0, 0) => 24,
        (14.0, 2) => 25, (14.0, 1) => 26, (14.0, 0) => 27,
        _ => -1,
    };

    internal static int ExactSpacerClass(int rdis, int len)
    {
        if (rdis < 5 && len < 5) return 2;
        if (rdis < 5 && len >= 5) return 1;
        if (rdis > 10 && rdis <= 12 && len < 5) return 1;
        if (rdis > 10 && rdis <= 12 && len >= 5) return 2;
        if (rdis >= 13) return 3;
        return 0;
    }

    private static int ShineDalgarnoExact(byte[] s, int pos, int start, double[] rwt)
    {
        int limit = Math.Min(6, start - 4 - pos);
        var match = ExactMatchVector(s, pos, limit);
        int maxVal = 0, curVal = 0;
        for (int i = limit; i >= 3; i--)
        {
            for (int j = 0; j <= limit - i; j++)
            {
                double ctr = -2.0; int mism = 0;
                for (int k = j; k < j + i; k++) { ctr += match[k]; if (match[k] < 0.0) mism++; }
                if (mism > 0) continue;
                int rdis = start - (pos + j + i);
                int dis = ExactSpacerClass(rdis, i);
                if (rdis > 15 || ctr < 6.0) continue;
                int b = ExactBinOf(ctr, dis);
                if (b >= 0) curVal = b;
                if (rwt[curVal] < rwt[maxVal]) continue;
                if (rwt[curVal] == rwt[maxVal] && curVal < maxVal) continue;
                maxVal = curVal;
            }
        }
        return maxVal;
    }

    private static int ShineDalgarnoMm(byte[] s, int pos, int start, double[] rwt)
    {
        int limit = Math.Min(6, start - 4 - pos);
        var match = new double[6];
        for (int i = 0; i < 6; i++) match[i] = -10.0;
        for (int i = 0; i < limit; i++)
        {
            if (pos + i < 0) continue;
            byte b = Seq.B(s, pos + i);
            if (i % 3 == 0) match[i] = b == 0 ? 2.0 : -3.0;
            else match[i] = b == 1 ? 3.0 : -2.0;
        }
        int maxVal = 0, curVal = 0;
        for (int i = limit; i >= 5; i--)
        {
            for (int j = 0; j <= limit - i; j++)
            {
                double ctr = -2.0; int mism = 0;
                for (int k = j; k < j + i; k++)
                {
                    ctr += match[k];
                    if (match[k] < 0.0) mism++;
                    if (match[k] < 0.0 && (k <= j + 1 || k >= j + i - 2)) ctr -= 10.0;
                }
                if (mism != 1) continue;
                int rdis = start - (pos + j + i);
                int dis = 0;
                if (rdis < 5) dis = 1;
                else if (rdis > 10 && rdis <= 12) dis = 2;
                else if (rdis >= 13) dis = 3;
                if (rdis > 15 || ctr < 6.0) continue;
                int bin = (ctr, dis) switch
                {
                    (6.0, 3) => 2, (7.0, 3) => 2, (9.0, 3) => 3, (6.0, 2) => 4, (6.0, 1) => 5, (6.0, 0) => 9,
                    (7.0, 2) => 7, (7.0, 1) => 8, (7.0, 0) => 14, (9.0, 2) => 17, (9.0, 1) => 18, (9.0, 0) => 19,
                    _ => -1,
                };
                if (bin >= 0) curVal = bin;
                if (rwt[curVal] < rwt[maxVal]) continue;
                if (rwt[curVal] == rwt[maxVal] && curVal < maxVal) continue;
                maxVal = curVal;
            }
        }
        return maxVal;
    }

    private static void RbsScore(Seq q, List<Node> nod, Training t)
    {
        foreach (var n in nod)
        {
            if (n.Type == STOP || n.Edge == 1) continue;
            ScanRbs(q, n, t.RbsWt, out int e, out int m);
            n.Rbs[0] = e; n.Rbs[1] = m;
        }
    }

    private static void ScanRbs(Seq q, Node n, double[] w, out int exact, out int mm)
    {
        int slen = q.Len; exact = 0; mm = 0;
        if (n.Strand == 1)
        {
            for (int j = n.Ndx - 20; j <= n.Ndx - 6; j++)
            {
                if (j < 0) continue;
                int c0 = ShineDalgarnoExact(q.F, j, n.Ndx, w), c1 = ShineDalgarnoMm(q.F, j, n.Ndx, w);
                if (c0 > exact) exact = c0;
                if (c1 > mm) mm = c1;
            }
        }
        else
        {
            for (int j = slen - n.Ndx - 21; j <= slen - n.Ndx - 7; j++)
            {
                if (j > slen - 1) continue;
                int c0 = ShineDalgarnoExact(q.R, j, slen - 1 - n.Ndx, w), c1 = ShineDalgarnoMm(q.R, j, slen - 1 - n.Ndx, w);
                if (c0 > exact) exact = c0;
                if (c1 > mm) mm = c1;
            }
        }
    }

    private static int PickBin(double[] w, int e, int m)
    {
        if (w[e] > w[m] + 1.0 || m == 0) return e;
        if (w[e] < w[m] - 1.0 || e == 0) return m;
        return Math.Max(e, m);
    }

    internal static RbsStartClass ClassifyStart(string seq, int startIndex, int strand, RbsModel? model)
    {
        var q = Encode(seq);
        var w = model?.RbsWt ?? new double[28];
        var n = new Node { Ndx = startIndex, Strand = strand, Type = 0 };
        ScanRbs(q, n, w, out int e, out int m);
        int bin = PickBin(w, e, m);
        return new RbsStartClass(e, m, bin, w[bin]);
    }

    private static void CountUpstream(byte[] s, int slen, int strand, int pos, Training t)
    {
        int start = strand == 1 ? pos : slen - 1 - pos, count = 0;
        for (int i = 1; i < 45; i++)
        {
            if (i > 2 && i < 15) continue;
            if (start - i >= 0) t.UpsComp[count, MerNdx(1, s, start - i)]++;
            count++;
        }
    }

    private static void ScoreUpstream(byte[] s, int slen, Node n, Training t)
    {
        int start = n.Strand == 1 ? n.Ndx : slen - 1 - n.Ndx, count = 0;
        n.Uscore = 0.0;
        for (int i = 1; i < 45; i++)
        {
            if (i > 2 && i < 15) continue;
            if (start - i < 0) continue;
            n.Uscore += 0.4 * t.StWt * t.UpsComp[count, MerNdx(1, s, start - i)];
            count++;
        }
    }

    // returns (rreal, rbg, treal) of the last iteration
    private static (double[] rreal, double[] rbg, double[] treal) TrainStartsSd(Seq q, List<Node> nod, Training t)
    {
        int nn = nod.Count, slen = q.Len;
        double wt = t.StWt, sthresh = 35.0;
        Array.Clear(t.TypeWt); Array.Clear(t.RbsWt); Array.Clear(t.UpsComp);
        var tbg = new double[3];
        foreach (var n in nod) if (n.Type != STOP) tbg[n.Type] += 1.0;
        double sum = tbg.Sum();
        for (int i = 0; i < 3; i++) tbg[i] /= sum;
        double[] rbg = new double[28], rreal = new double[28], treal = new double[3];
        var best = new double[3]; var bndx = new int[3]; var rbs = new int[3]; var type = new int[3];
        for (int it = 0; it < 10; it++)
        {
            Array.Clear(rbg);
            foreach (var n in nod)
            {
                if (n.Type == STOP || n.Edge == 1) continue;
                rbg[PickBin(t.RbsWt, n.Rbs[0], n.Rbs[1])] += 1.0;
            }
            sum = rbg.Sum();
            for (int j = 0; j < 28; j++) rbg[j] /= sum;
            Array.Clear(rreal); Array.Clear(treal);

            for (int pass = 0; pass < 2; pass++)
            {
                int strand = pass == 0 ? 1 : -1;
                for (int j = 0; j < 3; j++) { best[j] = 0.0; bndx[j] = -1; rbs[j] = 0; type[j] = 0; }
                for (int jj = 0; jj < nn; jj++)
                {
                    int j = pass == 0 ? jj : nn - 1 - jj;
                    var n = nod[j];
                    if (n.Type != STOP && n.Edge == 1) continue;
                    int fr = n.Ndx % 3;
                    if (n.Type == STOP && n.Strand == strand)
                    {
                        if (best[fr] >= sthresh && nod[bndx[fr]].Ndx % 3 == fr)
                        {
                            rreal[rbs[fr]] += 1.0; treal[type[fr]] += 1.0;
                            if (it == 9) CountUpstream(strand == 1 ? q.F : q.R, slen, strand, nod[bndx[fr]].Ndx, t);
                        }
                        best[fr] = 0.0; bndx[fr] = -1; rbs[fr] = 0; type[fr] = 0;
                    }
                    else if (n.Strand == strand && n.Type != STOP)
                    {
                        int maxRb = PickBin(t.RbsWt, n.Rbs[0], n.Rbs[1]);
                        if (n.Cscore + wt * t.RbsWt[maxRb] + wt * t.TypeWt[n.Type] >= best[fr])
                        {
                            best[fr] = n.Cscore + wt * t.RbsWt[maxRb];
                            best[fr] += wt * t.TypeWt[n.Type];
                            bndx[fr] = j; type[fr] = n.Type; rbs[fr] = maxRb;
                        }
                    }
                }
            }
            sum = rreal.Sum();
            if (sum == 0.0) Array.Clear(t.RbsWt);
            else
                for (int j = 0; j < 28; j++)
                {
                    rreal[j] /= sum;
                    t.RbsWt[j] = rbg[j] != 0 ? Math.Clamp(Math.Log(rreal[j] / rbg[j]), -4.0, 4.0) : -4.0;
                }
            sum = treal.Sum();
            if (sum == 0.0) Array.Clear(t.TypeWt);
            else
                for (int j = 0; j < 3; j++)
                {
                    treal[j] /= sum;
                    t.TypeWt[j] = tbg[j] != 0 ? Math.Clamp(Math.Log(treal[j] / tbg[j]), -4.0, 4.0) : -4.0;
                }
            if (sum <= nn / 2000.0) sthresh /= 2.0;
        }
        for (int i = 0; i < 32; i++)
        {
            sum = 0.0;
            for (int j = 0; j < 4; j++) sum += t.UpsComp[i, j];
            if (sum == 0.0)
            {
                for (int j = 0; j < 4; j++) t.UpsComp[i, j] = 0.0;
                continue;
            }
            for (int j = 0; j < 4; j++)
            {
                t.UpsComp[i, j] /= sum;
                double v;
                bool at = j == 0 || j == 3;
                if (t.Gc > 0.1 && t.Gc < 0.9) v = Math.Log(t.UpsComp[i, j] * 2.0 / (at ? 1.0 - t.Gc : t.Gc));
                else if (t.Gc <= 0.1) v = Math.Log(t.UpsComp[i, j] * 2.0 / (at ? 0.90 : 0.10));
                else v = Math.Log(t.UpsComp[i, j] * 2.0 / (at ? 0.10 : 0.90));
                t.UpsComp[i, j] = Math.Clamp(v, -4.0, 4.0);
            }
        }
        return (rreal, rbg, treal);
    }

    private static bool DetermineSdUsage(Training t)
    {
        var w = t.RbsWt;
        bool uses = true;
        if (w[0] >= 0.0) uses = false;
        if (w[16] < 1.0 && w[13] < 1.0 && w[15] < 1.0 && (w[0] >= -0.5 || (w[22] < 2.0 && w[24] < 2.0 && w[27] < 2.0)))
            uses = false;
        return uses;
    }

    // ------------------------------------------------------------------ node scoring
    private static void ScoreNodes(Seq q, List<Node> nod, Training t)
    {
        int nn = nod.Count, slen = q.Len;
        RawCodingScore(q, nod, t);
        RbsScore(q, nod, t);
        for (int i = 0; i < nn; i++)
        {
            var n = nod[i];
            if (n.Type == STOP) continue;
            int edgeGene = 0;
            if (n.Edge == 1) edgeGene++;
            if ((n.Strand == 1 && !Seq.Stop(q.F, n.StopVal)) || (n.Strand == -1 && !Seq.Stop(q.R, slen - 1 - n.StopVal))) edgeGene++;

            if (n.Edge == 1)
            {
                n.Tscore = EdgeBonus * t.StWt / edgeGene; n.Uscore = 0.0; n.Rscore = 0.0;
            }
            else
            {
                n.Tscore = t.TypeWt[n.Type] * t.StWt;
                double sd = Math.Max(t.RbsWt[n.Rbs[0]], t.RbsWt[n.Rbs[1]]) * t.StWt;
                n.Rscore = sd;
                ScoreUpstream(n.Strand == 1 ? q.F : q.R, slen, n, t);
                if (n.Ndx <= 2 && n.Strand == 1) n.Uscore += EdgeUps * t.StWt;
                else if (n.Ndx >= slen - 3 && n.Strand == -1) n.Uscore += EdgeUps * t.StWt;
                else if (i < 500 && n.Strand == 1)
                {
                    for (int j = i - 1; j >= 0; j--)
                        if (nod[j].Edge == 1 && n.StopVal == nod[j].StopVal) { n.Uscore += EdgeUps * t.StWt; break; }
                }
                else if (i >= nn - 500 && n.Strand == -1)
                {
                    for (int j = i + 1; j < nn; j++)
                        if (nod[j].Edge == 1 && n.StopVal == nod[j].StopVal) { n.Uscore += EdgeUps * t.StWt; break; }
                }
            }

            if (((n.Ndx <= 2 && n.Strand == 1) || (n.Ndx >= slen - 3 && n.Strand == -1)) && n.Edge == 0)
            {
                edgeGene++; n.Edge = 1; n.Tscore = 0.0; n.Uscore = EdgeBonus * t.StWt / edgeGene; n.Rscore = 0.0;
            }
            if (n.Edge == 0 && edgeGene == 1) n.Uscore -= 0.5 * EdgeBonus * t.StWt;

            if (edgeGene == 0 && Math.Abs(n.Ndx - n.StopVal) < 250)
            {
                double negf = 250.0 / (float)Math.Abs(n.Ndx - n.StopVal);
                double posf = (float)Math.Abs(n.Ndx - n.StopVal) / 250.0;
                if (n.Rscore < 0) n.Rscore *= negf;
                if (n.Uscore < 0) n.Uscore *= negf;
                if (n.Tscore < 0) n.Tscore *= negf;
                if (n.Rscore > 0) n.Rscore *= posf;
                if (n.Uscore > 0) n.Uscore *= posf;
                if (n.Tscore > 0) n.Tscore *= posf;
            }
            n.Sscore = n.Tscore + n.Rscore + n.Uscore;
            if (n.Cscore < 0.0)
            {
                if (edgeGene > 0 && n.Edge == 0) n.Sscore -= t.StWt;
                else n.Sscore -= 0.5;
            }
        }
    }

    private static void EliminateBadGenes(List<Node> nod, int dbeg, Training t)
    {
        if (dbeg == -1) return;
        int path = dbeg;
        while (nod[path].Traceb != -1) path = nod[path].Traceb;
        while (nod[path].Tracef != -1)
        {
            var p = nod[path];
            if (p.Strand == 1 && p.Type == STOP) nod[p.Tracef].Sscore += IntergenicMod(p, nod[p.Tracef], t);
            if (p.Strand == -1 && p.Type != STOP) p.Sscore += IntergenicMod(p, nod[p.Tracef], t);
            path = p.Tracef;
        }
        path = dbeg;
        while (nod[path].Traceb != -1) path = nod[path].Traceb;
        while (nod[path].Tracef != -1)
        {
            var p = nod[path];
            if (p.Strand == 1 && p.Type != STOP && p.Cscore + p.Sscore < 0) { p.Elim = 1; nod[p.Tracef].Elim = 1; }
            if (p.Strand == -1 && p.Type == STOP && nod[p.Tracef].Cscore + nod[p.Tracef].Sscore < 0) { p.Elim = 1; nod[p.Tracef].Elim = 1; }
            path = p.Tracef;
        }
    }

    internal sealed class PGene { public int StartNdx, StopNdx; }

    private static List<PGene> AddGenes(List<Node> nod, int dbeg)
    {
        var genes = new List<PGene>();
        if (dbeg == -1) return genes;
        int path = dbeg;
        while (nod[path].Traceb != -1) path = nod[path].Traceb;
        PGene cur = new();
        while (path != -1)
        {
            var n = nod[path];
            if (n.Elim == 1) { path = n.Tracef; continue; }
            if (n.Strand == 1 && n.Type != STOP) cur.StartNdx = path;
            if (n.Strand == -1 && n.Type == STOP) cur.StopNdx = path;
            if (n.Strand == 1 && n.Type == STOP) { cur.StopNdx = path; genes.Add(cur); cur = new(); }
            if (n.Strand == -1 && n.Type != STOP) { cur.StartNdx = path; genes.Add(cur); cur = new(); }
            path = n.Tracef;
        }
        return genes;
    }

    private static void TweakFinalStarts(List<PGene> genes, List<Node> nod, Training t)
    {
        int ng = genes.Count, nn = nod.Count;
        for (int i = 0; i < ng; i++)
        {
            int ndx = genes[i].StartNdx; var cn = nod[ndx];
            double sc = cn.Sscore + cn.Cscore, igm = 0.0;
            if (i > 0 && cn.Strand == 1 && nod[genes[i - 1].StartNdx].Strand == 1) igm = IntergenicMod(nod[genes[i - 1].StopNdx], cn, t);
            if (i > 0 && cn.Strand == 1 && nod[genes[i - 1].StartNdx].Strand == -1) igm = IntergenicMod(nod[genes[i - 1].StartNdx], cn, t);
            if (i < ng - 1 && cn.Strand == -1 && nod[genes[i + 1].StartNdx].Strand == 1) igm = IntergenicMod(cn, nod[genes[i + 1].StartNdx], t);
            if (i < ng - 1 && cn.Strand == -1 && nod[genes[i + 1].StartNdx].Strand == -1) igm = IntergenicMod(cn, nod[genes[i + 1].StopNdx], t);

            int[] maxndx = { -1, -1 }; double[] maxsc = { 0, 0 }, maxigm = { 0, 0 };
            for (int j = ndx - 100; j < ndx + 100; j++)
            {
                if (j < 0 || j >= nn || j == ndx) continue;
                var nj = nod[j];
                if (nj.Type == STOP || nj.StopVal != cn.StopVal) continue;
                double tigm = 0.0;
                if (i > 0 && nj.Strand == 1 && nod[genes[i - 1].StartNdx].Strand == 1)
                {
                    if (nod[genes[i - 1].StopNdx].Ndx - nj.Ndx > MaxSamOvlp) continue;
                    tigm = IntergenicMod(nod[genes[i - 1].StopNdx], nj, t);
                }
                if (i > 0 && nj.Strand == 1 && nod[genes[i - 1].StartNdx].Strand == -1)
                {
                    if (nod[genes[i - 1].StartNdx].Ndx - nj.Ndx >= 0) continue;
                    tigm = IntergenicMod(nod[genes[i - 1].StartNdx], nj, t);
                }
                if (i < ng - 1 && nj.Strand == -1 && nod[genes[i + 1].StartNdx].Strand == 1)
                {
                    if (nj.Ndx - nod[genes[i + 1].StartNdx].Ndx >= 0) continue;
                    tigm = IntergenicMod(nj, nod[genes[i + 1].StartNdx], t);
                }
                if (i < ng - 1 && nj.Strand == -1 && nod[genes[i + 1].StartNdx].Strand == -1)
                {
                    if (nj.Ndx - nod[genes[i + 1].StopNdx].Ndx > MaxSamOvlp) continue;
                    tigm = IntergenicMod(nj, nod[genes[i + 1].StopNdx], t);
                }
                if (maxndx[0] == -1) { maxndx[0] = j; maxsc[0] = nj.Cscore + nj.Sscore; maxigm[0] = tigm; }
                else if (nj.Cscore + nj.Sscore + tigm > maxsc[0])
                {
                    maxndx[1] = maxndx[0]; maxsc[1] = maxsc[0]; maxigm[1] = maxigm[0];
                    maxndx[0] = j; maxsc[0] = nj.Cscore + nj.Sscore; maxigm[0] = tigm;
                }
                else if (maxndx[1] == -1 || nj.Cscore + nj.Sscore + tigm > maxsc[1])
                {
                    maxndx[1] = j; maxsc[1] = nj.Cscore + nj.Sscore; maxigm[1] = tigm;
                }
            }
            for (int j = 0; j < 2; j++)
            {
                int mndx = maxndx[j]; if (mndx == -1) continue;
                var m = nod[mndx];
                if (m.Tscore < cn.Tscore && maxsc[j] - m.Tscore >= sc - cn.Tscore + t.StWt && m.Rscore > cn.Rscore &&
                    m.Uscore > cn.Uscore && m.Cscore > cn.Cscore && Math.Abs(m.Ndx - cn.Ndx) > 15)
                    maxsc[j] += cn.Tscore - m.Tscore;
                else if (Math.Abs(m.Ndx - cn.Ndx) <= 15 && m.Rscore + m.Tscore > cn.Rscore + cn.Tscore && cn.Edge == 0 && m.Edge == 0)
                {
                    if (cn.Cscore > m.Cscore) maxsc[j] += cn.Cscore - m.Cscore;
                    if (cn.Uscore > m.Uscore) maxsc[j] += cn.Uscore - m.Uscore;
                    if (igm > maxigm[j]) maxsc[j] += igm - maxigm[j];
                }
                else maxsc[j] = -1000.0;
            }
            int sel = -1;
            for (int j = 0; j < 2; j++)
            {
                if (maxndx[j] == -1) continue;
                if (sel == -1 && maxsc[j] + maxigm[j] > sc + igm) sel = j;
                else if (sel >= 0 && maxsc[j] + maxigm[j] > maxsc[sel] + maxigm[sel]) sel = j;
            }
            if (sel != -1) genes[i].StartNdx = maxndx[sel];
        }
    }

    // ------------------------------------------------------------------ drivers
    private sealed class Trained
    {
        public required Seq Q;
        public required Training T;
        public required List<Node> Nodes;
        public required double[] RReal, Rbg, TReal;
    }

    private static Trained TrainCore(string genome)
    {
        var q = Encode(genome);
        var t = new Training();
        int gcCount = 0;
        for (int i = 0; i < q.Len; i++) if (q.F[i] == 1 || q.F[i] == 2) gcCount++;
        t.Gc = (double)gcCount / q.Len;
        var nodes = AddNodes(q, closed: false);
        var gcf = GcFramePlot(q.F);
        RecordGcBias(gcf, nodes, t);
        RecordOverlappingStarts(nodes, t, flag: false);
        int ipath = DProg(nodes, t, flag: false);
        CalcDicodonGene(t, q, nodes, ipath);
        RawCodingScore(q, nodes, t);
        RbsScore(q, nodes, t);
        var (rreal, rbg, treal) = TrainStartsSd(q, nodes, t);
        t.UsesSd = DetermineSdUsage(t);
        return new Trained { Q = q, T = t, Nodes = nodes, RReal = rreal, Rbg = rbg, TReal = treal };
    }

    internal static RbsModel Train(string genome)
    {
        var tr = TrainCore(genome);
        if (!tr.T.UsesSd)
            throw new NotSupportedException("The genome does not use Shine-Dalgarno motifs (Prodigal's non-SD motif trainer is not ported).");
        return ToModel(tr.T, tr.RReal, tr.Rbg, tr.TReal);
    }

    private static RbsModel ToModel(Training t, double[] rreal, double[] rbg, double[] treal) =>
        new(t.Gc, t.Bias, t.TypeWt, t.RbsWt, t.UpsComp, t.GeneDc, t.UsesSd, rreal, rbg, treal);

    private static Training FromModel(RbsModel m) => new()
    {
        Gc = m.GcContent, Bias = m.Bias, TypeWt = m.TypeWt, RbsWt = m.RbsWt, UpsComp = m.UpsComp,
        GeneDc = m.GeneDc ?? new double[4096], UsesSd = m.UsesShineDalgarno,
    };

    internal static RbsModel TrainFromGenes(string genome, IEnumerable<GenomeAnnotator.GeneAnnotation> genes)
    {
        var q = Encode(genome);
        int slen = q.Len, gcCount = 0;
        for (int i = 0; i < slen; i++) if (q.F[i] == 1 || q.F[i] == 2) gcCount++;
        var t = new Training { Gc = (double)gcCount / slen };
        var nodes = AddNodes(q, closed: false);
        RbsScore(q, nodes, t);
        var byKey = new Dictionary<(int, int), Node>();
        foreach (var n in nodes) if (n.Type != STOP) byKey[(n.Strand, n.Ndx)] = n;
        var starts = new List<Node>();
        foreach (var g in genes)
        {
            int strand = g.Strand == '-' ? -1 : 1;
            int ndx = strand == 1 ? g.Start : g.End - 1;
            if (byKey.TryGetValue((strand, ndx), out var nd)) starts.Add(nd);
        }
        var tbg = new double[3];
        foreach (var n in nodes) if (n.Type != STOP) tbg[n.Type] += 1.0;
        double sum = tbg.Sum(); for (int i = 0; i < 3; i++) tbg[i] /= sum;
        double[] rreal = new double[28], rbg = new double[28], treal = new double[3];
        for (int it = 0; it < 10; it++)
        {
            Array.Clear(rbg); Array.Clear(rreal); Array.Clear(treal);
            foreach (var n in nodes)
            {
                if (n.Type == STOP || n.Edge == 1) continue;
                rbg[PickBin(t.RbsWt, n.Rbs[0], n.Rbs[1])] += 1.0;
            }
            sum = rbg.Sum(); for (int j = 0; j < 28; j++) rbg[j] /= sum;
            foreach (var n in starts)
            {
                if (n.Edge == 1) continue;
                rreal[PickBin(t.RbsWt, n.Rbs[0], n.Rbs[1])] += 1.0; treal[n.Type] += 1.0;
            }
            sum = rreal.Sum();
            if (sum == 0.0) Array.Clear(t.RbsWt);
            else
            {
                for (int j = 0; j < 28; j++)
                {
                    rreal[j] /= sum;
                    t.RbsWt[j] = rbg[j] != 0 ? Math.Clamp(Math.Log(rreal[j] / rbg[j]), -4.0, 4.0) : -4.0;
                }
            }
            sum = treal.Sum();
            if (sum == 0.0) Array.Clear(t.TypeWt);
            else
            {
                for (int j = 0; j < 3; j++)
                {
                    treal[j] /= sum;
                    t.TypeWt[j] = tbg[j] != 0 ? Math.Clamp(Math.Log(treal[j] / tbg[j]), -4.0, 4.0) : -4.0;
                }
            }
        }
        t.UsesSd = DetermineSdUsage(t);
        return new RbsModel(t.Gc, new double[3], t.TypeWt, t.RbsWt, new double[32, 4], null, t.UsesSd, rreal, rbg, treal);
    }

    /// <summary>Scored genome: every start node with Prodigal scores under a trained model.</summary>
    internal sealed class ScoredGenome
    {
        public required Seq Q;
        public required Training T;
        public required List<Node> Nodes;
        public required Dictionary<(int strand, int ndx), Node> StartByKey;
    }

    internal static ScoredGenome ScoreGenome(string genome, RbsModel? model, out RbsModel usedModel)
    {
        Trained? tr = null;
        Training t;
        if (model == null || model.GeneDc == null)
        {
            if (genome.Length < 20000)
                throw new ArgumentException("Prodigal start selection trains on the input and needs at least 20000 nt (or pass a model trained on a larger genome).", nameof(genome));
            tr = TrainCore(genome);
            if (!tr.T.UsesSd) throw new NotSupportedException("The genome does not use Shine-Dalgarno motifs (Prodigal's non-SD motif trainer is not ported).");
            t = tr.T;
            if (model != null)
            {   // keep the caller's RBS/type weights but Prodigal's coding model of this genome
                Array.Copy(model.RbsWt, t.RbsWt, 28); Array.Copy(model.TypeWt, t.TypeWt, 3);
            }
            usedModel = model ?? ToModel(t, tr.RReal, tr.Rbg, tr.TReal);
        }
        else { t = FromModel(model); usedModel = model; }
        var q = tr?.Q ?? Encode(genome);
        var nodes = AddNodes(q, closed: false);
        ScoreNodes(q, nodes, t);
        var map = new Dictionary<(int, int), Node>();
        foreach (var n in nodes) if (n.Type != STOP) map[(n.Strand, n.Ndx)] = n;
        return new ScoredGenome { Q = q, T = t, Nodes = nodes, StartByKey = map };
    }

    internal readonly record struct DpGene(int Start, int End, char Strand, int StartType, double Cscore, double Sscore, double Rscore, double Tscore, double Uscore, int RbsBin);

    internal static List<DpGene> RunDp(string genome, RbsModel? model)
    {
        var sg = ScoreGenome(genome, model, out _);
        RecordOverlappingStarts(sg.Nodes, sg.T, flag: true);
        int ipath = DProg(sg.Nodes, sg.T, flag: true);
        EliminateBadGenes(sg.Nodes, ipath, sg.T);
        var genes = AddGenes(sg.Nodes, ipath);
        TweakFinalStarts(genes, sg.Nodes, sg.T);
        var res = new List<DpGene>();
        foreach (var g in genes)
        {
            var s = sg.Nodes[g.StartNdx]; var e = sg.Nodes[g.StopNdx];
            int bin = PickBin(sg.T.RbsWt, s.Rbs[0], s.Rbs[1]);
            if (s.Strand == 1) res.Add(new DpGene(s.Ndx, e.Ndx + 2, '+', s.Type, s.Cscore, s.Sscore, s.Rscore, s.Tscore, s.Uscore, bin));
            else res.Add(new DpGene(e.Ndx - 2, s.Ndx, '-', s.Type, s.Cscore, s.Sscore, s.Rscore, s.Tscore, s.Uscore, bin));
        }
        return res;
    }

    /// <summary>Returns the Prodigal start score (cscore + sscore) of the start at (strand, ndx), or null if it is not a candidate node.</summary>
    internal static double? StartScore(ScoredGenome sg, int strand, int ndx) =>
        sg.StartByKey.TryGetValue((strand, ndx), out var n) ? n.Cscore + n.Sscore : null;
}
