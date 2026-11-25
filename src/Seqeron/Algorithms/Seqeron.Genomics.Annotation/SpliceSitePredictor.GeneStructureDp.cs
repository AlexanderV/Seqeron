using Seqeron.Genomics.Core;

namespace Seqeron.Genomics.Annotation;

public static partial class SpliceSitePredictor
{
    #region Gene Structure Prediction (splice-site + open-reading-frame dynamic program)

    /// <summary>
    /// Result of <see cref="PredictGeneStructureDp"/>.
    /// </summary>
    /// <param name="Structure">Exons (tiling the input), selected introns, spliced sequence and
    /// <see cref="GeneStructure.OverallScore"/> = <see cref="IntronScore"/> + <see cref="CodingScore"/>.</param>
    /// <param name="CdsStart">0-based genomic index of the first coding nucleotide (the start codon).</param>
    /// <param name="CodingSequence">Spliced coding sequence from the start codon to the first in-frame
    /// stop codon (stop included), RNA alphabet.</param>
    /// <param name="IntronScore">Sum over the selected introns of MaxEntScan score5 + score3 - intronPenalty (bits).</param>
    /// <param name="CodingScore">Sum of the optional user-supplied exon scores (0 without a scorer).</param>
    public readonly record struct GeneStructureDpResult(
        GeneStructure Structure,
        int CdsStart,
        string CodingSequence,
        double IntronScore,
        double CodingScore);

    // Per-intron cost (bits) of the objective. Without it the sum of positive log-odds adds every
    // frame-preserving GT..AG pair (CFTR exons 2-5: 121 predicted introns for 3 annotated). No published
    // value exists (GENSCAN / GeneSplicer use trained length and transition priors); 12 bits was chosen on
    // 12 annotated human genes that are not among the evaluation genes (ACTB, GAPDH, HRAS, INS, GH1, TPI1,
    // LDHA, SOD1, APOE, TTR, CYCS, GCK, ALB; 32-38 of 65 introns recovered at 12 bits versus 15-30 at 0-10, see Evidence).
    private const double DefaultDpIntronPenalty = 12.0;

    /// <summary>
    /// Predicts the exon/intron structure of a gene as the globally optimal chain of introns under a
    /// splice-site log-odds score and an open-reading-frame constraint (dynamic program).
    /// </summary>
    /// <remarks>
    /// <para>
    /// <b>Model</b> (in the spirit of GENSCAN, Burge &amp; Karlin 1997 J Mol Biol 268:78, and
    /// GeneSplicer, Pertea, Lin &amp; Salzberg 2001 Nucleic Acids Res 29:1185): the structure is
    /// 5' UTR | initial exon starting at the start codon | (intron, internal exon)* | final exon holding
    /// the stop codon | 3' UTR. Splicing removes introns only (Gilbert 1978), the coding frame
    /// carries over every junction (the phase states 0/1/2 of GENSCAN: a codon may be split across an
    /// intron, also the stop codon), and the spliced coding sequence must be an open reading frame
    /// from the start codon to the first in-frame stop with no premature stop. UTRs are modelled as
    /// intron-free, as in GENSCAN's single-gene model.
    /// </para>
    /// <para>
    /// <b>Objective.</b> Maximise the sum over the selected introns of (the MaxEntScan (Yeo &amp; Burge
    /// 2004) 5' score (<see cref="ScoreDonorMaxEnt"/>, window -3..+6 around a GT or GC) plus 3' score
    /// (<see cref="ScoreAcceptorMaxEnt"/>, window -20..+3 around an AG) - <paramref name="intronPenalty"/>);
    /// both scores are log2 likelihood ratios (bits), the sum is the log-odds of the donor and acceptor
    /// being true sites, the penalty the prior cost of an intron (unpublished, tunable). Optionally
    /// a user-supplied exon coding score is added (e.g. a CPAT hexamer log-likelihood, see
    /// <see cref="CpatModel"/>). The default coding term is the open-reading-frame constraint only.
    /// The optimum is exact (Bellman recursion over (acceptor, partial codon) states with sliding-window
    /// maxima over donors, O(acceptors x donors-in-reading-frame-window) time, O(acceptors) memory); it
    /// equals exhaustive search over all intron sets.
    /// </para>
    /// <para>
    /// <b>Constraints.</b> Introns are GT-AG or GC-AG (U2-type), non-overlapping, of length
    /// <paramref name="minIntronLength"/>..<paramref name="maxIntronLength"/> (at least 4 nt); every
    /// exon (the first one counted from the start codon, the last one up to the end of the input) has at
    /// least <paramref name="minExonLength"/> nt, effectively at least 2 so that a codon spans at most
    /// two exons. No intron lies in the 5' UTR (before the start codon) or the 3' UTR (after the stop
    /// codon). Donor / acceptor sites below <paramref name="minDonorScore"/> /
    /// <paramref name="minAcceptorScore"/> bits are not considered. U12-type introns are not modelled.
    /// </para>
    /// <para>
    /// The start codon is <paramref name="startPosition"/> when given (any 0-based position, e.g. the
    /// first codon of a window into a gene) and otherwise the first AUG of the sequence. Exons tile the
    /// whole input like <see cref="PredictGeneStructure"/>; <see cref="Exon.Phase"/> is the number of
    /// coding nucleotides in the preceding exons modulo 3. <see cref="Intron.Score"/> is donor + acceptor
    /// bits; <see cref="Intron.BranchPoint"/> is not used by the model and is null. Returns null when no
    /// start exists or no open reading frame is reachable.
    /// </para>
    /// </remarks>
    /// <param name="sequence">Genomic sequence of the gene (DNA or RNA).</param>
    /// <param name="startPosition">0-based first coding nucleotide; null = first AUG.</param>
    /// <param name="minExonLength">Minimum exon length (nt), effectively at least 2.</param>
    /// <param name="minIntronLength">Minimum intron length (nt), effectively at least 4.</param>
    /// <param name="maxIntronLength">Maximum intron length (nt).</param>
    /// <param name="minDonorScore">Minimum MaxEntScan 5' score (bits) of a candidate donor.</param>
    /// <param name="minAcceptorScore">Minimum MaxEntScan 3' score (bits) of a candidate acceptor.</param>
    /// <param name="exonCodingScore">Optional score of a coding exon: (exon sequence in DNA alphabet,
    /// offset of its first complete codon, 0..2) -> score added to the objective.</param>
    /// <param name="intronPenalty">Bits subtracted per intron (the prior cost of an intron, see remarks);
    /// default 12. 0 maximises the plain sum of log-odds and over-predicts introns.</param>
    public static GeneStructureDpResult? PredictGeneStructureDp(
        string sequence,
        int? startPosition = null,
        int minExonLength = 3,
        int minIntronLength = 60,
        int maxIntronLength = 100000,
        double minDonorScore = 0.0,
        double minAcceptorScore = 0.0,
        Func<string, int, double>? exonCodingScore = null,
        double intronPenalty = DefaultDpIntronPenalty)
    {
        if (minExonLength < 1)
            throw new ArgumentOutOfRangeException(nameof(minExonLength), minExonLength, "Must be at least 1.");
        if (minIntronLength < 1)
            throw new ArgumentOutOfRangeException(nameof(minIntronLength), minIntronLength, "Must be at least 1.");
        if (maxIntronLength < minIntronLength)
            throw new ArgumentOutOfRangeException(nameof(maxIntronLength), maxIntronLength, "Must be >= minIntronLength.");
        if (string.IsNullOrEmpty(sequence))
            return null;
        if (startPosition is int sp && (sp < 0 || sp + 3 > sequence.Length))
            throw new ArgumentOutOfRangeException(nameof(startPosition), sp, "Start codon must lie inside the sequence.");

        string upper = sequence.ToUpperInvariant().Replace('T', 'U');
        int start = startPosition ?? upper.IndexOf("AUG", StringComparison.Ordinal);
        if (start < 0)
            return null;

        var donorSites = FindDonorSites(upper, minDonorScore, true, DonorScoringModel.MaxEntScan)
            .Where(s => s.Position > start).OrderBy(s => s.Position).ToList();
        var acceptorSites = FindAcceptorSites(upper, minAcceptorScore, false, AcceptorScoringModel.MaxEntScan)
            .Where(s => s.Position > start).OrderBy(s => s.Position).ToList();

        var solution = SolveSpliceDp(upper, start,
            donorSites.Select(s => (s.Position, s.Score)).ToList(),
            acceptorSites.Select(s => (s.Position, s.Score - intronPenalty)).ToList(),
            minExonLength, minIntronLength, maxIntronLength, exonCodingScore);
        if (solution is null)
            return null;

        var donorByPos = donorSites.ToDictionary(s => s.Position);
        var acceptorByPos = acceptorSites.ToDictionary(s => s.Position);
        var introns = solution.Value.Introns.Select(p =>
        {
            int length = p.Acceptor - p.Donor + 1;
            var donor = donorByPos[p.Donor];
            var acceptor = acceptorByPos[p.Acceptor];
            return new Intron(p.Donor, p.Acceptor, length, donor, acceptor, null,
                upper.Substring(p.Donor, length),
                upper[p.Donor + 1] == 'C' ? IntronType.GcAg : IntronType.U2,
                donor.Score + acceptor.Score);
        }).ToList();

        var tiled = DeriveExons(upper, introns);
        var exons = new List<Exon>(tiled.Count);
        int codingBefore = 0;
        for (int k = 0; k < tiled.Count; k++)
        {
            exons.Add(tiled[k] with { Phase = codingBefore % 3 });
            codingBefore += k == 0 ? tiled[k].End - start + 1 : tiled[k].Length;
        }

        string spliced = GenerateSplicedSequence(exons);
        var code = GeneticCode.Standard;
        int stopEnd = -1;
        for (int i = start; i + 3 <= spliced.Length; i += 3)
        {
            if (code.IsStopCodon(spliced.Substring(i, 3)))
            {
                stopEnd = i + 3;
                break;
            }
        }

        var structure = new GeneStructure(exons, introns, spliced, solution.Value.Score);
        return new GeneStructureDpResult(structure, start, spliced.Substring(start, stopEnd - start),
            solution.Value.IntronScore, solution.Value.CodingScore);
    }

    internal readonly record struct SpliceDpSolution(
        double Score, double IntronScore, double CodingScore, List<(int Donor, int Acceptor)> Introns);

    // Partial-codon states carried across an intron: 0 = none, 1..5 = one nucleotide (A,C,G,U,other),
    // 6..30 = two nucleotides. Nucleotide codes 0..4 (4 = anything but A/C/G/U).
    private const int PartialStates = 31;
    private const int NoStop = int.MaxValue;

    private static readonly bool[] StopCodonTable = BuildStopCodonTable();

    private static bool[] BuildStopCodonTable()
    {
        const string bases = "ACGU";
        var table = new bool[125];
        for (int x = 0; x < 4; x++)
            for (int y = 0; y < 4; y++)
                for (int z = 0; z < 4; z++)
                    table[25 * x + 5 * y + z] =
                        GeneticCode.Standard.IsStopCodon(string.Concat(bases[x], bases[y], bases[z]));
        return table;
    }

    private static int BaseCode(char c) => c switch { 'A' => 0, 'C' => 1, 'G' => 2, 'U' => 3, _ => 4 };

    private static int PartialPhase(int p)
    {
        if (p == 0) return 0;
        return p < 6 ? 1 : 2;
    }

    // State of the 'phase' coding nucleotides immediately before genomic index pos.
    private static int PartialCode(byte[] code, int pos, int phase)
    {
        if (phase == 0) return 0;
        if (phase == 1) return 1 + code[pos - 1];
        return 6 + 5 * code[pos - 2] + code[pos - 1];
    }

    /// <summary>
    /// Exact optimum of the splice + open-reading-frame model over the given donor / acceptor
    /// candidates (positions = first / last intron nucleotide, ascending, scores in bits). See
    /// <see cref="PredictGeneStructureDp"/> for the constraints. Returns null when no open reading
    /// frame from <paramref name="start"/> to a stop exists under the constraints.
    /// </summary>
    internal static SpliceDpSolution? SolveSpliceDp(
        string upper,
        int start,
        IReadOnlyList<(int Position, double Score)> donors,
        IReadOnlyList<(int Position, double Score)> acceptors,
        int minExonLength,
        int minIntronLength,
        int maxIntronLength,
        Func<string, int, double>? exonCodingScore)
    {
        int len = upper.Length;
        int minExon = Math.Max(2, minExonLength);
        int minIntron = Math.Max(4, minIntronLength);
        int nD = donors.Count, nA = acceptors.Count;
        string dna = exonCodingScore is null ? "" : upper.Replace('U', 'T');

        var code = new byte[len];
        for (int i = 0; i < len; i++) code[i] = (byte)BaseCode(upper[i]);

        // stopNext[x] = start of the first in-frame (x, x+3, ...) stop codon at or after x.
        var stopNext = new int[len + 3];
        Array.Fill(stopNext, NoStop);
        for (int x = len - 1; x >= 0; x--)
        {
            bool stop = x + 2 < len && StopCodonTable[25 * code[x] + 5 * code[x + 1] + code[x + 2]];
            stopNext[x] = stop ? x : stopNext[x + 3];
        }

        var dPos = new int[nD];
        var dScore = new double[nD];
        for (int j = 0; j < nD; j++) (dPos[j], dScore[j]) = donors[j];
        var aPos = new int[nA];
        var aScore = new double[nA];
        for (int i = 0; i < nA; i++) (aPos[i], aScore[i]) = acceptors[i];

        // G[j*3+ph]: best score of a prefix whose last exon ends right before donor j with ph
        // leftover coding nucleotides; back pointers (acceptor index or -1 = start, partial state).
        var gVal = new double[nD * 3];
        Array.Fill(gVal, double.NegativeInfinity);
        var gBackA = new int[nD * 3];
        var gBackP = new int[nD * 3];
        // F[i*31+p]: best score of a prefix ending with an intron at acceptor i, partial state p.
        var fVal = new double[nA * PartialStates];
        Array.Fill(fVal, double.NegativeInfinity);
        var fBack = new int[nA * PartialStates];

        double bestTotal = double.NegativeInfinity;
        int bestI = -2, bestP = 0;

        double Coding(int from, int toExclusive, int offset) =>
            exonCodingScore is null ? 0.0 : exonCodingScore(dna.Substring(from, toExclusive - from), offset);

        // Handles one acceptor (or the virtual start, accIdx = -1) whose F values are final: evaluates
        // the terminal exon and relaxes all donors reachable through one exon.
        void Expand(int accIdx, int a, ReadOnlySpan<double> f, int firstExonMin)
        {
            var v = new double[3];
            var vp = new int[3];
            Array.Fill(v, double.NegativeInfinity);
            for (int p = 0; p < PartialStates; p++)
            {
                double fp = f[p];
                if (double.IsNegativeInfinity(fp)) continue;
                int ph = PartialPhase(p);
                bool junctionStop = false;
                if (ph == 1)
                {
                    if (a + 2 >= len) continue;
                    junctionStop = StopCodonTable[25 * (p - 1) + 5 * code[a + 1] + code[a + 2]];
                }
                else if (ph == 2)
                {
                    if (a + 1 >= len) continue;
                    int q = p - 6;
                    junctionStop = StopCodonTable[25 * (q / 5) + 5 * (q % 5) + code[a + 1]];
                }

                int c0 = ph == 0 ? a + 1 : a + 1 + (3 - ph);
                if (len - 1 - a >= firstExonMin)
                {
                    // Terminal exon: the first in-frame stop must exist; for a junction stop it is the junction codon.
                    int stopEnd = -1;
                    if (junctionStop) stopEnd = c0 - 1;
                    else if (stopNext[c0] != NoStop) stopEnd = stopNext[c0] + 2;
                    if (stopEnd >= 0)
                    {
                        int from = accIdx < 0 ? start : a + 1;
                        double coding = Coding(from, stopEnd + 1, ph == 0 ? 0 : 3 - ph);
                        double total = fp + coding;
                        if (total > bestTotal)
                        {
                            bestTotal = total; bestI = accIdx; bestP = p;
                        }
                    }
                }

                if (!junctionStop && fp > v[ph])
                {
                    v[ph] = fp;
                    vp[ph] = p;
                }
            }

            for (int ph = 0; ph < 3; ph++)
            {
                if (double.IsNegativeInfinity(v[ph])) continue;
                int c0 = ph == 0 ? a + 1 : a + 1 + (3 - ph);
                int lo = a + 1 + firstExonMin;
                int qStop = stopNext[c0];
                int hi = qStop == NoStop ? len - 1 : qStop + 2;
                int idx = LowerBound(dPos, lo);
                for (; idx < nD && dPos[idx] <= hi; idx++)
                {
                    int d = dPos[idx];
                    int phNew = (d - c0) % 3;
                    double cand = v[ph];
                    if (exonCodingScore is not null)
                        cand += Coding(accIdx < 0 ? start : a + 1, d, ph == 0 ? 0 : 3 - ph);
                    int g = idx * 3 + phNew;
                    if (cand > gVal[g])
                    {
                        gVal[g] = cand;
                        gBackA[g] = accIdx;
                        gBackP[g] = vp[ph];
                    }
                }
            }
        }

        // Start state: coding begins at 'start'; the first exon must hold at least the start codon.
        var startF = new double[PartialStates];
        Array.Fill(startF, double.NegativeInfinity);
        startF[0] = 0.0;
        Expand(-1, start - 1, startF, Math.Max(3, minExon));

        // One max-queue (monotonic deque) per partial state over the donors inside the intron-length window.
        var dqVal = new List<double>[PartialStates];
        var dqIdx = new List<int>[PartialStates];
        var dqHead = new int[PartialStates];
        for (int p = 0; p < PartialStates; p++) { dqVal[p] = new List<double>(); dqIdx[p] = new List<int>(); }
        int inserted = 0;

        for (int i = 0; i < nA; i++)
        {
            int a = aPos[i];
            while (inserted < nD && dPos[inserted] <= a - minIntron + 1)
            {
                int j = inserted++;
                for (int ph = 0; ph < 3; ph++)
                {
                    double gv = gVal[j * 3 + ph];
                    if (double.IsNegativeInfinity(gv) || dPos[j] - ph < 0) continue;
                    int p = PartialCode(code, dPos[j], ph);
                    double h = gv + dScore[j];
                    var vals = dqVal[p];
                    var idxs = dqIdx[p];
                    while (vals.Count > dqHead[p] && vals[^1] < h)
                    {
                        vals.RemoveAt(vals.Count - 1);
                        idxs.RemoveAt(idxs.Count - 1);
                    }
                    vals.Add(h);
                    idxs.Add(j);
                }
            }

            int windowStart = a - maxIntronLength + 1;
            bool any = false;
            for (int p = 0; p < PartialStates; p++)
            {
                var idxs = dqIdx[p];
                while (idxs.Count > dqHead[p] && dPos[idxs[dqHead[p]]] < windowStart)
                    dqHead[p]++;
                if (idxs.Count > dqHead[p])
                {
                    fVal[i * PartialStates + p] = aScore[i] + dqVal[p][dqHead[p]];
                    fBack[i * PartialStates + p] = idxs[dqHead[p]];
                    any = true;
                }
            }

            if (any)
                Expand(i, a, new ReadOnlySpan<double>(fVal, i * PartialStates, PartialStates), minExon);
        }

        if (bestI == -2)
            return null;

        var introns = new List<(int, int)>();
        int cur = bestI, cp = bestP;
        double intronScore = 0;
        while (cur >= 0)
        {
            int j = fBack[cur * PartialStates + cp];
            introns.Add((dPos[j], aPos[cur]));
            intronScore += dScore[j] + aScore[cur];
            int g = j * 3 + PartialPhase(cp);
            cur = gBackA[g];
            cp = gBackP[g];
        }

        introns.Reverse();
        return new SpliceDpSolution(bestTotal, intronScore, exonCodingScore is null ? 0.0 : bestTotal - intronScore, introns);
    }

    private static int LowerBound(int[] sorted, int value)
    {
        int lo = 0, hi = sorted.Length;
        while (lo < hi)
        {
            int mid = (lo + hi) >>> 1;
            if (sorted[mid] < value) lo = mid + 1; else hi = mid;
        }
        return lo;
    }

    #endregion
}
