#pragma warning disable S3776, S107, S1854, S1066, S3358, S1125, S127
using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Analysis;

internal static partial class EmbossSequencePrograms
{
    /// <summary>
    /// <c>plotcon</c> of EMBOSS 6.6.0 (B09 WP67, F308, R14-13): the similarity of an alignment along its length. Every column scores the sum over all pairs of rows of the matrix score of the two residues (the integer matrix of
    /// <c>-scorefile</c>: EBLOSUM62 for protein, EDNAFULL for nucleic sets; a residue that the matrix does not contain, a gap, scores nothing) times the weight of the sequence; the score of the window of <c>-winsize</c> columns that ends at
    /// column k is the sum of the column scores of every row over the window, put at the middle of the window (<c>binmid = k - winsize + 1 + winsize / 2</c>), and divided by <c>numseq (numseq - 1) winsize^2</c> (unsigned 32-bit
    /// arithmetic, single precision as the program computes). The result is the graph of <c>-graph</c> (the devices that write text; the others need plplot and are refused).
    /// </summary>
    internal static EmbossRunResult Plotcon(EmbossAcdArgs a, IReadOnlyDictionary<string, string> files)
    {
        var ctx = EmbossProgramBase.Context("plotcon");
        ctx.InputType ??= "gapany";
        EmbossSeqset? set;
        try
        {
            set = EmbossSeqset.Read("plotcon", a.Raw("sequences") ?? "", files, a.Sequence("sequences"), ctx, "gapany");
        }
        catch (EmbossAcdDieException e)
        {
            return EmbossProgramBase.Died(ctx, e);
        }

        if (set is null)
            return EmbossProgramBase.Died(ctx, "sequences");
        bool protein = !set.IsNuc;
        SubstitutionMatrix matrix;
        int winsizeInt;
        EmbossGraphXy? graph;
        try
        {
            string matrixName = a.Raw("scorefile") ?? (protein ? "EBLOSUM62" : "EDNAFULL");
            if (matrixName.Length == 0)
            {
                ctx.Error("Matrix is required");
                throw a.Die("scorefile");
            }

            var read = ReadMatrix(matrixName, files);
            if (read is null)
            {
                ctx.Error($"Unable to read matrix '{matrixName}'");
                throw a.Die("scorefile");
            }

            matrix = read;
            winsizeInt = a.Int("winsize", 4);
            string baseName = EmbossOutPath.DefaultBaseName("", a.Raw("sequences") ?? "");
            graph = EmbossGraphArgs.Open(a, ctx, "plotcon", "graph", "", 1, baseName, files);
            if (graph is not null)
                graph.Xaxis = a.RawAssoc("graph", "gxtitle") ?? "Relative Residue Position";
        }
        catch (EmbossAcdDieException e)
        {
            return EmbossProgramBase.Died(ctx, e);
        }

        if (graph is null)
            return EmbossProgramBase.Died(ctx, "graph");

        // ajSeqsetFill pads every sequence with gaps to the longest; ajSeqsetFmtUpper
        var seqs = set.Seqs;
        int numseq = seqs.Count;
        var texts = seqs.Select(s => s.SeqText.ToUpperInvariant()).ToList();
        int lenseq = texts.Count == 0 ? 0 : texts.Max(t => t.Length);
        for (int i = 0; i < numseq; i++)
            texts[i] = texts[i].PadRight(lenseq, '-');
        uint winsize = unchecked((uint)winsizeInt);
        int bin2 = (int)((float)winsize / 2.0);
        string alphabet = matrix.Alphabet;
        int Code(char c)
        {
            int k = alphabet.IndexOf(char.ToUpperInvariant(c), StringComparison.Ordinal);
            return k < 0 ? 0 : k + 1;
        }

        float Mij(int m1, int m2) => matrix.TryGetScore(alphabet[m1 - 1], alphabet[m2 - 1], out int sc) ? sc : 0;

        var score = new float[numseq, lenseq];
        var sumscore = new float[lenseq];
        for (int k = 0; k < lenseq; k++)
        {
            for (int i = 0; i < numseq; i++)
            {
                int m1 = Code(texts[i][k]);
                for (int j = i + 1; j < numseq; j++)
                {
                    int m2 = Code(texts[j][k]);
                    if (m1 != 0 && m2 != 0)
                    {
                        float m = Mij(m1, m2);
                        float contrj = (m * seqs[j].Weight) + score[i, k];
                        float contri = (m * seqs[i].Weight) + score[j, k];
                        score[i, k] = contrj;
                        score[j, k] = contri;
                    }
                }
            }

            if ((uint)k >= unchecked(winsize - 1))
            {
                int binup = k + 1;
                int binlo = unchecked((int)((uint)k - winsize + 1));
                int binmid = binlo + bin2;
                for (int i = 0; i < numseq; i++)
                    for (int bin = binlo; bin < binup; bin++)
                        sumscore[binmid] += score[i, bin];
            }
        }

        float denominator = unchecked((uint)numseq * ((uint)numseq - 1) * winsize * winsize);
        for (int bin = 0; bin < lenseq; bin++)
            sumscore[bin] = sumscore[bin] / denominator;

        var data = new EmbossGraphData { Type = "2D Plot", X = new float[lenseq], Y = new float[lenseq] };
        for (int bin = 0; bin < lenseq; bin++)
        {
            data.X[bin] = bin + 1;
            data.Y[bin] = sumscore[bin];
        }

        // ajGraphicsCalcRange
        float ymin = float.MaxValue, ymax = -float.MaxValue;
        for (int i = 0; i < lenseq; i++)
        {
            ymin = ymin < data.Y[i] ? ymin : data.Y[i];
            ymax = ymax > data.Y[i] ? ymax : data.Y[i];
        }

        if (ymin == ymax)
        {
            ymin *= 0.9f;
            ymax *= 1.1f;
        }

        data.SetTruescale(0, lenseq, ymin, ymax);
        graph.Graphs.Add(data);
        graph.Title = "Similarity Plot of Aligned Sequences";
        graph.Yaxis = "Similarity";
        graph.Display(ctx, true);
        string created = graph.Close();                                         // ajGraphicsClose again: a display that could not open its file left its name on the list
        return new EmbossRunResult(graph.Stdout + created, EmbossRunResult.Join(ctx.Messages), 0, new Dictionary<string, string>(graph.Files));
    }
}
