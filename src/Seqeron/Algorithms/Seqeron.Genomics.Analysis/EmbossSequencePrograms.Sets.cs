#pragma warning disable S3776, S107, S1854, S1066, S3358, S1125
using System.Globalization;
using System.Text;
using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Analysis;

/// <summary>The programs of EMBOSS 6.6.0 that read sets of sequences (<c>seqset</c> / <c>seqsetall</c>): <c>listor</c>, <c>nthseqset</c> (B09 WP62, F281).</summary>
internal static partial class EmbossSequencePrograms
{
    // ---------------------------------------------------------------------------------------------------------------- nthseqset

    /// <summary><c>nthseqset</c>: the <c>-number</c>-th set of the sets of the input (a file of several alignments) written as a set.</summary>
    internal static EmbossRunResult Nthseqset(EmbossAcdArgs a, IReadOnlyDictionary<string, string> files)
    {
        var ctx = EmbossProgramBase.Context("nthseqset");
        bool features;
        try
        {
            features = a.Bool("feature", false);
        }
        catch (EmbossAcdDieException e)
        {
            return EmbossProgramBase.Died(ctx, e);
        }

        ctx.InputType ??= "gapany";
        var sets = EmbossSeqset.ReadAll("nthseqset", a.Raw("sequence") ?? "", files, a.Sequence("sequence"), ctx, "gapany", features);
        if (sets is null)
            return EmbossProgramBase.Died(ctx);
        int n;
        EmbossOutseqOptions outOptions;
        try
        {
            n = a.IntLimited("number", 1, 1, int.MaxValue, ctx);
            outOptions = a.Seqout("outseq");
        }
        catch (EmbossAcdDieException e)
        {
            return EmbossProgramBase.Died(ctx, e);
        }

        var stream = new EmbossSeqoutStream("nthseqset", ctx, outOptions, null);
        if (stream.OpenFailed)
            return EmbossRunResult.From(stream.DiedOutput());
        int count = sets.Count;
        if (count < n)
        {
            stream.Close();                                                    // ajSeqoutClose: the warning for an empty output, then the fatal error
            return FatalAt(ctx, "nthseqset.c", 55, $"No such sequence set - only {count} sequence sets were input.");
        }

        stream.WriteSet(sets[n - 1].Seqs);
        return Finish(stream);
    }

    // ---------------------------------------------------------------------------------------------------------------- listor

    private const string ListorValues = "OR:OR - merger of both sets, AND:AND - only those in both sets, XOR:XOR - only those not in both sets, NOT:NOT - those of the first set that are not in the second";

    /// <summary><c>listor</c>: a list file of the USAs of the logical OR / AND / XOR / NOT of two sets of sequences (sequences are equal when their residues are equal, ignoring the case).</summary>
    internal static EmbossRunResult Listor(EmbossAcdArgs a, IReadOnlyDictionary<string, string> files)
    {
        var ctx = EmbossProgramBase.Context("listor");
        ctx.InputType ??= "gapany";
        var seq1 = EmbossSeqset.Read("listor", a.Raw("firstsequences") ?? "", files, a.Sequence("firstsequences"), ctx, "gapany");
        if (seq1 is null)
            return EmbossProgramBase.Died(ctx, "firstsequences");
        var seq2 = EmbossSeqset.Read("listor", a.Raw("secondsequences") ?? "", files, a.Sequence("secondsequences"), ctx, "gapany");
        if (seq2 is null)
            return EmbossProgramBase.Died(ctx, "secondsequences");
        string op;
        string outfile;
        try
        {
            op = a.List("operator", "OR", ListorValues, 1, 1, ",")[0];
            outfile = a.Raw("outfile") ?? DefaultOutfile(a, files, "listor", "firstsequences");
            if (outfile.Length == 0)
                throw a.Die("outfile", "Error: Output file is required");
        }
        catch (EmbossAcdDieException e)
        {
            return EmbossProgramBase.Died(ctx, e);
        }

        int n1 = seq1.Seqs.Count, n2 = seq2.Seqs.Count;
        var lengths1 = new int[n1];
        var hits1 = new int[n1];
        var order1 = new uint[n1];
        for (int i = 0; i < n1; i++)
        {
            lengths1[i] = seq1.Seqs[i].Seq.Length;
            order1[i] = (uint)i;
            hits1[i] = -1;
        }

        EmbossRandom.SortIntIncI(lengths1, order1, n1);
        var lengths2 = new int[n2];
        var hits2 = new int[n2];
        var order2 = new uint[n2];
        for (int i = 0; i < n2; i++)
        {
            lengths2[i] = seq2.Seqs[i].Seq.Length;
            order2[i] = (uint)i;
            hits2[i] = -1;
        }

        EmbossRandom.SortIntIncI(lengths2, order2, n2);
        int curr1 = 0, curr2 = 0;
        while (curr1 < n1 && curr2 < n2)
        {
            if (lengths1[order1[curr1]] < lengths2[order2[curr2]])
                curr1++;
            else if (lengths1[order1[curr1]] > lengths2[order2[curr2]])
                curr2++;
            else
            {
                int tmp1, tmp2 = 0;
                for (tmp1 = curr1; tmp1 < n1 && lengths1[order1[tmp1]] == lengths2[order2[curr2]]; tmp1++)
                {
                    for (tmp2 = curr2; tmp2 < n2 && lengths2[order2[tmp2]] == lengths2[order2[curr2]]; tmp2++)
                    {
                        if (string.Equals(seq1.Seqs[(int)order1[tmp1]].SeqText, seq2.Seqs[(int)order2[tmp2]].SeqText, StringComparison.OrdinalIgnoreCase))
                        {
                            hits1[order1[tmp1]] = (int)order2[tmp2];
                            hits2[order2[tmp2]] = (int)order1[tmp1];
                        }
                    }
                }

                curr1 = tmp1;
                curr2 = tmp2;
            }
        }

        var sb = new StringBuilder();
        for (int i = 0; i < n1; i++)
        {
            bool write = op switch
            {
                "OR" => true,
                "AND" => hits1[i] != -1,
                _ => hits1[i] == -1,                  // XOR, NOT
            };
            if (write)
                sb.Append(UsaOf(seq1.Seqs[i])).Append('\n');
        }

        if (op is "OR" or "XOR")
        {
            for (int i = 0; i < n2; i++)
            {
                if (hits2[i] == -1)
                    sb.Append(UsaOf(seq2.Seqs[i])).Append('\n');
            }
        }

        return OutfileResult(ctx, outfile, sb.ToString());
    }
}
