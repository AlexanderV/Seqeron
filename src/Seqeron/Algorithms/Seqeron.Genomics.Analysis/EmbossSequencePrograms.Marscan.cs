#pragma warning disable S3776, S107, S1854, S1066, S3358, S1125, S127, S1905
using System.Globalization;
using Seqeron.Genomics.Alignment;
using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Analysis;

/// <summary><c>marscan</c> of EMBOSS 6.6.0 (<c>emboss/marscan.c</c>): the MAR/SAR recognition signatures of a nucleotide sequence (B09 WP63, F282).</summary>
internal static partial class EmbossSequencePrograms
{
    private const int MarscanMaxDist = 200;

    // marscan_getpos: takes the next position of a list while it is not past the other list's last position and stores or outputs the matches within MAXDIST of each other
    private sealed class MarscanState
    {
        public bool StoredMatch;
        public int StoredDist, Stored8Pos, Stored16Pos, StoredLastPos;
        public int EndOfLastOutput = -(MarscanMaxDist + 1);
        public readonly List<EmbossFeature> Features = [];
    }

    private static bool MarscanGetPos(Queue<int> l, ref int thisprev, int otherprev, MarscanState st, bool thisIs8)
    {
        while (thisprev <= otherprev)
        {
            if (l.Count == 0)
                return false;
            thisprev = l.Dequeue();
            int s8 = thisIs8 ? thisprev : otherprev;
            int s16 = thisIs8 ? otherprev : thisprev;
            int e8 = s8 + 7, e16 = s16 + 15;
            int lastpos = e8 > e16 ? e8 : e16;
            int dist = e8 < s16 ? s16 - e8 : e16 < s8 ? s8 - e16 : 0;
            int firstpos = s8 > s16 ? s16 : s8;
            if (otherprev == -1)
            {
                dist = MarscanMaxDist + 1;
                firstpos = MarscanMaxDist + 1;
            }

            if (dist < 0)
                throw new InvalidOperationException("marscan: a negative distance (ajFatal).");
            if (dist <= MarscanMaxDist)
            {
                if (st.StoredMatch)
                {
                    if (dist < st.StoredDist)
                    {
                        if (firstpos - st.EndOfLastOutput > MarscanMaxDist)
                            MarscanStore(st, dist, s8, s16, lastpos);
                    }
                    else
                    {
                        MarscanOutput(st);
                        st.StoredMatch = false;
                        st.EndOfLastOutput = st.StoredLastPos;
                    }
                }
                else if (firstpos - st.EndOfLastOutput > MarscanMaxDist)
                {
                    MarscanStore(st, dist, s8, s16, lastpos);
                }
            }
            else if (st.StoredMatch)
            {
                MarscanOutput(st);
                st.StoredMatch = false;
                st.EndOfLastOutput = st.StoredLastPos;
            }
        }

        return true;
    }

    private static void MarscanStore(MarscanState st, int dist, int s8, int s16, int lastpos)
    {
        st.StoredMatch = true;
        st.StoredDist = dist;
        st.Stored8Pos = s8;
        st.Stored16Pos = s16;
        st.StoredLastPos = lastpos;
    }

    // marscan_output_stored_match
    private static void MarscanOutput(MarscanState st)
    {
        if (!st.StoredMatch)
            return;
        int s8 = st.Stored8Pos, s16 = st.Stored16Pos;
        int e8 = s8 + 7, e16 = s16 + 15;
        int start = s8 < s16 ? s8 : s16;
        int end = e8 > e16 ? e8 : e16;
        string I(int v) => v.ToString(CultureInfo.InvariantCulture);
        st.Features.Add(new EmbossFeature("SO:0005836", "marscan", start, end, 0.0, '+', 0,
            [new("note", "*type MAR/SAR recognition site (MRS)"), new("note", "*start8bp " + I(s8)), new("note", "*end8bp " + I(e8)), new("note", "*start16bp " + I(s16)), new("note", "*end16bp " + I(e16))]));
    }

    // marscan_stepdown
    private static void MarscanStepDown(Queue<int> l16, Queue<int> l8, MarscanState st)
    {
        int prev16 = -1, prev8 = -1;
        bool notend16 = true, notend8 = true;
        while (notend16 || notend8)
        {
            notend16 = MarscanGetPos(l16, ref prev16, prev8, st, false);
            if (prev16 >= prev8 && !notend8)
                notend16 = false;
            notend8 = MarscanGetPos(l8, ref prev8, prev16, st, true);
            if (prev8 >= prev16 && !notend16)
                notend8 = false;
        }

        if (st.Stored8Pos - st.EndOfLastOutput > MarscanMaxDist || st.Stored16Pos - st.EndOfLastOutput > MarscanMaxDist)
            MarscanOutput(st);
    }

    /// <summary>
    /// <c>marscan</c>: the MRS signatures: a 16 bp site (awwrtaannwwgnnnc, one mismatch allowed, or its reverse gnnncwwnnttaywwt) and an 8 bp site (aataayaa or ttrttatt, exact) within 200 bp of each other, one signature for each cluster
    /// (the pair with the least distance), reported as GFF (default) or any report format.
    /// </summary>
    internal static EmbossRunResult Marscan(EmbossAcdArgs a, IReadOnlyDictionary<string, string> files)
    {
        var ctx = EmbossProgramBase.Context("marscan");
        var all = EmbossProgramBase.Open("marscan", a, files, ctx, "dna");
        if (all is null)
            return EmbossProgramBase.Died(ctx);
        EmbossReportArgs.Opened? report;
        try
        {
            report = EmbossReportArgs.Open(a, ctx, "outfile", "marscan", "gff", [new("str", "type", "type"), new("int", "start8bp", "start8bp"), new("int", "end8bp", "end8bp"), new("int", "start16bp", "start16bp"), new("int", "end16bp", "end16bp")], 3, "", a.Raw("outfile"),
                EmbossOutPath.DefaultBaseName(all.Peek().Name, a.Raw("sequence") ?? ""), EmbossRunEnvironment.Current.RunDate ?? DateTime.Now, files);
        }
        catch (EmbossAcdDieException e)
        {
            return EmbossProgramBase.Died(ctx, e);
        }

        if (report is null)
            return EmbossProgramBase.Died(ctx, "outfile");
        var tables = new List<EmbossReport>();
        long count = 0, total = 0;
        while (EmbossSeqretPrograms.TryNext(all, ctx, out var seq))
        {
            count++;
            total += EmbossSeqRange.GetEnd(seq) - EmbossSeqRange.GetBegin(seq) + 1;
            int begin = all.GetSeqBegin(seq), end = all.GetSeqEnd(seq);
            string text = EmbossProgramBase.SubS(seq.SeqText, begin - 1, end - 1).ToUpperInvariant();
            Queue<int> Search(string pattern, int mismatches)
            {
                var hits = GenomicAnalyzer.FindPatternHits(text, pattern, mismatches, false, default, true);
                var starts = hits.Select(h => h.Position + begin).OrderBy(x => x).ToList();
                return new Queue<int>(starts);
            }

            var l16 = Search("AWWRTAANNWWGNNNC", 1);
            var l16rev = Search("GNNNCWWNNTTAYWWT", 1);
            var l8 = Search("AATAAYAA", 0);
            var l8rev = Search("TTRTTATT", 0);
            var m8 = new Queue<int>(l8.Concat(l8rev).OrderBy(x => x));
            var m16 = new Queue<int>(l16.Concat(l16rev).OrderBy(x => x));
            var st = new MarscanState();
            MarscanStepDown(m16, m8, st);
            var sequence = new EmbossReportSequence(seq.Name, EmbossFeatureReports.MakeUsa(seq), seq.Acc, seq.Desc, seq.Seq.ToString(), EmbossSeqRange.GetBegin(seq), EmbossSeqRange.GetEnd(seq), seq.Offset, seq.Reversed, true);
            long length = st.Features.Count == 0 ? 0 : st.Features.Max(f => (long)Math.Max(f.Start, f.End));
            tables.Add(new EmbossReport(sequence, st.Features, false, "N", TableLength: length, Tags: report.Options.Tags));
        }

        string textOut = EmbossReportArgs.Write(report, tables, ctx, count, total);
        return EmbossReportArgs.Result(ctx, report.File, textOut);
    }
}
