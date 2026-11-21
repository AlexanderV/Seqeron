// The control flow of the ported C functions is kept literally (loop variables moved inside the body, empty search loops, the C compound conditions).
#pragma warning disable S127, S108, S1994, S3358, S1066, S1854

namespace Seqeron.Genomics.Alignment;

public static partial class Minimap2Chaining
{
    // the parts of mm_idx_t the base-level stage reads: the 4-bit sequences (seq_nt4_table: A C G T = 0..3, everything else 4), the names, the ALT flags (--alt) and the junction intervals (--junc-bed)
    internal sealed class MmIndex
    {
        internal required byte[][] Seq { get; init; }
        internal required int[] Length { get; init; }
        internal string?[] Names { get; init; } = [];
        internal int K { get; init; }
        internal bool Hpc { get; init; }
        internal bool[] IsAlt { get; set; } = [];
        internal int NAlt { get; set; }
        internal List<(int St, int En, int Strand)>[]? Junctions { get; set; }

        internal int Count => Seq.Length;

        internal static MmIndex Create(IReadOnlyList<string> targets, IReadOnlyList<string?>? names, int k, bool hpc)
        {
            var seq = new byte[targets.Count][];
            for (int i = 0; i < seq.Length; i++)
            {
                var s = targets[i];
                var b = new byte[s.Length];
                for (int j = 0; j < b.Length; j++)
                    b[j] = (byte)Nt4(s[j]);
                seq[i] = b;
            }

            return new MmIndex
            {
                Seq = seq, Length = [.. targets.Select(t => t.Length)], Names = names is null ? new string?[targets.Count] : [.. names], K = k, Hpc = hpc, IsAlt = new bool[targets.Count],
            };
        }

        // mm_idx_getseq: the target bases [st, en) (en clipped to the sequence length), or null when st is outside the sequence
        internal byte[]? GetSeq(int rid, int st, int en)
        {
            if (rid < 0 || rid >= Count || st < 0 || st >= Length[rid])
                return null;
            if (en > Length[rid])
                en = Length[rid];
            var r = new byte[Math.Max(0, en - st)];
            Array.Copy(Seq[rid], st, r, 0, r.Length);
            return r;
        }

        // mm_idx_getseq_rev: the reverse complement of the bases [st, en)
        internal byte[]? GetSeqRev(int rid, int st, int en)
        {
            if (rid < 0 || rid >= Count || st < 0 || st >= Length[rid])
                return null;
            if (en > Length[rid])
                en = Length[rid];
            int st1 = Length[rid] - en, en1 = Length[rid] - st;
            var r = new byte[Math.Max(0, en - st)];
            for (int i = st1; i < en1; i++)
            {
                byte c = Seq[rid][i];
                r[en1 - i - 1] = c < 4 ? (byte)(3 - c) : c;
            }

            return r;
        }

        // mm_idx_getseq2
        internal byte[] GetSeq2(bool rev, int rid, int st, int en) => (rev ? GetSeqRev(rid, st, en) : GetSeq(rid, st, en)) ?? [];

        // mm_idx_bed_junc: the junction flags of the positions [st, en): bit 0 / 1 = donor / acceptor of a + junction, bit 3 / 2 = donor / acceptor of a - junction
        internal byte[] BedJunction(int ctg, int st, int en)
        {
            var s = new byte[Math.Max(0, en - st)];
            if (Junctions is null || ctg < 0 || ctg >= Count)
                return s;
            foreach (var (jst, jen, strand) in Junctions[ctg])
            {
                if (jst < st || !(st <= jst && en >= jen) || strand == 0)
                    continue;
                if (strand > 0)
                {
                    s[jst - st] |= 1;
                    s[jen - 1 - st] |= 2;
                }
                else
                {
                    s[jst - st] |= 8;
                    s[jen - 1 - st] |= 4;
                }
            }

            return s;
        }

        // mm_idx_alt_read: the first word of every line is a target name; matching targets become ALT contigs
        internal void ReadAlt(string? text)
        {
            if (text is null)
                return;
            int n = 0;
            foreach (string line in text.Split('\n'))
            {
                string l = line.TrimEnd('\r');
                int e = 0;
                while (e < l.Length && !char.IsWhiteSpace(l[e]))
                    e++;
                string name = l[..e];
                int id = Array.IndexOf(Names, name);
                if (id >= 0)
                {
                    IsAlt[id] = true;
                    ++n;
                }
            }

            NAlt = n;
        }

        // mm_idx_bed_read with read_junc = 1 (the BED12 blocks give the introns; a BED6 line is one interval with its strand)
        internal void ReadJunctionBed(string text)
        {
            Junctions = new List<(int, int, int)>[Count];
            for (int i = 0; i < Count; i++)
                Junctions[i] = [];
            foreach (string line in text.Split('\n'))
            {
                string l = line.TrimEnd('\r');
                var f = l.Split('\t');
                int id = f.Length > 0 ? Array.IndexOf(Names, f[0]) : -1;
                if (id < 0 || f[0].Length == 0)
                    continue;
                long stL = -1, enL = -1;
                int strand = 0;
                int nBlk = 0;
                string? bl = null, bs = null;
                bool ok = true;
                bool stop = false;
                for (int i = 1; i < f.Length && ok && !stop; i++)
                {
                    switch (i)
                    {
                        case 1:
                            stL = CAtol(f[i]);
                            ok = stL >= 0;
                            break;
                        case 2:
                            enL = CAtol(f[i]);
                            ok = enL >= 0;
                            break;
                        case 5:
                            strand = f[i].Length > 0 && f[i][0] == '+' ? 1 : f[i].Length > 0 && f[i][0] == '-' ? -1 : 0;
                            break;
                        case 9:
                            if (f[i].Length == 0 || !char.IsAsciiDigit(f[i][0]))
                                stop = true;
                            else
                                nBlk = (int)CAtol(f[i]);
                            break;
                        case 10:
                            bl = f[i];
                            break;
                        case 11:
                            bs = f[i];
                            break;
                    }

                    if (i == 11)
                        stop = true;
                }

                if (!ok || stL < 0 || enL < 0 || stL >= enL)
                    continue;
                int st = (int)stL, en = (int)enL;
                if (bs is not null && bl is not null)
                {
                    var starts = bs.Split(',', StringSplitOptions.None);
                    var sizes = bl.Split(',', StringSplitOptions.None);
                    int e0 = st + (int)CAtol(starts[0]) + (int)CAtol(sizes[0]);
                    for (int b = 1; b < nBlk && b < starts.Length && b < sizes.Length; b++)
                    {
                        int jst = e0, jen = st + (int)CAtol(starts[b]);
                        e0 = st + (int)CAtol(starts[b]) + (int)CAtol(sizes[b]);
                        if (jen > jst)
                            Junctions[id].Add((jst, jen, strand));
                    }
                }
                else
                {
                    Junctions[id].Add((st, en, strand));
                }
            }

            for (int i = 0; i < Count; i++)
                Junctions[i] = [.. Junctions[i].OrderBy(x => x.St)];
        }

        // atol: the leading decimal integer (optional sign), 0 when there is none
        private static long CAtol(string s)
        {
            int i = 0;
            while (i < s.Length && char.IsWhiteSpace(s[i]))
                i++;
            bool neg = false;
            if (i < s.Length && (s[i] == '+' || s[i] == '-'))
                neg = s[i++] == '-';
            long v = 0;
            while (i < s.Length && char.IsAsciiDigit(s[i]))
                v = (v * 10) + (s[i++] - '0');
            return neg ? -v : v;
        }
    }
}
