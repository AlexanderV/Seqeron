#pragma warning disable S127, S3358, S1066, S3776, S107, S1854, IDE0059, S2325
using System.Globalization;
using System.Text;

namespace Seqeron.Genomics.Alignment;

public static partial class SequenceAligner
{
    #region MAFFT 7.505: adding sequences to an alignment (addfunctions.c, disttbfast.c with -I)

    /// <summary>The state of <c>disttbfast -I n</c> (<c>--add</c>, <c>--addprofile</c>): which step of the guide tree does what (<c>mergeoralign</c>) and which sequences are already aligned.</summary>
    internal sealed class MafftAddState
    {
        internal MafftAddState(char[] mergeOrAlign, bool[] alreadyAligned)
        {
            MergeOrAlign = mergeOrAlign;
            AlreadyAligned = alreadyAligned;
        }

        /// <summary><c>mergeoralign[l]</c>: 'n' nothing to do, '1' / '2' the first / second group is an alignment that must keep its columns, 'w' / 'a' ordinary profile alignment.</summary>
        internal char[] MergeOrAlign { get; }

        /// <summary><c>alreadyaligned</c>.</summary>
        internal bool[] AlreadyAligned { get; }
    }

    internal sealed partial class MafftRun
    {
        /// <summary><c>strncpy0</c> + C string semantics: at most <paramref name="n"/> characters of <paramref name="s"/> from <paramref name="start"/> (a NUL, the end of the string, ends the copy).</summary>
        private static string CopyN(string s, int start, int n)
        {
            if (start >= s.Length || n <= 0)
                return "";
            return s.Substring(start, Math.Min(n, s.Length - start));
        }

        /// <summary><c>findcommongaps</c>: for every position of the alignment without its common gaps the number of common-gap columns in front of it (index = number of kept columns before; length + 1 entries).</summary>
        private static int[] FindCommonGaps(int n, IReadOnlyList<string> seq)
        {
            int len = seq[0].Length;
            var gapmap = new int[len + 1];
            int pos = 0;
            for (int i = 0; i < len; i++)
            {
                int j;
                for (j = 0; j < n; j++)
                    if (seq[j][i] != '-')
                        break;
                if (j == n)
                    gapmap[pos]++;
                else
                    pos++;
            }
            return gapmap;
        }

        /// <summary><c>commongappick( n, seq )</c>: the columns in which the first <paramref name="n"/> rows have gaps only are removed (from all <paramref name="n"/> rows).</summary>
        private static void CommonGapPickN(int n, string[] seq)
        {
            int len = seq[0].Length;
            var keep = new List<int>(len);
            for (int i = 0; i < len; i++)
            {
                int j;
                for (j = 0; j < n; j++)
                    if (seq[j][i] != '-')
                        break;
                if (j != n)
                    keep.Add(i);
            }
            var sb = new StringBuilder(keep.Count);
            for (int j = 0; j < n; j++)
            {
                sb.Clear();
                foreach (int i in keep)
                    sb.Append(seq[j][i]);
                seq[j] = sb.ToString();
            }
        }

        /// <summary><c>adjustgapmap</c>: the gap map for the alignment <paramref name="seq"/> (a '=' marks a column that the new alignment inserted: nothing was removed in front of it).</summary>
        private static int[] AdjustGapMap(int newlen, int[] gapmap, string seq)
        {
            var tmpmap = new int[newlen + 2];
            int j = 0, pos = 0;
            foreach (char c in seq)
            {
                if (c == '=')
                    tmpmap[j++] = 0;
                else
                    tmpmap[j++] = pos < gapmap.Length ? gapmap[pos++] : 0;
            }
            tmpmap[j++] = pos < gapmap.Length ? gapmap[pos] : 0;
            var res = new int[newlen + 1];
            Array.Copy(tmpmap, res, newlen + 1);
            return res;
        }

        /// <summary><c>restorecommongaps</c>: puts the removed common-gap columns back into the rows of the two groups; returns the gap map in the coordinates of the restored rows (every gap column is followed by a zero entry per column).</summary>
        private int[] RestoreCommonGaps(int njob, int n0, string[] seq, int[] ex1, int[] ex2, int[] gapmap, char gapchar)
        {
            if (n0 == 0)
                return gapmap;
            var mem = new List<int>(ex1);
            mem.AddRange(ex2);
            int rep = mem[0];
            int len = seq[rep].Length;
            int len1 = len + 1;
            var sb = new StringBuilder();
            foreach (int i in mem)
            {
                sb.Clear();
                for (int j = 0; j < len1; j++)
                {
                    int klim = j < gapmap.Length ? gapmap[j] : 0;
                    sb.Append(gapchar, klim);
                    if (j < len)
                        sb.Append(seq[i][j]);
                }
                seq[i] = sb.ToString();
            }
            var tmpgapmap = new List<int>();
            for (int j = 0; j < len1; j++)
            {
                int g = j < gapmap.Length ? gapmap[j] : 0;
                tmpgapmap.Add(g);
                for (int k = 0; k < g; k++)
                    tmpgapmap.Add(0);
            }
            return tmpgapmap.ToArray();
        }

        /// <summary><c>findnewgaps</c>: the number of inserted ('=') columns in front of every old column of the row <paramref name="rep"/>.</summary>
        private static int[] FindNewGaps(string[] seq, int rep)
        {
            int len = seq[0].Length;
            var gaplen = new int[len + 1];
            int pos = 0;
            for (int i = 0; i < len; i++)
            {
                if (seq[rep][i] == '=')
                    gaplen[pos]++;
                else
                    pos++;
            }
            return gaplen;
        }

        /// <summary><c>eq2dashmatometehayaku</c>: the columns that are '=' in the first row become '-' in all <paramref name="n"/> rows.</summary>
        private static void Eq2DashMatome(string[] s, int n)
        {
            var cols = new List<int>();
            for (int j = 0; j < s[0].Length; j++)
                if (s[0][j] == '=')
                    cols.Add(j);
            if (cols.Count == 0)
                return;
            for (int i = 0; i < n; i++)
            {
                var chars = s[i].ToCharArray();
                foreach (int j in cols)
                    chars[j] = '-';
                s[i] = new string(chars);
            }
        }

        /// <summary><c>recordoriginalgaps</c>: '-' for a column in which all of the first <paramref name="n"/> rows have a gap, 'o' for the others.</summary>
        private static string RecordOriginalGaps(int n, IReadOnlyList<string> s)
        {
            int len = s[0].Length;
            var sb = new StringBuilder(len);
            for (int i = 0; i < len; i++)
            {
                int j;
                for (j = 0; j < n; j++)
                    if (s[j][i] != '-')
                        break;
                sb.Append(j == n ? '-' : 'o');
            }
            return sb.ToString();
        }

        /// <summary><c>restoreoriginalgaps</c>: the columns that only gaps of the original alignment had are put back into the first <paramref name="n"/> rows.</summary>
        private static void RestoreOriginalGaps(int n, string[] seq, string originalgaps)
        {
            for (int i = 0; i < n; i++)
            {
                var sb = new StringBuilder(originalgaps.Length);
                int p = 0;
                foreach (char og in originalgaps)
                    sb.Append(og == '-' ? '-' : seq[i][p++]);
                seq[i] = sb.ToString();
            }
        }

        /// <summary>The runs of residues that <c>--keeplength</c> deleted from one added sequence: position in the sequence (0-based, residues only) and length.</summary>
        internal readonly record struct GapPos(int Pos, int Len);

        /// <summary><c>deletenewinsertions_whole</c>: the columns in which only added sequences have residues are removed (<c>--keeplength</c>); returns the deleted runs of every added sequence.</summary>
        private static List<GapPos>[] DeleteNewInsertions(int on, int an, string[] oseq, string[] aseq)
        {
            int len = oseq[0].Length;
            var eq = new bool[len + 1]; // eqseq[len] = 0
            for (int i = 0; i < len; i++)
            {
                int j;
                for (j = 0; j < on; j++)
                    if (oseq[j][i] != '-')
                        break;
                eq[i] = j == on;
            }
            var deleted = new List<GapPos>[an];
            for (int j = 0; j < an; j++)
            {
                deleted[j] = new List<GapPos>();
                int q = 0;
                string a = aseq[j];
                for (int i = 0; i < len; i++)
                {
                    char tmpc = a[i];
                    int qstart = q;
                    while (tmpc != '-' && eq[i])
                    {
                        i++;
                        q++;
                        tmpc = i < a.Length ? a[i] : '\0';
                    }
                    if (qstart != q)
                        deleted[j].Add(new GapPos(qstart, q - qstart));
                    if (tmpc != '-')
                        q++;
                }
            }
            for (int j = 0; j < on; j++)
                oseq[j] = new string(oseq[j].Where((_, i) => !eq[i]).ToArray());
            for (int j = 0; j < an; j++)
                aseq[j] = new string(aseq[j].Where((_, i) => i >= eq.Length || !eq[i]).ToArray());
            return deleted;
        }

        /// <summary>The name of an added sequence as the map files print it (<c>strstr( nameptr, "_oe_" )</c> and the outnumber prefix removed).</summary>
        private string MapName(string programName)
        {
            string nameptr = programName[1..];
            if (G.OutNumber != 0)
            {
                int c = nameptr.IndexOf("_numo_e", StringComparison.Ordinal);
                nameptr = c >= 0 ? nameptr[(c + 8)..] : "";
            }
            int oe = nameptr.IndexOf("_oe_", StringComparison.Ordinal);
            return oe >= 0 ? nameptr[(oe + 4)..] : nameptr;
        }

        /// <summary><c>reconstructdeletemap</c> (<c>--mapout</c>): for every added sequence the letters, their positions in the original sequence and in the reference alignment ('-' for a deleted letter).</summary>
        private string ReconstructDeleteMap(int nadd, string[] addbk, List<GapPos>[] deletelist, string[] realn, IReadOnlyList<string> names)
        {
            var sb = new StringBuilder();
            for (int i = 0; i < nadd; i++)
            {
                int len = addbk[i].Length;
                var gapped = new bool[len + 1];
                sb.Append('>').Append(MapName(names[i])).Append('\n');
                sb.Append("# letter, position in the original sequence, position in the reference alignment\n");
                foreach (var d in deletelist[i])
                    for (int g = 0; g < d.Len; g++)
                        gapped[d.Pos + g] = true;
                int p = 0;
                for (int j = 0; j < len; j++)
                {
                    while (p < realn[i].Length && realn[i][p] == '-')
                        p++;
                    if (gapped[j])
                    {
                        sb.Append(addbk[i][j]).Append(", ").Append(j + 1).Append(", -\n");
                    }
                    else
                    {
                        sb.Append(addbk[i][j]).Append(", ").Append(j + 1).Append(", ").Append(p + 1).Append('\n');
                        p++;
                    }
                }
            }
            return sb.ToString();
        }

        /// <summary><c>reconstructdeletemap_compact</c> (<c>--compactmapout</c>): the insertions as "9A - 11C &gt; 8v9"; the title of a sequence with insertions is changed to "Nins:...|title" (the new program names are in <paramref name="newNames"/>).</summary>
        private string ReconstructDeleteMapCompact(int nadd, string[] addbk, List<GapPos>[] deletelist, string[] realn, string[] names, string[] newNames)
        {
            var sb = new StringBuilder();
            sb.Append("# Insertion in added sequence > Position in reference\n");
            for (int i = 0; i < nadd; i++)
            {
                int len = addbk[i].Length;
                var gapped = new bool[len + 1];
                string programName = names[i];
                string nameptr = programName[1..];
                int prefix = 1;
                if (G.OutNumber != 0)
                {
                    int c = nameptr.IndexOf("_numo_e", StringComparison.Ordinal);
                    nameptr = c >= 0 ? nameptr[(c + 8)..] : "";
                    prefix = programName.Length - nameptr.Length;
                }
                int oe = nameptr.IndexOf("_oe_", StringComparison.Ordinal);
                if (oe >= 0)
                {
                    nameptr = nameptr[(oe + 4)..];
                    prefix = programName.Length - nameptr.Length;
                }
                bool any = false;
                foreach (var d in deletelist[i])
                {
                    for (int g = 0; g < d.Len; g++)
                        gapped[d.Pos + g] = true;
                    any = true;
                }
                if (!any)
                    continue;
                sb.Append('>').Append(nameptr).Append('\n');
                int status = -1;
                var insstr = new StringBuilder();
                int nins = 0;
                int p = 0;
                int j;
                for (j = 0; j < len; j++)
                {
                    while (p < realn[i].Length && realn[i][p] == '-')
                        p++;
                    if (gapped[j])
                    {
                        if (status != 1)
                        {
                            status = 1;
                            sb.Append(j + 1).Append(addbk[i][j]).Append(" - ");
                            if (nins == 0)
                            {
                                insstr.Clear();
                                insstr.Append(j + 1).Append(addbk[i][j]).Append('-');
                            }
                            else if (nins == 1)
                            {
                                insstr.Append("etc,");
                            }
                            nins++;
                        }
                    }
                    else
                    {
                        if (status == 1)
                        {
                            sb.Append(j).Append(addbk[i][j - 1]).Append(" > ").Append(p).Append('v').Append(p + 1).Append('\n');
                            if (nins == 1)
                                insstr.Append(j).Append(addbk[i][j - 1]).Append(',');
                        }
                        status = 0;
                        p++;
                    }
                }
                if (status == 1)
                {
                    sb.Append(j).Append(addbk[i][j - 1]).Append(" > ").Append(p).Append('v').Append(p + 1).Append('\n');
                    if (nins == 1)
                        insstr.Append(j).Append(addbk[i][j - 1]).Append(',');
                }
                if (insstr.Length > 0)
                    insstr.Length--;
                string newname = programName[..prefix] + nins.ToString(CultureInfo.InvariantCulture) + "ins:" + insstr + "|" + nameptr;
                newNames[i] = newname.Length > 256 ? newname[..256] : newname;
            }
            return sb.ToString();
        }

        /// <summary><c>profilealignment</c> (<c>addfunctions.c</c>): the sequences of group 0 (aligned before, not in the two groups) and group 2 (a group of the step, with the gaps of the new columns) are aligned as profiles; group 1 follows them with '=' in the columns that only group 2 has.</summary>
        private void ProfileAlignment(int n0, int n1, int n2, string[] aln0, string[] aln1, string[] aln2, int alloclen)
        {
            if (aln0[0].Length == 1 && aln2[0].Length == 1)
                return;
            CommonGapPickN(n0, aln0);
            CommonGapPickN(n2, aln2);
            if (aln2[0].Length == 0)
            {
                int newlen = aln0[0].Length;
                for (int i = 0; i < n2; i++)
                    aln2[i] = new string('-', newlen);
                return;
            }
            var effarr0 = new double[n0];
            var effarr2 = new double[n2];
            var allgap0 = new bool[n0];
            var allgap2 = new bool[n2];
            int alcount0 = 0, alcount2 = 0;
            for (int i = 0; i < n0; i++)
            {
                allgap0[i] = aln0[i].All(c => c == '-');
                if (!allgap0[i])
                    alcount0++;
            }
            for (int i = 0; i < n2; i++)
            {
                allgap2[i] = aln2[i].All(c => c == '-');
                if (!allgap2[i])
                    alcount2++;
            }
            for (int i = 0; i < n0; i++)
                effarr0[i] = allgap0[i] ? 0.0 : 1.0 / alcount0;
            for (int i = 0; i < n2; i++)
                effarr2[i] = allgap2[i] ? 0.0 : 1.0 / alcount2;
            char saved = G.NewGapChar;
            G.NewGapChar = '-';
            if (G.Alg == 'M')
                MsAlign(G.NDisConsweightMulti, aln0, aln2, effarr0, effarr2, n0, n2, alloclen, null, null, null, null, 1, 1, null, -1, -1, 0, 0.0, 0.0);
            else
                AAlign(G.NDisConsweightMulti, G.Penalty, G.PenaltyEx, aln0, aln2, effarr0, effarr2, n0, n2, alloclen, null, null, null, null, 1, 1, null, -1, -1, 0, 0.0, 0.0);
            G.NewGapChar = saved;
            int nl = aln0[0].Length;
            for (int i = 0; i < n1; i++)
                aln1[i] = new string('-', nl);
            var chars1 = aln1.Select(x => x.ToCharArray()).ToArray();
            for (int j = 0; j < nl; j++)
            {
                int i;
                for (i = 0; i < n0; i++)
                    if (aln0[i][j] != '-')
                        break;
                if (i == n0)
                {
                    for (i = 0; i < n1; i++)
                        if (chars1[i][j] != '-')
                            break;
                }
                else
                {
                    i = -1;
                }
                if (i == n1)
                    for (i = 0; i < n1; i++)
                        chars1[i][j] = '=';
            }
            for (int i = 0; i < n1; i++)
                aln1[i] = new string(chars1[i]);
        }

        /// <summary>
        /// <c>insertnewgaps</c> (<c>addfunctions.c</c>): the columns that the alignment of the step inserted ('=' in the group that was aligned) are put into the sequences that are already aligned and do not take part in the step
        /// (<paramref name="seq"/> are the rows of all sequences; the two groups are <paramref name="ex1"/> and <paramref name="ex2"/>).
        /// </summary>
        private void InsertNewGaps(int njob, bool[] alreadyaligned, string[] seq, int[] ex1, int[] ex2, int[] gaplen, int[] gapmap, int alloclen, char gapchar)
        {
            var mar = new int[njob];
            for (int i = 0; i < njob; i++)
                if (!alreadyaligned[i])
                    mar[i] = 3;
            foreach (int k in ex1)
                mar[k] = 1;
            foreach (int k in ex2)
                mar[k] = 2;
            var list0 = new List<int>();
            var list1 = new List<int>();
            var list2 = new List<int>();
            for (int i = 0; i < njob; i++)
            {
                if (mar[i] == 2)
                    list2.Add(i);
                if (mar[i] == 1)
                    list1.Add(i);
                if (mar[i] == 0)
                    list0.Add(i);
            }
            int ngroup0 = list0.Count, ngroup1 = list1.Count, ngroup2 = list2.Count;
            if (ngroup0 == 0)
                return;
            int rep = list0[0];
            int len = seq[rep].Length;
            int len0 = len + 1;
            var aseq = new StringBuilder[njob];
            for (int i = 0; i < njob; i++)
                aseq[i] = new StringBuilder();
            int posin12 = 0;
            for (int j = 0; j < len0; j++)
            {
                if (j < gaplen.Length && gaplen[j] != 0)
                {
                    var mseq0 = new string[ngroup0];
                    var mseq1 = new string[ngroup1];
                    var mseq2 = new string[ngroup2];
                    int gapshift = gaplen[j];
                    string gaps = new string(gapchar, gapshift);
                    for (int i = 0; i < ngroup0; i++)
                        mseq0[i] = gaps;
                    for (int i = 0; i < ngroup1; i++)
                        mseq1[i] = CopyN(seq[list1[i]], posin12, gapshift);
                    for (int i = 0; i < ngroup2; i++)
                        mseq2[i] = CopyN(seq[list2[i]], posin12, gapshift);
                    posin12 += gapshift;
                    gapshift = posin12 < gapmap.Length ? gapmap[posin12] : 0;
                    for (int i = 0; i < ngroup0; i++)
                        mseq0[i] += CopyN(seq[list0[i]], j, gapshift);
                    for (int i = 0; i < ngroup1; i++)
                        mseq1[i] += CopyN(seq[list1[i]], posin12, gapshift);
                    for (int i = 0; i < ngroup2; i++)
                        mseq2[i] += CopyN(seq[list2[i]], posin12, gapshift);
                    if (gapshift != 0)
                        ProfileAlignment(ngroup0, ngroup1, ngroup2, mseq0, mseq1, mseq2, alloclen);
                    j += gapshift;
                    posin12 += gapshift;
                    for (int i = 0; i < ngroup0; i++)
                        aseq[list0[i]].Append(mseq0[i]);
                    for (int i = 0; i < ngroup1; i++)
                        aseq[list1[i]].Append(mseq1[i]);
                    for (int i = 0; i < ngroup2; i++)
                        aseq[list2[i]].Append(mseq2[i]);
                }
                int blocklen = 1;
                for (int q = j + 1; q < len0 && q < gaplen.Length && gaplen[q] == 0; q++)
                    blocklen++;
                for (int i = 0; i < ngroup0; i++)
                    aseq[list0[i]].Append(CopyN(seq[list0[i]], j, blocklen));
                for (int i = 0; i < ngroup1; i++)
                    aseq[list1[i]].Append(CopyN(seq[list1[i]], posin12, blocklen));
                for (int i = 0; i < ngroup2; i++)
                    aseq[list2[i]].Append(CopyN(seq[list2[i]], posin12, blocklen));
                j += blocklen - 1;
                posin12 += blocklen - 1;
                posin12++;
            }
            foreach (int i in list0)
                seq[i] = aseq[i].ToString();
            foreach (int i in list1)
                seq[i] = aseq[i].ToString();
            foreach (int i in list2)
                seq[i] = aseq[i].ToString();
        }
    }

    #endregion
}
