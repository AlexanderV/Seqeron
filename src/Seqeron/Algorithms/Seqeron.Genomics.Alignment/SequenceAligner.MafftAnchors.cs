#pragma warning disable S127, S3358, S1066, S3776, S107, S1854, IDE0059, S2325, S138, S1541
using System.Globalization;
using System.Text;

namespace Seqeron.Genomics.Alignment;

public static partial class SequenceAligner
{
    #region MAFFT 7.505: --anchors (disttbfast.c: readexternalanchors, pickpairanch, checkanchors_strongestfirst; Falign.c: Falign_givenanchors)

    /// <summary><c>ExtAnch</c> of <c>disttbfast.c</c>: a stretch of two sequences (0-based residue positions) that is aligned without gaps, with its score; <see cref="StartI"/> of -1 marks an anchor that was dropped.</summary>
    internal sealed class MafftExtAnch
    {
        internal int I, J, StartI, EndI, StartJ, EndJ, Score;

        internal MafftExtAnch Copy() => (MafftExtAnch)MemberwiseClone();
    }

    internal sealed partial class MafftGlobals
    {
        /// <summary><c>useexternalanchors</c> (<c>-l</c>), <c>maxanchorseparation</c> (<c>-x n</c>) of <c>disttbfast</c> and the text of the file <c>_externalanchors</c>.</summary>
        internal bool UseExternalAnchors;
        internal double MaxAnchorSeparation = 1000.0;
        internal string? ExternalAnchorsText;
    }

    internal sealed partial class MafftRun
    {
        private const string AnchorBehindEndMessage =
            "An anchor ends behind the last residue of its sequence: recountpositions reads the table of positions (a calloc of one int more than the alignment length) behind its end, which is undefined behaviour.";

        private List<MafftExtAnch>? _extAnch;
        private int[][]? _anchIndex;

        /// <summary>The lines of <c>disttbfast</c> that read and index the anchors (before <c>constants</c>): "reading anchors", <c>readexternalanchors</c>, "sorting anchors", <c>indexanchors</c>.</summary>
        private void PrepareExternalAnchors(int njob, string[] seq)
        {
            var nogaplen = new int[njob];
            for (int i = 0; i < njob; i++)
                nogaplen[i] = GapPick0(seq[i]).Length;
            Report("reading anchors\n");
            var list = ReadExternalAnchors(njob, nogaplen);
            _anchIndex = new int[njob][];
            for (int i = 0; i < njob; i++)
            {
                _anchIndex[i] = new int[njob];
                Array.Fill(_anchIndex[i], -1);
            }

            Report("sorting anchors\n");
            var sorted = list.ToArray();
            MergeSortLikeGlibc(sorted, (p, q) => p.I != q.I ? p.I - q.I : p.J - q.J);
            _extAnch = sorted.ToList();
            for (int n = 0; n < _extAnch.Count; n++)
            {
                if (_anchIndex[_extAnch[n].I][_extAnch[n].J] == -1)
                    _anchIndex[_extAnch[n].I][_extAnch[n].J] = n;
            }
        }

        /// <summary><c>readexternalanchors</c> (<c>io.c</c>): one anchor per line of the file <c>_externalanchors</c> ("i j starti endi startj endj score", all counted from 1, i &lt; j).</summary>
        private List<MafftExtAnch> ReadExternalAnchors(int nseq, int[] nogaplen)
        {
            var list = new List<MafftExtAnch>();
            string text = G.ExternalAnchorsText ?? throw new MafftProgramFailureException("Cannot open _externalanchors\n");
            int lineno = 0;
            foreach (string buf in SplitFgetsLines(text, 10001))
            {
                lineno++;
                if (buf.Length == 0 || buf[^1] != '\n')
                    break; // fgets read the last line without its line break and set the end of file: the program stops reading there
                if (buf[0] == '#')
                    continue;
                var sc = new CScanner(buf);
                if (!sc.Int(out int i) || !sc.Int(out int j) || !sc.Int(out int starti) || !sc.Int(out int endi) || !sc.Int(out int startj) || !sc.Int(out int endj) || !sc.Int(out int score))
                    throw new NotSupportedException("A line of the file of anchors has fewer than seven numbers: the program keeps the values of the fields that sscanf did not set (uninitialised memory).");
                var a = new MafftExtAnch { I = i - 1, J = j - 1, StartI = starti - 1, EndI = endi - 1, StartJ = startj - 1, EndJ = endj - 1, Score = score };
                if (a.I >= nseq || a.J >= nseq)
                    throw new MafftProgramFailureException("\nOut of range?  The input file has " + nseq.ToString(CultureInfo.InvariantCulture) + " sequences but pair " + (a.I + 1).ToString(CultureInfo.InvariantCulture) + "-" + (a.J + 1).ToString(CultureInfo.InvariantCulture)
                        + " was specified in line " + lineno.ToString(CultureInfo.InvariantCulture) + ".\nNote that sequence IDs are counted from 1.\n");
                if (a.I >= a.J)
                    throw new MafftProgramFailureException("\nFormat problem?  \"" + (a.I + 1).ToString(CultureInfo.InvariantCulture) + " " + (a.J + 1).ToString(CultureInfo.InvariantCulture) + "\" in line " + lineno.ToString(CultureInfo.InvariantCulture)
                        + ".\nThe sequence id of the first column must be less than the second.\n");
                if (a.I < 0 || a.J < 0)
                    throw new NotSupportedException("An anchor names the sequence 0: the program reads nogaplen[-1].");
                if (a.StartI > nogaplen[a.I])
                    throw new MafftProgramFailureException("\nOut of range?  len(seq" + (a.I + 1).ToString(CultureInfo.InvariantCulture) + ")=" + nogaplen[a.I].ToString(CultureInfo.InvariantCulture) + ", but anchor=" + (a.StartI + 1).ToString(CultureInfo.InvariantCulture)
                        + " in line " + lineno.ToString(CultureInfo.InvariantCulture) + ".\nNote that position is counted from 1.\n");
                if (a.StartJ > nogaplen[a.J])
                    throw new MafftProgramFailureException("\nOut of range?  len(seq" + a.J.ToString(CultureInfo.InvariantCulture) + ")=" + (nogaplen[a.J] + 1).ToString(CultureInfo.InvariantCulture) + ", but anchor=" + (a.StartJ + 1).ToString(CultureInfo.InvariantCulture)
                        + " in line " + lineno.ToString(CultureInfo.InvariantCulture) + ".\nNote that position is counted from 1.\n");
                list.Add(a);
            }

            return list;
        }

        /// <summary><c>pickpairanch</c>: the anchors of the pairs (member of group 1, member of group 2) in the coordinates of the two groups' rows (the position of the residue in the row), the weakest ones given up (<c>checkanchors_strongestfirst</c>).</summary>
        private List<MafftExtAnch> PickPairAnch(int n1, int n2, int[] m1, int[] m2, string[] seq1, string[] seq2)
        {
            var pair = new List<MafftExtAnch>();
            var ext = _extAnch!;
            for (int i = 0; i < n1; i++)
            {
                for (int j = 0; j < n2; j++)
                {
                    if (m1[i] < m2[j])
                    {
                        int k = _anchIndex![m1[i]][m2[j]];
                        while (k != -1 && k < ext.Count && ext[k].I == m1[i] && ext[k].J == m2[j])
                        {
                            pair.Add(new MafftExtAnch { I = i, J = j, StartI = ext[k].StartI, EndI = ext[k].EndI, StartJ = ext[k].StartJ, EndJ = ext[k].EndJ, Score = ext[k].Score });
                            k++;
                        }
                    }
                    else
                    {
                        int k = _anchIndex![m2[j]][m1[i]];
                        while (k != -1 && k < ext.Count && ext[k].I == m2[j] && ext[k].J == m1[i])
                        {
                            pair.Add(new MafftExtAnch { I = i, J = j, StartI = ext[k].StartJ, EndI = ext[k].EndJ, StartJ = ext[k].StartI, EndJ = ext[k].EndI, Score = ext[k].Score });
                            k++;
                        }
                    }
                }
            }

            RecountPositions(pair, n1, n2, seq1, seq2);
            return CheckAnchorsStrongestFirst(pair, GapNonGapRatio(n1, seq1), GapNonGapRatio(n2, seq2));
        }

        /// <summary><c>recountpositions</c>: the residue positions of the anchors become the columns of the rows of the groups.</summary>
        private static void RecountPositions(List<MafftExtAnch> pair, int n1, int n2, string[] seq1, string[] seq2)
        {
            var map = new int[seq1[0].Length + 1];
            for (int k = 0; k < n1; k++)
            {
                int pos = 0;
                string s = seq1[k];
                for (int i = 0; i < seq1[0].Length + 1; i++)
                {
                    if (i == s.Length || s[i] != '-')
                    {
                        map[pos] = i;
                        pos++;
                    }
                }

                foreach (var a in pair)
                {
                    if (a.I == k)
                    {
                        if (a.StartI >= map.Length || a.EndI >= map.Length || a.StartI < 0 || a.EndI < 0)
                            throw new NotSupportedException(AnchorBehindEndMessage);
                        a.StartI = map[a.StartI];
                        a.EndI = map[a.EndI];
                    }
                }
            }

            map = new int[seq2[0].Length + 1];
            for (int k = 0; k < n2; k++)
            {
                int pos = 0;
                string s = seq2[k];
                for (int i = 0; i < seq2[0].Length + 1; i++)
                {
                    if (i == s.Length || s[i] != '-')
                    {
                        map[pos] = i;
                        pos++;
                    }
                }

                foreach (var a in pair)
                {
                    if (a.J == k)
                    {
                        if (a.StartJ >= map.Length || a.EndJ >= map.Length || a.StartJ < 0 || a.EndJ < 0)
                            throw new NotSupportedException(AnchorBehindEndMessage);
                        a.StartJ = map[a.StartJ];
                        a.EndJ = map[a.EndJ];
                    }
                }
            }
        }

        /// <summary><c>gapnongapratio</c>: the number of gaps between the first and the last residue of every row over the number of residues.</summary>
        private static double GapNonGapRatio(int n, string[] s)
        {
            double fv = 0.0, ng = 0.0;
            for (int i = 0; i < n; i++)
            {
                string t = s[i].Trim('-');
                foreach (char c in t)
                {
                    if (c == '-')
                        fv += 1.0;
                    else
                        ng += 1.0;
                }
            }

            return fv / ng;
        }

        /// <summary><c>checkanchors_strongestfirst</c>: the anchors in the order of their score (the earlier one of equal score last: glibc's merge sort puts the right one first when the tie is broken by the address), an anchor that is too far from the nearest stronger one or crosses a stronger one is dropped (<c>starti</c> = -1), anchors that touch lose one position.</summary>
        private List<MafftExtAnch> CheckAnchorsStrongestFirst(List<MafftExtAnch> list, double gapratio1, double gapratio2)
        {
            var a = list.ToArray();
            MergeSortLikeGlibc(a, (p, q) => p.Score != q.Score ? q.Score - p.Score : 1);
            int s = a.Length;
            for (int p = 0; p < s; p++)
            {
                if (a[p].StartI == -1)
                    continue;
                double zurei, zurej;
                if (p != 0)
                {
                    int nearest = -1, mindist = 999999999;
                    for (int q = 0; q < p; q++)
                    {
                        if (a[q].StartI == -1)
                            continue;
                        if (Math.Abs(a[p].StartI - a[q].StartI) < mindist)
                        {
                            nearest = q;
                            mindist = Math.Abs(a[p].StartI - a[q].StartI);
                        }
                    }

                    if (a[nearest].StartI < a[p].StartI)
                    {
                        zurei = (double)(a[p].StartI - a[nearest].EndI) / (1.0 + gapratio1);
                        zurej = (double)(a[p].StartJ - a[nearest].EndJ) / (1.0 + gapratio2);
                    }
                    else
                    {
                        zurei = (double)(a[nearest].StartI - a[p].EndI) / (1.0 + gapratio1);
                        zurej = (double)(a[nearest].StartJ - a[p].EndJ) / (1.0 + gapratio2);
                    }
                }
                else
                {
                    zurei = zurej = 0.0;
                }

                if (Math.Abs(zurei - zurej) > G.MaxAnchorSeparation)
                {
                    a[p].StartI = a[p].StartJ = a[p].EndJ = -1;
                    continue;
                }

                for (int q = p + 1; q < s; q++)
                {
                    if (a[q].StartI == -1)
                        continue;
                    if (a[p].EndI < a[q].StartI && a[p].EndJ < a[q].StartJ)
                        continue; // q lies behind the stronger anchor p
                    if (a[p].EndI == a[q].StartI && a[p].EndJ < a[q].StartJ && a[q].StartI < a[q].EndI)
                        a[q].StartI += 1;
                    else if (a[p].EndI < a[q].StartI && a[p].EndJ == a[q].StartJ && a[q].StartJ < a[q].EndJ)
                        a[q].StartJ += 1;
                    else if (!(a[q].EndI < a[p].StartI && a[q].EndJ < a[p].StartJ)) // q lies in front of p: nothing to do
                    {
                        if (a[q].EndI == a[p].StartI && a[q].EndJ < a[p].StartJ && a[q].StartI < a[q].EndI)
                            a[q].EndI -= 1;
                        else if (a[q].EndI < a[p].StartI && a[q].EndJ == a[p].StartJ && a[q].StartJ < a[q].EndJ)
                            a[q].EndJ -= 1;
                        else
                            a[q].StartI = a[q].StartJ = a[q].EndJ = -1;
                    }
                }
            }

            MergeSortLikeGlibc(a, (p, q) => p.StartI != q.StartI ? p.StartI - q.StartI : -1);
            return a.ToList();
        }

        private static double EstimateGapFreq(int n, string[] s)
        {
            double fv = 0.0;
            for (int i = 0; i < n; i++)
            {
                string t = s[i].Trim('-');
                double f = t.Count(c => c == '-');
                fv += f / t.Length;
            }

            return fv / n;
        }

        private static int TerminalMargin(int lshorter, double groupsizefac) => (int)((lshorter * 2.0 + 100) * groupsizefac);

        /// <summary><c>Falign_givenanchors</c>: the groups are cut at the anchors (an anchor of the same length in both is copied, the stretches between them are aligned by <c>A__align</c> with the boundary gaps of the neighbours); stretches of an anchor longer than 100 are cut every 100 columns.</summary>
        private void FalignGivenAnchors(List<MafftExtAnch> pairanch, double[][] dyn, string[] seq1, string[] seq2, double[] eff1, double[] eff2, int clus1, int clus2, int alloclen)
        {
            int len1 = seq1[0].Length, len2 = seq2[0].Length;
            double marginfac1 = 1.0 + EstimateGapFreq(clus1, seq1);
            double marginfac2 = 1.0 + EstimateGapFreq(clus2, seq2);
            int starttermcut1 = 0, starttermcut2 = 0, endtermcut1 = 0, endtermcut2 = 0;
            var cut1 = new List<int> { 0 };
            var cut2 = new List<int> { 0 };
            var alignorcopy = new List<char> { 'a' };
            void Put(int idx, int c1, int c2, char ac)
            {
                while (cut1.Count <= idx)
                {
                    cut1.Add(0);
                    cut2.Add(0);
                    alignorcopy.Add('\0');
                }

                cut1[idx] = c1;
                cut2[idx] = c2;
                alignorcopy[idx] = ac;
            }

            int count = 0;
            int count0 = 0;
            while (count0 < pairanch.Count)
            {
                var a = pairanch[count0];
                if (a.StartI == -1)
                {
                    count0++;
                    continue;
                }

                if (count + 2 > 100000 - 3)
                    throw new MafftProgramFailureException("TOO MANY SEGMENTS.\n");
                if (count == 0)
                {
                    if (a.StartI > TerminalMargin(a.StartJ, marginfac1))
                    {
                        int cutadd = a.StartI - TerminalMargin(a.StartJ, marginfac1);
                        Put(1, cutadd, 0, 'A');
                        count += 1;
                        starttermcut1 = 1;
                    }
                    else if (a.StartJ > TerminalMargin(a.StartI, marginfac2))
                    {
                        int cutadd = a.StartJ - TerminalMargin(a.StartI, marginfac2);
                        Put(1, 0, cutadd, 'A');
                        count += 1;
                        starttermcut2 = 1;
                    }
                }

                Put(count + 1, a.StartI, a.StartJ, 'c');
                count += 1;
                if (a.EndI - cut1[count] == a.EndJ - cut2[count])
                {
                    while (a.EndI + 1 - cut1[count] > 100 && a.EndJ + 1 - cut2[count] > 100)
                    {
                        Report("added an anchor, because the length is " + (a.EndI + 1 - cut1[count]).ToString(CultureInfo.InvariantCulture) + "," + (a.EndJ + 1 - cut2[count]).ToString(CultureInfo.InvariantCulture) + " > 100\n");
                        Put(count + 1, cut1[count] + 100, cut2[count] + 100, 'c');
                        count += 1;
                    }
                }

                Put(count + 1, a.EndI + 1, a.EndJ + 1, 'a');
                count += 1;
                count0++;
            }

            Put(count, cut1[count], cut2[count], 'a');
            if (count > 1 && len1 - cut1[count] > TerminalMargin(len2 - cut2[count], marginfac1))
            {
                int cutadd = len1 - 1 - (len1 - cut1[count] - TerminalMargin(len2 - cut2[count], marginfac1));
                alignorcopy[count] = 'A';
                Put(count + 1, cutadd, len2, 'a');
                Put(count + 2, len1, len2, 'c');
                count += 1;
                endtermcut1 = 1;
            }
            else if (count > 1 && len2 - cut2[count] > TerminalMargin(len1 - cut1[count], marginfac2))
            {
                int cutadd = len2 - 1 - (len2 - cut2[count] - TerminalMargin(len1 - cut1[count], marginfac2));
                alignorcopy[count] = 'A';
                Put(count + 1, len1, cutadd, 'a');
                Put(count + 2, len1, len2, 'c');
                count += 1;
                endtermcut2 = 1;
            }

            if (cut1[count] != len1 || cut2[count] != len2)
            {
                Put(count + 1, len1, len2, 'c');
                count += 1;
            }

            count += 1;
            var result1 = new StringBuilder[clus1];
            var result2 = new StringBuilder[clus2];
            for (int j = 0; j < clus1; j++)
                result1[j] = new StringBuilder();
            for (int j = 0; j < clus2; j++)
                result2[j] = new StringBuilder();
            int totallen = 0;
            for (int i = 0; i < count - 1; i++)
            {
                int headgp = i == 0 || (i == 1 && alignorcopy[1] == 'A') ? G.Outgap : 1;
                int tailgp = i == count - 2 || (i == count - 3 && alignorcopy[count - 3] == 'A') ? G.Outgap : 1;
                var sgap1 = new char[clus1];
                var egap1 = new char[clus1];
                var sgap2 = new char[clus2];
                var egap2 = new char[clus2];
                if (cut1[i] != 0)
                    for (int j = 0; j < clus1; j++)
                        sgap1[j] = seq1[j][cut1[i] - 1];
                else
                    Array.Fill(sgap1, 'o');
                if (cut2[i] != 0)
                    for (int j = 0; j < clus2; j++)
                        sgap2[j] = seq2[j][cut2[i] - 1];
                else
                    Array.Fill(sgap2, 'o');
                if (cut1[i + 1] != len1)
                    for (int j = 0; j < clus1; j++)
                        egap1[j] = seq1[j][cut1[i + 1]];
                else
                    Array.Fill(egap1, 'o');
                if (cut2[i + 1] != len2)
                    for (int j = 0; j < clus2; j++)
                        egap2[j] = seq2[j][cut2[i + 1]];
                else
                    Array.Fill(egap2, 'o');
                Report("DP " + (i + 1).ToString("D5", CultureInfo.InvariantCulture) + " / " + (count - 1).ToString("D5", CultureInfo.InvariantCulture) + " " + new string('\b', 17));
                var tmpres1 = new string[clus1];
                var tmpres2 = new string[clus2];
                for (int j = 0; j < clus1; j++)
                    tmpres1[j] = seq1[j].Substring(cut1[i], cut1[i + 1] - cut1[i]);
                int orilen1 = tmpres1[0].Length;
                for (int j = 0; j < clus2; j++)
                    tmpres2[j] = seq2[j].Substring(cut2[i], cut2[i + 1] - cut2[i]);
                int orilen2 = tmpres2[0].Length;
                bool copy = alignorcopy[i] == 'c' && orilen1 == orilen2 && orilen1 == tmpres1[0].Length && string.Equals(tmpres1[0], tmpres2[0], StringComparison.Ordinal);
                if (copy)
                {
                    // the two stretches are identical: copied, not aligned
                }
                else if (alignorcopy[i] == 'A')
                {
                    int penaltyExx = G.PenaltyEx == 0 ? -179 : G.PenaltyEx;
                    _ = AAlign(dyn, G.Penalty, penaltyExx, tmpres1, tmpres2, eff1, eff2, clus1, clus2, alloclen, new string(sgap1), new string(sgap2), new string(egap1), new string(egap2), headgp, tailgp, null, -1, -1, 0, 1.0, 1.0);
                }
                else if (G.Alg == 'A')
                {
                    _ = AAlign(dyn, G.Penalty, G.PenaltyEx, tmpres1, tmpres2, eff1, eff2, clus1, clus2, alloclen, new string(sgap1), new string(sgap2), new string(egap1), new string(egap2), headgp, tailgp, null, -1, -1, 0, 1.0, 1.0);
                }
                else if (G.Alg == 'M')
                {
                    _ = MsAlign(dyn, tmpres1, tmpres2, eff1, eff2, clus1, clus2, alloclen, new string(sgap1), new string(sgap2), new string(egap1), new string(egap2), headgp, tailgp, null, -1, -1, 0, 1.0, 1.0);
                }
                else
                {
                    throw new MafftProgramFailureException("alg = " + G.Alg + "\nERROR IN SOURCE FILE Falign.c\n");
                }

                if (i == 1 && alignorcopy[i] == 'A')
                {
                    if (starttermcut1 != 0)
                    {
                        if (Enumerable.Range(0, clus2).Any(j => tmpres2[j][0] != '-'))
                            Report("There may be a problem in a hard-coded parameter (1).  Please contact katoh@ifrec.osaka-u.ac.jp\n");
                    }
                    else if (starttermcut2 != 0)
                    {
                        if (Enumerable.Range(0, clus1).Any(j => tmpres1[j][0] != '-'))
                            Report("There may be a problem in a hard-coded parameter (2).  Please contact katoh@ifrec.osaka-u.ac.jp\n");
                    }
                    else
                    {
                        break;
                    }
                }
                else if (alignorcopy[i] == 'A')
                {
                    if (endtermcut1 != 0)
                    {
                        int tmplen = tmpres2[0].Length;
                        if (Enumerable.Range(0, clus2).Any(j => tmpres2[j][tmplen - 1] != '-'))
                            Report("There may be a problem in a hard-coded parameter (3).  Please contact katoh@ifrec.osaka-u.ac.jp\n");
                    }
                    else if (endtermcut2 != 0)
                    {
                        int tmplen = tmpres1[0].Length;
                        if (Enumerable.Range(0, clus1).Any(j => tmpres1[j][tmplen - 1] != '-'))
                            Report("There may be a problem in a hard-coded parameter (4).  Please contact katoh@ifrec.osaka-u.ac.jp\n");
                    }
                    else
                    {
                        break;
                    }
                }

                int nlen = tmpres1[0].Length;
                if (totallen + nlen > alloclen)
                    throw new MafftProgramFailureException("totallen=" + totallen.ToString(CultureInfo.InvariantCulture) + " +  nlen=" + nlen.ToString(CultureInfo.InvariantCulture) + " > alloclen = " + alloclen.ToString(CultureInfo.InvariantCulture) + "\nLENGTH OVER in Falign\n \n");
                for (int j = 0; j < clus1; j++)
                    result1[j].Append(tmpres1[j]);
                for (int j = 0; j < clus2; j++)
                    result2[j].Append(tmpres2[j]);
                totallen += nlen;
            }

            for (int j = 0; j < clus1; j++)
                seq1[j] = result1[j].ToString();
            for (int j = 0; j < clus2; j++)
                seq2[j] = result2[j].ToString();
        }
    }

    #endregion
}
