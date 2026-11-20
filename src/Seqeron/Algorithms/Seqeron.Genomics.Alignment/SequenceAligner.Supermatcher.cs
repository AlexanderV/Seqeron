using System.Globalization;
using System.Text;
using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Alignment;

// EMBOSS 6.6 supermatcher (audit A2-28c, B09): the Rabin-Karp word seeds of nucleus/embword.c (embWordGetTable, embWordRabinKarpInit, embWordRabinKarpSearch, embWordMatchFirstMax), supermatcher.c
// (supermatcher_findendpoints) and the banded Smith-Waterman of nucleus/embaln.c (embAlignPathCalcSWFast, embAlignWalkSWMatrixFast, embAlignReportLocal), with the programs' single-precision arithmetic.
public static partial class SequenceAligner
{
    private const ulong RkModulus = 1073741789UL;
    private const ulong RkRadix = 256UL;

    /// <summary>
    /// EMBOSS <c>supermatcher</c> (local alignment of large sequences): the words of <paramref name="wordLength"/> residues of the target sequences (<c>-bsequence</c>) are looked up in every query sequence
    /// (<c>-asequence</c>) by Rabin-Karp multi-pattern search; the longest word match (the first of equal length) fixes the diagonal, and a banded Smith-Waterman of <paramref name="width"/> cells (<c>embAlignPathCalcSWFast</c>)
    /// along that diagonal, from the sequence starts to the first sequence end, gives the alignment, reported only when its score exceeds <paramref name="minScore"/>. Rows, offsets, scores and the error-file
    /// lines are those of the installed supermatcher 6.6.0 (the program is a heuristic: its alignment is the best one within the band, not necessarily the optimum).
    /// </summary>
    /// <param name="queries">The <c>-asequence</c> sequences (first member of every output pair).</param>
    /// <param name="targets">The <c>-bsequence</c> sequence set (second member of every output pair).</param>
    /// <param name="matrix">Residue-pair scores (<c>-datafile</c>); default EDNAFULL for nucleotide input, EBLOSUM62 otherwise.</param>
    /// <param name="gapOpen">Cost of a gap of length one (<c>-gapopen</c>, default 10).</param>
    /// <param name="gapExtend">Cost of each further gap position (<c>-gapextend</c>, default 0.5).</param>
    /// <param name="wordLength">Word length for the initial matching (<c>-wordlen</c>, at least 3, default 6).</param>
    /// <param name="width">Alignment band width (<c>-width</c>, at least 1, default 16).</param>
    /// <param name="minScore">Minimum score to report an alignment (<c>-minscore</c>, default 0).</param>
    /// <param name="outputFormat">The <c>-aformat</c> the result will be written in: the formats that do not show sequences (<see cref="EmbossAlignFormat.Match"/>, <see cref="EmbossAlignFormat.Score"/>) get the program's plain-rows alignment (<c>ajAlignDefineCC</c>: the TARGET is the first sequence, no offsets) instead of the sequence-range report.</param>
    public static EmbossProgramResult Supermatcher(IReadOnlyList<EmbossSequence> queries, IReadOnlyList<EmbossSequence> targets, SubstitutionMatrix? matrix = null,
        double gapOpen = 10, double gapExtend = 0.5, int wordLength = 6, int width = 16, double minScore = 0, EmbossAlignFormat? outputFormat = null)
    {
        ArgumentNullException.ThrowIfNull(queries);
        ArgumentNullException.ThrowIfNull(targets);
        queries = EmbossSequenceReader.Clean(queries);
        targets = EmbossSequenceReader.Clean(targets);
        RequireEmbossRange(gapOpen, 0, 100, nameof(gapOpen), "-gapopen", "supermatcher");
        RequireEmbossRange(gapExtend, 0, 10, nameof(gapExtend), "-gapextend", "supermatcher");
        if (wordLength < 3)
            throw new ArgumentOutOfRangeException(nameof(wordLength), "-wordlen is at least 3 in EMBOSS supermatcher.");
        if (width < 1)
            throw new ArgumentOutOfRangeException(nameof(width), "-width is at least 1 in EMBOSS supermatcher.");
        if (!(minScore >= 0 && double.IsFinite(minScore)))
            throw new ArgumentOutOfRangeException(nameof(minScore), minScore, "-minscore is a finite number, 0 or more, in EMBOSS supermatcher.");
        if (targets.Count == 0)
            throw new ArgumentException("At least one target sequence is required.", nameof(targets));
        matrix ??= queries.Concat(targets).All(s => EmbossSequenceReader.IsNucleotide(s.Residues)) ? SubstitutionMatrix.EdnaFull : SubstitutionMatrix.Blosum62;
        var table = new EmbossScoreTable(matrix);
        float gapopen = RoundFloat((float)gapOpen, 8);
        float gapextend = RoundFloat((float)gapExtend, 8);
        float minscore = (float)minScore;

        var words = new SupermatcherWords(targets, wordLength);
        var alignments = new List<EmbossPairAlignment>();
        var errors = new StringBuilder();
        foreach (var query in queries)
        {
            AlignmentCancellation.ThrowIfRequested();
            var matchLists = words.Search(query.Residues);
            for (int k = 0; k < targets.Count; k++)
            {
                var target = targets[k];
                if (matchLists[k].Count == 0)
                {
                    errors.Append(CultureInfo.InvariantCulture, $"No wordmatch start points for {query.Name} vs {target.Name}. No alignment\n");
                    continue;
                }

                var max = FirstMaxWordMatch(matchLists[k]);
                int targetstart = max.TargetStart;
                int querystart = max.QueryStart;
                int offset = targetstart - querystart;
                if (offset > 0)
                {
                    targetstart = offset;
                    querystart = 0;
                }
                else
                {
                    querystart = 0 - offset;
                    targetstart = 0;
                }

                int amax = target.Residues.Length - 1;
                int bmax = query.Residues.Length - 1;
                int targetend = targetstart;
                int queryend = querystart;
                while (targetend < amax && queryend < bmax)
                {
                    targetend++;
                    queryend++;
                }

                int lena = targetend - targetstart + 1;
                int lenb = queryend - querystart + 1;
                string a = target.Residues.Substring(targetstart, lena);
                string b = query.Residues.Substring(querystart, lenb);
                int pathwidth = width;
                RequireCells(((long)lena + 1) * pathwidth, SupermatcherMaxPathCells, nameof(SupermatcherMaxPathCells), "EMBOSS supermatcher", "keeps a banded float path matrix and an int compass matrix (8 bytes per cell)",
                    $"{lena:N0} target residues x width {pathwidth:N0} = {((long)lena + 1) * pathwidth:N0} cells", "The program has no fixed limit and dies when the allocation fails", "use a smaller -width");
                var path = new float[((long)lena + 1) * pathwidth];
                var compass = new int[path.Length];
                float score = SwPathCalcFast(a, b, lena, lenb, 0, pathwidth, gapopen, gapextend, path, compass, table);
                if (!(score > minscore))
                    continue;
                var (targetRow, queryRow) = SwWalkFast(path, compass, gapopen, gapextend, target.Residues, query.Residues, lena, lenb, 0, pathwidth, ref targetstart, ref querystart);
                if (outputFormat is EmbossAlignFormat.Match or EmbossAlignFormat.Score)
                    alignments.Add(new EmbossPairAlignment(new EmbossAlignedSequence(target.Name, targetRow.Replace('.', '-')), new EmbossAlignedSequence(query.Name, queryRow.Replace('.', '-')), score));
                else
                    alignments.Add(ReportLocal(query, target, queryRow, targetRow, querystart, targetstart, score));
            }
        }

        return new EmbossProgramResult(alignments, errors.ToString());
    }

    // embAlignReportLocal: the rows carry '.' for gaps (standardised to '-'), the residues outside the aligned part are the offsets (before) and offends (after)
    private static EmbossPairAlignment ReportLocal(EmbossSequence first, EmbossSequence second, string rowFirst, string rowSecond, int startFirst, int startSecond, float score)
    {
        static int Gaps(string row) => row.Count(c => c is '.' or '-' or ' ');
        int offend1 = first.Residues.Length - rowFirst.Length - startFirst + Gaps(rowFirst);
        int offend2 = second.Residues.Length - rowSecond.Length - startSecond + Gaps(rowSecond);
        return new EmbossPairAlignment(
            new EmbossAlignedSequence(first.Name, rowFirst.Replace('.', '-'), startFirst, offend1, Description: first.Description, Usa: first.Name),
            new EmbossAlignedSequence(second.Name, rowSecond.Replace('.', '-'), startSecond, offend2, Description: second.Description, Usa: second.Name),
            score);
    }

    private readonly record struct WordMatch(int TargetStart, int QueryStart, int Length);

    // embWordMatchFirstMax: the longest match; of equal length the one with the smaller target start (the first listed otherwise)
    private static WordMatch FirstMaxWordMatch(List<WordMatch> list)
    {
        int maxLength = 0;
        WordMatch max = default;
        foreach (var p in list)
        {
            if (p.Length > maxLength)
            {
                max = p;
                maxLength = p.Length;
            }
            else if (p.Length == maxLength && p.TargetStart < max.TargetStart)
            {
                max = p;
            }
        }

        return max;
    }

    // embWordGetTable + embWordRabinKarpInit + embWordRabinKarpSearch
    private sealed class SupermatcherWords
    {
        private sealed class Word
        {
            internal string Key = "";
            // seqlocs of the word keyed by sequence NAME (ajTablestrNew): the index of the first sequence of that name and the positions of every sequence of the name
            internal readonly Dictionary<string, (int SeqIndex, List<int> Locs)> BySequence = new(StringComparer.Ordinal);
            internal ulong Hash;
        }

        private readonly int _wl;
        private readonly IReadOnlyList<EmbossSequence> _targets;
        private readonly Word[] _sorted;

        internal SupermatcherWords(IReadOnlyList<EmbossSequence> targets, int wordLength)
        {
            _wl = wordLength;
            _targets = targets;
            var table = new Dictionary<string, Word>(StringComparer.Ordinal);
            var order = new List<Word>();
            for (int k = 0; k < targets.Count; k++)
                AddSequence(table, order, targets[k], k);
            foreach (var w in order)
            {
                ulong h = 0;
                foreach (char c in w.Key)
                    h = ((RkRadix * h) + char.ToUpperInvariant(c)) % RkModulus;
                w.Hash = h;
            }

            _sorted = [.. order.OrderBy(w => w.Hash).ThenBy(w => w.Key, StringComparer.Ordinal)];
        }

        private void AddSequence(Dictionary<string, Word> table, List<Word> order, EmbossSequence seq, int seqIndex)
        {
            string s = seq.Residues;
            int len = s.Length;
            if (len < _wl)
                return;
            char skip = EmbossSequenceReader.IsNucleotide(s) ? 'N' : 'X';
            char At(int idx) => idx >= 0 && idx < len ? char.ToUpperInvariant(s[idx]) : '\0';
            int ilast = len - _wl;
            int i = 0;
            int start = 0;
            int j = 0;
            while (j < _wl)
            {
                if (At(start + j) == skip)
                {
                    while (At(start + j) == skip)
                    {
                        i++;
                        start++;
                    }

                    j = 0;
                    if (i > ilast)
                        return;
                }
                else
                {
                    j++;
                }
            }

            j = _wl - 1;
            while (i <= ilast)
            {
                if (At(start + j) == skip)
                {
                    while (At(start + j) == skip)
                    {
                        i++;
                        start++;
                    }

                    i += j;
                    start += j;
                    continue;
                }

                string key = s.Substring(start, _wl);
                string upper = key.ToUpperInvariant();
                if (!table.TryGetValue(upper, out var rec))
                {
                    rec = new Word { Key = key };
                    table[upper] = rec;
                    order.Add(rec);
                }

                if (!rec.BySequence.TryGetValue(seq.Name, out var locs))
                {
                    locs = (seqIndex, []);
                    rec.BySequence[seq.Name] = locs;
                }

                locs.Locs.Add(i);
                start++;
                i++;
            }
        }

        // embWordRabinKarpSearch for one query: the match list of every target sequence
        internal List<WordMatch>[] Search(string text)
        {
            var lists = new List<WordMatch>[_targets.Count];
            for (int k = 0; k < lists.Length; k++)
                lists[k] = [];
            var lastLocation = new uint[_targets.Count];
            int plen = _wl;
            int tlen = text.Length;
            ulong rm = 1;
            for (int i = 1; i <= plen - 1; i++)
                rm = (RkRadix * rm) % RkModulus;
            int TextAt(int idx) => idx >= 0 && idx < tlen ? char.ToUpperInvariant(text[idx]) : 0;
            ulong textHash = 0;
            for (int i = 0; i < plen; i++)
                textHash = ((RkRadix * textHash) + (ulong)TextAt(i)) % RkModulus;
            for (int i = plen; i <= tlen; i++)
            {
                var word = FindByHash(textHash);
                if (word is not null)
                {
                    int seq2start = i - plen;
                    foreach (var (seqIndex, locs) in word.BySequence.Values)
                    {
                        string seq = _targets[seqIndex].Residues;
                        if (lastLocation[seqIndex] < (uint)i)
                        {
                            int maxloc = 0;
                            foreach (int pos in locs)
                            {
                                int matchlen = 0;
                                int ii = seq2start;
                                while (matchlen < plen)
                                {
                                    if (char.ToUpperInvariant(seq[pos + matchlen]) != TextAt(ii++))
                                        break;
                                    matchlen++;
                                }

                                if (matchlen < plen)
                                    continue;
                                while (ii < tlen && pos + matchlen < seq.Length)
                                {
                                    if (char.ToUpperInvariant(seq[pos + matchlen]) != TextAt(ii++))
                                        break;
                                    ++matchlen;
                                }

                                lists[seqIndex].Add(new WordMatch(pos, seq2start, matchlen));
                                if (ii > maxloc)
                                    maxloc = ii;
                            }

                            if (maxloc > 0)
                                lastLocation[seqIndex] = (uint)maxloc;
                        }
                    }
                }

                textHash = (((textHash + ((ulong)TextAt(i - plen) * (RkModulus - rm))) * RkRadix) + (ulong)TextAt(i)) % RkModulus;
            }

            return lists;
        }

        // bsearch over the words ordered by hash
        private Word? FindByHash(ulong hash)
        {
            int lo = 0;
            int hi = _sorted.Length - 1;
            while (lo <= hi)
            {
                int mid = (lo + hi) >> 1;
                ulong h = _sorted[mid].Hash;
                if (h == hash)
                    return _sorted[mid];
                if (h < hash)
                    lo = mid + 1;
                else
                    hi = mid - 1;
            }

            return null;
        }
    }

    // embAlignPathCalcSWFast: the single-precision banded Smith-Waterman path matrix (DIAG = 0, LEFT = 1, DOWN = 2)
    private static float SwPathCalcFast(string a, string b, int lena, int lenb, int offset, int pathwidth, float gapopen, float gapextend, float[] path, int[] compass, EmbossScoreTable table)
    {
        const int Left = 1;
        const int Down = 2;
        float ret = -float.MaxValue;
        int width = pathwidth;
        if (width < 1)
            width = 1;
        if (lena < width)
            width = lena;
        if (lenb < width)
            width = lenb;
        int leftwidth = width / 2;
        int rightwidth = width - leftwidth - 1;
        int ymin = offset >= 0 ? 0 : -offset;
        int ymax = lenb + leftwidth - offset;
        if (ymax > lena)
            ymax = lena;
        int xmin = offset <= 0 ? -leftwidth : offset - leftwidth;
        int xmax = lena + rightwidth + offset;
        if (xmax > lenb)
            xmax = lenb;
        var maxa = new float[Math.Max(xmax + 1, width)];
        float bx = 0f;
        float fnew;
        int ip;
        float Penalty(int previous, int gapKind) => previous == gapKind ? gapextend : gapopen;

        // the first column
        for (int i = 0; i <= lena; ++i)
        {
            ip = (i * width) + leftwidth - i;
            if (i < leftwidth)
            {
                path[ip] = table.Score(a[i], b[0]);
                if (path[ip] < 0)
                    path[ip] = 0;
                fnew = i == 0 ? 0f : path[ip - width - 1] - Penalty(compass[ip - width - 1], Down);
                if (fnew > path[ip])
                {
                    path[ip] = fnew;
                    compass[ip] = Down;
                }
                else
                {
                    compass[ip] = 0;
                }
            }
            else
            {
                path[ip] = 0;
                compass[ip] = 0;
            }
        }

        // the first row
        for (int j = 0; j < width; ++j)
        {
            if (j >= leftwidth)
            {
                path[j] = table.Score(a[0], b[j - leftwidth]);
                if (path[j] < 0)
                    path[j] = 0;
                fnew = j == 0 ? 0f : path[j - 1] - Penalty(compass[j - 1], Left);
                if (fnew > path[j])
                {
                    path[j] = fnew;
                    compass[j] = Left;
                }
                else
                {
                    compass[j] = 0;
                }
            }
            else
            {
                path[j] = 0;
                compass[j] = 0;
            }

            maxa[j] = path[j] - (j == 0 ? gapopen : Penalty(compass[j - 1], Left));
        }

        int irow = ymin;
        while (++irow < ymax)
        {
            AlignmentCancellation.ThrowIfRequested();       // A3-5: one poll per row of the band
            ip = (irow * width) - 1;
            int icol = xmin++;
            for (int i = 0; i < width; i++)
            {
                compass[++ip] = 0;
                if (icol++ < 0)
                    continue;
                if (icol >= xmax)
                    break;
                float match = table.Score(a[irow], b[icol]);
                float mscore = path[ip - width] + match;
                if (mscore < 0.0)
                    mscore = 0.0f;
                path[ip] = mscore;
                if (i > 0 && i != width - 1)
                {
                    fnew = path[ip - width + 1];
                    fnew -= gapopen;
                    if (i == width - 2)
                    {
                        maxa[icol] = fnew;
                    }
                    else
                    {
                        maxa[icol] -= gapextend;
                        if (fnew > maxa[icol])
                            maxa[icol] = fnew;
                    }

                    if (maxa[icol] > mscore)
                    {
                        mscore = maxa[icol];
                        path[ip] = mscore;
                        compass[ip] = Down;
                    }
                }

                if (irow > 0 && i != 0)
                {
                    if (i == 1 || icol == 1)
                    {
                        bx = path[ip - 1] - gapopen;
                    }
                    else
                    {
                        bx -= gapextend;
                        if (compass[ip - 1] != 1)
                        {
                            fnew = path[ip - 1];
                            fnew -= gapopen;
                            if (fnew > bx)
                                bx = fnew;
                        }
                    }

                    if (bx > mscore)
                    {
                        mscore = bx;
                        path[ip] = mscore;
                        compass[ip] = Left;
                    }
                }

                if (mscore > ret)
                    ret = mscore;
                if (path[ip] < 0)
                    path[ip] = 0f;
            }
        }

        return ret;
    }

    // embAlignWalkSWMatrixFast: the walk back from the first cell of the maximum; returns the two rows ('.' for gaps) and moves start1 / start2 to the first aligned residue
    private static (string Row1, string Row2) SwWalkFast(float[] path, int[] compass, float gapopen, float gapextend, string seqA, string seqB, int lena, int lenb, int offset, int width, ref int start1, ref int start2)
    {
        const int Left = 1;
        const int Down = 2;
        const float Feps = 1.192e-6F;
        static bool FpEq(double a, double b) => (b - Feps) < a && a < (b + Feps);
        if (lena < width)
            width = lena;
        if (lenb < width)
            width = lenb;
        int leftwidth = width / 2;
        int ymin = offset >= 0 ? 0 : -offset;
        int ymax = lenb + leftwidth - offset;
        if (ymax > lena)
            ymax = lena;
        double pmax = -float.MaxValue;
        int k = 0;
        int xpos = 0;
        int ypos = 0;
        for (int i = ymin; i < ymax; ++i)
        {
            for (int j = 0; j < width; ++j)
            {
                if (path[k++] > pmax)
                {
                    pmax = path[k - 1];
                    xpos = j;
                    ypos = i;
                }
            }
        }

        var m = new StringBuilder();
        var n = new StringBuilder();
        int pBase = start1;
        int qBase = start2;
        int xpos2 = xpos + ypos - leftwidth;
        while (xpos2 >= 0 && ypos >= 0 && path[(ypos * width) + xpos] > 0.0)
        {
            int ip = (ypos * width) + xpos;
            if (compass[ip] == 0)
            {
                m.Append(seqA[pBase + ypos--]);
                n.Append(seqB[qBase + xpos2--]);
                if (xpos2 >= 0 && ypos >= 0 && ip >= width && path[ip - width] <= 0.0)
                    break;
                continue;
            }

            if (compass[ip] == Left)
            {
                double score = path[ip];
                int gapcnt = 0;
                int ix = xpos - 1;
                double bimble = 0.0;
                while (ip > 0)
                {
                    bimble = path[--ip] - gapopen - (gapcnt * gapextend);
                    if (ix == 0 || FpEq(score, bimble))
                        break;
                    --ix;
                    ++gapcnt;
                }

                if (bimble <= 0.0)
                    break;
                for (int ic = 0; ic <= gapcnt; ++ic)
                {
                    m.Append('.');
                    n.Append(seqB[qBase + xpos2--]);
                    xpos--;
                }

                continue;
            }

            if (compass[ip] == Down)
            {
                double score = path[ip];
                int gapcnt = 0;
                int iy = ypos - 1;
                double bimble = 0.0;
                while (iy >= 0)
                {
                    ip -= width;
                    bimble = path[++ip] - gapopen - (gapcnt * gapextend);
                    if (FpEq(score, bimble))
                        break;
                    --iy;
                    ++gapcnt;
                }

                if (bimble <= 0.0)
                    break;
                for (int ic = 0; ic <= gapcnt; ++ic)
                {
                    m.Append(seqA[pBase + ypos--]);
                    n.Append('.');
                    xpos++;
                }

                continue;
            }

            throw new InvalidOperationException("Walk Error in SW");
        }

        start1 += ypos + 1;
        start2 += xpos2 + 1;
        return (Reverse(m), Reverse(n));
    }
}
