using System.Globalization;
using System.Text;
using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Alignment;

/// <summary>The switches of <c>combineMUMs</c> (MUMmer 3.23).</summary>
public sealed record MummerCombineMumsOptions
{
    /// <summary><c>-D</c>: write only the positions and characters of the differences to standard output instead of the alignments.</summary>
    public bool OnlyDifferencePositions { get; init; }

    /// <summary><c>-e rate</c>: the error rate cutoff of the extensions (e.g. 0.02 for two percent); the branch point score of a match becomes <c>1 - (1 - rate)</c> in single precision. Null = the built-in 0.272.</summary>
    public double? ErrorRate { get; init; }

    /// <summary><c>-n</c>: matches only between a, c, g and t (every other character of the reference becomes '2', of the query '1').</summary>
    public bool NucleotidesOnly { get; init; }

    /// <summary><c>-N num</c>: break matches at this many consecutive characters that are not a, c, g or t (0 = never).</summary>
    public int ConsecutiveNonAcgt { get; init; }

    /// <summary><c>-q tag</c>: suffix of the query tag (default <c>Query</c>).</summary>
    public string QuerySuffix { get; init; } = "Query";

    /// <summary><c>-r tag</c>: suffix of the reference tag (default <c>Ref</c>).</summary>
    public string ReferenceSuffix { get; init; } = "Ref";

    /// <summary><c>-S</c>: write every difference of the extensions and gaps.</summary>
    public bool ShowDifferences { get; init; }

    /// <summary><c>-t</c>: tag the matches with the first word of the FASTA headers.</summary>
    public bool TagFromFastaLine { get; init; }

    /// <summary><c>-v num</c>: verbose level (above 0 every extension is written with its alignment).</summary>
    public int Verbose { get; init; }

    /// <summary><c>-x</c>: do not write the <c>.cover</c> files (nor their summary).</summary>
    public bool NoCoverFiles { get; init; }

    /// <summary><c>-W file</c>: name of the gaps file with the number of errors (default <c>witherrors.gaps</c>).</summary>
    public string ErrorsFileName { get; init; } = "witherrors.gaps";

    /// <summary>The file name of the reference as given on the command line (the reference tag is its base name without extension).</summary>
    public string ReferenceName { get; init; } = "reference.fa";

    /// <summary>The file name of the queries as given on the command line (the query tag is its base name without extension).</summary>
    public string QueryName { get; init; } = "query.fa";
}

/// <summary>Output of <see cref="MummerCombineMums.Run"/>.</summary>
/// <param name="Text">What the program writes to standard output: the gaps file with the extensions, alignments and <c>Region:</c> lines.</param>
/// <param name="ErrorsFile">The file <c>-W</c> names (<c>witherrors.gaps</c>): the gap lines with the number of errors appended.</param>
/// <param name="Log">What the program writes to standard error: the reading warnings, the summary of the matches and, without <c>-x</c>, the coverage of both sequences.</param>
/// <param name="Files">The files the program leaves behind by name: the errors file and, without <c>-x</c>, <c>&lt;ref tag&gt;.cover</c> and <c>&lt;query tag&gt;.cover</c>.</param>
public sealed record MummerCombineMumsResult(string Text, string ErrorsFile, string Log, IReadOnlyDictionary<string, string> Files);

/// <summary>
/// <c>combineMUMs</c> of MUMmer 3.23 (<c>src/tigr/combineMUMs.cc</c>, Delcher): takes the clusters of exact matches of <c>mgaps</c> and extends every match off its ends and fills the gaps between the matches by edit-distance alignments
/// (diagonal bands of the Myers / Landau-Vishkin kind with a branch point score of 0.272 per match and -0.728 per error, a give-up length of 200, at most 500 errors, extensions in steps of 10 000 bases), prints the regions
/// (<c>Region:  start .. end  ...  errors / length  percent</c>) and the coverage of both sequences.
/// </summary>
/// <remarks>
/// The program reads the reference (its first record) and the query records in the order the <c>&gt;</c> lines of the gaps file ask for them (forward block, reverse block, ...); the reverse blocks use the reverse complement of the
/// query. The scoring uses single precision floats (<c>float Branch_Pt_Match_Value</c>), reproduced. A failed assertion of the program (no record, a line that is not three numbers, a gaps line before any <c>&gt;</c> line is undefined) is a
/// <see cref="FormatException"/>; so is a coordinate beyond a sequence (the program reads outside its buffers).
/// </remarks>
public static class MummerCombineMums
{
    /// <summary>Runs <c>combineMUMs [options] reference query gaps</c>.</summary>
    /// <param name="referenceFasta">The reference (only its first record is read).</param>
    /// <param name="queryFasta">The multi-FASTA file of the query records.</param>
    /// <param name="gaps">The gaps file (output of <c>mgaps</c>).</param>
    /// <param name="options">The switches (null = defaults).</param>
    /// <exception cref="FormatException">The program fails (assertion, missing record, undefined memory access).</exception>
    public static MummerCombineMumsResult Run(string referenceFasta, string queryFasta, string gaps, MummerCombineMumsOptions? options = null)
    {
        ArgumentNullException.ThrowIfNull(referenceFasta);
        ArgumentNullException.ThrowIfNull(queryFasta);
        ArgumentNullException.ThrowIfNull(gaps);
        return new Engine(options ?? new MummerCombineMumsOptions(), referenceFasta, queryFasta).Execute(gaps);
    }

    private sealed class Cover(long lo, long hi)
    {
        public long Lo = lo, Hi = hi;
        public Cover? Next;
    }

    private const int InitialBufferSize = 100000;
    private const int MaxErrors = 500, MaxExtension = 10000, GiveUpLen = 200, DefaultPad = 10, DisplayWidth = 60, MaxLine = 1024;

    private sealed class Engine
    {
        private readonly MummerCombineMumsOptions _o;
        private readonly string _referenceFasta, _queryFasta;
        private readonly StringBuilder _out = new(), _errorFile = new(), _err = new();
        private readonly int[] _editSpace = new int[(MaxErrors + 4) * MaxErrors];
        private readonly int[] _rowOffset = new int[MaxErrors];
        private float _branchPtMatchValue = 0.272f;
        private byte[] _ref = [0, 0], _query = new byte[InitialBufferSize];
        private long _refLen, _queryLen, _querySize = InitialBufferSize;
        private int _fillCt;
        private bool _padUndefined;

        public Engine(MummerCombineMumsOptions o, string referenceFasta, string queryFasta)
        {
            _o = o;
            _referenceFasta = referenceFasta;
            _queryFasta = queryFasta;
            int offset = 2, del = 6;
            for (int i = 0; i < MaxErrors; i++)
            {
                _rowOffset[i] = offset;
                offset += del;
                del += 2;
            }
        }

        private static FormatException Abort(string what) => new("combineMUMs: " + what);

        public MummerCombineMumsResult Execute(string gaps)
        {
            var ci = CultureInfo.InvariantCulture;
            if (_o.ErrorRate is double rate)
            {
                double percentId = 1.0 - rate;
                _branchPtMatchValue = (float)(1.0 - percentId);
            }

            string refTag = BaseTag(_o.ReferenceName) + _o.ReferenceSuffix;
            string queryTag = BaseTag(_o.QueryName) + _o.QuerySuffix;

            var refReader = new MummerFastaReader(_referenceFasta, ownReader: true);
            if (!ReadRecord(refReader, out string refHeader, out _ref, out _refLen))
                throw Abort("Assertion `Read_String (fp, & Ref, & Ref_Size, ref_header)' failed: the reference has no sequence record.");
            if (_o.TagFromFastaLine)
                refTag = FirstToken(refHeader, " \t\n>") ?? throw Abort("segmentation fault: the reference header has no tag (strtok returns NULL).");
            for (long i = 1; i <= _refLen; i++)
            {
                switch ((char)_ref[i])
                {
                    case 'a' or 'c' or 'g' or 't':
                        break;
                    default:
                        if (_o.NucleotidesOnly)
                            _ref[i] = (byte)'2';
                        _fillCt++;
                        break;
                }
            }

            var queryReader = new MummerFastaReader(_queryFasta, ownReader: true);
            var errorFile = _errorFile;
            long refLo = 0, refHi = 0, queryLo = 0, queryHi = 0, priorMatchLen = 0, totalErrors = 0;
            int matchCt = 0;
            double matchTotal = 0.0;
            bool isForward = false, firstMatch = true;
            Cover? refCover = null, queryCover = null;
            bool queryRead = false;

            void EndRegion(bool resetErrors)
            {
                totalErrors += ExtendForward(ref refHi, ref queryHi);
                long maxLen = 1 + refHi - refLo;
                if (1 + Math.Abs(queryHi - queryLo) > maxLen)
                    maxLen = 1 + Math.Abs(queryHi - queryLo);
                PrintRegion(refLo, refHi, queryLo, queryHi, totalErrors, maxLen, isForward);
                matchCt++;
                matchTotal += 1 + refHi - refLo;
                if (resetErrors)
                    totalErrors = 0;
                AddCoverage(ref refCover, refLo, refHi);
                if (isForward)
                    AddCoverage(ref queryCover, queryLo, queryHi);
                else
                    AddCoverage(ref queryCover, 1 + _queryLen - queryHi, 1 + _queryLen - queryLo);
            }

            foreach (string line in MummerCText.Fgets(gaps, MaxLine))
            {
                int lineLen = line.Length;
                if (line[0] == '>')
                {
                    if (!firstMatch)
                        EndRegion(true);
                    if (_o.TagFromFastaLine)
                    {
                        string? p = FirstToken(line[1..], " >\t\n");
                        if (p is null)
                            throw Abort("segmentation fault: a '>' line without a tag (strlen of NULL).");
                        queryTag = p;
                    }

                    isForward = !isForward;
                    if (isForward)
                    {
                        if (!ReadRecord(queryReader, out _, out byte[] record, out _queryLen))
                            throw Abort("Assertion `Read_String (match_fp, & Query, & query_size, header)' failed: fewer query records than forward blocks in the gaps file.");
                        StoreQuery(record);
                        queryRead = true;
                        if (_o.NucleotidesOnly)
                        {
                            for (long i = 1; i <= _queryLen; i++)
                            {
                                if (_query[i] is not ((byte)'a' or (byte)'c' or (byte)'g' or (byte)'t'))
                                    _query[i] = (byte)'1';
                            }
                        }
                    }
                    else
                    {
                        RevComplement(_query, _queryLen);
                    }

                    firstMatch = true;
                    _out.Append(line);
                    errorFile.Append(line);
                }
                else if (line[0] == '#')
                {
                    if (!firstMatch)
                        EndRegion(true);
                    firstMatch = true;
                    _out.Append(line);
                    errorFile.Append(line);
                }
                else
                {
                    var v = new long[3];
                    if (MummerCText.ScanLongs(line, v) != 3)
                        throw Abort("Assertion `sscanf (line, \"%ld %ld %ld\", & ref_pos, & query_pos, & match_len) == 3' failed: the gaps file has a line that is not three numbers.");
                    if (!queryRead)
                        throw Abort("a match line before the first '>' line: the program reads a query it never read (undefined).");
                    long refPos = v[0], queryPos = v[1], matchLen = v[2];
                    string lineNoNl = lineLen > 0 ? line[..^1] : line;
                    if (firstMatch)
                    {
                        refLo = refPos;
                        queryLo = queryPos;
                        totalErrors += ExtendBackward(ref refLo, ref queryLo);
                        if (!_o.OnlyDifferencePositions)
                            _out.Append(line);
                        errorFile.Append(lineNoNl).Append(' ').Append("-".PadLeft(7)).Append('\n');
                    }
                    else
                    {
                        int aLen = (int)(refPos - refHi - 1), bLen = (int)(queryPos - queryHi - 1);
                        var delta = new int[MaxErrors];
                        int errors = PrefixEditDist(false, refHi + 1, aLen, queryHi + 1, bLen, MaxErrors - 1, out int aEnd, out int bEnd, out bool matchToEnd, delta, out int deltaLen, false);
                        if (_o.ShowDifferences)
                            ShowDiffs(refHi + 1, aEnd, queryHi + 1, bEnd, delta, deltaLen, refHi + 1, queryHi + 1);
                        if (!_o.OnlyDifferencePositions)
                            _out.Append(line);
                        var leftPad = SetPad(true, refHi, queryHi, (int)priorMatchLen);
                        var rightPad = SetPad(false, refPos, queryPos, (int)matchLen);
                        if (aEnd == 0 && bEnd == 0)
                        {
                            if (!_o.OnlyDifferencePositions)
                                _out.Append("   No extension\n");
                            errorFile.Append(lineNoNl).Append(' ').Append("-".PadLeft(7)).Append('\n');
                        }
                        else if (_o.OnlyDifferencePositions)
                        {
                            long qStart = isForward ? queryHi + 1 : _queryLen - queryHi;
                            int bIncr = isForward ? 1 : -1;
                            DisplayDifferencePositions(refHi + 1, refHi + 1, aEnd, queryHi + 1, qStart, bEnd, bIncr, delta, deltaLen, refTag, queryTag);
                        }
                        else
                        {
                            _out.Append("     Errors = ").Append(errors.ToString(ci)).Append('\n');
                            DisplayAlignmentWithPad(refHi + 1, aEnd, queryHi + 1, bEnd, delta, deltaLen, leftPad, rightPad);
                            errorFile.Append(lineNoNl).Append(' ').Append(errors.ToString(ci).PadLeft(7)).Append('\n');
                        }

                        totalErrors += errors;
                        if (!matchToEnd)
                        {
                            refHi += aEnd;
                            queryHi += bEnd;
                            long maxLen = 1 + refHi - refLo;
                            if (1 + Math.Abs(queryHi - queryLo) > maxLen)
                                maxLen = 1 + Math.Abs(queryHi - queryLo);
                            PrintRegion(refLo, refHi, queryLo, queryHi, totalErrors, maxLen, isForward);
                            matchCt++;
                            matchTotal += 1 + refHi - refLo;
                            totalErrors = 0;
                            AddCoverage(ref refCover, refLo, refHi);
                            if (isForward)
                                AddCoverage(ref queryCover, queryLo, queryHi);
                            else
                                AddCoverage(ref queryCover, 1 + _queryLen - queryHi, 1 + _queryLen - queryLo);
                            refLo = refPos;
                            queryLo = queryPos;
                            totalErrors += ExtendBackward(ref refLo, ref queryLo);
                        }
                    }

                    refHi = refPos + matchLen - 1;
                    queryHi = queryPos + matchLen - 1;
                    priorMatchLen = matchLen;
                    firstMatch = false;
                }
            }

            if (!firstMatch)
                EndRegion(false);

            Err("           Ref len = " + _refLen.ToString(ci) + "\n");
            Err("            acgt's = " + (_refLen - _fillCt).ToString(ci) + "\n");
            Err("        Non acgt's = " + _fillCt.ToString(ci) + "\n");
            Err(" Number of matches = " + matchCt.ToString(ci) + "\n");
            Err("Sum of match bases = " + CFormat.Fixed(matchTotal, 0) + "\n");
            Err("   Avg match bases = " + CFormat.Fixed(matchCt == 0 ? 0.0 : matchTotal / matchCt, 0) + "\n");

            var files = new Dictionary<string, string> { [_o.ErrorsFileName] = errorFile.ToString() };
            if (!_o.NoCoverFiles)
            {
                files[refTag + ".cover"] = ShowCoverage(refCover, refTag, _ref);
                if (!queryRead)
                    throw Abort("no query was read: the coverage of the query reads memory the program never wrote (undefined).");
                RevComplement(_query, _queryLen);
                files[queryTag + ".cover"] = ShowCoverage(queryCover, queryTag, _query);
            }

            return new MummerCombineMumsResult(_out.ToString(), errorFile.ToString(), _err.ToString(), files);
        }

        private void Err(string s) => _err.Append(s);

        /// <summary>
        /// The query buffer of the program is one block that every record is read into again (<c>Read_String (match_fp, &amp; Query, ...)</c>): the bases from index 1, a zero byte after them, and whatever a longer earlier record left behind beyond
        /// that (the coverage of the query counts those bytes when the covered regions of earlier records reach beyond the last record). The block grows by 1.4 when a record does not fit.
        /// </summary>
        private void StoreQuery(byte[] record)
        {
            while (_queryLen >= _querySize - 1)
                _querySize = (long)(_querySize * 1.4);
            if (_querySize > _query.Length)
                Array.Resize(ref _query, (int)_querySize);
            Array.Copy(record, 1, _query, 1, _queryLen);
            _query[_queryLen + 1] = 0;
        }

        private static string BaseTag(string path)
        {
            int slash = path.LastIndexOf('/');
            string tag = slash < 0 ? path : path[(slash + 1)..];
            int dot = tag.LastIndexOf('.');
            return dot < 0 ? tag : tag[..dot];
        }

        /// <summary><c>strtok (s, delimiters)</c>: the first token, null when there is none.</summary>
        private static string? FirstToken(string s, string delimiters)
        {
            int i = 0;
            while (i < s.Length && delimiters.Contains(s[i]))
                i++;
            if (i >= s.Length)
                return null;
            int j = i;
            while (j < s.Length && !delimiters.Contains(s[j]))
                j++;
            return s[i..j];
        }

        private bool ReadRecord(MummerFastaReader reader, out string header, out byte[] seq, out long len)
        {
            var warn = new StringBuilder();
            bool ok = reader.Read(out header, out seq, out len, warn);
            _err.Append(warn);
            return ok;
        }

        private void PrintRegion(long refLo, long refHi, long queryLo, long queryHi, long totalErrors, long maxLen, bool isForward)
        {
            if (_o.OnlyDifferencePositions)
                return;
            var ci = CultureInfo.InvariantCulture;
            double error = maxLen == 0 ? 0.0 : (double)totalErrors / maxLen;
            long q1 = isForward ? queryLo : 1 + _queryLen - queryLo, q2 = isForward ? queryHi : 1 + _queryLen - queryHi;
            _out.Append("Region:  ").Append(refLo.ToString(ci).PadLeft(9)).Append(" .. ").Append(refHi.ToString(ci).PadRight(9)).Append("  ")
                .Append(q1.ToString(ci).PadLeft(9)).Append(" .. ").Append(q2.ToString(ci).PadRight(9)).Append("  ")
                .Append(totalErrors.ToString(ci).PadLeft(7)).Append(" / ").Append(maxLen.ToString(ci).PadRight(9)).Append("  ")
                .Append(CFormat.Fixed(100.0 * error, 2).PadLeft(5)).Append("%\n");
        }

        private static void RevComplement(byte[] s, long len)
        {
            long i = 1, j = len;
            for (; i < j; i++, j--)
            {
                byte ch = MummerSequence.Complement(s[i]);
                s[i] = MummerSequence.Complement(s[j]);
                s[j] = ch;
            }

            if (i == j)
                s[i] = MummerSequence.Complement(s[i]);
        }

        private static void AddCoverage(ref Cover? list, long lo, long hi)
        {
            if (list is null || hi + 1 < list.Lo)
            {
                list = new Cover(lo, hi) { Next = list };
                return;
            }

            Cover? prev = null, p = list;
            for (; p is not null && lo - 1 > p.Hi; p = p.Next)
                prev = p;
            if (p is null || hi + 1 < p.Lo)
            {
                prev!.Next = new Cover(lo, hi) { Next = p };
                return;
            }

            if (lo < p.Lo)
                p.Lo = lo;
            while (p.Next is not null && hi + 1 >= p.Next.Lo)
            {
                p.Hi = p.Next.Hi;
                p.Next = p.Next.Next;
            }

            if (hi > p.Hi)
                p.Hi = hi;
        }

        // ------------------------------------------------------------------ the edit distance

        private int ExtendBackward(ref long refLo, ref long queryLo)
        {
            int ct = 0, sum = 0;
            int aEnd, bEnd;
            do
            {
                long aLen = refLo - 1;
                if (aLen > MaxExtension)
                    aLen = MaxExtension;
                long bLen = queryLo - 1;
                if (bLen > MaxExtension)
                    bLen = MaxExtension;
                var delta = new int[MaxErrors];
                int errors = PrefixEditDist(true, refLo - 1, (int)aLen, queryLo - 1, (int)bLen, MaxErrors - 1, out aEnd, out bEnd, out bool matchToEnd, delta, out int deltaLen, true);
                if (_o.ShowDifferences)
                    RevShowDiffs(refLo - 1, aEnd, queryLo - 1, bEnd, delta, deltaLen, refLo - 1, queryLo - 1);
                if (_o.Verbose > 0)
                {
                    _out.Append(Verbose("revextend#", ++ct, refLo - 1, queryLo - 1, (int)aLen, (int)bLen, errors, matchToEnd, aEnd, bEnd));
                    DisplayAlignment(-1, refLo - 1, aEnd, queryLo - 1, bEnd, delta, deltaLen);
                }

                refLo -= aEnd;
                queryLo -= bEnd;
                sum += errors;
            }
            while (aEnd > 0.9 * MaxExtension || bEnd > MaxExtension);
            return sum;
        }

        private int ExtendForward(ref long refHi, ref long queryHi)
        {
            int ct = 0, sum = 0;
            int aEnd, bEnd;
            do
            {
                long aLen = _refLen - refHi;
                if (aLen > MaxExtension)
                    aLen = MaxExtension;
                long bLen = _queryLen - queryHi;
                if (bLen > MaxExtension)
                    bLen = MaxExtension;
                var delta = new int[MaxErrors];
                int errors = PrefixEditDist(false, refHi + 1, (int)aLen, queryHi + 1, (int)bLen, MaxErrors - 1, out aEnd, out bEnd, out bool matchToEnd, delta, out int deltaLen, true);
                if (_o.ShowDifferences)
                    ShowDiffs(refHi + 1, aEnd, queryHi + 1, bEnd, delta, deltaLen, refHi + 1, queryHi + 1);
                if (_o.Verbose > 0)
                {
                    _out.Append(Verbose("extend#", ++ct, refHi + 1, queryHi + 1, (int)aLen, (int)bLen, errors, matchToEnd, aEnd, bEnd));
                    DisplayAlignment(1, refHi + 1, aEnd, queryHi + 1, bEnd, delta, deltaLen);
                }

                refHi += aEnd;
                queryHi += bEnd;
                sum += errors;
            }
            while (aEnd > 0.9 * MaxExtension || bEnd > MaxExtension);
            return sum;
        }

        private static string Verbose(string label, int ct, long a, long b, int aLen, int bLen, int errors, bool matchToEnd, int aEnd, int bEnd)
        {
            var ci = CultureInfo.InvariantCulture;
            return label + ct.ToString(ci) + "  " + a.ToString(ci).PadLeft(9) + " " + b.ToString(ci).PadLeft(9) + "  " + aLen.ToString(ci).PadLeft(5) + " " + bLen.ToString(ci).PadLeft(5) + "  "
                + errors.ToString(ci).PadLeft(3) + "  " + (matchToEnd ? 'T' : 'F') + "  " + aEnd.ToString(ci).PadLeft(4) + " " + bEnd.ToString(ci).PadLeft(4) + "\n";
        }

        private int E(int e, int d) => _editSpace[_rowOffset[e] + d];

        private void SetE(int e, int d, int v) => _editSpace[_rowOffset[e] + d] = v;

        /// <summary>
        /// <c>Prefix_Edit_Dist</c> (forward) and <c>Rev_Prefix_Edit_Dist</c> (reverse: the strings are read from the base downwards): the minimum number of changes to match <c>A [0 .. m-1]</c> with a prefix of <c>T [0 .. n-1]</c>
        /// (<c>A [i]</c> is <c>Ref [aBase + dir * i]</c>), <paramref name="extending"/> as in the program.
        /// </summary>
        private int PrefixEditDist(bool reverse, long aBase, int m, long tBase, int n, int errorLimit, out int aEnd, out int tEnd, out bool matchToEnd, int[] delta, out int deltaLen, bool extending)
        {
            int dir = reverse ? -1 : 1;
            CheckRange(_refLen, aBase, dir, m, "reference");
            CheckRange(_queryLen, tBase, dir, n, "query");
            byte A(long i) => _ref[aBase + dir * i];
            byte T(long i) => _query[tBase + dir * i];

            double score, maxScore;
            int maxScoreLen, maxScoreBestD, maxScoreBestE;
            int bestD = 0, bestE = 0, longest = 0, row;
            int left, right;
            deltaLen = 0;

            if (_o.ConsecutiveNonAcgt > 0)
            {
                int ct = 0, i;
                for (i = 0; i < m && ct < _o.ConsecutiveNonAcgt; i++)
                    ct = A(i) is (byte)'a' or (byte)'c' or (byte)'g' or (byte)'t' ? 0 : ct + 1;
                if (ct >= _o.ConsecutiveNonAcgt)
                {
                    m = i - ct;
                    extending = true;
                }

                ct = 0;
                for (i = 0; i < n && ct < _o.ConsecutiveNonAcgt; i++)
                    ct = T(i) is (byte)'a' or (byte)'c' or (byte)'g' or (byte)'t' ? 0 : ct + 1;
                if (ct >= _o.ConsecutiveNonAcgt)
                {
                    n = i - ct;
                    extending = true;
                }
            }

            int shorter = Math.Min(m, n);
            row = 0;
            while (row < shorter && A(row) == T(row))
                row++;

            SetE(0, 0, row);
            if (row == m && row == n)
            {
                deltaLen = 0;
                aEnd = tEnd = row;
                matchToEnd = !extending;
                return 0;
            }

            left = right = 0;
            maxScore = row * _branchPtMatchValue;
            maxScoreLen = row;
            maxScoreBestD = 0;
            maxScoreBestE = 0;

            for (int e = 1; e <= errorLimit; e++)
            {
                left = Math.Max(left - 1, -e);
                right = Math.Min(right + 1, e);
                SetE(e - 1, left, -2);
                SetE(e - 1, left - 1, -2);
                SetE(e - 1, right, -2);
                SetE(e - 1, right + 1, -2);

                for (int d = left; d <= right; d++)
                {
                    row = 1 + E(e - 1, d);
                    int j = E(e - 1, d - 1);
                    if (j > row)
                        row = j;
                    j = 1 + E(e - 1, d + 1);
                    if (j > row)
                        row = j;
                    while (row < m && row + d < n && A(row) == T(row + d))
                        row++;

                    SetE(e, d, row);

                    if (row == m && row + d == n)
                    {
                        if (extending)
                        {
                            for (j = left; j <= d; j++)
                            {
                                if (E(e, j) > longest)
                                    longest = E(e, j);
                            }

                            score = longest * _branchPtMatchValue - e;
                            if (score < maxScore)
                            {
                                aEnd = maxScoreLen;
                                tEnd = maxScoreLen + maxScoreBestD;
                                SetDeltas(delta, out deltaLen, maxScoreLen, maxScoreBestD, maxScoreBestE);
                                matchToEnd = false;
                                return maxScoreBestE;
                            }
                        }

                        aEnd = row;
                        tEnd = row + d;
                        SetDeltas(delta, out deltaLen, row, d, e);
                        matchToEnd = !extending;
                        return e;
                    }
                }

                int cutoff = longest - GiveUpLen;
                while (left <= right && E(e, left) < cutoff)
                    left++;
                if (left > right)
                    break;
                while (right > left && E(e, right) < cutoff)
                    right--;

                for (int d = left; d <= right; d++)
                {
                    if (E(e, d) > longest)
                    {
                        bestD = d;
                        bestE = e;
                        longest = E(e, d);
                    }
                }

                score = longest * _branchPtMatchValue - e;
                if (score > maxScore)
                {
                    maxScore = score;
                    maxScoreLen = longest;
                    maxScoreBestD = bestD;
                    maxScoreBestE = bestE;
                }

                if (longest - maxScoreLen >= GiveUpLen)
                    break;
            }

            aEnd = maxScoreLen;
            tEnd = maxScoreLen + maxScoreBestD;
            SetDeltas(delta, out deltaLen, maxScoreLen, maxScoreBestD, maxScoreBestE);
            matchToEnd = false;
            return maxScoreBestE;
        }

        private static void CheckRange(long len, long start, int dir, int count, string name)
        {
            if (count <= 0)
                return;
            long last = start + dir * (count - 1L);
            if (start < 1 || start > len || last < 1 || last > len)
                throw new FormatException("combineMUMs: a match or extension reaches outside the " + name + " sequence; the program reads memory beyond its buffer (undefined).");
        }

        private void SetDeltas(int[] delta, out int deltaLen, int row, int d, int e)
        {
            var stack = new int[MaxErrors];
            int last = row, count = 0;
            for (int k = e; k > 0; k--)
            {
                int from = d;
                int max = 1 + E(k - 1, d);
                int j = E(k - 1, d - 1);
                if (j > max)
                {
                    from = d - 1;
                    max = j;
                }

                j = 1 + E(k - 1, d + 1);
                if (j > max)
                {
                    from = d + 1;
                    max = j;
                }

                if (from == d - 1)
                {
                    stack[count++] = max - last - 1;
                    d--;
                    last = E(k - 1, from);
                }
                else if (from == d + 1)
                {
                    stack[count++] = last - (max - 1);
                    d++;
                    last = E(k - 1, from);
                }
            }

            stack[count++] = last + 1;
            int n = 0;
            for (int i = count - 1; i > 0; i--)
                delta[n++] = Math.Abs(stack[i]) * Math.Sign(stack[i - 1]);
            deltaLen = count - 1;
        }

        // ------------------------------------------------------------------ the display

        /// <summary>The two rows of an alignment: the gaps coded in <paramref name="delta"/> over a (<c>A [i]</c> = <c>Ref [aBase + dir * i]</c>) and b.</summary>
        private (StringBuilder Top, StringBuilder Bottom) Rows(int dir, long aBase, int aLen, long bBase, int bLen, int[] delta, int deltaCt)
        {
            var top = new StringBuilder();
            var bottom = new StringBuilder();
            int i = 0, j = 0;
            for (int k = 0; k < deltaCt; k++)
            {
                for (int m = 1; m < Math.Abs(delta[k]); m++)
                {
                    top.Append((char)_ref[aBase + dir * i++]);
                    j++;
                }

                if (delta[k] < 0)
                {
                    top.Append('-');
                    j++;
                }
                else
                {
                    top.Append((char)_ref[aBase + dir * i++]);
                }
            }

            while (i < aLen && j < bLen)
            {
                top.Append((char)_ref[aBase + dir * i++]);
                j++;
            }

            i = j = 0;
            for (int k = 0; k < deltaCt; k++)
            {
                for (int m = 1; m < Math.Abs(delta[k]); m++)
                {
                    bottom.Append((char)_query[bBase + dir * j++]);
                    i++;
                }

                if (delta[k] > 0)
                {
                    bottom.Append('-');
                    i++;
                }
                else
                {
                    bottom.Append((char)_query[bBase + dir * j++]);
                }
            }

            while (j < bLen && i < aLen)
            {
                bottom.Append((char)_query[bBase + dir * j++]);
                i++;
            }

            return (top, bottom);
        }

        /// <summary><c>Display_Alignment</c> (dir 1) and <c>Rev_Display_Alignment</c> (dir -1, with an empty line before every block).</summary>
        private void DisplayAlignment(int dir, long aBase, int aLen, long bBase, int bLen, int[] delta, int deltaCt)
        {
            var (top, bottom) = Rows(dir, aBase, aLen, bBase, bLen, delta, deltaCt);
            int topLen = top.Length, bottomLen = bottom.Length;
            for (int i = 0; i < topLen || i < bottomLen; i += DisplayWidth)
            {
                if (dir < 0)
                    _out.Append('\n');
                _out.Append("A: ");
                for (int j = 0; j < DisplayWidth && i + j < topLen; j++)
                    _out.Append(top[i + j]);
                _out.Append("\nB: ");
                for (int j = 0; j < DisplayWidth && i + j < bottomLen; j++)
                    _out.Append(bottom[i + j]);
                _out.Append("\n   ");
                for (int j = 0; j < DisplayWidth && i + j < bottomLen && i + j < topLen; j++)
                    _out.Append(top[i + j] != ' ' && bottom[i + j] != ' ' && top[i + j] != bottom[i + j] ? '^' : ' ');
                _out.Append('\n');
            }
        }

        private void DisplayAlignmentWithPad(long aBase, int aLen, long bBase, int bLen, int[] delta, int deltaCt, byte[][] leftPad, byte[][] rightPad)
        {
            if (_padUndefined)
                throw new FormatException("combineMUMs: a match reaches outside a sequence; the program reads memory beyond its buffer (undefined).");
            var (top, bottom) = Rows(1, aBase, aLen, bBase, bLen, delta, deltaCt);
            for (int k = DefaultPad - 1; k >= 0; k--)
            {
                top.Insert(0, (char)leftPad[0][k]);
                bottom.Insert(0, (char)leftPad[1][k]);
            }

            for (int k = 0; k < DefaultPad; k++)
            {
                top.Append((char)rightPad[0][k]);
                bottom.Append((char)rightPad[1][k]);
            }

            int topLen = top.Length, bottomLen = bottom.Length;
            if (topLen != bottomLen)
                throw new FormatException("combineMUMs: ERROR:  Alignment length mismatch  top = " + topLen.ToString(CultureInfo.InvariantCulture) + "  bottom = " + bottomLen.ToString(CultureInfo.InvariantCulture)
                    + " (the program exits with status 255)");
            for (int i = 0; i < topLen || i < bottomLen; i += DisplayWidth)
            {
                _out.Append("A: ");
                for (int j = 0; j < DisplayWidth && i + j < topLen; j++)
                    _out.Append(top[i + j]);
                _out.Append("\nB: ");
                for (int j = 0; j < DisplayWidth && i + j < bottomLen; j++)
                    _out.Append(bottom[i + j]);
                _out.Append("\n   ");
                for (int j = 0; j < DisplayWidth && i + j < bottomLen && i + j < topLen; j++)
                {
                    if (i + j < DefaultPad)
                        _out.Append((char)leftPad[2][i + j]);
                    else if (topLen - i - j <= DefaultPad)
                        _out.Append((char)rightPad[2][i + j - (topLen - DefaultPad)]);
                    else
                        _out.Append(top[i + j] != ' ' && bottom[i + j] != ' ' && top[i + j] != bottom[i + j] ? '^' : ' ');
                }

                _out.Append('\n');
            }
        }

        /// <summary><c>Set_Left_Pad</c> / <c>Set_Right_Pad</c>: the pad rows (reference, query, '=' marks) of a match of length <paramref name="mumLen"/> ending at (left) / beginning at (right) the given positions.</summary>
        private byte[][] SetPad(bool left, long refPos, long queryPos, int mumLen)
        {
            var pad = new[] { new byte[DefaultPad], new byte[DefaultPad], new byte[DefaultPad] };
            for (int i = 0; i < DefaultPad; i++)
            {
                int at = left ? DefaultPad - i - 1 : i;
                long r = left ? refPos - i : refPos + i, q = left ? queryPos - i : queryPos + i;
                bool refOk = left ? r > 0 : r <= _refLen, queryOk = left ? q > 0 : q <= _queryLen;
                pad[0][at] = refOk ? PadByte(_ref, _refLen, r) : (byte)'.';
                pad[1][at] = queryOk ? PadByte(_query, _queryLen, q) : (byte)'.';
                pad[2][at] = i < mumLen ? (byte)'=' : (byte)' ';
            }

            return pad;
        }

        private byte PadByte(byte[] seq, long len, long index)
        {
            if (index < 1 || index > len)
            {
                _padUndefined = true;
                return 0;
            }

            return seq[index];
        }

        private void DisplayDifferencePositions(long aBase, long aStart, int aLen, long bBase, long bStart, int bLen, int bIncr, int[] delta, int deltaCt, string aTag, string bTag)
        {
            var ci = CultureInfo.InvariantCulture;
            void Line(long ap, long bp, char x, char y) =>
                _out.Append(aTag.PadRight(15)).Append(' ').Append(bTag.PadRight(15)).Append(' ').Append(((int)ap).ToString(ci).PadLeft(8)).Append(' ').Append(((int)bp).ToString(ci).PadLeft(8)).Append("  ").Append(x).Append("  ").Append(y).Append('\n');
            int i = 0, j = 0;
            for (int k = 0; k < deltaCt; k++)
            {
                for (int m = 1; m < Math.Abs(delta[k]); m++)
                {
                    char x = (char)_ref[aBase + i], y = (char)_query[bBase + j];
                    if (x != y)
                        Line(aStart + i, bStart + j * bIncr, x, y);
                    i++;
                    j++;
                }

                if (delta[k] < 0)
                {
                    Line(aStart + i, bStart + j * bIncr, '-', (char)_query[bBase + j]);
                    j++;
                }
                else
                {
                    Line(aStart + i, bIncr > 0 ? bStart + j : bStart - j + 1, (char)_ref[aBase + i], '-');
                    i++;
                }
            }

            while (i < aLen && j < bLen)
            {
                char x = (char)_ref[aBase + i], y = (char)_query[bBase + j];
                if (x != y)
                    Line(aStart + i, bStart + j * bIncr, x, y);
                i++;
                j++;
            }
        }

        private void ShowDiffs(long aBase, int aLen, long bBase, int bLen, int[] delta, int deltaCt, long aRef, long bRef)
        {
            var ci = CultureInfo.InvariantCulture;
            void Line(long x, long y, string what) => _out.Append(x.ToString(ci).PadLeft(8)).Append(' ').Append(y.ToString(ci).PadLeft(8)).Append(':').Append(what).Append('\n');
            int i = 0, j = 0;
            for (int k = 0; k < deltaCt; k++)
            {
                for (int m = 1; m < Math.Abs(delta[k]); m++)
                {
                    if (_ref[aBase + i] != _query[bBase + j])
                        Line(aRef + i, bRef + j, "  " + (char)_ref[aBase + i] + "  " + (char)_query[bBase + j]);
                    i++;
                    j++;
                }

                if (delta[k] < 0)
                {
                    Line(aRef + i - 1, bRef + j, " ins " + (char)_query[bBase + j]);
                    j++;
                }
                else
                {
                    Line(aRef + i, bRef + j - 1, " del " + (char)_ref[aBase + i]);
                    i++;
                }
            }

            while (i < aLen && j < bLen)
            {
                if (_ref[aBase + i] != _query[bBase + j])
                    Line(aRef + i, bRef + j, "  " + (char)_ref[aBase + i] + "  " + (char)_query[bBase + j]);
                i++;
                j++;
            }
        }

        private void RevShowDiffs(long aBase, int aLen, long bBase, int bLen, int[] delta, int deltaCt, long aRef, long bRef)
        {
            var ci = CultureInfo.InvariantCulture;
            void Line(long x, long y, string what) => _out.Append(x.ToString(ci).PadLeft(8)).Append(' ').Append(y.ToString(ci).PadLeft(8)).Append(':').Append(what).Append('\n');
            int i = 0, j = 0;
            for (int k = 0; k < deltaCt; k++)
            {
                for (int m = 1; m < Math.Abs(delta[k]); m++)
                {
                    if (_ref[aBase - i] != _query[bBase - j])
                        Line(aRef - i, bRef - j, "  " + (char)_ref[aBase - i] + "  " + (char)_query[bBase - j]);
                    i++;
                    j++;
                }

                if (delta[k] < 0)
                {
                    Line(aRef - i, bRef - j, " ins " + (char)_query[bBase - j]);
                    j++;
                }
                else
                {
                    Line(aRef - i, bRef - j, " del " + (char)_ref[aBase - i]);
                    i++;
                }
            }

            while (i < aLen && j < bLen)
            {
                if (_ref[aBase - i] != _query[bBase - j])
                    Line(aRef - i, bRef - j, "  " + (char)_ref[aBase - i] + "  " + (char)_query[bBase - j]);
                i++;
                j++;
            }
        }

        // ------------------------------------------------------------------ the coverage

        /// <summary><c>Show_Coverage</c>: the text of the <c>.cover</c> file; the summary goes to standard error.</summary>
        private string ShowCoverage(Cover? list, string tag, byte[] seq)
        {
            var ci = CultureInfo.InvariantCulture;
            var fp = new StringBuilder();
            int ct = 0;
            long prevHi = 0, covNs = 0, gapNs = 0;
            double covTotal = 0.0, gapTotal = 0.0;
            long seqLen = 0;
            while (seq[seqLen + 1] != 0)
                seqLen++;
            if (ReferenceEquals(seq, _ref) && _refLen % 16 == 7)
                throw new FormatException("combineMUMs: the coverage of a reference whose length is 7 modulo 16 reads the first byte after its realloc'd buffer (strlen of the reference): the length printed depends on the heap of the process (undefined behaviour; valgrind: invalid read).");
            long CountNs(long from, long to)
            {
                long n = 0;
                for (long i = from; i <= to; i++)
                {
                    if (i < 1 || i >= seq.Length)
                        throw new FormatException("combineMUMs: a covered region lies outside its sequence buffer; the program reads memory beyond it (undefined).");
                    // strchr ("acgt", seq [i]) == NULL; the terminator (0) is found in the string
                    if (!(seq[i] is (byte)'a' or (byte)'c' or (byte)'g' or (byte)'t' or (byte)0))
                        n++;
                }

                return n;
            }

            string Pct(long n, long len) => CFormat.Fixed(len == 0 ? 0.0 : 100.0 * n / len, 1).PadLeft(5);
            fp.Append("Region".PadRight(9)).Append(' ').Append("Start".PadLeft(9)).Append(' ').Append("End".PadLeft(9)).Append(' ').Append("Len".PadLeft(8)).Append(' ').Append("N's".PadLeft(8)).Append("   ").Append("%N's".PadLeft(6)).Append('\n');
            void Row(string name, long start, long end, long len, long nCt) =>
                fp.Append(name).Append(' ').Append(start.ToString(ci).PadLeft(9)).Append(' ').Append(end.ToString(ci).PadLeft(9)).Append(' ').Append(len.ToString(ci).PadLeft(8)).Append(' ').Append(nCt.ToString(ci).PadLeft(8)).Append("  ").Append(Pct(nCt, len)).Append("%\n");

            for (var p = list; p is not null; p = p.Next)
            {
                long nCt = CountNs(prevHi + 1, p.Lo - 1);
                long len = p.Lo - prevHi - 1;
                Row("gap" + ct.ToString(ci).PadRight(6), prevHi + 1, p.Lo - 1, len, nCt);
                gapTotal += len;
                gapNs += nCt;
                ct++;
                nCt = CountNs(p.Lo, p.Hi);
                len = 1 + p.Hi - p.Lo;
                Row("cov" + ct.ToString(ci).PadRight(6), p.Lo, p.Hi, len, nCt);
                covTotal += len;
                covNs += nCt;
                prevHi = p.Hi;
            }

            long lastNs = CountNs(prevHi + 1, seqLen);
            long lastLen = seqLen - prevHi;
            Row("gap" + ct.ToString(ci).PadRight(6), prevHi + 1, seqLen, lastLen, lastNs);
            gapTotal += lastLen;
            gapNs += lastNs;
            long totalNs = covNs + gapNs;

            string F0(double x) => CFormat.Fixed(x, 0);
            string F1(double x) => CFormat.Fixed(x, 1);
            double sl = seqLen, acgt = seqLen - totalNs;
            Err("\n" + tag + " Sequence Coverage:\n");
            Err("   Sequence length = " + seqLen.ToString(ci) + "\n");
            Err("            acgt's = " + (seqLen - totalNs).ToString(ci) + "\n");
            Err("        Non acgt's = " + totalNs.ToString(ci) + "\n");
            Err(" Number of regions = " + ct.ToString(ci) + "\n");
            Err("     Matched bases = " + F0(covTotal) + "  (" + F1(sl == 0.0 ? 0.0 : 100.0 * covTotal / sl) + "%, " + F1(acgt == 0.0 ? 0.0 : 100.0 * (covTotal - covNs) / acgt) + "% of acgt's)\n");
            Err("     Avg match len = " + F0(ct == 0 ? 0.0 : covTotal / ct) + "\n");
            Err("    N's in matches = " + covNs.ToString(ci) + "  (" + F1(covTotal == 0.0 ? 0.0 : 100.0 * covNs / covTotal) + "%)\n");
            Err("   Unmatched bases = " + F0(gapTotal) + "  (" + F1(sl == 0.0 ? 0.0 : 100.0 * gapTotal / sl) + "%, " + F1(acgt == 0.0 ? 0.0 : 100.0 * (gapTotal - gapNs) / acgt) + "% of acgt's)\n");
            Err("       Avg gap len = " + F0(gapTotal / (1.0 + ct)) + "\n");
            Err("       N's in gaps = " + gapNs.ToString(ci) + "  (" + F1(gapTotal == 0.0 ? 0.0 : 100.0 * gapNs / gapTotal) + "%)\n");
            return fp.ToString();
        }
    }
}
