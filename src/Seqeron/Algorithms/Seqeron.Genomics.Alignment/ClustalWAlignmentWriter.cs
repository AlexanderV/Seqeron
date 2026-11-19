using System.Globalization;
using System.Text;
using System.Text.RegularExpressions;
using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Alignment;

/// <summary>
/// The alignment files of ClustalW 2.1 (<c>AlignmentOutput::createAlignmentOutput</c>: <c>clustalOut</c>, <c>gcgOut</c>, <c>gdeOut</c>, <c>phylipOut</c>, <c>nbrfOut</c>, <c>nexusOut</c>, <c>fastaOut</c>), byte for byte: the
/// layout, the line widths (<c>PAGEWIDTH</c> 80, <c>LINELENGTH</c> 60, <c>GCG_LINELENGTH</c> 50), the conservation line of the CLUSTAL format, the GCG checksums and weights, the <c>-RANGE</c> / <c>-SEQNO_RANGE</c> /
/// <c>-SEQNOS</c> / <c>-CASE</c> options and the order of the rows (<see cref="ClustalWAlignment.OutputIndex"/>).
/// </summary>
public static class ClustalWAlignmentWriter
{
    private const int PageWidth = 80;
    private const int LineLength = 60;
    private const int GcgLineLength = 50;
    private const int MaxNamesToDisplay = 30;
    private const int MinNamesToDisplay = 10;
    private const int ClusSequenceOffset = 16;   // clusSequenceOffset(15 + 1)
    private const int ClusSecStructOffset = 9;   // clusSecStructOffset

    private static readonly string[] StrongGroups = ["STA", "NEQK", "NHQK", "NDEQ", "QHRK", "MILV", "MILF", "HY", "FYW"];
    private static readonly string[] WeakGroups = ["CSA", "ATV", "SAG", "STNK", "STPA", "SGND", "SNDEQK", "NDEQHK", "NEQHRK", "FVLIM", "HFY"];

    /// <summary>The alignment file of <paramref name="alignment"/> in the format of <paramref name="options"/>.</summary>
    /// <param name="alignment">The alignment (rows in input order, <see cref="ClustalWAlignment.OutputIndex"/> = the order of the file).</param>
    /// <param name="options">Format and format options.</param>
    /// <param name="structure">The structure parameters (<c>-HELIXENDIN</c>, <c>-STRANDENDIN</c> shape the structure lines).</param>
    public static string Write(ClustalWAlignment alignment, ClustalWOutputOptions? options = null, ClustalWStructureOptions? structure = null)
    {
        ArgumentNullException.ThrowIfNull(alignment);
        options ??= new ClustalWOutputOptions();
        var ctx = new Context(alignment, options, structure ?? new ClustalWStructureOptions());
        return options.Format switch
        {
            ClustalWOutputFormat.Clustal => ctx.Clustal(),
            ClustalWOutputFormat.Gcg => ctx.Gcg(),
            ClustalWOutputFormat.Gde => ctx.Gde(),
            ClustalWOutputFormat.Phylip => ctx.Phylip(),
            ClustalWOutputFormat.Pir => ctx.Pir(),
            ClustalWOutputFormat.Nexus => ctx.Nexus(),
            ClustalWOutputFormat.Fasta => ctx.Fasta(),
            _ => throw new ArgumentOutOfRangeException(nameof(options), options.Format, "Unknown alignment format."),
        };
    }

    /// <summary>The default file extension of a format when <c>-OUTFILE</c> is not given (<c>AlignmentOutput::openAlignmentOutput</c>).</summary>
    public static string DefaultExtension(ClustalWOutputFormat format) => format switch
    {
        ClustalWOutputFormat.Clustal => "aln",
        ClustalWOutputFormat.Gcg => "msf",
        ClustalWOutputFormat.Gde => "gde",
        ClustalWOutputFormat.Phylip => "phy",
        ClustalWOutputFormat.Pir => "pir",
        ClustalWOutputFormat.Nexus => "nxs",
        ClustalWOutputFormat.Fasta => "fasta",
        _ => throw new ArgumentOutOfRangeException(nameof(format)),
    };

    private sealed class Context
    {
        private readonly ClustalWAlignment _a;
        private readonly ClustalWOutputOptions _o;
        private readonly int _n;
        private readonly int _length;
        private readonly int _firstRes;
        private readonly int _lastRes;
        private readonly int _maxNames;
        private readonly string[] _rows;
        private readonly ClustalWStructureOptions _st;

        public Context(ClustalWAlignment a, ClustalWOutputOptions o, ClustalWStructureOptions st)
        {
            _a = a;
            _o = o;
            _st = st;
            _n = a.Rows.Count;
            if (_n == 0)
                throw new ArgumentException("The alignment has no sequences.", nameof(a));
            if (a.OutputIndex.Count != _n || a.Names.Count != _n)
                throw new ArgumentException("Names, rows and output index must have one entry per sequence.", nameof(a));
            _rows = [.. a.Rows];
            _length = _rows.Max(r => r.Length);
            _maxNames = Math.Max(MinNamesToDisplay, a.Names.Max(x => x.Length));

            // AlignmentOutput::createAlignmentOutput: -RANGE
            _firstRes = 1;
            _lastRes = _length;
            if (o.Range is { } range)
            {
                if (range.From > range.To || range.From == -1 || range.To == -1)
                {
                    // "seqrange numbers are not set properly, using default...."
                }
                else
                {
                    if (range.From < 1)
                        throw new ArgumentException("A range starting before the first column reads outside the sequences (undefined in ClustalW).", nameof(o));
                    _firstRes = range.From;
                    _lastRes = Math.Min(range.To, _length);
                    if (_firstRes > _lastRes)
                        throw new ArgumentException("The range lies behind the end of the alignment.", nameof(o));
                }
            }
        }

        private int Index(int ii) => _a.OutputIndex[ii];

        // the residue letter of row i (0-based) at column j (1-based), '-' for a gap; null beyond the end of the row (val == -3)
        private char? Cell(int i, int j) => j <= _rows[i].Length ? _rows[i][j - 1] : null;

        private static char Letter(char? c, char gap = '-') => c is null || c == '-' ? gap : c.Value;

        private bool Dna => !_a.IsProtein;

        private static string NameOnly(string s)
        {
            int slash = s.IndexOf('/', StringComparison.Ordinal);
            return slash < 0 ? s : s[..slash];
        }

        // AlignmentOutput::findRangeValues for the sequence at output position ii (0-based): the first and the last residue number of the written columns
        private (int Start, int End) FindRangeValues(int ii)
        {
            int i = Index(ii);
            int iStart = 0, tmpStart = 0, ngaps = 0, pregaps = 0, found = 0;
            var m = Regex.Match(_a.Names[i], @"^([^/]+)/\s*([+-]?\d+)-\s*([+-]?\d+)");
            bool isRange = m.Success;
            if (isRange)
                tmpStart = int.Parse(m.Groups[2].Value, CultureInfo.InvariantCulture);
            for (int k = 1; k < _firstRes; k++)
            {
                // a column behind the end of a (converted, unaligned) row is read from the unused capacity of the program's vector: zero, i.e. not a gap
                if (Cell(i, k) == '-')
                    pregaps++;
            }

            // len = lastRes: the loop runs over firstRes .. firstRes + lastRes - 1 and stops at the end of the row
            for (int j = _firstRes; j < _firstRes + _lastRes; j++)
            {
                var c = Cell(i, j);
                if (c is null)
                    break;
                if (c == '-')
                    ngaps++;
                else
                    found = j;
                if (found != 0 && iStart == 0)
                    iStart = found;
            }

            int formula = 1;
            if (_o.SequenceNumberRange)
            {
                formula = iStart - pregaps;
                if (isRange)
                    formula += tmpStart == 1 ? 0 : tmpStart - 1;
            }

            if (pregaps == _firstRes - 1)
                formula = tmpStart;
            formula = formula <= 0 ? 1 : formula;
            if (pregaps == 0 && tmpStart == 0)
                formula = _firstRes;
            return (formula, formula + _lastRes - ngaps - 1);
        }

        // the name part is cut to maxName characters first when the format cuts the names (CLUSTAL), then to the text before the first '/'
        private string RangeName(int ii, int maxName = int.MaxValue)
        {
            var (start, end) = FindRangeValues(ii);
            string name = _a.Names[Index(ii)];
            if (name.Length > maxName)
                name = name[..maxName];
            return $"{NameOnly(name)}/{start.ToString(CultureInfo.InvariantCulture)}-{end.ToString(CultureInfo.InvariantCulture)}";
        }

        private static string PadRight(string s, int width) => s.Length >= width ? s : s + new string(' ', width - s.Length);

        private static string PadLeft(string s, int width) => s.Length >= width ? s : new string(' ', width - s.Length) + s;

        private static int LineWidth(int namesWidth)
        {
            int w = PageWidth - namesWidth;
            w -= w % 10;
            return w > LineLength || w <= 0 ? LineLength : w;
        }

        // the residues of row i from firstRes to lastRes: a gap is '-', the row ends where its data end
        private string RowSlice(int i, char gap = '-')
        {
            var sb = new StringBuilder();
            for (int j = _firstRes; j <= _lastRes; j++)
            {
                var c = Cell(i, j);
                if (c is null)
                    break;
                sb.Append(Letter(c, gap));
            }

            return sb.ToString();
        }

        public string Fasta()
        {
            var sb = new StringBuilder();
            int lineLength = LineWidth(_maxNames);
            for (int ii = 0; ii < _n; ii++)
            {
                int i = Index(ii);
                sb.Append('>').Append(NameOnly(_a.Names[i]));
                if (_o.SequenceNumberRange)
                {
                    var (s, e) = FindRangeValues(ii);
                    sb.Append('/').Append(s.ToString(CultureInfo.InvariantCulture)).Append('-').Append(e.ToString(CultureInfo.InvariantCulture));
                }

                sb.Append('\n');
                AppendWrapped(sb, RowSlice(i), lineLength);
            }

            return sb.ToString();
        }

        private static void AppendWrapped(StringBuilder sb, string text, int lineLength)
        {
            for (int j = 1; j <= text.Length; j++)
            {
                sb.Append(text[j - 1]);
                if (j % lineLength == 0 || j == text.Length)
                    sb.Append('\n');
            }
        }

        public string Pir()
        {
            var sb = new StringBuilder();
            int lineLength = LineWidth(_maxNames);
            string prefix = Dna ? ">DL;" : ">P1;";
            for (int ii = 0; ii < _n; ii++)
            {
                int i = Index(ii);
                sb.Append(prefix);
                sb.Append(_o.SequenceNumberRange ? RangeName(ii) : _a.Names[i]).Append('\n').Append('\n');   // the title of a FASTA entry is empty
                AppendWrapped(sb, RowSlice(i), lineLength);
                sb.Append("*\n");
            }

            return sb.ToString();
        }

        public string Gde()
        {
            var sb = new StringBuilder();
            int lineLength = LineWidth(_maxNames);
            string prefix = Dna ? "#" : "%";
            for (int ii = 0; ii < _n; ii++)
            {
                int i = Index(ii);
                sb.Append(prefix).Append(_o.SequenceNumberRange ? RangeName(ii) : _a.Names[i]).Append('\n');
                string slice = RowSlice(i);
                if (_o.Lowercase)
                    slice = slice.ToLowerInvariant();
                AppendWrapped(sb, slice, lineLength);
            }

            AppendStructureGde(sb, lineLength);
            return sb.ToString();
        }

        // AlignmentOutput::printSecStructMask: the structure line as it is written (the first helixEndMinus / strandEndMinus residues of a run and the residues before its end are lower case)
        private char[] PrintedStructure(ClustalWAlignmentStructure st, int prfLength)
        {
            var mask = st.SecondaryStructure;
            char M(int idx) => idx >= 0 && idx < mask.Length ? mask[idx] : '\0';
            var ss = new char[prfLength + 10];
            for (int k = 0; k < prfLength; k++)
                ss[k] = M(k);
            int helix = _st.HelixEndInside, strand = _st.StrandEndInside;
            int i = 0;
            while (i < prfLength)
            {
                if (char.ToLowerInvariant(M(i)) == 'a' || M(i) == '$')
                {
                    int j;
                    for (j = 0; j < helix; j++)
                    {
                        if (i + j >= prfLength || (char.ToLowerInvariant(M(i + j)) != 'a' && M(i + j) != '$'))
                            break;
                        ss[i + j] = 'a';
                    }

                    i += j;
                    while (char.ToLowerInvariant(M(i)) == 'a' || M(i) == '$')
                    {
                        if (i >= prfLength)
                            break;
                        if (M(i) == '$')
                        {
                            ss[i] = 'A';
                            i++;
                            break;
                        }

                        ss[i] = M(i);
                        i++;
                    }

                    for (j = 0; j < helix; j++)
                    {
                        if (i - j - 1 >= 0 && (char.ToLowerInvariant(M(i - j - 1)) == 'a' || M(i - j - 1) == '$'))
                            ss[i - j - 1] = 'a';
                    }
                }
                else if (char.ToLowerInvariant(M(i)) == 'b' || M(i) == '%')
                {
                    int j;
                    for (j = 0; j < strand; j++)
                    {
                        if (i + j >= prfLength || (char.ToLowerInvariant(M(i + j)) != 'b' && M(i + j) != '%'))
                            break;
                        ss[i + j] = 'b';
                    }

                    i += j;
                    while (char.ToLowerInvariant(M(i)) == 'b' || M(i) == '%')
                    {
                        if (i >= prfLength)
                            break;
                        if (M(i) == '%')
                        {
                            ss[i] = 'B';
                            i++;
                            break;
                        }

                        ss[i] = M(i);
                        i++;
                    }

                    for (j = 0; j < strand; j++)
                    {
                        if (i - j - 1 >= 0 && (char.ToLowerInvariant(M(i - j - 1)) == 'b' || M(i - j - 1) == '%'))
                            ss[i - j - 1] = 'b';
                    }
                }
                else
                {
                    i++;
                }
            }

            return ss;
        }

        private int ProfileFirstRowLength(int profile) => _rows[profile == 1 ? 0 : _a.Profile1Count].Length;

        private static bool SecStructUsed(ClustalWAlignmentStructure? st) => st is { HasSecondaryStructure: true, Used: true };

        private void AppendStructureGde(StringBuilder sb, int lineLength)
        {
            int length = _lastRes - _firstRes + 1;
            var buf = new char[Math.Max(_length, 1) + 1];
            Array.Fill(buf, '0');
            // the buffer holds the characters of the last sequence written (it is reused for the structure lines)
            if (_n > 0)
            {
                string last = RowSlice(Index(_n - 1));
                if (_o.Lowercase)
                    last = last.ToLowerInvariant();
                for (int k = 0; k < last.Length && k < buf.Length; k++)
                    buf[k] = last[k];
            }

            void Block(bool useLastRes)
            {
                int count = useLastRes ? _lastRes : length;
                for (int i = 1; i <= count; i++)
                {
                    sb.Append(buf[i - 1]);
                    if (i % lineLength == 0 || i == count)
                        sb.Append('\n');
                }
            }

            static char MaskChar(char c) => c;
            int outStruct = (int)_o.Structure;
            var s1 = _a.Structure1;
            var s2 = _a.Structure2;
            if (outStruct is 0 or 2)
            {
                if (SecStructUsed(s1))
                {
                    sb.Append("\"SS_").Append(PadRight(s1!.Name, _maxNames)).Append('\n');
                    var ss = PrintedStructure(s1, ProfileFirstRowLength(1));
                    for (int i = _firstRes; i <= _lastRes; i++)
                        buf[i - _firstRes] = i - 1 < ss.Length ? MaskChar(ss[i - 1]) : '\0';
                    Block(useLastRes: true);
                }

                if (SecStructUsed(s2))
                {
                    sb.Append("\"SS_").Append(PadRight(s2!.Name, _maxNames)).Append('\n');
                    var ss = PrintedStructure(s2, ProfileFirstRowLength(2));
                    for (int i = _firstRes; i <= _lastRes; i++)
                        buf[i - _firstRes] = i - 1 < ss.Length ? MaskChar(ss[i - 1]) : '\0';
                    Block(useLastRes: false);
                }
            }

            if (outStruct is 1 or 2)
            {
                if (s1 is { Used: true })
                {
                    sb.Append("\"GM_").Append(PadRight(s1.Name, _maxNames)).Append('\n');
                    for (int i = _firstRes; i <= _lastRes; i++)
                        buf[i - _firstRes] = i - 1 < s1.GapPenaltyMask.Length ? MaskChar(s1.GapPenaltyMask[i - 1]) : '\0';
                    Block(useLastRes: false);
                }

                if (s2 is { Used: true })
                {
                    sb.Append("\"GM_").Append(PadRight(s2.Name, _maxNames)).Append('\n');
                    // the loop of the program stops one entry short (i < length): the last entry keeps what the buffer held
                    for (int i = _firstRes; i < length; i++)
                        buf[i - _firstRes] = i - 1 < s2.GapPenaltyMask.Length ? MaskChar(s2.GapPenaltyMask[i - 1]) : '\0';
                    Block(useLastRes: false);
                }
            }
        }

        public string Phylip()
        {
            var sb = new StringBuilder();
            int length = _lastRes - _firstRes + 1;
            int chunks = length / GcgLineLength;
            if (length % GcgLineLength != 0)
                chunks++;
            sb.Append(PadLeft(_n.ToString(CultureInfo.InvariantCulture), 6)).Append(' ').Append(PadLeft(length.ToString(CultureInfo.InvariantCulture), 6));
            for (int block = 1; block <= chunks; block++)
            {
                int pos1 = ((block - 1) * GcgLineLength) + 1;
                int pos2 = length < pos1 + GcgLineLength - 1 ? length : pos1 + GcgLineLength - 1;
                for (int ii = 0; ii < _n; ii++)
                {
                    int i = Index(ii);
                    if (block == 1)
                    {
                        sb.Append('\n');
                        if (!_o.SequenceNumberRange)
                        {
                            string name = _a.Names[i];
                            sb.Append(PadRight(name.Length > 10 ? name[..10] : name, 10)).Append(' ');
                        }
                        else
                        {
                            sb.Append(PadRight(RangeName(ii), _maxNames + 15));
                        }
                    }
                    else
                    {
                        sb.Append("\n           ");
                    }

                    for (int j = pos1; j <= pos2; j++)
                    {
                        var c = Cell(i, j + _firstRes - 1);
                        if (c is null)
                            break;
                        sb.Append(Letter(c));
                        if (j % 10 == 0)
                            sb.Append(' ');
                    }
                }

                sb.Append('\n');
            }

            return sb.ToString();
        }

        public string Nexus()
        {
            var sb = new StringBuilder();
            int length = _lastRes - _firstRes + 1;
            int chunks = length / GcgLineLength;
            if (length % GcgLineLength != 0)
                chunks++;
            sb.Append("#NEXUS\n");
            sb.Append("BEGIN DATA;\n");
            sb.Append("dimensions ntax=").Append(_n.ToString(CultureInfo.InvariantCulture)).Append(" nchar=").Append(length.ToString(CultureInfo.InvariantCulture)).Append(";\n");
            sb.Append("format missing=?\n");
            sb.Append("interleave datatype=").Append(Dna ? "DNA " : "PROTEIN ").Append("gap= -;\n");
            sb.Append("\nmatrix");
            for (int block = 1; block <= chunks; block++)
            {
                int pos1 = ((block - 1) * GcgLineLength) + 1;
                int pos2 = length < pos1 + GcgLineLength - 1 ? length : pos1 + GcgLineLength - 1;
                for (int ii = 0; ii < _n; ii++)
                {
                    int i = Index(ii);
                    sb.Append('\n');
                    if (!_o.SequenceNumberRange)
                        sb.Append(PadRight(_a.Names[i], _maxNames + 1)).Append(' ');
                    else
                        sb.Append(PadRight(RangeName(ii), _maxNames + 15));
                    for (int j = pos1; j <= pos2; j++)
                    {
                        var c = Cell(i, j + _firstRes - 1);
                        if (c is null)
                            break;
                        sb.Append(Letter(c));
                    }
                }

                sb.Append('\n');
            }

            sb.Append(";\nend;\n");
            return sb.ToString();
        }

        public string Gcg()
        {
            var sb = new StringBuilder();
            int length = _lastRes - _firstRes + 1;
            var checks = new int[_n];
            for (int i = 0; i < _n; i++)
            {
                // sequence[1 .. lastRes]: the residues of the range, then '.' up to lastRes (gaps are '.')
                var buf = new List<char>();
                for (int j = _firstRes; j <= _lastRes; j++)
                {
                    var c = Cell(i, j);
                    if (c is null)
                        break;
                    buf.Add(Letter(c, '.'));
                }

                while (buf.Count < _lastRes)
                    buf.Add('.');
                long check = 0;
                for (int k = 0; k < _lastRes; k++)
                    check += ((k % 57) + 1) * char.ToUpperInvariant(buf[k]);
                checks[i] = (int)(check % 10000);
            }

            long grand = 0;
            for (int ii = 0; ii < _n; ii++)
                grand += checks[Index(ii)];
            grand %= 10000;

            sb.Append("PileUp\n\n");
            sb.Append("\n\n   MSF:").Append(PadLeft(_lastRes.ToString(CultureInfo.InvariantCulture), 5)).Append("  Type: ").Append(Dna ? 'N' : 'P');
            sb.Append("    Check:").Append(PadLeft(grand.ToString(CultureInfo.InvariantCulture), 6)).Append("   .. \n\n");
            for (int ii = 0; ii < _n; ii++)
            {
                int i = Index(ii);
                float weight = (_a.SequenceWeights[i] * 100) / (float)1000;
                sb.Append(" Name: ").Append(_a.Names[i]).Append(" oo  Len:").Append(PadLeft(length.ToString(CultureInfo.InvariantCulture), 5))
                    .Append("  Check:").Append(PadLeft(checks[i].ToString(CultureInfo.InvariantCulture), 6))
                    .Append("  Weight:  ").Append(SequenceAligner.ClustalWFixed(weight, 1)).Append('\n');
            }

            sb.Append("\n//\n");
            int chunks = length / GcgLineLength;
            if (length % GcgLineLength != 0)
                chunks++;
            for (int block = 1; block <= chunks; block++)
            {
                sb.Append("\n\n");
                int pos1 = ((block - 1) * GcgLineLength) + 1;
                int pos2 = length < pos1 + GcgLineLength - 1 ? length : pos1 + GcgLineLength - 1;
                for (int ii = 0; ii < _n; ii++)
                {
                    int i = Index(ii);
                    sb.Append('\n');
                    if (!_o.SequenceNumberRange)
                        sb.Append(PadRight(_a.Names[i], _maxNames + 5)).Append(' ');
                    else
                        sb.Append(PadRight(RangeName(ii), _maxNames + 15));
                    for (int j = pos1; j <= pos2; j++)
                    {
                        sb.Append(Letter(Cell(i, j + _firstRes - 1), '.'));
                        if (j % 10 == 0)
                            sb.Append(' ');
                    }
                }
            }

            sb.Append("\n\n");
            return sb.ToString();
        }

        public string Clustal()
        {
            var sb = new StringBuilder();
            int namesWidth = Math.Max(Math.Min(MaxNamesToDisplay, _maxNames), MinNamesToDisplay);
            var seqNo = new int[_n];
            var printSeqNo = new bool[_n];
            for (int i = 0; i < _n; i++)
            {
                // for j < firstRes: `if ((val >= 0) || (val <= maxAA))` is always true: every column before the range counts, gaps included
                for (int j = 1; j < _firstRes; j++)
                    seqNo[i]++;
            }

            sb.Append("CLUSTAL 2.1 multiple sequence alignment\n\n");
            int lineLength = LineWidth(namesWidth);
            int length = _lastRes - _firstRes + 1;
            int chunks = length / lineLength;
            if (length % lineLength != 0)
                chunks++;
            for (int lv1 = 1; lv1 <= chunks; lv1++)
            {
                int pos = ((lv1 - 1) * lineLength) + 1;
                int ptr = length < pos + lineLength - 1 ? length : pos + lineLength - 1;
                sb.Append('\n');
                AppendStructureClustal(sb, pos, ptr, namesWidth);
                for (int ii = 0; ii < _n; ii++)
                {
                    int i = Index(ii);
                    printSeqNo[i] = false;
                    var line = new StringBuilder();
                    for (int j = pos; j <= ptr; j++)
                    {
                        var c = Cell(i, j + _firstRes - 1);
                        if (c is null)
                        {
                            while (line.Length < ptr - pos + 1)
                                line.Append('-');
                            break;
                        }

                        if (c == '-')
                        {
                            line.Append('-');
                        }
                        else
                        {
                            line.Append(c.Value);
                            seqNo[i]++;
                            printSeqNo[i] = true;
                        }
                    }

                    string seqName = _a.Names[i];
                    if (seqName.Length > namesWidth)
                        seqName = seqName[..namesWidth];
                    if (!_o.SequenceNumberRange)
                        sb.Append(PadRight(seqName, namesWidth + 6));
                    else
                        sb.Append(PadRight(RangeName(ii, namesWidth), namesWidth + ClusSequenceOffset));
                    sb.Append(line);
                    if (_o.SequenceNumbers && printSeqNo[i])
                        sb.Append(' ').Append(seqNo[i].ToString(CultureInfo.InvariantCulture));
                    sb.Append('\n');
                }

                // conservation line: all sequences identical to the first one (*), the same strong (:) or weak (.) group (protein)
                var cons = new StringBuilder();
                for (int p = pos; p <= ptr; p++)
                    cons.Append(Conservation(p + _firstRes - 1));
                sb.Append(new string(' ', namesWidth + 6));
                if (_o.SequenceNumberRange)
                    sb.Append("          ");
                sb.Append(cons).Append('\n');
            }

            return sb.ToString();
        }

        private char Conservation(int col)
        {
            char first = col <= _rows[0].Length ? _rows[0][col - 1] : '-';
            if (first == '-')
                return ' ';
            int ident = 0;
            var strong = new int[StrongGroups.Length];
            var weak = new int[WeakGroups.Length];
            for (int j = 0; j < _n; j++)
            {
                if (col > _rows[j].Length || col > _rows[0].Length)
                    continue;
                char c = _rows[j][col - 1];
                if (c == first)
                    ident++;
                if (c == '-')
                    continue;
                for (int k = 0; k < StrongGroups.Length; k++)
                {
                    if (StrongGroups[k].Contains(c, StringComparison.Ordinal))
                        strong[k]++;
                }

                for (int k = 0; k < WeakGroups.Length; k++)
                {
                    if (WeakGroups[k].Contains(c, StringComparison.Ordinal))
                        weak[k]++;
                }
            }

            if (ident == _n)
                return '*';
            if (!Dna)
            {
                if (strong.Any(x => x == _n))
                    return ':';
                if (weak.Any(x => x == _n))
                    return '.';
            }

            return ' ';
        }

        private void AppendStructureClustal(StringBuilder sb, int pos, int ptr, int namesWidth)
        {
            int outStruct = (int)_o.Structure;
            AppendProfileStructureClustal(sb, pos, ptr, namesWidth, _a.Structure1, 1, outStruct);
            AppendProfileStructureClustal(sb, pos, ptr, namesWidth, _a.Structure2, 2, outStruct);
        }

        private void AppendProfileStructureClustal(StringBuilder sb, int pos, int ptr, int namesWidth, ClustalWAlignmentStructure? st, int profile, int outStruct)
        {
            if (st == null)
                return;
            string name = st.Name;
            // profile 1 compares with MAXNAMESTODISPLAY, profile 2 with the names width (both cut to the names width)
            if (name.Length > (profile == 1 ? MaxNamesToDisplay : namesWidth))
                name = name[..namesWidth];
            string Line(char[] source)
            {
                var sbl = new StringBuilder();
                for (int i = pos; i <= ptr; i++)
                {
                    int idx = i + _firstRes - 2;
                    if (source.Length <= idx)
                        break;
                    sbl.Append(source[idx]);
                }

                string text = sbl.ToString();
                int nul = text.IndexOf('\0', StringComparison.Ordinal);
                return nul >= 0 ? text[..nul] : text;
            }

            if (outStruct is 0 or 2 && st is { HasSecondaryStructure: true, Used: true })
            {
                var ss = PrintedStructure(st, ProfileFirstRowLength(profile));
                // the copy of the program has length lengthSeq + 10 (zero filled)
                sb.Append("!SS_").Append(PadRight(name, _o.SequenceNumberRange ? namesWidth + ClusSecStructOffset : namesWidth)).Append("  ").Append(Line(ss)).Append('\n');
            }

            if (outStruct is 1 or 2 && st is { Used: true })
                sb.Append("!GM_").Append(PadRight(name, namesWidth)).Append("  ").Append(Line(st.GapPenaltyMask.ToCharArray())).Append('\n');
        }
    }
}
