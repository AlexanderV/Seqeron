using System.Globalization;
using System.Text;
using System.Text.Json.Nodes;
using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Alignment;

/// <summary>Which row of an <see cref="AlignmentResult"/> is the query and what happens to the free end gaps in the bridge to the coordinate formats (<see cref="AlignmentResultCoordinates"/>).</summary>
public sealed record AlignmentBridgeOptions
{
    /// <summary>
    /// False (default): sequence 1 is the target (reference) and sequence 2 the query, as Biopython's <c>aligner.align(target, query)</c> and the <see cref="BioCoordAlignment"/> order; a residue of row 2 against a gap of row 1
    /// is an insertion (<c>I</c>), a residue of row 1 against a gap of row 2 a deletion (<c>D</c>). True: sequence 1 is the query (the result of <c>SemiGlobalAlign(query, reference)</c>: sequence 2 is the reference).
    /// </summary>
    public bool FirstIsQuery { get; init; }

    /// <summary>
    /// False (default): the coordinates are those of Biopython's <c>PairwiseAligner</c>: the free end gaps of a semi-global (ends-free) result stay in the alignment as leading / trailing gap blocks (<c>2D4M2D</c>).
    /// True: the leading and trailing columns that have a gap in either row are dropped (Biopython <c>alignment[:, first:last]</c>), so that the alignment starts and ends with an aligned pair as a SAM line, a PSL or a BED record does (<c>4M</c> at POS 3).
    /// </summary>
    public bool TrimEndGaps { get; init; }
}

/// <summary>
/// The bridge from the rows of an <see cref="AlignmentResult"/> to the coordinate world (R14-9): the block coordinates of Biopython's <c>Bio.Align.Alignment</c> (<see cref="BioCoordAlignment"/>, written by
/// <see cref="AlignmentCoordinateFormats"/> as BED / PSL / chain / MAF / SAM / Exonerate), the standard and the extended CIGAR, the SAM <c>MD</c> tag and the minimap2 <c>cs</c> string. The coordinates are derived from the
/// rows alone (the residues before and after the rows come from <see cref="EmbossAlignedSequence.Flanks"/>), so they do not depend on the <c>End</c> convention of the result.
/// </summary>
public static class AlignmentResultCoordinates
{
    public const char Gap = '-';

    /// <summary>The column walk of a result: runs of aligned pairs (<c>M</c>), insertions (<c>I</c>, query only) and deletions (<c>D</c>, target only) with the positions they start at.</summary>
    private sealed class Walk
    {
        public required string TargetRow { get; init; }
        public required string QueryRow { get; init; }
        public required long TargetStart { get; init; }
        public required long QueryStart { get; init; }
        public required long TargetLength { get; init; }
        public required long QueryLength { get; init; }
        public required List<(char Op, int Length)> Runs { get; init; }
        public required long TrimmedTargetStart { get; init; }
        public required long TrimmedQueryStart { get; init; }
        public required long TrimmedTargetEnd { get; init; }
        public required long TrimmedQueryEnd { get; init; }
    }

    // the residues of the sequence that precede the aligned part: the start of the row, clamped to what the alignment leaves of the sequence
    private static (int Offset, int Offend) LeadingFlank(string row, int start, int length)
    {
        int consumed = 0;
        foreach (char c in row)
        {
            if (c is not ('-' or '.' or '~' or '?' or ' '))
                consumed++;
        }

        if (consumed >= length)
            return (0, 0);
        int offset = Math.Clamp(start, 0, length - consumed);
        return (offset, length - offset - consumed);
    }

    private static Walk Analyse(AlignmentResult r, long length1, long length2, AlignmentBridgeOptions? options)
    {
        ArgumentNullException.ThrowIfNull(r);
        if (r.AlignedSequence1.Length != r.AlignedSequence2.Length)
            throw new ArgumentException("The two rows of an alignment have the same length.", nameof(r));
        if (length1 < 0 || length2 < 0)
            throw new ArgumentOutOfRangeException(nameof(length1), "A sequence length is not negative.");
        options ??= new AlignmentBridgeOptions();
        var (off1, _) = LeadingFlank(r.AlignedSequence1, r.StartPosition1, (int)length1);
        var (off2, _) = LeadingFlank(r.AlignedSequence2, r.StartPosition2, (int)length2);
        bool swap = options.FirstIsQuery;
        string tRow = swap ? r.AlignedSequence2 : r.AlignedSequence1;
        string qRow = swap ? r.AlignedSequence1 : r.AlignedSequence2;
        long tOff = swap ? off2 : off1, qOff = swap ? off1 : off2;
        long tLen = swap ? length2 : length1, qLen = swap ? length1 : length2;

        var runs = new List<(char, int)>();
        for (int i = 0; i < tRow.Length; i++)
        {
            bool tGap = tRow[i] == Gap, qGap = qRow[i] == Gap;
            if (tGap && qGap)
                throw new ArgumentException("A column of the alignment has a gap in both rows.", nameof(r));
            char op = 'M';
            if (tGap)
                op = 'I';
            else if (qGap)
                op = 'D';

            if (runs.Count > 0 && runs[^1].Item1 == op)
                runs[^1] = (op, runs[^1].Item2 + 1);
            else
                runs.Add((op, 1));
        }

        long tStart = tOff, qStart = qOff;
        long tEnd = tOff + tRow.Count(c => c != Gap), qEnd = qOff + qRow.Count(c => c != Gap);
        if (options.TrimEndGaps)
        {
            while (runs.Count > 0 && runs[0].Item1 != 'M')
            {
                if (runs[0].Item1 == 'D')
                    tStart += runs[0].Item2;
                else
                    qStart += runs[0].Item2;
                runs.RemoveAt(0);
            }

            while (runs.Count > 0 && runs[^1].Item1 != 'M')
            {
                if (runs[^1].Item1 == 'D')
                    tEnd -= runs[^1].Item2;
                else
                    qEnd -= runs[^1].Item2;
                runs.RemoveAt(runs.Count - 1);
            }
        }

        return new Walk
        {
            TargetRow = tRow,
            QueryRow = qRow,
            TargetStart = tOff,
            QueryStart = qOff,
            TargetLength = tLen,
            QueryLength = qLen,
            Runs = runs,
            TrimmedTargetStart = tStart,
            TrimmedQueryStart = qStart,
            TrimmedTargetEnd = tEnd,
            TrimmedQueryEnd = qEnd,
        };
    }

    private static BioCoordSequence Record(string id, long length, string? full, string row, long offset)
    {
        if (full is not null)
        {
            string aligned = row.Replace("-", "", StringComparison.Ordinal);
            if (full.Length != length || offset + aligned.Length > full.Length || !string.Equals(full.Substring((int)offset, aligned.Length), aligned, StringComparison.OrdinalIgnoreCase))
                throw new ArgumentException($"The sequence '{id}' is not the sequence the alignment rows were made from.");
            return new BioCoordSequence { Id = id, Length = length, Parts = AlignmentCoordinateFormats.SeqParts([new KeyValuePair<long, string>(0, full)], length) };
        }

        string residues = row.Replace("-", "", StringComparison.Ordinal);
        var parts = residues.Length == 0
            ? AlignmentCoordinateFormats.UndefinedParts(length)
            : AlignmentCoordinateFormats.SeqParts([new KeyValuePair<long, string>(offset, residues)], length);
        return new BioCoordSequence { Id = id, Length = length, Parts = parts };
    }

    private static BioCoordAlignment Build(AlignmentResult r, string name1, long length1, string? sequence1, string name2, long length2, string? sequence2, AlignmentBridgeOptions? options)
    {
        ArgumentNullException.ThrowIfNull(name1);
        ArgumentNullException.ThrowIfNull(name2);
        var w = Analyse(r, length1, length2, options);
        if (w.Runs.Count == 0)
            throw new InvalidOperationException("The alignment has no aligned column: there are no block coordinates to describe.");
        bool swap = options?.FirstIsQuery ?? false;
        var target = swap ? Record(name2, length2, sequence2, w.TargetRow, w.TargetStart) : Record(name1, length1, sequence1, w.TargetRow, w.TargetStart);
        var query = swap ? Record(name1, length1, sequence1, w.QueryRow, w.QueryStart) : Record(name2, length2, sequence2, w.QueryRow, w.QueryStart);
        var t = new List<long> { w.TrimmedTargetStart };
        var q = new List<long> { w.TrimmedQueryStart };
        foreach (var (op, n) in w.Runs)
        {
            t.Add(t[^1] + (op == 'I' ? 0 : n));
            q.Add(q[^1] + (op == 'D' ? 0 : n));
        }

        return new BioCoordAlignment
        {
            Sequences = [target, query],
            Coordinates = [[.. t], [.. q]],
            Attributes = new JsonObject { ["score"] = JsonValue.Create(r.ExactScore) },
        };
    }

    /// <summary>
    /// The alignment as Biopython's <c>Alignment</c> (target first, query second, <c>score</c>): the block coordinates are the runs of aligned pairs and gaps of the rows. The residues of the sequences are known only where the rows carry them
    /// (the aligned part: a SAM line of a local result then has <c>*</c> as its SEQ and no <c>MD</c>); use the overload with the sequences for the complete records.
    /// </summary>
    /// <param name="result">The result of any pairwise aligner of the library.</param>
    /// <param name="name1">The identifier of sequence 1.</param>
    /// <param name="length1">The length of sequence 1.</param>
    /// <param name="name2">The identifier of sequence 2.</param>
    /// <param name="length2">The length of sequence 2.</param>
    /// <param name="options">The roles of the sequences and the end gaps (default: sequence 1 is the target, end gaps stay).</param>
    /// <exception cref="InvalidOperationException">The result has no column (an empty local result).</exception>
    public static BioCoordAlignment ToBioCoordAlignment(this AlignmentResult result, string name1, int length1, string name2, int length2, AlignmentBridgeOptions? options = null)
        => Build(result, name1, length1, null, name2, length2, null, options);

    /// <summary>The alignment as Biopython's <c>Alignment</c> with the complete sequences (the strings the aligner was called with; the aligners compare letters case-insensitively and write them upper case).</summary>
    /// <param name="result">The result of any pairwise aligner of the library.</param>
    /// <param name="name1">The identifier of sequence 1.</param>
    /// <param name="sequence1">Sequence 1.</param>
    /// <param name="name2">The identifier of sequence 2.</param>
    /// <param name="sequence2">Sequence 2.</param>
    /// <param name="options">The roles of the sequences and the end gaps.</param>
    /// <exception cref="InvalidOperationException">The result has no column.</exception>
    public static BioCoordAlignment ToBioCoordAlignment(this AlignmentResult result, string name1, string sequence1, string name2, string sequence2, AlignmentBridgeOptions? options = null)
    {
        ArgumentNullException.ThrowIfNull(sequence1);
        ArgumentNullException.ThrowIfNull(sequence2);
        return Build(result, name1, sequence1.Length, sequence1, name2, sequence2.Length, sequence2, options);
    }

    private static string RunText(IEnumerable<(char Op, int Length)> runs)
    {
        var sb = new StringBuilder();
        foreach (var (op, n) in runs)
            sb.Append(n.ToString(CultureInfo.InvariantCulture)).Append(op);
        return sb.ToString();
    }

    // M runs split into = / X runs (a letter pair is equal when the letters are equal ignoring case: the aligners compare case-insensitively)
    private static List<(char Op, int Length)> Operations(Walk w, bool extended)
    {
        if (!extended)
            return w.Runs.ToList();
        var list = new List<(char, int)>();
        var (tRes, qRes) = AlignedResidues(w);
        int tIdx = 0, qIdx = 0;
        foreach (var (op, n) in w.Runs)
        {
            if (op == 'I')
            {
                Add(list, 'I', n);
                qIdx += n;
            }
            else if (op == 'D')
            {
                Add(list, 'D', n);
                tIdx += n;
            }
            else
            {
                for (int k = 0; k < n; k++)
                    Add(list, char.ToUpperInvariant(tRes[tIdx + k]) == char.ToUpperInvariant(qRes[qIdx + k]) ? '=' : 'X', 1);
                tIdx += n;
                qIdx += n;
            }
        }

        return list;

        static void Add(List<(char, int)> l, char op, int n)
        {
            if (l.Count > 0 && l[^1].Item1 == op)
                l[^1] = (op, l[^1].Item2 + n);
            else
                l.Add((op, n));
        }
    }

    /// <summary>
    /// The CIGAR of the aligned columns: <c>M</c> for an aligned pair (standard) or <c>=</c> / <c>X</c> for an identical / different pair (<paramref name="extended"/>, the SAM extended CIGAR of minimap2 <c>--eqx</c> and parasail),
    /// <c>I</c> for a residue of the query against a gap, <c>D</c> for a residue of the target against a gap. An empty alignment gives <c>*</c>.
    /// </summary>
    /// <param name="result">The result of any pairwise aligner.</param>
    /// <param name="extended">Write <c>=</c> / <c>X</c> instead of <c>M</c>.</param>
    /// <param name="options">The roles of the sequences (default: row 1 is the target, row 2 the query) and the end gaps.</param>
    public static string ToCigar(this AlignmentResult result, bool extended = false, AlignmentBridgeOptions? options = null)
    {
        var w = Analyse(result, int.MaxValue, int.MaxValue, options);
        return w.Runs.Count == 0 ? "*" : RunText(Operations(w, extended));
    }

    /// <summary>
    /// The CIGAR of the SAM line <c>Bio.Align</c> writes for the alignment: the aligned columns as <see cref="ToCigar(AlignmentResult, bool, AlignmentBridgeOptions?)"/> and the query residues outside the alignment as soft clips
    /// (<c>S</c>) before and after (a local result: the query residues around the aligned part).
    /// </summary>
    /// <param name="result">The result of any pairwise aligner.</param>
    /// <param name="length1">The length of sequence 1.</param>
    /// <param name="length2">The length of sequence 2.</param>
    /// <param name="extended">Write <c>=</c> / <c>X</c> instead of <c>M</c>.</param>
    /// <param name="options">The roles of the sequences and the end gaps.</param>
    public static string ToCigar(this AlignmentResult result, int length1, int length2, bool extended = false, AlignmentBridgeOptions? options = null)
    {
        var w = Analyse(result, length1, length2, options);
        if (w.Runs.Count == 0)
            return "*";
        var sb = new StringBuilder();
        if (w.TrimmedQueryStart > 0)
            sb.Append(w.TrimmedQueryStart.ToString(CultureInfo.InvariantCulture)).Append('S');
        sb.Append(RunText(Operations(w, extended)));
        if (w.TrimmedQueryEnd < w.QueryLength)
            sb.Append((w.QueryLength - w.TrimmedQueryEnd).ToString(CultureInfo.InvariantCulture)).Append('S');
        return sb.ToString();
    }

    /// <summary>
    /// The SAM <c>MD</c> tag value of the alignment (mismatching and deleted target residues; <c>alignment.format("sam", md=True)</c> of Biopython: the number before a leading deletion is omitted when it is zero).
    /// Letters are compared as they stand in the rows (case-sensitive, as Biopython).
    /// </summary>
    /// <param name="result">The result of any pairwise aligner.</param>
    /// <param name="options">The roles of the sequences (row 1 is the target by default) and the end gaps.</param>
    /// <exception cref="InvalidOperationException">The result has no column.</exception>
    public static string ToMd(this AlignmentResult result, AlignmentBridgeOptions? options = null)
    {
        var w = Analyse(result, int.MaxValue, int.MaxValue, options);
        if (w.Runs.Count == 0)
            throw new InvalidOperationException("The alignment has no aligned column.");
        // the MD tag depends on the residues of the columns only: describe them as two sequences that start at 0
        var (tRes, qRes) = AlignedResidues(w);
        var t = new List<long> { 0 };
        var q = new List<long> { 0 };
        foreach (var (op, n) in w.Runs)
        {
            t.Add(t[^1] + (op == 'I' ? 0 : n));
            q.Add(q[^1] + (op == 'D' ? 0 : n));
        }

        var target = new BioCoordSequence { Id = "target", Length = tRes.Length, Parts = new SortedDictionary<long, string> { [0] = tRes } };
        var query = new BioCoordSequence { Id = "query", Length = qRes.Length, Parts = new SortedDictionary<long, string> { [0] = qRes } };
        return AlignmentCoordinateFormats.MdTag(new BioCoordAlignment { Sequences = [target, query], Coordinates = [[.. t], [.. q]] });
    }

    /// <summary>
    /// The minimap2 <c>cs</c> difference string of the alignment (<c>:N</c> for N matching bases or, <paramref name="longForm"/>, <c>=ACGT</c>; <c>*ag</c> a substitution target-&gt;query; <c>+acg</c> an insertion in the query; <c>-acg</c> a deletion),
    /// by the <c>mm_gen_cs</c> port of <see cref="Minimap2Chaining"/>: nucleotides only, every other letter is <c>n</c>.
    /// </summary>
    /// <param name="result">The result of any pairwise aligner.</param>
    /// <param name="longForm">The long form (<c>--cs=long</c>).</param>
    /// <param name="options">The roles of the sequences and the end gaps.</param>
    public static string ToCs(this AlignmentResult result, bool longForm = false, AlignmentBridgeOptions? options = null)
    {
        var w = Analyse(result, int.MaxValue, int.MaxValue, options);
        if (w.Runs.Count == 0)
            return "";
        var (tRes, qRes) = AlignedResidues(w);
        return Minimap2Chaining.CsOfAlignment(tRes, qRes, w.Runs, longForm);
    }

    // the residues of the columns that are kept after the trim
    private static (string Target, string Query) AlignedResidues(Walk w)
    {
        string tAll = w.TargetRow.Replace("-", "", StringComparison.Ordinal), qAll = w.QueryRow.Replace("-", "", StringComparison.Ordinal);
        int t0 = (int)(w.TrimmedTargetStart - w.TargetStart), q0 = (int)(w.TrimmedQueryStart - w.QueryStart);
        return (tAll.Substring(t0, (int)(w.TrimmedTargetEnd - w.TrimmedTargetStart)), qAll.Substring(q0, (int)(w.TrimmedQueryEnd - w.TrimmedQueryStart)));
    }

    /// <summary>True when <paramref name="format"/> names one of the coordinate / CIGAR formats of <see cref="Format"/>.</summary>
    public static bool IsFormatName(string? format) => TryParseSpec(format, out _);

    private enum Kind { Bed, Psl, Chain, Maf, Sam, Exonerate, Cigar, Md, Cs }

    private readonly record struct Spec(Kind Kind, int BedN = 12, bool SamMd = false, bool ExonerateCigar = false, bool Extended = false, bool LongCs = false);

    private static bool TryParseSpec(string? format, out Spec spec)
    {
        spec = default;
        string f = (format ?? "").Trim().ToLowerInvariant();
        string option = "";
        int colon = f.IndexOf(':');
        if (colon >= 0)
        {
            option = f[(colon + 1)..];
            f = f[..colon];
        }

        if (f.StartsWith("bed", StringComparison.Ordinal) && option.Length == 0)
        {
            string digits = f[3..];
            int n = 12;
            if (digits.Length > 0 && (!int.TryParse(digits, NumberStyles.None, CultureInfo.InvariantCulture, out n) || n is < 3 or > 12))
                return false;
            spec = new Spec(Kind.Bed, n);
            return true;
        }

        switch (f)
        {
            case "psl" when option.Length == 0: spec = new Spec(Kind.Psl); return true;
            case "chain" when option.Length == 0: spec = new Spec(Kind.Chain); return true;
            case "maf" when option.Length == 0: spec = new Spec(Kind.Maf); return true;
            case "sam" when option is "" or "md": spec = new Spec(Kind.Sam, SamMd: option == "md"); return true;
            case "exonerate" when option is "" or "vulgar" or "cigar": spec = new Spec(Kind.Exonerate, ExonerateCigar: option == "cigar"); return true;
            case "cigar" when option is "" or "eqx" or "extended": spec = new Spec(Kind.Cigar, Extended: option != ""); return true;
            case "md" when option.Length == 0: spec = new Spec(Kind.Md); return true;
            case "cs" when option is "" or "long": spec = new Spec(Kind.Cs, LongCs: option == "long"); return true;
            default: return false;
        }
    }

    /// <summary>
    /// The alignment as text in one of the coordinate formats, byte for byte as Biopython 1.88 <c>alignment.format(...)</c> writes the coordinates of <see cref="ToBioCoordAlignment(AlignmentResult, string, string, string, string, AlignmentBridgeOptions?)"/>
    /// (no file header): <c>bed</c> (or <c>bed3</c> ... <c>bed12</c>), <c>psl</c>, <c>chain</c>, <c>maf</c>, <c>sam</c> (<c>sam:md</c> with the MD tag), <c>exonerate</c> (<c>exonerate:cigar</c>), or the strings
    /// <c>cigar</c> (<c>cigar:eqx</c> extended; with soft clips), <c>md</c>, <c>cs</c> (<c>cs:long</c>).
    /// </summary>
    /// <param name="result">The result of any pairwise aligner.</param>
    /// <param name="format">The format name.</param>
    /// <param name="name1">The identifier of sequence 1.</param>
    /// <param name="sequence1">Sequence 1.</param>
    /// <param name="name2">The identifier of sequence 2.</param>
    /// <param name="sequence2">Sequence 2.</param>
    /// <param name="options">The roles of the sequences and the end gaps.</param>
    /// <exception cref="ArgumentException">The format name is not one of the above.</exception>
    /// <exception cref="InvalidOperationException">The result has no column.</exception>
    public static string Format(this AlignmentResult result, string format, string name1, string sequence1, string name2, string sequence2, AlignmentBridgeOptions? options = null)
    {
        if (!TryParseSpec(format, out var spec))
            throw new ArgumentException($"Unknown coordinate format '{format}'. Names: bed, bed3 ... bed12, psl, chain, maf, sam, sam:md, exonerate, exonerate:cigar, cigar, cigar:eqx, md, cs, cs:long.", nameof(format));
        switch (spec.Kind)
        {
            case Kind.Cigar:
                return ToCigar(result, sequence1.Length, sequence2.Length, spec.Extended, options);
            case Kind.Md:
                return ToMd(result, options);
            case Kind.Cs:
                return ToCs(result, spec.LongCs, options);
        }

        var a = ToBioCoordAlignment(result, name1, sequence1, name2, sequence2, options);
        var write = new CoordinateWriteOptions { BedN = spec.BedN, Md = spec.SamMd, ExonerateFormat = spec.ExonerateCigar ? "cigar" : "vulgar" };
        var fmt = spec.Kind switch
        {
            Kind.Bed => CoordinateAlignmentFormat.Bed,
            Kind.Psl => CoordinateAlignmentFormat.Psl,
            Kind.Chain => CoordinateAlignmentFormat.Chain,
            Kind.Maf => CoordinateAlignmentFormat.Maf,
            Kind.Sam => CoordinateAlignmentFormat.Sam,
            _ => CoordinateAlignmentFormat.Exonerate,
        };
        return AlignmentCoordinateFormats.Format(a, fmt, write);
    }
}
