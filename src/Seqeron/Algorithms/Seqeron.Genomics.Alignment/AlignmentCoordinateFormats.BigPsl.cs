using System.Globalization;
using System.Text;
using System.Text.Json.Nodes;
using System.Text.RegularExpressions;
using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Alignment;

public static partial class AlignmentCoordinateFormats
{
    /// <summary>The autoSQL declaration of <c>Bio/Align/bigpsl.py</c> (<c>bigPsl</c>: BED12 plus the 13 PSL fields).</summary>
    public static BigBedDeclaration BigPslDeclaration { get; } = new(
        "bigPsl",
        "bigPsl pairwise alignment",
        [
            new BigBedField("string", "chrom", "Reference sequence chromosome or scaffold"),
            new BigBedField("uint", "chromStart", "Start position in chromosome"),
            new BigBedField("uint", "chromEnd", "End position in chromosome"),
            new BigBedField("string", "name", "Name or ID of item, ideally both human readable and unique"),
            new BigBedField("uint", "score", "Score (0-1000)"),
            new BigBedField("char[1]", "strand", "+ or - indicates whether the query aligns to the + or - strand on the reference"),
            new BigBedField("uint", "thickStart", "Start of where display should be thick (start codon)"),
            new BigBedField("uint", "thickEnd", "End of where display should be thick (stop codon)"),
            new BigBedField("uint", "reserved", "RGB value (use R,G,B string in input file)"),
            new BigBedField("int", "blockCount", "Number of blocks"),
            new BigBedField("int[blockCount]", "blockSizes", "Comma separated list of block sizes"),
            new BigBedField("int[blockCount]", "chromStarts", "Start positions relative to chromStart"),
            new BigBedField("uint", "oChromStart", "Start position in other chromosome"),
            new BigBedField("uint", "oChromEnd", "End position in other chromosome"),
            new BigBedField("char[1]", "oStrand", "+ or -, - means that psl was reversed into BED-compatible coordinates"),
            new BigBedField("uint", "oChromSize", "Size of other chromosome."),
            new BigBedField("int[blockCount]", "oChromStarts", "Start positions relative to oChromStart or from oChromStart+oChromSize depending on strand"),
            new BigBedField("lstring", "oSequence", "Sequence on other chrom (or edit list, or empty)"),
            new BigBedField("string", "oCDS", "CDS in NCBI format"),
            new BigBedField("uint", "chromSize", "Size of target chromosome"),
            new BigBedField("uint", "match", "Number of bases matched."),
            new BigBedField("uint", "misMatch", "Number of bases that don't match"),
            new BigBedField("uint", "repMatch", "Number of bases that match but are part of repeats"),
            new BigBedField("uint", "nCount", "Number of 'N' bases"),
            new BigBedField("uint", "seqType", "0=empty, 1=nucleotide, 2=amino_acid"),
        ]);

    /// <summary>
    /// Writes the pairwise alignments as <c>Bio.Align.write(alignments, handle, "bigpsl", ...)</c> does: a bigBed file of the BED12 columns plus the PSL columns of <c>bigPsl.as</c> (the match counts from the sequences, the other-chromosome fields of the query, the CDS in NCBI style,
    /// the sequence). The <c>extraIndex</c> argument of the Biopython writer is never used by it, so it is not an option here either.
    /// </summary>
    /// <param name="alignments">Pairwise alignments (target first); every alignment needs an <c>annotations</c> object attribute.</param>
    /// <param name="targets">The chromosomes (<c>alignments.targets</c>).</param>
    /// <param name="options">Compress / Compressor, Cds, Fa, Mask, Wildcard.</param>
    public static byte[] WriteBigPsl(IReadOnlyList<BioCoordAlignment> alignments, IReadOnlyList<BioCoordSequence> targets, BigBedWriteOptions? options = null)
    {
        ArgumentNullException.ThrowIfNull(alignments);
        ArgumentNullException.ThrowIfNull(targets);
        try
        {
            return WriteBigPslCore(alignments, targets, options ?? new BigBedWriteOptions());
        }
        catch (Exception e) when (e is IndexOutOfRangeException or ArgumentOutOfRangeException or KeyNotFoundException or InvalidCastException or NullReferenceException or FormatException or OverflowException or InvalidOperationException)
        {
            throw new ArgumentException(e.Message, e);
        }
    }

    private static char AsciiUpper(char c) => c is >= 'a' and <= 'z' ? (char)(c - 32) : c;

    private static char AsciiLower(char c) => c is >= 'A' and <= 'Z' ? (char)(c + 32) : c;

    /// <summary><c>bytes(seq[start:end])</c>: null when a residue of the slice is undefined (<c>UndefinedSequenceError</c>); the slice is clamped like a Python slice.</summary>
    private static string? PySliceBytes(BioCoordSequence seq, long start, long end)
    {
        if (seq.Parts is null)
            return null;
        if (seq.Parts.Count == 1 && seq.Parts.TryGetValue(0, out var whole) && whole.Length == seq.Length)
        {
            long s = Math.Min(start, whole.Length), e = Math.Min(end, whole.Length);
            return e <= s ? "" : whole.Substring((int)s, (int)(e - s));
        }

        return seq.Slice(start, end);
    }

    private static (long Matches, long MisMatches, long RepMatches, long NCount) CountBlock(string tSeq, string qSeq, string? mask, int? wildcard)
    {
        long matches = 0, mis = 0, rep = 0, n = 0;
        int count = Math.Min(tSeq.Length, qSeq.Length);
        for (int i = 0; i < count; i++)
        {
            char c1 = tSeq[i];
            char u1, u2;
            if (mask == "upper")
            {
                u1 = AsciiLower(c1);
                u2 = AsciiLower(qSeq[i]);
            }
            else
            {
                u1 = AsciiUpper(c1);
                u2 = AsciiUpper(qSeq[i]);
            }

            if (wildcard is int w && (u1 == w || u2 == w))
                n++;
            else if (u1 != u2)
                mis++;
            else if (mask is "lower" or "upper")
            {
                if (u1 == c1)
                    matches++;
                else
                    rep++;
            }
            else
                matches++;
        }

        return (matches, mis, rep, n);
    }

    private static byte[] WriteBigPslCore(IReadOnlyList<BioCoordAlignment> alignments, IReadOnlyList<BioCoordSequence> targets, BigBedWriteOptions o)
    {
        var fixedAlignments = new List<BioCoordAlignment>();
        foreach (var alignment in alignments)
        {
            if (alignment.Coordinates.Length == 0 || alignment.Coordinates[0].Length == 0)
                continue;
            if (alignment.Sequences.Count != 2)
                throw new ArgumentException("ValueError: too many values to unpack (expected 2)");
            var targetRecord = alignment.Sequences[0];
            var queryRecord = alignment.Sequences[1];
            var target = targetRecord;
            var query = queryRecord;
            long tSize = target.Length, qSize = query.Length;
            bool? dnax = null;
            var c = Copy(alignment.Coordinates);
            string strand;
            if (c[1][0] > c[1][^1])
            {
                strand = "-";
                query = ReverseComplement(query);
                for (int i = 0; i < c[1].Length; i++)
                    c[1][i] = qSize - c[1][i];
            }
            else if (c[0][0] > c[0][^1])
            {
                strand = "-";
                target = ReverseComplement(target);
                for (int i = 0; i < c[0].Length; i++)
                    c[0][i] = tSize - c[0][i];
                dnax = true;
            }
            else
            {
                strand = "+";
            }

            long matches = 0, misMatches = 0, repMatches = 0, nCount = 0;
            var blockSizes = new List<long>();
            var qStarts = new List<long>();
            long tStart = c[0][0], qStart = c[1][0], qCount = 0;
            for (int col = 1; col < c[0].Length; col++)
            {
                long tEnd = c[0][col], qEnd = c[1][col];
                if (tStart == tEnd)
                {
                    qStart = qEnd;
                }
                else if (qStart == qEnd)
                {
                    tStart = tEnd;
                }
                else
                {
                    long tCount = tEnd - tStart;
                    qCount = qEnd - qStart;
                    qStarts.Add(qStart);
                    blockSizes.Add(qCount);
                    if (tCount == qCount)
                    {
                        if (dnax == true)
                            throw new ArgumentException("AssertionError");
                        dnax = false;
                    }
                    else
                    {
                        if (tCount != 3 * qCount || dnax == false)
                            throw new ArgumentException("AssertionError");
                        dnax = true;
                    }

                    string? tSeq = PySliceBytes(target, tStart, tEnd);
                    string? qSeq = PySliceBytes(query, qStart, qEnd);
                    if (tSeq is null || qSeq is null)
                    {
                        matches += qCount;
                    }
                    else
                    {
                        var (m, mm, rm, nc) = CountBlock(tSeq, qSeq, o.Mask, o.Wildcard);
                        matches += m;
                        misMatches += mm;
                        repMatches += rm;
                        nCount += nc;
                    }

                    tStart = tEnd;
                    qStart = qEnd;
                }
            }

            string Override(string name, long counted) => alignment.Attributes.TryGetPropertyValue(name, out var v) && v is not null ? PyStr(v) : Str(counted);
            string matchesText = Override("matches", matches), misMatchText = Override("misMatches", misMatches), repMatchText = Override("repMatches", repMatches), nCountText = Override("nCount", nCount);
            if (qStarts.Count == 0)
                throw new ArgumentException("IndexError: index 0 is out of bounds for axis 0 with size 0");
            long oStart = qStarts[0], oEnd = qStarts[^1] + qCount;
            string oStrand = "+";
            long[][] coordinates = alignment.Coordinates;
            if (strand == "-")
            {
                if (dnax == true)
                {
                    oStrand = "-";
                    var flipped = new List<long>();
                    for (int i = 0; i < qStarts.Count; i++)
                        flipped.Add(qSize - (qStarts[i] + blockSizes[i]));
                    flipped.Reverse();
                    qStarts = flipped;
                    coordinates = [.. alignment.Coordinates.Select(r => r.Reverse().ToArray())];
                }
                else
                {
                    (oStart, oEnd) = (qSize - oEnd, qSize - oStart);
                }
            }

            string oSequence = "";
            if (o.Fa)
            {
                if (queryRecord.Parts is null || queryRecord.Parts.Count != 1 || !queryRecord.Parts.TryGetValue(0, out var text) || text.Length != queryRecord.Length)
                    throw new ArgumentException("UndefinedSequenceError: Sequence content is undefined");
                oSequence = text;
            }

            string oCds = "";
            if (o.Cds)
            {
                oCds = "n/a";
                foreach (var feature in queryRecord.Features)
                {
                    if (feature is JsonObject f && f["type"]?.GetValue<string>() == "CDS")
                    {
                        oCds = InsdcLocationString(f["location"]!, queryRecord.Length);
                        break;
                    }
                }
            }

            string molecule = queryRecord.Annotations.TryGetPropertyValue("molecule_type", out var mt) && mt is JsonValue mv && mv.TryGetValue<string>(out var mts) ? mts : "";
            string seqType = "0";
            if (molecule == "DNA")
                seqType = "1";
            else if (molecule == "protein")
                seqType = "2";
            if (!alignment.Attributes.TryGetPropertyValue("annotations", out var annNode) || annNode is not JsonObject)
                throw new ArgumentException("AttributeError: 'Alignment' object has no attribute 'annotations'");
            var annotations = (JsonObject)annNode.DeepClone();
            annotations["oChromStart"] = Str(oStart);
            annotations["oChromEnd"] = Str(oEnd);
            annotations["oStrand"] = oStrand;
            annotations["oChromSize"] = Str(qSize);
            annotations["oChromStarts"] = string.Join(",", qStarts.Select(Str));
            annotations["oSequence"] = oSequence;
            annotations["oCDS"] = oCds;
            annotations["chromSize"] = Str(tSize);
            annotations["match"] = matchesText;
            annotations["misMatch"] = misMatchText;
            annotations["repMatch"] = repMatchText;
            annotations["nCount"] = nCountText;
            annotations["seqType"] = seqType;
            var attributes = new JsonObject();
            foreach (var kv in alignment.Attributes)
                attributes[kv.Key] = kv.Key == "annotations" ? annotations : kv.Value?.DeepClone();
            fixedAlignments.Add(new BioCoordAlignment { Sequences = alignment.Sequences, Coordinates = coordinates, Attributes = attributes });
        }

        var sorted = fixedAlignments.OrderBy(a => a.Sequences[0].Id, StringComparer.Ordinal).ThenBy(a => a.Coordinates[0][0]).ToList();
        var inner = new BigBedWriteOptions
        {
            BedN = 12,
            Compress = o.Compress,
            Compressor = o.Compressor,
            DeclarationFields = BigPslDeclaration.Fields,
            DeclarationName = BigPslDeclaration.Name,
            DeclarationComment = BigPslDeclaration.Comment,
        };
        return WriteBigBedCore(sorted, targets, inner);
    }

    /// <summary><c>Bio.SeqIO.InsdcIO._insdc_location_string</c> of an exact-position location ({start, end, strand} or {parts, operator}).</summary>
    internal static string InsdcLocationString(JsonNode location, long recLength)
    {
        if (location is JsonObject obj && obj["parts"] is JsonArray parts)
        {
            string op = obj["operator"]!.GetValue<string>();
            int? strand = CommonStrand(parts);
            if (strand == -1)
                return $"complement({op}({string.Join(",", parts.Reverse().Select(p => InsdcSimpleLocation(p!, recLength)))}))";
            return $"{op}({string.Join(",", parts.Select(p => InsdcLocationString(p!, recLength)))})";
        }

        string loc = InsdcSimpleLocation(location, recLength);
        return location["strand"]?.GetValue<int>() == -1 ? $"complement({loc})" : loc;
    }

    private static int? CommonStrand(JsonArray parts)
    {
        int? first = null;
        for (int i = 0; i < parts.Count; i++)
        {
            int? s = parts[i]!["strand"]?.GetValue<int>();
            if (i == 0)
                first = s;
            else if (s != first)
                return null;
        }

        return first;
    }

    private static string InsdcSimpleLocation(JsonNode location, long recLength)
    {
        long start = location["start"]!.GetValue<long>(), end = location["end"]!.GetValue<long>();
        if (start == end)
            return end == recLength ? $"{recLength}^1" : $"{end}^{end + 1}";
        if (start + 1 == end)
            return Str(end);
        return $"{start + 1}..{end}";
    }

    private const string LocOneOfPosition = @"one\-of\(\d+[,\d+]+\)";
    private const string LocReference = @"(?:[a-zA-Z][a-zA-Z0-9_\.\|]*[a-zA-Z0-9]?\:)";
    private static readonly Regex LocPair = new(@"^[<>]?-?\d+\.\.[<>]?-?\d+", RegexOptions.ECMAScript | RegexOptions.CultureInvariant);
    private static readonly Regex LocBetween = new(@"^\d+\^\d+", RegexOptions.ECMAScript | RegexOptions.CultureInvariant);
    private static readonly Regex LocWithin = new(@"^([<>]?\d+|\(\d+\.\d+\))\.\.([<>]?\d+|\(\d+\.\d+\))", RegexOptions.ECMAScript | RegexOptions.CultureInvariant);
    private static readonly Regex LocOneOf = new(@"^([<>]?\d+|" + LocOneOfPosition + @")\.\.([<>]?\d+|" + LocOneOfPosition + ")", RegexOptions.ECMAScript | RegexOptions.CultureInvariant);
    private static readonly Regex LocBond = new(@"^bond\([<>]?\d+\)", RegexOptions.ECMAScript | RegexOptions.CultureInvariant);
    private static readonly Regex LocSolo = new(@"^[<>]?\d+", RegexOptions.ECMAScript | RegexOptions.CultureInvariant);
    private static readonly Regex LocWithinPosition = new(@"^\((\d+)\.(\d+)\)$", RegexOptions.ECMAScript | RegexOptions.CultureInvariant);
    private static readonly Regex LocOneOfChoices = new(@"^one\-of\((\d+[,\d+]+)\)$", RegexOptions.ECMAScript | RegexOptions.CultureInvariant);

    // the splitter of the parts of join() / order() / bond(): Bio.SeqFeature._split (_any_location)
    private const string LocOneOfLocation = @"[<>]?(?:\d+|" + LocOneOfPosition + @")\.\.[<>]?(?:\d+|" + LocOneOfPosition + ")";
    private static readonly Regex LocAnyLocation = new("(" + LocReference + "?" + LocOneOfLocation + @"|complement\(" + LocOneOfLocation + @"\)|[^,]+|complement\([^,]+\))", RegexOptions.ECMAScript | RegexOptions.CultureInvariant);

    /// <summary><c>Position.fromstring(text, offset)</c> as an integer: the position of a start (offset -1: the lowest choice / the left bound) or of an end (the highest choice / the right bound).</summary>
    private static long PositionValue(string text, bool start)
    {
        int offset = start ? -1 : 0;
        if (text.Length > 0 && text[0] is '<' or '>')
            return PyInt(text[1..]) + offset;
        var within = LocWithinPosition.Match(text);
        if (within.Success)
            return start ? PyInt(within.Groups[1].Value) + offset : PyInt(within.Groups[2].Value);
        var oneOf = LocOneOfChoices.Match(text);
        if (oneOf.Success)
        {
            var choices = oneOf.Groups[1].Value.Split(',').Select(c => PyInt(c) + offset).ToList();
            return start ? choices.Min() : choices.Max();
        }

        return PyInt(text) + offset;
    }

    private static JsonObject SimpleLocationSpec(long start, long end, int? strand)
        => new() { ["start"] = start, ["end"] = end, ["strand"] = strand };

    /// <summary><c>Location.fromstring</c>: exact, fuzzy (<c>&lt;</c>, <c>&gt;</c>), within <c>(a.b)</c> and <c>one-of(...)</c> positions, <c>a..b</c>, <c>a</c>, <c>a^b</c>, <c>ref:</c> prefixes (dropped: the description holds start, end and strand), <c>complement</c>, <c>join</c>, <c>order</c>, <c>bond</c>.</summary>
    internal static JsonObject LocationFromString(string text)
    {
        int strand = 1;
        if (text.StartsWith("complement(", StringComparison.Ordinal))
        {
            if (text.Length == 0 || text[^1] != ')')
                throw new ArgumentException($"closing bracket missing in '{text}'");
            text = text[11..^1];
            strand = -1;
        }

        string? op = Array.Find(["join", "order", "bond"], name => text.StartsWith(name + "(", StringComparison.Ordinal));
        if (op is null)
        {
            var loc = SimpleLocationFromString(text);
            loc["strand"] = strand;
            return loc;
        }

        string inner = text[(op.Length + 1)..Math.Max(op.Length + 1, text.Length - 1)];
        var locs = new List<JsonObject>();
        JsonObject last = [];
        foreach (Match part in LocAnyLocation.Matches(inner))
        {
            last = SimpleLocationFromString(part.Groups[1].Value);
            if (last["strand"]?.GetValue<int>() == -1)
            {
                if (strand == -1)
                    throw new ArgumentException($"double complement in '{text}'?");
            }
            else
            {
                last["strand"] = strand;
            }

            locs.Add((JsonObject)last.DeepClone());
        }

        if (locs.Count == 0)
            throw new ArgumentException($"failed to parse feature location '{text}'");
        if (locs.Count == 1)
            return last;
        if (strand == -1)
        {
            if (locs.Any(l => l["strand"]?.GetValue<int>() != -1))
                throw new ArgumentException("AssertionError");
            locs.Reverse();
        }

        return new JsonObject { ["parts"] = new JsonArray([.. locs.Select(l => (JsonNode)l)]), ["operator"] = op };
    }

    private static JsonObject SimpleLocationFromString(string text)
    {
        int? strand = null;
        if (text.StartsWith("complement(", StringComparison.Ordinal))
        {
            text = text.Length >= 12 ? text[11..^1] : "";
            strand = -1;
        }

        var halves = text.Split("..");
        if (halves.Length == 2 && long.TryParse(halves[0].Trim(), NumberStyles.AllowLeadingSign, CultureInfo.InvariantCulture, out long s0) && long.TryParse(halves[1].Trim(), NumberStyles.AllowLeadingSign, CultureInfo.InvariantCulture, out long e0)
            && 0 <= s0 - 1 && s0 - 1 < e0)
            return SimpleLocationSpec(s0 - 1, e0, strand);
        var colon = text.Split(':');
        if (colon.Length == 2)
            text = colon[1];
        (Regex Pattern, string Kind)[] categories = [(LocPair, "range"), (LocBetween, "between"), (LocWithin, "range"), (LocOneOf, "range"), (LocBond, "bond"), (LocSolo, "solo")];
        Match? match = null;
        string kind = "";
        foreach (var (pattern, k) in categories)
        {
            var candidate = pattern.Match(text);
            if (candidate.Success)
            {
                match = candidate;
                kind = k;
                break;
            }
        }

        if (match is null)
            throw new ArgumentException($"Could not parse feature location '{text}'");
        RequireWhole(match, text);
        long start, end;
        if (kind == "range")
        {
            (start, end) = PositionRange(text);
        }
        else if (kind == "between")
        {
            var p = text.Split('^');
            long s = PyInt(p[0]), e = PyInt(p[1]);
            if (s + 1 != e)
                throw new ArgumentException($"invalid feature location '{text}'");
            start = end = s;
        }
        else
        {
            string number = kind == "bond" ? text[5..^1] : text;
            start = PositionValue(number, true);
            end = PositionValue(number, false);
        }

        if (start < 0)
            throw new ArgumentException($"negative starting position in feature location '{text}'");
        if (end < start)
            throw new ArgumentException("End location must be greater than or equal to start location");
        return SimpleLocationSpec(start, end, strand);
    }

    private static void RequireWhole(Match m, string text)
    {
        if (m.Value != text)
            throw new ArgumentException("AssertionError");
    }

    private static (long Start, long End) PositionRange(string text)
    {
        var p = text.Split("..");
        return (PositionValue(p[0], true), PositionValue(p[1], false));
    }

    private static readonly string[] BigPslFieldNames =
    [
        "chrom", "chromStart", "chromEnd", "name", "score", "strand", "thickStart", "thickEnd", "reserved", "blockCount", "blockSizes", "chromStarts",
        "oChromStart", "oChromEnd", "oStrand", "oChromSize", "oChromStarts", "oSequence", "oCDS", "chromSize", "match", "misMatch", "repMatch", "nCount", "seqType",
    ];

    private static void BbCheckPslFields(IReadOnlyList<BigBedField> fields)
    {
        for (int i = 0; i < BigPslFieldNames.Length; i++)
        {
            if (fields[i].Name != BigPslFieldNames[i])
                throw new ArgumentException($"Expected field name '{BigPslFieldNames[i]}'; found '{fields[i].Name}'");
        }
    }

    /// <summary><c>bigpsl.AlignmentIterator._create_alignment</c>.</summary>
    private static BioCoordAlignment BbCreatePslAlignment(long chromId, long tStart, long tEnd, byte[] restBytes, List<BioCoordSequence> targets)
    {
        var words = Encoding.UTF8.GetString(restBytes).Split('\t');
        if (words.Length != 22)
            throw new ArgumentException($"Unexpected number of fields ({words.Length}, expected 22)");
        var targetRecord = targets[(int)chromId];
        long tSize = PyInt(words[16]);
        if (targetRecord.Length != tSize)
            throw new ArgumentException($"Unexpected chromosome size {tSize} (expected {targetRecord.Length})");
        string strand = words[2];
        string qName = words[0];
        long qSize = PyInt(words[12]);
        long blockCount = PyInt(words[6]);
        long[] blockSizes = [.. words[7].TrimEnd(',').Split(',').Select(PyInt)];
        long[] tStarts = [.. words[8].TrimEnd(',').Split(',').Select(PyInt)];
        long[] qStarts = [.. words[13].TrimEnd(',').Split(',').Select(PyInt)];
        if (blockSizes.Length != blockCount)
            throw new ArgumentException($"Inconsistent number of blocks ({blockSizes.Length} found, expected {blockCount})");
        if (qStarts.Length != blockCount)
            throw new ArgumentException("Inconsistent number of query start positions");
        if (tStarts.Length != blockCount)
            throw new ArgumentException("Inconsistent number of target start positions");
        long[] tBlockSizes = [.. blockSizes];
        string querySequence = words[14];
        SortedDictionary<long, string>? queryParts;
        if (querySequence.Length == 0)
        {
            queryParts = UndefinedParts(qSize);
        }
        else
        {
            if (querySequence.Length != qSize)
                throw new ArgumentException($"Inconsistent query sequence length ({querySequence.Length}, expected {qSize})");
            queryParts = new SortedDictionary<long, string> { [0] = querySequence };
        }

        var queryAnnotations = new JsonObject();
        var features = new JsonArray();
        string cds = words[15];
        if (cds.Length > 0 && cds != "n/a")
            features.Add(new JsonObject { ["type"] = "CDS", ["location"] = LocationFromString(cds), ["qualifiers"] = new JsonObject() });
        long[] qBlockSizes;
        switch (words[21])
        {
            case "0":
                qBlockSizes = tBlockSizes;
                break;
            case "1":
                queryAnnotations["molecule_type"] = "DNA";
                qBlockSizes = tBlockSizes;
                break;
            case "2":
                queryAnnotations["molecule_type"] = "protein";
                qBlockSizes = [.. tBlockSizes.Select(x => (long)Math.Floor(x / 3.0))];
                break;
            default:
                throw new ArgumentException($"Unexpected sequence type '{words[21]}'");
        }

        for (int i = 0; i < tStarts.Length; i++)
            tStarts[i] += tStart;
        string qStrand = words[11];
        if (qStrand == "-" && strand == "-")
        {
            (tStart, tEnd) = (tEnd, tStart);
            qStarts = [.. qStarts.Select((x, i) => qSize - x - qBlockSizes[i]).Reverse()];
            tStarts = [.. tStarts.Select((x, i) => tSize - x - tBlockSizes[i]).Reverse()];
            qBlockSizes = [.. qBlockSizes.Reverse()];
            tBlockSizes = [.. tBlockSizes.Reverse()];
        }

        long qPosition = qStarts[0], tPosition = tStarts[0];
        var t = new List<long> { tPosition };
        var q = new List<long> { qPosition };
        for (int i = 0; i < tBlockSizes.Length; i++)
        {
            if (tStarts[i] != tPosition)
            {
                t.Add(tStarts[i]);
                q.Add(qPosition);
                tPosition = tStarts[i];
            }

            if (qStarts[i] != qPosition)
            {
                t.Add(tPosition);
                q.Add(qStarts[i]);
                qPosition = qStarts[i];
            }

            tPosition += tBlockSizes[i];
            qPosition += qBlockSizes[i];
            t.Add(tPosition);
            q.Add(qPosition);
        }

        long[][] coordinates = [[.. t], [.. q]];
        long qStart = PyInt(words[9]), qEnd = PyInt(words[10]);
        if (strand == "-")
        {
            if (qStrand == "-")
            {
                for (int i = 0; i < coordinates[0].Length; i++)
                    coordinates[0][i] = tSize - coordinates[0][i];
            }
            else
            {
                (qStart, qEnd) = (qEnd, qStart);
                for (int i = 0; i < coordinates[1].Length; i++)
                    coordinates[1][i] = qSize - coordinates[1][i];
            }
        }

        if (tStart != coordinates[0][0])
            throw new ArgumentException($"Inconsistent tStart found ({tStart}, expected {coordinates[0][0]})");
        if (tEnd != coordinates[0][^1])
            throw new ArgumentException($"Inconsistent tEnd found ({tEnd}, expected {coordinates[0][^1]})");
        if (qStart != coordinates[1][0])
            throw new ArgumentException($"Inconsistent qStart found ({qStart}, expected {coordinates[1][0]})");
        if (qEnd != coordinates[1][^1])
            throw new ArgumentException($"Inconsistent qEnd found ({qEnd}, expected {coordinates[1][^1]})");
        var tgt = new BioCoordSequence { Id = targetRecord.Id, Length = targetRecord.Length, Parts = targetRecord.Parts, Description = targetRecord.Description };
        var query = new BioCoordSequence { Id = qName, Length = qSize, Parts = queryParts, Description = BiopythonUnknownDescription, Annotations = queryAnnotations, Features = features };
        var attrs = new JsonObject { ["annotations"] = new JsonObject() };
        string scoreText = words[1];
        if (double.TryParse(scoreText, NumberStyles.Float, CultureInfo.InvariantCulture, out double score) && !scoreText.Contains('_'))
            attrs["score"] = score == Math.Floor(score) && !double.IsInfinity(score) ? (long)score : score;
        else
            attrs["score"] = scoreText;
        attrs["thickStart"] = PyInt(words[3]);
        attrs["thickEnd"] = PyInt(words[4]);
        attrs["itemRgb"] = words[5];
        attrs["matches"] = PyInt(words[17]);
        attrs["misMatches"] = PyInt(words[18]);
        attrs["repMatches"] = PyInt(words[19]);
        attrs["nCount"] = PyInt(words[20]);
        return new BioCoordAlignment { Sequences = [tgt, query], Coordinates = coordinates, Attributes = attrs };
    }

    /// <summary>Reads a bigPsl file as <c>Bio.Align.parse(handle, "bigpsl")</c>.</summary>
    public static CoordinateAlignmentFile ParseBigPsl(byte[] data)
    {
        ArgumentNullException.ThrowIfNull(data);
        try
        {
            return ParseBigBedCore(data, false, true);
        }
        catch (Exception e) when (e is IndexOutOfRangeException or ArgumentOutOfRangeException or KeyNotFoundException or InvalidCastException or NullReferenceException or FormatException or OverflowException or InvalidOperationException or InvalidDataException)
        {
            throw new ArgumentException(e.Message, e);
        }
    }
}
