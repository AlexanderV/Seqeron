using System.Buffers.Binary;
using System.Globalization;
using System.IO.Compression;
using System.Text;
using System.Text.Json.Nodes;
using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Alignment;

public static partial class AlignmentCoordinateFormats
{
    /// <summary>The default declaration (<c>AutoSQLTable.default</c>): the 12 fields of BED.</summary>
    public static IReadOnlyList<BigBedField> DefaultBigBedFields { get; } =
    [
        new("string", "chrom", "Reference sequence chromosome or scaffold"),
        new("uint", "chromStart", "Start position in chromosome"),
        new("uint", "chromEnd", "End position in chromosome"),
        new("string", "name", "Name of item."),
        new("uint", "score", "Score (0-1000)"),
        new("char[1]", "strand", "+ or - for strand"),
        new("uint", "thickStart", "Start of where display should be thick (start codon)"),
        new("uint", "thickEnd", "End of where display should be thick (stop codon)"),
        new("uint", "reserved", "Used as itemRgb as of 2004-11-22"),
        new("int", "blockCount", "Number of blocks"),
        new("int[blockCount]", "blockSizes", "Comma separated list of block sizes"),
        new("int[blockCount]", "chromStarts", "Start positions relative to chromStart"),
    ];

    /// <summary><c>str(AutoSQLTable)</c>: the autoSql text of a declaration.</summary>
    public static string BigBedDeclarationText(string name, string comment, IReadOnlyList<BigBedField> fields)
    {
        int typeWidth = fields.Max(f => f.AsType.Length);
        int nameWidth = fields.Max(f => f.Name.Length) + 1;
        var sb = new StringBuilder();
        sb.Append("table ").Append(name).Append('\n');
        sb.Append('"').Append(comment).Append("\"\n");
        sb.Append("(\n");
        foreach (var f in fields)
            sb.Append("   ").Append(f.AsType.PadRight(typeWidth)).Append(' ').Append((f.Name + ";").PadRight(nameWidth)).Append("    \"").Append(f.Comment).Append("\"\n");
        sb.Append(")\n");
        return sb.ToString();
    }

    private sealed class ByteSink
    {
        public readonly MemoryStream Stream = new();

        public long Position => Stream.Position;

        public void Write(byte[] data) => Stream.Write(data, 0, data.Length);

        public void U8(int v) => Stream.WriteByte((byte)v);

        public void U16(int v) => Write(BitConverter.GetBytes((ushort)v));

        public void U32(long v) => Write(BitConverter.GetBytes(unchecked((uint)v)));

        public void I32(long v) => Write(BitConverter.GetBytes(unchecked((int)v)));

        public void U64(long v) => Write(BitConverter.GetBytes(unchecked((ulong)v)));

        public void F32(float v) => Write(BitConverter.GetBytes(v));

        public void F64(double v) => Write(BitConverter.GetBytes(v));

        public void Zeros(long n)
        {
            if (n > 0)
                Write(new byte[n]);
        }

        public void Seek(long position) => Stream.Position = position;
    }

    private sealed class BbRegion
    {
        public long ChromId, Start, End, Offset;
    }

    private sealed class BbSummary
    {
        public long ChromId, Start, End, Offset;
        public long ValidCount;
        public float MinVal, MaxVal, SumData, SumSquares;

        public BbSummary(long chromId, long start, long end, long value)
        {
            ChromId = chromId;
            Start = start;
            End = end;
            MinVal = value;
            MaxVal = value;
        }

        public BbSummary Copy() => (BbSummary)MemberwiseClone();

        public void Add(BbSummary other)
        {
            End = other.End;
            ValidCount += other.ValidCount;
            MinVal = Math.Min(MinVal, other.MinVal);
            MaxVal = Math.Max(MaxVal, other.MaxVal);
            SumData += other.SumData;
            SumSquares += other.SumSquares;
        }

        public void Update(long overlap, long val)
        {
            ValidCount += overlap;
            if (MinVal > val)
                MinVal = val;
            if (MaxVal < val)
                MaxVal = val;
            SumData += (float)(val * overlap);
            SumSquares += (float)(val * val * overlap);
        }

        public byte[] ToBytes()
        {
            var s = new ByteSink();
            s.U32(ChromId);
            s.U32(Start);
            s.U32(End);
            s.U32(ValidCount);
            s.F32(MinVal);
            s.F32(MaxVal);
            s.F32(SumData);
            s.F32(SumSquares);
            return s.Stream.ToArray();
        }
    }

    private sealed class BbRange
    {
        public long Start, End, Val;
    }

    /// <summary>The red-black interval tree of <c>_RangeTree</c> as a sorted list (the in-order traversal is all that the summaries use).</summary>
    private sealed class BbRangeTree
    {
        public readonly List<BbRange> Ranges = [];
        public readonly long ChromId;
        public readonly long ChromSize;

        public BbRangeTree(long chromId, long chromSize)
        {
            ChromId = chromId;
            ChromSize = chromSize;
        }

        private BbRange? Find(long start, long end) => Ranges.FirstOrDefault(r => end > r.Start && r.End > start);

        private void Add(BbRange item)
        {
            if (Ranges.Any(r => item.End > r.Start && r.End > item.Start))
                return;
            int i = Ranges.FindIndex(r => r.Start > item.Start);
            Ranges.Insert(i < 0 ? Ranges.Count : i, item);
        }

        public void AddToCoverageDepth(BioCoordAlignment alignment)
        {
            long start = alignment.Coordinates[0][0];
            long end = alignment.Coordinates[0][^1];
            if (start > end)
                (start, end) = (end, start);
            var existing = Find(start, end);
            if (existing is null)
            {
                Add(new BbRange { Start = start, End = end, Val = 1 });
                return;
            }

            if (existing.Start <= start && existing.End >= end)
            {
                if (existing.Start < start)
                {
                    var r = new BbRange { Start = existing.Start, End = start, Val = existing.Val };
                    existing.Start = start;
                    Add(r);
                }

                if (existing.End > end)
                {
                    var r = new BbRange { Start = end, End = existing.End, Val = existing.Val };
                    existing.End = end;
                    Add(r);
                }

                existing.Val++;
                return;
            }

            var items = Ranges.Where(r => r.End > start && end > r.Start).ToList();
            long s = start;
            long e = end;
            foreach (var item in items)
            {
                if (s < item.Start)
                {
                    var r = new BbRange { Start = s, End = item.Start, Val = 1 };
                    s = item.Start;
                    Add(r);
                }
                else if (s > item.Start)
                {
                    var r = new BbRange { Start = item.Start, End = s, Val = item.Val };
                    item.Start = s;
                    Add(r);
                }

                if (item.Start < end && item.End > end)
                {
                    var r = new BbRange { Start = end, End = item.End, Val = item.Val };
                    item.End = end;
                    Add(r);
                }

                item.Val++;
                s = item.End;
            }

            if (s < e)
                Add(new BbRange { Start = s, End = e, Val = 1 });
        }

        /// <summary><c>_RangeTree.generate_summaries</c>.</summary>
        public IEnumerable<BbSummary> GenerateSummaries(long scale, BbTotalSum totalSum)
        {
            if (Ranges.Count == 0)
                throw new InvalidOperationException("StopIteration");
            int idx = 0;
            long start = Ranges[0].Start, end = Ranges[0].End, val = Ranges[0].Val;
            var summary = new BbSummary(ChromId, start, Math.Min(start + scale, ChromSize), val);
            while (true)
            {
                long size = Math.Max(end - start, 1);
                totalSum.Update(size, val);
                if (summary.End <= start)
                    summary = new BbSummary(ChromId, start, Math.Min(start + scale, ChromSize), val);
                while (end > summary.End)
                {
                    long overlap = Math.Min(end, summary.End) - Math.Max(start, summary.Start);
                    if (overlap <= 0)
                        throw new InvalidOperationException("AssertionError");
                    summary.Update(overlap, val);
                    size -= overlap;
                    start = summary.End;
                    yield return summary;
                    summary = new BbSummary(ChromId, start, Math.Min(start + scale, ChromSize), val);
                }

                summary.Update(size, val);
                idx++;
                if (idx >= Ranges.Count)
                    break;
                start = Ranges[idx].Start;
                end = Ranges[idx].End;
                val = Ranges[idx].Val;
                if (summary.End <= start)
                    yield return summary;
            }

            yield return summary;
        }
    }

    private sealed class BbTotalSum
    {
        public long ValidCount;
        public double MinVal = long.MaxValue;
        public double MaxVal = -long.MaxValue;
        public double SumData;
        public double SumSquares;

        public void Update(long size, long val)
        {
            ValidCount += size;
            if (val < MinVal)
                MinVal = val;
            if (val > MaxVal)
                MaxVal = val;
            SumData += (double)val * size;
            SumSquares += (double)val * val * size;
        }

        public byte[] ToBytes()
        {
            var s = new ByteSink();
            s.U64(ValidCount);
            s.F64(MinVal);
            s.F64(MaxVal);
            s.F64(SumData);
            s.F64(SumSquares);
            return s.Stream.ToArray();
        }
    }

    private sealed class BbRNode
    {
        public readonly List<BbRNode> Children = [];
        public BbRNode? Parent;
        public long StartChromId, StartBase, EndChromId, EndBase, StartFileOffset, EndFileOffset;

        public void CalcLevelSizes(long[] levelSizes, int level)
        {
            levelSizes[level]++;
            level++;
            if (level == levelSizes.Length)
                return;
            foreach (var child in Children)
                child.CalcLevelSizes(levelSizes, level);
        }
    }

    private const uint RTreeSignature = 0x2468ACE0;
    private const uint BPlusSignature = 0x78CA8C91;
    private const uint BigBedSignature = 0x8789F2EB;

    private sealed record ChromUsage(string Name, long Id, long Size);

    /// <summary>Python's <c>0 &lt;= int(text)</c> of a color component.</summary>
    private static long PyIntStrict(string s) => PyInt(s);

    /// <summary><c>bigbed.AlignmentWriter._extract_fields</c>.</summary>
    private static (string, long, long, byte[]) BbExtractFields(BioCoordAlignment alignment, BigBedWriteOptions o, IReadOnlyList<BigBedField> declaration)
    {
        int bedN = o.BedN;
        var row = new List<string>();
        string chrom = alignment.Sequences[0].Id ?? throw new ArgumentException("TypeError: object of type 'NoneType' has no len()");
        if (chrom.Length >= 255)
            throw new ArgumentException($"alignment target name '{chrom}' is too long (must not exceed 254 characters)");
        if (chrom.Length < 1)
            throw new ArgumentException("alignment target name cannot be blank or empty");
        long chromStart = alignment.Coordinates[0][0];
        long chromEnd = alignment.Coordinates[0][^1];
        if (chromEnd < chromStart)
            throw new ArgumentException($"chromStart after chromEnd ({chromEnd} > {chromStart})");
        string name = "";
        if (bedN > 3)
        {
            name = alignment.Sequences[1].Id ?? throw new ArgumentException("TypeError: object of type 'NoneType' has no len()");
            if (name == "")
                name = ".";
            else if (name.Length > 255)
                throw new ArgumentException($"alignment query name '{name}' is too long (must not exceed 255 characters");
            row.Add(name);
        }

        if (bedN > 4)
        {
            string score = ".";
            if (TryAttr(alignment, "score", out var sc))
            {
                var v = sc.AsValue();
                double d = v.TryGetValue<long>(out var l) ? l : v.GetValue<double>();
                if (d < 0 || d > 1000)
                    throw new ArgumentException($"score ({PyStr(sc)}) must be between 0 and 1000");
                score = PyStr(sc);
            }

            row.Add(score);
        }

        if (bedN > 5)
            row.Add(alignment.Coordinates[1][0] <= alignment.Coordinates[1][^1] ? "+" : "-");
        long thickStart = chromStart;
        if (bedN > 6)
        {
            if (TryAttr(alignment, "thickStart", out var ts))
                thickStart = Int(ts);
            row.Add(Str(thickStart));
        }

        if (bedN > 7)
        {
            long thickEnd = chromEnd;
            if (TryAttr(alignment, "thickEnd", out var te))
            {
                thickEnd = Int(te);
                if (thickEnd < thickStart)
                    throw new ArgumentException($"thickStart ({thickStart}) after thickEnd ({thickEnd})");
                if (thickStart != 0 && (thickStart < chromStart || thickStart > chromEnd))
                    throw new ArgumentException($"thickStart out of range for {name}:{chromStart}-{chromEnd}, thick:{thickStart}-{thickEnd}");
                if (thickEnd != 0 && (thickEnd < chromStart || thickEnd > chromEnd))
                    throw new ArgumentException($"thickEnd out of range for {name}:{chromStart}-{chromEnd}, thick:{thickStart}-{thickEnd}");
            }

            row.Add(Str(thickEnd));
        }

        if (bedN > 8)
        {
            string itemRgb = ".";
            if (TryAttr(alignment, "itemRgb", out var rgb))
            {
                itemRgb = rgb.GetValue<string>();
                var colors = itemRgb.TrimEnd(',').Split(',');
                bool ok = colors.Length == 3 && colors.All(c =>
                {
                    long x = PyIntStrict(c);
                    return x >= 0 && x < 256;
                });
                if (!ok)
                {
                    long x = PyIntStrict(itemRgb);
                    if (!(x >= 0 && x < (2L << 32)))
                        throw new ArgumentException($"Expecting color to consist of r,g,b values from 0 to 255. Got '{itemRgb}'");
                }
            }

            row.Add(itemRgb);
        }

        var blockSizes = new List<long>();
        var chromStarts = new List<long>();
        if (bedN > 9)
        {
            int m = alignment.Coordinates[0].Length;
            for (int c = 0; c + 1 < m; c++)
            {
                long ts = alignment.Coordinates[0][c + 1] - alignment.Coordinates[0][c];
                long qs = alignment.Coordinates[1][c + 1] - alignment.Coordinates[1][c];
                if (ts != 0 && qs != 0)
                {
                    blockSizes.Add(Math.Max(ts, qs));
                    chromStarts.Add(alignment.Coordinates[0][c] - chromStart);
                }
            }

            row.Add(Str(blockSizes.Count));
        }

        if (bedN > 10)
            row.Add(string.Join(",", blockSizes.Select(Str)) + ",");
        if (bedN > 11)
            row.Add(string.Join(",", chromStarts.Select(Str)) + ",");
        alignment.Attributes.TryGetPropertyValue("annotations", out var annotationsNode);
        foreach (var field in declaration.Skip(bedN))
        {
            if (annotationsNode is not JsonObject ann || !ann.TryGetPropertyValue(field.Name, out var value) || value is null)
                throw new ArgumentException($"KeyError: '{field.Name}'");
            if (value is JsonValue jv && jv.TryGetValue<string>(out var sv))
                row.Add(sv);
            else if (value is JsonValue)
                row.Add(PyStr(value));
            else if (value is JsonArray arr)
                row.Add(string.Join(",", arr.Select(x => PyValueStr(x))));
            else
                throw new ArgumentException("TypeError");
        }

        return (chrom, chromStart, chromEnd, Encoding.UTF8.GetBytes(string.Join("\t", row)));
    }

    /// <summary>Writes the alignments as <c>Bio.Align.write(alignments, handle, "bigbed", ...)</c> does (the bytes of the file).</summary>
    /// <param name="alignments">The alignments, sorted by target name and start.</param>
    /// <param name="targets">The chromosomes in the order of the alignments (<c>alignments.targets</c>): identifier and length.</param>
    /// <param name="options">The writer options.</param>
    public static byte[] WriteBigBed(IReadOnlyList<BioCoordAlignment> alignments, IReadOnlyList<BioCoordSequence> targets, BigBedWriteOptions? options = null)
    {
        ArgumentNullException.ThrowIfNull(alignments);
        ArgumentNullException.ThrowIfNull(targets);
        options ??= new BigBedWriteOptions();
        try
        {
            return WriteBigBedCore(alignments, targets, options);
        }
        catch (Exception e) when (e is IndexOutOfRangeException or ArgumentOutOfRangeException or KeyNotFoundException or InvalidCastException or NullReferenceException or FormatException or OverflowException or InvalidOperationException)
        {
            throw new ArgumentException(e.Message, e);
        }
    }

    private static byte[] DefaultCompress(byte[] data)
    {
        using var ms = new MemoryStream();
        using (var z = new ZLibStream(ms, CompressionLevel.Optimal, true))
            z.Write(data, 0, data.Length);
        return ms.ToArray();
    }

    private static byte[] WriteBigBedCore(IReadOnlyList<BioCoordAlignment> alignments, IReadOnlyList<BioCoordSequence> targets, BigBedWriteOptions o)
    {
        if (o.BedN < 3 || o.BedN > 12)
            throw new ArgumentException("bedN must be between 3 and 12");
        var compress = o.Compressor ?? DefaultCompress;
        var fields = o.DeclarationFields is not null ? [.. o.DeclarationFields] : DefaultBigBedFields.Take(o.BedN).ToList();
        string declName = o.DeclarationFields is not null ? o.DeclarationName : "bed";
        string declComment = o.DeclarationFields is not null ? o.DeclarationComment : "Browser Extensible Data";
        var extra = o.ExtraIndex.Select(n =>
        {
            int index = fields.FindIndex(f => f.Name == n);
            if (index < 0)
                throw new ArgumentException($"extraIndex field {n} not a standard bed field or found in 'as' file.");
            if (fields[index].AsType != "string")
                throw new ArgumentException("Sorry for now can only index string fields.");
            return (Name: n, Index: index);
        }).ToList();

        // _get_chrom_usage
        var usage = new List<ChromUsage>();
        long bedCount = 0, totalBases = 0, chromId = 0, chromSize = -1;
        string name = "";
        int keySize = 0;
        long lastStart = -1;
        var maxFieldSize = new long[extra.Count];
        foreach (var alignment in alignments)
        {
            string chrom = alignment.Sequences[0].Id ?? throw new ArgumentException("TypeError");
            long start = alignment.Coordinates[0][0];
            long end = alignment.Coordinates[0][^1];
            for (int x = 0; x < extra.Count; x++)
                maxFieldSize[x] = Math.Max(maxFieldSize[x], BbExtraValue(alignment, extra[x].Name).Length);
            if (start > end)
                throw new ArgumentException($"end ({end}) before start ({start}) in alignment [{bedCount}]");
            bedCount++;
            totalBases += end - start;
            if (name != chrom)
            {
                if (string.CompareOrdinal(name, chrom) > 0)
                    throw new ArgumentException($"alignments are not sorted by target name at alignment [{bedCount}]");
                if (name.Length > 0)
                {
                    usage.Add(new ChromUsage(name, chromId, chromSize));
                    chromId++;
                }

                var target = targets.FirstOrDefault(t => t.Id == chrom) ?? throw new ArgumentException($"failed to find target '{chrom}' in target list at alignment [{bedCount}]");
                name = chrom;
                keySize = Math.Max(keySize, Encoding.UTF8.GetByteCount(chrom));
                chromSize = target.Length;
                lastStart = -1;
            }

            if (end > chromSize)
                throw new ArgumentException($"end coordinate {end} bigger than {chrom} size of {chromSize} at alignment [{bedCount}]'");
            if (lastStart >= 0)
            {
                long diff = start - lastStart;
                if (diff < 0)
                    throw new ArgumentException($"alignments are not sorted at alignment [{bedCount}]");
            }

            lastStart = start;
        }

        if (name.Length > 0)
            usage.Add(new ChromUsage(name, chromId, chromSize));
        double aveSize = bedCount > 0 ? (double)totalBases / bedCount : 0;
        if (alignments.Count == 0)
            throw new ArgumentException("an empty alignment list cannot be written as bigBed");

        var output = new ByteSink();
        int headerSize = 4 + 2 + 2 + 8 + 8 + 8 + 2 + 2 + 8 + 8 + 4 + 8;
        output.Zeros(headerSize);
        output.Zeros(10 * 24);
        long autoSqlOffset = output.Position;
        output.Write(Encoding.UTF8.GetBytes(BigBedDeclarationText(declName, declComment, fields) + "\0"));
        long totalSummaryOffset = output.Position;
        output.Zeros(8 + 8 * 4);
        long extraIndicesOffset = output.Position;
        output.Zeros(64 + (extra.Count * 20));
        long chromosomeTreeOffset = output.Position;
        BbWriteBPlusTree(output, usage.Select(u => (Encoding.UTF8.GetBytes(u.Name), u.Id, u.Size)).ToList(), keySize, Math.Min(o.BlockSize, usage.Count), 8);
        long fullDataOffset = output.Position;
        var reductions = BbCalculateReductions(aveSize);
        output.U64(alignments.Count);

        // extra index chunks
        var chunks = extra.Select(_ => new (byte[] Name, long Offset, long Size)[alignments.Count]).ToList();
        for (int x = 0; x < extra.Count; x++)
            for (int i = 0; i < chunks[x].Length; i++)
                chunks[x][i] = (new byte[maxFieldSize[x]], 0, 0);

        // write_alignments
        int itemsPerSlot = o.ItemsPerSlot;
        long curChromId = -1;
        long itemIx = 0, sectionStartIx = 0, sectionEndIx = 0;
        string? currentChrom = null;
        var regions = new List<BbRegion>();
        var buffer = new ByteSink();
        long maxBlockSize = 0;
        bool done = false;
        BbRegion? region = null;
        int pos = 0;
        long start0 = 0, end0 = 0;
        while (true)
        {
            byte[] rest = [];
            BioCoordAlignment? alignment = null;
            if (pos < alignments.Count)
            {
                alignment = alignments[pos++];
                var (chrom, st, en, rs) = BbExtractFields(alignment, o, fields);
                start0 = st;
                end0 = en;
                rest = rs;
                if (chrom != currentChrom)
                {
                    if (currentChrom is not null)
                        itemIx = itemsPerSlot;
                    currentChrom = chrom;
                    curChromId++;
                    foreach (var r in reductions)
                        r.End = 0;
                }
            }
            else
            {
                itemIx = itemsPerSlot;
                done = true;
            }

            if (itemIx == itemsPerSlot)
            {
                long blockStartOffset = output.Position;
                long size = buffer.Position;
                if (size > maxBlockSize)
                    maxBlockSize = size;
                byte[] data = buffer.Stream.ToArray();
                if (o.Compress)
                    data = compress(data);
                output.Write(data);
                buffer.Stream.SetLength(0);
                buffer.Seek(0);
                if (extra.Count > 0)
                {
                    long blockEndOffset = output.Position;
                    long blockSize = blockEndOffset - blockStartOffset;
                    for (int x = 0; x < extra.Count; x++)
                    {
                        for (long k = sectionStartIx; k < sectionEndIx; k++)
                            chunks[x][k] = (chunks[x][k].Name, blockStartOffset, blockSize);
                    }

                    sectionStartIx = sectionEndIx;
                }

                region!.Offset = blockStartOffset;
                if (done)
                    break;
                itemIx = 0;
            }

            if (itemIx == 0)
            {
                region = new BbRegion { ChromId = curChromId, Start = start0, End = end0 };
                regions.Add(region);
            }
            else if (end0 > region!.End)
            {
                region.End = end0;
            }

            itemIx++;
            foreach (var row in reductions)
            {
                if (start0 >= row.End)
                {
                    row.Size++;
                    row.End = (int)(start0 + row.Scale);
                }

                while (end0 > row.End)
                {
                    row.Size++;
                    row.End += row.Scale;
                }
            }

            if (extra.Count > 0)
            {
                for (int x = 0; x < extra.Count; x++)
                {
                    byte[] v = BbExtraValue(alignment!, extra[x].Name);
                    var nameBytes = new byte[maxFieldSize[x]];
                    Array.Copy(v, nameBytes, Math.Min(v.Length, nameBytes.Length));
                    chunks[x][sectionEndIx] = (nameBytes, chunks[x][sectionEndIx].Offset, chunks[x][sectionEndIx].Size);
                }

                sectionEndIx++;
            }

            buffer.U32(curChromId);
            buffer.U32(start0);
            buffer.U32(end0);
            buffer.Write(rest);
            buffer.U8(0);
        }

        long uncompressBufSize = o.Compress ? Math.Max(maxBlockSize, (long)itemsPerSlot * 32) : 0;
        long fullIndexOffset = output.Position;
        BbWriteRTree(output, regions.Select(r => r).ToList(), o.BlockSize, 1, fullIndexOffset);
        var (zoomList, totalSum) = BbWriteZoomLevels(alignments, output, fullIndexOffset - fullDataOffset, usage, reductions, o, compress);
        var extraOffsets = new long[extra.Count];
        for (int x = 0; x < extra.Count; x++)
        {
            extraOffsets[x] = output.Position;
            var sorted = chunks[x].OrderBy(c => c.Name, Comparer<byte[]>.Create(CompareBytes)).ThenBy(c => c.Offset).ThenBy(c => c.Size).ToList();
            BbWriteBPlusTreeChunks(output, sorted, (int)maxFieldSize[x], o.BlockSize);
        }

        output.Seek(0);
        output.U32(BigBedSignature);
        output.U16(4);
        output.U16(zoomList.Count);
        output.U64(chromosomeTreeOffset);
        output.U64(fullDataOffset);
        output.U64(fullIndexOffset);
        output.U16(fields.Count);
        output.U16(o.BedN);
        output.U64(autoSqlOffset);
        output.U64(totalSummaryOffset);
        output.U32(uncompressBufSize);
        output.U64(extraIndicesOffset);
        foreach (var z in zoomList)
        {
            output.U32(z.ReductionLevel);
            output.Zeros(4);
            output.U64(z.DataOffset);
            output.U64(z.IndexOffset);
        }

        output.Seek(totalSummaryOffset);
        output.Write(totalSum.ToBytes());
        if (output.Position != extraIndicesOffset)
            throw new ArgumentException("AssertionError");
        if (extra.Count > 0)
        {
            long off = output.Position + 64;
            output.U16(64);
            output.U16(extra.Count);
            output.U64(off);
            output.Zeros(52);
            for (int x = 0; x < extra.Count; x++)
            {
                output.Zeros(2);
                output.U16(1);
                output.U64(extraOffsets[x]);
                output.Zeros(4);
                output.U16(extra[x].Index);
                output.Zeros(2);
            }
        }
        else
        {
            output.U16(64);
            output.U16(0);
            output.U64(0);
            output.Zeros(52);
        }

        output.Seek(output.Stream.Length);
        output.U32(BigBedSignature);
        return output.Stream.ToArray();
    }

    private static int CompareBytes(byte[]? a, byte[]? b)
    {
        for (int i = 0; i < Math.Min(a!.Length, b!.Length); i++)
        {
            if (a[i] != b[i])
                return a[i].CompareTo(b[i]);
        }

        return a.Length.CompareTo(b.Length);
    }

    private static byte[] BbExtraValue(BioCoordAlignment alignment, string name)
    {
        string v;
        if (name == "chrom")
            v = alignment.Sequences[0].Id ?? "";
        else if (name == "name")
            v = alignment.Sequences[1].Id ?? "";
        else if (alignment.Attributes.TryGetPropertyValue("annotations", out var an) && an is JsonObject ao && ao.TryGetPropertyValue(name, out var val) && val is not null)
            v = val.GetValue<string>();
        else
            throw new ArgumentException($"KeyError: '{name}'");
        return Encoding.UTF8.GetBytes(v);
    }

    private sealed class BbReduction
    {
        public int Scale, Size, End;
    }

    private static List<BbReduction> BbCalculateReductions(double aveSize)
    {
        const int maxLevels = 10;
        var all = new List<BbReduction>();
        for (int i = 0; i < maxLevels; i++)
            all.Add(new BbReduction());
        long res = Math.Max((long)aveSize, 10);
        int resTry = 0;
        for (; resTry < maxLevels; resTry++)
        {
            if (res > int.MaxValue)
                break;
            all[resTry].Scale = (int)res;
            res *= 4;
        }

        if (resTry == maxLevels)
            resTry = maxLevels - 1;
        return all.Take(resTry).ToList();
    }

    private sealed class BbZoom
    {
        public long ReductionLevel, DataOffset, IndexOffset;
    }

    private static (List<BbZoom> Zooms, BbTotalSum Total) BbWriteZoomLevels(IReadOnlyList<BioCoordAlignment> alignments, ByteSink output, long dataSize, List<ChromUsage> usage, List<BbReduction> reductions, BigBedWriteOptions o, Func<byte[], byte[]> compress)
    {
        var zoomList = new List<BbZoom> { new() };
        var totalSum = new BbTotalSum();
        double maxReducedSize = dataSize / 2.0;
        zoomList[0].DataOffset = output.Position;
        BbReduction? initialReduction = null;
        foreach (var r in reductions)
        {
            initialReduction = r;
            double reducedSize = r.Size * 32.0;
            if (o.Compress)
                reducedSize /= 2;
            if (reducedSize <= maxReducedSize)
                break;
            initialReduction = null;
        }

        initialReduction ??= reductions[0];
        output.I32(initialReduction.Size);
        int blockBytes = o.ItemsPerSlot * 32;
        var buf = new ByteSink();
        void Flush(ByteSink sink, bool zip)
        {
            byte[] data = buf.Stream.ToArray();
            sink.Write(zip ? compress(data) : data);
            buf.Stream.SetLength(0);
            buf.Seek(0);
        }

        void WriteItem(ByteSink sink, BbSummary item, int size, bool zip)
        {
            item.Offset = sink.Position;
            buf.Write(item.ToBytes());
            if (buf.Position == size)
                Flush(sink, zip);
        }

        var regions = new List<BbSummary>();
        var rezoomedList = new List<BbSummary>();
        long scale = initialReduction.Scale;
        long doubleReductionSize = scale * 4;
        // the trees: one per chromosome, built from the alignments in order
        var trees = new List<BbRangeTree>();
        int pos = 0;
        foreach (var u in usage)
        {
            var tree = new BbRangeTree(u.Id, u.Size);
            while (pos < alignments.Count && alignments[pos].Sequences[0].Id == u.Name)
                tree.AddToCoverageDepth(alignments[pos++]);
            trees.Add(tree);
        }

        BbSummary? rezoomed = null;
        foreach (var tree in trees)
        {
            long start = -long.MaxValue;
            foreach (var summary in tree.GenerateSummaries(scale, totalSum))
            {
                WriteItem(output, summary, blockBytes, o.Compress);
                regions.Add(summary);
                if (start + doubleReductionSize < summary.End)
                {
                    rezoomed = summary.Copy();
                    start = rezoomed.Start;
                    rezoomedList.Add(rezoomed);
                }
                else
                {
                    rezoomed!.Add(summary);
                }
            }
        }

        Flush(output, o.Compress);
        if (regions.Count != initialReduction.Size)
            throw new ArgumentException("AssertionError");
        zoomList[0].ReductionLevel = initialReduction.Scale;
        long indexOffset = output.Position;
        zoomList[0].IndexOffset = indexOffset;
        BbWriteRTree(output, regions.Select(r => new BbRegion { ChromId = r.ChromId, Start = r.Start, End = r.End, Offset = r.Offset }).ToList(), o.BlockSize, o.ItemsPerSlot, indexOffset);

        // zoomList.reduce
        long zoomCount = initialReduction.Size;
        long reduction = (long)initialReduction.Scale * 4;
        var summaries = rezoomedList;
        bool broke = false;
        for (int level = 1; level < 10; level++)
        {
            long rezoomCount = summaries.Count;
            if (rezoomCount >= zoomCount)
            {
                broke = true;
                break;
            }

            zoomCount = rezoomCount;
            var z = new BbZoom { DataOffset = output.Position };
            zoomList.Add(z);
            output.I32(zoomCount);
            int itemSize = o.Compress ? blockBytes : 32;
            foreach (var s in summaries)
                WriteItem(output, s, itemSize, o.Compress);
            Flush(output, o.Compress);
            long idxOffset = output.Position;
            BbWriteRTree(output, summaries.Select(r => new BbRegion { ChromId = r.ChromId, Start = r.Start, End = r.End, Offset = r.Offset }).ToList(), o.BlockSize, o.ItemsPerSlot, idxOffset);
            z.IndexOffset = idxOffset;
            z.ReductionLevel = reduction;
            reduction *= 4;
            int i = 0;
            long chromIdCur = -1;
            long end = 0;
            BbSummary? currentSummary = null;
            foreach (var summary in summaries.ToList())
            {
                if (summary.ChromId != chromIdCur || summary.End > end)
                {
                    end = summary.Start + reduction;
                    chromIdCur = summary.ChromId;
                    currentSummary = summary;
                    summaries[i] = currentSummary;
                    i++;
                }
                else
                {
                    currentSummary!.Add(summary);
                }
            }

            summaries.RemoveRange(i, summaries.Count - i);
        }

        if (!broke && zoomList.Count == 10)
            zoomList.RemoveAt(9);
        return (zoomList, totalSum);
    }

    // ------------------------------------------------------------------ trees

    private static void BbWriteBPlusTree(ByteSink output, List<(byte[] Name, long Id, long Size)> items, int keySize, int blockSize, int valSize)
    {
        var records = items.Select(it =>
        {
            var b = new byte[keySize + valSize];
            Array.Copy(it.Name, b, Math.Min(it.Name.Length, keySize));
            BitConverter.GetBytes((int)it.Id).CopyTo(b, keySize);
            BitConverter.GetBytes((int)it.Size).CopyTo(b, keySize + 4);
            return (Key: b[..keySize], Record: b);
        }).ToList();
        BbWriteBPlusTreeRecords(output, records, keySize, keySize + valSize, blockSize);
    }

    private static void BbWriteBPlusTreeChunks(ByteSink output, List<(byte[] Name, long Offset, long Size)> items, int keySize, int blockSize)
    {
        var records = items.Select(it =>
        {
            var b = new byte[keySize + 16];
            Array.Copy(it.Name, b, Math.Min(it.Name.Length, keySize));
            BitConverter.GetBytes((ulong)it.Offset).CopyTo(b, keySize);
            BitConverter.GetBytes((ulong)it.Size).CopyTo(b, keySize + 8);
            return (Key: b[..keySize], Record: b);
        }).ToList();
        BbWriteBPlusTreeRecords(output, records, keySize, keySize + 16, blockSize);
    }

    private static void BbWriteBPlusTreeRecords(ByteSink output, List<(byte[] Key, byte[] Record)> items, int keySize, int itemSize, int blockSize)
    {
        int valSize = itemSize - keySize;
        long itemCount = items.Count;
        output.U32(BPlusSignature);
        output.U32(blockSize);
        output.U32(keySize);
        output.U32(valSize);
        output.U64(itemCount);
        output.Zeros(8);
        int levels = 1;
        long n = itemCount;
        while (n > blockSize)
        {
            n = (n + blockSize - 1) / blockSize;
            levels++;
        }

        int nonLeafSize = keySize + 8;
        long bytesInIndexBlock = 4 + ((long)blockSize * nonLeafSize);
        long bytesInLeafBlock = 4 + ((long)blockSize * itemSize);
        long indexOffset = output.Position;
        for (int level = levels - 1; level > 0; level--)
        {
            long slotSizePer = (long)Math.Pow(blockSize, level);
            long nodeSizePer = slotSizePer * blockSize;
            long bytesInNextLevelBlock = level == 1 ? bytesInLeafBlock : bytesInIndexBlock;
            long indices = (itemCount + nodeSizePer - 1) / nodeSizePer;
            long levelSize = indices * bytesInIndexBlock;
            long endLevel = indexOffset + levelSize;
            long nextChild = endLevel;
            for (long index = 0; index < itemCount; index += nodeSizePer)
            {
                var block = new List<(byte[] Key, byte[] Record)>();
                for (long k = index; k < Math.Min(itemCount, index + nodeSizePer); k += slotSizePer)
                    block.Add(items[(int)k]);
                output.U8(0);
                output.U8(0);
                output.U16(block.Count);
                foreach (var item in block)
                {
                    output.Write(item.Key);
                    output.U64(nextChild);
                    nextChild += bytesInNextLevelBlock;
                }

                output.Zeros((long)(blockSize - block.Count) * nonLeafSize);
            }

            indexOffset = endLevel;
        }

        for (long index = 0; index < itemCount; index += blockSize)
        {
            var block = items.Skip((int)index).Take(blockSize).ToList();
            if (block.Count == 0)
                break;
            output.U8(1);
            output.U8(0);
            output.U16(block.Count);
            foreach (var item in block)
                output.Write(item.Record);
            output.Zeros((long)(blockSize - block.Count) * itemSize);
        }
    }

    private static void BbWriteRTree(ByteSink output, List<BbRegion> items, int blockSize, int itemsPerSlot, long endFileOffset)
    {
        var (root, levelCount) = BbBuildRTree(blockSize, items, endFileOffset);
        output.U32(RTreeSignature);
        output.U32(blockSize);
        output.U64(items.Count);
        output.U32(root.StartChromId);
        output.U32(root.StartBase);
        output.U32(root.EndChromId);
        output.U32(root.EndBase);
        output.U64(endFileOffset);
        output.U32(itemsPerSlot);
        output.Zeros(4);
        var levelSizes = new long[levelCount];
        root.CalcLevelSizes(levelSizes, 0);
        long size = 4 + (24L * blockSize);
        long levelOffset = output.Position;
        for (int i = 0; i < levelCount - 2; i++)
        {
            levelOffset += levelSizes[i] * size;
            if (i == levelCount - 3)
                size = 4 + (32L * blockSize);
            BbWriteIndexLevel(output, root, blockSize, size, 0, i, levelOffset);
        }

        int leafLevel = levelCount - 2;
        BbWriteLeaves(output, blockSize, size, root, 0, leafLevel);
    }

    private static (BbRNode Root, int Levels) BbBuildRTree(int blockSize, List<BbRegion> items, long endFileOffset)
    {
        int itemCount = items.Count;
        if (itemCount == 0)
            throw new ArgumentException("TypeError: cannot unpack non-iterable NoneType object");
        var children = new List<BbRNode>();
        long nextOffset = items[0].Offset;
        int i = 0;
        while (i < itemCount)
        {
            var child = new BbRNode();
            children.Add(child);
            var startItem = items[i];
            child.StartChromId = child.EndChromId = startItem.ChromId;
            child.StartBase = startItem.Start;
            child.EndBase = startItem.End;
            child.StartFileOffset = nextOffset;
            int oneSize = 1;
            bool broke = false;
            for (int j = i + 1; j < itemCount; j++)
            {
                var endItem = items[j];
                nextOffset = endItem.Offset;
                if (nextOffset != child.StartFileOffset)
                {
                    broke = true;
                    break;
                }

                oneSize++;
            }

            if (!broke)
                nextOffset = endFileOffset;
            child.EndFileOffset = nextOffset;
            for (int k = i + 1; k < Math.Min(itemCount, i + oneSize); k++)
            {
                var item = items[k];
                if (item.ChromId < child.StartChromId)
                {
                    child.StartChromId = item.ChromId;
                    child.StartBase = item.Start;
                }
                else if (item.ChromId == child.StartChromId && item.Start < child.StartBase)
                {
                    child.StartBase = item.Start;
                }

                if (item.ChromId > child.EndChromId)
                {
                    child.EndChromId = item.ChromId;
                    child.EndBase = item.End;
                }
                else if (item.ChromId == child.EndChromId && item.End > child.EndBase)
                {
                    child.EndBase = item.End;
                }
            }

            i += oneSize;
        }

        int levelCount = 1;
        BbRNode parent;
        while (true)
        {
            var parents = new List<BbRNode>();
            int slotsUsed = blockSize;
            parent = null!;
            foreach (var child in children)
            {
                if (slotsUsed >= blockSize)
                {
                    slotsUsed = 1;
                    parent = new BbRNode
                    {
                        Parent = child.Parent,
                        StartChromId = child.StartChromId,
                        StartBase = child.StartBase,
                        EndChromId = child.EndChromId,
                        EndBase = child.EndBase,
                        StartFileOffset = child.StartFileOffset,
                        EndFileOffset = child.EndFileOffset,
                    };
                    parents.Add(parent);
                }
                else
                {
                    slotsUsed++;
                    if (child.StartChromId < parent.StartChromId)
                    {
                        parent.StartChromId = child.StartChromId;
                        parent.StartBase = child.StartBase;
                    }
                    else if (child.StartChromId == parent.StartChromId && child.StartBase < parent.StartBase)
                    {
                        parent.StartBase = child.StartBase;
                    }

                    if (child.EndChromId > parent.EndChromId)
                    {
                        parent.EndChromId = child.EndChromId;
                        parent.EndBase = child.EndBase;
                    }
                    else if (child.EndChromId == parent.EndChromId && child.EndBase > parent.EndBase)
                    {
                        parent.EndBase = child.EndBase;
                    }
                }

                parent.Children.Add(child);
                child.Parent = parent;
            }

            levelCount++;
            if (parents.Count == 1)
                break;
            children = parents;
        }

        return (parent, levelCount);
    }

    private static void BbWriteLeaves(ByteSink output, int itemsPerSlot, long lNodeSize, BbRNode tree, int curLevel, int leafLevel)
    {
        if (curLevel == leafLevel)
        {
            output.U8(1);
            output.U8(0);
            output.U16(tree.Children.Count);
            foreach (var child in tree.Children)
            {
                output.U32(child.StartChromId);
                output.U32(child.StartBase);
                output.U32(child.EndChromId);
                output.U32(child.EndBase);
                output.U64(child.StartFileOffset);
                output.U64(child.EndFileOffset - child.StartFileOffset);
            }

            output.Zeros((long)(itemsPerSlot - tree.Children.Count) * 24);
        }
        else
        {
            foreach (var child in tree.Children)
                BbWriteLeaves(output, itemsPerSlot, lNodeSize, child, curLevel + 1, leafLevel);
        }
    }

    private static long BbWriteIndexLevel(ByteSink output, BbRNode parent, int blockSize, long childNodeSize, int curLevel, int destLevel, long offset)
    {
        long previousOffset = offset;
        if (curLevel == destLevel)
        {
            output.U8(0);
            output.U8(0);
            output.U16(parent.Children.Count);
            foreach (var child in parent.Children)
            {
                output.U32(child.StartChromId);
                output.U32(child.StartBase);
                output.U32(child.EndChromId);
                output.U32(child.EndBase);
                output.U64(offset);
                offset += childNodeSize;
            }

            output.Zeros((long)(blockSize - parent.Children.Count) * 24);
        }
        else
        {
            foreach (var child in parent.Children)
                offset = BbWriteIndexLevel(output, child, blockSize, childNodeSize, curLevel + 1, destLevel, offset);
        }

        if (output.Position != previousOffset)
            throw new ArgumentException($"Internal error: offset mismatch ({output.Position} vs {previousOffset})");
        return offset;
    }
}
