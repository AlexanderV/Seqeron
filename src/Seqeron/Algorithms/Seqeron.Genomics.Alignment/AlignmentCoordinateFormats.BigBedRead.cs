using System.IO.Compression;
using System.Text;
using System.Text.Json.Nodes;
using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Alignment;

public static partial class AlignmentCoordinateFormats
{
    /// <summary>A little / big endian reader over the bytes of a file with Python's <c>stream.read</c> semantics (a short read returns the bytes that exist).</summary>
    private sealed class ByteSource(byte[] data)
    {
        public long Pos;
        public bool BigEndian;

        public byte[] Read(long n)
        {
            long end = Math.Min(data.Length, Pos + n);
            long start = Math.Min(data.Length, Pos);
            var r = data.AsSpan((int)start, (int)(end - start)).ToArray();
            Pos += n;
            return r;
        }

        private byte[] Exact(int n)
        {
            var r = Read(n);
            if (r.Length != n)
                throw new ArgumentException($"unpack requires a buffer of {n} bytes");
            if (BigEndian)
                Array.Reverse(r);
            return r;
        }

        public int U8() => Exact(1)[0];

        public int U16() => BitConverter.ToUInt16(Exact(2));

        public long U32() => BitConverter.ToUInt32(Exact(4));

        public long U64() => (long)BitConverter.ToUInt64(Exact(8));

        public void Skip(int n) => Exact(n);
    }

    /// <summary>The fields of an autoSql declaration read from a bigBed file (<c>AutoSQLTable.from_bytes</c>).</summary>
    private static (string Name, string Comment, List<BigBedField> Fields) BbParseAutoSql(byte[] data)
    {
        if (data.Length == 0 || data[^1] != 0)
            throw new ArgumentException("AssertionError");
        string text = Encoding.UTF8.GetString(data, 0, data.Length - 1);
        var (word, rest) = SplitWs1(text);
        if (word != "table")
            throw new ArgumentException("AssertionError");
        var (name, rest2) = SplitWs1(rest);
        text = rest2;
        if (!text.StartsWith('"'))
            throw new ArgumentException("AssertionError");
        int i = text.IndexOf('"', 1);
        string comment = text[1..i];
        text = text[(i + 1)..].Trim();
        if (!text.StartsWith('(') || !text.EndsWith(')'))
            throw new ArgumentException("AssertionError");
        text = text[1..^1].Trim();
        var fields = new List<BigBedField>();
        while (text.Length > 0)
        {
            i = text.IndexOf('"');
            if (i < 0)
                throw new ArgumentException("ValueError: substring not found");
            int j = text.IndexOf('"', i + 1);
            if (j < 0)
                throw new ArgumentException("ValueError: substring not found");
            string fieldComment = text[(i + 1)..j];
            string definition = text[..i].Trim();
            if (!definition.EndsWith(';'))
                throw new ArgumentException("AssertionError");
            string inner = definition[..^1].TrimEnd();
            int sp = inner.LastIndexOfAny([' ', '\t', '\n', '\r']);
            if (sp < 0)
                throw new ArgumentException("not enough values to unpack");
            string fieldType = inner[..sp].TrimEnd(), fieldName = inner[(sp + 1)..];
            string dataType = fieldType.EndsWith(']') ? fieldType[..fieldType.IndexOf('[')] : fieldType;
            if (dataType is not ("int" or "uint" or "short" or "ushort" or "byte" or "ubyte" or "float" or "char" or "string" or "lstring"))
                throw new ArgumentException("AssertionError");
            fields.Add(new BigBedField(fieldType, fieldName, fieldComment));
            text = text[(j + 1)..].Trim();
        }

        return (name, comment, fields);
    }

    private sealed class BbReadNode
    {
        public BbReadNode? Parent;
        public List<BbReadNode> Children = [];
        public long DataOffset, DataSize;
        public bool IsLeaf;
    }

    /// <summary><c>_RTreeFormatter.read</c>: the tree of a file whose stream is positioned at the R-tree header.</summary>
    private static BbReadNode BbReadRTree(ByteSource s)
    {
        long magic = s.U32();
        s.U32();
        long itemCount = s.U64();
        s.Skip(16);
        s.U64();
        s.U32();
        s.Skip(4);
        if (magic != RTreeSignature)
            throw new ArgumentException("AssertionError");
        var root = new BbReadNode();
        var node = root;
        long itemsCounted = 0;
        while (true)
        {
            bool isLeaf = s.U8() != 0;
            s.Skip(1);
            int count = s.U16();
            if (isLeaf)
            {
                for (int i = 0; i < count; i++)
                {
                    s.Skip(16);
                    long dataOffset = s.U64();
                    long dataSize = s.U64();
                    node.Children.Add(new BbReadNode { Parent = node, IsLeaf = true, DataOffset = dataOffset, DataSize = dataSize });
                }

                itemsCounted += count;
                while (true)
                {
                    var parent = node.Parent;
                    if (parent is null)
                    {
                        if (itemsCounted != itemCount)
                            throw new ArgumentException("AssertionError");
                        return node;
                    }

                    int index = parent.Children.IndexOf(node);
                    if (index + 1 < parent.Children.Count)
                    {
                        node = parent.Children[index + 1];
                        break;
                    }

                    node = parent;
                }
            }
            else
            {
                for (int i = 0; i < count; i++)
                {
                    s.Skip(16);
                    long dataOffset = s.U64();
                    node.Children.Add(new BbReadNode { Parent = node, DataOffset = dataOffset });
                }

                node = node.Children[0];
            }

            s.Pos = node.DataOffset;
        }
    }

    private static List<BioCoordSequence> BbReadBPlusTree(ByteSource s)
    {
        long magic = s.U32();
        s.U32();
        long keySize = s.U32();
        long valSize = s.U32();
        long itemCount = s.U64();
        s.Skip(8);
        if (magic != BPlusSignature)
            throw new ArgumentException("AssertionError");
        if (valSize != 8)
            throw new ArgumentException("AssertionError");
        var targets = new List<BioCoordSequence>();
        var parents = new Stack<List<long>>();
        List<long>? current = null;
        while (true)
        {
            bool isLeaf = s.U8() != 0;
            s.Skip(1);
            int count = s.U16();
            if (isLeaf)
            {
                for (int i = 0; i < count; i++)
                {
                    var key = s.Read(keySize);
                    long chromId = s.U32();
                    long chromSize = s.U32();
                    int end = key.Length;
                    while (end > 0 && key[end - 1] == 0)
                        end--;
                    string name = Encoding.UTF8.GetString(key, 0, end);
                    if (chromId != targets.Count)
                        throw new ArgumentException("AssertionError");
                    targets.Add(new BioCoordSequence { Id = name, Length = chromSize, Parts = UndefinedParts(chromSize), Description = BiopythonUnknownDescription });
                }
            }
            else
            {
                var children = new List<long>();
                for (int i = 0; i < count; i++)
                {
                    s.Read(keySize);
                    children.Add(s.U64());
                }

                if (current is not null)
                    parents.Push(current);
                current = children;
            }

            long pos;
            while (true)
            {
                if (current is null)
                {
                    if (targets.Count != itemCount)
                        throw new ArgumentException("AssertionError");
                    return targets;
                }

                if (current.Count > 0)
                {
                    pos = current[0];
                    current.RemoveAt(0);
                    break;
                }

                current = parents.Count > 0 ? parents.Pop() : null;
            }

            s.Pos = pos;
        }
    }

    private static Func<string, object> BbConverter(string asType)
    {
        string t = asType;
        bool array = t.Contains('[') && t.Contains(']');
        if (array)
            t = t.Split('[')[0].Trim();
        Func<string, object> convert;
        if (t is "int" or "uint" or "short" or "ushort")
            convert = w => PyInt(w);
        else if (t is "byte" or "ubyte")
            convert = _ => throw new ArgumentException("TypeError: string argument without an encoding");
        else if (t == "float")
            convert = w => PyFloatParse(w);
        else if (t is "char" or "string" or "lstring")
            convert = w => w;
        else
            throw new ArgumentException($"Unknown field type {asType}");
        if (array)
            return w => w.TrimEnd(',').Split(',').Select(convert).ToList();
        return convert;
    }

    private static JsonNode BbToJson(object v) => v switch
    {
        long l => JsonValue.Create(l),
        double d => JsonValue.Create(d),
        string s => JsonValue.Create(s),
        List<object> list => new JsonArray(list.Select(BbToJson).ToArray()),
        _ => throw new ArgumentException("TypeError"),
    };

    /// <summary>Reads a bigBed file as <c>Bio.Align.parse(handle, "bigbed")</c> does.</summary>
    public static CoordinateAlignmentFile ParseBigBed(byte[] data)
    {
        ArgumentNullException.ThrowIfNull(data);
        try
        {
            return ParseBigBedCore(data);
        }
        catch (Exception e) when (e is IndexOutOfRangeException or ArgumentOutOfRangeException or KeyNotFoundException or InvalidCastException or NullReferenceException or FormatException or OverflowException or InvalidOperationException or InvalidDataException)
        {
            throw new ArgumentException(e.Message, e);
        }
    }

    private static CoordinateAlignmentFile ParseBigBedCore(byte[] data, bool maf = false)
    {
        var s = new ByteSource(data);
        var magic = s.Read(4);
        long little = magic.Length == 4 ? BitConverter.ToUInt32(magic) : -1;
        long big = magic.Length == 4 ? System.Buffers.Binary.BinaryPrimitives.ReadUInt32BigEndian(magic) : -1;
        if (little == BigBedSignature)
            s.BigEndian = false;
        else if (big == BigBedSignature)
            s.BigEndian = true;
        else
            throw new ArgumentException("not a bigBed file");
        int version = s.U16();
        int zoomLevels = s.U16();
        long chromosomeTreeOffset = s.U64();
        long fullDataOffset = s.U64();
        long fullIndexOffset = s.U64();
        int fieldCount = s.U16();
        int definedFieldCount = s.U16();
        long autoSqlOffset = s.U64();
        long totalSummaryOffset = s.U64();
        long uncompressBufSize = s.U32();
        s.U64();
        if (version != 4)
            throw new ArgumentException("AssertionError");
        _ = zoomLevels;
        if (definedFieldCount < 3 || definedFieldCount > 12)
            throw new ArgumentException($"expected between 3 and 12 columns, found {definedFieldCount}");
        var zoomIndexOffsets = new List<long>();
        for (int z = 0; z < 10; z++)
        {
            long reductionLevel = s.U32();
            s.Skip(4);
            s.U64();
            long indexOffset = s.U64();
            if (reductionLevel == 0)
                break;
            zoomIndexOffsets.Add(indexOffset);
        }

        foreach (long off in zoomIndexOffsets)
        {
            s.Pos = off;
            long sig = s.U32();
            s.U32();
            s.U64();
            s.Skip(16);
            s.U64();
            s.U32();
            if (sig != RTreeSignature)
                throw new ArgumentException("AssertionError");
        }

        // autoSql
        long autoSqlSize = totalSummaryOffset - autoSqlOffset;
        s.Pos = autoSqlOffset;
        var (declName, declComment, declFields) = BbParseAutoSql(s.Read(autoSqlSize));
        int bedN = definedFieldCount;
        string[] names = ["chrom", "chromStart", "chromEnd", "name", "score", "strand", "thickStart", "thickEnd", "reserved", "blockCount", "blockSizes", "chromStarts"];
        for (int i = 0; i < bedN; i++)
        {
            string fname = declFields[i].Name;
            if (fname != names[i])
                throw new ArgumentException($"Expected field name '{names[i]}'; found '{fname}'");
        }

        var customFields = new List<(string Name, string AsType)>();
        for (int i = definedFieldCount; i < fieldCount; i++)
            customFields.Add((declFields[i].Name, declFields[i].AsType));
        var converters = customFields.Select(f => BbConverter(f.AsType)).ToList();

        s.Pos = fullDataOffset;
        _ = s.U64();
        bool compressed = uncompressBufSize > 0;
        s.Pos = chromosomeTreeOffset;
        var targets = BbReadBPlusTree(s);
        s.Pos = fullIndexOffset;
        var tree = BbReadRTree(s);

        var alignments = new List<BioCoordAlignment>();
        // depth-first walk over the leaves
        var leaves = new List<BbReadNode>();
        void Walk(BbReadNode node)
        {
            if (node.IsLeaf)
            {
                leaves.Add(node);
                return;
            }

            foreach (var c in node.Children)
                Walk(c);
        }

        Walk(tree);
        string? mafReference = null;
        foreach (var leaf in leaves)
        {
            s.Pos = leaf.DataOffset;
            byte[] block = s.Read(leaf.DataSize);
            if (compressed)
            {
                using var input = new MemoryStream(block);
                using var z = new ZLibStream(input, CompressionMode.Decompress);
                using var output = new MemoryStream();
                z.CopyTo(output);
                block = output.ToArray();
            }

            if (maf && mafReference is null)
            {
                mafReference = BigMafReference(block);
                targets[0] = new BioCoordSequence { Id = mafReference + "." + targets[0].Id, Length = targets[0].Length, Parts = targets[0].Parts, Description = targets[0].Description, Annotations = targets[0].Annotations };
            }

            int pos = 0;
            while (pos < block.Length)
            {
                if (pos + 12 > block.Length)
                    throw new ArgumentException("unpack requires a buffer of 12 bytes");
                var bs = new ByteSource(block[pos..(pos + 12)]) { BigEndian = s.BigEndian };
                long chromId = bs.U32(), chromStart = bs.U32(), chromEnd = bs.U32();
                int zero = Array.IndexOf(block, (byte)0, pos + 12);
                if (zero < 0)
                    throw new ArgumentException("ValueError: subsection not found");
                byte[] rest = block[(pos + 12)..zero];
                pos = zero + 1;
                alignments.Add(maf ? BbCreateMafAlignment(rest) : BbCreateAlignment(chromId, chromStart, chromEnd, rest, targets, bedN, customFields, converters));
            }
        }

        var decl = new BigBedDeclaration(declName, declComment, declFields);
        return new CoordinateAlignmentFile { Alignments = alignments, Metadata = [], Targets = targets, Declaration = decl };
    }

    private static BioCoordAlignment BbCreateAlignment(long chromId, long chromStart, long chromEnd, byte[] restBytes, List<BioCoordSequence> targets, int bedN, List<(string Name, string AsType)> custom, List<Func<string, object>> converters)
    {
        string text = Encoding.UTF8.GetString(restBytes);
        var words = text.Length > 0 ? text.Split('\t') : [];
        var targetRecord = targets[(int)chromId];
        string? name = bedN > 3 ? words[0] : null;
        string strand = bedN > 5 ? words[2] : "+";
        long[][] coordinates;
        long qSize;
        if (bedN > 9)
        {
            long blockCount = PyInt(words[6]);
            var blockSizes = words[7].TrimEnd(',').Split(',').Select(PyInt).ToArray();
            var blockStarts = words[8].TrimEnd(',').Split(',').Select(PyInt).ToArray();
            if (blockSizes.Length != blockCount)
                throw new ArgumentException($"Inconsistent number of block sizes ({blockSizes.Length} found, expected {blockCount})");
            if (blockStarts.Length != blockCount)
                throw new ArgumentException($"Inconsistent number of block start positions ({blockStarts.Length} found, expected {blockCount})");
            long tPosition = 0, qPosition = 0;
            var t = new List<long> { tPosition };
            var q = new List<long> { qPosition };
            for (int i = 0; i < blockSizes.Length; i++)
            {
                if (blockStarts[i] != tPosition)
                {
                    t.Add(blockStarts[i]);
                    q.Add(qPosition);
                    tPosition = blockStarts[i];
                }

                tPosition += blockSizes[i];
                qPosition += blockSizes[i];
                t.Add(tPosition);
                q.Add(qPosition);
            }

            coordinates = [[.. t], [.. q]];
            qSize = blockSizes.Sum();
        }
        else
        {
            long blockSize = chromEnd - chromStart;
            coordinates = [[0, blockSize], [0, blockSize]];
            qSize = blockSize;
        }

        for (int i = 0; i < coordinates[0].Length; i++)
            coordinates[0][i] += chromStart;
        var querySeq = new BioCoordSequence { Id = name, Length = qSize, Parts = UndefinedParts(qSize), Description = BiopythonUnknownDescription };
        if (strand == "-")
        {
            for (int i = 0; i < coordinates[1].Length; i++)
                coordinates[1][i] = qSize - coordinates[1][i];
        }

        if (chromStart != coordinates[0][0])
            throw new ArgumentException($"Inconsistent chromStart found ({chromStart}, expected {coordinates[0][0]})");
        if (chromEnd != coordinates[0][^1])
            throw new ArgumentException($"Inconsistent chromEnd found ({chromEnd}, expected {coordinates[0][^1]})");
        var attrs = new JsonObject();
        if (words.Length > bedN - 3)
        {
            var ann = new JsonObject();
            for (int i = 0; i < Math.Min(words.Length - (bedN - 3), custom.Count); i++)
                ann[custom[i].Name] = BbToJson(converters[i](words[(bedN - 3) + i]));
            attrs["annotations"] = ann;
        }

        var tgt = new BioCoordSequence { Id = targetRecord.Id, Length = targetRecord.Length, Parts = targetRecord.Parts, Description = BiopythonUnknownDescription };
        var alignment = new BioCoordAlignment { Sequences = [tgt, querySeq], Coordinates = coordinates, Attributes = attrs };
        if (bedN <= 4)
            return alignment;
        string scoreText = words[1];
        if (double.TryParse(scoreText, System.Globalization.NumberStyles.Float, System.Globalization.CultureInfo.InvariantCulture, out double score) && !scoreText.Contains('_'))
        {
            if (score == Math.Floor(score) && !double.IsInfinity(score))
                attrs["score"] = (long)score;
            else
                attrs["score"] = score;
        }
        else
        {
            attrs["score"] = scoreText;
        }

        if (bedN <= 6)
            return alignment;
        attrs["thickStart"] = PyInt(words[3]);
        if (bedN <= 7)
            return alignment;
        attrs["thickEnd"] = PyInt(words[4]);
        if (bedN <= 8)
            return alignment;
        attrs["itemRgb"] = words[5];
        return alignment;
    }
}
