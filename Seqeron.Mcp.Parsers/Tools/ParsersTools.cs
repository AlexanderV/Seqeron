using System.ComponentModel;
using ModelContextProtocol.Server;
using Seqeron.Genomics;

namespace Seqeron.Mcp.Parsers.Tools;

[McpServerToolType]
public static class ParsersTools
{
    // ========================
    // FASTA Tools
    // ========================

    [McpServerTool(Name = "fasta_parse")]
    [Description("Parse FASTA format string into sequence entries. Returns list of sequences with their IDs and descriptions.")]
    public static FastaParseResult FastaParse(
        [Description("FASTA format content to parse")] string content)
    {
        if (string.IsNullOrEmpty(content))
            throw new ArgumentException("Content cannot be null or empty", nameof(content));

        var entries = FastaParser.Parse(content).ToList();
        var results = entries.Select(e => new FastaEntryResult(
            e.Id,
            e.Description,
            e.Sequence.Sequence,
            e.Sequence.Length
        )).ToList();

        return new FastaParseResult(results, results.Count);
    }

    [McpServerTool(Name = "fasta_format")]
    [Description("Format sequence(s) to FASTA string. Accepts ID, optional description, and sequence.")]
    public static FastaFormatResult FastaFormat(
        [Description("Sequence identifier")] string id,
        [Description("DNA sequence")] string sequence,
        [Description("Optional sequence description")] string? description = null,
        [Description("Line width for sequence wrapping (default: 80)")] int lineWidth = 80)
    {
        if (string.IsNullOrEmpty(id))
            throw new ArgumentException("ID cannot be null or empty", nameof(id));
        if (string.IsNullOrEmpty(sequence))
            throw new ArgumentException("Sequence cannot be null or empty", nameof(sequence));
        if (lineWidth < 1)
            throw new ArgumentException("Line width must be at least 1", nameof(lineWidth));

        var dnaSeq = new DnaSequence(sequence);
        var entry = new FastaEntry(id, description, dnaSeq);
        var fasta = FastaParser.ToFasta(new[] { entry }, lineWidth);

        return new FastaFormatResult(fasta.TrimEnd());
    }

    [McpServerTool(Name = "fasta_write")]
    [Description("Write sequence(s) to a FASTA file. Creates or overwrites the file at specified path.")]
    public static FastaWriteResult FastaWrite(
        [Description("File path to write FASTA output")] string filePath,
        [Description("Sequence identifier")] string id,
        [Description("DNA sequence")] string sequence,
        [Description("Optional sequence description")] string? description = null,
        [Description("Line width for sequence wrapping (default: 80)")] int lineWidth = 80)
    {
        if (string.IsNullOrEmpty(filePath))
            throw new ArgumentException("File path cannot be null or empty", nameof(filePath));
        if (string.IsNullOrEmpty(id))
            throw new ArgumentException("ID cannot be null or empty", nameof(id));
        if (string.IsNullOrEmpty(sequence))
            throw new ArgumentException("Sequence cannot be null or empty", nameof(sequence));
        if (lineWidth < 1)
            throw new ArgumentException("Line width must be at least 1", nameof(lineWidth));

        var dnaSeq = new DnaSequence(sequence);
        var entry = new FastaEntry(id, description, dnaSeq);
        FastaParser.WriteFile(filePath, new[] { entry }, lineWidth);

        return new FastaWriteResult(filePath, 1, sequence.Length);
    }

    // ========================
    // FASTQ Tools
    // ========================

    [McpServerTool(Name = "fastq_parse")]
    [Description("Parse FASTQ format string into sequence entries with quality scores. Supports Phred+33 and Phred+64 encodings.")]
    public static FastqParseResult FastqParse(
        [Description("FASTQ format content to parse")] string content,
        [Description("Quality encoding: 'phred33' (Sanger/Illumina 1.8+), 'phred64' (Illumina 1.3-1.7), or 'auto' (default)")] string encoding = "auto")
    {
        if (string.IsNullOrEmpty(content))
            throw new ArgumentException("Content cannot be null or empty", nameof(content));

        var qualityEncoding = encoding.ToLowerInvariant() switch
        {
            "phred33" => FastqParser.QualityEncoding.Phred33,
            "phred64" => FastqParser.QualityEncoding.Phred64,
            "auto" => FastqParser.QualityEncoding.Auto,
            _ => throw new ArgumentException($"Invalid encoding: {encoding}. Use 'phred33', 'phred64', or 'auto'", nameof(encoding))
        };

        var records = FastqParser.Parse(content, qualityEncoding).ToList();
        var results = records.Select(r => new FastqRecordResult(
            r.Id,
            r.Description,
            r.Sequence,
            r.QualityString,
            r.QualityScores.ToList(),
            r.Sequence.Length
        )).ToList();

        return new FastqParseResult(results, results.Count);
    }

    [McpServerTool(Name = "fastq_statistics")]
    [Description("Calculate quality statistics for FASTQ data. Returns read counts, quality metrics, and GC content.")]
    public static FastqStatisticsResult FastqStatistics(
        [Description("FASTQ format content to analyze")] string content,
        [Description("Quality encoding: 'phred33' (Sanger/Illumina 1.8+), 'phred64' (Illumina 1.3-1.7), or 'auto' (default)")] string encoding = "auto")
    {
        if (string.IsNullOrEmpty(content))
            throw new ArgumentException("Content cannot be null or empty", nameof(content));

        var qualityEncoding = encoding.ToLowerInvariant() switch
        {
            "phred33" => FastqParser.QualityEncoding.Phred33,
            "phred64" => FastqParser.QualityEncoding.Phred64,
            "auto" => FastqParser.QualityEncoding.Auto,
            _ => throw new ArgumentException($"Invalid encoding: {encoding}. Use 'phred33', 'phred64', or 'auto'", nameof(encoding))
        };

        var records = FastqParser.Parse(content, qualityEncoding).ToList();
        var stats = FastqParser.CalculateStatistics(records);

        return new FastqStatisticsResult(
            stats.TotalReads,
            stats.TotalBases,
            stats.MeanReadLength,
            stats.MeanQuality,
            stats.MinReadLength,
            stats.MaxReadLength,
            stats.Q20Percentage,
            stats.Q30Percentage,
            stats.GcContent);
    }

    [McpServerTool(Name = "fastq_filter")]
    [Description("Filter FASTQ reads by minimum average quality score. Returns reads meeting the quality threshold.")]
    public static FastqFilterResult FastqFilter(
        [Description("FASTQ format content to filter")] string content,
        [Description("Minimum average quality score threshold (e.g., 20 or 30)")] double minQuality,
        [Description("Quality encoding: 'phred33' (Sanger/Illumina 1.8+), 'phred64' (Illumina 1.3-1.7), or 'auto' (default)")] string encoding = "auto")
    {
        if (string.IsNullOrEmpty(content))
            throw new ArgumentException("Content cannot be null or empty", nameof(content));
        if (minQuality < 0)
            throw new ArgumentException("Minimum quality must be non-negative", nameof(minQuality));

        var qualityEncoding = encoding.ToLowerInvariant() switch
        {
            "phred33" => FastqParser.QualityEncoding.Phred33,
            "phred64" => FastqParser.QualityEncoding.Phred64,
            "auto" => FastqParser.QualityEncoding.Auto,
            _ => throw new ArgumentException($"Invalid encoding: {encoding}. Use 'phred33', 'phred64', or 'auto'", nameof(encoding))
        };

        var records = FastqParser.Parse(content, qualityEncoding).ToList();
        var filtered = FastqParser.FilterByQuality(records, minQuality).ToList();

        var results = filtered.Select(r => new FastqRecordResult(
            r.Id,
            r.Description,
            r.Sequence,
            r.QualityString,
            r.QualityScores.ToList(),
            r.Sequence.Length
        )).ToList();

        return new FastqFilterResult(
            results,
            results.Count,
            records.Count,
            records.Count > 0 ? (double)results.Count / records.Count * 100 : 0);
    }

    [McpServerTool(Name = "fastq_detect_encoding")]
    [Description("Detect quality score encoding from FASTQ quality string. Returns 'Phred33' (Sanger/Illumina 1.8+) or 'Phred64' (Illumina 1.3-1.7).")]
    public static FastqDetectEncodingResult FastqDetectEncoding(
        [Description("Quality string from FASTQ record")] string qualityString)
    {
        if (string.IsNullOrEmpty(qualityString))
            throw new ArgumentException("Quality string cannot be null or empty", nameof(qualityString));

        var encoding = FastqParser.DetectEncoding(qualityString);
        return new FastqDetectEncodingResult(
            encoding.ToString(),
            encoding == FastqParser.QualityEncoding.Phred33 ? 33 : 64);
    }

    [McpServerTool(Name = "fastq_encode_quality")]
    [Description("Encode Phred quality scores to a quality string.")]
    public static FastqEncodeQualityResult FastqEncodeQuality(
        [Description("List of Phred quality scores (0-41)")] List<int> scores,
        [Description("Quality encoding: 'phred33' (default) or 'phred64'")] string encoding = "phred33")
    {
        if (scores == null || scores.Count == 0)
            throw new ArgumentException("Scores cannot be null or empty", nameof(scores));

        var qualityEncoding = encoding.ToLowerInvariant() switch
        {
            "phred33" => FastqParser.QualityEncoding.Phred33,
            "phred64" => FastqParser.QualityEncoding.Phred64,
            _ => throw new ArgumentException($"Invalid encoding: {encoding}. Use 'phred33' or 'phred64'", nameof(encoding))
        };

        var qualityString = FastqParser.EncodeQualityScores(scores, qualityEncoding);
        return new FastqEncodeQualityResult(qualityString, scores.Count);
    }

    [McpServerTool(Name = "fastq_phred_to_error")]
    [Description("Convert Phred quality score to error probability. Formula: P = 10^(-Q/10).")]
    public static FastqPhredToErrorResult FastqPhredToError(
        [Description("Phred quality score")] int phredScore)
    {
        if (phredScore < 0)
            throw new ArgumentException("Phred score cannot be negative", nameof(phredScore));

        var errorProbability = FastqParser.PhredToErrorProbability(phredScore);
        return new FastqPhredToErrorResult(phredScore, errorProbability);
    }

    [McpServerTool(Name = "fastq_error_to_phred")]
    [Description("Convert error probability to Phred quality score. Formula: Q = -10 * log10(P).")]
    public static FastqErrorToPhredResult FastqErrorToPhred(
        [Description("Error probability (0-1)")] double errorProbability)
    {
        if (errorProbability < 0 || errorProbability > 1)
            throw new ArgumentException("Error probability must be between 0 and 1", nameof(errorProbability));

        var phredScore = FastqParser.ErrorProbabilityToPhred(errorProbability);
        return new FastqErrorToPhredResult(errorProbability, phredScore);
    }

    [McpServerTool(Name = "fastq_trim_quality")]
    [Description("Trim low-quality bases from both ends of FASTQ reads. Removes bases below the quality threshold.")]
    public static FastqTrimQualityResult FastqTrimQuality(
        [Description("FASTQ format content to trim")] string content,
        [Description("Minimum quality score threshold (default: 20)")] int minQuality = 20,
        [Description("Quality encoding: 'phred33' (default), 'phred64', or 'auto'")] string encoding = "auto")
    {
        if (string.IsNullOrEmpty(content))
            throw new ArgumentException("Content cannot be null or empty", nameof(content));
        if (minQuality < 0)
            throw new ArgumentException("Minimum quality must be non-negative", nameof(minQuality));

        var qualityEncoding = encoding.ToLowerInvariant() switch
        {
            "phred33" => FastqParser.QualityEncoding.Phred33,
            "phred64" => FastqParser.QualityEncoding.Phred64,
            "auto" => FastqParser.QualityEncoding.Auto,
            _ => throw new ArgumentException($"Invalid encoding: {encoding}. Use 'phred33', 'phred64', or 'auto'", nameof(encoding))
        };

        var records = FastqParser.Parse(content, qualityEncoding).ToList();
        var trimmed = records.Select(r => FastqParser.TrimByQuality(r, minQuality)).ToList();

        var results = trimmed.Select(r => new FastqRecordResult(
            r.Id,
            r.Description,
            r.Sequence,
            r.QualityString,
            r.QualityScores.ToList(),
            r.Sequence.Length
        )).ToList();

        var originalBases = records.Sum(r => r.Sequence.Length);
        var trimmedBases = trimmed.Sum(r => r.Sequence.Length);

        return new FastqTrimQualityResult(
            results,
            results.Count,
            originalBases,
            trimmedBases,
            originalBases > 0 ? (double)(originalBases - trimmedBases) / originalBases * 100 : 0);
    }

    [McpServerTool(Name = "fastq_trim_adapter")]
    [Description("Trim adapter sequences from FASTQ reads. Searches for adapter at the 3' end and within the sequence.")]
    public static FastqTrimAdapterResult FastqTrimAdapter(
        [Description("FASTQ format content to trim")] string content,
        [Description("Adapter sequence to remove")] string adapter,
        [Description("Minimum overlap length to consider a match (default: 5)")] int minOverlap = 5,
        [Description("Quality encoding: 'phred33' (default), 'phred64', or 'auto'")] string encoding = "auto")
    {
        if (string.IsNullOrEmpty(content))
            throw new ArgumentException("Content cannot be null or empty", nameof(content));
        if (string.IsNullOrEmpty(adapter))
            throw new ArgumentException("Adapter cannot be null or empty", nameof(adapter));
        if (minOverlap < 1)
            throw new ArgumentException("Minimum overlap must be at least 1", nameof(minOverlap));

        var qualityEncoding = encoding.ToLowerInvariant() switch
        {
            "phred33" => FastqParser.QualityEncoding.Phred33,
            "phred64" => FastqParser.QualityEncoding.Phred64,
            "auto" => FastqParser.QualityEncoding.Auto,
            _ => throw new ArgumentException($"Invalid encoding: {encoding}. Use 'phred33', 'phred64', or 'auto'", nameof(encoding))
        };

        var records = FastqParser.Parse(content, qualityEncoding).ToList();
        var trimmed = records.Select(r => FastqParser.TrimAdapter(r, adapter, minOverlap)).ToList();

        var results = trimmed.Select(r => new FastqRecordResult(
            r.Id,
            r.Description,
            r.Sequence,
            r.QualityString,
            r.QualityScores.ToList(),
            r.Sequence.Length
        )).ToList();

        var readsWithAdapter = records.Zip(trimmed, (original, trim) => original.Sequence.Length != trim.Sequence.Length).Count(x => x);
        var originalBases = records.Sum(r => r.Sequence.Length);
        var trimmedBases = trimmed.Sum(r => r.Sequence.Length);

        return new FastqTrimAdapterResult(
            results,
            results.Count,
            readsWithAdapter,
            originalBases,
            trimmedBases);
    }

    [McpServerTool(Name = "fastq_format")]
    [Description("Format a single FASTQ record to string format.")]
    public static FastqFormatResult FastqFormat(
        [Description("Sequence identifier")] string id,
        [Description("DNA sequence")] string sequence,
        [Description("Quality string (same length as sequence)")] string qualityString,
        [Description("Optional sequence description")] string? description = null)
    {
        if (string.IsNullOrEmpty(id))
            throw new ArgumentException("ID cannot be null or empty", nameof(id));
        if (string.IsNullOrEmpty(sequence))
            throw new ArgumentException("Sequence cannot be null or empty", nameof(sequence));
        if (string.IsNullOrEmpty(qualityString))
            throw new ArgumentException("Quality string cannot be null or empty", nameof(qualityString));
        if (sequence.Length != qualityString.Length)
            throw new ArgumentException("Sequence and quality string must have the same length");

        var record = new FastqParser.FastqRecord(
            id,
            description ?? "",
            sequence,
            qualityString,
            FastqParser.DecodeQualityScores(qualityString));

        var fastq = FastqParser.ToFastqString(record);
        return new FastqFormatResult(fastq.TrimEnd());
    }

    [McpServerTool(Name = "fastq_write")]
    [Description("Write FASTQ records to a file. Creates or overwrites the file at specified path.")]
    public static FastqWriteResult FastqWrite(
        [Description("File path to write FASTQ output")] string filePath,
        [Description("FASTQ format content to write")] string content,
        [Description("Quality encoding for parsing: 'phred33' (default), 'phred64', or 'auto'")] string encoding = "auto")
    {
        if (string.IsNullOrEmpty(filePath))
            throw new ArgumentException("File path cannot be null or empty", nameof(filePath));
        if (string.IsNullOrEmpty(content))
            throw new ArgumentException("Content cannot be null or empty", nameof(content));

        var qualityEncoding = encoding.ToLowerInvariant() switch
        {
            "phred33" => FastqParser.QualityEncoding.Phred33,
            "phred64" => FastqParser.QualityEncoding.Phred64,
            "auto" => FastqParser.QualityEncoding.Auto,
            _ => throw new ArgumentException($"Invalid encoding: {encoding}. Use 'phred33', 'phred64', or 'auto'", nameof(encoding))
        };

        var records = FastqParser.Parse(content, qualityEncoding).ToList();
        FastqParser.WriteToFile(filePath, records);

        var totalBases = records.Sum(r => r.Sequence.Length);
        return new FastqWriteResult(filePath, records.Count, totalBases);
    }

    // ========================
    // BED Tools
    // ========================

    [McpServerTool(Name = "bed_parse")]
    [Description("Parse BED format content into genomic region records. Supports BED3, BED6, and BED12 formats.")]
    public static BedParseResult BedParse(
        [Description("BED format content to parse")] string content,
        [Description("BED format: 'bed3', 'bed6', 'bed12', or 'auto' (default)")] string format = "auto")
    {
        if (string.IsNullOrEmpty(content))
            throw new ArgumentException("Content cannot be null or empty", nameof(content));

        var bedFormat = format.ToLowerInvariant() switch
        {
            "bed3" => BedParser.BedFormat.BED3,
            "bed6" => BedParser.BedFormat.BED6,
            "bed12" => BedParser.BedFormat.BED12,
            "auto" => BedParser.BedFormat.Auto,
            _ => throw new ArgumentException($"Invalid format: {format}. Use 'bed3', 'bed6', 'bed12', or 'auto'", nameof(format))
        };

        var records = BedParser.Parse(content, bedFormat).ToList();
        var results = records.Select(r => new BedRecordResult(
            r.Chrom,
            r.ChromStart,
            r.ChromEnd,
            r.Length,
            r.Name,
            r.Score,
            r.Strand?.ToString(),
            r.ThickStart,
            r.ThickEnd,
            r.ItemRgb,
            r.BlockCount,
            r.BlockSizes?.ToList(),
            r.BlockStarts?.ToList()
        )).ToList();

        return new BedParseResult(results, results.Count);
    }

    [McpServerTool(Name = "bed_filter")]
    [Description("Filter BED records by chromosome, region, strand, length, or score. All filters are optional and can be combined.")]
    public static BedFilterResult BedFilter(
        [Description("BED format content to filter")] string content,
        [Description("Filter by chromosome name (e.g., 'chr1')")] string? chrom = null,
        [Description("Filter by region start position (requires chrom and regionEnd)")] int? regionStart = null,
        [Description("Filter by region end position (requires chrom and regionStart)")] int? regionEnd = null,
        [Description("Filter by strand: '+' or '-'")] string? strand = null,
        [Description("Filter by minimum feature length")] int? minLength = null,
        [Description("Filter by maximum feature length")] int? maxLength = null,
        [Description("Filter by minimum score")] int? minScore = null,
        [Description("Filter by maximum score")] int? maxScore = null)
    {
        if (string.IsNullOrEmpty(content))
            throw new ArgumentException("Content cannot be null or empty", nameof(content));

        var records = BedParser.Parse(content).ToList();
        IEnumerable<BedParser.BedRecord> filtered = records;

        // Apply chromosome filter
        if (!string.IsNullOrEmpty(chrom))
        {
            filtered = BedParser.FilterByChrom(filtered, chrom);
        }

        // Apply region filter (requires both start and end)
        if (regionStart.HasValue && regionEnd.HasValue && !string.IsNullOrEmpty(chrom))
        {
            filtered = BedParser.FilterByRegion(filtered, chrom, regionStart.Value, regionEnd.Value);
        }

        // Apply strand filter
        if (!string.IsNullOrEmpty(strand))
        {
            if (strand != "+" && strand != "-")
                throw new ArgumentException("Strand must be '+' or '-'", nameof(strand));
            filtered = BedParser.FilterByStrand(filtered, strand[0]);
        }

        // Apply length filter
        if (minLength.HasValue || maxLength.HasValue)
        {
            filtered = BedParser.FilterByLength(filtered, minLength ?? 0, maxLength);
        }

        // Apply score filter
        if (minScore.HasValue || maxScore.HasValue)
        {
            filtered = BedParser.FilterByScore(filtered, minScore ?? 0, maxScore);
        }

        var filteredList = filtered.ToList();
        var results = filteredList.Select(r => new BedRecordResult(
            r.Chrom,
            r.ChromStart,
            r.ChromEnd,
            r.Length,
            r.Name,
            r.Score,
            r.Strand?.ToString(),
            r.ThickStart,
            r.ThickEnd,
            r.ItemRgb,
            r.BlockCount,
            r.BlockSizes?.ToList(),
            r.BlockStarts?.ToList()
        )).ToList();

        return new BedFilterResult(
            results,
            results.Count,
            records.Count,
            records.Count > 0 ? (double)results.Count / records.Count * 100 : 0);
    }

    [McpServerTool(Name = "bed_merge")]
    [Description("Merge overlapping BED records into single intervals. Adjacent or overlapping features on the same chromosome are combined.")]
    public static BedMergeResult BedMerge(
        [Description("BED format content to merge")] string content)
    {
        if (string.IsNullOrEmpty(content))
            throw new ArgumentException("Content cannot be null or empty", nameof(content));

        var records = BedParser.Parse(content).ToList();
        var merged = BedParser.MergeOverlapping(records).ToList();

        var results = merged.Select(r => new BedRecordResult(
            r.Chrom,
            r.ChromStart,
            r.ChromEnd,
            r.Length,
            r.Name,
            r.Score,
            r.Strand?.ToString(),
            r.ThickStart,
            r.ThickEnd,
            r.ItemRgb,
            r.BlockCount,
            r.BlockSizes?.ToList(),
            r.BlockStarts?.ToList()
        )).ToList();

        return new BedMergeResult(results, results.Count, records.Count);
    }

    [McpServerTool(Name = "bed_intersect")]
    [Description("Find intersecting regions between two BED datasets. Returns overlapping portions of features.")]
    public static BedIntersectResult BedIntersect(
        [Description("First BED format content (features to intersect)")] string contentA,
        [Description("Second BED format content (reference features)")] string contentB)
    {
        if (string.IsNullOrEmpty(contentA))
            throw new ArgumentException("Content A cannot be null or empty", nameof(contentA));
        if (string.IsNullOrEmpty(contentB))
            throw new ArgumentException("Content B cannot be null or empty", nameof(contentB));

        var recordsA = BedParser.Parse(contentA).ToList();
        var recordsB = BedParser.Parse(contentB).ToList();
        var intersected = BedParser.Intersect(recordsA, recordsB).ToList();

        var results = intersected.Select(r => new BedRecordResult(
            r.Chrom,
            r.ChromStart,
            r.ChromEnd,
            r.Length,
            r.Name,
            r.Score,
            r.Strand?.ToString(),
            r.ThickStart,
            r.ThickEnd,
            r.ItemRgb,
            r.BlockCount,
            r.BlockSizes?.ToList(),
            r.BlockStarts?.ToList()
        )).ToList();

        return new BedIntersectResult(results, results.Count, recordsA.Count, recordsB.Count);
    }

    // ========================
    // VCF Tools
    // ========================

    [McpServerTool(Name = "vcf_parse")]
    [Description("Parse VCF (Variant Call Format) content into variant records. Returns chromosome, position, reference/alternate alleles, quality, and filter status.")]
    public static VcfParseResult VcfParse(
        [Description("VCF format content to parse")] string content)
    {
        if (string.IsNullOrEmpty(content))
            throw new ArgumentException("Content cannot be null or empty", nameof(content));

        var records = VcfParser.Parse(content).ToList();
        var results = records.Select(r => new VcfRecordResult(
            r.Chrom,
            r.Pos,
            r.Id,
            r.Ref,
            r.Alt.ToList(),
            r.Qual,
            r.Filter.ToList(),
            r.Info.ToDictionary(kv => kv.Key, kv => kv.Value),
            VcfParser.ClassifyVariant(r).ToString()
        )).ToList();

        return new VcfParseResult(results, results.Count);
    }

    [McpServerTool(Name = "vcf_statistics")]
    [Description("Calculate statistics for VCF variants. Returns counts by variant type, chromosome distribution, and quality metrics.")]
    public static VcfStatisticsResult VcfStatistics(
        [Description("VCF format content to analyze")] string content)
    {
        if (string.IsNullOrEmpty(content))
            throw new ArgumentException("Content cannot be null or empty", nameof(content));

        var records = VcfParser.Parse(content).ToList();
        var stats = VcfParser.CalculateStatistics(records);

        return new VcfStatisticsResult(
            stats.TotalVariants,
            stats.SnpCount,
            stats.IndelCount,
            stats.ComplexCount,
            stats.PassingCount,
            stats.ChromosomeCounts.ToDictionary(kv => kv.Key, kv => kv.Value),
            stats.MeanQuality);
    }

    [McpServerTool(Name = "vcf_filter")]
    [Description("Filter VCF variants by type, quality, chromosome, or PASS status.")]
    public static VcfFilterResult VcfFilter(
        [Description("VCF format content to filter")] string content,
        [Description("Filter by variant type: 'snp', 'indel', 'insertion', 'deletion', 'complex'")] string? variantType = null,
        [Description("Filter by chromosome name")] string? chrom = null,
        [Description("Minimum quality score")] double? minQuality = null,
        [Description("Only include PASS variants")] bool passOnly = false)
    {
        if (string.IsNullOrEmpty(content))
            throw new ArgumentException("Content cannot be null or empty", nameof(content));

        var records = VcfParser.Parse(content).ToList();
        IEnumerable<VcfParser.VcfRecord> filtered = records;

        // Filter by variant type
        if (!string.IsNullOrEmpty(variantType))
        {
            filtered = variantType.ToLowerInvariant() switch
            {
                "snp" => VcfParser.FilterSNPs(filtered),
                "indel" => VcfParser.FilterIndels(filtered),
                "insertion" => VcfParser.FilterByType(filtered, VcfParser.VariantType.Insertion),
                "deletion" => VcfParser.FilterByType(filtered, VcfParser.VariantType.Deletion),
                "complex" => VcfParser.FilterByType(filtered, VcfParser.VariantType.Complex),
                _ => throw new ArgumentException($"Invalid variant type: {variantType}. Use 'snp', 'indel', 'insertion', 'deletion', or 'complex'", nameof(variantType))
            };
        }

        // Filter by chromosome
        if (!string.IsNullOrEmpty(chrom))
        {
            filtered = VcfParser.FilterByChrom(filtered, chrom);
        }

        // Filter by quality
        if (minQuality.HasValue)
        {
            filtered = VcfParser.FilterByQuality(filtered, minQuality.Value);
        }

        // Filter PASS only
        if (passOnly)
        {
            filtered = VcfParser.FilterPassing(filtered);
        }

        var filteredList = filtered.ToList();
        var results = filteredList.Select(r => new VcfRecordResult(
            r.Chrom,
            r.Pos,
            r.Id,
            r.Ref,
            r.Alt.ToList(),
            r.Qual,
            r.Filter.ToList(),
            r.Info.ToDictionary(kv => kv.Key, kv => kv.Value),
            VcfParser.ClassifyVariant(r).ToString()
        )).ToList();

        return new VcfFilterResult(
            results,
            results.Count,
            records.Count,
            records.Count > 0 ? (double)results.Count / records.Count * 100 : 0);
    }

    // ========================
    // GFF/GTF Tools
    // ========================

    [McpServerTool(Name = "gff_parse")]
    [Description("Parse GFF3/GTF format content into feature records. Supports GFF3, GTF, and GFF2 formats for gene annotations.")]
    public static GffParseResult GffParse(
        [Description("GFF/GTF format content to parse")] string content,
        [Description("Format: 'gff3', 'gtf', 'gff2', or 'auto' (default)")] string format = "auto")
    {
        if (string.IsNullOrEmpty(content))
            throw new ArgumentException("Content cannot be null or empty", nameof(content));

        var gffFormat = format.ToLowerInvariant() switch
        {
            "gff3" => GffParser.GffFormat.GFF3,
            "gtf" => GffParser.GffFormat.GTF,
            "gff2" => GffParser.GffFormat.GFF2,
            "auto" => GffParser.GffFormat.Auto,
            _ => throw new ArgumentException($"Invalid format: {format}. Use 'gff3', 'gtf', 'gff2', or 'auto'", nameof(format))
        };

        var records = GffParser.Parse(content, gffFormat).ToList();
        var results = records.Select(r => new GffRecordResult(
            r.Seqid,
            r.Source,
            r.Type,
            r.Start,
            r.End,
            r.End - r.Start + 1,
            r.Score,
            r.Strand.ToString(),
            r.Phase,
            r.Attributes.ToDictionary(kv => kv.Key, kv => kv.Value),
            GffParser.GetGeneName(r)
        )).ToList();

        return new GffParseResult(results, results.Count);
    }

    [McpServerTool(Name = "gff_statistics")]
    [Description("Calculate statistics for GFF/GTF annotations. Returns feature type counts, sequence IDs, and gene/exon counts.")]
    public static GffStatisticsResult GffStatistics(
        [Description("GFF/GTF format content to analyze")] string content,
        [Description("Format: 'gff3', 'gtf', 'gff2', or 'auto' (default)")] string format = "auto")
    {
        if (string.IsNullOrEmpty(content))
            throw new ArgumentException("Content cannot be null or empty", nameof(content));

        var gffFormat = format.ToLowerInvariant() switch
        {
            "gff3" => GffParser.GffFormat.GFF3,
            "gtf" => GffParser.GffFormat.GTF,
            "gff2" => GffParser.GffFormat.GFF2,
            "auto" => GffParser.GffFormat.Auto,
            _ => throw new ArgumentException($"Invalid format: {format}. Use 'gff3', 'gtf', 'gff2', or 'auto'", nameof(format))
        };

        var records = GffParser.Parse(content, gffFormat).ToList();
        var stats = GffParser.CalculateStatistics(records);

        return new GffStatisticsResult(
            stats.TotalFeatures,
            stats.FeatureTypeCounts.ToDictionary(kv => kv.Key, kv => kv.Value),
            stats.SequenceIds.ToList(),
            stats.Sources.ToList(),
            stats.GeneCount,
            stats.ExonCount);
    }

    [McpServerTool(Name = "gff_filter")]
    [Description("Filter GFF/GTF records by feature type, sequence ID, or genomic region.")]
    public static GffFilterResult GffFilter(
        [Description("GFF/GTF format content to filter")] string content,
        [Description("Filter by feature type (e.g., 'gene', 'exon', 'CDS')")] string? featureType = null,
        [Description("Filter by sequence ID (chromosome/contig name)")] string? seqid = null,
        [Description("Filter by region start position (requires seqid and regionEnd)")] int? regionStart = null,
        [Description("Filter by region end position (requires seqid and regionStart)")] int? regionEnd = null,
        [Description("Format: 'gff3', 'gtf', 'gff2', or 'auto' (default)")] string format = "auto")
    {
        if (string.IsNullOrEmpty(content))
            throw new ArgumentException("Content cannot be null or empty", nameof(content));

        var gffFormat = format.ToLowerInvariant() switch
        {
            "gff3" => GffParser.GffFormat.GFF3,
            "gtf" => GffParser.GffFormat.GTF,
            "gff2" => GffParser.GffFormat.GFF2,
            "auto" => GffParser.GffFormat.Auto,
            _ => throw new ArgumentException($"Invalid format: {format}. Use 'gff3', 'gtf', 'gff2', or 'auto'", nameof(format))
        };

        var records = GffParser.Parse(content, gffFormat).ToList();
        IEnumerable<GffParser.GffRecord> filtered = records;

        // Filter by feature type
        if (!string.IsNullOrEmpty(featureType))
        {
            filtered = GffParser.FilterByType(filtered, featureType);
        }

        // Filter by sequence ID
        if (!string.IsNullOrEmpty(seqid))
        {
            filtered = GffParser.FilterBySeqid(filtered, seqid);
        }

        // Filter by region (requires seqid and both start/end)
        if (regionStart.HasValue && regionEnd.HasValue && !string.IsNullOrEmpty(seqid))
        {
            filtered = GffParser.FilterByRegion(filtered, seqid, regionStart.Value, regionEnd.Value);
        }

        var filteredList = filtered.ToList();
        var results = filteredList.Select(r => new GffRecordResult(
            r.Seqid,
            r.Source,
            r.Type,
            r.Start,
            r.End,
            r.End - r.Start + 1,
            r.Score,
            r.Strand.ToString(),
            r.Phase,
            r.Attributes.ToDictionary(kv => kv.Key, kv => kv.Value),
            GffParser.GetGeneName(r)
        )).ToList();

        return new GffFilterResult(
            results,
            results.Count,
            records.Count,
            records.Count > 0 ? (double)results.Count / records.Count * 100 : 0);
    }

    // ========================
    // GenBank Tools
    // ========================

    [McpServerTool(Name = "genbank_parse")]
    [Description("Parse GenBank flat file format into structured records. Returns locus info, definition, accession, features, and sequence.")]
    public static GenBankParseResult GenBankParse(
        [Description("GenBank format content to parse")] string content)
    {
        if (string.IsNullOrEmpty(content))
            throw new ArgumentException("Content cannot be null or empty", nameof(content));

        var records = GenBankParser.Parse(content).ToList();
        var results = records.Select(r => new GenBankRecordResult(
            r.Locus,
            r.SequenceLength,
            r.MoleculeType,
            r.Topology,
            r.Division,
            r.Date?.ToString("yyyy-MM-dd"),
            r.Definition,
            r.Accession,
            r.Version,
            r.Keywords.ToList(),
            r.Organism,
            r.Taxonomy,
            r.Features.Select(f => new GenBankFeatureResult(
                f.Key,
                f.Location.Start,
                f.Location.End,
                f.Location.IsComplement,
                f.Location.IsJoin,
                f.Location.RawLocation,
                f.Qualifiers.ToDictionary(kv => kv.Key, kv => kv.Value)
            )).ToList(),
            r.Sequence.Length,
            r.Sequence
        )).ToList();

        return new GenBankParseResult(results, results.Count);
    }

    [McpServerTool(Name = "genbank_features")]
    [Description("Extract features from GenBank records by feature type (gene, CDS, mRNA, etc.).")]
    public static GenBankFeaturesResult GenBankFeatures(
        [Description("GenBank format content to parse")] string content,
        [Description("Feature type to extract (e.g., 'gene', 'CDS', 'mRNA', 'exon')")] string? featureType = null)
    {
        if (string.IsNullOrEmpty(content))
            throw new ArgumentException("Content cannot be null or empty", nameof(content));

        var records = GenBankParser.Parse(content).ToList();
        var allFeatures = new List<GenBankFeatureResult>();

        foreach (var record in records)
        {
            IEnumerable<GenBankParser.Feature> features = record.Features;

            if (!string.IsNullOrEmpty(featureType))
            {
                features = GenBankParser.GetFeatures(record, featureType);
            }

            allFeatures.AddRange(features.Select(f => new GenBankFeatureResult(
                f.Key,
                f.Location.Start,
                f.Location.End,
                f.Location.IsComplement,
                f.Location.IsJoin,
                f.Location.RawLocation,
                f.Qualifiers.ToDictionary(kv => kv.Key, kv => kv.Value)
            )));
        }

        // Count by feature type
        var featureTypeCounts = allFeatures
            .GroupBy(f => f.Key)
            .ToDictionary(g => g.Key, g => g.Count());

        return new GenBankFeaturesResult(allFeatures, allFeatures.Count, featureTypeCounts);
    }

    [McpServerTool(Name = "genbank_statistics")]
    [Description("Calculate statistics for GenBank records. Returns counts of records, features, and sequence lengths.")]
    public static GenBankStatisticsResult GenBankStatistics(
        [Description("GenBank format content to analyze")] string content)
    {
        if (string.IsNullOrEmpty(content))
            throw new ArgumentException("Content cannot be null or empty", nameof(content));

        var records = GenBankParser.Parse(content).ToList();

        var featureTypeCounts = records
            .SelectMany(r => r.Features)
            .GroupBy(f => f.Key)
            .ToDictionary(g => g.Key, g => g.Count());

        var moleculeTypes = records
            .Select(r => r.MoleculeType)
            .Where(m => !string.IsNullOrEmpty(m))
            .Distinct()
            .ToList();

        var divisions = records
            .Select(r => r.Division)
            .Where(d => !string.IsNullOrEmpty(d))
            .Distinct()
            .ToList();

        var totalSequenceLength = records.Sum(r => r.SequenceLength);
        var totalFeatures = records.Sum(r => r.Features.Count);
        var geneCount = records.Sum(r => GenBankParser.GetGenes(r).Count());
        var cdsCount = records.Sum(r => GenBankParser.GetCDS(r).Count());

        return new GenBankStatisticsResult(
            records.Count,
            totalFeatures,
            totalSequenceLength,
            featureTypeCounts,
            moleculeTypes,
            divisions,
            geneCount,
            cdsCount);
    }

    // ========================
    // EMBL Tools
    // ========================

    [McpServerTool(Name = "embl_parse")]
    [Description("Parse EMBL flat file format into structured records. Returns accession, description, features, and sequence.")]
    public static EmblParseResult EmblParse(
        [Description("EMBL format content to parse")] string content)
    {
        if (string.IsNullOrEmpty(content))
            throw new ArgumentException("Content cannot be null or empty", nameof(content));

        var records = EmblParser.Parse(content).ToList();
        var results = records.Select(r => new EmblRecordResult(
            r.Accession,
            r.SequenceVersion,
            r.DataClass,
            r.MoleculeType,
            r.Topology,
            r.TaxonomicDivision,
            r.SequenceLength,
            r.Description,
            r.Keywords.ToList(),
            r.Organism,
            r.OrganismClassification.ToList(),
            r.Features.Select(f => new EmblFeatureResult(
                f.Key,
                f.Location.Start,
                f.Location.End,
                f.Location.IsComplement,
                f.Location.IsJoin,
                f.Location.RawLocation,
                f.Qualifiers.ToDictionary(kv => kv.Key, kv => kv.Value)
            )).ToList(),
            r.Sequence.Length,
            r.Sequence
        )).ToList();

        return new EmblParseResult(results, results.Count);
    }

    [McpServerTool(Name = "embl_features")]
    [Description("Extract features from EMBL records by feature type (gene, CDS, mRNA, etc.).")]
    public static EmblFeaturesResult EmblFeatures(
        [Description("EMBL format content to parse")] string content,
        [Description("Feature type to extract (e.g., 'gene', 'CDS', 'mRNA', 'exon')")] string? featureType = null)
    {
        if (string.IsNullOrEmpty(content))
            throw new ArgumentException("Content cannot be null or empty", nameof(content));

        var records = EmblParser.Parse(content).ToList();
        var allFeatures = new List<EmblFeatureResult>();

        foreach (var record in records)
        {
            IEnumerable<EmblParser.Feature> features = record.Features;

            if (!string.IsNullOrEmpty(featureType))
            {
                features = EmblParser.GetFeatures(record, featureType);
            }

            allFeatures.AddRange(features.Select(f => new EmblFeatureResult(
                f.Key,
                f.Location.Start,
                f.Location.End,
                f.Location.IsComplement,
                f.Location.IsJoin,
                f.Location.RawLocation,
                f.Qualifiers.ToDictionary(kv => kv.Key, kv => kv.Value)
            )));
        }

        var featureTypeCounts = allFeatures
            .GroupBy(f => f.Key)
            .ToDictionary(g => g.Key, g => g.Count());

        return new EmblFeaturesResult(allFeatures, allFeatures.Count, featureTypeCounts);
    }

    [McpServerTool(Name = "embl_statistics")]
    [Description("Calculate statistics for EMBL records. Returns counts of records, features, and sequence lengths.")]
    public static EmblStatisticsResult EmblStatistics(
        [Description("EMBL format content to analyze")] string content)
    {
        if (string.IsNullOrEmpty(content))
            throw new ArgumentException("Content cannot be null or empty", nameof(content));

        var records = EmblParser.Parse(content).ToList();

        var featureTypeCounts = records
            .SelectMany(r => r.Features)
            .GroupBy(f => f.Key)
            .ToDictionary(g => g.Key, g => g.Count());

        var moleculeTypes = records
            .Select(r => r.MoleculeType)
            .Where(m => !string.IsNullOrEmpty(m))
            .Distinct()
            .ToList();

        var divisions = records
            .Select(r => r.TaxonomicDivision)
            .Where(d => !string.IsNullOrEmpty(d))
            .Distinct()
            .ToList();

        var totalSequenceLength = records.Sum(r => r.SequenceLength);
        var totalFeatures = records.Sum(r => r.Features.Count);
        var geneCount = records.Sum(r => EmblParser.GetGenes(r).Count());
        var cdsCount = records.Sum(r => EmblParser.GetCDS(r).Count());

        return new EmblStatisticsResult(
            records.Count,
            totalFeatures,
            totalSequenceLength,
            featureTypeCounts,
            moleculeTypes,
            divisions,
            geneCount,
            cdsCount);
    }
}

// ========================
// Result Records
// ========================

public record FastaEntryResult(string Id, string? Description, string Sequence, int Length);
public record FastaParseResult(List<FastaEntryResult> Entries, int Count);
public record FastaFormatResult(string Fasta);
public record FastaWriteResult(string FilePath, int EntriesWritten, int TotalBases);
public record FastqRecordResult(string Id, string? Description, string Sequence, string QualityString, List<int> QualityScores, int Length);
public record FastqParseResult(List<FastqRecordResult> Entries, int Count);
public record FastqStatisticsResult(
    int TotalReads,
    long TotalBases,
    double MeanReadLength,
    double MeanQuality,
    int MinReadLength,
    int MaxReadLength,
    double Q20Percentage,
    double Q30Percentage,
    double GcContent);
public record FastqFilterResult(
    List<FastqRecordResult> Entries,
    int PassedCount,
    int TotalCount,
    double PassedPercentage);
public record FastqDetectEncodingResult(string Encoding, int Offset);
public record FastqEncodeQualityResult(string QualityString, int Length);
public record FastqPhredToErrorResult(int PhredScore, double ErrorProbability);
public record FastqErrorToPhredResult(double ErrorProbability, int PhredScore);
public record FastqTrimQualityResult(
    List<FastqRecordResult> Entries,
    int Count,
    long OriginalBases,
    long TrimmedBases,
    double TrimmedPercentage);
public record FastqTrimAdapterResult(
    List<FastqRecordResult> Entries,
    int Count,
    int ReadsWithAdapter,
    long OriginalBases,
    long TrimmedBases);
public record FastqFormatResult(string Fastq);
public record FastqWriteResult(string FilePath, int RecordsWritten, long TotalBases);
public record BedRecordResult(
    string Chrom,
    int ChromStart,
    int ChromEnd,
    int Length,
    string? Name = null,
    int? Score = null,
    string? Strand = null,
    int? ThickStart = null,
    int? ThickEnd = null,
    string? ItemRgb = null,
    int? BlockCount = null,
    List<int>? BlockSizes = null,
    List<int>? BlockStarts = null);
public record BedParseResult(List<BedRecordResult> Records, int Count);
public record BedFilterResult(
    List<BedRecordResult> Records,
    int PassedCount,
    int TotalCount,
    double PassedPercentage);
public record BedMergeResult(List<BedRecordResult> Records, int MergedCount, int OriginalCount);
public record BedIntersectResult(List<BedRecordResult> Records, int IntersectionCount, int CountA, int CountB);
public record VcfRecordResult(
    string Chrom,
    int Pos,
    string Id,
    string Ref,
    List<string> Alt,
    double? Qual,
    List<string> Filter,
    Dictionary<string, string> Info,
    string VariantType);
public record VcfParseResult(List<VcfRecordResult> Records, int Count);
public record VcfStatisticsResult(
    int TotalVariants,
    int SnpCount,
    int IndelCount,
    int ComplexCount,
    int PassingCount,
    Dictionary<string, int> ChromosomeCounts,
    double? MeanQuality);
public record VcfFilterResult(
    List<VcfRecordResult> Records,
    int PassedCount,
    int TotalCount,
    double PassedPercentage);
public record GffRecordResult(
    string Seqid,
    string Source,
    string Type,
    int Start,
    int End,
    int Length,
    double? Score,
    string Strand,
    int? Phase,
    Dictionary<string, string> Attributes,
    string? GeneName);
public record GffParseResult(List<GffRecordResult> Records, int Count);
public record GffStatisticsResult(
    int TotalFeatures,
    Dictionary<string, int> FeatureTypeCounts,
    List<string> SequenceIds,
    List<string> Sources,
    int GeneCount,
    int ExonCount);
public record GffFilterResult(
    List<GffRecordResult> Records,
    int PassedCount,
    int TotalCount,
    double PassedPercentage);
public record GenBankFeatureResult(
    string Key,
    int Start,
    int End,
    bool IsComplement,
    bool IsJoin,
    string RawLocation,
    Dictionary<string, string> Qualifiers);
public record GenBankRecordResult(
    string Locus,
    int SequenceLength,
    string MoleculeType,
    string Topology,
    string Division,
    string? Date,
    string Definition,
    string Accession,
    string Version,
    List<string> Keywords,
    string Organism,
    string Taxonomy,
    List<GenBankFeatureResult> Features,
    int ActualSequenceLength,
    string Sequence);
public record GenBankParseResult(List<GenBankRecordResult> Records, int Count);
public record GenBankFeaturesResult(
    List<GenBankFeatureResult> Features,
    int Count,
    Dictionary<string, int> FeatureTypeCounts);
public record GenBankStatisticsResult(
    int RecordCount,
    int TotalFeatures,
    int TotalSequenceLength,
    Dictionary<string, int> FeatureTypeCounts,
    List<string> MoleculeTypes,
    List<string> Divisions,
    int GeneCount,
    int CdsCount);
public record EmblFeatureResult(
    string Key,
    int Start,
    int End,
    bool IsComplement,
    bool IsJoin,
    string RawLocation,
    Dictionary<string, string> Qualifiers);
public record EmblRecordResult(
    string Accession,
    string SequenceVersion,
    string DataClass,
    string MoleculeType,
    string Topology,
    string TaxonomicDivision,
    int SequenceLength,
    string Description,
    List<string> Keywords,
    string Organism,
    List<string> OrganismClassification,
    List<EmblFeatureResult> Features,
    int ActualSequenceLength,
    string Sequence);
public record EmblParseResult(List<EmblRecordResult> Records, int Count);
public record EmblFeaturesResult(
    List<EmblFeatureResult> Features,
    int Count,
    Dictionary<string, int> FeatureTypeCounts);
public record EmblStatisticsResult(
    int RecordCount,
    int TotalFeatures,
    int TotalSequenceLength,
    Dictionary<string, int> FeatureTypeCounts,
    List<string> MoleculeTypes,
    List<string> Divisions,
    int GeneCount,
    int CdsCount);
