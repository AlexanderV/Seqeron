using System.ComponentModel;
using ModelContextProtocol.Server;
using Seqeron.Genomics;

namespace Seqeron.Mcp.Parsers.Tools;

[McpServerToolType]
public class ParsersTools
{
    // ========================
    // FASTA Tools
    // ========================

    [McpServerTool(Name = "fasta_parse", Title = "FASTA — Parse", ReadOnly = true)]
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

    [McpServerTool(Name = "fasta_format", Title = "FASTA — Format", ReadOnly = true)]
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

    [McpServerTool(Name = "fasta_write", Title = "FASTA — Write to File")]
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

    [McpServerTool(Name = "fastq_parse", Title = "FASTQ — Parse", ReadOnly = true)]
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

    [McpServerTool(Name = "fastq_statistics", Title = "FASTQ — Statistics", ReadOnly = true)]
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

    [McpServerTool(Name = "fastq_filter", Title = "FASTQ — Filter by Quality", ReadOnly = true)]
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

    [McpServerTool(Name = "fastq_detect_encoding", Title = "FASTQ — Detect Encoding", ReadOnly = true)]
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

    [McpServerTool(Name = "fastq_encode_quality", Title = "FASTQ — Encode Quality", ReadOnly = true)]
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

    [McpServerTool(Name = "fastq_phred_to_error", Title = "FASTQ — Phred to Error Rate", ReadOnly = true)]
    [Description("Convert Phred quality score to error probability. Formula: P = 10^(-Q/10).")]
    public static FastqPhredToErrorResult FastqPhredToError(
        [Description("Phred quality score")] int phredScore)
    {
        if (phredScore < 0)
            throw new ArgumentException("Phred score cannot be negative", nameof(phredScore));

        var errorProbability = FastqParser.PhredToErrorProbability(phredScore);
        return new FastqPhredToErrorResult(phredScore, errorProbability);
    }

    [McpServerTool(Name = "fastq_error_to_phred", Title = "FASTQ — Error to Phred Score", ReadOnly = true)]
    [Description("Convert error probability to Phred quality score. Formula: Q = -10 * log10(P).")]
    public static FastqErrorToPhredResult FastqErrorToPhred(
        [Description("Error probability (0-1)")] double errorProbability)
    {
        if (errorProbability < 0 || errorProbability > 1)
            throw new ArgumentException("Error probability must be between 0 and 1", nameof(errorProbability));

        var phredScore = FastqParser.ErrorProbabilityToPhred(errorProbability);
        return new FastqErrorToPhredResult(errorProbability, phredScore);
    }

    [McpServerTool(Name = "fastq_trim_quality", Title = "FASTQ — Trim by Quality", ReadOnly = true)]
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

    [McpServerTool(Name = "fastq_trim_adapter", Title = "FASTQ — Trim Adapter", ReadOnly = true)]
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

    [McpServerTool(Name = "fastq_format", Title = "FASTQ — Format", ReadOnly = true)]
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

    [McpServerTool(Name = "fastq_write", Title = "FASTQ — Write to File")]
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

    [McpServerTool(Name = "bed_parse", Title = "BED — Parse", ReadOnly = true)]
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

    [McpServerTool(Name = "bed_filter", Title = "BED — Filter", ReadOnly = true)]
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

    [McpServerTool(Name = "bed_merge", Title = "BED — Merge Overlapping", ReadOnly = true)]
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

    [McpServerTool(Name = "bed_intersect", Title = "BED — Intersect", ReadOnly = true)]
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

    [McpServerTool(Name = "vcf_parse", Title = "VCF — Parse", ReadOnly = true)]
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

    [McpServerTool(Name = "vcf_statistics", Title = "VCF — Statistics", ReadOnly = true)]
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

    [McpServerTool(Name = "vcf_filter", Title = "VCF — Filter", ReadOnly = true)]
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

    [McpServerTool(Name = "vcf_classify")]
    [Description("Classify variant type for a VCF record. Returns SNP, MNP, Insertion, Deletion, Complex, Symbolic, or Unknown.")]
    public static VcfClassifyResult VcfClassify(
        [Description("Reference allele")] string refAllele,
        [Description("Alternate allele")] string altAllele)
    {
        if (string.IsNullOrEmpty(refAllele))
            throw new ArgumentException("Reference allele cannot be null or empty", nameof(refAllele));
        if (string.IsNullOrEmpty(altAllele))
            throw new ArgumentException("Alternate allele cannot be null or empty", nameof(altAllele));

        // Create a minimal record for classification
        var record = new VcfParser.VcfRecord(
            "chr1", 1, ".", refAllele, new[] { altAllele }, null, Array.Empty<string>(),
            new Dictionary<string, string>());

        var variantType = VcfParser.ClassifyVariant(record);
        var length = VcfParser.GetVariantLength(record);

        return new VcfClassifyResult(variantType.ToString(), refAllele.Length, altAllele.Length, length);
    }

    [McpServerTool(Name = "vcf_is_snp", Title = "VCF — Is SNP", ReadOnly = true)]
    [Description("Check if a variant is a SNP (Single Nucleotide Polymorphism).")]
    public static VcfIsSNPResult VcfIsSNP(
        [Description("Reference allele")] string refAllele,
        [Description("Alternate allele")] string altAllele)
    {
        if (string.IsNullOrEmpty(refAllele))
            throw new ArgumentException("Reference allele cannot be null or empty", nameof(refAllele));
        if (string.IsNullOrEmpty(altAllele))
            throw new ArgumentException("Alternate allele cannot be null or empty", nameof(altAllele));

        var record = new VcfParser.VcfRecord(
            "chr1", 1, ".", refAllele, new[] { altAllele }, null, Array.Empty<string>(),
            new Dictionary<string, string>());

        var isSNP = VcfParser.IsSNP(record);
        return new VcfIsSNPResult(isSNP, refAllele, altAllele);
    }

    [McpServerTool(Name = "vcf_is_indel", Title = "VCF — Is Indel", ReadOnly = true)]
    [Description("Check if a variant is an indel (insertion or deletion).")]
    public static VcfIsIndelResult VcfIsIndel(
        [Description("Reference allele")] string refAllele,
        [Description("Alternate allele")] string altAllele)
    {
        if (string.IsNullOrEmpty(refAllele))
            throw new ArgumentException("Reference allele cannot be null or empty", nameof(refAllele));
        if (string.IsNullOrEmpty(altAllele))
            throw new ArgumentException("Alternate allele cannot be null or empty", nameof(altAllele));

        var record = new VcfParser.VcfRecord(
            "chr1", 1, ".", refAllele, new[] { altAllele }, null, Array.Empty<string>(),
            new Dictionary<string, string>());

        var isIndel = VcfParser.IsIndel(record);
        var variantType = VcfParser.ClassifyVariant(record);
        var isInsertion = variantType == VcfParser.VariantType.Insertion;
        var isDeletion = variantType == VcfParser.VariantType.Deletion;

        return new VcfIsIndelResult(isIndel, isInsertion, isDeletion, refAllele, altAllele);
    }

    [McpServerTool(Name = "vcf_variant_length", Title = "VCF — Variant Length", ReadOnly = true)]
    [Description("Get the length difference of a variant (absolute difference between ref and alt lengths).")]
    public static VcfVariantLengthResult VcfVariantLength(
        [Description("Reference allele")] string refAllele,
        [Description("Alternate allele")] string altAllele)
    {
        if (string.IsNullOrEmpty(refAllele))
            throw new ArgumentException("Reference allele cannot be null or empty", nameof(refAllele));
        if (string.IsNullOrEmpty(altAllele))
            throw new ArgumentException("Alternate allele cannot be null or empty", nameof(altAllele));

        var record = new VcfParser.VcfRecord(
            "chr1", 1, ".", refAllele, new[] { altAllele }, null, Array.Empty<string>(),
            new Dictionary<string, string>());

        var length = VcfParser.GetVariantLength(record);
        return new VcfVariantLengthResult(length, refAllele.Length, altAllele.Length);
    }

    [McpServerTool(Name = "vcf_is_hom_ref")]
    [Description("Check if a genotype is homozygous reference (0/0 or 0|0).")]
    public static VcfGenotypeCheckResult VcfIsHomRef(
        [Description("Genotype string (e.g., '0/0', '0/1', '1/1')")] string genotype)
    {
        if (string.IsNullOrEmpty(genotype))
            throw new ArgumentException("Genotype cannot be null or empty", nameof(genotype));

        var isHomRef = genotype == "0/0" || genotype == "0|0";
        return new VcfGenotypeCheckResult(isHomRef, genotype, "HomozygousReference");
    }

    [McpServerTool(Name = "vcf_is_hom_alt")]
    [Description("Check if a genotype is homozygous alternate (e.g., 1/1, 2/2).")]
    public static VcfGenotypeCheckResult VcfIsHomAlt(
        [Description("Genotype string (e.g., '0/0', '0/1', '1/1')")] string genotype)
    {
        if (string.IsNullOrEmpty(genotype))
            throw new ArgumentException("Genotype cannot be null or empty", nameof(genotype));

        var alleles = genotype.Replace('|', '/').Split('/');
        var isHomAlt = alleles.Length == 2 &&
                       alleles[0] == alleles[1] &&
                       alleles[0] != "0" &&
                       alleles[0] != ".";

        return new VcfGenotypeCheckResult(isHomAlt, genotype, "HomozygousAlternate");
    }

    [McpServerTool(Name = "vcf_is_het")]
    [Description("Check if a genotype is heterozygous (different alleles).")]
    public static VcfGenotypeCheckResult VcfIsHet(
        [Description("Genotype string (e.g., '0/0', '0/1', '1/1')")] string genotype)
    {
        if (string.IsNullOrEmpty(genotype))
            throw new ArgumentException("Genotype cannot be null or empty", nameof(genotype));

        var alleles = genotype.Replace('|', '/').Split('/');
        var isHet = alleles.Length == 2 && alleles[0] != alleles[1];

        return new VcfGenotypeCheckResult(isHet, genotype, "Heterozygous");
    }

    [McpServerTool(Name = "vcf_has_flag", Title = "VCF — Has Info Flag", ReadOnly = true)]
    [Description("Check if a VCF INFO field flag is present.")]
    public static VcfHasFlagResult VcfHasFlag(
        [Description("VCF format content")] string content,
        [Description("INFO field flag name to check")] string flag)
    {
        if (string.IsNullOrEmpty(content))
            throw new ArgumentException("Content cannot be null or empty", nameof(content));
        if (string.IsNullOrEmpty(flag))
            throw new ArgumentException("Flag cannot be null or empty", nameof(flag));

        var records = VcfParser.Parse(content).ToList();
        var recordsWithFlag = records.Where(r => VcfParser.HasInfoFlag(r, flag)).ToList();

        return new VcfHasFlagResult(
            flag,
            recordsWithFlag.Count,
            records.Count,
            records.Count > 0 ? (double)recordsWithFlag.Count / records.Count * 100 : 0);
    }

    [McpServerTool(Name = "vcf_write", Title = "VCF — Write to File")]
    [Description("Write VCF records to a file. Creates or overwrites the file at specified path.")]
    public static VcfWriteResult VcfWrite(
        [Description("File path to write VCF output")] string filePath,
        [Description("VCF format content to write")] string content)
    {
        if (string.IsNullOrEmpty(filePath))
            throw new ArgumentException("File path cannot be null or empty", nameof(filePath));
        if (string.IsNullOrEmpty(content))
            throw new ArgumentException("Content cannot be null or empty", nameof(content));

        var (header, records) = VcfParser.ParseWithHeader(content);
        var recordsList = records.ToList();
        VcfParser.WriteToFile(filePath, recordsList, header);

        return new VcfWriteResult(filePath, recordsList.Count);
    }

    // ========================
    // GFF/GTF Tools
    // ========================

    [McpServerTool(Name = "gff_parse", Title = "GFF — Parse", ReadOnly = true)]
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

    [McpServerTool(Name = "gff_statistics", Title = "GFF — Statistics", ReadOnly = true)]
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

    [McpServerTool(Name = "gff_filter", Title = "GFF — Filter", ReadOnly = true)]
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

    [McpServerTool(Name = "genbank_parse", Title = "GenBank — Parse", ReadOnly = true)]
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

    [McpServerTool(Name = "genbank_features", Title = "GenBank — Extract Features", ReadOnly = true)]
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

    [McpServerTool(Name = "genbank_statistics", Title = "GenBank — Statistics", ReadOnly = true)]
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

    [McpServerTool(Name = "genbank_parse_location", Title = "GenBank — Parse Location", ReadOnly = true)]
    [Description("Parse a GenBank feature location string into its components. Handles simple ranges, complement, join, and complex locations.")]
    public static GenBankLocationResult GenBankParseLocation(
        [Description("Feature location string (e.g., '100..200', 'complement(100..200)', 'join(100..200,300..400)')")] string locationString)
    {
        if (string.IsNullOrEmpty(locationString))
            throw new ArgumentException("Location string cannot be null or empty", nameof(locationString));

        var location = GenBankParser.ParseLocation(locationString);

        var parts = location.Parts.Select(p => new LocationPartResult(p.Start, p.End, p.End - p.Start + 1)).ToList();

        return new GenBankLocationResult(
            location.Start,
            location.End,
            location.End - location.Start + 1,
            location.IsComplement,
            location.IsJoin,
            parts,
            location.RawLocation);
    }

    [McpServerTool(Name = "genbank_extract_sequence", Title = "GenBank — Extract Sequence", ReadOnly = true)]
    [Description("Extract a subsequence from a GenBank record based on a feature location string. Handles complement and join locations.")]
    public static GenBankExtractSequenceResult GenBankExtractSequence(
        [Description("GenBank format content")] string content,
        [Description("Feature location string (e.g., '100..200', 'complement(100..200)', 'join(100..200,300..400)')")] string locationString)
    {
        if (string.IsNullOrEmpty(content))
            throw new ArgumentException("Content cannot be null or empty", nameof(content));
        if (string.IsNullOrEmpty(locationString))
            throw new ArgumentException("Location string cannot be null or empty", nameof(locationString));

        var records = GenBankParser.Parse(content).ToList();
        if (records.Count == 0)
            throw new ArgumentException("No GenBank records found in content", nameof(content));

        var record = records[0];
        var location = GenBankParser.ParseLocation(locationString);
        var sequence = GenBankParser.ExtractSequence(record, location);

        return new GenBankExtractSequenceResult(
            sequence,
            sequence.Length,
            location.IsComplement,
            location.IsJoin,
            locationString);
    }

    // ========================
    // EMBL Tools
    // ========================

    [McpServerTool(Name = "embl_parse", Title = "EMBL — Parse", ReadOnly = true)]
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

    [McpServerTool(Name = "embl_features", Title = "EMBL — Extract Features", ReadOnly = true)]
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

    [McpServerTool(Name = "embl_statistics", Title = "EMBL — Statistics", ReadOnly = true)]
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
