namespace Seqeron.Mcp.Parsers.Tools;

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
public record VcfClassifyResult(string VariantType, int RefLength, int AltLength, int LengthDifference);
public record VcfIsSNPResult(bool IsSNP, string RefAllele, string AltAllele);
public record VcfIsIndelResult(bool IsIndel, bool IsInsertion, bool IsDeletion, string RefAllele, string AltAllele);
public record VcfVariantLengthResult(int Length, int RefLength, int AltLength);
public record VcfGenotypeCheckResult(bool Result, string Genotype, string CheckType);
public record VcfHasFlagResult(string Flag, int RecordsWithFlag, int TotalRecords, double Percentage);
public record VcfWriteResult(string FilePath, int RecordsWritten);
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
public record LocationPartResult(int Start, int End, int Length);
public record GenBankLocationResult(
    int Start,
    int End,
    int Length,
    bool IsComplement,
    bool IsJoin,
    List<LocationPartResult> Parts,
    string RawLocation);
public record GenBankExtractSequenceResult(
    string Sequence,
    int Length,
    bool IsComplement,
    bool IsJoin,
    string Location);
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
