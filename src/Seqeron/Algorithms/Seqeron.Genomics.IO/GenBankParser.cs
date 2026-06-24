using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Text;
using System.Text.RegularExpressions;

namespace Seqeron.Genomics.IO;

/// <summary>
/// Parser for GenBank flat file format (.gb, .gbk, .genbank).
/// Supports parsing of sequence records with annotations, features, and references.
/// </summary>
public static partial class GenBankParser
{
    #region Records

    /// <summary>Represents a complete GenBank record</summary>
    public readonly record struct GenBankRecord(
        string Locus,
        int SequenceLength,
        string MoleculeType,
        string Topology,
        string Division,
        DateTime? Date,
        string Definition,
        string Accession,
        string Version,
        IReadOnlyList<string> Keywords,
        string Organism,
        string Taxonomy,
        IReadOnlyList<Reference> References,
        IReadOnlyList<Feature> Features,
        string Sequence,
        IReadOnlyDictionary<string, string> AdditionalFields);

    /// <summary>Literature reference</summary>
    public readonly record struct Reference(
        int Number,
        string Authors,
        string Title,
        string Journal,
        string PubMed,
        int? BaseFrom,
        int? BaseTo);

    /// <summary>Sequence feature with location and qualifiers</summary>
    public readonly record struct Feature(
        string Key,
        Location Location,
        IReadOnlyDictionary<string, string> Qualifiers);

    /// <summary>Feature location (supports joins, complements, orders, and partial indicators)</summary>
    public readonly record struct Location(
        int Start,
        int End,
        bool IsComplement,
        bool IsJoin,
        bool IsOrder,
        bool Is5PrimePartial,
        bool Is3PrimePartial,
        IReadOnlyList<(int Start, int End)> Parts,
        string RawLocation);

    #endregion

    #region Main Parsing Methods

    /// <summary>
    /// Parses GenBank records from a file.
    /// </summary>
    public static IEnumerable<GenBankRecord> ParseFile(string filePath)
    {
        if (string.IsNullOrEmpty(filePath))
            yield break;

        if (!File.Exists(filePath))
            yield break;

        foreach (var record in Parse(File.ReadAllText(filePath)))
        {
            yield return record;
        }
    }

    /// <summary>
    /// Parses GenBank records from text content.
    /// Multiple records are separated by // delimiter.
    /// </summary>
    public static IEnumerable<GenBankRecord> Parse(string content)
    {
        if (string.IsNullOrEmpty(content))
            yield break;

        // Split by record delimiter
        var recordTexts = content.Split(new[] { "\n//" }, StringSplitOptions.RemoveEmptyEntries);

        foreach (var recordText in recordTexts)
        {
            var trimmed = recordText.Trim();
            if (!string.IsNullOrEmpty(trimmed) && trimmed.StartsWith("LOCUS", StringComparison.Ordinal))
            {
                var record = ParseRecord(trimmed);
                if (record.HasValue)
                    yield return record.Value;
            }
        }
    }

    /// <summary>
    /// Parses a single GenBank record.
    /// </summary>
    private static GenBankRecord? ParseRecord(string text)
    {
        var lines = text.Split('\n').Select(l => l.TrimEnd('\r')).ToList();

        // Parse LOCUS line
        var locusInfo = ParseLocusLine(lines.FirstOrDefault(l => l.StartsWith("LOCUS", StringComparison.Ordinal)) ?? "");

        // Collect sections
        var sections = ExtractSections(lines);

        var definition = sections.GetValueOrDefault("DEFINITION", "").Trim();
        var accession = sections.GetValueOrDefault("ACCESSION", "").Trim().Split(' ')[0];
        var version = sections.GetValueOrDefault("VERSION", "").Trim();

        // Keywords
        var keywordsText = sections.GetValueOrDefault("KEYWORDS", "");
        var keywords = ParseKeywords(keywordsText);

        // Organism and taxonomy from SOURCE section
        var sourceSection = sections.GetValueOrDefault("SOURCE", "");
        var (organism, taxonomy) = ParseSource(sourceSection);

        // References
        var references = ParseReferences(sections.GetValueOrDefault("REFERENCE", ""));

        // Features
        var features = ParseFeatures(sections.GetValueOrDefault("FEATURES", ""));

        // Sequence
        var sequence = ParseSequence(sections.GetValueOrDefault("ORIGIN", ""));

        // Additional fields
        var additionalFields = new Dictionary<string, string>();
        foreach (var (key, value) in sections)
        {
            if (!IsStandardField(key))
            {
                additionalFields[key] = value.Trim();
            }
        }

        return new GenBankRecord(
            locusInfo.Locus,
            locusInfo.Length,
            locusInfo.MoleculeType,
            locusInfo.Topology,
            locusInfo.Division,
            locusInfo.Date,
            definition,
            accession,
            version,
            keywords,
            organism,
            taxonomy,
            references,
            features,
            sequence,
            additionalFields
        );
    }

    #endregion

    #region Section Parsing

    private static Dictionary<string, string> ExtractSections(List<string> lines)
    {
        var sections = new Dictionary<string, string>();
        string currentSection = "";
        var currentContent = new StringBuilder();

        foreach (var line in lines)
        {
            // Check for section headers: standard 12-char column or short ORIGIN/special sections
            bool isNewSection = (line.Length >= 12 && !char.IsWhiteSpace(line[0])) ||
                                line.TrimEnd() == "ORIGIN" ||
                                (line.StartsWith("ORIGIN", StringComparison.Ordinal) &&
                                 (line.Length == 6 || char.IsWhiteSpace(line[6])));

            if (isNewSection)
            {
                // New section - save previous
                if (!string.IsNullOrEmpty(currentSection))
                {
                    SaveSection(sections, currentSection, currentContent.ToString());
                }

                var match = SectionHeaderRegex().Match(line);
                if (match.Success)
                {
                    currentSection = match.Groups[1].Value;
                    currentContent = new StringBuilder(line.Length > 12 ? line[12..] : "");
                }
                else
                {
                    // Handle short section names like ORIGIN without trailing content
                    currentSection = line.Trim();
                    currentContent = new StringBuilder();
                }
            }
            else if (!string.IsNullOrEmpty(currentSection))
            {
                // Continuation of current section
                currentContent.AppendLine();
                // FEATURES section needs to preserve leading spaces for proper parsing
                if (currentSection == "FEATURES")
                {
                    currentContent.Append(line);
                }
                else
                {
                    currentContent.Append(line.TrimStart());
                }
            }
        }

        if (!string.IsNullOrEmpty(currentSection))
        {
            SaveSection(sections, currentSection, currentContent.ToString());
        }

        return sections;
    }

    /// <summary>
    /// Saves a section to the dictionary, appending for REFERENCE sections (which can appear multiple times).
    /// </summary>
    private static void SaveSection(Dictionary<string, string> sections, string sectionName, string content)
    {
        if (sectionName == "REFERENCE")
        {
            // REFERENCE sections can appear multiple times - concatenate with newline + section marker
            if (sections.TryGetValue(sectionName, out var existing))
            {
                sections[sectionName] = existing + "\nREFERENCE   " + content;
            }
            else
            {
                sections[sectionName] = content;
            }
        }
        else
        {
            sections[sectionName] = content;
        }
    }

    // NCBI GenBank division codes
    // Source: https://www.ncbi.nlm.nih.gov/genbank/ and INSDC standards
    // "UNK" is not an official INSDC code but appears in real-world GenBank files for unclassified entries.
    private static readonly HashSet<string> KnownDivisionCodes = new(StringComparer.Ordinal)
    {
        "PRI", "ROD", "MAM", "VRT", "INV", "PLN", "BCT", "VRL", "PHG",
        "SYN", "UNA", "EST", "PAT", "STS", "GSS", "HTG", "HTC", "ENV",
        "CON", "TSA", "UNK"
    };

    private static (string Locus, int Length, string MoleculeType, string Topology, string Division, DateTime? Date)
        ParseLocusLine(string line)
    {
        if (string.IsNullOrEmpty(line))
            return ("", 0, "", "", "", null);

        // LOCUS format: LOCUS name length bp type topology division date
        var parts = line.Split(new[] { ' ' }, StringSplitOptions.RemoveEmptyEntries);

        string locus = parts.Length > 1 ? parts[1] : "";
        int length = 0;
        string moleculeType = "";
        string topology = "";
        string division = "";
        DateTime? date = null;

        if (parts.Length > 2 && int.TryParse(parts[2], out var len))
        {
            length = len;
        }

        // Find molecule type, topology, division, date
        foreach (var part in parts.Skip(3))
        {
            if (part is "DNA" or "RNA" or "mRNA" or "rRNA" or "tRNA"
                    or "ss-DNA" or "ds-DNA" or "ss-RNA" or "ds-RNA" or "cRNA")
            {
                moleculeType = part;
            }
            else if (part is "linear" or "circular")
            {
                topology = part;
            }
            else if (KnownDivisionCodes.Contains(part))
            {
                division = part;
            }
            else if (DateTime.TryParseExact(part, new[] { "dd-MMM-yyyy", "dd-MMM-yy" },
                     CultureInfo.InvariantCulture, DateTimeStyles.None, out var d))
            {
                date = d;
            }
        }

        return (locus, length, moleculeType, topology, division, date);
    }

    private static IReadOnlyList<string> ParseKeywords(string text)
    {
        if (string.IsNullOrWhiteSpace(text) || text.Trim() == ".")
            return Array.Empty<string>();

        return text
            .Replace("\n", " ")
            .Split(new[] { ';', ',' }, StringSplitOptions.RemoveEmptyEntries)
            .Select(k => k.Trim().TrimEnd('.'))
            .Where(k => !string.IsNullOrEmpty(k))
            .ToList();
    }

    private static (string Organism, string Taxonomy) ParseSource(string text)
    {
        if (string.IsNullOrWhiteSpace(text))
            return ("", "");

        var lines = text.Split('\n');
        var organism = "";
        var taxonomy = new StringBuilder();
        var inTaxonomy = false;

        foreach (var line in lines)
        {
            var trimmed = line.Trim();
            if (trimmed.StartsWith("ORGANISM", StringComparison.Ordinal))
            {
                organism = trimmed.Length > 9 ? trimmed[9..].Trim() : "";
                inTaxonomy = true;
            }
            else if (inTaxonomy && !string.IsNullOrEmpty(trimmed))
            {
                if (taxonomy.Length > 0) taxonomy.Append(' ');
                taxonomy.Append(trimmed);
            }
        }

        return (organism, taxonomy.ToString().TrimEnd('.'));
    }

    private static IReadOnlyList<Reference> ParseReferences(string text)
    {
        if (string.IsNullOrWhiteSpace(text))
            return Array.Empty<Reference>();

        var references = new List<Reference>();
        var refBlocks = text.Split(new[] { "\nREFERENCE" }, StringSplitOptions.None);

        // First block doesn't have "REFERENCE" prefix (it was stripped by ExtractSections)
        // Add it back so all blocks have consistent format
        var blocks = new List<string> { "REFERENCE   " + refBlocks[0] };
        blocks.AddRange(refBlocks.Skip(1).Select(b => "REFERENCE" + b));

        foreach (var block in blocks.Where(b => !string.IsNullOrWhiteSpace(b)))
        {
            var refNum = 0;
            string authors = "", title = "", journal = "", pubmed = "";
            int? baseFrom = null, baseTo = null;

            var lines = block.Split('\n');
            string currentField = "";
            var currentValue = new StringBuilder();

            foreach (var line in lines)
            {
                var trimmed = line.Trim();
                if (line.Length >= 12 && !char.IsWhiteSpace(line[0]) && line.Length > 2)
                {
                    // Save previous field
                    SaveRefField(ref authors, ref title, ref journal, ref pubmed, currentField, currentValue.ToString());

                    var spaceIdx = trimmed.IndexOf(' ');
                    currentField = spaceIdx > 0 ? trimmed[..spaceIdx] : trimmed;
                    currentValue = new StringBuilder(spaceIdx > 0 ? trimmed[(spaceIdx + 1)..] : "");

                    // Parse reference number and bases
                    if (currentField == "REFERENCE" || currentField.StartsWith("REFERENCE", StringComparison.Ordinal))
                    {
                        var refMatch = ReferenceNumberRegex().Match(currentValue.ToString());
                        if (refMatch.Success)
                        {
                            refNum = int.Parse(refMatch.Groups[1].Value, CultureInfo.InvariantCulture);
                            if (refMatch.Groups[2].Success)
                                baseFrom = int.Parse(refMatch.Groups[2].Value, CultureInfo.InvariantCulture);
                            if (refMatch.Groups[3].Success)
                                baseTo = int.Parse(refMatch.Groups[3].Value, CultureInfo.InvariantCulture);
                        }
                    }
                }
                else if (!string.IsNullOrEmpty(currentField))
                {
                    currentValue.Append(' ');
                    currentValue.Append(trimmed);
                }
            }

            SaveRefField(ref authors, ref title, ref journal, ref pubmed, currentField, currentValue.ToString());

            if (refNum > 0 || !string.IsNullOrEmpty(title))
            {
                references.Add(new Reference(refNum, authors, title, journal, pubmed, baseFrom, baseTo));
            }
        }

        return references;
    }

    private static void SaveRefField(ref string authors, ref string title, ref string journal,
        ref string pubmed, string field, string value)
    {
        switch (field)
        {
            case "AUTHORS": authors = value.Trim(); break;
            case "TITLE": title = value.Trim(); break;
            case "JOURNAL": journal = value.Trim(); break;
            case "PUBMED": pubmed = value.Trim(); break;
        }
    }

    private static IReadOnlyList<Feature> ParseFeatures(string text)
    {
        if (string.IsNullOrWhiteSpace(text))
            return Array.Empty<Feature>();

        var features = new List<Feature>();
        var lines = text.Split('\n').Where(l => l.Length > 0).ToList();

        string currentKey = "";
        string currentLocation = "";
        var qualifiers = new Dictionary<string, string>();

        // Pending (possibly multi-line) qualifier. Continuation lines are buffered
        // raw, separated by '\n', and finalized only when the qualifier ends so the
        // outer quotes / embedded "" escapes are processed across the whole value.
        string? pendingQualName = null;
        var pendingQualValue = new StringBuilder();

        void FlushQualifier()
        {
            if (pendingQualName is null)
                return;
            qualifiers[pendingQualName] = FinalizeQualifierValue(pendingQualName, pendingQualValue.ToString());
            pendingQualName = null;
            pendingQualValue.Clear();
        }

        foreach (var line in lines)
        {
            // Skip header line
            if (line.TrimStart().StartsWith("Location/Qualifiers", StringComparison.OrdinalIgnoreCase))
                continue;

            // New feature (starts at column 5)
            if (line.Length > 5 && !char.IsWhiteSpace(line[5]) && char.IsWhiteSpace(line[0]))
            {
                FlushQualifier();

                // Save previous feature
                if (!string.IsNullOrEmpty(currentKey))
                {
                    features.Add(CreateFeature(currentKey, currentLocation, qualifiers));
                }

                var trimmed = line.Trim();
                var spaceIdx = trimmed.IndexOf(' ');
                if (spaceIdx > 0)
                {
                    currentKey = trimmed[..spaceIdx];
                    currentLocation = trimmed[(spaceIdx + 1)..].Trim();
                }
                else
                {
                    currentKey = trimmed;
                    currentLocation = "";
                }
                qualifiers = new Dictionary<string, string>();
            }
            // Qualifier (starts at column 21, with /)
            else if (line.Length > 21 && line.TrimStart().StartsWith('/'))
            {
                FlushQualifier();

                var qualLine = line.Trim()[1..]; // Remove leading /
                var eqIdx = qualLine.IndexOf('=');
                if (eqIdx > 0)
                {
                    pendingQualName = qualLine[..eqIdx];
                    pendingQualValue.Append(qualLine[(eqIdx + 1)..].Trim());
                }
                else
                {
                    qualifiers[qualLine] = "true";
                }
            }
            // Continuation of the current qualifier value
            else if (pendingQualName is not null && line.Length > 21)
            {
                // INSDC: continuation lines begin at column 22. Join with '\n' to mark
                // the wrap point; FinalizeQualifierValue converts it to a single space
                // (mirrors Biopython feature_qualifier newline->space handling).
                pendingQualValue.Append('\n').Append(line.Trim());
            }
        }

        FlushQualifier();

        // Save last feature
        if (!string.IsNullOrEmpty(currentKey))
        {
            features.Add(CreateFeature(currentKey, currentLocation, qualifiers));
        }

        return features;
    }

    // Qualifiers whose value is a contiguous biological sequence with no internal
    // whitespace: line wraps must be collapsed entirely (no inserted space), per
    // Biopython's _BaseGenBankConsumer.remove_space_keys.
    private static readonly HashSet<string> NoSpaceQualifierKeys =
        new(StringComparer.Ordinal) { "translation" };

    /// <summary>
    /// Finalizes a (possibly multi-line) raw qualifier value. Continuation wraps
    /// (encoded as '\n') become a single space; the outer quote pair is removed and
    /// INSDC-escaped embedded quotes ("") collapse to one; for sequence qualifiers
    /// such as /translation all whitespace is removed so the value reassembles exactly.
    /// </summary>
    private static string FinalizeQualifierValue(string qualName, string rawValue)
    {
        // Wrap points -> single space (Biopython feature_qualifier: newline -> space).
        var joined = rawValue.Replace("\n", " ");
        var unquoted = UnquoteQualifierValue(joined);

        if (NoSpaceQualifierKeys.Contains(qualName))
        {
            // /translation etc.: strip every whitespace char so the amino-acid string
            // is contiguous regardless of where the flat file wrapped it.
            return new string(unquoted.Where(c => !char.IsWhiteSpace(c)).ToArray());
        }

        return unquoted;
    }

    /// <summary>
    /// Unquotes a free-text qualifier value. When the value is enclosed in double
    /// quotes, the single outer pair is removed and INSDC-escaped embedded quotes
    /// (a literal '"' encoded as two consecutive quotes "") are collapsed to one.
    /// Unquoted values are returned unchanged.
    /// </summary>
    /// <remarks>
    /// Per the INSDC Feature Table Definition (qualifier value format): free text is
    /// surrounded by double quotation marks and an embedded quotation mark is
    /// represented by a pair of adjacent quotation marks.
    /// </remarks>
    private static string UnquoteQualifierValue(string value)
    {
        if (value.Length >= 2 && value[0] == '"' && value[^1] == '"')
        {
            return value[1..^1].Replace("\"\"", "\"");
        }

        return value;
    }

    private static Feature CreateFeature(string key, string rawLocation, Dictionary<string, string> qualifiers)
    {
        var location = ParseLocation(rawLocation);
        return new Feature(key, location, qualifiers);
    }

    /// <summary>
    /// Parses a feature location string.
    /// </summary>
    public static Location ParseLocation(string locationStr)
    {
        if (string.IsNullOrEmpty(locationStr))
            return new Location(0, 0, false, false, false, false, false,
                Array.Empty<(int, int)>(), locationStr);

        var (start, end, isComplement, isJoin, isOrder, is5PrimePartial, is3PrimePartial, parts) =
            SequenceFormatHelper.ParseLocationParts(locationStr, useStartsWithForComplement: true);
        return new Location(start, end, isComplement, isJoin, isOrder,
            is5PrimePartial, is3PrimePartial, parts, locationStr);
    }

    private static string ParseSequence(string text)
    {
        if (string.IsNullOrWhiteSpace(text))
            return "";

        var sb = new StringBuilder();
        foreach (var line in text.Split('\n'))
        {
            foreach (var c in line)
            {
                if (char.IsLetter(c))
                {
                    sb.Append(char.ToUpperInvariant(c));
                }
            }
        }
        return sb.ToString();
    }

    private static bool IsStandardField(string field)
    {
        return field is "LOCUS" or "DEFINITION" or "ACCESSION" or "VERSION" or
               "KEYWORDS" or "SOURCE" or "REFERENCE" or "FEATURES" or "ORIGIN" or
               "DBLINK" or "DBSOURCE" or "COMMENT";
    }

    #endregion

    #region Utility Methods

    /// <summary>
    /// Extracts a specific feature type from a GenBank record.
    /// </summary>
    public static IEnumerable<Feature> GetFeatures(GenBankRecord record, string featureKey)
    {
        return record.Features.Where(f =>
            f.Key.Equals(featureKey, StringComparison.OrdinalIgnoreCase));
    }

    /// <summary>
    /// Gets all CDS (coding sequence) features.
    /// </summary>
    public static IEnumerable<Feature> GetCDS(GenBankRecord record)
    {
        return GetFeatures(record, "CDS");
    }

    /// <summary>
    /// Gets all gene features.
    /// </summary>
    public static IEnumerable<Feature> GetGenes(GenBankRecord record)
    {
        return GetFeatures(record, "gene");
    }

    /// <summary>
    /// Extracts a subsequence based on a feature location.
    /// </summary>
    public static string ExtractSequence(GenBankRecord record, Location location)
        => FeatureLocationHelper.ExtractSequence(record.Sequence, location);

    /// <summary>
    /// Gets qualifier value from a feature.
    /// </summary>
    public static string? GetQualifier(Feature feature, string qualifierName)
    {
        return feature.Qualifiers.TryGetValue(qualifierName, out var value) ? value : null;
    }

    /// <summary>
    /// Translates a CDS feature to protein sequence.
    /// </summary>
    public static string? TranslateCDS(GenBankRecord record, Feature cds)
    {
        // Check if translation is already in qualifiers
        if (GetQualifier(cds, "translation") is { } existingTranslation)
            return existingTranslation;

        var dnaSeq = ExtractSequence(record, cds.Location);
        if (string.IsNullOrEmpty(dnaSeq))
            return null;

        // Simple translation
        return Translate(dnaSeq);
    }

    private static string Translate(string dna)
    {
        var sb = new StringBuilder();
        for (int i = 0; i + 2 < dna.Length; i += 3)
        {
            var codon = dna.Substring(i, 3);
            sb.Append(CodonToAminoAcid(codon));
        }
        return sb.ToString();
    }

    private static char CodonToAminoAcid(string codon)
    {
        return codon.ToUpperInvariant() switch
        {
            "TTT" or "TTC" => 'F',
            "TTA" or "TTG" or "CTT" or "CTC" or "CTA" or "CTG" => 'L',
            "ATT" or "ATC" or "ATA" => 'I',
            "ATG" => 'M',
            "GTT" or "GTC" or "GTA" or "GTG" => 'V',
            "TCT" or "TCC" or "TCA" or "TCG" or "AGT" or "AGC" => 'S',
            "CCT" or "CCC" or "CCA" or "CCG" => 'P',
            "ACT" or "ACC" or "ACA" or "ACG" => 'T',
            "GCT" or "GCC" or "GCA" or "GCG" => 'A',
            "TAT" or "TAC" => 'Y',
            "TAA" or "TAG" or "TGA" => '*',
            "CAT" or "CAC" => 'H',
            "CAA" or "CAG" => 'Q',
            "AAT" or "AAC" => 'N',
            "AAA" or "AAG" => 'K',
            "GAT" or "GAC" => 'D',
            "GAA" or "GAG" => 'E',
            "TGT" or "TGC" => 'C',
            "TGG" => 'W',
            "CGT" or "CGC" or "CGA" or "CGG" or "AGA" or "AGG" => 'R',
            "GGT" or "GGC" or "GGA" or "GGG" => 'G',
            _ => 'X'
        };
    }

    #endregion

    #region Regex Patterns

    [GeneratedRegex(@"^([A-Z]+)\s+")]
    private static partial Regex SectionHeaderRegex();

    [GeneratedRegex(@"(\d+)(?:\s+\(bases\s+(\d+)\s+to\s+(\d+)\))?")]
    private static partial Regex ReferenceNumberRegex();

    #endregion
}

