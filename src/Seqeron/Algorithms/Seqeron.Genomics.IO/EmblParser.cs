using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Text;
using System.Text.RegularExpressions;

namespace Seqeron.Genomics.IO;

/// <summary>
/// Parser for EMBL flat file format (.embl, .dat).
/// EMBL format uses two-letter line type prefixes.
/// </summary>
public static partial class EmblParser
{
    #region Records

    /// <summary>Represents a complete EMBL record</summary>
    public readonly record struct EmblRecord(
        string Accession,
        string SequenceVersion,
        string DataClass,
        string MoleculeType,
        string Topology,
        string TaxonomicDivision,
        int SequenceLength,
        string Description,
        IReadOnlyList<string> Keywords,
        string Organism,
        IReadOnlyList<string> OrganismClassification,
        IReadOnlyList<Reference> References,
        IReadOnlyList<Feature> Features,
        string Sequence,
        IReadOnlyDictionary<string, string> AdditionalFields);

    /// <summary>Literature reference</summary>
    public readonly record struct Reference(
        int Number,
        string Citation,
        string Authors,
        string Title,
        string Journal,
        string CrossReference,
        string Comment,
        string Positions,
        string Group);

    /// <summary>Sequence feature with location and qualifiers</summary>
    public readonly record struct Feature(
        string Key,
        Location Location,
        IReadOnlyDictionary<string, string> Qualifiers);

    /// <summary>
    /// A single remote-entry segment of a feature location: a span that points into
    /// another INSDC entry's sequence (INSDC Feature Table 3.4.2.1(e),
    /// <c>accession[.version]:descriptor</c>). Used when a remote reference appears
    /// nested inside a <c>join</c>/<c>order</c>/<c>complement</c> operator, e.g. the span
    /// <c>J00194.1:100..202</c> inside <c>join(1..100,J00194.1:100..202)</c>.
    /// </summary>
    public readonly record struct RemotePart(string Accession, string? Version, int Start, int End);

    /// <summary>Feature location</summary>
    /// <remarks>
    /// The trailing members (<see cref="IsBetween"/>, <see cref="IsSingleBaseFromRange"/>,
    /// <see cref="RemoteAccession"/>, <see cref="RemoteVersion"/>, <see cref="RemoteParts"/>)
    /// support the rarer INSDC Feature Table location forms (section 3.4.2.1 / 3.4.3): the
    /// site-between operator <c>n^m</c>, the deprecated single-base-from-range <c>n.m</c>,
    /// top-level remote references <c>accession[.version]:descriptor</c>, and remote
    /// references nested inside a <c>join</c>/<c>order</c>/<c>complement</c> operator
    /// (captured per-segment in <see cref="RemoteParts"/>). They default to
    /// <c>false</c>/<c>null</c>/empty so ordinary local spans are reported exactly as before.
    /// </remarks>
    public readonly record struct Location(
        int Start,
        int End,
        bool IsComplement,
        bool IsJoin,
        bool IsOrder,
        bool Is5PrimePartial,
        bool Is3PrimePartial,
        IReadOnlyList<(int Start, int End)> Parts,
        string RawLocation,
        bool IsBetween = false,
        bool IsSingleBaseFromRange = false,
        string? RemoteAccession = null,
        string? RemoteVersion = null,
        IReadOnlyList<RemotePart>? RemoteParts = null)
    {
        /// <summary>
        /// True when this location points into another INSDC entry — either a top-level
        /// remote reference (such as <c>J00194.1:100..202</c>) or a remote reference
        /// nested inside an operator (captured in <see cref="RemoteParts"/>).
        /// </summary>
        public bool IsRemote => RemoteAccession is not null || (RemoteParts is { Count: > 0 });
    }

    #endregion

    #region Line Type Prefixes

    // Standard EMBL line types
    private const string ID = "ID";   // Identification
    private const string AC = "AC";   // Accession number
    private const string SV = "SV";   // Sequence version
    private const string DT = "DT";   // Date
    private const string DE = "DE";   // Description
    private const string KW = "KW";   // Keywords
    private const string OS = "OS";   // Organism species
    private const string OC = "OC";   // Organism classification
    private const string OG = "OG";   // Organelle
    private const string RN = "RN";   // Reference number
    private const string RC = "RC";   // Reference comment
    private const string RP = "RP";   // Reference positions
    private const string RX = "RX";   // Reference cross-reference
    private const string RG = "RG";   // Reference group
    private const string RA = "RA";   // Reference authors
    private const string RT = "RT";   // Reference title
    private const string RL = "RL";   // Reference location
    private const string DR = "DR";   // Database cross-reference
    private const string CC = "CC";   // Comments
    private const string FH = "FH";   // Feature header
    private const string FT = "FT";   // Feature table
    private const string SQ = "SQ";   // Sequence header
    private const string XX = "XX";   // Spacer line
    private const string END = "//";  // Entry terminator

    #endregion

    #region Controlled Vocabularies

    // INSDC Feature Table v11.3 /mol_type qualifier controlled vocabulary.
    private static readonly HashSet<string> ValidMolTypes = new(StringComparer.Ordinal)
    {
        "genomic DNA", "genomic RNA", "mRNA", "tRNA", "rRNA",
        "other RNA", "other DNA", "transcribed RNA", "viral cRNA",
        "unassigned DNA", "unassigned RNA"
    };

    // EBI EMBL User Manual Release 143, Section 3.1.
    private static readonly HashSet<string> ValidDataClasses = new(StringComparer.Ordinal)
    {
        "CON", "PAT", "EST", "GSS", "HTC", "HTG", "WGS", "TSA", "STS", "STD"
    };

    // EBI EMBL User Manual Release 143, Section 3.2.
    private static readonly HashSet<string> ValidDivisions = new(StringComparer.Ordinal)
    {
        "PHG", "ENV", "FUN", "HUM", "INV", "MAM", "VRT",
        "MUS", "PLN", "PRO", "ROD", "SYN", "TGN", "UNC", "VRL"
    };

    #endregion

    #region Main Parsing Methods

    /// <summary>
    /// Parses EMBL records from a file.
    /// </summary>
    public static IEnumerable<EmblRecord> ParseFile(string filePath)
    {
        if (string.IsNullOrEmpty(filePath) || !File.Exists(filePath))
            yield break;

        foreach (var record in Parse(File.ReadAllText(filePath)))
        {
            yield return record;
        }
    }

    /// <summary>
    /// Parses EMBL records from text content.
    /// Multiple records are separated by // delimiter.
    /// </summary>
    public static IEnumerable<EmblRecord> Parse(string content)
    {
        if (string.IsNullOrEmpty(content))
            yield break;

        // Split by record delimiter
        var recordTexts = content.Split(new[] { "\n//" }, StringSplitOptions.RemoveEmptyEntries);

        foreach (var recordText in recordTexts)
        {
            var trimmed = recordText.Trim();
            if (!string.IsNullOrEmpty(trimmed) && trimmed.StartsWith("ID", StringComparison.Ordinal))
            {
                var record = ParseRecord(trimmed);
                if (record.HasValue)
                    yield return record.Value;
            }
        }
    }

    /// <summary>
    /// Parses a single EMBL record.
    /// </summary>
    private static EmblRecord? ParseRecord(string text)
    {
        var lines = text.Split('\n')
            .Select(l => l.TrimEnd('\r'))
            .Where(l => !string.IsNullOrEmpty(l) && l.Length >= 2)
            .ToList();

        // Extract line groups by prefix
        var lineGroups = GroupLinesByPrefix(lines);

        // Parse ID line
        var (accession, version, dataClass, moleculeType, topology, division, seqLength) =
            ParseIdLine(GetFirstLine(lineGroups, ID));

        // Parse version from separate SV line if not in ID line
        if (string.IsNullOrEmpty(version))
        {
            version = ParseAccessionLine(GetFirstLine(lineGroups, SV));
        }
        if (string.IsNullOrEmpty(accession))
        {
            accession = ParseAccessionLine(GetLines(lineGroups, AC).FirstOrDefault() ?? "");
        }

        // Description
        var description = JoinLines(GetLines(lineGroups, DE));

        // Keywords
        var keywords = ParseKeywords(JoinLines(GetLines(lineGroups, KW)));

        // Organism
        var organism = JoinLines(GetLines(lineGroups, OS)).TrimEnd('.');
        var classification = ParseClassification(JoinLines(GetLines(lineGroups, OC)));

        // References
        var references = ParseReferences(lines);

        // Features - parse directly from lines, not via GroupLinesByPrefix
        var features = ParseFeaturesFromLines(lines);

        // Sequence
        var sequence = ParseSequence(lines);

        // Additional fields — store both non-standard prefixes and standard prefixes
        // that are not consumed by dedicated parsers above (DT, DR, CC, OG).
        var consumedPrefixes = new HashSet<string>(StringComparer.Ordinal)
        {
            ID, AC, SV, DE, KW, OS, OC, RN, RC, RP, RX, RG, RA, RT, RL, FH, FT, SQ, XX
        };
        var additionalFields = new Dictionary<string, string>();
        foreach (var (prefix, content) in lineGroups)
        {
            if (!consumedPrefixes.Contains(prefix))
            {
                additionalFields[prefix] = content;
            }
        }

        return new EmblRecord(
            accession,
            version,
            dataClass,
            moleculeType,
            topology,
            division,
            seqLength,
            description,
            keywords,
            organism,
            classification,
            references,
            features,
            sequence,
            additionalFields
        );
    }

    #endregion

    #region Line Parsing Helpers

    private static Dictionary<string, string> GroupLinesByPrefix(List<string> lines)
    {
        var groups = new Dictionary<string, StringBuilder>();

        foreach (var line in lines)
        {
            if (line.Length < 2) continue;

            var prefix = line.Length >= 2 ? line[..2] : line;
            var content = line.Length > 5 ? line[5..] : "";

            if (!groups.ContainsKey(prefix))
                groups[prefix] = new StringBuilder();

            if (groups[prefix].Length > 0)
                groups[prefix].Append(' ');
            groups[prefix].Append(content.Trim());
        }

        return groups.ToDictionary(kv => kv.Key, kv => kv.Value.ToString());
    }

    private static string GetFirstLine(Dictionary<string, string> groups, string prefix)
    {
        return groups.GetValueOrDefault(prefix, "");
    }

    private static IEnumerable<string> GetLines(Dictionary<string, string> groups, string prefix)
    {
        if (groups.TryGetValue(prefix, out var content))
            yield return content;
    }

    private static string JoinLines(IEnumerable<string> lines)
    {
        return string.Join(" ", lines).Trim();
    }

    private static (string Accession, string Version, string DataClass, string MoleculeType, string Topology,
        string Division, int Length) ParseIdLine(string line)
    {
        if (string.IsNullOrEmpty(line))
            return ("", "", "", "", "", "", 0);

        // Format: ACCESSION; SV VERSION; TOPOLOGY; MOLECULE; DATA_CLASS; DIVISION; LENGTH BP.
        var parts = line.Split(';').Select(p => p.Trim()).ToArray();

        string accession = parts.Length > 0 ? parts[0].Split(' ')[0] : "";
        string version = "";
        string dataClass = "";
        string moleculeType = "";
        string topology = "";
        string division = "";
        int length = 0;

        foreach (var part in parts.Skip(1))
        {
            var trimmed = part.Trim();
            if (trimmed.EndsWith("BP", StringComparison.OrdinalIgnoreCase) ||
                trimmed.EndsWith("BP.", StringComparison.OrdinalIgnoreCase))
            {
                var lengthMatch = LengthRegex().Match(trimmed);
                if (lengthMatch.Success)
                    length = int.Parse(lengthMatch.Groups[1].Value, CultureInfo.InvariantCulture);
            }
            else if (trimmed.StartsWith("SV ", StringComparison.OrdinalIgnoreCase))
            {
                version = trimmed[3..].Trim();
            }
            else if (trimmed is "linear" or "circular")
            {
                topology = trimmed;
            }
            else if (ValidMolTypes.Contains(trimmed))
            {
                moleculeType = trimmed;
            }
            else if (ValidDataClasses.Contains(trimmed))
            {
                dataClass = trimmed;
            }
            else if (ValidDivisions.Contains(trimmed))
            {
                division = trimmed;
            }
        }

        return (accession, version, dataClass, moleculeType, topology, division, length);
    }

    private static string ParseAccessionLine(string line)
    {
        if (string.IsNullOrEmpty(line))
            return "";

        // Take first accession number (primary)
        return line.Split(new[] { ';', ' ' }, StringSplitOptions.RemoveEmptyEntries)
            .FirstOrDefault()?.Trim() ?? "";
    }

    private static IReadOnlyList<string> ParseKeywords(string text)
    {
        if (string.IsNullOrWhiteSpace(text) || text.Trim() == ".")
            return Array.Empty<string>();

        return text
            .Split(new[] { ';' }, StringSplitOptions.RemoveEmptyEntries)
            .Select(k => k.Trim().TrimEnd('.'))
            .Where(k => !string.IsNullOrEmpty(k))
            .ToList();
    }

    private static IReadOnlyList<string> ParseClassification(string text)
    {
        if (string.IsNullOrWhiteSpace(text))
            return Array.Empty<string>();

        return text
            .Split(new[] { ';' }, StringSplitOptions.RemoveEmptyEntries)
            .Select(c => c.Trim().TrimEnd('.'))
            .Where(c => !string.IsNullOrEmpty(c))
            .ToList();
    }

    private static IReadOnlyList<Reference> ParseReferences(List<string> lines)
    {
        var references = new List<Reference>();

        int currentRefNum = 0;
        string citation = "", authors = "", title = "", journal = "", xref = "", comment = "";
        string positions = "", group = "";

        // Normalize and append a finished reference. Applied identically whether the block is
        // flushed on encountering the next RN line or at end-of-record, so every reference in a
        // multi-reference entry has the same trailing-separator trimming.
        void AddReference()
        {
            references.Add(new Reference(currentRefNum, citation, authors.TrimEnd(';', ' '),
                title.Trim('"', ';', ' '), journal.TrimEnd(';', ' '), xref, comment,
                positions, group.TrimEnd(';', ' ')));
        }

        foreach (var line in lines)
        {
            if (line.Length < 5) continue;

            var prefix = line[..2];
            var content = line.Length > 5 ? line[5..].Trim() : "";

            switch (prefix)
            {
                case RN:
                    // Save previous reference
                    if (currentRefNum > 0)
                    {
                        AddReference();
                    }
                    // Parse new reference number
                    var numMatch = ReferenceNumberRegex().Match(content);
                    currentRefNum = numMatch.Success ? int.Parse(numMatch.Groups[1].Value, CultureInfo.InvariantCulture) : 0;
                    citation = authors = title = journal = xref = comment = "";
                    positions = group = "";
                    break;
                case RC:
                    comment = string.IsNullOrEmpty(comment) ? content : comment + " " + content;
                    break;
                case RA:
                    authors = string.IsNullOrEmpty(authors) ? content : authors + " " + content;
                    break;
                case RT:
                    title = string.IsNullOrEmpty(title) ? content : title + " " + content;
                    break;
                case RL:
                    journal = string.IsNullOrEmpty(journal) ? content : journal + " " + content;
                    break;
                case RX:
                    xref = string.IsNullOrEmpty(xref) ? content : xref + "; " + content;
                    break;
                case RP:
                    positions = string.IsNullOrEmpty(positions) ? content : positions + " " + content;
                    break;
                case RG:
                    group = string.IsNullOrEmpty(group) ? content : group + " " + content;
                    break;
            }
        }

        // Save last reference
        if (currentRefNum > 0)
        {
            AddReference();
        }

        return references;
    }

    /// <summary>
    /// Parses features directly from FT lines, preserving whitespace structure.
    /// EMBL FT format:
    /// FT   key             location
    /// FT                   /qualifier="value"
    /// </summary>
    private static IReadOnlyList<Feature> ParseFeaturesFromLines(List<string> lines)
    {
        var features = new List<Feature>();

        string? currentKey = null;
        string currentLocation = "";
        var currentQualifiers = new Dictionary<string, string>();
        var qualifierBuilder = new StringBuilder();
        string? currentQualifierName = null;

        foreach (var line in lines)
        {
            if (!line.StartsWith("FT", StringComparison.Ordinal) || line.Length < 5)
                continue;

            // FT lines: positions 0-1 = "FT", 2-4 = spaces, 5+ = content
            // Feature key starts at position 5, qualifiers are indented further
            var content = line.Length > 5 ? line[5..] : "";

            if (string.IsNullOrWhiteSpace(content))
                continue;

            // Check if this is a new feature (key starts at column 5, not indented)
            // or a continuation (qualifiers/location continuation are indented)
            bool isNewFeature = content.Length > 0 && !char.IsWhiteSpace(content[0]) && !content.TrimStart().StartsWith('/');
            bool isQualifier = content.TrimStart().StartsWith('/');

            if (isNewFeature)
            {
                // Save previous feature
                if (currentKey != null)
                {
                    // Finish any pending qualifier
                    FinishQualifier(currentQualifierName, qualifierBuilder, currentQualifiers);
                    features.Add(CreateFeature(currentKey, currentLocation, currentQualifiers));
                }

                // Parse new feature: key and location
                var trimmed = content.Trim();
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
                currentQualifiers = new Dictionary<string, string>();
                qualifierBuilder.Clear();
                currentQualifierName = null;
            }
            else if (isQualifier)
            {
                // Finish previous qualifier if any
                FinishQualifier(currentQualifierName, qualifierBuilder, currentQualifiers);

                // Parse new qualifier: /name=value or /name (boolean)
                // Manual parsing instead of regex to correctly handle values containing '/'.
                var qualContent = content.TrimStart();
                if (qualContent.StartsWith('/'))
                {
                    var eqIdx = qualContent.IndexOf('=', 1);
                    if (eqIdx > 1)
                    {
                        currentQualifierName = qualContent[1..eqIdx];
                        qualifierBuilder.Clear();
                        qualifierBuilder.Append(qualContent[(eqIdx + 1)..]);
                    }
                    else
                    {
                        currentQualifierName = qualContent[1..].Trim();
                        qualifierBuilder.Clear();
                        qualifierBuilder.Append("true");
                    }
                }
            }
            else
            {
                // Continuation line - could be location or qualifier value continuation
                var trimmed = content.Trim();
                if (currentQualifierName != null)
                {
                    // Qualifier value continuation
                    qualifierBuilder.Append(trimmed);
                }
                else if (currentKey != null)
                {
                    // Location continuation
                    currentLocation += trimmed;
                }
            }
        }

        // Save last feature
        if (currentKey != null)
        {
            FinishQualifier(currentQualifierName, qualifierBuilder, currentQualifiers);
            features.Add(CreateFeature(currentKey, currentLocation, currentQualifiers));
        }

        return features;
    }

    private static void FinishQualifier(string? name, StringBuilder value, Dictionary<string, string> qualifiers)
    {
        if (name != null)
        {
            qualifiers[name] = UnquoteQualifierValue(value.ToString().Trim());
        }
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

        // INSDC FT 3.4.2.1(e): a remote reference is a remote entry identifier
        // ("accession[.version]") followed by ':' and a local location descriptor that
        // applies to that entry's sequence (e.g. "J00194.1:100..202"). Strip and capture
        // the prefix so the accession-version digits do not leak into the parsed span;
        // the trailing descriptor is then parsed exactly like a local location.
        string? remoteAccession = null;
        string? remoteVersion = null;
        string descriptor = locationStr;
        var remoteMatch = RemoteReferenceRegex().Match(locationStr);
        if (remoteMatch.Success)
        {
            remoteAccession = remoteMatch.Groups["acc"].Value;
            remoteVersion = remoteMatch.Groups["ver"].Success ? remoteMatch.Groups["ver"].Value : null;
            descriptor = locationStr[remoteMatch.Length..];
        }

        // INSDC FT 3.4.2.1(b)/(c): the site-between operator "n^m" (a site between two
        // adjoining bases) and the deprecated single-base-from-range "n.m" (a single
        // base somewhere in [n,m]) are both bare two-number forms with a distinct
        // separator. Detect them before delegating, since the shared range regex only
        // understands the two-period span "n..m".
        var betweenMatch = SiteBetweenRegex().Match(descriptor);
        if (betweenMatch.Success)
        {
            int bStart = int.Parse(betweenMatch.Groups[1].Value, CultureInfo.InvariantCulture);
            int bEnd = int.Parse(betweenMatch.Groups[2].Value, CultureInfo.InvariantCulture);
            return new Location(bStart, bEnd, false, false, false, false, false,
                new[] { (bStart, bEnd) }, locationStr,
                IsBetween: true, RemoteAccession: remoteAccession, RemoteVersion: remoteVersion);
        }

        var singleDotMatch = SingleBaseFromRangeRegex().Match(descriptor);
        if (singleDotMatch.Success)
        {
            int sStart = int.Parse(singleDotMatch.Groups[1].Value, CultureInfo.InvariantCulture);
            int sEnd = int.Parse(singleDotMatch.Groups[2].Value, CultureInfo.InvariantCulture);
            return new Location(sStart, sEnd, false, false, false, false, false,
                new[] { (sStart, sEnd) }, locationStr,
                IsSingleBaseFromRange: true, RemoteAccession: remoteAccession, RemoteVersion: remoteVersion);
        }

        // INSDC FT 3.4.2.1(e) + 3.4.3: a remote reference may appear NESTED inside a
        // join/order/complement operator, e.g. "join(1..100,J00194.1:100..202)" — "Joins
        // region 1..100 of the existing entry with the region 100..202 of remote entry
        // J00194". The top-level prefix strip above is anchored, so it does not catch these.
        // Scan the descriptor for embedded "accession[.version]:" prefixes, capture each
        // segment's remote entry per-part, and REMOVE the prefixes before delegating so the
        // accession-version digits cannot leak into the local numeric span (the shared range
        // regex would otherwise read the version, e.g. ".1", as a spurious single-base part).
        var (cleanedDescriptor, remoteParts) = ExtractNestedRemoteReferences(descriptor);

        var (start, end, isComplement, isJoin, isOrder, is5PrimePartial, is3PrimePartial, parts) =
            SequenceFormatHelper.ParseLocationParts(cleanedDescriptor, useStartsWithForComplement: false);
        return new Location(start, end, isComplement, isJoin, isOrder,
            is5PrimePartial, is3PrimePartial, parts, locationStr,
            RemoteAccession: remoteAccession, RemoteVersion: remoteVersion,
            RemoteParts: remoteParts);
    }

    /// <summary>
    /// Finds remote-entry references nested inside an operator location descriptor
    /// (INSDC FT 3.4.2.1(e) / 3.4.3, e.g. <c>J00194.1:100..202</c> inside
    /// <c>join(1..100,J00194.1:100..202)</c>). For each embedded
    /// <c>accession[.version]:span</c> it records a <see cref="RemotePart"/> with the
    /// remote span's bounds, and returns the descriptor with the remote prefixes removed so
    /// the accession-version digits do not leak into the local numeric parse.
    /// </summary>
    /// <returns>
    /// The descriptor with embedded <c>accession[.version]:</c> prefixes stripped, and the
    /// list of captured remote segments (<c>null</c> when none were found, so ordinary
    /// local operator forms are reported exactly as before).
    /// </returns>
    private static (string Cleaned, IReadOnlyList<RemotePart>? RemoteParts)
        ExtractNestedRemoteReferences(string descriptor)
    {
        var matches = NestedRemoteReferenceRegex().Matches(descriptor);
        if (matches.Count == 0)
            return (descriptor, null);

        var remoteParts = new List<RemotePart>();
        var cleaned = new StringBuilder(descriptor.Length);
        int cursor = 0;

        foreach (Match match in matches)
        {
            // Append the text before this prefix unchanged.
            cleaned.Append(descriptor, cursor, match.Index - cursor);
            cursor = match.Index + match.Length;

            string accession = match.Groups["acc"].Value;
            string? version = match.Groups["ver"].Success ? match.Groups["ver"].Value : null;

            // The remote span is the bare range descriptor that immediately follows the
            // stripped "accession[.version]:" prefix, up to the next separator (',' / ')').
            int spanEnd = cursor;
            while (spanEnd < descriptor.Length && descriptor[spanEnd] != ',' && descriptor[spanEnd] != ')')
                spanEnd++;
            string span = descriptor[cursor..spanEnd];

            var spanMatch = SequenceFormatHelper.LocationRangeRegex().Match(span);
            int rStart = 0, rEnd = 0;
            if (spanMatch.Success)
            {
                rStart = int.Parse(spanMatch.Groups[1].Value, CultureInfo.InvariantCulture);
                rEnd = spanMatch.Groups[2].Success
                    ? int.Parse(spanMatch.Groups[2].Value, CultureInfo.InvariantCulture)
                    : rStart;
            }

            remoteParts.Add(new RemotePart(accession, version, rStart, rEnd));
        }

        // Append the remaining tail after the last prefix.
        cleaned.Append(descriptor, cursor, descriptor.Length - cursor);

        return (cleaned.ToString(), remoteParts);
    }

    private static string ParseSequence(List<string> lines)
    {
        var sb = new StringBuilder();
        bool inSequence = false;

        foreach (var line in lines)
        {
            if (line.StartsWith("SQ", StringComparison.Ordinal))
            {
                inSequence = true;
                continue;
            }

            if (inSequence && !line.StartsWith("//", StringComparison.Ordinal))
            {
                foreach (var c in line)
                {
                    if (char.IsLetter(c))
                    {
                        sb.Append(char.ToUpperInvariant(c));
                    }
                }
            }
        }

        return sb.ToString();
    }

    #endregion

    #region Utility Methods

    /// <summary>
    /// Converts an EMBL record to GenBank format.
    /// </summary>
    public static GenBankParser.GenBankRecord ToGenBank(EmblRecord embl)
    {
        var gbFeatures = embl.Features.Select(f =>
            new GenBankParser.Feature(
                f.Key,
                new GenBankParser.Location(f.Location.Start, f.Location.End, f.Location.IsComplement,
                    f.Location.IsJoin, f.Location.IsOrder,
                    f.Location.Is5PrimePartial, f.Location.Is3PrimePartial,
                    f.Location.Parts, f.Location.RawLocation),
                f.Qualifiers
            )).ToList();

        var gbReferences = embl.References.Select(r =>
            new GenBankParser.Reference(r.Number, r.Authors, r.Title, r.Journal, r.CrossReference, null, null)
        ).ToList();

        return new GenBankParser.GenBankRecord(
            embl.Accession,
            embl.SequenceLength,
            embl.MoleculeType,
            embl.Topology,
            embl.TaxonomicDivision,
            null,
            embl.Description,
            embl.Accession,
            embl.SequenceVersion,
            embl.Keywords,
            embl.Organism,
            string.Join("; ", embl.OrganismClassification),
            gbReferences,
            gbFeatures,
            embl.Sequence,
            embl.AdditionalFields
        );
    }

    /// <summary>
    /// Extracts features by type.
    /// </summary>
    public static IEnumerable<Feature> GetFeatures(EmblRecord record, string featureKey)
    {
        return record.Features.Where(f =>
            f.Key.Equals(featureKey, StringComparison.OrdinalIgnoreCase));
    }

    /// <summary>
    /// Gets CDS features.
    /// </summary>
    public static IEnumerable<Feature> GetCDS(EmblRecord record)
    {
        return GetFeatures(record, "CDS");
    }

    /// <summary>
    /// Gets gene features.
    /// </summary>
    public static IEnumerable<Feature> GetGenes(EmblRecord record)
    {
        return GetFeatures(record, "gene");
    }

    /// <summary>
    /// Extracts subsequence based on location.
    /// </summary>
    public static string ExtractSequence(EmblRecord record, Location location)
        => FeatureLocationHelper.ExtractSequence(record.Sequence, location);

    #endregion

    #region Regex Patterns

    [GeneratedRegex(@"(\d+)\s*BP", RegexOptions.IgnoreCase)]
    private static partial Regex LengthRegex();

    [GeneratedRegex(@"\[(\d+)\]")]
    private static partial Regex ReferenceNumberRegex();

    // INSDC FT 3.4.2.1(e): remote entry identifier "accession[.version]" then ':'.
    [GeneratedRegex(@"^(?<acc>[A-Za-z][A-Za-z0-9_]*)(?:\.(?<ver>\d+))?:")]
    private static partial Regex RemoteReferenceRegex();

    // INSDC FT 3.4.2.1(e) / 3.4.3: a remote entry identifier "accession[.version]:" that
    // appears NESTED inside a join/order/complement operator (e.g.
    // "join(1..100,J00194.1:100..202)"). Anchored via a lookbehind to an operator boundary —
    // the accession must be immediately preceded by '(' (operator open paren, including a
    // 'complement(' that wraps the individual segment) or ',' (a join/order element
    // separator) — so it only matches segment-leading remote prefixes and never a position
    // number inside a span. The lookbehind is NOT consumed, so any wrapping 'complement(' is
    // left intact in the descriptor for the shared complement detector.
    [GeneratedRegex(@"(?<=[(,])(?<acc>[A-Za-z][A-Za-z0-9_]*)(?:\.(?<ver>\d+))?:")]
    private static partial Regex NestedRemoteReferenceRegex();

    // INSDC FT 3.4.2.1(b): site between two adjoining bases, "n^m" (e.g. 123^124).
    [GeneratedRegex(@"^(\d+)\^(\d+)$")]
    private static partial Regex SiteBetweenRegex();

    // INSDC FT 3.4.2.1(c): deprecated single base from a range, "n.m" (single period,
    // e.g. 102.110) — distinct from the two-period sequence span "n..m".
    [GeneratedRegex(@"^(\d+)\.(\d+)$")]
    private static partial Regex SingleBaseFromRangeRegex();

    #endregion
}

