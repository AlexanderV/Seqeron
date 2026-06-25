using System;
using System.Collections.Generic;
using System.Globalization;
using System.Text;
using System.Text.RegularExpressions;

namespace Seqeron.Genomics.IO
{
    /// <summary>
    /// Helper methods for extracting feature sequences from genomic records.
    /// </summary>
    public static partial class FeatureLocationHelper
    {
        /// <summary>
        /// Extracts the sequence for a feature based on its location.
        /// Handles join locations and complement strands.
        /// </summary>
        /// <param name="fullSequence">The full genomic sequence.</param>
        /// <param name="location">The feature location.</param>
        /// <returns>The extracted sequence, reverse complemented if needed.</returns>
        private static string ExtractSequenceInternal(
            string fullSequence,
            IReadOnlyList<(int Start, int End)> parts,
            int start,
            int end,
            bool isComplement)
        {
            if (string.IsNullOrEmpty(fullSequence))
                return "";

            var sb = new StringBuilder();

            if (parts.Count > 0)
            {
                foreach (var (partStart, partEnd) in parts)
                {
                    var realStart = Math.Max(0, partStart - 1); // GenBank/EMBL are 1-based
                    var realEnd = Math.Min(fullSequence.Length, partEnd);
                    if (realStart < realEnd)
                    {
                        sb.Append(fullSequence[realStart..realEnd]);
                    }
                }
            }
            else if (start > 0 && end > 0)
            {
                var realStart = Math.Max(0, start - 1);
                var realEnd = Math.Min(fullSequence.Length, end);
                if (realStart < realEnd)
                {
                    sb.Append(fullSequence[realStart..realEnd]);
                }
            }

            var seq = sb.ToString();
            // Use static string-based reverse complement to support all IUPAC ambiguity codes.
            // DnaSequence constructor rejects non-ACGT characters, but real GenBank sequences
            // can contain IUPAC codes (N, R, Y, etc.) per INSDC 7.4.1.
            return isComplement ? DnaSequence.GetReverseComplementString(seq) : seq;
        }

        /// <summary>
        /// Extracts the sequence for a GenBank feature.
        /// </summary>
        public static string ExtractSequence(string fullSequence, GenBankParser.Location location)
            => ExtractSequenceInternal(fullSequence, location.Parts, location.Start, location.End, location.IsComplement);

        /// <summary>
        /// Extracts the sequence for an EMBL feature.
        /// </summary>
        public static string ExtractSequence(string fullSequence, EmblParser.Location location)
            => ExtractSequenceInternal(fullSequence, location.Parts, location.Start, location.End, location.IsComplement);

        #region Remote-aware assembly (INSDC FT 3.4 / 3.5)

        /// <summary>
        /// Resolves a remote-entry span for <see cref="ResolveLocationSequence"/>.
        /// </summary>
        /// <param name="accession">
        /// The remote entry accession (the part before the first '.'/':' in
        /// "accession[.version]:span"), e.g. <c>J00194</c>.
        /// </param>
        /// <param name="version">
        /// The sequence version following the accession ('.N'), or <c>null</c> when the
        /// reference carries no explicit version (INSDC FT 3.4.2.1(e)).
        /// </param>
        /// <returns>
        /// The FULL sequence (5'-to-3', base 1 first) of the referenced remote entry, or
        /// <c>null</c> when the caller cannot supply it. The library performs the 1-based
        /// inclusive slicing and any strand handling itself; the resolver only fetches.
        /// </returns>
        public delegate string? RemoteSequenceResolver(string accession, int? version);

        /// <summary>
        /// Assembles the FULL feature sequence for a parsed location, splicing local spans
        /// (taken from <paramref name="localSequence"/>) together with remote spans (fetched
        /// through the caller-supplied <paramref name="remoteResolver"/>), honouring the
        /// INSDC Feature Table location semantics: segment order under <c>join</c>/<c>order</c>,
        /// reverse-complement under <c>complement(...)</c>, 1-based inclusive coordinates, and
        /// the <c>&lt;</c>/<c>&gt;</c> partial markers.
        /// </summary>
        /// <remarks>
        /// <para>
        /// The library is offline-first and performs NO network I/O itself. Resolving the
        /// sequence of a remote entry referenced in a location is delegated to the caller via
        /// <paramref name="remoteResolver"/> (the caller does any network/database access);
        /// this method does the assembly. When the location contains no remote reference the
        /// resolver is never invoked and the result is identical to a local extraction.
        /// </para>
        /// <para>
        /// Assembly rules (INSDC FT Definition v11.x, §3.4–3.5):
        /// </para>
        /// <list type="bullet">
        /// <item><description>
        /// Base 1 is the first base (5' end); a span <c>n..m</c> selects bases n..m inclusive
        /// (1-based). A single base number <c>n</c> selects base n.
        /// </description></item>
        /// <item><description>
        /// <c>join(a,b,...)</c> places the elements end-to-end in the listed order;
        /// <c>order(a,b,...)</c> assembles the same ordered concatenation.
        /// </description></item>
        /// <item><description>
        /// <c>complement(location)</c> reads the complement of the span in its 5'-to-3'
        /// direction, i.e. the reverse complement. <c>complement(join(a,b))</c> joins the
        /// elements first and then reverse-complements the whole result — equivalent to
        /// <c>join(complement(b),complement(a))</c> (segment order is reversed).
        /// </description></item>
        /// <item><description>
        /// A remote element <c>accession[.version]:span</c> is sliced from the resolver-supplied
        /// remote sequence with the same 1-based inclusive convention.
        /// </description></item>
        /// <item><description>
        /// The <c>&lt;</c>/<c>&gt;</c> partial markers mean the true endpoint lies beyond the
        /// stated number; the stated number is the only available coordinate, so assembly uses
        /// it as written.
        /// </description></item>
        /// </list>
        /// <para>
        /// Missing-resolver / resolver-returns-null behaviour: when a remote element is present
        /// and <paramref name="remoteResolver"/> is <c>null</c>, or the resolver returns
        /// <c>null</c>/empty for that accession, the remote element contributes nothing (empty
        /// string) to the assembly — the local spans are still spliced in their correct
        /// positions. This mirrors the existing local extraction, which clamps out-of-range
        /// spans rather than throwing.
        /// </para>
        /// </remarks>
        /// <param name="rawLocation">
        /// The raw location string as parsed from the feature table, e.g.
        /// <c>join(1..10,J00194.1:5..14)</c> or <c>complement(join(1..5,X.1:1..4))</c>.
        /// </param>
        /// <param name="localSequence">The full sequence of the local entry (base 1 first).</param>
        /// <param name="remoteResolver">
        /// Caller-supplied delegate that returns the full sequence of a remote entry, or
        /// <c>null</c>. May be <c>null</c> when the caller knows the location is local-only.
        /// </param>
        /// <returns>The assembled feature sequence.</returns>
        public static string ResolveLocationSequence(
            string rawLocation,
            string localSequence,
            RemoteSequenceResolver? remoteResolver)
        {
            if (string.IsNullOrEmpty(rawLocation))
                return string.Empty;

            // INSDC FT 3.5.3: complement(...) applies to the WHOLE enclosed span. Detect the
            // top-level outer complement, then assemble the inner ordered concatenation and
            // reverse-complement the assembled string as a single unit. This realises
            // complement(join(a,b)) == join(complement(b),complement(a)) without having to
            // reverse the segment list explicitly — reverse-complementing the concatenation
            // both reverses order and complements each segment.
            string descriptor = rawLocation.Trim();
            bool outerComplement = false;
            var outer = OuterComplementRegex().Match(descriptor);
            if (outer.Success)
            {
                outerComplement = true;
                descriptor = outer.Groups["inner"].Value;
            }

            var assembled = new StringBuilder();
            foreach (var segment in EnumerateSegments(descriptor))
            {
                assembled.Append(ResolveSegment(segment, localSequence, remoteResolver));
            }

            var result = assembled.ToString();
            return outerComplement ? DnaSequence.GetReverseComplementString(result) : result;
        }

        /// <summary>
        /// Convenience overload that assembles directly from a parsed
        /// <see cref="EmblParser.Location"/> using its preserved <c>RawLocation</c>.
        /// </summary>
        public static string ResolveLocationSequence(
            EmblParser.Location location,
            string localSequence,
            RemoteSequenceResolver? remoteResolver)
            => ResolveLocationSequence(location.RawLocation, localSequence, remoteResolver);

        /// <summary>
        /// Splits a location descriptor (with any outer complement already stripped) into its
        /// top-level segments in listed order. A <c>join(...)</c>/<c>order(...)</c> wrapper is
        /// removed and its comma-separated elements are yielded individually; a bare span is
        /// yielded as the single segment. Splitting respects parenthesis nesting so a
        /// per-segment <c>complement(...)</c> stays intact.
        /// </summary>
        private static IEnumerable<string> EnumerateSegments(string descriptor)
        {
            descriptor = descriptor.Trim();

            var wrapper = JoinOrderWrapperRegex().Match(descriptor);
            string inner = wrapper.Success ? wrapper.Groups["inner"].Value : descriptor;

            int depth = 0, start = 0;
            for (int i = 0; i < inner.Length; i++)
            {
                char c = inner[i];
                if (c == '(') depth++;
                else if (c == ')') depth--;
                else if (c == ',' && depth == 0)
                {
                    yield return inner[start..i];
                    start = i + 1;
                }
            }

            if (start <= inner.Length)
                yield return inner[start..];
        }

        /// <summary>
        /// Resolves a single segment to its assembled sequence. A segment may itself be a
        /// per-element <c>complement(...)</c>, a remote <c>accession[.version]:span</c>, or a
        /// local span. Per-element complement reverse-complements only that element's slice.
        /// </summary>
        private static string ResolveSegment(
            string segment,
            string localSequence,
            RemoteSequenceResolver? remoteResolver)
        {
            segment = segment.Trim();
            if (segment.Length == 0)
                return string.Empty;

            bool complement = false;
            var seg = SegmentComplementRegex().Match(segment);
            if (seg.Success)
            {
                complement = true;
                segment = seg.Groups["inner"].Value.Trim();
            }

            // Remote element: "accession[.version]:span" — fetch via resolver, then slice.
            var remote = SegmentRemoteRegex().Match(segment);
            string source;
            string spanText;
            if (remote.Success)
            {
                string accession = remote.Groups["acc"].Value;
                int? version = remote.Groups["ver"].Success
                    ? int.Parse(remote.Groups["ver"].Value, CultureInfo.InvariantCulture)
                    : null;
                source = remoteResolver?.Invoke(accession, version) ?? string.Empty;
                spanText = segment[remote.Length..];
            }
            else
            {
                source = localSequence ?? string.Empty;
                spanText = segment;
            }

            string slice = SliceSpan(source, spanText);
            return complement ? DnaSequence.GetReverseComplementString(slice) : slice;
        }

        /// <summary>
        /// Slices a 1-based inclusive span ("n..m", "n", or single-base "n") from
        /// <paramref name="source"/>, clamping to the available bounds and ignoring the
        /// <c>&lt;</c>/<c>&gt;</c> partial markers (the stated number is used as written).
        /// </summary>
        private static string SliceSpan(string source, string spanText)
        {
            if (string.IsNullOrEmpty(source))
                return string.Empty;

            var match = SequenceFormatHelper.LocationRangeRegex().Match(spanText);
            if (!match.Success)
                return string.Empty;

            int start = int.Parse(match.Groups[1].Value, CultureInfo.InvariantCulture);
            int end = match.Groups[2].Success
                ? int.Parse(match.Groups[2].Value, CultureInfo.InvariantCulture)
                : start;

            int realStart = Math.Max(0, start - 1);      // 1-based → 0-based inclusive
            int realEnd = Math.Min(source.Length, end);   // inclusive end → exclusive slice bound
            if (realStart >= realEnd)
                return string.Empty;

            return source[realStart..realEnd];
        }

        // INSDC FT 3.5.3: an outer complement wrapping the whole location, e.g.
        // "complement(join(1..5,X.1:1..4))". The inner group is balanced by the trailing ')'.
        [GeneratedRegex(@"^complement\((?<inner>.*)\)$", RegexOptions.IgnoreCase)]
        private static partial Regex OuterComplementRegex();

        // INSDC FT 3.5.3: a join/order wrapper around the comma-separated elements.
        [GeneratedRegex(@"^(?:join|order)\((?<inner>.*)\)$", RegexOptions.IgnoreCase)]
        private static partial Regex JoinOrderWrapperRegex();

        // A per-element complement(...) inside a join/order, e.g. "complement(4918..5163)".
        [GeneratedRegex(@"^complement\((?<inner>.*)\)$", RegexOptions.IgnoreCase)]
        private static partial Regex SegmentComplementRegex();

        // INSDC FT 3.4.2.1(e): a remote-entry prefix "accession[.version]:" on a segment.
        [GeneratedRegex(@"^(?<acc>[A-Za-z][A-Za-z0-9_]*)(?:\.(?<ver>\d+))?:")]
        private static partial Regex SegmentRemoteRegex();

        #endregion
    }
}

