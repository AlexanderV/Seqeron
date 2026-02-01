using System;
using System.Collections.Generic;
using System.Text;

namespace Seqeron.Genomics
{
    /// <summary>
    /// Helper methods for extracting feature sequences from genomic records.
    /// </summary>
    public static class FeatureLocationHelper
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
            return isComplement ? new DnaSequence(seq).ReverseComplement().Sequence : seq;
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
    }
}
