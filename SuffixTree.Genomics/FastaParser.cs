using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;

namespace SuffixTree.Genomics
{
    /// <summary>
    /// Parser for FASTA format - the standard format for biological sequences.
    /// </summary>
    public static class FastaParser
    {
        /// <summary>
        /// Parses a FASTA string into DNA sequences.
        /// </summary>
        public static IEnumerable<FastaEntry> Parse(string fastaContent)
        {
            if (string.IsNullOrWhiteSpace(fastaContent))
                yield break;

            using var reader = new StringReader(fastaContent);
            foreach (var entry in ParseReader(reader))
            {
                yield return entry;
            }
        }

        /// <summary>
        /// Parses a FASTA file into DNA sequences.
        /// </summary>
        public static IEnumerable<FastaEntry> ParseFile(string filePath)
        {
            using var reader = new StreamReader(filePath);
            foreach (var entry in ParseReader(reader))
            {
                yield return entry;
            }
        }

        /// <summary>
        /// Reads sequences from a FASTA file asynchronously.
        /// </summary>
        public static async IAsyncEnumerable<FastaEntry> ParseFileAsync(string filePath)
        {
            using var reader = new StreamReader(filePath);
            string? header = null;
            var sequenceBuilder = new StringBuilder();

            string? line;
            while ((line = await reader.ReadLineAsync()) != null)
            {
                if (line.StartsWith('>'))
                {
                    if (header != null && sequenceBuilder.Length > 0)
                    {
                        yield return CreateEntry(header, sequenceBuilder.ToString());
                    }

                    header = line.Substring(1).Trim();
                    sequenceBuilder.Clear();
                }
                else if (!string.IsNullOrWhiteSpace(line))
                {
                    sequenceBuilder.Append(line.Trim());
                }
            }

            if (header != null && sequenceBuilder.Length > 0)
            {
                yield return CreateEntry(header, sequenceBuilder.ToString());
            }
        }

        /// <summary>
        /// Writes sequences to FASTA format.
        /// </summary>
        public static string ToFasta(IEnumerable<FastaEntry> entries, int lineWidth = 80)
        {
            var sb = new StringBuilder();
            foreach (var entry in entries)
            {
                sb.Append('>').AppendLine(entry.Header);

                string seq = entry.Sequence.Sequence;
                for (int i = 0; i < seq.Length; i += lineWidth)
                {
                    int len = Math.Min(lineWidth, seq.Length - i);
                    sb.AppendLine(seq.Substring(i, len));
                }
            }
            return sb.ToString();
        }

        /// <summary>
        /// Writes sequences to a FASTA file.
        /// </summary>
        public static void WriteFile(string filePath, IEnumerable<FastaEntry> entries, int lineWidth = 80)
        {
            File.WriteAllText(filePath, ToFasta(entries, lineWidth));
        }

        private static IEnumerable<FastaEntry> ParseReader(TextReader reader)
        {
            string? header = null;
            var sequenceBuilder = new StringBuilder();

            string? line;
            while ((line = reader.ReadLine()) != null)
            {
                if (line.StartsWith('>'))
                {
                    if (header != null && sequenceBuilder.Length > 0)
                    {
                        yield return CreateEntry(header, sequenceBuilder.ToString());
                    }

                    header = line.Substring(1).Trim();
                    sequenceBuilder.Clear();
                }
                else if (!string.IsNullOrWhiteSpace(line))
                {
                    sequenceBuilder.Append(line.Trim());
                }
            }

            if (header != null && sequenceBuilder.Length > 0)
            {
                yield return CreateEntry(header, sequenceBuilder.ToString());
            }
        }

        private static FastaEntry CreateEntry(string header, string sequence)
        {
            // Parse header: typically "ID description"
            var parts = header.Split(new[] { ' ', '\t' }, 2);
            string id = parts[0];
            string? description = parts.Length > 1 ? parts[1] : null;

            return new FastaEntry(id, description, new DnaSequence(sequence));
        }
    }

    /// <summary>
    /// Represents a single entry in a FASTA file.
    /// </summary>
    public sealed class FastaEntry
    {
        public FastaEntry(string id, string? description, DnaSequence sequence)
        {
            Id = id;
            Description = description;
            Sequence = sequence;
        }

        /// <summary>
        /// The sequence identifier (first word of header).
        /// </summary>
        public string Id { get; }

        /// <summary>
        /// Optional description (rest of header after ID).
        /// </summary>
        public string? Description { get; }

        /// <summary>
        /// The DNA sequence.
        /// </summary>
        public DnaSequence Sequence { get; }

        /// <summary>
        /// Full header line (ID + Description).
        /// </summary>
        public string Header => Description != null ? $"{Id} {Description}" : Id;

        public override string ToString() => $"{Id} ({Sequence.Length} bp)";
    }
}
