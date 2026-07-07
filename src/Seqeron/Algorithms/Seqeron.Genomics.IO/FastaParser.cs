using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using Seqeron.Genomics.Core;

namespace Seqeron.Genomics.IO
{
    /// <summary>
    /// Selectable sequence alphabet for FASTA parsing.
    /// </summary>
    /// <remarks>
    /// Alphabets are defined by their authoritative single-letter code sets:
    /// <list type="bullet">
    /// <item><see cref="StrictDna"/> — A/C/G/T only (default, unchanged historical behaviour).</item>
    /// <item><see cref="IupacNucleotide"/> — NC-IUB (1985) "Nomenclature for Incompletely Specified
    /// Bases in Nucleic Acid Sequences", Nucleic Acids Research 13(9):3021–3030: A C G T U R Y S W K M
    /// B D H V N plus gap '-'.</item>
    /// <item><see cref="Rna"/> — A/C/G/U only.</item>
    /// <item><see cref="Protein"/> — IUPAC amino-acid one-letter codes: 20 standard residues plus
    /// ambiguity codes B (Asx), Z (Glx), J (Xle), X (Xaa), the rare residues U (Sec) and O (Pyl),
    /// and stop '*'.</item>
    /// </list>
    /// </remarks>
    public enum SequenceAlphabet
    {
        /// <summary>Strict DNA: only A, C, G, T (default; matches the historical DNA-only parser).</summary>
        StrictDna = 0,

        /// <summary>IUPAC nucleotide codes (A C G T U R Y S W K M B D H V N) plus the gap symbol '-'.</summary>
        IupacNucleotide = 1,

        /// <summary>RNA: only A, C, G, U.</summary>
        Rna = 2,

        /// <summary>IUPAC amino-acid one-letter codes (20 standard + B Z J X U O) plus stop '*'.</summary>
        Protein = 3
    }

    /// <summary>
    /// Parser for FASTA format - the standard format for biological sequences.
    /// </summary>
    /// <remarks>
    /// The default <see cref="Parse(string)"/> / <see cref="ParseFile(string)"/> /
    /// <see cref="ParseFileAsync(string)"/> path is strict DNA-only and unchanged: it returns
    /// <see cref="FastaEntry"/> instances backed by <see cref="DnaSequence"/> (A/C/G/T only).
    /// The opt-in <see cref="SequenceAlphabet"/> overloads return <see cref="FastaRecord"/> and
    /// accept RNA, protein, or IUPAC-ambiguous nucleotide FASTA.
    /// </remarks>
    public static class FastaParser
    {
        // --- Authoritative alphabet code sets ---------------------------------------------------

        // IUPAC nucleotide codes per NC-IUB (1985), Nucleic Acids Research 13(9):3021–3030:
        // bases A C G T U; ambiguity R Y S W K M B D H V N; gap '-'.
        private static readonly HashSet<char> IupacNucleotideCodes = new()
        {
            'A', 'C', 'G', 'T', 'U',
            'R', 'Y', 'S', 'W', 'K', 'M',
            'B', 'D', 'H', 'V', 'N',
            '-'
        };

        // RNA bases.
        private static readonly HashSet<char> RnaCodes = new() { 'A', 'C', 'G', 'U' };

        // IUPAC amino-acid one-letter codes: 20 standard residues + ambiguity B (Asx), Z (Glx),
        // J (Xle), X (Xaa) + rare U (Sec), O (Pyl) + stop '*'.
        private static readonly HashSet<char> ProteinCodes = new()
        {
            'A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I',
            'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V',
            'B', 'Z', 'J', 'X', 'U', 'O',
            '*'
        };

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
                else
                {
                    foreach (char c in line)
                    {
                        if (!char.IsWhiteSpace(c))
                            sequenceBuilder.Append(c);
                    }
                }
            }

            if (header != null && sequenceBuilder.Length > 0)
            {
                yield return CreateEntry(header, sequenceBuilder.ToString());
            }
        }

        // --- Opt-in alphabet-aware overloads ---------------------------------------------------

        /// <summary>
        /// Parses a FASTA string under a selectable <paramref name="alphabet"/>, returning
        /// <see cref="FastaRecord"/> instances whose raw (uppercased) sequence string is validated
        /// against the chosen alphabet. Out-of-alphabet characters throw <see cref="ArgumentException"/>.
        /// </summary>
        /// <remarks>
        /// <see cref="SequenceAlphabet.StrictDna"/> reproduces the default DNA-only behaviour
        /// (A/C/G/T only). The other alphabets are opt-in.
        /// </remarks>
        public static IEnumerable<FastaRecord> Parse(string fastaContent, SequenceAlphabet alphabet)
        {
            if (string.IsNullOrWhiteSpace(fastaContent))
                yield break;

            using var reader = new StringReader(fastaContent);
            foreach (var record in ParseReaderTyped(reader, alphabet))
            {
                yield return record;
            }
        }

        /// <summary>
        /// Parses a FASTA file under a selectable <paramref name="alphabet"/>.
        /// </summary>
        public static IEnumerable<FastaRecord> ParseFile(string filePath, SequenceAlphabet alphabet)
        {
            using var reader = new StreamReader(filePath);
            foreach (var record in ParseReaderTyped(reader, alphabet))
            {
                yield return record;
            }
        }

        /// <summary>
        /// Reads sequences from a FASTA file asynchronously under a selectable <paramref name="alphabet"/>.
        /// </summary>
        public static async IAsyncEnumerable<FastaRecord> ParseFileAsync(string filePath, SequenceAlphabet alphabet)
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
                        yield return CreateRecord(header, sequenceBuilder.ToString(), alphabet);
                    }

                    header = line.Substring(1).Trim();
                    sequenceBuilder.Clear();
                }
                else
                {
                    foreach (char c in line)
                    {
                        if (!char.IsWhiteSpace(c))
                            sequenceBuilder.Append(c);
                    }
                }
            }

            if (header != null && sequenceBuilder.Length > 0)
            {
                yield return CreateRecord(header, sequenceBuilder.ToString(), alphabet);
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
                // Emit '\n' explicitly (not AppendLine/Environment.NewLine) so FASTA output is
                // byte-identical across platforms — sequence files must not carry OS-dependent CRLF.
                sb.Append('>').Append(entry.Header).Append('\n');

                string seq = entry.Sequence.Sequence;
                for (int i = 0; i < seq.Length; i += lineWidth)
                {
                    int len = Math.Min(lineWidth, seq.Length - i);
                    sb.Append(seq.AsSpan(i, len)).Append('\n');
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
                else
                {
                    foreach (char c in line)
                    {
                        if (!char.IsWhiteSpace(c))
                            sequenceBuilder.Append(c);
                    }
                }
            }

            if (header != null && sequenceBuilder.Length > 0)
            {
                yield return CreateEntry(header, sequenceBuilder.ToString());
            }
        }

        private static IEnumerable<FastaRecord> ParseReaderTyped(TextReader reader, SequenceAlphabet alphabet)
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
                        yield return CreateRecord(header, sequenceBuilder.ToString(), alphabet);
                    }

                    header = line.Substring(1).Trim();
                    sequenceBuilder.Clear();
                }
                else
                {
                    foreach (char c in line)
                    {
                        if (!char.IsWhiteSpace(c))
                            sequenceBuilder.Append(c);
                    }
                }
            }

            if (header != null && sequenceBuilder.Length > 0)
            {
                yield return CreateRecord(header, sequenceBuilder.ToString(), alphabet);
            }
        }

        private static (string Id, string? Description) SplitHeader(string header)
        {
            // Parse header: typically "ID description". Split on the first space or tab.
            var parts = header.Split(new[] { ' ', '\t' }, 2);
            string id = parts[0];
            string? description = parts.Length > 1 ? parts[1] : null;
            return (id, description);
        }

        private static FastaEntry CreateEntry(string header, string sequence)
        {
            var (id, description) = SplitHeader(header);
            return new FastaEntry(id, description, new DnaSequence(sequence));
        }

        private static FastaRecord CreateRecord(string header, string sequence, SequenceAlphabet alphabet)
        {
            var (id, description) = SplitHeader(header);
            // Lower-case letters are accepted and mapped to upper-case (NCBI/Wikipedia FASTA),
            // matching the DnaSequence/RnaSequence/ProteinSequence constructors.
            string normalized = sequence.ToUpperInvariant();
            ValidateAgainstAlphabet(normalized, alphabet);
            return new FastaRecord(id, description, normalized, alphabet);
        }

        private static void ValidateAgainstAlphabet(string sequence, SequenceAlphabet alphabet)
        {
            for (int i = 0; i < sequence.Length; i++)
            {
                char c = sequence[i];
                bool valid = alphabet switch
                {
                    SequenceAlphabet.StrictDna => c is 'A' or 'C' or 'G' or 'T',
                    SequenceAlphabet.IupacNucleotide => IupacNucleotideCodes.Contains(c),
                    SequenceAlphabet.Rna => RnaCodes.Contains(c),
                    SequenceAlphabet.Protein => ProteinCodes.Contains(c),
                    _ => throw new ArgumentOutOfRangeException(nameof(alphabet), alphabet, "Unknown sequence alphabet.")
                };

                if (!valid)
                {
                    throw new ArgumentException(
                        $"Invalid character '{c}' at position {i} for alphabet {alphabet}.",
                        nameof(sequence));
                }
            }
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

    /// <summary>
    /// Represents a single FASTA record parsed under a selectable <see cref="SequenceAlphabet"/>.
    /// The raw (uppercased) sequence string is preserved verbatim, so IUPAC ambiguity codes,
    /// RNA (U), gaps, and protein residues are kept rather than coerced into a DNA type.
    /// </summary>
    public sealed class FastaRecord
    {
        public FastaRecord(string id, string? description, string sequence, SequenceAlphabet alphabet)
        {
            Id = id;
            Description = description;
            Sequence = sequence;
            Alphabet = alphabet;
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
        /// The raw (uppercased) sequence string, validated against <see cref="Alphabet"/>.
        /// </summary>
        public string Sequence { get; }

        /// <summary>
        /// The alphabet the sequence was validated against.
        /// </summary>
        public SequenceAlphabet Alphabet { get; }

        /// <summary>
        /// Full header line (ID + Description).
        /// </summary>
        public string Header => Description != null ? $"{Id} {Description}" : Id;

        public override string ToString() => $"{Id} ({Sequence.Length} residues, {Alphabet})";
    }
}
