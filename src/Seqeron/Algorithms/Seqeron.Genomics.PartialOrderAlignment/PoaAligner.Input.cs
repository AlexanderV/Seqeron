// SPDX-License-Identifier: GPL-2.0-or-later
// Derived from POA 2.0 (Partial Order Alignment, https://sourceforge.net/projects/poamsa/), Copyright (C) 2004-2006 Christopher Lee; GPL-2.0-or-later.
// Part of Seqeron.Genomics.PartialOrderAlignment, licensed under the GNU General Public License, version 2 or (at your option) any later version (see LICENSE in the project directory).
// This file is NOT covered by the MIT licence of the other Seqeron assemblies.
using System.Globalization;
using System.Text;

namespace Seqeron.Genomics.PartialOrderAlignment;

using Seqeron.Genomics.Alignment;
using static Seqeron.Genomics.Alignment.SequenceAligner;

/// <summary>One alignment / partial-order file of the <c>poa</c> command line (<c>-read_msa</c>, <c>-read_msa2</c>, an entry of <c>-read_msa_list</c>): its name, text and the optional <c>-subset</c> / <c>-remove</c> list.</summary>
/// <param name="Name">The file name (for the program's messages).</param>
/// <param name="Text">null = the file cannot be opened. The text of the file: a PO file (<c>VERSION=LPO.0.1</c>), a FASTA-PIR alignment (<c>&gt;name title</c>) or a CLUSTAL alignment.</param>
/// <param name="SubsetText">The text of the <c>-subset</c> / <c>-remove</c> file (lines <c>SOURCENAME=name</c>); null = none.</param>
/// <param name="RemoveListed"><c>true</c> for <c>-remove</c> (keep every sequence except the listed ones), <c>false</c> for <c>-subset</c>.</param>
internal sealed record PoaMsaFile(string Name, string? Text, string? SubsetName = null, string? SubsetText = null, bool RemoveListed = false);

/// <summary>The inputs of a <c>poa</c> run in the program's order: the first MSA, the second MSA, the MSA list entries, the FASTA sequences.</summary>
internal sealed record PoaProgramInputs(PoaMsaFile? Msa, PoaMsaFile? Msa2, IReadOnlyList<PoaMsaFile> List, IReadOnlyList<PoaSequence> Fasta);

/// <summary>The program stops with a message (<c>exit(1)</c>) while reading an input.</summary>
public sealed class PoaInputException(string message, int exitCode = 1) : Exception(message)
{
    /// <summary>The exit status of the program (1, or 255 = <c>exit(-1)</c> of an invalid memory request).</summary>
    public int ExitCode { get; } = exitCode;
}

public static partial class PoaAligner
{
    private sealed partial class PoaEngine
    {
        /// <summary>The C stdio calls of <c>msa_format.c</c> / <c>lpo_format.c</c> over a string (<c>fgets</c>, <c>getc</c> and the <c>fscanf</c> conversions they use).</summary>
        private sealed class CText(string text)
        {
            private int _p;

            public int Getc() => _p < text.Length ? text[_p++] : -1;

            /// <summary><c>fgets(buf, size, f)</c>: at most <c>size - 1</c> characters up to and including the line feed; null at the end of the text.</summary>
            public string? Fgets(int size)
            {
                if (_p >= text.Length) return null;
                int start = _p;
                while (_p < text.Length && _p - start < size - 1)
                {
                    if (text[_p++] == '\n') break;
                }

                return text[start.._p];
            }

            private static bool IsSpace(char c) => c == ' ' || (c >= '\t' && c <= '\r');

            public void SkipWhite()
            {
                while (_p < text.Length && IsSpace(text[_p])) _p++;
            }

            /// <summary>A literal of the format: false at the first character that differs (it stays unread).</summary>
            public bool Literal(string lit)
            {
                foreach (char c in lit)
                {
                    if (_p >= text.Length || text[_p] != c) return false;
                    _p++;
                }

                return true;
            }

            /// <summary><c>%s</c>: skips white space, reads the characters up to the next white space.</summary>
            public string? Word()
            {
                SkipWhite();
                int start = _p;
                while (_p < text.Length && !IsSpace(text[_p])) _p++;
                return _p > start ? text[start.._p] : null;
            }

            /// <summary><c>%[^\n]</c> (or <c>%[^\n\r]</c>): at least one character up to the line end.</summary>
            public string? Rest(bool stopAtCr = false)
            {
                int start = _p;
                while (_p < text.Length && text[_p] != '\n' && !(stopAtCr && text[_p] == '\r')) _p++;
                return _p > start ? text[start.._p] : null;
            }

            /// <summary><c>%d</c>.</summary>
            public bool Int(out int value)
            {
                value = 0;
                SkipWhite();
                int q = _p;
                bool neg = false;
                if (q < text.Length && (text[q] == '+' || text[q] == '-'))
                {
                    neg = text[q] == '-';
                    q++;
                }

                int ds = q;
                long v = 0;
                while (q < text.Length && text[q] >= '0' && text[q] <= '9')
                {
                    v = Math.Min(v * 10 + (text[q] - '0'), int.MaxValue);
                    q++;
                }

                if (q == ds) return false;
                _p = q;
                value = (int)(neg ? -v : v);
                return true;
            }

            /// <summary><c>%*[\r\n]</c>.</summary>
            public void SkipLineEnds()
            {
                while (_p < text.Length && (text[_p] == '\r' || text[_p] == '\n')) _p++;
            }

            /// <summary><c>%*[ \t]</c>.</summary>
            public void SkipBlanks()
            {
                while (_p < text.Length && (text[_p] == ' ' || text[_p] == '\t')) _p++;
            }
        }

        // ------------------------------------------------------------------ lpo_index_symbols

        /// <summary><c>lpo_index_symbols</c>: letters read from a PO file are ASCII characters; the first letter decides whether they still need the translation to symbol indices.</summary>
        private void IndexSymbols(PoaLpo lpo)
        {
            if (lpo.Letters.Length == 0 || lpo.Letters[0].Letter < _m.Symbols.Length) return;
            foreach (var l in lpo.Letters)
            {
                int k = _m.Symbols.Length - 1;
                for (int j = 0; j < _m.Symbols.Length; j++)
                {
                    if (_m.Symbols[j] == (char)l.Letter)
                    {
                        k = j;
                        break;
                    }
                }

                l.Letter = k;
            }
        }

        // ------------------------------------------------------------------ the program's input loop

        /// <summary>Reads every input of the command line in the order of <c>main</c> and runs the alignment (<c>buildup_progressive_lpo</c>) and the writers.</summary>
        /// <param name="inputs">The input files and sequences.</param>
        /// <param name="err">Receives the program's progress lines of standard error.</param>
        /// <exception cref="PoaInputException">The program stops with this message.</exception>
        internal PoaRunResult RunProgram(PoaProgramInputs inputs, StringBuilder err, string? fastaFile)
        {
            Header();
            var all = new List<PoaLpo>();
            int inputIndex = 0;
            void Take(PoaLpo lpo)
            {
                foreach (var s in lpo.Sources) s.Input = inputIndex++;
                all.Add(lpo);
            }

            if (inputs.Msa != null)
            {
                var lpo = ReadPartialOrderFile(inputs.Msa);
                err.Append("...Read ").Append(lpo.Sources.Count.ToString(CultureInfo.InvariantCulture)).Append(" sequences from MSA file ").Append(inputs.Msa.Name).Append("...\n");
                Take(lpo);
            }

            if (inputs.Msa2 != null)
            {
                var lpo = ReadPartialOrderFile(inputs.Msa2);
                err.Append("...Read ").Append(lpo.Sources.Count.ToString(CultureInfo.InvariantCulture)).Append(" sequences from second MSA file ").Append(inputs.Msa2.Name).Append("...\n");
                Take(lpo);
            }

            foreach (var entry in inputs.List)
            {
                var lpo = ReadPartialOrderFile(entry);
                err.Append("...Read ").Append(lpo.Sources.Count.ToString(CultureInfo.InvariantCulture)).Append(" sequences from PO list entry ").Append(entry.Name).Append("...\n");
                Take(lpo);
            }

            if (fastaFile != null)
                err.Append("...Read ").Append(inputs.Fasta.Count.ToString(CultureInfo.InvariantCulture)).Append(" sequences from sequence file ").Append(fastaFile).Append("...\n");
            foreach (var rec in inputs.Fasta)
            {
                var lpo = MakeSequence(rec);
                if (lpo.Sequence.Length == 0) continue;
                Initialize(lpo, inputIndex++);
                all.Add(lpo);
            }

            if (all.Count == 0) throw new PoaInputException("No input sequences were provided; use one of the -read_ flags.\nExiting.");
            return Finish(all);
        }

        /// <summary><c>read_partial_order_file</c> of <c>main.c</c>.</summary>
        private PoaLpo ReadPartialOrderFile(PoaMsaFile f)
        {
            if (f.Text == null) throw new PoaInputException($"Couldn't open MSA file {f.Name}.\nExiting.");
            if (f.SubsetName != null && f.SubsetText == null) throw new PoaInputException($"Couldn't open subset file {f.SubsetName}.\nExiting.");
            var lpo = ReadMsaSelect(f.Text, f.SubsetText, false, f.RemoveListed);
            if (lpo == null || lpo.Sources.Count == 0)
            {
                throw new PoaInputException(f.SubsetText != null
                    ? $"MSA file {f.Name}, filtered with subset file, couldn't be read or contains no sequences.\nExiting."
                    : $"MSA file {f.Name} couldn't be read or contains no sequences.\nExiting.");
            }

            return lpo;
        }

        // ------------------------------------------------------------------ msa_format.c

        /// <summary><c>read_msa_select</c> with the format determined from the first line(s).</summary>
        internal PoaLpo? ReadMsaSelect(string text, string? subset, bool keepAllLinks, bool removeListed)
        {
            var f = new CText(text);
            string? line;
            while ((line = f.Fgets(511)) != null)
            {
                if (line[0] == '>') return ReadPir(f, line, subset, removeListed);
                if (line.StartsWith("VERSION=", StringComparison.Ordinal))
                    return subset != null ? ReadLpoSelect(f, subset, keepAllLinks, removeListed) : ReadLpo(f);
                if (line.StartsWith("CLUSTAL", StringComparison.Ordinal)) return ReadClustal(f, line, subset, removeListed);
                if (line[0] is '#' or '*' or ' ' or '\t' or '\n' or '\r') continue;
                return ReadClustal(f, line, subset, removeListed); // "Unable to determine MSA file type... trying CLUSTAL."
            }

            return null;
        }

        private static bool IsResidueChar(char ch) => (ch >= 'a' && ch <= 'z') || (ch >= 'A' && ch <= 'Z') || ch == '?' || ch == '[' || ch == ']';

        private static bool IsGapChar(char ch) => ch == '.' || ch == '-';

        private static bool IsNameFirstChar(char ch) => !(ch == '\n' || ch == '\r' || ch == ' ' || ch == '\t' || ch == '#' || ch == '*');

        private sealed class AlnRows
        {
            public readonly List<string> Names = new(), Titles = new();
            public readonly List<StringBuilder> Rows = new();
        }

        private static void AppendRow(StringBuilder row, string aln)
        {
            foreach (char ch in aln)
                if (IsResidueChar(ch) || IsGapChar(ch)) row.Append(ch);
        }

        /// <summary><c>read_clustal</c>.</summary>
        private PoaLpo? ReadClustal(CText fp, string firstLine, string? subset, bool removeListed)
        {
            var rows = new AlnRows();
            int currSeq = 0, lineNum = 0;
            bool expectRepeats = false, expectHeader = true;
            string? line = null;
            string aln = "";
            while ((firstLine != null && lineNum == 0) || (line = fp.Fgets(511)) != null)
            {
                if (firstLine != null && lineNum == 0) line = firstLine;
                lineNum++;
                if (expectHeader)
                {
                    if (line!.StartsWith("CLUSTAL", StringComparison.Ordinal))
                    {
                        expectHeader = false;
                        continue;
                    }

                    if (!IsNameFirstChar(line[0])) continue;
                    expectHeader = false;
                }

                if (!IsNameFirstChar(line![0]))
                {
                    currSeq = 0;
                    if (rows.Names.Count > 0) expectRepeats = true;
                    continue;
                }

                // sscanf(line, "%s %[^\n\r]", name, aln) < 2
                var ls = new CText(line);
                string? name = ls.Word();
                ls.SkipWhite();
                string? rest = ls.Rest(true);
                if (name == null || rest == null)
                    throw new PoaInputException($"Error: Trouble reading CLUSTAL-formatted file near line {lineNum}: \n>>>\n{line}<<<\nBailing out.\n(the program then dereferences a null alignment and crashes)");
                aln = rest;
                if (!expectRepeats)
                {
                    rows.Names.Add(name);
                    rows.Titles.Add("");
                    rows.Rows.Add(new StringBuilder());
                }
                else if (currSeq >= rows.Names.Count || name != rows.Names[currSeq])
                {
                    throw new PoaInputException($"Error: Trouble reading CLUSTAL-formatted file at line {lineNum}: \n>>>\n{line}<<<\nSequence name ({name}) does not match expected sequence name.  Bailing out.\n(the program then dereferences a null alignment and crashes)");
                }

                AppendRow(rows.Rows[currSeq], aln);
                currSeq++;
            }

            return FromRows(rows, subset, removeListed);
        }

        /// <summary><c>read_pir</c>. A line that holds no text (<c>sscanf("%[^\n\r]")</c> does not match) re-uses the previous row: <c>aln</c> keeps its old value.</summary>
        private PoaLpo? ReadPir(CText fp, string firstLine, string? subset, bool removeListed)
        {
            var rows = new AlnRows();
            int currSeq = 0, lineNum = 0;
            string? line = null;
            string aln = "";
            while ((firstLine != null && lineNum == 0) || (line = fp.Fgets(511)) != null)
            {
                if (firstLine != null && lineNum == 0) line = firstLine;
                lineNum++;
                if (line![0] == '>')
                {
                    var ls = new CText(line);
                    ls.Literal(">");
                    string? name = ls.Word();
                    string title = "";
                    if (name != null)
                    {
                        ls.SkipWhite();
                        title = ls.Rest(true) ?? "";
                    }

                    if (name == null)
                        throw new PoaInputException($"Error: Trouble reading PIR-formatted file near line {lineNum} (no sequence name?):\n>>>\n{line}<<<\nBailing out.\n" + (rows.Names.Count > 0 ? "(the program then dereferences a null alignment and crashes)" : string.Empty));
                    rows.Names.Add(name);
                    rows.Titles.Add(title);
                    rows.Rows.Add(new StringBuilder());
                    currSeq = rows.Names.Count - 1;
                }
                else if (line[0] == '#' || line[0] == '*')
                {
                    continue;
                }
                else
                {
                    var ls = new CText(line);
                    string? rest = ls.Rest(true);
                    if (rest != null) aln = rest;
                    if (rows.Names.Count == 0)
                        return null; // "no preceding '>seqname' line?": the program frees nothing and returns NULL
                    AppendRow(rows.Rows[currSeq], aln);
                }
            }

            return FromRows(rows, subset, removeListed);
        }

        private PoaLpo? FromRows(AlnRows rows, string? subset, bool removeListed)
        {
            FilterSequenceSet(rows, subset, removeListed);
            if (rows.Names.Count == 0) return null;
            var lpo = LpoFromAlnMat(rows);
            lpo.Name = lpo.Sources[0].Name;
            lpo.Title = lpo.Sources[0].Title;
            return lpo;
        }

        /// <summary><c>filter_sequence_set</c>: <c>SOURCENAME=name</c> lines of the select file.</summary>
        private static void FilterSequenceSet(AlnRows rows, string? subset, bool removeListed)
        {
            int n = rows.Names.Count;
            var keep = new bool[n];
            for (int i = 0; i < n; i++) keep[i] = subset == null || removeListed;
            if (subset != null)
            {
                var sf = new CText(subset);
                while (sf.Literal("SOURCENAME="))
                {
                    string? name = sf.Rest(true);
                    if (name == null) break;
                    sf.SkipLineEnds(); // %*[\r\n]
                    for (int i = 0; i < n; i++)
                        if (rows.Names[i] == name) keep[i] = !removeListed;
                }
            }

            for (int i = n; i-- > 0;)
            {
                if (keep[i]) continue;
                rows.Names.RemoveAt(i);
                rows.Titles.RemoveAt(i);
                rows.Rows.RemoveAt(i);
            }
        }

        /// <summary><c>lpo_from_aln_mat</c>: the alignment rows become a partial order by fusing every row into a consensus row (the first residue of each column) and stripping the consensus again.</summary>
        private PoaLpo LpoFromAlnMat(AlnRows rows)
        {
            int nSeqs = rows.Names.Count;
            int maxAln = rows.Rows.Max(r => r.Length);
            var cons = new char[maxAln];
            for (int col = 0; col < maxAln; col++)
            {
                cons[col] = '-';
                for (int i = 0; i < nSeqs; i++)
                {
                    if (col < rows.Rows[i].Length && IsResidueChar(rows.Rows[i][col]))
                    {
                        cons[col] = rows.Rows[i][col];
                        break;
                    }
                }
            }

            char Case(char ch)
            {
                if (_o.Case == PoaCase.Lower) return char.ToLowerInvariant(ch);
                return _o.Case == PoaCase.Upper ? char.ToUpperInvariant(ch) : ch;
            }

            var columnIds0 = new List<int>();
            var resIds0 = new int[maxAln];
            var sb = new StringBuilder();
            for (int col = 0; col < maxAln; col++)
            {
                char ch = Case(cons[col]);
                if (IsResidueChar(ch))
                {
                    sb.Append(ch);
                    columnIds0.Add(col);
                    resIds0[col] = sb.Length - 1;
                }
            }

            var lposeq = new PoaLpo { Name = "consens_row", Title = string.Empty, Sequence = sb.ToString() };
            Initialize(lposeq, -1);
            for (int i = 0; i < nSeqs; i++)
            {
                string rawName = rows.Names[i];
                if (rawName.Length > 31) throw new PoaInputException($"The sequence name '{rawName}' is longer than 31 characters: the program copies it with strcpy into a 32-byte field and aborts (glibc buffer overflow detection).");
                var curr = new PoaLpo { Name = Trunc31(rawName), Title = rows.Titles[i] };
                var sbr = new StringBuilder();
                var colIds = new List<int>();
                var row = rows.Rows[i];
                for (int col = 0; col < row.Length; col++)
                {
                    char ch = Case(row[col]);
                    if (IsResidueChar(ch))
                    {
                        sbr.Append(ch);
                        colIds.Add(col);
                    }
                }

                curr.Sequence = sbr.ToString();
                if (curr.Sequence.Length == 0)
                    throw new PoaInputException($"...*** invalid memory request: al_x[0].\nThe sequence '{rawName}' of the alignment holds no residue (CALLOC of size 0: the program exits with status 255).", 255);
                Initialize(curr, -1);
                int len = curr.Length;
                var toCons = new int[len];
                var fromCons = new int[lposeq.Length];
                Array.Fill(fromCons, -1);
                var seqToPo = SeqToPo(lposeq, 0);
                for (int j = 0; j < len; j++)
                {
                    int resId = resIds0[colIds[j]];
                    int letterId = seqToPo[resId];
                    toCons[j] = letterId;
                    fromCons[letterId] = j;
                }

                FuseRingIdentities(lposeq.Letters, curr.Letters, fromCons, toCons);
                FuseLpo(lposeq, curr, fromCons, toCons);
            }

            StripSeq(lposeq, 0, false);
            lposeq.Sources.ForEach(s => s.Input = -1);
            return lposeq;
        }

        /// <summary><c>build_seq_to_po_index</c> for one source: the letter of every position of the sequence.</summary>
        private static int[] SeqToPo(PoaLpo lpo, int iseq)
        {
            var map = new int[lpo.Sources[iseq].Length];
            for (int i = 0; i < lpo.Letters.Length; i++)
                foreach (var (s, ipos) in lpo.Letters[i].Source)
                    if (s == iseq) map[ipos] = i;
            return map;
        }

        /// <summary><c>strip_seq_from_lpo</c>: removes a source and, unless links are kept, rebuilds the links from the remaining sequences.</summary>
        private static void StripSeq(PoaLpo lpo, int bad, bool keepAllLinks)
        {
            lpo.Sources.RemoveAt(bad);
            foreach (var l in lpo.Letters)
            {
                for (int k = l.Source.Count; k-- > 0;)
                {
                    if (l.Source[k].Iseq == bad) l.Source.RemoveAt(k);
                    else if (l.Source[k].Iseq >= bad) l.Source[k] = (l.Source[k].Iseq - 1, l.Source[k].Ipos);
                }
            }

            if (keepAllLinks) return;
            var toPo = Enumerable.Range(0, lpo.Sources.Count).Select(i => SeqToPo(lpo, i)).ToArray();
            foreach (var l in lpo.Letters)
            {
                l.Left.Clear();
                l.Right.Clear();
            }

            for (int i = 0; i < lpo.Sources.Count; i++)
            {
                int len = lpo.Sources[i].Length;
                for (int j = 0; j < len - 1; j++)
                {
                    int left = toPo[i][j], right = toPo[i][j + 1];
                    AddLink(lpo.Letters[left].Right, right);
                    AddLink(lpo.Letters[right].Left, left);
                }
            }
        }

        // ------------------------------------------------------------------ lpo_format.c: read_lpo / read_lpo_select

        private static void ReadLpoHeader(CText f, out string name, out string title, out int length, out int nsource)
        {
            name = "";
            title = "";
            if (f.Literal("VERSION=")) f.Word(); // fscanf stops at the first literal that does not match: the version word is not read then
            f.SkipWhite();
            if (f.Literal("NAME=")) name = f.Rest() ?? "";
            f.SkipWhite();
            if (f.Literal("TITLE=")) title = f.Rest() ?? "";
            f.SkipWhite();
            if (!f.Literal("LENGTH=") || !f.Int(out length))
                throw new PoaInputException("LENGTH / SOURCECOUNT of the PO file not found.");
            f.SkipWhite();
            nsource = 0;
            if (!f.Literal("SOURCECOUNT=") || !f.Int(out nsource))
                throw new PoaInputException("LENGTH / SOURCECOUNT of the PO file not found.");
        }

        private static PoaSourceInfo ReadSourceInfo(CText f)
        {
            f.SkipWhite();
            if (!f.Literal("SOURCENAME=")) throw new PoaInputException("SOURCENAME of the PO file not found.");
            string? name = f.Rest();
            f.SkipWhite();
            if (name == null || !f.Literal("SOURCEINFO=") || !f.Int(out int length) || !f.Int(out int istart) || !f.Int(out int weight) || !f.Int(out int bundle))
                throw new PoaInputException("SOURCEINFO of the PO file not found.");
            f.SkipBlanks();
            string title = f.Rest() ?? "";
            return new PoaSourceInfo { Name = Trunc31(name), Title = title, Length = length, Istart = istart, Weight = weight, BundleId = bundle };
        }

        /// <summary><c>read_lpo</c> (the first line, <c>VERSION=...</c>, was consumed by <c>fgets</c>; a PO file that does not begin with it still reads as the program reads it).</summary>
        private PoaLpo? ReadLpo(CText f)
        {
            ReadLpoHeader(f, out string name, out string title, out int length, out int nsource);
            var lpo = new PoaLpo { Name = Trunc31(name), Title = title };
            if (nsource <= 0 || length < 0) return null;
            var letters = new PoaLetter[length];
            for (int i = 0; i < length; i++) letters[i] = new PoaLetter { AlignRing = i, RingId = -1 };
            lpo.Letters = letters;
            for (int i = 0; i < nsource; i++) lpo.Sources.Add(ReadSourceInfo(f));
            var posCount = new int[nsource];
            for (int i = 0; i < length; i++)
            {
                f.SkipWhite();
                int c = f.Getc();
                if (c < 0) throw new PoaInputException("The PO file ends inside the letter list.");
                f.Literal(":");
                letters[i].Letter = c;
                int fieldId;
                while ((fieldId = f.Getc()) >= 0 && fieldId != '\n')
                {
                    if (!f.Int(out int value)) throw new PoaInputException("A letter field of the PO file is not a number.");
                    switch (fieldId)
                    {
                        case 'L':
                            if (value < 0 || value >= length) throw new PoaInputException("A link of the PO file points outside the letters (the program reads outside its array).");
                            AddLink(letters[i].Left, value);
                            AddLink(letters[value].Right, i);
                            break;
                        case 'S':
                            if (value < 0 || value >= nsource) throw new PoaInputException("A source number of the PO file is outside the sources (the program writes outside its array).");
                            letters[i].Source.Add((value, posCount[value]++));
                            break;
                        case 'A':
                            letters[i].AlignRing = value;
                            break;
                    }
                }
            }

            for (int i = 0; i < length; i++)
            {
                if (letters[i].RingId >= 0) continue;
                int j = i;
                int guard = 0;
                do
                {
                    letters[j].RingId = i;
                    j = letters[j].AlignRing;
                    if (j < 0 || j >= length || ++guard > length + 1) throw new PoaInputException("An alignment ring of the PO file is not closed (the program loops or reads outside its array).");
                }
                while (j != i);
            }

            return lpo;
        }

        /// <summary><c>read_lpo_select</c>: keeps (<c>-subset</c>) or drops (<c>-remove</c>) the listed sequences, compacts the letters that remain and rebuilds links and rings.</summary>
        private PoaLpo? ReadLpoSelect(CText f, string subset, bool keepAllLinks, bool removeListed)
        {
            const int invalidLink = -99;
            ReadLpoHeader(f, out string name, out string title, out int length, out int nsourceFile);
            var lpo = new PoaLpo { Name = Trunc31(name), Title = title };
            if (nsourceFile <= 0 || length < 0) return null;
            var letters = new PoaLetter[length];
            for (int i = 0; i < length; i++) letters[i] = new PoaLetter { AlignRing = i, RingId = -1 };
            var sources = new List<PoaSourceInfo>();
            for (int i = 0; i < nsourceFile; i++) sources.Add(ReadSourceInfo(f));
            var posCount = new int[nsourceFile];
            var iseqCompact = new int[nsourceFile];
            var lastPos = new int[nsourceFile];
            var matchPos = new int[nsourceFile];
            int retentionMode = removeListed ? 0 : 1; // default_retention_mode = 0, default_no_retention_mode = 1
            for (int i = 0; i < nsourceFile; i++) iseqCompact[i] = -retentionMode;
            var sf = new CText(subset);
            while (sf.Literal("SOURCENAME="))
            {
                string? sname = sf.Rest();
                if (sname == null) break;
                sf.SkipWhite();
                for (int i = 0; i < nsourceFile; i++)
                {
                    if (sources[i].Name == sname)
                    {
                        iseqCompact[i] = retentionMode - 1;
                        break;
                    }
                }
            }

            int jj = 0;
            for (int i = 0; i < nsourceFile; i++)
            {
                lastPos[i] = -1;
                if (iseqCompact[i] >= 0)
                {
                    iseqCompact[i] = jj;
                    lpo.Sources.Add(sources[i]);
                    matchPos[jj] = invalidLink;
                    jj++;
                }
                else
                {
                    iseqCompact[i] = -1;
                }
            }

            int nsource = jj;
            var linkList = new List<int>();
            var posCompact = new int[length];
            var ringOld = new int[length];
            int npos = 0;
            for (int i = 0; i < length; i++)
            {
                f.SkipWhite();
                int c = f.Getc();
                if (c < 0) throw new PoaInputException("The PO file ends inside the letter list.");
                f.Literal(":");
                letters[npos].Letter = c;
                linkList.Clear();
                bool keepThis = false;
                int fieldId;
                while ((fieldId = f.Getc()) >= 0 && fieldId != '\n')
                {
                    if (!f.Int(out int value)) throw new PoaInputException("A letter field of the PO file is not a number.");
                    switch (fieldId)
                    {
                        case 'L':
                            if (value < 0 || value >= length) throw new PoaInputException("A link of the PO file points outside the letters (the program reads outside its array).");
                            if (posCompact[value] >= 0) linkList.Add(posCompact[value]);
                            break;
                        case 'S':
                            if (value < 0 || value >= nsourceFile) throw new PoaInputException("A source number of the PO file is outside the sources (the program writes outside its array).");
                            if (iseqCompact[value] >= 0)
                            {
                                keepThis = true;
                                value = iseqCompact[value];
                                if (lastPos[value] >= 0 && !linkList.Contains(lastPos[value])) linkList.Add(lastPos[value]);
                                lastPos[value] = npos;
                                letters[npos].Source.Add((value, posCount[value]++));
                                letters[npos].RingId = -1;
                                letters[npos].AlignRing = i;
                                ringOld[i] = i;
                            }

                            break;
                        case 'A':
                            ringOld[i] = value;
                            if (keepThis) letters[npos].AlignRing = i;
                            break;
                    }
                }

                if (keepThis)
                {
                    foreach (var (iseq, ipos) in letters[npos].Source) matchPos[iseq] = ipos - 1;
                    foreach (int lk in linkList)
                    {
                        foreach (var (iseq, ipos) in letters[lk].Source)
                        {
                            if (keepAllLinks || ipos == matchPos[iseq])
                            {
                                AddLink(letters[npos].Left, lk);
                                AddLink(letters[lk].Right, npos);
                                break;
                            }
                        }
                    }

                    foreach (var (iseq, _) in letters[npos].Source) matchPos[iseq] = invalidLink;
                    posCompact[i] = npos++;
                }
                else
                {
                    posCompact[i] = -1;
                }
            }

            var kept = new PoaLetter[npos];
            Array.Copy(letters, kept, npos);
            lpo.Letters = kept;
            for (int i = 0; i < npos; i++)
            {
                if (kept[i].RingId >= 0) continue;
                int j = kept[i].AlignRing;
                int k = 0;
                int guard = 0;
                do
                {
                    if (j < 0 || j >= length || ++guard > 4 * length + 4) throw new PoaInputException("An alignment ring of the PO file is not closed (the program loops or reads outside its array).");
                    if (posCompact[j] >= 0)
                    {
                        k = posCompact[j];
                        kept[k].RingId = i;
                    }

                    j = ringOld[j];
                    if (j < 0 || j >= length) throw new PoaInputException("An alignment ring of the PO file is not closed (the program loops or reads outside its array).");
                    if (posCompact[j] >= 0) kept[posCompact[j]].AlignRing = k;
                }
                while (posCompact[j] != i);
            }

            _ = nsource;
            return lpo;
        }
    }
}
