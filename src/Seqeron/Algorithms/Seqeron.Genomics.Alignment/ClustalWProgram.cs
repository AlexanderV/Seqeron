using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Alignment;

/// <summary>What a ClustalW command wrote: the output files by name, and what the program printed to the error stream while reading its command line.</summary>
/// <param name="Files">Output file name to its text (alignment file, guide-tree files).</param>
/// <param name="Messages">Warnings of the command line.</param>
public sealed record ClustalWProgramResult(IReadOnlyDictionary<string, string> Files, IReadOnlyList<string> Messages);

/// <summary>
/// ClustalW 2.1 as a function: a command line (<see cref="ClustalWCommandLine"/>) plus the text of the files it names gives the text of the files the program writes. Covers the alignment commands
/// (<c>-ALIGN</c>, <c>-NEWTREE</c>, <c>-USETREE</c>, <c>-PROFILE1/2 -PROFILE</c>, <c>-SEQUENCES</c>, <c>-CONVERT</c>) with every alignment, guide-tree, row-order, output-format (<see cref="ClustalWAlignmentWriter"/>) and
/// secondary-structure option; input files in FASTA, CLUSTAL or GDE format (the formats of the verified oracles).
/// </summary>
public static class ClustalWProgram
{
    /// <summary>Runs a ClustalW command line.</summary>
    /// <param name="arguments">The arguments, as on the command line (<c>-INFILE=x.fa -TYPE=DNA -ALIGN -OUTPUT=GCG</c>).</param>
    /// <param name="files">The files named on the command line: name to text.</param>
    /// <exception cref="NotSupportedException">The verb or the input format is not implemented (<c>-TREE</c>, <c>-BOOTSTRAP</c>, PIR / MSF / RSF / EMBL / PHYLIP / NEXUS input).</exception>
    public static ClustalWProgramResult Run(IReadOnlyList<string> arguments, IReadOnlyDictionary<string, string> files)
    {
        ArgumentNullException.ThrowIfNull(arguments);
        ArgumentNullException.ThrowIfNull(files);
        var cmd = ClustalWCommandLine.Parse(arguments, files);
        var output = new Dictionary<string, string>(StringComparer.Ordinal);
        string Read(string name) => ReadInput(files, name);

        switch (cmd.Verb)
        {
            case ClustalWVerb.Align:
            case ClustalWVerb.AlignUsingTree:
            case ClustalWVerb.GuideTree:
            {
                string inFile = cmd.InFile ?? throw new ArgumentException("No input file(s) specified", nameof(arguments));
                var seqs = ReadFile(Read(inFile), cmd.Type, out bool dna, out _);
                var run = cmd.Run;
                var result = dna ? SequenceAligner.RunClustalW(seqs, cmd.Dna, run) : SequenceAligner.RunClustalW(seqs, cmd.Protein, run);
                string baseName = BaseName(inFile);
                if (result.GuideTree != null)
                    output[cmd.NewTreeFile ?? baseName + "dnd"] = result.GuideTree;
                if (result.Alignment != null)
                    output[cmd.OutFile ?? baseName + ClustalWAlignmentWriter.DefaultExtension(cmd.Run.Output.Format)] = ClustalWAlignmentWriter.Write(result.Alignment, cmd.Run.Output, cmd.Run.Structure);
                break;
            }

            case ClustalWVerb.Profile:
            case ClustalWVerb.Sequences:
            {
                string p1Name = cmd.Profile1 ?? throw new ArgumentException("No input file(s) specified", nameof(arguments));
                string p2Name = cmd.Profile2 ?? throw new ArgumentException("Only 1 profile file (profile 2) specified.", nameof(arguments));
                var p1 = ReadFile(Read(p1Name), cmd.Type, out bool dna, out var s1);
                var p2 = ReadFile(Read(p2Name), dna ? ClustalWSequenceType.Dna : ClustalWSequenceType.Protein, out _, out var s2);
                var run = cmd.Run with { ProfileStructure1 = s1, ProfileStructure2 = s2 };
                ClustalWRunResult result;
                if (cmd.Verb == ClustalWVerb.Profile)
                    result = dna ? SequenceAligner.RunClustalWProfile(p1, p2, cmd.Dna, run) : SequenceAligner.RunClustalWProfile(p1, p2, cmd.Protein, run);
                else
                    result = dna ? SequenceAligner.RunClustalWSequences(p1, p2, cmd.Dna, run) : SequenceAligner.RunClustalWSequences(p1, p2, cmd.Protein, run);
                if (result.GuideTree1 != null)
                    output[BaseName(p1Name) + "dnd"] = result.GuideTree1;
                if (result.GuideTree2 != null)
                    output[BaseName(p2Name) + "dnd"] = result.GuideTree2;
                if (result.GuideTree != null)
                    output[cmd.NewTreeFile ?? BaseName(p2Name) + "dnd"] = result.GuideTree;
                if (result.Alignment != null)
                {
                    string baseName = BaseName(cmd.Verb == ClustalWVerb.Profile ? p1Name : p2Name);
                    output[cmd.OutFile ?? baseName + ClustalWAlignmentWriter.DefaultExtension(cmd.Run.Output.Format)] = ClustalWAlignmentWriter.Write(result.Alignment, cmd.Run.Output, cmd.Run.Structure);
                }

                break;
            }

            case ClustalWVerb.Convert:
            {
                // Clustal::outputNow: the sequences as they are (no alignment: rows of different lengths, weights 100, input order) in the chosen format
                string inFile = cmd.InFile ?? throw new ArgumentException("No input file(s) specified", nameof(arguments));
                var seqs = ReadFile(Read(inFile), cmd.Type, out bool dna, out _);
                int n = seqs.Count;
                var alignment = new ClustalWAlignment(seqs.Select(s => s.Name).ToList(), seqs.Select(s => s.Residues.ToUpperInvariant()).ToList(), Enumerable.Range(0, n).ToList(), !dna, Enumerable.Repeat(100, n).ToList());
                output[cmd.OutFile ?? BaseName(inFile) + ClustalWAlignmentWriter.DefaultExtension(cmd.Run.Output.Format)] = ClustalWAlignmentWriter.Write(alignment, cmd.Run.Output, cmd.Run.Structure);
                break;
            }

            default:
                throw new NotSupportedException($"ClustalW verb {cmd.Verb} is not implemented in this library (-TREE and -BOOTSTRAP are registered as audit item A2-14f).");
        }

        return new ClustalWProgramResult(output, cmd.Messages);
    }

    private static string ReadInput(IReadOnlyDictionary<string, string> files, string name)
    {
        if (!files.TryGetValue(name, out var text))
            throw new ArgumentException($"Cannot open sequence file [{name}]: its text is not among the files given.", nameof(files));
        return text;
    }

    // Utility::getPath: the name without its extension, keeping the dot ("x.fa" -> "x."; no extension: the whole name and a dot)
    private static string BaseName(string name)
    {
        int dot = name.LastIndexOf('.');
        int slash = Math.Max(name.LastIndexOf('/'), name.LastIndexOf('\\'));
        return dot > slash && dot >= 0 ? name[..(dot + 1)] : name + ".";
    }

    // FileReader::checkInfile: the format by its first characters: CLUSTAL, GDE ('#' or '%' lines), else FASTA (Pearson)
    private static List<ClustalWSequence> ReadFile(string text, ClustalWSequenceType? type, out bool dna, out ClustalWProfileStructure? structure)
    {
        structure = null;
        string first = text.TrimStart().Split('\n', 2)[0];
        List<ClustalWSequence> seqs;
        if (first.StartsWith("CLUSTAL", StringComparison.Ordinal))
        {
            var f = ClustalWProfileFileReader.ParseClustal(text);
            seqs = [.. f.Sequences];
            structure = f.Structure;
        }
        else if (first.StartsWith('>') && !System.Text.RegularExpressions.Regex.IsMatch(first, @"^>(P1|F1|DL|DC|RL|RC|N3|N1);"))
        {
            seqs = ParseFasta(text);
        }
        else if (first.Length > 0 && first[0] is '#' or '%')
        {
            bool gdeDna = first[0] == '#';
            var f = ClustalWProfileFileReader.ParseGde(text, gdeDna);
            seqs = [.. f.Sequences];
            structure = f.Structure;
            type ??= gdeDna ? ClustalWSequenceType.Dna : ClustalWSequenceType.Protein;
        }
        else
        {
            throw new NotSupportedException("Only FASTA, CLUSTAL and GDE input files are implemented (PIR, MSF, RSF, EMBL, PHYLIP and NEXUS input are registered as audit item A2-14f).");
        }

        dna = (type ?? SequenceAligner.GuessClustalWSequenceType(seqs.Select(s => s.Residues))) == ClustalWSequenceType.Dna;
        return seqs;
    }

    // PearsonFileParser: '>' starts an entry, its name is the first word (<= 150 characters, blanks become '_'), the residues are the letters of the lines below
    private static List<ClustalWSequence> ParseFasta(string text)
    {
        var list = new List<ClustalWSequence>();
        string? name = null;
        var residues = new System.Text.StringBuilder();
        foreach (string raw in text.Split('\n'))
        {
            string line = raw.TrimEnd('\r');
            if (line.StartsWith('>'))
            {
                if (name != null)
                    list.Add(new ClustalWSequence(name, residues.ToString()));
                string header = line[1..].TrimStart(' ');
                int space = header.IndexOf(' ', StringComparison.Ordinal);
                name = (space >= 0 ? header[..space] : header).Replace(' ', '_');
                if (name.Length > 150)
                    name = name[..150];
                residues.Clear();
            }
            else if (name != null)
            {
                foreach (char c in line)
                {
                    char u = char.ToUpperInvariant(c);
                    if (ClustalWUserMatrix.ResidueCodes.Contains(u, StringComparison.Ordinal))
                        residues.Append(u);
                }
            }
        }

        if (name != null)
            list.Add(new ClustalWSequence(name, residues.ToString()));
        return list;
    }
}
