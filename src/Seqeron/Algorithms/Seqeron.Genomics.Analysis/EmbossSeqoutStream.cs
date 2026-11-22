#pragma warning disable S3776, S107, S1854
using System.Text;
using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Analysis;

/// <summary>The sequence output qualifiers of an ACD <c>seqout</c> / <c>seqoutall</c> (B09 WP61, F280).</summary>
/// <param name="OutSeq">The value of the parameter (<c>-outseq</c>): a file name, optionally with a <c>format::</c> prefix; <c>stdout</c> writes to the output of the program.</param>
/// <param name="Format"><c>-osformat</c> (default <c>fasta</c>).</param>
/// <param name="Single"><c>-ossingle</c>: every sequence is written to a file of its own.</param>
/// <param name="Directory"><c>-osdirectory</c>.</param>
/// <param name="Extension"><c>-osextension</c> (the extension of the files of <c>-ossingle</c>).</param>
/// <param name="DbName"><c>-osdbname</c>.</param>
/// <param name="Name"><c>-osname</c> (no effect on the names of the sequences, measured on 6.6.0).</param>
/// <param name="ExistingDirectories">The directories that exist where the program runs (null: all).</param>
public sealed record EmbossOutseqOptions(string OutSeq = "stdout", string? Format = null, bool Single = false, string? Directory = null, string? Extension = null, string? DbName = null, string? Name = null, IReadOnlyCollection<string>? ExistingDirectories = null)
{
    /// <summary>The date of the run (the headers of some writers); null = now.</summary>
    public DateTime? RunDate { get; init; }

    /// <summary>The command line of the run (the SAM writer prints it).</summary>
    public string CommandLine { get; init; } = "";
}

/// <summary>
/// The sequence output of a program (<c>AjPSeqout</c> of an ACD <c>seqoutall</c>: <c>ajSeqoutWriteSeq</c> + <c>ajSeqoutClose</c> of <c>ajseqwrite.c</c>; B09 WP61, F280): the format found by name, the type check of the
/// ACD type of the output (<c>ajSeqTypeCheckS</c>), the default names of unnamed sequences, the formats that write a set when the output is closed, the files of <c>-ossingle</c>, the footer of the DAS formats, the
/// BGZF blocks of BAM, the warning "No sequences written to output file". The code is that of <c>EmbossSeqretPrograms</c> for the part that does not need a feature file.
/// </summary>
internal sealed class EmbossSeqoutStream
{
    private readonly string _program;
    private readonly EmbossFeatureReadContext _context;
    private readonly EmbossOutseqOptions _o;
    private readonly string _outType;
    private readonly StringBuilder _seqOut = new();
    private readonly Dictionary<string, string> _singles = [];
    private readonly EmbossSeqWriteRun _run;
    private readonly EmbossSeqOutFormat? _format;
    private readonly string _formatName, _outName, _mainPath, _outDir;
    private string? _lastSingle;
    private int _unnamed, _count, _written;

    /// <summary>The output file could not be opened (the program dies).</summary>
    public bool OpenFailed { get; }

    /// <param name="program">The program name.</param>
    /// <param name="context">The messages of the run.</param>
    /// <param name="options">The qualifiers.</param>
    /// <param name="outType">The ACD <c>type</c> of the output (<c>stopprotein</c> ...); null = the type of the first sequence input.</param>
    /// <param name="defaultFormat">The format of the output when no qualifier gives one (<c>fasta</c>).</param>
    public EmbossSeqoutStream(string program, EmbossFeatureReadContext context, EmbossOutseqOptions options, string? outType, string defaultFormat = "fasta")
    {
        _program = program;
        _context = context;
        _o = options;
        _outType = outType ?? context.InputType ?? "";                // acdInTypeSeq: an output with no type of its own takes the type of the first input
        _formatName = options.Format ?? defaultFormat;
        _outName = options.OutSeq;
        // ajSeqoutUsaProcess: "format::file" names the output format (an empty format is "unknown")
        var um = System.Text.RegularExpressions.Regex.Match(_outName, "^([A-Za-z0-9-]*)::?(.*)$", System.Text.RegularExpressions.RegexOptions.CultureInvariant | System.Text.RegularExpressions.RegexOptions.Singleline);
        if (um.Success)
        {
            _formatName = um.Groups[1].Value.Length == 0 ? "unknown" : um.Groups[1].Value;
            _outName = um.Groups[2].Value;
        }

        _run = new EmbossSeqWriteRun { RunDate = options.RunDate ?? DateTime.Now, Program = program, CommandLine = options.CommandLine, Setoutdb = options.DbName ?? "", FileName = _outName };
        // ajSeqoutOpen (at the start of the run, by acdSetSeqout): -osdirectory is the directory of the output file, too (not only of the files of -ossingle)
        _outDir = EmbossOutPath.FixDirectory(options.Directory);
        _mainPath = options.Single ? _outName : EmbossOutPath.Combine(_outDir, _outName);
        if (!options.Single && !EmbossOutPath.DirectoryExists(options.ExistingDirectories, _mainPath))
        {
            EmbossOutPath.OpenFailed(context, _outName, _outDir, options.OutSeq);
            OpenFailed = true;
            return;
        }

        _format = EmbossSeqWrite.Find(_formatName);
    }

    /// <summary>The program dies: the message of <c>acdSetSeqout</c>.</summary>
    public EmbossSeqretOutput DiedOutput(string parameter = "outseq")
    {
        _context.Messages.Add($"Died: {_program} terminated: Bad value for '-{parameter}' with -auto defined");
        return new EmbossSeqretOutput("", "", 0, [.. _context.Messages]) { Died = true };
    }

    /// <summary><c>ajSeqoutWriteSeq</c>. A sequence that the format cannot write is reported and skipped.</summary>
    /// <exception cref="InvalidOperationException">The output is one on which the reference program dies by a fault.</exception>
    public void Write(EmbossSeq seq)
    {
        if (_format is null)
        {
            // ajSeqoutWriteSeq: the format is looked up when the first sequence is written (the message is printed once)
            if (_count++ == 0)
                _context.Error($"unknown output sequence format '{_formatName}'");
            return;
        }

        bool isNuc = seq.Type == "N", isProt = seq.Type == "P";
        if (!(isNuc && _format.Nucleotide) && !(isProt && _format.Protein))
        {
            string kind = isNuc ? "not supported for nucleotide sequences" : "not supported for protein sequences";
            if (!isNuc && !isProt)
                kind = "failed: unknown sequence type";
            _context.Error($"Sequence format '{_formatName}' {kind}");
            return;
        }

        bool outSingle = _o.Single;
        if (_format.Save)
        {
            // seqClone then seqWriteListAppend: the sequence (a copy cut to its range) is kept and the set is written when the output is closed
            var copy = EmbossSeqRange.CopyForList(seq);
            EmbossSeqRange.Trim(copy);
            if (outSingle && EmbossSeqWrite.WritesToNullFileDies(_format, Math.Max(copy.Seq.Length, _run.Saved.Count == 0 ? 0 : _run.Saved.Max(x => x.Seq.Length))))
                throw new InvalidOperationException("A set format with -ossingle writes the saved list to a file that is not open (the reference program dies by a fault).");
            if (_format.Writer is "Nexus" or "Msf" && copy.Seq.Length == 0)
                throw new InvalidOperationException("The nexus and msf writers die by a segmentation fault on a set that holds an empty sequence (measured on 6.6.0).");
            if (seq.Type.Length > 0)
                _run.OutType = seq.Type;
            EmbossSeqWrite.ListAppend(_run, new EmbossSeqSaved { Name = copy.Name, Desc = copy.Desc, Type = copy.Type, Seq = copy.SeqText, Weight = copy.Weight, Gene = seq.Genelist.Count > 0 ? seq.Genelist[0].Name : null }, seq.Entryname);
            _count++;
            _written++;
            return;
        }

        if (!EmbossSeqWrite.IsPorted(_format))
            throw new NotSupportedException($"The sequence output format '{_format.Name}' is not written by this port.");
        _run.UfoLocal = false;
        // ajSeqoutSetNameDefaultS: an unnamed sequence is called EMBOSS_001, EMBOSS_002 ... (EMBOSS with -ossingle: the file of every unnamed sequence has the same name)
        string name = seq.Name;
        if (name.Length == 0)
            name = outSingle ? "EMBOSS" : "EMBOSS_" + (++_unnamed).ToString("000", System.Globalization.CultureInfo.InvariantCulture);
        string residues = EmbossSeqRange.Cut(seq, out int ibegin, out _, out _);                  // seqClone: the range of the sequence
        residues = EmbossSeqType.CheckS(residues, _outType);                                       // seqClone: ajSeqTypeCheckS of the output type
        var target = outSingle ? new StringBuilder() : _seqOut;
        var item = EmbossSeqWrite.Clone(seq, residues, ibegin, _run);
        item.Name = name;
        _run.Count = outSingle ? 0 : _count;
        _run.HeaderCount = _count;
        if (_format.Writer == "Debug")
        {
            string ext = _o.Extension ?? _o.Format ?? "fasta";
            string dirName = _o.Directory ?? "";
            if (dirName.Length > 0 && !dirName.EndsWith('/'))
                dirName += "/";
            _run.Debug = new EmbossSeqDebugInfo
            {
                Usa = _o.OutSeq, InputFormat = seq.Formatstr, OutputFormat = _formatName, Filename = _outName, Directory = dirName, Entryname = seq.Entryname, Single = outSingle, Features = false,
                FileNameOpen = outSingle ? (name + "." + ext).ToLowerInvariant() : _outName, Extension = ext, OutputType = _outType.Length > 0 ? _outType : "gapany",
            };
        }

        EmbossSeqWrite.Write(target, _format, item, _run);
        if (outSingle)
        {
            string file = (name + "." + (_o.Extension ?? _o.Format ?? "fasta")).ToLowerInvariant();
            string path = EmbossOutPath.Combine(_outDir, file);
            if (!EmbossOutPath.DirectoryExists(_o.ExistingDirectories, path))
            {
                if (!EmbossOutPath.SingleOpenSurvives(_format))
                    throw new InvalidOperationException("The directory of a file of -ossingle does not exist: the library warns and the reference program dies by a segmentation fault.");
                EmbossOutPath.SingleOpenFailed(_context, file, _outDir);
            }
            else
            {
                _lastSingle = path;
                _singles[path] = target.ToString();
            }
        }

        _count++;
        _written++;
    }

    /// <summary><c>ajSeqoutClose</c>: the sets, the BAM blocks and the footer are written, the warning for an empty output is printed.</summary>
    public EmbossSeqretOutput Close()
    {
        bool outSingle = _o.Single;
        if (_format is { Save: true } && _run.Saved.Count > 0 && !outSingle)
        {
            _run.Count = _count;                                       // ajSeqoutClose writes the saved list
            EmbossSeqWrite.WriteSaved(_seqOut, _format, _run);
        }

        // seqCleanBam: the BGZF blocks of a BAM file are written when the output is closed (with -ossingle into the only file the writer can use)
        if (_format is { Writer: "Bam" } && EmbossSeqWrite.FinishBam(_run) is { } bam)
        {
            if (_lastSingle is not null)
                _singles[_lastSingle] += bam;
            else
                _seqOut.Append(bam);
        }

        // seqCleanDasdna / seqCleanDasseq: the footer of a DAS file is written when the file that is open at the end of the run is closed (with -ossingle only the last file has it)
        if (_written > 0 && _format is not null && EmbossSeqWrite.Footer(_format) is { Length: > 0 } footer)
        {
            if (_lastSingle is not null)
                _singles[_lastSingle] += footer;
            else if (!outSingle)
                _seqOut.Append(footer);
        }

        if (_context.Fatal)
            return new EmbossSeqretOutput(_seqOut.ToString(), "", _written, [.. _context.Messages]) { Died = true, SingleFiles = _singles, SequenceFileName = _mainPath };
        if (_count == 0 && !outSingle)                                 // with -ossingle there is no output file that stays empty: the files are opened for the sequences
            _context.Warn($"No sequences written to output file '{_outName}'");
        return new EmbossSeqretOutput(_seqOut.ToString(), "", _written, [.. _context.Messages]) { SingleFiles = _singles, SequenceFileName = _mainPath };
    }
}
