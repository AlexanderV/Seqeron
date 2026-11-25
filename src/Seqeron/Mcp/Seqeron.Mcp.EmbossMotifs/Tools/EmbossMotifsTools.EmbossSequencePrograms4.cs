// SPDX-License-Identifier: GPL-2.0-or-later
// This file is part of Seqeron.Mcp.EmbossMotifs, licensed under the terms stated in LICENSE (this directory).
using System.ComponentModel;
using System.Text;
using ModelContextProtocol.Server;
using Seqeron.Mcp.Analysis.Tools;
using static Seqeron.Mcp.Analysis.Tools.AnalysisTools;

namespace Seqeron.Mcp.EmbossMotifs.Tools;

public partial class EmbossMotifsTools
{
    [McpServerTool(Name = "patmatmotifs", Title = "Protein — EMBOSS patmatmotifs (scan a protein with the PROSITE database)", ReadOnly = true)]
    [Description("EMBOSS 6.6.0 patmatmotifs: the protein sequence is searched with every pattern of the PROSITE database in the form prosextract makes (prositeLines: the text of PROSITE/prosite.lines: identifier, accession, description word, pattern and regular expression of every PATTERN entry; the tool prosextract makes it from prosite.dat), the hits are the report out.txt (default format dbmotif: every report format of EMBOSS, the report qualifiers of reportQualifiers); the simple patterns of post-translational modification sites (myristyl, asn_glycosylation, camp_phospho_site, pkc_phospho_site, ck2_phospho_site, tyr_phospho_site) are skipped unless noprune; with full the documentation of every pattern that has a hit is appended (documentation: accession -> the text of the file PROSITE/PSnnnnn that prosextract makes). Equals the installed program on 300 runs (every input format of seqret, ranges, the report options, the messages of the library; standard output, error stream, exit status and written files byte for byte). B09 WP67, F306.")]
    public static EmbossProgramToolResult Patmatmotifs(
        [Description("The text of the protein sequence file (any input format of seqret), or a sequence as letters.")] string sequence,
        [Description("The text of the PROSITE data file prosite.lines made by prosextract (PROSITE/prosite.lines): required, the sandbox of the tool has no PROSITE database.")] string? prositeLines = null,
        [Description("The documentation files of the PROSITE entries made by prosextract: accession (PS00001) -> text; read with full.")] IReadOnlyDictionary<string, string>? documentation = null,
        [Description("Provide the full documentation of the patterns that match (-full).")] bool full = false,
        [Description("Ignore the simple patterns (-prune, default true; false = -noprune).")] bool prune = true,
        [Description("The report format (-rformat, default dbmotif): simple, embl, genbank, gff, gff3, pir, swiss, dasgff, debug, listfile, dbmotif, diffseq, draw, restrict, excel, feattable, motif, nametable, regions, seqtable, srs, summary, table, tagseq.")] string? reportFormat = null,
        [Description("Report qualifiers as typed on the command line: -raccshow -rdesshow -rusashow -norscoreshow -norstrandshow -rmaxall N -rmaxseq N (-rformat is the parameter reportFormat).")] string? reportQualifiers = null,
        [Description("The input format (-sformat1); default: the program finds the format.")] string? sequenceFormat = null,
        [Description("Sequence input qualifiers as typed after the sequence on the command line: -sbegin1 -send1 -sreverse1, -slower1 -supper1 -snucleotide1 -sprotein1 -scircular1, -sid1 name, -sdbname1 db, -squick1, -iquery1 query, -ioffset1 bytes.")] string? qualifiers = null,
        CancellationToken cancellationToken = default,
        [Description("More files by name (file name -> text) that the program may read (a file named by a USA, @file, the PROSITE data files of the data directory).")] IReadOnlyDictionary<string, string>? files = null)
    {
        if (string.IsNullOrWhiteSpace(sequence))
            throw new ArgumentException("The sequence cannot be empty", nameof(sequence));
        cancellationToken.ThrowIfCancellationRequested();
        var all = new Dictionary<string, string>(files ?? new Dictionary<string, string>());
        var args = new List<string>();
        all["sequence.dat"] = sequence;
        args.AddRange(["-sequence", SequenceSpec(null, sequence, "sequence.dat")]);
        if (!string.IsNullOrWhiteSpace(sequenceFormat))
            args.AddRange(["-sformat1", sequenceFormat]);
        AddQualifierTokens(args, qualifiers);
        if (!string.IsNullOrEmpty(prositeLines))
            all["PROSITE/prosite.lines"] = prositeLines;
        if (documentation is not null)
            foreach (var kv in documentation)
                all["PROSITE/" + kv.Key] = kv.Value;
        if (full)
            args.Add("-full");
        if (!prune)
            args.Add("-noprune");
        if (!string.IsNullOrEmpty(reportFormat))
            args.AddRange(["-rformat", reportFormat]);
        AddQualifierTokens(args, reportQualifiers);
        args.AddRange(["-outfile", "out.txt"]);
        return RunEmbossMultiProgram("patmatmotifs", args, all, nameof(sequence), "out.txt");
    }
    [McpServerTool(Name = "patmatdb", Title = "Protein — EMBOSS patmatdb (search protein sequences with a PROSITE pattern)", ReadOnly = true)]
    [Description("EMBOSS 6.6.0 patmatdb: every protein sequence is searched with the PROSITE pattern motif (case-independent; < > anchor, x any residue, [ ] a class, { } excluded residues, (n) (n,m) counts) turned into a regular expression of the library (PCRE 7); the hits are the report out.txt (default format dbmotif: every report format of EMBOSS, the report qualifiers of reportQualifiers); a text that the PROSITE syntax refuses (an unmatched bracket, a character that is no residue) stops the program with the message of the library, a pattern that the regular expression library refuses (a range out of order, a count above 65535, nothing to repeat) makes the reference program die by a segmentation fault and is refused with the message of the library. Equals the installed program on 300 runs (every input format of seqret, ranges, the report options, the messages of the library; standard output, error stream, exit status and written files byte for byte). B09 WP67, F306.")]
    public static EmbossProgramToolResult Patmatdb(
        [Description("The text of the protein sequence file (any input format of seqret), or a sequence as letters.")] string sequence,
        [Description("The PROSITE pattern to search for (-motif, at least 2 characters), e.g. [DE](2)HS{P}X(2)PX(2,4)C.")] string motif,
        [Description("The report format (-rformat, default dbmotif): simple, embl, genbank, gff, gff3, pir, swiss, dasgff, debug, listfile, dbmotif, diffseq, draw, restrict, excel, feattable, motif, nametable, regions, seqtable, srs, summary, table, tagseq.")] string? reportFormat = null,
        [Description("Report qualifiers as typed on the command line: -raccshow -rdesshow -rusashow -norscoreshow -norstrandshow -rmaxall N -rmaxseq N (-rformat is the parameter reportFormat).")] string? reportQualifiers = null,
        [Description("The input format (-sformat1); default: the program finds the format.")] string? sequenceFormat = null,
        [Description("Sequence input qualifiers as typed after the sequence on the command line: -sbegin1 -send1 -sreverse1, -slower1 -supper1 -snucleotide1 -sprotein1 -scircular1, -sid1 name, -sdbname1 db, -squick1, -iquery1 query, -ioffset1 bytes.")] string? qualifiers = null,
        CancellationToken cancellationToken = default,
        [Description("More files by name (file name -> text) that the program may read (a file named by a USA, @file, the PROSITE data files of the data directory).")] IReadOnlyDictionary<string, string>? files = null)
    {
        if (string.IsNullOrWhiteSpace(sequence))
            throw new ArgumentException("The sequence cannot be empty", nameof(sequence));
        if (string.IsNullOrWhiteSpace(motif))
            throw new ArgumentException("The motif cannot be empty", nameof(motif));
        cancellationToken.ThrowIfCancellationRequested();
        var all = new Dictionary<string, string>(files ?? new Dictionary<string, string>());
        var args = new List<string>();
        all["sequence.dat"] = sequence;
        args.AddRange(["-sequence", SequenceSpec(null, sequence, "sequence.dat")]);
        args.AddRange(["-motif", motif]);
        if (!string.IsNullOrWhiteSpace(sequenceFormat))
            args.AddRange(["-sformat1", sequenceFormat]);
        AddQualifierTokens(args, qualifiers);
        if (!string.IsNullOrEmpty(reportFormat))
            args.AddRange(["-rformat", reportFormat]);
        AddQualifierTokens(args, reportQualifiers);
        args.AddRange(["-outfile", "out.txt"]);
        return RunEmbossMultiProgram("patmatdb", args, all, nameof(sequence), "out.txt");
    }
    [McpServerTool(Name = "prosextract", Title = "Protein — EMBOSS prosextract (build the PROSITE database for patmatmotifs)", ReadOnly = true)]
    [Description("EMBOSS 6.6.0 prosextract: for every PATTERN entry of prosite.dat (the text of the SIB file) the file PROSITE/prosite.lines (returned as the text) holds the identifier, the accession, the first word of the description, the pattern and its regular expression; for every {PSnnnnn; NAME} entry of prosite.doc one documentation file PROSITE/PSnnnnn (returned in files, by the path of the data directory /usr/share/EMBOSS/data/PROSITE/PSnnnnn) holds the text between {BEGIN} and {END}. Equals the installed program on 300 runs (every input format of seqret, ranges, the report options, the messages of the library; standard output, error stream, exit status and written files byte for byte). B09 WP67, F306.")]
    public static EmbossProgramToolResult Prosextract(
        [Description("The text of the PROSITE data file prosite.dat (ftp.expasy.org/databases/prosite).")] string prositeDat,
        [Description("The text of the PROSITE documentation file prosite.doc.")] string prositeDoc,
        CancellationToken cancellationToken = default,
        [Description("More files by name (file name -> text) that the program may read (a file named by a USA, @file, the PROSITE data files of the data directory).")] IReadOnlyDictionary<string, string>? files = null)
    {
        if (string.IsNullOrWhiteSpace(prositeDat))
            throw new ArgumentException("The PROSITE data file cannot be empty", nameof(prositeDat));
        cancellationToken.ThrowIfCancellationRequested();
        var all = new Dictionary<string, string>(files ?? new Dictionary<string, string>()) { ["prosite.dat"] = prositeDat, ["prosite.doc"] = prositeDoc ?? "" };
        var args = new List<string> { "-prositedir", "./" };
        return RunEmbossMultiProgram("prosextract", args, all, nameof(prositeDat), "/usr/share/EMBOSS/data/PROSITE/prosite.lines");
    }
}
