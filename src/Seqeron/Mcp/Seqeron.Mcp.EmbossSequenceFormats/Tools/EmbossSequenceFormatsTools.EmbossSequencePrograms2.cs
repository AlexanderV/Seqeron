// SPDX-License-Identifier: GPL-2.0-or-later
// This file is part of Seqeron.Mcp.EmbossSequenceFormats, licensed under the terms stated in LICENSE (this directory).
using System.ComponentModel;
using System.Text;
using ModelContextProtocol.Server;
using Seqeron.Mcp.Analysis.Tools;
using static Seqeron.Mcp.Analysis.Tools.AnalysisTools;

namespace Seqeron.Mcp.EmbossSequenceFormats.Tools;

public partial class EmbossSequenceFormatsTools
{
    [McpServerTool(Name = "shuffleseq", Title = "Genome — EMBOSS shuffleseq (shuffle sequences keeping the composition)", ReadOnly = true)]
    [Description("EMBOSS 6.6.0 shuffleseq: every sequence is written shuffle times with its residues in random order (a random number for every position, then the Shell sort of the positions by it: ajStrRandom of the EMBOSS library); the random numbers are ajRandomDouble (Wichmann-Hill with the 101-term lagged generator of ajmath.c), seeded from the microseconds of the clock. Equals the installed program on 300 runs (every input format of seqret, ranges, the output formats, the messages of the library; byte for byte). B09 WP62, F281.")]
    public static EmbossProgramToolResult Shuffleseq(
        [Description("The text of the sequence file (any input format of seqret), or a sequence as letters.")] string sequence,
        [Description("The number of shuffles of every sequence (-shuffle, default 1).")] int? shuffle = null,
        [Description("The microseconds of the clock when the random number generator is seeded (the seed is microseconds % 9999 + 1; the program takes them from the clock of the machine): the same value gives the same sequences. Default: the clock.")] int? microseconds = null,
        [Description("The input format (-sformat1); default: the program finds the format.")] string? sequenceFormat = null,
        [Description("Sequence input qualifiers as typed after the sequence on the command line: -sbegin1 -send1 -sreverse1, -slower1 -supper1 -snucleotide1 -sprotein1 -scircular1, -sid1 name, -sdbname1 db, -squick1, -iquery1 query, -ioffset1 bytes.")] string? qualifiers = null,
        [Description("The sequence output format (-osformat): fasta (default), embl, genbank, gcg, nbrf, ig, raw, msf, clustal, phylip, nexus, ... every format of seqret.")] string outputFormat = "fasta",
        CancellationToken cancellationToken = default,
        [Description("More files by name (file name -> text) that the program may read (a file named by a USA, @file).")] IReadOnlyDictionary<string, string>? files = null,
        [Description("The USA of the sequence input (-sequence) as the program reads it, instead of the text: [format::]file[begin:end:r], file:id, @list, asis::SEQUENCE.")] string? sequenceValue = null)
    {
        if (string.IsNullOrWhiteSpace(sequence))
            throw new ArgumentException("The sequence cannot be empty", nameof(sequence));
        cancellationToken.ThrowIfCancellationRequested();
        var all = new Dictionary<string, string>(files ?? new Dictionary<string, string>());
        var args = new List<string>();
        all["sequence.dat"] = sequence;
        args.Add("-sequence");
        args.Add(SequenceSpec(sequenceValue, sequence, "sequence.dat"));
        if (!string.IsNullOrWhiteSpace(sequenceFormat))
            args.AddRange(["-sformat1", sequenceFormat]);
        AddQualifierTokens(args, qualifiers);
        if (shuffle.HasValue)
            args.AddRange(["-shuffle", shuffle.Value.ToString(System.Globalization.CultureInfo.InvariantCulture)]);
        return RunEmbossSequenceProgram("shuffleseq", args, all, outputFormat, nameof(sequence), microseconds.HasValue ? new EmbossRunEnvironment(Microseconds: microseconds) : null);
    }
    [McpServerTool(Name = "msbar", Title = "Genome — EMBOSS msbar (mutate sequences)", ReadOnly = true)]
    [Description("EMBOSS 6.6.0 msbar: count times the block, codon (nucleic sequences) and point mutations of the lists are done at random places: insertions, deletions, changes, duplications and moves (1 = a random one of them, 0 = none; a comma list of at most 4); the result is tried again (10 times) when it equals a sequence of otherSequence. Equals the installed program on 300 runs (every input format of seqret, ranges, the output formats, the messages of the library; byte for byte). B09 WP62, F281.")]
    public static EmbossProgramToolResult Msbar(
        [Description("The text of the sequence file (any input format of seqret), or a sequence as letters.")] string sequence,
        [Description("The number of times the mutation operations are done (-count, default 1).")] int? count = null,
        [Description("The point mutations (-point): 0 none, 1 any, 2 insertions, 3 deletions, 4 changes, 5 duplications, 6 moves; a comma list of at most 4 (default 0).")] string? point = null,
        [Description("The block mutations (-block), as point (default 0).")] string? block = null,
        [Description("The codon mutations (-codon), as point (default 0; nucleic sequences only).")] string? codon = null,
        [Description("Do the codon and block operations in frame (-inframe).")] bool inframe = false,
        [Description("The minimum size of a block mutation (-minimum, default 1).")] int? minimum = null,
        [Description("The maximum size of a block mutation (-maximum, default 10).")] int? maximum = null,
        [Description("Sequences the result must not equal (-othersequence: the text of a sequence file, or a USA such as asis:ACGT).")] string? otherSequence = null,
        [Description("The microseconds of the clock when the random number generator is seeded (the seed is microseconds % 9999 + 1; the program takes them from the clock of the machine): the same value gives the same sequences. Default: the clock.")] int? microseconds = null,
        [Description("The input format (-sformat1); default: the program finds the format.")] string? sequenceFormat = null,
        [Description("Sequence input qualifiers as typed after the sequence on the command line: -sbegin1 -send1 -sreverse1, -slower1 -supper1 -snucleotide1 -sprotein1 -scircular1, -sid1 name, -sdbname1 db, -squick1, -iquery1 query, -ioffset1 bytes.")] string? qualifiers = null,
        [Description("The sequence output format (-osformat): fasta (default), embl, genbank, gcg, nbrf, ig, raw, msf, clustal, phylip, nexus, ... every format of seqret.")] string outputFormat = "fasta",
        CancellationToken cancellationToken = default,
        [Description("More files by name (file name -> text) that the program may read (a file named by a USA, @file).")] IReadOnlyDictionary<string, string>? files = null,
        [Description("The USA of the sequence input (-sequence) as the program reads it, instead of the text: [format::]file[begin:end:r], file:id, @list, asis::SEQUENCE.")] string? sequenceValue = null)
    {
        if (string.IsNullOrWhiteSpace(sequence))
            throw new ArgumentException("The sequence cannot be empty", nameof(sequence));
        cancellationToken.ThrowIfCancellationRequested();
        var all = new Dictionary<string, string>(files ?? new Dictionary<string, string>());
        var args = new List<string>();
        all["sequence.dat"] = sequence;
        args.Add("-sequence");
        args.Add(SequenceSpec(sequenceValue, sequence, "sequence.dat"));
        if (!string.IsNullOrWhiteSpace(sequenceFormat))
            args.AddRange(["-sformat1", sequenceFormat]);
        AddQualifierTokens(args, qualifiers);
        if (count.HasValue)
            args.AddRange(["-count", count.Value.ToString(System.Globalization.CultureInfo.InvariantCulture)]);
        if (!string.IsNullOrEmpty(point))
            args.AddRange(["-point", point]);
        if (!string.IsNullOrEmpty(block))
            args.AddRange(["-block", block]);
        if (!string.IsNullOrEmpty(codon))
            args.AddRange(["-codon", codon]);
        if (inframe)
            args.Add("-inframe");
        if (minimum.HasValue)
            args.AddRange(["-minimum", minimum.Value.ToString(System.Globalization.CultureInfo.InvariantCulture)]);
        if (maximum.HasValue)
            args.AddRange(["-maximum", maximum.Value.ToString(System.Globalization.CultureInfo.InvariantCulture)]);
        if (!string.IsNullOrEmpty(otherSequence))
        {
            all["other.dat"] = otherSequence;
            args.AddRange(["-othersequence", otherSequence.StartsWith("asis:", StringComparison.Ordinal) ? otherSequence : "other.dat"]);
        }
        return RunEmbossSequenceProgram("msbar", args, all, outputFormat, nameof(sequence), microseconds.HasValue ? new EmbossRunEnvironment(Microseconds: microseconds) : null);
    }
    [McpServerTool(Name = "makenucseq", Title = "Genome — EMBOSS makenucseq (random nucleotide sequences)", ReadOnly = true)]
    [Description("EMBOSS 6.6.0 makenucseq: amount sequences of length lower case nucleotides, each residue picked at random from ACGT (or from the codons of the codon usage file with the frequencies of the file, the end cut to the length), with the insert put at start; the names are EMBOSS_001, EMBOSS_002 .... Equals the installed program on 300 runs (every input format of seqret, ranges, the output formats, the messages of the library; byte for byte). B09 WP62, F281.")]
    public static EmbossProgramToolResult Makenucseq(
        [Description("The number of sequences (-amount, default 100).")] int? amount = null,
        [Description("The length of every sequence (-length, default 100).")] int? length = null,
        [Description("A string put into every sequence (-insert).")] string? insert = null,
        [Description("The position of the insert (-start, default 1).")] int? start = null,
        [Description("The text of a codon usage file in the EMBOSS (cusp) format (-codonfile).")] string? codonFile = null,
        [Description("The microseconds of the clock when the random number generator is seeded (the seed is microseconds % 9999 + 1; the program takes them from the clock of the machine): the same value gives the same sequences. Default: the clock.")] int? microseconds = null,
        [Description("The sequence output format (-osformat): fasta (default), embl, genbank, gcg, nbrf, ig, raw, msf, clustal, phylip, nexus, ... every format of seqret.")] string outputFormat = "fasta",
        CancellationToken cancellationToken = default)
    {
        cancellationToken.ThrowIfCancellationRequested();
        var all = new Dictionary<string, string>();
        var args = new List<string>();
        if (amount.HasValue)
            args.AddRange(["-amount", amount.Value.ToString(System.Globalization.CultureInfo.InvariantCulture)]);
        if (length.HasValue)
            args.AddRange(["-length", length.Value.ToString(System.Globalization.CultureInfo.InvariantCulture)]);
        if (!string.IsNullOrEmpty(insert))
            args.AddRange(["-insert", insert]);
        if (start.HasValue)
            args.AddRange(["-start", start.Value.ToString(System.Globalization.CultureInfo.InvariantCulture)]);
        if (!string.IsNullOrEmpty(codonFile))
        {
            all["codon.cut"] = codonFile;
            args.AddRange(["-codonfile", "codon.cut"]);
        }
        return RunEmbossSequenceProgram("makenucseq", args, all, outputFormat, nameof(amount), microseconds.HasValue ? new EmbossRunEnvironment(Microseconds: microseconds) : null);
    }
    [McpServerTool(Name = "makeprotseq", Title = "Genome — EMBOSS makeprotseq (random protein sequences)", ReadOnly = true)]
    [Description("EMBOSS 6.6.0 makeprotseq: amount sequences of length lower case residues, each picked at random from the 20 amino acids (or with the composition of the first sequence of a pepstats file), with the insert put at start; the names are EMBOSS_001, EMBOSS_002 .... Equals the installed program on 300 runs (every input format of seqret, ranges, the output formats, the messages of the library; byte for byte). B09 WP62, F281.")]
    public static EmbossProgramToolResult Makeprotseq(
        [Description("The number of sequences (-amount, default 100).")] int? amount = null,
        [Description("The length of every sequence (-length, default 100).")] int? length = null,
        [Description("A string put into every sequence (-insert).")] string? insert = null,
        [Description("The position of the insert (-start, default 1).")] int? start = null,
        [Description("The text of a pepstats output file (-pepstatsfile).")] string? pepstatsFile = null,
        [Description("The microseconds of the clock when the random number generator is seeded (the seed is microseconds % 9999 + 1; the program takes them from the clock of the machine): the same value gives the same sequences. Default: the clock.")] int? microseconds = null,
        [Description("The sequence output format (-osformat): fasta (default), embl, genbank, gcg, nbrf, ig, raw, msf, clustal, phylip, nexus, ... every format of seqret.")] string outputFormat = "fasta",
        CancellationToken cancellationToken = default)
    {
        cancellationToken.ThrowIfCancellationRequested();
        var all = new Dictionary<string, string>();
        var args = new List<string>();
        if (amount.HasValue)
            args.AddRange(["-amount", amount.Value.ToString(System.Globalization.CultureInfo.InvariantCulture)]);
        if (length.HasValue)
            args.AddRange(["-length", length.Value.ToString(System.Globalization.CultureInfo.InvariantCulture)]);
        if (!string.IsNullOrEmpty(insert))
            args.AddRange(["-insert", insert]);
        if (start.HasValue)
            args.AddRange(["-start", start.Value.ToString(System.Globalization.CultureInfo.InvariantCulture)]);
        if (!string.IsNullOrEmpty(pepstatsFile))
        {
            all["pepstats.txt"] = pepstatsFile;
            args.AddRange(["-pepstatsfile", "pepstats.txt"]);
        }
        return RunEmbossSequenceProgram("makeprotseq", args, all, outputFormat, nameof(amount), microseconds.HasValue ? new EmbossRunEnvironment(Microseconds: microseconds) : null);
    }
    [McpServerTool(Name = "nthseqset", Title = "Genome — EMBOSS nthseqset (write one set of sequences)", ReadOnly = true)]
    [Description("EMBOSS 6.6.0 nthseqset: the sets of the input (a file of several alignments: phylip, nexus ... or one set of all the sequences of the file) are counted and the number-th is written; with fewer sets the library reports 'No such sequence set - only N sequence sets were input.'. Equals the installed program on 300 runs (every input format of seqret, ranges, the output formats, the messages of the library; byte for byte). B09 WP62, F281.")]
    public static EmbossProgramToolResult Nthseqset(
        [Description("The text of the sequence file (any input format of seqret), or a sequence as letters.")] string sequence,
        [Description("The number of the set (-number, default 1).")] int? number = null,
        [Description("The input format (-sformat1); default: the program finds the format.")] string? sequenceFormat = null,
        [Description("Sequence input qualifiers as typed after the sequence on the command line: -sbegin1 -send1 -sreverse1, -slower1 -supper1 -snucleotide1 -sprotein1 -scircular1, -sid1 name, -sdbname1 db, -squick1, -iquery1 query, -ioffset1 bytes.")] string? qualifiers = null,
        [Description("The sequence output format (-osformat): fasta (default), embl, genbank, gcg, nbrf, ig, raw, msf, clustal, phylip, nexus, ... every format of seqret.")] string outputFormat = "fasta",
        CancellationToken cancellationToken = default,
        [Description("More files by name (file name -> text) that the program may read (a file named by a USA, @file).")] IReadOnlyDictionary<string, string>? files = null,
        [Description("The USA of the sequence input (-sequence) as the program reads it, instead of the text: [format::]file[begin:end:r], file:id, @list, asis::SEQUENCE.")] string? sequenceValue = null)
    {
        if (string.IsNullOrWhiteSpace(sequence))
            throw new ArgumentException("The sequence cannot be empty", nameof(sequence));
        cancellationToken.ThrowIfCancellationRequested();
        var all = new Dictionary<string, string>(files ?? new Dictionary<string, string>());
        var args = new List<string>();
        all["sequence.dat"] = sequence;
        args.Add("-sequence");
        args.Add(SequenceSpec(sequenceValue, sequence, "sequence.dat"));
        if (!string.IsNullOrWhiteSpace(sequenceFormat))
            args.AddRange(["-sformat1", sequenceFormat]);
        AddQualifierTokens(args, qualifiers);
        if (number.HasValue)
            args.AddRange(["-number", number.Value.ToString(System.Globalization.CultureInfo.InvariantCulture)]);
        return RunEmbossSequenceProgram("nthseqset", args, all, outputFormat, nameof(sequence), null);
    }
    [McpServerTool(Name = "listor", Title = "Genome — EMBOSS listor (logical OR of two sets of sequences)", ReadOnly = true)]
    [Description("EMBOSS 6.6.0 listor: the USAs of the sequences of the first set and the second set that are (OR: all, AND: in both, XOR: in one only, NOT: in the first and not in the second); sequences are equal when their residues are equal ignoring the case. Equals the installed program on 300 runs (every input format of seqret, ranges, the output formats, the messages of the library; byte for byte). B09 WP62, F281.")]
    public static EmbossProgramToolResult Listor(
        [Description("The text of the sequence file (any input format of seqret), or a sequence as letters.")] string firstSequences,
        [Description("The text of the second sequence file (any input format of seqret), or a sequence as letters.")] string secondSequences = "",
        [Description("OR, AND, XOR or NOT (-operator, default OR; a unique abbreviation is accepted).")] string? @operator = null,
        [Description("The input format (-sformat1); default: the program finds the format.")] string? sequenceFormat = null,
        [Description("Sequence input qualifiers as typed after the sequence on the command line: -sbegin1 -send1 -sreverse1, -slower1 -supper1 -snucleotide1 -sprotein1 -scircular1, -sid1 name, -sdbname1 db, -squick1, -iquery1 query, -ioffset1 bytes.")] string? qualifiers = null,
        CancellationToken cancellationToken = default,
        [Description("More files by name (file name -> text) that the program may read (a file named by a USA, @file).")] IReadOnlyDictionary<string, string>? files = null)
    {
        if (string.IsNullOrWhiteSpace(firstSequences))
            throw new ArgumentException("The sequence cannot be empty", nameof(firstSequences));
        cancellationToken.ThrowIfCancellationRequested();
        var all = new Dictionary<string, string>(files ?? new Dictionary<string, string>());
        var args = new List<string>();
        all["firstsequences.dat"] = firstSequences;
        args.Add("-firstsequences");
        args.Add(SequenceSpec(null, firstSequences, "firstsequences.dat"));
        all["secondsequences.dat"] = secondSequences;
        args.Add("-secondsequences");
        args.Add(SequenceSpec(null, secondSequences, "secondsequences.dat"));
        if (!string.IsNullOrWhiteSpace(sequenceFormat))
            args.AddRange(["-sformat1", sequenceFormat]);
        AddQualifierTokens(args, qualifiers);
        if (!string.IsNullOrEmpty(@operator))
            args.AddRange(["-operator", @operator]);
        return RunEmbossOutfileProgram("listor", args, all, nameof(firstSequences), null);
    }
    [McpServerTool(Name = "infoseq", Title = "Genome — EMBOSS infoseq (basic information about sequences)", ReadOnly = true)]
    [Description("EMBOSS 6.6.0 infoseq: for every sequence the USA, database, name, accession, GI, version, type, length, percent GC (nucleic), organism and description, in columns (or delimited, or as an HTML table), with a heading line. Equals the installed program on 300 runs (every input format of seqret, ranges, the output formats, the messages of the library; byte for byte). B09 WP62, F281.")]
    public static EmbossProgramToolResult Infoseq(
        [Description("The text of the sequence file (any input format of seqret), or a sequence as letters.")] string sequence,
        [Description("Print the information in aligned columns (-columns, default true).")] bool columns = true,
        [Description("The delimiter of the records without columns (-delimiter, default |).")] string? delimiter = null,
        [Description("Format the output as an HTML table (-html).")] bool html = false,
        [Description("Display only the columns that are asked for (-only).")] bool only = false,
        [Description("Display the column headings (-heading, default true).")] bool heading = true,
        [Description("The USA column (-usa).")] bool usa = true,
        [Description("The database column (-database).")] bool database = true,
        [Description("The name column (-name).")] bool name = true,
        [Description("The accession column (-accession).")] bool accession = true,
        [Description("The GI column (-gi).")] bool gi = false,
        [Description("The version column (-seqversion).")] bool seqversion = false,
        [Description("The type column (-type).")] bool @type = true,
        [Description("The length column (-length).")] bool length = true,
        [Description("The percent GC column (-pgc).")] bool pgc = true,
        [Description("The organism column (-organism).")] bool organism = true,
        [Description("The description column (-description).")] bool description = true,
        [Description("The input format (-sformat1); default: the program finds the format.")] string? sequenceFormat = null,
        [Description("Sequence input qualifiers as typed after the sequence on the command line: -sbegin1 -send1 -sreverse1, -slower1 -supper1 -snucleotide1 -sprotein1 -scircular1, -sid1 name, -sdbname1 db, -squick1, -iquery1 query, -ioffset1 bytes.")] string? qualifiers = null,
        CancellationToken cancellationToken = default,
        [Description("More files by name (file name -> text) that the program may read (a file named by a USA, @file).")] IReadOnlyDictionary<string, string>? files = null,
        [Description("The USA of the sequence input (-sequence) as the program reads it, instead of the text: [format::]file[begin:end:r], file:id, @list, asis::SEQUENCE.")] string? sequenceValue = null)
    {
        if (string.IsNullOrWhiteSpace(sequence))
            throw new ArgumentException("The sequence cannot be empty", nameof(sequence));
        cancellationToken.ThrowIfCancellationRequested();
        var all = new Dictionary<string, string>(files ?? new Dictionary<string, string>());
        var args = new List<string>();
        all["sequence.dat"] = sequence;
        args.Add("-sequence");
        args.Add(SequenceSpec(sequenceValue, sequence, "sequence.dat"));
        if (!string.IsNullOrWhiteSpace(sequenceFormat))
            args.AddRange(["-sformat1", sequenceFormat]);
        AddQualifierTokens(args, qualifiers);
        if (!columns)
            args.Add("-nocolumns");
        if (!string.IsNullOrEmpty(delimiter))
            args.AddRange(["-delimiter", delimiter]);
        if (html)
            args.Add("-html");
        if (only)
            args.Add("-only");
        if (!heading)
            args.Add("-noheading");
        if (!usa)
            args.Add("-nousa");
        if (!database)
            args.Add("-nodatabase");
        if (!name)
            args.Add("-noname");
        if (!accession)
            args.Add("-noaccession");
        if (gi)
            args.Add("-gi");
        if (seqversion)
            args.Add("-seqversion");
        if (!@type)
            args.Add("-notype");
        if (!length)
            args.Add("-nolength");
        if (!pgc)
            args.Add("-nopgc");
        if (!organism)
            args.Add("-noorganism");
        if (!description)
            args.Add("-nodescription");
        return RunEmbossOutfileProgram("infoseq", args, all, nameof(sequence), null);
    }
    [McpServerTool(Name = "wordcount", Title = "Genome — EMBOSS wordcount (count the words of sequences)", ReadOnly = true)]
    [Description("EMBOSS 6.6.0 wordcount: the words of wordSize residues (default 4 for nucleic, 2 for protein sequences; words with the ambiguity code N / X are skipped) of all the sequences with their counts, the most frequent first (equal counts in the order of the hash table of the library), at least minCount. Equals the installed program on 300 runs (every input format of seqret, ranges, the output formats, the messages of the library; byte for byte). B09 WP62, F281.")]
    public static EmbossProgramToolResult Wordcount(
        [Description("The text of the sequence file (any input format of seqret), or a sequence as letters.")] string sequence,
        [Description("The size of a word (-wordsize).")] int? wordSize = null,
        [Description("The minimum count to report (-mincount, default 1).")] int? minCount = null,
        [Description("The input format (-sformat1); default: the program finds the format.")] string? sequenceFormat = null,
        [Description("Sequence input qualifiers as typed after the sequence on the command line: -sbegin1 -send1 -sreverse1, -slower1 -supper1 -snucleotide1 -sprotein1 -scircular1, -sid1 name, -sdbname1 db, -squick1, -iquery1 query, -ioffset1 bytes.")] string? qualifiers = null,
        CancellationToken cancellationToken = default,
        [Description("More files by name (file name -> text) that the program may read (a file named by a USA, @file).")] IReadOnlyDictionary<string, string>? files = null,
        [Description("The USA of the sequence input (-sequence) as the program reads it, instead of the text: [format::]file[begin:end:r], file:id, @list, asis::SEQUENCE.")] string? sequenceValue = null)
    {
        if (string.IsNullOrWhiteSpace(sequence))
            throw new ArgumentException("The sequence cannot be empty", nameof(sequence));
        cancellationToken.ThrowIfCancellationRequested();
        var all = new Dictionary<string, string>(files ?? new Dictionary<string, string>());
        var args = new List<string>();
        all["sequence.dat"] = sequence;
        args.Add("-sequence");
        args.Add(SequenceSpec(sequenceValue, sequence, "sequence.dat"));
        if (!string.IsNullOrWhiteSpace(sequenceFormat))
            args.AddRange(["-sformat1", sequenceFormat]);
        AddQualifierTokens(args, qualifiers);
        if (wordSize.HasValue)
            args.AddRange(["-wordsize", wordSize.Value.ToString(System.Globalization.CultureInfo.InvariantCulture)]);
        if (minCount.HasValue)
            args.AddRange(["-mincount", minCount.Value.ToString(System.Globalization.CultureInfo.InvariantCulture)]);
        return RunEmbossOutfileProgram("wordcount", args, all, nameof(sequence), null);
    }
    [McpServerTool(Name = "seqxref", Title = "Genome — EMBOSS seqxref (database cross-references of sequence entries)", ReadOnly = true)]
    [Description("EMBOSS 6.6.0 seqxref: for every entry '#USA: count' and the cross-references db:id (the DR lines of EMBL and Swiss-Prot entries ...) sorted by database name and identifier, equal ones once. Equals the installed program on 300 runs (every input format of seqret, ranges, the output formats, the messages of the library; byte for byte). B09 WP62, F281.")]
    public static EmbossProgramToolResult Seqxref(
        [Description("The text of the sequence file (any input format of seqret), or a sequence as letters.")] string sequence,
        [Description("The input format (-sformat1); default: the program finds the format.")] string? sequenceFormat = null,
        [Description("Sequence input qualifiers as typed after the sequence on the command line: -sbegin1 -send1 -sreverse1, -slower1 -supper1 -snucleotide1 -sprotein1 -scircular1, -sid1 name, -sdbname1 db, -squick1, -iquery1 query, -ioffset1 bytes.")] string? qualifiers = null,
        CancellationToken cancellationToken = default,
        [Description("More files by name (file name -> text) that the program may read (a file named by a USA, @file).")] IReadOnlyDictionary<string, string>? files = null,
        [Description("The USA of the sequence input (-sequence) as the program reads it, instead of the text: [format::]file[begin:end:r], file:id, @list, asis::SEQUENCE.")] string? sequenceValue = null)
    {
        if (string.IsNullOrWhiteSpace(sequence))
            throw new ArgumentException("The sequence cannot be empty", nameof(sequence));
        cancellationToken.ThrowIfCancellationRequested();
        var all = new Dictionary<string, string>(files ?? new Dictionary<string, string>());
        var args = new List<string>();
        all["sequence.dat"] = sequence;
        args.Add("-sequence");
        args.Add(SequenceSpec(sequenceValue, sequence, "sequence.dat"));
        if (!string.IsNullOrWhiteSpace(sequenceFormat))
            args.AddRange(["-sformat1", sequenceFormat]);
        AddQualifierTokens(args, qualifiers);
        return RunEmbossOutfileProgram("seqxref", args, all, nameof(sequence), null);
    }
    [McpServerTool(Name = "showorf", Title = "Genome — EMBOSS showorf (nucleotide sequence and its translations)", ReadOnly = true)]
    [Description("EMBOSS 6.6.0 showorf: the sequence in blocks of width bases with a ruler and the numbers of the bases, and under it the translations of the frames (F1 F2 F3 forward, R1 R2 R3 reverse: the amino acid under the first base of every codon, the numbers of the first and last residue of the block); the genetic code is the EMBOSS table of the code number. Equals the installed program on 300 runs (every input format of seqret, ranges, the output formats, the messages of the library; byte for byte). B09 WP62, F281.")]
    public static EmbossProgramToolResult Showorf(
        [Description("The text of the sequence file (the first sequence is used), or a sequence as letters.")] string sequence,
        [Description("The frames to translate (-frames): a comma list of 0 none, 1 F1, 2 F2, 3 F3, 4 R1, 5 R2, 6 R3 (default 1,2,3,4,5,6).")] string? frames = null,
        [Description("The genetic code (-table): 0 .. 6, 9 .. 16, 21 .. 23.")] string? table = null,
        [Description("The width of the screen (-width, default 50; at least 10).")] int? width = null,
        [Description("Add a ruler (-ruler, default true).")] bool ruler = true,
        [Description("Number the translations (-plabel, default true).")] bool plabel = true,
        [Description("Number the sequence (-nlabel, default true).")] bool nlabel = true,
        [Description("The input format (-sformat1); default: the program finds the format.")] string? sequenceFormat = null,
        [Description("Sequence input qualifiers as typed after the sequence on the command line: -sbegin1 -send1 -sreverse1, -slower1 -supper1 -snucleotide1 -sprotein1 -scircular1, -sid1 name, -sdbname1 db, -squick1, -iquery1 query, -ioffset1 bytes.")] string? qualifiers = null,
        CancellationToken cancellationToken = default,
        [Description("More files by name (file name -> text) that the program may read (a file named by a USA, @file).")] IReadOnlyDictionary<string, string>? files = null,
        [Description("The USA of the sequence input (-sequence) as the program reads it, instead of the text: [format::]file[begin:end:r], file:id, @list, asis::SEQUENCE.")] string? sequenceValue = null)
    {
        if (string.IsNullOrWhiteSpace(sequence))
            throw new ArgumentException("The sequence cannot be empty", nameof(sequence));
        cancellationToken.ThrowIfCancellationRequested();
        var all = new Dictionary<string, string>(files ?? new Dictionary<string, string>());
        var args = new List<string>();
        all["sequence.dat"] = sequence;
        args.Add("-sequence");
        args.Add(SequenceSpec(sequenceValue, sequence, "sequence.dat"));
        if (!string.IsNullOrWhiteSpace(sequenceFormat))
            args.AddRange(["-sformat1", sequenceFormat]);
        AddQualifierTokens(args, qualifiers);
        if (!string.IsNullOrEmpty(frames))
            args.AddRange(["-frames", frames]);
        if (!string.IsNullOrEmpty(table))
            args.AddRange(["-table", table]);
        if (width.HasValue)
            args.AddRange(["-width", width.Value.ToString(System.Globalization.CultureInfo.InvariantCulture)]);
        if (!ruler)
            args.Add("-noruler");
        if (!plabel)
            args.Add("-noplabel");
        if (!nlabel)
            args.Add("-nonlabel");
        return RunEmbossOutfileProgram("showorf", args, all, nameof(sequence), null);
    }
    [McpServerTool(Name = "prettyseq", Title = "Genome — EMBOSS prettyseq (nucleotide sequence and its translation)", ReadOnly = true)]
    [Description("EMBOSS 6.6.0 prettyseq: the sequence in blocks of width bases with a ruler and numbers; the bases of the ranges are in lower case and the amino acid is written under the first base of each of their codons (the rest of the sequence in upper case with blanks under it). Equals the installed program on 300 runs (every input format of seqret, ranges, the output formats, the messages of the library; byte for byte). B09 WP62, F281.")]
    public static EmbossProgramToolResult Prettyseq(
        [Description("The text of the sequence file (the first sequence is used), or a sequence as letters.")] string sequence,
        [Description("The ranges to translate (-range '4-57,78-94'; default the whole sequence).")] string? range = null,
        [Description("The genetic code (-table): 0 .. 6, 9 .. 16, 21 .. 23.")] string? table = null,
        [Description("Add a ruler (-ruler, default true).")] bool ruler = true,
        [Description("Number the translations (-plabel, default true).")] bool plabel = true,
        [Description("Number the sequence (-nlabel, default true).")] bool nlabel = true,
        [Description("The width of the screen (-width, default 60; at least 10).")] int? width = null,
        [Description("The input format (-sformat1); default: the program finds the format.")] string? sequenceFormat = null,
        [Description("Sequence input qualifiers as typed after the sequence on the command line: -sbegin1 -send1 -sreverse1, -slower1 -supper1 -snucleotide1 -sprotein1 -scircular1, -sid1 name, -sdbname1 db, -squick1, -iquery1 query, -ioffset1 bytes.")] string? qualifiers = null,
        CancellationToken cancellationToken = default,
        [Description("More files by name (file name -> text) that the program may read (a file named by a USA, @file).")] IReadOnlyDictionary<string, string>? files = null,
        [Description("The USA of the sequence input (-sequence) as the program reads it, instead of the text: [format::]file[begin:end:r], file:id, @list, asis::SEQUENCE.")] string? sequenceValue = null)
    {
        if (string.IsNullOrWhiteSpace(sequence))
            throw new ArgumentException("The sequence cannot be empty", nameof(sequence));
        cancellationToken.ThrowIfCancellationRequested();
        var all = new Dictionary<string, string>(files ?? new Dictionary<string, string>());
        var args = new List<string>();
        all["sequence.dat"] = sequence;
        args.Add("-sequence");
        args.Add(SequenceSpec(sequenceValue, sequence, "sequence.dat"));
        if (!string.IsNullOrWhiteSpace(sequenceFormat))
            args.AddRange(["-sformat1", sequenceFormat]);
        AddQualifierTokens(args, qualifiers);
        if (!string.IsNullOrEmpty(range))
            args.AddRange(["-range", range]);
        if (!string.IsNullOrEmpty(table))
            args.AddRange(["-table", table]);
        if (!ruler)
            args.Add("-noruler");
        if (!plabel)
            args.Add("-noplabel");
        if (!nlabel)
            args.Add("-nonlabel");
        if (width.HasValue)
            args.AddRange(["-width", width.Value.ToString(System.Globalization.CultureInfo.InvariantCulture)]);
        return RunEmbossOutfileProgram("prettyseq", args, all, nameof(sequence), null);
    }
}
