// SPDX-License-Identifier: GPL-2.0-or-later
// This file is part of Seqeron.Mcp.EmbossFeatureTools, licensed under the terms stated in LICENSE (this directory).
using System.ComponentModel;
using System.Text;
using ModelContextProtocol.Server;
using Seqeron.Mcp.Analysis.Tools;
using static Seqeron.Mcp.Analysis.Tools.AnalysisTools;

namespace Seqeron.Mcp.EmbossFeatureTools.Tools;

public partial class EmbossFeatureToolsTools
{
    [McpServerTool(Name = "featcopy", Title = "Alignment — EMBOSS featcopy (read and write feature tables)", ReadOnly = true)]
    [Description("EMBOSS 6.6.0 featcopy: every feature table of a file in ANY input format is read by the feature readers of the EMBOSS library and written in an output format. Input formats: embl, genbank, refseq (nucleotide, feature lines of five blanks), genpept, refseqp (protein), gff2 (versions 1 / 2, ##Type, ##sequence-region, group fields, featflags), gff3 (ID / Parent grouping, several lines of one ID, percent-encoded values, several tables), the gff3 EMBOSS wrote before 6.4, swiss, pir, dasgff (XML); with no format the library tries embl, refseq, refseqp, genpept, gff3, gff2, swiss, pir and dasgff in turn on the start of the file and keeps the first that finds features, and a 'format::' prefix of features forces one (an unknown name is an error message and the formats are tried); the value may also be a LIST: '@listfile' or 'list::listfile' (files given in the files argument; a word of a line is a file, 'format::file', 'file[begin:end]' or a list again, '#' starts a comment, a line of several words is 'Bad record', a missing file or list is an error message; the queries of a nested list are read before the pending ones; 500 runs equal the installed program, the lists that name themselves, on which the program never stops, are refused), and 'asis::' is refused for features as the program refuses it. Output formats written for a feature table: gff3 (gff, unknown input = refused as the program does), gff2, embl (em), genbank (gb, ddbj, refseq), bed (track line, score 0, blocks of the sub-features), refseqp, pir (nbrf), swiss (sw, swissprot, uniprot: FT lines wrapped at 75 with the /FTId= rule), dasgff (XML, type category from the ontology table), draw (cirdna / lindna group), debug (the full internal trace: type, source, score, strand, positions, flags, group, exon, tags, sub-features); a format the program refuses for the kind of sequence (embl for a protein table, pir for a nucleotide one) is an error message and nothing is written. The first table is cut to the range featureBegin / featureEnd (a table with a range is trimmed and its features shifted: a feature cut at the end is flagged). The reader and writer pairs equal the installed program on 3 800 input files (embl, genbank, refseq, genpept, refseqp, gff2, gff3, gff3emboss, swiss, pir, dasgff: 236 or more each, damaged and mixed files, the format:: prefix, ranges, the files on which the program never stops are refused with an error) and 1 800 output files in gff3 / gff2 / embl / genbank and their aliases and 2 750 output files in bed, refseq, refseqp, pir / nbrf, swiss / sw / swissprot / uniprot, dasgff and draw (a Swiss-Prot feature line the library's wrapper cannot break - a word as long as the line - never ends in the program and is refused with an error); the messages (Warning: bad /protein_id value ..., Error: Unknown input format ...) are returned as the program prints them.")]
    public static FeatcopyToolResult Featcopy(
        [Description("The text of the feature file (a feature table of any input format, or an EMBL / GenBank entry; several tables may follow one another).")] string features,
        [Description("The value of -features: 'format::' prefix and [begin:end] suffix are read as the program reads them; default the file itself.")] string? featuresValue = null,
        [Description("The output feature format (-offormat): gff3 (default), gff, gff2, bed, embl, em, genbank, gb, ddbj, refseq, refseqp, pir, nbrf, swiss, sw, swissprot, uniprot, dasgff, draw, debug.")] string outputFormat = "gff3",
        [Description("-fbegin: the first table is cut to this range start (default 1; negative counts from the end).")] int featureBegin = 1,
        [Description("-fend: the first table is cut to this range end (default: its length).")] int? featureEnd = null,
        CancellationToken cancellationToken = default,
        [Description("More files by name (file name -> text) that featuresValue may name: a list file ('@list' or 'list::list': one query per word of a line, '#' starts a comment, the queries are files with an optional 'format::' prefix and [begin:end] suffix, or lists again), the files the list names; a file that is not here does not exist. 'asis::' is refused for features as the program refuses it.")] IReadOnlyDictionary<string, string>? files = null)
    {
        if (string.IsNullOrWhiteSpace(features))
            throw new ArgumentException("The feature file cannot be empty", nameof(features));
        cancellationToken.ThrowIfCancellationRequested();
        string name = "features.dat";
        string value = string.IsNullOrWhiteSpace(featuresValue) ? name : featuresValue;
        var all = new Dictionary<string, string>(files ?? new Dictionary<string, string>()) { [name] = features };
        try
        {
            var o = EmbossFeaturePrograms.Featcopy(new EmbossFeatcopyRequest(value, all, featureBegin, featureEnd, "unknown.gff", outputFormat, null, DateTime.Now));
            return new FeatcopyToolResult(o.FileName, o.Text, o.Format, o.Tables, [.. o.Messages]);
        }
        catch (EmbossFeatureDiedException e)
        {
            throw new ArgumentException(e.Message + (e.Errors.Count > 0 ? " (" + string.Join("; ", e.Errors.Where(m => m.Length > 0)) + ")" : ""), nameof(features));
        }
        catch (Exception e) when (e is NotSupportedException or InvalidOperationException)
        {
            throw new ArgumentException(e.Message, nameof(outputFormat));
        }
    }
    [McpServerTool(Name = "featmerge", Title = "Alignment — EMBOSS featmerge (combine two feature tables)", ReadOnly = true)]
    [Description("EMBOSS 6.6.0 featmerge: the first feature table of featuresA and the first of featuresB (each in any input format of the EMBOSS feature readers: embl, genbank, refseq, genpept, refseqp, gff2, gff3, swiss, pir, dasgff; with no format every format is tried; 'format::' prefix and [begin:end] suffix in the *Value arguments) are read, a copy of every feature of B (every tag checked again, the group numbers kept) is added to A, and A is written in an output format: gff3 (default), gff, gff2, bed, embl, em, genbank, gb, ddbj, refseq, refseqp, pir, nbrf, swiss, sw, swissprot, uniprot, dasgff, draw, debug. The range of A (-fbegin1 / -fend1) cuts and shifts the merged table as it does in the program; the range of B is read and not used. Equals the installed program on 450 runs (every pair of input formats, ranges of both inputs, every output format, missing files, the inputs on which the program never stops refused with an error); the messages (Warning: ..., Error: Failed to open filename ...) are returned as the program prints them.")]
    public static FeatcopyToolResult Featmerge(
        [Description("The text of the first feature file (-afeatures).")] string featuresA,
        [Description("The text of the second feature file (-bfeatures).")] string featuresB,
        [Description("The value of -afeatures: 'format::' prefix and [begin:end] suffix; default the file itself.")] string? featuresAValue = null,
        [Description("The value of -bfeatures.")] string? featuresBValue = null,
        [Description("The output feature format (-offormat): gff3 (default), gff, gff2, bed, embl, em, genbank, gb, ddbj, refseq, refseqp, pir, nbrf, swiss, sw, swissprot, uniprot, dasgff, draw, debug.")] string outputFormat = "gff3",
        [Description("-fbegin1: the range start of the first table (default 1).")] int featureBegin = 1,
        [Description("-fend1: the range end of the first table (default: its length).")] int? featureEnd = null,
        CancellationToken cancellationToken = default,
        [Description("More files by name (file name -> text) that featuresAValue / featuresBValue may name, list files included ('@list').")] IReadOnlyDictionary<string, string>? files = null)
    {
        if (string.IsNullOrWhiteSpace(featuresA))
            throw new ArgumentException("The feature file cannot be empty", nameof(featuresA));
        if (string.IsNullOrWhiteSpace(featuresB))
            throw new ArgumentException("The feature file cannot be empty", nameof(featuresB));
        cancellationToken.ThrowIfCancellationRequested();
        try
        {
            var o = EmbossFeaturePrograms.Featmerge(new EmbossFeatmergeRequest(string.IsNullOrWhiteSpace(featuresAValue) ? "a.dat" : featuresAValue, string.IsNullOrWhiteSpace(featuresBValue) ? "b.dat" : featuresBValue,
                new Dictionary<string, string>(files ?? new Dictionary<string, string>()) { ["a.dat"] = featuresA, ["b.dat"] = featuresB }, featureBegin, featureEnd, 1, null, "unknown.gff", outputFormat, null, DateTime.Now));
            return new FeatcopyToolResult(o.FileName, o.Text, o.Format, o.Tables, [.. o.Messages]);
        }
        catch (EmbossFeatureDiedException e)
        {
            throw new ArgumentException(e.Message + (e.Errors.Count > 0 ? " (" + string.Join("; ", e.Errors.Where(m => m.Length > 0)) + ")" : ""), nameof(featuresA));
        }
        catch (Exception e) when (e is NotSupportedException or InvalidOperationException)
        {
            throw new ArgumentException(e.Message, nameof(outputFormat));
        }
    }
    [McpServerTool(Name = "feattext", Title = "Alignment — EMBOSS feattext (original text of a feature table)", ReadOnly = true)]
    [Description("EMBOSS 6.6.0 feattext: every feature table of a file (any input format of the EMBOSS feature readers; with no format every format is tried; 'format::' prefix and [begin:end] suffix in featuresValue) is read with the input storing the text it reads, and the text of each entry is returned as the library stored it: the EMBL / GenBank reader and the GFF3 reader keep the lines they consumed (a GFF3 table stops at the next ##gff-version line, which is given back with the library's off-by-the-line-ending cut of the stored text), the readers of the other formats keep none (nothing is returned for their tables). Equals the installed program on 400 runs (every input format, damaged files, the format:: prefix, the runs on which the program never stops refused with an error); the messages are returned as the program prints them.")]
    public static FeattextToolResult Feattext(
        [Description("The text of the feature file (a feature table of any input format, or an EMBL / GenBank entry; several tables may follow one another).")] string features,
        [Description("The value of -features: 'format::' prefix and [begin:end] suffix are read as the program reads them; default the file itself.")] string? featuresValue = null,
        CancellationToken cancellationToken = default,
        [Description("More files by name (file name -> text) that featuresValue may name, list files included ('@list').")] IReadOnlyDictionary<string, string>? files = null)
    {
        if (string.IsNullOrWhiteSpace(features))
            throw new ArgumentException("The feature file cannot be empty", nameof(features));
        cancellationToken.ThrowIfCancellationRequested();
        string name = "features.dat";
        try
        {
            var o = EmbossFeaturePrograms.Feattext(new EmbossFeattextRequest(string.IsNullOrWhiteSpace(featuresValue) ? name : featuresValue, new Dictionary<string, string>(files ?? new Dictionary<string, string>()) { [name] = features }));
            return new FeattextToolResult(o.Text, o.Tables, [.. o.Messages]);
        }
        catch (EmbossFeatureDiedException e)
        {
            throw new ArgumentException(e.Message + (e.Errors.Count > 0 ? " (" + string.Join("; ", e.Errors.Where(m => m.Length > 0)) + ")" : ""), nameof(features));
        }
        catch (InvalidOperationException e)
        {
            throw new ArgumentException(e.Message, nameof(features));
        }
    }
    [McpServerTool(Name = "maskfeat", Title = "Alignment — EMBOSS maskfeat (mask the regions of features)", ReadOnly = true)]
    [Description("EMBOSS 6.6.0 maskfeat: the regions of the features of the wanted types are masked in the sequences of a file (any input format of seqret that carries features: embl, genbank, refseq, refseqp, swiss, nbrf, gff2, gff3; the feature types are those of the vocabulary, wildcards * and ? allowed, several types separated by blanks, commas, semicolons or bars; the default repeat* masks repeat_region and repeat_unit): the residues become the mask character (default N for nucleotide and X for protein sequences) or lower case (toLower, or a blank or no mask character), and the sequences are returned in FASTA (the sequence type of the program is any: every character that is no letter and no * is dropped). Equals the installed program on 1 200 runs (every input format with features, damaged files, the runs on which the program refuses an input or dies by a fault). A5-8c3g, B09 F249.")]
    public static MaskfeatToolResult Maskfeat(
        [Description("The text of the sequence file with features (see the formats above), or a sequence as letters (nothing to mask).")] string sequence,
        [Description("The input format (-sformat); default: the program finds the format.")] string? sequenceFormat = null,
        [Description("The feature types to mask (-type): names separated by blanks, commas, semicolons or bars; * and ? wildcards; default repeat*.")] string type = "repeat*",
        [Description("Make the masked regions lower case (-tolower).")] bool toLower = false,
        [Description("The mask character (-maskchar, one character; none or a blank = lower case); default N (nucleotide) / X (protein).")] string? maskChar = null,
        [Description("Sequence input qualifiers as typed after the sequence on the command line: -sbegin1 -send1 -sreverse1, -slower1 -supper1 -snucleotide1 -sprotein1 -scircular1, -sformat1 F, -sid1 name, -sdbname1 db, -squick1, -iquery1 query, -ioffset1 bytes.")] string? qualifiers = null,
        CancellationToken cancellationToken = default,
        [Description("More files by name (file name -> text) that sequence may name with a 'format::' prefix.")] IReadOnlyDictionary<string, string>? files = null,
        [Description("The USA of the sequence input (-sequence) as the program reads it, instead of the file the text is put in: [format::]file[begin:end:r] (file = the name of an entry of files; the text of sequence is the file 'sequence.dat'), file:id / file:field:query / file{a|b} (id, acc, sv, gi, des, key, org; wildcards * and ?), @list / list::list (a list file of files; # comments; every entry has its own range and format), asis::SEQUENCE; the range qualifiers -sbegin1 -send1 -sreverse1 of qualifiers are applied to the sequences and to their features.")] string? sequenceValue = null)
    {
        if (string.IsNullOrWhiteSpace(sequence))
            throw new ArgumentException("The sequence cannot be empty", nameof(sequence));
        cancellationToken.ThrowIfCancellationRequested();
        string name = "sequence.dat";
        var q = EmbossSequenceQualifiers.ParseFile(qualifiers) with { Format = sequenceFormat };
        string spec = SequenceSpec(sequenceValue, sequence, name);
        try
        {
            var all = new Dictionary<string, string>(files ?? new Dictionary<string, string>()) { [name] = sequence };
            var o = EmbossMaskfeatPrograms.Maskfeat(new EmbossMaskfeatRequest(spec, all, q, type, toLower, maskChar));
            if (o.Died)
                throw new ArgumentException(string.Join("; ", o.Messages.Where(m => m.Length > 0)), nameof(sequence));
            return new MaskfeatToolResult(o.SequenceText, o.Sequences, [.. o.Messages]);
        }
        catch (Exception e) when (e is NotSupportedException or InvalidOperationException)
        {
            throw new ArgumentException(e.Message, nameof(sequence));
        }
    }
    [McpServerTool(Name = "extractfeat", Title = "Alignment — EMBOSS extractfeat (extract features as sequences)", ReadOnly = true)]
    [Description("EMBOSS 6.6.0 extractfeat: the features of the sequences of a file (any input format of seqret that carries features: embl, genbank, refseq, refseqp, swiss, nbrf, gff2, gff3) that pass the filters (source, type, sense, score, tag, value; wildcards * and ?, several patterns separated by blanks, commas, semicolons or bars) are extracted as sequences, named SEQUENCE_START_END (with featInName also _TYPE) and described by [type] and the tags of describe, optionally with the sequence before and after the feature (before / after: negative values count from the other end of the feature; the missing sequence at the ends is padded with N or X) and with the parts of a join as one sequence (join). The sequences are returned in FASTA or in the output format of outputFormat (fasta unknown pearson, gcg gcg8, nbrf pir, codata, strider, staden, text plain raw, ig iguspto, gde, fitch, acedb, mase, ncbi, gifasta, nib (nibble), dasdna, das, asn1, embl (em, emblnew, emold, emblold), genbank (gb, ddbj), refseq, genpept, refseqp, swiss (sw, swissprot, uniprot, swissnew, swnew, swold, swissold, swissprotold, swissprotnew), gff2, gff3 (gff), fastq (fastq-sanger, fastq-illumina, fastq-solexa: the quality values of the sequences are written), sam, scf, and the sequence-set formats msf, clustal (aln), selex, phylip, phylipnon (phylip3), hennig86, mega, meganon, nexus (paup), nexusnon (paupnon), jackknifer, jackknifernon, treecon (written when the output is closed), experiment, debug and bam (BGZF blocks deflated as zlib 1.3 does)). A part of a join that names another entry (J00194.1:100..103) is read from the same file or database. Equals the installed program on 1 000 runs (every input format with features, every filter, before / after / join, 29 output format names, damaged files). A5-8c3h, B09 F258.")]
    public static ExtractfeatToolResult Extractfeat(
        [Description("The text of the sequence file with features (see the formats above), or a sequence as letters (no features: nothing is extracted).")] string sequence,
        [Description("The input format (-sformat); default: the program finds the format.")] string? sequenceFormat = null,
        [Description("Sequence before each feature (-before): bases / residues; negative: the start is that far before the END of the feature.")] int before = 0,
        [Description("Sequence after each feature (-after); negative: the end is that far after the START of the feature.")] int after = 0,
        [Description("Write the parts of a join (CDS, mRNA, ...) as one sequence (-join).")] bool join = false,
        [Description("Add the feature type to the sequence name (-featinname).")] bool featInName = false,
        [Description("Tag names to add with their values to the description (-describe): names separated by blanks, commas, semicolons or bars; wildcards.")] string describe = "",
        [Description("Sources of the features to extract (-source); default *.")] string source = "*",
        [Description("Types of the features to extract (-type); default *.")] string type = "*",
        [Description("Sense of the features (-sense): 0 any, 1 forward, -1 reverse.")] int sense = 0,
        [Description("Minimum score (-minscore); with maxScore both 0 the score is not tested.")] float minScore = 0f,
        [Description("Maximum score (-maxscore).")] float maxScore = 0f,
        [Description("Tags a feature must have (-tag); default *.")] string tag = "*",
        [Description("Values the tag must have (-value); default *.")] string value = "*",
        [Description("The sequence output format (-osformat); default fasta.")] string? outputFormat = null,
        [Description("Sequence input qualifiers as typed after the sequence on the command line: -sbegin1 -send1 -sreverse1, -slower1 -supper1 -snucleotide1 -sprotein1 -scircular1, -sformat1 F, -sid1 name, -sdbname1 db, -squick1, -iquery1 query, -ioffset1 bytes.")] string? qualifiers = null,
        CancellationToken cancellationToken = default,
        [Description("More files by name (file name -> text) that sequenceValue may name: list files, other sequence files, the data files of the databases of definitions.")] IReadOnlyDictionary<string, string>? files = null,
        [Description("The USA of the sequence input (-sequence): [format::]file[begin:end:r], file:id, file:field:value, file{a|b}, @list, list::list, asis::TEXT, db:id; the text of sequence is the file 'sequence.dat'.")] string? sequenceValue = null,
        [Description("The text of an EMBOSS configuration file (.embossrc) that defines databases for sequenceValue: DBNAME name [ type: N format: fasta method: direct directory: \"dir\" filename: \"*.fa\" ] (methods direct, emboss and emblcd), RES / RESOURCE, SERVER, INCLUDE (the included file is an entry of files), IF / IFDEF / ELSE / ENDIF, SET name value; the index files of the methods emboss (dbxfasta / dbxflat: NAME.ent, NAME.xid, NAME.pxid ...) and emblcd (dbifasta / dbiflat: division.lkp, entrynam.idx, acnum.trg ...) are entries of files by path, as latin-1 text.")] string? definitions = null,
        [Description("Write every extracted sequence to a file of its own, name.extension in lower case (-ossingle2); the files are returned in singleFiles.")] bool outputSingle = false,
        [Description("The directory of the files of outputSingle (-osdirectory2).")] string? outputDirectory = null,
        [Description("The extension of the files of outputSingle (-osextension2); default: the output format as typed.")] string? outputExtension = null,
        [Description("The database name of the output sequences (-osdbname2): the writers put it in front of the names or into the identifiers.")] string? outputDbName = null)
    {
        if (string.IsNullOrWhiteSpace(sequence))
            throw new ArgumentException("The sequence cannot be empty", nameof(sequence));
        cancellationToken.ThrowIfCancellationRequested();
        string file = "sequence.dat";
        var q = EmbossSequenceQualifiers.ParseFile(qualifiers) with { Format = sequenceFormat };
        string spec = SequenceSpec(sequenceValue, sequence, file);
        try
        {
            var all = new Dictionary<string, string>(files ?? new Dictionary<string, string>()) { [file] = sequence };
            var o = EmbossExtractfeatProgram.Extractfeat(new EmbossExtractfeatRequest(spec, all, q, before, after, join, featInName, describe, source, type, sense, minScore, maxScore, tag, value, "stdout", outputFormat) { Definitions = string.IsNullOrWhiteSpace(definitions) ? null : [definitions], OutSingle = outputSingle, OutDirectory = outputDirectory, OutExtension = outputExtension, OutDbName = outputDbName });
            if (o.Died)
                throw new ArgumentException(string.Join("; ", o.Messages.Where(m => m.Length > 0)), nameof(sequence));
            return new ExtractfeatToolResult(o.SequenceText, o.Sequences, [.. o.Messages], o.SingleFiles.Count > 0 ? o.SingleFiles : null);
        }
        catch (Exception e) when (e is NotSupportedException or InvalidOperationException)
        {
            throw new ArgumentException(e.Message, nameof(sequence));
        }
    }
    [McpServerTool(Name = "coderet", Title = "Alignment — EMBOSS coderet (CDS, mRNA, translations, non-coding)", ReadOnly = true)]
    [Description("EMBOSS 6.6.0 coderet: from the nucleotide sequences of a file with feature tables (embl, genbank, refseq, gff2, gff3, ...) the coding sequences (the CDS features, joins as one sequence; files NAME.cds, sequences ACCESSION_cds_N), the messenger RNAs (mRNA features; NAME.mrna, _mrna_N), the translations (the /translation tag of the CDS features; NAME.prot, _pro_N) and the non-coding regions (what the CDS, mRNA and exon features leave; NAME.noncoding, _noncoding_START) are written, NAME being the first sequence in lower case. The tool returns the log of the program (the counts of every sequence), the files by name and the messages; cds, mrna, translation and nonCoding select the outputs and outputFormat is the sequence output format (fasta unknown pearson, gcg gcg8, nbrf pir, codata, strider, staden, text plain raw, ig iguspto, gde, fitch, acedb, mase, ncbi, gifasta, nib (nibble), dasdna, das, asn1, embl (em, emblnew, emold, emblold), genbank (gb, ddbj), refseq, genpept, refseqp, swiss (sw, swissprot, uniprot, swissnew, swnew, swold, swissold, swissprotold, swissprotnew), gff2, gff3 (gff), fastq (fastq-sanger, fastq-illumina, fastq-solexa: the quality values of the sequences are written), sam, scf, and the sequence-set formats msf, clustal (aln), selex, phylip, phylipnon (phylip3), hennig86, mega, meganon, nexus (paup), nexusnon (paupnon), jackknifer, jackknifernon, treecon (written when the output is closed), experiment, debug and bam (BGZF blocks deflated as zlib 1.3 does)). Equals the installed program on 1 000 runs (EMBL entries with joins, complements, orders, partial ends and translations, the reader inputs with features, damaged files, one to four outputs in 29 format names). A5-8c3h, B09 F259.")]
    public static CoderetToolResult Coderet(
        [Description("The text of the sequence file with features (nucleotide sequences), see the formats above.")] string sequence,
        [Description("The input format (-sformat); default: the program finds the format.")] string? sequenceFormat = null,
        [Description("Write the coding sequences (-cdsoutseq).")] bool cds = true,
        [Description("Write the messenger RNAs (-mrnaoutseq).")] bool mrna = true,
        [Description("Write the translations (-translationoutseq).")] bool translation = true,
        [Description("Write the non-coding regions (-restoutseq).")] bool nonCoding = true,
        [Description("The sequence output format of the four outputs (-osformat3 .. -osformat6); default fasta.")] string? outputFormat = null,
        [Description("Sequence input qualifiers as typed after the sequence on the command line: -sbegin1 -send1 -sreverse1, -slower1 -supper1, -sformat1 F, -sid1 name, -sdbname1 db, -squick1, -iquery1 query, -ioffset1 bytes.")] string? qualifiers = null,
        CancellationToken cancellationToken = default,
        [Description("More files by name (file name -> text) that sequenceValue may name.")] IReadOnlyDictionary<string, string>? files = null,
        [Description("The USA of the sequence input (-seqall): [format::]file[begin:end:r], file:id, file:field:value, file{a|b}, @list, list::list, db:id; the text of sequence is the file 'sequence.dat'.")] string? sequenceValue = null,
        [Description("The text of an EMBOSS configuration file (.embossrc) that defines databases for sequenceValue: DBNAME definitions (methods direct, emboss and emblcd), RES / RESOURCE, SERVER, INCLUDE (the included file is an entry of files), IF / IFDEF / ELSE / ENDIF and SET blocks; the files of a database, and for the methods emboss (dbxfasta / dbxflat index files NAME.ent, .xid, .pxid ...) and emblcd (dbifasta / dbiflat files division.lkp, entrynam.idx, acnum.trg ...) its index files as latin-1 text, are entries of files by path.")] string? definitions = null,
        [Description("Write every sequence of the four outputs to a file of its own, name.extension in lower case (-ossingle3 .. -ossingle6); the files are in files.")] bool outputSingle = false,
        [Description("The directory of the output files (-osdirectory3 .. -osdirectory6).")] string? outputDirectory = null,
        [Description("The extension of the output files (-osextension3 .. -osextension6); default cds, mrna, prot and noncoding.")] string? outputExtension = null,
        [Description("The base name of the output files (-osname3 .. -osname6); default: the name of the first sequence in lower case.")] string? outputName = null)
    {
        if (string.IsNullOrWhiteSpace(sequence))
            throw new ArgumentException("The sequence cannot be empty", nameof(sequence));
        cancellationToken.ThrowIfCancellationRequested();
        string file = "sequence.dat";
        var q = EmbossSequenceQualifiers.ParseFile(qualifiers) with { Format = sequenceFormat };
        string spec = SequenceSpec(sequenceValue, sequence, file);
        try
        {
            var all = new Dictionary<string, string>(files ?? new Dictionary<string, string>()) { [file] = sequence };
            var o = EmbossCoderetProgram.Coderet(new EmbossCoderetRequest(spec, all, q, cds ? null : "", mrna ? null : "", translation ? null : "", nonCoding ? null : "", outputFormat, outputFormat, outputFormat, outputFormat)
            {
                Definitions = string.IsNullOrWhiteSpace(definitions) ? null : [definitions],
                CdsOptions = new EmbossSeqoutOptions(outputName, outputExtension, outputDirectory, outputSingle),
                MrnaOptions = new EmbossSeqoutOptions(outputName, outputExtension, outputDirectory, outputSingle),
                TranslationOptions = new EmbossSeqoutOptions(outputName, outputExtension, outputDirectory, outputSingle),
                RestOptions = new EmbossSeqoutOptions(outputName, outputExtension, outputDirectory, outputSingle),
            });
            if (o.Died)
                throw new ArgumentException(string.Join("; ", o.Messages.Where(m => m.Length > 0)), nameof(sequence));
            return new CoderetToolResult(o.Report, o.Files, [.. o.Messages]);
        }
        catch (Exception e) when (e is NotSupportedException or InvalidOperationException)
        {
            throw new ArgumentException(e.Message, nameof(sequence));
        }
    }
    [McpServerTool(Name = "showfeat", Title = "Alignment — EMBOSS showfeat (display features)", ReadOnly = true)]
    [Description("EMBOSS 6.6.0 showfeat: the features of the sequences of a file (any input format of seqret that carries features) as lines of ASCII graphics under a scale: one line per feature (or per part of a join; with joinFeatures the parts of a join on one line, with collapse the features on the same line) with the strand (|---> for forward, <---| for reverse), optionally the strand, source, positions, type and tags with their values. The features are filtered by sourceMatch, typeMatch, tagMatch and valueMatch (wildcards * and ?, several patterns separated by blanks, commas, semicolons or bars; stricttags shows only the tags that match), sorted (sort: start, source, type or nosort; by strand first), the sense is selected by forward / reverse / unknown, and annotation adds ranges with a text (\"5-10 foo,20-30 bar\") as lines of their own. width is the width of the graphics (below 4 no graphics and no scale), html wraps the output in HTML. Equals the installed program on 1 000 runs (every input format with features, every option, damaged files, the runs on which the program dies by a fault). A5-8c3h, B09 F260.")]
    public static ShowfeatToolResult Showfeat(
        [Description("The text of the sequence file with features (see the formats above), or a sequence as letters (no features).")] string sequence,
        [Description("The input format (-sformat); default: the program finds the format.")] string? sequenceFormat = null,
        [Description("Sources of the features to show (-sourcematch); default *.")] string sourceMatch = "*",
        [Description("Types of the features to show (-typematch); default *.")] string typeMatch = "*",
        [Description("Tags a feature must have (-tagmatch); default *.")] string tagMatch = "*",
        [Description("Values the tag must have (-valuematch); default *.")] string valueMatch = "*",
        [Description("Sort order (-sort): start (default), source, type or nosort.")] string sort = "start",
        [Description("Show the parts of a join on one line (-joinfeatures).")] bool joinFeatures = false,
        [Description("Output as HTML (-html).")] bool html = false,
        [Description("Show the name of the sequence (-id); default true.")] bool id = true,
        [Description("Show the description of the sequence (-description); default true.")] bool description = true,
        [Description("Show the scale line (-scale); default true.")] bool scale = true,
        [Description("Width of the graphics lines (-width); default 60; below 4 there are no graphics and no scale; negative values are reset to 0 with a warning.")] int width = 60,
        [Description("Show all features of a type on one line (-collapse).")] bool collapse = false,
        [Description("Show forward features (-forward); default true.")] bool forward = true,
        [Description("Show reverse features (-reverse); default true.")] bool reverse = true,
        [Description("Show features of unknown sense (-unknown); default true.")] bool unknown = true,
        [Description("Show the strand of the features (-strand).")] bool showStrand = false,
        [Description("Show the source (-sourceshow).")] bool showSource = false,
        [Description("Show the positions (-position).")] bool showPosition = false,
        [Description("Show the type (-typeshow); default true.")] bool showType = true,
        [Description("Show the tags (-tagshow).")] bool showTags = false,
        [Description("Show the values of the tags (-valueshow); default true.")] bool showValues = true,
        [Description("Show only the tags that match tagMatch and valueMatch (-stricttags).")] bool strictTags = false,
        [Description("Ranges to show as annotation (-annotation): \"5-10 text,20-30 text\".")] string annotation = "",
        [Description("Sequence input qualifiers as typed after the sequence on the command line: -sbegin1 -send1 -sreverse1, -slower1 -supper1 -snucleotide1 -sprotein1 -scircular1, -sformat1 F, -sid1 name, -sdbname1 db, -squick1, -iquery1 query, -ioffset1 bytes.")] string? qualifiers = null,
        CancellationToken cancellationToken = default,
        [Description("More files by name (file name -> text) that sequenceValue may name.")] IReadOnlyDictionary<string, string>? files = null,
        [Description("The USA of the sequence input (-sequence): [format::]file[begin:end:r], file:id, file:field:value, file{a|b}, @list, list::list, asis::TEXT, db:id; the text of sequence is the file 'sequence.dat'.")] string? sequenceValue = null,
        [Description("The text of an EMBOSS configuration file (.embossrc) that defines databases for sequenceValue: DBNAME definitions (methods direct, emboss and emblcd), RES / RESOURCE, SERVER, INCLUDE (the included file is an entry of files), IF / IFDEF / ELSE / ENDIF and SET blocks; the files of a database, and for the methods emboss (dbxfasta / dbxflat index files NAME.ent, .xid, .pxid ...) and emblcd (dbifasta / dbiflat files division.lkp, entrynam.idx, acnum.trg ...) its index files as latin-1 text, are entries of files by path.")] string? definitions = null)
    {
        if (string.IsNullOrWhiteSpace(sequence))
            throw new ArgumentException("The sequence cannot be empty", nameof(sequence));
        cancellationToken.ThrowIfCancellationRequested();
        string file = "sequence.dat";
        var q = EmbossSequenceQualifiers.ParseFile(qualifiers) with { Format = sequenceFormat };
        string spec = SequenceSpec(sequenceValue, sequence, file);
        try
        {
            var all = new Dictionary<string, string>(files ?? new Dictionary<string, string>()) { [file] = sequence };
            var o = EmbossShowfeatProgram.Showfeat(new EmbossShowfeatRequest(spec, all, q)
            {
                SourceMatch = sourceMatch, TypeMatch = typeMatch, TagMatch = tagMatch, ValueMatch = valueMatch, Sort = sort, JoinFeatures = joinFeatures, Html = html, Id = id, Description = description, Scale = scale, Width = width,
                Collapse = collapse, Forward = forward, Reverse = reverse, Unknown = unknown, Strand = showStrand, SourceShow = showSource, Position = showPosition, TypeShow = showType, TagShow = showTags, ValueShow = showValues,
                StrictTags = strictTags, Annotation = annotation, Definitions = string.IsNullOrWhiteSpace(definitions) ? null : [definitions],
            });
            if (o.Died)
                throw new ArgumentException(string.Join("; ", o.Messages.Where(m => m.Length > 0)), nameof(sequence));
            return new ShowfeatToolResult(o.Text, o.Sequences, [.. o.Messages]);
        }
        catch (Exception e) when (e is NotSupportedException or InvalidOperationException or ArgumentException)
        {
            throw new ArgumentException(e.Message, nameof(sequence));
        }
    }
    [McpServerTool(Name = "featreport", Title = "Alignment — EMBOSS featreport (report a feature table)", ReadOnly = true)]
    [Description("EMBOSS 6.6.0 featreport: the first feature table of a feature file (embl, genbank, refseq, refseqp, swiss, pir, gff2, gff3, dasgff, or any: the library tries every reader) is reported for a sequence in a report format (reportFormat: gff default, gff3, embl, genbank, pir, swiss, dasgff, debug, listfile, dbmotif, diffseq, draw, restrict, excel, feattable, motif, nametable, regions, seqtable, simple, srs, summary, table, tagseq, or the report qualifiers as typed: -rformat table -raccshow -rdesshow -rscoreshow N -rstrandshow N -rusashow -rmaxall 3 -rmaxseq 2). The formats gff, gff3, embl, genbank, pir, swiss and dasgff are written by the feature writers of the library (joins, tags, flags); the others list the start, end, strand and score of the features. The features come with an optional range (featureBegin, featureEnd). Equals the installed program on 1 000 runs (every feature format, every report format, ranges, damaged files). A5-8c3h, B09 F262. The report file name qualifiers of reportFormat (-rname -rextension -rdirectory: the default file name is the name of the sequence + .featreport, a directory is put in front of a name that has none, a directory that does not exist stops the program) and the qualifiers of the features input in qualifiers (-fformat2 -fopenfile2 -iquery2 -ioffset2 -fask2 -fbegin2 -fend2 -freverse2 -fcircular2; only -fbegin2 / -fend2 change the report, as in the program) are applied as the program applies them; the file name is returned in reportFile. A5-8c3h, B09 F278.")]
    public static FeatreportToolResult Featreport(
        [Description("The text of the sequence file (any input format of seqret), or a sequence as letters.")] string sequence,
        [Description("The text of the feature file (see the formats above).")] string features,
        [Description("The input format of the sequence (-sformat); default: the program finds the format.")] string? sequenceFormat = null,
        [Description("The report format name or the report qualifiers as typed (see above); default gff.")] string reportFormat = "gff",
        [Description("Sequence input qualifiers as typed after the sequence on the command line: -sbegin1 -send1 -sreverse1, -slower1 -supper1 -snucleotide1 -sprotein1 -scircular1, -sformat1 F, -sid1 name, -sdbname1 db, -squick1, -iquery1 query, -ioffset1 bytes.")] string? qualifiers = null,
        [Description("The first position of the features to report (-fbegin2); default 1.")] int featureBegin = 1,
        [Description("The last position of the features to report (-fend2); default the length of the table.")] int? featureEnd = null,
        CancellationToken cancellationToken = default,
        [Description("More files by name (file name -> text) that sequenceValue and featuresValue may name.")] IReadOnlyDictionary<string, string>? files = null,
        [Description("The USA of the sequence input (-sequence): [format::]file[begin:end:r], file:id, file:field:value, file{a|b}, @list, list::list, asis::TEXT, db:id; the text of sequence is the file 'sequence.dat'.")] string? sequenceValue = null,
        [Description("The value of -features: [format::]file[begin:end:r]; the text of features is the file 'features.dat'.")] string? featuresValue = null,
        [Description("The text of an EMBOSS configuration file (.embossrc) that defines databases for sequenceValue: DBNAME definitions (methods direct, emboss and emblcd), RES / RESOURCE, SERVER, INCLUDE (the included file is an entry of files), IF / IFDEF / ELSE / ENDIF and SET blocks; the files of a database, and for the methods emboss (dbxfasta / dbxflat index files NAME.ent, .xid, .pxid ...) and emblcd (dbifasta / dbiflat files division.lkp, entrynam.idx, acnum.trg ...) its index files as latin-1 text, are entries of files by path.")] string? definitions = null,
        [Description("The name of the report file (-outfile); default: standard output, or - when reportFormat has -rname / -rextension / -rdirectory - the file named after the sequence (name.featreport) in that directory; it is returned in reportFile.")] string? outputFile = null)
    {
        if (string.IsNullOrWhiteSpace(sequence))
            throw new ArgumentException("The sequence cannot be empty", nameof(sequence));
        if (string.IsNullOrWhiteSpace(features))
            throw new ArgumentException("The features cannot be empty", nameof(features));
        cancellationToken.ThrowIfCancellationRequested();
        string file = "sequence.dat", featFile = "features.dat";
        string format = (reportFormat ?? "").Trim();
        string reportText = format.StartsWith('-') ? format : "-rformat " + format;
        string spec = SequenceSpec(sequenceValue, sequence, file);
        try
        {
            var (inputs, rq) = EmbossSequenceQualifiers.ParseInputsAndReport("featreport", (qualifiers ?? "") + " " + reportText);
            bool named = rq.Name is not null || rq.Extension is not null || rq.Directory is not null;
            bool defaultName = outputFile is null && named;
            var commandLine = new List<string> { EmbossCommandLineEcho.Parameter("sequence", spec), EmbossCommandLineEcho.Parameter("features", featuresValue ?? featFile) };
            if (!defaultName)
                commandLine.Add(EmbossCommandLineEcho.Parameter("outfile", outputFile ?? "stdout"));
            commandLine.AddRange(NumberedReportEntries(reportText, 3));
            commandLine.AddRange(HeaderEntries(qualifiers));
            var all = new Dictionary<string, string>(files ?? new Dictionary<string, string>()) { [file] = sequence, [featFile] = features };
            var o = EmbossFeatreportProgram.Featreport(new EmbossFeatreportRequest(spec, featuresValue ?? featFile, all, inputs[0] with { Format = sequenceFormat ?? inputs[0].Format }, rq, featureBegin, featureEnd, commandLine, DateTime.Now, outputFile ?? "stdout") { Definitions = string.IsNullOrWhiteSpace(definitions) ? null : [definitions], DefaultReportFile = defaultName });
            if (o.Died)
                throw new ArgumentException(string.Join("; ", o.Messages.Where(m => m.Length > 0)), nameof(sequence));
            return new FeatreportToolResult(o.Report, [.. o.Messages], o.ReportFile);
        }
        catch (Exception e) when (e is NotSupportedException or InvalidOperationException)
        {
            throw new ArgumentException(e.Message, nameof(sequence));
        }
    }
    [McpServerTool(Name = "twofeat", Title = "Alignment — EMBOSS twofeat (neighbouring pairs of features)", ReadOnly = true)]
    [Description("EMBOSS 6.6.0 twofeat: the features of the sequences of a file (any input format of seqret that carries features) are tested against the criteria of feature A (aSource, aType, aSense 0 / + / -, aMinScore, aMaxScore, aTag, aValue) and of feature B (bSource ... bValue); wildcards * and ?, several patterns separated by blanks, commas, semicolons or bars; the scores are cut to integers) and every pair of an A and a B feature in the required relation is reported: overlap (A any, O overlap, NO none, NW overlap but not within, AW A within B, BW B within A), distance between minRange and maxRange (ignored unless min is less than max; measured from the nearest, left, right or furthest ends: rangeType N / L / R / F), sense (A any, S same, O opposite: the opposite sense stops the program with a fatal error of the library), as a new feature (type typeOut, source twofeat, covering both, tags startA endA startB endB) or, with twoOut, as the two features. The report is in reportFormat (table default; the names of featreport or the report qualifiers as typed). Equals the installed program on 1 000 runs (every input format with features, every option, every report format). A5-8c3h, B09 F263. The report file name qualifiers of reportFormat (-rname -rextension -rdirectory: the default file name is the name of the sequence + .twofeat) are applied as the program applies them; the file name is returned in reportFile. A5-8c3h, B09 F278.")]
    public static TwofeatToolResult Twofeat(
        [Description("The text of the sequence file with features (any input format of seqret), or a sequence as letters (no features).")] string sequence,
        [Description("The input format (-sformat); default: the program finds the format.")] string? sequenceFormat = null,
        [Description("Source of feature A (-asource); default *.")] string aSource = "*",
        [Description("Type of feature A (-atype); default *.")] string aType = "*",
        [Description("Sense of feature A (-asense): 0 any, + forward ('-' is read as any by the program).")] string aSense = "0",
        [Description("Minimum score of feature A (-aminscore).")] float aMinScore = 0f,
        [Description("Maximum score of feature A (-amaxscore).")] float aMaxScore = 0f,
        [Description("Tag of feature A (-atag); default *.")] string aTag = "*",
        [Description("Value of the tags of feature A (-avalue); default *.")] string aValue = "*",
        [Description("Source of feature B (-bsource); default *.")] string bSource = "*",
        [Description("Type of feature B (-btype); default *.")] string bType = "*",
        [Description("Sense of feature B (-bsense).")] string bSense = "0",
        [Description("Minimum score of feature B (-bminscore).")] float bMinScore = 0f,
        [Description("Maximum score of feature B (-bmaxscore).")] float bMaxScore = 0f,
        [Description("Tag of feature B (-btag); default *.")] string bTag = "*",
        [Description("Value of the tags of feature B (-bvalue); default *.")] string bValue = "*",
        [Description("Overlap (-overlap): A any, O overlap, NO none, NW overlap but not within, AW A within B, BW B within A.")] string overlap = "A",
        [Description("Minimum distance between the features (-minrange).")] int minRange = 0,
        [Description("Maximum distance between the features (-maxrange).")] int maxRange = 0,
        [Description("Where the distance is measured (-rangetype): N nearest ends, L left ends, R right ends, F furthest ends.")] string rangeType = "N",
        [Description("Sense of the pair (-sense): A any, S same, O opposite (fatal in the program).")] string sense = "A",
        [Description("Order of the features (-order): accepted, not used by the program.")] string order = "A",
        [Description("Report the two features of a pair instead of one feature that covers them (-twoout).")] bool twoOut = false,
        [Description("The type of the new feature (-typeout); default misc_feature.")] string typeOut = "misc_feature",
        [Description("The report format name or the report qualifiers as typed (see featreport); default table.")] string reportFormat = "table",
        [Description("Sequence input qualifiers as typed after the sequence on the command line: -sbegin1 -send1 -sreverse1, -slower1 -supper1 -snucleotide1 -sprotein1 -scircular1, -sformat1 F, -sid1 name, -sdbname1 db, -squick1, -iquery1 query, -ioffset1 bytes.")] string? qualifiers = null,
        CancellationToken cancellationToken = default,
        [Description("More files by name (file name -> text) that sequenceValue may name.")] IReadOnlyDictionary<string, string>? files = null,
        [Description("The USA of the sequence input (-sequence): [format::]file[begin:end:r], file:id, file:field:value, file{a|b}, @list, list::list, asis::TEXT, db:id; the text of sequence is the file 'sequence.dat'.")] string? sequenceValue = null,
        [Description("The text of an EMBOSS configuration file (.embossrc) that defines databases for sequenceValue: DBNAME definitions (methods direct, emboss and emblcd), RES / RESOURCE, SERVER, INCLUDE (the included file is an entry of files), IF / IFDEF / ELSE / ENDIF and SET blocks; the files of a database, and for the methods emboss (dbxfasta / dbxflat index files NAME.ent, .xid, .pxid ...) and emblcd (dbifasta / dbiflat files division.lkp, entrynam.idx, acnum.trg ...) its index files as latin-1 text, are entries of files by path.")] string? definitions = null,
        [Description("The name of the report file (-outfile); default: standard output, or - when reportFormat has -rname / -rextension / -rdirectory - the file named after the sequence (name.twofeat) in that directory; it is returned in reportFile.")] string? outputFile = null)
    {
        if (string.IsNullOrWhiteSpace(sequence))
            throw new ArgumentException("The sequence cannot be empty", nameof(sequence));
        cancellationToken.ThrowIfCancellationRequested();
        string file = "sequence.dat";
        string format = (reportFormat ?? "").Trim();
        string reportText = format.StartsWith('-') ? format : "-rformat " + format;
        string spec = SequenceSpec(sequenceValue, sequence, file);
        try
        {
            var (inputs, rq) = EmbossSequenceQualifiers.ParseInputsAndReport("twofeat", (qualifiers ?? "") + " " + reportText);
            bool named = rq.Name is not null || rq.Extension is not null || rq.Directory is not null;
            bool defaultName = outputFile is null && named;
            var commandLine = new List<string> { EmbossCommandLineEcho.Parameter("sequence", spec) };
            if (!defaultName)
                commandLine.Add(EmbossCommandLineEcho.Parameter("outfile", outputFile ?? "stdout"));
            commandLine.AddRange(NumberedReportEntries(reportText, 2));
            commandLine.AddRange(HeaderEntries(qualifiers));
            var all = new Dictionary<string, string>(files ?? new Dictionary<string, string>()) { [file] = sequence };
            var o = EmbossTwofeatProgram.Twofeat(new EmbossTwofeatRequest(spec, all, inputs[0] with { Format = sequenceFormat ?? inputs[0].Format }, rq, commandLine, DateTime.Now, outputFile ?? "stdout")
            {
                DefaultReportFile = defaultName,
                ASource = aSource, AType = aType, ASense = aSense, AMinScore = aMinScore, AMaxScore = aMaxScore, ATag = aTag, AValue = aValue,
                BSource = bSource, BType = bType, BSense = bSense, BMinScore = bMinScore, BMaxScore = bMaxScore, BTag = bTag, BValue = bValue,
                Overlap = overlap, MinRange = minRange, MaxRange = maxRange, RangeType = rangeType, Sense = sense, Order = order, TwoOut = twoOut, TypeOut = typeOut,
                Definitions = string.IsNullOrWhiteSpace(definitions) ? null : [definitions],
            });
            if (o.Died)
                throw new ArgumentException(string.Join("; ", o.Messages.Where(m => m.Length > 0)), nameof(sequence));
            return new TwofeatToolResult(o.Report, [.. o.Messages], o.ReportFile);
        }
        catch (Exception e) when (e is NotSupportedException or InvalidOperationException)
        {
            throw new ArgumentException(e.Message, nameof(sequence));
        }
    }
}
