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
    [McpServerTool(Name = "seqret", Title = "Alignment — EMBOSS seqret (read and write sequences)", ReadOnly = true)]
    [Description("EMBOSS 6.6.0 seqret: the sequences of a file are read by the sequence readers of the EMBOSS library (ajseqread.c) and written in FASTA or, with outputFormat, in another sequence format (with the GFF3 or another feature format of the features the file carries). Input formats read byte for byte as the installed program reads them (200 inputs of each of 66 names, with the library's warnings and errors): fasta (ncbi, pearson, gifasta, dbid), gcg (gcg8), embl (em), swiss (sw, swissprot, uniprot; the DE description of RecName / AltName / SubName / Includes / Contains / Flags is assembled as the library does), nbrf (pir), genbank (gb, ddbj, refseq), refseqp, genpept, codata, strider, ig, igstrict, iguspto, staden, textonly (plain, asis, unknown), fastq, fastq-sanger, fastq-illumina, fastq-solexa, sam, gff2, gff3 (gff; the sequence after ##FASTA / the ##DNA lines), gde, raw, and the alignment formats clustal (aln), msf, fitch, mase, biomart, treecon, acedb, phylip (interleaved; falls back to phylipnon as the library does), phylipnon, nexus (paup), mega, selex, stockholm (pfam), jackknifer, hennig86 (the library loops over memory on a sequence line: refused), the entries of experiment, pdb (the ATOM residues), pdbseq (SEQRES), pdbnuc, pdbnucseq (the members of the alignment are returned one after the other). ace / consed (the library reads the AS line and then dies by a segmentation fault) are refused, as are the few inputs on which the library reads freed or foreign memory (a repeated taxon label of a nexus file, ...). With no format every format the program tries is tried in turn (a text no format reads is read as plain text). The features of embl, genbank, refseq, refseqp, swiss, nbrf, gff2 and gff3 files are returned with featureFormat when feature is true. Output: fasta unknown pearson, gcg gcg8, nbrf pir, codata, strider, staden, text plain raw, ig iguspto, gde, fitch, acedb, mase, ncbi, gifasta, nib (nibble), dasdna, das, asn1, embl (em, emblnew, emold, emblold), genbank (gb, ddbj), refseq, genpept, refseqp, swiss (sw, swissprot, uniprot, swissnew, swnew, swold, swissold, swissprotold, swissprotnew), gff2, gff3 (gff), fastq (fastq-sanger, fastq-illumina, fastq-solexa: the quality values of the sequences are written), sam, scf, and the sequence-set formats msf, clustal (aln), selex, phylip, phylipnon (phylip3), hennig86, mega, meganon, nexus (paup), nexusnon (paupnon), jackknifer, jackknifernon, treecon (written when the output is closed), experiment, debug and bam (BGZF blocks deflated as zlib 1.3 does) - written byte for byte as the installed program writes them (200 or more runs of each name: the format as -osformat2 or the 'format::' prefix of the output name, outputSingle with outputDirectory / outputExtension, outputDbName (-osdbname), unknown formats, the format-for-sequence-type errors; the entry formats print the references, keywords, taxonomy, comments, cross-references and dates the readers keep, and with feature true the formats that keep their features (embl, genbank, refseq, refseqp, swiss, nbrf / pir, gff2, gff3) write them inside the sequence file as the installed program does and no feature file is made; a set format with outputSingle dies in the installed program on most inputs and is refused). Not read yet (the binary abi, bam, nibble, scf, das, ensembl and the CO contig lines of embl): the tool refuses such an input. A5-8c3g / A5-8c3h, B09 F248 / F251 - F252 / F265 - F271.")]
    public static SeqretToolResult Seqret(
        [Description("The text of the sequence file (any input format above), or a sequence as letters. A value 'asis::ACGT' is the sequence itself.")] string sequence,
        [Description("The input format (-sformat): the names above, any case, or an EDAM term; default: the program finds the format.")] string? sequenceFormat = null,
        [Description("Read and write the features of the sequences (-feature).")] bool feature = false,
        [Description("Read one sequence and stop (-firstonly).")] bool firstOnly = false,
        [Description("The feature output format (-offormat2): gff3 (default), gff, gff2, bed, embl, em, genbank, gb, ddbj, refseq, refseqp, pir, nbrf, swiss, sw, swissprot, uniprot, dasgff, draw, debug.")] string featureFormat = "gff3",
        [Description("Sequence input qualifiers as typed after the sequence on the command line: -sbegin1 -send1 -sreverse1 (the range and the strand: applied to the sequences and to their features), -slower1 -supper1 -snucleotide1 -sprotein1 -scircular1, -sformat1 F, -sid1 name, -sdbname1 db, -squick1, -iquery1 query, -ioffset1 bytes.")] string? qualifiers = null,
        CancellationToken cancellationToken = default,
        [Description("More files by name (file name -> text) that sequence may name with a 'format::' prefix.")] IReadOnlyDictionary<string, string>? files = null,
        [Description("The sequence output format (-osformat2): fasta (default; unknown, pearson), gcg (gcg8), nbrf (pir), codata, strider, staden, text (plain, raw), ig, gde, fitch, acedb, mase, ncbi, gifasta, nib (nibble), dasdna, das, asn1, embl (em, emblnew, emold, emblold), genbank (gb, ddbj), refseq, genpept, refseqp, swiss (sw, swissprot, uniprot, swissnew, swnew, swold, swissold, swissprotold, swissprotnew), gff2, gff3 (gff), fastq, fastq-sanger, fastq-illumina, fastq-solexa, sam, scf, msf, clustal (aln), selex, phylip, phylipnon (phylip3), hennig86, mega, meganon, nexus (paup), nexusnon (paupnon), jackknifer, jackknifernon, treecon. iguspto is accepted by the library but its writer dies on every sequence that has no full description (all but swiss and iguspto input), so such a run is refused. experiment, debug and bam are not written yet.")] string outputFormat = "fasta",
        [Description("Write every sequence to a file of its own, name.extension in lower case (-ossingle2); the files are returned in singleFiles.")] bool outputSingle = false,
        [Description("The directory of the files of outputSingle (-osdirectory2).")] string? outputDirectory = null,
        [Description("The extension of the files of outputSingle (-osextension2); default: the output format as typed.")] string? outputExtension = null,
        [Description("The database name of the output sequences (-osdbname2): the writers put it in front of the names (db:name: fasta, fastq, gde) or into the identifiers (gnl|db|name: ncbi, gifasta).")] string? outputDbName = null,
        [Description("The USA of the sequence input (-sequence) as the program reads it, instead of the file the text is put in: [format::]file[begin:end:r] (file = the name of an entry of files; the text of sequence is the file 'sequence.dat'), file:id / file:field:query / file{a|b} (id, acc, sv, gi, des, key, org; wildcards * and ?), @list / list::list (a list file of files; # comments; every entry has its own range and format), asis::SEQUENCE; the range qualifiers -sbegin1 -send1 -sreverse1 of qualifiers are applied to the sequences and to their features.")] string? sequenceValue = null,
        [Description("The text of an EMBOSS configuration file (.embossrc) that defines databases for sequenceValue: DBNAME name [ type: N format: fasta method: direct directory: \"dir\" filename: \"*.fa\" ], SET name value; the data files are in files.")] string? definitions = null)
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
            var o = EmbossSeqretPrograms.Seqret(new EmbossSeqretRequest(spec, all, q, feature, firstOnly, "stdout", outputFormat, "outfeat.gff", featureFormat, DateTime.Now) { OutSingle = outputSingle, OutDirectory = outputDirectory, OutExtension = outputExtension, OutDbName = outputDbName, CommandLine = $"-auto -sequence {spec} -outseq stdout -osformat2 {outputFormat}" + (feature ? " -feature" : ""), Definitions = string.IsNullOrWhiteSpace(definitions) ? null : [definitions] });
            if (o.Died)
                throw new ArgumentException(string.Join("; ", o.Messages.Where(m => m.Length > 0)), nameof(sequence));
            return new SeqretToolResult(o.SequenceText, o.FeatureText, o.Sequences, [.. o.Messages], o.SingleFiles.Count > 0 ? o.SingleFiles : null);
        }
        catch (Exception e) when (e is NotSupportedException or InvalidOperationException)
        {
            throw new ArgumentException(e.Message, nameof(sequence));
        }
    }
    [McpServerTool(Name = "union", Title = "Alignment — EMBOSS union (join sequences into one)", ReadOnly = true)]
    [Description("EMBOSS 6.6.0 union: the sequences of a file (any input format of seqret) are joined into ONE sequence that carries the name, description, accession and type of the first; with findOverlap the bases by which a sequence repeats the end of the one before (the longest end of the previous sequence that is the start of this one) are joined once, and with overlapFile the overlap of every sequence is returned (one number per sequence, the first is 0). The result is written in FASTA or in the output format of outputFormat (fasta unknown pearson, gcg gcg8, nbrf pir, codata, strider, staden, text plain raw, ig iguspto, gde, fitch, acedb, mase, ncbi, gifasta, nib (nibble), dasdna, das, asn1, embl (em, emblnew, emold, emblold), genbank (gb, ddbj), refseq, genpept, refseqp, swiss (sw, swissprot, uniprot, swissnew, swnew, swold, swissold, swissprotold, swissprotnew), gff2, gff3 (gff), fastq (fastq-sanger, fastq-illumina, fastq-solexa: the quality values of the sequences are written), sam, scf, and the sequence-set formats msf, clustal (aln), selex, phylip, phylipnon (phylip3), hennig86, mega, meganon, nexus (paup), nexusnon (paupnon), jackknifer, jackknifernon, treecon (written when the output is closed), experiment, debug and bam (BGZF blocks deflated as zlib 1.3 does)). With feature the features of the sequences are joined into the features of the result (the positions moved by the length of the sequences before; with source one feature of the type source, tag origid = the name of the sequence, for every sequence joined) and returned as featureText in the feature format of featureFormat (gff3 default, gff2, embl, genbank, refseq, bed, dasgff, debug; an unknown name is written as gff3 with the library's errors). Equals the installed program on 200 runs (fragments that overlap by random lengths, the reader inputs, damaged files, every output format, the format-for-sequence-type errors). A5-8c3h, B09 F253. The joined sequence can be written to a file of its own (outputSingle with outputDirectory / outputExtension: -ossingle2 -osdirectory2 -osextension2; returned in singleFiles), as seqret does; equals the installed program on 1 200 runs (every output format including experiment, debug and bam, the output qualifiers, directories that do not exist). A5-8c3h, B09 F277.")]
    public static UnionToolResult Union(
        [Description("The text of the sequence file (any input format of seqret), or a sequence as letters.")] string sequence,
        [Description("The input format (-sformat); default: the program finds the format.")] string? sequenceFormat = null,
        [Description("Look for the overlaps of consecutive sequences when joining (-findoverlap).")] bool findOverlap = false,
        [Description("Return the overlap of every sequence (-overlapfile).")] bool overlapFile = false,
        [Description("The sequence output format (-osformat); default fasta.")] string outputFormat = "fasta",
        [Description("Sequence input qualifiers as typed after the sequence on the command line: -sbegin1 -send1 -sreverse1, -slower1 -supper1 -snucleotide1 -sprotein1 -scircular1, -sformat1 F, -sid1 name, -sdbname1 db, -squick1, -iquery1 query, -ioffset1 bytes.")] string? qualifiers = null,
        CancellationToken cancellationToken = default,
        [Description("More files by name (file name -> text) that sequence may name with a 'format::' prefix.")] IReadOnlyDictionary<string, string>? files = null,
        [Description("The USA of the sequence input (-sequence) as the program reads it, instead of the file the text is put in: [format::]file[begin:end:r] (file = the name of an entry of files; the text of sequence is the file 'sequence.dat'), file:id / file:field:query / file{a|b} (id, acc, sv, gi, des, key, org; wildcards * and ?), @list / list::list (a list file of files; # comments; every entry has its own range and format), asis::SEQUENCE; the range qualifiers -sbegin1 -send1 -sreverse1 of qualifiers are applied to the sequences and to their features.")] string? sequenceValue = null,
        [Description("Join the features of the sequences (-feature); the features are returned as featureText.")] bool feature = false,
        [Description("Create a source feature for every sequence joined (-source; with feature).")] bool source = false,
        [Description("The format of the feature output (-offormat2); default gff3 (gff).")] string? featureFormat = null,
        [Description("Write the joined sequence to a file of its own, name.extension in lower case (-ossingle2); the file is returned in singleFiles.")] bool outputSingle = false,
        [Description("The directory of the file of outputSingle (-osdirectory2).")] string? outputDirectory = null,
        [Description("The extension of the file of outputSingle (-osextension2); default: the output format as typed.")] string? outputExtension = null)
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
            var o = EmbossUnionProgram.Union(new EmbossUnionRequest(spec, all, q, findOverlap, overlapFile, "stdout", outputFormat, DateTime.Now) { Feature = feature, Source = source, OutFeatFormat = featureFormat, OutSingle = outputSingle, OutDirectory = outputDirectory, OutExtension = outputExtension });
            if (o.Died)
                throw new ArgumentException(string.Join("; ", o.Messages.Where(m => m.Length > 0)), nameof(sequence));
            return new UnionToolResult(o.SequenceText, o.OverlapText, [.. o.Messages], o.FeatureText, o.SingleFiles);
        }
        catch (Exception e) when (e is NotSupportedException or InvalidOperationException)
        {
            throw new ArgumentException(e.Message, nameof(sequence));
        }
    }
    [McpServerTool(Name = "seqretsplit", Title = "Alignment — EMBOSS seqretsplit (write every sequence to its own file)", ReadOnly = true)]
    [Description("EMBOSS 6.6.0 seqretsplit: the sequences of a file (any input format of seqret) are written each to a file of its own, name.extension in lower case (the extension is the output format as typed, or outputExtension; outputDirectory is the directory), in FASTA or in the output format of outputFormat (the formats of the union tool). A sequence whose name an earlier sequence used (in any case) gets the name name.002, name.003 ... and the library warns 'Duplicate name ... changed to ...'. With outputSingle false the sequences go to one file (sequenceText), as seqret writes them. Equals the installed program on 200 runs (duplicate names in any case, the reader inputs, damaged files, every output format, directories and extensions). A5-8c3h, B09 F254.")]
    public static SeqretsplitToolResult Seqretsplit(
        [Description("The text of the sequence file (any input format of seqret), or a sequence as letters.")] string sequence,
        [Description("The input format (-sformat); default: the program finds the format.")] string? sequenceFormat = null,
        [Description("Read one sequence and stop (-firstonly).")] bool firstOnly = false,
        [Description("The sequence output format (-osformat); default fasta.")] string outputFormat = "fasta",
        [Description("Write every sequence to a file of its own (-ossingle2; the default of the program).")] bool outputSingle = true,
        [Description("The directory of the files (-osdirectory2).")] string? outputDirectory = null,
        [Description("The extension of the files (-osextension2); default: the output format as typed.")] string? outputExtension = null,
        [Description("Sequence input qualifiers as typed after the sequence on the command line: -sbegin1 -send1 -sreverse1, -slower1 -supper1 -snucleotide1 -sprotein1 -scircular1, -sformat1 F, -sid1 name, -sdbname1 db, -squick1, -iquery1 query, -ioffset1 bytes.")] string? qualifiers = null,
        CancellationToken cancellationToken = default,
        [Description("More files by name (file name -> text) that sequence may name with a 'format::' prefix.")] IReadOnlyDictionary<string, string>? files = null,
        [Description("The USA of the sequence input (-sequence) as the program reads it, instead of the file the text is put in: [format::]file[begin:end:r] (file = the name of an entry of files; the text of sequence is the file 'sequence.dat'), file:id / file:field:query / file{a|b} (id, acc, sv, gi, des, key, org; wildcards * and ?), @list / list::list (a list file of files; # comments; every entry has its own range and format), asis::SEQUENCE; the range qualifiers -sbegin1 -send1 -sreverse1 of qualifiers are applied to the sequences and to their features.")] string? sequenceValue = null,
        [Description("The text of an EMBOSS configuration file (.embossrc) that defines databases for sequenceValue: DBNAME name [ type: N format: fasta method: direct directory: \"dir\" filename: \"*.fa\" ], SET name value; the data files are in files.")] string? definitions = null)
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
            var o = EmbossSeqretPrograms.Seqretsplit(new EmbossSeqretRequest(spec, all, q, false, firstOnly, "stdout", outputFormat, null, null, DateTime.Now) { OutSingle = outputSingle, OutDirectory = outputDirectory, OutExtension = outputExtension, Definitions = string.IsNullOrWhiteSpace(definitions) ? null : [definitions] });
            if (o.Died)
                throw new ArgumentException(string.Join("; ", o.Messages.Where(m => m.Length > 0)), nameof(sequence));
            return new SeqretsplitToolResult(o.SequenceText, o.SingleFiles, o.Sequences, [.. o.Messages]);
        }
        catch (Exception e) when (e is NotSupportedException or InvalidOperationException)
        {
            throw new ArgumentException(e.Message, nameof(sequence));
        }
    }
    [McpServerTool(Name = "entret", Title = "Alignment — EMBOSS entret (retrieve flatfile entries)", ReadOnly = true)]
    [Description("EMBOSS 6.6.0 entret: the entries of a sequence file (any input format of seqret) are read with the text of every entry stored, and the text of each entry is returned exactly as the file has it (the lines the sequence reader read: the whole alignment for the first sequence of an alignment format, nothing for the later ones). sequenceValue is the USA of the program: a query on the file (file:id, file:acc:value, file:wild*, file{a|b}), a range, a list file or a database defined in definitions (db:id); firstOnly stops after the first entry. The embl, genbank, swiss, gff and the other readers store the lines as the library does (a gff3 entry keeps its version line without the line end). Equals the installed program on 756 runs (every input format, damaged files, queries). A5-8c3h, B09 F257.")]
    public static EntretToolResult Entret(
        [Description("The text of the sequence file (any input format of seqret), or a sequence as letters.")] string sequence,
        [Description("The input format (-sformat); default: the program finds the format.")] string? sequenceFormat = null,
        [Description("Read one entry and stop (-firstonly).")] bool firstOnly = false,
        [Description("Sequence input qualifiers as typed after the sequence on the command line: -sformat1 F, -sid1 name, -sdbname1 db, -iquery1 query, -ioffset1 bytes, -sbegin1 -send1 -sreverse1 (the range does not change the text of an entry).")] string? qualifiers = null,
        CancellationToken cancellationToken = default,
        [Description("More files by name (file name -> text) that sequenceValue may name: list files, other sequence files, the data files of the databases of definitions.")] IReadOnlyDictionary<string, string>? files = null,
        [Description("The USA of the sequence input (-sequence): [format::]file[begin:end:r], file:id, file:field:value, file{a|b}, @list, list::list, asis::TEXT, db:id; the text of sequence is the file 'sequence.dat'.")] string? sequenceValue = null,
        [Description("The text of an EMBOSS configuration file (.embossrc) that defines databases for sequenceValue: DBNAME name [ type: N format: fasta method: direct directory: \"dir\" filename: \"*.fa\" ] (methods direct, emboss and emblcd), RES / RESOURCE, SERVER, INCLUDE (the included file is an entry of files), IF / IFDEF / ELSE / ENDIF, SET name value; the index files of the methods emboss (dbxfasta / dbxflat: NAME.ent, NAME.xid, NAME.pxid ...) and emblcd (dbifasta / dbiflat: division.lkp, entrynam.idx, acnum.trg ...) are entries of files by path, as latin-1 text.")] string? definitions = null)
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
            var o = EmbossEntretProgram.Entret(new EmbossEntretRequest(spec, all, q, firstOnly) { Definitions = string.IsNullOrWhiteSpace(definitions) ? null : [definitions] });
            if (o.Died)
                throw new ArgumentException(string.Join("; ", o.Messages.Where(m => m.Length > 0)), nameof(sequence));
            return new EntretToolResult(o.Text, o.Entries, [.. o.Messages]);
        }
        catch (Exception e) when (e is NotSupportedException or InvalidOperationException)
        {
            throw new ArgumentException(e.Message, nameof(sequence));
        }
    }
    [McpServerTool(Name = "descseq", Title = "Alignment — EMBOSS descseq (alter the name or description of a sequence)", ReadOnly = true)]
    [Description("EMBOSS 6.6.0 descseq: the name and / or the description of the FIRST sequence of a file (any input format of seqret) is replaced by name / description, or appended to the existing ones (append), and the sequence is returned in FASTA (the header is the name, the sequence version or accession, and the description). Equals the installed program on 1 000 runs (19 input formats, names and descriptions with blanks, tabs and punctuation, damaged files). A5-8c3g, B09 F250.")]
    public static DescseqToolResult Descseq(
        [Description("The text of the sequence file (any input format of seqret), or a sequence as letters.")] string sequence,
        [Description("The input format (-sformat); default: the program finds the format.")] string? sequenceFormat = null,
        [Description("The new name (-name); empty = the name is kept.")] string name = "",
        [Description("The new description (-description); empty = the description is kept.")] string description = "",
        [Description("Append the name and the description to the existing ones (-append).")] bool append = false,
        [Description("Sequence input qualifiers as typed after the sequence on the command line: -sbegin1 -send1 -sreverse1, -slower1 -supper1 -snucleotide1 -sprotein1 -scircular1, -sformat1 F, -sid1 name, -sdbname1 db, -squick1, -iquery1 query, -ioffset1 bytes.")] string? qualifiers = null,
        CancellationToken cancellationToken = default,
        [Description("More files by name (file name -> text) that sequence may name with a 'format::' prefix.")] IReadOnlyDictionary<string, string>? files = null,
        [Description("The USA of the sequence input (-sequence) as the program reads it, instead of the file the text is put in: [format::]file[begin:end:r] (file = the name of an entry of files; the text of sequence is the file 'sequence.dat'), file:id / file:field:query / file{a|b} (id, acc, sv, gi, des, key, org; wildcards * and ?), @list / list::list (a list file of files; # comments; every entry has its own range and format), asis::SEQUENCE; the range qualifiers -sbegin1 -send1 -sreverse1 of qualifiers are applied to the sequences and to their features.")] string? sequenceValue = null)
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
            var o = EmbossDescseqPrograms.Descseq(new EmbossDescseqRequest(spec, all, q, name ?? "", description ?? "", append));
            if (o.Died)
                throw new ArgumentException(string.Join("; ", o.Messages.Where(m => m.Length > 0)), nameof(sequence));
            return new DescseqToolResult(o.SequenceText, o.Sequences, [.. o.Messages]);
        }
        catch (Exception e) when (e is NotSupportedException or InvalidOperationException)
        {
            throw new ArgumentException(e.Message, nameof(sequence));
        }
    }
}
