#!/usr/bin/env python3
"""HISTORIC (B09 WP80): the tools this script generated now live in the EMBOSS MCP servers by area (Seqeron.Mcp.Emboss*; copyleft, docs/Validation/review-2026-09/LICENSES.md) and the tool files are maintained by hand; the paths and class names below are those of WP61-WP67.
B09 WP61 (F280): generates the MCP tools of the EMBOSS sequence programs of work package 61 on the Analysis server
(src/Seqeron/Mcp/Seqeron.Mcp.Analysis/Tools/AnalysisTools.EmbossSequencePrograms.cs: every tool builds the command line of the program from its typed parameters and calls EmbossProgramRunner.Run) and the spec
(scripts/docs/emboss_wp61_tool_spec.json) the document generators read (gen_mummer_est_tool_docs.py, register_mcp_tools.py).
Usage (repository root): python3 scripts/docs/gen_emboss_wp61_tools.py ; python3 scripts/docs/gen_mummer_est_tool_docs.py scripts/docs/emboss_wp61_tool_spec.json ;
python3 scripts/docs/register_mcp_tools.py scripts/docs/emboss_wp61_tool_spec.json ; python3 scripts/docs/sync_mcp_tool_docs.py ; python3 scripts/skills/gen-catalog.py"""
import gzip, json, os

ORACLE_DIR = 'tests/Seqeron/Seqeron.Genomics.Tests/TestData/Alignment'

# parameter kinds: flag (bool; default true -> -noNAME when false, default false -> -NAME when true), int / float / string / menu (-NAME value, only when given)
# tool: name, title, program, summary, desc, inputs [(acd, param, description)], opts [(param, kind, default, acd, description)], after, example
SEQ_IN = "The text of the sequence file (any input format of seqret), or a sequence as letters."
TOOLS = [
    dict(name='getorf', title='Genome — EMBOSS getorf (open reading frames)', program='getorf', after='entret',
         summary="Finds and extracts the open reading frames (ORFs) of nucleotide sequences, as getorf of EMBOSS 6.6.0.",
         desc="the regions between STOP codons or between START and STOP codons of the six frames (the circular sequence and the nucleotides around the codons too) are returned as protein or nucleotide sequences named name_1, name_2, ... with the header [start - end] (REVERSE SENSE) and the description of the sequence; the scan is Translator.FindOrfsGetorfProgram (the model of getorf_FindORFs over the EMBOSS genetic code files)",
         inputs=[('sequence', 'sequence', SEQ_IN)],
         opts=[('table', 'menu', '0', 'table', 'The genetic code (-table): 0 standard with AUG start only, 1 .. 6, 9 .. 16, 21 .. 23 (the EMBOSS EGC files).'),
               ('minSize', 'int', 30, 'minsize', 'The minimum size of an ORF in nucleotides (-minsize, default 30).'),
               ('maxSize', 'int', 1000000, 'maxsize', 'The maximum size of an ORF in nucleotides (-maxsize, default 1000000).'),
               ('find', 'menu', '0', 'find', 'The type of output (-find): 0 protein between STOP codons, 1 protein between START and STOP, 2 nucleotides between STOP codons, 3 nucleotides between START and STOP, 4 nucleotides around START codons, 5 around initial STOP codons, 6 around ending STOP codons.'),
               ('methionine', 'flag', True, 'methionine', 'A START codon at the beginning of a protein is an M (-methionine, default true).'),
               ('circular', 'flag', False, 'circular', 'The sequence is circular (-circular).'),
               ('reverse', 'flag', True, 'reverse', 'ORFs of the reverse complement are found as well (-reverse, default true).'),
               ('flanking', 'int', 100, 'flanking', 'The number of nucleotides on each side of the codon for find 4 .. 6 (-flanking, default 100).')],
         example=(">s1 desc here\nATGAAATTTGGGCCCTAGATGCCCTAAGGGATGAAATAA\n", {'minSize': 6})),
    dict(name='revseq', title='Genome — EMBOSS revseq (reverse and complement sequences)', program='revseq', after='getorf',
         summary="Reverses and complements nucleotide sequences, as revseq of EMBOSS 6.6.0.",
         desc="every sequence is reverse-complemented (or only reversed, or only complemented); with both and tag the description starts with 'Reversed:'",
         inputs=[('sequence', 'sequence', SEQ_IN)],
         opts=[('reverse', 'flag', True, 'reverse', 'Reverse the sequence (-reverse, default true; false = -noreverse).'),
               ('complement', 'flag', True, 'complement', 'Complement the sequence (-complement, default true).'),
               ('tag', 'flag', True, 'tag', "Put 'Reversed:' in front of the description of a reverse-complemented sequence (-tag, default true).")],
         example=(">s1 d\nAACCGGTT\n", {})),
    dict(name='transeq', title='Genome — EMBOSS transeq (translate nucleotide sequences)', program='transeq', after='revseq',
         summary="Translates nucleotide sequences to protein in the chosen frames, as transeq of EMBOSS 6.6.0.",
         desc="the frames 1 2 3 -1 -2 -3 (F = forward three, R = reverse three, 6 = all; a comma list), the genetic code, regions (translated joined together), trim of trailing X and *, clean (stops to X), the alternative definition of the reverse frames, and the first START codon as M for regions; the translation is the ajtrn.c model (EmbossTrn: ambiguity codes resolved as the EMBOSS genetic code tables do)",
         inputs=[('sequence', 'sequence', SEQ_IN)],
         opts=[('frame', 'string', None, 'frame', "The frames (-frame): 1 2 3 F -1 -2 -3 R 6, separated by commas; default 1."),
               ('table', 'menu', '0', 'table', 'The genetic code (-table): 0 .. 6, 9 .. 16, 21 .. 23.'),
               ('regions', 'string', None, 'regions', "The regions to translate (-regions): '4-57,78-94'; default the whole sequence."),
               ('trim', 'flag', False, 'trim', "Remove the 'X' and '*' at the end of the translation (-trim)."),
               ('clean', 'flag', False, 'clean', "Change the '*' of the STOP codons to 'X' (-clean)."),
               ('alternative', 'flag', False, 'alternative', "Frame -1 starts in the last codon (-alternative)."),
               ('methionine', 'flag', True, 'methionine', 'A START codon at the beginning of a region translation is an M (-methionine, default true).')],
         example=(">s1 desc here\nATGAAATTTGGGCCCTAGATGCCCTAAGGGATGAAATAA\n", {'frame': '1'})),
    dict(name='extractseq', title='Genome — EMBOSS extractseq (extract regions of a sequence)', program='extractseq', after='transeq', single=True,
         summary="Extracts regions of a sequence, as extractseq of EMBOSS 6.6.0.",
         desc="the regions (default the whole range of the sequence) are joined into one sequence, or with separate written one sequence per region named name_start_end",
         inputs=[('sequence', 'sequence', "The text of the sequence file (the first sequence is used), or a sequence as letters.")],
         opts=[('regions', 'string', None, 'regions', "The regions (-regions): '4-57,78-94', '1:45 67=99'; default the range of the sequence."),
               ('separate', 'flag', False, 'separate', 'Write every region as a separate sequence (-separate).')],
         example=(">s1 d\nACGTACGTACGTACGTACGT\n", {'regions': '3-5,8-10', 'separate': True})),
    dict(name='cutseq', title='Genome — EMBOSS cutseq (remove a section from a sequence)', program='cutseq', after='extractseq', single=True,
         summary="Removes a section from a sequence, as cutseq of EMBOSS 6.6.0.",
         desc="the residues from 'from' to 'to' are removed (default: from the first to the last residue, so the whole sequence goes; a value outside the limits is reset to the limit with the warning of the library)",
         inputs=[('sequence', 'sequence', "The text of the sequence file (the first sequence is used), or a sequence as letters.")],
         opts=[('from', 'int', None, 'from', 'The first residue to remove (-from).'), ('to', 'int', None, 'to', 'The last residue to remove (-to).')],
         example=(">s1 d\nACGTACGTACGTACGTACGT\n", {'from': 3, 'to': 5})),
    dict(name='degapseq', title='Genome — EMBOSS degapseq (remove gap characters)', program='degapseq', after='cutseq',
         summary="Removes the gap characters from sequences, as degapseq of EMBOSS 6.6.0.",
         desc="every character that is no letter and no star (the gaps '-', '.', '~' and blanks) is removed",
         inputs=[('sequence', 'sequence', SEQ_IN)], opts=[], example=(">s1 d\nAC--GT..AC\n", {})),
    dict(name='trimseq', title='Genome — EMBOSS trimseq (trim ambiguous ends)', program='trimseq', after='degapseq',
         summary="Trims the ambiguity characters and gaps from the ends of sequences, as trimseq of EMBOSS 6.6.0.",
         desc="a window moves in from each end; the sequence is cut where the share of bad characters (ambiguity codes, with strict also B and Z or N, X; gaps always) in the window falls below percent",
         inputs=[('sequence', 'sequence', SEQ_IN)],
         opts=[('window', 'int', None, 'window', 'The window size (-window, default 1; not longer than the first sequence).'),
               ('percent', 'float', None, 'percent', 'The percentage of bad characters that makes a window to be trimmed (-percent, default 100).'),
               ('strict', 'flag', False, 'strict', 'Trim the ambiguity codes too (-strict).'),
               ('star', 'flag', False, 'star', 'Protein: an asterisk is a good character (-star).'),
               ('left', 'flag', True, 'left', 'Trim the start (-left, default true).'),
               ('right', 'flag', True, 'right', 'Trim the end (-right, default true).')],
         example=(">s1 d\nNNNNACGTACGTNN\n", {})),
    dict(name='skipseq', title='Genome — EMBOSS skipseq (skip the first sequences)', program='skipseq', after='trimseq',
         summary="Reads and writes sequences, skipping the first n, as skipseq of EMBOSS 6.6.0.",
         desc="the first 'skip' sequences are left out",
         inputs=[('sequence', 'sequence', SEQ_IN)], opts=[('skip', 'int', None, 'skip', 'The number of sequences to skip (-skip, default 0).')],
         example=(">a\nAC\n>b\nGT\n", {'skip': 1})),
    dict(name='nthseq', title='Genome — EMBOSS nthseq (write one sequence)', program='nthseq', after='skipseq',
         summary="Writes the n-th sequence of a file, as nthseq of EMBOSS 6.6.0.",
         desc="the sequence with the number 'number'; with fewer sequences the library reports 'No such sequence - only N sequences were input.'",
         inputs=[('sequence', 'sequence', SEQ_IN)], opts=[('number', 'int', None, 'number', 'The number of the sequence (-number, default 1).')],
         example=(">a\nAC\n>b\nGT\n", {'number': 2})),
    dict(name='notseq', title='Genome — EMBOSS notseq (exclude sequences)', program='notseq', after='nthseq',
         summary="Writes the sequences that are not on a list, as notseq of EMBOSS 6.6.0.",
         desc="the sequences whose name or accession matches one of the patterns of exclude (wildcards * and ?, separated by commas or semicolons; '@file' reads the patterns from a file of the files dictionary) go to junkOutputFile (default: dropped), the others are the output",
         inputs=[('sequence', 'sequence', SEQ_IN)],
         opts=[('exclude', 'string', None, 'exclude', 'The patterns of the names or accessions to exclude (-exclude), or @file.')],
         example=(">a1 x\nAC\n>b2 y\nGT\n", {'exclude': 'a*'})),
    dict(name='pasteseq', title='Genome — EMBOSS pasteseq (insert a sequence into another)', program='pasteseq', after='notseq',
         summary="Inserts one sequence into another, as pasteseq of EMBOSS 6.6.0.",
         desc="the insert sequence is put into the sequence after position pos (0 = at the start; default the end)",
         inputs=[('asequence', 'sequence', 'The text of the sequence file (the first sequence is used), or a sequence as letters.'), ('bsequence', 'insertSequence', 'The sequence to insert (the same kind).')],
         opts=[('pos', 'int', None, 'pos', 'The position after which to insert (-pos, 0 .. length; default the end).')],
         example=(">a x\nAAAATTTT\n", {'insertSequence': ">b\nCCCC\n", 'pos': 4})),
dict(name='trimest', title='Genome — EMBOSS trimest (trim poly-A tails of ESTs)', program='trimest', after='pasteseq',
         summary="Removes the poly-A tails of EST sequences, as trimest of EMBOSS 6.6.0.",
         desc="the longer of the 3' poly-A tail and the 5' poly-T head (with mismatches allowed; N is ignored) is removed from every sequence, or written in lower case with tolower; a removed poly-T head makes the sequence reverse-complemented unless reverse is false; the description tells what was done",
         inputs=[('sequence', 'sequence', SEQ_IN)],
         opts=[('minLength', 'int', None, 'minlength', 'The minimum length of a tail (-minlength, default 4).'),
               ('mismatches', 'int', None, 'mismatches', 'The number of contiguous mismatches allowed in a tail (-mismatches, default 1).'),
               ('reverse', 'flag', True, 'reverse', 'Reverse-complement a sequence whose poly-T head was removed (-reverse, default true).'),
               ('toLower', 'flag', False, 'tolower', 'Change the tail to lower case instead of removing it (-tolower).'),
               ('fivePrime', 'flag', True, 'fiveprime', "Also look for a 5' poly-T head (-fiveprime, default true).")],
         example=(">e1 est\nTTTTTTTTACGTACGTNNACGTAAAAAAAAAAAA\n", {})),
    dict(name='newseq', title='Genome — EMBOSS newseq (type in a short sequence)', program='newseq', after='trimest', kind='noinput',
         summary="Makes a sequence from a name, a description, a type and the residues, as newseq of EMBOSS 6.6.0.",
         desc="the sequence is written in the chosen format with the gaps and the ambiguity codes of its type (nucleic: gapnucleotide; protein: gapprotein) as the output sequence check leaves them",
         inputs=[], opts=[('name', 'string', None, 'name', 'The name of the sequence (-name).'), ('description', 'string', None, 'description', 'The description of the sequence (-description).'),
                          ('type', 'menu', 'N', 'type', 'N nucleic or P protein (-type).'), ('residues', 'string', None, 'sequence', 'The residues (-sequence).')],
         example=(None, {'name': 'nn', 'description': 'd d', 'type': 'N', 'residues': 'ACGTRXYZ'})),
    dict(name='yank', title='Genome — EMBOSS yank (the USA of a sequence)', program='yank', after='newseq', single=True, kind='outfile',
         summary="Writes the USA (the universal sequence address) of a sequence, as yank of EMBOSS 6.6.0.",
         desc="the text is the USA of the first sequence as the library makes it: format::file:name with the range [begin:end] or [end:begin:r]",
         inputs=[('sequence', 'sequence', "The text of the sequence file (the first sequence is used), or a sequence as letters.")], opts=[],
         example=(">s1 d\nACGT\n", {})),
    dict(name='seqcount', title='Genome — EMBOSS seqcount (count the sequences)', program='seqcount', after='yank', kind='outfile',
         summary="Counts the sequences of a file, as seqcount of EMBOSS 6.6.0.",
         desc="the number of sequences the file holds (with the range qualifiers, the sequences read)",
         inputs=[('sequence', 'sequence', SEQ_IN)], opts=[], example=(">a\nAC\n>b\nGT\n", {})),
    dict(name='splitter', title='Genome — EMBOSS splitter (split sequences into pieces)', program='splitter', after='seqcount',
         summary="Splits sequences into smaller sequences, as splitter of EMBOSS 6.6.0.",
         desc="the pieces of size residues (overlapping by overlap; addOverlap puts the overlap on top of the size) are named name_begin-end and keep the description; the features of the pieces (-feature) are not written",
         inputs=[('sequence', 'sequence', SEQ_IN)],
         opts=[('size', 'int', None, 'size', 'The size of the pieces (-size, default 10000).'), ('overlap', 'int', None, 'overlap', 'The overlap of the pieces (-overlap, default 0; at most size - 1).'),
               ('addOverlap', 'flag', False, 'addoverlap', 'The overlap is added to the pieces (-addoverlap).')],
         example=(">s1 d\nACGTACGTACGTACGTACGT\n", {'size': 8})),
    dict(name='maskseq', title='Genome — EMBOSS maskseq (mask regions of a sequence)', program='maskseq', after='splitter', single=True,
         summary="Masks regions of a sequence with a character or lower case, as maskseq of EMBOSS 6.6.0.",
         desc="the regions (-regions '4-57,78-94') are replaced by the mask character (N, X for protein) or written in lower case with toLower",
         inputs=[('sequence', 'sequence', "The text of the sequence file (the first sequence is used), or a sequence as letters.")],
         opts=[('regions', 'string', None, 'regions', 'The regions to mask (-regions).'), ('toLower', 'flag', False, 'tolower', 'Write the regions in lower case instead of masking them (-tolower).'),
               ('maskChar', 'string', None, 'maskchar', 'The mask character (-maskchar, one character).')],
         example=(">s1 d\nACGTACGTACGTACGTACGT\n", {'regions': '3-5'})),
    dict(name='maskambignuc', title='Genome — EMBOSS maskambignuc (mask ambiguous nucleotides)', program='maskambignuc', after='maskseq',
         summary="Masks all ambiguity nucleotide codes with N, as maskambignuc of EMBOSS 6.6.0.",
         desc="the input type gapnucleotidesimple turns every nucleotide ambiguity code into N (gaps are kept)",
         inputs=[('sequence', 'sequence', SEQ_IN)], opts=[], example=(">a\nACGTRYKMNBDHV-\n", {})),
    dict(name='maskambigprot', title='Genome — EMBOSS maskambigprot (mask ambiguous amino acids)', program='maskambigprot', after='maskambignuc',
         summary="Masks all ambiguity amino acid codes with X, as maskambigprot of EMBOSS 6.6.0.",
         desc="the input type gapproteinsimple turns B, Z, U, J, O and * into X (gaps are kept)",
         inputs=[('sequence', 'sequence', SEQ_IN)], opts=[], example=(">a\nACDEFBZXUJO*-\n", {})),

dict(name='trimspace', title='Genome — EMBOSS trimspace (remove extra whitespace)', program='trimspace', after='maskambigprot', kind='text',
         summary="Cleans up whitespace in a text file, as trimspace of EMBOSS 6.6.0.",
         desc="every line loses its leading and trailing white space and every run of blanks becomes one (tabs and carriage returns are blanks)",
         inputs=[], opts=[], example=("  a \t b   c  \n\n x\ty\n", {})),
    dict(name='nospace', title='Genome — EMBOSS nospace (remove whitespace)', program='nospace', after='trimspace', kind='text',
         summary="Removes whitespace from a text file, as nospace of EMBOSS 6.6.0.",
         desc="menu all removes all white space, end the trailing white space, excess the leading and trailing white space and the extra blanks (a unique abbreviation of the menu value is accepted)",
         inputs=[], opts=[('menu', 'menu', 'all', 'menu', 'all, end or excess (-menu, default all).')], example=("  a \t b   c  \n\n x\ty\n", {})),
    dict(name='notab', title='Genome — EMBOSS notab (replace tabs by blanks)', program='notab', after='nospace', kind='text',
         summary="Replaces the tabs of a text file with blanks, as notab of EMBOSS 6.6.0.",
         desc="every tab becomes the blanks up to the next multiple of eight columns",
         inputs=[], opts=[], example=("a\tb\n12345678\tc\n", {})),
    dict(name='noreturn', title='Genome — EMBOSS noreturn (change the line ends)', program='noreturn', after='notab', kind='text',
         summary="Changes the line ends of a text file, as noreturn of EMBOSS 6.6.0.",
         desc="the line ends become those of the system (unix LF, pc CRLF, mac CR); a carriage return inside a line is a line end unless the output is mac",
         inputs=[], opts=[('system', 'menu', 'unix', 'system', 'unix, pc or mac (-system, default unix).')], example=("a\rb\r\nc\n", {})),

]

CS_TYPES = {'flag': 'bool', 'int': 'int?', 'float': 'float?', 'string': 'string?', 'menu': 'string?'}


def cs_str(s):
    return '"' + s.replace('\\', '\\\\').replace('"', '\\"').replace('\n', '\\n').replace('\r', '\\r').replace('\t', '\\t') + '"'


def oracle_rows(prog):
    p = os.path.join(ORACLE_DIR, 'emboss_prog_%s_oracle.jsonl.gz' % prog)
    return sum(1 for _ in gzip.open(p, 'rt')) if os.path.exists(p) else 0


def method_name(t):
    return t['name'][0].upper() + t['name'][1:]


def gen_cs():
    out = ['using System.ComponentModel;', 'using ModelContextProtocol.Server;', 'using global::Seqeron.Genomics.Analysis;', '', 'namespace Seqeron.Mcp.Analysis.Tools;', '',
           '// <auto-generated by scripts/docs/gen_emboss_wp61_tools.py: do not edit> EMBOSS 6.6.0 sequence programs (B09 WP61, F280): every tool builds the command line of the program and runs it with EmbossProgramRunner.', 'public partial class AnalysisTools', '{']
    for t in TOOLS:
        kind = t.get('kind', 'seqout')
        n = oracle_rows(t['program'])
        desc = "EMBOSS 6.6.0 %s: %s. Equals the installed program on %d runs (every input format of seqret, ranges, the output formats, the messages of the library; byte for byte). B09 WP61, F280." % (t['program'], t['desc'], n)
        t['_desc'] = desc
        t['_n'] = n
        out.append('    [McpServerTool(Name = %s, Title = %s, ReadOnly = true)]' % (cs_str(t['name']), cs_str(t['title'])))
        out.append('    [Description(%s)]' % cs_str(desc))
        params = []
        first = None
        if kind == 'text':
            first = 'text'
            params.append('        [Description("The text of the file to process.")] string text')
        for acd, p, d in t['inputs']:
            first = first or p
            params.append('        [Description(%s)] string %s' % (cs_str(d), p) + ('' if p == t['inputs'][0][1] else ' = ""'))
        for p, k, dflt, acd, d in t['opts']:
            pn = p if p not in ('from', 'to', 'type') else '@' + p
            if kind == 'noinput':
                first = first or p
                params.append('        [Description(%s)] string %s' % (cs_str(d), pn) + ('' if pn == first else (' = "N"' if k == 'menu' else ' = ""')))
                continue
            dv = {'flag': 'true' if dflt else 'false', 'int': 'null', 'float': 'null', 'string': 'null', 'menu': 'null'}[k]
            params.append('        [Description(%s)] %s %s = %s' % (cs_str(d), CS_TYPES[k], pn, dv))
        if kind not in ('noinput', 'text'):
            params.append('        [Description("The input format (-sformat1); default: the program finds the format.")] string? sequenceFormat = null')
            params.append('        [Description(%s)] string? qualifiers = null' % cs_str("Sequence input qualifiers as typed after the sequence on the command line: -sbegin1 -send1 -sreverse1, -slower1 -supper1 -snucleotide1 -sprotein1 -scircular1, -sid1 name, -sdbname1 db, -squick1, -iquery1 query, -ioffset1 bytes."))
        if kind not in ('outfile', 'text'):
            params.append('        [Description("The sequence output format (-osformat): fasta (default), embl, genbank, gcg, nbrf, ig, raw, msf, clustal, phylip, nexus, ... every format of seqret.")] string outputFormat = "fasta"')
        params.append('        CancellationToken cancellationToken = default')
        if kind not in ('noinput', 'text'):
            params.append('        [Description("More files by name (file name -> text) that the program may read (a file named by a USA, @file).")] IReadOnlyDictionary<string, string>? files = null')
            params.append('        [Description("The USA of the sequence input (-sequence) as the program reads it, instead of the text: [format::]file[begin:end:r], file:id, @list, asis::SEQUENCE.")] string? sequenceValue = null')
        out.append('    public static EmbossProgramToolResult %s(' % method_name(t))
        out.append(',\n'.join(params) + ')')
        out.append('    {')
        if kind != 'noinput':
            out.append('        if (string.IsNullOrWhiteSpace(%s))' % first)
            out.append('            throw new ArgumentException("The sequence cannot be empty", nameof(%s));' % first)
        out.append('        cancellationToken.ThrowIfCancellationRequested();')
        out.append('        var all = new Dictionary<string, string>(%s);' % ('files ?? new Dictionary<string, string>()' if kind not in ('noinput', 'text') else ''))
        if kind == 'text':
            out.append('        all["in.txt"] = text;')
            out.append('        var args = new List<string> { "-infile", "in.txt" };')
        if kind != 'text':
            out.append('        var args = new List<string>();')
        for i, (acd, p, d) in enumerate(t['inputs']):
            f = '%s.dat' % acd
            out.append('        all[%s] = %s;' % (cs_str(f), p))
            out.append('        args.Add("-%s");' % acd)
            out.append('        args.Add(SequenceSpec(%s, %s, %s));' % ('sequenceValue' if i == 0 else 'null', p, cs_str(f)))
        if kind not in ('noinput', 'text'):
            out.append('        if (!string.IsNullOrWhiteSpace(sequenceFormat))')
            out.append('            args.AddRange(["-sformat1", sequenceFormat]);')
            out.append('        AddQualifierTokens(args, qualifiers);')
        for p, k, dflt, acd, d in t['opts']:
            pn = p if p not in ('from', 'to', 'type') else '@' + p
            if kind == 'noinput':
                out.append('        args.AddRange(["-%s", %s]);' % (acd, pn))
            elif kind == 'text' and k == 'menu':
                out.append('        if (!string.IsNullOrEmpty(%s))\n            args.AddRange(["-%s", %s]);' % (pn, acd, pn))
            elif k == 'flag':
                if dflt:
                    out.append('        if (!%s)\n            args.Add("-no%s");' % (pn, acd))
                else:
                    out.append('        if (%s)\n            args.Add("-%s");' % (pn, acd))
            elif k in ('int', 'float'):
                out.append('        if (%s.HasValue)\n            args.AddRange(["-%s", %s.Value.ToString(System.Globalization.CultureInfo.InvariantCulture)]);' % (pn, acd, pn))
            else:
                out.append('        if (!string.IsNullOrEmpty(%s))\n            args.AddRange(["-%s", %s]);' % (pn, acd, pn))
        if t['program'] == 'notseq':
            out.append('        args.AddRange(["-junkoutseq", "junk.out"]);')
        param = first if kind != 'noinput' else 'name'
        if kind in ('outfile', 'text'):
            out.append('        return RunEmbossOutfileProgram(%s, args, all, nameof(%s));' % (cs_str(t['program']), param))
        else:
            out.append('        return RunEmbossSequenceProgram(%s, args, all, outputFormat, nameof(%s));' % (cs_str(t['program']), param))
        out.append('    }')
        out.append('')
    out.append('}')
    open('src/Seqeron/Mcp/Seqeron.Mcp.Analysis/Tools/AnalysisTools.EmbossSequencePrograms.cs', 'w').write('\n'.join(out) + '\n')


def example_output(t, seq, opts):
    """the output of the installed program for the example of the tool"""
    import subprocess, tempfile
    kind = t.get('kind', 'seqout')
    with tempfile.TemporaryDirectory() as d:
        args = []
        if kind == 'text':
            open(os.path.join(d, 'in.txt'), 'w', newline='').write(seq)
            args += ['-infile', 'in.txt']
        for i, (acd, p, desc) in enumerate(t['inputs']):
            f = '%s.dat' % acd
            open(os.path.join(d, f), 'w').write(seq if i == 0 else opts.get(p, ''))
            args += ['-' + acd, f]
        for p, k, dflt, acd, desc in t['opts']:
            if p not in opts and k != 'flag':
                continue
            v = opts.get(p, dflt)
            if k == 'flag':
                if bool(v) != bool(dflt): args.append(('-no' if dflt else '-') + acd)
            else:
                args += ['-' + acd, str(v)]
        if t['program'] == 'notseq': args += ['-junkoutseq', 'junk.out']
        args += ['-outfile', 'out.txt'] if kind in ('outfile', 'text') else ['-outseq', 'stdout']
        r = subprocess.run([t['program'], '-auto'] + args, cwd=d, capture_output=True, text=True)
        text = open(os.path.join(d, 'out.txt'), newline='').read() if kind in ('outfile', 'text') and os.path.exists(os.path.join(d, 'out.txt')) else r.stdout
        return {'text': text, 'files': {}, 'messages': [m for m in r.stderr.split('\n') if m]}


def gen_spec():
    spec = []
    for t in TOOLS:
        seq, opts = t['example']
        ex_in = {t['inputs'][0][1]: seq} if t['inputs'] else ({'text': seq} if t.get('kind') == 'text' else {})
        ex_in.update(opts)
        ex_in = {k: v for k, v in ex_in.items() if v is not None}
        spec.append({
            'name': t['name'], 'summary': "EMBOSS 6.6.0 %s: %s" % (t['program'], t['summary'].rstrip('.') + '.'),
            'methodId': 'EmbossProgramRunner.Run', 'sourceFile': 'EmbossSequencePrograms.cs' if t['program'] != 'getorf' else 'EmbossGetorfProgram.cs',
            'sourcePath': 'Seqeron.Genomics.Analysis/' + ('EmbossSequencePrograms.cs' if t['program'] != 'getorf' else 'EmbossGetorfProgram.cs'),
            'algorithmDoc': 'Alignment/Pairwise_Output_Formats_And_EMBOSS_Programs.md', 'oracleScript': 'generate_emboss_prog_getorf_oracle.py' if t['program'] == 'getorf' else 'generate_emboss_prog_seq_oracle.py',
            'short': "EMBOSS %s: %s; equals %s 6.6.0 on %d runs." % (t['program'], t['summary'].split(', as ')[0].rstrip('.').lower().replace('finds and ', 'finds and '), t['program'], t['_n']),
            'after': t['after'], 'file': 'AnalysisTools.EmbossSequencePrograms.cs', 'description': t['_desc'],
            'output': {'text': ['string', 'The standard output of the program: the sequences written to the standard output.'], 'files': ['object', 'The files the program wrote, by name.'], 'messages': ['string[]', 'What the library printed on the error stream (warnings and errors).']},
            'errors': ['The sequence cannot be empty', 'The program died: the messages of the library (a value it refuses, an input no reader reads, a sequence of the wrong kind)', 'An input on which the reference program never stops or dies by a fault is refused'],
            'example': {'title': 'A run on a small input', 'prompt': 'Run %s on this sequence.' % t['program'], 'input': ex_in, 'output': example_output(t, seq, opts), 'note': 'The output is the text the installed program writes.'},
        })
    json.dump(spec, open('scripts/docs/emboss_wp61_tool_spec.json', 'w'), indent=1)


def gen_tests():
    out = ['using NUnit.Framework;', 'using Seqeron.Mcp.Analysis.Tools;', '', 'namespace Seqeron.Mcp.Analysis.Tests;', '',
           '// <auto-generated by scripts/docs/gen_emboss_wp61_tools.py: do not edit> B09 WP61 (F280): the MCP tools of the EMBOSS sequence programs; the expected text is what the installed EMBOSS 6.6.0 program writes for the example of the tool document.',
           '[TestFixture]', 'public class EmbossSequenceProgramToolsTests', '{']
    for t in TOOLS:
        kind_t = t.get('kind', 'seqout')
        seq, opts = t['example']
        ex = example_output(t, seq, opts)
        args = []
        if t['inputs'] or kind_t == 'text':
            args.append(cs_str(seq))
        for p, k, dflt, acd, d in t['opts']:
            if p not in opts:
                continue
            pn = p if p not in ('from', 'to', 'type') else '@' + p
            v = opts[p]
            lit = ('true' if v else 'false') if isinstance(v, bool) else cs_str(v) if isinstance(v, str) else str(v)
            args.append('%s: %s' % (pn, lit))
        if t['program'] == 'pasteseq':
            args = [a for a in args if not a.startswith('insertSequence')] 
            args.insert(1, 'insertSequence: ' + cs_str(opts['insertSequence']))
        out.append('    [Test]')
        out.append('    public void %s_WritesWhatTheInstalledProgramWrites()' % method_name(t))
        out.append('    {')
        out.append('        var r = AnalysisTools.%s(%s);' % (method_name(t), ', '.join(args)))
        out.append('        Assert.That(r.Text, Is.EqualTo(%s));' % cs_str(ex['text']))
        out.append('    }')
        out.append('')
    out.append('    [Test]')
    out.append('    public void Tools_RefuseAnEmptySequenceAndABadValueWithTheMessagesOfTheLibrary()')
    out.append('    {')
    out.append('        Assert.Throws<ArgumentException>(() => AnalysisTools.Revseq(""));')
    out.append('        var e = Assert.Throws<ArgumentException>(() => AnalysisTools.Getorf(">s1\\nACGTACGTACGT\\n", table: "7"));')
    out.append('        Assert.That(e!.Message, Does.Contain("Bad menu option \'7\'"));')
    out.append('        Assert.That(AnalysisTools.Nthseq(">a\\nAC\\n", number: 2).Messages, Does.Contain("Error: No such sequence - only 1 sequences were input."));')
    out.append('    }')
    out.append('}')
    open('tests/Seqeron/Seqeron.Mcp.EmbossSequenceFormats.Tests/EmbossSequenceProgramToolsTests.cs', 'w').write('\n'.join(out) + '\n')


if __name__ == '__main__':
    gen_cs()
    gen_spec()
    gen_tests()
