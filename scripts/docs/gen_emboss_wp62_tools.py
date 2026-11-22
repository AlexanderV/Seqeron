#!/usr/bin/env python3
"""B09 WP62 (F281): generates the MCP tools of the EMBOSS programs of work package 62 on the Analysis server (shuffleseq msbar makenucseq makeprotseq nthseqset listor infoseq wordcount seqxref; the programs of
later parts of the package are added to TOOLS): src/Seqeron/Mcp/Seqeron.Mcp.Analysis/Tools/AnalysisTools.EmbossSequencePrograms2.cs (every tool builds the command line of the program and calls EmbossProgramRunner.Run), the spec
(scripts/docs/emboss_wp62_tool_spec.json) the document generators read, and the tests of the tools (the expected text is what the installed program writes for the example, run under the fixed clock of the example).
Usage (repository root): python3 scripts/docs/gen_emboss_wp62_tools.py ; python3 scripts/docs/gen_mummer_est_tool_docs.py scripts/docs/emboss_wp62_tool_spec.json ;
python3 scripts/docs/register_mcp_tools.py scripts/docs/emboss_wp62_tool_spec.json ; python3 scripts/docs/sync_mcp_tool_docs.py ; python3 scripts/skills/gen-catalog.py"""
import gzip, json, os, sys, subprocess, tempfile
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from gen_emboss_wp61_tools import cs_str, CS_TYPES, method_name  # noqa: E402

ORACLE_DIR = 'tests/Seqeron/Seqeron.Genomics.Tests/TestData/Alignment'
SEQ_IN = "The text of the sequence file (any input format of seqret), or a sequence as letters."
USEC = ('microseconds', 'int', None, None, "The microseconds of the clock when the random number generator is seeded (the seed is microseconds % 9999 + 1; the program takes them from the clock of the machine): the same value gives the same sequences. Default: the clock.")
# tool: name, title, program, kind (seqout | outfile | noinput | outfile2), inputs [(acd, param, description)], opts [(param, kind, default, acd, description[, file])]
TOOLS = [
    dict(name='shuffleseq', title='Genome — EMBOSS shuffleseq (shuffle sequences keeping the composition)', program='shuffleseq', after='seqcount', random=True,
         summary="Shuffles sequences keeping their composition, as shuffleseq of EMBOSS 6.6.0.",
         desc="every sequence is written shuffle times with its residues in random order (a random number for every position, then the Shell sort of the positions by it: ajStrRandom of the EMBOSS library); the random numbers are ajRandomDouble (Wichmann-Hill with the 101-term lagged generator of ajmath.c), seeded from the microseconds of the clock",
         inputs=[('sequence', 'sequence', SEQ_IN)], opts=[('shuffle', 'int', None, 'shuffle', 'The number of shuffles of every sequence (-shuffle, default 1).')],
         example=(">s1 d\nACGTACGTACGTAAAACCCGGGTTT\n", {'shuffle': 2, 'microseconds': 7})),
    dict(name='msbar', title='Genome — EMBOSS msbar (mutate sequences)', program='msbar', after='shuffleseq', random=True,
         summary="Mutates sequences with point, codon and block mutations, as msbar of EMBOSS 6.6.0.",
         desc="count times the block, codon (nucleic sequences) and point mutations of the lists are done at random places: insertions, deletions, changes, duplications and moves (1 = a random one of them, 0 = none; a comma list of at most 4); the result is tried again (10 times) when it equals a sequence of otherSequence",
         inputs=[('sequence', 'sequence', SEQ_IN)],
         opts=[('count', 'int', None, 'count', 'The number of times the mutation operations are done (-count, default 1).'),
               ('point', 'string', None, 'point', "The point mutations (-point): 0 none, 1 any, 2 insertions, 3 deletions, 4 changes, 5 duplications, 6 moves; a comma list of at most 4 (default 0)."),
               ('block', 'string', None, 'block', 'The block mutations (-block), as point (default 0).'),
               ('codon', 'string', None, 'codon', 'The codon mutations (-codon), as point (default 0; nucleic sequences only).'),
               ('inframe', 'flag', False, 'inframe', 'Do the codon and block operations in frame (-inframe).'),
               ('minimum', 'int', None, 'minimum', 'The minimum size of a block mutation (-minimum, default 1).'),
               ('maximum', 'int', None, 'maximum', 'The maximum size of a block mutation (-maximum, default 10).'),
               ('otherSequence', 'string', None, 'othersequence', 'Sequences the result must not equal (-othersequence: the text of a sequence file, or a USA such as asis:ACGT).', 'other.dat')],
         example=(">s1 d\nACGTACGTACGTAAAACCCGGGTTT\n", {'point': '2,3', 'count': 3, 'microseconds': 7})),
    dict(name='makenucseq', title='Genome — EMBOSS makenucseq (random nucleotide sequences)', program='makenucseq', after='msbar', kind='noinput', random=True,
         summary="Creates random nucleotide sequences, as makenucseq of EMBOSS 6.6.0.",
         desc="amount sequences of length lower case nucleotides, each residue picked at random from ACGT (or from the codons of the codon usage file with the frequencies of the file, the end cut to the length), with the insert put at start; the names are EMBOSS_001, EMBOSS_002 ...",
         inputs=[],
         opts=[('amount', 'int', None, 'amount', 'The number of sequences (-amount, default 100).'), ('length', 'int', None, 'length', 'The length of every sequence (-length, default 100).'),
               ('insert', 'string', None, 'insert', 'A string put into every sequence (-insert).'), ('start', 'int', None, 'start', 'The position of the insert (-start, default 1).'),
               ('codonFile', 'string', None, 'codonfile', 'The text of a codon usage file in the EMBOSS (cusp) format (-codonfile).', 'codon.cut')],
         example=(None, {'amount': 2, 'length': 20, 'microseconds': 7})),
    dict(name='makeprotseq', title='Genome — EMBOSS makeprotseq (random protein sequences)', program='makeprotseq', after='makenucseq', kind='noinput', random=True,
         summary="Creates random protein sequences, as makeprotseq of EMBOSS 6.6.0.",
         desc="amount sequences of length lower case residues, each picked at random from the 20 amino acids (or with the composition of the first sequence of a pepstats file), with the insert put at start; the names are EMBOSS_001, EMBOSS_002 ...",
         inputs=[],
         opts=[('amount', 'int', None, 'amount', 'The number of sequences (-amount, default 100).'), ('length', 'int', None, 'length', 'The length of every sequence (-length, default 100).'),
               ('insert', 'string', None, 'insert', 'A string put into every sequence (-insert).'), ('start', 'int', None, 'start', 'The position of the insert (-start, default 1).'),
               ('pepstatsFile', 'string', None, 'pepstatsfile', 'The text of a pepstats output file (-pepstatsfile).', 'pepstats.txt')],
         example=(None, {'amount': 2, 'length': 12, 'microseconds': 7})),
    dict(name='nthseqset', title='Genome — EMBOSS nthseqset (write one set of sequences)', program='nthseqset', after='msbar', single=False,
         summary="Writes the n-th set of sequences of an input of several sets, as nthseqset of EMBOSS 6.6.0.",
         desc="the sets of the input (a file of several alignments: phylip, nexus ... or one set of all the sequences of the file) are counted and the number-th is written; with fewer sets the library reports 'No such sequence set - only N sequence sets were input.'",
         inputs=[('sequence', 'sequence', SEQ_IN)], opts=[('number', 'int', None, 'number', 'The number of the set (-number, default 1).')],
         example=(">a\nAC\n>b\nGT\n", {'number': 1})),
    dict(name='listor', title='Genome — EMBOSS listor (logical OR of two sets of sequences)', program='listor', after='nthseqset', kind='outfile2',
         summary="Writes a list file of the logical OR / AND / XOR / NOT of two sets of sequences, as listor of EMBOSS 6.6.0.",
         desc="the USAs of the sequences of the first set and the second set that are (OR: all, AND: in both, XOR: in one only, NOT: in the first and not in the second); sequences are equal when their residues are equal ignoring the case",
         inputs=[('firstsequences', 'firstSequences', SEQ_IN), ('secondsequences', 'secondSequences', 'The text of the second sequence file (any input format of seqret), or a sequence as letters.')],
         opts=[('operator', 'menu', 'OR', 'operator', 'OR, AND, XOR or NOT (-operator, default OR; a unique abbreviation is accepted).')],
         example=(">a\nAC\n>b\nGT\n", {'secondSequences': ">c\nac\n>d\nTT\n", 'operator': 'AND'})),
    dict(name='infoseq', title='Genome — EMBOSS infoseq (basic information about sequences)', program='infoseq', after='listor', kind='outfile',
         summary="Displays basic information about sequences, as infoseq of EMBOSS 6.6.0.",
         desc="for every sequence the USA, database, name, accession, GI, version, type, length, percent GC (nucleic), organism and description, in columns (or delimited, or as an HTML table), with a heading line",
         inputs=[('sequence', 'sequence', SEQ_IN)],
         opts=[('columns', 'flag', True, 'columns', 'Print the information in aligned columns (-columns, default true).'), ('delimiter', 'string', None, 'delimiter', 'The delimiter of the records without columns (-delimiter, default |).'),
               ('html', 'flag', False, 'html', 'Format the output as an HTML table (-html).'), ('only', 'flag', False, 'only', 'Display only the columns that are asked for (-only).'),
               ('heading', 'flag', True, 'heading', 'Display the column headings (-heading, default true).'), ('usa', 'flag', True, 'usa', 'The USA column (-usa).'), ('database', 'flag', True, 'database', 'The database column (-database).'),
               ('name', 'flag', True, 'name', 'The name column (-name).'), ('accession', 'flag', True, 'accession', 'The accession column (-accession).'), ('gi', 'flag', False, 'gi', 'The GI column (-gi).'),
               ('seqversion', 'flag', False, 'seqversion', 'The version column (-seqversion).'), ('type', 'flag', True, 'type', 'The type column (-type).'), ('length', 'flag', True, 'length', 'The length column (-length).'),
               ('pgc', 'flag', True, 'pgc', 'The percent GC column (-pgc).'), ('organism', 'flag', True, 'organism', 'The organism column (-organism).'), ('description', 'flag', True, 'description', 'The description column (-description).')],
         example=(">a desc here\nACGTACGTAC\n>b\nMKVLAAGIV\n", {})),
    dict(name='wordcount', title='Genome — EMBOSS wordcount (count the words of sequences)', program='wordcount', after='infoseq', kind='outfile',
         summary="Counts the words of a given size in sequences, as wordcount of EMBOSS 6.6.0.",
         desc="the words of wordSize residues (default 4 for nucleic, 2 for protein sequences; words with the ambiguity code N / X are skipped) of all the sequences with their counts, the most frequent first (equal counts in the order of the hash table of the library), at least minCount",
         inputs=[('sequence', 'sequence', SEQ_IN)], opts=[('wordSize', 'int', None, 'wordsize', 'The size of a word (-wordsize).'), ('minCount', 'int', None, 'mincount', 'The minimum count to report (-mincount, default 1).')],
         example=(">a\nACGTACGTACGTTTTT\n", {'wordSize': 3})),
    dict(name='seqxref', title='Genome — EMBOSS seqxref (database cross-references of sequence entries)', program='seqxref', after='wordcount', kind='outfile',
         summary="Reports the database cross-references of sequence entries, as seqxref of EMBOSS 6.6.0.",
         desc="for every entry '#USA: count' and the cross-references db:id (the DR lines of EMBL and Swiss-Prot entries ...) sorted by database name and identifier, equal ones once",
         inputs=[('sequence', 'sequence', SEQ_IN)], opts=[],
         example=("ID   TEST_HUMAN    Reviewed;    10 AA.\nAC   P12345;\nDR   PDB; 1ABC; X-ray; 1.9 A; A=1-10.\nDR   EMBL; X1234; -.\nSQ   SEQUENCE   10 AA;  1000 MW;  ABCD CRC64;\n     MKVLAAGIVG\n//\n", {})),
    dict(name='showorf', title='Genome — EMBOSS showorf (nucleotide sequence and its translations)', program='showorf', after='seqxref', kind='outfile', single=True,
         summary="Displays a nucleotide sequence and its translations in pretty format, as showorf of EMBOSS 6.6.0.",
         desc="the sequence in blocks of width bases with a ruler and the numbers of the bases, and under it the translations of the frames (F1 F2 F3 forward, R1 R2 R3 reverse: the amino acid under the first base of every codon, the numbers of the first and last residue of the block); the genetic code is the EMBOSS table of the code number",
         inputs=[('sequence', 'sequence', "The text of the sequence file (the first sequence is used), or a sequence as letters.")],
         opts=[('frames', 'string', None, 'frames', "The frames to translate (-frames): a comma list of 0 none, 1 F1, 2 F2, 3 F3, 4 R1, 5 R2, 6 R3 (default 1,2,3,4,5,6)."),
               ('table', 'string', None, 'table', 'The genetic code (-table): 0 .. 6, 9 .. 16, 21 .. 23.'), ('width', 'int', None, 'width', 'The width of the screen (-width, default 50; at least 10).'),
               ('ruler', 'flag', True, 'ruler', 'Add a ruler (-ruler, default true).'), ('plabel', 'flag', True, 'plabel', 'Number the translations (-plabel, default true).'), ('nlabel', 'flag', True, 'nlabel', 'Number the sequence (-nlabel, default true).')],
         example=(">s1 d\nATGAAATTTGGGCCCTAGATGCCCAAATTT\n", {'frames': '1,4', 'width': 30})),
    dict(name='prettyseq', title='Genome — EMBOSS prettyseq (nucleotide sequence and its translation)', program='prettyseq', after='showorf', kind='outfile', single=True,
         summary="Writes a nucleotide sequence and its translation, as prettyseq of EMBOSS 6.6.0.",
         desc="the sequence in blocks of width bases with a ruler and numbers; the bases of the ranges are in lower case and the amino acid is written under the first base of each of their codons (the rest of the sequence in upper case with blanks under it)",
         inputs=[('sequence', 'sequence', "The text of the sequence file (the first sequence is used), or a sequence as letters.")],
         opts=[('range', 'string', None, 'range', "The ranges to translate (-range '4-57,78-94'; default the whole sequence)."), ('table', 'string', None, 'table', 'The genetic code (-table): 0 .. 6, 9 .. 16, 21 .. 23.'),
               ('ruler', 'flag', True, 'ruler', 'Add a ruler (-ruler, default true).'), ('plabel', 'flag', True, 'plabel', 'Number the translations (-plabel, default true).'), ('nlabel', 'flag', True, 'nlabel', 'Number the sequence (-nlabel, default true).'),
               ('width', 'int', None, 'width', 'The width of the screen (-width, default 60; at least 10).')],
         example=(">s1 d\nATGAAATTTGGGCCCTAGATGCCCAAATTT\n", {'range': '4-15'})),
]


def oracle_rows(prog):
    p = os.path.join(ORACLE_DIR, 'emboss_prog_%s_oracle.jsonl.gz' % prog)
    return sum(1 for _ in gzip.open(p, 'rt')) if os.path.exists(p) else 0


def pn_of(p):
    return p if p not in ('from', 'to', 'type', 'operator') else '@' + p


def gen_cs():
    out = ['using System.ComponentModel;', 'using ModelContextProtocol.Server;', 'using global::Seqeron.Genomics.Analysis;', '', 'namespace Seqeron.Mcp.Analysis.Tools;', '',
           '// <auto-generated by scripts/docs/gen_emboss_wp62_tools.py: do not edit> EMBOSS 6.6.0 programs (B09 WP62, F281): every tool builds the command line of the program and runs it with EmbossProgramRunner.', 'public partial class AnalysisTools', '{']
    for t in TOOLS:
        kind = t.get('kind', 'seqout')
        n = oracle_rows(t['program'])
        desc = "EMBOSS 6.6.0 %s: %s. Equals the installed program on %d runs (every input format of seqret, ranges, the output formats, the messages of the library; byte for byte). B09 WP62, F281." % (t['program'], t['desc'], n)
        t['_desc'] = desc
        t['_n'] = n
        out.append('    [McpServerTool(Name = %s, Title = %s, ReadOnly = true)]' % (cs_str(t['name']), cs_str(t['title'])))
        out.append('    [Description(%s)]' % cs_str(desc))
        params = []
        first = None
        for acd, p, d in t['inputs']:
            first = first or p
            params.append('        [Description(%s)] string %s' % (cs_str(d), p) + ('' if p == t['inputs'][0][1] else ' = ""'))
        for o in t['opts']:
            p, k, dflt, acd, d = o[:5]
            pn = pn_of(p)
            dv = {'flag': 'true' if dflt else 'false', 'int': 'null', 'float': 'null', 'string': 'null', 'menu': 'null'}[k]
            params.append('        [Description(%s)] %s %s = %s' % (cs_str(d), CS_TYPES[k], pn, dv))
        if t.get('random'):
            params.append('        [Description(%s)] int? microseconds = null' % cs_str(USEC[4]))
        if kind in ('seqout', 'outfile', 'outfile2'):
            params.append('        [Description("The input format (-sformat1); default: the program finds the format.")] string? sequenceFormat = null')
            params.append('        [Description(%s)] string? qualifiers = null' % cs_str("Sequence input qualifiers as typed after the sequence on the command line: -sbegin1 -send1 -sreverse1, -slower1 -supper1 -snucleotide1 -sprotein1 -scircular1, -sid1 name, -sdbname1 db, -squick1, -iquery1 query, -ioffset1 bytes."))
        if kind in ('seqout', 'noinput'):
            params.append('        [Description("The sequence output format (-osformat): fasta (default), embl, genbank, gcg, nbrf, ig, raw, msf, clustal, phylip, nexus, ... every format of seqret.")] string outputFormat = "fasta"')
        params.append('        CancellationToken cancellationToken = default')
        if kind != 'noinput':
            params.append('        [Description("More files by name (file name -> text) that the program may read (a file named by a USA, @file).")] IReadOnlyDictionary<string, string>? files = null')
        if kind in ('seqout', 'outfile'):
            params.append('        [Description("The USA of the sequence input (-sequence) as the program reads it, instead of the text: [format::]file[begin:end:r], file:id, @list, asis::SEQUENCE.")] string? sequenceValue = null')
        out.append('    public static EmbossProgramToolResult %s(' % method_name(t))
        out.append(',\n'.join(params) + ')')
        out.append('    {')
        if t['inputs']:
            out.append('        if (string.IsNullOrWhiteSpace(%s))' % first)
            out.append('            throw new ArgumentException("The sequence cannot be empty", nameof(%s));' % first)
        out.append('        cancellationToken.ThrowIfCancellationRequested();')
        out.append('        var all = new Dictionary<string, string>(%s);' % ('files ?? new Dictionary<string, string>()' if kind != 'noinput' else ''))
        out.append('        var args = new List<string>();')
        for i, (acd, p, d) in enumerate(t['inputs']):
            f = '%s.dat' % acd
            out.append('        all[%s] = %s;' % (cs_str(f), p))
            out.append('        args.Add("-%s");' % acd)
            out.append('        args.Add(SequenceSpec(%s, %s, %s));' % ('sequenceValue' if i == 0 and kind != 'outfile2' else 'null', p, cs_str(f)))
        if kind in ('seqout', 'outfile', 'outfile2'):
            out.append('        if (!string.IsNullOrWhiteSpace(sequenceFormat))')
            out.append('            args.AddRange(["-sformat1", sequenceFormat]);')
            out.append('        AddQualifierTokens(args, qualifiers);')
        for o in t['opts']:
            p, k, dflt, acd, d = o[:5]
            pn = pn_of(p)
            if len(o) > 5:
                out.append('        if (!string.IsNullOrEmpty(%s))\n        {\n            all[%s] = %s;\n            args.AddRange(["-%s", %s]);\n        }' % (pn, cs_str(o[5]), pn, acd, cs_str(o[5])) if acd != 'othersequence' else
                           '        if (!string.IsNullOrEmpty(%s))\n        {\n            all[%s] = %s;\n            args.AddRange(["-%s", %s.StartsWith("asis:", StringComparison.Ordinal) ? %s : %s]);\n        }' % (pn, cs_str(o[5]), pn, acd, pn, pn, cs_str(o[5])))
            elif k == 'flag':
                if dflt:
                    out.append('        if (!%s)\n            args.Add("-no%s");' % (pn, acd))
                else:
                    out.append('        if (%s)\n            args.Add("-%s");' % (pn, acd))
            elif k in ('int', 'float'):
                out.append('        if (%s.HasValue)\n            args.AddRange(["-%s", %s.Value.ToString(System.Globalization.CultureInfo.InvariantCulture)]);' % (pn, acd, pn))
            else:
                out.append('        if (!string.IsNullOrEmpty(%s))\n            args.AddRange(["-%s", %s]);' % (pn, acd, pn))
        param = first if t['inputs'] else 'amount'
        env = 'microseconds.HasValue ? new EmbossRunEnvironment(Microseconds: microseconds) : null' if t.get('random') else 'null'
        if kind in ('outfile', 'outfile2'):
            out.append('        return RunEmbossOutfileProgram(%s, args, all, nameof(%s), %s);' % (cs_str(t['program']), param, env))
        else:
            out.append('        return RunEmbossSequenceProgram(%s, args, all, outputFormat, nameof(%s), %s);' % (cs_str(t['program']), param, env))
        out.append('    }')
        out.append('')
    out.append('}')
    open('src/Seqeron/Mcp/Seqeron.Mcp.Analysis/Tools/AnalysisTools.EmbossSequencePrograms2.cs', 'w').write('\n'.join(out) + '\n')


def example_output(t, seq, opts):
    """the output of the installed program for the example of the tool (the clock of the example is fixed by the shim)"""
    kind = t.get('kind', 'seqout')
    shim_dir = tempfile.mkdtemp()
    shim = os.path.join(shim_dir, 'faketime.so')
    subprocess.run(['gcc', '-shared', '-fPIC', '-o', shim, 'scripts/alignment/faketime.c', '-ldl'], check=True)
    with tempfile.TemporaryDirectory() as d:
        args = []
        for i, (acd, p, desc) in enumerate(t['inputs']):
            f = '%s.dat' % acd
            open(os.path.join(d, f), 'w').write(seq if i == 0 else opts.get(p, ''))
            args += ['-' + acd, f]
        for o in t['opts']:
            p, k, dflt, acd, desc = o[:5]
            if p not in opts and k != 'flag':
                continue
            v = opts.get(p, dflt)
            if len(o) > 5:
                open(os.path.join(d, o[5]), 'w').write(v)
                v = o[5]
            if k == 'flag':
                if bool(v) != bool(dflt): args.append(('-no' if dflt else '-') + acd)
            else:
                args += ['-' + acd, str(v)]
        args += ['-outfile', 'out.txt'] if kind in ('outfile', 'outfile2') else ['-outseq', 'stdout']
        env = dict(os.environ, LD_PRELOAD=shim, FAKE_EPOCH='1700000000', TZ='UTC')
        if 'microseconds' in opts: env['FAKE_USEC'] = str(opts['microseconds'])
        r = subprocess.run([t['program'], '-auto'] + args, cwd=d, capture_output=True, text=True, env=env)
        text = open(os.path.join(d, 'out.txt'), newline='').read() if kind in ('outfile', 'outfile2') and os.path.exists(os.path.join(d, 'out.txt')) else r.stdout
        return {'text': text, 'files': {}, 'messages': [m for m in r.stderr.split('\n') if m]}


def gen_spec():
    spec = []
    for t in TOOLS:
        seq, opts = t['example']
        ex_in = {t['inputs'][0][1]: seq} if t['inputs'] else {}
        ex_in.update(opts)
        ex_in = {k: v for k, v in ex_in.items() if v is not None}
        src = {'shuffleseq': 'EmbossSequencePrograms.Random.cs', 'msbar': 'EmbossSequencePrograms.Random.cs', 'makenucseq': 'EmbossSequencePrograms.Random.cs', 'makeprotseq': 'EmbossSequencePrograms.Random.cs',
               'nthseqset': 'EmbossSequencePrograms.Sets.cs', 'listor': 'EmbossSequencePrograms.Sets.cs'}.get(t['program'], 'EmbossSequencePrograms.Report.cs')
        gen = 'generate_emboss_prog_random_oracle.py' if t.get('random') else 'generate_emboss_prog_report_oracle.py'
        spec.append({
            'name': t['name'], 'summary': "EMBOSS 6.6.0 %s: %s" % (t['program'], t['summary'].rstrip('.') + '.'),
            'methodId': 'EmbossProgramRunner.Run', 'sourceFile': src, 'sourcePath': 'Seqeron.Genomics.Analysis/' + src,
            'algorithmDoc': 'Alignment/Pairwise_Output_Formats_And_EMBOSS_Programs.md', 'oracleScript': gen,
            'short': "EMBOSS %s: %s; equals %s 6.6.0 on %d runs." % (t['program'], t['summary'].split(', as ')[0].rstrip('.').lower(), t['program'], t['_n']),
            'after': t['after'], 'file': 'AnalysisTools.EmbossSequencePrograms2.cs', 'description': t['_desc'],
            'output': {'text': ['string', 'The standard output of the program (the sequences) or the text of its output file.'], 'files': ['object', 'The files the program wrote, by name.'], 'messages': ['string[]', 'What the library printed on the error stream (warnings and errors).']},
            'errors': ['The sequence cannot be empty', 'The program died: the messages of the library (a value it refuses, an input no reader reads, a sequence of the wrong kind)', 'An input on which the reference program never stops or dies by a fault is refused'],
            'example': {'title': 'A run on a small input', 'prompt': 'Run %s on this input.' % t['program'], 'input': ex_in, 'output': example_output(t, seq, opts), 'note': 'The output is the text the installed program writes.'},
        })
    json.dump(spec, open('scripts/docs/emboss_wp62_tool_spec.json', 'w'), indent=1)


def gen_tests():
    out = ['using NUnit.Framework;', 'using Seqeron.Mcp.Analysis.Tools;', '', 'namespace Seqeron.Mcp.Analysis.Tests;', '',
           '// <auto-generated by scripts/docs/gen_emboss_wp62_tools.py: do not edit> B09 WP62 (F281): the MCP tools of the EMBOSS programs; the expected text is what the installed EMBOSS 6.6.0 program writes for the example of the tool document (the clock of the random programs fixed).',
           '[TestFixture]', 'public class EmbossProgramTools2Tests', '{']
    for t in TOOLS:
        seq, opts = t['example']
        ex = example_output(t, seq, opts)
        args = []
        if t['inputs']:
            args.append(cs_str(seq))
        for o in t['opts'] + ([USEC] if t.get('random') else []):
            p = o[0]
            if p not in opts:
                continue
            v = opts[p]
            lit = ('true' if v else 'false') if isinstance(v, bool) else cs_str(v) if isinstance(v, str) else str(v)
            args.append('%s: %s' % (pn_of(p), lit))
        if t['program'] == 'listor':
            args = [cs_str(seq), 'secondSequences: ' + cs_str(opts['secondSequences']), '@operator: ' + cs_str(opts['operator'])]
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
    out.append('        Assert.Throws<ArgumentException>(() => AnalysisTools.Shuffleseq(""));')
    out.append('        var e = Assert.Throws<ArgumentException>(() => AnalysisTools.Msbar(">s1\\nACGTACGT\\n", point: "7"));')
    out.append('        Assert.That(e!.Message, Does.Contain("Bad menu option \'7\'"));')
    out.append('        var f = Assert.Throws<ArgumentException>(() => AnalysisTools.Makenucseq(amount: 1, length: 2, insert: "ACGTA"));')
    out.append('        Assert.That(f!.Message, Does.Contain("Sequence smaller than inserted part"));')
    out.append('    }')
    out.append('}')
    open('tests/Seqeron/Seqeron.Mcp.Analysis.Tests/EmbossProgramTools2Tests.cs', 'w').write('\n'.join(out) + '\n')


if __name__ == '__main__':
    gen_cs()
    gen_spec()
    gen_tests()
