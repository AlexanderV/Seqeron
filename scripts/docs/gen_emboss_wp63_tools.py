#!/usr/bin/env python3
"""HISTORIC (B09 WP80): the tools this script generated now live in the EMBOSS MCP servers by area (Seqeron.Mcp.Emboss*; copyleft, docs/Validation/review-2026-09/LICENSES.md) and the tool files are maintained by hand; the paths and class names below are those of WP61-WP67.
B09 WP63 (F282): generates the MCP tools of the EMBOSS programs of work package 63 on the Analysis server (tcode checktrans plotorf; the programs of later parts of the package are added to TOOLS):
src/Seqeron/Mcp/Seqeron.Mcp.Analysis/Tools/AnalysisTools.EmbossSequencePrograms3.cs (every tool builds the command line of the program and calls EmbossProgramRunner.Run), the spec
(scripts/docs/emboss_wp63_tool_spec.json) the document generators read, and the tests of the tools (the expected text is what the installed program writes for the example, run under the fixed clock of the example).
Usage (repository root): python3 scripts/docs/gen_emboss_wp63_tools.py ; python3 scripts/docs/gen_mummer_est_tool_docs.py scripts/docs/emboss_wp63_tool_spec.json ;
python3 scripts/docs/register_mcp_tools.py scripts/docs/emboss_wp63_tool_spec.json ; python3 scripts/docs/sync_mcp_tool_docs.py ; python3 scripts/skills/gen-catalog.py"""
import gzip, json, os, sys, subprocess, tempfile
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from gen_emboss_wp61_tools import cs_str, CS_TYPES, method_name  # noqa: E402

ORACLE_DIR = 'tests/Seqeron/Seqeron.Genomics.Tests/TestData/Alignment'
SEQ_IN = "The text of the sequence file (any input format of seqret), or a sequence as letters."
QUALS = "Sequence input qualifiers as typed after the sequence on the command line: -sbegin1 -send1 -sreverse1, -slower1 -supper1 -snucleotide1 -sprotein1 -scircular1, -sid1 name, -sdbname1 db, -squick1, -iquery1 query, -ioffset1 bytes."
RQUALS = "Report qualifiers as typed on the command line: -raccshow -rdesshow -rusashow -norscoreshow -norstrandshow -rmaxall N -rmaxseq N (-rformat is the parameter reportFormat)."
GQUALS = "Graph qualifiers as typed on the command line: -gtitle T -gsubtitle T -gxtitle T -gytitle T -gdesc T -goutfile NAME -gdirectory DIR."
# tool: name, title, program, kind (report | multi | graph), inputs, opts [(param, kind, default, acd, description[, file])], fixed = the output arguments, main = the file that is the text of the result
TOOLS = [
    dict(name='tcode', title='Genome — EMBOSS tcode (Fickett TESTCODE coding statistic)', program='tcode', after='showorf', kind='report', main='out.txt',
         summary="Identifies protein-coding regions of a nucleotide sequence with the Fickett TESTCODE statistic, as tcode of EMBOSS 6.6.0.",
         desc="the Fickett TESTCODE statistic of the windows of window bases (at least 200) every step bases, from the base compositions and the position asymmetries of the three codon positions with the probability tables of the data file Etcode.dat (or the datafile given); every window is reported with its score and the estimation Non-coding (below 0.74), No opinion or Coding (0.95 and above) as the report out.txt in the report format reportFormat (default table: every report format of EMBOSS, the report qualifiers of reportQualifiers); with plot the values are plotted against the window mid positions (graph: data writes the data file tcode1.dat, bedgraph the track file tcode.bedgraph, none draws nothing; the devices that draw with plplot (ps, png, svg ...) are refused)",
         inputs=[('sequence', 'sequence', SEQ_IN)],
         opts=[('window', 'int', None, 'window', 'The length of the sliding window (-window, default 200: at least 200 and at most the length of the first sequence).'),
               ('step', 'int', None, 'step', 'The step of the window (-step, default 3).'),
               ('dataFile', 'string', None, 'datafile', 'The text of a testcode data file in the format of Etcode.dat (-datafile; default the file of EMBOSS).', 'etcode.dat'),
               ('plot', 'flag', False, 'plot', 'Plot the statistic (-plot).'),
               ('graph', 'string', None, 'graph', 'The graph device of -plot (-graph): data, bedgraph or none.'),
               ('reportFormat', 'string', None, 'rformat', 'The report format (-rformat, default table): simple, embl, genbank, gff, gff3, pir, swiss, dasgff, debug, listfile, dbmotif, diffseq, draw, restrict, excel, feattable, motif, nametable, regions, seqtable, srs, summary, table, tagseq.'),
               ('reportQualifiers', 'qualifiers', None, '', RQUALS), ('graphQualifiers', 'qualifiers', None, '', GQUALS)],
         example=(">s1 test\n" + "ATGAAATTTGGGCCCGATGAACTGACCGGCTTAGGTACC" * 8 + "\n", {'window': 200, 'step': 30})),
    dict(name='checktrans', title='Genome — EMBOSS checktrans (STOP codons and ORFs of a protein)', program='checktrans', after='tcode', kind='multi', main='out.txt',
         summary="Reports the STOP codons and the open reading frames of a protein translation, as checktrans of EMBOSS 6.6.0.",
         desc="for every protein sequence the stretches between its stops (*) that are at least orfml residues long: the report out.txt lists their positions and the number of stops, the sequences name_N are returned in the file seqout.txt (sequence format outputFormat) and the protein feature table in the file features.txt (feature format featureFormat, default gff3); with addlast (default) a stop is forced at the end (it is not counted)",
         inputs=[('sequence', 'sequence', SEQ_IN)], fixed=['-outfile', 'out.txt', '-outseq', 'seqout.txt', '-outfeat', 'features.txt'],
         opts=[('orfml', 'int', None, 'orfml', 'The minimum length of an ORF to report (-orfml, default 100).'),
               ('addlast', 'flag', True, 'addlast', 'Force the sequence to end with a stop (-addlast, default true).'),
               ('outputFormat', 'string', None, 'osformat', 'The sequence output format (-osformat, default fasta).'),
               ('featureFormat', 'string', None, 'offormat', 'The feature output format (-offormat, default gff3): gff, gff2, gff3, bed, pir, swiss, dasgff, draw, debug, refseqp.')],
         example=(">p1 test\nMKVLAAGIV*ACDEFGHIK*PQRS\n", {'orfml': 3})),
    dict(name='wordfinder', title='Alignment — EMBOSS wordfinder (local alignments of large sequences by word matches)', program='wordfinder', after='checktrans', kind='multi', main='out.txt',
         summary="Matches large sequences against one or more other sequences by their words, as wordfinder of EMBOSS 6.6.0.",
         desc="the target sequences (bsequence) are indexed by their words of wordlen residues, every query sequence (asequence) is looked up in the index, the diagonal with the most word residues fixes the start points and a banded Smith-Waterman alignment of width cells (gap penalties gapopen / gapextend, matrix datafile: an EMBOSS matrix name such as EDNAFULL or EBLOSUM62, or the text of a matrix file) is written in the alignment format alignFormat (every -aformat of EMBOSS, default simple, with the sub-header lines Wordscore and Alignlength) as out.txt; the matches and alignments excluded by limitmatch / lowmatch / limitalign / lowalign and the matches per target are written to err.txt",
         inputs=[('asequence', 'querySequence', "The text of the query sequence file (-asequence: any input format of seqret; one or more sequences), or a sequence as letters."), ('bsequence', 'targetSequence', "The text of the target sequence file (-bsequence), or a sequence as letters.")],
         fixed=['-outfile', 'out.txt', '-errfile', 'err.txt'],
         opts=[('gapOpen', 'float', None, 'gapopen', 'The gap opening penalty (-gapopen, 0 to 1000, default 30).'), ('gapExtend', 'float', None, 'gapextend', 'The gap extension penalty (-gapextend, 0 to 10, default 1.5).'),
               ('dataFile', 'string', None, 'datafile', 'The scoring matrix (-datafile): the text of a matrix file.', 'matrix.txt'), ('matrix', 'string', None, 'datafile', 'The scoring matrix as an EMBOSS data file name (-datafile): EDNAFULL (nucleotides, default), EBLOSUM62 (proteins, default), EPAM250 ...'),
               ('width', 'int', None, 'width', 'The width of the alignment band (-width, at least 1, default 16).'), ('wordLength', 'int', None, 'wordlen', 'The word length of the initial matching (-wordlen, at least 3, default 6).'),
               ('limitMatch', 'int', None, 'limitmatch', 'The maximum match score, 0 for no limit (-limitmatch).'), ('limitAlign', 'int', None, 'limitalign', 'The maximum alignment length, 0 for no limit (-limitalign).'),
               ('lowMatch', 'int', None, 'lowmatch', 'The minimum match score (-lowmatch).'), ('lowAlign', 'int', None, 'lowalign', 'The minimum alignment length (-lowalign).'),
               ('alignFormat', 'string', None, 'aformat', 'The alignment format (-aformat, default simple): simple, pair, fasta, msf, clustal, mega, nexus, phylip, selex, treecon, markx0 ... markx10, match, multiple, score, srs, tcoffee, unknown, sam, ...'),
               ('alignWidth', 'int', None, 'awidth', 'The residues per line of the alignment (-awidth).'), ('alignQualifiers', 'qualifiers', None, '', 'Alignment qualifiers as typed on the command line: -aaccshow -adesshow -ausashow -aglobal.')],
         example=(">q1\nACGTACGTACGGGTTTTAAACCCGGGATATATCGCGCGTAGCTAGCTAGCTAGGGGATCCC\n", {'targetSequence': ">t1\nGGGACGTACGTACGGGTTTTAAACCCGGGATATATCGCGCGTAGCTAGCTAGCTAGGGGATCCCAAAAAA\n"})),
    dict(name='marscan', title='Genome — EMBOSS marscan (matrix/scaffold recognition signatures)', program='marscan', after='checktrans', kind='report', main='out.txt',
         summary="Finds matrix/scaffold recognition (MRS) signatures in DNA sequences, as marscan of EMBOSS 6.6.0.",
         desc="the MAR/SAR recognition signatures: the 16 bp site awwrtaannwwgnnnc (one mismatch allowed) or its reverse gnnncwwnnttaywwt and the 8 bp site aataayaa or ttrttatt within 200 bp of each other; of a cluster of signatures the pair with the least distance is reported as the report out.txt (default format gff: the feature misc_signal with the positions of the two sites; every report format of EMBOSS, the report qualifiers of reportQualifiers)",
         inputs=[('sequence', 'sequence', SEQ_IN)],
         opts=[('reportFormat', 'string', None, 'rformat', 'The report format (-rformat, default gff): simple, embl, genbank, gff, gff3, pir, swiss, dasgff, debug, listfile, dbmotif, diffseq, draw, restrict, excel, feattable, motif, nametable, regions, seqtable, srs, summary, table, tagseq.'),
               ('reportQualifiers', 'qualifiers', None, '', RQUALS)],
         example=(">s1 test\nCCCCCCCCCCAATAACAACCCCCCCCCCCCCCAATGTAAGGTTGCCCCTTTTTTTTTT\n", {})),
    dict(name='plotorf', title='Genome — EMBOSS plotorf (potential ORFs of a nucleotide sequence)', program='plotorf', after='checktrans', kind='graph', main='stdout',
         summary="Plots the potential open reading frames of a nucleotide sequence, as plotorf of EMBOSS 6.6.0.",
         desc="for each of the six frames the stretches from a start codon (start, default ATG) to the next stop codon of the frame (stop, default TAA,TAG,TGA; a stretch that reaches the end of the sequence ends there) as the filled rectangles of six plots; graph: data writes the six data files plotorf1.dat ... plotorf6.dat (returned in files), bedgraph one track file, none draws nothing (the devices that draw with plplot are refused)",
         inputs=[('sequence', 'sequence', SEQ_IN)],
         opts=[('start', 'string', None, 'start', 'The start codons, a comma list (-start, default ATG).'), ('stop', 'string', None, 'stop', 'The stop codons, a comma list (-stop, default TAA,TAG,TGA).'),
               ('graph', 'string', 'data', 'graph', 'The graph device (-graph): data (default), bedgraph or none.'), ('graphQualifiers', 'qualifiers', None, '', GQUALS)],
         example=(">o1\nATGAAATTTGGGCCCTAGATGCCCAAATTTTAAATGGGGGGGGTGATAGCATGAAACCCTGAGGATAA\n", {})),
]


def oracle_rows(prog):
    p = os.path.join(ORACLE_DIR, 'emboss_prog_%s_oracle.jsonl.gz' % prog)
    return sum(1 for _ in gzip.open(p, 'rt')) if os.path.exists(p) else 0


def gen_cs():
    out = ['using System.ComponentModel;', 'using ModelContextProtocol.Server;', 'using global::Seqeron.Genomics.Analysis;', '', 'namespace Seqeron.Mcp.Analysis.Tools;', '',
           '// <auto-generated by scripts/docs/gen_emboss_wp63_tools.py: do not edit> EMBOSS 6.6.0 programs (B09 WP63, F282): every tool builds the command line of the program and runs it with EmbossProgramRunner.', 'public partial class AnalysisTools', '{']
    for t in TOOLS:
        n = oracle_rows(t['program'])
        desc = "EMBOSS 6.6.0 %s: %s. Equals the installed program on %d runs (every input format of seqret, ranges, the report and graph options, the messages of the library; standard output, error stream, exit status and written files byte for byte). B09 WP63, F282." % (t['program'], t['desc'], n)
        t['_desc'], t['_n'] = desc, n
        out.append('    [McpServerTool(Name = %s, Title = %s, ReadOnly = true)]' % (cs_str(t['name']), cs_str(t['title'])))
        out.append('    [Description(%s)]' % cs_str(desc))
        params = []
        first = t['inputs'][0][1]
        for acd, p, d in t['inputs']:
            params.append('        [Description(%s)] string %s' % (cs_str(d), p))
        for o in t['opts']:
            p, k, dflt, acd, d = o[:5]
            if k == 'qualifiers':
                params.append('        [Description(%s)] string? %s = null' % (cs_str(d), p))
            else:
                dv = 'true' if (k == 'flag' and dflt) else 'false' if k == 'flag' else cs_str(dflt) if (k == 'string' and dflt) else 'null'
                ty = 'string' if k == 'string' and dflt else CS_TYPES[k]
                params.append('        [Description(%s)] %s %s = %s' % (cs_str(d), ty, p, dv))
        params.append('        [Description("The input format (-sformat1); default: the program finds the format.")] string? sequenceFormat = null')
        params.append('        [Description(%s)] string? qualifiers = null' % cs_str(QUALS))
        params.append('        CancellationToken cancellationToken = default')
        params.append('        [Description("More files by name (file name -> text) that the program may read (a file named by a USA, @file).")] IReadOnlyDictionary<string, string>? files = null')
        out.append('    public static EmbossProgramToolResult %s(' % method_name(t))
        out.append(',\n'.join(params) + ')')
        out.append('    {')
        out.append('        if (string.IsNullOrWhiteSpace(%s))' % first)
        out.append('            throw new ArgumentException("The sequence cannot be empty", nameof(%s));' % first)
        out.append('        cancellationToken.ThrowIfCancellationRequested();')
        out.append('        var all = new Dictionary<string, string>(files ?? new Dictionary<string, string>());')
        out.append('        var args = new List<string>();')
        for acd, p_, d_ in t['inputs']:
            out.append('        all["%s.dat"] = %s;' % (acd, p_))
            out.append('        args.AddRange(["-%s", SequenceSpec(null, %s, "%s.dat")]);' % (acd, p_, acd))
        out.append('        if (!string.IsNullOrWhiteSpace(sequenceFormat))')
        out.append('            args.AddRange(["-sformat1", sequenceFormat]);')
        out.append('        AddQualifierTokens(args, qualifiers);')
        for o in t['opts']:
            p, k, dflt, acd, d = o[:5]
            if k == 'qualifiers':
                out.append('        AddQualifierTokens(args, %s);' % p)
            elif len(o) > 5:
                out.append('        if (!string.IsNullOrEmpty(%s))\n        {\n            all[%s] = %s;\n            args.AddRange(["-%s", %s]);\n        }' % (p, cs_str(o[5]), p, acd, cs_str(o[5])))
            elif k == 'flag':
                out.append('        if (%s%s)\n            args.Add("-%s%s");' % ('!' if dflt else '', p, 'no' if dflt else '', acd))
            elif k in ('int', 'float'):
                out.append('        if (%s.HasValue)\n            args.AddRange(["-%s", %s.Value.ToString(System.Globalization.CultureInfo.InvariantCulture)]);' % (p, acd, p))
            else:
                out.append('        if (!string.IsNullOrEmpty(%s))\n            args.AddRange(["-%s", %s]);' % (p, acd, p))
        if t['kind'] == 'report' and any(o[0] == 'plot' for o in t['opts']):
            out.append('        if (plot && string.IsNullOrEmpty(graph))\n            args.AddRange(["-graph", "data"]);')
        fixed = t.get('fixed', ['-outfile', 'out.txt'] if t['kind'] == 'report' else [])
        if fixed:
            out.append('        args.AddRange([%s]);' % ', '.join(cs_str(x) for x in fixed))
        out.append('        return RunEmbossMultiProgram(%s, args, all, nameof(%s), %s);' % (cs_str(t['program']), first, cs_str(t['main'])))
        out.append('    }')
        out.append('')
    out.append('}')
    open('src/Seqeron/Mcp/Seqeron.Mcp.Analysis/Tools/AnalysisTools.EmbossSequencePrograms3.cs', 'w').write('\n'.join(out) + '\n')


def command(t, seq, opts, d):
    args = []
    for i, (acd, p_, desc) in enumerate(t['inputs']):
        open(os.path.join(d, acd + '.dat'), 'w').write(seq if i == 0 else opts[p_])
        args += ['-' + acd, acd + '.dat']
    for o in t['opts']:
        p, k, dflt, acd, desc = o[:5]
        if p not in opts and not (k == 'string' and dflt):
            continue
        v = opts.get(p, dflt)
        if k == 'qualifiers':
            args += v.split()
        elif len(o) > 5:
            open(os.path.join(d, o[5]), 'w').write(v)
            args += ['-' + acd, o[5]]
        elif k == 'flag':
            if bool(v) != bool(dflt): args.append(('-no' if dflt else '-') + acd)
        else:
            args += ['-' + acd, str(v)]
    if t['kind'] == 'report' and opts.get('plot') and not opts.get('graph'):
        args += ['-graph', 'data']
    args += t.get('fixed', ['-outfile', 'out.txt'] if t['kind'] == 'report' else [])
    return args


def example_output(t, seq, opts):
    """the output of the installed program for the example of the tool (the clock of the example is fixed by the shim)"""
    shim_dir = tempfile.mkdtemp()
    shim = os.path.join(shim_dir, 'faketime.so')
    subprocess.run(['gcc', '-shared', '-fPIC', '-o', shim, 'scripts/alignment/faketime.c', '-ldl'], check=True)
    with tempfile.TemporaryDirectory() as d:
        args = command(t, seq, opts, d)
        env = dict(os.environ, LD_PRELOAD=shim, FAKE_EPOCH='1700000000', TZ='UTC')
        r = subprocess.run([t['program'], '-auto'] + args, cwd=d, capture_output=True, text=True, env=env)
        files = {}
        for f in sorted(os.listdir(d)):
            if f not in [acd + '.dat' for acd, _, _ in t['inputs']] + ['faketime.so', 'etcode.dat'] and f != t['main']:
                files[f] = open(os.path.join(d, f), newline='').read()
        text = open(os.path.join(d, t['main']), newline='').read() if t['main'] != 'stdout' else r.stdout
        return {'text': text, 'files': files, 'messages': [m for m in r.stderr.split('\n') if m]}


def gen_spec():
    spec = []
    for t in TOOLS:
        seq, opts = t['example']
        ex_in = {t['inputs'][0][1]: seq}
        ex_in.update(opts)
        src = 'EmbossSequencePrograms.%s.cs' % t['program'].capitalize()
        spec.append({
            'name': t['name'], 'summary': "EMBOSS 6.6.0 %s: %s" % (t['program'], t['summary'].rstrip('.') + '.'),
            'methodId': 'EmbossProgramRunner.Run', 'sourceFile': src, 'sourcePath': 'Seqeron.Genomics.Analysis/' + src,
            'algorithmDoc': 'Alignment/Pairwise_Output_Formats_And_EMBOSS_Programs.md', 'oracleScript': 'generate_emboss_prog_wp63_oracle.py',
            'short': "EMBOSS %s: %s; equals %s 6.6.0 on %d runs." % (t['program'], t['summary'].split(', as ')[0].rstrip('.').lower(), t['program'], t['_n']),
            'after': t['after'], 'file': 'AnalysisTools.EmbossSequencePrograms3.cs', 'description': t['_desc'],
            'output': {'text': ['string', 'The text of the main output of the program (the report, or its standard output).'], 'files': ['object', 'The other files the program wrote (sequences, features, graph data), by name.'], 'messages': ['string[]', 'What the library printed on the error stream (warnings and errors).']},
            'errors': ['The sequence cannot be empty', 'The program died: the messages of the library (a value it refuses, an input no reader reads, a sequence of the wrong kind)', 'An input on which the reference program never stops or dies by a fault is refused'],
            'example': {'title': 'A run on a small input', 'prompt': 'Run %s on this input.' % t['program'], 'input': ex_in, 'output': example_output(t, seq, opts), 'note': 'The output is the text the installed program writes.'},
        })
    json.dump(spec, open('scripts/docs/emboss_wp63_tool_spec.json', 'w'), indent=1)


def gen_tests():
    out = ['using NUnit.Framework;', 'using Seqeron.Mcp.Analysis.Tools;', '', 'namespace Seqeron.Mcp.Analysis.Tests;', '',
           '// <auto-generated by scripts/docs/gen_emboss_wp63_tools.py: do not edit> B09 WP63 (F282): the MCP tools of the EMBOSS programs; the expected text is what the installed EMBOSS 6.6.0 program writes for the example of the tool document.',
           '[TestFixture]', 'public class EmbossProgramTools3Tests', '{',
           '    // the programs print the date of the run: the tools take it from the clock, so the date lines are compared as a placeholder',
           '    private static string Norm(string s) => System.Text.RegularExpressions.Regex.Replace(System.Text.RegularExpressions.Regex.Replace(s, "# Rundate: .*", "# Rundate: X"), "(#!Date |##date )\\\\d{4}-\\\\d{2}-\\\\d{2}", "$1X");', '']
    for t in TOOLS:
        seq, opts = t['example']
        ex = example_output(t, seq, opts)
        args = [cs_str(seq)]
        for acd, p_, d_ in t['inputs'][1:]:
            args.append('%s: %s' % (p_, cs_str(opts[p_])))
        for o in t['opts']:
            p = o[0]
            if p not in opts:
                continue
            v = opts[p]
            lit = ('true' if v else 'false') if isinstance(v, bool) else cs_str(v) if isinstance(v, str) else str(v)
            args.append('%s: %s' % (p, lit))
        out.append('    [Test]')
        out.append('    public void %s_WritesWhatTheInstalledProgramWrites()' % method_name(t))
        out.append('    {')
        out.append('        var r = AnalysisTools.%s(%s);' % (method_name(t), ', '.join(args)))
        norm = lambda x: __import__('re').sub(r'(#!Date |##date )\d{4}-\d{2}-\d{2}', r'\1X', __import__('re').sub(r'# Rundate: .*', '# Rundate: X', x))
        out.append('        Assert.That(Norm(r.Text), Is.EqualTo(%s));' % cs_str(norm(ex['text'])))
        for f, v in ex['files'].items():
            out.append('        Assert.That(Norm(r.Files[%s]), Is.EqualTo(%s));' % (cs_str(f), cs_str(norm(v))))
        out.append('    }')
        out.append('')
    out.append('    [Test]')
    out.append('    public void Tools_RefuseAnEmptySequenceAndABadValueWithTheMessagesOfTheLibrary()')
    out.append('    {')
    out.append('        Assert.Throws<ArgumentException>(() => AnalysisTools.Tcode(""));')
    out.append('        var e = Assert.Throws<ArgumentException>(() => AnalysisTools.Plotorf(">o1\\nATGAAATAG\\n", graph: "bogus"));')
    out.append('        Assert.That(e!.Message, Does.Contain("Invalid XY graph value \'bogus\'"));')
    out.append('        var f = Assert.Throws<ArgumentException>(() => AnalysisTools.Tcode(">s1\\n" + new string(\'A\', 300) + "\\n", plot: true, graph: "png"));')
    out.append('        Assert.That(f!.Message, Does.Contain("plplot"));')
    out.append('    }')
    out.append('}')
    open('tests/Seqeron/Seqeron.Mcp.EmbossOpenReadingFrames.Tests/EmbossProgramTools3Tests.cs', 'w').write('\n'.join(out) + '\n')


if __name__ == '__main__':
    gen_cs()
    gen_spec()
    gen_tests()
