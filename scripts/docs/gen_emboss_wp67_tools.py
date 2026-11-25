#!/usr/bin/env python3
"""HISTORIC (B09 WP80): the tools this script generated now live in the EMBOSS MCP servers by area (Seqeron.Mcp.Emboss*; copyleft, docs/Validation/review-2026-09/LICENSES.md) and the tool files are maintained by hand; the paths and class names below are those of WP61-WP67.
B09 WP67 (F306): generates the MCP tools of the PROSITE programs of EMBOSS 6.6.0 on the Analysis server (patmatmotifs, patmatdb, prosextract):
src/Seqeron/Mcp/Seqeron.Mcp.Analysis/Tools/AnalysisTools.EmbossSequencePrograms4.cs (every tool builds the command line of the program and calls EmbossProgramRunner.Run), the spec
(scripts/docs/emboss_wp67_tool_spec.json) the document generators read, and the tests of the tools (the expected text is what the installed program writes for the example, run under the fixed clock of the example).
Usage (repository root): python3 scripts/docs/gen_emboss_wp67_tools.py ; python3 scripts/docs/gen_mummer_est_tool_docs.py scripts/docs/emboss_wp67_tool_spec.json ;
python3 scripts/docs/register_mcp_tools.py scripts/docs/emboss_wp67_tool_spec.json ; python3 scripts/docs/sync_mcp_tool_docs.py ; python3 scripts/skills/gen-catalog.py"""
import gzip, json, os, sys, subprocess, tempfile
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from gen_emboss_wp61_tools import cs_str, CS_TYPES, method_name  # noqa: E402

ORACLE_DIR = 'tests/Seqeron/Seqeron.Genomics.Tests/TestData/Alignment'
SEQ_IN = "The text of the sequence file (any input format of seqret), or a sequence as letters."
QUALS = "Sequence input qualifiers as typed after the sequence on the command line: -sbegin1 -send1 -sreverse1, -slower1 -supper1 -snucleotide1 -sprotein1 -scircular1, -sid1 name, -sdbname1 db, -squick1, -iquery1 query, -ioffset1 bytes."
RQUALS = "Report qualifiers as typed on the command line: -raccshow -rdesshow -rusashow -norscoreshow -norstrandshow -rmaxall N -rmaxseq N (-rformat is the parameter reportFormat)."
GQUALS = "Graph qualifiers as typed on the command line: -gtitle T -gsubtitle T -gxtitle T -gytitle T -gdesc T -goutfile NAME -gdirectory DIR."
# tool: name, title, program, kind (report | multi | graph | prosextract), inputs, required, opts [(param, kind, default, acd, description[, file])], fixed = the output arguments, main = the file that is the text of the result
PROT_IN = "The text of the protein sequence file (any input format of seqret), or a sequence as letters."
TOOLS = [
    dict(name='patmatmotifs', title='Protein — EMBOSS patmatmotifs (scan a protein with the PROSITE database)', program='patmatmotifs', after='marscan', kind='report', main='out.txt',
         summary="Scans a protein sequence with the patterns of the PROSITE database, as patmatmotifs of EMBOSS 6.6.0.",
         desc="the protein sequence is searched with every pattern of the PROSITE database in the form prosextract makes (prositeLines: the text of PROSITE/prosite.lines: identifier, accession, description word, pattern and regular expression of every PATTERN entry; the tool prosextract makes it from prosite.dat), the hits are the report out.txt (default format dbmotif: every report format of EMBOSS, the report qualifiers of reportQualifiers); the simple patterns of post-translational modification sites (myristyl, asn_glycosylation, camp_phospho_site, pkc_phospho_site, ck2_phospho_site, tyr_phospho_site) are skipped unless noprune; with full the documentation of every pattern that has a hit is appended (documentation: accession -> the text of the file PROSITE/PSnnnnn that prosextract makes)",
         inputs=[('sequence', 'sequence', PROT_IN)],
         opts=[('prositeLines', 'datafile', None, '', 'The text of the PROSITE data file prosite.lines made by prosextract (PROSITE/prosite.lines): required, the sandbox of the tool has no PROSITE database.', 'PROSITE/prosite.lines'),
               ('documentation', 'docs', None, '', 'The documentation files of the PROSITE entries made by prosextract: accession (PS00001) -> text; read with full.'),
               ('full', 'flag', False, 'full', 'Provide the full documentation of the patterns that match (-full).'),
               ('prune', 'flag', True, 'prune', 'Ignore the simple patterns (-prune, default true; false = -noprune).'),
               ('reportFormat', 'string', None, 'rformat', 'The report format (-rformat, default dbmotif): simple, embl, genbank, gff, gff3, pir, swiss, dasgff, debug, listfile, dbmotif, diffseq, draw, restrict, excel, feattable, motif, nametable, regions, seqtable, srs, summary, table, tagseq.'),
               ('reportQualifiers', 'qualifiers', None, '', RQUALS)],
         example=(">s1 test\nMAKKCAGGCLLAKKCDEFGHIKLMNPQRSTVWY\n", {'prositeLines': "MYPAT PS99999\n Test\n A-x(2)-C\n^A[^BJOUXZ]{2}C\n\n", 'full': True, 'documentation': {'PS99999': 'MYPAT documentation line 1\nline 2\n'}})),
    dict(name='patmatdb', title='Protein — EMBOSS patmatdb (search protein sequences with a PROSITE pattern)', program='patmatdb', after='patmatmotifs', kind='report', main='out.txt',
         summary="Searches protein sequences with a PROSITE sequence motif, as patmatdb of EMBOSS 6.6.0.",
         desc="every protein sequence is searched with the PROSITE pattern motif (case-independent; < > anchor, x any residue, [ ] a class, { } excluded residues, (n) (n,m) counts) turned into a regular expression of the library (PCRE 7); the hits are the report out.txt (default format dbmotif: every report format of EMBOSS, the report qualifiers of reportQualifiers); a text that the PROSITE syntax refuses (an unmatched bracket, a character that is no residue) stops the program with the message of the library, a pattern that the regular expression library refuses (a range out of order, a count above 65535, nothing to repeat) makes the reference program die by a segmentation fault and is refused with the message of the library",
         inputs=[('sequence', 'sequence', PROT_IN)], required=[('motif', 'The PROSITE pattern to search for (-motif, at least 2 characters), e.g. [DE](2)HS{P}X(2)PX(2,4)C.')],
         opts=[('reportFormat', 'string', None, 'rformat', 'The report format (-rformat, default dbmotif): simple, embl, genbank, gff, gff3, pir, swiss, dasgff, debug, listfile, dbmotif, diffseq, draw, restrict, excel, feattable, motif, nametable, regions, seqtable, srs, summary, table, tagseq.'),
               ('reportQualifiers', 'qualifiers', None, '', RQUALS)],
         example=(">p1 test\nMAAAAKKCAGGCLLAKKC\n", {'motif': 'A-x(2)-C'})),
    dict(name='prosextract', title='Protein — EMBOSS prosextract (build the PROSITE database for patmatmotifs)', program='prosextract', after='patmatdb', kind='prosextract', main='/usr/share/EMBOSS/data/PROSITE/prosite.lines',
         summary="Builds the PROSITE motif database for patmatmotifs from prosite.dat and prosite.doc, as prosextract of EMBOSS 6.6.0.",
         desc="for every PATTERN entry of prosite.dat (the text of the SIB file) the file PROSITE/prosite.lines (returned as the text) holds the identifier, the accession, the first word of the description, the pattern and its regular expression; for every {PSnnnnn; NAME} entry of prosite.doc one documentation file PROSITE/PSnnnnn (returned in files, by the path of the data directory /usr/share/EMBOSS/data/PROSITE/PSnnnnn) holds the text between {BEGIN} and {END}",
         inputs=[('prosite.dat', 'prositeDat', 'The text of the PROSITE data file prosite.dat (ftp.expasy.org/databases/prosite).'), ('prosite.doc', 'prositeDoc', 'The text of the PROSITE documentation file prosite.doc.')],
         opts=[],
         example=("ID   A_PAT; PATTERN.\nAC   PS99999;\nDT   01-APR-1990 CREATED.\nDE   Test pattern site.\nPA   A-x(2)-C.\n//\n", {'prositeDoc': "{PDOC99999}\n{PS99999; A_PAT}\n{BEGIN}\nDocumentation of A_PAT.\n{END}\n"})),
]


def oracle_rows(prog):
    p = os.path.join(ORACLE_DIR, 'emboss_prog_%s_oracle.jsonl.gz' % prog)
    return sum(1 for _ in gzip.open(p, 'rt')) if os.path.exists(p) else 0


def gen_cs():
    out = ['using System.ComponentModel;', 'using ModelContextProtocol.Server;', 'using global::Seqeron.Genomics.Analysis;', '', 'namespace Seqeron.Mcp.Analysis.Tools;', '',
           '// <auto-generated by scripts/docs/gen_emboss_wp67_tools.py: do not edit> PROSITE programs of EMBOSS 6.6.0 (B09 WP67, F306): every tool builds the command line of the program and runs it with EmbossProgramRunner.', 'public partial class AnalysisTools', '{']
    for t in TOOLS:
        n = oracle_rows(t['program'])
        desc = "EMBOSS 6.6.0 %s: %s. Equals the installed program on %d runs (every input format of seqret, ranges, the report options, the messages of the library; standard output, error stream, exit status and written files byte for byte). B09 WP67, F306." % (t['program'], t['desc'], n)
        t['_desc'], t['_n'] = desc, n
        out.append('    [McpServerTool(Name = %s, Title = %s, ReadOnly = true)]' % (cs_str(t['name']), cs_str(t['title'])))
        out.append('    [Description(%s)]' % cs_str(desc))
        params = []
        first = t['inputs'][0][1]
        for acd, p, d in t['inputs']:
            params.append('        [Description(%s)] string %s' % (cs_str(d), p))
        for p, d in t.get('required', []):
            params.append('        [Description(%s)] string %s' % (cs_str(d), p))
        for o in t['opts']:
            p, k, dflt, acd, d = o[:5]
            if k in ('qualifiers', 'datafile'):
                params.append('        [Description(%s)] string? %s = null' % (cs_str(d), p))
            elif k == 'docs':
                params.append('        [Description(%s)] IReadOnlyDictionary<string, string>? %s = null' % (cs_str(d), p))
            else:
                dv = 'true' if (k == 'flag' and dflt) else 'false' if k == 'flag' else cs_str(dflt) if (k == 'string' and dflt) else 'null'
                ty = 'string' if k == 'string' and dflt else CS_TYPES[k]
                params.append('        [Description(%s)] %s %s = %s' % (cs_str(d), ty, p, dv))
        if t['kind'] != 'prosextract':
            params.append('        [Description("The input format (-sformat1); default: the program finds the format.")] string? sequenceFormat = null')
            params.append('        [Description(%s)] string? qualifiers = null' % cs_str(QUALS))
        params.append('        CancellationToken cancellationToken = default')
        params.append('        [Description("More files by name (file name -> text) that the program may read (a file named by a USA, @file, the PROSITE data files of the data directory).")] IReadOnlyDictionary<string, string>? files = null')
        out.append('    public static EmbossProgramToolResult %s(' % method_name(t))
        out.append(',\n'.join(params) + ')')
        out.append('    {')
        if t['kind'] == 'prosextract':
            out.append('        if (string.IsNullOrWhiteSpace(prositeDat))\n            throw new ArgumentException("The PROSITE data file cannot be empty", nameof(prositeDat));')
            out.append('        cancellationToken.ThrowIfCancellationRequested();')
            out.append('        var all = new Dictionary<string, string>(files ?? new Dictionary<string, string>())')
            out[-1] += ' { ["prosite.dat"] = prositeDat, ["prosite.doc"] = prositeDoc ?? "" };'
            out.append('        var args = new List<string> { "-prositedir", "./" };')
            out.append('        return RunEmbossMultiProgram("prosextract", args, all, nameof(prositeDat), %s);' % cs_str(t['main']))
            out.append('    }')
            out.append('')
            continue
        out.append('        if (string.IsNullOrWhiteSpace(%s))' % first)
        out.append('            throw new ArgumentException("The sequence cannot be empty", nameof(%s));' % first)
        for p, d in t.get('required', []):
            out.append('        if (string.IsNullOrWhiteSpace(%s))' % p)
            out.append('            throw new ArgumentException("The %s cannot be empty", nameof(%s));' % (p, p))
        out.append('        cancellationToken.ThrowIfCancellationRequested();')
        out.append('        var all = new Dictionary<string, string>(files ?? new Dictionary<string, string>());')
        out.append('        var args = new List<string>();')
        for acd, p_, d_ in t['inputs']:
            out.append('        all["%s.dat"] = %s;' % (acd, p_))
            out.append('        args.AddRange(["-%s", SequenceSpec(null, %s, "%s.dat")]);' % (acd, p_, acd))
        for p, d in t.get('required', []):
            out.append('        args.AddRange(["-%s", %s]);' % (p, p))
        out.append('        if (!string.IsNullOrWhiteSpace(sequenceFormat))')
        out.append('            args.AddRange(["-sformat1", sequenceFormat]);')
        out.append('        AddQualifierTokens(args, qualifiers);')
        for o in t['opts']:
            p, k, dflt, acd, d = o[:5]
            if k == 'qualifiers':
                out.append('        AddQualifierTokens(args, %s);' % p)
            elif k == 'datafile':
                out.append('        if (!string.IsNullOrEmpty(%s))\n            all[%s] = %s;' % (p, cs_str(o[5]), p))
            elif k == 'docs':
                out.append('        if (%s is not null)\n            foreach (var kv in %s)\n                all["PROSITE/" + kv.Key] = kv.Value;' % (p, p))
            elif k == 'flag':
                out.append('        if (%s%s)\n            args.Add("-%s%s");' % ('!' if dflt else '', p, 'no' if dflt else '', acd))
            elif k in ('int', 'float'):
                out.append('        if (%s.HasValue)\n            args.AddRange(["-%s", %s.Value.ToString(System.Globalization.CultureInfo.InvariantCulture)]);' % (p, acd, p))
            else:
                out.append('        if (!string.IsNullOrEmpty(%s))\n            args.AddRange(["-%s", %s]);' % (p, acd, p))
        fixed = t.get('fixed', ['-outfile', 'out.txt'] if t['kind'] == 'report' else [])
        if fixed:
            out.append('        args.AddRange([%s]);' % ', '.join(cs_str(x) for x in fixed))
        out.append('        return RunEmbossMultiProgram(%s, args, all, nameof(%s), %s);' % (cs_str(t['program']), first, cs_str(t['main'])))
        out.append('    }')
        out.append('')
    out.append('}')
    open('src/Seqeron/Mcp/Seqeron.Mcp.Analysis/Tools/AnalysisTools.EmbossSequencePrograms4.cs', 'w').write('\n'.join(out) + '\n')


def command(t, seq, opts, d):
    """the command line and the files of the example in the directory d"""
    args = []
    for i, (acd, p_, desc) in enumerate(t['inputs']):
        text = seq if i == 0 else opts.get(p_, '')
        open(os.path.join(d, acd if t['kind'] == 'prosextract' else acd + '.dat'), 'w').write(text)
        if t['kind'] != 'prosextract':
            args += ['-' + acd, acd + '.dat']
    for p, desc in t.get('required', []):
        args += ['-' + p, opts[p]]
    if t['kind'] == 'prosextract':
        return ['-prositedir', './']
    for o in t['opts']:
        p, k, dflt, acd, desc = o[:5]
        if p not in opts and not (k == 'string' and dflt):
            continue
        v = opts.get(p, dflt)
        if k == 'qualifiers':
            args += v.split()
        elif k == 'datafile':
            os.makedirs(os.path.dirname(os.path.join(d, o[5])), exist_ok=True)
            open(os.path.join(d, o[5]), 'w').write(v)
        elif k == 'docs':
            for acc, text in v.items():
                open(os.path.join(d, 'PROSITE', acc), 'w').write(text)
        elif k == 'flag':
            if bool(v) != bool(dflt): args.append(('-no' if dflt else '-') + acd)
        else:
            args += ['-' + acd, str(v)]
    args += t.get('fixed', ['-outfile', 'out.txt'] if t['kind'] == 'report' else [])
    return args


def example_output(t, seq, opts):
    """the output of the installed program for the example of the tool (the clock of the example is fixed by the shim; the PROSITE data of the installation are hidden: the data are the files of the example)"""
    shim_dir = tempfile.mkdtemp()
    shim = os.path.join(shim_dir, 'faketime.so')
    subprocess.run(['gcc', '-shared', '-fPIC', '-o', shim, 'scripts/alignment/faketime.c', '-ldl'], check=True)
    with tempfile.TemporaryDirectory() as d:
        os.makedirs(os.path.join(d, 'PROSITE'), exist_ok=True)
        args = command(t, seq, opts, d)
        env = dict(os.environ, LD_PRELOAD=shim, FAKE_EPOCH='1700000000', TZ='UTC')
        inputs = {acd if t['kind'] == 'prosextract' else acd + '.dat' for acd, _, _ in t['inputs']}
        before = set()
        for root, _, names in os.walk(d):
            before.update(os.path.relpath(os.path.join(root, n), d) for n in names)
        if t['kind'] == 'prosextract':
            os.makedirs(os.path.join(d, 'edata', 'PROSITE'))
            env['EMBOSS_DATA'] = os.path.join(d, 'edata') + '/'
        cmd = ['unshare', '-m', 'sh', '-c', 'mount -t tmpfs none /usr/share/EMBOSS/data/PROSITE && exec "$@"', 'sh', t['program'], '-auto'] + args
        r = subprocess.run(cmd, cwd=d, capture_output=True, text=True, env=env)
        files = {}
        for root, _, names in os.walk(d):
            for n in sorted(names):
                rel = os.path.relpath(os.path.join(root, n), d)
                if rel in before or rel == 'faketime.so':
                    continue
                files[rel] = open(os.path.join(root, n), newline='').read()
        if t['kind'] == 'prosextract':
            files = {'/usr/share/EMBOSS/data/PROSITE/' + k[len('edata/PROSITE/'):]: v for k, v in files.items()}
            text = files.pop(t['main'], '')
        else:
            text = files.pop(t['main'], '') if t['main'] != 'stdout' else r.stdout
        return {'text': text, 'files': files, 'messages': [m for m in r.stderr.split('\n') if m]}


def gen_spec():
    spec = []
    for t in TOOLS:
        seq, opts = t['example']
        ex_in = {t['inputs'][0][1]: seq}
        ex_in.update(opts)
        src = 'EmbossSequencePrograms.Prosite.cs'
        spec.append({
            'name': t['name'], 'summary': "EMBOSS 6.6.0 %s: %s" % (t['program'], t['summary'].rstrip('.') + '.'),
            'methodId': 'EmbossProgramRunner.Run', 'sourceFile': src, 'sourcePath': 'Seqeron.Genomics.Analysis/' + src,
            'algorithmDoc': 'Alignment/Pairwise_Output_Formats_And_EMBOSS_Programs.md', 'oracleScript': 'generate_emboss_prog_prosite_oracle.py',
            'short': "EMBOSS %s: %s; equals %s 6.6.0 on %d runs." % (t['program'], t['summary'].split(', as ')[0].rstrip('.').lower(), t['program'], t['_n']),
            'after': t['after'], 'file': 'AnalysisTools.EmbossSequencePrograms4.cs', 'description': t['_desc'],
            'output': {'text': ['string', 'The text of the main output of the program (the report, or its standard output).'], 'files': ['object', 'The other files the program wrote (sequences, features, graph data), by name.'], 'messages': ['string[]', 'What the library printed on the error stream (warnings and errors).']},
            'errors': ['The sequence cannot be empty', 'The program died: the messages of the library (a value it refuses, an input no reader reads, a sequence of the wrong kind)', 'An input on which the reference program never stops or dies by a fault is refused'],
            'example': {'title': 'A run on a small input', 'prompt': 'Run %s on this input.' % t['program'], 'input': ex_in, 'output': example_output(t, seq, opts), 'note': 'The output is the text the installed program writes.'},
        })
    json.dump(spec, open('scripts/docs/emboss_wp67_tool_spec.json', 'w'), indent=1)


def cs_lit(v):
    if isinstance(v, bool):
        return 'true' if v else 'false'
    if isinstance(v, str):
        return cs_str(v)
    if isinstance(v, dict):
        return 'new Dictionary<string, string> { ' + ', '.join('[%s] = %s' % (cs_str(k), cs_str(x)) for k, x in v.items()) + ' }'
    return str(v)


def gen_tests():
    out = ['using NUnit.Framework;', 'using Seqeron.Mcp.Analysis.Tools;', '', 'namespace Seqeron.Mcp.Analysis.Tests;', '',
           '// <auto-generated by scripts/docs/gen_emboss_wp67_tools.py: do not edit> B09 WP67 (F306): the MCP tools of the PROSITE programs; the expected text is what the installed EMBOSS 6.6.0 program writes for the example of the tool document.',
           '[TestFixture]', 'public class EmbossProgramTools4Tests', '{',
           '    // the programs print the date of the run: the tools take it from the clock, so the date lines are compared as a placeholder',
           '    private static string Norm(string s) => System.Text.RegularExpressions.Regex.Replace(System.Text.RegularExpressions.Regex.Replace(s, "# Rundate: .*", "# Rundate: X"), "(#!Date |##date )\\\\d{4}-\\\\d{2}-\\\\d{2}", "$1X");', '']
    for t in TOOLS:
        seq, opts = t['example']
        ex = example_output(t, seq, opts)
        args = [cs_str(seq)]
        for acd, p_, d_ in t['inputs'][1:]:
            args.append('%s: %s' % (p_, cs_str(opts[p_])))
        for p, d_ in t.get('required', []):
            args.append('%s: %s' % (p, cs_str(opts[p])))
        for o in t['opts']:
            p = o[0]
            if p in opts:
                args.append('%s: %s' % (p, cs_lit(opts[p])))
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
    out.append('    public void Tools_RefuseAnEmptyInputAndABadPatternWithTheMessagesOfTheLibrary()')
    out.append('    {')
    out.append('        Assert.Throws<ArgumentException>(() => AnalysisTools.Patmatmotifs("", prositeLines: "x"));')
    out.append('        Assert.Throws<ArgumentException>(() => AnalysisTools.Patmatdb(">p\\nMAAAK\\n", ""));')
    out.append('        Assert.Throws<ArgumentException>(() => AnalysisTools.Prosextract("", ""));')
    out.append('        var e = Assert.Throws<ArgumentException>(() => AnalysisTools.Patmatdb(">p\\nMAAAK\\n", "A-x(2"));')
    out.append('        Assert.That(e!.Message, Does.Contain("Unmatched \'(\' in A-X(2"));')
    out.append('        var f = Assert.Throws<ArgumentException>(() => AnalysisTools.Patmatdb(">p\\nMAAAK\\n", "A(2,1)"));')
    out.append('        Assert.That(f!.Message, Does.Contain("numbers out of order in {} quantifier"));')
    out.append('        var g = Assert.Throws<ArgumentException>(() => AnalysisTools.Patmatmotifs(">p\\nMAAAK\\n"));')
    out.append('        Assert.That(g!.Message, Does.Contain("PROSEXTRACT needs running"));')
    out.append('    }')
    out.append('}')
    open('tests/Seqeron/Seqeron.Mcp.EmbossMotifs.Tests/EmbossProgramTools4Tests.cs', 'w').write('\n'.join(out) + '\n')


if __name__ == '__main__':
    gen_cs()
    gen_spec()
    gen_tests()
