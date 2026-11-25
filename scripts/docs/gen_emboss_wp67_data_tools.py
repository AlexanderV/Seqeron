#!/usr/bin/env python3
"""HISTORIC (B09 WP80): the tools this script generated now live in the EMBOSS MCP servers by area (Seqeron.Mcp.Emboss*; copyleft, docs/Validation/review-2026-09/LICENSES.md) and the tool files are maintained by hand; the paths and class names below are those of WP61-WP67.
B09 WP67 (F309): the spec (scripts/docs/emboss_wp67_data_tool_spec.json) of the MCP tools drfinddata, dbtell and showdb of the Analysis server (src/Seqeron/Mcp/Seqeron.Mcp.Analysis/Tools/AnalysisTools.EmbossDataPrograms.cs, written by hand):
the example output of every tool is what the installed EMBOSS 6.6.0 program writes. Usage (repository root): python3 scripts/docs/gen_emboss_wp67_data_tools.py ; python3 scripts/docs/gen_mummer_est_tool_docs.py scripts/docs/emboss_wp67_data_tool_spec.json ;
python3 scripts/docs/register_mcp_tools.py scripts/docs/emboss_wp67_data_tool_spec.json ; python3 scripts/docs/sync_mcp_tool_docs.py ; python3 scripts/skills/gen-catalog.py"""
import json, os, re, subprocess, tempfile

CS = 'src/Seqeron/Mcp/Seqeron.Mcp.Analysis/Tools/AnalysisTools.EmbossDataPrograms.cs'
EXAMPLES = [
    ('drfinddata', {'query': 'microarray', 'subclasses': False, 'outputFormat': 'list'}, ['-query', 'microarray', '-nosubclasses', '-oformat', 'list', '-outfile', 'stdout'], 'Find the data resources whose output is a microarray data type'),
    ('dbtell', {'database': 'go_process'}, ['-database', 'go_process', '-outfile', 'stdout'], 'Tell me the definition of the database go_process'),
    ('showdb', {'database': 'taxon'}, ['-database', 'taxon', '-outfile', 'stdout'], 'Show the database taxon'),
]
INFO = {
    'drfinddata': ('EMBOSS drfinddata: finds the data resources of the DRCAT catalogue by the EDAM data term of their output; equals drfinddata 6.6.0 on 233 runs.', 'prosextract'),
    'dbtell': ('EMBOSS dbtell: the definition of a database of the standard EMBOSS configuration; equals dbtell 6.6.0 on 302 runs.', 'drfinddata'),
    'showdb': ('EMBOSS showdb: the table of the databases of the standard EMBOSS configuration; equals showdb 6.6.0 on 301 runs.', 'dbtell'),
}


def run(prog, args):
    with tempfile.TemporaryDirectory() as d:
        return subprocess.run([prog, '-auto'] + args, cwd=d, capture_output=True, text=True, env=dict(os.environ, HOME=d)).stdout


def main():
    cs = open(CS).read()
    desc = {m.group(1): m.group(2).replace('\\"', '"') for m in re.finditer(r'Name = "(\w+)".*?\n\s*\[Description\("((?:[^"\\]|\\.)*)"\)\]', cs)}
    spec = []
    for prog, inp, args, prompt in EXAMPLES:
        src = 'EmbossSequencePrograms.Drfinddata.cs' if prog == 'drfinddata' else 'EmbossSequencePrograms.Dbtell.cs'
        spec.append({
            'name': prog, 'summary': 'EMBOSS 6.6.0 %s: %s.' % (prog, INFO[prog][0].split(': ', 1)[1].split(';')[0].rstrip('.')), 'methodId': 'EmbossProgramRunner.Run', 'sourceFile': src, 'sourcePath': 'Seqeron.Genomics.Analysis/' + src,
            'algorithmDoc': 'Alignment/Pairwise_Output_Formats_And_EMBOSS_Programs.md', 'oracleScript': 'generate_emboss_prog_db_oracle.py', 'short': INFO[prog][0], 'after': INFO[prog][1], 'file': 'AnalysisTools.EmbossDataPrograms.cs', 'description': desc[prog],
            'output': {'text': ['string', 'The text the program writes to standard output.'], 'files': ['object', 'The other files the program wrote, by name (none).'], 'messages': ['string[]', 'What the library printed on the error stream (warnings and errors).']},
            'errors': ['The program died: the messages of the library (a value it refuses)'] + ([] if prog == 'dbtell' else ['An input on which the reference program dies by a fault is refused']),
            'example': {'title': 'A run on the standard configuration', 'prompt': prompt + '.', 'input': inp, 'output': {'text': run(prog, args), 'files': {}, 'messages': []}, 'note': 'The output is the text the installed program writes.'},
        })
    json.dump(spec, open('scripts/docs/emboss_wp67_data_tool_spec.json', 'w'), indent=1)


if __name__ == '__main__':
    main()
