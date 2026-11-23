#!/usr/bin/env python3
"""Regenerates (or checks) the MCP tool documents of B09's tools from the code's tool metadata (audit A5-5, B09 F170).

SOURCE OF TRUTH: the `[McpServerTool(Name)]` + `[Description]` attributes of `src/Seqeron/Mcp/**/Tools/*.cs`, read by REFLECTION (scripts/docs/McpToolDump: tool name, description, every parameter with its
type, default value and `[Description]`), never by text matching.

WRITES (only what the code defines; everything else in the documents - output schema, examples, errors, worked example - is hand text and stays):
  docs/mcp/tools/<server>/<tool>.mcp.json : `description`, `inputSchema.properties` (names, order, type, `default`, `description`; `items` of a new array parameter), `inputSchema.required`
  docs/mcp/tools/<server>/<tool>.md       : the tool description (under the title and under "## Description") and the rows of the "## Input Schema" table (a row of a parameter the code no longer
                                            has is removed, a new one is added, the description cell is the code's text + "(default: X)")

SCOPE: every tool of the Alignment server and the Analysis-server tools whose document names a B09 class in `methodId` (SequenceAligner, GenomicAnalyzer, MultipleAlignment*, MotifFinder consambig ...);
`--all-owned` is the default, `--server NAME` / `--tool NAME` narrow it.

  python3 scripts/docs/sync_mcp_tool_docs.py            regenerate
  python3 scripts/docs/sync_mcp_tool_docs.py --check    exit 1 and list every difference (CI / the contract test of Seqeron.Mcp.*.Tests does the same over reflection)
  python3 scripts/docs/sync_mcp_tool_docs.py --dump f   use a metadata dump written earlier (skips the dotnet run)
Needs the MCP projects built (`dotnet build src/Seqeron/Mcp/Seqeron.Mcp.Alignment src/Seqeron/Mcp/Seqeron.Mcp.Analysis -c Debug`); python3 stdlib only.
"""
import argparse
import glob
import json
import os
import re
import subprocess
import sys
import tempfile

ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir, os.pardir))
TOOLS = os.path.join(ROOT, 'docs', 'mcp', 'tools')
DUMPER = os.path.join(ROOT, 'scripts', 'docs', 'McpToolDump')
B09_CLASSES = ('EmbossReportPrograms.', 'SequenceAligner.', 'GenomicAnalyzer.', 'MultipleAlignment', 'AnchorBasedAligner.', 'DotMatcher.', 'Minimap2Chaining.', 'NucmerClusters.', 'Ksw2Aligner.', 'AlignmentFormats.', 'MashCommands.', 'MashTaxScreen.', 'MotifFinder.GenerateEmbossAmbiguousConsensus', 'Mummer', 'PromerProgram.', 'EmbossEst2genomePrograms.', 'EmbossFeaturePrograms.', 'EmbossSeqretPrograms.', 'EmbossMaskfeatPrograms.', 'EmbossDescseqPrograms.', 'EmbossUnionProgram.', 'EmbossEntretProgram.', 'EmbossExtractfeatProgram.', 'EmbossCoderetProgram.', 'EmbossShowfeatProgram.', 'EmbossFeatreportProgram.', 'EmbossTwofeatProgram.', 'EmbossProgramRunner.', 'SummaryInfo', 'ConservationScores', 'CodonAligner.', 'MultipleSeqAlignment')


def run_dump(servers):
    bins = ':'.join(os.path.join(ROOT, 'src', 'Seqeron', 'Mcp', 'Seqeron.Mcp.' + s, 'bin', 'Debug', 'net10.0') for s in servers)
    out = tempfile.mkdtemp(prefix='mcp-tool-dump-')
    subprocess.run(['dotnet', 'build', DUMPER, '-c', 'Release', '-o', out, '-v', 'q', '-nologo'], check=True, stdout=subprocess.DEVNULL)
    r = subprocess.run(['dotnet', os.path.join(out, 'McpToolDump.dll'), bins] + ['Seqeron.Mcp.' + s for s in servers], check=True, capture_output=True, text=True, encoding='utf-8')
    return json.loads(r.stdout)


def style_of(raw, data):
    for indent in (2, 1, 4):
        for ascii_ in (False, True):
            for nl in ('\n', ''):
                if json.dumps(data, indent=indent, ensure_ascii=ascii_) + nl == raw:
                    return indent, ascii_, nl
    return 2, False, '\n'


def same(a, b):
    if isinstance(a, bool) or isinstance(b, bool):
        return isinstance(a, bool) and isinstance(b, bool) and a == b
    if isinstance(a, (int, float)) and isinstance(b, (int, float)):
        return float(a) == float(b)
    return a == b


def md_type(json_type, items):
    return 'array<%s>' % items if json_type == 'array' and items else json_type


def update_json(doc, tool):
    """Returns the list of differences and mutates doc in place."""
    diffs = []
    if tool['description'] is not None and doc.get('description') != tool['description']:
        diffs.append('description')
        doc['description'] = tool['description']
    schema = doc.setdefault('inputSchema', {'type': 'object', 'properties': {}})
    props = schema.setdefault('properties', {})
    new_props = {}
    required = []
    for p in tool['parameters']:
        name = p['name']
        existing = props.get(name)
        if existing is None:
            diffs.append('missing parameter ' + name)
            node = {'type': p['json']}
            if p.get('items'):
                node['items'] = {'type': p['items']}
            node['description'] = p['description'] or ''
            if 'default' in p:
                node['default'] = p['default']
            new_props[name] = node
        else:
            node = existing
            if node.get('type') != p['json']:
                diffs.append('type of ' + name)
                node['type'] = p['json']
            if p['description'] is not None and node.get('description') != p['description']:
                diffs.append('description of ' + name)
                node['description'] = p['description']
            if 'default' in p:
                if 'default' not in node or not same(node['default'], p['default']):
                    diffs.append('default of ' + name)
                    node['default'] = p['default']
            elif 'default' in node:
                diffs.append('default of ' + name + ' (the code has none)')
                del node['default']
            new_props[name] = node
        if p['required']:
            required.append(name)
    for name in props:
        if name not in new_props:
            diffs.append('extra parameter ' + name)
    if list(props.keys()) != list(new_props.keys()) and not any(d.startswith(('missing', 'extra')) for d in diffs):
        diffs.append('parameter order')
    schema['properties'] = new_props
    if required or 'required' in schema:
        if schema.get('required') != required:
            diffs.append('required')
        if required:
            schema['required'] = required
        else:
            schema.pop('required', None)
    return diffs


def md_rows(doc, tool, old_rows):
    rows = []
    props = doc['inputSchema']['properties']
    by_name = {p['name']: p for p in tool['parameters']}
    for name, node in props.items():
        p = by_name[name]
        old = old_rows.get(name)
        type_cell = old[0] if old else md_type(node.get('type', 'string'), (node.get('items') or {}).get('type'))
        desc = node.get('description', '')
        if 'default' in node and '(default' not in desc:
            dv = node['default']
            desc += ' (default: %s)' % ('true' if dv is True else 'false' if dv is False else dv)
        rows.append('| `%s` | %s | %s | %s |' % (name, type_cell, 'Yes' if p['required'] else 'No', desc.replace('|', '\\|')))
    return rows


def generated_from(old_json_description, md_text):
    """True when the md paragraph is (a prefix / extension of) the tool description of the old JSON: text that the generator owns; a hand-written paragraph is left alone."""
    o = (old_json_description or '').strip()
    t = md_text.strip()
    return bool(o) and bool(t) and (t == o or o.startswith(t) or t.startswith(o))


def update_md(path, doc, tool, old_json_description):
    text = open(path, encoding='utf-8', newline='').read()
    orig = text
    desc = tool['description']
    # 1. the paragraph under the title and 2. "## Description": replaced when they are the old tool description (stale), a hand-written summary stays
    m = re.match(r'(# [^\n]*\n\n)(.*?)(\n\n## )', text, re.S)
    if m and desc is not None and generated_from(old_json_description, m.group(2)):
        text = text[:m.start(2)] + desc + text[m.end(2):]
    m = re.search(r'(## Description\n\n)(.*?)(\n\n## )', text, re.S)
    if m and desc is not None and generated_from(old_json_description, m.group(2)):
        text = text[:m.start(2)] + desc + text[m.end(2):]
    # 3. the Input Schema table
    m = re.search(r'(## Input Schema\n\n)(.*?)(\n\n## |\n*$)', text, re.S)
    if m:
        table = m.group(2).split('\n')
        old_rows = {}
        for line in table[2:]:
            cells = [c.strip() for c in re.split(r'(?<!\\)\|', line.strip().strip('|'))]
            if len(cells) >= 4:
                old_rows[cells[0].strip('`')] = (cells[1], cells[2], cells[3])
        if tool['parameters'] or old_rows:
            new = table[:2] if len(table) >= 2 and table[0].startswith('| Parameter') else ['| Parameter | Type | Required | Description |', '|-----------|------|----------|-------------|']
            new = new + md_rows(doc, tool, old_rows)
            text = text[:m.start(2)] + '\n'.join(new) + text[m.end(2):]
    return orig, text


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--check', action='store_true')
    ap.add_argument('--dump', help='metadata dump (JSON) written earlier')
    ap.add_argument('--server', action='append', help='alignment | analysis (default: both)')
    ap.add_argument('--tool', action='append')
    ap.add_argument('--all-owned', action='store_true', help='(default) every Alignment tool + the Analysis tools of the B09 classes')
    ap.add_argument('--all', action='store_true', help='every tool of the selected servers, whatever its owner')
    a = ap.parse_args()
    servers = a.server or ['alignment', 'analysis']
    tools = json.load(open(a.dump, encoding='utf-8')) if a.dump else run_dump([s.capitalize() for s in servers])
    changed, stale = 0, []
    for tool in sorted(tools, key=lambda t: (t['server'], t['name'])):
        server = tool['server'].lower()
        if server not in servers or (a.tool and tool['name'] not in a.tool):
            continue
        jpath = os.path.join(TOOLS, server, tool['name'] + '.mcp.json')
        mpath = os.path.join(TOOLS, server, tool['name'] + '.md')
        if not os.path.exists(jpath):
            stale.append((tool['name'], ['no document ' + os.path.relpath(jpath, ROOT)]))
            continue
        raw = open(jpath, encoding='utf-8', newline='').read()
        doc = json.loads(raw)
        if server != 'alignment' and not a.all and not str(doc.get('methodId', '')).startswith(B09_CLASSES):
            continue
        style = style_of(raw, doc)
        old_description = doc.get('description')
        diffs = update_json(doc, tool)
        md_old = md_new = None
        if os.path.exists(mpath):
            md_old, md_new = update_md(mpath, doc, tool, old_description)
            if md_old != md_new:
                diffs.append('md')
        if diffs:
            stale.append((tool['name'], diffs))
            if not a.check:
                indent, ascii_, nl = style
                with open(jpath, 'w', encoding='utf-8', newline='') as f:
                    f.write(json.dumps(doc, indent=indent, ensure_ascii=ascii_) + nl)
                if md_old is not None and md_old != md_new:
                    with open(mpath, 'w', encoding='utf-8', newline='') as f:
                        f.write(md_new)
                changed += 1
    for name, diffs in stale:
        print('%s %s: %s' % ('STALE' if a.check else 'updated', name, ', '.join(diffs)))
    if a.check:
        print('%d stale tool document(s)' % len(stale))
        return 1 if stale else 0
    print('%d tool document(s) regenerated' % changed)
    return 0


if __name__ == '__main__':
    sys.exit(main())
