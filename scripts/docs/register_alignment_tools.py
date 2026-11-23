#!/usr/bin/env python3
"""Registers new Alignment-server MCP tools in the hand-kept registries (B09 WP67, F307+): creates their documents (docs/mcp/tools/alignment/<tool>.md / .mcp.json; the input-schema rows are filled by sync_mcp_tool_docs.py from the code),
MCP_STATUS.md (counts + rows), the server README (count + rows), traceability.md, scripts/skills/gen-catalog.py EXPECTED_COUNTS / EXPECTED_TOTAL, the quoted tool totals of the README / skill texts (old total -> new total); then run
`python3 scripts/docs/sync_mcp_tool_docs.py`, `python3 scripts/skills/gen-catalog.py`, `python3 scripts/skills/sync-github-skills.py` (catalogue + mirror).
Usage: python3 scripts/docs/register_alignment_tools.py spec.json
spec.json: [{"name", "short", "methodId", "sourceFile", "sourcePath", "algorithmDoc", "oracleScript", "description", "output": {field: [type, description]}, "errors": [...], "example": {title, prompt, input, output, note}, "after": "<existing tool the README / traceability rows follow>"}]. Idempotent."""
import json, os, re, sys

spec = json.load(open(sys.argv[1]))
def rd(p): return open(p).read()
def wr(p, s): open(p, 'w').write(s)

ROOT = 'docs/mcp/tools/alignment'
for t in spec:
    name = t['name']
    md, js = '%s/%s.md' % (ROOT, name), '%s/%s.mcp.json' % (ROOT, name)
    if os.path.exists(md) and os.path.exists(js):
        continue
    out_rows = ''.join('| `%s` | %s | %s |\n' % (k, v[0], v[1]) for k, v in t['output'].items())
    err_rows = ''.join('| %d | %s |\n' % (1001 + i, e) for i, e in enumerate(t['errors']))
    ex = t['example']
    text = f"""# {name}

{t.get('summary', t['short'])}

## Overview

| Property | Value |
|----------|-------|
| **Server** | Alignment |
| **Tool Name** | `{name}` |
| **Method ID** | `{t['methodId']}` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

{t['description']}

## Core Documentation Reference

- Source: [{t['sourceFile']}](../../../../src/Seqeron/Algorithms/{t['sourcePath']})
- Algorithm: [{t['algorithmDoc'].split('/')[-1]}](../../../algorithms/{t['algorithmDoc']})
- Oracle: [{t['oracleScript']}](../../../../scripts/alignment/{t['oracleScript']})

## Input Schema

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|

## Output Schema

| Field | Type | Description |
|-------|------|-------------|
{out_rows}
## Errors

| Code | Message |
|------|---------|
{err_rows}
## Examples

### Example 1: {ex['title']}

**User Prompt:**
> {ex['prompt']}

**Expected Tool Call:**
```json
{json.dumps({'tool': name, 'arguments': ex['input']})}
```

**Response:**
```json
{json.dumps(ex['output'])}
```

{ex['note']}

## See Also

- [pairwise_aligner](pairwise_aligner.md)
- [alignment_operation](alignment_operation.md)
"""
    wr(md, text)
    doc = {
        '$schema': 'https://json-schema.org/draft/2020-12/schema',
        'toolName': name, 'serverName': 'Alignment', 'version': '1.0.0', 'stability': 'stable', 'methodId': t['methodId'],
        'docRef': t['sourcePath'],
        'description': t['description'],
        'inputSchema': {'type': 'object', 'properties': {}, 'required': []},
        'outputSchema': {'type': 'object', 'properties': {k: {'type': v[0], 'description': v[1]} for k, v in t['output'].items()}, 'required': list(t['output'].keys())},
        'errors': [{'code': 1001 + i, 'message': e} for i, e in enumerate(t['errors'])],
        'examples': [{'name': ex['title'], 'input': ex['input'], 'output': ex['output']}],
    }
    wr(js, json.dumps(doc, indent=2) + '\n')
    print('created', name)

status = rd('docs/mcp/MCP_STATUS.md')
readme = rd('src/Seqeron/Mcp/Seqeron.Mcp.Alignment/README.md')
trace = rd('docs/mcp/traceability.md')
cat = rd('scripts/skills/gen-catalog.py')
contract = rd('tests/Seqeron/Seqeron.Mcp.Alignment.Tests/Contracts/ToolRegistrationContractTests.cs')

new = [t for t in spec if '| `%s` |' % t['name'] not in status]
if not new:
    print('nothing to register'); sys.exit(0)
total_old = int(re.search(r'EXPECTED_TOTAL = (\d+)', cat).group(1))
al_old = int(re.search(r'"alignment": (\d+)', cat).group(1))
k = len(new)
total_new, al_new = total_old + k, al_old + k
status = status.replace('%d MCP tools across 11 servers' % total_old, '%d MCP tools across 11 servers' % total_new)
for lbl in ('Gold-standard binding', 'Schema+Binding tests', 'Docs (.md + .mcp.json)'):
    status = status.replace('| %s | %d | %d |' % (lbl, total_old, total_old), '| %s | %d | %d |' % (lbl, total_new, total_new))
status = status.replace('| Alignment | Seqeron.Mcp.Alignment | %d | %d/%d | %d/%d | %d/%d | ✅ done |' % ((al_old,) * 7), '| Alignment | Seqeron.Mcp.Alignment | %d | %d/%d | %d/%d | %d/%d | ✅ done |' % ((al_new,) * 7))
status = status.replace('### Alignment (Seqeron.Mcp.Alignment) — %d tools' % al_old, '### Alignment (Seqeron.Mcp.Alignment) — %d tools' % al_new)
m = re.search(r'(^\| (\d+) \| `[^`]+` \| ☑ \| ☑ \| ☑ \|\n)(?=\n?### Analysis)', status, re.M)
num = int(m.group(2))
rows = ''.join('| %d | `%s` | ☑ | ☑ | ☑ |\n' % (num + i + 1, t['name']) for i, t in enumerate(new))
status = status[:m.end()] + rows + status[m.end():]
wr('docs/mcp/MCP_STATUS.md', status)
readme = readme.replace('Exposes **%d tools**' % al_old, 'Exposes **%d tools**' % al_new)
for t in new:
    anchor = re.search(r'^\| `%s` \|[^\n]*\n' % re.escape(t['after']), readme, re.M)
    readme = readme[:anchor.end()] + '| `%s` | %s |\n' % (t['name'], t['short']) + readme[anchor.end():]
    anchor = re.search(r'^\| %s \| Alignment \|[^\n]*\n' % re.escape(t['after']), trace, re.M)
    trace = trace[:anchor.end()] + '| %s | Alignment | %s | Alignment/%s | 0/0 | – | – | Build-Ready |\n' % (t['name'], t['methodId'], t['sourceFile']) + trace[anchor.end():]
    if '        "%s",\n' % t['name'] not in contract:
        contract = contract.replace('        "%s",\n' % t['after'], '        "%s",\n        "%s",\n' % (t['after'], t['name']), 1)
wr('src/Seqeron/Mcp/Seqeron.Mcp.Alignment/README.md', readme)
wr('docs/mcp/traceability.md', trace)
wr('tests/Seqeron/Seqeron.Mcp.Alignment.Tests/Contracts/ToolRegistrationContractTests.cs', contract)
cat = cat.replace('"alignment": %d,' % al_old, '"alignment": %d,' % al_new).replace('EXPECTED_TOTAL = %d' % total_old, 'EXPECTED_TOTAL = %d' % total_new)
cat = re.sub(r'\((\d+) documented;', '(%d documented;' % total_new, cat)
wr('scripts/skills/gen-catalog.py', cat)
# the quoted totals of the README / skill texts
for p in ['README.md', 'docs/mcp/README.md', '.claude/skills/README.md', '.claude/skills/marketplace.json', '.claude/skills/seqeron-dev/reference/api-conventions.md', '.claude/skills/seqeron-discovery/SKILL.md', '.claude/skills/seqeron-setup/SKILL.md', '.github/skills/README.md', '.github/copilot-instructions.md']:
    s = rd(p)
    s2 = re.sub(r'(?<![\d.])%d(?= tools| MCP tool| tool| schemas|%%20tools| strict-schema|\)|\*\* tools)' % total_old, str(total_new), s)
    s2 = s2.replace('= %d tools' % total_old, '= %d tools' % total_new).replace('among %d' % total_old, 'among %d' % total_new).replace('MCP-%d%%20tools' % total_old, 'MCP-%d%%20tools' % total_new)
    if s2 != s:
        wr(p, s2)
print('registered', [t['name'] for t in new], 'alignment', al_new, 'total', total_new)
