#!/usr/bin/env python3
"""Creates the per-tool documents (docs/mcp/tools/analysis/<tool>.md / .mcp.json) of the B09 A3-9b / A2-28d tools (nucmer, delta-filter, show-*, dnadiff, est2genome) from the SPEC below, then
`python3 scripts/docs/sync_mcp_tool_docs.py` fills the input-schema rows from the code. Run from the repository root: python3 scripts/docs/gen_mummer_est_tool_docs.py
Only creates a document that does not exist yet (hand text of an existing one is never overwritten)."""
import json, os, sys

ROOT = 'docs/mcp/tools/analysis'
SPEC = json.load(open(sys.argv[1] if len(sys.argv) > 1 else os.path.join(os.path.dirname(__file__), 'mummer_est_tool_spec.json')))

for t in SPEC:
    name = t['name']
    md = '%s/%s.md' % (ROOT, name)
    js = '%s/%s.mcp.json' % (ROOT, name)
    if os.path.exists(md) and os.path.exists(js):
        continue
    out_rows = ''.join('| `%s` | %s | %s |\n' % (k, v[0], v[1]) for k, v in t['output'].items())
    err_rows = ''.join('| 1001 | %s |\n' % e for e in t['errors'])
    ex = t['example']
    text = f"""# {name}

{t['summary']}

## Overview

| Property | Value |
|----------|-------|
| **Server** | Analysis |
| **Tool Name** | `{name}` |
| **Method ID** | `{t['methodId']}` |
| **Version** | 1.0.0 |
| **Stability** | Stable |

## Description

{t['description']}

## Core Documentation Reference

- Source: [{t['sourceFile']}](../../../../src/Seqeron/Algorithms/{t['sourcePath']})
- Algorithm doc: [{t['algorithmDoc'].split('/')[-1]}](../../../algorithms/{t['algorithmDoc']})
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

- [nucmer_clusters](nucmer_clusters.md) - the anchor-clustering stage of nucmer
"""
    open(md, 'w').write(text)
    doc = {
        '$schema': 'https://json-schema.org/draft/2020-12/schema',
        'toolName': name, 'serverName': 'Analysis', 'version': '1.0.0', 'stability': 'stable', 'methodId': t['methodId'],
        'docRef': t['sourcePath'].split('/', 1)[1] if '/' in t['sourcePath'] else t['sourcePath'],
        'description': t['description'],
        'inputSchema': {'type': 'object', 'properties': {}, 'required': []},
        'outputSchema': {'type': 'object', 'properties': {k: {'type': v[0], 'description': v[1]} for k, v in t['output'].items()}, 'required': list(t['output'].keys())},
        'errors': [{'code': 1001, 'message': e} for e in t['errors']],
        'examples': [{'name': ex['title'], 'input': ex['input'], 'output': ex['output']}],
    }
    json.dump(doc, open(js, 'w'), indent=2)
    open(js, 'a').write('\n')
    print('created', name)
