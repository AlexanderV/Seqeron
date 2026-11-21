#!/usr/bin/env python3
"""Registers new Analysis-server MCP tools in the hand-kept registries (B09 A3-9b / A2-28d): MCP_STATUS.md (counts + row), the server README, traceability.md, the bio-annotation tool-map, scripts/skills/gen-catalog.py
EXPECTED_COUNTS / EXPECTED_TOTAL and the contract list of the tool-registration test; then `python3 scripts/skills/gen-catalog.py` regenerates the catalog.  Usage: python3 scripts/docs/register_mcp_tools.py spec.json
spec.json: [{"name", "short", "methodId", "file", "after": "<existing tool the new rows follow>"}]. Idempotent: a tool already present is skipped."""
import json, re, sys

spec = json.load(open(sys.argv[1]))
def rd(p): return open(p).read()
def wr(p, s): open(p, 'w').write(s)

status = rd('docs/mcp/MCP_STATUS.md')
readme = rd('src/Seqeron/Mcp/Seqeron.Mcp.Analysis/README.md')
trace = rd('docs/mcp/traceability.md')
tmap = rd('.claude/skills/bio-annotation/reference/tool-map.md')
cat = rd('scripts/skills/gen-catalog.py')
contract = rd('tests/Seqeron/Seqeron.Mcp.Analysis.Tests/Contracts/ToolRegistrationContractTests.cs')

new = [t for t in spec if '| `%s` |' % t['name'] not in status]
if not new:
    print('nothing to register'); sys.exit(0)
total_old = int(re.search(r'EXPECTED_TOTAL = (\d+)', cat).group(1))
an_old = int(re.search(r'"analysis": (\d+)', cat).group(1))
k = len(new)
total_new, an_new = total_old + k, an_old + k
status = status.replace('%d MCP tools across 11 servers' % total_old, '%d MCP tools across 11 servers' % total_new)
for lbl in ('Gold-standard binding', 'Schema+Binding tests', 'Docs (.md + .mcp.json)'):
    status = status.replace('| %s | %d | %d |' % (lbl, total_old, total_old), '| %s | %d | %d |' % (lbl, total_new, total_new))
status = status.replace('| Analysis | Seqeron.Mcp.Analysis | %d | %d/%d | %d/%d | %d/%d | ✅ done |' % ((an_old,) * 7), '| Analysis | Seqeron.Mcp.Analysis | %d | %d/%d | %d/%d | %d/%d | ✅ done |' % ((an_new,) * 7))
status = status.replace('### Analysis (Seqeron.Mcp.Analysis) — %d tools' % an_old, '### Analysis (Seqeron.Mcp.Analysis) — %d tools' % an_new)
last = re.findall(r'^\| (\d+) \|[^\n]*\n(?=\n?### Annotation)', status, re.M)
m = re.search(r'(^\| (\d+) \| `[^`]+` \| ☑ \| ☑ \| ☑ \|\n)(?=\n?### Annotation)', status, re.M)
num = int(m.group(2))
rows = ''.join('| %d | `%s` | ☑ | ☑ | ☑ |\n' % (num + i + 1, t['name']) for i, t in enumerate(new))
status = status[:m.end()] + rows + status[m.end():]
wr('docs/mcp/MCP_STATUS.md', status)

readme = readme.replace('Exposes **%d tools**' % an_old, 'Exposes **%d tools**' % an_new)
for t in new:
    anchor = re.search(r'^\| `%s` \|[^\n]*\n' % re.escape(t['after']), readme, re.M)
    readme = readme[:anchor.end()] + '| `%s` | %s |\n' % (t['name'], t['short']) + readme[anchor.end():]
    anchor = re.search(r'^\| %s \|[^\n]*\n' % re.escape(t['after']), trace, re.M)
    trace = trace[:anchor.end()] + '| %s | Analysis | %s | %s | 0/0 | – | – | Build-Ready |\n' % (t['name'], t['methodId'], t['file']) + trace[anchor.end():]
    anchor = re.search(r'^\| `%s` \| X \|[^\n]*\n' % re.escape(t['after']), tmap, re.M) or re.search(r'^\| `%s` \| X \|[^\n]*\n' % re.escape(t.get('mapAfter', t['after'])), tmap, re.M)
    if anchor is not None:      # the tool-map lists only some tools (the Mash tools are not in it): a tool whose anchor is absent is left out of it
        tmap = tmap[:anchor.end()] + '| `%s` | X | %s | `%s` |\n' % (t['name'], t['short'], t['methodId']) + tmap[anchor.end():]
    contract = contract.replace('        "%s",\n' % t['after'], '        "%s",\n        "%s",\n' % (t['after'], t['name']), 1) if '        "%s",\n' % t['name'] not in contract else contract
wr('src/Seqeron/Mcp/Seqeron.Mcp.Analysis/README.md', readme)
wr('docs/mcp/traceability.md', trace)
wr('.claude/skills/bio-annotation/reference/tool-map.md', tmap)
wr('tests/Seqeron/Seqeron.Mcp.Analysis.Tests/Contracts/ToolRegistrationContractTests.cs', contract)
cat = cat.replace('"analysis": %d,' % an_old, '"analysis": %d,' % an_new).replace('EXPECTED_TOTAL = %d' % total_old, 'EXPECTED_TOTAL = %d' % total_new)
cat = re.sub(r'\((\d+) documented;', '(%d documented;' % total_new, cat)
wr('scripts/skills/gen-catalog.py', cat)
print('registered', [t['name'] for t in new], 'analysis', an_new, 'total', total_new)
