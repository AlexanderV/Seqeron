#!/usr/bin/env python3
"""Keep the MolTools tool docs' `docRef` line anchors on the declaration of the method they name.

Every `docs/mcp/tools/moltools/<tool>.mcp.json` carries `"docRef": "Seqeron.Genomics.MolTools/<File>.cs#L<n>"` and
its .md links the same anchor; `docs/mcp/traceability.md` repeats it. Line numbers rot with every edit of the
source (B08 round 4, H29: 18 of 55 were stale), so this script recomputes them:

  python3 scripts/docs/sync_moltools_docref.py            # report drift, exit 1 when there is any
  python3 scripts/docs/sync_moltools_docref.py --fix      # rewrite json, md and traceability.md
  python3 scripts/docs/sync_moltools_docref.py --fix --files CrisprDesigner RestrictionAnalyzer

After `--fix` run `python3 scripts/skills/gen-catalog.py` (the skills catalog reads the md link).

Which declaration: the one of the method named by `methodId` (`Class.Method`) in the docRef file (or, when it has
moved, another file whose name starts with the same stem) with `public static`; with overloads, the one whose
parameter names overlap most with the tool's input properties (snake_case -> camelCase), ties to the one nearest
the old anchor. The test `B07ToolDocsContractTests.DocRef_LineAnchor_PointsAtMethodDeclaration` enforces the same
rule for the files listed in its `AuditedSourceStems`.
"""
import argparse
import glob
import json
import os
import re
import sys

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
SRC = os.path.join(ROOT, 'src', 'Seqeron', 'Algorithms', 'Seqeron.Genomics.MolTools')
DOCS = os.path.join(ROOT, 'docs', 'mcp', 'tools', 'moltools')
TRACE = os.path.join(ROOT, 'docs', 'mcp', 'traceability.md')
ANCHOR = re.compile(r'^Seqeron\.Genomics\.MolTools/(?P<file>[\w.]+\.cs)#L(?P<line>\d+)\b')


def camel(snake):
    head, *rest = snake.split('_')
    return head + ''.join(p[:1].upper() + p[1:] for p in rest)


def parameters(lines, index, method):
    """Parameter names of the declaration starting on lines[index] (may span lines)."""
    text = '\n'.join(lines[index:index + 40])
    start = text.index(method + '(') + len(method) + 1
    depth, angle, current, names = 1, 0, '', []
    for ch in text[start:]:
        if ch == '(':
            depth += 1
        elif ch == ')':
            depth -= 1
            if depth == 0:
                break
        elif ch == '<':
            angle += 1
        elif ch == '>':
            angle -= 1
        if ch == ',' and depth == 1 and angle == 0:
            names.append(current)
            current = ''
        else:
            current += ch
    names.append(current)
    out = []
    for segment in names:
        segment = segment.split('=')[0].strip()
        if segment:
            out.append(re.split(r'\s+', segment)[-1])
    return out


def declarations(path, method):
    lines = open(path, encoding='utf-8').read().split('\n')
    strict, lenient = [], []
    name = re.escape(method)
    for i, line in enumerate(lines):
        if line.lstrip().startswith('//') or not re.search(r'\b' + name + r'\(', line):
            continue
        # a declaration: `public static <return type> Method(` - no call, no `=>` before the name
        if re.search(r'public static\s+[^=(]*\b' + name + r'\(', line):
            strict.append((i + 1, parameters(lines, i, method)))
        elif 'public static' in line or (i > 0 and 'public static' in lines[i - 1]):
            lenient.append((i + 1, parameters(lines, i, method)))
    return strict or lenient


def sync_secondary(md, primary_file, primary_line):
    """Other `File.cs#Ln` links in the md (e.g. 'Overhang (`GetOverhangSequence`): [..#L1708](..#L1708)') follow the
    method named in backticks on the same line; a line without a resolvable name is reported, not guessed."""
    notes, out = [], []
    primary = re.compile(r'Source: ')
    for line in md.split('\n'):
        anchors = list(re.finditer(r'([\w.]+\.cs)#L(\d+)', line))
        if not anchors or primary.search(line):
            out.append(line)
            continue
        names = re.findall(r'`(\w+)(?:\([^`]*\))?`', line)
        for m in anchors:
            f, n = m[1], int(m[2])
            path = os.path.join(SRC, f)
            if not os.path.exists(path):
                continue
            target = None
            for name in names:
                decls = declarations(path, name)
                if decls:
                    keep = n in [d[0] for d in decls]
                    target = n if keep else min(decls, key=lambda d: abs(d[0] - n))[0]
                    break
            if target is None:
                notes.append(f'UNRESOLVED secondary anchor {f}#L{n}: name the method in backticks on that line')
            elif target != n:
                notes.append(f'secondary {f}#L{n} -> L{target}')
                line = re.sub(re.escape(f'{f}#L{n}') + r'(?!\d)', f'{f}#L{target}', line)
        out.append(line)
    return '\n'.join(out), notes


def best(decls, tool_params, old_line):
    def score(d):
        overlap = len(set(d[1]) & tool_params)
        return (-overlap, len(set(d[1]) - tool_params), abs(d[0] - old_line))
    return min(decls, key=score)[0]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--fix', action='store_true')
    ap.add_argument('--files', nargs='*', default=None, help='only docRefs whose file name starts with one of these stems')
    args = ap.parse_args()

    trace = open(TRACE, encoding='utf-8').read()
    drift = 0
    for jpath in sorted(glob.glob(os.path.join(DOCS, '*.mcp.json'))):
        text = open(jpath, encoding='utf-8').read()
        doc = json.loads(text)
        tool, method = doc['toolName'], doc['methodId'].split('.')[-1]
        m = ANCHOR.match(doc['docRef'])
        if not m:
            print(f'UNPARSABLE {tool}: {doc["docRef"]}')
            drift += 1
            continue
        file, old = m['file'], int(m['line'])
        if args.files is not None and not any(file.startswith(s) for s in args.files):
            continue
        stem = file.split('.')[0]
        tool_params = {camel(p) for p in doc.get('inputSchema', {}).get('properties', {})}
        candidates = [file] + sorted(os.path.basename(p) for p in glob.glob(os.path.join(SRC, stem + '*.cs'))
                                     if os.path.basename(p) != file)
        target = None
        for cand in candidates:
            decls = declarations(os.path.join(SRC, cand), method)
            if decls:
                # an anchor that already sits on a declaration of the method stays (a manual overload choice is kept)
                keep = cand == file and old in [d[0] for d in decls]
                target = (cand, old if keep else best(decls, tool_params, old))
                break
        if target is None:
            print(f'NOT FOUND {tool}: {method} in {file}')
            drift += 1
            continue
        new_file, new = target
        mpath = jpath[:-len('.mcp.json')] + '.md'
        md = open(mpath, encoding='utf-8').read()
        want = f'{new_file}#L{new}'
        # the md `Source: [File.cs#Ln](...#Ln)` link and the traceability row must carry the same anchor as the json
        md_link = re.compile(r'(Source: \[)[\w.]+\.cs#L\d+(\]\([^)]*/)[\w.]+\.cs#L\d+(\))')
        trace_row = re.compile(r'(?m)^(\| ' + re.escape(tool) + r' \|.*?\| )[\w.]+\.cs#L\d+( \|)')
        md2 = md_link.sub(lambda t: f'{t[1]}{want}{t[2]}{want}{t[3]}', md, count=1)
        if (new_file, new) != (file, old):
            md2 = re.sub(re.escape(f'{file}#L{old}') + r'(?!\d)', want, md2)
        trace2 = trace_row.sub(lambda t: f'{t[1]}{want}{t[2]}', trace, count=1)
        md2, secondary = sync_secondary(md2, new_file, new)
        if (new_file, new) == (file, old) and md2 == md and trace2 == trace and not secondary:
            continue
        drift += 1
        print(f'{tool}: {file}#L{old} -> {want}' + ('' if (new_file, new) != (file, old) else '  (md / traceability only)'))
        for note in secondary:
            print('    ' + note)
        if not args.fix:
            continue
        if (new_file, new) != (file, old):
            text = text.replace(f'{file}#L{old}"', f'{want}"', 1)
            open(jpath, 'w', encoding='utf-8', newline='').write(text)
        open(mpath, 'w', encoding='utf-8', newline='').write(md2)
        trace = trace2
    if args.fix:
        open(TRACE, 'w', encoding='utf-8', newline='').write(trace)
    print(f'{drift} anchor(s) {"rewritten" if args.fix else "drifted"}')
    return 1 if drift and not args.fix else 0


if __name__ == '__main__':
    sys.exit(main())
