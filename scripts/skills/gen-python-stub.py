#!/usr/bin/env python3
"""Offline codegen: turn a Seqeron MCP tool schema into a typed Python wrapper.

Reads a tool's machine-readable ``.mcp.json`` schema (located via the generated
catalog) and emits a typed Python function that calls the Seqeron MCP client
helper::

    from seqeron_mcp_client import call_tool
    call_tool("<ServerName>", "<tool>", { ...args... }, project_root=project_root)

Pure Python stdlib. Runs fully OFFLINE — no dotnet, no ``mcp`` package; it only
parses JSON schema files.

Usage:
    python3 scripts/skills/gen-python-stub.py <tool> [<tool> ...] \\
        [--out FILE] [--project-root PATH]
"""

from __future__ import annotations

import argparse
import json
import keyword
import os
import re
import sys

# Repo root = two levels up from this file (scripts/skills/gen-python-stub.py).
_THIS_DIR = os.path.dirname(os.path.abspath(__file__))
_REPO_ROOT = os.path.abspath(os.path.join(_THIS_DIR, os.pardir, os.pardir))

CATALOG_REL = os.path.join("docs", "skills", "_generated", "catalog.json")

# JSON-Schema primitive type -> Python type hint.
_TYPE_MAP = {
    "string": "str",
    "integer": "int",
    "number": "float",
    "boolean": "bool",
    "array": "list",
    "object": "dict",
    "null": "None",
}


def _fail(msg):
    sys.stderr.write("gen-python-stub.py: error: " + msg + "\n")
    sys.exit(1)


def _load_catalog(project_root):
    path = os.path.join(project_root, CATALOG_REL)
    if not os.path.isfile(path):
        _fail("catalog not found: " + path)
    with open(path, "r", encoding="utf-8") as fh:
        return json.load(fh)


def _schema_path_for(entry, project_root):
    """Resolve a catalog entry's .mcp.json path (docPath with .md -> .mcp.json)."""
    doc_path = entry.get("docPath")
    if not doc_path:
        _fail("catalog entry for tool '%s' has no docPath" % entry.get("tool"))
    if doc_path.endswith(".md"):
        mcp_rel = doc_path[: -len(".md")] + ".mcp.json"
    else:
        mcp_rel = doc_path + ".mcp.json"
    return os.path.join(project_root, mcp_rel)


def _load_schema(tool, catalog, project_root):
    entry = None
    for item in catalog:
        if item.get("tool") == tool:
            entry = item
            break
    if entry is None:
        known = ", ".join(sorted(i.get("tool", "") for i in catalog)[:10])
        _fail(
            "unknown tool '%s' (not in catalog). Examples of known tools: %s ..."
            % (tool, known)
        )
    path = _schema_path_for(entry, project_root)
    if not os.path.isfile(path):
        _fail("schema file not found for tool '%s': %s" % (tool, path))
    with open(path, "r", encoding="utf-8") as fh:
        return json.load(fh)


def _sanitize_identifier(name):
    """Turn an arbitrary tool/param name into a valid Python identifier."""
    ident = re.sub(r"\W", "_", str(name))
    if not ident:
        ident = "_"
    if ident[0].isdigit():
        ident = "_" + ident
    if keyword.iskeyword(ident):
        ident = ident + "_"
    return ident


def _py_type(prop_schema):
    """Map a JSON-Schema property to a Python type-hint string."""
    t = prop_schema.get("type")
    if isinstance(t, list):
        # Union type: pick the first non-null primitive, else fall back.
        for candidate in t:
            if candidate != "null" and candidate in _TYPE_MAP:
                return _TYPE_MAP[candidate]
        return "object"
    if t in _TYPE_MAP:
        return _TYPE_MAP[t]
    return "dict"


def _default_literal(prop_schema):
    """Return a Python literal for a schema 'default', or None sentinel."""
    if "default" not in prop_schema:
        return "None"
    return repr(prop_schema["default"])


def _wrap(text, width, indent):
    """Simple word-wrap for docstring lines (stdlib textwrap avoided for determinism)."""
    words = str(text).split()
    lines = []
    cur = ""
    for w in words:
        if not cur:
            cur = w
        elif len(cur) + 1 + len(w) <= width:
            cur += " " + w
        else:
            lines.append(cur)
            cur = w
    if cur:
        lines.append(cur)
    if not lines:
        lines = [""]
    return [indent + ln for ln in lines]


def _build_docstring(schema, props, required):
    lines = []
    desc = schema.get("description")
    if desc:
        lines.extend(_wrap(desc, 76, ""))
    else:
        lines.append(schema.get("toolName", "MCP tool wrapper."))

    # Parameters.
    if props:
        lines.append("")
        lines.append("Args:")
        for name in props:  # props already ordered
            pschema = props[name]
            pdesc = pschema.get("description", "")
            req = "required" if name in required else "optional"
            head = "    %s (%s, %s):" % (
                _sanitize_identifier(name),
                _py_type(pschema),
                req,
            )
            if pdesc:
                wrapped = _wrap(pdesc, 68, "        ")
                lines.append(head)
                lines.extend(wrapped)
            else:
                lines.append(head)
        lines.append("    project_root (str, optional):")
        lines.append("        Path to the Seqeron repository root. Defaults to '.'.")

    # Output fields.
    out = schema.get("outputSchema") or {}
    out_props = out.get("properties") or {}
    if out_props:
        lines.append("")
        lines.append("Returns:")
        lines.append("    dict: Structured tool output with fields:")
        for name in out_props:
            fschema = out_props[name]
            fdesc = fschema.get("description", "")
            head = "        %s (%s):" % (name, _py_type(fschema))
            if fdesc:
                lines.append(head)
                lines.extend(_wrap(fdesc, 60, "            "))
            else:
                lines.append(head)

    # One example.
    examples = schema.get("examples") or []
    if examples:
        ex = examples[0]
        lines.append("")
        lines.append("Example:")
        ex_name = ex.get("name")
        if ex_name:
            lines.extend(_wrap(ex_name, 72, "    "))
        if "input" in ex:
            lines.append("    input:  " + json.dumps(ex["input"], sort_keys=True))
        if "output" in ex:
            lines.append("    output: " + json.dumps(ex["output"], sort_keys=True))

    return lines


def generate_function(schema):
    tool = schema.get("toolName")
    server = schema.get("serverName")
    if not tool or not server:
        _fail("schema missing toolName/serverName")

    func_name = _sanitize_identifier(tool)

    in_schema = schema.get("inputSchema") or {}
    props = in_schema.get("properties") or {}
    # Preserve insertion order from JSON (Python dicts keep it).
    required = list(in_schema.get("required") or [])

    # Order: required params first (positional), then optional (with defaults).
    ordered_names = [n for n in props if n in required] + [
        n for n in props if n not in required
    ]

    # Build the signature.
    sig_parts = []
    for name in ordered_names:
        pschema = props[name]
        ident = _sanitize_identifier(name)
        hint = _py_type(pschema)
        if name in required:
            sig_parts.append("%s: %s" % (ident, hint))
        else:
            sig_parts.append(
                "%s: %s = %s" % (ident, hint, _default_literal(pschema))
            )
    # project_root is always keyword-only.
    sig_parts.append("*")
    sig_parts.append('project_root: str = "."')

    lines = []
    if len(sig_parts) <= 3 and all("\n" not in p for p in sig_parts):
        # Compact single-line signature when short.
        joined = ", ".join(sig_parts)
        if len("def %s(%s):" % (func_name, joined)) <= 88:
            lines.append("def %s(%s):" % (func_name, joined))
        else:
            lines.append("def %s(" % func_name)
            for p in sig_parts:
                lines.append("    %s," % p)
            lines.append("):")
    else:
        lines.append("def %s(" % func_name)
        for p in sig_parts:
            lines.append("    %s," % p)
        lines.append("):")

    # Docstring.
    doc_lines = _build_docstring(schema, props, required)
    lines.append('    """' + (doc_lines[0] if doc_lines else ""))
    for dl in doc_lines[1:]:
        lines.append("    " + dl if dl else "")
    lines.append('    """')

    # Body: build arguments dict.
    lines.append("    arguments = {}")
    for name in ordered_names:
        ident = _sanitize_identifier(name)
        if name in required:
            lines.append("    arguments[%r] = %s" % (name, ident))
        else:
            lines.append("    if %s is not None:" % ident)
            lines.append("        arguments[%r] = %s" % (name, ident))
    lines.append(
        '    return call_tool(%r, %r, arguments, project_root=project_root)'
        % (server, tool)
    )

    return "\n".join(lines)


def generate_module(tools, project_root):
    catalog = _load_catalog(project_root)
    header = [
        "# generated by gen-python-stub.py — do not edit by hand",
        "from seqeron_mcp_client import call_tool",
        "",
        "",
    ]
    funcs = []
    for tool in tools:
        schema = _load_schema(tool, catalog, project_root)
        funcs.append(generate_function(schema))
    body = "\n\n\n".join(funcs)
    return "\n".join(header) + body + "\n"


def main(argv=None):
    parser = argparse.ArgumentParser(
        description="Generate typed Python MCP wrapper stubs from Seqeron tool schemas."
    )
    parser.add_argument("tools", nargs="+", help="tool name(s) to generate")
    parser.add_argument("--out", default=None, help="output file (default: stdout)")
    parser.add_argument(
        "--project-root",
        default=_REPO_ROOT,
        help="path to the Seqeron repo root (default: inferred from script location)",
    )
    args = parser.parse_args(argv)

    module_text = generate_module(args.tools, args.project_root)

    if args.out:
        with open(args.out, "w", encoding="utf-8") as fh:
            fh.write(module_text)
    else:
        sys.stdout.write(module_text)
    return 0


if __name__ == "__main__":
    sys.exit(main())
