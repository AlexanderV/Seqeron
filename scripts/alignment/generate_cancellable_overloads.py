#!/usr/bin/env python3
"""Generates src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/SequenceAligner.Cancellable.cs (audit R14-15, B09 F294).

Every public pairwise aligner of `SequenceAligner` (a static method returning AlignmentResult or IReadOnlyList<AlignmentResult> whose first parameters are the two sequences:
GlobalAlign / LocalAlign / SemiGlobalAlign, every Affine, LinearSpace, GeneralGap, CoOptimal, TopN, TopNLinearSpace, DualAffine, NeedleAlign, WaterAlign, StretcherAlign) gets an overload with the SAME parameters
(defaults dropped) and a trailing `CancellationToken cancellationToken` (required: an optional one would overlap the original, which the S3427 analyzer rejects; the existing
GlobalAlign token overloads are the same shape). The overload enters the ambient cancellation scope (`AlignmentCancellation.Enter`, the one mechanism the dynamic programs
poll) and calls the original. A call that omits the token still binds to the original, so nothing existing changes.
Methods that already have a token overload (the linear GlobalAlign string / DnaSequence overloads) are skipped.

  python3 scripts/alignment/generate_cancellable_overloads.py            write the file
  python3 scripts/alignment/generate_cancellable_overloads.py --check    exit 1 when the committed file is stale
"""
import glob
import os
import re
import sys

ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir, os.pardir))
DIR = os.path.join(ROOT, "src", "Seqeron", "Algorithms", "Seqeron.Genomics.Alignment")
OUT = os.path.join(DIR, "SequenceAligner.Cancellable.cs")
RX = re.compile(r"public static\s+(AlignmentResult|IReadOnlyList<AlignmentResult>)\s+(\w+)\(")


def balanced(text, i):
    depth, j = 1, i
    while depth:
        c = text[j]
        if c == "(":
            depth += 1
        elif c == ")":
            depth -= 1
        j += 1
    return text[i:j - 1]


def split_params(s):
    s = re.sub(r"//[^\n]*", "", s)
    out, depth, cur = [], 0, ""
    for c in s:
        if c in "<(":
            depth += 1
        elif c in ">)":
            depth -= 1
        if c == "," and depth == 0:
            out.append(cur.strip())
            cur = ""
        else:
            cur += c
    if cur.strip():
        out.append(cur.strip())
    return out


def parse(p):
    p = " ".join(p.split())
    head = p.split(" = ")[0] if " = " in p else p
    typ, name = head.rsplit(" ", 1)
    return typ.strip(), name.strip()


def collect():
    methods = []
    for f in sorted(glob.glob(os.path.join(DIR, "SequenceAligner*.cs"))):
        if os.path.abspath(f) == OUT:
            continue
        text = open(f, encoding="utf-8").read()
        for m in RX.finditer(text):
            params = [parse(p) for p in split_params(balanced(text, m.end()))]
            if not params or params[0][1] != "sequence1":
                continue
            # attributes directly above the declaration
            line_start = text.rfind("\n", 0, m.start()) + 1
            above = text[:line_start].rstrip("\n").split("\n")
            attrs = []
            while above and above[-1].strip().startswith("["):
                attrs.insert(0, above.pop().strip())
            methods.append({"ret": m.group(1), "name": m.group(2), "params": params, "attrs": attrs, "file": os.path.basename(f)})
    return methods


def cref(name, params):
    return "%s(%s)" % (name, ", ".join(t.replace("<", "{").replace(">", "}") for t, _ in params))


def render():
    methods = collect()
    have = {(m["name"], tuple(t for t, _ in m["params"])) for m in methods}
    out = []
    seen = set()
    for m in methods:
        types = tuple(t for t, _ in m["params"])
        if any(t == "CancellationToken" for t in types):
            continue
        # an existing overload that takes the token right after these parameters (GlobalAlign(string, string, ScoringMatrix?, CancellationToken, IProgress<double>? = null))
        if any(n == m["name"] and ts[:len(types)] == types and len(ts) > len(types) and ts[len(types)] == "CancellationToken" for n, ts in have):
            continue
        key = (m["name"], types)
        if key in seen:
            continue
        seen.add(key)
        decl = ", ".join("%s %s" % (t, n) for t, n in m["params"])
        args = ", ".join(n for _, n in m["params"])
        c = cref(m["name"], m["params"])
        out.append("    /// <summary><see cref=\"%s\"/> that observes <paramref name=\"cancellationToken\"/> (audit R14-15): the dynamic programs poll it at their row / iteration boundaries and the call throws <see cref=\"OperationCanceledException\"/> without a partial result.</summary>" % c)
        out.append("    /// <inheritdoc cref=\"%s\"/>" % c)
        out.append("    /// <param name=\"cancellationToken\">The token to observe; <c>CancellationToken.None</c> = never cancelled (the ambient scope of an outer call still applies).</param>")
        attrs = m["attrs"]
        for a in attrs:
            out.append("    " + a)
        out.append("    public static %s %s(%s, CancellationToken cancellationToken)" % (m["ret"], m["name"], decl))
        out.append("    {")
        out.append("        using var cancelScope = AlignmentCancellation.Enter(cancellationToken);")
        out.append("        return %s(%s);" % (m["name"], args))
        out.append("    }")
        out.append("")
    head = [
        "using Seqeron.Genomics.Core;",
        "using Seqeron.Genomics.Infrastructure;",
        "",
        "namespace Seqeron.Genomics.Alignment;",
        "",
        "// Generated by scripts/alignment/generate_cancellable_overloads.py (audit R14-15, B09 F294): one overload with a trailing CancellationToken for every public pairwise aligner.",
        "// Do not edit by hand; run the script. The overload enters the ambient cancellation scope (SequenceAligner.Cancellation.cs) and calls the original, which polls it at its row / iteration boundaries.",
        "public static partial class SequenceAligner",
        "{",
    ]
    return "\n".join(head + out).rstrip("\n") + "\n}\n"


def main(argv):
    text = render()
    if "--check" in argv:
        current = open(OUT, encoding="utf-8").read() if os.path.isfile(OUT) else ""
        if current != text:
            print("SequenceAligner.Cancellable.cs is stale: run python3 scripts/alignment/generate_cancellable_overloads.py")
            return 1
        print("SequenceAligner.Cancellable.cs is fresh (%d overloads)" % text.count("public static"))
        return 0
    with open(OUT, "w", encoding="utf-8", newline="\n") as f:
        f.write(text)
    print("wrote %s (%d overloads)" % (OUT, text.count("public static")))
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
