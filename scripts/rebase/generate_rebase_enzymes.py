#!/usr/bin/env python3
"""Reproducible generator for the vendored REBASE restriction-enzyme table (review 2026-09 B08, D4/F59).

This is a *build-time* developer tool. Its outputs are committed; nothing runs it at build or test time.

Regenerate (from the repo root):

    python3 -I scripts/rebase/generate_rebase_enzymes.py table \
        --biopython-dict <site-packages>/Bio/Restriction/Restriction_Dictionary.py \
        --allenz <crisprBase>/inst/enzymes/allenz.txt

    <venv with biopython==1.88>/bin/python -I scripts/rebase/generate_rebase_enzymes.py oracle

`table` needs only the two pinned input files (parsed with `ast.literal_eval`, nothing is imported or
executed from them); `oracle` imports Biopython 1.88 itself (`pip install biopython==1.88`) and reads the
table written by `table`.

Inputs (pinned; the script refuses to run on any other bytes):

  * Biopython 1.88 `Bio/Restriction/Restriction_Dictionary.py` -- "Used REBASE emboss files version 404
    (2024)", 1088 enzymes. sha256 2a79099295dbad6061ea67a11e053787c591fcb2eb10fc8c0f89ead908dfa02b
    (byte-identical to tag `biopython-188` on GitHub; docs/Evidence/B08-Source-Provenance.md section 2).
    Licence: the file is "part of the Biopython distribution and governed by its license" (Biopython
    License Agreement: permission to use, copy, modify and distribute ... provided that copyright notices
    and the permission notice appear in supporting documentation).
  * REBASE `allenz.104` (Copyright (c) Dr. Richard J. Roberts, 2021) as vendored in the Bioconductor
    package crisprBase (`crisprVerse/crisprBase` `inst/enzymes/allenz.txt`, package licence MIT).
    Redistribution of these two fields is tracked in docs/Validation/review-2026-09/REBASE-permission-request.md.
    sha256 43b12432ed3c1a19e961d019e3a131336e8ee690a4568a68ddf55ed0e3869388
    (B08-Source-Provenance.md section 3a).

What is taken from where (B08 rule): recognition site and cut offsets from Biopython (REBASE 404);
organism (field <3>) and cognate-methylation site (field <6>) from allenz.104, joined by exact enzyme
name. An enzyme absent from allenz.104 gets empty organism and methylation fields (counted in the header).

Which enzymes: every Biopython enzyme with a defined cut (OneCut: fst5 is not None; TwoCuts: scd5 is not
None as well). The NoCut/NotDefined enzymes (fst5 is None -- REBASE knows the site but not where it is
cut) are EXCLUDED: a RestrictionEnzyme record must carry integer cut offsets, and inventing them would
be a false statement. Their count is written to the header.

Mapping (RestrictionEnzyme convention, F54/F56): Forward = fst5, Reverse = size + fst3;
SecondForward = scd5, SecondReverse = size + scd3 (empty for one-cut enzymes).

Outputs:
  table  -> src/Seqeron/Algorithms/Seqeron.Genomics.MolTools/Resources/rebase_enzymes.tsv
  oracle -> tests/Seqeron/Seqeron.Genomics.Tests/TestData/Restriction/biopython_rebase_search.tsv
            One line per vendored enzyme: a fixed pseudo-random A/C/G/T molecule (random.Random(SEED)) with
            the site planted twice -- once as written, once reverse-complemented, ambiguity codes
            resolved to random concrete bases -- and Biopython 1.88 `enzyme.search(seq, linear=True)` and
            `search(seq, linear=False)` (1-based top-strand cuts, as Biopython returns them).
"""
import argparse
import ast
import hashlib
import os
import random
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))
TABLE = os.path.join(REPO, "src", "Seqeron", "Algorithms", "Seqeron.Genomics.MolTools", "Resources",
                     "rebase_enzymes.tsv")
ORACLE = os.path.join(REPO, "tests", "Seqeron", "Seqeron.Genomics.Tests", "TestData", "Restriction",
                      "biopython_rebase_search.tsv")

BIOPYTHON_DICT_SHA256 = "2a79099295dbad6061ea67a11e053787c591fcb2eb10fc8c0f89ead908dfa02b"
ALLENZ_SHA256 = "43b12432ed3c1a19e961d019e3a131336e8ee690a4568a68ddf55ed0e3869388"
SEED = 20261015

IUPAC = {"A": "A", "C": "C", "G": "G", "T": "T", "R": "AG", "Y": "CT", "S": "CG", "W": "AT", "K": "GT",
         "M": "AC", "B": "CGT", "D": "AGT", "H": "ACT", "V": "ACG", "N": "ACGT"}
COMPLEMENT = str.maketrans("ACGTRYSWKMBDHVN", "TGCAYRSWMKVHDBN")


def sha256_of(path):
    with open(path, "rb") as f:
        return hashlib.sha256(f.read()).hexdigest()


def require_pin(path, expected):
    actual = sha256_of(path)
    if actual != expected:
        sys.exit(f"{path}: sha256 {actual} != pinned {expected}")


def read_biopython_dict(path):
    """rest_dict[...] = {...} assignments, evaluated as literals only."""
    tree = ast.parse(open(path, encoding="utf-8").read())
    enzymes = {}
    for node in tree.body:
        if (isinstance(node, ast.Assign) and len(node.targets) == 1
                and isinstance(node.targets[0], ast.Subscript)
                and isinstance(node.targets[0].value, ast.Name)
                and node.targets[0].value.id == "rest_dict"):
            name = ast.literal_eval(node.targets[0].slice)
            enzymes[name] = ast.literal_eval(node.value)
    return enzymes


def read_allenz(path):
    records, current = {}, None
    with open(path, encoding="latin-1") as f:
        for line in f:
            m = re.match(r"^<(\d)>(.*)$", line.rstrip("\n"))
            if not m:
                continue
            key, value = int(m.group(1)), m.group(2)
            if key == 1:
                if value in records:
                    sys.exit(f"duplicate allenz record {value}")
                current = {}
                records[value] = current
            current[key] = value
    return records


def cmd_table(args):
    require_pin(args.biopython_dict, BIOPYTHON_DICT_SHA256)
    require_pin(args.allenz, ALLENZ_SHA256)
    bio = read_biopython_dict(args.biopython_dict)
    allenz = read_allenz(args.allenz)

    rows, undefined, missing = [], 0, 0
    one_cut = two_cut = 0
    for name in sorted(bio, key=lambda n: (n.upper(), n)):
        e = bio[name]
        if e["fst5"] is None:
            undefined += 1
            continue
        site, size = e["site"], e["size"]
        if size != len(site) or any(c not in IUPAC for c in site):
            sys.exit(f"{name}: unexpected site {site!r} / size {size}")
        fwd, rev = e["fst5"], size + e["fst3"]
        if e["scd5"] is not None:
            two_cut += 1
            sfwd, srev = str(e["scd5"]), str(size + e["scd3"])
        else:
            one_cut += 1
            sfwd = srev = ""
        rec = allenz.get(name)
        if rec is None:
            missing += 1
            organism = methylation = ""
        else:
            organism, methylation = rec.get(3, "").strip(), rec.get(6, "").strip()
        for field in (organism, methylation):
            if "\t" in field or "\n" in field:
                sys.exit(f"{name}: tab/newline in a field")
        rows.append((name, site, str(fwd), str(rev), sfwd, srev, organism, methylation))

    os.makedirs(os.path.dirname(TABLE), exist_ok=True)
    with open(TABLE, "w", encoding="utf-8", newline="\n") as out:
        out.write("# REBASE restriction enzymes with a defined cut, vendored for Seqeron.Genomics.MolTools.RestrictionAnalyzer.\n")
        out.write("# GENERATED by scripts/rebase/generate_rebase_enzymes.py -- do not edit by hand; see that script's header.\n")
        out.write(f"# Sites and cuts: Biopython 1.88 Bio/Restriction/Restriction_Dictionary.py (REBASE emboss 404, 2024),"
                  f" sha256 {BIOPYTHON_DICT_SHA256}; Biopython License Agreement.\n")
        out.write(f"# Organism (<3>) and cognate methylation site (<6>): REBASE allenz.104, Copyright (c) Dr. Richard J. Roberts, 2021,"
                  f" as redistributed in crisprBase (MIT); sha256 {ALLENZ_SHA256}."
                  f" Redistribution permission: see docs/Validation/review-2026-09/REBASE-permission-request.md.\n")
        out.write(f"# Census: {len(bio)} Biopython enzymes = {one_cut} one-cut + {two_cut} two-cut (vendored, {len(rows)})"
                  f" + {undefined} with no defined cut (excluded). Not in allenz.104 (organism/methylation empty): {missing}.\n")
        out.write("# Columns: Name, Site, Forward (fst5), Reverse (size+fst3), SecondForward (scd5), SecondReverse (size+scd3), Organism, MethylationField\n")
        for row in rows:
            out.write("\t".join(row) + "\n")
    print(f"wrote {TABLE}: {len(rows)} enzymes ({one_cut} one-cut, {two_cut} two-cut); "
          f"excluded {undefined} undefined; {missing} not in allenz.104")


def read_table():
    rows = []
    with open(TABLE, encoding="utf-8") as f:
        for line in f:
            if line.startswith("#"):
                continue
            rows.append(line.rstrip("\n").split("\t"))
    return rows


def cmd_oracle(_args):
    import Bio
    from Bio.Seq import Seq
    from Bio import Restriction
    if Bio.__version__ != "1.88":
        sys.exit(f"Biopython {Bio.__version__}; 1.88 required")
    require_pin(os.path.join(os.path.dirname(Restriction.__file__), "Restriction_Dictionary.py"),
                BIOPYTHON_DICT_SHA256)

    rng = random.Random(SEED)
    lines = []
    for row in read_table():
        name, site = row[0], row[1]
        enzyme = getattr(Restriction, name)
        concrete = lambda: "".join(rng.choice(IUPAC[c]) for c in site)
        rand = lambda n: "".join(rng.choice("ACGT") for _ in range(n))
        forward = concrete()
        reverse = concrete().translate(COMPLEMENT)[::-1]
        seq = rand(rng.randint(0, 30)) + forward + rand(rng.randint(5, 40)) + reverse + rand(rng.randint(0, 30))
        linear = enzyme.search(Seq(seq), linear=True)
        circular = enzyme.search(Seq(seq), linear=False)
        lines.append("\t".join([name, seq, ",".join(map(str, linear)), ",".join(map(str, circular))]))

    os.makedirs(os.path.dirname(ORACLE), exist_ok=True)
    with open(ORACLE, "w", encoding="utf-8", newline="\n") as out:
        out.write("# Biopython 1.88 Bio.Restriction <enzyme>.search(Seq(seq), linear=True|False) for every enzyme of rebase_enzymes.tsv.\n")
        out.write(f"# GENERATED by scripts/rebase/generate_rebase_enzymes.py oracle (random.Random({SEED})); do not edit by hand.\n")
        out.write("# Each molecule carries the site once as written and once reverse-complemented (ambiguity codes resolved at random).\n")
        out.write("# Columns: Name, Sequence, LinearCuts (1-based top strand, Biopython order), CircularCuts (as Biopython returns them)\n")
        out.write("\n".join(lines) + "\n")
    print(f"wrote {ORACLE}: {len(lines)} enzymes")


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = parser.add_subparsers(dest="cmd", required=True)
    t = sub.add_parser("table")
    t.add_argument("--biopython-dict", required=True)
    t.add_argument("--allenz", required=True)
    t.set_defaults(func=cmd_table)
    o = sub.add_parser("oracle")
    o.set_defaults(func=cmd_oracle)
    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
