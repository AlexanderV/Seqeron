#!/usr/bin/env python3
"""Generate the EMBOSS 6.6.0 `restrict` oracle for the enzymes where EMBOSS and Biopython disagree.

Background (review 2026-09 B08 H7). A recognition window that reads as the site on BOTH strands
admits two binding orientations. Biopython reports one reading per window
(`Palindromic._search` applies `_modify` only; `NonPalindromic._search` takes the forward branch
whenever the forward named group matched), so when the two readings cut different bonds it loses
one physically real cut. EMBOSS `embPatRestrictScan` (nucleus/embpat.c) instead runs two
independent scan loops -- forward over `substr`, then `forward = ajFalse` over `revstr` -- and
afterwards keeps one hit per distinct `cut1`, preferring the forward hit. Seqeron follows EMBOSS.

Exactly five of the 754 REBASE enzymes with a defined cut are affected: HauII, SgeI, MspJI, LpnPI
and SgrTI (see `ReverseReadingCutsAlike` in RestrictionAnalyzer.cs for the derivation).

This script writes `tests/.../TestData/Restriction/emboss_restrict_bothstrand.tsv`, using the same
molecules as the Biopython fixture `biopython_rebase_search.tsv` so the two oracles are directly
comparable. It does NOT touch the frozen `rebase_enzymes.tsv` or
`generate_rebase_enzymes.py`.

Requirements: `emboss` (6.6.0) and `biopython==1.88`. EMBOSS ships no `embossre.enz` because
REBASE is licensed, so the enzyme file is generated here from Biopython's
`Restriction_Dictionary` and placed at `$EMBOSS_DATA/REBASE/embossre.enz`.

Usage:  python3 scripts/rebase/generate_emboss_restrict_oracle.py [output.tsv]
"""

from __future__ import annotations

import os
import re
import subprocess
import sys
import tempfile
from pathlib import Path

from Bio.Restriction.Restriction_Dictionary import rest_dict as REST_DICT

AFFECTED = ("HauII", "SgeI", "MspJI", "LpnPI", "SgrTI")
REPO = Path(__file__).resolve().parents[2]
BIOPYTHON_FIXTURE = (
    REPO / "tests/Seqeron/Seqeron.Genomics.Tests/TestData/Restriction/biopython_rebase_search.tsv"
)
DEFAULT_OUT = (
    REPO / "tests/Seqeron/Seqeron.Genomics.Tests/TestData/Restriction/emboss_restrict_bothstrand.tsv"
)


def write_enzyme_file(path: Path) -> int:
    """EMBOSS embossre.enz: name, pattern, len, ncuts, blunt, c1, c2, c3, c4 (tab separated).

    c1/c2 are the first 5'/3' cut positions counted from the first base of the recognition site,
    the cut lying AFTER base c; c3/c4 are the second pair of a two-sided (Type IIB) enzyme.
    Biopython stores c1 as `fst5` and c2 as `size + fst3` (likewise `scd5` / `size + scd3`).
    `blunt` is read by EMBOSS as an INTEGER -- a 'Y'/'N' column makes `restrict` load no enzyme at
    all and then report every requested name as absent from the database.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    written = 0
    with path.open("w") as fh:
        fh.write("# Generated from Biopython 1.88 Restriction_Dictionary by\n")
        fh.write("# scripts/rebase/generate_emboss_restrict_oracle.py -- not REBASE-derived data.\n")
        for name in sorted(REST_DICT):
            d = REST_DICT[name]
            if d.get("substrat") != "DNA" or d.get("fst5") is None or d.get("fst3") is None:
                continue
            size = d["size"]
            c1, c2 = d["fst5"], size + d["fst3"]
            if d.get("scd5") is not None and d.get("scd3") is not None:
                c3, c4, ncuts = d["scd5"], size + d["scd3"], 4
            else:
                c3, c4, ncuts = 0, 0, 2
            blunt = 1 if c1 == c2 else 0
            fh.write(
                f"{name.replace('_', '-')}\t{d['site']}\t{size}\t{ncuts}\t{blunt}"
                f"\t{c1}\t{c2}\t{c3}\t{c4}\n"
            )
            written += 1
    # restrict also opens REBASE/embossre.equ; an empty file is a valid "no equivalences".
    (path.parent / "embossre.equ").write_text("")
    return written


def restrict_cuts(data_dir: Path, sequence: str, enzyme: str, circular: bool) -> list[int]:
    """The 1-based top-strand cut positions EMBOSS reports, deduped and sorted."""
    env = dict(os.environ, EMBOSS_DATA=str(data_dir))
    with tempfile.TemporaryDirectory() as tmp:
        fasta, out = Path(tmp) / "in.fa", Path(tmp) / "out.restrict"
        fasta.write_text(f">oracle\n{sequence}\n")
        cmd = [
            "restrict", "-sequence", str(fasta), "-enzymes", enzyme,
            "-sitelen", "2", "-outfile", str(out), "-auto",
        ]
        if circular:
            cmd.append("-plasmid")
        subprocess.run(cmd, env=env, capture_output=True, check=False)
        if not out.exists():
            return []
        cuts = set()
        for line in out.read_text().splitlines():
            line = line.strip()
            if not line or line.startswith("#") or line.startswith("Start"):
                continue
            fields = re.split(r"\s+", line)
            if len(fields) >= 7 and fields[3].lower() == enzyme.lower():
                # column 6 is `5prime`, the 0-based top-strand bond
                cuts.add(int(fields[5]) + 1)
        return sorted(cuts)


def main() -> int:
    out_path = Path(sys.argv[1]) if len(sys.argv) > 1 else DEFAULT_OUT
    molecules: dict[str, str] = {}
    for line in BIOPYTHON_FIXTURE.read_text().splitlines():
        if not line or line.startswith("#"):
            continue
        name, sequence = line.split("\t")[:2]
        if name in AFFECTED:
            molecules[name] = sequence
    missing = set(AFFECTED) - set(molecules)
    if missing:
        raise SystemExit(f"molecules missing from {BIOPYTHON_FIXTURE.name}: {sorted(missing)}")

    with tempfile.TemporaryDirectory() as tmp:
        data_dir = Path(tmp) / "data"
        count = write_enzyme_file(data_dir / "REBASE" / "embossre.enz")
        print(f"enzyme file: {count} enzymes", file=sys.stderr)

        rows = []
        for name in AFFECTED:
            seq = molecules[name]
            rows.append((
                name, seq,
                restrict_cuts(data_dir, seq, name, circular=False),
                restrict_cuts(data_dir, seq, name, circular=True),
            ))

    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w") as fh:
        fh.write("# EMBOSS 6.6.0 `restrict` cut positions for the five enzymes where EMBOSS and\n")
        fh.write("# Biopython disagree: a window reading as the site on both strands admits two\n")
        fh.write("# cleavages, which EMBOSS reports and Biopython collapses (B08 H7).\n")
        fh.write("# GENERATED by scripts/rebase/generate_emboss_restrict_oracle.py; do not edit by hand.\n")
        fh.write("# Molecules are the same rows as biopython_rebase_search.tsv.\n")
        fh.write("# Columns: Name, Sequence, LinearCuts, CircularCuts (1-based top strand, sorted, deduped)\n")
        for name, seq, lin, circ in rows:
            fh.write(f"{name}\t{seq}\t{','.join(map(str, lin))}\t{','.join(map(str, circ))}\n")
    print(f"wrote {out_path}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
