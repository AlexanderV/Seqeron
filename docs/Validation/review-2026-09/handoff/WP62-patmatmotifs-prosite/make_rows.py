#!/usr/bin/env python3
"""The argument lines of the runs (inputs/rows.txt): the input files x the options (-full, -prune, the report formats of the report layer, the report qualifiers)."""
import itertools, sys
files = ['one.fa', 'two.fa', 'five.fa', 'many.fa', 'damaged.fa', 'empty.fa']
rows = []
for f in files:
    for full in ['', '-full']:
        for prune in ['', '-noprune']:
            rows.append('-sequence %s %s %s -outfile out.rpt' % (f, full, prune))
for f in ['two.fa', 'many.fa']:
    for rf in ['dbmotif', 'gff', 'simple', 'table', 'excel', 'seqtable', 'motif', 'regions', 'srs', 'listfile', 'tagseq', 'embl', 'swiss', 'pir', 'nametable', 'draw', 'restrict', 'feattable', 'summary', 'debug', 'dasgff', 'genbank']:
        rows.append('-sequence %s -rformat %s -outfile out.rpt' % (f, rf))
for q in ['-raccshow', '-rdesshow', '-rscoreshow', '-rusashow', '-rstrandshow', '-rmaxall 5', '-rmaxseq 2', '-rdirectory outdir', '-rname myreport', '-rextension rpt', '-rformat dbmotif -rmaxall 1']:
    rows.append('-sequence many.fa %s -outfile out.rpt' % q)
rows += ['-sequence many.fa -sbegin 10 -send 120 -outfile out.rpt', '-sequence many.fa -sreverse -outfile out.rpt', '-sequence many.fa -slower -outfile out.rpt', '-sequence nothere.fa -outfile out.rpt', '-sequence many.fa -outfile stdout', '-sequence many.fa']
open(sys.argv[1], 'w').write('\n'.join(rows) + '\n')
