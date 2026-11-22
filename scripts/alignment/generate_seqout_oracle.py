#!/usr/bin/env python3
"""A5-8c3h (B09 F252) oracle: the sequence WRITERS of EMBOSS 6.6.0 (ajseqwrite.c, `seqOutFormat`) and the output options of seqret
    seqret -auto -sequence in.dat [-sformat1 F] -outseq OUT -osformat2 FMT [-ossingle2 [-osextension2 E] [-osdirectory2 sub]] [-osname2 N] [-firstonly]
(the output files are collected: the sequence file and, with -ossingle, one file per sequence). The input sequences are the inputs of the reader oracle (generate_seqret_oracle.py: writer outputs of seqret, entries with feature
tables, damaged texts) so that names, descriptions, accessions, versions, nucleotide / protein types, gaps and lower case reach the writer in every way the readers give them.
JSONL.GZ rows: id, fmt (the output format), text, spec, args, osformat, outseq, single, osext, osdir, osname, firstonly, epoch, rc, fa (the sequence file), files ({path: text} of -ossingle), err.
Usage: python3 scripts/alignment/generate_seqout_oracle.py OUT.jsonl.gz PER_FORMAT [seed] [format ...]      (needs EMBOSS 6.6.0 and gcc)
"""
import gzip, json, os, random, shutil, subprocess, sys, tempfile, textwrap

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import generate_seqret_oracle as G  # noqa: E402

FORMATS = ['unknown', 'fasta', 'pearson', 'gcg', 'gcg8', 'nbrf', 'pir', 'codata', 'strider', 'staden', 'text', 'plain', 'raw', 'ig', 'iguspto', 'gde', 'fitch', 'acedb', 'mase', 'ncbi', 'gifasta', 'nib', 'nibble', 'dasdna', 'das', 'asn1',
           'fastq-sanger', 'fastq-illumina', 'fastq-solexa', 'fastq', 'sam', 'gff2', 'gff3', 'gff',
           'msf', 'clustal', 'aln', 'selex', 'phylip', 'phylipnon', 'phylip3', 'hennig86', 'mega', 'meganon', 'nexus', 'paup', 'nexusnon', 'paupnon', 'jackknifer', 'jackknifernon', 'treecon', 'scf', 'embl', 'em', 'emblnew', 'emold', 'emblold']
# the writers that print the accuracy values (quality) of the sequences: their inputs are mostly texts with random qualities, and every run is checked with valgrind (the buffer of the values keeps what it had)
ENTRY_FORMATS = {'embl', 'em', 'emblnew', 'emold', 'emblold', 'genbank', 'gb', 'ddbj', 'refseq', 'genpept', 'refseqp', 'swiss', 'sw', 'swissprot', 'uniprot', 'swissnew', 'swnew', 'swold', 'swissold', 'swissprotold', 'swissprotnew', 'experiment', 'debug'}
QUAL_FORMATS = {'fastq-sanger', 'fastq-illumina', 'fastq-solexa', 'fastq', 'sam', 'scf', 'bam', 'experiment', 'debug'}
QUAL_INPUTS = ['qsanger', 'qsanger', 'qillumina', 'qsolexa', 'qsam', 'qexperiment', 'fastq-sanger', 'fastq-illumina', 'fastq-solexa', 'sam', 'experiment']
# the inputs of the readers whose reference runs are defined (the nexus and mega readers of the library read freed or foreign memory on some inputs, see generate_seqret_oracle.py: memerr)
INPUTS = ['ncbi', 'ncbi', 'ncbi', 'fasta', 'fasta', 'fasta', 'gcg', 'embl', 'genbank', 'nbrf', 'swiss', 'ig', 'codata', 'strider', 'staden', 'gde', 'fitch', 'treecon', 'pir', 'clustal', 'msf', 'phylip', 'raw', 'textonly', 'asis', 'fastq', 'sam']
# the entry formats whose readers keep the annotation (references, keywords, taxonomy, comments, dates) that the entry writers print: SEQOUT_INPUTS=annot makes the inputs of the entry writers mostly such texts
ANNOT_INPUTS = ['embl', 'embl', 'embl', 'embl', 'genbank', 'genbank', 'swiss', 'swiss', 'refseq', 'genpept', 'refseqp', 'experiment', 'nbrf', 'fasta', 'ncbi', 'gcg', 'pir', 'sam']
PROT_ANNOT_INPUTS = ['genpept', 'genpept', 'genpept', 'refseqp', 'refseqp', 'swiss', 'swiss', 'swiss', 'fasta', 'ncbi', 'pir', 'gcg']
ANNOT_READY = set(os.environ.get('SEQOUT_READY', 'embl,nbrf,experiment,fasta,ncbi,gcg,pir,sam,genbank,refseq,genpept,refseqp,swiss').split(','))


DBNAME_FORCE = os.environ.get('SEQOUT_OSDBNAME') == '1'     # a run that gives -osdbname to most inputs (the file seqout_oracle_osdbname.jsonl.gz)


def ncbi_text(rnd):
    """Entries with the identifiers of the NCBI FASTA headers (gi, database, accession.version, name, description) in every form the NCBI reader knows, so that Gi, Db, Setdb, Sv and Acc reach the ncbi and gifasta writers."""
    dbs = ['gb', 'emb', 'dbj', 'sp', 'ref', 'pir', 'prf', 'tpd', 'tpe', 'tpg', 'lcl', 'bbs', 'gnl|mydb', 'gnl|BL_ORD_ID', 'pdb', 'gnl|Seq', 'unk', 'GB', 'LCL']
    out = []
    for i in range(rnd.randint(1, 4)):
        seq = ''.join(rnd.choice('ACGT') for _ in range(rnd.choice([5, 30, 60, 61, 130])))
        acc = rnd.choice(['AB%06d' % rnd.randint(1, 999999), 'NM_%06d' % rnd.randint(1, 999999), 'P%05d' % rnd.randint(1, 99999), 'X%d' % rnd.randint(1, 99)])
        ver = acc + rnd.choice(['', '.1', '.2', '.12'])
        db = rnd.choice(dbs)
        gi = str(rnd.randint(1, 99999999))
        name = rnd.choice(['NAME_%d' % rnd.randint(1, 9), 'abc', ver, acc, gi])
        desc = rnd.choice(['', 'some protein', 'x', 'alpha beta gamma'])
        k = rnd.random()
        if k < 0.3: head = 'gi|%s|%s|%s|%s' % (gi, db if '|' not in db else 'gnl', ver, name) if '|' not in db else 'gi|%s|%s|%s' % (gi, db, name)
        elif k < 0.45: head = 'gi|%s' % gi
        elif k < 0.6: head = '%s|%s|%s' % (db if '|' not in db else db.replace('|', '|'), ver, name) if '|' not in db else '%s|%s' % (db, name)
        elif k < 0.7: head = 'lcl|%s' % name
        elif k < 0.8: head = 'gnl|BL_ORD_ID|%d' % rnd.randint(0, 20)
        else: head = name
        out.append('>%s%s\n%s\n' % (head, (' ' + desc) if desc else '', seq))
    return ''.join(out)


def qual_text(rnd, kind):
    """Entries of the quality formats with random quality characters (fastq-sanger, fastq-illumina, fastq-solexa, sam, experiment); a few have a quality line of the wrong length, qualities below the minimum, or no quality at all (sam '*')."""
    out = []
    qmin = {'qsanger': 33, 'qillumina': 64, 'qsolexa': 59, 'qsam': 33}.get(kind, 33)
    for i in range(rnd.choice([1, 1, 2, 3, 5])):
        n = rnd.choice([1, 4, 10, 20, 21, 40, 61])
        seq = ''.join(rnd.choice('ACGTN') if rnd.random() < 0.95 else rnd.choice('acgtRYK') for _ in range(n))
        name = rnd.choice(['r%d' % i, 'read_%d' % i, 'HWI-ST:%d:1' % i, 'A%d' % (i + 7)])
        desc = rnd.choice(['', '', 'x y', 'length=%d' % n, 'some description'])
        qn = n
        if rnd.random() < 0.1: qn = max(0, n + rnd.choice([-3, -1, 1, 4]))
        lo = qmin if rnd.random() < 0.85 else max(33, qmin - rnd.choice([1, 5, 20]))
        hi = rnd.choice([qmin + 10, qmin + 40, 126, 90])
        qual = ''.join(chr(rnd.randint(lo, min(hi, 126))) for _ in range(qn))
        if kind == 'qexperiment':
            av = ' '.join(rnd.choice(['%d' % rnd.randint(0, 99), '%d,%d' % (rnd.randint(0, 9), rnd.randint(0, 99)), '%d' % rnd.randint(0, 9)]) for _ in range(qn))
            body = ''.join('SQ   Sequence %d BP; 0 A; 0 C; 0 G; 0 T; 0 other;\n' % n if False else '') or ''
            lines = ['ID   %s' % name, 'EX   %s' % (desc or 'experiment')]
            for j in range(0, len(av), 70): lines.append('AV   ' + av[j:j + 70])
            lines.append('SQ   Sequence %d BP; 0 A; 0 C; 0 G; 0 T; 0 other;' % n)
            for j in range(0, n, 60): lines.append('     ' + seq[j:j + 60].lower())
            lines.append('//')
            out.append('\n'.join(lines) + '\n')
        elif kind == 'qsam':
            q = qual if rnd.random() < 0.8 and qn == n else '*'
            out.append('%s\t4\t*\t0\t0\t*\t*\t0\t0\t%s\t%s\n' % (name, seq, q))
        else:
            out.append('@%s%s\n%s\n+%s\n%s\n' % (name, (' ' + desc) if desc else '', seq, name if rnd.random() < 0.3 else '', qual))
    text = ''.join(out)
    if kind == 'qsam' and rnd.random() < 0.8: text = '@HD\tVN:1.0\n' + text
    return text


MONTHS = ['JAN', 'FEB', 'MAR', 'APR', 'MAY', 'JUN', 'JUL', 'AUG', 'SEP', 'OCT', 'NOV', 'DEC']
TAXA = [('Homo sapiens', 'human', ['Eukaryota', 'Metazoa', 'Chordata', 'Craniata', 'Vertebrata', 'Mammalia']), ('Escherichia coli', '', ['Bacteria', 'Proteobacteria', 'Gammaproteobacteria']),
        ('Mus musculus', 'house mouse', ['Eukaryota', 'Metazoa', 'Chordata']), ('Arabidopsis thaliana', 'thale cress', ['Eukaryota', 'Viridiplantae', 'Streptophyta', 'Embryophyta', 'Tracheophyta', 'Spermatophyta'])]
TITLES = ['A title', 'Cloning and sequencing of the gene for a very long protein kinase that regulates the cell cycle in many different organisms and tissues',
          'Short title', 'Direct Submission', 'Structure of the complex between two proteins at 2.5 Angstrom resolution; implications for function']
JOURNALS = ['Nature 1:1-2(2000).', 'Submitted (12-JAN-2000) to the INSDC. Some Institute, 1 Some Street, Some City, SOME COUNTRY.', 'J. Mol. Biol. 100:200-210(1990).', '(er) Thesis (1999), University of Somewhere.', 'Patent number US 1234567-A/1, 01-JAN-1990.']
AUTHORS = ['Smith J.', 'Smith J., Jones K.', 'Smith J., Jones K., Lee M.', 'Smith J., Jones K., Lee M., Brown A.B., White C.D., Black E.F., Green G.H., Gray I.J., Blue K.L.', 'The Consortium']


def rdate(rnd):
    return '%02d-%s-%d' % (rnd.randint(1, 31), rnd.choice(MONTHS), rnd.randint(1980, 2024))


def rseq_lines(rnd, seq, fmt):
    out = []
    if fmt == 'embl':
        for i in range(0, len(seq), 60):
            chunk = seq[i:i + 60]
            blocks = ' '.join(chunk[j:j + 10] for j in range(0, len(chunk), 10))
            out.append('     %-66s %9d' % (blocks, min(i + 60, len(seq))))
    else:
        for i in range(0, len(seq), 60):
            chunk = seq[i:i + 60]
            out.append('%9d %s' % (i + 1, ' '.join(chunk[j:j + 10] for j in range(0, len(chunk), 10))))
    return '\n'.join(out) + '\n'


def embl_rich(rnd, prot=False):
    """An EMBL entry with every kind of annotation line: references (RN RC RP RX RG RA RT RL), DR cross-references, KW, OS OC OG, CC comments, DT dates, in the new (release 87) and the old ID line forms."""
    out = ''
    for k in range(rnd.choice([1, 1, 2])):
        n = rnd.choice([12, 60, 61, 130, 301])
        seq = ''.join(rnd.choice('acgt') for _ in range(n))
        name = rnd.choice(['AB000001', 'X56789', 'SEQ1', 'HSA_1', 'U12345'])
        acc = [name] + rnd.sample(['AB000002', 'X1', 'Y22222', 'Z333333'], rnd.choice([0, 0, 1, 2, 3]))
        if rnd.random() < 0.8:
            out += 'ID   %s; SV %d; %s; %s; %s; %s; %d BP.\nXX\n' % (name, rnd.randint(1, 5), rnd.choice(['linear', 'circular']), rnd.choice(['genomic DNA', 'mRNA', 'genomic RNA', 'tRNA', 'other DNA', 'unassigned DNA', 'bogus']), rnd.choice(['STD', 'CON', 'EST', 'HTG', 'PAT', 'XXX']), rnd.choice(['HUM', 'PRO', 'VRL', 'UNC', 'MUS', 'ENV', 'FOO']), n)
        else:
            out += 'ID   %s  standard; DNA; UNC; %d BP.\n' % (name, n)
        if rnd.random() < 0.9: out += 'AC   %s;\nXX\n' % '; '.join(acc)
        if rnd.random() < 0.5: out += 'SV   %s.%d\nXX\n' % (name, rnd.randint(1, 3))
        if rnd.random() < 0.6:
            out += 'DT   %s (Rel. %d, Created)\nDT   %s (Rel. %d, Last updated, Version %d)\nXX\n' % (rdate(rnd), rnd.randint(1, 120), rdate(rnd), rnd.randint(1, 120), rnd.randint(1, 9))
        if rnd.random() < 0.9:
            d = rnd.choice(['A short description.', 'Homo sapiens mRNA for a very long named protein kinase that is involved in the regulation of the cell cycle, complete cds.', 'x'])
            out += ''.join('DE   %s\n' % d[i:i + 70] for i in range(0, len(d), 70)) + 'XX\n'
        if rnd.random() < 0.6: out += 'KW   %s.\nXX\n' % '; '.join(rnd.sample(['protein kinase', 'cell cycle', 'EST', 'mRNA', 'HTG', 'complete cds'], rnd.choice([1, 2, 3])))
        if rnd.random() < 0.7:
            sp, common, tax = rnd.choice(TAXA)
            out += 'OS   %s%s\n' % (sp, (' (%s)' % common) if common else '')
            oc = '; '.join(tax) + '.'
            out += ''.join('OC   %s\n' % oc[i:i + 70] for i in range(0, len(oc), 70))
            if rnd.random() < 0.3: out += 'OG   %s\n' % rnd.choice(['Mitochondrion', 'Plastid; Chloroplast'])
            out += 'XX\n'
        for r in range(rnd.choice([0, 1, 2, 3])):
            out += 'RN   [%d]\n' % (r + 1)
            if rnd.random() < 0.3: out += 'RC   %s\n' % rnd.choice(['revised', 'a comment that runs over more than one line because it is very long indeed and keeps going'])
            if rnd.random() < 0.8: out += 'RP   %d-%d\n' % (rnd.randint(1, 5), rnd.randint(6, 300))
            if rnd.random() < 0.5: out += 'RX   %s; %d.\n' % (rnd.choice(['PUBMED', 'DOI', 'MEDLINE']), rnd.randint(100, 99999999))
            if rnd.random() < 0.15: out += 'RG   The Gene Consortium\n'
            a = rnd.choice(AUTHORS)
            out += ''.join('RA   %s\n' % l for l in (textwrap.wrap(a + ';', 70) or [a + ';']))
            t = rnd.choice(TITLES)
            out += ''.join('RT   %s\n' % l for l in textwrap.wrap('"' + t + '";', 70)) if rnd.random() < 0.85 else 'RT   ;\n'
            out += 'RL   %s\n' % rnd.choice(JOURNALS)
            out += 'XX\n'
        if rnd.random() < 0.5:
            for d in range(rnd.choice([1, 2, 3])):
                out += 'DR   %s\n' % rnd.choice(['EMBL; AB000001; BAA00001.1.', 'UniProtKB/Swiss-Prot; P12345; ABC_HUMAN.', 'GO; GO:0000001; C:cell; IDA:UniProtKB.', 'PDB; 1ABC; X-ray; 2.50 A; A=1-100.', 'MIM; 100100; gene.'])
            out += 'XX\n'
        if rnd.random() < 0.5:
            c = rnd.choice(['a comment', 'Some comment, with commas, that is long enough to wrap onto more than one line of the entry when it is written out again by the writer', 'first line\nsecond line'])
            out += ''.join('CC   %s\n' % l for l in c.split('\n')) + 'XX\n'
        out += 'SQ   Sequence %d BP; %d A; %d C; %d G; %d T; 0 other;\n' % (n, seq.count('a'), seq.count('c'), seq.count('g'), seq.count('t'))
        out += rseq_lines(rnd, seq, 'embl') + '//\n'
    return out


def genbank_rich(rnd, kind='genbank'):
    """A GenBank entry with DEFINITION, ACCESSION, VERSION, KEYWORDS, SOURCE / ORGANISM, REFERENCE blocks (AUTHORS TITLE JOURNAL PUBMED), COMMENT and ORIGIN."""
    out = ''
    for k in range(rnd.choice([1, 1, 2])):
        n = rnd.choice([12, 60, 61, 130, 301])
        seq = ''.join(rnd.choice('acgt') for _ in range(n))
        name = rnd.choice(['AB000001', 'X56789', 'NM_000001', 'HSA_1', 'U12345'])
        mol = rnd.choice(['DNA', 'mRNA', 'RNA', 'ss-RNA', 'tRNA'])
        out += 'LOCUS       %-16s %6d bp    %-7s %s %s %s\n' % (name, n, mol, rnd.choice(['linear', 'circular']), rnd.choice(['PRI', 'BCT', 'VRL', 'EST', 'CON', 'UNA', 'ROD']), rdate(rnd))
        d = rnd.choice(['A short description.', 'Homo sapiens mRNA for a very long named protein kinase that is involved in the regulation of the cell cycle, complete cds.', 'x'])
        out += 'DEFINITION  %s\n' % '\n            '.join(textwrap.wrap(d, 67))
        out += 'ACCESSION   %s\n' % ' '.join([name] + rnd.sample(['AB000002', 'X1', 'Y22222'], rnd.choice([0, 1, 2])))
        out += 'VERSION     %s.%d%s\n' % (name, rnd.randint(1, 3), rnd.choice(['', '  GI:%d' % rnd.randint(1, 99999999)]))
        if rnd.random() < 0.6: out += 'KEYWORDS    %s.\n' % '; '.join(rnd.sample(['protein kinase', 'cell cycle', 'EST', 'mRNA', 'HTG', 'complete cds'], rnd.choice([1, 2, 3])))
        if rnd.random() < 0.8:
            sp, common, tax = rnd.choice(TAXA)
            out += 'SOURCE      %s%s\n  ORGANISM  %s\n' % (sp, (' (%s)' % common) if common else '', sp)
            lines = textwrap.wrap('; '.join(tax) + '.', 67)
            out += ''.join('            %s\n' % l for l in lines)
        for r in range(rnd.choice([0, 1, 2, 3])):
            out += 'REFERENCE   %d' % (r + 1) + rnd.choice(['', '  (bases 1 to %d)' % n, '  (bases 1 to 50; 60 to 100)', '  (residues 1 to 3)']) + '\n'
            if rnd.random() < 0.9:
                a = rnd.choice(AUTHORS).replace('. ', '.').replace(', ', ',')
                lines = textwrap.wrap(a, 67)
                out += '  AUTHORS   %s\n' % lines[0] + ''.join('            %s\n' % l for l in lines[1:])
            if rnd.random() < 0.85:
                lines = textwrap.wrap(rnd.choice(TITLES), 67)
                out += '  TITLE     %s\n' % lines[0] + ''.join('            %s\n' % l for l in lines[1:])
            out += '  JOURNAL   %s\n' % rnd.choice(['Nature 1 (1), 1-2 (2000)', 'Submitted (12-JAN-2000) to the INSDC. Some Institute', 'J. Mol. Biol. 100 (2), 200-210 (1990)'])
            if rnd.random() < 0.5: out += '   PUBMED   %d\n' % rnd.randint(100, 99999999)
        if rnd.random() < 0.5:
            c = rnd.choice(['a comment', 'Some comment, with commas, that is long enough to wrap onto more than one line of the entry', 'first paragraph\n\nsecond paragraph'])
            ls = c.split('\n')
            out += 'COMMENT     %s\n' % ls[0] + ''.join('            %s\n' % l for l in ls[1:])
        out += 'ORIGIN\n' + rseq_lines(rnd, seq, 'genbank') + '//\n'
    return out


def genpept_rich(rnd):
    """A GenPept / RefSeq protein entry: LOCUS (7 fields), DEFINITION, ACCESSION, VERSION, KEYWORDS, SOURCE with ORGANISM and the lineage, REFERENCE blocks with AUTHORS TITLE JOURNAL PUBMED, COMMENT, ORIGIN."""
    out = ''
    for k in range(rnd.choice([1, 1, 2])):
        n = rnd.choice([12, 60, 61, 130, 301])
        seq = ''.join(rnd.choice('acdefghiklmnpqrstvwy') for _ in range(n))
        name = rnd.choice(['NP_000001', 'AAA12345', 'P12345', 'XP_1'])
        out += 'LOCUS       %-12s %6d aa    linear   %s %s\n' % (name, n, rnd.choice(['PRI', 'BCT', 'VRL', 'ROD', 'UNA', 'FOO']), rdate(rnd))
        d = rnd.choice(['A short description.', 'protein kinase C alpha type [Homo sapiens] with a long description that wraps over lines of the entry', 'x'])
        out += 'DEFINITION  %s\n' % '\n            '.join(textwrap.wrap(d, 67))
        out += 'ACCESSION   %s\n' % ' '.join([name] + rnd.sample(['AB000002', 'X1', 'Y22222'], rnd.choice([0, 1, 2])))
        out += 'VERSION     %s.%d%s\n' % (name, rnd.randint(1, 3), rnd.choice(['', '  GI:%d' % rnd.randint(1, 99999999)]))
        if rnd.random() < 0.6: out += 'KEYWORDS    %s.\n' % '; '.join(rnd.sample(['RefSeq', 'kinase', 'cell cycle'], rnd.choice([1, 2, 3])))
        if rnd.random() < 0.8:
            sp, common, tax = rnd.choice(TAXA)
            out += 'SOURCE      %s%s\n  ORGANISM  %s\n' % (sp, (' (%s)' % common) if common else '', sp)
            out += ''.join('            %s\n' % l for l in textwrap.wrap('; '.join(tax) + '.', 67))
        for r in range(rnd.choice([0, 1, 2, 3])):
            out += 'REFERENCE   %d' % (r + 1) + rnd.choice(['', '  (residues 1 to %d)' % n, '  (residues 1 to 50; 60 to 100)']) + '\n'
            if rnd.random() < 0.9:
                a = rnd.choice(AUTHORS).replace('. ', '.').replace(', ', ',')
                lines = textwrap.wrap(a, 67)
                out += '  AUTHORS   %s\n' % lines[0] + ''.join('            %s\n' % l for l in lines[1:])
            if rnd.random() < 0.85:
                lines = textwrap.wrap(rnd.choice(TITLES), 67)
                out += '  TITLE     %s\n' % lines[0] + ''.join('            %s\n' % l for l in lines[1:])
            out += '  JOURNAL   %s\n' % rnd.choice(['Nature 1 (1), 1-2 (2000)', 'Submitted (12-JAN-2000) to the INSDC. Some Institute', 'J. Mol. Biol. 100 (2), 200-210 (1990)'])
            if rnd.random() < 0.5: out += '   PUBMED   %d\n' % rnd.randint(100, 99999999)
        if rnd.random() < 0.5:
            c = rnd.choice(['a comment', 'Some comment, with commas, that is long enough to wrap onto more than one line of the entry', 'first paragraph\n\nsecond paragraph'])
            ls = c.split('\n')
            out += 'COMMENT     %s\n' % ls[0] + ''.join('            %s\n' % l for l in ls[1:])
        out += 'ORIGIN\n' + rseq_lines(rnd, seq, 'genbank') + '//\n'
    return out


def swiss_rich(rnd):
    """A Swiss-Prot entry for a protein with ID, AC, DT, DE, GN, OS, OC, OG, OX, RN / RP / RC / RX / RG / RA / RT / RL, CC, DR, KW, FT and SQ lines."""
    out = ''
    for k in range(rnd.choice([1, 1, 2])):
        n = rnd.choice([12, 60, 61, 130, 301])
        seq = ''.join(rnd.choice('ACDEFGHIKLMNPQRSTVWY') for _ in range(n))
        name = rnd.choice(['ABC_HUMAN', 'KPCA_MOUSE', 'Q12345_ECOLI', 'TEST_1'])
        out += 'ID   %-24s Reviewed;         %d AA.\n' % (name, n)
        out += 'AC   %s;\n' % '; '.join(['P12345'] + rnd.sample(['Q99999', 'O00001', 'A0A000'], rnd.choice([0, 1, 2])))
        out += 'DT   %s, integrated into UniProtKB/Swiss-Prot.\nDT   %s, sequence version %d.\nDT   %s, entry version %d.\n' % (rdate(rnd), rdate(rnd), rnd.randint(1, 4), rdate(rnd), rnd.randint(1, 150))
        dl = rnd.choice([['RecName: Full=Protein kinase C alpha type;', '         Short=PKC-alpha;', '         EC=2.7.11.13;'], ['SubName: Full=Putative uncharacterized protein;'],
                         ['RecName: Full=A very long protein name that goes over the line;', 'AltName: Full=Other name;', 'Flags: Precursor;']])
        out += ''.join('DE   %s\n' % l for l in dl) if rnd.random() < 0.9 else ''
        if rnd.random() < 0.7: out += 'GN   Name=%s; Synonyms=ABC1, XYZ; ORFNames=b0001;\n' % rnd.choice(['PRKCA', 'abc', 'gene1'])
        sp, common, tax = rnd.choice(TAXA)
        out += 'OS   %s%s.\n' % (sp, (' (%s)' % common) if common else '')
        lines = textwrap.wrap('; '.join(tax) + '.', 70)
        out += ''.join('OC   %s\n' % l for l in lines)
        if rnd.random() < 0.7: out += 'OX   NCBI_TaxID=%d;\n' % rnd.choice([9606, 10090, 562, 3702])
        for r in range(rnd.choice([0, 1, 2, 3])):
            out += 'RN   [%d]\nRP   %s\n' % (r + 1, rnd.choice(['NUCLEOTIDE SEQUENCE [MRNA].', 'X-RAY CRYSTALLOGRAPHY (2.5 ANGSTROMS).']))
            if rnd.random() < 0.5: out += 'RC   STRAIN=C57BL/6J; TISSUE=Brain;\n'
            if rnd.random() < 0.6: out += 'RX   PubMed=%d; DOI=10.1000/xyz;\n' % rnd.randint(100, 99999999)
            if rnd.random() < 0.15: out += 'RG   The Consortium;\n'
            a = rnd.choice(AUTHORS)
            out += ''.join('RA   %s\n' % l for l in textwrap.wrap(a + ';', 70))
            out += ''.join('RT   %s\n' % l for l in textwrap.wrap('"' + rnd.choice(TITLES) + '";', 70)) if rnd.random() < 0.8 else ''
            out += 'RL   %s\n' % rnd.choice(['Nature 1:1-2(2000).', 'Submitted (JAN-2000) to the EMBL/GenBank/DDBJ databases.'])
        if rnd.random() < 0.6:
            c = rnd.choice(['-!- FUNCTION: A function that is described here.', '-!- SUBCELLULAR LOCATION: Cytoplasm.\nCC       Nucleus.', '-!- SIMILARITY: Belongs to the protein kinase superfamily.'])
            out += ''.join('CC   %s\n' % l for l in c.split('\n'))
        if rnd.random() < 0.5: out += 'DR   EMBL; AB000001; BAA00001.1; -; mRNA.\nDR   PDB; 1ABC; X-ray; 2.50 A; A=1-100.\n'
        if rnd.random() < 0.6: out += 'KW   %s.\n' % '; '.join(rnd.sample(['Kinase', 'Transferase', 'ATP-binding', 'Cytoplasm'], rnd.choice([1, 2, 3])))
        mw = rnd.randint(1000, 90000)
        out += 'SQ   SEQUENCE %5d AA; %6d MW;  %016X CRC64;\n' % (n, mw, rnd.getrandbits(64))
        for i in range(0, n, 60):
            chunk = seq[i:i + 60]
            out += '     ' + ' '.join(chunk[j:j + 10] for j in range(0, len(chunk), 10)) + '\n'
        out += '//\n'
    return out


def annot_text(rnd, src):
    """(text, how) of an annotated entry for the source name, or None when the generic reader inputs are to be used"""
    if src == 'embl' and rnd.random() < 0.85: return embl_rich(rnd), 'rich'
    if src in ('genbank', 'refseq') and rnd.random() < 0.85: return genbank_rich(rnd), 'rich'
    if src == 'swiss' and rnd.random() < 0.85: return swiss_rich(rnd), 'rich'
    if src in ('genpept', 'refseqp') and rnd.random() < 0.85: return genpept_rich(rnd), 'rich'
    return None


def collect(work):
    files = {}
    for root, _, names in os.walk(work):
        for n in names:
            p = os.path.join(root, n)
            rel = os.path.relpath(p, work)
            if rel == 'in.dat':
                continue
            files[rel] = open(p, encoding='latin-1', newline='').read()
    return files


def main(out, per, seed, formats):
    rnd = random.Random(seed)
    rows = []
    with tempfile.TemporaryDirectory() as d:
        shim = os.path.join(d, 'faketime.so')
        subprocess.run(['gcc', '-shared', '-fPIC', '-o', shim, os.path.join(os.path.dirname(os.path.abspath(__file__)), 'faketime.c')], check=True)
        work = os.path.join(d, 'w'); os.mkdir(work)
        gen = os.path.join(d, 'g'); os.mkdir(gen)
        for fmt in formats:
            made = 0
            while made < per:
                if fmt in ENTRY_FORMATS:
                    # the writers of protein entries are given protein sources (the others reject nucleotide sequences)
                    pool = PROT_ANNOT_INPUTS if fmt in ('refseqp', 'swiss', 'sw', 'swissprot', 'uniprot', 'swissnew', 'swnew', 'swissprotnew', 'swold', 'swissold', 'swissprotold') else ANNOT_INPUTS
                    src = rnd.choice([x for x in pool if x in ANNOT_READY])
                else:
                    src = rnd.choice(QUAL_INPUTS if fmt in QUAL_FORMATS and rnd.random() < 0.8 else INPUTS)
                rich = annot_text(rnd, src) if fmt in ENTRY_FORMATS else None
                if rich:
                    text, how = rich
                elif src.startswith('q'):
                    text, how = qual_text(rnd, src), 'qual'
                    src = {'qsanger': 'fastq-sanger', 'qillumina': 'fastq-illumina', 'qsolexa': 'fastq-solexa', 'qsam': 'sam', 'qexperiment': 'experiment'}[src]
                else:
                    text, how = (ncbi_text(rnd), 'ncbi') if src == 'ncbi' else G.make_input(rnd, gen, src)
                if rnd.random() < 0.15: text, how = G.damaged(rnd, text), how + '+damaged'
                if rnd.random() < 0.05: text = text.rstrip('\n')
                spec, args = 'in.dat', []
                if rnd.random() < 0.7: args = ['-sformat1', 'fasta' if src == 'ncbi' else src]
                if rnd.random() < 0.08: args += ['-firstonly']
                if rnd.random() < 0.08: args += [rnd.choice(['-slower1', '-supper1', '-sprotein1', '-snucleotide1'])]
                osformat, outseq = fmt, rnd.choice(['o.out', 'out', 'o.txt'])
                if rnd.random() < 0.3:
                    osformat, outseq = None, '%s::%s' % (rnd.choice([fmt, fmt.upper()]), outseq)
                single = rnd.random() < 0.25
                osext = rnd.choice(['abc', 'X1', 'Txt']) if single and rnd.random() < 0.4 else None
                osdir = rnd.choice(['sub', 'sub/']) if single and rnd.random() < 0.4 else None
                osname = rnd.choice(['nm', 'N1']) if rnd.random() < 0.2 else None
                osdbname = rnd.choice(['MYDB', 'db1', 'sp', 'gb', 'lcl', 'BL_ORD_ID']) if DBNAME_FORCE and rnd.random() < 0.9 else None
                if rnd.random() < 0.03: osformat, outseq = 'bogus', 'o.out'
                for nme in os.listdir(work):
                    p = os.path.join(work, nme)
                    if os.path.isdir(p):
                        for f in os.listdir(p): os.remove(os.path.join(p, f))
                        os.rmdir(p)
                    else:
                        os.remove(p)
                os.mkdir(os.path.join(work, 'sub'))
                open(os.path.join(work, 'in.dat'), 'w', newline='\n', encoding='latin-1').write(text)
                epoch = rnd.randint(978307200, 2051222400)
                env = dict(os.environ, LD_PRELOAD=shim, FAKE_EPOCH=str(epoch))
                cmd = ['-auto', '-sequence', spec] + args + ['-outseq', outseq]
                if osformat: cmd += ['-osformat2', osformat]
                if single: cmd += ['-ossingle2']
                if osext: cmd += ['-osextension2', osext]
                if osdir: cmd += ['-osdirectory2', osdir]
                if osname: cmd += ['-osname2', osname]
                if osdbname: cmd += ['-osdbname2', osdbname]
                rc, err, hang = G.seqret(work, cmd, env=env)
                files = collect(work) if not hang else {}
                memerr = False
                if fmt in QUAL_FORMATS and not hang and shutil.which('valgrind'):
                    # a run that reads memory the library never wrote has no defined output: valgrind -q with an exit code (the rows are refused by the port)
                    import re
                    try:
                        vg = subprocess.run(['valgrind', '-q', '--error-exitcode=99', 'seqret'] + cmd, cwd=work, capture_output=True, timeout=120, env=dict(os.environ))
                        memerr = vg.returncode == 99 or re.search(rb'^==\d+== (Invalid|Conditional|Use of|Syscall|Mismatched|Source and destination|Process terminating)', vg.stderr, re.M) is not None
                    except subprocess.TimeoutExpired:
                        memerr = True
                if err.startswith('Read and write (return) sequences\n'): err = err[len('Read and write (return) sequences\n'):]
                rows.append({'id': len(rows), 'fmt': fmt, 'how': how + '/' + src, 'text': text, 'spec': spec, 'args': args, 'osformat': osformat, 'outseq': outseq, 'single': single, 'osext': osext, 'osdir': osdir, 'osname': osname, 'osdbname': osdbname,
                             'epoch': epoch, 'rc': rc, 'files': files, 'err': err, 'hang': hang, 'memerr': memerr, 'cmd': cmd})
                made += 1
    with gzip.open(out, 'wt', encoding='utf-8', newline='\n', compresslevel=9) as fh:
        for r in rows: fh.write(json.dumps(r, ensure_ascii=True) + '\n')
    from collections import Counter
    c = Counter((r['fmt'], 'files' if r['files'] else 'none') for r in rows)
    print(len(rows), 'rows;', sorted(c.items()))


if __name__ == '__main__':
    fm = sys.argv[4:] or FORMATS
    main(sys.argv[1], int(sys.argv[2]), int(sys.argv[3]) if len(sys.argv) > 3 else 20265711, fm)
