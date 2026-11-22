#!/usr/bin/env python3
"""A5-8c3h (B09 F258, WP58) oracle: EMBOSS 6.6.0 `extractfeat` (emboss/extractfeat.c: the sequence of every feature that matches the options)
    extractfeat -auto -sequence in.dat [-sformat1 F] [-before N] [-after N] [-join] [-featinname] [-describe TAGS] [-source S] [-type T] [-sense N] [-minscore X] [-maxscore X] [-tag T] [-value V]
                -outseq o.out -osformat2 OF
Inputs: entries with feature tables in the formats whose readers carry the features (embl, genbank, swiss, nbrf, gff3, refseqp, gff2) made as in generate_seqret_oracle.py, FASTA files without features, damaged files.
JSONL.GZ rows: id, fmt, how, text, spec, args, epoch, rc, out, err, hang, memerr
Usage: python3 scripts/alignment/generate_extractfeat_oracle.py OUT.jsonl.gz [n=1000] [seed]      (needs EMBOSS 6.6.0, valgrind)
"""
import gzip, json, os, random, re, shutil, subprocess, sys, tempfile

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import generate_seqret_oracle as S  # noqa: E402

FMTS = ['embl', 'genbank', 'swiss', 'nbrf', 'gff3', 'refseqp', 'gff2', 'embl', 'genbank', 'gff3', 'swiss', 'embl', 'genbank', 'embl', 'genbank', 'gff3', 'fasta']
TYPES = ['repeat*', 'gene', 'CDS', 'mRNA', 'exon', 'intron', 'misc_feature', 'source', '*UTR', 'repeat_region', 'tRNA rRNA', 'gene,CDS', 'exon;intron', 'gene|CDS', 'misc*', 'g?ne', '*', 'binding_site', 'variation',
         'polypeptide_region', 'polypeptide_domain', 'domain', 'region', 'Region', 'CDS gene mRNA', 'x', 'repeat_unit', 'databank_entry', '*_region', 'sequence_feature', '']
SOURCES = ['*', 'EMBL', 'embl', 'GenBank', 'genbank', 'gene*', 'embl|gff', 'SWISSPROT', 'swiss*', 'pir', 'x', 'GFF*', 'gff3', 'gff*|embl', '']
TAGS = ['*', 'gene', 'note', 'product', 'label', 'db_xref', 'gene|label', 'number', 'codon_start', 'translation', 'x', '', 'ge*', '*note*', 'organism', 'mol_type']
VALUES = ['*', 'a*', 'x', 'pax*', '1', '10', '*1', 'abc', 'a|b', 'a b', 'gene*', 'hypothetical*', '*protein', '?']
DESCRIBE = ['gene', 'note', 'product', 'gene|label', 'gene note', '*', 'x', 'db_xref', 'label', 'ge*', 'codon_start', 'number', '']
OSFORMATS = ['fasta', 'fasta', 'fasta', 'fasta', 'gcg', 'gcg8', 'pir', 'nbrf', 'raw', 'text', 'plain', 'staden', 'ig', 'iguspto', 'codata', 'fitch', 'mase', 'strider', 'gde', 'acedb', 'ncbi', 'gifasta', 'dasdna', 'das', 'asn1', 'nib', 'nibble', 'unknown', 'pearson', 'embl', 'em', 'emblnew', 'emold', 'emblold', 'genbank', 'gb', 'ddbj', 'refseq', 'genpept', 'refseqp', 'swiss', 'sw', 'swissprot', 'uniprot', 'swissnew', 'swold', 'fastq', 'fastq-sanger', 'fastq-illumina', 'fastq-solexa', 'sam', 'gff', 'gff2', 'gff3', 'scf', 'msf', 'clustal', 'aln', 'selex', 'phylip', 'phylipnon', 'phylip3', 'hennig86', 'mega', 'meganon', 'nexus', 'paup', 'nexusnon', 'paupnon', 'jackknifer', 'jackknifernon', 'treecon']   # the writers ported so far (every writer but experiment, debug and bam: B09 WP59)


def main(out, n, seed):
    rnd = random.Random(seed)
    rows = []
    with tempfile.TemporaryDirectory() as d:
        work = os.path.join(d, 'w'); os.mkdir(work)
        gen = os.path.join(d, 'g'); os.mkdir(gen)
        shim = os.path.join(d, 'faketime.so')
        subprocess.run(['gcc', '-shared', '-fPIC', '-o', shim, os.path.join(os.path.dirname(os.path.abspath(__file__)), 'faketime.c')], check=True)
        while len(rows) < n:
            fmt = FMTS[len(rows) % len(FMTS)]
            text, how = S.make_input(rnd, gen, fmt)
            if rnd.random() < 0.12: text, how = S.damaged(rnd, text), how + '+damaged'
            spec = 'in.dat'
            args = ['-sformat1', fmt] if rnd.random() < 0.5 else []
            if rnd.random() < 0.3: args += ['-type', rnd.choice(TYPES)]
            if rnd.random() < 0.3: args += ['-before', str(rnd.choice([0, 1, 3, 10, 25, -2, -5, 100]))]
            if rnd.random() < 0.3: args += ['-after', str(rnd.choice([0, 1, 3, 10, 25, -2, -5, 100]))]
            if rnd.random() < 0.3: args.append('-join')
            if rnd.random() < 0.2: args.append('-featinname')
            if rnd.random() < 0.3: args += ['-describe', rnd.choice(DESCRIBE)]
            if rnd.random() < 0.07: args += ['-source', rnd.choice(SOURCES)]
            if rnd.random() < 0.15: args += ['-sense', str(rnd.choice([0, 1, -1, 2]))]
            if rnd.random() < 0.05: args += ['-minscore', str(rnd.choice([0, 1, 0.5, 100, -1]))]
            if rnd.random() < 0.05: args += ['-maxscore', str(rnd.choice([0, 1, 0.5, 100, -1, 1000]))]
            if rnd.random() < 0.1: args += ['-tag', rnd.choice(TAGS)]
            if rnd.random() < 0.06: args += ['-value', rnd.choice(VALUES)]
            if rnd.random() < 0.05: args += ['-supper1' if rnd.random() < 0.5 else '-slower1']
            if rnd.random() < 0.05: args += ['-sbegin1', str(rnd.randint(1, 20))]
            if rnd.random() < 0.05: args += ['-firstonly'] if False else []
            osf = rnd.choice(OSFORMATS)
            for nme in os.listdir(work): os.remove(os.path.join(work, nme))
            open(os.path.join(work, 'in.dat'), 'w', newline='\n', encoding='latin-1').write(text)
            cmd = ['-auto', '-sequence', spec] + args + ['-outseq', 'o.out', '-osformat2', osf]
            hang = False
            epoch = rnd.randint(978307200, 2051222400)
            env = dict(os.environ, LD_PRELOAD=shim, FAKE_EPOCH=str(epoch))
            try:
                import resource
                def limit():
                    resource.setrlimit(resource.RLIMIT_FSIZE, (4 << 20, 4 << 20))
                p = subprocess.run(['extractfeat'] + cmd, cwd=work, capture_output=True, timeout=5, preexec_fn=limit, env=env)
                rc, err = p.returncode, p.stderr.decode('latin-1')
                if rc == -25: hang = True
            except subprocess.TimeoutExpired:
                rc, err, hang = -9, '', True
                subprocess.run(['pkill', '-9', 'extractfeat'])
            q = os.path.join(work, 'o.out')
            outtext = open(q, encoding='latin-1', newline='').read() if os.path.exists(q) and not hang else ''
            if err.startswith('Extract features from sequence(s)\n'): err = err[len('Extract features from sequence(s)\n'):]
            if not outtext and not hang and rnd.random() < 0.65:
                continue                                    # keep the share of runs without output small
            memerr = False
            if not hang and shutil.which('valgrind') and (rc != 0 or len(rows) % 6 == 0):
                vgargs = [a.replace('o.out', 'ovg.out') for a in cmd]
                try:
                    vg = subprocess.run(['valgrind', '-q', '--error-exitcode=99', 'extractfeat'] + vgargs, cwd=work, capture_output=True, timeout=120)
                    memerr = vg.returncode == 99 or re.search(rb'^==\d+== (Invalid|Conditional|Use of|Syscall|Mismatched|Source and destination|Process terminating)', vg.stderr, re.M) is not None
                except subprocess.TimeoutExpired:
                    memerr = True
            rows.append({'id': len(rows), 'fmt': fmt, 'how': how, 'text': text, 'spec': spec, 'args': args + ['-osformat2', osf], 'epoch': epoch, 'rc': rc, 'out': outtext, 'err': err, 'hang': hang, 'memerr': memerr})
    with gzip.open(out, 'wt', encoding='utf-8', newline='\n', compresslevel=9) as fh:
        for r in rows: fh.write(json.dumps(r, ensure_ascii=True) + '\n')
    from collections import Counter
    print(len(rows), 'rows;', sorted(Counter((r['fmt'], 'hang' if r['hang'] else 'out' if r['out'] else 'none') for r in rows).items()))


if __name__ == '__main__':
    main(sys.argv[1], int(sys.argv[2]) if len(sys.argv) > 2 else 1000, int(sys.argv[3]) if len(sys.argv) > 3 else 20265804)
