#!/usr/bin/env python3
"""Tests of the streaming generator s1b_gen.py (work package 1b.9f): the text of every kind is byte-identical to what the benchmark's own generators make in memory (gen_inputs.py for the random DNA and
the repetitive text, make_1b9e_inputs.py for the protein-like text), for any chunk size; the sums file says what was written; a file is renamed into place only when complete; an existing file is not
overwritten without --force. usage: python3 -I test_s1b_gen.py"""
import hashlib
import importlib.util
import os
import subprocess
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
GEN = os.path.join(HERE, 's1b_gen.py')
EM = os.path.join(HERE, '..', 'em')
BENCH_GEN = os.path.join(HERE, '..', 'bench', 'gen')


def load(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def generate(kind, n, extra=()):
    d = tempfile.mkdtemp(prefix='s1bgen')
    out = os.path.join(d, f'{kind}.txt')
    p = subprocess.run([sys.executable, '-I', GEN, '--kind', kind, '--n', str(n), '--out', out, *extra], capture_output=True, text=True)
    return p, out


def read(path):
    with open(path, 'rb') as fh:
        return fh.read()


class GeneratorTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.gi = load(os.path.join(BENCH_GEN, 'gen_inputs.py'), 'gi')
        cls.em = load(os.path.join(EM, 'make_1b9e_inputs.py'), 'em')

    def test_dna_is_the_prefix_of_the_benchmark_stream(self):
        for n in (1, 3, 4, 5, 1000, 2_000_003):
            p, out = generate('dna', n)
            self.assertEqual(p.returncode, 0, p.stderr)
            # the first n symbols of the word stream: for a multiple of four this is rand_dna(rng, n) itself, otherwise the prefix of the next multiple (as dna_10m is the prefix of the 50 Mbp stream)
            self.assertEqual(read(out), self.gi.rand_dna(self.gi.rng('dna50'), n + (-n) % 4)[:n], n)
            if n % 4 == 0:
                self.assertEqual(read(out), self.gi.rand_dna(self.gi.rng('dna50'), n))

    def test_dna_does_not_depend_on_the_chunk_size(self):
        a = generate('dna', 1_000_001, ('--chunk-mb', '1'))[1]
        b = generate('dna', 1_000_001, ('--chunk-mb', '4'))[1]
        self.assertEqual(read(a), read(b))

    def test_rep_is_the_repetitive_text_of_the_benchmark(self):
        for n in (1, 99_999, 100_000, 250_017):
            p, out = generate('rep', n)
            self.assertEqual(p.returncode, 0, p.stderr)
            rr = self.gi.rng('rep100')
            block = self.gi.rand_dna(rr, 100_000)
            expected = b''.join(self.gi.mutate(rr, block, 0.01, self.gi.DNA) for _ in range(-(-n // 100_000)))[:n]
            self.assertEqual(read(out), expected, n)

    def test_prot_is_the_protein_like_text_of_1b9e(self):
        for n in (1, 5000, 123_457, 600_000):
            p, out = generate('prot', n)
            self.assertEqual(p.returncode, 0, p.stderr)
            self.assertEqual(read(out), self.em.protein_like(n), n)

    def test_the_sums_file_describes_the_text(self):
        for kind in ('dna', 'rep', 'prot'):
            p, out = generate(kind, 30_000)
            data = read(out)
            sums = dict(line.split('=', 1) for line in read(out + '.sums').decode().splitlines())
            self.assertEqual(sums['n'], '30000')
            self.assertEqual(sums['sha256'], hashlib.sha256(data).hexdigest())
            self.assertEqual(sums['kind'], kind)
            letters = sums['alphabet']
            self.assertEqual(sorted(set(data.decode())), sorted(letters))
            self.assertEqual(sums['counts'], ','.join(str(data.count(c.encode())) for c in letters))

    def test_an_existing_file_is_kept_without_force_and_replaced_with_it(self):
        p, out = generate('dna', 100)
        q = subprocess.run([sys.executable, '-I', GEN, '--kind', 'dna', '--n', '200', '--out', out], capture_output=True, text=True)
        self.assertNotEqual(q.returncode, 0)
        self.assertEqual(len(read(out)), 100)
        q = subprocess.run([sys.executable, '-I', GEN, '--kind', 'dna', '--n', '200', '--out', out, '--force'], capture_output=True, text=True)
        self.assertEqual(q.returncode, 0, q.stderr)
        self.assertEqual(len(read(out)), 200)
        self.assertFalse(os.path.exists(out + '.part'))

    def test_a_bad_kind_or_size_is_refused(self):
        self.assertNotEqual(generate('nonsense', 10)[0].returncode, 0)
        self.assertNotEqual(generate('dna', 0)[0].returncode, 0)

    def test_the_size_accepts_decimal_suffixes(self):
        p, out = generate('dna', '2k')
        self.assertEqual(len(read(out)), 2000)


if __name__ == '__main__':
    unittest.main(verbosity=1)
