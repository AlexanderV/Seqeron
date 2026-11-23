#!/usr/bin/env python3
"""Deterministic cases for Bio.Align.AlignInfo.SummaryInfo of Biopython 1.79 (B09 WP67, R14-13): writes cases.jsonl, one JSON object per line:
{"id", "records": [{"id", "seq", "weight" (optional)}], "calls": [{"method", "kwargs"}]}.
Python 3.9 standard library only. Usage: python3 -I make_cases.py OUT.jsonl"""
import json, random, sys

AA = 'ACDEFGHIKLMNPQRSTVWY'


def mutate(rnd, base, alphabet, rate, gap_rate):
    out = []
    for c in base:
        r = rnd.random()
        if r < gap_rate:
            out.append('-')
        elif r < gap_rate + rate:
            out.append(rnd.choice(alphabet))
        else:
            out.append(c)
    return ''.join(out)


def uniform(letters):
    return {c: 1.0 / len(letters) for c in letters}


def calls_for(rnd, seqs, alphabet):
    letters = ''.join(sorted(set(''.join(seqs)) - {'-'}))
    n = len(seqs[0])
    calls = []
    for thr in (0.7, 0.5, 1.0, 0.0, 0.34):
        for amb in ('X', 'N'):
            calls.append({'method': 'dumb_consensus', 'kwargs': {'threshold': thr, 'ambiguous': amb}})
            calls.append({'method': 'gap_consensus', 'kwargs': {'threshold': thr, 'ambiguous': amb}})
    calls.append({'method': 'dumb_consensus', 'kwargs': {'require_multiple': True}})
    calls.append({'method': 'gap_consensus', 'kwargs': {'require_multiple': True, 'ambiguous': '?'}})
    calls.append({'method': 'replacement_dictionary', 'kwargs': {}})
    calls.append({'method': 'replacement_dictionary', 'kwargs': {'letters': letters}})
    calls.append({'method': 'replacement_dictionary', 'kwargs': {'letters': letters[:3]}})
    calls.append({'method': 'replacement_dictionary', 'kwargs': {'letters': letters + '-'}})
    calls.append({'method': 'replacement_dictionary', 'kwargs': {'letters': letters, 'skip_chars': ['-']}})
    for ignore in (None, [], ['-'], ['-', '*'], ['A'], 'notalist'):
        kw = {}
        if ignore is not None:
            kw['chars_to_ignore'] = ignore
        calls.append({'method': 'pos_specific_score_matrix', 'kwargs': kw})
        calls.append({'method': 'pos_specific_score_matrix', 'kwargs': dict(kw, axis_seq='dumb')})
    calls.append({'method': 'pos_specific_score_matrix', 'kwargs': {'axis_seq': 'K' * n, 'chars_to_ignore': ['-']}})
    calls.append({'method': 'pos_specific_score_matrix', 'kwargs': {'axis_seq': 'K' * (n + 1)}})
    calls.append({'method': 'information_content', 'kwargs': {}})
    for table in (uniform(letters) if letters else {}, uniform(letters + 'Z') if letters else {}, {k: v for k, v in list(uniform(letters).items())[:-1]} if letters else {}, {}):
        for ignore in ([], ['-'], ['-', '*']):
            for base in (2, 2.718281828459045, 10):
                for pc in (0, 1, 0.5):
                    if rnd.random() < 0.35:
                        calls.append({'method': 'information_content', 'kwargs': {'e_freq_table': table, 'chars_to_ignore': ignore, 'log_base': base, 'pseudo_count': pc}})
    t = uniform(letters) if letters else {}
    calls.append({'method': 'information_content', 'kwargs': {'e_freq_table': t, 'chars_to_ignore': ['-'], 'start': 1, 'end': max(2, n - 1)}})
    calls.append({'method': 'information_content', 'kwargs': {'e_freq_table': t, 'chars_to_ignore': ['-'], 'start': n}})
    calls.append({'method': 'information_content', 'kwargs': {'e_freq_table': t, 'chars_to_ignore': ['-'], 'start': -1}})
    calls.append({'method': 'information_content', 'kwargs': {'e_freq_table': t, 'chars_to_ignore': ['-'], 'end': n + 1}})
    calls.append({'method': 'information_content', 'kwargs': {'e_freq_table': t, 'chars_to_ignore': ['-'], 'pseudo_count': -1}})
    calls.append({'method': 'information_content', 'kwargs': {'e_freq_table': t, 'chars_to_ignore': ['-'], 'log_base': 1}})
    calls.append({'method': 'information_content', 'kwargs': {'e_freq_table': t, 'chars_to_ignore': ['-'], 'log_base': 0}})
    if letters:
        calls.append({'method': 'information_content', 'kwargs': {'e_freq_table': dict(t, **{letters[0]: 0.0}), 'chars_to_ignore': ['-']}})
    calls.append({'method': 'ic_vector_after', 'kwargs': {'e_freq_table': t, 'chars_to_ignore': ['-'], 'start': 1, 'end': max(2, n - 1)}})
    calls.append({'method': 'print_info_content', 'kwargs': {'e_freq_table': t, 'chars_to_ignore': ['-']}})
    calls.append({'method': 'print_info_content', 'kwargs': {}})
    for col in (0, n // 2, n - 1, n, -1):
        calls.append({'method': 'get_column', 'kwargs': {'col': col}})
    return calls


def main():
    out = sys.argv[1]
    rnd = random.Random(20262503)
    cases = []
    # the cases of Tests/test_AlignInfo.py of Biopython 1.79
    cases.append({'id': 'bp-nucleotides', 'records': [{'id': 'test%d' % (i + 1), 'seq': s} for i, s in enumerate(['ACGTCGCG', 'GGGGCCCC', 'AAACACAC'])]})      # GFF/multi.fna of the Biopython tests (ic 7.32029999423075 with the uniform 0.25 table)
    cases.append({'id': 'bp-proteins', 'records': [{'id': 'ID001', 'seq': 'MHQAIFIYQIGYP*LKSGYIQSIRSPEYDNW-'}, {'id': 'ID002', 'seq': 'MH--IFIYQIGYAYLKSGYIQSIRSPEY-NW*'}, {'id': 'ID003', 'seq': 'MHQAIFIYQIGYPYLKSGYIQSIRSPEYDNW*'}]})
    cases.append({'id': 'bp-pseudo', 'records': [{'id': 'ID00%d' % (i + 1), 'seq': s} for i, s in enumerate(['AACCACGTTTAA', 'CACCACGTGGGT', 'CACCACGTTCGC', 'GCGCACGTGGGG', 'TCGCACGTTGTG', 'TGGCACGTGTTT', 'TGACACGTGGGA', 'TTACACGTGCGC'])]})
    cases.append({'id': 'bp-weights', 'records': [{'id': 'a', 'seq': 'GTATC', 'weight': 0.5}, {'id': 'b', 'seq': 'AT--C', 'weight': 0.8}, {'id': 'c', 'seq': 'CTGTC', 'weight': 1.0}]})
    for k in range(60):
        kind = rnd.choice(['n', 'p'])
        alphabet = 'ACGT' if kind == 'n' else AA
        nseq = rnd.choice([1, 2, 3, 4, 6])
        length = rnd.choice([1, 3, 6, 10, 25])
        base = ''.join(rnd.choice(alphabet) for _ in range(length))
        seqs = [mutate(rnd, base, alphabet + ('*' if kind == 'p' and rnd.random() < 0.3 else ''), rnd.choice([0.05, 0.3, 0.6]), rnd.choice([0, 0.05, 0.25])) for _ in range(nseq)]
        recs = [{'id': 's%d' % i, 'seq': s} for i, s in enumerate(seqs)]
        if rnd.random() < 0.4:
            for r in recs:
                r['weight'] = rnd.choice([0.25, 0.5, 0.75, 1.0, 1.5, 2.0])
        cases.append({'id': 'rnd-%02d' % k, 'records': recs})
    with open(out, 'w') as fh:
        for c in cases:
            c['calls'] = calls_for(rnd, [r['seq'] for r in c['records']], None)
            fh.write(json.dumps(c) + '\n')
    print(len(cases), 'cases', sum(len(c['calls']) for c in cases), 'calls')


if __name__ == '__main__':
    main()
