#!/usr/bin/env python3
"""Tests of the result checker s1b_check.py (work package 1b.9f): a clean form passes, and each of the faults that the owner's run could show, or that a doctored form would hide, makes it print FAIL:
RssAnon above the gate (in process, sampled outside), a container hash that differs from the reference or from the other run of the same text, a missing section, a verifier check that is not a pass,
a query mismatch, a damaged build, an emulated process, a test set with failures, a structure level that is "not applicable" below the reader's limit. Also what must NOT fail: the structure level and
the engine above the reader's limit, a set that was not run (INCOMPLETE, not FAIL), an unreadable file (exit code 2).
usage: python3 -I test_s1b_check.py"""
import copy
import json
import os
import subprocess
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
CHECK = os.path.join(HERE, 's1b_check.py')

CHECKS = ['files', 'meta', 'text', 'input', 'sa_permutation', 'sa_order', 'lcp_values', 'lcp_width', 'stats', 'child', 'sha256']


def run_entry(n=100_000_000, kind='genome', sha='a' * 64, structure='pass', engine='ok'):
    return {
        'kind': kind, 'n': n, 'threads': 4, 'psi': True, 'resume': True, 'fsync': True,
        'input': {'sha256': 'b' * 64, 'generator': 'e8-gen'},
        'build': {
            'exit_code': 0, 'attempts': [{'exit_code': 0, 'total_s': 30.0, 'wall_s': 31.0, 'skipped': '', 'rssanon_in_process_mb': 296.0, 'rssanon_external_mb': 297.0, 'samples': 300}],
            'phases_s': {'scan': 1.0, 'pack': 2.0, 'ph1': 3.0, 'ph23': 9.0, 'ph45': 12.0, 'bwt': 1.5, 'psi': 1.0, 'crc': 0.3, 'seal': 0.1},
            'total_s_sum': 30.0, 'wall_s_sum': 31.0, 'container_bytes': 1_348_644_448, 'scratch_peak_bytes': 1_200_000_000, 'arena_peak_bytes': 160_000_000,
            'lcp_kind': 'U16', 'overflow_pairs': 1_390_028, 'counted_read_bytes': 7_000_000_000, 'counted_write_bytes': 6_000_000_000, 'crc_hardware': 'yes'},
        'rss': {'gate_mb': 320, 'in_process_peak_mb': 296.0, 'external_peak_mb': 297.0, 'external_source': 'proc', 'verdict': 'Within'},
        'disk': {'scratch_peak_bytes_sampled': 1_190_000_000, 'total_peak_bytes_sampled': 2_700_000_000, 'free_min_bytes': 5_000_000_000},
        'container': {'bytes': 1_348_644_448, 'sha256': sha, 'reference_key': f'{kind}:{n}:psi', 'reference_sha256': None},
        'verify': {'status': 'pass', 'exit_code': 0, 'seconds': 6.2, 'all': 'pass', 'checks': {c: 'pass' for c in CHECKS}},
        'structure': {'status': structure, 'seconds': 27.0, 'rssanon_peak_mb': 109.0, 'message': ''},
        'queries': {'status': 'pass', 'total': 1500, 'sa_mismatches': 0, 'engine_status': engine, 'engine_mismatches': 0, 'set_sha256': 'c' * 64, 'oracle_present_without_hit': 0},
    }


def form(**kw):
    return {
        'form': 's1b-handoff', 'form_version': 1, 'generated_utc': '2026-11-22T00:00:00Z', 'package_commit': 'deadbeef',
        'environment': {'arch': 'aarch64', 'image_arch': 'arm64', 'process_arch': 'Arm64', 'os_arch': 'Arm64', 'emulation': 'no', 'cores': 18, 'mem_total_bytes': 8_317_202_432,
                        'cgroup_version': '2', 'fs_type': 'ext4', 'dotnet': '10.0.15', 'crc32c_hardware': True, 'crc32c_instruction': 'Crc32.Arm64.ComputeCrc32C',
                        'crc32c_check_value': 'ok', 'crc32c_table_equals_hardware': True, 'preflight_result': 'OK', 'rssanon_source': 'proc'},
        'runs': {'genome-n100000000': run_entry(**kw)},
        'tests': {'persistent': {'status': 'pass', 'passed': 672, 'failed': 0, 'skipped': 0, 'seconds': 35.0, 'configuration': 'Debug'}},
        'notes': [],
    }


def check(f, *extra, raw=None):
    with tempfile.TemporaryDirectory() as d:
        path = os.path.join(d, 'form.json')
        with open(path, 'w') as fh:
            fh.write(raw if raw is not None else json.dumps(f))
        ref = os.path.join(d, 'ref.tsv')
        with open(ref, 'w') as fh:
            fh.write('kind\tn\tpsi\tcontainer_sha256\tcontainer_bytes\tinput_sha256\tsource\n')
            fh.write('genome\t100000000\tpsi\t' + 'a' * 64 + '\t1348644448\t' + 'b' * 64 + '\ttest\n')
        p = subprocess.run([sys.executable, '-I', CHECK, path, '--reference', ref, *extra], capture_output=True, text=True)
        return p.returncode, p.stdout


def lines(out, tag):
    return [x for x in out.splitlines() if x.lstrip().startswith(tag)]


class CheckerTests(unittest.TestCase):
    def test_a_clean_form_passes(self):
        code, out = check(form())
        self.assertEqual(code, 0, out)
        self.assertIn('RESULT: PASS', out)
        self.assertEqual(lines(out, 'FAIL'), [], out)
        for gate in ('R1', 'R2', 'B1', 'I1', 'V1', 'S1', 'Q1', 'P1'):
            self.assertTrue([x for x in lines(out, 'PASS') if f' {gate} ' in x], gate + '\n' + out)

    def test_rssanon_above_the_gate_fails_in_process_and_sampled(self):
        f = form()
        f['runs']['genome-n100000000']['rss']['in_process_peak_mb'] = 321.0
        code, out = check(f)
        self.assertEqual(code, 1, out)
        self.assertTrue([x for x in lines(out, 'FAIL') if ' R1 ' in x], out)
        f = form()
        f['runs']['genome-n100000000']['rss']['external_peak_mb'] = 400.0
        code, out = check(f)
        self.assertEqual(code, 1, out)
        self.assertTrue([x for x in lines(out, 'FAIL') if ' R2 ' in x], out)

    def test_the_gate_is_the_checkers_own_not_the_forms(self):
        f = form()
        f['runs']['genome-n100000000']['rss']['gate_mb'] = 1000
        f['runs']['genome-n100000000']['rss']['external_peak_mb'] = 500.0
        code, out = check(f)
        self.assertEqual(code, 1, out)

    def test_rssanon_exactly_at_the_gate_passes_and_just_above_fails(self):
        f = form()
        f['runs']['genome-n100000000']['rss'].update(in_process_peak_mb=320.0, external_peak_mb=320.0)
        self.assertEqual(check(f)[0], 0)
        f['runs']['genome-n100000000']['rss'].update(in_process_peak_mb=320.001)
        self.assertEqual(check(f)[0], 1)

    def test_a_peak_taken_from_ps_is_not_judged(self):
        f = form()
        f['runs']['genome-n100000000']['rss'].update(external_source='ps-rss', external_peak_mb=900.0)
        code, out = check(f)
        self.assertEqual(code, 0, out)
        self.assertTrue([x for x in lines(out, 'N/A') + lines(out, 'NOT JUDGED') if ' R2 ' in x], out)

    def test_a_hash_that_differs_from_the_reference_fails(self):
        f = form(sha='d' * 64)
        code, out = check(f)
        self.assertEqual(code, 1, out)
        self.assertTrue([x for x in lines(out, 'FAIL') if ' I1 ' in x], out)

    def test_two_runs_of_one_text_must_give_one_container(self):
        f = form(n=3_000_000_000, structure='not_applicable', engine='not_applicable')
        first = next(iter(f['runs']))
        f['runs']['second'] = run_entry(n=3_000_000_000, sha='e' * 64, structure='not_applicable', engine='not_applicable')
        f['runs'][first]['container']['sha256'] = 'f' * 64
        code, out = check(f)
        self.assertEqual(code, 1, out)
        self.assertTrue([x for x in lines(out, 'FAIL') if ' I1 ' in x], out)
        f['runs']['second']['container']['sha256'] = 'f' * 64
        code, out = check(f)
        self.assertEqual(code, 0, out)
        self.assertTrue([x for x in lines(out, 'PASS') if ' I1 ' in x and 'other run' in x], out)

    def test_no_reference_at_a_new_size_is_not_a_failure(self):
        code, out = check(form(n=3_000_000_000, structure='not_applicable', engine='not_applicable'))
        self.assertEqual(code, 0, out)
        self.assertTrue([x for x in lines(out, 'N/A') + lines(out, 'INFO') if ' I1 ' in x], out)

    def test_a_missing_section_fails(self):
        for section in ('build', 'rss', 'container', 'verify', 'structure', 'queries', 'input'):
            f = form()
            del f['runs']['genome-n100000000'][section]
            code, out = check(f)
            self.assertEqual(code, 1, section + '\n' + out)
            self.assertTrue([x for x in lines(out, 'FAIL') if 'missing' in x and section in x], section + '\n' + out)

    def test_a_missing_field_inside_a_section_fails_not_crashes(self):
        f = form()
        del f['runs']['genome-n100000000']['rss']['in_process_peak_mb']
        code, out = check(f)
        self.assertEqual(code, 1, out)
        self.assertNotIn('Traceback', out)

    def test_a_verifier_check_that_is_not_a_pass_fails_and_names_it(self):
        for state in ('FAIL', 'skipped'):
            f = form()
            f['runs']['genome-n100000000']['verify']['checks']['child'] = state
            f['runs']['genome-n100000000']['verify']['all'] = 'pass'
            code, out = check(f)
            self.assertEqual(code, 1, out)
            self.assertTrue([x for x in lines(out, 'FAIL') if ' V1 ' in x and 'child' in x], out)

    def test_a_verifier_with_a_bad_exit_code_fails(self):
        f = form()
        f['runs']['genome-n100000000']['verify'].update(exit_code=1)
        self.assertEqual(check(f)[0], 1)

    def test_a_verifier_that_was_not_run_is_incomplete_not_a_pass(self):
        f = form()
        f['runs']['genome-n100000000']['verify'] = {'status': 'not_run'}
        code, out = check(f)
        self.assertEqual(code, 3, out)
        self.assertIn('INCOMPLETE', out)

    def test_a_query_mismatch_fails(self):
        f = form()
        f['runs']['genome-n100000000']['queries'].update(sa_mismatches=2)
        code, out = check(f)
        self.assertEqual(code, 1, out)
        self.assertTrue([x for x in lines(out, 'FAIL') if ' Q1 ' in x], out)
        f = form()
        f['runs']['genome-n100000000']['queries'].update(engine_mismatches=1)
        self.assertEqual(check(f)[0], 1)

    def test_a_query_set_that_is_too_small_fails(self):
        f = form()
        f['runs']['genome-n100000000']['queries'].update(total=10)
        self.assertEqual(check(f)[0], 1)

    def test_an_oracle_that_does_not_find_a_cut_query_fails(self):
        f = form()
        f['runs']['genome-n100000000']['queries'].update(oracle_present_without_hit=3)
        self.assertEqual(check(f)[0], 1)

    def test_above_the_readers_limit_the_engine_and_the_structure_level_are_not_applicable(self):
        code, out = check(form(n=3_000_000_000, structure='not_applicable', engine='not_applicable'))
        self.assertEqual(code, 0, out)
        self.assertTrue([x for x in lines(out, 'N/A') if ' S1 ' in x], out)

    def test_below_the_readers_limit_not_applicable_is_a_failure(self):
        code, out = check(form(n=2_000_000_000, structure='not_applicable'))
        self.assertEqual(code, 1, out)
        code, out = check(form(n=2_000_000_000, engine='not_applicable'))
        self.assertEqual(code, 1, out)

    def test_an_undocumented_exception_of_the_reader_is_a_warning_above_its_limit_and_a_failure_below(self):
        code, out = check(form(n=3_000_000_000, structure='error', engine='not_applicable'))
        self.assertEqual(code, 0, out)
        self.assertTrue([x for x in lines(out, 'WARN') if ' S1 ' in x], out)
        code, out = check(form(n=2_000_000_000, structure='error'))
        self.assertEqual(code, 1, out)
        self.assertTrue([x for x in lines(out, 'FAIL') if ' S1 ' in x], out)

    def test_the_structure_level_is_optional(self):
        code, out = check(form(structure='not_run'))
        self.assertEqual(code, 0, out)
        self.assertTrue([x for x in lines(out, 'NOT RUN') if ' S1 ' in x], out)

    def test_a_failed_structure_level_fails(self):
        self.assertEqual(check(form(structure='fail'))[0], 1)

    def test_a_build_that_did_not_finish_fails(self):
        f = form()
        f['runs']['genome-n100000000']['build']['exit_code'] = 137
        code, out = check(f)
        self.assertEqual(code, 1, out)
        self.assertTrue([x for x in lines(out, 'FAIL') if ' B1 ' in x], out)

    def test_a_killed_attempt_followed_by_a_resumed_one_is_fine_and_each_attempt_is_gated(self):
        f = form()
        b = f['runs']['genome-n100000000']['build']
        b['attempts'] = [{'exit_code': 137, 'total_s': 0, 'wall_s': 12.0, 'skipped': '', 'rssanon_in_process_mb': None, 'rssanon_external_mb': 250.0, 'samples': 100},
                         {'exit_code': 0, 'total_s': 20.0, 'wall_s': 21.0, 'skipped': 'packed,ranks', 'rssanon_in_process_mb': 290.0, 'rssanon_external_mb': 291.0, 'samples': 200}]
        code, out = check(f)
        self.assertEqual(code, 0, out)
        b['attempts'][0]['rssanon_external_mb'] = 330.0
        self.assertEqual(check(f)[0], 1, 'a killed attempt that was over the gate is still over the gate')

    def test_the_container_size_must_equal_what_the_build_reported(self):
        f = form()
        f['runs']['genome-n100000000']['container']['bytes'] = 1_000
        self.assertEqual(check(f)[0], 1)

    def test_an_emulated_process_fails(self):
        f = form()
        f['environment'].update(emulation='yes: an x86-64 process on an ARM CPU')
        code, out = check(f)
        self.assertEqual(code, 1, out)
        self.assertTrue([x for x in lines(out, 'FAIL') if ' P1 ' in x], out)

    def test_a_native_run_without_an_image_is_not_an_emulation(self):
        f = form()
        f['environment'].update(image_arch='unknown', arch='x86_64', process_arch='X64', os_arch='X64')
        code, out = check(f)
        self.assertEqual(code, 0, out)
        f['environment'].update(image_arch='arm64')
        self.assertEqual(check(f)[0], 1, 'an arm64 image on an x86-64 kernel is an emulation')

    def test_an_arm64_cpu_without_the_hardware_crc_is_a_warning_and_a_table_path_that_differs_is_a_failure(self):
        f = form()
        f['environment'].update(crc32c_hardware=False)
        code, out = check(f)
        self.assertEqual(code, 0, out)
        self.assertTrue([x for x in lines(out, 'WARN') if ' X1 ' in x], out)
        f['environment'].update(crc32c_table_equals_hardware=False)
        self.assertEqual(check(f)[0], 1)
        f = form()
        f['environment'].update(crc32c_check_value='WRONG')
        self.assertEqual(check(f)[0], 1)

    def test_a_failed_test_set_fails_and_a_different_count_is_a_warning(self):
        f = form()
        f['tests']['persistent'].update(failed=1, passed=671)
        code, out = check(f)
        self.assertEqual(code, 1, out)
        self.assertTrue([x for x in lines(out, 'FAIL') if ' T1 ' in x and 'persistent' in x], out)
        f = form()
        f['tests']['persistent'].update(passed=600)
        code, out = check(f)
        self.assertEqual(code, 0, out)
        self.assertTrue([x for x in lines(out, 'WARN') if ' T1 ' in x], out)

    def test_a_filtered_test_run_is_a_warning_not_a_pass_of_the_set(self):
        f = form()
        f['tests']['persistent'].update(filter='FullyQualifiedName~Handoff', passed=20, skipped=0)
        code, out = check(f)
        self.assertEqual(code, 0, out)
        self.assertTrue([x for x in lines(out, 'WARN') if ' T1 ' in x and 'filtered' in x], out)
        f['tests']['persistent'].update(failed=1)
        self.assertEqual(check(f)[0], 1)

    def test_a_form_with_no_tests_is_not_a_failure_unless_they_are_required(self):
        f = form()
        f['tests'] = {}
        self.assertEqual(check(f)[0], 0)
        code, out = check(f, '--require-tests')
        self.assertEqual(code, 3, out)

    def test_a_form_that_is_not_json_or_not_a_form_exits_with_two(self):
        self.assertEqual(check(None, raw='{not json')[0], 2)
        self.assertEqual(check(None, raw='[]')[0], 2)
        self.assertEqual(check({'form': 'something-else'})[0], 2)
        f = form()
        f['form_version'] = 99
        self.assertEqual(check(f)[0], 2)

    def test_run_selects_the_runs_that_are_judged_and_the_others_still_serve_as_references(self):
        f = form()
        f['runs']['second'] = run_entry(sha='a' * 64)
        f['runs']['second']['rss']['in_process_peak_mb'] = 500.0
        self.assertEqual(check(f)[0], 1, 'without --run every run is judged')
        code, out = check(f, '--run', 'genome-n100000000')
        self.assertEqual(code, 0, out)
        self.assertNotIn('run second (', out)
        self.assertEqual(check(f, '--run', 'second')[0], 1)
        code, out = check(f, '--run', 'nonexistent')
        self.assertEqual(code, 2, out)
        f['runs']['second']['container']['sha256'] = 'f' * 64
        f['runs']['second']['rss']['in_process_peak_mb'] = 100.0
        code, out = check(f, '--run', 'genome-n100000000')
        self.assertEqual(code, 1, 'a run that is judged is compared with the other runs of the same text, judged or not')
        self.assertTrue([x for x in lines(out, 'FAIL') if ' I1 ' in x], out)

    def test_the_gate_numbers_of_the_checker_are_the_ones_of_the_owner_decision(self):
        out = subprocess.run([sys.executable, '-I', CHECK, '--gates'], capture_output=True, text=True).stdout
        self.assertIn('RSSANON_GATE_MB = 320', out)


if __name__ == '__main__':
    unittest.main(verbosity=1)
