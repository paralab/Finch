import math
import os
import pickle
import shutil
import subprocess
import tarfile
import tempfile
import time
import typing
import struct
import argparse
import sys
import zlib
from fnmatch import fnmatch

from gen_folders import gen_folders, get_varied_keys
from file_generator import FileGenerator, FileCopyGenerator, FileLinkGenerator, FolderCopyGenerator


class ReferenceFileMissing:
    def __init__(self, src_test, path: str):
        if src_test and hasattr(src_test, 'name'):
            self.msg = '{} - reference file missing ({})'.format(src_test.name, path)
        else:
            self.msg = "Reference file missing ({})".format(path)

    def __str__(self):
        return self.msg


class ResultFileMissing:
    def __init__(self, src_test, path: str):
        if src_test and hasattr(src_test, 'name'):
            self.msg = '{} - result file missing ({})'.format(src_test.name, path)
        else:
            self.msg = "Result file missing ({})".format(path)

    def __str__(self):
        return self.msg


class TestPassed:
    def __init__(self, src_test, msg: str):
        prefix = 'Ok'
        if src_test and hasattr(src_test, 'name'):
            prefix = src_test.name + ' ok'
        self.msg = "{} - {}".format(prefix, msg)

    def __str__(self):
        return self.msg


class TestFailed:
    def __init__(self, src_test, msg: str):
        prefix = 'Failed'
        if src_test and hasattr(src_test, 'name'):
            prefix = src_test.name + ' failed'
        self.msg = "{} - {}".format(prefix, msg)

    def __str__(self):
        return self.msg


def _get_precision(tol: float):
    """
    :return: number of figures required to represent tol in decimal notation
    """
    return int(math.ceil(-math.log10(tol))) + 1


class RegexDiffMetric:
    def __init__(self, test_name: str, path: str, exp: typing.Pattern, tol: float = 1e-12, required: bool = True):
        """
        Parse the file at path and pull out lines with exp.
        Checks that the fields of exp match within tol.
        Note this only supports non-overlapping regex matches.
        :param test_name: name of test, printed in results
        :param path: path to file
        :param exp: regular expression
        :param tol: tolerance
        :param required: require at least one match (if no matches found in either file, give error)
        """
        self.name = test_name
        self.path = path
        self.exp = exp
        self.tol = tol
        self.required = required

    def files(self):
        return [self.path]

    def test(self, ref_root, check_root):
        ref_path = os.path.join(ref_root, self.path)
        check_path = os.path.join(check_root, self.path)
        if not os.path.exists(ref_path):
            return ReferenceFileMissing(self, ref_path)
        if not os.path.exists(check_path):
            return ResultFileMissing(self, check_path)

        # note: reads entire file into memory, may fail for giant files
        with open(ref_path, 'r') as ref, open(check_path, 'r') as check:
            ref_matches = list(self.exp.finditer(ref.read()))
            check_matches = list(self.exp.finditer(check.read()))

        if len(ref_matches) == len(check_matches) and len(ref_matches) == 0 and self.required:
            return TestFailed(self, "pattern '{}' not found in reference file".format(self.exp.pattern))

        # if len(ref_matches) != len(check_matches), zip stops after the shorter iterable is exhausted
        max_error = -math.inf
        precision = _get_precision(self.tol)  # how many decimals to print
        for ref, check in zip(ref_matches, check_matches):
            for ref_val, check_val in zip(ref.groups(), check.groups()):
                ref_val = float(ref_val)
                check_val = float(check_val)
                err = abs(ref_val - check_val)
                if err > self.tol:
                    return TestFailed(self, "reference {:.{prec}f} vs result {:.{prec}f} exceeds tol {} ('{}' vs '{}')"
                                      .format(ref_val, check_val, self.tol, ref.groups(0), check.groups(0), prec=precision))
                max_error = max(max_error, err)

        # intentionally after zip, as differences between what matches there are may explain why
        # the total number of matches differs
        if len(ref_matches) != len(check_matches):
            return TestFailed(self, "different number of matches (expected {} matches, found {})"
                              .format(len(ref_matches), len(check_matches)))

        return TestPassed(self, "max difference = {:.{prec}f}".format(max_error, prec=precision))


class VecDiffMetric:
    def __init__(self, path: str, tol: float = 1e-12):
        """
        Parse the PETSc vector dumped to path, then calculates the L2 norm.
        Checks that the norm is within tol.
        NOTE: Assumes PetscScalar is double (8 bytes) and PetscInt is int (4 bytes).
              Will break if PETSc is compiled with the --with-64bit-indices option.
        :param path: path to petsc vector file
        :param tol: tolerance
        """
        self.name = os.path.splitext(path)[0]
        self.path = path
        self.tol = tol

    def files(self):
        return [self.path]

    def test(self, ref_root, check_root):
        ref_path = os.path.join(ref_root, self.path)
        check_path = os.path.join(check_root, self.path)
        if not os.path.exists(ref_path):
            return ReferenceFileMissing(self, ref_path)
        if not os.path.exists(check_path):
            return ResultFileMissing(self, check_path)

        max_error = -math.inf
        precision = _get_precision(self.tol)  # how many decimals to print
        with open(ref_path, 'rb') as ref, open(check_path, 'rb') as check:
            ref_vec_class, ref_n_rows = struct.unpack('>ii', ref.read(8))
            check_vec_class, check_n_rows = struct.unpack('>ii', check.read(8))

            if ref_n_rows != check_n_rows:
                return TestFailed(self, "different number of rows in vector (expected {}, got {})"
                                  .format(ref_n_rows, check_n_rows))

            failed_rows = []
            for row in range(ref_n_rows):
                ref_val = struct.unpack('>d', ref.read(8))[0]
                check_val = struct.unpack('>d', check.read(8))[0]
                err = abs(ref_val - check_val)
                max_error = max(max_error, err)
                if err > self.tol:
                    failed_rows.append((row, ref_val, check_val, err))

            if len(failed_rows) > 0:
                if len(failed_rows) < 20:
                    msg = 'Found {} rows exceeding tolerance {}:\n'.format(len(failed_rows), self.tol)
                    row_msg = "  on row {}, reference {:.{prec}f} vs result {:.{prec}f}, err = {:.{prec}f}"
                    msg += '\n'.join([row_msg.format(r[0], r[1], r[2], r[3], prec=precision) for r in failed_rows])
                    return TestFailed(self, msg)
                else:
                    min_row = min([r[0] for r in failed_rows])
                    max_row = max([r[0] for r in failed_rows])
                    return TestFailed(self, "found differences on {} rows (ranging from {} to {}), max difference was {}"
                                      .format(len(failed_rows), min_row, max_row, max_error))

        return TestPassed(self, "max difference = {:.{prec}f}".format(max_error, prec=precision))


def _filter_cases(cases, filters):
    """
    filter cases by names
    :param cases:
    :param filters: Todo Only filter by name at this moment
    :return: list of filtered cases
    """
    if filters is None:
        return list(cases)
    else:
        filtered_case = []
        for case in cases:
            for filter in filters:
                if case.get("name") and case.get("name") in filter:
                    filtered_case.append(case)
                    break
    return filtered_case


def _case_hash(case: dict) -> int:
    return zlib.adler32(repr(tuple(sorted(case.items()))).encode())


class ResultsArchive:
    def __init__(self, path: str):
        """
        :param path: base directory for storing results
        """
        self.path = path

        if not os.path.exists(path):
            os.makedirs(self.path)

    def _case_path(self, case: dict) -> str:
        """
        :param case: case configuration
        :return: archive path for reference data for case (may or may not exist)
        """
        return os.path.join(self.path, str(_case_hash(case)) + '.tar.gz')

    def add(self, case: dict, add_dir: str, exclude_patterns: typing.List[str] = list()) -> None:
        """
        Compress add_dir as _case_path(case) (a tgz archive), overwriting it if it already exists,
        excluding any files that match any of the Unix shell-style patterns in exclude_patterns.
        It can be retrieved later through extract_reference().
        :param case: case configuration
        :param add_dir: directory to store (should already have results in it)
        :param exclude_patterns: list of Unix shell-style patterns for files to exclude from the archive
        """
        with tarfile.open(self._case_path(case), mode='w:gz') as tf:
            for root, dirs, files in os.walk(add_dir):
                for file in files:
                    if any([fnmatch(file, e) for e in exclude_patterns]):
                        continue

                    tf.add(os.path.join(root, file), arcname=file, recursive=False)

    def has_reference(self, case: dict) -> bool:
        """
        :param case: case configuration
        :return: True if case has been add'd before
        """
        return os.path.exists(self._case_path(case))

    def extract_reference(self, case: dict):
        """
        Extract previously add'd data for case to a temporary directory.
        Multiple references may be extracted simultaneously.
        :param case: case configuration
        :return: A tempfile context manager for the extracted file (use it with the 'with' keyword).
                 If you use 'with', the temporary directory will be automatically deleted at
                 the end of the with's scope.
        """
        tempdir = tempfile.TemporaryDirectory()
        with tarfile.open(self._case_path(case), mode='r') as f:
            f.extractall(tempdir.name)
        return tempdir

    def remove_all_except(self, cases: typing.List[dict]) -> None:
        """
        Removes all files in self.path that do not match "_case_hash(cases).tar.gz".
        Used for cleaning stale reference data (i.e. after cases gains a new key, which changes all hashes).
        :param cases: list of configurations
        :return:
        """
        hashes = [str(_case_hash(c)) for c in cases]
        for root, dirs, files in os.walk(self.path):
            for f in files:
                name = os.path.basename(f)[:-7]
                if not name in hashes:
                    print("Removing", f)
                    os.remove(os.path.join(root, f))


class RegressionTester:
    def __init__(self):
        self.cases = []  # type: typing.List[dict]
        self.generators = []  # type: typing.List[typing.Callable[[dict, str], None], list]
        self.metrics = []
        self.run_cmd = None  # type: typing.List[str]
        self.run_dir = 'regression_tests'
        self.results_archive = ResultsArchive('regression_reference_data')
        self.consistent_keys = ['ntasks', 'mfree']
        self.results_exclude_patterns = []  # type: typing.List[str]

    def add_metric(self, metric):
        """
        Add a metric used to compare cases.
        :param metric: should implement test(ref_root: str, check_root: str)
        :return:
        """
        self.metrics.append(metric)

    def add_generator(self, generator, condition=None) -> None:
        self.generators.append([generator, condition])

    def add_file_generator(self, filename: str, template: str) -> None:
        self.add_generator(FileGenerator(filename, template))

    def set_run_cmd(self, run_cmd: str) -> None:
        """
        String for the program run command. May contain placeholders for case parameters.
        (e.g. {refine_lvl} will be replaced by case['refine_lvl'])
        Should not contain 'mpirun' or similar (the runner will be automatically prefixed to the command for you).
        :param run_cmd: string for program run command, may contain format placeholders matching case keys
        """
        self.run_cmd = run_cmd

    def add_cases(self, base_cfg, cases: typing.List[dict]=list()) -> None:
        """
        Add cases, using base_cfg as a common base, then merging base_cfg + cases[i] (i.e. keys in cases
        overwrite or append to base_cfg). May be called multiple times to add multiple sets of cases.
        :param base_cfg: common configuration parameters that all cases should start with
        :param cases: list of cases, which are added to (or overwrite) base_cfg to get the final list of cases
        """
        if len(cases) == 0:
            cases = [{}]

        for case in cases:
            case = {**base_cfg, **case}
            for key, val in case.items():
                if callable(val):
                    case[key] = val(case)
            self.cases.append(case)

    def run_case(self, runner: str, case: dict) -> typing.Tuple[int, str]:
        """
        Run a case in a fixed case-specific location inside self.run_dir.
        If data already exists in this folder (i.e. from a previous run), the case directory is first rm -rf'd.
        :param runner: mpirun-like parallel run command (must support the '-n' argument),
                       may be prefixed with arguments (like -q for aprun)
        :param case: case configuration
        :return: the program return code and the directory the case was run in (as a tuple)
        """
        run_dir = os.path.join(self.run_dir, str(_case_hash(case)))

        # (re)create run_dir
        if os.path.exists(run_dir):
            shutil.rmtree(run_dir)
        os.makedirs(run_dir)

        for generator in self.generators:
            # generator[0] is the callable function, generator[1] is the condition for the generator
            gen = generator[0]
            cond = generator[1]
            if cond is None:
                gen(case, run_dir)
            elif all(item in case.items() for item in cond.items()):
                gen(case, run_dir)

        # save case data as a Python pickle file
        with open(os.path.join(run_dir, 'case.pkl'), 'wb') as f:
            pickle.dump(case, f)

        # run the case
        start_time = time.perf_counter()
        with open(os.path.join(run_dir, 'output.txt'), 'wb') as outf:
            cmd = '{} -n {} {}'.format(runner, case.get('ntasks', 1), self.run_cmd.format(**case))
            print("Running {}...".format(self._case_name(case)))
            ok = subprocess.call(cmd, shell=True, stdout=outf, stderr=subprocess.STDOUT, cwd=run_dir)
        end_time = time.perf_counter()
        print("  (completed in {:.2f}s)".format(end_time - start_time))
        # move the files in root directory
        for root, dirs, files in os.walk(run_dir):
            for file in files:
                try:
                    shutil.move(os.path.join(root, file), run_dir)
                except OSError:
                    pass
        return ok, run_dir

    def set_consistent_keys(self, keys: typing.List[str]) -> None:
        """
        Set a list of "consistent keys." Cases that vary only by a "consistent key" are checked to have the same values.
        For example, consider the following case list:
            [ {'refine_lvl': 4, 'ntasks': 1}, {'refine_lvl': 4, 'ntasks': 8}, {'refine_lvl': 5, 'ntasks': 1} ]
        If 'ntasks' is a consistent key, the first two cases will be grouped together
        (as they have the same refine_lvl, and only ntasks varies) and checked that their results match.
        This checking is performed after data generation in "_ensure_consistency()".
        :param keys: list of "consistent" keys
        """
        self.consistent_keys = keys

    def set_exclude_patterns(self, exclude_patterns: typing.List[str]) -> None:
        """
        List of patterns describing files to leave out of the reference data archives (i.e. VTK files).
        This is done to make reference data archives smaller.
        :param exclude_patterns: list of Unix shell-style patterns for what files to exclude from reference data archive
                                 (see fnmatch for supported syntax - wildcards like '*.vtk' work)
        """
        self.results_exclude_patterns = exclude_patterns

    def _ensure_consistency(self):
        """
        Ensure that test cases that vary only by self.consistent_keys have the same reference data.
        This is used to guarantee that parallel runs give the same results as serial runs.
        See "set_consistent_keys" for more information.
        :return:
        """
        # TODO
        pass

    def _case_name(self, case: dict) -> str:
        """
        Returns a human-readable name for a case. Hides keys that are the same for all test cases.
        Do not rely on this name staying the same.
        :param case: case configuration
        :return: currently, a name in the form of "varied_key = value, ..."
        """
        varied_keys = get_varied_keys(self.cases)  # TODO cache this
        items = filter(lambda x: x[0] in varied_keys or len(varied_keys) == 0, case.items())
        items = sorted(items)
        return ', '.join(['{} = {}'.format(*x) for x in items])

    def run_main(self):
        parser = argparse.ArgumentParser()
        parser.add_argument('--clean', action='store_true', default=False,
                            help="Remove reference data belonging to cases no longer in the database. May be combined with --generate.")
        parser.add_argument('--generate', action='store_true', default=False,
                            help="Regenerate and save reference data.")
        parser.add_argument('--filter', nargs='+', default=None,
                            help="Act on only the subset of tests whose configurations match this Python expression.")
        parser.add_argument('--full', action='store_true', default=False,
                            help="Print all results (including passed tests).")
        parser.add_argument('--keep-data', action='store_true', default=False, help="Don't clean up results.")
        parser.add_argument('--runner', type=str, default='mpirun', help="Specify an alternative to mpirun.")
        args = parser.parse_args()

        # remove reference data from cases that no longer exist (common if keys are added, which changes all hashes)
        if args.clean:
            self.results_archive.remove_all_except(self.cases)
            if not args.generate:
                sys.exit(0)

        # generate reference data
        generate_cases = _filter_cases(self.cases, args.filter) if args.generate else []
        if len(generate_cases) > 0:
            failed = False
            for case in generate_cases:
                rc, result_dir = self.run_case(args.runner, case)
                if rc != 0:
                    color = "\033[93m"
                    end_color = "\033[0m"
                    print("{}ERROR: Case '{}' gave non-zero return code.{} Will not archive results."
                          .format(color, self._case_name(case), end_color))

                    # print last 30 lines of output.txt, indented with 2 spaces
                    log_path = os.path.join(result_dir, 'output.txt')
                    if os.path.exists(log_path):
                        print("Last 30 lines of output:")
                        with open(log_path, 'r') as f:
                            print('  \n'.join(f.read().split('\n')[-30:]))

                    failed = True
                    continue

                self.results_archive.add(case, result_dir, self.results_exclude_patterns)

            self._ensure_consistency()

            # clean up data
            if not args.keep_data and not failed and os.path.exists(self.run_dir):
                shutil.rmtree(self.run_dir)

            return (1 if failed else 0)

        # perform tests
        test_cases = _filter_cases(self.cases, args.filter)
        results = [list() for _ in test_cases]
        if len(test_cases) > 0:
            for case_idx, case in enumerate(test_cases):
                if not self.results_archive.has_reference(case):
                    results[case_idx].append(TestFailed(None, 'Missing reference data (run --generate)'))
                    continue

                rc, result_dir = self.run_case(args.runner, case)
                if rc != 0:
                    # print last 30 lines of output.txt, indented with 2 spaces
                    log_path = os.path.join(result_dir, 'output.txt')
                    msg = 'Non-zero return code {}'.format(rc)
                    if os.path.exists(log_path):
                        msg += "\n  Last 30 lines of output:"
                        with open(log_path, 'r') as f:
                            msg += '    \n'.join(f.read().split('\n')[-30:])
                    results[case_idx].append(TestFailed(None, msg))

                with self.results_archive.extract_reference(case) as ref_path:
                    for test in self.metrics:
                        result = test.test(ref_path, result_dir)
                        results[case_idx].append(result)

            # print results
            print("\n-----------------\n")
            for case, this_results in zip(test_cases, results):
                n_passed = sum([1 if type(r) is TestPassed else 0 for r in this_results])
                this_failed = n_passed != len(this_results)
                color = "\033[91m" if this_failed else ""
                end_color = "\033[0m"
                print("{}{} ({}) - {} ({} of {} passed){}"
                      .format(color, self._case_name(case), _case_hash(case), "failed" if this_failed else "OK",
                              n_passed, len(this_results), end_color))

                if this_failed or args.full:
                    for test in this_results:
                        color = "" if type(test) is TestPassed else "\033[91m"
                        end_color = "\033[0m"
                        print("  {}{}{}".format(color, str(test), end_color))

            failed = any([type(r) is not TestPassed for rr in results for r in rr])

            # clean up data
            if not args.keep_data and not failed and os.path.exists(self.run_dir):
                shutil.rmtree(self.run_dir)

            return(1 if failed else 0)

        parser.print_help()
        return(0)
