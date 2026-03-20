import reframe as rfm
import reframe.utility.sanity as sn


class X3D2BuildMixin:
    """Common build configuration for x3d2 tests."""
    valid_systems = ['*']
    valid_prog_environs = ['*']
    sourcesdir = '../../'
    build_system = 'CMake'
    build_config = {
        'config_opts': ['-DWITH_2DECOMPFFT=ON', '-DBUILD_TESTING=ON'],
        'build_type': 'Release'
    }


@rfm.simple_test
class X3D2UnitTests(X3D2BuildMixin, rfm.RegressionTest):
    """Runs ctest -L unit"""
    executable = 'ctest'
    executable_opts = ['-L', 'unit', '--output-on-failure', '--timeout', '120']

    @rfm.run_before('sanity')
    def set_sanity_patterns(self):
        self.sanity_patterns = sn.assert_found(
            r'100% tests passed', self.stdout
        )


@rfm.simple_test
class X3D2VerificationTests(X3D2BuildMixin, rfm.RegressionTest):
    """Runs ctest -L verification"""
    executable = 'ctest'
    executable_opts = [
        '-L', 'verification', '--output-on-failure', '--timeout', '300'
    ]

    @rfm.run_before('sanity')
    def set_sanity_patterns(self):
        self.sanity_patterns = sn.assert_found(
            r'100% tests passed', self.stdout
        )


@rfm.simple_test
class X3D2PerfThom(X3D2BuildMixin, rfm.RegressionTest):
    """Thomas algorithm bandwidth benchmark"""
    executable = 'bin/perf_omp_thom_omp_1'

    @rfm.run_before('sanity')
    def set_sanity_patterns(self):
        self.sanity_patterns = sn.assert_found(
            r'PERF_METRIC: omp_thom_periodic', self.stdout
        )

    @rfm.run_before('performance')
    def set_perf_patterns(self):
        self.perf_patterns = {
            'periodic_bw': sn.extractsingle(
                r'PERF_METRIC: omp_thom_periodic time=\S+s bw=(\S+)',
                self.stdout, 1, float
            ),
            'dirichlet_bw': sn.extractsingle(
                r'PERF_METRIC: omp_thom_dirichlet time=\S+s bw=(\S+)',
                self.stdout, 1, float
            ),
        }

        self.reference = {
            '*': {
                'periodic_bw': (1.0, -0.2, None, 'GiB/s'),
                'dirichlet_bw': (1.0, -0.2, None, 'GiB/s'),
            }
        }


@rfm.simple_test
class X3D2PerfTridiag(X3D2BuildMixin, rfm.RegressionTest):
    """Distributed tridiagonal solver benchmark"""
    num_tasks = 4
    executable = 'bin/perf_omp_tridiag_omp_4'

    @rfm.run_before('sanity')
    def set_sanity_patterns(self):
        self.sanity_patterns = sn.assert_found(
            r'PERF_METRIC: omp_tridiag_periodic', self.stdout
        )

    @rfm.run_before('performance')
    def set_perf_patterns(self):
        self.perf_patterns = {
            'bw_min': sn.extractsingle(
                r'PERF_METRIC: omp_tridiag_bw_min bw=(\S+)',
                self.stdout, 1, float
            ),
            'bw_max': sn.extractsingle(
                r'PERF_METRIC: omp_tridiag_bw_max bw=(\S+)',
                self.stdout, 1, float
            ),
        }

        self.reference = {
            '*': {
                'bw_min': (1.0, -0.2, None, 'GiB/s'),
                'bw_max': (1.0, -0.2, None, 'GiB/s'),
            }
        }
