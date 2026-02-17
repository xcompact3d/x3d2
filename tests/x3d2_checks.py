import reframe as rfm
import reframe.utility.sanity as sn

@rfm.simple_test
class X3D2UnitTest(rfm.RegressionTest):
    """
    Runs the unit tests via CTest.
    Useful for checking if the build works on the HPC environment.
    """
    valid_systems = ['*']
    valid_prog_environs = ['*']
    
    # ReFrame will copy the source to a stage directory
    sourcesdir = '../'
    
    build_system = 'CMake'
    build_config = {
        'config_opts': ['-DWITH_2DECOMPFFT=ON', '-DBUILD_TESTING=ON'],
        'build_type': 'Release'
    }

    # We use ctest to drive the execution of unit tests
    executable = 'ctest'
    executable_opts = ['-L', 'unit', '--output-on-failure']

    @run_before('sanity')
    def set_sanity_patterns(self):
        # Check that 100% of tests passed
        self.sanity_patterns = sn.assert_found(r'100% tests passed', self.stdout)

@rfm.simple_test
class X3D2PerformanceBenchmark(rfm.RegressionTest):
    """
    Runs a specific performance test and tracks the timing.
    """
    valid_systems = ['*']
    valid_prog_environs = ['*']
    sourcesdir = '../'
    
    build_system = 'CMake'
    build_config = {
        'config_opts': ['-DWITH_2DECOMPFFT=ON', '-DBUILD_TESTING=ON'],
        'build_type': 'Release'
    }

    # ReFrame handles the MPI launch (srun/mpirun), so we point to the binary directly.
    # Note: The binary name depends on your CMake logic: test_name + _ + backend + _ + np
    executable = 'tests/bin/test_mesh_omp_4'
    num_tasks = 4

    @run_before('sanity')
    def set_sanity_patterns(self):
        self.sanity_patterns = sn.assert_found(r'Test Passed', self.stdout)

    @run_before('performance')
    def set_perf_patterns(self):
        # Example: Extract a value like "Time: 0.123s" from stdout
        # You will need to ensure your Fortran test prints this format.
        self.perf_patterns = {
            'execution_time': sn.extractsingle(r'Time:\s+(\S+)s', self.stdout, 1, float)
        }
        
        self.reference = {
            '*': {'execution_time': (0.5, -0.1, 0.1, 's')} # ref, -10%, +10% tol
        }