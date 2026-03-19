# ReFrame site-specific settings template for x3d2
#
# Copy this file and customise for your HPC system.
# See https://reframe-hpc.readthedocs.io/en/stable/configure.html

site_configuration = {
    'systems': [
        {
            'name': 'local',
            'descr': 'Local workstation',
            'hostnames': ['localhost'],
            'partitions': [
                {
                    'name': 'default',
                    'scheduler': 'local',
                    'launcher': 'mpirun',
                    'environs': ['builtin'],
                }
            ]
        },
        # Add your HPC system here, e.g.:
        # {
        #     'name': 'archer2',
        #     'descr': 'ARCHER2 UK National Supercomputer',
        #     'hostnames': ['ln\d+\.archer2\.ac\.uk'],
        #     'partitions': [
        #         {
        #             'name': 'compute',
        #             'scheduler': 'slurm',
        #             'launcher': 'srun',
        #             'access': ['--partition=standard', '--qos=standard'],
        #             'environs': ['gnu', 'cray'],
        #         }
        #     ]
        # },
    ],
    'environments': [
        {
            'name': 'builtin',
            'cc': 'gcc',
            'ftn': 'mpif90',
        },
    ],
}
