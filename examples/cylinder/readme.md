# Generate the IBM file

PYTHONPATH=~/opt/adios2/build/local/lib/python3.12/dist-packages/ python3 ~/opt/py4x3d2/tests/run.py -v --iibm 1 --dx 1 --dy 1 --dz 1 --nx 257 --ny 128 --nz 32 --save --cyl 6.4 64 64 16 0 0 1

This produces `ibm.bp` in the current directory. Before running the case, it must be moved to the `x3d2` directory and renamed according to the Poisson solver used.
 
The filename convention is `ibm_XYZ.bp` where `X`, `Y`, `Z` reflect the parity of `nx`, `ny`, `nz` respectively — `0` if even, `1` if odd.
 
## Example scenarios

| nx  | ny  | nz | Target filename |
|-----|-----|----|-----------------|
| 256 | 128 | 32 | `ibm_000.bp`    |
| 257 | 128 | 32 | `ibm_100.bp`    |
| 256 | 129 | 32 | `ibm_010.bp`    |
| 256 | 128 | 33 | `ibm_001.bp`    |
| 257 | 129 | 33 | `ibm_111.bp`    |
 

# Run the case

LD_LIBRARY_PATH=~/opt/adios2/build/lib/:$LD_LIBRARY_PATH mpirun -n 1 ../build/bin/xcompact input.x3d
