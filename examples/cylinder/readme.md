# Generate the IBM file

PYTHONPATH=~/opt/adios2/build/local/lib/python3.12/dist-packages/ python3 ~/opt/py4x3d2/tests/run.py -v --iibm 1 --dx 1 --dy 1 --dz 1 --nx 256 --ny 128 --nz 32 --save --cyl 6.4 64 64 16 0 0 1

# Run the case

LD_LIBRARY_PATH=~/opt/adios2/build/lib/:$LD_LIBRARY_PATH mpirun -n 1 ../build/bin/xcompact input.x3d
