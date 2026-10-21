# Generate the IBM file

The mask is generated with [x3d2-tools](https://github.com/xcompact3d/x3d2-tools):

PYTHONPATH=~/opt/adios2/build/local/lib/python3.12/dist-packages/ python3 ~/opt/x3d2-tools/tests/run.py -v --iibm 1 --dx 1 --dy 1 --dz 1 --nx 257 --ny 128 --nz 32 --save --cyl 6.4 64 64 16 0 0 1

This produces `ibm.bp` in the current directory. Move it next to `input.x3d`
and point `ibm_file` in `&solver_params` at it (the default is `ibm.bp`):

```
&solver_params
...
ibm_on = T
ibm_file = 'ibm.bp'
/End
```

The mask holds `ep1` on the mesh vertices, so it must have the same size as
`dims_global` (`nx`, `ny`, `nz` above). x3d2 checks this when it reads the
file and stops if the sizes differ; regenerate the mask after changing the
grid.

If you write the mask yourself from Python, do as x3d2-tools does: store a
C-ordered array of shape `(nz, ny, nx)` and declare the ADIOS2 variable with
`shape=[nz, ny, nx]`.

# Run the case

LD_LIBRARY_PATH=~/opt/adios2/build/lib/:$LD_LIBRARY_PATH mpirun -n 1 ../build/bin/xcompact input.x3d
