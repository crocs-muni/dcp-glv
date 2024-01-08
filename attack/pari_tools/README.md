## Implementation of DCP solver and baby-step giant-step in pari/gp

Run `make` to compile.

Run `make test` and `./test` to run tests.

Description of files:

`dcp_pari.py`
- Wrapper for the dcp solvers (below) in Python. Pari binaries are run using os.system.

`dcp_solver.c`
- implementation of dcp solver simple scalar multiplication.

`dcp_glv_solver.c`
- implementation of dcp solver for 2-dimensional scalar decomposition.

`bsgs.cpp`
- implementation of aby-step Giant-step for 2-dim scalar decomposition. Should be run from `bsgssolver.py`.


Note: pari stack can be changed in the `c/cpp` files.