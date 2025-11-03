# Computing critical Rayleigh number using Kore

This benchmark uses a slightly modified version of the [`jones2000`](https://github.com/repepo/kore/tree/main/tests/jones2000) test of [kore](https://github.com/repepo/kore).
The `find_Rac.py` script has been adapted to include azimuthal wave numbers `m=1,2,3`, and a different `sed` command to work on `MacOS 15`. It also only saves the values of `Ra_c`, `m` and `ω_c` to the `critical_params.dat` file.

The `parameters.py` file has been changed w.r.t the benchmark in `kore` to use constant temperature BC at the top of the core.


## To run this benchmark
1) Install a working version of kore (this was computed with commit [`2306d1d7e2c397e32c1ed9c07b70c989a81e91c0`](https://github.com/repepo/kore/tree/2306d1d7e2c397e32c1ed9c07b70c989a81e91c0)).
2) Copy `find_Rac.py` into kore's root folder and copy the `parameters.py` file into `bin/parameters.py`.
3) Run the benchmark using `./find_Rac.py`.

This will generate the file `critical_params.dat` containing
```
8.25010970e+06 1 -5.70072886e-03
6.61663621e+06 2 -9.24957284e-03
5.79607663e+06 3 -1.18798295e-02
```

