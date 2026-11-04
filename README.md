# NestedLogit

https://github.com/tothbog/NestedLogit

Interval branch-and-bound solver for the continuous competitive single facility
location and design problem with the nested logit choice rule, together with the
test instances used in the article

> J. Fernández, B. G.-Tóth, G. Méndez-Vogel, V. Marianov:
> *A continuous competitive single facility location and design problem with the
> nested logit choice rule.* Submitted to Omega.

The solver finds the location (x1, x2) and the quality (alpha) of a new facility
of the locating chain that maximise its profit, with guaranteed bounds on the
global optimum, using interval arithmetic and automatic differentiation.

## Requirements

- A C++ compiler (tested with g++ on Linux, x86-64).
- [PROFIL/BIAS](https://www.tuhh.de/ti3/keil/profil/) interval library, version
  2.0.8. Build it following its own instructions; it provides the libraries
  `libProfil`, `libBias` and `liblr`.

## Building

```bash
make PROFIL_DIR=/path/to/Profil-2.0.8
```

This creates the executable `NestedLogit` in the top directory. If `PROFIL_DIR`
is not given, `$HOME/programs/Profil-2.0.8` is used.

## Running

```bash
./NestedLogit -f instances/e1-12 -ep 1e-6 -mt -fv
```

Main options:

| Option      | Meaning                                                        |
|-------------|----------------------------------------------------------------|
| `-f file`   | instance file (required)                                       |
| `-ep eps`   | tolerance of the stopping criterion (default 1e-3)             |
| `-mt`       | use the monotonicity test                                      |
| `-fv`       | fix fi(i,2) = 1.25 and v(i) = 1 for all demand points          |
| `-ef`       | bound the subboxes with the gradient of the father box         |
| `-ub`       | update the upper bound at the Baumann point                    |
| `-vi`       | use the probabilistic choice rule instead of the nested logit  |
| `-h`        | print the list of all options                                  |

The options `-ef` and `-ub` need `-mt`. The `-f` option must be the
first one, because the names of the output files are built from the remaining
options.

The solver writes a summary to standard error and creates two files next to the
instance: `<instance>-ep<eps>-<options>.res` (summary of the run) and
`<instance>-ep<eps>-<options>_res.tclin` (the final list of boxes).

To solve all instances listed in a file:

```bash
./run.sh p-e1 1e-6 -fv
```

The article's settings were `-ep 1e-6 -mt -fv` for the Murcia instances and
`-ep 1e-4 -mt -fv` for the random instances (`./run.sh rand_probs 1e-4 -fv`).

## Instances

`instances/` holds the instances based on the Murcia data (71 demand points).
Instance names follow `<case>-<l1><lc>[-a|-x]`:

- `e0`, `e1`, `e2`: 0, 1 or 2 existing facilities of the locating chain.
- `<l1><lc>` gives (lambda_1, lambda_c): `1` = 1.0, `2` = 2.0, `5` = 1.5
  (for example `12` means lambda_1 = 1, lambda_c = 2).
- `-a`: quality fixed (alpha = 3); `-x`: location fixed to the optimal one,
  only the quality is optimised.
- `b` (e.g. `e1b`): variant with reduced demand for Murcia (from 10 to 4).

The files `p-*` list the instances of each experiment, `rand_probs` lists the
450 random instances in `instances/random/`. Random instance names are
`r-<n>-<m>-<k>-<id>-<l1><lc>`, where n is the number of demand points, m the
number of existing facilities, k the number of those belonging to the locating
chain, and `<l1><lc>` is as above.

### Instance file format

Plain text, with a label line before each block of values: problem name,
number of demand points (n), of existing facilities (m) and of facilities owned
by the locating chain (k); bounds of x1, x2 and alpha; parameters of the
distance function, of the G(alpha) cost function, of the nested logit model
(lambda_1, lambda_c) and of F; the Psi parameters, locations of the existing
facilities and of the demand points, buying power, and the attractiveness
values of the existing facilities. The last two fields (`Minima`, `Minimo`) are
not used by the solver, except that `-ft` evaluates the objective at `Minima`.

## Random instance generator

`generator/` contains the program that created the random instances:

```bash
cd generator
make PROFIL_DIR=/path/to/Profil-2.0.8
./run_generator.sh params.txt      # creates r-<n>-<m>-<k>-<id>-11 files
./lambda.sh r-100-6-3-0-11         # creates the -12, -21, -22, -51, -52 variants
```

The generator is seeded with the current time, so a new run produces new
instances; the instances used in the article are the ones in
`instances/random/`.

## Third-party code

`src/HESS_ari.*`, `src/XI_ari.*`, `src/I_util.*` and `src/Global.h` are based
on the C++ Toolbox for Verified Computing I (R. Hammer, M. Hocks, U. Kulisch,
D. Ratz, Springer, 1995), adapted to PROFIL/BIAS by M. Cs. Markót. They keep
their original copyright notice and terms (free for non-commercial use).

## License

Except for the third-party files listed above, the code is distributed under the
GNU General Public License v3.0 or later, see `LICENSE`.

> J. Fernández, B. G.-Tóth, G. Méndez-Vogel, V. Marianov:
> *A continuous competitive single facility location and design problem with the
> nested logit choice rule.* Submitted to Omega.
