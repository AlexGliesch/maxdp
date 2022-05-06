# Supplementary material for *A hybrid heuristic for the Maximum Dispersion Problem*

This repository holds the source code, instance generator and detailed result tables for the following paper:

> **A hybrid heuristic for the Maximum Dispersion Problem**<br>
> Alex Gliesch, Marcus Ritt (2021) <br>
> European Journal of Operational Research 288.3, pp. 721-735 <br>
> https://doi.org/10.1016/j.ejor.2020.06.011

**Bibtex**

```bibtex
@article{GlieschRitt/2021,
  title     = {A hybrid heuristic for the Maximum Dispersion Problem},
  author    = {Gliesch, Alex and Ritt, Marcus},
  journal   = {European Journal of Operational Research},
  volume    = {288},
  number    = {3},
  pages     = {721--735},
  year      = {2021},
  publisher = {Elsevier},
  doi       = {10.1016/j.ejor.2020.06.011}
}
```

Please use the reference above if you use this material in your research.

## Dependencies

- CPLEX 20.1.0
- [boost 1.60.0](https://www.boost.org/users/history/version_1_60_0.html)
- [fmtlib](https://github.com/fmtlib/fmt)
- [GNU parallel](https://www.gnu.org/software/parallel/)

## Generating instances
1. Under `src`, run `cmake .` then `make`.
1. Generate the test instances used in the paper with `./generateinstances.sh`. The instances will be placed in `src/inst`.
1. To generate more instances as you like, see `./generateinstance --help`.

## Running the code 
1. Under `src`, run `cmake .` then `make`.
1. Run using `./maxdp -i {instance} -t {timeLimit} -o {outFile}`. For more options, see `--help`.
1. For our implementation of Fernandez et al. (2013)'s solver, under `src/julia` run `./fkn.sh {instance} {timeLimit} {alpha}`.

