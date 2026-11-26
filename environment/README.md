# Computing environment

The saved results were checked with R 4.3.1 and Python 3.12. `R-dependencies.csv` records the R package versions and their dependencies; `R-session.txt` records the R environment. Python versions are pinned in the root `requirements.txt`.

Check an installed R library:

```sh
Rscript --vanilla scripts/check_R_environment.R
```

For a separate R library, download and install the recorded sources explicitly:

```sh
Rscript --vanilla scripts/check_R_environment.R --download --source-dir environment/sources
Rscript --vanilla scripts/check_R_environment.R --install --source-dir environment/sources --library environment/R-library
R_LIBS_USER="$PWD/environment/R-library" Rscript --vanilla scripts/check_R_environment.R
```

Keep `R_LIBS_USER` set for subsequent R commands if you use that library. Source installation needs the usual C/C++/Fortran tools and libraries for GMP/MPFR. The full R source installation has not been tested on a new machine. A different BLAS, solver, or package version can change numerical results; each calculation saves its session information.

The main R dependencies include MASS, Matrix, glmnet, randomForest, ggplot2, scales, digest, jsonlite, hdi, SIHR, lpSolve, gmp, Rmpfr, CVXR, ECOSolveR, rags2ridges, bestNormalize, VGAM, igraph, and glasso. The connected-target and independent-target power examples primarily use base R and jsonlite. SIHR 2.0.1 is the recorded version for the CGM comparator.
