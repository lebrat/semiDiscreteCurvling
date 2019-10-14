# semi-discrete optimal transport for Curvling
______________________
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3484184.svg)](https://doi.org/10.5281/zenodo.3484184)

This codes are the implementation of the following [paper](https://arxiv.org/abs/1804.08356). It allows to compute the exact *L²*  optimal transport between a continuous measure (with density) and measure carried by a set of Dirac masses. This toolbox allows the calculation in **2D** of the optimal Transport distance. Moreover it includes Curvling constraints, that is imposing to the Diracs masses to be taken along a curve with a bounded speed and acceleration. We provide hands-on tutorials on the [wiki](https://github.com/lebrat/semiDiscreteCurvling/wiki).

The codes are released for Linux platforms (tested for *Mint 18 Cinnamon 64-bit* and *Ubuntu 19.04 Disco Dingo*). 

The back-end computations are coded in `C++` and make use of the computational geometry library [CGAL](www.cgal.org) and the linear algebra library [Eigen3](http://eigen.tuxfamily.org). We provide a `python 3.7` interface by using the wrapper [swig](http://www.swig.org/).



## Authors
This software was developed by:
* Frédéric de Gournay
* Jonas Kahn
* Alban Gossard
* Pierre Weiss
* [Léo Lebrat](lebrat.org)

All members of the [Toulouse institute of Mathematics](https://www.math.univ-toulouse.fr/?lang=en) France, granted by the **ANR-17-CE23-0013**.

## License

This software is open source distributed under license MIT.
