# bskit

`bskit` is a Python package for measuring density bispectra from snapshots of cosmological N-body or hydrodynamical simulations. It can measure auto or cross bispectra in a user-specified set of triangle bins (that is, triplets of 3-vector wavenumbers (k_1,k_2,k_3) such that k_1, k_2, and k_3 fall into separately specified vector magnitude bins and k_1+k_2+k_3=0). Several common sets of bins are also implemented:
- all triangle bins with min(|k_1|,|k_2|,|k_3|) > k_min and max(|k_1|,|k_2|,|k_3|) < k_max for specified k_min and k_max;
- equilateral triangles between specified k_min and k_max;
- isosceles triangles defined by m*|k_1| ~ |k_2| ~ |k_3| for specified multiplier m;
- squeezed isosceles triangles with |k_1| fixed and |k_2| ~ |k_3|.

`bskit` is built upon the [nbodykit](github.com/bccp/nbodykit) simulation analysis package, and users should familiarize themselves with the central `nbodykit` concepts in the [documentation](https://nbodykit.readthedocs.io/en/latest/) before getting started with `bskit`.

This package uses the FFT-based bispectrum measurement algorithm presented e.g. in Sec. 2 of

> Tomlinson et al., *Efficient parallel algorithm for estimating higher-order polyspectra*, [Astrophys. J., 158, 116 (2019)](10.3847/1538-3881/ab3223), arXiv:[1904.11055](https://arxiv.org/abs/1904.11055)

(see below for further references), and was first used for the bispectrum measurements in the following paper:

> Foreman et al., *Baryonic effects on the matter bispectrum*, 2019, arXiv:19mm.xxxxx

## Installation

`bskit`’s main dependency is on `nbodykit`. After following the `nbodykit` [installation instructions](https://nbodykit.readthedocs.io/en/latest/getting-started/install.html) (preferably via `conda`), all `bskit` dependencies will also be installed. At this point, simply clone this repository and ensure that the root `bskit` directory is in your `PYTHONPATH`, e.g. via 

``
export PYTHONPATH=/path/to/bskit:$PYTHONPATH
``

## Usage

Usage instructions and a guide to the included examples can be found [here]()

## References

If `bskit` is used in original research, please cite the associated paper:

> Foreman et al., *Baryonic effects on the matter bispectrum*, 2019, arXiv:19mm.xxxxx

In addition, please cite the `nbodykit` paper,

> Hand et al., *nbodykit: an open-source, massively parallel toolkit for large-scale structure*, [Astron. J., 156, 160 (2018)](https://dx.doi.org/10.3847/1538-3881/aadae0), arXiv:[1712.05834](https://arxiv.org/abs/1712.05834)

and the following standard references for the FFT-based bispectrum estimator that the package implements:

> Scoccimarro, *The Bispectrum: From Theory to Observations*, [Astrophys. J., 544, 597 (2000)](https://dx.doi.org/10.1086/317248), arXiv:[astro-ph/0004086](https://arxiv.org/abs/astro-ph/0004086)

> Sefusatti et al., *Accurate estimators of correlation functions in Fourier space*, [Mon. Not. Roy. Astron. Soc., 460, 3624 (2016)](https://dx.doi.org/10.1093/mnras/stw1229), arXiv:[1512.07295](https://arxiv.org/abs/1512.07295) 

