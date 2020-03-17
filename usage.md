# Usage

`bskit` is documented with Python docstrings, and most included scripts have headers and command-line help. Here is a general guide on how to get started.

## Included files

This repository is organized as follows:

`bskit/`: source for module

`scripts/`: a variety of Python scripts that should cover the most common use cases, mostly accepting command-line inputs. They are sorted into scripts for processing simulation snapshots and density grids (`grids/`), measuring power spectra and bispectra (`measure/`), and post-processing the measurements (`process/`).

`examples/batch/`: examples of submitting the Python scripts to a cluster (specifically, a cluster using the `TORQUE` batch system), used for the measurements in Foreman et al. 2019.

`examples/tests/scripts/`: ready-to-run tests of various `bskit` functions, executable from the command line, which take the density grids in `examples/tests/grids/` as inputs and write outputs to `examples/tests/output_user/`.These outputs should be almost identical to the files in `examples/tests/output_ref/`.

## Getting started

For reasons of computational efficiency, measuring the bispectrum from a given snapshot is a multi-step process, split into computing binning information, computing unnormalized binned bispectrum values, and combining the two into a final file with all relevant information. The example scripts in `examples/tests/scripts/` will get you acquainted with the syntax for each step. Before running these examples, use `set_locations.sh` to set environment variables with various locations:

``
source set_locations.sh
``

And also run `examples/tests/grids/generate_gaussian_test_grids.py` to generate the test grids used for these examples.

1. Computing binning information: based on a specified simulation box size, Fourier grid resolution, and bispectrum binning scheme, compute the number of triangles within each bin, along with the mean k vector magnitudes within each bin, and save to disk. In cases where you will be measuring the bispectrum from multiple snapshots with the same binning scheme, this step only needs to be completed once.

> Example: `examples/tests/scripts/test_compute_bs_gridinfo.sh`

2. Computing unnormalized binned bispectrum values: for the binning scheme used in step 1, measure the bispectrum in each bin, and save to disk. The output bispectrum values will *not* be normalized by the number of triangles in each bin.

> Example: `examples/tests/scripts/test_measure_bs_faster.sh`

3. Combine the output files from steps 1 and 2: this will assemble a single file containing the following information for each bin: k magnitude means, k magnitude boundaries, correctly normalized bispectrum values, and number of triangles.

> Example: `examples/tests/scripts/test_process_bs_output.sh`

### Other examples in `examples/tests/scripts/`:

`test_compute_bs_gridinfo_slow.sh`: Compute the binning information using an algorithm that contains a lot of redundant computation, and is therefore slower than `test_compute_bs_gridinfo.sh`, but is much less memory-intensive.

`test_compute_bs_gridinfo_subbox.sh`: Compute the binning information with a binning scheme specified over the entire simulation box, but only corresponding to a sub-volume of the full volume.

`test_measure_bs_faster_subboxes.sh`: Measure the unnormalized bispectrum from sub-volumes of the simulation box, using the same binning scheme specified for the full box.

`test_measure_bs_slower.sh`: Compute binning information *and* the normalized bispectrum, using the slower, less memory-intensive algorithm.

`test_measure_cross_bs_eq_slower.sh`: Compute the cross bispectrum (and binning information) between two density grids on equilateral triangles only, using the slower algorithm

The final 3 examples don’t use anything from `bskit` specifically, instead using only standard `nbodykit` routines:

`test_measure_cross_ps_from_bigfile.sh`: Measure the cross power spectrum between two density grids, along with the mean k and number of modes per bin.

`test_measure_ps_from_bigfile.sh`: Measure the (auto) power spectrum of a given snapshot, along with the mean k and number of modes per bin.

`test_measure_ps_from_bigfile_subboxes.sh`: Measure the (auto) power spectrum of sub-volumes of a given snapshot, along with the mean k and number of modes per bin, using the k bins defined for the full simulation box.


## Example batch scripts

Many of the scripts in `examples/tests/scripts/` largely overlap with those in `examples/tests/scripts/`, but also include examples of converting particle snapshots into density grids that are saved to disk, or downsampling density grids to a lower resolution. Converting particles to density as a separate step, instead of directly loading the particle snapshots into `bskit` routines, can save disk space (since the density grids typically take up less) or computation time (if you expect to make multiple measurements on the same snapshot, it is more efficient to load the density grid than the raw particles each time).


## Included general scripts

The above examples mostly call one of the scripts located in the main `scripts/` directory. Each script also contains command-line help, if your desired use is not covered by one of the examples.

The names of the scripts in `scripts/grids/` and `scripts/measure/` are self-explanatory, but those in `scripts/process/` are described below.

`transform_bininfo_file.py`: Takes an existing file with binning info for a given binning scheme (i.e. bin boundaries), box size, and grid resolution, and transforms it to correspond to a different box size, with the same grid resolution and binning scheme in terms of the fundamental mode of the new box size.

`process_fast_bs_measurement.py`: Takes files with binning info and unnormalized bispectrum values (along with a specified k_max), and combines them into a single file containing bin info and properly normalized bispectrum values.

`process_ps_subboxes_into_errorbars.py`: Takes files (generated by `scripts/measure/measure_ps_from_bigfile_subboxes.py`) with power spectrum measurements from sub-volumes of two different snapshots, and computes an estimate of the sample variance uncertainty on the *ratio* of power spectra from the snapshots, based on the spread in the ratios from each sub-volume. See Foreman et al. 2019 for more details.

`process_bs_subboxes_into_errorbars.py`: As above, but for bispectrum measurements from sub-volumes (generated by `scripts/measure/measure_bs_from_bigfile_subboxes.py`).

`make_all_…`: Example shell scripts for processing batches of power spectrum or bispectrum measurements to produce sample variance errorbars or combined bispectrum files. (These scripts have everything hard-coded within them, and will need to be modified directly before use).


## `bskit` details

Most of the action happens within the `FFTBispectrum` class, based on the `FFTPower` class in `nbodykit`. When initialized, an instance of `FFTBispectrum` takes a `CatalogSource` or `MeshSource` object corresponding to the snapshot of interest (or more than one object, for cross-correlations), along with details about bispectrum binning. The `measure_bispectrum` routine measures the bispectrum (along with binning info) over a range of bins, while the `measure_bispectrum_faster` and `measure_gridinfo_faster` routines compute unnormalized bispectrum and binning info values separately, using a faster but more memory-intensive algorithm than `measure_bispectrum`. Calling these routines with an output file name causes computations to be written to disk intermittently (which is recommended, to guard against crashes or timeouts), but the `save_bispectrum` routine can also save the computations to disk separately. Further details, as well as other routines that can be directly applied to a density grid without initializing an `FFTBispectrum`, can be found be exploring the `bskit/main.py` file directly.
