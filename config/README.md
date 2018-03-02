# Scripts for running code
There are many command line arguments for running `h5group-sorter`.
You can check them out by running
```sh
./h5group-sorter -h
```
after compiling the code. We provide a few scripts to run the code.

## Scripts to load modules
* `module.sh`: for LANL clusters
* `module_cray.sh`: for Cray systems at NERSC

## Scripts to check HDF5
* `check_hdf5.sh`: check the information of each tracer variables.
You are going to the indices of `Ux` and `q` for sorting.
    - The position each variable may not be the same as `h5dump -H` gives.

## Tracer sorting scripts
These scripts will use `h5group-sorter` to sort particles based on different
variables, such as particles tags and energy.

## Scripts to get particle trajectories
* `trajectory.sh`: get particle trajectories from tag-sorted tracers. It requires
that the tracers at each time frames have sorted by theirs tags.
    - `ratio_emax`: starting maximum energy / maximum energy of all particles
