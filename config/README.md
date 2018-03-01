# Scripts for running code

## Module loading scripts
* `module.sh`: for LANL clusters
* `module_cray.sh`: for Cray systems at NERSC

## Tracer sorting scripts
These scripts will use `h5group-sorter` to sort particles based on different
variables, such as particles tags and energy. In the root directory,
```sh
./h5group-sorter -h
```
to check out the commandline arguments to run the program.

sorter.sh: sort particles with different key
* change filepath to the directory containing the tracked trajectories
* particle: ion or electron
* the rest of the script loops over all time frames
  -- The configuration information is explained in configuration.c
  -- -k changes the key:
    0: Ux
    1: Uy
    2: Uz
    3: dX
    4: dY
    5: dZ
    6: i
    7: q

energy_sorter.sh: sort particles with energy as the key
* change tstep to the last time step
* the rest is similar as sorter.sh
* -k argument must be 8, so it will sort the particles with energy as the key

trajectory.sh: get high energy particle trajectories
* -n 1000 indicates number of particles
* -r argument can be modified. 10 indicates 1/10 of the maximum energy, while
  100 indicates 1/100 of the maximum energy

Use one one script to get the particle trajectory
time_steps.sh:
* sort the particle in the step by energy
* get particle tags for the particles you are interested in
  --nptl_traj: number of particles you need
  --ratio_emax: ration of th maximum energy of all particles over
                that of tracked particles
  -w: won't save the sorted particle in the middle steps. Removing it will save
      all the sorted particles.
