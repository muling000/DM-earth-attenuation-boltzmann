This is an example of the attenuation of the Solar DM.

First, run the simulation

   $ mpiexec -n N darkprop solardm.toml

where N should be replaced with the number of CPU cores.
After waiting about 10 core hours, 10^7 samples will be collected for each cross section.

Then the DM flux can be reconstructed using the Python script ``analysis.py``

   $ python3 analysis.py

which requires Python packages, `numpy`, `scipy`, `h5py`, and `matplotlib`.
The figure ``solardm-flux.png`` will be generated.

Note that in order to save storage space, the grid points of the tables in the ``data``
directory used for interpolation are very sparse.
