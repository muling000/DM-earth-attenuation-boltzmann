This is an example of the Earth attenuation of isotropic CRDM flux.

First, run the simulation

   $ mpiexec -n N darkprop config_crdm.toml

where N should be replaced with the number of CPU cores.
After waiting about 20 core hours, 10^7 samples will be collected for each
cross section.

Then the underground flux can be reconstructed using the Python script
`analysis.py`

   $ python3 analysis.py

which requires Python packages, `numpy`, `scipy`, `h5py`, and `matplotlib`.
The figure ``crdm-flux.png`` will be generated, which corresponds to Fig. 5 in
the paper arXiv: 2111:05559.
