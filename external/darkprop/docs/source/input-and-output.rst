Input and Output
================

A complete example of the configuration file is shown below, which can be found in the
`app` directory of the source code.

.. literalinclude:: ../../app/config.toml
   :caption: complete DarkProp configuration file
   :language: toml

This is a TOML format configuration file consisting of parameter-value pairs.
Some parameters can be missing and the default values will be used. The parameters are
grouped into three main sections (tables) ``input``, ``output``, and ``numerical`` which
are written in square brackets on a line by themselves. A section may have subsections
defined as the ``[section.subsection]`` form. Comments are started with a ``#`` symbol.

Using the configuration file above, DarkProp will simulate halo DM with mass
:math:`m_\chi = 0.1~\text{GeV}` and cross section
:math:`\sigma_{\chi p} = 10^{-30}~\text{cm}^2` in velocity range from 1 cm/s to 784 km/s
(in the Earth frame). The events of particles crossing a 1.4 km depth underground surface
will be stored in the file :file:`./out/example/example_m0_s0.hdf5`.

The Galactic coordinate system is used in :program:`darkprop` with the Galactic longitude
:math:`l` from :math:`0` to :math:`2\pi` and the Galactic latitude :math:`b` from
:math:`-\pi/2` to :math:`\pi/2`. The transformation between Cartesian coordinates and
Galactic coordinates is

.. math::

   x &= r \cos b \cos l, \\
   y &= r \cos b \sin l, \\
   z &= r \sin b.

The following is a detailed explanation of the input parameters and the output files.

Configuration parameters
------------------------

In this section we use the ``section.subsection.parameter`` notation for parameters. The
type of each parameter is marked after it and the unit of the dimensional parameters are
also indicated in square brackets behind. The default value, if defined, is written in
parentheses.

- id ``String``
   The ``id`` parameter should be defined before any other sections.
   Its value will be the prefix of the output files.

[input]
^^^^^^^

- input.type ``String``
   The type of the initial conditions. Options:

   :"halo": Halo DM simulation.
   :"flux": Isotropic arbitrary incident flux in tabular form.
   :"intensity": Anisotropic arbitrary incident flux in HDF5 format.
   :"sample": Samples of initial conditions in HDF5 format.
   :"source": Production of DM from internal of the medium.

- input.dm_model ``String``
   Choosing the DM model implemented in DarkProp. Options (case insensitive):

   :"SI": Spin independent constant cross section (class :ref:`SIDM`).
   :"SIHelm": Spin independent constant cross section with the Helm nuclear form factor
              (class :ref:`SIHelmDM`).
   :"DMElectron": DM-electron vector coupling (class :ref:`DMElectron`).
   :"SolarDM": The solar dark matter model (class :ref:`SolarDM`).

- input.medium_model ``String``
   Choosing the medium model implemented in DarkProp. Options (case insensitive):

   :"HomoEarth": A homogeneous Earth model composed of 8 components
                 (class :ref:`HomoEarth`).
   :"PREM": The Preliminary Reference Earth Model (class :ref:`PREMEarth`).
   :"HomoElectronEarth": A homogeneous Earth model composed of electrons
                         (class :ref:`HomoElectronEarth`).
   :"Sun": A Sun model (class :ref:`Sun`).

- input.mchis ``Array`` [GeV]
   List of DM particle masses.

- input.sigmas ``Array`` [:math:`\text{cm}^2`]
   List of DM-nucleon scattering cross sections :math:`\sigma_{\chi p}`.

- input.mchi_file ``String``
   Path of a text file containing a list of DM particle masses in GeV.
   Its effect is the same as the ``input.mchis`` parameter. If both parameters exist, it
   will override ``input.mchis``.

   The file should contain only one column of numbers, or a single line of numbers
   separated by spaces.

- input.sigma_file ``String``
   Path of a text file containing a list of cross sections in :math:`\text{cm}^2`.
   Its effect is the same as the ``input.sigmas`` parameter. If both parameters exist, it
   will override ``input.sigmas``.

   This file should contain only one column of numbers, or a single line of numbers
   separated by spaces.

- input.mmed ``Float`` [GeV] (default 0.0)
   The mediator mass.

- input.mediator_limit ``String`` (default "")
   The limitation case of the mediator mass :math:`m_\phi`. Options:

   :"light": :math:`m_\phi \to 0`.
   :"heavy": :math:`m_\phi \to \infty`.
   :"": An empty string means no limit is taken.

- input.Tcut ``Float`` [GeV]
   The termination energy of the simulation.
   The simulation of the trajectory will stop when the kinetic energy of the particle is
   less than ``input.Tcut``.

   If ``input.type`` is ``"halo"``, this parameter will be overwritten by the
   ``input.halo.v_cut`` parameter as :math:`T_\text{cut}= \frac{1}{2}m_\chi v_\text{cut}^2`.

   If ``input.type`` is ``"flux"``, this parameter can be absent, and :math:`T_\text{cut}`
   will be the lowest energy of the input flux.

- [input.halo]
   The parameters in this subsection only take effect if ``input.type`` is ``"halo"``.

   * input.halo.v_esc ``Float`` [km/s]
      The escape speed of the standard halo model (SHM) velocity distribution.

   * input.halo.v_0 ``Float`` [km/s]
      The most probable speed of the SHM velocity distribution.

   * input.halo.v_earth ``Array`` [km/s]
      The 3-D velocity of the Earth. The direction will be used to
      determine the isodetection rings.

   * input.halo.v_min ``Float`` [km/s]
      The lower bound of sampling velocity in the Earth frame.

   * input.halo.v_max ``Float`` [km/s]
      The upper bound of sampling velocity in the Earth frame.
      It will be :math:`v_\oplus + v_\text{esc}` if a negative value is given.

   * input.halo.v_cut ``Float`` [km/s]
      The termination speed of the simulation. It overwrites the ``input.Tcut``
      parameter.

- [input.flux]
   The parameters in this subsection only take effect if ``input.type`` is ``"flux"``.

   * input.flux.point_source ``Boolean`` (default ``false``)
      Choosing if injecting from a point source. Options:

      :``true``: DM particles come from the same direction.
      :``false``: DM particles are injected isotropically.

   * input.flux.flux_file ``String`` [:math:`\text{GeV}`, :math:`\text{GeV}^{-1}\text{cm}^{-2}\text{s}^{-1}\text{sr}^{-1}`]
      The path to a file containing the incident fluxes.

      The file should be a space separated CSV file with multiple columns.
      The first column (column 0) is the kinetic energy of the DM particle.
      The following columns are the corresponding DM fluxes of different masses and cross
      sections which are set by the ``input.mchis`` (or ``input.mchi_file``) and
      ``input.sigmas`` (or ``input.sigma_file``) parameters.
      These masses and cross sections constitute 2-dimensional grid points. The flux of the
      :math:`i`-th mass and the :math:`j`-th cross section is at the
      :math:`(1 + j + i \times N_s)`-th column, where :math:`N_s` is the number of cross
      sections. The unit of energy is GeV and the unit of flux is
      :math:`\text{GeV}^{-1}\text{cm}^{-2}\text{s}^{-1}\text{sr}^{-1}`.

   * input.flux.Tmin ``Float`` (default 0.0)
      The lower boundary of the sampling kinetic energy.
      If :math:`\leq 0`, it will be set to be the smallest energy in the `flux_file`.

   * input.flux.Tmax ``Float`` (default 0.0)
      The upper boundary of the sampling kinetic energy.
      If :math:`\leq 0`, it will be set to be the largest energy in the `flux_file`.

- [input.intensity]
   The parameters in this subsection only take effect if ``input.type`` is ``"intensity"``.

   * input.intensity.intensity_file ``String``
      The path to a file containing anisotropic incident fluxes.

      The file must be a HDF5 file with the following definitive structure:

      .. code-block:: text

         file
         |--mchi          {Nm}         DM particle mass [GeV]
         |--sigma         {SCALAR}     DM-nucleon cross section [cm^2]
         |--m_0           Group        the first DM particle mass
         |  |--mchi       {SCALAR}     DM particle mass in this group [GeV]
         |  |--b          {Nb}         Galactic latitude [rad]
         |  |--l          {Nl}         Galactic longitude [rad]
         |  |--T          {NT}         DM kinetic energy [GeV]
         |  |--intensity  {Nb, Nl, NT} DM intensity (flux) [GeV^-1 cm^-2 s^-1 sr^-1]
         |
         |--m_1           Group
         |--m_2           Group
         |  ...
         |--m_{Nm - 1}    Group

      :program:`darkprop` will sample initial conditions based on the interpolation of
      the 3-dimensional `intensity` grid for each DM mass.

   * input.intensity.select_mass ``Array`` (default [])
      An array of **indices** to select a subset of DM particle masses from the
      `intensity_file` for simulation.
      If empty, all masses are selected.

- [input.sample]
   The parameters in this subsection only take effect if ``input.type`` is ``"sample"``.

   * input.sample.sample_file ``String``
      The path to a file containing individual initial conditions.

      The file must be a HDF5 file with the following definitive structure:

      .. code-block:: text

         file
         |--mchi          {Nm}      DM particle mass [GeV]
         |--ran_E         {NE}      DM kinetic energy [GeV]
         |--ran_bl        {NE, 2}   latitude, longitude [rad]
         |--weight        {Nm, NE}  weights of the samples

- [input.source]
   The parameters in this subsection only take effect if ``input.type`` is ``"source"``.

   * input.source.inv_cdf_r_file ``String``
      The path to a file containing the inverse CDF of radius distribution of the source.

      The file should be a space separated CSV file with **two** columns.
      The first column is the CDF ranges from 0 to 1.
      The second column is the radius r, in unit of **cm**.

   * input.source.inv_cdf_T_file ``String``
      The path to a file containing the inverse CDF of kinetic energy distribution of the
      source.

      The file should be a space separated CSV file with **two** columns.
      The first column is the CDF ranges from 0 to 1.
      The second column is the kinetic energy T, in unit of **MeV**.

- [input.solardm]
   The parameters in this subsection only take effect if ``input.dm_model`` is
   ``"SolarDM"``.

   * input.solardm.solar_density_integral_file ``String``
      The path to a file containing the normalized density integral of the Sun.

      The file should be a space separated CSV file with **three** columns.
      The first column is :math:`y`, the second column is :math:`x`, and the third column
      is :math:`z`. The definitions are given in the documentation of the :ref:`Sun` model.

      The file should be grouped in blocks that separated by an empty line.
      In each block the numbers (:math:`y`) in the first column are the same.

   * input.solardm.temperature_scale ``Float`` (default 1.0)
      A scale factor multiplied to the solar temperature.

[output]
^^^^^^^^
- outdir ``String``
   The directory saving all outputs.
   The directory will be created if it does not exist.

- overwrite ``Boolean`` (default ``false``)
   Choosing whether to overwrite existing output files.

- file_name_style ``String`` (default "long")
   Naming style of the output files. Options:

   :"long": "{outdir}/{id}/{id}_mchi{mchi:.3e}GeV_sigma{sigma:.3e}cm2.hdf5"
   :"short": "{outdir}/{id}/{id}_m{mi}_s{si}.hdf5" where ``mi`` and ``si`` are the index of
             mchis and sigmas.

- maximal_simulation ``Integer`` (default 1000000000)
   The maximum number of simulated particles.
   Forcibly terminate the simulation even if the number of samples has not reached the
   requirement, and the rest will be ``NaN``.

- sample_size ``Integer``
   The number of samples **at each depth on each isodetection ring**.
   The storage space will be allocated from the beginning based on this value, and filled
   with ``NaN``.

- depth ``Array`` [km]
   The list of depth to collect samples.

- rings ``Integer`` (default 1)
   Number of isodetection rings.
   Any number :math:`\leq 1` will be set to 1.

- cos_theta_lower ``Float`` (default -1.0)
   Lower boundary of the cosine of the ring angle :math:`\Theta`.

- cos_theta_upper ``Float`` (default 1.0)
   Upper boundary of the cosine of the ring angle :math:`\Theta`.

- ring_space ``String`` (default "linear")
   The method of separating the isodetection ring angle :math:`\Theta`. Options:

   :"linear": :math:`\Theta` is uniformly separated from 0 to :math:`2\pi`.
   :"cos_linear": :math:`\cos\Theta` is uniformly separated from 1 to -1.

- axis ``Array`` (default [1.0, 0.0, 0.0])
   The axis of the isodetection ring.
   It takes no effect if ``input.type`` is "halo", in which case the axis will always be
   ``input.halo.v_earth``.

- chunk_size ``Integer`` (default 10000)
   The chunk size for HDF5 output.

- store_time ``Boolean`` (default ``false``)
   Whether to store the time of the events.

- full_tracks ``Integer`` (default 0)
   Number of full trajectories additionally stored in the
   ``./<output.outdir>/<id>/fulltracks`` directory.

- [output.log]
   * screen_log_lag ``Integer`` (default 10000)
      Control the frequency of screen log output.

   * warning_long_track ``Integer`` (default 1000000)
      A warning will be printed once a trajectory is longer than this value.

[numerical]
^^^^^^^^^^^
The parameters in this section control some details of the numerical
calculation, generally using the default values.

- random_seed ``Integer`` (default -1)
   The seed of the random number generator.

   If the seed :math:`\geq 0`, the result will be reproducible.

- interp_Tmin ``Float`` [GeV] (default 1e-100)
    Lower boundary of the interpolation for :ref:`SIHelmDM`.

- interp_num ``Integer`` (default 100000)
    Interpolation grid points for :ref:`SIHelmDM`.

- sample_buffer_size ``Integer`` (default 10000)
   The buffer size used when ``input.type`` is ``"sample"``.


Output file format
------------------

Crossing Events
^^^^^^^^^^^^^^^
The crossing events will be stored in HDF5 files with file names specified in the
``file_name_style`` parameter. The HDF5 files contain datasets (scalar or 1D array or 2D
array) arranged in groups like directories in a file system. The data file is organized as
follows

.. code-block:: text

   file
   |--mchi    {SCALAR}  DM particle mass
   |--sigma   {SCALAR}  DM-nucleon cross section
   |--R0      {SCALAR}  surface radius
   |--Tcut    {SCALAR}  termination kinetic energy of the simulation
   |--Tmin    {SCALAR}  minimal incident energy
   |--Tmax    {SCALAR}  maximal incident energy
   |--V0      {SCALAR}  (halo DM) most probable velocity
   |--Vesc    {SCALAR}  (halo DM) escape velocity
   |--Vcut    {SCALAR}  (halo DM) termination velocity of the simulation
   |--Vmin    {SCALAR}  (halo DM) minimal incident veloctiy
   |--Vmax    {SCALAR}  (halo DM) maximal incident velocity
   |--Vearth  {3}       (halo DM) Earth's velocity
   |--depth   {Ndepth}  array of the depth
   |--depth_0 Group     the first depth
   |  |--R          {SCALAR}  radius of the surface at this depth
   |  |--depth      {SCALAR}  the detph
   |  |--ring_bins  {Nrings}  array of the bins of the isodetection ring
   |  |--ring_0     Group     the first ring
   |  |  |--Nsample {SCALAR}  number of samples
   |  |  |--Nsim    {SCALAR}  number of simulated particles
   |  |  |--weight  {output.sample_size/Inf}    array of weights from IS
   |  |  |--T       {output.sample_size/Inf}    array of kinetic energies
   |  |  |--p       {output.sample_size/Inf, 2} array of momentum directions (b, l)
   |  |  |--r       {output.sample_size/Inf, 2} array of position directions (b, l)
   |  |  |--t       {output.sample_size/Inf}    array of times
   |  |
   |  |--ring_1            Group
   |  |--ring_2            Group
   |  |  ...
   |  |--ring_{Nrings - 1} Group
   |
   |--depth_1              Group
   |--depth_2              Group
   |  ...
   |--depth_{Ndepth - 1}   Group

.. note::

   - If ``output.rings`` <= 1, there will **not** be a ``ring_0`` group but all
     the datasets will be in the ``depth_i`` group directly.
   - ``Nsample`` should be equal to the number of non-NaN data in the datasets, otherwise
     there is a bug or the data is corrupted.
   - Units are stored with each dataset as an attribute.


Full Trajectories
^^^^^^^^^^^^^^^^^
Currently, the full trajectories are stored in plain CSV (space separated) files. The
quantities in each column are :math:`r_x, r_y, r_z, p_x, p_y, p_z, T`, and :math:`t`
respectively. The units are in km, GeV, and second.
