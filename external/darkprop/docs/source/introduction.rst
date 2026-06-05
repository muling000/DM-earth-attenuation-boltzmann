Introduction
============

DarkProp is a Monte Carlo simulation code for the propagation of dark matter particles in a
medium. It was originally developed for simulating the Earth attenuation effect of the
cosmic ray boosted dark matter (CRDM).

Overview of the simulation
--------------------------

The simulation of dark matter (DM) propagation includes the following procedures:

1. Prepare a DM particle with momentum and position sampled from specific initial
   conditions.

2. Repeat propagation and scattering processes

    * propagation:
       Sample a length :math:`l` according to the mean free path. Then the DM particle
       travel the distance :math:`l` in a straight line along the momentum direction.
    * scattering:
       Sample the momentum of the DM particle after scattering using the differential
       cross section.

3. Terminate the simulation based on some conditions, for example

    * Particle escapes the medium.
    * Particle's kinetic energy is below a certain threshold.

Once the trajectories of the particles are obtained, we can do any interesting analysis.

The two core steps of the simulation are propagation and scattering, which are explained in
the next two sections.


Propagation
^^^^^^^^^^^

Suppose we have a beam of DM particles with intensity :math:`I` traveling inside the Earth.
The number density of the crustal component nucleus :math:`N` is :math:`n_N`. We are
interested in how many particles can **freely** travel a distance longer than :math:`L`.
The change in intensity :math:`\mathrm{d}I` when the particles travel a small distance
:math:`\mathrm{d}z` is

.. math::
   \mathrm{d}I = - I \sum_N n_N \sigma_{\chi N}^\mathrm{tot} \mathrm{d}z,
   :label: eq:dIdz

where

.. math::
   \sigma_{\chi N}^{\mathrm{tot}}
   \equiv \int_0^{q^2_{\max}}
   \frac{\mathrm{d}\sigma_{\chi N}}{\mathrm{d}q^{2}}\mathrm{d}q^2
   :label: eq:sigma_tot

is the total cross section of DM-nucleus (or DM-electron) scattering, where :math:`q` is
the momentum transfer. The cross section may depend on the energy of the incident DM
particle, and :math:`n_N` may vary with location. The solution of :eq:`eq:dIdz` is

.. math::
   I = I_0 \exp(-\int_0^L \frac{\mathrm{d}z}{\lambda}),
   :label: eq:II0

where :math:`I_0` is the intensity at :math:`L = 0`, and we defined the **mean free path**
:math:`\lambda` as

.. math::
   \lambda^{-1} = \sum_N \lambda_N^{-1}
   \equiv \sum_N n_N \sigma_{\chi N}^{\mathrm{tot}}.
   :label: eq:mean_free_path


The key point is that :math:`I / I_0` can be considered as the probability that a single
particle freely passes a distance greater than :math:`L`, therefore

.. math::
   P(\hat{l} > L) = \exp\left(-\int_0^L \frac{\mathrm{d}z}{\lambda}\right),
   :label: eq:prob_l

where :math:`\hat{l}` denotes the random variable of distance.
Then the cumulative distribution function (CDF) is

.. math::
   F(L) \equiv P(\hat{l} \leqslant L)
   = 1 - \exp\left(-\int_0^L \frac{\mathrm{d}z}{\lambda}\right).
   :label: eq:cdf_l

We can draw a sample of :math:`\hat{l}` by the inverse CDF method

.. math::
   \hat{l} = F^{-1}(\xi),
   :label: eq:sample_free_path_general

where :math:`\xi` is a random number drawn from a uniform distribution between 0 and 1.
So the question of sampling a free path :math:`L` reduces to solving :math:`L` from the
equation

.. math::
   \int_0^L \frac{\mathrm{d}z}{\lambda} = -\ln(1-\xi).
   :label: eq:sample_free_path

The specific form of :math:`\lambda` depends on the cross section and the Earth model.

Scattering
^^^^^^^^^^

Before scattering, we need to determine what kind of nucleus the DM particle collides with.
This can be sampled using the probability of scattering on species :math:`N`

.. math::
   P_N = \frac{n_N \sigma_{\chi N}^{\mathrm{tot}}}{\sum_{N'}n_{N'}
   \sigma_{\chi N'}^{\mathrm{tot}}}
   = \frac{\lambda}{\lambda_N}.
   :label: eq:prob_scatter

The probability density function (PDF) of momentum transfer squared :math:`q^2` is the
normalized differential cross section

.. math::
    p(q^2) = \frac{1}{\sigma_{\chi N}^{\mathrm{tot}}}
   \frac{\mathrm{d}\sigma_{\chi N}}{\mathrm{d}q^{2}}.
   :label: eq:prob_q2

The :math:`q^2` and the corresponding recoil energy :math:`T_r` of the target nucleus can
be sampled from this distribution.

For elastic scattering, the scattering angle (between :math:`\vec{p}_i` and
:math:`\vec{p}_f` in laboratory frame) can be obtained from kinematics

.. math::
   \cos\theta = \frac{p_i^2 - T_r (T_i + m_\chi + m_N)}{p_i p_f},
   :label: eq:scattering_angle

where subscripts :math:`i` and :math:`j` stand for the initial and final states,
respectively.

Implementation of the Code
--------------------------

We designed two abstract base classes to implement above processes, the
:ref:`darkprop::Particle <class Particle>` class and the
:ref:`darkprop::Medium <class Medium>` class. General speaking, the ``Particle`` class
implements the total cross section in :eq:`eq:sigma_tot`, the sampling of the recoil energy
using :eq:`eq:prob_q2`, and the calculation of the scattering angle from the sampled recoil
energy (e.g. :eq:`eq:scattering_angle` for elastic scattering). The ``Medium`` class
implements the sampling of the free path in :eq:`eq:sample_free_path_general` and the
scattering probability in :eq:`eq:prob_scatter`. Using such an abstraction, various DM and
medium models can be implemented by inheriting these two base classes.
See :doc:`basic-interface` and :doc:`implemented-models` sections for details.
