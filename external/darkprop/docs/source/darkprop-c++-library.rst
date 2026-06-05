darkprop C++ Library
====================

The core of DarkProp is a header only template library which can be used by including the
header files. The advantage of using the library directly is that you can flexibly
customize the simulation and extend the physical models through inheritance.

The template parameters are ``Vector`` and ``Value``. The internal ``Vector`` type
``darkprop::Vector3d`` and the ``Value`` type ``double`` will be used by default.
Alternatively, ``Eigen::Vector3d`` in the `Eigen 3 <https://eigen.tuxfamily.org/>`_
library has also been tested.

.. warning::
   Since DarkProp uses the natural unit system with :math:`\text{GeV}\equiv 1`, the
   ``Value`` type *cannot* be ``float``. In general use ``double``. ``long double`` is also
   supported, but it is generally only used in test situations with precision requirements
   due to the lower efficiency.

.. _usage-library:

Usage
-----

The DarkProp library should be available after the cmake installation.
To use the library simply include the needed headers in your source code like this

.. code-block:: cpp

   #include <darkprop/core.hpp>
   #include <darkprop/IO.hpp>  // optional hdf5 io

Below is an working example.

.. _example-library:

Example
^^^^^^^

.. literalinclude:: ../../examples/trajectory-csv/trajectory.cpp
   :linenos:
   :caption: examples/trajectory-csv/trajectory.cpp
   :language: cpp


Simulation Outline
------------------

The input of a simulation is the initial state of the incident particle (initial position,
momentum, kinetic energy, and time, etc.). They can be set using the
:ref:`basic interface<class Particle>` defined in the ``Particle`` base class.
The :ref:`Initialization.hpp` file provides some functions to produce DM with specific
distributions. The :ref:`Injector.hpp` file provides some ``Injectors`` for typical
particle injection scenarios.

The output of DarkProp's built-in :ref:`simulation functions<Simulation.hpp>` is a
``std::vector`` of ``Event`` (``std::vector<Event>``). Users can freely write their
functions to store or analyze these data. DarkProp provides tool functions in :ref:`IO.hpp`
for storing data in HDF5 format. The HDF5 format is suitable for storing large data sets in
scientific computing and widely supported by various programming languages.

To perform a simulation, first create an :ref:`Medium<Medium>` instance and a :ref:`Particle<Particle>` instance (their concrete subclasses), then set the initial state of
the Particle, and then call
:ref:`propagate and scatter<Interface Functions: propagate and scatter>` functions in a
loop, or use the implemented :ref:`simulation functions<Simulation.hpp>`.

DarkProp is designed to change the state of the particle **in place**, instead of storing
the entire trajectory in memory. Therefore, necessary analysis needs to be written in the
simulation function (e.g. to find out the events passing through a sphere) and the output
is also determined by the user.


Class Inheritance
-----------------

.. image:: _static/images/inherit_graph_medium.png
.. image:: _static/images/inherit_graph_particle.png
.. image:: _static/images/inherit_graph_injector.png
