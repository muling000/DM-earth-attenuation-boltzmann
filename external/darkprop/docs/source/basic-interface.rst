Basic Interface
===============
The core of this library is two abstract classes, ``Particle`` and ``Medium``. Users who
want to implement their own dark matter and medium models only need to inherit them and
implement all virtual functions. We show all the abstract classes below.

class Medium
------------

.. doxygenclass:: darkprop::Medium
   :members:

.. doxygenclass:: darkprop::MediumBall
   :members:

.. doxygenclass:: darkprop::Earth
   :members:

class Particle
--------------

.. doxygenclass:: darkprop::Particle
   :members:


Interface Functions: propagate and scatter
------------------------------------------

.. doxygenfunction:: darkprop::propagate

.. doxygenfunction:: darkprop::scatter

Other Ingredients
-----------------

.. doxygenstruct:: darkprop::Target
   :members:

.. doxygenclass:: darkprop::RandomNumber
   :members:

.. doxygenstruct:: darkprop::Event
   :members:
