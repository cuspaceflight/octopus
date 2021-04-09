Octopus Documentation
=====================

Octopus is a software package designed primarily to model the flow of nitrous oxide (N2O) through a rocket engine
injector. It has a built-in EOS based on a Helmholz energy method outlined in `[3]<main ref 3>`_, and includes the
required coefficients to use it with nitrous oxide. 19 other chemicals are also available, but CO2 is not.
:class:`octopus.main.Nist` contains methods that may be used, along with :mod:`scipy.optimize` and
:func:`octopus.utils.derivative` to produce the required coefficients for the modelling of CO2, if required.

.. toctree::
    installation
    usage
    _autosummary/octopus.main
    _autosummary/octopus.utils

.. autosummary::
    :recursive:

    octopus.main
    octopus.utils

Octopus.main module
===================

.. automodule:: octopus.main
   :members:
   :undoc-members:
   :show-inheritance:
   :member-order: bysource

   .. rubric:: Classes

   .. autosummary::

      Element
      Fluid
      Manifold
      Orifice

Octopus.utils module
====================

.. automodule:: octopus.utils
   :members:
   :undoc-members:
   :show-inheritance:
   :member-order: bysource

   .. rubric:: Functions

   .. autosummary::

      derivative

   .. rubric:: Classes

   .. autosummary::

      Nist

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`