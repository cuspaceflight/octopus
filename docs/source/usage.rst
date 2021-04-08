Usage
=====

Octopus was designed primarily to model nitrous oxide (N2O) flow through a rocket injector. To start working with the
module, a :class:`octopus.main.Fluid` must be initiated. Throughout this example we will be using nitrous oxide, with a Helmholz
EOS, and IPA, using :mod:`thermo`'s EOS as it is a better-behaved fluid.

>>> from octopus import Fluid
>>> nitrous_oxide = Fluid('N2O',method='helmholz')
>>> isopropanol = Fluid('isopropanol',P=18e5,T=293)

If using the helmholz EOS, P and T may not be specified (they are not used, as they are thermo-inherited properties,
however the ``method`` parameter must be specified as "helmholz" as "thermo" is the default. Currently, "helmholz" is
only available for N2O, but it will not throw an error if you add your own data for another fluid.

Next, we can define a :class:`octopus.main.PropertySource` object to represent the initial properties of ``nitrous_oxide``. A
PropertySource is simply an object that has ``p`` and ``T`` properties (note: lowercase p), which the :class:`octopus.main.Manifold`
class uses to define its inlet conditions. The :class:`octopus.main.Fluid` object itself is not used, as its properties change
throughout the system, yet the same object is passed to all classes.

>>> from octopus import PropertySource, Manifold
>>> nitrous_manifold = Manifold(nitrous_oxide,PropertySource(p=18e5,T=250))
>>> ipa_manifold = Manifold(isopropanol,PropertySource(p=18e5,T=293))

Next we can define :class:`octopus.main.Orifice`s that draw fluid from the manifolds defined above, and an Element to
connect the two.

>>> from octopus import Orifice,Element
>>> nitrous_orifice = Orifice(nitrous_manifold,L=1e-2,D=2e-3,orifice_type=Orifice.STRAIGHT)
>>> ipa_orifice = Orifice(ipa_manifold,L=1e-2,D=1e-3,orifice_type=Orifice.STRAIGHT)
>>> n2o_ipa_element = Element(nitrous_orifice,ipa_orifice,downstream=PropertySource(p=15e5)) # in progress




