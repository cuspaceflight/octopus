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
>>> n2o_ipa_element = Element([nitrous_orifice],[ipa_orifice]) # in progress

If we simply want to use the EOS from the :class:`Fluid` class, there are many methods available if
``method='helmholz'`` is set. The next example will show how to get the saturation properties at 250K. We assume an all-
liquid state, and never use the bare property methods (:meth:`Fluid.p()`) as they may return invalid results.

To get one of ``p,chi,h,s`` as a function of ``[rho,T]``:

>>> T = 250
>>> rho = nitrous.rho_l(T)
>>> properties = nitrous.get_properties(rho,T)
>>> p = properties['p']

To access a property, we can use the properties dictionary to return it using its name as a key. In order to compute rho
or T, we need a non-linear solver:

>>> from scipy.optimize import least_squares

To get ``[rho,T]`` as a funtion of two of ``p,chi,h,s`` (see :meth:`octopus.main.Fluid.fun_ps`):

>>> p=18e5          # pressure we want
>>> chi=0           # vapour fraction we want (all liquid)
>>> x0 = [800,250]  # initial values of rho and T near to the answer
>>> u = ['p','chi'] # names of dependent variables we know
>>> y = [p,chi]     # values of dependent variables we know (same order)
>>> properties = least_squares(nitrous.fun_ps,x0,args=(u,y))
>>> rho,T = properties.x

The above code is used in :meth:`octopus.main.Orifice.m_dot_HEM` to calculate the inital and final conditions, and hence
that is a good example to look to for further context.

DOES NOT WORK - SHOWN AS A TODO
===============================
To get density and vapour fraction as a function of T and one of ``p,chi,h,s``, i.e. average density and vapour fraction
in a 5MPa tank at ambient temperature:

>>> T = 293         # temperature we want
>>> p = 50e5      # vapour fraction we want
>>> x0 = [800,T]    # guess of density, and known temperature
>>> u = ['T','p'] # T and chi are known
>>> y = [T,p]
>>> properties = least_squares(nitrous.fun_ps,x0,args=(u,y))
>>> rho = properties.x[0]
>>> chi = nitrous.get_properties(rho,T)['p']




