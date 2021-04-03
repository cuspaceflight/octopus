# ***OCTOPUS***

<img src="img/logo.png" alt="logo" width="600"/>

## A repository for m5 injector design engineers to collate code performing analysis on injector flow and efficacy.

### Current features:
*   Choice of SPI, HEM and Dyer models returning mass flow rate
*   Full Helmholz rho/T EOS implementation in the Fluid class for all relevant properties
*   Orifice types - basic

### Eventual planned features:
*   Orifice types - Waxman cavitating
*   Estimation of orifice discharge coefficients (from empirical data) (?)
*   Estimation of atomisation performance (from empirical data) (?)
*   Determine element O/F, overall O/F
*   Determine film cooling mass flow proportion (~20% is said to be sufficient to ignore the core O/F for temperature calculations, need to find the source for this - possibly NASA SP-8089)
*   Warnings for combustion stability criteria (from empirical data)
*   Document everything - kill two birds with one stone here and write up all the injector notes in preparation for the new wiki?

### Currently in progress:
*   Setting up class heirarchy
  * Current intended structure: Injector plate --> Manifolds --> Injector elements --> Orifices

### Things that need fixing:
*   Not sure the HEM graph produced matches Solomon's

SI base units used everywhere, except for possibly some final output and graphing.
