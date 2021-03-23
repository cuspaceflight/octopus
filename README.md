# ***OCTOPUS***

<img src="img/logo.png" alt="logo" width="600"/>

## A repository for m5 injector design engineers to collate code performing analysis on injector flow and efficacy.

### Eventual planned features:
*   Choice of SPI, HEM and Dyer (with Solomon correction) models
*   Compressibility corrections
*   Orifice types - Waxman cavitating, basic
*   Estimation of orifice discharge coefficients (from empirical data) (?)
*   Estimation of atomisation performance (from empirical data) (?)
*   Determine pressure drop, orifice mass flow, element O/F, overall O/F, etc with the above
*   Determine film cooling mass flow proportion (~20% is said to be sufficient to ignore the core O/F for temperature calculations, need to find the source for this - possibly NASA SP-8089)
*   Warnings for combustion stability criteria (from empirical data)
*   Document everything - kill two birds with one stone here and write up all the injector notes in preparation for the new wiki?

### Currently in progress:
*   Finding specific enthalpy and any other important chemical properties, transport properties
  * Possibly with multiple source data options / averaging / sense checks
*   Decisions on classes / where everything should live
  * Current intended structure: Injector plate --> Manifolds --> Injector elements --> Orifices

### Things that need fixing:
*   Sure this won't stay empty long...

SI base units used everywhere, except for possibly some final output and graphing.