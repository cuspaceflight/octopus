# ***OCTOPUS***

<img src="img/logo.png" alt="logo" width="600"/>

## A repository for m5 injector design engineers to collate code performing analysis on injector flow and efficacy.

### Install with pip:
`pip install git+https://github.com/cuspaceflight/octopus.git`

### Documentation:
[Octopus Docs](https://cuspaceflight.github.io/octopus)

### Current features:
*   Model flow through a simple straight orifice with a selection of 1- and 2-phase models.
*   Models supported:
 *  Single Phase Incompressible: `python Orifice.m_dot_SPI(Pcc)`
 *  Homogeneous Equilibrium Model: `Orifice.m_dot_HEM(Pcc)`
 *  Dyer Correction Model: `Orifice.m_dot_dyer(Pcc)`
*   Accurate Equations of State for multiple fluids, coefficients included for N2O.
*   Example code showing some features of the `octopus` module, most documentation is in docstrings on source code.

### Currently in progress:
- [ ] Improving class hierarchy, including a Plate class
- [ ] Adding venturi calculations
- [ ] Making nitrous data easier to access
- [ ] Improving documentation
- [ ] GUI for creating candidate injector patterns and analysing their mass flow distribution, including a graphical visualisation

### Eventual planned features:
*   Orifice types - Waxman cavitating
*   Estimation of orifice discharge coefficients (from empirical data) (?)
*   Estimation of atomisation performance (from empirical data) (?)
*   Determine element O/F, overall O/F
*   Determine film cooling mass flow proportion (~20% is said to be sufficient to ignore the core O/F for temperature calculations, need to find the source for this - possibly NASA SP-8089)
*   Warnings for combustion stability criteria (from empirical data)
