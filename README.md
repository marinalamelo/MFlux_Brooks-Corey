# MFlux_Brooks-Corey
This repository refers to the analytical solution of the microscopic root water uptake model developed by Quirijn de Jong van Lier for the special case of Brooks and Corey soils. 
Link to the capsule on Code Ocean: https://doi.org/10.24433/CO.2960439.v1

The software calculates root water uptake and soil and crop potentials according to a process-based model, De Jong van Lier et al. (2013), http://dx.doi.org/10.2136/vzj2013.02.0039 for the special case of Brooks and Corey soils. Scenario parameters are defined in the file input.qvl. The simulation scenario includes root and xylem radii, longitudinal and radial root conductance, and minimum leaf water potential.

Soil hydraulic properties are described by the Brooks-Corey parameters, as well as the root length density. The model output (file Info.out) includes actual (Ta) transpiration rate (as defined in de Jong van Lier et al., 2013), the root water uptake from the soil layer, the pressure head at the surface of the roots in the soil layer, and the xylem and leaf water potentials.
