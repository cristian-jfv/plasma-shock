Three numerical MHD schemes are implemented to model a schock inside a plasma:
- Lax-Friedrichs (LxF)
- Total Variation Diminishing Lax-Friedrichs (TVD-LxF)
- First Order Arithmetic Mean Matrix (FOAMM)

The conservative form of the MHD equations is used

$$
\frac{\partial \vec{u}}{\partial t} + \frac{\partial \vec{f}}{\partial x} = \vec{0},
$$

taking the discretization, where the flux $\vec{f}$ is approximated at the cells interfaces by the schemes.
