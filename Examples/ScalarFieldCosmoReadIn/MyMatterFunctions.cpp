/* GRTresna
 * Copyright 2024 The GRTL collaboration.
 * Please refer to LICENSE in GRTresna's root directory.
 */

#include "ScalarField.hpp"

Real ScalarField::my_potential_function(const Real &phi_here) const
{
    return 0.5 * pow(m_matter_params.scalar_mass * phi_here, 2.0);
}

Real ScalarField::my_phi_function(const RealVect &loc) const
{
    Real rr = sqrt(loc[0] * loc[0] + loc[1] * loc[1] + loc[2] * loc[2]);
    Real L = domainLength[0];
    Real dphi_value = m_matter_params.dphi / 3. *
                      (sin(2 * M_PI * loc[0] / L) + sin(2 * M_PI * loc[1] / L) +
                       sin(2 * M_PI * loc[2] / L));
    return m_matter_params.phi_0 + dphi_value;
}

Real ScalarField::my_Pi_function(const RealVect &loc) const
{
    Real rr = sqrt(loc[0] * loc[0] + loc[1] * loc[1] + loc[2] * loc[2]);
    Real L = domainLength[0];
    Real dpi_value = m_matter_params.dpi / 3. *
                     (sin(2 * M_PI * loc[0] / L) + sin(2 * M_PI * loc[1] / L) +
                      sin(2 * M_PI * loc[2] / L));
    return m_matter_params.pi_0 + dpi_value;
}
