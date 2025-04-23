/* GRTresna
 * Copyright 2024 The GRTL collaboration.
 * Please refer to LICENSE in GRTresna's root directory.
 */

#include "ScalarField.hpp"

Real ScalarField::my_potential_function(const Real &phi_here) const
{
    return 0.5 * pow(m_matter_params.scalar_mass * phi_here, 2.0);
}

Real ScalarField::my_phi_function(const RealVect &loc,
                                  const RealVect &a_dx) const
{
    // Just for readability
    RealVect L = domainLength;
    int num_lines = m_matter_params.lines;
    Real spacing = m_matter_params.spacing;
    Real *input_dphi = m_matter_params.input_dphi;

    // where am i?
    Real xx = loc[0] + L[0] / 2;
    Real yy = loc[1] + L[1] / 2;
    Real zz = loc[2] + L[2] / 2;

    // Periodic BCs:
    if (xx < 0)
    {
        xx += L[0];
    }
    if (yy < 0)
    {
        yy += L[1];
    }
    if (zz < 0)
    {
        zz += L[2];
    }

    if (xx > L[0])
    {
        xx += -L[0];
    }
    if (yy > L[1])
    {
        yy += -L[1];
    }
    if (zz > L[2])
    {
        zz += -L[2];
    }

    int x_H = 0;
    int y_H = 0;
    int z_H = 0;
    if (xx > L[0] - a_dx[0] && xx < L[0])
    {
        // xx -= L[0];
        x_H = num_lines;
    }
    if (yy > L[1] - a_dx[1] && xx < L[1])
    {
        // yy -= L[1];
        y_H = num_lines;
    }
    if (zz > L[2] - a_dx[2] && xx < L[2])
    {
        // zz -= L[2];
        z_H = num_lines;
    }

    // https://www.wikiwand.com/en/Trilinear_interpolation#/Method

    int i_L = static_cast<int>(floor(xx / spacing));
    int i_H = static_cast<int>(ceil(xx / spacing));
    int j_L = static_cast<int>(floor(yy / spacing));
    int j_H = static_cast<int>(ceil(yy / spacing));
    int k_L = static_cast<int>(floor(zz / spacing));
    int k_H = static_cast<int>(ceil(zz / spacing));

    Real x0 = i_L * spacing;
    Real x1 = i_H * spacing;
    Real y0 = j_L * spacing;
    Real y1 = j_H * spacing;
    Real z0 = k_L * spacing;
    Real z1 = k_H * spacing;

    Real xd = (xx - x0) / (x1 - x0);
    Real yd = (yy - y0) / (y1 - y0);
    Real zd = (zz - z0) / (z1 - z0);

    if (xx > L[0] - a_dx[0] && xx < L[0])
    {
        i_H += -num_lines;
    }
    if (yy > L[1] - a_dx[1] && xx < L[1])
    {
        j_H += -num_lines;
    }
    if (zz > L[2] - a_dx[2] && xx < L[2])
    {
        k_H += -num_lines;
    }
    Real c000 =
        *(input_dphi + num_lines * num_lines * i_L + num_lines * j_L + k_L);
    Real c100 =
        *(input_dphi + num_lines * num_lines * i_H + num_lines * j_L + k_L);
    Real c010 =
        *(input_dphi + num_lines * num_lines * i_L + num_lines * j_H + k_L);
    Real c110 =
        *(input_dphi + num_lines * num_lines * i_H + num_lines * j_H + k_L);
    Real c001 =
        *(input_dphi + num_lines * num_lines * i_L + num_lines * j_L + k_H);
    Real c101 =
        *(input_dphi + num_lines * num_lines * i_H + num_lines * j_L + k_H);
    Real c011 =
        *(input_dphi + num_lines * num_lines * i_L + num_lines * j_H + k_H);
    Real c111 =
        *(input_dphi + num_lines * num_lines * i_H + num_lines * j_H + k_H);

    // Interpolate in x
    Real c00 = c000 * (1 - xd) + c100 * xd;
    Real c01 = c001 * (1 - xd) + c101 * xd;
    Real c10 = c010 * (1 - xd) + c110 * xd;
    Real c11 = c011 * (1 - xd) + c111 * xd;

    // Interpolate in y
    Real c0 = c00 * (1 - yd) + c10 * yd;
    Real c1 = c01 * (1 - yd) + c11 * yd;

    // Interpolate in z
    Real dphi_value = c0 * (1 - zd) + c1 * zd;

    return m_matter_params.phi_0 + m_matter_params.dphi * dphi_value;
}

Real ScalarField::my_Pi_function(const RealVect &loc,
                                 const RealVect &a_dx) const
{
    // Just for readability
    RealVect L = domainLength;
    int num_lines = m_matter_params.lines;
    Real spacing = m_matter_params.spacing;
    Real *input_dpi = m_matter_params.input_dpi;

    // where am i?
    Real xx = loc[0] + L[0] / 2;
    Real yy = loc[1] + L[1] / 2;
    Real zz = loc[2] + L[2] / 2;

    // Periodic BCs:
    if (xx < 0)
    {
        xx += L[0];
    }
    if (yy < 0)
    {
        yy += L[1];
    }
    if (zz < 0)
    {
        zz += L[2];
    }

    if (xx > L[0])
    {
        xx += -L[0];
    }
    if (yy > L[1])
    {
        yy += -L[1];
    }
    if (zz > L[2])
    {
        zz += -L[2];
    }

    int x_H = 0;
    int y_H = 0;
    int z_H = 0;
    if (xx > L[0] - a_dx[0] && xx < L[0])
    {
        // xx -= L[0];
        x_H = num_lines;
    }
    if (yy > L[1] - a_dx[1] && xx < L[1])
    {
        // yy -= L[1];
        y_H = num_lines;
    }
    if (zz > L[2] - a_dx[2] && xx < L[2])
    {
        // zz -= L[2];
        z_H = num_lines;
    }

    // https://www.wikiwand.com/en/Trilinear_interpolation#/Method

    int i_L = static_cast<int>(floor(xx / spacing));
    int i_H = static_cast<int>(ceil(xx / spacing));
    int j_L = static_cast<int>(floor(yy / spacing));
    int j_H = static_cast<int>(ceil(yy / spacing));
    int k_L = static_cast<int>(floor(zz / spacing));
    int k_H = static_cast<int>(ceil(zz / spacing));

    Real x0 = i_L * spacing;
    Real x1 = i_H * spacing;
    Real y0 = j_L * spacing;
    Real y1 = j_H * spacing;
    Real z0 = k_L * spacing;
    Real z1 = k_H * spacing;

    Real xd = (xx - x0) / (x1 - x0);
    Real yd = (yy - y0) / (y1 - y0);
    Real zd = (zz - z0) / (z1 - z0);

    if (xx > L[0] - a_dx[0] && xx < L[0])
    {
        i_H += -num_lines;
    }
    if (yy > L[1] - a_dx[1] && xx < L[1])
    {
        j_H += -num_lines;
    }
    if (zz > L[2] - a_dx[2] && xx < L[2])
    {
        k_H += -num_lines;
    }
    Real c000 =
        *(input_dpi + num_lines * num_lines * i_L + num_lines * j_L + k_L);
    Real c100 =
        *(input_dpi + num_lines * num_lines * i_H + num_lines * j_L + k_L);
    Real c010 =
        *(input_dpi + num_lines * num_lines * i_L + num_lines * j_H + k_L);
    Real c110 =
        *(input_dpi + num_lines * num_lines * i_H + num_lines * j_H + k_L);
    Real c001 =
        *(input_dpi + num_lines * num_lines * i_L + num_lines * j_L + k_H);
    Real c101 =
        *(input_dpi + num_lines * num_lines * i_H + num_lines * j_L + k_H);
    Real c011 =
        *(input_dpi + num_lines * num_lines * i_L + num_lines * j_H + k_H);
    Real c111 =
        *(input_dpi + num_lines * num_lines * i_H + num_lines * j_H + k_H);

    // Interpolate in x
    Real c00 = c000 * (1 - xd) + c100 * xd;
    Real c01 = c001 * (1 - xd) + c101 * xd;
    Real c10 = c010 * (1 - xd) + c110 * xd;
    Real c11 = c011 * (1 - xd) + c111 * xd;

    // Interpolate in y
    Real c0 = c00 * (1 - yd) + c10 * yd;
    Real c1 = c01 * (1 - yd) + c11 * yd;

    // Interpolate in z
    Real dpi_value = c0 * (1 - zd) + c1 * zd;

    return m_matter_params.pi_0 + m_matter_params.dpi * dpi_value;
}
