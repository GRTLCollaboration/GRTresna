/* GRTresna
 * Copyright 2024 The GRTL Collaboration.
 * Please refer to LICENSE in GRTresna's root directory.
 */

#include "PsiAndAijFunctions.hpp"
#include "REAL.H"
#include "RealVect.H"
#include "Tensor.hpp"
#include "TensorAlgebra.hpp"

void PsiAndAijFunctions::read_params(GRParmParse &pp,
                                     params_t &a_psi_and_Aij_params)
{
    // Initial conditions for the black holes
    pp.get("bh1_bare_mass", a_psi_and_Aij_params.bh1_bare_mass);
    pp.get("bh2_bare_mass", a_psi_and_Aij_params.bh2_bare_mass);
    std::vector<double> temp_spin1(SpaceDim);
    std::vector<double> temp_spin2(SpaceDim);
    std::vector<double> temp_offset1(SpaceDim);
    std::vector<double> temp_offset2(SpaceDim);
    std::vector<double> temp_mom1(SpaceDim);
    std::vector<double> temp_mom2(SpaceDim);
    pp.getarr("bh1_spin", temp_spin1, 0, SpaceDim);
    pp.getarr("bh2_spin", temp_spin2, 0, SpaceDim);
    pp.getarr("bh1_offset", temp_offset1, 0, SpaceDim);
    pp.getarr("bh2_offset", temp_offset2, 0, SpaceDim);
    pp.getarr("bh1_momentum", temp_mom1, 0, SpaceDim);
    pp.getarr("bh2_momentum", temp_mom2, 0, SpaceDim);
    for (int idir = 0; idir < SpaceDim; idir++)
    {
        a_psi_and_Aij_params.bh1_spin[idir] = temp_spin1[idir];
        a_psi_and_Aij_params.bh2_spin[idir] = temp_spin2[idir];
        a_psi_and_Aij_params.bh1_offset[idir] = temp_offset1[idir];
        a_psi_and_Aij_params.bh2_offset[idir] = temp_offset2[idir];
        a_psi_and_Aij_params.bh1_momentum[idir] = temp_mom1[idir];
        a_psi_and_Aij_params.bh2_momentum[idir] = temp_mom2[idir];
    }

    if (abs(a_psi_and_Aij_params.bh1_bare_mass) > 0.0 ||
        abs(a_psi_and_Aij_params.bh2_bare_mass) > 0.0)
    {
        pout() << "Spacetime contains black holes with bare masses "
               << a_psi_and_Aij_params.bh1_bare_mass << " and "
               << a_psi_and_Aij_params.bh2_bare_mass << endl;
    }

    pp.get("use_compact_Vi_ansatz", a_psi_and_Aij_params.use_compact_Vi_ansatz,
           false);
}

void PsiAndAijFunctions::get_bh_coords(Real &bh_radius, RealVect &loc_bh,
                                       const RealVect &loc,
                                       const RealVect &bh_offset)
{
    // set coords
    loc_bh = loc - bh_offset;

    // set radius
    Real bh_radius_squared = 0.0;
    FOR1(i) { bh_radius_squared += loc_bh[i] * loc_bh[i]; }
    bh_radius = sqrt(bh_radius_squared);
}

Real PsiAndAijFunctions::compute_bowenyork_psi(const RealVect &loc)
{
    // the Bowen York params
    Real m1 = m_psi_and_Aij_params.bh1_bare_mass;
    Real m2 = m_psi_and_Aij_params.bh2_bare_mass;
    RealVect J1 = m_psi_and_Aij_params.bh1_spin;
    RealVect J2 = m_psi_and_Aij_params.bh2_spin;

    Real J1_squared = 0.0;
    Real J2_squared = 0.0;
    FOR1(i)
    {
        J1_squared += J1[i] * J1[i];
        J2_squared += J2[i] * J2[i];
    }
    Real a1_squared = 0.0;
    Real a2_squared = 0.0;
    if (m1 > 0.0)
    {
        a1_squared = J1_squared / m1 / m1;
    }
    if (m2 > 0.0)
    {
        a2_squared = J2_squared / m2 / m2;
    }

    // set the BH values - location
    RealVect loc_bh1;
    Real rbh1;
    get_bh_coords(rbh1, loc_bh1, loc, m_psi_and_Aij_params.bh1_offset);

    RealVect loc_bh2;
    Real rbh2;
    get_bh_coords(rbh2, loc_bh2, loc, m_psi_and_Aij_params.bh2_offset);

    // equation (22) in arXiv:gr-qc/0612001

    return 0.5 * (sqrt(m1 * m1 - a1_squared) / rbh1 +
                  sqrt(m2 * m2 - a2_squared) / rbh2);
}

// Set Aij Bowen York data
// see Alcubierre pg 110 eqn (3.4.22)
void PsiAndAijFunctions::compute_bowenyork_Aij(
    Tensor<2, Real> &Aij, // const IntVect &iv,
    const RealVect &loc)
{
    // set the BH values - location
    RealVect loc_bh1;
    Real rbh1;
    get_bh_coords(rbh1, loc_bh1, loc, m_psi_and_Aij_params.bh1_offset);

    RealVect loc_bh2;
    Real rbh2;
    get_bh_coords(rbh2, loc_bh2, loc, m_psi_and_Aij_params.bh2_offset);

    RealVect n1 = {loc_bh1[0] / rbh1, loc_bh1[1] / rbh1, loc_bh1[2] / rbh1};
    RealVect n2 = {loc_bh2[0] / rbh2, loc_bh2[1] / rbh2, loc_bh2[2] / rbh2};

    // the Bowen York params
    RealVect J1 = m_psi_and_Aij_params.bh1_spin;
    RealVect J2 = m_psi_and_Aij_params.bh2_spin;
    RealVect P1 = m_psi_and_Aij_params.bh1_momentum;
    RealVect P2 = m_psi_and_Aij_params.bh2_momentum;

    using namespace TensorAlgebra;
    Tensor<3, Real> epsilon = TensorAlgebra::epsilon();

    FOR2(i, j)
    {
        Aij[i][j] = 1.5 / rbh1 / rbh1 * (n1[i] * P1[j] + n1[j] * P1[i]) +
                    1.5 / rbh2 / rbh2 * (n2[i] * P2[j] + n2[j] * P2[i]);

        FOR1(k)
        {
            Aij[i][j] += 1.5 / rbh1 / rbh1 * (n1[i] * n1[j] - delta(i, j)) *
                             P1[k] * n1[k] +
                         1.5 / rbh2 / rbh2 * (n2[i] * n2[j] - delta(i, j)) *
                             P2[k] * n2[k];

            FOR1(l)
            {
                Aij[i][j] +=
                    -3.0 / rbh1 / rbh1 / rbh1 *
                        (epsilon[i][l][k] * n1[j] + epsilon[j][l][k] * n1[i]) *
                        n1[l] * J1[k] -
                    3.0 / rbh2 / rbh2 / rbh2 *
                        (epsilon[i][l][k] * n2[j] + epsilon[j][l][k] * n2[i]) *
                        n2[l] * J2[k];
            }
        }
    }
}

// The part of Aij excluding the Brill Lindquist BH Aij
// Using ansatz in B&S Appendix B Eq B.5
void PsiAndAijFunctions::compute_ctt_Aij(Tensor<2, Real> &Aij,
                                         const FArrayBox &multigrid_vars_box,
                                         const IntVect &iv,
                                         const RealVect &a_dx,
                                         const RealVect &loc) const
{
    DerivativeOperators derivs(a_dx);

    // get the derivs
    Tensor<2, Real, SpaceDim> d2_U;
    derivs.get_d2(d2_U, iv, multigrid_vars_box, c_U_0);
    Tensor<2, Real, SpaceDim> d1_Vi;
    derivs.get_d1_vector(d1_Vi, iv, multigrid_vars_box,
                         Interval(c_V1_0, c_V3_0));
    Tensor<3, Real, SpaceDim> d2_Vi;
    derivs.get_d2_vector(d2_Vi, iv, multigrid_vars_box,
                         Interval(c_V1_0, c_V3_0));

    // Periodic: Use ansatz B.3 in B&S (p547)
    // Non-periodic: Compact ansatz B.7 in B&S (p547)
    Real trace = 0.0;
    if (!m_psi_and_Aij_params.use_compact_Vi_ansatz)
    {
        FOR1(i) { trace += d1_Vi[i][i] + d2_U[i][i]; }

        // set the values of Aij
        FOR2(i, j)
        {
            Aij[i][j] = d1_Vi[i][j] + d2_U[i][j] + d1_Vi[j][i] + d2_U[j][i] -
                        2.0 / 3.0 * TensorAlgebra::delta(i, j) * trace;
        }
    }
    else
    {
        FOR1(i)
        {
            trace +=
                0.75 * d1_Vi[i][i] -
                0.125 * (d2_U[i][i] + loc[0] * d2_Vi[0][i][i] +
                         loc[1] * d2_Vi[1][i][i] + loc[2] * d2_Vi[2][i][i]);
        }
        // set the values of Aij
        FOR2(i, j)
        {
            Aij[i][j] = 0.75 * d1_Vi[i][j] - 0.125 * d2_U[i][j] +
                        0.75 * d1_Vi[j][i] - 0.125 * d2_U[j][i] -
                        2.0 / 3.0 * TensorAlgebra::delta(i, j) * trace;
            FOR1(k)
            {
                Aij[i][j] -=
                    0.125 * (loc[k] * (d2_Vi[k][i][j] + d2_Vi[k][j][i]));
            }
        }
    }
}
