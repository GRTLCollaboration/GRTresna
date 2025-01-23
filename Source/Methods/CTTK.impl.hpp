/* GRTresna
 * Copyright 2024 The GRTL Collaboration.
 * Please refer to LICENSE in GRTresna's root directory.
 */

#ifndef CTTK_HPP_
#error "This file should only be included through CTTK.hpp"
#endif

#include "DimensionDefinitions.hpp"
#include "GRParmParse.hpp"
#include "Tensor.hpp"
#include "TensorAlgebra.hpp"

template <typename matter_t> struct CTTK<matter_t>::params_t
{
    int sign_of_K;
    bool use_compact_Vi_ansatz;
    Real regularised_part_psi;
    bool deactivate_zero_mode;
};

template <typename matter_t>
CTTK<matter_t>::CTTK(params_t a_method_params, matter_t *a_matter,
                     PsiAndAijFunctions *a_psi_and_Aij_functions,
                     int a_numLevels,
                     const std::array<double, SpaceDim> a_center,
                     Real a_G_Newton)
    : m_method_params(a_method_params), matter(a_matter), G_Newton(a_G_Newton),
      numLevels(a_numLevels), psi_and_Aij_functions(a_psi_and_Aij_functions),
      center(a_center)
{
}

template <typename matter_t>
void CTTK<matter_t>::read_params(GRParmParse &pp, params_t &a_method_params)
{
    pp.load("sign_of_K", a_method_params.sign_of_K, 1);
    pp.load("use_compact_Vi_ansatz", a_method_params.use_compact_Vi_ansatz,
            false);
    pp.load("regularised_part_psi", a_method_params.regularised_part_psi, 1.0);
    pp.load("deactivate_zero_mode", a_method_params.deactivate_zero_mode,
            false);
}

template <typename matter_t>
void CTTK<matter_t>::solve_analytic(LevelData<FArrayBox> *a_multigrid_vars,
                                    LevelData<FArrayBox> *a_rhs,
                                    const RealVect &a_dx)
{
    DerivativeOperators derivs(a_dx);
    // Iterate through the boxes in turn
    DataIterator dit = a_rhs->dataIterator();
    for (dit.begin(); dit.ok(); ++dit)
    {
        FArrayBox &multigrid_vars_box = (*a_multigrid_vars)[dit()];
        FArrayBox &rhs_box = (*a_rhs)[dit()];
        Box unghosted_box = rhs_box.box();

        // Iterate through the interior of boxes
        // (ghosts need to be filled later due to gradient terms)
        BoxIterator bit(unghosted_box);
        for (bit.begin(); bit.ok(); ++bit)
        {
            // work out location on the grid
            IntVect iv = bit();
            RealVect loc;
            Grids::get_loc(loc, iv, a_dx, center);

            // Calculate the actual value of psi including BH part
            Real psi_reg = multigrid_vars_box(iv, c_psi_reg);
            Real psi_bh = psi_and_Aij_functions->compute_bowenyork_psi(loc);
            Real psi_0 = psi_reg + psi_bh;
            Real laplacian_psi_reg;
            derivs.scalar_Laplacian(laplacian_psi_reg, iv, multigrid_vars_box,
                                    c_psi_reg);

            // Assign values of Aij
            Tensor<2, Real> Aij_reg;
            psi_and_Aij_functions->compute_ctt_Aij(Aij_reg, multigrid_vars_box,
                                                   iv, a_dx, loc);
            Tensor<2, Real> Aij_bh;
            psi_and_Aij_functions->compute_bowenyork_Aij(Aij_bh, loc);
            // This is \bar  A_ij \bar A^ij
            Real A2_0 = 0.0;
            FOR2(i, j)
            {
                A2_0 += (Aij_reg[i][j] + Aij_bh[i][j]) *
                        (Aij_reg[i][j] + Aij_bh[i][j]);
            }

            // Compute emtensor components
            const auto emtensor =
                matter->compute_emtensor(iv, a_dx, multigrid_vars_box);

            // Now work out K using ansatz which sets it to (roughly)
            // the FRW value based on the local densities
            Real K_0_squared = 24.0 * M_PI * G_Newton * emtensor.rho +
                               1.5 * A2_0 * pow(psi_0, -12.0) +
                               12.0 * laplacian_psi_reg * pow(psi_0, -5.0);

            // Set value for K
            // be careful if at a point K = 0, may have discontinuity
            multigrid_vars_box(iv, c_K_0) =
                m_method_params.sign_of_K * sqrt(K_0_squared);

            // set values for \bar Aij_0
            multigrid_vars_box(iv, c_A11_0) = Aij_reg[0][0] + Aij_bh[0][0];
            multigrid_vars_box(iv, c_A22_0) = Aij_reg[1][1] + Aij_bh[1][1];
            multigrid_vars_box(iv, c_A33_0) = Aij_reg[2][2] + Aij_bh[2][2];
            multigrid_vars_box(iv, c_A12_0) = Aij_reg[0][1] + Aij_bh[0][1];
            multigrid_vars_box(iv, c_A13_0) = Aij_reg[0][2] + Aij_bh[0][2];
            multigrid_vars_box(iv, c_A23_0) = Aij_reg[1][2] + Aij_bh[1][2];
        }
    }
}

template <typename matter_t>
void CTTK<matter_t>::set_elliptic_terms(
    LevelData<FArrayBox> *a_multigrid_vars, LevelData<FArrayBox> *a_rhs,
    RefCountedPtr<LevelData<FArrayBox>> a_aCoef,
    RefCountedPtr<LevelData<FArrayBox>> a_bCoef, const RealVect &a_dx)
{
    DerivativeOperators derivs(a_dx);
    DataIterator dit = a_rhs->dataIterator();
    for (dit.begin(); dit.ok(); ++dit)
    {
        FArrayBox &multigrid_vars_box = (*a_multigrid_vars)[dit()];
        FArrayBox &rhs_box = (*a_rhs)[dit()];
        FArrayBox &aCoef_box = (*a_aCoef)[dit()];
        FArrayBox &bCoef_box = (*a_bCoef)[dit()];
        // JCAurre: Initialise rhs=0, aCoef=0 and bCoef=1 for all constraint
        // variables
        for (int comp = 0; comp < NUM_CONSTRAINT_VARS; comp++)
        {
            rhs_box.setVal(0.0, comp);
            aCoef_box.setVal(0.0, comp);
            bCoef_box.setVal(1.0, comp);

            // this prevents small amounts of noise in the sources
            // activating the zero modes - (Garfinkle trick) see 2207.03125
            if (m_method_params.deactivate_zero_mode)
            {
                Real small_number = 1e-10;
                aCoef_box.setVal(-small_number, comp);
            }
        }
        Box unghosted_box = rhs_box.box();
        BoxIterator bit(unghosted_box);
        for (bit.begin(); bit.ok(); ++bit)
        {
            // work out location on the grid
            IntVect iv = bit();
            RealVect loc;
            Grids::get_loc(loc, iv, a_dx, center);

            // Calculate the actual value of psi including BH part
            Real psi_reg = multigrid_vars_box(iv, c_psi_reg);
            Real psi_bh = psi_and_Aij_functions->compute_bowenyork_psi(loc);
            Real psi_0 = psi_reg + psi_bh;
            Real laplacian_psi_reg;
            derivs.scalar_Laplacian(laplacian_psi_reg, iv, multigrid_vars_box,
                                    c_psi_reg);

            // Get values of Aij
            Tensor<2, Real> Aij_reg;
            psi_and_Aij_functions->compute_ctt_Aij(Aij_reg, multigrid_vars_box,
                                                   iv, a_dx, loc);
            Tensor<2, Real> Aij_bh;
            psi_and_Aij_functions->compute_bowenyork_Aij(Aij_bh, loc);
            // This is \bar  A_ij \bar A^ij
            Real A2_0 = 0.0;
            FOR2(i, j)
            {
                A2_0 += (Aij_reg[i][j] + Aij_bh[i][j]) *
                        (Aij_reg[i][j] + Aij_bh[i][j]);
            }

            // Compute emtensor components
            const auto emtensor =
                matter->compute_emtensor(iv, a_dx, multigrid_vars_box);

            Tensor<1, Real, SpaceDim> d1_K;
            derivs.get_d1(d1_K, iv, multigrid_vars_box, c_K_0);

            // Get d_i V_i and laplacians
            Tensor<2, Real, SpaceDim> d1_Vi;
            derivs.get_d1_vector(d1_Vi, iv, multigrid_vars_box,
                                 Interval(c_V1_0, c_V3_0));

            Tensor<1, Real, SpaceDim> laplacian_Vi;
            derivs.vector_Laplacian(laplacian_Vi, iv, multigrid_vars_box,
                                    Interval(c_V1_0, c_V3_0));

            Real laplacian_U;
            derivs.scalar_Laplacian(laplacian_U, iv, multigrid_vars_box, c_U_0);

            // now set the rhs values in the box
            rhs_box(iv, c_psi) = 0.0; // K cancels all terms in CTTK
            rhs_box(iv, c_V1) =
                pow(psi_0, 6.0) *
                (2.0 / 3.0 * d1_K[0] + 8.0 * M_PI * G_Newton * emtensor.Si[0]);
            rhs_box(iv, c_V2) =
                pow(psi_0, 6.0) *
                (2.0 / 3.0 * d1_K[1] + 8.0 * M_PI * G_Newton * emtensor.Si[1]);
            rhs_box(iv, c_V3) =
                pow(psi_0, 6.0) *
                (2.0 / 3.0 * d1_K[2] + 8.0 * M_PI * G_Newton * emtensor.Si[2]);

            // Periodic: Use ansatz B.3 in B&S (p547)
            // Non-periodic: Compact ansatz B.7 in B&S (p547)
            if (!m_method_params.use_compact_Vi_ansatz)
            {
                rhs_box(iv, c_U) =
                    -0.25 * (d1_Vi[0][0] + d1_Vi[1][1] + d1_Vi[2][2]);
            }
            else
            {
                rhs_box(iv, c_U) = 0.0;
                FOR1(i)
                {
                    rhs_box(iv, c_U) +=
                        -loc[i] * (pow(psi_0, 6.0) *
                                   (2.0 / 3.0 * d1_K[i] +
                                    8.0 * M_PI * G_Newton * emtensor.Si[i]));
                }
            }

            if (m_method_params.deactivate_zero_mode)
            {
                rhs_box(iv, c_V1) += -laplacian_Vi[0];
                rhs_box(iv, c_V2) += -laplacian_Vi[1];
                rhs_box(iv, c_V3) += -laplacian_Vi[2];
                rhs_box(iv, c_U) += -laplacian_U;
            }
        }
    }
}

template <typename matter_t>
void CTTK<matter_t>::initialise_method_vars(
    LevelData<FArrayBox> &a_multigrid_vars, const RealVect &a_dx) const
{
    CH_assert(a_multigrid_vars.nComp() == NUM_MULTIGRID_VARS);

    DataIterator dit = a_multigrid_vars.dataIterator();
    for (dit.begin(); dit.ok(); ++dit)
    {
        // These contain the vars in the boxes, set them all to zero
        FArrayBox &multigrid_vars_box = a_multigrid_vars[dit()];
        for (int comp = 0; comp < NUM_MULTIGRID_VARS; comp++)
        {
            multigrid_vars_box.setVal(0.0, comp);
        }

        // Iterate over the box and set non zero comps
        Box ghosted_box = multigrid_vars_box.box();
        BoxIterator bit(ghosted_box);
        for (bit.begin(); bit.ok(); ++bit)
        {

            // work out location on the grid
            IntVect iv = bit();
            RealVect loc;
            Grids::get_loc(loc, iv, a_dx, center);

            // note that we don't include the singular part of psi
            // for the BHs - this is added at the output data stage
            // and when we calculate psi_reg in the rhs etc
            // as it already satisfies Laplacian(psi) = 0
            multigrid_vars_box(iv, c_psi_reg) =
                m_method_params.regularised_part_psi;
        }
    }
}

template <typename matter_t>
void CTTK<matter_t>::initialise_constraint_vars(
    LevelData<FArrayBox> &a_constraint_vars, const RealVect &a_dx) const
{

    DataIterator dit = a_constraint_vars.dataIterator();
    for (dit.begin(); dit.ok(); ++dit)
    {
        // These contain the vars in the boxes, set them all to zero
        FArrayBox &constraint_vars_box = a_constraint_vars[dit()];

        for (int comp = 0; comp < NUM_CONSTRAINT_VARS; comp++)
        {
            constraint_vars_box.setVal(0.0, comp);
        }
    }
}
