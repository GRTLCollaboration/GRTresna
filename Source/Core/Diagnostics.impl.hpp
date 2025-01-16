/* GRTresna
 * Copyright 2024 The GRTL Collaboration.
 * Please refer to LICENSE in GRTresna's root directory.
 */

#include "DerivativeOperators.hpp"
#include "DiagnosticVariables.hpp"

template <typename method_t, typename matter_t>
Diagnostics<method_t, matter_t>::Diagnostics(
    method_t *a_method, matter_t *a_matter,
    PsiAndAijFunctions *a_psi_and_Aij_functions, const Real a_G_Newton,
    const std::array<double, SpaceDim> a_center)
    : method(a_method), matter(a_matter),
      psi_and_Aij_functions(a_psi_and_Aij_functions), G_Newton(a_G_Newton),
      center(a_center)
{
}

template <typename method_t, typename matter_t>
void Diagnostics<method_t, matter_t>::compute_constraint_terms(
    LevelData<FArrayBox> *a_multigrid_vars,
    LevelData<FArrayBox> *a_diagnostic_vars, LevelData<FArrayBox> *a_rhs,
    const RealVect &a_dx) const
{
    DerivativeOperators derivs(a_dx);
    // Iterate through the boxes in turn
    DataIterator dit = a_rhs->dataIterator();
    for (dit.begin(); dit.ok(); ++dit)
    {
        FArrayBox &multigrid_vars_box = (*a_multigrid_vars)[dit()];
        FArrayBox &diagnostic_vars_box = (*a_diagnostic_vars)[dit()];
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
            const Real psim6 = 1.0 / pow(psi_0, 6.0);

            Real laplacian_psi_reg;
            derivs.scalar_Laplacian(laplacian_psi_reg, iv, multigrid_vars_box,
                                    c_psi_reg);
            Tensor<1, Real, SpaceDim> d1_K;
            derivs.get_d1(d1_K, iv, multigrid_vars_box, c_K_0);
            Tensor<3, Real, SpaceDim> d2_Vi;
            derivs.get_d2_vector(d2_Vi, iv, multigrid_vars_box,
                                 Interval(c_V1_0, c_V3_0));

            // Assign values of Aij
            Tensor<2, Real> Aij_reg;
            method->psi_and_Aij_functions->compute_ctt_Aij(
                Aij_reg, multigrid_vars_box, iv, a_dx, loc);
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

            diagnostic_vars_box(iv, c_rho) = emtensor.rho;
            diagnostic_vars_box(iv, c_S1) = emtensor.Si[0];
            diagnostic_vars_box(iv, c_S2) = emtensor.Si[1];
            diagnostic_vars_box(iv, c_S3) = emtensor.Si[2];

            // Set value for K
            Real K = multigrid_vars_box(iv, c_K_0);
            Real K_0_squared = K * K;

            diagnostic_vars_box(iv, c_Ham) =
                K_0_squared - 24.0 * M_PI * G_Newton * emtensor.rho -
                1.5 * A2_0 * pow(psi_0, -12.0) -
                12.0 * laplacian_psi_reg * pow(psi_0, -5.0);
            diagnostic_vars_box(iv, c_Ham_abs) =
                K_0_squared + 24.0 * M_PI * G_Newton * emtensor.rho +
                1.5 * abs(A2_0) * pow(psi_0, -12.0) +
                12.0 * abs(laplacian_psi_reg) * pow(psi_0, -5.0);

            Real Mom1 =
                -2.0 / 3.0 * d1_K[0] - 8.0 * M_PI * G_Newton * emtensor.Si[0];
            Real Mom2 =
                -2.0 / 3.0 * d1_K[1] - 8.0 * M_PI * G_Newton * emtensor.Si[1];
            Real Mom3 =
                -2.0 / 3.0 * d1_K[2] - 8.0 * M_PI * G_Newton * emtensor.Si[2];

            Real Mom1_abs = 2.0 / 3.0 * abs(d1_K[0]) +
                            8.0 * M_PI * G_Newton * abs(emtensor.Si[0]);
            Real Mom2_abs = 2.0 / 3.0 * abs(d1_K[1]) +
                            8.0 * M_PI * G_Newton * abs(emtensor.Si[1]);
            Real Mom3_abs = 2.0 / 3.0 * abs(d1_K[2]) +
                            8.0 * M_PI * G_Newton * abs(emtensor.Si[2]);

            FOR(i)
            {
                Mom1 += psim6 * d2_Vi[0][i][i];
                Mom2 += psim6 * d2_Vi[1][i][i];
                Mom3 += psim6 * d2_Vi[2][i][i];

                Mom1_abs += abs(psim6 * d2_Vi[0][i][i]);
                Mom2_abs += abs(psim6 * d2_Vi[1][i][i]);
                Mom3_abs += abs(psim6 * d2_Vi[2][i][i]);
            }

            Real Mom = sqrt(Mom1 * Mom1 + Mom2 * Mom2 + Mom3 * Mom3);

            diagnostic_vars_box(iv, c_Mom1) = Mom1;
            diagnostic_vars_box(iv, c_Mom2) = Mom2;
            diagnostic_vars_box(iv, c_Mom3) = Mom3;
            diagnostic_vars_box(iv, c_Mom) = Mom;
            diagnostic_vars_box(iv, c_Mom1_abs) = Mom1_abs;
            diagnostic_vars_box(iv, c_Mom2_abs) = Mom2_abs;
            diagnostic_vars_box(iv, c_Mom3_abs) = Mom3_abs;
            diagnostic_vars_box(iv, c_Mom_abs) =
                sqrt(Mom1_abs * Mom1_abs + Mom2_abs * Mom2_abs +
                     Mom3_abs * Mom3_abs);
        }
    }
}

template <typename method_t, typename matter_t>
void Diagnostics<method_t, matter_t>::normalise_constraints(
    LevelData<FArrayBox> *a_multigrid_vars,
    LevelData<FArrayBox> *a_diagnostic_vars, LevelData<FArrayBox> *a_rhs,
    const RealVect &a_dx, IntVect &nCells) const
{
    DerivativeOperators derivs(a_dx);
    // Iterate through the boxes in turn
    DataIterator dit = a_rhs->dataIterator();
    for (dit.begin(); dit.ok(); ++dit)
    {
        FArrayBox &multigrid_vars_box = (*a_multigrid_vars)[dit()];
        FArrayBox &diagnostic_vars_box = (*a_diagnostic_vars)[dit()];
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
            Real Ham = diagnostic_vars_box(iv, c_Ham);
            Real Mom1 = diagnostic_vars_box(iv, c_Mom1);
            Real Mom2 = diagnostic_vars_box(iv, c_Mom2);
            Real Mom3 = diagnostic_vars_box(iv, c_Mom3);
            Real Mom = diagnostic_vars_box(iv, c_Mom);

            Real Mom1_abs = diagnostic_vars_box(iv, c_Mom1_abs);
            Real Mom2_abs = diagnostic_vars_box(iv, c_Mom2_abs);
            Real Mom3_abs = diagnostic_vars_box(iv, c_Mom3_abs);

            Real Ham_abs = diagnostic_vars_box(iv, c_Ham_abs);
            Real Mom_abs = diagnostic_vars_box(iv, c_Mom_abs);

            diagnostic_vars_box(iv, c_Ham_norm) = abs(Ham / Ham_abs);
            diagnostic_vars_box(iv, c_Mom_norm) = abs(Mom / Mom_abs);

            // This sets the error norms at the boundary cells to zero
            IntVect lo = IntVect::Zero;
            IntVect hi = nCells - IntVect::Unit;
            if (iv[0] == lo[0] || iv[1] == lo[1] || iv[2] == lo[2] ||
                iv[0] == hi[0] || iv[1] == hi[1] || iv[2] == hi[2])
            {
                diagnostic_vars_box(iv, c_Ham_norm) = 0;
                diagnostic_vars_box(iv, c_Mom_norm) = 0;
            }
        }
    }
}

template <typename method_t, typename matter_t>
void Diagnostics<method_t, matter_t>::initialise_diagnostic_vars(
    LevelData<FArrayBox> &a_diagnostic_vars, const RealVect &a_dx) const
{
    DataIterator dit = a_diagnostic_vars.dataIterator();
    for (dit.begin(); dit.ok(); ++dit)
    {
        // These contain the vars in the boxes, set them all to zero
        FArrayBox &diagnostic_vars_box = a_diagnostic_vars[dit()];

        for (int comp = 0; comp < NUM_DIAGNOSTIC_VARS; comp++)
        {
            diagnostic_vars_box.setVal(0.0, comp);
        }
    }
}