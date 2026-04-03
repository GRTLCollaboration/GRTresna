/* GRTresna
 * Copyright 2024 The GRTL Collaboration.
 * Please refer to LICENSE in GRTresna's root directory.
 */

#include "ScalarField.hpp"
#include "DerivativeOperators.hpp"
#include "EMTensor.hpp"
#include "FArrayBox.H"
#include "GRParmParse.hpp"
#include "Grids.hpp"
#include "IntVect.H"
#include "LevelData.H"
#include "MultigridVariables.hpp"
#include "PsiAndAijFunctions.hpp"
#include "REAL.H"
#include "RealVect.H"
#include "Tensor.hpp"

void ScalarField::initialise_matter_vars(LevelData<FArrayBox> &a_multigrid_vars,
                                         const RealVect &a_dx) const
{
    CH_assert(a_multigrid_vars.nComp() == NUM_MULTIGRID_VARS);

    DataIterator dit = a_multigrid_vars.dataIterator();
    for (dit.begin(); dit.ok(); ++dit)
    {
        // These contain the vars in the boxes, set them all to zero
        FArrayBox &multigrid_vars_box = a_multigrid_vars[dit()];

        // Iterate over the box and set non zero comps
        Box ghosted_box = multigrid_vars_box.box();
        BoxIterator bit(ghosted_box);
        for (bit.begin(); bit.ok(); ++bit)
        {

            // work out location on the grid
            IntVect iv = bit();
            RealVect loc;
            Grids::get_loc(loc, iv, a_dx, center);

            multigrid_vars_box(iv, c_phi_0) = my_phi_function(loc);
            multigrid_vars_box(iv, c_Pi_0) = my_Pi_function(loc);
        }
    }
}

// template <class data_t>
emtensor_t ScalarField::compute_emtensor(const IntVect a_iv,
                                         const RealVect &a_dx,
                                         FArrayBox &a_multigrid_vars_box) const
{
    emtensor_t out;

    DerivativeOperators derivs(a_dx);
    RealVect loc;
    Grids::get_loc(loc, a_iv, a_dx, center);

    Real psi_reg = a_multigrid_vars_box(a_iv, c_psi_reg);
    Real psi_bh = psi_and_Aij_functions->compute_bowenyork_psi(loc);
    Real psi_0 = psi_reg + psi_bh;
    Real Pi_0 = a_multigrid_vars_box(a_iv, c_Pi_0);
    Real phi_0 = a_multigrid_vars_box(a_iv, c_phi_0);

    Tensor<1, Real, SpaceDim> d1_phi;
    derivs.get_d1(d1_phi, a_iv, a_multigrid_vars_box, c_phi_0);
    Real d1_phi_squared = 0;
    FOR1(i) { d1_phi_squared += d1_phi[i] * d1_phi[i]; }

    Real V_of_phi = my_potential_function(phi_0);

    out.rho =
        0.5 * pow(psi_0, -4.0) * d1_phi_squared + 0.5 * Pi_0 * Pi_0 + V_of_phi;
    FOR1(i) { out.Si[i] = -Pi_0 * d1_phi[i]; }

    return out;
}
