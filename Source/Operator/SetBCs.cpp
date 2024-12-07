/* GRTresna
 * Copyright 2024 The GRTL Collaboration.
 * Please refer to LICENSE in GRTresna's root directory.
 */

#include "SetBCs.H"
#include "AMRIO.H"
#include "BCFunc.H"
#include "BRMeshRefine.H"
#include "BiCGStabSolver.H"
#include "BoxIterator.H"
#include "CONSTANTS.H"
#include "CoarseAverage.H"
#include "ConstraintVariables.hpp"
#include "GRParmParse.hpp"
#include "LoadBalance.H"
#include "computeNorm.H"
#include "parstream.H"
#include <cmath>

// Global BCRS definitions
bool GlobalBCRS::s_areBCsParsed = false;
BoundaryConditions::params_t GlobalBCRS::s_boundary_params;

void ParseBC(FArrayBox &a_state, const Box &a_valid,
             const ProblemDomain &a_domain, Real a_dx, bool a_homogeneous)
{
    if (!a_domain.domainBox().contains(a_state.box()))
    {
        if (!GlobalBCRS::s_areBCsParsed)
        {
            GRParmParse pp;
            GlobalBCRS::s_boundary_params.read_params(pp);
            GlobalBCRS::s_areBCsParsed = true;
        }

        // fill the boundary cells and ghosts
        BoundaryConditions solver_boundaries;
        int num_ghosts = 1;
        std::array<double, 3> center = {0.0, 0.0, 0.0};
        solver_boundaries.define(a_dx, GlobalBCRS::s_boundary_params, a_domain,
                                 num_ghosts);

        // this will populate the multigrid boundaries according to the BCs
        // in particular it will fill cells for Aij, and updated K
        solver_boundaries.fill_constraint_box(
            Side::Lo, a_state, Interval(0, NUM_CONSTRAINT_VARS - 1));
        solver_boundaries.fill_constraint_box(
            Side::Hi, a_state, Interval(0, NUM_CONSTRAINT_VARS - 1));
        // Do it twice to catch corners... must be a better way
        // but for now this works
        solver_boundaries.fill_constraint_box(
            Side::Lo, a_state, Interval(0, NUM_CONSTRAINT_VARS - 1));
        solver_boundaries.fill_constraint_box(
            Side::Hi, a_state, Interval(0, NUM_CONSTRAINT_VARS - 1));
    }
}
