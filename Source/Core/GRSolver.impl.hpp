/* GRTresna
 * Copyright 2024 The GRTL Collaboration.
 * Please refer to LICENSE in GRTresna's root directory.
 */

#ifndef GRSOLVER_HPP_
#error "This file should only be included through GRSolver.hpp"
#endif

#include "Diagnostics.hpp"
#include "GRParmParse.hpp"
#include "GRSolver.hpp"
#include "PsiAndAijFunctions.hpp"
#include "RHSTagging.hpp"
#include "WriteFile.hpp"
#include "WriteOutput.H"

template <class method_t, class matter_t>
GRSolver<method_t, matter_t>::GRSolver(GRParmParse &a_pp)
    : pp(a_pp), params(pp), numLevels(params.grid_params.numLevels),
      multigrid_vars(numLevels, NULL), constraint_vars(numLevels, NULL),
      rhs(numLevels, NULL), aCoef(numLevels), bCoef(numLevels),
      diagnostic_vars(numLevels, NULL), Ham_error(0.), Mom_error(0.)
{
    psi_and_Aij_functions = new PsiAndAijFunctions(params.psi_and_Aij_params);
    matter = new matter_t(params.matter_params, psi_and_Aij_functions,
                          params.grid_params.center,
                          params.grid_params.domainLength);
    method =
        new method_t(params.method_params, matter, psi_and_Aij_functions,
                     params.grid_params.numLevels, params.grid_params.center,
                     params.base_params.G_Newton);
    tagging_criterion = new RHSTagging<method_t, matter_t>(
        method, matter, params.base_params.G_Newton);
    grids = new Grids(params.grid_params, tagging_criterion,
                      params.base_params.readin_matter_data);
    diagnostics = new Diagnostics<method_t, matter_t>(
        method, matter, psi_and_Aij_functions, params.base_params.G_Newton,
        params.grid_params.center);
}

template <class method_t, class matter_t>
void GRSolver<method_t, matter_t>::setup()
{
    // set up the grids, using the rhs for tagging to decide
    // where needs additional levels
    if (params.base_params.readin_matter_data)
    {
        grids->read_grids(params.base_params.input_filename);
    }
    else
    {
        grids->set_grids();
    }

    create_vars();

    mlOp.m_num_mg_iterations = params.base_params.numMGIter;
    mlOp.m_num_mg_smooth = params.base_params.numMGSmooth;
    mlOp.m_preCondSolverDepth = params.base_params.preCondSolverDepth;

    // define the multi level operator
    grids->define_operator(mlOp, aCoef, bCoef, params.base_params.alpha,
                           params.base_params.beta);

    // set the solver params
    bool homogeneousBC = false;
    solver.define(&mlOp, homogeneousBC);
    solver.m_verbosity = params.base_params.verbosity;
    solver.m_normType = 2;
    solver.m_eps = params.base_params.iter_tolerance;
    solver.m_imax = params.base_params.max_iter;

    output_solver_data(constraint_vars, multigrid_vars, diagnostic_vars,
                       grids->grids_data, params, 0);
}

template <class method_t, class matter_t>
int GRSolver<method_t, matter_t>::run()
{
    // Iterate linearised Poisson eqn for NL solution

    bool filling_solver_vars = false;
    grids->fill_ghosts_correct_coarse(multigrid_vars, filling_solver_vars);

    openFile(params.base_params.error_filename);
    for (int NL_iter = 0; NL_iter < params.base_params.max_NL_iter; NL_iter++)
    {
        pout() << "Main Loop Iteration " << (NL_iter + 1) << " out of "
               << params.base_params.max_NL_iter << endl;

        for (int ilev = 0; ilev < numLevels; ilev++)
        {
            RealVect dxLevel = grids->vectDx[ilev];
            method->solve_analytic(multigrid_vars[ilev], rhs[ilev],
                                   grids->vectDx[ilev]);
        }
        filling_solver_vars = false;
        grids->fill_ghosts_correct_coarse(multigrid_vars, filling_solver_vars);

        for (int ilev = 0; ilev < numLevels; ilev++)
        {
            RealVect dxLevel = grids->vectDx[ilev];
            method->set_elliptic_terms(multigrid_vars[ilev], rhs[ilev],
                                       aCoef[ilev], bCoef[ilev],
                                       grids->vectDx[ilev]);
        }

        calculate_diagnostics(NL_iter);

        grids->define_operator(mlOp, aCoef, bCoef, params.base_params.alpha,
                               params.base_params.beta);
        bool homogeneousBC = false;
        solver.define(&mlOp, homogeneousBC);

        solver.solve(constraint_vars, rhs);

        grids->update_psi0(multigrid_vars, constraint_vars,
                           params.method_params.deactivate_zero_mode);

        bool filling_solver_vars = true;
        grids->fill_ghosts_correct_coarse(multigrid_vars, filling_solver_vars);

        // Only write out at requested intervals
        bool at_diagnostic_interval =
            ((NL_iter + 1) % params.base_params.diagnostic_interval == 0);
        if (params.base_params.write_diagnostics && at_diagnostic_interval)
        {
            output_solver_data(constraint_vars, multigrid_vars, diagnostic_vars,
                               grids->grids_data, params, NL_iter + 1);
        }
    }

    pout() << "Converged!" << endl
           << "Ham relative error: " << Ham_error << " %" << endl
           << "Mom relative error: " << Mom_error << " %" << endl;

    // Mayday if result not converged (> 100% error)
    if (Ham_error > 1e2 || Mom_error > 1e2)
    {
        MayDay::Error(
            "NL iterations did not converge - may need a better initial guess");
    }

    output_final_data(multigrid_vars, grids->grids_data, grids->vectDx,
                      grids->vectDomain, params,
                      params.base_params.output_filename);

    int exitStatus = solver.m_exitStatus;
    // note that for AMRMultiGrid, success = 1.
    exitStatus -= 1;
    return exitStatus;
}

template <typename method_t, typename matter_t>
void GRSolver<method_t, matter_t>::calculate_diagnostics(const int NL_iter)
{
    if (params.grid_params.periodic_directions_exist)
    {
        // Calculate values for integrand here
        pout() << "Computing integrability of rhs for periodic domain... "
               << endl;
        for (int icomp = 0; icomp < NUM_CONSTRAINT_VARS; icomp++)
        {
            Real integral = grids->compute_sum(rhs, Interval(icomp, icomp));

            pout() << "Integral of rhs " << icomp << " is " << integral << endl;
            pout() << "(This should be a small number)" << endl;
        }
    }

    for (int ilev = 0; ilev < numLevels; ilev++)
    {
        RealVect dxLevel = grids->vectDx[ilev];
        diagnostics->compute_constraint_terms(
            multigrid_vars[ilev], diagnostic_vars[ilev], rhs[ilev], dxLevel);
    }

    for (int ilev = 0; ilev < numLevels; ilev++)
    {
        RealVect dxLevel = grids->vectDx[ilev];
        diagnostics->normalise_constraints(multigrid_vars[ilev],
                                           diagnostic_vars[ilev], rhs[ilev],
                                           dxLevel, params.grid_params.nCells);
    }
    Real Ham_norm =
        grids->compute_norm(diagnostic_vars, Interval(c_Ham, c_Ham));
    Real Mom_norm =
        grids->compute_norm(diagnostic_vars, Interval(c_Mom, c_Mom));
    Real Ham_abs_norm =
        grids->compute_norm(diagnostic_vars, Interval(c_Ham_abs, c_Ham_abs));
    Real Mom_abs_norm =
        grids->compute_norm(diagnostic_vars, Interval(c_Mom_abs, c_Mom_abs));

    Ham_error = 100 * Ham_norm / Ham_abs_norm;
    Mom_error = 100 * Mom_norm / Mom_abs_norm;

    pout() << "The relative error of Ham before step " << NL_iter << " is "
           << Ham_error << " %" << endl;
    pout() << "The relative error of Mom before step " << NL_iter << " is "
           << Mom_error << " %" << endl;
    writeFile(params.base_params.error_filename, NL_iter, Ham_error, Mom_error);
}

template <typename method_t, typename matter_t>
void GRSolver<method_t, matter_t>::create_vars()
{
    IntVect ghosts = params.grid_params.num_ghosts * IntVect::Unit;
    IntVect no_ghosts = IntVect::Zero;

    // Initialise core vars
    for (int ilev = 0; ilev < numLevels; ilev++)
    {
        multigrid_vars[ilev] = new LevelData<FArrayBox>(
            grids->grids_data[ilev], NUM_MULTIGRID_VARS, ghosts);

        constraint_vars[ilev] = new LevelData<FArrayBox>(
            grids->grids_data[ilev], NUM_CONSTRAINT_VARS, ghosts);
        rhs[ilev] = new LevelData<FArrayBox>(grids->grids_data[ilev],
                                             NUM_CONSTRAINT_VARS, no_ghosts);
        aCoef[ilev] =
            RefCountedPtr<LevelData<FArrayBox>>(new LevelData<FArrayBox>(
                grids->grids_data[ilev], NUM_CONSTRAINT_VARS, no_ghosts));
        bCoef[ilev] =
            RefCountedPtr<LevelData<FArrayBox>>(new LevelData<FArrayBox>(
                grids->grids_data[ilev], NUM_CONSTRAINT_VARS, no_ghosts));
        diagnostic_vars[ilev] = new LevelData<FArrayBox>(
            grids->grids_data[ilev], NUM_DIAGNOSTIC_VARS, no_ghosts);
    }

    for (int ilev = 0; ilev < numLevels; ilev++)
    {
        RealVect dxLevel = grids->vectDx[ilev];

        method->initialise_method_vars(*multigrid_vars[ilev], dxLevel);
        matter->initialise_matter_vars(*multigrid_vars[ilev], dxLevel);

        method->initialise_constraint_vars(*constraint_vars[ilev], dxLevel);
        method->initialise_constraint_vars(*rhs[ilev], dxLevel);
        method->initialise_constraint_vars(*aCoef[ilev], dxLevel);
        method->initialise_constraint_vars(*bCoef[ilev], dxLevel);

        diagnostics->initialise_diagnostic_vars(*diagnostic_vars[ilev],
                                                dxLevel);
    }
}

template <class method_t, class matter_t>
GRSolver<method_t, matter_t>::~GRSolver()
{

    delete grids;
    delete psi_and_Aij_functions;
    delete diagnostics;

    for (int ilev = 0; ilev < numLevels; ilev++)
    {
        delete multigrid_vars[ilev];
        delete constraint_vars[ilev];
        delete rhs[ilev];
        delete diagnostic_vars[ilev];
    }
}
