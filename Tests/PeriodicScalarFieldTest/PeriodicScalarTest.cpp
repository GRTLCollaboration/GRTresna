/* GRTresna
 * Copyright 2024 The GRTL Collaboration.
 * Please refer to LICENSE in GRTresna's root directory.
 */

#ifdef CH_MPI
#include "mpi.h"
#endif

#include <iostream>

#include "CTTK.hpp"
#include "CTTKHybrid.hpp"
#include "Diagnostics.hpp"
#include "GRSolver.hpp"
#include "RHSTagging.hpp"
#include "ScalarField.hpp"
#include "SimulationParameters.hpp"

using namespace std;

int main(int argc, char *argv[])
{
    int failed = 0;

#ifdef CH_MPI
    MPI_Init(&argc, &argv);
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank == 0)
        cout << "Running with MPI" << endl;
#endif

    std::string in_string = argv[argc - 1];
    GRParmParse pp(0, argv + argc, NULL, in_string.c_str());

    Real Ham_norm = 0.;
    Real Mom_norm = 0.;

    SimulationParameters<CTTK<ScalarField>, ScalarField> params(pp);

    int numLevels = params.grid_params.numLevels;

    PsiAndAijFunctions *psi_and_Aij_functions =
        new PsiAndAijFunctions(params.psi_and_Aij_params);
    ScalarField *matter = new ScalarField(
        params.matter_params, psi_and_Aij_functions, params.grid_params.center,
        params.grid_params.domainLength);
    CTTK<ScalarField> *method = new CTTK<ScalarField>(
        params.method_params, matter, psi_and_Aij_functions,
        params.grid_params.numLevels, params.grid_params.center,
        params.base_params.G_Newton);
    TaggingCriterion *tagging_criterion =
        new RHSTagging<CTTK<ScalarField>, ScalarField>(
            method, matter, params.base_params.G_Newton);
    Grids *grids = new Grids(params.grid_params, tagging_criterion,
                             params.base_params.readin_matter_data);
    Diagnostics<CTTK<ScalarField>, ScalarField> *diagnostics =
        new Diagnostics<CTTK<ScalarField>, ScalarField>(
            method, matter, psi_and_Aij_functions, params.base_params.G_Newton,
            params.grid_params.center);

    Vector<LevelData<FArrayBox> *> multigrid_vars(numLevels, NULL);
    Vector<LevelData<FArrayBox> *> constraint_vars(numLevels, NULL);
    Vector<LevelData<FArrayBox> *> rhs(numLevels, NULL);
    Vector<LevelData<FArrayBox> *> diagnostic_vars(numLevels, NULL);
    Vector<RefCountedPtr<LevelData<FArrayBox>>> aCoef(numLevels);
    Vector<RefCountedPtr<LevelData<FArrayBox>>> bCoef(numLevels);

    grids->set_grids();

    IntVect ghosts = params.grid_params.num_ghosts * IntVect::Unit;
    IntVect no_ghosts = IntVect::Zero;

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

    MultilevelLinearOp<FArrayBox> mlOp;
    BiCGStabSolver<Vector<LevelData<FArrayBox> *>> solver;

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

    /* ------------- RUN SOLVER -------------- */

    bool filling_solver_vars = false;
    grids->fill_ghosts_correct_coarse(multigrid_vars, filling_solver_vars);

    openFile(params.base_params.error_filename);
    for (int NL_iter = 0; NL_iter < params.base_params.max_NL_iter; NL_iter++)
    {
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

        for (int ilev = 0; ilev < numLevels; ilev++)
        {
            RealVect dxLevel = grids->vectDx[ilev];
            diagnostics->compute_constraint_terms(multigrid_vars[ilev],
                                                  diagnostic_vars[ilev],
                                                  rhs[ilev], dxLevel);
        }

        for (int ilev = 0; ilev < numLevels; ilev++)
        {
            RealVect dxLevel = grids->vectDx[ilev];
            diagnostics->normalise_constraints(
                multigrid_vars[ilev], diagnostic_vars[ilev], rhs[ilev], dxLevel,
                params.grid_params.nCells);
        }
        Ham_norm = grids->compute_norm(diagnostic_vars, Interval(c_Ham, c_Ham));
        Mom_norm =
            grids->compute_norm(diagnostic_vars, Interval(c_Mom1, c_Mom3));

        grids->define_operator(mlOp, aCoef, bCoef, params.base_params.alpha,
                               params.base_params.beta);
        bool homogeneousBC = false;
        solver.define(&mlOp, homogeneousBC);

        solver.solve(constraint_vars, rhs);

        grids->update_psi0(multigrid_vars, constraint_vars,
                           params.method_params.deactivate_zero_mode);

        bool filling_solver_vars = true;
        grids->fill_ghosts_correct_coarse(multigrid_vars, filling_solver_vars);
    }

    if (abs(Ham_norm) > 1.e-12 || abs(Ham_norm) > 0.01)
    {
        failed = -1;
        pout() << "Tests failed, norms at the end are Ham: " << Ham_norm
               << " and Mom: " << Mom_norm << endl;
    }

    if (failed == 0)
        std::cout << "PeriodicScalar test passed..." << std::endl;
    else
        std::cout << "PeriodicScalar test failed..." << std::endl;

#ifdef CH_MPI
    MPI_Finalize();
#endif

    return failed;
}
