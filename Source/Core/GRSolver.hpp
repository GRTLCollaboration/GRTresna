/* GRTresna
 * Copyright 2024 The GRTL Collaboration.
 * Please refer to LICENSE in GRTresna's root directory.
 */

#ifndef GRSOLVER_HPP_
#define GRSOLVER_HPP_

#include "Diagnostics.hpp"
#include "GRParmParse.hpp"
#include "PsiAndAijFunctions.hpp"
#include "SimulationParameters.hpp"
#include "TaggingCriterion.hpp"

/*
Class that manages high-level solver functionality, independent of specific
method
*/
template <class method_t, class matter_t, class tagging_t> class GRSolver
{
  public:
    GRSolver(GRParmParse &pp);

    void setup();

    int run();

    ~GRSolver();

  private:
    void create_vars();

    void calculate_diagnostics(const int NL_iter);

    Real Ham_error;
    Real Mom_error;

    GRParmParse pp;
    SimulationParameters<method_t, matter_t> params;

    int numLevels; // Must be after params

    method_t *method;
    matter_t *matter;

    Diagnostics<method_t, matter_t> *diagnostics;

    PsiAndAijFunctions *psi_and_Aij_functions;

    Grids *grids;
    TaggingCriterion *tagging_criterion;

    MultilevelLinearOp<FArrayBox> mlOp;
    BiCGStabSolver<Vector<LevelData<FArrayBox> *>> solver;

    Vector<LevelData<FArrayBox> *> multigrid_vars;

    Vector<LevelData<FArrayBox> *> constraint_vars;
    Vector<LevelData<FArrayBox> *> rhs;
    Vector<LevelData<FArrayBox> *> diagnostic_vars;
    Vector<RefCountedPtr<LevelData<FArrayBox>>> aCoef;
    Vector<RefCountedPtr<LevelData<FArrayBox>>> bCoef;
};

#include "GRSolver.impl.hpp"

#endif /* GRSOLVER_HPP_ */