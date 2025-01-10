/* GRTresna
 * Copyright 2024 The GRTL Collaboration.
 * Please refer to LICENSE in GRTresna's root directory.
 */

#ifndef BOUNDARYCONDITIONS_HPP_
#define BOUNDARYCONDITIONS_HPP_

// Chombo includes
#include "BoxIterator.H"
#include "Copier.H"
#include "FourthOrderInterpStencil.H"
#include "Interval.H"
#include "RealVect.H"

// Our includes
#include "ConstraintVariables.hpp"
#include "DimensionDefinitions.hpp"
#include "GRChomboVariables.hpp"
#include "GRParmParse.hpp"
#include "LevelData.H"
#include "MultigridVariables.hpp"
#include "ParityDefinitions.hpp"
#include "VariableType.hpp"

// Chombo namespace
#include "UsingNamespace.H"

/// Class which deals with the boundaries at the edge of the physical domain in
/// cases where they are not periodic. Currently only options are extrapolating
/// and reflective
class BoundaryConditions
{
  public:
    /// enum for possible boundary states
    enum
    {
        EXTRAPOLATING_BC,
        REFLECTIVE_BC
    };

    /// Structure containing the boundary condition params
    struct params_t
    {
        std::array<int, CH_SPACEDIM> hi_boundary;
        std::array<int, CH_SPACEDIM> lo_boundary;
        std::array<bool, CH_SPACEDIM> is_periodic;
        bool nonperiodic_boundaries_exist;
        bool reflective_boundaries_exist;
        bool extrapolating_boundaries_exist;
        bool boundary_conditions_written;
        std::array<int, NUM_METRIC_VARS> vars_parity_metric;
        std::array<int, NUM_MULTIGRID_VARS - NUM_METRIC_VARS>
            vars_parity_matter;
        std::array<int, NUM_GRCHOMBO_VARS> vars_parity_grchombo;
        std::array<int, NUM_CONSTRAINT_VARS> vars_parity_constraint;
        int extrapolation_order;
        bool Vi_extrapolated_at_boundary;
        params_t(); // sets the defaults
        void
        set_is_periodic(const std::array<bool, CH_SPACEDIM> &a_is_periodic);
        void set_hi_boundary(const std::array<int, CH_SPACEDIM> &a_hi_boundary);
        void set_lo_boundary(const std::array<int, CH_SPACEDIM> &a_lo_boundary);
        void read_params(GRParmParse &pp);
    };

  protected:
    // Member values
    double m_dx;            // The grid spacing
    int m_num_ghosts;       // the number of ghosts (usually 3)
    params_t m_params;      // the boundary params
    ProblemDomain m_domain; // the problem domain (excludes boundary cells)
    Box m_domain_box;       // The box representing the domain
    bool is_defined; // whether the BoundaryConditions class members are defined

  public:
    /// Default constructor - need to call define afterwards
    BoundaryConditions() { is_defined = false; }

    /// define function sets members and is_defined set to true
    void define(double a_dx, const params_t &a_params, ProblemDomain a_domain,
                int a_num_ghosts);

    /// write out boundary params (used during setup for debugging)
    static void write_boundary_conditions(const params_t &a_params);

    /// The function which returns the parity of each of the vars in
    int
    get_var_parity(int a_comp, int a_dir,
                   const VariableType var_type = VariableType::multigrid) const;

    /// static version used for initial output of boundary values
    static int
    get_var_parity(int a_comp, int a_dir, const params_t &a_params,
                   const VariableType var_type = VariableType::multigrid);

    /// enforce solution boundary conditions on multigrid vars
    void fill_multigrid_boundaries(const Side::LoHiSide a_side,
                                   LevelData<FArrayBox> &a_state,
                                   const Interval &a_comps,
                                   const bool filling_solver_vars);

    /// fill grchombo boundaries - used to fill output ghosts
    void fill_grchombo_boundaries(
        const Side::LoHiSide a_side, LevelData<FArrayBox> &a_state,
        const Interval &a_comps = Interval(0, NUM_GRCHOMBO_VARS - 1));

    /// fill constraint box - used to fill output ghosts
    void fill_constraint_box(
        const Side::LoHiSide a_side, FArrayBox &a_state,
        const Interval &a_comps = Interval(0, NUM_CONSTRAINT_VARS - 1));

    /// Fill the boundary values appropriately based on the params set
    /// in the direction dir
    void fill_boundary_cells_dir(const Side::LoHiSide a_side,
                                 const LevelData<FArrayBox> &a_soln,
                                 LevelData<FArrayBox> &a_out, const int dir,
                                 const int boundary_condition,
                                 const Interval &a_comps,
                                 const VariableType var_type,
                                 const bool filling_solver_vars = false);

    /// Get the boundary condition for a_dir and a_side
    int get_boundary_condition(const Side::LoHiSide a_side, const int a_dir);

    /// get the boundary box to fill if we are at a boundary
    Box get_boundary_box(const Side::LoHiSide a_side, const int a_dir,
                         const IntVect &offset_lo, const IntVect &offset_hi,
                         Box &this_ghostless_box, int shrink_for_coarse = 0);

    /// This function takes a default constructed open DisjointBoxLayout and
    /// grows the boxes lying along the boundary to include the boundaries if
    /// necessary (i.e. in the Sommerfeld BC case). It is used to define the
    /// correct DisjointBoxLayout for the exchange copier so that shared
    /// boundary ghosts are exchanged correctly.
    void expand_grids_to_boundaries(DisjointBoxLayout &a_out_grids,
                                    const DisjointBoxLayout &a_in_grids);

    friend class ExpandGridsToBoundaries;

  private:
    /// write out reflective conditions
    static void write_reflective_conditions(int idir, const params_t &a_params);

    void fill_extrapolating_cell(FArrayBox &out_box, const IntVect iv,
                                 const Side::LoHiSide a_side, const int dir,
                                 const std::vector<int> &extrapolating_comps,
                                 const int order = 1) const;

    void fill_constant_cell(FArrayBox &out_box, const IntVect iv,
                            const Side::LoHiSide a_side, const int dir,
                            const std::vector<int> &extrapolating_comps,
                            const double a_value = 0.0) const;

    void fill_reflective_cell(
        FArrayBox &out_box, const IntVect iv, const Side::LoHiSide a_side,
        const int dir, const std::vector<int> &reflective_comps,
        const VariableType var_type = VariableType::multigrid) const;
};

/// This derived class is used by expand_grids_to_boundaries to grow the
/// boxes along the Sommerfeld BC boundaries
class ExpandGridsToBoundaries : public BaseTransform
{
  public:
    ExpandGridsToBoundaries(BoundaryConditions &a_boundaries)
        : m_boundaries(a_boundaries)
    {
    }

    /// Operator called by transform to grow the boxes where required
    Box operator()(const Box &a_in_box) override;

  protected:
    BoundaryConditions &m_boundaries;
};

#endif /* BOUNDARYCONDITIONS_HPP_ */