/* GRTresna
 * Copyright 2024 The GRTL Collaboration.
 * Please refer to LICENSE in GRTresna's root directory.
 */

#ifndef SCALARFIELDVARIABLES_HPP_
#define SCALARFIELDVARIABLES_HPP_

#include "MetricVariables.hpp"
#include "ParityDefinitions.hpp"

// Matter Vars
enum
{
    c_phi_0 = NUM_METRIC_VARS,
    c_Pi_0,
    NUM_MULTIGRID_VARS
};

namespace MatterVariables
{

static const std::array<std::string, NUM_MULTIGRID_VARS - NUM_METRIC_VARS>
    variable_names = {"phi_0", "Pi_0"};

static constexpr std::array<int, NUM_MULTIGRID_VARS - NUM_METRIC_VARS> const
    vars_parity = {EVEN, EVEN};

} // namespace MatterVariables

#endif // SCALARFIELDVARIABLES_HPP_