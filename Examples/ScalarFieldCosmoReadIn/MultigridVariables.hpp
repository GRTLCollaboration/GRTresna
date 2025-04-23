/* GRTresna
 * Copyright 2024 The GRTL collaboration.
 * Please refer to LICENSE in GRTresna's root directory.
 */

#ifndef MULTIGRIDVARIABLES_HPP
#define MULTIGRIDVARIABLES_HPP

#include "MetricVariables.hpp"
#include "ScalarFieldVariables.hpp"

namespace MultigridVariables
{
static const std::array<std::string, NUM_METRIC_VARS> metric_variable_names =
    MetricVariables::variable_names;
static const std::array<std::string, NUM_MULTIGRID_VARS - NUM_METRIC_VARS>
    matter_variable_names = MatterVariables::variable_names;
} // namespace MultigridVariables

#endif /* MULTIGRIDVARIABLES_HPP */