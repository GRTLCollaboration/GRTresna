/* GRTresna
 * Copyright 2024 The GRTL Collaboration.
 * Please refer to LICENSE in GRTresna's root directory.
 */

#ifndef MATTERPARAMS_HPP_
#define MATTERPARAMS_HPP_

#include "GRParmParse.hpp"
#include "REAL.H"
#include <fstream>
#include <iostream>
#include <string>

namespace MatterParams
{

struct params_t
{
    Real phi_0;
    Real dphi;
    Real pi_0;
    Real dpi;
    Real scalar_mass;

    // params for read in
    std::string read_from_data_dphi;
    std::string read_from_data_dpi;
    Real *input_dphi;
    Real *input_dpi;
    int lines;
    Real spacing;
};

inline void read_params(GRParmParse &pp, params_t &matter_params)
{
    pp.get("phi_0", matter_params.phi_0);
    pp.get("dphi", matter_params.dphi);
    pp.get("pi_0", matter_params.pi_0);
    pp.get("dpi", matter_params.dpi);
    pp.get("scalar_mass", matter_params.scalar_mass);

    // read in of data
    pp.get("read_from_data_dphi", matter_params.read_from_data_dphi);
    pp.get("read_from_data_dpi", matter_params.read_from_data_dpi);
    pp.get("data_lines", matter_params.lines);
    pp.get("data_spacing", matter_params.spacing);
#define num 256
    static Real inputdata_dphi[num * num * num];
    static Real inputdata_dpi[num * num * num];
    Real tmp_dphi = 0.0;
    Real tmp_dpi = 0.0;
    ifstream ifspsi_dphi(matter_params.read_from_data_dphi);
    ifstream ifspsi_dpi(matter_params.read_from_data_dpi);
    for (int i = 0;
         i < matter_params.lines * matter_params.lines * matter_params.lines;
         ++i)
    {
        ifspsi_dphi >> tmp_dphi;
        ifspsi_dpi >> tmp_dpi;
        inputdata_dphi[i] = tmp_dphi;
        inputdata_dpi[i] = tmp_dpi;
    }
    matter_params.input_dphi = inputdata_dphi;
    matter_params.input_dpi = inputdata_dpi;
}

}; // namespace MatterParams

#endif