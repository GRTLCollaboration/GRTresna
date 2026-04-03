/* GRTresna
 * Copyright 2024 The GRTL Collaboration.
 * Please refer to LICENSE in GRTresna's root directory.
 */

#ifndef SCALARFIELD_HPP_
#define SCALARFIELD_HPP_

#include "EMTensor.hpp"
#include "FArrayBox.H"
#include "GRParmParse.hpp"
#include "IntVect.H"
#include "LevelData.H"
#include "MatterParams.hpp"
#include "PsiAndAijFunctions.hpp"
#include "REAL.H"
#include "RealVect.H"
#include "Tensor.hpp"

class ScalarField
{
  public:
    using params_t = MatterParams::params_t;

    ScalarField(params_t a_matter_params,
                PsiAndAijFunctions *a_psi_and_Aij_functions,
                const std::array<double, SpaceDim> a_center,
                RealVect a_domainLength)
        : m_matter_params(a_matter_params),
          psi_and_Aij_functions(a_psi_and_Aij_functions), center(a_center),
          domainLength(a_domainLength)
    {
    }

    //! The function which calculates the EM Tensor, given the vars and
    //! derivatives, including the potential
    // template <class data_t>
    emtensor_t compute_emtensor(const IntVect a_iv, const RealVect &a_dx,
                                FArrayBox &a_multigrid_vars_box) const;

    static void read_params(GRParmParse &pp, params_t &matter_params)
    {
        MatterParams::read_params(pp, matter_params);
    }

    void initialise_matter_vars(LevelData<FArrayBox> &a_multigrid_vars,
                                const RealVect &a_dx) const;

    Real my_potential_function(const Real &phi_here) const;

    Real my_phi_function(const RealVect &locr) const;

    Real my_Pi_function(const RealVect &loc) const;

    params_t m_matter_params;

    PsiAndAijFunctions *psi_and_Aij_functions;

    ~ScalarField() {}

  private:
    const std::array<double, SpaceDim> center;
    RealVect domainLength;
};

#endif /* SCALARFIELD_HPP_ */
