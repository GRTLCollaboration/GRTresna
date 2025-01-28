---
title: 'GRTresna: An open-source code to solve the initial data constraints in numerical relativity'
tags:
  - C++
  - MPI
  - Open MP
  - gravity
  - general relativity
  - numerical relativity
  - initial data
authors:
- name: Josu C. Aurrekoetxea
  orcid: 0000-0001-9584-5791
  affiliation: 1
- name: Sam E. Brady
  orcid: 0009-0000-5568-839X
  affiliation: 2
- name: Llibert Areste Salo
  orcid: 0000-0002-3812-8523
  affiliation: 3
- name: Jamie Bamber
  orcid: 0000-0001-7181-3365
  affiliation: 4
- name: Liina Chung-Jukko
  orcid: 0000-0002-2559-7734
  affiliation: 5
- name: Katy Clough
  orcid: 0000-0001-8841-1522
  affiliation: 2
- name: Eloy de Jong
  orcid: 0000-0002-4505-0808
  affiliation: 5
- name: Matthew Elley
  orcid: 0000-0003-3167-5601
  affiliation: 6
- name: Pau Figueras
  orcid: 0000-0001-6438-315X
  affiliation: 2
- name: Thomas Helfer
  orcid: 0000-0001-6880-1005
  affiliation: 7
- name: Eugene A. Lim
  orcid: 0000-0002-6227-9540
  affiliation: 5
- name: Miren Radia
  orcid: 0000-0001-8861-2025
  affiliation: 8
- name: Areef Waeming
  orcid: 0009-0001-6700-4116
  affiliation: 2
- name: Zipeng Wang
  orcid: 0000-0002-4745-8209
  affiliation: 9
affiliations:
- name: Center for Theoretical Physics, Massachusetts Institute of Technology, 77 Massachusetts Ave, Cambridge, MA 02139, USA
  index: 1
- name: Geometry, Analysis and Gravitation, School of Mathematical Sciences, Queen Mary University of London, Mile End Road, London E1 4NS, United Kingdom
  index: 2
- name: Institute for Theoretical Physics, KU Leuven, Celestijnenlaan 200D, B-3001 Leuven, Belgium
  index: 3  
- name: Department of Physics, University of Illinois Urbana-Champaign, Urbana, IL 61801, USA
  index: 4
- name: Theoretical Particle Physics and Cosmology, King's College London, Strand, London, WC2R 2LS, United Kingdom
  index: 5
- name: Department of Physics, University of Basque Country, UPV/EHU, 48080, Bilbao, Spain
  index: 6
- name: Institute for Advanced Computational Science, Stony Brook University, NY 11794 USA
  index: 7
- name: Information Services, University of Cambridge, Roger Needham Building, 7 JJ Thomson Avenue, Cambridge, CB3 0WA, United Kingdom
  index: 8
- name: Institute for Advanced Computational Science, Stony Brook University, NY 11794 USA
  index: 9
- name: Department of Physics and Astronomy, Johns Hopkins University, Baltimore, MD 21218, USA
  index: 10
date: 22 Jan 2025
bibliography: paper.bib

---

GRTresna is a multigrid solver designed to solve the constraint equations for the initial data required in numerical relativity simulations. In particular it is focussed on scenarios with fundamental fields around black holes and inhomogeneous cosmological spacetimes. The following overview has been prepared as part of the submission of the code to the Journal of Open Source Software. The code is based on the formalism in Aurrekoetxea, Clough \& Lim [@Aurrekoetxea:2022mpw] and can be found at https://github.com/GRTLCollaboration/GRTresna

# Summary

Numerical relativity (NR) is a tool for the solution of the Einstein Equations, which describe gravity in strong field regimes. The equations can be expressed as a set of coupled partial differential equations (PDEs) for the 10 metric quantities $g_{\mu\nu}$ and their time derivatives $\partial_t g_{\mu\nu}$. NR is primarily focussed on the hyperbolic PDEs that describe their time evolution from an initial data set, but the initial data itself must satisfy a set of four coupled non-linear elliptic PDEs known as the Hamiltonian and momentum constraints. Whilst these constraints can be solved more straightforwardly by making certain assumptions, this significantly restricts the range of physical scenarios that can be studied. A general solver therefore expands the physics that NR evolutions can be used to probe. 

In the ADM form of the Einstein Equations [@Arnowitt:1962hi], we slice the spacetime into 3-dimensional hypersurfaces

$$ds^2 = -(\alpha^2 - \beta_i\beta^i) dt^2 + 2\beta_i dx^i dt + \gamma_{ij} dx^i dx^j$$

and the elliptic constraints are expressed as

$$H \equiv R + K^2-K_{ij}K^{ij}-16\pi \rho = 0,$$

$$M_i \equiv D^j (K_{ij}- \gamma_{ij} K) - 8\pi S_i = 0.$$

Here, $\gamma_{ij}$ is the 3-metric of the hypersurface, $R$ is the Ricci scalar associated to this metric, and $K_{ij}\sim \partial_t \gamma_{ij}$ is the extrinsic curvature tensor, with $K=\gamma^{ij}K_{ij}$ its trace. The decomposed components of the stress-energy tensor of matter (measured by normal observers) are defined as $\rho = n_\mu n_\nu T^{\mu\nu}$ and $S_i = -\gamma_{i\mu} n_\nu T^{\mu\nu}$, where $n_\mu = (-\alpha,0,0,0)$. These equations constitute the set of four PDEs to be solved. There are 16 unknowns: 6 in $\gamma_{ij}$, 6 in $K_{ij}$, 1 in $\rho$ and 3 in $S_i$. Usually the matter configuration is set by the physical scenario, which determines $\rho$ and $S_i$. The constraints only determine 4 quantities, and the remaining 8 (4 of which are physical degrees of freedom, and 4 gauge choices) must be chosen according to physical principles or knowledge about the system. In this short paper, we introduce GRTresna, an open-source code to solve these equations.


The two main methods for finding initial conditions in numerical relativity are the *conformal traverse-traceless* (CTT) and the *conformal thin sandwich* (CTS) approaches.  We refer the reader to the standard NR texts [@Alcubierre:2008co;@Gourgoulhon:2007ue;@Baumgarte:2010ndz;@Baumgarte:2021skc;@Shibata_book] for more details about these. GRTresna implements two variations of the CTT method recently introduced in Aurrekoetxea, Clough \& Lim [@Aurrekoetxea:2022mpw]: the CTTK and CTTK-Hybrid methods, which are particularly well-suited to cases with fundamental fields in the matter content. Documentation about using and modifying GRTresna can be found in the code wiki https://github.com/GRTLCollaboration/GRTresna/wiki.

![plot_grtresna](plot.png)
*Some highlights of work using GRTresna to date: (Left:) Dark matter around binary black holes, from [@Bamber:2022pbs;@Aurrekoetxea:2023jwk;@Aurrekoetxea:2024cqd] (Middle:) Evolution of inflationary perturbations during preheating, from [@Aurrekoetxea:2023jwd]. (Right:) Scalar fields around black holes in* $4\partial ST$ *gravity, from [@Brady:2023dgu]*

# Key features of GRTresna

The key features of GRTresna are as follows

- Flexibility: GRTresna is designed to be extended to various physical scenarios, including different matter types and gravitational theories beyond GR. It currently supports cosmological-type periodic spacetimes and a superposition of two boosted and/or spinning black holes (Bowen-York initial data), with fully general scalar field matter source configurations and the flexibility to adapt to other setups. While scalar fields are the only matter sources included in the current version of the code, the templated methods allow users to easily replace them with other matter types by copying the scalar field implementation and modifying the methods to compute the corresponding energy and momentum densities.
    
- Methods: GRTresna incorporates the CTTK and CTTK-Hybrid methods to solve the Hamiltonian and momentum constraints. These methods offer several advantages when dealing with fundamental fields, as discussed in [@Aurrekoetxea:2022mpw]. The method code is also templated, so users can easily implement their preferred methods.

- Initial conditions: The code supports analytical initial data for the matter fields, as well as the option to read grids and data from an existing HDF5 file. This functionality is especially useful when combined with our code that evolves matter on fixed metric backgrounds GRDzhadzha [@Aurrekoetxea:2023fhl], meaning that we can upgrade the resulting matter configurations to full NR simulations with backreaction.
    
- Boundary conditions: The code implements extrapolating, reflective, and periodic boundary conditions, compatible with those in the NR evolution code GRChombo [@Andrade:2021rbd;@Clough:2015sqa].

- Diagnostics: The code computes the Hamiltonian and momentum constraint errors at each iteration step, and outputs the norm of these values across the grid to a text file.

- Compatibility: As GRTresna is developed on top of Chombo, the solver is primarily designed to be compatible with GRChombo [@Andrade:2021rbd;@Clough:2015sqa] and the family of codes developed by the GRTL Collaboration. We provide two examples that integrate directly with existing examples in the GRChombo evolution code (via the output of a checkpoint file for restart at $t=0$), and provide guidance and tools to validate the results. However, the code outputs data in the standard HDF5 data format, which should be straightforward to adapt to other NR codes that support HDF5 input or can be accessed using Python.

Other features that are inherited from Chombo include

- C++ class structure: GRTresna is written in the C++ language, and makes heavy use of object-oriented programming (OOP) and templating.
    
- Parallelism: GRTresna uses hybrid OpenMP/MPI parallelism.
    
- Adaptive Mesh Refinement: The code inherits the flexible AMR grid structure of Chombo, with block-structured Berger-Rigoutsos grid generation [@Berger:1991]. The tagging of refinement regions is fully flexible and while it is based on the sources of the elliptic equations by default, other user-defined measures can be defined [@Radia:2021smk].

- Fast: The code uses a multigrid method to efficiently reduce errors across a hierarchy of discretizations, enabling the solver to achieve rapid convergence while minimizing computational costs. This makes GRTresna highly optimized for handling the demanding computations of initial data in the presence of AMR.

Forthcoming features currently under development include the addition of other methods, in particular the Extended Conformal Thin Sandwich (XCTS) method, non-conformally flat metric data, new matter types including vector fields, the modified scalar-tensor gravity formalism of Brady et. al. [@Brady:2023dgu] and dimensional reduction to 2D using the cartoon formalism [@Alcubierre:1999ab;@Cook:2016soy].


# Statement of Need


There are a number of existing initial data solvers for numerical relativity, most of which are primarily designed to solve for initial conditions in compact object mergers (i.e. neutron stars and black holes). These include TwoPunctures [@Ansorg:2004ds], SGRID [@Tichy:2009yr], BAM [@Bruegmann:2006ulg], LORENE [@LORENE;@Gourgoulhon:2000nn], Spells [@Pfeiffer:2002wt], SpECTRE [@Vu:2021coj;@Vu:2024cgf] ([@Nee:2024bur] for modified gravity) COCAL [@Uryu:2011ky;@Tsokaros:2012kp;@Tsokaros:2015fea], PCOCAL [@Boukas:2023ckb], Elliptica [@Rashti:2021ihv], NRPyElliptic [@Assumpcao:2021fhq], KADATH/FUKA [@FUKA;@Grandclement:2009ju;@Papenfort:2021hod}], SPHINCS\_ID [@SPHINCSID;@Diener:2022hui;@Rosswog:2023nnl], and the solver of East *et al.* [@East:2012zn]. Many of these codes, particularly those using spectral methods like TwoPunctures and SpECTRE, provide a higher accuracy in the solution compared to GRTresna, which is limited to second order accuracy by the multigrid method used. They are therefore better suited to initial data for waveform generation where precision is key. GRTresna is, however, designed to be more flexible and general purpose, tackling both cosmological and black hole spacetimes in a range of scenarios beyond GR and the Standard Model.

In particular, to the best of our knowledge, there is no fully general, publicly available initial condition solver for inhomogeneous cosmological spacetimes. One exception is FLRWSolver, developed by Macpherson \emph{et al.} [@Macpherson:2016ict] as part of the Einstein Toolkit [@Loffler:2011ay], which specializes in initializing data for cosmological perturbations arising from inflation for studies of late-time cosmology. However, this is limited to only weakly non-linear initial data. GRTresna aims to provide an open-source tool that not only incorporates the general features of existing initial data solvers for compact objects in GR but also extends their capabilities to cosmological spacetimes (see [@Aurrekoetxea:2024mdy] for a review of the application of numerical relativity in cosmology). GRTresna is particularly well-suited for fundamental field matter types, such as scalar and vector fields. Its flexible design allows users to implement new solver methods, additional matter types, or extend the code to study theories beyond GR. It is fully compatible with the GRTL Collaboration's ecosystem of codes but can also serve as a complementary tool for generating constraint-satisfying initial data for other numerical relativity codes.

# Key research projects using GRChombo

The code has already been used successfully to study a range of problems in fundamental physics, including:
- The robustness of inflation to inhomogeneities in the scalar field [@Aurrekoetxea:2019fhr;@Elley:2024alx].
- The formation of oscillons during inflationary preheating [@Aurrekoetxea:2023jwd].
- Formation of spinning primordial black holes [@deJong:2023gsx].
- The effect of scalar dark matter environments around binary black holes [@Bamber:2022pbs;@Aurrekoetxea:2023jwk;@Aurrekoetxea:2024cqd].
- The general relativistic evolution of polarized Proca stars [@Wang:2023tly].
- Solving the initial conditions problem for modified gravity theories [@Brady:2023dgu].


# Acknowledgements

We thank the GRTL collaboration (www.grtlcollaboration.org) for their support and code development work.
JCA acknowledges funding from the Department of Physics at MIT. SEB is supported by a QMUL Principal studentship. JB acknowledges funding by the National Science Foundation (NSF) Grants PHY-2308242, OAC-2310548 and PHY-2006066 to the
University of Illinois at Urbana-Champaign. KC acknowledges funding from the UKRI Ernest Rutherford Fellowship (grant number ST/V003240/1). KC and PF acknowledge funding from STFC Research Grant ST/X000931/1 (Astronomy at Queen Mary 2023-2026). ME has been supported in part by the PID2021-123703NB-C21 grant funded by MCIN/AEI/10.13039/501100011033/and by ERDF: ``A way of making Europe''; the Basque Government grant (IT-1628-22). EAL acknowledges support from a Leverhulme Trust Research Project Grant. AW acknowledges the support of the Development and Promotion of Science and Technology Talents Project (DPST), the Institute for the Promotion of Teaching Science and Technology (IPST), Thailand. ZW is supported by NSF Grants No. AST-2307146, PHY-2207502, PHY-090003 and PHY-20043, by NASA Grant No. 21-ATP21-0010, by the John Templeton Foundation Grant 62840, by the Simons Foundation, and by the Italian Ministry of Foreign Affairs and International Cooperation grant No. PGR01167.

Development of this code used the DiRAC Memory Intensive services Cosma8 and Cosma7 at Durham University, managed by the Institute for Computational Cosmology on behalf of the STFC DiRAC HPC Facility (www.dirac.ac.uk). The DiRAC service at Durham was funded by BEIS, UKRI and STFC capital funding, Durham University and STFC operations grants. DiRAC is part of the UKRI Digital Research Infrastructure.

# References

