/* GRTresna
 * Copyright 2024 The GRTL Collaboration.
 * Please refer to LICENSE in GRTresna's root directory.
 */

#include "mpi.h"
#include <iostream>

#include "CTTK.hpp"
#include "CTTKHybrid.hpp"
#include "GRParmParse.hpp"
#include "GRSolver.hpp"
#include "ScalarField.hpp"

using namespace std;

int main(int argc, char *argv[])
{
    int status = 0;
#ifdef _OPENMP
    std::cout << "#threads = " << omp_get_max_threads() << std::endl;
#endif
#ifdef CH_MPI
    MPI_Init(&argc, &argv);
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank == 0)
        cout << "Running with MPI" << endl;
#endif

    if (argc < 1)
    {
        cerr << " usage " << argv[0] << " <input_file_name> " << endl;
        exit(0);
    }

    // Read params input file
    char *inFile = argv[1];
    GRParmParse pp(argc - 2, argv + 2, NULL, inFile);

    // Create solver, specifying method and matter
    // Need to define these in the environment - see
    // https://github.com/GRTLCollaboration/GRTresna/wiki/Solver-methods
    // for more details
#if defined(USE_CTTK)
    GRSolver<CTTK<ScalarField>, ScalarField> solver(pp);
    pout() << "Using CTTK" << endl;
#elif defined(USE_CTTKHybrid)
    GRSolver<CTTKHybrid<ScalarField>, ScalarField> solver(pp);
    pout() << "Using CTTK Hybrid" << endl;
#else
#error                                                                         \
    "No valid GRSolver method defined. Please define one of the methods in the GNUMakefile. See the wiki page on solver methods for more details."
#endif

    solver.setup();
    status = solver.run();

#ifdef CH_MPI
    MPI_Finalize();
#endif
    return status;
}