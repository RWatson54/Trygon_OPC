/*
!********************************************************************!
!********************************************************************!
!                                                                    !
!   Trygon -- Main program for fully featured FR solver              !
!                                                                    !
!********************************************************************!
!                                                                    !
!   Version history:                                                 !
!                    Program created: 12Apr23                - raw54 !
!                                                                    !
!********************************************************************!
!                                                                    !
!   Known issues:                                                    !
!                    No known issues.                        - raw54 !
!                                                                    !
!********************************************************************!
!********************************************************************!
! 
!   Doxygen section: 
! 
!>  @author 
!>  Rob Watson 
! 
!>  @brief 
!>  Trygon is a flux reconstruction based solver which is designed to
!>  make use of the OP2 unstructured parallelisation libraries.
!>
!>  Key dependencies are (to be installed in this order):
!>    -- (Parallel-enabled) HDF5
!>
!>    -- ParMetis 
!>       -- Because of classic stupidity, GKlib isn't added to OP2 search path
!>       -- Install GKlib, Metis, ParMetis all into the same prefix
!>       -- Then merge GKlib and Metis libraries into just Metis:
!>          -- ar -x libGKlib.a
!>          -- ar -x libmetis.a
!>          -- ar -qc libmetis.a *.o
!>          -- rm *.o; rm libGKlib.a
!>    --     OR
!>    -- Scotch and PT-Scotch
!>       -- bison and flex
!>       -- (Change the Make.inc makefile with -DIDXSIZE32 and without DSCOTCH_PTHREAD)
!>       -- make -j scotch
!>       -- (Change the Make.inc makefile with CCD -> CCD = mpicc)
!>       -- make pt scotch
!>    --     OR
!>    -- KaHIP
!>       -- Not tried this yet, so no instructions forthcoming
!>
!>    -- HDF5 (with and without --enable-parallel)
!>       -- Try and compile only the static versions - or delete all the *.so* files
!>       -- CC=gcc FC=gfortran ./configure --enable-fortran --prefix=/home/raw54/Source/hdf5-1.14.0/hdf5_seq
!>       -- make
!>       -- make check
!>       -- make install
!>       -- make check-install
!>       -- CC=mpicc FC=mpifort ./configure --enable-fortran --enable-parallel --prefix=/home/raw54/Source/hdf5-1.14.0/hdf5_par
!>       -- make
!>       -- make check
!>       -- make install
!>       -- make check-install
!>
!>    -- OP2-DSL
!>       -- export OP2_COMPILER=gnu
!>       -- export PTSCOTCH_INSTALL_PATH=/home/raw54/Source/Scotch
!>       -- export HDF5_SEQ_INSTALL_PATH=/home/raw54/Source/hdf5-1.14.0/hdf5_seq
!>       -- export HDF5_PAR_INSTALL_PATH=/home/raw54/Source/hdf5-1.14.0/hdf5_par
!>       -- make -C op2 config
!>       -- If using the GCC compiler, this will probably cause an error related to a bug in the KaHIP extension
!>          -- either install KaHIP, or take the real_t declaration out of the kahip ifdef 
!>             in op2/src/mpi/op_mpi_part_core.cpp, or switch to intel compiler (?)
!>       -- make -C op2
! 
!********************************************************************!
!********************************************************************!
*/
/*
// standard headers
*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*
// global constants
*/

double gam, gm1, cfl, eps, mach, alpha, qinf[4];

/*
// OP header file
*/

#include "op_seq.h"

/*
// kernel routines for parallel loops
*/

#include "qts2q0s.h"
#include "qts2qtf.h"
#include "qts2fhs.h"
#include "fhs2dts.h"
#include "qtf2itf.h"
#include "fhs2ftf.h"
#include "dts2cts.h"
#include "cts2qts.h"

/*
// include some other header files
*/
#include "welcome.h"

/*
// main program
*/

int main(int argc, char **argv) {

  /* 
  // Define a couple of compiler-level constants for dimensioning
  */ 

  #define nDim 2
  #define nPde 1

  /*
  // OP initialisation
  */

  op_init(argc, argv, 2);

  /*
  // Say hello
  */

  WL_hello();

  // timer and counter
  int niter;
  double cpu_t1, cpu_t2, wall_t1, wall_t2;

  /* 
  // Set the filename for the grid
  */

  char file[] = "new_grid.h5";

  /*
  // Initialise the data sets from the HDF file
  */
  
  op_printf(" Initialising data sets and reading from flow file... \n");

  // Sets

  op_set elements   = op_decl_set_hdf5(file, "elements");
  op_set edges      = op_decl_set_hdf5(file, "edges");
  op_set fluxpoints = op_decl_set_hdf5(file, "fluxpoints");
  op_set solnpoints = op_decl_set_hdf5(file, "solnpoints");
  op_set intfces    = op_decl_set_hdf5(file, "intfces");

  // Maps

  op_map pEdgesElements      = op_decl_map_hdf5(edges,      elements,    2, file, "pEdgesElements");
  op_map pElementsSolnPoints = op_decl_map_hdf5(elements,   solnpoints, 25, file, "pElementsSolnPoints");
  op_map pElementsFluxPoints = op_decl_map_hdf5(elements,   fluxpoints, 20, file, "pElementsFluxPoints");
  op_map pInfceFluxPoints    = op_decl_map_hdf5(intfces,    fluxpoints,  2, file, "pInfceFluxPoints");

  // Data to be read from file

  op_dat p_xCentrs     = op_decl_dat_hdf5(elements,       nDim, "double", file, "p_xCentrs");
  op_dat p_xSolnPoints = op_decl_dat_hdf5(solnpoints,     nDim, "double", file, "p_xSolnPoints");
  op_dat p_qSolnPoints = op_decl_dat_hdf5(solnpoints,     nPde, "double", file, "p_qSolnPoints");
  op_dat p_mSolnPoints = op_decl_dat_hdf5(solnpoints,nDim*nDim, "double", file, "p_mSolnPoints");
  op_dat p_jSolnPoints = op_decl_dat_hdf5(solnpoints,        1, "double", file, "p_jSolnPoints");
  op_dat p_nIntfces    = op_decl_dat_hdf5(intfces,        nDim, "double", file, "p_nIntfces");

  /*
  // Data to be generated locally
  */ 

  double *oSolnPoints; oSolnPoints = (double *)malloc(nPde * op_get_size(solnpoints) * sizeof(double));
  op_dat p_oSolnPoints = op_decl_dat(solnpoints,      nPde, "double", oSolnPoints, "p_oSolnPoints");

  double *qFluxPoints; qFluxPoints = (double *)malloc(nPde * op_get_size(fluxpoints) * sizeof(double));
  op_dat p_qFluxPoints = op_decl_dat(fluxpoints,      nPde, "double", qFluxPoints, "p_qFluxPoints");

  double *fSolnPoints; fSolnPoints = (double *)malloc(nDim*nPde * op_get_size(solnpoints) * sizeof(double));
  op_dat p_fSolnPoints = op_decl_dat(solnpoints, nDim*nPde, "double", fSolnPoints, "p_fSolnPoints");

  double *dSolnPoints; dSolnPoints = (double *)malloc(nPde * op_get_size(solnpoints) * sizeof(double));
  op_dat p_dSolnPoints = op_decl_dat(solnpoints,      nPde, "double", dSolnPoints, "p_dSolnPoints");

  double *iFluxPoints; iFluxPoints = (double *)malloc(nPde * op_get_size(fluxpoints) * sizeof(double));
  op_dat p_iFluxPoints = op_decl_dat(fluxpoints,      nPde, "double", iFluxPoints, "p_iFluxPoints");

  double *fFluxPoints; fFluxPoints = (double *)malloc(nPde * op_get_size(fluxpoints) * sizeof(double));
  op_dat p_fFluxPoints = op_decl_dat(fluxpoints,      nPde, "double", fFluxPoints, "p_fFluxPoints");

  double *cSolnPoints; cSolnPoints = (double *)malloc(nPde * op_get_size(solnpoints) * sizeof(double));
  op_dat p_cSolnPoints = op_decl_dat(solnpoints,      nPde, "double", cSolnPoints, "p_cSolnPoints");

  op_printf(" ...done. \n");

  /*
  // Print out some information from OP2 about the setup
  */

  op_diagnostic_output();

  /* 
  // Do the partitioning of the mesh using the appropriate method
  */

  //op_partition("PARMETIS", "KWAY", elements, pEdgesElements, p_xPoints);
  //op_partition("BLOCK", "ANY", elements, pEdgesElements, p_xCentrs);
  op_partition("PTSCOTCH", "KWAY", elements, pEdgesElements, p_xCentrs);

  /*
  // initialise timers for total execution wall time
  */ 

  op_timers(&cpu_t1, &wall_t1);

  /* 
  // The main time-marching loop
  */ 

  niter = 20000;

  for (int iter = 1; iter <= niter; iter++) {

    /*
    // Step 0 --- Save the old conserved variables at the solution points
    */

    op_par_loop(qts2q0s, "qts2q0s", solnpoints, 
    		          op_arg_dat(p_qSolnPoints, -1, OP_ID,      nPde, "double", OP_READ),
    		          op_arg_dat(p_oSolnPoints, -1, OP_ID,      nPde, "double", OP_WRITE) );

    /*
    // Step 1 --- Project the solution to the flux points
    */

    op_par_loop(qts2qtf, "qts2qtf", elements, 
    		          op_arg_dat(p_qSolnPoints, -25, pElementsSolnPoints, nPde, "double", OP_READ),
    		          op_arg_dat(p_qFluxPoints, -20, pElementsFluxPoints, nPde, "double", OP_WRITE) );

    /*
    // Step 2 --- Compute the fluxes at the solution points
    */

    op_par_loop(qts2fhs, "qts2fhs", solnpoints, 
    		          op_arg_dat(p_qSolnPoints, -1, OP_ID,      nPde, "double", OP_READ),
    		          op_arg_dat(p_mSolnPoints, -1, OP_ID, nDim*nDim, "double", OP_READ),
    		          op_arg_dat(p_fSolnPoints, -1, OP_ID, nDim*nPde, "double", OP_WRITE) );

    /*
    // Step 3 --- Compute the divergence of the fluxes at the solution points
    */

    op_par_loop(fhs2dts, "fhs2dts", elements, 
    		          op_arg_dat(p_fSolnPoints, -25, pElementsSolnPoints, nDim*nPde, "double", OP_READ),
    		          op_arg_dat(p_dSolnPoints, -25, pElementsSolnPoints,      nPde, "double", OP_WRITE) );

    /*
    // Step 4 --- Compute the common interface fluxes at the flux points
    */

    op_par_loop(qtf2itf, "qtf2itf", intfces,
    		          op_arg_dat(p_qFluxPoints, -2, pInfceFluxPoints, nPde, "double", OP_READ),
    		          op_arg_dat(p_nIntfces,    -1, OP_ID,            nDim, "double", OP_READ),
    		          op_arg_dat(p_iFluxPoints, -2, pInfceFluxPoints, nPde, "double", OP_WRITE) );

    /*
    // Step 5 --- Compute the projection of the normal fluxes to the flux points
    */

    op_par_loop(fhs2ftf, "fhs2ftf", elements, 
    		          op_arg_dat(p_fSolnPoints, -25, pElementsSolnPoints, nDim*nPde, "double", OP_READ),
    		          op_arg_dat(p_fFluxPoints, -20, pElementsFluxPoints,      nPde, "double", OP_WRITE) );

    /*
    // Step 6 --- Correct the divergence of the fluxes at the solution points
    */

    op_par_loop(dts2cts, "dts2cts", elements,
    		          op_arg_dat(p_iFluxPoints, -20, pElementsFluxPoints, nPde, "double", OP_READ),
    		          op_arg_dat(p_fFluxPoints, -20, pElementsFluxPoints, nPde, "double", OP_READ), 
    		          op_arg_dat(p_dSolnPoints, -25, pElementsSolnPoints, nPde, "double", OP_READ),
    		          op_arg_dat(p_cSolnPoints, -25, pElementsSolnPoints, nPde, "double", OP_WRITE) );

    /*
    // Step 7 --- Update the solution at the solution points
    */

    op_par_loop(cts2qts, "cts2qts", solnpoints, 
    		          op_arg_dat(p_cSolnPoints, -1, OP_ID,      nPde, "double", OP_READ),
    		          op_arg_dat(p_jSolnPoints, -1, OP_ID,         1, "double", OP_READ),
    		          op_arg_dat(p_oSolnPoints, -1, OP_ID,      nPde, "double", OP_READ),
    		          op_arg_dat(p_qSolnPoints, -1, OP_ID,      nPde, "double", OP_WRITE) );

  }

  /*
  // write out to the HDF file
  */
  op_fetch_data_hdf5_file(p_xSolnPoints, "new_data.h5");
  op_fetch_data_hdf5_file(p_qSolnPoints, "new_data.h5");

  /*
  // check the timers for total execution wall time
  */ 

  op_timers(&cpu_t2, &wall_t2);

  /*
  // Print out the timers
  */ 

  op_timing_output();

  /*
  // Say hello
  */

  WL_goodbye();

  /*
  // And exit OP2
  */

  op_exit();
}
