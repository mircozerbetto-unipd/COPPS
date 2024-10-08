---------------------------------------------------------
| C++OPPS v 2.2 - Installation notes                    |
|                                                       |
| Mirco Zerbetto - Dipartimento di Scienze Chimiche     |
| Universita' degli Studi di Padova, 35131, Padova (IT) |
| email: mirco.zerbetto@unipd.it                        |
---------------------------------------------------------

A. INTRODUCTION

  C++OPPS is a C++ software package for the calculation of
  NMR relaxation data of molecules in liquid phases, based
  on a stochastic modeling of the dynamics of the system.
 
  These stochastic models are implemented:

  - SRLS: the Slowly Relaxing Local Structure model. This is a
         phenomenological model where the dynamics of the molecule
         is described in terms of two coupled rotators. One rotator
         describes the global tumbling of the molecule, while the
         other rotator collects local internal motions around the NMR
         probe. Coupling is originated by a potential of mean force,
         while no hydrodynamic coupling is included in the model.

  - TS-SRLS: the Two States SRLS is a SRLS model augmented by a jump motion
             between 2 discrete states that differ in the orientation of
             the magnetic tensors.

  -FB1: the Flexible Body model with 1 internal degree of freedom. This is
        a molecular model, in which the stochastic coordiantes are the
        global tumbling and the rotation about one dihedral angle.

  -FB2: same as FB1, but 2 dihedral angles are included in the relevant
        dynamics.

        NB: the implementation of this model has not been optimized yet.
            Therefore, calculations may take long time.

  - 3S-FB: 3 Sites-Flexible Body is an augmented FB1 model coupled to a
           random jump among 3 discrete states that differ for the
           orientation of the magnetic tensors.

B. COMPILING C++OPPS

  The copps copps program needs the following libraries to run

  - libcqp++.a
  - libminpack_jac.a
  - libmv.a, libsparse.a, libspblas.a
  - liblevmar.a
  - libbessel.a
  - libjemdi.a, libcmrg.a
  - libcubature.a

  Sources are located in the lib/src/ directory of the C+OPPS-v2.2 package.

  Optionally, it is possible to compile copps program to run in parallel in
  multiprocessor workstations and multicore new-generation CPUs. Copps is
  parallelized under the MPI paradgm. A MPI library implementation is needed.

   E.g.: MPICH (http://www-unix.mcs.anl.gov/mpi/mpich1/download.html)
         OPENMPI (http://www.open-mpi.org/)

  To compile copps it should be sufficient to run:

   sh ./build.sh

  from the C++OPPS-v2.2 software package directory, and follow the instructions.

  If building the parallelized version, the build script assumes that the
  mpicxx and mpif77 compilers can be called from the bash shell.

  Lapack and Blas libraries are necessary and the build script asks for their
  location (the default the path is /usr/lib64/).
  
C. RUNNING C++OPPS

  The latest version of C++OPPS is run only by command line. The binary
  executable, named 'copps', is found in the src/ directory if the compilation
  ends succefully. This executable can be moved wherever it is preferred.

  The program expects that in the working directory (say its absolute path is <WDIR>)
  it is present at least a file named <PROJECT>_copps.input, where PROJECT is
  a label that is used to distringuish one calculation from another (e.g., if one
  is calcultating NMR data for two different proteins). The calculation is run
  by calling

   <PATH_TO>/copp WDIR PROJECT

  and the results will be placed in the file <WDIR>/<PROJECT>_copps.output.

  If fitting has to be done, experimental data need to be specified in a
  file called <WDIR>/<PROJECT>_copps.exp.

  Please, see the examples for writing the input and exp files.
  
  If C++OPPS was compiled with MPI, than the call changes to

   mpirun -np <P> <PATH_TO>/copp WDIR PROJECT

  where P is the number of cores on which the calculation will be parallelized.


