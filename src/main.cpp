#define PRINT_PARALLELIZATION_INFO

#include "pp/range.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include "harness.h"
#include "shared_data.h"
#include "include/timer.h"
#include "axisRegistry.h"
#include "variableRegistry.h"
#include "mpiManager.h"
#include "welcome.h"

#include "builtInBoundaryConditions.h"

#include "runner.h"

int main(int argc, char *argv[]){

    //Initialize MPI
    SAMS::MPI::initialize(argc, argv);

    SAMS::harness harness;

    SAMS::runner<
        LARE::LARE3D> runner;

    runner.initialize(argc, argv);
    runner.activateSimulation("LARE3D");
    runner.initializeSimulations();

    SAMS::printWelcomeMessage();
    //MPI auto decomposition
    auto& mpi = harness.MPIManager;
    auto& axisReg = harness.axisRegistry;
    auto& varReg = harness.variableRegistry;
    mpi.autoDecomposition({false,false,false});
    //Initialize portable wrapper
    portableWrapper::initialize(argc, argv);
    SAMS::finishWelcomeMessage();

    //Create the simulation (LARE) and data objects
    LARE::LARE3D S(harness);
    LARE::simulationData data;

    //Setup control variables
    S.controlvariables(data);
    data.visc2_norm=data.visc2;

    //Register axes and attach them to MPI dimensions
    axisReg.registerAxis("X", SAMS::MPIAxis(0));
    axisReg.registerAxis("Y", SAMS::MPIAxis(1));
    axisReg.registerAxis("Z", SAMS::MPIAxis(2));
    //Tell LARE to register its variables
    S.registerVars();
    //Other simulations would register their variables here too

    //Set the axis domains and decompose them
    axisReg.setDomain("X", data.nx, data.x_min, data.x_max);
    axisReg.setDomain("Y", data.ny, data.y_min, data.y_max);
    axisReg.setDomain("Z", data.nz, data.z_min, data.z_max);

    mpi.decomposeAllAxes();

    varReg.allocateAll();

    //Tell LARE to grab the shared allocated variables
		S.allocate(data);
    //Tell LARE to set up its grid
    S.grid(data);

		portableWrapper::fence();
    S.initial_conditions(data);
		portableWrapper::fence();
    S.boundary_conditions(data);
    portableWrapper::fence();
    timer t;
    t.begin("Main Loop");
    data.step=0;

    while (true)
    {
      SAMS::cout << data.step << " " << data.time << std::endl;      
      //if (data.step%10==0) S.output(data);
      if ((data.step >= data.nsteps && data.nsteps >= 0) || (data.time >= data.t_end))
        break;

      S.lagrangian_step(data);    // lagran.cpp
      S.eulerian_remap(data); // remap.cpp
      data.step++;
      if (data.rke) S.energy_correction(data); // diagnostics.cpp
      S.eta_calc(data);            // lagran.cpp
    }
    t.end();

		S.output(data);

		S.manager.clear();
    axisReg.finalize();
    varReg.finalize();
    portableWrapper::finalize();
    SAMS::MPI::finalize();

}
