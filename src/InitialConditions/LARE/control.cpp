#include "lareic.h"
#include "lareBoundaryClass.h"

namespace LARE{

    void LARE3DInitialConditions::control_variables([[maybe_unused]] SAMS::harness &harnessRef, [[maybe_unused]] simulationData &data){
        std::cout << "Setting control variables for LARE3D simulation" << std::endl;
        data.nx = 128; // Number of cells in the x-direction
        data.ny = 128; // Number of cells in the y-direction
        data.nz = 128; // Number of cells in the z-direction

        data.dt_multiplier = 0.8; // Default multiplier for time step
        data.dt = 0.0;

        // Maximum number of iterations; if nsteps < 0, run until t_end
        data.nsteps = 10;
        data.t_end = 60.0 * 60.0 * 24.0 * 24.0; // One day in seconds

        // Geometry options: cartesian, cylindrical, spherical
        data.geometry = geometryType::Cartesian;

        // Domain limits
        data.x_min = -1.0e6;
        data.x_max = 1.0e6;
        data.y_min = -1.0e6;
        data.y_max = 1.0e6;
        data.z_min = -1.0e6;
        data.z_max = 1.0e6;

        // Boundary conditions
        data.xbc_min = BCType::BC_OTHER;
        data.xbc_max = BCType::BC_OTHER;
        data.ybc_min = BCType::BC_OTHER;
        data.ybc_max = BCType::BC_OTHER;
        data.zbc_min = BCType::BC_OTHER;
        data.zbc_max = BCType::BC_OTHER;

        // Grid stretching
        data.x_stretch = false;
        data.y_stretch = false;
        data.z_stretch = false;

        // Shock viscosity coefficients
        data.visc1 = 0.1;
        data.visc2 = 1.0;

        // Ratio of specific heat capacities
        data.gas_gamma = 1.4;

        // Average mass of an ion in proton masses
        data.mf = 1.2;

        // Resistive MHD options
        data.resistiveMHD = false;
        data.eta_background = 1.e-10;
        data.j_max = 1.0;
        data.eta0 = 2.e-10;

        // Remap kinetic energy correction
        data.rke = true;

        // Output frequency and directory
        data.dt_snapshots = 10.0;
    }

 void LARE3DInitialConditions::initial_conditions([[maybe_unused]] SAMS::harness &harness, [[maybe_unused]] simulationData &data)
  {

    SAMS::cout << "Setting up initial conditions" << std::endl;
    // Set initial conditions for the LARE3D
    pw::assign(data.vx, 0.0);
    pw::assign(data.vy, 0.0);
    pw::assign(data.vz, 0.0);

    pw::assign(data.rho, 1.0);
    pw::assign(data.energy_electron, 1.0);
    pw::assign(data.energy_ion, 1.0);

    pw::assign(data.bx, 0.0);
    pw::assign(data.by, 0.0);
    pw::assign(data.bz, 0.0);

    return;

    /*T_dataType v0 = 0.e3;
    T_dataType a0 = 1.0e5;
    T_dataType a2 = a0 * a0;
    T_dataType amp = 0.5;

    T_dataType xcentre = 0.0;
    T_dataType ycentre = 0.5e6;
    T_dataType zcentre = 0.0;

    if (harness.MPIManager.getRank() == 0 && harness.MPIManager.getSize() > 1)
    {
      amp = 0.0;
    }

    // Set the initial thermal energy of electrons and ions
    T_dataType T0 = 1.e6;
    T_dataType energy = 0.5 * kb_si * T0 / mh_si / (data.gas_gamma - 1.0);
    pw::assign(data.energy_electron, energy);
    pw::assign(data.energy_ion, energy);

    pw::applyKernel(
        LAMBDA(T_indexType ix, T_indexType iy, T_indexType iz) {
          T_dataType x2 = (data.xc(ix) - xcentre) * (data.xc(ix) - xcentre);
          T_dataType y2 = (data.yc(iy) - ycentre) * (data.yc(iy) - ycentre);
          T_dataType z2 = (data.zc(iz) - zcentre) * (data.zc(iz) - zcentre);
          T_dataType r2 = x2 + y2 + z2;
          T_dataType v = v0 * std::exp(-r2 / a2);
          data.energy_electron(ix, iy, iz) *= (1.0 + amp * std::exp(-r2 / a2));
          data.energy_ion(ix, iy, iz) *= (1.0 + amp * std::exp(-r2 / a2));
        },
        pw::Range(-1, data.nx + 2),
        pw::Range(-1, data.ny + 2),
        pw::Range(-1, data.nz + 2));

    T_dataType bmult = 000.0;
    pw::assign(data.bx, 0.01 * bmult);
    pw::assign(data.by, 0.00 * bmult);
    pw::assign(data.bz, 0.00 * bmult);
    // Set the initial density field in kg/m^3
    pw::assign(data.rho, 1.0e-6);

    if (data.rke)
      pw::assign(data.delta_ke, 0.0);*/
  }

}//namespace LARE