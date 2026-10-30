//#include "lareic.h"
//#include "lareBoundaryClass.h"
#ifdef USE_KOKKOS
#include "shared_data.h"

namespace LARE{

    void LARE3DKokkos::defaultValues(simulationData &data){

        data.dt_multiplier = 0.8; // Default multiplier for time step

        // Geometry options: cartesian, cylindrical, spherical
        data.geometry = geometryType::Cartesian;

        // Shock viscosity coefficients
        data.visc1 = 0.1;
        data.visc2 = 1.0;

        // Ratio of specific heat capacities
        data.gas_gamma = 1.4;

        // Average mass of an ion in proton masses
        data.mf = 1.2;

        data.rke = false; // Remap kinetic energy correction off by default

        // Physical constants
        data.mu0 = mu0_si;

        // Initialize time=0
        data.time = 0.0;
    }

 void LARE3DKokkos::defaultVariables([[maybe_unused]] simulationData &data)
  {

    // Set default variable values
    Kokkos::deep_copy(data.vx, 0.0);
    Kokkos::deep_copy(data.vy, 0.0);
    Kokkos::deep_copy(data.vz, 0.0);

    Kokkos::deep_copy(data.rho, 0.0);
    Kokkos::deep_copy(data.energy_electron, 0.0);
    Kokkos::deep_copy(data.energy_ion, 0.0);

    Kokkos::deep_copy(data.bx, 0.0);
    Kokkos::deep_copy(data.by, 0.0);
    Kokkos::deep_copy(data.bz, 0.0);

    if (data.rke)
    {
        Kokkos::deep_copy(data.delta_ke, 0.0);
    }

  }

}//namespace LARE

#endif // USE_KOKKOS