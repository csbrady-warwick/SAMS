/*
 *    Copyright 2025 SAMS Team
   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.
*/
#ifdef USE_KOKKOS
#include "shared_data.h"
#include "variableRegistry.h"
#include "axisRegistry.h"

namespace LARE
{

    namespace pw = portableWrapper;

    /**
     * Register LARE's axes with the axis registry.
     */
    void LARE3DKokkos::registerAxes(SAMS::harness &harness)
    {
        auto &axisReg = harness.axisRegistry;
        axisReg.registerAxis("X", SAMS::MPIAxis(0));
        axisReg.registerAxis("Y", SAMS::MPIAxis(1));
        axisReg.registerAxis("Z", SAMS::MPIAxis(2));
    }

    /**
     * Register variables with the portable array manager.
     */
    void LARE3DKokkos::registerVariables(SAMS::harness &harness)
    {

        auto &varRegistry = harness.variableRegistry;

        const int ghosts = 2; // 2 Ghost cells at top and bottom of each dimension

        varRegistry.registerVariable<T_dataType>("energy_electron", pw::arrayTags::accelerated, SAMS::dimension("X", ghosts), SAMS::dimension("Y", ghosts), SAMS::dimension("Z", ghosts));

        varRegistry.registerVariable<T_dataType>("energy_ion", pw::arrayTags::accelerated, SAMS::dimension("X", ghosts), SAMS::dimension("Y", ghosts), SAMS::dimension("Z", ghosts));

        varRegistry.registerVariable<T_dataType>("rho", pw::arrayTags::accelerated, SAMS::dimension("X", ghosts), SAMS::dimension("Y", ghosts), SAMS::dimension("Z", ghosts));

        varRegistry.registerVariable<T_dataType>("vx", pw::arrayTags::accelerated, SAMS::dimension("X", ghosts, SAMS::staggerType::HALF_CELL), SAMS::dimension("Y", ghosts, SAMS::staggerType::HALF_CELL), SAMS::dimension("Z", ghosts, SAMS::staggerType::HALF_CELL));

        varRegistry.registerVariable<T_dataType>("vy", pw::arrayTags::accelerated, SAMS::dimension("X", ghosts, SAMS::staggerType::HALF_CELL), SAMS::dimension("Y", ghosts, SAMS::staggerType::HALF_CELL), SAMS::dimension("Z", ghosts, SAMS::staggerType::HALF_CELL));

        varRegistry.registerVariable<T_dataType>("vz", pw::arrayTags::accelerated, SAMS::dimension("X", ghosts, SAMS::staggerType::HALF_CELL), SAMS::dimension("Y", ghosts, SAMS::staggerType::HALF_CELL), SAMS::dimension("Z", ghosts, SAMS::staggerType::HALF_CELL));

        varRegistry.registerVariable<T_dataType>("LARE/vx1", pw::arrayTags::accelerated, SAMS::dimension("X", ghosts, SAMS::staggerType::HALF_CELL), SAMS::dimension("Y", ghosts, SAMS::staggerType::HALF_CELL), SAMS::dimension("Z", ghosts, SAMS::staggerType::HALF_CELL));

        varRegistry.registerVariable<T_dataType>("LARE/vy1", pw::arrayTags::accelerated, SAMS::dimension("X", ghosts, SAMS::staggerType::HALF_CELL), SAMS::dimension("Y", ghosts, SAMS::staggerType::HALF_CELL), SAMS::dimension("Z", ghosts, SAMS::staggerType::HALF_CELL));

        varRegistry.registerVariable<T_dataType>("LARE/vz1", pw::arrayTags::accelerated, SAMS::dimension("X", ghosts, SAMS::staggerType::HALF_CELL), SAMS::dimension("Y", ghosts, SAMS::staggerType::HALF_CELL), SAMS::dimension("Z", ghosts, SAMS::staggerType::HALF_CELL));

        varRegistry.registerVariable<T_dataType>("bx", pw::arrayTags::accelerated, SAMS::dimension("X", ghosts, SAMS::staggerType::HALF_CELL), SAMS::dimension("Y", ghosts), SAMS::dimension("Z", ghosts));

        varRegistry.registerVariable<T_dataType>("by", pw::arrayTags::accelerated, SAMS::dimension("X", ghosts), SAMS::dimension("Y", ghosts, SAMS::staggerType::HALF_CELL), SAMS::dimension("Z", ghosts));

        varRegistry.registerVariable<T_dataType>("bz", pw::arrayTags::accelerated, SAMS::dimension("X", ghosts), SAMS::dimension("Y", ghosts), SAMS::dimension("Z", ghosts, SAMS::staggerType::HALF_CELL));

        varRegistry.registerVariable<T_dataType>("LARE/dm", pw::arrayTags::accelerated, SAMS::dimension("X", ghosts, SAMS::staggerType::CENTRED), SAMS::dimension("Y", ghosts, SAMS::staggerType::CENTRED), SAMS::dimension("Z", ghosts, SAMS::staggerType::CENTRED));
    }

    /**
     * Allocate the data arrays for the LARE3DKokkos.
     * This allocates the permanent state arrays that are used throughout the LARE3DKokkos.
     */
    void LARE3DKokkos::allocate(SAMS::harness &harness, simulationData &data)
    {
        T_sizeType nx, ny, nz;

        auto &axRegistry = harness.axisRegistry;
        auto &varRegistry = harness.variableRegistry;
        // Centred since LARE thinks in terms of cell centres for nx, ny, nz
        nx = axRegistry.getLocalDomainElements("X", SAMS::staggerType::CENTRED);
        ny = axRegistry.getLocalDomainElements("Y", SAMS::staggerType::CENTRED);
        nz = axRegistry.getLocalDomainElements("Z", SAMS::staggerType::CENTRED);

        // Get the ranges for the whole local domain
        data.xcLocalRange = axRegistry.getLocalRange("X", SAMS::staggerType::CENTRED);
        data.ycLocalRange = axRegistry.getLocalRange("Y", SAMS::staggerType::CENTRED);
        data.zcLocalRange = axRegistry.getLocalRange("Z", SAMS::staggerType::CENTRED);
        data.xbLocalRange = axRegistry.getLocalRange("X", SAMS::staggerType::HALF_CELL);
        data.ybLocalRange = axRegistry.getLocalRange("Y", SAMS::staggerType::HALF_CELL);
        data.zbLocalRange = axRegistry.getLocalRange("Z", SAMS::staggerType::HALF_CELL);

        // Get the ranges for the ghost cells for centred axes
        data.xcminBCRange = axRegistry.getLocalNonDomainRange("X", SAMS::staggerType::CENTRED, SAMS::domain::edges::lower );
        data.xcmaxBCRange = axRegistry.getLocalNonDomainRange("X", SAMS::staggerType::CENTRED, SAMS::domain::edges::upper);
        data.ycminBCRange = axRegistry.getLocalNonDomainRange("Y", SAMS::staggerType::CENTRED, SAMS::domain::edges::lower);
        data.ycmaxBCRange = axRegistry.getLocalNonDomainRange("Y", SAMS::staggerType::CENTRED, SAMS::domain::edges::upper);
        data.zcminBCRange = axRegistry.getLocalNonDomainRange("Z", SAMS::staggerType::CENTRED, SAMS::domain::edges::lower);
        data.zcmaxBCRange = axRegistry.getLocalNonDomainRange("Z", SAMS::staggerType::CENTRED, SAMS::domain::edges::upper);

        // Get the ranges for the ghost cells for half cell staggered axes
        data.xbminBCRange = axRegistry.getLocalNonDomainRange("X", SAMS::staggerType::HALF_CELL, SAMS::domain::edges::lower);
        data.xbmaxBCRange = axRegistry.getLocalNonDomainRange("X", SAMS::staggerType::HALF_CELL, SAMS::domain::edges::upper);
        data.ybminBCRange = axRegistry.getLocalNonDomainRange("Y", SAMS::staggerType::HALF_CELL, SAMS::domain::edges::lower);
        data.ybmaxBCRange = axRegistry.getLocalNonDomainRange("Y", SAMS::staggerType::HALF_CELL, SAMS::domain::edges::upper);
        data.zbminBCRange = axRegistry.getLocalNonDomainRange("Z", SAMS::staggerType::HALF_CELL, SAMS::domain::edges::lower);
        data.zbmaxBCRange = axRegistry.getLocalNonDomainRange("Z", SAMS::staggerType::HALF_CELL, SAMS::domain::edges::upper);


        // Get the ranges for the actual domain (no ghost cells)
        data.xcLocalDomainRange = axRegistry.getLocalDomainRange("X", SAMS::staggerType::CENTRED);
        data.ycLocalDomainRange = axRegistry.getLocalDomainRange("Y", SAMS::staggerType::CENTRED);
        data.zcLocalDomainRange = axRegistry.getLocalDomainRange("Z", SAMS::staggerType::CENTRED);
        data.xbLocalDomainRange = axRegistry.getLocalDomainRange("X", SAMS::staggerType::HALF_CELL);
        data.ybLocalDomainRange = axRegistry.getLocalDomainRange("Y", SAMS::staggerType::HALF_CELL);
        data.zbLocalDomainRange = axRegistry.getLocalDomainRange("Z", SAMS::staggerType::HALF_CELL);

        data.nx = nx;
        data.ny = ny;
        data.nz = nz;

        data.nx_global = axRegistry.getDimension("X").getGlobalDomainCount(SAMS::staggerType::CENTRED);
        data.ny_global = axRegistry.getDimension("Y").getGlobalDomainCount(SAMS::staggerType::CENTRED);
        data.nz_global = axRegistry.getDimension("Z").getGlobalDomainCount(SAMS::staggerType::CENTRED);

        using Range = pw::Range;
        // Grab the final variable sizes from the registry and wrap the arrays
        varRegistry.fillKOVArray("energy_electron", data.energy_electron);
        Kokkos::deep_copy(data.energy_electron, 0.0);
        varRegistry.fillKOVArray("energy_ion", data.energy_ion);
        Kokkos::deep_copy(data.energy_ion, 0.0);
        varRegistry.fillKOVArray("rho", data.rho);
        Kokkos::deep_copy(data.rho, 0.0);
        varRegistry.fillKOVArray("vx", data.vx);
        Kokkos::deep_copy(data.vx, 0.0);
        varRegistry.fillKOVArray("vy", data.vy);
        Kokkos::deep_copy(data.vy, 0.0);
        varRegistry.fillKOVArray("vz", data.vz);
        Kokkos::deep_copy(data.vz, 0.0);
        varRegistry.fillKOVArray("bx", data.bx);
        Kokkos::deep_copy(data.bx, 0.0);
        varRegistry.fillKOVArray("by", data.by);
        Kokkos::deep_copy(data.by, 0.0);
        varRegistry.fillKOVArray("bz", data.bz);
        Kokkos::deep_copy(data.bz, 0.0);
        varRegistry.fillKOVArray("LARE/vx1", data.vx1);
        Kokkos::deep_copy(data.vx1, 0.0);
        varRegistry.fillKOVArray("LARE/vy1", data.vy1);
        Kokkos::deep_copy(data.vy1, 0.0);
        varRegistry.fillKOVArray("LARE/vz1", data.vz1);
        Kokkos::deep_copy(data.vz1, 0.0);
        varRegistry.fillKOVArray("LARE/dm", data.dm);
        Kokkos::deep_copy(data.dm, 0.0);

        data.isxLB = harness.MPIManager.isEdge(0, SAMS::domain::edges::lower);
        data.isxUB = harness.MPIManager.isEdge(0, SAMS::domain::edges::upper);
        data.isyLB = harness.MPIManager.isEdge(1, SAMS::domain::edges::lower);
        data.isyUB = harness.MPIManager.isEdge(1, SAMS::domain::edges::upper);
        data.iszLB = harness.MPIManager.isEdge(2, SAMS::domain::edges::lower);
        data.iszUB = harness.MPIManager.isEdge(2, SAMS::domain::edges::upper);

        SAMS::debugAll3 << "Edge detection: "
                        << " XLB: " << data.isxLB << " XUB: " << data.isxUB
                        << " YLB: " << data.isyLB << " YUB: " << data.isyUB
                        << " ZLB: " << data.iszLB << " ZUB: " << data.iszUB
                        << std::endl;

        manager.allocateKOV(data.p_visc, data.xcLocalRange, data.ycLocalRange, data.zcLocalRange);
        manager.allocateKOV(data.eta, data.xcLocalRange, data.ycLocalRange, data.zcLocalRange);
        manager.allocateKOV(data.dxab, data.xbLocalRange, data.ycLocalRange, data.zcLocalRange);
        manager.allocateKOV(data.dyab, data.xcLocalRange, data.ybLocalRange, data.zcLocalRange);
        manager.allocateKOV(data.dzab, data.xcLocalRange, data.ycLocalRange, data.zbLocalRange);
        manager.allocateKOV(data.dxac, data.xcLocalRange, data.ycLocalRange, data.zcLocalRange);
        manager.allocateKOV(data.dyac, data.xcLocalRange, data.ycLocalRange, data.zcLocalRange);
        manager.allocateKOV(data.dzac, data.xcLocalRange, data.ycLocalRange, data.zcLocalRange);
        manager.allocateKOV(data.cv, data.xcLocalRange, data.ycLocalRange, data.zcLocalRange);
        manager.allocateKOV(data.cv1, data.xcLocalRange, data.ycLocalRange, data.zcLocalRange);
        manager.allocateKOV(data.cvc, data.xcLocalRange, data.ycLocalRange, data.zcLocalRange);

        Range xcp = pw::Range(0, data.nx + 1);
        Range ycp = pw::Range(0, data.ny + 1);
        Range zcp = pw::Range(0, data.nz + 1);
        Range ycpp = pw::Range(0, data.ny + 2);
        Range zcpp = pw::Range(0, data.nz + 2);
        Range xbp = pw::Range(-1, data.nx + 1);
        Range ybp = pw::Range(-1, data.ny + 1);
        Range zbp = pw::Range(-1, data.nz + 1);
        // Allocate arrays using the portableArrayManager
        manager.allocateKOV(data.bx1, data.xcLocalRange, data.ycLocalRange, data.zcLocalRange);
        manager.allocateKOV(data.by1, data.xcLocalRange, data.ycLocalRange, data.zcLocalRange);
        manager.allocateKOV(data.bz1, data.xcLocalRange, data.ycLocalRange, data.zcLocalRange);
        manager.allocateKOV(data.alpha1, xcp, ycpp, zcpp);
        manager.allocateKOV(data.alpha2, xbp, ycp, zcpp);
        manager.allocateKOV(data.alpha3, data.xcLocalRange, data.ycLocalRange, zcp);
        manager.allocateKOV(data.visc_heat, xcp, ycp, zcp);
        manager.allocateKOV(data.pressure, data.xcLocalRange, data.ycLocalRange, data.zcLocalRange);
        manager.allocateKOV(data.p_e, data.xcLocalRange, data.ycLocalRange, data.zcLocalRange);
        manager.allocateKOV(data.p_i, data.xcLocalRange, data.ycLocalRange, data.zcLocalRange);
        manager.allocateKOV(data.rho_v, xbp, ybp, zbp);
        manager.allocateKOV(data.cv_v, xbp, ybp, zbp);
        manager.allocateKOV(data.fx, data.xbLocalDomainRange, data.ybLocalDomainRange, data.zbLocalDomainRange);
        manager.allocateKOV(data.fy, data.xbLocalDomainRange, data.ybLocalDomainRange, data.zbLocalDomainRange);
        manager.allocateKOV(data.fz, data.xbLocalDomainRange, data.ybLocalDomainRange, data.zbLocalDomainRange);
        manager.allocateKOV(data.fx_visc, data.xbLocalDomainRange, data.ybLocalDomainRange, data.zbLocalDomainRange);
        manager.allocateKOV(data.fy_visc, data.xbLocalDomainRange, data.ybLocalDomainRange, data.zbLocalDomainRange);
        manager.allocateKOV(data.fz_visc, data.xbLocalDomainRange, data.ybLocalDomainRange, data.zbLocalDomainRange);
        manager.allocateKOV(data.flux_x, data.xbLocalDomainRange, data.ybLocalDomainRange, data.zbLocalDomainRange);
        manager.allocateKOV(data.flux_y, data.xbLocalDomainRange, data.ybLocalDomainRange, data.zbLocalDomainRange);
        manager.allocateKOV(data.flux_z, data.xbLocalDomainRange, data.ybLocalDomainRange, data.zbLocalDomainRange);
        manager.allocateKOV(data.curlb, data.xbLocalDomainRange, data.ybLocalDomainRange, data.zbLocalDomainRange);

        axRegistry.fillKOVLocalAxis("X", data.xc, SAMS::staggerType::CENTRED);
        axRegistry.fillKOVLocalAxis("Y", data.yc, SAMS::staggerType::CENTRED);
        axRegistry.fillKOVLocalAxis("Z", data.zc, SAMS::staggerType::CENTRED);
        axRegistry.fillKOVLocalAxis("X", data.xb, SAMS::staggerType::HALF_CELL);
        axRegistry.fillKOVLocalAxis("Y", data.yb, SAMS::staggerType::HALF_CELL);
        axRegistry.fillKOVLocalAxis("Z", data.zb, SAMS::staggerType::HALF_CELL);
        axRegistry.fillPPLocalAxis("X", data.xb_host, SAMS::staggerType::HALF_CELL);
        axRegistry.fillPPLocalAxis("Y", data.yb_host, SAMS::staggerType::HALF_CELL);
        axRegistry.fillPPLocalAxis("Z", data.zb_host, SAMS::staggerType::HALF_CELL);
        axRegistry.fillPPLocalAxis("X", data.xc_host, SAMS::staggerType::CENTRED);
        axRegistry.fillPPLocalAxis("Y", data.yc_host, SAMS::staggerType::CENTRED);
        axRegistry.fillPPLocalAxis("Z", data.zc_host, SAMS::staggerType::CENTRED);
        axRegistry.fillKOVAxis("X", data.xb_global, SAMS::staggerType::CENTRED);
        axRegistry.fillKOVAxis("Y", data.yb_global, SAMS::staggerType::CENTRED);
        axRegistry.fillKOVAxis("Z", data.zb_global, SAMS::staggerType::CENTRED);

        axRegistry.fillKOVLocalDelta("X", data.dxc, SAMS::staggerType::CENTRED);
        axRegistry.fillKOVLocalDelta("Y", data.dyc, SAMS::staggerType::CENTRED);
        axRegistry.fillKOVLocalDelta("Z", data.dzc, SAMS::staggerType::CENTRED);
        axRegistry.fillKOVLocalDelta("X", data.dxb, SAMS::staggerType::HALF_CELL);
        axRegistry.fillKOVLocalDelta("Y", data.dyb, SAMS::staggerType::HALF_CELL);
        axRegistry.fillKOVLocalDelta("Z", data.dzb, SAMS::staggerType::HALF_CELL);

        manager.allocateKOV(data.hy, Range(-2, nx + 2));
        manager.allocateKOV(data.hz, Range(-2, nx + 2), Range(-2, ny + 2));
        manager.allocateKOV(data.hyc, Range(-1, nx + 2));
        manager.allocateKOV(data.hzc, Range(-1, nx + 2), Range(-1, ny + 2));
        manager.allocateKOV(data.hz1, Range(-2, nx + 2), Range(-1, ny + 2));
        manager.allocateKOV(data.hz2, Range(-1, nx + 2), Range(-2, ny + 2));
        manager.allocateKOV(data.x, Range(-2, nx + 2), Range(-2, ny + 2), Range(-2, nz + 2));
        manager.allocateKOV(data.y, Range(-2, nx + 2), Range(-2, ny + 2), Range(-2, nz + 2));
        manager.allocateKOV(data.z, Range(-2, nx + 2), Range(-2, ny + 2), Range(-2, nz + 2));
        manager.allocateKOV(data.xp, Range(-2, nx + 2), Range(-2, ny + 2), Range(-2, nz + 2));
        manager.allocateKOV(data.yp, Range(-2, nx + 2), Range(-2, ny + 2), Range(-2, nz + 2));
        manager.allocateKOV(data.zp, Range(-2, nx + 2), Range(-2, ny + 2), Range(-2, nz + 2));
        if (data.rke)
        {
            manager.allocateKOV(data.delta_ke, Range(-1, nx + 2), Range(-1, ny + 2), Range(-1, nz + 2));
        }

        data.mpiType = SAMS::gettypeRegistry().getMPIType<T_dataType>();
    }

    /**
     * Copy the parameters that the initial conditions set in the data of the portable LARE3D.
     * The number of cells is not copied, because allocate sets it from the decomposed axes.
     */
    void LARE3DKokkos::propagateInitialConditions(simulationData &data, LARE::LARE3D<>::simulationData &lare)
    {
        data.t_end = lare.t_end;
        data.dt_snapshots = lare.dt_snapshots;
        data.x_min = lare.x_min;
        data.x_max = lare.x_max;
        data.y_min = lare.y_min;
        data.y_max = lare.y_max;
        data.z_min = lare.z_min;
        data.z_max = lare.z_max;
        data.dt_multiplier = lare.dt_multiplier;
        data.geometry = lare.geometry;
        data.visc1 = lare.visc1;
        data.visc2 = lare.visc2;
        data.gas_gamma = lare.gas_gamma;
        data.mf = lare.mf;
        data.mu0 = lare.mu0;
        data.resistiveMHD = lare.resistiveMHD;
        data.eta_background = lare.eta_background;
        data.j_max = lare.j_max;
        data.eta0 = lare.eta0;
        data.rke = lare.rke;

        grid(data);
    }

    /**
     * Setup the basic LARE3DKokkos parameters like grid points etc.
     */
    void LARE3DKokkos::grid(simulationData &data)
    {

        pw::portableArrayManager localManager;

        auto hyv = localManager.createKOV<double>(pw::Range(-2, data.nx + 2));
        auto hzv = localManager.createKOV<double>(pw::Range(-2, data.nx + 2), pw::Range(-2, data.ny + 2));

        data.length_x = data.x_max - data.x_min;
        data.length_y = data.y_max - data.y_min;
        data.length_z = data.z_max - data.z_min;

        data.dx = data.length_x / static_cast<T_dataType>(data.nx);
        data.dy = data.length_y / static_cast<T_dataType>(data.ny);
        data.dz = data.length_z / static_cast<T_dataType>(data.nz);

        // The grid axes are already filled, just need to set up the metric terms
        if (data.geometry == geometryType::Cartesian)
        {
            Kokkos::deep_copy(data.hy, 1.0);
            Kokkos::deep_copy(data.hyc, 1.0);
            Kokkos::deep_copy(hyv, 1.0);
            Kokkos::fence();
        }
        else if (data.geometry == geometryType::Cylindrical || data.geometry == geometryType::Spherical)
        {
            auto l1 = KOKKOS_LAMBDA(T_indexType ix)
            {
                data.hy(ix) = std::abs(data.xb(ix));
                data.hyc(ix) = std::abs(data.xc(ix));
                hyv(ix) = (std::abs(data.xc(0)) < none_zero && ix == 1) ? std::abs(data.xc(ix)) : 0.25 * std::abs(data.xc(ix));
            };
            Kokkos::parallel_for("Parallel", Kokkos::RangePolicy<KOKKOS_EXECUTION_SPACE, Kokkos::IndexType<SIGNED_INDEX_TYPE>>(static_cast<SIGNED_INDEX_TYPE>(-2), static_cast<SIGNED_INDEX_TYPE>(data.nx + 2) + 1), l1);
        }
        Kokkos::fence();

        if (data.geometry == geometryType::Cartesian || data.geometry == geometryType::Cylindrical)
        {
            Kokkos::deep_copy(data.hz, 1.0);
            Kokkos::deep_copy(data.hzc, 1.0);
            Kokkos::deep_copy(hzv, 1.0);
            Kokkos::deep_copy(data.hz1, 1.0);
            Kokkos::deep_copy(data.hz2, 1.0);
            Kokkos::fence();
        }
        else if (data.geometry == geometryType::Spherical)
        {
            auto l1 = KOKKOS_LAMBDA(T_indexType iy)
            {
                T_dataType s = std::abs(std::sin(data.yb(iy)));
                T_indexType ix = -2;
                data.hz(ix, iy) = data.hy(ix) * s;
                hzv(ix + 2, iy + 2) = hyv(ix) * s;
                for (ix = -1; ix <= data.nx + 2; ++ix)
                {
                    data.hz(ix, iy) = data.hy(ix) * s;
                    hzv(ix, iy) = hyv(ix) * s;
                    data.hz2(ix, iy) = data.hyc(ix) * s;
                }
            };
            Kokkos::parallel_for("Parallel", Kokkos::RangePolicy<KOKKOS_EXECUTION_SPACE, Kokkos::IndexType<SIGNED_INDEX_TYPE>>(static_cast<SIGNED_INDEX_TYPE>(-2), static_cast<SIGNED_INDEX_TYPE>(data.ny + 2) + 1), l1);

            auto l2 = KOKKOS_LAMBDA(T_indexType iy)
            {
                T_dataType sc = std::abs(std::sin(data.yc(iy)));
                T_indexType ix = -2;
                data.hz1(ix, iy) = data.hy(ix) * sc;
                for (ix = -1; ix <= data.nx + 2; ++ix)
                {
                    data.hz1(ix, iy) = data.hy(ix) * sc;
                    data.hzc(ix, iy) = data.hyc(ix) * sc;
                }
            };
            Kokkos::parallel_for("Parallel", Kokkos::RangePolicy<KOKKOS_EXECUTION_SPACE, Kokkos::IndexType<SIGNED_INDEX_TYPE>>(static_cast<SIGNED_INDEX_TYPE>(-1), static_cast<SIGNED_INDEX_TYPE>(data.ny + 2) + 1), l2);
            Kokkos::fence();
        }

        {
            // Cell centred areas and volumes
            auto l = KOKKOS_LAMBDA(T_indexType iy, T_indexType iz)
            {
                T_dataType dy = data.dyb(iy);
                T_dataType dz = data.dzb(iz);
                T_dataType dydz = dy * dz;
                data.dxab(-2, iy, iz) = dydz * data.hy(-2) * data.hz1(-2, iy);
                for (SIGNED_INDEX_TYPE ix = -1; ix <= data.nx + 2; ++ix)
                {
                    T_dataType dx = data.dxb(ix);
                    data.dxab(ix, iy, iz) = dydz * data.hy(ix) * data.hz1(ix, iy);
                    data.dyab(ix, iy, iz) = dx * dz * data.hz2(ix, iy);
                    data.dzab(ix, iy, iz) = dx * dy * data.hy(ix);
                    data.cv(ix, iy, iz) = dx * dydz * data.hyc(ix) * data.hzc(ix, iy);
                }
            };
            Kokkos::parallel_for("Parallel", Kokkos::MDRangePolicy<KOKKOS_EXECUTION_SPACE, Kokkos::Rank<2, LARE_ITERATE, LARE_ITERATE>, Kokkos::IndexType<SIGNED_INDEX_TYPE>>({static_cast<SIGNED_INDEX_TYPE>(-1), static_cast<SIGNED_INDEX_TYPE>(-1)}, {static_cast<SIGNED_INDEX_TYPE>(data.ny + 2) + 1, static_cast<SIGNED_INDEX_TYPE>(data.nz + 2) + 1}), l);
            Kokkos::fence();
        }

        {
            // Fix negative y boundaries
            Kokkos::parallel_for("Parallel", Kokkos::MDRangePolicy<KOKKOS_EXECUTION_SPACE, Kokkos::Rank<2, LARE_ITERATE, LARE_ITERATE>, Kokkos::IndexType<SIGNED_INDEX_TYPE>>({static_cast<SIGNED_INDEX_TYPE>(-1), static_cast<SIGNED_INDEX_TYPE>(-1)}, {static_cast<SIGNED_INDEX_TYPE>(data.nx + 2) + 1, static_cast<SIGNED_INDEX_TYPE>(data.nz + 2) + 1}), KOKKOS_LAMBDA(T_indexType ix, T_indexType iz) { data.dyab(ix, -2, iz) = data.dxb(ix) * data.dzb(iz) * data.hz2(ix, -2); });
            Kokkos::fence();
        }
        {
            // Fix negative z boundaries
            auto l = KOKKOS_LAMBDA(T_indexType ix, T_indexType iy)
            {
                T_dataType dx = data.dxb(ix);
                T_dataType dy = data.dyc(iy);
                data.dzab(ix, iy, -2) = dx * dy * data.hyc(ix);
            };
            Kokkos::parallel_for("Parallel", Kokkos::MDRangePolicy<KOKKOS_EXECUTION_SPACE, Kokkos::Rank<2, LARE_ITERATE, LARE_ITERATE>, Kokkos::IndexType<SIGNED_INDEX_TYPE>>({static_cast<SIGNED_INDEX_TYPE>(-1), static_cast<SIGNED_INDEX_TYPE>(-1)}, {static_cast<SIGNED_INDEX_TYPE>(data.nx + 2) + 1, static_cast<SIGNED_INDEX_TYPE>(data.ny + 2) + 1}), l);
            Kokkos::parallel_for("Parallel", Kokkos::MDRangePolicy<KOKKOS_EXECUTION_SPACE, Kokkos::Rank<2, LARE_ITERATE, LARE_ITERATE>, Kokkos::IndexType<SIGNED_INDEX_TYPE>>({static_cast<SIGNED_INDEX_TYPE>(-1), static_cast<SIGNED_INDEX_TYPE>(-1)}, {static_cast<SIGNED_INDEX_TYPE>(data.nx + 2) + 1, static_cast<SIGNED_INDEX_TYPE>(data.ny + 2) + 1}), KOKKOS_LAMBDA(T_indexType ix, T_indexType iy) { data.dzab(ix, iy, -2) = data.dxb(ix) * data.dyb(iy) * data.hyc(ix); });
            Kokkos::fence();
        }
        {
            // Node centred areas and volumes
            auto l = KOKKOS_LAMBDA(T_indexType ix, T_indexType iy, T_indexType iz)
            {
                T_dataType dx = data.dxb(ix);
                T_dataType dy = data.dyb(iy);
                T_dataType dz = data.dzc(iz);
                T_dataType dydz = dy * dz;
                data.dxac(ix, iy, iz) = dydz * data.hyc(ix) * data.hz2(ix, iy);
                data.dyac(ix, iy, iz) = dx * dz * data.hz1(ix, iy);
                data.dzac(ix, iy, iz) = dx * dy * data.hy(ix);
                data.cvc(ix, iy, iz) = dx * dydz * hyv(ix) * data.hz(ix, iy);
            };
            Kokkos::parallel_for("Parallel", Kokkos::MDRangePolicy<KOKKOS_EXECUTION_SPACE, Kokkos::Rank<3, LARE_ITERATE, LARE_ITERATE>, Kokkos::IndexType<SIGNED_INDEX_TYPE>>({static_cast<SIGNED_INDEX_TYPE>(-1), static_cast<SIGNED_INDEX_TYPE>(-1), static_cast<SIGNED_INDEX_TYPE>(-1)}, {static_cast<SIGNED_INDEX_TYPE>(data.nx + 2) + 1, static_cast<SIGNED_INDEX_TYPE>(data.ny + 2) + 1, static_cast<SIGNED_INDEX_TYPE>(data.nz + 2) + 1}), l);
            Kokkos::fence();
        }
        // Set up the cartesian coordinates array
        if (data.geometry == geometryType::Cartesian)
        {
            auto l = KOKKOS_LAMBDA(T_indexType ix, T_indexType iy, T_indexType iz)
            {
                data.x(ix, iy, iz) = data.xc(ix);
                data.y(ix, iy, iz) = data.yc(iy);
                data.z(ix, iy, iz) = data.zc(iz);
            };
            Kokkos::parallel_for("Parallel", Kokkos::MDRangePolicy<KOKKOS_EXECUTION_SPACE, Kokkos::Rank<3, LARE_ITERATE, LARE_ITERATE>, Kokkos::IndexType<SIGNED_INDEX_TYPE>>({static_cast<SIGNED_INDEX_TYPE>(-1), static_cast<SIGNED_INDEX_TYPE>(-1), static_cast<SIGNED_INDEX_TYPE>(-1)}, {static_cast<SIGNED_INDEX_TYPE>(data.nx + 2) + 1, static_cast<SIGNED_INDEX_TYPE>(data.ny + 2) + 1, static_cast<SIGNED_INDEX_TYPE>(data.nz + 2) + 1}), l);
            Kokkos::fence();
        }
        else if (data.geometry == geometryType::Cylindrical)
        {
            auto l = KOKKOS_LAMBDA(T_indexType ix, T_indexType iy, T_indexType iz)
            {
                T_dataType r = data.xb(ix);
                T_dataType theta = data.yb(iy);
                T_dataType zz = data.zb(iz);
                data.x(ix, iy, iz) = r;
                data.y(ix, iy, iz) = theta;
                data.z(ix, iy, iz) = zz;
                data.xp(ix, iy, iz) = r * std::cos(theta);
                data.yp(ix, iy, iz) = r * std::sin(theta);
                data.zp(ix, iy, iz) = zz;
            };
            Kokkos::parallel_for("Parallel", Kokkos::MDRangePolicy<KOKKOS_EXECUTION_SPACE, Kokkos::Rank<3, LARE_ITERATE, LARE_ITERATE>, Kokkos::IndexType<SIGNED_INDEX_TYPE>>({static_cast<SIGNED_INDEX_TYPE>(-2), static_cast<SIGNED_INDEX_TYPE>(-2), static_cast<SIGNED_INDEX_TYPE>(-2)}, {static_cast<SIGNED_INDEX_TYPE>(data.nx + 2) + 1, static_cast<SIGNED_INDEX_TYPE>(data.ny + 2) + 1, static_cast<SIGNED_INDEX_TYPE>(data.nz + 2) + 1}), l);
            Kokkos::fence();
        }
        else if (data.geometry == geometryType::Spherical)
        {
            auto l = KOKKOS_LAMBDA(T_indexType ix, T_indexType iy, T_indexType iz)
            {
                T_dataType r = data.xb(ix);
                T_dataType theta = data.yb(iy);
                T_dataType phi = data.zb(iz);
                data.x(ix, iy, iz) = r;
                data.y(ix, iy, iz) = theta;
                data.z(ix, iy, iz) = phi;
                data.xp(ix, iy, iz) = r * std::sin(theta) * std::cos(phi);
                data.yp(ix, iy, iz) = r * std::sin(theta) * std::sin(phi);
                data.zp(ix, iy, iz) = r * std::cos(theta);
            };
            Kokkos::parallel_for("Parallel", Kokkos::MDRangePolicy<KOKKOS_EXECUTION_SPACE, Kokkos::Rank<3, LARE_ITERATE, LARE_ITERATE>, Kokkos::IndexType<SIGNED_INDEX_TYPE>>({static_cast<SIGNED_INDEX_TYPE>(-2), static_cast<SIGNED_INDEX_TYPE>(-2), static_cast<SIGNED_INDEX_TYPE>(-2)}, {static_cast<SIGNED_INDEX_TYPE>(data.nx + 2) + 1, static_cast<SIGNED_INDEX_TYPE>(data.ny + 2) + 1, static_cast<SIGNED_INDEX_TYPE>(data.nz + 2) + 1}), l);
            Kokkos::fence();
        }
    }
}

#endif // USE_KOKKOS