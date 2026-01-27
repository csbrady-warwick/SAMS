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
#include "shared_data.h"
#include "variableRegistry.h"

#define DISABLE // DISABLE

namespace LARE
{

    namespace pw = portableWrapper;

    void LARE3D::boundary_conditions(simulationData &data)
    {
        DISABLE
        bfield_bcs(data);
        energy_bcs(data);
        density_bcs(data);
        velocity_bcs(data);
        pw::fence(); // Ensure all operations are complete before returning
    }

    void LARE3D::bfield_bcs(simulationData &data)
    {
        DISABLE
        harness.variableRegistry.haloExchange("bx");
        harness.variableRegistry.haloExchange("by");
        harness.variableRegistry.haloExchange("bz");
        if (data.xbc_min == BCType::BC_OTHER && data.isxLB)
        {
            pw::assign(
                data.bx(-2, pw::Range(), pw::Range()),
                data.bx(2, pw::Range(), pw::Range()));
            pw::assign(
                data.bx(pw::Range(-1, -1), pw::Range(), pw::Range()),
                data.bx(pw::Range(1, 1), pw::Range(), pw::Range()));
            pw::assign(
                data.by(pw::Range(-1, -1), pw::Range(), pw::Range()),
                data.by(pw::Range(2, 2), pw::Range(), pw::Range()));
            pw::assign(
                data.by(pw::Range(0, 0), pw::Range(), pw::Range()),
                data.by(pw::Range(1, 1), pw::Range(), pw::Range()));
            pw::assign(
                data.bz(pw::Range(-1, -1), pw::Range(), pw::Range()),
                data.bz(pw::Range(2, 2), pw::Range(), pw::Range()));
            pw::assign(
                data.bz(pw::Range(0, 0), pw::Range(), pw::Range()),
                data.bz(pw::Range(1, 1), pw::Range(), pw::Range()));
        }

        if (data.xbc_max == BCType::BC_OTHER && data.isxUB)
        {
            pw::assign(
                data.bx(data.nx + 1, pw::Range(), pw::Range()),
                data.bx(data.nx - 1, pw::Range(), pw::Range()));
            pw::assign(
                data.bx(data.nx + 2, pw::Range(), pw::Range()),
                data.bx(data.nx - 2, pw::Range(), pw::Range()));
            pw::assign(
                data.by(data.nx + 1, pw::Range(), pw::Range()),
                data.by(data.nx, pw::Range(), pw::Range()));
            pw::assign(
                data.by(data.nx + 2, pw::Range(), pw::Range()),
                data.by(data.nx - 1, pw::Range(), pw::Range()));
            pw::assign(
                data.bz(data.nx + 1, pw::Range(), pw::Range()),
                data.bz(data.nx, pw::Range(), pw::Range()));
            pw::assign(
                data.bz(data.nx + 2, pw::Range(), pw::Range()),
                data.bz(data.nx - 1, pw::Range(), pw::Range()));
        }

        if (data.ybc_min == BCType::BC_OTHER && data.isyLB)
        {
            pw::assign(
                data.by(pw::Range(), -2, pw::Range()),
                data.by(pw::Range(), 2, pw::Range()));
            pw::assign(
                data.by(pw::Range(), -1, pw::Range()),
                data.by(pw::Range(), 1, pw::Range()));
            pw::assign(
                data.bx(pw::Range(), -1, pw::Range()),
                data.bx(pw::Range(), 2, pw::Range()));
            pw::assign(
                data.bx(pw::Range(), 0, pw::Range()),
                data.bx(pw::Range(), 1, pw::Range()));
            pw::assign(
                data.bz(pw::Range(), -1, pw::Range()),
                data.bz(pw::Range(), 2, pw::Range()));
            pw::assign(
                data.bz(pw::Range(), 0, pw::Range()),
                data.bz(pw::Range(), 1, pw::Range()));
        }

        if (data.ybc_max == BCType::BC_OTHER && data.isyUB)
        {
            pw::assign(
                data.by(pw::Range(), data.ny + 1, pw::Range()),
                data.by(pw::Range(), data.ny - 1, pw::Range()));
            pw::assign(
                data.by(pw::Range(), data.ny + 2, pw::Range()),
                data.by(pw::Range(), data.ny - 2, pw::Range()));
            pw::assign(
                data.bx(pw::Range(), data.ny + 1, pw::Range()),
                data.bx(pw::Range(), data.ny, pw::Range()));
            pw::assign(
                data.bx(pw::Range(), data.ny + 2, pw::Range()),
                data.bx(pw::Range(), data.ny - 1, pw::Range()));
            pw::assign(
                data.bz(pw::Range(), data.ny + 1, pw::Range()),
                data.bz(pw::Range(), data.ny, pw::Range()));
            pw::assign(
                data.bz(pw::Range(), data.ny + 2, pw::Range()),
                data.bz(pw::Range(), data.ny - 1, pw::Range()));
        }

        if (data.zbc_min == BCType::BC_OTHER && data.iszLB)
        {
            pw::assign(
                data.bz(pw::Range(), pw::Range(), -2),
                data.bz(pw::Range(), pw::Range(), 2));
            pw::assign(
                data.bz(pw::Range(), pw::Range(), -1),
                data.bz(pw::Range(), pw::Range(), 1));
            pw::assign(
                data.by(pw::Range(), pw::Range(), -1),
                data.by(pw::Range(), pw::Range(), 2));
            pw::assign(
                data.by(pw::Range(), pw::Range(), 0),
                data.by(pw::Range(), pw::Range(), 1));
            pw::assign(
                data.bx(pw::Range(), pw::Range(), -1),
                data.bx(pw::Range(), pw::Range(), 2));
            pw::assign(
                data.bx(pw::Range(), pw::Range(), 0),
                data.bx(pw::Range(), pw::Range(), 1));
        }

        if (data.zbc_max == BCType::BC_OTHER && data.iszUB)
        {
            pw::assign(
                data.bz(pw::Range(), pw::Range(), data.nz + 1),
                data.bz(pw::Range(), pw::Range(), data.nz - 1));
            pw::assign(
                data.bz(pw::Range(), pw::Range(), data.nz + 2),
                data.bz(pw::Range(), pw::Range(), data.nz - 2));
            pw::assign(
                data.by(pw::Range(), pw::Range(), data.nz + 1),
                data.by(pw::Range(), pw::Range(), data.nz));
            pw::assign(
                data.by(pw::Range(), pw::Range(), data.nz + 2),
                data.by(pw::Range(), pw::Range(), data.nz - 1));
            pw::assign(
                data.bx(pw::Range(), pw::Range(), data.nz + 1),
                data.bx(pw::Range(), pw::Range(), data.nz));
            pw::assign(
                data.bx(pw::Range(), pw::Range(), data.nz + 2),
                data.bx(pw::Range(), pw::Range(), data.nz - 1));
        }
    }

    void LARE3D::energy_bcs(simulationData &data)
    {
        DISABLE
        harness.variableRegistry.haloExchange("energy_electron");
        harness.variableRegistry.haloExchange("energy_ion");
        if (data.xbc_min == BCType::BC_OTHER && data.isxLB)
        {
            pw::assign(
                data.energy_electron(-1, pw::Range(), pw::Range()),
                data.energy_electron(2, pw::Range(), pw::Range()));
            pw::assign(
                data.energy_electron(0, pw::Range(), pw::Range()),
                data.energy_electron(1, pw::Range(), pw::Range()));
            pw::assign(
                data.energy_ion(-1, pw::Range(), pw::Range()),
                data.energy_ion(2, pw::Range(), pw::Range()));
            pw::assign(
                data.energy_ion(0, pw::Range(), pw::Range()),
                data.energy_ion(1, pw::Range(), pw::Range()));
        }

        if (data.xbc_max == BCType::BC_OTHER && data.isxUB)
        {
            pw::assign(
                data.energy_electron(data.nx + 1, pw::Range(), pw::Range()),
                data.energy_electron(data.nx, pw::Range(), pw::Range()));
            pw::assign(
                data.energy_electron(data.nx + 2, pw::Range(), pw::Range()),
                data.energy_electron(data.nx - 1, pw::Range(), pw::Range()));
            pw::assign(
                data.energy_ion(data.nx + 1, pw::Range(), pw::Range()),
                data.energy_ion(data.nx, pw::Range(), pw::Range()));
            pw::assign(
                data.energy_ion(data.nx + 2, pw::Range(), pw::Range()),
                data.energy_ion(data.nx - 1, pw::Range(), pw::Range()));
        }

        if (data.ybc_min == BCType::BC_OTHER && data.isyLB)
        {
            pw::assign(
                data.energy_electron(pw::Range(), -1, pw::Range()),
                data.energy_electron(pw::Range(), 2, pw::Range()));
            pw::assign(
                data.energy_electron(pw::Range(), 0, pw::Range()),
                data.energy_electron(pw::Range(), 1, pw::Range()));
            pw::assign(
                data.energy_ion(pw::Range(), -1, pw::Range()),
                data.energy_ion(pw::Range(), 2, pw::Range()));
            pw::assign(
                data.energy_ion(pw::Range(), 0, pw::Range()),
                data.energy_ion(pw::Range(), 1, pw::Range()));
        }

        if (data.ybc_max == BCType::BC_OTHER && data.isyUB)
        {
            pw::assign(
                data.energy_electron(pw::Range(), data.ny + 1, pw::Range()),
                data.energy_electron(pw::Range(), data.ny, pw::Range()));
            pw::assign(
                data.energy_electron(pw::Range(), data.ny + 2, pw::Range()),
                data.energy_electron(pw::Range(), data.ny - 1, pw::Range()));
            pw::assign(
                data.energy_ion(pw::Range(), data.ny + 1, pw::Range()),
                data.energy_ion(pw::Range(), data.ny, pw::Range()));
            pw::assign(
                data.energy_ion(pw::Range(), data.ny + 2, pw::Range()),
                data.energy_ion(pw::Range(), data.ny - 1, pw::Range()));
        }
        if (data.zbc_min == BCType::BC_OTHER && data.iszLB)
        {
            pw::assign(
                data.energy_electron(pw::Range(), pw::Range(), -1),
                data.energy_electron(pw::Range(), pw::Range(), 2));
            pw::assign(
                data.energy_electron(pw::Range(), pw::Range(), 0),
                data.energy_electron(pw::Range(), pw::Range(), 1));
            pw::assign(
                data.energy_ion(pw::Range(), pw::Range(), -1),
                data.energy_ion(pw::Range(), pw::Range(), 2));
            pw::assign(
                data.energy_ion(pw::Range(), pw::Range(), 0),
                data.energy_ion(pw::Range(), pw::Range(), 1));
        }
        if (data.zbc_max == BCType::BC_OTHER && data.iszUB)
        {
            pw::assign(
                data.energy_electron(pw::Range(), pw::Range(), data.nz + 1),
                data.energy_electron(pw::Range(), pw::Range(), data.nz));
            pw::assign(
                data.energy_electron(pw::Range(), pw::Range(), data.nz + 2),
                data.energy_electron(pw::Range(), pw::Range(), data.nz - 1));
            pw::assign(
                data.energy_ion(pw::Range(), pw::Range(), data.nz + 1),
                data.energy_ion(pw::Range(), pw::Range(), data.nz));
            pw::assign(
                data.energy_ion(pw::Range(), pw::Range(), data.nz + 2),
                data.energy_ion(pw::Range(), pw::Range(), data.nz - 1));
        }
    }

    void LARE3D::density_bcs(simulationData &data)
    {
        DISABLE
        harness.variableRegistry.haloExchange("rho");
        if (data.xbc_min == BCType::BC_OTHER && data.isxLB)
        {
            pw::assign(
                data.rho(-1, pw::Range(), pw::Range()),
                data.rho(2, pw::Range(), pw::Range()));
            pw::assign(
                data.rho(0, pw::Range(), pw::Range()),
                data.rho(1, pw::Range(), pw::Range()));
        }

        if (data.xbc_max == BCType::BC_OTHER && data.isxUB)
        {
            pw::assign(
                data.rho(data.nx + 1, pw::Range(), pw::Range()),
                data.rho(data.nx - 1, pw::Range(), pw::Range()));
            pw::assign(
                data.rho(data.nx + 2, pw::Range(), pw::Range()),
                data.rho(data.nx - 2, pw::Range(), pw::Range()));
        }

        if (data.ybc_min == BCType::BC_OTHER && data.isyLB)
        {
            pw::assign(
                data.rho(pw::Range(), -1, pw::Range()),
                data.rho(pw::Range(), 2, pw::Range()));
            pw::assign(
                data.rho(pw::Range(), 0, pw::Range()),
                data.rho(pw::Range(), 1, pw::Range()));
        }

        if (data.ybc_max == BCType::BC_OTHER && data.isyUB)
        {
            pw::assign(
                data.rho(pw::Range(), data.ny + 1, pw::Range()),
                data.rho(pw::Range(), data.ny - 1, pw::Range()));
            pw::assign(
                data.rho(pw::Range(), data.ny + 2, pw::Range()),
                data.rho(pw::Range(), data.ny - 2, pw::Range()));
        }

        if (data.zbc_min == BCType::BC_OTHER && data.iszLB)
        {
            pw::assign(
                data.rho(pw::Range(), pw::Range(), -1),
                data.rho(pw::Range(), pw::Range(), 2));
            pw::assign(
                data.rho(pw::Range(), pw::Range(), 0),
                data.rho(pw::Range(), pw::Range(), 1));
        }
    }

    void LARE3D::velocity_bcs(simulationData &data)
    {
        DISABLE
        harness.variableRegistry.haloExchange("vx");
        harness.variableRegistry.haloExchange("vy");
        harness.variableRegistry.haloExchange("vz");
        // Other boundaries clamp v=0
        if (data.xbc_min == BCType::BC_OTHER && data.isxLB)
        {
            pw::assign(
                data.vx(pw::Range(-2, 0), pw::Range(), pw::Range()),
                0.0);
            pw::assign(
                data.vy(-1, pw::Range(), pw::Range()),
                data.vy(1, pw::Range(), pw::Range()));
            pw::assign(
                data.vy(-2, pw::Range(), pw::Range()),
                data.vy(2, pw::Range(), pw::Range()));

            pw::assign(
                data.vz(-1, pw::Range(), pw::Range()),
                data.vz(1, pw::Range(), pw::Range()));
            pw::assign(
                data.vz(-2, pw::Range(), pw::Range()),
                data.vz(2, pw::Range(), pw::Range()));
        }

        if (data.xbc_max == BCType::BC_OTHER && data.isxUB)
        {
            pw::assign(
                data.vx(pw::Range(data.nx, data.nx + 2), pw::Range(), pw::Range()),
                0.0);

            pw::assign(
                data.vy(data.nx + 1, pw::Range(), pw::Range()),
                data.vy(data.nx - 1, pw::Range(), pw::Range()));
            pw::assign(
                data.vy(data.nx + 2, pw::Range(), pw::Range()),
                data.vy(data.nx - 2, pw::Range(), pw::Range()));
            pw::assign(
                data.vz(data.nx + 1, pw::Range(), pw::Range()),
                data.vz(data.nx, pw::Range(), pw::Range()));
            pw::assign(
                data.vz(data.nx + 2, pw::Range(), pw::Range()),
                data.vz(data.nx - 1, pw::Range(), pw::Range()));
        }

        if (data.ybc_min == BCType::BC_OTHER && data.isyLB)
        {
            pw::assign(
                data.vx(pw::Range(), -1, pw::Range()),
                data.vx(pw::Range(), 1, pw::Range()));
            pw::assign(
                data.vx(pw::Range(), -2, pw::Range()),
                data.vx(pw::Range(), 2, pw::Range()));
            pw::assign(
                data.vy(pw::Range(), pw::Range(-2, 0), pw::Range()),
                0.0);
            pw::assign(
                data.vz(pw::Range(), -1, pw::Range()),
                data.vz(pw::Range(), 1, pw::Range()));
            pw::assign(
                data.vz(pw::Range(), -2, pw::Range()),
                data.vz(pw::Range(), 2, pw::Range()));
        }

        if (data.ybc_max == BCType::BC_OTHER && data.isyUB)
        {
            pw::assign(
                data.vx(pw::Range(), data.ny + 1, pw::Range()),
                data.vx(pw::Range(), data.ny - 1, pw::Range()));
            pw::assign(
                data.vx(pw::Range(), data.ny + 2, pw::Range()),
                data.vx(pw::Range(), data.ny - 2, pw::Range()));
            pw::assign(
                data.vy(pw::Range(), pw::Range(data.ny, data.ny + 2), pw::Range()),
                0.0);
            pw::assign(
                data.vz(pw::Range(), data.ny + 1, pw::Range()),
                data.vz(pw::Range(), data.ny, pw::Range()));
            pw::assign(
                data.vz(pw::Range(), data.ny + 2, pw::Range()),
                data.vz(pw::Range(), data.ny - 1, pw::Range()));
        }

        if (data.zbc_min == BCType::BC_OTHER && data.iszLB)
        {
            pw::assign(
                data.vx(pw::Range(), pw::Range(), -2),
                data.vx(pw::Range(), pw::Range(), 2));
            pw::assign(
                data.vx(pw::Range(), pw::Range(), -1),
                data.vx(pw::Range(), pw::Range(), 1));
            pw::assign(
                data.vy(pw::Range(), pw::Range(), -1),
                data.vy(pw::Range(), pw::Range(), 1));
            pw::assign(
                data.vy(pw::Range(), pw::Range(), -2),
                data.vy(pw::Range(), pw::Range(), 2));
            pw::assign(
                data.vz(pw::Range(), pw::Range(), pw::Range(-2, 0)),
                0.0);
        }

        if (data.zbc_max == BCType::BC_OTHER && data.iszUB)
        {
            pw::assign(
                data.vx(pw::Range(), pw::Range(), data.nz + 1),
                data.vx(pw::Range(), pw::Range(), data.nz - 1));
            pw::assign(
                data.vy(pw::Range(), pw::Range(), data.nz + 1),
                data.vy(pw::Range(), pw::Range(), data.nz - 1));
            pw::assign(
                data.vx(pw::Range(), pw::Range(), data.nz + 2),
                data.vx(pw::Range(), pw::Range(), data.nz - 2));
            pw::assign(
                data.vy(pw::Range(), pw::Range(), data.nz + 2),
                data.vy(pw::Range(), pw::Range(), data.nz - 2));
            pw::assign(
                data.vz(pw::Range(), pw::Range(), pw::Range(data.nz, data.nz + 2)),
                0.0);
        }
    }

    void LARE3D::remap_v_bcs(simulationData &data)
    {
        DISABLE
        harness.variableRegistry.haloExchange("LARE/vx1");
        harness.variableRegistry.haloExchange("LARE/vy1");
        harness.variableRegistry.haloExchange("LARE/vz1");
        // Other boundaries clamp v=0
        if (data.xbc_min == BCType::BC_OTHER && data.isxLB)
        {
            pw::assign(
                data.vx1(pw::Range(-2, 0), pw::Range(), pw::Range()),
                0.0);

            pw::assign(
                data.vy1(-2, pw::Range(), pw::Range()),
                data.vy1(2, pw::Range(), pw::Range()));

            pw::assign(
                data.vy1(-1, pw::Range(), pw::Range()),
                data.vy1(1, pw::Range(), pw::Range()));

            pw::assign(
                data.vz1(-2, pw::Range(), pw::Range()),
                data.vz1(2, pw::Range(), pw::Range()));

            pw::assign(
                data.vz1(-1, pw::Range(), pw::Range()),
                data.vz1(1, pw::Range(), pw::Range()));
        }

        if (data.xbc_max == BCType::BC_OTHER && data.isxUB)
        {
            pw::assign(
                data.vx1(pw::Range(data.nx, data.nx + 2), pw::Range(), pw::Range()),
                0.0);
            pw::assign(
                data.vy1(data.nx + 1, pw::Range(), pw::Range()),
                data.vy1(data.nx - 1, pw::Range(), pw::Range()));
            pw::assign(
                data.vy1(data.nx + 2, pw::Range(), pw::Range()),
                data.vy1(data.nx - 2, pw::Range(), pw::Range()));
            pw::assign(
                data.vz1(data.nx + 1, pw::Range(), pw::Range()),
                data.vz1(data.nx - 1, pw::Range(), pw::Range()));
            pw::assign(
                data.vz1(data.nx + 2, pw::Range(), pw::Range()),
                data.vz1(data.nx - 2, pw::Range(), pw::Range()));
        }

        if (data.ybc_min == BCType::BC_OTHER && data.isyLB)
        {
            pw::assign(
                data.vy1(pw::Range(), pw::Range(-2, 0), pw::Range()),
                0.0);
            pw::assign(
                data.vx1(pw::Range(), -2, pw::Range()),
                data.vx1(pw::Range(), 2, pw::Range()));
            pw::assign(
                data.vx1(pw::Range(), -1, pw::Range()),
                data.vx1(pw::Range(), 1, pw::Range()));
            pw::assign(
                data.vz1(pw::Range(), -2, pw::Range()),
                data.vz1(pw::Range(), 2, pw::Range()));
            pw::assign(
                data.vz1(pw::Range(), -1, pw::Range()),
                data.vz1(pw::Range(), 1, pw::Range()));
        }

        if (data.ybc_max == BCType::BC_OTHER && data.isyUB)
        {
            pw::assign(
                data.vy1(pw::Range(), pw::Range(data.ny, data.ny + 2), pw::Range()),
                0.0);
            pw::assign(
                data.vx1(pw::Range(), data.ny + 1, pw::Range()),
                data.vx1(pw::Range(), data.ny - 1, pw::Range()));
            pw::assign(
                data.vx1(pw::Range(), data.ny + 2, pw::Range()),
                data.vx1(pw::Range(), data.ny - 2, pw::Range()));
            pw::assign(
                data.vz1(pw::Range(), data.ny + 1, pw::Range()),
                data.vz1(pw::Range(), data.ny - 1, pw::Range()));
            pw::assign(
                data.vz1(pw::Range(), data.ny + 2, pw::Range()),
                data.vz1(pw::Range(), data.ny - 2, pw::Range()));
        }

        if (data.zbc_min == BCType::BC_OTHER && data.iszLB)
        {
            pw::assign(
                data.vz1(pw::Range(), pw::Range(), pw::Range(-2, 0)),
                0.0);
            pw::assign(
                data.vx1(pw::Range(), pw::Range(), -2),
                data.vx1(pw::Range(), pw::Range(), 2));
            pw::assign(
                data.vx1(pw::Range(), pw::Range(), -1),
                data.vx1(pw::Range(), pw::Range(), 1));
            pw::assign(
                data.vy1(pw::Range(), pw::Range(), -2),
                data.vy1(pw::Range(), pw::Range(), 2));
            pw::assign(
                data.vy1(pw::Range(), pw::Range(), -1),
                data.vy1(pw::Range(), pw::Range(), 1));
        }

        if (data.zbc_max == BCType::BC_OTHER && data.iszUB)
        {
            pw::assign(
                data.vz1(pw::Range(), pw::Range(), pw::Range(data.nz, data.nz + 2)),
                0.0);
            pw::assign(
                data.vx1(pw::Range(), pw::Range(), data.nz + 1),
                data.vx1(pw::Range(), pw::Range(), data.nz - 1));
            pw::assign(
                data.vx1(pw::Range(), pw::Range(), data.nz + 2),
                data.vx1(pw::Range(), pw::Range(), data.nz - 2));
            pw::assign(
                data.vy1(pw::Range(), pw::Range(), data.nz + 1),
                data.vy1(pw::Range(), pw::Range(), data.nz - 1));
            pw::assign(
                data.vy1(pw::Range(), pw::Range(), data.nz + 2),
                data.vy1(pw::Range(), pw::Range(), data.nz - 2));
        }
    }

    void LARE3D::dm_x_bcs(simulationData &data, remapData &remap_data)
    {
        DISABLE
        harness.variableRegistry.haloExchange("LARE/dm");
        // Continous flux at x-max boundary
        if (data.isxUB)
        {
            pw::assign(
                remap_data.dm(data.nx + 1, pw::Range(), pw::Range()),
                remap_data.dm(data.nx, pw::Range(), pw::Range()));
        }

        // Continous flux at x-min boundary
        if (data.isxLB)
        {
            pw::assign(
                remap_data.dm(-1, pw::Range(), pw::Range()),
                remap_data.dm(0, pw::Range(), pw::Range()));
        }
        pw::fence();
    }

    void LARE3D::dm_y_bcs(simulationData &data, remapData &remap_data)
    {
        DISABLE
        harness.variableRegistry.haloExchange("LARE/dm");
        // Continous flux at y-max boundary
        if (data.isyUB)
        {
            pw::assign(
                remap_data.dm(pw::Range(), data.ny + 1, pw::Range()),
                remap_data.dm(pw::Range(), data.ny, pw::Range()));
        }
        // Continous flux at y-min boundary
        if (data.isyLB)
        {
            pw::assign(
                remap_data.dm(pw::Range(), -1, pw::Range()),
                remap_data.dm(pw::Range(), 0, pw::Range()));
        }

        pw::fence();
    }

    void LARE3D::dm_z_bcs(simulationData &data, remapData &remap_data)
    {
        DISABLE
        harness.variableRegistry.haloExchange("LARE/dm");
        // Continous flux at z-max boundary
        if (data.iszUB)
        {
            pw::assign(
                remap_data.dm(pw::Range(), pw::Range(), data.nz + 1),
                remap_data.dm(pw::Range(), pw::Range(), data.nz));
        }

        // Continous flux at z-min boundary
        if (data.iszLB)
        {
            pw::assign(
                remap_data.dm(pw::Range(), pw::Range(), -1),
                remap_data.dm(pw::Range(), pw::Range(), 0));
        }
        pw::fence();
    }
}
