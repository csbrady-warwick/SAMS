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
#include "remapData.h"

namespace LARE
{

    namespace pw = portableWrapper;

    void LARE3DKokkos::eulerian_remap(simulationData &data, remapData &remap_data)
    {
        using Range = pw::Range;
        int case_test;
        pw::portableArrayManager remapManager;

        // We can allocate everything other than flux here
        remapManager.allocateKOV(remap_data.rho1, Range(-1, data.nx + 2), Range(-1, data.ny + 2), Range(-1, data.nz + 2));
        remapManager.allocateKOV(remap_data.cv2, Range(-1, data.nx + 2), Range(-1, data.ny + 2), Range(-1, data.nz + 2));
        remapManager.allocateKOV(remap_data.cvc1, Range(-1, data.nx + 2), Range(-1, data.ny + 2), Range(-1, data.nz + 2));
        remapManager.allocateKOV(remap_data.db1, Range(-1, data.nx + 2), Range(-1, data.ny + 2), Range(-1, data.nz + 2));
        remapManager.allocateKOV(remap_data.rho_v, Range(-1, data.nx + 2), Range(-1, data.ny + 2), Range(-1, data.nz + 2));
        remapManager.allocateKOV(remap_data.rho_v1, Range(-1, data.nx + 2), Range(-1, data.ny + 2), Range(-1, data.nz + 2));
        // Flux is one element larger in the direction of remap, so it is allocated in each remap function

        if (data.rke)
        {
            Kokkos::deep_copy(data.delta_ke, 0.0);
        }
        remap_data.xpass = 1.0;
        remap_data.ypass = 1.0;
        remap_data.zpass = 1.0;

        Kokkos::parallel_for("Parallel", Kokkos::MDRangePolicy<KOKKOS_EXECUTION_SPACE, Kokkos::Rank<3, LARE_ITERATE, LARE_ITERATE>, Kokkos::IndexType<SIGNED_INDEX_TYPE>>({static_cast<SIGNED_INDEX_TYPE>(-2), static_cast<SIGNED_INDEX_TYPE>(-1), static_cast<SIGNED_INDEX_TYPE>(-1)}, {static_cast<SIGNED_INDEX_TYPE>(data.nx + 2) + 1, static_cast<SIGNED_INDEX_TYPE>(data.ny + 2) + 1, static_cast<SIGNED_INDEX_TYPE>(data.nz + 2) + 1}), KOKKOS_LAMBDA(T_indexType ix, T_indexType iy, T_indexType iz) { data.bx(ix, iy, iz) *= data.dxab(ix, iy, iz); });

        Kokkos::parallel_for("Parallel", Kokkos::MDRangePolicy<KOKKOS_EXECUTION_SPACE, Kokkos::Rank<3, LARE_ITERATE, LARE_ITERATE>, Kokkos::IndexType<SIGNED_INDEX_TYPE>>({static_cast<SIGNED_INDEX_TYPE>(-1), static_cast<SIGNED_INDEX_TYPE>(-2), static_cast<SIGNED_INDEX_TYPE>(-1)}, {static_cast<SIGNED_INDEX_TYPE>(data.nx + 2) + 1, static_cast<SIGNED_INDEX_TYPE>(data.ny + 2) + 1, static_cast<SIGNED_INDEX_TYPE>(data.nz + 2) + 1}), KOKKOS_LAMBDA(T_indexType ix, T_indexType iy, T_indexType iz) { data.by(ix, iy, iz) *= data.dyab(ix, iy, iz); });

        Kokkos::parallel_for("Parallel", Kokkos::MDRangePolicy<KOKKOS_EXECUTION_SPACE, Kokkos::Rank<3, LARE_ITERATE, LARE_ITERATE>, Kokkos::IndexType<SIGNED_INDEX_TYPE>>({static_cast<SIGNED_INDEX_TYPE>(-1), static_cast<SIGNED_INDEX_TYPE>(-1), static_cast<SIGNED_INDEX_TYPE>(-2)}, {static_cast<SIGNED_INDEX_TYPE>(data.nx + 2) + 1, static_cast<SIGNED_INDEX_TYPE>(data.ny + 2) + 1, static_cast<SIGNED_INDEX_TYPE>(data.nz + 2) + 1}), KOKKOS_LAMBDA(T_indexType ix, T_indexType iy, T_indexType iz) { data.bz(ix, iy, iz) *= data.dzab(ix, iy, iz); });

        case_test = data.step % 6;

        // Strang ordering
        switch (case_test)
        {
        case 0:
            remap_x(data, remap_data);
            remap_y(data, remap_data);
            remap_z(data, remap_data);
            break;
        case 1:
            remap_y(data, remap_data);
            remap_z(data, remap_data);
            remap_x(data, remap_data);
            break;
        case 2:
            remap_z(data, remap_data);
            remap_x(data, remap_data);
            remap_y(data, remap_data);
            break;
        case 3:
            remap_x(data, remap_data);
            remap_z(data, remap_data);
            remap_y(data, remap_data);
            break;
        case 4:
            remap_z(data, remap_data);
            remap_y(data, remap_data);
            remap_x(data, remap_data);
            break;
        case 5:
            remap_y(data, remap_data);
            remap_x(data, remap_data);
            remap_z(data, remap_data);
            break;
        }

        // bx = bx / (dxab + none_zero), etc.

        Kokkos::parallel_for("Parallel", Kokkos::MDRangePolicy<KOKKOS_EXECUTION_SPACE, Kokkos::Rank<3, LARE_ITERATE, LARE_ITERATE>, Kokkos::IndexType<SIGNED_INDEX_TYPE>>({static_cast<SIGNED_INDEX_TYPE>(-2), static_cast<SIGNED_INDEX_TYPE>(-1), static_cast<SIGNED_INDEX_TYPE>(-1)}, {static_cast<SIGNED_INDEX_TYPE>(data.nx + 2) + 1, static_cast<SIGNED_INDEX_TYPE>(data.ny + 2) + 1, static_cast<SIGNED_INDEX_TYPE>(data.nz + 2) + 1}), KOKKOS_LAMBDA(T_indexType ix, T_indexType iy, T_indexType iz) { data.bx(ix, iy, iz) /= (data.dxab(ix, iy, iz) + data.none_zero); });

        Kokkos::parallel_for("Parallel", Kokkos::MDRangePolicy<KOKKOS_EXECUTION_SPACE, Kokkos::Rank<3, LARE_ITERATE, LARE_ITERATE>, Kokkos::IndexType<SIGNED_INDEX_TYPE>>({static_cast<SIGNED_INDEX_TYPE>(-1), static_cast<SIGNED_INDEX_TYPE>(-2), static_cast<SIGNED_INDEX_TYPE>(-1)}, {static_cast<SIGNED_INDEX_TYPE>(data.nx + 2) + 1, static_cast<SIGNED_INDEX_TYPE>(data.ny + 2) + 1, static_cast<SIGNED_INDEX_TYPE>(data.nz + 2) + 1}), KOKKOS_LAMBDA(T_indexType ix, T_indexType iy, T_indexType iz) { data.by(ix, iy, iz) /= (data.dyab(ix, iy, iz) + data.none_zero); });

        Kokkos::parallel_for("Parallel", Kokkos::MDRangePolicy<KOKKOS_EXECUTION_SPACE, Kokkos::Rank<3, LARE_ITERATE, LARE_ITERATE>, Kokkos::IndexType<SIGNED_INDEX_TYPE>>({static_cast<SIGNED_INDEX_TYPE>(-1), static_cast<SIGNED_INDEX_TYPE>(-1), static_cast<SIGNED_INDEX_TYPE>(-2)}, {static_cast<SIGNED_INDEX_TYPE>(data.nx + 2) + 1, static_cast<SIGNED_INDEX_TYPE>(data.ny + 2) + 1, static_cast<SIGNED_INDEX_TYPE>(data.nz + 2) + 1}), KOKKOS_LAMBDA(T_indexType ix, T_indexType iy, T_indexType iz) { data.bz(ix, iy, iz) /= (data.dzab(ix, iy, iz) + data.none_zero); });

        bfield_bcs();

        // Set the grid positions back to their default value
        Kokkos::parallel_for("Parallel", Kokkos::MDRangePolicy<KOKKOS_EXECUTION_SPACE, Kokkos::Rank<3, LARE_ITERATE, LARE_ITERATE>, Kokkos::IndexType<SIGNED_INDEX_TYPE>>({static_cast<SIGNED_INDEX_TYPE>(-2), static_cast<SIGNED_INDEX_TYPE>(-1), static_cast<SIGNED_INDEX_TYPE>(-1)}, {static_cast<SIGNED_INDEX_TYPE>(data.nx + 2) + 1, static_cast<SIGNED_INDEX_TYPE>(data.ny + 2) + 1, static_cast<SIGNED_INDEX_TYPE>(data.nz + 2) + 1}), KOKKOS_LAMBDA(T_indexType ix, T_indexType iy, T_indexType iz) {
        data.x(ix, iy, iz) = data.xb(ix);
        data.y(ix, iy, iz) = data.yb(iy);
        data.z(ix, iy, iz) = data.zb(iz); });

        remapManager.deallocateKOV(remap_data.rho1);
        remapManager.deallocateKOV(remap_data.cv2);
        remapManager.deallocateKOV(remap_data.cvc1);
        remapManager.deallocateKOV(remap_data.db1);
        remapManager.deallocateKOV(remap_data.rho_v);
        remapManager.deallocateKOV(remap_data.rho_v1);
    }
}

#endif // USE_KOKKOS