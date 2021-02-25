//---- GPL ---------------------------------------------------------------------
//
// Copyright (C)  Stichting Deltares, 2011-2021.
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation version 3.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
// contact: delft3d.support@deltares.nl
// Stichting Deltares
// P.O. Box 177
// 2600 MH Delft, The Netherlands
//
// All indications and logos of, and references to, "Delft3D" and "Deltares"
// are registered trademarks of Stichting Deltares, and remain the property of
// Stichting Deltares. All rights reserved.
//
//------------------------------------------------------------------------------

#pragma once

#include <memory>

#include <MeshKernel/Contacts.hpp>
#include <MeshKernel/CurvilinearGridFromSplines.hpp>
#include <MeshKernel/Mesh1D.hpp>
#include <MeshKernel/Mesh2D.hpp>
#include <MeshKernel/OrthogonalizationAndSmoothing.hpp>

namespace meshkernelapi
{
    /// @brief The class holding the state of the C API library
    class MeshKernelState
    {
    public:
        /// @brief Constructor initializing mesh and contacts classes
        MeshKernelState()
        {
            m_mesh1d = std::make_shared<meshkernel::Mesh1D>();
            m_mesh2d = std::make_shared<meshkernel::Mesh2D>();
            m_contacts = std::make_shared<meshkernel::Contacts>(m_mesh1d, m_mesh2d);
        }

        // Mesh instances
        std::shared_ptr<meshkernel::Mesh1D> m_mesh1d;     ///< Shared pointer to meshkernel::Mesh1D instance
        std::shared_ptr<meshkernel::Mesh2D> m_mesh2d;     ///< Shared pointer to meshkernel::Mesh2D instance
        std::shared_ptr<meshkernel::Contacts> m_contacts; ///< Shared pointer to meshkernel::Contacts instance

        // For interactivity
        std::shared_ptr<meshkernel::OrthogonalizationAndSmoothing> m_orthogonalization;       ///< Shared pointer to meshkernel::OrthogonalizationAndSmoothing instance
        std::shared_ptr<meshkernel::CurvilinearGridFromSplines> m_curvilinearGridFromSplines; ///< Shared pointer to meshkernel::CurvilinearGridFromSplines instance
    };

} // namespace meshkernelapi
