//---- GPL ---------------------------------------------------------------------
//
// Copyright (C)  Stichting Deltares, 2011-2020.
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
#include <vector>

#include <MeshKernel/Entities.hpp>
#include <MeshKernel/Mesh.hpp>
#include <MeshKernel/Network1D.hpp>

/// \namespace meshkernel
/// @brief Contains the logic of the C++ static library
namespace meshkernel
{
    /// @brief A class derived from Mesh, which describes 1d meshes.
    ///
    /// A 1d mesh is composed of a series of connected edges
    /// representing 1d real word features, such as pipes or a sewage network.
    class Mesh1D : public Mesh
    {

    public:
        /// @brief Default constructor
        Mesh1D() = default;

        /// @brief Construct a mesh1d starting from the edges and nodes
        /// @param[in] edges The input edges
        /// @param[in] nodes The input nodes
        /// @param[in] projection  The projection to use
        explicit Mesh1D(std::vector<Edge> const& edges,
                        std::vector<Point> const& nodes,
                        Projection projection);

        /// @brief Constructs a mesh 1d from a network 1d. The network contains the chainages where the discratization points will be computed.
        /// @param[in] Network1D The input network
        /// @param[in] minFaceSize The minimum face size below which two nodes will be merged
        explicit Mesh1D(Network1D& network1d, double minFaceSize);

        /// @brief Inquire if a mesh 1d-node is on boundary
        /// @param[in] node The node index
        /// @return If the node is on boundary
        [[nodiscard]] bool IsNodeOnBoundary(size_t node) const { return m_nodesNumEdges[node] == 1; }
    };

} // namespace meshkernel
