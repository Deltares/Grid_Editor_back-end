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

#include <vector>

#include <MeshKernel/Entities.hpp>

namespace meshkernel
{
    /// @brief A class used for triangulation interpolation
    ///
    /// As for averaging, the triangle interpolation operates at three specific
    /// \ref MeshLocations - Faces, Nodes, and Edges.
    /// The idea is to triangulate the samples and identify for each location
    /// the triangle that fully contains
    /// Only the values at the nodes of the identified triangle are used in
    /// the computation of each location.
    /// The algorithm operates as follow:
    ///
    /// - The triangulation of the samples is computed
    ///
    /// - For each triangle, the circumcentre is computed
    ///
    /// - The triangle circumcentres are ordered in an RTree for a fast search
    ///
    /// - For each location, the closest circumcentre is found
    ///
    /// - If the corresponding triangle contains the location then the linear
    ///   interpolation is performed, otherwise, the next neighbouring
    ///   triangle is searched. The next neighbouring triangle is the first
    ///   triangle that satisfies these two conditions:
    ///
    ///   1.  shares one of the current triangle edges.
    ///
    ///   2.  the crossing of the edge and the line connecting the location
    ///       the current triangle circumcentre exists.
    ///
    /// - If the next triangle does not contain the location, repeat the step
    ///   above for a maximum number of times equal to two times the number of
    ///   triangles.
    ///
    /// When a triangle enclosing a specific location is not found, the
    /// interpolated value at that location is invalid. The handling of
    /// spherical accurate projection occurs at low-level geometrical functions
    /// (\ref IsPointInPolygonNodes, \ref AreSegmentsCrossing). Therefore, the algorithm is
    /// independent of the implementation details that occur at the level of the
    /// geometrical functions.
    class TriangulationInterpolation
    {

    public:
        /// @brief Constructor
        /// @param[in] locations interpolation points (where the values should be computed)
        /// @param[in] samples  Values to use for the interpolation
        /// @param[in] projection Projection to use (\ref Projection)
        TriangulationInterpolation(const std::vector<Point>& locations,
                                   const std::vector<Sample>& samples,
                                   Projection projection);

        /// @brief Compute results on the interpolation points
        void Compute();

        /// @brief Get the results
        /// @return
        [[nodiscard]] const auto& GetResults() const
        {
            return m_results;
        }

    private:
        const std::vector<Point>& m_locations;
        const std::vector<Sample>& m_samples;
        Projection m_projection;
        std::vector<double> m_results;
    };

}; // namespace meshkernel