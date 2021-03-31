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

#include <MeshKernel/CurvilinearGrid.hpp>
#include <MeshKernel/Entities.hpp>

namespace meshkernel
{

    /// @brief A class implementing the curvilinear grid smoothing algorithm
    class CurvilinearGridSmoothing
    {
    public:
        /// @brief Class constructor
        /// @param[in] grid                        The input curvilinear grid
        /// @param[in] smoothingIterations         The number of smoothing iterations to perform
        CurvilinearGridSmoothing(std::shared_ptr<CurvilinearGrid> grid,
                                 size_t smoothingIterations);

        /// @brief Compute curvilinear grid block smoothing (modifies the m_grid nodal values)
        void Compute();

        /// @brief Compute curvilinear grid line smoothing. The algorithm smooths the grid along the direction specified by the line.
        /// The line must be an m or n grid line of the curvilinear grid.
        /// The grid is smoothed in the area specified by the line in one direction and by the two corner points on the other direction.
        /// @param[in] firstSegmentNode The first point of the line
        /// @param[in] secondSegmentNode The second point of the line
        /// @param[in] lowerLeftCornerSmoothingArea The left point of the smoothing area
        /// @param[in] upperRightCornerSmootingArea The right point of the smoothing area
        void ComputedDirectionalSmooth(Point const& firstSegmentNode,
                                       Point const& secondSegmentNode,
                                       Point const& lowerLeftCornerSmoothingArea,
                                       Point const& upperRightCornerSmootingArea);

        /// @brief Sets the smoothing block
        /// @param[in] firstCornerPoint The first point defining the smoothing area
        /// @param[in] secondCornerPoint The second point defining the smoothing area
        void SetBlock(Point const& firstCornerPoint, Point const& secondCornerPoint);

    private:
        /// @brief Solve one iteration of block smoothing
        void Solve();

        /// @brief Solve one iteration of directional smoothing
        /// @param[in] isSmoothingAlongM True if smoothing is along M, false otherwise
        /// @param[in] pointOnLineIndices A point on the line defining the directional smooth
        /// @param[in] lowerLeftCornerRegion The lower left corner of the smoothing region
        /// @param[in] upperRightCornerSmoothingRegion The upper right corner of the smoothing region
        void SolveDirectionalSmooth(bool isSmoothingAlongM,
                                    CurvilinearGrid::NodeIndices const& pointOnLineIndices,
                                    CurvilinearGrid::NodeIndices const& lowerLeftCornerRegion,
                                    CurvilinearGrid::NodeIndices const& upperRightCornerSmoothingRegion);

        /// @brief Function for computing the smoothing factors at the current location given a line and a zone of influence (SMEERFUNCTIE)
        /// The smoothing factor is maximum at the line and 0 at the boundary of the smoothing zone.
        /// @param[in] currentPointIndices The indices of the current point
        /// @param[in] pointOnSmoothingLineIndices The indices of a point on the smoothing line
        /// @param[in] lowerLeftIndices The lower left indices of the smoothing area
        /// @param[in] upperRightIndices The upper right indices of the smoothing area
        /// @return A tuple containing the horizontal, the vertical and mixed smoothing factors
        std::tuple<double, double, double> ComputeDirectionalSmoothingFactors(CurvilinearGrid::NodeIndices const& currentPointIndices,
                                                                              CurvilinearGrid::NodeIndices const& pointOnSmoothingLineIndices,
                                                                              CurvilinearGrid::NodeIndices const& lowerLeftIndices,
                                                                              CurvilinearGrid::NodeIndices const& upperRightIndices) const;

        /// @brief Projects a point on the closest grid boundary
        /// @param[in] point The point to project
        /// @param[in] m The current m coordinate on the boundary of the curvilinear grid
        /// @param[in] n The current n coordinate on the boundary of the curvilinear grid
        void ProjectPointOnClosestGridBoundary(Point const& point, size_t m, size_t n);

        std::shared_ptr<CurvilinearGrid> m_grid; ///< A pointer to the curvilinear grid to modify
        size_t m_smoothingIterations;            ///< The orthogonalization parameters

        CurvilinearGrid::NodeIndices m_lowerLeft;  ///< The lower left corner of the smoothing block, used in grid block smoothing
        CurvilinearGrid::NodeIndices m_upperRight; ///< The upper right corner of the smoothing block, used in grid block smoothing

        std::vector<std::vector<Point>> m_gridNodesCache; ///< A cache for storing current iteration node positions
    };
} // namespace meshkernel