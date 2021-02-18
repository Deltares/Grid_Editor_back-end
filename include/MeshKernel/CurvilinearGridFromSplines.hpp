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
#include <vector>

#include <MeshKernel/Entities.hpp>
#include <MeshKernelApi/CurvilinearParameters.hpp>
#include <MeshKernelApi/SplinesToCurvilinearParameters.hpp>

namespace meshkernel
{
    class CurvilinearGrid;
    class Splines;

    /// @brief A class used to create a curvilinear grid from a central spline
    /// Usually this class should be preferred over CurvilinearGridFromSplinesTransfinite.
    ///
    /// In this class, the algorithm to gradually develop a mesh from the
    /// centreline of the channel towards the boundaries is implemented. It is
    /// the most complex algorithm in MeshKernel. The curvilinear mesh is
    /// developed from the center spline by the following steps:
    ///
    /// - Initialization step
    ///
    ///     - The splines are labeled (central or transversal spline) based
    ///       on the number of corner points and the intersecting angles.
    ///
    ///     - The grid heights at a different position along the central
    ///       spline are computed from the crossing splines.
    ///
    ///     - The normal vectors of each edge are computed, as these
    ///       determine the growing front directions.
    ///
    ///     - The edge velocities to apply to each normal direction are
    ///       computed.
    ///
    /// - Iteration step, where the mesh is grown of one layer at the time
    ///   from the left and right sides of the central spline:
    ///
    ///     - Compute the node velocities from the edge velocities.
    ///
    ///     - Find the nodes at the front (the front might miss some faces and
    ///       be irregular).
    ///
    ///     - Compute the maximum growth time to avoid faces with intersecting
    ///       edges.
    ///
    ///     - Grow the grid by translating the nodes at the front by an amount
    ///       equal to the product of the nodal velocity by the maximum grow
    ///       time.
    ///
    /// - Post-processing
    ///
    ///     - Remove the skewed faces whose aspect ratio exceeds a prescribed
    ///       value.
    ///
    ///     - Compute the resulting CurvilinearGrid from the internal table of
    ///       computed nodes (m_gridPoints).
    class CurvilinearGridFromSplines
    {
    public:
        /// @brief Ctor
        /// @param splines Input splines
        /// @param curvilinearParameters The parameters for OrthogonalCurvilinearGridFromSplines algorithm
        /// @param splinesToCurvilinearParameters The parameters for OrthogonalCurvilinearGridFromSplines algorithm
        CurvilinearGridFromSplines(std::shared_ptr<Splines> splines,
                                   const meshkernelapi::CurvilinearParameters& curvilinearParameters,
                                   const meshkernelapi::SplinesToCurvilinearParameters& splinesToCurvilinearParameters);

        /// Computes the spline properties, such as cross splines (get_splineprops)
        /// @param restoreOriginalProperties
        void ComputeSplineProperties(bool restoreOriginalProperties);

        /// @brief Computes a curvilinear grid using the growing front method (spline2curvi).
        ///
        /// OrthogonalCurvilinearGrid algorithm.
        /// 1. Compute spline properties (the crossings)
        /// 2. Make all grid lines of the central spline
        /// 3. Adds artificial splines
        /// 4. Compute properties with artificial splines added
        /// 5. Compute the edge velocities
        /// 6. Grow layers
        /// @returns The computed curvilinear grid
        CurvilinearGrid Compute();

        /// @brief Computes the curvilinear grid from the mesh points
        /// @returns The computed curvilinear grid
        CurvilinearGrid ComputeCurvilinearGridFromGridPoints();

        /// @brief Initialize the OrthogonalCurvilinearGrid algorithm.
        void Initialize();

        /// @brief Performs one iteration for generating another layer on the advancing fronts
        /// @param[in] layer The index of the layer to be generated
        void Iterate(size_t layer);

        /// @brief For the central spline, computes the spline subdivisions along the spline (make_wholegridline)
        void MakeAllGridLines();

        std::vector<Point> m_gridLine;                        ///< Coordinates of the first gridline (xg1, yg1)
        std::vector<double> m_gridLineDimensionalCoordinates; ///< Center spline coordinates of the first gridline (sg1)
        std::vector<double> m_maximumGridHeights;             ///< Maximum transversal grid height ()
        size_t m_numM = 0;                                    ///< Number of columns in the curvilinear grid

        std::shared_ptr<Splines> m_splines; ///< A pointer to spline class instance

    private:
        /// @brief From the layerIndex index gets the next grid layerIndex and the transversal sublayer index (get_isub)
        /// @param[in] layerIndex The current grid layerIndex index
        /// @returns The next grid layerIndex and the sub layerIndex index
        std::tuple<size_t, size_t> ComputeGridLayerAndSubLayer(size_t layerIndex);

        /// @brief Grow a layer starting from a given layer index
        void GrowLayer(size_t layerIndex);

        /// @brief Compute the maximum allowable time step when growing a grid (comp_tmax_self)
        /// @param[in] coordinates The starting point coordinates
        /// @param[in] velocities The velocities at the coordinates
        /// @returns The maximum allowable time step for each edge
        std::vector<double> ComputeMaximumEdgeGrowTime(const std::vector<Point>& coordinates,
                                                       const std::vector<Point>& velocities) const;

        /// @brief Copy growth velocities to the advancing front, add points at front corners corners (copy_vel_to_front)
        /// @brief[in] layerIndex The current grid layerIndex index
        /// @brief[in] previousFrontVelocities The previous front velocities
        /// @returns The front velocities for the next front
        std::vector<Point> CopyVelocitiesToFront(size_t layerIndex, const std::vector<Point>& previousFrontVelocities);

        /// @brief Computes the points at front, which have to be moved.
        /// @returns The indices of the grid points, the front grid points and the number of front points
        std::tuple<std::vector<std::vector<size_t>>, std::vector<Point>, size_t> FindFront();

        /// @brief Compute growth velocity vectors at grid points (comp_vel)
        /// @param[in] layerIndex The current grid layerIndex index
        /// @returns The velocities at the grid point
        std::vector<Point> ComputeVelocitiesAtGridPoints(size_t layerIndex);

        /// @brief Get left and right points at given layer for a given index (get_LR)
        /// @brief[in] gridPoints The front grid points
        /// @brief[in] index The front grid point index from which the neighbors should be searched
        /// @returns The left and the right index for the current front grid point index
        std::tuple<size_t, size_t> GetNeighbours(const std::vector<Point>& gridPoints, size_t index) const;

        /// @brief Compute the edge grow velocities (comp_edgevel)
        /// TODO: can this be split in compute heights and computeGrowFactors
        void ComputeEdgeVelocities();

        /// @brief Compute the grid grow factor for a given total grid height, first grid layer height and number of grid layers (comp_dgrow)
        /// @param[in] The current spline index
        /// @returns The computed grow factor
        double ComputeGrowFactor(size_t splineIndex) const;

        /// @brief Computes the exponential grid height from the aspect ratio, the heights and the number of grid layers
        /// @brief[in] aspectRatio The aspect ratio
        /// @brief[in] height The current height
        /// @brief[in] numLayers The number of layers
        /// @returns The new grid height
        [[nodiscard]] double ComputeTotalExponentialHeight(double aspectRatio,
                                                           double height,
                                                           size_t numLayers) const;

        /// @brief Compute the number of exponential layers for a given heightRatio
        /// @brief[in] heightRatio The max height ratio between layers
        /// @returns the number of exponential layers
        [[nodiscard]] size_t ComputeNumberExponentialLayers(double heightRatio) const;

        /// @brief Computes the sub-interval velocities (left and right)
        /// @param[in] s
        /// @param[in] startGridLineIndex
        /// @param[in] endGridLineIndex
        /// @param[in] numHeights
        /// @param[in] numOtherSideHeights
        /// @param[in] firstHeight
        /// @param[in] gridLineIndex
        /// @param[in] otherGridLineIndex
        /// @param[out] numPerpendicularFacesOnSubintervalAndEdge
        /// @param[out] edgeVelocities
        /// @param[out] hh0MaxRatio
        void ComputeVelocitiesSubIntervals(size_t s,
                                           size_t startGridLineIndex,
                                           size_t endGridLineIndex,
                                           size_t numHeights,
                                           size_t numOtherSideHeights,
                                           double firstHeight,
                                           const std::vector<size_t>& gridLineIndex,
                                           const std::vector<size_t>& otherGridLineIndex,
                                           std::vector<std::vector<size_t>>& numPerpendicularFacesOnSubintervalAndEdge,
                                           std::vector<double>& edgeVelocities,
                                           double& hh0MaxRatio);

        /// @brief Compute the grid heights at grid edges on the center spline (comp_gridheights)
        void ComputeGridHeights();

        /// @brief Find the nearest cross spline
        /// @param[in] s
        /// @param[in] j
        /// @param[in] numHeightsLeft
        /// @param[in] edgesCenterPoints
        /// @param[in] crossSplineLeftHeights
        /// @param[in] localValidSplineIndices
        /// @param[out] localSplineDerivatives
        /// @param[out] crossingSplinesDimensionalCoordinates
        /// @param[out] heights
        void FindNearestCrossSplines(size_t s,
                                     size_t j,
                                     const std::vector<size_t>& numHeightsLeft,
                                     const std::vector<std::vector<double>>& crossSplineLeftHeights,
                                     const std::vector<double>& edgesCenterPoints,
                                     std::vector<size_t>& localValidSplineIndices,
                                     std::vector<double>& localSplineDerivatives,
                                     std::vector<double>& crossingSplinesDimensionalCoordinates,
                                     std::vector<std::vector<double>>& heights);

        /// @brief Computes the intersections on a given spline (get_crosssplines)
        /// @brief[in] splineIndex The current spline index
        void GetSplineIntersections(size_t splineIndex);

        /// @brief Generate a gridline on a spline with a prescribed maximum mesh width (make_gridline)
        /// @param[in] splineIndex  The current spline index
        /// @param[in] startingIndex The start index in the gridline
        /// @returns The number of m subdivisions
        size_t MakeGridLine(size_t splineIndex,
                            size_t startingIndex);

        /// @brief Compute the grid heights using ComputeSubHeights and calculates the maximum sub height (get_heights)
        void ComputeHeights();

        /// @brief Compute the heights at left and right of the center spline index.
        ///
        /// Heights are determined by the lengths of the left and right parts of the crossing splines (comp_subheights)
        ///
        /// @param[in] centerSplineIndex
        /// @param[in] crossingSplineLocalIndex
        void ComputeSubHeights(size_t centerSplineIndex, size_t crossingSplineLocalIndex);

        /// @brief Delete skewed cells and cells whose aspect ratio exceeds a prescribed value (postgrid)
        void DeleteSkinnyTriangles();

        /// @brief Allocate spline properties arrays
        void AllocateSplinesProperties();

        /// @brief The spline type
        enum class SplineTypes
        {
            central,
            crossing,
            artificial,
            lateral
        };

        // algorithm parameters
        meshkernelapi::CurvilinearParameters m_curvilinearParameters;                   ///< Curvilinear parameters
        meshkernelapi::SplinesToCurvilinearParameters m_splinesToCurvilinearParameters; ///< Splines to curvilinear parameters

        const size_t m_maxNumCenterSplineHeights = 10; ///< Nsubmax, naz number of different heights a cross spline can have (is determined by how many crossing spline the user can input)
        const size_t m_maxNUniformPart = 5;            ///< Maximum number of layers in the uniform part
        double m_onTopOfEachOtherSquaredTolerance;     ///< On top of each other squared tolerance
        size_t m_numOriginalSplines = 0;               ///< The original number of splines

        // Spline properties (first index is the spline number)
        std::vector<SplineTypes> m_type;                           ///< Spline type
        std::vector<int> m_centralSplineIndex;                     ///< For each spline the index to its central
        std::vector<size_t> m_numCrossingSplines;                  ///< Number of crossing splines
        std::vector<std::vector<size_t>> m_crossingSplinesIndices; ///< Indices for each cross spline, the indices of the center splines

        std::vector<std::vector<bool>> m_isLeftOriented;                         ///< isLeftOriented cross spline is left to right(.true.) or not (.false.) w.r.t.center spline
        std::vector<std::vector<double>> m_crossSplineCoordinates;               ///< t center spline coordinates of cross splines
        std::vector<std::vector<double>> m_cosCrossingAngle;                     ///< cosPhi cosine of crossing angle
        std::vector<std::vector<std::vector<double>>> m_crossSplineLeftHeights;  ///< hL left - hand side grid heights at cross spline locations for each grid layer subinterval, hL(1, :) being the height of the first subinterval, etc.
        std::vector<std::vector<std::vector<double>>> m_crossSplineRightHeights; ///< hR right - hand side grid heights at cross spline locations for each grid layer subinterval, hR(1, :) being the height of the first subinterval, etc.
        std::vector<std::vector<size_t>> m_numCrossSplineLeftHeights;            ///< NsubL number of subintervals of grid layers at cross spline locations at the left - hand side of the spline, each having their own exponential grow factor
        std::vector<std::vector<size_t>> m_numCrossSplineRightHeights;           ///< NsubR number of subintervals of grid layers at cross spline locations at the right - hand side of the spline, each having their own exponential grow factor
        std::vector<size_t> m_numMSplines;                                       ///< mfac number of grid intervals on the spline
        std::vector<size_t> m_leftGridLineIndex;                                 ///< iL index in the whole gridline array of the first grid point on the left - hand side of the spline
        std::vector<size_t> m_rightGridLineIndex;                                ///< iR index in the whole gridline array of the first grid point on the right - hand side of the spline
        std::vector<std::vector<double>> m_gridHeights;                          ///< Heights of all grid elements

        std::vector<size_t> m_leftGridLineIndexOriginal;  ///< Original index of the left grid line
        std::vector<size_t> m_rightGridLineIndexOriginal; ///< Original index of the right grid line
        std::vector<size_t> m_mfacOriginal;               ///< Original mfac number
        std::vector<double> m_maximumGridHeightsOriginal; ///< Original maximum transversal grid height
        std::vector<SplineTypes> m_originalTypes;         ///< Original types

        //cache variables during iterations
        std::vector<double> m_edgeVelocities;                                         ///< Edge velocities
        std::vector<size_t> m_validFrontNodes;                                        ///< Valid front nodes
        std::vector<std::vector<Point>> m_gridPoints;                                 ///< Grid points
        double m_timeStep = 1.0;                                                      ///< Time step
        std::vector<size_t> m_subLayerGridPoints;                                     ///< Sublayer grid points
        std::vector<std::vector<size_t>> m_numPerpendicularFacesOnSubintervalAndEdge; ///< Perpendicular faces on subinterval and edge
        std::vector<std::vector<double>> m_growFactorOnSubintervalAndEdge;            ///< Grow factor on subinterval and edge
    };
} // namespace meshkernel
