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

#include "Operations.cpp"
#include "SpatialTrees.hpp"

namespace meshkernel
{
    // Forward declarations
    class Polygons;
    class Mesh;

    class Averaging
    {
    public:
        /// @brief Interpolation based on averaging
        /// @param mesh The input mesh.
        /// @param samples The samples with x,y locations and values.
        /// @param averagingMethod The averaging method to use.
        /// @param locationType The location type (faces, edges, nodes).
        /// @param relativeSearchRadius The relative search radius, used to enlarge the search area when looking for samples.
        /// @param useClosestSampleIfNoneAvailable If no sample are found simply use the closest one.
        /// @param subtractSampleValues For some algorithms (e.g. refinement based on levels) we need to subtract 1 the sample value.
        explicit Averaging(std::shared_ptr<Mesh> mesh,
                           std::vector<Sample>& samples,
                           AveragingMethod averagingMethod,
                           InterpolationLocation locationType,
                           double relativeSearchRadius,
                           bool useClosestSampleIfNoneAvailable,
                           bool subtractSampleValues);

        /// @brief Compute interpolation
        /// @return If the method succeded
        bool Compute();

        /// @brief Get the result values
        /// @return
        const auto& GetResults() const
        {
            return m_results;
        }

    private:
        /// @brief Compute The averaging results in polygon.
        /// @param polygon The bounding polygon where the samples are included.
        /// @param interpolationPoint The interpolation point.
        /// @param result The resulting value
        /// @return
        bool ComputeOnPolygon(const std::vector<Point>& polygon,
                              Point interpolationPoint,
                              double& result);

        std::shared_ptr<Mesh> m_mesh;
        std::vector<Sample>& m_samples;
        AveragingMethod m_averagingMethod;
        InterpolationLocation m_interpolationLocation;
        double m_relativeSearchRadius;
        bool m_useClosestSampleIfNoneAvailable = false;
        bool m_transformSamples = false;

        SpatialTrees::RTree m_samplesRtree;
        std::vector<double> m_results;
        std::vector<bool> m_visitedSamples;
    };
} // namespace meshkernel