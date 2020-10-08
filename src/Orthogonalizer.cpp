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
#include <algorithm>
#include <numeric>
#include <string>

#include "Operations.cpp"
#include "Constants.cpp"
#include "Mesh.hpp"
#include "Entities.hpp"
#include "Orthogonalizer.hpp"

meshkernel::Orthogonalizer::Orthogonalizer(std::shared_ptr<Mesh> mesh) : m_mesh(mesh)
{
}

bool meshkernel::Orthogonalizer::Compute()
{
    m_mesh->ComputeNodeNeighbours();
    m_weights.resize(m_mesh->GetNumNodes(), std::vector<double>(m_mesh->m_maxNumNeighbours, 0.0));
    m_rhs.resize(m_mesh->GetNumNodes(), std::vector<double>(2, 0.0));
    std::fill(m_rhs.begin(), m_rhs.end(), std::vector<double>(2, 0.0));

    // Compute mesh aspect ratios
    m_mesh->GetAspectRatios(m_aspectRatios);

    for (auto n = 0; n < m_mesh->GetNumNodes(); n++)
    {
        if (m_mesh->m_nodesTypes[n] != 1 && m_mesh->m_nodesTypes[n] != 2)
        {
            continue;
        }

        for (auto nn = 0; nn < m_mesh->m_nodesNumEdges[n]; nn++)
        {

            const auto edgeIndex = m_mesh->m_nodesEdges[n][nn];
            double aspectRatio = m_aspectRatios[edgeIndex];
            m_weights[n][nn] = 0.0;

            if (aspectRatio != doubleMissingValue)
            {
                // internal nodes
                m_weights[n][nn] = aspectRatio;

                if (m_mesh->m_edgesNumFaces[edgeIndex] == 1)
                {
                    // boundary nodes
                    m_weights[n][nn] = 0.5 * aspectRatio;

                    // compute the edge length
                    Point neighbouringNode = m_mesh->m_nodes[m_mesh->m_nodesNodes[n][nn]];
                    double neighbouringNodeDistance = Distance(neighbouringNode, m_mesh->m_nodes[n], m_mesh->m_projection);

                    auto leftFace = m_mesh->m_edgesFaces[edgeIndex][0];
                    bool flippedNormal;
                    Point normal;
                    NormalVectorInside(m_mesh->m_nodes[n], neighbouringNode, m_mesh->m_facesMassCenters[leftFace], normal, flippedNormal, m_mesh->m_projection);

                    if (m_mesh->m_projection == Projections::spherical && m_mesh->m_projection != Projections::sphericalAccurate)
                    {
                        normal.x = normal.x * std::cos(degrad_hp * 0.5 * (m_mesh->m_nodes[n].y + neighbouringNode.y));
                    }

                    m_rhs[n][0] += neighbouringNodeDistance * normal.x * 0.5;
                    m_rhs[n][1] += neighbouringNodeDistance * normal.y * 0.5;
                }
            }
        }

        // normalize
        double factor = std::accumulate(m_weights[n].begin(), m_weights[n].end(), 0.0);
        if (std::abs(factor) > 1e-14)
        {
            factor = 1.0 / factor;
            for (auto& w : m_weights[n])
                w = w * factor;
            m_rhs[n][0] = factor * m_rhs[n][0];
            m_rhs[n][1] = factor * m_rhs[n][1];
        }
    }
    return true;
}