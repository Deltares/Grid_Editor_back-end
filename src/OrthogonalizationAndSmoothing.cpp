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

#include <vector>
#include <algorithm>
#include <numeric>
#include <string>

#include "Operations.cpp"
#include "Smoother.hpp"
#include "Orthogonalizer.hpp"
#include "OrthogonalizationAndSmoothing.hpp"
#include "Entities.hpp"
#include "Mesh.hpp"
#include "LandBoundaries.hpp"
#include "Polygons.hpp"

meshkernel::OrthogonalizationAndSmoothing::OrthogonalizationAndSmoothing(std::shared_ptr<Mesh> mesh,
                                                                         std::shared_ptr<Smoother> smoother,
                                                                         std::shared_ptr<Orthogonalizer> orthogonalizer,
                                                                         std::shared_ptr<Polygons> polygon,
                                                                         std::shared_ptr<LandBoundaries> landBoundaries,
                                                                         int isTriangulationRequired,
                                                                         int isAccountingForLandBoundariesRequired,
                                                                         int projectToLandBoundaryOption,
                                                                         const meshkernelapi::OrthogonalizationParametersNative& orthogonalizationParametersNative) : m_mesh(mesh),
                                                                                                                                                                      m_smoother(smoother),
                                                                                                                                                                      m_orthogonalizer(orthogonalizer),
                                                                                                                                                                      m_polygons(polygon),
                                                                                                                                                                      m_landBoundaries(landBoundaries),
                                                                                                                                                                      m_projectToLandBoundaryOption(projectToLandBoundaryOption)
{
    m_orthogonalizationToSmoothingFactor = orthogonalizationParametersNative.OrthogonalizationToSmoothingFactor;
    m_orthogonalizationToSmoothingFactorBoundary = orthogonalizationParametersNative.OrthogonalizationToSmoothingFactorBoundary;
    m_smoothorarea = orthogonalizationParametersNative.Smoothorarea;
    m_orthogonalizationOuterIterations = orthogonalizationParametersNative.OuterIterations;
    m_orthogonalizationBoundaryIterations = orthogonalizationParametersNative.BoundaryIterations;
    m_orthogonalizationInnerIterations = orthogonalizationParametersNative.InnerIterations;
}

bool meshkernel::OrthogonalizationAndSmoothing::Initialize()
{
    // Sets the node mask
    m_mesh->Administrate(Mesh::AdministrationOptions::AdministrateMeshEdgesAndFaces);
    m_mesh->MaskNodesInPolygons(*m_polygons, true);

    // Flag nodes outside the polygon as corner points
    for (auto n = 0; n < m_mesh->GetNumNodes(); n++)
    {
        if (m_mesh->m_nodeMask[n] == 0)
        {
            m_mesh->m_nodesTypes[n] = 3;
        }
    }

    // TODO: calculate volume weights for areal smoother
    m_mumax = (1.0 - m_smoothorarea) * 0.5;
    m_mu = std::min(1e-2, m_mumax);
    m_orthogonalCoordinates.resize(m_mesh->GetNumNodes());

    // back-up original nodes, for projection on original mesh boundary
    m_originalNodes = m_mesh->m_nodes;
    m_orthogonalCoordinates = m_mesh->m_nodes;

    // project on land boundary
    if (m_projectToLandBoundaryOption >= 1)
    {
        // account for enclosing polygon
        m_landBoundaries->FindNearestMeshBoundary(m_projectToLandBoundaryOption);
    }

    // for spherical accurate computations we need to call PrapareOuterIteration (orthonet_comp_ops)
    if (m_mesh->m_projection == Projections::sphericalAccurate)
    {
        if (m_orthogonalizationToSmoothingFactor < 1.0)
        {
            bool successful = PrapareOuterIteration();
            if (!successful)
            {
                return false;
            }
        }

        m_localCoordinatesIndexes.resize(m_mesh->GetNumNodes() + 1);
        m_localCoordinatesIndexes[0] = 1;
        for (int n = 0; n < m_mesh->GetNumNodes(); ++n)
        {
            m_localCoordinatesIndexes[n + 1] = m_localCoordinatesIndexes[n] + std::max(m_mesh->m_nodesNumEdges[n] + 1, m_smoother->GetNumConnectedNodes(n));
        }

        m_localCoordinates.resize(m_localCoordinatesIndexes.back() - 1, {doubleMissingValue, doubleMissingValue});
    }

    return true;
}

bool meshkernel::OrthogonalizationAndSmoothing::Compute()
{
    bool successful = true;

    for (auto outerIter = 0; outerIter < m_orthogonalizationOuterIterations; outerIter++)
    {
        if (successful)
        {
            successful = PrapareOuterIteration();
        }
        for (auto boundaryIter = 0; boundaryIter < m_orthogonalizationBoundaryIterations; boundaryIter++)
        {
            for (auto innerIter = 0; innerIter < m_orthogonalizationInnerIterations; innerIter++)
            {
                if (successful)
                {
                    successful = InnerIteration();
                }
            } // inner iteration
        }     // boundary iter

        //update mu
        if (successful)
        {
            successful = FinalizeOuterIteration();
        }
    } // outer iter

    DeallocateLinearSystem();

    return true;
}

bool meshkernel::OrthogonalizationAndSmoothing::PrapareOuterIteration()
{

    bool successful = true;

    // compute weights and rhs of orthogonalizer
    if (successful)
    {
        successful = m_orthogonalizer->Compute();
    }

    // computes the smoother weights
    if (successful)
    {
        successful = m_smoother->Compute();
    }

    // allocate linear system for smoother and orthogonalizer
    if (successful)
    {
        successful = AllocateLinearSystem();
    }

    // compute linear system terms for smoother and orthogonalizer
    if (successful)
    {
        successful = ComputeLinearSystemTerms();
    }
    return successful;
}

bool meshkernel::OrthogonalizationAndSmoothing::AllocateLinearSystem()
{
    bool successful = true;
    // reallocate caches
    if (successful && m_nodeCacheSize == 0)
    {
        m_compressedRhs.resize(m_mesh->GetNumNodes() * 2);
        std::fill(m_compressedRhs.begin(), m_compressedRhs.end(), 0.0);

        m_compressedEndNodeIndex.resize(m_mesh->GetNumNodes());
        std::fill(m_compressedEndNodeIndex.begin(), m_compressedEndNodeIndex.end(), 0.0);

        m_compressedStartNodeIndex.resize(m_mesh->GetNumNodes());
        std::fill(m_compressedStartNodeIndex.begin(), m_compressedStartNodeIndex.end(), 0.0);

        for (int n = 0; n < m_mesh->GetNumNodes(); n++)
        {
            m_compressedEndNodeIndex[n] = m_nodeCacheSize;
            m_nodeCacheSize += std::max(m_mesh->m_nodesNumEdges[n] + 1, m_smoother->GetNumConnectedNodes(n));
            m_compressedStartNodeIndex[n] = m_nodeCacheSize;
        }

        m_compressedNodesNodes.resize(m_nodeCacheSize);
        m_compressedWeightX.resize(m_nodeCacheSize);
        m_compressedWeightY.resize(m_nodeCacheSize);
    }
    return successful;
}

bool meshkernel::OrthogonalizationAndSmoothing::DeallocateLinearSystem()
{
    m_compressedRhs.resize(0);
    m_compressedEndNodeIndex.resize(0);
    m_compressedStartNodeIndex.resize(0);
    m_compressedNodesNodes.resize(0);
    m_compressedWeightX.resize(0);
    m_compressedWeightY.resize(0);
    m_nodeCacheSize = 0;

    return true;
}

bool meshkernel::OrthogonalizationAndSmoothing::FinalizeOuterIteration()
{
    m_mu = std::min(2.0 * m_mu, m_mumax);

    //compute new faces circumcenters
    if (!m_keepCircumcentersAndMassCenters)
    {
        m_mesh->ComputeFaceCircumcentersMassCentersAndAreas();
    }

    return true;
}

bool meshkernel::OrthogonalizationAndSmoothing::ComputeLinearSystemTerms()
{
    double max_aptf = std::max(m_orthogonalizationToSmoothingFactorBoundary, m_orthogonalizationToSmoothingFactor);
#pragma omp parallel for
    for (int n = 0; n < m_mesh->GetNumNodes(); n++)
    {
        if ((m_mesh->m_nodesTypes[n] != 1 && m_mesh->m_nodesTypes[n] != 2) || m_mesh->m_nodesNumEdges[n] < 2)
        {
            continue;
        }
        if (m_keepCircumcentersAndMassCenters != false && (m_mesh->m_nodesNumEdges[n] != 3 || m_mesh->m_nodesNumEdges[n] != 1))
        {
            continue;
        }

        const double atpfLoc = m_mesh->m_nodesTypes[n] == 2 ? max_aptf : m_orthogonalizationToSmoothingFactor;
        const double atpf1Loc = 1.0 - atpfLoc;
        int maxnn = m_compressedStartNodeIndex[n] - m_compressedEndNodeIndex[n];
        auto cacheIndex = m_compressedEndNodeIndex[n];
        for (auto nn = 1; nn < maxnn; nn++)
        {
            double wwx = 0.0;
            double wwy = 0.0;

            // Smoother
            if (atpf1Loc > 0.0 && m_mesh->m_nodesTypes[n] == 1)
            {
                wwx = atpf1Loc * m_smoother->GetWeight(n, nn);
                wwy = atpf1Loc * m_smoother->GetWeight(n, nn);
            }

            // Orthogonalizer
            if (nn < m_mesh->m_nodesNumEdges[n] + 1)
            {
                wwx += atpfLoc * m_orthogonalizer->GetWeight(n, nn - 1);
                wwy += atpfLoc * m_orthogonalizer->GetWeight(n, nn - 1);
                m_compressedNodesNodes[cacheIndex] = m_mesh->m_nodesNodes[n][nn - 1];
            }
            else
            {
                m_compressedNodesNodes[cacheIndex] = m_smoother->GetCoonectedNodeIndex(n, nn);
            }

            m_compressedWeightX[cacheIndex] = wwx;
            m_compressedWeightY[cacheIndex] = wwy;
            cacheIndex++;
        }
        int firstCacheIndex = n * 2;
        m_compressedRhs[firstCacheIndex] = atpfLoc * m_orthogonalizer->GetRightHandSide(n, 0);
        m_compressedRhs[firstCacheIndex + 1] = atpfLoc * m_orthogonalizer->GetRightHandSide(n, 1);
    }

    return true;
}

bool meshkernel::OrthogonalizationAndSmoothing::InnerIteration()
{
#pragma omp parallel for
    for (int n = 0; n < m_mesh->GetNumNodes(); n++)
    {
        UpdateNodeCoordinates(n);
    }

    // update mesh node coordinates
    m_mesh->m_nodes = m_orthogonalCoordinates;

    // project on the original net boundary
    ProjectOnOriginalMeshBoundary();

    // compute local coordinates
    ComputeCoordinates();

    // project on land boundary
    if (m_projectToLandBoundaryOption >= 1)
    {
        m_landBoundaries->SnapMeshToLandBoundaries();
    }

    return true;
}

bool meshkernel::OrthogonalizationAndSmoothing::ProjectOnOriginalMeshBoundary()
{
    Point normalSecondPoint{doubleMissingValue, doubleMissingValue};
    Point normalThirdPoint{doubleMissingValue, doubleMissingValue};

    // in this case the nearest point is the point itself
    std::vector<int> nearestPoints(m_mesh->GetNumNodes(), 0);
    std::iota(nearestPoints.begin(), nearestPoints.end(), 0);

    for (auto n = 0; n < m_mesh->GetNumNodes(); n++)
    {
        int nearestPointIndex = nearestPoints[n];
        if (m_mesh->m_nodesTypes[n] == 2 && m_mesh->m_nodesNumEdges[n] > 0 && m_mesh->m_nodesNumEdges[nearestPointIndex] > 0)
        {
            Point firstPoint = m_mesh->m_nodes[n];
            if (!firstPoint.IsValid())
            {
                continue;
            }

            int numEdges = m_mesh->m_nodesNumEdges[nearestPointIndex];
            int numNodes = 0;
            int leftNode = -1;
            int rightNode = -1;
            Point secondPoint{doubleMissingValue, doubleMissingValue};
            Point thirdPoint{doubleMissingValue, doubleMissingValue};
            for (auto nn = 0; nn < numEdges; nn++)
            {
                auto edgeIndex = m_mesh->m_nodesEdges[nearestPointIndex][nn];
                if (m_mesh->m_edgesNumFaces[edgeIndex] == 1)
                {
                    numNodes++;
                    if (numNodes == 1)
                    {
                        leftNode = m_mesh->m_nodesNodes[n][nn];
                        if (leftNode == intMissingValue)
                        {
                            return false;
                        }
                        secondPoint = m_originalNodes[leftNode];
                    }
                    else if (numNodes == 2)
                    {
                        rightNode = m_mesh->m_nodesNodes[n][nn];
                        if (rightNode == intMissingValue)
                        {
                            return false;
                        }
                        thirdPoint = m_originalNodes[rightNode];
                    }
                }
            }

            if (!secondPoint.IsValid() || !thirdPoint.IsValid())
            {
                continue;
            }

            //Project the moved boundary point back onto the closest original edge (either between 0 and 2 or 0 and 3)
            double rl2 = 0.0;
            double dis2 = DistanceFromLine(firstPoint, m_originalNodes[nearestPointIndex], secondPoint, normalSecondPoint, rl2, m_mesh->m_projection);

            double rl3 = 0.0;
            double dis3 = DistanceFromLine(firstPoint, m_originalNodes[nearestPointIndex], thirdPoint, normalThirdPoint, rl3, m_mesh->m_projection);

            if (dis2 < dis3)
            {
                m_mesh->m_nodes[n] = normalSecondPoint;
                if (rl2 > 0.5 && m_mesh->m_nodesTypes[n] != 3)
                {
                    nearestPoints[n] = leftNode;
                }
            }
            else
            {
                m_mesh->m_nodes[n] = normalThirdPoint;
                if (rl3 > 0.5 && m_mesh->m_nodesTypes[n] != 3)
                {
                    nearestPoints[n] = rightNode;
                }
            }
        }
    }
    return true;
}

bool meshkernel::OrthogonalizationAndSmoothing::ComputeCoordinates() const
{
    if (m_mesh->m_projection == Projections::sphericalAccurate)
    {
        //TODO : missing implementation
        return true;
    }

    return true;
}

bool meshkernel::OrthogonalizationAndSmoothing::UpdateNodeCoordinates(int nodeIndex)
{

    double dx0 = 0.0;
    double dy0 = 0.0;
    double increments[2]{0.0, 0.0};
    ComputeLocalIncrements(nodeIndex, dx0, dy0, increments);

    if (increments[0] <= 1e-8 || increments[1] <= 1e-8)
    {
        return true;
    }

    int firstCacheIndex = nodeIndex * 2;
    dx0 = (dx0 + m_compressedRhs[firstCacheIndex]) / increments[0];
    dy0 = (dy0 + m_compressedRhs[firstCacheIndex + 1]) / increments[1];
    constexpr double relaxationFactor = 0.75;
    if (m_mesh->m_projection == Projections::cartesian || m_mesh->m_projection == Projections::spherical)
    {
        double x0 = m_mesh->m_nodes[nodeIndex].x + dx0;
        double y0 = m_mesh->m_nodes[nodeIndex].y + dy0;
        static constexpr double relaxationFactorCoordinates = 1.0 - relaxationFactor;

        m_orthogonalCoordinates[nodeIndex].x = relaxationFactor * x0 + relaxationFactorCoordinates * m_mesh->m_nodes[nodeIndex].x;
        m_orthogonalCoordinates[nodeIndex].y = relaxationFactor * y0 + relaxationFactorCoordinates * m_mesh->m_nodes[nodeIndex].y;
    }
    if (m_mesh->m_projection == Projections::sphericalAccurate)
    {
        Point localPoint{relaxationFactor * dx0, relaxationFactor * dy0};

        double exxp[3];
        double eyyp[3];
        double ezzp[3];
        ComputeThreeBaseComponents(m_mesh->m_nodes[nodeIndex], exxp, eyyp, ezzp);

        //get 3D-coordinates in rotated frame
        Cartesian3DPoint cartesianLocalPoint;
        SphericalToCartesian(localPoint, cartesianLocalPoint);

        //project to fixed frame
        Cartesian3DPoint transformedCartesianLocalPoint;
        transformedCartesianLocalPoint.x = exxp[0] * cartesianLocalPoint.x + eyyp[0] * cartesianLocalPoint.y + ezzp[0] * cartesianLocalPoint.z;
        transformedCartesianLocalPoint.y = exxp[1] * cartesianLocalPoint.x + eyyp[1] * cartesianLocalPoint.y + ezzp[1] * cartesianLocalPoint.z;
        transformedCartesianLocalPoint.z = exxp[2] * cartesianLocalPoint.x + eyyp[2] * cartesianLocalPoint.y + ezzp[2] * cartesianLocalPoint.z;

        //tranform to spherical coordinates
        CartesianToSpherical(transformedCartesianLocalPoint, m_mesh->m_nodes[nodeIndex].x, m_orthogonalCoordinates[nodeIndex]);
    }
    return true;
}

bool meshkernel::OrthogonalizationAndSmoothing::ComputeLocalIncrements(int nodeIndex, double& dx0, double& dy0, double* weightsSum)
{
    int numConnectedNodes = m_compressedStartNodeIndex[nodeIndex] - m_compressedEndNodeIndex[nodeIndex];
    auto cacheIndex = m_compressedEndNodeIndex[nodeIndex];
    for (int nn = 1; nn < numConnectedNodes; nn++)
    {
        const auto wwx = m_compressedWeightX[cacheIndex];
        const auto wwy = m_compressedWeightY[cacheIndex];
        const auto currentNode = m_compressedNodesNodes[cacheIndex];
        cacheIndex++;

        if (m_mesh->m_projection == Projections::cartesian)
        {
            const double wwxTransformed = wwx;
            const double wwyTransformed = wwy;
            dx0 = dx0 + wwxTransformed * (m_mesh->m_nodes[currentNode].x - m_mesh->m_nodes[nodeIndex].x);
            dy0 = dy0 + wwyTransformed * (m_mesh->m_nodes[currentNode].y - m_mesh->m_nodes[nodeIndex].y);
            weightsSum[0] += wwxTransformed;
            weightsSum[1] += wwyTransformed;
        }

        if (m_mesh->m_projection == Projections::spherical)
        {
            const double wwxTransformed = wwx * earth_radius * degrad_hp *
                                          std::cos(0.5 * (m_mesh->m_nodes[nodeIndex].y + m_mesh->m_nodes[currentNode].y) * degrad_hp);
            const double wwyTransformed = wwy * earth_radius * degrad_hp;

            dx0 = dx0 + wwxTransformed * (m_mesh->m_nodes[currentNode].x - m_mesh->m_nodes[nodeIndex].x);
            dy0 = dy0 + wwyTransformed * (m_mesh->m_nodes[currentNode].y - m_mesh->m_nodes[nodeIndex].y);
            weightsSum[0] += wwxTransformed;
            weightsSum[1] += wwyTransformed;
        }
        if (m_mesh->m_projection == Projections::sphericalAccurate)
        {
            const double wwxTransformed = wwx * earth_radius * degrad_hp;
            const double wwyTransformed = wwy * earth_radius * degrad_hp;

            dx0 = dx0 + wwxTransformed * m_localCoordinates[m_localCoordinatesIndexes[nodeIndex] + currentNode - 1].x;
            dy0 = dy0 + wwyTransformed * m_localCoordinates[m_localCoordinatesIndexes[nodeIndex] + currentNode - 1].y;
            weightsSum[0] += wwxTransformed;
            weightsSum[1] += wwyTransformed;
        }
    }

    return true;
}
