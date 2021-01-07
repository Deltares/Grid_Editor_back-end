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

namespace meshkernel
{
    class Mesh;

    /// <summary>
    /// Orthogonalizion (optimize the aspect ratios) and mesh smoothing (optimize internal face angles or area).
    /// </summary>
    class Smoother
    {

    public:
        /// @brief Mesh ctor
        /// @brief mesh
        /// @returns
        explicit Smoother(std::shared_ptr<Mesh> mesh);

        /// @brief Computes the smoother weights
        void Compute();

        /// @brief Gets the weight for a certain node and connected node
        /// @brief node
        /// @brief connectedNode
        /// @returns
        [[nodiscard]] auto GetWeight(size_t node, int connectedNode)
        {
            return m_weights[node][connectedNode];
        }

        /// @brief Get the index of the connected node as assigned by the smoother administration
        /// @brief node
        /// @brief connectedNode
        /// @returns
        [[nodiscard]] auto GetConnectedNodeIndex(size_t node, int connectedNode)
        {
            return m_connectedNodes[node][connectedNode];
        }

        /// @brief Get number of connected nodes
        /// @brief node
        /// @returns
        [[nodiscard]] auto GetNumConnectedNodes(size_t node)
        {
            return m_numConnectedNodes[node];
        }

    private:
        /// @brief Initialize smoother topologies. A topology is determined by how many nodes are connected to the current node.
        ///        There are at maximum mesh.m_numNodes topologies, most likely much less
        void Initialize();

        /// @brief Computes all topologies of the elliptic smoother
        void ComputeTopologies();

        /// @brief Computes all operators of the elliptic smoother
        void ComputeOperators();

        /// @brief Compute nodes local coordinates, sice-effects only for sphericalAccurate projection (comp_local_coords)
        /// @brief mesh
        /// @returns
        bool ComputeCoordinates();

        /// @brief Computes the smoother weights from the operators (orthonet_compweights_smooth)
        void ComputeWeights();

        /// Computes operators of the elliptic smoother by node (orthonet_comp_operators)
        /// @param[in] currentNode
        void ComputeOperatorsNode(size_t currentNode);

        /// @brief Computes m_faceNodeMappingCache, m_sharedFacesCache, m_connectedNodes for the current node, required before computing xi and eta
        /// @param[in] currentNode
        /// @param[out] numSharedFaces
        /// @param[out] numConnectedNodes
        void NodeAdministration(size_t currentNode,
                                size_t& numSharedFaces,
                                size_t& numConnectedNodes);

        /// @brief Compute compute current node xi and eta (orthonet_assign_xieta)
        /// @param[in] currentNode
        /// @param[in] numSharedFaces
        /// @param[in] numConnectedNodes
        void ComputeNodeXiEta(size_t currentNode,
                              size_t numSharedFaces,
                              size_t numConnectedNodes);

        /// @brief Compute optimal edge angle
        /// @brief numFaceNodes
        /// @brief theta1
        /// @brief theta2
        /// @brief isBoundaryEdge
        /// @returns
        [[nodiscard]] double OptimalEdgeAngle(size_t numFaceNodes,
                                              double theta1 = -1.0,
                                              double theta2 = -1.0,
                                              bool isBoundaryEdge = false) const;

        /// @brief Allocate smoother operators
        /// @param[in] topologyIndex
        void AllocateNodeOperators(size_t topologyIndex);

        /// @brief If it is a new topology, save it
        /// @param[in] currentNode
        /// @param[in] numSharedFaces
        /// @param[in] numConnectedNodes
        void SaveNodeTopologyIfNeeded(size_t currentNode,
                                      size_t numSharedFaces,
                                      size_t numConnectedNodes);

        /// @brief Computes local coordinates jacobian from the mapped jacobians m_Jxi and m_Jeta
        /// @param[in] currentNode
        /// @param[out] J
        void ComputeJacobian(size_t currentNode, std::vector<double>& J) const;

        /// @brief Compute the matrix norm
        /// @brief x
        /// @brief y
        /// @brief matCoefficents
        /// @returns
        [[nodiscard]] double MatrixNorm(const std::vector<double>& x,
                                        const std::vector<double>& y,
                                        const std::vector<double>& matCoefficents) const;

        // The mesh to smooth
        std::shared_ptr<Mesh> m_mesh;

        // Smoother weights
        std::vector<std::vector<double>> m_weights;

        // Smoother operators
        std::vector<std::vector<std::vector<double>>> m_Gxi;  // Node to edge xi derivative
        std::vector<std::vector<std::vector<double>>> m_Geta; // Node to edge etha derivative
        std::vector<std::vector<double>> m_Divxi;             // Edge to node xi derivative
        std::vector<std::vector<double>> m_Diveta;            // Edge to node etha derivative
        std::vector<std::vector<std::vector<double>>> m_Az;   // Coefficients to estimate values at cell circumcenters
        std::vector<std::vector<double>> m_Jxi;               // Node to node xi derivative (Jacobian)
        std::vector<std::vector<double>> m_Jeta;              // Node to node eta derivative (Jacobian)
        std::vector<std::vector<double>> m_ww2;               // weights

        // Smoother local caches
        std::vector<size_t> m_sharedFacesCache;
        std::vector<size_t> m_connectedNodesCache;
        std::vector<std::vector<size_t>> m_faceNodeMappingCache;
        std::vector<double> m_xiCache;
        std::vector<double> m_etaCache;
        std::vector<size_t> m_boundaryEdgesCache;
        std::vector<double> m_leftXFaceCenterCache;
        std::vector<double> m_leftYFaceCenterCache;
        std::vector<double> m_rightXFaceCenterCache;
        std::vector<double> m_rightYFaceCenterCache;
        std::vector<double> m_xisCache;
        std::vector<double> m_etasCache;

        // Smoother topologies
        size_t m_numTopologies = 0;
        std::vector<size_t> m_nodeTopologyMapping;
        std::vector<size_t> m_numTopologyNodes;
        std::vector<size_t> m_numTopologyFaces;
        std::vector<std::vector<double>> m_topologyXi;
        std::vector<std::vector<double>> m_topologyEta;
        std::vector<std::vector<size_t>> m_topologySharedFaces;
        std::vector<std::vector<std::vector<size_t>>> m_topologyFaceNodeMapping;
        std::vector<std::vector<size_t>> m_topologyConnectedNodes;

        std::vector<size_t> m_numConnectedNodes;           // (nmk2)
        std::vector<std::vector<size_t>> m_connectedNodes; // (kk2)

        // Class variables
        size_t m_maximumNumConnectedNodes = 0;
        size_t m_maximumNumSharedFaces = 0;

        static constexpr int m_topologyInitialSize = 10;
        static constexpr double m_thetaTolerance = 1e-4;
    };
} // namespace meshkernel
