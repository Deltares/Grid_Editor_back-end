#include "TestUtils/MakeMeshes.hpp"

#include <gtest/gtest.h>

#include <MeshKernel/Contacts.hpp>
#include <MeshKernel/Entities.hpp>
#include <MeshKernel/Mesh1D.hpp>
#include <MeshKernel/Mesh2D.hpp>
#include <MeshKernel/Polygons.hpp>

TEST(ComputeSingleConnections, ComputeContacts)
{
    // Create 1d mesh
    std::vector<meshkernel::Point> nodes{{-20, 0}, {-20, 10}, {-20, 20}, {-20, 30}, {-20, 40}, {-10, 40}, {0, 40}, {10, 40}};
    std::vector<meshkernel::Edge> edges{{0, 1}, {1, 2}, {2, 3}, {3, 4}, {4, 5}, {5, 6}, {6, 7}};
    const auto mesh1d = std::make_shared<meshkernel::Mesh1D>(edges, nodes, meshkernel::Projection::cartesian);

    // Create 2d mesh
    const auto mesh2d = MakeRectangularMeshForTesting(3, 3, 10, meshkernel::Projection::cartesian, {0.0, 0.0});

    // Create contacts
    std::vector<bool> onedNodeMask(nodes.size(), true);
    meshkernel::Contacts contacts(mesh1d, mesh2d, onedNodeMask);

    // Set the polygon where to generate the contacts
    std::vector<meshkernel::Point> polygonPoints{{-30, -20}, {40, -20}, {40, 50}, {-40, 50}, {-30, -20}};
    meshkernel::Polygons polygon(polygonPoints, meshkernel::Projection::cartesian);

    // Execute
    contacts.ComputeSingleConnections(polygon);

    //Assert
    constexpr double tolerance = 1e-6;
    ASSERT_EQ(0, contacts.m_mesh2dIndices[0]);
}