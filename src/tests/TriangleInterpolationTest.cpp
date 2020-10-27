#include "../Mesh.hpp"
#include "../Entities.hpp"
#include "../TriangulationInterpolation.cpp"
#include "SampleFileReader.cpp"
#include "MakeMeshes.cpp"
#include <gtest/gtest.h>

TEST(TriangleInterpolation, InterpolateOnNodes)
{
    // Set up
    std::vector<meshkernel::Sample> samples = ReadSampleFile("..\\..\\tests\\TriangleInterpolationTests\\inTestTriangleInterpolation.xyz");
    auto mesh = ReadLegacyMeshFromFile("..\\..\\tests\\TriangleInterpolationTests\\simple_grid_net.nc");
    ASSERT_GT(mesh->GetNumNodes(), 0);

    meshkernel::TriangulationInterpolation triangulationInterpolation(mesh->m_nodes, samples, meshkernel::Projections::cartesian);
    triangulationInterpolation.Compute();

    const auto results = triangulationInterpolation.GetResults();

    // test internal results
    constexpr double tolerance = 1e-9;
    ASSERT_NEAR(2.0, results[20], tolerance);
    ASSERT_NEAR(2.0, results[21], tolerance);
    ASSERT_NEAR(2.0, results[22], tolerance);
    ASSERT_NEAR(2.0, results[23], tolerance);
    ASSERT_NEAR(2.0, results[24], tolerance);
    ASSERT_NEAR(1.0, results[25], tolerance);
    ASSERT_NEAR(1.0, results[26], tolerance);
    ASSERT_NEAR(1.0, results[27], tolerance);
    ASSERT_NEAR(1.0, results[28], tolerance);
    ASSERT_NEAR(1.0, results[29], tolerance);
}

TEST(TriangleInterpolation, InterpolateOnEdges)
{
    // Set up
    std::vector<meshkernel::Sample> samples = ReadSampleFile("..\\..\\tests\\TriangleInterpolationTests\\inTestTriangleInterpolation.xyz");
    auto mesh = ReadLegacyMeshFromFile("..\\..\\tests\\TriangleInterpolationTests\\simple_grid_net.nc");
    ASSERT_GT(mesh->GetNumNodes(), 0);

    mesh->ComputeEdgesCenters();

    meshkernel::TriangulationInterpolation triangulationInterpolation(mesh->m_edgesCenters, samples, meshkernel::Projections::cartesian);
    triangulationInterpolation.Compute();

    const auto results = triangulationInterpolation.GetResults();

    // test internal results
    constexpr double tolerance = 1e-9;
    ASSERT_NEAR(2.0, results[14], tolerance);
    ASSERT_NEAR(2.0, results[15], tolerance);
    ASSERT_NEAR(2.0, results[16], tolerance);
    ASSERT_NEAR(2.0, results[17], tolerance);
    ASSERT_NEAR(2.0, results[18], tolerance);
    ASSERT_NEAR(2.0, results[19], tolerance);
    ASSERT_NEAR(2.0, results[20], tolerance);
    ASSERT_NEAR(2.0, results[21], tolerance);
    ASSERT_NEAR(2.0, results[22], tolerance);
    ASSERT_NEAR(2.0, results[23], tolerance);
}
TEST(TriangleInterpolation, InterpolateOnFaces)
{
    // Set up
    std::vector<meshkernel::Sample> samples = ReadSampleFile("..\\..\\tests\\TriangleInterpolationTests\\inTestTriangleInterpolation.xyz");
    auto mesh = ReadLegacyMeshFromFile("..\\..\\tests\\TriangleInterpolationTests\\simple_grid_net.nc");
    ASSERT_GT(mesh->GetNumNodes(), 0);

    meshkernel::TriangulationInterpolation triangulationInterpolation(mesh->m_facesMassCenters, samples, meshkernel::Projections::cartesian);
    triangulationInterpolation.Compute();

    const auto results = triangulationInterpolation.GetResults();

    // test internal results
    constexpr double tolerance = 1e-9;
    ASSERT_NEAR(2.0, results[15], tolerance);
    ASSERT_NEAR(2.0, results[16], tolerance);
    ASSERT_NEAR(2.0, results[17], tolerance);
    ASSERT_NEAR(2.0, results[18], tolerance);
    ASSERT_NEAR(2.0, results[19], tolerance);
    ASSERT_NEAR(0.0, results[20], tolerance);
    ASSERT_NEAR(0.0, results[21], tolerance);
    ASSERT_NEAR(0.0, results[22], tolerance);
    ASSERT_NEAR(0.0, results[23], tolerance);
    ASSERT_NEAR(0.0, results[24], tolerance);
    ASSERT_NEAR(0.0, results[25], tolerance);
}