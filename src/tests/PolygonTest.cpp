#include "../Mesh.hpp"
#include "../Entities.hpp"
#include "../Polygons.cpp"
#include "../Orthogonalization.hpp"
#include "../Splines.hpp"
#include <gtest/gtest.h>

TEST(Polygons, MeshBoundaryToPolygon)
{
    //One gets the edges
    std::vector<GridGeom::Point> nodes;

    nodes.push_back(GridGeom::Point{ 322.252624511719,454.880187988281 });
    nodes.push_back(GridGeom::Point{ 227.002044677734,360.379241943359 });
    nodes.push_back(GridGeom::Point{ 259.252227783203,241.878051757813 });
    nodes.push_back(GridGeom::Point{ 428.003295898438,210.377746582031 });
    nodes.push_back(GridGeom::Point{ 536.003967285156,310.878753662109 });
    nodes.push_back(GridGeom::Point{ 503.753784179688,432.379974365234 });
    nodes.push_back(GridGeom::Point{ 350.752807617188,458.630249023438 });
    nodes.push_back(GridGeom::Point{ 343.15053976393,406.232256102912 });
    nodes.push_back(GridGeom::Point{ 310.300984548069,319.41005739802 });
    nodes.push_back(GridGeom::Point{ 423.569603308318,326.17986967523 });

    std::vector<GridGeom::Edge> edges;
    // Local edges
    edges.push_back({ 3, 9 });
    edges.push_back({ 9, 2 });
    edges.push_back({ 2, 3 });
    edges.push_back({ 3, 4 });
    edges.push_back({ 4, 9 });
    edges.push_back({ 2, 8 });
    edges.push_back({ 8, 1 });
    edges.push_back({ 1, 2 });
    edges.push_back({ 9, 8 });
    edges.push_back({ 8, 7 });
    edges.push_back({ 7, 1 });
    edges.push_back({ 9, 10 });
    edges.push_back({ 10, 8 });
    edges.push_back({ 4, 5 });
    edges.push_back({ 5, 10 });
    edges.push_back({ 10, 4 });
    edges.push_back({ 8, 6 });
    edges.push_back({ 6, 7 });
    edges.push_back({ 10, 6 });
    edges.push_back({ 5, 6 });

    for (int i = 0; i < edges.size(); i++)
    {
        edges[i].first -= 1;
        edges[i].second -= 1;
    }

    // now build node-edge mapping
    GridGeom::Mesh mesh;
    mesh.Set(edges, nodes, GridGeom::Projections::cartesian);

    GridGeom::Polygons polygons;
    const std::vector<GridGeom::Point> polygon;
    std::vector<GridGeom::Point> meshBoundaryPolygon;
    int numNodesBoundaryPolygons;
    polygons.Set(polygon, GridGeom::Projections::cartesian);
    polygons.MeshBoundaryToPolygon(mesh, 0, 1, meshBoundaryPolygon, numNodesBoundaryPolygons);


    //constexpr double tolerance = 1e-2;

    //ASSERT_NEAR(325.590101919525, mesh.m_nodes[0].x, tolerance);
    //ASSERT_NEAR(229.213730481198, mesh.m_nodes[1].x, tolerance);
    //ASSERT_NEAR(263.439319753147, mesh.m_nodes[2].x, tolerance);
    //ASSERT_NEAR(429.191105834504, mesh.m_nodes[3].x, tolerance);
    //ASSERT_NEAR(535.865215426468, mesh.m_nodes[4].x, tolerance);
    //ASSERT_NEAR(503.753784179688, mesh.m_nodes[5].x, tolerance);
    //ASSERT_NEAR(354.048340705929, mesh.m_nodes[6].x, tolerance);
    //ASSERT_NEAR(346.790050854504, mesh.m_nodes[7].x, tolerance);
    //ASSERT_NEAR(315.030130405285, mesh.m_nodes[8].x, tolerance);
    //ASSERT_NEAR(424.314957449766, mesh.m_nodes[9].x, tolerance);

    //ASSERT_NEAR(455.319334078551, mesh.m_nodes[0].y, tolerance);
    //ASSERT_NEAR(362.573521507281, mesh.m_nodes[1].y, tolerance);
    //ASSERT_NEAR(241.096458631763, mesh.m_nodes[2].y, tolerance);
    //ASSERT_NEAR(211.483073921775, mesh.m_nodes[3].y, tolerance);
    //ASSERT_NEAR(311.401495506714, mesh.m_nodes[4].y, tolerance);
    //ASSERT_NEAR(432.379974365234, mesh.m_nodes[5].y, tolerance);
    //ASSERT_NEAR(458.064836627594, mesh.m_nodes[6].y, tolerance);
    //ASSERT_NEAR(405.311585650679, mesh.m_nodes[7].y, tolerance);
    //ASSERT_NEAR(319.612138503550, mesh.m_nodes[8].y, tolerance);
    //ASSERT_NEAR(327.102805172725, mesh.m_nodes[9].y, tolerance);
}

TEST(Polygons, CreatePointsInPolygons)
{
    // Prepare
    GridGeom::Polygons polygons;
    std::vector<GridGeom::Point> nodes;

    nodes.push_back(GridGeom::Point{ 302.002502,472.130371 });
    nodes.push_back(GridGeom::Point{ 144.501526, 253.128174 });
    nodes.push_back(GridGeom::Point{ 368.752930, 112.876755 });
    nodes.push_back(GridGeom::Point{ 707.755005, 358.879242 });
    nodes.push_back(GridGeom::Point{ 301.252502, 471.380371 });
    nodes.push_back(GridGeom::Point{ 302.002502, 472.130371 });

    polygons.Set(nodes, GridGeom::Projections::cartesian);

    // Execute
    std::vector<std::vector<GridGeom::Point>> generatedPoints;
    bool success = polygons.CreatePointsInPolygons(generatedPoints);
    ASSERT_TRUE(success);

    // Assert
    constexpr double tolerance = 1e-5;

    ASSERT_NEAR(302.00250199999999, generatedPoints[0][0].x, tolerance);
    ASSERT_NEAR(472.13037100000003, generatedPoints[0][0].y, tolerance);

    ASSERT_NEAR(144.50152600000001, generatedPoints[0][1].x, tolerance);
    ASSERT_NEAR(253.12817400000000, generatedPoints[0][1].y, tolerance);

    ASSERT_NEAR(368.75292999999999, generatedPoints[0][2].x, tolerance);
    ASSERT_NEAR(112.87675500000000, generatedPoints[0][2].y, tolerance);

    ASSERT_NEAR(707.75500499999998, generatedPoints[0][3].x, tolerance);
    ASSERT_NEAR(358.87924199999998, generatedPoints[0][3].y, tolerance);

    ASSERT_NEAR(301.25250199999999, generatedPoints[0][4].x, tolerance);
    ASSERT_NEAR(471.38037100000003, generatedPoints[0][4].y, tolerance);
}

TEST(Polygons, RefinePolygon)
{
    // Prepare
    GridGeom::Polygons polygons;
    std::vector<GridGeom::Point> nodes;

    nodes.push_back(GridGeom::Point{ 0,0 });
    nodes.push_back(GridGeom::Point{ 3, 0 });
    nodes.push_back(GridGeom::Point{ 3, 3 });
    nodes.push_back(GridGeom::Point{ 0, 3 });
    nodes.push_back(GridGeom::Point{ 0, 0 });

    polygons.Set(nodes, GridGeom::Projections::cartesian);

    // Execute
    std::vector<std::vector<GridGeom::Point>> generatedPoints;
    std::vector<GridGeom::Point> refinedPolygon;
    bool success = polygons.RefinePart(0, 0, 1.0, refinedPolygon);
    ASSERT_TRUE(success);

    ASSERT_EQ(13, refinedPolygon.size());
    constexpr double tolerance = 1e-5;

    ASSERT_NEAR(0, refinedPolygon[0].x, tolerance);
    ASSERT_NEAR(1, refinedPolygon[1].x, tolerance);
    ASSERT_NEAR(2, refinedPolygon[2].x, tolerance);
    ASSERT_NEAR(3, refinedPolygon[3].x, tolerance);
    ASSERT_NEAR(3, refinedPolygon[4].x, tolerance);
    ASSERT_NEAR(3, refinedPolygon[5].x, tolerance);
    ASSERT_NEAR(3, refinedPolygon[6].x, tolerance);
    ASSERT_NEAR(2, refinedPolygon[7].x, tolerance);
    ASSERT_NEAR(1, refinedPolygon[8].x, tolerance);
    ASSERT_NEAR(0, refinedPolygon[9].x, tolerance);
    ASSERT_NEAR(0, refinedPolygon[10].x, tolerance);
    ASSERT_NEAR(0, refinedPolygon[11].x, tolerance);
    ASSERT_NEAR(0, refinedPolygon[12].x, tolerance);

    ASSERT_NEAR(0, refinedPolygon[0].y, tolerance);
    ASSERT_NEAR(0, refinedPolygon[1].y, tolerance);
    ASSERT_NEAR(0, refinedPolygon[2].y, tolerance);
    ASSERT_NEAR(0, refinedPolygon[3].y, tolerance);
    ASSERT_NEAR(1, refinedPolygon[4].y, tolerance);
    ASSERT_NEAR(2, refinedPolygon[5].y, tolerance);
    ASSERT_NEAR(3, refinedPolygon[6].y, tolerance);
    ASSERT_NEAR(3, refinedPolygon[7].y, tolerance);
    ASSERT_NEAR(3, refinedPolygon[8].y, tolerance);
    ASSERT_NEAR(3, refinedPolygon[9].y, tolerance);
    ASSERT_NEAR(2, refinedPolygon[10].y, tolerance);
    ASSERT_NEAR(1, refinedPolygon[11].y, tolerance);
    ASSERT_NEAR(0, refinedPolygon[12].y, tolerance);
}

TEST(Polygons, RefinePolygonLongerSquare)
{
    // Prepare
    GridGeom::Polygons polygons;
    std::vector<GridGeom::Point> nodes;

    nodes.push_back(GridGeom::Point{ 0,0 });
    nodes.push_back(GridGeom::Point{ 3, 0 });
    nodes.push_back(GridGeom::Point{ 3, 3 });
    nodes.push_back(GridGeom::Point{ 0, 3.5 });
    nodes.push_back(GridGeom::Point{ 0, 0 });

    polygons.Set(nodes, GridGeom::Projections::cartesian);

    // Execute
    std::vector<std::vector<GridGeom::Point>> generatedPoints;
    std::vector<GridGeom::Point> refinedPolygon;
    bool success = polygons.RefinePart(0, 0, 1.0, refinedPolygon);
    ASSERT_TRUE(success);

    ASSERT_EQ(13, refinedPolygon.size());
    constexpr double tolerance = 1e-5;

    ASSERT_NEAR(0, refinedPolygon[0].x, tolerance);
    ASSERT_NEAR(1, refinedPolygon[1].x, tolerance);
    ASSERT_NEAR(2, refinedPolygon[2].x, tolerance);
    ASSERT_NEAR(3, refinedPolygon[3].x, tolerance);
    ASSERT_NEAR(3, refinedPolygon[4].x, tolerance);
    ASSERT_NEAR(3, refinedPolygon[5].x, tolerance);
    ASSERT_NEAR(3, refinedPolygon[6].x, tolerance);
    ASSERT_NEAR(2.0136060761678563, refinedPolygon[7].x, tolerance);
    ASSERT_NEAR(1.0272121523357125, refinedPolygon[8].x, tolerance);
    ASSERT_NEAR(0.040818228503568754, refinedPolygon[9].x, tolerance);
    ASSERT_NEAR(0, refinedPolygon[10].x, tolerance);
    ASSERT_NEAR(0, refinedPolygon[11].x, tolerance);
    ASSERT_NEAR(0, refinedPolygon[12].x, tolerance);

    ASSERT_NEAR(0, refinedPolygon[0].y, tolerance);
    ASSERT_NEAR(0, refinedPolygon[1].y, tolerance);
    ASSERT_NEAR(0, refinedPolygon[2].y, tolerance);
    ASSERT_NEAR(0, refinedPolygon[3].y, tolerance);
    ASSERT_NEAR(1, refinedPolygon[4].y, tolerance);
    ASSERT_NEAR(2, refinedPolygon[5].y, tolerance);
    ASSERT_NEAR(3, refinedPolygon[6].y, tolerance);
    ASSERT_NEAR(3.1643989873053573, refinedPolygon[7].y, tolerance);
    ASSERT_NEAR(3.3287979746107146, refinedPolygon[8].y, tolerance);
    ASSERT_NEAR(3.4931969619160719, refinedPolygon[9].y, tolerance);
    ASSERT_NEAR(2.5413812651491092, refinedPolygon[10].y, tolerance);
    ASSERT_NEAR(1.5413812651491092, refinedPolygon[11].y, tolerance);
    ASSERT_NEAR(0.54138126514910923, refinedPolygon[12].y, tolerance);
}