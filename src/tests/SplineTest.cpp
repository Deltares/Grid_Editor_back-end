#include "GridGeomTest.hpp"
#include "../Splines.hpp"
#include "../Entities.hpp"

TEST(SplineTests, SetSpline)
{
    //One gets the edges
    std::vector<GridGeom::Point> splineNodes;
    splineNodes.push_back(GridGeom::Point{GridGeom::doubleMissingValue, GridGeom::doubleMissingValue });
    splineNodes.push_back(GridGeom::Point{ GridGeom::doubleMissingValue, GridGeom::doubleMissingValue });
    splineNodes.push_back(GridGeom::Point{ GridGeom::doubleMissingValue, GridGeom::doubleMissingValue });
    splineNodes.push_back(GridGeom::Point{ GridGeom::doubleMissingValue, GridGeom::doubleMissingValue });

    GridGeom::Splines splines;
    bool success = splines.Set(splineNodes);

    ASSERT_TRUE(success);
    EXPECT_EQ(splines.m_numSplines, 1);
    EXPECT_EQ(splines.m_splines[0].size(), splineNodes.size());
    EXPECT_EQ(splines.m_numAllocatedSplines, 5);
    EXPECT_EQ(splines.m_numAllocatedSplineNodes[0], 10);
}

TEST(SplineTests, CubicSplineInterpolation)
{
    //One gets the edges
    std::vector<GridGeom::Point> splineNodes;
    
    splineNodes.push_back(GridGeom::Point{ 212.001953125000, 155.627197265625 });
    splineNodes.push_back(GridGeom::Point{ 529.253906250000, 432.379974365234 });
    splineNodes.push_back(GridGeom::Point{ 930.506469726562, 453.380187988281 });

    int pointsBetweenVertices = 20;
    std::vector<GridGeom::Point> coordinatesDerivatives(splineNodes.size());
    GridGeom::Splines::SecondOrderDerivative(splineNodes, splineNodes.size(), coordinatesDerivatives);
    std::vector<GridGeom::Point> splineCoordinates;

    for (int n = 0; n < splineNodes.size() - 1; n++)
    {
        for (int p = 0; p <= pointsBetweenVertices; p++)
        {
            const double pointAdimensionalCoordinate = n + double(p) / double(pointsBetweenVertices);
            GridGeom::Point pointCoordinate;
            GridGeom::Splines::Interpolate(splineNodes, coordinatesDerivatives, pointAdimensionalCoordinate, pointCoordinate);
            splineCoordinates.push_back({ pointCoordinate.x, pointCoordinate.y });
        }
    }

    const double tolerance = 1e-3;
    ASSERT_NEAR(226.817168170929, splineCoordinates[1].x, tolerance);
    ASSERT_NEAR(241.648133331299, splineCoordinates[2].x, tolerance);
    ASSERT_NEAR(256.510598720551, splineCoordinates[3].x, tolerance);
    ASSERT_NEAR(271.420314453125, splineCoordinates[4].x, tolerance);
    ASSERT_NEAR(286.393030643463, splineCoordinates[5].x, tolerance);
    ASSERT_NEAR(930.506469726562, splineCoordinates.back().x, tolerance);

    ASSERT_NEAR(172.653750896454, splineCoordinates[1].y, tolerance);
    ASSERT_NEAR(189.632350921631, splineCoordinates[2].y, tolerance);
    ASSERT_NEAR(206.515043735504, splineCoordinates[3].y, tolerance);
    ASSERT_NEAR(223.253875732422, splineCoordinates[4].y, tolerance);
    ASSERT_NEAR(453.380187988281, splineCoordinates.back().y, tolerance);
}

TEST(SplineTests, SplineIntersection)
{
    std::vector<GridGeom::Point> firstSpline;
    firstSpline.push_back(GridGeom::Point{ 152.001571655273, 86.6264953613281 });
    firstSpline.push_back(GridGeom::Point{ 374.752960205078, 336.378997802734 });
    firstSpline.push_back(GridGeom::Point{ 850.255920410156, 499.130676269531 });
    GridGeom::Splines splines;
    bool success = splines.Set(firstSpline);

    std::vector<GridGeom::Point> secondSpline;
    secondSpline.push_back(GridGeom::Point{ 72.5010681152344,391.129577636719 });
    secondSpline.push_back(GridGeom::Point{ 462.503479003906, 90.3765411376953 });
    success = splines.Set(secondSpline);

    double crossProductIntersection;
    GridGeom::Point dimensionalIntersection;
    GridGeom::Point adimensionalIntersection;

    splines.GetSplinesIntersection(0, 1, GridGeom::Projections::cartesian,
        crossProductIntersection, dimensionalIntersection, adimensionalIntersection);

    const double tolerance = 1e-5;
    ASSERT_NEAR(261.736770097059, dimensionalIntersection.x, tolerance);
    ASSERT_NEAR(245.199166962145, dimensionalIntersection.y, tolerance);
    ASSERT_NEAR(0.601498208554790, adimensionalIntersection.x, tolerance);
    ASSERT_NEAR(0.485216749175026, adimensionalIntersection.y, tolerance);
    ASSERT_NEAR(-0.996215079635043, crossProductIntersection, tolerance);

}