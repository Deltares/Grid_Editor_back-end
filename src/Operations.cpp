#pragma once

#include <cmath>
#include "Entities.hpp"
#include "Constants.cpp"

#include <boost/geometry/geometries/segment.hpp> 
#include <boost/geometry/geometries/register/point.hpp>
#include <boost/geometry/geometries/adapted/boost_tuple.hpp>

#include <boost/geometry.hpp>
// register node so we can use boost geometry algorithms
BOOST_GEOMETRY_REGISTER_POINT_2D(GridGeom::Point, double, boost::geometry::cs::cartesian, x, y);

namespace GridGeom
{

    // functions that depends of the point type
    template<OperationTypes operationType>
    struct Operations;

    // coordinate reference indipendent operations
    template<typename T>
    static double dotProduct(const T& dx1, const T& dx2)
    {
        return dx1 * dx2;
    }

    template<typename T, typename... Args>
    static T dotProduct(T dx1, T dx2, Args... args) 
    {
        return dx1 * dx2 + dotProduct(args...);
    }

    // transform 2d spherical to 3d cartesian
    static void sphericalToCartesian(const Point& sphericalPoint, cartesian3DPoint& cartesianPoint)
    {
        cartesianPoint.z = earth_radius * sin(sphericalPoint.y * degrad_hp);
        double rr = earth_radius * cos(sphericalPoint.y * degrad_hp);
        cartesianPoint.x = rr * cos(sphericalPoint.x * degrad_hp);
        cartesianPoint.y = rr * sin(sphericalPoint.x * degrad_hp);
    }

    //  transform 3d cartesian coordinates to 2d spherical
    static void cartesianToSpherical(const cartesian3DPoint& cartesianPoint, const double referenceLongitude, Point& sphericalPoint)
    {
        double angle = atan2(cartesianPoint.y, cartesianPoint.x) * raddeg_hp;
        sphericalPoint.y = atan2(cartesianPoint.z, sqrt(cartesianPoint.x * cartesianPoint.x + cartesianPoint.y * cartesianPoint.y)) * raddeg_hp;
        sphericalPoint.x = angle + std::lround((referenceLongitude - angle) / 360.0) * 360.0;
    }

    // isLeft(): tests if a point is Left|On|Right of an infinite line.
    //    Input:  three points leftPoint, rightPoint, and point
    //    Return: >0 for point left of the line through leftPoint and rightPoint
    //            =0 for point  on the line
    //            <0 for point  right of the line
    static double isLeft(const Point& leftPoint, const Point& rightPoint, const Point& point)
    {
        double left = (rightPoint.x - leftPoint.x) * (point.y - leftPoint.y) - (point.x - leftPoint.x) * (rightPoint.y - leftPoint.y);
        return left;
    }

    // check if a point is in polygon using the winding number method
    // polygon: vector of points in counter clockwise order
    static bool pointInPolygon(const Point& point, const std::vector<Point>& polygon, const int numberOfPolygonPoints)
    {
        int windingNumber = 0;
        for (int n = 0; n < numberOfPolygonPoints; n++)
        {
            if (polygon[n].y <= point.y) // an upward crossing
            {
                if (polygon[n + 1].y > point.y)
                {
                    if (isLeft(polygon[n], polygon[n + 1], point) > 0.0)
                    {
                        ++windingNumber; // have  a valid up intersect
                    }
                }
            }
            else
            {
                if (polygon[n + 1].y <= point.y) // a downward crossing
                {
                    if (isLeft(polygon[n], polygon[n + 1], point) < 0.0)
                    {
                        --windingNumber; // have  a valid down intersect
                    }
                }
            }
        }
        return windingNumber == 0 ? false : true;
    }


    //faceAreaAndCenterOfMass: for cartesian, spherical point and spherical3dPoint
    template<OperationTypes operationType>
    static bool faceAreaAndCenterOfMass(std::vector<Point>& polygon, const int numberOfPolygonPoints, double& area, Point& centerOfMass)
    {
        if (numberOfPolygonPoints < 1)
        {
            return false;
        }
        
        double minX = std::numeric_limits<double>::max();
        double minY = std::numeric_limits<double>::max();
        Operations<operationType>::referencePoint(polygon, numberOfPolygonPoints, minX, minY);

        Point reference{ minX, minY };
        area = 0.0;
        double xCenterOfMass = 0.0;
        double yCenterOfMass = 0.0;
        for (int p = 0; p < numberOfPolygonPoints; p++)
        {
            double dx0 = Operations<operationType>::getDx(reference, polygon[p]);
            double dy0 = Operations<operationType>::getDy(reference, polygon[p]);
            double dx1 = Operations<operationType>::getDx(reference, polygon[p + 1]);
            double dy1 = Operations<operationType>::getDy(reference, polygon[p + 1]);

            double xc = 0.5 * (dx0 + dx1);
            double yc = 0.5 * (dy0 + dy1);

            dx0 = Operations<operationType>::getDx(polygon[p], polygon[p + 1]);
            dy0 = Operations<operationType>::getDy(polygon[p], polygon[p + 1]);
            double dsx = dy0;
            double dsy = -dx0;
            double xds = xc * dsx + yc * dsy;
            area = area + 0.5 * xds;

            xCenterOfMass = xCenterOfMass + xds * xc;
            yCenterOfMass = yCenterOfMass + xds * yc;
        }

        double fac = 1.0 / (3.0 * area);
        xCenterOfMass = fac * xCenterOfMass;
        yCenterOfMass = fac * yCenterOfMass;

        if constexpr (operationType == sphericalOperations)
        {
            yCenterOfMass = yCenterOfMass / (earth_radius * degrad_hp);
            xCenterOfMass = xCenterOfMass / (earth_radius * degrad_hp * cos((yCenterOfMass + minY) * degrad_hp));
        }

        centerOfMass.x = xCenterOfMass + minX;
        centerOfMass.y = yCenterOfMass + minY;

        area = std::abs(area);

        return true;
    }

    static bool lineCrossing(const Point& firstSegmentFistPoint, const Point& firstSegmentSecondPoint, const Point& secondSegmentFistPoint, const Point& secondSegmentSecondPoint, Point& intersection)
    {
        typedef boost::geometry::model::segment<Point> Segment;
        Segment firstSegment(firstSegmentFistPoint, firstSegmentSecondPoint);
        Segment secondSegment(secondSegmentFistPoint, secondSegmentSecondPoint);

        std::vector<Point> intersections;
        boost::geometry::intersection(firstSegment, secondSegment, intersections);

        if (!intersections.empty())
        {
            intersection = intersections[0];
            return true;
        }
        return false;
    }

    // cartesian points
    template <>
    struct Operations<cartesianOperations>
    {
        //normalout, Creates the relative unit normal vector to edge 1->2
        static void normalVector(const Point& firstPoint, const Point& secondPoint, const Point& insidePoint, Point& result)
        {
            double dx = getDx(firstPoint, secondPoint);
            double dy = getDy(firstPoint, secondPoint);
            const double squaredDistance = dx * dx + dy * dy;
            if (squaredDistance != 0.0)
            {
                const double distance = sqrt(squaredDistance);
                result.x = dy / distance;
                result.y = -dx / distance;
            }
        }

        //normaloutchk
        static void normalVectorInside(const Point& firstPoint, const Point& secondPoint, const Point& insidePoint, Point& result, bool& flippedNormal)
        {
            normalVector(firstPoint, secondPoint, insidePoint, result);
            flippedNormal = false;
            Point thirdPoint{ firstPoint.x + result.x, firstPoint.y + result.y };

            if (outerProductTwoSegments(firstPoint, thirdPoint, firstPoint, secondPoint) * outerProductTwoSegments(firstPoint, insidePoint, firstPoint, secondPoint) > 0.0)
            {
                result.x = -result.x;
                result.y = -result.y;
                flippedNormal = true;
            }
            else
            {
                flippedNormal = false;
            }
        }

        static double getDx(const Point& firstPoint, const Point& secondPoint)
        {
            return secondPoint.x - firstPoint.x;
        }

        static double getDy(const Point& firstPoint, const Point& secondPoint)
        {
            return secondPoint.y - firstPoint.y;
        }

        static void add(Point& point, const Point& normal, const double increment)
        {
            point.x = point.x + normal.x * increment;
            point.y = point.y + normal.y * increment;
        }

        static void referencePoint(std::vector<Point>& polygon, const int numPoints, double& minX, double& minY)
        {
            minX = std::numeric_limits<double>::max();
            minY = std::numeric_limits<double>::max();
            for (int i=0; i< numPoints;i++)
            {
                if (polygon[i].x < minX)
                {
                    minX = polygon[i].x;
                }
                if (abs(polygon[i].y) < abs(minY))
                {
                    minY = polygon[i].y;
                }
            }
        }

        //dbdistance
        static double distance(const Point& firstPoint, const Point& secondPoint)
        {
            double dx = getDx(firstPoint, secondPoint);
            double dy = getDy(firstPoint, secondPoint);
            const double squaredDistance = dx * dx + dy * dy;
            double distance = 0.0;
            if (squaredDistance != 0.0)
            {
                distance = sqrt(squaredDistance);
            }
            return distance;
        }

        //dLINEDIS3
        static double distanceFromLine(const Point& p3, const Point& p1, const Point& p2, Point& normalPoint, double& rlout)
        {
            double dis = 0.0;
            double r2 = distance(p2, p1);
            if (r2 != 0.0)
            {
                double rl = (getDx(p1, p3) * getDx(p1, p2) + getDy(p1, p3) * getDy(p1, p2)) / (r2 * r2);
                rlout = std::max(std::min(1.0, rl), 0.0);
                normalPoint.x = p1.x + rlout * (p2.x - p1.x);
                normalPoint.y = p1.y + rlout * (p2.y - p1.y);
                dis = distance(p3, normalPoint);
            }
            return dis;
        }

        //out product of two segments
        static double outerProductTwoSegments(const Point& firstPointFirstSegment, const Point& secondPointFirstSegment, const Point& firstPointSecondSegment, const Point& secondPointSecondSegment)
        {
            double dx1 = getDx(firstPointFirstSegment, secondPointFirstSegment);
            double dx2 = getDx(firstPointSecondSegment, secondPointSecondSegment);

            double dy1 = getDy(firstPointFirstSegment, secondPointFirstSegment);
            double dy2 = getDy(firstPointSecondSegment, secondPointSecondSegment);
         
            return dx1 * dy2 - dy1 * dx2;
        }

        //inner product of two segments
        static double innerProductTwoSegments(const Point& firstPointFirstSegment, const Point& secondPointFirstSegment, const Point& firstPointSecondSegment, const Point& secondPointSecondSegment)
        {
            double dx1 = getDx(firstPointFirstSegment, secondPointFirstSegment);
            double dx2 = getDx(firstPointSecondSegment, secondPointSecondSegment);

            double dy1 = getDy(firstPointFirstSegment, secondPointFirstSegment);
            double dy2 = getDy(firstPointSecondSegment, secondPointSecondSegment);

            return dx1 * dx2 + dy1 * dy2;
        }


        static bool orthogonalizationComputeLocalCoordinates(const std::vector<size_t>& m_nodesNumEdges, const std::vector<size_t>& numConnectedNodes, std::vector<int>& localCoordinates)
        {
            //do nothing
            return true;
        }


        static inline bool orthogonalizationComputeJacobian(const int currentNode, const std::vector<double>& Jxi, const std::vector<double>& Jeta, const std::vector<size_t>& connectedNodes, const int numNodes, const std::vector<Point>& nodes, std::vector<double>& J)
        {
            J[0] = 0.0; 
            J[1] = 0.0;
            J[2] = 0.0;
            J[3] = 0.0;
            for (int i = 0; i < numNodes; i++)
            {
                J[0] += Jxi[i] * nodes[connectedNodes[i]].x;
                J[1] += Jxi[i] * nodes[connectedNodes[i]].y;
                J[2] += Jeta[i] * nodes[connectedNodes[i]].x;
                J[3] += Jeta[i] * nodes[connectedNodes[i]].y;
            }
            return true;
        }

        static  bool orthogonalizationComputeDeltas(int firstNode, int secondNode, double wwx, double wwy, const std::vector<Point>& nodes, double& dx0, double& dy0, std::vector<double>& increments)
        {

            increments[0] += wwx;
            increments[1] += wwy;

            dx0 = dx0 + wwx * (nodes[firstNode].x - nodes[secondNode].x);
            dy0 = dy0 + wwy * (nodes[firstNode].y - nodes[secondNode].y);
            return true;
        }

        static bool orthogonalizationComputeCoordinates(double dx0, double dy0, const Point& point, Point& updatedPoint)
        {
            double x0 = point.x + dx0;
            double y0 = point.y + dy0;
            static constexpr double relaxationFactorCoordinates = 1.0 - relaxationFactorOrthogonalizationUpdate;

            updatedPoint.x = relaxationFactorOrthogonalizationUpdate * x0 + relaxationFactorCoordinates * point.x;
            updatedPoint.y = relaxationFactorOrthogonalizationUpdate * y0 + relaxationFactorCoordinates * point.y;

            return true;
        }

        static bool circumcenterOfTriangle(const Point& p1, const Point& p2, const Point& p3, Point& circumcenter)
        {
            double dx2 = getDx(p1, p2);
            double dy2 = getDy(p1, p2);

            double dx3 = getDx(p1, p3);
            double dy3 = getDy(p1, p3);

            double den = dy2 * dx3 - dy3 * dx2;
            double z = 0.0;
            if (den != 0.0)
            {
                z = (dx2 * (dx2 - dx3) + dy2 * (dy2 - dy3)) / den;
            }

            circumcenter.x = p1.x + 0.5 * (dx3 - z * dy3);
            circumcenter.y = p1.y + 0.5 * (dy3 + z * dx3);
            return true;
        }

    };
    
    // spherical point
    template <>
    struct Operations<sphericalOperations>
    {
        static void normalVector(const Point& firstPoint, const Point& secondPoint, const Point& insidePoint, Point& result)
        {
            cartesian3DPoint firstPointCartesianCoordinates;
            cartesian3DPoint secondPointCartesianCoordinates;
            sphericalToCartesian(firstPoint, firstPointCartesianCoordinates);
            sphericalToCartesian(secondPoint, secondPointCartesianCoordinates);

            double lambda = insidePoint.x * degrad_hp;
            double phi = insidePoint.y * degrad_hp;
            double elambda[3] = { -sin(lambda), cos(lambda), 0.0 };
            double ephi[3] = { -sin(lambda), cos(lambda), 0.0 };

            double dx = (secondPointCartesianCoordinates.x - firstPointCartesianCoordinates.x) * elambda[0] +
                (secondPointCartesianCoordinates.y - firstPointCartesianCoordinates.y) * elambda[1] +
                (secondPointCartesianCoordinates.z - firstPointCartesianCoordinates.z) * elambda[2];

            double dy = (secondPointCartesianCoordinates.x - firstPointCartesianCoordinates.x) * ephi[0] +
                (secondPointCartesianCoordinates.y - firstPointCartesianCoordinates.y) * ephi[1] +
                (secondPointCartesianCoordinates.z - firstPointCartesianCoordinates.z) * ephi[2];

            double squaredDistance = dx * dx + dy * dy;
            if (squaredDistance != 0.0)
            {
                const double distance = sqrt(squaredDistance);
                result.x = dx / distance;
                result.y = dy / distance;
            }
        }

        static void normalVectorInside(const Point& firstPoint, const Point& secondPoint, const Point& insidePoint, Point& result, bool& flippedNormal)
        {

            //if (JSFERIC.eq.1 . and .jasfer3D.eq.0) xn = xn * cos(dg2rd * 0.5d0 * (y0 + y1)) !normal vector needs to be in Cartesian coordinates
        }

        static double getDx(const Point& firstPoint, const Point& secondPoint)
        {
            double  firstPointYDiff = abs(abs(firstPoint.y) - 90.0);
            double  secondPointYDiff = abs(abs(secondPoint.y) - 90.0);
            if (firstPointYDiff <= dtol_pole && secondPointYDiff > dtol_pole || firstPointYDiff > dtol_pole && secondPointYDiff <= dtol_pole)
            {
                return 0.0;
            }
            double firstPointX = firstPoint.x;
            double secondPointX = secondPoint.x;
            if (firstPointX - secondPointX > 180.0)
            {
                firstPointX -= 360.0;
            }
            else if (firstPointX - secondPointX < -180.0)
            {
                firstPointX += 360.0;
            }

            firstPointX = firstPointX * degrad_hp;
            secondPointX = secondPointX * degrad_hp;
            double firstPointY = firstPoint.y * degrad_hp;
            double secondPointY = secondPoint.y * degrad_hp;
            double cosPhi = cos(0.5 * (firstPointY + secondPointY));
            double dx = earth_radius * cosPhi * (secondPointX - firstPointX);
            return dx;
        }

        static double getDy(const Point& firstPoint, const Point& secondPoint)
        {
            double firstPointY = firstPoint.y * degrad_hp;
            double secondPointY = secondPoint.y * degrad_hp;
            double dy = earth_radius * (secondPointY - firstPointY);
            return dy;
        }

        static void add(Point& point, const Point& normal, const double increment)
        {
            double convertedIncrement = raddeg_hp * increment / earth_radius;
            double xf = 1.0 / cos(degrad_hp * point.y);
            point.x = point.x + normal.x * convertedIncrement * xf;
            point.y = point.y + normal.y * convertedIncrement;
        }

        static void referencePoint(std::vector<Point>& polygon, const int numPoints, double& minX, double& minY)
        {
            minX = std::numeric_limits<double>::max();
            minY = std::numeric_limits<double>::max();
            double xmax = std::numeric_limits<double>::min();
            for (int i=0;i< numPoints;i++)
            {
                if (polygon[i].x < minX)
                {
                    minX =polygon[i].x;
                }
                if (abs(polygon[i].y) < abs(minY))
                {
                    minY = polygon[i].y;
                }
                if (polygon[i].x> xmax)
                {
                    xmax = polygon[i].x;
                }
            }

            if (xmax - minX > 180.0)
            {
                double deltaX = xmax - 180.0;
                for (int i = 0; i < numPoints; i++)
                {
                    if(polygon[i].x < deltaX)
                    {
                        polygon[i].x = polygon[i].x + 360.0;
                    } 
                }
                minX = minX + 360.0;
            }
            //TODO: check result
            minX = std::min_element(polygon.begin(), polygon.end(), [](const Point& p1, const Point& p2) { return p1.x < p2.x; })->x;
        }

        static double distance(const Point& firstPoint, const Point& secondPoint)
        {
            return -1.0;
        }

        //dLINEDIS3
        static double distanceFromLine(const Point& p3, const Point& p1, const Point& p2, Point& normalPoint, double& rlout)
        {
            //TODO: implement me
            return -1.0;
        }


        //out product of two segments
        static double outerProductTwoSegments(const Point& firstPointFirstSegment, const Point& secondPointFirstSegment, const Point& firstPointSecondSegment, const Point& secondPointSecondSegment)
        {
            //TODO: IMPLEMENTATION IS MISSING

            return 0.0;
        }

        static double innerProductTwoSegments(const Point& firstPointFirstSegment, const Point& secondPointFirstSegment, const Point& firstPointSecondSegment, const Point& secondPointSecondSegment)
        {
            cartesian3DPoint firstPointFirstSegment3D;
            cartesian3DPoint secondPointFirstSegment3D;
            cartesian3DPoint firstPointSecondSegment3D;
            cartesian3DPoint secondPointSecondSegment3D;

            sphericalToCartesian(firstPointFirstSegment, firstPointFirstSegment3D);
            sphericalToCartesian(secondPointFirstSegment, secondPointFirstSegment3D);
            sphericalToCartesian(firstPointSecondSegment, firstPointSecondSegment3D);
            sphericalToCartesian(secondPointSecondSegment, secondPointSecondSegment3D);

            double dx1 = secondPointFirstSegment3D.x - firstPointFirstSegment3D.x;
            double dy1 = secondPointFirstSegment3D.y - firstPointFirstSegment3D.y;
            double dz1 = secondPointFirstSegment3D.z - firstPointFirstSegment3D.z;

            double dx2 = secondPointSecondSegment3D.x - firstPointSecondSegment3D.x;
            double dy2 = secondPointSecondSegment3D.y - firstPointSecondSegment3D.y;
            double dz2 = secondPointSecondSegment3D.z - firstPointSecondSegment3D.z;

            return dotProduct(dx1, dx2, dy1, dy2, dz1, dz2);
        }

        //TODO:: comp_local_coords
        static bool orthogonalizationComputeLocalCoordinates(const std::vector<size_t>& m_nodesNumEdges, const std::vector<size_t>& numConnectedNodes, std::vector<int>& localCoordinates)
        {
            localCoordinates.resize(m_nodesNumEdges.size(), 0);
            localCoordinates[0]=1;
            for (int i = 0; i < m_nodesNumEdges.size(); i++)
            {
                localCoordinates[i+1]= localCoordinates[i]+std::max(m_nodesNumEdges[i]+1, numConnectedNodes[i]);
            }
            return true;
        }


        static bool orthogonalizationComputeJacobian(const int currentNode, const std::vector<double>& Jxi, const std::vector<double>& Jeta, const std::vector<size_t>& connectedNodes, const int numNodes, const std::vector<Point>& nodes, std::vector<double>& J)
        {
            double factor = std::cos(nodes[currentNode].y) * degrad_hp;
            J[0] = 0.0;
            J[1] = 0.0;
            J[2] = 0.0;
            J[3] = 0.0;
            for (int i = 0; i < numNodes; i++)
            {
                J[0] += Jxi[i] * nodes[connectedNodes[i]].x;
                J[1] += Jxi[i] * nodes[connectedNodes[i]].y;
                J[2] += Jeta[i] * nodes[connectedNodes[i]].x;
                J[3] += Jeta[i] * nodes[connectedNodes[i]].y;
            }
            return true;
        }

        static bool orthogonalizationComputeDeltas(int firstNode, int secondNode, double wwx, double wwy, const std::vector<Point>& nodes, double& dx0, double& dy0, std::vector<double>& increments)
        {
            double wwxTransformed = wwx * earth_radius * degrad_hp;
            double wwyTransformed = wwy * earth_radius * degrad_hp;

            increments[0] += wwxTransformed;
            increments[1] += wwyTransformed;

            dx0 = dx0 + wwxTransformed * (nodes[firstNode].x - nodes[secondNode].x);
            dy0 = dy0 + wwxTransformed * (nodes[firstNode].y - nodes[secondNode].y);

            return true;
        }

        static bool orthogonalizationComputeCoordinates(double dx0, double dy0, const Point& point, Point& updatedPoint)
        {
            //TODO: implement
            //if (jsferic.eq.1 . and .jasfer3D.eq.1) then
            //    dumx(1) = relaxin * Dx0
            //    dumy(1) = relaxin * Dy0
            //    call loc2spher(xk(k), yk(k), 1, dumx, dumy, xk1(k), yk1(k))
            //else

            return true;
        }

        static bool circumcenterOfTriangle(const Point& p1, const Point& p2, const Point& p3, Point& circumcenter)
        {
            double dx2 = getDx(p1, p2);
            double dy2 = getDy(p1, p2);

            double dx3 = getDx(p1, p3);
            double dy3 = getDy(p1, p3);

            double den = dy2 * dx3 - dy3 * dx2;
            double z = 0.0;
            if (den > 1e-16)
            {
                z = (dx2 * (dx2 - dx3) + dy2 * (dy2 - dy3)) / den;
            }

            //TODO circumcenter3 FINISH
            //phi = (y(1) + y(2) + y(3)) / 3d0
            //    xf = 1d0 / dcos(degrad_hp * phi)
            //    xz = x(1) + xf * 0.5d0 * (dx3 - z * dy3) * raddeg_hp / earth_radius
            //    yz = y(1) + 0.5d0 * (dy3 + z * dx3) * raddeg_hp / earth_radius
            return true;
        }
    };  

}