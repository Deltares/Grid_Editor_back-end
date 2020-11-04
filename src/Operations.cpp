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

#include <cmath>
#include <algorithm>
#include <numeric>
#include "Entities.hpp"
#include "Constants.cpp"
#include "SpatialTrees.hpp"
#include "Exceptions.hpp"

namespace meshkernel
{
    // coordinate reference independent operations
    template <typename T>
    [[nodiscard]] T DotProduct(const T& dx1, const T& dx2)
    {
        return dx1 * dx2;
    }

    template <typename T, typename... Args>
    [[nodiscard]] T DotProduct(const T& dx1, const T& dx2, Args&... args)
    {
        return dx1 * dx2 + DotProduct(args...);
    }

    template <typename T>
    void ResizeVectorIfNeeded(int newSize, std::vector<T>& vectorToResize, T fillValue = T())
    {
        const int currentSize = int(vectorToResize.size());
        if (newSize > currentSize)
        {
            newSize = std::max(newSize, int(currentSize * 1.2));
            vectorToResize.resize(newSize, fillValue);
        }
    }

    template <typename T>
    void ResizeVectorIfNeededWithMinimumSize(int newSize, std::vector<T>& vectorToResize, int minSize, T fillValue = T())
    {
        const int currentSize = int(vectorToResize.size());
        if (newSize > currentSize)
        {
            newSize = std::max(minSize, int(5 * newSize));
            vectorToResize.resize(newSize, fillValue);
        }
    }

    template <typename T>
    [[nodiscard]] int FindIndex(const std::vector<T>& vec, T el)
    {
        int index = 0;
        for (int n = 0; n < vec.size(); n++)
        {
            if (vec[n] == el)
            {
                index = n;
                break;
            }
        }
        return index;
    }

    static int FindIndexes(const std::vector<Point>& vec,
                           const size_t start,
                           const size_t end,
                           const double& separator,
                           std::vector<std::vector<size_t>>& indexes)
    {
        if (vec.empty())
        {
            return 0;
        }
        // set an invalid index
        auto invalidIndex = std::numeric_limits<size_t>::max();
        for (int n = 0; n < indexes.size(); n++)
        {
            indexes[n][0] = invalidIndex;
            indexes[n][1] = invalidIndex;
        }

        if (start > vec.size() || end > vec.size() || indexes.empty())
        {
            return -1;
        }

        int pos = 0;
        for (auto n = start; n < end; n++)
        {

            if (vec[n].x != separator && indexes[pos][0] == invalidIndex)
            {
                indexes[pos][0] = n;
            }
            else if (vec[n].x == separator && indexes[pos][1] == invalidIndex)
            {
                indexes[pos][1] = n - 1;
                pos++;
            }
            if (pos >= indexes.size())
            {
                return -1;
            }
        }

        if (pos < indexes.size() && indexes[pos][1] == invalidIndex)
        {
            indexes[pos][1] = int(vec.size()) - 1;
            pos++;
        }

        return pos;
    }

    template <typename T>
    [[nodiscard]] std::vector<int> SortedIndexes(const std::vector<T>& v)
    {
        std::vector<int> idx(v.size());
        iota(idx.begin(), idx.end(), 0);
        std::stable_sort(idx.begin(), idx.end(), [&v](size_t i1, size_t i2) { return v[i1] < v[i2]; });
        return idx;
    }

    //chmike's algorithm
    template <class T>
    void ReorderVector(std::vector<T>& v, std::vector<int> const& order)
    {
        std::vector<T> ordered(v.size());
        for (int i = 0; i < order.size(); ++i)
        {
            ordered[i] = v[order[i]];
        }
        v = ordered;
    }

    template <typename T>
    void MakeMonothonic(std::vector<T>& v)
    {

        bool isMonotonic = false;
        int maxIter = 10;
        int iter = 0;
        while (!isMonotonic && iter < maxIter)
        {
            isMonotonic = true;
            maxIter++;
            for (int n = 0; n < v.size(); ++n)
            {
                if (v[n + 1] - v[n] < 0.0)
                {
                    isMonotonic = false;
                    break;
                }
            }
            if (!isMonotonic)
            {
                for (int n = 1; n < v.size() - 1; ++n)
                {
                    v[n] = 0.5 * (v[n - 1] + v[n + 1]);
                }
            }
        }
    }

    template <typename T>
    void AddValueToVector(std::vector<T>& vec, const T value)
    {
        for (auto& val : vec)
        {
            val += value;
        }
    }

    // algorithm performing the zero's search using the golden section algorithm's
    template <typename F>
    [[nodiscard]] double FindFunctionRootWithGoldenSectionSearch(F func, double min, double max)
    {
        //golden distance factors
        const double c = 0.38196602;
        const double r = 0.61803399;
        const double tolerance = 0.00001;

        double left = min;
        double middle = (min + max) * 0.5;
        double right = max;

        double x0 = left;
        double x1 = middle - c * (middle - left);
        double x2 = middle;
        double x3 = right;
        if (std::abs(right - middle) > std::abs(middle - left))
        {
            x1 = middle;
            x2 = middle + c * (middle - left);
        }

        double f1 = func(x1);
        double f2 = func(x2);

        while (std::abs(x3 - x0) > tolerance * std::max(std::abs(x1) + std::abs(x2), 1e-8))
        {
            if (f2 < f1)
            {
                x0 = x1;
                x1 = x2;
                x2 = r * x1 + c * x3;

                f1 = f2;
                f2 = func(x2);
            }
            else
            {
                x3 = x2;
                x2 = x1;
                x1 = r * x2 + c * x0;

                f2 = f1;
                f1 = func(x1);
            }
        }

        return f1 < f2 ? x1 : x2;
    }

    [[nodiscard]] static bool IsPointOnPole(const Point& point)
    {
        return std::abs(std::abs(point.y) - 90.0) < absLatitudeAtPoles;
    }

    ///sphertocart3D transform 2d spherical to 3d cartesian
    static void SphericalToCartesian(const Point& sphericalPoint, Cartesian3DPoint& cartesianPoint)
    {
        cartesianPoint.z = earth_radius * sin(sphericalPoint.y * degrad_hp);
        double rr = earth_radius * cos(sphericalPoint.y * degrad_hp);
        cartesianPoint.x = rr * cos(sphericalPoint.x * degrad_hp);
        cartesianPoint.y = rr * sin(sphericalPoint.x * degrad_hp);
    }

    ///Cart3Dtospher Transform 3d cartesian coordinates to 2d spherical
    static void CartesianToSpherical(const Cartesian3DPoint& cartesianPoint, double referenceLongitude, Point& sphericalPoint)
    {
        double angle = atan2(cartesianPoint.y, cartesianPoint.x) * raddeg_hp;
        sphericalPoint.y = atan2(cartesianPoint.z, sqrt(cartesianPoint.x * cartesianPoint.x + cartesianPoint.y * cartesianPoint.y)) * raddeg_hp;
        sphericalPoint.x = angle + std::lround((referenceLongitude - angle) / 360.0) * 360.0;
    }

    /// IsLeft(): tests if a point is Left|On|Right of an infinite line.
    ///    Input:  three points leftPoint, rightPoint, and point
    ///    Return: >0 for point left of the line through leftPoint and rightPoint
    ///            =0 for point  on the line
    ///            <0 for point  right of the line
    [[nodiscard]] static double IsLeft(const Point& leftPoint, const Point& rightPoint, const Point& point)
    {
        double left = (rightPoint.x - leftPoint.x) * (point.y - leftPoint.y) - (point.x - leftPoint.x) * (rightPoint.y - leftPoint.y);
        return left;
    }

    /// Check if a point is in polygonNodes using the winding number method
    /// polygonNodes: a closed polygonNodes consisting f a vector of numberOfPolygonPoints + 1 in counter clockwise order
    [[nodiscard]] static bool IsPointInPolygonNodes(const Point& point, const std::vector<Point>& polygonNodes, int startNode, int endNode)
    {
        if (endNode <= startNode)
        {
            return true;
        }
        const int currentPolygonSize = endNode - startNode + 1;
        if (polygonNodes.size() < currentPolygonSize)
        {
            return false;
        }
        if (polygonNodes[startNode] != polygonNodes[endNode])
        {
            return false;
        }

        int windingNumber = 0;
        for (int n = startNode; n < endNode; n++)
        {
            if (polygonNodes[n].y <= point.y) // an upward crossing
            {
                if (polygonNodes[n + 1].y >= point.y && IsLeft(polygonNodes[n], polygonNodes[n + 1], point) >= 0.0)
                {
                    ++windingNumber; // have  a valid up intersect
                }
            }
            else
            {
                if (polygonNodes[n + 1].y <= point.y && IsLeft(polygonNodes[n], polygonNodes[n + 1], point) <= 0.0) // a downward crossing
                {
                    --windingNumber; // have  a valid down intersect
                }
            }
        }
        return windingNumber == 0 ? false : true;
    }

    static void ComputeThreeBaseComponents(const Point& point, double (&exxp)[3], double (&eyyp)[3], double (&ezzp)[3])
    {
        double phi0 = point.y * degrad_hp;
        double lambda0 = point.x * degrad_hp;

        exxp[0] = cos(phi0) * cos(lambda0);
        exxp[1] = cos(phi0) * sin(lambda0);
        exxp[2] = sin(phi0);

        eyyp[0] = -sin(lambda0);
        eyyp[1] = cos(lambda0);
        eyyp[2] = 0.0;

        ezzp[0] = -sin(phi0) * cos(lambda0);
        ezzp[1] = -sin(phi0) * sin(lambda0);
        ezzp[2] = cos(phi0);
    };

    static void ComputeTwoBaseComponents(const Point& point, double (&elambda)[3], double (&ephi)[3])
    {
        double phi0 = point.y * degrad_hp;
        double lambda0 = point.x * degrad_hp;

        elambda[0] = -sin(lambda0);
        elambda[1] = cos(lambda0);
        elambda[2] = 0.0;

        ephi[0] = -sin(phi0) * cos(lambda0);
        ephi[1] = -sin(phi0) * sin(lambda0);
        ephi[2] = cos(phi0);
    };

    [[nodiscard]] static double GetDx(const Point& firstPoint, const Point& secondPoint, const Projections& projection)
    {
        double delta = secondPoint.x - firstPoint.x;
        if (std::abs(delta) <= nearlyZero)
        {
            return 0.0;
        }

        if (projection == Projections::cartesian)
        {
            return delta;
        }
        if (projection == Projections::spherical || projection == Projections::sphericalAccurate)
        {
            bool isFirstPointOnPole = IsPointOnPole(firstPoint);
            bool isSecondPointOnPole = IsPointOnPole(secondPoint);
            if (isFirstPointOnPole && !isSecondPointOnPole || !isFirstPointOnPole && isSecondPointOnPole)
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
        return doubleMissingValue;
    }

    [[nodiscard]] static double GetDy(const Point& firstPoint, const Point& secondPoint, const Projections& projection)
    {
        double delta = secondPoint.y - firstPoint.y;
        if (std::abs(delta) <= nearlyZero)
        {
            return 0.0;
        }

        if (projection == Projections::cartesian)
        {
            return delta;
        }
        if (projection == Projections::spherical || projection == Projections::sphericalAccurate)
        {
            double firstPointY = firstPoint.y * degrad_hp;
            double secondPointY = secondPoint.y * degrad_hp;
            double dy = earth_radius * (secondPointY - firstPointY);
            return dy;
        }
        return doubleMissingValue;
    }

    ///dprodout: out product of two segments
    [[nodiscard]] static double OuterProductTwoSegments(const Point& firstPointFirstSegment, const Point& secondPointFirstSegment,
                                                        const Point& firstPointSecondSegment, const Point& secondPointSecondSegment, const Projections& projection)
    {
        if (projection == Projections::sphericalAccurate)
        {
            Cartesian3DPoint firstPointFirstSegmentCartesian;
            SphericalToCartesian(firstPointFirstSegment, firstPointFirstSegmentCartesian);
            auto xx1 = firstPointFirstSegmentCartesian.x;
            auto yy1 = firstPointFirstSegmentCartesian.y;
            auto zz1 = firstPointFirstSegmentCartesian.z;

            Cartesian3DPoint secondPointFirstSegmentCartesian;
            SphericalToCartesian(secondPointFirstSegment, secondPointFirstSegmentCartesian);
            auto xx2 = secondPointFirstSegmentCartesian.x;
            auto yy2 = secondPointFirstSegmentCartesian.y;
            auto zz2 = secondPointFirstSegmentCartesian.z;

            Cartesian3DPoint firstPointSecondSegmentCartesian;
            SphericalToCartesian(firstPointSecondSegment, firstPointSecondSegmentCartesian);
            auto xx3 = firstPointSecondSegmentCartesian.x;
            auto yy3 = firstPointSecondSegmentCartesian.y;
            auto zz3 = firstPointSecondSegmentCartesian.z;

            Cartesian3DPoint secondPointSecondSegmentCartesian;
            SphericalToCartesian(secondPointSecondSegment, secondPointSecondSegmentCartesian);
            auto xx4 = secondPointSecondSegmentCartesian.x;
            auto yy4 = secondPointSecondSegmentCartesian.y;
            auto zz4 = secondPointSecondSegmentCartesian.z;

            double vxx = (yy2 - yy1) * (zz4 - zz3) - (zz2 - zz1) * (yy4 - yy3);
            double vyy = (zz2 - zz1) * (xx4 - xx3) - (xx2 - xx1) * (zz4 - zz3);
            double vzz = (xx2 - xx1) * (yy4 - yy3) - (yy2 - yy1) * (xx4 - xx3);

            double result = std::sqrt(vxx * vxx + vyy * vyy + vzz * vzz);

            //check if vector is pointing outwards of earth
            if (vxx * xx1 + vyy * yy1 + vzz * zz1 < 0.0)
            {
                result = -result;
            }
            return result;
        }

        // cartesian and spherical
        if (projection == Projections::cartesian || projection == Projections::spherical)
        {
            double dx1 = GetDx(firstPointFirstSegment, secondPointFirstSegment, projection);
            double dx2 = GetDx(firstPointSecondSegment, secondPointSecondSegment, projection);

            double dy1 = GetDy(firstPointFirstSegment, secondPointFirstSegment, projection);
            double dy2 = GetDy(firstPointSecondSegment, secondPointSecondSegment, projection);

            return dx1 * dy2 - dy1 * dx2;
        }

        return doubleMissingValue;
    }

    /// half
    static void MiddlePoint(const Point& firstPoint, const Point& secondPoint, Point& result, const Projections& projection)
    {
        result.x = doubleMissingValue;
        result.y = doubleMissingValue;
        if (!firstPoint.IsValid() || !secondPoint.IsValid())
        {
            return;
        }
        if (projection == Projections::sphericalAccurate)
        {
            Cartesian3DPoint firstPointCartesianCoordinates{doubleMissingValue, doubleMissingValue};
            SphericalToCartesian(firstPoint, firstPointCartesianCoordinates);
            Cartesian3DPoint secondPointCartesianCoordinates{doubleMissingValue, doubleMissingValue};
            SphericalToCartesian(secondPoint, secondPointCartesianCoordinates);

            Cartesian3DPoint middleCartesianPointCoordinate{doubleMissingValue, doubleMissingValue};
            middleCartesianPointCoordinate.x = 0.5 * (firstPointCartesianCoordinates.x + secondPointCartesianCoordinates.x);
            middleCartesianPointCoordinate.y = 0.5 * (firstPointCartesianCoordinates.y + secondPointCartesianCoordinates.y);
            double referenceLongitude = std::max(firstPoint.x, secondPoint.x);

            CartesianToSpherical(middleCartesianPointCoordinate, referenceLongitude, result);
        }

        // cartesian and spherical
        if (projection == Projections::cartesian || projection == Projections::spherical)
        {
            result = (firstPoint + secondPoint) * 0.5;
        }
    }

    static void ComputeMiddlePoint(const Point& firstPoint, const Point& secondPoint, const Projections& projection, Point& centre)
    {
        centre.x = doubleMissingValue;
        centre.y = doubleMissingValue;
        if (projection == Projections::spherical || projection == Projections::sphericalAccurate)
        {
            centre.y = (firstPoint.y + secondPoint.y) / 2.0;
            const auto isFirstNodeOnPole = IsPointOnPole(firstPoint);
            const auto isSecondNodeOnPole = IsPointOnPole(secondPoint);

            if (isFirstNodeOnPole && !isSecondNodeOnPole)
            {
                centre.x = secondPoint.x;
            }
            else if (!isFirstNodeOnPole && isSecondNodeOnPole)
            {
                centre.x = firstPoint.x;
            }
            else
            {
                const auto maxx = std::max(firstPoint.x, secondPoint.x);
                const auto minx = std::min(firstPoint.x, secondPoint.x);

                if (std::abs(maxx - minx) > std::numeric_limits<double>::epsilon())
                {
                    centre.x = centre.x + 180.0;
                }
            }
        }

        // cartesian
        if (projection == Projections::cartesian)
        {
            centre = (firstPoint + secondPoint) * 0.5;
        }
    }

    ///normalin, Normalized vector in direction 1 -> 2, in the orientation of (xu,yu)
    static void NormalVector(const Point& firstPoint, const Point& secondPoint, const Point& insidePoint, Point& result, const Projections& projection)
    {
        if (projection == Projections::sphericalAccurate)
        {
            Cartesian3DPoint firstPointCartesianCoordinates;
            Cartesian3DPoint secondPointCartesianCoordinates;
            SphericalToCartesian(firstPoint, firstPointCartesianCoordinates);
            SphericalToCartesian(secondPoint, secondPointCartesianCoordinates);

            double elambda[3];
            double ephi[3];
            ComputeTwoBaseComponents(insidePoint, elambda, ephi);

            double dx = (secondPointCartesianCoordinates.x - firstPointCartesianCoordinates.x) * elambda[0] +
                        (secondPointCartesianCoordinates.y - firstPointCartesianCoordinates.y) * elambda[1] +
                        (secondPointCartesianCoordinates.z - firstPointCartesianCoordinates.z) * elambda[2];

            double dy = (secondPointCartesianCoordinates.x - firstPointCartesianCoordinates.x) * ephi[0] +
                        (secondPointCartesianCoordinates.y - firstPointCartesianCoordinates.y) * ephi[1] +
                        (secondPointCartesianCoordinates.z - firstPointCartesianCoordinates.z) * ephi[2];

            double squaredDistance = dx * dx + dy * dy;
            if (squaredDistance > 0.0)
            {
                const double distance = sqrt(squaredDistance);
                result.x = dx / distance;
                result.y = dy / distance;
            }
        }

        // cartesian and spherical
        if (projection == Projections::cartesian || projection == Projections::spherical)
        {
            double dx = GetDx(firstPoint, secondPoint, projection);
            double dy = GetDy(firstPoint, secondPoint, projection);
            const double squaredDistance = dx * dx + dy * dy;
            if (squaredDistance > 0.0)
            {
                const double distance = sqrt(squaredDistance);
                result.x = dx / distance;
                result.y = dy / distance;
            }
        }
    }

    //spher2locvec, transforms vector with componentis in global spherical coordinate directions(xglob, yglob)
    ///to local coordinate directions(xloc, yloc) around reference point(xref, yref)
    static void TransformGlobalVectorToLocal(const Point& reference, const Point& globalCoordinates, const Point& globalComponents, Projections projection, Point& localComponents)
    {
        if (projection == Projections::sphericalAccurate)
        {
            double exxp[3];
            double eyyp[3];
            double ezzp[3];
            ComputeThreeBaseComponents(reference, exxp, eyyp, ezzp);

            // get the 3D coordinate
            Cartesian3DPoint globalCoordinatesCartesian;
            SphericalToCartesian(globalCoordinates, globalCoordinatesCartesian);

            //project to rotated frame
            Cartesian3DPoint globalCoordinatesCartesianRotated;
            globalCoordinatesCartesianRotated.x = exxp[0] * globalCoordinatesCartesian.x + exxp[1] * globalCoordinatesCartesian.y + exxp[2] * globalCoordinatesCartesian.z;
            globalCoordinatesCartesianRotated.y = eyyp[0] * globalCoordinatesCartesian.x + eyyp[1] * globalCoordinatesCartesian.y + eyyp[2] * globalCoordinatesCartesian.z;
            globalCoordinatesCartesianRotated.z = ezzp[0] * globalCoordinatesCartesian.x + ezzp[1] * globalCoordinatesCartesian.y + ezzp[2] * globalCoordinatesCartesian.z;

            // Compute global base vectors at other point in 3D(xx, yy, zz) frame
            double elambda[3];
            double ephi[3];
            ComputeTwoBaseComponents(globalCoordinates, elambda, ephi);

            double vxx = globalComponents.x * elambda[0] + globalComponents.y * ephi[0];
            double vyy = globalComponents.x * elambda[1] + globalComponents.y * ephi[1];
            double vzz = globalComponents.x * elambda[2] + globalComponents.y * ephi[2];

            //tranform to local spherical coordinates
            Point globalCoordinatesToLocal;
            CartesianToSpherical(globalCoordinatesCartesianRotated, reference.x, globalCoordinatesToLocal);

            //compute base vectors at other point in rotated 3D(xxp, yyp, zzp) frame
            double elambdap[3];
            double ephip[3];
            ComputeTwoBaseComponents(globalCoordinatesToLocal, elambdap, ephip);

            //compute local base vectors in(xx, yy, zz) frame
            double elambdaloc[3];
            elambdaloc[0] = exxp[0] * elambdap[0] + eyyp[0] * elambdap[1] + ezzp[0] * elambda[2];
            elambdaloc[1] = exxp[1] * elambdap[0] + eyyp[1] * elambdap[1] + ezzp[1] * elambda[2];
            elambdaloc[2] = exxp[2] * elambdap[0] + eyyp[2] * elambdap[1] + ezzp[2] * elambda[2];

            double ephiloc[3];
            ephiloc[0] = exxp[0] * ephip[0] + eyyp[0] * ephip[1] + ezzp[0] * ephip[2];
            ephiloc[1] = exxp[1] * ephip[0] + eyyp[1] * ephip[1] + ezzp[1] * ephip[2];
            ephiloc[2] = exxp[2] * ephip[0] + eyyp[2] * ephip[1] + ezzp[2] * ephip[2];

            //compute vectors in other point in local base(elambdaloc, ephiloc)
            localComponents.x = elambdaloc[0] * vxx + elambdaloc[1] * vyy + elambdaloc[2] * vzz;
            localComponents.y = ephiloc[0] * vxx + ephiloc[1] * vyy + ephiloc[2] * vzz;
        }
        else
        {
            // cartesian and spherical
            if (projection == Projections::cartesian || projection == Projections::spherical)
            {
                localComponents = globalComponents;
            }
        }
    }

    ///normalout
    static void NormalVectorOutside(const Point& firstPoint, const Point& secondPoint, Point& result, const Projections& projection)
    {
        result.x = doubleMissingValue;
        result.y = doubleMissingValue;
        if (projection == Projections::sphericalAccurate)
        {
            Point middlePoint;
            MiddlePoint(firstPoint, secondPoint, middlePoint, projection);

            Cartesian3DPoint firstPointCartesianCoordinates;
            SphericalToCartesian(firstPoint, firstPointCartesianCoordinates);
            Cartesian3DPoint secondPointCartesianCoordinates;
            SphericalToCartesian(secondPoint, secondPointCartesianCoordinates);

            //compute the base vectors at middle point
            double elambda[3];
            double ephi[3];
            ComputeTwoBaseComponents(middlePoint, elambda, ephi);

            // project vector in local base
            double dx = (secondPointCartesianCoordinates.x - firstPointCartesianCoordinates.x) * elambda[0] +
                        (secondPointCartesianCoordinates.y - firstPointCartesianCoordinates.y) * elambda[1] +
                        (secondPointCartesianCoordinates.z - firstPointCartesianCoordinates.z) * elambda[2];

            double dy = (secondPointCartesianCoordinates.x - firstPointCartesianCoordinates.x) * ephi[0] +
                        (secondPointCartesianCoordinates.y - firstPointCartesianCoordinates.y) * ephi[1] +
                        (secondPointCartesianCoordinates.z - firstPointCartesianCoordinates.z) * ephi[2];

            const double squaredDistance = dx * dx + dy * dy;
            if (squaredDistance > 0.0)
            {
                const double distance = sqrt(squaredDistance);
                result.x = dy / distance;
                result.y = -dx / distance;
            }
        }

        // cartesian and spherical
        if (projection == Projections::cartesian || projection == Projections::spherical)
        {
            double dx = GetDx(firstPoint, secondPoint, projection);
            double dy = GetDy(firstPoint, secondPoint, projection);

            const double squaredDistance = dx * dx + dy * dy;
            if (squaredDistance > 0.0)
            {
                const double distance = sqrt(squaredDistance);
                result.x = dy / distance;
                result.y = -dx / distance;
            }

            if (projection == Projections::spherical)
            {
                Point middlePoint;
                MiddlePoint(firstPoint, secondPoint, middlePoint, projection);
                result.x = result.x / cos(degrad_hp * 0.5 * (firstPoint.y + secondPoint.y));
                result.y = result.y;
            }
        }
    }

    ///normaloutchk
    ///Computes the normal vector to a line 1-2, which is *outward* w.r.t.
    ///an 'inside' point 3. Similar to normalout, except that the normal
    ///vector may be flipped based on the 'inside' point.
    ///TODO:test me
    static void NormalVectorInside(const Point& firstPoint, const Point& secondPoint, const Point& insidePoint, Point& normal, bool& flippedNormal, Projections projection)
    {
        NormalVectorOutside(firstPoint, secondPoint, normal, projection);
        Point thirdPoint;
        if (projection == Projections::cartesian || projection == Projections::spherical)
        {
            flippedNormal = false;
            thirdPoint.x = firstPoint.x + normal.x;
            thirdPoint.y = firstPoint.y + normal.y;
        }
        if (projection == Projections::sphericalAccurate)
        {
            Point middle;
            MiddlePoint(firstPoint, secondPoint, middle, projection);
            Point localComponents;
            TransformGlobalVectorToLocal(firstPoint, middle, normal, projection, localComponents);
            double elambda[3];
            double ephi[3];
            ComputeTwoBaseComponents(firstPoint, elambda, ephi);

            double vxx = localComponents.x * elambda[0] + localComponents.y * ephi[0];
            double vyy = localComponents.x * elambda[1] + localComponents.y * ephi[1];
            double vzz = localComponents.x * elambda[2] + localComponents.y * ephi[2];

            Cartesian3DPoint firstPointCartesian;
            SphericalToCartesian(firstPoint, firstPointCartesian);

            Cartesian3DPoint rotatedPoint;
            double alpha = 0.0;
            rotatedPoint.x = firstPointCartesian.x + alpha * vxx;
            rotatedPoint.y = firstPointCartesian.y + alpha * vyy;
            rotatedPoint.z = firstPointCartesian.z + alpha * vzz;

            CartesianToSpherical(rotatedPoint, firstPoint.x, thirdPoint);
        }

        if (OuterProductTwoSegments(firstPoint, thirdPoint, firstPoint, secondPoint, projection) * OuterProductTwoSegments(firstPoint, insidePoint, firstPoint, secondPoint, projection) > 0.0)
        {
            normal.x = -normal.x;
            normal.y = -normal.y;
            flippedNormal = true;
        }
        else
        {
            flippedNormal = false;
        }
    }

    static void Add(Point& point, const Point& normal, double increment, double xf, const Projections& projection)
    {
        if (projection == Projections::cartesian)
        {
            point.x = point.x + normal.x * increment;
            point.y = point.y + normal.y * increment;
        }
        if (projection == Projections::spherical || projection == Projections::sphericalAccurate)
        {
            double convertedIncrement = raddeg_hp * increment / earth_radius;
            point.x = point.x + normal.x * convertedIncrement * xf;
            point.y = point.y + normal.y * convertedIncrement;
        }
    }

    static void ReferencePoint(std::vector<Point>& polygon, const int numPoints, double& minX, double& minY, const Projections& projection)
    {
        minX = std::numeric_limits<double>::max();
        minY = std::numeric_limits<double>::max();
        for (int i = 0; i < numPoints; i++)
        {
            minX = std::min(polygon[i].x, minX);
            if (abs(polygon[i].y) < abs(minY))
            {
                minY = polygon[i].y;
            }
        }

        if (projection == Projections::spherical)
        {
            double maxX = std::numeric_limits<double>::min();
            for (int i = 0; i < numPoints; i++)
            {
                maxX = std::max(polygon[i].x, maxX);
            }

            if (maxX - minX > 180.0)
            {
                double deltaX = maxX - 180.0;
                for (int i = 0; i < numPoints; i++)
                {
                    if (polygon[i].x < deltaX)
                    {
                        polygon[i].x = polygon[i].x + 360.0;
                    }
                }
                minX = minX + 360.0;
            }
        }
    }

    [[nodiscard]] static double ComputeSquaredDistance(const Point& firstPoint, const Point& secondPoint, const Projections& projection)
    {

        if (!firstPoint.IsValid() || !secondPoint.IsValid())
        {
            return 0.0;
        }

        if (projection == Projections::sphericalAccurate)
        {
            Cartesian3DPoint firstPointCartesian{doubleMissingValue, doubleMissingValue};
            SphericalToCartesian(firstPoint, firstPointCartesian);
            auto xx1 = firstPointCartesian.x;
            auto yy1 = firstPointCartesian.y;
            auto zz1 = firstPointCartesian.z;

            Cartesian3DPoint secondPointCartesian{doubleMissingValue, doubleMissingValue};
            SphericalToCartesian(secondPoint, secondPointCartesian);
            auto xx2 = secondPointCartesian.x;
            auto yy2 = secondPointCartesian.y;
            auto zz2 = secondPointCartesian.z;

            return (xx2 - xx1) * (xx2 - xx1) + (yy2 - yy1) * (yy2 - yy1) + (zz2 - zz1) * (zz2 - zz1);
        }

        //cartesian and spherical
        if (projection == Projections::cartesian || projection == Projections::spherical)
        {
            double dx = GetDx(firstPoint, secondPoint, projection);
            double dy = GetDy(firstPoint, secondPoint, projection);
            return dx * dx + dy * dy;
        }

        return doubleMissingValue;
    }

    //dbdistance
    [[nodiscard]] static double Distance(const Point& firstPoint, const Point& secondPoint, const Projections& projection)
    {
        double distance = ComputeSquaredDistance(firstPoint, secondPoint, projection);
        if (distance >= 0.0)
        {
            distance = sqrt(distance);
        }
        return distance;
    }

    // dLINEDIS3
    // Computes the perpendicular distance from point to a line firstNode - secondNode.
    // normalPoint: coordinates of the projected point from point onto the line
    [[nodiscard]] static double DistanceFromLine(const Point& point, const Point& firstNode, const Point& secondNode, Point& normalPoint, double& ratio, const Projections& projection)
    {
        if (projection == Projections::cartesian || projection == Projections::spherical)
        {
            double dis = 0.0;
            double squaredDistance = ComputeSquaredDistance(secondNode, firstNode, projection);
            if (squaredDistance != 0.0)
            {
                ratio = (GetDx(firstNode, point, projection) * GetDx(firstNode, secondNode, projection) +
                         GetDy(firstNode, point, projection) * GetDy(firstNode, secondNode, projection)) /
                        squaredDistance;
                double correctedRatio = std::max(std::min(1.0, ratio), 0.0);
                normalPoint.x = firstNode.x + correctedRatio * (secondNode.x - firstNode.x);
                normalPoint.y = firstNode.y + correctedRatio * (secondNode.y - firstNode.y);
                dis = Distance(point, normalPoint, projection);
            }
            return dis;
        }

        if (projection == Projections::sphericalAccurate)
        {
            Cartesian3DPoint firstNodeCartesian;
            SphericalToCartesian(firstNode, firstNodeCartesian);
            auto xx1 = firstNodeCartesian.x;
            auto yy1 = firstNodeCartesian.y;
            auto zz1 = firstNodeCartesian.z;

            Cartesian3DPoint secondNodeCartesian;
            SphericalToCartesian(secondNode, secondNodeCartesian);
            auto xx2 = secondNodeCartesian.x;
            auto yy2 = secondNodeCartesian.y;
            auto zz2 = secondNodeCartesian.z;

            Cartesian3DPoint pointCartesian;
            SphericalToCartesian(point, pointCartesian);
            auto xx3 = pointCartesian.x;
            auto yy3 = pointCartesian.y;
            auto zz3 = pointCartesian.z;

            double x21 = xx2 - xx1;
            double y21 = yy2 - yy1;
            double z21 = zz2 - zz1;
            double x31 = xx3 - xx1;
            double y31 = yy3 - yy1;
            double z31 = zz3 - zz1;

            double r2 = x21 * x21 + y21 * y21 + z21 * z21;

            ratio = 0.0;
            if (r2 >= 0.0)
            {

                ratio = (x31 * x21 + y31 * y21 + z31 * z21) / r2;
                double correctedRatio = std::max(std::min(1.0, ratio), 0.0);

                Cartesian3DPoint cartesianNormal3DPoint;
                cartesianNormal3DPoint.x = firstNodeCartesian.x + correctedRatio * x21;
                cartesianNormal3DPoint.y = firstNodeCartesian.y + correctedRatio * y21;
                cartesianNormal3DPoint.z = firstNodeCartesian.z + correctedRatio * z21;

                cartesianNormal3DPoint.x = cartesianNormal3DPoint.x - xx3;
                cartesianNormal3DPoint.y = cartesianNormal3DPoint.y - yy3;
                cartesianNormal3DPoint.z = cartesianNormal3DPoint.z - zz3;

                double dis = std::sqrt(cartesianNormal3DPoint.x * cartesianNormal3DPoint.x +
                                       cartesianNormal3DPoint.y * cartesianNormal3DPoint.y +
                                       cartesianNormal3DPoint.z * cartesianNormal3DPoint.z);

                double referenceLongitude = std::max({firstNode.x, secondNode.x, point.x});
                CartesianToSpherical(cartesianNormal3DPoint, referenceLongitude, normalPoint);

                return dis;
            }
        }

        return -1.0;
    }

    /// dprodin inner product of two segments
    [[nodiscard]] static double InnerProductTwoSegments(const Point& firstPointFirstSegment, const Point& secondPointFirstSegment, const Point& firstPointSecondSegment, const Point& secondPointSecondSegment, const Projections& projection)
    {
        if (projection == Projections::sphericalAccurate)
        {
            Cartesian3DPoint firstPointFirstSegment3D;
            Cartesian3DPoint secondPointFirstSegment3D;
            Cartesian3DPoint firstPointSecondSegment3D;
            Cartesian3DPoint secondPointSecondSegment3D;

            SphericalToCartesian(firstPointFirstSegment, firstPointFirstSegment3D);
            SphericalToCartesian(secondPointFirstSegment, secondPointFirstSegment3D);
            SphericalToCartesian(firstPointSecondSegment, firstPointSecondSegment3D);
            SphericalToCartesian(secondPointSecondSegment, secondPointSecondSegment3D);

            double dx1 = secondPointFirstSegment3D.x - firstPointFirstSegment3D.x;
            double dy1 = secondPointFirstSegment3D.y - firstPointFirstSegment3D.y;
            double dz1 = secondPointFirstSegment3D.z - firstPointFirstSegment3D.z;

            double dx2 = secondPointSecondSegment3D.x - firstPointSecondSegment3D.x;
            double dy2 = secondPointSecondSegment3D.y - firstPointSecondSegment3D.y;
            double dz2 = secondPointSecondSegment3D.z - firstPointSecondSegment3D.z;

            return DotProduct(dx1, dx2, dy1, dy2, dz1, dz2);
        }

        // cartesian and spherical
        if (projection == Projections::cartesian || projection == Projections::spherical)
        {
            double dx1 = GetDx(firstPointFirstSegment, secondPointFirstSegment, projection);
            double dx2 = GetDx(firstPointSecondSegment, secondPointSecondSegment, projection);

            double dy1 = GetDy(firstPointFirstSegment, secondPointFirstSegment, projection);
            double dy2 = GetDy(firstPointSecondSegment, secondPointSecondSegment, projection);

            return dx1 * dx2 + dy1 * dy2;
        }

        return doubleMissingValue;
    }

    // dcosphi
    [[nodiscard]] static double NormalizedInnerProductTwoSegments(const Point& firstPointFirstSegment, const Point& secondPointFirstSegment, const Point& firstPointSecondSegment, const Point& secondPointSecondSegment, const Projections& projection)
    {
        if (projection == Projections::sphericalAccurate)
        {
            Cartesian3DPoint firstPointFirstSegmentCartesian;
            SphericalToCartesian(firstPointFirstSegment, firstPointFirstSegmentCartesian);
            auto xx1 = firstPointFirstSegmentCartesian.x;
            auto yy1 = firstPointFirstSegmentCartesian.y;
            auto zz1 = firstPointFirstSegmentCartesian.z;

            Cartesian3DPoint secondPointFirstSegmentCartesian;
            SphericalToCartesian(secondPointFirstSegment, secondPointFirstSegmentCartesian);
            auto xx2 = secondPointFirstSegmentCartesian.x;
            auto yy2 = secondPointFirstSegmentCartesian.y;
            auto zz2 = secondPointFirstSegmentCartesian.z;

            Cartesian3DPoint firstPointSecondSegmentCartesian;
            SphericalToCartesian(firstPointSecondSegment, firstPointSecondSegmentCartesian);
            auto xx3 = firstPointSecondSegmentCartesian.x;
            auto yy3 = firstPointSecondSegmentCartesian.y;
            auto zz3 = firstPointSecondSegmentCartesian.z;

            Cartesian3DPoint secondPointSecondSegmentCartesian;
            SphericalToCartesian(secondPointSecondSegment, secondPointSecondSegmentCartesian);
            auto xx4 = secondPointSecondSegmentCartesian.x;
            auto yy4 = secondPointSecondSegmentCartesian.y;
            auto zz4 = secondPointSecondSegmentCartesian.z;

            auto dx1 = xx2 - xx1;
            auto dy1 = yy2 - yy1;
            auto dz1 = zz2 - zz1;
            auto firstSegmentDistance = dx1 * dx1 + dy1 * dy1 + dz1 * dz1;

            auto dx2 = xx4 - xx3;
            auto dy2 = yy4 - yy3;
            auto dz2 = zz4 - zz3;
            auto secondSegmentDistance = dx2 * dx2 + dy2 * dy2 + dz2 * dz2;

            double cosphi;
            if (firstSegmentDistance <= 0.0 || secondSegmentDistance <= 0.0)
            {
                cosphi = doubleMissingValue;
            }
            else
            {
                cosphi = (dx1 * dx2 + dy1 * dy2 + dz1 * dz2) / sqrt(firstSegmentDistance * secondSegmentDistance);
            }
            return cosphi;
        }

        // cartesian and spherical
        if (projection == Projections::cartesian || projection == Projections::spherical)
        {
            const auto dx1 = GetDx(firstPointFirstSegment, secondPointFirstSegment, projection);
            const auto dx2 = GetDx(firstPointSecondSegment, secondPointSecondSegment, projection);

            const auto dy1 = GetDy(firstPointFirstSegment, secondPointFirstSegment, projection);
            const auto dy2 = GetDy(firstPointSecondSegment, secondPointSecondSegment, projection);

            const auto r1 = dx1 * dx1 + dy1 * dy1;
            const auto r2 = dx2 * dx2 + dy2 * dy2;

            double cosphi;
            if (r1 <= 0.0 || r2 <= 0.0)
            {
                return doubleMissingValue;
            }

            cosphi = (dx1 * dx2 + dy1 * dy2) / std::sqrt(r1 * r2);
            cosphi = std::max(std::min(cosphi, 1.0), -1.0);
            return cosphi;
        }

        return doubleMissingValue;
    }

    static void CircumcenterOfTriangle(const Point& p1, const Point& p2, const Point& p3, const Projections projection, Point& circumcenter)
    {

        double dx2 = GetDx(p1, p2, projection);
        double dy2 = GetDy(p1, p2, projection);

        double dx3 = GetDx(p1, p3, projection);
        double dy3 = GetDy(p1, p3, projection);

        double den = dy2 * dx3 - dy3 * dx2;
        double z = 0.0;
        if (std::abs(den) > 0.0)
        {
            z = (dx2 * (dx2 - dx3) + dy2 * (dy2 - dy3)) / den;
        }

        if (projection == Projections::cartesian)
        {
            circumcenter.x = p1.x + 0.5 * (dx3 - z * dy3);
            circumcenter.y = p1.y + 0.5 * (dy3 + z * dx3);
        }
        if (projection == Projections::spherical)
        {
            double phi = (p1.y + p2.y + p3.y) / 3.0;
            double xf = 1.0 / cos(degrad_hp * phi);
            circumcenter.x = p1.x + xf * 0.5 * (dx3 - z * dy3) * raddeg_hp / earth_radius;
            circumcenter.y = p1.y + 0.5 * (dy3 + z * dx3) * raddeg_hp / earth_radius;
        }
        if (projection == Projections::sphericalAccurate)
        {
            //TODO: compute in case of spherical accurate
        }
    }

    /// (CROSS)
    [[nodiscard]] static bool AreLinesCrossing(const Point& firstSegmentFistPoint,
                                               const Point& firstSegmentSecondPoint,
                                               const Point& secondSegmentFistPoint,
                                               const Point& secondSegmentSecondPoint,
                                               bool adimensional,
                                               Point& intersection,
                                               double& crossProduct,
                                               double& firstRatio,
                                               double& secondRatio,
                                               const Projections& projection)
    {
        bool isCrossing = false;

        firstRatio = doubleMissingValue;
        secondRatio = doubleMissingValue;
        auto const x21 = GetDx(firstSegmentFistPoint, firstSegmentSecondPoint, projection);
        auto const y21 = GetDy(firstSegmentFistPoint, firstSegmentSecondPoint, projection);

        auto const x43 = GetDx(secondSegmentFistPoint, secondSegmentSecondPoint, projection);
        auto const y43 = GetDy(secondSegmentFistPoint, secondSegmentSecondPoint, projection);

        auto const x31 = GetDx(firstSegmentFistPoint, secondSegmentFistPoint, projection);
        auto const y31 = GetDy(firstSegmentFistPoint, secondSegmentFistPoint, projection);

        auto const det = x43 * y21 - y43 * x21;

        std::vector<double> values{x21, y21, x43, y43};
        double eps = std::max(0.00001 * (*std::max_element(values.begin(), values.end())), std::numeric_limits<double>::denorm_min());

        if (std::abs(det) < eps)
        {
            return isCrossing;
        }

        secondRatio = (y31 * x21 - x31 * y21) / det;
        firstRatio = (y31 * x43 - x31 * y43) / det;
        if (firstRatio >= 0.0 && firstRatio <= 1.0 && secondRatio >= 0.0 && secondRatio <= 1.0)
        {
            isCrossing = true;
        }
        intersection.x = firstSegmentFistPoint.x + firstRatio * (firstSegmentSecondPoint.x - firstSegmentFistPoint.x);
        intersection.y = firstSegmentFistPoint.y + firstRatio * (firstSegmentSecondPoint.y - firstSegmentFistPoint.y);
        crossProduct = -det;
        if (adimensional)
        {
            crossProduct = -det / (std::sqrt(x21 * x21 + y21 * y21) * std::sqrt(x43 * x43 + y43 * y43) + 1e-8);
        }

        return isCrossing;
    }

    //faceAreaAndCenterOfMass: for cartesian, spherical point and spherical3dPoint
    static void FaceAreaAndCenterOfMass(std::vector<Point>& polygon, int numberOfPolygonPoints, Projections projection, double& area, Point& centerOfMass)
    {
        if (numberOfPolygonPoints <= 0)
        {
            throw std::invalid_argument("FaceAreaAndCenterOfMass: The polygon contains no nodes.");
        }

        double minX;
        double minY;
        ReferencePoint(polygon, numberOfPolygonPoints, minX, minY, projection);

        Point reference{minX, minY};
        area = 0.0;
        double xCenterOfMass = 0.0;
        double yCenterOfMass = 0.0;
        const double minArea = 1e-8;
        for (int p = 0; p < numberOfPolygonPoints; p++)
        {
            double dx0 = GetDx(reference, polygon[p], projection);
            double dy0 = GetDy(reference, polygon[p], projection);
            double dx1 = GetDx(reference, polygon[p + 1], projection);
            double dy1 = GetDy(reference, polygon[p + 1], projection);

            double xc = 0.5 * (dx0 + dx1);
            double yc = 0.5 * (dy0 + dy1);

            dx0 = GetDx(polygon[p], polygon[p + 1], projection);
            dy0 = GetDy(polygon[p], polygon[p + 1], projection);
            double dsx = dy0;
            double dsy = -dx0;
            double xds = xc * dsx + yc * dsy;
            area = area + 0.5 * xds;

            xCenterOfMass = xCenterOfMass + xds * xc;
            yCenterOfMass = yCenterOfMass + xds * yc;
        }
        area = std::abs(area) < minArea ? minArea : area;

        const double fac = 1.0 / (3.0 * area);
        xCenterOfMass = fac * xCenterOfMass;
        yCenterOfMass = fac * yCenterOfMass;

        if (projection == Projections::spherical)
        {
            yCenterOfMass = yCenterOfMass / (earth_radius * degrad_hp);
            xCenterOfMass = xCenterOfMass / (earth_radius * degrad_hp * std::cos((yCenterOfMass + minY) * degrad_hp));
        }

        centerOfMass.x = xCenterOfMass + minX;
        centerOfMass.y = yCenterOfMass + minY;

        area = std::abs(area);
    }

    [[nodiscard]] static double Averaging(const std::vector<Sample>& samples,
                                          int numPolygonNodes,
                                          const std::vector<Point>& polygon,
                                          const Point centerOfMass,
                                          const Projections& projection,
                                          SpatialTrees::RTree& rtree,
                                          int averagingMethod)
    {

        std::vector<Point> searchPolygon(numPolygonNodes);

        // averaging settings
        const double relativeFaceSearchSize = 1.01;
        double minx = std::numeric_limits<double>::max();
        double maxx = std::numeric_limits<double>::min();
        double miny = std::numeric_limits<double>::max();
        double maxy = std::numeric_limits<double>::min();

        for (int i = 0; i < numPolygonNodes; i++)
        {
            searchPolygon[i] = polygon[i] * relativeFaceSearchSize + centerOfMass * (1 - relativeFaceSearchSize);
            minx = std::min(minx, searchPolygon[i].x);
            maxx = std::max(maxx, searchPolygon[i].x);
            miny = std::min(miny, searchPolygon[i].y);
            maxy = std::max(maxy, searchPolygon[i].y);
        }

        if (projection == Projections::spherical && maxx - minx > 180.0)
        {

            double xmean = 0.5 * (maxx + minx);
            minx = std::numeric_limits<double>::max();
            maxx = std::numeric_limits<double>::min();
            for (int i = 0; i < numPolygonNodes; i++)
            {
                if (searchPolygon[i].x < xmean)
                {
                    searchPolygon[i].x = searchPolygon[i].x + 360.0;
                    minx = std::min(minx, searchPolygon[i].x);
                    maxx = std::max(maxx, searchPolygon[i].x);
                }
            }
        }

        double result = doubleMissingValue;
        double searchRadiusSquared = std::numeric_limits<double>::min();
        for (int i = 0; i < numPolygonNodes; i++)
        {
            double squaredDistance = ComputeSquaredDistance(centerOfMass, searchPolygon[i], projection);
            searchRadiusSquared = std::max(searchRadiusSquared, squaredDistance);
        }
        if (searchRadiusSquared <= 0.0)
        {
            throw std::invalid_argument("Averaging: The search radius is <= 0.");
        }

        rtree.NearestNeighboursOnSquaredDistance(centerOfMass, searchRadiusSquared);
        if (rtree.GetQueryResultSize() == 0)
        {
            return result;
        }

        int numValidSamplesInPolygon = 0;
        double wall = 0;
        bool firstValidSampleFound = false;

        for (int i = 0; i < rtree.GetQueryResultSize(); i++)
        {
            //do stuff based on the averaging method
            auto sampleIndex = rtree.GetQuerySampleIndex(i);
            auto sampleValue = samples[sampleIndex].value;
            if (sampleValue <= doubleMissingValue)
            {
                continue;
            }

            Point samplePoint{samples[sampleIndex].x, samples[sampleIndex].y};
            // assume here polygon has a size equal to numPolygonNodes + 1
            bool isInPolygon = IsPointInPolygonNodes(samplePoint, polygon, 0, numPolygonNodes);
            if (isInPolygon)
            {
                if (averagingMethod == SimpleAveraging)
                {
                    if (!firstValidSampleFound)
                    {
                        firstValidSampleFound = true;
                        result = 0.0;
                    }
                    result += sampleValue;
                    numValidSamplesInPolygon++;
                }
                if (averagingMethod == KdTree)
                {
                    if (!firstValidSampleFound)
                    {
                        firstValidSampleFound = true;
                        result = sampleValue;
                    }
                    result = std::min(std::abs(result), std::abs(sampleValue));
                }
                if (averagingMethod == Max)
                {
                    if (!firstValidSampleFound)
                    {
                        firstValidSampleFound = true;
                        result = -std::numeric_limits<double>::max();
                    }
                    result = std::max(result, sampleValue);
                }
                if (averagingMethod == InverseWeightDistance)
                {
                    double distance = std::max(0.01, Distance(centerOfMass, samplePoint, projection));
                    double weight = 1.0 / distance;
                    wall += weight;
                    numValidSamplesInPolygon++;
                    result += weight * sampleValue;
                }
            }
        }

        if (averagingMethod == SimpleAveraging && numValidSamplesInPolygon > 0)
        {
            if (result > doubleMissingValue)
            {
                result /= numValidSamplesInPolygon;
            }
        }

        if (averagingMethod == InverseWeightDistance && numValidSamplesInPolygon > 0)
        {
            result /= wall;
        }

        return result;
    }

    [[nodiscard]] static int NextCircularForwardIndex(int currentIndex, int size)
    {
        int index = currentIndex + 1;
        if (index >= size)
        {
            index = index - size;
        }
        return index;
    }

    [[nodiscard]] static int NextCircularBackwardIndex(int currentIndex, int size)
    {
        int index = currentIndex - 1;
        if (index < 0)
        {
            index = index + size;
        }
        return index;
    }

    template <typename T>
    [[nodiscard]] bool InterpolateSplinePoint(const std::vector<T>& coordinates,
                                              const std::vector<T>& coordinatesDerivatives,
                                              double pointAdimensionalCoordinate,
                                              T& pointCoordinate)
    {
        if (pointAdimensionalCoordinate < 0)
        {
            return false;
        }

        const double eps = 1e-5;
        const double splFac = 1.0;
        const auto intCoordinate = int(std::floor(pointAdimensionalCoordinate));
        if (pointAdimensionalCoordinate - intCoordinate < eps)
        {
            pointCoordinate = coordinates[intCoordinate];
            return true;
        }

        int low = intCoordinate;
        int high = low + 1;
        double a = high - pointAdimensionalCoordinate;
        double b = pointAdimensionalCoordinate - low;

        pointCoordinate = coordinates[low] * a + coordinates[high] * b +
                          (coordinatesDerivatives[low] * (pow(a, 3) - a) + coordinatesDerivatives[high] * (pow(b, 3) - b)) / 6.0 * splFac;

        return true;
    }

    template <class T>
    static void SwapVectorElements(std::vector<T>& v, int numElements)
    {
        if (numElements > v.size())
        {
            return;
        }

        for (int i = 0; i < numElements / 2; i++)
        {
            const auto a = v[i];
            v[i] = v[i + 1];
            v[i + 1] = a;
        }
    }

    static void ComputeAdimensionalDistancesFromPointSerie(const std::vector<Point>& v, Projections projection, std::vector<double>& result, double& totalDistance)
    {
        result[0] = 0;
        for (int i = 1; i < v.size(); i++)
        {
            result[i] = result[i - 1] + Distance(v[i - 1], v[i], projection);
        }
        totalDistance = result.back();
        const double inverseTotalDistance = 1.0 / totalDistance;
        // normalize
        for (int i = 1; i < v.size(); i++)
        {
            result[i] = result[i] * inverseTotalDistance;
        }
    }

    // get the sign
    template <typename T>
    [[nodiscard]] static int sgn(T val)
    {
        return (T(0) < val ? 1 : 0) - (val < T(0) ? 1 : 0);
    }

    //(DUITPL)
    [[nodiscard]] static int TwoSegmentsSign(const Point& p1, const Point& p2, const Point& p3, const Point& p4, Projections projection)
    {

        auto dx1 = GetDx(p1, p2, projection);
        auto dy1 = GetDy(p1, p2, projection);
        auto dx2 = GetDx(p3, p4, projection);
        auto dy2 = GetDy(p3, p4, projection);
        auto val = dx1 * dy2 - dy1 * dx2;
        return sgn(val);
    }

    //(TRANFN2)
    static void InterpolateTransfinite(const std::vector<Point>& sideOne,
                                       const std::vector<Point>& sideTwo,
                                       const std::vector<Point>& sideThree,
                                       const std::vector<Point>& sideFour,
                                       Projections projections,
                                       int numM,
                                       int numN,
                                       std::vector<std::vector<Point>>& result)
    {
        double totalLengthOne;
        std::vector<double> sideOneAdimensional(sideOne.size());
        ComputeAdimensionalDistancesFromPointSerie(sideOne, projections, sideOneAdimensional, totalLengthOne);

        double totalLengthTwo;
        std::vector<double> sideTwoAdimensional(sideTwo.size());
        ComputeAdimensionalDistancesFromPointSerie(sideTwo, projections, sideTwoAdimensional, totalLengthTwo);

        double totalLengthThree;
        std::vector<double> sideThreeAdimensional(sideThree.size());
        ComputeAdimensionalDistancesFromPointSerie(sideThree, projections, sideThreeAdimensional, totalLengthThree);

        double totalLengthFour;
        std::vector<double> sideFourAdimensional(sideFour.size());
        ComputeAdimensionalDistancesFromPointSerie(sideFour, projections, sideFourAdimensional, totalLengthFour);

        // now compute the adimensional distance of each point to be filled
        const int numMPoints = numM + 1;
        const int numNPoints = numN + 1;

        std::vector<std::vector<double>> iWeightFactor(numMPoints, std::vector<double>(numNPoints));
        std::vector<std::vector<double>> jWeightFactor(numMPoints, std::vector<double>(numNPoints));
        for (int i = 0; i < numMPoints; i++)
        {
            for (int j = 0; j < numNPoints; j++)
            {
                const double mWeight = double(i) / double(numM);
                const double nWeight = double(j) / double(numN);

                iWeightFactor[i][j] = (1.0 - nWeight) * sideThreeAdimensional[i] + nWeight * sideFourAdimensional[i];
                jWeightFactor[i][j] = (1.0 - mWeight) * sideOneAdimensional[j] + mWeight * sideTwoAdimensional[j];
            }
        }

        std::vector<std::vector<double>> weightOne(numMPoints, std::vector<double>(numNPoints));
        std::vector<std::vector<double>> weightTwo(numMPoints, std::vector<double>(numNPoints));
        std::vector<std::vector<double>> weightThree(numMPoints, std::vector<double>(numNPoints));
        std::vector<std::vector<double>> weightFour(numMPoints, std::vector<double>(numNPoints));
        for (int i = 0; i < numMPoints; i++)
        {
            for (int j = 0; j < numNPoints; j++)
            {

                weightOne[i][j] = (1.0 - jWeightFactor[i][j]) * totalLengthThree + jWeightFactor[i][j] * totalLengthFour;
                weightTwo[i][j] = (1.0 - iWeightFactor[i][j]) * totalLengthOne + iWeightFactor[i][j] * totalLengthTwo;
                weightThree[i][j] = weightTwo[i][j] / weightOne[i][j];
                weightFour[i][j] = weightOne[i][j] / weightTwo[i][j];
                const double wa = 1.0 / (weightThree[i][j] + weightFour[i][j]);
                weightOne[i][j] = wa * weightThree[i][j];
                weightTwo[i][j] = wa * weightFour[i][j];
            }
        }

        //border points
        result.resize(numMPoints, std::vector<Point>(numNPoints));
        for (int i = 0; i < numMPoints; i++)
        {
            result[i][0] = sideThree[i];
            result[i][numN] = sideFour[i];
        }
        for (int i = 0; i < numNPoints; i++)
        {
            result[0][i] = sideOne[i];
            result[numM][i] = sideTwo[i];
        }

        // first interpolation
        for (int i = 1; i < numM; i++)
        {
            for (int j = 1; j < numN; j++)
            {

                result[i][j].x = (sideOne[j].x * (1.0 - iWeightFactor[i][j]) + sideTwo[j].x * iWeightFactor[i][j]) * weightOne[i][j] +
                                 (sideThree[i].x * (1.0 - jWeightFactor[i][j]) + sideFour[i].x * jWeightFactor[i][j]) * weightTwo[i][j];

                result[i][j].y = (sideOne[j].y * (1.0 - iWeightFactor[i][j]) + sideTwo[j].y * iWeightFactor[i][j]) * weightOne[i][j] +
                                 (sideThree[i].y * (1.0 - jWeightFactor[i][j]) + sideFour[i].y * jWeightFactor[i][j]) * weightTwo[i][j];
            }
        }

        // update weights
        for (int i = 0; i < numMPoints; i++)
        {
            for (int j = 0; j < numNPoints; j++)
            {
                weightOne[i][j] = (1.0 - jWeightFactor[i][j]) * sideThreeAdimensional[i] * totalLengthThree +
                                  jWeightFactor[i][j] * sideFourAdimensional[i] * totalLengthFour;
                weightTwo[i][j] = (1.0 - iWeightFactor[i][j]) * sideOneAdimensional[j] * totalLengthOne +
                                  iWeightFactor[i][j] * sideTwoAdimensional[j] * totalLengthTwo;
            }
        }

        for (int i = 1; i < numMPoints; i++)
        {
            for (int j = 0; j < numNPoints; j++)
            {
                weightThree[i][j] = weightOne[i][j] - weightOne[i - 1][j];
            }
        }

        for (int i = 0; i < numMPoints; i++)
        {
            for (int j = 1; j < numNPoints; j++)
            {
                weightFour[i][j] = weightTwo[i][j] - weightTwo[i][j - 1];
            }
        }

        for (int i = 1; i < numMPoints; i++)
        {
            for (int j = 1; j < numNPoints - 1; j++)
            {
                weightOne[i][j] = 0.25 * (weightFour[i][j] + weightFour[i][j + 1] + weightFour[i - 1][j] + weightFour[i - 1][j + 1]) / weightThree[i][j];
            }
        }

        for (int i = 1; i < numMPoints - 1; i++)
        {
            for (int j = 1; j < numNPoints; j++)
            {
                weightTwo[i][j] = 0.25 * (weightThree[i][j] + weightThree[i][j - 1] + weightThree[i + 1][j] + weightThree[i + 1][j - 1]) / weightFour[i][j];
            }
        }

        // Iterate several times over
        const int numIterations = 25;
        for (int iter = 0; iter < numIterations; iter++)
        {
            // re-assign the weights
            for (int i = 0; i < numMPoints; i++)
            {
                for (int j = 0; j < numNPoints; j++)
                {
                    weightThree[i][j] = result[i][j].x;
                    weightFour[i][j] = result[i][j].y;
                }
            }

            for (int i = 1; i < numM; i++)
            {
                for (int j = 1; j < numN; j++)
                {

                    const double wa = 1.0 / (weightOne[i][j] + weightOne[i + 1][j] + weightTwo[i][j] + weightTwo[i][j + 1]);

                    result[i][j].x = wa * (weightThree[i - 1][j] * weightOne[i][j] + weightThree[i + 1][j] * weightOne[i + 1][j] +
                                           weightThree[i][j - 1] * weightTwo[i][j] + weightThree[i][j + 1] * weightTwo[i][j + 1]);

                    result[i][j].y = wa * (weightFour[i - 1][j] * weightOne[i][j] + weightFour[i + 1][j] * weightOne[i + 1][j] +
                                           weightFour[i][j - 1] * weightTwo[i][j] + weightFour[i][j + 1] * weightTwo[i][j + 1]);
                }
            }
        }
    }

} // namespace meshkernel
