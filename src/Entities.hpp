#pragma once

#include <utility>
#include <vector>

namespace GridGeom
{
    enum OperationTypes 
    {
        cartesianOperations,
        sphericalOperations
    };

    struct Point
    {
        double x;
        double y;
    };

    struct cartesian3DPoint
    {
        double x;
        double y;
        double z;
    };

    struct Nodes
    {
        std::vector<double> x;
        std::vector<double> y;
    };

    typedef std::pair<std::size_t, std::size_t> Edge;

    enum class Projections
    {
        cartesian,
        spherical
    };

}