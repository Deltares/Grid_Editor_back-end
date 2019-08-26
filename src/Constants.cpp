#ifndef CONSTANTS_CPP
#define CONSTANTS_CPP

#ifndef M_PI
#define M_PI (3.14159265358979323846)
#endif

namespace GridGeom
{
    static constexpr double degrad_hp = M_PI / 180.0; // conversion factor from degrees to radians(pi / 180)
    static constexpr double earth_radius = 6378137.0; // earth radius(m)
    static constexpr double dtol_pole = 0.0001;       // pole tolerance in degrees
    static constexpr double raddeg_hp = 180.0 / M_PI; // conversion factor from radians to degrees(180 / pi)
}

#endif