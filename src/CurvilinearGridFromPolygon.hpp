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

#include <vector>
#include <memory>
#include "Entities.hpp"


namespace MeshKernel
{
    // Forward declarations
    class Mesh;
    class Polygons;
    class CurvilinearGrid;

    class CurvilinearGridFromPolygon
    {
    public:

        /// <summary>
        /// Ctor
        /// </summary>
        /// <returns></returns>
        CurvilinearGridFromPolygon();

        /// <summary>
        /// 
        /// </summary>
        /// <param name="mesh">The input mesh</param>
        /// <param name="polygon">The input polygon</param>
        /// <returns></returns>
        CurvilinearGridFromPolygon(std::shared_ptr<Polygons> polygon);

        /// <summary>
        /// Compute curvilinear in a quad (pol2curvi)
        /// </summary>
        /// <returns>If the method succeeded</returns>
        bool Compute(int firstNode, int secondNode, int thirdNode, bool useFourthSide, CurvilinearGrid& curvilinearGrid) const;

    private:

        std::shared_ptr<Polygons> m_polygon;                  // A pointer to Polygons

    };

}


