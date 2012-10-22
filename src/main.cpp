// ====================================================================
//
// File:        main.cpp
// Author:      Alexander Voigt
// Description: main() routine
//
// This file is part of FlexibleSUSY.
//
// FlexibleSUSY is free software: you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License
// as published by the Free Software Foundation, either version 3 of
// the License, or (at your option) any later version.
//
// FlexibleSUSY is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with FlexibleSUSY.  If not, see
// <http://www.gnu.org/licenses/>.
// ====================================================================

#include "two_scale_solver.hpp"

int main()
{
   Two_scale_solver ds(0);
   ds.solve(100);

   return 0;
}
