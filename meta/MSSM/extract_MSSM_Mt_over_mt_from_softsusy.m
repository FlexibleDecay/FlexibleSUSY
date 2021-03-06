(* :Copyright:

   ====================================================================
   This file is part of FlexibleSUSY.

   FlexibleSUSY is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published
   by the Free Software Foundation, either version 3 of the License,
   or (at your option) any later version.

   FlexibleSUSY is distributed in the hope that it will be useful, but
   WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with FlexibleSUSY.  If not, see
   <http://www.gnu.org/licenses/>.
   ====================================================================

*)

(* Converts the 2-loop SUSY-QCD corrections O(alpha_s^2) to the DR-bar
   top Yukawa coupling in the MSSM from GiNaC form to Mathematica
   form.

   The GiNaC expression in "dmtas2.expr" has been extracted from
   SOFTSUSY 3.7.4 by adding the following C++ code snippet into the
   file
   src/two_loop_thresholds/tquark_corrections.cpp
   right after the cache has been filled:

   if (cache.size() == 1) {
      ofstream out("dmtas2.expr");
      out << cache[0] << endl;
   }

   Note: The expression does not include the 2-loop factor 1/(4 Pi)^4 .
 *)

str = Import["dmtas2.expr", "String"];
ex  = ToExpression[str, TraditionalForm];
ex >> "dmtas2.m"
