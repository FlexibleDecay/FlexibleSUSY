// ====================================================================
// This file is part of FlexibleSUSY.
//
// FlexibleSUSY is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published
// by the Free Software Foundation, either version 3 of the License,
// or (at your option) any later version.
//
// FlexibleSUSY is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with FlexibleSUSY.  If not, see
// <http://www.gnu.org/licenses/>.
// ====================================================================

#ifndef MSSM_TWO_SCALE_CONVERGENCE_TESTER_H
#define MSSM_TWO_SCALE_CONVERGENCE_TESTER_H

#include "two_scale_convergence_tester_skeleton.hpp"
#include "mssm_two_scale.hpp"

class Mssm_convergence_tester : public Convergence_tester_skeleton<Mssm<Two_scale> > {
public:
   Mssm_convergence_tester(Mssm<Two_scale>*, double, unsigned long);
   virtual ~Mssm_convergence_tester();
   virtual unsigned int max_iterations() const;

protected:
   virtual double max_rel_diff() const;

private:
   unsigned long maximum_iterations;

   double sumTol(const Mssm<Two_scale>&, const Mssm<Two_scale>&) const;
};

#endif
