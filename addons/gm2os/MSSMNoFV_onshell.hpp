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

#ifndef GM2_CALCULATOR_H
#define GM2_CALCULATOR_H

#include "MSSMNoFV_mass_eigenstates.hpp"
#include <Eigen/Core>

namespace flexiblesusy {
namespace gm2os {

class MSSMNoFV_onshell : public MSSMNoFV_mass_eigenstates {
public:
   MSSMNoFV_onshell(const MSSMNoFV_mass_eigenstates&);
   virtual ~MSSMNoFV_onshell() {}

   double get_MUDIM() const {return get_scale();}
   double get_EL0() const {return EL0;}
   double get_gY() const {return sqrt(0.6) * get_g1();}
   double get_EL() const;
   double get_TB() const {return get_vu() / get_vd();}
   Eigen::Matrix<double,3,3> get_Ae() const;
   Eigen::Matrix<double,3,3> get_Au() const;
   Eigen::Matrix<double,3,3> get_Ad() const;

   // (g-2) wrappers for pole mass getters
   double get_MW() const { return get_physical().MVWm; }
   double get_MZ() const { return get_physical().MVZ; }
   double get_ME() const { return get_physical().MFe; }
   double get_MM() const { return get_physical().MFm; }
   double get_ML() const { return get_physical().MFtau; }
   double get_MU() const { return get_physical().MFu; }
   double get_MC() const { return get_physical().MFc; }
   double get_MT() const { return get_physical().MFt; }
   double get_MD() const { return get_physical().MFd; }
   double get_MS() const { return get_physical().MFs; }
   double get_MB() const { return get_physical().MFb; }
   double get_MA0() const { return get_physical().MAh(1); }
   Eigen::Array<double,2,1> get_MSmu() const { return get_physical().MSm; }
   Eigen::Matrix<double,2,2> get_USmu() const { return get_physical().ZM; }
   Eigen::Array<double,2,1> get_MStau() const { return get_physical().MStau; }
   Eigen::Matrix<double,2,2> get_UStau() const { return get_physical().ZTau; }
   Eigen::Array<double,2,1> get_MSbot() const { return get_physical().MSb; }
   Eigen::Matrix<double,2,2> get_USbot() const { return get_physical().ZB; }
   Eigen::Array<double,2,1> get_MStop() const { return get_physical().MSt; }
   Eigen::Matrix<double,2,2> get_UStop() const { return get_physical().ZT; }

   void convert_to_onshell();

private:
   double EL0; ///< electromagnetic gauge coupling in the Thompson limit

   static bool is_equal(double, double, double);
   template <class Derived>
   static bool is_equal(const Eigen::ArrayBase<Derived>&,
                        const Eigen::ArrayBase<Derived>&, double);
   unsigned find_bino_like_neutralino();

   void convert_gauge_couplings();
   void convert_BMu();
   void convert_Mu_M1_M2(double precision_goal = 1e-5,
                         unsigned max_iterations = 100);
   void convert_vevs();
   void convert_yukawa_couplings();
};

} // namespace gm2os
} // namespace flexiblesusy

#endif