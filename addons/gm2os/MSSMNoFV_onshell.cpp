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

#include "MSSMNoFV_onshell.hpp"
#include "wrappers.hpp"
#include "logger.hpp"
#include "numerics2.hpp"

namespace flexiblesusy {
namespace gm2os {

MSSMNoFV_onshell::MSSMNoFV_onshell(const MSSMNoFV_mass_eigenstates& model_)
   : MSSMNoFV_mass_eigenstates(model_)
   , EL0(0.30282212)
{}

/**
 * Calculates the Ae parameter from the lepton trilinear and Yukawa
 * couplings.
 */
Eigen::Matrix<double,3,3> MSSMNoFV_onshell::get_Ae() const {
   Eigen::Matrix<double,3,3> Ae(get_TYe());
   Eigen::Matrix<double,3,3> Ye(get_Ye());

   for (int i = 0; i < 3; i++)
      if (!is_zero(Ye(i,i)))
         Ae(i,i) /= Ye(i,i);

   return Ae;
}

/**
 * Calculates the Au parameter from the up-type quark trilinear and
 * Yukawa couplings.
 */
Eigen::Matrix<double,3,3> MSSMNoFV_onshell::get_Au() const {
   Eigen::Matrix<double,3,3> Au(get_TYu());
   Eigen::Matrix<double,3,3> Yu(get_Yu());

   for (int i = 0; i < 3; i++)
      if (!is_zero(Yu(i,i)))
         Au(i,i) /= Yu(i,i);

   return Au;
}

/**
 * Calculates the Ad parameter from the down-type quark trilinear and
 * Yukawa couplings.
 */
Eigen::Matrix<double,3,3> MSSMNoFV_onshell::get_Ad() const {
   Eigen::Matrix<double,3,3> Ad(get_TYd());
   Eigen::Matrix<double,3,3> Yd(get_Yd());

   for (int i = 0; i < 3; i++)
      if (!is_zero(Yd(i,i)))
         Ad(i,i) /= Yd(i,i);

   return Ad;
}

/**
 * Returns the electromagnetig gauge coupling, calculated from gY and
 * g2.
 */
double MSSMNoFV_onshell::get_EL() const {
   const double gY = get_gY();
   const double g2 = get_g2();
   return gY * g2 / sqrt(sqr(gY) + sqr(g2));
}

/**
 * Converts the model parameters from the DR-bar scheme to the
 * on-shell scheme.
 *
 * The function assumes that the physical struct is filled with pole
 * masses and corresponding mixing matrices.  From these quantities,
 * the on-shell model parameters are calculated.
 *
 * @todo implement and check this function
 */
void MSSMNoFV_onshell::convert_to_onshell() {

   convert_gauge_couplings();
   convert_BMu();
   convert_vevs();
   convert_yukawa_couplings();
   convert_Mu_M1_M2();
}

void MSSMNoFV_onshell::convert_gauge_couplings()
{
   const double EL = get_EL(); // DR-bar
   const double MW = get_MW(); // pole mass
   const double MZ = get_MZ(); // pole mass
   const double cW = MW / MZ;  // on-shell weak mixing angle

   set_g1(sqrt(5. / 3.) * EL / cW);
   set_g2(EL / sqrt(1. - sqr(cW)));
}

void MSSMNoFV_onshell::convert_BMu()
{
   const double TB = get_TB(); // vu / vd
   const double MA = get_MA0(); // pole mass
   const double tan2b = 2. * TB / (1. - sqr(TB));

   set_BMu(0.5 * sqr(MA) * (tan2b / sqrt(1. + sqr(tan2b))));
}

void MSSMNoFV_onshell::convert_vevs()
{
   const double TB = get_TB(); // vu / vd
   const double MW = get_MW(); // pole mass
   const double vev = 2. * MW / get_g2();

   set_vu(vev / sqrt(1. + 1. / sqr(TB)));
   set_vd(get_vu() / TB);
}

void MSSMNoFV_onshell::convert_yukawa_couplings()
{
   Eigen::Matrix<double,3,3> Ye_neu(get_Ye());
   Ye_neu(0, 0) = sqrt(2.) * get_ME() / get_vd();
   Ye_neu(1, 1) = sqrt(2.) * get_MM() / get_vd();
   Ye_neu(2, 2) = sqrt(2.) * get_ML() / get_vd();
   set_Ye(Ye_neu);

   Eigen::Matrix<double,3,3> Yu_neu(get_Yu());
   Yu_neu(0, 0) = sqrt(2.) * get_MU() / get_vu();
   Yu_neu(1, 1) = sqrt(2.) * get_MC() / get_vu();
   Yu_neu(2, 2) = sqrt(2.) * get_MT() / get_vu();
   set_Yu(Yu_neu);

   Eigen::Matrix<double,3,3> Yd_neu(get_Yd());
   Yd_neu(0, 0) = sqrt(2.) * get_MD() / get_vd();
   Yd_neu(1, 1) = sqrt(2.) * get_MS() / get_vd();
   Yd_neu(2, 2) = sqrt(2.) * get_MB() / get_vd();
   set_Yd(Yd_neu);

   set_TYe(Ye_neu * get_Ae());
   set_TYu(Yu_neu * get_Au());
   set_TYd(Yd_neu * get_Ad());
}

template <class Derived>
bool MSSMNoFV_onshell::is_equal(const Eigen::ArrayBase<Derived>& a,
                                const Eigen::ArrayBase<Derived>& b,
                                double precision_goal)
{
   return MaxRelDiff(a,b) < precision_goal;
}

bool MSSMNoFV_onshell::is_equal(double a, double b, double precision_goal)
{
   return MaxRelDiff(a,b) < precision_goal;
}

/**
 * Returns index of most bino-like neutralino.  The function extracts
 * this information from the neutralino pole mass mixing matrix.
 */
unsigned MSSMNoFV_onshell::find_bino_like_neutralino()
{
   unsigned max_bino;
   get_physical().ZN.col(0).cwiseAbs().maxCoeff(&max_bino);

   return max_bino;
}

void MSSMNoFV_onshell::convert_Mu_M1_M2(
   double precision_goal,
   unsigned max_iterations)
{
   // find neutralino, which is most bino like
   const unsigned max_bino = find_bino_like_neutralino();

   const auto MCha_goal(get_physical().MCha);
   auto MChi_goal(get_MChi());
   MChi_goal(max_bino) = get_physical().MChi(max_bino);

   bool accuracy_goal_reached =
      MSSMNoFV_onshell::is_equal(MCha_goal, get_MCha(), precision_goal) &&
      MSSMNoFV_onshell::is_equal(MChi_goal(max_bino), get_MChi(max_bino), precision_goal);
   unsigned it = 0;

   while (!accuracy_goal_reached && it < max_iterations) {

      const auto U(get_UM()); // neg. chargino mixing matrix
      const auto V(get_UP()); // pos. chargino mixing matrix
      const auto N(get_ZN()); // neutralino mixing matrix
      const auto X(U.transpose() * MCha_goal.matrix().asDiagonal() * V);
      const auto Y(N.transpose() * MChi_goal.matrix().asDiagonal() * N);

      set_MassB(Re(Y(0,0)));
      set_MassWB(Re(X(0,0)));
      set_Mu(Re(X(1,1)));

      calculate_DRbar_masses();

      MChi_goal = get_MChi();
      MChi_goal(max_bino) = get_physical().MChi(max_bino);

      accuracy_goal_reached =
         MSSMNoFV_onshell::is_equal(MCha_goal, get_MCha(), precision_goal) &&
         MSSMNoFV_onshell::is_equal(MChi_goal(max_bino), get_MChi(max_bino), precision_goal);
      it++;
   }

   if (it == max_iterations)
      WARNING("DR-bar to on-shell conversion did not converge.");
}

} // gm2os
} // namespace flexiblesusy