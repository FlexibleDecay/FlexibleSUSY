// ====================================================================
// This file is part of GM2Calc.
//
// GM2Calc is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published
// by the Free Software Foundation, either version 3 of the License,
// or (at your option) any later version.
//
// GM2Calc is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with GM2Calc.  If not, see
// <http://www.gnu.org/licenses/>.
// ====================================================================

#include "MSSMNoFV_onshell.hpp"
#include "numerics2.hpp"
#include "gm2_error.hpp"
#include "gm2_1loop.hpp"
#include "ffunctions.hpp"

#include <cmath>
#include <complex>
#include <iostream>

#define WARNING(message) std::cerr << "Warning: " << message << '\n';

namespace {
   static const double ALPHA_EM_THOMPSON = 1./137.035999074;
   static const double DELTA_ALPHA_EM_MZ =
      + 0.031498 /*leptonic*/
      - 0.0000728 /*top*/
      + 0.027626 /*hadronic, arXiv:1105.3149v2 */;
   static const double ALPHA_EM_MZ =
      ALPHA_EM_THOMPSON / (1. - DELTA_ALPHA_EM_MZ);

   double calculate_e(double alpha) {
      return std::sqrt(4. * M_PI * alpha);
   }
   double calculate_alpha(double e) {
      return e * e / (4. * M_PI);
   }
}

namespace flexiblesusy {
namespace gm2calc {

MSSMNoFV_onshell::MSSMNoFV_onshell()
   : MSSMNoFV_onshell_mass_eigenstates()
   , EL(calculate_e(ALPHA_EM_MZ))
   , EL0(calculate_e(ALPHA_EM_THOMPSON))
   , Ae(Eigen::Matrix<double,3,3>::Zero())
   , Au(Eigen::Matrix<double,3,3>::Zero())
   , Ad(Eigen::Matrix<double,3,3>::Zero())
{
   get_physical().MFm  = 0.1056583715;
   get_physical().MVWm = 80.385;
   get_physical().MVZ  = 91.1876;
   set_scale(get_physical().MVZ);
}

MSSMNoFV_onshell::MSSMNoFV_onshell(const MSSMNoFV_onshell_mass_eigenstates& model_)
   : MSSMNoFV_onshell_mass_eigenstates(model_)
   , EL(calculate_e(ALPHA_EM_MZ))
   , EL0(calculate_e(ALPHA_EM_THOMPSON))
   , Ae(get_Ae())
   , Au(get_Au())
   , Ad(get_Ad())
{
}

void MSSMNoFV_onshell::set_alpha_MZ(double alpha)
{
   EL = calculate_e(alpha);
}

void MSSMNoFV_onshell::set_alpha_thompson(double alpha)
{
   EL0 = calculate_e(alpha);
}

/**
 * Returns the electromagnetig gauge coupling in the Thompson limit.
 */
double MSSMNoFV_onshell::get_EL() const {
   return EL;
}

/**
 * Converts the model parameters from the DR-bar scheme to the
 * on-shell scheme.
 *
 * The function assumes that the physical struct is filled with pole
 * masses and corresponding mixing matrices.  From these quantities,
 * the on-shell model parameters are calculated.
 */
void MSSMNoFV_onshell::convert_to_onshell(double precision) {
   calculate_DRbar_masses();
   check_input();

   convert_gauge_couplings();
   convert_BMu();
   convert_vev();
   convert_yukawa_couplings_treelevel();
   convert_Mu_M1_M2(precision, 1000);
   convert_yukawa_couplings(); // first guess of resummed yukawas
   convert_mf2(precision, 1000);
   convert_yukawa_couplings();
   calculate_DRbar_masses();
}

void MSSMNoFV_onshell::check_input()
{
   if (is_zero(get_MW()))
      throw EInvalidInput("W mass is zero");
   if (is_zero(get_MZ()))
      throw EInvalidInput("Z mass is zero");

   // if pole masses are zero, interpret the tree-level masses as pole
   // masses

#define COPY_IF_ZERO_0(m)                                               \
   if (MSSMNoFV_onshell::is_zero(get_physical().m))                     \
      get_physical().m = get_##m();
#define COPY_IF_ZERO_1(m,z)                                             \
   if (MSSMNoFV_onshell::is_zero(get_physical().m)) {                   \
      get_physical().m = get_##m();                                     \
      get_physical().z = get_##z();                                     \
   }
#define COPY_IF_ZERO_2(m,u,v)                                           \
   if (MSSMNoFV_onshell::is_zero(get_physical().m)) {                   \
      get_physical().m = get_##m();                                     \
      get_physical().u = get_##u();                                     \
      get_physical().v = get_##v();                                     \
   }

   COPY_IF_ZERO_1(MChi, ZN);
   COPY_IF_ZERO_2(MCha, UM, UP);
   COPY_IF_ZERO_0(MSveL);
   COPY_IF_ZERO_0(MSvmL);
   COPY_IF_ZERO_0(MSvtL);
   COPY_IF_ZERO_1(MSd, ZD);
   COPY_IF_ZERO_1(MSu, ZU);
   COPY_IF_ZERO_1(MSe, ZE);
   COPY_IF_ZERO_1(MSm, ZM);
   COPY_IF_ZERO_1(MStau, ZTau);
   COPY_IF_ZERO_1(MSs, ZS);
   COPY_IF_ZERO_1(MSc, ZC);
   COPY_IF_ZERO_1(MSb, ZB);
   COPY_IF_ZERO_1(MSt, ZT);

#undef COPY_IF_ZERO_0
#undef COPY_IF_ZERO_1
#undef COPY_IF_ZERO_2
}

void MSSMNoFV_onshell::convert_gauge_couplings()
{
   const double MW = get_MW(); // pole mass
   const double MZ = get_MZ(); // pole mass
   const double cW = MW / MZ;  // on-shell weak mixing angle

   set_g1(sqrt(5. / 3.) * EL / cW);
   set_g2(EL / sqrt(1. - sqr(cW)));
}

void MSSMNoFV_onshell::convert_BMu()
{
   const double TB = get_TB(); // DR-bar
   const double MA = get_MA0(); // pole mass
   const double sin2b = 2. * TB / (1. + sqr(TB));

   set_BMu(0.5 * sqr(MA) * sin2b);
}

void MSSMNoFV_onshell::convert_vev()
{
   const double TB = get_TB(); // DR-bar
   const double MW = get_MW(); // pole mass
   const double vev = 2. * MW / get_g2();

   set_vu(vev / sqrt(1. + 1. / sqr(TB)));
   set_vd(get_vu() / TB);
}

void MSSMNoFV_onshell::convert_yukawa_couplings_treelevel()
{
   Eigen::Matrix<double,3,3> Ye_neu(Eigen::Matrix<double,3,3>::Zero());
   Ye_neu(0, 0) = sqrt(2.) * get_ME() / get_vd();
   Ye_neu(1, 1) = sqrt(2.) * get_MM() / get_vd();
   Ye_neu(2, 2) = sqrt(2.) * get_ML() / get_vd();
   set_Ye(Ye_neu);

   Eigen::Matrix<double,3,3> Yu_neu(Eigen::Matrix<double,3,3>::Zero());
   Yu_neu(0, 0) = sqrt(2.) * get_MU() / get_vu();
   Yu_neu(1, 1) = sqrt(2.) * get_MC() / get_vu();
   Yu_neu(2, 2) = sqrt(2.) * get_MT() / get_vu();
   set_Yu(Yu_neu);

   Eigen::Matrix<double,3,3> Yd_neu(Eigen::Matrix<double,3,3>::Zero());
   Yd_neu(0, 0) = sqrt(2.) * get_MD() / get_vd();
   Yd_neu(1, 1) = sqrt(2.) * get_MS() / get_vd();
   Yd_neu(2, 2) = sqrt(2.) * get_MB() / get_vd();
   set_Yd(Yd_neu);

   // recalculate trilinear couplings with new Yukawas
   set_TYe(Ye_neu * Ae);
   set_TYu(Yu_neu * Au);
   set_TYd(Yd_neu * Ad);
}

/**
 * Calculate Yukawa couplings with resummed tan(beta) corrections
 *
 * @note The resummations depend on the model parameters ml2, me2, Mu,
 * MassB, MassWB, MassG.  Therefore, this routine needs to be called
 * after all these model parameters have been determined.
 */
void MSSMNoFV_onshell::convert_yukawa_couplings()
{
   Eigen::Matrix<double,3,3> Ye_neu(Eigen::Matrix<double,3,3>::Zero());
   Ye_neu(0, 0) = sqrt(2.) * get_ME() / get_vd();
   Ye_neu(1, 1) = sqrt(2.) * get_MM() / get_vd() / (1 + delta_mu_correction(*this));
   Ye_neu(2, 2) = sqrt(2.) * get_ML() / get_vd() / (1 + delta_tau_correction(*this));
   set_Ye(Ye_neu);

   Eigen::Matrix<double,3,3> Yu_neu(Eigen::Matrix<double,3,3>::Zero());
   Yu_neu(0, 0) = sqrt(2.) * get_MU() / get_vu();
   Yu_neu(1, 1) = sqrt(2.) * get_MC() / get_vu();
   Yu_neu(2, 2) = sqrt(2.) * get_MT() / get_vu();
   set_Yu(Yu_neu);

   Eigen::Matrix<double,3,3> Yd_neu(Eigen::Matrix<double,3,3>::Zero());
   Yd_neu(0, 0) = sqrt(2.) * get_MD() / get_vd();
   Yd_neu(1, 1) = sqrt(2.) * get_MS() / get_vd();
   Yd_neu(2, 2) = sqrt(2.) * get_MB() / get_vd() / (1 + delta_bottom_correction(*this));
   set_Yd(Yd_neu);

   // recalculate trilinear couplings with new Yukawas
   set_TYe(Ye_neu * Ae);
   set_TYu(Yu_neu * Au);
   set_TYd(Yd_neu * Ad);
}


template <class Derived>
bool MSSMNoFV_onshell::is_equal(const Eigen::ArrayBase<Derived>& a,
                                const Eigen::ArrayBase<Derived>& b,
                                double precision_goal)
{
   return (a - b).cwiseAbs().maxCoeff() < precision_goal;
}

bool MSSMNoFV_onshell::is_equal(double a, double b, double precision_goal)
{
   return std::abs(a - b) < precision_goal;
}

template <class Derived>
bool MSSMNoFV_onshell::is_zero(const Eigen::ArrayBase<Derived>& a,
                               double eps)
{
   return a.cwiseAbs().minCoeff() < eps;
}

template <class Derived>
bool MSSMNoFV_onshell::is_zero(const Eigen::MatrixBase<Derived>& a,
                               double eps)
{
   return a.cwiseAbs().minCoeff() < eps;
}

bool MSSMNoFV_onshell::is_zero(double a, double eps)
{
   return std::abs(a) < eps;
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

      set_MassB(std::real(Y(0,0)));
      set_MassWB(std::real(X(0,0)));
      set_Mu(std::real(X(1,1)));

      calculate_DRbar_masses();

      MChi_goal = get_MChi();
      MChi_goal(max_bino) = get_physical().MChi(max_bino);

      accuracy_goal_reached =
         MSSMNoFV_onshell::is_equal(MCha_goal, get_MCha(), precision_goal) &&
         MSSMNoFV_onshell::is_equal(MChi_goal(max_bino), get_MChi(max_bino), precision_goal);
      it++;
   }

   if (it == max_iterations)
      WARNING("DR-bar to on-shell conversion for Mu, M1 and M2 did not converge.");
}

void MSSMNoFV_onshell::convert_mf2(
   double precision_goal,
   unsigned max_iterations)
{
   const Eigen::Array<double,2,1> MSm_pole(get_physical().MSm);
   Eigen::Array<double,2,1> MSm(get_MSm());

   bool accuracy_goal_reached =
      MSSMNoFV_onshell::is_equal(MSm, MSm_pole, precision_goal);
   unsigned it = 0;

   while (!accuracy_goal_reached && it < max_iterations) {
      const Eigen::Matrix<double,2,2> ZM(get_ZM()); // smuon mixing matrix
      const Eigen::Matrix<double,2,2> M(ZM.adjoint() * MSm_pole.square().matrix().asDiagonal() * ZM);

      const double vd2 = sqr(get_vd());
      const double vu2 = sqr(get_vu());
      const double g12 = sqr(get_g1());
      const double g22 = sqr(get_g2());
      const double ymu2 = std::norm(Ye(1,1));

      const double ml211 = M(0,0)
         - (0.5*ymu2*vd2 + 0.075*g12*vd2 - 0.125*g22*vd2
            - 0.075*g12*vu2 + 0.125*g22*vu2);

      const double me211 = M(1,1)
         - (0.5*ymu2*vd2 - 0.15*g12*vd2 + 0.15*g12*vu2);

      set_ml2(1,1,ml211);
      set_me2(1,1,me211);

      calculate_DRbar_masses();

      accuracy_goal_reached =
         MSSMNoFV_onshell::is_equal(get_MSm(), MSm_pole, precision_goal);

      it++;
   }

   if (it == max_iterations)
      WARNING("DR-bar to on-shell conversion for ml2 and me2 did not converge.");
}

std::ostream& operator<<(std::ostream& os, const MSSMNoFV_onshell& model)
{
   os <<
      "======================================\n"
      " (g-2) parameters \n"
      "======================================\n"
      << "1/alpha(MZ) = " << 1./calculate_alpha(model.get_EL()) << '\n'
      << "1/alpha(0)  = " << 1./calculate_alpha(model.get_EL0()) << '\n'
      <<
      "--------------------------------------\n"
      " on-shell parameters \n"
      "--------------------------------------\n"
      "MM          = " << model.get_MM() << '\n' <<
      "MSm         = " << model.get_MSmu().transpose() << '\n' <<
      "USm         = " << model.get_USmu().row(0)
                       << model.get_USmu().row(1) << '\n' <<
      "MT          = " << model.get_MT() << '\n' <<
      "MW          = " << model.get_MW() << '\n' <<
      "MZ          = " << model.get_MZ() << '\n' <<
      "tan(beta)   = " << model.get_TB() << '\n' <<
      "yu          = " << model.get_Yu().diagonal().transpose() << '\n' <<
      "yd          = " << model.get_Yd().diagonal().transpose() << '\n' <<
      "ye          = " << model.get_Ye().diagonal().transpose() << '\n' <<
      "BMu         = " << model.get_BMu() << '\n' <<
      "Mu          = " << model.get_Mu() << '\n' <<
      "M1          = " << model.get_MassB() << '\n' <<
      "M2          = " << model.get_MassWB() << '\n' <<
      "msl(2,2)    = " << signed_abs_sqrt(model.get_ml2(1,1)) << '\n' <<
      "mse(2,2)    = " << signed_abs_sqrt(model.get_me2(1,1)) << '\n'
      ;

   return os;
}

} // gm2calc
} // namespace flexiblesusy