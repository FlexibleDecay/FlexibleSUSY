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

#include "threshold_loop_functions.hpp"

#include <cmath>
#include <limits>

namespace flexiblesusy {
namespace threshold_loop_functions {

namespace {
   const double EPS = 1.e-8;
   double sqr(double x) { return x*x; }

   template <typename T>
   bool is_zero(T a, T prec = std::numeric_limits<T>::epsilon())
   {
      return std::fabs(a) < prec;
   }

   template <typename T>
   bool is_equal(T a, T b, T prec = std::numeric_limits<T>::epsilon())
   {
      return is_zero(a - b, prec);
   }
}

double F1(double x)
{
   if (is_equal(x, 1., EPS))
      return 1.;

   const double x2 = sqr(x);

   return x*std::log(x2)/(x2-1);
}

double F2(double x)
{
   if (is_equal(x, 1., EPS))
      return 1.;

   const double x2 = sqr(x);

   return 6*x2*(2-2*x2+(1+x2)*std::log(x2))/std::pow(x2-1,3);
}

double F3(double x)
{
   if (is_equal(x, 1., EPS))
      return 1.;

   const double x2 = sqr(x);

   return 2*x*(5*(1-x2)+(1+4*x2)*std::log(x2))/(3*sqr(x2-1));
}

double F4(double x)
{
   if (is_equal(x, 1., EPS))
      return 1.;

   const double x2 = sqr(x);

   return 2*x*(x2-1-std::log(x2))/sqr(x2-1);
}

double F5(double x)
{
   if (is_equal(x, 1., EPS))
      return 1.;

   const double x2 = sqr(x);
   const double x4 = std::pow(x,4);

   return 3*x*(1-x4+2*x2*std::log(x2))/std::pow(1-x2,3);
}

double F6(double x)
{
   if (is_equal(x, 1., EPS))
      return 0.;

   const double x2 = sqr(x);

   return (x2-3)/(4*(1-x2)) + x2*(x2-2)/(2*sqr(1.-x2))*std::log(x2);
}

double F7(double x)
{
   if (is_equal(x, 1., EPS))
      return 1.;

   const double x2 = sqr(x);
   const double x4 = std::pow(x,4);

   return (-3*(x4-6*x2+1.))/(2*sqr(x2-1))
      + (3*x4*(x2-3.))/(std::pow(x2-1.,3))*std::log(x2);
}

double F8(double x1, double x2)
{
   if (is_equal(x1, 1., EPS) && is_equal(x2, 1., EPS))
      return 1.;

   const double x12 = sqr(x1);
   const double x22 = sqr(x2);

   return -2. + 2./(x12-x22)
      *(std::pow(x1,4)/(x12-1.)*std::log(x12)
        -std::pow(x2,4)/(x22-1.)*std::log(x22));
}

double F9(double x1, double x2)
{
   if (is_equal(x1, 1., EPS) && is_equal(x2, 1., EPS))
      return 1.;

   const double x12 = sqr(x1);
   const double x22 = sqr(x2);

   return 2./(x12-x22)*(x12/(x12-1.)*std::log(x12)-x22/(x22-1.)*std::log(x22));
}

double f(double r)
{
   return F5(r);
}

double g(double r)
{
   return F7(r);
}

double f1(double r)
{
   if (is_equal(r, 1., EPS))
      return 1.;

   const double r2 = sqr(r);

   return (6*(r2+3)*r2)/(7*sqr(r2-1))
      + (6*(r2-5)*std::pow(r,4)*std::log(r2))/(7*std::pow(r2-1,3));
}

double f2(double r)
{
   if (is_equal(r, 1., EPS))
      return 1.;

   const double r2 = sqr(r);

   return (2*(r2+11)*r2)/(9*sqr(r2-1))
      + (2*(5*r2-17)*std::pow(r,4)*std::log(r2))/(9*std::pow(r2-1,3));
}

double f3(double r)
{
   if (is_equal(r, 1., EPS))
      return 1.;

   const double r2 = sqr(r);
   const double r4 = std::pow(r,4);

   return (2*(r4+9*r2+2))/(3*sqr(r2-1))
      + (2*(r4-7*r2-6)*r2*std::log(r2))/(3*std::pow(r2-1,3));
}

double f4(double r)
{
   if (is_equal(r, 1., EPS))
      return 1.;

   const double r2 = sqr(r);
   const double r4 = std::pow(r,4);

   return (2*(5*r4+25*r2+6))/(7*sqr(r2-1))
      + (2*(r4-19*r2-18)*r2*std::log(r2))/(7*std::pow(r2-1,3));
}

/// f5(r1,r2) in the limit r1 -> 1
static double f5_1_1(double r1, double r2)
{
   return 0.772943722943723 - 0.5524891774891774*r2 +
      0.7870670995670994*std::pow(r2,2) - 0.3316558441558441*std::pow(r2,3) +
      0.056277056277056266*std::pow(r2,4) + r1*(-0.5524891774891774 +
      1.0700757575757573*r2 - 0.6625541125541123*std::pow(r2,2) +
      0.22483766233766228*std::pow(r2,3) - 0.03344155844155843*std::pow(r2,4)) +
      std::pow(r1,3)*(-0.33165584415584404 + 0.22483766233766223*r2 -
      0.08755411255411245*std::pow(r2,2) + 0.01650432900432896*std::pow(r2,3) -
      0.0007034632034631958*std::pow(r2,4)) + std::pow(r1,4)*(0.05627705627705626
      - 0.03344155844155841*r2 + 0.010281385281385256*std::pow(r2,2) -
      0.0007034632034631921*std::pow(r2,3) -
      0.0002705627705627725*std::pow(r2,4)) + std::pow(r1,2)*(0.7870670995670994
      - 0.6625541125541123*r2 + 0.32061688311688297*std::pow(r2,2) -
      0.08755411255411248*std::pow(r2,3) + 0.01028138528138527*std::pow(r2,4));
}

/// f5(r1,r2) in the limit r1 -> 1
static double f5_1_r2(double r1, double r2)
{
   return (-0.025*std::pow(-1. + r1,3)*(4. - 17.*r2 + 4.*std::pow(r2,2) - 25.*std::pow(r2,3) -
        20.*std::pow(r2,4) + 41.*std::pow(r2,5) + 12.*std::pow(r2,6) + std::pow(r2,7) +
        (-30.*std::pow(r2,3) - 30.*std::pow(r2,5))*std::log(std::pow(r2,2))))/
    (std::pow(-1. + r2,6)*std::pow(1. + r2,2)) -
   (0.125*std::pow(-1. + r1,2)*(1. - 4.*r2 + std::pow(r2,2) - 4.*std::pow(r2,3) -
        5.*std::pow(r2,4) + 8.*std::pow(r2,5) + 3.*std::pow(r2,6) +
        (-6.*std::pow(r2,3) - 6.*std::pow(r2,5))*std::log(std::pow(r2,2))))/
    (std::pow(-1. + r2,5)*std::pow(1. + r2,2)) +
   (0.75*(-1 + r2 + 2*std::pow(r2,2) - std::pow(r2,4) - std::pow(r2,5) +
        (std::pow(r2,3) + std::pow(r2,5))*std::log(std::pow(r2,2))))/
    (std::pow(-1 + r2,3)*std::pow(1 + r2,2)) +
   (0.25*(-1. + r1)*(1. - 1.*r2 - 2.*std::pow(r2,2) + 8.*std::pow(r2,3) + std::pow(r2,4) -
        7.*std::pow(r2,5) + (3.*std::pow(r2,3) + 3.*std::pow(r2,5))*std::log(std::pow(r2,2))))/
    (std::pow(-1. + r2,4)*std::pow(1. + r2,2)) +
   (0.05*std::pow(-1. + r1,4)*(-1. + 4.5*r2 + 2.*std::pow(r2,2) + 16.5*std::pow(r2,3) -
        16.5*std::pow(r2,5) - 2.*std::pow(r2,6) - 4.5*std::pow(r2,7) + 1.*std::pow(r2,8) +
        (15.*std::pow(r2,3) + 15.*std::pow(r2,5))*std::log(std::pow(r2,2))))/
      (std::pow(-1. + r2,7)*std::pow(1. + r2,2));
}

double f5(double r1, double r2)
{
   if (is_equal(r1, 1., 0.01) && is_equal(r2, 1., 0.01))
      return f5_1_1(r1, r2);

   if (is_equal(r1, 1., 0.01))
      return f5_1_r2(r1, r2);

   if (is_equal(r2, 1., 0.01))
      return f5_1_r2(r2, r1);

   const double r12 = sqr(r1);

   const double result
      = (1+sqr(r1+r2)-r12*sqr(r2))/((r12-1)*(sqr(r2)-1))
      + (std::pow(r1,3)*(r12+1)*std::log(r12))/(sqr(r12-1)*(r1-r2))
      - (std::pow(r2,3)*(sqr(r2)+1)*std::log(sqr(r2)))/((r1-r2)*sqr(sqr(r2)-1));

   return 0.75 * result;
}

double f6(double r1, double r2)
{
   if (is_equal(r1, 1., EPS) && is_equal(r2, 1., EPS))
      return 1.;

   const double r12 = sqr(r1);
   const double r22 = sqr(r2);

   const double result
      = (r12+r22+r1*r2-r12*r22)/((r12-1)*(r22-1))
      + (std::pow(r1,5)*std::log(r12))/(sqr(r12-1)*(r1-r2))
      - (std::pow(r2,5)*std::log(r22))/((r1-r2)*sqr(r22-1));

   return 6./7. * result;
}

double f7(double r1, double r2)
{
   if (is_equal(r1, 1., EPS) && is_equal(r2, 1., EPS))
      return 1.;

   const double r12 = sqr(r1);
   const double r22 = sqr(r2);

   const double result
      = (1+r1*r2)/((r12-1)*(r22-1))
      + (std::pow(r1,3)*std::log(r12))/(sqr(r12-1)*(r1-r2))
      - (std::pow(r2,3)*std::log(r22))/((r1-r2)*sqr(r22-1));

   return 6. * result;
}

double f8(double r1, double r2)
{
   if (is_equal(r1, 1., EPS) && is_equal(r2, 1., EPS))
      return 1.;

   const double r12 = sqr(r1);
   const double r22 = sqr(r2);

   const double result
      = (r1+r2)/((r12-1)*(r22-1))
      + (std::pow(r1,4)*std::log(r12))/(sqr(r12-1)*(r1-r2))
      - (std::pow(r2,4)*std::log(r22))/((r1-r2)*sqr(r22-1));

   return 1.5 * result;
}

} // namespace threshold_loop_functions
} // namespace flexiblesusy
