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

#include "dilog.hpp"
#include <cmath>
#include <limits>

namespace flexiblesusy {

namespace {
   template <typename T>
   T sqr(T x) noexcept { return x*x; }
} // namespace

/**
 * @brief Real dilogarithm \f$\mathrm{Li}_2(z)\f$
 * @param x real argument
 * @note Implementation translated by R.Brun from CERNLIB DILOG function C332
 * @return \f$\mathrm{Li}_2(z)\f$
 *
 * Implemented as a truncated series expansion in terms of Chebyshev
 * polynomials, see [Yudell L. Luke: Mathematical functions and their
 * approximations, Academic Press Inc., New York 1975, p.67].
 */
double dilog(double x) noexcept {
   const double PI = M_PI;
   const double HF  = 0.5;
   const double PI2 = PI*PI;
   const double PI3 = PI2/3;
   const double PI6 = PI2/6;
   const double PI12 = PI2/12;
   const double C[20] = {0.42996693560813697, 0.40975987533077105,
     -0.01858843665014592, 0.00145751084062268,-0.00014304184442340,
      0.00001588415541880,-0.00000190784959387, 0.00000024195180854,
     -0.00000003193341274, 0.00000000434545063,-0.00000000060578480,
      0.00000000008612098,-0.00000000001244332, 0.00000000000182256,
     -0.00000000000027007, 0.00000000000004042,-0.00000000000000610,
      0.00000000000000093,-0.00000000000000014, 0.00000000000000002};

   double T,H,Y,S,A,ALFA,B1,B2,B0;

   if (x == 1) {
       H = PI6;
   } else if (x == -1) {
       H = -PI12;
   } else {
       T = -x;
       if (T <= -2) {
           Y = -1/(1+T);
           S = 1;
           B1= log(-T);
           B2= log(1+1/T);
           A = -PI3+HF*(B1*B1-B2*B2);
       } else if (T < -1) {
           Y = -1-T;
           S = -1;
           A = log(-T);
           A = -PI6+A*(A+log(1+1/T));
       } else if (T <= -0.5) {
           Y = -(1+T)/T;
           S = 1;
           A = log(-T);
           A = -PI6+A*(-HF*A+log(1+T));
       } else if (T < 0) {
           Y = -T/(1+T);
           S = -1;
           B1= log(1+T);
           A = HF*B1*B1;
       } else if (T <= 1) {
           Y = T;
           S = 1;
           A = 0;
       } else {
           Y = 1/T;
           S = -1;
           B1= log(T);
           A = PI6+HF*B1*B1;
       }
       H    = Y+Y-1;
       ALFA = H+H;
       B1   = 0;
       B2   = 0;
       for (int i=19;i>=0;i--){
          B0 = C[i] + ALFA*B1-B2;
          B2 = B1;
          B1 = B0;
       }
       H = -(S*(B0-H*B2)+A);
    }
    return H;
}

/**
 * @brief Real dilogarithm \f$\mathrm{Li}_2(z)\f$ with long double precision
 * @param x real argument
 * @note Implementation translated by R.Brun from CERNLIB DILOG function C332
 * @return \f$\mathrm{Li}_2(z)\f$
 *
 * Implemented as a truncated series expansion in terms of Chebyshev
 * polynomials, see [Yudell L. Luke: Mathematical functions and their
 * approximations, Academic Press Inc., New York 1975, p.67].
 */
long double dilog(long double x) noexcept {
   const long double PI  = 3.141592653589793238463l;
   const long double HF  = 0.5l;
   const long double PI2 = PI*PI;
   const long double PI3 = PI2/3;
   const long double PI6 = PI2/6;
   const long double PI12 = PI2/12;
   const long double C[24] = { 0.42996693560813697204l, 0.40975987533077105847l,
     -0.01858843665014591965l, 0.00145751084062267855l,-0.00014304184442340049l,
      0.00001588415541879553l,-0.00000190784959386583l, 0.00000024195180854165l,
     -0.00000003193341274252l, 0.00000000434545062677l,-0.00000000060578480118l,
      0.00000000008612097799l,-0.00000000001244331660l, 0.00000000000182255696l,
     -0.00000000000027006766l, 0.00000000000004042209l,-0.00000000000000610325l,
      0.00000000000000092863l,-0.00000000000000014226l, 0.00000000000000002193l,
      0.00000000000000000340l, 0.00000000000000000053l, 0.00000000000000000008l,
      0.00000000000000000001l};

   long double T,H,Y,S,A,ALFA,B1,B2,B0;

   if (x == 1) {
       H = PI6;
   } else if (x == -1) {
       H = -PI12;
   } else {
       T = -x;
       if (T <= -2) {
           Y = -1/(1+T);
           S = 1;
           B1= log(-T);
           B2= log(1+1/T);
           A = -PI3+HF*(B1*B1-B2*B2);
       } else if (T < -1) {
           Y = -1-T;
           S = -1;
           A = log(-T);
           A = -PI6+A*(A+log(1+1/T));
       } else if (T <= -0.5l) {
           Y = -(1+T)/T;
           S = 1;
           A = log(-T);
           A = -PI6+A*(-HF*A+log(1+T));
       } else if (T < 0) {
           Y = -T/(1+T);
           S = -1;
           B1= log(1+T);
           A = HF*B1*B1;
       } else if (T <= 1) {
           Y = T;
           S = 1;
           A = 0;
       } else {
           Y = 1/T;
           S = -1;
           B1= log(T);
           A = PI6+HF*B1*B1;
       }
       H    = Y+Y-1;
       ALFA = H+H;
       B1   = 0;
       B2   = 0;
       for (int i = 23; i >= 0; i--) {
          B0 = C[i] + ALFA*B1-B2;
          B2 = B1;
          B1 = B0;
       }
       H = -(S*(B0-H*B2)+A);
    }
    return H;
}

/**
 * @brief Complex dilogarithm \f$\mathrm{Li}_2(z)\f$
 * @param z complex argument
 * @note Implementation translated from SPheno to C++
 * @return \f$\mathrm{Li}_2(z)\f$
 */
std::complex<long double> dilog(const std::complex<long double>& z) noexcept
{
   const long double PI = 3.141592653589793238l;
   static const int N = 12;

   // bf[1..N-1] are the even Bernoulli numbers / (2 n + 1)!
   // generated by: Table[BernoulliB[2 n]/(2 n + 1)!, {n, 1, 11}]
   const long double bf[N] = {
      - 1.l/4.l                  , 1.l/36.l                 ,
      - 1.l/36.e2l               , 1.l/21168.e1l            ,
      - 1.l/108864.e2l           , 1.l/52690176.e1l         ,
      - 4.064761645144225527e-11l, 8.921691020456452555e-13l,
      - 1.993929586072107569e-14l, 4.518980029619918192e-16l,
      - 1.035651761218124701e-17l, 2.395218621026186746e-19l
   };

   const long double rz = std::real(z);
   const long double iz = std::imag(z);
   const long double nz = sqr(rz) + sqr(iz);

   // special cases
   if (iz == 0.l) {
      if (rz <= 1.l)
         return std::complex<long double>(dilog(rz), 0.l);
      if (rz > 1.l)
         return std::complex<long double>(dilog(rz), -PI*std::log(rz));
   } else if (nz < std::numeric_limits<long double>::epsilon()) {
      return z;
   }

   std::complex<long double> cy, cz;
   int jsgn, ipi12;

   // transformation to |z|<1, Re(z)<=0.5
   if (rz <= 0.5l) {
      if (nz > 1.l) {
         cy = -0.5l * sqr(std::log(-z));
         cz = -std::log(1.l - 1.l / z);
         jsgn = -1;
         ipi12 = -2;
      } else { // nz <= 1.l
         cy = 0;
         cz = -std::log(1.l - z);
         jsgn = 1;
         ipi12 = 0;
      }
   } else { // rz > 0.5
      if (nz <= 2*rz) {
         cz = -std::log(z);
         cy = cz * std::log(1.l - z);
         jsgn = -1;
         ipi12 = 2;
      } else { // nz > 2*rz
         cy = -0.5l * sqr(std::log(-z));
         cz = -std::log(1.l - 1.l / z);
         jsgn = -1;
         ipi12 = -2;
      }
   }

   // the dilogarithm
   const std::complex<long double> cz2(sqr(cz));
   const std::complex<long double> sum =
      cz +
      cz2 * (bf[0] +
      cz  * (bf[1] +
      cz2 * (bf[2] +
      cz2 * (bf[3] +
      cz2 * (bf[4] +
      cz2 * (bf[5] +
      cz2 * (bf[6] +
      cz2 * (bf[7] +
      cz2 * (bf[8] +
      cz2 * (bf[9] +
      cz2 * (bf[10] +
      cz2 * (bf[11]))))))))))));

   return static_cast<long double>(jsgn) * sum + cy + ipi12 * PI * PI / 12.l;
}

/**
 * @brief Complex dilogarithm \f$\mathrm{Li}_2(z)\f$
 * @param z complex argument
 * @note Implementation translated from SPheno to C++
 * @return \f$\mathrm{Li}_2(z)\f$
 */
std::complex<double> dilog(const std::complex<double>& z) noexcept
{
   const double PI = 3.141592653589793;
   static const int N = 10;

   // bf[1..N-1] are the even Bernoulli numbers / (2 n + 1)!
   // generated by: Table[BernoulliB[2 n]/(2 n + 1)!, {n, 1, 9}]
   const double bf[N] = {
      - 1./4.,
      + 1./36.,
      - 1./3600.,
      + 1./211680.,
      - 1./10886400.,
      + 1./526901760.,
      - 4.064761645144226e-11,
      + 8.921691020456453e-13,
      - 1.993929586072108e-14,
      + 4.518980029619918e-16
   };

   const double rz = std::real(z);
   const double iz = std::imag(z);
   const double nz = sqr(rz) + sqr(iz);

   // special cases
   if (iz == 0.) {
      if (rz <= 1.)
         return {dilog(rz), 0.};
      else // (rz > 1.)
         return {dilog(rz), -PI*std::log(rz)};
   } else if (nz < std::numeric_limits<double>::epsilon()) {
      return z;
   }

   std::complex<double> cy, cz;
   int jsgn, ipi12;

   // transformation to |z|<1, Re(z)<=0.5
   if (rz <= 0.5) {
      if (nz > 1.) {
         cy = -0.5 * sqr(std::log(-z));
         cz = -std::log(1. - 1. / z);
         jsgn = -1;
         ipi12 = -2;
      } else { // nz <= 1.
         cy = 0;
         cz = -std::log(1. - z);
         jsgn = 1;
         ipi12 = 0;
      }
   } else { // rz > 0.5
      if (nz <= 2*rz) {
         cz = -std::log(z);
         cy = cz * std::log(1. - z);
         jsgn = -1;
         ipi12 = 2;
      } else { // nz > 2*rz
         cy = -0.5 * sqr(std::log(-z));
         cz = -std::log(1. - 1. / z);
         jsgn = -1;
         ipi12 = -2;
      }
   }

   // the dilogarithm
   const std::complex<double> cz2(sqr(cz));
   const std::complex<double> sum =
      cz +
      cz2 * (bf[0] +
      cz  * (bf[1] +
      cz2 * (bf[2] +
      cz2 * (bf[3] +
      cz2 * (bf[4] +
      cz2 * (bf[5] +
      cz2 * (bf[6] +
      cz2 * (bf[7] +
      cz2 * (bf[8] +
      cz2 * (bf[9]))))))))));

   return double(jsgn) * sum + cy + ipi12 * PI * PI / 12.;
}

/**
 * @brief Clausen function \f$\mathrm{Cl}_2(\theta)\f$
 * @param x real angle
 * @return \f$\mathrm{Cl}_2(\theta)\f$
 */
double clausen_2(double x) noexcept
{
   using std::exp;
   const double PI = 3.141592653589793;
   const std::complex<double> i(0.,1.);

   while (x >= 2*PI)
      x -= 2*PI;

   while (x < 0.)
      x += 2*PI;

   if (std::abs(x) < std::numeric_limits<double>::epsilon() ||
       std::abs(x - PI) < std::numeric_limits<double>::epsilon() ||
       std::abs(x - 2*PI) < std::numeric_limits<double>::epsilon())
      return 0.;

   return std::imag(dilog(exp(i*x)));
}

} // namespace flexiblesusy
