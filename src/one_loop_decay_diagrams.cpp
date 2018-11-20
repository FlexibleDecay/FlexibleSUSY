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

#include "one_loop_decay_diagrams.hpp"
#include "pv.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

/**
 * @brief Evaluates T1G1N1 diagram for process S -> SS
 *
 * @param[in] mext1 mass of external field S[1]
 * @param[in] mext2 mass of external field S[2]
 * @param[in] mext3 mass of external field S[3]
 * @param[in] mF4 mass of internal field F[4]
 * @param[in] mF5 mass of internal field F[5]
 * @param[in] mF6 mass of internal field F[6]
 * @param[in] CpcF4cF5S1PL coupling Cp[-F[4], -F[5], S[1]][PL]
 * @param[in] CpcF4cF5S1PR coupling Cp[-F[4], -F[5], S[1]][PR]
 * @param[in] CpF5F6S3PL coupling Cp[F[5], F[6], S[3]][PL]
 * @param[in] CpF5F6S3PR coupling Cp[F[5], F[6], S[3]][PR]
 * @param[in] CpcF6F4S2PL coupling Cp[-F[6], F[4], S[2]][PL]
 * @param[in] CpcF6F4S2PR coupling Cp[-F[6], F[4], S[2]][PR]
 *
 * @return value of the one-loop diagram
 */
Decay_amplitude_SSS calculate_diagram_SSS_t1g1n1_FFF(
   double mext1, double mext2, double mext3,
   double mF4, double mF5, double mF6,
   const std::complex<double>& CpcF4cF5S1PL, const std::complex<double>&
      CpcF4cF5S1PR, const std::complex<double>& CpF5F6S3PL, const std::complex<
      double>& CpF5F6S3PR, const std::complex<double>& CpcF6F4S2PL, const std::
      complex<double>& CpcF6F4S2PR,
   double scale)
{
   const auto b0tmp1 = passarino_veltman::B0(mext3*mext3, mF5*mF5, mF6*mF6,
      scale*scale);
   const auto c0tmp2 = passarino_veltman::C0(mext2*mext2, mext3*mext3, mext1*
      mext1, mF4*mF4, mF6*mF6, mF5*mF5, scale*scale);
   const auto c1tmp3 = passarino_veltman::C1(mext2*mext2, mext3*mext3, mext1*
      mext1, mF4*mF4, mF6*mF6, mF5*mF5, scale*scale);
   const auto c2tmp4 = passarino_veltman::C2(mext2*mext2, mext3*mext3, mext1*
      mext1, mF4*mF4, mF6*mF6, mF5*mF5, scale*scale);

   Decay_amplitude_SSS result;

   result.m_decay = mext1;
   result.m_out_1 = mext2;
   result.m_out_2 = mext3;

   result.form_factor = oneOver16PiSqr*(std::complex<double>(0,1)*(2*b0tmp1*(
      CpcF4cF5S1PR*(CpcF6F4S2PR*CpF5F6S3PL*mF4 + CpcF6F4S2PL*CpF5F6S3PR*mF5 +
      CpcF6F4S2PL*CpF5F6S3PL*mF6) + CpcF4cF5S1PL*(CpcF6F4S2PL*CpF5F6S3PR*mF4 +
      CpcF6F4S2PR*CpF5F6S3PL*mF5 + CpcF6F4S2PR*CpF5F6S3PR*mF6)) + c1tmp3*
      CpcF4cF5S1PR*CpcF6F4S2PR*CpF5F6S3PL*mF4*Sqr(mext1) + 3*c2tmp4*
      CpcF4cF5S1PR*CpcF6F4S2PR*CpF5F6S3PL*mF4*Sqr(mext1) + c1tmp3*CpcF4cF5S1PL*
      CpcF6F4S2PL*CpF5F6S3PR*mF4*Sqr(mext1) + 3*c2tmp4*CpcF4cF5S1PL*CpcF6F4S2PL
      *CpF5F6S3PR*mF4*Sqr(mext1) + c2tmp4*CpcF4cF5S1PL*CpcF6F4S2PR*CpF5F6S3PL*
      mF5*Sqr(mext1) + c2tmp4*CpcF4cF5S1PR*CpcF6F4S2PL*CpF5F6S3PR*mF5*Sqr(mext1
      ) + c1tmp3*CpcF4cF5S1PR*CpcF6F4S2PL*CpF5F6S3PL*mF6*Sqr(mext1) + 2*c2tmp4*
      CpcF4cF5S1PR*CpcF6F4S2PL*CpF5F6S3PL*mF6*Sqr(mext1) + c1tmp3*CpcF4cF5S1PL*
      CpcF6F4S2PR*CpF5F6S3PR*mF6*Sqr(mext1) + 2*c2tmp4*CpcF4cF5S1PL*CpcF6F4S2PR
      *CpF5F6S3PR*mF6*Sqr(mext1) + 3*c1tmp3*CpcF4cF5S1PR*CpcF6F4S2PR*CpF5F6S3PL
      *mF4*Sqr(mext2) + c2tmp4*CpcF4cF5S1PR*CpcF6F4S2PR*CpF5F6S3PL*mF4*Sqr(
      mext2) + 3*c1tmp3*CpcF4cF5S1PL*CpcF6F4S2PL*CpF5F6S3PR*mF4*Sqr(mext2) +
      c2tmp4*CpcF4cF5S1PL*CpcF6F4S2PL*CpF5F6S3PR*mF4*Sqr(mext2) + 2*c1tmp3*
      CpcF4cF5S1PL*CpcF6F4S2PR*CpF5F6S3PL*mF5*Sqr(mext2) + c2tmp4*CpcF4cF5S1PL*
      CpcF6F4S2PR*CpF5F6S3PL*mF5*Sqr(mext2) + 2*c1tmp3*CpcF4cF5S1PR*CpcF6F4S2PL
      *CpF5F6S3PR*mF5*Sqr(mext2) + c2tmp4*CpcF4cF5S1PR*CpcF6F4S2PL*CpF5F6S3PR*
      mF5*Sqr(mext2) + c1tmp3*CpcF4cF5S1PR*CpcF6F4S2PL*CpF5F6S3PL*mF6*Sqr(mext2
      ) + c1tmp3*CpcF4cF5S1PL*CpcF6F4S2PR*CpF5F6S3PR*mF6*Sqr(mext2) - c1tmp3*
      CpcF4cF5S1PR*CpcF6F4S2PR*CpF5F6S3PL*mF4*Sqr(mext3) - c2tmp4*CpcF4cF5S1PR*
      CpcF6F4S2PR*CpF5F6S3PL*mF4*Sqr(mext3) - c1tmp3*CpcF4cF5S1PL*CpcF6F4S2PL*
      CpF5F6S3PR*mF4*Sqr(mext3) - c2tmp4*CpcF4cF5S1PL*CpcF6F4S2PL*CpF5F6S3PR*
      mF4*Sqr(mext3) - c2tmp4*CpcF4cF5S1PL*CpcF6F4S2PR*CpF5F6S3PL*mF5*Sqr(mext3
      ) - c2tmp4*CpcF4cF5S1PR*CpcF6F4S2PL*CpF5F6S3PR*mF5*Sqr(mext3) - c1tmp3*
      CpcF4cF5S1PR*CpcF6F4S2PL*CpF5F6S3PL*mF6*Sqr(mext3) - c1tmp3*CpcF4cF5S1PL*
      CpcF6F4S2PR*CpF5F6S3PR*mF6*Sqr(mext3) + c0tmp2*mF4*(2*(CpcF4cF5S1PR*
      CpcF6F4S2PL*CpF5F6S3PL + CpcF4cF5S1PL*CpcF6F4S2PR*CpF5F6S3PR)*mF4*mF6 + 2
      *mF5*(CpcF4cF5S1PL*CpF5F6S3PL*(CpcF6F4S2PR*mF4 + CpcF6F4S2PL*mF6) +
      CpcF4cF5S1PR*CpF5F6S3PR*(CpcF6F4S2PL*mF4 + CpcF6F4S2PR*mF6)) + (
      CpcF4cF5S1PR*CpcF6F4S2PR*CpF5F6S3PL + CpcF4cF5S1PL*CpcF6F4S2PL*CpF5F6S3PR
      )*(Sqr(mext1) + Sqr(mext2) - Sqr(mext3)) + 2*(CpcF4cF5S1PR*CpcF6F4S2PR*
      CpF5F6S3PL + CpcF4cF5S1PL*CpcF6F4S2PL*CpF5F6S3PR)*Sqr(mF4))));

   return result;
}

/**
 * @brief Evaluates T1G2N2 diagram for process S -> SS
 *
 * @param[in] mext1 mass of external field S[1]
 * @param[in] mext2 mass of external field S[2]
 * @param[in] mext3 mass of external field S[3]
 * @param[in] mS4 mass of internal field S[4]
 * @param[in] mS5 mass of internal field S[5]
 * @param[in] mS6 mass of internal field S[6]
 * @param[in] CpS1cS4cS5 coupling Cp[S[1], -S[4], -S[5]]
 * @param[in] CpS2S4cS6 coupling Cp[S[2], S[4], -S[6]]
 * @param[in] CpS3S5S6 coupling Cp[S[3], S[5], S[6]]
 *
 * @return value of the one-loop diagram
 */
Decay_amplitude_SSS calculate_diagram_SSS_t1g2n2_SSS(
   double mext1, double mext2, double mext3,
   double mS4, double mS5, double mS6,
   const std::complex<double>& CpS1cS4cS5, const std::complex<double>&
      CpS2S4cS6, const std::complex<double>& CpS3S5S6,
   double scale)
{
   const auto c0tmp1 = passarino_veltman::C0(mext2*mext2, mext3*mext3, mext1*
      mext1, mS4*mS4, mS6*mS6, mS5*mS5, scale*scale);

   Decay_amplitude_SSS result;

   result.m_decay = mext1;
   result.m_out_1 = mext2;
   result.m_out_2 = mext3;

   result.form_factor = oneOver16PiSqr*(std::complex<double>(0,-1)*c0tmp1*
      CpS1cS4cS5*CpS2S4cS6*CpS3S5S6);

   return result;
}

/**
 * @brief Evaluates T1G3N3 diagram for process S -> SS
 *
 * @param[in] mext1 mass of external field S[1]
 * @param[in] mext2 mass of external field S[2]
 * @param[in] mext3 mass of external field S[3]
 * @param[in] mU4 mass of internal field U[4]
 * @param[in] mU5 mass of internal field U[5]
 * @param[in] mU6 mass of internal field U[6]
 * @param[in] CpS1cU4cU5 coupling Cp[S[1], -U[4], -U[5]]
 * @param[in] CpS2cU6U4 coupling Cp[S[2], -U[6], U[4]]
 * @param[in] CpS3U5U6 coupling Cp[S[3], U[5], U[6]]
 *
 * @return value of the one-loop diagram
 */
Decay_amplitude_SSS calculate_diagram_SSS_t1g3n3_UUU(
   double mext1, double mext2, double mext3,
   double mU4, double mU5, double mU6,
   const std::complex<double>& CpS1cU4cU5, const std::complex<double>&
      CpS2cU6U4, const std::complex<double>& CpS3U5U6,
   double scale)
{
   const auto c0tmp1 = passarino_veltman::C0(mext2*mext2, mext3*mext3, mext1*
      mext1, mU4*mU4, mU6*mU6, mU5*mU5, scale*scale);

   Decay_amplitude_SSS result;

   result.m_decay = mext1;
   result.m_out_1 = mext2;
   result.m_out_2 = mext3;

   result.form_factor = oneOver16PiSqr*(std::complex<double>(0,1)*c0tmp1*
      CpS1cU4cU5*CpS2cU6U4*CpS3U5U6);

   return result;
}

/**
 * @brief Evaluates T1G4N4 diagram for process S -> SS
 *
 * @param[in] mext1 mass of external field S[1]
 * @param[in] mext2 mass of external field S[2]
 * @param[in] mext3 mass of external field S[3]
 * @param[in] mS4 mass of internal field S[4]
 * @param[in] mS5 mass of internal field S[5]
 * @param[in] mV6 mass of internal field V[6]
 * @param[in] CpS1cS4cS5 coupling Cp[S[1], -S[4], -S[5]]
 * @param[in] CpS2S4cV6 coupling Cp[S[2], S[4], -V[6]][Mom[S[2]] - Mom[S[4]]]
 * @param[in] CpS3S5V6 coupling Cp[S[3], S[5], V[6]][Mom[S[3]] - Mom[S[5]]]
 *
 * @return value of the one-loop diagram
 */
Decay_amplitude_SSS calculate_diagram_SSS_t1g4n4_SSV(
   double mext1, double mext2, double mext3,
   double mS4, double mS5, double mV6,
   const std::complex<double>& CpS1cS4cS5, const std::complex<double>&
      CpS2S4cV6, const std::complex<double>& CpS3S5V6,
   double scale)
{
   const auto b0tmp1 = passarino_veltman::B0(mext3*mext3, mS5*mS5, mV6*mV6,
      scale*scale);
   const auto c0tmp2 = passarino_veltman::C0(mext2*mext2, mext3*mext3, mext1*
      mext1, mS4*mS4, mV6*mV6, mS5*mS5, scale*scale);
   const auto c1tmp3 = passarino_veltman::C1(mext2*mext2, mext3*mext3, mext1*
      mext1, mS4*mS4, mV6*mV6, mS5*mS5, scale*scale);
   const auto c2tmp4 = passarino_veltman::C2(mext2*mext2, mext3*mext3, mext1*
      mext1, mS4*mS4, mV6*mV6, mS5*mS5, scale*scale);

   Decay_amplitude_SSS result;

   result.m_decay = mext1;
   result.m_out_1 = mext2;
   result.m_out_2 = mext3;

   result.form_factor = oneOver16PiSqr*(std::complex<double>(0,-1)*CpS1cS4cS5*
      CpS2S4cV6*CpS3S5V6*(b0tmp1 + c1tmp3*Sqr(mext1) + c2tmp4*Sqr(mext1) -
      c1tmp3*Sqr(mext2) - c2tmp4*Sqr(mext2) - c1tmp3*Sqr(mext3) + c2tmp4*Sqr(
      mext3) + c0tmp2*(-Sqr(mext1) + Sqr(mext3) + Sqr(mS4))));

   return result;
}

/**
 * @brief Evaluates T1G5N5 diagram for process S -> SS
 *
 * @param[in] mext1 mass of external field S[1]
 * @param[in] mext2 mass of external field S[2]
 * @param[in] mext3 mass of external field S[3]
 * @param[in] mS4 mass of internal field S[4]
 * @param[in] mV5 mass of internal field V[5]
 * @param[in] mS6 mass of internal field S[6]
 * @param[in] CpS1cS4cV5 coupling Cp[S[1], -S[4], -V[5]][Mom[S[1]] - Mom[-S[4]]
    ]
 * @param[in] CpS2S4cS6 coupling Cp[S[2], S[4], -S[6]]
 * @param[in] CpS3S6V5 coupling Cp[S[3], S[6], V[5]][Mom[S[3]] - Mom[S[6]]]
 *
 * @return value of the one-loop diagram
 */
Decay_amplitude_SSS calculate_diagram_SSS_t1g5n5_SVS(
   double mext1, double mext2, double mext3,
   double mS4, double mV5, double mS6,
   const std::complex<double>& CpS1cS4cV5, const std::complex<double>&
      CpS2S4cS6, const std::complex<double>& CpS3S6V5,
   double scale)
{
   const auto b0tmp1 = passarino_veltman::B0(mext3*mext3, mS6*mS6, mV5*mV5,
      scale*scale);
   const auto c0tmp2 = passarino_veltman::C0(mext2*mext2, mext3*mext3, mext1*
      mext1, mS4*mS4, mS6*mS6, mV5*mV5, scale*scale);
   const auto c1tmp3 = passarino_veltman::C1(mext2*mext2, mext3*mext3, mext1*
      mext1, mS4*mS4, mS6*mS6, mV5*mV5, scale*scale);
   const auto c2tmp4 = passarino_veltman::C2(mext2*mext2, mext3*mext3, mext1*
      mext1, mS4*mS4, mS6*mS6, mV5*mV5, scale*scale);

   Decay_amplitude_SSS result;

   result.m_decay = mext1;
   result.m_out_1 = mext2;
   result.m_out_2 = mext3;

   result.form_factor = oneOver16PiSqr*(std::complex<double>(0,-1)*CpS1cS4cV5*
      CpS2S4cS6*CpS3S6V5*(b0tmp1 - c2tmp4*Sqr(mext1) - c0tmp2*Sqr(mext2) +
      c2tmp4*Sqr(mext2) + c0tmp2*Sqr(mext3) - c2tmp4*Sqr(mext3) + c1tmp3*(-Sqr(
      mext1) + Sqr(mext2) + Sqr(mext3)) + c0tmp2*Sqr(mS4)));

   return result;
}

/**
 * @brief Evaluates T1G6N6 diagram for process S -> SS
 *
 * @param[in] mext1 mass of external field S[1]
 * @param[in] mext2 mass of external field S[2]
 * @param[in] mext3 mass of external field S[3]
 * @param[in] mV4 mass of internal field V[4]
 * @param[in] mS5 mass of internal field S[5]
 * @param[in] mS6 mass of internal field S[6]
 * @param[in] CpS1cS5cV4 coupling Cp[S[1], -S[5], -V[4]][Mom[S[1]] - Mom[-S[5]]
    ]
 * @param[in] CpS2cS6V4 coupling Cp[S[2], -S[6], V[4]][Mom[S[2]] - Mom[-S[6]]]
 * @param[in] CpS3S5S6 coupling Cp[S[3], S[5], S[6]]
 *
 * @return value of the one-loop diagram
 */
Decay_amplitude_SSS calculate_diagram_SSS_t1g6n6_VSS(
   double mext1, double mext2, double mext3,
   double mV4, double mS5, double mS6,
   const std::complex<double>& CpS1cS5cV4, const std::complex<double>&
      CpS2cS6V4, const std::complex<double>& CpS3S5S6,
   double scale)
{
   const auto b0tmp1 = passarino_veltman::B0(mext3*mext3, mS5*mS5, mS6*mS6,
      scale*scale);
   const auto c0tmp2 = passarino_veltman::C0(mext2*mext2, mext3*mext3, mext1*
      mext1, mV4*mV4, mS6*mS6, mS5*mS5, scale*scale);
   const auto c1tmp3 = passarino_veltman::C1(mext2*mext2, mext3*mext3, mext1*
      mext1, mV4*mV4, mS6*mS6, mS5*mS5, scale*scale);
   const auto c2tmp4 = passarino_veltman::C2(mext2*mext2, mext3*mext3, mext1*
      mext1, mV4*mV4, mS6*mS6, mS5*mS5, scale*scale);

   Decay_amplitude_SSS result;

   result.m_decay = mext1;
   result.m_out_1 = mext2;
   result.m_out_2 = mext3;

   result.form_factor = oneOver16PiSqr*(std::complex<double>(0,-1)*CpS1cS5cV4*
      CpS2cS6V4*CpS3S5S6*(b0tmp1 + c1tmp3*Sqr(mext1) + 3*c2tmp4*Sqr(mext1) + 3*
      c1tmp3*Sqr(mext2) + c2tmp4*Sqr(mext2) - c1tmp3*Sqr(mext3) - c2tmp4*Sqr(
      mext3) + c0tmp2*(2*(Sqr(mext1) + Sqr(mext2) - Sqr(mext3)) + Sqr(mV4))));

   return result;
}

/**
 * @brief Evaluates T1G7N7 diagram for process S -> SS
 *
 * @param[in] mext1 mass of external field S[1]
 * @param[in] mext2 mass of external field S[2]
 * @param[in] mext3 mass of external field S[3]
 * @param[in] mS4 mass of internal field S[4]
 * @param[in] mV5 mass of internal field V[5]
 * @param[in] mV6 mass of internal field V[6]
 * @param[in] CpS1cS4cV5 coupling Cp[S[1], -S[4], -V[5]][Mom[S[1]] - Mom[-S[4]]
    ]
 * @param[in] CpS2S4cV6 coupling Cp[S[2], S[4], -V[6]][Mom[S[2]] - Mom[S[4]]]
 * @param[in] CpS3V5V6 coupling Cp[S[3], V[5], V[6]][g[lt2, lt3]]
 *
 * @return value of the one-loop diagram
 */
Decay_amplitude_SSS calculate_diagram_SSS_t1g7n7_SVV(
   double mext1, double mext2, double mext3,
   double mS4, double mV5, double mV6,
   const std::complex<double>& CpS1cS4cV5, const std::complex<double>&
      CpS2S4cV6, const std::complex<double>& CpS3V5V6,
   double scale)
{
   const auto b0tmp1 = passarino_veltman::B0(mext3*mext3, mV5*mV5, mV6*mV6,
      scale*scale);
   const auto c0tmp2 = passarino_veltman::C0(mext2*mext2, mext3*mext3, mext1*
      mext1, mS4*mS4, mV6*mV6, mV5*mV5, scale*scale);
   const auto c1tmp3 = passarino_veltman::C1(mext2*mext2, mext3*mext3, mext1*
      mext1, mS4*mS4, mV6*mV6, mV5*mV5, scale*scale);
   const auto c2tmp4 = passarino_veltman::C2(mext2*mext2, mext3*mext3, mext1*
      mext1, mS4*mS4, mV6*mV6, mV5*mV5, scale*scale);

   Decay_amplitude_SSS result;

   result.m_decay = mext1;
   result.m_out_1 = mext2;
   result.m_out_2 = mext3;

   result.form_factor = oneOver16PiSqr*(std::complex<double>(0,0.5)*CpS1cS4cV5*
      CpS2S4cV6*CpS3V5V6*(2*b0tmp1 - c1tmp3*Sqr(mext1) - 3*c2tmp4*Sqr(mext1) -
      3*c1tmp3*Sqr(mext2) - c2tmp4*Sqr(mext2) + c1tmp3*Sqr(mext3) + c2tmp4*Sqr(
      mext3) + c0tmp2*(Sqr(mext1) + Sqr(mext2) - Sqr(mext3) + 2*Sqr(mS4))));

   return result;
}

/**
 * @brief Evaluates T1G8N8 diagram for process S -> SS
 *
 * @param[in] mext1 mass of external field S[1]
 * @param[in] mext2 mass of external field S[2]
 * @param[in] mext3 mass of external field S[3]
 * @param[in] mV4 mass of internal field V[4]
 * @param[in] mS5 mass of internal field S[5]
 * @param[in] mV6 mass of internal field V[6]
 * @param[in] CpS1cS5cV4 coupling Cp[S[1], -S[5], -V[4]][Mom[S[1]] - Mom[-S[5]]
    ]
 * @param[in] CpS2V4cV6 coupling Cp[S[2], V[4], -V[6]][g[lt2, lt3]]
 * @param[in] CpS3S5V6 coupling Cp[S[3], S[5], V[6]][Mom[S[3]] - Mom[S[5]]]
 *
 * @return value of the one-loop diagram
 */
Decay_amplitude_SSS calculate_diagram_SSS_t1g8n8_VSV(
   double mext1, double mext2, double mext3,
   double mV4, double mS5, double mV6,
   const std::complex<double>& CpS1cS5cV4, const std::complex<double>&
      CpS2V4cV6, const std::complex<double>& CpS3S5V6,
   double scale)
{
   const auto b0tmp1 = passarino_veltman::B0(mext3*mext3, mS5*mS5, mV6*mV6,
      scale*scale);
   const auto c0tmp2 = passarino_veltman::C0(mext2*mext2, mext3*mext3, mext1*
      mext1, mV4*mV4, mV6*mV6, mS5*mS5, scale*scale);
   const auto c1tmp3 = passarino_veltman::C1(mext2*mext2, mext3*mext3, mext1*
      mext1, mV4*mV4, mV6*mV6, mS5*mS5, scale*scale);
   const auto c2tmp4 = passarino_veltman::C2(mext2*mext2, mext3*mext3, mext1*
      mext1, mV4*mV4, mV6*mV6, mS5*mS5, scale*scale);

   Decay_amplitude_SSS result;

   result.m_decay = mext1;
   result.m_out_1 = mext2;
   result.m_out_2 = mext3;

   result.form_factor = oneOver16PiSqr*(std::complex<double>(0,0.5)*CpS1cS5cV4*
      CpS2V4cV6*CpS3S5V6*(2*b0tmp1 + 4*c1tmp3*Sqr(mext1) + 7*c2tmp4*Sqr(mext1)
      + 2*c1tmp3*Sqr(mext2) - c2tmp4*Sqr(mext2) - 4*c1tmp3*Sqr(mext3) + c2tmp4*
      Sqr(mext3) + 2*c0tmp2*(3*Sqr(mext1) - Sqr(mext2) + Sqr(mext3) + Sqr(mV4))
      ));

   return result;
}

/**
 * @brief Evaluates T1G9N9 diagram for process S -> SS
 *
 * @param[in] mext1 mass of external field S[1]
 * @param[in] mext2 mass of external field S[2]
 * @param[in] mext3 mass of external field S[3]
 * @param[in] mV4 mass of internal field V[4]
 * @param[in] mV5 mass of internal field V[5]
 * @param[in] mS6 mass of internal field S[6]
 * @param[in] CpS1cV4cV5 coupling Cp[S[1], -V[4], -V[5]][g[lt2, lt3]]
 * @param[in] CpS2cS6V4 coupling Cp[S[2], -S[6], V[4]][Mom[S[2]] - Mom[-S[6]]]
 * @param[in] CpS3S6V5 coupling Cp[S[3], S[6], V[5]][Mom[S[3]] - Mom[S[6]]]
 *
 * @return value of the one-loop diagram
 */
Decay_amplitude_SSS calculate_diagram_SSS_t1g9n9_VVS(
   double mext1, double mext2, double mext3,
   double mV4, double mV5, double mS6,
   const std::complex<double>& CpS1cV4cV5, const std::complex<double>&
      CpS2cS6V4, const std::complex<double>& CpS3S6V5,
   double scale)
{
   const auto b0tmp1 = passarino_veltman::B0(mext3*mext3, mS6*mS6, mV5*mV5,
      scale*scale);
   const auto c0tmp2 = passarino_veltman::C0(mext2*mext2, mext3*mext3, mext1*
      mext1, mV4*mV4, mS6*mS6, mV5*mV5, scale*scale);
   const auto c1tmp3 = passarino_veltman::C1(mext2*mext2, mext3*mext3, mext1*
      mext1, mV4*mV4, mS6*mS6, mV5*mV5, scale*scale);
   const auto c2tmp4 = passarino_veltman::C2(mext2*mext2, mext3*mext3, mext1*
      mext1, mV4*mV4, mS6*mS6, mV5*mV5, scale*scale);

   Decay_amplitude_SSS result;

   result.m_decay = mext1;
   result.m_out_1 = mext2;
   result.m_out_2 = mext3;

   result.form_factor = oneOver16PiSqr*(std::complex<double>(0,0.5)*CpS1cV4cV5*
      CpS2cS6V4*CpS3S6V5*(2*b0tmp1 - c1tmp3*Sqr(mext1) + 2*c2tmp4*Sqr(mext1) +
      7*c1tmp3*Sqr(mext2) + 4*c2tmp4*Sqr(mext2) + c1tmp3*Sqr(mext3) - 4*c2tmp4*
      Sqr(mext3) + 2*c0tmp2*(-Sqr(mext1) + 3*Sqr(mext2) + Sqr(mext3) + Sqr(mV4)
      )));

   return result;
}

/**
 * @brief Evaluates T1G10N10 diagram for process S -> SS
 *
 * @param[in] mext1 mass of external field S[1]
 * @param[in] mext2 mass of external field S[2]
 * @param[in] mext3 mass of external field S[3]
 * @param[in] mV4 mass of internal field V[4]
 * @param[in] mV5 mass of internal field V[5]
 * @param[in] mV6 mass of internal field V[6]
 * @param[in] CpS1cV4cV5 coupling Cp[S[1], -V[4], -V[5]][g[lt2, lt3]]
 * @param[in] CpS2V4cV6 coupling Cp[S[2], V[4], -V[6]][g[lt2, lt3]]
 * @param[in] CpS3V5V6 coupling Cp[S[3], V[5], V[6]][g[lt2, lt3]]
 *
 * @return value of the one-loop diagram
 */
Decay_amplitude_SSS calculate_diagram_SSS_t1g10n10_VVV(
   double mext1, double mext2, double mext3,
   double mV4, double mV5, double mV6,
   const std::complex<double>& CpS1cV4cV5, const std::complex<double>&
      CpS2V4cV6, const std::complex<double>& CpS3V5V6,
   double scale)
{
   const auto c0tmp1 = passarino_veltman::C0(mext2*mext2, mext3*mext3, mext1*
      mext1, mV4*mV4, mV6*mV6, mV5*mV5, scale*scale);

   Decay_amplitude_SSS result;

   result.m_decay = mext1;
   result.m_out_1 = mext2;
   result.m_out_2 = mext3;

   result.form_factor = oneOver16PiSqr*(std::complex<double>(0,4)*c0tmp1*
      CpS1cV4cV5*CpS2V4cV6*CpS3V5V6);

   return result;
}

/**
 * @brief Evaluates T2G1N11 diagram for process S -> SS
 *
 * @param[in] mext1 mass of external field S[1]
 * @param[in] mext2 mass of external field S[2]
 * @param[in] mext3 mass of external field S[3]
 * @param[in] mS4 mass of internal field S[4]
 * @param[in] mS5 mass of internal field S[5]
 * @param[in] CpS1cS4cS5 coupling Cp[S[1], -S[4], -S[5]]
 * @param[in] CpS2S3S4S5 coupling Cp[S[2], S[3], S[4], S[5]]
 *
 * @return value of the one-loop diagram
 */
Decay_amplitude_SSS calculate_diagram_SSS_t2g1n11_SS(
   double mext1, double mext2, double mext3,
   double mS4, double mS5,
   const std::complex<double>& CpS1cS4cS5, const std::complex<double>&
      CpS2S3S4S5,
   double scale)
{
   const auto b0tmp1 = passarino_veltman::B0(mext1*mext1, mS4*mS4, mS5*mS5,
      scale*scale);

   Decay_amplitude_SSS result;

   result.m_decay = mext1;
   result.m_out_1 = mext2;
   result.m_out_2 = mext3;

   result.form_factor = oneOver16PiSqr*(-0.5*b0tmp1*CpS1cS4cS5*CpS2S3S4S5);

   return result;
}

/**
 * @brief Evaluates T2G2N12 diagram for process S -> SS
 *
 * @param[in] mext1 mass of external field S[1]
 * @param[in] mext2 mass of external field S[2]
 * @param[in] mext3 mass of external field S[3]
 * @param[in] mV4 mass of internal field V[4]
 * @param[in] mV5 mass of internal field V[5]
 * @param[in] CpS1cV4cV5 coupling Cp[S[1], -V[4], -V[5]][g[lt2, lt3]]
 * @param[in] CpS2S3V4V5 coupling Cp[S[2], S[3], V[4], V[5]][g[lt3, lt4]]
 *
 * @return value of the one-loop diagram
 */
Decay_amplitude_SSS calculate_diagram_SSS_t2g2n12_VV(
   double mext1, double mext2, double mext3,
   double mV4, double mV5,
   const std::complex<double>& CpS1cV4cV5, const std::complex<double>&
      CpS2S3V4V5,
   double scale, double finite)
{
   const auto b0tmp1 = passarino_veltman::B0(mext1*mext1, mV4*mV4, mV5*mV5,
      scale*scale);

   Decay_amplitude_SSS result;

   result.m_decay = mext1;
   result.m_out_1 = mext2;
   result.m_out_2 = mext3;

   result.form_factor = oneOver16PiSqr*(CpS1cV4cV5*CpS2S3V4V5*(-2*b0tmp1 +
      finite));

   return result;
}

/**
 * @brief Evaluates T3G1N13 diagram for process S -> SS
 *
 * @param[in] mext1 mass of external field S[1]
 * @param[in] mext2 mass of external field S[2]
 * @param[in] mext3 mass of external field S[3]
 * @param[in] mS4 mass of internal field S[4]
 * @param[in] mS5 mass of internal field S[5]
 * @param[in] CpS2cS4cS5 coupling Cp[S[2], -S[4], -S[5]]
 * @param[in] CpS1S3S4S5 coupling Cp[S[1], S[3], S[4], S[5]]
 *
 * @return value of the one-loop diagram
 */
Decay_amplitude_SSS calculate_diagram_SSS_t3g1n13_SS(
   double mext1, double mext2, double mext3,
   double mS4, double mS5,
   const std::complex<double>& CpS2cS4cS5, const std::complex<double>&
      CpS1S3S4S5,
   double scale)
{
   const auto b0tmp1 = passarino_veltman::B0(mext2*mext2, mS4*mS4, mS5*mS5,
      scale*scale);

   Decay_amplitude_SSS result;

   result.m_decay = mext1;
   result.m_out_1 = mext2;
   result.m_out_2 = mext3;

   result.form_factor = oneOver16PiSqr*(-0.5*b0tmp1*CpS1S3S4S5*CpS2cS4cS5);

   return result;
}

/**
 * @brief Evaluates T3G2N14 diagram for process S -> SS
 *
 * @param[in] mext1 mass of external field S[1]
 * @param[in] mext2 mass of external field S[2]
 * @param[in] mext3 mass of external field S[3]
 * @param[in] mV4 mass of internal field V[4]
 * @param[in] mV5 mass of internal field V[5]
 * @param[in] CpS2cV4cV5 coupling Cp[S[2], -V[4], -V[5]][g[lt2, lt3]]
 * @param[in] CpS1S3V4V5 coupling Cp[S[1], S[3], V[4], V[5]][g[lt3, lt4]]
 *
 * @return value of the one-loop diagram
 */
Decay_amplitude_SSS calculate_diagram_SSS_t3g2n14_VV(
   double mext1, double mext2, double mext3,
   double mV4, double mV5,
   const std::complex<double>& CpS2cV4cV5, const std::complex<double>&
      CpS1S3V4V5,
   double scale, double finite)
{
   const auto b0tmp1 = passarino_veltman::B0(mext2*mext2, mV4*mV4, mV5*mV5,
      scale*scale);

   Decay_amplitude_SSS result;

   result.m_decay = mext1;
   result.m_out_1 = mext2;
   result.m_out_2 = mext3;

   result.form_factor = oneOver16PiSqr*(CpS1S3V4V5*CpS2cV4cV5*(-2*b0tmp1 +
      finite));

   return result;
}

/**
 * @brief Evaluates T4G1N15 diagram for process S -> SS
 *
 * @param[in] mext1 mass of external field S[1]
 * @param[in] mext2 mass of external field S[2]
 * @param[in] mext3 mass of external field S[3]
 * @param[in] mS4 mass of internal field S[4]
 * @param[in] mS5 mass of internal field S[5]
 * @param[in] CpS3cS4cS5 coupling Cp[S[3], -S[4], -S[5]]
 * @param[in] CpS1S2S4S5 coupling Cp[S[1], S[2], S[4], S[5]]
 *
 * @return value of the one-loop diagram
 */
Decay_amplitude_SSS calculate_diagram_SSS_t4g1n15_SS(
   double mext1, double mext2, double mext3,
   double mS4, double mS5,
   const std::complex<double>& CpS3cS4cS5, const std::complex<double>&
      CpS1S2S4S5,
   double scale)
{
   const auto b0tmp1 = passarino_veltman::B0(mext3*mext3, mS4*mS4, mS5*mS5,
      scale*scale);

   Decay_amplitude_SSS result;

   result.m_decay = mext1;
   result.m_out_1 = mext2;
   result.m_out_2 = mext3;

   result.form_factor = oneOver16PiSqr*(-0.5*b0tmp1*CpS1S2S4S5*CpS3cS4cS5);

   return result;
}

/**
 * @brief Evaluates T4G2N16 diagram for process S -> SS
 *
 * @param[in] mext1 mass of external field S[1]
 * @param[in] mext2 mass of external field S[2]
 * @param[in] mext3 mass of external field S[3]
 * @param[in] mV4 mass of internal field V[4]
 * @param[in] mV5 mass of internal field V[5]
 * @param[in] CpS3cV4cV5 coupling Cp[S[3], -V[4], -V[5]][g[lt2, lt3]]
 * @param[in] CpS1S2V4V5 coupling Cp[S[1], S[2], V[4], V[5]][g[lt3, lt4]]
 *
 * @return value of the one-loop diagram
 */
Decay_amplitude_SSS calculate_diagram_SSS_t4g2n16_VV(
   double mext1, double mext2, double mext3,
   double mV4, double mV5,
   const std::complex<double>& CpS3cV4cV5, const std::complex<double>&
      CpS1S2V4V5,
   double scale, double finite)
{
   const auto b0tmp1 = passarino_veltman::B0(mext3*mext3, mV4*mV4, mV5*mV5,
      scale*scale);

   Decay_amplitude_SSS result;

   result.m_decay = mext1;
   result.m_out_1 = mext2;
   result.m_out_2 = mext3;

   result.form_factor = oneOver16PiSqr*(CpS1S2V4V5*CpS3cV4cV5*(-2*b0tmp1 +
      finite));

   return result;
}

/**
 * @brief Evaluates T1G1N1 diagram for process S -> VV
 *
 * @param[in] mext1 mass of external field S[1]
 * @param[in] mext2 mass of external field V[2]
 * @param[in] mext3 mass of external field V[3]
 * @param[in] mF4 mass of internal field F[4]
 * @param[in] mF5 mass of internal field F[5]
 * @param[in] mF6 mass of internal field F[6]
 * @param[in] CpcF4cF5S1PL coupling Cp[-F[4], -F[5], S[1]][PL]
 * @param[in] CpcF4cF5S1PR coupling Cp[-F[4], -F[5], S[1]][PR]
 * @param[in] CpF5F6V3PL coupling Cp[F[5], F[6], V[3]][LorentzProduct[gamma[lt3
    ], PL]]
 * @param[in] CpF5F6V3PR coupling Cp[F[5], F[6], V[3]][LorentzProduct[gamma[lt3
    ], PR]]
 * @param[in] CpcF6F4V2PL coupling Cp[-F[6], F[4], V[2]][LorentzProduct[gamma[
    lt3], PL]]
 * @param[in] CpcF6F4V2PR coupling Cp[-F[6], F[4], V[2]][LorentzProduct[gamma[
    lt3], PR]]
 *
 * @return value of the one-loop diagram
 */
Decay_amplitude_SVV calculate_diagram_SVV_t1g1n1_FFF(
   double mext1, double mext2, double mext3,
   double mF4, double mF5, double mF6,
   const std::complex<double>& CpcF4cF5S1PL, const std::complex<double>&
      CpcF4cF5S1PR, const std::complex<double>& CpF5F6V3PL, const std::complex<
      double>& CpF5F6V3PR, const std::complex<double>& CpcF6F4V2PL, const std::
      complex<double>& CpcF6F4V2PR,
   double scale)
{
   const auto b0tmp1 = passarino_veltman::B0(mext3*mext3, mF5*mF5, mF6*mF6,
      scale*scale);
   const auto c0tmp2 = passarino_veltman::C0(mext2*mext2, mext3*mext3, mext1*
      mext1, mF4*mF4, mF6*mF6, mF5*mF5, scale*scale);
   const auto c00tmp3 = passarino_veltman::C00(mext2*mext2, mext3*mext3, mext1*
      mext1, mF4*mF4, mF6*mF6, mF5*mF5, scale*scale);
   const auto c1tmp4 = passarino_veltman::C1(mext2*mext2, mext3*mext3, mext1*
      mext1, mF4*mF4, mF6*mF6, mF5*mF5, scale*scale);
   const auto c12tmp5 = passarino_veltman::C12(mext2*mext2, mext3*mext3, mext1*
      mext1, mF4*mF4, mF6*mF6, mF5*mF5, scale*scale);
   const auto c2tmp6 = passarino_veltman::C2(mext2*mext2, mext3*mext3, mext1*
      mext1, mF4*mF4, mF6*mF6, mF5*mF5, scale*scale);
   const auto c22tmp7 = passarino_veltman::C22(mext2*mext2, mext3*mext3, mext1*
      mext1, mF4*mF4, mF6*mF6, mF5*mF5, scale*scale);

   Decay_amplitude_SVV result;

   result.m_decay = mext1;
   result.m_out_1 = mext2;
   result.m_out_2 = mext3;

   result.form_factor_eps = oneOver16PiSqr*(-2*(c0tmp2*(CpcF4cF5S1PL*
      CpcF6F4V2PL*CpF5F6V3PL - CpcF4cF5S1PR*CpcF6F4V2PR*CpF5F6V3PR)*mF4 +
      c2tmp6*(CpcF4cF5S1PL*CpcF6F4V2PL*CpF5F6V3PL*mF4 - CpcF4cF5S1PR*
      CpcF6F4V2PR*CpF5F6V3PR*mF4 + CpcF4cF5S1PR*CpcF6F4V2PL*CpF5F6V3PL*mF5 -
      CpcF4cF5S1PL*CpcF6F4V2PR*CpF5F6V3PR*mF5) + c1tmp4*(CpcF4cF5S1PL*
      CpcF6F4V2PL*CpF5F6V3PL*mF4 - CpcF4cF5S1PR*CpcF6F4V2PR*CpF5F6V3PR*mF4 -
      CpcF4cF5S1PL*CpcF6F4V2PR*CpF5F6V3PL*mF6 + CpcF4cF5S1PR*CpcF6F4V2PL*
      CpF5F6V3PR*mF6)));
   result.form_factor_g = oneOver16PiSqr*(std::complex<double>(0,-1)*(-4*
      c00tmp3*(CpcF4cF5S1PL*CpcF6F4V2PL*CpF5F6V3PL*mF4 + CpcF4cF5S1PR*
      CpcF6F4V2PR*CpF5F6V3PR*mF4 + CpcF4cF5S1PR*CpcF6F4V2PL*CpF5F6V3PL*mF5 +
      CpcF4cF5S1PL*CpcF6F4V2PR*CpF5F6V3PR*mF5) - 2*c0tmp2*CpcF4cF5S1PR*
      CpcF6F4V2PR*CpF5F6V3PL*mF4*mF5*mF6 - 2*c0tmp2*CpcF4cF5S1PL*CpcF6F4V2PL*
      CpF5F6V3PR*mF4*mF5*mF6 + 2*b0tmp1*(CpcF4cF5S1PL*(CpcF6F4V2PL*CpF5F6V3PL*
      mF4 + CpcF6F4V2PR*CpF5F6V3PR*mF5 - CpcF6F4V2PR*CpF5F6V3PL*mF6) +
      CpcF4cF5S1PR*(CpcF6F4V2PR*CpF5F6V3PR*mF4 + CpcF6F4V2PL*CpF5F6V3PL*mF5 -
      CpcF6F4V2PL*CpF5F6V3PR*mF6)) + 2*c0tmp2*CpcF4cF5S1PL*CpcF6F4V2PL*
      CpF5F6V3PL*Cube(mF4) + 2*c0tmp2*CpcF4cF5S1PR*CpcF6F4V2PR*CpF5F6V3PR*Cube(
      mF4) + c0tmp2*CpcF4cF5S1PL*CpcF6F4V2PL*CpF5F6V3PL*mF4*Sqr(mext1) + c1tmp4
      *CpcF4cF5S1PL*CpcF6F4V2PL*CpF5F6V3PL*mF4*Sqr(mext1) + 3*c2tmp6*
      CpcF4cF5S1PL*CpcF6F4V2PL*CpF5F6V3PL*mF4*Sqr(mext1) + c0tmp2*CpcF4cF5S1PR*
      CpcF6F4V2PR*CpF5F6V3PR*mF4*Sqr(mext1) + c1tmp4*CpcF4cF5S1PR*CpcF6F4V2PR*
      CpF5F6V3PR*mF4*Sqr(mext1) + 3*c2tmp6*CpcF4cF5S1PR*CpcF6F4V2PR*CpF5F6V3PR*
      mF4*Sqr(mext1) + c2tmp6*CpcF4cF5S1PR*CpcF6F4V2PL*CpF5F6V3PL*mF5*Sqr(mext1
      ) + c2tmp6*CpcF4cF5S1PL*CpcF6F4V2PR*CpF5F6V3PR*mF5*Sqr(mext1) - c1tmp4*
      CpcF4cF5S1PL*CpcF6F4V2PR*CpF5F6V3PL*mF6*Sqr(mext1) - 2*c2tmp6*
      CpcF4cF5S1PL*CpcF6F4V2PR*CpF5F6V3PL*mF6*Sqr(mext1) - c1tmp4*CpcF4cF5S1PR*
      CpcF6F4V2PL*CpF5F6V3PR*mF6*Sqr(mext1) - 2*c2tmp6*CpcF4cF5S1PR*CpcF6F4V2PL
      *CpF5F6V3PR*mF6*Sqr(mext1) + c0tmp2*CpcF4cF5S1PL*CpcF6F4V2PL*CpF5F6V3PL*
      mF4*Sqr(mext2) + 3*c1tmp4*CpcF4cF5S1PL*CpcF6F4V2PL*CpF5F6V3PL*mF4*Sqr(
      mext2) + c2tmp6*CpcF4cF5S1PL*CpcF6F4V2PL*CpF5F6V3PL*mF4*Sqr(mext2) +
      c0tmp2*CpcF4cF5S1PR*CpcF6F4V2PR*CpF5F6V3PR*mF4*Sqr(mext2) + 3*c1tmp4*
      CpcF4cF5S1PR*CpcF6F4V2PR*CpF5F6V3PR*mF4*Sqr(mext2) + c2tmp6*CpcF4cF5S1PR*
      CpcF6F4V2PR*CpF5F6V3PR*mF4*Sqr(mext2) + 2*c1tmp4*CpcF4cF5S1PR*CpcF6F4V2PL
      *CpF5F6V3PL*mF5*Sqr(mext2) + c2tmp6*CpcF4cF5S1PR*CpcF6F4V2PL*CpF5F6V3PL*
      mF5*Sqr(mext2) + 2*c1tmp4*CpcF4cF5S1PL*CpcF6F4V2PR*CpF5F6V3PR*mF5*Sqr(
      mext2) + c2tmp6*CpcF4cF5S1PL*CpcF6F4V2PR*CpF5F6V3PR*mF5*Sqr(mext2) -
      c1tmp4*CpcF4cF5S1PL*CpcF6F4V2PR*CpF5F6V3PL*mF6*Sqr(mext2) - c1tmp4*
      CpcF4cF5S1PR*CpcF6F4V2PL*CpF5F6V3PR*mF6*Sqr(mext2) - c0tmp2*CpcF4cF5S1PL*
      CpcF6F4V2PL*CpF5F6V3PL*mF4*Sqr(mext3) - c1tmp4*CpcF4cF5S1PL*CpcF6F4V2PL*
      CpF5F6V3PL*mF4*Sqr(mext3) - c2tmp6*CpcF4cF5S1PL*CpcF6F4V2PL*CpF5F6V3PL*
      mF4*Sqr(mext3) - c0tmp2*CpcF4cF5S1PR*CpcF6F4V2PR*CpF5F6V3PR*mF4*Sqr(mext3
      ) - c1tmp4*CpcF4cF5S1PR*CpcF6F4V2PR*CpF5F6V3PR*mF4*Sqr(mext3) - c2tmp6*
      CpcF4cF5S1PR*CpcF6F4V2PR*CpF5F6V3PR*mF4*Sqr(mext3) - c2tmp6*CpcF4cF5S1PR*
      CpcF6F4V2PL*CpF5F6V3PL*mF5*Sqr(mext3) - c2tmp6*CpcF4cF5S1PL*CpcF6F4V2PR*
      CpF5F6V3PR*mF5*Sqr(mext3) + c1tmp4*CpcF4cF5S1PL*CpcF6F4V2PR*CpF5F6V3PL*
      mF6*Sqr(mext3) + c1tmp4*CpcF4cF5S1PR*CpcF6F4V2PL*CpF5F6V3PR*mF6*Sqr(mext3
      ) + 2*c0tmp2*CpcF4cF5S1PR*CpcF6F4V2PL*CpF5F6V3PL*mF5*Sqr(mF4) + 2*c0tmp2*
      CpcF4cF5S1PL*CpcF6F4V2PR*CpF5F6V3PR*mF5*Sqr(mF4) - 2*c0tmp2*CpcF4cF5S1PL*
      CpcF6F4V2PR*CpF5F6V3PL*mF6*Sqr(mF4) - 2*c0tmp2*CpcF4cF5S1PR*CpcF6F4V2PL*
      CpF5F6V3PR*mF6*Sqr(mF4)));
   result.form_factor_11 = oneOver16PiSqr*(std::complex<double>(0,2)*(c1tmp4*
      CpcF4cF5S1PL*CpcF6F4V2PL*CpF5F6V3PL*mF4 + 2*c22tmp7*CpcF4cF5S1PL*
      CpcF6F4V2PL*CpF5F6V3PL*mF4 + 3*c2tmp6*CpcF4cF5S1PL*CpcF6F4V2PL*CpF5F6V3PL
      *mF4 + c1tmp4*CpcF4cF5S1PR*CpcF6F4V2PR*CpF5F6V3PR*mF4 + 2*c22tmp7*
      CpcF4cF5S1PR*CpcF6F4V2PR*CpF5F6V3PR*mF4 + 3*c2tmp6*CpcF4cF5S1PR*
      CpcF6F4V2PR*CpF5F6V3PR*mF4 + c0tmp2*(CpcF4cF5S1PL*CpcF6F4V2PL*CpF5F6V3PL
      + CpcF4cF5S1PR*CpcF6F4V2PR*CpF5F6V3PR)*mF4 + 2*c22tmp7*CpcF4cF5S1PR*
      CpcF6F4V2PL*CpF5F6V3PL*mF5 + c2tmp6*CpcF4cF5S1PR*CpcF6F4V2PL*CpF5F6V3PL*
      mF5 + 2*c22tmp7*CpcF4cF5S1PL*CpcF6F4V2PR*CpF5F6V3PR*mF5 + c2tmp6*
      CpcF4cF5S1PL*CpcF6F4V2PR*CpF5F6V3PR*mF5 + 2*c12tmp5*(CpcF4cF5S1PL*
      CpcF6F4V2PL*CpF5F6V3PL*mF4 + CpcF4cF5S1PR*CpcF6F4V2PR*CpF5F6V3PR*mF4 +
      CpcF4cF5S1PR*CpcF6F4V2PL*CpF5F6V3PL*mF5 + CpcF4cF5S1PL*CpcF6F4V2PR*
      CpF5F6V3PR*mF5) - c1tmp4*CpcF4cF5S1PL*CpcF6F4V2PR*CpF5F6V3PL*mF6 - c1tmp4
      *CpcF4cF5S1PR*CpcF6F4V2PL*CpF5F6V3PR*mF6));
   result.form_factor_21 = oneOver16PiSqr*(std::complex<double>(0,2)*(c1tmp4*
      CpcF4cF5S1PL*CpcF6F4V2PL*CpF5F6V3PL*mF4 + 2*c22tmp7*CpcF4cF5S1PL*
      CpcF6F4V2PL*CpF5F6V3PL*mF4 + 3*c2tmp6*CpcF4cF5S1PL*CpcF6F4V2PL*CpF5F6V3PL
      *mF4 + c1tmp4*CpcF4cF5S1PR*CpcF6F4V2PR*CpF5F6V3PR*mF4 + 2*c22tmp7*
      CpcF4cF5S1PR*CpcF6F4V2PR*CpF5F6V3PR*mF4 + 3*c2tmp6*CpcF4cF5S1PR*
      CpcF6F4V2PR*CpF5F6V3PR*mF4 + c0tmp2*(CpcF4cF5S1PL*CpcF6F4V2PL*CpF5F6V3PL
      + CpcF4cF5S1PR*CpcF6F4V2PR*CpF5F6V3PR)*mF4 + 2*c22tmp7*CpcF4cF5S1PR*
      CpcF6F4V2PL*CpF5F6V3PL*mF5 + c2tmp6*CpcF4cF5S1PR*CpcF6F4V2PL*CpF5F6V3PL*
      mF5 + 2*c22tmp7*CpcF4cF5S1PL*CpcF6F4V2PR*CpF5F6V3PR*mF5 + c2tmp6*
      CpcF4cF5S1PL*CpcF6F4V2PR*CpF5F6V3PR*mF5 + 2*c12tmp5*(CpcF4cF5S1PL*
      CpcF6F4V2PL*CpF5F6V3PL*mF4 + CpcF4cF5S1PR*CpcF6F4V2PR*CpF5F6V3PR*mF4 +
      CpcF4cF5S1PR*CpcF6F4V2PL*CpF5F6V3PL*mF5 + CpcF4cF5S1PL*CpcF6F4V2PR*
      CpF5F6V3PR*mF5) - c1tmp4*CpcF4cF5S1PL*CpcF6F4V2PR*CpF5F6V3PL*mF6 - c1tmp4
      *CpcF4cF5S1PR*CpcF6F4V2PL*CpF5F6V3PR*mF6));
   result.form_factor_12 = oneOver16PiSqr*(std::complex<double>(0,2)*(c1tmp4*
      CpcF4cF5S1PL*CpcF6F4V2PL*CpF5F6V3PL*mF4 + 2*c22tmp7*CpcF4cF5S1PL*
      CpcF6F4V2PL*CpF5F6V3PL*mF4 + 3*c2tmp6*CpcF4cF5S1PL*CpcF6F4V2PL*CpF5F6V3PL
      *mF4 + c1tmp4*CpcF4cF5S1PR*CpcF6F4V2PR*CpF5F6V3PR*mF4 + 2*c22tmp7*
      CpcF4cF5S1PR*CpcF6F4V2PR*CpF5F6V3PR*mF4 + 3*c2tmp6*CpcF4cF5S1PR*
      CpcF6F4V2PR*CpF5F6V3PR*mF4 + c0tmp2*(CpcF4cF5S1PL*CpcF6F4V2PL*CpF5F6V3PL
      + CpcF4cF5S1PR*CpcF6F4V2PR*CpF5F6V3PR)*mF4 + 2*c22tmp7*CpcF4cF5S1PR*
      CpcF6F4V2PL*CpF5F6V3PL*mF5 + c2tmp6*CpcF4cF5S1PR*CpcF6F4V2PL*CpF5F6V3PL*
      mF5 + 2*c22tmp7*CpcF4cF5S1PL*CpcF6F4V2PR*CpF5F6V3PR*mF5 + c2tmp6*
      CpcF4cF5S1PL*CpcF6F4V2PR*CpF5F6V3PR*mF5 + 2*c12tmp5*(CpcF4cF5S1PL*
      CpcF6F4V2PL*CpF5F6V3PL*mF4 + CpcF4cF5S1PR*CpcF6F4V2PR*CpF5F6V3PR*mF4 +
      CpcF4cF5S1PR*CpcF6F4V2PL*CpF5F6V3PL*mF5 + CpcF4cF5S1PL*CpcF6F4V2PR*
      CpF5F6V3PR*mF5) - c1tmp4*CpcF4cF5S1PL*CpcF6F4V2PR*CpF5F6V3PL*mF6 - c1tmp4
      *CpcF4cF5S1PR*CpcF6F4V2PL*CpF5F6V3PR*mF6));
   result.form_factor_22 = oneOver16PiSqr*(std::complex<double>(0,2)*(c1tmp4*
      CpcF4cF5S1PL*CpcF6F4V2PL*CpF5F6V3PL*mF4 + 2*c22tmp7*CpcF4cF5S1PL*
      CpcF6F4V2PL*CpF5F6V3PL*mF4 + 3*c2tmp6*CpcF4cF5S1PL*CpcF6F4V2PL*CpF5F6V3PL
      *mF4 + c1tmp4*CpcF4cF5S1PR*CpcF6F4V2PR*CpF5F6V3PR*mF4 + 2*c22tmp7*
      CpcF4cF5S1PR*CpcF6F4V2PR*CpF5F6V3PR*mF4 + 3*c2tmp6*CpcF4cF5S1PR*
      CpcF6F4V2PR*CpF5F6V3PR*mF4 + c0tmp2*(CpcF4cF5S1PL*CpcF6F4V2PL*CpF5F6V3PL
      + CpcF4cF5S1PR*CpcF6F4V2PR*CpF5F6V3PR)*mF4 + 2*c22tmp7*CpcF4cF5S1PR*
      CpcF6F4V2PL*CpF5F6V3PL*mF5 + c2tmp6*CpcF4cF5S1PR*CpcF6F4V2PL*CpF5F6V3PL*
      mF5 + 2*c22tmp7*CpcF4cF5S1PL*CpcF6F4V2PR*CpF5F6V3PR*mF5 + c2tmp6*
      CpcF4cF5S1PL*CpcF6F4V2PR*CpF5F6V3PR*mF5 + 2*c12tmp5*(CpcF4cF5S1PL*
      CpcF6F4V2PL*CpF5F6V3PL*mF4 + CpcF4cF5S1PR*CpcF6F4V2PR*CpF5F6V3PR*mF4 +
      CpcF4cF5S1PR*CpcF6F4V2PL*CpF5F6V3PL*mF5 + CpcF4cF5S1PL*CpcF6F4V2PR*
      CpF5F6V3PR*mF5) - c1tmp4*CpcF4cF5S1PL*CpcF6F4V2PR*CpF5F6V3PL*mF6 - c1tmp4
      *CpcF4cF5S1PR*CpcF6F4V2PL*CpF5F6V3PR*mF6));

   return result;
}

/**
 * @brief Evaluates T1G2N2 diagram for process S -> VV
 *
 * @param[in] mext1 mass of external field S[1]
 * @param[in] mext2 mass of external field V[2]
 * @param[in] mext3 mass of external field V[3]
 * @param[in] mS4 mass of internal field S[4]
 * @param[in] mS5 mass of internal field S[5]
 * @param[in] mS6 mass of internal field S[6]
 * @param[in] CpS1cS4cS5 coupling Cp[S[1], -S[4], -S[5]]
 * @param[in] CpS4cS6V2 coupling Cp[S[4], -S[6], V[2]][Mom[S[4]] - Mom[-S[6]]]
 * @param[in] CpS5S6V3 coupling Cp[S[5], S[6], V[3]][Mom[S[5]] - Mom[S[6]]]
 *
 * @return value of the one-loop diagram
 */
Decay_amplitude_SVV calculate_diagram_SVV_t1g2n2_SSS(
   double mext1, double mext2, double mext3,
   double mS4, double mS5, double mS6,
   const std::complex<double>& CpS1cS4cS5, const std::complex<double>&
      CpS4cS6V2, const std::complex<double>& CpS5S6V3,
   double scale)
{
   const auto c00tmp1 = passarino_veltman::C00(mext2*mext2, mext3*mext3, mext1*
      mext1, mS4*mS4, mS6*mS6, mS5*mS5, scale*scale);
   const auto c12tmp2 = passarino_veltman::C12(mext2*mext2, mext3*mext3, mext1*
      mext1, mS4*mS4, mS6*mS6, mS5*mS5, scale*scale);
   const auto c2tmp3 = passarino_veltman::C2(mext2*mext2, mext3*mext3, mext1*
      mext1, mS4*mS4, mS6*mS6, mS5*mS5, scale*scale);
   const auto c22tmp4 = passarino_veltman::C22(mext2*mext2, mext3*mext3, mext1*
      mext1, mS4*mS4, mS6*mS6, mS5*mS5, scale*scale);

   Decay_amplitude_SVV result;

   result.m_decay = mext1;
   result.m_out_1 = mext2;
   result.m_out_2 = mext3;

   result.form_factor_g = oneOver16PiSqr*(std::complex<double>(0,4)*c00tmp1*
      CpS1cS4cS5*CpS4cS6V2*CpS5S6V3);
   result.form_factor_11 = oneOver16PiSqr*(std::complex<double>(0,4)*(c12tmp2 +
      c22tmp4 + c2tmp3)*CpS1cS4cS5*CpS4cS6V2*CpS5S6V3);
   result.form_factor_21 = oneOver16PiSqr*(std::complex<double>(0,4)*(c12tmp2 +
      c22tmp4 + c2tmp3)*CpS1cS4cS5*CpS4cS6V2*CpS5S6V3);
   result.form_factor_12 = oneOver16PiSqr*(std::complex<double>(0,4)*(c12tmp2 +
      c22tmp4 + c2tmp3)*CpS1cS4cS5*CpS4cS6V2*CpS5S6V3);
   result.form_factor_22 = oneOver16PiSqr*(std::complex<double>(0,4)*(c12tmp2 +
      c22tmp4 + c2tmp3)*CpS1cS4cS5*CpS4cS6V2*CpS5S6V3);

   return result;
}

/**
 * @brief Evaluates T1G3N3 diagram for process S -> VV
 *
 * @param[in] mext1 mass of external field S[1]
 * @param[in] mext2 mass of external field V[2]
 * @param[in] mext3 mass of external field V[3]
 * @param[in] mU4 mass of internal field U[4]
 * @param[in] mU5 mass of internal field U[5]
 * @param[in] mU6 mass of internal field U[6]
 * @param[in] CpS1cU4cU5 coupling Cp[S[1], -U[4], -U[5]]
 * @param[in] CpU5U6V3 coupling Cp[U[5], U[6], V[3]][Mom[U[5]]]
 * @param[in] CpcU6U4V2 coupling Cp[-U[6], U[4], V[2]][Mom[-U[6]]]
 *
 * @return value of the one-loop diagram
 */
Decay_amplitude_SVV calculate_diagram_SVV_t1g3n3_UUU(
   double mext1, double mext2, double mext3,
   double mU4, double mU5, double mU6,
   const std::complex<double>& CpS1cU4cU5, const std::complex<double>& CpU5U6V3
      , const std::complex<double>& CpcU6U4V2,
   double scale)
{
   const auto c00tmp1 = passarino_veltman::C00(mext2*mext2, mext3*mext3, mext1*
      mext1, mU4*mU4, mU6*mU6, mU5*mU5, scale*scale);
   const auto c12tmp2 = passarino_veltman::C12(mext2*mext2, mext3*mext3, mext1*
      mext1, mU4*mU4, mU6*mU6, mU5*mU5, scale*scale);
   const auto c2tmp3 = passarino_veltman::C2(mext2*mext2, mext3*mext3, mext1*
      mext1, mU4*mU4, mU6*mU6, mU5*mU5, scale*scale);
   const auto c22tmp4 = passarino_veltman::C22(mext2*mext2, mext3*mext3, mext1*
      mext1, mU4*mU4, mU6*mU6, mU5*mU5, scale*scale);

   Decay_amplitude_SVV result;

   result.m_decay = mext1;
   result.m_out_1 = mext2;
   result.m_out_2 = mext3;

   result.form_factor_g = oneOver16PiSqr*(std::complex<double>(0,1)*c00tmp1*
      CpcU6U4V2*CpS1cU4cU5*CpU5U6V3);
   result.form_factor_11 = oneOver16PiSqr*(std::complex<double>(0,1)*(c12tmp2 +
      c22tmp4 + c2tmp3)*CpcU6U4V2*CpS1cU4cU5*CpU5U6V3);
   result.form_factor_21 = oneOver16PiSqr*(std::complex<double>(0,1)*(c12tmp2 +
      c22tmp4 + c2tmp3)*CpcU6U4V2*CpS1cU4cU5*CpU5U6V3);
   result.form_factor_12 = oneOver16PiSqr*(std::complex<double>(0,1)*(c12tmp2 +
      c22tmp4 + c2tmp3)*CpcU6U4V2*CpS1cU4cU5*CpU5U6V3);
   result.form_factor_22 = oneOver16PiSqr*(std::complex<double>(0,1)*(c12tmp2 +
      c22tmp4 + c2tmp3)*CpcU6U4V2*CpS1cU4cU5*CpU5U6V3);

   return result;
}

/**
 * @brief Evaluates T1G4N4 diagram for process S -> VV
 *
 * @param[in] mext1 mass of external field S[1]
 * @param[in] mext2 mass of external field V[2]
 * @param[in] mext3 mass of external field V[3]
 * @param[in] mS4 mass of internal field S[4]
 * @param[in] mS5 mass of internal field S[5]
 * @param[in] mV6 mass of internal field V[6]
 * @param[in] CpS1cS4cS5 coupling Cp[S[1], -S[4], -S[5]]
 * @param[in] CpS4V2cV6 coupling Cp[S[4], V[2], -V[6]][g[lt2, lt3]]
 * @param[in] CpS5V3V6 coupling Cp[S[5], V[3], V[6]][g[lt2, lt3]]
 *
 * @return value of the one-loop diagram
 */
Decay_amplitude_SVV calculate_diagram_SVV_t1g4n4_SSV(
   double mext1, double mext2, double mext3,
   double mS4, double mS5, double mV6,
   const std::complex<double>& CpS1cS4cS5, const std::complex<double>&
      CpS4V2cV6, const std::complex<double>& CpS5V3V6,
   double scale)
{
   const auto c0tmp1 = passarino_veltman::C0(mext2*mext2, mext3*mext3, mext1*
      mext1, mS4*mS4, mV6*mV6, mS5*mS5, scale*scale);

   Decay_amplitude_SVV result;

   result.m_decay = mext1;
   result.m_out_1 = mext2;
   result.m_out_2 = mext3;

   result.form_factor_g = oneOver16PiSqr*(std::complex<double>(0,1)*c0tmp1*
      CpS1cS4cS5*CpS4V2cV6*CpS5V3V6);

   return result;
}

/**
 * @brief Evaluates T1G5N5 diagram for process S -> VV
 *
 * @param[in] mext1 mass of external field S[1]
 * @param[in] mext2 mass of external field V[2]
 * @param[in] mext3 mass of external field V[3]
 * @param[in] mS4 mass of internal field S[4]
 * @param[in] mV5 mass of internal field V[5]
 * @param[in] mS6 mass of internal field S[6]
 * @param[in] CpS1cS4cV5 coupling Cp[S[1], -S[4], -V[5]][Mom[S[1]] - Mom[-S[4]]
    ]
 * @param[in] CpS4cS6V2 coupling Cp[S[4], -S[6], V[2]][Mom[S[4]] - Mom[-S[6]]]
 * @param[in] CpS6V3V5 coupling Cp[S[6], V[3], V[5]][g[lt2, lt3]]
 *
 * @return value of the one-loop diagram
 */
Decay_amplitude_SVV calculate_diagram_SVV_t1g5n5_SVS(
   double mext1, double mext2, double mext3,
   double mS4, double mV5, double mS6,
   const std::complex<double>& CpS1cS4cV5, const std::complex<double>&
      CpS4cS6V2, const std::complex<double>& CpS6V3V5,
   double scale)
{
   const auto c00tmp1 = passarino_veltman::C00(mext2*mext2, mext3*mext3, mext1*
      mext1, mS4*mS4, mS6*mS6, mV5*mV5, scale*scale);
   const auto c12tmp2 = passarino_veltman::C12(mext2*mext2, mext3*mext3, mext1*
      mext1, mS4*mS4, mS6*mS6, mV5*mV5, scale*scale);
   const auto c2tmp3 = passarino_veltman::C2(mext2*mext2, mext3*mext3, mext1*
      mext1, mS4*mS4, mS6*mS6, mV5*mV5, scale*scale);
   const auto c22tmp4 = passarino_veltman::C22(mext2*mext2, mext3*mext3, mext1*
      mext1, mS4*mS4, mS6*mS6, mV5*mV5, scale*scale);

   Decay_amplitude_SVV result;

   result.m_decay = mext1;
   result.m_out_1 = mext2;
   result.m_out_2 = mext3;

   result.form_factor_g = oneOver16PiSqr*(std::complex<double>(0,2)*c00tmp1*
      CpS1cS4cV5*CpS4cS6V2*CpS6V3V5);
   result.form_factor_11 = oneOver16PiSqr*(std::complex<double>(0,2)*(c12tmp2 +
      c22tmp4 - c2tmp3)*CpS1cS4cV5*CpS4cS6V2*CpS6V3V5);
   result.form_factor_21 = oneOver16PiSqr*(std::complex<double>(0,2)*(c12tmp2 +
      c22tmp4 - c2tmp3)*CpS1cS4cV5*CpS4cS6V2*CpS6V3V5);
   result.form_factor_12 = oneOver16PiSqr*(std::complex<double>(0,2)*(c12tmp2 +
      c22tmp4 - c2tmp3)*CpS1cS4cV5*CpS4cS6V2*CpS6V3V5);
   result.form_factor_22 = oneOver16PiSqr*(std::complex<double>(0,2)*(c12tmp2 +
      c22tmp4 - c2tmp3)*CpS1cS4cV5*CpS4cS6V2*CpS6V3V5);

   return result;
}

/**
 * @brief Evaluates T1G6N6 diagram for process S -> VV
 *
 * @param[in] mext1 mass of external field S[1]
 * @param[in] mext2 mass of external field V[2]
 * @param[in] mext3 mass of external field V[3]
 * @param[in] mV4 mass of internal field V[4]
 * @param[in] mS5 mass of internal field S[5]
 * @param[in] mS6 mass of internal field S[6]
 * @param[in] CpS1cS5cV4 coupling Cp[S[1], -S[5], -V[4]][Mom[S[1]] - Mom[-S[5]]
    ]
 * @param[in] CpS5S6V3 coupling Cp[S[5], S[6], V[3]][Mom[S[5]] - Mom[S[6]]]
 * @param[in] CpcS6V2V4 coupling Cp[-S[6], V[2], V[4]][g[lt2, lt3]]
 *
 * @return value of the one-loop diagram
 */
Decay_amplitude_SVV calculate_diagram_SVV_t1g6n6_VSS(
   double mext1, double mext2, double mext3,
   double mV4, double mS5, double mS6,
   const std::complex<double>& CpS1cS5cV4, const std::complex<double>& CpS5S6V3
      , const std::complex<double>& CpcS6V2V4,
   double scale)
{
   const auto c0tmp1 = passarino_veltman::C0(mext2*mext2, mext3*mext3, mext1*
      mext1, mV4*mV4, mS6*mS6, mS5*mS5, scale*scale);
   const auto c00tmp2 = passarino_veltman::C00(mext2*mext2, mext3*mext3, mext1*
      mext1, mV4*mV4, mS6*mS6, mS5*mS5, scale*scale);
   const auto c1tmp3 = passarino_veltman::C1(mext2*mext2, mext3*mext3, mext1*
      mext1, mV4*mV4, mS6*mS6, mS5*mS5, scale*scale);
   const auto c12tmp4 = passarino_veltman::C12(mext2*mext2, mext3*mext3, mext1*
      mext1, mV4*mV4, mS6*mS6, mS5*mS5, scale*scale);
   const auto c2tmp5 = passarino_veltman::C2(mext2*mext2, mext3*mext3, mext1*
      mext1, mV4*mV4, mS6*mS6, mS5*mS5, scale*scale);
   const auto c22tmp6 = passarino_veltman::C22(mext2*mext2, mext3*mext3, mext1*
      mext1, mV4*mV4, mS6*mS6, mS5*mS5, scale*scale);

   Decay_amplitude_SVV result;

   result.m_decay = mext1;
   result.m_out_1 = mext2;
   result.m_out_2 = mext3;

   result.form_factor_g = oneOver16PiSqr*(std::complex<double>(0,2)*c00tmp2*
      CpcS6V2V4*CpS1cS5cV4*CpS5S6V3);
   result.form_factor_11 = oneOver16PiSqr*(std::complex<double>(0,2)*(2*c0tmp1
      + c12tmp4 + 2*c1tmp3 + c22tmp6 + 3*c2tmp5)*CpcS6V2V4*CpS1cS5cV4*CpS5S6V3)
      ;
   result.form_factor_21 = oneOver16PiSqr*(std::complex<double>(0,2)*(2*c0tmp1
      + c12tmp4 + 2*c1tmp3 + c22tmp6 + 3*c2tmp5)*CpcS6V2V4*CpS1cS5cV4*CpS5S6V3)
      ;
   result.form_factor_12 = oneOver16PiSqr*(std::complex<double>(0,2)*(2*c0tmp1
      + c12tmp4 + 2*c1tmp3 + c22tmp6 + 3*c2tmp5)*CpcS6V2V4*CpS1cS5cV4*CpS5S6V3)
      ;
   result.form_factor_22 = oneOver16PiSqr*(std::complex<double>(0,2)*(2*c0tmp1
      + c12tmp4 + 2*c1tmp3 + c22tmp6 + 3*c2tmp5)*CpcS6V2V4*CpS1cS5cV4*CpS5S6V3)
      ;

   return result;
}

/**
 * @brief Evaluates T1G7N7 diagram for process S -> VV
 *
 * @param[in] mext1 mass of external field S[1]
 * @param[in] mext2 mass of external field V[2]
 * @param[in] mext3 mass of external field V[3]
 * @param[in] mS4 mass of internal field S[4]
 * @param[in] mV5 mass of internal field V[5]
 * @param[in] mV6 mass of internal field V[6]
 * @param[in] CpS1cS4cV5 coupling Cp[S[1], -S[4], -V[5]][Mom[S[1]] - Mom[-S[4]]
    ]
 * @param[in] CpS4V2cV6 coupling Cp[S[4], V[2], -V[6]][g[lt2, lt3]]
 * @param[in] CpV3V5V6 coupling Cp[V[3], V[5], V[6]][g[lt1, lt2] (-Mom[V[3]] +
    Mom[V[5]]) + g[lt1, lt3] (Mom[V[3]] - Mom[V[6]]) + g[lt2, lt3] (-Mom[V[5]]
    + Mom[V[6]])]
 *
 * @return value of the one-loop diagram
 */
Decay_amplitude_SVV calculate_diagram_SVV_t1g7n7_SVV(
   double mext1, double mext2, double mext3,
   double mS4, double mV5, double mV6,
   const std::complex<double>& CpS1cS4cV5, const std::complex<double>&
      CpS4V2cV6, const std::complex<double>& CpV3V5V6,
   double scale)
{
   const auto b0tmp1 = passarino_veltman::B0(mext3*mext3, mV5*mV5, mV6*mV6,
      scale*scale);
   const auto c0tmp2 = passarino_veltman::C0(mext2*mext2, mext3*mext3, mext1*
      mext1, mS4*mS4, mV6*mV6, mV5*mV5, scale*scale);
   const auto c00tmp3 = passarino_veltman::C00(mext2*mext2, mext3*mext3, mext1*
      mext1, mS4*mS4, mV6*mV6, mV5*mV5, scale*scale);
   const auto c1tmp4 = passarino_veltman::C1(mext2*mext2, mext3*mext3, mext1*
      mext1, mS4*mS4, mV6*mV6, mV5*mV5, scale*scale);
   const auto c12tmp5 = passarino_veltman::C12(mext2*mext2, mext3*mext3, mext1*
      mext1, mS4*mS4, mV6*mV6, mV5*mV5, scale*scale);
   const auto c2tmp6 = passarino_veltman::C2(mext2*mext2, mext3*mext3, mext1*
      mext1, mS4*mS4, mV6*mV6, mV5*mV5, scale*scale);
   const auto c22tmp7 = passarino_veltman::C22(mext2*mext2, mext3*mext3, mext1*
      mext1, mS4*mS4, mV6*mV6, mV5*mV5, scale*scale);

   Decay_amplitude_SVV result;

   result.m_decay = mext1;
   result.m_out_1 = mext2;
   result.m_out_2 = mext3;

   result.form_factor_g = oneOver16PiSqr*(std::complex<double>(0,1)*CpS1cS4cV5*
      CpS4V2cV6*CpV3V5V6*(b0tmp1 - c00tmp3 - c1tmp4*Sqr(mext1) - c2tmp6*Sqr(
      mext1) - c0tmp2*Sqr(mext2) + c1tmp4*Sqr(mext2) + c2tmp6*Sqr(mext2) +
      c0tmp2*Sqr(mext3) + c1tmp4*Sqr(mext3) - c2tmp6*Sqr(mext3) + c0tmp2*Sqr(
      mS4)));
   result.form_factor_11 = oneOver16PiSqr*(std::complex<double>(0,-1)*(c12tmp5
      - 4*c1tmp4 + c22tmp7 - c2tmp6)*CpS1cS4cV5*CpS4V2cV6*CpV3V5V6);
   result.form_factor_21 = oneOver16PiSqr*(std::complex<double>(0,-1)*(c12tmp5
      - 4*c1tmp4 + c22tmp7 - c2tmp6)*CpS1cS4cV5*CpS4V2cV6*CpV3V5V6);
   result.form_factor_12 = oneOver16PiSqr*(std::complex<double>(0,-1)*(c12tmp5
      - 4*c1tmp4 + c22tmp7 - c2tmp6)*CpS1cS4cV5*CpS4V2cV6*CpV3V5V6);
   result.form_factor_22 = oneOver16PiSqr*(std::complex<double>(0,-1)*(c12tmp5
      - 4*c1tmp4 + c22tmp7 - c2tmp6)*CpS1cS4cV5*CpS4V2cV6*CpV3V5V6);

   return result;
}

/**
 * @brief Evaluates T1G8N8 diagram for process S -> VV
 *
 * @param[in] mext1 mass of external field S[1]
 * @param[in] mext2 mass of external field V[2]
 * @param[in] mext3 mass of external field V[3]
 * @param[in] mV4 mass of internal field V[4]
 * @param[in] mS5 mass of internal field S[5]
 * @param[in] mV6 mass of internal field V[6]
 * @param[in] CpS1cS5cV4 coupling Cp[S[1], -S[5], -V[4]][Mom[S[1]] - Mom[-S[5]]
    ]
 * @param[in] CpS5V3V6 coupling Cp[S[5], V[3], V[6]][g[lt2, lt3]]
 * @param[in] CpV2V4cV6 coupling Cp[V[2], V[4], -V[6]][g[lt1, lt2] (-Mom[V[2]]
    + Mom[V[4]]) + g[lt1, lt3] (Mom[V[2]] - Mom[-V[6]]) + g[lt2, lt3] (-Mom[V[4
    ]] + Mom[-V[6]])]
 *
 * @return value of the one-loop diagram
 */
Decay_amplitude_SVV calculate_diagram_SVV_t1g8n8_VSV(
   double mext1, double mext2, double mext3,
   double mV4, double mS5, double mV6,
   const std::complex<double>& CpS1cS5cV4, const std::complex<double>& CpS5V3V6
      , const std::complex<double>& CpV2V4cV6,
   double scale)
{
   const auto b0tmp1 = passarino_veltman::B0(mext3*mext3, mS5*mS5, mV6*mV6,
      scale*scale);
   const auto c0tmp2 = passarino_veltman::C0(mext2*mext2, mext3*mext3, mext1*
      mext1, mV4*mV4, mV6*mV6, mS5*mS5, scale*scale);
   const auto c00tmp3 = passarino_veltman::C00(mext2*mext2, mext3*mext3, mext1*
      mext1, mV4*mV4, mV6*mV6, mS5*mS5, scale*scale);
   const auto c1tmp4 = passarino_veltman::C1(mext2*mext2, mext3*mext3, mext1*
      mext1, mV4*mV4, mV6*mV6, mS5*mS5, scale*scale);
   const auto c12tmp5 = passarino_veltman::C12(mext2*mext2, mext3*mext3, mext1*
      mext1, mV4*mV4, mV6*mV6, mS5*mS5, scale*scale);
   const auto c2tmp6 = passarino_veltman::C2(mext2*mext2, mext3*mext3, mext1*
      mext1, mV4*mV4, mV6*mV6, mS5*mS5, scale*scale);
   const auto c22tmp7 = passarino_veltman::C22(mext2*mext2, mext3*mext3, mext1*
      mext1, mV4*mV4, mV6*mV6, mS5*mS5, scale*scale);

   Decay_amplitude_SVV result;

   result.m_decay = mext1;
   result.m_out_1 = mext2;
   result.m_out_2 = mext3;

   result.form_factor_g = oneOver16PiSqr*(std::complex<double>(0,1)*CpS1cS5cV4*
      CpS5V3V6*CpV2V4cV6*(b0tmp1 - c00tmp3 + 2*c0tmp2*Sqr(mext1) + c1tmp4*Sqr(
      mext1) + 3*c2tmp6*Sqr(mext1) + 2*c0tmp2*Sqr(mext2) + 3*c1tmp4*Sqr(mext2)
      + c2tmp6*Sqr(mext2) - 2*c0tmp2*Sqr(mext3) - c1tmp4*Sqr(mext3) - c2tmp6*
      Sqr(mext3) + c0tmp2*Sqr(mV4)));
   result.form_factor_11 = oneOver16PiSqr*(std::complex<double>(0,-1)*(2*c0tmp2
       + c12tmp5 - 2*c1tmp4 + c22tmp7 + 3*c2tmp6)*CpS1cS5cV4*CpS5V3V6*CpV2V4cV6
      );
   result.form_factor_21 = oneOver16PiSqr*(std::complex<double>(0,-1)*(2*c0tmp2
       + c12tmp5 - 2*c1tmp4 + c22tmp7 + 3*c2tmp6)*CpS1cS5cV4*CpS5V3V6*CpV2V4cV6
      );
   result.form_factor_12 = oneOver16PiSqr*(std::complex<double>(0,-1)*(2*c0tmp2
       + c12tmp5 - 2*c1tmp4 + c22tmp7 + 3*c2tmp6)*CpS1cS5cV4*CpS5V3V6*CpV2V4cV6
      );
   result.form_factor_22 = oneOver16PiSqr*(std::complex<double>(0,-1)*(2*c0tmp2
       + c12tmp5 - 2*c1tmp4 + c22tmp7 + 3*c2tmp6)*CpS1cS5cV4*CpS5V3V6*CpV2V4cV6
      );

   return result;
}

/**
 * @brief Evaluates T1G9N9 diagram for process S -> VV
 *
 * @param[in] mext1 mass of external field S[1]
 * @param[in] mext2 mass of external field V[2]
 * @param[in] mext3 mass of external field V[3]
 * @param[in] mV4 mass of internal field V[4]
 * @param[in] mV5 mass of internal field V[5]
 * @param[in] mS6 mass of internal field S[6]
 * @param[in] CpS1cV4cV5 coupling Cp[S[1], -V[4], -V[5]][g[lt2, lt3]]
 * @param[in] CpcS6V2V4 coupling Cp[-S[6], V[2], V[4]][g[lt2, lt3]]
 * @param[in] CpS6V3V5 coupling Cp[S[6], V[3], V[5]][g[lt2, lt3]]
 *
 * @return value of the one-loop diagram
 */
Decay_amplitude_SVV calculate_diagram_SVV_t1g9n9_VVS(
   double mext1, double mext2, double mext3,
   double mV4, double mV5, double mS6,
   const std::complex<double>& CpS1cV4cV5, const std::complex<double>&
      CpcS6V2V4, const std::complex<double>& CpS6V3V5,
   double scale)
{
   const auto c0tmp1 = passarino_veltman::C0(mext2*mext2, mext3*mext3, mext1*
      mext1, mV4*mV4, mS6*mS6, mV5*mV5, scale*scale);

   Decay_amplitude_SVV result;

   result.m_decay = mext1;
   result.m_out_1 = mext2;
   result.m_out_2 = mext3;

   result.form_factor_g = oneOver16PiSqr*(std::complex<double>(0,-1)*c0tmp1*
      CpcS6V2V4*CpS1cV4cV5*CpS6V3V5);

   return result;
}

/**
 * @brief Evaluates T1G10N10 diagram for process S -> VV
 *
 * @param[in] mext1 mass of external field S[1]
 * @param[in] mext2 mass of external field V[2]
 * @param[in] mext3 mass of external field V[3]
 * @param[in] mV4 mass of internal field V[4]
 * @param[in] mV5 mass of internal field V[5]
 * @param[in] mV6 mass of internal field V[6]
 * @param[in] CpS1cV4cV5 coupling Cp[S[1], -V[4], -V[5]][g[lt2, lt3]]
 * @param[in] CpV2V4cV6 coupling Cp[V[2], V[4], -V[6]][g[lt1, lt2] (-Mom[V[2]]
    + Mom[V[4]]) + g[lt1, lt3] (Mom[V[2]] - Mom[-V[6]]) + g[lt2, lt3] (-Mom[V[4
    ]] + Mom[-V[6]])]
 * @param[in] CpV3V5V6 coupling Cp[V[3], V[5], V[6]][g[lt1, lt2] (-Mom[V[3]] +
    Mom[V[5]]) + g[lt1, lt3] (Mom[V[3]] - Mom[V[6]]) + g[lt2, lt3] (-Mom[V[5]]
    + Mom[V[6]])]
 *
 * @return value of the one-loop diagram
 */
Decay_amplitude_SVV calculate_diagram_SVV_t1g10n10_VVV(
   double mext1, double mext2, double mext3,
   double mV4, double mV5, double mV6,
   const std::complex<double>& CpS1cV4cV5, const std::complex<double>&
      CpV2V4cV6, const std::complex<double>& CpV3V5V6,
   double scale, double finite)
{
   const auto b0tmp1 = passarino_veltman::B0(mext3*mext3, mV5*mV5, mV6*mV6,
      scale*scale);
   const auto c0tmp2 = passarino_veltman::C0(mext2*mext2, mext3*mext3, mext1*
      mext1, mV4*mV4, mV6*mV6, mV5*mV5, scale*scale);
   const auto c00tmp3 = passarino_veltman::C00(mext2*mext2, mext3*mext3, mext1*
      mext1, mV4*mV4, mV6*mV6, mV5*mV5, scale*scale);
   const auto c1tmp4 = passarino_veltman::C1(mext2*mext2, mext3*mext3, mext1*
      mext1, mV4*mV4, mV6*mV6, mV5*mV5, scale*scale);
   const auto c12tmp5 = passarino_veltman::C12(mext2*mext2, mext3*mext3, mext1*
      mext1, mV4*mV4, mV6*mV6, mV5*mV5, scale*scale);
   const auto c2tmp6 = passarino_veltman::C2(mext2*mext2, mext3*mext3, mext1*
      mext1, mV4*mV4, mV6*mV6, mV5*mV5, scale*scale);
   const auto c22tmp7 = passarino_veltman::C22(mext2*mext2, mext3*mext3, mext1*
      mext1, mV4*mV4, mV6*mV6, mV5*mV5, scale*scale);

   Decay_amplitude_SVV result;

   result.m_decay = mext1;
   result.m_out_1 = mext2;
   result.m_out_2 = mext3;

   result.form_factor_g = oneOver16PiSqr*(std::complex<double>(0,-0.5)*
      CpS1cV4cV5*CpV2V4cV6*CpV3V5V6*(4*b0tmp1 + 20*c00tmp3 - 4*finite - 4*
      c0tmp2*Sqr(mext1) + c1tmp4*Sqr(mext1) + 4*c2tmp6*Sqr(mext1) + 6*c0tmp2*
      Sqr(mext2) + 5*c1tmp4*Sqr(mext2) + 2*c2tmp6*Sqr(mext2) + 4*c0tmp2*Sqr(
      mext3) - c1tmp4*Sqr(mext3) - 2*c2tmp6*Sqr(mext3) + 4*c0tmp2*Sqr(mV4)));
   result.form_factor_11 = oneOver16PiSqr*(std::complex<double>(0,-1)*(5*c0tmp2
       + c1tmp4 + 10*(c12tmp5 + c22tmp7 + c2tmp6))*CpS1cV4cV5*CpV2V4cV6*
      CpV3V5V6);
   result.form_factor_21 = oneOver16PiSqr*(std::complex<double>(0,-1)*(5*c0tmp2
       + c1tmp4 + 10*(c12tmp5 + c22tmp7 + c2tmp6))*CpS1cV4cV5*CpV2V4cV6*
      CpV3V5V6);
   result.form_factor_12 = oneOver16PiSqr*(std::complex<double>(0,-1)*(5*c0tmp2
       + c1tmp4 + 10*(c12tmp5 + c22tmp7 + c2tmp6))*CpS1cV4cV5*CpV2V4cV6*
      CpV3V5V6);
   result.form_factor_22 = oneOver16PiSqr*(std::complex<double>(0,-1)*(5*c0tmp2
       + c1tmp4 + 10*(c12tmp5 + c22tmp7 + c2tmp6))*CpS1cV4cV5*CpV2V4cV6*
      CpV3V5V6);

   return result;
}

/**
 * @brief Evaluates T2G1N11 diagram for process S -> VV
 *
 * @param[in] mext1 mass of external field S[1]
 * @param[in] mext2 mass of external field V[2]
 * @param[in] mext3 mass of external field V[3]
 * @param[in] mS4 mass of internal field S[4]
 * @param[in] mS5 mass of internal field S[5]
 * @param[in] CpS1cS4cS5 coupling Cp[S[1], -S[4], -S[5]]
 * @param[in] CpS4S5V2V3 coupling Cp[S[4], S[5], V[2], V[3]][g[lt3, lt4]]
 *
 * @return value of the one-loop diagram
 */
Decay_amplitude_SVV calculate_diagram_SVV_t2g1n11_SS(
   double mext1, double mext2, double mext3,
   double mS4, double mS5,
   const std::complex<double>& CpS1cS4cS5, const std::complex<double>&
      CpS4S5V2V3,
   double scale)
{
   const auto b0tmp1 = passarino_veltman::B0(mext1*mext1, mS4*mS4, mS5*mS5,
      scale*scale);

   Decay_amplitude_SVV result;

   result.m_decay = mext1;
   result.m_out_1 = mext2;
   result.m_out_2 = mext3;

   result.form_factor_g = oneOver16PiSqr*(-0.5*b0tmp1*CpS1cS4cS5*CpS4S5V2V3);

   return result;
}

/**
 * @brief Evaluates T2G2N12 diagram for process S -> VV
 *
 * @param[in] mext1 mass of external field S[1]
 * @param[in] mext2 mass of external field V[2]
 * @param[in] mext3 mass of external field V[3]
 * @param[in] mV4 mass of internal field V[4]
 * @param[in] mV5 mass of internal field V[5]
 * @param[in] CpS1cV4cV5 coupling Cp[S[1], -V[4], -V[5]][g[lt2, lt3]]
 * @param[in] CpV2V3V4V51 coupling Cp[V[2], V[3], V[4], V[5]][g[lt1, lt2] g[lt3
    , lt4]]
 * @param[in] CpV2V3V4V52 coupling Cp[V[2], V[3], V[4], V[5]][g[lt1, lt4] g[lt2
    , lt3]]
 * @param[in] CpV2V3V4V53 coupling Cp[V[2], V[3], V[4], V[5]][g[lt1, lt3] g[lt2
    , lt4]]
 *
 * @return value of the one-loop diagram
 */
Decay_amplitude_SVV calculate_diagram_SVV_t2g2n12_VV(
   double mext1, double mext2, double mext3,
   double mV4, double mV5,
   const std::complex<double>& CpS1cV4cV5, const std::complex<double>&
      CpV2V3V4V51, const std::complex<double>& CpV2V3V4V52, const std::complex<
      double>& CpV2V3V4V53,
   double scale, double finite)
{
   const auto b0tmp1 = passarino_veltman::B0(mext1*mext1, mV4*mV4, mV5*mV5,
      scale*scale);

   Decay_amplitude_SVV result;

   result.m_decay = mext1;
   result.m_out_1 = mext2;
   result.m_out_2 = mext3;

   result.form_factor_g = oneOver16PiSqr*(-0.5*b0tmp1*CpS1cV4cV5*(4*CpV2V3V4V51
       + CpV2V3V4V52 + CpV2V3V4V53) + CpS1cV4cV5*CpV2V3V4V51*finite);

   return result;
}

/**
 * @brief Evaluates T3G1N13 diagram for process S -> VV
 *
 * @param[in] mext1 mass of external field S[1]
 * @param[in] mext2 mass of external field V[2]
 * @param[in] mext3 mass of external field V[3]
 * @param[in] mS4 mass of internal field S[4]
 * @param[in] mV5 mass of internal field V[5]
 * @param[in] CpcS4V2cV5 coupling Cp[-S[4], V[2], -V[5]][g[lt2, lt3]]
 * @param[in] CpS1S4V3V5 coupling Cp[S[1], S[4], V[3], V[5]][g[lt3, lt4]]
 *
 * @return value of the one-loop diagram
 */
Decay_amplitude_SVV calculate_diagram_SVV_t3g1n13_SV(
   double mext1, double mext2, double mext3,
   double mS4, double mV5,
   const std::complex<double>& CpcS4V2cV5, const std::complex<double>&
      CpS1S4V3V5,
   double scale)
{
   const auto b0tmp1 = passarino_veltman::B0(mext2*mext2, mS4*mS4, mV5*mV5,
      scale*scale);

   Decay_amplitude_SVV result;

   result.m_decay = mext1;
   result.m_out_1 = mext2;
   result.m_out_2 = mext3;

   result.form_factor_g = oneOver16PiSqr*(b0tmp1*CpcS4V2cV5*CpS1S4V3V5);

   return result;
}

/**
 * @brief Evaluates T4G1N14 diagram for process S -> VV
 *
 * @param[in] mext1 mass of external field S[1]
 * @param[in] mext2 mass of external field V[2]
 * @param[in] mext3 mass of external field V[3]
 * @param[in] mS4 mass of internal field S[4]
 * @param[in] mV5 mass of internal field V[5]
 * @param[in] CpcS4V3cV5 coupling Cp[-S[4], V[3], -V[5]][g[lt2, lt3]]
 * @param[in] CpS1S4V2V5 coupling Cp[S[1], S[4], V[2], V[5]][g[lt3, lt4]]
 *
 * @return value of the one-loop diagram
 */
Decay_amplitude_SVV calculate_diagram_SVV_t4g1n14_SV(
   double mext1, double mext2, double mext3,
   double mS4, double mV5,
   const std::complex<double>& CpcS4V3cV5, const std::complex<double>&
      CpS1S4V2V5,
   double scale)
{
   const auto b0tmp1 = passarino_veltman::B0(mext3*mext3, mS4*mS4, mV5*mV5,
      scale*scale);

   Decay_amplitude_SVV result;

   result.m_decay = mext1;
   result.m_out_1 = mext2;
   result.m_out_2 = mext3;

   result.form_factor_g = oneOver16PiSqr*(b0tmp1*CpcS4V3cV5*CpS1S4V2V5);

   return result;
}

/**
 * @brief Evaluates T1G1N1 diagram for process S -> SV
 *
 * @param[in] mext1 mass of external field S[1]
 * @param[in] mext2 mass of external field S[2]
 * @param[in] mext3 mass of external field V[3]
 * @param[in] mF4 mass of internal field F[4]
 * @param[in] mF5 mass of internal field F[5]
 * @param[in] mF6 mass of internal field F[6]
 * @param[in] CpcF4cF5S1PL coupling Cp[-F[4], -F[5], S[1]][PL]
 * @param[in] CpcF4cF5S1PR coupling Cp[-F[4], -F[5], S[1]][PR]
 * @param[in] CpF5F6V3PL coupling Cp[F[5], F[6], V[3]][LorentzProduct[gamma[lt3
    ], PL]]
 * @param[in] CpF5F6V3PR coupling Cp[F[5], F[6], V[3]][LorentzProduct[gamma[lt3
    ], PR]]
 * @param[in] CpcF6F4S2PL coupling Cp[-F[6], F[4], S[2]][PL]
 * @param[in] CpcF6F4S2PR coupling Cp[-F[6], F[4], S[2]][PR]
 *
 * @return value of the one-loop diagram
 */
Decay_amplitude_SSV calculate_diagram_SSV_t1g1n1_FFF(
   double mext1, double mext2, double mext3,
   double mF4, double mF5, double mF6,
   const std::complex<double>& CpcF4cF5S1PL, const std::complex<double>&
      CpcF4cF5S1PR, const std::complex<double>& CpF5F6V3PL, const std::complex<
      double>& CpF5F6V3PR, const std::complex<double>& CpcF6F4S2PL, const std::
      complex<double>& CpcF6F4S2PR,
   double scale)
{
   const auto b0tmp1 = passarino_veltman::B0(mext3*mext3, mF5*mF5, mF6*mF6,
      scale*scale);
   const auto c0tmp2 = passarino_veltman::C0(mext2*mext2, mext3*mext3, mext1*
      mext1, mF4*mF4, mF6*mF6, mF5*mF5, scale*scale);
   const auto c1tmp3 = passarino_veltman::C1(mext2*mext2, mext3*mext3, mext1*
      mext1, mF4*mF4, mF6*mF6, mF5*mF5, scale*scale);
   const auto c2tmp4 = passarino_veltman::C2(mext2*mext2, mext3*mext3, mext1*
      mext1, mF4*mF4, mF6*mF6, mF5*mF5, scale*scale);

   Decay_amplitude_SSV result;

   result.m_decay = mext1;
   result.m_scalar = mext2;
   result.m_vector = mext3;

   result.form_factor = oneOver16PiSqr*(std::complex<double>(0,-2)*(b0tmp1*
      CpcF4cF5S1PL*CpcF6F4S2PR*CpF5F6V3PL + b0tmp1*CpcF4cF5S1PR*CpcF6F4S2PL*
      CpF5F6V3PR + c2tmp4*CpcF4cF5S1PR*CpcF6F4S2PR*CpF5F6V3PL*mF4*mF5 + c2tmp4*
      CpcF4cF5S1PL*CpcF6F4S2PL*CpF5F6V3PR*mF4*mF5 + c2tmp4*CpcF4cF5S1PL*
      CpcF6F4S2PL*CpF5F6V3PL*mF4*mF6 + c2tmp4*CpcF4cF5S1PR*CpcF6F4S2PR*
      CpF5F6V3PR*mF4*mF6 + c2tmp4*CpcF4cF5S1PR*CpcF6F4S2PL*CpF5F6V3PL*mF5*mF6 +
      c2tmp4*CpcF4cF5S1PL*CpcF6F4S2PR*CpF5F6V3PR*mF5*mF6 + c0tmp2*mF4*(
      CpcF4cF5S1PL*(2*CpcF6F4S2PR*CpF5F6V3PL*mF4 + CpcF6F4S2PL*CpF5F6V3PR*mF5 +
      CpcF6F4S2PL*CpF5F6V3PL*mF6) + CpcF4cF5S1PR*(2*CpcF6F4S2PL*CpF5F6V3PR*mF4
      + CpcF6F4S2PR*CpF5F6V3PL*mF5 + CpcF6F4S2PR*CpF5F6V3PR*mF6)) + c2tmp4*
      CpcF4cF5S1PL*CpcF6F4S2PR*CpF5F6V3PL*Sqr(mext1) + c2tmp4*CpcF4cF5S1PR*
      CpcF6F4S2PL*CpF5F6V3PR*Sqr(mext1) + c2tmp4*CpcF4cF5S1PL*CpcF6F4S2PR*
      CpF5F6V3PL*Sqr(mF4) + c2tmp4*CpcF4cF5S1PR*CpcF6F4S2PL*CpF5F6V3PR*Sqr(mF4)
      + c1tmp3*(CpcF4cF5S1PL*CpcF6F4S2PL*mF4*(CpF5F6V3PR*mF5 + CpF5F6V3PL*mF6)
      + CpcF4cF5S1PR*CpcF6F4S2PR*mF4*(CpF5F6V3PL*mF5 + CpF5F6V3PR*mF6) +
      CpcF4cF5S1PL*CpcF6F4S2PR*(CpF5F6V3PR*mF5*mF6 + CpF5F6V3PL*(Sqr(mext2) +
      Sqr(mF4))) + CpcF4cF5S1PR*CpcF6F4S2PL*(CpF5F6V3PL*mF5*mF6 + CpF5F6V3PR*(
      Sqr(mext2) + Sqr(mF4))))));

   return result;
}

/**
 * @brief Evaluates T1G2N2 diagram for process S -> SV
 *
 * @param[in] mext1 mass of external field S[1]
 * @param[in] mext2 mass of external field S[2]
 * @param[in] mext3 mass of external field V[3]
 * @param[in] mS4 mass of internal field S[4]
 * @param[in] mS5 mass of internal field S[5]
 * @param[in] mS6 mass of internal field S[6]
 * @param[in] CpS1cS4cS5 coupling Cp[S[1], -S[4], -S[5]]
 * @param[in] CpS2S4cS6 coupling Cp[S[2], S[4], -S[6]]
 * @param[in] CpS5S6V3 coupling Cp[S[5], S[6], V[3]][Mom[S[5]] - Mom[S[6]]]
 *
 * @return value of the one-loop diagram
 */
Decay_amplitude_SSV calculate_diagram_SSV_t1g2n2_SSS(
   double mext1, double mext2, double mext3,
   double mS4, double mS5, double mS6,
   const std::complex<double>& CpS1cS4cS5, const std::complex<double>&
      CpS2S4cS6, const std::complex<double>& CpS5S6V3,
   double scale)
{
   const auto c0tmp1 = passarino_veltman::C0(mext2*mext2, mext3*mext3, mext1*
      mext1, mS4*mS4, mS6*mS6, mS5*mS5, scale*scale);
   const auto c1tmp2 = passarino_veltman::C1(mext2*mext2, mext3*mext3, mext1*
      mext1, mS4*mS4, mS6*mS6, mS5*mS5, scale*scale);
   const auto c2tmp3 = passarino_veltman::C2(mext2*mext2, mext3*mext3, mext1*
      mext1, mS4*mS4, mS6*mS6, mS5*mS5, scale*scale);

   Decay_amplitude_SSV result;

   result.m_decay = mext1;
   result.m_scalar = mext2;
   result.m_vector = mext3;

   result.form_factor = oneOver16PiSqr*(std::complex<double>(0,-2)*(c0tmp1 +
      c1tmp2 + c2tmp3)*CpS1cS4cS5*CpS2S4cS6*CpS5S6V3);

   return result;
}

/**
 * @brief Evaluates T1G3N3 diagram for process S -> SV
 *
 * @param[in] mext1 mass of external field S[1]
 * @param[in] mext2 mass of external field S[2]
 * @param[in] mext3 mass of external field V[3]
 * @param[in] mU4 mass of internal field U[4]
 * @param[in] mU5 mass of internal field U[5]
 * @param[in] mU6 mass of internal field U[6]
 * @param[in] CpS1cU4cU5 coupling Cp[S[1], -U[4], -U[5]]
 * @param[in] CpS2cU6U4 coupling Cp[S[2], -U[6], U[4]]
 * @param[in] CpU5U6V3 coupling Cp[U[5], U[6], V[3]][Mom[U[5]]]
 *
 * @return value of the one-loop diagram
 */
Decay_amplitude_SSV calculate_diagram_SSV_t1g3n3_UUU(
   double mext1, double mext2, double mext3,
   double mU4, double mU5, double mU6,
   const std::complex<double>& CpS1cU4cU5, const std::complex<double>&
      CpS2cU6U4, const std::complex<double>& CpU5U6V3,
   double scale)
{
   const auto c0tmp1 = passarino_veltman::C0(mext2*mext2, mext3*mext3, mext1*
      mext1, mU4*mU4, mU6*mU6, mU5*mU5, scale*scale);
   const auto c1tmp2 = passarino_veltman::C1(mext2*mext2, mext3*mext3, mext1*
      mext1, mU4*mU4, mU6*mU6, mU5*mU5, scale*scale);
   const auto c2tmp3 = passarino_veltman::C2(mext2*mext2, mext3*mext3, mext1*
      mext1, mU4*mU4, mU6*mU6, mU5*mU5, scale*scale);

   Decay_amplitude_SSV result;

   result.m_decay = mext1;
   result.m_scalar = mext2;
   result.m_vector = mext3;

   result.form_factor = oneOver16PiSqr*(std::complex<double>(0,1)*(c0tmp1 +
      c1tmp2 + c2tmp3)*CpS1cU4cU5*CpS2cU6U4*CpU5U6V3);

   return result;
}

/**
 * @brief Evaluates T1G4N4 diagram for process S -> SV
 *
 * @param[in] mext1 mass of external field S[1]
 * @param[in] mext2 mass of external field S[2]
 * @param[in] mext3 mass of external field V[3]
 * @param[in] mS4 mass of internal field S[4]
 * @param[in] mS5 mass of internal field S[5]
 * @param[in] mV6 mass of internal field V[6]
 * @param[in] CpS1cS4cS5 coupling Cp[S[1], -S[4], -S[5]]
 * @param[in] CpS2S4cV6 coupling Cp[S[2], S[4], -V[6]][Mom[S[2]] - Mom[S[4]]]
 * @param[in] CpS5V3V6 coupling Cp[S[5], V[3], V[6]][g[lt2, lt3]]
 *
 * @return value of the one-loop diagram
 */
Decay_amplitude_SSV calculate_diagram_SSV_t1g4n4_SSV(
   double mext1, double mext2, double mext3,
   double mS4, double mS5, double mV6,
   const std::complex<double>& CpS1cS4cS5, const std::complex<double>&
      CpS2S4cV6, const std::complex<double>& CpS5V3V6,
   double scale)
{
   const auto c0tmp1 = passarino_veltman::C0(mext2*mext2, mext3*mext3, mext1*
      mext1, mS4*mS4, mV6*mV6, mS5*mS5, scale*scale);
   const auto c1tmp2 = passarino_veltman::C1(mext2*mext2, mext3*mext3, mext1*
      mext1, mS4*mS4, mV6*mV6, mS5*mS5, scale*scale);
   const auto c2tmp3 = passarino_veltman::C2(mext2*mext2, mext3*mext3, mext1*
      mext1, mS4*mS4, mV6*mV6, mS5*mS5, scale*scale);

   Decay_amplitude_SSV result;

   result.m_decay = mext1;
   result.m_scalar = mext2;
   result.m_vector = mext3;

   result.form_factor = oneOver16PiSqr*(std::complex<double>(0,-1)*(c0tmp1 -
      c1tmp2 - c2tmp3)*CpS1cS4cS5*CpS2S4cV6*CpS5V3V6);

   return result;
}

/**
 * @brief Evaluates T1G5N5 diagram for process S -> SV
 *
 * @param[in] mext1 mass of external field S[1]
 * @param[in] mext2 mass of external field S[2]
 * @param[in] mext3 mass of external field V[3]
 * @param[in] mS4 mass of internal field S[4]
 * @param[in] mV5 mass of internal field V[5]
 * @param[in] mS6 mass of internal field S[6]
 * @param[in] CpS1cS4cV5 coupling Cp[S[1], -S[4], -V[5]][Mom[S[1]] - Mom[-S[4]]
    ]
 * @param[in] CpS2S4cS6 coupling Cp[S[2], S[4], -S[6]]
 * @param[in] CpS6V3V5 coupling Cp[S[6], V[3], V[5]][g[lt2, lt3]]
 *
 * @return value of the one-loop diagram
 */
Decay_amplitude_SSV calculate_diagram_SSV_t1g5n5_SVS(
   double mext1, double mext2, double mext3,
   double mS4, double mV5, double mS6,
   const std::complex<double>& CpS1cS4cV5, const std::complex<double>&
      CpS2S4cS6, const std::complex<double>& CpS6V3V5,
   double scale)
{
   const auto c0tmp1 = passarino_veltman::C0(mext2*mext2, mext3*mext3, mext1*
      mext1, mS4*mS4, mS6*mS6, mV5*mV5, scale*scale);
   const auto c1tmp2 = passarino_veltman::C1(mext2*mext2, mext3*mext3, mext1*
      mext1, mS4*mS4, mS6*mS6, mV5*mV5, scale*scale);
   const auto c2tmp3 = passarino_veltman::C2(mext2*mext2, mext3*mext3, mext1*
      mext1, mS4*mS4, mS6*mS6, mV5*mV5, scale*scale);

   Decay_amplitude_SSV result;

   result.m_decay = mext1;
   result.m_scalar = mext2;
   result.m_vector = mext3;

   result.form_factor = oneOver16PiSqr*(std::complex<double>(0,1)*(c0tmp1 -
      c1tmp2 - c2tmp3)*CpS1cS4cV5*CpS2S4cS6*CpS6V3V5);

   return result;
}

/**
 * @brief Evaluates T1G6N6 diagram for process S -> SV
 *
 * @param[in] mext1 mass of external field S[1]
 * @param[in] mext2 mass of external field S[2]
 * @param[in] mext3 mass of external field V[3]
 * @param[in] mV4 mass of internal field V[4]
 * @param[in] mS5 mass of internal field S[5]
 * @param[in] mS6 mass of internal field S[6]
 * @param[in] CpS1cS5cV4 coupling Cp[S[1], -S[5], -V[4]][Mom[S[1]] - Mom[-S[5]]
    ]
 * @param[in] CpS2cS6V4 coupling Cp[S[2], -S[6], V[4]][Mom[S[2]] - Mom[-S[6]]]
 * @param[in] CpS5S6V3 coupling Cp[S[5], S[6], V[3]][Mom[S[5]] - Mom[S[6]]]
 *
 * @return value of the one-loop diagram
 */
Decay_amplitude_SSV calculate_diagram_SSV_t1g6n6_VSS(
   double mext1, double mext2, double mext3,
   double mV4, double mS5, double mS6,
   const std::complex<double>& CpS1cS5cV4, const std::complex<double>&
      CpS2cS6V4, const std::complex<double>& CpS5S6V3,
   double scale)
{
   const auto c0tmp1 = passarino_veltman::C0(mext2*mext2, mext3*mext3, mext1*
      mext1, mV4*mV4, mS6*mS6, mS5*mS5, scale*scale);
   const auto c00tmp2 = passarino_veltman::C00(mext2*mext2, mext3*mext3, mext1*
      mext1, mV4*mV4, mS6*mS6, mS5*mS5, scale*scale);
   const auto c1tmp3 = passarino_veltman::C1(mext2*mext2, mext3*mext3, mext1*
      mext1, mV4*mV4, mS6*mS6, mS5*mS5, scale*scale);
   const auto c11tmp4 = passarino_veltman::C11(mext2*mext2, mext3*mext3, mext1*
      mext1, mV4*mV4, mS6*mS6, mS5*mS5, scale*scale);
   const auto c12tmp5 = passarino_veltman::C12(mext2*mext2, mext3*mext3, mext1*
      mext1, mV4*mV4, mS6*mS6, mS5*mS5, scale*scale);
   const auto c2tmp6 = passarino_veltman::C2(mext2*mext2, mext3*mext3, mext1*
      mext1, mV4*mV4, mS6*mS6, mS5*mS5, scale*scale);
   const auto c22tmp7 = passarino_veltman::C22(mext2*mext2, mext3*mext3, mext1*
      mext1, mV4*mV4, mS6*mS6, mS5*mS5, scale*scale);

   Decay_amplitude_SSV result;

   result.m_decay = mext1;
   result.m_scalar = mext2;
   result.m_vector = mext3;

   result.form_factor = oneOver16PiSqr*(std::complex<double>(0,-2)*CpS1cS5cV4*
      CpS2cS6V4*CpS5S6V3*(4*c00tmp2 + c11tmp4*Sqr(mext1) + 4*c12tmp5*Sqr(mext1)
      + 3*c1tmp3*Sqr(mext1) + 3*c22tmp7*Sqr(mext1) + 5*c2tmp6*Sqr(mext1) + 3*
      c11tmp4*Sqr(mext2) + 4*c12tmp5*Sqr(mext2) + 5*c1tmp3*Sqr(mext2) + c22tmp7
      *Sqr(mext2) + 3*c2tmp6*Sqr(mext2) - c11tmp4*Sqr(mext3) - 2*c12tmp5*Sqr(
      mext3) - 3*c1tmp3*Sqr(mext3) - c22tmp7*Sqr(mext3) - 3*c2tmp6*Sqr(mext3) +
      c1tmp3*Sqr(mV4) + c2tmp6*Sqr(mV4) + c0tmp1*(2*(Sqr(mext1) + Sqr(mext2) -
      Sqr(mext3)) + Sqr(mV4))));

   return result;
}

/**
 * @brief Evaluates T1G7N7 diagram for process S -> SV
 *
 * @param[in] mext1 mass of external field S[1]
 * @param[in] mext2 mass of external field S[2]
 * @param[in] mext3 mass of external field V[3]
 * @param[in] mS4 mass of internal field S[4]
 * @param[in] mV5 mass of internal field V[5]
 * @param[in] mV6 mass of internal field V[6]
 * @param[in] CpS1cS4cV5 coupling Cp[S[1], -S[4], -V[5]][Mom[S[1]] - Mom[-S[4]]
    ]
 * @param[in] CpS2S4cV6 coupling Cp[S[2], S[4], -V[6]][Mom[S[2]] - Mom[S[4]]]
 * @param[in] CpV3V5V6 coupling Cp[V[3], V[5], V[6]][g[lt1, lt2] (-Mom[V[3]] +
    Mom[V[5]]) + g[lt1, lt3] (Mom[V[3]] - Mom[V[6]]) + g[lt2, lt3] (-Mom[V[5]]
    + Mom[V[6]])]
 *
 * @return value of the one-loop diagram
 */
Decay_amplitude_SSV calculate_diagram_SSV_t1g7n7_SVV(
   double mext1, double mext2, double mext3,
   double mS4, double mV5, double mV6,
   const std::complex<double>& CpS1cS4cV5, const std::complex<double>&
      CpS2S4cV6, const std::complex<double>& CpV3V5V6,
   double scale)
{
   const auto b0tmp1 = passarino_veltman::B0(mext3*mext3, mV5*mV5, mV6*mV6,
      scale*scale);
   const auto c0tmp2 = passarino_veltman::C0(mext2*mext2, mext3*mext3, mext1*
      mext1, mS4*mS4, mV6*mV6, mV5*mV5, scale*scale);
   const auto c00tmp3 = passarino_veltman::C00(mext2*mext2, mext3*mext3, mext1*
      mext1, mS4*mS4, mV6*mV6, mV5*mV5, scale*scale);
   const auto c1tmp4 = passarino_veltman::C1(mext2*mext2, mext3*mext3, mext1*
      mext1, mS4*mS4, mV6*mV6, mV5*mV5, scale*scale);
   const auto c11tmp5 = passarino_veltman::C11(mext2*mext2, mext3*mext3, mext1*
      mext1, mS4*mS4, mV6*mV6, mV5*mV5, scale*scale);
   const auto c12tmp6 = passarino_veltman::C12(mext2*mext2, mext3*mext3, mext1*
      mext1, mS4*mS4, mV6*mV6, mV5*mV5, scale*scale);
   const auto c2tmp7 = passarino_veltman::C2(mext2*mext2, mext3*mext3, mext1*
      mext1, mS4*mS4, mV6*mV6, mV5*mV5, scale*scale);
   const auto c22tmp8 = passarino_veltman::C22(mext2*mext2, mext3*mext3, mext1*
      mext1, mS4*mS4, mV6*mV6, mV5*mV5, scale*scale);

   Decay_amplitude_SSV result;

   result.m_decay = mext1;
   result.m_scalar = mext2;
   result.m_vector = mext3;

   result.form_factor = oneOver16PiSqr*(std::complex<double>(0,-1)*CpS1cS4cV5*
      CpS2S4cV6*CpV3V5V6*(4*b0tmp1 - 4*c00tmp3 - c11tmp5*Sqr(mext1) - 4*c12tmp6
      *Sqr(mext1) + c1tmp4*Sqr(mext1) - 3*c22tmp8*Sqr(mext1) - c2tmp7*Sqr(mext1
      ) - 3*c11tmp5*Sqr(mext2) - 4*c12tmp6*Sqr(mext2) - c1tmp4*Sqr(mext2) -
      c22tmp8*Sqr(mext2) + c2tmp7*Sqr(mext2) + c0tmp2*Sqr(mext3) + c11tmp5*Sqr(
      mext3) + 2*c12tmp6*Sqr(mext3) - 2*c1tmp4*Sqr(mext3) + c22tmp8*Sqr(mext3)
      - 2*c2tmp7*Sqr(mext3) + 4*c0tmp2*Sqr(mS4)));

   return result;
}

/**
 * @brief Evaluates T1G8N8 diagram for process S -> SV
 *
 * @param[in] mext1 mass of external field S[1]
 * @param[in] mext2 mass of external field S[2]
 * @param[in] mext3 mass of external field V[3]
 * @param[in] mV4 mass of internal field V[4]
 * @param[in] mS5 mass of internal field S[5]
 * @param[in] mV6 mass of internal field V[6]
 * @param[in] CpS1cS5cV4 coupling Cp[S[1], -S[5], -V[4]][Mom[S[1]] - Mom[-S[5]]
    ]
 * @param[in] CpS2V4cV6 coupling Cp[S[2], V[4], -V[6]][g[lt2, lt3]]
 * @param[in] CpS5V3V6 coupling Cp[S[5], V[3], V[6]][g[lt2, lt3]]
 *
 * @return value of the one-loop diagram
 */
Decay_amplitude_SSV calculate_diagram_SSV_t1g8n8_VSV(
   double mext1, double mext2, double mext3,
   double mV4, double mS5, double mV6,
   const std::complex<double>& CpS1cS5cV4, const std::complex<double>&
      CpS2V4cV6, const std::complex<double>& CpS5V3V6,
   double scale)
{
   const auto c0tmp1 = passarino_veltman::C0(mext2*mext2, mext3*mext3, mext1*
      mext1, mV4*mV4, mV6*mV6, mS5*mS5, scale*scale);
   const auto c1tmp2 = passarino_veltman::C1(mext2*mext2, mext3*mext3, mext1*
      mext1, mV4*mV4, mV6*mV6, mS5*mS5, scale*scale);
   const auto c2tmp3 = passarino_veltman::C2(mext2*mext2, mext3*mext3, mext1*
      mext1, mV4*mV4, mV6*mV6, mS5*mS5, scale*scale);

   Decay_amplitude_SSV result;

   result.m_decay = mext1;
   result.m_scalar = mext2;
   result.m_vector = mext3;

   result.form_factor = oneOver16PiSqr*(std::complex<double>(0,-1)*(2*c0tmp1 +
      c1tmp2 + c2tmp3)*CpS1cS5cV4*CpS2V4cV6*CpS5V3V6);

   return result;
}

/**
 * @brief Evaluates T1G9N9 diagram for process S -> SV
 *
 * @param[in] mext1 mass of external field S[1]
 * @param[in] mext2 mass of external field S[2]
 * @param[in] mext3 mass of external field V[3]
 * @param[in] mV4 mass of internal field V[4]
 * @param[in] mV5 mass of internal field V[5]
 * @param[in] mS6 mass of internal field S[6]
 * @param[in] CpS1cV4cV5 coupling Cp[S[1], -V[4], -V[5]][g[lt2, lt3]]
 * @param[in] CpS2cS6V4 coupling Cp[S[2], -S[6], V[4]][Mom[S[2]] - Mom[-S[6]]]
 * @param[in] CpS6V3V5 coupling Cp[S[6], V[3], V[5]][g[lt2, lt3]]
 *
 * @return value of the one-loop diagram
 */
Decay_amplitude_SSV calculate_diagram_SSV_t1g9n9_VVS(
   double mext1, double mext2, double mext3,
   double mV4, double mV5, double mS6,
   const std::complex<double>& CpS1cV4cV5, const std::complex<double>&
      CpS2cS6V4, const std::complex<double>& CpS6V3V5,
   double scale)
{
   const auto c0tmp1 = passarino_veltman::C0(mext2*mext2, mext3*mext3, mext1*
      mext1, mV4*mV4, mS6*mS6, mV5*mV5, scale*scale);
   const auto c1tmp2 = passarino_veltman::C1(mext2*mext2, mext3*mext3, mext1*
      mext1, mV4*mV4, mS6*mS6, mV5*mV5, scale*scale);
   const auto c2tmp3 = passarino_veltman::C2(mext2*mext2, mext3*mext3, mext1*
      mext1, mV4*mV4, mS6*mS6, mV5*mV5, scale*scale);

   Decay_amplitude_SSV result;

   result.m_decay = mext1;
   result.m_scalar = mext2;
   result.m_vector = mext3;

   result.form_factor = oneOver16PiSqr*(std::complex<double>(0,1)*(2*c0tmp1 +
      c1tmp2 + c2tmp3)*CpS1cV4cV5*CpS2cS6V4*CpS6V3V5);

   return result;
}

/**
 * @brief Evaluates T1G10N10 diagram for process S -> SV
 *
 * @param[in] mext1 mass of external field S[1]
 * @param[in] mext2 mass of external field S[2]
 * @param[in] mext3 mass of external field V[3]
 * @param[in] mV4 mass of internal field V[4]
 * @param[in] mV5 mass of internal field V[5]
 * @param[in] mV6 mass of internal field V[6]
 * @param[in] CpS1cV4cV5 coupling Cp[S[1], -V[4], -V[5]][g[lt2, lt3]]
 * @param[in] CpS2V4cV6 coupling Cp[S[2], V[4], -V[6]][g[lt2, lt3]]
 * @param[in] CpV3V5V6 coupling Cp[V[3], V[5], V[6]][g[lt1, lt2] (-Mom[V[3]] +
    Mom[V[5]]) + g[lt1, lt3] (Mom[V[3]] - Mom[V[6]]) + g[lt2, lt3] (-Mom[V[5]]
    + Mom[V[6]])]
 *
 * @return value of the one-loop diagram
 */
Decay_amplitude_SSV calculate_diagram_SSV_t1g10n10_VVV(
   double mext1, double mext2, double mext3,
   double mV4, double mV5, double mV6,
   const std::complex<double>& CpS1cV4cV5, const std::complex<double>&
      CpS2V4cV6, const std::complex<double>& CpV3V5V6,
   double scale)
{
   const auto c0tmp1 = passarino_veltman::C0(mext2*mext2, mext3*mext3, mext1*
      mext1, mV4*mV4, mV6*mV6, mV5*mV5, scale*scale);
   const auto c1tmp2 = passarino_veltman::C1(mext2*mext2, mext3*mext3, mext1*
      mext1, mV4*mV4, mV6*mV6, mV5*mV5, scale*scale);
   const auto c2tmp3 = passarino_veltman::C2(mext2*mext2, mext3*mext3, mext1*
      mext1, mV4*mV4, mV6*mV6, mV5*mV5, scale*scale);

   Decay_amplitude_SSV result;

   result.m_decay = mext1;
   result.m_scalar = mext2;
   result.m_vector = mext3;

   result.form_factor = oneOver16PiSqr*(std::complex<double>(0,-6)*(c0tmp1 +
      c1tmp2 + c2tmp3)*CpS1cV4cV5*CpS2V4cV6*CpV3V5V6);

   return result;
}

/**
 * @brief Evaluates T2G1N11 diagram for process S -> SV
 *
 * @param[in] mext1 mass of external field S[1]
 * @param[in] mext2 mass of external field S[2]
 * @param[in] mext3 mass of external field V[3]
 * @param[in] mS4 mass of internal field S[4]
 * @param[in] mV5 mass of internal field V[5]
 * @param[in] CpS1cS4cV5 coupling Cp[S[1], -S[4], -V[5]][Mom[S[1]] - Mom[-S[4]]
    ]
 * @param[in] CpS2S4V3V5 coupling Cp[S[2], S[4], V[3], V[5]][g[lt3, lt4]]
 *
 * @return value of the one-loop diagram
 */
Decay_amplitude_SSV calculate_diagram_SSV_t2g1n11_SV(
   double mext1, double mext2, double mext3,
   double mS4, double mV5,
   const std::complex<double>& CpS1cS4cV5, const std::complex<double>&
      CpS2S4V3V5,
   double scale)
{
   const auto b0tmp1 = passarino_veltman::B0(mext1*mext1, mS4*mS4, mV5*mV5,
      scale*scale);
   const auto b1tmp2 = passarino_veltman::B1(mext1*mext1, mS4*mS4, mV5*mV5,
      scale*scale);

   Decay_amplitude_SSV result;

   result.m_decay = mext1;
   result.m_scalar = mext2;
   result.m_vector = mext3;

   result.form_factor = oneOver16PiSqr*((b0tmp1 - b1tmp2)*CpS1cS4cV5*CpS2S4V3V5
      );

   return result;
}

/**
 * @brief Evaluates T3G1N12 diagram for process S -> SV
 *
 * @param[in] mext1 mass of external field S[1]
 * @param[in] mext2 mass of external field S[2]
 * @param[in] mext3 mass of external field V[3]
 * @param[in] mS4 mass of internal field S[4]
 * @param[in] mV5 mass of internal field V[5]
 * @param[in] CpS2cS4cV5 coupling Cp[S[2], -S[4], -V[5]][Mom[S[2]] - Mom[-S[4]]
    ]
 * @param[in] CpS1S4V3V5 coupling Cp[S[1], S[4], V[3], V[5]][g[lt3, lt4]]
 *
 * @return value of the one-loop diagram
 */
Decay_amplitude_SSV calculate_diagram_SSV_t3g1n12_SV(
   double mext1, double mext2, double mext3,
   double mS4, double mV5,
   const std::complex<double>& CpS2cS4cV5, const std::complex<double>&
      CpS1S4V3V5,
   double scale)
{
   const auto b0tmp1 = passarino_veltman::B0(mext2*mext2, mS4*mS4, mV5*mV5,
      scale*scale);
   const auto b1tmp2 = passarino_veltman::B1(mext2*mext2, mS4*mS4, mV5*mV5,
      scale*scale);

   Decay_amplitude_SSV result;

   result.m_decay = mext1;
   result.m_scalar = mext2;
   result.m_vector = mext3;

   result.form_factor = oneOver16PiSqr*((-b0tmp1 + b1tmp2)*CpS1S4V3V5*
      CpS2cS4cV5);

   return result;
}

/**
 * @brief Evaluates T1G1N1 diagram for process S -> FF
 *
 * @param[in] mext1 mass of external field S[1]
 * @param[in] mext2 mass of external field F[2]
 * @param[in] mext3 mass of external field F[3]
 * @param[in] mF4 mass of internal field F[4]
 * @param[in] mF5 mass of internal field F[5]
 * @param[in] mS6 mass of internal field S[6]
 * @param[in] CpF2F4cS6PL coupling Cp[F[2], F[4], -S[6]][PL]
 * @param[in] CpF2F4cS6PR coupling Cp[F[2], F[4], -S[6]][PR]
 * @param[in] CpcF4cF5S1PL coupling Cp[-F[4], -F[5], S[1]][PL]
 * @param[in] CpcF4cF5S1PR coupling Cp[-F[4], -F[5], S[1]][PR]
 * @param[in] CpF5F3S6PL coupling Cp[F[5], F[3], S[6]][PL]
 * @param[in] CpF5F3S6PR coupling Cp[F[5], F[3], S[6]][PR]
 *
 * @return value of the one-loop diagram
 */
Decay_amplitude_SFF calculate_diagram_SFF_t1g1n1_FFS(
   double mext1, double mext2, double mext3,
   double mF4, double mF5, double mS6,
   const std::complex<double>& CpF2F4cS6PL, const std::complex<double>&
      CpF2F4cS6PR, const std::complex<double>& CpcF4cF5S1PL, const std::complex
      <double>& CpcF4cF5S1PR, const std::complex<double>& CpF5F3S6PL, const std
      ::complex<double>& CpF5F3S6PR,
   double scale)
{
   const auto b0tmp1 = passarino_veltman::B0(mext3*mext3, mF5*mF5, mS6*mS6,
      scale*scale);
   const auto c0tmp2 = passarino_veltman::C0(mext1*mext1, mext3*mext3, mext2*
      mext2, mF4*mF4, mF5*mF5, mS6*mS6, scale*scale);
   const auto c1tmp3 = passarino_veltman::C1(mext1*mext1, mext3*mext3, mext2*
      mext2, mF4*mF4, mF5*mF5, mS6*mS6, scale*scale);
   const auto c2tmp4 = passarino_veltman::C2(mext1*mext1, mext3*mext3, mext2*
      mext2, mF4*mF4, mF5*mF5, mS6*mS6, scale*scale);

   Decay_amplitude_SFF result;

   result.m_decay = mext1;
   result.m_out_1 = mext2;
   result.m_out_2 = mext3;

   result.form_factor_right = oneOver16PiSqr*(std::complex<double>(0,1)*(b0tmp1
      *CpcF4cF5S1PL*CpF2F4cS6PR*CpF5F3S6PR - c0tmp2*CpcF4cF5S1PL*CpF2F4cS6PL*
      CpF5F3S6PR*mext2*mF4 + c0tmp2*CpcF4cF5S1PR*CpF2F4cS6PR*CpF5F3S6PL*mext3*
      mF4 + c0tmp2*CpcF4cF5S1PR*CpF2F4cS6PR*CpF5F3S6PR*mF4*mF5 + c1tmp3*
      CpF2F4cS6PR*CpF5F3S6PL*mext3*(CpcF4cF5S1PR*mF4 + CpcF4cF5S1PL*mF5) -
      c2tmp4*CpF2F4cS6PL*mext2*(CpcF4cF5S1PR*CpF5F3S6PL*mext3 + CpcF4cF5S1PL*
      CpF5F3S6PR*mF4 + CpcF4cF5S1PR*CpF5F3S6PR*mF5) + c1tmp3*CpF5F3S6PR*(-(
      CpF2F4cS6PL*mext2*(CpcF4cF5S1PL*mF4 + CpcF4cF5S1PR*mF5)) + CpcF4cF5S1PL*
      CpF2F4cS6PR*Sqr(mext1)) + c2tmp4*CpcF4cF5S1PL*CpF2F4cS6PR*CpF5F3S6PR*Sqr(
      mext2) + c0tmp2*CpcF4cF5S1PL*CpF2F4cS6PR*CpF5F3S6PR*Sqr(mF4)));
   result.form_factor_left = oneOver16PiSqr*(std::complex<double>(0,1)*(b0tmp1*
      CpcF4cF5S1PR*CpF2F4cS6PL*CpF5F3S6PL - c2tmp4*CpcF4cF5S1PL*CpF2F4cS6PR*
      CpF5F3S6PR*mext2*mext3 - c2tmp4*CpcF4cF5S1PR*CpF2F4cS6PR*CpF5F3S6PL*mext2
      *mF4 - c2tmp4*CpcF4cF5S1PL*CpF2F4cS6PR*CpF5F3S6PL*mext2*mF5 - c1tmp3*
      CpF2F4cS6PR*CpF5F3S6PL*mext2*(CpcF4cF5S1PR*mF4 + CpcF4cF5S1PL*mF5) +
      c0tmp2*mF4*(-(CpcF4cF5S1PR*CpF2F4cS6PR*CpF5F3S6PL*mext2) + CpcF4cF5S1PL*
      CpF2F4cS6PL*CpF5F3S6PR*mext3 + CpcF4cF5S1PR*CpF2F4cS6PL*CpF5F3S6PL*mF4 +
      CpcF4cF5S1PL*CpF2F4cS6PL*CpF5F3S6PL*mF5) + c1tmp3*CpF2F4cS6PL*(
      CpcF4cF5S1PL*CpF5F3S6PR*mext3*mF4 + CpcF4cF5S1PR*CpF5F3S6PR*mext3*mF5 +
      CpcF4cF5S1PR*CpF5F3S6PL*Sqr(mext1)) + c2tmp4*CpcF4cF5S1PR*CpF2F4cS6PL*
      CpF5F3S6PL*Sqr(mext2)));

   return result;
}

/**
 * @brief Evaluates T1G2N2 diagram for process S -> FF
 *
 * @param[in] mext1 mass of external field S[1]
 * @param[in] mext2 mass of external field F[2]
 * @param[in] mext3 mass of external field F[3]
 * @param[in] mS4 mass of internal field S[4]
 * @param[in] mS5 mass of internal field S[5]
 * @param[in] mF6 mass of internal field F[6]
 * @param[in] CpF2cF6S4PL coupling Cp[F[2], -F[6], S[4]][PL]
 * @param[in] CpF2cF6S4PR coupling Cp[F[2], -F[6], S[4]][PR]
 * @param[in] CpF6F3S5PL coupling Cp[F[6], F[3], S[5]][PL]
 * @param[in] CpF6F3S5PR coupling Cp[F[6], F[3], S[5]][PR]
 * @param[in] CpS1cS4cS5 coupling Cp[S[1], -S[4], -S[5]]
 *
 * @return value of the one-loop diagram
 */
Decay_amplitude_SFF calculate_diagram_SFF_t1g2n2_SSF(
   double mext1, double mext2, double mext3,
   double mS4, double mS5, double mF6,
   const std::complex<double>& CpF2cF6S4PL, const std::complex<double>&
      CpF2cF6S4PR, const std::complex<double>& CpF6F3S5PL, const std::complex<
      double>& CpF6F3S5PR, const std::complex<double>& CpS1cS4cS5,
   double scale)
{
   const auto c0tmp1 = passarino_veltman::C0(mext2*mext2, mext1*mext1, mext3*
      mext3, mF6*mF6, mS4*mS4, mS5*mS5, scale*scale);
   const auto c1tmp2 = passarino_veltman::C1(mext2*mext2, mext1*mext1, mext3*
      mext3, mF6*mF6, mS4*mS4, mS5*mS5, scale*scale);
   const auto c2tmp3 = passarino_veltman::C2(mext2*mext2, mext1*mext1, mext3*
      mext3, mF6*mF6, mS4*mS4, mS5*mS5, scale*scale);

   Decay_amplitude_SFF result;

   result.m_decay = mext1;
   result.m_out_1 = mext2;
   result.m_out_2 = mext3;

   result.form_factor_right = oneOver16PiSqr*(std::complex<double>(0,-1)*
      CpS1cS4cS5*(c1tmp2*CpF2cF6S4PL*CpF6F3S5PR*mext2 + c2tmp3*CpF2cF6S4PR*
      CpF6F3S5PL*mext3 - c0tmp1*CpF2cF6S4PR*CpF6F3S5PR*mF6));
   result.form_factor_left = oneOver16PiSqr*(std::complex<double>(0,-1)*
      CpS1cS4cS5*(c1tmp2*CpF2cF6S4PR*CpF6F3S5PL*mext2 + c2tmp3*CpF2cF6S4PL*
      CpF6F3S5PR*mext3 - c0tmp1*CpF2cF6S4PL*CpF6F3S5PL*mF6));

   return result;
}

/**
 * @brief Evaluates T1G3N3 diagram for process S -> FF
 *
 * @param[in] mext1 mass of external field S[1]
 * @param[in] mext2 mass of external field F[2]
 * @param[in] mext3 mass of external field F[3]
 * @param[in] mF4 mass of internal field F[4]
 * @param[in] mF5 mass of internal field F[5]
 * @param[in] mV6 mass of internal field V[6]
 * @param[in] CpF2F4cV6PL coupling Cp[F[2], F[4], -V[6]][LorentzProduct[gamma[
    lt3], PL]]
 * @param[in] CpF2F4cV6PR coupling Cp[F[2], F[4], -V[6]][LorentzProduct[gamma[
    lt3], PR]]
 * @param[in] CpcF4cF5S1PL coupling Cp[-F[4], -F[5], S[1]][PL]
 * @param[in] CpcF4cF5S1PR coupling Cp[-F[4], -F[5], S[1]][PR]
 * @param[in] CpF5F3V6PR coupling Cp[F[5], F[3], V[6]][LorentzProduct[gamma[lt3
    ], PR]]
 * @param[in] CpF5F3V6PL coupling Cp[F[5], F[3], V[6]][LorentzProduct[gamma[lt3
    ], PL]]
 *
 * @return value of the one-loop diagram
 */
Decay_amplitude_SFF calculate_diagram_SFF_t1g3n3_FFV(
   double mext1, double mext2, double mext3,
   double mF4, double mF5, double mV6,
   const std::complex<double>& CpF2F4cV6PL, const std::complex<double>&
      CpF2F4cV6PR, const std::complex<double>& CpcF4cF5S1PL, const std::complex
      <double>& CpcF4cF5S1PR, const std::complex<double>& CpF5F3V6PR, const std
      ::complex<double>& CpF5F3V6PL,
   double scale, double finite)
{
   const auto b0tmp1 = passarino_veltman::B0(mext3*mext3, mF5*mF5, mV6*mV6,
      scale*scale);
   const auto c0tmp2 = passarino_veltman::C0(mext1*mext1, mext3*mext3, mext2*
      mext2, mF4*mF4, mF5*mF5, mV6*mV6, scale*scale);
   const auto c1tmp3 = passarino_veltman::C1(mext1*mext1, mext3*mext3, mext2*
      mext2, mF4*mF4, mF5*mF5, mV6*mV6, scale*scale);
   const auto c2tmp4 = passarino_veltman::C2(mext1*mext1, mext3*mext3, mext2*
      mext2, mF4*mF4, mF5*mF5, mV6*mV6, scale*scale);

   Decay_amplitude_SFF result;

   result.m_decay = mext1;
   result.m_out_1 = mext2;
   result.m_out_2 = mext3;

   result.form_factor_right = oneOver16PiSqr*(std::complex<double>(0,2)*(-2*
      b0tmp1*CpcF4cF5S1PR*CpF2F4cV6PL*CpF5F3V6PR + CpcF4cF5S1PL*(c1tmp3*
      CpF2F4cV6PL*CpF5F3V6PL*mext3*mF4 - c1tmp3*CpF2F4cV6PR*CpF5F3V6PR*mext2*
      mF5 - c2tmp4*CpF2F4cV6PR*CpF5F3V6PR*mext2*mF5 + c0tmp2*CpF2F4cV6PL*mF4*(
      CpF5F3V6PL*mext3 - 2*CpF5F3V6PR*mF5)) + CpcF4cF5S1PR*(-((c0tmp2 + c1tmp3
      + c2tmp4)*CpF2F4cV6PR*CpF5F3V6PR*mext2*mF4) + c1tmp3*CpF2F4cV6PL*
      CpF5F3V6PL*mext3*mF5 + CpF2F4cV6PL*CpF5F3V6PR*(finite - 2*c1tmp3*Sqr(
      mext1) - c2tmp4*Sqr(mext1) - c2tmp4*Sqr(mext2) + c2tmp4*Sqr(mext3) - 2*
      c0tmp2*Sqr(mF4)))));
   result.form_factor_left = oneOver16PiSqr*(std::complex<double>(0,-2)*(2*
      b0tmp1*CpcF4cF5S1PL*CpF2F4cV6PR*CpF5F3V6PL + CpcF4cF5S1PR*(-(c1tmp3*
      CpF2F4cV6PR*CpF5F3V6PR*mext3*mF4) + c1tmp3*CpF2F4cV6PL*CpF5F3V6PL*mext2*
      mF5 + c2tmp4*CpF2F4cV6PL*CpF5F3V6PL*mext2*mF5 + c0tmp2*CpF2F4cV6PR*mF4*(-
      (CpF5F3V6PR*mext3) + 2*CpF5F3V6PL*mF5)) + CpcF4cF5S1PL*((c0tmp2 + c1tmp3
      + c2tmp4)*CpF2F4cV6PL*CpF5F3V6PL*mext2*mF4 - c1tmp3*CpF2F4cV6PR*
      CpF5F3V6PR*mext3*mF5 + CpF2F4cV6PR*CpF5F3V6PL*(-finite + 2*c1tmp3*Sqr(
      mext1) + c2tmp4*Sqr(mext1) + c2tmp4*Sqr(mext2) - c2tmp4*Sqr(mext3) + 2*
      c0tmp2*Sqr(mF4)))));

   return result;
}

/**
 * @brief Evaluates T1G4N4 diagram for process S -> FF
 *
 * @param[in] mext1 mass of external field S[1]
 * @param[in] mext2 mass of external field F[2]
 * @param[in] mext3 mass of external field F[3]
 * @param[in] mS4 mass of internal field S[4]
 * @param[in] mV5 mass of internal field V[5]
 * @param[in] mF6 mass of internal field F[6]
 * @param[in] CpF2cF6S4PL coupling Cp[F[2], -F[6], S[4]][PL]
 * @param[in] CpF2cF6S4PR coupling Cp[F[2], -F[6], S[4]][PR]
 * @param[in] CpF6F3V5PL coupling Cp[F[6], F[3], V[5]][LorentzProduct[gamma[lt3
    ], PL]]
 * @param[in] CpF6F3V5PR coupling Cp[F[6], F[3], V[5]][LorentzProduct[gamma[lt3
    ], PR]]
 * @param[in] CpS1cS4cV5 coupling Cp[S[1], -S[4], -V[5]][Mom[S[1]] - Mom[-S[4]]
    ]
 *
 * @return value of the one-loop diagram
 */
Decay_amplitude_SFF calculate_diagram_SFF_t1g4n4_SVF(
   double mext1, double mext2, double mext3,
   double mS4, double mV5, double mF6,
   const std::complex<double>& CpF2cF6S4PL, const std::complex<double>&
      CpF2cF6S4PR, const std::complex<double>& CpF6F3V5PL, const std::complex<
      double>& CpF6F3V5PR, const std::complex<double>& CpS1cS4cV5,
   double scale)
{
   const auto b0tmp1 = passarino_veltman::B0(mext3*mext3, mF6*mF6, mV5*mV5,
      scale*scale);
   const auto c0tmp2 = passarino_veltman::C0(mext2*mext2, mext1*mext1, mext3*
      mext3, mF6*mF6, mS4*mS4, mV5*mV5, scale*scale);
   const auto c1tmp3 = passarino_veltman::C1(mext2*mext2, mext1*mext1, mext3*
      mext3, mF6*mF6, mS4*mS4, mV5*mV5, scale*scale);
   const auto c2tmp4 = passarino_veltman::C2(mext2*mext2, mext1*mext1, mext3*
      mext3, mF6*mF6, mS4*mS4, mV5*mV5, scale*scale);

   Decay_amplitude_SFF result;

   result.m_decay = mext1;
   result.m_out_1 = mext2;
   result.m_out_2 = mext3;

   result.form_factor_right = oneOver16PiSqr*(std::complex<double>(0,-1)*
      CpS1cS4cV5*(-(b0tmp1*CpF2cF6S4PR*CpF6F3V5PR) + 2*c2tmp4*CpF2cF6S4PL*
      CpF6F3V5PL*mext2*mext3 + 2*c0tmp2*CpF2cF6S4PL*CpF6F3V5PR*mext2*mF6 -
      c0tmp2*CpF2cF6S4PR*CpF6F3V5PL*mext3*mF6 + c1tmp3*CpF2cF6S4PL*mext2*(
      CpF6F3V5PL*mext3 + CpF6F3V5PR*mF6) + c2tmp4*CpF2cF6S4PR*(CpF6F3V5PL*mext3
      *mF6 + CpF6F3V5PR*(Sqr(mext1) - Sqr(mext2))) + c0tmp2*CpF2cF6S4PR*
      CpF6F3V5PR*Sqr(mext2) - c0tmp2*CpF2cF6S4PR*CpF6F3V5PR*Sqr(mS4)));
   result.form_factor_left = oneOver16PiSqr*(std::complex<double>(0,-1)*
      CpS1cS4cV5*(-(b0tmp1*CpF2cF6S4PL*CpF6F3V5PL) + 2*c2tmp4*CpF2cF6S4PR*
      CpF6F3V5PR*mext2*mext3 + 2*c0tmp2*CpF2cF6S4PR*CpF6F3V5PL*mext2*mF6 -
      c0tmp2*CpF2cF6S4PL*CpF6F3V5PR*mext3*mF6 + c1tmp3*CpF2cF6S4PR*mext2*(
      CpF6F3V5PR*mext3 + CpF6F3V5PL*mF6) + c2tmp4*CpF2cF6S4PL*(CpF6F3V5PR*mext3
      *mF6 + CpF6F3V5PL*(Sqr(mext1) - Sqr(mext2))) + c0tmp2*CpF2cF6S4PL*
      CpF6F3V5PL*Sqr(mext2) - c0tmp2*CpF2cF6S4PL*CpF6F3V5PL*Sqr(mS4)));

   return result;
}

/**
 * @brief Evaluates T1G5N5 diagram for process S -> FF
 *
 * @param[in] mext1 mass of external field S[1]
 * @param[in] mext2 mass of external field F[2]
 * @param[in] mext3 mass of external field F[3]
 * @param[in] mV4 mass of internal field V[4]
 * @param[in] mS5 mass of internal field S[5]
 * @param[in] mF6 mass of internal field F[6]
 * @param[in] CpF2cF6V4PR coupling Cp[F[2], -F[6], V[4]][LorentzProduct[gamma[
    lt3], PR]]
 * @param[in] CpF2cF6V4PL coupling Cp[F[2], -F[6], V[4]][LorentzProduct[gamma[
    lt3], PL]]
 * @param[in] CpF6F3S5PL coupling Cp[F[6], F[3], S[5]][PL]
 * @param[in] CpF6F3S5PR coupling Cp[F[6], F[3], S[5]][PR]
 * @param[in] CpS1cS5cV4 coupling Cp[S[1], -S[5], -V[4]][Mom[S[1]] - Mom[-S[5]]
    ]
 *
 * @return value of the one-loop diagram
 */
Decay_amplitude_SFF calculate_diagram_SFF_t1g5n5_VSF(
   double mext1, double mext2, double mext3,
   double mV4, double mS5, double mF6,
   const std::complex<double>& CpF2cF6V4PR, const std::complex<double>&
      CpF2cF6V4PL, const std::complex<double>& CpF6F3S5PL, const std::complex<
      double>& CpF6F3S5PR, const std::complex<double>& CpS1cS5cV4,
   double scale)
{
   const auto b0tmp1 = passarino_veltman::B0(mext3*mext3, mF6*mF6, mS5*mS5,
      scale*scale);
   const auto c0tmp2 = passarino_veltman::C0(mext2*mext2, mext1*mext1, mext3*
      mext3, mF6*mF6, mV4*mV4, mS5*mS5, scale*scale);
   const auto c1tmp3 = passarino_veltman::C1(mext2*mext2, mext1*mext1, mext3*
      mext3, mF6*mF6, mV4*mV4, mS5*mS5, scale*scale);
   const auto c2tmp4 = passarino_veltman::C2(mext2*mext2, mext1*mext1, mext3*
      mext3, mF6*mF6, mV4*mV4, mS5*mS5, scale*scale);

   Decay_amplitude_SFF result;

   result.m_decay = mext1;
   result.m_out_1 = mext2;
   result.m_out_2 = mext3;

   result.form_factor_right = oneOver16PiSqr*(std::complex<double>(0,1)*
      CpS1cS5cV4*(-(b0tmp1*CpF2cF6V4PL*CpF6F3S5PR) + 2*c1tmp3*CpF2cF6V4PR*
      CpF6F3S5PL*mext2*mext3 + c1tmp3*CpF2cF6V4PR*CpF6F3S5PR*mext2*mF6 + 2*
      c0tmp2*CpF2cF6V4PL*CpF6F3S5PL*mext3*mF6 + c2tmp4*CpF6F3S5PL*mext3*(
      CpF2cF6V4PR*mext2 + CpF2cF6V4PL*mF6) + 2*c1tmp3*CpF2cF6V4PL*CpF6F3S5PR*
      Sqr(mext1) + c1tmp3*CpF2cF6V4PL*CpF6F3S5PR*Sqr(mext2) - 2*c1tmp3*
      CpF2cF6V4PL*CpF6F3S5PR*Sqr(mext3) - c2tmp4*CpF2cF6V4PL*CpF6F3S5PR*(Sqr(
      mext1) - Sqr(mext2) + Sqr(mext3)) - c0tmp2*CpF6F3S5PR*(CpF2cF6V4PR*mext2*
      mF6 + CpF2cF6V4PL*(-Sqr(mext2) + Sqr(mV4)))));
   result.form_factor_left = oneOver16PiSqr*(std::complex<double>(0,1)*
      CpS1cS5cV4*(-(b0tmp1*CpF2cF6V4PR*CpF6F3S5PL) + c2tmp4*CpF2cF6V4PL*
      CpF6F3S5PR*mext2*mext3 - c0tmp2*CpF2cF6V4PL*CpF6F3S5PL*mext2*mF6 + c2tmp4
      *CpF2cF6V4PR*CpF6F3S5PR*mext3*mF6 + c1tmp3*CpF2cF6V4PL*mext2*(2*
      CpF6F3S5PR*mext3 + CpF6F3S5PL*mF6) - c2tmp4*CpF2cF6V4PR*CpF6F3S5PL*Sqr(
      mext1) + c2tmp4*CpF2cF6V4PR*CpF6F3S5PL*Sqr(mext2) + c1tmp3*CpF2cF6V4PR*
      CpF6F3S5PL*(2*Sqr(mext1) + Sqr(mext2) - 2*Sqr(mext3)) - c2tmp4*
      CpF2cF6V4PR*CpF6F3S5PL*Sqr(mext3) + c0tmp2*CpF2cF6V4PR*(2*CpF6F3S5PR*
      mext3*mF6 + CpF6F3S5PL*(Sqr(mext2) - Sqr(mV4)))));

   return result;
}

/**
 * @brief Evaluates T1G6N6 diagram for process S -> FF
 *
 * @param[in] mext1 mass of external field S[1]
 * @param[in] mext2 mass of external field F[2]
 * @param[in] mext3 mass of external field F[3]
 * @param[in] mV4 mass of internal field V[4]
 * @param[in] mV5 mass of internal field V[5]
 * @param[in] mF6 mass of internal field F[6]
 * @param[in] CpF2cF6V4PR coupling Cp[F[2], -F[6], V[4]][LorentzProduct[gamma[
    lt3], PR]]
 * @param[in] CpF2cF6V4PL coupling Cp[F[2], -F[6], V[4]][LorentzProduct[gamma[
    lt3], PL]]
 * @param[in] CpF6F3V5PR coupling Cp[F[6], F[3], V[5]][LorentzProduct[gamma[lt3
    ], PR]]
 * @param[in] CpF6F3V5PL coupling Cp[F[6], F[3], V[5]][LorentzProduct[gamma[lt3
    ], PL]]
 * @param[in] CpS1cV4cV5 coupling Cp[S[1], -V[4], -V[5]][g[lt2, lt3]]
 *
 * @return value of the one-loop diagram
 */
Decay_amplitude_SFF calculate_diagram_SFF_t1g6n6_VVF(
   double mext1, double mext2, double mext3,
   double mV4, double mV5, double mF6,
   const std::complex<double>& CpF2cF6V4PR, const std::complex<double>&
      CpF2cF6V4PL, const std::complex<double>& CpF6F3V5PR, const std::complex<
      double>& CpF6F3V5PL, const std::complex<double>& CpS1cV4cV5,
   double scale)
{
   const auto c0tmp1 = passarino_veltman::C0(mext2*mext2, mext1*mext1, mext3*
      mext3, mF6*mF6, mV4*mV4, mV5*mV5, scale*scale);
   const auto c1tmp2 = passarino_veltman::C1(mext2*mext2, mext1*mext1, mext3*
      mext3, mF6*mF6, mV4*mV4, mV5*mV5, scale*scale);
   const auto c2tmp3 = passarino_veltman::C2(mext2*mext2, mext1*mext1, mext3*
      mext3, mF6*mF6, mV4*mV4, mV5*mV5, scale*scale);

   Decay_amplitude_SFF result;

   result.m_decay = mext1;
   result.m_out_1 = mext2;
   result.m_out_2 = mext3;

   result.form_factor_right = oneOver16PiSqr*(std::complex<double>(0,2)*
      CpS1cV4cV5*(c1tmp2*CpF2cF6V4PR*CpF6F3V5PR*mext2 + c2tmp3*CpF2cF6V4PL*
      CpF6F3V5PL*mext3 + 2*c0tmp1*CpF2cF6V4PL*CpF6F3V5PR*mF6));
   result.form_factor_left = oneOver16PiSqr*(std::complex<double>(0,2)*
      CpS1cV4cV5*(c1tmp2*CpF2cF6V4PL*CpF6F3V5PL*mext2 + c2tmp3*CpF2cF6V4PR*
      CpF6F3V5PR*mext3 + 2*c0tmp1*CpF2cF6V4PR*CpF6F3V5PL*mF6));

   return result;
}

} // namespace flexiblesusy
