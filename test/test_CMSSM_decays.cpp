
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_CMSSM_decays

#include <boost/test/unit_test.hpp>

#include "wrappers.hpp"

#include "CMSSM_decays.hpp"
#include "CMSSM_mass_eigenstates.hpp"
#include "test_CMSSM.hpp"

namespace flexiblesusy {

double hh_to_barFuFu_leading_order(
   const CMSSM_mass_eigenstates& model, int gI1, int gO1)
{
   const int N_c = 3;

   const auto vu = model.get_vu();
   const auto vd = model.get_vd();
   const auto v = Sqrt(vu * vu + vd * vd);
   const auto tb = vu / vd;
   const auto sb = Sin(ArcTan(tb));

   const auto G_mu = 1. / (Sqrt(2.) * v * v);

   const auto ca = Cos(model.Alpha());

   const auto m_f = model.get_MFu(gO1);
   const auto M_phi = model.get_Mhh(gI1);

   const auto coup_SM = Sqrt(Sqrt(2.) * G_mu) * m_f;
   const auto coup = ca * coup_SM / sb;
   const auto beta_f = Sqrt(1. - 4. * m_f * m_f / (M_phi * M_phi));

   return N_c * G_mu * Sqr(m_f) * Sqr(coup) * M_phi * Cube(beta_f) /
      (4. * Sqrt(2.) * Pi);
}

double hh_to_barFdFd_leading_order(
   const CMSSM_mass_eigenstates& model, int gI1, int gO1)
{
   const int N_c = 3;

   const auto vu = model.get_vu();
   const auto vd = model.get_vd();
   const auto v = Sqrt(vu * vu + vd * vd);
   const auto tb = vu / vd;
   const auto cb = Cos(ArcTan(tb));

   const auto G_mu = 1. / (Sqrt(2.) * v * v);

   const auto sa = Sin(model.Alpha());

   const auto m_f = model.get_MFd(gO1);
   const auto M_phi = model.get_Mhh(gI1);

   const auto coup_SM = Sqrt(Sqrt(2.) * G_mu) * m_f;
   const auto coup = -sa * coup_SM / cb;
   const auto beta_f = Sqrt(1. - 4. * m_f * m_f / (M_phi * M_phi));

   return N_c * G_mu * Sqr(m_f) * Sqr(coup) * M_phi * Cube(beta_f) /
      (4. * Sqrt(2.) * Pi);
}

double hh_to_barFeFe_leading_order(
   const CMSSM_mass_eigenstates& model, int gI1, int gO1)
{
   const int N_c = 1;

   const auto vu = model.get_vu();
   const auto vd = model.get_vd();
   const auto v = Sqrt(vu * vu + vd * vd);
   const auto tb = vu / vd;
   const auto cb = Cos(ArcTan(tb));

   const auto G_mu = 1. / (Sqrt(2.) * v * v);

   const auto sa = Sin(model.Alpha());

   const auto m_f = model.get_MFd(gO1);
   const auto M_phi = model.get_Mhh(gI1);

   const auto coup_SM = Sqrt(Sqrt(2.) * G_mu) * m_f;
   const auto coup = -sa * coup_SM / cb;
   const auto beta_f = Sqrt(1. - 4. * m_f * m_f / (M_phi * M_phi));

   return N_c * G_mu * Sqr(m_f) * Sqr(coup) * M_phi * Cube(beta_f) /
      (4. * Sqrt(2.) * Pi);
}

} // namespace flexiblesusy

using namespace flexiblesusy;

BOOST_AUTO_TEST_CASE( test_CMSSM_hh_to_SM_fermions_leading_order )
{
   const softsusy::QedQcd qedqcd;

   CMSSM_input_parameters input;
   input.m0 = 500.;
   input.m12 = 125.;
   input.TanBeta = 10.;
   input.SignMu = 1.;
   input.Azero = 50.;

   CMSSM_mass_eigenstates model(setup_CMSSM(input));
   model.do_calculate_sm_pole_masses(true);
   model.set_pole_mass_loop_order(0);
   model.calculate_pole_masses();

   const auto hh1_to_barFuFu_expected = hh_to_barFuFu_leading_order(
      model, 0, 2);
   const auto hh2_to_barFuFu_expected = hh_to_barFuFu_leading_order(
      model, 1, 2);
   const auto hh1_to_barFdFd_expected = hh_to_barFdFd_leading_order(
      model, 0, 2);
   const auto hh2_to_barFdFd_expected = hh_to_barFdFd_leading_order(
      model, 1, 2);
   const auto hh1_to_barFeFe_expected = hh_to_barFeFe_leading_order(
      model, 0, 2);
   const auto hh2_to_barFeFe_expected = hh_to_barFeFe_leading_order(
      model, 1, 2);

   CMSSM_decays decays(model, qedqcd, input, HigherOrderSMCorrections::disable);

   const auto hh1_to_barFuFu = decays.partial_width_hh_to_barFuFu(
      &model, 0, 2, 2);
   const auto hh2_to_barFuFu = decays.partial_width_hh_to_barFuFu(
      &model, 1, 2, 2);
   const auto hh1_to_barFdFd = decays.partial_width_hh_to_barFdFd(
      &model, 0, 2, 2);
   const auto hh2_to_barFdFd = decays.partial_width_hh_to_barFdFd(
      &model, 1, 2, 2);
   const auto hh1_to_barFeFe = decays.partial_width_hh_to_barFeFe(
      &model, 0, 2, 2);
   const auto hh2_to_barFeFe = decays.partial_width_hh_to_barFeFe(
      &model, 1, 2, 2);

   BOOST_CHECK_CLOSE_FRACTION(hh1_to_barFuFu, hh1_to_barFuFu_expected, 1e-15);
   BOOST_CHECK_CLOSE_FRACTION(hh2_to_barFuFu, hh2_to_barFuFu_expected, 1e-15);
   BOOST_CHECK_CLOSE_FRACTION(hh1_to_barFdFd, hh1_to_barFdFd_expected, 1e-15);
   BOOST_CHECK_CLOSE_FRACTION(hh2_to_barFdFd, hh2_to_barFdFd_expected, 1e-15);
   BOOST_CHECK_CLOSE_FRACTION(hh1_to_barFeFe, hh1_to_barFeFe_expected, 1e-15);
   BOOST_CHECK_CLOSE_FRACTION(hh2_to_barFeFe, hh2_to_barFeFe_expected, 1e-15);

}

