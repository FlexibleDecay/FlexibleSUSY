
// special case for H -> H H
template <>
double CLASSNAME::get_partial_width<H, H, H>(
   const ContextName& context,
   typename field_indices<H>::type const& indexIn,
   typename field_indices<H>::type const& indexOut1,
   typename field_indices<H>::type const& indexOut2) const
{
   const double mHIn = context.mass<H>(indexIn);
   const double mH1 = context.mass<H>(indexOut1);
   const double mH2 = context.mass<H>(indexOut2);

   if (mHIn < mH1 + mH2) {
         return 0.;
   }

   const double flux = 1. / (2. * mHIn);
   const double ps = beta(mHIn, mH1, mH2) / (8. * Pi);
   const auto ps_symmetry =
      final_state_symmetry_factor<H, H>(indexOut1, indexOut2);

   return flux * ps * ps_symmetry *
      effective_coupling<H, H, H>(context, indexIn, indexOut1, indexOut2).square();
}