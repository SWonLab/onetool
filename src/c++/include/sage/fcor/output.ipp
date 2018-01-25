// ---------------------------------------------------------------------------
// Inline Implementation of FcorViewBase
// ---------------------------------------------------------------------------

inline size_t
FcorView::trait_count() const
{
  return my_pairsetdata->get_mped()->info().trait_count();
}

inline string
FcorView::get_trait_name(size_t t) const
{
  return my_pairsetdata->get_mped()->info().trait_info(t).name();
}

inline double
FcorView::eff_count(double cor, double se) const
{
  if( SAGE::isnan(cor) || SAGE::isnan(se) )
    return SAGE::QNAN;

  if( fabs(se) < std::numeric_limits<double>::epsilon() )
    return SAGE::QNAN;

  double N = 1.0 + ( 1.0 - cor * cor ) * ( 1.0 - cor * cor ) / ( se * se );

  return ( N + sqrt( N * N + 22.0 * cor * cor * ( 1.0 - cor * cor ) / ( se * se ) ) ) / 2.0;
}

inline const analysis_option_type&
FcorView::get_analysis_options() const
{
  return my_options;
}
