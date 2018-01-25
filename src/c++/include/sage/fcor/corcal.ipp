// ---------------------------------------------------------------------------
// Inline Implementation of CorrelationCal
// ---------------------------------------------------------------------------

inline void
CorrelationCal::set_mped(const RefMultiPedigree* mp)
{
  my_mped = mp;
}

inline bool
CorrelationCal::is_valid() const
{
  return my_valid;
}

inline const RefMultiPedigree*
CorrelationCal::get_mped() const
{
  return my_mped;
}

inline const pairset_vector*
CorrelationCal::get_pairset() const
{
  return my_pairset;
}

inline const pairset_info_vector*
CorrelationCal::get_pairset_info() const
{
  return my_pairset_info;
}

inline const corinfo_vector&
CorrelationCal::get_corinfo() const
{
  return my_corinfo;
}

inline const pair_to_corinfo_map&
CorrelationCal::get_corinfo_map() const
{
  return my_corinfo_map;
}

// end of CorrelationCal Implementation
