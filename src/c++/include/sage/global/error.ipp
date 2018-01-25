namespace SAGE {

//=============================================
// basic_errorstream
//=============================================

template<class charT,class traits>
inline
basic_errorstream<charT,traits>::basic_errorstream(const sb_type_ptr& sbuffer)
  : std::basic_ostream<charT,traits>(sbuffer.get()), sbuf(sbuffer)
{
  DEBUG(std::cout << "basic_errorstream constructor(sb=" << (void*)sbuffer.get() << ")"<< std::endl;)
  rdbuf(rdbuf());
  DEBUG(std::cout << "basic_errorstream constructor(sb=" << (void*)rdbuf().get() << ")"<< std::endl;)
}

template<class charT,class traits>
basic_errorstream<charT,traits>::basic_errorstream(const basic_errorstream<charT,traits> &s)
  : std::basic_ostream<charT,traits>(s.sbuf.get()), sbuf(s.sbuf)
{
  DEBUG(std::cout << "basic_errorstream constructor(sb=" << (void*)sb << ")"<< std::endl;)
}

template<class charT,class traits>
basic_errorstream<charT,traits>::~basic_errorstream()
{
  DEBUG(std::cout << "deleting basic_errorstream at " << (void*) this << std::endl;)
}

template<class charT,class traits>
basic_errorstream<charT,traits>&
basic_errorstream<charT,traits>::operator=(const basic_errorstream<charT,traits>& s)
{
  if(&s != this)
    rdbuf(s.rdbuf());

  return *this;
}

template<class charT,class traits>
void basic_errorstream<charT,traits>::rdbuf(const sb_type_ptr& sbuffer)
{
  if(rdbuf() == sbuffer && sbuf == sbuffer) return;

  DEBUG(std::cout << "basic_errorstream::rdbuf(sb=" << (void*)sbuffer.get() << ")"<< std::endl;)
  sbuf = sbuffer;

  std::basic_ostream<charT,traits>::rdbuf(sbuffer.get());
}

template<class charT, class traits, class T>
basic_errorstream<char,traits>& operator<<(basic_errorstream<char,traits>& out, const T& t)
{
  static_cast<std::basic_ostream<charT,traits>&>(out) << t;
  return out;
}

template<class charT,class traits>
inline typename basic_errorstream<charT,traits>::sb_type_ptr 
basic_errorstream<charT,traits>::rdbuf() const 
{ return sbuf; }

template<class charT,class traits>
inline typename basic_errorstream<charT,traits>::sb_error_info & 
basic_errorstream<charT,traits>::get_info() 
{ return rdbuf()->get_info(); }

template<class charT,class traits> 
inline typename basic_errorstream<charT,traits>::string_type 
basic_errorstream<charT,traits>::str() const 
{ return rdbuf()->str(); }

template<class charT,class traits>
inline void
basic_errorstream<charT,traits>::str(const string_type& str_arg)
{ rdbuf()->str(str_arg); }

template<class charT,class traits> 
inline void
basic_errorstream<charT,traits>::prefix(const string_type& str_arg)
{ rdbuf()->prefix(str_arg); }

template<class charT,class traits> 
inline typename basic_errorstream<charT,traits>::string_type 
basic_errorstream<charT,traits>::prefix() const
{ return rdbuf()->prefix(); }

template<class charT,class traits> 
inline void
basic_errorstream<charT,traits>::suffix(const string_type& str_arg)
{ rdbuf()->suffix(str_arg); }

template<class charT,class traits> 
inline typename basic_errorstream<charT,traits>::string_type
basic_errorstream<charT,traits>::suffix() const
{ return rdbuf()->suffix(); }

template<class charT,class traits> 
inline void
basic_errorstream<charT,traits>::location(const string_type& file, 
                                                     size_t line)
{ rdbuf()->location(file,line); }

template<class charT,class traits> 
inline void
basic_errorstream<charT,traits>::filename(const string_type &f)
{ rdbuf()->filename(f); }

template<class charT,class traits> 
inline typename basic_errorstream<charT,traits>::string_type
basic_errorstream<charT,traits>::filename() const
{return rdbuf()->filename(); }

template<class charT,class traits>
inline void
basic_errorstream<charT,traits>::linenumber(size_t l)
{ rdbuf()->linenumber(l); }

template<class charT,class traits>
inline size_t
basic_errorstream<charT,traits>::linenumber() const
{ return rdbuf()->linenumber(); }

template<class charT,class traits> 
inline void
basic_errorstream<charT,traits>::priority(error_priority_type pr)
{ rdbuf()->priority(pr); }

template<class charT,class traits> 
inline typename basic_errorstream<charT,traits>::error_priority_type
basic_errorstream<charT,traits>::priority() const
{ return rdbuf()->priority(); }

template<class charT,class traits> 
inline void
basic_errorstream<charT,traits>::offset(size_t o)
{ rdbuf()->offset(o); }

template<class charT,class traits>
inline size_t
basic_errorstream<charT,traits>::offset() const
{ return rdbuf()->offset(); }

template<class charT,class traits>
inline void
basic_errorstream<charT,traits>::indent(size_t i)
{ rdbuf()->indent(i); }

template<class charT,class traits>
inline size_t
basic_errorstream<charT,traits>::indent() const
{ return rdbuf()->indent(); }

template<class charT,class traits>
inline void
basic_errorstream<charT,traits>::line_width(size_t i)
{ rdbuf()->line_width(i); }

template<class charT,class traits>
inline size_t
basic_errorstream<charT,traits>::line_width() const
{ return rdbuf()->line_width(); }

template<class charT,class traits>
inline bool
basic_errorstream<charT,traits>::raw_mode() const
{ return rdbuf()->raw_mode(); }

template<class charT,class traits>
inline void
basic_errorstream<charT,traits>::set_raw_mode()
{ rdbuf()->set_raw_mode(); }

template<class charT,class traits>
inline bool
basic_errorstream<charT,traits>::cooked_mode() const
{ return rdbuf()->cooked_mode(); }

template<class charT,class traits>
inline void
basic_errorstream<charT,traits>::set_cooked_mode()
{ rdbuf()->set_cooked_mode(); }

//=============================================
// errorstream
//=============================================

template<class charT,class traits>
errorstream<charT,traits>::errorstream(const ErrorNoInit& )
  :  basic_errorstream<charT,traits>( basic_errorstream<charT,traits>::rdbuf() )
{ }

template<class charT,class traits>
errorstream<charT,traits>::errorstream(std::ostream &ostr)
  : basic_errorstream<charT,traits>(typename basic_errorstream<charT,traits>::sb_type_ptr(new errstream_type(ostr)))
{ }

template<class charT,class traits>
errorstream<charT,traits>::errorstream(const base_sb_type_ptr& sbuffer)
  : basic_errorstream<charT,traits>(sbuffer)
{
  DEBUG((std::cout << "errorstream non-init constructor(sb=" << (void*)basic_errorbuf<charT,traits>::rdbuf().get() << ")"<< std::endl;))
}

template<class charT,class traits>
errorstream<charT,traits>::errorstream(std::ostream&      ostr,
                                       const string_type& st)
  : basic_errorstream<charT,traits>(sb_type_ptr(new errstream_type(ostr,st)))
{ }

template<class charT,class traits>
errorstream<charT,traits>::errorstream(const errorstream<charT,traits>& o)
                         : basic_errorstream<charT,traits>(o)
{ }

template<class charT,class traits>
errorstream<charT,traits>::~errorstream()
{ }

//=============================================
// errormultistream
//=============================================

template<class charT,class traits>
errormultistream<charT,traits>::errormultistream()
  : errorstream<charT,traits>(typename errorstream<charT,traits>::base_sb_type_ptr(new multibuf_type()))
{ }

template<class charT,class traits>
errormultistream<charT,traits>::errormultistream(const errormultistream<charT,traits>& o)
  : errorstream<charT,traits>(o.rdbuf())
{ }

template<class charT,class traits>
errormultistream<charT,traits>::~errormultistream()
{ }

template<class charT,class traits> 
inline void
errormultistream<charT,traits>::insert(const basic_errorstream<charT,traits> &o)
{ static_cast<multibuf_type*>(basic_errorstream<charT,traits>::rdbuf().get())->insert(o); }

template<class charT,class traits>
inline void
errormultistream<charT,traits>::insert(std::ostream &o)
{ static_cast<multibuf_type*>(basic_errorbuf<charT,traits>::rdbuf().get())->insert(o); }

template<class charT,class traits>
inline void
errormultistream<charT,traits>::insert(std::ostream &o, const sb_error_info &e)
{ static_cast<multibuf_type*>(basic_errorbuf<charT,traits>::rdbuf().get())->insert(o,e); }

template<class charT,class traits>
inline void
errormultistream<charT,traits>::restrict(error_restriction r, error_priority p)
{ static_cast<multibuf_type*>(basic_errorstream<charT,traits>::rdbuf().get())->restrict(r,p); }

template<class charT,class traits>
inline void
errormultistream<charT,traits>::restrict_range(error_priority p1, error_priority p2)
{ static_cast<multibuf_type*>(basic_errorbuf<charT,traits>::rdbuf().get())->restrict_range(p1,p2); }

  
} // End namesapce SAGE
