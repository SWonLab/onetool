//==========================================================================
// File:      string_util.cpp
//
// Author:    Yeunjoo E. Song
//
// Notes:     Implements the string utility function in string_util.h
//
// Copyright (c) 2014 Sungho Won
//   All Rights Reserved
//==========================================================================

#include "sage/global/string_util.h"

#ifdef _WIN32
#	ifdef max
#		undef max
#		undef min
#	endif
#	define finite std::isfinite
#endif

namespace SAGE{
/*
int kill_ws(FILE *i, const char *ws, const char *eol)
{
  if(!i || ferror(i) || feof(i)) return 0;

  int l = 0;
  unsigned int c;

  for(c=getc(i); !ferror(i) && !feof(i); c=getc(i))
  {
    if( strchr(ws, c) ) continue;
    else if( strchr(eol,c) ) { ++l; continue; }
    else break;
  }

  if(!feof(i))
    ungetc(c, i);

  return l;
}

string getString(FILE *i, const char *delim, const char *ws)
{
  const int bsize = 512;
  char buffer[512];
  int bcount = 0;
  string s;

  if( !i || feof(i) || ferror(i) )
    return s;

  char        c;
  const char* sel = NULL;

  for( c=getc(i); !ferror(i) && !feof(i); c=getc(i) )
  {
    if( (sel=strchr(delim,  c)) != NULL ||
             strchr(ws,     c)  != NULL )
      break;

    if( bcount == bsize )
    {
      s.append(buffer, bcount);
      bcount = 0;
    }

    buffer[bcount++] = c;
  }

  if( bcount != 0 )
    s.append(buffer, bcount);

  if( sel )
    ungetc(c,i);

  return s;
}
*/
int kill_ws(istream &i, const char *ws, const char *eol)
{
  if(!i || i.eof()) return 0;
  int l = 0;
  char c;
  for(c=i.get(); i; c=i.get())
  {
    if( strchr(ws, c) ) continue;
    else if( strchr(eol,c) ) { ++l; continue; }
    else break;
  }

  if(!i.eof())
    i.putback(c);

  return l;
}

string getStringSAGE(istream &i, const char *delim, const char *ws )
{
  char buffer[100];
  char *stop  = buffer;
  string s = "";

  if(!i || i.eof())
    return s;

  char        c;
  const char* sel = NULL;
  //while(i.get(c))
  for(c=i.get(); i; c=i.get() )
  {
    if( (sel=strchr(delim,  c)) != NULL ||
             strchr(ws,     c)  != NULL)
       break;
    *stop++ = c;
    if(stop == buffer + 98)
    {
      *stop = '\0';
      s+=buffer;
      stop = buffer;
    }
  }

  if(stop != buffer)
  {
    *stop = '\0';
    s+=buffer;
  }

  if(sel)
    i.putback(c);

  return s;
}

string getQStringSAGE(istream &i, const char *delim, const char *ws, const char *eol, int *lines)
{
  char buffer[100];
  char *stop  = buffer;
  string s;

  if(!i || i.eof()) return s;

  bool literal = false;
  bool quoted  = false;
  char        c;
  const char* sel = NULL;

  for(c=i.get(); i; c=i.get())
  {
    sel = NULL;
    if( !quoted && !literal && ( (sel=strchr(delim,  c)) != NULL) )
      break;

    if(stop == buffer + 98)
    {
      *stop = '\0';
      s+=buffer;
      stop = buffer;
    }

    if(literal)
    {
      literal = false;
      switch(c)
      {
        case 'n' : *stop++ = '\n'; break;
        case 't' : *stop++ = '\t'; break;
        default: if(!strchr(eol, c)) *stop++ = c;
                 else if( lines ) (*lines)++;
                    break;
      }
      continue;
    }
    else if(c == '\\')            literal = true;
    else if(c == '\"' && quoted)  quoted  = false;
    else if(c == '\"')            quoted  = true;
    else if(strchr(eol, c))       break;
    else if(!quoted && strchr( ws, c )) /*skip white space*/ ;
    else *stop++ = c;
  }

  if(stop != buffer)
  {
    *stop = '\0';
    s+=buffer;
  }

  if(sel)
    i.putback(c);

  return s;
}

int printQuoted(ostream &o, string& s, int pos, int mark, int m)
{
  if(!o) return pos;

  int len = s.length();
  int stop  = len;

  while (m < len)
  {
    if( m && pos != mark && len-m+pos > 75 && len-m > 1)
    {
      o << endl;
      for(pos=0; pos<mark; ++pos) o << ' ';
    }

    stop = len;
    if( pos + len + 2 - m > 76 )
    {
      stop = m+76-pos;
      if(stop > len)
        stop = len;
      for(int p = 5+m; p < len && pos + len - p + 2 < 77; ++p)
        if(strchr(" \t\n;.,*-+/|?!", s[p]))
          stop = p;
      if(stop <= m) stop = m+1;
    }

    o << '\"';
    for(unsigned int i=0; i < s.size(); ++i)
    {
      switch( s[i] )
      {
        case '\r' : break;
        case '\n' : o << "\\n";  break;
        case '\t' : o << "\\t";  break;
        case '\\' : o << "\\\\"; break;
        case '\"' : o << "\\\""; break;
        default   : o <<  s[i];  break;
      }
    }
    o << '\"';
    pos += 2+stop-m;
    m = stop;

    if (m < len)
    {
      o << " \\";
      pos += 2;
    }
  }

  return pos;
}

string to_upper(const string &str)
{
  string s(str);
  std::transform(s.begin(), s.end(), s.begin(), ::toupper);

  return s;
}

string to_lower(const string &str)
{
  string s(str);
  std::transform(s.begin(), s.end(), s.begin(), ::tolower);

  return s;
}

string strip_ws(const string& s, const char *ws)
{
  size_t b, e;

  if( s.length() == 0 ) return string();

  if( ws )
  {
    b = s.find_first_not_of(ws);
    e = s.find_last_not_of(ws);
  }
  else
  {
    for( b=0; b < s.length() && isspace(s[b]); ++b );
    for( e=s.length()-1; isspace(s[e]); --e );
  }

  if( e < b || b == (size_t)-1 ) return string();

  return s.substr(b,e-b+1);
}

string long2str (long l, int width, long flags, char fch)
{
  //return std::to_string(l);

  std::ostringstream o;

  if(flags) o << std::setiosflags( (std::ios_base::fmtflags)flags );
  if(width) o << std::setw(width);
  if(fch)   o << std::setfill(fch);

  o << l;

  return o.str();
}

string ptr2str (const void* v, int width, long flags, char fch)
{
  ostringstream o;

  if(flags) o << std::setiosflags( (std::ios_base::fmtflags)flags);
  if(width) o << std::setw(width);
  if(fch)   o << std::setfill(fch);

  o << v;

  return o.str();
}


string toUpper(const string &s) { return to_upper(s); }
string toLower(const string &s) { return to_lower(s); }

string doub2str(double d, int width, int pres, long flags, char fch)
{
/*
  if(width < 0) width = 0;

  bool b = (pres < 0);

  if(pres < 0 && width) pres = width;

  while(true)
  {
    std::ostringstream o;

    if(flags) o << std::setiosflags( (std::ios_base::fmtflags)flags );
    if(width) o << std::setw(width);
    if(fch)   o << std::setfill(fch);

    o << d;
    int size = o.str().size();
    if(b && width &&  size > width)
    {
      if(pres) pres -= (size - width);
      else
        return o.str();
    }
    else
      return o.str();

    if(pres < 0) pres = 0;
  }
*/
  std::ostringstream o;
  o << std::fixed << std::setw( width ) << std::setprecision( pres )
    << std::setfill( fch ) << d;

  return o.str();
}

double str2doub(const string& val)
{
  char *end_ptr = NULL;
  const char *begin = val.c_str();

  double v = strtod(begin, &end_ptr);
  if( end_ptr == begin )   // Invalid Input -- cannot create anything
    return SAGE::QNAN;

  return v;
}

long str2long(const string& val)
{
  char *end_ptr = NULL;
  const char *begin = val.c_str();

  long v = strtol(begin, &end_ptr, 10);
  if( end_ptr == begin )   // Invalid Input -- cannot create anything
    return 0;

  return v;
}

string fp(double d, size_t w, size_t p, char in)
{
  string invalid(in, 45);
  if( !finite(d) )
    return invalid.substr(0,w);

  string n = doub2str(d, w, p, std::ios::showpoint | std::ios::fixed, ' ');
//  if( n.size() > w )
//    n = doub2str(d, w, p, std::ios::scientific | std::ios::showpoint);

  if( n.size() > w )
    n = n.substr(0,w);

  return n;
}

string pval(double p, size_t w, int prec, int max_stars)
{
  // NOTE: Work in progress.  Should justify output to max_stars in all
  //       cases.  Right now it works only when max_stars == 2.

  if(prec < 0)
    prec = w - 3 - 2;

  string pv = fp(p,w-3,prec).substr(0,w);
  pv += " ";

  if( !finite(p) )
  {
    for( int i = 0; i < max_stars; ++i )
      pv += " ";
  }
  else if( p < 0.01 )
  {
    int stars = max_stars;

    if( p > 0 )
      stars = std::min( (int)-log10(p), max_stars );

    for( int i = 0; i < stars; ++i )
      pv += "*";

    if( stars < max_stars )
    {
      for( int i = stars; i < max_stars; ++i )
        pv += " ";
    }
  }
  else if( p < 0.05 )
  {
    pv += "*";

    for( int i = 1; i < max_stars; ++i )
      pv += " ";
  }
  else
  {
    for( int i = 0; i < max_stars; ++i )
      pv += " ";
  }

  return pv;
}

string fp_scientific(double d, size_t w, size_t p, char in)
{
  string invalid(in, 45);

  if( !finite(d) )
    return invalid.substr(0,w);

  string n = doub2str(d, w, p, std::ios::scientific | std::ios::showpoint, ' ');

  return n;
}

string pval_scientific(double p, size_t w, int prec, int max_stars)
{
  // NOTE: Work in progress.  Should justify output to max_stars in all
  //       cases.  Right now it works only when max_stars == 2.

  if(prec < 0)
    prec = w - 3 - 2;

  string pv = fp_scientific(p,w-3,prec).substr(0,w);

  pv += " ";

  if( !finite(p) )
  {
    for( int i = 0; i < max_stars; ++i )
      pv += " ";
  }
  else if( p < 0.01 )
  {
    int stars = max_stars;

    if( p > 0 )
      stars = std::min( (int)-log10(p), max_stars );

    for( int i = 0; i < stars; ++i )
      pv += "*";

    if( stars < max_stars )
    {
      for( int i = stars; i < max_stars; ++i )
        pv += " ";
    }
  }
  else if( p < 0.05 )
  {
    pv += "*";

    for( int i = 1; i < max_stars; ++i )
      pv += " ";
  }
  else
  {
    for( int i = 0; i < max_stars; ++i )
      pv += " ";
  }

  return pv;
}

}
