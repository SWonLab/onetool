#include "sage/global/string_tokenizer.h"

namespace SAGE {

string_tokenizer::iterator::iterator(const string_tokenizer& s, size_t i)
{
  st = &s;
  pos = i;

  value.resize(0);
  first = last = false;
  last_delim = '\0';

  find_start( true );
  find_end();
}

string_tokenizer::iterator& string_tokenizer::iterator::operator++()
{
  value.resize(0);
  last_delim = '\0';

  find_start();
  find_end();
  return *this;
}

void string_tokenizer::iterator::find_start(bool first_field)
{
  const string &s=st->st;

  // Catch out of bounds state
  if(pos >= s.size())
    return;

  // Construct shorthand notation for the following logic
  bool l = st->skip_leading_delimiters();
  bool t = st->skip_trailing_delimiters();
  bool c = st->skip_consecutive_delimiters();

  // See state table below for justification
  // Intuitively, no delimiter eliding is required so we return.
  if( !l && !t && !c )
    return;

  // Scan string for last consecutive delimiter
  size_t i = pos;
  for(; i < s.size() && strchr(st->delim.c_str(), s[i]) != NULL; ++i);

  // If no characters after the current are delimiters, then return
  if(i==pos) return;

  bool f = first_field;
  bool e = (i == s.size());

  // The next section of code implements the rules for when delimiters are
  // elided based on a set of state variables.  The logic looks fairly
  // obscure, but is based on the following state table and the previously
  // defined 5 variables (f,e,l,t,c).  Some cases are non-trivial to derive
  // since they rely on related behaviors of the 'first' and 'last'
  // variables.
  //
  //        \  L| 0 0 0 0 1 1 1 1
  //         \ T| 0 0 1 1 1 1 0 0
  //        FE\C| 0 1 1 0 0 1 1 0
  //        ----|-----------------    Coverage of R:
  //        0 0 | R     R R     R       1) !T !C !F
  //        0 1 | R             R       2)  T !C !F !E
  //        1 1 | R R R R     R R       3) !L  F
  //        1 0 | R R R R               4)  F  E !L !T
  //        
  //        R     = return without action
  //        blank = skip consecutive delimiters         

  if( (!t && !c && !f) || (t && !c && !f && !e) || (!l && f) || (f && e && !l && !t))
    return;

  // Update the current position if in bounds
  if(i <= s.size())
    pos=i;

  // Invoke the trailing field flag if we do not elide trailing fields and
  // a delimiter character is at the end of the string.
  if( !t && e )
  {
    last = true;
    pos=i-1;
  }
}

void string_tokenizer::iterator::find_end() 
{
  const string &s = st->st;

  last_delim = '\0';

  if(pos >= s.size())
  {
    pos = s.size();
    first = last = true;
    return;
  }

  bool literal    = false;
  bool quoted     = false;
  bool last_quote = false;
  char c;

  const int bsize = 512;
  char buffer[512];
  int bcount = 0;
  size_t size = s.size();  

  while(pos < size)
  {
    c=s[pos++];

    if(    !quoted 
        && !(c == '"' && last_quote) 
        && !literal 
        && strchr(st->delim.c_str(),  c) != NULL )
    {
      last_delim = c;

      if(last)
        last_delim = '\0';

      // Invoke last field rule (last character is a delimiter)
      if(!last && pos == s.size())
      {
        --pos;
        last = true;
      }

      break;
    }

    if(bcount == bsize)
    {
      value.append(buffer, bcount);
      bcount = 0;
    }

    if(literal)                 
    { 
      literal = false;
      switch(c)
      {
        case 'n' : buffer[bcount++] = '\n'; break;
        case 't' : buffer[bcount++] = '\t'; break;
        default  : buffer[bcount++] = c;    break;
      }
      continue;
    }
    else if(c == '\\')
    {
      last_quote       = false;
      literal          = true;
    }
    else if(c == '"')
    {
      if(last_quote)
      { 
        last_quote = false; 
        quoted     = !quoted; 
        buffer[bcount++] = c; 
      }
      else if(quoted)
      {
        quoted     = false;
        last_quote = true;
      }
      else
      {
        // last_quote = true; // the first quote is never an internal quote
        quoted     = true;
      }
    }
    else 
    {
      last_quote       = false;
      buffer[bcount++] = c;
    }
  }

  if(literal)
    last_delim = '\\';

#if REPORT_ERROR
  if(quoted)
  // report error in file format
#endif

  if(bcount != 0)
    value.append(buffer, bcount);

  if(st->ws.size())
    value = strip_ws(value, st->ws.c_str());
  else
    value = strip_ws(value);
}

}

