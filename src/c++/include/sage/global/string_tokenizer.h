#ifndef __STRING_TOK_H
#define __STRING_TOK_H

//==========================================================================
// File:      string_tokenizer.h
//
// Author:    S.A.G.E. crew
//
// History:   From S.A.G.E. to oneTool                          yes May 2015
//
// Notes:     string_tokenizer class
//            - provides const forward iteration through delimited strings.
//
// Copyright (c) 2014 Sungho Won
//   All Rights Reserved
//==========================================================================

#include "sage/global/output_streams.h"

namespace SAGE {

class string_tokenizer
{
  public:

    class iterator;
    friend class iterator;

    typedef iterator const_iterator;

    string_tokenizer(const string& s = "", const string &d = ",", 
                                           const string &w =" \t\n\r")
      : st(s), delim(d), ws(w), skip_consecutive_delim(false),
        skip_leading_delim(false), skip_trailing_delim(false) { }

    ~string_tokenizer() { }
    
    // Iterators  

    iterator begin() const { return iterator(*this, 0);           }
    iterator   end() const { return iterator(*this, st.size()+1); }

    const string &str()        const { return st;    }
    const string &whitespace() const { return ws;    }
    const string &delimiters() const { return delim; }

    void set_str(const string &s)        { st = s;    }
    void set_whitespace(const string &w) { ws = w;    }
    void set_delimiters(const string &d) { delim = d; }

    bool skip_consecutive_delimiters() const
    { 
      return skip_consecutive_delim;
    }
    bool skip_leading_delimiters() const
    { 
      return skip_leading_delim;
    }
    bool skip_trailing_delimiters() const
    { 
      return skip_trailing_delim;
    }
    void set_skip_consecutive_delimiters(bool skip = true) 
    { 
      skip_consecutive_delim = skip;
    }
    void set_skip_leading_delimiters(bool skip = true) 
    { 
      skip_leading_delim = skip;
    }
    void set_skip_trailing_delimiters(bool skip = true) 
    { 
      skip_trailing_delim = skip;
    }
    
  protected:

    string st, delim, ws;
    bool skip_consecutive_delim;
    bool skip_leading_delim;
    bool skip_trailing_delim;

  public:

    class iterator
    {
      public:

        iterator() : st(NULL), pos(0), first(false), last(false), last_delim('\0') { }
        iterator(const string_tokenizer&, size_t);
        iterator(const iterator& s)
           : st(s.st), pos(s.pos), first(s.first), last(s.last), 
             last_delim(s.last_delim), value(s.value)              { }
        
        ~iterator() { st = NULL; }
        
        iterator& operator= (const iterator& rhs)
        { 
          if(&rhs == this) return *this;
        
          st         = rhs.st; 
          value      = rhs.value; 
          first      = rhs.first;
          last       = rhs.last;
          pos        = rhs.pos;
          last_delim = rhs.last_delim;
          return *this; 
        }
        
        char last_delimiter() const { return last_delim; }

        // Dereference

        const string &operator*() const  { return value;  }
        const string *operator->() const { return &value; }

        // Comparisons

        bool operator==(const iterator& rhs) const
        { return st == rhs.st && pos == rhs.pos && first == rhs.first 
            && last == rhs.last; }

        bool operator!= (const iterator& rhs) const
        { return !(*this == rhs); }
        
        // iteration

        iterator& operator++();
        iterator  operator++(int) { iterator i(*this); ++(*this); return i; }

        iterator& operator--();
        iterator  operator--(int) { iterator i(*this); --(*this); return i; }
        
      protected:

        void find_end();
        void find_start(bool first_field = false);
        
        const string_tokenizer* st;
        size_t pos;
        bool first;
        bool last;
        char last_delim;
        string value;
    };
};

}

#endif
