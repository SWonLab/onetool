#include "global/app_log.h"
//#include "global/lib_include.h"
//#include "global/errormanip.h"

namespace ONETOOL {

string to_upper(const string &str)
{
  return SAGE::to_upper(str);
}

string to_lower(const string &str)
{
  return SAGE::to_lower(str);
}

string strip_ws(const string& s, const char *ws)
{
  return SAGE::strip_ws(s, ws);
}
/*
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

void app_log::init_output_streams()
{
  cerrorstream Screen = cerrorstream(cout);

  //   screen      -> Screen
  my_screen_stream.insert(Screen);

  if( !log_file )
  {
    infoqc_cerr << priority(fatal)
                << "Cannot open information log file.  Exiting..." << endl;
    exit(EXIT_FAILURE);
  }

  log_file.open((my_program_name + ".log").c_str());

  cerrorstream Info(log_file);

  //   info         -> Screen + Info (+ Debug)
  my_info_stream.insert(Screen);
  my_info_stream.insert(Info);

  // error stream: set the prefix for all error messages produced.
  // Each stream added to the error stream must be prefixed with this code.

  string prefix = "%%" + to_upper(my_program_name) + "-%p: ";

  //   errors      -> infoqc_cerr (if priority >= error)
  infoqc_cerr.prefix(prefix);

  my_error_stream.insert(infoqc_cerr);
  my_error_stream.restrict(r_ge, error);

  //   errors      -> Info (if priority >= information)
  Info.prefix(prefix);

  my_error_stream.insert(Info);
  my_error_stream.restrict(r_ge, information);

  // if debug is true we construct the debug stream and add it to the others
  if( my_debug )
  {
    //   create Debug stream
    dbg_file.open((my_program_name + ".dbg").c_str());

    cerrorstream Debug(dbg_file);

    //   add Debug to info and errors (priority >= debug)
    my_info_stream.insert(Debug);
    my_debug_stream.insert(Debug);

    Debug.prefix(prefix);

    my_error_stream.insert(Debug);
    my_error_stream.restrict(r_ge, debug);

    //   add Screen to errors (priority == debug)
    Screen.prefix(prefix);

    my_error_stream.insert(Screen);
    my_error_stream.restrict(r_eq, debug);
  }

  // Classify the non-error streams as raw output.
  my_screen_stream.set_raw_mode();
  my_info_stream.set_raw_mode();
}
*/
} // End namespace ONETOOL
