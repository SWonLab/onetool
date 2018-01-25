//==========================================================================
// File:      app_base.cpp
//
// Author:    Yeunjoo E. Song
//
// History:   Initial implementation                            yes Nov 2014
//
// Notes:     This source file implements app_base.h.
//
// Copyright (c) 2014 Sungho Won
//   All Rights Reserved
//==========================================================================

#include "app/app_base.h"

using namespace std;

namespace ONETOOL {
namespace APP    {

#ifndef CXXFLAGS
  const char *app_base::cxxflags = NULL;
#else
  const char *app_base::cxxflags = CXXFLAGS;
#endif

#ifndef BUILD
  const char *app_base::build = NULL;
#else
  const char *app_base::build = BUILD;
#endif

#ifndef LDFLAGS
  const char *app_base::ldflags = NULL;
#else
  const char *app_base::ldflags = LDFLAGS;
#endif

#ifndef LDLIBS
  const char *app_base::ldlibs = NULL;
#else
  const char *app_base::ldlibs = LDLIBS;
#endif

#ifndef BUILD_DATE
  const char *app_base::build_date = "";
#else
  const char *app_base::build_date = BUILD_DATE;
#endif

string app_base::release_string;

string app_base::get_release_string() { return release_string; }

//================================================
//
//  CONSTRUCTOR
//
//================================================
app_base::app_base(std::string name, int argc, char **argv)
{
  free(malloc(1));

  my_name    = name;
  this->argc = argc;
  this->argv = argv;

  dbg = 0;

  calculate_release_string();

  print_title(cout);
}

//===============================================
//
//  calculateReleaseString()
//
//===============================================
void
app_base::calculate_release_string() const
{
  // Process current time:

  time_t cur_time = time(NULL);

  struct tm * now = localtime(&cur_time);

  std::ostringstream cur_time_str;

  cur_time_str << now->tm_mday        << " ";

       if(now->tm_mon ==  0) cur_time_str << "Jan";
  else if(now->tm_mon ==  1) cur_time_str << "Feb";
  else if(now->tm_mon ==  2) cur_time_str << "Mar";
  else if(now->tm_mon ==  3) cur_time_str << "Apr";
  else if(now->tm_mon ==  4) cur_time_str << "May";
  else if(now->tm_mon ==  5) cur_time_str << "Jun";
  else if(now->tm_mon ==  6) cur_time_str << "Jul";
  else if(now->tm_mon ==  7) cur_time_str << "Aug";
  else if(now->tm_mon ==  8) cur_time_str << "Sep";
  else if(now->tm_mon ==  9) cur_time_str << "Oct";
  else if(now->tm_mon == 10) cur_time_str << "Nov";
  else if(now->tm_mon == 11) cur_time_str << "Dec";

  cur_time_str << " "
               << 1900 + now->tm_year << " "
               << (now->tm_hour < 10 ? "0" : "") << now->tm_hour        << ":"
               << (now->tm_min  < 10 ? "0" : "") << now->tm_min         << ":"
               << (now->tm_sec  < 10 ? "0" : "") << now->tm_sec;

  // Generate release string:

  release_string = "[" + SAGE::toUpper(my_name) + " " + version.print_pretty() + "; " + build_date + "]" +
                   " -- " + cur_time_str.str() + "\n";
}

//===================================
//
//  getItselfPath(...)
//
//===================================
char getItselfPath(char *Sp_path)
{
#ifdef _WIN32
    Sp_path[0] = '\0';
    DWORD r = GetModuleFileName(NULL, Sp_path, MAX_PATH);
    /* Remove file itself */
    for (char *q = Sp_path + r - 1; q >= Sp_path && *q != '\\'; q--)
        *q = '\0';
    if (!Sp_path[0]) return 0;
#else
    char    S_tmp[PATH_MAX];
    pid_t   pid = getpid();
#   if defined(__FreeBSD__) || defined(__APPLE__)
    getcmdline(pid, S_tmp);
    strcpy(Sp_path, S_tmp);
#   elif defined(__sun)
    strcpy(Sp_path, getexecname());
#   else
    sprintf(S_tmp, "/proc/%d/exe", pid);
    if (readlink(S_tmp, Sp_path, PATH_MAX) == -1) {
        perror("readlink");
        return 0;
    }
    /* Remove file itself */
    for (char *q=Sp_path+strlen(Sp_path)-1 ; q>=Sp_path && *q != '/' ; q--)
        *q = '\0';
#   endif
#endif
    return 1;
}

//===================================
//
//  print_title(...)
//
//===================================
void app_base::print_title(ostream& o)
{
  o << "###############################################################################" << endl
    << "###############################################################################" << endl
    << "#" << endl
    << "#              " << get_release_string()
    << "#" << endl;
  o << "# Remember you have agreed to add an appropriate statement (including the" << endl
    << "# NIH grant number) under \"acknowledgments\" in any publication of results" << endl
    << "# obtained by using this program. Suggested wording is:" << endl
    << "#" << endl
    << "# \"(Some of)The results of this paper were obtained by using the software" << endl
    << "# package " << SAGE::toUpper(my_name) << ", which was supported by the National Research Foundation" << endl
    << "# of Korea Grant funded by Korean Government (NRF-2014S1A2A2028559).\"" << endl
    << "#" << endl;
  o << "###############################################################################" << endl
    << "###############################################################################" << endl << endl;

  return;
}
 
//===================================
//
//  print_exe_info(...)
//
//===================================
void app_base::print_exe_info(ostream& o)
{
  char S_exePath[512] = { 0, };
  if (S_exePath[0] == '\0' && getItselfPath(S_exePath) == 0)
      return;

  o << "exe path: " << S_exePath << endl;
  o << "exe command: " << argv[0];
  for( int i = 1 ; i < argc ; i++ )
    o << " " << argv[i];
  o << endl << endl;

  return;
}

//===================================
//
//  print_debug(...)
//
//===================================
void app_base::print_debug(ostream& o)
{
  o << "Application Information: ";

  if( !build_date && !cxxflags && !ldflags && !ldlibs && !build )
  {
    cout << "none. " << endl;
    return;
  }
  cout << endl;

                                       o << "    Version: " << version.print_pretty() << endl;
  if(build_date && strlen(build_date)) o << "      Built: " << build_date << endl;
  if(build      && strlen(build))      o << "      BUILD: " << build << endl;
  if(cxxflags   && strlen(cxxflags))   o << "   CXXFLAGS: " << cxxflags << endl;
  if(ldflags    && strlen(ldflags))    o << "    LDFLAGS: " << ldflags << endl;
  if(ldlibs     && strlen(ldlibs))     o << "     LDLIBS: " << ldlibs << endl;

  o << endl
    << "-------------------------------------------------------------------------------"
    << endl << endl;

  return;
}

//===================================
//
//  print_help(...)
//
//===================================
void app_base::print_help(ostream& o)
{
  o //<< " Usage 1: " << my_name << " --script script_file_name" << endl
    //<< " Usage 2: " << my_name << " --vfile file_root_name [options]" << endl
    //<< " Usage 3: " << my_name << " --bfile file_root_name [options]" << endl
    << endl
    << " Usage 1: " << my_name << " --vcf vcf_file_name --fam fam_file_name [options]" << endl
    << " Usage 2: " << my_name << " --bed bed_file_name --bim bim_file_name --fam fam_file_name [options]" << endl;
  o << endl
    << "Command line options for input/output:" << endl;
  //o << "\nSpecify the input files with the same root name with the expected file extensions" << endl
  //  << "  --bfile   PLINK binary ped files, e.g. test.bed, test.bim and test.fam" << endl;
  o << "\nRequired input files:" << endl
    << "  --vcf     VCF file name" << endl
    << "  --fam     PLINK family file name" << endl
    << "  or" << endl
    << "  --bed     PLINK binary ped file name" << endl
    << "  --bim     PLINK bim file name" << endl
    << "  --fam     PLINK family file name" << endl
    << "  or" << endl
    << "  --ped     PLINK ped file name" << endl
    << "  --map     PLINK map file name" << endl;
  o << "\n[Optional input files]" << endl
    << "  --pheno   Phenotype/Covariate file name" << endl
    << "  --typ     Type probability file from SEGREG to run LODLINK analysis" << endl
    << "  --map     Map file to run MERLIN analysis" << endl;
  o << "\n[option for output file]" << endl
    << "  --out     root name for all output files, 'onetool' by default" << endl;
  o << "\n\n[more options - please refer to the manual]" << endl;
                                                                      
  o << endl
    << "-------------------------------------------------------------------------------"
    << endl << endl;
  return;
}

} // End namespace APP
} // End namespace ONETOOL
