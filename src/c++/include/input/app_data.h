#ifndef __APP_DATA_H
#define __APP_DATA_H
//==========================================================================
// File:      app_data.h
//
// Author:    Yeunjoo E. Song
//
// History:   Initial implementation                            yes Jan 2015
//
// Notes:     This header file defines the input data for infoqc.
//
// Copyright (c) 2015 Sungho Won
//   All Rights Reserved
//==========================================================================

#include "sage/rped/rped.h"
#include "input/vcf.h"

namespace ONETOOL {

class app_data
{
  public:

    app_data(string name, bool debug);

    ~app_data();

    cIO* input_wisard();
    void input_sage(cIO* vio);

          SAGE::RPED::RefMultiPedigree& pedigrees();
    const SAGE::RPED::RefMultiPedigree& pedigrees() const;

    SAGE::Output_Streams& get_ostreams() const;
    SAGE::cerrorstream&   errors() const;
    ostream&              info() const;
    ostream&              screen() const;
    ostream&              messages() const;

  protected:
  
    bool read_pedigree_file(cIO* vio, char* fname);
    bool read_vcf_file(cIO* vio, char* fname);
    bool read_pheno_file(cIO* vio, char* fname);
    bool read_locus_description_file(cIO* vio, char* fname);

    string                                my_program_name;
    SAGE::RPED::RefMultiPedigree          my_pedigrees;

    mutable SAGE::Output_Streams          my_output;
};

} // End namespace ONETOOL

#include "input/app_data.ipp"

#endif
