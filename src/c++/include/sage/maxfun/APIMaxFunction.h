#ifndef APIMAXFUNCTION_H
#define APIMAXFUNCTION_H

#include "sage/maxfun/Results.h"

namespace SAGE   {
namespace MAXFUN {

class APIMaxFunction : public MaxFunction
{
  public:

    APIMaxFunction(MaxFunction &, ParameterMgr &, const DebugCfg &);
    
    APIMaxFunction(const APIMaxFunction &);

    virtual double evaluate      (vector<double> & params);
    virtual int    update_bounds (vector<double> & params);

    ParameterMgr & getParameterMgr() const;

    // Tells this object to which OUTPUT::Table it should direct runtime output.
    //void setTable(OUTPUT::Table *);
    
  private:

    APIMaxFunction & operator= (const APIMaxFunction &);

    MaxFunction&        my_max_function;
    ParameterMgr&       my_parameter_mgr;
    const DebugCfg&     my_debug_cfg;
    //OUTPUT::Table*      my_table;
    double lastvalue; // due to JA
};

}} // End namespace

#endif
