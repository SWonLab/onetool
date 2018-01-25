#ifndef SUBMODEL_H
#define SUBMODEL_H

#include "sage/maxfun/ParameterInput.h"

namespace SAGE   {
namespace MAXFUN {

class Submodel
{
  public:

    Submodel(cerrorstream& errors = sage_cerr);
    Submodel(const Submodel& other);
    Submodel&  operator=(const Submodel& other);
  
    virtual ~Submodel();

    bool isLinked() const;

  protected:

    friend class ParameterMgr;

    virtual int update() = 0;
    virtual int finalizeConfiguration();
    virtual void setAdvancedParameterOptions();

    double& getParam(int local_id);
    double  getParam(int local_id) const;

    int addParametersToParameterMgr();

    void linkToParameterMgr(ParameterMgr* mi);
    void unlinkFromParameterMgr();
    
          ParameterMgr* getParameterMgr();
    const ParameterMgr* getParameterMgr() const;

    vector<ParameterInput> my_parameters;

    mutable cerrorstream  my_errors;

  private:

    ParameterMgr* my_info;
};

} // End of namespace MAXFUN
} // End of namespace SAGE

#include "sage/maxfun/Submodel.ipp"

#endif
