#ifndef PARAMETERMGR_H
#define PARAMETERMGR_H

#include "sage/maxfun/TransformationSubmodel.h"

namespace SAGE   {
namespace MAXFUN {

typedef std::shared_ptr<const ParamCalculator> ParamCalculatorShCstPtr;
typedef std::shared_ptr<ParamCalculator>       ParamCalculatorShPtr;

//===============================================================================
// class ParamCalculator, AdditiveParamCalculator : public ParamCalculator
//===============================================================================

class ParamCalculator
{
  public:
    virtual ~ParamCalculator() { }
    virtual double calculateParam (const ParameterMgr * mgr) const = 0;
};

class AdditiveParamCalculator : public ParamCalculator
{
  public:
    AdditiveParamCalculator(const vector<int> & ids);
    virtual double calculateParam (const ParameterMgr * mgr) const;
  private:
    vector<int> my_ids;
};

//===============================================================================
// class Transformer, SimpleTransformer : public Transformer
//===============================================================================

class Transformer
{
  public:

    // Required to make compiler happy.
    virtual ~Transformer() { }

    Transformer            () {};
    Transformer            (const Transformer &);
    Transformer& operator= (const Transformer &);

    virtual double transformToMaximizationScale (const ParameterMgr * mgr, double val) const = 0;
    virtual double transformToReportedScale    (const ParameterMgr * mgr, double val) const = 0;

    virtual double getLowerBound() const { return -MF_INFINITY; }
    virtual double getUpperBound() const { return MF_INFINITY; }
};

class SimpleTransformer : public Transformer
{
  public:
    virtual double transformToMaximizationScale (const ParameterMgr * mgr, double val) const { return val; }
    virtual double transformToReportedScale     (const ParameterMgr * mgr, double val) const { return val; }
};

typedef std::shared_ptr<const Transformer> TransformerShCstPtr;
typedef std::shared_ptr<Transformer> TransformerShPtr;

struct TransformerInfo
{
  TransformerInfo();
  TransformerInfo(const TransformerInfo & other);
  TransformerInfo& operator=(const TransformerInfo & other);

  int reported_param_id;
  int maximization_param_id;
  TransformerShCstPtr transformer;
};

inline TransformerInfo::TransformerInfo() {}

inline TransformerInfo::TransformerInfo(const TransformerInfo & other)
{
  reported_param_id = other.reported_param_id; maximization_param_id = other.maximization_param_id; transformer = other.transformer;
}

inline TransformerInfo& TransformerInfo::operator=(const TransformerInfo & other)
{
  reported_param_id = other.reported_param_id; maximization_param_id = other.maximization_param_id; transformer = other.transformer; return *this;
}

//===============================================================================
// class NewSubmodel stuff 
//===============================================================================

template<class SM>
struct SMType
{
  typedef SM sm_type;
};

typedef std::shared_ptr<NewSubmodel> NewSubmodelShPtr;
typedef std::shared_ptr<const NewSubmodel> NewSubmodelShCstPtr;

class NewSubmodel
{
  friend class ParameterMgr;

  public:

    inline NewSubmodel(const string & name, cerrorstream & errors = sage_cerr);
    inline NewSubmodel(const NewSubmodel & other);
    inline NewSubmodel& operator=(const NewSubmodel & other);
    virtual inline ~NewSubmodel() { }

    virtual int update() = 0;
    virtual NewSubmodelShPtr clone() = 0;
    virtual inline int finalizeConfiguration();
    virtual inline void setAdvancedParameterOptions();
    inline void setName(const string & name);
    inline const string & getName() const;
    inline const ParameterMgr & getParameterMgr() const;
    inline bool isFinalized() const;
    inline double & getParam(int local_id);
    inline double getParam(int local_id) const;

  protected:

    vector<ParameterInput> my_parameters;
    mutable cerrorstream  my_errors;

  private:

    string        my_name;
    bool          my_finalized;
    ParameterMgr* my_mgr;

    inline void setParameterMgr(ParameterMgr * mgr);
    inline int addParametersToMgr();
    inline void setFinalized(bool finalized);
};    

//===============================================================================
// class ParameterMgr
//===============================================================================

class ParameterMgr
{
  public:

    friend class APIMaxFunction;
    friend class Function;
    friend class OutputFormatter;
    friend class Results;
    friend class Submodel;
    friend class Maximizer;

    MEM_FRIEND(MAXFUN::ParameterMgr);

    ParameterMgr();
    ParameterMgr(const ParameterMgr & other);
    ParameterMgr& operator= (const ParameterMgr & other);

    ~ParameterMgr();

    void reset();

    int addParameter(
      string                   group_name,
      string                   param_name, 
      Parameter::ParamTypeEnum initial_type     =  MAXFUN::Parameter::INDEPENDENT, 
      double                   initial_estimate =  0.0, 
      double                   lower_bound      = -MAXFUN::MF_INFINITY, 
      double                   upper_bound      =  MAXFUN::MF_INFINITY,
      double                   init_stepsize    =  0.1);
      
    Parameter&  addParameterAlt(
      string                   group_name,
      string                   param_name, 
      Parameter::ParamTypeEnum initial_type     =  MAXFUN::Parameter::INDEPENDENT, 
      double                   initial_estimate =  0.0, 
      double                   lower_bound      = -MAXFUN::MF_INFINITY, 
      double                   upper_bound      =  MAXFUN::MF_INFINITY,
      double                   init_stepsize    =  0.1);      

    int addParameter(ParameterInput& param);

    Parameter& addParameterAlt(ParameterInput& param);    

    int addParamCalculator(int idx, ParamCalculatorShPtr calculator);

    int addTransformer(TransformerShCstPtr transformer, 
                       const string& maximization_group_name, 
                       const string& maximization_param_name,
                       const string& reported_group_name,
                       const string& reported_param_name);
                
    bool isInUse() const;

    void addGroup(string group_name) const;
    //void dumpConfiguration() const;

    int addSubModel(Submodel* sub_mod);

    void removeSubmodels();

    bool hasLinkedSubmodels() const;

    template<class SMType> std::shared_ptr<typename SMType::sm_type> createSubmodel(SMType t);

    NewSubmodelShPtr    getSubmodel(const string& name);
    NewSubmodelShCstPtr getSubmodel(const string& name) const;

    template<class SM> const typename SM::sm_type& getSubmodel(const string & name, SM t) const;
    template<class SM>       typename SM::sm_type& getSubmodel(const string & name, SM t);

    bool doesParamExist(string group_name, string param_name) const;

    int getParamID(string group_name, string param_name) const;
    int getParamID(string group_name, int group_id) const;

    int getParamCount() const;
    int getParamCount(string group_name) const;

    int getEstimatedParamCount() const;
                
    int getGroupCount() const;

    double getEst(int param_id) const;

    double& operator() (string group_name, string param_name);
    double  operator() (string group_name, string param_name) const;
    double& operator() (int param_id);
    double  operator() (int param_id) const;

          Parameter& getParameter (string group_name, string param_name);
    const Parameter& getParameter (string group_name, string param_name) const;

          Parameter& getParameter (int param_id);
    const Parameter& getParameter (int param_id) const;

          Parameter& getParameter(string group_name, int group_id);
    const Parameter& getParameter(string group_name, int group_id) const;

    ParameterIterator      getParamBegin();
    ParameterConstIterator getParamBegin() const;

    ParameterIterator      getParamEnd();
    ParameterConstIterator getParamEnd() const;

    ParameterIterator      getParamBegin(string group_name);
    ParameterConstIterator getParamBegin(string group_name) const;

    ParameterIterator      getParamEnd(string group_name);
    ParameterConstIterator getParamEnd (string group_name) const;
    
    // Added 8-29-7. djb
    //
    void  dumpParameterLookupTable() const;    

    // Added JA July -09 
    // Are all parameters fixed ?
    bool allfixed();

  protected:

    int finalizeSubmodels();

    void copy(const ParameterMgr &);
    void copyInUseParameterMgr(const ParameterMgr &);

    int update(vector<double> & params);

    void transformParameters();

    void calculateInternalInitialEstimates();

    int updateDependents(vector<double> & params);

    //void dump(DebugCfg & debug) const;

    vector<string> getGroupNames() const;

  public:

    vector<string> getOrderedGroupNames() const;

  protected:

    bool doesGroupExist(string group_name) const;

          param_group& getGroup(string group_name);
    const param_group& getGroup(string group_name) const;
    
    void removeSubmodel(Submodel* sm);

  private:

    // Master list of parameters:
    param_vector params;

    // Lookup table for parameter names:
    std::map<pair<string, string>, int> param_lookup_table;

    // Sub-lists of parameters:
    mutable std::map<string, param_group> groups;

    // Group ordering:
    mutable vector<string> group_ordering;

    // List of Submodels
    list<Submodel*> my_submodels;

    // List of NewSubmodels
    vector<NewSubmodelShPtr> my_new_submodels;

    // Parameter transformation
    vector<TransformerInfo> my_transformer_infos;

    // Parameter calculators:
    std::map<int, ParamCalculatorShPtr> my_param_calcs;
};

//================================================================
//  NewSubmodel INLINES
//================================================================

inline NewSubmodel::NewSubmodel(const string & name, cerrorstream & errors)
{ 
  my_finalized  = false;
  my_errors     = errors;
  my_name       = name;  
  my_parameters . clear();
}
 
inline NewSubmodel::NewSubmodel(const NewSubmodel & other)
{
  my_finalized  = other.my_finalized;
  my_name       = other.my_name;
  my_errors     = other.my_errors;
  my_parameters = other.my_parameters;
}

inline NewSubmodel&
NewSubmodel::operator=(const NewSubmodel & other)
{
  my_finalized  = other.my_finalized;
  my_name       = other.my_name;
  my_errors     = other.my_errors;
  my_parameters = other.my_parameters;

  return *this;
}

inline void NewSubmodel::setName(const string & name) { my_name = name; }

inline const string& NewSubmodel::getName() const { return my_name; }

inline bool NewSubmodel::isFinalized() const { return my_finalized; }
inline void NewSubmodel::setFinalized(bool finalized) { my_finalized = finalized; }

inline double&  NewSubmodel::getParam(int local_id)       { return my_mgr->operator()(my_parameters[local_id].index); }
inline double   NewSubmodel::getParam(int local_id) const { return my_mgr->operator()(my_parameters[local_id].index); }

inline void NewSubmodel::setParameterMgr(ParameterMgr* mgr) { my_mgr = mgr; }

inline int
NewSubmodel::addParametersToMgr()
{
  for(size_t i = 0; i < my_parameters.size(); i++)
    my_mgr->addParameter(my_parameters[i]);

  return 0;
}

inline int
NewSubmodel::finalizeConfiguration()
{
  return 0;
}

inline void
NewSubmodel::setAdvancedParameterOptions()
{
}

inline const ParameterMgr &
NewSubmodel::getParameterMgr() const
{
  return *my_mgr;
}


class FooSubmodel : public NewSubmodel
{
  friend class ParameterMgr;

  public:

    FooSubmodel() : NewSubmodel("Foo")
    {
      my_parameters.push_back(ParameterInput("Foo", "tempy", Parameter::FIXED, 1.0, -1.0, 1.0));
    }

    FooSubmodel(const FooSubmodel & other) : NewSubmodel(other) {}

    virtual int update() { return 0; }

    virtual NewSubmodelShPtr clone() { return NewSubmodelShPtr(new FooSubmodel(*this)); }

};

//=======================================================================
//  hasLinkedSubmodels()
//=======================================================================

inline bool
ParameterMgr::hasLinkedSubmodels() const
{
  return !my_submodels.empty();
}

//=======================================================================
//  isInUse()
//=======================================================================
inline bool
ParameterMgr::isInUse() const
{
  return hasLinkedSubmodels();
}

//=======================================================================
//  getParamID(...) #1
//=======================================================================
inline int
ParameterMgr::getParamID(string group_name, string param_name) const
{
  std::map<pair<string, string>, int>::const_iterator i = param_lookup_table.find(make_pair(group_name, param_name));

  if(i == param_lookup_table.end())
  {
    cout << "Error: Parameter '" << param_name << "' in group '" << group_name << "' does not exist." << endl;
    exit(1);
  }

  return i->second;
}

//=======================================================================
//  getParamID(...) #2
//=======================================================================
inline int
ParameterMgr::getParamID(string group_name, int group_id) const
{
  if(!doesGroupExist(group_name))
  {
    cout << "Error: Group " << group_name << " does not exist!" << endl;
    exit(1);
  }

  const param_group & g = getGroup(group_name);

  if((size_t)group_id >= g.size())
  {
    cout << "Error: Group id " << group_id << " too large for group " << group_name << endl;
    exit(1);
  }

  return g[group_id];
}

//=====================================================================
//  createSubmodel(...)
//=====================================================================

template<class SMType> 
inline std::shared_ptr<typename SMType::sm_type> 
ParameterMgr::createSubmodel(SMType t) 
{
  std::shared_ptr<typename SMType::sm_type> sm_sh_ptr = std::shared_ptr<typename SMType::sm_type>(new typename SMType::sm_type);

  sm_sh_ptr->setParameterMgr(this);

  my_new_submodels.push_back(sm_sh_ptr);

  return sm_sh_ptr;
}

//=====================================================================
//  getSubmodel(...) TEMPLATIZED CONST
//=====================================================================
template<class SM> 
inline const typename SM::sm_type& 
ParameterMgr::getSubmodel(const string & name, SM t) const
{
  return *(static_cast<const typename SM::sm_type *> (getSubmodel(name).get()));
}
                
//=====================================================================
//  getSubmodel(...) TEMPLATIZED NON-CONST
//=====================================================================
template<class SM> 
inline typename SM::sm_type& 
ParameterMgr::getSubmodel(const string & name, SM t)
{
  return *(static_cast<typename SM::sm_type *> (getSubmodel(name).get()));
}

} // End of namespace MAXFUN

MEM_COUNT_BEGIN(MAXFUN::ParameterMgr)
{
  size_t x = 0;
  
  for(size_t i = 0; i < t.params.size(); ++i)
    x += get_mem(t.params[i]);
    
  for(std::map<pair<string, string>, int>::const_iterator i = t.param_lookup_table.begin(); i != t.param_lookup_table.end(); ++i)
  {
    x += get_mem(i->first.first) + get_mem(i->first.second) + get_mem(i->second);
  }
 
  for(std::map<string, MAXFUN::param_group>::const_iterator i = t.groups.begin(); i != t.groups.end(); ++i)
  {
    x += get_mem(i->first) + get_mem(i->second);
  }

  x += get_mem(t.group_ordering);

  x += sizeof(MAXFUN::ParameterMgr);

  return x;
}
MEM_COUNT_END

} // End of namespace SAGE

#endif
