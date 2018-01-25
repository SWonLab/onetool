#include "sage/global/LSF.h"

#ifdef NAMESPACE_STD
NAMESPACE_STD
#endif

class NameManager;
NameManager* LSFNameMgr;

long int LSFBase::existing_components = 0L;

#ifdef VERBOSE_MEM_DEBUG

#include <set>

typedef set< LSFBase *, less<LSFBase *> > LSFmemset;
static LSFmemset *LSFmemmap;

static void addMemMap(LSFBase *b)
{
  if( !LSFmemmap ) 
    LSFmemmap = new LSFmemset();
  LSFmemmap->insert( b );
}

static void removeMemMap(LSFBase *b)
{
  if(LSFmemmap) LSFmemmap->erase( b );
}

LSFBase *LSFBase::getLSFObjs()
{
  LSFBase *b = new LSFBase("memory list");
  LSFmemset::iterator i;
  for( i = LSFmemmap->begin(); i != LSFmemmap->end(); i++)
    if( (*i) != b )
      b->List(TRUE)->add( (*i) );
  return b;
}

bool LSFBase::isValidLSF(LSFBase *b)
{
  return ( LSFmemmap->find(b) != LSFmemmap->end() );
}

#else

LSFBase *LSFBase::getLSFObjs()      { return NULL; }
bool LSFBase::isValidLSF(LSFBase *) { return true; }

#endif // VERBOSE_MEM_DEBUG

LSFBase::~LSFBase()
{
  free_attrs();
  free_list();
  --existing_components;
#ifdef VERBOSE_MEM_DEBUG
  removeMemMap(this);
#endif
}

void LSFBase::init(const string &n, lsf_t t)
{ 
  _name = n;
  _type = LSF_BASE; _local_type = 0;
  if(t) type(t);
  _attrs = NULL; _list = NULL; 
  _ref = 1;
  _selfish = false; 
  _state = 0;
#ifdef VERBOSE_MEM_DEBUG
  addMemMap(this);
#endif
  ++existing_components; 
}

LSFBase::AList *LSFBase::attrs(bool create) 
{ if(create && !_attrs) _attrs = new AList; return _attrs; }

LSFList *LSFBase::List(bool create)
{ if(create && !_list) _list = new LSFList; return _list; }

void LSFBase::free_attrs() { if(_attrs) { delete _attrs; _attrs = NULL; } }
void LSFBase::free_list()  { if(_list)  { delete _list;  _list  = NULL; } }

void LSFBase::AcceptAll(LSFVisitor *v) 
{ 
  if(!_list) return;
  LSFList::iterator i = _list->begin();
  for( ; i != _list->end(); i++ )
    (*i)->Accept(v); 
}

LSFIterator::LSFIterator()                      { }
LSFIterator::LSFIterator(LSFBase *x)            { it = x->local_iterator(); }
LSFIterator::LSFIterator(const LSFIterator &x)  { it = x.it;                }
LSFIterator::LSFIterator(LSF_list_iterator *x)  { it = x;                   }

LSF_list_iterator::LSF_list_iterator(LSFBase *c) 
{ 
  l = NULL; 
  if(c && c->List())
  {
    l=c->List();
    i=c->List()->begin(); 
  }
}

LSF_list_iterator::LSF_list_iterator(LSFBase *c, iterator &x) 
{ 
  l = NULL;
  if(c && c->List() )
  {
    l = c->List(); 
    i = x; 
  }
}

void LSFVisitor::Process(LSFBase      *p) { def(p); }

bool LSFmap::set(const string &s, LSFBase *v, bool b)
{ 
  LSFBase* l = find(s);

  if (!b && !l) return false;
  if(l == v) return true;
  p_map[s] = v;  
  if(v == this)
    v->release();
  return true;
}

void LSFmap::remove(const string &s)
{ 
  LSFmap::iterator i=p_map.find(s);
  if( i == p_map.end() ) return;
  
  if( (*i).second == this) 
      (*i).second->grab();
  p_map.erase(i);
}

//
//
//
LSFFactory::LSFFactory(const char *n = NULL)
          : LSFBase(n)
{
  _type = LSF_FACTORY;
  add_factory(NULL, lsf_mappings);
}

LSFFactory::~LSFFactory()
{
  for(factorymap::iterator i=_factory_map.begin(); i!=_factory_map.end(); i++)
    if( (*i).second && (*i).second != this )
      (*i).second->release();
}

const string &LSFFactory::type_map( lsf_t t )
{
  static string nilString = "";
  // First check our local object mappings
  LSFFactory::typemap::iterator n = _type_map.begin();
  for( ; n != _type_map.end(); n++)
    if( (*n).second == t )
      return (*n).first;

  // If no valid binding exists then return nothing
  return nilString;
}

void LSFFactory::add_factory( LSFFactory *f, Type_mapping *m )
{
  if( !m ) return;
  for(int i=0; (int) m[i].value != -1; i++)
  {
    _type_map [ m[i].name ] = m[i].value; // add the name binding
    _factory_map[ m[i].value ] = f;
    if(f && f != this) f->grab();
  }
}

void LSFFactory::add_factory( LSFFactory *f, LSF_mapping *m )
{
  if( !m ) return;
  for(int i=0; m[i].name && (int) m[i].value != -1; i++)
  {
    _type_map [ m[i].name ] = m[i].value << LSF_SHIFT; // add the name binding
    _factory_map[ m[i].value << LSF_SHIFT ] = f;
    if(f && f != this) f->grab();
  }
}

void LSFFactory::add_factory( LSFFactory *f, lsf_t t, 
                                  Local_mapping *m )
{
  for(int i=0; (int) m[i].value != -1; i++)
  {
    _type_map [ m[i].name ] = t << LSF_SHIFT | m[i].value;
    _factory_map[ t<< LSF_SHIFT | m[i].value ] = f;
    if(f && f != this) f->grab();
  }
}

void LSFFactory::dump_map(ostream &o)
{
  LSFFactory::typemap::iterator n = _type_map.begin();
  for( ; n != _type_map.end(); n++)
    o << (*n).first << '\t' << (*n).second << endl;
}

LSFBase *LSFFactory::build( lsf_t t, const char *name, const char *type, 
                                AList *s, LSFBase *p)
{
  // Note: We don't check if we have the mapping since this is
  //       a generic and therefor _ABSTRACT_

  lsf_t  obtype  = t;

  if ( (long)obtype == -1 && type && *type ) 
  {
    LSFFactory::typemap::iterator n = _type_map.find(type);
    if( n != _type_map.end() )
      obtype = (*n).second;
  }

  // If the type is still unknown at this point then lets try a local
  // build and hope for the best.

  if( (long)obtype == -1 )
  {
    if((name && *name) || (type && *type) || (s && s->size()))
      return build_local( obtype, name, type, s, p);
    else
      return NULL;
  }

  // Find the factory which makes the object
  LSFFactory::factorymap::iterator f;
  f = _factory_map.find( obtype );

  if( f != _factory_map.end() )
  {
    if ( (*f).second )      // foreward the request if needed
      return (*f).second->build((*f).first, name, type, s, p);
    else                    // if not build locally
      return        build_local((*f).first, name, type, s, p);
  }
  return NULL;              // if all else fails, Give up!
}

LSFFactory* Factory;

LSFBase *LSFFactory::build_local( lsf_t t, const char* name,
                                  const char *, AList *l, LSFBase*)
{
  LSFBase *c = NULL;

  // Do storage type mapping
  if((int) t == -1 && l)
  {
    AList::iterator i = l->find(0);
    if( i != l->end() )
      switch( (*i).second.Type() )
      {
        case AttrVal::Attr_Integer: t = LSF_INT    << LSF_SHIFT; break;
        case AttrVal::Attr_Real:    t = LSF_REAL   << LSF_SHIFT; break;
        case AttrVal::Attr_String:  t = LSF_STRING << LSF_SHIFT; break;
        case AttrVal::Attr_Ptr:   
        case AttrVal::Attr_None:
        default:                  break;
      }
  }

  // Create the object
  switch (t>>LSF_SHIFT)
  {
    case LSF_MAP       :     c = new LSFmap ();       break;
    case LSF_REF       :   //c = new LSFref ();       break;
    case LSF_ITER      :
    case LSF_COMPOSITE :     
    case LSF_COMPONENT :     
    case LSF_INT       :
    case LSF_REAL      :
    case LSF_STRING    :
    case LSF_EXPR      :
    case LSF_SWITCH    :
    case LSF_GUARD     :     c = new LSFBase();
                             if(c) c->type(t);        break;
    case LSF_BASE      :     c = new LSFBase();       break;
    default            :     if(l || name) 
                               c = new LSFBase();     break;
  }

  if(!c) return NULL;

  if(name && *name)
    c->name( name );

  if(l)
    for(AList::iterator i = l->begin(); i != l->end(); i++)
      c->attrs(TRUE)->set( (*i) );

  return c;
}
