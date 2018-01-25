#include "sage/global/LSFtypes.h"

LSF_mapping lsf_mappings[] =
{ { LSF_NONE,      (char *)"None"      },  // LSF objects: formal names
  { LSF_BASE,      (char *)"base"      },
  { LSF_COMPONENT, (char *)"item"      },
  { LSF_COMPOSITE, (char *)"list"      },
  { LSF_REF,       (char *)"ref"       },
  { LSF_MAP,       (char *)"map"       },
  { LSF_FACTORY,   (char *)"factory"   },
  { LSF_INT,       (char *)"int"       },
  { LSF_REAL,      (char *)"real"      },
  { LSF_STRING,    (char *)"string"    },  // LSF storage types
  { LSF_ITER,      (char *)"iterator"  },
  { LSF_GUARD,     (char *)"guard"     },  // LSF interpreted types
  { LSF_EXPR,      (char *)"expr"      },
  { LSF_EXPR,      (char *)"expression"},
  { LSF_SWITCH,    (char *)"switch"    },
  { (unsigned long) -1,     (char *)0  }  };

