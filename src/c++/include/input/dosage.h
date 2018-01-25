#pragma once
#include "global/io.h"

namespace ONETOOL {

typedef enum _xDosageType {
	DSG_UNKNOWN,
	DSG_BEAGLE,
	DSG_MINIMAC,
	DSG_MERLIN,
	DSG_IMPUTE2,
	DSG_RAW,
	DSG_OTHER
} xDosageType;

class cDosageIO : public cIO
{
	xDosageType	X_type;
	void		_exportDosageDist(wsReal R_minDsg, wsReal R_maxDsg);
	wsStrCst		getFormat() { return "io.dosage"; }
	void		_init(wsStrCst S_fn);
public:
	cDosageIO(wsStrCst S_fn, char B_inDry);
	~cDosageIO();
};

} // End namespace ONETOOL
