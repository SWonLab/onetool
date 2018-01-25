#include <stdio.h>
#include <math.h>
#include "utils/util.h"
#include "utils/timer.h"

namespace ONETOOL {

#if TOOLSET_TYPE == TOOLSET_ONETOOL
cTimer t;
#endif

void cTimer::start()
{
	gettimeofday(&X_tiStart, NULL);
}

wsReal cTimer::get()
{
	gettimeofday(&X_tiEnd, NULL);
	D_tiElapsed = (X_tiEnd.tv_sec - X_tiStart.tv_sec) * REAL_CONST(1000.0);      // sec to ms
	D_tiElapsed += (X_tiEnd.tv_usec - X_tiStart.tv_usec) / REAL_CONST(1000.0);   // us to ms

	return D_tiElapsed;
}

const char* cTimer::getReadable()
{
	wsReal		R_time = get();

	if (R_time < REAL_CONST(1000.0))
		sprintf(S_buf, "%gmsec.", R_time);
	else if (R_time < REAL_CONST(60000.0))
		sprintf(S_buf, "%gsec.", R_time/REAL_CONST(1000.0));
	else
		sprintf(S_buf, "%dmin.%gsec.", (int)floor(R_time/REAL_CONST(60000.0)),
			(R_time-floor(R_time/REAL_CONST(60000.0))*
				REAL_CONST(60000.0))/REAL_CONST(1000.0));

	return S_buf;
}

const char* cTimer::fmt(wsReal R_time)
{
	static char S_buf[64];

	if (R_time < 1000.0f)
		sprintf(S_buf, "%g msec.", R_time);
	else if (R_time < 60000.0f)
		sprintf(S_buf, "%g sec.", R_time/1000.0f);
	else
		sprintf(S_buf, "%d min. %g sec.", (int)floor(R_time/60000.0f),
			(R_time-floor(R_time/60000.0f)*60000.0f)/1000.0f);

	return S_buf;
}

} // End namespace ONETOOL
