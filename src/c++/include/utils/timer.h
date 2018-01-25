#pragma once
#ifndef __WISARD_TIMER_H__
#define __WISARD_TIMER_H__

#include <time.h>
#ifdef _WIN32
#	include <WinSock.h>
#endif
#include "global/common.h"

namespace ONETOOL {

class cTimer
{
	char		S_buf[256];
	timeval		X_tiStart, X_tiEnd;
	wsReal		D_tiElapsed;
public:
	cTimer(int B_start=0) { if (B_start) start(); }
	void		start();				/* 타이머를 시작한다 */
	wsReal		get();					/* 시작 시점으로부터 경과한 초를 반환 */
	const char*	getReadable();
	static const char*
				fmt(wsReal R_val);
};

} // End namespace ONETOOL

#endif

