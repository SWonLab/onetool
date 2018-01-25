#pragma once
#ifndef __WISARD_ANALYSIS_H__
#define __WISARD_ANALYSIS_H__

#include <map>
#include <math.h>
#include "global/worker.h"
#include "global/io.h"
#include "utils/util.h"
#include "global/option.h"
using namespace std;

namespace ONETOOL {

int forVariant_equal(int N_mode, int N_thread, void *Vp_thrData,
	wsUint L_thrData, void *Vp_shareData, wsUint *Na_waitQueue);
int forAllSNP_ana(int N_mode, int N_thread, void *Vp_thrData,
	wsUint L_thrData, void *Vp_shareData, wsUint *Na_waitQueue);
int forAllSample_ana(int N_mode, int N_thread, void *Vp_thrData,
	wsUint L_thrData, void *Vp_shareData, wsUint *Na_waitQueue);
int forAllSample_anaXthr(int N_mode, int N_thread, void *Vp_thrData,
	wsUint L_thrData, void *Vp_shareData, wsUint *Na_waitQueue);

typedef struct _xAnaThread
{
	cIO		*Cp_IO;
	void	*Vp_data;
} xAnaThread;

typedef struct _xAnaPermThread
{
	cIO		*Cp_IO;
	void	*Vp_data;
} xAnaPermThread;

class cSetManagerAnalysis;
typedef struct _xAnaSetThread
{
	cIO*					Cp_IO;
	cSetManagerAnalysis*	Cp_gs;
	void*					Vp_data;
} xAnaSetThread;

class cAnalysis
{
protected:
	cIO	*Cp_IO;
public:
	cAnalysis(cIO *Cp_inpIO) {
		Cp_IO = Cp_inpIO;
	}
	virtual ~cAnalysis() {
	}
	cIO* getIO() const { return Cp_IO; }
	virtual void run() = 0;
};

} // End namespace ONETOOL

#endif
