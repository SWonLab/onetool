#pragma once
#ifndef __WISARD_WORKER_H__
#define __WISARD_WORKER_H__

#ifdef _WIN32
#	include <windows.h>
#endif
#include "global/common.h"
#include "utils/util.h"

#ifdef _WIN32
	#define myThreadRet DWORD WINAPI
	#define myThread HANDLE
	#define myEvent HANDLE
	#define myEventSet SetEvent
	#define myEventReset ResetEvent
	#define myEventWait(a) WaitForSingleObject(a, INFINITE)
	#define myEventCreate(a, n) a = CreateEvent(NULL, TRUE, FALSE, n)
#else
	#define myThreadRet void*
	#include <semaphore.h>
	#include <unistd.h>
	#include <signal.h>
	#include <errno.h>
	#define myThread pthread_t
	#define myEventReset(a)
#	ifdef __APPLE__
#		define myEvent sem_t*
#		define myEventWait(a) while (sem_trywait(a) == -1)
#		define myEventSet(a) if (sem_post(a) == -1) halt("SYSERR: SEMPOST failed [%d]", errno);
#		define myEventCreate(a, n) { char name[64]; \
			struct timeval m; \
			gettimeofday(&m, NULL); \
			sprintf(name, "wisard%d.%d", getpid(), m.tv_usec); \
			sem_t *ss = sem_open(name, O_CREAT, 0777, 0); \
			if (ss == SEM_FAILED) halt("SYSERR: SEMOPEN[%s] failed [%d]", name, errno); a = ss; }
#	else
#		define myEvent sem_t
#		define myEventReset(a)
#		define myEventSet(a) if (sem_post(&a) == -1) halt("SYSERR: SEMPOST failed [%d]", errno)
#		define myEventWait(a) while (sem_trywait(&a) == -1)
#		define myEventCreate(a, n) sem_init(&a, 0, 0)
#	endif
#endif

#define DISTFUN_INIT		0
#define DISTFUN_DISTRIBUTE	1
#define DISTFUN_UNINIT		2
#define DISTFUN_CHECKEND	3
#define DISTFUN_DISTASYNC	4
#define DISTFUN_AFTERLOOP	5

namespace ONETOOL {

typedef struct _xThread {
	wsUint		N_size;
} xThread;

typedef struct _xRangedThread {
	wsUint		N_start;
	wsUint		N_end;
} xRangedThread;

typedef int (*WORKERPROC)(int N_idx, void *Vp_shareData, void *Vp_data,
	void *Vp_result);
typedef int (*VARIANTPROC)(int N_idx, void *Vp_shareData, void *Vp_data,
	void *Vp_result);
typedef int (*DISTRIBUTEPROC)(int N_mode, int N_thread, void *Vp_thrData,
	wsUint L_thrData, void *Vp_shareData, wsUint *Na_waitQueue);
wsUint	getCoreCount();
int		distEqual(int N_mode, int N_thread, void *Vp_thrData,
	wsUint L_thrData, void *Vp_shareData, wsUint *Na_waitQueue);

class cWorker;
cWorker& WORKER();

struct gData2
{
	cWorker *pWorker;
	int runnable;
	int tid;
	WORKERPROC func;
	void *data, *sdata, *res;
#ifdef USE_CUDA
	cudaDeviceProp X_propGPU;
#endif
	gData2(cWorker* _pWorker, int _tid) : pWorker(_pWorker), tid(_tid), func(NULL) {}

	void set(WORKERPROC _func, void *_shareData, void *_data, void *_res) {
		func	= _func;
		sdata	= _shareData;
		data	= _data;
		res		= _res;
	}
};

#ifdef _WIN32
DWORD WINAPI thrS(void *Vp);
#endif

class cWorker
{
	pthread_mutex_t		H_workerSync;

	/* Indicates the number of threads used in parallelization */
	wsUint		*Na_waitQueue;
	int			N_thread;
	myThread	*Ha_thrHandle;
	gData2		**Xa_thrData;
public:
	myEvent		*Ha_evtHandle;
	myEvent		*Ha_evtHandleToHost;

	cWorker();
	~cWorker();
	void	clear();
	int		init(int N_thread);
	//WORKER().run(calcMAF, forAllSNP, this, Na_missGeno, sizeof(MYuint))
	void	run(WORKERPROC thrFunc, DISTRIBUTEPROC distFunc,
		void *shareData, void *resArray, wsUint N_szData=sizeof(int));
	void	runAsync(WORKERPROC thrFunc, DISTRIBUTEPROC distFunc,
		void *shareData, void *resArray, wsUint N_szData=sizeof(int));
	void	runVariant(VARIANTPROC thrFunc, void *shareData, void *resArray,
		wsUint N_szData=sizeof(int));
	const gData2*
			getThrData(wsUint N_thrIdx) { return Xa_thrData[N_thrIdx]; }
	void	setWait(wsUint N_thrIdx);
};

} // End namespace ONETOOL

#endif
