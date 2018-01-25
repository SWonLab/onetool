#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "utils/util.h"
#include "global/worker.h"

namespace ONETOOL {

/**< The one and only instance of cWorker */
cWorker C_worker;

/**
 * Getter function to access to worker instance
 *
 * @param Sp_str		A pointer that points to the beginning of string
 * @param Sp_nextStr	A pointer that be pointed to the beginning of
 * @return				The reference of one and only instance of cWorker
 */
cWorker& WORKER()
{
	return C_worker;
}

myThreadRet workerThread(void *lpdwThreadParam)
{
	gData2	&X_workerData = *((gData2 *)lpdwThreadParam);
	cWorker	*pWorker = X_workerData.pWorker;

#ifdef USE_CUDA
	/* Initialize CUDA device */
	cutilSafeCall(cudaSetDevice(X_workerData.tid));
	/* Fetch device info */
	cutilSafeCall(cudaGetDeviceProperties(&(X_workerData.X_propGPU),
		X_workerData.tid));
#endif

	//	printf("[%d] init\n", GetCurrentThread());
	while (pWorker && pWorker->Ha_evtHandle) {
		/* Wait until event is awaken */
		//fprintf(stderr, "Waiting thread %d...\n", X_workerData.tid );
		myEventWait(pWorker->Ha_evtHandle[X_workerData.tid]);
		//fprintf(stderr, "Go thread %d...\n", X_workerData.tid);

		/* Immediately terminate if do not work */
		if (X_workerData.runnable == 0) {
			myEventSet(pWorker->Ha_evtHandleToHost[X_workerData.tid]);
			myEventReset(pWorker->Ha_evtHandle[X_workerData.tid]);
			break;
		}

		/* Function check */
		if (X_workerData.func == NULL) {
			LOG("SYSERR: Tried to execute NULL function");
			exit(1);
		}


		/* Execute */
		int N_ret = X_workerData.func(X_workerData.tid,
			X_workerData.sdata, X_workerData.data, X_workerData.res);

		if (!pWorker->Ha_evtHandleToHost)
			goto _halt;
		myEventSet(pWorker->Ha_evtHandleToHost[X_workerData.tid]);
		myEventReset(pWorker->Ha_evtHandle[X_workerData.tid]);

		/* Return value check */
		switch (N_ret) {
		case -9:
			/* Immediate exit */
			goto _halt;
			break;
		case 0:
			/* Normal exit */
			break;
		default:
			/* Invalid exit */
			LOG("Invalid return value\n");
			exit(1);
		}
	}
_halt:
	return 0;
}

cWorker::cWorker() {
	Ha_thrHandle		= NULL;
	Xa_thrData			= NULL;
	Ha_evtHandle		= NULL;
	Ha_evtHandleToHost	= NULL;
	pthread_mutex_init(H_workerSync);
}

cWorker::~cWorker()
{
//	clear();
}

int do_terminate(int N_idx, void *Vp_shareData, void *Vp_data,
	void *Vp_result)
{
	return -9;
}

int dist_terminate(int N_mode, int N_thread, void *Vp_thrData,
	wsUint L_thrData,  void *Vp_shareData, wsUint *Na_waitQueue)
{
	static char B_term = 0;
	char B_retThread = 0;

	switch (N_mode) {
	case DISTFUN_INIT:
		break;
	case DISTFUN_DISTRIBUTE:
		if (B_term == 0)
			B_retThread = 1;
		B_term = 1;
		break;
	case DISTFUN_UNINIT:
		break;
	case DISTFUN_AFTERLOOP:
		break;
	default:
		halt("Unsupported distFunc command `%d`", N_mode);
		break;
	}

	return B_retThread?N_thread:0;
}

void cWorker::clear()
{
	//fprintf(stderr, "AA");
	run(do_terminate, dist_terminate, NULL, NULL);
	//fprintf(stderr, "BB");
	DEALLOC(Ha_thrHandle);
	for (int i=0 ; i<N_thread ; i++) {
		delete Xa_thrData[i];
	}
	DEALLOC(Xa_thrData);
	DEALLOC(Ha_evtHandle);
	DEALLOC(Ha_evtHandleToHost);
	DEALLOC(Na_waitQueue);

	N_thread = 0;
}

int cWorker::init(int N_inpThread)
{
	timeval X_uniq;
	gettimeofday(&X_uniq, NULL);

	this->N_thread		= N_inpThread;
	wsAlloc(Ha_evtHandle, myEvent, N_thread);
	wsAlloc(Ha_evtHandleToHost, myEvent, N_thread);
	wsAlloc(Ha_thrHandle, myThread, N_thread);
	wsAlloc(Xa_thrData, gData2*, N_thread);
	wsAlloc(Na_waitQueue, wsUint, N_thread);

	int N_core = (int)getCoreCount();
//	LOG("Maximum number of available cores is [%d]\n", N_core);
	if (N_thread > N_core) {
		LOGnote("Number of threads adjusted [%d] to [%d] for performance\n",
			N_thread, N_core);
		N_thread = N_core;
	}
	LOG("Initialized with [%d] out of maximum [%d] threads\n", N_thread, N_core);
	report(WISARD_REP_TOTNCORE, N_core);
	report(WISARD_REP_USEDNCORE, N_core);

	for (int i=0 ; i<N_thread ; i++)  {
		char S_evtName[32];
		Xa_thrData[i]	= new gData2(this, i);
		Xa_thrData[i]->runnable = 1;

		sprintf(S_evtName, "evt%d_%d", i, (int)(X_uniq.tv_usec));
		myEventCreate(Ha_evtHandle[i], S_evtName);
		sprintf(S_evtName, "evtH%d_%d", i, (int)(X_uniq.tv_usec));
		myEventCreate(Ha_evtHandleToHost[i], S_evtName);

#ifdef _WIN32
		Ha_thrHandle[i] = CreateThread(NULL, NULL, workerThread,
			Xa_thrData[i], NULL, NULL);
#else
		pthread_create(&Ha_thrHandle[i], NULL, workerThread, Xa_thrData[i]);
#endif
	}

	return N_thread;
}

#ifdef _WIN32
DWORD WINAPI thrS(void *Vp)
{
	void **Vpp = (void **)Vp;
	myEvent *Hp = (myEvent *)Vpp[0];
	unsigned int N = *((unsigned int *)Vpp[1]);

	for (unsigned int i=0 ; i<N ; i++)
		myEventSet(Hp[i]);
	return 0;
}
#endif

wsUint getCoreCount()
{
#ifdef _WIN32
	SYSTEM_INFO sysinfo;
	GetSystemInfo(&sysinfo);

	return (wsUint)(sysinfo.dwNumberOfProcessors);
#elif defined (CTL_HW)
	int mib[4];
	wsUint numCPU;
	size_t len = sizeof(numCPU); 

	/* set the mib for hw.ncpu */
	mib[0] = CTL_HW;
	mib[1] = HW_AVAILCPU;  // alternatively, try HW_NCPU;

	/* get the number of CPUs from the system */
	sysctl(mib, 2, &numCPU, &len, NULL, 0);

	if( numCPU < 1 ) 
	{
		mib[1] = HW_NCPU;
		sysctl( mib, 2, &numCPU, &len, NULL, 0 );

		if( numCPU < 1 )
		{
			numCPU = 1;
		}
	}
	return numCPU;
#else
	return sysconf( _SC_NPROCESSORS_ONLN );
#endif
}

int distEqual(int N_mode, int N_thread, void *Vp_thrData, wsUint L_thrData,
	void *Vp_shareData, wsUint *Na_waitQueue)
{
	unsigned char*
			Np_curThrData	= (unsigned char *)Vp_thrData;
	wsUint	N_szTarget		= ((xThread *)Vp_shareData)->N_size - 1;
	UINT_s	N_proc			= 0;
	int		N_ret			= N_proc != 0 ? 0 : N_thread;
	int		i=0;

	pverbose("mode [%d] thread [%d]\n", N_mode, N_thread);

	switch (N_mode) {
	case DISTFUN_INIT:
		i = 0;
		N_proc = 0;
		break;
	case DISTFUN_DISTRIBUTE:
		memset(Np_curThrData, 0x00, L_thrData*N_thread);
		for ( ; i<N_thread ; i++,Np_curThrData+=L_thrData) {
			xRangedThread* Xp_rt = (xRangedThread *)Np_curThrData;
			Xp_rt->N_start	= i;
			Xp_rt->N_end	= i + (wsUint)((wsReal)(N_szTarget - i) / (wsReal)N_thread) * N_thread;
		}
		N_proc = N_thread;
		break;
	case DISTFUN_UNINIT:
		LOG("%d/%d samples processed\n", N_szTarget, N_szTarget);
		break;
	case DISTFUN_AFTERLOOP:
		break;
	default:
		halt_fmt(WISARD_SYST_INVL_DISTFUNC_CMD, N_mode);
		//		halt("Unsupported distFunc command [%d]", N_mode);
		break;
	}
	return N_ret;
}

void cWorker::run(WORKERPROC thrFunc, DISTRIBUTEPROC distFunc,
	void *Vp_shareData, void *resArray, wsUint L_thrData)
{
	static char *Va_thrData = NULL;

	/* Init distFunc */
	distFunc(DISTFUN_INIT, N_thread, NULL, 0, Vp_shareData, NULL);

	if (Va_thrData)
		halt("SYSERR: Va_data should be NULL, maybe it called"
			" 'WITHIN' thread?");

	/* While distFunc is active */
	Va_thrData = (char *)malloc(L_thrData * N_thread);
	int N_allowedThr;
	while ((N_allowedThr = distFunc(DISTFUN_DISTRIBUTE, N_thread,
		Va_thrData, L_thrData, Vp_shareData, NULL)) != 0) {
		/* Wake up threads */
		for (int i=0 ; i<N_allowedThr ; i++)
			Xa_thrData[i]->set(thrFunc, Vp_shareData, Va_thrData + (L_thrData*i),
				resArray);
		for (int i=0 ; i<N_allowedThr ; i++)
			myEventReset(Ha_evtHandleToHost[i]);

		//fprintf(stderr, "Q1 %d", N_allowedThr);
#ifdef _WIN32
		void *Vp[] = { Ha_evtHandle, &N_allowedThr };
		HANDLE h = CreateThread(NULL, 0, thrS, Vp, 0, NULL);
#else
		for (int i=0 ; i<N_allowedThr ; i++) {
			myEventSet(Ha_evtHandle[i]);
		}
#endif
		//fprintf(stderr, "Q2");

		/* Wait */
#ifdef _WIN32
		//			Sleep(1);
		//			fprintf(stderr, "SigWait\n");
		WaitForMultipleObjects(N_allowedThr, Ha_evtHandleToHost, TRUE,
			INFINITE);
		CloseHandle(h);
#else
		for (int i=0 ; i<N_allowedThr ; i++) {
			//fprintf(stderr, "Release %d\n", i);
			myEventWait(Ha_evtHandleToHost[i]);
		}
#endif
		//fprintf(stderr, "X");
		distFunc(DISTFUN_AFTERLOOP, N_allowedThr, Va_thrData, L_thrData,
			Vp_shareData, NULL);
	}
	free(Va_thrData);
	Va_thrData = NULL;
	//fprintf(stderr, "QU");
	/* Uninit distFunc */
	distFunc(DISTFUN_UNINIT, N_thread, NULL, 0, Vp_shareData, NULL);
	//fprintf(stderr, "XU");
}

void cWorker::runAsync(WORKERPROC thrFunc, DISTRIBUTEPROC distFunc,
	void *Vp_shareData, void *resArray, wsUint L_thrData)
{
	static char *Va_thrData = NULL;

	/* Init distFunc */
	distFunc(DISTFUN_INIT, N_thread, NULL, 0, Vp_shareData, NULL);

	if (Va_thrData)
		halt("SYSERR: Va_data should be NULL, maybe it called"
		" 'WITHIN' thread?");

	/* While distFunc is active */
	/* Fill up wait queue */
	Va_thrData = (char *)malloc(L_thrData * N_thread);
for (int i=0 ; i<N_thread ; i++)
		Na_waitQueue[i] = 1;

	/* Until this run is complete */
	while (distFunc(DISTFUN_CHECKEND, N_thread, NULL, 0, NULL, NULL)) {
		/* This area should protect Na_distIdx */

		/* Give work to waiting threads */
		pthread_mutex_lock(&H_workerSync);
		distFunc(DISTFUN_DISTASYNC, N_thread, Va_thrData, L_thrData,
			Vp_shareData, Na_waitQueue);

		/* Wake up threads */
		for (int i=0 ; i<N_thread ; i++) {
			if (Na_waitQueue[i])
				Xa_thrData[i]->set(thrFunc, Vp_shareData, Va_thrData + (L_thrData*i),
					resArray);
		}
		for (int i=0 ; i<N_thread ; i++) {
			if (Na_waitQueue[i])
				myEventReset(Ha_evtHandleToHost[i]);
		}
		for (int i=0 ; i<N_thread ; i++)
			if (Na_waitQueue[i])
				myEventSet(Ha_evtHandle[i]);
		pthread_mutex_unlock(&H_workerSync);

		/* Waiting for update, using first only */
#ifdef _WIN32
		WaitForSingleObject(Ha_evtHandleToHost[0], INFINITE);
#else
		myEventWait(Ha_evtHandleToHost[0]);
#endif
	}
	free(Va_thrData);
	Va_thrData = NULL;

	/* Uninit distFunc */
	distFunc(DISTFUN_UNINIT, N_thread, NULL, 0, Vp_shareData, NULL);
}

void cWorker::setWait( wsUint N_thrIdx )
{
	pthread_mutex_lock(&H_workerSync);
	Na_waitQueue[N_thrIdx] = 1;
	pthread_mutex_unlock(&H_workerSync);
}

} // End namespace ONETOOL
