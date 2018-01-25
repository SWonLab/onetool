#pragma once
#include <stdio.h>
#include <stdarg.h>
#ifdef USE_GZ
#	include "input/bgzf.h"
#endif

namespace ONETOOL {

typedef enum _xExportType {
	ET_AUTO,
	ET_PLAIN,
	ET_BIN,
	ET_GZIP,
	ET_BGZF
} xExportType;

class cIO;
struct xVariant;
class cExporter
{
protected:
	char			S_fn[512];
public:
	cExporter();
	virtual ~cExporter();
	static cExporter*	summon(const char *S_ext, char B_append=0,
		xExportType X_eType=ET_AUTO);
	virtual void*		open(const char *S_ext)=0;
	virtual void		fmt(const char *S_fmt, ...)=0;
	virtual void		put(const char *S_fmt)=0;
	virtual void		put(const char *S_fmt, wsUint N_times)=0;
	virtual void		write(const void *S_buf, wsUint N_len)=0;
	virtual void		sampleWise(cIO *Cp_IO, wsVecCst Ra_v)=0;
};

typedef enum _xTableType {
	TT_TABULAR,
	TT_XML,
	TT_JSON,
	TT_CSV,
} xTableType;

class cTableExporter
{
	xTableType	X_tt;
	char		S_fmt[512];
	cExporter*	Cp_exp;
	vStr		Xv_headers;
	int			N_elem;
	int			N_intPos;
	void		_init(const char *S_ext, const char *S_desc);
public:
	cTableExporter(const char *S_ext, vStr Xv_inpHeaders,
		const char *S_inpFmt, const char *S_desc=NULL, char B_variant=0);
	cTableExporter(const char *S_ext, const char *S_inpFmt,
		const char *S_desc, char B_variant, wsUint N_sz, ...);
	cTableExporter(const char *S_ext, const char *S_inpFmt,
		const char *S_desc, char B_variant, vStr& Sa_header);
	cTableExporter(const char *S_ext, const char** Sa_headers, wsUint N_sz,
		const char *S_inpFmt, const char *S_desc=NULL, char B_variant=0);
	~cTableExporter();
//	FILE*	getHandle() { return H_fp; }
	void	close();
	void	write(wsUint N_sz, ...);
	void	put(wsUint N_sz, ...);
	void	next();
	void	writeVariant(xVariant* Xp_vrt, ...);
	void	sampleWise(cIO *Cp_IO, wsVecCst Ra_v);
};

class cPlainExporter : public cExporter
{
	FILE *H_fp;
#ifdef _DEBUG
	size_t N_written;
#endif
public:
	cPlainExporter() {
		H_fp = NULL;
#ifdef _DEBUG
		N_written = 0;
#endif
	}
	cPlainExporter(const char *S_ext, char B_append=0, char B_bin=0);
	~cPlainExporter();
	FILE*	getHandle() { return H_fp; }
	void*	open(const char *S_ext);
	void	fmt(const char *S_fmt, ...);
	void	put(const char *S_fmt);
	void	put(const char *S_fmt, wsUint N_times);
	void	write(const void *S_buf, wsUint N_len);
	void	sampleWise(cIO *Cp_IO, wsVecCst Ra_v);
};

class cGzipExporter : public cExporter
{
	myGzHandler H_fp;
public:
	cGzipExporter() {
		H_fp = NULL;
	}
	cGzipExporter(const char *S_ext, char B_append=0, char B_bin=0);
	~cGzipExporter();
	myGzHandler	getHandle() { return H_fp; }
	void*		open(const char *S_ext);
	void		fmt(const char *S_fmt, ...);
	void		put(const char *S_fmt);
	void		put(const char *S_fmt, wsUint N_times);
	void		write(const void *S_buf, wsUint N_len) {
		myGzWrite(S_buf, sizeof(char), N_len, H_fp);
	}
	void		sampleWise(cIO *Cp_IO, wsVecCst Ra_v);
};

#ifdef USE_GZ
class cBgzfExporter : public cExporter
{
	BGZF* H_fp;
public:
	cBgzfExporter() {
		H_fp = NULL;
	}
	cBgzfExporter(const char *S_ext, char B_append=0, char B_bin=0);
	~cBgzfExporter();
	BGZF*	getHandle() { return H_fp; }
	void*	open(const char *S_ext);
	void	fmt(const char *S_fmt, ...);
	void	put(const char *S_fmt);
	void	put(const char *S_fmt, wsUint N_times);
	void	write(const void *S_buf, wsUint N_len) {
		bgzf_write(H_fp, S_buf, N_len);
	}
	void	sampleWise(cIO *Cp_IO, wsVecCst Ra_v);
};
#endif

} // End namespace ONETOOL
