#pragma once
#include "global/common.h"
#include "utils/util.h"

#define STR_OPEN_FAILED	1

namespace ONETOOL {

class cStream
{
protected:
	wsStrCst	S_desc;
	char*	S_inpName;
	char*	S_path;
	void	*H_stream;
	char	B_gzip;
	char	B_fail;
	char	B_binary;
public:
	cStream();
	virtual ~cStream();
	wsStrCst	getPath() { return S_path; }
	/* Open a stream,
	 * 0 if no error, otherwise error */
	virtual char	open(wsStrCst S_path, wsStrCst S_desc, char B_dontHalt,
		char B_binary)=0;
	/* Close a stream,
	 * 0 if no error, otherwise error */
	virtual char	close()=0;
	/* 1 if end, 0 otherwise */
	virtual char	end()=0;
	/* Read a line from given stream */
	virtual char*	gets(char *S_buf, wsUint N_len)=0;
	/* Tell current offset of given stream */
	virtual size_t	tell()=0;
	/* Seek specific position */
	virtual void	seek(size_t N_pos, int X_type)=0;
	/* Read byte */
	virtual size_t	read(void *Bp, wsUint N_sz)=0;
	virtual size_t	read(unsigned char *Bp, wsUint N_sz=1)=0;
	virtual size_t	read(char *Bp, wsUint N_sz=1)=0;
	virtual size_t	read(unsigned int *Bp, wsUint N_sz=1)=0;
	virtual size_t	read(double *Bp, wsUint N_sz=1)=0;
	/* Rewind given stream
	 * 0 if no error, otherwise error */
	virtual char	rewind()=0;
	/* Get the name */
	virtual wsStrCst	getName() {
		return S_inpName;
	}
};

class cStrFile : public cStream
{
public:
	cStrFile();
	cStrFile(FILE* H_fp);
	cStrFile(wsStrCst S_inpPath, wsStrCst S_inpDesc=NULL, char B_dontHalt=0,
		char B_inpBinary=0);
	~cStrFile();

	char	isFailed() { return B_fail; }
	virtual char	open(wsStrCst S_inpPath, wsStrCst S_inpDesc=NULL, char B_dontHalt=0,
		char B_inpBinary=0);
	char	close();
	char	end();
	char*	gets(char *S_buf, wsUint N_len);
	char*	ugets(char B_draw=0);
	char	rewind();
	size_t	tell();
	int		isEnd();
	void	seek(size_t N_pos, int X_type);
	size_t	read(float *Bp, wsUint N_sz);
	size_t	read(void *Bp, wsUint N_sz);
	size_t	read(unsigned char *Bp, wsUint N_sz=1);
	size_t	read(char *Bp, wsUint N_sz=1);
	size_t	read(unsigned int *Bp, wsUint N_sz=1);
	size_t	read(double *Bp, wsUint N_sz=1);
};

class cNetFile : public cStream
{
protected:
	void	_init();
	char*	_getHost(wsStrCst S_url, wsUint *Np_port=NULL);
public:
	cNetFile();
	cNetFile(wsStrCst S_url, wsStrCst S_inpDesc=NULL);
	~cNetFile();

	char open(wsStrCst S_url, wsStrCst S_inpDesc=NULL, char B_dontHalt=0,
		char B_inpBinary=0);
	char close();
};

} // End namespace ONETOOL
