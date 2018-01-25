#include <stdarg.h>
#include "input/stream.h"
#include "global/option.h"

#ifdef USE_NET
#	ifdef _WIN32
#		ifndef VLD_RPTHOOK_INSTALL /* Resolve conflict with VLD */
#			include <WinSock.h>
#		endif
#		pragma comment(lib, "ws2_32.lib")
#	else
#		include <sys/socket.h>
#		include <netinet/in.h>
#		include <arpa/inet.h>
#		define SOCKET		int 
#		define SOCKADDR_IN	sockaddr_in
#		define SOCKADDR		sockaddr 
#	endif
#endif

#ifndef _WIN32
#	include <unistd.h>
#	include <sys/types.h>
#	include <pwd.h>
#endif

#ifndef USE_GZ
#	define	gzclose	fclose
#	define	gzeof	feof
#	define	gzopen	fopen
#	define	gztell	ftell
#	define	gzseek	fseek
#	define	gzgets(a,b,c)	fgets(b,c,a)
#	define	gzread(a,b,c)	fread(b,c,a)
#	define	gzrewind	::rewind
#	define	gzgetc	fgetc
#	define	z_off_t	size_t
typedef FILE* gzFile;
#endif

namespace ONETOOL {

cStream::cStream()
{
	H_stream	= NULL;
}

cStream::~cStream()
{
}

cStrFile::cStrFile() : cStream()
{
	/* In default, set to fail */
	H_stream	= NULL;
	S_path		= NULL;
	S_inpName	= NULL;
	B_fail = 1;
}

cStrFile::cStrFile(FILE* H_fp)
{
	H_stream	= (void *)H_fp;
	S_path		= NULL;
	B_fail		= 0;
	B_binary	= 0;
}

cStrFile::cStrFile(wsStrCst S_inpPath, wsStrCst S_inpDesc/*=NULL*/,
	char B_dontHalt/*=0*/, char B_binary/*=0*/) : cStream()
{
	H_stream	= NULL;
	S_path		= NULL;
	B_fail = open(S_inpPath, S_inpDesc, B_dontHalt, B_binary);
}

cStrFile::~cStrFile()
{
	if (S_path)		delete [] S_path;
	if (S_inpName)	delete [] S_inpName;
	if (!H_stream)	return;

	B_gzip ? gzclose((gzFile)H_stream) : fclose((FILE *)H_stream);
}

char cStrFile::open(wsStrCst S_inpPath, wsStrCst S_inpDesc/*=NULL*/,
	char B_dontHalt/*=0*/, char B_inpBinary/*=0*/)
{
	if (!S_inpPath) halt("SYSERR : Empty path value");
	close();
#ifndef _WIN32
	/* Process home directory */
	if (S_inpPath[0] == '~') {
		struct passwd *pw = NULL;
		int k = 0;
		if (S_inpPath[1] == '/') {
			pw = getpwuid(getuid());
			k = 1;
		} else {
			/* Find first '/' */
			char *Sp = strchr((char *)S_inpPath, '/');
			if (Sp) {
				char S_nametmp[256] = { 0, };
				memcpy(S_nametmp, S_inpPath+1, Sp-S_inpPath-1);
				pw = getpwnam((char *)S_nametmp);
				k = Sp-S_inpPath;
			} else {
				pw = getpwnam(S_inpPath+1);
				k = strlen(S_inpPath);
			}
			/* Get home dir of that user */
			if (pw == NULL) {
				S_path = new char[1];
				S_path[0] = '\0';
				if (B_dontHalt)
					return STR_OPEN_FAILED;
				else
					S_desc ? halt_fmt(WISARD_FAIL_OPEN_FILEWITHDESC, S_desc, S_path) :
					halt_fmt(WISARD_FAIL_OPEN_FILE, S_path);
			}
		}
		const char *homedir = pw->pw_dir;
		size_t N_len = strlen(S_inpPath) - k + strlen(homedir) + 1;
		S_path = new char[N_len];
		sprintf(S_path, "%s%s", pw->pw_dir, S_inpPath+k);
	} else {
#endif
		wsAlloc(S_path, char, strlen(S_inpPath) + 1);
		strcpy(S_path, S_inpPath);
#ifndef _WIN32
	}
#endif
	S_inpName	= NULL;
	S_desc		= S_inpDesc;
	B_binary	= B_inpBinary;

	/* Get the name from the path */ {
		char *Sp_end = S_path + strlen(S_path) - 1;
		while (Sp_end>=S_path && *Sp_end != WISARD_DIR_CHAR) Sp_end--;

		wsAlloc(S_inpName, char, strlen(S_path));
		strcpy(S_inpName, S_path + 1);
		//S_inpName = strdup(S_path+1);
	}

	FILE *H_fpTest = S_path[0]=='.'&&!S_path[1]?stdin : fopen(S_path, "rb");
	if (H_fpTest == NULL) {
		if (B_dontHalt)
			return STR_OPEN_FAILED;
		else
			S_desc ? halt_fmt(WISARD_FAIL_OPEN_FILEWITHDESC, S_desc, S_path) :
				halt_fmt(WISARD_FAIL_OPEN_FILE, S_path);
	}

	/* Open  */
	USHORT_t N_sig;
	if (fread(&N_sig, sizeof(USHORT_t), 1, H_fpTest)==1 &&
		(N_sig==0x8b1f || N_sig==0x8b1f)) {
		fclose(H_fpTest);
#ifndef USE_GZ
		S_desc ? halt_fmt(WISARD_FAIL_OPEN_GZWITHDESC, S_desc, S_path) :
			halt_fmt(WISARD_FAIL_OPEN_GZ, S_path);
#endif
		/* Gzip open from stdin requires its filename to '-' */
		H_stream	= gzopen(S_path[0]=='.'&&!S_path[1]?"-":S_path,
			B_binary?"rb":"r");

		B_gzip		= 1;
	} else {
		fclose(H_fpTest);
		H_stream	= S_path[0]=='.'&&!S_path[1] ? stdin :
			fopen(S_path, B_binary?"rb":"r");

		B_gzip		= 0;
	}
	if (H_stream == NULL) {
		if (B_dontHalt)
			return STR_OPEN_FAILED;
		else
			S_desc ? halt_fmt(WISARD_FAIL_OPEN_FILEWITHDESC, S_desc,
				S_path[0]=='.'&&!S_path[1]?"<<STANDARD INPUT>>":S_path) :
				halt_fmt(WISARD_FAIL_OPEN_FILE, S_path);
	}
	//	return STR_OPEN_FAILED;

	/* Successfully initialized */
	B_fail = 0;
	return 0;
}

size_t cStrFile::tell()
{
	return B_gzip ? gztell((gzFile)H_stream) : ftell((FILE *)H_stream);
}

int cStrFile::isEnd()
{
	return B_gzip ? gzeof((gzFile)H_stream) : feof((FILE *)H_stream);
}

char* cStrFile::gets(char *S_buf, wsUint N_len)
{
	char *ret = NULL;

	char B_isEOF = char(B_gzip ? gzeof((gzFile)H_stream) : feof((FILE *)H_stream));
	if (B_isEOF) return NULL;
	if (OPT_ENABLED(passemptyline)) {
		size_t N_lastPos = ftell((FILE *)H_stream);
		if (B_gzip)
			do {
				S_buf[0] = '\0';
				ret = gzgets((gzFile)H_stream, S_buf, N_len);
			} while (S_buf[0] == '\r' || S_buf[0] == '\n' || S_buf[0] == '\t' || S_buf[0] == ' ');
		else
			do {
				S_buf[0] = '\0';
				ret = fgets(S_buf, N_len, (FILE *)H_stream);
				// FIx malformed
				if (N_lastPos == ftell((FILE *)H_stream)) {
					fseek((FILE *)H_stream, 0, SEEK_END);
					return NULL;
				}
			} while (S_buf[0] == '\r' || S_buf[0] == '\n' || S_buf[0] == '\t' || S_buf[0] == ' ');
	} else {
		if (B_gzip)
			ret = gzgets((gzFile)H_stream, S_buf, N_len);
		else
			ret = fgets(S_buf, N_len, (FILE *)H_stream);
	}
	if (ret && (!ret[0] || ret[0] == '\r' || ret[0] == '\n')) {
		if (B_isEOF)
			return NULL;
		else
			pverbose("Unclear end found\n");
	}

	return ret;
}

char* cStrFile::ugets(char B_draw/*=0*/)
{
	int size = 1024;
	const int maxlen = 1024;
	char* buffer = (char*)malloc(maxlen);

	if (buffer) { /* NULL if malloc() fails */
		int ch = EOF;
		int pos = 0;

		/* Read input one character at a time, resizing the buffer as necessary */
		if (B_gzip) {
			gzFile H_zFp = (gzFile)H_stream;
			while((ch = gzgetc(H_zFp)) != '\n' && ch != EOF && !gzeof(H_zFp)) {
				buffer[pos++] = (char)ch;
				if (pos == size) { /* Next character to be inserted needs more memory */
					size = pos + maxlen;
					buffer = (char*)realloc(buffer, size);
				}
			}
		} else {
			FILE *H_fp = (FILE *)H_stream;
			while((ch = fgetc(H_fp)) != '\n' && ch != EOF && !feof(H_fp)) {
				buffer[pos++] = (char)ch;
				if (pos == size) { /* Next character to be inserted needs more memory */
					size = pos + maxlen;
					buffer = (char*)realloc(buffer, size);
				}
			}
		}
		buffer[pos] = '\0'; /* Null-terminate the completed string */
	}
	if (buffer && (!buffer[0] || buffer[0] == '\r' || buffer[0] == '\n')) {
		DEALLOC(buffer);
		return NULL;
	}
	return buffer;
}

char cStrFile::close()
{
	if (S_path) {
		free(S_path);
		S_path = NULL;
	}
	if (!H_stream) return 0;
	int ret = B_gzip ? gzclose((gzFile)H_stream) : fclose((FILE *)H_stream);
	H_stream = NULL;
	return (char)ret;
}

char cStrFile::end()
{
	int ret = B_gzip ? gzeof((gzFile)H_stream) : feof((FILE *)H_stream);
	return (char)ret;
}

char cStrFile::rewind()
{
	/* STDIN can't use rewind */
	if (S_path[0]=='.'&&!S_path[1])
		halt("Input from STDIN does not supports rewind()");
// ::rewind does not have a return value, hence this part should be separated
#ifdef USE_GZ
	if (B_gzip) {
		int ret = gzrewind((gzFile)H_stream);
		if (ret == -1) halt("GZIP library error : REWIND failed\n");
		return (char)ret;
	}
#endif
	/* Normal file stream does not have return, so always return 0 */
	::rewind((FILE *)H_stream);
	return 0;
}

void cStrFile::seek(size_t N_pos, int X_type)
{
	if (B_gzip) {
		z_off_t N_res = gzseek((gzFile)H_stream, (z_off_t)N_pos, X_type);
		if (N_res == -1)
			halt("Seek to [%d] failed", N_pos);
	} else {
		if (fseek((FILE *)H_stream, (long)N_pos, X_type) != 0)
			halt("Seek to [%d] failed", N_pos);
	}
}

size_t cStrFile::read(float *Bp, wsUint N_sz/*=1*/)
{
#ifdef USE_GZ
	return B_gzip ? gzread((gzFile)H_stream, Bp, sizeof(float)*N_sz)/sizeof(float) :
		fread(Bp, sizeof(float), N_sz, (FILE *)H_stream);
#else
	return fread(Bp, sizeof(float), N_sz, (FILE *)H_stream);
#endif
}

size_t cStrFile::read(void *Bp, wsUint N_sz/*=1*/)
{
#ifdef USE_GZ
	return B_gzip ? gzread((gzFile)H_stream, Bp, sizeof(unsigned char)*N_sz) :
		fread(Bp, sizeof(unsigned char), N_sz, (FILE *)H_stream);
#else
	return fread(Bp, sizeof(unsigned char), N_sz, (FILE *)H_stream);
#endif
}

size_t cStrFile::read(unsigned char *Bp, wsUint N_sz/*=1*/)
{
#ifdef USE_GZ
	return B_gzip ? gzread((gzFile)H_stream, Bp, sizeof(unsigned char)*N_sz) :
		fread(Bp, sizeof(unsigned char), N_sz, (FILE *)H_stream);
#else
	return fread(Bp, sizeof(unsigned char), N_sz, (FILE *)H_stream);
#endif
}

size_t cStrFile::read(char *Bp, wsUint N_sz/*=1*/)
{
#ifdef USE_GZ
	return B_gzip ? gzread((gzFile)H_stream, Bp, sizeof(char)*N_sz) :
		fread(Bp, sizeof(char), N_sz, (FILE *)H_stream);
#else
	return fread(Bp, sizeof(char), N_sz, (FILE *)H_stream);
#endif
}

size_t cStrFile::read(unsigned int *Bp, wsUint N_sz/*=1*/)
{
#ifdef USE_GZ
	return B_gzip ? gzread((gzFile)H_stream, Bp, sizeof(unsigned int)*N_sz) :
		fread(Bp, sizeof(unsigned int), N_sz, (FILE *)H_stream);
#else
	return fread(Bp, sizeof(unsigned int), N_sz, (FILE *)H_stream);
#endif
}

size_t cStrFile::read(double *Bp, wsUint N_sz/*=1*/)
{
#ifdef USE_GZ
	return B_gzip ? gzread((gzFile)H_stream, Bp, sizeof(double)*N_sz) :
		fread(Bp, sizeof(double), N_sz, (FILE *)H_stream);
#else
	return fread(Bp, sizeof(double), N_sz, (FILE *)H_stream);
#endif
}

#ifdef USE_NET
cNetFile::cNetFile(wsStrCst S_url, wsStrCst S_inpDesc/*=NULL*/)
{
	S_desc = S_inpDesc;

	if (open(S_url))
		S_desc ? halt_fmt(WISARD_FAIL_OPEN_URLWITHDESC, S_desc, S_url) :
			halt_fmt(WISARD_FAIL_OPEN_URL, S_url);
}

cNetFile::~cNetFile()
{
#ifdef _WIN32
	WSACleanup();
#endif
}

void cNetFile::_init()
{
#ifdef _WIN32
	/* Initialize socket interface */
	WSADATA	X_wsaData;
	WORD	N_wsVer	= MAKEWORD(2, 2);
	int		N_stat	= WSAStartup(N_wsVer, &X_wsaData);
	/* Initialization check */
	if (N_stat != 0)
		halt_fmt(WISARD_FAIL_INIT_WINSOCK, HIBYTE(N_wsVer), LOBYTE(N_wsVer));
	/* Version check */
	if (LOBYTE(X_wsaData.wVersion) != LOBYTE(N_wsVer) ||
		HIBYTE(X_wsaData.wVersion) != HIBYTE(N_wsVer))
		halt_fmt(WISARD_INVL_WSOCKVER, HIBYTE(X_wsaData.wVersion),
			LOBYTE(X_wsaData.wVersion), HIBYTE(N_wsVer), LOBYTE(N_wsVer));
#endif
}

char* cNetFile::_getHost(wsStrCst S_url, wsUint *Np_port/*=NULL*/)
{
	char	*S_ret	= NULL;
	wsStrCst	S_http	= "http://";

	/* URL must starts with http:// */
	if (memcmp(S_url, S_http, strlen(S_http)))
		halt_fmt(WISARD_INVL_URLFORMAT, S_url);

	/* Find the first '/' after http:// */
	wsStrCst Sp_endHost = strchr(S_url+strlen(S_http), '/');

	/* Make return buffer and copy */
	wsUint N_lenHost;
	if (Sp_endHost)
		N_lenHost = (wsUint)(Sp_endHost-S_url-strlen(S_http));
	else
		N_lenHost = (wsUint)(strlen(S_url)-strlen(S_http));
	if (N_lenHost == 0)
		halt("No host URL was found on given url [%s]", S_url);
	wsCalloc(S_ret, char, N_lenHost+1);
	memcpy(S_ret, S_url+strlen(S_http), N_lenHost);

	/* Check out there is a port suggestion */
	char *Sp_port = strchr(S_ret, ':');
	if (Sp_port != NULL) {
		*Sp_port = '\0';
		if (Np_port) *Np_port = atoi(Sp_port+1);
	} else if (Np_port) *Np_port = 80;

	/* This URL is valid */
	return S_ret;
}

char cNetFile::open(wsStrCst S_url, wsStrCst S_inpDesc, char B_dontHalt,
	char B_inpBinary)
{
	/* Find out the IP address of requested URL */
	wsUint	N_port	= 80;
	char	*S_host	= _getHost(S_url, &N_port);

	/* Get host address */
#if 1
	unsigned long N_addrHost = inet_addr(S_host);
#else
	struct in_addr iaddr;
	unsigned long N_addrHost;
	if (inet_pton(AF_INET, S_host, &iaddr) == SUCCESS)
		N_addrHost = iaddr;
#endif

	/* Make socket */
	SOCKET X_sock = socket(PF_INET, SOCK_STREAM, 0);
	if (X_sock < 0)
		halt_fmt(WISARD_FAIL_SOCKCREATE, S_url);

	/* Connect to host */
	SOCKADDR_IN	X_soin;
	memset(&X_soin, 0x00, sizeof(SOCKADDR_IN));
	X_soin.sin_family		= AF_INET;
	X_soin.sin_addr.s_addr	= htonl(N_addrHost);
#ifdef _WIN32
	X_soin.sin_port			= htons((u_short)N_port);
#else
	X_soin.sin_port			= htons(N_port);
#endif
	if (connect(X_sock, (SOCKADDR *)&X_soin, sizeof(SOCKADDR_IN)) < 0)
		halt_fmt(WISARD_FAIL_SOCKCONN, S_url);
		/* Should be implemented */

	H_stream = (void *)X_sock;

	/* Make actual HTTP request */
//	send(X_sock, 

	/* Successfully connected */
	return 0;
}

char cNetFile::close()
{
#ifdef _WIN32
	if (::closesocket((SOCKET)H_stream) != 0)
		return 1;
#else
	if (::close((intptr_t)H_stream) != 0)
		return 1;
#endif

	/* Successfully closed */
	return 0;
}

#endif

} // End namespace ONETOOL
