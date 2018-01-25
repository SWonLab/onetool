#pragma once
#ifndef __WISARD_FTP_H__
#define __WISARD_FTP_H__
//#include <io.h>

#define DEFAULT_PORT_NUM 21
#define PASSWORD_LENGTH 256
#ifdef _WIN32
#	undef EINTR
#	define EINTR WSAEINTR 
#else
#	define SOCKET int
#endif
#ifndef bzero
#	define bzero(x,y) memset(x,0,y)
#	define bcopy(x,y,z) memcpy(y,x,z)
#endif


#ifdef __FreeBSD__
#	include <netinet/in.h>
#endif

namespace ONETOOL {

enum {
	LS = 0,
	BINARY,
	ASCII,
	PWD,
	CD,
	OPEN,
	CLOSE,
	QUIT,
	LCD,
	LLS,
	LDIR,
	USER,
	SHELL,
	IGNORE_COMMAND,
	GET,
	PUT,
	HELP,
	RHELP,

	FTP_COMPLETE=1, /* these will be used later */
	FTP_CONTINUE,
	FTP_ERROR
};


class cFtp
{
public:
	cFtp();
	cFtp(wsStrCst S_name, wsStrCst S_path);
	~cFtp();

	void DoOpen(wsStrCst);
	void DoList(char *);
	void DoCD(wsStrCst);
	void DoShellCommand(char *);
	void DoLogin(wsStrCst);
	void DoClose(void);
	void DoLCD(wsStrCst);
	void DoPut(wsStrCst);
	void DoGet(wsStrCst);
	void DoLLS(char * );
	void DoBinary();
	void DoRhelp( char *);
	void DoAscii();
	void DoPWD();

	int  CheckFds(char *);

private:
	char szBuffer[1025];  /* buffer used to read/write */
	char szUser[20];          /* stores username */
	char szPass[256];         /* stores user password */
	int Connected;     /* flag for connect status */


	SOCKET hListenSocket;
	SOCKET hControlSocket;
	SOCKET hDataSocket;	
	int bSendPort;
	int ReadCommand;
	int bMode;

	int		GetReply(char* Sp_ret=NULL);
	int		GetLine();
	void	CleanUp();
	int		SendControlMsg(const char *, size_t);
	int		SendDataMsg( char *szBuffer, size_t len);
	SOCKET	ConnectToServer(char *name, wsStrCst port="80");
	SOCKET	GetListenSocket();
	int		InitWinsock();
	SOCKET	AcceptConnection();
	void	CloseControlConnection( void );
	void	CloseDataConnection(SOCKET hDataSocket);
	void	CloseListenSocket();
	int		ReadDataMsg( char *szBuffer, size_t len);
	int		GetUnixInput(char *command);
	int GetWin32Input( char *command);
	void GetFile(wsStrCst fname);
	void PutFile(wsStrCst fname);
	int ReadControlMsg( char *szBuffer, int len);
	int CheckControlMsg( char *szPtr, int len);
	int CheckInput();

};

} // End namespace ONETOOL

#endif
