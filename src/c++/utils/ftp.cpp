#include "global/common.h"
#include "utils/util.h"
#include "utils/ftp.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#ifdef USE_FTP
#if defined(_WIN32)
#	include <winsock.h>
#	include <direct.h>
#	include <conio.h> 
#	include <io.h>
#else
#	include <unistd.h>
#	include <sys/types.h>			/* See NOTES */
#	include <arpa/inet.h>
#	include <sys/socket.h>
#	include <ctype.h>
#	include <netdb.h>
#	include <errno.h>
#	define INVALID_SOCKET -1
#	define closesocket close
#endif

namespace ONETOOL {

cFtp::cFtp() {
	Connected = 0;
	hListenSocket	= INVALID_SOCKET;
	hControlSocket	= INVALID_SOCKET;
	hDataSocket		= INVALID_SOCKET;
	bSendPort		= 1;
	ReadCommand = 0;
	bMode=ASCII;

	InitWinsock();
};

cFtp::cFtp(wsStrCst S_name, wsStrCst S_path)
{
	Connected		= 0;
	hListenSocket	= INVALID_SOCKET;
	hControlSocket	= INVALID_SOCKET;
	hDataSocket		= INVALID_SOCKET;
	bSendPort		= 1;
	ReadCommand		= 0;
	bMode			= ASCII;

	InitWinsock();
	DoOpen(S_name);
	DoLogin("user anonymous:biznok@gmail.com");
	char S_cmdCD[256];
	sprintf(S_cmdCD, "cd %s", S_path);
	DoCD(S_cmdCD);
	//	F.DoList("ls");
	//DoBinary();
	//DoGet("rs_chAltOnly.fas.gz");
};

cFtp::~cFtp() {
	
};


/*
 * DoOpen
 * this function is called when the o,open command is issued.
 * it connects to the requested host, and logs the user in
 *
 */
void cFtp::DoOpen(wsStrCst command)
{
  char *szHost=NULL;  /* remote host */

	/* 
	* do not do anything if we are already connected.
	*/
	if (Connected ) {
		printf("Already connected.  Close connection first.\n");
		fflush(stdout);
		return;
	}
	
	/*
	* extract the remote host from command line (if given ).
	* make a copy of the host name with strdup.
	*/
	if (!strcmp(command,"open") || !strcmp(command,"o")) {
		printf("Host:");
		if (!fgets(szBuffer,1024,stdin)) halt("STDERR: fgets() failed");
		(void)strtok(szBuffer,"\n");
		szHost = (char *)strdup(szBuffer);;
	}
	else if (!strncmp(command,"open ",5))
	szHost = strdup(&command[5]);
	else if (!strncmp(command,"o ",2) )
	szHost = strdup(&command[2]);
	else
	szHost = strdup(command);


	printf("Connecting to %s\n",szHost);
	hControlSocket = ConnectToServer(szHost, (char *)"21");

#if (defined(WIN32) || defined(_WIN32) )
	Sleep(1);
#else
	sleep(1);
#endif

	if (hControlSocket > 0)  {
	 printf("Connected to %s\n",szHost);
/*#if (defined(WIN32) || defined(_WIN32) )
	
	sprintf(command,"dftp: Connected to %s ", szHost);

	SetConsoleTitle(command);  // set console window title 
#endif*/

	 Connected = 1;			/* we ar now connected */
	 GetReply();				/* get reply (welcome message) from server */
	 //DoLogin((char *)NULL); /* prompt for username and password */
	 //DoBinary();				/* default binary mode */
	}
	free(szHost); /* free the strdupped string */
}



/*
 * DoLogin
 * this function logs the user into the remote host.
 * prompts for username and password
 * or parses username and password from command (user <username>:<password>)
 */
void cFtp::DoLogin(wsStrCst command)
{
 char *User=NULL;
 char *Password=NULL;
 
  if (command && *command)

	 User=strdup(&command[5]);

	if (Connected )  {

	 /* 
		* ignore leading whitespace
		*/
	 while(User && (*User == ' ' || *User == '\t') && *User)
	User++;
	Password = strchr(User, ':');
	if (Password) {
		printf("password given\n");
		*Password = '\0';
		Password++;
	};
	 /*
		* if user name was not provided via command line, read it in.
		*/
	 if (!User || !(*User) ) {
	printf("Login:");
	if (!fgets(szUser,20,stdin)) halt("STDERR: fgets() failed");
	User = szUser;
	(void)strtok(szUser,"\n");	/* remove '\n' */
	 }

	 /*
		* send user name & password to server  & get reply message
		*/
	 sprintf(szBuffer,"USER %s\r\n",User);
	 SendControlMsg(szBuffer, strlen(szBuffer));
	 GetReply();
	 if (!Password)
			halt("Password was not provided!");

	 sprintf(szBuffer,"PASS %s\r\n",Password);
	 SendControlMsg(szBuffer,strlen(szBuffer));
	 GetReply();
	}
	else
		printf("Not Connected.\n");

	free(User); /* free the strdupped string */
}

/*
 * DoClose
 * closes connection to the ftp server
 */
void cFtp::DoClose( void )
{
	if (!Connected  ) {
	 printf("Not Connected.\n");
	
	}
	else {
		SendControlMsg((char *)"quit\r\n",6);
		GetReply();
		CloseControlConnection();
		hControlSocket = (SOCKET)-1;
#if (defined(WIN32) || defined(_WIN32) )
	
		SetConsoleTitle("dftp: Connection closed");
#endif

		Connected = 0;
	}

}

/*
 * DoList
 * perform directory listing i.e: ls
 */
void cFtp::DoList( char *command)
{
	if (!Connected ) {
		printf("Not Connected.\n");
		return;
	}

	/*
	* obtain a listening socket
	*/
	if (GetListenSocket() < 0) {
		printf("Cannot obtain a listen socket.\n");
		return;
	}
	
	/*
	* parse command
	*/
	if (!strcmp(command,"ls") )  {
		sprintf(szBuffer,"NLST\r\n");
	}
	else if (!strcmp(command,"dir") ) 
		sprintf(szBuffer,"LIST\r\n");

	else if (!strncmp(command, "ls ",3)) {
		while( *command == ' ') command++;
		sprintf(szBuffer,"LIST %s\r\n",&command[3]);
	}
	/*
	* send command to server and get response
	*/
	SendControlMsg(szBuffer,strlen(szBuffer));
	memset(szBuffer,0,1024);
	GetReply();

	/*
	* accept server's connection
	*/
	if (AcceptConnection() < 0) {
		printf("Cannot accept connection.\n");
		return;
	}
	CloseListenSocket();		/* close listening socket */

	/*
	* display directory listing.
	*/
	while( ReadDataMsg(szBuffer,1024) > 0) {
		fflush(stdout);
		puts(szBuffer);
		memset(szBuffer,0,1024);
	}
	/*
	* read response
	*/
	(void)GetReply();

}

/*
 * DoCD
 * chang to another directory on the remote system
 */
void cFtp::DoCD(wsStrCst command)
{
	const char *dir=&command[2];

	if (!Connected ) {
		printf("Not Connected.\n");
		return;
	}

	/*
	* ignore leading whitespace
	*/
	while( *dir && (*dir == ' ' || *dir == '\t') ) 
		dir++;

	/*
	* if dir is not specified, read it in
	*/
	if (! (*dir) ) {
		printf("Remote directory:");
		if (!fgets(szBuffer,1024,stdin)) halt("STDERR: fgets() failed");
		(void)strtok(szBuffer,"\n");
		dir = (char *)strdup(szBuffer);
		while( *dir && (*dir) == ' ') 
		dir++;
		if (!(*dir) ) {
	printf("Usage: cd remote-directory\n");
	return;
		}
	}
	
	/*
	* send command to server and read response
	*/
	sprintf(szBuffer, "CWD %s\r\n",dir);
	SendControlMsg(szBuffer,strlen(szBuffer));
	(void)GetReply();
}

/*
 * DoLCD
 * change directory on the local system
 */
void cFtp::DoLCD(wsStrCst command)
{
	char *dir = strdup(command+3);

	while(*dir && (*dir == ' ' || *dir == '\t') ) dir++;

	/*
	* if dir is not specified, then print the current dir
	*/
	if (! *dir ) {
		dir = getcwd((char *)NULL,256);
		if (!dir)
	perror("getcwd");
		else
	printf("Current directory is: %s\n",dir);
	}
	else {
		if (chdir(dir) < 0) 
	perror("chdir");
		else {
	dir = getcwd((char *)NULL,256);
	if (!dir)
		perror("getcwd");
	else
	printf("Current directory is: %s\n",dir);
		}
	}

	free(dir);
}

/*
 * DoLLS
 * perform local directory listing.  winqvt implements this, but
 * this is not supported by many systems.  it's not really needed
 * since we already have '!' command.  you can just do !ls or !dir.
 */
void cFtp::DoLLS( char *command )
{
#if (!defined(_WIN32) || !defined(WIN32) )
	if (system("ls")) halt("SYSERR: List command have non-zero exit code");
#else
	system("dir");
#endif
}

/*
 * function to pass commands to the system, ie: dir.
 * this function gets called when '!' is encountered
 */
void DoShellCommand(char *S_cmd)
{
	S_cmd++;  /* ignore '!' */
#if (!defined(_WIN32) || !defined(WIN32) )
	if (system(S_cmd)) halt("SYSERR: Shell command [%s] have non-zero exit code", S_cmd);
#else
	if (!S_cmd || !*S_cmd)
		system("cmd");	/* have to fix this for win95 */
	else							/* maybe i should use 'start' */
		system(S_cmd); /* programatically determine which */ 
#endif									/* system we are on, then make the  */
				/* appropriate call */
}


/*
 * DoBinary
 * set file transfer mode to binary
 */
void cFtp::DoBinary()
{
  if (!Connected ) {
		printf("Not Connected.\n");
		return;
  }
	sprintf(szBuffer, "TYPE I\r\n");
	SendControlMsg(szBuffer,strlen(szBuffer));
	GetReply();
	printf("File transfer modes set to binary.\n");
	bMode = BINARY;
}

/*
 * DoAscii
 * set file transfer mode to ascii text
 */
void cFtp::DoAscii()
{
  if (!Connected ) {
		printf("Not Connected.\n");
		return;
  }
	sprintf(szBuffer, "TYPE A\r\n");
	SendControlMsg(szBuffer,strlen(szBuffer));
	GetReply();
	printf("File transfer modes set to ascii.\n");
	bMode = ASCII;

}

/*
 * DoGet
 * retrieve a file from the remote host.  calls GetFile(..)
 */
void cFtp::DoGet(wsStrCst command)
{

  if (!Connected ) {
		printf("Not Connected.\n");
		return;
  }
  //(void)strtok(command," ");
  //GetFile(strtok((char *)NULL, " "));
  GetFile(command);
}

/*
 * DoPut
 * send a file to the remote host.  calls PutFile(..)
 */
void cFtp::DoPut(wsStrCst command )
{
  if (!Connected ) {
		printf("Not Connected.\n");
		return;
  }
  //(void)strtok(command," ");
  //PutFile(strtok((char *)NULL, " "));
  PutFile(command);
}

/*
 * DoRhelp
 * sends a help command to the server.
 */
void cFtp::DoRhelp( char *command )
{
  char *szCommand;

  if (!Connected ) {
		printf("Not Connected.\n");
		return;
  }
  (void)strtok(command," ");
  szCommand=strtok((char *)NULL, " ");

  if (szCommand && *szCommand )
	 sprintf(szBuffer,"HELP %s\r\n",szCommand);
  else 
	 sprintf(szBuffer, "HELP\r\n");
  
  SendControlMsg(szBuffer, strlen(szBuffer));
  GetReply();
}

/*
 * retrieves the current directory on the remote host
 */
void cFtp::DoPWD()
{

  if (!Connected ) {
		printf("Not Connected.\n");
		return;
  }
  sprintf(szBuffer, "PWD\r\n");
  SendControlMsg(szBuffer,strlen(szBuffer));
  GetReply();
}


#ifdef _WIN32
#	define socklen_t int
#endif
SOCKET cFtp::GetListenSocket()
{
	SOCKET		sockfd;
	int			flag=1;
	socklen_t	len;
	struct sockaddr_in  serv_addr, TempAddr;
	char *port,*ipaddr;
	char szBuffer[64]={0};



	/*
	 * Open a TCP socket (an Internet stream socket).
	 */

	if ((sockfd = socket(AF_INET, SOCK_STREAM, 0)) < 0) {
		perror("socket");
		return INVALID_SOCKET;
	}



	/*
	 * Fill in structure fields for binding
	 */

	if (bSendPort) {
		bzero((char *) &serv_addr, sizeof(serv_addr));
		serv_addr.sin_family		= AF_INET;
		serv_addr.sin_addr.s_addr	= htonl(INADDR_ANY);
		serv_addr.sin_port			= htons(0); /* let system choose */
	} else {
		/* reuse the control socket then */
		if (setsockopt(sockfd,SOL_SOCKET,SO_REUSEADDR,
			(char *)&flag,sizeof(flag)) < 0) {
			perror("setsockopt");
			closesocket(sockfd);
			return INVALID_SOCKET;
		}
	}

	/*
	 * bind the address to the socket
	 */

	if (::bind(sockfd,(struct sockaddr *)&serv_addr, sizeof(serv_addr)) < 0) {
		perror("bind");
		closesocket(sockfd);
		return INVALID_SOCKET;
	}

	

	len = sizeof(serv_addr);
	if (getsockname(sockfd, (struct sockaddr *)&serv_addr, &len) < 0) {
		perror("getsockname");
		closesocket(sockfd);
		return INVALID_SOCKET;
	}

	len = sizeof(TempAddr);
	if (getsockname(hControlSocket, (struct sockaddr *)&TempAddr, &len) < 0) {
		perror("getsockname");
		closesocket(sockfd);
		return INVALID_SOCKET;
	}



	ipaddr	= (char *)&TempAddr.sin_addr;
	port	= (char *)&serv_addr.sin_port;

#define  UC(b)  (((int)b)&0xff)
	sprintf(szBuffer,"PORT %d,%d,%d,%d,%d,%d\r\n",
			UC(ipaddr[0]), UC(ipaddr[1]), UC(ipaddr[2]), UC(ipaddr[3]),
			UC(port[0]), UC(port[1]));
	/*
	 * allow ftp server to connect
	 * allow only one server
	 */
	if (listen(sockfd, 1) < 0) {
		perror("listen");
		closesocket(sockfd);
		return INVALID_SOCKET;
	}

	SendControlMsg(szBuffer,strlen(szBuffer));
	if ((wsUint)(GetReply()) == 4) {
		LOGwarn("GetListenSock failed\n");
		closesocket(sockfd);
		return INVALID_SOCKET;
	}
	hListenSocket = sockfd;

	return sockfd;

}



SOCKET cFtp::AcceptConnection()
{
	struct sockaddr_in cli_addr;
	socklen_t		clilen	= sizeof(cli_addr);
	SOCKET	sockfd	= accept(hListenSocket, (struct sockaddr *) &cli_addr,
				&clilen);

	if (sockfd < 0) {
		 perror("accept");
		return INVALID_SOCKET;
	}

	hDataSocket = sockfd;
	closesocket(hListenSocket);

	return sockfd;
}

SOCKET cFtp::ConnectToServer(char *name, wsStrCst port/*="80"*/)
{
	SOCKET s;
	unsigned int portnum;

	struct sockaddr_in server;
	struct hostent *hp;



  while( name && *name == ' ') name++;

  if (!name || ! (*name) )

		return INVALID_SOCKET;



  portnum = atoi(port);

  bzero((char *) &server, sizeof(server));



  if (isdigit(name[0])) {

	server.sin_family		= AF_INET;

	server.sin_addr.s_addr = inet_addr(name);

	server.sin_port		 = (u_short)htons(portnum);

  }

  else{ 

	if ((hp = gethostbyname(name)) == NULL)

	{
#if (!defined( _WIN32 ) || !defined( WIN32 ))

		perror("gethostbyname");
#else 
	 printf("gethostbyname: error code: %d\r\n", WSAGetLastError());

#endif
	 return INVALID_SOCKET;

	}

	bcopy(hp->h_addr,(char *) &server.sin_addr,hp->h_length);

	server.sin_family = hp->h_addrtype;

	server.sin_port = (u_short)htons(portnum);

	}/* else */



  /* create socket */

  if ((s = socket(AF_INET, SOCK_STREAM, 0)) < 1) {
#if (!defined( _WIN32 ) || !defined( WIN32 ))
	 perror("socket");
#else
	 printf("socket: error code: %d\r\n", WSAGetLastError());
#endif

	return INVALID_SOCKET;

  }

  

  if (connect(s,(struct sockaddr *)&server, sizeof(server))< 0) {

#if (!defined( _WIN32 ) || !defined( WIN32 ))
	 perror("connect");
#else
	 printf("connect: error code: %d\r\n", WSAGetLastError());
#endif	
	 

	return INVALID_SOCKET;

  }

  

  setsockopt(s,SOL_SOCKET,SO_LINGER,0,0);

  setsockopt(s,SOL_SOCKET,SO_REUSEADDR,0,0);

  setsockopt(s,SOL_SOCKET,SO_KEEPALIVE,0,0);

  hDataSocket = s;

  

  return s;

}



int cFtp::ReadControlMsg( char *szBuffer, int len)

{

  int ret;

	

	if ((ret=recv(hControlSocket,szBuffer,len,0)) <= 0)

		return 0;

	

	return ret;

}



int cFtp::SendControlMsg(const char *szBuffer, size_t len)
{
	if (send(hControlSocket, szBuffer, (int)len, 0) <= 0)
		return 0;
	return 1;
}



int cFtp::ReadDataMsg(char *szBuffer, size_t len)
{

	int ret;



	if ((ret=recv(hDataSocket,szBuffer,(int)len,0)) <= 0)

		return 0;



	return ret;
}



int cFtp::SendDataMsg( char *szBuffer, size_t len)

{

	

	if (send(hDataSocket,szBuffer,(int)len,0) <= 0)

		return 0;

	

	return 1;

}



void cFtp::CloseDataConnection(SOCKET hDataSocket)

{

#if (!defined( _WIN32 ) || !defined( WIN32 )) 

		close(hDataSocket);

#else

		closesocket(hDataSocket);

#endif



		hDataSocket = INVALID_SOCKET;

}



void cFtp::CloseControlConnection( void )

{

#if (!defined( _WIN32 ) || !defined( WIN32 )) 

		close(hControlSocket);

#else

		closesocket(hControlSocket);

#endif



		hControlSocket = INVALID_SOCKET;

}


void cFtp::CloseListenSocket( void )
{
#if (!defined( _WIN32 ) || !defined( WIN32 )) 
		close(hListenSocket);
#else
		closesocket(hListenSocket);
#endif

		hListenSocket = INVALID_SOCKET;
}


int cFtp::CheckControlMsg( char *szPtr, int len)

{

	

	return recv(hControlSocket,szPtr,len,MSG_PEEK);

}







int cFtp::GetLine()
{
	//int done=0, iRetCode =0;
	int iLen, iBuffLen=0;
	char *szPtr = szBuffer, nCode[3]={0},ch=0;

	while( (iBuffLen < 1024) && (CheckControlMsg(&ch,1)  > 0) ){

		 iLen = ReadControlMsg(&ch,1);
		iBuffLen += iLen;

		*szPtr = ch;
		szPtr += iLen;

		if (ch == '\n' )
			break;	// we have a line: return
	}

	*(szPtr+1) = (char)0;
	strncpy(nCode, szBuffer, 3);

	return (atoi(nCode));
}



int cFtp::GetReply(char* Sp_ret/*=NULL*/)
{
	int done = 0, iRetCode = 0;
	memset(szBuffer,0,1024);

	while (!done) {
		iRetCode = GetLine();
		(void)strtok(szBuffer,"\r\n");
		pverbose("FTP response [%s]\n", szBuffer);
		if (szBuffer[3] != '-' && iRetCode > 0)
			done = 1;
		if (Sp_ret) strcpy(Sp_ret, szBuffer);
		memset(szBuffer,0,1024);
	}

	return iRetCode;
}



int cFtp::CheckInput()
{
#if (!defined(WIN32) || !defined(_WIN32) )
	//int rval;
	int i;
	struct timeval timeout;

	fd_set readfds, writefds, exceptfds;
	FD_ZERO(&readfds);
	FD_ZERO(&writefds);
	FD_ZERO(&exceptfds);
	FD_CLR (fileno(stdin),&readfds);

	timeout.tv_sec = 0 ;					/* 0-second timeout. */
	timeout.tv_usec = 0 ;					/* 0 microseconds.  */
	FD_SET(fileno(stdin),&readfds);

	i=select ( fileno(stdin)+1,
		&readfds,
		&writefds,
		&exceptfds,
		&timeout);
	/* SELECT interrupted by signal - try again.  */
	if (errno == EINTR && i ==-1)  {
		return 0;

	}
	return ( FD_ISSET(fileno(stdin),&readfds) );
#else
	return (kbhit() );
#endif
}

int  cFtp::CheckFds( char *command)
{
 int rval, i;
 fd_set readfds, writefds, exceptfds;
 struct timeval timeout;
// char *szInput=command;


/*  memset(command,0,1024); */
  memset(szBuffer,0,1024);
  FD_ZERO(&readfds);
  FD_ZERO(&writefds);
  FD_ZERO(&exceptfds);
#if (!defined(WIN32) || !defined(_WIN32) )
  FD_CLR (fileno(stdin),&readfds);
#endif
  if (hControlSocket > 0) 
		FD_CLR (hControlSocket,&readfds);
  timeout.tv_sec = 0 ;					/* 1-second timeout. */
  timeout.tv_usec = 0 ;					/* 0 microseconds.  */

#if (!defined( _WIN32 ) || !defined( WIN32 ) )
  FD_SET(fileno(stdin),&readfds);
#endif

  if (hControlSocket > 0) 
		FD_SET(hControlSocket,&readfds);

  i=select ((hControlSocket > 0) ? (int)(hControlSocket+1) : 1,
				&readfds,
			&writefds,
			&exceptfds,
			&timeout);
	 /* SELECT interrupted by signal - try again.  */
	if (errno == EINTR && i ==-1)  {
		/*memset(command,0,1024);*/
		return 0;
	}


	if ((hControlSocket > 0) && FD_ISSET(hControlSocket, &readfds) ) {
		if ((rval = ReadControlMsg(szBuffer,1024)) > 0) {
			printf("%sftp>", szBuffer);
			fflush(stdout);
			return 0;
		} else {
			printf("\r\nConnection closed by server\r\n");
#if (defined(WIN32) || defined(_WIN32) )
			SetConsoleTitle("dftp: Connection closed");
#endif
			CloseControlConnection();
			hControlSocket = (SOCKET)-1;
			return 0;
		}
	}

#if (!defined(_WIN32) || !defined(WIN32) )
	if (FD_ISSET(fileno(stdin),&readfds) )  
		return (ReadCommand = GetUnixInput(command));
	return (ReadCommand = 0);
#else
	return (ReadCommand = GetWin32Input(command));
#endif
}


int cFtp::InitWinsock()
{
#ifdef _WIN32
	WSADATA WsaData;

	if (!WSAStartup(0x0101,&WsaData) ) {
		SetConsoleTitle("dftp: Connection closed");
		return 1;
	}
	else {
	printf("Winsock cannot be started\r\n");
	return 0;
	}
#else
	return 1;
#endif
}

void cFtp::CleanUp()
{
#ifdef _WIN32
	WSACleanup();
#endif
}

/*
 * GetUnixInput
 * function called to get user input on the UNIX side.  
 */

int cFtp::GetUnixInput(char *command)
{
#if (!defined( _WIN32 ) || !defined( WIN32 )) 
	char ch;
	int i;
	
	/*
	* ignore leading whitespace
	*/
	while( (ch=getchar()) == ' ' || ch == '\t') ;
	
	if (ch != '\n') {
		command[0] = ch;
		if (!fgets(&command[1],1024,stdin)) halt("STDERR: fgets() failed");
		strtok(command,"\n");
		i = strlen(command) - 1;

		while( i>0 && isspace(command[i]))
			i--;
		if (i>= 0)
			command[i+1] = 0;

	}
	return 1;
#else
		return 0;
#endif
}

/*
 * GetWin32Input
 * called to get input on the WIN32 side.
 */
int cFtp::GetWin32Input( char *command)
{
#if (defined( _WIN32 ) || defined( WIN32 )) 
	char ch;
	static int i=0;	
	
	/*
	* i, above, is static because this function can be called
	* many times before '\r' is pressed.  we set i=0 when
	* '\r' is encountered.
	*/


	while( kbhit() ) {
		if ((ch=(char)getch()) == 0)
	 getch();			/* ignore */
		else {
	if (ch == '\r') {
		 command[i] = 0;
		 i = 0;
		 printf("\r\n");
		 i = (int)strlen(command) - 1;

		 /*
			* ignore trailing whitespace
			*/
		 while( i>0 && isspace(command[i]))
		i--;
		if (i>= 0)
			command[i+1] = 0;
		i = 0;

		 return 1;
	}
	else if (ch == (char)8 && i > 0) {	/* process backspace */
		 printf("%c %c",8,8);
		 command[i] = 0;
		 i--;
	}
	else if (( ch >= 32) && ( ch <= 126) ) {
		 printf("%c", ch);
		 command[i++] = ch;
	}
		}
	}
	return 0;
#else
	return 0; 
#endif
}


/*
 * GetFile
 * called to retrive a file from remote host
 */
void cFtp::GetFile(wsStrCst fname)
{
	FILE *fp=NULL;
//	int retval;
	int	nTotal=0, nBytesRead=0, aborted=0;
	const char *abortstr = (const char *)"ABOR\r\n";
	char ch;

	/*/void (*OldHandler)(int); */

	/*
	* did we get a filename?
	*/
	if (!fname || ! (*fname)) {
		printf("No file specified.\n");
		return;
	}

	/*
	* open the file with current mode
	*/
	if (! (fp=fopen(fname,(bMode==ASCII) ? "wt" : "wb"))) {
		perror("file open");
		return;
	}

	/*
	* send command to server and read response
	*/
	if (1) {
		char S_pasv[1024];
		strcpy(szBuffer, "PASV\r\n");
		if (!SendControlMsg(szBuffer,strlen(szBuffer))) {
			fclose(fp);
			return;
		}
		(void)GetReply(S_pasv);
		int code, v1, v2, v3, v4, p1, p2;
		int Nscan = sscanf(S_pasv, "%d Entering Passive Mode (%d,%d,%d,%d,%d,%d)",
			&code, &v1, &v2, &v3, &v4, &p1, &p2);
		/* If it is failed to interpret */
		if (Nscan != 7) {
			LOGwarn("Failed to parse PASV message [%s]\n", S_pasv);
			fclose(fp);
			return;
		}

		wsUint N_port = p1*256 + p2;
		/* Establish new port */
		char S_host[256];
		char S_port[256];
		sprintf(S_host, "%d.%d.%d.%d", v1, v2, v3, v4);
		sprintf(S_port, "%d", N_port);
		ConnectToServer(S_host, S_port);
	} else {
		/*
		* obtain a listen socket
		*/
		if (GetListenSocket() < 0) {
			fclose(fp);
			return;
		}
	}
	sprintf(szBuffer,"RETR %s\r\n",fname);
	if (!SendControlMsg(szBuffer,strlen(szBuffer))) {
		fclose(fp);
		return;
	}
	if ((wsUint)(GetReply()/100) == 4) {
		LOGwarn("RETR failed\n");
		fclose(fp);
		return;
	}
	
	/*
	* accept server connection
	*/
	if (1) {
		
	} else {
		if (AcceptConnection() <= 0) {
			fclose(fp);
			return;
		}
	}
	/* 
	* now get file and store
	*/

	
	/*/ this did not work :((
	*  so i am taking it out.  well, i am commenting it out
	*  if you figure it out, then use it.  the way i did it
	*  now works well.
	*/
/*/	OldHandler = signal(SIGINT,SigHandler);
	 retval = setjmp(abortfile);

	if (retval != 0 ) {
		printf("\r\nAborting\r\n");
		fflush(stdout);
		close(fd);
		CloseDataConnection(hDataSocket);
		if (OldHandler)
	 (void)signal(SIGINT,OldHandler);
		ControlCHit = 0;
		GetReply();
		return;
	}
*/

#if (defined(_WIN32) || defined(WIN32))
	
	printf("Press ESC to abort\r\n");
#else
	printf("Type q and hit return to abort\r\n");
#endif
	while( (nBytesRead=ReadDataMsg(szBuffer,1024)) > 0) {
		fwrite(szBuffer, 1, nBytesRead, fp);	
//		write(fd,szBuffer,nBytesRead);
		nTotal+=nBytesRead;
		printf("%s : %d received\r",fname,nTotal);
		if (CheckInput() ) {
#if (defined(_WIN32) || defined(WIN32))

			if (!(ch = (char)getch()) ) getch();
			if (ch == 27 )
				aborted = 1;
#else
			ch = getchar();
			if (ch != '\n') {
				while( getchar() != '\n') ;		/* read 'til new line */
			}
			if (ch == 'q') 
				aborted = 1;
#endif
		}
		
		/*
		* did we abort?
		*/
		if (aborted ) {
		 
			printf("\r\nAbort: Waiting for server to finish.");
			SendControlMsg(abortstr,strlen(abortstr));
			break;
		}
	}

	if (aborted ) {			// ignore everything if aborted.

		while( (nBytesRead=ReadDataMsg(szBuffer,1024)) > 0);
		GetReply();
	}

 /*  (void)signal(SIGINT,OldHandler); */

	printf("\r\n");
	fclose(fp);
//	close(fd);
	CloseDataConnection(hDataSocket);
	GetReply();
}

/*
 * PutFile
 * called to transfer a file to the remote host using the current
 * file transfer mode.  it's just like GetFile.
 *
 * i have commented out lines that would have helped with trapping
 * ctrl-c because longjmp did not work :((
 * if you figure it out, then uncomment them
 */
void cFtp::PutFile(wsStrCst fname)
{
	FILE *fp=NULL;
//	int retval;
	int fd, nTotal=0, nBytesRead=0, aborted=0;
	const char *abortstr = (const char *)"ABOR\r\n";
	char ch;

  /* void (*OldHandler)(int); */

	if (!fname || ! (*fname)) {
		printf("No file specified.\n");
		return;
	}
	if (! (fp=fopen(fname,(bMode==ASCII) ? "rt" : "rb"))) {
		perror("file open");
		return;
	}

	if (GetListenSocket() < 0) {
		fclose(fp);
		return;
	}
	
	/*
	* send command to server & read reply
	*/
	sprintf(szBuffer,"STOR %s\r\n",fname);
	if (!SendControlMsg(szBuffer,strlen(szBuffer))) {
		fclose(fp);
		return;
	}
	GetReply();
	
	/*
	* accept server connection
	*/
	if (AcceptConnection() <= 0) {
		fclose(fp);
		return;
	}
	/* 
	* now send file
	*/
	
	fd = fileno(fp);
/*
	OldHandler = signal(SIGINT,SigHandler);
	retval = setjmp(abortfile);

	if (retval != 0 ) {
		printf("Aborting\r\n");
		fflush(stdout);
		close(fd);
		CloseDataConnection(hDataSocket);
		if (OldHandler)
	 (void)signal(SIGINT,OldHandler);
		ControlCHit = 0;
		GetReply();
		return;
	}
  */
#if (defined(_WIN32) || defined(WIN32))
	
	printf("Press ESC to abort\r\n");
#else
	printf("Type q and hit return to abort\r\n");
#endif
	while( (nBytesRead=read(fd,szBuffer,1024)) > 0) {
		SendDataMsg(szBuffer,nBytesRead);
		nTotal+=nBytesRead;
		printf("%s : %d sent\r",fname,nTotal);
		if (CheckInput() ) {
#if (defined(_WIN32) || defined(WIN32))

			if (!(ch = (char)getch()) ) getch();
			if (ch == 27 )
				aborted = 1;
#else
			ch = getchar();
			if (ch != '\n') {
				while( getchar() != '\n') ;		/* read 'til new line */
			}
			if (ch == 'q') 
				aborted = 1;
#endif
		}

		/*
		* send an abort command to server if we aborted.
		*/
		if (aborted ) {
			LOGwarn("Abort: Waiting for server to finish.\n");
			SendControlMsg(abortstr,strlen(abortstr));
			break;
		}
	}

	/*(void)signal(SIGINT,OldHandler); */

	printf("\r\n");
	/*
	* close data connection
	*/
	CloseDataConnection(hDataSocket);
	close(fd);
	GetReply();
	fclose(fp);
}

} // End namespace ONETOOL

#endif
