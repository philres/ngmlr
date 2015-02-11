/*
 * UpdateCheck.cpp
 *
 *  Created on: October 2, 2014
 *      Author: moritz
 */

#include "UpdateCheck.h"
#include "Version.h"

#include "Log.h"

#include <unistd.h>
#include <sys/socket.h>
#include <arpa/inet.h>
#include <stdlib.h>
#include <netdb.h>
#include <errno.h>
#include <cstring>
#include <time.h>

#include <string>
#include <vector>
using std::vector;
using std::string;

////////////////////////////////////////////////
// Update Check HTTP Response Body (Line End: \n)
// (Headers are ignored)
////////////////////////////////////////////////
//[UPDATE_MAGIC]
//[UP_TO_DATE|OUT_OF_DATE]
//[CURRENT_VERSION_MAJOR]
//[CURRENT_VERSION_MINOR]
//[CURRENT_VERSION_BUILD]
//[DOWNLOAD URL]
//[MESSAGE]

//#define UPDATE_DEBUG

#define UPDATE_MAGIC        "NGMUPDATE_1337_MAGIC"
#define UP_TO_DATE          "UP_TO_DATE"
#define OUT_OF_DATE         "OUT_OF_DATE"

#define UPDATE_HOST_DEFAULT "www.cibiv.at"
#define UPDATE_PATH_DEFAULT "/software/ngm/version.php"

#define UPDATE_REMIND_AFTER_MONTHS 6

#ifdef UPDATE_DEBUG
	#define UPDATE_LOG(x) Log.Message(string(x).c_str())
#else
	#define UPDATE_LOG(x) NULL
#endif

class UpdateCheckException
{
	string m_what;

public:
	UpdateCheckException(const string& what="God knows what") { m_what = what; };
	const string& what() { return m_what; }
};

#ifndef _WIN32

class SocketClient
{
protected:
	int m_socket;

public:
	 SocketClient(const string host="wwww.google.at",const string& get="/")
	 {
	 	m_socket = -1;
	 }

	~SocketClient()
	{

	}

	void create()
	{
		UPDATE_LOG("Create socket");

		//Create socket
		m_socket = socket( AF_INET, SOCK_STREAM, IPPROTO_TCP );

		if( m_socket < 0 )
			throw UpdateCheckException("Failed to create socket");
	}

	void connect(const string& host,int port)
	{
		struct hostent* ent;

		char ip[255];
		memset(ip,0,255);

		UPDATE_LOG("Resolve IP");

		ent = gethostbyname( host.c_str() );
		if( ent == 0 )
			throw UpdateCheckException("Failed to get IP for '" + host + "'");

		struct in_addr ** addr_list = (struct in_addr **) ent->h_addr_list;

		for(int i = 0; addr_list[i] != NULL; i++) 
		{
			//Return the first one;
			strcpy(ip, inet_ntoa(*addr_list[i]));
		}

		sockaddr_in* remote = new sockaddr_in();
		remote->sin_family = AF_INET;

		int set_res = inet_pton( AF_INET, ip, (void *)(&(remote->sin_addr.s_addr)) );
		if( set_res < 0 )
			throw UpdateCheckException("Cannot set IP to socket");
		else if( set_res == 0 )
			throw UpdateCheckException("IP not valid");

		remote->sin_port = htons(port);

		UPDATE_LOG("Connecting...");

		int connect_res = ::connect( m_socket, (sockaddr*) remote, sizeof(sockaddr_in) );
		if( connect_res < 0 )
		{
			throw UpdateCheckException("Failed to connect to " + string(ip));
		}

		UPDATE_LOG("Connected!");

		delete remote;
	}

	void disconnect()
	{
		close(m_socket);
	}

	void write(const string& data)
	{
		if( send( m_socket, data.c_str(), data.size(), 0 ) < 0 )
			throw UpdateCheckException("Failed to write to socket");
	}

	string read()
	{
		char buffer[8192];
		memset(buffer,0,sizeof(buffer));

		int bytecount = recv( m_socket, buffer, sizeof(buffer), 0 );

		if( bytecount < 0 )
			throw UpdateCheckException("Failed to read from socket");

		return string( buffer, bytecount );
	}
};

#endif

//No HTTP/Socket library, reinvent wheel instead
class HTTPRequest : public SocketClient
{
	string m_host;
	string m_get;

public:
	 HTTPRequest(const string host="www.google.at",const string& get="/")
	 {
	 	m_host = host;
	 	m_get = get;
	 }

	~HTTPRequest()
	{

	}

	string buildRequestString()
	{
		return "GET " + m_get + " HTTP/1.1\r\nHost: " + m_host + "\r\n\r\n";
	}

	vector<string> processReponseString(const string& resp)
	{
		//Remove HTTP headers

		int header_line_count = 0;
		string line = "";

		int i = 0;

		for( i = 0; i < resp.size(); ++ i )
		{
			if( resp[ i ] == '\r' )
			{
				i ++;
				if( i < resp.size() && resp[ i ] == '\n' )
				{
					if( line == "")
					{
						//End of HTTP header
						i++;
						break;
					}
					line = "";
					header_line_count ++;
					continue;
				}
			}

			line.push_back( resp[ i ] );
		}

		if( header_line_count == 0 )
			throw UpdateCheckException("No HTTP header in response");

		vector<string> result;
		line = "";

		for( ; i < resp.size(); ++ i )
		{
			if( resp[ i ] == '\n' )
			{
				result.push_back( line );
				line = "";
				continue;
			}

			line.push_back( resp[ i ] );
		}

		if( line != "" )
			result.push_back( line );

		return result;
	}

	vector<string> execute()
	{
		create();
		connect(m_host, 80);
		write( buildRequestString() );
		string response = read();
		disconnect();
		return processReponseString( response );
	}

};

//TODO: Warn to enable version check option when distance of curr time - build time > treshold
//TODO: Only connect to update check server when option enabled
//TODO: Allow override of update check server by option
struct tm GetBuildTime()
{
	const char mon_name[][4] = {
		"Jan", "Feb", "Mar", "Apr", "May", "Jun",
		"Jul", "Aug", "Sep", "Oct", "Nov", "Dec"
	};

	struct tm tm_build;
	memset( (void*) &tm_build, 0, sizeof(tm_build) );

	string build_stamp( __DATE__ ); //Apr  4 2014


	vector<string> parts;
	string buffer;
	for( int i = 0; i < build_stamp.size(); ++ i )
	{
		if( build_stamp[i] == ' ' )
		{
			if( buffer != "" )
			{
				parts.push_back( buffer );
				buffer = "";
			}
		} else {
			buffer.push_back( build_stamp[i] );
		}
	}

	if( buffer.size() > 0 )
		parts.push_back( buffer );

	if( parts.size() < 3 )
		throw UpdateCheckException("Unrecognized __DATE__ stamp in build");

	int month = -1;
	for( int i = 0; i < 12; i ++ )
	{
		if( string(mon_name[ i ]) == parts[0] )
		{
			month = i;
		}
	}

	string year_text = parts[2];
	int year = 0;
	int power = 1;
	for( int i = year_text.size() - 1; i >= 0; -- i )
	{
		year += power * (year_text[i]-'0');
		power *= 10;
	}

	tm_build.tm_mon = month;
	tm_build.tm_year = year - 1900;

	//Log.Message("%s//mon:%s, year:%s -ddk mon: %d, year: %d",__DATE__,parts[0].c_str(),parts[2].c_str(),month,year);

	return tm_build;
}


void UpdateCheckInterface::reminder()
{
	try
	{
		time_t time_raw;
		time(&time_raw);

		struct tm * tm_now_ptr;
		tm_now_ptr = localtime(&time_raw);

		struct tm tm_now = *tm_now_ptr;

		struct tm tm_build;
		tm_build = GetBuildTime();

		int months_now = tm_now.tm_year * 12 + tm_now.tm_mon;
		int months_build = tm_build.tm_year * 12 + tm_build.tm_mon;

		if( months_now - months_build > UPDATE_REMIND_AFTER_MONTHS )
		{
			Log.Message("[UPDATE_CHECK] Your version of NGM is more than %d months old - a newer version may be available. (For performing an automatic check use --update-check)",UPDATE_REMIND_AFTER_MONTHS);
		}

	} catch( UpdateCheckException ex ) {
		Log.Message(("[UPDATE_CHECK] Failure: " + ex.what()).c_str());
	}
}

void UpdateCheckInterface::remoteCheck()
{
	try
	{
		string params = string("?major=") + VERSION_MAJOR + string("&minor=") + VERSION_MINOR + string("&build=") + VERSION_BUILD;
		HTTPRequest req = HTTPRequest( UPDATE_HOST_DEFAULT, UPDATE_PATH_DEFAULT + params );

		vector<string> lines = req.execute();

		if( lines.size() < 7 )
			throw UpdateCheckException("Response too short");

		if( lines[0] != UPDATE_MAGIC )
			throw UpdateCheckException("Response header magic invalid");

		//If we reach this point, the update check was successful!

		Log.Message((string("[UPDATE_CHECK] Requesting version information from ") + UPDATE_HOST_DEFAULT).c_str());

		if( lines[1] == UP_TO_DATE )
		{
			Log.Message(string("[UPDATE_CHECK] You are running the latest version of NGM.").c_str());

		} else if( lines[1] == OUT_OF_DATE ) {
			Log.Message(string("[UPDATE_CHECK] >>>>>>> Your version of NGM is out of date! <<<<<<<").c_str());
			Log.Message(string("[UPDATE_CHECK] You are running  : " + string(VERSION_MAJOR) + "-" + string(VERSION_MINOR) + "-" + string(VERSION_BUILD)).c_str());
			Log.Message(string("[UPDATE_CHECK] Available version: >>> " + string(lines[2]) + "-" + string(lines[3]) + "-" + string(lines[4]) + " <<<").c_str());
			Log.Message(string("[UPDATE_CHECK] Download URL: " + lines[5]).c_str());
			Log.Message(string(lines[6]).c_str());

		} else {
			throw UpdateCheckException("Status value invalid");
		}

	} catch( UpdateCheckException ex ) {
		Log.Message(("[UPDATE_CHECK] Failure: " + ex.what()).c_str());
	}
}
