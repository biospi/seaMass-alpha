#include "NetCDFile.hpp"

NetCDFile::NetCDFile(const string _fileName, int omode) : fileName(_fileName)
{
	if (omode == NC_NOWRITE)
	{
		if ((retval = nc_open(fileName.c_str(), omode, &ncid)))
			ERR(retval);
		fileStatus = true;
	}
}

void NetCDFile::open(const string _fileName, int omode)
{
	if (fileStatus == false)
	{
		fileName=_fileName;
		if (omode == NC_NOWRITE)
		{
			if ((retval = nc_open(fileName.c_str(), omode, &ncid)))
				ERR(retval);
			fileStatus = true;
		}
	}
	else
	{
		cout<<"File already opened"<<endl;
		exit(ERRCODE);
	}
}

void NetCDFile::close(void)
{
	if (fileStatus == true)
	{
		if ((retval = nc_close(ncid)))
			ERR(retval);
		fileStatus = false;
	}
	else
	{
		cout<<"No file to close"<<endl;
		exit(ERRCODE);
	}
}

NetCDFile::~NetCDFile()
{
	if (fileStatus == true)
	{
		if ((retval = nc_close(ncid)))
			ERR(retval);
	}
}

void NetCDFile::err(int e)
{
	cout<<"Error: "<<nc_strerror(e)<<endl;
	exit(ERRCODE);
}

