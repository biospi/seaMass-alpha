#include "NetCDFile.hpp"

NetCDFile::NetCDFile(const string _fileName, int omode) : fileName(_fileName)
{
	switch(omode)
	{
	case NC_NOWRITE:
		if ((retval = nc_open(fileName.c_str(), omode, &ncid)))
			err(retval);
		fileStatus = true;
		break;
	case NC_WRITE:
		if ((retval = nc_open(fileName.c_str(), omode, &ncid)))
			err(retval);
		fileStatus = true;
		break;
	case NC_NETCDF4:
	   if ((retval = nc_create(fileName.c_str(), omode|NC_CLOBBER, &ncid)))
	      err(retval);
	   fileStatus = true;
	   break;
	}
}

void NetCDFile::open(const string _fileName, int omode)
{
	if (fileStatus == false)
	{
		fileName=_fileName;
		switch(omode)
		{
		case NC_NOWRITE:
			if ((retval = nc_open(fileName.c_str(), omode, &ncid)))
				err(retval);
			fileStatus = true;
			break;
		case NC_WRITE:
			if ((retval = nc_open(fileName.c_str(), omode, &ncid)))
				err(retval);
			fileStatus = true;
			break;
		case NC_NETCDF4:
			if ((retval = nc_create(fileName.c_str(), omode|NC_CLOBBER, &ncid)))
				err(retval);
			fileStatus = true;
			break;
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
			err(retval);
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
			err(retval);
	}
}

void NetCDFile::err(int e)
{
	cout<<"Error: "<<nc_strerror(e)<<endl;
	exit(ERRCODE);
}

