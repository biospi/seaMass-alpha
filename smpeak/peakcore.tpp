//
// $Id$
//
//
// Original author: Ranjeet Bhamber <ranjeet <a.t> liverpool.ac.uk>
//
// Copyright (C) 2015  Biospi Laboratory for Medical Bioinformatics, University of Liverpool, UK
//
// This file is part of seaMass.
//
// seaMass is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// seaMass is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with seaMass.  If not, see <http://www.gnu.org/licenses/>.
//

#ifndef SMPEAK_PEAKCORE_TPP_
#define SMPEAK_PEAKCORE_TPP_

template<typename T>
VecMat<T>::VecMat(void):row(0),col(0),m(NULL){}

template<typename T>
VecMat<T>::VecMat(hsize_t _r, hsize_t _c, vector<T> &_vec):row(_r), col(_c),v(_vec)
{
	matIdx.resize(row,0);
	for(hsize_t i=0; i < row; ++i)
		matIdx[i]=&v[i*col];
	m=&matIdx[0];
}

template<typename T>
VecMat<T>::VecMat(hsize_t _r, hsize_t _c) : row(_r), col(_c)
{
	v.resize(row*col);
	matIdx.resize(row,0);
	for(hsize_t i=0; i < row; ++i)
		matIdx[i]=&v[i*col];
	m=&matIdx[0];
}

template<typename T>
void VecMat<T>::set(hsize_t _r, hsize_t _c, vector<T> &_vec)
{
	row=_r;
	col=_c;
	v=_vec;
	matIdx.resize(row,0);
	for(hsize_t i=0; i < row; ++i)
		matIdx[i]=&v[i*col];
	m=&matIdx[0];
}

template<typename T>
void VecMat<T>::set(hsize_t _r, hsize_t _c)
{
	row=_r;
	col=_c;
	v.resize(row*col,0);
	matIdx.resize(row,0);
	for(hsize_t i=0; i < row; ++i)
		matIdx[i]=&v[i*col];
	m=&matIdx[0];
}

template<typename T>
void VecMat<T>::getDims(hsize_t dims[])
{
	dims[0]=row;
	dims[1]=col;
}

template<typename T>
void VecMat<T>::clear(void)
{
	row=0;
	col=0;
	vector<T>().swap(this->v);
	vector<T*>().swap(this->matIdx);
}

template<typename T>
void findVecString(vector<char> &vecStr, vector<T> &vec,
		const string subStr, const string endSubStr)
{
	size_t nSub = subStr.length();
	string str(&vecStr[0]);
	if(vec.size()>0) vec.resize(0);

	for(size_t i=0; i < (vecStr.size() - nSub);)
	{
		size_t pos=str.find(subStr,i);
		if(pos == string::npos)
		{
			pos = str.find(endSubStr,i);
			vec.push_back(pos);
			i=vecStr.size();
		}
		else{
			vec.push_back(pos);
			i=pos+nSub;
		}
	}
}



template<typename T>
vector<size_t> findSize(VecMat<T> data)
{
	hsize_t dims[2];
	data.getDims(dims);
	vector<size_t> dsize;

	for(size_t idx=0; idx < dims[0]; ++idx)
	{
		size_t len=0;
		while(len < dims[1] && data.m[idx][len] > 0) ++len;
		dsize.push_back(len);
	}
	return dsize;
}

template<typename T>
vector<size_t> findSizeT(VecMat<T> data)
{
	hsize_t dims[2];
	data.getDims(dims);
	vector<size_t> dsize;

	for(size_t idx=0; idx < dims[1]; ++idx)
	{
		size_t len=0;
		while(len < dims[0] && data.m[len][idx] > 0) ++len;
		dsize.push_back(len);
	}
	return dsize;
}


template<typename T>
void repackPeakData(VecMat<T> &peak, VecMat<T> &raw, vector<int> msl,
				vector<size_t> psize, vector<size_t> rsize)
{
	vector<vector<T> > dataBuff;

	hsize_t pdim[2];
	hsize_t rdim[2];
	hsize_t r;
	hsize_t c;

	peak.getDims(pdim);
	raw.getDims(rdim);

	(pdim[0] > rdim[0]) ?  r = pdim[0]: r = rdim[0];
	(pdim[1] > rdim[1]) ?  c = pdim[1]: c = rdim[1];

	dataBuff.resize(r);

	for(hsize_t i = 0; i < r; ++i)
		dataBuff[i].resize(c,0);

	for(size_t idx=0, rdx=0; idx < msl.size(); ++idx)
	{
		if((msl[idx] == 1) && (rdx < psize.size() ))
		{
			for(int i = 0; i < psize[rdx]; ++i)
			{
				dataBuff[idx][i]=peak.m[rdx][i];
			}
			++rdx;
		}
		else if(msl[idx] == 2)
		{
			for(int i = 0; i < rsize[idx]; ++i)
			{
				dataBuff[idx][i]=raw.m[idx][i];
			}
		}
		else if((msl[idx] == 1) && (rdx >= psize.size() ))
		{
			cout << "Repack old data in missing ms-level!!!"<<endl;
			for(int i = 0; i < rsize[idx]; ++i)
			{
				dataBuff[idx][i]=raw.m[idx][i];
			}
		}
		else
		{
			cout << "Error in Repack!!! Unknown ms-level!!!"<<endl;
		}
	}

	raw.clear();
	raw.set(r,c);

	for(int i=0; i < r; ++i)
		for(int j=0; j < c; ++j)
			raw.m[i][j]=dataBuff[i][j];
}

template<typename T>
void repackPeakDataT(VecMat<T> &peak, VecMat<T> &raw, vector<int> msl,
				vector<size_t> psize, vector<size_t> rsize)
{
	vector<vector<T> > dataBuff;

	hsize_t pdim[2];
	hsize_t rdim[2];
	hsize_t r;
	hsize_t c;

	peak.getDims(pdim);
	raw.getDims(rdim);

	(pdim[0] > rdim[0]) ?  r = pdim[0]: r = rdim[0];
	(pdim[1] > rdim[1]) ?  c = pdim[1]: c = rdim[1];

	dataBuff.resize(r);

	for(hsize_t i = 0; i < r; ++i)
		dataBuff[i].resize(c,0);

	for(size_t idx=0, rdx=0; idx < msl.size(); ++idx)
	{
		if((msl[idx] == 1) && (rdx < psize.size() ))
		{
			for(int i = 0; i < psize[rdx]; ++i)
			{
				dataBuff[i][idx]=peak.m[i][rdx];
			}
			++rdx;
		}
		else if(msl[idx] == 2)
		{
			for(int i = 0; i < rsize[idx]; ++i)
			{
				dataBuff[i][idx]=raw.m[i][idx];
			}
		}
		else if((msl[idx] == 1) && (rdx >= psize.size() ))
		{
			cout << "Repack old data in missing ms-level!!!"<<endl;
			for(int i = 0; i < rsize[idx]; ++i)
			{
				dataBuff[i][idx]=raw.m[i][idx];
			}
		}
		else
		{
			cout << "Error in Repack!!! Unknown ms-level!!!"<<endl;
		}
	}

	raw.clear();
	raw.set(r,c);

	for(int i=0; i < r; ++i)
		for(int j=0; j < c; ++j)
			raw.m[i][j]=dataBuff[i][j];
}

#endif /* SMPEAK_PEAKCORE_TPP_ */
