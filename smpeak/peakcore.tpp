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
	v.resize(row*col);
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


#endif /* SMPEAK_PEAKCORE_TPP_ */
