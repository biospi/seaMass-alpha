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
void VecMat<T>::getDims(hsize_t dims[])
{
	dims[0]=row;
	dims[1]=col;
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

#endif /* SMPEAK_PEAKCORE_TPP_ */
