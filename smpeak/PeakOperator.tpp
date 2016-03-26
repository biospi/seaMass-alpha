//
// $Id$
//
//
// Author: Ranjeet Bhamber <ranjeet <a.t> liverpool.ac.uk>
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

#ifndef SMPEAK_PEAKOPERATOR_TPP_
#define SMPEAK_PEAKOPERATOR_TPP_

template<typename T, typename R>
void MathOp<T,R>::calMidPoint(lli rtIdx, lli mzIdx, T** alpha,
		const vector<double> &mza, double &mz1, double &a1)
{
	mz1=0.5*(mza[mzIdx]+mza[mzIdx+1]);
	a1=double(0.5*(alpha[rtIdx][mzIdx]+alpha[rtIdx][mzIdx+1]));
}

template<typename T, typename R>
double MathOp<T,R>::calT(const double y0, const double y1, const double y2)
{
	double a=y1-y0;
	double b=y0-2.0*y1+y2;
	double g=a/b;
	double r=sqrt(g*g-y0/b);
	double tp=-g+r;
	double tn=-g-r;

	if(tp >= 0 && tp <= 1) return tp;
	else if(tn >= 0 && tn <= 1)	return tn;
	else
	{
		//cout<<"Insignificant False Peak Detected - Peak Ignored"<<endl;
		return -1.0;
	}
}

template<typename T, typename R>
double MathOp<T,R>::calX(double t, double x0, double x1, double x2)
{
	double c=(1.0-t);
	return c*c*x0+2.0*c*t*x1+t*t*x2;
}

template<typename T, typename R>
vector<T> MathOp<T,R>::cal3rdMidPoint(lli row, lli col, T **P)
{
	vector<T> ry(4,0.0);
	T p0=P[row][col-1];
	T p1=P[row][col];
	T p2=P[row][col+1];
	T p3=P[row][col+2];
	T u=1.0/3.0;
	T v=2.0/3.0;
	T w=1.0/2.0;
	T up=v; // (1-u);
	T vp=u; // (1-v);
	T wp=w; // (1-w);

	ry[0]=wp*(vp*p0+v*p1)+w*(up*p1+u*p2);
	ry[1]=up*p1+u*p2;
	ry[2]=vp*p1+v*p2;
	ry[3]=wp*(vp*p1+v*p2)+w*(up*p2+u*p3);

	return ry;
}

template<typename T, typename R>
T MathOp<T,R>::calPeakCount(vector<T> &ry, double t)
{
	double tp=1-t;
	T p0=tp*ry[0]+t*ry[1];
	T p1=tp*ry[1]+t*ry[2];
	T p2=tp*ry[2]+t*ry[3];
	T q0=tp*p0+t*p1;
	T q1=tp*p1+t*p2;
	return tp*q0+t*q1;
}

template<typename T, typename R>
void MathOp<T,R>::calPeakWidth(lli rtIdx,lli mzIdx, T** alpha, const vector<double> d2mz,
		double &mzlhs, double &mzrhs)
{
	lli n=d2mz.size();
	lli lhsIdx=mzIdx+1;
	lli rhsIdx=mzIdx+1;

	do
	{
		--lhsIdx;
		if(lhsIdx < 0) break;
	} while( alpha[rtIdx][lhsIdx] < 0);
	do
	{
		++rhsIdx;
		if(rhsIdx > n) break;
	} while( alpha[rtIdx][rhsIdx] < 0);

	if(lhsIdx > 0 && rhsIdx < n ){
		double x1,x2,y1,y2,m;

		x1=d2mz[lhsIdx];
		y1=alpha[rtIdx][lhsIdx];
		x2=d2mz[lhsIdx+1];
		y2=alpha[rtIdx][lhsIdx+1];
		m=(y2-y1)/(x2-x1);
		mzlhs=(-y1/m)+x1;

		x1=d2mz[rhsIdx];
		y1=alpha[rtIdx][rhsIdx];
		x2=d2mz[rhsIdx-1];
		y2=alpha[rtIdx][rhsIdx-1];
		m=(y2-y1)/(x2-x1);
		mzrhs=(-y1/m)+x1;

		/*
		if (mzlhs < d2mz[0])
		{
			cout<<"WTF LHS is out of range!!"<<endl;
		}
		if (mzrhs > d2mz[n-1])
		{
			cout<<"WTF RHS is out of range!!"<<endl;
		}
		*/
	}
	else
	{
		mzlhs=-1;
		mzrhs=-1;
	}

}

template<typename T, typename R>
void MathOp<T,R>::calPeakMZ(
		DataAxis<T,R> const *bs,DataAxis<T,R> const *dbs,DataAxis<T,R> const *d2bs,
		lli i, lli j, double &mzPeak, T &countMax, double &mzlhs, double &mzrhs,
		double &t0, lli &falsePeak)
{
	if((dbs->alpha->m[i][j] > 0) && (dbs->alpha->m[i][j+1] < 0) )
	{
		double pa1=0.0;
		double pmz1=0.0;
		calMidPoint(i,j,dbs->alpha->m,dbs->mz,pmz1,pa1);
		if(pa1 < 0){
			double pa0=0.0;
			double pmz0=0.0;
			vector<T> ry;
			calMidPoint(i,j-1,dbs->alpha->m,dbs->mz,pmz0,pa0);
			t0= calT(pa0,double(dbs->alpha->m[i][j]),pa1);
			if(t0>=0)
			{
				mzPeak=calX(t0,pmz0,dbs->mz[j],pmz1);
				mzlhs=0.0;
				mzrhs=0.0;
				ry = cal3rdMidPoint(i,j-1,bs->alpha->m);
				countMax = calPeakCount(ry,t0);
				calPeakWidth(i,j,d2bs->alpha->m,d2bs->mz,mzlhs,mzrhs);
			}
			else
			{
				++falsePeak;
			}
		}
		else
		{
			double pa2=0.0;
			double pmz2=0.0;
			vector<T> ry;
			calMidPoint(i,j+1,dbs->alpha->m,dbs->mz,pmz2,pa2);
			t0 = calT(pa1,double(dbs->alpha->m[i][j+1]),pa2);
			if (t0>=0)
			{
				mzPeak= calX(t0,pmz1,dbs->mz[j+1],pmz2);
				mzlhs=0.0;
				mzrhs=0.0;
				ry = cal3rdMidPoint(i,j,bs->alpha->m);
				countMax = calPeakCount(ry,t0);
				calPeakWidth(i,j,d2bs->alpha->m,d2bs->mz,mzlhs,mzrhs);
			}
			else
			{
				++falsePeak;
			}
		}
	}
}


template<template<class> class pPeak, template<class,class> class pData, typename T, typename R>
void Centroid1D<pPeak,pData,T,R>::calculate(pPeak<T> *peak, pData<R,T> *data, T threshold)
{
	vector<DataAxis<T,R>* > bsData =  data->get();

	DataAxis<T,R> const *bs=bsData[0];
	DataAxis<T,R> const *dbs=bsData[1];
	DataAxis<T,R> const *d2bs=bsData[2];

	hsize_t dims[2];
	bs->alpha->getDims(dims);
	lli i = 0;
	lli col = dims[1];

	lli falsePeak=0;
	lli falseWidth=0;

	// Find Peaks and exact MZ values.
	for(lli j = 2; j < col-2; ++j)
	{
		if((dbs->alpha->m[i][j] > 0) && (dbs->alpha->m[i][j+1] < 0) &&
			bs->alpha->m[i][j] > threshold)
		{
			double mzPeak=-1.0;
			T countMax=-1.0;
			double mzlhs=-1.0;
			double mzrhs=-1.0;
			double t0=-1.0;
			this->calPeakMZ(bs,dbs,d2bs,i,j,mzPeak,countMax,mzlhs,mzrhs,t0,falsePeak);

			if (t0 >=0)
			{
				if(mzlhs >= 0 && mzrhs >= 0 && countMax >= threshold)
				{
					{
						peak->addPeak(mzPeak,bs->rt[i].second,countMax,
								make_pair(mzlhs,mzrhs),
								make_pair(bs->rt[i].second,bs->rt[i].second),t0,j,
								bs->rt[i].first);
					}
				}
				else
				{
					++falseWidth;
				}
			}
		}
	}
	peak->updateFalseData(falsePeak,falseWidth);
}


template<template<class> class pPeak, template<class,class> class pData, typename T, typename R>
void Centroid2D<pPeak,pData,T,R>::calculate(pPeak<T> *peak, pData<R,T> *data, T threshold)
{
	vector<DataAxis<T,R>* > bsData =  data->get();

	DataAxis<T,R> const *bs=bsData[0];
	DataAxis<T,R> const *dbs=bsData[1];
	DataAxis<T,R> const *d2bs=bsData[2];

	hsize_t dims[2];
	bs->alpha->getDims(dims);
	lli row = dims[0];
	lli col = dims[1];

	lli falsePeak=0;
	lli falseWidth=0;
	// Find Peaks and exact MZ values.
	cout<<"Extract Peaks from Mass Spec Data"<<endl;

	#pragma omp parallel for reduction(+:falsePeak,falseWidth)
	for(lli i = 0; i < row; ++i)
	{
		for(lli j = 2; j < col-2; ++j)
		{
			if((dbs->alpha->m[i][j] > 0) && (dbs->alpha->m[i][j+1] < 0) &&
				bs->alpha->m[i][j] > threshold)
			{
				double mzPeak=-1.0;
				T countMax=-1.0;
				double mzlhs=-1.0;
				double mzrhs=-1.0;
				double t0=-1.0;
				this->calPeakMZ(bs,dbs,d2bs,i,j,mzPeak,countMax,mzlhs,mzrhs,t0,falsePeak);

				if (t0 >=0)
				{
					if(mzlhs >= 0 && mzrhs >= 0 && countMax >= threshold)
					{
						#pragma omp critical(peak)
						{
							peak->addPeak(mzPeak,bs->rt[i],countMax,make_pair(mzlhs,mzrhs),
								make_pair(bs->rt[i],bs->rt[i]),t0,j,i);
						}
					}
					else
					{
						++falseWidth;
					}
				}
			}
		}
		//cout<<"Found 2D ["<<peak->numOfPeaks()<<"] Peaks."<<endl;
	}
	cout<<"Found ["<<peak->numOfPeaks()<<"] Peaks."<<endl;
	if(falsePeak > 0)
		cout<<"Warning !!! Found ["<<falsePeak<<"] Insignificant False Peaks Detected - Peaks Ignored"<<endl;
	if(falseWidth > 0)
		cout<<"Warning !!! Found ["<<falseWidth<<"] False Peaks Detected with Incorrect Peak Widths - Peaks Ignored"<<endl;
}


template<template<class> class pPeak, template<class,class> class pData, typename T, typename R>
void ExtractPeak<pPeak,pData,T,R>::calculate(pPeak<T> *peak, pData<R,T> *data, T threshold)
{
	// Patch used for calculations 5x4, 5x7, 7x7 Coefficients,
	// 3 extra points -2/3,-1/3,1/3 - 4 extra points -2/3,-1/3,1/3,2/3
	//     0  1  2  3
	// 0   *  *  *  *
	//                                0  1  2  3
	// 1   *  *  *  *              0  *  *  *  *
	//     -  -  -  -              1  -  -  -  -
	//     -  -  -  -              2  -  -  -  -
	// 2   *  X  *  * mul 5x5 -->  3  *  X  *  *
	//     -  -  -  -              4  -  -  -  -
	// 3   *  *  *  *
	//
	// 4   *  *  *  *
	//
	//    0  1  2  3  4  5  6
	// 0  *  *  *  *  *  *  *
	//                                          0  1  2  3  4  5  6
	// 1  *  *  *  *  *  *  *                0  *  *  *  *  *  *  *
	//    -  -  -  -  -  -  -                1  -  -  -  -  -  -  -
	//    -  -  -  -  -  -  -                2  -  -  -  -  -  -  -
	// 2  *  *  *  X  *  *  *   mul 5x5 -- > 3  *  *  *  X  *  *  *
	//    -  -  -  -  -  -  -                4  -  -  -  -  -  -  -
	// 3  *  *  *  *  *  *  *
	//
	// 4  *  *  *  *  *  *  *

	//    0  1  2  3  4  5  6
	// 0  *  *  *  *  *  *  *
	//                                          0  1  2  3  4  5  6
	// 1  *  *  *  *  *  *  *                0  *  *  *  *  *  *  *
	//    -  -  -  -  -  -  -                1  -  -  -  -  -  -  -
	//    -  -  -  -  -  -  -                2  -  -  -  -  -  -  -
	// 2  *  *  *  X  *  *  *   mul 5x5 -- > 3  *  *  *  X  *  *  *
	//    -  -  -  -  -  -  -                4  -  -  -  -  -  -  -
	//    -  -  -  -  -  -  -                5  -  -  -  -  -  -  -
	// 3  *  *  *  *  *  *  *                6  *  *  *  *  *  *  *
	//
	// 4  *  *  *  *  *  *  *
	//

	vector<DataAxis<T,R>* > bsData;
	BasisPatch<T> *pMat;

	data->get(bsData,pMat);

	DataAxis<T,R> const *bs=bsData[0];
	DataAxis<T,R> const *dhbs=bsData[1];
	DataAxis<T,R> const *dh2bs=bsData[2];
	DataAxis<T,R> const *dvbs=bsData[3];
	DataAxis<T,R> const *dv2bs=bsData[4];

	hsize_t dims[2];
	bs->alpha->getDims(dims);
	lli row = dims[0];
	lli col = dims[1];

	lli falsePeak=0;
	lli falseWidth=0;
	lli falseInnerPeak=0;
	// Find Peaks and exact MZ values.
	cout<<"Extract Peaks from Mass Spec Data"<<endl;

	int nthrd=omp_get_num_threads();
	int run;
	#pragma omp parallel
	{
		int thrdid=omp_get_thread_num();

		VecMat<T> csPat(7,7);
		vector<int> offset(2,0);

		csPat.getDims(dims);
		lli pr = dims[0];
		lli pc = dims[1];

		run = 0;
		#pragma omp for reduction(+:falsePeak,falseWidth,falseInnerPeak) schedule(dynamic)
		for(lli i = 2; i < row-3; ++i)
		{
			if(thrdid == 0)
				cout<<"\r"<<"Processing Scan: "<<row-5<<"/"<<run<<flush;

			for(lli j = 3; j < col-4; ++j)
			{
				if( (dhbs->alpha->m[i][j] > 0) && (dhbs->alpha->m[i][j+1] < 0) &&
					(dvbs->alpha->m[i][j] > 0) && (dvbs->alpha->m[i+1][j] < 0) &&
					bs->alpha->m[i][j] > threshold)
				{
					vector<DataPoint> spt(pr); // Surface Points
					vector<DataPoint> p(4); // Surface Points

					double rtPeak=-1.0;
					T countMaxRT=-1.0;
					double t0rt=-1.0;
					double rtrhs = -1.0;
					double rtlhs = -1.0;
					double mzrhs = -1.0;
					double mzlhs = -1.0;

					mulVecMat(*(bs->alpha),pMat->b,csPat,i,j);

					SMData2D<OpInterface> A(dims,&offset[0],1,1,csPat.v);
					SMData2D<OpNablaHS> dhpA(dims,&offset[0],1,1,csPat.v);

					for(int k=0 ; k < pc; ++k)
					{
						A.mz[k]=bs->mz[j-3+k];
						dhpA.mz[k]=dhbs->mz[j-3+k];
					}

					double drt=(bs->rt[i]-bs->rt[i-1])/3.0;
					double ddrt=(dhbs->rt[i]-dhbs->rt[i-1])/3.0;
					for(int k=0 ; k < pr; ++k)
					{
						A.rt[k]=bs->rt[i-1]+k*drt;
						dhpA.rt[k]=dhbs->rt[i-1]+k*ddrt;
						spt[k].rt=A.rt[k];
					}

					vector<int> u(7);
					u[0]=calInnerPeakMZ(&A, &dhpA, 0, 3, spt[0].mz, spt[0].pk, spt[0].t0, falseInnerPeak);
					u[1]=calInnerPeakMZ(&A, &dhpA, 1, 3, spt[1].mz, spt[1].pk, spt[1].t0, falseInnerPeak);
					u[2]=calInnerPeakMZ(&A, &dhpA, 2, 3, spt[2].mz, spt[2].pk, spt[2].t0, falseInnerPeak);
					u[3]=calInnerPeakMZ(&A, &dhpA, 3, 3, spt[3].mz, spt[3].pk, spt[3].t0, falseInnerPeak);

					p[0].t0=-1;
					vector<DataPoint> ospt(4);
					for(int v=0; v <4; ++v)
					{
						if(v > 0)
						{
							int vidx=v+3;
							u[vidx]=calInnerPeakMZ(&A,&dhpA,vidx,3,spt[vidx].mz,spt[vidx].pk,spt[vidx].t0,falseInnerPeak);
						}
						if(u[v]+u[v+1]+u[v+2]+u[v+3] == 4)
						{
							ospt[0]=spt[v];
							ospt[1]=spt[v+1];
							ospt[2]=spt[v+2];
							ospt[3]=spt[v+3];
							bezierP(ospt,p);
							bezierT(p);
						}
						if(p[0].t0 != -1) break;
					}

					if(	p[0].mz < bs->mz[0] || p[0].mz > bs->mz[col-1] ||
						p[1].mz < bs->mz[0] || p[1].mz > bs->mz[col-1] ||
						p[2].mz < bs->mz[0] || p[2].mz > bs->mz[col-1] ||
						p[3].mz < bs->mz[0] || p[3].mz > bs->mz[col-1]) p[0].t0=-1;

					if(p[0].t0 >= 0)
					{
						this->calPeakLenRT(i,j,dv2bs->alpha->m,dv2bs->rt,rtlhs,rtrhs);
						this->calPeakWidth(i,j,dh2bs->alpha->m,dh2bs->mz,mzlhs,mzrhs);
						if(rtlhs > 0 && rtrhs > 0 && mzlhs > bs->mz[0] && mzrhs > bs->mz[0])
						{
							DataPoint B = bezierO3(p);
							if(rtlhs < B.rt && B.rt < rtrhs &&
							   mzlhs < bs->mz[col-1] && mzrhs < bs->mz[col-1] &&
							   B.pk >= threshold)
							{
								#pragma omp critical(peak)
								{
									peak->addPeak(B.mz,B.rt,B.pk,make_pair(mzlhs,mzrhs),
											make_pair(rtlhs,rtrhs),B.t0,j,i);
								}
							}
							else
							{
								++falsePeak;
							}
						}
						else
						{
							++falseWidth;
						}
					}
					else
					{
						++falsePeak;
					}
				}
			}
			#pragma omp atomic
				++run;
		}
	}
	cout<<"\r"<<"Processing Scan: "<<row-5<<"/"<<row-5<<endl;
	peak->updateFalseData(falsePeak,falseWidth);
}

template<template<class> class pPeak, template<class,class> class pData, typename T, typename R>
void ExtractPeak<pPeak,pData,T,R>::mulVecMat(const VecMat<T> &cs, const VecMat<T> &bp, VecMat<T> &pcs, lli i, lli j)
{
	// (0,0) at patch
	lli rdx=i-2;
	lli cdx=j-3;

	lli cols = 7;
	lli rows = 5;

	/*
	hsize_t dim[2];
	cout<<"CS Peak at ("<<i<<","<<j<<")"<<endl;
	cout<<"CS Patch from ("<<rdx<<","<<cdx<<")"<<endl;
	for(int i=rdx; i < rdx+rows; ++i){
		for(int j=cdx; j < cdx+cols; ++j){
			cout<<setw(12);
			cout<<cs.m[i][j]<<"  ";
		}
		cout<<endl;
	}

	cout<<"Result Patch"<<endl;
	pcs.getDims(dim);
	for(int i=0; i < dim[0]; ++i){
		for(int j=0; j < dim[1]; ++j){
			cout<<setw(12);
			cout<<pcs.m[i][j]<<"  ";
		}
		cout<<endl;
	}

	cout<<"B-Spline"<<endl;
	for(int i=0; i < rows; ++i){
		for(int j=0; j < cols; ++j){
			cout<<setw(12);
			cout<<bp.m[i][j]<<"  ";
		}
		cout<<endl;
	}
	*/

	lli rpcs=0;
	lli cpcs=0;
	for(int n = 0; n < cols; ++n)
	{
		rpcs=0;
		for(int m = 0; m < 7; ++m)
		{
			pcs.m[rpcs][cpcs]=
					cs.m[rdx][cdx]  *bp.m[0][m]+
					cs.m[rdx+1][cdx]*bp.m[1][m]+
					cs.m[rdx+2][cdx]*bp.m[2][m]+
					cs.m[rdx+3][cdx]*bp.m[3][m]+
					cs.m[rdx+4][cdx]*bp.m[4][m];
			++rpcs;
		}
		++cdx;
		++cpcs;
	}

	/*
	cout<<"Result Patch"<<endl;
	for(int i=0; i < dim[0]; ++i){
		for(int j=0; j < dim[1]; ++j){
			cout<<setw(12);
			cout<<pcs.m[i][j]<<"  ";
		}
		cout<<endl;
	}
	*/
}

template<template<class> class pPeak, template<class,class> class pData, typename T, typename R>
int ExtractPeak<pPeak,pData,T,R>::calInnerPeakMZ(
		DataAxis<T,R> const *bs,DataAxis<T,R> const *dhbs,
		lli i, lli j, double &mzPeak, T &countMax,
		double &t0, lli &falsePeak)
{
	if((dhbs->alpha->m[i][j] > 0) && (dhbs->alpha->m[i][j+1] < 0) )
	{
		lli u=i;
		lli v=j;
		double pa1=0.0;
		double pmz1=0.0;
		MathOp<T,R>::calMidPoint(u,v,dhbs->alpha->m,dhbs->mz,pmz1,pa1);
		if(pa1 < 0){
			double pa0=0.0;
			double pmz0=0.0;
			vector<T> ry;
			MathOp<T,R>::calMidPoint(u,v-1,dhbs->alpha->m,dhbs->mz,pmz0,pa0);
			t0= MathOp<T,R>::calT(pa0,double(dhbs->alpha->m[u][v]),pa1);
			if(t0>=0)
			{
				mzPeak=MathOp<T,R>::calX(t0,pmz0,dhbs->mz[v],pmz1);
				ry = MathOp<T,R>::cal3rdMidPoint(u,v-1,bs->alpha->m);
				countMax = MathOp<T,R>::calPeakCount(ry,t0);
				return 1;
			}
			else
			{
				++falsePeak;
				return 0;
			}
		}
		else
		{
			double pa2=0.0;
			double pmz2=0.0;
			vector<T> ry;
			MathOp<T,R>::calMidPoint(u,v+1,dhbs->alpha->m,dhbs->mz,pmz2,pa2);
			t0 = MathOp<T,R>::calT(pa1,double(dhbs->alpha->m[u][v+1]),pa2);
			if (t0>=0)
			{
				mzPeak= MathOp<T,R>::calX(t0,pmz1,dhbs->mz[v+1],pmz2);
				ry = MathOp<T,R>::cal3rdMidPoint(u,v,bs->alpha->m);
				countMax = MathOp<T,R>::calPeakCount(ry,t0);
				return 1;
			}
			else
			{
				++falsePeak;
				return 0;
			}
		}
	}
	else if((dhbs->alpha->m[i][j-1] > 0) && (dhbs->alpha->m[i][j] < 0) )
	{
		lli u=i;
		lli v=j-1;
		double pa1=0.0;
		double pmz1=0.0;
		MathOp<T,R>::calMidPoint(u,v,dhbs->alpha->m,dhbs->mz,pmz1,pa1);
		if(pa1 < 0){
			double pa0=0.0;
			double pmz0=0.0;
			vector<T> ry;
			MathOp<T,R>::calMidPoint(u,v-1,dhbs->alpha->m,dhbs->mz,pmz0,pa0);
			t0= MathOp<T,R>::calT(pa0,double(dhbs->alpha->m[u][v]),pa1);
			if(t0>=0)
			{
				mzPeak=MathOp<T,R>::calX(t0,pmz0,dhbs->mz[v],pmz1);
				ry = MathOp<T,R>::cal3rdMidPoint(u,v-1,bs->alpha->m);
				countMax = MathOp<T,R>::calPeakCount(ry,t0);
				return 1;
			}
			else
			{
				++falsePeak;
				return 0;
			}
		}
		else
		{
			double pa2=0.0;
			double pmz2=0.0;
			vector<T> ry;
			MathOp<T,R>::calMidPoint(u,v+1,dhbs->alpha->m,dhbs->mz,pmz2,pa2);
			t0 = MathOp<T,R>::calT(pa1,double(dhbs->alpha->m[u][v+1]),pa2);
			if (t0>=0)
			{
				mzPeak= MathOp<T,R>::calX(t0,pmz1,dhbs->mz[v+1],pmz2);
				ry = MathOp<T,R>::cal3rdMidPoint(u,v,bs->alpha->m);
				countMax = MathOp<T,R>::calPeakCount(ry,t0);
				return 1;
			}
			else
			{
				++falsePeak;
				return 0;
			}
		}
	}
	else
	{
		return 0;
	}
}

template<template<class> class pPeak, template<class,class> class pData, typename T, typename R>
void ExtractPeak<pPeak,pData,T,R>::bezierP(vector<DataPoint> const &spt, vector<DataPoint> &p)
{
	p[0].mz=spt[0].mz;
	p[0].rt=spt[0].rt;
	p[0].pk=spt[0].pk;

	p[3].mz=spt[3].mz;
	p[3].rt=spt[3].rt;
	p[3].pk=spt[3].pk;

	p[1].mz=-(5.0/6.0)*p[0].mz + p[3].mz/3.0 + 3.0*spt[1].mz - (3.0/2.0)*spt[2].mz;
	p[1].rt=-(5.0/6.0)*p[0].rt + p[3].rt/3.0 + 3.0*spt[1].rt - (3.0/2.0)*spt[2].rt;
	p[1].pk=-(5.0/6.0)*p[0].pk + p[3].pk/3.0 + 3.0*spt[1].pk - (3.0/2.0)*spt[2].pk;

	p[2].mz=p[0].mz/3.0 - (5.0/6.0)*p[3].mz - (3.0/2.0)*spt[1].mz + 3.0*spt[2].mz;
	p[2].rt=p[0].rt/3.0 - (5.0/6.0)*p[3].rt - (3.0/2.0)*spt[1].rt + 3.0*spt[2].rt;
	p[2].pk=p[0].pk/3.0 - (5.0/6.0)*p[3].pk - (3.0/2.0)*spt[1].pk + 3.0*spt[2].pk;
}

template<template<class> class pPeak, template<class,class> class pData, typename T, typename R>
void ExtractPeak<pPeak,pData,T,R>::bezierT(vector<DataPoint> &p)
{
	double tn=0.0;
	double tp=0.0;

	T p0=p[0].pk;
	T p1=p[1].pk;
	T p2=p[2].pk;
	T p3=p[3].pk;

	T c0=p0-2.0*p1+p2;
	T c1=-p0*p2+p0*p3+p1*p1-p1*p2-p1*p3+p2*p2;
	T c2=p0-3.0*p1+3.0*p2-p3;

	tn=(c0-sqrt(c1))/c2;
	tp=(c0+sqrt(c1))/c2;

	if(tp >= 0 && tp <= 1) p[0].t0 = tp;
	else if(tn >= 0 && tn <= 1)	p[0].t0 = tn;
	else p[0].t0 = -1.0;
}

template<template<class> class pPeak, template<class,class> class pData, typename T, typename R>
typename ExtractPeak<pPeak,pData,T,R>::DataPoint ExtractPeak<pPeak,pData,T,R>::bezierO3(vector<DataPoint> &p)
{
	double t = p[0].t0;
	double tp = (1-t);
	DataPoint b;
	double c0=tp*tp*tp;
	double c1=3.0*tp*tp*t;
	double c2=3.0*tp*t*t;
	double c3=t*t*t;

	b.mz= c0*p[0].mz + c1*p[1].mz + c2*p[2].mz + c3*p[3].mz;
	b.rt= c0*p[0].rt + c1*p[1].rt + c2*p[2].rt + c3*p[3].rt;
	b.pk= c0*p[0].pk + c1*p[1].pk + c2*p[2].pk + c3*p[3].pk;
	b.t0 = t;

	return b;
}

template<template<class> class pPeak, template<class,class> class pData, typename T, typename R>
void ExtractPeak<pPeak,pData,T,R>::calMidPointRT(lli rtIdx, lli mzIdx, T** alpha,
		const vector<double> &rta, double &rt1, double &a1)
{
	rt1=0.5*(rta[rtIdx]+rta[rtIdx+1]);
	a1=double(0.5*(alpha[rtIdx][mzIdx]+alpha[rtIdx+1][mzIdx]));
}

template<template<class> class pPeak, template<class,class> class pData, typename T, typename R>
vector<T> ExtractPeak<pPeak,pData,T,R>::cal3rdMidPointRT(lli row, lli col, T **P)
{
	vector<T> ry(4,0.0);
	T p0=P[row-1][col];
	T p1=P[row][col];
	T p2=P[row+1][col];
	T p3=P[row+2][col];
	T u=1.0/3.0;
	T v=2.0/3.0;
	T w=1.0/2.0;
	T up=v; // (1-u);
	T vp=u; // (1-v);
	T wp=w; // (1-w);

	ry[0]=wp*(vp*p0+v*p1)+w*(up*p1+u*p2);
	ry[1]=up*p1+u*p2;
	ry[2]=vp*p1+v*p2;
	ry[3]=wp*(vp*p1+v*p2)+w*(up*p2+u*p3);

	return ry;
}

template<template<class> class pPeak, template<class,class> class pData, typename T, typename R>
void ExtractPeak<pPeak,pData,T,R>::calPeakWidthRT(lli rtIdx,lli mzIdx, T** alpha, const vector<double> d2rt,
		double &rtlhs, double &rtrhs)
{
	lli n=d2rt.size();
	lli lhsIdx=rtIdx+1;
	lli rhsIdx=rtIdx+1;

	do
	{
		--lhsIdx;
		if(lhsIdx < 0) break;
	} while( alpha[lhsIdx][mzIdx] < 0);
	do
	{
		++rhsIdx;
		if(rhsIdx > n) break;
	} while( alpha[rhsIdx][mzIdx] < 0);

	if(lhsIdx > 0 && rhsIdx < n ){
		double x1,x2,y1,y2,m;

		x1=d2rt[lhsIdx];
		y1=alpha[lhsIdx][mzIdx];
		x2=d2rt[lhsIdx+1];
		y2=alpha[lhsIdx+1][mzIdx];
		m=(y2-y1)/(x2-x1);
		rtlhs=(-y1/m)+x1;

		x1=d2rt[rhsIdx];
		y1=alpha[rhsIdx][mzIdx];
		x2=d2rt[rhsIdx-1];
		y2=alpha[rhsIdx-1][mzIdx];
		m=(y2-y1)/(x2-x1);
		rtrhs=(-y1/m)+x1;

		/*
		if (mzlhs < d2mz[0])
		{
			cout<<"WTF LHS is out of range!!"<<endl;
		}
		if (mzrhs > d2mz[n-1])
		{
			cout<<"WTF RHS is out of range!!"<<endl;
		}
		*/
	}
	else
	{
		rtlhs=-1;
		rtrhs=-1;
	}
}



template<template<class> class pPeak, template<class,class> class pData, typename T, typename R>
void ExtractPeak<pPeak,pData,T,R>::calPeakRT(
		DataAxis<T,R> const *bs,DataAxis<T,R> const *dbs,DataAxis<T,R> const *d2bs,
		lli i, lli j, double &mzPeak, T &countMax, double &rtlhs, double &rtrhs,
		double &t0, lli &falsePeak)
{
	if((dbs->alpha->m[i][j] > 0) && (dbs->alpha->m[i+1][j] < 0) )
	{
		double pa1=0.0;
		double prt1=0.0;
		calMidPointRT(i,j,dbs->alpha->m,dbs->rt,prt1,pa1);
		if(pa1 < 0){
			double pa0=0.0;
			double prt0=0.0;
			vector<T> ry;
			calMidPointRT(i-1,j,dbs->alpha->m,dbs->rt,prt0,pa0);
			t0=MathOp<T,R>::calT(pa0,double(dbs->alpha->m[i][j]),pa1);
			if(t0>=0)
			{
				mzPeak=MathOp<T,R>::calX(t0,prt0,dbs->rt[j],prt1);
				rtlhs=0.0;
				rtrhs=0.0;
				ry = cal3rdMidPointRT(i-1,j,bs->alpha->m);
				countMax = MathOp<T,R>::calPeakCount(ry,t0);
				calPeakWidthRT(i,j,d2bs->alpha->m,d2bs->rt,rtlhs,rtrhs);
			}
			else
			{
				++falsePeak;
			}
		}
		else
		{
			double pa2=0.0;
			double pmz2=0.0;
			vector<T> ry;
			calMidPointRT(i+1,j,dbs->alpha->m,dbs->rt,pmz2,pa2);
			t0 = MathOp<T,R>::calT(pa1,double(dbs->alpha->m[i+1][j]),pa2);
			if (t0>=0)
			{
				mzPeak= MathOp<T,R>::calX(t0,prt1,dbs->rt[j+1],pmz2);
				rtlhs=0.0;
				rtrhs=0.0;
				ry = cal3rdMidPointRT(i,j,bs->alpha->m);
				countMax = MathOp<T,R>::calPeakCount(ry,t0);
				calPeakWidthRT(i,j,d2bs->alpha->m,d2bs->rt,rtlhs,rtrhs);
			}
			else
			{
				++falsePeak;
			}
		}
	}
}

template<template<class> class pPeak, template<class,class> class pData, typename T, typename R>
void ExtractPeak<pPeak,pData,T,R>::calPeakLenRT(lli rtIdx,lli mzIdx, T** alpha,
		const vector<double> d2rt, double &rtlhs, double &rtrhs)
{
	lli n=d2rt.size()-1;
	lli lhsIdx=rtIdx+1;
	lli rhsIdx=rtIdx+1;

	do
	{
		--lhsIdx;
		if(lhsIdx == 0) break;
	} while(alpha[lhsIdx][mzIdx] < alpha[lhsIdx-1][mzIdx]);
	do
	{
		++rhsIdx;
		if(rhsIdx == n) break;
	} while(alpha[rhsIdx][mzIdx] < alpha[rhsIdx+1][mzIdx]);

	if(lhsIdx >= 0)
	{
		rtlhs=d2rt[lhsIdx];
	}
	else
	{
		rtlhs=-1;
	}
	if(rhsIdx <= n)
	{
		rtrhs=d2rt[rhsIdx];
	}
	else
	{

		rtrhs=-1;
	}

	/*
	if (mzlhs < d2mz[0])
	{
		cout<<"WTF LHS is out of range!!"<<endl;
	}
	if (mzrhs > d2mz[n-1])
	{
		cout<<"WTF RHS is out of range!!"<<endl;
	}
	*/
}

#endif /* SMPEAK_PEAKOPERATOR_TPP_ */
