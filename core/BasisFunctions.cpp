//
// $Id$
//
//
// Original author: Andrew Dowsey <andrew.dowsey <a.t> bristol.ac.uk>
//
// Copyright (C) 2016  biospi Laboratory, University of Bristol, UK
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


#include "BasisFunctions.hpp"
#include "BSpline.hpp"
#include <iostream>
#include <limits>
#include <iomanip>


using namespace std;


CoeffsMetadata::
CoeffsMetadata(ii _dc) :
d(_dc),
l(_dc),
o(_dc),
n(_dc)
{
}


void
CoeffsMetadata::
operator=(const CoeffsMetadata& cm)
{
	d = cm.d;
	l = cm.l;
	o = cm.o;
	n = cm.n;
}


CoeffsMetadata::
~CoeffsMetadata()
{
}


li
CoeffsMetadata::
size() const
{
	li size = 1;
	for (ii i = 0; i < d; i++) size *= n[i];
	return size;
}


void
CoeffsMetadata::
print(ostream& out) const
{
	out << "lc=[";
	for (ii i = 0; i < d; i++)
	{
		if (l[i] == numeric_limits<ii>::min()) out << "?"; else out << l[i];
		if (i < d - 1) out << ",";
	}
	out << "] oc=[";
	for (ii i = 0; i < d; i++)
	{
		if (o[i] == numeric_limits<ii>::min()) out << "?"; else out << o[i];
		if (i < d - 1) out << ",";
	}
	out << "] nc=[";
	for (ii i = 0; i < d; i++)
	{
		if (n[i] == numeric_limits<ii>::min()) out << "?"; else out << n[i];
		if (i < d - 1) out << ",";
	}
	out << "]:" << size();
}


////////////////////////////////////////////////////////////////////////////////

Basis::
Basis(vector<Basis*>& bases, ii dimensions, Basis* _parent, bool _transient) :
    parent(_parent),
    transient(_transient),
    cm(dimensions),
    volume(0.0),
    discrep(0.0)
{
    index = bases.size();
    bases.push_back(this);
    if (parent) parent->child_count++;
}


void
Basis::
error(vector<fp>& fs, const vector<fp>& gs)
{
    ii size_d = 0;
    double dis = 0.0;
    double err = 0.0;
    double sum_g = 0.0;
    double sum_f = 0.0;
    maxerr = 0;
    for (li i = 0; i < gs.size(); i++)
    {
        double v = fabs(gs[i]-fs[i]);
        maxerr = maxerr > v ? maxerr : v;
    }
    // bug in this openmp section at present for err and discrep
    //#pragma omp parallel for simd reduction(+:dis,err,size_d,sum_g,sum_f)
    for (li i = 0; i < gs.size(); i++)
    {
        sum_g += gs[i];
        sum_f += fs[i];
        
        if (fs[i] > 0.0 && gs[i] >= 0.0)
        {
			if (gs[i] > 0.0)
			{
				dis += ((gs[i] - ceil(fs[i]))*(gs[i] - ceil(fs[i]))) / ceil(fs[i]);
				err += fabs(gs[i] - fs[i]);
				size_d++;
			}
            
            fs[i] = gs[i]/fs[i];
        }
        else
        {
            fs[i] = 1.0;
        }
    }
    volume = sum_f / sum_g;
    discrep = dis / size_d;
    erro = err / size_d;
}


////////////////////////////////////////////////////////////////////////////////


BasisResampleMZ::
BasisResampleMZ(std::vector<Basis*>& bases, const std::vector<fp>& bin_counts, const std::vector<li>& spectrum_index, const std::vector<double>& bin_edges, short resolution, ii order, bool transient) :
    Basis(bases, 2, 0, transient)
{
	if (spectrum_index.size() > 0)
	{
		is = spectrum_index;
	}
	else
	{
		is.push_back(0);
		is.push_back(bin_counts.size());
	}

    ///////////////////////////////////////////////////////////////////////
    // create A as a temporary COO matrix
    
    // calculate indicies of non-empty spectra
    cm.l[1] = numeric_limits<ii>::min();
    cm.o[1] = numeric_limits<ii>::min();
    cm.n[1] = is.size() - 1;

    // init arrays
    a.resize(cm.n[1]);
    ia.resize(cm.n[1]);
    ja.resize(cm.n[1]);
    m.resize(cm.n[1]);
    nnz.resize(cm.n[1], 0);

    // find min and max m/z across spectra, and narrowest bin
    mz_min = numeric_limits<double>::max();
    mz_max = 0.0;
	double mz_diff = numeric_limits<double>::max();
    for (ii j = 0; j < cm.n[1]; j++)
    {
        m[j] = (ii) (is[j+1] - is[j]);
		mz_min = bin_edges[is[j] + j] < mz_min ? bin_edges[is[j] + j] : mz_min;
		mz_max = bin_edges[is[j + 1] + j] > mz_max ? bin_edges[is[j + 1] + j] : mz_max;

		for (ii i = 0; i < m[j]; i++)
		{
			double diff = bin_edges[is[j] + j + i + 1] - bin_edges[is[j] + j + i];
			mz_diff = diff < mz_diff ? diff : mz_diff;
		}
	}

	ii resolution_auto = floor(log2(1.0 / mz_diff / 60.0 / 1.0033548378));
	if (resolution == numeric_limits<short>::max())
	{
		resolution = resolution_auto;
	}
	// Bases per 1.0033548378Th (difference between carbon12 and carbon13)
	double bpi = pow(2.0, (double)resolution) * 60 / 1.0033548378;

	cm.l[0] = resolution;
    cm.o[0] = (ii) floor(mz_min * bpi);
    cm.n[0] = ((ii) ceil(mz_max * bpi)) + order - cm.o[0];

    // figure out nnz
    #pragma omp parallel for
	for (ii j = 0; j < cm.n[1]; ++j)
	{
		for (ii i = 0; i < m[j]; i++)
		{
			if (bin_counts[is[j] + i] >= 0.0)
			{
				double cf0 = bin_edges[is[j] + j + i] * bpi;
				double cf1 = bin_edges[is[j] + j + i + 1] * bpi;

				ii ci0 = (ii)floor(cf0);
				ii ci1 = ((ii)ceil(cf1)) + order;

				nnz[j] += ci1 - ci0;
			}
		}
	}

    // populate coo matrix
    // should also implement b-splines as a lookup table, possibly faster
    ii done = 0;
    #pragma omp parallel for
    for(ii j = 0; j < cm.n[1]; ++j)
    {
        vector<fp> acoo(nnz[j]);
        vector<ii> rowind(nnz[j]);
        vector<ii> colind(nnz[j]);
        
        ii k = 0;
        for (ii i = 0; i < m[j]; i++)
        if (bin_counts[is[j]+i] >= 0.0)
        {
			double cf0 = bin_edges[is[j] + j + i] * bpi;
			double cf1 = bin_edges[is[j] + j + i + 1] * bpi;

            ii ci0 = (ii) floor(cf0);
            ii ci1 = ((ii) ceil(cf1)) + order;
            
            // work out basis coefficients
            for (ii ci = ci0; ci < ci1; ci++)
            {
                double bf0 = ci - order;
                double bf1 = ci + 1;
                
                // intersection of bin and basis, between 0 and order+1
                double b0 = cf0 > bf0 ? cf0 - bf0 : 0.0;
                double b1 = cf1 < bf1 ? cf1 - bf0 : bf1 - bf0;
                
                // basis coefficient b is _integral_ of area under b-spline basis
				double b = BSpline::im(b1, order + 1) - BSpline::im(b0, order + 1);
                if (b <= 0.000001) b = 0.0; // saves computation and problems with l2norm later
                
                acoo[k] = b;
                rowind[k] = i;
                colind[k] = ci - cm.o[0];
                k++;
            }
        }
        
        // create A and free coo
        a[j].resize(nnz[j]);
        ia[j].resize(m[j]+1);
        ja[j].resize(nnz[j]);
        
        ii job[] = {2, 0, 0, 0, nnz[j], 0}; ii info;
        mkl_scsrcoo(job, &m[j], a[j].data(), ja[j].data(), ia[j].data(), &(nnz[j]), acoo.data(), rowind.data(), colind.data(), &info);
        
        // display progress update
        #pragma omp critical
        {
            done++;
			if (done % 100 == 0)
			{
				for (int i = 0; i < 256; ++i) cout << '\b';
				cout << index << " BasisResampleMZ " << setw(1+(int)(log10((float)cm.n[1]))) << done << "/" << cm.n[1] << " " << flush;
			}
        }
    }
    for (int i = 0; i < 256; ++i) cout << '\b';
    
    cout << index << " BasisResampleMZ ";
    cm.print(cout);
    li size = is.back();
    for (ii j = 0; j < a.size(); j++) size += 2*nnz[j]+1;
    cout << " mem=" << setprecision(2) << fixed << (sizeof(this)+size*sizeof(fp))/(1024.0*1024.0) << "Mb";
    if (transient) cout << " (t)";
    cout << endl;
	cout << "   resolution=" << resolution << " (" << bpi << " bases per 1.0033548378Th)" << " range=" << setprecision(3) << mz_min << ":" << mz_diff << ":" << mz_max << "Th" << endl;
	if (resolution_auto != resolution)
	{
		cerr << endl << "WARNING: resolution is not the suggested value of " << resolution_auto << ". Continue at your own risk!" << endl << endl;
	}
}


BasisResampleMZ::~BasisResampleMZ()
{
}


void
BasisResampleMZ::
synthesis(vector<fp>& fs, const vector<fp>& cs, bool accum)
{
    static fp alpha = 1.0;
    fp beta = accum ? 1.0 : 0.0;
    # pragma omp parallel for
    for (li j = 0; j < cm.n[1]; j++)
    {
        //cout << "synthesis " << j << "," << cm.n[0] << "," << cm.n[1] << "," << m[j] << endl;
        fp* c = const_cast<fp*>(&(cs.data()[j*cm.n[0]]));
        mkl_scsrmv("N", &(m[j]), &(cm.n[0]), &alpha, "G**C", a[j].data(), ja[j].data(), ia[j].data(), &(ia[j].data()[1]), c, &beta, &(fs.data()[is[j]]));
    }
}


void
BasisResampleMZ::
analysis(vector<fp>& es, const vector<fp>& fs)
{
    static fp alpha = 1.0, beta = 0.0;
    //# pragma omp parallel for
    for (li j = 0; j < cm.n[1]; j++)
    {
        //cout << "analysis " << j << "," << cm.n[0] << "," << cm.n[1] << "," << m[j] << ":" << j*cm.n[0] << ":" << j*cm.n[0] << endl;
        fp* f = const_cast<fp*>(&(fs.data()[is[j]]));
        //mkl_scsrmv("N", &(cm.n[0]), &(m[j]), &alpha, "G**C", at[j].data(), jat[j].data(), iat[j].data(), &(iat[j].data()[1]), f, &beta, &(es.data()[j*cm.n[0]]));
        mkl_scsrmv("T", &(m[j]), &(cm.n[0]), &alpha, "G**C", a[j].data(), ja[j].data(), ia[j].data(), &(ia[j].data()[1]), f, &beta, &(es.data()[j*cm.n[0]]));
    }
}


void
BasisResampleMZ::
l2norm(vector<fp>& es, const vector<fp>& fs)
{
    static fp alpha = 1.0, beta = 0.0;
    //# pragma omp parallel for
    for (li j = 0; j < cm.n[1]; j++)
    {
        //cout << "l2norm " << j << "," << cm.n[0] << "," << cm.n[1] << "," << m[j] << "," << nnz[j] << endl;
        vector<fp> a2(nnz[j]);
        vsSqr(nnz[j], a[j].data(), a2.data());
        fp* f = const_cast<fp*>(&(fs.data()[is[j]]));
        mkl_scsrmv("T", &(m[j]), &(cm.n[0]), &alpha, "G**C", a2.data(), ja[j].data(), ia[j].data(), &(ia[j].data()[1]), f, &beta, &(es.data()[j*cm.n[0]]));
    }
}


////////////////////////////////////////////////////////////////////////////////


BasisResampleRT::
BasisResampleRT(std::vector<Basis*>& bases, Basis* parent, const std::vector<double>& start_times, const std::vector<double>& finish_times, const std::vector<fp>& exposures, short resolution, ii order, bool transient) :
    Basis(bases, parent->get_cm().d, parent, transient)
{
	rt_min = start_times.front();
	rt_max = finish_times.back();

	double rt_diff = numeric_limits<double>::max();
	for (ii j = 0; j < exposures.size(); j++)
	{
		double diff = finish_times[j] - start_times[j];
		rt_diff = diff < rt_diff ? diff : rt_diff;
	}
	ii resolution_auto = floor(log2(1.0 / rt_diff));
	if (resolution == numeric_limits<short>::max())
	{
		resolution = resolution_auto;
	}
	double bpi = pow(2.0, (double) resolution);
    double spacing = 1.0 / bpi;

    
    // create A as a temporary COO matrix
    m = exposures.size();
    nnz = 0;
    for (ii j = 0; j < exposures.size(); j++)
    {
		double cf0 = start_times[j] / spacing;
		double cf1 = finish_times[j] / spacing;
        
        ii ci0 = (ii) floor(cf0);
        ii ci1 = ((ii) ceil(cf1)) + order;
        
        nnz += ci1 - ci0;
    }
    
    cm.l[0] = parent->get_cm().l[0];
    cm.o[0] = parent->get_cm().o[0];
    cm.n[0] = parent->get_cm().n[0];
    cm.l[1] = resolution;
    cm.o[1] = (ii) floor(start_times.front() / spacing);
    cm.n[1] = ((ii) ceil(finish_times.back() / spacing)) + order - cm.o[1];

    // populate coo matrix
    vector<fp> acoo(nnz);
    vector<ii> rowind(nnz);
    vector<ii> colind(nnz);
    
    // create A as a temporary COO matrix
    ii mi = 0;
    for (ii i = 0, j = 0; j < exposures.size(); j++)
    {
		double cf0 = start_times[j] / spacing;
		double cf1 = finish_times[j] / spacing;

        ii ci0 = (ii) floor(cf0);
        ii ci1 = ((ii) ceil(cf1)) + order;
        
        for (int ci = ci0; ci < ci1; ci++)
        {
            double bf0 = ci - order;
            double bf1 = ci + 1;
            
            // intersection of bin and basis, between 0 and order+1
            double b0 = cf0 > bf0 ? cf0 - bf0 : 0.0;
            double b1 = cf1 < bf1 ? cf1 - bf0 : bf1 - bf0;
            
            // basis coefficient b is _integral_ of area under b-spline basis
			double b = BSpline::im(b1, order + 1) - BSpline::im(b0, order + 1);
            if (b <= FLT_MIN) b = FLT_MIN; // for numerical instability in im()
            
			acoo[mi] = b * exposures[j];
            rowind[mi] = i;
            colind[mi] = ci - cm.o[1];
            mi++;
        }
        i++;
    }
    
    // create A and free coo
    a.resize(nnz);
    ia.resize(m+1);
    ja.resize(nnz);
    at.resize(nnz);
    iat.resize(cm.n[1]+1);
    jat.resize(nnz);
    
    ii job[] = {2, 0, 0, 0, nnz, 0}; ii info;
    mkl_scsrcoo(job, &m, a.data(), ja.data(), ia.data(), &nnz, acoo.data(), rowind.data(), colind.data(), &info);
    mkl_scsrcoo(job, &cm.n[1], at.data(), jat.data(), iat.data(), &nnz, acoo.data(), colind.data(), rowind.data(), &info);
    
    cout << index << " BasisResampleRT parent=" << parent->get_index() << " ";
    cm.print(cout);
    cout << " A=[" << m << "," << cm.n[1] << "]:" << nnz;
    if (transient) cout << " (t)";
    cout << endl;
	cout << "   resolution=" << resolution << " (" << bpi << " bases per second)" << " range=" << setprecision(3) << rt_min << ":" << rt_diff << ":" << rt_max << "seconds" << endl;
	if (resolution_auto != resolution)
	{
		cerr << endl << "WARNING: resolution is not the suggested value of " << resolution_auto << ". Continue at your own risk!" << endl << endl;
	}
}


BasisResampleRT::~BasisResampleRT()
{
}


void
BasisResampleRT::
synthesis(vector<fp>& fs, const vector<fp>& cs, bool accum)
{
    static fp alpha = 1.0;
    fp beta = accum ? 1.0 : 0.0;
    fp* c = const_cast<fp*>(cs.data());
    mkl_scsrmm("N", &m, &(cm.n[0]), &(cm.n[1]), &alpha, "G**C", a.data(), ja.data(), ia.data(), &(ia.data()[1]), c, &(cm.n[0]),  &beta, fs.data(), &(cm.n[0]));
}


void
BasisResampleRT::
analysis(vector<fp>& es, const vector<fp>& fs)
{
    static fp alpha = 1.0, beta = 0.0;
    fp* f = const_cast<fp*>(fs.data());
    //mkl_scsrmm("T", &m, &(cm.n[0]), &(cm.n[1]), &alpha, "G**C", a.data(), ja.data(), ia.data(), &(ia.data()[1]), f, &(cm.n[0]),  &beta, es.data(), &(cm.n[0]));
    mkl_scsrmm("N", &(cm.n[1]), &(cm.n[0]), &m, &alpha, "G**C", at.data(), jat.data(), iat.data(), &(iat.data()[1]), f, &(cm.n[0]),  &beta, es.data(), &(cm.n[0]));
}


void
BasisResampleRT::
l2norm(vector<fp>& es, const vector<fp>& fs)
{
    vector<fp> a2(nnz);
    vsSqr(nnz, a.data(), a2.data());
    static fp alpha = 1.0, beta = 0.0;
    fp* f = const_cast<fp*>(fs.data());
    mkl_scsrmm("T", &m, &(cm.n[0]), &(cm.n[1]), &alpha, "G**C", a2.data(), ja.data(), ia.data(), &(ia.data()[1]), f, &(cm.n[0]),  &beta, es.data(), &(cm.n[0]));
}


////////////////////////////////////////////////////////////////////////////////


BasisDyadicScale::
BasisDyadicScale(vector<Basis*>& bases,
                 Basis* parent,
                 ii _dim,
                 ii order,
                 bool transient) :
    Basis(bases, parent->get_cm().d, parent, transient),
    dim(_dim)
{
    cm = parent->get_cm();
    cm.l[dim] = parent->get_cm().l[dim] - 1;
    cm.o[dim] = parent->get_cm().o[dim]/2;
    cm.n[dim] = (parent->get_cm().o[dim] + parent->get_cm().n[dim] - 1 - order)/2 + order + 1 - cm.o[dim];
    m = parent->get_cm().n[dim];
    n = cm.n[dim];
    
    ii stride = 1;
    for (ii j = 0; j < dim; j++) stride *= cm.n[j];

    // create our kernel
    ii nh = order + 2;
    vector<fp> hs(nh);
    double sum = 0.0;
    for (ii i = 0; i < nh; i++)
    {
        hs[i] = 1.0/pow(2.0, (double)order) /**
        bspline::factorial(order+1)/(double)(bspline::factorial(i)*bspline::factorial(order+1-i))*/;
        sum += hs[i];
    }
    for (ii i = 0; i < nh; i++)
    {
        hs[i] /= sum;
    }
    
    // create A as a temporary COO matrix
    vector<fp> acoo(nh * n);
    vector<ii> rowind(nh * n);
    vector<ii> colind(nh * n);
    
    nnz = 0;
    ii offset = order + ((parent->get_cm().o[dim] + 1) % 2);
    for (ii j = 0; j < n; j++)
    {
        for (ii i = 0; i < nh; i++)
        {
            rowind[nnz] = 2 * j + i - offset;
            if (rowind[nnz] < 0 || rowind[nnz] >= m) continue;
            acoo[nnz] = hs[i];
            colind[nnz] = j;

            nnz++;
        }
    }
    
    // create A and free coo
    a.resize(nnz);
    ia.resize(m+1);
    ja.resize(nnz);
    at.resize(nnz);
    iat.resize(n+1);
    jat.resize(nnz);
    
    //ii major = (dim == 0) ? 1 : 0;
    ii job[] = {2, 0, 0, 0, nnz, 0}; ii info;
    mkl_scsrcoo(job, &m, a.data(), ja.data(), ia.data(), &nnz, acoo.data(), rowind.data(), colind.data(), &info);
    mkl_scsrcoo(job, &n, at.data(), jat.data(), iat.data(), &nnz, acoo.data(), colind.data(), rowind.data(), &info);

    cout << index << " BasisDyadicScale parent=" << parent->get_index() << " dim=" << dim << " ";
    cm.print(cout);
    cout << " A=[" << m << "," << n << "]:" << nnz;
    if (transient) cout << " (t)";
    cout << endl;
}


BasisDyadicScale::~BasisDyadicScale()
{
}


void
BasisDyadicScale::
synthesis(vector<fp>& fs, const vector<fp>& cs, bool accum)
{
    static fp alpha = 1.0;
    fp beta = accum ? 1.0 : 0.0;
    
    if (dim == 0)
    {
        #pragma omp parallel for
        for (ii i = 0; i < cm.n[1]; i++)
        {
            fp* c = const_cast<fp*>(&(cs.data()[i*cm.n[0]]));
            mkl_scsrmv("N", &m, &n, &alpha, "G**C", a.data(), ja.data(), ia.data(), &(ia.data()[1]), c, &beta, &(fs.data()[i*parent->get_cm().n[0]]));
        }
    }
    else
    {
        fp* c = const_cast<fp*>(cs.data());
        mkl_scsrmm("N", &m, &(cm.n[0]), &n, &alpha, "G**C", a.data(), ja.data(), ia.data(), &(ia.data()[1]), c, &(cm.n[0]),  &beta, fs.data(), &(cm.n[0]));
    }
}


void
BasisDyadicScale::
analysis(vector<fp>& es, const vector<fp>& fs)
{
    static fp alpha = 1.0, beta = 0.0;

    if (dim == 0)
    {
        #pragma omp parallel for
        for (ii i = 0; i < cm.n[1]; i++)
        {
            fp* f = const_cast<fp*>(&(fs.data()[i*parent->get_cm().n[0]]));
            mkl_scsrmv("T", &m, &n, &alpha, "G**C", a.data(), ja.data(), ia.data(), &(ia.data()[1]), f, &beta, &(es.data()[i*cm.n[0]]));
        }
    }
    else
    {
        fp* f = const_cast<fp*>(fs.data());
        //mkl_scsrmm("T", &m, &(cm.n[0]), &n, &alpha, "G**C", a.data(), ja.data(), ia.data(), &(ia.data()[1]), f, &(cm.n[0]),  &beta, es.data(), &(cm.n[0]));
        mkl_scsrmm("N", &n, &(cm.n[0]), &m, &alpha, "G**C", at.data(), jat.data(), iat.data(), &(iat.data()[1]), f, &(cm.n[0]),  &beta, es.data(), &(cm.n[0]));
    }
}


void
BasisDyadicScale::
l2norm(vector<fp>& es, const vector<fp>& fs)
{
    
    static fp alpha = 1.0, beta = 0.0;

    vector<fp> a2(nnz);
    vsSqr(nnz, a.data(), a2.data());
    
    if (dim == 0)
    {
        static ii one = 1;
        for (ii i = 0; i < cm.n[1]; i++)
        {
            fp* f = const_cast<fp*>(&(fs.data()[i*parent->get_cm().n[0]]));
            mkl_scsrmm("T", &m, &one, &n, &alpha, "G**C", a2.data(), ja.data(), ia.data(), &(ia.data()[1]), f, &one,  &beta, &(es.data()[i*cm.n[0]]), &one);
        }
    }
    else
    {
        fp* f = const_cast<fp*>(fs.data());
        mkl_scsrmm("T", &m, &(cm.n[0]), &n, &alpha, "G**C", a2.data(), ja.data(), ia.data(), &(ia.data()[1]), f, &(cm.n[0]),  &beta, es.data(), &(cm.n[0]));
    }
}


