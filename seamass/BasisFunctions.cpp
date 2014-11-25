//
// $Id$
//
//
// Original author: Andrew Dowsey <andrew.dowsey <a.t> manchester.ac.uk>
//
// Copyright (C) 2013  CADET Laboratory for Medical Bioinformatics, University of Manchester, UK
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
#include <limits>


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
    for (ii i = 0; i < gs.size(); i++)
    {
        double v = fabs(gs[i]-fs[i]);
        maxerr = maxerr > v ? maxerr : v;
    }
    // bug in this openmp section at present for err and discrep
    //#pragma omp parallel for simd reduction(+:dis,err,size_d,sum_g,sum_f)
    for (ii i = 0; i < gs.size(); i++)
    {
        sum_g += gs[i];
        sum_f += fs[i];
        
        if (fs[i] > 0.0 && gs[i] >= 0.0)
        {
            dis += ((gs[i]-ceil(fs[i]))*(gs[i]-ceil(fs[i])))/ceil(fs[i]);
            err += fabs(gs[i]-fs[i]);
            size_d++;
            
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
BasisResampleMZ(vector<Basis*>& bases,
                const vector< vector<double> >& mzs,
                const vector<fp>& gs,
                const vector<li>& _is,
                const vector<ii>& js,
                ii rci,
                ii order,
                bool transient) :
    Basis(bases, 2, 0, transient),
    is(_is)
{
	double rc = pow(2.0, (double) rci) * 60 / 1.0033548378;

    ///////////////////////////////////////////////////////////////////////
    // create A as a temporary COO matrix
    
    // calculate indicies of non-empty spectra
    cm.l[1] = numeric_limits<ii>::min();
    cm.o[1] = numeric_limits<ii>::min();
    cm.n[1] = js.size();

    // init arrays
    a.resize(cm.n[1]);
    ia.resize(cm.n[1]);
    ja.resize(cm.n[1]);
    //at.resize(cm.n[1]);
    //iat.resize(cm.n[1]);
    //jat.resize(cm.n[1]);
    m.resize(cm.n[1]);
    nnz.resize(cm.n[1], 0);

    // find min and max m/z across spectra
    double mz_min = DBL_MAX;
    double mz_max = 0.0;
    for (ii j = 0; j < cm.n[1]; j++)
    {
        m[j] = (ii) (is[j+1] - is[j]);
        mz_min = mzs[js[j]].front() < mz_min ? mzs[js[j]].front() : mz_min;
        mz_max = mzs[js[j]].back() > mz_max ? mzs[js[j]].back() : mz_max;
    }
    cm.l[0] = rci;
    cm.o[0] = (ii) floor(mz_min * rc);
    cm.n[0] = ((ii) ceil(mz_max * rc)) + order - cm.o[0];

    // figure out nnz
    #pragma omp parallel for
    for(ii j = 0; j < cm.n[1]; ++j)
    for (ii i = 0; i < m[j]; i++)
    if (gs[is[j]+i] >= 0.0)
    {
        double cf0 = mzs[js[j]][i] * rc;
        double cf1 = mzs[js[j]][i+1] * rc;
        
        ii ci0 = (ii) floor(cf0);
        ii ci1 = ((ii) ceil(cf1)) + order;
        
        nnz[j] += ci1 - ci0;
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
        if (gs[is[j]+i] >= 0.0)
        {
            double cf0 = mzs[js[j]][i] * rc;
            double cf1 = mzs[js[j]][i+1] * rc;
            
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
                double b = bspline::im(b1, order+1) - bspline::im(b0, order+1);
                if (b <= FLT_MIN) b = FLT_MIN; // for numerical instability in im()
                
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
        //at[j].resize(nnz[j]);
        //iat[j].resize(cm.n[0]+1);
        //jat[j].resize(nnz[j]);
        
        ii job[] = {2, 0, 0, 0, nnz[j], 0}; ii info;
        mkl_scsrcoo(job, &m[j], a[j].data(), ja[j].data(), ia[j].data(), &(nnz[j]), acoo.data(), rowind.data(), colind.data(), &info);
        //mkl_scsrcoo(job, &cm.n[0], at[j].data(), jat[j].data(), iat[j].data(), &(nnz[j]), acoo.data(), colind.data(), rowind.data(), &info);
        
        // display progress update
        #pragma omp critical
        {
			if (done % 100 == 0)
			{
				for (int i = 0; i < 256; ++i) cout << '\b';
				cout << index << " BasisResampleMZ " << setw(1+(int)(log10((float)cm.n[1]))) << ++done << "/" << cm.n[1] << " " << flush;
			}
        }
    }
    for (int i = 0; i < 256; ++i) cout << '\b';
    
    cout << index << " BasisResampleMZ ";
    cm.print(cout);
    ii size = is.back();
    for (ii j = 0; j < a.size(); j++) size += 2*nnz[j]+1;
    cout << " mem=" << setprecision(2) << fixed << (sizeof(this)+size*sizeof(fp))/(1024.0*1024.0) << "Mb";
    if (transient) cout << " (t)";
    cout << endl;
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
    for (ii j = 0; j < cm.n[1]; j++)
    {
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
    for (ii j = 0; j < cm.n[1]; j++)
    {
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
    for (ii j = 0; j < cm.n[1]; j++)
    {
        vector<fp> a2(nnz[j]);
        vsSqr(nnz[j], a[j].data(), a2.data());
        fp* f = const_cast<fp*>(&(fs.data()[is[j]]));
        mkl_scsrmv("T", &(m[j]), &(cm.n[0]), &alpha, "G**C", a2.data(), ja[j].data(), ia[j].data(), &(ia[j].data()[1]), f, &beta, &(es.data()[j*cm.n[0]]));
    }
}


////////////////////////////////////////////////////////////////////////////////


BasisResampleRT::
BasisResampleRT(vector<Basis*>& bases,
                Basis* parent,
                const vector<double>& rts,
                const vector<ii>& js,
                ii rci,
                ii order,
                bool transient) :
    Basis(bases, parent->get_cm().d, parent, transient)
{
	double rc = pow(2.0, (double) rci);
    double sp = 1.0 / rc;
    
    // create A as a temporary COO matrix
    m = js.size();
    nnz = 0;
    for (ii j = 0; j < js.size(); j++)
    {
        double cf0 = rts[js[j]] / sp;
        double cf1 = rts[js[j]+1] / sp;
        
        ii ci0 = (ii) floor(cf0);
        ii ci1 = ((ii) ceil(cf1)) + order;
        
        nnz += ci1 - ci0;
    }
    
    cm.l[0] = parent->get_cm().l[0];
    cm.o[0] = parent->get_cm().o[0];
    cm.n[0] = parent->get_cm().n[0];
    cm.l[1] = rci;
    cm.o[1] = (ii) floor(rts.front() / sp);
    cm.n[1] = ((ii) ceil(rts.back() / sp)) + order - cm.o[1];

    // populate coo matrix
    vector<fp> acoo(nnz);
    vector<ii> rowind(nnz);
    vector<ii> colind(nnz);
    
    // create A as a temporary COO matrix
    ii mi = 0;
    for (ii i = 0, j = 0; j < js.size(); j++)
    {
        double cf0 = rts[js[j]] / sp;
        double cf1 = rts[js[j]+1] / sp;
        
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
            double b = bspline::im(b1, order+1) - bspline::im(b0, order+1);
            if (b <= FLT_MIN) b = FLT_MIN; // for numerical instability in im()
            
            acoo[mi] = b;
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
        hs[i] = 1.0/pow(2.0, (double)order) *
        bspline::factorial(order+1)/(double)(bspline::factorial(i)*bspline::factorial(order+1-i));
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


