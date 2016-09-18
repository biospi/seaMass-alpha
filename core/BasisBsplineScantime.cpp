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


#include "BasisBsplineScantime.hpp"

#include "Bspline.hpp"

#include <limits>
#include <iomanip>


using namespace std;


/*BasisBsplineScantime::BasisBsplineScantime(std::vector<Basis*>& bases, ii parentIndex, const std::vector<double>& start_times, const std::vector<double>& finish_times, const std::vector<fp>& exposures, short resolution, ii order = 3, bool isTransient)
	: BsplineBasis(bases, parent->get_cm().d, parent, transient)
{
	rtMin = start_times.front();
	rtMax = finish_times.back();

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
			double b = Bspline::im(b1, order + 1) - Bspline::im(b0, order + 1);
            if (b <= numeric_limits<float>::min()) b = numeric_limits<float>::min(); // for numerical instability in im()

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
    
    cout << index << " BasisBsplineScantime parent=" << parent->get_index() << " ";
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


BasisBsplineScantime::~BasisBsplineScantime()
{
}


void
BasisBsplineScantime::
synthesis(Matrix& f, const Matrix& c, bool accum)
{
    static fp alpha = 1.0;
    fp beta = accum ? 1.0 : 0.0;
    fp* cs = const_cast<fp*>(c.vs_);
    mkl_scsrmm("N", &m, &(cm.n[0]), &(cm.n[1]), &alpha, "G**C", a.data(), ja.data(), ia.data(), &(ia.data()[1]), cs, &(cm.n[0]),  &beta, f.vs_, &(cm.n[0]));
}


void
BasisBsplineScantime::
analysis(Matrix& e, const Matrix& f, bool a_sqrd)
{
    static fp alpha = 1.0, beta = 0.0;

	if (a_sqrd)
	{
		vector<fp> a2(nnz);
		vsSqr(nnz, a.data(), a2.data());
		fp* fs = const_cast<fp*>(f.vs_);
		mkl_scsrmm("T", &m, &(cm.n[0]), &(cm.n[1]), &alpha, "G**C", a2.data(), ja.data(), ia.data(), &(ia.data()[1]), fs, &(cm.n[0]), &beta, e.vs_, &(cm.n[0]));
	}
	else
	{
		fp* fs = const_cast<fp*>(f.vs_);
		mkl_scsrmm("N", &(cm.n[1]), &(cm.n[0]), &m, &alpha, "G**C", at.data(), jat.data(), iat.data(), &(iat.data()[1]), fs, &(cm.n[0]), &beta, e.vs_, &(cm.n[0]));
	}
}


double BasisBsplineScantime::getMin() const
{
	return rtMin;
}


double BasisBsplineScantime::getMax() const
{
	return rtMax;
}*/


