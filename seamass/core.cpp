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


#include "core.hpp"
#include <limits>

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


ii
CoeffsMetadata::
size() const
{
    ii size = 1;
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
        if (i < d-1) out << ",";
    }
    out << "] oc=[";
    for (ii i = 0; i < d; i++)
    {
        if (o[i] == numeric_limits<ii>::min()) out << "?"; else out << o[i];  
        if (i < d-1) out << ",";
    }
    out << "] nc=[";
    for (ii i = 0; i < d; i++)
    {
        if (n[i] == numeric_limits<ii>::min()) out << "?"; else out << n[i];  
        if (i < d-1) out << ",";
    }
    out << "]:" << size();
}


////////////////////////////////////////////////////////////////////////////////


void bin_mzs_intensities(vector< vector<double> >& mzs,
                         vector< vector<double> >& intensities,
                         ii instrument_type)
{
    // This modifies the raw data for some limitations of the mzML spec and makes
    // sure the intensities are treated as binned between m/z datapoints.
    //
    // For ToF data, it interpolates the mzs to represent the edges of the bins rather
    //   than the centres
    // For FT data, it converts the sampled data to binned counts by treating the
    //   mz values as the bin edges, and using trapezoid rule to integrate intensity values
    //
    // This is all a bit rough at the moment, should be fitting splines to the data

    // if more than one spectrum, ignore last as we do not know its scan end time
    if (mzs.size() > 1) mzs.resize(mzs.size() - 1);
    if (intensities.size() > 1) intensities.resize(intensities.size() - 1);
    
    // Looks like calculation for centrioded Orbitrap as it is creating MZ on the edges of the bins.
    if (instrument_type == 2)
    {
        #pragma omp parallel for
        for (ii j = 0; j < mzs.size(); j++)
        if (mzs[j].size() >= 2)
        {
            // we drop the first and last m/z datapoint as we don't know both their bin edges
            for (ii i = 1; i < mzs[j].size(); i++)
            {
                // linear interpolation of mz extent (probably should be cubic)
                mzs[j][i-1] = 0.5 * (mzs[j][i-1] + mzs[j][i]);
                intensities[j][i-1] = intensities[j][i];
            }
            mzs[j].resize(mzs[j].size() - 1);
            intensities[j].resize(intensities[j].size() - 2);
        }
        else
        {
            mzs[j].resize(0);
            intensities[j].resize(0);
        }
    }
    /*else if (instrument_type == 1) // Orbitrap - but disabled for now as doesn't play nicely with merge_bins function
    {
        #pragma omp parallel for
        for (ii j = 0; j < mzs.size(); j++)
        if (mzs[j].size() >= 2)
        {
            // for Orbitrap only, mark zeros as missing data
            for (ii i = 1; i < mzs[j].size(); i++)
            if (intensities[j][i-1] <= 0.0 || intensities[j][i] <= 0.0)
            {
                intensities[j][i-1] = -1.0;
            }
            else
            {
                intensities[j][i-1] = (mzs[j][i] - mzs[j][i-1]) * 0.5 * (intensities[j][i] + intensities[j][i-1]);
            }
            intensities[j].resize(intensities[j].size() - 1);
        }
        else
        {
            mzs[j].resize(0);
            intensities[j].resize(0);
        }
    }*/
    else // This looks like the Calculation for Time of Flight (ToF) as it is integrating the Intensities
    {
        #pragma omp parallel for
        for (ii j = 0; j < mzs.size(); j++)
        if (mzs[j].size() >= 2)
        {
            // dividing by minimum to get back to ion counts for SWATH data which appears to be scaled to correct for dynamic range restrictions (hack!)
            double minimum = std::numeric_limits<double>::max();
            for (ii i = 1; i < mzs[j].size(); i++)
            {
                if (intensities[j][i-1] > 0) minimum = minimum < intensities[j][i-1] ? minimum : intensities[j][i-1];
                intensities[j][i-1] = (mzs[j][i] - mzs[j][i-1]) * 0.5 * (intensities[j][i] + intensities[j][i-1]);
            }
            intensities[j].resize(intensities[j].size() - 1);
            for (ii i = 0; i < intensities[j].size(); i++)
            {
                intensities[j][i] /= minimum;
            }
            cout << minimum << endl;
        }
        else
        {
            mzs[j].resize(0);
            intensities[j].resize(0);
        }
    }
}


////////////////////////////////////////////////////////////////////////////////


void merge_bins(vector< vector<double> >& mzs,
                vector< vector<double> >& intensities,
                double width)
{
    #pragma omp parallel for
    for (ii j = 0; j < mzs.size(); j++)
    {
        double w = 0;
        double v = 0;
        ii n = 0;
        ii k = 1;
        for (ii i = 1; i < mzs[j].size();)
        {
            w += mzs[j][i] - mzs[j][i-1];
            n++;
            
            if (w > width || i == mzs[j].size()-1)
            {
                if (n == 1)
                {
                    mzs[j][k] = mzs[j][i];
                    intensities[j][k-1] = intensities[j][i-1];
                    i++;
                }
                else
                {
                    mzs[j][k] = mzs[j][i-1];
                    intensities[j][k-1] = v;
                }
                
                w = 0;
                v = 0;
                n = 0;
                k++;
            }
            else
            {
                i++;
            }
            
            if (i < intensities[j].size()) v += intensities[j][i-1];
        }
        
        mzs[j].resize(k);
        intensities[j].resize(k-1);
    }
}


void create_gs(vector<fp>& gs,
               vector<li>& is,
               vector<ii>& js,
               const vector< vector<double> >& intensities)
{
    ii nj = 0;
    for (ii j = 0; j < intensities.size(); j++)
    if (intensities[j].size() > 0) nj++;
    
    js.resize(nj);
    is.resize(nj+1);
    is.front() = 0;
    for (ii i = 0, j = 0; j < intensities.size(); j++)
    if (intensities[j].size() > 0)
    {
        js[i] = j;
        is[i+1] = is[i] + intensities[j].size();
        i++;
    }
    
    gs.resize(is.back());
    #pragma omp parallel for
    for (ii j = 0; j < nj; j++)
    for (ii i = 0; i < intensities[js[j]].size(); i++)
    {
         gs[is[j]+i] = intensities[js[j]][i];
    }
    
    cout << "Raw gs=[" << js.size() << "/" << intensities.size() << "]:" << gs.size() << " mem=" << fixed << setprecision(2) << (is.size()*sizeof(li) + js.size()*sizeof(ii) + gs.size()*sizeof(fp))/(1024.0*1024.0) << "Mb";
    cout << endl;
}


////////////////////////////////////////////////////////////////////////////////


namespace bspline {

double
m(double x, int k, int i, vector<fp>& ks)
{
    if (k == 1)
    {
        if (ks[i] <= x && x < ks[i+1])
        {
            return 1.0 / (ks[i+1] - ks[i]);
        }
        else
        {
            return 0.0;
        }
    }
    else
    {
        if (ks[i+k]-ks[i] == 0.0)
        {
            return 0.0;
        }
        else
        {
            return (k*((x-ks[i]) * m(x,k-1,i,ks) +
                       (ks[i+k]-x) * m(x,k-1,i+1,ks)))/((k-1)*(ks[i+k]-ks[i]));
        }
    }
}


double
m(double x, int k, int i)
{
    if (k == 1)
    {
        if (i <= x && x < i+1)
        {
            return 1.0;
        }
        else
        {
            return 0.0;
        }
    }
    else
    {
        return (k*((x-i) * m(x,k-1,i) +
                   ((i+k)-x) * m(x,k-1,i+1)))/((k-1)*k);
    }
}


double
im(double x, int k)
{
    double v = 0.0;
    for (int i = 0; i < k+1; ++i)
    {
        v += m(x,k+1,i);
    }
    return v;
}
 
    
int
factorial(int n)
{
    return (n == 1 || n == 0) ? 1 : factorial(n - 1) * n;
}

}



