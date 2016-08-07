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




double CatmullRomInterpolate(
	double y0, double y1,
	double y2, double y3,
	double mu)
{
	double a0, a1, a2, a3, mu2;

	mu2 = mu*mu;
	a0 = -0.5*y0 + 1.5*y1 - 1.5*y2 + 0.5*y3;
	a1 = y0 - 2.5*y1 + 2 * y2 - 0.5*y3;
	a2 = -0.5*y0 + 0.5*y2;
	a3 = y1;

	return a0*mu*mu2 + a1*mu2 + a2*mu + a3;
}


void remove_zeros(vector< vector<double> >& mzs, vector< vector<double> >& intensities)
{
	for (ii j = 0; j < intensities.size(); j++)
	{
		ii shift = 0;
		for (ii i = 0; i < intensities[j].size(); i++)
		{
			if (intensities[j][i] <= 0.0)
			{
				shift++;
			}
			else if (shift > 0 )
			{
				intensities[j][i - shift] = intensities[j][i];
				mzs[j][i - shift] = mzs[j][i];
			}
		}
		intensities[j].resize(intensities[j].size() - shift);
		mzs[j].resize(mzs[j].size() - shift);
	}
}


////////////////////////////////////////////////////////////////////////////////


void merge_bins(vector< vector<fp> >& mzs,
	vector< vector<fp> >& intensities,
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
        
        if(k != 1) mzs[j].resize(k);
        intensities[j].resize(k-1);
    }
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



