//
// $Id$
//
//
// Author: Ranjeet Bhamber <ranjeet <a.t> bristol.ac.uk>
//
// Copyright (C) 2015  Biospi Laboratory for Medical Bioinformatics, University of Bristol, UK
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

#include "imagecore.hpp"


MassSpecData::MassSpecData(void)
{
	mzRange.first=0.0;
	mzRange.second=0.0;
	N=0;
}

MassSpecData::MassSpecData(int mz, int rt)
{
	N=0;
}

void MassSpecData::calMZi(void)
{
	mzi.resize(sci.size(),0);

	for(lli i=0; i < mzi.size(); ++i)
	{
		mzi[i]=sci[i]+i;
	}
}

void MassSpecData::calRange(void)
{
	mzRange.first=mz[mzi[0]];
	mzRange.second=mz[mzi[1]-1];
	vector<pair<double,double> > delta;
	delta.push_back(mzRange);
	for(size_t i=1; i < mzi.size()-1; ++i)
	{
		double minMZ=mz[mzi[i]];
		double maxMZ=mz[mzi[i+1]-1];

		// We want largest MZ range
		if (minMZ < mzRange.first) mzRange.first=minMZ;
		if (maxMZ > mzRange.second) mzRange.second=maxMZ;
	}
	cout<<"Range Min MZ: "<<mzRange.first<<"   Range Max MZ: "<<mzRange.second<<endl;
	rtRange.first=rt[0];
	lli idx=0;
	do{
		rtRange.first=rt[idx];
		idx++;
	} while(rt[idx-1]<0);
	rtRange.second=rt.back();
	cout<<"Range Min RT(min): "<<rtRange.first/60.0<<"   Range Max RT(min): "<<rtRange.second/60.0<<endl;
	cout<<"Range Min RT(s):   "<<rtRange.first<<"   Range Max RT(s):   "<<rtRange.second<<endl;
}


MassSpecImage::MassSpecImage(pair<ui,ui> _xyview) : xypxl(_xyview)
{
	mz.resize(xypxl.first+1,0);// Number of pixels along MZ axis, MZ are edges of bins
	rt.resize(xypxl.second+1,0); // Number of pixels along RT axis
	sc.resize(xypxl.first*xypxl.second,0.0);
	drt=0.0;
	dmz=0.0;
}

void scanMZ(MassSpecData &raw, MassSpecImage &imgbox, lli idx, lli row, const double scaleRT)
{
	for(lli i=raw.mzi[idx]; i < raw.mzi[idx+1]; ++i)
	{
		double dmz=fabs(raw.mz[i+1]-raw.mz[i]);
		double xn=(raw.mz[i]-imgbox.mz[0])/imgbox.dmz;
		double xp=(raw.mz[i+1]-imgbox.mz[0])/imgbox.dmz;
		lli n = floor(xn);
		lli p = floor(xp);
		lli rowOffSet= (row > 0) ? row*(imgbox.xypxl.first): 0;

		if(p > imgbox.xypxl.first-1) p=imgbox.xypxl.first-1;
		if(n >= imgbox.xypxl.first-1) break;
		// Overlapping image 1st box on boundary
		if(n < 0 && p >= 0)
		{
			n=0;
			for(lli col=n; col <= p; ++col)
			{
				double boxBegin=0.0;
				double boxEnd=0.0;
				double scaleMZ=0.0;
				// Cut top overlap else its contained
				(imgbox.mz[col] < raw.mz[i]) ? boxBegin=raw.mz[i] : boxBegin=imgbox.mz[col];
				// Cut bottom overlap else its contained
				(imgbox.mz[col+1] < raw.mz[i+1]) ? boxEnd=imgbox.mz[col+1] : boxEnd=raw.mz[i+1];
				// Calculate MZ scaling
				scaleMZ=(boxEnd-boxBegin)/dmz;
				// Update Spectrum Count
				imgbox.sc[col+rowOffSet]=imgbox.sc[col+rowOffSet]+raw.sc[i-idx]*float(scaleRT*scaleMZ);
			}
		}
		// Multiple smaller image boxes within single raw data scanline
		else if(n < p && n >=0)
		{
			for(lli col=n; col <= p; ++col)
			{
				double boxBegin=0.0;
				double boxEnd=0.0;
				double scaleMZ=0.0;
				// Cut top overlap else its contained
				(imgbox.mz[col] < raw.mz[i]) ? boxBegin=raw.mz[i] : boxBegin=imgbox.mz[col];
				// Cut bottom overlap else its contained
				(imgbox.mz[col+1] < raw.mz[i+1]) ? boxEnd=imgbox.mz[col+1] : boxEnd=raw.mz[i+1];
				// Calculate MZ scaling
				scaleMZ=(boxEnd-boxBegin)/dmz;
				// Update Spectrum Count
				imgbox.sc[col+rowOffSet]=imgbox.sc[col+rowOffSet]+raw.sc[i-idx]*float(scaleRT*scaleMZ);
			}
		}
		// Raw data scanline contained within a Image box hence no scale needed.
		else if(n == p && n > 0)
		{
			lli col =n;
			double scaleMZ=1.0;
			// Update Sectural Count
			imgbox.sc[col+rowOffSet]=imgbox.sc[col+rowOffSet]+raw.sc[i-idx]*float(scaleRT*scaleMZ);
		}
	}
}
