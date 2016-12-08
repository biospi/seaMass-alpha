//
// Original author: Andrew Dowsey <andrew.dowsey <a.t> liverpool.ac.uk>
//
// Copyright (C) 2015  biospi Laboratory, EEE, University of Liverpool, UK
//
// This file is part of seaMass-TD.
//
// seaMass-TD is free software: you can redistribute it and/or modify
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
// along with seaMass-TD.  If not, see <http://www.gnu.org/licenses/>.
//


#include "BasisIsotopeDistribution.hpp"
#include <limits>
#include <map>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <cmath>
#include <cfloat>
using namespace std;


/*struct Factor
{
	ii gi;
	fp ave_mass;
	vector<fp> p;

	Factor(ii _gi) : gi(_gi), ave_mass(0.0) {}
};


BasisIsotopeDistribution::
BasisIsotopeDistribution(vector<Basis*>& bases, BasisChargeDistribution* _parent,
	                     ii out_res, ii factor_res, ii max_z, bool transient) :
	Basis(bases, _parent, transient),
	parent(_parent),
	ns(parent->get_ns()),
	nc(parent->get_nc())
{
	///////////////////////////////////////////////////////////////////////
	// load learnt factors into multimap

	ifstream ifs("seamass-td.factors");
	multimap<ii, Factor>* factors = new multimap<ii, Factor>;
	ii fi0 = parent->get_ci0s().front() / (1 << (out_res - factor_res));
	ii fi1 = parent->get_ci1s().back() / (1 << (out_res - factor_res));

	cout << get_index() << " BasisIsotopeDistribution " << flush;
	while (ifs.good())
	{
		char comma;
		ii fi, z, fgi;
		ifs >> fi;
		if (fi < fi0) { ifs.ignore(65535, '\n'); continue; }
		if (fi > fi1) break;

		ifs >> comma >> z >> comma >> fgi;
		map<ii,Factor>::iterator elem = factors->insert(pair<ii,Factor>(fi, Factor(fgi)));
		for (ii j = 0; j < 101; j++)
		{
			fp isotope;
			ifs >> comma >> isotope;
			if (elem->second.p.size() > 0 && isotope <= elem->second.p.back() && isotope < 0.00001)
			{
				ifs.ignore(65535, '\n');
				break;
			}
			else
			{
				elem->second.ave_mass += j * isotope;
				elem->second.p.push_back(isotope);
			}
		}
		//ofs << ci << " " << z << " " << id << " " << elem->second.size() << endl;

		if ((fi-fi0) % 1000 == 0)
		{
			for (int i = 0; i < 256; ++i) cout << '\b';
			cout << get_index() << " BasisIsotopeDistribution " << setw(1 + (int)(log10((float)fi1))) << fi-fi0 << "/" << fi1-fi0 << " " << flush;
		}
	}
	for (int i = 0; i < 256; ++i) cout << '\b';
	ifs.close();

	///////////////////////////////////////////////////////////////////////
	// create A as a temporary COO matrix

	ii ni = 0;
	ii ng = 0;
	ii nnz = 0;
	for (ii oi = parent->get_ci0s().front(); oi <= parent->get_ci1s().back(); oi++)
	{
		ii fi = oi / (1 << (out_res - factor_res));
		ii last_fgi = -1;
		for (pair<multimap<ii, Factor>::iterator, multimap<ii, Factor>::iterator> fs = factors->equal_range(fi); fs.first != fs.second; ++fs.first)
		{
			ii fgi = fs.second->first;
			if (fgi != last_fgi)
			{
				ng++;
				last_fgi = fgi;
			}

			Factor& fac = fs.second->second;
			for (ii z = 0; z < parent->get_ci0s().size(); z++)
			if (oi >= parent->get_ci0s()[z] && oi <= parent->get_ci1s()[z])
			{
				for (ii i = 0; i < fs.first->second.p.size(); i++)
				{
					ii pi = oi + i * (1 << out_res);
					if (pi <= parent->get_ci1s()[z])
					{
						nnz++;
					}
				}
				ni++;
			}
		}
	}

	vector<fp> acoo(nnz);
	vector<ii> rowind(nnz);
	vector<ii> colind(nnz);
	gis.resize(ng+1);
	gis[0] = 0;
	ois.resize(parent->get_ci1s().back() - parent->get_ci0s().front() + 1);
	ois[0] = 0;
	ave_masses.resize(ni);

	ii ci = 0;
	ii gi = 0;
	ii k = 0;
	cout << get_index() << " BasisIsotopeDistribution " << flush;
	for (ii oi = parent->get_ci0s().front(); oi <= parent->get_ci1s().back(); oi++)
	{
		ii fi = oi / (1 << (out_res - factor_res));
		ii last_fgi = -1;
		for (pair<multimap<ii, Factor>::iterator, multimap<ii, Factor>::iterator> fs = factors->equal_range(fi); fs.first != fs.second; ++fs.first)
		{
			ii fgi = fs.second->first;
			if (fgi != last_fgi)
			{
				gi++;
				gis[gi] = ci;
				last_fgi = fgi;
			}

			Factor& fac = fs.second->second;
			for (ii z = 0; z < parent->get_ci0s().size(); z++)
			if (oi >= parent->get_ci0s()[z] && oi <= parent->get_ci1s()[z])
			{
				for (ii i = 0; i < fs.first->second.p.size(); i++)
				{
					ii pi = oi + i * (1 << out_res);
					if (pi <= parent->get_ci1s()[z])
					{
						acoo[k] = fs.first->second.p[i];
						rowind[k] = pi - parent->get_ci0s()[z] + parent->get_cos()[z];
						colind[k] = ci;
						k++;
					}
				}
				//ave_masses[ci] = fs.first->second.ave_mass;
				ci++;
			}
		}
		ois[oi - parent->get_ci0s().front() + 1] = ci;

		if (oi % 1000 == 0)
		{
			for (int i = 0; i < 256; ++i) cout << '\b';
			cout << get_index() << " BasisIsotopeDistribution " << setw(1 + (int)(log10((float)parent->get_ci1s().back() - parent->get_ci0s().front() + 1))) << oi - parent->get_ci0s().front() << "/" << parent->get_ci1s().back() - parent->get_ci0s().front() + 1 << " " << flush;
		}
	}
	for (int i = 0; i < 256; ++i) cout << '\b';
	nc = ns * ci;

	// create A
	delete factors;
	a.init(parent->get_cos().back(), nc / ns, acoo, rowind, colind);

	cout << get_index() << " BasisIsotopeDistribution ";
	cout << " A=";
	a.print(cout);
	cout << " ng=" << ng;
	if (transient) cout << " (t)" << endl; else cout << " nc=" << nc << endl;
}


BasisIsotopeDistribution::~BasisIsotopeDistribution()
{
}


void
BasisIsotopeDistribution::
synthesis(vector<fp>& fs, const vector<fp>& cs, bool accum) const
{
	for (ii j = 0; j < ns; j++)
	{
		a.mult(&(fs.data()[j*a.get_m()]), &(cs.data()[j*a.get_n()]), false, accum);
	}
}


void
BasisIsotopeDistribution::
analysis(vector<fp>& es, const vector<fp>& fs) const
{
	for (ii j = 0; j < ns; j++)
	{
		a.mult(&(es.data()[j*a.get_n()]), &(fs.data()[j*a.get_m()]), true);
	}
}


void
BasisIsotopeDistribution::
l2norm(vector<fp>& es, const vector<fp>& fs) const
{
	for (ii j = 0; j < ns; j++)
	{
		a.sqr_mult(&(es.data()[j*a.get_n()]), &(fs.data()[j*a.get_m()]), true);
	}
}


void
BasisIsotopeDistribution::
shrink(std::vector<fp>& es, const std::vector<fp>& cs, const std::vector<fp>& l2, const std::vector<fp>& wcs, double shrinkage) const
{
	// GROUP-WISE SHRINKAGE! (note - intuitive implemention, not mathematically verified yet)
	ii n = nc / ns;
	#pragma omp parallel for
	for (ii g = 0; g < gis.size() - 1; g++)
	{
		// sum up coefficients per group
		fp sum = 0.0;
		for (ii j = 0; j < ns; j++)
		for (ii i = gis[g]; i < gis[g + 1]; i++)
		{
			sum += cs[j*n + i];
		}

		// scale the shrinkage to be proportional to the contribution of this coefficient to the group total
		for (ii j = 0; j < ns; j++)
		for (ii i = gis[g]; i < gis[g + 1]; i++)
		{
			double scale = cs[j*n + i] / sum;
			es[j*n + i] *= cs[j*n + i] / (scale * shrinkage * l2[j*n + i] + wcs[j*n + i]);
		}
	}
}


void
BasisIsotopeDistribution::
write_cs(const std::vector<fp>& cs) const
{
	ii n = nc / ns;
	vector<fp> sums(ois.size() - 1, 0.0);
	for (ii o = 0; o < ois.size() - 1; o++)
	{
		for (ii j = 0; j < ns; j++)
		for (ii i = ois[o]; i < ois[o + 1]; i++)
		{
			sums[o] += cs[j*n + i];
		}
	}

	ostringstream oss; oss << "profile" << get_index() << ".csv";
	ofstream ofs(oss.str().c_str());
	ofs << "mass,intensity" << setprecision(10) << endl;
	for (ii o = 0; o < ois.size() - 1; o++)
	{
		if (sums[o] > 0.0 || o > 0 && sums[o - 1] > 0.0 || o < sums.size() - 1 && sums[o + 1] > 0.0)
		{
			ofs << (parent->get_ci0s().front() + o) * parent->get_mass_interval() << "," << sums[o] << endl;
		}
	}

	ostringstream oss2; oss2 << "all" << get_index() << ".csv";
	ofstream ofs2(oss2.str().c_str());
	ofs2 << "mass,intensity";
	for (ii f = 0; f < 3; f++)
	for (ii z = 0; z < parent->get_ci0s().size(); z++)
	{
		ofs2 << "," << "z" << z+1 << "." << f;
	} 
	ofs2 << setprecision(10) << endl;
	for (ii o = 0; o < ois.size() - 1; o++)
	{
		if (sums[o] > 0.0 || o > 0 && sums[o - 1] > 0.0 || o < sums.size() - 1 && sums[o + 1] > 0.0)
		{
			ofs2 << (parent->get_ci0s().front() + o) * parent->get_mass_interval() << "," << sums[o];

			for (ii i = ois[o]; i < ois[o + 1];)
			{
				for (ii z = 0; z < parent->get_ci0s().size(); z++)
				{
					ofs2 << ",";
					if (parent->get_ci0s().front() + o >= parent->get_ci0s()[z] && parent->get_ci0s().front() + o <= parent->get_ci1s()[z])
					{
						fp sum = 0.0;
						for (ii j = 0; j < ns; j++) sum += cs[j*n + i];
						ofs2 << sum;
						i++;
					}						
				}
			}

			ofs2 << endl;
		}
	}
}

void
BasisIsotopeDistribution::
restrict_range(std::vector<fp>& cs, double mass0, double mass1) const
{
	ii n = nc / ns;
	vector<fp> sums(ois.size() - 1, 0.0);
	for (ii o = 0; o < ois.size() - 1; o++)
	{
		if ((parent->get_ci0s().front() + o) * parent->get_mass_interval() < mass0 || (parent->get_ci0s().front() + o) * parent->get_mass_interval()  > mass1)
		for (ii j = 0; j < ns; j++)
		for (ii i = ois[o]; i < ois[o + 1]; i++)
		{
			cs[j*n + i] = 0.0;
		}
	}
}*/