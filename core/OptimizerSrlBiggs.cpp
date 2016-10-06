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


#include "OptimizerSrlBiggs.hpp"
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <omp.h>


using namespace std;


OptimizerSrlBiggs::OptimizerSrlBiggs(const vector<Basis*>& bases, const Matrix& g, fp pruneThreshold) : OptimizerSrl(bases, g, pruneThreshold)
{
#ifndef NDEBUG
	cout << "Init Optimizer Biggs-Andrews accelerated SRL..." << endl;
#endif

	// temporaries required for acceleration
	c0s_.resize(bases.size());
	u0s_.resize(bases.size());
	for (ii i = 0; i < (ii)bases.size(); i++)
	{
		if (!bases[i]->isTransient())
		{
			c0s_[i].init(bases[i]->getM(), bases[i]->getN());
			u0s_[i].init(bases[i]->getM(), bases[i]->getN());
		}
	}

	// how much memory are we using?
	li mem = 0;
	for (ii i = 0; i < c0s_.size(); i++) mem += c0s_[i].mem();
	for (ii i = 0; i < u0s_.size(); i++) mem += u0s_[i].mem();
	cout << "Biggs-Andrews Acceleration mem=" << defaultfloat << setprecision(3) << mem / (1024.0*1024.0) << "Mb" << endl;
}


OptimizerSrlBiggs::~OptimizerSrlBiggs()
{
}


void OptimizerSrlBiggs::update(std::vector<Matrix>& cs, std::vector<Matrix>& c1s, ii iteration)
{
	if (iteration == 1)
	{
		for (ii i = 0; i < (ii)cs.size(); i++)
		{
			if (!!cs[i])
			{
				// save for next iteration
				u0s_[i].elementwiseDiv(c1s[i], cs[i]); 
				c0s_[i].copy(c1s[i]); 

				cs[i].copy(c1s[i]); 
				c1s[i].free();
			}
		}
	}
	else
	{
		double sumU0U = 0.0;
		double sumU0U0 = 0.0;
		for (ii i = 0; i < (ii)cs.size(); i++)
		{
			if (!!cs[i])
			{
				Matrix oldU0;
				oldU0.elementwiseLn(u0s_[i]);
				oldU0.elementwiseMul(c0s_[i], oldU0);
				sumU0U0 += oldU0.sumSqrs();

				// save for next iteration
				u0s_[i].elementwiseDiv(c1s[i], cs[i]);

				Matrix newU0;
				newU0.elementwiseLn(u0s_[i]);
				newU0.elementwiseMul(c1s[i], newU0);
				newU0.elementwiseMul(oldU0, newU0);
				sumU0U += newU0.sum();
			}
		}

		// linear vector extrapolation
		fp a = (fp)sqrt(sumU0U / sumU0U0);
		a = a > 0.0f ? a : 0.0f;
		a = a < 1.0f ? a : 1.0f;

		cout << "a=" << a << endl;

		for (ii i = 0; i < (ii)cs.size(); i++)
		{
			if (!!cs[i])
			{
				Matrix c;
				c.elementwiseDiv(c1s[i], c0s_[i]);
				c.elementwisePow(c, a);
				c.elementwiseMul(c, c1s[i]);

				cs[i].copy(c); // for this iteration
				c0s_[i].copy(c1s[i]); // save for next iteration
				c1s[i].free();
			}
		}
	}
}

// quadratic vector extrapolation below (todo)
/*
		else if (iteration_ == 1) // linear vector extrapolation this time, but save the qs
		{
			for (ii j = 0; j < (ii)bases_.size(); j++)
			{
				if (!bases_[j]->isTransient())
				{
					#pragma omp parallel for reduction(+:sum,sumd)
					for (li i = 0; i < cs_[j].size(); i++)
					{
						if (cs_[j].vs_[i] > 0.0)
						{
							fp q = cEs[j].vs_[i] / c0s_[j].vs_[i];
							fp c1 = cEs[j].vs_[i] * powf(cEs[j].vs_[i] / c0s_[j].vs_[i], a);

							sum += cs_[j].vs_[i] * cs_[j].vs_[i];
							sumd += (c1 - cs_[j].vs_[i])*(c1 - cs_[j].vs_[i]);

							// for this itteration
							cs_[j].vs_[i] = c1;

							// for next itteration
							c0s_[j].vs_[i] = cEs[j].vs_[i];
							q0s_[j].vs_[i] = q;
						}
					}
					cEs[j].free();
				}
			}
		}
		else // quadratic vector extrapolation
		{
			for (ii j = 0; j < (ii)bases_.size(); j++)
			{
				if (!bases_[j]->isTransient())
				{
					#pragma omp parallel for reduction(+:sum,sumd)
					for (li i = 0; i < cs_[j].size(); i++)
					{
						if (cs_[j].vs_[i] > 0.0)
						{
							fp q = cEs[j].vs_[i] / c0s_[j].vs_[i];
							fp c1 = cEs[j].vs_[i] * powf(q, a) * powf(q / q0s_[j].vs_[i], 0.5f*a*a);

							sum += cs_[j].vs_[i] * cs_[j].vs_[i];
							sumd += (c1 - cs_[j].vs_[i])*(c1 - cs_[j].vs_[i]);

							// for this itteration
							cs_[j].vs_[i] = c1;

							// for next itteration
							c0s_[j].vs_[i] = cEs[j].vs_[i];
							q0s_[j].vs_[i] = q;
						}
					}
					cEs[j].free();
				}
			}
		}
	}

	return sqrt(sumd) / sqrt(sum);
*/



