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


#include "OptimizerAsrl.hpp"
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <omp.h>


using namespace std;


OptimizerAsrl::OptimizerAsrl(const vector<Basis*>& bases, const Matrix& g, ii accelleration) : bases_(bases), g_(g), accelleration_(accelleration), iteration_(0)
{
#ifndef NDEBUG
	cout << "Init Optimizer ASRL..." << endl;
#endif

	// compute L1 norm of each basis function and store in 'l1s'
#ifndef NDEBUG
	cout << " Init L1s..." << endl;
#endif
	l1s_.resize(bases_.size());
	for (ii i = 0; i < (ii)bases_.size(); i++)
	{
		if (i == 0)
		{
			Matrix t; t.init(g_.m(), g_.n(), 1.0);
			bases_[i]->analysis(l1s_[0], t, false);
		}
		else
		{
			bases_[i]->analysis(l1s_[i], l1s_[bases_[i]->getParentIndex()], false);
		}
	}
	for (ii i = 0; i < (ii)bases_.size(); i++)
	{
		if (bases_[i]->isTransient())
		{
			l1s_[i].free();
		}
	}

	// compute L2 norm of each basis function and store in 'l2s'
#ifndef NDEBUG
	cout << " Init L2s..." << endl;
#endif
	l2s_.resize(bases_.size());
	for (ii i = 0; i < (ii)bases_.size(); i++)
	{
		if (i == 0)
		{
			Matrix t; t.init(g_.m(), g_.n(), 1.0);
			bases_[i]->analysis(l2s_[0], t, true);
		}
		else
		{
			bases_[i]->analysis(l2s_[i], l2s_[bases_[i]->getParentIndex()], true);
		}
	}
	for (ii i = 0; i < (ii)bases_.size(); i++)
	{
		if (bases_[i]->isTransient())
		{
			l2s_[i].free();
		}
		else
		{
			l2s_[i].elementwiseSqrt(l2s_[i]);
		}
	}

	// initialise starting estimate of 'c' from analysis of 'g'
#ifndef NDEBUG
	cout << " Init Cs from G" << endl;
#endif
	cs_.resize(bases_.size());
	for (ii i = 0; i < (ii)bases_.size(); i++)
	{
		if (i == 0)
		{
			bases_[i]->analysis(cs_[0], g_, false);
		}
		else
		{
			bases_[i]->analysis(cs_[i], cs_[bases_[i]->getParentIndex()], false);
		}
	}
	ii notTransient = 0;
	for (ii i = 0; i < (ii)bases_.size(); i++)
	{
		if (bases_[i]->isTransient())
		{
			cs_[i].free();
		}
		else
		{
			notTransient++;
		}
	}
	for (ii i = 0; i < (ii)bases_.size(); i++)
	{
		if (!bases_[i]->isTransient())
		{
			cs_[i].elementwiseDiv(cs_[i], l1s_[i], (fp)0.0);
			cs_[i].elementwiseMul((fp)1.0 / notTransient, cs_[i]);
		}
	}

	// temporaries required for acceleration
	if (accelleration_ >= 1)
	{
		c0s_.resize(bases_.size());
		u0s_.resize(bases_.size());
		for (ii i = 0; i < (ii)bases.size(); i++)
		{
			if (!bases[i]->isTransient())
			{
				c0s_[i].init(bases_[i]->getM(), bases_[i]->getN());
				u0s_[i].init(bases_[i]->getM(), bases_[i]->getN());
			}
		}
	}
	if (accelleration_ >= 2)
	{
		q0s_.resize(bases.size());
		for (ii i = 0; i < (ii)bases.size(); i++)
		{
			if (!bases[i]->isTransient())
			{
				q0s_[i].init(bases_[i]->getM(), bases_[i]->getN());
			}
		}
	}

	// how much memory are we using?
	li mem = 0;
	for (ii i = 0; i < cs_.size(); i++) mem += cs_[i].mem();
	for (ii i = 0; i < l1s_.size(); i++) mem += l1s_[i].mem();
	for (ii i = 0; i < l2s_.size(); i++) mem += l2s_[i].mem();
	for (ii i = 0; i < c0s_.size(); i++) mem += c0s_[i].mem();
	for (ii i = 0; i < u0s_.size(); i++) mem += u0s_[i].mem();
	for (ii i = 0; i < q0s_.size(); i++) mem += q0s_[i].mem();
	cout << "Vector Extrapolated Sparse Richardson-Lucy mem=" << defaultfloat << setprecision(3) << mem / (1024.0*1024.0) << "Mb" << endl;
}


OptimizerAsrl::~OptimizerAsrl()
{
}


double
OptimizerAsrl::
step(fp lambda)
{
	iteration_++;

	// temporary variables
	Matrix f;
	vector<Matrix> cEs(bases_.size());
	Basis::ErrorInfo errorInfo;
	double sumSqrs = 0.0;
	double sumSqrDiffs = 0.0;

	// SYNTHESIS
#ifndef NDEBUG
	cout << iteration_ << " synthesis" << endl;
#endif
	double synthesisStart = omp_get_wtime();
	{
		synthesis(f);
	}
	double synthesisDuration = omp_get_wtime() - synthesisStart;

    // ERROR
#ifndef NDEBUG
	cout << iteration_ << " error" << endl;
#endif
	double errorStart = omp_get_wtime();
	{
		errorInfo = bases_[0]->error(f, f, g_);
	}
	double errorDuration = omp_get_wtime() - errorStart;

    // ANALYSIS
#ifndef NDEBUG
	cout << iteration_ << " analysis" << endl;
#endif
	double analysisStart = omp_get_wtime();
	{
		for (ii i = 0; i < (ii)bases_.size(); i++)
		{
			if (i == 0)
			{
				bases_.front()->analysis(cEs[0], f, false);
				f.free();
			}
			else
			{
				bases_[i]->analysis(cEs[i], cEs[bases_[i]->getParentIndex()], false);
			}
		}
	}
	double analysisDuration = omp_get_wtime() - analysisStart;

	// SHRINKAGE
#ifndef NDEBUG
	cout << iteration_ << " shrinkage" << endl;
#endif
	double shrinkageStart = omp_get_wtime();
	{
		for (ii i = 0; i < (ii)bases_.size(); i++)
		{
			if (!bases_[i]->isTransient())
			{
				bases_[i]->shrinkage(cEs[i], cEs[i], cs_[i], l1s_[i], l2s_[i], lambda);
			}
			else
			{
				cEs[i].free();
			}
		}
	}
    double shrinkageDuration = omp_get_wtime() - shrinkageStart;
   
    // ACCELERATION (currently turned off)
#ifndef NDEBUG
	cout << iteration_ << " acceleration" << endl;
#endif
	double accelerationStart = omp_get_wtime();
	{
		//if (accelleration_ == 0)
		{
			for (ii i = 0; i < (ii)bases_.size(); i++)
			{
				if (!bases_[i]->isTransient())
				{
					sumSqrs += cs_[i].sumSqrs();
					sumSqrDiffs += cEs[i].sumSqrDiffs(cs_[i]);

					cs_[i].elementwiseMul((fp)1.0, cEs[i]);
				}
			}
		}
	}
	double accelerationDuration = omp_get_wtime() - accelerationStart;

#ifndef NDEBUG
	cout << "Durations: synthesis=" << defaultfloat << setprecision(3) << synthesisDuration;
	cout << " error=" << errorDuration;
	cout << " analysis=" << analysisDuration;
	cout << " shrinkage=" << shrinkageDuration;
	cout << " acceleration=" << accelerationDuration;
	cout << " all=" << synthesisDuration + errorDuration + analysisDuration + shrinkageDuration + accelerationDuration << endl;
#endif

	return sqrt(sumSqrDiffs) / sqrt(sumSqrs);

	/*fp a = 0.0;

	else // accelerated update
	{
	// init/update u0s and compute acceleration factor a
	if (iteration_ == 0)
	{
	for (ii j = 0; j < (ii)bases_.size(); j++)
	{
	if (!bases_[j]->isTransient())
	{
	#pragma omp parallel for
	for (li i = 0; i < cs_[j].size(); i++)
	{
	if (cs_[j].vs_[i] > 0.0)
	{
	u0s_[j].vs_[i] = cEs[j].vs_[i] / cs_[j].vs_[i];
	}
	}
	}
	}
	a = 0.0;
	}
	else
	{
	double sum_u0u = 0.0;
	double sum_u0u0 = 0.0;
	for (ii j = 0; j < (ii)bases_.size(); j++)
	{
	if (!bases_[j]->isTransient())
	{
	#pragma omp parallel for reduction(+:sum_u0u,sum_u0u0)
	for (li i = 0; i < cs_[j].size(); i++)
	{
	if (cs_[j].vs_[i] > 0.0)
	{
	double old_u0 = u0s_[j].vs_[i];
	u0s_[j].vs_[i] = cEs[j].vs_[i] / cs_[j].vs_[i];

	old_u0 = old_u0 > 0.0 ? c0s_[j].vs_[i] * log(old_u0) : 0.0;
	sum_u0u += old_u0 * (u0s_[j].vs_[i] > 0.0 ? cEs[j].vs_[i] * log(u0s_[j].vs_[i]) : 0.0);
	sum_u0u0 += old_u0 * old_u0;
	}
	}
	}
	}

	a = (fp) sqrt(sum_u0u / sum_u0u0);
	a = a > 0.0f ? a : 0.0f;
	a = a < 1.0f ? a : 1.0f;
	}

	if (iteration_ == 0) // unaccelerated this time, but save es as c0s
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
	sum += cs_[j].vs_[i] * cs_[j].vs_[i];
	sumd += (cEs[j].vs_[i] - cs_[j].vs_[i])*(cEs[j].vs_[i] - cs_[j].vs_[i]);

	// for this itteration
	cs_[j].vs_[i] = cEs[j].vs_[i];
	// for next itteration
	c0s_[j].vs_[i] = cEs[j].vs_[i];
	}
	}
	cEs[j].free();
	}
	}
	}
	else if (accelleration_ == 1) // linear vector extrapolation
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
	fp c1 = cEs[j].vs_[i] * powf(cEs[j].vs_[i] / c0s_[j].vs_[i], a);

	sum += cs_[j].vs_[i] * cs_[j].vs_[i];
	sumd += (c1 - cs_[j].vs_[i])*(c1 - cs_[j].vs_[i]);

	// for this itteration
	cs_[j].vs_[i] = c1;

	// for next itteration
	c0s_[j].vs_[i] = cEs[j].vs_[i];
	}
	}
	cEs[j].free();
	}
	}
	}
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

	return sqrt(sumd) / sqrt(sum);*/
}


void
OptimizerAsrl::
synthesis(Matrix& f, ii basis) const
{
	vector<Matrix> ts(bases_.size());
	for (ii i = (ii)bases_.size() - 1; i >= 0; i--)
	{
		if (!ts[i] && !bases_[i]->isTransient())
		{
			ts[i].init(bases_[i]->getM(), bases_[i]->getN());
			ts[i].elementwiseDiv(cs_[i], l2s_[i], (fp)0.0);
		}

		if (basis == i) // return with B-spline control points
		{
			f.copy(ts[i]);

			break;
		}

		if (i > 0)
		{
			ii pi = bases_[i]->getParentIndex();
			if (!ts[pi] && !bases_[pi]->isTransient())
			{
				ts[pi].init(bases_[pi]->getM(), bases_[pi]->getN());
				ts[pi].elementwiseDiv(cs_[pi], l2s_[pi], (fp)0.0);
			}

			bases_[i]->synthesis(ts[pi], ts[i], true);
		}
		else
		{
			bases_[0]->synthesis(f, ts[0], false);
		}

		ts[i].free();
	}
}


const std::vector<Matrix>& OptimizerAsrl::getCs() const
{ 
	return cs_;
}


/*void
OptimizerAsrl::
threshold(fp thresh)
{
	for (ii j = 0; j < (ii)bases.size(); j++)
	{
		if (!bases[j]->is_transient())
		{
			#pragma omp parallel for
			for (li i = 0; i < (li)bases[j]->get_cm().size(); i++)
			{
				if (cs[j][i] < thresh) cs[j][i] = 0.0;
			}
		}
	}*/

	/*double volume = bases.front()->get_volume();

	for (ii j = 0; j < (ii)bases.size(); j++)
	{
		if (!bases[j]->is_transient())
		{
			#pragma omp parallel for
			for (li i = 0; i < (li)bases[j]->get_cm().size(); i++)
			{
				cs[j][i] /= volume;
			}
		}
	}
}*/
