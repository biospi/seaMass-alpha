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


#include "OptimizerSrl.hpp"
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <omp.h>


using namespace std;


OptimizerSrl::OptimizerSrl(const vector<Basis*>& bases, const Matrix& g, fp pruneThreshold) : bases_(bases), g_(g), pruneThreshold_(pruneThreshold), lambda_(0.0), iteration_(0)
{
#ifndef NDEBUG
	cout << "Init Optimizer SRL..." << endl;
#endif

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

	// compute L1 norm of each L2 normalised basis function and store in 'l1l2s'
#ifndef NDEBUG
	cout << " Init L1s..." << endl;
#endif
	l1l2s_.resize(bases_.size());
	for (ii i = 0; i < (ii)bases_.size(); i++)
	{
		if (i == 0)
		{
			Matrix t; t.init(g_.m(), g_.n(), 1.0);
			bases_[i]->analysis(l1l2s_[0], t, false);
		}
		else
		{
			bases_[i]->analysis(l1l2s_[i], l1l2s_[bases_[i]->getParentIndex()], false);
		}
	}
	for (ii i = 0; i < (ii)bases_.size(); i++)
	{
		if (!bases_[i]->isTransient())
		{
			l1l2s_[i].elementwiseDiv(l1l2s_[i], l2s_[i]);
		}
		else
		{
			l1l2s_[i].free();
		}
	}

	// initialise starting estimate of 'c' from analysis of 'g'
#ifndef NDEBUG
	cout << " Init Cs from G" << endl;
#endif
	double sumC = 0.0;
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
		sumC += cs_[i].sum();
	}
	double sumG = g_.sum();
	for (ii i = 0; i < (ii)bases_.size(); i++)
	{
		if (!bases_[i]->isTransient())
		{
			Matrix l1; l1.elementwiseMul(l1l2s_[i], l2s_[i]);
			cs_[i].elementwiseDiv(cs_[i], l1);
			cs_[i].elementwiseMul((fp)(sumG / sumC), cs_[i]);
			cs_[i].elementwiseMul(cs_[i], l2s_[i]);
			cs_[i].prune(pruneThreshold);
		}
		else
		{
			cs_[i].free();
		}
	}

	// how much memory are we using?
	li mem = 0;
	for (ii i = 0; i < cs_.size(); i++) mem += cs_[i].mem();
	for (ii i = 0; i < l2s_.size(); i++) mem += l2s_[i].mem();
	for (ii i = 0; i < l1l2s_.size(); i++) mem += l1l2s_[i].mem();
	cout << "Sparse Richardson-Lucy mem=";
    cout << cout.unsetf(std::ios::floatfield);
    cout << setprecision(3) << mem / (1024.0*1024.0) << "Mb" << endl;
    //cout << "Sparse Richardson-Lucy mem=" << defaultfloat << setprecision(3) << mem / (1024.0*1024.0) << "Mb" << endl;
}


OptimizerSrl::~OptimizerSrl()
{
}


void OptimizerSrl::init(fp lambda)
{
	lambda_ = lambda;
	iteration_ = 0;
}


double OptimizerSrl::step()
{
	iteration_++;

	// temporary variables
	Matrix f;
	vector<Matrix> cEs(bases_.size());
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

	double volG = g_.sum() / g_.size();
	double volF = f.sum() / f.size();
	cout << "  volF=" << volF << " volG=" << volG << " vol=" << volF / volG << endl;

    // ERROR
#ifndef NDEBUG
	cout << iteration_ << " error" << endl;
#endif
	double errorStart = omp_get_wtime();
	{
		f.elementwiseDiv(g_, f);
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
	for (ii i = 0; i < (ii)bases_.size(); i++)
	{
		if (!bases_[i]->isTransient())
		{
			cEs[i].elementwiseDiv(cEs[i], l2s_[i]);
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
				bases_[i]->shrinkage(cEs[i], cEs[i], cs_[i], l1l2s_[i], lambda_);
			}
			else
			{
				cEs[i].free();
			}
		}
	}
    double shrinkageDuration = omp_get_wtime() - shrinkageStart;
   
    // UPDATE
#ifndef NDEBUG
	cout << iteration_ << " acceleration" << endl;
#endif
	double accelerationStart = omp_get_wtime();
	{
		// gradient calculation is based on unaccelerated update so we can terminate at the same point regardless 
		for (ii i = 0; i < (ii)bases_.size(); i++)
		{
			if (!bases_[i]->isTransient())
			{
				sumSqrs += cs_[i].sumSqrs();
				sumSqrDiffs += cEs[i].sumSqrDiffs(cs_[i]);
			}
		}

		// possible accelerated update
		update(cs_, cEs, iteration_);

		// prune small coefficients
		for (ii i = 0; i < (ii)bases_.size(); i++)
		{
			if (!bases_[i]->isTransient())
			{
				cs_[i].prune(pruneThreshold_);
			}
		}
	}
	double accelerationDuration = omp_get_wtime() - accelerationStart;

#ifndef NDEBUG
	cout << "Durations: synthesis=";
    cout << cout.unsetf(std::ios::floatfield);
    cout << setprecision(3) << synthesisDuration;
    //cout << "Durations: synthesis=" << defaultfloat << setprecision(3) << synthesisDuration;
	cout << " error=" << errorDuration;
	cout << " analysis=" << analysisDuration;
	cout << " shrinkage=" << shrinkageDuration;
	cout << " acceleration=" << accelerationDuration;
	cout << " all=" << synthesisDuration + errorDuration + analysisDuration + shrinkageDuration + accelerationDuration << endl;
#endif

	return sqrt(sumSqrDiffs) / sqrt(sumSqrs);
}


// unaccelerated update
void OptimizerSrl::update(std::vector<Matrix>& cs, std::vector<Matrix>& c1s, ii iteration)
{
	for (ii i = 0; i < (ii)bases_.size(); i++)
	{
		if (!bases_[i]->isTransient())
		{
			cs[i].copy(c1s[i]);
			c1s[i].free();
		}
	}
}


void OptimizerSrl::synthesis(Matrix& f, ii basis) const
{
	vector<Matrix> ts(bases_.size());
	for (ii i = (ii)bases_.size() - 1; i >= 0; i--)
	{
		if (!ts[i] && !bases_[i]->isTransient())
		{
			ts[i].elementwiseDiv(cs_[i], l2s_[i]);
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
				ts[pi].elementwiseDiv(cs_[pi], l2s_[pi]);
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


const std::vector<Matrix>& OptimizerSrl::getCs() const
{ 
	return cs_;
}

