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
#include <math.h>


using namespace std;


OptimizerSrl::OptimizerSrl(const vector<Basis*>& bases, const MatrixSparse& b, ii debugLevel, fp pruneThreshold) : bases_(bases), b_(b), pruneThreshold_(pruneThreshold), lambda_(0.0), iteration_(0)
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
			MatrixSparse t;
			t.init(b_.m(), b_.n(), (fp)1.0);
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
			l2s_[i].elementwiseSqrt();
		}
	}

	// compute L1 norm of each L2 normalised basis function and store in 'l1l2s'
#ifndef NDEBUG
	cout << " Init L2L1s..." << endl;
#endif
	l1l2s_.resize(bases_.size());
	for (ii i = 0; i < (ii)bases_.size(); i++)
	{
		if (i == 0)
		{
			MatrixSparse t;
			t.init(b_.m(), b_.n(), (fp)1.0);
			bases_[i]->analysis(l1l2s_[0], t, false);
		}
		else
		{
			bases_[i]->analysis(l1l2s_[i], l1l2s_[bases_[i]->getParentIndex()], false);
		}
	}
	for (ii i = 0; i < (ii)bases_.size(); i++)
	{
		if (bases_[i]->isTransient())
		{
			l1l2s_[i].free();
		}
		else
		{
			l1l2s_[i].elementwiseDiv(l2s_[i]);
		}
	}

	// initialise starting estimate of 'x' from analysis of 'b'
#ifndef NDEBUG
	cout << " Init Cs from G" << endl;
#endif
	double sumX = 0.0;
	xs_.resize(bases_.size());
	for (ii i = 0; i < (ii)bases_.size(); i++)
	{
		if (i == 0)
		{
			bases_[i]->analysis(xs_[0], b_, false);
		}
		else
		{
			bases_[i]->analysis(xs_[i], xs_[bases_[i]->getParentIndex()], false);
		}
		sumX += xs_[i].sum();
	}
	double sumB = b_.sum();
	for (ii i = 0; i < (ii)bases_.size(); i++)
	{
		if (bases_[i]->isTransient())
		{
			xs_[i].free();
		}
		else
		{
			MatrixSparse l1;
			l1.init(l1l2s_[i]);
			l1.elementwiseMul(l2s_[i]);
			MatrixSparse x;
			x.init(xs_[i]);
			x.elementwiseDiv(l1);
			x.elementwiseMul((fp)(sumB / sumX));
			x.elementwiseMul(l2s_[i]);
			xs_[i].init(x, pruneThreshold);
		}
	}

#ifndef NDEBUG
	// how much memory are we using?
	li mem = 0;
	for (ii i = 0; i < xs_.size(); i++) mem += xs_[i].mem();
	for (ii i = 0; i < l2s_.size(); i++) mem += l2s_[i].mem();
	for (ii i = 0; i < l1l2s_.size(); i++) mem += l1l2s_[i].mem();
	cout << "Sparse Richardson-Lucy mem="; 
	cout.unsetf(ios::floatfield);
	cout << setprecision(3) << mem / (1024.0*1024.0) << "Mb" << endl;
#endif
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
	MatrixSparse f, e;
	vector<MatrixSparse> xEs(bases_.size());
	double sumSqrs = 0.0;
	double sumSqrDiffs = 0.0;
	double sumSqrs2 = 0.0;

/*	// SYNTHESIS
#ifndef NDEBUG
	cout << iteration_ << " synthesis" << endl;
#endif
	double synthesisStart = omp_get_wtime();
	{
		synthesis(f);
	}
	double synthesisDuration = omp_get_wtime() - synthesisStart;
	*/
/*	// ERROR
#ifndef NDEBUG
	cout << iteration_ << " error" << endl;
#endif
	double errorStart = omp_get_wtime();
	{
		e.init(b_);
		e.elementwiseDiv(f);
		f.free();
	}
	double errorDuration = omp_get_wtime() - errorStart;
	*/
	// ANALYSIS

	
	e.init(b_.m(), b_.n(), 1.0);


#ifndef NDEBUG
	cout << iteration_ << " analysis" << endl;
#endif
	double analysisStart = omp_get_wtime();
	{
		for (ii i = 0; i < (ii)bases_.size(); i++)
		{
			if (i == 0)
			{
				bases_.front()->analysis(xEs[0], e, false);
				e.free();
			}
			else
			{
				bases_[i]->analysis(xEs[i], xEs[bases_[i]->getParentIndex()], false);
				/*double sum1 = xEs[i].sum();
				double sum1b = 0.0;
				for (ii nz = 0; nz < xEs[i].nnz(); nz++)
				{
					cout << xEs[i].js_[nz] << "," << xEs[i].vs_[nz] << endl;
					sum1b += xEs[i].vs_[nz];
				}
				bases_[i]->analysis(xEs[i], xEs[bases_[i]->getParentIndex()], false);
				double sum2 = xEs[i].sum();
				double sum2b = 0.0;
				for (ii nz = 0; nz < xEs[i].nnz(); nz++)
				{
					cout << xEs[i].js_[nz] << "," << xEs[i].vs_[nz] << endl;
					sum2b += xEs[i].vs_[nz];
				}
				cout << sum1 << "," << sum2 << endl;
				cout << sum1b << "," << sum2b << endl;
				if (sum1 != sum2) exit(0);*/
			}
			cout << i << ":" << (xEs[i].sum()) << endl;
		}
	}
	for (ii i = 0; i < (ii)bases_.size(); i++)
	{
		if (!bases_[i]->isTransient())
		{
			xEs[i].elementwiseDiv(l2s_[i]);
		}
	}
	double analysisDuration = omp_get_wtime() - analysisStart;

/*	// SHRINKAGE
#ifndef NDEBUG
	cout << iteration_ << " shrinkage" << endl;
#endif
	double shrinkageStart = omp_get_wtime();
	{
		for (ii i = 0; i < (ii)bases_.size(); i++)
		{
			if (!bases_[i]->isTransient())
			{
				bases_[i]->shrinkage(xEs[i], xs_[i], l1l2s_[i], lambda_);
			}
			else
			{
				xEs[i].free();
			}
		}
	}
	double shrinkageDuration = omp_get_wtime() - shrinkageStart;
	*/
	// UPDATE
#ifndef NDEBUG
	cout << iteration_ << " termination check" << endl;
#endif
	double updateStart = omp_get_wtime();
	{
		// termination check
		for (ii i = 0; i < (ii)bases_.size(); i++)
		{
			if (!bases_[i]->isTransient())
			{
				sumSqrs += xs_[i].sumSqrs();
				sumSqrDiffs += xs_[i].sumSqrDiffs(xEs[i]);
				sumSqrs2 += xEs[i].sumSqrs();
			}
		}

		// copy into xs_, pruning small coefficients
		for (ii i = 0; i < (ii)bases_.size(); i++)
		{
			if (!bases_[i]->isTransient())
			{
				//xs_[i].init(xEs[i], pruneThreshold_);
				xs_[i].init(xEs[i]);
				xEs[i].free();
			}
		}
	}
	double updateDuration = omp_get_wtime() - updateStart;

/*#ifndef NDEBUG
	cout << "Durations: synthesis=";
	cout.unsetf(ios::floatfield);
	cout << setprecision(3) << synthesisDuration;
	cout << " error=" << errorDuration;
	cout << " analysis=" << analysisDuration;
	cout << " shrinkage=" << shrinkageDuration;
	cout << " update=" << updateDuration;
	cout << " all=" << synthesisDuration + errorDuration + analysisDuration + shrinkageDuration + updateDuration << endl;
#endif*/

	cout << sumSqrs << "," << sumSqrs2 << "," << sumSqrDiffs << endl;
	return sqrt(sumSqrDiffs) / sqrt(sumSqrs);
}


void OptimizerSrl::synthesis(MatrixSparse& f, ii basis) const
{
	vector<MatrixSparse> ts(bases_.size());
	for (ii i = (ii)bases_.size() - 1; i >= 0; i--)
	{
		if (!ts[i] && !bases_[i]->isTransient())
		{
			ts[i].init(xs_[i]);
			ts[i].elementwiseDiv(l2s_[i]);
		}

		if (basis == i) // return with B-spline control points
		{
			f.init(ts[i]);

			break;
		}

		if (i > 0)
		{
			ii pi = bases_[i]->getParentIndex();
			if (!ts[pi] && !bases_[pi]->isTransient())
			{
				ts[pi].init(xs_[pi]);
				ts[pi].elementwiseDiv(l2s_[pi]);
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


ii OptimizerSrl::getIteration() const
{
	return iteration_;
}


const std::vector<Basis*>& OptimizerSrl::getBases() const
{
	return bases_;
}

std::vector<MatrixSparse>& OptimizerSrl::xs()
{
	return xs_;
}
