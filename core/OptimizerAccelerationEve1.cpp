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


#include "OptimizerAccelerationEve1.hpp"
#include <algorithm>
#include <iostream>
#include <iomanip>


using namespace std;


/*OptimizerAccelerationEve1::OptimizerAccelerationEve1(Optimizer* optimizer, ii debugLevel) : optimizer_(optimizer)
{
#ifndef NDEBUG
	cout << "Init Optimizer Biggs-Andrews Acceleration (EVE1)..." << endl;
#endif

	// temporaries required for acceleration
	x0s_.resize(getBases().size());
	y0s_.resize(getBases().size());
	u0s_.resize(getBases().size());
	for (ii i = 0; i < (ii)getBases().size(); i++)
	{
		if (!getBases()[i]->isTransient())
		{
			x0s_[i].init(getBases()[i]->getM(), getBases()[i]->getN());
			y0s_[i].init(getBases()[i]->getM(), getBases()[i]->getN());
			u0s_[i].init(getBases()[i]->getM(), getBases()[i]->getN());
		}
	}

#ifndef NDEBUG
	// how much memory are we using?
	li mem = 0;
	for (ii i = 0; i < x0s_.size(); i++) mem += x0s_[i].mem();
	for (ii i = 0; i < y0s_.size(); i++) mem += y0s_[i].mem();
	for (ii i = 0; i < u0s_.size(); i++) mem += u0s_[i].mem();

	cout << "Eve1-Andrews Acceleration mem=";
	cout.unsetf(ios::floatfield); 
	cout << setprecision(3) << mem / (1024.0*1024.0) << "Mb" << endl;
#endif
}


OptimizerAccelerationEve1::~OptimizerAccelerationEve1()
{
}


void OptimizerAccelerationEve1::init(fp lambda)
{
	optimizer_->init(lambda);
}


double OptimizerAccelerationEve1::step()
{
#ifndef NDEBUG
	cout << getIteration() << " acceleration" << endl;
	fp aDebug = 0.0;
#endif
	double start = omp_get_wtime();

	if (getIteration() == 0)
	{
		for (ii i = 0; i < (ii)getBases().size(); i++)
		{
			if (!getBases()[i]->isTransient())
			{
				// no extrapolation this iteration, just save 'xs'
				y0s_[i].copy(xs()[i]);
			}
		}
	}
	else if (getIteration() == 1)
	{
		for (ii i = 0; i < (ii)getBases().size(); i++)
		{
			if (!getBases()[i]->isTransient())
			{
				// can now calcaulte first gradient vector 'u0s'
				u0s_[i].elementwiseDiv(xs()[i], y0s_[i]);
				// no extrapolation this iteration, just save 'xs'
				x0s_[i].copy(xs()[i]);
				y0s_[i].copy(xs()[i]);
			}
		}
	}
	else
	{
		// calculate acceleration paramater 'a'
		double numerator = 0.0;
		double denominator = 0.0;
		for (ii i = 0; i < (ii)getBases().size(); i++)
		{
			if (!getBases()[i]->isTransient())
			{
				// using old gradient vector 'u0s'
				Matrix cLogU0;
				cLogU0.elementwiseLn(u0s_[i]);
				cLogU0.elementwiseMul(x0s_[i], cLogU0); // (x[k-1] . log u[k-2])
				denominator += cLogU0.sumSqrs();  // (x[k-1] . log u[k-2]) T (x[k-1] . log u[k-2])

				// update to new gradient vector 'u0s'
				u0s_[i].elementwiseDiv(xs()[i], y0s_[i]);

				// using new gradient vector 'u0s'
				Matrix c1LogU;
				c1LogU.elementwiseLn(u0s_[i]);
				c1LogU.elementwiseMul(xs()[i], c1LogU); // (x[k] . log u[k-1])
				c1LogU.elementwiseMul(c1LogU, cLogU0); // (x[k] . log u[k-1]) . (x[k-1] . log u[k-2])
				numerator += c1LogU.sum(); // (x[k] . log u[k-1]) T (x[k-1] . log u[k-2])
			}
		}
		fp a = (fp)(numerator / denominator);
#ifndef NDEBUG
		aDebug = a;
#endif
		a = a > 0.0f ? a : 0.0f;
		a = a < 1.0f ? a : 1.0f;

		// linear extrapolation of 'xs'
		for (ii i = 0; i < (ii)getBases().size(); i++)
		{
			if (!getBases()[i]->isTransient())
			{
				// extrapolate 'xs' and save for next iteration as 'y0s'
				y0s_[i].elementwiseDiv(xs()[i], x0s_[i]);
				y0s_[i].elementwisePow(y0s_[i], a);
				y0s_[i].elementwiseMul(y0s_[i], xs()[i]); // x[k] . (x[k] / x[k-1])^a

				x0s_[i].copy(xs()[i]); // previous 'xs' saved as 'x0s' for next iteration
				xs()[i].copy(y0s_[i]); // extrapolated 'xs' for this iteration
			}
		}
	}

	double duration = omp_get_wtime() - start;

#ifndef NDEBUG
	cout << "Acceleration: a=";
	cout.unsetf(ios::floatfield); 
	cout << setprecision(3) << aDebug << " duration=" << duration << endl;
#endif

	// now perform the optimizer iteration on the extrapolated 'xs'
	return optimizer_->step();
}


void OptimizerAccelerationEve1::synthesis(Matrix& f, ii basis) const
{
	optimizer_->synthesis(f, basis);
}


ii OptimizerAccelerationEve1::getIteration() const
{
	return optimizer_->getIteration();
}


const std::vector<Basis*>& OptimizerAccelerationEve1::getBases() const
{
	return optimizer_->getBases();
}


std::vector<Matrix>& OptimizerAccelerationEve1::xs()
{
	return optimizer_->xs();
}


// quadratic vector extrapolation below (todo)

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



