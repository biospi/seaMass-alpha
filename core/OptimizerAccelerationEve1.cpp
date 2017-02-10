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


OptimizerAccelerationEve1::OptimizerAccelerationEve1(Optimizer* optimizer) : optimizer_(optimizer), accelerationDuration_(0.0)
{
    if (getDebugLevel() % 10 >= 1)
    {
        cout << getTimeStamp() << "Initialising Biggs-Andrews Acceleration (EVE1)..." << endl;
    }

	// temporaries required for acceleration
	x0s_.resize(getBases().size());
	y0s_.resize(getBases().size());
	u0s_.resize(getBases().size());
}


OptimizerAccelerationEve1::~OptimizerAccelerationEve1()
{
}


void OptimizerAccelerationEve1::init(fp lambda)
{
	optimizer_->init(lambda);
}


fp OptimizerAccelerationEve1::step()
{
    if (getDebugLevel() % 10 >= 2)
    {
        cout << getTimeStamp() << "  Acceleration ..." << endl;
    }
    
    double accelerationStart = getElapsedTime();

    fp a = 0.0;
	if (getIteration() == 0)
	{
		for (ii i = 0; i < (ii)getBases().size(); i++)
		{
			if (getBases()[i]->getTransient() == Basis::Transient::NO)
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
			if (getBases()[i]->getTransient() == Basis::Transient::NO)
			{
				// can now calcaulte first gradient vector 'u0s'
                u0s_[i].copy(xs()[i]);
				u0s_[i].divCorrespondingNonzeros(y0s_[i]);
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
			if (getBases()[i]->getTransient() == Basis::Transient::NO)
			{
				// using old gradient vector 'u0s'
				MatrixSparse cLogU0;
                cLogU0.copy(u0s_[i]);
				cLogU0.lnNonzeros();
				cLogU0.mul(x0s_[i]); // (x[k-1] . log u[k-2])
				denominator += cLogU0.sumSqrs();  // (x[k-1] . log u[k-2]) T (x[k-1] . log u[k-2])

				// update to new gradient vector 'u0s'
                u0s_[i].copy(xs()[i]);
                MatrixSparse t;
                t.copy(xs()[i]);
                t.subsetElementwiseCopy(y0s_[i]);
				u0s_[i].divCorrespondingNonzeros(t);
                
                t.copy(xs()[i]);
                t.subsetElementwiseCopy(cLogU0);

				// using new gradient vector 'u0s'
				Matrix c1LogU;
                c1LogU.copy(u0s_[i]);
				c1LogU.lnNonzeros();
				c1LogU.mul(xs()[i]); // (x[k] . log u[k-1])
                c1LogU.mul(t); // (x[k] . log u[k-1]) . (x[k-1] . log u[k-2])
				numerator += c1LogU.sum(); // (x[k] . log u[k-1]) T (x[k-1] . log u[k-2])
			}
		}
		a = (fp)(numerator / denominator);
		fp aThresh = a > 0.0f ? a : 0.0f;
		aThresh = aThresh < 1.0f ? aThresh : 1.0f;

		// linear extrapolation of 'xs'
		for (ii i = 0; i < (ii)getBases().size(); i++)
		{
			if (getBases()[i]->getTransient() == Basis::Transient::NO)
			{
				// extrapolate 'xs' and save for next iteration as 'y0s'
                y0s_[i].copy(xs()[i]);
                
                MatrixSparse t;
                t.copy(xs()[i]);
                t.subsetElementwiseCopy(x0s_[i]);
                
				y0s_[i].divCorrespondingNonzeros(t);
				y0s_[i].pow(aThresh);
				y0s_[i].mul(xs()[i]); // x[k] . (x[k] / x[k-1])^a

				x0s_[i].copy(xs()[i]); // previous 'xs' saved as 'x0s' for next iteration
				xs()[i].copy(y0s_[i]); // extrapolated 'xs' for this iteration
			}
		}
	}
    
    double accelerationDuration = getElapsedTime() - accelerationStart;
    
    if (getDebugLevel() % 10 >= 2 && getElapsedTime() != 0.0)
    {
        cout << getTimeStamp() << "    acceleration=" << a << endl;
        
        cout << getTimeStamp() << "    duration=";
        cout.unsetf(ios::floatfield);
        cout << setprecision(3) << accelerationDuration << endl;
        
        accelerationDuration_ += accelerationDuration;
        
        cout << getTimeStamp();
        cout << "    total=";
        cout.unsetf(ios::floatfield);
        cout << setprecision(3) << accelerationDuration_ << endl;
    }

	// now perform the optimizer iteration on the extrapolated 'xs'
	return optimizer_->step();
}

void OptimizerAccelerationEve1::synthesis(MatrixSparse& f, ii basis)
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



