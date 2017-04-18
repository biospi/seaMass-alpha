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


#include "SeamassTopdown.hpp"

#include "../asrl/OptimizerSrl.hpp"
#include "../asrl/OptimizerAccelerationEve1.hpp"

#include "BasisChargeDistribution.hpp"

#include <iostream>
#include <iomanip>


using namespace std;


SeamassTopdown::SeamassTopdown(Input& input, ii maxMass, ii binsPerDalton, double shrinkage, double tolerance, ii debugLevel) : shrinkage_(shrinkage), tolerance_(tolerance), iteration_(0), debugLevel_(debugLevel)
{
	init(input, maxMass, binsPerDalton);
}


SeamassTopdown::~SeamassTopdown()
{
	delete optimizer_;
	delete inner_optimizer_;

	for (ii i = 0; i < (ii)bases_.size(); i++)
	{
		delete bases_[i];
	}
}


void SeamassTopdown::notice()
{
	cout << endl;
	cout << "seaMass-TD - Copyright (C) 2016 - biospi Laboratory, University of Bristol, UK" << endl;
	cout << "This program comes with ABSOLUTELY NO WARRANTY." << endl;
	cout << "This is free software, and you are welcome to redistribute it under certain conditions." << endl;
	cout << endl;
}


void SeamassTopdown::init(Input& input, ii maxMass, ii binsPerDalton)
{
	/*BasisChargeDistribution* basisCharge = new BasisChargeDistribution(bases_, input.counts, input.scale, input.offset, maxMass, binsPerDalton);

	b_.init((li)input.counts.size(), 1, input.counts.data());
	innerOptimizer_ = new OptimizerSrl(bases_, b_, debugLevel_);
	optimizer_ = new OptimizerAccelerationEve1(innerOptimizer_, debugLevel_);
	optimizer_->init((fp)shrinkage_);*/
}


bool SeamassTopdown::step()
{
	/*if (iteration_ == 0 && debugLevel_ > 0)
	{
		li nnz = 0;
		li nx = 0;
		for (ii j = 0; j < (ii)bases_.size(); j++)
		{
			if (!static_cast<BasisBspline*>(bases_[j])->isTransient())
			{
				nnz += optimizer_->xs()[j].nnz();
				nx += optimizer_->xs()[j].size();
			}
		}

		cout << " it:     0 nx: " << setw(10) << nx << " nnz: " << setw(10) << nnz << " tol:  " << fixed << setprecision(8) << setw(10) << tolerance_ << endl;
	}

	iteration_++;
	double grad = optimizer_->step();

	if (debugLevel_ > 0)
	{
		li nnz = 0;
		for (ii j = 0; j < (ii)bases_.size(); j++)
		{
			if (!static_cast<BasisBspline*>(bases_[j])->isTransient())
			{
				nnz += optimizer_->xs()[j].nnz();
			}
		}
		cout << " it: " << setw(5) << iteration_;
		cout << " shrink: ";
		cout.unsetf(ios::floatfield); 
		cout << setprecision(4) << setw(6) << shrinkage_;
		cout << " nnz: " << setw(10) << nnz;
		cout << " grad: " << fixed << setprecision(8) << setw(10) << grad << endl;
	}

	if (grad <= tolerance_)
	{
		if (shrinkage_ == 0)
		{
			if (debugLevel_ == 0) cout << "o" << endl;
			return false;
		}
		else
		{
			if (debugLevel_ == 0) cout << "o" << flush;
			shrinkage_ *= (shrinkage_ > 0.0625 ? 0.5 : 0.0);
			optimizer_->init((fp)shrinkage_);
		}
	}
	else
	{
		if (debugLevel_ == 0) cout << "." << flush;
	}*/

	return true;
}


ii SeamassTopdown::getIteration() const
{
	return iteration_;
}
