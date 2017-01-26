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
#include <math.h>


using namespace std;


OptimizerSrl::OptimizerSrl(const vector<Basis*>& bases, const MatrixSparse& b, int debugLevel, fp pruneThreshold) : bases_(bases), pruneThreshold_(pruneThreshold), lambda_(0.0), iteration_(0), debugLevel_(debugLevel), xs_(bases_.size()), l2s_(bases_.size()), l1l2sPlusLambda_(bases_.size()), synthesisDuration_(0.0), errorDuration_(0.0), analysisDuration_(0.0), shrinkageDuration_(0.0), updateDuration_(0.0)
{
    cout << getTimeStamp() << "Initialising Optimizer SRL ..." << endl;

	// compute L2 norm of each basis function and store in 'l2s'
    cout << getTimeStamp() << "Initialising L2 norms ..." << endl;
	for (ii i = 0; i < (ii)bases_.size(); i++)
	{
		if (i == 0)
		{
			MatrixSparse t;
            t.copy(b);
            t.set((fp)1.0);
			bases_[i]->analysis(l2s_[0], t, true);
		}
		else
		{
			bases_[i]->analysis(l2s_[i], l2s_[bases_[i]->getParentIndex()], true);
		}
	}
	for (ii i = 0; i < (ii)bases_.size(); i++)
	{
        if (bases_[i]->getTransient() == Basis::Transient::NO)
		{
            MatrixSparse t;
            l2s_[i].elementwiseSqrt();
            t.prune(l2s_[i], pruneThreshold);
            l2s_[i].copy(t);
		}
		else
		{
            l2s_[i].free();
		}
	}

	// compute L1 norm of each L2 normalised basis function and store in 'l1l2s'
    cout << getTimeStamp() << "Initialising L1 norms of L2 norms ..." << endl;
	for (ii i = 0; i < (ii)bases_.size(); i++)
	{
		if (i == 0)
		{
			MatrixSparse t;
            t.copy(b);
            t.set((fp)1.0);
			bases_[i]->analysis(l1l2sPlusLambda_[0], t, false);
		}
		else
		{
			bases_[i]->analysis(l1l2sPlusLambda_[i], l1l2sPlusLambda_[bases_[i]->getParentIndex()], false);
		}
	}
	for (ii i = 0; i < (ii)bases_.size(); i++)
	{
		if (bases_[i]->getTransient() == Basis::Transient::NO)
		{
            MatrixSparse t;
            t.copy(l2s_[i]);
            t.subsetElementwiseCopy(l1l2sPlusLambda_[i]);
            l1l2sPlusLambda_[i].elementwiseDiv(l2s_[i]);
		}
		else
		{
            l1l2sPlusLambda_[i].free();
		}
	}
    
    // now make 'b_' 'b' but with zeros removed
    b_.prune(b, 0.0);

	// initialise starting estimate of 'x' from analysis of 'b'
    cout << getTimeStamp() << "Seeding from analysis of input ..." << endl;
	//double sumX = 0.0;
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
  
        /*for (ii nz = 0; nz < 1000000; nz++)
        {
            cout << xs_[i].vs_[nz] << ",";
        }
        cout << i << "ARGH1" << endl;*/
        
		//sumX += xs_[i].sum();
	}
	//double sumB = b_.sum();
	for (ii i = 0; i < (ii)bases_.size(); i++)
	{
		if (bases_[i]->getTransient() == Basis::Transient::NO)
		{
            // remove unneeded l1l2sPlusLambda
            MatrixSparse l1l2PlusLambda;
            l1l2PlusLambda.copy(xs_[i]);
            l1l2PlusLambda.subsetElementwiseCopy(l1l2sPlusLambda_[i]);
            
            // normalise and prune xs
            MatrixSparse x;
            x.copy(xs_[i]);
            x.elementwiseDiv(l1l2PlusLambda);
            //x.elementwiseMul((fp)(sumB / sumX));
            xs_[i].prune(x, pruneThreshold);
            
            // remove unneeded l12sPlusLambda again (after pruning)
            MatrixSparse t;
            t.copy(xs_[i]);
            t.subsetElementwiseCopy(l1l2PlusLambda);
            l1l2sPlusLambda_[i].copy(t);
            
            // remove unneeded l2s
            t.copy(xs_[i]);
            t.subsetElementwiseCopy(l2s_[i]);
            l2s_[i].copy(t);
        
            /*for (ii nz = 0; nz < 1000000; nz++)
            {
                cout << xs_[i].vs_[nz] << ",";
            }
            cout << i << "ARGH2" << endl;*/
		}
		else
		{
            xs_[i].free();
		}
	}
    
#ifndef NDEBUG
	cout << getTimeStamp() << "Sparse Richardson-Lucy ..." << endl;
#endif
}


OptimizerSrl::~OptimizerSrl()
{
}


void OptimizerSrl::init(fp lambda)
{
#ifndef NDEBUG
    cout << getTimeStamp() << iteration_ << " lambda=" << lambda << endl;
#endif

    for (ii i = 0; i < (ii)bases_.size(); i++)
    {
        if (bases_[i]->getTransient() == Basis::Transient::NO)
        {
            l1l2sPlusLambda_[i].elementwiseAdd(lambda - lambda_);
        }
    }
    
	lambda_ = lambda;
	iteration_ = 0;
}


fp OptimizerSrl::step()
{
	iteration_++;

	// SYNTHESIS
#ifndef NDEBUG
	cout << getTimeStamp() << iteration_ << " synthesis" << endl;
#endif
    MatrixSparse f;
	double synthesisStart = getElapsedTime();
	{
		synthesis(f);
	}
	double synthesisDuration = getElapsedTime() - synthesisStart;
	
	// ERROR
#ifndef NDEBUG
	cout << getTimeStamp() << iteration_ << " error" << endl;
#endif
    MatrixSparse fE;
	double errorStart = getElapsedTime();
	{
		fE.copy(b_);
		fE.subsetElementwiseDiv(f); // make this dense mm
		f.free();
	}
	double errorDuration = getElapsedTime() - errorStart;
	
	// ANALYSIS
#ifndef NDEBUG
	cout << getTimeStamp() << iteration_ << " analysis" << endl;
#endif
    vector<MatrixSparse> xEs(bases_.size());
	double analysisStart = getElapsedTime();
	{
		for (ii i = 0; i < (ii)bases_.size(); i++)
		{
			if (i == 0)
			{
				bases_.front()->analysis(xEs[0], fE, false);
				fE.free();
			}
			else
			{
				bases_[i]->analysis(xEs[i], xEs[bases_[i]->getParentIndex()], false);
			}
		}
	}
	for (ii i = 0; i < (ii)bases_.size(); i++)
	{
		if (bases_[i]->getTransient() == Basis::Transient::NO)
		{
            MatrixSparse t;
            t.copy(xs_[i]);
            t.subsetElementwiseCopy(xEs[i]);
            xEs[i].copy(t);
            
			xEs[i].elementwiseDiv(l2s_[i]);
		}
	}
	double analysisDuration = getElapsedTime() - analysisStart;

	// SHRINKAGE
#ifndef NDEBUG
	cout << getTimeStamp() << iteration_ << " shrinkage" << endl;
#endif
    vector<MatrixSparse> ys(bases_.size());
	double shrinkageStart = getElapsedTime();
	{
		for (ii i = 0; i < (ii)bases_.size(); i++)
		{
			if (bases_[i]->getTransient() == Basis::Transient::NO)
			{
				bases_[i]->shrinkage(ys[i], xs_[i], xEs[i], l1l2sPlusLambda_[i]);
			}
			
             xEs[i].free();
		}
	}
	double shrinkageDuration = getElapsedTime() - shrinkageStart;
	
	// UPDATE
#ifndef NDEBUG
	cout << getTimeStamp() << iteration_ << " termination check" << endl;
#endif
    fp sumSqrs = 0.0;
    fp sumSqrDiffs = 0.0;
	double updateStart = getElapsedTime();
	{
		// termination check
		for (ii i = 0; i < (ii)bases_.size(); i++)
		{
			if (bases_[i]->getTransient() == Basis::Transient::NO)
			{
				sumSqrs += xs_[i].sumSqrs();
                ys[i].sumSqrs();
				sumSqrDiffs += xs_[i].sumSqrDiffs(ys[i]);
			}
		}

		// copy into xs_, pruning small coefficients
		for (ii i = 0; i < (ii)bases_.size(); i++)
		{
			if (bases_[i]->getTransient() == Basis::Transient::NO)
			{
				xs_[i].prune(ys[i], pruneThreshold_);
                ys[i].free();
                
                MatrixSparse t;
                t.copy(xs_[i]);
                t.subsetElementwiseCopy(l1l2sPlusLambda_[i]);
                l1l2sPlusLambda_[i].copy(t);
                
                t.copy(xs_[i]);
                t.subsetElementwiseCopy(l2s_[i]);
                l2s_[i].copy(t);
			}
		}
	}
	double updateDuration = getElapsedTime() - updateStart;
    
    if (debugLevel_ >= 2 && getElapsedTime() != 0.0)
    {
        cout << getTimeStamp();
        cout << "  Durations: synthesis=";
        cout.unsetf(ios::floatfield);
        cout << setprecision(3) << synthesisDuration;
        cout << " error=" << errorDuration;
        cout << " analysis=" << analysisDuration;
        cout << " shrinkage=" << shrinkageDuration;
        cout << " update=" << updateDuration;
        cout << " all=" << synthesisDuration + errorDuration + analysisDuration + shrinkageDuration + updateDuration << endl;
        
        synthesisDuration_ += synthesisDuration;
        errorDuration_ += errorDuration;
        analysisDuration_ += analysisDuration;
        shrinkageDuration_ += shrinkageDuration;
        updateDuration_ += updateDuration;
        
        cout << getTimeStamp();
        cout << "  Total Durations: synthesis=";
        cout.unsetf(ios::floatfield);
        cout << setprecision(3) << synthesisDuration_;
        cout << " error=" << errorDuration_;
        cout << " analysis=" << analysisDuration_;
        cout << " shrinkage=" << shrinkageDuration_;
        cout << " update=" << updateDuration_;
        cout << " all=" << synthesisDuration_ + errorDuration_ + analysisDuration_ + shrinkageDuration_ + updateDuration_ << endl;
    }

	return sqrt(sumSqrDiffs) / sqrt(sumSqrs);
}


void OptimizerSrl::synthesis(MatrixSparse& f, ii basis)
{
	vector<MatrixSparse> ts(bases_.size());
	for (ii i = (ii)bases_.size() - 1; i >= 0; i--)
	{
		if (!ts[i].m() && bases_[i]->getTransient() == Basis::Transient::NO)
		{
            ts[i].copy(xs_[i]);
			if (l2s_[i].m()) ts[i].elementwiseDiv(l2s_[i]);
		}

		if (basis == i) // return with B-spline control points
		{
			f.copy(ts[i]);

			break;
		}

		if (i > 0)
		{
			ii pi = bases_[i]->getParentIndex();
			if (!ts[pi].m() && bases_[pi]->getTransient() == Basis::Transient::NO)
			{
                ts[pi].copy(xs_[pi]);
				if (l2s_[pi].m()) ts[pi].elementwiseDiv(l2s_[pi]);
			}

            bases_[i]->deleteBasisFunctions(ts[i], 100000);
			bases_[i]->synthesis(ts[pi], ts[i], bases_[pi]->getTransient() == Basis::Transient::NO);
		}
		else
		{
            bases_[0]->deleteBasisFunctions(ts[0], 100000);
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
