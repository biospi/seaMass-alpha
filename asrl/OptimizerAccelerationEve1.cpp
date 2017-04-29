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
#include <iomanip>
using namespace std;
using namespace kernel;


OptimizerAccelerationEve1::OptimizerAccelerationEve1(Optimizer* optimizer) : optimizer_(optimizer), accelerationDuration_(0.0)
{
    if (getDebugLevel() % 10 >= 1)
        cout << getTimeStamp() << "  Initialising Biggs-Andrews Acceleration (EVE1) ..." << endl;

    // temporaries required for acceleration
    x0s_.resize(xs().size());
    y0s_.resize(xs().size());
    u0s_.resize(xs().size());
    for (size_t i = 0; i < xs().size(); i++)
    {
        x0s_[i].resize(xs()[i].size());
        y0s_[i].resize(xs()[i].size());
        u0s_[i].resize(xs()[i].size());
    }
}


OptimizerAccelerationEve1::~OptimizerAccelerationEve1()
{
}


void OptimizerAccelerationEve1::setLambda(fp lambda,fp lambdaGroup)
{
    optimizer_->setLambda(lambda, lambdaGroup);
}


fp OptimizerAccelerationEve1::step()
{
    if (getDebugLevel() % 10 >= 3)
        cout << getTimeStamp() << "    Acceleration ..." << endl;

    double accelerationStart = getElapsedTime();

    fp a = 0.0;
    if (getIteration() == 0)
    {
        for (ii i = 0; i < (ii)getBases().size(); i++)
        {
            if (!getBases()[i]->isTransient())
            {
                // no extrapolation this iteration, just save 'xs'
                for (size_t k = 0; k < xs()[i].size(); k++)
                {
                    y0s_[i][k].copy(xs()[i][k]);
               }
            }
        }
    }
    else if (getIteration() == 1)
    {
        for (ii i = 0; i < (ii)getBases().size(); i++)
        {
            if (!getBases()[i]->isTransient())
            {
                for (size_t k = 0; k < xs()[i].size(); k++)
                {
                    // can now calcaulte first gradient vector 'u0s'
                    MatrixSparse t;
                    t.copy(xs()[i][k]);
                    t.copySubset(y0s_[i][k]);
                    
                    u0s_[i][k].copy(xs()[i][k]);
                    u0s_[i][k].divNonzeros(t);
                    // no extrapolation this iteration, just save 'xs'
                    x0s_[i][k].copy(xs()[i][k]);
                    y0s_[i][k].copy(xs()[i][k]);
                }
            }
        }
    }
    else
    {
        // calculate acceleration parameter 'a'
        double numerator = 0.0;
        double denominator = 0.0;
        for (ii i = 0; i < (ii)getBases().size(); i++)
        {
            if (!getBases()[i]->isTransient())
            {
                for (size_t k = 0; k < xs()[i].size(); k++)
                {
                    // using old gradient vector 'u0s'
                    MatrixSparse cLogU0;
                    cLogU0.copy(u0s_[i][k]);
                    cLogU0.lnNonzeros();
                    cLogU0.mul(x0s_[i][k]); // (x[k-1] . log u[k-2])
                    denominator += cLogU0.sumSqrs();  // (x[k-1] . log u[k-2]) T (x[k-1] . log u[k-2])
                    
                    // update to new gradient vector 'u0s'
                    u0s_[i][k].copy(xs()[i][k]);
                    MatrixSparse t;
                    t.copy(xs()[i][k]);
                    t.copySubset(y0s_[i][k]);
                    u0s_[i][k].divNonzeros(t);
                    
                    t.copy(xs()[i][k]);
                    t.copySubset(cLogU0);
                    
                    // using new gradient vector 'u0s'
                    MatrixSparse c1LogU;
                    c1LogU.copy(u0s_[i][k]);
                    c1LogU.lnNonzeros();
                    c1LogU.mul(xs()[i][k]); // (x[k] . log u[k-1])
                    c1LogU.mul(t); // (x[k] . log u[k-1]) . (x[k-1] . log u[k-2])
                    numerator += c1LogU.sum(); // (x[k] . log u[k-1]) T (x[k-1] . log u[k-2])
                }
            }
        }
        a = (fp)(numerator / denominator);
        fp aThresh = a > 0.0f ? a : 0.0f;
        aThresh = aThresh < 1.0f ? aThresh : 1.0f;

        // linear extrapolation of 'xs'
        for (ii i = 0; i < (ii)getBases().size(); i++)
        {
            if (!getBases()[i]->isTransient())
            {
                for (size_t k = 0; k < xs()[i].size(); k++)
                {
                    // extrapolate 'xs' and save for next iteration as 'y0s'
                    y0s_[i][k].copy(xs()[i][k]);
                    
                    MatrixSparse t;
                    t.copy(xs()[i][k]);
                    t.copySubset(x0s_[i][k]);
                    
                    y0s_[i][k].divNonzeros(t);
                    y0s_[i][k].pow(aThresh);
                    y0s_[i][k].mul(xs()[i][k]); // x[k] . (x[k] / x[k-1])^a
                    
                    x0s_[i][k].copy(xs()[i][k]); // previous 'xs' saved as 'x0s' for next iteration
                    xs()[i][k].copy(y0s_[i][k]); // extrapolated 'xs' for this iteration
               }
            }
        }
    }
    
    double accelerationDuration = getElapsedTime() - accelerationStart;

    if (getDebugLevel() % 10 >= 3)
        cout << getTimeStamp() << fixed << setprecision(4) <<  "    acceleration       = " << a << endl;


        if (getDebugLevel() % 10 >= 2 && getElapsedTime() != 0.0)
    {
        accelerationDuration_ += accelerationDuration;

        cout << getTimeStamp()  << "    duration_acceleration = " << fixed << setprecision(6) << setw(12) << accelerationDuration << "  total = " << setprecision(4) << setw(12) << accelerationDuration_ << endl;
    }

    // now perform the optimizer iteration on the extrapolated 'xs'
    return optimizer_->step();
}


void OptimizerAccelerationEve1::synthesize(vector<MatrixSparse>& f, vector< vector<MatrixSparse> >& xEs, ii basis)
{
    optimizer_->synthesize(f, xEs, basis);
}


void
OptimizerAccelerationEve1::analyze(std::vector< std::vector<MatrixSparse> > &xEs, std::vector<MatrixSparse> &fE, bool l2, bool l2Normalize)
{
    optimizer_->analyze(xEs, fE, l2, l2Normalize);
}


ii OptimizerAccelerationEve1::getIteration() const
{
    return optimizer_->getIteration();
}


const std::vector<Basis*>& OptimizerAccelerationEve1::getBases() const
{
    return optimizer_->getBases();
}


std::vector< std::vector<MatrixSparse> >& OptimizerAccelerationEve1::xs()
{
    return optimizer_->xs();
}


std::vector< std::vector<MatrixSparse> >& OptimizerAccelerationEve1::l2s()
{
    return optimizer_->l2s();
}


std::vector< std::vector<MatrixSparse> >& OptimizerAccelerationEve1::l1l2s()
{
    return optimizer_->l1l2s();
}


// TODO: quadratic vector extrapolation below
/*
        else if (iteration_ == 1) // linear vector extrapolation this time, but save the qs
        {
            for (ii j = 0; j < (ii)bases_.size(); j++)
            {
                if (!bases_[j]->isTransient())
                {
                    //#pragma omp parallel for reduction(+:sum,sumd)
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
                    //#pragma omp parallel for reduction(+:sum,sumd)
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



