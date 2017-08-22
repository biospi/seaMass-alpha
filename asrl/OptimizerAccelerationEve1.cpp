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
#include <kernel.hpp>
#include <algorithm>
#include <iomanip>
#include <sstream>
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
    for (ii k = 0; k < ii(xs().size()); k++)
    {
        x0s_[k].resize(xs()[k].size());
        y0s_[k].resize(xs()[k].size());
        u0s_[k].resize(xs()[k].size());
    }
}


OptimizerAccelerationEve1::~OptimizerAccelerationEve1()
{
}


void OptimizerAccelerationEve1::setLambda(fp lambda, fp lambdaGroup)
{
    optimizer_->setLambda(lambda, lambdaGroup);
}


fp OptimizerAccelerationEve1::getLambda() const
{
    return optimizer_->getLambda();
}

fp OptimizerAccelerationEve1::getLambdaGroup() const
{
    return optimizer_->getLambdaGroup();
}


fp OptimizerAccelerationEve1::step()
{
    if (getDebugLevel() % 10 >= 3)
        cout << getTimeStamp() << "    Acceleration ..." << endl;

    fp a = 0.0;

    double accelerationStart = getElapsedTime();
    {
        if (getIteration() == 0)
        {
            for (ii l = 0; l < (ii)getBases().size(); l++)
            {
                if (!getBases()[l]->isTransient())
                {
                    if (getDebugLevel() % 10 >= 3)
                    {
                        ostringstream oss;
                        oss << getTimeStamp() << "     " << l << " OptimizerAccelerationEve1::acceleration0";
                        info(oss.str());
                    }

                    // no extrapolation this iteration, just save 'xs'
                    for (ii k = 0; k < ii(xs()[l].size()); k++)
                        y0s_[l][k].copy(xs()[l][k]);
                }
            }
        }
        else if (getIteration() == 1)
        {
            for (ii l = 0; l < (ii)getBases().size(); l++)
            {
                if (!getBases()[l]->isTransient())
                {
                    for (ii k = 0; k < ii(xs()[l].size()); k++)
                    {
                        if (getDebugLevel() % 10 >= 3)
                        {
                            ostringstream oss;
                            oss << getTimeStamp() << "     " << l << " OptimizerAccelerationEve1::acceleration1";
                            info(oss.str());
                        }
                        // can now calcaulte first gradient vector 'u0s'
                        MatrixSparse t;
                        t.copyAatB(y0s_[l][k], xs()[l][k]);

                        u0s_[l][k].divNonzeros(xs()[l][k], t);
                        // no extrapolation this iteration, just save 'xs'
                        x0s_[l][k].copy(xs()[l][k]);
                        y0s_[l][k].copy(xs()[l][k]);
                    }
                }
            }
        }
        else
        {
            // calculate acceleration parameter 'a'
            double numerator = 0.0;
            double denominator = 0.0;
            for (ii l = 0; l < ii(getBases().size()); l++)
            {
                if (!getBases()[l]->isTransient())
                {
                    if (getDebugLevel() % 10 >= 3)
                    {
                        ostringstream oss;
                        oss << getTimeStamp() << "     " << l << " OptimizerAccelerationEve1::accelerationCalcA";
                        info(oss.str());
                    }

                    for (ii k = 0; k < ii(xs()[l].size()); k++)
                    {
                        // using old gradient vector 'u0s'
                        MatrixSparse cLogU0;
                        cLogU0.lnNonzeros(u0s_[l][k]);
                        cLogU0.mul(cLogU0, x0s_[l][k]); // (x[k-1] . log u[k-2])
                        denominator += cLogU0.sumSqrs();  // (x[k-1] . log u[k-2]) T (x[k-1] . log u[k-2])

                        // update to new gradient vector 'u0s'
                        MatrixSparse t;
                        t.copyAatB(y0s_[l][k], xs()[l][k]);
                        u0s_[l][k].divNonzeros(xs()[l][k], t);

                        t.copyAatB(cLogU0, xs()[l][k]);

                        // using new gradient vector 'u0s'
                        MatrixSparse c1LogU;
                        c1LogU.lnNonzeros(u0s_[l][k]);
                        c1LogU.mul(c1LogU, xs()[l][k]); // (x[k] . log u[k-1])
                        c1LogU.mul(c1LogU, t); // (x[k] . log u[k-1]) . (x[k-1] . log u[k-2])
                        numerator += c1LogU.sum(); // (x[k] . log u[k-1]) T (x[k-1] . log u[k-2])
                    }
                }
            }
            a = fp(numerator / denominator);
            fp aThresh = a > 0.0f ? a : 0.0f;
            aThresh = aThresh < 1.0f ? aThresh : 1.0f;

            // linear extrapolation of 'xs'
            for (ii l = 0; l < ii(getBases().size()); l++)
            {
                if (!getBases()[l]->isTransient())
                {
                    if (getDebugLevel() % 10 >= 3)
                    {
                        ostringstream oss;
                        oss << getTimeStamp() << "     " << l << " OptimizerAccelerationEve1::acceleration2+";
                        info(oss.str());
                    }

                    for (ii k = 0; k < ii(xs()[l].size()); k++)
                    {
                        // extrapolate 'xs' and save for next iteration as 'y0s'
                        y0s_[l][k].copy(xs()[l][k]);

                        MatrixSparse t;
                        t.copyAatB(x0s_[l][k], xs()[l][k]);

                        y0s_[l][k].divNonzeros(y0s_[l][k], t);
                        y0s_[l][k].pow(y0s_[l][k], aThresh);
                        y0s_[l][k].mul(y0s_[l][k], xs()[l][k]); // x[k] . (x[k] / x[k-1])^a

                        x0s_[l][k].copy(xs()[l][k]); // previous 'xs' saved as 'x0s' for next iteration
                        xs()[l][k].copy(y0s_[l][k]); // extrapolated 'xs' for this iteration
                    }
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



