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
#include <kernel.hpp>
#include <iomanip>
#include <cmath>
#include <sstream>
#include <limits>
using namespace std;
using namespace kernel;


OptimizerSrl::OptimizerSrl(const vector<Basis*>& bases, const std::vector<Matrix>& b, bool seed, fp pruneThreshold) : bases_(bases), b_(b), pruneThreshold_(pruneThreshold), lambda_(0.0), lambdaGroup_(0.0), iteration_(0), synthesisDuration_(0.0), errorDuration_(0.0), analysisDuration_(0.0), shrinkageDuration_(0.0), updateDuration_(0.0)
{
    if (getDebugLevel() % 10 >= 1)
        cout << getTimeStamp() << "  Creating optimizer SRL ..." << endl;

    if (seed)
    {
        {   // compute L2 and L1 norm of each basis function and store in 'l2s' and 'l1l2s'
            if (getDebugLevel() % 10 >= 1)
                cout << getTimeStamp() << "  Initialising L2 norms ..." << endl;

            vector<MatrixSparse> t(b_.size());
            for (ii k = 0; k < ii(t.size()); k++)
                t[k].initDense(1, b_[k].n(), fp(1.0));

            analyze(l2s_, t, true, false);

            if (getDebugLevel() % 10 >= 1)
                cout << getTimeStamp() << "  Initialising L1 norms of L2 norms ..." << endl;

            analyze(l1l2sPlusLambda_, t, false);
        }

        {   // initialise starting estimate of 'x' from analysis of 'b'
            if (getDebugLevel() % 10 >= 1)
                cout << getTimeStamp() << "  Seeding from analysis of input ..." << endl;

            vector<MatrixSparse> t(b_.size());
            for (ii k = 0; k < ii(t.size()); k++)
                t[k].importFromMatrix(b_[k]);

            analyze(xs_, t, false, false);

            double sumB = 0.0;
            for (ii k = 0; k < ii(b_.size()); k++)
                sumB += b_[k].sum();

            if (getDebugLevel() % 10 >= 2)
                cout << getTimeStamp() << "    volume_b=" << fixed << sumB << endl;

            double sumX = 0.0;
            for (ii l = 0; l < ii(bases_.size()); l++)
            {
                for (ii k = 0; k < ii(xs_[l].size()); k++)
                    sumX += xs_[l][k].sum();
            }

            for (ii l = 0; l < ii(bases_.size()); l++)
            {
                if (!bases_[l]->isTransient())
                {
                    for (ii k = 0; k < ii(xs_[l].size()); k++)
                    {
                        // remove unneeded l1l2sPlusLambda
                        MatrixSparse l1l2PlusLambda;
                        l1l2PlusLambda.copySubset(l1l2sPlusLambda_[l][k], xs_[l][k]);

                        // normalise and prune xs
                        MatrixSparse x;
                        x.divNonzeros(xs_[l][k], l1l2PlusLambda);
                        l1l2PlusLambda.empty();
                        x.mul((fp) (sumB / sumX));
                        xs_[l][k].pruneCells(x, pruneThreshold);
                        x.empty();

                        // remove unneeded l2s
                        MatrixSparse t;
                        t.copySubset(l2s_[l][k], xs_[l][k]);
                        l2s_[l][k].swap(t);

                        // remove unneeded l1l2s
                        t.copySubset(l1l2sPlusLambda_[l][k], xs_[l][k]);
                        l1l2sPlusLambda_[l][k].swap(t);
                    }
                }
            }
        }
    }
}


OptimizerSrl::~OptimizerSrl()
{
}


void OptimizerSrl::setLambda(fp lambda, fp lambdaGroup)
{
    if (getDebugLevel() % 10 >= 3)
        cout << getTimeStamp() << "   lambda=" << lambda << endl;

    for (ii l = 0; l < ii(bases_.size()); l++)
    {
        if (!bases_[l]->isTransient())
        {
            for (ii k = 0; k < ii(l1l2sPlusLambda_[l].size()); k++)
                l1l2sPlusLambda_[l][k].addNonzeros(lambda - lambda_);
        }
    }

    lambda_ = lambda;
    lambdaGroup_ = lambdaGroup;
    iteration_ = 0;
}


fp OptimizerSrl::step()
{
    iteration_++;

    // SYNTHESISE
    if (getDebugLevel() % 10 >= 3)
        cout << getTimeStamp() << "    Synthesis ..." << endl;

    vector<MatrixSparse> f_fE;
    vector< vector<MatrixSparse> > xEs_ys;
    double synthesisStart = getElapsedTime();
    {
        synthesize(f_fE, xEs_ys);

        if (getDebugLevel() % 10 >= 2)
        {
            double sumF = 0.0;
            for (ii k = 0; k < ii(b_.size()); k++)
                sumF += f_fE[k].sum();
            cout << getTimeStamp() << "    volume_f=" << fixed << sumF << endl;
        }
    }
    double synthesisDuration = getElapsedTime() - synthesisStart;

    // ERROR
    if (getDebugLevel() % 10 >= 3)
        cout << getTimeStamp() << "    Error ..." << endl;

    double errorStart = getElapsedTime();
    {
        for (ii k = 0; k < ii(f_fE.size()); k++)
        {
            // any zeros in f_fE are due to underflow. We need to do this to avoid divide by zero error
            f_fE[k].censorLeft(f_fE[k], numeric_limits<fp>::min());
            f_fE[k].div2Dense(b_[k]);

            MatrixSparse t;
            t.pruneCells(f_fE[k]);
            f_fE[k].swap(t);
         }
    }
    double errorDuration = getElapsedTime() - errorStart;

    // ANALYSIS
    if (getDebugLevel() % 10 >= 3)
        cout << getTimeStamp() << "    Analysis ..." << endl;

    double analysisStart = getElapsedTime();
    {
        analyze(xEs_ys, f_fE, false);

        for (ii k = 0; k < ii(f_fE.size()); k++)
            f_fE[k].empty();
    }
    double analysisDuration = getElapsedTime() - analysisStart;

    // SHRINKAGE
    if (getDebugLevel() % 10 >= 3)
        cout << getTimeStamp() << "    Shrinkage ..." << endl;

    double shrinkageStart = getElapsedTime();
    {
        for (ii l = 0; l < ii(bases_.size()); l++)
        {
            if (!bases_[l]->isTransient())
            {
                if (lambdaGroup_ > 0.0)
                {
                    if (getDebugLevel() % 10 >= 3)
                    {
                        ostringstream oss;
                        oss << getTimeStamp() << "     " << l << " OptimizerSrl::shrinkageGroup";
                        info(oss.str());
                    }

                    vector<MatrixSparse>* g = bases_[l]->getGroups(false);
                    vector<MatrixSparse>* gT = bases_[l]->getGroups(true);

                    // group and individual shrinkage
                    for (ii k = 0; k < ii(xEs_ys[l].size()); k++)
                    {
                        // y = groupNorm(x)
                        MatrixSparse t;
                        t.sqr(xs_[l][k]);
                        MatrixSparse y;
                        y.matmul(false, t, (*gT)[k], false);
                        t.matmul(false, y, (*g)[k], false);
                        y.copySubset(t, xs_[l][k]);
                        t.empty();
                        y.sqrt(y);

                        // y = x * groupNorm(x)^-1)
                        y.divNonzeros(xs_[l][k], y);

                        // y = lambdaGroup * x * groupNorm(x)^-1
                        y.mul(lambdaGroup_);

                        // y = l1l2 + lambda + lambdaGroup * x * groupNorm(x)^-1
                        y.addNonzeros(y, l1l2sPlusLambda_[l][k]);

                        // y = x / (l1l2 + lambda + lambdaGroup * x * groupNorm(x)^-1)
                        y.divNonzeros(xs_[l][k], y);

                        // y = xE * x / (l1l2 + lambda + lambdaGroup * x * groupNorm(x)^-1)
                        xEs_ys[l][k].mul(xEs_ys[l][k], y);
                    }
                }
                else
                {
                    if (getDebugLevel() % 10 >= 3)
                    {
                        ostringstream oss;
                        oss << getTimeStamp() << "     " << l << " OptimizerSrl::shrinkage";
                        info(oss.str());
                    }

                    // individual shrinkage only
                    for (ii k = 0; k < ii(xEs_ys[l].size()); k++)
                    {
                        // y = x / (l1l2 + lambda)
                        MatrixSparse y;
                        y.divNonzeros(xs_[l][k], l1l2sPlusLambda_[l][k]);

                        // y = xE * x / (l1l2 + lambda)
                        xEs_ys[l][k].mul(xEs_ys[l][k], y);
                    }
                }
            }
        }
    }
    double shrinkageDuration = getElapsedTime() - shrinkageStart;

    // UPDATE
    if (getDebugLevel() % 10 >= 3)
        cout << getTimeStamp() << "    Termination Check and Pruning..." << endl;

    fp sumSqrs = 0.0;
    fp sumSqrs2 = 0.0;
    fp sumSqrDiffs = 0.0;
    double updateStart = getElapsedTime();
    {
        // termination check
        for (ii l = 0; l < ii(bases_.size()); l++)
        {
            if (!bases_[l]->isTransient())
            {
                if (getDebugLevel() % 10 >= 3)
                {
                    ostringstream oss;
                    oss << getTimeStamp() << "     " << l << " OptimizerSrl::updateGradient";
                    info(oss.str());
                }

                for (ii k = 0; k < ii(xs_[l].size()); k++)
                {
                    sumSqrs += xs_[l][k].sumSqrs();
                    sumSqrs2 += xEs_ys[l][k].sumSqrs();
                    sumSqrDiffs += xs_[l][k].sumSqrDiffsNonzeros(xEs_ys[l][k]);
                }
            }
        }

        // copy into xs_, pruning small coefficients
        for (ii l = 0; l < ii(bases_.size()); l++)
        {
            if (!bases_[l]->isTransient())
            {
                if (getDebugLevel() % 10 >= 3)
                {
                    ostringstream oss;
                    oss << getTimeStamp() << "     " << l << " OptimizerSrl::updatePrune";
                    info(oss.str());
                }

                for (ii k = 0; k < ii(xs_[l].size()); k++)
                {
                    xs_[l][k].pruneCells(xEs_ys[l][k], pruneThreshold_);
                    xEs_ys[l][k].empty();

                    // prune l1l2s
                    MatrixSparse t;
                    t.copySubset(l1l2sPlusLambda_[l][k], xs_[l][k]);
                    l1l2sPlusLambda_[l][k].swap(t);

                    // prune l2s
                    t.copySubset(l2s_[l][k], xs_[l][k]);
                    l2s_[l][k].swap(t);
                }
            }
        }
    }
    double updateDuration = getElapsedTime() - updateStart;
    
    if (getDebugLevel() % 10 >= 2 && getElapsedTime() != 0.0)
    {
        synthesisDuration_ += synthesisDuration;
        errorDuration_ += errorDuration;
        analysisDuration_ += analysisDuration;
        shrinkageDuration_ += shrinkageDuration;
        updateDuration_ += updateDuration;

        cout << getTimeStamp()  << "    duration_synthesis    = " << fixed << setprecision(6) << setw(12) << synthesisDuration << "  total = " << setprecision(4) << setw(12) << synthesisDuration_ << endl;
        cout << getTimeStamp()  << "    duration_error        = " << fixed << setprecision(6) << setw(12) << errorDuration << "  total = " << setprecision(4) << setw(12) << errorDuration_ << endl;
        cout << getTimeStamp()  << "    duration_analysis     = " << fixed << setprecision(6) << setw(12) << analysisDuration << "  total = " << setprecision(4) << setw(12) << analysisDuration_ << endl;
        cout << getTimeStamp()  << "    duration_shrinkage    = " << fixed << setprecision(6) << setw(12) << shrinkageDuration << "  total = " << setprecision(4) << setw(12) << shrinkageDuration_ << endl;
        cout << getTimeStamp()  << "    duration_update       = " << fixed << setprecision(6) << setw(12) << updateDuration << "  total = " << setprecision(4) << setw(12) << updateDuration_ << endl;
        cout << getTimeStamp()  << "                                     total_sort = " << setprecision(4) << setw(12) << MatrixSparse::sortElapsed_ << endl;
     }

    return sqrt(sumSqrDiffs) / sqrt(sumSqrs);
}


void OptimizerSrl::synthesize(vector<MatrixSparse>& f, vector< vector<MatrixSparse> >& xEs, ii basis)
{
    if (xEs.size() != bases_.size())
        xEs.resize(bases_.size());

    for (ii l = ii(bases_.size()) - 1; l >= 0; l--)
    {
        if (!xEs[l].size() && !bases_[l]->isTransient())
        {
            xEs[l].resize(xs_[l].size());

            if (l2s_.size())
            {
                for (ii k = 0; k < ii(xEs[l].size()); k++)
                    xEs[l][k].divNonzeros(xs_[l][k], l2s_[l][k]);
            }
        }

        if (basis == l) // return with B-spline control points
        {
            for (ii k = 0; k < ii(xEs[l].size()); k++)
                f[k].swap(xEs[l][k]);

            break;
        }

        if (l > 0)
        {
            ii pi = bases_[l]->getParentIndex();
            
            if (!xEs[pi].size() && !bases_[pi]->isTransient())
            {
                xEs[pi].resize(xs_[pi].size());

                if (l2s_.size())
                {
                    for (ii k = 0; k < ii(xEs[pi].size()); k++)
                        xEs[pi][k].divNonzeros(xs_[pi][k], l2s_[pi][k]);
                }
            }

            bases_[l]->synthesize(xEs[pi], xEs[l], !bases_[pi]->isTransient());
        }
        else
        {
            bases_[0]->synthesize(f, xEs[0], false);
        }
    }
}


void OptimizerSrl::analyze(std::vector<std::vector<MatrixSparse> > &xEs, std::vector<MatrixSparse> &fE, bool l2, bool l2Normalize) const
{
    if (xEs.size() != bases_.size())
        xEs.resize(bases_.size());

    vector<MatrixSparse> t;
    bases_.front()->analyze(t, fE, l2);

    if (xEs[0].size() != t.size())
    {
        xEs[0].resize(t.size());
        for (ii k = 0; k < ii(xEs[0].size()); k++)
            xEs[0][k].swap(t[k]);
    }
    else
    {
        for (ii k = 0; k < ii(xEs[0].size()); k++)
        {
            xEs[0][k].swap(t[k]);
            /*xEs[0][k].copySubset(t[k]);
            if (t[k].nnz() != xEs[0][k].nnz())
            {
                cout << "input  0[" << k << "]: " << t[k] << endl;
                cout << "output 0[" << k << "]: " << xEs[0][k] << endl << endl;
            }*/
       }
    }

    for (ii l = 1; l < ii(bases_.size()); l++)
    {
        vector<MatrixSparse> t;
        bases_[l]->analyze(t, xEs[bases_[l]->getParentIndex()], l2);

        if (xEs[l].size() != t.size())
        {
            xEs[l].resize(t.size());
            for (ii k = 0; k < ii(xEs[l].size()); k++)
                xEs[l][k].swap(t[k]);
        }
        else
        {
            for (ii k = 0; k < ii(xEs[l].size()); k++)
            {
                xEs[l][k].swap(t[k]);
                /*xEs[l][k].copySubset(t[k]);
                if (t[k].nnz() != xEs[l][k].nnz())
                {
                    cout << "input  " << l << "[" << k << "]: " << t[k] << endl;
                    cout << "output " << l << "[" << k << "]: " << xEs[l][k] << endl << endl;
                }*/
            }
        }
   }

    for (ii l = 0; l < ii(bases_.size()); l++)
    {
        if (!bases_[l]->isTransient())
        {
            if (xs_.size() > 0)
            {
                for (ii k = 0; k < ii(xEs[l].size()); k++)
                {
                    if (xEs[l][k].nnz() > xs_[l][k].nnz())
                    {
                        MatrixSparse t;
                        t.copySubset(xEs[l][k], xs_[l][k]);
                        xEs[l][k].swap(t);
                    }
                }
            }

            if(l2)
            {
                for (ii k = 0; k < ii(xEs[l].size()); k++)
                    xEs[l][k].sqrt(xEs[l][k]);
            }

            if (l2Normalize)
            {
                for (ii k = 0; k < ii(xEs[l].size()); k++)
                    xEs[l][k].divNonzeros(xEs[l][k], l2s_[l][k]);
            }
        }
        else
        {
            vector<MatrixSparse>().swap(xEs[l]);
        }
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


std::vector< std::vector<MatrixSparse> >& OptimizerSrl::xs()
{
    return xs_;
}


std::vector< std::vector<MatrixSparse> >& OptimizerSrl::l2s()
{
    return l2s_;
}


std::vector< std::vector<MatrixSparse> >& OptimizerSrl::l1l2s()
{
    return l1l2sPlusLambda_;
}
