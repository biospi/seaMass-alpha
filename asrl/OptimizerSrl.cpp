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
#include <iomanip>
#include <cmath>
using namespace std;
using namespace kernel;


OptimizerSrl::OptimizerSrl(const vector<Basis*>& bases, const std::vector<Matrix>& b, bool seed, fp pruneThreshold) : bases_(bases), b_(b), pruneThreshold_(pruneThreshold), lambda_(0.0), lambdaGroup_(0.0), iteration_(0), synthesisDuration_(0.0), errorDuration_(0.0), analysisDuration_(0.0), shrinkageDuration_(0.0), updateDuration_(0.0)
{
    if (getDebugLevel() % 10 >= 1)
        cout << getTimeStamp() << "  Creating optimizer SRL ..." << endl;

    if (seed)
    {
        // compute L2 norm of each basis function and store in 'l2s'
        if (getDebugLevel() % 10 >= 1)
            cout << getTimeStamp() << "  Initialising L2 norms ..." << endl;

        l2s_.resize(bases_.size());
        for (ii i = 0; i < (ii)bases_.size(); i++)
        {
            if (i == 0)
            {
                vector<MatrixSparse> t(b_.size());
                for (size_t k = 0; k < t.size(); k++)
                {
                    t[k].copy(1, b_[k].n(), (fp) 1.0);
                }
                bases_[i]->analyse(l2s_[0], t, true);
            }
            else
            {
                bases_[i]->analyse(l2s_[i], l2s_[bases_[i]->getParentIndex()], true);
            }
        }

        for (ii i = 0; i < (ii)bases_.size(); i++)
        {
            if (!bases_[i]->isTransient())
            {
                for (size_t k = 0; k < l2s_[i].size(); k++)
                {
                    l2s_[i][k].sort();
                    l2s_[i][k].sqrt();
                }
            }
        }

        // compute L1 norm of each L2 normalised basis function and store in 'l1l2s'
        if (getDebugLevel() % 10 >= 1)
            cout << getTimeStamp() << "    Initialising L1 norms of L2 norms ..." << endl;

        l1l2s_.resize(bases_.size());
        for (ii i = 0; i < (ii)bases_.size(); i++)
        {
            if (i == 0)
            {
                vector<MatrixSparse> t(b_.size());
                for (size_t k = 0; k < t.size(); k++)
                {
                    t[k].copy(1, b_[k].n(), (fp) 1.0);
                }
                bases_[i]->analyse(l1l2s_[0], t, false);
            }
            else
            {
                bases_[i]->analyse(l1l2s_[i], l1l2s_[bases_[i]->getParentIndex()], false);
            }
        }

        for (ii i = 0; i < (ii)bases_.size(); i++)
        {
            if (!bases_[i]->isTransient())
            {
                for (size_t k = 0; k < l1l2s_[i].size(); k++)
                {
                    l1l2s_[i][k].sort();
                    l1l2s_[i][k].divNonzeros(l2s_[i][k].vs());
                }
            }
        }


        // initialise starting estimate of 'x' from analysis of 'b'
        if (getDebugLevel() % 10 >= 1)
            cout << getTimeStamp() << "  Seeding from analysis of input ..." << endl;

        xs_.resize(bases_.size());
        for (ii i = 0; i < (ii)bases_.size(); i++)
        {
            if (i == 0)
            {
                vector<MatrixSparse> t(b_.size());
                for (size_t k = 0; k < t.size(); k++)
                    t[k].copy(b_[k]);

                bases_[i]->analyse(xs_[0], t, false);
            }
            else
            {
                bases_[i]->analyse(xs_[i], xs_[bases_[i]->getParentIndex()], false);
            }
        }

        for (ii i = 0; i < (ii)bases_.size(); i++)
        {
            if (!bases_[i]->isTransient())
            {
                for (size_t k = 0; k < xs_[i].size(); k++)
                {
                    // need to sort xs after matmul
                    xs_[i][k].sort();

                    MatrixSparse x;
                    x.copy(xs_[i][k]);
                    xs_[i][k].prune(x, pruneThreshold);

                    // remove unneeded l2s
                    MatrixSparse t;
                    t.copy(xs_[i][k]);
                    t.copySubset(l2s_[i][k]);
                    l2s_[i][k].copy(t);

                    // remove unneeded l1l2s
                    t.copy(xs_[i][k]);
                    t.copySubset(l1l2s_[i][k]);
                    l1l2s_[i][k].copy(t);
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
    vector< vector<MatrixSparse> > xEs;
    double synthesisStart = getElapsedTime();
    {
        synthesise(f_fE, xEs);
    }
    double synthesisDuration = getElapsedTime() - synthesisStart;

    // ERROR
    if (getDebugLevel() % 10 >= 3)
        cout << getTimeStamp() << "    Error ..." << endl;

    double errorStart = getElapsedTime();
    {
        if (getDebugLevel() % 10 >= 3)
            cout << getTimeStamp() << "     OptimizerSrl::error" << endl;

        for (size_t k = 0; k < f_fE.size(); k++)
            f_fE[k].div2Nonzeros(b_[k].vs());
    }
    double errorDuration = getElapsedTime() - errorStart;

    // init l1s_ and l1l2s_
    /*if (l2s_.size() != bases_.size())
    {
     // TODO
    }*/

    // ANALYSIS
    if (getDebugLevel() % 10 >= 3)
        cout << getTimeStamp() << "    Analysis ..." << endl;

    double analysisStart = getElapsedTime();
    {
        for (ii i = 0; i < (ii)bases_.size(); i++)
        {
            if (i == 0)
            {
                vector<MatrixSparse> t;
                bases_.front()->analyse(t, f_fE, false);

                for (size_t k = 0; k < t.size(); k++)
                {
                    xEs[0][k].copy(t[k]);

                    //t[k].sort();
                    //xEs[0][k].copySubset(t[k]);
                }
            }
            else
            {
                vector<MatrixSparse> t;
                bases_[i]->analyse(t, xEs[bases_[i]->getParentIndex()], false);

                for (size_t k = 0; k < t.size(); k++)
                {
                    xEs[i][k].copy(t[k]);

                    //t[k].sort();
                    //xEs[i][k].copySubset(t[k]);
                }
            }
        }
    }

    for (ii i = 0; i < (ii)bases_.size(); i++)
    {
        if (!bases_[i]->isTransient())
        {
            for (size_t k = 0; k < xEs[i].size(); k++)
            {
                xEs[i][k].sort();

                MatrixSparse t;
                t.copy(xs_[i][k]);
                t.copySubset(xEs[i][k]);
                xEs[i][k].copy(t);

                xEs[i][k].divNonzeros(l2s_[i][k].vs());
            }
        }
    }
    double analysisDuration = getElapsedTime() - analysisStart;

    // SHRINKAGE
    if (getDebugLevel() % 10 >= 3)
        cout << getTimeStamp() << "    Shrinkage ..." << endl;

    vector< vector<MatrixSparse> > ys(bases_.size());
    double shrinkageStart = getElapsedTime();
    {
        for (ii i = 0; i < (ii)bases_.size(); i++)
        {
            if (!bases_[i]->isTransient())
            {
                ys[i].resize(xs_[i].size());
                if (lambdaGroup_ > 0.0)
                {
                    const vector<MatrixSparse>* g = bases_[i]->getGroups(false);
                    const vector<MatrixSparse>* gT = bases_[i]->getGroups(true);

                    // group and individual shrinkage
                    for (size_t k = 0; k < ys[i].size(); k++)
                    {

                        // y = groupNorm(x)
                        MatrixSparse t;
                        t.copy(xs_[i][k]);
                        t.sqr();
                        ys[i][k].matmul(false, t, (*gT)[k], false);
                        t.matmul(false, ys[i][k], (*g)[k], false);
                        t.sort();
                        ys[i][k].copy(xs_[i][k]);
                        ys[i][k].copySubset(t);
                        ys[i][k].sqrt();

                        // y = x * groupNorm(x)^-1)
                        ys[i][k].div2Nonzeros(xs_[i][k].vs());

                        // y = lambdaGroup * x * groupNorm(x)^-1
                        ys[i][k].mul(lambdaGroup_);

                        // y = lambda + lambdaGroup * x * groupNorm(x)^-1
                        ys[i][k].addNonzeros(lambda_);

                        // y = l1l2 + lambda + lambdaGroup * x * groupNorm(x)^-1
                        ys[i][k].addNonzeros(l1l2s_[i][k].vs());

                        // y = x / (l1l2 + lambda + lambdaGroup * x * groupNorm(x)^-1)
                        ys[i][k].div2Nonzeros(xs_[i][k].vs());

                        // y = xE * x / (l1l2 + lambda + lambdaGroup * x * groupNorm(x)^-1)
                        ys[i][k].mul(xEs[i][k].vs());
                    }
                }
                else
                {
                    // individual shrinkage only
                    for (size_t k = 0; k < ys[i].size(); k++)
                    {
                        // y = l1l2
                        ys[i][k].copy(l1l2s_[i][k]);

                        // y = l1l2 + lambda
                        ys[i][k].addNonzeros(lambda_);

                        // y = x / (l1l2 + lambda)
                        ys[i][k].div2Nonzeros(xs_[i][k].vs());

                        // y = xE * x / (l1l2 + lambda)
                        ys[i][k].mul(xEs[i][k].vs());
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
    fp sumSqrDiffs = 0.0;
    double updateStart = getElapsedTime();
    {
        // termination check
        for (ii i = 0; i < (ii)bases_.size(); i++)
        {
            if (!bases_[i]->isTransient())
            {
                if (getDebugLevel() % 10 >= 3)
                     cout << getTimeStamp() << "     " << i << " OptimizerSrl::grad" << endl;

                for (size_t k = 0; k < xs_[i].size(); k++)
                {
                    sumSqrs += xs_[i][k].sumSqrs();
                    ys[i][k].sumSqrs();
                    sumSqrDiffs += xs_[i][k].sumSqrDiffsNonzeros(ys[i][k].vs());
                }
            }
        }

        // copy into xs_, pruning small coefficients
        for (ii i = 0; i < (ii)bases_.size(); i++)
        {
            if (!bases_[i]->isTransient())
            {
                for (size_t k = 0; k < xs_[i].size(); k++)
                {
                    xs_[i][k].prune(ys[i][k], pruneThreshold_);
                    ys[i][k].init();
                    
                    MatrixSparse t;
                    t.copy(xs_[i][k]);
                    t.copySubset(l1l2s_[i][k]);
                    l1l2s_[i][k].copy(t);
                    
                    t.copy(xs_[i][k]);
                    t.copySubset(l2s_[i][k]);
                    l2s_[i][k].copy(t);
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
     }

    return sqrt(sumSqrDiffs) / sqrt(sumSqrs);
}


void OptimizerSrl::synthesise(vector<MatrixSparse>& f, vector< vector<MatrixSparse> >& cs, ii basis)
{
    cs.resize(bases_.size());
    for (ii i = (ii)bases_.size() - 1; i >= 0; i--)
    {
        if (!cs[i].size() && !bases_[i]->isTransient())
        {
            cs[i].resize(xs_[i].size());
            for (size_t k = 0; k < cs[i].size(); k++)
            {
                cs[i][k].copy(xs_[i][k]);
                
                if (l2s_.size() == xs_.size() && l2s_[i][k].m())
                    cs[i][k].divNonzeros(l2s_[i][k].vs());
            }
        }

        if (basis == i) // return with B-spline control points
        {
            for (size_t k = 0; k < cs[i].size(); k++)
                f[k].copy(cs[i][k]);

            break;
        }

        if (i > 0)
        {
            ii pi = bases_[i]->getParentIndex();
            
            if (!cs[pi].size() && !bases_[pi]->isTransient())
            {
                cs[pi].resize(xs_[pi].size());
                for (size_t k = 0; k < cs[pi].size(); k++)
                {
                    cs[pi][k].copy(xs_[pi][k]);
                    
                    if (l2s_.size() == xs_.size() && l2s_[pi][k].m())
                        cs[pi][k].divNonzeros(l2s_[pi][k].vs());
                }
            }

            bases_[i]->synthesise(cs[pi], cs[i], !bases_[pi]->isTransient());
        }
        else
        {
            bases_[0]->synthesise(f, cs[0], false);
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
    return l1l2s_;
}
