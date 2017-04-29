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
        {   // initialise starting estimate of 'x' from analysis of 'b'
            if (getDebugLevel() % 10 >= 1)
                cout << getTimeStamp() << "  Seeding from analysis of input ..." << endl;

            vector<MatrixSparse> t(b_.size());
            for (ii k = 0; k < ii(t.size()); k++)
                t[k].copy(b_[k]);

            analyze(xs_, t, false, false);
        }

        {   // compute L2 norm of each basis function and store in 'l2s'
            if (getDebugLevel() % 10 >= 1)
                cout << getTimeStamp() << "  Initialising L2 norms ..." << endl;

            vector<MatrixSparse> t(b_.size());
            for (ii k = 0; k < ii(t.size()); k++)
                t[k].copy(1, b_[k].n(), (fp) 1.0);

            analyze(l2s_, t, true, false);
        }

        {   // compute L1 norm of each L2 normalised basis function and store in 'l1l2s'
            if (getDebugLevel() % 10 >= 1)
                cout << getTimeStamp() << "  Initialising L1 norms of L2 norms ..." << endl;

            vector<MatrixSparse> t(b_.size());
            for (size_t k = 0; k < t.size(); k++)
                t[k].copy(1, b_[k].n(), (fp) 1.0);

            l1l2s_.resize(bases_.size());

            analyze(l1l2s_, t, false);
        }

        // initialise starting estimate of 'x' from analysis of 'b'
        if (getDebugLevel() % 10 >= 1)
            cout << getTimeStamp() << "  Seeding from analysis of input 2 ..." << endl;

        double sumB = 0.0;
        double sumX = 0.0;
        xs_.resize(bases_.size());
        for (ii i = 0; i < (ii)bases_.size(); i++)
        {
            if (i == 0)
            {
                vector<MatrixSparse> t(b_.size());
                for (size_t k = 0; k < t.size(); k++)
                {
                    t[k].copy(b_[k]);
                    sumB += t[k].sum();
                }

                bases_[i]->analyze(xs_[0], t, false);
            }
            else
            {
                bases_[i]->analyze(xs_[i], xs_[bases_[i]->getParentIndex()], false);
            }

            for (size_t k = 0; k < xs_[i].size(); k++)
                sumX += xs_[i][k].sum();
        }

        for (ii i = 0; i < (ii)bases_.size(); i++)
        {
            if (!bases_[i]->isTransient())
            {
                for (size_t k = 0; k < xs_[i].size(); k++)
                {
                   // remove unneeded l1l2sPlusLambda
                    MatrixSparse l1l2PlusLambda;
                    l1l2PlusLambda.copy(xs_[i][k]);
                    l1l2PlusLambda.copySubset(l1l2s_[i][k]);

                    // normalise and prune xs
                    MatrixSparse x;
                    x.copy(xs_[i][k]);
                    x.divNonzeros(l1l2PlusLambda);
                    x.mul((fp)(sumB / sumX));
                    xs_[i][k].copyPrune(x, pruneThreshold);

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
        synthesize(f_fE, xEs);
    }
    double synthesisDuration = getElapsedTime() - synthesisStart;

    // ERROR
    if (getDebugLevel() % 10 >= 3)
        cout << getTimeStamp() << "    Error ..." << endl;

    double errorStart = getElapsedTime();
    {
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
        analyze(xEs, f_fE, false);
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
                    vector<MatrixSparse>* g = bases_[i]->getGroups(false);
                    vector<MatrixSparse>* gT = bases_[i]->getGroups(true);

                    // group and individual shrinkage
                    for (size_t k = 0; k < ys[i].size(); k++)
                    {

                        // y = groupNorm(x)
                        MatrixSparse t;
                        t.copy(xs_[i][k]);
                        t.sqr();
                        ys[i][k].matmul(false, t, (*gT)[k], false);
                        t.matmul(false, ys[i][k], (*g)[k], false);
                        ys[i][k].copy(xs_[i][k]);
                        ys[i][k].copySubset(t);
                        ys[i][k].sqrt();

                        // y = x * groupNorm(x)^-1)
                        ys[i][k].div2Nonzeros(xs_[i][k]);

                        // y = lambdaGroup * x * groupNorm(x)^-1
                        ys[i][k].mul(lambdaGroup_);

                        // y = lambda + lambdaGroup * x * groupNorm(x)^-1
                        ys[i][k].addNonzeros(lambda_);

                        // y = l1l2 + lambda + lambdaGroup * x * groupNorm(x)^-1
                        ys[i][k].addNonzeros(l1l2s_[i][k]);

                        // y = x / (l1l2 + lambda + lambdaGroup * x * groupNorm(x)^-1)
                        ys[i][k].div2Nonzeros(xs_[i][k]);

                        // y = xE * x / (l1l2 + lambda + lambdaGroup * x * groupNorm(x)^-1)
                        ys[i][k].mul(xEs[i][k]);
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
                        ys[i][k].div2Nonzeros(xs_[i][k]);

                        // y = xE * x / (l1l2 + lambda)
                        ys[i][k].mul(xEs[i][k]);
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
                    sumSqrDiffs += xs_[i][k].sumSqrDiffsNonzeros(ys[i][k]);
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
                    xs_[i][k].copyPrune(ys[i][k], pruneThreshold_);
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
        cout << getTimeStamp()  << "                                     total_sort = " << setprecision(4) << setw(12) << MatrixSparse::sortElapsed_ << endl;
     }

    return sqrt(sumSqrDiffs) / sqrt(sumSqrs);
}


void OptimizerSrl::synthesize(vector<MatrixSparse>& f, vector< vector<MatrixSparse> >& xEs, ii basis)
{
    if (xEs.size() != bases_.size())
        xEs.resize(bases_.size());

    for (ii i = ii(bases_.size()) - 1; i >= 0; i--)
    {
        if (!xEs[i].size() && !bases_[i]->isTransient())
        {
            xEs[i].resize(xs_[i].size());
            for (size_t k = 0; k < xEs[i].size(); k++)
            {
                xEs[i][k].copy(xs_[i][k]);
                
                if (l2s_.size() == xs_.size() && l2s_[i][k].m())
                    xEs[i][k].divNonzeros(l2s_[i][k]);
            }
        }

        if (basis == i) // return with B-spline control points
        {
            for (size_t k = 0; k < xEs[i].size(); k++)
                f[k].copy(xEs[i][k]);

            break;
        }

        if (i > 0)
        {
            ii pi = bases_[i]->getParentIndex();
            
            if (!xEs[pi].size() && !bases_[pi]->isTransient())
            {
                xEs[pi].resize(xs_[pi].size());
                for (size_t k = 0; k < xEs[pi].size(); k++)
                {
                    xEs[pi][k].copy(xs_[pi][k]);
                    
                    if (l2s_.size() == xs_.size() && l2s_[pi][k].m())
                        xEs[pi][k].divNonzeros(l2s_[pi][k]);
                }
            }

            bases_[i]->synthesize(xEs[pi], xEs[i], !bases_[pi]->isTransient());
        }
        else
        {
            bases_[0]->synthesize(f, xEs[0], false);
        }
    }
}


void OptimizerSrl::analyze(std::vector<std::vector<MatrixSparse> > &xEs, std::vector<MatrixSparse> &fE, bool l2, bool l2Normalize)
{
    if (xEs.size() != bases_.size())
        xEs.resize(bases_.size());

    vector<MatrixSparse> t;
    bases_.front()->analyze(t, fE, l2);

    if (xEs[0].size() != t.size())
        xEs[0].resize(t.size());

    for (ii k = 0; k < ii(xEs[0].size()); k++)
        //xEs[0][k].copySubset(t[k]);
        xEs[0][k].copy(t[k]);

    for (ii l = 1; l < ii(bases_.size()); l++)
    {
        vector<MatrixSparse> t;
        bases_[l]->analyze(t, xEs[bases_[l]->getParentIndex()], l2);

        if (xEs[l].size() != t.size())
            xEs[l].resize(t.size());

        for (ii k = 0; k < ii(xEs[l].size()); k++)
            //xEs[l][k].copySubset(t[k]);
            xEs[l][k].copy(t[k]);
    }

    for (ii l = 0; l < ii(bases_.size()); l++)
    {
        if (!bases_[l]->isTransient())
        {
            for (ii k = 0; k < ii(xEs[l].size()); k++)
            {
                if (xEs[l][k].nnz() > xs_[l][k].nnz())
                {
                    MatrixSparse t;
                    t.copy(xs_[l][k]);
                    t.copySubset(xEs[l][k]);
                    xEs[l][k].copy(t);
                }
            }

            if(l2)
            {
                for (ii k = 0; k < ii(xEs[l].size()); k++)
                    l2s_[l][k].sqrt();
            }

            if (l2Normalize)
            {
                for (ii k = 0; k < ii(xEs[l].size()); k++)
                    xEs[l][k].divNonzeros(l2s_[l][k]);
            }
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
