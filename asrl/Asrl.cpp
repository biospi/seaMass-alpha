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


#include "Asrl.hpp"
#include "BasisMatrix.hpp"
#include "OptimizerSrl.hpp"
#include "OptimizerAccelerationEve1.hpp"
#include <iomanip>
using namespace std;
using namespace kernel;


void Asrl::notice()
{
    cout << "seaMass-ASRL : Copyright (C) 2016 - biospi Laboratory, University of Bristol, UK" << endl;
    cout << "This program comes with ABSOLUTELY NO WARRANTY." << endl;
    cout << "This is free software, and you are welcome to redistribute it under certain conditions." << endl;
}


Asrl::Asrl(Input &input, fp lambda, fp lambdaGroup, bool taperShrinkage, fp tolerance) : bT_(input.bT), lambda_(lambda), lambdaGroup_(lambdaGroup), lambdaGroupStart_(lambdaGroup), taperShrinkage_(taperShrinkage), tolerance_(tolerance), iteration_(0)
{
    if (getDebugLevel() % 10 >= 1)
    {
        cout << getTimeStamp() << "  Initialising basis functions ..." << endl;
    }

    new BasisMatrix(bases_, input.aT, input.gT.size() > 0 ? &input.gT : 0, false);

    innerOptimizer_ = new OptimizerSrl(bases_, bT_);
    optimizer_ = new OptimizerAccelerationEve1(innerOptimizer_);
    optimizer_->setLambda(fp(lambda_), fp(lambdaGroup_));
}


Asrl::~Asrl()
{
    delete optimizer_;
    delete innerOptimizer_;

    for (ii i = 0; i < ii(bases_.size()); i++)
        delete bases_[i];
}


bool Asrl::step()
{
    if (iteration_ == 0 && getDebugLevel() % 10 >= 1)
    {
        li nnz = 0;
        li nx = 0;
        for (ii j = 0; j < (ii)bases_.size(); j++)
        {
            if (!static_cast<Basis *>(bases_[j])->isTransient())
            {
                for (size_t k = 0; k < optimizer_->xs()[j].size(); k++)
                {
                    nnz += optimizer_->xs()[j][k].nnz();
                    nx += optimizer_->xs()[j][k].size();
                }
            }
        }

        cout << getTimeStamp();
        cout << "   it:     0 nx: " << setw(10) << nx << " nnz: " << setw(10) << nnz;
        cout << " tol:  " << fixed << setprecision(8) << setw(10) << tolerance_ << endl;
    }

    iteration_++;
    double grad = optimizer_->step();

    if (getDebugLevel() % 10 >= 1)
    {
        li nnz = 0;
        for (ii j = 0; j < (ii)bases_.size(); j++)
        {
            if (!static_cast<Basis *>(bases_[j])->isTransient())
            {
                for (size_t k = 0; k < optimizer_->xs()[j].size(); k++)
                {
                    nnz += optimizer_->xs()[j][k].nnz();
                }
            }
        }
        cout << getTimeStamp();
        cout << "   it: " << setw(5) << iteration_;
        cout << " shrink: ";
        cout.unsetf(ios::floatfield);
        cout << setprecision(4) << setw(6) << lambda_;
        cout << " nnz: " << setw(10) << nnz;
        cout << " grad: " << fixed << setprecision(8) << setw(10) << grad << endl;
    }

    if (grad <= tolerance_)
    {
        if (lambda_ == 0.0 || !taperShrinkage_)
        {
            if (getDebugLevel() % 10 == 0) cout << "o" << endl;

            return false;
        }
        else
        {
            if (getDebugLevel() % 10 == 0) cout << "o" << flush;

            lambda_ *= 0.5;
            lambdaGroup_ *= 0.5;

            if (lambda_ < 0.0625 && lambdaGroup_ < 0.0625)
            {
                lambda_ = 0.0;
                lambdaGroup_ = 0.0;
            }

            optimizer_->setLambda(fp(lambda_), fp(lambdaGroup_));
        }
    }
    else
    {
        if (getDebugLevel() % 10 == 0) cout << "." << flush;
    }
    
    if (grad != grad)
    {
        throw runtime_error("BUG: grad turned into NAN.");
    }
    
    return true;
}


ii Asrl::getIteration() const
{
    return iteration_;
}


void Asrl::getOutput(Output& output) const
{
    if (getDebugLevel() % 10 >= 1)
         cout << getTimeStamp() << "  Getting output ..." << endl;

    output.xT.resize(optimizer_->xs()[0].size());
    for (ii i = 0; i < ii(output.xT.size()); i++)
        output.xT[i].copy(optimizer_->xs()[0][i]);

    vector<MatrixSparse> f;
    {
        vector<vector<MatrixSparse> > cs;
        optimizer_->synthesise(f, cs);
    }
    output.xTaT.resize(f.size());
    for (ii i = 0; i < ii(output.xTaT.size()); i++)
    {
        output.xTaT[i].alloc(f[i].m(), f[i].n());
        f[i].exportTo(output.xTaT[i].vs());
    }

    if (lambdaGroupStart_ > 0.0)
        bases_[0]->synthesiseGroups(output.xTgT, optimizer_->xs()[0], false);
}
