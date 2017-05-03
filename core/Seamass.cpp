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


#include "Seamass.hpp"

#include "../asrl/OptimizerAccelerationEve1.hpp"

#include "BasisBsplineMz.hpp"
#include "BasisBsplineScale.hpp"
#include "BasisBsplineScantime.hpp"

#include <cstring>
#include <iomanip>


using namespace std;
using namespace kernel;


void Seamass::notice()
{
    cout << "seaMass - Copyright (C) 2016 - biospi Laboratory, University of Bristol, UK" << endl;
    cout << "This program comes with ABSOLUTELY NO WARRANTY." << endl;
    cout << "This is free software, and you are welcome to redistribute it under certain conditions." << endl;
}


Seamass::Seamass(const Input& input, const std::vector<char>& scale, fp lambda, fp tolerance) : lambda_(lambda), lambdaStart_(lambda), tolerance_(tolerance), iteration_(0)
{
    init(input, scale, true);
}


Seamass::Seamass(const Input& input, const Output& seed) : lambda_(seed.shrinkage), lambdaStart_(seed.shrinkage), tolerance_(seed.tolerance), iteration_(0)
{
    init(input, seed.scale, false);

    // import seed
    for (ii k = 0; k < (ii)bases_.size(); k++)
    {
        if (!bases_[k]->isTransient())
        {
            optimizer_->xs()[k].resize(1);
            optimizer_->xs()[k][0].copy(seed.xs[k]);

            optimizer_->l2s()[k].resize(1);
            optimizer_->l2s()[k][0].copy(seed.l2s[k]);

            optimizer_->l1l2s()[k].resize(1);
            optimizer_->l1l2s()[k][0].copy(seed.l1l2s[k]);
        }
   }
}


Seamass::~Seamass()
{
    delete optimizer_;
    delete innerOptimizer_;

    for (ii i = 0; i < (ii)bases_.size(); i++)
        delete bases_[i];
}


void Seamass::init(const Input& input, const std::vector<char>& scales, bool seed)
{
    // for speed only, merge bins if rc_mz is set more than 8 times higher than the bin width
    // this is conservative, 4 times might be ok, but 2 times isn't enough
    //merge_bins(mzs, intensities, 0.125 / rc_mz);
    
    // INIT BASIS FUNCTIONS
    // Create our tree of bases
    if (getDebugLevel() % 10 >= 1)
    {
        cout << getTimeStamp() << "  Initialising overcomplete tree of basis functions ..." << endl;
    }
    if (input.countsIndex.size() <= 2)
    {
        dimensions_ = 1;

        new BasisBsplineMz(bases_, input.counts, input.countsIndex, input.locations, scales[0], false);
        
        while (static_cast<BasisBspline*>(bases_.back())->getGridInfo().scale[0] > -6)
        {
            new BasisBsplineScale(bases_, bases_.back()->getIndex(), 0, false);
        }
    }
    else
    {
        dimensions_ = 2;

        new BasisBsplineMz(bases_, input.counts, input.countsIndex, input.locations, scales[0], true);
        Basis* previousBasis = new BasisBsplineScantime(bases_, bases_.back()->getIndex(), input.startTimes, input.finishTimes, input.exposures, scales[1], false);
 
        for (ii i = 0; static_cast<BasisBspline*>(bases_.back())->getGridInfo().scale[0] > -6; i++)
        {
            if (i > 0)
            {
                previousBasis = new BasisBsplineScale(bases_, previousBasis->getIndex(), 0, false);
            }
            
            while (static_cast<BasisBspline*>(bases_.back())->getGridInfo().extent[1] > 4)
            {
                new BasisBsplineScale(bases_, bases_.back()->getIndex(), 1, false);
            }
        }
    }

    // INIT B
    if (input.countsIndex.size() == 0)
    {
        b_.resize(1);
        b_[0].copy(1, input.counts.size(), input.counts.data());
    }
    else
    {
        b_.resize(input.countsIndex.size() - 1);
        for (ii i = 0; i < input.countsIndex.size() - 1; i++)
        {
            b_[i].copy(1, input.countsIndex[i + 1] - input.countsIndex[i], &input.counts.data()[input.countsIndex[i]]);
        }
    }

    // INIT OPTIMISER
    innerOptimizer_ = new OptimizerSrl(bases_, b_, seed);
    optimizer_ = new OptimizerAccelerationEve1(innerOptimizer_);
    optimizer_->setLambda((fp) lambda_);
}


bool Seamass::step()
{
    if (iteration_ == 0 && getDebugLevel() % 10 >= 1)
    {
        li nnz = 0;
        li nx = 0;
        for (ii j = 0; j < (ii)bases_.size(); j++)
        {
            if (!static_cast<BasisBspline *>(bases_[j])->isTransient())
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
            if (!static_cast<BasisBspline *>(bases_[j])->isTransient())
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
        if (lambda_ == 0)
        {
            if (getDebugLevel() % 10 == 0) cout << "o" << endl;
            return false;
        }
        else
        {
            if (getDebugLevel() % 10 == 0) cout << "o" << flush;
            lambda_ *= (lambda_ > 0.0625 ? 0.5 : 0.0);
            optimizer_->setLambda((fp) lambda_);
        }
    }
    else
    {
        if (getDebugLevel() % 10 == 0) cout << "." << flush;
    }
    
    if (grad != grad)
    {
        cout << "ARGGH!" << endl;
        return false;
    }
    
    return true;
}


ii Seamass::getIteration() const
{
    return iteration_;
}


void Seamass::getOutput(Output& output) const
{
    if (getDebugLevel() % 10 >= 1)
    {
        cout << getTimeStamp() << "  Getting output ..." << endl;
    }

    output = Output();

    const BasisBspline::GridInfo& meshInfo = static_cast<BasisBspline*>(bases_[dimensions_ - 1])->getGridInfo();
    output.scale = meshInfo.scale;

    output.shrinkage = lambdaStart_;
    output.tolerance = tolerance_;

    output.xs.resize(bases_.size());
    output.l2s.resize(bases_.size());
    output.l1l2s.resize(bases_.size());
    for (ii k = 0; k < (ii)bases_.size(); k++)
    {
        output.xs[k].copy(optimizer_->xs()[k][0]);
        output.l2s[k].copy(optimizer_->l2s()[k][0]);
        output.l1l2s[k].copy(optimizer_->l1l2s()[k][0]);
    }

    /*output.baselineScale.resize(dimensions_);
    output.baselineOffset.resize(dimensions_);
    output.baselineExtent.resize(dimensions_);
    for (ii i = 0; i < dimensions_; i++)
    {
        const BasisBspline::GridInfo& meshInfo = static_cast<BasisBspline*>(bases_[dimensions_ - 1])->getGridInfo();

        output.baselineScale[i] = meshInfo.scale[i];
        output.baselineOffset[i] = meshInfo.offset[i];
        output.baselineExtent[i] = meshInfo.extent[i];
    }

    output.shrinkage = lambdaStart_;
    output.tolerance = tolerance_;

    output.scales.resize(dimensions_);
    output.offsets.resize(dimensions_);
    for (ii k = 0; k < (ii)bases_.size(); k++)
    {
        if (bases_[k]->isTransient() == Basis::Transient::NO)
        {
            const BasisBspline::GridInfo& meshInfo = static_cast<BasisBspline*>(bases_[k])->getGridInfo();

            vector<fp> vs;
            vector<ii> is;
            vector<ii> js;
            optimizer_->xs()[k][0].exportTo(is, js, vs);

            for (ii nz = 0; nz < vs.size(); nz++)
            {
                cout << "scale=" << (int) meshInfo.scale[0] << " offset=" << meshInfo.offset[0] + js[nz] << " value=" << vs[nz] << endl;

                output.weights.push_back(vs[nz]);
                output.scales[0].push_back(meshInfo.scale[0]);
                output.offsets[0].push_back(meshInfo.offset[0] + js[nz]);
                //cout << acoo[nz] << ":" << rowind[nz] << "," << colind[nz] << endl;
            }
        }
    }*/
}


void Seamass::getOutputBinCounts(std::vector<fp>& binCounts) const
{
    if (getDebugLevel() % 10 >= 1)
        cout << getTimeStamp() << "  Deriving restored bin counts ..." << endl;

    vector<MatrixSparse> f;
    optimizer_->synthesise(f);

    li i = 0;
    for (ii k = 0; k < ii(f.size()); k++)
    {
        f[k].exportTo(&binCounts.data()[i]);
        i += f[k].size();
    }
}


void Seamass::getOutputControlPoints(ControlPoints& controlPoints) const
{
    if (getDebugLevel() % 10 >= 1)
        cout << getTimeStamp() << "  Deriving control points ..." << endl;

    const BasisBspline::GridInfo& meshInfo = static_cast<BasisBspline*>(bases_[dimensions_ - 1])->getGridInfo();

    vector<MatrixSparse> c(1);
    optimizer_->synthesise(c, dimensions_ - 1);
    vector<fp>(meshInfo.size()).swap(controlPoints.coeffs);
    c[0].exportTo(controlPoints.coeffs.data());

    controlPoints.scale = meshInfo.scale;
    controlPoints.offset = meshInfo.offset;
    controlPoints.extent = meshInfo.extent;
}


void Seamass::getOutputControlPoints1d(ControlPoints& controlPoints) const
{
    if (getDebugLevel() % 10 >= 1)
        cout << getTimeStamp() << "  Deriving 1D control points ..." << endl;

    const BasisBspline::GridInfo& meshInfo = static_cast<BasisBspline*>(bases_[0])->getGridInfo();

    vector<MatrixSparse> c(1);
    optimizer_->synthesise(c, 0);
    vector<fp>(meshInfo.size()).swap(controlPoints.coeffs);
    c[0].exportTo(controlPoints.coeffs.data());

    controlPoints.scale = meshInfo.scale;
    controlPoints.offset = meshInfo.offset;
    controlPoints.extent = meshInfo.extent;
    controlPoints.extent.push_back(meshInfo.count);
}
