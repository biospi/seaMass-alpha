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
#include "DatasetSeamass.hpp"
#include "BasisBsplineMz.hpp"
#include "BasisBsplinePeak.hpp"
#include "BasisBsplineScale.hpp"
#include "BasisBsplineScantime.hpp"
#include "../asrl/OptimizerAccelerationEve1.hpp"
#include <kernel.hpp>
#include <cstring>
#include <cmath>
#include <iomanip>
#include <sstream>
using namespace std;
using namespace kernel;


void Seamass::notice()
{
    cout << "seaMass - Copyright (C) 2016 - biospi Laboratory, University of Bristol, UK" << endl;
    cout << "This program comes with ABSOLUTELY NO WARRANTY." << endl;
    cout << "This is free software, and you are welcome to redistribute it under certain conditions." << endl;
}


Seamass::Seamass(Input& input, const std::vector<char>& scale, fp lambda, bool taperShrinkage, fp tolerance,
                 double peakFwhm) : innerOptimizer_(0), lambda_(lambda), lambdaStart_(lambda), taperShrinkage_(taperShrinkage), tolerance_(tolerance), peakFwhm_(peakFwhm), iteration_(0)
{
    init(input, scale, true);
    optimizer_->setLambda(fp(lambda_));
}


Seamass::Seamass(Input& input, const Output& seed) : innerOptimizer_(0), lambda_(seed.shrinkage), lambdaStart_(seed.shrinkage), tolerance_(seed.tolerance), peakFwhm_(seed.peakFwhm), iteration_(0)
{
    init(input, seed.scale, false);

    // import seed
    optimizer_->xs().resize(seed.xs.size());
    for (ii k = 0; k < ii(optimizer_->xs().size()); k++)
    {
        if (!bases_[k]->isTransient())
        {
            optimizer_->xs()[k].resize(1);
            optimizer_->xs()[k][0].copy(seed.xs[k]);
        }
   }

    optimizer_->l2s().resize(seed.l2s.size());
    for (ii k = 0; k < ii(optimizer_->l2s().size()); k++)
    {
        if (!bases_[k]->isTransient())
        {
            optimizer_->l2s()[k].resize(1);
            optimizer_->l2s()[k][0].copy(seed.l2s[k]);
        }
    }

    optimizer_->l1l2s().resize(seed.l1l2s.size());
    for (ii k = 0; k < ii(optimizer_->l1l2s().size()); k++)
    {
        if (!bases_[k]->isTransient())
        {
            optimizer_->l1l2s()[k].resize(1);
            optimizer_->l1l2s()[k][0].copy(seed.l1l2s[k]);
        }
    }

    optimizer_->setLambda((fp) lambda_);
}


Seamass::~Seamass()
{
    delete optimizer_;
    if (innerOptimizer_)
        delete innerOptimizer_;

    for (ii i = 0; i < (ii)bases_.size(); i++)
        delete bases_[i];
}


void Seamass::init(Input& input, const std::vector<char>& scales, bool seed)
{
    // INIT BASIS FUNCTIONS
    // Create our tree of bases
    if (getDebugLevel() % 10 >= 1)
    {
        ostringstream oss;
        oss << getTimeStamp() << "  Initialising overcomplete tree of basis functions ...";
        info(oss.str());
    }
    if (input.countsIndex.size() <= 2)
    {
        dimensions_ = 1;

        new BasisBsplineMz(bases_, b_, input.counts, input.countsIndex, input.locations,scales[0], true);

        if (peakFwhm_ > 0.0)
            new BasisBsplinePeak(bases_, bases_.back()->getIndex(), peakFwhm_, false);

        while (static_cast<BasisBspline*>(bases_.back())->getGridInfo().scale[0] > 10)
            new BasisBsplineScale(bases_, bases_.back()->getIndex(), 0, false);
    }
    else
    {
        dimensions_ = 2;

        new BasisBsplineMz(bases_, b_, input.counts, input.countsIndex, input.locations, scales[0], true);

        if (peakFwhm_ > 0.0)
            new BasisBsplinePeak(bases_, bases_.back()->getIndex(), peakFwhm_, true);

        Basis* previousBasis = new BasisBsplineScantime(bases_, bases_.back()->getIndex(), input.startTimes,
                                                        input.finishTimes, input.exposures, scales[1], false);

        for (ii i = 0; static_cast<BasisBspline*>(bases_.back())->getGridInfo().scale[0] > 10; i++)
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
    /*if (input.countsIndex.size() == 0)
    {
        Matrix b;
        b.importFromArray(1, input.counts.size(), input.counts.data());
        b_[0].importFromMatrix(b);
    }
    else
    {
        b_.resize(input.countsIndex.size() - 1);
        for (ii i = 0; i < input.countsIndex.size() - 1; i++)
        {
            Matrix b;
            b.importFromArray(1, input.countsIndex[i + 1] - input.countsIndex[i],
                                  &input.counts.data()[input.countsIndex[i]]);
            b_[i].importFromMatrix(b);
        }
    }*/

    // INIT OPTIMISER
    innerOptimizer_ = new OptimizerSrl(bases_, b_, seed);
    optimizer_ = new OptimizerAccelerationEve1(innerOptimizer_);
    //optimizer_ = new OptimizerSrl(bases_, b_, seed);
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

        ostringstream oss;
        oss << getTimeStamp();
        oss << "   it:     0 nx: " << setw(10) << nx << " nnz: " << setw(10) << nnz;
        oss << " tol:  " << fixed << setprecision(8) << setw(10) << tolerance_;
        info(oss.str());
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
                     nnz += optimizer_->xs()[j][k].nnz();
             }
        }

        ostringstream oss;
        oss << getTimeStamp();
        oss << "   it: " << setw(5) << iteration_;
        oss << " shrink: ";
        oss.unsetf(ios::floatfield);
        oss << setprecision(4) << setw(6) << lambda_;
        oss << " nnz: " << setw(10) << nnz;
        oss << " grad: " << fixed << setprecision(8) << setw(10) << grad;
        info(oss.str());
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
        ostringstream oss;
        oss << "BUG: convergence failed!";
        error(oss.str());
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
        ostringstream oss;
        oss << getTimeStamp() << "  Getting output ...";
        info(oss.str());
    }

    output = Output();

    const BasisBspline::GridInfo& meshInfo =
            static_cast<BasisBspline*>(bases_[dimensions_ - (peakFwhm_ > 0.0 ? 0 : 1)])->getGridInfo();
    output.scale = meshInfo.scale;

    output.shrinkage = lambdaStart_;
    output.tolerance = tolerance_;
    output.peakFwhm = peakFwhm_;

    output.xs.resize(optimizer_->xs().size());
    for (ii k = 0; k < ii(optimizer_->xs().size()); k++)
        if (!bases_[k]->isTransient())
            output.xs[k].copy(optimizer_->xs()[k][0]);

    output.l2s.resize(optimizer_->l2s().size());
    for (ii k = 0; k < ii(optimizer_->l2s().size()); k++)
        if (!bases_[k]->isTransient())
            output.l2s[k].copy(optimizer_->l2s()[k][0]);

    output.l1l2s.resize(optimizer_->l1l2s().size());
    for (ii k = 0; k < ii(optimizer_->l1l2s().size()); k++)
        if (!bases_[k]->isTransient())
            output.l1l2s[k].copy(optimizer_->l1l2s()[k][0]);

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


void Seamass::getInput(Input &input, bool reconstruct) const
{
    if (getDebugLevel() % 10 >= 1)
    {
        ostringstream oss;
        oss << getTimeStamp() << (reconstruct ? "  Deriving reconstructed output ..." : "  Deriving blurred input ...");
        info(oss.str());
    }

    const BasisBspline::GridInfo& meshInfo = static_cast<BasisBsplineMz*>(bases_[0])->getBGridInfo();
    vector<fp>(meshInfo.size()).swap(input.counts);

    if (reconstruct)
    {
        vector<MatrixSparse> f;
        {
            vector<vector<MatrixSparse> > cs;
            optimizer_->synthesize(f, cs);
        }

        f[0].exportToDense(input.counts.data());
    }
    else
    {
        b_[0].exportToDense(input.counts.data());
    }

    vector<double>(meshInfo.size() + meshInfo.m()).swap(input.locations);
    for (ii i = 0; i < meshInfo.m(); i++)
    {
        input.countsIndex[i + 1] = (i + 1) * meshInfo.n();

        for (ii j = 0; j <= meshInfo.n(); j++)
        {
            double mz = pow(2.0, (meshInfo.offset[0] + j - 0.5) / double(1L << meshInfo.scale[0])) +
                        BasisBsplineMz::PROTON_MASS;
            input.locations[i * (meshInfo.n() + 1) + j] = mz;
        }
    }
}


void Seamass::getOutputControlPoints(ControlPoints& controlPoints, bool deconvolve) const
{
    if (getDebugLevel() % 10 >= 1)
    {
        ostringstream oss;
        oss << getTimeStamp() << "  Deriving control points ...";
        info(oss.str());
    }

    const BasisBspline::GridInfo& meshInfo =
            static_cast<BasisBspline*>(bases_[dimensions_ - (peakFwhm_ > 0.0 && deconvolve ? 0 : 1)])->getGridInfo();

    vector<MatrixSparse> c(1);
    {
        vector<vector<MatrixSparse> > cs;
        optimizer_->synthesize(c, cs, dimensions_ - (peakFwhm_ > 0.0 && deconvolve ? 0 : 1));
    }

    vector<fp>(meshInfo.size()).swap(controlPoints.coeffs);
    c[0].exportToDense(controlPoints.coeffs.data());

    controlPoints.scale = meshInfo.scale;
    controlPoints.offset = meshInfo.offset;
    controlPoints.extent = meshInfo.extent;
}


void Seamass::getOutputControlPoints1d(ControlPoints& controlPoints, bool deconvolve) const
{
    if (getDebugLevel() % 10 >= 1)
    {
        ostringstream oss;
        oss << getTimeStamp() << "  Deriving 1D control points ...";
        info(oss.str());
    }

    const BasisBspline::GridInfo& meshInfo = static_cast<BasisBspline*>(bases_[peakFwhm_ > 0.0 && deconvolve ? 1 : 0])->getGridInfo();

    vector<MatrixSparse> c(1);
    {
        vector<vector<MatrixSparse> > cs;
        optimizer_->synthesize(c, cs, peakFwhm_ > 0.0 && deconvolve ? 1 : 0);
    }

    // temporary hack (differentiate instead)
    for (ii nz = 0; nz < c[0].nnz(); nz++)
    {
        double mz0 = pow(2.0, (meshInfo.offset[0] + c[0].js_[nz] - 0.5) / double(1L << meshInfo.scale[0])) +
                BasisBsplineMz::PROTON_MASS;
        double mz1 = pow(2.0, (meshInfo.offset[0] + c[0].js_[nz] + 0.5) / double(1L << meshInfo.scale[0])) +
                BasisBsplineMz::PROTON_MASS;

        c[0].vs_[nz] /= mz1 - mz0;
    }

    vector<fp>(meshInfo.size()).swap(controlPoints.coeffs);
    c[0].exportToDense(controlPoints.coeffs.data());

    controlPoints.scale = meshInfo.scale;
    controlPoints.offset = meshInfo.offset;
    controlPoints.extent = meshInfo.extent;
    controlPoints.extent.push_back(meshInfo.count);
}
