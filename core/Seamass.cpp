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
#include "BasisBsplineCharge.hpp"
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


Seamass::Seamass(Input& input, const string& isotopesFilename, const std::vector<short>& scale,
                 fp lambda, fp lambdaGroup, bool taperShrinkage, fp tolerance, double peakFwhm, short chargeStates) :
        innerOptimizer_(0), scale_(scale), isotopesFilename_(isotopesFilename), lambda_(lambda),
        lambdaStart_(lambda), lambdaGroup_(lambdaGroup), lambdaGroupStart_(lambdaGroup),
        taperShrinkage_(taperShrinkage), tolerance_(tolerance),
        iteration_(0), peakFwhm_(peakFwhm), chargeStates_(chargeStates)
{
    init(input, true);
    optimizer_->setLambda(lambda_, lambdaGroup_);
}


Seamass::Seamass(Input& input, const Output& output) :
        innerOptimizer_(0), scale_(output.scale), isotopesFilename_(output.isotopesFilename), lambda_(output.lambda),
        lambdaStart_(output.lambda), lambdaGroup_(output.lambdaGroup), lambdaGroupStart_(output.lambdaGroup),
        tolerance_(output.tolerance), iteration_(0), peakFwhm_(output.peakFwhm), chargeStates_(output.chargeStates)
{
    init(input, false);

    // import seed
    optimizer_->xs().resize(output.xs.size());
    for (ii k = 0; k < ii(optimizer_->xs().size()); k++)
    {
        if (!bases_[k]->isTransient())
        {
            optimizer_->xs()[k].resize(1);
            optimizer_->xs()[k][0].copy(output.xs[k]);
        }
   }

    optimizer_->l2s().resize(output.l2s.size());
    for (ii k = 0; k < ii(optimizer_->l2s().size()); k++)
    {
        if (!bases_[k]->isTransient())
        {
            optimizer_->l2s()[k].resize(1);
            optimizer_->l2s()[k][0].copy(output.l2s[k]);
        }
    }

    optimizer_->l1l2s().resize(output.l1l2s.size());
    for (ii k = 0; k < ii(optimizer_->l1l2s().size()); k++)
    {
        if (!bases_[k]->isTransient())
        {
            optimizer_->l1l2s()[k].resize(1);
            optimizer_->l1l2s()[k][0].copy(output.l1l2s[k]);
        }
    }

    optimizer_->setLambda(lambda_, lambdaGroup_);
}


Seamass::~Seamass()
{
    delete optimizer_;
    if (innerOptimizer_)
        delete innerOptimizer_;

    for (ii i = 0; i < (ii)bases_.size(); i++)
        delete bases_[i];
}


void Seamass::init(Input& input, bool seed)
{
    // INIT BASIS FUNCTIONS
    // Create our tree of bases
    if (getDebugLevel() % 10 >= 1)
    {
        ostringstream oss;
        oss << getTimeStamp() << "  Initialising overcomplete tree of basis functions ...";
        info(oss.str());
    }

    /*if (input.countsIndex.size() <= 2)
    {
        dimensions_ = 1;

        //if (peakFwhm_ > 0.0)
        //    new BasisBsplinePeak(bases_, bases_.back()->getIndex(), peakFwhm_, true);

        outputBasis_ = new BasisBsplineMz(bases_, b_, isotopesFilename_, input.counts, input.countsIndex,
                                          input.locations, scale_[0], chargeStates_, false);

        for (ii i = 0; static_cast<BasisBspline*>(bases_.back())->getGridInfo().colScale[1] > 8; i++)
            new BasisBsplineScale(bases_,  bases_.back()->getIndex(), 1, 1, true, false);
    }
    else
    {
        dimensions_ = 2;

        //if (peakFwhm_ > 0.0)
        //    new BasisBsplinePeak(bases_, bases_.back()->getIndex(), peakFwhm_, true);

        new BasisBsplineMz(bases_, b_, isotopesFilename_, input.counts, input.countsIndex, input.locations,
                           scale_[0], chargeStates_, true);

        outputBasis_ = new BasisBsplineScantime(bases_, bases_.back()->getIndex(), input.startTimes,
                                                input.finishTimes, input.exposures, scale_[1], false);

        Basis* previousBasis = outputBasis_;
        for (ii i = 0; static_cast<BasisBspline*>(bases_.back())->getGridInfo().colScale[1] > 8; i++)
        {
            if (i > 0)
                previousBasis = new BasisBsplineScale(bases_, previousBasis->getIndex(), 1, 1, true, false);

            while (static_cast<BasisBspline*>(bases_.back())->getGridInfo().rowExtent[0] > 4)
                new BasisBsplineScale(bases_, bases_.back()->getIndex(), 0, 0, true, false);
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
            optimizer_->setLambda(lambda_, lambdaGroup_);
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


void Seamass::getOutput(Output& output, bool synthesize) const
{
    if (getDebugLevel() % 10 >= 1)
    {
        ostringstream oss;
        oss << getTimeStamp() << "  Getting " << (synthesize ? "synthesized " : "") << "output ...";
        info(oss.str());
    }

    output = Output();

    output.scale = scale_;
    output.lambda = lambdaStart_;
    output.lambdaGroup = lambdaGroupStart_;
    output.tolerance = tolerance_;
    output.peakFwhm = peakFwhm_;
    output.chargeStates = chargeStates_;
    output.isotopesFilename = isotopesFilename_;

    output.gridInfos.resize(bases_.size());
    for (ii k = 0; k < ii(bases_.size()); k++)
        output.gridInfos[k] = static_cast<BasisBspline*>(bases_[k])->getGridInfo();

    if (synthesize)
    {
        vector<vector<MatrixSparse> > xs;
        {
            vector<MatrixSparse> f;
            optimizer_->synthesize(f, xs);
        }

        output.xs.resize(xs.size());
        for (ii k = 0; k < ii(xs.size()); k++)
            output.xs[k].copy(xs[k][0]);

        output.l2s.resize(optimizer_->l2s().size());
        for (ii k = 0; k < ii(optimizer_->l2s().size()); k++)
                output.l2s[k].copy(optimizer_->l2s()[k][0]);

        output.l1l2s.resize(optimizer_->l1l2s().size());
        for (ii k = 0; k < ii(optimizer_->l1l2s().size()); k++)
        {
            output.l1l2s[k].copy(optimizer_->l1l2s()[k][0]);

            if (!bases_[k]->isTransient())
                output.l1l2s[k].addNonzeros(-optimizer_->getLambda());
        }
    }
    else
    {
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
        {
            if (!bases_[k]->isTransient())
            {
                output.l1l2s[k].copy(optimizer_->l1l2s()[k][0]);
                output.l1l2s[k].addNonzeros(-optimizer_->getLambda());
            }
        }
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


void Seamass::getInput(Input &input, bool reconstruct) const
{
    if (getDebugLevel() % 10 >= 1)
    {
        ostringstream oss;
        oss << getTimeStamp() << (reconstruct ? "  Deriving reconstructed output ..." : "  Deriving blurred input ...");
        info(oss.str());
    }

    /*const BasisBspline::GridInfo& meshInfo = static_cast<BasisBsplineMz*>(bases_[0])->getBGridInfo();
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
        if (input.countsIndex.size() > 0)
            input.countsIndex[i + 1] = (i + 1) * meshInfo.n();

        for (ii j = 0; j <= meshInfo.n(); j++)
        {
            double mz = pow(2.0, (meshInfo.colOffset[0] + j) / double(1L << meshInfo.colScale[0])) +
                        BasisBsplineMz::PROTON_MASS;
            input.locations[i * (meshInfo.n() + 1) + j] = mz;
        }
    }*/
}


void Seamass::getOutputControlPoints(ControlPoints& controlPoints) const
{
    if (getDebugLevel() % 10 >= 1)
    {
        ostringstream oss;
        oss << getTimeStamp() << "  Deriving control points ...";
        info(oss.str());
    }

    const BasisBspline::GridInfo& meshInfo =
            static_cast<BasisBspline*>(bases_[dimensions_])->getGridInfo();

    vector<MatrixSparse> c(1);
    {
        vector<vector<MatrixSparse> > cs;
        optimizer_->synthesize(c, cs, dimensions_);
    }

    vector<fp>(meshInfo.size()).swap(controlPoints.coeffs);
    c[0].exportToDense(controlPoints.coeffs.data());

    controlPoints.scale = meshInfo.colScale;
    controlPoints.offset = meshInfo.colOffset;
    controlPoints.extent = meshInfo.colExtent;
}


void Seamass::getOutputControlPoints1d(ControlPoints& controlPoints, bool density) const
{
    // todo: move this code to Basis class

    if (getDebugLevel() % 10 >= 1)
    {
        ostringstream oss;
        oss << getTimeStamp() << "  Deriving 1D control points ...";
        info(oss.str());
    }

    const BasisBspline::GridInfo& meshInfo = static_cast<BasisBspline*>(bases_[1])->getGridInfo();
    controlPoints.scale.resize(1);
    controlPoints.scale[0] = meshInfo.colScale[1];

    controlPoints.offset.resize(1);
    controlPoints.offset[0] = meshInfo.colOffset[1];

    controlPoints.extent.resize(2);
    controlPoints.extent[0] = meshInfo.colExtent[1]
                              + ii(log2(double(meshInfo.colExtent[0])) * double(1L << controlPoints.scale[0]));
    controlPoints.extent[1] = meshInfo.rowExtent[0];

    vector<MatrixSparse> c(1);
    {
        vector<vector<MatrixSparse> > cs;
        optimizer_->synthesize(c, cs, 1);
    }

    // sum across charge states
    vector<fp>(controlPoints.extent[0] * controlPoints.extent[1], 0.0f).swap(controlPoints.coeffs);
    {
        vector<fp> t(meshInfo.size());
        c[0].exportToDense(t.data());

        for (ii i = 0; i < meshInfo.rowExtent[0]; i++)
        {
            for (ii j = 0; j < meshInfo.colExtent[1]; j++)
            {
                for (ii z = 0; z < meshInfo.colExtent[0]; z++)
                {
                    controlPoints.coeffs[i * controlPoints.extent[0] +
                            ii(log2(double(z+1)) * double(1L << controlPoints.scale[0])) + j] +=
                            t[i * meshInfo.colExtent[0] * meshInfo.colExtent[1] + z * meshInfo.colExtent[1] + j];
                }
            }
        }
    }

    // temporary hack (integrate instead)
    if (density)
    {
        for (ii x = 0; x < controlPoints.extent[0]; x++)
        {
            double mz0 = pow(2.0, (controlPoints.offset[0] + x - 0.5) / double(1L << controlPoints.scale[0])) +
                         BasisBsplineMz::PROTON_MASS;
            double mz1 = pow(2.0, (controlPoints.offset[0] + x + 0.5) / double(1L << controlPoints.scale[0])) +
                         BasisBsplineMz::PROTON_MASS;

            for (ii y = 0; y < controlPoints.extent[1]; y++)
                controlPoints.coeffs[x + y * controlPoints.extent[0]] /= mz1 - mz0;
        }
    }
}
