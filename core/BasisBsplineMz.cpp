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


#include "BasisBsplineMz.hpp"
#include "Bspline.hpp"
#include <limits>
#include <iomanip>
#include <cmath>
#include <sstream>
#include "../io/FileNetcdf.hpp"
using namespace std;
using namespace kernel;


double BasisBsplineMz::PROTON_MASS = 1.007276466879;


void convolution(vector<double>& x, const vector<double>& a, const vector<double>& b)
{
    x.resize(a.size() + b.size() - 1, 0);
    for (ii i = 0; i < ii(a.size()); i++)
    {
        for (ii j = 0; j < ii(b.size()); j++)
        {
            x[i + j] += a[i] * b[j];
        }
    }
}


BasisBsplineMz::BasisBsplineMz(std::vector<Basis*>& bases, vector<MatrixSparse>& b, const string& isotopesFilename,
                               const std::vector<fp>& binCounts, const std::vector<li>& binCountsIndex_,
                               const std::vector<double>& binEdges, short scale, short chargeStates, bool transient) :
        BasisBspline(bases, 1, 2, transient), bGridInfo_(1, 1), chargeDeconvolution_(chargeStates > 0), gTs_(1), gs_(1)
{
    if (getDebugLevel() % 10 >= 2)
    {
        ostringstream oss;
        oss << getTimeStamp();
        if (getDebugLevel() % 10 >= 2)
            oss << "   " << getIndex() << " BasisBsplineMz";
        else
            oss << "   BasisBsplineMz";
        if (isTransient()) oss << " (transient)";
        oss << " ...";
        info(oss.str());
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // blur the input data by the B-spline kernel and place in 'b'
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // bspline basis function lookup table
    Bspline bspline(3, 65536);
    {
        std::vector<li> binCountsIndex;
        std::vector<li> binEdgesIndex;
        if (binCountsIndex_.size() > 0)
        {
            binCountsIndex = binCountsIndex_;
            binEdgesIndex = binCountsIndex_;
            for (ii i = 0; i < ii(binEdgesIndex.size()); i++) binEdgesIndex[i] += i;
        }
        else
        {
            binCountsIndex.push_back(0);
            binCountsIndex.push_back(binCounts.size());
            binEdgesIndex.push_back(0);
            binEdgesIndex.push_back(binEdges.size());
        }

        // find min and max m/z across spectra, m for each A
        double mz0 = numeric_limits<double>::max();
        double mz1 = 0.0;
        double xDiff = 0.0;
        li n = 0;
        for (ii k = 0; k < ii(binEdgesIndex.size()) - 1; k++)
        {
            mz0 = binEdges[binEdgesIndex[k]] < mz0 ? binEdges[binEdgesIndex[k]] : mz0;
            mz1 = binEdges[binEdgesIndex[k + 1] - 1] > mz1 ? binEdges[binEdgesIndex[k + 1] - 1] : mz1;

            // find mean difference in index between edges, ignoring first, last and zeros
            for (ii i = 1; i < binCountsIndex[k + 1] - binCountsIndex[k] - 1; i++)
            {
                if (binCounts[binCountsIndex[k] + i] != 0)
                {
                    li ei = binEdgesIndex[k] + i;
                    xDiff += log2(binEdges[ei + 1] - PROTON_MASS) - log2(binEdges[ei] - PROTON_MASS);
                    n++;
                }
            }
        }
        xDiff /= double(n);

        ii scaleAuto = ii(ceil(log2(1.0 / xDiff))) + 1;
        if (scale == numeric_limits<short>::max())
        {
            scale = scaleAuto;

            if (getDebugLevel() % 10 >= 1)
            {
                ostringstream oss;
                oss << getTimeStamp() << "     autodetected_mz_scale=" << fixed << setprecision(1) << int(scale);
                info(oss.str());
            }
        }

        double scale2 = pow(2.0, scale);
        bGridInfo_.rowScale[0] = numeric_limits<short>::min();
        bGridInfo_.rowOffset[0] = 0;
        bGridInfo_.rowExtent[0] = ii(binCountsIndex.size()) - 1;
        bGridInfo_.colScale[0] = scale;
        bGridInfo_.colOffset[0] = ii(floor(log2(mz0 - PROTON_MASS) * scale2));
        bGridInfo_.colExtent[0] = (ii(ceil(log2(mz1 - PROTON_MASS) * scale2))) - bGridInfo_.colOffset[0] + 1;

        if (getDebugLevel() % 10 >= 2)
        {
            ostringstream oss;
            oss << getTimeStamp() << "     range=" << fixed << setprecision(3) << mz0 << ":" << mz1 << "Th";
            info(oss.str());
            ostringstream oss2;
            oss2 << getTimeStamp() << "     charge_deconvolution=" << chargeDeconvolution_;
            info(oss2.str());
            ostringstream oss3;
            oss3 << getTimeStamp() << "     charge_states=" << chargeStates;
            info(oss3.str());
            ostringstream oss4;
            oss4 << getTimeStamp() << "     b_" << bGridInfo_;
            info(oss4.str());
        }

        vector<MatrixSparse> bs(bGridInfo_.rowExtent[0]);
        for (ii k = 0; k < bGridInfo_.rowExtent[0]; k++)
        {
            vector<ii> rowind;
            vector<ii> colind;
            vector<fp> acoo;

            // create transformation matrix
            for (ii i = 0; i < binEdgesIndex[k + 1] - binEdgesIndex[k] - 1; i++)
            {
                auto startNz = ii(acoo.size());
                double rowSum = 0.0;

                if (binCounts[binCountsIndex[k] + i] >= 0.0)
                {
                    li ei = binEdgesIndex[k] + i;
                    double xfMin = log2(binEdges[ei] - PROTON_MASS) * scale2;
                    double xfMax = log2(binEdges[ei + 1] - PROTON_MASS) * scale2;

                    auto xMin = ii(floor(xfMin)) - 2;
                    auto xMax = ii(ceil(xfMax)) + 2;

                    // work out basis coefficients
                    for (ii x = xMin; x <= xMax; x++)
                    {
                        double bfMin = x - 1.5;
                        double bfMax = x + 2.5;

                        // intersection of bin and basis, between 0.0 and 4.0
                        double bMin = xfMin > bfMin ? xfMin - bfMin : 0.0;
                        double bMax = xfMax < bfMax ? xfMax - bfMin : bfMax - bfMin;

                        // basis coefficient b is _integral_ of area under b-spline basis
                        auto bc = fp(bspline.ibasis(bMax) - bspline.ibasis(bMin));

                        ii j = x - bGridInfo_.colOffset[0];
                        if (j >= 0 && j < bGridInfo_.colExtent[0] && bc > 0.0)
                        {
                            rowSum += bc;
                            acoo.push_back(bc);
                            rowind.push_back(i);
                            colind.push_back(j);
                       }
                    }
                }

                // normalise column
                if (rowSum > 0.0)
                {
                    for (ii nz = startNz; nz < ii(acoo.size()); nz++)
                        acoo[nz] /= rowSum;
                }
            }

            // create b
            MatrixSparse a;
            a.importFromCoo(ii(binCountsIndex[k + 1] - binCountsIndex[k]), bGridInfo_.n(), acoo.size(),
                            rowind.data(), colind.data(), acoo.data());

            Matrix t;
            t.importFromArray(1, ii(binCountsIndex[k + 1] - binCountsIndex[k]), &binCounts.data()[binCountsIndex[k]]);
            MatrixSparse t2, t3;
            t2.importFromMatrix(t);
            t3.matmul(false, t2, a, false);
            bs[k].pruneCells(t3);
        }

        b.resize(1);
        b[0].concatenateRows(bs);

        if (scaleAuto != scale && getDebugLevel() % 10 >= 2)
        {
            ostringstream oss;
            oss << "WARNING: mz_scale is not the suggested value of " << scaleAuto << ". Continue at your own risk!";
            warning(oss.str());
        }
    }

    if (getDebugLevel() % 10 >= 2)
    {
        ostringstream oss;
        oss << getTimeStamp() << "     input=B" << b[0];
        info(oss.str());
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Load 'A'
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    if (getDebugLevel() % 10 >= 2)
    {
        ostringstream oss;
        oss << getTimeStamp() << "    Loading " << isotopesFilename << "...";
        info(oss.str());
    }

    FileNetcdf fileIn(isotopesFilename);
    ostringstream oss2;
    oss2 << "s=" << setfill('0') << setw(2) << ii(bGridInfo_.colScale[0]);
    int groupId = fileIn.openGroup(oss2.str());

    vector<ii> offset;
    fileIn.readAttribute(offset, "offset", "", groupId);

    {
        vector<short> zs;
        vector<ii> is;
        vector<ii> js;
        vector<fp> vs;

        ii iMin = numeric_limits<ii>::max();
        ii iMax = 0;
        ii jMin = bGridInfo_.colOffset[0];
        ii jMax = bGridInfo_.colOffset[0] + bGridInfo_.colExtent[0] - 1;

        for (short z = 0; z < chargeStates; z++)
        {
            if (getDebugLevel() % 10 >= 2)
            {
                ostringstream oss;
                oss << getTimeStamp() << "     z" << (z + 1);
                info(oss.str());
            }

            ostringstream oss;
            oss << "z=" << setfill('0') << setw(4) << (z + 1);
            MatrixSparse aTz;
            fileIn.readMatrixSparseCsr(aTz, oss.str(), groupId);

            // this should be in MatrixSparse
            for (ii _i = 0; _i < aTz.m(); _i++)
            {
                for (ii nz = aTz.ijs()[_i]; nz < aTz.ijs()[_i + 1]; nz++)
                {
                    ii i = offset[0] + _i;
                    ii j = offset[1] + aTz.js()[nz];

                    if (jMin <= j && j <= jMax)
                    {
                        iMin = iMin < i ? iMin : i;
                        iMax = iMax > i ? iMax : i;

                        zs.push_back(z);
                        is.push_back(i);
                        js.push_back(j - jMin);
                        vs.push_back(aTz.vs()[nz]);
                    }
                }
            }
        }

        for (ii nz = 0; nz < ii(is.size()); nz++)
        {
            is[nz] = (is[nz] - iMin) + zs[nz] * (iMax - iMin + 1);

            //cout << iNs[nz] << "," << jNs[nz] << "=" << vNs[nz] << endl;
        }

        ii mN = chargeStates * (iMax - iMin + 1);
        ii nN = jMax - jMin + 1;

        aT_.importFromCoo(mN, nN, vs.size(), is.data(), js.data(), vs.data());
        a_.transpose(aT_);

        // Set up 'A'

        gridInfo().rowScale[0] = bGridInfo_.rowScale[0];
        gridInfo().rowOffset[0] = bGridInfo_.rowOffset[0];
        gridInfo().rowExtent[0] = bGridInfo_.rowExtent[0];

        gridInfo().colScale[0] = numeric_limits<short>::min();
        gridInfo().colOffset[0] = 0;
        gridInfo().colExtent[0] = chargeStates > 0 ? chargeStates : 1;

        gridInfo().colScale[1] = bGridInfo_.colScale[0];
        gridInfo().colOffset[1] = iMin;
        gridInfo().colExtent[1] = iMax - iMin + 1;
    }

    if (getDebugLevel() % 10 >= 2)
    {
        ostringstream oss;
        oss << getTimeStamp() << "     a_" << gridInfo();
        info(oss.str());
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Gt = m x n matrix where m are the coefficients and n are the groups (monoisotope centroid mass).
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    {
        vector<ii> is;
        vector<ii> js;
        vector<fp> vs;

        ii m = getGridInfo().colExtent[0] * getGridInfo().colExtent[1];
        ii n = getGridInfo().colExtent[1] +  ii(round(log2(double(chargeStates)) * (1L << getGridInfo().colScale[1])));

        vector<ii> gSizes(n, 0);
        for (ii z = 0; z < gridInfo().colExtent[0]; z++)
        {
            auto g0 = ii(round(log2(double(z + 1)) * (1L << getGridInfo().colScale[1])));

            for (ii x = 0; x < gridInfo().colExtent[1]; x++)
            {
                ii g = g0 + x;
                gSizes[g]++;

                double mass = pow(2.0, (gridInfo().colOffset[1] + g) / double(1L << gridInfo().colScale[1]));
                //cout << mass << endl;

                is.push_back(x + z * gridInfo().colExtent[1]);
                js.push_back(g);
                vs.push_back(1.0);
                //vs.push_back(1.0 / sqrt(sqrt(mass)));
            }
        }

        /*for (ii nz = 0; nz < ii(vs.size()); nz++)
        {
            vs[nz] /= sqrt(fp(gSizes[js[nz]]));
        }*/

        gTs_[0].importFromCoo(m, n, vs.size(), is.data(), js.data(), vs.data());
        gs_[0].transpose(gTs_[0]);
    }
}


BasisBsplineMz::~BasisBsplineMz()
{
}


void
BasisBsplineMz::
synthesize(vector<MatrixSparse> &f, const vector<MatrixSparse> &x, bool accumulate)
{
    if (getDebugLevel() % 10 >= 3)
    {
        ostringstream oss;
        oss << getTimeStamp() << "     " << getIndex() << " BasisBsplineMz::synthesize";
        info(oss.str());
    }

    if (!f.size())
        f.resize(1);

    // zero basis functions that are no longer needed
    MatrixSparse t;
    ii rowsPruned = t.pruneRows(aT_, x[0], false, 0.75);
    if (rowsPruned > 0)
    {
        aT_.swap(t);
        a_.transpose(aT_);

        if (getDebugLevel() % 10 >= 3)
        {
            ostringstream oss;
            oss << getTimeStamp() << "      " << getIndex() << " pruned " << rowsPruned << " basis functions";
            info(oss.str());
        }
    }

    // synthesise
    f[0].matmul(false, x[0], aT_, accumulate);

    if (getDebugLevel() % 10 >= 3)
    {
        ostringstream oss;
        oss << getTimeStamp() << "       " << f[0];
        info(oss.str());
    }
}


void BasisBsplineMz::analyze(vector<MatrixSparse> &xE, const vector<MatrixSparse> &fE, bool sqrA)
{
    if (getDebugLevel() % 10 >= 3)
    {
        ostringstream oss;
        oss << getTimeStamp() << "     " << getIndex() << " BasisBsplineMz::analyse";
        info(oss.str());
    }

    if (!xE.size())
        xE.resize(1);

    if (sqrA)
    {
        MatrixSparse t;
        t.sqr(a_);
        xE[0].matmul(false, fE[0], t, false);
    }
    else
    {
        xE[0].matmul(false, fE[0], a_, false);
    }

    if (getDebugLevel() % 10 >= 3)
    {
        ostringstream oss;
        oss << getTimeStamp() << "       " << xE[0];
        info(oss.str());
    }
}


const vector<MatrixSparse> * BasisBsplineMz::getColGroups(bool transpose) const
{
    if (transpose)
        return &gTs_;
    else
        return &gs_;
}


const BasisBsplineMz::GridInfo& BasisBsplineMz::getBGridInfo() const
{
    return bGridInfo_;
}