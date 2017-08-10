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


BasisBsplineMz::BasisBsplineMz(std::vector<Basis*>& bases, vector<MatrixSparse>& b, const std::vector<fp>& binCounts,
                               const std::vector<li>& binCountsIndex_, const std::vector<double>& binEdges,
                               char scale, short chargeStates, bool transient)
        : BasisBspline(bases, chargeStates <= 1 ? 1 : 2, transient), bGridInfo_(1),
          chargeDeconvolution_(chargeStates > 0), chargeStates_(chargeStates == 0 ? 1 : chargeStates)
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

    // blur the input data by the B-spline kernel and place in 'b'
    Bspline bspline(3, 65536); // bspline basis function lookup table
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
        double mzMin = numeric_limits<double>::max();
        double mzMax = 0.0;
        double xDiff = 0.0;
        li n = 0;
        for (ii k = 0; k < ii(binEdgesIndex.size()) - 1; k++)
        {
            mzMin = binEdges[binEdgesIndex[k]] < mzMin ? binEdges[binEdgesIndex[k]] : mzMin;
            mzMax = binEdges[binEdgesIndex[k + 1] - 1] > mzMax ? binEdges[binEdgesIndex[k + 1] - 1] : mzMax;

            // find mean difference in index between edges, ignoring first, last and zeros
            for (ii i = 1; i < binCountsIndex[k + 1] - binCountsIndex[k] - 1; i++)
            {
                if (binCounts[binCountsIndex[k] + i] != 0)
                {
                    ii ei = binEdgesIndex[k] + i;
                    xDiff += log2(binEdges[ei + 1] - PROTON_MASS) - log2(binEdges[ei] - PROTON_MASS);
                    n++;
                }
            }
        }
        xDiff /= double(n);

        ii scaleAuto = ii(ceil(log2(1.0 / xDiff))) + 1;
        if (scale == numeric_limits<char>::max())
        {
            scale = scaleAuto;

            if (getDebugLevel() % 10 >= 1)
            {
                ostringstream oss;
                oss << getTimeStamp() << "     autodetected_mz_scale=" << fixed << setprecision(1) << int(scale);
                info(oss.str());
            }
        }
        if (getDebugLevel() % 10 >= 2)
        {
            ostringstream oss;
            oss << getTimeStamp() << "     range=" << fixed << setprecision(3) << mzMin << ":" << mzMax << "Th";
            info(oss.str());
            ostringstream oss2;
            oss2 << getTimeStamp() << "     charge_deconvolution=" << chargeDeconvolution_;
            info(oss2.str());
            ostringstream oss3;
            oss3 << getTimeStamp() << "     charge_states=" << chargeStates_;
            info(oss3.str());
        }

        bGridInfo_.count = ii(binCountsIndex.size()) - 1;
        bGridInfo_.scale[0] = scale;
        bGridInfo_.offset[0] = ii(floor(log2(mzMin - PROTON_MASS) * (1L << scale))) - 1;
        bGridInfo_.extent[0] = (ii(ceil(log2(mzMax - PROTON_MASS) * (1L << scale))) + 1) - bGridInfo_.offset[0] + 1;

        vector<MatrixSparse> bs(bGridInfo_.count);
        for (ii k = 0; k < bGridInfo_.count; k++)
        {
            vector<ii> rowind;
            vector<ii> colind;
            vector<fp> acoo;

            // create transformation matrix
            for (ii i = 0; i < binEdgesIndex[k + 1] - binEdgesIndex[k] - 1; i++)
            {
                ii startNz = ii(acoo.size());
                double rowSum = 0.0;

                if (binCounts[binCountsIndex[k] + i] >= 0.0)
                {
                    ii ei = binEdgesIndex[k] + i;
                    double xfMin = log2(binEdges[ei] - PROTON_MASS) * (1L << scale);
                    double xfMax = log2(binEdges[ei + 1] - PROTON_MASS) * (1L << scale);

                    ii xMin = ii(floor(xfMin)) - 1;
                    ii xMax = ii(ceil(xfMax)) + 1;

                    // work out basis coefficients
                    for (ii x = xMin; x <= xMax; x++)
                    {
                        double bfMin = double(x - 2);
                        double bfMax = double(x + 2);

                        // intersection of bin and basis, between 0.0 and 4.0
                        double bMin = xfMin > bfMin ? xfMin - bfMin : 0.0;
                        double bMax = xfMax < bfMax ? xfMax - bfMin : bfMax - bfMin;

                        // basis coefficient b is _integral_ of area under b-spline basis
                        fp bc = fp(bspline.ibasis(bMax) - bspline.ibasis(bMin));
                        if (bc > 0.0)
                        {
                            rowSum += bc;
                            acoo.push_back(bc);
                            rowind.push_back(i);
                            colind.push_back(x - bGridInfo_.offset[0]);
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
            a.importFromCoo(binCountsIndex[k + 1] - binCountsIndex[k], bGridInfo_.n(), acoo.size(),
                            rowind.data(), colind.data(), acoo.data());

            Matrix t;
            t.importFromArray(1, binCountsIndex[k + 1] - binCountsIndex[k], &binCounts.data()[binCountsIndex[k]]);
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

    // create our kernel
    vector<fp> hs(5);
    gridInfo().count = bGridInfo_.count;
    gridInfo().scale[0] = bGridInfo_.scale[0];
    gridInfo().offset[0] = bGridInfo_.offset[0] - ii(hs.size() - 1) / 2;
    gridInfo().extent[0] = (bGridInfo_.extent[0] + ii(hs.size() - 1)) + ii(log2(fp(chargeStates_)) * (1L << scale));
    if (gridInfo().dimensions > 1)
    {
        gridInfo().scale[1] = 0;
        gridInfo().offset[1] = 0;
        gridInfo().extent[1] = chargeStates_;
    }

    if (getDebugLevel() % 10 >= 2)
    {
        ostringstream oss;
        oss << getTimeStamp() << "     input=B" << b[0];
        info(oss.str());

        /*for (ii k = 0; k < ii(hs.size()); k++)
        {
            ostringstream oss2;
            oss2 << getTimeStamp() << "     kernel=" << hs[k];
            info(oss2.str());
        }*/

        ostringstream oss3;
        oss3 << getTimeStamp() << "     " << gridInfo();
        info(oss3.str());
    }

    // isotopes
    double mz1 = pow(2.0, (gridInfo().offset[0] + gridInfo().extent[0]) / double(1L << scale)) + PROTON_MASS;
    ii maxCarbons = ii(ceil(mz1 * 0.008));
    cout << "maxCarbons=" << maxCarbons;
    vector< vector<double> > carbons(maxCarbons);
    carbons[0].resize(2);
    carbons[0][0] = 0.9893;
    carbons[0][1] = 0.0107;

    for (ii i = 1; i < ii(carbons.size()); i++)
        convolution(carbons[i], carbons[i - 1], carbons[0]);

    ii maxIsotopes = 0;
    for (ii i = 1; i < ii(carbons.size()); i++)
    {
        double sum = 0.0;
        for (ii j = 0; j < ii(carbons[i].size()); j++)
        {
            sum += carbons[i][j];

            if (sum > 0.99999)
            {
                carbons[i].resize(j + 1);
                maxIsotopes = maxIsotopes > carbons[i].size() ? maxIsotopes : carbons[i].size();
                break;
            }
        }
    }

    // create A as a temporary COO matrix
    ii m = bGridInfo_.extent[0];
    ii n = gridInfo().extent[0];

    aTs_.resize(chargeStates_);
    as_.resize(chargeStates_);
    for (short z = 0; z < chargeStates_; z++)
    {
        aTs_[z].init(n, m);
        fp iz = log2(fp(z + 1)) * (1L << scale);

        //for (ii c = 0; c < maxIsotopes; c++)
        for (ii c = 0; c < 1; c++)
        {
            cout << z << ":" << (c+1) << "/" << maxIsotopes << endl;

            vector<ii> is;
            vector<ii> js;
            vector<fp> vs;
            for (ii j = 0; j < n; j++)
            {
                double iMono = j - iz;
                double mz0 = pow(2.0, (gridInfo().offset[0] + iMono + 0.5) / double(1L << scale)) + PROTON_MASS;

                double mass0 = pow(2.0, (gridInfo().offset[0] + iMono + 0.5) / double(1L << scale) + log2(z+1));

                ii nCarbons = ii(ceil(mass0 * 0.008));
                if (c < carbons[nCarbons - 1].size())
                {
                    //cout << z << ":" << mass0 << ":" << nCarbons << ":" << carbons[nCarbons - 1].size() << endl;

                    double mz = mz0 + c * 1.0033548378 / (z + 1);
                    double iIsotope = log2(mz - PROTON_MASS) * (1L << scale) - bGridInfo_.offset[0] - 0.5;

                    double t;
                    double hc = modf(iIsotope, &t);
                    ii i0 = ii(iIsotope);
                    if (hc >= 0.5)
                    {
                        i0++;
                        hc -= 1.0;
                    }

                    hs[0] = bspline.ibasis(0.5 - hc) - bspline.ibasis(0.0);
                    hs[1] = bspline.ibasis(1.5 - hc) - bspline.ibasis(0.5 - hc);
                    hs[2] = bspline.ibasis(2.5 - hc) - bspline.ibasis(1.5 - hc);
                    hs[3] = bspline.ibasis(3.5 - hc) - bspline.ibasis(2.5 - hc);
                    hs[4] = bspline.ibasis(4.0) - bspline.ibasis(3.5 - hc);

                    for (ii k = 0; k < ii(hs.size()); k++)
                    {
                        ii i = i0 + k - ii(hs.size() - 1) / 2;

                        if (i >= 0 && i < m && hs[k] > 0.0)
                        {
                            is.push_back(i);
                            js.push_back(j);
                            //vs.push_back(hs[k] * carbons[nCarbons - 1][c]);
                            vs.push_back(hs[k]);
                        }
                    }
                }
            }

            // add to A
            MatrixSparse a;
            a.importFromCoo(n, m, vs.size(), js.data(), is.data(), vs.data());

            MatrixSparse t;
            t.add(1.0, false, a, aTs_[z]);
            aTs_[z].swap(t);
        }

        //ostringstream oss; oss << z << ".coo";
        //FileNetcdf fileOut(oss.str(), NC_NETCDF4);
        //fileOut.write(aTs_[z], "At");

        as_[z].transpose(aTs_[z]);
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

    for (short c = 0; c < chargeStates_; c++)
    {
        MatrixSparseView row(x[0], c);

        // zero basis functions that are no longer needed
        /*MatrixSparse t;
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
        }*/

        // synthesise
        f[0].matmul(false, row, aTs_[c], (accumulate || c > 0) ? true : false);
    }

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

    vector<MatrixSparse> xEs(chargeStates_);
    for (short c = 0; c < chargeStates_; c++)
    {
        if (sqrA)
        {
            MatrixSparse t;
            t.sqr(as_[c]);

            xEs[c].matmul(false, fE[0], t, false);
        }
        else
        {
            xEs[c].matmul(false, fE[0], as_[c], false);
        }
    }

    xE.resize(1);
    xE[0].concatenateRows(xEs);

    if (getDebugLevel() % 10 >= 3)
    {
        ostringstream oss;
        oss << getTimeStamp() << "       " << xE[0];
        info(oss.str());
    }
}


const BasisBsplineMz::GridInfo& BasisBsplineMz::getBGridInfo() const
{
    return bGridInfo_;
}