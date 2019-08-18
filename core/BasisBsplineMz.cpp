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


BasisBsplineMz::BasisBsplineMz(std::vector<Basis*>& bases, vector<MatrixSparse>& b, bool transient) :
        BasisBspline(bases, 1, 2, transient), bGridInfo_(1, 1), gTs_(1), gs_(1)
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

    bGridInfo_.rowScale[0] = 0;
    bGridInfo_.rowOffset[0] = 0;
    bGridInfo_.rowExtent[0] = b[0].m();
    bGridInfo_.colScale[0] = 0;
    bGridInfo_.colOffset[0] = 0;
    bGridInfo_.colExtent[0] = b[0].n();

    if (getDebugLevel() % 10 >= 2)
    {
        ostringstream oss4;
        oss4 << getTimeStamp() << "     b_" << bGridInfo_;
        info(oss4.str());
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

    /*if (getDebugLevel() % 10 >= 2)
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
    }*/

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Gt = m x n matrix where m are the coefficients and n are the groups (monoisotope centroid mass).
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    /*{
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

                //double mass = pow(2.0, (gridInfo().colOffset[1] + g) / double(1L << gridInfo().colScale[1]));
                //cout << mass << endl;

                is.push_back(x + z * gridInfo().colExtent[1]);
                js.push_back(g);
                vs.push_back(1.0);
                //vs.push_back(1.0 / sqrt(mass); // this does not work
                //vs.push_back(1.0 / pow(300.0*mass, 1.0/4.0));
                //vs.push_back(1.0 / pow(6.0*mass, 1.0/3.0)); //vs.push_back(1.0 / sqrt(pow(6.0*mass, 2.0/3.0)));
            }
        }

        //for (ii nz = 0; nz < ii(vs.size()); nz++)
        //{
        //    vs[nz] /= sqrt(fp(gSizes[js[nz]]));
        //}

        gTs_[0].importFromCoo(m, n, vs.size(), is.data(), js.data(), vs.data());
        gs_[0].transpose(gTs_[0]);
    }*/
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
