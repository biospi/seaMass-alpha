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


#include "BasisLibrary.hpp"
#include "Bspline.hpp"
#include "../io/FileNetcdf.hpp"
#include <limits>
#include <iomanip>
#include <cmath>
#include <sstream>
using namespace std;
using namespace kernel;


BasisLibrary::BasisLibrary(std::vector<Basis*>& bases, const BasisGrid::GridInfo& parentGridInfo,
                           const std::string& dbFilename, bool transient) :
    BasisGrid(bases, parentGridInfo, transient), gTs_(1), gs_(1)
{
    ostringstream oss3;
    oss3 << "Library filename=" << dbFilename;
    type() = oss3.str();

    if (getDebugLevel() % 10 >= 2)
    {
        ostringstream oss;
        oss << getTimeStamp() << "   " << getIndex() << " " << oss3.str() << " parent = " << getParentIndex() << " ...";
        info(oss.str());
    }

    if (getDebugLevel() % 10 >= 1)
    {
        ostringstream oss;
        oss << getTimeStamp() << "     Loading " << dbFilename << " ...";
        info(oss.str());
    }

    {
        FileNetcdf fileIn(dbFilename);
        ostringstream oss2;
        oss2 << "m" << setfill('0') << setw(2) << ii(parentGridInfo.colScale[0]);

        MatrixSparse db;
        fileIn.readMatrixSparseCsr(db, oss2.str());

        //  read in offset
        ii db_offset = fileIn.readAttribute<ii>("offset", "", fileIn.openGroup(oss2.str()));

        // extract relevant submatrix (not efficient atm and ought to be moved to SparseMatrix)
        vector<ii> is1;
        vector<ii> js1;
        vector<fp> vs1;
        ii i = -1;
        ii j_min = parentGridInfo.colOffset[0];
        ii j_max = parentGridInfo.colOffset[0] + parentGridInfo.colExtent[0] - 1;
        bool new_row;
        for (ii i0 = 0; i0 < db.m(); ++i0)
        {
            new_row = true;
            for (ii k = db.ijs()[i0]; k < db.ijs()[i0 + 1]; ++k)
            {
                ii j = db_offset + db.js()[k];

                if (j >= j_min && j <= j_max) {
                    if (new_row)
                    {
                        i++;
                        new_row = false;
                        ids_.push_back(i0);
                    }

                    is1.push_back(i);
                    js1.push_back(j - j_min);
                    vs1.push_back(db.vs()[k]);
                }
            }
        }
        ii m = i + 1;
        ii n = parentGridInfo.colExtent[0];

        gridInfo().colOffset[0] = -1;
        gridInfo().colExtent[0] = m;

        aT_.importFromCoo(m, n, is1.size(), is1.data(), js1.data(), vs1.data());
    }
    a_.transpose(aT_);
 
    if (getDebugLevel() % 10 >= 2)
    {
        ostringstream oss3;
        oss3 << getTimeStamp() << "     " << gridInfo();
        info(oss3.str());
    }
}


BasisLibrary::~BasisLibrary()
{
}


void
BasisLibrary::
synthesize(vector<MatrixSparse> &f, const vector<MatrixSparse> &x, bool accumulate)
{
    if (getDebugLevel() % 10 >= 3)
    {
        ostringstream oss;
        oss << getTimeStamp() << "     " << getIndex() << " BasisGridPeak::synthesise";
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


void BasisLibrary::analyze(vector<MatrixSparse> &xE, const vector<MatrixSparse> &fE, bool sqrA)
{
    if (getDebugLevel() % 10 >= 3)
    {
        ostringstream oss;
        oss << getTimeStamp() << "     " << getIndex() << " BasisGridPeak::analyse";
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



