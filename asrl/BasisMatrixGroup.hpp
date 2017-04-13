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

#ifndef SEAMASS_ASRL_BASISMATRIXGROUP_HPP
#define SEAMASS_ASRL_BASISMATRIXGROUP_HPP


#include "BasisMatrix.hpp"


class BasisMatrixGroup : public BasisMatrix
{
public:
    BasisMatrixGroup(std::vector<Basis*>& bases, ii aM, ii aN, std::vector<fp>& aV, std::vector<ii>& aI, std::vector<ii>& aJ,
                     ii gM, std::vector<fp>& gV, std::vector<ii>& gI, std::vector<ii>& gJ, Transient transient);
    virtual ~BasisMatrixGroup();

    void groupSynthesis(std::vector<MatrixSparse>& f, const std::vector<MatrixSparse>& x, bool accumulate) const;

    void shrinkage(std::vector<MatrixSparse>& y, const std::vector<MatrixSparse>& x, const std::vector<MatrixSparse>& xE, const std::vector<MatrixSparse>& l1l2, fp lambda) const;


private:
    MatrixSparse gT_;
    MatrixSparse g_;
};


#endif //SEAMASS_BASISMATRIXGROUP_HPP
