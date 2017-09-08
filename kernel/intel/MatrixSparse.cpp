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


#include "MatrixSparse.hpp"
#include "kernel.hpp"
#include <iomanip>
#include <sstream>
#include <cmath>
#include <cassert>
#include <cstring>
#include <algorithm>
#include <ippcore.h>
#include <ipps.h>
#if defined(_OPENMP)
  #include <omp.h>
#endif
using namespace std;
using namespace kernel;


MatrixSparse::MatrixSparse() : m_(0), n_(0), mat_(0), ijs1_(0), isCsrOwned_(false), isSorted_(true)
{
}


MatrixSparse::~MatrixSparse()
{
    empty();
}


void MatrixSparse::init(ii m, ii n)
{
    empty();

    m_ = m;
    n_ = n;
}


void MatrixSparse::empty()
{
    if (mat_)
    {
        sparse_status_t status = mkl_sparse_destroy(mat_);
        assert(!status);
        mat_ = 0;
    }

    if (ijs1_)
    {
        if (isCsrOwned_)
        {
            mkl_free(ijs_);
            mkl_free(js_);
            mkl_free(vs_);
        }
        ijs1_ = 0;
    }

    isCsrOwned_ = false;
}


void MatrixSparse::swap(MatrixSparse& a)
{
    MatrixSparse t = *this;
    *this = a;
    a = t;
    t.ijs1_ = 0;
    t.mat_ = 0;
}


ii MatrixSparse::m() const
{
    return m_;
}


ii MatrixSparse::n() const
{
    return n_;
}


li MatrixSparse::size() const
{
    return li(m_) * n_;
}


ii MatrixSparse::nnz() const
{
    if (ijs1_)
        return ijs1_[m_ - 1];
    else
        return 0;
}


ii MatrixSparse::nnzActual() const
{
    ii count = 0;

    for (ii nz = 0; nz < nnz(); nz++)
    {
        if (vs_[nz] != 0.0)
            count++;
    }

    return count;
}


void MatrixSparse::copy(const MatrixSparse& a)
{
    if (getDebugLevel() % 10 >= 4)
    {
        ostringstream oss;
        oss << getTimeStamp() << "       A" << a << " := ...";
        info(oss.str());
    }

    assert(this != &a);

    if (createCsr(a.m(), a.n(), a.nnz()))
    {
        //IppStatus status = ippsCopy_32s(a.ijs_, ijs_, a.m_ + 1);
        IppStatus status = vectorCopy(a.ijs_, ijs_, a.m_ + 1);
        assert(!status);

        //status = ippsCopy_32s(a.js_, js_, a.nnz());
        status = vectorCopy(a.js_, js_, a.nnz());
        assert(!status);

        //status = ippsCopy_32f(a.vs_, vs_, a.nnz());
        status = vectorCopy(a.vs_, vs_, a.nnz());
        assert(!status);

        commitCsr(a.isSorted_);
    }

    if (getDebugLevel() % 10 >= 4)
    {
        ostringstream oss;
        oss << getTimeStamp() << "       ... X" << *this;
        info(oss.str(), this);
    }
}


void MatrixSparse::transpose(const MatrixSparse& a)
{
    if (getDebugLevel() % 10 >= 4)
    {
        ostringstream oss;
        oss << getTimeStamp() << "       t(A" << a << ") := ...";
        info(oss.str());
    }

    assert(this != &a);

    if (createMkl(a.n(), a.m(), a.nnz() > 0))
    {
        sparse_status_t status = mkl_sparse_convert_csr(a.mat_, SPARSE_OPERATION_TRANSPOSE, &mat_);
        assert(!status);

        commitMkl(true);
    }

    if (getDebugLevel() % 10 >= 4)
    {
        ostringstream oss;
        oss << getTimeStamp() << "       ... X" << *this;
        info(oss.str(), this);
    }
}


void MatrixSparse::importFromCoo(ii a_m, ii a_n, ii a_nnz, const ii* a_is, const ii* a_js, const fp* a_vs)
{
    if (getDebugLevel() % 10 >= 4)
    {
        ostringstream oss;
        oss << getTimeStamp() << "       importFrom(COO) := ...";
        info(oss.str());
    }

    if (createMkl(a_m, a_n, a_nnz > 0))
    {
        sparse_matrix_t t;
        ii* _a_is = const_cast<ii*>(a_is);
        ii* _a_js = const_cast<ii*>(a_js);
        fp* _a_vs = const_cast<fp*>(a_vs);

        sparse_status_t status = mkl_sparse_s_create_coo(&t, SPARSE_INDEX_BASE_ZERO, m_, n_, a_nnz, _a_is, _a_js, _a_vs);
        assert(!status);

        status = mkl_sparse_convert_csr(t, SPARSE_OPERATION_NON_TRANSPOSE, &mat_);
        assert(!status);

        status = mkl_sparse_destroy(t);
        assert(!status);

        commitMkl(true);
    }

    if (getDebugLevel() % 10 >= 4)
    {
        ostringstream oss;
        oss << getTimeStamp() << "       ... X" << *this;
        info(oss.str(), this);
    }
}


// todo: optimize
void MatrixSparse::importFromMatrix(const Matrix &a)
{
    empty();

    vector<ii> is;
    vector<ii> js;
    vector<fp> vs;

    for (ii i = 0; i < a.m(); i++)
    {
        for (ii j = 0; j < a.n(); j++)
        {
            if (a.vs()[j + i * a.n()] != 0.0)
            {
                is.push_back(i);
                js.push_back(j);
                vs.push_back(a.vs()[j + i * a.n()]);
            }
        }
    }

    importFromCoo(a.m(), a.n(), vs.size(), is.data(), js.data(), vs.data());
}


// todo: optimize
void MatrixSparse::initDense(ii m, ii n, fp v)
{
    empty();

    vector<ii> is;
    vector<ii> js;
    vector<fp> vs;

    for (ii i = 0; i < m; i++)
    {
        for (ii j = 0; j < n; j++)
        {
            is.push_back(i);
            js.push_back(j);
            vs.push_back(v);
        }
    }

    importFromCoo(m, n, vs.size(), is.data(), js.data(), vs.data());
}


void MatrixSparse::concatenateRows(const std::vector<MatrixSparse> &as)
{
    if (as.size() > 0)
    {
        if (getDebugLevel() % 10 >= 4)
        {
            ostringstream oss;
            oss << getTimeStamp() << "       " << "concatenateRows(A" << as.front() << " x " << as.size() << ")" << " := ...";
            info(oss.str());
        }

        assert(as.size() > 0);

        for (size_t k = 0; k < as.size(); k++)
            assert(as[k].m_ == 1);

        for (size_t k = 0; k < as.size() - 1; k++)
            assert(as[k].n_ == as[k + 1].n_);

        ii as_nnz = 0;
        for (ii i = 0; i < ii(as.size()); i++)
            as_nnz += as[i].nnz();

        if (createCsr(ii(as.size()), as[0].n(), as_nnz))
        {
            ijs_[0] = 0;
            for (ii i = 0; i < m_; i++)
                ijs1_[i] = ijs_[i] + as[i].nnz();

            //#pragma omp parallel
            for (ii i = 0; i < m_; i++)
            {
                if (as[i].ijs1_)
                {
                    //IppStatus status = ippsCopy_32s(as[i].js_, &js_[ijs_[i]], as[i].nnz());
                    IppStatus status = vectorCopy(as[i].js_, &js_[ijs_[i]], as[i].nnz());
                    assert(!status);
                }
             }

            //#pragma omp parallel
            for (ii i = 0; i < m_; i++)
            {
                if (as[i].ijs1_)
                {
                    //IppStatus status = ippsCopy_32f(as[i].vs_, &vs_[ijs_[i]], as[i].nnz());
                    IppStatus status = vectorCopy(as[i].vs_, &vs_[ijs_[i]], as[i].nnz());
                    assert(!status);
                }
            }

            bool isSorted = true;
            for (ii k = 0; k < ii(as.size()); k++)
            {
                if (!as[k].isSorted_)
                {
                    isSorted = false;
                    break;
                }
            }
            commitCsr(isSorted);
        }

        if (getDebugLevel() % 10 >= 4)
        {
            ostringstream oss;
            oss << getTimeStamp() << "       ... X" << *this;
            info(oss.str(), this);
        }
    }
}


void MatrixSparse::copySubset(const MatrixSparse &a, const MatrixSparse &b)
{
    if (getDebugLevel() % 10 >= 4)
    {
        ostringstream oss;
        oss << getTimeStamp() << "       copySubset(A" << a << " within B" << b << ") := ...";
        info(oss.str());
    }

    assert(this != &a);
    assert(this != &b);
    assert(a.m() == b.m());
    assert(a.n() == b.n());

    if (createCsr(b.m(), b.n(), b.nnz()))
    {
        a.sort();
        b.sort();

        //IppStatus status = ippsCopy_32s(b.ijs_, ijs_, m_ + 1);
        IppStatus status = vectorCopy(b.ijs_, ijs_, m_ + 1);
        assert(!status);

        //status = ippsCopy_32s(b.js_, js_, nnz());
        status = vectorCopy(b.js_, js_, nnz());
        assert(!status);

        for (ii i = 0; i < m_; i++)
        {
            ii a_nz = a.ijs_[i];
            for (ii nz = ijs_[i]; nz < ijs1_[i]; nz++)
            {
                for (; a_nz < a.ijs1_[i]; a_nz++)
                {
                    if (b.js_[nz] == a.js_[a_nz])
                    {
                        vs_[nz] = a.vs_[a_nz];
                        break;
                    }
                }
            }
        }

        commitCsr(true);
    }

    if (getDebugLevel() % 10 >= 4)
    {
        ostringstream oss;
        oss << getTimeStamp() << "       ... X" << *this;
        info(oss.str(), this);
    }
}


ii MatrixSparse::pruneCells(const MatrixSparse &a, fp threshold)
{
    if (getDebugLevel() % 10 >= 4)
    {
        ostringstream oss;
        oss << getTimeStamp() << "       pruneCells(A" << a << " <= ";
        oss.unsetf(ios::floatfield);
        oss << setprecision(8) << threshold << ") := ...";
        info(oss.str());
    }

    assert(this != &a);

    ii nnzCells = 0;
    vector<ii> nnzs;
    if (a.nnz() > 0)
    {
        nnzs.resize(a.m_, 0);
        for (ii i = 0; i < a.m_; i++)
        {
            for (ii a_nz = a.ijs_[i]; a_nz < a.ijs1_[i]; a_nz++)
            {
                if (a.vs_[a_nz] > threshold)
                    nnzs[i]++;
            }
            nnzCells += nnzs[i];
        }
    }

    if (createCsr(a.m_, a.n_, nnzCells))
    {
        ijs_[0] = 0;
        for (ii i = 0; i < a.m_; i++)
            ijs1_[i] = ijs_[i] + nnzs[i];

        //#pragma omp parallel
        for (ii i = 0; i < a.m_; i++)
        {
            ii nz = ijs_[i];
            for (ii a_nz = a.ijs_[i]; a_nz < a.ijs1_[i]; a_nz++)
            {
                if (a.vs_[a_nz] > threshold)
                {
                    js_[nz] = a.js_[a_nz];
                    vs_[nz] = a.vs_[a_nz];
                    nz++;
                }
            }
        }

        commitCsr(a.isSorted_);
    }

    if (getDebugLevel() % 10 >= 4)
    {
        ostringstream oss;
        oss << getTimeStamp() << "       ... X" << *this;
        info(oss.str(), this);
    }

    return nnzCells;
}


// THIS WOULD BE BETTER DONE FAKE-INPLACE SO NO PRUNING IS A NULL OP
ii MatrixSparse::pruneRows(const MatrixSparse &a, const MatrixSparse &b, bool bRows, fp threshold)
{
    if (getDebugLevel() % 10 >= 4)
    {
        ostringstream oss;
        oss << getTimeStamp() << "       pruneRows(" << a << ",";
        oss << (bRows ? "rows(" : "columns(") << b << ")) where ";
        oss << (bRows ? "nRows" : "nColumns") << " to prune > " << fixed << setprecision(1) << threshold * 100.0 << "% := ...";
        info(oss.str());
    }

    assert(this != &a);
    assert(this != &b);
    assert(bRows ? a.m_ == b.m_ : a.m_ == b.n_);

    init(a.m_, a.n_);

    ii rowsPruned = 0;
    if (a.nnz() > 0)
    {
        ii aNnzRows = 0;
        for (ii i = 0; i < a.m_; i++)
            if (a.ijs1_[i] - a.ijs_[i] > 0)
                aNnzRows++;

        if (b.nnz() > 0)
        {
            vector<ii> rowOrColNnzs(a.m_, 0);
            if (bRows)
            {
                for (ii i = 0; i < b.m_; i++)
                    rowOrColNnzs[i] = b.ijs1_[i] - b.ijs_[i];
            }
            else
            {
                for (ii nz = 0; nz < b.nnz(); nz++)
                    rowOrColNnzs[b.js_[nz]]++;
            }

            ii bNnzRowsOrCols = 0;
            for (ii i = 0; i < a.m_; i++)
            {
                if (rowOrColNnzs[i] > 0)
                    bNnzRowsOrCols++;
            }

            if (bNnzRowsOrCols / fp(aNnzRows) < threshold)
            {
                ii bNnzRowsOrCols = 0;
                ii outNnz = 0;
                for (ii i = 0; i < a.m_; i++)
                {
                    if (rowOrColNnzs[i] > 0)
                        outNnz += a.ijs1_[i] - a.ijs_[i];
                }

                createCsr(m_, n_, outNnz);

                ijs_[0] = 0;
                for (ii i = 0; i < m_; i++)
                    ijs1_[i] = ijs_[i] + (rowOrColNnzs[i] > 0 ? a.ijs1_[i] - a.ijs_[i] : 0);

                //#pragma omp parallel
                for (ii i = 0; i < m_; i++)
                {
                    if (ijs1_[i] - ijs_[i] > 0)
                    {
                        //IppStatus status = ippsCopy_32s(&a.js_[a.ijs_[i]], &js_[ijs_[i]], ijs1_[i] - ijs_[i]);
                        IppStatus status = vectorCopy(&a.js_[a.ijs_[i]], &js_[ijs_[i]], ijs1_[i] - ijs_[i]);
                        assert(!status);
                    }
                }

                //#pragma omp parallel
                for (ii i = 0; i < m_; i++)
                {
                    if (ijs1_[i] - ijs_[i] > 0)
                    {
                        //IppStatus status = ippsCopy_32f(&a.vs_[a.ijs_[i]], &vs_[ijs_[i]], ijs1_[i] - ijs_[i]);
                        IppStatus status = vectorCopy(&a.vs_[a.ijs_[i]], &vs_[ijs_[i]], ijs1_[i] - ijs_[i]);
                        assert(!status);
                    }
                }

                commitCsr(a.isSorted_);

                rowsPruned = aNnzRows- bNnzRowsOrCols;
            }
        }
        else
        {
            rowsPruned = aNnzRows;
        }
    }

    if (getDebugLevel() % 10 >= 4)
    {
        ostringstream oss;
        oss << getTimeStamp() << "       ... X" << *this << " (" << rowsPruned << " rows pruned)";
        info(oss.str(), this);
    }

    return rowsPruned;
}


void MatrixSparse::exportToCoo(ii *is, ii *js, fp *vs) const
{
    if (getDebugLevel() % 10 >= 4)
    {
        ostringstream oss;
        oss << getTimeStamp() << "       X" << *this << " := ...";
        info(oss.str());
    }

    if (initCsr(nnz() > 0))
    {
        ii job[] = { 0, 0, 0, 0 , ijs1_[m_ - 1], 3 };
        ii info;
        mkl_scsrcoo(job, &m_, vs_, js_, ijs_, &ijs1_[m_ - 1], vs, is, js, &info);
    }

    if (getDebugLevel() % 10 >= 4)
    {
        ostringstream oss;
        oss << getTimeStamp() << "       ... exportTo(COO)";
        info(oss.str(), this);
    }
}


void MatrixSparse::exportToDense(fp *vs) const
{
    if (getDebugLevel() % 10 >= 4)
    {
        ostringstream oss;
        oss << getTimeStamp() << "       X" << *this << " := ...";
        info(oss.str());
    }

    if (initCsr(nnz() > 0))
    {
        ii job[] = { 1, 0, 0, 2, 0, 0 };
        ii info;
        mkl_sdnscsr(job, &m_, &n_, vs, &n_, vs_, js_, ijs_, &info);
    }
    else
    {
        memset(vs, 0, sizeof(fp) * size());
    }

    if (getDebugLevel() % 10 >= 4)
    {
        ostringstream oss;
        oss << getTimeStamp() << "       ... exportTo(DENSE)";
        info(oss.str(), this);
    }
}


void MatrixSparse::add(fp alpha, bool transposeA, const MatrixSparse& a, const MatrixSparse& b)
{
    if (getDebugLevel() % 10 >= 4)
    {
        ostringstream oss;
        oss << getTimeStamp() << "       " << alpha << " * " << (transposeA ? "t(" : "") << "A" << a << (transposeA ? ")" : "") << " + B" << b;
        oss << " := ...";
        info(oss.str());
    }

    assert(this != &a);
    assert(this != &b);
    assert((transposeA ? a.n() : a.m()) == b.m());
    assert((transposeA ? a.m() : a.n()) == b.n());

    if (a.nnz() > 0 && b.nnz() == 0)
    {
        if (transposeA)
            transpose(a);
        else
            copy(a);
    }
    else if (a.nnz() == 0 && b.nnz() > 0)
    {
        copy(b);
    }
    else if (createMkl(transposeA ? a.n() : a.m(), b.n(), a.nnz() > 0 && b.nnz() > 0))
    {
        sparse_status_t status = mkl_sparse_s_add(transposeA ? SPARSE_OPERATION_TRANSPOSE : SPARSE_OPERATION_NON_TRANSPOSE, a.mat_, alpha, b.mat_, &mat_); assert(!status);

        commitMkl(true);
    }

    if (getDebugLevel() % 10 >= 4)
    {
        ostringstream oss;
        oss << getTimeStamp() << "       ... X" << *this;
        info(oss.str(), this);
    }
}


void MatrixSparse::matmul(bool transposeA, const MatrixSparse& a, const MatrixSparse& b, bool accumulate)
{
    if (getDebugLevel() % 10 >= 4)
    {
        ostringstream oss;
        oss << getTimeStamp() << "       " << (transposeA ? "t(" : "") << "A" << a << (transposeA ? ")" : "") << " %*% B" << b;
        if (accumulate) oss << " + X" << *this;
        oss << " := ...";
        info(oss.str());
    }

    assert(this != &a);
    assert(this != &b);
    assert((transposeA ? a.m() : a.n()) == b.m());

    if (nnz() == 0)
        accumulate = false;

    if (accumulate)
    {
        if (initMkl(a.nnz() > 0 && b.nnz() > 0))
        {
            sparse_matrix_t t;
            sparse_status_t status = mkl_sparse_spmm(transposeA ? SPARSE_OPERATION_TRANSPOSE : SPARSE_OPERATION_NON_TRANSPOSE, a.mat_, b.mat_, &t);
            assert(!status);

            sparse_matrix_t y;
            status = mkl_sparse_s_add(SPARSE_OPERATION_NON_TRANSPOSE, mat_, 1.0, t, &y);
            assert(!status);

            status = mkl_sparse_destroy(t);
            assert(!status);

            empty();
            mat_ = y;

            commitMkl(false);
        }
    }
    else
    {
        if (createMkl(transposeA ? a.n() : a.m(), b.n(), a.nnz() > 0 && b.nnz() > 0))
        {
             sparse_status_t status = mkl_sparse_spmm(transposeA ? SPARSE_OPERATION_TRANSPOSE : SPARSE_OPERATION_NON_TRANSPOSE, a.mat_, b.mat_, &mat_);
            assert(!status);

            commitMkl(false);
        }
    }

    if (getDebugLevel() % 10 >= 4)
    {
        ostringstream oss;
        oss << getTimeStamp() << "       ... X" << *this;
        oss;
        info(oss.str(), this);
    }
}


void MatrixSparse::matmulDense(bool transposeA, const MatrixSparse &a, const MatrixSparse &b)
{
    if (getDebugLevel() % 10 >= 4)
    {
        ostringstream oss;
        oss << getTimeStamp() << "       " << (transposeA ? "t(" : "") << "A" << a << (transposeA ? ")" : "") << " %*% B" << b;
        oss << " := ...";
        info(oss.str());
    }

    assert(this != &a);
    assert(this != &b);
    assert((transposeA ? a.m() : a.n()) == b.m());

    if (createCsr(transposeA ? a.n() : a.m(), b.n(), (transposeA ? a.n() : a.m()) * b.n()))
    {
        ijs_[0] = 0;
        for (ii i = 0; i < m_; i++)
            ijs1_[i] = ijs_[i] + n_;

        for (ii nz = 0; nz < nnz(); nz++)
            js_[nz] = nz % n_;

        if (a.nnz() > 0 && b.nnz() > 0)
        {
            sparse_status_t status = mkl_sparse_s_spmmd(transposeA ? SPARSE_OPERATION_TRANSPOSE : SPARSE_OPERATION_NON_TRANSPOSE, a.mat_, b.mat_, SPARSE_LAYOUT_ROW_MAJOR, vs_, n_);
            assert(!status);
        }
        else
        {
            memset(vs_, 0, sizeof(fp) * size());
        }

        commitCsr(true);
    }


    if (getDebugLevel() % 10 >= 4)
    {
        ostringstream oss;
        oss << getTimeStamp() << "       ... X" << *this << " (DENSE)";
        oss;
        info(oss.str(), this);
    }
}


void MatrixSparse::mul(fp beta)
{
    if (getDebugLevel() % 10 >= 4)
    {
        ostringstream oss;
        oss << getTimeStamp() << "       X" << *this << " * ";
        oss.unsetf(ios::floatfield);
        oss << setprecision(8) << beta << " := ...";
        info(oss.str());
    }

    if (initCsr(nnz() > 0))
    {
        IppStatus status;
        status = ippsMulC_32f_I(beta, vs_, ijs1_[m_ - 1]);
        assert(!status);

        commitCsr(isSorted_);
    }

    if (getDebugLevel() % 10 >= 4)
    {
        ostringstream oss;
        oss << getTimeStamp() << "       ... X" << *this;
        info(oss.str(), this);
    }
}


void MatrixSparse::mul(const MatrixSparse& a, const MatrixSparse& b)
{
    if (getDebugLevel() % 10 >= 4)
    {
        ostringstream oss;
        oss << getTimeStamp() << "       " << (this == &a ? "X" : "A") << a << " * " << (this == &b ? "X" : "B") << b << " := ...";
        info(oss.str());
    }

    if (this != &a && this != &b && createCsr(a.m(), a.n(), a.nnz()))
    {
        //IppStatus status = ippsCopy_32s(a.ijs_, ijs_, m_ + 1);
        IppStatus status = vectorCopy(a.ijs_, ijs_, m_ + 1);
        assert(!status);

        //status = ippsCopy_32s(a.js_, js_, nnz());
        status = vectorCopy(a.js_, js_, nnz());
        assert(!status);
    }

    if (initCsr(nnz() > 0))
    {
        if (&a != &b)
        {
            a.sort();
            b.sort();
        }

        assert(a.nnz() == b.nnz());
        for (ii i = 0; i < m_; i++)
            assert(a.ijs_[i] == b.ijs_[i]);
        for (ii nz = 0; nz < ijs1_[m_ - 1]; nz++)
            assert(a.js_[nz] == b.js_[nz]);

        vsMul(nnz(), a.vs_, b.vs_, vs_);
        int err = vmlGetErrStatus();
        assert(err == VML_STATUS_OK);

        commitCsr(true);
    }

    if (getDebugLevel() % 10 >= 4)
    {
        ostringstream oss;
        oss << getTimeStamp() << "       ... X" << *this;
        info(oss.str(), this);
    }
}


void MatrixSparse::sqr(const MatrixSparse& a)
{
    if (getDebugLevel() % 10 >= 4)
    {
        ostringstream oss;
        oss << getTimeStamp() << "       sqr(" << (this == &a ? "X" : "A") << a << ") := ...";
        info(oss.str());
    }

    if (this != &a && createCsr(a.m(), a.n(), a.nnz()))
    {
        //IppStatus status = ippsCopy_32s(a.ijs_, ijs_, m_ + 1);
        IppStatus status = vectorCopy(a.ijs_, ijs_, m_ + 1);
        assert(!status);

        //status = ippsCopy_32s(a.js_, js_, nnz());
        status = vectorCopy(a.js_, js_, nnz());
        assert(!status);
    }

    if (initCsr(nnz() > 0))
    {
        vsSqr(nnz(), a.vs_, vs_);
        int err = vmlGetErrStatus();
        assert(err == VML_STATUS_OK);

        commitCsr(a.isSorted_);
    }

    if (getDebugLevel() % 10 >= 4)
    {
        ostringstream oss;
        oss << getTimeStamp() << "       ... X" << *this;
        info(oss.str(), this);
    }
}


void MatrixSparse::sqrt(const MatrixSparse& a)
{
    if (getDebugLevel() % 10 >= 4)
    {
        ostringstream oss;
        oss << getTimeStamp() << "       sqrt(" << (this == &a ? "X" : "A") << a << ") := ...";
        info(oss.str());
    }

    if (this != &a && createCsr(a.m(), a.n(), a.nnz()))
    {
        //IppStatus status = ippsCopy_32s(a.ijs_, ijs_, m_ + 1);
        IppStatus status = vectorCopy(a.ijs_, ijs_, m_ + 1);
        assert(!status);

        //status = ippsCopy_32s(a.js_, js_, nnz());
        status = vectorCopy(a.js_, js_, nnz());
        assert(!status);
    }

    if (initCsr(nnz() > 0))
    {
        vsSqrt(nnz(), a.vs_, vs_);
        int err = vmlGetErrStatus();
        assert(err == VML_STATUS_OK);

        commitCsr(a.isSorted_);
    }

    if (getDebugLevel() % 10 >= 4)
    {
        ostringstream oss;
        oss << getTimeStamp() << "       ... X" << *this;
        info(oss.str(), this);
    }
}


void MatrixSparse::pow(const MatrixSparse& a, fp power)
{
    if (getDebugLevel() % 10 >= 4)
    {
        ostringstream oss;
        oss << getTimeStamp() << "       pow(" << (this == &a ? "X" : "A") << a << ", " << power << ") := ...";
        info(oss.str());
    }

    if (this != &a && createCsr(a.m(), a.n(), a.nnz()))
    {
        //IppStatus status = ippsCopy_32s(a.ijs_, ijs_, m_ + 1);
        IppStatus status = vectorCopy(a.ijs_, ijs_, m_ + 1);
        assert(!status);

        //status = ippsCopy_32s(a.js_, js_, ijs1_[m_ - 1]);
        status = vectorCopy(a.js_, js_, ijs1_[m_ - 1]);
        assert(!status);
    }

    if (initCsr(nnz() > 0))
    {
        vsPowx(nnz(), a.vs_, power, vs_);
        int err = vmlGetErrStatus();
        assert(err == VML_STATUS_OK);

        commitCsr(a.isSorted_);
    }

    if (getDebugLevel() % 10 >= 4)
    {
        ostringstream oss;
        oss << getTimeStamp() << "       ... X" << *this;
        info(oss.str(), this);
    }
}


void MatrixSparse::censorLeft(const MatrixSparse& a, fp threshold)
{
    if (getDebugLevel() % 10 >= 4)
    {
        ostringstream oss;
        oss << getTimeStamp() << "       censorLeft(" << (this == &a ? "X" : "A") << a << ", ";
        oss.unsetf(ios::floatfield);
        oss << setprecision(8) << threshold << ") := ...";
        info(oss.str());
    }

    if (this != &a && createCsr(a.m(), a.n(), a.nnz()))
    {
        //IppStatus status = ippsCopy_32s(a.ijs_, ijs_, m_ + 1);
        IppStatus status = vectorCopy(a.ijs_, ijs_, m_ + 1);
        assert(!status);

        //status = ippsCopy_32s(a.js_, js_, nnz());
        status = vectorCopy(a.js_, js_, nnz());
        assert(!status);
    }

    if (initCsr(nnz() > 0))
    {
        IppStatus status = ippsThreshold_LT_32f(a.vs_, vs_, nnz(), threshold);
        assert(!status);

        commitCsr(a.isSorted_);
    }

    if (getDebugLevel() % 10 >= 4)
    {
        ostringstream oss;
        oss << getTimeStamp() << "       ... X" << *this;
        info(oss.str(), this);
    }
}


void MatrixSparse::addNonzeros(fp beta)
{
    if (getDebugLevel() % 10 >= 4)
    {
        ostringstream oss;
        oss << getTimeStamp() << "       X" << *this << " + ";
        oss.unsetf(ios::floatfield);
        oss << setprecision(8) << beta << " := ...";
        info(oss.str());
    }

    if (initCsr(nnz() > 0))
    {
        IppStatus status = ippsAddC_32f_I(beta, vs_, nnz());
        assert(!status);
    }

    if (getDebugLevel() % 10 >= 4)
    {
        ostringstream oss;
        oss << getTimeStamp() << "       ... X" << *this;
        info(oss.str(), this);
    }
}


void MatrixSparse::addNonzeros(const MatrixSparse& a, const MatrixSparse& b)
{
    if (getDebugLevel() % 10 >= 4)
    {
        ostringstream oss;
        oss << getTimeStamp() << "       " << (this == &a ? "X" : "A") << a << " / " << (this == &b ? "X" : "B") << b << " := ...";
        info(oss.str());
    }

    if (this != &a && this != &b && createCsr(a.m(), a.n(), a.nnz()))
    {
        //IppStatus status = ippsCopy_32s(a.ijs_, ijs_, m_ + 1);
        IppStatus status = vectorCopy(a.ijs_, ijs_, m_ + 1);
        assert(!status);

        //status = ippsCopy_32s(a.js_, js_, ijs1_[m_ - 1]);
        status = vectorCopy(a.js_, js_, ijs1_[m_ - 1]);
        assert(!status);
    }

    if (initCsr(nnz() > 0))
    {
        if (&a != &b)
        {
            a.sort();
            b.sort();
        }

        assert(a.nnz() == b.nnz());
        for (ii i = 0; i < m_; i++)
            assert(a.ijs_[i] == b.ijs_[i]);
        for (ii nz = 0; nz < ijs1_[m_ - 1]; nz++)
            assert(a.js_[nz] == b.js_[nz]);

        vsAdd(nnz(), a.vs_, b.vs_, vs_);
        int err = vmlGetErrStatus();
        assert(err == VML_STATUS_OK);

        commitCsr(true);
    }

    if (getDebugLevel() % 10 >= 4)
    {
        ostringstream oss;
        oss << getTimeStamp() << "       ... X" << *this;
        info(oss.str(), this);
    }
}


void MatrixSparse::lnNonzeros(const MatrixSparse& a)
{
    if (getDebugLevel() % 10 >= 4)
    {
        ostringstream oss;
        oss << getTimeStamp() << "       ln(" << (this == &a ? "X" : "A") << a << ") := ...";
        info(oss.str());
    }

    if (this != &a && createCsr(a.m(), a.n(), a.nnz()))
    {
        //IppStatus status = ippsCopy_32s(a.ijs_, ijs_, m_ + 1);
        IppStatus status = vectorCopy(a.ijs_, ijs_, m_ + 1);
        assert(!status);

        //status = ippsCopy_32s(a.js_, js_, ijs1_[m_ - 1]);
        status = vectorCopy(a.js_, js_, ijs1_[m_ - 1]);
        assert(!status);
    }

    if (initCsr(nnz() > 0))
    {
//#ifdef NDEBUG // for some reason this causes valgrind to crash
        vsLn(nnz(), a.vs_, vs_);
        int err = vmlGetErrStatus();
        assert(err == VML_STATUS_OK);
//#else
//        for (ii i = 0; i < ijs1_[m_ - 1]; i++)
//            vs_[i] = log(a.vs_[i]);
//#endif

        commitCsr(a.isSorted_);
    }

    if (getDebugLevel() % 10 >= 4)
    {
        ostringstream oss;
        oss << getTimeStamp() << "       ... X" << *this;
        info(oss.str(), this);
    }
}


void MatrixSparse::divNonzeros(const MatrixSparse& a, const MatrixSparse& b)
{
    if (getDebugLevel() % 10 >= 4)
    {
        ostringstream oss;
        oss << getTimeStamp() << "       " << (this == &a ? "X" : "A") << a << " / " << (this == &a ? "X" : "A") << a << " := ...";
        info(oss.str());
    }

    if (this != &a && this != &b && createCsr(a.m(), a.n(), a.nnz()))
    {
        //IppStatus status = ippsCopy_32s(a.ijs_, ijs_, m_ + 1);
        IppStatus status = vectorCopy(a.ijs_, ijs_, m_ + 1);
        assert(!status);

        //status = ippsCopy_32s(a.js_, js_, ijs1_[m_ - 1]);
        status = vectorCopy(a.js_, js_, ijs1_[m_ - 1]);
        assert(!status);
    }

    if (initCsr(nnz() > 0))
    {
        if (&a != &b)
        {
            a.sort();
            b.sort();
        }

        assert(a.nnz() == b.nnz());
        for (ii i = 0; i < m_; i++)
            assert(a.ijs_[i] == b.ijs_[i]);
        for (ii nz = 0; nz < ijs1_[m_ - 1]; nz++)
            assert(a.js_[nz] == b.js_[nz]);

        vsDiv(nnz(), a.vs_, b.vs_, vs_);
        int err = vmlGetErrStatus();
        assert(err == VML_STATUS_OK);

        commitCsr(true);
    }

    if (getDebugLevel() % 10 >= 4)
    {
        ostringstream oss;
        oss << getTimeStamp() << "       ... X" << *this;
        info(oss.str(), this);
    }
}


void MatrixSparse::div2Dense(const Matrix &a)
{
    if (getDebugLevel() % 10 >= 4)
    {
        ostringstream oss;
        oss << getTimeStamp() << "       A" << a << " / X" << *this << " := ...";
        info(oss.str());
    }

    sort();
    assert(nnz() == size());
    assert(m_ = a.m());
    assert(n_ = a.n());

    if (initCsr(nnz() > 0))
    {
        vsDiv(nnz(), a.vs(), vs_, vs_);
        int err = vmlGetErrStatus();
        assert(err == VML_STATUS_OK);

        commitCsr(true);
    }

    if (getDebugLevel() % 10 >= 4)
    {
        ostringstream oss;
        oss << getTimeStamp() << "       ... X" << *this;
        info(oss.str(), this);
    }
}


fp MatrixSparse::sum() const
{
    if (getDebugLevel() % 10 >= 4)
    {
        ostringstream oss;
        oss << getTimeStamp() << "       sum(X" << *this << ") := ...";
        info(oss.str());
    }

    fp sum = 0.0;
    if (initCsr(nnz() > 0))
    {
        IppStatus status = ippsSum_32f(vs_, nnz(), &sum, ippAlgHintFast);
        assert(!status);
    }

    if (getDebugLevel() % 10 >= 4)
    {
        ostringstream oss;
        oss << getTimeStamp() << "       ... ";
        oss.unsetf(ios::floatfield);
        oss << setprecision(8) << sum;
        info(oss.str(), this);
    }

    return sum;
}


fp MatrixSparse::sumSqrs() const
{
   if (getDebugLevel() % 10 >= 4)
    {
        ostringstream oss;
        oss << getTimeStamp() << "       sumSqrs(X" << *this << ") := ...";
        info(oss.str());
    }

    fp sum = 0.0;
    if (initCsr(nnz() > 0))
    {
        IppStatus status = ippsNorm_L2_32f(vs_, nnz(), &sum);
        assert(!status);

        sum *= sum;
    }

    if (getDebugLevel() % 10 >= 4)
    {
        ostringstream oss;
        oss << getTimeStamp() << "       ... ";
        oss.unsetf(ios::floatfield);
        oss << setprecision(8) << sum;
        info(oss.str(), this);
    }

    return sum;
}


fp MatrixSparse::sumSqrDiffsNonzeros(const MatrixSparse& a) const
{
    if (getDebugLevel() % 10 >= 4)
    {
        ostringstream oss;
        oss << getTimeStamp() << "       sumSqrDiffs(X" << *this << ", A" << a << ") := ...";
        info(oss.str());
    }

    fp sum = 0.0;
    if (initCsr(nnz() > 0))
    {
        sort();
        a.sort();

        assert(nnz() == a.nnz());
        for (ii i = 0; i < m_; i++)
            assert(ijs_[i] == a.ijs_[i]);
        for (ii nz = 0; nz < ijs1_[m_ - 1]; nz++)
            assert(js_[nz] == a.js_[nz]);

        IppStatus status = ippsNormDiff_L2_32f(vs_, a.vs_, nnz(), &sum);
        assert(!status);

        sum *= sum;
    }

    if (getDebugLevel() % 10 >= 4)
    {
        ostringstream oss;
        oss << getTimeStamp() << "       ... ";
        oss.unsetf(ios::floatfield);
        oss << setprecision(8) << sum;
        info(oss.str(), this);
    }
    
    return sum;
}


double MatrixSparse::sortElapsed_ = 0.0;


bool MatrixSparse::createCsr(ii m, ii n, ii nnz)
{
    init(m, n);

    if (nnz > 0)
    {
        ijs_ = static_cast<ii*>(mkl_malloc(sizeof(ii) * (m_ + 1), 64));
        ijs_[m_] = nnz;
        ijs1_ = ijs_ + 1;

        js_ = static_cast<ii*>(mkl_malloc(sizeof(ii) * nnz, 64));
        vs_ = static_cast<fp*>(mkl_malloc(sizeof(fp) * nnz, 64));

        isCsrOwned_ = true;

        // todo: defer this MKL setup to initMkl
        sparse_status_t status = mkl_sparse_s_create_csr(&mat_, SPARSE_INDEX_BASE_ZERO, m_, n_, ijs_, ijs1_, js_, vs_);
        assert(!status);

        return true;
    }
    else
    {
        return false;
    }
}


bool MatrixSparse::initCsr(bool notEmpty) const
{
    // todo: make sure Csr is setup
    return notEmpty;
}


void MatrixSparse::commitCsr(bool isSorted) const
{
    // this function sorts in place and I'd like everything else to think the object hasn't changed
    bool& _isSorted_ = const_cast<bool&>(isSorted_);
    _isSorted_ = isSorted;
}


bool MatrixSparse::createMkl(ii m, ii n, bool notEmpty)
{
    init(m, n);

    isCsrOwned_ = false;

    return notEmpty;
}


bool MatrixSparse::initMkl(bool notEmpty) const
{
    // todo: make sure Mkl is setup
    return notEmpty;
}


void MatrixSparse::commitMkl(bool isSorted) const
{
    // this function sorts in place and I'd like everything else to think the object hasn't changed
    bool& _isSorted_ = const_cast<bool&>(isSorted_);
    _isSorted_ = isSorted;

    // todo: defer this CSR setup to initCsr
    ii& _m_  = const_cast<ii&>(m_);
    ii& _n_  = const_cast<ii&>(n_);
    ii*& _ijs_  = const_cast<ii*&>(ijs_);
    ii*& _ijs1_  = const_cast<ii*&>(ijs1_);
    ii*& _js_  = const_cast<ii*&>(js_);
    fp*& _vs_  = const_cast<fp*&>(vs_);

    sparse_index_base_t indexing;
    sparse_status_t status;
    status = mkl_sparse_s_export_csr(mat_, &indexing, &_m_, &_n_, &_ijs_, &_ijs1_, &_js_, &_vs_);
    assert(!status);
    assert(!indexing);
}


// I've tried sorting with mkl_scsradd, ippsSortIndexAscend_32s_I & std::sort, this is fastest
void MatrixSparse::sort() const
{
    if (!isSorted_)
    {
        if (initCsr(nnz() > 0))
        {
            // check if we really need to sort
            bool isSorted = true;
            for (ii i = 0; i < m_; i++)
            {
                for (ii nz = ijs_[i]; nz < ijs1_[i] - 1; nz++)
                {
                    if (js_[nz] > js_[nz + 1])
                    {
                        isSorted = false;
                        break;
                    }
                }

                if (!isSorted)
                    break;
            }

            if (!isSorted)
            {
                if (getDebugLevel() % 10 >= 4)
                {
                    ostringstream oss;
                    oss << getTimeStamp() << "       sort(X" << *this << ") := ...";
                    info(oss.str());
                }

                double sortStart = getElapsedTime();
                {
                    //#pragma omp for
                    for (ii i = 0; i < m_; i++)
                    {
                        if (ijs1_[i] > ijs_[i])
                        {
                            ii bufSize;
                            IppStatus status = ippsSortRadixIndexGetBufferSize(ijs1_[i] - ijs_[i], ipp32s, &bufSize);
                            assert(!status);

                            Ipp8u* buffer = static_cast<Ipp8u*>(mkl_malloc(sizeof(Ipp8u) * bufSize, 64));

                            ii* idxs = static_cast<ii*>(mkl_malloc(sizeof(ii) * (ijs1_[i] - ijs_[i]), 64));
                            status = ippsSortRadixIndexAscend_32s(&js_[ijs_[i]], sizeof(ii), idxs, ijs1_[i] - ijs_[i], buffer);
                            assert(!status);
                            mkl_free(buffer);

                            ii* newJs = static_cast<ii*>(mkl_malloc(sizeof(fp) * (ijs1_[i] - ijs_[i]), 64));
                            for (ii nz = 0; nz < ijs1_[i] - ijs_[i]; nz++)
                                newJs[nz] = js_[ijs_[i] + idxs[nz]];

                            fp* newVs = static_cast<fp*>(mkl_malloc(sizeof(fp) * (ijs1_[i] - ijs_[i]), 64));
                            vsPackV(ijs1_[i] - ijs_[i], &vs_[ijs_[i]], idxs, newVs);
                            int err = vmlGetErrStatus();
                            assert(err == VML_STATUS_OK);
                            mkl_free(idxs);

                            //status = ippsCopy_32s(newJs, &js_[ijs_[i]], ijs1_[i] - ijs_[i]);
                            status = vectorCopy(newJs, &js_[ijs_[i]], ijs1_[i] - ijs_[i]);
                            assert(!status);
                            mkl_free(newJs);

                            //status = ippsCopy_32f(newVs, &vs_[ijs_[i]], ijs1_[i] - ijs_[i]);
                            status = vectorCopy(newVs, &vs_[ijs_[i]], ijs1_[i] - ijs_[i]);
                            assert(!status);
                            mkl_free(newVs);
                        }
                    }
                }
                sortElapsed_ += getElapsedTime() - sortStart;

                commitCsr(true);

                if (getDebugLevel() % 10 >= 4)
                {
                    ostringstream oss;
                    oss << getTimeStamp() << "       ... X" << *this;
                    info(oss.str(), this);
                }
            }
        }
        else
        {
            // this function sorts in place and I'd like everything else to think the object hasn't changed
            bool& _isSorted_ = const_cast<bool&>(isSorted_);
            _isSorted_ = true;
        }
    }
}


ostream& operator<<(ostream& os, const MatrixSparse& a)
{
    if (a.m() == 0)
    {
        os << "[]";
    }
    else
    {
        os << "[" << a.m_ << "," << a.n_ << "]:" << a.nnz();
        
        if (MatrixSparse::getDebugLevel()% 10 >= 4)
        {
            os << "(" << a.nnzActual() << ")";
        }

        os << "/" << a.size() << ":";
        os.unsetf(ios::floatfield);
        os << setprecision(3) << 100.0 * a.nnz() / (double)a.size() << "%";
    }

    return  os;
}


MatrixSparseView::MatrixSparseView(const MatrixSparse &a, ii row) : isOwned_(false)
{
    if (getDebugLevel() % 10 >= 4)
    {
        ostringstream oss;
        oss << getTimeStamp() << "       A" << a << "[" << row << "] := ...";
        info(oss.str());
    }

    assert(row >= 0 && row < a.m_);

    init(1, a.n_);

    if (initCsr(a.nnz() > 0))
    {
        ii newNnz = a.ijs1_[row] - a.ijs_[row];

        if (newNnz > 0)
        {
            ijs_ = static_cast<ii*>(mkl_malloc(sizeof(ii) * 2, 64));
            ijs_[0] = 0;
            ijs1_ = ijs_ + 1;

            ijs1_[0] = newNnz;
            js_ = &a.js_[a.ijs_[row]];
            vs_ = &a.vs_[a.ijs_[row]];

            sparse_status_t status;
            status = mkl_sparse_s_create_csr(&mat_, SPARSE_INDEX_BASE_ZERO, m_, n_, ijs_, ijs1_, js_, vs_);
            assert(!status);

            isOwned_ = true;

            commitCsr(a.isSorted_);
        }
    }

    if (getDebugLevel() % 10 >= 4)
    {
        ostringstream oss;
        oss << getTimeStamp() << "       ... X" << *this;
        info(oss.str(), this);
    }
}


MatrixSparseView::~MatrixSparseView()
{
    if (isOwned_)
    {
        mkl_free(ijs_);
    }
}
