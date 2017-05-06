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
    free();
}


void MatrixSparse::init(ii m, ii n)
{
    free();

    m_ = m;
    n_ = n;
}


bool MatrixSparse::allocCsr(ii nnz)
{
    free();

    if (nnz > 0)
    {
        ijs_ = static_cast<ii*>(mkl_malloc(sizeof(ii) * (m_ + 1), 64));
        ijs1_ = ijs_ + 1;
        ijs1_[m_ - 1] = nnz;
        js_ = static_cast<ii*>(mkl_malloc(sizeof(ii) * ijs1_[m_ - 1], 64));
        vs_ = static_cast<fp*>(mkl_malloc(sizeof(fp) * ijs1_[m_ - 1], 64));
        isCsrOwned_ = true;

        return true;
    }
    else
    {
        return false;
    }
 }


void MatrixSparse::free()
{
    if (mat_)
    {
        status_ = mkl_sparse_destroy(mat_);
        assert(!status_);
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


bool MatrixSparse::getRidOfCsr() const
{
    if (ijs1_)
    {
        return true;
    }
    else if (mat_)
    {
        sparse_index_base_t indexing;
        ii m;
        ii n;
        ii* _ijs_ = const_cast<ii*>(ijs_);
        ii* _ijs1_ = const_cast<ii*>(ijs1_);
        ii* _js_ = const_cast<ii*>(js_);
        fp* _vs_ = const_cast<fp*>(vs_);

        sparse_status_t status = mkl_sparse_s_export_csr(mat_, &indexing, &m, &n, &_ijs_, &_ijs1_, &_js_, &_vs_);

        assert(!status_);
        assert(!indexing);
        assert(m_ == m);
        assert(n_ == n);

        return true;
    }
    else
    {
        return false;
    }
}


bool MatrixSparse::initCsr(ii m, ii n, ii nnz)
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

        return true;
    }
    else
    {
        return false;
    }
}


bool MatrixSparse::initMkl(ii m, ii n, ii nnz)
{
    init(m, n);

    return nnz > 0;
}


void MatrixSparse::commitCsr(bool isSorted)
{
    isSorted_ = isSorted;

    sparse_status_t status = mkl_sparse_s_create_csr(&mat_, SPARSE_INDEX_BASE_ZERO, m_, n_, ijs_, ijs1_, js_, vs_);
    assert(!status);
}


void MatrixSparse::commitMkl(bool isSorted)
{
    isSorted_ = isSorted;

    sparse_index_base_t indexing;
    status_ = mkl_sparse_s_export_csr(mat_, &indexing, &m_, &n_, &ijs_, &ijs1_, &js_, &vs_);
    assert(!status_);
    assert(!indexing);
}


void MatrixSparse::copy(const MatrixSparse& a)
{
    if (getDebugLevel() % 10 >= 4)
    {
        ostringstream oss;
        oss << getTimeStamp() << "       A" << a << " := ...";
        info(oss.str());
    }

    if (initCsr(a.m_, a.n_, a.nnz()))
    {
        ippsCopy_32s(a.ijs_, ijs_, a.m_);
        ippsCopy_32s(a.js_, js_, a.nnz());
        ippsCopy_32f(a.vs_, vs_, a.nnz());

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

    if (initMkl(a.m(), a.n(), a.nnz()))
    {
        mkl_sparse_convert_csr(a.mat_, SPARSE_OPERATION_TRANSPOSE, &mat_);

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

    if (initMkl(a_m, a_n, a_nnz))
    {
        sparse_matrix_t t;
        ii* _a_is = const_cast<ii*>(a_is);
        ii* _a_js = const_cast<ii*>(a_js);
        fp* _a_vs = const_cast<fp*>(a_vs);
        status_ = mkl_sparse_s_create_coo(&t, SPARSE_INDEX_BASE_ZERO, m(), n(), a_nnz, _a_is, _a_js, _a_vs);
        assert(!status_);

        status_ = mkl_sparse_convert_csr(t, SPARSE_OPERATION_NON_TRANSPOSE, &mat_);
        assert(!status_);

        status_ = mkl_sparse_destroy(t);
        assert(!status_);

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
    free();

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
void MatrixSparse::importFromMatrix(ii m, ii n, fp v)
{
    free();

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


void MatrixSparse::concatenateSparseVectors(const std::vector<MatrixSparse> &as)
{
    if (as.size() > 0)
    {
        if (getDebugLevel() % 10 >= 4)
        {
            ostringstream oss;
            oss << getTimeStamp() << "       " << "concatenateSparseVectors(A" << as.front() << " x " << as.size() << ")" << " := ...";
            info(oss.str());
        }

        for (size_t k = 0; k < as.size(); k++)
            assert(as[k].m_ == 1);

        for (size_t k = 0; k < as.size() - 1; k++)
            assert(as[k].n_ == as[k + 1].n_);

        init(ii(as.size()), as[0].n());

        ii as_nnz = 0;
        for (ii i = 0; i < m_; i++)
            as_nnz += as[i].nnz();

        if (initCsr(ii(as.size()), as[0].n(), as_nnz))
        {
            ijs_[0] = 0;
            for (ii i = 0; i < m_ - 1; i++)
                ijs_[i + 1] = ijs_[i] + as[i].nnz();

            //#pragma omp parallel
            for (ii i = 0; i < m_; i++)
            {
                if (as[i].ijs1_)
                    ippsCopy_32s(as[i].js_, &js_[ijs_[i]], as[i].ijs1_[0]);
             }

            //#pragma omp parallel
            for (ii i = 0; i < m_; i++)
            {
                if (as[i].ijs1_)
                    ippsCopy_32f(as[i].vs_, &vs_[ijs_[i]], as[i].ijs1_[0]);
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


// SEEMS OPTIMAL
void MatrixSparse::copySubset(const MatrixSparse &a)
{
    if (getDebugLevel() % 10 >= 4)
    {
        ostringstream oss;
        oss << getTimeStamp() << "       copySubset(A" << a << ") within X" << *this << " := ...";
        info(oss.str());
    }

    assert(m_ == a.m_);
    assert(n_ == a.n_);

    if (ijs1_ && a.ijs1_)
    {
        sort();
        a.sort();

        //#pragma omp parallel
        for (ii i = 0; i < m_; i++)
        {
            ii a_nz = a.ijs_[i];
            for (ii nz = ijs_[i]; nz < ijs1_[i]; nz++)
            {
                for (; a_nz < a.ijs1_[i]; a_nz++)
                {
                    if (js_[nz] == a.js_[a_nz])
                    {
                        vs_[nz] = a.vs_[a_nz];
                        break;
                    }
                }
             }
        }

        isSorted_ = true;
    }

    if (getDebugLevel() % 10 >= 4)
    {
        ostringstream oss;
        oss << getTimeStamp() << "       ... X" << *this;
        info(oss.str(), this);
    }
}


// SEEMS OPTIMAL
void MatrixSparse::copySubset(const MatrixSparse &a, const MatrixSparse &b)
{
    if (getDebugLevel() % 10 >= 4)
    {
        ostringstream oss;
        oss << getTimeStamp() << "       copySubset(A" << a << " within B" << b << ") := ...";
        info(oss.str());
    }

    assert(a.m() == b.m());
    assert(a.n() == b.n());

    init(a.m(), a.n());

    if (a.getRidOfCsr() && b.getRidOfCsr())
    {
        a.sort();
        b.sort();

        allocCsr(b.nnz());

        ippsCopy_32s(b.ijs_, ijs_, m_);
        ippsCopy_32s(b.js_, js_, nnz());

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

        isSorted_ = true;

        status_ = mkl_sparse_s_create_csr(&mat_, SPARSE_INDEX_BASE_ZERO, m_, n_, ijs_, ijs1_, js_, vs_);
        assert(!status_);
    }

    if (getDebugLevel() % 10 >= 4)
    {
        ostringstream oss;
        oss << getTimeStamp() << "       ... X" << *this;
        info(oss.str(), this);
    }
}


ii MatrixSparse::copyPrune(const MatrixSparse &a, fp threshold)
{
    if (getDebugLevel() % 10 >= 4)
    {
        ostringstream oss;
        oss << getTimeStamp() << "       copyPrune(A" << a << " <= ";
        oss.unsetf(ios::floatfield);
        oss << setprecision(8) << threshold << ") := ...";
        info(oss.str());
    }

    init(a.m_, a.n_);

    ii nnzCells = 0;
    if (a.ijs1_)
    {
        vector<ii> nnzs(a.m_, 0);
        for (ii i = 0; i < a.m_; i++)
        {
            for (ii a_nz = a.ijs_[i]; a_nz < a.ijs1_[i]; a_nz++)
            {
                if (a.vs_[a_nz] > threshold)
                    nnzs[i]++;
            }
            nnzCells += nnzs[i];
        }

        if (nnzCells > 0)
        {
            ijs_ = static_cast<ii*>(mkl_malloc(sizeof(ii) * (a.m_ + 1), 64));
            ijs_[0] = 0;
            ijs1_ = ijs_ + 1;
            for (ii i = 0; i < a.m_; i++)
                ijs1_[i] = ijs_[i] + nnzs[i];

            js_ = static_cast<ii*>(mkl_malloc(sizeof(ii) * ijs1_[a.m_ - 1], 64));
            vs_ = static_cast<fp*>(mkl_malloc(sizeof(fp) * ijs1_[a.m_ - 1], 64));
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

            status_ = mkl_sparse_s_create_csr(&mat_, SPARSE_INDEX_BASE_ZERO, m_, n_, ijs_, ijs1_, js_, vs_);
            assert(!status_);

            isCsrOwned_ = true;
            isSorted_ = a.isSorted_;
        }
    }

    if (getDebugLevel() % 10 >= 4)
    {
        ostringstream oss;
        oss << getTimeStamp() << "       ... X" << *this;
        info(oss.str(), this);
    }

    return nnzCells;
}


// MOSTLY OPTIMAL, MIGHT BE BETTER IF OMP LOOP DIDN'T INCLUDE PRUNED ROWS
ii MatrixSparse::copyPruneRows(const MatrixSparse &a, const MatrixSparse &b, bool bRows, fp threshold)
{
    if (getDebugLevel() % 10 >= 4)
    {
        ostringstream oss;
        oss << getTimeStamp() << "       copyPruneRows(" << a << ",";
        oss << (bRows ? "rows(" : "columns(") << b << ")) where ";
        oss << (bRows ? "nRows" : "nColumns") << " to prune > " << fixed << setprecision(1) << threshold * 100.0 << "% := ...";
        info(oss.str());
    }

    assert(bRows ? a.m_ == b.m_ : a.m_ == b.n_);

    init(a.m_, a.n_);

    ii rowsPruned = 0;
    if (a.ijs1_ && b.ijs1_)
    {
        ii aNnzRows = 0;
        for (ii i = 0; i < m_; i++)
            if (a.ijs1_[i] - a.ijs_[i] > 0)
                aNnzRows++;
        //oss << "aNnzRows=" << aNnzRows;

        vector<ii> rowOrColNnzs(m_, 0);
        if (bRows)
        {
            for (ii i = 0; i < m_; i++)
                rowOrColNnzs[i] = b.ijs1_[i] - b.ijs_[i];
        }
        else
        {
            for (ii nz = 0; nz < b.nnz(); nz++)
                rowOrColNnzs[b.js_[nz]]++;
        }

        ii bNnzRowsOrCols = 0;
        for (ii i = 0; i < m_; i++)
            if (rowOrColNnzs[i] > 0)
                bNnzRowsOrCols++;
        //oss << "bNnzRowsOrCols=" << bNnzRowsOrCols;

        if (bNnzRowsOrCols / (fp) aNnzRows < threshold)
        {
            ijs_ = static_cast<ii*>(mkl_malloc(sizeof(ii) * (m_ + 1), 64));
            ijs_[0] = 0;
            ijs1_ = ijs_ + 1;
            for (ii i = 0; i < m_; i++)
                ijs1_[i] = ijs_[i] + (rowOrColNnzs[i] > 0 ? a.ijs1_[i] - a.ijs_[i] : 0);

            js_ = static_cast<ii*>(mkl_malloc(sizeof(ii) * ijs1_[m_ - 1], 64));
            vs_ = static_cast<fp*>(mkl_malloc(sizeof(fp) * ijs1_[m_ - 1], 64));
            //#pragma omp parallel
            for (ii i = 0; i < m_; i++)
            {
                ippsCopy_32s(&a.js_[a.ijs_[i]], &js_[ijs_[i]], ijs1_[i] - ijs_[i]);
                ippsCopy_32f(&a.vs_[a.ijs_[i]], &vs_[ijs_[i]], ijs1_[i] - ijs_[i]);
            }

            status_ = mkl_sparse_s_create_csr(&mat_, SPARSE_INDEX_BASE_ZERO, m_, n_, ijs_, ijs1_, js_, vs_);
            assert(!status_);

            isCsrOwned_ = true;
            isSorted_ = a.isSorted_;

            rowsPruned = aNnzRows- bNnzRowsOrCols;
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


void MatrixSparse::exportTo(ii* rowind, ii* colind, fp* acoo) const
{
    if (getDebugLevel() % 10 >= 4)
    {
        ostringstream oss;
        oss << getTimeStamp() << "       X" << *this << " := ...";
        info(oss.str());
    }

    ii length = nnz();
    if (length > 0)
    {
        ii job[] = { 0, 0, 0, 0 , length, 3 };
        ii info;
        mkl_scsrcoo(job, &m_, vs_, js_, ijs_, &length, acoo, rowind, colind, &info);
    }

    if (getDebugLevel() % 10 >= 4)
    {
        ostringstream oss;
        oss << getTimeStamp() << "       ... exportTo(COO)";
        info(oss.str(), this);
    }
}


void MatrixSparse::exportTo(fp *vs) const
{
    if (getDebugLevel() % 10 >= 4)
    {
        ostringstream oss;
        oss << getTimeStamp() << "       X" << *this << " := ...";
        info(oss.str());
    }

    if (nnz() > 0)
    {
        ii job[] = { 1, 0, 0, 2, 0, 0 };
        ii info;
        mkl_sdnscsr(job, &m_, &n_, vs, &n_, vs_, js_, ijs_, &info);
    }
    else
    {
        for (li x = 0; x < size(); x++)
            vs[x] = 0.0;
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

    assert((transposeA ? a.n() : a.m()) == b.m());
    assert((transposeA ? a.m() : a.n()) == b.n());

    if (!a.ijs1_ || !b.ijs1_)
    {
        init(transposeA ? a.n() : a.m(), b.n());
    }
    else if (!a.ijs1_)
    {
        copy(b);
    }
    else if (!b.ijs1_)
    {
        if (transposeA)
            transpose(a);
        else
            copy(a);
    }
    else
    {
        init(transposeA ? a.n() : a.m(), b.n());

        status_ = mkl_sparse_s_add(transposeA ? SPARSE_OPERATION_TRANSPOSE : SPARSE_OPERATION_NON_TRANSPOSE, a.mat_, alpha, b.mat_, &mat_); assert(!status_);
        
        sparse_index_base_t indexing;
        status_ = mkl_sparse_s_export_csr(mat_, &indexing, &m_, &n_, &ijs_, &ijs1_, &js_, &vs_);
        assert(!status_);

        isCsrOwned_ = false;
        isSorted_ = true;
    }

    if (getDebugLevel() % 10 >= 4)
    {
        ostringstream oss;
        oss << getTimeStamp() << "       ... X" << *this;
        info(oss.str(), this);
    }
}


void MatrixSparse::matmul(bool transposeA, const MatrixSparse& a, const MatrixSparse& b, bool accumulate, bool denseOutput)
{
    if (getDebugLevel() % 10 >= 4)
    {
        ostringstream oss;
        oss << getTimeStamp() << "       " << (transposeA ? "t(" : "") << "A" << a << (transposeA ? ")" : "") << " %*% B" << b;
        if (accumulate) oss << " + X" << *this;
        oss << " := ...";
        info(oss.str());
    }

    assert((transposeA ? a.m() : a.n()) == b.m());

    if (!ijs1_)
        accumulate = false;

    if (a.ijs1_ && b.ijs1_)
    {
        if (!accumulate)
            init(transposeA ? a.n() : a.m(), b.n());

        if (denseOutput)
        {
            ii m = transposeA ? a.n() : a.m();
            ii n = b.n();

            ii* is0 = static_cast<ii*>(mkl_malloc(sizeof(ii) * (m + 1), 64));
            is0[0] = 0;
            ii* is1 = is0 + 1;
            for (ii i = 0; i < m; i++)
                is1[i] = is0[i] + n;

            ii* js = static_cast<ii*>(mkl_malloc(sizeof(ii) * is1[m - 1], 64));
            for (ii nz = 0; nz < is1[m - 1]; nz++)
                js[nz] = nz % n;

            fp* vs = static_cast<fp*>(mkl_malloc(sizeof(fp) * is1[m - 1], 64));
            status_ = mkl_sparse_s_spmmd(transposeA ? SPARSE_OPERATION_TRANSPOSE : SPARSE_OPERATION_NON_TRANSPOSE, a.mat_, b.mat_, SPARSE_LAYOUT_ROW_MAJOR, vs, n); assert(!status_);

            sparse_matrix_t t;
            status_ = mkl_sparse_s_create_csr(&t, SPARSE_INDEX_BASE_ZERO, m, n, is0, is1, js, vs);
            assert(!status_);

            if (accumulate)
            {
                sparse_matrix_t y;
                status_ = mkl_sparse_s_add(SPARSE_OPERATION_NON_TRANSPOSE, mat_, 1.0, t, &y);
                assert(!status_);

                status_ = mkl_sparse_destroy(t);
                assert(!status_);

                free();
                mat_ = y;

                sparse_index_base_t indexing;
                status_ = mkl_sparse_s_export_csr(mat_, &indexing, &m_, &n_, &ijs_, &ijs1_, &js_, &vs_);
                assert(!status_);

                isCsrOwned_ = false;
            }
            else
            {
                m_ = m;
                n_ = n;
                ijs_ = is0;
                ijs1_ = is1;
                js_ = js;
                vs_ = vs;
                mat_ = t;
                isCsrOwned_ = true;
            }

            isSorted_ = true;
        }
        else
        {
            if (accumulate)
            {
                sparse_matrix_t t;
                status_ = mkl_sparse_spmm(transposeA ? SPARSE_OPERATION_TRANSPOSE : SPARSE_OPERATION_NON_TRANSPOSE,
                                          a.mat_, b.mat_, &t);

                // annoying hack: it fails with this error when the output has no non-zeros
                if (status_ != SPARSE_STATUS_ALLOC_FAILED)
                {
                    assert(!status_);

                    sparse_matrix_t y;
                    status_ = mkl_sparse_s_add(SPARSE_OPERATION_NON_TRANSPOSE, mat_, 1.0, t, &y);
                    assert(!status_);

                    status_ = mkl_sparse_destroy(t);
                    assert(!status_);

                    free();
                    mat_ = y;

                    sparse_index_base_t indexing;
                    status_ = mkl_sparse_s_export_csr(mat_, &indexing, &m_, &n_, &ijs_, &ijs1_, &js_, &vs_);
                    assert(!status_);
                }
            }
            else
            {
                status_ = mkl_sparse_spmm(transposeA ? SPARSE_OPERATION_TRANSPOSE : SPARSE_OPERATION_NON_TRANSPOSE,
                                          a.mat_, b.mat_, &mat_);

                // annoying hack: it fails with this error when the output has no non-zeros
                if(status_ == SPARSE_STATUS_ALLOC_FAILED)
                {
                    init(transposeA ? a.n() : a.m(), b.n());
                }
                else
                {
                    assert(!status_);

                    sparse_index_base_t indexing;
                    status_ = mkl_sparse_s_export_csr(mat_, &indexing, &m_, &n_, &ijs_, &ijs1_, &js_, &vs_);
                    assert(!status_);
                }
            }

            isCsrOwned_ = false;
            isSorted_ = false;
        }
    }
    else
    {
        if (!accumulate)
        {
            if (denseOutput)
                importFromMatrix(transposeA ? a.n() : a.m(), b.n(), fp(0.0));
            else
                init(transposeA ? a.n() : a.m(), b.n());
        }
    }

    if (getDebugLevel() % 10 >= 4)
    {
        ostringstream oss;
        oss << getTimeStamp() << "       ... X" << *this;
        if (denseOutput) oss << " (DENSE)";
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

    if (ijs1_)
        ippsMulC_32f_I(beta, vs_, ijs1_[m_ - 1]);

    if (getDebugLevel() % 10 >= 4)
    {
        ostringstream oss;
        oss << getTimeStamp() << "       ... X" << *this;
        info(oss.str(), this);
    }
}


void MatrixSparse::mul(const MatrixSparse& a)
{
    if (getDebugLevel() % 10 >= 4)
    {
        ostringstream oss;
        oss << getTimeStamp() << "       X" << *this << " * A := ...";
        info(oss.str());
    }

    if (ijs1_)
    {
        sort();
        a.sort();

        assert(ijs1_[m_ - 1] == a.ijs1_[m_ - 1]);
        for (ii i = 0; i < m_; i++)
            assert(ijs_[i] == a.ijs_[i]);
        for (ii nz = 0; nz < ijs1_[m_ - 1]; nz++)
            assert(js_[nz] == a.js_[nz]);

        vsMul(ijs1_[m_ - 1], vs_, a.vs_, vs_);
    }

    if (getDebugLevel() % 10 >= 4)
    {
        ostringstream oss;
        oss << getTimeStamp() << "       ... X" << *this;
        info(oss.str(), this);
    }
}


void MatrixSparse::sqr()
{
    if (getDebugLevel() % 10 >= 4)
    {
        ostringstream oss;
        oss << getTimeStamp() << "       sqr(X" << *this << ") := ...";
        info(oss.str());
    }

    if (ijs1_)
        vsSqr(ijs1_[m_ - 1], vs_, vs_);

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
        oss << getTimeStamp() << "       sqr(A" << a << ") := ...";
        info(oss.str());
    }

    init(a.m_, a.n_);

    if (a.ijs1_)
    {
        ijs_ = static_cast<ii*>(mkl_malloc(sizeof(ii) * (m_ + 1), 64));
        ippsCopy_32s(a.ijs_, ijs_, m_ + 1);
        ijs1_ = ijs_ + 1;
        js_ = static_cast<ii*>(mkl_malloc(sizeof(ii) * ijs1_[m_ - 1], 64));
        ippsCopy_32s(a.js_, js_, ijs1_[m_ - 1]);
        vs_ = static_cast<fp*>(mkl_malloc(sizeof(fp) * ijs1_[m_ - 1], 64));

        vsSqr(ijs1_[m_ - 1], a.vs_, vs_);

        status_ = mkl_sparse_s_create_csr(&mat_, SPARSE_INDEX_BASE_ZERO, m_, n_, ijs_, ijs1_, js_, vs_);
        assert(!status_);

        isCsrOwned_ = true;
        isSorted_ = a.isSorted_;
    }

    if (getDebugLevel() % 10 >= 4)
    {
        ostringstream oss;
        oss << getTimeStamp() << "       ... X" << *this;
        info(oss.str(), this);
    }
}


void MatrixSparse::sqrt()
{
    if (getDebugLevel() % 10 >= 4)
    {
        ostringstream oss;
        oss << getTimeStamp() << "       sqrt(X" << *this << ") := ...";
        info(oss.str());
    }

    if (ijs1_)
        vsSqrt(ijs1_[m_ - 1], vs_, vs_);

    if (getDebugLevel() % 10 >= 4)
    {
        ostringstream oss;
        oss << getTimeStamp() << "       ... X" << *this;
        info(oss.str(), this);
    }
}


void MatrixSparse::pow(fp power)
{
    if (getDebugLevel() % 10 >= 4)
    {
        ostringstream oss;
        oss << getTimeStamp() << "       pow(X" << *this << ", " << power << ") := ...";
        info(oss.str());
    }

    if (ijs1_)
        vsPowx(ijs1_[m_ - 1], vs_, power, vs_);

    if (getDebugLevel() % 10 >= 4)
    {
        ostringstream oss;
        oss << getTimeStamp() << "       ... X" << *this;
        info(oss.str(), this);
    }
}


void MatrixSparse::censorLeft(fp threshold)
{
    if (getDebugLevel() % 10 >= 4)
    {
        ostringstream oss;
        oss << getTimeStamp() << "       censorLeft(X" << *this << ", ";
        oss.unsetf(ios::floatfield);
        oss << setprecision(8) << threshold << ") := ...";
        info(oss.str());
    }

    if (ijs1_)
        ippsThreshold_LT_32f(vs_, vs_, ijs1_[m_ - 1], threshold);

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

    if (ijs1_)
        ippsAddC_32f_I(beta, vs_, ijs1_[m_ - 1]);

    if (getDebugLevel() % 10 >= 4)
    {
        ostringstream oss;
        oss << getTimeStamp() << "       ... X" << *this;
        info(oss.str(), this);
    }
}


void MatrixSparse::addNonzeros(const MatrixSparse& a)
{
    if (getDebugLevel() % 10 >= 4)
    {
        ostringstream oss;
        oss << getTimeStamp() << "       X" << *this << " / A" << a << " := ...";
        info(oss.str());
    }

    if (ijs1_)
    {
        sort();
        a.sort();

        assert(ijs1_[m_ - 1] == a.ijs1_[m_ - 1]);
        for (ii i = 0; i < m_; i++)
            assert(ijs_[i] == a.ijs_[i]);
        for (ii nz = 0; nz < ijs1_[m_ - 1]; nz++)
            assert(js_[nz] == a.js_[nz]);

        vsAdd(ijs1_[m_ - 1], vs_, a.vs_, vs_);
    }

    if (getDebugLevel() % 10 >= 4)
    {
        ostringstream oss;
        oss << getTimeStamp() << "       ... X" << *this;
        info(oss.str(), this);
    }
}



void MatrixSparse::lnNonzeros()
{
    if (getDebugLevel() % 10 >= 4)
    {
        ostringstream oss;
        oss << getTimeStamp() << "       ln(X" << *this << ") := ...";
        info(oss.str());
    }

    if (ijs1_)
    {
#ifdef NDEBUG
        vsLn(is1_[m_ - 1], vs_, vs_); // for some reason this causes valgrind to crash
#else
        for (ii i = 0; i < ijs1_[m_ - 1]; i++)
            vs_[i] = log(vs_[i]);
#endif
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
        oss << getTimeStamp() << "       ln(X" << a << ") := ...";
        info(oss.str());
    }

    init(a.m_, a.n_);

    if (a.ijs1_)
    {
        ijs_ = static_cast<ii*>(mkl_malloc(sizeof(ii) * (m_ + 1), 64));
        ippsCopy_32s(a.ijs_, ijs_, m_ + 1);
        ijs1_ = ijs_ + 1;
        js_ = static_cast<ii*>(mkl_malloc(sizeof(ii) * ijs1_[m_ - 1], 64));
        ippsCopy_32s(a.js_, js_, ijs1_[m_ - 1]);
        vs_ = static_cast<fp*>(mkl_malloc(sizeof(fp) * ijs1_[m_ - 1], 64));

#ifdef NDEBUG
        vsLn(is1_[m_ - 1], a.vs_, vs_); // for some reason this causes valgrind to crash
#else
        for (ii i = 0; i < ijs1_[m_ - 1]; i++)
            vs_[i] = log(a.vs_[i]);
#endif

        status_ = mkl_sparse_s_create_csr(&mat_, SPARSE_INDEX_BASE_ZERO, m_, n_, ijs_, ijs1_, js_, vs_);
        assert(!status_);

        isCsrOwned_ = true;
        isSorted_ = a.isSorted_;
    }

    if (getDebugLevel() % 10 >= 4)
    {
        ostringstream oss;
        oss << getTimeStamp() << "       ... X" << *this;
        info(oss.str(), this);
    }
}


void MatrixSparse::expNonzeros()
{
    if (getDebugLevel() % 10 >= 4)
    {
        ostringstream oss;
        oss << getTimeStamp() << "       exp(X" << *this << ") := ...";
        info(oss.str());
    }

    if (ijs1_)
        vsExp(ijs1_[m_ - 1], vs_, vs_);

    if (getDebugLevel() % 10 >= 4)
    {
        ostringstream oss;
        oss << getTimeStamp() << "       ... X" << *this;
        info(oss.str(), this);
    }
}


void MatrixSparse::divNonzeros(const MatrixSparse& a)
{
    if (getDebugLevel() % 10 >= 4)
    {
        ostringstream oss;
        oss << getTimeStamp() << "       X" << *this << " / A" << a << " := ...";
        info(oss.str());
    }

    if (ijs1_)
    {
        sort();
        a.sort();

        assert(ijs1_[m_ - 1] == a.ijs1_[m_ - 1]);
        for (ii i = 0; i < m_; i++)
            assert(ijs_[i] == a.ijs_[i]);
        for (ii nz = 0; nz < ijs1_[m_ - 1]; nz++)
            assert(js_[nz] == a.js_[nz]);

        vsDiv(ijs1_[m_ - 1], vs_, a.vs_, vs_);
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
        oss << getTimeStamp() << "       A" << a << " / B" << b << " := ...";
        info(oss.str());
    }

    init(a.m_, a.n_);

    if (a.ijs1_ && b.ijs1_)
    {
        a.sort();
        b.sort();

        assert(a.ijs1_[m_ - 1] == b.ijs1_[m_ - 1]);
        for (ii i = 0; i < m_; i++)
            assert(a.ijs_[i] == b.ijs_[i]);
        for (ii nz = 0; nz < a.ijs1_[m_ - 1]; nz++)
            assert(a.js_[nz] == b.js_[nz]);

        ijs_ = static_cast<ii*>(mkl_malloc(sizeof(ii) * (m_ + 1), 64));
        ippsCopy_32s(a.ijs_, ijs_, m_ + 1);
        ijs1_ = ijs_ + 1;
        js_ = static_cast<ii*>(mkl_malloc(sizeof(ii) * ijs1_[m_ - 1], 64));
        ippsCopy_32s(a.js_, js_, ijs1_[m_ - 1]);
        vs_ = static_cast<fp*>(mkl_malloc(sizeof(fp) * ijs1_[m_ - 1], 64));

        vsDiv(ijs1_[m_ - 1], a.vs_, b.vs_, vs_);

        status_ = mkl_sparse_s_create_csr(&mat_, SPARSE_INDEX_BASE_ZERO, m_, n_, ijs_, ijs1_, js_, vs_);
        assert(!status_);

        isCsrOwned_ = true;
        isSorted_ = true;
    }

    if (getDebugLevel() % 10 >= 4)
    {
        ostringstream oss;
        oss << getTimeStamp() << "       ... X" << *this;
        info(oss.str(), this);
    }
}


void MatrixSparse::div2Nonzeros(const MatrixSparse& a)
{
    if (getDebugLevel() % 10 >= 4)
    {
        ostringstream oss;
        oss << getTimeStamp() << "       A" << a << " / X" << *this << " := ...";
        info(oss.str());
    }

    if (ijs1_)
    {
        sort();
        a.sort();

        assert(ijs1_[m_ - 1] == a.ijs1_[m_ - 1]);
        for (ii i = 0; i < m_; i++)
            assert(ijs_[i] == a.ijs_[i]);
        for (ii nz = 0; nz < ijs1_[m_ - 1]; nz++)
            assert(js_[nz] == a.js_[nz]);

        vsDiv(ijs1_[m_ - 1], a.vs_, vs_, vs_);
    }

    if (getDebugLevel() % 10 >= 4)
    {
        ostringstream oss;
        oss << getTimeStamp() << "       ... X" << *this;
        info(oss.str(), this);
    }
}


void MatrixSparse::div2(const Matrix &a)
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

    if (ijs1_)
    {
        vsDiv(ijs1_[m_ - 1], a.vs(), vs_, vs_);

        //for (ii i = 0; i < is1_[m_ - 1]; i++)
        //    vs_[i] = vs_[i] > 0.0 ? a.vs()[i] / vs_[i] : 0.0;
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
    if (ijs1_)
        ippsSum_32f(vs_, ijs1_[m_ - 1], &sum, ippAlgHintFast);

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
    if (ijs1_)
    {
        ippsNorm_L2_32f(vs_, ijs1_[m_ - 1], &sum);
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
    if (ijs1_)
    {
        sort();
        a.sort();

        assert(ijs1_[m_ - 1] == a.ijs1_[m_ - 1]);
        for (ii i = 0; i < m_; i++)
            assert(ijs_[i] == a.ijs_[i]);
        for (ii nz = 0; nz < ijs1_[m_ - 1]; nz++)
            assert(js_[nz] == a.js_[nz]);

        ippsNormDiff_L2_32f(vs_, a.vs_, ijs1_[m_ - 1], &sum);
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


// I've tried sorting with mkl_scsradd, ippsSortIndexAscend_32s_I & std::sort, this is fastest
void MatrixSparse::sort() const
{
    // this function sorts in place and I'd like everything else to think the object hasn't changed
    bool& _isSorted_ = const_cast<bool&>(isSorted_);

    if (ijs1_)
    {
        if (!isSorted_ && ijs1_)
        {
            // check if we really need to sort
            _isSorted_ = true;
            for (ii i = 0; i < m_; i++)
            {
                for (ii nz = ijs_[i]; nz < ijs1_[i] - 1; nz++)
                {
                    if (js_[nz] > js_[nz + 1])
                    {
                        _isSorted_ = false;
                        break;
                    }
                }

                if (!isSorted_)
                    break;
            }

            if (!isSorted_)
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
                            ippsSortRadixIndexGetBufferSize(ijs1_[i] - ijs_[i], ipp32s, &bufSize);
                            Ipp8u* buffer = static_cast<Ipp8u*>(mkl_malloc(sizeof(Ipp8u) * bufSize, 64));

                            ii* idxs = static_cast<ii*>(mkl_malloc(sizeof(ii) * (ijs1_[i] - ijs_[i]), 64));
                            ippsSortRadixIndexAscend_32s(&js_[ijs_[i]], sizeof(ii), idxs, ijs1_[i] - ijs_[i], buffer);
                            mkl_free(buffer);

                            ii* newJs = static_cast<ii*>(mkl_malloc(sizeof(fp) * (ijs1_[i] - ijs_[i]), 64));
                            for (ii nz = 0; nz < ijs1_[i] - ijs_[i]; nz++)
                                newJs[nz] = js_[ijs_[i] + idxs[nz]];

                            fp* newVs = static_cast<fp*>(mkl_malloc(sizeof(fp) * (ijs1_[i] - ijs_[i]), 64));
                            vsPackV(ijs1_[i] - ijs_[i], &vs_[ijs_[i]], idxs, newVs);
                            mkl_free(idxs);

                            ippsCopy_32s(newJs, &js_[ijs_[i]], ijs1_[i] - ijs_[i]);
                            mkl_free(newJs);

                            ippsCopy_32f(newVs, &vs_[ijs_[i]], ijs1_[i] - ijs_[i]);
                            mkl_free(newVs);
                        }
                    }
                }
                sortElapsed_ += getElapsedTime() - sortStart;

                _isSorted_ = true;

                if (getDebugLevel() % 10 >= 4)
                {
                    ostringstream oss;
                    oss << getTimeStamp() << "       ... X" << *this;
                    info(oss.str(), this);
                }
            }
        }
    }
    else
    {
        _isSorted_ = true;
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

    if (a.nnz() > 0)
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

            status_ = mkl_sparse_s_create_csr(&mat_, SPARSE_INDEX_BASE_ZERO, m_, n_, ijs_, ijs1_, js_, vs_);
            assert(!status_);

            isOwned_ = true;
            isSorted_ = a.isSorted_;
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
