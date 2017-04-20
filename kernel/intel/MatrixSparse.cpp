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
#include "../FileNetcdf.hpp"
#include <iomanip>
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


MatrixSparse::MatrixSparse() : m_(0), n_(0), is1_(0)
{
}


MatrixSparse::~MatrixSparse()
{
    clear();
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
    return (li)m_ * n_;
}


ii MatrixSparse::nnz() const
{
    if (is1_)
    {
        return is1_[m_ - 1];
    }
    else
    {
        return 0;
    }
}


ii MatrixSparse::nnzActual() const
{
    if (is1_)
    {
        ii count = 0;

        for (ii nz = 0; nz < nnz(); nz++)
        {
            if (vs_[nz] != 0.0) count++;
        }

        return count;
    }
    else
    {
        return 0;
    }
}


fp* MatrixSparse::vs() const
{
    return vs_;
}


void MatrixSparse::init(ii m, ii n)
{
    clear();
    m_ = m;
    n_ = n;
}


void MatrixSparse::init(const MatrixSparse &a, ii row)
{
    assert(a.is1_);
    assert(a.nnz() > 0);

    init();

    if (getDebugLevel() % 10 >= 4)
        cout << getTimeStamp() << "       A" << a << "[" << row << "] := ..." << endl;

    m_ = 1;
    n_ = a.n();
    // we are not deallocating this at the moment... (and valgrind spots the small leak)
    is0_ = static_cast<ii*>(mkl_malloc(sizeof(ii) * 2, 64));
    is0_[0] = 0;
    is1_ = is0_ + 1;

    ii newNnz = a.is1_[row] - a.is0_[row];
    is1_[0] = newNnz;
    js_ = &a.js_[a.is0_[row]];
    vs_ = &a.vs_[a.is0_[row]];

    status_ = mkl_sparse_s_create_csr(&mat_, SPARSE_INDEX_BASE_ZERO, m_, n_, is0_, is1_, js_, vs_); assert(!status_);
    isOwned_ = false;

    if (getDebugLevel() % 10 >= 4)
        cout << getTimeStamp() << "       ... X" << *this << endl;
}


// SEEMS OPTIMAL
void MatrixSparse::copy(const MatrixSparse& a, bool transpose)
{
    init();

    if (getDebugLevel() % 10 >= 4)
        cout << getTimeStamp() << "       " << (transpose ? "t(" : "") << "A" << a << (transpose ? ")" : "") << " := ..." << endl;

    if (a.is1_)
    {
        if (transpose)
        {
            if (a.is1_)
            {
                mkl_sparse_convert_csr(a.mat_, SPARSE_OPERATION_TRANSPOSE, &mat_);

                sparse_index_base_t indexing;
                mkl_sparse_s_export_csr(mat_, &indexing, &m_, &n_, &is0_, &is1_, &js_, &vs_);
                isOwned_ = false;
            }
            else
            {
                m_ = a.n_;
                n_ = a.m_;
            }
        }
        else
        {
            m_ = a.m_;
            n_ = a.n_;

            if (a.is1_)
            {
                is0_ = static_cast<ii*>(mkl_malloc(sizeof(ii) * (a.m_ + 1), 64));
                is1_ = is0_ + 1;
                js_ = static_cast<ii*>(mkl_malloc(sizeof(ii) * a.is1_[m_ - 1], 64));
                vs_ = static_cast<fp*>(mkl_malloc(sizeof(fp) * a.is1_[m_ - 1], 64));

                memcpy(is0_, a.is0_, sizeof(ii) * (a.m_ + 1));
                memcpy(vs_, a.vs_, sizeof(fp) * a.is1_[m_ - 1]);
                memcpy(js_, a.js_, sizeof(ii) * a.is1_[m_ - 1]);

                status_ = mkl_sparse_s_create_csr(&mat_, SPARSE_INDEX_BASE_ZERO, m_, n_, is0_, is1_, js_, vs_); assert(!status_);
                isOwned_ = true;
            }
        }
    }
    else
    {
        if (transpose)
            init(a.n_, a.m_);
        else
            init(a.m_, a.n_);
    }


    if (getDebugLevel() % 10 >= 4)
        cout << getTimeStamp() << "       ... X" << *this << endl;
}


void MatrixSparse::copy(const std::vector<MatrixSparse> &as)
{
    clear();

    for (size_t k = 0; k < as.size(); k++) assert(as[k].m_ == 1);
    for (size_t k = 0; k < as.size() - 1; k++) assert(as[k].n_ == as[k + 1].n_);

    m_ = (ii)as.size();
    n_ = as[0].n();

    ii nnz = 0;
    for (ii i = 0; i < m_; i++)
        nnz += as[i].nnz();

    if (nnz > 0)
    {
        is0_ = static_cast<ii*>(mkl_malloc(sizeof(ii) * (m_ + 1), 64));
        is0_[0] = 0;
        is1_ = is0_ + 1;
        for (ii i = 0; i < m_; i++)
             is1_[i] = is0_[i] + as[i].nnz();

        js_ = static_cast<ii*>(mkl_malloc(sizeof(ii) * is1_[m_ - 1], 64));
        vs_ = static_cast<fp*>(mkl_malloc(sizeof(fp) * is1_[m_ - 1], 64));
        #pragma omp parallel
        for (ii i = 0; i < m_; i++)
        {
            if (as[i].is1_)
            {
                memcpy(&js_[is0_[i]], as[i].js_, sizeof(ii) * as[i].is1_[0]);
                memcpy(&vs_[is0_[i]], as[i].vs_, sizeof(fp) * as[i].is1_[0]);
            }
        }

        status_ = mkl_sparse_s_create_csr(&mat_, SPARSE_INDEX_BASE_ZERO, m_, n_, is0_, is1_, js_, vs_); assert(!status_);
        isOwned_ = true;
    }
}


// SEEMS OPTIMAL
void MatrixSparse::copy(ii m, ii n, const std::vector<ii> &is, const std::vector<ii> &js,
                        const std::vector<fp> &vs)
{
    clear();
    
    if (getDebugLevel() % 10 >= 4)
        cout << getTimeStamp() << "       COO := ..." << endl;

    m_ = m;
    n_ = n;

    if (vs.size() > 0)
    {
        fp* c_acoo = const_cast<fp*>(vs.data());
        ii* c_rowind = const_cast<ii*>(is.data());
        ii* c_colind = const_cast<ii*>(js.data());

        sparse_matrix_t t;
        status_ = mkl_sparse_s_create_coo(&t, SPARSE_INDEX_BASE_ZERO, m_, n_, vs.size(), c_rowind, c_colind, c_acoo); assert(!status_);
        status_ = mkl_sparse_convert_csr(t, SPARSE_OPERATION_NON_TRANSPOSE, &mat_); assert(!status_);
        status_ = mkl_sparse_destroy(t); assert(!status_);

        sparse_index_base_t indexing;
        status_ = mkl_sparse_s_export_csr(mat_, &indexing, &m_, &n_, &is0_, &is1_, &js_, &vs_); assert(!status_);
        isOwned_ = false;
    }
    
    if (getDebugLevel() % 10 >= 4)
        cout << getTimeStamp() << "       ... X" << *this << endl;
}


void MatrixSparse::alloc(ii m, ii n, fp v)
{
    init();

    vector<fp> vs;
    vector<ii> is;
    vector<ii> js;

    for (ii i = 0; i < m; i++)
    {
        for (ii j = 0; j < n; j++)
        {
            is.push_back(i);
            js.push_back(j);
            vs.push_back(v);
        }
    }

    copy(m, n, is, js, vs);
}


void MatrixSparse::copy(const Matrix &a)
{
    clear();

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

    copy(a.m(), a.n(), is, js, vs);
}


// SEEMS OPTIMAL
void MatrixSparse::copySubset(const MatrixSparse &a)
{
    if (getDebugLevel() % 10 >= 4)
        cout << getTimeStamp() << "       " << a << " within " << *this << " := ..." << endl;

    assert(m_ == a.m_);
    assert(n_ == a.n_);

    if (is1_ && a.is1_)
    {
        for (ii i = 0; i < m_; i++)
            for (ii nz = is0_[i]; nz < is1_[i] - 1; nz++)
                assert(js_[nz] <= js_[nz + 1]);
        for (ii i = 0; i < a.m_; i++)
            for (ii a_nz = a.is0_[i]; a_nz < a.is1_[i] - 1; a_nz++)
                assert(a.js_[a_nz] <= a.js_[a_nz + 1]);

#pragma omp parallel
        for (ii i = 0; i < m_; i++)
        {
            ii a_nz = a.is0_[i];
            for (ii nz = is0_[i]; nz < is1_[i]; nz++)
            {
                while(a_nz < a.is1_[i])
                {
                    if (js_[nz] == a.js_[a_nz])
                    {
                        vs_[nz] = a.vs_[a_nz];
                        a_nz++;
                        break;
                    }
                    else
                        a_nz++;
                }
            }
        }
    }

    if (getDebugLevel() % 10 >= 4)
        cout << getTimeStamp() << "       ... X" << *this << endl;
}


void MatrixSparse::exportTo(fp *vs) const
{
    if (getDebugLevel() % 10 >= 4)
        cout << getTimeStamp() << "       X" << *this << " := ..." << endl;

    if (nnz() > 0)
    {
        ii job[] = { 1, 0, 0, 2, 0, 0 };
        ii info;
        mkl_sdnscsr(job, &m_, &n_, vs, &n_, vs_, js_, is0_, &info);
    }
    else
    {
        for (li x = 0; x < size(); x++)
            vs[x] = 0.0;
    }

    if (getDebugLevel() % 10 >= 4)
        cout << getTimeStamp() << "       ... dense" << endl;
 }


void MatrixSparse::exportTo(vector<ii> &is, vector<ii> &js, vector<fp> &vs) const
{
    if (getDebugLevel() % 10 >= 4)
        cout << getTimeStamp() << "       X" << *this << " := ..." << endl;

    ii _nnz = nnz();
    vector<ii>(_nnz).swap(is);
    vector<ii>(_nnz).swap(js);
    vector<fp>(_nnz).swap(vs);

    if (_nnz > 0)
    {
        ii job[] = { 0, 0, 0, 0 , _nnz, 3 };
        ii info;
        mkl_scsrcoo(job, &m_, vs_, js_, is0_, &_nnz, vs.data(), is.data(), js.data(), &info);
    }

    if (getDebugLevel() % 10 >= 4)
        cout << getTimeStamp() << "       ... COO" << endl;
}


void MatrixSparse::clear()
{
    if (is1_)
    {
        status_ = mkl_sparse_destroy(mat_); assert(!status_);

        if (isOwned_)
        {
            mkl_free(is0_);
            mkl_free(js_);
            mkl_free(vs_);
        }

        is1_ = 0;
    }
}


struct MyComparator
{
    ii* js_;
    
    MyComparator(ii* js): js_(js) {}
    
    bool operator()(ii i1, ii i2)
    {
        return js_[i1] < js_[i2];
    }
};


// SHOULD BE IN-PLACE SORTING WITHOUT CREATION OF INDICIES (LOOK AT IPP SORT)
void MatrixSparse::sort()
{
    if (getDebugLevel() % 10 >= 4)
        cout << getTimeStamp() << "       sort(X" << *this << ") := ..." << endl;

    if (is1_)
    {
        //#pragma omp parallel // why does omp make it crash here?
        for (ii i = 0; i < m_; i++)
        {
            vector<ii> indicies(is1_[i] - is0_[i]);
            for (ii nz = 0; nz < indicies.size(); nz++)
                indicies[nz] = nz;
            std::sort(indicies.begin(), indicies.end(), MyComparator(&js_[is0_[i]]));
            
            {
                vector<ii> js(is1_[i] - is0_[i]);
                for (ii nz = 0; nz < js.size(); nz++)
                    js[nz] = js_[is0_[i] + indicies[nz]];
                for (ii nz = 0; nz < js.size(); nz++)
                    js_[is0_[i] + nz] = js[nz];
            }
            
            {
                vector<fp> vs(is1_[i] - is0_[i]);
                for (ii nz = 0; nz < vs.size(); nz++)
                    vs[nz] = vs_[is0_[i] + indicies[nz]];
                for (ii nz = 0; nz < vs.size(); nz++)
                    vs_[is0_[i] + nz] = vs[nz];
            }
        }
    }
    
    if (getDebugLevel() % 10 >= 4)
         cout << getTimeStamp() << "       ... X" << *this << endl;
}


ii MatrixSparse::prune(const MatrixSparse& a, fp pruneThreshold)
{
    if (getDebugLevel() % 10 >= 4)
    {
        cout << getTimeStamp() << "       (A" << a << " > ";
        cout.unsetf(ios::floatfield);
        cout << setprecision(8) << pruneThreshold << ") ? A : 0.0 := ..." << endl;
    }

    ii nnzCells = 0;
    if (a.is1_)
    {
        vector<ii> nnzs(a.m_, 0);
        for (ii i = 0; i < a.m_; i++)
        {
            for (ii a_nz = a.is0_[i]; a_nz < a.is1_[i]; a_nz++)
            {
                if (a.vs_[a_nz] > pruneThreshold)
                    nnzs[i]++;
            }
            nnzCells += nnzs[i];
        }
        
        if (nnzCells > 0)
        {
            ii* is0 = static_cast<ii*>(mkl_malloc(sizeof(ii) * (a.m_ + 1), 64));
            is0[0] = 0;
            ii* is1 = is0 + 1;
            for (ii i = 0; i < a.m_; i++)
                is1[i] = is0[i] + nnzs[i];

            ii* js = static_cast<ii*>(mkl_malloc(sizeof(ii) * is1[a.m_ - 1], 64));
            fp* vs = static_cast<fp*>(mkl_malloc(sizeof(fp) * is1[a.m_ - 1], 64));
            #pragma omp parallel
            for (ii i = 0; i < a.m_; i++)
            {
                ii nz = is0[i];
                for (ii a_nz = a.is0_[i]; a_nz < a.is1_[i]; a_nz++)
                {
                    if (a.vs_[a_nz] > pruneThreshold)
                    {
                        js[nz] = a.js_[a_nz];
                        vs[nz] = a.vs_[a_nz];
                        nz++;
                    }
                }
            }
            
            clear();
            m_ = a.m_;
            n_ = a.n_;
            is0_ = is0;
            is1_ = is1;
            js_ = js;
            vs_ = vs;
            status_ = mkl_sparse_s_create_csr(&mat_, SPARSE_INDEX_BASE_ZERO, m_, n_, is0_, is1_, js_, vs_); assert(!status_);
            isOwned_ = true;
        }
        else
        {
            init(a.m_, a.n_);
        }
    }
    else
    {
        init(a.m_, a.n_);
   }

    if (getDebugLevel() % 10 >= 4)
         cout << getTimeStamp() << "       ... X" << *this << endl;

    return nnzCells;
}


// MOSTLY OPTIMAL, MIGHT BE BETTER IF OMP LOOP DIDN'T INCLUDE PRUNED ROWS
ii MatrixSparse::pruneRows(const MatrixSparse& a, ii aNnzRows, const MatrixSparse& b, bool bRows, fp threshold)
{
    if (getDebugLevel() % 10 >= 4)
    {
        cout << getTimeStamp() << "       pruneRows(" << a << ",";
        cout << (bRows ? "rows(" : "cols(") << b << ")) where ";
        cout << (bRows ? "nRows" : "nCols") << " to prune > " << fixed << setprecision(1) << threshold * 100.0 << "% := ..." << endl;
    }

    if (b.is1_ && a.is1_)
    {
        vector<ii> rowOrColNnzs(a.m_, 0);
        if (bRows)
        {
            for (ii i = 0; i < b.m(); i++)
                rowOrColNnzs[i] += b.is1_[i] - b.is0_[i];
        }
        else
        {
            for (ii nz = 0; nz < b.nnz(); nz++)
                rowOrColNnzs[b.js_[nz]]++;
        }
        
        ii b_nnzRowsOrCols = 0;
        for (ii i = 0; i < a.m_; i++)
            if (rowOrColNnzs[i] > 0) b_nnzRowsOrCols++;

        if (b_nnzRowsOrCols / (fp) aNnzRows <= threshold)
        {
            ii* is0 = static_cast<ii*>(mkl_malloc(sizeof(ii) * (a.m_ + 1), 64));
            is0[0] = 0;
            ii* is1 = is0 + 1;
            for (ii i = 0; i < a.m_; i++)
                is1[i] = is0[i] + (rowOrColNnzs[i] > 0 ? a.is1_[i] - a.is0_[i] : 0);

            ii* js = static_cast<ii*>(mkl_malloc(sizeof(ii) * is1[a.m_ - 1], 64));
            fp* vs = static_cast<fp*>(mkl_malloc(sizeof(fp) * is1[a.m_ - 1], 64));
            #pragma omp parallel
            for (ii i = 0; i < a.m_; i++)
            {
                memcpy(&js[is0[i]], &a.js_[a.is0_[i]], sizeof(ii) * (is1[i] - is0[i]));
                memcpy(&vs[is0[i]], &a.vs_[a.is0_[i]], sizeof(fp) * (is1[i] - is0[i]));
            }
            
            clear();
            m_ = a.m_;
            n_ = a.n_;
            is0_ = is0;
            is1_ = is1;
            js_ = js;
            vs_ = vs;
            status_ = mkl_sparse_s_create_csr(&mat_, SPARSE_INDEX_BASE_ZERO, m_, n_, is0_, is1_, js_, vs_); assert(!status_);
            isOwned_ = true;
            
            aNnzRows = b_nnzRowsOrCols;
        }
    }
    else
    {
        clear();
        m_ = a.m_;
        n_ = a.n_;
        
        aNnzRows = 0;
    }
    
    if (getDebugLevel() % 10 >= 4)
        cout << getTimeStamp() << "       ... X" << *this << " (" << aNnzRows << " " << (bRows ? "rows" : "cols") << " pruned)" << endl;

    return aNnzRows;
}


void MatrixSparse::add(fp alpha, bool transposeA, const MatrixSparse& a, const MatrixSparse& b)
{
    init();
    
    if (getDebugLevel() % 10 >= 4)
    {
        cout << getTimeStamp() << "       " << alpha << " * " << (transposeA ? "t(" : "") << "A" << a << (transposeA ? ")" : "") << " + B" << b;
        cout << " := ..." << endl;
    }
    
    assert((transposeA ? a.n() : a.m()) == b.m());
    assert((transposeA ? a.m() : a.n()) == b.n());
    
    if (!a.is1_ || !b.is1_)
    {
        m_ = transposeA ? a.n() : a.m();
        n_ = b.n();
    }
    else if (!a.is1_)
    {
        copy(b);
    }
    else if (!b.is1_)
    {
        copy(a, transposeA);
        mul(alpha);
    }
    else
    {
        status_ = mkl_sparse_s_add(transposeA ? SPARSE_OPERATION_TRANSPOSE : SPARSE_OPERATION_NON_TRANSPOSE, a.mat_, alpha, b.mat_, &mat_); assert(!status_);
        
        sparse_index_base_t indexing;
        status_ = mkl_sparse_s_export_csr(mat_, &indexing, &m_, &n_, &is0_, &is1_, &js_, &vs_); assert(!status_);
        isOwned_ = false;
    }
    
    if (getDebugLevel() % 10 >= 4)
        cout << getTimeStamp() << "       ... X" << *this << endl;
 }


void MatrixSparse::matmul(bool transposeA, const MatrixSparse& a, const MatrixSparse& b, bool accumulate, bool denseOutput)
{
    if (!accumulate) init();
    
    if (getDebugLevel() % 10 >= 4)
    {
        cout << getTimeStamp() << "       " << (transposeA ? "t(" : "") << "A" << a << (transposeA ? ")" : "") << " %*% B" << b;
        if (accumulate) cout << " + X" << *this;
        cout << " := ..." << endl;
    }
    
    assert((transposeA ? a.m() : a.n()) == b.m());
    
    if (!is1_)
        accumulate = false;

    if (!a.is1_ || !b.is1_)
    {
        m_ = transposeA ? a.n() : a.m();
        n_ = b.n();
    }
    else
    {
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
            status_ = mkl_sparse_s_create_csr(&t, SPARSE_INDEX_BASE_ZERO, m, n, is0, is1, js, vs); assert(!status_);

            if (accumulate)
            {
                sparse_matrix_t y;
                status_ = mkl_sparse_s_add(SPARSE_OPERATION_NON_TRANSPOSE, mat_, 1.0, t, &y); assert(!status_);
                status_ = mkl_sparse_destroy(t); assert(!status_);

                init();
                mat_ = y;
                
                sparse_index_base_t indexing;
                status_ = mkl_sparse_s_export_csr(mat_, &indexing, &m_, &n_, &is0_, &is1_, &js_, &vs_); assert(!status_);
                isOwned_ = false;
            }
            else
            {
                m_ = m;
                n_ = n;
                is0_ = is0;
                is1_ = is1;
                js_ = js;
                vs_ = vs;
                mat_ = t;
                isOwned_ = true;
            }
        }
        else
        {
            if (accumulate)
            {
                sparse_matrix_t t;
                status_ = mkl_sparse_spmm(transposeA ? SPARSE_OPERATION_TRANSPOSE : SPARSE_OPERATION_NON_TRANSPOSE, a.mat_, b.mat_, &t); assert(!status_);
                
                sparse_matrix_t y;
                status_ = mkl_sparse_s_add(SPARSE_OPERATION_NON_TRANSPOSE, mat_, 1.0, t, &y); assert(!status_);
                status_ = mkl_sparse_destroy(t); assert(!status_);

                init();
                mat_ = y;
                
                sparse_index_base_t indexing;
                status_ = mkl_sparse_s_export_csr(mat_, &indexing, &m_, &n_, &is0_, &is1_, &js_, &vs_); assert(!status_);
                isOwned_ = false;
            }
            else
            {
                status_ = mkl_sparse_spmm(transposeA ? SPARSE_OPERATION_TRANSPOSE : SPARSE_OPERATION_NON_TRANSPOSE, a.mat_, b.mat_, &mat_); assert(!status_);
                
                sparse_index_base_t indexing;
                status_ = mkl_sparse_s_export_csr(mat_, &indexing, &m_, &n_, &is0_, &is1_, &js_, &vs_); assert(!status_);
                isOwned_ = false;
            }
        }
    }

    if (getDebugLevel() % 10 >= 4)
    {
        cout << getTimeStamp() << "       ... X" << *this;
        if (denseOutput) cout << " (DENSE)";
        cout << endl;
    }
}


// WARNING: IPP IS NOT 64BIT LENGTH
void MatrixSparse::mul(fp beta)
{
    if (getDebugLevel() % 10 >= 4)
    {
        cout << getTimeStamp() << "       X" << *this << " * ";
        cout.unsetf(ios::floatfield);
        cout << setprecision(8) << beta << " := ..." << endl;
    }

    if (is1_)
        ippsMulC_32f_I(beta, vs_, is1_[m_ - 1]);

    if (getDebugLevel() % 10 >= 4)
        cout << getTimeStamp() << "       ... X" << *this << endl;
}


void MatrixSparse::mul(const fp* a_vs)
{
    if (getDebugLevel() % 10 >= 4)
        cout << getTimeStamp() << "       X" << *this << " * A := ..." << endl;

    if (is1_)
        vsMul(is1_[m_ - 1], vs_, a_vs, vs_);

    if (getDebugLevel() % 10 >= 4)
        cout << getTimeStamp() << "       ... X" << *this << endl;
 }


void MatrixSparse::sqr()
{
    if (getDebugLevel() % 10 >= 4)
        cout << getTimeStamp() << "       (X" << *this << ")^2 := ..." << endl;

    if (is1_)
        vsSqr(is1_[m_ - 1], vs_, vs_);

    if (getDebugLevel() % 10 >= 4)
        cout << getTimeStamp() << "       ... X" << *this << endl;
}


void MatrixSparse::sqrt()
{
    if (getDebugLevel() % 10 >= 4)
        cout << getTimeStamp() << "       sqrt(X" << *this << ") := ..." << endl;

    if (is1_)
        vsSqrt(is1_[m_ - 1], vs_, vs_);

    if (getDebugLevel() % 10 >= 4)
        cout << getTimeStamp() << "       ... X" << *this << endl;
 }


void MatrixSparse::pow(fp power)
{
    if (getDebugLevel() % 10 >= 4)
        cout << getTimeStamp() << "       X" << *this << "^" << power << " := ..." << endl;

    if (is1_)
        vsPowx(is1_[m_ - 1], vs_, power, vs_);

    if (getDebugLevel() % 10 >= 4)
        cout << getTimeStamp() << "       ... X" << *this << endl;
 }


// WARNING: IPP IS NOT 64BIT LENGTH
void MatrixSparse::addNonzeros(fp beta)
{
    if (getDebugLevel() % 10 >= 4)
    {
        cout << getTimeStamp() << "       X" << *this << " + ";
        cout.unsetf(ios::floatfield);
        cout << setprecision(8) << beta << " := ..." << endl;
    }

    if (is1_)
        ippsAddC_32f_I(beta, vs_, is1_[m_ - 1]);

    if (getDebugLevel() % 10 >= 4)
        cout << getTimeStamp() << "       ... X" << *this << endl;
}


void MatrixSparse::addNonzeros(const fp* a_vs)
{
    if (getDebugLevel() % 10 >= 4)
        cout << getTimeStamp() << "       X" << *this << " / A := ..." << endl;

    if (is1_)
        vsAdd(is1_[m_ - 1], vs_, a_vs, vs_);

    if (getDebugLevel() % 10 >= 4)
        cout << getTimeStamp() << "       ... X" << *this << endl;
}


void MatrixSparse::lnNonzeros()
{
    if (getDebugLevel() % 10 >= 4)
        cout << getTimeStamp() << "       ln(X" << *this << ") := ..." << endl;

    if (is1_)
    {
#ifdef NDEBUG
        vsLn(is1_[m_ - 1], vs_, vs_); // for some reason this causes valgrind to crash
#else
        for (ii i = 0; i < is1_[m_ - 1]; i++)
            vs_[i] = log(vs_[i]);
#endif
    }

    if (getDebugLevel() % 10 >= 4)
        cout << getTimeStamp() << "       ... X" << *this << endl;
}


void MatrixSparse::expNonzeros()
{
    if (getDebugLevel() % 10 >= 4)
        cout << getTimeStamp() << "       exp(X" << *this << ") := ..." << endl;

    if (is1_)
        vsExp(is1_[m_ - 1], vs_, vs_);

    if (getDebugLevel() % 10 >= 4)
        cout << getTimeStamp() << "       ... X" << *this << endl;
}

void MatrixSparse::divNonzeros(const fp* a_vs)
{
    if (getDebugLevel() % 10 >= 4)
        cout << getTimeStamp() << "       X" << *this << " / A := ..." << endl;

    if (is1_)
        vsDiv(is1_[m_ - 1], vs_, a_vs, vs_);

    if (getDebugLevel() % 10 >= 4)
         cout << getTimeStamp() << "       ... X" << *this << endl;
 }


void MatrixSparse::div2Nonzeros(const fp* a_vs)
{
    if (getDebugLevel() % 10 >= 4)
         cout << getTimeStamp() << "       X" << *this << " / A := ..." << endl;

    if (is1_)
        vsDiv(is1_[m_ - 1], a_vs, vs_, vs_);

    if (getDebugLevel() % 10 >= 4)
        cout << getTimeStamp() << "       ... X" << *this << endl;
 }


fp MatrixSparse::sum() const
{
    if (getDebugLevel() % 10 >= 4)
        cout << getTimeStamp() << "       sum(X" << *this << ") := ..." << endl;

    fp sum = 0.0;
    if (is1_)
        ippsSum_32f(vs_, is1_[m_ - 1], &sum, ippAlgHintFast);

    if (getDebugLevel() % 10 >= 4)
    {
        cout << getTimeStamp() << "       ... ";
        cout.unsetf(ios::floatfield);
        cout << setprecision(8) << sum << endl;
    }

    return sum;
}


// WARNING: IPP IS NOT 32BIT LENGTH
fp MatrixSparse::sumSqrs() const
{
    if (getDebugLevel() % 10 >= 4)
        cout << getTimeStamp() << "       sum((X" << *this << ")^2) := ..." << endl;

    fp sum = 0.0;
    if (is1_)
    {
        ippsNorm_L2_32f(vs_, is1_[m_ - 1], &sum);
        sum *= sum;
    }

    if (getDebugLevel() % 10 >= 4)
    {
        cout << getTimeStamp() << "       ... ";
        cout.unsetf(ios::floatfield);
        cout << setprecision(8) << sum << endl;
    }

    return sum;
}


// WARNING: IPP IS NOT 32BIT LENGTH
fp MatrixSparse::sumSqrDiffsNonzeros(const fp* a_vs) const
{
    if (getDebugLevel() % 10 >= 4)
        cout << getTimeStamp() << "       sum((X" << *this << " - A)^2) := ..." << endl;

    fp sum = 0.0;
    if (is1_)
    {
        ippsNormDiff_L2_32f(vs_, a_vs, is1_[m_ - 1], &sum);
        sum *= sum;
    }
    
    if (getDebugLevel() % 10 >= 4)
    {
        cout << getTimeStamp() << "       ... ";
        cout.unsetf(ios::floatfield);
        cout << setprecision(8) << sum << endl;
    }
    
    return sum;
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
        
        if (getDebugLevel() % 10 >= 4)
        {
            os << "(" << a.nnzActual() << ")";
        }
        
        os << "/" << a.size() << ":";
        os.unsetf(ios::floatfield);
        os << setprecision(3) << 100.0 * a.nnz() / (double)a.size() << "%";
    }

    return  os;
}
