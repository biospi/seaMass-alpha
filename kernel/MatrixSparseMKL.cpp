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


#include "MatrixSparseMKL.hpp"

#include <iomanip>
//#include <iostream>
#include <cassert>
//#include <cstring>
//#include <cmath>
//#include <sstream>
//#include <algorithm>

#if defined(_OPENMP)
  #include <omp.h>
#endif

#include <ippcore.h>
#include <ipps.h>
#include <ippi.h>

#include "../kernel/NetcdfFile.hpp"


using namespace std;


MatrixSparseMKL::MatrixSparseMKL() : m_(0), n_(0), is1_(0)
{
}


MatrixSparseMKL::~MatrixSparseMKL()
{
	free();
}


// SEEMS OPTIMAL
void MatrixSparseMKL::init(ii m, ii n, ii nnz, const fp* acoo, const ii* rowind, const ii* colind)
{
    free();
    
    if (getDebugLevel() % 10 >= 4)
    {
        cout << getTimeStamp() << "     COO := ..." << endl;
    }
    
    m_ = m;
    n_ = n;
    
    fp* c_acoo = const_cast<fp*>(acoo);
    ii* c_rowind = const_cast<ii*>(rowind);
    ii* c_colind = const_cast<ii*>(colind);
    
    sparse_matrix_t t;
    status_ = mkl_sparse_s_create_coo(&t, SPARSE_INDEX_BASE_ZERO, m_, n_, nnz, c_rowind, c_colind, c_acoo); assert(!status_);
    status_ = mkl_sparse_convert_csr(t, SPARSE_OPERATION_NON_TRANSPOSE, &mat_); assert(!status_);
    status_ = mkl_sparse_destroy(t); assert(!status_);
    
    sparse_index_base_t indexing;
    status_ = mkl_sparse_s_export_csr(mat_, &indexing, &m_, &n_, &is0_, &is1_, &js_, &vs_); assert(!status_);
    isOwned_ = false;
    
    if (getDebugLevel() % 10 >= 4)
    {
        cout << getTimeStamp() << "     ... X" << *this << endl;
        
        if (getDebugLevel() >= 14)
        {
            ostringstream oss; oss << setw(9) << setfill('0') << getId() << ".csr";
            write(oss.str());
        }
    }
}


void MatrixSparseMKL::init(ii m, ii n, fp v)
{
    free();
    
    vector<fp> acoo;
    vector<ii> rowind;
    vector<ii> colind;
    
    for (ii i = 0; i < m; i++)
    {
        for (ii j = 0; j < n; j++)
        {
            acoo.push_back(v);
            rowind.push_back(i);
            colind.push_back(j);
        }
    }
    
    init(m, n, (ii)acoo.size(), acoo.data(), rowind.data(), colind.data());
}


void MatrixSparseMKL::init(const std::vector<MatrixSparseMKL>& as)
{
    free();
    
    for (size_t k = 0; k < as.size(); k++) assert(as[k].m_ == 1);
    for (size_t k = 0; k < as.size() - 1; k++) assert(as[k].n_ == as[k + 1].n_);
    
    m_ = (ii)as.size();
    n_ = as[0].n();
    
    ii nnz = 0;
    for (ii i = 0; i < m_; i++)
    {
        nnz += as[i].nnz();
    }
    
    if (nnz > 0)
    {
        is0_ = static_cast<ii*>(mkl_malloc(sizeof(ii) * (m_ + 1), 64));
        is0_[0] = 0;
        is1_ = is0_ + 1;
        for (ii i = 0; i < m_; i++)
        {
            is1_[i] = is0_[i] + as[i].is1_[0];
        }
        
        js_ = static_cast<ii*>(mkl_malloc(sizeof(ii) * is1_[m_ - 1], 64));
        vs_ = static_cast<fp*>(mkl_malloc(sizeof(fp) * is1_[m_ - 1], 64));
        #pragma omp parallel
        for (ii i = 0; i < m_; i++)
        {
            memcpy(&js_[is0_[i]], as[i].js_, sizeof(ii) * as[i].is1_[0]);
            memcpy(&vs_[is0_[i]], as[i].vs_, sizeof(fp) * as[i].is1_[0]);
        }
        
        status_ = mkl_sparse_s_create_csr(&mat_, SPARSE_INDEX_BASE_ZERO, m_, n_, is0_, is1_, js_, vs_); assert(!status_);
        isOwned_ = true;
    }
}


void MatrixSparseMKL::init(const MatrixMKL& a)
{
    free();
    
    vector<fp> acoo;
    vector<ii> rowind;
    vector<ii> colind;
    
    for (ii i = 0; i < a.m(); i++)
    {
        for (ii j = 0; j < a.n(); j++)
        {
            if (a.vs()[j + i * a.n()] != 0.0)
            {
                acoo.push_back(a.vs()[j + i * a.n()]);
                rowind.push_back(i);
                colind.push_back(j);
            }
        }
    }
    
    init(a.m(), a.n(), (ii)acoo.size(), acoo.data(), rowind.data(), colind.data());
}


// we are not deallocating is0_ at the moment...
void MatrixSparseMKL::wrap(const MatrixSparseMKL& a, ii row)
{
    free();
    
    if (getDebugLevel() % 10 >= 4)
    {
        cout << getTimeStamp() << "     A" << a << "[" << row << "] := ..." << endl;
    }
    
    m_ = 1;
    n_ = a.n();
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
    {
        cout << getTimeStamp() << "     ... X" << *this << endl;
        
        if (getDebugLevel() >= 14)
        {
            ostringstream oss; oss << setw(9) << setfill('0') << getId() << ".csr";
            write(oss.str());
        }
    }
}


void MatrixSparseMKL::clear()
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


void MatrixSparseMKL::free()
{
    clear();
    m_ = 0;
    n_ = 0;
}


ii MatrixSparseMKL::m() const
{
    return m_;
}


ii MatrixSparseMKL::n() const
{
    return n_;
}


li MatrixSparseMKL::size() const
{
    return (li)m_ * n_;
}


ii MatrixSparseMKL::nnz() const
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


ii MatrixSparseMKL::nnzActual() const
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


ii* MatrixSparseMKL::js() const
{
    return js_;
}


fp* MatrixSparseMKL::vs() const
{
    return vs_;
}


// SEEMS OPTIMAL
void MatrixSparseMKL::copy(const MatrixSparseMKL& a, bool transpose)
{
    free();
    
    if (getDebugLevel() % 10 >= 4)
    {
        cout << getTimeStamp() << "     " << (transpose ? "t(" : "") << "A" << a << (transpose ? ")" : "") << " := ..." << endl;
    }
    
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

    if (getDebugLevel() % 10 >= 4)
    {
        cout << getTimeStamp() << "     ... X" << *this << endl;
        
        if (getDebugLevel() >= 14)
        {
            ostringstream oss; oss << setw(9) << setfill('0') << getId() << ".csr";
            write(oss.str());
        }
    }
}


class MatrixSparseMKL::MyComparator
{
public:
    
    ii* js_;
    
    MyComparator(ii* js): js_(js) {}
    
    bool operator()(ii i1, ii i2)
    {
        return js_[i1] < js_[i2];
    }
};


// SHOULD BE IN-PLACE SORTING WITHOUT CREATION OF INDICIES (LOOK AT IPP SORT)
void MatrixSparseMKL::sort()
{
    if (getDebugLevel() % 10 >= 4)
    {
        cout << getTimeStamp() << "     sort(X" << *this << ") := ..." << endl;
    }
    
    if (is1_)
    {
        //#pragma omp parallel // why does omp make it crash here?
        for (ii i = 0; i < m_; i++)
        {
            vector<ii> indicies(is1_[i] - is0_[i]);
            for (ii nz = 0; nz < indicies.size(); nz++) indicies[nz] = nz;
            std::sort(indicies.begin(), indicies.end(), MyComparator(&js_[is0_[i]]));
            
            {
                vector<ii> js(is1_[i] - is0_[i]);
                for (ii nz = 0; nz < js.size(); nz++) js[nz] = js_[is0_[i] + indicies[nz]];
                for (ii nz = 0; nz < js.size(); nz++) js_[is0_[i] + nz] = js[nz];
            }
            
            {
                vector<fp> vs(is1_[i] - is0_[i]);
                for (ii nz = 0; nz < vs.size(); nz++) vs[nz] = vs_[is0_[i] + indicies[nz]];
                for (ii nz = 0; nz < vs.size(); nz++) vs_[is0_[i] + nz] = vs[nz];
            }
        }
    }
    
    if (getDebugLevel() % 10 >= 4)
    {
        cout << getTimeStamp() << "     ... X" << *this << endl;
        
        if (getDebugLevel() >= 14)
        {
            ostringstream oss; oss << setw(9) << setfill('0') << getId() << ".csr";
            write(oss.str());
        }
    }
}


ii MatrixSparseMKL::prune(const MatrixSparseMKL& a, fp pruneThreshold)
{
    if (getDebugLevel() % 10 >= 4)
    {
        cout << getTimeStamp() << "     (A" << a << " > ";
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
                if (a.vs_[a_nz] > pruneThreshold) nnzs[i]++;
            }
            nnzCells += nnzs[i];
        }
        
        if (nnzCells > 0)
        {
            ii* is0 = static_cast<ii*>(mkl_malloc(sizeof(ii) * (a.m_ + 1), 64));
            is0[0] = 0;
            ii* is1 = is0 + 1;
            for (ii i = 0; i < a.m_; i++)
            {
                is1[i] = is0[i] + nnzs[i];
            }
            
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
            clear();
            m_ = a.m_;
            n_ = a.n_;
        }
    }
    else
    {
        clear();
        m_ = a.m_;
        n_ = a.n_;
    }

    if (getDebugLevel() % 10 >= 4)
    {
        cout << getTimeStamp() << "     ... X" << *this << endl;
        
        if (getDebugLevel() >= 14)
        {
            ostringstream oss; oss << setw(9) << setfill('0') << getId() << ".csr";
            write(oss.str());
        }
    }
    
    return nnzCells;
}


// MOSTLY OPTIMAL, MIGHT BE BETTER IF OMP LOOP DIDN'T INCLUDE PRUNED ROWS
ii MatrixSparseMKL::pruneRows(const MatrixSparseMKL& a, ii aNnzRows, const MatrixSparseMKL& b, bool bRows, fp threshold)
{
    if (getDebugLevel() % 10 >= 4)
    {
        cout << getTimeStamp() << "     pruneRows(" << a << ",";
        cout << (bRows ? "rows(" : "cols(") << b << ")) where ";
        cout << (bRows ? "nRows" : "nCols") << " to prune > " << fixed << setprecision(1) << threshold * 100.0 << "% := ..." << endl;
    }

    if (b.is1_ && a.is1_)
    {
        vector<ii> rowOrColNnzs(a.m_, 0);
        if (bRows)
        {
            for (ii i = 0; i < b.m(); i++)
            {
                rowOrColNnzs[i] += b.is1_[i] - b.is0_[i];
            }
        }
        else
        {
            for (ii nz = 0; nz < b.nnz(); nz++)
            {
                rowOrColNnzs[b.js_[nz]]++;
            }
        }
        
        ii b_nnzRowsOrCols = 0;
        for (ii i = 0; i < a.m_; i++)
        {
            if (rowOrColNnzs[i] > 0) b_nnzRowsOrCols++;
        }
        
        if (b_nnzRowsOrCols / (fp) aNnzRows <= threshold)
        {
            ii* is0 = static_cast<ii*>(mkl_malloc(sizeof(ii) * (a.m_ + 1), 64));
            is0[0] = 0;
            ii* is1 = is0 + 1;
            for (ii i = 0; i < a.m_; i++)
            {
                is1[i] = is0[i] + (rowOrColNnzs[i] > 0 ? a.is1_[i] - a.is0_[i] : 0);
            }
            
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
    {
        cout << getTimeStamp() << "     ... X" << *this << " (" << aNnzRows << " " << (bRows ? "rows" : "cols") << " pruned)" << endl;
        
        if (getDebugLevel() >= 14)
        {
            ostringstream oss; oss << setw(9) << setfill('0') << getId() << ".csr";
            write(oss.str());
        }
    }
    
    return aNnzRows;
}



// USE MKL FUNCTION THATS CONVERTS TO DENSE
void MatrixSparseMKL::output(fp* vs) const
{
    if (getDebugLevel() % 10 >= 4)
    {
        cout << getTimeStamp() << "     X" << *this << " := ..." << endl;
    }
    
    for (li x = 0; x < size(); x++) vs[x] = 0.0;
    
    ii i = 0;
    for (ii nz = 0; nz < nnz(); nz++)
    {
        while (nz >= is1_[i]) i++; // row of nz'th non-zero
        ii j = js_[nz]; // column of nz'th non-zero
        
        vs[j + i * n_] = vs_[nz];
    }
    
    if (getDebugLevel() % 10 >= 4)
    {
        cout << getTimeStamp() << "     ... out" << endl;
    }
    
    /*ii nnz = 0;
    for (ii j = 0; j < n_; j++)
    {
        for (ii i = 0; i < m_; i++)
        {
            cout << fixed << setprecision(1) << vs[j + i * n_] << ",";
            if (vs[j + i * n_] > 0.0) nnz++;
        }
        cout << endl ;
    }
    cout << nnz << endl;*/
}


void MatrixSparseMKL::write(const string& filename) const
{
    if (m_ > 0)
    {
        NetCDFile outFile(filename, NC_NETCDF4);
        if (is1_)
        {
            outFile.write_VecNC("ia", is0_, m_ + 1, NC_INT64);
            outFile.write_VecNC("ja", js_, nnz(), NC_INT64);
            outFile.write_VecNC("a", vs_, nnz(), NC_FLOAT);
        }
        else
        {
            vector<ii> is(m_ + 1, 0);
            outFile.write_VecNC("ia", is, NC_INT64);
        }
        vector<ii> n;
        n.push_back(n_);
        outFile.write_AttNC("ia", "n", n, NC_INT64);
    }
}


void MatrixSparseMKL::add(fp alpha, bool transposeA, const MatrixSparseMKL& a, const MatrixSparseMKL& b)
{
    free();
    
    if (getDebugLevel() % 10 >= 4)
    {
        cout << getTimeStamp() << "     " << alpha << " * " << (transposeA ? "t(" : "") << "A" << a << (transposeA ? ")" : "") << " + B" << b;
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
    {
        cout << getTimeStamp() << "     ... X" << *this << endl;
        
        if (getDebugLevel() >= 14)
        {
            ostringstream oss; oss << setw(9) << setfill('0') << getId() << ".csr";
            write(oss.str());
        }
    }
}


void MatrixSparseMKL::matmul(bool transposeA, const MatrixSparseMKL& a, const MatrixSparseMKL& b, bool accumulate, bool denseOutput)
{
    if (!accumulate) free();
    
    if (getDebugLevel() % 10 >= 4)
    {
        cout << getTimeStamp() << "     " << (transposeA ? "t(" : "") << "A" << a << (transposeA ? ")" : "") << " %*% B" << b;
        if (accumulate) cout << " + X" << *this;
        cout << " := ..." << endl;
    }
    
    assert((transposeA ? a.m() : a.n()) == b.m());
    
    if (!is1_)
    {
        accumulate = false;
    }
    
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
            {
                is1[i] = is0[i] + n;
            }
            
            ii* js = static_cast<ii*>(mkl_malloc(sizeof(ii) * is1[m - 1], 64));
            for (ii nz = 0; nz < is1[m - 1]; nz++)
            {
                js[nz] = nz % n;
            }
            
            fp* vs = static_cast<fp*>(mkl_malloc(sizeof(fp) * is1[m - 1], 64));
            status_ = mkl_sparse_s_spmmd(transposeA ? SPARSE_OPERATION_TRANSPOSE : SPARSE_OPERATION_NON_TRANSPOSE, a.mat_, b.mat_, SPARSE_LAYOUT_ROW_MAJOR, vs, n); assert(!status_);

            sparse_matrix_t t;
            status_ = mkl_sparse_s_create_csr(&t, SPARSE_INDEX_BASE_ZERO, m, n, is0, is1, js, vs); assert(!status_);

            if (accumulate)
            {
                sparse_matrix_t y;
                status_ = mkl_sparse_s_add(SPARSE_OPERATION_NON_TRANSPOSE, mat_, 1.0, t, &y); assert(!status_);
                status_ = mkl_sparse_destroy(t); assert(!status_);
                
                free();
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
                
                free();
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
        cout << getTimeStamp() << "     ... X" << *this;
        if (denseOutput) cout << " (DENSE)";
        cout << endl;
        
        if (getDebugLevel() >= 14)
        {
            ostringstream oss; oss << setw(9) << setfill('0') << getId() << ".csr";
            write(oss.str());
        }
    }
}


// WARNING: IPP IS NOT 32BIT LENGTH
void MatrixSparseMKL::mul(fp beta)
{
    if (getDebugLevel() % 10 >= 4)
    {
        cout << getTimeStamp() << "     X" << *this << " * ";
        cout.unsetf(ios::floatfield);
        cout << setprecision(8) << beta << " := ..." << endl;
    }
    
    ippsMulC_32f_I(beta, vs_, nnz());
    
    if (getDebugLevel() % 10 >= 4)
    {
        cout << getTimeStamp() << "     ... X" << *this << endl;
        
        if (getDebugLevel() >= 14)
        {
            ostringstream oss; oss << setw(9) << setfill('0') << getId() << ".csr";
            write(oss.str());
        }
    }
}


void MatrixSparseMKL::mul(const MatrixSparseMKL& a)
{
    if (getDebugLevel() % 10 >= 4)
    {
        cout << getTimeStamp() << "     X" << *this << " * A" << a << " := ..." << endl;
    }
    
    assert(m_ == a.m_);
    assert(n_ == a.n_);
    assert(nnz() == a.nnz());
    for (ii nz = 0; nz < nnz(); nz++) assert(js_[nz] == a.js_[nz]);
    
    vsMul(nnz(), vs_, a.vs_, vs_);
    
    if (getDebugLevel() % 10 >= 4)
    {
        cout << getTimeStamp() << "     ... X" << *this << endl;
        
        if (getDebugLevel()   >= 14)
        {
            ostringstream oss; oss << setw(9) << setfill('0') << getId() << ".csr";
            write(oss.str());
        }
    }
}


void MatrixSparseMKL::sqr()
{
    if (getDebugLevel() % 10 >= 4)
    {
        cout << getTimeStamp() << "     (X" << *this << ")^2 := ..." << endl;
    }
    
    vsSqr(nnz(), vs_, vs_);
    
    if (getDebugLevel() % 10 >= 4)
    {
        cout << getTimeStamp() << "     ... X" << *this << endl;
        
        if (getDebugLevel() >= 14)
        {
            ostringstream oss; oss << setw(9) << setfill('0') << getId() << ".csr";
            write(oss.str());
        }
    }
}


void MatrixSparseMKL::sqrt()
{
    if (getDebugLevel() % 10 >= 4)
    {
        cout << getTimeStamp() << "     sqrt(X" << *this << ") := ..." << endl;
    }
    
    vsSqrt(nnz(), vs_, vs_);
    
    if (getDebugLevel() % 10 >= 4)
    {
        cout << getTimeStamp() << "     ... X" << *this << endl;
        
        if (getDebugLevel() >= 14)
        {
            ostringstream oss; oss << setw(9) << setfill('0') << getId() << ".csr";
            write(oss.str());
        }
    }
}


// WARNING: IPP IS NOT 32BIT LENGTH
void MatrixSparseMKL::setNonzeros(fp v)
{
    if (getDebugLevel() % 10 >= 4)
    {
        cout << getTimeStamp() << "     (X" << *this << " != 0.0) ? " << v << " : 0.0 := ..." << endl;
    }
    
    ippsSet_32f(v, vs_, nnz());
    
    if (getDebugLevel() % 10 >= 4)
    {
        cout << getTimeStamp() << "     ... X" << *this << endl;
        
        if (getDebugLevel() >= 14)
        {
            ostringstream oss; oss << setw(9) << setfill('0') << getId() << ".csr";
            write(oss.str());
        }
    }
}


// WARNING: IPP IS NOT 32BIT LENGTH
void MatrixSparseMKL::addNonzeros(fp beta)
{
    if (getDebugLevel() % 10 >= 4)
    {
        cout << getTimeStamp() << "     X" << *this << " + ";
        cout.unsetf(ios::floatfield);
        cout << setprecision(8) << beta << " := ..." << endl;
    }
    
    ippsAddC_32f_I(beta, vs_, nnz());
    
    if (getDebugLevel() % 10 >= 4)
    {
        cout << getTimeStamp() << "     ... X" << *this << endl;
        
        if (getDebugLevel() >= 14)
        {
            ostringstream oss; oss << setw(9) << setfill('0') << getId() << ".csr";
            write(oss.str());
        }
    }
}


void MatrixSparseMKL::lnNonzeros()
{
    if (getDebugLevel() % 10 >= 4)
    {
        cout << getTimeStamp() << "     ln(X" << *this << ") := ..." << endl;
    }
    
    vsLn(nnz(), vs_, vs_);
    
    if (getDebugLevel() % 10 >= 4)
    {
        cout << getTimeStamp() << "     ... X" << *this << endl;
        
        if (getDebugLevel() >= 14)
        {
            ostringstream oss; oss << setw(9) << setfill('0') << getId() << ".csr";
            write(oss.str());
        }
    }
}


void MatrixSparseMKL::expNonzeros()
{
    if (getDebugLevel() % 10 >= 4)
    {
        cout << getTimeStamp() << "     exp(X" << *this << ") := ..." << endl;
    }
    
    vsExp(nnz(), vs_, vs_);
    
    if (getDebugLevel() % 10 >= 4)
    {
        cout << getTimeStamp() << "     ... X" << *this << endl;
        
        if (getDebugLevel() >= 14)
        {
            ostringstream oss; oss << setw(9) << setfill('0') << getId() << ".csr";
            write(oss.str());
        }
    }
}


void MatrixSparseMKL::pow(fp power)
{
    if (getDebugLevel() % 10 >= 4)
    {
        cout << getTimeStamp() << "     X" << *this << "^" << power << " := ..." << endl;
    }
    
    vsPowx(nnz(), vs_, power, vs_);
    
    if (getDebugLevel() % 10 >= 4)
    {
        cout << getTimeStamp() << "     ... X" << *this << endl;
        
        if (getDebugLevel() >= 14)
        {
            ostringstream oss; oss << setw(9) << setfill('0') << getId() << ".csr";
            write(oss.str());
        }
    }
}


fp MatrixSparseMKL::sum() const
{
    if (getDebugLevel() % 10 >= 4)
    {
        cout << getTimeStamp() << "     sum(X" << *this << ") := ..." << endl;
    }
    
    fp sum = 0.0;
    ippsSum_32f(vs_, nnz(), &sum, ippAlgHintFast);
    
    if (getDebugLevel() % 10 >= 4)
    {
        cout << getTimeStamp() << "     ... ";
        cout.unsetf(ios::floatfield);
        cout << setprecision(8) << sum << endl;
    }

	return sum;
}


void MatrixSparseMKL::divNonzeros(const fp* a_vs)
{
    if (getDebugLevel() % 10 >= 4)
    {
        cout << getTimeStamp() << "     X" << *this << " / A := ..." << endl;
    }
    
    vsDiv(nnz(), vs_, a_vs, vs_);
    
    if (getDebugLevel() % 10 >= 4)
    {
        cout << getTimeStamp() << "     ... X" << *this << endl;
        
        if (getDebugLevel() >= 14)
        {
            ostringstream oss; oss << setw(9) << setfill('0') << getId() << ".csr";
            write(oss.str());
        }
    }
}


void MatrixSparseMKL::div2Nonzeros(const fp* a_vs)
{
    if (getDebugLevel() % 10 >= 4)
    {
        cout << getTimeStamp() << "     X" << *this << " / A := ..." << endl;
    }
    
    vsDiv(nnz(), a_vs, vs_, vs_);
    
    if (getDebugLevel() % 10 >= 4)
    {
        cout << getTimeStamp() << "     ... X" << *this << endl;
        
        if (getDebugLevel() >= 14)
        {
            ostringstream oss; oss << setw(9) << setfill('0') << getId() << ".csr";
            write(oss.str());
        }
    }
}


// WARNING: IPP IS NOT 32BIT LENGTH
fp MatrixSparseMKL::sumSqrs() const
{
    if (getDebugLevel() % 10 >= 4)
    {
        cout << getTimeStamp() << "     sum((X" << *this << ")^2) := ..." << endl;
    }
    
    fp sum = 0.0;
    ippsNorm_L2_32f(vs_, nnz(), &sum);
    sum *= sum;

    if (getDebugLevel() % 10 >= 4)
    {
        cout << getTimeStamp() << "     ... ";
        cout.unsetf(ios::floatfield);
        cout << setprecision(8) << sum << endl;
    }

	return sum;
}


// WARNING: IPP IS NOT 32BIT LENGTH
fp MatrixSparseMKL::sumSqrDiffsNonzeros(const fp* a_vs) const
{
    if (getDebugLevel() % 10 >= 4)
    {
        cout << getTimeStamp() << "     sum((X" << *this << " - A)^2) := ..." << endl;
    }
    
    fp sum = 0.0;
    ippsNormDiff_L2_32f(vs_, a_vs, nnz(), &sum);
    sum *= sum;
    
    if (getDebugLevel() % 10 >= 4)
    {
        cout << getTimeStamp() << "     ... ";
        cout.unsetf(ios::floatfield);
        cout << setprecision(8) << sum << endl;
    }
    
    return sum;
}


// SEEMS OPTIMAL
void MatrixSparseMKL::subsetCopy(const MatrixSparseMKL& a)
{
    if (getDebugLevel() % 10 >= 4)
    {
        cout << getTimeStamp() << "     " << a << " within " << *this << " := ..." << endl;
    }
    
    assert(m_ == a.m_);
    assert(n_ == a.n_);
    
    if (is1_ && a.is1_)
    {
        for (ii i = 0; i < m_; i++) for (ii nz = is0_[i]; nz < is1_[i] - 1; nz++) assert(js_[nz] <= js_[nz + 1]);
        for (ii i = 0; i < a.m_; i++) for (ii a_nz = a.is0_[i]; a_nz < a.is1_[i] - 1; a_nz++) assert(a.js_[a_nz] <= a.js_[a_nz + 1]);
        
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
                    {
                        a_nz++;
                    }
                }
            }
        }
    }
    
    if (getDebugLevel() % 10 >= 4)
    {
        cout << getTimeStamp() << "     ... X" << *this << endl;
        
        if (getDebugLevel() >= 14)
        {
            ostringstream oss; oss << getId() << ".csr";
            write(oss.str());
        }
    }
}


ostream& operator<<(ostream& os, const MatrixSparseMKL& a)
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
