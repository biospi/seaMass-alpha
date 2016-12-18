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

#include <iomanip>
#include <iostream>
#include <cassert>
#include <cstring>

#if defined(_OPENMP)
  #include <omp.h>
#endif

#include <ippcore.h>
#include <ipps.h>


using namespace std;


void setNumThreads(int threads)
{
    cout << "Config: " << 8 * sizeof(ii) << "bit MKL addressing, ";

#if defined(_OPENMP)
    if (threads == 0)
    {
        threads = mkl_get_max_threads();
    }
    else
    {
        mkl_set_num_threads(threads);
    }

    ippSetNumThreads(threads);
    omp_set_num_threads(threads);
    
    cout << threads << " MKL/IPP threads, " << threads << " OpenMP threads" << endl;
#else
    cout << "non-OpenMP build (hence no thread control)" << endl;
#endif
}


double getWallTime()
{
    return dsecnd();
}


MatrixSparse::MatrixSparse() : nnz_(-1), data_(Data::NONE)
{
}


MatrixSparse::~MatrixSparse()
{
	free();
}


void MatrixSparse::free()
{
    if (data_ != Data::NONE)
    {
        mkl_sparse_destroy(*mat_);
        delete mat_;
        
        if (data_ == Data::EXTERNAL)
        {
            mkl_free(is_);
            mkl_free(js_);
            mkl_free(vs_);
        }
        else if (data_ == Data::EXTERNAL_MKL)
        {
            mkl_sparse_destroy(*mkl_);
            delete mkl_;
        }
        
        data_ = Data::NONE;
    }
}


bool MatrixSparse::operator!() const
{
    return data_ == Data::NONE;
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
    if (nnz_ == -1)
    {
        sparse_index_base_t indexing; ii m; ii n; ii* is0; ii* is1; ii* js; fp* vs;
        mkl_sparse_s_export_csr(*mat_, &indexing, &m, &n, &is0, &is1, &js, &vs);
    
        return is1[m - 1];
    }
    else
    {
        return nnz_;
    }
}


li MatrixSparse::mem() const
{
    return sizeof(*this) + sizeof(fp) * (li)nnz() + (li)sizeof(ii) * ((li)nnz() + m_ + 1);
}


void MatrixSparse::init(ii m, ii n)
{
#ifndef NDEBUG
    cout << "  0 := " << flush;
#endif
    
    free();
    
    m_ = m;
    n_ = n;
    nnz_ = 0; // HAVE TO SPECIFY THIS TO AVOID MKL PROCESSING EMPTY MATRICES
    is_ = static_cast<ii*>(mkl_calloc(m_ + 1, sizeof(ii), 64));
    js_ = static_cast<ii*>(mkl_malloc(sizeof(ii), 64)); // MKL doesn't like the pointer being null
    vs_ = static_cast<fp*>(mkl_malloc(sizeof(fp), 64)); // MKL doesn't like the pointer being null
    
    mat_ = new sparse_matrix_t;
    mkl_sparse_s_create_csr(mat_, SPARSE_INDEX_BASE_ZERO, m_, n_, is_, &is_[1], js_, vs_);
    data_ = Data::EXTERNAL;
    
#ifndef NDEBUG
    cout << "X" << *this << endl;
#endif
}


void MatrixSparse::init(ii m, ii n, fp v)
{
#ifndef NDEBUG
    cout << "  " << v << " := " << flush;
#endif
    
    free();
    
    m_ = m;
    n_ = n;
    is_ = static_cast<ii*>(mkl_malloc(sizeof(ii) * (m_ + 1), 64));
    for (ii i = 0; i <= m_; i++) is_[i] = i * n;
    js_ = static_cast<ii*>(mkl_malloc(sizeof(ii) * is_[m_], 64));
    for (ii i = 0; i < m; i++) for (ii j = 0; j < n; j++) js_[j + i * n] = j;
    vs_ = static_cast<fp*>(mkl_malloc(sizeof(fp) * is_[m_], 64));
    for (ii nz = 0; nz < is_[m_]; nz++) vs_[nz] = v;
    
    mat_ = new sparse_matrix_t;
    mkl_sparse_s_create_csr(mat_, SPARSE_INDEX_BASE_ZERO, m_, n_, is_, &is_[1], js_, vs_);
    data_ = Data::EXTERNAL;
    
#ifndef NDEBUG
    cout << "X" << *this << endl;
#endif
}


void MatrixSparse::init(ii m, ii n, const fp* vs)
{
    vector<fp> acoo;
    vector<ii> rowind;
    vector<ii> colind;
    for (ii i = 0; i < m; i++)
    {
        for (ii j = 0; j < n; j++)
        {
            if (vs[j + i*n] != 0.0)
            {
                acoo.push_back(vs[j + i*n]);
                rowind.push_back(i);
                colind.push_back(j);
            }
        }
    }
    
    init(m, n, (ii)acoo.size(), acoo.data(), rowind.data(), colind.data());
}


void MatrixSparse::init(ii m, ii n, ii nnz, const fp* acoo, const ii* rowind, const ii* colind)
{
#ifndef NDEBUG
    cout << "  COO := " << flush;
#endif
    
    free();
    
    m_ = m;
    n_ = n;
    fp* c_acoo = const_cast<fp*>(acoo);
    ii* c_rowind = const_cast<ii*>(rowind);
    ii* c_colind = const_cast<ii*>(colind);
    
    sparse_matrix_t t;
    mkl_sparse_s_create_coo(&t, SPARSE_INDEX_BASE_ZERO, m_, n_, nnz, c_rowind, c_colind, c_acoo);
    mat_ = new sparse_matrix_t;
    mkl_sparse_convert_csr(t, SPARSE_OPERATION_NON_TRANSPOSE, mat_);
    data_ = Data::MKL;
    
#ifndef NDEBUG
    cout << "X" << *this << endl;
#endif
}


void MatrixSparse::copy(const MatrixSparse& a, Transpose transpose)
{
#ifndef NDEBUG
    cout << "  A" << (transpose == Transpose::NO ? "" : "t") << a << " := " << flush;
#endif
    
    free();

    if (transpose == Transpose::NO)
    {
        if (a.nnz_ == 0)
        {
            init(a.m_, a.n_);
        }
        else
        {
            sparse_index_base_t a_indexing; ii a_m; ii a_n; ii* a_is0; ii* a_is1; ii* a_js; fp* a_vs;
            mkl_sparse_s_export_csr(*a.mat_, &a_indexing, &a_m, &a_n, &a_is0, &a_is1, &a_js, &a_vs);
            
            m_ = a.m_;
            n_ = a.n_;
            is_ = static_cast<ii*>(mkl_malloc(sizeof(ii) * (m_ + 1), 64));
            memcpy(is_, a_is0, sizeof(ii) * (m_ + 1));
            js_ = static_cast<ii*>(mkl_malloc(sizeof(ii) * is_[m_], 64));
            memcpy(js_, a_js, sizeof(ii) * is_[m_]);
            vs_ = static_cast<fp*>(mkl_malloc(sizeof(fp) * is_[m_], 64));
            memcpy(vs_, a_vs, sizeof(fp) * is_[m_]);
            
            mat_ = new sparse_matrix_t;
            mkl_sparse_s_create_csr(mat_, SPARSE_INDEX_BASE_ZERO, m_, n_, is_, &is_[1], js_, vs_);
            data_ = Data::EXTERNAL;
        }
    }
    else
    {
        if (a.nnz_ == 0)
        {
            init(a.n_, a.m_);
        }
        else
        {
            mat_ = new sparse_matrix_t;
            mkl_sparse_convert_csr(*a.mat_, SPARSE_OPERATION_TRANSPOSE, mat_);
            data_ = Data::MKL;
        }
    }

#ifndef NDEBUG
	cout << "X" << *this << endl;
#endif
}


// todo: optimise
void MatrixSparse::copy(const MatrixSparse& a, fp pruneThreshold)
{
#ifndef NDEBUG
	cout << "  (A" << a << " >= ";
	cout.unsetf(ios::floatfield);
	cout << setprecision(8) << pruneThreshold << ") ? A : 0.0 := " << flush;
#endif
    
    free();
    
    sparse_index_base_t a_indexing; ii a_m; ii a_n; ii* a_is0; ii* a_is1; ii* a_js; fp* a_vs;
    mkl_sparse_s_export_csr(*a.mat_, &a_indexing, &a_m, &a_n, &a_is0, &a_is1, &a_js, &a_vs);

	m_ = a.m_;
	n_ = a.n_;
    is_ = static_cast<ii*>(mkl_calloc(m_ + 1, sizeof(ii), 64));

	ii nnz = 0;
	for (ii a_nz = 0; a_nz < a_is1[a_m - 1]; a_nz++)
	{
		if (a_vs[a_nz] >= pruneThreshold)
		{
			nnz++;
		}
	}
    nnz_ = nnz; // important if zero
    js_ = static_cast<ii*>(mkl_malloc(sizeof(ii) * (nnz > 0 ? nnz : 1), 64));
    vs_ = static_cast<fp*>(mkl_malloc(sizeof(fp) * (nnz > 0 ? nnz : 1), 64));

	ii nz = 0;
	ii a_i = 0;
	for (ii a_nz = 0; a_nz < a_is1[a_m - 1]; a_nz++)
	{
		while (a_nz >= a_is1[a_i]) a_i++; // row of nz'th non-zero in a

		if (a_vs[a_nz] >= pruneThreshold)
		{
			is_[a_i + 1]++;
			js_[nz] = a_js[a_nz];
			vs_[nz] = a_vs[a_nz];
			nz++;
		}
	}
	for (ii i = 0; i < m_; i++)
	{
		is_[i + 1] += is_[i];
	}

    mat_ = new sparse_matrix_t;
	mkl_sparse_s_create_csr(mat_, SPARSE_INDEX_BASE_ZERO, m_, n_, is_, &is_[1], js_, vs_);
    data_ = Data::EXTERNAL;

#ifndef NDEBUG
	cout << "X" << *this << endl;
#endif
}


void MatrixSparse::output(fp* vs) const
{
#ifndef NDEBUG
    cout << "  X" << *this << " := " << flush;
#endif
    
    sparse_index_base_t indexing; ii m; ii n; ii* is0; ii* is1; ii* js; fp* vsIn;
    mkl_sparse_s_export_csr(*mat_, &indexing, &m, &n, &is0, &is1, &js, &vsIn);
    
    for (ii x = 0; x < m_ * n_; x++) vs[x] = 0.0;
    
    ii i = 0;
    for (ii nz = 0; nz < is1[m - 1]; nz++)
    {
        while (nz >= is1[i]) i++; // row of nz'th non-zero
        ii j = js[nz]; // column of nz'th non-zero
        
        vs[j + i * n_] = vsIn[nz];
    }
    
#ifndef NDEBUG
    cout << "out" << flush;
#endif
}


void MatrixSparse::deleteRows(const MatrixSparse& a)
{
#ifndef NDEBUG
    cout << "  X" << *this << ":deleteRows(A" << a << ") := " << flush;
#endif

    sparse_index_base_t indexing; ii m; ii n; ii* is0; ii* is1; ii* js; fp* vs;
    mkl_sparse_s_export_csr(*mat_, &indexing, &m, &n, &is0, &is1, &js, &vs);
    
    sparse_index_base_t a_indexing; ii a_m; ii a_n; ii* a_is0; ii* a_is1; ii* a_js; fp* a_vs;
    mkl_sparse_s_export_csr(*a.mat_, &a_indexing, &a_m, &a_n, &a_is0, &a_is1, &a_js, &a_vs);

    if (a.nnz() == 0)
    {
        // delete all of a
        init(m_, n_);
    }
    else
    {
        ii cellsDeleted = 0;
        ii a_nz = 0;
        for (ii i = 0; i < m_; i++)
        {
            // is this row to be deleted?
            bool del = (is1[i] - is0[i]) > 0;
            if (del)
            {
                //cout << "in" << endl;
                for (; a_nz < a_is1[a_m - 1] && a_js[a_nz] <= i; a_nz++)
                {
                    if (a_js[a_nz] == i)
                    {
                        del = false;
                        break;
                    }
                }
                //cout << "out" << endl;
            }
            
            ii cellsDeleted0 = cellsDeleted;
            
            if (del)
            {
                cellsDeleted += is1[i] - is0[i];
            }
            else if (cellsDeleted0 > 0)
            {
                // this row not deleted but has to be shifted
                for (ii nz = is0[i]; nz < is1[i]; nz++)
                {
                    js[nz - cellsDeleted0] = js[nz];
                    vs[nz - cellsDeleted0] = vs[nz];
                }
            }
            
            // ensure row index up to date
            if (cellsDeleted0 > 0) is0[i] -= cellsDeleted0;
        }
        if (cellsDeleted > 0) is1[m_ - 1] -= cellsDeleted;
        
        mkl_realloc(js, sizeof(ii) * is1[m_ - 1]); // not sure this does anything
        mkl_realloc(vs, sizeof(fp) * is1[m_ - 1]); // not sure this does anything
        
        if (data_ == Data::EXTERNAL || data_ == Data::EXTERNAL_MKL)
        {
            mkl_sparse_destroy(*mat_);
            mkl_sparse_s_create_csr(mat_, indexing, m, n, is0, is1, js, vs);
        }
        else // data_ = MKL
        {
            mkl_ = mat_;
            mat_ = new sparse_matrix_t;
            mkl_sparse_s_create_csr(mat_, indexing, m, n, is0, is1, js, vs);
            data_= Data::EXTERNAL_MKL;
        }
        
        nnz_ = is1[m_ - 1]; // incase it is ZERO
    }
 
#ifndef NDEBUG
    cout << "X" << *this << endl;
#endif
}


void MatrixSparse::mul(Accumulate accumulate, const MatrixSparse& a, Transpose transposeX, const MatrixSparse& b)
{
#ifndef NDEBUG
    cout << "  A" << (transposeX== Transpose::NO ? "" : "t") << a << " %*% B" << b;
    if (accumulate == Accumulate::YES) cout << " + X" << *this;
    cout << " := " << flush;
#endif

	if (!*this)
    {
        accumulate = Accumulate::NO;
    }

	// mkl can't handle completely empty csr matrix...
	if (b.nnz() == 0 || a.nnz() == 0)
	{
		if (accumulate == Accumulate::NO)
		{
			init(transposeX == Transpose::NO ? a.m_ : a.n_, b.n_);
		}
		else
		{
			// noop
		}
	}
	else
	{
		if (accumulate == Accumulate::NO)
		{
			free();
            mat_ = new sparse_matrix_t;
            mkl_sparse_spmm(transposeX == Transpose::NO ? SPARSE_OPERATION_NON_TRANSPOSE : SPARSE_OPERATION_TRANSPOSE, *a.mat_, *b.mat_, mat_);
		}
		else
		{
            sparse_matrix_t t;
            mkl_sparse_spmm(transposeX == Transpose::NO ? SPARSE_OPERATION_NON_TRANSPOSE : SPARSE_OPERATION_TRANSPOSE, *a.mat_, *b.mat_, &t);
            
            sparse_matrix_t* y = new sparse_matrix_t;
			mkl_sparse_s_add(SPARSE_OPERATION_NON_TRANSPOSE, *mat_, 1.0, t, y);
            mkl_sparse_destroy(t);
			free();
            
			mat_ = y;
		}
        
        data_ = Data::MKL;
        m_ = transposeX == Transpose::NO ? a.m_ : a.n_;
        n_ = b.n_;
	}

#ifndef NDEBUG
	cout << "X" << *this << endl;
#endif
}


void MatrixSparse::elementwiseSqr()
{
#ifndef NDEBUG
	cout << "  (X" << *this << ")^2 := " << flush;
#endif

    sparse_index_base_t indexing; ii m; ii n; ii* is0; ii* is1; ii* js; fp* vs;
    mkl_sparse_s_export_csr(*mat_, &indexing, &m, &n, &is0, &is1, &js, &vs);
    
	vsSqr(is1[m - 1], vs, vs);
    
    if (data_ == Data::EXTERNAL || data_ == Data::EXTERNAL_MKL)
    {
        mkl_sparse_destroy(*mat_);
        mkl_sparse_s_create_csr(mat_, indexing, m, n, is0, is1, js, vs);
    }
    else // data_ = MKL
    {
        mkl_ = mat_;
        mat_ = new sparse_matrix_t;
        mkl_sparse_s_create_csr(mat_, indexing, m, n, is0, is1, js, vs);
        data_= Data::EXTERNAL_MKL;
    }
    
#ifndef NDEBUG
	cout << "X" << *this << endl;
#endif
}


void MatrixSparse::elementwiseSqrt()
{
#ifndef NDEBUG
	cout << "  sqrt(X" << *this << ") := " << flush;
#endif

    sparse_index_base_t indexing; ii m; ii n; ii* is0; ii* is1; ii* js; fp* vs;
    mkl_sparse_s_export_csr(*mat_, &indexing, &m, &n, &is0, &is1, &js, &vs);

	vsSqrt(is1[m - 1], vs, vs);
    
    if (data_ == Data::EXTERNAL || data_ == Data::EXTERNAL_MKL)
    {
        mkl_sparse_destroy(*mat_);
        mkl_sparse_s_create_csr(mat_, indexing, m, n, is0, is1, js, vs);
    }
    else // data_ = MKL
    {
        mkl_ = mat_;
        mat_ = new sparse_matrix_t;
        mkl_sparse_s_create_csr(mat_, indexing, m, n, is0, is1, js, vs);
        data_= Data::EXTERNAL_MKL;
    }

#ifndef NDEBUG
	cout << "X" << *this << endl;
#endif
}


void MatrixSparse::elementwiseAdd(fp beta)
{
#ifndef NDEBUG
	cout << "  X" << *this << " + ";
	cout.unsetf(ios::floatfield);
	cout << setprecision(8) << beta << " := " << flush;
#endif

    sparse_index_base_t indexing; ii m; ii n; ii* is0; ii* is1; ii* js; fp* vs;
    mkl_sparse_s_export_csr(*mat_, &indexing, &m, &n, &is0, &is1, &js, &vs);
    
    ippsAddC_32f_I(beta, vs, is1[m - 1]);
    
    if (data_ == Data::EXTERNAL || data_ == Data::EXTERNAL_MKL)
    {
        mkl_sparse_destroy(*mat_);
        mkl_sparse_s_create_csr(mat_, indexing, m, n, is0, is1, js, vs);
    }
    else // data_ = MKL
    {
        mkl_ = mat_;
        mat_ = new sparse_matrix_t;
        mkl_sparse_s_create_csr(mat_, indexing, m, n, is0, is1, js, vs);
        data_= Data::EXTERNAL_MKL;
    }

#ifndef NDEBUG
	cout << "X" << *this << endl;
#endif
}


// todo: use IPP library
void MatrixSparse::elementwiseMul(fp beta)
{
#ifndef NDEBUG
	cout << "  X" << *this << " * ";
	cout.unsetf(ios::floatfield);
	cout << setprecision(8) << beta << " := " << flush;
#endif

    sparse_index_base_t indexing; ii m; ii n; ii* is0; ii* is1; ii* js; fp* vs;
    mkl_sparse_s_export_csr(*mat_, &indexing, &m, &n, &is0, &is1, &js, &vs);

    ippsMulC_32f_I(beta, vs, is1[m - 1]);

    if (data_ == Data::EXTERNAL || data_ == Data::EXTERNAL_MKL)
    {
        mkl_sparse_destroy(*mat_);
        mkl_sparse_s_create_csr(mat_, indexing, m, n, is0, is1, js, vs);
    }
    else // data_ = MKL
    {
        mkl_ = mat_;
        mat_ = new sparse_matrix_t;
        mkl_sparse_s_create_csr(mat_, indexing, m, n, is0, is1, js, vs);
        data_= Data::EXTERNAL_MKL;
    }

#ifndef NDEBUG
	cout << "X" << *this << endl;
#endif
}


void MatrixSparse::elementwiseMul(const MatrixSparse& a)
{
#ifndef NDEBUG
    cout << "  X" << *this << " / A" << a << " := " << flush;
#endif
    
    sparse_index_base_t indexing; ii m; ii n; ii* is0; ii* is1; ii* js; fp* vs;
    mkl_sparse_s_export_csr(*mat_, &indexing, &m, &n, &is0, &is1, &js, &vs);
    
    sparse_index_base_t a_indexing; ii a_m; ii a_n; ii* a_is0; ii* a_is1; ii* a_js; fp* a_vs;
    mkl_sparse_s_export_csr(*a.mat_, &a_indexing, &a_m, &a_n, &a_is0, &a_is1, &a_js, &a_vs);
    
    vsMul(is1[m - 1], vs, a_vs, vs);
    
    if (data_ == Data::EXTERNAL || data_ == Data::EXTERNAL_MKL)
    {
        mkl_sparse_destroy(*mat_);
        mkl_sparse_s_create_csr(mat_, indexing, m, n, is0, is1, js, vs);
    }
    else // data_ = MKL
    {
        mkl_ = mat_;
        mat_ = new sparse_matrix_t;
        mkl_sparse_s_create_csr(mat_, indexing, m, n, is0, is1, js, vs);
        data_= Data::EXTERNAL_MKL;
    }
    
#ifndef NDEBUG
    cout << "X" << *this << endl;
#endif
}


void MatrixSparse::elementwiseDiv(const MatrixSparse& a)
{
#ifndef NDEBUG
    cout << "  X" << *this << " / A" << a << " := " << flush;
#endif
    
    sparse_index_base_t indexing; ii m; ii n; ii* is0; ii* is1; ii* js; fp* vs;
    mkl_sparse_s_export_csr(*mat_, &indexing, &m, &n, &is0, &is1, &js, &vs);
    
    sparse_index_base_t a_indexing; ii a_m; ii a_n; ii* a_is0; ii* a_is1; ii* a_js; fp* a_vs;
    mkl_sparse_s_export_csr(*a.mat_, &a_indexing, &a_m, &a_n, &a_is0, &a_is1, &a_js, &a_vs);
    
    vsDiv(is1[m - 1], vs, a_vs, vs);
    
    if (data_ == Data::EXTERNAL || data_ == Data::EXTERNAL_MKL)
    {
        mkl_sparse_destroy(*mat_);
        mkl_sparse_s_create_csr(mat_, indexing, m, n, is0, is1, js, vs);
    }
    else // data_ = MKL
    {
        mkl_ = mat_;
        mat_ = new sparse_matrix_t;
        mkl_sparse_s_create_csr(mat_, indexing, m, n, is0, is1, js, vs);
        data_= Data::EXTERNAL_MKL;
    }
    
#ifndef NDEBUG
    cout << "X" << *this << endl;
#endif
}


fp MatrixSparse::sum() const
{
#ifndef NDEBUG
	cout << "  sum(X" << *this << ") := " << flush;
#endif

    sparse_index_base_t indexing; ii m; ii n; ii* is0; ii* is1; ii* js; fp* vs;
    mkl_sparse_s_export_csr(*mat_, &indexing, &m, &n, &is0, &is1, &js, &vs);
    
    fp sum = 0.0;
    ippsSum_32f(vs, is1[m - 1], &sum, ippAlgHintFast);

#ifndef NDEBUG
	cout.unsetf(ios::floatfield);
	cout << setprecision(8) << sum << endl;
#endif

	return sum;
}


fp MatrixSparse::sumSqrs() const
{
#ifndef NDEBUG
	cout << "  sum((X" << *this << ")^2) := " << flush;
#endif

    sparse_index_base_t indexing; ii m; ii n; ii* is0; ii* is1; ii* js; fp* vs;
    mkl_sparse_s_export_csr(*mat_, &indexing, &m, &n, &is0, &is1, &js, &vs);
    
    fp sum = 0.0;
	ippsNorm_L2_32f(vs, is1[m - 1], &sum);
    sum *= sum;

#ifndef NDEBUG
	cout.unsetf(ios::floatfield);
	cout << setprecision(8) << sum << endl;
#endif

	return sum;
}


fp MatrixSparse::sumSqrDiffs(const MatrixSparse& a) const
{
    
#ifndef NDEBUG
    cout << "  sum((X" << *this << " - A" << a << ")^2) = " << flush;
#endif
    
    sparse_index_base_t indexing; ii m; ii n; ii* is0; ii* is1; ii* js; fp* vs;
    mkl_sparse_s_export_csr(*mat_, &indexing, &m, &n, &is0, &is1, &js, &vs);
 
    sparse_index_base_t a_indexing; ii a_m; ii a_n; ii* a_is0; ii* a_is1; ii* a_js; fp* a_vs;
    mkl_sparse_s_export_csr(*a.mat_, &a_indexing, &a_m, &a_n, &a_is0, &a_is1, &a_js, &a_vs);

    fp sum = 0.0;
    ippsNormDiff_L2_32f(vs, a_vs, is1[m - 1], &sum);
    sum *= sum;
    
#ifndef NDEBUG
    cout.unsetf(ios::floatfield);
    cout << setprecision(8) << sum << endl;
#endif
    
    return sum;
}


// todo: use OpenMP
void MatrixSparse::subsetElementwiseCopy(const MatrixSparse& a)
{
#ifndef NDEBUG
    cout << "  subset(A" << a << ") := " << flush;
#endif
    
    sparse_index_base_t indexing; ii m; ii n; ii* is0; ii* is1; ii* js; fp* vs;
    mkl_sparse_s_export_csr(*mat_, &indexing, &m, &n, &is0, &is1, &js, &vs);
    
    sparse_index_base_t a_indexing; ii a_m; ii a_n; ii* a_is0; ii* a_is1; ii* a_js; fp* a_vs;
    mkl_sparse_s_export_csr(*a.mat_, &a_indexing, &a_m, &a_n, &a_is0, &a_is1, &a_js, &a_vs);
    
    ii i = 0;
    ii a_i = 0;
    ii a_nz = 0;
    for (ii nz = 0; nz < is1[m_ - 1]; nz++)
    {
        while (nz >= is1[i]) i++; // row of nz'th non-zero
        
        for (; a_nz < a_is1[m - 1]; a_nz++)
        {
            while (a_nz >= a_is1[a_i]) a_i++; // row of nz'th non-zero of a
            
            if (a_js[a_nz] == js[nz] && a_i == i)
            {
                vs[nz] = a_vs[a_nz];
                break;
            }
        }
    }
    
    if (data_ == Data::EXTERNAL || data_ == Data::EXTERNAL_MKL)
    {
        mkl_sparse_destroy(*mat_);
        mkl_sparse_s_create_csr(mat_, indexing, m, n, is0, is1, js, vs);
    }
    else // data_ = MKL
    {
        mkl_ = mat_;
        mat_ = new sparse_matrix_t;
        mkl_sparse_s_create_csr(mat_, indexing, m, n, is0, is1, js, vs);
        data_= Data::EXTERNAL_MKL;
    }
    
#ifndef NDEBUG
    cout << "X" << *this << endl;
#endif
}


// todo: use OpenMP
/*void MatrixSparse::subsetElementwiseMul(const MatrixSparse& a)
{
#ifndef NDEBUG
    cout << "  X" << *this << " * subset(" << a << ") := " << flush;
#endif
    
    sparse_index_base_t indexing; ii m; ii n; ii* is0; ii* is1; ii* js; fp* vs;
    mkl_sparse_s_export_csr(*mat_, &indexing, &m, &n, &is0, &is1, &js, &vs);
    
    sparse_index_base_t a_indexing; ii a_m; ii a_n; ii* a_is0; ii* a_is1; ii* a_js; fp* a_vs;
    mkl_sparse_s_export_csr(*a.mat_, &a_indexing, &a_m, &a_n, &a_is0, &a_is1, &a_js, &a_vs);
    
    ii i = 0;
    ii a_i = 0;
    ii a_nz = 0;
    for (ii nz = 0; nz < is1[m_ - 1]; nz++)
    {
        while (nz >= is1[i]) i++; // row of nz'th non-zero
        
        for (; a_nz < a_is1[m - 1]; a_nz++)
        {
            while (a_nz >= a_is1[a_i]) a_i++; // row of nz'th non-zero of a
            
            if (a_js[a_nz] == js[nz] && a_i == i)
            {
                vs[nz] *= a_vs[a_nz];
                break;
            }
        }
    }
    
    if (data_ == Data::EXTERNAL || data_ == Data::EXTERNAL_MKL)
    {
        mkl_sparse_destroy(*mat_);
        mkl_sparse_s_create_csr(mat_, indexing, m, n, is0, is1, js, vs);
    }
    else // data_ = MKL
    {
        mkl_ = mat_;
        mat_ = new sparse_matrix_t;
        mkl_sparse_s_create_csr(mat_, indexing, m, n, is0, is1, js, vs);
        data_= Data::EXTERNAL_MKL;
    }
    
#ifndef NDEBUG
    cout << "X" << *this << endl;
#endif
}*/


// todo: use OpenMP
void MatrixSparse::subsetElementwiseDiv(const MatrixSparse& a)
{
#ifndef NDEBUG
    cout << "  X" << *this << " / subset(" << a << ") := " << flush;
#endif
    
    sparse_index_base_t indexing; ii m; ii n; ii* is0; ii* is1; ii* js; fp* vs;
    mkl_sparse_s_export_csr(*mat_, &indexing, &m, &n, &is0, &is1, &js, &vs);
    
    sparse_index_base_t a_indexing; ii a_m; ii a_n; ii* a_is0; ii* a_is1; ii* a_js; fp* a_vs;
    mkl_sparse_s_export_csr(*a.mat_, &a_indexing, &a_m, &a_n, &a_is0, &a_is1, &a_js, &a_vs);
    
    ii i = 0;
    ii a_i = 0;
    ii a_nz = 0;
    for (ii nz = 0; nz < is1[m_ - 1]; nz++)
    {
        while (nz >= is1[i]) i++; // row of nz'th non-zero
        
        for (; a_nz < a_is1[m - 1]; a_nz++)
        {
            while (a_nz >= a_is1[a_i]) a_i++; // row of nz'th non-zero of a
            
            if (a_js[a_nz] == js[nz] && a_i == i)
            {
                vs[nz] /= a_vs[a_nz];
                break;
            }
        }
    }
    
    if (data_ == Data::EXTERNAL || data_ == Data::EXTERNAL_MKL)
    {
        mkl_sparse_destroy(*mat_);
        mkl_sparse_s_create_csr(mat_, indexing, m, n, is0, is1, js, vs);
    }
    else // data_ = MKL
    {
        mkl_ = mat_;
        mat_ = new sparse_matrix_t;
        mkl_sparse_s_create_csr(mat_, indexing, m, n, is0, is1, js, vs);
        data_= Data::EXTERNAL_MKL;
    }
    
#ifndef NDEBUG
    cout << "X" << *this << endl;
#endif
}


// todo: use OpenMP
/*double MatrixSparse::subsetSumSqrDiffs(const MatrixSparse& a) const
{
	double sum = 0.0;

#ifndef NDEBUG
	cout << "  sum((X" << *this << " - subset(A" << a << "))^2) = " << flush;
#endif
    
    sparse_index_base_t indexing; ii m; ii n; ii* is0; ii* is1; ii* js; fp* vs;
    mkl_sparse_s_export_csr(*mat_, &indexing, &m, &n, &is0, &is1, &js, &vs);
    
    sparse_index_base_t a_indexing; ii a_m; ii a_n; ii* a_is0; ii* a_is1; ii* a_js; fp* a_vs;
    mkl_sparse_s_export_csr(*a.mat_, &a_indexing, &a_m, &a_n, &a_is0, &a_is1, &a_js, &a_vs);
    
    ii i = 0;
    ii a_i = 0;
    ii a_nz = 0;
    for (ii nz = 0; nz < is1[m_ - 1]; nz++)
    {
        while (nz >= is1[i]) i++; // row of nz'th non-zero
        
        for (; a_nz < a_is1[m - 1]; a_nz++)
        {
            while (a_nz >= a_is1[a_i]) a_i++; // row of nz'th non-zero of a
            
            if (a_js[a_nz] == js[nz] && a_i == i)
            {
                sum += (vs[nz] - a.vs_[a_nz]) * (vs[nz] - a.vs_[a_nz]);
                break;
            }
        }
    }

#ifndef NDEBUG
	cout.unsetf(ios::floatfield);
	cout << setprecision(8) << sum << endl;
#endif

	return sum;
}*/


ostream& operator<<(ostream& os, const MatrixSparse& a)
{
	if (a.m() == 0)
	{
		os << "[]";
	}
	else
	{
		os << "{" << a.m() << "," << a.n() << "}:" << a.nnz() << "/" << a.size() << ":";
        os.unsetf(ios::floatfield);
        os << setprecision(3) << 100.0 * a.nnz() / (double)a.size() << "%";
	}

	return  os;
}
