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


void printNumThreads()
{
    cout << "Config: " << 8 * sizeof(ii) << "bit MKL addressing, " << mkl_get_max_threads() << " MKL threads, ";
#if defined(_OPENMP)
    cout << omp_get_max_threads() << " OpenMP threads" << endl;
#else
    cout << "non-OpenMP build" << endl;
#endif
}


double getWallTime()
{
    return dsecnd();
}


MatrixSparse::MatrixSparse() : m_(0), n_(0), isEmpty_(true), isTransposed_(false)
{
}


MatrixSparse::~MatrixSparse()
{
	free();
}


void MatrixSparse::free()
{
    if (!isEmpty_)
    {
        mkl_sparse_destroy(mat_);
        
        if (!isMklData_)
        {
            mkl_free(is0_);
            mkl_free(js_);
            mkl_free(vs_);
        }
        
        isEmpty_ = true;
    }
}


/*void MatrixSparse::transpose()
{
#ifndef NDEBUG
    cout << "  t(X" << *this << ") := " << flush;
#endif
    
    if (isTransposed_)
    {
        isTransposed_ = false;
    }
    else
    {
        isTransposed_ = true;
    }
    
#ifndef NDEBUG
    cout << "X" << *this << endl;
#endif
}*/


ii MatrixSparse::m() const
{
    return (isTransposed_) ? n_ : m_;
}


ii MatrixSparse::n() const
{
    return (isTransposed_) ? m_ : n_;
}


li MatrixSparse::size() const
{
    return (li)m_ * n_;
}


bool MatrixSparse::isTransposed() const
{
    return isTransposed_;
}


ii MatrixSparse::nnz() const
{
    if (isEmpty_)
    {
        return 0;
    }
    else
    {
        return is0_[m_];
    }
}


li MatrixSparse::mem() const
{
    return sizeof(*this) + sizeof(fp) * (li)nnz() + (li)sizeof(ii) * ((li)nnz() + m_ + 1);
}


void MatrixSparse::init(ii m, ii n)
{
    free();
    
    m_ = m;
    n_ = n;
    isEmpty_ = true;
    isTransposed_ = false;
}


void MatrixSparse::init(ii m, ii n, ii nnz, const fp* acoo, const ii* rowind, const ii* colind)
{
#ifndef NDEBUG
    cout << "  COO := " << flush;
#endif
    
    free();
    
    m_ = m;
    n_ = n;
    isEmpty_ = false;
    isTransposed_ = false;
    
    fp* c_acoo = const_cast<fp*>(acoo);
    ii* c_rowind = const_cast<ii*>(rowind);
    ii* c_colind = const_cast<ii*>(colind);
   
    sparse_matrix_t t;
    status_ = mkl_sparse_s_create_coo(&t, SPARSE_INDEX_BASE_ZERO, m_, n_, nnz, c_rowind, c_colind, c_acoo); assert(!status_);
    status_ = mkl_sparse_convert_csr(t, SPARSE_OPERATION_NON_TRANSPOSE, &mat_); assert(!status_);
    
    sparse_index_base_t indexing;
    status_ = mkl_sparse_s_export_csr(mat_, &indexing, &m_, &n_, &is0_, &is1_, &js_, &vs_); assert(!status_);
    isMklData_ = true;
    
#ifndef NDEBUG
    cout << "X" << *this << endl;
#endif
}


void MatrixSparse::set(fp v)
{
#ifndef NDEBUG
    cout << "  (X" << *this << " != 0.0) ? " << v << " : 0.0 := " << flush;
#endif
    
    for (ii nz = 0; nz < nnz(); nz++)
    {
        vs_[nz] = v;
    }
    
#ifndef NDEBUG
    cout << "X" << *this << endl;
#endif
}


void MatrixSparse::copy(const MatrixSparse& a, Operation operation)
{
#ifndef NDEBUG
    cout << "  " << (operation == Operation::PACK_ROWS ? "pack(" : "") << (operation == Operation::UNPACK_ROWS ? "unpack(" : "");
    cout << "A" << a;
    cout << (operation == Operation::PACK_ROWS || operation == Operation::UNPACK_ROWS ? ")" : "") << " := " << flush;
#endif
    
    free();
    
    if (operation == Operation::TRANSPOSE)
    {
        if (a.isEmpty_)
        {
            init(a.n(), a.m());
        }
        else
        {
            mkl_sparse_convert_csr(a.mat_, SPARSE_OPERATION_TRANSPOSE, &mat_);
            
            sparse_index_base_t indexing;
            mkl_sparse_s_export_csr(mat_, &indexing, &m_, &n_, &is0_, &is1_, &js_, &vs_);
            isMklData_ = true;
            
            isEmpty_ = false;
            isTransposed_ = !a.isTransposed_;
        }
    }
    else
    {
        if (a.isEmpty_)
        {
            init(a.m(), a.n());
        }
        else
        {
            m_ = a.m_;
            isEmpty_ = false;
            isTransposed_ = a.isTransposed_;
            
            is0_ = static_cast<ii*>(mkl_malloc(sizeof(ii) * (a.m_ + 1), 64));
            is1_ = is0_ + 1;
            js_ = static_cast<ii*>(mkl_malloc(sizeof(ii) * a.is0_[m_], 64));
            vs_ = static_cast<fp*>(mkl_malloc(sizeof(fp) * a.is0_[m_], 64));
            
            memcpy(is0_, a.is0_, sizeof(ii) * (a.m_ + 1));
            memcpy(vs_, a.vs_, sizeof(fp) * a.is0_[m_]);

            if (operation == Operation::PACK_ROWS)
            {
                assert(!a.isTransposed_);
                
                n_ = a.n_ / a.m_;
                for (ii i = 0; i < a.m_; i++)
                {
                    for (ii nnz = a.is0_[i]; nnz < a.is1_[i]; nnz++)
                    {
                        js_[nnz] = a.js_[nnz] % n_;
                    }
                }
            }
            else if (operation == Operation::UNPACK_ROWS)
            {
                assert(!a.isTransposed_);
                
                n_ = a.n_ * a.m_;
                for (ii i = 0; i < a.m_; i++)
                {
                    for (ii nnz = a.is0_[i]; nnz < a.is1_[i]; nnz++)
                    {
                        //cout << i << "," << a.js_[nnz] << ":" << a.vs_[nnz] << " -> " << flush;
                        js_[nnz] = a.js_[nnz] + i * a.n_;
                        //cout << i << "," << js_[nnz] << ":" << vs_[nnz] << endl;
                    }
                }
            }
            else
            {
                n_ = a.n_;
                memcpy(js_, a.js_, sizeof(ii) * a.is0_[m_]);
            }

            status_ = mkl_sparse_s_create_csr(&mat_, SPARSE_INDEX_BASE_ZERO, m_, n_, is0_, is1_, js_, vs_); assert(!status_);
            isMklData_ = false;
        }
        
    }

#ifndef NDEBUG
	cout << "X" << *this << endl;
#endif
}


// todo: make this mkl
void MatrixSparse::copy(const MatrixSparse& a, fp pruneThreshold)
{
#ifndef NDEBUG
	cout << "  (A" << a << " > ";
	cout.unsetf(ios::floatfield);
	cout << setprecision(8) << pruneThreshold << ") ? A : 0.0 := " << flush;
#endif
    
    free();
    
	m_ = a.m_;
	n_ = a.n_;
    isTransposed_ = a.isTransposed_;

	ii nnz = 0;
    if (!a.isEmpty_)
    {
        for (ii a_nz = 0; a_nz < a.is0_[a.m_]; a_nz++)
        {
            if (a.vs_[a_nz] >= pruneThreshold)
            {
                nnz++;
            }
        }
    }
    isEmpty_ = (nnz == 0);
    
    if (!isEmpty_)
    {
        is0_ = static_cast<ii*>(mkl_calloc(m_ + 1, sizeof(ii), 64));
        is1_ = is0_ + 1;
        js_ = static_cast<ii*>(mkl_malloc(sizeof(ii) * nnz, 64));
        vs_ = static_cast<fp*>(mkl_malloc(sizeof(fp) * nnz, 64));
        assert(nnz > 0);
        
        ii nz = 0;
        ii a_i = 0;
        for (ii a_nz = 0; a_nz < a.is0_[a.m_]; a_nz++)
        {
            while (a_nz >= a.is1_[a_i]) a_i++; // row of nz'th non-zero in a
            
            if (a.vs_[a_nz] > pruneThreshold)
            {
                is1_[a_i]++;
                js_[nz] = a.js_[a_nz];
                vs_[nz] = a.vs_[a_nz];
                nz++;
            }
        }
        for (ii i = 0; i < m_; i++)
        {
            is1_[i] += is0_[i];
        }
        
        status_ = mkl_sparse_s_create_csr(&mat_, SPARSE_INDEX_BASE_ZERO, m_, n_, is0_, is1_, js_, vs_); assert(!status_);
        isMklData_ = false;
    }

#ifndef NDEBUG
	cout << "X" << *this << endl;
#endif
}


void MatrixSparse::output(fp* vs) const
{
#ifndef NDEBUG
    cout << "  X" << *this << " := " << flush;
#endif
    
    if (!isEmpty_)
    {
        for (li x = 0; x < size(); x++) vs[x] = 0.0;
        
        ii i = 0;
        for (ii nz = 0; nz < is0_[m_]; nz++)
        {
            while (nz >= is1_[i]) i++; // row of nz'th non-zero
            ii j = js_[nz]; // column of nz'th non-zero
            
            vs[j + i * n_] = vs_[nz];
        }
    }
    
#ifndef NDEBUG
    cout << "out" << endl;
#endif
}


void MatrixSparse::deleteRows(const MatrixSparse& a)
{
#ifndef NDEBUG
    cout << "  X" << *this << ":deleteRows(A" << a << ") := " << flush;
#endif

    if (!isEmpty_)
    {
        if (a.isEmpty_)
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
                bool del = (is1_[i] - is0_[i]) > 0;
                if (del)
                {
                    //cout << "in" << endl;
                    for (; a_nz < a.is0_[a.m_] && a.js_[a_nz] <= i; a_nz++)
                    {
                        if (a.js_[a_nz] == i)
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
                    cellsDeleted += is1_[i] - is0_[i];
                }
                else if (cellsDeleted0 > 0)
                {
                    // this row not deleted but has to be shifted
                    for (ii nz = is0_[i]; nz < is1_[i]; nz++)
                    {
                        js_[nz - cellsDeleted0] = js_[nz];
                        vs_[nz - cellsDeleted0] = vs_[nz];
                    }
                }
                
                // ensure row index up to date
                if (cellsDeleted0 > 0) is0_[i] -= cellsDeleted0;
            }
            if (cellsDeleted > 0) is1_[m_ - 1] -= cellsDeleted;
            
            mkl_realloc(js_, sizeof(ii) * is0_[m_]); // not sure this does anything
            mkl_realloc(vs_, sizeof(fp) * is0_[m_]); // not sure this does anything
        }
    }
 
#ifndef NDEBUG
    cout << "X" << *this << endl;
#endif
}


void MatrixSparse::mul(bool transposeA, const MatrixSparse& a, const MatrixSparse& b, bool accumulate, bool transpose)
{
#ifndef NDEBUG
    cout << "  " << (transpose ? "t(" : "") << (transposeA ? "t(" : "") << "A" << a << (transposeA ? ")" : "") << " %*% B" << b;
    if (transpose) cout << ")" << flush;
    if (accumulate) cout << " + X" << *this;
    cout << " := " << flush;
    
    assert((transposeA ? a.m() : a.n()) == b.m());
#endif
    
    if (isEmpty_)
    {
        accumulate = false;
    }
    
    transposeA = a.isTransposed_ ? !transposeA : transposeA;

	if (a.isEmpty_ || b.isEmpty_)
	{
		if (!accumulate)
		{
			init(transposeA ? a.n() : a.m(), b.n());
		}
	}
	else
	{
		if (accumulate)
        {
            sparse_matrix_t t;
            status_ = mkl_sparse_spmm(transposeA ? SPARSE_OPERATION_TRANSPOSE : SPARSE_OPERATION_NON_TRANSPOSE, a.mat_, b.mat_, &t); assert(!status_);
            
            sparse_matrix_t y;
            status_ = mkl_sparse_s_add(transpose == isTransposed_ ? SPARSE_OPERATION_NON_TRANSPOSE : SPARSE_OPERATION_TRANSPOSE, mat_, 1.0, t, &y); assert(!status_);
            status_ = mkl_sparse_destroy(t); assert(!status_);
            free();
            
            mat_ = y;
            
            isTransposed_ = transpose;
        }
        else
		{
			free();
            status_ = mkl_sparse_spmm(transposeA ? SPARSE_OPERATION_TRANSPOSE : SPARSE_OPERATION_NON_TRANSPOSE, a.mat_, b.mat_, &mat_); assert(!status_);
            
            isTransposed_ = transpose;
		}

        isEmpty_ = false;
        
        sparse_index_base_t indexing;
        status_ = mkl_sparse_s_export_csr(mat_, &indexing, &m_, &n_, &is0_, &is1_, &js_, &vs_); assert(!status_);
        isMklData_ = true;
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

    if (!isEmpty_)
    {
        vsSqr(is0_[m_], vs_, vs_);
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

    if (!isEmpty_)
    {
        vsSqrt(is0_[m_], vs_, vs_);
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
    
    if (!isEmpty_)
    {
        ippsAddC_32f_I(beta, vs_, is0_[m_]);
    }

#ifndef NDEBUG
	cout << "X" << *this << endl;
#endif
}


void MatrixSparse::elementwiseMul(fp beta)
{
#ifndef NDEBUG
	cout << "  X" << *this << " * ";
	cout.unsetf(ios::floatfield);
	cout << setprecision(8) << beta << " := " << flush;
#endif

    if (!isEmpty_)
    {
        ippsMulC_32f_I(beta, vs_, is0_[m_]);
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
    
    if (!isEmpty_)
    {
        vsMul(is0_[m_], vs_, a.vs_, vs_);
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
    
    if (!isEmpty_)
    {
        vsDiv(is0_[m_], vs_, a.vs_, vs_);
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

    fp sum = 0.0;
    
    if (!isEmpty_)
    {
        ippsSum_32f(vs_, is0_[m_], &sum, ippAlgHintFast);
    }
    
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

    fp sum = 0.0;
    
    if (!isEmpty_)
    {
        ippsNorm_L2_32f(vs_, is0_[m_], &sum);
        sum *= sum;
    }

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
    
    fp sum = 0.0;
    
    if (!isEmpty_)
    {
        ippsNormDiff_L2_32f(vs_, a.vs_, is0_[m_], &sum);
        sum *= sum;
    }
    
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
    
    if (!isEmpty_)
    {
        ii i = 0;
        ii a_i = 0;
        ii a_nz = 0;
        for (ii nz = 0; nz < is0_[m_]; nz++)
        {
            while (nz >= is1_[i]) i++; // row of nz'th non-zero
            
            for (; a_nz < a.is0_[m_]; a_nz++)
            {
                while (a_nz >= a.is1_[a_i]) a_i++; // row of nz'th non-zero of a
                
                if (a.js_[a_nz] == js_[nz] && a_i == i)
                {
                    vs_[nz] = a.vs_[a_nz];
                    break;
                }
            }
        }
    }
    
#ifndef NDEBUG
    cout << "X" << *this << endl;
#endif
}


// todo: use OpenMP
void MatrixSparse::subsetElementwiseDiv(const MatrixSparse& a)
{
#ifndef NDEBUG
    cout << "  X" << *this << " / subset(A" << a << ") := " << flush;
#endif
    
    if (!isEmpty_ && !a.isEmpty_)
    {
        ii i = 0;
        ii a_i = 0;
        ii a_nz = 0;
        for (ii nz = 0; nz < is0_[m_]; nz++)
        {
            while (nz >= is1_[i]) i++; // row of nz'th non-zero
            
            for (; a_nz < a.is0_[m_]; a_nz++)
            {
                while (a_nz >= a.is1_[a_i]) a_i++; // row of nz'th non-zero of a
                
                if (a.js_[a_nz] == js_[nz] && a_i == i)
                {
                    vs_[nz] /= a.vs_[a_nz];
                    break;
                }
            }
        }
    }
    
#ifndef NDEBUG
    cout << "X" << *this << endl;
#endif
}


ostream& operator<<(ostream& os, const MatrixSparse& a)
{
	if (a.m() == 0)
	{
		os << "[]";
	}
	else
	{
        os << (a.isTransposed_ ? "t" : "") << "[" << a.m_ << "," << a.n_ << "]:" << a.nnz() << "/" << a.size() << ":";
        os.unsetf(ios::floatfield);
        os << setprecision(3) << 100.0 * a.nnz() / (double)a.size() << "%";
	}

	return  os;
}
