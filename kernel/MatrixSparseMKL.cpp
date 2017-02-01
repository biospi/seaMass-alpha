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
#include <iostream>
#include <cassert>
#include <cstring>
#include <sstream>

#if defined(_OPENMP)
  #include <omp.h>
#endif

#include <ippcore.h>
#include <ipps.h>


#include "../kernel/NetcdfFile.hpp"


using namespace std;


std::string getThreadInfo()
{
    ostringstream out;
    out << "Config: " << 8 * sizeof(ii) << "bit MKL addressing, " << mkl_get_max_threads() << " MKL threads, ";
#if defined(_OPENMP)
    out << omp_get_max_threads() << " OpenMP threads" << endl;
#else
    out << "non-OpenMP build" << endl;
#endif
    return out.str();
}


static li id_ = -1;


li getId()
{
    return id_;
}


static double startTime_ = dsecnd();


void resetElapsedTime()
{
    startTime_ = dsecnd();
}


double getElapsedTime()
{
    return dsecnd() - startTime_;
}


li getUsedMemory()
{
    int allocatedBuffers;
    return mkl_mem_stat(&allocatedBuffers);
}


string getTimeStamp()
{
    ostringstream out;
    out << "[" << setw(9) << ++id_ << "," << fixed << internal << setw(9) << std::setprecision(3) << getElapsedTime() << "," << setw(9) << getUsedMemory()/1024.0/1024.0 << "] ";
    return out.str();
}


MatrixSparseMKL::MatrixSparseMKL() : m_(0), n_(0), is1_(0)
{
}


MatrixSparseMKL::~MatrixSparseMKL()
{
	free();
}


void MatrixSparseMKL::free()
{
    if (is1_)
    {
#ifndef NDEBUG
        cout << getTimeStamp() << "   X" << *this << " := ..." << endl;
#endif
        
        status_ = mkl_sparse_destroy(mat_); assert(!status_);
        
        if (!isMklData_)
        {
            mkl_free(is0_);
            mkl_free(js_);
            mkl_free(vs_);
        }
        
        m_ = 0;
        n_ = 0;
        is1_ = 0;
        
#ifndef NDEBUG
        cout << getTimeStamp() << "   ... X" << *this << endl;
#endif
    }
    else
    {
        m_ = 0;
        n_ = 0;
    }
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


void MatrixSparseMKL::init(ii m, ii n, ii nnz, const fp* acoo, const ii* rowind, const ii* colind)
{
    free();
    
#ifndef NDEBUG
    cout << getTimeStamp() << "   COO := ..." << endl;
#endif
    
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
    isMklData_ = true;
    
#ifndef NDEBUG
    cout << getTimeStamp() << "   ... X" << *this << endl;
#endif
}


void MatrixSparseMKL::set(fp v)
{
#ifndef NDEBUG
    cout << getTimeStamp() << "   (X" << *this << " != 0.0) ? " << v << " : 0.0 := ..." << endl;
#endif
    
    for (ii nz = 0; nz < nnz(); nz++)
    {
        vs_[nz] = v;
    }
    
#ifndef NDEBUG
    cout << getTimeStamp() << "   ... X" << *this << endl;
#endif
}


void MatrixSparseMKL::copy(const MatrixSparseMKL& a, Operation operation)
{
    free();
    
#ifndef NDEBUG
    cout << getTimeStamp() << "   " << (operation == Operation::TRANSPOSE ? "t(" : "");
    cout << (operation == Operation::PACK_ROWS ? "pack(" : "") << (operation == Operation::UNPACK_ROWS ? "unpack(" : "");
    cout << "A" << a << (operation == Operation::PACK_ROWS || operation == Operation::UNPACK_ROWS ? ")" : "");
    cout << (operation == Operation::TRANSPOSE ? ")" : "") << " := ..." << endl;
#endif
    
    if (operation == Operation::TRANSPOSE)
    {
        mkl_sparse_convert_csr(a.mat_, SPARSE_OPERATION_TRANSPOSE, &mat_);
            
        sparse_index_base_t indexing;
        mkl_sparse_s_export_csr(mat_, &indexing, &m_, &n_, &is0_, &is1_, &js_, &vs_);
        isMklData_ = true;
    }
    else
    {
        m_ = a.m_;
        if (operation == Operation::PACK_ROWS)
        {
            n_ = a.n_ / a.m_;
        }
        else if (operation == Operation::UNPACK_ROWS)
        {
            n_ = a.n_ * a.m_;
        }
        else
        {
            n_ = a.n_;
        }
        
        if (a.is1_)
        {
            is0_ = static_cast<ii*>(mkl_malloc(sizeof(ii) * (a.m_ + 1), 64));
            is1_ = is0_ + 1;
            js_ = static_cast<ii*>(mkl_malloc(sizeof(ii) * a.is1_[m_ - 1], 64));
            vs_ = static_cast<fp*>(mkl_malloc(sizeof(fp) * a.is1_[m_ - 1], 64));
            
            memcpy(is0_, a.is0_, sizeof(ii) * (a.m_ + 1));
            memcpy(vs_, a.vs_, sizeof(fp) * a.is1_[m_ - 1]);

            if (operation == Operation::PACK_ROWS)
            {
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
                for (ii i = 0; i < a.m_; i++)
                {
                    for (ii nnz = a.is0_[i]; nnz < a.is1_[i]; nnz++)
                    {
                        js_[nnz] = a.js_[nnz] + i * a.n_;
                    }
                }
            }
            else
            {
                memcpy(js_, a.js_, sizeof(ii) * a.is0_[m_]);
            }

            status_ = mkl_sparse_s_create_csr(&mat_, SPARSE_INDEX_BASE_ZERO, m_, n_, is0_, is1_, js_, vs_); assert(!status_);
            isMklData_ = false;
        }
    }

#ifndef NDEBUG
    cout << getTimeStamp() << "   ... X" << *this << endl;
#endif
}


// todo: make this mkl
void MatrixSparseMKL::prune(const MatrixSparseMKL& a, fp pruneThreshold)
{
    free();
    
#ifndef NDEBUG
	cout << getTimeStamp() << "   (A" << a << " > ";
	cout.unsetf(ios::floatfield);
	cout << setprecision(8) << pruneThreshold << ") ? A : 0.0 := ..." << endl;
#endif
    
	m_ = a.m_;
	n_ = a.n_;

	ii nnz = 0;
    for (ii a_nz = 0; a_nz < a.nnz(); a_nz++)
    {
        if (a.vs_[a_nz] >= pruneThreshold)
        {
            nnz++;
        }
    }
    
    if (nnz > 0)
    {
        is0_ = static_cast<ii*>(mkl_calloc(m_ + 1, sizeof(ii), 64));
        is1_ = is0_ + 1;
        js_ = static_cast<ii*>(mkl_malloc(sizeof(ii) * (nnz == 0 ? 1 : nnz), 64));
        vs_ = static_cast<fp*>(mkl_malloc(sizeof(fp) * (nnz == 0 ? 1 : nnz), 64));
        
        ii nz = 0;
        ii a_i = 0;
        for (ii a_nz = 0; a_nz < a.nnz(); a_nz++)
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
    cout << getTimeStamp() << "   ... X" << *this << endl;
#endif
}


void MatrixSparseMKL::output(fp* vs) const
{
#ifndef NDEBUG
    cout << getTimeStamp() << "   X" << *this << " := ..." << endl;
#endif
    
    for (li x = 0; x < size(); x++) vs[x] = 0.0;
    
    ii i = 0;
    for (ii nz = 0; nz < nnz(); nz++)
    {
        while (nz >= is1_[i]) i++; // row of nz'th non-zero
        ii j = js_[nz]; // column of nz'th non-zero
        
        vs[i + j * m_] = vs_[nz];
    }
    
#ifndef NDEBUG
    cout << getTimeStamp() << "   ... out" << endl;
#endif
}


void MatrixSparseMKL::zeroRowsOfZeroColumns(const MatrixSparseMKL& a, const MatrixSparseMKL& x)
{
    free();
    
#ifndef NDEBUG
    cout << getTimeStamp() << "   zeroRowsOfZeroColumns(A" << a << ", X" << x << ") := ..." << endl;
#endif
    
    if (x.is1_)
    {
        // forces row vector to be column vector and hence have all row indicies (DOES NOT WORK IN 2D YET)
        MatrixSparseMKL t;
        t.copy(x, MatrixSparseMKL::Operation::TRANSPOSE);
        
        // make into unit diagonal matrix
        t.n_ = t.m_;
        t.set((fp)1.0);
        for (ii i = 0; i < t.m_; i++)
        {
            for (ii nz = t.is0_[i]; nz < t.is1_[i]; nz++)
            {
                t.js_[nz] = i;
            }
        }
        
        sparse_matrix_t diag;
        status_ = mkl_sparse_s_create_csr(&diag, SPARSE_INDEX_BASE_ZERO, t.m_, t.n_, t.is0_, t.is1_, t.js_, t.vs_); assert(!status_);
        status_ = mkl_sparse_spmm(SPARSE_OPERATION_NON_TRANSPOSE, diag, a.mat_, &mat_); assert(!status_);
        status_ = mkl_sparse_destroy(diag); assert(!status_);
        
        sparse_index_base_t indexing;
        status_ = mkl_sparse_s_export_csr(mat_, &indexing, &m_, &n_, &is0_, &is1_, &js_, &vs_); assert(!status_);
        isMklData_ = true;
    }
    else
    {
        m_ = a.m_;
        n_ = a.n_;
    }
    
#ifndef NDEBUG
    cout << getTimeStamp() << "   ... X" << *this << endl;
#endif
}


void MatrixSparseMKL::mul(bool transposeA, const MatrixSparseMKL& a, const MatrixSparseMKL& b, bool accumulate, bool transpose)
{
    if (!accumulate) free();
    
#ifndef NDEBUG
    cout << getTimeStamp() << "   " << (transpose ? "t(" : "") << (transposeA ? "t(" : "") << "A" << a << (transposeA ? ")" : "") << " %*% B" << b;
    if (transpose) cout << ")";
    if (accumulate) cout << " + X" << *this;
    cout << " := ..." << endl;
    
    assert((transposeA ? a.m() : a.n()) == b.m());
#endif
    
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
		if (accumulate)
        {
            sparse_matrix_t t;
            status_ = mkl_sparse_spmm(transposeA ? SPARSE_OPERATION_TRANSPOSE : SPARSE_OPERATION_NON_TRANSPOSE, a.mat_, b.mat_, &t); assert(!status_);
 
            sparse_matrix_t y;
            status_ = mkl_sparse_s_add(SPARSE_OPERATION_NON_TRANSPOSE, mat_, 1.0, t, &y); assert(!status_);
            status_ = mkl_sparse_destroy(t); assert(!status_);
            
            free();
            mat_ = y;
        }
        else
		{
            status_ = mkl_sparse_spmm(transposeA ? SPARSE_OPERATION_TRANSPOSE : SPARSE_OPERATION_NON_TRANSPOSE, a.mat_, b.mat_, &mat_); assert(!status_);
		}
        
        sparse_index_base_t indexing;
        status_ = mkl_sparse_s_export_csr(mat_, &indexing, &m_, &n_, &is0_, &is1_, &js_, &vs_); assert(!status_);
        isMklData_ = true;
	}

#ifndef NDEBUG
    cout << getTimeStamp() << "   ... X" << *this << endl;
    
    if (m_ > 0)
    {
        ostringstream oss; oss << getId() << ".csr";
        NetCDFile outFile(oss.str(), NC_NETCDF4);
        
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

#endif
}


void MatrixSparseMKL::elementwiseSqr()
{
#ifndef NDEBUG
	cout << getTimeStamp() << "   (X" << *this << ")^2 := ..." << endl;
#endif

    vsSqr(nnz(), vs_, vs_);
    
#ifndef NDEBUG
    cout << getTimeStamp() << "   ... X" << *this << endl;
#endif
}


void MatrixSparseMKL::elementwiseSqrt()
{
#ifndef NDEBUG
	cout << getTimeStamp() << "   sqrt(X" << *this << ") := ..." << endl;
#endif

    vsSqrt(nnz(), vs_, vs_);

#ifndef NDEBUG
    cout << getTimeStamp() << "   ... X" << *this << endl;
#endif
}


void MatrixSparseMKL::elementwiseAdd(fp beta)
{
#ifndef NDEBUG
	cout << getTimeStamp() << "   X" << *this << " + ";
	cout.unsetf(ios::floatfield);
	cout << setprecision(8) << beta << " := ..." << endl;
#endif
    
    ippsAddC_32f_I(beta, vs_, nnz());

#ifndef NDEBUG
    cout << getTimeStamp() << "   ... X" << *this << endl;
#endif
}


void MatrixSparseMKL::elementwiseMul(fp beta)
{
#ifndef NDEBUG
	cout << getTimeStamp() << "   X" << *this << " * ";
	cout.unsetf(ios::floatfield);
	cout << setprecision(8) << beta << " := ..." << endl;
#endif

    ippsMulC_32f_I(beta, vs_, nnz());

#ifndef NDEBUG
    cout << getTimeStamp() << "   ... X" << *this << endl;
#endif
}


void MatrixSparseMKL::elementwiseMul(const MatrixSparseMKL& a)
{
#ifndef NDEBUG
    cout << getTimeStamp() << "   X" << *this << " * A" << a << " := ..." << endl;
#endif
    
    vsMul(nnz(), vs_, a.vs_, vs_);
    
#ifndef NDEBUG
    cout << getTimeStamp() << "   ... X" << *this << endl;
#endif
}


void MatrixSparseMKL::elementwiseDiv(const MatrixSparseMKL& a)
{
#ifndef NDEBUG
    cout << getTimeStamp() << "   X" << *this << " / A" << a << " := ..." << endl;
#endif
    
    vsDiv(nnz(), vs_, a.vs_, vs_);
    
#ifndef NDEBUG
    cout << getTimeStamp() << "   ... X" << *this << endl;
#endif
}


fp MatrixSparseMKL::sum() const
{
#ifndef NDEBUG
	cout << getTimeStamp() << "   sum(X" << *this << ") := ..." << endl;
#endif

    fp sum = 0.0;
    ippsSum_32f(vs_, nnz(), &sum, ippAlgHintFast);
    
#ifndef NDEBUG
    cout << getTimeStamp() << "   ... ";
	cout.unsetf(ios::floatfield);
	cout << setprecision(8) << sum << endl;
#endif

	return sum;
}


#include <math.h>


fp MatrixSparseMKL::sumSqrs() const
{
#ifndef NDEBUG
	cout << getTimeStamp() << "   sum((X" << *this << ")^2) := ..." << endl;
#endif

    fp sum = 0.0;
    ippsNorm_L2_32f(vs_, nnz(), &sum);
    sum *= sum;

#ifndef NDEBUG
    cout << getTimeStamp() << "   ... ";
	cout.unsetf(ios::floatfield);
	cout << setprecision(8) << sum << endl;
#endif

	return sum;
}


fp MatrixSparseMKL::sumSqrDiffs(const MatrixSparseMKL& a) const
{
#ifndef NDEBUG
    cout << getTimeStamp() << "   sum((X" << *this << " - A" << a << ")^2) := ..." << endl;
#endif
    
    fp sum = 0.0;
    ippsNormDiff_L2_32f(vs_, a.vs_, nnz(), &sum);
    sum *= sum;
    
#ifndef NDEBUG
    cout << getTimeStamp() << "   ... ";
    cout.unsetf(ios::floatfield);
    cout << setprecision(8) << sum << endl;
#endif
    
    return sum;
}


// todo: use OpenMP
void MatrixSparseMKL::subsetElementwiseCopy(const MatrixSparseMKL& a)
{
#ifndef NDEBUG
    cout << getTimeStamp() << "   " << *this << ":subset(A" << a << ") := ..." << endl;
#endif
    
    if (is1_)
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
    cout << getTimeStamp() << "   ... X" << *this << endl;
#endif
}


// todo: use OpenMP
void MatrixSparseMKL::subsetElementwiseDiv(const MatrixSparseMKL& a)
{
#ifndef NDEBUG
    cout << getTimeStamp() << "   X" << *this << " / subset(A" << a << ", " << *this << ") := ..." << endl;
#endif
    
    if (is1_ && a.is1_)
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
    cout << getTimeStamp() << "   ... X" << *this << endl;
#endif
}


ostream& operator<<(ostream& os, const MatrixSparseMKL& a)
{
	if (a.m() == 0)
	{
		os << "[]";
	}
	else
	{
        os << "[" << a.m_ << "," << a.n_ << "]:" << a.nnz() << "/" << a.size() << ":";
        os.unsetf(ios::floatfield);
        os << setprecision(3) << 100.0 * a.nnz() / (double)a.size() << "%";
	}

	return  os;
}
