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
#include <cmath>
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


static int debugLevel_ = 0;


void setDebugLevel(int debugLevel)
{
    debugLevel_ = debugLevel;
}


int getDebugLevel()
{
    return debugLevel_;
}


MatrixSparseMKL::MatrixSparseMKL() : m_(0), n_(0), is1_(0)
{
}


MatrixSparseMKL::~MatrixSparseMKL()
{
	free();
}


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
    isMklData_ = true;
    
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


void MatrixSparseMKL::free()
{
    if (is1_)
    {
        if (getDebugLevel() % 10 >= 4)
        {
            cout << getTimeStamp() << "     X" << *this << " := ..." << endl;
        }
        
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
        
        if (getDebugLevel() % 10 >= 4)
        {
            cout << getTimeStamp() << "     ... X" << *this << endl;
        }
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


void MatrixSparseMKL::copy(const MatrixSparseMKL& a, Operation operation)
{
    free();
    
    if (getDebugLevel() % 10 >= 4)
    {
        cout << getTimeStamp() << "     " << (operation == Operation::TRANSPOSE ? "t(" : "");
        cout << (operation == Operation::PACK_ROWS ? "pack(" : "") << (operation == Operation::UNPACK_ROWS ? "unpack(" : "");
        cout << "A" << a << (operation == Operation::PACK_ROWS || operation == Operation::UNPACK_ROWS ? ")" : "");
        cout << (operation == Operation::TRANSPOSE ? ")" : "") << " := ..." << endl;
    }
    
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


// todo: make this mkl
void MatrixSparseMKL::prune(const MatrixSparseMKL& a, fp pruneThreshold)
{
    free();
    
    if (getDebugLevel() % 10 >= 4)
    {
        cout << getTimeStamp() << "     (A" << a << " > ";
        cout.unsetf(ios::floatfield);
        cout << setprecision(8) << pruneThreshold << ") ? A : 0.0 := ..." << endl;
    }
    
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


void MatrixSparseMKL::zeroRowsOfZeroColumns(const MatrixSparseMKL& a, const MatrixSparseMKL& x)
{
    free();
    
    if (getDebugLevel() % 10 >= 4)
    {
        cout << getTimeStamp() << "     zeroRowsOfZeroColumns(A" << a << ", X" << x << ") := ..." << endl;
    }
    
    if (x.is1_)
    {
        // forces row vector to be column vector and hence have all row indicies (DOES NOT WORK IN 2D YET)
        MatrixSparseMKL t;
        t.copy(x, MatrixSparseMKL::Operation::TRANSPOSE);
        
        // make into unit diagonal matrix
        t.n_ = t.m_;
        t.setNonzeros((fp)1.0);
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
    
    if (getDebugLevel() % 10 >= 4)
    {
        cout << getTimeStamp() << "     ... X" << *this << endl;
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
        copy(a, transposeA ? Operation::TRANSPOSE : Operation::NONE);
        mul(alpha);
    }
    else
    {
        status_ = mkl_sparse_s_add(transposeA ? SPARSE_OPERATION_TRANSPOSE : SPARSE_OPERATION_NON_TRANSPOSE, a.mat_, alpha, b.mat_, &mat_); assert(!status_);
        
        sparse_index_base_t indexing;
        status_ = mkl_sparse_s_export_csr(mat_, &indexing, &m_, &n_, &is0_, &is1_, &js_, &vs_); assert(!status_);
        isMklData_ = true;
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


void MatrixSparseMKL::matmul(bool transposeA, const MatrixSparseMKL& a, const MatrixSparseMKL& b, bool accumulate)
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
    for (ii nz = 0; nz < nnz(); nz++) assert(js_[nz] == a.js_[nz]);
    
    vsMul(nnz(), vs_, a.vs_, vs_);
    
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


void MatrixSparseMKL::setNonzeros(fp v)
{
    if (getDebugLevel() % 10 >= 4)
    {
        cout << getTimeStamp() << "     (X" << *this << " != 0.0) ? " << v << " : 0.0 := ..." << endl;
    }
    
    for (ii nz = 0; nz < nnz(); nz++)
    {
        vs_[nz] = v;
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


void MatrixSparseMKL::powNonzeros(fp power)
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


void MatrixSparseMKL::divCorrespondingNonzeros(const MatrixSparseMKL& a)
{
    if (getDebugLevel() % 10 >= 4)
    {
        cout << getTimeStamp() << "     X" << *this << " / A" << a << " := ..." << endl;
    }
    
    assert(m_ == a.m_);
    assert(n_ == a.n_);
    for (ii nz = 0; nz < nnz(); nz++) assert(js_[nz] == a.js_[nz]);
    
    vsDiv(nnz(), vs_, a.vs_, vs_);
    
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


fp MatrixSparseMKL::sumSqrDiffsCorrespondingNonzeros(const MatrixSparseMKL& a) const
{
    if (getDebugLevel() % 10 >= 4)
    {
        cout << getTimeStamp() << "     sum((X" << *this << " - A" << a << ")^2) := ..." << endl;
    }
 
    assert(m_ == a.m_);
    assert(n_ == a.n_);
    
    fp sum = 0.0;
    ippsNormDiff_L2_32f(vs_, a.vs_, nnz(), &sum);
    sum *= sum;
    
    if (getDebugLevel() % 10 >= 4)
    {
        cout << getTimeStamp() << "     ... ";
        cout.unsetf(ios::floatfield);
        cout << setprecision(8) << sum << endl;
    }
    
    return sum;
}


// todo: use OpenMP
// todo: argh, mul output not sorted so a may not be
void MatrixSparseMKL::subsetElementwiseCopy(const MatrixSparseMKL& a)
{
    if (getDebugLevel() % 10 >= 4)
    {
        cout << getTimeStamp() << "     " << a << " within " << *this << " := ..." << endl;
    }
    
    assert(m_ == a.m_);
    assert(n_ == a.n_);
    
    ii nnz = 0;
    if (is1_ && a.is1_)
    {
        for (ii i = 0; i < m_; i++)
        {
            ii a_nz = a.is0_[i];
            for (ii nz = is0_[i]; nz < is1_[i]; nz++)
            {
                bool found = false;

                while(a_nz < a.is1_[i])
                {
                     if (js_[nz] == a.js_[a_nz])
                    {
                        vs_[nz] = a.vs_[a_nz];
                        a_nz++;
                        nnz++;
                        found = true;
                        break;
                    }
                    else
                    {
                        a_nz++;
                    }
                }
                
                if (!found) // go back to beginning because we might have missed it
                {
                    a_nz = a.is0_[i];
                    while(a_nz < a.is1_[i])
                    {
                        if (js_[nz] == a.js_[a_nz])
                        {
                            vs_[nz] = a.vs_[a_nz];
                            a_nz++;
                            nnz++;
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
