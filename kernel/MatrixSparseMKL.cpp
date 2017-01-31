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


static li id_ = 0;


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
    out << "[" << setw(9) << id_++ << "," << fixed << internal << setw(9) << std::setprecision(3) << getElapsedTime() << "," << setw(9) << getUsedMemory()/1024.0/1024.0 << "] ";
    return out.str();
}


MatrixSparseMKL::MatrixSparseMKL() : m_(0), n_(0), isEmpty_(true), isTransposed_(false)
{
}


MatrixSparseMKL::~MatrixSparseMKL()
{
	free();
}


void MatrixSparseMKL::free()
{
    if (!isEmpty_)
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
        
        isEmpty_ = true;
        isTransposed_ = false;
        m_ = 0;
        n_ = 0;
        
#ifndef NDEBUG
        cout << getTimeStamp() << "   ... X" << *this << endl;
#endif
    }
}


/*void MatrixSparseMKL::transpose()
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


ii MatrixSparseMKL::m() const
{
    return (isTransposed_) ? n_ : m_;
}


ii MatrixSparseMKL::n() const
{
    return (isTransposed_) ? m_ : n_;
}


li MatrixSparseMKL::size() const
{
    return (li)m_ * n_;
}


bool MatrixSparseMKL::isTransposed() const
{
    return isTransposed_;
}


ii MatrixSparseMKL::nnz() const
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


/*li MatrixSparseMKL::mem() const
{
    return sizeof(*this) + sizeof(fp) * (li)nnz() + (li)sizeof(ii) * ((li)nnz() + m_ + 1);
}*/


void MatrixSparseMKL::init(ii m, ii n)
{
    free();
    
    m_ = m;
    n_ = n;
    isEmpty_ = true;
    isTransposed_ = false;
}


void MatrixSparseMKL::init(ii m, ii n, ii nnz, const fp* acoo, const ii* rowind, const ii* colind)
{
    free();
    
#ifndef NDEBUG
    cout << getTimeStamp() << "   COO := ..." << endl;
#endif
    
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
            isTransposed_ = a.isTransposed_;
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
                        js_[nnz] = a.js_[nnz] + i * a.n_;
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
    cout << getTimeStamp() << "   ... X" << *this << endl;
#endif
}


void MatrixSparseMKL::output(fp* vs) const
{
#ifndef NDEBUG
    cout << getTimeStamp() << "   X" << *this << " := ..." << endl;
#endif
    
    if (!isEmpty_)
    {
        for (li x = 0; x < size(); x++) vs[x] = 0.0;
        
        ii i = 0;
        for (ii nz = 0; nz < is0_[m_]; nz++)
        {
            while (nz >= is1_[i]) i++; // row of nz'th non-zero
            ii j = js_[nz]; // column of nz'th non-zero
            
            vs[i + j * m_] = vs_[nz];
            //cout << vs_[nz] << endl;
        }
    }
    
#ifndef NDEBUG
    cout << getTimeStamp() << "   ... out" << endl;
#endif
}


/*void MatrixSparseMKL::deleteRows(const MatrixSparseMKL& a)
{
//#ifndef NDEBUG
    cout << "  X" << *this << ":deleteRows(A" << a << ") := " << flush;
//#endif

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
 
//#ifndef NDEBUG
    cout << "X" << *this << endl;
    cout << "   mem=" << fixed << setprecision(2) << getUsedMemory()/1024.0/1024.0 << "Mb time=" << setprecision(10) << getWallTime << endl;
//#endif
}*/


void MatrixSparseMKL::zeroRowsOfZeroColumns(const MatrixSparseMKL& a, const MatrixSparseMKL& x)
{
    free();
    
#ifndef NDEBUG
    cout << getTimeStamp() << "   zeroRowsOfZeroColumns(A" << a << ", X" << x << ") := ..." << endl;
#endif
    
    if (x.isEmpty_)
    {
        // delete all
        init(m_, n_);
    }
    else
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
        
        isTransposed_ = false;
        isEmpty_ = false;
        
        sparse_index_base_t indexing;
        status_ = mkl_sparse_s_export_csr(mat_, &indexing, &m_, &n_, &is0_, &is1_, &js_, &vs_); assert(!status_);
        isMklData_ = true;
    }
    
#ifndef NDEBUG
    cout << getTimeStamp() << "   ... X" << *this << endl;
#endif
}


void MatrixSparseMKL::mul(bool transposeA, const MatrixSparseMKL& a, const MatrixSparseMKL& b, bool accumulate, bool transpose) {
    if (!accumulate) free();

#ifndef NDEBUG
    cout << getTimeStamp() << "   " << (transpose ? "t(" : "") << (transposeA ? "t(" : "")
         << "A" << a << (transposeA ? ")" : "") << " %*% B" << b;
    if (transpose) cout << ")";
    if (accumulate) cout << " + X" << *this;
    cout << " := ..." << endl;

    assert((transposeA ? a.m() : a.n()) == b.m());
#endif

    if (isEmpty_) {
        accumulate = false;
    }

    transposeA = a.isTransposed_ ? !transposeA : transposeA;

    if (a.isEmpty_ || b.isEmpty_) {
        if (!accumulate) {
            init(transposeA ? a.n() : a.m(), b.n());
        }
    } else {
        if (accumulate) {
            sparse_matrix_t t;
            status_ = mkl_sparse_spmm(transposeA ? SPARSE_OPERATION_TRANSPOSE
                                                 : SPARSE_OPERATION_NON_TRANSPOSE, a.mat_,
                                      b.mat_, &t);
            assert(!status_);

            sparse_matrix_t y;
            status_ = mkl_sparse_s_add(
                    transpose == isTransposed_ ? SPARSE_OPERATION_NON_TRANSPOSE
                                               : SPARSE_OPERATION_TRANSPOSE, mat_, 1.0, t,
                    &y);
            assert(!status_);
            status_ = mkl_sparse_destroy(t);
            assert(!status_);

            free();
            mat_ = y;
        } else {
            status_ = mkl_sparse_spmm(transposeA ? SPARSE_OPERATION_TRANSPOSE
                                                 : SPARSE_OPERATION_NON_TRANSPOSE, a.mat_,
                                      b.mat_, &mat_);
            assert(!status_);
        }

        isTransposed_ = transpose;
        isEmpty_ = false;

        sparse_index_base_t indexing;
        status_ = mkl_sparse_s_export_csr(mat_, &indexing, &m_, &n_, &is0_, &is1_, &js_,
                                          &vs_);
        assert(!status_);
        isMklData_ = true;
    }

#ifndef NDEBUG
    cout << getTimeStamp() << "   ... X" << *this << endl;
    vector<ii> dims;
    dims.push_back(m_);
    dims.push_back(n_);

    ostringstream oss;
    oss << getId() << ".csr";
    NetCDFile outFile(oss.str(), NC_NETCDF4);
    try{
        if(this->is0_==NULL){
           throw "is0_ is pointing to NULL Array";
        }
        else{
            outFile.write_VecNC("is", this->is0_, this->m_ + 1, NC_INT64);
        }
        if(this->js_==NULL){
           throw "js_ is pointing to NULL Array";
        }
        else{
            outFile.write_VecNC("js", this->js_, this->is0_[m_], NC_INT64);
        }
        if(this->vs_==NULL){
           throw "vs_ is pointing to NULL Array";
        }
        else{
            outFile.write_VecNC("vs", this->vs_, this->is0_[m_], NC_FLOAT);
        }
        if(dims.empty() == true){
           throw "Error in dims in array";
        }
        else{
            outFile.write_AttNC("js", "dims", dims, NC_INT64);
        }
    }
    catch(const char* msg)
    {
        cerr<<"Matirx error: "<<msg<<endl;
        exit(2);
    };

#endif
}


void MatrixSparseMKL::elementwiseSqr()
{
#ifndef NDEBUG
	cout << getTimeStamp() << "   (X" << *this << ")^2 := ..." << endl;
#endif

    if (!isEmpty_)
    {
        vsSqr(is0_[m_], vs_, vs_);
    }
    
#ifndef NDEBUG
    cout << getTimeStamp() << "   ... X" << *this << endl;
#endif
}


void MatrixSparseMKL::elementwiseSqrt()
{
#ifndef NDEBUG
	cout << getTimeStamp() << "   sqrt(X" << *this << ") := ..." << endl;
#endif

    if (!isEmpty_)
    {
        vsSqrt(is0_[m_], vs_, vs_);
    }

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
    
    if (!isEmpty_)
    {
        ippsAddC_32f_I(beta, vs_, is0_[m_]);
    }

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

    if (!isEmpty_)
    {
        ippsMulC_32f_I(beta, vs_, is0_[m_]);
    }

#ifndef NDEBUG
    cout << getTimeStamp() << "   ... X" << *this << endl;
#endif
}


void MatrixSparseMKL::elementwiseMul(const MatrixSparseMKL& a)
{
#ifndef NDEBUG
    cout << getTimeStamp() << "   X" << *this << " * A" << a << " := ..." << endl;
#endif
    
    if (!isEmpty_)
    {
        vsMul(is0_[m_], vs_, a.vs_, vs_);
    }
    
#ifndef NDEBUG
    cout << getTimeStamp() << "   ... X" << *this << endl;
#endif
}


void MatrixSparseMKL::elementwiseDiv(const MatrixSparseMKL& a)
{
#ifndef NDEBUG
    cout << getTimeStamp() << "   X" << *this << " / A" << a << " := ..." << endl;
#endif
    
    if (!isEmpty_)
    {
        vsDiv(is0_[m_], vs_, a.vs_, vs_);
    }
    
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
    
    if (!isEmpty_)
    {
        ippsSum_32f(vs_, is0_[m_], &sum, ippAlgHintFast);
    }
    
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
    
    if (!isEmpty_)
    {
        ippsNorm_L2_32f(vs_, is0_[m_], &sum);
        sum *= sum;
    }
    
    /*if (isinf(sum))
    {
        for (ii nz = 0; nz < nnz(); nz++)
        {
            cout << fixed << vs_[nz] << endl;
        }
        cout << "ARGH" << endl;
        exit(0);
    }*/

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
    
    if (!isEmpty_)
    {
        ippsNormDiff_L2_32f(vs_, a.vs_, is0_[m_], &sum);
        sum *= sum;
    }
    
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
    cout << getTimeStamp() << "   subset(A" << a << ") := ..." << endl;
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
    cout << getTimeStamp() << "   ... X" << *this << endl;
#endif
}


// todo: use OpenMP
void MatrixSparseMKL::subsetElementwiseDiv(const MatrixSparseMKL& a)
{
#ifndef NDEBUG
    cout << getTimeStamp() << "   X" << *this << " / subset(A" << a << ") := ..." << endl;
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
        os << (a.isTransposed_ ? "t" : "") << "[" << a.m_ << "," << a.n_ << "]:" << a.nnz() << "/" << a.size() << ":";
        os.unsetf(ios::floatfield);
        os << setprecision(3) << 100.0 * a.nnz() / (double)a.size() << "%";
	}

	return  os;
}
