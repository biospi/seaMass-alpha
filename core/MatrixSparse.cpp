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


using namespace std;


MatrixSparse::MatrixSparse()
	: m_(0), n_(0), nnz_(0), vs_(0), is_(0), js_(0), isIsJsOwned_(false), isVsOwned_(false), mat_(0)
{
}


MatrixSparse::~MatrixSparse()
{
	free();
}


void MatrixSparse::init(ii m, ii n, ii nnz)
{
	free();

	m_ = m;
	n_ = n;
	nnz_ = nnz;
	is_ = new ii[m_ + 1];
	js_ = new ii[nnz_];
	isIsJsOwned_ = true;
	vs_ = new fp[nnz_];
	isVsOwned_ = true;

	mkl_sparse_s_create_csr(&mat_, SPARSE_INDEX_BASE_ZERO, m_, n_, is_, &(is_[1]), js_, vs_);
}


void MatrixSparse::init(const MatrixSparse& a)
{
	free();

	m_ = a.m_;
	n_ = a.n_;
	nnz_ = a.nnz_;
	is_ = a.is_;
	js_ = a.js_;
	isIsJsOwned_ = false;
	vs_ = new fp[nnz_];
	isVsOwned_ = true;

	mkl_sparse_s_create_csr(&mat_, SPARSE_INDEX_BASE_ZERO, m_, n_, is_, &(is_[1]), js_, vs_);
}


void MatrixSparse::init(ii m, ii n, ii nnz, const fp* acoo, const ii* rowind, const ii* colind)
{
	init(m, n, nnz);

	fp* c_acoo = const_cast<fp*>(acoo);
	ii* c_rowind = const_cast<ii*>(rowind);
	ii* c_colind = const_cast<ii*>(colind);
	ii job[] = { 2, 0, 0, 0, nnz_, 0 }; ii info;
	mkl_scsrcoo(job, &m_, vs_, js_, is_, &nnz, c_acoo, c_rowind, c_colind, &info);

	mkl_sparse_s_create_csr(&mat_, SPARSE_INDEX_BASE_ZERO, m_, n_, is_, &(is_[1]), js_, vs_);
}


void MatrixSparse::free()
{
	m_ = 0;
	n_ = 0;
	nnz_ = 0;

	if (isIsJsOwned_)
	{
		delete[] is_;
		delete[] js_;
		isIsJsOwned_ = false;
	}
	if (isVsOwned_)
	{
		delete[] vs_;
		isVsOwned_ = false;
	}

	if (mat_) mkl_sparse_destroy(mat_);
}


ii MatrixSparse::m() const 
{ 
	return m_; 
}


ii MatrixSparse::n() const 
{ 
	return n_; 
}


ii MatrixSparse::nnz() const  
{ 
	return nnz_; 
}


li MatrixSparse::size() const
{
	return (li)m_ * n_;
}

bool MatrixSparse::operator!() const 
{ 
	return vs_ == 0; 
}


li MatrixSparse::mem() const
{
	return sizeof(*this) + (li)nnz_ * sizeof(fp) + (li)(m_ + 1) * sizeof(ii) + (li)nnz_ * sizeof(ii);
}


void MatrixSparse::elementwiseSqr(const MatrixSparse& a)
{
#ifndef NDEBUG
	cout << "  (A" << a << ")^2 = " << flush;
#endif

	if (!*this) init(a);

	vsSqr(nnz_, a.vs_, vs_);

#ifndef NDEBUG
	cout << "Y" << *this << endl;
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
		os << "{" << a.m() << "," << a.n() << "}:" << a.nnz() << "/" << a.size() << ":";
        os.unsetf(ios::floatfield);
        os << setprecision(3) << 100.0 * a.nnz() / (double)a.size() << "%";
	}

	return  os;
}
