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


using namespace std;


MatrixSparse::MatrixSparse()
	: mat_(0), isOwned_(false)
{
}


MatrixSparse::~MatrixSparse()
{
	free();
}


void MatrixSparse::init(const MatrixSparse& a)
{
	free();

	m_ = a.m_;
	n_ = a.n_;
	is0_ = new ii[m_ + 1];
	is1_ = is0_ + 1;
	for (ii i = 0; i < m_; i++) is0_[i] = a.is0_[i];
	is1_[m_ - 1] = a.is1_[m_ - 1];
	js_ = new ii[is1_[m_ - 1]];
	for (ii nz = 0; nz < is1_[m_ - 1]; nz++) js_[nz] = a.js_[nz];
	vs_ = new fp[is1_[m_ - 1]];
	for (ii nz = 0; nz < is1_[m_ - 1]; nz++) vs_[nz] = a.vs_[nz];
	isOwned_ = true;

	mkl_sparse_s_create_csr(&mat_, SPARSE_INDEX_BASE_ZERO, m_, n_, is0_, is1_, js_, vs_);

#ifndef NDEBUG
	cout << "  A" << a << " = Y" << *this << endl;
#endif
}


void MatrixSparse::init(const MatrixSparse& a, fp pruneThreshold)
{
	free();

#ifndef NDEBUG
	cout << "  (A" << a << " >= ";
	cout.unsetf(ios::floatfield);
	cout << setprecision(8) << pruneThreshold << ") ? A : 0.0 = " << flush;
#endif

	m_ = a.m_;
	n_ = a.n_;
	is0_ = new ii[m_ + 1];
	is0_[0] = 0;
	is1_ = is0_ + 1;
	for (ii i = 0; i < m_; i++) is1_[i] = 0;

	ii nnz = 0;
	for (ii nza = 0; nza < a.nnz(); nza++)
	{
		if (a.vs_[nza] >= pruneThreshold)
		{
			nnz++;
		}
	}
	js_ = new ii[nnz];
	vs_ = new fp[nnz];

	ii nz = 0;
	ii ia = 0;
	for (ii nza = 0; nza < a.nnz(); nza++)
	{
		while (nza >= a.is1_[ia]) ia++; // row of nz'th non-zero in a 

		if (a.vs_[nza] >= pruneThreshold)
		{
			is1_[ia]++;
			js_[nz] = a.js_[nza];
			vs_[nz] = a.vs_[nza];
			nz++;
		}
	}
	for (ii i = 0; i < m_; i++)
	{
		is0_[i + 1] += is0_[i];
	}

	isOwned_ = true;

	mkl_sparse_s_create_csr(&mat_, SPARSE_INDEX_BASE_ZERO, m_, n_, is0_, is1_, js_, vs_);

#ifndef NDEBUG
	cout << "Y" << *this << endl;
#endif
}


void MatrixSparse::init(ii m, ii n)
{
	free();

	m_ = m;
	n_ = n;
	is0_ = new ii[m_ + 1];
	is0_[0] = 0;
	is1_ = is0_ + 1;
	for (ii i = 0; i < m_; i++) is1_[i] = 0;
	js_ = 0;
	vs_ = 0;
	isOwned_ = true;

	mkl_sparse_s_create_csr(&mat_, SPARSE_INDEX_BASE_ZERO, m_, n_, is0_, is1_, js_, vs_);

#ifndef NDEBUG
	cout << "  = Y" << *this << endl;
#endif
}


void MatrixSparse::init(ii m, ii n, ii nnz, const fp* acoo, const ii* rowind, const ii* colind)
{
	free();

	m_ = m;
	n_ = n;
	is0_ = new ii[m_ + 1];
	is1_ = is0_ + 1;
	js_ = new ii[nnz];
	vs_ = new fp[nnz];
	isOwned_ = true;

	ii job[] = { 2, 0, 0, 0, nnz, 0 };
	ii info;
	fp* c_acoo = const_cast<fp*>(acoo);
	ii* c_rowind = const_cast<ii*>(rowind);
	ii* c_colind = const_cast<ii*>(colind);
	mkl_scsrcoo(job, &m_, vs_, js_, is0_, &nnz, c_acoo, c_rowind, c_colind, &info);

	mkl_sparse_s_create_csr(&mat_, SPARSE_INDEX_BASE_ZERO, m_, n_, is0_, is1_, js_, vs_);

#ifndef NDEBUG
	cout << "  = Y" << *this << endl;
#endif
}


void MatrixSparse::init(ii m, ii n, fp v)
{
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


void MatrixSparse::convertFromDense(ii m, ii n, const fp* vs)
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


void MatrixSparse::convertToDense(fp* vs) const
{
	for (ii x = 0; x < m_ * n_; x++)
	{
		vs[x] = 0.0;
	}

	ii i = 0;
	for (ii nz = 0; nz < nnz(); nz++)
	{
		while (nz >= is1_[i]) i++; // row of nz'th non-zero 
		ii j = js_[nz]; // column of nz'th non-zero 

		vs[j + i * n_] = vs_[nz];
	}
}


void MatrixSparse::free()
{
	if (mat_)
	{
		mkl_sparse_destroy(mat_);
		mat_ = 0;

		if (isOwned_)
		{
			delete[] is0_;
			delete[] js_;
			delete[] vs_;
			isOwned_ = false;
		}
	}
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


bool MatrixSparse::operator!() const 
{ 
	return mat_ == 0; 
}


ii MatrixSparse::nnz(bool actual) const
{	
	if (actual)
	{
		ii nnz = 0;
		for (ii x = 0; x < is1_[m_ - 1]; x++)
		{
			if (vs_[x] != 0.0) nnz++;
		}
		return nnz;
	}
	else
	{
		return is1_[m_ - 1];
	}
}


li MatrixSparse::mem() const
{
	return sizeof(*this) + (li)nnz() * sizeof(fp) + (li)(m_ + 1) * sizeof(ii) + (li)nnz() * sizeof(ii);
}


void MatrixSparse::mul(const MatrixSparse& a, Transpose transpose, const MatrixSparse& x, Accumulate accumulate)
{
#ifndef NDEBUG
	cout << "  A" << (transpose == Transpose::NO ? "" : "t") << a << " %*% X" << x << (accumulate == Accumulate::NO ? " = " : " =+ ") << flush;
#endif

	if (!*this) accumulate = Accumulate::NO;

	// mkl can't handle completely empty csr matrix...
	if (x.nnz() == 0 || a.nnz() == 0)
	{
		if (accumulate == Accumulate::NO)
		{
			init(transpose == Transpose::NO ? a.m_ : a.n_, x.n_);
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

			mkl_sparse_spmm(transpose == Transpose::NO ? SPARSE_OPERATION_NON_TRANSPOSE : SPARSE_OPERATION_TRANSPOSE, a.mat_, x.mat_, &mat_);
			sparse_index_base_t indexing = SPARSE_INDEX_BASE_ZERO;
			mkl_sparse_s_export_csr(mat_, &indexing, &m_, &n_, &is0_, &is1_, &js_, &vs_);

			sparse_matrix_t t2;
			mkl_sparse_spmm(transpose == Transpose::NO ? SPARSE_OPERATION_NON_TRANSPOSE : SPARSE_OPERATION_TRANSPOSE, a.mat_, x.mat_, &t2);
			ii m; ii n; ii* is0; ii* is1; ii* js; fp* vs;
			indexing = SPARSE_INDEX_BASE_ZERO;
			mkl_sparse_s_export_csr(t2, &indexing, &m, &n, &is0, &is1, &js, &vs);

			for (ii nz = 0; nz < is1_[m_ - 1]; nz++)
			{
				if (vs_[nz] != vs[nz])
				{
					cout << "assert " << nz << ": " << vs_[nz] << " != " << vs[nz] << " transpose:" << (int)transpose << " accumulate:" << (int)accumulate << endl;
					exit(1);
				}
			}

		}
		else
		{
			sparse_matrix_t t, y;
			mkl_sparse_spmm(transpose == Transpose::NO ? SPARSE_OPERATION_NON_TRANSPOSE : SPARSE_OPERATION_TRANSPOSE, a.mat_, x.mat_, &t);
			ii m; ii n; ii* is0; ii* is1; ii* js; fp* vs;
			sparse_index_base_t indexing = SPARSE_INDEX_BASE_ZERO;
			mkl_sparse_s_export_csr(t, &indexing, &m, &n, &is0, &is1, &js, &vs);

			sparse_matrix_t t2;
			mkl_sparse_spmm(transpose == Transpose::NO ? SPARSE_OPERATION_NON_TRANSPOSE : SPARSE_OPERATION_TRANSPOSE, a.mat_, x.mat_, &t2);
			fp* vs2;
			indexing = SPARSE_INDEX_BASE_ZERO;
			mkl_sparse_s_export_csr(t2, &indexing, &m, &n, &is0, &is1, &js, &vs2);

			for (ii nz = 0; nz < is1[m - 1]; nz++)
			{
				if (vs[nz] != vs2[nz])
				{
					cout << "assert " << nz << ": " << vs[nz] << " != " << vs2[nz] << " transpose:" << (int)transpose << " accumulate:" << (int)accumulate << endl;
					exit(1);
				}
			}

			mkl_sparse_s_add(SPARSE_OPERATION_NON_TRANSPOSE, mat_, 1.0, t, &y);
			free();
			mkl_sparse_destroy(t);
			mat_ = y;

			indexing = SPARSE_INDEX_BASE_ZERO;
			mkl_sparse_s_export_csr(mat_, &indexing, &m_, &n_, &is0_, &is1_, &js_, &vs_);
		}
	}

#ifndef NDEBUG
	cout << "Y" << *this << endl;
#endif
}


void MatrixSparse::elementwiseSqr()
{
#ifndef NDEBUG
	cout << "  (Y" << *this << ")^2 = " << flush;
#endif

	vsSqr(nnz(), vs_, vs_);

#ifndef NDEBUG
	cout << "Y" << *this << endl;
#endif
}


void MatrixSparse::elementwiseSqrt()
{
#ifndef NDEBUG
	cout << "  sqrt(Y" << *this << ") = " << flush;
#endif

	vsSqrt(nnz(), vs_, vs_);

#ifndef NDEBUG
	cout << "Y" << *this << endl;
#endif
}


void MatrixSparse::elementwiseAdd(fp beta)
{
#ifndef NDEBUG
	cout << "  Y" << *this << " + ";
	cout.unsetf(ios::floatfield);
	cout << setprecision(8) << beta << " = " << flush;
#endif

	#pragma omp parallel for
	for (ii nz = 0; nz < nnz(); nz++)
	{
		vs_[nz] += beta;
	}

#ifndef NDEBUG
	cout << "Y" << *this << endl;
#endif
}


void MatrixSparse::elementwiseMul(fp beta)
{
#ifndef NDEBUG
	cout << "  Y" << *this << " * ";
	cout.unsetf(ios::floatfield);
	cout << setprecision(8) << beta << " = " << flush;
#endif

	#pragma omp parallel for
	for (ii nz = 0; nz < nnz(); nz++)
	{
		vs_[nz] *= beta;
	}

#ifndef NDEBUG
	cout << "Y" << *this << endl;
#endif
} 


void MatrixSparse::elementwiseMul(const MatrixSparse& a)
{
#ifndef NDEBUG
	cout << "  A" << a << " =* " << flush;
#endif

	ii i = 0;
	ii ia = 0;
	ii nza = 0;
	for (ii nz = 0; nz < is1_[m_ - 1]; nz++)
	{
		while (nz >= is1_[i]) i++; // row of nz'th non-zero

		for (; nza < a.is1_[m_ - 1]; nza++)
		{
			while (nza >= a.is1_[ia]) ia++; // row of nz'th non-zero

			if (a.js_[nza] == js_[nz] && ia == i)
			{
				vs_[nz] *= a.vs_[nza];
				break;
			}
		}
	}

#ifndef NDEBUG
	cout << "Y" << *this << endl;
#endif
}


void MatrixSparse::elementwiseDiv(const MatrixSparse& a)
{
#ifndef NDEBUG
	cout << "  A" << a << " =/ " << flush;
#endif

	ii i = 0;
	ii ia = 0;
	ii nza = 0;
	for (ii nz = 0; nz < is1_[m_ - 1]; nz++)
	{
		while (nz >= is1_[i]) i++; // row of nz'th non-zero

		for (; nza < a.is1_[m_ - 1]; nza++)
		{
			while (nza >= a.is1_[ia]) ia++; // row of nz'th non-zero

			if (a.js_[nza] == js_[nz] && ia == i)
			{
				vs_[nz] /= a.vs_[nza];
				break;
			}
		}
	}

#ifndef NDEBUG
	cout << "Y" << *this << endl;
#endif
}


double MatrixSparse::sum() const
{
	double sum = 0.0;

#ifndef NDEBUG
	cout << "  sum(Y" << *this << ") = " << flush;
#endif

	#pragma omp parallel for reduction(+:sum)
	for (ii nz = 0; nz < nnz(); nz++)
	{
		sum += vs_[nz];
	}

#ifndef NDEBUG
	cout.unsetf(ios::floatfield);
	cout << setprecision(8) << sum << endl;
#endif

	return sum;
}


double MatrixSparse::sumSqrs() const
{
	double sum = 0.0;

#ifndef NDEBUG
	cout << "  sum((Y" << *this << ")^2) = " << flush;
#endif

	#pragma omp parallel for reduction(+:sum)
	for (ii nz = 0; nz < nnz(); nz++)
	{
		sum += vs_[nz] * vs_[nz];
	}

#ifndef NDEBUG
	cout.unsetf(ios::floatfield);
	cout << setprecision(8) << sum << endl;
#endif

	return sum;
}


double MatrixSparse::sumSqrDiffs(const MatrixSparse& a) const
{
	double sum = 0.0;

#ifndef NDEBUG
	cout << "  sum((Y" << *this << " - A" << a << ")^2) = " << flush;
#endif

	ii i = 0;
	ii ia = 0;
	ii nza = 0;
	for (ii nz = 0; nz < is1_[m_ - 1]; nz++)
	{
		while (nz >= is1_[i]) i++; // row of nz'th non-zero

		for (; nza < a.is1_[m_ - 1]; nza++)
		{
			while (nza >= a.is1_[ia]) ia++; // row of nz'th non-zero

			if (a.js_[nza] == js_[nz] && ia == i)
			{
				sum += (vs_[nz] - a.vs_[nza]) * (vs_[nz] - a.vs_[nza]);
				break;
			}
		}
	}

#ifndef NDEBUG
	cout.unsetf(ios::floatfield);
	cout << setprecision(8) << sum << endl;
#endif

	return sum;
}


const fp* MatrixSparse::getVs() const
{
	return vs_;
}


ostream& operator<<(ostream& os, const MatrixSparse& a)
{
	if (a.m() == 0)
	{
		os << "[]";
	}
	else
	{
		os << "{" << a.m() << "," << a.n() << "}:(" << a.nnz(true) << "," << a.nnz() << ")/" << a.size() << ":";
        os.unsetf(ios::floatfield);
        os << setprecision(3) << 100.0 * a.nnz() / (double)a.size() << "%";
	}

	return  os;
}
